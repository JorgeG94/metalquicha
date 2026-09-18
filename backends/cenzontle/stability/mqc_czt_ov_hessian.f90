!! The electronic Hessian on the rotations that can move a converged SCF
module mqc_czt_ov_hessian
   !! `(A+B)` over the non-redundant occupied-virtual rotations, as an operator.
   !!
   !! A converged SCF is a stationary point. Whether it is a *minimum* is the
   !! question of whether the second derivative of the energy with respect to
   !! the orbital rotations is positive definite, and the matrix of those second
   !! derivatives is the electronic Hessian. Nothing here diagonalises it --
   !! `mqc_czt_stability` hands this operator to OpenTrustRegion, which does --
   !! and nothing here needs OpenTrustRegion, which is why the two are separate
   !! modules: the physics is compiled into every build with CPU integrals and
   !! only the eigensolver is optional.
   !!
   !! ## Where the Hessian comes from
   !!
   !! It is already in the program. `mqc_czt_response`'s operator drives the
   !! coupled-perturbed equations, and the matrix those equations are solved
   !! against *is* this Hessian. Recovering it is arithmetic on what
   !! `nuclear_apply` already computes, and the arithmetic is worth writing out
   !! because getting the factor or the sign wrong gives a curvature that is
   !! plausible and useless.
   !!
   !! Write `Delta_ai = e_a - e_i` for the orbital energy differences and
   !!
   !!     M_{ai,bj} = 4(ai|bj) - (ab|ij) - (aj|ib)
   !!
   !! for the two-electron part. `mqc_czt_cphf`'s `assemble_hessian` builds
   !! `(A+B) = Delta + M` densely from MO integrals, so that is the convention
   !! this program already uses and the one kept here.
   !!
   !! `nuclear_apply` forms the symmetrised trial density at double occupancy,
   !! contracts it as an SCF contracts a density, transforms back to the MO
   !! basis and then divides the virtual rows by `e_i - e_a`. Those steps
   !! amount to
   !!
   !!     apply(x)_ai = (M x)_ai / (e_i - e_a) = -(M x)_ai / Delta_ai
   !!
   !! which is what makes the coupled-perturbed equations a fixed point:
   !! `Delta x + M x = -h` rearranges to `x = -h/Delta + apply(x)`, and
   !! `solve_response` iterates exactly that. Multiplying back through by
   !! `Delta` recovers the Hessian:
   !!
   !!     (A+B) x |_ai = Delta_ai * (x_ai - apply(x)_ai)
   !!
   !! One `apply` -- one Fock build over the whole basis -- per Hessian-vector
   !! product, which is the cost a matrix-free eigensolver was chosen for.
   !!
   !! ## Why the vector is shorter than the response operator's
   !!
   !! `nuclear_response_t` works on `n_mo` by `n_occ`: a displaced nucleus drags
   !! its basis functions along, so orthonormality has to be maintained and the
   !! occupied-occupied block of the response is fixed rather than solved for.
   !! `nuclear_apply` zeroes those rows in its image and `nuclear_denominators`
   !! reports zero for them.
   !!
   !! For a stability analysis that layout is actively wrong. Rotating occupied
   !! orbitals among themselves does not change a closed-shell determinant, so
   !! those directions are *redundant*: the Hessian is identically zero on them.
   !! Handed the full layout an eigensolver would find `n_occ^2` exact zeros,
   !! report the lowest eigenvalue of the spectrum as zero, and call every
   !! reference marginally stable. So this operator is parameterised on the
   !! virtual-occupied rectangle alone -- `n_vir` by `n_occ`, the non-redundant
   !! rotations, the same space `mqc_czt_mcscf`'s `is_redundant` picks out --
   !! and embeds into the response operator's longer vector on the way in.
   !!
   !! Because the parameterisation carries no redundancy, OpenTrustRegion's
   !! `project` hook is not needed and is left null; see `mqc_czt_stability`.
   !!
   !! ## What this does and does not decide
   !!
   !! `(A+B)` is the **real singlet** orbital-rotation Hessian: it answers
   !! whether the reference is a minimum among real closed-shell determinants,
   !! the RHF-to-RHF question. A triplet instability -- the RHF-to-UHF question
   !! -- lives in the triplet Hessian, which combines the same integrals
   !! differently, and a complex instability lives in `(A-B)`. Neither is
   !! reached from here, and a caller must not read a stable verdict as more
   !! than it is.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response, only: response_operator_t
   use mqc_czt_hessian, only: nuclear_response_t
   implicit none
   private

   public :: ov_hessian_t
   public :: stability_result_t
   public :: build_scf_ov_hessian

   type :: stability_result_t
      !! What a stability analysis found, whichever eigensolver ran it
      logical :: ran = .false.
         !! Whether the analysis was attempted at all. False leaves every other
         !! field at its default and means the deck did not ask, or the build
         !! could not.
      logical :: stable = .true.
         !! Whether the converged SCF is a minimum with respect to real
         !! closed-shell orbital rotations. See the module header for what that
         !! does *not* cover.
      logical :: has_curvature = .false.
         !! Whether `lowest_curvature` was recovered. It is not always
         !! available; `mqc_czt_stability` says why.
      real(dp) :: lowest_curvature = 0.0_dp
         !! The smallest eigenvalue of `(A+B)`, in the convention of this
         !! module's header. Meaningful only when `has_curvature`.
      real(dp), allocatable :: rotation(:)
         !! The `n_vir` by `n_occ` rotation that lowers the energy, flattened
         !! virtual fastest. Allocated only when the reference is unstable --
         !! there is no downhill direction from a minimum.
      integer :: n_parameters = 0
         !! `n_vir*n_occ`, the size of the space that was searched.
      integer :: n_products = 0
         !! Hessian-vector products spent, each one a Fock build.
   end type stability_result_t

   type :: ov_hessian_t
      !! `(A+B)` as something that can multiply a vector
      !!
      !! Deliberately not an extension of `response_operator_t`: that interface
      !! is `Delta^-1 K` on the response operator's own longer layout, and this
      !! is the Hessian itself on the shorter one. Extending it would put two
      !! different meanings behind one `apply`.
      class(response_operator_t), pointer :: response => null()
         !! Where the Fock build happens. A pointer rather than a component so
         !! that a test can put an operator with a known matrix behind it; the
         !! target must outlive this object.
      real(dp), allocatable :: gaps(:)
         !! `e_a - e_i` over the occupied-virtual rectangle, virtual fastest.
         !! Also the Hessian diagonal this hands an eigensolver -- see
         !! `diagonal`.
      integer :: n_occ = 0
      integer :: n_vir = 0
      integer :: n_mo = 0
      type(error_t) :: error
         !! What `apply` could not do. The eigensolvers this is handed report
         !! integer codes and carry no room for a message, so the message is
         !! kept on the object that the callback reached in the first place.
      integer :: n_apply = 0
   contains
      procedure :: apply => ov_hessian_apply
      procedure :: diagonal => ov_hessian_diagonal
      procedure :: length => ov_hessian_length
   end type ov_hessian_t

contains

   pure function ov_hessian_length(this) result(n)
      !! How long a rotation vector is: the non-redundant rotations and no more
      class(ov_hessian_t), intent(in) :: this
      integer :: n

      n = this%n_vir*this%n_occ
   end function ov_hessian_length

   function ov_hessian_diagonal(this) result(diag)
      !! The diagonal an eigensolver preconditions and starts on
      !!
      !! `Delta` rather than `diag(Delta + M)`: the exact diagonal costs one
      !! Fock build per element, which is the whole matrix. Nothing in the
      !! answer depends on it -- a preconditioner and a starting guess move the
      !! iteration count and not the eigenvalue -- and `Delta` is the right
      !! approximation anyway, since its smallest element is the HOMO-LUMO
      !! rotation and that is where an instability usually is.
      class(ov_hessian_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      diag = this%gaps
   end function ov_hessian_diagonal

   subroutine ov_hessian_apply(this, x, hx)
      !! `(A+B) x`, from one application of the response operator
      !!
      !! `Delta * (x - apply(x))`, which the module header derives. The
      !! embedding into and out of the response operator's longer vector is the
      !! rest of it: occupied rows go in zero, because a rotation among
      !! occupied orbitals is not one of the parameters, and come back ignored,
      !! because the response operator zeroes them.
      class(ov_hessian_t), intent(inout) :: this
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hx(:)

      real(dp), allocatable :: full(:), image(:)
      integer :: i, a, n_ov, wide

      n_ov = this%length()
      if (size(x) /= n_ov .or. size(hx) /= n_ov) then
         call this%error%set(ERROR_VALIDATION, "a trial rotation handed to the "// &
                             "electronic Hessian is not the length of the "// &
                             "non-redundant occupied-virtual space")
         hx = 0.0_dp
         return
      end if
      if (.not. associated(this%response)) then
         call this%error%set(ERROR_VALIDATION, "the electronic Hessian has no "// &
                             "response operator behind it, so there is nothing "// &
                             "to build a Fock matrix with")
         hx = 0.0_dp
         return
      end if

      wide = this%n_mo*this%n_occ
      allocate (full(wide), image(wide))
      full = 0.0_dp
      do i = 1, this%n_occ
         do a = 1, this%n_vir
            full(this%n_occ + a + (i - 1)*this%n_mo) = x(a + (i - 1)*this%n_vir)
         end do
      end do

      call this%response%apply(full, image, this%error)
      this%n_apply = this%n_apply + 1
      if (this%error%has_error()) then
         hx = 0.0_dp
         return
      end if

      do i = 1, this%n_occ
         do a = 1, this%n_vir
            hx(a + (i - 1)*this%n_vir) = this%gaps(a + (i - 1)*this%n_vir) &
                                         *(x(a + (i - 1)*this%n_vir) &
                                           - image(this%n_occ + a + (i - 1)*this%n_mo))
         end do
      end do
      deallocate (full, image)
   end subroutine ov_hessian_apply

   subroutine build_scf_ov_hessian(mol, orbitals, energies, n_occ, response, hessian, &
                                   error, xc, reference, k_scale, rs_k_lr, rs_omega, bounds)
      !! Point an `ov_hessian_t` at a converged SCF
      !!
      !! `response` is an argument rather than a local because the operator has
      !! to outlive this call: `hessian%response` points at it, and a local
      !! would be gone before the first Hessian-vector product. The caller owns
      !! it and has to keep it in scope for as long as it uses `hessian`.
      !!
      !! `xc` and `reference` are one argument in two halves, refused
      !! separately for the reason `solve_mo1_batch` refuses them: the kernel
      !! has to be evaluated at a density, and a context without one reaches
      !! the Fock build with nothing allocated.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      type(nuclear_response_t), intent(out), target :: response
      type(ov_hessian_t), intent(out) :: hessian
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
         !! Exact exchange in the response operator: one for Hartree-Fock, the
         !! mixing fraction for a hybrid, zero for a pure functional.
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      real(dp), intent(in), optional :: bounds(:, :)

      integer :: n_ao, n_mo, n_vir, i, a

      if (error%has_error()) return

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ

      if (n_occ < 1 .or. n_vir < 1) then
         call error%set(ERROR_VALIDATION, "a stability analysis needs at least one "// &
                        "occupied and one virtual orbital; there are no rotations "// &
                        "to examine otherwise")
         return
      end if
      if (size(energies) < n_mo) then
         call error%set(ERROR_VALIDATION, "there are fewer orbital energies than "// &
                        "orbitals, so the electronic Hessian cannot be scaled")
         return
      end if
      if (present(xc) .neqv. present(reference)) then
         call error%set(ERROR_VALIDATION, "the electronic Hessian was given an "// &
                        "exchange-correlation context without the reference density "// &
                        "its kernel is evaluated at, or the reverse; it needs both "// &
                        "or neither")
         return
      end if

      response%mol => mol
      response%orbitals = orbitals
      response%c_occ = orbitals(:, 1:n_occ)
      response%energies = energies
      response%n_occ = n_occ
      response%n_mo = n_mo
      response%n_pert = 1
      allocate (response%zero_h(n_ao, n_ao))
      response%zero_h = 0.0_dp
      if (present(xc)) response%xc => xc
      if (present(reference)) response%reference = reference
      if (present(k_scale)) response%k_scale = k_scale
      if (present(rs_k_lr)) response%rs_k_lr = rs_k_lr
      if (present(rs_omega)) response%rs_omega = rs_omega
      if (present(bounds)) then
         response%bounds = bounds
      else
         call schwarz_bounds(mol, response%bounds, error)
         if (error%has_error()) return
      end if

      hessian%response => response
      hessian%n_occ = n_occ
      hessian%n_vir = n_vir
      hessian%n_mo = n_mo
      allocate (hessian%gaps(n_vir*n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            hessian%gaps(a + (i - 1)*n_vir) = energies(n_occ + a) - energies(i)
         end do
      end do

      if (any(hessian%gaps <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "an occupied orbital lies above a virtual "// &
                        "one, so the orbitals are not the aufbau solution and the "// &
                        "orbital-energy denominators a stability analysis "// &
                        "preconditions on are not positive")
         return
      end if
   end subroutine build_scf_ov_hessian

end module mqc_czt_ov_hessian
