!! Transition moments, oscillator strengths and natural transition orbitals
module mqc_czt_tddft_properties
   !! What a linear-response spectrum is read for, once the energies are in
   !! hand: how bright each root is, which way it is polarised, and which
   !! orbital pair it actually is.
   !!
   !! Everything here is arithmetic over the amplitudes `mqc_czt_tddft`
   !! returns and two one-electron integral sets. No Fock build, no
   !! quadrature and no iteration, so the whole module costs a fraction of
   !! one matrix-vector product of the solve that produced its input.
   !!
   !! ## The amplitude convention this is written for
   !!
   !! `|X|^2 - |Y|^2 = 1/2`, every route, which is what `response_excitations`
   !! hands back. A closed-shell excitation is two spin-orbital excitations of
   !! equal weight, the spatial amplitude carries both, and the factor two in
   !! every moment below is that spin sum. Feeding unit-normalised amplitudes
   !! in instead would multiply every moment by the square root of two and
   !! every oscillator strength by two, silently.
   !!
   !! ## The two gauges
   !!
   !! **Length**, about the nuclear charge centroid `R0 = sum Z_A R_A / sum Z_A`:
   !!
   !!     mu = 2 sum_ia <i| r - R0 |a> (X+Y)_ia
   !!     f_len = (2/3) w |mu|^2
   !!
   !! The origin matters only through the charge of the transition density,
   !! which is zero -- `<i|a>` vanishes by orthogonality -- so `mu` is origin
   !! independent to round-off and `R0` is a convention rather than a choice.
   !! It is PySCF's, so a cross-code comparison is of the same number.
   !!
   !! **Velocity**, from the same amplitudes' difference:
   !!
   !!     v = - [ 2 sum_ia <grad i|a> (X-Y)_ia ]
   !!     f_vel = (2/3) |v|^2 / w
   !!
   !! and the sign is the only thing here that is not obvious, so it is
   !! written out. libcint's `int1e_ipovlp` is `(grad mu | nu)`, the gradient
   !! on the **bra**, which is `-<mu| grad |nu>`: the operator is
   !! anti-Hermitian over real functions and integration by parts moves it
   !! across at the cost of a sign. PySCF contracts that matrix as it comes
   !! and negates the result (`transition_velocity_dipole`, and the momentum
   !! operator is `p = -i grad` so what both codes report is the imaginary
   !! part of `<0|p|n>`). This follows it exactly, which makes `v` comparable
   !! to PySCF's component by component and not merely in magnitude.
   !!
   !! The two gauges agree only in a complete basis. In cc-pVDZ they differ
   !! by a factor of three on the lowest root of water, which is a statement
   !! about the basis and not about either implementation; both are reported
   !! because the gap between them is the diagnostic.
   !!
   !! ## Natural transition orbitals
   !!
   !! The singular value decomposition of the `(n_occ, n_vir)` excitation
   !! amplitude renormalised to a unit vector -- `Y` is dropped, which is
   !! Martin's definition and PySCF's. The weights are the squared singular
   !! values and sum to one, so the leading one says how nearly the root is a
   !! single orbital pair. The orbitals themselves are `C_occ U` and
   !! `C_vir V`, with each column's largest component made positive so the
   !! phase is reproducible.
   !!
   !! ## Triplets
   !!
   !! A spin-forbidden transition has no dipole at all: the spatial integral
   !! is multiplied by an overlap of orthogonal spin functions. So a triplet
   !! root's moments are set to zero rather than computed, and the zero is
   !! exact rather than the result of a cancellation that a finite basis
   !! would leave at 1e-16.
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_gesvd
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_multipole, only: multipole_matrices, DIPOLE_COMPONENTS
   use mqc_czt_gradient, only: one_electron_deriv, DERIV_OVLP
   use mqc_result_types, only: STATE_SPIN_TRIPLET
   implicit none
   private

   public :: excited_properties_t
   public :: excited_properties
   public :: natural_transition_orbitals
   public :: nuclear_charge_centroid
   public :: log_property_table

   real(dp), parameter :: SPIN_SUM = 2.0_dp
      !! The closed-shell spin sum in every transition moment.
      !!
      !! Two, because a spatial-orbital amplitude at `|X|^2-|Y|^2 = 1/2`
      !! stands for the alpha and beta spin-orbital excitations together. The
      !! factor and the normalisation are one convention in two places and
      !! move together or not at all.

   real(dp), parameter :: OSCILLATOR_PREFACTOR = 2.0_dp/3.0_dp
      !! The isotropic average in `f = (2/3) w |mu|^2`: one third from
      !! averaging over the three polarisation directions, two from the
      !! definition.

   type :: excited_properties_t
      !! One transition property set per root, in the solver's own order
      !!
      !! Every array runs over the same states as the excitation energies it
      !! was built from, and they are allocated together. A triplet row is
      !! exactly zero in all four moment arrays; `nto_weights` is real for a
      !! triplet, because the orbital content of a spin-forbidden excitation
      !! is perfectly well defined.
      real(dp), allocatable :: transition_dipole(:, :)
         !! (3, n_states) length gauge, atomic units, about `origin`
      real(dp), allocatable :: velocity_moment(:, :)
         !! (3, n_states) the imaginary part of `<0|p|n>`, atomic units
      real(dp), allocatable :: f_length(:)         !! (n_states) dimensionless
      real(dp), allocatable :: f_velocity(:)       !! (n_states) dimensionless
      real(dp), allocatable :: nto_weights(:, :)
         !! (n_pairs, n_states) descending, summing to one down each column.
         !! `n_pairs` is `min(n_occ, n_vir)`, which is `n_occ` for any basis
         !! worth running.
      real(dp), allocatable :: total_energy(:)
         !! (n_states) `E_SCF + w`, the excited state's own total energy
      real(dp) :: origin(3) = 0.0_dp
         !! Where the length-gauge dipole was measured from, Bohr
   contains
      procedure :: destroy => properties_destroy
   end type excited_properties_t

contains

   subroutine properties_destroy(this)
      !! Release everything this holds and reset the origin
      class(excited_properties_t), intent(inout) :: this

      if (allocated(this%transition_dipole)) deallocate (this%transition_dipole)
      if (allocated(this%velocity_moment)) deallocate (this%velocity_moment)
      if (allocated(this%f_length)) deallocate (this%f_length)
      if (allocated(this%f_velocity)) deallocate (this%f_velocity)
      if (allocated(this%nto_weights)) deallocate (this%nto_weights)
      if (allocated(this%total_energy)) deallocate (this%total_energy)
      this%origin = 0.0_dp
   end subroutine properties_destroy

   pure function nuclear_charge_centroid(charges, coords) result(origin)
      !! `sum Z_A R_A / sum Z_A`, in Bohr
      !!
      !! The charges are the ones the molecule presents, so an atom behind an
      !! effective core potential contributes its reduced charge. That is
      !! PySCF's `_charge_center` as well, which is the point: the transition
      !! dipole of a neutral transition density does not depend on the origin,
      !! but any disagreement in it would show up here first and be blamed on
      !! the integrals.
      real(dp), intent(in) :: charges(:)      !! (natm)
      real(dp), intent(in) :: coords(:, :)    !! (3, natm)
      real(dp) :: origin(3)

      real(dp) :: total
      integer :: atom

      origin = 0.0_dp
      total = sum(charges)
      if (abs(total) <= 0.0_dp) return
      do atom = 1, size(charges)
         origin = origin + charges(atom)*coords(:, atom)
      end do
      origin = origin/total
   end function nuclear_charge_centroid

   subroutine mo_block(matrices, c_occ, c_vir, blocks)
      !! `<i| O_c |a>` for each component of a one-electron operator
      !!
      !! Two GEMMs per component rather than a quadruple loop: the AO matrix
      !! is `n_ao` squared and the occupied-virtual block is much smaller, so
      !! transforming is cheaper than anything that touches the AO pair list
      !! per amplitude.
      real(dp), intent(in) :: matrices(:, :, :)   !! (n_ao, n_ao, n_comp)
      real(dp), intent(in) :: c_occ(:, :)         !! (n_ao, n_occ)
      real(dp), intent(in) :: c_vir(:, :)         !! (n_ao, n_vir)
      real(dp), allocatable, intent(out) :: blocks(:, :, :)
         !! (n_occ, n_vir, n_comp)

      real(dp), allocatable :: half(:, :)
      integer :: n_ao, n_occ, n_vir, n_comp, comp

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      n_vir = size(c_vir, 2)
      n_comp = size(matrices, 3)

      allocate (blocks(n_occ, n_vir, n_comp), half(n_ao, n_vir))
      do comp = 1, n_comp
         call pic_gemm(matrices(:, :, comp), c_vir, half)
         call pic_gemm(c_occ, half, blocks(:, :, comp), transa="T")
      end do
      deallocate (half)
   end subroutine mo_block

   pure function moment_of(blocks, amplitude, n_occ, n_vir) result(moment)
      !! `2 sum_ia <i|O_c|a> t_ia` for the three Cartesian components
      !!
      !! The amplitude arrives flat at `idx = (i-1)*n_vir + a`, virtual
      !! fastest, which is the layout the solver and `response_product`
      !! share; the transformed integrals are `(i, a)`, so the two indices
      !! are read out rather than reshaped.
      real(dp), intent(in) :: blocks(:, :, :)   !! (n_occ, n_vir, 3)
      real(dp), intent(in) :: amplitude(:)      !! (n_occ*n_vir) flat
      integer, intent(in) :: n_occ, n_vir
      real(dp) :: moment(3)

      integer :: comp, i, a

      moment = 0.0_dp
      do comp = 1, 3
         do i = 1, n_occ
            do a = 1, n_vir
               moment(comp) = moment(comp) + &
                              blocks(i, a, comp)*amplitude((i - 1)*n_vir + a)
            end do
         end do
      end do
      moment = SPIN_SUM*moment
   end function moment_of

   subroutine natural_transition_orbitals(amplitude, c_occ, c_vir, weights, &
                                          nto_occ, nto_vir, error)
      !! The natural transition orbital pairs of one root, by SVD of `X`
      !!
      !! `X` is renormalised to a unit vector first, which is what makes the
      !! weights sum to one: the stored amplitude is at `|X|^2 = 1/2` for a
      !! Tamm-Dancoff root and at something between that and a half for a
      !! paired one, and neither is the CIS coefficient Martin's construction
      !! is written for. `Y` is not used -- the transition density matrix
      !! that would include it is a different object and breaks the point
      !! group symmetry of the orbitals; PySCF makes the same choice and says
      !! the same thing about it.
      !!
      !! The phase of each column is fixed by making its largest-magnitude
      !! component positive. Without that the orbitals are only defined up to
      !! a sign per pair and no two runs need agree.
      real(dp), intent(in) :: amplitude(:)     !! (n_occ*n_vir) flat `X`
      real(dp), intent(in) :: c_occ(:, :)      !! (n_ao, n_occ)
      real(dp), intent(in) :: c_vir(:, :)      !! (n_ao, n_vir)
      real(dp), allocatable, intent(out) :: weights(:)
         !! (min(n_occ, n_vir)) squared singular values, descending, sum one
      real(dp), allocatable, intent(out) :: nto_occ(:, :)
         !! (n_ao, n_pairs) the occupied natural transition orbitals
      real(dp), allocatable, intent(out) :: nto_vir(:, :)
         !! (n_ao, n_pairs) their virtual partners, pair by pair
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: t(:, :), u(:, :), vt(:, :), v(:, :), sigma(:)
      real(dp) :: norm
      integer :: n_ao, n_occ, n_vir, n_pairs, i, a, k, info

      if (error%has_error()) return

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      n_vir = size(c_vir, 2)
      n_pairs = min(n_occ, n_vir)

      if (size(amplitude) /= n_occ*n_vir) then
         call error%set(ERROR_VALIDATION, "the amplitude handed to the natural "// &
                        "transition orbitals is not the length of the "// &
                        "occupied-virtual space")
         return
      end if

      norm = sqrt(dot_product(amplitude, amplitude))
      if (norm <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "a root's excitation amplitude is "// &
                        "identically zero, so it has no natural transition orbitals")
         return
      end if

      allocate (t(n_occ, n_vir))
      do a = 1, n_vir
         do i = 1, n_occ
            t(i, a) = amplitude((i - 1)*n_vir + a)/norm
         end do
      end do

      allocate (u(n_occ, n_pairs), vt(n_pairs, n_vir), sigma(n_pairs))
      ! `t` is destroyed here, which is why it was built rather than aliased.
      call pic_gesvd(t, sigma, u, vt, info=info)
      deallocate (t)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "the singular value decomposition of an "// &
                        "excitation amplitude failed with LAPACK info "// &
                        to_char(info))
         return
      end if

      allocate (v(n_vir, n_pairs))
      do k = 1, n_pairs
         v(:, k) = vt(k, :)
      end do
      deallocate (vt)

      call fix_column_phases(u)
      call fix_column_phases(v)

      allocate (weights(n_pairs), nto_occ(n_ao, n_pairs), nto_vir(n_ao, n_pairs))
      weights = sigma*sigma
      call pic_gemm(c_occ, u, nto_occ)
      call pic_gemm(c_vir, v, nto_vir)
      deallocate (u, v, sigma)
   end subroutine natural_transition_orbitals

   subroutine fix_column_phases(columns)
      !! Make each column's largest-magnitude entry positive
      real(dp), intent(inout) :: columns(:, :)

      integer :: k, row

      do k = 1, size(columns, 2)
         row = maxloc(abs(columns(:, k)), dim=1)
         if (columns(row, k) < 0.0_dp) columns(:, k) = -columns(:, k)
      end do
   end subroutine fix_column_phases

   subroutine excited_properties(mol, orbitals, n_occ, excitations, state_spin, &
                                 x_amplitudes, y_amplitudes, scf_energy, props, error)
      !! Every transition property of a converged spectrum
      !!
      !! One pass over the dipole integrals and one over `int1e_ipovlp`,
      !! transformed once and contracted against every root -- the integrals
      !! do not depend on the state, and evaluating them per state is the
      !! obvious way to make a cheap analysis cost more than the solve.
      !!
      !! A triplet root short-circuits: its moments are written as zeros and
      !! no contraction is performed for it. Its natural transition orbitals
      !! are computed like any other root's.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: excitations(:)      !! (n_states) Hartree
      integer, intent(in) :: state_spin(:)        !! (n_states) `STATE_SPIN_*`
      real(dp), intent(in) :: x_amplitudes(:, :)  !! (n_ov, n_states)
      real(dp), intent(in) :: y_amplitudes(:, :)  !! (n_ov, n_states)
      real(dp), intent(in) :: scf_energy
         !! The reference total energy the excitations sit above, Hartree.
      type(excited_properties_t), intent(out) :: props
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dipole_ao(:, :, :), nabla_ao(:, :, :)
      real(dp), allocatable :: dipole_mo(:, :, :), nabla_mo(:, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      real(dp), allocatable :: sum_amplitude(:), diff_amplitude(:)
      real(dp), allocatable :: weights(:), nto_occ(:, :), nto_vir(:, :)
      integer :: n_mo, n_vir, n_ov, n_states, n_pairs, k

      if (error%has_error()) return

      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      n_ov = n_occ*n_vir
      n_states = size(excitations)
      n_pairs = min(n_occ, n_vir)

      if (n_occ < 1 .or. n_vir < 1) then
         call error%set(ERROR_VALIDATION, "transition properties need at least one "// &
                        "occupied and one virtual orbital")
         return
      end if
      if (size(x_amplitudes, 1) /= n_ov .or. size(y_amplitudes, 1) /= n_ov) then
         call error%set(ERROR_VALIDATION, "the amplitudes handed to the transition "// &
                        "properties are not the length of the occupied-virtual space")
         return
      end if
      if (size(x_amplitudes, 2) < n_states .or. size(y_amplitudes, 2) < n_states &
          .or. size(state_spin) < n_states) then
         call error%set(ERROR_VALIDATION, "there are fewer amplitude columns or spin "// &
                        "labels than there are excitation energies")
         return
      end if

      allocate (props%transition_dipole(3, n_states), props%velocity_moment(3, n_states))
      allocate (props%f_length(n_states), props%f_velocity(n_states))
      allocate (props%nto_weights(n_pairs, n_states), props%total_energy(n_states))
      props%transition_dipole = 0.0_dp
      props%velocity_moment = 0.0_dp
      props%f_length = 0.0_dp
      props%f_velocity = 0.0_dp
      props%nto_weights = 0.0_dp
      props%total_energy = scf_energy + excitations
      if (n_states < 1) return

      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)

      props%origin = nuclear_charge_centroid(mol%charges, mol%coords)
      call multipole_matrices(mol, props%origin, 1, dipole_ao, error)
      if (error%has_error()) return
      if (size(dipole_ao, 3) /= DIPOLE_COMPONENTS) then
         call error%set(ERROR_GENERIC, "the dipole integrals came back with a "// &
                        "component count that is not three")
         return
      end if
      call one_electron_deriv(mol, nabla_ao, DERIV_OVLP)

      call mo_block(dipole_ao, c_occ, c_vir, dipole_mo)
      call mo_block(nabla_ao, c_occ, c_vir, nabla_mo)
      deallocate (dipole_ao, nabla_ao)

      allocate (sum_amplitude(n_ov), diff_amplitude(n_ov))
      do k = 1, n_states
         if (state_spin(k) /= STATE_SPIN_TRIPLET) then
            sum_amplitude = x_amplitudes(:, k) + y_amplitudes(:, k)
            diff_amplitude = x_amplitudes(:, k) - y_amplitudes(:, k)
            props%transition_dipole(:, k) = moment_of(dipole_mo, sum_amplitude, &
                                                      n_occ, n_vir)
            ! The negation is libcint's gradient sitting on the bra; see the
            ! module header, which is where the whole sign argument lives.
            props%velocity_moment(:, k) = -moment_of(nabla_mo, diff_amplitude, &
                                                     n_occ, n_vir)
            props%f_length(k) = OSCILLATOR_PREFACTOR*excitations(k)* &
                                sum(props%transition_dipole(:, k)**2)
            if (excitations(k) > 0.0_dp) then
               props%f_velocity(k) = OSCILLATOR_PREFACTOR* &
                                     sum(props%velocity_moment(:, k)**2)/excitations(k)
            end if
         end if

         call natural_transition_orbitals(x_amplitudes(:, k), c_occ, c_vir, &
                                          weights, nto_occ, nto_vir, error)
         if (error%has_error()) return
         props%nto_weights(:, k) = weights
         deallocate (weights, nto_occ, nto_vir)
      end do
      deallocate (sum_amplitude, diff_amplitude, dipole_mo, nabla_mo, c_occ, c_vir)
   end subroutine excited_properties

   subroutine log_property_table(props, excitations, state_spin)
      !! The spectrum as a brightness table, one line per root
      !!
      !! Both gauges, because the gap between them is the diagnostic: they
      !! agree only in a complete basis, and a length-gauge strength three
      !! times its velocity-gauge partner is a statement about the basis
      !! rather than about the solve. The leading natural transition orbital
      !! weight is the last column because it is what says whether the
      !! dominant-amplitude list printed beside the energies is the whole
      !! story.
      type(excited_properties_t), intent(in) :: props
      real(dp), intent(in) :: excitations(:)
      integer, intent(in) :: state_spin(:)

      character(len=MAX_LINE_LENGTH) :: line
      integer :: k

      if (size(excitations) < 1) return
      if (.not. allocated(props%f_length)) return

      call logger%info("  transition moments are in atomic units, about the "// &
                       "nuclear charge centroid; a triplet's are exactly zero")
      call logger%info("   state     f(length)  f(velocity)        mu_x        "// &
                       "mu_y        mu_z    lead NTO")
      do k = 1, size(excitations)
         write (line, "(a,i4,2f13.6,3f12.6,f12.6)") "   ", k, props%f_length(k), &
            props%f_velocity(k), props%transition_dipole(1, k), &
            props%transition_dipole(2, k), props%transition_dipole(3, k), &
            leading_weight(props, k)
         if (state_spin(k) == STATE_SPIN_TRIPLET) then
            line = trim(line)//"   (spin forbidden)"
         end if
         call logger%info(trim(line))
      end do
   end subroutine log_property_table

   pure function leading_weight(props, k) result(weight)
      !! The largest natural transition orbital weight of root `k`, or zero
      type(excited_properties_t), intent(in) :: props
      integer, intent(in) :: k
      real(dp) :: weight

      weight = 0.0_dp
      if (.not. allocated(props%nto_weights)) return
      if (size(props%nto_weights, 1) < 1) return
      weight = maxval(props%nto_weights(:, k))
   end function leading_weight

end module mqc_czt_tddft_properties
