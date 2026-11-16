!! The potential integral contracted over a batch of points, for seminumerical exchange
module mqc_czt_esp_batch
   !! `G(g, nu) = sum_lambda <nu| 1/|r - r_g| |lambda> F(g, lambda)` for a batch
   !! of points, without the integral library.
   !!
   !! The same quantity `cosx_esp_contract` gets from `int1e_grids`, by the
   !! McMurchie-Davidson scheme laid out so that every inner loop runs over the
   !! points of the batch:
   !!
   !!     <a| 1/|r - C| |b> = (2 pi / p) K_ab sum_tuv E_tuv R_tuv(p, P - C)
   !!
   !! per primitive pair, with `E` the Hermite expansion of the product (no
   !! dependence on the point) and `R` the Hermite Coulomb integrals (all of
   !! the dependence on it), built from the Boys function by recurrence.
   !!
   !! Everything between the two transforms is Cartesian. `F` is taken to the
   !! Cartesian functions of each shell once per batch, the shell pairs
   !! accumulate into a Cartesian `G`, and that is taken back once at the end,
   !! so no shell pair pays for a spherical transform.
   !!
   !! Conventions are libcint's, because the basis functions are: the Cartesian
   !! order within a shell, the s and p factor it keeps outside its transform
   !! table, and that table itself, all through `mqc_czt_ao_data`.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: PI
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_ao_data, only: C2S_LMAX, c2s_block, common_fac_sp
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
   implicit none
   private

   public :: esp_shells_t
   public :: build_esp_shells
   public :: esp_batch_contract
   public :: boys_values
   public :: ESP_BATCH_LMAX

   integer, parameter :: ESP_BATCH_LMAX = C2S_LMAX
      !! Highest angular momentum of a shell this evaluates

   integer, parameter :: L_PAIR_MAX = 2*ESP_BATCH_LMAX
      !! Highest total angular momentum of a shell pair, and so the highest
      !! Boys order and Hermite degree
   integer, parameter :: N_CART_MAX = (ESP_BATCH_LMAX + 1)*(ESP_BATCH_LMAX + 2)/2
   integer, parameter :: N_HERM_MAX = (L_PAIR_MAX + 1)*(L_PAIR_MAX + 2)*(L_PAIR_MAX + 3)/6

   real(dp), parameter :: TWO_PI = 2.0_dp*PI

   ! The Boys function is tabulated on a uniform mesh and expanded about the
   ! nearest node. Half a step is 1/128, so five terms leave (1/128)^5/5! ~
   ! 2e-13 of the leading one; past the table the asymptotic form is used,
   ! whose neglected part is of order exp(-x). The step is a power of two so
   ! that a node's abscissa is exact.
   real(dp), parameter :: BOYS_STEP = 1.0_dp/64.0_dp
   real(dp), parameter :: BOYS_INV_STEP = 64.0_dp
   real(dp), parameter :: BOYS_X_MAX = 36.0_dp
   integer, parameter :: BOYS_N_NODE = 2304       !! BOYS_X_MAX/BOYS_STEP
   integer, parameter :: BOYS_TAYLOR = 4          !! Highest power kept
   integer, parameter :: BOYS_SERIES_TERMS = 400  !! For the table itself

   real(dp), parameter :: PRIM_EXPONENT_CUT = 50.0_dp
      !! A primitive pair is dropped when its Gaussian product prefactor is
      !! below exp(-this)

   type :: esp_shells_t
      !! A molecule's shells laid out for the batch evaluator
      !!
      !! One entry per contraction: a generally contracted shell appears once
      !! per contraction column, sharing its exponents.
      integer :: n_shell = 0
      integer :: n_cart = 0                    !! Cartesian functions in all
      integer :: n_ao = 0                      !! Functions in the molecule's own form
      integer :: l_max = 0
      logical :: cartesian = .false.           !! The molecule's angular form
      integer, allocatable :: l(:)             !! (n_shell)
      integer, allocatable :: parent(:)        !! (n_shell) shell of the molecule, 1-based
      integer, allocatable :: n_prim(:)        !! (n_shell)
      integer, allocatable :: prim_first(:)    !! (n_shell) 0-based start in `expo`/`coef`
      integer, allocatable :: cart_first(:)    !! (n_shell) 0-based first Cartesian function
      integer, allocatable :: ao_first(:)      !! (n_shell) 0-based first function, own form
      real(dp), allocatable :: centre(:, :)    !! (3, n_shell), Bohr
      real(dp), allocatable :: expo(:)
      real(dp), allocatable :: coef(:)         !! Contraction times libcint's s/p factor
      real(dp), allocatable :: c2s(:, :, :)
         !! (2l+1, n_cart(l), 0:l_max) spherical rows over Cartesian columns
      integer, allocatable :: powers(:, :, :)
         !! (3, n_cart(l), 0:l_max) the exponents of x, y, z for each Cartesian
         !! function of a shell, libcint's order
      integer, allocatable :: h_of(:, :, :)
         !! (0:2*l_max, 0:2*l_max, 0:2*l_max) position of the Hermite function
         !! (t, u, v) in a flat list ordered by degree, so one table serves
         !! every shell pair
      real(dp), allocatable :: boys(:, :)
         !! (0:L_PAIR_MAX + BOYS_TAYLOR, 0:BOYS_N_NODE) the Boys function at the
         !! mesh nodes, every order the expansion reaches. Order first, so the
         !! orders one point needs sit together.
   end type esp_shells_t

contains

   subroutine build_esp_shells(mol, shells, error)
      !! Lay a molecule's shells out for `esp_batch_contract`
      !!
      !! Refuses a basis with a shell above `ESP_BATCH_LMAX`; the caller keeps
      !! the integral library for that molecule.
      type(czt_molecule_t), intent(in) :: mol
      type(esp_shells_t), intent(out) :: shells
      type(error_t), intent(inout) :: error

      integer :: ish, ic, ip, l, nprim, nctr, exp_ptr, coef_ptr, iatom
      integer :: n_shell, n_prim_total, is, ipr, n_here, ao_off

      n_shell = 0
      n_prim_total = 0
      shells%l_max = 0
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         if (l > ESP_BATCH_LMAX) then
            call error%set(ERROR_VALIDATION, "the batch potential integral stops at "// &
                           "g functions and this basis goes higher")
            return
         end if
         shells%l_max = max(shells%l_max, l)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         n_shell = n_shell + nctr
         n_prim_total = n_prim_total + nctr*mol%bas(LIBCINT_NPRIM_OF, ish)
      end do

      shells%n_shell = n_shell
      shells%cartesian = mol%cartesian
      shells%n_ao = mol%nao
      allocate (shells%l(n_shell), shells%parent(n_shell), shells%n_prim(n_shell), &
                shells%prim_first(n_shell), shells%cart_first(n_shell), &
                shells%ao_first(n_shell), shells%centre(3, n_shell))
      allocate (shells%expo(n_prim_total), shells%coef(n_prim_total))

      is = 0
      ipr = 0
      shells%n_cart = 0
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         coef_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1
         if (mol%cartesian) then
            n_here = (l + 1)*(l + 2)/2
         else
            n_here = 2*l + 1
         end if
         ao_off = mol%shell_offset(ish)
         ! libcint lays a shell's functions out with the contraction index
         ! outermost.
         do ic = 1, nctr
            is = is + 1
            shells%l(is) = l
            shells%parent(is) = ish
            shells%n_prim(is) = nprim
            shells%prim_first(is) = ipr
            shells%cart_first(is) = shells%n_cart
            shells%ao_first(is) = ao_off + (ic - 1)*n_here
            shells%centre(:, is) = mol%coords(:, iatom)
            do ip = 1, nprim
               shells%expo(ipr + ip) = mol%env(exp_ptr + ip)
               shells%coef(ipr + ip) = common_fac_sp(l)*mol%env(coef_ptr + (ic - 1)*nprim + ip)
            end do
            ipr = ipr + nprim
            shells%n_cart = shells%n_cart + (l + 1)*(l + 2)/2
         end do
      end do

      allocate (shells%c2s(2*shells%l_max + 1, (shells%l_max + 1)*(shells%l_max + 2)/2, &
                           0:shells%l_max))
      shells%c2s = 0.0_dp
      do l = 0, shells%l_max
         call c2s_block(l, shells%c2s(1:2*l + 1, 1:(l + 1)*(l + 2)/2, l))
      end do

      allocate (shells%powers(3, (shells%l_max + 1)*(shells%l_max + 2)/2, 0:shells%l_max))
      shells%powers = 0
      do l = 0, shells%l_max
         call cartesian_powers(l, shells%powers(:, 1:(l + 1)*(l + 2)/2, l))
      end do
      allocate (shells%h_of(0:2*shells%l_max, 0:2*shells%l_max, 0:2*shells%l_max))
      call hermite_index(2*shells%l_max, shells%h_of)

      call tabulate_boys(shells%boys)
   end subroutine build_esp_shells

   pure subroutine tabulate_boys(table)
      !! The Boys function at the mesh nodes, for every order the expansion needs
      !!
      !! From `F_n(x) = exp(-x) sum_k (2x)^k / ((2n+1)(2n+3)...(2n+2k+1))`, whose
      !! terms are all positive, so nothing cancels at large `x`.
      real(dp), allocatable, intent(out) :: table(:, :)

      integer :: i, n, k
      real(dp) :: x, term, total

      allocate (table(0:L_PAIR_MAX + BOYS_TAYLOR, 0:BOYS_N_NODE))
      do i = 0, BOYS_N_NODE
         x = BOYS_STEP*real(i, dp)
         do n = 0, L_PAIR_MAX + BOYS_TAYLOR
            term = 1.0_dp/real(2*n + 1, dp)
            total = term
            do k = 1, BOYS_SERIES_TERMS
               term = term*2.0_dp*x/real(2*n + 2*k + 1, dp)
               total = total + term
               if (term < 1.0e-17_dp*total) exit
            end do
            table(n, i) = exp(-x)*total
         end do
      end do
   end subroutine tabulate_boys

   pure subroutine boys_values(table, n_max, x, values)
      !! `F_n(x)` for `n = 0 .. n_max` at each of a set of arguments
      real(dp), intent(in) :: table(0:, 0:)   !! From `build_esp_shells`
      integer, intent(in) :: n_max            !! At most `2*ESP_BATCH_LMAX`
      real(dp), intent(in) :: x(:)            !! Non-negative
      real(dp), intent(out) :: values(:, 0:)  !! (size(x), 0:n_max)

      real(dp), parameter :: HALF = 0.5_dp, THIRD = 1.0_dp/3.0_dp, QUARTER = 0.25_dp
      integer :: ig, n, np
      integer :: node(size(x))
      real(dp) :: d(size(x)), inv_2x(size(x)), far(size(x))

      logical :: any_near, any_far

      ! A batch is compact, so for one primitive pair its points are usually
      ! all inside the table or all beyond it, and only that branch is
      ! computed. Where both are needed, both are computed for every point and
      ! selected afterwards: no branch sits inside a loop over the points, and
      ! the arguments are clamped so the branch a point does not take still
      ! computes something finite.
      np = size(x)
      any_near = .false.
      any_far = .false.
      do ig = 1, np
         any_near = any_near .or. x(ig) < BOYS_X_MAX
         any_far = any_far .or. x(ig) >= BOYS_X_MAX
      end do

      if (any_near) then
         do ig = 1, np
            node(ig) = int(min(x(ig), BOYS_X_MAX)*BOYS_INV_STEP + HALF)
            d(ig) = BOYS_STEP*real(node(ig), dp) - x(ig)
         end do
         do n = 0, n_max
            do ig = 1, np
               ! F_n(x) = sum_k F_(n+k)(x0) (x0 - x)^k / k!, by Horner.
               values(ig, n) = table(n, node(ig)) + d(ig)*(table(n + 1, node(ig)) &
                                                           + HALF*d(ig)*(table(n + 2, node(ig)) &
                                                                         + THIRD*d(ig)*(table(n + 3, node(ig)) &
                                                                                 + QUARTER*d(ig)*table(n + 4, node(ig)))))
            end do
         end do
      end if
      if (.not. any_far) return

      ! (2n-1)!! / (2x)^n * sqrt(pi/x) / 2
      do ig = 1, np
         inv_2x(ig) = HALF/max(x(ig), BOYS_X_MAX)
         far(ig) = sqrt(HALF*PI*inv_2x(ig))
      end do
      do n = 0, n_max
         if (n > 0) then
            do ig = 1, np
               far(ig) = far(ig)*real(2*n - 1, dp)*inv_2x(ig)
            end do
         end if
         if (any_near) then
            do ig = 1, np
               if (x(ig) >= BOYS_X_MAX) values(ig, n) = far(ig)
            end do
         else
            do ig = 1, np
               values(ig, n) = far(ig)
            end do
         end if
      end do
   end subroutine boys_values

   subroutine esp_batch_contract(shells, coords, f, g, error, bounds, screen_tol, n_skipped)
      !! `G(g, nu) = sum_lambda <nu| 1/|r - r_g| |lambda> F(g, lambda)`
      !!
      !! Point `g` of `f` meets the potential integral at point `g` of `coords`
      !! and no other. Without `bounds` this is `cosx_esp_contract` to
      !! rounding. With them the same halves of shell pairs are skipped, and
      !! within a kept pair a primitive pair is dropped too when a bound on its
      !! integral times the `F` it multiplies is under the threshold, so the
      !! two then agree to the threshold and not to rounding.
      type(esp_shells_t), intent(in) :: shells
      real(dp), intent(in) :: coords(:, :)     !! (3, n_batch), Bohr
      real(dp), intent(in) :: f(:, :)          !! (n_batch, n_ao)
      real(dp), allocatable, intent(out) :: g(:, :)   !! (n_batch, n_ao)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: bounds(:, :)
         !! Schwarz bounds over the molecule's shells; absent keeps every pair
      real(dp), intent(in), optional :: screen_tol
         !! Threshold on bound times the largest `|F|` of a shell; zero keeps
         !! everything. Read only with `bounds`.
      integer, intent(out), optional :: n_skipped
         !! Shell pairs, of the lower triangle, not evaluated for this batch

      real(dp), allocatable :: fc(:, :), gc(:, :), f_max(:)
      real(dp), allocatable :: block(:, :, :), r(:, :, :), boys(:, :), pc(:, :), x(:)
      real(dp), allocatable :: points(:, :)
      real(dp) :: tol, f_pair
      logical :: screening, into_i, into_j
      integer :: n, is, js, l, nc, ns, a, b, cs, as, skipped, l_pair

      n = size(coords, 2)
      if (size(coords, 1) /= 3 .or. n < 1) then
         call error%set(ERROR_VALIDATION, "ESP batch: grid points must be (3, n) with n > 0")
         return
      end if
      if (size(f, 1) /= n .or. size(f, 2) /= shells%n_ao) then
         call error%set(ERROR_VALIDATION, "ESP batch: F must be (n_points, n_ao)")
         return
      end if

      allocate (fc(n, shells%n_cart), gc(n, shells%n_cart), f_max(shells%n_shell))
      gc = 0.0_dp

      ! F in the Cartesian functions of each shell, and the largest element
      ! each shell carries over the batch.
      do is = 1, shells%n_shell
         l = shells%l(is)
         nc = (l + 1)*(l + 2)/2
         cs = shells%cart_first(is)
         as = shells%ao_first(is)
         if (shells%cartesian) then
            fc(:, cs + 1:cs + nc) = f(:, as + 1:as + nc)
            f_max(is) = maxval(abs(f(:, as + 1:as + nc)))
         else
            ns = 2*l + 1
            f_max(is) = maxval(abs(f(:, as + 1:as + ns)))
            do b = 1, nc
               fc(:, cs + b) = 0.0_dp
               do a = 1, ns
                  if (shells%c2s(a, b, l) == 0.0_dp) cycle
                  fc(:, cs + b) = fc(:, cs + b) + shells%c2s(a, b, l)*f(:, as + a)
               end do
            end do
         end if
      end do

      screening = present(bounds)
      tol = 0.0_dp
      if (present(screen_tol)) tol = screen_tol
      skipped = 0

      l_pair = 2*shells%l_max
      allocate (block(n, N_CART_MAX, N_CART_MAX), r(n, N_HERM_MAX, 0:l_pair), &
                boys(n, 0:l_pair), pc(n, 3), x(n))
      ! One column per direction, so the loops over the points run down
      ! contiguous memory.
      points = transpose(coords)

      do is = 1, shells%n_shell
         do js = 1, is
            into_i = .true.
            into_j = is /= js
            if (screening) then
               into_i = bounds(shells%parent(is), shells%parent(js))*f_max(js) >= tol
               into_j = into_j .and. &
                        bounds(shells%parent(is), shells%parent(js))*f_max(is) >= tol
               if (.not. (into_i .or. into_j)) then
                  skipped = skipped + 1
                  cycle
               end if
            end if
            f_pair = 0.0_dp
            if (into_i) f_pair = f_max(js)
            if (into_j) f_pair = max(f_pair, f_max(is))
            call pair_contract(shells, is, js, n, points, fc, gc, into_i, into_j, &
                               merge(tol, 0.0_dp, screening), f_pair, block, r, boys, pc, x)
         end do
      end do
      if (present(n_skipped)) n_skipped = skipped

      ! Back to the molecule's own functions.
      allocate (g(n, shells%n_ao))
      do is = 1, shells%n_shell
         l = shells%l(is)
         nc = (l + 1)*(l + 2)/2
         cs = shells%cart_first(is)
         as = shells%ao_first(is)
         if (shells%cartesian) then
            g(:, as + 1:as + nc) = gc(:, cs + 1:cs + nc)
         else
            do a = 1, 2*l + 1
               g(:, as + a) = 0.0_dp
               do b = 1, nc
                  if (shells%c2s(a, b, l) == 0.0_dp) cycle
                  g(:, as + a) = g(:, as + a) + shells%c2s(a, b, l)*gc(:, cs + b)
               end do
            end do
         end if
      end do
   end subroutine esp_batch_contract

   subroutine pair_contract(shells, is, js, n, points, fc, gc, into_i, into_j, &
                            tol, f_pair, block, r, boys, pc, x)
      !! One shell pair: its Cartesian integral block over the batch, then both halves
      type(esp_shells_t), intent(in) :: shells
      integer, intent(in) :: is, js, n
      real(dp), intent(in) :: points(:, :)   !! (n, 3): one column per direction
      real(dp), intent(in) :: fc(:, :)
      real(dp), intent(inout) :: gc(:, :)
      logical, intent(in) :: into_i, into_j
      real(dp), intent(in) :: tol
         !! A primitive pair is dropped when a bound on its integral, times
         !! `f_pair`, is below this. Zero keeps every one.
      real(dp), intent(in) :: f_pair
         !! The largest `|F|` this pair's integral will multiply
      real(dp), intent(inout) :: block(:, :, :)   !! (n, N_CART_MAX, N_CART_MAX) scratch
      real(dp), intent(inout) :: r(:, :, 0:)      !! (n, N_HERM_MAX, 0:) scratch
      real(dp), intent(inout) :: boys(:, 0:)      !! (n, 0:) scratch
      real(dp), intent(inout) :: pc(:, :)         !! (n, 3) scratch
      real(dp), intent(inout) :: x(:)             !! (n) scratch

      integer :: la, lb, l_pair, nca, ncb, ka, kb, ia, ib, ig, t, u, v, h
      integer :: ax, ay, az, bx, by, bz, csa, csb
      real(dp) :: ex(0:L_PAIR_MAX, 0:ESP_BATCH_LMAX, 0:ESP_BATCH_LMAX)
      real(dp) :: ey(0:L_PAIR_MAX, 0:ESP_BATCH_LMAX, 0:ESP_BATCH_LMAX)
      real(dp) :: ez(0:L_PAIR_MAX, 0:ESP_BATCH_LMAX, 0:ESP_BATCH_LMAX)
      real(dp) :: a_centre(3), b_centre(3), p_centre(3)
      real(dp) :: ab2, alpha, beta, p, mu
      real(dp) :: prefactor, coefficient, bound
      logical :: any_primitive

      la = shells%l(is)
      lb = shells%l(js)
      l_pair = la + lb
      nca = (la + 1)*(la + 2)/2
      ncb = (lb + 1)*(lb + 2)/2
      csa = shells%cart_first(is)
      csb = shells%cart_first(js)
      a_centre = shells%centre(:, is)
      b_centre = shells%centre(:, js)
      ab2 = sum((a_centre - b_centre)**2)

      block(:, 1:nca, 1:ncb) = 0.0_dp
      any_primitive = .false.

      do ka = 1, shells%n_prim(is)
         alpha = shells%expo(shells%prim_first(is) + ka)
         do kb = 1, shells%n_prim(js)
            beta = shells%expo(shells%prim_first(js) + kb)
            p = alpha + beta
            mu = alpha*beta/p
            if (mu*ab2 > PRIM_EXPONENT_CUT) cycle
            any_primitive = .true.
            p_centre = (alpha*a_centre + beta*b_centre)/p
            prefactor = TWO_PI/p*exp(-mu*ab2) &
                        *shells%coef(shells%prim_first(is) + ka) &
                        *shells%coef(shells%prim_first(js) + kb)

            if (l_pair == 0) then
               ! Two s functions: the prefactor times F_0 and nothing else.
               if (tol > 0.0_dp) then
                  if (abs(prefactor)*f_pair < tol) cycle
               end if
               do ig = 1, n
                  x(ig) = p*((p_centre(1) - points(ig, 1))**2 + (p_centre(2) - points(ig, 2))**2 &
                             + (p_centre(3) - points(ig, 3))**2)
               end do
               call boys_values(shells%boys, 0, x(1:n), boys(1:n, 0:0))
               do ig = 1, n
                  block(ig, 1, 1) = block(ig, 1, 1) + prefactor*boys(ig, 0)
               end do
               cycle
            end if

            call hermite_expansion(la, lb, p, p_centre(1) - a_centre(1), &
                                   p_centre(1) - b_centre(1), ex)
            call hermite_expansion(la, lb, p, p_centre(2) - a_centre(2), &
                                   p_centre(2) - b_centre(2), ey)
            call hermite_expansion(la, lb, p, p_centre(3) - a_centre(3), &
                                   p_centre(3) - b_centre(3), ez)

            if (tol > 0.0_dp) then
               ! |R_tuv| is at most of order (2p)^(degree/2), and F_0 at most
               ! one, so this bounds every element of the primitive block.
               bound = abs(prefactor)*maxval(abs(ex(0:l_pair, 0:la, 0:lb))) &
                       *maxval(abs(ey(0:l_pair, 0:la, 0:lb))) &
                       *maxval(abs(ez(0:l_pair, 0:la, 0:lb))) &
                       *max(1.0_dp, 2.0_dp*p)**(0.5_dp*real(l_pair, dp))
               if (bound*f_pair < tol) then
                  cycle
               end if
            end if
            do ig = 1, n
               pc(ig, 1) = p_centre(1) - points(ig, 1)
               pc(ig, 2) = p_centre(2) - points(ig, 2)
               pc(ig, 3) = p_centre(3) - points(ig, 3)
               x(ig) = p*(pc(ig, 1)**2 + pc(ig, 2)**2 + pc(ig, 3)**2)
            end do
            call boys_values(shells%boys, l_pair, x(1:n), boys(1:n, 0:l_pair))
            call hermite_coulomb(n, l_pair, p, pc, boys, shells%h_of, r)

            do ib = 1, ncb
               bx = shells%powers(1, ib, lb)
               by = shells%powers(2, ib, lb)
               bz = shells%powers(3, ib, lb)
               do ia = 1, nca
                  ax = shells%powers(1, ia, la)
                  ay = shells%powers(2, ia, la)
                  az = shells%powers(3, ia, la)
                  do v = 0, az + bz
                     do u = 0, ay + by
                        do t = 0, ax + bx
                           coefficient = prefactor*ex(t, ax, bx)*ey(u, ay, by)*ez(v, az, bz)
                           if (coefficient == 0.0_dp) cycle
                           h = shells%h_of(t, u, v)
                           do ig = 1, n
                              block(ig, ia, ib) = block(ig, ia, ib) + coefficient*r(ig, h, 0)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (.not. any_primitive) return

      do ib = 1, ncb
         do ia = 1, nca
            if (into_i) then
               do ig = 1, n
                  gc(ig, csa + ia) = gc(ig, csa + ia) + block(ig, ia, ib)*fc(ig, csb + ib)
               end do
            end if
            if (into_j) then
               do ig = 1, n
                  gc(ig, csb + ib) = gc(ig, csb + ib) + block(ig, ia, ib)*fc(ig, csa + ia)
               end do
            end if
         end do
      end do
   end subroutine pair_contract

   pure subroutine cartesian_powers(l, powers)
      !! The exponents of x, y, z for each Cartesian function of a shell, libcint's order
      integer, intent(in) :: l
      integer, intent(out) :: powers(:, :)   !! (3, n_cart)

      integer :: i, j, comp

      comp = 0
      do i = 0, l
         do j = 0, i
            comp = comp + 1
            powers(1, comp) = l - i
            powers(2, comp) = i - j
            powers(3, comp) = j
         end do
      end do
   end subroutine cartesian_powers

   pure subroutine hermite_index(l_pair, h_of)
      !! Position of the Hermite function (t, u, v) in a flat list, lowest degree first
      integer, intent(in) :: l_pair
      integer, intent(out) :: h_of(0:, 0:, 0:)

      integer :: degree, t, u, v, h

      h_of = 0
      h = 0
      do degree = 0, l_pair
         do t = degree, 0, -1
            do u = degree - t, 0, -1
               v = degree - t - u
               h = h + 1
               h_of(t, u, v) = h
            end do
         end do
      end do
   end subroutine hermite_index

   pure subroutine hermite_expansion(la, lb, p, pa, pb, e)
      !! One Cartesian direction of the Hermite expansion of a Gaussian product
      !!
      !! `e(t, i, j)`, with the product's own prefactor left out:
      !! `e(0, 0, 0) = 1`.
      integer, intent(in) :: la, lb
      real(dp), intent(in) :: p    !! Sum of the two exponents
      real(dp), intent(in) :: pa   !! Product centre minus the first centre
      real(dp), intent(in) :: pb   !! Product centre minus the second
      real(dp), intent(out) :: e(0:, 0:, 0:)

      integer :: i, j, t
      real(dp) :: half_inv_p

      half_inv_p = 0.5_dp/p
      e(:, 0:la, 0:lb) = 0.0_dp
      e(0, 0, 0) = 1.0_dp
      do i = 1, la
         do t = 0, i
            e(t, i, 0) = pa*e(t, i - 1, 0)
            if (t > 0) e(t, i, 0) = e(t, i, 0) + half_inv_p*e(t - 1, i - 1, 0)
            if (t < i - 1) e(t, i, 0) = e(t, i, 0) + real(t + 1, dp)*e(t + 1, i - 1, 0)
         end do
      end do
      do j = 1, lb
         do i = 0, la
            do t = 0, i + j
               e(t, i, j) = pb*e(t, i, j - 1)
               if (t > 0) e(t, i, j) = e(t, i, j) + half_inv_p*e(t - 1, i, j - 1)
               if (t < i + j - 1) e(t, i, j) = e(t, i, j) + real(t + 1, dp)*e(t + 1, i, j - 1)
            end do
         end do
      end do
   end subroutine hermite_expansion

   pure subroutine hermite_coulomb(n, l_pair, p, pc, boys, h_of, r)
      !! The Hermite Coulomb integrals `R_tuv` at every point of the batch
      !!
      !! `r(:, h, 0)` on return is `R_tuv` for `h = h_of(t, u, v)`; the higher
      !! auxiliary orders in the last index are scratch.
      integer, intent(in) :: n, l_pair
      real(dp), intent(in) :: p
      real(dp), intent(in) :: pc(:, :)        !! (n, 3) product centre minus the point
      real(dp), intent(in) :: boys(:, 0:)     !! (n, 0:l_pair)
      integer, intent(in) :: h_of(0:, 0:, 0:)
      real(dp), intent(inout) :: r(:, :, 0:)  !! (n, n_herm, 0:l_pair)

      integer :: m, degree, t, u, v, h, h1, h2, ig, axis, lower
      real(dp) :: scale, weight

      ! R^m_000 = (-2p)^m F_m
      scale = 1.0_dp
      do m = 0, l_pair
         do ig = 1, n
            r(ig, 1, m) = scale*boys(ig, m)
         end do
         scale = -2.0_dp*p*scale
      end do

      do degree = 1, l_pair
         do t = degree, 0, -1
            do u = degree - t, 0, -1
               v = degree - t - u
               h = h_of(t, u, v)
               ! Step down along the first direction with a nonzero index.
               if (t > 0) then
                  axis = 1
                  lower = t - 1
                  h1 = h_of(t - 1, u, v)
                  h2 = 0
                  if (t > 1) h2 = h_of(t - 2, u, v)
               else if (u > 0) then
                  axis = 2
                  lower = u - 1
                  h1 = h_of(t, u - 1, v)
                  h2 = 0
                  if (u > 1) h2 = h_of(t, u - 2, v)
               else
                  axis = 3
                  lower = v - 1
                  h1 = h_of(t, u, v - 1)
                  h2 = 0
                  if (v > 1) h2 = h_of(t, u, v - 2)
               end if
               weight = real(lower, dp)
               do m = 0, l_pair - degree
                  if (h2 > 0) then
                     do ig = 1, n
                        r(ig, h, m) = weight*r(ig, h2, m + 1) + pc(ig, axis)*r(ig, h1, m + 1)
                     end do
                  else
                     do ig = 1, n
                        r(ig, h, m) = pc(ig, axis)*r(ig, h1, m + 1)
                     end do
                  end if
               end do
            end do
         end do
      end do
   end subroutine hermite_coulomb

end module mqc_czt_esp_batch
