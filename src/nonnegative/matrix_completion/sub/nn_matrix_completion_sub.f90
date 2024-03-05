submodule (nn_matrix_completion_mod) nn_matrix_completion_sub
implicit none (type, external)

contains
!----------------------------------------------------------------------------------------------------------------------------

module subroutine rgd_completion_step(m, n, r, S, U, ldu, VT, ldvt, &
        sample_size, sample, sample_indices, &
        work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP => DP
use maria_constants_mod,only:   &
    ONE => D_ONE,               &
    ZERO => D_ZERO
use maria_comparison_mod,   only:   &
    safe_eq
use maria_la_core_mod,  only:   &
    MM => dmatmul,              &
    droundup_lwork,             &
    ddot, ddgmm, dlacpy, dgemv
use maria_lr_la_mod,    only:   &
    dlrval
use maria_lr_geom_mod,  only:   &
    dlrproj_tangent,            &
    dlrdotf_tangent,            &
    dlrretr_tangent
use maria_reports_mod,          only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    integer,    intent(in   )               :: r
    real(WP),   intent(inout),  contiguous  :: S(:)
    real(WP),   intent(inout),  contiguous  :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout),  contiguous  :: VT(:)
    integer,    intent(in   )               :: ldvt
    integer,    intent(in   )               :: sample_size
    real(WP),   intent(in   ),  contiguous  :: sample(:)
    integer,    intent(in   ),  contiguous  :: sample_indices(:)
    real(WP),   intent(  out),  contiguous  :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out),  contiguous  :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr

!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'RGD_COMPLETION_STEP'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(16)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i, mem_U, mem_x, ifoo(1), &
                    mem_grad, i_grad, i_wrk, lwrk, mem_pU, i_pU, mem_pVT, i_pVT,  &
                    mem_C, i_C
    real(WP)    ::  foo(1), step, rgrad_nrm
    procedure(MM), pointer  :: B2BA, B2AB

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1)  = (m < 0)
    bad_arg(2)  = (n < 0)
    bad_arg(3)  = (r < 0) .or. (r > min(m,n))
    bad_arg(6)  = (ldu < max(1,m))
    bad_arg(8)  = (ldvt < max(1,r))
    bad_arg(9)  = (sample_size < 0) .or. (sample_size > m*n)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]:
    mem_U = m*r
    mem_x = r
    mem_grad = sample_size
    mem_pU = m*r
    mem_pVT = r*n
    mem_C = r*r
    
    ! Query procedures
    call dlrretr_tangent(m, n, r, foo, foo, ldu, foo, ldvt, ONE, &
        foo, m, foo, r, foo, r, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)

    ! Total
    req_lw = mem_pu + mem_pvt + mem_c + max(mem_U + mem_grad, lw_proc)
    req_liw = liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(13) = (lwork < req_lw)
    bad_arg(15) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------
i_pU = 1
i_pVT = i_pU + mem_pu
i_c = i_pvt + mem_pvt
i_grad = i_c + mem_c
i_wrk = i_grad + mem_grad
lwrk = lwork - i_wrk + 1
associate(pu => work(i_pu:), pvt => work(i_pvt:), c => work(i_c:),  gradient => work(i_grad:), wrk => work(i_wrk:))
    ! Form Euclidean gradient
    call dlacpy('a', m, r, U, ldu, wrk, m)
    call ddgmm('r', m, r, wrk, m, S, 1, ierr)
    do i = 1, sample_size
    associate( irow => sample_indices(2*i-1), icol => sample_indices(2*i) )
        gradient(i) = dlrval(m, n, irow, icol, r, wrk, m, VT, ldvt, ierr)
    end associate
        gradient(i) = gradient(i) - sample(i)
    end do

    ! Compute Riemannian gradient
    B2AB => proj_grad_r
    B2BA => proj_grad_l
    call dlrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, m, pVT, r, C, r, ierr)
end associate

associate(pu => work(i_pu:), pvt => work(i_pvt:), c => work(i_c:), sample_rgrad => work(i_grad:), wrk => work(i_wrk:))
    ! Evaluate Riemannian gradient on the sampling set
    do i = 1, sample_size
    associate( irow => sample_indices(2*i-1), icol => sample_indices(2*i) )
        sample_rgrad(i) = dlrval(m, n, irow, icol, r, U, ldu, pVT, r, ierr)
        sample_rgrad(i) = sample_rgrad(i) &
            + dlrval(m, n, irow, icol, r, pU, m, VT, ldvt, ierr)
        associate( row => U(irow:), col => VT(1 + (icol-1)*ldvt:) )
            call dgemv('n', r, r, ONE, C, r, col, 1, ZERO, wrk, 1)
            sample_rgrad(i) = sample_rgrad(i) +  ddot(r, row, 1, wrk, 1)
        end associate
    end associate
    end do

    ! Compute the step of steepest descent
    step = sqrt(dlrdotf_tangent(m, n, r, pU, m, pVT, r, C, r, pU, m, pVT, r, C, r, ierr))
    rgrad_nrm = sqrt(ddot(sample_size, sample_rgrad, 1, sample_rgrad, 1))
    if (safe_eq(rgrad_nrm, ZERO)) return
    step = (step / rgrad_nrm)**2
end associate

i_wrk = i_c + mem_C
lwrk = lwork - i_wrk + 1
associate(pu => work(i_pu:), pvt => work(i_pvt:), c => work(i_c:), wrk => work(i_wrk:))
    ! Do a step and retract
    call dlrretr_tangent(m, n, r, S, U, ldu, VT, ldvt, -step, pU, m, pVT, r, C, r, wrk, lwrk, iwork, liwork, ierr)
end associate
contains
    subroutine proj_grad_r &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgescal, dgemm
    implicit none
        character(1), intent(in)                :: transB
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: B(:)
        integer,      intent(in)                :: ldB
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: ierr
    
        integer :: i, j

!        associate(grad => work(i_grad:))
!            call dgemm('n', transB, m, n, k, alpha, grad, m, B, ldb, beta, C, ldc)
!        end associate
!        return

        call dgescal(m, n, beta, C, ldc, ierr)
        do i = 1, sample_size
        associate( irow => sample_indices(2*i-1), icol => sample_indices(2*i), grad => work(i_grad:) )
            do j = 1, n
                C(irow + (j-1)*ldc) = C(irow + (j-1)*ldc) &
                    + alpha * grad(i) * B(j + (icol-1)*ldb)
            end do
        end associate
        end do
    end subroutine proj_grad_r

    subroutine proj_grad_l &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgescal, dgemm
    implicit none
        character(1), intent(in)                :: transB
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: B(:)
        integer,      intent(in)                :: ldB
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: ierr
    
        integer :: i, j

!        associate(grad => work(i_grad:))
!            call dgemm(transB, 'n', m, n, k, alpha, B, ldb, grad, k, beta, C, ldc)
!        end associate
!        return

        call dgescal(m, n, beta, C, ldc, ierr)
        do i = 1, sample_size
        associate( irow => sample_indices(2*i-1), icol => sample_indices(2*i), grad => work(i_grad:) )
            do j = 1, m
                C(j + (icol-1)*ldc) = C(j + (icol-1)*ldc) &
                    + alpha * grad(i) * B(irow + (j-1)*ldb)
            end do
        end associate
        end do
    end subroutine proj_grad_l
end subroutine rgd_completion_step

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol &
(method, m, n, r, U, ldu, VT, ldvt, irow, icol, niter, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, neg_pos)
use maria_kinds_mod,            only:   &
    WP      => DP
use maria_constants_mod,        only:   &
    ZERO    => D_ZERO,                  &
    ONE     => D_ONE
use maria_comparison_mod,       only:   &
    safe_less
use maria_la_core_mod,          only:   &
    dgemv,                              &
    dlacpy,                             &
    ddgmm,                              &
    droundup_lwork
use maria_access_matrix_mod,    only:   &
    MS      => dmatslc
use maria_lr_la_mod,            only:   &
    dlrval,                             &
    dlrort,                             &
    dlrsvd_ort
use maria_lr_maxvol_mod,        only:   &
    dgemaxvol
use maria_lr_cross_mod,         only:   &
    dmatcross
use maria_reports_mod,          only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    integer,    intent(in   )               :: r
    real(WP),   intent(inout), contiguous   :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout), contiguous   :: VT(:)
    integer,    intent(in   )               :: ldvt
    integer,    intent(inout), contiguous   :: irow(:)
    integer,    intent(inout), contiguous   :: icol(:)
    integer,    intent(in   )               :: niter
    real(WP),   intent(in   )               :: thresh
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
    integer,    intent(inout), optional     :: neg_pos(2)
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_MAXVOL'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(19)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_UU,    &
                    i_UU, mem_VVT, i_VVT, ifoo(1)
    real(WP)    ::  foo(1), lowest, ratio

    procedure(MS), pointer :: fun

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1)  = (method < 1) .or. (method > 4)
    bad_arg(2)  = (m < 0)
    bad_arg(3)  = (n < 0)
    bad_arg(4)  = (r < 0) .or. (r > min(m,n))
    bad_arg(6)  = (ldu < max(1,m))
    bad_arg(8)  = (ldvt < max(1,r))
    bad_arg(11) = (niter < 0)
    bad_arg(12) = safe_less(thresh, ZERO)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: UU, VVT
    mem_UU = m * r
    mem_VVT = r * n

    ! Query procedures
    call dgemaxvol(m, n, fun, r, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dmatcross(m, n, fun, r, ifoo, r, ifoo, r, foo, m, foo, r, foo, -1, ifoo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    liw_proc = max(liw_proc, ifoo(1))

    ! Total
    req_lw = mem_UU + mem_VVT + lw_proc
    req_liw = liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(14) = (lwork < req_lw)
    bad_arg(16) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..UU..|..VVT..|..wrk..|
i_UU = 1
i_VVT = i_UU + mem_UU
i_wrk = i_VVT + mem_VVT
lwrk = lwork - i_wrk + 1
associate (UU => work(i_UU:), VVT => work(i_VVT:), wrk => work(i_wrk:))
    lowest = ZERO
    if (method >= 3) then
        fun => dlrslc_min
        call dgemaxvol(m, n, fun, 1, neg_pos(1:), neg_pos(2:), niter, ONE, wrk, lwrk, iwork, liwork, ierr)
        lowest = dlrval(m, n, neg_pos(1), neg_pos(2), r, U, ldu, VT, ldvt, ierr)
        lowest = min(lowest, ZERO)
        lowest = -lowest * ratio
    end if    

    ! Low-rank and nonnegative steps combined
    fun => dlrslc_nn
    call dgemaxvol(m, n, fun, r, irow, icol, niter, thresh, wrk, lwrk, iwork, liwork, ierr)
    if (ierr > 0) return
    call dmatcross(m, n, fun, r, icol, r, irow, r, UU, m, VVT, r, wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, UU, m, U, ldu)
    call dlacpy('a', r, n, VVT, r, VT, ldvt)
end associate

contains
    function dlrval_min(m, n, i, j, ierr)
    use maria_constants_mod,    only:   &
        PI => D_PI
    use maria_lr_la_mod,    only:   &
        dlrval
        integer,    intent(in   )   :: m
        integer,    intent(in   )   :: n
        integer,    intent(in   )   :: i
        integer,    intent(in   )   :: j
        integer,    intent(  out)   :: ierr
        real(WP)                    :: dlrval_min

        dlrval_min = dlrval(m, n, i, j, r, U, ldu, VT, ldvt, ierr)
        dlrval_min = PI / 2 - atan(dlrval_min - lowest)
    end function dlrval_min

    subroutine dlrslc_min(m, n, mode, ind, x, incx, ierr)
    use maria_access_matrix_mod,    only:   &
        dmatval2slc,                        &
        MV  => dmatval
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        integer,    intent(in   )               :: mode
        integer,    intent(in   )               :: ind
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr

        procedure(MV),  pointer :: fun

        fun => dlrval_min
        call dmatval2slc(fun, m, n, mode, ind, x, incx, ierr)
    end subroutine dlrslc_min

    function dlrval_nn(m, n, i, j, ierr)
    use maria_lr_la_mod,    only:   &
        dlrval
        integer,    intent(in   )   :: m
        integer,    intent(in   )   :: n
        integer,    intent(in   )   :: i
        integer,    intent(in   )   :: j
        integer,    intent(  out)   :: ierr
        real(WP)                    :: dlrval_nn

        dlrval_nn = dlrval(m, n, i, j, r, U, ldu, VT, ldvt, ierr)
            if (method == 1 .or. method == 3) then
                ! Projection
                dlrval_nn = max(lowest, dlrval_nn)
            elseif (method == 2) then
                ! Reflection
                dlrval_nn = abs(dlrval_nn)
            else
                dlrval_nn = dlrval_nn + lowest
            end if
    end function dlrval_nn

    subroutine dlrslc_nn(m, n, mode, ind, x, incx, ierr)
    use maria_access_matrix_mod,    only:   &
        dmatval2slc,                        &
        MV  => dmatval
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        integer,    intent(in   )               :: mode
        integer,    intent(in   )               :: ind
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr

        procedure(MV),  pointer :: fun

        fun => dlrval_nn
        call dmatval2slc(fun, m, n, mode, ind, x, incx, ierr)
    end subroutine dlrslc_nn
end subroutine nn_lrma_maxvol

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol_proj &
(method, m, n, r, U, ldu, VT, ldvt, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, &
    work, lwork, iwork, liwork, ierr, shift_ratio, neg_pos)
use maria_kinds_mod,            only:   &
    WP      => DP
use maria_constants_mod,        only:   &
    ZERO    => D_ZERO,                  &
    ONE     => D_ONE
use maria_comparison_mod,       only:   &
    safe_less
use maria_la_core_mod,          only:   &
    dgemv,                              &
    dlacpy,                             &
    droundup_lwork
use maria_access_matrix_mod,    only:   &
    MS      => dmatslc
use maria_lr_la_mod,            only:   &
    dlrval
use maria_lr_maxvol_mod,        only:   &
    dgemaxvol_proj,                     &
    dgemaxvol
use maria_lr_cross_mod,         only:   &
    dmatcross
use maria_reports_mod,          only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: method                       !  1
    integer,    intent(in   )               :: m                            !  2
    integer,    intent(in   )               :: n                            !  3
    integer,    intent(in   )               :: r                            !  4
    real(WP),   intent(inout), contiguous   :: U(:)                         !  5
    integer,    intent(in   )               :: ldu                          !  6
    real(WP),   intent(inout), contiguous   :: VT(:)                        !  7
    integer,    intent(in   )               :: ldvt                         !  8
    integer,    intent(in   )               :: kr                           !  9
    integer,    intent(inout), contiguous   :: irow(:)                      ! 10
    integer,    intent(inout), contiguous   :: icol_short(:)                ! 11
    integer,    intent(in   )               :: kc                           ! 12
    integer,    intent(inout), contiguous   :: icol(:)                      ! 13
    integer,    intent(inout), contiguous   :: irow_short(:)                ! 14
    integer,    intent(in   )               :: niter                        ! 15
    real(WP),   intent(in   )               :: thresh                       ! 16
    real(WP),   intent(  out), contiguous   :: work(:)                      ! 17
    integer,    intent(in   )               :: lwork                        ! 18
    integer,    intent(  out), contiguous   :: iwork(:)                     ! 19
    integer,    intent(in   )               :: liwork                       ! 20
    integer,    intent(  out)               :: ierr                         ! 21
    real(WP),   intent(in   ), optional     :: shift_ratio                  ! 22
    integer,    intent(inout), optional     :: neg_pos(2)

!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_MAXVOL_PROJ'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(21)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_UU,    &
                    i_UU, mem_VVT, i_VVT, ifoo(1)
    real(WP)    ::  foo(1), lowest, ratio

    procedure(MS), pointer :: fun

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1)  = (method < 1) .or. (method > 4)
    bad_arg(2)  = (m < 0)
    bad_arg(3)  = (n < 0)
    bad_arg(4)  = (r < 0) .or. (r > min(m,n))
    bad_arg(6)  = (ldu < max(1,m))
    bad_arg(8)  = (ldvt < max(1,r))
    bad_arg(11) = (niter < 0)
    bad_arg(12) = safe_less(thresh, ZERO)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: UU, VVT
    mem_UU = m * r
    mem_VVT = r * n
    ! Query procedures
    call dgemaxvol_proj(m, n, fun, r, kr, ifoo, ifoo, kc, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dmatcross(m, n, fun, kc, ifoo, kr, ifoo, r, foo, m, foo, r, foo, -1, ifoo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    liw_proc = max(liw_proc, ifoo(1))
    ! Total
    req_lw = mem_UU + mem_VVT + lw_proc
    req_liw = liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(14) = (lwork < req_lw)
    bad_arg(16) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..UU..|..VVT..|..wrk..|
i_UU = 1
i_VVT = i_UU + mem_UU
i_wrk = i_VVT + mem_VVT
lwrk = lwork - i_wrk + 1
associate (UU => work(i_UU:), VVT => work(i_VVT:), wrk => work(i_wrk:))   
    lowest = ZERO
    if (method >= 3) then
        fun => dlrslc_min
        call dgemaxvol(m, n, fun, 1, neg_pos(1:), neg_pos(2:), niter, ONE, wrk, lwrk, iwork, liwork, ierr)
        lowest = dlrval(m, n, neg_pos(1), neg_pos(2), r, U, ldu, VT, ldvt, ierr)
        lowest = min(lowest, ZERO)
        lowest = -lowest * ratio
    end if    

    ! Low-rank and nonnegative steps combined
    fun => dlrslc_nn
    call dgemaxvol_proj(m, n, fun, r, kr, irow, icol_short, kc, icol, irow_short, &
        niter, thresh, wrk, lwrk, iwork, liwork, ierr)
    if (ierr > 0) return
    call dmatcross(m, n, fun, kc, icol, kr, irow, r, UU, m, VVT, r, wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, UU, m, U, ldu)
    call dlacpy('a', r, n, VVT, r, VT, ldvt)
end associate

contains
    function dlrval_min(m, n, i, j, ierr)
    use maria_constants_mod,    only:   &
        PI => D_PI
    use maria_lr_la_mod,    only:   &
        dlrval
        integer,    intent(in   )   :: m
        integer,    intent(in   )   :: n
        integer,    intent(in   )   :: i
        integer,    intent(in   )   :: j
        integer,    intent(  out)   :: ierr
        real(WP)                    :: dlrval_min

        dlrval_min = dlrval(m, n, i, j, r, U, ldu, VT, ldvt, ierr)
        dlrval_min = PI / 2 - atan(dlrval_min - lowest)
    end function dlrval_min

    subroutine dlrslc_min(m, n, mode, ind, x, incx, ierr)
    use maria_access_matrix_mod,    only:   &
        dmatval2slc,                        &
        MV  => dmatval
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        integer,    intent(in   )               :: mode
        integer,    intent(in   )               :: ind
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr

        procedure(MV),  pointer :: fun

        fun => dlrval_min
        call dmatval2slc(fun, m, n, mode, ind, x, incx, ierr)
    end subroutine dlrslc_min

    function dlrval_nn(m, n, i, j, ierr)
    use maria_lr_la_mod,    only:   &
        dlrval
        integer,    intent(in   )   :: m
        integer,    intent(in   )   :: n
        integer,    intent(in   )   :: i
        integer,    intent(in   )   :: j
        integer,    intent(  out)   :: ierr
        real(WP)                    :: dlrval_nn

        dlrval_nn = dlrval(m, n, i, j, r, U, ldu, VT, ldvt, ierr)
            if (method == 1 .or. method == 3) then
                ! Projection
                dlrval_nn = max(lowest, dlrval_nn)
            elseif (method == 2) then
                ! Reflection
                dlrval_nn = abs(dlrval_nn)
            else
                dlrval_nn = dlrval_nn + lowest
            end if
    end function dlrval_nn

    subroutine dlrslc_nn(m, n, mode, ind, x, incx, ierr)
    use maria_access_matrix_mod,    only:   &
        dmatval2slc,                        &
        MV  => dmatval
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        integer,    intent(in   )               :: mode
        integer,    intent(in   )               :: ind
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr

        procedure(MV),  pointer :: fun

        fun => dlrval_nn
        call dmatval2slc(fun, m, n, mode, ind, x, incx, ierr)
    end subroutine dlrslc_nn
end subroutine nn_lrma_maxvol_proj

!----------------------------------------------------------------------------------------------------------------------------
end submodule nn_matrix_completion_sub
