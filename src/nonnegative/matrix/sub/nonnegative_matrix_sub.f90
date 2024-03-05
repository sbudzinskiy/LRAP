submodule (nonnegative_matrix_mod) nonnegative_matrix_sub
implicit none (type, external)

contains
!----------------------------------------------------------------------------------------------------------------------------

module subroutine distance_to_nn &
(m, n, A, lda, distc, distf)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_comparison_mod,   only:   &
    safe_less
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ),  contiguous  :: A(:)
    integer,    intent(in   )               :: lda
    real(WP),   intent(  out)               :: distc
    real(WP),   intent(  out)               :: distf

!-- Variables ------------------------------------------------------------------
    integer     :: i, j

!-- Executable section ---------------------------------------------------------
    distc = ZERO
    distf = ZERO

    do j = 1, n
        do i = 1, m
            associate (val => A(i + (j-1)*lda))
                if (safe_less(val, ZERO)) then
                    distf = distf + val**2
                    distc = max(distc, -val)
                end if
            end associate
        end do
    end do

    distf = sqrt(distf)
end subroutine distance_to_nn

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_svd &
(method, m, n, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, ierr, shift_ratio)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_comparison_mod,   only:   &
    safe_less
use maria_la_core_mod,      only:   &
    dgemm,                          &
    dlacpy,                         &
    dgesdd_q,                       &
    ddgmm,                          &
    droundup_lwork
use maria_reports_mod,      only:   &
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
    real(WP),   intent(  out), contiguous   :: work(:)                      !  9
    integer,    intent(in   )               :: lwork                        ! 10
    integer,    intent(  out), contiguous   :: iwork(:)                     ! 11
    integer,    intent(in   )               :: liwork                       ! 12
    integer,    intent(  out)               :: ierr                         ! 13
    real(WP),   intent(in   ), optional     :: shift_ratio                  ! 14
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_SVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(14)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_A,     &
                    i_A, mem_bigU, i_bigU, mem_bigVT, i_bigVT, mem_bigS,        &
                    i_bigS, i, j, ifoo(1)
    real(WP)    ::  foo(1), lowest, ratio

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1) = (method < 1) .or. (method > 4)
    bad_arg(2) = (m < 0)
    bad_arg(3) = (n < 0)
    bad_arg(4) = (r < 0) .or. (r > min(m,n))
    bad_arg(6) = (ldu < max(1,m))
    bad_arg(8) = (ldvt < max(1,r))

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: A, bigU, bigS, bigVT
    mem_A = m * n
    mem_bigU = m * min(m,n)
    mem_bigS = min(m,n)
    mem_bigVT = min(m,n) * n
    ! Query procedures
    call dgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    ! Total
    req_lw = mem_A + mem_bigU + mem_bigS + mem_bigVT + lw_proc
    req_liw = liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(10) = (lwork < req_lw)
    bad_arg(12) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..A..|..bigS..|..bigU..|..bigVT..|..wrk..|
i_A = 1
i_bigS = i_A + mem_A
i_bigU = i_bigS + mem_bigS
i_bigVT = i_bigU + mem_bigU
i_wrk = i_bigVT + mem_bigVT
lwrk = lwork - i_wrk + 1
associate (A => work(i_A:), bigS => work(i_bigS:), bigU => work(i_bigU:),       &
bigVT => work(i_bigVT:), wrk => work(i_wrk:))
    call dgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)
    
    if (method == 1) then
        lowest = ZERO
    else
        ! Choose shift as the absolute value of the largest negative entry
        call distance_to_nn(m, n, A, m, lowest, foo(1))
        lowest = lowest * ratio
    end if

    ! Nonnegative step
    do j = 1, n
        do i = 1, m
            if (method == 1 .or. method == 3) then
                ! Projection
                A(i + (j-1)*m) = max(lowest, A(i + (j-1)*m))
            elseif (method == 2) then
                ! Reflection
                A(i + (j-1)*m) = abs(A(i + (j-1)*m))
            else
                A(i + (j-1)*m) = A(i + (j-1)*m) + lowest
            end if
        end do
    end do

    ! Low-rank step
    call dgesdd_q('s', m, n, A, m, bigS, bigU, m, bigVT, min(m,n), wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, bigU, m, U, ldu)
    call ddgmm('r', m, r, U, ldu, bigS, 1, ierr)
    call dlacpy('a', r, n, bigVT, min(m,n), VT, ldvt)
end associate
end subroutine nn_lrma_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_rsvd &
(method, m, n, r, U, ldu, VT, ldvt, rng, p, q, work, lwork, iwork, liwork, ierr, shift_ratio)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_prng_mod,     only:   &
    prng
use maria_comparison_mod,   only:   &
    safe_less
use maria_la_core_mod,      only:   &
    dgemm,                          &
    dlacpy,                         &
    ddgmm,                          &
    droundup_lwork,                 &
    MM => dmatmul
use maria_lr_tsvd_mod,      only:   &
    dgersvd1
use maria_reports_mod,      only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,        intent(in   )               :: method                   !  1
    integer,        intent(in   )               :: m                        !  2
    integer,        intent(in   )               :: n                        !  3
    integer,        intent(in   )               :: r                        !  4
    real(WP),       intent(inout), contiguous   :: U(:)                     !  5
    integer,        intent(in   )               :: ldu                      !  6
    real(WP),       intent(inout), contiguous   :: VT(:)                    !  7
    integer,        intent(in   )               :: ldvt                     !  8
    class(prng),    intent(in   )               :: rng                      !  9
    integer,        intent(in   )               :: p                        ! 10
    integer,        intent(in   )               :: q                        ! 11
    real(WP),       intent(  out), contiguous   :: work(:)                  ! 12
    integer,        intent(in   )               :: lwork                    ! 13
    integer,        intent(  out), contiguous   :: iwork(:)                 ! 14
    integer,        intent(in   )               :: liwork                   ! 15
    integer,        intent(  out)               :: ierr                     ! 16
    real(WP),       intent(in   ), optional     :: shift_ratio              ! 17
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_RSVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(17)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, lda, mem_A,&
                    i_A, mem_bigU, i_bigU, mem_bigVT, i_bigVT, mem_bigS,        &
                    i_bigS, mem_gauss, i_gauss, mem_csk, i_csk, i, j, ifoo(1)
    real(WP)    ::  foo(1), lowest, ratio
    procedure(MM), pointer :: B2AB, B2BA

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1) = (method < 1) .or. (method > 4)
    bad_arg(2) = (m < 0)
    bad_arg(3) = (n < 0)
    bad_arg(4) = (r < 0) .or. (r > min(m,n))
    bad_arg(6) = (ldu < max(1,m))
    bad_arg(8) = (ldvt < max(1,r))
    bad_arg(10) = (p < 0)
    bad_arg(11) = (q < 0)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: A, gauss, csk, bigU, bigS, bigVT
    lda = m
    mem_A = lda * n
    mem_gauss = n * (r+p)
    mem_csk = m * (r+p)
    mem_bigU = m * (r+p)
    mem_bigS = r+p
    mem_bigVT = (r+p) * n
    ! Query procedures
    call dgersvd1(m, n, B2AB, B2BA, r+p, foo, m, q, foo, foo, m, foo, r+p, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    ! Total
    req_lw = mem_A + mem_gauss + mem_csk + mem_bigU + mem_bigS + mem_bigVT + lw_proc
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

B2AB => mm_A_left
B2BA => mm_A_right

! Slice work:   |..A..|..gauss..|..csk..|..bigS..|..bigU..|..bigVT..|..wrk..|
i_A = 1
i_gauss = i_A + mem_A
i_csk = i_gauss + mem_gauss
i_bigS = i_csk + mem_csk
i_bigU = i_bigS + mem_bigS
i_bigVT = i_bigU + mem_bigU
i_wrk = i_bigVT + mem_bigVT
lwrk = lwork - i_wrk + 1
associate (A => work(i_A:), gauss => work(i_gauss:), csk => work(i_csk:), &
        bigS => work(i_bigS:), bigU => work(i_bigU:), bigVT => work(i_bigVT:), wrk => work(i_wrk:))
    call dgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, lda)
    
    if (method == 1) then
        lowest = ZERO
    else
        ! Choose shift as the absolute value of the largest negative entry
        call distance_to_nn(m, n, A, lda, lowest, foo(1))
        lowest = lowest * ratio
    end if

    ! Nonnegative step
    do j = 1, n
        do i = 1, m
            if (method == 1 .or. method == 3) then
                ! Projection
                A(i + (j-1)*lda) = max(lowest, A(i + (j-1)*lda))
            elseif (method == 2) then
                ! Reflection
                A(i + (j-1)*lda) = abs(A(i + (j-1)*lda))
            else
                A(i + (j-1)*lda) = A(i + (j-1)*lda) + lowest
            end if
        end do
    end do

    ! Low-rank step
    call rng%dnormal(n*(r+p), gauss, ZERO, ONE, ierr)
    call dgemm('n', 'n', m, r+p, n, ONE, A, lda, gauss, n, ZERO, csk, m)
    call dgersvd1(m, n, B2AB, B2BA, r+p, csk, m, q, bigS, bigU, m, bigVT, r+p, wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, bigU, m, U, ldu)
    call ddgmm('r', m, r, U, ldu, bigS, 1, ierr)
    call dlacpy('a', r, n, bigVT, r+p, VT, ldvt)
end associate
contains
    subroutine mm_A_left &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgemm
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

    associate( A => work(i_A:) )
        call dgemm('n', transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        ierr = 0
    end associate
    end subroutine mm_A_left

    subroutine mm_A_right &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgemm
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

    associate( A => work(i_A:) )
        call dgemm(transB, 'n', m, n, k, alpha, B, ldb, A, lda, beta, C, ldc)
        ierr = 0
    end associate
    end subroutine mm_A_right
end subroutine nn_lrma_rsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol &
(method, m, n, r, U, ldu, VT, ldvt, p, irow, icol, niter, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, i1, j1)
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
    integer,    intent(in   )               :: method                       !  1
    integer,    intent(in   )               :: m                            !  2
    integer,    intent(in   )               :: n                            !  3
    integer,    intent(in   )               :: r                            !  4
    real(WP),   intent(inout), contiguous   :: U(:)                         !  5
    integer,    intent(in   )               :: ldu                          !  6
    real(WP),   intent(inout), contiguous   :: VT(:)                        !  7
    integer,    intent(in   )               :: ldvt                         !  8
    integer,    intent(in   )               :: p                            !  9
    integer,    intent(inout), contiguous   :: irow(:)                      ! 10
    integer,    intent(inout), contiguous   :: icol(:)                      ! 11
    integer,    intent(in   )               :: niter                        ! 12
    real(WP),   intent(in   )               :: thresh                       ! 13
    real(WP),   intent(  out), contiguous   :: work(:)                      ! 14
    integer,    intent(in   )               :: lwork                        ! 15
    integer,    intent(  out), contiguous   :: iwork(:)                     ! 16
    integer,    intent(in   )               :: liwork                       ! 17
    integer,    intent(  out)               :: ierr                         ! 18
    real(WP),   intent(in   ), optional     :: shift_ratio                  ! 19
    integer,    intent(in   ), optional     :: i1
    integer,    intent(in   ), optional     :: j1
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_MAXVOL'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(19)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_UU,    &
                    i_UU, mem_VVT, i_VVT, ifoo(1), max_pos(2),  &
                    mem_tau, i_tau, mem_SS, i_SS,       &
                    mem_UUU, i_UUU, mem_VVVT, i_VVVT
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
    bad_arg(9)  = (p < 0)
    bad_arg(12) = (niter < 0)
    bad_arg(13) = safe_less(thresh, ZERO)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: UU, VVT, tau, SS, UUU, VVVT
    mem_UU = m * (r+p)
    mem_VVT = (r+p) * n
    if (p > 0) then
        mem_tau = r+p
        mem_SS = r+p
        mem_UUU = mem_UU
        mem_VVVT = mem_VVT
    else
        mem_tau = 0
        mem_SS = 0
        mem_UUU = 0
        mem_VVVT = 0
    end if
    ! Query procedures
    call dgemaxvol(m, n, fun, r+p, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dmatcross(m, n, fun, r+p, ifoo, r+p, ifoo, r+p, foo, m, foo, r+p, foo, -1, ifoo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    liw_proc = max(liw_proc, ifoo(1))
    if (p > 0) then
        call dlrort('l', m, n, r+p, foo, m, foo, r+p, foo, foo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        call dlrsvd_ort('s', 'l', m, n, r+p, foo, m, foo, r+p, foo, foo, foo, m, foo, r+p, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
    end if

    ! Total
    req_lw = mem_UU + mem_VVT + lw_proc + mem_tau + mem_SS + mem_UUU + mem_VVVT
    req_liw = liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(15) = (lwork < req_lw)
    bad_arg(17) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..UU..|..VVT..|..tau..|..SS..|..UUU..|..VVVT..|..wrk..|
i_UU = 1
i_VVT = i_UU + mem_UU
i_tau = i_VVT + mem_VVT
i_SS = i_tau + mem_tau
i_UUU = i_SS + mem_SS
i_VVVT = i_UUU + mem_UUU
i_wrk = i_VVVT + mem_VVVT
lwrk = lwork - i_wrk + 1
associate (UU => work(i_UU:), VVT => work(i_VVT:), tau => work(i_tau:), SS => work(i_SS:), &
        UUU => work(i_UUU:), VVVT => work(i_VVVT:), wrk => work(i_wrk:))
    lowest = ZERO
    if (method >= 3) then
        fun => dlrslc_min
        if (present(i1)) then
            max_pos(1) = i1
        else
            max_pos(1) = 1
        end if
        if (present(j1)) then
            max_pos(2) = j1
        else
            max_pos(2) = 1
        end if
        call dgemaxvol(m, n, fun, 1, max_pos(1:), max_pos(2:), niter, ONE, wrk, lwrk, iwork, liwork, ierr)
        lowest = dlrval(m, n, max_pos(1), max_pos(2), r, U, ldu, VT, ldvt, ierr)
        lowest = min(lowest, ZERO)
        lowest = -lowest * ratio
    end if    

    ! Low-rank and nonnegative steps combined
    fun => dlrslc_nn
    call dgemaxvol(m, n, fun, r+p, irow, icol, niter, thresh, wrk, lwrk, iwork, liwork, ierr)
    if (ierr > 0) return
    call dmatcross(m, n, fun, r+p, icol, r+p, irow, r+p, UU, m, VVT, r+p, wrk, lwrk, iwork, liwork, ierr)
    if (p > 0) then
        call dlacpy('a', m, r+p, UU, m, UUU, m)
        call dlacpy('a', r+p, n, VVT, r+p, VVVT, r+p)
        call dlrort('l', m, n, r+p, UUU, m, VVVT, r+p, tau, wrk, lwrk, ierr)
        call dlrsvd_ort('s', 'l', m, n, r+p, UUU, m, VVVT, r+p, tau, SS, UU, m, VVT, r+p, wrk, lwrk, iwork, liwork, ierr)
        call ddgmm('r', m, r, UU, m, SS, 1, ierr)
    end if
    call dlacpy('a', m, r, UU, m, U, ldu)
    call dlacpy('a', r, n, VVT, r+p, VT, ldvt)
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
    work, lwork, iwork, liwork, ierr, shift_ratio, i1, j1)
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
    integer,    intent(in   ), optional     :: i1
    integer,    intent(in   ), optional     :: j1

!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRMA_MAXVOL_PROJ'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(21)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_UU,    &
                    i_UU, mem_VVT, i_VVT, ifoo(1), max_pos(2)
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
        if (present(i1)) then
            max_pos(1) = i1
        else
            max_pos(1) = 1
        end if
        if (present(j1)) then
            max_pos(2) = j1
        else
            max_pos(2) = 1
        end if
        call dgemaxvol(m, n, fun, 1, max_pos(1:), max_pos(2:), niter, ONE, wrk, lwrk, iwork, liwork, ierr)
        lowest = dlrval(m, n, max_pos(1), max_pos(2), r, U, ldu, VT, ldvt, ierr)
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
end submodule nonnegative_matrix_sub
