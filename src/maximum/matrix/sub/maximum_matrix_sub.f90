submodule (maximum_matrix_mod) maximum_matrix_sub
implicit none (type, external)

contains
!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrma_svd &
(m, n, A, lda, r, U, ldu, VT, ldvt, tol, work, lwork, iwork, liwork, ierr)
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
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   )               :: lda
    integer,    intent(in   )               :: r                            !  4
    real(WP),   intent(inout), contiguous   :: U(:)                         !  5
    integer,    intent(in   )               :: ldu                          !  6
    real(WP),   intent(inout), contiguous   :: VT(:)                        !  7
    integer,    intent(in   )               :: ldvt
    real(WP),   intent(in   )               :: tol
    real(WP),   intent(  out), contiguous   :: work(:)                      !  9
    integer,    intent(in   )               :: lwork                        ! 10
    integer,    intent(  out), contiguous   :: iwork(:)                     ! 11
    integer,    intent(in   )               :: liwork                       ! 12
    integer,    intent(  out)               :: ierr                         ! 13
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'MAX_LRMA_SVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(14)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_B,     &
                    i_B, mem_bigU, i_bigU, mem_bigVT, i_bigVT, mem_bigS,        &
                    i_bigS, i, j, ifoo(1)
    real(WP)    ::  foo(1)

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: B, bigU, bigS, bigVT
    mem_B = m * n
    mem_bigU = m * min(m,n)
    mem_bigS = min(m,n)
    mem_bigVT = min(m,n) * n
    ! Query procedures
    call dgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    ! Total
    req_lw = mem_B + mem_bigU + mem_bigS + mem_bigVT + lw_proc
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

! Slice work:   |..B..|..bigS..|..bigU..|..bigVT..|..wrk..|
i_B = 1
i_bigS = i_B + mem_B
i_bigU = i_bigS + mem_bigS
i_bigVT = i_bigU + mem_bigU
i_wrk = i_bigVT + mem_bigVT
lwrk = lwork - i_wrk + 1
associate (B => work(i_B:), bigS => work(i_bigS:), bigU => work(i_bigU:),       &
bigVT => work(i_bigVT:), wrk => work(i_wrk:))
    call dgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, B, m)
    
    ! Max-ball step
    do j = 1, n
        do i = 1, m
        associate( val0 => A(i+(j-1)*lda),  val1 => B(i + (j-1)*m) )
            val1 = val1 - val0
            val1 = min(val1, tol)
            val1 = max(val1, -tol)
            val1 = val0 + val1
        end associate
        end do
    end do

    ! Low-rank step
    call dgesdd_q('s', m, n, B, m, bigS, bigU, m, bigVT, min(m,n), wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, bigU, m, U, ldu)
    call ddgmm('r', m, r, U, ldu, bigS, 1, ierr)
    call dlacpy('a', r, n, bigVT, min(m,n), VT, ldvt)
end associate
end subroutine max_lrma_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_svd &
(m, n, A, r, rng, init_randort, c1, c2, c3, c4)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_prng_mod,         only:   &
    prng
use maria_la_core_mod,      only:   &
    dgenrmc,                        &
    daxpy,                          &
    dscal
use maria_la_utils_mod,     only:   &
    drandort
use maria_lr_la_mod,        only:   &
    dlr2full
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in)              :: m
    integer,    intent(in)              :: n
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in)              :: r
    class(prng),intent(in)              :: rng
    logical,    intent(in)              :: init_randort
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4

!-- Variables ------------------------------------------------------------------
    integer     :: ifoo(1), lwork, liwork, ierr
    real(WP)    :: foo(1), eps_min, eps_max, eps_cur, err_pre, err_cur, nrmA, eps_len

    integer,    allocatable :: iwork(:)
    real(WP),   allocatable :: AA(:), U(:), VT(:), work(:)

!-- Estimate workspace ---------------------------------------------------------
    call max_lrma_svd(m, n, foo, m, r, foo, m, foo, r, ONE, foo, -1, ifoo, -1, ierr)
    lwork = int(foo(1))
    liwork = ifoo(1)
    if (init_randort) then
        call drandort(rng, m, r, foo, m, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call drandort(rng, r, n, foo, r, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
    end if

    allocate(AA(m*n), U(m*r), VT(r*n), work(lwork), iwork(liwork))

! Random initial approximation
    if (init_randort) then
        call drandort(rng, m, r, U, m, work, lwork, ierr)
        call drandort(rng, r, n, VT, r, work, lwork, ierr)
    else
        call rng%dnormal(m*r, U, ZERO, ONE, ierr)
        call rng%dnormal(r*n, VT, ZERO, ONE, ierr)
        call dscal(m*r, ONE / r, U, 1)
    end if

! Norm of A
    nrmA = dgenrmc(m, n, A, m, ierr)
    print '( (1pe12.5) )', nrmA

! Initial performance
    eps_min = ZERO
    call dlr2full(m, n, r, U, m, VT, r, AA, m, ierr)
    call daxpy(m*n, -ONE, A, 1, AA, 1)
    err_pre = dgenrmc(m, n, AA, m, ierr)
    eps_max = err_pre
    print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"

! Iterations
    binsearch: do
        eps_len = eps_max - eps_min
        if ( eps_len < c1 * eps_max ) then
            print '( (1pe12.5) )', eps_max
            return
        end if

        eps_cur = (eps_min + eps_max) / 2
        print '( (a), (1pe12.5) )', "Trying: ", eps_cur
        ap: do
            call max_lrma_svd(m, n, A, m, r, U, m, VT, r, eps_cur, work, lwork, iwork, liwork, ierr)
            call dlr2full(m, n, r, U, m, VT, r, AA, m, ierr)
            call daxpy(m*n, -ONE, A, 1, AA, 1)
            err_cur = dgenrmc(m, n, AA, m, ierr)
            eps_max = min(err_cur, eps_max)
            print '( (1pe12.5) )', 1 - err_cur / err_pre
            if ( err_cur > err_pre * (ONE - c2) ) then
                ! Convergence of AP slowed down
                if ( err_cur > eps_min + eps_len * c3 ) then
                    ! Too far, increase lower bound
                    eps_min = eps_min + eps_len * c4
                end if
                print '( (a), (1pe12.5) )', "Reached: ", err_cur
                exit ap
            end if
            err_pre = err_cur
        end do ap

        print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"
    end do binsearch
end subroutine run_ap_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrma_rsvd &
(m, n, A, lda, r, U, ldu, VT, ldvt, tol, rng, p, q, work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
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
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   )               :: lda
    integer,    intent(in   )               :: r                            !  4
    real(WP),   intent(inout), contiguous   :: U(:)                         !  5
    integer,    intent(in   )               :: ldu                          !  6
    real(WP),   intent(inout), contiguous   :: VT(:)                        !  7
    integer,    intent(in   )               :: ldvt
    real(WP),   intent(in   )               :: tol
    class(prng),intent(in   )               :: rng
    integer,    intent(in   )               :: p
    integer,    intent(in   )               :: q
    real(WP),   intent(  out), contiguous   :: work(:)                      !  9
    integer,    intent(in   )               :: lwork                        ! 10
    integer,    intent(  out), contiguous   :: iwork(:)                     ! 11
    integer,    intent(in   )               :: liwork                       ! 12
    integer,    intent(  out)               :: ierr                         ! 13
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'MAX_LRMA_RSVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(18)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_B,     &
                    i_B, mem_gauss, i_gauss, mem_csk, i_csk, mem_bigU, i_bigU,  &
                    mem_bigVT, i_bigVT, mem_bigS, i_bigS, i, j, ifoo(1), ldb
    real(WP)    ::  foo(1)
    procedure(MM), pointer :: C2BC, C2CB

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick

    if (m == 0 .or. n == 0 .or. r == 0) return
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    ! Storage [real]: B, gauss, csk, bigU, bigS, bigVT
    ldb = m
    mem_B = ldb * n
    mem_gauss = n * (r+p)
    mem_csk = m * (r+p)
    mem_bigU = m * (r+p)
    mem_bigS = r+p
    mem_bigVT = (r+p) * n
    ! Query procedures
    call dgersvd1(m, n, C2BC, C2CB, r+p, foo, m, q, foo, foo, m, foo, r+p, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    ! Total
    req_lw = mem_B + mem_gauss + mem_csk + mem_bigU + mem_bigS + mem_bigVT + lw_proc
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

C2BC => mm_B_left
C2CB => mm_B_right

! Slice work:   |..B..|..gauss..|..csk..|..bigS..|..bigU..|..bigVT..|..wrk..|
i_B = 1
i_gauss = i_B + mem_B
i_csk = i_gauss + mem_gauss
i_bigS = i_csk + mem_csk
i_bigU = i_bigS + mem_bigS
i_bigVT = i_bigU + mem_bigU
i_wrk = i_bigVT + mem_bigVT
lwrk = lwork - i_wrk + 1
associate (B => work(i_B:), gauss => work(i_gauss:), csk => work(i_csk:), &
        bigS => work(i_bigS:), bigU => work(i_bigU:), bigVT => work(i_bigVT:), wrk => work(i_wrk:))
    call dgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, B, m)
    
    ! Max-ball step
    do j = 1, n
        do i = 1, m
        associate( val0 => A(i+(j-1)*lda),  val1 => B(i + (j-1)*m) )
            val1 = val1 - val0
            val1 = min(val1, tol)
            val1 = max(val1, -tol)
            val1 = val0 + val1
        end associate
        end do
    end do

    ! Low-rank step
    call rng%dnormal(n*(r+p), gauss, ZERO, ONE, ierr)
    call dgemm('n', 'n', m, r+p, n, ONE, B, ldb, gauss, n, ZERO, csk, m)
    call dgersvd1(m, n, C2BC, C2CB, r+p, csk, m, q, bigS, bigU, m, bigVT, r+p, wrk, lwrk, iwork, liwork, ierr)
    call dlacpy('a', m, r, bigU, m, U, ldu)
    call ddgmm('r', m, r, U, ldu, bigS, 1, ierr)
    call dlacpy('a', r, n, bigVT, r+p, VT, ldvt)
end associate
contains
    subroutine mm_B_left &
    (transX, m, n, k, alpha, X, ldx, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: ierr

    associate( B => work(i_B:) )
        call dgemm('n', transX, m, n, k, alpha, B, ldb, X, ldx, beta, C, ldc)
        ierr = 0
    end associate
    end subroutine mm_B_left

    subroutine mm_B_right &
    (transX, m, n, k, alpha, X, ldx, beta, C, ldc, ierr)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_la_core_mod,  only:   &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: ierr

    associate( B => work(i_B:) )
        call dgemm(transX, 'n', m, n, k, alpha, X, ldx, B, ldb, beta, C, ldc)
        ierr = 0
    end associate
    end subroutine mm_B_right
end subroutine max_lrma_rsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_rsvd &
(m, n, A, r, p, q, rng, init_randort, c1, c2, c3, c4)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_prng_mod,     only:   &
    prng
use maria_la_core_mod,      only:   &
    dgenrmc,                        &
    daxpy,                          &
    dscal
use maria_la_utils_mod,     only:   &
    drandort
use maria_lr_la_mod,        only:   &
    dlr2full
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in)              :: m
    integer,    intent(in)              :: n
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in)              :: r
    integer,    intent(in)              :: p
    integer,    intent(in)              :: q
    class(prng),intent(in)              :: rng
    logical,    intent(in)              :: init_randort
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4

!-- Variables ------------------------------------------------------------------
    integer     :: ifoo(1), lwork, liwork, ierr
    real(WP)    :: foo(1), eps_min, eps_max, eps_cur, err_pre, err_cur, nrmA, eps_len

    integer,    allocatable :: iwork(:)
    real(WP),   allocatable :: AA(:), U(:), VT(:), work(:)

!-- Estimate workspace ---------------------------------------------------------
    call max_lrma_rsvd(m, n, foo, m, r, foo, m, foo, r, ONE, rng, p, q, foo, -1, ifoo, -1, ierr)
    lwork = int(foo(1))
    liwork = ifoo(1)
    if (init_randort) then
        call drandort(rng, m, r, foo, m, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call drandort(rng, r, n, foo, r, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
    end if

    allocate(AA(m*n), U(m*r), VT(r*n), work(lwork), iwork(liwork))

! Random initial approximation
    if (init_randort) then
        call drandort(rng, m, r, U, m, work, lwork, ierr)
        call drandort(rng, r, n, VT, r, work, lwork, ierr)
    else
        call rng%dnormal(m*r, U, ZERO, ONE, ierr)
        call rng%dnormal(r*n, VT, ZERO, ONE, ierr)
        call dscal(m*r, ONE / r, U, 1)
    end if

! Norm of A
    nrmA = dgenrmc(m, n, A, m, ierr)
    print '( (1pe12.5) )', nrmA

! Initial performance
    eps_min = ZERO
    call dlr2full(m, n, r, U, m, VT, r, AA, m, ierr)
    call daxpy(m*n, -ONE, A, 1, AA, 1)
    err_pre = dgenrmc(m, n, AA, m, ierr)
    eps_max = err_pre
    print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"

! Iterations
    binsearch: do
        eps_len = eps_max - eps_min
        if ( eps_len < c1 * eps_max ) then
            print '( (1pe12.5) )', eps_max
            return
        end if

        eps_cur = (eps_min + eps_max) / 2
        print '( (a), (1pe12.5) )', "Trying: ", eps_cur
        ap: do
            call max_lrma_rsvd(m, n, A, m, r, U, m, VT, r, eps_cur, rng, p, q, work, lwork, iwork, liwork, ierr)
            call dlr2full(m, n, r, U, m, VT, r, AA, m, ierr)
            call daxpy(m*n, -ONE, A, 1, AA, 1)
            err_cur = dgenrmc(m, n, AA, m, ierr)
            eps_max = min(err_cur, eps_max)
            print '( (1pe12.5) )', 1 - err_cur / err_pre
            if ( err_cur > err_pre * (ONE - c2) ) then
                ! Convergence of AP slowed down
                if ( err_cur > eps_min + eps_len * c3 ) then
                    ! Too far, increase lower bound
                    eps_min = eps_min + eps_len * c4
                end if
                print '( (a), (1pe12.5) )', "Reached: ", err_cur
                exit ap
            end if
            err_pre = err_cur
        end do ap

        print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"
    end do binsearch
end subroutine run_ap_rsvd

!----------------------------------------------------------------------------------------------------------------------------
end submodule maximum_matrix_sub
