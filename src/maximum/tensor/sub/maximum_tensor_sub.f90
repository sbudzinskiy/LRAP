submodule (maximum_tensor_mod) maximum_tensor_sub
implicit none (type, external)

contains
!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrta_ttsvd &
(d, n, A, r, cores, tol, work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP  =>  DP
use maria_arr_mod,      only:   & 
    AR  =>  darr
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_comparison_mod,   only:   &
    safe_less
use maria_la_core_mod,      only:   &
    droundup_lwork
use maria_tt_utils_mod,     only:   &
    dtt2full
use maria_tt_tsvd_mod,      only:   &
    dttsvd
use maria_reports_mod,      only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: d
    integer,    intent(in   ), contiguous   :: n(:)
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   ), contiguous   :: r(:)
    type(AR),   intent(inout)               :: cores(:)
    real(WP),   intent(in   )               :: tol
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'MAX_LRTA_TTSVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(14)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_B,     &
                    i_B, mem_rr, i_rr, i_iwrk, liwrk, i, numel, ifoo(1)
    real(WP)    ::  foo(1)
    type(AR)    ::  arfoo(1)

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    numel = product( n(1:d) )

    ! Storage [real]: B
    mem_B = numel
    ! Storage [int]: rr
    mem_rr = d-1
    ! Query procedures
    call dttsvd('l', d, n, foo, ifoo, arfoo, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dtt2full(d, n, r, arfoo, foo, foo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    ! Total
    req_lw = mem_B  + lw_proc
    req_liw = mem_rr + liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(8) = (lwork < req_lw)
    bad_arg(10) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..B..|..wrk..|
i_B = 1
i_wrk = i_B + mem_B
lwrk = lwork - i_wrk + 1
! Slice IWORK: |..rr..|..iwrk..|
i_rr = 1
i_iwrk = i_rr + mem_rr
liwrk = liwork - i_iwrk + 1
associate (B => work(i_B:), wrk => work(i_wrk:), rr => iwork(i_rr:), iwrk => iwork(i_iwrk:))
    call dtt2full(d, n, r, cores, B, wrk, lwrk, ierr)
    
    ! Max-ball step
    do i = 1, numel
    associate( val0 => A(i),  val1 => B(i) )
        val1 = val1 - val0
        val1 = min(val1, tol)
        val1 = max(val1, -tol)
        val1 = val0 + val1
    end associate
    end do

    ! Low-rank step
    rr(1:d-1) = r(1:d-1)
    call dttsvd('l', d, n, B, rr, cores, wrk, lwrk, iwork, liwork, ierr, r)
end associate
end subroutine max_lrta_ttsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_ttsvd &
(d, n, A, r, rng, c1, c2, c3, c4)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_arr_mod,      only:   & 
    AR  =>  darr
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_prng_mod,     only:   &
    prng
use maria_la_core_mod,  only:   &
    dgenrmc,                    &
    dcopy,                      &
    daxpy
use maria_tt_utils_mod, only:   &
    dtt2full, ttrank
use maria_tt_tsvd_mod,  only:   &
    dttsvd
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in)              :: d
    integer,    intent(in), contiguous  :: n(:)
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in), contiguous  :: r(:)
    class(prng),intent(in)              :: rng
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4

!-- Variables ------------------------------------------------------------------
    integer     :: ifoo(1), lwork, liwork, ierr, rl, rr, i, numel
    real(WP)    :: foo(1), eps_min, eps_max, eps_cur, err_pre, err_cur, eps_len
    type(AR)    :: arfoo(1)

    integer,    allocatable :: iwork(:), tmpr(:)
    real(WP),   allocatable :: B(:), work(:)
    type(AR),   allocatable :: cores(:)

!-- Estimate workspace ---------------------------------------------------------
    call max_lrta_ttsvd(d, n, foo, r, arfoo, eps_cur, foo, -1, ifoo, -1, ierr)
    lwork = int(foo(1))
    liwork = ifoo(1)
    call dttsvd('l', d, n, foo, ifoo, arfoo, foo, -1, ifoo, -1, ierr, r)
    lwork = max(lwork, int(foo(1)))
    liwork = max(liwork, ifoo(1))
    call dtt2full(d, n, r, arfoo, foo, foo, -1, ierr)
    lwork = max(lwork, int(foo(1)))

    numel = product(n(1:d))
    allocate(B(numel), tmpr(d-1), cores(d))
    tmpr(1:d-1) = r(1:d-1)
    allocate(work(lwork), iwork(liwork))

    ! Random initial approximation
    call dcopy(numel, A, 1, B, 1)
    do i = 1, d
        call ttrank(d, r, i, rl, rr, ierr)
        allocate( cores(i)%arr(rl * n(i) * rr) )
        call rng%dnormal(rl * n(i) * rr, cores(i)%arr, ZERO, sqrt(ONE / rl / rr), ierr)
    end do

    ! Initial performance
    eps_min = ZERO
    call dtt2full(d, n, r, cores, B, work, lwork, ierr)
    call daxpy(numel, -ONE, A, 1, B, 1)
    err_pre = dgenrmc(1, numel, B, 1, ierr)
    eps_max = err_pre
!    print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"

    binsearch: do
        eps_len = eps_max - eps_min
        if ( eps_len < c1 * eps_max ) then
            print '( (1pe12.5) )', eps_max
            return
        end if

        eps_cur = (eps_min + eps_max) / 2
!            print '( (a), (1pe12.5) )', "Trying: ", eps_cur
        ap: do
            call max_lrta_ttsvd(d, n, A, r, cores, eps_cur, work, lwork, iwork, liwork, ierr)
            call dtt2full(d, n, r, cores, B, work, lwork, ierr)
            call daxpy(numel, -ONE, A, 1, B, 1)
            err_cur = dgenrmc(1, numel, B, 1, ierr)
            eps_max = min(err_cur, eps_max)
!                print '( (1pe12.5) )', 1 - err_cur / err_pre
            if ( err_cur > err_pre * (ONE - c2) ) then
                ! Convergence of AP slowed down
                if ( err_cur > eps_min + eps_len * c3 ) then
                    ! Too far, increase lower bound
                    eps_min = eps_min + eps_len * c4
                end if
!                print '( (a), (1pe12.5) )', "Reached: ", err_cur
                exit ap
            end if
            err_pre = err_cur
        end do ap

!            print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"
    end do binsearch
end subroutine run_ap_ttsvd

!----------------------------------------------------------------------------------------------------------------------------
end submodule maximum_tensor_sub
