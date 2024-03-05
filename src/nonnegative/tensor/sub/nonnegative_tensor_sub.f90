submodule (nonnegative_tensor_mod) nonnegative_tensor_sub
implicit none (type, external)

contains
!----------------------------------------------------------------------------------------------------------------------------

module subroutine distance_to_nn &
(d, n, A, distc, distf)
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_comparison_mod,   only:   &
    safe_less
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: d
    integer,    intent(in   ),  contiguous  :: n(:)
    real(WP),   intent(in   ),  contiguous  :: A(:)
    real(WP),   intent(  out)               :: distc
    real(WP),   intent(  out)               :: distf

!-- Variables ------------------------------------------------------------------
    integer     :: i, numel

!-- Executable section ---------------------------------------------------------
    distc = ZERO
    distf = ZERO

    numel = product( n(1:d) )

    do i = 1, numel
    associate( val => A(i) )
        if ( safe_less(val, ZERO) ) then
            distf = distf + val**2
            distc = max(distc, -val)
        end if
    end associate
    end do

    distf = sqrt(distf)
end subroutine distance_to_nn

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrta_ttsvd &
(method, d, n, r, cores, work, lwork, iwork, liwork, ierr, shift_ratio)
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
use maria_tt_utils_mod,         only:   &
    dtt2full
use maria_tt_tsvd_mod,          only:   &
    dttsvd
use maria_reports_mod,      only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: d
    integer,    intent(in   ), contiguous   :: n(:)
    integer,    intent(in   ), contiguous   :: r(:)
    type(AR),   intent(inout)               :: cores(:)
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRTA_TTSVD'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(14)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_A,     &
                    i_A, mem_rr, i_rr, i_iwrk, liwrk, ifoo(1), i, numel
    real(WP)    ::  foo(1), lowest, ratio
    type(AR)    ::  arfoo(1)

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1) = (method < 1) .or. (method > 4)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    numel = product( n(1:d) )

    ! Storage [real]: A
    mem_A = numel
    ! Storage [int]: rr
    mem_rr = d-1
    ! Query procedures
    call dttsvd('l', d, n, foo, ifoo, arfoo, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dtt2full(d, n, r, arfoo, foo, foo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    ! Total
    req_lw = mem_A  + lw_proc
    req_liw = mem_rr + liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(7) = (lwork < req_lw)
    bad_arg(9) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..A..|..wrk..|
i_A = 1
i_wrk = i_A + mem_A
lwrk = lwork - i_wrk + 1
! Slice IWORK: |..rr..|..iwrk..|
i_rr = 1
i_iwrk = i_rr + mem_rr
liwrk = liwork - i_iwrk + 1
associate (A => work(i_A:), wrk => work(i_wrk:), rr => iwork(i_rr:), iwrk => iwork(i_iwrk:))
    call dtt2full(d, n, r, cores, A, wrk, lwrk, ierr)

    if (method == 1) then
        lowest = ZERO
    else
        ! Choose shift as the absolute value of the largest negative entry
        call distance_to_nn(d, n, A, lowest, foo(1))
        lowest = lowest * ratio
    end if

    ! Nonnegative step
    do i = 1, numel
        if (method == 1 .or. method == 3) then
            A(i) = max(lowest, A(i))
        elseif (method == 2) then
            A(i) = abs(A(i))
        else
            A(i) = A(i) + lowest
        end if
    end do

    ! Low-rank step
    rr(1:d-1) = r(1:d-1)
    call dttsvd('l', d, n, A, rr, cores, wrk, lwrk, iwork, liwork, ierr, r)
end associate
end subroutine nn_lrta_ttsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrta_ttcross &
(method, d, n, r, cores, row_inds, col_inds, niter, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, mi1)
use maria_kinds_mod,    only:   &
    WP  =>  DP
use maria_arr_mod,      only:   & 
    AR  =>  darr,               &
    iarr
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
use maria_comparison_mod,   only:   &
    safe_less
use maria_la_core_mod,      only:   &
    droundup_lwork
use maria_tt_utils_mod,         only:   &
    dtt2full
use maria_tt_cross_mod,          only:   &
    dttcross
use maria_access_tensor_mod,    only:   &
    dtenfib
use maria_reports_mod,      only:   &
    report_bad_arg
!-- Arguments ------------------------------------------------------------------
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: d
    integer,    intent(in   ), contiguous   :: n(:)
    integer,    intent(in   ), contiguous   :: r(:)
    type(AR),   intent(inout)               :: cores(:)
    type(iarr), intent(inout)               :: row_inds(:)
    type(iarr), intent(inout)               :: col_inds(:)
    integer,    intent(in   )               :: niter
    real(WP),   intent(in   )               :: thresh
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
    integer,    intent(in   ), optional     :: mi1(:)
   
!-- Parameters -----------------------------------------------------------------
    character(*), parameter :: SRNAME = 'NN_LRTA_TTCROSS'

!-- Variables ------------------------------------------------------------------
    logical     ::  lw_query, bad_arg(14)
    integer     ::  req_lw, req_liw, lw_proc, liw_proc, i_wrk, lwrk, mem_A,     &
                    i_A, mem_iwrk1, i_iwrk1, i_iwrk, liwrk, ifoo(1), numel
    real(WP)    ::  foo(1), lowest, ratio
    type(AR)    ::  arfoo(1)
    type(iarr)  :: iarfoo(1)
    procedure(dtenfib), pointer :: fun

    if (present(shift_ratio)) then
        ratio = shift_ratio
    else
        ratio = ONE
    end if

    if (present(mi1)) continue

!-- Sanity check ---------------------------------------------------------------
    ierr = 0
    bad_arg = .false.

    bad_arg(1) = (method < 1) .or. (method > 4)

!-- Quick return if possible ---------------------------------------------------
quick: block
    if (any(bad_arg)) exit quick
end block quick

!-- Estimate workspace ---------------------------------------------------------
workspace: block
    if (any(bad_arg)) exit workspace

    numel = product( n(1:d) )

    ! Storage [real]: A
    mem_A = numel
    ! Storage [int]: iwrk1
    mem_iwrk1 = d
    ! Query procedures
    call dttcross(fun, d, n, r, iarfoo, iarfoo, arfoo, niter, thresh, foo, -1, ifoo, -1, ierr)
    lw_proc = int(foo(1))
    liw_proc = ifoo(1)
    call dtt2full(d, n, r, arfoo, foo, foo, -1, ierr)
    lw_proc = max(lw_proc, int(foo(1)))
    ! Total
    req_lw = mem_A  + lw_proc
    req_liw = mem_iwrk1 + liw_proc

    lw_query = (lwork == -1 .or. liwork == -1)
    if (lw_query) then
        work(1) = droundup_lwork(req_lw)
        iwork(1) = req_liw
        return
    end if

    bad_arg(11) = (lwork < req_lw)
    bad_arg(13) = (liwork < req_liw) 
end block workspace

!-- Report incorrect arguments -------------------------------------------------
    if (any(bad_arg)) then
        ierr = -findloc(bad_arg, .true., dim=1)
        call report_bad_arg(SRNAME, -ierr)
        return
    end if

!-- Executable section ---------------------------------------------------------

! Slice work:   |..A..|..wrk..|
i_A = 1
i_wrk = i_A + mem_A
lwrk = lwork - i_wrk + 1
! Slice IWORK: |..iwrk1..|..iwrk..|
i_iwrk1 = 1
i_iwrk = i_iwrk1 + mem_iwrk1
liwrk = liwork - i_iwrk + 1
associate (A => work(i_A:), wrk => work(i_wrk:), iwrk1 => iwork(i_iwrk1:), iwrk => iwork(i_iwrk:))
    call dtt2full(d, n, r, cores, A, wrk, lwrk, ierr)

!    lowest = ZERO
!    if (method == 2) then
!        fun => dttfib_min
!                
!        lowest = dlrval(m, n, max_pos(1), max_pos(2), r, U, ldu, VT, ldvt, ierr)
!        lowest = min(lowest, ZERO)
!        lowest = -lowest * ratio
!    end if

    if (method == 1) then
        lowest = ZERO
    else
        ! Choose shift as the absolute value of the largest negative entry
        call distance_to_nn(d, n, A, lowest, foo(1))
        lowest = lowest * ratio
    end if

    ! Low-rank and nonnegative steps combined
    fun => dttfib_nn
    call dttcross(fun, d, n, r, row_inds, col_inds, cores, niter, thresh, wrk, lwrk, iwrk, liwrk, ierr)
end associate
contains
    function dttval_nn(d, n, mi, ierr)
    use maria_access_tensor_mod, only: &
        mi2i

        integer, intent(in) :: d
        integer, intent(in), contiguous :: n(:)
        integer, intent(in), contiguous :: mi(:)
        integer, intent(out) :: ierr
        real(WP)            :: dttval_nn

        integer :: i

        call mi2i(d, n, mi, i, ierr)
        dttval_nn = work(i_A + i - 1)
        if (method == 1 .or. method == 3) then
            dttval_nn = max(lowest, dttval_nn)
        elseif (method == 2) then
            dttval_nn = abs(dttval_nn)
        else
            dttval_nn = dttval_nn + lowest
        end if
    end function dttval_nn

    subroutine dttfib_nn(d, n, k, li, ri, x, incx, ierr)
    use maria_access_tensor_mod, only: &
        dtenval, dtenval2fib

        integer, intent(in) :: d
        integer, intent(in), contiguous :: n(:)
        integer, intent(in) :: k
        integer, intent(in), contiguous :: li(:)
        integer, intent(in), contiguous :: ri(:)
        real(WP), intent(out), contiguous :: x(:)
        integer, intent(in) :: incx
        integer, intent(out) :: ierr

        procedure(dtenval), pointer :: fun

        fun => dttval_nn
        call dtenval2fib(fun, d, n, k, li, ri, x, incx, iwork, mem_iwrk1, ierr)
    end subroutine dttfib_nn
end subroutine nn_lrta_ttcross

!----------------------------------------------------------------------------------------------------------------------------
end submodule nonnegative_tensor_sub
