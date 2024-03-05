program tensor_integer
use maria_kinds_mod,    only:   &
    WP  =>  DP
use maria_arr_mod,      only:   & 
    AR  =>  darr
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE, D_MACHTOL
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only:   &
    prng    => prng_mkl
#else
use maria_prng_builtin_mod, only:   &
    prng    => prng_builtin
#endif
use maria_comparison_mod,   only:   &
    safe_eq, safe_leq
use maria_utils_mod,        only:   &
    are_close
use maria_la_utils_mod,     only:   &
    drandort
use maria_la_core_mod,      only:   &
    dscal
use maria_tt_utils_mod,     only:   &
    ttrank, dtt2full
implicit none (type, external)

character(100)  :: seed_str
integer         :: d, maxn, maxr, ierr, seed, i, lwork, numel, rl, rr
type(prng)      :: rng
real(WP)        :: foo(1), maxerr
type(AR)        :: arfoo(1)

integer,    allocatable :: n(:), r(:)
real(WP),   allocatable :: A(:), work(:)
type(AR),   allocatable :: cores(:)

! Read command-line arguments
i = 1
call get_command_argument(i, seed_str)
read(seed_str, *) seed
i = i + 1
! Data description -------------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) d
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) maxn
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) maxr
i = i + 1

numel = maxn**d
allocate(n(d), r(d-1), A(numel), cores(d))
n = maxn
r = maxr
call rng%init(seed, ierr)

! Generate data ----------------------------------------------------------
do i = 1, d
    call ttrank(d, r, i, rl, rr, ierr)
    allocate( cores(i)%arr(rl * n(i) * rr) )
    call rng%dnormal(rl * n(i) * rr, cores(i)%arr, ZERO, ONE, ierr)
end do

call dtt2full(d, n, r, arfoo, foo, foo, -1, ierr)
lwork = int(foo(1))
allocate(work(lwork), stat=ierr)
call dtt2full(d, n, r, cores, A, work, lwork, ierr)
maxerr = ZERO
do i = 1, numel
    maxerr = max( maxerr, abs(A(i) - nint(A(i))) )
    A(i) = ONE * nint(A(i))
end do
deallocate(work)
print '( (1pe12.5) )', 0.5_WP - maxerr

! Compute its MAX_LRTA
call run_ap_ttsvd(d, n, A, r)

contains

    subroutine run_ap_ttsvd &
    (d, n, A, r)
    use maria_la_core_mod,      only:   &
        dgenrmc,                        &
        dcopy,                         &
        daxpy
    use maria_tt_utils_mod,     only:   &
        dtt2full
    use maria_tt_tsvd_mod,      only:   &
        dttsvd
    use maximum_tensor_mod,        only:   &
        max_lrta_ttsvd

        integer,    intent(in   )               :: d
        integer,    intent(in   ), contiguous   :: n(:)
        real(WP),   intent(in   ), contiguous   :: A(:)
        integer,    intent(in   ), contiguous   :: r(:)

        integer     :: ifoo(1), lwork, liwork, ierr
        real(WP)    :: foo(1), eps_min, eps_max, eps_cur, err_pre, err_cur
        type(AR)    :: arfoo(1)

        integer,    allocatable :: iwork(:), tmpr(:)
        real(WP),   allocatable :: B(:), work(:)
        type(AR),   allocatable :: cores(:)

        ! Estimate workspace
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

        ! Initialization
        call dcopy(numel, A, 1, B, 1)
        call dttsvd('l', d, n, B, tmpr, cores, work, lwork, iwork, liwork, ierr, r)

        ! Initial performance
        eps_min = ZERO
        call dtt2full(d, n, r, cores, B, work, lwork, ierr)
        call daxpy(numel, -ONE, A, 1, B, 1)
        err_pre = dgenrmc(1, numel, B, 1, ierr)
        eps_max = 0.5_WP
        print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"

        binsearch: do
            if ( are_close(eps_min, eps_max, ierr, rtol=1e-8_WP) ) then
                print '( (1pe12.5) )', max(err_pre, eps_max)
                return
            end if

            eps_cur = (eps_min + eps_max) / 2
            print '( (a), (1pe12.5) )', "Trying: ", eps_cur
            ap: do
                call max_lrta_ttsvd(d, n, A, r, cores, eps_cur, work, lwork, iwork, liwork, ierr)
                call dtt2full(d, n, r, cores, B, work, lwork, ierr)
                call daxpy(numel, -ONE, A, 1, B, 1)
                err_cur = dgenrmc(1, numel, B, 1, ierr)
                eps_max = min(err_cur, eps_max)
                print '( (1pe12.5) )', 1 - err_cur / err_pre
                if ( safe_leq(err_cur, eps_cur) ) then
                    ! Reached desired max-error
                    exit ap
                elseif ( err_cur > err_pre * (ONE-5e-6_WP) ) then 
                    ! Alternating projections converged but max-error not reached
                    eps_min = eps_min + (eps_max - eps_min) / 4
                    print '( (a), (1pe12.5) )', "Reached: ", err_cur
                    exit ap
                end if
                err_pre = err_cur
            end do ap

            print '( (a), (1pe12.5), (a), (1pe12.5), (a) )', "Interval: (", eps_min, ", ", eps_max, "]"
        end do binsearch
    end subroutine run_ap_ttsvd

end program tensor_integer
