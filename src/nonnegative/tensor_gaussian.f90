program tensor_gaussian
use maria_kinds_mod,        only:   &
    WP      => DP
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
implicit none (type, external)

character(100)  :: seed_str
integer         :: d, method_ap, method_lrproj, ierr, seed, niter_ap, niter_maxvol, i, numel
type(prng)      :: rng
real(WP)        :: rtolf, thresh, shift_ratio

integer,    allocatable :: n(:)
real(WP),   allocatable :: A(:)

! Read command-line arguments
i = 1
call get_command_argument(i, seed_str)
read(seed_str, *) seed
i = i + 1
! Data description -------------------------------------------------------------
d = 10
allocate(n(d))
n = 2
call get_command_argument(i, seed_str)
read(seed_str, *) rtolf
i = i + 1
! Alternating projections ------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) method_ap
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) niter_ap
i = i + 1
if (method_ap >= 3) then
    call get_command_argument(i, seed_str)
    read(seed_str, *) shift_ratio
    i = i + 1
end if
! Low-rank projection ----------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) method_lrproj
i = i + 1
if (method_lrproj == 2) then
    call get_command_argument(i, seed_str)
    read(seed_str, *) niter_maxvol
    i = i + 1
    call get_command_argument(i, seed_str)
    read(seed_str, *) thresh
    i = i + 1
end if

numel = product( n(1:d) )

allocate(A(numel))
call rng%init(seed, ierr)

! Read data from file ----------------------------------------------------------
open(unit=42, file="mixture_gaussian_2**10.data")
read(42, *) A
close(42)

! Compute its NN_LRTA
if (method_lrproj == 1) then
    call run_ap_ttsvd(d, n, A, rtolf, method_ap, niter_ap, shift_ratio)
end if

contains

    subroutine run_ap_ttsvd &
    (d, n, A, rtolf, method_ap, niter_ap, shift_ratio)
    use maria_arr_mod,      only:   & 
        AR  =>  darr
    use maria_la_core_mod,      only:   &
        dgenrmf,                        &
        dgenrmc,                        &
        dcopy
    use maria_tt_utils_mod,        only:   &
        dtt2full
    use maria_tt_tsvd_mod,          only:   &
        dttsvd
    use nonnegative_tensor_mod,        only:   &
        nn_lrta_ttsvd

        integer,    intent(in)              :: d
        integer,    intent(in), contiguous  :: n(:)
        real(WP),   intent(in), contiguous  :: A(:)
        real(WP),   intent(in)              :: rtolf
        integer,    intent(in)              :: method_ap
        integer,    intent(in)              :: niter_ap
        real(WP),   intent(in)              :: shift_ratio

        integer     :: ifoo(1), lwork, liwork, ierr, i
        real(WP)    :: foo(1), nrmf, nrmc

        type(AR),   allocatable :: cores(:)
        integer,    allocatable :: iwork(:), r(:)
        real(WP),   allocatable :: AA(:), work(:)

        allocate(r(d-1), cores(d))
        r = 50
        ! Estimate workspace
        call dttsvd('l', d, n, foo, ifoo, cores, foo, -1, ifoo, -1, ierr, rtolf=rtolf)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call nn_lrta_ttsvd(method_ap, d, n, r, cores, foo, -1, ifoo, -1, ierr, shift_ratio)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        numel = product( n(1:d) )
        allocate(AA(numel), work(lwork), iwork(liwork))

        ! Norm of A
        nrmf = dgenrmf(1, numel, A, 1, ierr)
        nrmc = dgenrmc(1, numel, A, 1, ierr)

        ! Best rank-r approximation
        call dcopy(numel, A, 1, AA, 1)
        call dttsvd('l', d, n, AA, r, cores, work, lwork, iwork, liwork, ierr, rtolf=rtolf)
        ! Initial performance
        call dtt2full(d, n, r, cores, AA, work, lwork, ierr)
        open(unit=42, file="initial.data", status="replace")
        write(42, '( (1pe15.8) )') AA
        close(42)
        call log_performance(0, d, n, A, AA, nrmf, nrmc)

        do i = 1, niter_ap
            call nn_lrta_ttsvd(method_ap, d, n, r, cores, work, lwork, iwork, liwork, ierr, shift_ratio)
            call dtt2full(d, n, r, cores, AA, work, lwork, ierr)
            call log_performance(i, d, n, A, AA, nrmf, nrmc)
        end do
        call dtt2full(d, n, r, cores, AA, work, lwork, ierr)
        open(unit=42, file="final.data", status="replace")
        write(42, '( (1pe15.8) )') AA
        close(42)
    end subroutine run_ap_ttsvd

    subroutine log_performance &
    (iter, d, n, A, AA, nrmf, nrmc)
    use maria_la_core_mod,  only:   &
        dgenrmf, dgenrmc, daxpy
    use nonnegative_tensor_mod,    only:   &
        distance_to_nn

        integer,    intent(in   )               :: iter
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        real(WP),   intent(in   ),  contiguous  :: A(:)
        real(WP),   intent(inout),  contiguous  :: AA(:)
        real(WP),   intent(in   )               :: nrmf
        real(WP),   intent(in   )               :: nrmc

        integer     :: numel
        real(WP)    :: rerrf, rerrc, distc, distf

        numel = product( n(1:d) )
        call distance_to_nn(d, n, AA, distc, distf)
        print '( (I0), (1pe12.5) )', iter, distf / nrmf
        print '( (I0), (1pe12.5) )', iter, distc / nrmc

        call daxpy(numel, -ONE, A, 1, AA, 1)
        rerrf = dgenrmf(1, numel, AA, 1, ierr) / nrmf
        rerrc = dgenrmc(1, numel, AA, 1, ierr) / nrmc
        print '( (I0), (1pe12.5) )', iter, rerrf
        print '( (I0), (1pe12.5) )', iter, rerrc
        flush(6)
    end subroutine log_performance
end program tensor_gaussian
