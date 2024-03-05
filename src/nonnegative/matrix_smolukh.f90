program matrix_smolukh
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
integer         :: m, n, r, k, p, q, method_ap, method_lrproj, ierr, seed, niter_ap, niter_maxvol, i
type(prng)      :: rng
real(WP)        :: thresh, shift_ratio

real(WP),   allocatable :: A(:)

! Read command-line arguments
i = 1
call get_command_argument(i, seed_str)
read(seed_str, *) seed
i = i + 1
! Data description -------------------------------------------------------------
m = 1024
n = 1024
call get_command_argument(i, seed_str)
read(seed_str, *) r
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
    read(seed_str, *) p
    i = i + 1
    call get_command_argument(i, seed_str)
    read(seed_str, *) q
    i = i + 1
elseif (method_lrproj >= 3) then
    call get_command_argument(i, seed_str)
    read(seed_str, *) niter_maxvol
    i = i + 1
    call get_command_argument(i, seed_str)
    read(seed_str, *) thresh
    i = i + 1
end if
if (method_lrproj == 4) then
    call get_command_argument(i, seed_str)
    read(seed_str, *) k
    i = i + 1
end if

allocate(A(m*n))
call rng%init(seed, ierr)

! Read data from file ----------------------------------------------------------
open(unit=42, file="smolukh_1024.data")
read(42, *) A
close(42)

! Compute its NN_LRMA
if (method_lrproj == 1) then
    call run_ap_svd(m, n, A, r, method_ap, niter_ap, shift_ratio)
elseif (method_lrproj == 2) then
    call run_ap_rsvd(m, n, A, r, method_ap, niter_ap, shift_ratio, rng, p, q)
elseif (method_lrproj == 3) then
    call run_ap_maxvol(rng, m, n, A, r, method_ap, niter_ap, niter_maxvol, thresh, shift_ratio, 0)
else
    call run_ap_maxvol_proj(rng, m, n, A, r, method_ap, niter_ap, k, k, niter_maxvol, thresh, shift_ratio)
end if

contains

    subroutine run_ap_svd &
    (m, n, A, r, method, niter_ap, shift_ratio)
    use maria_la_core_mod,      only:   &
        dgenrmf,                        &
        dgenrmc,                        &
        dlacpy,                         &
        ddgmm,                          &
        dgesdd_q
    use maria_lr_la_mod,        only:   &
        dlr2full
    use nonnegative_matrix_mod,        only:   &
        nn_lrma_svd

        integer,    intent(in)              :: m
        integer,    intent(in)              :: n
        real(WP),   intent(in), contiguous  :: A(:)
        integer,    intent(in)              :: r
        integer,    intent(in)              :: method
        integer,    intent(in)              :: niter_ap
        real(WP),   intent(in)              :: shift_ratio

        integer     :: ifoo(1), lwork, liwork, ierr, i
        real(WP)    :: foo(1), nrmf, nrmc

        integer,    allocatable :: iwork(:)
        real(WP),   allocatable :: AA(:), S(:), U(:), VT(:), work(:)

        ! Estimate workspace
        call dgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call nn_lrma_svd(method, m, n, r, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr, shift_ratio)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(AA(m*n), S(min(m,n)), U(m*min(m,n)), VT(min(m,n)*n), work(lwork), iwork(liwork))

        ! Norm of A
        nrmf = dgenrmf(m, n, A, m, ierr)
        nrmc = dgenrmc(m, n, A, m, ierr)

        ! Best rank-r approximation
        call dlacpy('a', m, n, A, m, AA, m)
        call dgesdd_q('s', m, n, AA, m, S, U, m, VT, min(m,n), work, lwork, iwork, liwork, ierr)
        call ddgmm('r', m, r, U, m, S, 1, ierr)

        ! Initial performance
        call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
        call log_performance(0, m, n, A, AA, nrmf, nrmc)

        do i = 1, niter_ap
            call nn_lrma_svd(method, m, n, r, U, m, VT, min(m,n), work, lwork, iwork, liwork, ierr, shift_ratio)
            call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
            call log_performance(i, m, n, A, AA, nrmf, nrmc)
        end do
    end subroutine run_ap_svd

    subroutine run_ap_rsvd &
    (m, n, A, r, method, niter_ap, shift_ratio, rng, p, q)
    use maria_la_core_mod,      only:   &
        dgenrmf,                        &
        dgenrmc,                        &
        dlacpy,                         &
        ddgmm,                          &
        dgesdd_q
    use maria_lr_la_mod,        only:   &
        dlr2full
    use nonnegative_matrix_mod,        only:   &
        nn_lrma_rsvd

        integer,    intent(in)              :: m
        integer,    intent(in)              :: n
        real(WP),   intent(in), contiguous  :: A(:)
        integer,    intent(in)              :: r
        integer,    intent(in)              :: method
        integer,    intent(in)              :: niter_ap
        real(WP),   intent(in)              :: shift_ratio
        type(prng), intent(in)              :: rng
        integer,    intent(in)              :: p
        integer,    intent(in)              :: q

        integer     :: ifoo(1), lwork, liwork, ierr, i
        real(WP)    :: foo(1), nrmf, nrmc

        integer,    allocatable :: iwork(:)
        real(WP),   allocatable :: AA(:), S(:), U(:), VT(:), work(:)

        ! Estimate workspace
        call dgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call nn_lrma_rsvd(method, m, n, r, foo, m, foo, min(m,n), rng, p, q, foo, -1, ifoo, -1, ierr, shift_ratio)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(AA(m*n), S(min(m,n)), U(m*min(m,n)), VT(min(m,n)*n), work(lwork), iwork(liwork))

        ! Norm of A
        nrmf = dgenrmf(m, n, A, m, ierr)
        nrmc = dgenrmc(m, n, A, m, ierr)

        ! Best rank-r approximation
        call dlacpy('a', m, n, A, m, AA, m)
        call dgesdd_q('s', m, n, AA, m, S, U, m, VT, min(m,n), work, lwork, iwork, liwork, ierr)
        call ddgmm('r', m, r, U, m, S, 1, ierr)

        ! Initial performance
        call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
        call log_performance(0, m, n, A, AA, nrmf, nrmc)

        do i = 1, niter_ap
            call nn_lrma_rsvd(method, m, n, r, U, m, VT, min(m,n), rng, p, q, work, lwork, iwork, liwork, ierr, shift_ratio)
            call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
            call log_performance(i, m, n, A, AA, nrmf, nrmc)
        end do
    end subroutine run_ap_rsvd

    subroutine run_ap_maxvol &
    (rng, m, n, A, r, method, niter_ap, niter_maxvol, thresh, shift_ratio, p)
    use maria_utils_mod,        only:   &
        arange,                         &
        permute
    use maria_la_core_mod,      only:   &
        dgenrmf,                        &
        dgenrmc,                        &
        dlacpy,                         &
        ddgmm,                          &
        dgesdd_q
    use maria_lr_la_mod,        only:   &
        dlr2full
    use nonnegative_matrix_mod,        only:   &
        nn_lrma_maxvol

        class(prng), intent(in)              :: rng
        integer,    intent(in)              :: m
        integer,    intent(in)              :: n
        real(WP),   intent(in), contiguous  :: A(:)
        integer,    intent(in)              :: r
        integer,    intent(in)              :: method
        integer,    intent(in)              :: niter_ap
        integer,    intent(in)              :: niter_maxvol
        real(WP),   intent(in)              :: thresh
        real(WP),   intent(in)              :: shift_ratio
        integer,    intent(in)              :: p

        integer     :: ifoo(2), lwork, liwork, ierr, i
        real(WP)    :: foo(1), nrmf, nrmc

        integer,    allocatable :: iwork(:), icol(:), irow(:)
        real(WP),   allocatable :: AA(:), S(:), U(:), VT(:), work(:)

        ! Estimate workspace
        call dgesdd_q('s', m, m, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call nn_lrma_maxvol(method, m, n, r, foo, m, foo, min(m,n), p, ifoo, ifoo, niter_maxvol,   &
            thresh, foo, -1, ifoo, -1, ierr, shift_ratio)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(AA(m*n), S(min(m,n)), U(m*min(m,n)), VT(min(m,n)*n), work(lwork), iwork(liwork), irow(m), icol(n))

        ! Norm of A
        nrmf = dgenrmf(m, n, A, m, ierr)
        nrmc = dgenrmc(m, n, A, m, ierr)

        ! Best rank-r approximation
        call dlacpy('a', m, n, A, m, AA, m)
        call dgesdd_q('s', m, n, AA, m, S, U, m, VT, min(m,n), work, lwork, iwork, liwork, ierr)
        call ddgmm('r', m, r, U, m, S, 1, ierr)

        ! Initial performance
        call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
        call log_performance(0, m, n, A, AA, nrmf, nrmc)

            call arange(m, irow, 1, 1, ierr)
            call arange(n, icol, 1, 1, ierr)
            call permute(rng, m, irow, ierr)
            call permute(rng, n, icol, ierr)

        do i = 1, niter_ap
            call rng%iuniform(2, ifoo, 1, m, ierr)
            call nn_lrma_maxvol(method, m, n, r, U, m, VT, min(m,n), p, irow, icol, niter_maxvol,  &
                thresh, work, lwork, iwork, liwork, ierr, shift_ratio, ifoo(1), ifoo(2))
            call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
            call log_performance(i, m, n, A, AA, nrmf, nrmc)
        end do
    end subroutine run_ap_maxvol

    subroutine run_ap_maxvol_proj &
    (rng, m, n, A, r, method, niter_ap, kr, kc, niter_maxvol, thresh, shift_ratio)
    use maria_utils_mod,        only:   &
        arange,                         &
        permute
    use maria_la_core_mod,      only:   &
        dgenrmf,                        &
        dgenrmc,                        &
        dlacpy,                         &
        ddgmm,                          &
        dgesdd_q
    use maria_lr_la_mod,        only:   &
        dlr2full
    use nonnegative_matrix_mod,        only:   &
        nn_lrma_maxvol_proj

        class(prng), intent(in)              :: rng
        integer,    intent(in)              :: m
        integer,    intent(in)              :: n
        real(WP),   intent(in), contiguous  :: A(:)
        integer,    intent(in)              :: r
        integer,    intent(in)              :: method
        integer,    intent(in)              :: niter_ap
        integer,    intent(in)              :: kr
        integer,    intent(in)              :: kc
        integer,    intent(in)              :: niter_maxvol
        real(WP),   intent(in)              :: thresh
        real(WP),   intent(in)              :: shift_ratio

        integer     :: ifoo(2), lwork, liwork, ierr, i
        real(WP)    :: foo(1), nrmf, nrmc

        integer,    allocatable :: iwork(:), icol(:), irow(:), icol_short(:), irow_short(:)
        real(WP),   allocatable :: AA(:), S(:), U(:), VT(:), work(:)

        ! Estimate workspace
        call dgesdd_q('s', m, m, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call nn_lrma_maxvol_proj(method, m, n, r, foo, m, foo, min(m,n), kr, ifoo, ifoo, kc, ifoo, ifoo, &
            niter_maxvol, thresh, foo, -1, ifoo, -1, ierr, shift_ratio)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(AA(m*n), S(min(m,n)), U(m*min(m,n)), VT(min(m,n)*n), work(lwork), iwork(liwork),   &
            irow(m), icol(n), irow_short(m), icol_short(n))

        ! Norm of A
        nrmf = dgenrmf(m, n, A, m, ierr)
        nrmc = dgenrmc(m, n, A, m, ierr)

        ! Best rank-r approximation
        call dlacpy('a', m, n, A, m, AA, m)
        call dgesdd_q('s', m, n, AA, m, S, U, m, VT, min(m,n), work, lwork, iwork, liwork, ierr)
        call ddgmm('r', m, r, U, m, S, 1, ierr)

        ! Initial performance
        call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
        call log_performance(0, m, n, A, AA, nrmf, nrmc)

            call arange(m, irow, 1, 1, ierr)
            call arange(n, icol, 1, 1, ierr)
            call permute(rng, m, irow, ierr)
            call permute(rng, n, icol, ierr)
            call arange(m, irow_short, 1, 1, ierr)
            call arange(n, icol_short, 1, 1, ierr)
            call permute(rng, m, irow_short, ierr)
            call permute(rng, n, icol_short, ierr)

        do i = 1, niter_ap

            call rng%iuniform(2, ifoo, 1, m, ierr)
            call nn_lrma_maxvol_proj(method, m, n, r, U, m, VT, min(m,n), kr, irow, icol_short, &
                kc, icol, irow_short, niter_maxvol, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, ifoo(1), ifoo(2))
            call dlr2full(m, n, r, U, m, VT, min(m,n), AA, m, ierr)
            call log_performance(i, m, n, A, AA, nrmf, nrmc)
        end do
    end subroutine run_ap_maxvol_proj

    subroutine log_performance &
    (iter, m, n, A, AA, nrmf, nrmc)
    use maria_la_core_mod,  only:   &
        dgenrmf, dgenrmc, daxpy
    use nonnegative_matrix_mod,    only:   &
        distance_to_nn

        integer,    intent(in   )               :: iter
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        real(WP),   intent(in   ),  contiguous  :: A(:)
        real(WP),   intent(inout),  contiguous  :: AA(:)
        real(WP),   intent(in   )               :: nrmf
        real(WP),   intent(in   )               :: nrmc

        real(WP)    :: rerrf, rerrc, distc, distf

        call distance_to_nn(m, n, AA, m, distc, distf)
        print '( (I0), (1pe12.5) )', iter, distf / nrmf
        print '( (I0), (1pe12.5) )', iter, distc / nrmc

        call daxpy(m*n, -ONE, A, 1, AA, 1)
        rerrf = dgenrmf(m, n, AA, m, ierr) / nrmf
        rerrc = dgenrmc(m, n, AA, m, ierr) / nrmc
        print '( (I0), (1pe12.5) )', iter, rerrf
        print '( (I0), (1pe12.5) )', iter, rerrc
        flush(6)
    end subroutine log_performance
end program matrix_smolukh
