program matrix_lr_randort
use maria_kinds_mod,        only:   &
    WP      => DP
use maria_constants_mod,    only:   &
    ZERO    => D_ZERO,              &
    ONE     => D_ONE
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only:   &
    prng    => prng_mkl
#else
use maria_prng_builtin_mod, only:   &
    prng    => prng_builtin
#endif
use maria_la_utils_mod,     only:   &
    drandort
use maria_la_core_mod,      only:   &
    dscal, dgenrmc
use maria_lr_la_mod,        only:   &
    dlr2full
use maximum_matrix_mod,     only:   &
    run_ap_svd, run_ap_rsvd

implicit none (type, external)

character(100)  :: seed_str
integer         :: n, r, rr, method_lrproj, ierr, seed, i, lwork, p, q
type(prng)      :: rng
real(WP)        :: foo(1), nrmA

real(WP),   allocatable :: A(:), U(:), VT(:), work(:)

! Read command-line arguments
i = 1
call get_command_argument(i, seed_str)
read(seed_str, *) seed
i = i + 1
! Data description -------------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) n
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) r
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) rr
i = i + 1
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
end if

allocate(U(n*rr), VT(rr*n), A(n*n))
call rng%init(seed+rr, ierr)

! Generate data ----------------------------------------------------------
call drandort(rng, n, rr, U, n, foo, -1, ierr)
lwork = int(foo(1))
call drandort(rng, rr, n, VT, rr, foo, -1, ierr)
lwork = max(lwork, int(foo(1)))
allocate( work(lwork) )
call drandort(rng, n, rr, U, n, work, lwork, ierr)
call drandort(rng, rr, n, VT, rr, work, lwork, ierr)
call dlr2full(n, n, rr, U, n, VT, rr, A, n, ierr)
nrmA = dgenrmc(n, n, A, n, ierr)
call dscal(n*n, ONE/nrmA, A, 1)
deallocate( work, U, VT )

! Compute its MAX_LRMA
if (method_lrproj == 1) then
    call run_ap_svd(n, n, A, r, rng, .true., 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)
elseif (method_lrproj == 2) then
    call run_ap_rsvd(n, n, A, r, p, q, rng, .true., 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)
end if

end program matrix_lr_randort
