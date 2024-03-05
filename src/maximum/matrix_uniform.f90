program matrix_uniform
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
use maximum_matrix_mod,     only:   &
    run_ap_svd, run_ap_rsvd

implicit none (type, external)

character(100)  :: seed_str
integer         :: n, r, method_lrproj, ierr, seed, i, p, q
type(prng)      :: rng

real(WP),   allocatable :: A(:)

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

allocate(A(n*n))
call rng%init(seed, ierr)

! Generate data ----------------------------------------------------------
call rng%duniform(n*n, A, -ONE, ONE, ierr)

! Compute its MAX_LRMA
if (method_lrproj == 1) then
    call run_ap_svd(n, n, A, r, rng, .false., 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)
elseif (method_lrproj == 2) then
    call run_ap_rsvd(n, n, A, r, p, q, rng, .false., 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)
end if

end program matrix_uniform
