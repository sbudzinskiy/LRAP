program tensor_uniform
use maria_kinds_mod,    only:   &
    WP  =>  DP
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
use maximum_tensor_mod,        only:   &
    run_ap_ttsvd
implicit none (type, external)

character(100)  :: seed_str
integer         :: d, maxn, maxr, ierr, seed, i, numel
type(prng)      :: rng

integer,    allocatable :: n(:), r(:)
real(WP),   allocatable :: A(:)

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
allocate(n(d), r(d-1), A(numel))
n = maxn
r = maxr
call rng%init(seed, ierr)

! Generate data ----------------------------------------------------------
call rng%duniform(numel, A, -ONE, ONE, ierr)

! Compute its MAX_LRTA
call run_ap_ttsvd(d, n, A, r, rng, 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)

end program tensor_uniform
