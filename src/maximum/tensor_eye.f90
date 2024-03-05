program tensor_eye
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
use maria_la_core_mod,      only:   &
    dscal
use maximum_tensor_mod,        only:   &
    run_ap_ttsvd
implicit none (type, external)

character(100)  :: seed_str
integer         :: d, maxn, maxr, ierr, seed, i, numel, ind, j
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
call dscal(numel, ZERO, A, 1)
ind = 1
do i = 1, maxn
    ind = i
    do j = 2, d
        ind = ind + (i-1) * maxn**(j-1)
    end do
    A(ind) = ONE
end do

! Compute its MAX_LRTA
call run_ap_ttsvd(d, n, A, r, rng, 1e-2_WP, 1e-3_WP, 2.0_WP/3, 2e-2_WP)

end program tensor_eye
