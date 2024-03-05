program matrix_cmplt
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
use maria_la_core_mod,  only:   &
    dscal, dlacpy
use maria_lr_la_mod,    only:   &
    dlrval
implicit none (type, external)

character(100)  :: seed_str
integer         :: m, n, r, k, method_ap, method_lrproj, ierr, seed, niter_ap, niter_maxvol, i, sample_size, niter_rgd
type(prng)      :: rng
real(WP)        :: thresh, shift_ratio, oversample

integer,    allocatable :: sample_indices(:)
real(WP),   allocatable :: U(:), VT(:), sample(:)

! Read command-line arguments
i = 1
call get_command_argument(i, seed_str)
read(seed_str, *) seed
i = i + 1
! Data description -------------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) m
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) n
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) r
i = i + 1
! Completion -------------------------------------------------------------------
call get_command_argument(i, seed_str)
read(seed_str, *) oversample
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) niter_rgd
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
call get_command_argument(i, seed_str)
read(seed_str, *) niter_maxvol
i = i + 1
call get_command_argument(i, seed_str)
read(seed_str, *) thresh
i = i + 1
if (method_lrproj == 2) then
    call get_command_argument(i, seed_str)
    read(seed_str, *) k
    i = i + 1
end if

seed = ceiling(seed * oversample)
call rng%init(seed, ierr)

! Generate nonnegative low-rank matrix  ----------------------------------------
allocate(U(m*r), VT(r*n))
call rng%duniform(m*r, U, ZERO, ONE, ierr)
call rng%duniform(r*n, VT, ZERO, ONE, ierr)
call dscal(m*r, sqrt(ONE / m), U, 1)
call dscal(r*n, sqrt(ONE / n), VT, 1)

! Generate sample --------------------------------------------------------------
sample_size = floor(r * (m + n - r) * oversample)
allocate( sample(sample_size), sample_indices(2*sample_size) )
do i = 1, sample_size
associate( irow => sample_indices(2*i-1:), icol => sample_indices(2*i:) )
    call rng%iuniform(1, irow, 1, m, ierr)
    call rng%iuniform(1, icol, 1, n, ierr)
!    irow(1) = mod(i-1, n) + 1
!    icol(1) = (i - 1) / n + 1
    sample(i) = dlrval(m, n, irow(1), icol(1), r, U, m, VT, r, ierr)
end associate
end do

! Compute its NN_LRMA
if (method_lrproj == 1) then
    call nn_rgd_maxvol(rng, m, n, r, U, VT, sample_size, sample, sample_indices, niter_rgd, method_ap, niter_ap, shift_ratio)
end if

contains
    
subroutine nn_rgd_maxvol &
(rng, m, n, r, U, VT, sample_size, sample, sample_indices, niter_rgd, method_ap, niter_ap, shift_ratio)
use maria_utils_mod,    only:   &
    arange,                     &
    permute
use maria_la_utils_mod, only:   &
    drandort
use maria_la_core_mod,  only:   &
    dlaset,                     &
    dlacpy,                     &
    ddgmm
use maria_lr_la_mod,    only:   &
    dlrnrmf,                    &
    dlrort,                     &
    dlrsvd_ort
use nn_matrix_completion_mod,   only:   &
    rgd_completion_step,                &
    nn_lrma_maxvol

    type(prng), intent(in)  :: rng
    integer,    intent(in)  :: m
    integer,    intent(in)  :: n
    integer,    intent(in)  :: r
    real(WP),   intent(in),  contiguous  :: U(:)
    real(WP),   intent(in),  contiguous  :: VT(:)
    integer,    intent(in)  :: sample_size
    real(WP),   intent(in), contiguous  :: sample(:)
    integer,    intent(in), contiguous  :: sample_indices(:)
    integer,    intent(in)  :: niter_rgd
    integer,    intent(in)  :: method_ap
    integer,    intent(in)  :: niter_ap
    real(WP),   intent(in)  :: shift_ratio

    integer     :: ifoo(1), iter, lwork, liwork, j, neg_pos(2)
    real(WP)    :: nrmA, foo(1)
    integer,    allocatable :: iwork(:), irow(:), icol(:)
    real(WP), allocatable :: S0(:), U0(:), VT0(:), work(:), UU(:), VVT(:), tau(:)

    call drandort(rng, m, r, foo, m, foo, -1, ierr)
    lwork = int(foo(1))
    call drandort(rng, r, n, foo, r, foo, -1, ierr)
    lwork = max(lwork, int(foo(1)))
    call rgd_completion_step(m, n, r, foo, foo, m, foo, r, sample_size, foo, ifoo, &
        foo, -1, ifoo, -1, ierr)
    lwork = max(lwork, int(foo(1)))
    liwork = ifoo(1)
    call nn_lrma_maxvol(method_ap, m, n, r, foo, m, foo, r, ifoo, ifoo, niter_maxvol, thresh, foo, -1, ifoo, -1, ierr, shift_ratio, neg_pos)
    lwork = max(lwork, int(foo(1)))
    liwork = max(liwork, ifoo(1))
    call dlrort('l', m, n, r, foo, m, foo, r, foo, foo, -1, ierr)
    lwork = max(lwork, int(foo(1)))
    call dlrsvd_ort('s', 'l', m, n, r, foo, m, foo, r, foo, foo, foo, m, foo, r, foo, -1, ifoo, -1, ierr)
    lwork = max(lwork, int(foo(1)))
    liwork = max(liwork, ifoo(1))

    allocate( work(lwork), iwork(liwork), irow(m), icol(n), tau(r) )

    nrmA = dlrnrmf(m, n, r, U, m, VT, r, work, lwork, ierr)
    
    ! Generate initial approximation as SVD
    allocate( S0(r), U0(m*r), VT0(r*n), UU(m*r), VVT(r*n) )
    call drandort(rng, m, r, U0, m, work, lwork, ierr)
    call drandort(rng, r, n, VT0, r, work, lwork, ierr)
    call dlaset('a', 1, r, ONE, ONE, S0, 1)
    call log_performance(0, m, n, r, U, VT, S0, U0, VT0, nrmA, work)
    
    ! Iterate
    call arange(m, irow, 1, 1, ierr)
    call arange(n, icol, 1, 1, ierr)
    call permute(rng, m, irow, ierr)
    call permute(rng, n, icol, ierr)
    do iter = 1, niter_rgd
        call rgd_completion_step(m, n, r, S0, U0, m, VT0, r, sample_size, sample, sample_indices, &
            work, lwork, iwork, liwork, ierr)

        ! Regularize
        call ddgmm('r', m, r, U0, m, S0, 1, ierr)
        do j = 1, niter_ap
            call rng%iuniform(2, neg_pos, 1, m, ierr)
            call nn_lrma_maxvol(method_ap, m, n, r, U0, m, VT0, r, irow, icol, niter_maxvol, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, neg_pos)
        end do
        call dlrort('l', m, n, r, U0, m, VT0, r, tau, work, lwork, ierr)
        call dlrsvd_ort('s', 'l', m, n, r, U0, m, VT0, r, tau, S0, UU, m, VVT, r, work, lwork, iwork, liwork, ierr)
        call dlacpy('a', m, r, UU, m, U0, m)
        call dlacpy('a', r, n, VVT, r, VT0, r)

        call log_performance(iter, m, n, r, U, VT, S0, U0, VT0, nrmA, work)
    end do
end subroutine nn_rgd_maxvol

    subroutine log_performance &
    (iter, m, n, r, U, VT, S0, U0, VT0, nrmf, work)
    use maria_la_core_mod,  only:   &
        dlacpy, ddgmm
    use maria_lr_la_mod,    only:   &
        dlrnrmf_diff

        integer,    intent(in   )               :: iter
        integer,    intent(in   )               :: m
        integer,    intent(in   )               :: n
        integer,    intent(in   )               :: r
        real(WP),   intent(in   ),  contiguous  :: U(:)
        real(WP),   intent(in   ),  contiguous  :: VT(:)
        real(WP),   intent(in   ),  contiguous  :: S0(:)
        real(WP),   intent(in   ),  contiguous  :: U0(:)
        real(WP),   intent(in   ),  contiguous  :: VT0(:)
        real(WP),   intent(in   )               :: nrmf
        real(WP),   intent(  out),  contiguous  :: work(:) ! m*r

        real(WP)    :: aerrf

        call dlacpy('a', m, r, U0, m, work, m)
        call ddgmm('r', m, r, work, m, S0, 1, ierr)
        associate(wrk => work(1 + m*r:))
            aerrf = dlrnrmf_diff(m, n, r, U, m, VT, r, r, work, m, VT0, r, wrk, size(work)-m*r, ierr)
        end associate
        print '( (I0), (1pe12.5) )', iter, aerrf / nrmf
        flush(6)
    end subroutine log_performance
end program matrix_cmplt
