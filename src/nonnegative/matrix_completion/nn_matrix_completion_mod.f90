module nn_matrix_completion_mod
implicit none (type, external)

interface
!----------------------------------------------------------------------------------------------------------------------------

module subroutine rgd_completion_step(m, n, r, S, U, ldu, VT, ldvt, &
        sample_size, sample, sample_indices, &
        work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP => DP
implicit none
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    integer,    intent(in   )               :: r
    real(WP),   intent(inout),  contiguous  :: S(:)
    real(WP),   intent(inout),  contiguous  :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout),  contiguous  :: VT(:)
    integer,    intent(in   )               :: ldvt
    integer,    intent(in   )               :: sample_size
    real(WP),   intent(in   ),  contiguous  :: sample(:)
    integer,    intent(in   ),  contiguous  :: sample_indices(:)
    real(WP),   intent(  out),  contiguous  :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out),  contiguous  :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
end subroutine rgd_completion_step

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol &
(method, m, n, r, U, ldu, VT, ldvt, irow, icol, niter, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, neg_pos)
use maria_kinds_mod,    only:   &
    WP  => DP
implicit none
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    integer,    intent(in   )               :: r
    real(WP),   intent(inout), contiguous   :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout), contiguous   :: VT(:)
    integer,    intent(in   )               :: ldvt
    integer,    intent(inout), contiguous   :: irow(:)
    integer,    intent(inout), contiguous   :: icol(:)
    integer,    intent(in   )               :: niter
    real(WP),   intent(in   )               :: thresh
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
    integer,    intent(inout), optional     :: neg_pos(2)
end subroutine nn_lrma_maxvol

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol_proj &
(method, m, n, r, U, ldu, VT, ldvt, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, &
    work, lwork, iwork, liwork, ierr, shift_ratio, neg_pos)
use maria_kinds_mod,    only:   &
    WP  => DP
implicit none
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    integer,    intent(in   )               :: r
    real(WP),   intent(inout), contiguous   :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout), contiguous   :: VT(:)
    integer,    intent(in   )               :: ldvt
    integer,    intent(in   )               :: kr
    integer,    intent(inout), contiguous   :: irow(:)
    integer,    intent(inout), contiguous   :: icol_short(:)
    integer,    intent(in   )               :: kc
    integer,    intent(inout), contiguous   :: icol(:)
    integer,    intent(inout), contiguous   :: irow_short(:)
    integer,    intent(in   )               :: niter
    real(WP),   intent(in   )               :: thresh
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
    integer,    intent(inout), optional     :: neg_pos(2)
end subroutine nn_lrma_maxvol_proj

!----------------------------------------------------------------------------------------------------------------------------
end interface
end module nn_matrix_completion_mod
