module nonnegative_matrix_mod
implicit none (type, external)

interface
!----------------------------------------------------------------------------------------------------------------------------

module subroutine distance_to_nn &
(m, n, A, lda, distc, distf)
use maria_kinds_mod,    only:   &
    WP  => DP
implicit none
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ),  contiguous  :: A(:)
    integer,    intent(in   )               :: lda
    real(WP),   intent(  out)               :: distc
    real(WP),   intent(  out)               :: distf
end subroutine distance_to_nn

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_svd &
(method, m, n, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, ierr, shift_ratio)
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
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
end subroutine nn_lrma_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_rsvd &
(method, m, n, r, U, ldu, VT, ldvt, rng, p, q, work, lwork, iwork, liwork, ierr, shift_ratio)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
implicit none
    integer,        intent(in   )               :: method
    integer,        intent(in   )               :: m
    integer,        intent(in   )               :: n
    integer,        intent(in   )               :: r
    real(WP),       intent(inout), contiguous   :: U(:)
    integer,        intent(in   )               :: ldu
    real(WP),       intent(inout), contiguous   :: VT(:)
    integer,        intent(in   )               :: ldvt
    class(prng),    intent(in   )               :: rng
    integer,        intent(in   )               :: p
    integer,        intent(in   )               :: q
    real(WP),       intent(  out), contiguous   :: work(:)
    integer,        intent(in   )               :: lwork
    integer,        intent(  out), contiguous   :: iwork(:)
    integer,        intent(in   )               :: liwork
    integer,        intent(  out)               :: ierr
    real(WP),       intent(in   ), optional     :: shift_ratio
end subroutine nn_lrma_rsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol &
(method, m, n, r, U, ldu, VT, ldvt, p, irow, icol, niter, thresh, work, lwork, iwork, liwork, ierr, shift_ratio, i1, j1)
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
    integer,    intent(in   )               :: p
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
    integer,    intent(in   ), optional     :: i1
    integer,    intent(in   ), optional     :: j1
end subroutine nn_lrma_maxvol

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrma_maxvol_proj &
(method, m, n, r, U, ldu, VT, ldvt, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, &
    work, lwork, iwork, liwork, ierr, shift_ratio, i1, j1)
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
    integer,    intent(in   ), optional     :: i1
    integer,    intent(in   ), optional     :: j1
end subroutine nn_lrma_maxvol_proj

!----------------------------------------------------------------------------------------------------------------------------
end interface
end module nonnegative_matrix_mod
