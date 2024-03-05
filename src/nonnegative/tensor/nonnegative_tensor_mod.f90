module nonnegative_tensor_mod
implicit none (type, external)

interface
!----------------------------------------------------------------------------------------------------------------------------

module subroutine distance_to_nn &
(d, n, A, distc, distf)
use maria_kinds_mod,    only:   &
    WP  => DP
implicit none
    integer,    intent(in   )               :: d
    integer,    intent(in   ),  contiguous  :: n(:)
    real(WP),   intent(in   ),  contiguous  :: A(:)
    real(WP),   intent(  out)               :: distc
    real(WP),   intent(  out)               :: distf
end subroutine distance_to_nn

!----------------------------------------------------------------------------------------------------------------------------

module subroutine nn_lrta_ttsvd &
(method, d, n, r, cores, work, lwork, iwork, liwork, ierr, shift_ratio)
use maria_kinds_mod,    only:   &
    WP  =>  DP
use maria_arr_mod,      only:   & 
    AR  =>  darr
implicit none
    integer,    intent(in   )               :: method
    integer,    intent(in   )               :: d
    integer,    intent(in   ), contiguous   :: n(:)
    integer,    intent(in   ), contiguous   :: r(:)
    type(AR),   intent(inout)               :: cores(:)
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
    real(WP),   intent(in   ), optional     :: shift_ratio
end subroutine nn_lrta_ttsvd

!----------------------------------------------------------------------------------------------------------------------------
end interface
end module nonnegative_tensor_mod
