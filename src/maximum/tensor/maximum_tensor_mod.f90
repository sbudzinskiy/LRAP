module maximum_tensor_mod
implicit none (type, external)

interface
!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrta_ttsvd &
(d, n, A, r, cores, tol, work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP  =>  DP
use maria_arr_mod,      only:   & 
    AR  =>  darr
implicit none
    integer,    intent(in   )               :: d
    integer,    intent(in   ), contiguous   :: n(:)
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   ), contiguous   :: r(:)
    type(AR),   intent(inout)               :: cores(:)
    real(WP),   intent(in   )               :: tol
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
end subroutine max_lrta_ttsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_ttsvd &
(d, n, A, r, rng, c1, c2, c3, c4)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
implicit none
    integer,    intent(in)              :: d
    integer,    intent(in), contiguous  :: n(:)
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in), contiguous  :: r(:)
    class(prng),intent(in)              :: rng
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4
end subroutine run_ap_ttsvd

!----------------------------------------------------------------------------------------------------------------------------
end interface
end module maximum_tensor_mod
