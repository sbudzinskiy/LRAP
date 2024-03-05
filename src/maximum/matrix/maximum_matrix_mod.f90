module maximum_matrix_mod
implicit none (type, external)

interface
!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrma_svd &
(m, n, A, lda, r, U, ldu, VT, ldvt, tol, work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP  => DP
implicit none
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   )               :: lda
    integer,    intent(in   )               :: r
    real(WP),   intent(inout), contiguous   :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout), contiguous   :: VT(:)
    integer,    intent(in   )               :: ldvt
    real(WP),   intent(in   )               :: tol
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
end subroutine max_lrma_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_svd &
(m, n, A, r, rng, init_randort, c1, c2, c3, c4)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
implicit none
    integer,    intent(in)              :: m
    integer,    intent(in)              :: n
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in)              :: r
    class(prng),intent(in)              :: rng
    logical,    intent(in)              :: init_randort
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4
end subroutine run_ap_svd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine max_lrma_rsvd &
(m, n, A, lda, r, U, ldu, VT, ldvt, tol, rng, p, q, work, lwork, iwork, liwork, ierr)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
implicit none
    integer,    intent(in   )               :: m
    integer,    intent(in   )               :: n
    real(WP),   intent(in   ), contiguous   :: A(:)
    integer,    intent(in   )               :: lda
    integer,    intent(in   )               :: r
    real(WP),   intent(inout), contiguous   :: U(:)
    integer,    intent(in   )               :: ldu
    real(WP),   intent(inout), contiguous   :: VT(:)
    integer,    intent(in   )               :: ldvt
    real(WP),   intent(in   )               :: tol
    class(prng),intent(in   )               :: rng
    integer,    intent(in   )               :: p
    integer,    intent(in   )               :: q
    real(WP),   intent(  out), contiguous   :: work(:)
    integer,    intent(in   )               :: lwork
    integer,    intent(  out), contiguous   :: iwork(:)
    integer,    intent(in   )               :: liwork
    integer,    intent(  out)               :: ierr
end subroutine max_lrma_rsvd

!----------------------------------------------------------------------------------------------------------------------------

module subroutine run_ap_rsvd &
(m, n, A, r, p, q, rng, init_randort, c1, c2, c3, c4)
use maria_kinds_mod,    only:   &
    WP  => DP
use maria_prng_mod,     only:   &
    prng
implicit none
    integer,    intent(in)              :: m
    integer,    intent(in)              :: n
    real(WP),   intent(in), contiguous  :: A(:)
    integer,    intent(in)              :: r
    integer,    intent(in)              :: p
    integer,    intent(in)              :: q
    class(prng),intent(in)              :: rng
    logical,    intent(in)              :: init_randort
    real(WP),   intent(in)              :: c1
    real(WP),   intent(in)              :: c2
    real(WP),   intent(in)              :: c3
    real(WP),   intent(in)              :: c4
end subroutine run_ap_rsvd

!----------------------------------------------------------------------------------------------------------------------------
end interface
end module maximum_matrix_mod
