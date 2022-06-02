!globals module
MODULE globals
  IMPLICIT NONE

  !customer locations
  REAL(8),ALLOCATABLE :: cust_locs(:,:)

  !number of customers
  INTEGER :: num_customers=40

  !traveling salesman problem dimensions
  INTEGER :: prob_dim=1

  !best from sa and sort
  REAL(8) :: sa_best,sort_best
CONTAINS
END MODULE globals
