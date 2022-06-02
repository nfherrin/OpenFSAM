!globals module
MODULE globals
  IMPLICIT NONE

  !customer locations
  REAL(8),ALLOCATABLE :: cust_locs(:,:)

  !number of customers
  INTEGER :: num_customers=10

  !traveling salesman problem dimensions
  INTEGER :: prob_dim=2
CONTAINS
END MODULE globals
