!globals module
MODULE globals
  IMPLICIT NONE

  !customer locations
  REAL(8),ALLOCATABLE :: cust_locs(:,:)

  !number of customers
  INTEGER :: num_customers=4
CONTAINS
END MODULE globals
