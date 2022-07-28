!globals module
MODULE globals
  USE OpenFSAM
  IMPLICIT NONE

  !customer locations
  REAL(8),ALLOCATABLE :: cust_locs(:,:)

  !number of customers
  INTEGER :: num_customers=10

  !traveling salesman problem dimensions
  INTEGER :: prob_dim=1

  !best from sa and sort
  REAL(8) :: sa_best,sort_best

  !simulated annealing object for the traveling salesman problem
  TYPE(sa_comb_type) :: ts_simanneal
CONTAINS
END MODULE globals
