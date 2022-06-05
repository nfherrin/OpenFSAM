!A Fortran based neural network accelerated simulated annealing software.
PROGRAM fortNNASA
  USE globals
  USE infuncs
  USE outfuncs
  USE sim_anneal
  USE travel_sales
  IMPLICIT NONE

  INTEGER :: i

  CALL get_cmdargs()

  CALL ts_init()

  !sets up the simulated annealing for the traveling salesman problem
  CALL setup_ts_sa()

  !CALL simulate_anneal()

  IF(sort_best .GT. 1.0E-10)WRITE(*,'(A,ES16.8)')'  Error: ',ABS(sort_best-sa_best)/sort_best

END PROGRAM fortNNASA
