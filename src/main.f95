!A Fortran based neural network accelerated simulated annealing software.
PROGRAM fortNNASA
  USE globals
  USE infuncs
  USE outfuncs
  USE sim_anneal
  USE travel_sales
  IMPLICIT NONE

  CALL ts_init()

END PROGRAM fortNNASA
