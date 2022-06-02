!A Fortran based neural network accelerated simulated annealing software.
PROGRAM fortNNASA
  USE globals
  USE infuncs
  USE outfuncs
  USE sim_anneal
  USE travel_sales
  IMPLICIT NONE

  INTEGER :: i

  DO i=1,1000
    CALL ts_init()

    CALL simulate_anneal()

    IF(sort_best .GT. 1.0E-10)WRITE(*,'(A,ES16.8)')'  Error: ',ABS(sort_best-sa_best)/sort_best

    DEALLOCATE(cust_locs)
  ENDDO

END PROGRAM fortNNASA
