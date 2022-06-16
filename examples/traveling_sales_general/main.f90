!Simulated annealing solution to the traveling salesman problem.
PROGRAM fortNNASA
  USE globals
  USE infuncs
  USE outfuncs
  USE travel_sales
  IMPLICIT NONE

  CALL get_cmdargs()

  CALL ts_init()

  !sets up the simulated annealing for the traveling salesman problem
  CALL setup_ts_sa()


  write(*,'(A)')'**********************************************************************************'
  write(*,'(A)')'**********************************************************************************'
  write(*,'(A)')'**********************************************************************************'
  write(*,'(A)')'**********************************************************************************'
  write(*,'(A)')'**********************************************************************************'
  write(*,'(A)')'                     performing simulated annealing'
  write(*,'(A)')'**********************************************************************************'
  CALL ts_simanneal%optimize()
  WRITE(*,'(A,ES16.8)')'final length: ',ts_simanneal%e_best
  WRITE(*,'(A,ES16.8)')'number of steps: ',ts_simanneal%total_steps*1.0
  sa_best=ts_simanneal%e_best

  IF(sort_best .GT. 1.0E-10)WRITE(*,'(A,ES16.8)')'Error: ',ABS(sort_best-sa_best)/sort_best

END PROGRAM fortNNASA
