!A simulated annealing for simulated annealing parameters.
PROGRAM sinannealing_sa
  USE globals
  USE infuncs
  USE outfuncs
  USE sa_ts_sa
  IMPLICIT NONE

  CALL get_cmdargs()
  !customer locations array, will change at each SA iteration
  ALLOCATE(cust_locs(num_customers,prob_dim))

  !sets up the simulated annealing for the traveling salesman problem
  CALL setup_sa_ts_sa()

  WRITE(*,'(A,40ES12.3)')'pre-optimized state',sa_ts_simanneal%state_curr(:)
  WRITE(*,'(A,40ES12.3)')'pre-optimized transformed state',state_transform(sa_ts_simanneal%state_curr(:))
  WRITE(*,'(A,40ES12.3)')'pre-optimized energy',sa_ts_simanneal%energy(sa_ts_simanneal%state_curr(:))

  WRITE(*,'(A)')'Annealing:'
  CALL sa_ts_simanneal%optimize()

  WRITE(*,'(A,40ES12.3)')'post-optimized state',sa_ts_simanneal%state_curr(:)
  WRITE(*,'(A,40ES12.3)')'post-optimized transformed state',state_transform(sa_ts_simanneal%state_curr(:))
  WRITE(*,'(A,40ES12.3)')'post-optimized energy',sa_ts_simanneal%e_best

ENDPROGRAM sinannealing_sa
