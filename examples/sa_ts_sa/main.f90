!A simulated annealing for simulated annealing parameters.
PROGRAM fortNNASA
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

  WRITE(*,*)sa_ts_simanneal%energy(state_transform(sa_ts_simanneal%state_curr(:)))
  WRITE(*,*)sa_ts_simanneal%state_curr(:)
  WRITE(*,*)state_transform(sa_ts_simanneal%state_curr(:))

  CALL sa_ts_simanneal%optimize()

  WRITE(*,*)sa_ts_simanneal%energy(state_transform(sa_ts_simanneal%state_curr(:)))
  WRITE(*,*)sa_ts_simanneal%state_curr(:)
  WRITE(*,*)state_transform(sa_ts_simanneal%state_curr(:))

END PROGRAM fortNNASA
