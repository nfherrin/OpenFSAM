!simulated annealing of optimal simulated annealing of traveling salesman problem
MODULE sa_ts_sa
  USE globals
  USE travel_sales
  IMPLICIT NONE
  PRIVATE
  PUBLIC setup_sa_ts_sa, state_transform
CONTAINS

  !sets up the traveling salesman problem
  SUBROUTINE setup_sa_ts_sa()
    sa_ts_simanneal%max_step=100
    sa_ts_simanneal%alpha=0.85
    sa_ts_simanneal%t_max=100
    sa_ts_simanneal%t_min=0
    sa_ts_simanneal%cool_opt='QuadMult'
    sa_ts_simanneal%mon_cool=.FALSE.
    ALLOCATE(sa_ts_simanneal%state_curr(3))
    !all state variables start at 1 and use a functional transform to actually get the values
    !this is so they can all use the same damping and max/mins
    sa_ts_simanneal%state_curr(1)=1.7
    sa_ts_simanneal%state_curr(2)=1.0
    sa_ts_simanneal%state_curr(3)=0.000001

    !the energy function
    sa_ts_simanneal%energy => sa_ts_eg
    !damping factor so changes are 0.5% of original values (1)
    sa_ts_simanneal%damping=0.01
    !max and min values of the state variables
    sa_ts_simanneal%smax=2.0D0-1.0D-12
    sa_ts_simanneal%smin=1.0D-12
  ENDSUBROUTINE setup_sa_ts_sa

  !this function runs a simulated annealing problem for the traveling salesman with state determined
  !temperatures and alpha values
  FUNCTION sa_ts_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: sa_ts_eg

    INTEGER :: i
    REAL(8) :: l2err,l2its

    l2err=0.0
    l2its=0.0
    !you want to do it multiple times so that your energy isn't too dependent on randomness
    DO i=1,10
      CALL ts_init()
      ALLOCATE(ts_simanneal)
      !setup the simulated annealing with the new parameters
      CALL setup_ts_sa(state_transform(sa_ts_simanneal%state_curr(:)))
      CALL ts_simanneal%optimize()
      sa_best=ts_simanneal%e_best
      !the L2 norm of the relative error is our first variable for the energy calculation
      l2err=l2err+(ABS(sa_best-sort_best)/sort_best)**2
      !the second variable is l2 norm of the number of iterations
      l2its=l2its+ts_simanneal%total_steps**2
      DEALLOCATE(ts_simanneal)
    ENDDO
    l2err=SQRT(l2err/(i-1))
    l2its=SQRT(l2its/(i-1))

    sa_ts_eg=l2err*1.0D+5+SQRT(l2its)

  ENDFUNCTION sa_ts_eg

  !state variable transformation to our actual parameters
  FUNCTION state_transform(state_vars)
    REAL(8),INTENT(IN) :: state_vars(3)
    REAL(8) :: state_transform(3)

    !alpha, starts at 0.5 and doesn't go above 1 or below 0
    state_transform(1)=0.5D0*state_vars(1)
    !maximum temp starts at 100 and doesn't go above 150 or below 50
    state_transform(2)=50.0D0*state_vars(2)+50.D0
    !minimum temp starts at 10 and doesn't go above 40 or below 0
    state_transform(3)=10.0D0*state_vars(3)**2
  ENDFUNCTION state_transform
END MODULE sa_ts_sa
