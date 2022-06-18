!simulated annealing of optimal simulated annealing of traveling salesman problem
MODULE sa_ts_sa
  USE globals
  USE travel_sales
  IMPLICIT NONE
  PRIVATE
  PUBLIC setup_sa_ts_sa, state_transform
CONTAINS

  !sets up the traveling salesman annealing problem
  SUBROUTINE setup_sa_ts_sa()
    sa_ts_simanneal%max_step=1000
    sa_ts_simanneal%t_max=100.0D0
    sa_ts_simanneal%t_min=1.0D-14
    sa_ts_simanneal%cool_opt='QuadAdd'
    sa_ts_simanneal%mon_cool=.FALSE.
    sa_ts_simanneal%prog_bar=.TRUE.
    !max and min values of the state variables
    sa_ts_simanneal%smax=2.0D0
    sa_ts_simanneal%smin=0.0D0
    sa_ts_simanneal%damping=0.0D0
    sa_ts_simanneal%resvar=50.0D0
    sa_ts_simanneal%damp_dyn=.TRUE.
    sa_ts_simanneal%num_perturb=1
    ALLOCATE(sa_ts_simanneal%state_curr(3))
    !all state variables use a functional transform to actually get
    !their values. This is so they can all use the same damping and max/min
    !start with ok initial guesses
    sa_ts_simanneal%state_curr(1)=1.0
    sa_ts_simanneal%state_curr(2)=1.0
    sa_ts_simanneal%state_curr(3)=1.0D-3

    !the energy function
    sa_ts_simanneal%energy => sa_ts_eg
  ENDSUBROUTINE setup_sa_ts_sa

  !this function runs a simulated annealing problem for the traveling salesman with state determined
  !temperatures and alpha values
  FUNCTION sa_ts_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: sa_ts_eg

    INTEGER :: i
    REAL(8) :: l2err,l2its

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)i=thisSA%size_states

    l2err=0.0
    l2its=0.0
    !you want to do it many times so that your energy isn't too dependent on randomness
    DO i=1,100
      CALL ts_init()
      ALLOCATE(ts_simanneal)
      !setup the simulated annealing with the new parameters
      CALL setup_ts_sa(state_transform(state_ord(:)))
      CALL ts_simanneal%optimize()
      sa_best=ts_simanneal%e_best
      !the L2 norm of the relative error is our first variable for the energy calculation
      l2err=l2err+(ABS(sa_best-sort_best)/sort_best)**2
      !the second variable is L2 norm of the number of iterations
      l2its=l2its+(ts_simanneal%total_steps*1.0D0)**2
      DEALLOCATE(ts_simanneal)
    ENDDO
    l2err=SQRT(l2err/(i-1))
    l2its=SQRT(l2its/(i-1))

    !heavily weight any error, we want the optimal solution. Less heavy weight for the
    ! number of iterations
    sa_ts_eg=(l2err*1.0D+8+l2its*1.0D-4)*1.0D-1

  ENDFUNCTION sa_ts_eg

  !state variable transformation to our actual parameters
  FUNCTION state_transform(state_vars)
    REAL(8),INTENT(IN) :: state_vars(3)
    REAL(8) :: state_transform(3)

    !num iterations doesn't go above 10000 or below 100
    state_transform(1)=4.95D+3*state_vars(1)+1.0D+2
    !maximum temp doesn't go above 200 or below 11
    state_transform(2)=94.5D0*state_vars(2)+11.D0
    !minimum temp doesn't go above 10 or below 1.0E-12
    state_transform(3)=10.0D0*state_vars(3)+1.0D-12
  ENDFUNCTION state_transform
END MODULE sa_ts_sa
