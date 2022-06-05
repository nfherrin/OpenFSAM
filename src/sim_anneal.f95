!simmulated annealing functions
MODULE sim_anneal
  USE globals
  USE travel_sales
  IMPLICIT NONE
  PRIVATE
  PUBLIC sa_comb_type,sa_cont_type

  REAL(8) :: alpha=0.0001,t_current,t_max=100,t_min=0,e_current,e_best
  INTEGER,ALLOCATABLE :: s_current(:),s_neigh(:)
  INTEGER :: step

!the base simulated annealing solver type
  TYPE,ABSTRACT :: sa_type_base
    !number of points in a state
    INTEGER :: size_states=0
    !maximum number of SA iterations
    INTEGER :: max_step=100
    !alpha value for cooling
    REAL(8) :: alpha=0.01
    !maximum temperature
    REAL(8) :: t_max=100
    !minimum temperature
    REAL(8) :: t_min=0
    !cooling option
    CHARACTER(64) :: cool_opt=''
    !Cooling schedule
    PROCEDURE(prototype_cooling),POINTER :: cool => NULL()
    CONTAINS
      !optimization subroutine using simulated annealing
      PROCEDURE :: optimize
  ENDTYPE sa_type_base

!combinatorial simulated annealing type
  TYPE,EXTENDS(sa_type_base) :: sa_comb_type
    !combinatorial problem state array (for perturbing combinatorial problems)
    INTEGER, POINTER, DIMENSION(:) :: state_cur
    !combinatorial problem neighbor array after pertubation
    INTEGER, ALLOCATABLE, DIMENSION(:) :: state_neigh
    !energy calculation
    PROCEDURE(prototype_eg_comb),POINTER :: energy => NULL()
  ENDTYPE sa_comb_type

!continuous simulated annealing type
  TYPE,EXTENDS(sa_type_base) :: sa_cont_type
    !continuous problem state array (for perturbing continuous problem values)
    REAL, POINTER, DIMENSION(:) :: state_cur
    !combinatorial problem neighbor array after pertubation
    REAL(8), ALLOCATABLE, DIMENSION(:) :: state_neigh
    !energy calculation
    PROCEDURE(prototype_eg_cont),POINTER :: energy => NULL()
  ENDTYPE sa_cont_type

!Simple abstract interface for a combinatorial energy computation subroutine
  ABSTRACT INTERFACE
    FUNCTION prototype_eg_comb(thisSA,state_val)
      IMPORT :: sa_comb_type
      CLASS(sa_comb_type),INTENT(INOUT) :: thisSA
      INTEGER,DIMENSION(:),INTENT(IN) :: state_val
      REAL(8) :: prototype_eg_comb
    ENDFUNCTION prototype_eg_comb
  ENDINTERFACE

!Simple abstract interface for a continuous energy computation subroutine
  ABSTRACT INTERFACE
    FUNCTION prototype_eg_cont(thisSA,state_val)
      IMPORT :: sa_cont_type
      CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
      REAL(8),DIMENSION(:),INTENT(IN) :: state_val
      REAL(8) :: prototype_eg_cont
    ENDFUNCTION prototype_eg_cont
  ENDINTERFACE

!Simple abstract interface for a continuous energy computation subroutine
  ABSTRACT INTERFACE
    FUNCTION prototype_cooling(thisSA,tmin,tmax,alpha,k,n)
      IMPORT :: sa_type_base
      CLASS(sa_type_base),INTENT(INOUT) :: thisSA
      REAL(8),INTENT(IN) :: tmin,tmax,alpha
      INTEGER,INTENT(IN) :: k,n
      REAL(8) :: prototype_cooling
    ENDFUNCTION prototype_cooling
  ENDINTERFACE

CONTAINS

  !simulate annealing for the traveling salesman
  SUBROUTINE optimize(thisSA)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8) :: e_neigh
    INTEGER,ALLOCATABLE :: s_best(:)
    INTEGER :: i,step_actual
    REAL(8) :: temp_r,start,finish

    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'performing simulated annealing'
    !write(*,*)'**********************************************************************************'

    CALL set_cooling(thisSA)

    !allocate the neighbor state variables and set the cooling
    SELECTTYPE(thisSA)
      TYPEIS(sa_comb_type)
        ALLOCATE(thisSA%state_neigh(thisSA%size_states))
        thisSA%state_neigh=thisSA%state_cur
      TYPEIS(sa_cont_type)
        ALLOCATE(thisSA%state_neigh(thisSA%size_states))
        thisSA%state_neigh=thisSA%state_cur
    ENDSELECT

    ALLOCATE(s_current(num_customers))
    ALLOCATE(s_neigh(num_customers))
    ALLOCATE(s_best(num_customers))

    !set the initial conditions
    t_current=t_max
    DO i=1,num_customers
      s_current(i)=i
    ENDDO
    s_best=s_current

    e_current=path_len(s_current)
    e_best=e_current

    alpha=0.85
    step=0
    step_actual=0
    CALL CPU_TIME(start)
    DO WHILE(step .LE. thisSA%max_step .AND. t_current .GE. t_min)
      step=step+1
      step_actual=step_actual+1
      !get a new neighbor and compute energy
      CALL get_neigh_comb()
      SELECTTYPE(thisSA)
        TYPEIS(sa_comb_type)
          e_neigh=thisSA%energy(thisSA%state_neigh)
        TYPEIS(sa_cont_type)
          e_neigh=thisSA%energy(thisSA%state_neigh)
      ENDSELECT
      e_neigh=path_len(s_neigh)
      !check and see if we accept the new temperature (lower temps alwasy accepted)
      CALL random_number(temp_r)
      IF(temp_r .LE. accept_prob(e_current,e_neigh))THEN
        s_current=s_neigh
        e_current=e_neigh
      ELSE
        !otherwise reject the neighbor and do nothing
      ENDIF
      !cool the temperature
      CALL temp_cool()
      !if it is the best energy, it's our new best value
      IF(e_current .LT. e_best)THEN
        e_best=e_current
        s_best=s_current
        !IF(step .GT. 100)step=step-100
        !write(*,*)step_actual,step,e_best
      ENDIF
      !IF(MOD(step,1000) .EQ. 0)THEN
      !  e_current=e_best
      !  s_current=s_best
      !ENDIF
      !WRITE(*,*)t_current
    ENDDO
    CALL CPU_TIME(finish)

    WRITE(*,'(A,ES16.8)',ADVANCE='NO')'iterations',(step_actual-1)*1.0
    !WRITE(*,'(A,ES16.8)')'SA Minimum path length:',e_best
    sa_best=e_best
    !WRITE(*,'(A,1000I6)')'Minimum path:',s_best
    !WRITE(*,'(A,ES16.8,A)')'Simulated annealing finished in: ',finish-start,' seconds'
    DEALLOCATE(s_current,s_neigh,s_best)
    thisSA%cool => NULL()
  ENDSUBROUTINE optimize

  !get a new neighbor state
  SUBROUTINE get_neigh_comb()
    INTEGER :: j1,j2,temp_i
    REAL(8) :: temp_r
    s_neigh=s_current
    !get switch indeces
    DO
      CALL random_number(temp_r)
      j1=1+FLOOR(num_customers*temp_r)
      CALL random_number(temp_r)
      j2=1+FLOOR(num_customers*temp_r)
      IF(j1 .NE. j2)EXIT
    ENDDO
    !set the new neighbor
    temp_i=s_neigh(j1)
    s_neigh(j1)=s_neigh(j2)
    s_neigh(j2)=temp_i
  ENDSUBROUTINE get_neigh_comb

  !function for the acceptance probability
  FUNCTION accept_prob(e_current,e_neigh)
    REAL(8),INTENT(IN) :: e_current,e_neigh
    REAL(8) :: accept_prob
    REAL(8) :: delta_e

    delta_e=e_neigh-e_current
    IF(-delta_e/t_current .LE. -700)THEN
      accept_prob=0.0
    ELSEIF(-delta_e/t_current .GE. 700)THEN
      accept_prob=10.0
    ELSE
      accept_prob=EXP(-delta_e/t_current)
    ENDIF
  ENDFUNCTION accept_prob

  !non-monotonic adaptive cooling schedule
  SUBROUTINE temp_cool()
    REAL(8) :: mu
    !compute the linear value
    t_current=t_max/(1.0+alpha*step**2)
    !now the non-monotonic factorization
    mu=1.0+(e_current-e_best)/e_current
    t_current=mu*t_current
  ENDSUBROUTINE temp_cool

  !set the cooling schedule
  SUBROUTINE set_cooling(thisSA)
    CLASS(sa_type_base) :: thisSA

    SELECTCASE(thisSA%cool_opt)
      CASE('LinMult')
        thisSA%cool => lin_mult_cool
      CASE('ExpMult')
        thisSA%cool => exp_mult_cool
      CASE('LogMult')
        thisSA%cool => log_mult_cool
      CASE('QuadMult')
        thisSA%cool => quad_mult_cool
    ENDSELECT
  ENDSUBROUTINE set_cooling

  !natural log exponential multiplicative cooling
  FUNCTION exp_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: exp_mult_cool

    exp_mult_cool=tmax*alpha**k
  ENDFUNCTION exp_mult_cool

  !linear multiplicative cooling
  FUNCTION lin_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: lin_mult_cool

    lin_mult_cool=tmax-alpha*k
  ENDFUNCTION lin_mult_cool

  !linear multiplicative cooling
  FUNCTION log_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: log_mult_cool

    log_mult_cool=tmax/(1.0+alpha*LOG(k+1.0))
  ENDFUNCTION log_mult_cool

  !linear multiplicative cooling
  FUNCTION quad_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: quad_mult_cool

    quad_mult_cool=tmax/(1.0+alpha*k**2)
  ENDFUNCTION quad_mult_cool
END MODULE sim_anneal
