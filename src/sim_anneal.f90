!simmulated annealing functions
MODULE sim_anneal
  IMPLICIT NONE
  PRIVATE
  PUBLIC sa_comb_type,sa_cont_type

  REAL(8),PARAMETER :: pi=4.D0*ATAN(1.D0)

!the base simulated annealing solver type
  TYPE,ABSTRACT :: sa_type_base
    !number of points in a state
    INTEGER :: size_states=0
    !maximum number of SA iterations
    INTEGER :: max_step=100
    !total number of steps it took
    INTEGER :: total_steps=0
    !alpha value for cooling
    REAL(8) :: alpha=0.01
    !maximum temperature
    REAL(8) :: t_max=100
    !minimum temperature
    REAL(8) :: t_min=0
    !best energy
    REAL(8) :: e_best=1.0D+307
    !cooling option
    CHARACTER(64) :: cool_opt='LinMult'
    !whether cooling is monotonic or not
    LOGICAL :: mon_cool=.TRUE.
    !Cooling schedule
    PROCEDURE(prototype_cooling),POINTER :: cool => NULL()
    CONTAINS
      !optimization subroutine using simulated annealing
      PROCEDURE :: optimize
  ENDTYPE sa_type_base

!combinatorial simulated annealing type
  TYPE,EXTENDS(sa_type_base) :: sa_comb_type
    !combinatorial problem state array (for perturbing combinatorial problems)
    INTEGER, POINTER, DIMENSION(:) :: state_curr
    !combinatorial problem neighbor array after pertubation
    INTEGER, ALLOCATABLE, DIMENSION(:) :: state_neigh
    !best energy state
    INTEGER, ALLOCATABLE, DIMENSION(:) :: state_best
    !energy calculation
    PROCEDURE(prototype_eg_comb),POINTER :: energy => NULL()
    CONTAINS
      !neighbor retrieval
      PROCEDURE,PASS :: get_neigh => get_neigh_comb
  ENDTYPE sa_comb_type

!continuous simulated annealing type
  TYPE,EXTENDS(sa_type_base) :: sa_cont_type
    !continuous problem state array (for perturbing continuous problem values)
    REAL(8), POINTER, DIMENSION(:) :: state_curr
    !continuous problem neighbor array after pertubation
    REAL(8), ALLOCATABLE, DIMENSION(:) :: state_neigh
    !best energy state
    REAL(8), ALLOCATABLE, DIMENSION(:) :: state_best
    !damping factor
    REAL(8) :: damping=1.0
    !upper and lower bounds, will be set to bounds of initial state if not changed
    REAL(8) :: smin=0.0,smax=0.0
    !energy calculation
    PROCEDURE(prototype_eg_cont),POINTER :: energy => NULL()
    CONTAINS
      !neighbor retrieval
      PROCEDURE,PASS :: get_neigh => get_neigh_cont
  ENDTYPE sa_cont_type

!Simple abstract interface for a combinatorial energy computation function
  ABSTRACT INTERFACE
    FUNCTION prototype_eg_comb(thisSA,state_val)
      IMPORT :: sa_comb_type
      CLASS(sa_comb_type),INTENT(INOUT) :: thisSA
      INTEGER,DIMENSION(:),INTENT(IN) :: state_val
      REAL(8) :: prototype_eg_comb
    ENDFUNCTION prototype_eg_comb
  ENDINTERFACE

!Simple abstract interface for a continuous energy computation function
  ABSTRACT INTERFACE
    FUNCTION prototype_eg_cont(thisSA,state_val)
      IMPORT :: sa_cont_type
      CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
      REAL(8),DIMENSION(:),INTENT(IN) :: state_val
      REAL(8) :: prototype_eg_cont
    ENDFUNCTION prototype_eg_cont
  ENDINTERFACE

!Simple abstract interface for a cooling function
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

  !optimize through simulated annealing
  SUBROUTINE optimize(thisSA)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8) :: e_neigh
    INTEGER :: step
    REAL(8) :: temp_r,start,finish,e_curr,t_curr

    CALL set_cooling(thisSA)

    !allocate the neighbor state variables and set the cooling, also set initial energy
    SELECTTYPE(thisSA)
      TYPEIS(sa_comb_type)
        thisSA%size_states=SIZE(thisSA%state_curr)
        IF(.NOT. ALLOCATED(thisSA%state_neigh))THEN
          ALLOCATE(thisSA%state_neigh(thisSA%size_states))
        ENDIF
        IF(.NOT. ALLOCATED(thisSA%state_best))THEN
          ALLOCATE(thisSA%state_best(thisSA%size_states))
        ENDIF
        thisSA%state_neigh=thisSA%state_curr
        thisSA%state_best=thisSA%state_curr
        !set energy to current energy
        e_curr=thisSA%energy(thisSA%state_curr)
        thisSA%e_best=e_curr
      TYPEIS(sa_cont_type)
        thisSA%size_states=SIZE(thisSA%state_curr)
        IF(.NOT. ALLOCATED(thisSA%state_neigh))THEN
          ALLOCATE(thisSA%state_neigh(thisSA%size_states))
        ENDIF
        IF(.NOT. ALLOCATED(thisSA%state_best))THEN
          ALLOCATE(thisSA%state_best(thisSA%size_states))
        ENDIF
        thisSA%state_neigh=thisSA%state_curr
        thisSA%state_best=thisSA%state_curr
        !set energy to current energy
        e_curr=thisSA%energy(thisSA%state_curr)
        thisSA%e_best=e_curr
        !set the bounds if not given
        IF(thisSA%smax-thisSA%smin .LT. 1.0D-13)THEN
          thisSA%smax=MAXVAL(thisSA%state_curr)
          thisSA%smin=MINVAL(thisSA%state_curr)
        ENDIF
    ENDSELECT

    t_curr=thisSA%t_max
    step=0
    CALL CPU_TIME(start)
    !actual simulated annealing happens here
    DO WHILE(step .LE. thisSA%max_step .AND. t_curr .GE. thisSA%t_min)
      step=step+1
      !get a new neighbor and compute energy
      SELECTTYPE(thisSA)
        TYPEIS(sa_comb_type)
          thisSA%state_neigh=thisSA%get_neigh(thisSA%state_curr)
          e_neigh=thisSA%energy(thisSA%state_neigh)
        TYPEIS(sa_cont_type)
          thisSA%state_neigh=thisSA%get_neigh(thisSA%state_curr,thisSA%damping,thisSA%smax,thisSA%smin)
          e_neigh=thisSA%energy(thisSA%state_neigh)
      ENDSELECT
      !check and see if we accept the new temperature (lower temps always accepted)
      CALL random_number(temp_r)
      IF(temp_r .LE. accept_prob(e_curr,e_neigh,t_curr))THEN
        SELECTTYPE(thisSA)
          TYPEIS(sa_comb_type)
            thisSA%state_curr=thisSA%state_neigh
          TYPEIS(sa_cont_type)
            thisSA%state_curr=thisSA%state_neigh
        ENDSELECT
        e_curr=e_neigh
      ELSE
        !otherwise reject the neighbor and do nothing
      ENDIF
      !cool the temperature
      t_curr=thisSA%cool(thisSA%t_min,thisSA%t_max,thisSA%alpha,step,thisSA%max_step)
      IF(.NOT. thisSA%mon_cool)t_curr=t_curr*(1.0+(e_curr-thisSA%e_best)/e_curr)
      !if it is the best energy, it's our new best value
      IF(e_curr .LT. thisSA%e_best)THEN
        thisSA%e_best=e_curr
        SELECTTYPE(thisSA)
          TYPEIS(sa_comb_type)
            thisSA%state_best=thisSA%state_neigh
          TYPEIS(sa_cont_type)
            thisSA%state_best=thisSA%state_neigh
        ENDSELECT
      ENDIF
    ENDDO
    CALL CPU_TIME(finish)

    thisSA%total_steps=step-1

    !set to the best state we ended up finding.
    SELECTTYPE(thisSA)
      TYPEIS(sa_comb_type)
        thisSA%state_curr=thisSA%state_best
      TYPEIS(sa_cont_type)
        thisSA%state_curr=thisSA%state_best
    ENDSELECT
  ENDSUBROUTINE optimize

  !function for the acceptance probability
  FUNCTION accept_prob(e_current,e_neigh,t_current)
    REAL(8),INTENT(IN) :: e_current,e_neigh,t_current
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

  !get a new neighbor state for a combinatorial problem
  FUNCTION get_neigh_comb(thisSA,s_curr)
    CLASS(sa_comb_type),INTENT(INOUT) :: thisSA
    INTEGER,INTENT(IN) :: s_curr(:)
    INTEGER,DIMENSION(SIZE(s_curr)) :: get_neigh_comb

    INTEGER :: j1,j2
    REAL(8) :: temp_r

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)temp_r=thisSA%e_best

    get_neigh_comb=s_curr
    !get switch indeces
    DO
      CALL random_number(temp_r)
      j1=1+FLOOR(SIZE(s_curr)*temp_r)
      CALL random_number(temp_r)
      j2=1+FLOOR(SIZE(s_curr)*temp_r)
      !make sure they're actually different
      IF(j1 .NE. j2)EXIT
    ENDDO
    !set the new neighbor by swapping those points
    get_neigh_comb(j1)=s_curr(j2)
    get_neigh_comb(j2)=s_curr(j1)
  ENDFUNCTION get_neigh_comb

  !get a new neighbor state for a continuous problem
  FUNCTION get_neigh_cont(thisSA,s_curr,damping,smin,smax)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: s_curr(:)
    REAL(8),INTENT(IN),OPTIONAL :: damping
    REAL(8),INTENT(IN),OPTIONAL :: smin,smax
    REAL(8),DIMENSION(SIZE(s_curr)) :: get_neigh_cont

    REAL(8) :: damp_app
    REAL(8) :: temp_r,max_ch,min_ch
    INTEGER :: i

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)temp_r=thisSA%e_best

    !set the damping factor and bounds
    IF(PRESENT(damping))THEN
      damp_app=damping
    ELSE
      damp_app=1.0
    ENDIF
    IF(PRESENT(smin))THEN
      min_ch=smin
    ELSE
      min_ch=MINVAL(s_curr)
    ENDIF
    IF(PRESENT(smax))THEN
      max_ch=smax
    ELSE
      max_ch=MAXVAL(s_curr)
    ENDIF

    DO i=1,SIZE(s_curr)
      CALL random_number(temp_r)
      !perturb the state
      get_neigh_cont(i)=s_curr(i)+(temp_r-0.5)*damp_app
      !make sure it doesn't go out of bounds
      get_neigh_cont(i)=MAX(MIN(get_neigh_cont(i),max_ch),min_ch)
    ENDDO
  ENDFUNCTION get_neigh_cont

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
      CASE('LinAdd')
        thisSA%cool => lin_add_cool
      CASE('QuadAdd')
        thisSA%cool => quad_add_cool
      CASE('ExpAdd')
        thisSA%cool => exp_add_cool
      CASE('TrigAdd')
        thisSA%cool => trig_add_cool
      CASE('custom')
        !do nothing, it is assumed the user already assgined a custom cooling schedule
        WRITE(*,*)'Using user specified cooling function'
      CASE DEFAULT
        WRITE(*,'(2A)')TRIM(ADJUSTL(thisSA%cool_opt)),' is not a valid cooling option'
        STOP 'The user MUST select a valid cooling option or give a custom cooling function'
    ENDSELECT
  ENDSUBROUTINE set_cooling

  !linear multiplicative cooling
  FUNCTION lin_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: lin_mult_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)lin_mult_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    lin_mult_cool=tmax-alpha*k
  ENDFUNCTION lin_mult_cool

  !natural log exponential multiplicative cooling
  FUNCTION exp_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: exp_mult_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)exp_mult_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    exp_mult_cool=tmax*alpha**k
  ENDFUNCTION exp_mult_cool

  !logarithmic multiplicative cooling
  FUNCTION log_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: log_mult_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)log_mult_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    log_mult_cool=tmax/(1.0+alpha*LOG10(k+1.0))
  ENDFUNCTION log_mult_cool

  !quadratic multiplicative cooling
  FUNCTION quad_mult_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: quad_mult_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)quad_mult_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    quad_mult_cool=tmax/(1.0+alpha*k**2)
  ENDFUNCTION quad_mult_cool

  !linear additive cooling
  FUNCTION lin_add_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: lin_add_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)lin_add_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    lin_add_cool=tmin+(tmax-tmin)*(n*1.0-k)/(n*1.0)
  ENDFUNCTION lin_add_cool

  !quadratic additive cooling
  FUNCTION quad_add_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: quad_add_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)quad_add_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    quad_add_cool=tmin+(tmax-tmin)*((n*1.0-k)/(n*1.0))**2
  ENDFUNCTION quad_add_cool

  !exponential additive cooling
  FUNCTION exp_add_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: exp_add_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)exp_add_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    exp_add_cool=tmin+(tmax-tmin)/(1.0+EXP(2.0*LOG(tmax-tmin)/(n*1.0))*(k-0.5*n))
  ENDFUNCTION exp_add_cool

  !trigonometric additive cooling
  FUNCTION trig_add_cool(thisSA,tmin,tmax,alpha,k,n)
    CLASS(sa_type_base),INTENT(INOUT) :: thisSA
    REAL(8),INTENT(IN) :: tmin,tmax,alpha
    INTEGER,INTENT(IN) :: k,n
    REAL(8) :: trig_add_cool

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)trig_add_cool=tmin+tmax+alpha+k+n+thisSA%e_best

    trig_add_cool=tmin+0.5*(tmax-tmin)*(1.0+COS(k*pi/n))
  ENDFUNCTION trig_add_cool
END MODULE sim_anneal
