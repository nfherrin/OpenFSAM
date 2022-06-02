!simmulated annealing functions
MODULE sim_anneal
  USE globals
  USE travel_sales
  IMPLICIT NONE
  PRIVATE
  PUBLIC simulate_anneal

  REAL(8) :: alpha=0.0001,t_current,t_max=100,t_min=0,e_current,e_best
  INTEGER,ALLOCATABLE :: s_current(:),s_neigh(:)
  INTEGER :: step
CONTAINS

  !simulate annealing for the traveling salesman
  SUBROUTINE simulate_anneal()
    REAL(8) :: e_neigh
    INTEGER,ALLOCATABLE :: s_best(:)
    INTEGER :: max_step,i,step_actual
    REAL(8) :: temp_r,start,finish

    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'**********************************************************************************'
    !write(*,*)'performing simulated annealing'
    !write(*,*)'**********************************************************************************'

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
    max_step=100000000
    step=0
    step_actual=0
    CALL CPU_TIME(start)
    DO WHILE(step .LE. max_step .AND. t_current .GE. t_min)
      step=step+1
      step_actual=step_actual+1
      !get a new neighbor and compute energy
      CALL get_neigh()
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
  ENDSUBROUTINE simulate_anneal

  !get a new neighbor state
  SUBROUTINE get_neigh()
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
  ENDSUBROUTINE get_neigh

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
END MODULE sim_anneal
