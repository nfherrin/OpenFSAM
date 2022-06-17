!traveling salesman initialization and energy functions
MODULE travel_sales
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC ts_init,path_len,dist,setup_ts_sa
CONTAINS
  !traveling salesman initialization problem
  SUBROUTINE ts_init()
    INTEGER :: i,j
    INTEGER,ALLOCATABLE :: min_ord(:)
    !each customer is given a random location on a 1 grid
    cust_locs=0.0
    DO i=1,num_customers
      DO j=1,prob_dim
        CALL random_number(cust_locs(i,prob_dim))
      ENDDO
    ENDDO

    ALLOCATE(min_ord(num_customers))

    !initialize the problem array (just ordered)
    DO i=1,num_customers
      min_ord(i)=i
    ENDDO

    !sort to get the exact solution
    CALL Bubble_Sort(min_ord)
    sort_best=path_len(min_ord)

    DEALLOCATE(min_ord)
  ENDSUBROUTINE ts_init

  FUNCTION path_len(state_ord)
    INTEGER,INTENT(IN) :: state_ord(num_customers)
    INTEGER :: i
    REAL(8) :: path_len
    path_len=0
    DO i=1,num_customers-1
      path_len=path_len+dist(cust_locs(state_ord(i),:),cust_locs(state_ord(i+1),:))
    ENDDO
  ENDFUNCTION path_len

  FUNCTION dist(loc1,loc2)
    REAL(8),INTENT(IN) :: loc1(:),loc2(:)
    INTEGER :: dim_prob,i
    REAL(8) :: dist
    dim_prob=SIZE(loc1)
    dist=0
    DO i=1,dim_prob
      dist=dist+(loc1(i)-loc2(i))**2
    ENDDO
    dist=SQRT(dist)
  ENDFUNCTION dist

  SUBROUTINE Bubble_Sort(a)
    INTEGER, INTENT(in out), DIMENSION(:) :: a
    INTEGER :: temp
    INTEGER :: i, j
    LOGICAL :: swapped

    DO j = SIZE(a)-1, 1, -1
      swapped = .FALSE.
      DO i = 1, j
        IF (cust_locs(a(i),1) > cust_locs(a(i+1),1)) THEN
          temp = a(i)
          a(i) = a(i+1)
          a(i+1) = temp
          swapped = .TRUE.
        END IF
      END DO
      IF (.NOT. swapped) EXIT
    END DO
  END SUBROUTINE Bubble_Sort

  SUBROUTINE setup_ts_sa(invars)
    REAL(8),INTENT(IN) :: invars(3)
    INTEGER :: i

    ts_simanneal%max_step=NINT(invars(1))
    ts_simanneal%t_max=invars(2)
    ts_simanneal%t_min=invars(3)
    ts_simanneal%cool_opt='QuadAdd'
    ts_simanneal%mon_cool=.FALSE.
    ALLOCATE(ts_simanneal%state_curr(num_customers))
    DO i=1,num_customers
      ts_simanneal%state_curr(i)=i
    ENDDO
    !point to a path length function that works with the SA type
    ts_simanneal%energy => path_len_eg
  ENDSUBROUTINE setup_ts_sa

  !path length function for the energy calculation
  FUNCTION path_len_eg(thisSA,state_ord)
    CLASS(sa_comb_type),INTENT(INOUT) :: thisSA
    INTEGER,DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: path_len_eg

    INTEGER :: i

    !this line is literally just to insure that it doesn't complain about not using the variables
    IF(.FALSE.)i=thisSA%size_states

    path_len_eg=0
    DO i=1,SIZE(state_ord)-1
      path_len_eg=path_len_eg+dist(cust_locs(state_ord(i),:),cust_locs(state_ord(i+1),:))
    ENDDO
  ENDFUNCTION path_len_eg
END MODULE travel_sales
