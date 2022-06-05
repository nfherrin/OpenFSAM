!traveling salesman initialization and energy functions
MODULE travel_sales
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC ts_init,path_len,dist,setup_ts_sa

  INTEGER,ALLOCATABLE :: glob_ord(:),min_ord(:),max_ord(:)
  REAL(8) :: min_len=1.0D+308,max_len=0.0,avg_len=0.0
CONTAINS
  !traveling salesman initialization problem
  SUBROUTINE ts_init()
    INTEGER :: i,j
    REAL(8) :: num_perms=1
    REAL(8) :: start,finish,est_time

    ALLOCATE(cust_locs(num_customers,prob_dim))
    !each customer is given a random 2D location on a 1x1 grid
    DO i=1,num_customers
      DO j=1,prob_dim
        CALL random_number(cust_locs(i,prob_dim))
      ENDDO
    ENDDO

    ALLOCATE(glob_ord(num_customers))
    ALLOCATE(min_ord(num_customers))
    ALLOCATE(max_ord(num_customers))

    !compute number of permutations and initialize the first permutation (just ordered)
    DO i=1,num_customers
      num_perms=num_perms*i
      glob_ord(i)=i
    ENDDO

    !WRITE(*,'(A,ES16.8)')'Number of possible paths: ',num_perms

    !find minimum path length brute force wise (only if estimated time is under 100 seconds)
    est_time=num_perms*num_customers*prob_dim*2.0E-09
    WRITE(*,'(A,ES16.8,A)')'Estimated brute force calculation time ',est_time,' seconds'
    sort_best=0
    IF(prob_dim .EQ. 1)THEN
      min_ord=glob_ord
      CALL Bubble_Sort(min_ord)
      min_len=path_len(min_ord)
      WRITE(*,'(A,ES16.8)')'Order Minimum path length: ',min_len
      sort_best=min_len
    ELSEIF(est_time .LE. 1.0E+03)THEN
      CALL CPU_TIME(start)
      CALL find_min(1)
      CALL CPU_TIME(finish)
      WRITE(*,'(A,ES16.8)')'Average path length: ',avg_len/num_perms
      WRITE(*,'(A,ES16.8)')'Minimum path length: ',min_len
      WRITE(*,'(A,ES16.8)')'Maximum path length: ',max_len
      WRITE(*,'(A,10000I6)')'Minimum path: ',min_ord
      !WRITE(*,'(A,13I2)')'Maximum path: ',max_ord
      WRITE(*,'(A,ES16.8,A)')'Brute force optimization finished in: ',finish-start,' seconds'
      sort_best=min_len
    ENDIF

    DEALLOCATE(glob_ord,min_ord,max_ord)
  ENDSUBROUTINE ts_init

  RECURSIVE SUBROUTINE find_min(i)
    INTEGER,INTENT(IN) :: i
    INTEGER :: j, t
    REAL(8) :: cur_len=0
    IF(i == num_customers)THEN
      !compute current permutations path length
      cur_len=path_len(glob_ord)
      avg_len=avg_len+cur_len
      !if it is less than the current minimum length, set the current minimum length to it.
      IF(cur_len .LE. min_len)THEN
        min_len=cur_len
        min_ord=glob_ord
      ENDIF
      IF(cur_len .GE. max_len)THEN
        max_len=cur_len
        max_ord=glob_ord
      ENDIF
    ELSE
      DO j = i, num_customers
        t = glob_ord(i)
        glob_ord(i) = glob_ord(j)
        glob_ord(j) = t
        call find_min(i + 1)
        t = glob_ord(i)
        glob_ord(i) = glob_ord(j)
        glob_ord(j) = t
      ENDDO
    ENDIF
  ENDSUBROUTINE find_min

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

  !sets up the traveling salesman problem
  SUBROUTINE setup_ts_sa()
    INTEGER :: i

    ts_simanneal%size_states=num_customers
    ts_simanneal%max_step=100000000
    ts_simanneal%alpha=0.85
    ts_simanneal%t_max=100
    ts_simanneal%t_min=0
    ts_simanneal%cool_opt='QuadMult'
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

    path_len_eg=0
    DO i=1,SIZE(state_ord)-1
      path_len_eg=path_len_eg+dist(cust_locs(state_ord(i),:),cust_locs(state_ord(i+1),:))
    ENDDO
  ENDFUNCTION path_len_eg
END MODULE travel_sales
