!traveling salesman initialization and energy functions
MODULE travel_sales
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ts_init
  INTEGER,ALLOCATABLE :: glob_ord(:),min_ord(:),max_ord(:)
  REAL(8) :: min_len=1.0D+308,max_len=0.0,avg_len=0.0
CONTAINS
  !traveling salesman initialization problem
  SUBROUTINE ts_init()
    INTEGER :: i
    INTEGER :: fact=1
    REAL(8) :: start,finish

    ALLOCATE(cust_locs(num_customers,2))
    !each customer is given a random 2D location on a 1x1 grid
    DO i=1,num_customers
      CALL random_number(cust_locs(i,1))
      CALL random_number(cust_locs(i,2))
    ENDDO

    ALLOCATE(glob_ord(num_customers))
    ALLOCATE(min_ord(num_customers))
    ALLOCATE(max_ord(num_customers))

    !compute number of permutations and initialize the first permutation (just ordered)
    DO i=1,num_customers
      fact=fact*i
      glob_ord(i)=i
    ENDDO

    WRITE(*,'(A,I0)')' Number of possible paths: ',fact

    !find minimum path length brute force wise (only if num_customers<=15)
    IF(num_customers .LE. 13)THEN
      CALL CPU_TIME(start)
      CALL find_min(1)
      CALL CPU_TIME(finish)
      WRITE(*,*)'Average path length: ',avg_len/(fact*1.0)
      WRITE(*,*)'Minimum length',min_len
      WRITE(*,*)'Minimum order',min_ord
      WRITE(*,*)'Maximum length',max_len
      WRITE(*,*)'Maximum order',max_ord
      WRITE(*,*)'Brute force optimization finished in: ',finish-start,' seconds'
    ENDIF
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
END MODULE travel_sales
