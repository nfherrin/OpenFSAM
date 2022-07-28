!Solution for 6 different non-monotonic functions and their combination using SA.
PROGRAM continuous_test
  USE OpenFSAM
  IMPLICIT NONE
  REAL(8) :: minvals(6),minlocs(6),t_max,t_min,smin,smax,damping,resvar,errcomb
  CHARACTER(64) :: cool_opt
  LOGICAL :: mon_cool,damp_dyn
  INTEGER :: ios,max_step

  TYPE(sa_cont_type) :: func_sa

  !minimum values and locations
  minvals(1)=1.5449192781044832D+01
  minvals(2)=1.4184540054199392D+01
  minvals(3)=1.5419753322096613D+01
  minvals(4)=1.1714693326018384D+01
  minvals(5)=7.7653742475768439D-01
  minvals(6)=1.4790147335588884D+01

  minlocs(1)=-1.1260671421437776D+00
  minlocs(2)= 2.8556531452530787D+00
  minlocs(3)= 8.5763080142082426D+00
  minlocs(4)= 8.9206430998425876D+00
  minlocs(5)=-8.3519922528629209D+00
  minlocs(6)=-2.9760703957656585D+00

  !SA parameters
  max_step=1000
  t_max=100.0D0
  t_min=1.0D-12
  cool_opt='QuadAdd'
  mon_cool=.FALSE.
  smin=-10.0D0
  smax=10.0D0
  damping=0.D0
  resvar=1.0D0
  damp_dyn=.TRUE.

  !function 1 solution
  WRITE(*,'(A)')'-----------------------------------function 1-----------------------------------'
  CALL setup_sa_func1()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(1),minlocs(1)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(1))/ABS(minvals(1)) &
    ,ABS(func_sa%state_best(1)-minlocs(1))/ABS(minlocs(1)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  !function 2 solution
  WRITE(*,'(A)')'-----------------------------------function 2-----------------------------------'
  CALL setup_sa_func2()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(2),minlocs(2)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(2))/ABS(minvals(2)) &
    ,ABS(func_sa%state_best(1)-minlocs(2))/ABS(minlocs(2)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  !function 3 solution
  WRITE(*,'(A)')'-----------------------------------function 3-----------------------------------'
  CALL setup_sa_func3()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(3),minlocs(3)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(3))/ABS(minvals(3)) &
    ,ABS(func_sa%state_best(1)-minlocs(3))/ABS(minlocs(3)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  !function 4 solution
  WRITE(*,'(A)')'-----------------------------------function 4-----------------------------------'
  CALL setup_sa_func4()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(4),minlocs(4)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(4))/ABS(minvals(4)) &
    ,ABS(func_sa%state_best(1)-minlocs(4))/ABS(minlocs(4)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  !function 5 solution
  WRITE(*,'(A)')'-----------------------------------function 5-----------------------------------'
  CALL setup_sa_func5()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(5),minlocs(5)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(5))/ABS(minvals(5)) &
    ,ABS(func_sa%state_best(1)-minlocs(5))/ABS(minlocs(5)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  !function 6 solution
  WRITE(*,'(A)')'-----------------------------------function 6-----------------------------------'
  CALL setup_sa_func6()
  CALL func_sa%optimize()
  WRITE(*,'(A)')'          min_erg     min_loc'
  WRITE(*,'(A,40ES12.4)')'ref: ',minvals(6),minlocs(6)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-minvals(6))/ABS(minvals(6)) &
    ,ABS(func_sa%state_best(1)-minlocs(6))/ABS(minlocs(6)),func_sa%total_steps
  WRITE(*,'(A)')
  DEALLOCATE(func_sa%state_curr)

  max_step=max_step*600
  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 6----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 6 variables
  func_sa%num_perturb=6
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 5----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 5 variables
  func_sa%num_perturb=5
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 4----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 4 variables
  func_sa%num_perturb=4
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 3----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 3 variables
  func_sa%num_perturb=3
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 2----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 2 variables
  func_sa%num_perturb=2
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

  !combined functions solution
  WRITE(*,'(A)')'-----------------------------combined funcs perturb 1----------------------------'
  WRITE(*,'(A)')'Annealing: '
  CALL setup_sa_func_comb()
  !perturb 1 variables
  func_sa%num_perturb=1
  CALL func_sa%optimize()
  errcomb=0.0D0
  DO ios=1,6
    errcomb=errcomb+((func_sa%state_best(ios)-minlocs(ios))/minlocs(ios))**2
  ENDDO
  errcomb=SQRT(errcomb)
  WRITE(*,'(A)')'          min_erg   min_loc_1   min_loc_2   min_loc_3   min_loc_4   min_loc_5   min_loc_6'
  WRITE(*,'(A,40ES12.4)')'ref: ',SUM(minvals(:)),minlocs(:)
  WRITE(*,'(A,40ES12.4)')'SA:  ',func_sa%e_best,func_sa%state_best
  WRITE(*,'(A)')' min_erg_err min_loc_err  num_steps'
  WRITE(*,'(2ES12.4,I11)')ABS(func_sa%e_best-SUM(minvals(:)))/SUM(ABS(minvals(:))) &
    ,errcomb,func_sa%total_steps
  DEALLOCATE(func_sa%state_curr)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func1()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func1_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func2()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func2_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func3()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func3_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func4()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func4_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func5()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func5_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func6()
    REAL(8) :: temp_r
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => func6_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_sa_func_comb()
    REAL(8) :: temp_r
    INTEGER :: i
    func_sa%max_step=max_step
    func_sa%t_max=t_max
    func_sa%t_min=t_min
    func_sa%cool_opt=cool_opt
    func_sa%mon_cool=mon_cool
    func_sa%prog_bar=.TRUE.
    func_sa%smin=smin
    func_sa%smax=smax
    func_sa%damping=damping
    func_sa%resvar=resvar
    func_sa%damp_dyn=damp_dyn
    !give energy function
    func_sa%energy => comb_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(6))
    DO i=1,6
      CALL random_number(temp_r)
      func_sa%state_curr(i)=temp_r*20.0d0-10.0d0
    ENDDO
  ENDSUBROUTINE setup_sa_func_comb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func1(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func1
    func1=10.0D0*SIN(x)-0.05D0*(x+2.0D0)+(x-1.0D0)**2+20.D0
  ENDFUNCTION func1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func2(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func2
    func2=1.78D-6*x**8+1.86D-5*x**7-3.75D-4*x**6-3.61D-3*x**5 &
      +2.55D-2*x**4+2.06D-1*x**3-4.85D-1*x**2-3.11D0*x+1.38D0+20.D0
  ENDFUNCTION func2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func3(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func3
    func3=1.11D-6*x**8+7.7D-6*x**7-1.32D-4*x**6-1.35D-3*x**5 &
      +3.59D-3*x**4+6.48D-2*x**3-6.49D-2*x**2-7.32D-1*x+2.57D0+20.D0
  ENDFUNCTION func3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func4(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func4
    func4=3.52D-6*x**8+1.2D-5*x**7-6.26D-4*x**6-1.99D-3*x**5 &
      +3.31D-2*x**4+9.51D-2*x**3-4.89D-1*x**2-1.7D0*x+1.45D0+20.D0
  ENDFUNCTION func4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func5(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func5
    func5=4.0D-6*x**8-1.0D-5*x**7-6.0D-4*x**6+2.0D-3*x**5+3.0D-2*x**4- &
      1.0D-1*x**3-5.0D-1*x**2+2.0D0*x-10.0D0+20.D0
  ENDFUNCTION func5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func6(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func6
    func6=1.57D-6*x**8+3.89D-6*x**7-2.8D-4*x**6-4.28D-4*x**5+ &
      1.39D-2*x**4-4.68D-3*x**3-7.05D-2*x**2+9.53D-1*x-2.87D0+20.D0
  ENDFUNCTION func6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func1_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func1_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func1_eg=func1(state_ord(1))
  ENDFUNCTION func1_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func2_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func2_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func2_eg=func2(state_ord(1))
  ENDFUNCTION func2_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func3_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func3_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func3_eg=func3(state_ord(1))
  ENDFUNCTION func3_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func4_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func4_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func4_eg=func4(state_ord(1))
  ENDFUNCTION func4_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func5_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func5_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func5_eg=func5(state_ord(1))
  ENDFUNCTION func5_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION func6_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func6_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    func6_eg=func6(state_ord(1))
  ENDFUNCTION func6_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION comb_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: comb_eg

    !this is to avoid the make catch complaints
    IF(.FALSE.)write(*,*)thisSA%cool_opt

    comb_eg=0.0D0
    comb_eg=comb_eg+func1(state_ord(1))
    comb_eg=comb_eg+func2(state_ord(2))
    comb_eg=comb_eg+func3(state_ord(3))
    comb_eg=comb_eg+func4(state_ord(4))
    comb_eg=comb_eg+func5(state_ord(5))
    comb_eg=comb_eg+func6(state_ord(6))
  ENDFUNCTION comb_eg
END PROGRAM continuous_test
