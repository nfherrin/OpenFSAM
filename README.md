# openSimAnn
A Fortran based open source simulated annealing utility.

This utility consists of a single module that can be generally assigned to solve a [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) optimization problem.
A user can easily add this module to any existing modern Fortran program since the module is self contained and sufficiently abstracted.

To use the simulated annealing module, point to it properly in your make file and add `USE sim_anneal` to the module/program that is using it.
The simulated annealing module has two public types `sa_comb_type` and `sa_cont_type`.
`sa_comb_type` is a combinatorial type simulated annealing optimizer and `sa_cont_type` is a continuous function simulated annealing optimizer.
To use, create a simulated annealing object with `TYPE(sa_comb_type) :: <sa_object>`.
This object now needs to be initialized.
To initialize, the user must specify the following variables
  1. Maximum number of steps in the annealing `<sa_object>%max_step` (`INT` default: `100`)
  2. Alpha parameter for cooling rate `<sa_object>%alpha` (`REAL(8)` default: `0.01`)
  3. Maximum (starting) temperature `<sa_object>%t_max` (`REAL(8)` default: `100`)
  4. Minimum (final) temperature `<sa_object>%t_min` (`REAL(8)` default: `0`)
  5. Cooling schedule option `<sa_object>%cool_opt` (`CHARACTER` default:`LinMult`)
  6. Cooling monotonic option `<sa_object>%mon_cool` (`LOGICAL` default:`.TRUE.`)
  7. The initial guess state variable `<sa_object>%state_curr` must be initialized and set.
    For `sa_comb_type`, this is a one dimensional integer array pointer.
    For `sa_cont_type`, this is a one dimensional double (`REAL(8)`) array pointer.
    This array has no default and MUST be both allocated and set by the user for the annealing to work.
    As a pointer, allocation can either occur through traditional allocation, or through pointing it to an already existing state variable array in the user's code.
  8. The energy function `<sa_object>%energy` must be assigned.
    This is a pointer function that points to a user defined energy function which takes in the state array and outputs an energy.
    This function has no default and MUST be assigned by the user for the annealing to work.
  9. \[For continuous annealing problems ONLY\] The damping factor `<sa_object>%damping` may be set (`REAL(8)` default: `1.0`).
  10. \[For continuous annealing problems ONLY\] The minimum state variable value, `<sa_object>%smin`, and the maximum state variable value, `<sa_object>%smax`, may be set (`REAL(8)` default: `0.0`).
    If they are not set or if the maximum is set below the minimum, the minimum and maximum values are taken from the minimum and maximum values of the initial state variables.

A brief explanation of the user initialized variables:
  1. The maximum number of steps, `n`, is the most steps the annealing counter will reach before stopping.
    If the cooling schedule is monotonic and additive, then this amount of steps will always be fully used before program completion.
    If the cooling schedule is non-monotonic and/or multiplicative, then the annealing may use less that `n` steps before the minimum temperature is reached.
  2. The alpha parameter is a user defined parameter for a multiplicative cooling schedule and the way in which it affects cooling depends on the schedule used.
  3. The maximum temperature is the initial temperature.
    For non-monotonic cooling the temperature may briefly go above this temperature.
  4. The minimum temperature is the temperature at which the annealing is considered complete.
    If the cooling schedule is monotonic and additive, this temperature will always be reached before stopping.
    If the cooling schedule is non-monotonic and/or multiplicative, then the annealing may reach the maximum number of steps and terminate before this temperature is reached.
  5. The cooling schedule is an equation used to determine how much to cool the temperature by each step.
    Available cooling schedules are taken from [A Comparison of Cooling Schedules for Simulated Annealing](http://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/) and are:
      1. `LinMult` - Linear multiplicative cooling
      2. `ExpMult` - Exponential multiplicative cooling
      3. `LogMult` - Logarithmic multiplicative cooling
      4. `QuadMult` - Quadratic multiplicative cooling
      5. `LinAdd` - Linear additive cooling
      6. `QuadAdd` - Quadratic additive cooling
      7. `ExpAdd` - Exponential additive cooling
      8. `TrigAdd` - Trigonometric additive cooling
    If they so desire, the user may also specify a custom cooling schedule beyond these by specifying the type as `custom` and pointing `<sa_object>%cool` to their desired cooling function.
    For specifying a custom cooling schedule, the generic cooling function interface takes in minimum temperature, maximum temperature, alpha, current step, and maximum number of steps.
  6. If the cooling is monotonic then the selected schedule alone is used.
    Otherwise the cooling result is multiplied by mu at each step, as described in the cooling reference.
  7. The initial guess state variable is the set and values of parameters that the annealing routine perturbs and then computes the energy of in order to optimize.
    In the example of the traveling salesman, it is the array of order of customer visitation by the salesman.
    Perturbation for combinatorial problems involves swapping two values.
    Perturbation for continuous problems involves perturbing each value by a random number from negative one half to positive one half times the damping factor.
  8. The energy function is a calculation of a given state variable's energy.
    Simulated annealing searches for the minimum energy, so the user should make certain that their energy function is minimum at the desired optimal result.
    In the example of the traveling salesman, the energy is the total path length the salesman must travel for the given order.
  9. \[For continuous annealing problems ONLY\] The damping factor is a factor determining how much state variables may be perturbed by for each neighbor perturbation.
    Typically this should be small relative to the magnitude of the state variables, but not so small that perturbations create essentially identical states.
    Half the damping factor is the maximum perturbation, positive or negative, that each state variable might be subjected to when finding a new neighbor.
  10. \[For continuous annealing problems ONLY\] The minimum and maximum state variables prevent the continuous state variables from being perturbed outside specified bounds of the problem.
    The user should ensure that only problem breaking values are excluded, otherwise the annealing may not be able to find the optimum state.

The user may now use the simulated annealing optimization in their code by calling `<sa_object>%optimize`.
This subroutine results in the optimal state array found stored in `<sa_object>%state_best` and the energy of that state is stored in `<sa_object>%e_best`.

An example of usage of this utility is given in `examples/traveling_sales_general` in which simulated annealing is used to optimize a traveling salesman problem.
The simulated annealing object is defined in `globals.f90` as
```
TYPE(sa_comb_type) :: ts_simanneal
```
And is initialized in `travel_sales.f90` with the following code:
```
  SUBROUTINE setup_ts_sa()
    INTEGER :: i

    ts_simanneal%max_step=1000000
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
```
where the energy function is defined as
```
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
```

With the traveling salesman problem setup and simulated annealing initialization complete, the optimization is then called in `main.f90` as
```
  CALL ts_simanneal%optimize()
```

Which can then view the optimal state array in the form of `ts_simanneal%state_best` with an energy of `ts_simanneal%e_best`.