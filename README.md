# OpenFSAM

---
## Description
---
OpenFSAM (**O**pen source **F**ortran **S**imulated **A**nnealing **M**odule).
A Fortran based open source simulated annealing module.

This module consists of a single module that can be generally assigned to solve a [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) optimization problem.
A user can easily add this module to any existing modern Fortran program since the module is self contained and sufficiently abstracted.

To use the simulated annealing module, point to it properly in your make file and add `USE OpenFSAM` to the module/program that is using it.
The simulated annealing module has three public types `sa_comb_type`, `sa_cont_type`, and `sa_disc_type`.
`sa_comb_type` is a combinatorial type simulated annealing optimizer, `sa_cont_type` is a continuous function simulated annealing optimizer, and `sa_disc_type` is a discrete function simulated annealing optimizer.
To use, create a simulated annealing object with `TYPE(sa_comb_type) :: <sa_object>`, `TYPE(sa_cont_type) :: <sa_object>`, or `TYPE(sa_disc_type)`.
This object now needs to be initialized.
To initialize, the user must specify the following variables
  1. Maximum number of steps in the annealing `<sa_object>%max_step` (`INT` default: `100`)
  2. Alpha parameter for cooling rate `<sa_object>%alpha` (`REAL(8)` default: `0.01`)
  3. Maximum (starting) temperature `<sa_object>%t_max` (`REAL(8)` default: `100`)
  4. Minimum (final) temperature `<sa_object>%t_min` (`REAL(8)` default: `0`)
  5. Cooling schedule option `<sa_object>%cool_opt` (`CHARACTER` default: `LinMult`)
  6. Cooling monotonic option `<sa_object>%mon_cool` (`LOGICAL` default: `.TRUE.`)
  7. Progress bar `<sa_object>%prog_bar` (`LOGICAL` default: `.FALSE.`)
  8. Restart Value `<sa_object>%resvar` (`REAL(8)` default: `0`)
  9. The initial guess state variable `<sa_object>%state_curr` must be initialized and set.
    For both `sa_comb_type` and `sa_disc_type`, this is a one dimensional integer array pointer.
    For `sa_cont_type`, this is a one dimensional double (`REAL(8)`) array pointer.
    This array has no default and MUST be both allocated and set by the user for the annealing to work.
    As a pointer, allocation can either occur through traditional allocation, or through pointing it to an already existing state variable array in the user's code.
  10. The energy function `<sa_object>%energy` must be assigned.
    This is a pointer function that points to a user defined energy function which takes in the state array and outputs an energy.
    This function has no default and MUST be assigned by the user for the annealing to work.
  11. \[For continuous annealing problems ONLY\] The damping factor `<sa_object>%damping` may be set (`REAL(8)` default: `1.0`).
  12. \[For continuous annealing problems ONLY\] The minimum state variable value, `<sa_object>%smin`, and the maximum state variable value, `<sa_object>%smax`, may be set (`REAL(8)` default: `0.0`).
    If they are not set or if the maximum is set below the minimum, the minimum and maximum values are taken from the minimum and maximum values of the initial state variables.
  13. \[For continuous annealing problems ONLY\] Dynamic damping option `<sa_object>%damp_dyn` (`LOGICAL` default: `.FALSE.`).
  14. \[For continuous annealing problems and discrete annealing problems ONLY\] Number of parameters to perturb for generating neighbors `<sa_object>%num_perturb` (`INTEGER` default: 0).
  15. \[For discrete annealing problems ONLY\] Acceptable integer values that each of the parameters to the energy function can adopt are stored in `<sa_object>%var_values`. (`INTEGER` default: array of all unique values in the initial state)

A brief explanation of the user initialized variables:
  1. The maximum number of steps, `n`, is the most steps the annealing counter will reach before stopping.
    If the cooling schedule is monotonic and additive, then this amount of steps will always be fully used before program completion.
    If the cooling schedule is non-monotonic and/or multiplicative, then the annealing may use less that `n` steps before the minimum temperature is reached.
    For this implementation of simulated annealing, the number of steps always increase if a neighbor state is accepted, otherwise if the neighbor state is rejected, the step count has a 1% chance of increasing.
    This feature allows for an even longer tail on simulated annealing cooling which been shown to be effective at increasing optimization.
  2. The alpha parameter is a user defined parameter for a multiplicative cooling schedule and the way in which it affects cooling depends on the schedule used.
  3. The maximum temperature is the initial temperature.
    For non-monotonic cooling the temperature may briefly go above this temperature.
  4. The minimum temperature is the temperature at which the annealing is considered complete.
    If the cooling schedule is monotonic and additive, this temperature will always be reached before stopping.
    If the cooling schedule is non-monotonic and/or multiplicative, then the annealing may reach the maximum number of steps and terminate before this temperature is reached.
  5. The cooling schedule is an equation used to determine how much to cool the temperature by each step.
    If they so desire, the user may also specify a custom cooling schedule beyond these by specifying the type as `custom` and pointing `<sa_object>%cool` to their desired cooling function.
    For specifying a custom cooling schedule, the generic cooling function interface takes in minimum temperature, maximum temperature, alpha, current step, and maximum number of steps.
    Available pre-made cooling schedules are taken from [*"A Comparison of Cooling Schedules for Simulated Annealing"*](http://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/) and are:
      1. `LinMult` - Linear multiplicative cooling
      2. `ExpMult` - Exponential multiplicative cooling
      3. `LogMult` - Logarithmic multiplicative cooling
      4. `QuadMult` - Quadratic multiplicative cooling
      5. `LinAdd` - Linear additive cooling
      6. `QuadAdd` - Quadratic additive cooling
      7. `ExpAdd` - Exponential additive cooling
      8. `TrigAdd` - Trigonometric additive cooling
  6. If the cooling is monotonic then the selected schedule alone is used.
    Otherwise the cooling result is multiplied by mu at each step, as described in the cooling reference.
  7. If the user specifies a progress bar, then as the annealing progresses the progress bar will be shown as it cools.
    Works best in additive cooling.
  8. If the user specifies a positive restart value, then upon reaching that temperature, the problem will set the current state as the best state found thus far.
    Additionally, the restart value will be halved so this will happen again later on.
    If used, it should be made small so that the annealing has had the opportunity to traverse much of the domain.
    This can be particularly useful in continuous annealing problems with multiple local minimums.
  9. The initial guess state variable is the set and values of parameters that the annealing routine perturbs and then computes the energy of in order to optimize.
    In the example of the traveling salesman, it is the array of order of customer visitation by the salesman.
    Perturbation for combinatorial problems involves swapping two values.
    Perturbation for continuous problems involves perturbing each value by a random number from negative one half to positive one half times the damping factor.
    Perturbation for discrete problems involves changing the value of one or more of the parameters.
  10. The energy function is a calculation of a given state variable's energy.
    Simulated annealing searches for the minimum energy, so the user should make certain that their energy function is minimum at the desired optimal result.
    In the example of the traveling salesman, the energy is the total path length the salesman must travel for the given order.
  11. \[For continuous annealing problems ONLY\] The damping factor is a factor determining how much state variables may be perturbed by for each neighbor perturbation.
    Typically this should be small relative to the magnitude of the state variables, but not so small that perturbations create essentially identical states.
    Half the damping factor is the maximum perturbation, positive or negative, that each state variable might be subjected to when finding a new neighbor.
    The default damping factor is half of the width of the problem bounds, i.e. `(<sa_object>%smax-<sa_object>%smin)/2`.
  12. \[For continuous annealing problems ONLY\] The minimum and maximum state variables prevent the continuous state variables from being perturbed outside specified bounds of the problem.
    The user should ensure that only problem breaking values are excluded, otherwise the annealing may not be able to find the optimum state.
  13. \[For continuous annealing problems ONLY\] Dynamic damping makes it so that each time the restart value is reached, the damping factor is halved.
    As such, dynamic damping ONLY comes into play if a non-zero restart value is given.
    This can be particularly useful for continuous problems where the user may wish to start with a large (often the default) damping factor, to traverse the whole domain.
    Coupled with dynamic damping and a low (but nonzero) restart value, this allows the user to gradually decrease the effective portion of the domain that is being searched towards the end of the annealing.
  14. \[For continuous annealing problems and discrete annealing problems ONLY\] When annealing continuous functions it can sometimes be useful to not perturb every parameter when generating neighbors.
    When annealing discrete functions, it could be desirable to perturb more than one of the input parameters at a time.
    To this end, the user can specify the number of parameters to perturb each iteration during generation of neighbors.
    This can be particularly useful when the energy is heavily dependent on some parameters and less dependent on others.
  15. \[For discrete annealing problems ONLY\] Since different discrete objective functions have different domains, the user is given the responsiblity to explicitly define the domain for their energy function. If they do not, then the domain will be assumed to span the unique values in the given initial state. This default is not a range of integers, but a list. So if the initial state contains \{1,3,6\} then the domain will just be the integers 1, 3, and 6 NOT all integers from 1 to 6.

The user may now use the simulated annealing optimization in their code by calling `<sa_object>%optimize`.
This subroutine results in the optimal state array found stored in `<sa_object>%state_best` and the energy of that state is stored in `<sa_object>%e_best`.

---
## Simulated Annealing for the Traveling Salesman Problem
---

An example of usage of this module for a combinatorial problem is given in `examples/traveling_sales_general/` in which simulated annealing is used to optimize a traveling salesman problem.
The simulated annealing object is defined in `globals.f90` as
```
  TYPE(sa_comb_type) :: ts_simanneal
```
And is initialized in `travel_sales.f90` with the following code:
```
  SUBROUTINE setup_ts_sa()
    INTEGER :: i

    ts_simanneal%max_step=10000*num_customers
    ts_simanneal%t_max=100
    ts_simanneal%t_min=0
    ts_simanneal%cool_opt='QuadAdd'
    ts_simanneal%mon_cool=.FALSE.
    ts_simanneal%prog_bar=.TRUE.
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

Which can then view the optimal state array in the form of `ts_simanneal%state_best` with an energy (path length) of `ts_simanneal%e_best`.

---
## Simulated Annealing for Continuous Function Optimization
---

An example of usage of this module for continuous problems is given in `examples/continuous_test/` in which simulated annealing is used to optimize six functions and a sum of all six functions.
The simulated annealing object is defined in `main.f90` as
```
  TYPE(sa_cont_type) :: func_sa
```
And is initialized in `main.f90` with the following code (shown here just for the first function):
```
  SUBROUTINE setup_sa_func1()
    REAL(8) :: temp_r
    func_sa%max_step=1000
    func_sa%t_max=100.0D0
    func_sa%t_min=1.0D-12
    func_sa%cool_opt='QuadAdd'
    func_sa%mon_cool=.FALSE.
    func_sa%smin=-10.0D0
    func_sa%smax=10.0D0
    func_sa%damping=0.D0
    func_sa%resvar=1.0D0
    func_sa%damp_dyn=.TRUE.
    !give energy function
    func_sa%energy => func1_eg
    !give random initial guess
    ALLOCATE(func_sa%state_curr(1))
    CALL random_number(temp_r)
    func_sa%state_curr=temp_r*20.0d0-10.0d0
  ENDSUBROUTINE setup_sa_func1
```
where the energy function is defined as
```
  FUNCTION func1_eg(thisSA,state_ord)
    CLASS(sa_cont_type),INTENT(INOUT) :: thisSA
    REAL(8),DIMENSION(:),INTENT(IN) :: state_ord
    REAL(8) :: func1_eg

    func1_eg=func1(state_ord(1))
  ENDFUNCTION func1_eg
```
for function:
```
  FUNCTION func1(x)
    REAL(8),INTENT(IN) :: x
    REAL(8) :: func1
    func1=10.0D0*SIN(x)-0.05D0*(x+2.0D0)+(x-1.0D0)**2+20.D0
  ENDFUNCTION func1
```

With the function minimization problem setup and simulated annealing initialization complete, the optimization is then called in `main.f90` as
```
  CALL func_sa%optimize()
```

Which can then view the optimal state array in the form of `func_sa%state_best` with an energy (function value) of `func_sa%e_best`.

This is done for all six functions and a sum of all six functions where each function is dependent upon a different variable.
It can be observed that the sum of all functions takes many more iterations to get a good solution and still the solution is not as good as the individual optimizations.
This is because continuous simulated annealing perturbs every parameter, so a decrease may occur in which one parameter moved in a way that increased the energy, while not increasing it as much as the other perturbations decreased it.
This demonstrates that for continuous simulated annealing, it is preferable to keep the number of parameters as low as possible.
In this case, since the energy is a sum of each function, then the energy does not rely on any cross parameter terms, so that each portion of the energy may be optimized independently.
Indeed, it can be observed in the example that reducing the number of perturbed parameters each iteration is effective in reaching a similar final energy while greatly reducing the number of total iterations used to find said energy.

---
## Simulated Annealing for Simulated Annealing
---

The traveling salesman problem is the problem for which simulated annealing was originally created (or at least, this was the application where it was first named as such).
As such, optimization of the traveling salesman problem is one of the common ways to demonstrate the effectiveness of simulated annealing.
Energy is typically expressed as some monotonically increasing function of the total path length, oftentimes simply the path length itself.
Since it is a combinatorial problem, damping factors and minimums/maximum state values are not a part of the annealing problem.
If we fix the energy function (to just the path length), then for a given additive cooling schedule, the effectiveness will be determined by the maximum temperature, the minimum temperature, and the number of iterations.
Well this is then a continuous optimization problem, figuring out which temperature bounds and how many steps provide the optimal simulated annealing for a given traveling salesman problem size.
Similarly a multiplicative cooling schedule will depend on alpha and the temperature bounds (but not on the number of steps if they are set sufficiently large).

One of the issues here is that determining how good the result of a simulated annealing calculation is can be difficult for large traveling salesman problems since a reference solution must be known.
However, it can be noticed that the only factor in the energy for a traveling salesman problem is total path length.
Since the path length is solely dependent upon the distance between each point regardless of the physical dimensions of the problem, then traveling salesman problems look identical when solved using simulated annealing, regardless of dimensionality.
As such, simulated annealing can be used on a 1D (and only on a 1D) traveling salesman problems to compare to the actual reference optimal solution energy found using a simple and inexpensive sort for an arbitrarily large traveling salesman problem.

The information in the two previous paragraphs allows the creation of a traveling salesman simulated annealing solver optimizer using simulated annealing.
That is to say, for any sized traveling salesman problem, we can search for an optimal number of iterations and temperature bounds using continuous simulated annealing.
This is done in `examples/sa_ts_sa/` where the continuous annealing problem is set up in `sa_ts_sa.f90` as:
```
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
```
where the energy function now solves the traveling salesman problem using simulated annealing with the current state parameters for the number of iterations as well as the temperature bounds.
This is again defined in `sa_ts_sa.f90` as:
```
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
```
The traveling salesman problem is setup each solve similar to the setup shown in the *Simulated Annealing for the Traveling Salesman Problem* example, this time in `travel_sales.f90`:
```
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
```
The actual state variables are converted into the temperature bounds and iterations using a functional transformation allowing them all to be valid in the minimum and maximum range as well as all use the same damping factor, this is defined in `sa_ts_sa.f90`:
```
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
```
This example is fairly tailored towards annealing of a traveling salesman problem of size 10.
For optimizing larger traveling salesman problems, the functional transform should likely be altered to allow for larger maximum iterations and potentially larger maximum temperatures.
Additionally, for very large traveling salesman problems it is likely acceptable for the annealing to not reach the exact solution but rather something sufficiently good.
As such, for very large problems it may be desirable to reduce the weight of the solution error so that states where the exact solution is not always found are acceptable.

---
## Simulated Annealing for the Ising Model
---

Finding the ground state of the [Ising model](https://web.stanford.edu/~jeffjar/statmech/intro4.html) is a discrete optimization problem. The Hamiltonian for such a system

$$H ( \vec{\sigma} ) = - \sum_{<i,j>} J_{ij} \sigma_i \sigma_j - \mu \sum_{j} h_j \sigma_{j}$$

is minimized assigning the binary values of $[-1, 1]$ to each of the $n$ elements in the vector $\vec{\sigma}$.

An example of ising model optimization is contained in `examples/ising/`. Initialization is straightforward:
```
    annealer%max_step = 100
    annealer%alpha = 0.01
    annealer%t_max = 100.0
    annealer%t_min = 0.0
    annealer%cool_opt = "QuadAdd"
    annealer%mon_cool = .true.
    annealer%prog_bar = .true.
    annealer%resvar = 0.0
    annealer%energy => ising_hamiltonian
    annealer%var_values = [1, -1]
    annealer%num_perturb = 10
    allocate(annealer%state_curr(n_spins))
    annealer%state_curr = state
```
The Hamiltonian objective function is given as:
```
    !> Ising Model Hamiltonian (with Magnetic Moment, mu = 1)
    function ising_hamiltonian(sa, state)
        class(sa_disc_type), intent(inout) :: sa
        integer, dimension(:), intent(in) :: state
        integer, dimension(:,:), allocatable :: state_mat
        real(8) :: ising_hamiltonian, r

        state_mat = reshape(spread(state, 2, size(state)), [size(state), size(state)])
        ising_hamiltonian = -1 * sum(j * state_mat * transpose(state_mat)) - dot_product(h, state)

    end function ising_hamiltonian
```
Most significantly, notice that the annealer was initialized with an array of discrete values which are considered valid, `annealer%var_values = [1, -1]`. You are free to choose any number of discrete integers to populate this parameter. For example, the [**Q**uadratic **U**nconstrained **B**inary **O**ptimization (QUBO) problem](https://arxiv.org/abs/1705.09844#) is very similar in structure to the Ising model, but it uses [0, 1] as valid parameter values.
