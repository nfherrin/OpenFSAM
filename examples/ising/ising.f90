program isingmodel
    use :: OpenFSAM, only : sa_disc_type
    implicit none

    type(sa_disc_type) :: annealer
    integer :: n_spins
    real(8) :: energy_initial
    real(8), dimension(:), allocatable :: h
    real(8), dimension(:,:), allocatable :: j
    integer, dimension(:), allocatable :: state

    n_spins = 1000
    allocate(j(n_spins,n_spins))
    allocate(h(n_spins))
    allocate(state(n_spins))

    call initialize_ising(n_spins, j, h, state)

    !> Initialize Annealer
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

    energy_initial = ising_hamiltonian(annealer, annealer%state_curr)

    write(*, '(A)')'----------------------------------------------------------------------------------'
    write(*, '(A)')'                    finding ising model ground state'
    write(*, '(A)')'----------------------------------------------------------------------------------'

    !> Find Ising Model's Ground State
    call annealer%optimize()

    write(*, '(A)', advance="no") "Random Initial Energy: "
    write(*, '(F0.2)')  energy_initial
    write(*, '(A)', advance="no") "Ground State Energy: "
    write(*, '(F0.2)')  annealer%e_best

    contains

        subroutine initialize_ising(n_spins, j, h, state)
            integer, intent(in) :: n_spins
            real(8), dimension(n_spins:n_spins), intent(inout) :: j
            real(8), dimension(n_spins), intent(inout) :: h
            integer, dimension(n_spins), intent(inout) :: state
            real(8), dimension(n_spins) :: temp_r

            call random_number(j)
            j = 2 * j - 0.5

            call random_number(h)
            h = 2 * h - 0.5

            call random_number(temp_r)
            state = merge(1, -1, (temp_r - 0.5) >= 0)
        end subroutine initialize_ising

        !> Ising Model Hamiltonian (with Magnetic Moment, mu = 1)
        function ising_hamiltonian(sa, state)
            class(sa_disc_type), intent(inout) :: sa
            integer, dimension(:), intent(in) :: state
            integer, dimension(:,:), allocatable :: state_mat
            real(8) :: ising_hamiltonian, r

            !this line is literally just to ensure that it doesn't complain about not using the variables
            if(.FALSE.)r=sa%e_best

            state_mat = reshape(spread(state, 2, size(state)), [size(state), size(state)])
            ising_hamiltonian = -1 * sum(j * state_mat * transpose(state_mat)) - dot_product(h, state)

        end function ising_hamiltonian

end program isingmodel
