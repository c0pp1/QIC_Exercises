program main
    use checkpoint_mod
    implicit none

    double precision, parameter :: PI = 4.D0*DATAN(1.D0)
    character(*), parameter :: ofdir = "./exercise06_"
    integer :: ii, jj, outunit = 100, d = 2, N = 10
    double complex, allocatable :: state(:, :)

    DB = 4

    ! test separable state normalization
    call get_separable_state(d, N, state)

    open(unit=outunit, file=ofdir // "prova.txt", action='write')
    do ii=1,N
        do jj=1,d
            write(outunit, '(ES12.4, ES12.4)', advance='no') real(state(jj, ii)), aimag(state(jj, ii))
        end do
        write(outunit, *)
    end do
    close(outunit)

    contains

        !> subroutine to get a state describing N-body, non interacting, separable
        !! pure state. Randomly initialized
        !!
        !! the state is represented as a dxN matrix, each column representing a different
        !! Hilbert space
        !!
        !! @param[in]   d   the local dimension
        !! @param[in]   N   #bodies in the system
        !! @param[out]  state   the randomly initialized state
        !!
        subroutine get_separable_state(d, N, state)
            integer :: d, N, info
            double complex, allocatable :: state(:, :)
            double precision, allocatable :: rel_phases(:, :), amps(:, :)
            intent(in) :: d, N
            intent(out) :: state

            call checkpoint(3, 'Allocating arrays for separable state computation...')
            if(.not. allocated(state)) allocate(state(d, N), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for separable state computation')
            allocate(rel_phases(d-1, N), amps(d-1, N), stat=info)

            ! generate random parameters: 
            ! angle between 0-2pi
            call random_number(rel_phases)
            rel_phases = rel_phases*2*PI
            ! amplitude between 0-1
            call random_number(amps)

            state(2:, :) = cdexp(cmplx(0., -rel_phases, kind=8))*amps

            ! calculate first coefficients given renormalization constraint
            state(1, :) = sqrt(1._8 - sum(state(2:, :)*conjg(state(2:, :)), dim=1))

        end subroutine

        !> subroutine to get a state describing N-body, interacting, non-separable
        !! pure state. Randomly initialized
        !!
        !! the local dimension d is fixed to 2, i.e. a qubit.
        !! the state is represented as a NxN matrix, each dimension representing 
        !! a different Hilbert space.
        !!
        !! @param[in]   N   #bodies in the system
        !! @param[out]  state   the randomly initialized state
        !!
        subroutine get_general_pure_state(N, state)
            integer :: d=2, N, info
            double complex, allocatable :: state(:, :)
            double precision, allocatable :: rel_phases(:, :), amps(:, :)
            intent(in) :: N
            intent(out) :: state

            call checkpoint(3, 'Allocating arrays for general pure state computation...')
            if(.not. allocated(state)) allocate(state(N, N), rel_phases(d-1, N), amps(d-1, N), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for general pure state computation')


        end subroutine
end program