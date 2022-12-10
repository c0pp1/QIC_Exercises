module time_independent
    use lapack_wrapper
    use checkpoint_mod
    implicit none

    double precision, parameter :: PI = 4.D0*DATAN(1.D0)

    contains

        subroutine get_groundstate(H, w, k, eigenvals, xlims, x, info, p)
            double precision, allocatable, intent(inout) :: H(:,:), eigenvals(:), p(:)
            double precision, intent(in) :: xlims(2), w
            intent(out) :: x
            optional :: p
            integer :: ii
            integer, intent(in) :: k
            integer, intent(inout) :: info
            double precision, allocatable :: x(:)
            double precision :: dx

            call checkpoint(3, 'Allocating H0, space, and eigenvalues vectors...')
            allocate(H(k,k), x(k), eigenvals(k), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate H0, space, and eigenvalues vectors!')

            dx = abs(xlims(2)-xlims(1))/real(k-1, kind=8)
            call checkpoint(3, 'spacing value is:', scalar=dx)
            x = (/(minval(xlims)+ii*dx, ii=0,k-1, 1)/)
            
            call get_Hamiltonian(H, w, x, dx)
            ! JOBZ='V' gives the eigenstates as columns in the input matrix
            call ev_wrap(JOBZ='V', A=H, W=eigenvals, INFO=info)
            eigenvals = eigenvals / 2._8

            if(present(p)) then
                if(.not. allocated(p)) allocate(p(k))
                do ii=0, floor(k/2.)
                    p(ii+1) = PI*ii/(x(k)-x(1))
                end do
                do ii=floor(k/2.)+1, k-1
                    p(ii+1) = PI*(ii-k)/(x(k)-x(1))
                end do
            end if

        end subroutine

        subroutine get_Hamiltonian(A, w, x, dx)
            double precision, intent(inout) :: A(:, :)
            double precision, intent(in) :: w, x(:), dx
            integer :: m, n, ii, jj

            m = ubound(A, 1)
            n = ubound(A, 2)
            if(m/=n) then
                call checkpoint(1, "matrix must be square")
                stop
            end if
            A = 0._8

            
            do jj = 1, n
                do ii=1, m
                    if(ii==jj) A(ii, jj) = 2._8/(dx*dx) + potential(x(ii), w)
                    if((ii==jj+1) .or. (ii==jj-1)) A(ii, jj) = -1._8/(dx*dx)
                end do
            end do

        end subroutine get_Hamiltonian

        function potential(xj, w_user) result(V)
            double precision :: xj
            double precision, intent(in), optional  :: w_user
            double precision :: V
            double precision :: wj = 1._8

            if(present(w_user)) wj = w_user
            
            V = wj**2*xj**2

        end function potential

end module time_independent

module trotter_suzuki
    use, intrinsic :: iso_c_binding
    use checkpoint_mod
    implicit none   
    include "fftw3.f" 
    double precision, parameter :: PI = 4.D0*DATAN(1.D0)
    
    contains

    subroutine propagate_T(T, nt, xlims)
        double precision, intent(in) :: T, xlims(2)
        integer, intent(in) :: nt
        double precision :: dt

        dt = T/(nt-1)

    end subroutine
    
    function potential_t(x, step, times, w_user) result(V)
        double precision, intent(in) :: x(:), times(:)
        double precision, intent(in), optional  :: w_user
        double precision :: w = 1._8, V(ubound(x, 1))
        integer, intent(in) :: step

        if(present(w_user)) w = w_user
        
        V = w**2*(x-times(step)/times(ubound(times, 1)))**2
        !call checkpoint(5, "x-t/T at step:", scalar=real(step, 8), vector=x-times(step)/times(ubound(times, 1)))
        !V = w**2*(x-2)**2

    end function

    subroutine propagate_step(plan_f, plan_b, times, dt, i, state, state_fft, x, p)
        intent(in) :: i, times, dt, p, x, state, plan_f, plan_b
        intent(inout) :: state_fft
        double precision :: times(:), dt, p(:), x(:)
        double complex :: state_fft(:), state(:)
        !pointer :: state_fft
        integer :: i, k
        integer*8 :: plan_f, plan_b

        call checkpoint(5, 'applying first V/2')
        state_fft = cdexp(cmplx(0., -1, kind=8)*0.5*dt*potential_t(x, i, times))*state

        call checkpoint(5, 'applying FFT')
        call dfftw_execute_dft(plan_f, state_fft, state_fft)

        call checkpoint(5, 'applying T')
        state_fft = cdexp(cmplx(0., -1, kind=8)*dt*p)*state_fft

        call checkpoint(5, 'applying reverse FFT')
        call dfftw_execute_dft(plan_b, state_fft, state_fft)
        call checkpoint(5, 'applying last V/2')
        state_fft = cdexp(cmplx(0., -1, kind=8)*0.5*dt*potential_t(x, i, times))*state_fft

        ! renormalize
        state_fft = state_fft / k

    end subroutine

    subroutine fftw_plans(k, plan_f, plan_b, state_fft)
        integer :: k
        integer*8 :: plan_f, plan_b
        double complex :: state_fft(:), state(:)
        allocatable :: state
        intent(in)  :: k
        intent(out) :: plan_f, plan_b
        optional :: state_fft

        if(present(state_fft)) then
            call checkpoint(3, 'Getting best fftw3 plan based on eigenvector size:', scalar=real(k, kind=8))
            call dfftw_plan_dft_1d(plan_f, k, state_fft, state_fft, FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_1d(plan_b, k, state_fft, state_fft, FFTW_BACKWARD, FFTW_MEASURE)

        else
            call checkpoint(2, 'For better performance pass the vector which will store the groundstate')
            allocate(state(k))
            call checkpoint(3, 'Getting best fftw3 plan based on eigenvector size:', scalar=real(k, kind=8))
            call dfftw_plan_dft_1d(plan_f, k, state, state, FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_1d(plan_b, k, state, state, FFTW_BACKWARD, FFTW_MEASURE)
        end if

    end subroutine

end module trotter_suzuki

program main
    use checkpoint_mod
    use time_independent
    use trotter_suzuki
    implicit none

    integer :: ii, mx=50, nt=50, info, outunit=100, E=1
    integer*8 :: plan_f, plan_b, recl
    double precision, allocatable :: H0(:, :), eigenvals(:), p(:), x(:), times(:)
    double precision :: xlims(2) = (/-1._8,1._8/), w = 1._8, T = 1., dt
    double complex, allocatable :: state_in(:), state_fft(:)
    logical :: skip_next = .false.
    character(len=64) :: arg, ofdir = './exercise05'

    DB = 5

    ! read command line arguments
    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-m')
                call get_command_argument(ii+1, arg)
                read (arg, *) mx
                skip_next = .true.

            case ('-e')
                call get_command_argument(ii+1, arg)
                read (arg, *) E
                E = E+1
                skip_next = .true.

            case ('-n')
                call get_command_argument(ii+1, arg)
                read (arg, *) nt
                skip_next = .true.

            case ('-t')
                call get_command_argument(ii+1, arg)
                read (arg, *) T
                skip_next = .true.

            case ('-w', '--frequency')
                call get_command_argument(ii+1, arg)
                read (arg, *) w
                skip_next = .true.

            case ('-o', '--odirname')
                call get_command_argument(ii+1, arg)
                read (arg, '(A)') ofdir
                skip_next = .true.

            case ('-x', '--xlims')
                call get_command_argument(ii+1, arg)
                read (arg, *) xlims
                skip_next = .true.

            case ('-h', '--help')
                call print_help()
                skip_next = .false.
                stop

            case default
                print '(2a, /)', 'unrecognised command-line option: ', arg
                call print_help()
                skip_next = .false.
                stop
        end select
    end do


    ! get the ground state of the harmonic oscillator
    call checkpoint(3, 'Getting initial state...')
    call get_groundstate(H0, w, mx, eigenvals, xlims, x, info, p)
    if(info/=0) call checkpoint(1, 'Failed to get groundstate!')
    call checkpoint(3, 'The requested eigenvalue is:', scalar=eigenvals(E))
    call checkpoint(4, 'Momentum discretization:', vector=p)

    allocate(state_fft(mx), state_in(mx), times(nt))

    call fftw_plans(mx, plan_f, plan_b, state_fft)
    call checkpoint(3, "The plans are:", vector=real((/plan_f, plan_b/), kind=8))

    dt = T/(nt)
    times = (/(ii*dt, ii=1,nt, 1)/)
    call checkpoint(3, 'The time step is:', scalar=dt)
    call checkpoint(4, 'The times are:', vector=times)

    state_in = cmplx(H0(:, E), 0., kind=8)
    call checkpoint(4, 'Input to time propagation:', vector=real(state_in), &
                    scalar=real(sum(state_in*conjg(state_in))))

    !> call propagate_step(plan_f, plan_b, dt, 1, state_in, state_fft_p, x, p)
    !! this subroutine doesn't work, gives unpredictable results

    arg = '(A, I'
    write(arg, '(A, I1, A)') trim(arg), floor(log10(real(E)))+1, ', A)'
    write(ofdir, arg) trim(ofdir), E, '.dat'   
    inquire(iolength=recl)abs(H0(:, 1))
    call checkpoint(3, 'Output file: ' // trim(ofdir))    
    open(unit=outunit, file=trim(ofdir), action='write', form='unformatted', access='direct', recl=recl)
    write(outunit, rec=1) abs(H0(:, 1))
        
    ! propagate 
    do ii=1, nt
        call checkpoint(5, 'applying first V/2')
        state_fft = cdexp(cmplx(0., -0.5_8*dt*potential_t(x, ii, times), kind=8))*state_in
        
        call checkpoint(5, 'applying FFT')
        call dfftw_execute_dft(plan_f, state_fft, state_fft)
        call checkpoint(5, 'result', vector=abs(state_fft))
        call checkpoint(5, 'applying T')
        state_fft = cdexp(cmplx(0., -1, kind=8)*dt*p)*state_fft

        call checkpoint(5, 'applying reverse FFT')
        call dfftw_execute_dft(plan_b, state_fft, state_fft)
        call checkpoint(5, 'applying last V/2')
        state_fft = cdexp(cmplx(0., -0.5_8*dt*potential_t(x, ii, times), kind=8))*state_fft

        ! renormalize
        state_fft = state_fft / mx
        state_in = state_fft
        write(outunit, rec=ii+1) abs(state_fft)
    end do


    contains
    subroutine print_help()
        print *, "Usage: ./exercise05 -m [#space steps (with prime factors 2,3,5,7)] -n [#time steps] " // &
                 "-x [xlims, comma separated without blank spaces] -t [time T] -w [frequency]"
    end subroutine print_help
end program