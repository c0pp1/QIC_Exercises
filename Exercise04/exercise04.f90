!> Harmonic oscillator
!! get the first k_user eigenvalues and eigenvectors from
!! a one dimensional harmonic oscillator
!!  H = p**2 + w**2*q**2
!!
!! with "--precision" you can set the relative precision of the output
!! eigenvalues. The relation between precision and the number of points
!! used in discretization k has been deduced through graphical results
!! studying the case k = k_user.
!! N.B. this relation is valid only if precision < 3%
!!
program main
    use lapack_wrapper

    implicit none

    integer :: ii, k, k_user, info, outunit=100
    double precision, allocatable :: x(:), H(:, :), eigenvals(:)
    double precision :: dx, xlims(2) = (/-1.,1./), w = 1., precision = 0.01
    logical :: skip_next = .false.
    character(len=32) :: arg, ofdir = './exercise04'

    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-k')
                call get_command_argument(ii+1, arg)
                read (arg, *) k_user
                skip_next = .true.

            case ('-w', '--frequency')
                call get_command_argument(ii+1, arg)
                read (arg, *) w
                skip_next = .true.

            case ('-p', '--precision')
                call get_command_argument(ii+1, arg)
                read (arg, *) precision
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

    ! this relation guarantees the precision the user wants
    k = ceiling(real(k_user, kind=8) / (10._8*precision))

    ! allocate memory for discretized space, hamiltonian, and eigenvalues
    allocate(x(k), H(k, k), eigenvals(k))
    dx = abs(xlims(2)-xlims(1))/real(k-1)
    x = (/(minval(xlims)+ii*dx, ii=0,k-1, 1)/)

    call get_Hamiltonian(H, w, x, dx)

    ! JOBZ='V' gives the eigenstates as columns in the input matrix
    call ev_wrap(JOBZ='V', A=H, W=eigenvals, INFO=info)
    if(info /= 0) then
        print *, "ERROR: failed to get eigenvalues"
        stop
    end if

    open(unit=outunit, file=trim(ofdir) // '.txt', action='write')
    do ii=1, k
        write(outunit, *) H(ii, :)
    end do
    write(outunit, *) eigenvals
    close(outunit)

    contains

        subroutine print_help()
            print *, "Usage: ./exercise04 -k [#eigen to calculate] -x [xlims, comma separated without blank spaces] -w [frequency]"
        end subroutine print_help

        subroutine get_Hamiltonian(A, w, x, dx)
            double precision, intent(inout) :: A(:, :)
            double precision, intent(in) :: w, x(:), dx
            integer :: m, n, ii, jj

            m = ubound(A, 1)
            n = ubound(A, 2)
            if(m/=n) then
                print *, "ERROR: matrix must be square"
                stop
            end if
            A = 0._8

            do ii=1, m
                do jj = 1, n
                    if(ii==jj) A(ii, jj) = 1./dx**2 + potential(x(ii), w)/2
                    if((ii==jj+1) .or. (ii==jj-1)) A(ii, jj) = -0.5/dx**2
                end do
            end do

        end subroutine get_Hamiltonian

        function potential(x, w_user) result(V)
            double precision :: x
            double precision, intent(in), optional  :: w_user
            double precision :: V
            double precision :: w = 1._8

            if(present(w_user)) w = w_user
            
            V = w**2*x**2

        end function potential

end program
