!> Harmonic oscillator
!! get the first k eigenvalues and eigenvectors from
!! a one dimensional harmonic oscillator
!!  H = p**2 + w**2*q**2
!!
program main
    use lapack_wrapper

    implicit none

    integer :: ii, k, k_user, info, outunit=100
    double precision, allocatable :: x(:), H(:, :), eigenvals(:)
    double precision :: dx, xlims(2) = (/-1.,1./), w = 1., hbar = 6.62d-34
    logical :: skip_next = .false.
    character(len=32) :: arg

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

    !! multiply the number of eigenvalues we want by ten in order
    !! to have enough evaliuation points
    k = k_user * 10

    allocate(x(k), H(k, k), eigenvals(k))
    dx = abs(xlims(2)-xlims(1))/real(k-1)
    x = (/(minval(xlims)+ii*dx, ii=0,k-1, 1)/)

    call get_Hamiltonian(H, w, x)

    call ev_wrap(JOBZ='V', A=H, W=eigenvals, INFO=info)
    if(info /= 0) then
        print *, "ERROR: failed to get eigenvalues"
        stop
    end if

    print *, eigenvals(:k_user)

    open(unit=outunit, file='./prova.txt', action='write')
    do ii=1, k
        write(outunit, *) H(ii, :)
    end do
    close(outunit)

    contains

        subroutine print_help()
            print *, "Usage: ./exercise04 -k [#eigen to calculate] -x [xlims, comma separated without blank spaces] -w [frequency]"
        end subroutine

        subroutine get_Hamiltonian(A, w, x)
            double precision, intent(inout) :: A(:, :)
            double precision, intent(in) :: w, x(:)
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
                    if(ii==jj) A(ii, jj) = 2. + w**2*x(ii)**2
                    if((ii==jj+1) .or. (ii==jj-1)) A(ii, jj) = -1.
                end do
            end do

            A = A*real(m-1)/(x(m) - x(1))
        end subroutine

end program
