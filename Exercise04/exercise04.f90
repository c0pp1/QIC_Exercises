!> Harmonic oscillator
!! get the first k_user eigenvalues and eigenvectors from
!! a one dimensional harmonic oscillator
!!  H = p**2 + w**2*q**2
!!
program main
    use lapack_wrapper

    implicit none

    integer :: ii, k, k_user, info, outunit=100
    double precision, allocatable :: x(:), H(:, :), eigenvals(:)
    double precision :: dx, xlims(2) = (/-1._8,1._8/), w = 1._8, precision = 0.01
    logical :: skip_next = .false.
    character(len=32) :: arg, ofdir = './exercise04'

    ! read command line arguments
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


    k = ceiling(real(k_user, kind=8) / (4._8*precision))

    ! allocate memory for discretized space, hamiltonian, and eigenvalues
    allocate(x(k), H(k, k), eigenvals(k))
    dx = abs(xlims(2)-xlims(1))/real(k-1, kind=8)
    !print *, dx
    x = (/(minval(xlims)+ii*dx, ii=0,k-1, 1)/)
    !print *, x
    call get_Hamiltonian(H, w, x, dx)
    !print *, H
    ! JOBZ='V' gives the eigenstates as columns in the input matrix
    call ev_wrap(JOBZ='V', A=H, W=eigenvals, INFO=info)
    if(info /= 0) then
        print *, "ERROR: failed to get eigenvalues"
        stop
    end if

    open(unit=outunit, file=trim(ofdir) // '.txt', action='write')
    do ii=1, k_user
        write(outunit, *) H(ii, :k_user)
    end do
    write(outunit, *) eigenvals(:k_user)
    close(outunit)

    contains

        subroutine print_help()
            print *, "Usage: ./exercise04 -k [#eigen to calculate] -x [xlims, comma separated without blank spaces] -w [frequency]"
        end subroutine print_help

        subroutine get_Hamiltonian(A, wi, xi, dxi)
            double precision, intent(inout) :: A(:, :)
            double precision, intent(in) :: wi, xi(:), dxi
            integer :: m, n, ij, jj

            m = ubound(A, 1)
            n = ubound(A, 2)
            if(m/=n) then
                print *, "ERROR: matrix must be square"
                stop
            end if
            A = 0._8

            
            do jj = 1, n
                do ij=1, m
                    if(ij==jj) A(ij, jj) = 2._8/(dxi*dxi) + potential(xi(ij), wi)
                    if((ij==jj+1) .or. (ij==jj-1)) A(ij, jj) = -1._8/(dxi*dxi)
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

end program
