program main
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    implicit none

    integer          :: sigma_z(2, 2) = reshape((/1, 0, 0, -1/), (/2,2/)), &
                        sigma_x(2, 2) = reshape((/0, 1, 1, 0/), (/2,2/)), &
                        N = 2, ii, info, outunit = 100, k
    double complex   :: H(:, :)
    double precision :: lambda = 1.d0, eigenvals(:)
    allocatable      :: H, eigenvals
    character(128)   :: arg, ofname, string
    character(1)     :: creturn = char(13)
    logical          :: skip_next = .false., exist

    DB = 4

    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-N')
                call get_command_argument(ii+1, arg)
                read (arg, *) N
                skip_next = .true.

            case ('-dl')
                call get_command_argument(ii+1, arg)
                read (arg, *) DB
                skip_next = .true.

            case ('-l', '--lambda')
                call get_command_argument(ii+1, arg)
                read (arg, *) lambda
                skip_next = .true.

            case ('-k')
                call get_command_argument(ii+1, arg)
                read (arg, *) k
                skip_next = .true.

            case ('-o', '--odirname')
                call get_command_argument(ii+1, arg)
                read (arg, '(A)') ofname
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


    call get_Hamiltonian(N, lambda, H)
    call checkpoint(4, 'Hamiltonian:', cmatrix=H)

    call checkpoint(3, 'Allocating space for eigenvals')
    allocate(eigenvals(2**N), stat=info)
    if(info/=0) then
        call checkpoint(0, 'Failed to allocate space for eigenvals')
        stop 1
    else
        call checkpoint(3, 'Done!')
    end if

    call ev_wrap(A=H, W=eigenvals, INFO=info)
    if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')

    ! save results to file, appending N, lambda, and eigenvals
    call checkpoint(3, "Writing to file '" // trim(ofname) // "'")
    inquire(file=trim(ofname), exist=exist)
    if (exist) then
        open(outunit, file=trim(ofname), status="old", position="append", action="write")
    else
        open(outunit, file=trim(ofname), status="new", action="write")
    end if

    write(outunit, *) N, lambda, eigenvals(:merge(k, 2**N, k<2**N))

    call checkpoint(3, 'All done!')
    contains

        subroutine get_Hamiltonian(N, lambda, H)
            integer :: N, ii, info
            double complex, allocatable, intent(out) :: H(:, :)
            double precision :: lambda
            intent(in) :: N, lambda

            call checkpoint(3, 'Allocating space for H')
            allocate(H(2**N, 2**N), stat=info)
            if(info/=0) then
                call checkpoint(0, 'Failed to allocate space for H')
                stop 1
            else
                call checkpoint(3, 'Done!')
            end if
        
            H = (0d0, 0d0)

            ! self_int part
            do ii=1, N
                H = H + get_op_on_ijN(N, sigma_z, ii)
                write (string, '(A, A, F6.2, A)') creturn, "Getting Hamiltonian... ", ii*100./N, "% done."
                call checkpoint(3, string, advance='no')
            end do
            call checkpoint(3, string='')
            H = lambda*H

            ! nn_int part
            do ii=1, N-1
                H = H - get_op_on_ijN(N, sigma_x, ii, ii+1)
            end do

        end subroutine

        function get_op_on_ijN(N, op, i, j) result(op_ij)
            integer     :: op(:, :), N, i, j, op_ij(:, :), ii, &
                           identity(ubound(op, 1), ubound(op, 2)), &
                           id_vect(ubound(op, 1)**2)
            optional    :: j
            allocatable :: op_ij
            intent(in)  :: N, op, i, j

            id_vect = 0
            id_vect(1::ubound(op, 1)+1) = 1
            identity = reshape(id_vect, shape(identity))
            allocate(op_ij(1,1))
            op_ij = 1
            if(present(j)) then
                do ii=1, N
                    if(ii==i .or. ii==j) then
                        op_ij = tensor_prod(op_ij, op)
                    else
                        op_ij = tensor_prod(op_ij, identity)
                    end if
                end do
            else
                do ii=1, N
                    if(ii==i) then
                        op_ij = tensor_prod(op_ij, op)
                    else
                        op_ij = tensor_prod(op_ij, identity)
                    end if
                end do
            end if

        end function

        subroutine print_help()
            print *, "Usage: ./exercise07 -N [#subsystems] -l [inter. strength] -k [#eigenval to save to file] -o [output fname]"

        end subroutine print_help
end program