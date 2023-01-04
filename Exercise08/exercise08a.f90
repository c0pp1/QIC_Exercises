module RSRG
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    use hamiltonian
    implicit none

    contains

        subroutine rsrg_step(HN, A, B)
            double complex   :: HN(:, :), H2N(:, :), P(:, :), &
                                A(:, :), B(:, :), id(ubound(HN, 1), ubound(HN, 2)), temp(:, :)
            integer          :: id_vect(ubound(HN, 1)**2), M, info
            double precision :: eigenvals(ubound(HN, 1)**2)
            intent(inout)    :: HN, A, B
            allocatable      :: H2N, P, temp

            call checkpoint(5, 'Input shape:', vector=real(shape(HN), 8))
            M = ubound(HN, 1)

            call checkpoint(5, 'Allocating space for computations...')
            allocate(H2N(M**2, M**2), P(M**2, M**2), temp(M**2, M), stat=info)
            if(info/=0) then
                call checkpoint(1, 'Failed to allocate space for two systems computations')
                stop
            end if

            if(ubound(A, 1) /= M .or. ubound(A, 1) /= ubound(B, 1)) then
                call checkpoint(1, 'Hamiltonian and interaction terms of subsystems must have same dimensions!')
                stop
            end if

            call checkpoint(5, 'Getting identity...')
            id_vect = 0
            id_vect(1::M+1) = 1
            id = reshape(id_vect, shape(HN))

            call checkpoint(5, 'Getting Hamiltonian of doubled system...')
            H2N = tensor_prod(HN, id)       !< hamiltonian on left system
            H2N = H2N + tensor_prod(id, HN) !< hamiltonian on right system

            H2N = H2N + tensor_prod(A, B)                  !< add interaction

            ! get projector
            call checkpoint(5, 'Getting projector...')
            P = H2N
            call ev_wrap(JOBZ='V', A=P, W=eigenvals, INFO=info)
            if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')
            call checkpoint(5, 'Minor eigenvals:', vector=eigenvals(1:M))

            ! get operators for next step
            call checkpoint(5, 'Projecting...')
            HN = matmul(.Adj. P(:, 1:M), matmul(H2N, P(:, 1:M))) / 2.
            call checkpoint(5, 'HN done')
            A  = matmul(.Adj. P(:, 1:M), matmul(tensor_prod(id, A), P(:, 1:M))) / sqrt(2.)
            call checkpoint(5, 'A done')
            B  = matmul(.Adj. P(:, 1:M), matmul(tensor_prod(B, id), P(:, 1:M))) / sqrt(2.)
            call checkpoint(5, 'B done')
        end subroutine

        subroutine rsrg_init(HN, A, B, N, int_op, lambda)
            double complex   :: HN(:, :), A(:, :), B(:, :)
            integer          :: int_op(:, :), N, &
                                id_vect(:), id(:, :)
            double precision :: lambda
            allocatable      :: HN, A, B, id_vect, id
            intent(in)       :: N, int_op, lambda
            intent(out)      :: HN, A, B

            allocate(id_vect(2**(2*(N-1))), id(2**(N-1), 2**(N-1)))
            if(.not. allocated(HN)) allocate(HN(2**N, 2**N))
            id_vect = 0
            id_vect(1::2**(N-1)+1) = 1
            id = reshape(id_vect, shape(id))

            ! get interaction terms
            A = tensor_prod(id, int_op)
            B = tensor_prod(int_op, id)

            call get_Hamiltonian(N, lambda, HN)

        end subroutine


end module RSRG

module hamiltonian
    use checkpoint_mod
    use complex_matrix_ops
    implicit none

    integer, parameter :: sigma_z(2, 2) = reshape((/1, 0, 0, -1/), (/2,2/)), &
                          sigma_x(2, 2) = reshape((/0, 1, 1, 0/), (/2,2/))
    character(1)       :: creturn = char(13)
    character(40)      :: string
    private           :: string, creturn

    contains

    subroutine get_Hamiltonian(N, lambda, H)
        integer :: N, ii, info
        double complex, allocatable, intent(out) :: H(:, :)
        double precision :: lambda
        intent(in) :: N, lambda
    
        call checkpoint(4, 'Allocating space for H')
        allocate(H(2**N, 2**N), stat=info)
        if(info/=0) then
            call checkpoint(0, 'Failed to allocate space for H')
            stop 1
        else
            call checkpoint(4, 'Done!')
        end if
    
        H = (0d0, 0d0)
    
        ! self_int part
        do ii=1, N
            H = H + get_op_on_ijN(N, sigma_z, ii)
            write (string, '( A, F6.2, A, A)') "Getting Hamiltonian... ", ii*100./N, "% done.", creturn
            call checkpoint(3, string, advance='no')
        end do
        call checkpoint(3, string='Getting Hamiltonian... 100.00% done.')
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
end module

program main
    use checkpoint_mod
    use complex_matrix_ops
    use RSRG
    implicit none

    integer          :: N = 2, ii, info, outunit = 100, iter
    double complex   :: HN(:, :), A(:, :), B(:, :)
    double precision :: lambda = 1.d0, eigenvals(:)
    allocatable      :: HN, A, B, eigenvals
    character(128)   :: arg, ofname='test', string
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

            case ('-i', '--iter')
                call get_command_argument(ii+1, arg)
                read (arg, *) iter
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

    call checkpoint(4, 'Initializing init Hamiltonian with #subsys:', scalar=real(N, 8))
    call rsrg_init(HN, A, B, N, sigma_x, lambda)
    call checkpoint(4, 'Done!')

    call checkpoint(4, 'Iterating...')
    do ii=1, iter
        call rsrg_step(HN, A, B)
        write (string, '(A, F6.2, A, A)') "Expanding system... ", ii*100./iter, "% done.", creturn
        call checkpoint(3, string, advance='no')
    end do
    call checkpoint(3, 'Expanding system... 100.00% done.')

    call checkpoint(4, 'Allocating space for eigenvals')
    allocate(eigenvals(2**N), stat=info)
    if(info/=0) then
        call checkpoint(0, 'Failed to allocate space for eigenvals')
        stop 1
    else
        call checkpoint(4, 'Done!')
    end if

    call ev_wrap(A=HN, W=eigenvals, INFO=info)
    if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')

    ! save results to file, appending N, lambda, and ground state energy
    call checkpoint(4, "Writing to file '" // trim(ofname) // "'")
    inquire(file=trim(ofname), exist=exist)
    if (exist) then
        open(outunit, file=trim(ofname), status="old", position="append", action="write")
    else
        open(outunit, file=trim(ofname), status="new", action="write")
    end if

    write(outunit, *) N, iter, lambda, eigenvals(1)

    call checkpoint(4, 'All done!')

contains
    
    subroutine print_help()
        print *, "Usage: ./exercise08a -N [#subsystems] -l [inter. strength] - [#iterations] -o [output fname]"

    end subroutine print_help

end program main


