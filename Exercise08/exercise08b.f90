module IDMRG
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    use hamiltonian
    implicit none

    contains

        subroutine idmrg_step(H1, H2, H4, H12, H23, H34)
            double complex   :: H1(:, :), H2(:, :), H4(:, :), H12(:, :), H23(:, :), H34(:, :), &
                                H_tmp(:, :), PDM1(:, :), PDM2(:, :), DM(:, :), id(ubound(H1, 1), ubound(H1, 2)), tmp(:,:)
            integer          :: id_vect(ubound(H1, 1)**2), M, info
            double precision :: eigenvals((ubound(H1, 1)*2)**2), eigenvals1((ubound(H1, 1)*2)), eigenvals2((ubound(H1, 1)*2))
            intent(inout)    :: H1, H2, H4, H12, H23, H34
            allocatable      :: PDM1, PDM2, DM, H_tmp, tmp

            M = ubound(H1, 1)

            call checkpoint(5, 'Allocating space for computations...')
            allocate(H_tmp((M*2)**2, (M*2)**2), DM((M*2)**2, (M*2)**2), stat=info)
            if(info/=0) then
                call checkpoint(1, 'Failed to allocate space for two systems computations')
                stop 1
            end if

            call checkpoint(5, 'Getting identity...')
            id_vect = 0
            id_vect(1::ubound(id, 1)+1) = 1
            id = reshape(id_vect, shape(id))

            call checkpoint(5, 'Getting total Hamiltonian of the system...')
            H_tmp = tot_hamiltonian(H1, H2, H4, H12, H23, H34)

            call checkpoint(5, 'Getting groundstate...')
            call ev_wrap(JOBZ='V', A=H_tmp, W=eigenvals, INFO=info)
            if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')

            call checkpoint(5, 'Getting DM...')
            call get_density_matrix(H_tmp(:,1), DM)
            call checkpoint(5, 'DM size:', vector=real((/ubound(DM, 1), ubound(DM, 2)/), 8))
            deallocate(H_tmp)

            call checkpoint(5, 'Getting partial DM for left part...')
            call partial_trace(DM, (/0,0,1,1/), (/M, 2, 2, M/), PDM1)
            call checkpoint(5, 'partial DM size:', vector=real((/ubound(PDM1, 1), ubound(PDM1, 2)/), 8))

            call checkpoint(5, 'Getting partial DM for right part...')
            call partial_trace(DM, (/1,1,0,0/), (/M, 2, 2, M/), PDM2)
            call checkpoint(5, 'partial DM size:', vector=real(shape(PDM2), 8))
            deallocate(DM)

            ! get projector
            call checkpoint(5, 'Getting projectors...')
            call ev_wrap(JOBZ='V', A=PDM1, W=eigenvals1, INFO=info)
            if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')
            call checkpoint(5, 'Major eigenvals L:', vector=eigenvals1(M*2:M*2-M+1:-1))

            call ev_wrap(JOBZ='V', A=PDM2, W=eigenvals2, INFO=info)
            if(info/=0) call checkpoint(1, 'Failed to get eigenvalues')
            call checkpoint(5, 'Major eigenvals R:', vector=eigenvals2(M*2:M*2-M+1:-1))

            ! get operators for next step
            call checkpoint(5, 'Projecting...')
            H1 = matmul(.Adj. PDM1(:, M*2:M*2-M+1:-1), &
                        matmul((H1 .tprod. id(:2, :2)) + (id .tprod. H2) + (H12 .tprod. cmplx(sigma_x, kind=8)), &
                               PDM1(:, M*2:M*2-M+1:-1)))
            call checkpoint(5, 'H1 done')

            call checkpoint(5, 'H4 shape:', vector=real(shape(H4), 8))
            tmp = matmul((id(:2, :2) .tprod. H4) + (H2 .tprod. id) + (cmplx(sigma_x, kind=8) .tprod. H34), &
                            PDM2(:, M*2:M*2-M+1:-1))
            call checkpoint(5, 'tmp shape:', vector=real(shape(tmp), 8))
            H4 = matmul(.Adj. PDM2(:, M*2:M*2-M+1:-1), tmp)
                        
            call checkpoint(5, 'H4 done')

            H12  = matmul(.Adj. PDM1(:, M*2:M*2-M+1:-1), &
                          matmul(H12 .tprod. id(:2,:2), PDM1(:, M*2:M*2-M+1:-1)))
            call checkpoint(5, 'H12 done')

            H34  = matmul(.Adj. PDM2(:, M*2:M*2-M+1:-1), &
                          matmul( id(:2,:2) .tprod. H34, PDM2(:, M*2:M*2-M+1:-1)))
            call checkpoint(5, 'H34 done')
        end subroutine

        subroutine idmrg_init(H1, H2, H12, H23, H34, N, int_op, lambda)
            double complex   :: H1(:, :), H2(:, :), H12(:, :), H23(:, :), H34(:, :)
            integer          :: int_op(:, :), N, info
            double precision :: lambda
            allocatable      :: H1, H2, H12, H23, H34
            intent(in)       :: N, int_op, lambda
            intent(out)      :: H1, H2, H12, H23, H34

            call checkpoint(4, 'Allocating space for initialization')
            allocate(H2(2, 2), H12(2**N, 2**N), &
                     H34(2**N, 2**N), H23(2**2, 2**2), stat=info)
            if(info/=0) then
                call checkpoint(0, 'Failed to allocate space for initialization')
                stop 1
            else
                call checkpoint(4, 'Done!')
            end if

            call get_Hamiltonian(N, lambda, H1, int_op)

            H2 = lambda*sigma_z
        
            H12 = get_op_on_ijN(N, int_op, N)   !< these are 'awaiting' interaction, remember to .tprod. sigma_x
            H34 = get_op_on_ijN(N, int_op, 1)   !< these are 'awaiting' interaction, remember to sigma_x .tprod.

            H23 = get_op_on_ijN(2, int_op, 1, 2)

        end subroutine

        function tot_hamiltonian(H1, H2, H4, H12, H23, H34) result(H_tot)
            double complex :: H1(:, :), H2(:, :), H4(:, :), H12(:, :), H23(:, :), H34(:, :), H_tot(:, :), id(:, :)
            integer        :: dim1, dim2, id_vect(:)
            intent(in)     :: H1, H2, H4, H12, H23, H34
            allocatable    :: H_tot, id_vect, id

            dim1 = ubound(H1, 1)
            dim2 = ubound(H2, 1)

            allocate(H_tot((dim1*dim2)**2, (dim1*dim2)**2), id(dim1*(dim2**2), dim1*(dim2**2)), id_vect((dim1*(dim2**2))**2))

            id_vect = 0
            id_vect(1::ubound(id, 1)+1) = 1
            id = reshape(id_vect, shape(id))

            H_tot = (H1 .tprod. id) + (id .tprod. H4) + &
                    (id(:dim1,:dim1) .tprod. (H2 .tprod. id(:dim1*dim2,:dim1*dim2))) + &
                    (id(:dim1*dim2,:dim1*dim2) .tprod. (H2 .tprod. id(:dim1,:dim1))) + &
                    (H12 .tprod. (cmplx(sigma_x, kind=8) .tprod. id(:dim1*dim2,:dim1*dim2))) + &
                    (id(:dim1*dim2,:dim1*dim2) .tprod. (cmplx(sigma_x, kind=8) .tprod. H34)) + &
                    (id(:dim1,:dim1) .tprod. (H23 .tprod. id(:dim1,:dim1)))

        end function tot_hamiltonian

end module IDMRG

module hamiltonian
    use checkpoint_mod
    use complex_matrix_ops
    implicit none

    integer, parameter :: sigma_z(2, 2) = reshape((/1, 0, 0, -1/), (/2,2/)), &
                          sigma_x(2, 2) = reshape((/0, 1, 1, 0/), (/2,2/))
    character(1)       :: creturn = char(13)
    character(40)      :: string
    private            :: string, creturn

    contains

    subroutine get_Hamiltonian(N, lambda, H, int_op)
        integer :: N, ii, info, int_op(:, :), op(:,:)
        double complex, allocatable, intent(out) :: H(:, :)
        double precision :: lambda
        intent(in) :: N, lambda, int_op
        optional :: int_op
        allocatable :: op

        if(present(int_op)) then
            op = int_op
        else
            op = sigma_x
        end if
    
        call checkpoint(4, 'Allocating space for H')
        if(.not. allocated(H)) then 
            allocate(H(2**N, 2**N), stat=info)
            if(info/=0) then
                call checkpoint(0, 'Failed to allocate space for H')
                stop 1
            else
                call checkpoint(4, 'Done!')
            end if
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
            H = H + get_op_on_ijN(N, op, ii, ii+1)
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
                    op_ij = op_ij .tprod. op
                else
                    op_ij = op_ij .tprod. identity
                end if
            end do
        else
            do ii=1, N
                if(ii==i) then
                    op_ij = op_ij .tprod. op
                else
                    op_ij = op_ij .tprod. identity
                end if
            end do
        end if
    
    end function
end module

program main
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    use IDMRG
    implicit none

    integer          :: N = 2, ii, info, outunit = 100, iter
    double complex   :: H1(:, :), H2(:, :), H4(:, :), H12(:, :), H23(:, :), H34(:, :)
    double precision :: lambda = 1.d0, eigenvals(:)
    allocatable      :: H1, H2, H4, H12, H23, H34, eigenvals
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
    call idmrg_init(H1, H2, H12, H23, H34, N, sigma_x, lambda)
    allocate(H4(2**N, 2**N))
    H4 = H1
    call checkpoint(5, 'Init Hamiltonian:', cmatrix=tot_hamiltonian(H1, H2, H4, H12, H23, H34))

    call checkpoint(4, 'Iterating...')
    do ii=1, iter
        call idmrg_step(H1, H2, H4, H12, H23, H34)
        write (string, '(A, F6.2, A, A)') "Expanding system... ", ii*100./iter, "% done.", creturn
        call checkpoint(3, string, advance='no')
    end do
    call checkpoint(3, 'Expanding system... 100.00% done.')

    call checkpoint(4, 'Allocating space for eigenvals')
    allocate(eigenvals(2**(N+1)**2), stat=info)
    if(info/=0) then
        call checkpoint(0, 'Failed to allocate space for eigenvals')
        stop 1
    else
        call checkpoint(4, 'Done!')
    end if

    call ev_wrap(A=tot_hamiltonian(H1, H2, H4, H12, H23, H34), W=eigenvals, INFO=info)
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
        print *, "Usage: ./exercise08b -N [#subsystems] -l [inter. strength] -i [#iterations] -o [output fname]"

    end subroutine print_help

end program main


