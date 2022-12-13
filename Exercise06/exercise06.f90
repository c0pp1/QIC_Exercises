program main
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    implicit none

    double precision, parameter :: PI = 4.D0*DATAN(1.D0)
    character(64) :: ofdir = "./exercise06_", arg, format
    character(1), parameter :: creturn = achar(13)
    integer :: ii, outunit = 100, d, N, iter = 1, M=2, info, which
    allocatable : d(:), N(:)
    double complex, allocatable :: sep_state(:, :), gen_state(:), density_matrix(:,:), partial_dm(:,:), user_state(:)
    double precision, allocatable :: perf_sep(:), perf_gen(:), eigenvals(:)
    double precision :: tstart, tstop, user_input(8)
    logical :: skip_next = .false., test = .false.
    double complex, pointer :: test_state(:)
    target :: user_state, gen_state

    
    DB = 4

    ! read command line arguments
    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-d')
                call get_command_argument(ii+1, arg)
                read (arg, *) d
                skip_next = .true.

            case ('-N')
                call get_command_argument(ii+1, arg)
                read (arg, *) N
                skip_next = .true.

            case ('-M')                             !< the system to trace out tot. dim
                call get_command_argument(ii+1, arg)
                read (arg, *) M
                skip_next = .true.

            case ('-i')
                call get_command_argument(ii+1, arg)
                read (arg, *) iter
                skip_next = .true.

            case ('-w', '--which')
                call get_command_argument(ii+1, arg)
                read (arg, *) which
                skip_next = .true.
                if(which/=1 .and. which/=2) then
                    call checkpoint(1, 'w must be 1 or 2, setting to 1')
                    which = 1
                end if

            case ('-t')
                allocate(user_state(4))
                call get_command_argument(ii+1, arg)
                read (arg, *) user_input
                user_state = cmplx(user_input(1::2), user_input(2::2), 8)
                test_state => user_state
                test = .true.
                skip_next = .true.

            case ('-o', '--odirname')
                call get_command_argument(ii+1, arg)
                read (arg, '(A)') ofdir
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

    allocate(perf_gen(iter), perf_sep(iter))

    do ii=1, iter
        ! test separable state normalization
        call cpu_time(tstart)
        call get_separable_state(d, N, sep_state)
        call cpu_time(tstop)
        perf_sep(ii) = tstop - tstart
        ! test general pure state
        call cpu_time(tstart)
        call get_general_pure_state(d, N, gen_state)
        call cpu_time(tstop)
        perf_gen(ii) = tstop - tstart

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Performance testing, ", ii*100./iter, "% done."
    end do
    write(*, *)

    write(format, '(A, I1, A, I1, A)') '(A, A, I', floor(log10(real(d)))+1, ', A, I', floor(log10(real(N)))+1, ')'
    write(ofdir, trim(format)) trim(ofdir), 'd', d, 'N', N
    open(unit=outunit, file=trim(ofdir) // "_perf.txt", action='write')
    write(outunit, *) perf_sep
    write(outunit, *) perf_gen
    close(outunit)

    open(unit=outunit, file=trim(ofdir) // "_sep.txt", action='write')
    call write_cmplx(outunit, sep_state)
    close(outunit)
    open(unit=outunit, file=trim(ofdir) // "_pure.txt", action='write')
    call write_cmplx(outunit, gen_state)
    close(outunit)

    if(.not. test) test_state => gen_state
    call checkpoint(3, 'Getting density matrix for ' // trim(merge('test   ', 'general', test)) // ' state')
    call get_density_matrix(test_state, density_matrix)
    call checkpoint(3, 'Trace of density matrix: ', scalar = real(.Tr. density_matrix, 8))
    call checkpoint(4, 'Density matrix:', cmatrix = density_matrix)

    allocate(eigenvals(size(density_matrix, 1)))
    call ev_wrap(A=density_matrix, W=eigenvals, INFO=info)
    call checkpoint(4, 'Density matrix eigenvalues:', vector = eigenvals)

    if((mod(d**N, M) /= 0) .or. (mod(4, M) /= 0)) then
        call checkpoint(1, 'M must be a factor of total state dimension')
        stop
    end if

    deallocate(eigenvals)
    allocate(eigenvals(merge(4/M,d**N/M, test)))

    call get_density_matrix(test_state, density_matrix) !< has to be recomputed as ev_wrap destroyed lower triangular part

    !call checkpoint(3, 'getting partial density matrices (trace out system 1 and then 2) with dimensions:', &
    !                vector = real((/merge(4/M,d**N/M, test), M/), 8))
    !call partial_trace(density_matrix, which, M, partial_dm)
    !call checkpoint(4, 'Partial density matrix of system:', scalar = real(mod(which, 2)+1, 8), cmatrix = partial_dm)
    !call ev_wrap(A=partial_dm, W=eigenvals, INFO=info)
    !if(info/=0) call checkpoint(0, 'Failed in getting partial_dm eigenvalues')
    !call checkpoint(4, 'partial density matrix eigenvalues:', vector = eigenvals)
    !call checkpoint(3, 'Trace of partial density matrix:', scalar = real(.Tr. partial_dm, 8))

    !deallocate(partial_dm)
    !deallocate(eigenvals)
    !allocate(eigenvals(M))
    !call partial_trace(density_matrix, mod(which, 2)+1, ubound(density_matrix, 1)/M, partial_dm)
    !call checkpoint(4, 'Partial density matrix of system:', scalar=real(which, 8), cmatrix = partial_dm)
    !call ev_wrap(A=partial_dm, W=eigenvals, INFO=info)
    !if(info/=0) call checkpoint(0, 'Failed in getting partial_dm eigenvalues')
    !call checkpoint(4, 'partial density matrix eigenvalues:', vector = eigenvals)
    !call checkpoint(3, 'Trace of partial density matrix:', scalar = real(.Tr. partial_dm, 8))

    call partial_trace2(density_matrix, (/1, 0/), (/4, 4/), partial_dm)
    call checkpoint(4, 'Partial2 density matrix:', scalar = real(mod(which, 2)+1, 8), cmatrix = partial_dm)
    call ev_wrap(A=partial_dm, W=eigenvals, INFO=info)
    if(info/=0) call checkpoint(0, 'Failed in getting partial_dm eigenvalues')
    call checkpoint(4, 'partial density matrix eigenvalues:', vector = eigenvals)
    call checkpoint(3, 'Trace of partial density matrix:', scalar = real(.Tr. partial_dm, 8))

    partial_dm = cmplx(get_trace_matrix((/0,1/), (/2, 2/), 2), kind=8)

    call checkpoint(4, 'Trace matrix:', cmatrix=partial_dm)
    call checkpoint(4, 'Tensor product:', cmatrix=cmplx(tensor_prod(reshape((/1,0/), (/1, 2/)), reshape((/1, 0, 0, 1/), (/2, 2/)))&
                    , kind=8))


                
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
            integer :: d, N, info=0
            double complex, allocatable :: state(:, :)
            double precision, allocatable :: rel_phases(:, :), amps(:, :)
            intent(in) :: d, N
            intent(out) :: state

            call checkpoint(3, 'Allocating arrays for separable state computation...')
            if(.not. allocated(state)) allocate(state(d, N), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for separable state computation')
            allocate(rel_phases(d-1, N), amps(d, N), stat=info)

            ! generate random parameters: 
            ! angle between 0-2pi
            call random_number(rel_phases)
            rel_phases = rel_phases*2*PI
            ! amplitude between 0-1
            call random_number(amps)

            state(2:, :) = cdexp(cmplx(0., -rel_phases, kind=8))*amps(2:, :)

            ! calculate first coefficients with zero phase
            state(1, :) = amps(1,:)
            ! renormalize
            state = state / sqrt(spread(sum(state*conjg(state), dim=1), 1, d))
            call checkpoint(4, 'State shape after renormalization:', vector=real(shape(state), kind=8))

        end subroutine

        !> subroutine to get a state describing N-body, interacting, non-separable
        !! pure state. Randomly initialized
        !!
        !! the state is represented as a d**N vector, where the i-th element is the
        !! coefficient of the state with i in basis d as local states
        !!
        !! @param[in]   N   #bodies in the system
        !! @param[out]  state   the randomly initialized state
        !!
        subroutine get_general_pure_state(d, N, state)
            integer :: d, N, info=0
            double complex, allocatable :: state(:)
            double precision, allocatable :: rel_phases(:), amps(:)
            intent(in) :: d, N
            intent(out) :: state

            call checkpoint(3, 'Allocating arrays for general pure state computation...')
            if(.not. allocated(state)) allocate(state(d**N), rel_phases(d**N-1), amps(d**N), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for general pure state computation')

            ! generate random parameters: 
            ! angle between 0-2pi
            call random_number(rel_phases)
            rel_phases = rel_phases*2*PI
            ! amplitude between 0-1
            call random_number(amps)

            ! calculate coefficients setting first phase to zero
            state = cdexp(cmplx(0., (/0._8, -rel_phases/), kind=8))*amps
            ! renormalize
            state = state / sqrt(sum(state*conjg(state)))

        end subroutine

        subroutine get_density_matrix(state, density_matrix)
            double complex, allocatable :: density_matrix(:,:)
            double complex :: state(:)
            integer :: info=0, m
            intent(in) :: state
            intent(out) :: density_matrix

            m = ubound(state, 1)
            call checkpoint(3, 'Allocating array for density matrix...')
            if(.not. allocated(density_matrix)) allocate(density_matrix(m, m), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for general pure state computation')

            call checkpoint(5, 'real state:', vector=real(state))
            
            density_matrix = spread(state, 2, m)* .Adj. spread(state, 2, m)
            call checkpoint(5, 'spread state:', vector=real(density_matrix(1,:)))

        end subroutine

        !> partial trace over the i-th system 
        !!
        !! given the density matrix suppose that it represents two systems
        !! with equal size s.t. if shape(density_matrix) = d**N, the output
        !! is always shape(output) = d**(N/2). Thus, i must be 1 or 2.
        !!
        !! @param[in]   density_matrix the matrix on which to perform the trace
        !! @param[in]   i   the system to trace out
        !! @param[in]   M   the system to trace out dimension
        !! @param[out]  partial_dm  the partial density matrix
        !!
        subroutine partial_trace(density_matrix, i, M, partial_dm)
            integer :: i, ii, jj, info=0, N, M
            double complex, allocatable :: density_matrix(:, :), partial_dm(:, :)
            intent(in) :: i, density_matrix, M
            intent(out) :: partial_dm

            N = ubound(density_matrix, 1) / M                                   !< dimension of the output
            if(.not. allocated(partial_dm)) allocate(partial_dm(N, N), stat=info)
            if(info/=0) call checkpoint(0, 'Failed to allocate space for partial trace computation')

            partial_dm = (0._8, 0._8)
            select case(i)
            case(1)
                do ii=1, N
                    do jj=1, N
                        partial_dm(ii, jj) = .Tr. density_matrix((ii-1)*M+1:ii*M, (jj-1)*M+1:jj*M)
                    end do
                end do
            case(2)
                do ii=1, N
                    do jj=1, N
                        partial_dm(ii, jj) = .Tr. density_matrix(mod(ii-1, M)+1::M, mod(jj-1, M)+1::M)
                    end do
                end do
            end select

        end subroutine

        subroutine partial_trace2(density_matrix, systems, dims, partial_dm)
            double complex :: density_matrix(:, :), partial_dm(:, :), trace_mat(:,:)
            integer :: systems(:), dims(:), ii, m
            allocatable :: partial_dm, trace_mat
            intent(in) :: density_matrix, systems, dims
            intent(out) :: partial_dm

            m = product(pack(dims, systems==0))
            if(allocated(partial_dm)) deallocate(partial_dm)
            allocate(partial_dm(m, m))
            partial_dm = cmplx(0._8, kind=8)
            do ii=1, product(pack(dims, systems==1))
        
                trace_mat = cmplx(get_trace_matrix(systems, dims, ii), kind=8)
                call checkpoint(5, 'Trace mat at iter:', scalar=real(ii, 8), cmatrix=cmplx(trace_mat, kind=8))
                call checkpoint(5, 'partial mat at iter:', scalar=real(ii, 8), cmatrix=cmplx(partial_dm, kind=8))
                partial_dm = partial_dm + matmul(trace_mat, matmul(density_matrix, transpose(trace_mat)))
        
            end do
            
        end subroutine

        function get_trace_matrix(systems, dims, el) result(tr_mat)
            integer, intent(in) :: systems(:), dims(:), el
            integer, allocatable :: tr_mat(:,:), identity(:), tmp(:,:), tmp_sys(:,:)
            integer :: ii, jj, kk
            logical :: mask(ubound(dims, 1))

            ! initialize tr_matrix
            allocate(tr_mat(1, 1))
            tr_mat = 1
            mask = (systems==1)

            do ii=1, size(systems)
                allocate(tmp(ubound(tr_mat, 1), ubound(tr_mat, 2)))
                tmp = tr_mat
                deallocate(tr_mat)

                if(systems(ii)==0) then
                    ! allocate space for result
                    allocate(tr_mat(ubound(tmp, 1)*dims(ii), ubound(tmp, 2)*dims(ii)))
                    ! create identity for the i-th system
                    allocate(identity(dims(ii)**2))
                    identity = 0
                    identity(1::dims(ii)+1) = 1
                    
                    tr_mat = tensor_prod(tmp, reshape(identity, (/dims(ii), dims(ii)/)))
                    deallocate(identity)

                else
                    mask(findloc(mask, .true.)) = .false.
                    ! allocate space for result
                    allocate(tr_mat(ubound(tmp, 1), ubound(tmp, 2)*dims(ii)))
                    allocate(tmp_sys(1, dims(ii)))
                    tmp_sys = 0
                    tmp_sys(1, mod((el-1)/product(pack(dims, mask)), dims(ii))+1) = 1

                    tr_mat = tensor_prod(tmp, tmp_sys)
                    deallocate(tmp_sys)
                end if

                deallocate(tmp)
            end do
        end function

        function tensor_prod(arr1, arr2) result(arr3)
            integer :: arr1(:,:), arr2(:,:), arr3(:,:), ii, jj, m1, n1, m2, n2
            allocatable :: arr3

            m1 = size(arr1, 1)
            n1 = size(arr1, 2)
            m2 = size(arr2, 1)
            n2 = size(arr2, 2)
            allocate(arr3(m1*m2, n1*n2))

            do jj=1, n1
                do ii=1, m1
                    arr3((ii-1)*m2+1:ii*m2, (jj-1)*n2+1:jj*n2) = arr1(ii, jj)*arr2
                end do
            end do

        end function

        subroutine print_help()
            print *, "Usage: ./exercise06 -m [#space steps (with prime factors 2,3,5,7)] -n [#time steps] " // &
                     "-x [xlims, comma separated without blank spaces] -t [time T] -w [frequency]"
        end subroutine print_help

end program