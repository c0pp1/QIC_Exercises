program main
    use checkpoint_mod
    use complex_matrix_ops
    use lapack_wrapper
    implicit none

    double precision, parameter :: PI = 4.D0*DATAN(1.D0)
    character(255) :: ofdir = "./exercise06_", arg, format
    character(1), parameter :: creturn = achar(13)
    integer :: ii, d(:), N=0, to_tr_out(:), outunit = 100, iter = 1, info
    allocatable :: d, to_tr_out, user_input
    double complex, allocatable :: sep_state(:, :), gen_state(:), density_matrix(:,:), partial_dm(:,:), user_state(:), test_dm(:,:)
    double precision, allocatable :: perf_sep(:), perf_gen(:), eigenvals(:)
    double precision :: tstart, tstop, user_input(:)
    logical :: skip_next = .false., test = .false., testisdm = .false.
    double complex, pointer :: test_state(:)
    target :: user_state, gen_state

    interface get_general_pure_state
        procedure get_general_pure_stateD
        procedure get_general_pure_stateDN
    end interface

    DB = 4

    ! read command line arguments
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

            case ('-d')
                if(N==0) then
                    call checkpoint(1, 'N must be the first argument')
                    stop
                else if(.not. allocated(d)) then
                    allocate(d(N), to_tr_out(N))
                end if
                call get_command_argument(ii+1, arg)
                read (arg, *, iostat=info) d
                if(info==-1 .or. info==-2) then
                    call checkpoint(0, "The user input doesn't match the system dimension")
                    stop
                end if
                skip_next = .true.

            case ('-tr')                             !< the systems to trace out, must be of same size as d
                if(N==0) then
                    call checkpoint(1, 'N must be the first argument')
                    stop
                else if(.not. allocated(to_tr_out)) then
                    allocate(d(N), to_tr_out(N))
                end if
                call get_command_argument(ii+1, arg)
                read (arg, *, iostat=info) to_tr_out
                if(info==-1 .or. info==-2) then
                    call checkpoint(0, "The user input doesn't match the system dimension")
                    stop
                end if

                skip_next = .true.

            case ('-i')
                call get_command_argument(ii+1, arg)
                read (arg, *) iter
                skip_next = .true.

            case ('-m')
                skip_next = .false.
                testisdm = .true.
                if(allocated(user_state)) then
                    call checkpoint(1, '-m must be put before -t argument')
                    stop
                end if

            case ('-t')
                if(testisdm) then
                    allocate( test_dm(product(d), product(d)), user_input(product(d)**2*2) )
                else 
                    allocate(user_state(product(d)), user_input(product(d)*2))
                end if
                call get_command_argument(ii+1, arg)
                read (arg, *, iostat=info) user_input
                if(info==-1 .or. info==-2) then
                    call checkpoint(0, "The user input doesn't match the system dimension")
                    stop
                end if

                if(testisdm) then
                    test_dm = reshape(cmplx(user_input(1::2), user_input(2::2), 8), shape(test_dm), order=(/2,1/))
                else 
                    user_state = cmplx(user_input(1::2), user_input(2::2), 8)
                    test_state => user_state
                end if
                
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
        call get_separable_state(d(1), N, sep_state)    !< in test mode simplify user inputs using only the first element of d
        call cpu_time(tstop)
        perf_sep(ii) = tstop - tstart
        ! test general pure state
        call cpu_time(tstart)
        call get_general_pure_state(d(1), N, gen_state) !< in test mode simplify user inputs using only the first element of d
        call cpu_time(tstop)
        perf_gen(ii) = tstop - tstart

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Performance testing, ", ii*100./iter, "% done."
    end do
    write(*, *)

    write(format, '(A, I1, A, I1, A)') '(A, A, I', floor(log10(real(d(1))))+1, ', A, I', floor(log10(real(N)))+1, ')'
    write(ofdir, trim(format)) trim(ofdir), 'd', d(1), 'N', N
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

    if(iter==1) then
        if(size(d) /= size(to_tr_out)) then
            call checkpoint(1, 'the arrays describing dimension of subsystems and which systems to trace out' //&
                                ' must have the same size! Got:', vector=real((/size(d), size(to_tr_out)/), 8))
            stop
        end if
        call get_general_pure_state(d, gen_state)   !< getting general pure state compatible with user inputs
        if(.not. test) test_state => gen_state
        call checkpoint(3, 'Getting density matrix for ' // trim(merge('test   ', 'general', test)) // ' state')
        if(testisdm) then
            if(allocated(density_matrix)) deallocate(density_matrix)
            allocate(density_matrix(ubound(test_dm, 1), ubound(test_dm, 2)))
            density_matrix = test_dm
        else
            call get_density_matrix(test_state, density_matrix)
        end if
        call checkpoint(3, 'Trace of density matrix: ', scalar = real(.Tr. density_matrix, 8))
        call checkpoint(4, 'Density matrix:', cmatrix = density_matrix)

        allocate(eigenvals(size(density_matrix, 1)))
        call ev_wrap(A=density_matrix, W=eigenvals, INFO=info)
        call checkpoint(4, 'Density matrix eigenvalues:', vector = eigenvals)

        if(size(test_state) /= product(d) .and. (.not. testisdm)) then
            call checkpoint(1, 'Dimension of test state is not coherent with user given dimensions', &
                            vector=real((/size(test_state), product(d)/), 8))
            stop
        end if
        deallocate(eigenvals)
        allocate(eigenvals(product(pack(d, to_tr_out==0))))

        if(testisdm) then
            if(allocated(density_matrix)) deallocate(density_matrix)
            allocate(density_matrix(ubound(test_dm, 1), ubound(test_dm, 2)))
            density_matrix = test_dm
        else
            call get_density_matrix(test_state, density_matrix)  !< has to be recomputed as ev_wrap destroyed lower triangular part
        end if


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

        call partial_trace2(density_matrix, to_tr_out, d, partial_dm)
        call checkpoint(4, 'Partial density matrix:', cmatrix = partial_dm)
        call ev_wrap(A=partial_dm, W=eigenvals, INFO=info)
        if(info/=0) call checkpoint(0, 'Failed in getting partial_dm eigenvalues')
        call checkpoint(4, 'partial density matrix eigenvalues:', vector = eigenvals)
        call checkpoint(3, 'Trace of partial density matrix:', scalar = real(.Tr. partial_dm, 8))

        partial_dm = cmplx(get_trace_matrix((/0,1/), (/2, 2/), 2), kind=8)

        call checkpoint(4, 'Trace matrix:', cmatrix=partial_dm)
        call checkpoint(4, 'Tensor product:', cmatrix=cmplx(tensor_prod(reshape((/1,0/), (/1, 2/)),&
                        reshape((/1, 0, 0, 1/), (/2, 2/))), kind=8))

    end if
                
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
        subroutine get_general_pure_stateDN(d, N, state)
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

        subroutine get_general_pure_stateD(d, state)
            integer :: d(:), N, info=0
            double complex, allocatable :: state(:)
            double precision, allocatable :: rel_phases(:), amps(:)
            intent(in) :: d
            intent(out) :: state

            N = product(d)
            call checkpoint(3, 'Allocating arrays for general pure state computation...')
            if(.not. allocated(state)) allocate(state(N), rel_phases(N-1), amps(N), stat=info)
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

        !> more general partial trace
        !!
        !! given the density matrix, a vector of zeros with ones on the indeces of the
        !! subsystems to be traced out and the vector with dimensions of the subsystems,
        !! calculate the partial trace
        !!
        !! @param[in]   density_matrix  the matrix on which to perform the partial trace
        !! @param[in]   systems a vector of zeros with ones on the indexes of the subsystems to be traced out
        !! @param[in]   dims    vector with dimensions of the subsystems
        !! @param[out]  partial_dm  the partial density matrix
        !!
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

        !> get the matrix to trace out the el-th state
        !!
        !! given the subsystems to trace out and the dimensions of all subsystems,
        !! computes the matrix for tracing out the el-th state, where el indicates
        !! the state among the possible state of a system built on the subsystems
        !! to trace out.
        !!
        !! @param[in]   systems a vector of zeros with ones on the indexes of the subsystems to be traced out
        !! @param[in]   dims    vector with dimensions of the subsystems
        !! @param[in]   el  the index of state to be traced out
        !! @return      tr_mat  the matrix to perform the trace
        !!
        function get_trace_matrix(systems, dims, el) result(tr_mat)
            integer, intent(in) :: systems(:), dims(:), el
            integer, allocatable :: tr_mat(:,:), identity(:), tmp(:,:), tmp_sys(:,:)
            integer :: ii
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

        !> tensor product of bidimensional arrays
        !!
        !! given two bidimensional arrays calculate the tensor product.
        !! The result is a bidiminesional array of size:
        !! (size(arr1, 1)*size(arr2, 1), size(arr1, 2)*size(arr2, 2))
        !!
        !! @param[in]   arr1    the first array
        !! @param[in]   arr2    the second array
        !! @result      arr3    the resulting tensor product
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
            print *, "Usage: ./exercise06 -N [# subsystems] -d [array of dims of subsys.]" //& 
                     " -tr [array of zeros with ones on the entryes of subsystems to be traced out]" // &
                     " -t [eventual state to be tested] -o [output directory] -i [#iterations for performance test]"
        end subroutine print_help

end program