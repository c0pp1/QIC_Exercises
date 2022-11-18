program matrix_mul_perf
    use advanced_ops

    implicit none

    integer*4 :: m, NN, ii, seed, outunit = 100
    real*8, allocatable :: A(:, :), B(:, :), result(:, :)
    real, allocatable :: perf_1(:), perf_2(:), perf_2_new(:), perf_default(:)
    real :: start, finish
    character(len=200) :: ofdir, ofname
    character*1 :: creturn = achar(13)
    character(len=32) :: arg, format
    logical           :: skip_next = .false.

    !> read command line arguments

    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-m', '--dimension')
                call get_command_argument(ii+1, arg)
                read (arg, *) m
                skip_next = .true.
                
            case ('-n', '--iter')
                call get_command_argument(ii+1, arg)
                read (arg, *) NN
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

    !performance test
    allocate(A(m, m), B(m, m), result(m, m), perf_1(NN), perf_2(NN), perf_2_new(NN), perf_default(NN))

    seed = 1
    call random_seed(seed)

    do ii = 1, NN
        call random_number(A)
        call random_number(B)

        call cpu_time(start)
        result = matmul1(A, B)
        call cpu_time(finish)
        perf_1(ii) = finish - start

        call cpu_time(start)
        result = matmul2(A, B)
        call cpu_time(finish)
        perf_2(ii) = finish - start

        call cpu_time(start)
        result = matmul2_new(A, B)
        call cpu_time(finish)
        perf_2_new(ii) = finish - start

        call cpu_time(start)
        result = matmul(A, B)
        call cpu_time(finish)
        perf_default(ii) = finish - start

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Performance testing, ", ii*100./NN, "% done."
    end do
    write (*, *)

    write (format, '(A, I2, A)') "(A, I", floor(log10(real(m)))+1, ", A)'"
    write (ofname, format) trim(ofdir), m, ".txt"
    open(unit = outunit, file=trim(ofname), action='write')
    write (outunit, *) perf_1
    write (outunit, *) perf_2
    write (outunit, *) perf_2_new
    write (outunit, *) perf_default
    close(outunit)

    contains

        subroutine print_help()
            print *, "you called print help"
        end subroutine

end program