program matrix_mul
    use advanced_ops
    implicit none
    integer*4 :: m, n, i, j, NN, ii, seed, outunit = 100
    real*8, allocatable :: A(:, :), B(:, :), result(:, :)
    real, allocatable :: perf_1(:), perf_2(:), perf_default(:)
    real :: start, finish
    character(len=100) :: ofname
    character*1 :: creturn = achar(13)

    !allocate simple matrices to test the correctness of calculations
    A = reshape( (/ 1, 2, 3, &
                    4, 5, 6 /), &
                (/2,3/), order=(/2,1/) )

    B = reshape( (/ 1., 2., 3., 1., &
                    4., 5., 6., 0.5, &
                    7., 8., 9., 2.5 /), &
                (/3,4/), order=(/2,1/) )

    write (*, '(A/,' // matrixformatter(A) // ')') "Matrix A:", transpose(A)
    write (*, '(A/,' // matrixformatter(B) // ')') "Matrix B:", transpose(B)
    result = matmul1(A, B)
    write (*, '(A/,'//matrixformatter(result)//')') "Result with matmul1:", transpose(result)
    result = matmul2(A, B)
    write (*, '(A/,'//matrixformatter(result)//')') "Result with matmul2:", transpose(result)
    result = matmul(A, B)
    write (*, '(A/,'//matrixformatter(result)//')') "Result with default matmul:", transpose(result)

    !ask for user input before testing the performances
    write (*, '(/A, /A)', advance='no') "Now let's test performances", &
                                       "Insert matrix A dimensions (space, comma or enter separated): "
    read  (*, *) m, n
    write (*, '(A)', advance='no') "Insert matrix B dimensions (space, comma or enter separated): "
    read  (*, *) i, j
    write (*, '(A)', advance='no') "How many operations?: "
    read  (*, *) NN
    write (*, '(A)', advance='no') "Output file name: "
    read  (*, '(A)') ofname
    write (*, *) ofname

    !performance test
    deallocate(A, B, result)
    allocate(A(m, n), B(i, j), result(m, j), perf_1(NN), perf_2(NN), perf_default(NN))

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
        result = matmul(A, B)
        call cpu_time(finish)
        perf_default(ii) = finish - start

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Performance testing, ", ii*100./NN, "% done."
    end do
    write (*, *)

    open(unit = outunit, file=trim(ofname), action='write', status='new')
    write (outunit, *) perf_1
    write (outunit, *) perf_2
    write (outunit, *) perf_default
    close(outunit)

end program