program main 
    use complex_matrix
    use checkpoint_mod
    use advanced_ops

    type(cmatrix)         :: matrix
    type(cmatrix)         :: adjoint
    integer, dimension(2) :: dims = (/1,2/), dims2 = (/4,3/)
    integer               :: counter1, unit = 100
    double complex, allocatable :: elems(:, :)
    logical               :: debug = .True.
    logical               :: error
    double complex        :: tr, det

    ! initialize matrix elements
    allocate(elems(dims2(1), dims2(2)))
    elems = reshape((/(2d0, 1.), (2d0, 3.), (4d0, 5.), (2.5d0, 3.2)/), shape=dims2, &
                    pad=(/(4d0, 3.), (7d0, 2.2), (3d0, 5.1)/), order=(/2, 1/))

    call checkpoint(debug, string="dummy initilization with size:" ,vector=real(dims))
    ! check dummy initialization
    call Init(matrix, dims, error)

    ! deallocate elems before initializing again
    deallocate(matrix%elems)
    ! check complete initialization
    call Init(matrix, elems, error=error)


    write (*, *) "Original matrix:"
    do counter1 = 1, ubound(elems, 1)
        write (*, *) matrix%elems(counter1, :)
    end do
    write (*, '(/)', advance='no')

    tr = .Tr. matrix
    det = Determinant(elems)

    write (*, *) "Trace from new data structure: ", matrix%tr
    print *, "Trace from operator: ", tr
    write (*, *) "Determinant: ", det

    write (*, *) "Adjugate matrix from new data structure:"
    do counter1 = 1, ubound(elems, 1)
        write (*, *) matrix%adj_elems(counter1, :)
    end do
    write (*, '(/)', advance='no')

    adjoint = .Adj. matrix

    write (*, *) "Adjugate matrix from operator:"
    do counter1 = 1, ubound(elems, 1)
        write (*, *) adjoint%elems(counter1, :)
    end do
    write (*, '(/)', advance='no')

    ! test writing on file
    open(unit = unit, file="./prova.txt", action='write')
    call writecmatrix(matrix, unit)
    close(unit=unit)
end program