program main 
    use complex_matrix
    use checkpoint_mod
    use advanced_ops

    type(cmatrix) :: matrix
    integer, dimension(2) :: dims = (/1,2/)
    double complex :: elems
    logical :: debug = .True.
    logical :: error

    call checkpoint(debug, vector=real(dims))
    call Init(matrix, dims, error)

    print *, matrix%MN, error, matrix%elems
end program