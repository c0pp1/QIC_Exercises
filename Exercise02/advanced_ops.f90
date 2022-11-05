module advanced_ops
    implicit none
    integer*8 ij    !< integer used as counter
      
    contains
    !> calculate the factorial of a positive integer
    !!
    !! @param[in]   n       must be a positive integer lesser than 20 due to max_int limitation
    !! @return      result  result of the factorial, initialized to 1
    function fact(n) result(result)
        integer*8, intent(in)  :: n
        integer*8              :: result
        result = 1

        
        if ( n > 20 ) then  !< check if n <= 20 as 21! ~= 5e19 > max_int_8 = 2^(8*8) ~= 2e19
            print *, 'The number is too big for the result to be representet in 8 bytes integer, returning 1'
        else if (n == 0 .or. n == 1) then   !< in this cases the factorial is 1, so immediately return
            return
        else    !< in other cases loop over a range from 2 to n
            do ij=2,n
                result = ij * result
            end do
            return 
        endif

    end function

    !> array multiplication
    !!
    !! this function takes two equal size arrays and performs their scalar multiplication
    !! if the lengths are different returns zero, printing an error
    !!
    !! @param[in]   A, B    the two arrays
    !! @param   dimA    the length of A array
    !! @param   dimB    the length of B array
    !! @param   ii      integer for loop over array length
    !! @result  C       resulting scalar
    function arrmul(A, B) result(C)
        real*8, intent(in) :: A(:)
        real*8, intent(in) :: B(:)
        real*8 :: C
        integer(kind=kind(ubound(A, 1))) :: dimA
        integer(kind=kind(ubound(B, 1))) :: dimB
        integer(kind=kind(dimA)) :: ii

        dimA = ubound(A, 1)
        dimB = ubound(B, 1)
        C = 0._8    !< initialize the result to real*8 zero

        if(dimA /= dimB) then   !< check for equal length arrays
            print *, "Array dimensions are incompatible! Returning zero value"
        else
            do ii = 1, dimA
                C = C + A(ii)*B(ii)
            end do
        end if
    end function arrmul

    !> matrix multiplication
    !!
    !! performs the matrix multiplication between A and B matrices, checking their shapes:
    !! if the second dimension of A and the first one of B are different returns a zero
    !! matrix with shape (ubound(A, 1), ubound(B, 2)).
    !! calculations are perfomed first on a row of matrix C (looping along matrix C column index)
    !! and then changing row (looping along matrix C row index) and looping again over column index
    !!
    !! @param[in]   A   first matrix
    !! @param[in]   B   second matrix
    !! @return      C   resulting matrix with shape (ubound(A, 1), ubound(B, 2))
    function matmul1(A, B) result(C)
        real*8, intent(in) :: A(:, :)
        real*8, intent(in) :: B(:, :)
        integer(kind=kind(ubound(A, 1))) :: dimA0   !< matrix A first dimension
        integer(kind=kind(ubound(B, 2))) :: dimB1   !< matrix B second dimension
        integer(kind=kind(ubound(B, 1))) :: dimB0   !< matrix B first dimension
        integer(kind=kind(ubound(A, 2))) :: dimA1   !< matrix A second dimension 
        real*8 :: C(ubound(A, 1), ubound(B, 2))

        integer(kind=kind(dimA0)) :: counter1 = 1   !< counter to loop along rows
        integer(kind=kind(dimB1)) :: counter2 = 1   !< counter to loop along columns

        dimA0 = ubound(A, 1)
        dimB1 = ubound(B, 2)
        dimA1 = ubound(A, 2)
        dimB0 = ubound(B, 1)

        if(dimA1 /= dimB0) then !< check if the second dimension of A and the first one of B are different, if so return zero
            print *, "Matrix dimensions are incompatible! Returning ", dimA0, "x", ubound(B, 2), &
                     " matrix with zero value"
            C = 0._8
        else
            do counter1 = 1, dimA0
                do counter2 = 1, dimB1
                    C(counter1, counter2) = arrmul(A(counter1, :), B(:, counter2))
                end do
            end do
        end if


    end function matmul1

    !> matrix multiplication
    !!
    !! performs the matrix multiplication between A and B matrices, checking their shapes:
    !! if the second dimension of A and the first one of B are different returns a zero
    !! matrix with shape (ubound(A, 1), ubound(B, 2)).
    !! calculations are perfomed first on a column of matrix C (looping along matrix C row index)
    !! and then changing column (looping along matrix C column index) and looping again over row index
    !!
    !! @param[in]   A   first matrix
    !! @param[in]   B   second matrix
    !! @return      C   resulting matrix with shape (ubound(A, 1), ubound(B, 2))
    function matmul2(A, B) result(C)
        real*8 :: A(:, :)
        real*8 :: B(:, :)
        integer(kind=kind(ubound(A, 1))) :: dimA0   !< matrix A first dimension
        integer(kind=kind(ubound(B, 2))) :: dimB1   !< matrix B second dimension
        integer(kind=kind(ubound(B, 1))) :: dimB0   !< matrix B first dimension
        integer(kind=kind(ubound(A, 2))) :: dimA1   !< matrix A second dimension
        real*8 :: C(ubound(A, 1), ubound(B, 2))

        integer(kind=kind(dimA0)) :: counter1 = 1   !< counter to loop along rows
        integer(kind=kind(dimB1)) :: counter2 = 1   !< counter to loop along columns

        dimA0 = ubound(A, 1)
        dimB1 = ubound(B, 2)
        dimA1 = ubound(A, 2)
        dimB0 = ubound(B, 1)

        if(dimA1 /= dimB0) then !< check if the second dimension of A and the first one of B are different, if so return zero
            print *, "Matrix dimensions are incompatible! Returning ", dimA0, "x", ubound(B, 2), &
                     " matrix with zero value"
            C = 0._8
        else
            do counter2 = 1, dimB1
                do counter1 = 1, dimA0
                    C(counter1, counter2) = arrmul(A(counter1, :), B(:, counter2))
                end do
            end do
        end if

    end function matmul2

    !> format string to print matrices
    !!
    !! takes in input the matrix you want to print and gives the format string
    !! to have a nice output. Numbers are printed in scientific representation
    !!
    !! @param[in]   A   the matrix to be printed, shape length 2
    !! @return      mformat the string with the format
    function matrixformatter(A) result(mformat)
        real*8, intent(in) :: A(:, :)
        integer(kind=kind(ubound(A, 1))) :: dimA0
        integer(kind=kind(ubound(A, 2))) :: dimA1
        integer(kind=kind(dimA0)) :: ii
        character(len=:), allocatable :: mformat

        dimA0 = ubound(A, 1)
        dimA1 = ubound(A, 2)

        mformat = ""
        do ii = 1, dimA0
            mformat = mformat // repeat("TR2ES8.2,", dimA1)
            mformat = mformat // "/"
        end do

    end function matrixformatter
end module advanced_ops