module complex_matrix_ops
    implicit none
    character(1), parameter :: tab = char(9)

    interface operator(.Tr.)
        module procedure Trace
    end interface
    interface operator(.Adj.)
        module procedure Adjoint
    end interface

    interface write_cmplx
        module procedure write_cmplx_1
        module procedure write_cmplx_2
    end interface write_cmplx

    interface tensor_prod
        module procedure tensor_prod_int
        module procedure tensor_prod_dc
    end interface tensor_prod

    contains

        function Trace(matrix) result(tr)
            double complex, intent(in) :: matrix(:,:)
            double complex :: tr
            integer :: m, ii

            m = ubound(matrix, 1)
            tr = (0._8, 0._8)
            do ii=1, m
                tr = tr + matrix(ii, ii)
            end do
        end function

        function Adjoint(matrix) result(amatrix)
            double complex :: matrix(:,:), amatrix(ubound(matrix, 2),ubound(matrix, 1))
            intent(in) :: matrix

            amatrix = conjg(transpose(matrix))
        end function


        function matrixformatter(A) result(mformat)
            double complex :: A(:, :)
            integer(kind=kind(ubound(A, 1))) :: dimA0
            integer(kind=kind(ubound(A, 2))) :: dimA1
            integer(kind=kind(dimA0)) :: ii
            character(len=:), allocatable :: mformat
    
            dimA0 = ubound(A, 1)
            dimA1 = ubound(A, 2)
    
            mformat = ""
            do ii = 1, dimA0
                mformat = mformat // repeat("TR2ES9.2','ES9.2,", dimA1)
                mformat = mformat // "/"
            end do
            mformat = trim(mformat(:len(mformat)-1 )//' ')
        end function matrixformatter

        subroutine write_cmplx_1(outunit, array)
            integer :: outunit, ii
            double complex, allocatable :: array(:)
            intent(inout) :: outunit, array

            do ii=1,ubound(array, 1)
                write(outunit, '(ES12.4E3, A, ES12.4E3, A)', advance='no') real(array(ii)), tab, aimag(array(ii)), tab
            end do

        end subroutine

        subroutine write_cmplx_2(outunit, array)
            integer :: outunit, ii, jj
            double complex :: array(:, :)
            intent(inout) :: outunit, array

            do ii=1, ubound(array, 1)
                do jj=1, ubound(array, 2)
                    write(outunit, '(ES12.4E3, A, ES12.4E3, A)', advance='no') real(array(ii, jj)), tab, aimag(array(ii, jj)), tab
                end do
                write(outunit, *)
            end do

        end subroutine

        !> tensor product of bidimensional arrays
        !!
        !! given two bidimensional arrays calculate the tensor product.
        !! The result is a bidiminesional array of size:
        !! (size(arr1, 1)*size(arr2, 1), size(arr1, 2)*size(arr2, 2))
        !!
        !! @param[in]   arr1    the first array
        !! @param[in]   arr2    the second array
        !! @result      arr3    the resulting tensor product
        function tensor_prod_int(arr1, arr2) result(arr3)
            integer :: arr1(:,:), arr2(:,:), arr3(:,:), ii, jj, m1, n1, m2, n2
            allocatable :: arr3
            intent(in) :: arr1, arr2

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

        function tensor_prod_dc(arr_dc1, arr_dc2) result(arr3)
            double complex :: arr_dc1(:,:), arr_dc2(:,:), arr3(:,:)
            integer :: ii, jj, m1, n1, m2, n2
            allocatable :: arr3
            intent(in) :: arr_dc1, arr_dc2

            m1 = size(arr_dc1, 1)
            n1 = size(arr_dc1, 2)
            m2 = size(arr_dc2, 1)
            n2 = size(arr_dc2, 2)


            allocate(arr3(m1*m2, n1*n2))

            do jj=1, n1
                do ii=1, m1
                    arr3((ii-1)*m2+1:ii*m2, (jj-1)*n2+1:jj*n2) = arr_dc1(ii, jj)*arr_dc2
                end do
            end do
            
        end function

end module