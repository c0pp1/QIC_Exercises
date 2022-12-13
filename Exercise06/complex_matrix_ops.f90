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
        procedure write_cmplx_1
        procedure write_cmplx_2
    end interface


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
            double complex :: matrix(:,:), amatrix(ubound(matrix, 1),ubound(matrix, 2))
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

end module