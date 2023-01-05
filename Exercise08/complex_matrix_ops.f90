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

    interface operator(.tprod.)
        module procedure tensor_prod_int
        module procedure tensor_prod_dc
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

        subroutine get_density_matrix(state, density_matrix)
            double complex, allocatable :: density_matrix(:,:)
            double complex :: state(:)
            integer :: info=0, m
            intent(in) :: state
            intent(out) :: density_matrix

            m = ubound(state, 1)
            if(.not. allocated(density_matrix)) allocate(density_matrix(m, m), stat=info)
            
            density_matrix = spread(state, 2, m)* .Adj. spread(state, 2, m)

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
        subroutine partial_trace(density_matrix, systems, dims, partial_dm)
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
                    
                    tr_mat = tmp .tprod. reshape(identity, (/dims(ii), dims(ii)/))
                    deallocate(identity)

                else
                    mask(findloc(mask, .true.)) = .false.
                    ! allocate space for result
                    allocate(tr_mat(ubound(tmp, 1), ubound(tmp, 2)*dims(ii)))
                    allocate(tmp_sys(1, dims(ii)))
                    tmp_sys = 0
                    tmp_sys(1, mod((el-1)/product(pack(dims, mask)), dims(ii))+1) = 1

                    tr_mat = tmp .tprod. tmp_sys
                    deallocate(tmp_sys)
                end if

                deallocate(tmp)
            end do
        end function


end module