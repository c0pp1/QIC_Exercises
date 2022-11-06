module complex_matrix
    use checkpoint_mod

    implicit none

    logical, private :: debug = .True.

    !< define new type complex matrix
    !!
    !! define also the main operations to handle complex matrix calculations.
    !! Two types of initialization are proposed: Init_size initializes only the 
    !! size of the matrix end set the elemnts to zero, Init_elems take care of
    !! initializing all the variables of the cmatrix type, including trace and 
    !! ajugate if possible.
    !!
    type :: cmatrix
        double complex, allocatable :: elems(:, :)
        integer, dimension(2)       :: MN = (/0,0/)
        double complex              :: tr
        double complex, allocatable :: adj_elems(:, :)
        
    end type cmatrix

    interface operator (.Tr.)
        module procedure Trace
    end interface

    interface operator (.Adj.)
        module procedure MatAdjoint
    end interface

    interface Init
        module procedure Init_size
        module procedure Init_elems 
    end interface

    contains  
    
        subroutine Init_size(self, MN, error)
            type(cmatrix)                     :: self
            integer, dimension(2), intent(in) :: MN
            logical, intent(out), optional    :: error
            
            if (present(error)) then
                error = .False.
            end if

            call checkpoint(debug, string="Called Init_size")
            if ((self%MN(1) /= 0) .and. (self%MN(2) /= 0)) then
                if (present(error)) then
                    error = .True.
                end if
                print *, "Dimensions already initialized!"
            else
                self%MN = MN
                allocate(self%elems(MN(1), MN(2)), self%adj_elems(MN(1), MN(2)))
                self%elems = (0._8, 0.)
            end if
             
        end subroutine

        subroutine Init_elems(self, elems, error)
            type(cmatrix), intent(inout)   :: self
            double complex, intent(in)     :: elems(:,:)
            logical, intent(out), optional :: error
            integer, dimension(2)          :: MN

            call checkpoint(debug, string="Called Init_elems")

            MN = (/ubound(elems, 1), ubound(elems, 2)/)

            if (present(error)) then
                error = .False.
            end if

            if (allocated(self%elems)) then
                print *, "Error: elements already initialized!"
                if (present(error)) then
                    error = .True.
                end if      
            else
                allocate(self%elems(MN(1), MN(2)))
                self%elems = elems
                self%MN = MN
                if (MN(1) == MN(2)) then
                    allocate(self%adj_elems(MN(1), MN(2)))
                    self%tr = Trace(self)
                    self%adj_elems = transpose(Cofactor(elems, MN))
                end if
            end if
        end subroutine

        function Trace(matrix) result(tr)
            type(cmatrix), intent(in) :: matrix
            double complex            :: tr
            integer                   :: ii

            if (matrix%MN(1) /= matrix%MN(2)) then
                print *, "Error: The trace of a non square matrix doesn't exist!"
                return
            else 
                tr = (0._8, 0.)
                do ii = 1, matrix%MN(1)
                    tr = tr + matrix%elems(ii, ii)
                end do
            end if
        end function Trace

        function slice(matrix_elems, i, j) result(sl_elems)
            double complex, intent(in)  :: matrix_elems(:, :)
            integer, intent(in)         :: i, j
            double complex, allocatable :: sl_elems(:, :)
            logical, allocatable        :: mask(:, :)
            integer :: m, n

            m = ubound(matrix_elems, 1)
            n = ubound(matrix_elems, 2)
            allocate(mask(m, n))
            allocate(sl_elems(m-1, n-1))

            mask = .true.
            mask(i, :) = .false.
            mask(:, j) = .false.

            sl_elems = reshape(pack(matrix_elems, mask), (/m-1, n-1/))
            deallocate(mask)
        end function

        recursive function Determinant(matrix_elems) result(det)
            double complex, intent(in) :: matrix_elems(:, :)
            double complex             :: det
            integer                    :: m, n, jj

            m = ubound(matrix_elems, 1)
            n = ubound(matrix_elems, 2)
            if (m /= n) then
                print *, "Error: The determinant of a non square matrix doesn't exist!"
                return
            end if

            if (n == 1) then
                det = matrix_elems(1, 1)
            else
                det = (0._8, 0)
                do jj = 1, n
                    det = det + (-1._8)**(1+jj)*matrix_elems(1, jj)*Determinant(slice(matrix_elems, 1, jj))
                end do
            end if
        end function

        function Cofactor(matrix_elems, MN) result(cofact)
            double complex, intent(in)  :: matrix_elems(:, :)
            integer, intent(in)         :: MN(2)
            integer                     :: ii, jj
            double complex, allocatable :: cofact(:, :)

            if (MN(1) /= MN(2)) then
                print *, "Error: The cofactor of a non square matrix doesn't exist!"
                return
            end if

            allocate(cofact(MN(1), MN(2)))

            do ii = 1, MN(1)
                do jj = 1, MN(2)
                    cofact(ii, jj) = (-1._8)**(ii+jj)*Determinant(slice(matrix_elems, ii, jj))
                end do
            end do

        end function

        function MatAdjoint(matrix) result(adjoint)
            type(cmatrix), intent(in):: matrix
            type(cmatrix) :: adjoint

            call Init(adjoint, matrix%MN)
            adjoint%elems = transpose(Cofactor(matrix%elems, matrix%MN))
        end function MatAdjoint

end module