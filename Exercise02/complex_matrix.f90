module complex_matrix
    use checkpoint_mod

    implicit none

    logical, private :: debug = .True.  !< private logical variable to enable debugging within this module
    private          :: Trace, MatAdjoint, Init_elems, Init_size

    !< define new type complex matrix
    !!
    !! define also the main operations to handle complex matrix calculations.
    !! Two types of initialization are proposed: Init_size initializes only the 
    !! size of the matrix end set the elemnts to zero, Init_elems take care of
    !! initializing all the variables of the cmatrix type, including trace and 
    !! ajugate if possible.
    !! Base operations provaded are .Adj. to calculate the adjoint and .Tr. to
    !! calculate the trace
    !!
    type :: cmatrix
        double complex, allocatable :: elems(:, :)      !< elements of the matrix
        integer, dimension(2)       :: MN = (/0,0/)     !< size of matrix
        double complex              :: tr               !< trace
        double complex, allocatable :: adj_elems(:, :)  !< elements of the adjugate matrix
        
    end type cmatrix

    interface operator (.Tr.)
        module procedure Trace
    end interface

    interface operator (.Adj.)
        module procedure MatAdjoint
    end interface

    !> Two possible initilizations:
    !! - Init_size is a dummy initialization, initialize MN, allocate space for
    !!   elems and initialize it to zero
    !! - Init_elems initializes all the data structure, given elems, calculating 
    !!   also the trace and the adjugate matrix
    !!
    !! If a variable is dummy initialized you have to manually deallocate elems
    !! before using Init_elems
    !!
    interface Init
        module procedure Init_size
        module procedure Init_elems 
    end interface

    contains  
    
        subroutine Init_size(self, MN, error)
            type(cmatrix)                     :: self
            integer, dimension(2), intent(in) :: MN
            logical, intent(out), optional    :: error  !< optional variable to track down errors
            
            if (present(error)) then
                error = .False.
            end if

            call checkpoint(debug, string="Called Init_size")

            ! check if already initialized, otherwise raise an error
            if ((self%MN(1) /= 0) .and. (self%MN(2) /= 0)) then
                if (present(error)) then
                    error = .True.
                end if
                print *, "Dimensions already initialized!"
            else
                self%MN = MN
                allocate(self%elems(MN(1), MN(2)))
                self%elems = (0._8, 0.)
            end if
             
        end subroutine

        subroutine Init_elems(self, elems, error)
            type(cmatrix), intent(inout)   :: self
            double complex, intent(in)     :: elems(:,:)
            logical, intent(out), optional :: error         !< optional variable to track down errors
            integer, dimension(2)          :: MN

            call checkpoint(debug, string="Called Init_elems")

            MN = (/ubound(elems, 1), ubound(elems, 2)/)

            if (present(error)) then
                error = .False.
            end if

            ! check if already initialized, otherwise raise an error
            if (allocated(self%elems)) then
                print *, "Error: elements already initialized!"
                if (present(error)) then
                    error = .True.
                end if      
            else
                allocate(self%elems(MN(1), MN(2)))
                self%elems = elems
                self%MN = MN

                ! only if possible also calculate trace and adjugate
                if (MN(1) == MN(2)) then
                    allocate(self%adj_elems(MN(1), MN(2)))
                    self%tr = Trace(self)
                    self%adj_elems = transpose(Cofactor(elems, MN))
                end if
            end if
        end subroutine

        !> calculate trace of square matrix
        !!
        !! due to IO limits cannot call this function directly in
        !! a IO operation; declare a variable to store its result 
        !! and use the variable for IO operations instead
        !!
        function Trace(matrix) result(tr)
            type(cmatrix), intent(in) :: matrix
            double complex            :: tr
            integer                   :: ii

            ! check for non-square matrix
            if (matrix%MN(1) /= matrix%MN(2)) then
                write (*,*) "Error: The trace of a non square matrix doesn't exist!"
                return
            end if
            
            tr = (0._8, 0.)
            do ii = 1, matrix%MN(1)
                tr = tr + matrix%elems(ii, ii)
            end do
        end function Trace

        !> slice input matrix elements
        !!
        !! takes matrix elements and returns matrix element without
        !! selected lines and rows
        !!
        function slice(matrix_elems, i, j) result(sl_elems)
            double complex, intent(in)  :: matrix_elems(:, :)
            integer, intent(in)         :: i, j
            double complex, allocatable :: sl_elems(:, :)
            logical, allocatable        :: mask(:, :)
            integer :: mm, nn

            mm = ubound(matrix_elems, 1)
            nn = ubound(matrix_elems, 2)
            allocate(mask(mm, nn))
            allocate(sl_elems(mm-1, nn-1))

            mask = .true.
            mask(i, :) = .false.
            mask(:, j) = .false.

            sl_elems = reshape(pack(matrix_elems, mask), (/mm-1, nn-1/))
            deallocate(mask)
        end function

        !> calculate the determinant recursively
        !!
        !! calculate determinant slicing along the first column.
        !! due to IO limits cannot call this function directly in
        !! a IO operation; declare a variable to store its result 
        !! and use the variable for IO operations instead.
        !!
        recursive function Determinant(matrix_elems) result(det)
            double complex, intent(in) :: matrix_elems(:, :)
            double complex             :: det
            integer                    :: MN(2), jj

            MN = (/ubound(matrix_elems, 1), ubound(matrix_elems, 2)/)
            if (MN(1) /= MN(2)) then
                print *, "Error: The determinant of a non square matrix doesn't exist!"
                return
            else if (MN(1) == 1) then
                det = matrix_elems(1, 1)
            else
                det = (0._8, 0)
                do jj = 1, MN(1)
                    det = det + (-1._8)**(1+jj)*matrix_elems(1, jj)*Determinant(slice(matrix_elems, 1, jj))
                end do
            end if
        end function

        !> get the cofactor matrix, useful for the adjoint
        !!
        !! due to IO limits cannot call this function directly in
        !! a IO operation; declare a variable to store its result 
        !! and use the variable for IO operations instead
        !!
        function Cofactor(matrix_elems, MN) result(cofact)
            double complex, intent(in)  :: matrix_elems(:, :)
            integer, intent(in)         :: MN(2)
            integer                     :: ii, jj
            double complex, allocatable :: cofact(:, :)

            ! check for non-square matrix
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

        subroutine writecmatrix(matrix, unit)
            type(cmatrix), intent(in) :: matrix
            integer, intent(in)       :: unit
            integer                   :: ii, jj

            do ii = 1, matrix%MN(1)
                do jj = 1, matrix%MN(2)
                    write (unit, '(ES12.6, A, ES12.6, TR4)', advance='no') &
                        real(matrix%elems(ii, jj)), ' +i', aimag(matrix%elems(ii, jj))
                end do
                write (unit, *)
            end do
        end subroutine

end module