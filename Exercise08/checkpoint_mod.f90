module checkpoint_mod
    use complex_matrix_ops
    implicit none
    integer :: DB = 5
    character(1), parameter :: ESC = achar(27)
    character(*), parameter :: BRED = ESC // '[31;1m', &
                                RED = ESC // '[31m', &
                              BLUE  = ESC // '[34m', &
                             YELLOW = ESC // '[33m', &
                              RESET = ESC // '[0m'

    contains
    !> used as checkpoint for debugging purposes
    !!
    !! subroutine to be called as a checkpoint: prints also optional string
    !! and variable to check if the execution is running as expected 
    !!
    !! @param[in]   debug   logical parameter to enable/disable checkpoints
    !! @param[in]   string  optional string to be printed in a new line
    !! @param[in]   variable    optional real*4 variable to be printed in a new line
    !! @param       t   variable to store the system time at which the subroutine was called
    subroutine checkpoint (debug_level, string, scalar, vector, cmatrix, advance)

        integer, intent(in)                    :: debug_level       !< enable/disable debugging
        character(len=*), intent(in), optional :: string, advance   !< optional strings
        real*8, optional                       :: scalar, vector(:) !< optional real*8 variable
        double complex, optional               :: cmatrix(:,:)
        intent(in)                             :: scalar, vector, cmatrix
        character*10                           :: t         !< time
        character(255)                         :: toprint
        if (debug_level <= DB) then
            call date_and_time(time=t)
            toprint = "[" // t(1:2) // ':' // t(3:4) // ':' // t(5:) // "]"
            select case(debug_level)
            case(0)
                toprint = trim(toprint) // '[' // BRED // 'FATAL' // RESET // ']'
            case(1)
                toprint = trim(toprint) // '[' // RED // 'ERROR' // RESET // ']'
            case(2)
                toprint = trim(toprint) // '[' // YELLOW // 'WARN' // RESET // ' ]'
            case default
                toprint = trim(toprint) // '[' // BLUE // 'INFO' // RESET // ' ]'
            end select
            if (present(string)) then
                toprint = trim(toprint) // ' ' // string
            end if
            if(present(advance)) then
                write(*, '(A)', advance=advance) trim(toprint)
            else
                write(*, *) trim(toprint)
            end if
            if (present(scalar)) then
                print *, scalar
            end if  
            if (present(vector)) then
                print *, vector
            end if    
            if (present(cmatrix)) then
                write(*, '(' // matrixformatter(cmatrix) // ')') transpose(cmatrix)
            end if         
        end if

    end subroutine
end module