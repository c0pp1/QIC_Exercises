module checkpoint_mod

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
    subroutine checkpoint (debug, string, scalar, vector)

        logical, intent(in)                    :: debug     !< enable/disable debugging
        character(len=*), intent(in), optional :: string    !< optional string
        real*4, intent(in), optional           :: scalar, vector(:)  !< optional real*4 variable
        character*10                           :: t         !< time
        if (debug) then
            call date_and_time(time=t)
            write (*, *) "checkpoint called at: ", t(1:2), ':', t(3:4), ':', t(5:)
            if (present(string)) then
                write (*, *) string
            end if
            if (present(scalar)) then
                print *, variable
            end if  
            if (present(vector)) then
                print *, vector
            end if            
        end if

    end subroutine
end module