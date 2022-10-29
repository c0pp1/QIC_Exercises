module advanced_ops
    implicit none
    integer*8 ii
    
    
    contains
    function fact(n) result(result)
        integer*8 n
        integer*8 result
        result = 1
        if (n == 0 .or. n == 1) then
            
            return
        else
            do ii=2,n
                result = ii * result
            end do
            return 
        endif

    end function
end module advanced_ops