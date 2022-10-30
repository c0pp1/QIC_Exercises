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

    function matmul1(A, B) result(C)
        real*8 :: A(:, :)
        real*8 :: B(:, :)
        real*8 :: C(ubound(A, 1), ubound(B, 2))

        integer(kind=kind(ubound(A, 1))) :: counter1 = 1
        integer(kind=kind(ubound(B, 2))) :: counter2 = 1

        if(ubound(A, 2) .ne. ubound(B, 1)) then

        end if


    end function matmul1


    function matmul2(A, B) result(C)
        real*8 :: A(:, :)
        real*8 :: B(:, :)
        real*8 :: C(ubound(A, 1), ubound(B, 2))


    end function matmul2
end module advanced_ops