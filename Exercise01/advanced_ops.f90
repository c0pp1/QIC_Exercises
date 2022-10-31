module advanced_ops
    implicit none
    integer*8 ij
    
    
    contains
    function fact(n) result(result)
        integer*8 n
        integer*8 result
        result = 1
        if (n == 0 .or. n == 1) then
            
            return
        else
            do ij=2,n
                result = ij * result
            end do
            return 
        endif

    end function

    function arrmul(A, B) result(C)
        real*8 :: A(:)
        real*8 :: B(:)
        real*8 :: C
        integer(kind=kind(ubound(A, 1))) :: dimA
        integer(kind=kind(ubound(B, 1))) :: dimB
        integer(kind=kind(dimA)) :: ii

        dimA = ubound(A, 1)
        dimB = ubound(B, 1)
        C = 0._8

        if(dimA .ne. dimB) then
            print *, "Array dimensions are incompatible! Returning zero value"
        else
            do ii = 1, dimA
                C = C + A(ii)*B(ii)
            end do
        end if
    end function arrmul

    function matmul1(A, B) result(C)
        real*8 :: A(:, :)
        real*8 :: B(:, :)
        integer(kind=kind(ubound(A, 1))) :: dimA0 
        integer(kind=kind(ubound(B, 2))) :: dimB1
        integer(kind=kind(ubound(B, 1))) :: dimB0 
        integer(kind=kind(ubound(A, 2))) :: dimA1
        real*8 :: C(ubound(A, 1), ubound(B, 2))

        integer(kind=kind(dimA0)) :: counter1 = 1
        integer(kind=kind(dimB1)) :: counter2 = 1

        dimA0 = ubound(A, 1)
        dimB1 = ubound(B, 2)
        dimA1 = ubound(A, 2)
        dimB0 = ubound(B, 1)

        if(dimA1 .ne. dimB0) then
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

    function matmul2(A, B) result(C)
        real*8 :: A(:, :)
        real*8 :: B(:, :)
        integer(kind=kind(ubound(A, 1))) :: dimA0 
        integer(kind=kind(ubound(B, 2))) :: dimB1
        integer(kind=kind(ubound(B, 1))) :: dimB0 
        integer(kind=kind(ubound(A, 2))) :: dimA1
        real*8 :: C(ubound(A, 1), ubound(B, 2))

        integer(kind=kind(dimA0)) :: counter1 = 1
        integer(kind=kind(dimB1)) :: counter2 = 1

        dimA0 = ubound(A, 1)
        dimB1 = ubound(B, 2)
        dimA1 = ubound(A, 2)
        dimB0 = ubound(B, 1)

        if(dimA1 .ne. dimB0) then
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

    function matrixformatter(A) result(mformat)
        real*8 :: A(:, :)
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