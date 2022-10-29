

program test
    use advanced_ops
    implicit none
    integer*8 n

    print *, "Insert an integer (at most 20, sorry):"
    read(*,*) n

    print *, "The result is:", fact(n)
stop
end program test