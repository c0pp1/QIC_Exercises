program main

    use lapack_wrapper

    integer*4 :: m, NN, ii, info, outunit = 100
    double complex, allocatable :: A(:, :)
    double precision, allocatable :: rA(:, :)
    double precision, allocatable :: eigenvals(:), spacings(:)
    character(len=200) :: ofdir, ofname
    character*1 :: creturn = achar(13)
    character(len=32) :: arg, format
    logical           :: skip_next = .false., neglect_1eig = .false., debug = .false., realA = .false.


    !> read command line arguments
    !!
    do ii = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-m', '--dimension')
                call get_command_argument(ii+1, arg)
                read (arg, *) m
                skip_next = .true.
                
            case ('-i', '--iter')
                call get_command_argument(ii+1, arg)
                read (arg, *) NN
                skip_next = .true.

            case ('-n', '--neglect1')
                neglect_1eig = .true.

            case ('-r', '--real')
                realA = .true.
                
            case ('-o', '--odirname')
                call get_command_argument(ii+1, arg)
                read (arg, '(A)') ofdir
                skip_next = .true.

            case ('-h', '--help')
                call print_help()
                skip_next = .false.
                stop

            case default
                print '(2a, /)', 'unrecognised command-line option: ', arg
                call print_help()
                skip_next = .false.
                stop
        end select
    end do

    !> allocate matrix and related arrays
    allocate(eigenvals(m), spacings(merge(m-2, m-1, neglect_1eig)))
    if(realA) then
        allocate(rA(m, m))
    else 
        allocate(A(m, m))
    end if

    !> open output file
    write (format, '(A, I2, A)') "(A, I", floor(log10(real(m)))+1, ", A)'"
    write (ofname, format) trim(ofdir), m, trim(merge("r.txt", ".txt ", realA))
    open(unit = outunit, file=trim(ofname), action='write')

    !> do calculations
    do ii=1, NN
        ! generate random hermitian matrix
        if(realA) then
            call latmr_wrap(A=rA, SYM='S', INFO=info)
        else
            call latmr_wrap(A=A, SYM='H', INFO=info)
        end if
        if(debug) print *, "done"
        if(info /= 0) then
            print *, "ERROR: failed to generate #", ii, " random hermitian matrix, skipping..."
            cycle
        end if

        if(debug) print *, A

        ! diagonalize and get eigenvalues
        if(realA) then
            call ev_wrap( A=rA, W=eigenvals, INFO=info)
        else
            call ev_wrap( A=A, W=eigenvals, INFO=info)
        end if
        if(info /= 0) then
            print *, "ERROR: failed to get eigenvalues for matrix ", ii, " skipping..."
            cycle
        end if

        if(neglect_1eig) then
            spacings = (eigenvals(3:) - eigenvals(2:m-1))
            spacings = spacings / (sum(spacings)/max(1,m-2))
        else
            spacings = (eigenvals(2:) - eigenvals(:-1))
            spacings = spacings / (sum(spacings)/max(1,m-1))
        end if
        write(outunit, *) spacings

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Spacings calculations, ", ii*100./NN, "% done."
    end do
    write(*, *)

    contains

        subroutine print_help()
            print *, "Usage: ./exercise03b -m [matrix size] -i [iterations]" // &
                     "-n [optional to neglect first eigenval] -o [output dir/first part of filename]"
        end subroutine

end program