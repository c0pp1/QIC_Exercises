program main

    external zheev

    integer*4 :: m, NN, ii, info, outunit = 100, lwork = -1, seed(4) = (/1,1,1,1/)
    double complex, allocatable :: A(:, :), work(:), D(:)
    double precision, allocatable :: eigenvals(:), spacings(:), rwork(:)
    character(len=200) :: ofdir, ofname
    character*1 :: creturn = achar(13), jobz = 'N', uplo = 'U'
    character(len=32) :: arg, format
    logical           :: skip_next = .false., debug = .false.

    integer :: MODE = 6
    double precision :: COND = 1.1d0
    complex*16 :: DMAX
    character :: RSIGN
    character :: GRADE = 'N'
    complex*16, allocatable :: DL(:)
    character :: PIVTNG = 'N'
    integer, allocatable :: IPIVOT(:)
    double precision :: SPARSE = 0.
    double precision :: ANORM = -1.
    character :: PACK = 'N'
    integer, allocatable :: IWORK(:)

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
                
            case ('-n', '--iter')
                call get_command_argument(ii+1, arg)
                read (arg, *) NN
                skip_next = .true.
                
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
    allocate( A(m, m), work(1), eigenvals(m), spacings(m-1), D(m), DL(m), IWORK(m), rwork(max(1, 3*m-2)))

    !> open output file
    write (format, '(A, I2, A)') "(A, I", floor(log10(real(m)))+1, ", A)'"
    write (ofname, format) trim(ofdir), m, ".txt"
    open(unit = outunit, file=trim(ofname), action='write')

    !> do calculations
    do ii=1, NN
        ! generate random hermitian matrix
        call zlatmr(m, m, 'U', seed, 'H', D, MODE, COND, DMAX, RSIGN, GRADE, D, MODE, COND, D, &
                    MODE, COND, PIVTNG, IPIVOT, m, m, SPARSE, ANORM, PACK, A, m, IWORK, info)
        if(info /= 0) then
            print *, "ERROR: failed to generate #", ii, " random hermitian matrix, skipping..."
            cycle
        end if

        if(debug) print *, A

        !> diagonalize and get eigenvalues
        !! the first call determines the optimal value for lwork 
        !! and stores it in work(1). Then we call again the subroutine
        !! with lwork optimal value
        !!
        call zheev(jobz, uplo, m, A, m, eigenvals, work, lwork, rwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call zheev(jobz, uplo, m, A, m, eigenvals, work, lwork, rwork, info)
        if(info /= 0) then
            print *, "ERROR: failed to get eigenvalues for matrix ", ii, " skipping..."
            cycle
        end if

        spacings = (eigenvals(2:) - eigenvals(:-1))/ merge((sum(eigenvals)/max(1,m)), 1.d0, (sum(eigenvals)/max(1,m))/=0.)
        write(outunit, *) spacings

        write (*, '(A, A, F6.2, A)', advance='no') creturn, "Spacings calculations, ", ii*100./NN, "% done."
    end do
    write(*, *)

    contains

        subroutine print_help()
            print *, "Usage: ./exercise03a -m [matrix size] -n [iterations] -o [output dir/first part of filename]"
        end subroutine

end program