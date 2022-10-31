program precision
    !#################################################################################################
    !exercise01b
    !compilation command:
    !gfortran -fno-range-check exercise01b.f90
    !the option "-fno-range-check" forces the compiler to ignore overflows in variable initializations
    !#################################################################################################
    
    implicit none

    !declare initial values, with different types
    !a)
    integer*2, parameter :: int2 = 2000000_2
    integer*4, parameter :: int4 = 2000000
    real*4, parameter :: PI4 = 4*ATAN(1.)*10.**32
    real*8, parameter :: PI8 = 4*DATAN(1.d0)*10._8**32
    !b)
    character(len=36), parameter :: table1_format = '(T15A, T20A, /A, I8, I8, /A, I8, I8)'
    character(len=81), parameter :: table2_format = '(T18A, TR10A, TR2A, /A, ES21.13, ES21.13, ES21.13, &
                                                    &/A, ES21.13, ES21.13, ES21.13)'

    !print quite nice output, showing differences between results and initial values
    write (*,'(A)')  "INTEGER precision"
    write (*, table1_format) "SUM", "START VALUE", "INTEGER*2: ", int2 + 1, int2, "INTEGER*4: ", int4 + 1, int4

    write (*,'(/A)') "REAL precision"
    write (*, table2_format) "VALUE", "FIRST ADDEND VALUE", "SECOND ADDEND VALUE", "REAL*4: ", PI4 + sqrt(2.)*10.**21, PI4, &
                             sqrt(2.)*10.**21, "REAL*8: ", PI8 + dsqrt(2._8)*10._8**21, PI8, dsqrt(2._8)*10._8**21
    
    stop
end program precision