module lapack_wrapper

    implicit none
    !> declare default values for arguments of lapack routines
    integer, parameter :: D_MODE = 6, D_SEED(4) = (/1,1,1,1/)
    character, parameter :: D_DIST = 'S', D_RSIGN = 'F', D_GRADE = 'N', D_PIVTNG = 'N', D_PACK = 'N'
    double precision, parameter :: D_COND = 1.1d0, D_SPARSE = 0., D_ANORM = -1.

    integer, parameter :: D_LWORK = -1
    character, parameter :: D_JOBZ = 'N', D_UPLO = 'U'

    logical, private, parameter :: debug = .false.

    interface latmr_wrap
        module procedure zlatmr_wrap
        module procedure dlatmr_wrap
    end interface

    interface ev_wrap
            module procedure zheev_wrap
            module procedure dsyev_wrap
    end interface

    contains 

        subroutine zlatmr_wrap(M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, &
                               CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO)
            
            complex*16, dimension(:, :) :: A
            integer, optional :: M, N
            character, optional :: DIST
            integer, dimension(4), optional :: ISEED
            character :: SYM
            complex*16, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: D
            integer, optional :: MODE
            double precision, optional :: COND
            complex*16, optional :: DMAX
            character, optional :: RSIGN
            character, optional :: GRADE
            complex*16, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: DL
            integer, optional :: MODEL
            double precision, optional :: CONDL
            complex*16, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: DR
            integer, optional :: MODER
            double precision, optional :: CONDR
            character, optional :: PIVTNG
            integer, optional, dimension(:) :: IPIVOT
            integer, optional :: KL
            integer, optional :: KU
            double precision, optional :: SPARSE
            double precision, optional :: ANORM
            character, optional :: PACK
            integer, optional :: LDA
            integer, optional, dimension(:) :: IWORK
            integer :: INFO 

            !> declare non optional variables to be passed to zlatmr in place of optional variables
            integer :: noM, noN, noISEED(4), noMODE, noMODEL, noMODER, noKU, noKL, noLDA
            character :: noDIST, noRSIGN, noGRADE, noPIVTNG, noPACK
            complex*16, dimension(min(ubound(A, 1), ubound(A, 2))) :: noD, noDL, noDR
            double precision :: noCOND, noCONDL, noCONDR, noSPARSE, noANORM
            complex*16 :: noDMAX
            integer, allocatable :: noIPIVOT(:), noIWORK(:)

            if(debug) print *, 'zlatmr was called'

            noM = merge(M, size(A, 1), present(M))
            noN = merge(N, size(A, 2), present(N))
            noDIST = merge(DIST, D_DIST, present(DIST))
            noISEED = merge(ISEED, D_SEED, present(ISEED))
            noMODE = merge(MODE, D_MODE, present(MODE))
            noCOND = merge(COND, D_COND, present(COND))
            noDMAX = merge(DMAX, (0.D0, 0.), present(DMAX))
            noRSIGN = merge(RSIGN, D_RSIGN, present(RSIGN))
            noGRADE = merge(GRADE, D_GRADE, present(GRADE))
            if(present(D)) then
                noD = D 
                noDL = merge(DL, noD, present(DL))
                noDR = merge(DR, noD, present(DR))
            end if
            noMODEL = merge(MODEL, D_MODE, present(MODEL))
            noMODER = merge(MODER, D_MODE, present(MODER))
            noCONDL = merge(CONDL, D_COND, present(CONDL))
            noCONDR = merge(CONDR, D_COND, present(CONDR))
            noPIVTNG = merge(PIVTNG, D_PIVTNG, present(PIVTNG))
            noKL = merge(KL, noM, present(KL))
            noKU = merge(KU, noN, present(KU))
            noSPARSE = merge(SPARSE, D_SPARSE, present(SPARSE))
            noANORM = merge(ANORM, D_ANORM, present(ANORM))
            noPACK = merge(PACK, D_PACK, present(PACK))
            noLDA = merge(LDA, noM, present(LDA))
            if(present(IWORK)) then
                allocate(noIWORK(size(IWORK)))
                noIWORK = IWORK
            end if
            if(present(IPIVOT)) then
                allocate(noIPIVOT(size(IPIVOT)))
                noIPIVOT = IPIVOT
            end if

            call zlatmr(noM, noN, noDIST, noISEED, SYM, noD, noMODE, noCOND, noDMAX, noRSIGN, noGRADE, noDL, noMODEL, noCONDL, &
                        noDR, noMODER, noCONDR, noPIVTNG, noIPIVOT, noKL, noKU, noSPARSE, noANORM, noPACK, A, noLDA, noIWORK, INFO)

        end subroutine

        subroutine dlatmr_wrap(M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, &
                               CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO)

            double precision, dimension(:, :) :: A
            integer, optional :: M, N
            character, optional :: DIST
            integer, dimension(4), optional :: ISEED
            character :: SYM
            double precision, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: D
            integer, optional :: MODE
            double precision, optional :: COND
            double precision, optional :: DMAX
            character, optional :: RSIGN
            character, optional :: GRADE
            double precision, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: DL
            integer, optional :: MODEL
            double precision, optional :: CONDL
            double precision, optional, dimension(min(ubound(A, 1), ubound(A, 2))) :: DR
            integer, optional :: MODER
            double precision, optional :: CONDR
            character, optional :: PIVTNG
            integer, optional, dimension(:) :: IPIVOT
            integer, optional :: KL
            integer, optional :: KU
            double precision, optional :: SPARSE
            double precision, optional :: ANORM
            character, optional :: PACK
            integer, optional :: LDA
            integer, optional, dimension(:) :: IWORK
            integer :: INFO 

            !> declare non optional variables to be passed to zlatmr in place of optional variables
            integer :: noM, noN, noISEED(4), noMODE, noMODEL, noMODER, noKU, noKL, noLDA
            character :: noDIST, noRSIGN, noGRADE, noPIVTNG, noPACK
            double precision, dimension(min(ubound(A, 1), ubound(A, 2))) :: noD, noDL, noDR
            double precision :: noCOND, noCONDL, noCONDR, noSPARSE, noANORM, noDMAX
            integer, allocatable :: noIPIVOT(:), noIWORK(:)

            if(debug) print *, 'dlatmr was called'

            noM = merge(M, size(A, 1), present(M))
            noN = merge(N, size(A, 2), present(N))
            noDIST = merge(DIST, D_DIST, present(DIST))
            noISEED = merge(ISEED, D_SEED, present(ISEED))
            noMODE = merge(MODE, D_MODE, present(MODE))
            noCOND = merge(COND, D_COND, present(COND))
            noDMAX = merge(DMAX, 0.D0, present(DMAX))
            noRSIGN = merge(RSIGN, D_RSIGN, present(RSIGN))
            noGRADE = merge(GRADE, D_GRADE, present(GRADE))
            if(present(D)) then
                noD = D 
                noDL = merge(DL, noD, present(DL))
                noDR = merge(DR, noD, present(DR))
            end if
            noMODEL = merge(MODEL, D_MODE, present(MODEL))
            noMODER = merge(MODER, D_MODE, present(MODER))
            noCONDL = merge(CONDL, D_COND, present(CONDL))
            noCONDR = merge(CONDR, D_COND, present(CONDR))
            noPIVTNG = merge(PIVTNG, D_PIVTNG, present(PIVTNG))
            noKL = merge(KL, noM, present(KL))
            noKU = merge(KU, noN, present(KU))
            noSPARSE = merge(SPARSE, D_SPARSE, present(SPARSE))
            noANORM = merge(ANORM, D_ANORM, present(ANORM))
            noPACK = merge(PACK, D_PACK, present(PACK))
            noLDA = merge(LDA, noM, present(LDA))
            if(present(IWORK)) then
                allocate(noIWORK(size(IWORK)))
                noIWORK = IWORK
            end if
            if(present(IPIVOT)) then
                allocate(noIPIVOT(size(IPIVOT)))
                noIPIVOT = IPIVOT
            end if

            call dlatmr(noM, noN, noDIST, noISEED, SYM, noD, noMODE, noCOND, noDMAX, noRSIGN, noGRADE, noDL, noMODEL, noCONDL, &
                        noDR, noMODER, noCONDR, noPIVTNG, noIPIVOT, noKL, noKU, noSPARSE, noANORM, noPACK, A, noLDA, noIWORK, INFO)

        end subroutine

        subroutine zheev_wrap(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)

            character, optional :: JOBZ, UPLO
            integer, optional :: N, LDA, LWORK
            complex*16, dimension(:, :) :: A
            double precision :: W(:)
            complex*16, allocatable, optional :: WORK(:)
            double precision, optional ::RWORK(:)
            integer :: INFO

            character :: noJOBZ, noUPLO
            integer :: noN, noLDA, noLWORK
            double precision, allocatable :: noRWORK(:)
            complex*16, allocatable :: noWORK(:)

            if(debug) print *, 'zheev was called'

            noJOBZ = merge(JOBZ, D_JOBZ, present(JOBZ))
            noUPLO = merge(UPLO, D_UPLO, present(UPLO))
            noN = merge(N, ubound(A, 1), present(N))
            noLDA = merge(LDA, ubound(A, 2), present(LDA))
            noLWORK = merge(LWORK, D_LWORK, present(LWORK))

            allocate(noRWORK(merge(size(RWORK), max(1, 3*noN-2), present(RWORK))))

            !> if lwork=-1, the first call determines the optimal value for
            !! lwork and stores it in work(1). Then we call again the subroutine
            !! with lwork optimal value
            !!
            if(noLWORK == -1) then
                allocate(noWORK(merge(size(WORK), 1, present(WORK))))
                call zheev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, noRWORK, INFO)
                noLWORK = int(noWORK(1))
                deallocate(noWORK)
                allocate(noWORK(noLWORK))

                call zheev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, noRWORK, INFO)
            else
                allocate(noWORK(noLWORK))
                call zheev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, noRWORK, INFO)
            end if

            if(present(WORK)) then
                deallocate(WORK)
                allocate(WORK(noLWORK))
                WORK = noWORK
            end if

            if(present(RWORK)) RWORK = noRWORK


        end subroutine

        subroutine dsyev_wrap(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

            character, optional :: JOBZ, UPLO
            integer, optional :: N, LDA, LWORK
            double precision, dimension(:, :) :: A
            double precision :: W(:)
            double precision, allocatable, optional :: WORK(:)
            integer :: INFO

            character :: noJOBZ, noUPLO
            integer :: noN, noLDA, noLWORK
            double precision, allocatable :: noWORK(:)

            if(debug) print *, 'dsyev was called'

            noJOBZ = merge(JOBZ, D_JOBZ, present(JOBZ))
            noUPLO = merge(UPLO, D_UPLO, present(UPLO))
            noN = merge(N, ubound(A, 1), present(N))
            noLDA = merge(LDA, ubound(A, 2), present(LDA))
            noLWORK = merge(LWORK, D_LWORK, present(LWORK))

            !> if lwork=-1, the first call determines the optimal value for
            !! lwork and stores it in work(1). Then we call again the subroutine
            !! with lwork optimal value
            !!
            if(noLWORK == -1) then
                allocate(noWORK(merge(size(WORK), 1, present(WORK))))
                call dsyev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, INFO)
                noLWORK = int(noWORK(1))
                deallocate(noWORK)
                allocate(noWORK(noLWORK))

                call dsyev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, INFO)
            else
                allocate(noWORK(noLWORK))
                call dsyev(noJOBZ, noUPLO, noN, A, noLDA, W, noWORK, noLWORK, INFO)
            end if

            if(present(WORK)) then
                deallocate(WORK)
                allocate(WORK(noLWORK))
                WORK = noWORK
            end if

        end subroutine

end module