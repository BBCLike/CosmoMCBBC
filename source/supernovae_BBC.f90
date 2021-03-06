MODULE BBC
    USE CosmologyTypes
    USE settings
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    IMPLICIT NONE

    type, extends(TCosmoCalcLikelihood) :: BBCLikelihood
        contains
        procedure :: LogLikeTheory => BBC_LnLike
    end type BBCLikelihood

    CHARACTER, PARAMETER, PRIVATE :: uplo = 'U' !For LAPACK

    !Supernova data type
    TYPE, PRIVATE :: supernova
        CHARACTER(LEN=12) :: name     !The name of the SN
        REAL(mcp) :: zhel             !The heliocentric frame redshift
        REAL(mcp) :: zcmb             !The CMB frame redshift
        REAL(mcp) :: mag              !Fitted m^b from JLA
        REAL(mcp) :: mag_err
        REAL(mcp) :: mu               !Distance Modulus (from BBC)
        REAL(mcp) :: mu_err
    END TYPE supernova

    INTEGER, PUBLIC :: nsn  !Number of supernovae
    TYPE(supernova), ALLOCATABLE, PRIVATE :: sndata(:) !Supernova data
    LOGICAL, PRIVATE :: FOUND_COVMAT_SYST = .FALSE.
    LOGICAL, PRIVATE :: FOUND_COVMAT_STAT = .FALSE.
    REAL(mcp), ALLOCATABLE, PRIVATE :: covmat_syst(:,:)
    REAL(mcp), ALLOCATABLE, PRIVATE :: covmat_stat(:,:)
    REAL(mcp), ALLOCATABLE, PRIVATE :: covmat_tot(:,:)
    REAL(mcp), ALLOCATABLE, PRIVATE :: covmat_inv_tot(:,:)

    ! THEORY LUM DIST IS A GLOBAL MODULE VARIABLE TO AVOID
    ! MALLOC/FREE AT EVERY CHAIN POINT
    REAL(mcp), ALLOCATABLE, PRIVATE :: lumdists(:)

    LOGICAL, PUBLIC :: BBC_read = .FALSE.
    LOGICAL, PUBLIC :: BBC_prepped = .FALSE.

    PRIVATE :: count_lines, read_BBC_data, read_cov_matrix, BBC_prep, read_BBC_dataset
    PUBLIC :: BBC_LnLike,  BBCLikelihood_Add

    CONTAINS

    ! SIMPLIFIED BACKGROUND EVOLUTION (TO COMPUTE H(Z)
    SUBROUTINE simplified_background(omm, grhok)
        REAL(mcp), intent(in) :: omm
        
    END SUBROUTINE

    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------

    !Counts the number of lines in an open file attached to lun,
    !returning the number of lines in lines and the number of
    !non-comment lines in noncommentlines, where a comment line
    !is defined to start with a #. The file is rewound on exit
    SUBROUTINE count_lines(lun, lines, noncommentlines)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: lun
        INTEGER, INTENT(out) :: lines, noncommentlines
        INTEGER, PARAMETER :: maxlines = 5000 !Maximum number allowed
        INTEGER :: i
        CHARACTER(LEN=80) :: inline, shiftline
        LOGICAL :: opened

        INTRINSIC ADJUSTL

        !Make sure the file is open
        INQUIRE( lun, OPENED=opened )
        IF (.NOT. opened) THEN
            WRITE(*,*) "File is not open in count_lines"
            STOP
        ENDIF

        !Now start reading
        lines = 0
        noncommentlines = 0
        DO i = 1, maxlines
            READ( lun, '(A)', ERR=2, END=100 ) inline
            lines = lines + 1
            shiftline = ADJUSTL( inline )
            IF ( shiftline(1:1) .NE. '#' ) noncommentlines = noncommentlines+1
        ENDDO
        GO TO 100

    2   WRITE(*,*) "Error reading input file in count_lines"
        STOP

    100 REWIND lun
    END SUBROUTINE count_lines

    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------

!------------------------------------------------------------
! Reads the covariance matrix from a file, given the filename
!------------------------------------------------------------
    SUBROUTINE read_cov_matrix(filename, mat, n)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: n
        REAL(mcp), INTENT(OUT) :: mat(n,n)
        INTEGER :: j, k, file_unit, nfile
        REAL(mcp) :: tmp

        OPEN(newunit=file_unit, FILE=TRIM(filename), FORM='formatted', STATUS='old', ERR=500)

        READ (file_unit, '(I5)', END=200, ERR=500) nfile
        IF (nfile /= n) THEN
            WRITE (*,'("File ",A," expected size ",I5," got ",I5)') TRIM(filename), n, nfile
            STOP
        ENDIF

        DO j=1,n
            DO k=1,n
                READ (file_unit, *, end=200) mat(j,k)
            END DO
        END DO

        READ (file_unit, *, err=150, end=150) tmp
        GOTO 200

    150 CLOSE(file_unit)
        RETURN

    200 WRITE (*,*) 'matrix file '//trim(filename)//' is the wrong size'
        WRITE (*,'("Expected: ",I5," by ",I5)') n,n
        STOP

    500 WRITE (*,*) 'Failed to open cov matrix file ' // TRIM(filename)
        STOP
    END SUBROUTINE read_cov_matrix

!------------------------------------------------------------
! Reads in a supernova data file, given knowledge of the
! number of lines to expect. Ignores lines that start with #.
!------------------------------------------------------------
    SUBROUTINE read_BBC_data(lun, nlines, nnoncommentlines, sndata)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: lun, nlines, nnoncommentlines
        TYPE(supernova), INTENT(out) :: sndata(nnoncommentlines)
        CHARACTER(LEN=80) :: inline, shiftline
        INTEGER:: i,count
        LOGICAL :: opened
        REAL(mcp) :: tmp
        INTRINSIC ADJUSTL

        INQUIRE(lun, OPENED=opened)
        IF (.NOT. opened) THEN
            WRITE(*,*) "File is not open in count_lines"
            STOP
        ENDIF

        count = 1
        DO i=1,nlines
            READ (lun, '(A)', ERR = 20, END = 20) inline
            shiftline = ADJUSTL(inline)
            IF (shiftline(1:1) .EQ. '#') THEN
                CYCLE
            ENDIF

    20      BACKSPACE lun

            READ (lun, *, ERR=60, END=50) sndata(count)%name, sndata(count)%zcmb, &
                sndata(count)%zhel, tmp, sndata(count)%mag, sndata(count)%mag_err
            sndata(count)%mu = sndata(count)%mag + 19.35 ! TODO: OFFSET
            sndata(count)%mu_err = sndata(count)%mag_err
            count = count+1
        END DO

        RETURN

    50  WRITE(*,'("File ended unexpectedly on line ",I3," expecting ",I3)') i,nlines
        STOP

    60  WRITE(*,*) 'Error reading in input data with: ',inline
        STOP
    END SUBROUTINE read_BBC_data

    SUBROUTINE BBC_prep
        IMPLICIT NONE
        integer :: i, j, k, status

        IF (.NOT. BBC_read) THEN
            STOP 'BBC data was not read in'
        ENDIF
        IF (nsn < 1) THEN
            STOP 'No BBC data read'
        ENDIF

        ALLOCATE(lumdists(nsn))

        IF (.NOT. FOUND_COVMAT_STAT) THEN
            DO j=1,nsn
                !BACKWARD COMPATIBILITY WITH LEGACY COSMOMC JLA SN CODE
                covmat_stat(j,j) = sndata(j)%mu_err*sndata(j)%mu_err
            END DO
        ENDIF
        
        DO i=1,nsn
            DO j=1,nsn
                covmat_tot(i,j) = covmat_stat(i,j) + covmat_syst(i,j)
                ! DPOTRF AND DPOTRI INV FUNCTIONS OVERWRITE covmat_inv_tot
                covmat_inv_tot(i,j) = covmat_tot(i,j)
            ENDDO
        ENDDO

        !Factor into Cholesky form, overwriting the input matrix
        CALL DPOTRF(uplo, nsn, covmat_inv_tot, nsn, status)
        IF (status .NE. 0) THEN
            STOP
        END IF
        !Note that DPOTRI only makes half of the matrix correct,
        !so we have to be careful in what follows
        CALL DPOTRI(uplo, nsn, covmat_inv_tot, nsn, status)
        IF (status .NE. 0) THEN
            STOP
        END IF

        DO k=1,nsn
            DO j=k+1,nsn
                covmat_inv_tot(j,k) = covmat_inv_tot(k,j)
            ENDDO
        ENDDO

        BBC_prepped = .TRUE.
        RETURN
    END SUBROUTINE BBC_prep

    SUBROUTINE read_BBC_dataset(filename)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(in) :: filename
        CHARACTER(LEN=:), allocatable :: covfile_syst
        CHARACTER(LEN=:), allocatable :: covfile_stat
        CHARACTER(LEN=:), allocatable :: data_file
        INTEGER :: i, j, nlines
        Type(TSettingIni) :: Ini
        integer file_unit

        IF (BBC_read) THEN
            STOP 'Error -- BBC data already read'
        ENDIF

        CALL Ini%Open(filename)
        data_file = Ini%Read_String('data_file')

        OPEN(newunit=file_unit, FILE=TRIM(data_file), FORM='formatted', STATUS='old', ERR=500)
        CALL count_lines(file_unit, nlines, nsn) !Find the number of lines

        ALLOCATE(sndata(nsn))
        CALL read_BBC_data(file_unit, nlines, nsn, sndata)
        CLOSE(file_unit)

        ALLOCATE(covmat_syst(nsn, nsn))     ! READ FROM FILE
        ALLOCATE(covmat_stat(nsn, nsn))     ! READ FROM FILE
        ALLOCATE(covmat_tot(nsn, nsn))      ! LOADED AT BBC PREP SUBROUTINE
        ALLOCATE(covmat_inv_tot(nsn, nsn))  ! LOADED AT BBC PREP SUBROUTINE
        DO i=1,nsn
            DO j=1,nsn
                covmat_syst(i,j) = 0.0
                covmat_stat(i,j) = 0.0
                covmat_tot(i,j) = 0.0
                covmat_inv_tot(i,j) = 0.0
            ENDDO
        ENDDO

        ! READ SYST COV
        covfile_syst = Ini%Read_String_Default('covmat_syst_file',default='',AllowBlank=.TRUE.,EnvDefault=.FALSE.)
        IF (covfile_syst .ne. '') THEN
            CALL read_cov_matrix(covfile_syst, covmat_syst, nsn)
            FOUND_COVMAT_SYST = .TRUE.
        ELSE
            write(*,*) 'BBC SN: SYS COV NOT FOUND, COV_SYST SET TO ZERO'
            FOUND_COVMAT_SYST = .FALSE.
        END IF

        ! READ STAT COV
        ! IF STAT FILE != EXISTS, THEN WE ASSUME STAT ERRORS = dm
        ! (dm = column on BBC DATA FILE - THIS ENSURES BACKWARDS COMPAT)
        covfile_stat = Ini%Read_String_Default('covmat_stat_file',default='',AllowBlank=.TRUE.,EnvDefault=.FALSE.)
        IF (covfile_stat .ne. '') THEN
            CALL read_cov_matrix(covfile_stat, covmat_stat, nsn)
            FOUND_COVMAT_STAT = .TRUE.
        ELSE
            write(*,*) 'BBC SN: STAT COV NOT FOUND, DIAG(COV_STAT) SET TO MU ERROR ON BBC DATA'
            FOUND_COVMAT_STAT = .FALSE.
        END IF

        IF(Ini%HasKey('mag_covmat_file')) THEN
            write(*,*) 'ERROR: OLD KEY mag_covmat_file IS NO LONGER VALID'
            write(*,*) 'USE THE KEYSET (covmat_stat_file, covmat_syst_file) INSTEAD'
            write(*,*) 'ABORT'
            STOP
        ENDIF

        CALL Ini%Close()

        BBC_read = .TRUE.
        RETURN

    500 WRITE(*,*) 'Error reading ' // data_file
        STOP
    END SUBROUTINE read_BBC_dataset

    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------

    subroutine BBCLikelihood_Add(LikeList, Ini)
        class(TLikelihoodList) :: LikeList
        class(TIniFile) :: ini
        Type(BBCLikelihood), pointer :: this
        character (LEN=:), allocatable:: BBC_filename

        if (.not. Ini%Read_Logical('use_BBC',.false.)) then
            return
        endif

        allocate(this)
        this%name='BBC'
        this%LikelihoodType = 'SN'
        this%needs_background_functions = .true.
        this%version = Ini%Read_String('BBC_version')
        call LikeList%Add(this)
        BBC_filename = Ini%Read_String('BBC_dataset')
        CALL read_BBC_dataset(BBC_filename)
        CALL BBC_prep
    end subroutine BBCLikelihood_Add

    FUNCTION BBC_LnLike(this, CMB)
        Class(BBCLikelihood) :: this
        Class(CMBParams) CMB
        real(mcp) :: ANALYTIC_MARG_A, ANALYTIC_MARG_B, ANALYTIC_MARG_C
        real(mcp) :: TMP_A, TMP_B, TMP_C, X
        real(mcp) :: zhel, zcmb, dl, BBC_LnLike, BBC_CHI2
        integer i, j

        BBC_LnLike = logZero

        !Make sure we're ready to actually do this
        IF (.NOT. BBC_read) THEN
            STOP 'BBC data not read in; must be by this point'
        ENDIF
        IF (.NOT. BBC_prepped ) THEN
            STOP 'BBC data not prepped; run BBC_prep'
        ENDIF

        DO i=1,nsn
            zhel = sndata(i)%zhel
            zcmb = sndata(i)%zcmb
            dl = (1.0+zhel)*(1.0+zcmb)*this%Calculator%AngularDiameterDistance(zcmb)
            lumdists(i) = 5.0*LOG10(dl/1e-5)
        ENDDO

        ANALYTIC_MARG_A = 0.0
        ANALYTIC_MARG_B = 0.0
        ANALYTIC_MARG_C = 0.0
        BBC_CHI2 = 0.0
        X = 5*LOG(CMB%H0/70.0)
        DO i=1,nsn
            TMP_A = 0
            TMP_B = 0
            TMP_C = 0
            DO j=1,nsn
                TMP_A = TMP_A + covmat_inv_tot(i,j)*((lumdists(j)+X)-sndata(j)%mu)
                TMP_B = TMP_B + covmat_inv_tot(i,j)
                TMP_C = TMP_C + covmat_inv_tot(i,j)
            ENDDO
            ANALYTIC_MARG_A = ANALYTIC_MARG_A + TMP_A*((lumdists(i)+X)-sndata(i)%mu)
            ANALYTIC_MARG_B = ANALYTIC_MARG_B + TMP_B*((lumdists(i)+X)-sndata(i)%mu)
            ANALYTIC_MARG_C = ANALYTIC_MARG_C + TMP_C
        ENDDO
        BBC_CHI2 = ANALYTIC_MARG_A + LOG(ANALYTIC_MARG_C/twopi) - ANALYTIC_MARG_B**2/ANALYTIC_MARG_C
        BBC_LnLike = 0.5*BBC_CHI2
    END FUNCTION BBC_LnLike

END MODULE BBC
