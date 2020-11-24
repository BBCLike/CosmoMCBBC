MODULE BBC
    USE CosmologyTypes
    USE settings
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    IMPLICIT NONE

    type, extends(TCosmoCalcLikelihood) :: BBCLikelihood
        contains
        procedure :: LogLike => BBC_LnLike
    end type BBCLikelihood

    character(LEN=*), parameter :: BBC_version =  'December_2020'
    REAL(mcp), PARAMETER, PRIVATE :: inv_twoPI = 1.0_mcp / twopi
    CHARACTER, PARAMETER, PRIVATE :: uplo = 'U' !For LAPACK
    INTEGER, PARAMETER, PRIVATE :: max_idisp_datasets = 10
    INTEGER, PARAMETER, PRIVATE :: snnamelen = 12
    REAL(mcp), PARAMETER, PRIVATE :: h0cfac = 5*LOG10( 100.0/299792.458 )

    !Variables we will try to get from the ini file
    CHARACTER(LEN=30), PRIVATE :: name !Name of data set
    
    !Supernova data type
    TYPE, PRIVATE :: supernova
        CHARACTER(LEN=snnamelen) :: name  !The name of the SN
        REAL(mcp) :: zhel, zcmb   !The heliocentric and CMB frame redshifts
        REAL(mcp) :: mag  !The K-corrected peak magnitude
    END TYPE supernova

    INTEGER, PUBLIC :: nsn  !Number of supernovae
    TYPE(supernova), ALLOCATABLE, PRIVATE :: sndata(:) !Supernova data
    REAL(mcp), ALLOCATABLE, PRIVATE :: covmat(:,:)
    REAL(mcp), ALLOCATABLE, PRIVATE :: lumdists(:) !Other convenience variables

    LOGICAL, PUBLIC :: BBC_read = .FALSE.
    LOGICAL, PUBLIC :: BBC_prepped = .FALSE.

    PRIVATE :: count_lines, read_BBC_data, read_cov_matrix
    PUBLIC :: BBC_prep, BBC_LnLike, BBC_cleanup, read_BBC_dataset, BBCLikelihood_Add

    CONTAINS

    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------

    !Counts the number of lines in an open file attached to lun,
    !returning the number of lines in lines and the number of
    !non-comment lines in noncommentlines, where a comment line
    !is defined to start with a #. The file is rewound on exit
    SUBROUTINE count_lines( lun, lines, noncommentlines )
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

    !Reads the covariance matrix from a file, given the filename
    !and the number of elements to expect
    !There are two possible formats supported
    !These are: as one big block, and then as n by n individual elements
    !The number of lines has to be the same as the number of SN, and
    !they have to be in the same order. Copied from settings::ReadMatrix
    SUBROUTINE read_cov_matrix(filename, mat, n)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: n
        REAL(mcp), INTENT(OUT) :: mat(n,n)
        INTEGER :: j,k, file_unit, nfile
        REAL(mcp) :: tmp

        IF (Feedback > 2) WRITE(*,*) 'reading: '//trim(filename)
        OPEN( newunit=file_unit, FILE=TRIM(filename), FORM='formatted', &
            STATUS='old', ERR = 500 )

        READ (file_unit, '(I5)', END=200, ERR=100) nfile
        IF (nfile /= n) THEN
            WRITE (*,'("For file ",A," expected size ",I5," got ",I5)') &
                TRIM(filename), n, nfile
            STOP
        ENDIF

        DO j=1,n
            READ (file_unit,*, end = 200, err=100) mat(j,1:n)
        ENDDO

        GOTO 120

    100 REWIND(file_unit)  !Try other possible format
        READ (file_unit, '(I5)', END=200, ERR=100) nfile

        DO j=1,n
            DO k=1,n
                READ (file_unit,*, end = 200) mat(j,k)
            END DO
        END DO

    120 READ (file_unit,*, err = 150, end =150) tmp
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
    !------------------------------------------------------------
    !------------------------------------------------------------
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
        this%LikelihoodType = 'SN'
        this%name='BBC'
        this%needs_background_functions = .true.
        this%version = Ini%Read_String_Default('BBC_version',BBC_version)
        call this%loadParamNames(trim(DataDir)//'BBC.paramnames')
        call LikeList%Add(this)
        BBC_filename = Ini%Read_String_Default('BBC_dataset',trim(DataDir)//'BBC.dataset')
        CALL read_BBC_dataset(BBC_filename)
        CALL BBC_prep
    end subroutine BBCLikelihood_Add

    !------------------------------------------------------------
    !Reads in a supernova data file, given knowledge of the number
    !of lines to expect.  Ignores lines that start with #.
    ! Input arguments:
    !  lun              The lun number of the file to read.  Must be already open
    !  nlines           The number of lines to expect in the file
    !  nnoncommentlines The number of non-comment lines in the file
    ! Output arguments:
    !  sndata           The returned SN data, of length nnoncommentlines
    !The file is not rewound on exit
    !------------------------------------------------------------
    SUBROUTINE read_BBC_data( lun, nlines, nnoncommentlines, sndata )
        IMPLICIT NONE
        INTEGER, INTENT(in) :: lun, nlines, nnoncommentlines
        TYPE(supernova), INTENT(out) :: sndata(nnoncommentlines)
        CHARACTER(LEN=80) :: inline, shiftline
        INTEGER:: i,count
        LOGICAL :: opened
        INTRINSIC ADJUSTL

        INQUIRE( lun, OPENED=opened )
        IF (.NOT. opened) THEN
            WRITE(*,*) "File is not open in count_lines"
            STOP
        ENDIF

        count = 1
        DO i=1,nlines
            !Read in line non-advancing
            READ (lun, '(A)', ERR = 20, END = 20) inline
            shiftline = ADJUSTL( inline )
            IF (shiftline(1:1) .EQ. '#') CYCLE

    20      BACKSPACE lun
  
            READ (lun, *, ERR=60, END=50) sndata(count)%name, sndata(count)%zcmb, &
                sndata(count)%zhel, sndata(count)%mag
            count = count+1
        END DO

        RETURN

    50  WRITE(*,'("File ended unexpectedly on line ",I3," expecting ",I3)') i,nlines
        STOP

    60  WRITE(*,*) 'Error reading in input data with: ',inline
        STOP

    END SUBROUTINE read_BBC_data

    !------------------------------------------------------------
    ! The public interface to reading data files
    ! This gets information from the .ini file and reads the data file
    ! Arguments:
    !  filename The name of the .ini file specifying the SN dataset
    !------------------------------------------------------------
    SUBROUTINE read_BBC_dataset( filename )
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(in) :: filename
        CHARACTER(LEN=:), allocatable :: covfile
        CHARACTER(LEN=:), allocatable :: data_file
        INTEGER :: nlines
        Type(TSettingIni) :: Ini
        integer file_unit

        IF (BBC_read) THEN
            STOP 'Error -- BBC data already read'
        ENDIF

        !Process the Ini file
        CALL Ini%Open(filename)

        name = Ini%Read_String( 'name', .FALSE. )
        data_file = Ini%Read_String_Default('data_file', trim(DataDir)//'BBC_lcparams.txt')
        !Now read the actual SN data
        OPEN( newunit=file_unit, FILE=TRIM(data_file), FORM='formatted', STATUS='old', ERR = 500 )
        CALL count_lines( file_unit, nlines, nsn ) !Find the number of lines
        ALLOCATE( sndata(nsn) )
        CALL read_BBC_data( file_unit, nlines, nsn, sndata )
        CLOSE( file_unit )

        covfile = Ini%Read_String( 'covmat_file',.TRUE. )
        ALLOCATE( covmat( nsn, nsn ) )
        CALL read_cov_matrix( covfile, covmat, nsn )
        CALL Ini%Close()

        BBC_read = .TRUE.
        RETURN

    500 WRITE(*,*) 'Error reading ' // data_file
        STOP
    END SUBROUTINE read_BBC_dataset

    SUBROUTINE invert_covariance_matrix(invcovmat, status)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: status
        REAL(mcp) :: invcovmat(:,:)

        invcovmat = covmat
        !Factor into Cholesky form, overwriting the input matrix
        CALL DPOTRF(uplo,nsn,invcovmat,nsn,status)
        IF ( status .NE. 0 ) THEN
            STOP
        END IF
        !Note that DPOTRI only makes half of the matrix correct,
        !so we have to be careful in what follows
        CALL DPOTRI(uplo, nsn, invcovmat, nsn, status)
        IF ( status .NE. 0 ) THEN
            STOP
        END IF
    END SUBROUTINE invert_covariance_matrix

    SUBROUTINE BBC_prep
        IMPLICIT NONE

        IF (.NOT. BBC_read) THEN
            STOP 'BBC data was not read in'
        ENDIF
        IF (nsn < 1) THEN
            STOP 'No BBC data read'
        ENDIF
        
        ALLOCATE(lumdists(nsn))

        BBC_prepped = .TRUE.

        RETURN

        STOP
    END SUBROUTINE BBC_prep

    SUBROUTINE BBC_cleanup
        IF ( ALLOCATED( sndata ) ) DEALLOCATE( sndata )
        IF ( ALLOCATED( lumdists ) ) DEALLOCATE( lumdists )
        IF ( ALLOCATED( covmat ) ) DEALLOCATE( covmat )
        BBC_prepped = .FALSE.
    END SUBROUTINE BBC_cleanup

    FUNCTION BBC_LnLike(this, CMB, Theory, DataParams)
        Class(BBCLikelihood) :: this
        Class(CMBParams) CMB
        Class(TCosmoTheoryPredictions), target :: Theory
        real(mcp) DataParams(:)
        REAL(mcp) :: BBC_LnLike
        real(mcp) zhel, zcmb
        integer i

        BBC_LnLike = logZero

        !Make sure we're ready to actually do this
        IF (.NOT. BBC_read) THEN
            STOP 'BBC data not read in; must be by this point'
        ENDIF
        IF (.NOT. BBC_prepped ) THEN
            STOP 'BBC data not prepped; run BBC_prep'
        ENDIF

        !Get the luminosity distances.
        DO i=1,nsn
            zhel = sndata(i)%zhel
            zcmb = sndata(i)%zcmb
            lumdists(i) = 5.0*LOG10((1.0+zhel)*(1.0+zcmb)*this%Calculator%AngularDiameterDistance(zcmb))
        ENDDO
        BBC_LnLike=0.0
    END FUNCTION BBC_LnLike

END MODULE BBC
