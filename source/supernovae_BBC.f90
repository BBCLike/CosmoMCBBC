    MODULE BBC
    USE CosmologyTypes
    USE settings
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    IMPLICIT NONE

    !Modified by AL to have option of internal alpha, beta marginalization
    logical :: BBC_marginalize = .false.
    REAL(mcp), allocatable :: BBC_marge_grid(:), alpha_grid(:),beta_grid(:)
    integer :: BBC_marge_steps = 0
    real(mcp) BBC_step_width_alpha, BBC_step_width_beta
    real(mcp), parameter :: BBC_alpha_center =  0.14
    real(mcp), parameter :: BBC_beta_center = 3.123
    integer :: BBC_int_points = 1

    type, extends(TCosmoCalcLikelihood) :: BBCLikelihood
    contains
    procedure :: LogLike => BBC_LnLike
    end type BBCLikelihood

    integer, parameter :: dl = mcp

    character(LEN=*), parameter :: BBC_version =  'December_2013'
    logical, parameter :: allow_inv_cache = .false. !AL inverse cache does not work.. have not checked why.

    !Constants
    REAL(dl), PARAMETER, PRIVATE :: inv_twoPI = 1.0_dl / twopi
    CHARACTER, PARAMETER, PRIVATE :: uplo = 'U' !For LAPACK
    INTEGER, PARAMETER, PRIVATE :: max_idisp_datasets = 10
    INTEGER, PARAMETER, PRIVATE :: snnamelen = 12
    REAL(dl), PARAMETER, PRIVATE :: h0cfac = 5*LOG10( 100.0/299792.458 )
    REAL(dl), PARAMETER, PRIVATE :: alphatol = 1E-10_dl, betatol = 1E-10_dl

    !Variables we will try to get from the ini file
    CHARACTER(LEN=30), PRIVATE :: name !Name of data set
    REAL(dl), PRIVATE :: pecz !Peculiar velocity error in z
    REAL(dl), DIMENSION( max_idisp_datasets ) :: intrinsicdisp !In magnitudes

    !Variables having to do with optional two-scripmt fit based
    ! on thirdvar cut
    LOGICAL, PRIVATE :: twoscriptmfit !Carry out two scriptm fit
    LOGICAL, PRIVATE :: has_thirdvar  !Data has third variable
    REAL(dl), PRIVATE :: scriptmcut !Cut in thirdvar between two scriptms

    !Supernova data type
    TYPE, PRIVATE :: supernova
        CHARACTER(LEN=snnamelen) :: name  !The name of the SN
        REAL(dl) :: zhel, zcmb    !The heliocentric and CMB frame redshifts
        REAL(dl) :: z_var         !The variance of the redshift
        REAL(dl) :: mag           !The K-corrected peak magnitude
        REAL(dl) :: mag_var       !The variance of mag
        REAL(dl) :: stretch       !The light-curve fit stretch parameter
        REAL(dl) :: stretch_var   !The variance in the stretch
        REAL(dl) :: colour        !The colour of the SN
        REAL(dl) :: colour_var    !The variance of colour
        REAL(dl) :: thirdvar      !Third variable for scripm split
        REAL(dl) :: thirdvar_var  !Variance in thirdvar
        REAL(dl) :: cov_mag_stretch !Covariance between mag and stretch
        REAL(dl) :: cov_mag_colour  !Covariance between mag and colour
        REAL(dl) :: cov_stretch_colour !Covariance between stretch and colour
        LOGICAL :: has_absdist    !This SN has an absolute distance
        INTEGER  :: dataset       !Subset identifier if subset dependent intrinsic disp is used
    END TYPE supernova

    INTEGER, PUBLIC :: nsn  !Number of supernovae
    TYPE( supernova ), ALLOCATABLE, PRIVATE :: sndata(:)  !Supernova data
    !Stores the parts of the error that can be pre-calculated
    REAL(dl), ALLOCATABLE, PRIVATE :: pre_vars(:)
    !Arrays which have 1 for SN in set 1 (A1) or 2 (A2).  For twoscriptm fit
    REAL(dl), ALLOCATABLE, PRIVATE :: A1(:), A2(:)

    !Covariance matrix stuff
    ! If we have no covariance matrix at all, diag_errors is .TRUE.
    LOGICAL, PRIVATE :: diag_errors =        .TRUE.

    !Which components of the covariance matrix do we have
    LOGICAL, PRIVATE :: has_mag_covmat =            .FALSE.
    LOGICAL, PRIVATE :: has_stretch_covmat =        .FALSE.
    LOGICAL, PRIVATE :: has_colour_covmat =         .FALSE.
    LOGICAL, PRIVATE :: has_mag_stretch_covmat =    .FALSE.
    LOGICAL, PRIVATE :: has_mag_colour_covmat =     .FALSE.
    LOGICAL, PRIVATE :: has_stretch_colour_covmat = .FALSE.
    LOGICAL, PRIVATE :: alphabeta_covmat =          .FALSE.
    REAL(dl), ALLOCATABLE, PRIVATE :: mag_covmat(:,:), stretch_covmat(:,:)
    REAL(dl), ALLOCATABLE, PRIVATE :: colour_covmat(:,:), mag_stretch_covmat(:,:)
    REAL(dl), ALLOCATABLE, PRIVATE :: mag_colour_covmat(:,:)
    REAL(dl), ALLOCATABLE, PRIVATE :: stretch_colour_covmat(:,:)

    !Structure for holding absolute distance information for SN
    LOGICAL, PRIVATE :: has_absdist =     .FALSE.
    INTEGER, PRIVATE :: nabsdist =         0
    TYPE, PRIVATE :: supernova_absdist
        CHARACTER(LEN=snnamelen) :: name  !The name of the SN
        REAL(dl) :: dl             !Distance in Mpc
        INTEGER :: index           !Index into sndata
    END TYPE supernova_absdist
    TYPE( supernova_absdist ), ALLOCATABLE, PRIVATE :: snabsdist(:)

    !Other convenience variables
    REAL(dl), ALLOCATABLE, PRIVATE :: lumdists(:)
    REAL(dl), PRIVATE :: alpha_prev, beta_prev

    LOGICAL, PRIVATE :: first_inversion
    LOGICAL, PUBLIC :: BBC_read = .FALSE.
    LOGICAL, PUBLIC :: BBC_prepped = .FALSE.

    PRIVATE :: count_lines, read_BBC_lc_data, read_cov_matrix
    PRIVATE :: read_BBC_absdist_data, match_BBC_absdist_indices
    PUBLIC :: BBC_prep, BBC_LnLike, BBC_cleanup, read_BBC_dataset,BBCLikelihood_Add

    CONTAINS


    subroutine BBCLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TIniFile) :: ini
    Type(BBCLikelihood), pointer :: this
    character (LEN=:), allocatable:: BBC_filename
    integer alpha_i, beta_i

    if (.not. Ini%Read_Logical('use_BBC',.false.)) return

    allocate(this)
    this%LikelihoodType = 'SN'
    this%name='BBC'
    this%needs_background_functions = .true.
    this%version = Ini%Read_String_Default('BBC_version',BBC_version)
    BBC_marginalize = Ini%Read_Logical('BBC_marginalize',.false.)
    if (BBC_marginalize) then
        BBC_marge_steps = Ini%Read_Int('BBC_marge_steps',7)
        BBC_step_width_alpha = Ini%Read_Double('BBC_step_width_alpha',0.003d0)
        BBC_step_width_beta = Ini%Read_Double('BBC_step_width_beta',0.04d0)
        BBC_int_points=0
        allocate(alpha_grid((2*BBC_marge_steps+1)**2))
        allocate(beta_grid((2*BBC_marge_steps+1)**2))
        do alpha_i = - BBC_marge_steps, BBC_marge_steps
            do beta_i = - BBC_marge_steps, BBC_marge_steps
                if (alpha_i**2 + beta_i**2 <= BBC_marge_steps**2) then
                    BBC_int_points=BBC_int_points+1
                    alpha_grid(BBC_int_points) = BBC_alpha_center + alpha_i* BBC_step_width_alpha
                    beta_grid(BBC_int_points)  = BBC_beta_center + beta_i* BBC_step_width_beta
                end if
            end do
        end do
        allocate(BBC_marge_grid(BBC_int_points))
    else
        call this%loadParamNames(trim(DataDir)//'BBC.paramnames')
    end if
    call LikeList%Add(this)
    BBC_filename = Ini%Read_String_Default('BBC_dataset',trim(DataDir)//'BBC.dataset')
    CALL read_BBC_dataset( BBC_filename )
    CALL BBC_prep
    If (Feedback>0) WRITE(*,*) 'read BBC dataset '//trim(BBC_filename)

    end subroutine BBCLikelihood_Add

    !Counts the number of lines in an open file attached to lun,
    ! returning the number of lines in lines and the number of
    ! non-comment lines in noncommentlines, where a comment line
    ! is defined to start with a #
    !The file is rewound on exit
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
    ! and the number of elements to expect
    !There are two possible formats supported
    ! These are: as one big block, and then as n by n individual elements
    ! The number of lines has to be the same as the number of SN, and
    ! they have to be in the same order
    !Copied from settings::ReadMatrix
    SUBROUTINE read_cov_matrix(filename, mat, n)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(OUT) :: mat(n,n)
    INTEGER :: j,k, file_unit, nfile
    REAL(dl) :: tmp

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
    ! Reads in a supernova data file, given knowledge of the number
    !  of lines to expect.  Ignores lines that start with #.
    ! Input arguments:
    !  lun              The lun number of the file to read.  Must be already open
    !  nlines           The number of lines to expect in the file
    !  nnoncommentlines The number of non-comment lines in the file
    ! Output arguments:
    !  sndata           The returned SN data, of length nnoncommentlines
    ! Notes:
    !  The file is not rewound on exit
    !------------------------------------------------------------
    SUBROUTINE read_BBC_lc_data( lun, nlines, nnoncommentlines, sndata )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: lun, nlines, nnoncommentlines
    TYPE(supernova), INTENT(out) :: sndata(nnoncommentlines)

    CHARACTER(LEN=80) :: inline, shiftline
    INTEGER:: i,count
    REAL :: dz, dm, ds, dc, dt
    LOGICAL :: opened

    INTRINSIC ADJUSTL

    INQUIRE( lun, OPENED=opened )
    IF (.NOT. opened) THEN
        WRITE(*,*) "File is not open in count_lines"
        STOP
    ENDIF

    count = 1
    has_thirdvar = .FALSE.
    sndata%has_absdist = .FALSE.
    DO i=1,nlines
        !Read in line non-advancing
        READ (lun, '(A)', ERR = 20, END = 20) inline
        shiftline = ADJUSTL( inline )
        IF (shiftline(1:1) .EQ. '#') CYCLE

        BACKSPACE lun

        !We have a few formats to try.  First, there is the very
        ! long format with thirdvar and dataset.  If that fails,
        ! try without data set.  If that fails, try without
        ! thirdvar but with dataset, and finally with neither

        !A further complication is that if one line has thirdvar,
        ! they had better all have them or else ugliness will probably
        ! result
        READ (lun, *, ERR=20, END=20) &
            sndata(count)%name, sndata(count)%zcmb, sndata(count)%zhel,&
            dz, sndata(count)%mag, dm, sndata(count)%stretch, ds, &
            sndata(count)%colour,dc,sndata(count)%thirdvar, dt,&
            sndata(count)%cov_mag_stretch,&
            sndata(count)%cov_mag_colour,sndata(count)%cov_stretch_colour,&
            sndata(count)%dataset
        IF ( (count .GT. 1) .AND. (.NOT. has_thirdvar) ) THEN
            WRITE(*,*) "Problem with third variable read"
            STOP
        ENDIF
        has_thirdvar = .TRUE.
        GOTO 10  !Success

        !That didn't work. Try without dataset.  First, undo the
        ! previous.  It should be 2 records out of place because
        ! we read over into the next line
20      BACKSPACE lun
        BACKSPACE lun
        READ (lun, *, ERR=30, END=30) &
            sndata(count)%name, sndata(count)%zcmb, sndata(count)%zhel,&
            dz, sndata(count)%mag, dm, sndata(count)%stretch, ds, &
            sndata(count)%colour,dc,sndata(count)%thirdvar,dt,&
            sndata(count)%cov_mag_stretch,&
            sndata(count)%cov_mag_colour,sndata(count)%cov_stretch_colour
        IF ( (count .GT. 1) .AND. (.NOT. has_thirdvar) ) THEN
            WRITE(*,*) "Problem with third variable read"
            STOP
        ENDIF
        has_thirdvar = .TRUE.
        GOTO 10  !Success

        !Ok, maybe there's no thirdvar
30      BACKSPACE lun
        BACKSPACE lun
        READ (lun, *, ERR=40, END=40) &
            sndata(count)%name, sndata(count)%zcmb, sndata(count)%zhel,&
            dz, sndata(count)%mag, dm, sndata(count)%stretch, ds, &
            sndata(count)%colour,dc,sndata(count)%cov_mag_stretch,&
            sndata(count)%cov_mag_colour,sndata(count)%cov_stretch_colour,&
            sndata(count)%dataset
        IF ( (count .GT. 1) .AND. (has_thirdvar) ) THEN
            WRITE(*,*) "Problem with third variable read"
            STOP
        ENDIF
        sndata(count)%thirdvar = 0.0
        dt = 0.0
        sndata(count)%dataset = 0

        !Still?
        !Ok, maybe there's no thirdvar and no dataset
40      BACKSPACE lun
        BACKSPACE lun
        READ (lun, *, ERR=60, END=50) &
            sndata(count)%name, sndata(count)%zcmb, sndata(count)%zhel,&
            dz, sndata(count)%mag, dm, sndata(count)%stretch, ds, &
            sndata(count)%colour,dc,sndata(count)%cov_mag_stretch,&
            sndata(count)%cov_mag_colour,sndata(count)%cov_stretch_colour,&
            sndata(count)%dataset
        IF ( (count .GT. 1) .AND. (has_thirdvar) ) THEN
            WRITE(*,*) "Problem with third variable read"
            STOP
        ENDIF
        sndata(count)%thirdvar = 0.0
        dt = 0.0
        sndata(count)%dataset = 0

10      sndata(count)%z_var = dz**2
        sndata(count)%mag_var = dm**2
        sndata(count)%stretch_var = ds**2
        sndata(count)%colour_var = dc**2
        sndata(count)%thirdvar_var = dt**2
        !sndata(count)%thirdvar = 6
        count = count+1
    END DO
    RETURN

50  WRITE(*,'("File ended unexpectedly on line ",I3," expecting ",I3)') i,nlines
    STOP

60  WRITE(*,*) 'Error reading in input data with: ',inline
    STOP

    END SUBROUTINE read_BBC_lc_data

    !------------------------------------------------------------
    ! Read in absolute distance info, given knowledge of the number
    !  of lines to expect.  Ignores lines that start with #.
    ! Input arguments:
    !  lun              The lun number of the file to read.  Must be already open
    !  nlines           The number of lines to expect in the file
    !  nnoncommentlines The number of non-comment lines in the file
    ! Output arguments:
    !  snabsdist        The absolute distance data, of length nnoncommentlines
    ! Notes:
    !  The file is not rewound on exit
    !------------------------------------------------------------
    SUBROUTINE read_BBC_absdist_data( lun, nlines, nnoncommentlines, snabsdist )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: lun, nlines, nnoncommentlines
    TYPE(supernova_absdist), INTENT(out) :: snabsdist(nnoncommentlines)

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
        !Read in line non-advancing mode
        READ (lun, '(A)', ERR = 140, END = 130) inline
        shiftline = ADJUSTL( inline )
        IF (shiftline(1:1) .EQ. '#') CYCLE

        BACKSPACE lun

        READ (lun, *, ERR=140, END=130) &
            snabsdist(count)%name, snabsdist(count)%dl
        count = count+1
    END DO
    RETURN

130 WRITE(*,'("File ended unexpectedly on line ",I3," expecting ",I3)') i,nlines
    STOP

140 WRITE(*,*) 'Error reading in input data with: ',inline
    STOP

    END SUBROUTINE read_BBC_absdist_data

    !------------------------------------------------------------
    ! The public interface to reading data files
    ! This gets information from the .ini file and reads the data file
    ! Arguments:
    !  filename        The name of the .ini file specifying the SN dataset
    !------------------------------------------------------------
    SUBROUTINE read_BBC_dataset(filename )
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=:), allocatable :: covfile
    CHARACTER(LEN=:), allocatable :: data_file
    INTEGER :: nlines, i
    REAL(dl) :: idisp_zero !Value for unspecified dataset numbers
    LOGICAL, DIMENSION( max_idisp_datasets ) :: idispdataset
    Type(TSettingIni) :: Ini
    integer file_unit

    IF (BBC_read) STOP 'Error -- BBC data already read'

    !Process the Ini file
    CALL Ini%Open(filename)

    name = Ini%Read_String( 'name', .FALSE. )
    data_file = Ini%Read_String_Default('data_file',trim(DataDir)//'BBC_lcparams.txt')

    !Now read the actual SN data
    OPEN( newunit=file_unit, FILE=TRIM(data_file), FORM='formatted', &
        STATUS='old', ERR = 500 )
    !Find the number of lines
    CALL count_lines( file_unit, nlines, nsn )
    ALLOCATE( sndata(nsn) )
    CALL read_BBC_lc_data( file_unit, nlines, nsn, sndata )
    CLOSE( file_unit )

    !Handle covariance matrix stuff
    has_mag_covmat=Ini%Read_Logical( 'has_mag_covmat', .TRUE. )
   
    !First test for covmat
    IF ( has_mag_covmat ) THEN
        diag_errors = .FALSE.

        covfile = Ini%Read_String('mag_covmat_file',.TRUE.)
        ALLOCATE( mag_covmat( nsn, nsn ) )
        CALL read_cov_matrix( covfile, mag_covmat, nsn )
    ELSE
        diag_errors = .TRUE.
    END IF

    CALL Ini%Close()

    IF (Feedback > 1) THEN
        WRITE(*,'(" BBC dataset name: ",A)') TRIM(name)
        WRITE(*,'(" BBC data file: ",A)') TRIM(data_file)
        WRITE(*,'(" Number of SN read: ",I4)') nsn
    ENDIF

    first_inversion = .true.
    BBC_read = .TRUE.
    BBC_prepped = .FALSE.
    RETURN

500 WRITE(*,*) 'Error reading ' // data_file
    STOP

    END SUBROUTINE read_BBC_dataset

    SUBROUTINE invert_covariance_matrix(invcovmat, status)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: status
        REAL(dl) :: invcovmat(:,:)

        invcovmat = mag_covmat

        !Factor into Cholesky form, overwriting the input matrix
        CALL DPOTRF(uplo,nsn,invcovmat,nsn,status)
        IF ( status .NE. 0 ) THEN
            STOP
        END IF

        !Note that DPOTRI only makes half of the matrix correct,
        ! so we have to be careful in what follows
        CALL DPOTRI(uplo, nsn, invcovmat, nsn, status)
        IF ( status .NE. 0 ) THEN
            STOP
        END IF
    END SUBROUTINE invert_covariance_matrix

    SUBROUTINE BBC_prep
        IMPLICIT NONE
        INTEGER ::  i

        IF (.NOT. BBC_read) STOP 'BBC data was not read in'
        IF (nsn < 1) STOP 'No BBC data read'

        IF ( MAXVAL( sndata%dataset ) .GE. max_idisp_datasets ) THEN
            WRITE(*,*) 'Invalid dataset number ',MAXVAL(sndata%dataset)
            WRITE(*,*) ' Maximum allowed is ',max_idisp_datasets
        END IF
        IF ( MINVAL( sndata%dataset ) .LT. 0 ) THEN
            WRITE(*,*) 'Invalid dataset number ',MINVAL(sndata%dataset)
            WRITE(*,*) ' Maximum allowed is 0'
        END IF

        ALLOCATE(lumdists(nsn))

        BBC_prepped = .TRUE.
        first_inversion = .TRUE.
        
        RETURN

    500 WRITE(*,*) 'Error reading ' // datafile
        STOP
    END SUBROUTINE BBC_prep

    SUBROUTINE BBC_cleanup
        IF ( ALLOCATED( sndata ) ) DEALLOCATE( sndata )
        IF ( ALLOCATED( lumdists ) ) DEALLOCATE( lumdists )
        IF ( ALLOCATED( mag_covmat ) ) DEALLOCATE( mag_covmat )
        BBC_prepped = .FALSE.
    END SUBROUTINE BBC_cleanup

    FUNCTION BBC_LnLike(this, CMB, Theory, DataParams)
    Class(BBCLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    REAL(mcp) :: BBC_LnLike
    real(dl) zhel, zcmb
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
        lumdists(i) = 5.0* LOG10( (1.0+zhel)*(1.0+zcmb) * this%Calculator%AngularDiameterDistance(zcmb) )
    ENDDO


    alpha=DataParams(1)
    beta=DataParams(2)

    BBC_LnLike=BBC_alpha_beta_like(alpha, beta, lumdists)

    END FUNCTION BBC_LnLike

    END MODULE BBC
