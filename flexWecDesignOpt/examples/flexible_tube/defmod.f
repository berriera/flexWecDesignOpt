      PROGRAM DEFMOD
!-----------------------------------------------------------------------
!
!     COPYRIGHT (C) 2007 WAMIT Incorporated
!     COPYRIGHT (C) 1992-1998 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
!
!-----------------------------------------------------------------------
!
!     Main program file :  DEFMOD.F
!
!-----------------------------------------------------------------------
!
!     WAMIT Version : 6.3
!
!-----------------------------------------------------------------------
!
!     DESCRIPTION : DEFMOD is a WAMIT pre-processor for defining special
!     modes of motion.  DEFMOD reads the input file gdf.PRE with coordin
!     of the panel centroids and normal vector, written by POTEN, and
!     constructs user-defined modes as specified by the vectors VEL,
!     and extended hydrostatic matric C(i,j).  Both are
!     output to the file gdf.MOD for input to WAMIT.
!
!     X, Y and Z in .PRE are panel centroid XCT in global coordinates 
!     system and VEL in .PRE is in body coordiantes system
!-----------------------------------------------------------------------
!
!     Subroutines called by the main program PREMOD
!     are listed below alphabetically. Indented subroutine names
!     indicate calls made by leading subroutine name
!
!   Name      File      Description
!   ------------------------------------------------------------------
!
!    RDNAMP   DEFMOD    Read FNAMES.WAM file
!    RBLANK     "       Removes leading blanks from filenames
!    NEWNAM     "       Derives output filename
!    DATTIM     "       Inputs system date and time
!    HEADRD     "       Outputs header
!    DEFINE     "       User-specified subroutine to evaluate VEL for
!                            each new mode
!-----------------------------------------------------------------------
!
!   MOD 1/98  removed parameters MAXBDY, MAXDFR.  Installed generic
!             headrd subroutine.  Modified fixed format source to
!             fixed/free format.  DATTIM revised for generic use.
!
!       1/07   hydrostatic pressure on dipole panels not added to C
!              Subroutines DEFILE1 to DEFINE4 are included in the
!              file and DEFINE2 is named as DEFINE to be called
!              by the main program. DEFINE0 is added.
!
!-----------------------------------------------------------------------
!     Strings for date and time
!-----------------------------------------------------------------------
      CHARACTER*11 ADATE
      CHARACTER*8 ATIME
!-----------------------------------------------------------------------
!     Strings for file headers and WAMIT version number
!-----------------------------------------------------------------------
      CHARACTER*72 HEAD
      CHARACTER*11 VERHDR
!-----------------------------------------------------------------------
!     Arrays for filenames, extensions, and device numbers
!-----------------------------------------------------------------------
      CHARACTER*20 FNAMES,FILEN(3)
      CHARACTER*4 EXTN(3)
      INTEGER IFILEN(3)
!-----------------------------------------------------------------------
!     Arrays for normal velocity of modes, hydrostatic matrix, and
!     vector displacement (U,V,W) of new modes, symmetry of modes
!-----------------------------------------------------------------------
      REAL VEL(:),C(:,:),U(:),V(:),W(:)
      INTEGER ISYM(:)
      INTEGER IDIPOL
!-----------------------------------------------------------------------
!     Array for products of symmetry factors
!     First pair of indices are for rigid-body modes (1-6)
!     Second pair of indices are for symmetry indices ISX,ISY
!-----------------------------------------------------------------------
      INTEGER MIJXY(6,6,0:1,0:1)
      ALLOCATABLE :: VEL,C,U,V,W,ISYM
!-----------------------------------------------------------------------
!     Assign version number for header
!-----------------------------------------------------------------------
      DATA VERHDR / '6.3' /
!-----------------------------------------------------------------------
!     Assign I/O default filenames and extensions
!-----------------------------------------------------------------------
      DATA FNAMES / 'fnames.wam' /
      DATA EXTN/ '.gdf','.pre','.mod' /
!-----------------------------------------------------------------------
!     Assign device numbers to output files:
!     INTIP=keyboard, INTOP=monitor
!     IFILEN are I/O files with names FILEN
!-----------------------------------------------------------------------
      DATA INTIP /5/, INTOP /6/
      DATA IFILEN / 1,2,3 /
      DATA ZERO/ 0.0 /, ONE / 1.0 /
!-----------------------------------------------------------------------
!     The following table gives the symmetry factors for integration of
!     products of two modes (i,j = 1,2,3,4,5,6).
!-----------------------------------------------------------------------
      DATA MIJXY /                                                      &
!-----------------------------------------------------------------------
!     ISX=0, ISY=0   (integration over complete body surface)
!-----------------------------------------------------------------------
     &     1,1,1,1,1,1,                                                 &
     &     1,1,1,1,1,1,                                                 &
     &     1,1,1,1,1,1,                                                 &
     &     1,1,1,1,1,1,                                                 &
     &     1,1,1,1,1,1,                                                 &
     &     1,1,1,1,1,1,                                                 &
!-----------------------------------------------------------------------
!     ISX=1, ISY=0   (integration over half of body, one side of x=0)
!-----------------------------------------------------------------------
     &     2,0,0,0,2,2,                                                 &
     &     0,2,2,2,0,0,                                                 &
     &     0,2,2,2,0,0,                                                 &
     &     0,2,2,2,0,0,                                                 &
     &     2,0,0,0,2,2,                                                 &
     &     2,0,0,0,2,2,                                                 &
!-----------------------------------------------------------------------
!     ISX=0, ISY=1   (integration over half of body, one side of y=0)
!-----------------------------------------------------------------------
     &     2,0,2,0,2,0,                                                 &
     &     0,2,0,2,0,2,                                                 &
     &     2,0,2,0,2,0,                                                 &
     &     0,2,0,2,0,2,                                                 &
     &     2,0,2,0,2,0,                                                 &
     &     0,2,0,2,0,2,                                                 &
!-----------------------------------------------------------------------
!     ISX=1, ISY=1   (integration over one quadrant of body)
!-----------------------------------------------------------------------
     &     4,0,0,0,4,0,                                                 &
     &     0,4,0,4,0,0,                                                 &
     &     0,0,4,0,0,0,                                                 &
     &     0,4,0,4,0,0,                                                 &
     &     4,0,0,0,4,0,                                                 &
     &     0,0,0,0,0,4/
!-----------------------------------------------------------------------
!     Write copyright header; read input files
!-----------------------------------------------------------------------
      CALL HEADRD(INTOP,VERHDR)
!-----------------------------------------------------------------------
!     Read file FNAMES.WAM with list of input filenames if it exists
!        If not, interactively input file name of PRE file
!-----------------------------------------------------------------------
      CALL RDNAMD(IFILEN,INTIP,INTOP,FILEN,FNAMES,EXTN)
!-----------------------------------------------------------------------
!     Call to system date and time routine
!-----------------------------------------------------------------------
      CALL DATTIM(ADATE,ATIME)
!-----------------------------------------------------------------------
!     Open I/O files and write header lines
!-----------------------------------------------------------------------
      OPEN (IFILEN(2), FILE=FILEN(2))
      OPEN (IFILEN(3), FILE=FILEN(3))
      READ (IFILEN(2),901,END=98) HEAD
901   FORMAT (A72)
      WRITE (INTOP,901) HEAD
      WRITE (IFILEN(3),901) HEAD
903   FORMAT (1X,A72)
      WRITE (INTOP,908) FILEN(3),ADATE,ATIME
      WRITE (IFILEN(3),908) FILEN(3),ADATE,ATIME
908   FORMAT (1X, 'DEFMOD Output, Filename: ',A,4X,A,2X,A)
!-----------------------------------------------------------------------
!     Read symmetry indices, number of panels, global NDFR (all bodies)
!     Allocate all arrays.
!-----------------------------------------------------------------------
      READ (IFILEN(2),*,END=98) ISX,ISY,NEQN,MDFR
      ALLOCATE (VEL(MDFR),C(MDFR,MDFR),U(MDFR),V(MDFR),W(MDFR),         &
     &              ISYM(MDFR))
!-----------------------------------------------------------------------
!     Initialize hydrostatic data for entire square matrix
!-----------------------------------------------------------------------
      C=ZERO
!-----------------------------------------------------------------------
!     Initialize indices ISYM for rigid-body modes
!     (ISYM denotes symmetry of new modes in terms of rigid-body modes)
!-----------------------------------------------------------------------
      DO I=1,6
        ISYM(I)=I
      END DO
!-----------------------------------------------------------------------
!     Initialize vertical velocities W for rigid-body modes
!     (pitch and roll are evaluated in loop over panels)
!-----------------------------------------------------------------------
      W(1)=ZERO
      W(2)=ZERO
      W(3)=ONE
      W(6)=ZERO
!-----------------------------------------------------------------------
!     Output symmetry indices, number of panels
!-----------------------------------------------------------------------
      WRITE (IFILEN(3),909) ISX,ISY,NEQN
909   FORMAT (3I8)
!-----------------------------------------------------------------------
!     If ISX,ISY are negative reset them to zero for DEFMOD
!-----------------------------------------------------------------------
      ISX=MAX(ISX,0)
      ISY=MAX(ISY,0)
!-----------------------------------------------------------------------
!     Loop over all panels, read PRE input data, call DEFINE. NEWMDS is
!     assigned in DEFINE.  After the first call evaluate NDFR, allocate
!     arrays, and write number of modes and symmetry indices,
!-----------------------------------------------------------------------
      DO 80 N=1,NEQN
        READ(IFILEN(2),*) X,Y,Z,AREA,(VEL(J),J=1,6),IDIPOL
        CALL DEFINE(X,Y,Z,U,V,W,ISYM,NEWMDS)
        IF (N.EQ.1) THEN
          NDFR=6+NEWMDS
          IF (NDFR.GT.MDFR) GO TO 97
          WRITE(IFILEN(3),951) NEWMDS
          WRITE(IFILEN(3),951) (ISYM(K),K=7,NDFR)
        ENDIF
!-----------------------------------------------------------------------
!     Evaluate new modes (VEL) and write to output file
!-----------------------------------------------------------------------
        DO K=7,NDFR
          VEL(K)=U(K)*VEL(1)+V(K)*VEL(2)+W(K)*VEL(3)
        END DO
        WRITE(IFILEN(3),952) (VEL(K),K=7,NDFR)
!-----------------------------------------------------------------------
!     Accumulate hydrostatic coefficients
!-----------------------------------------------------------------------
        W(4)=Y
        W(5)=-X
        IF(IDIPOL.NE.0) CYCLE
        DO 70 I=1,NDFR
          DO J=1,NDFR
            C(I,J)=C(I,J)+                                              &
     &             AREA*VEL(J)*W(I)*MIJXY(ISYM(I),ISYM(J),ISX,ISY)
          END DO
70      CONTINUE
80    CONTINUE
!-----------------------------------------------------------------------
!     Write hydrostatic coefficients to output file and monitor
!-----------------------------------------------------------------------
      DO I=1,NDFR
        WRITE(IFILEN(3),952) (C(I,J),J=1,NDFR)
        WRITE(INTOP,952) (C(I,J),J=1,NDFR)
      END DO
 951  FORMAT(10I6)
 952  FORMAT(1P,6E13.5)
      WRITE (INTOP,903) ' file gdf.MOD completed; now re-run POTEN '
      STOP
  97  WRITE (INTOP,903)                                                 &
     &  ' ERROR: too many NEWMDS in DEFMOD compared to POTEN input '
      STOP
  98  WRITE (INTOP,903) ' Error -- cannot read input file gdf.PRE '
      STOP
      END

      SUBROUTINE DEFINE(X,Y,Z,U,V,W,ISYM,NEWMDS)
C-----------------------------------------------------------------------
      REAL U(*),V(*),W(*)
      INTEGER ISYM(*)
      DATA ZERO/ 0.0 /
      DATA PI/ 3.141592653589793 /
      REAL wn(1:36),nrm(1:36),L,R,Di,Ts,rho,zs,Ss,lk,uk,S,dr

C-----------------------------------------------------------------------
C    FlexWEC Tube, Following Aurelien Babarit, 2017, Journal of
C    Journal of Fluids and Structures, analytic solution
C-----------------------------------------------------------------------

      ! Read in .txt file of design variables
      OPEN(unit = 15, file = 'defmod.txt', status = 'old')
      READ(15,*) NEWMDS
      READ(15,*) L
      READ(15,*) R
      READ(15,*) Di
      READ(15,*) Ts
      READ(15,*) rho
      READ(15,*) zs
      READ(15,*) (wn(i), i = 1, NEWMDS)    
      READ(15,*) (nrm(i), i = 1, NEWMDS)
      CLOSE(15)
      
      ! Calculate tube area
      Ss=PI*R*R
      DO J=7,NEWMDS+6
        N=J-6
        
        ! Constants for each mode shape
        lk=SQRT((2*PI/(Di*Ts))*(SQRT(1+Ts*rho*Di*Di*wn(N)*wn(N)/PI)-1))
        uk=SQRT((2*PI/(Di*Ts))*(SQRT(1+Ts*rho*Di*Di*wn(N)*wn(N)/PI)+1))
        
        ! Zero modal displacement along the x axis U
        U(J)=ZERO
        ISYM(J)=5
        
        ! First set of mode shapes
        IF (N.LE.(NEWMDS/2)) THEN
          S=Ss-Ss*nrm(N)*(lk*TANH(uk*L/2)*COS(lk*X)/COS(lk*L/2)-        &
     &      uk*TAN(lk*L/2)*COSH(uk*X)/COSH(uk*L/2))
          
        ! Second set of mode shapes
          ELSE
          S=Ss-Ss*nrm(N)*lk*uk*(-TANH(uk*L/2)*SIN(lk*X)/COS(lk*L/2)+    &
     &      TAN(lk*L/2)*SINH(uk*X)/COSH(uk*L/2))
        ENDIF 
        
        ! Calculate radial displacement in the tube and then break into components V and W
        dr=SQRT(S/PI)-R
        V(J)=(Y/R)*dr
        W(J)=((Z-zs)/R)*dr   
      END DO
   
      RETURN
      END

      SUBROUTINE RDNAMD(IFILEN,INTIP,INTOP,FILEN,FNAMES,EXTN)
!-----------------------------------------------------------------------
!
!     COPYRIGHT (C) 1992 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
!
!-----------------------------------------------------------------------
!
!     DESCRIPTION : Input filename of GDF file from the file
!                   FNAMES.WAM or from interactive input
!
!-----------------------------------------------------------------------
!
!     Input Arguments (in alphabetical order):
!
!     Symbolic  Where value  Description
!       name    is assigned
!     ------------------------------------------------------------------
!     IFILEN    DEFMOD      Device numbers for files
!     INTIP     DEFMOD      Device number for keyboard
!     INTOP     DEFMOD      Device number for monitor
!     FNAMES    DEFMOD      Name of input file FNAME.WAM
!     EXTN      DEFMOD      Extensions for file names
!
!     Output Arguments (in alphabetical order):
!
!     Symbolic    Description
!       name
!     ------------------------------------------------------------------
!     FILEN     Names of input files
!-----------------------------------------------------------------------
      CHARACTER*1 YN
      CHARACTER*4 EXTN(*)
      CHARACTER*20 FILEN(*),FNAMES
      INTEGER IFILEN(*)
      
      ! EDIT: Lines here are commented out in order to no longer interact
      ! with user. defmod now always reads file name from FNAMES.WAM 
      ! instead of interactively entering in .gdf file name or having the
      ! user type in 'Y' or 'y' to select reading from FNAMES.WAM
      
      !WRITE(INTOP,907) 
!-----------------------------------------------------------------------
!     Enter GDF interactively or read from FNAMES.WAM
!-----------------------------------------------------------------------
      !READ(INTIP,900) YN
      !IF(INDEX(YN,'Y')+INDEX(YN,'y').EQ.0) THEN
      !  WRITE(INTOP,908)
      !  READ (INTIP,900) FILEN(1)
      !  FILEN(1)=ADJUSTL(FILEN(1))
      !  GOTO 25
      !ENDIF
      OPEN (IFILEN(1), FILE=FNAMES, STATUS='OLD', ERR=11)
 10   READ (IFILEN(1),900,END=11) FILEN(1)
      DO I=20,1,-1
        IF(FILEN(1)(I:I).NE.' ') EXIT
      ENDDO
      IF((FILEN(1)(I-2:I-2)=='G'.OR.FILEN(1)(I-2:I-2)=='g').AND.
     &   (FILEN(1)(I-1:I-1)=='D'.OR.FILEN(1)(I-1:I-1)=='d').AND.
     &   (FILEN(1)(I:I)=='F'.OR.FILEN(1)(I:I)=='f')) THEN
        GOTO 20
      ENDIF
      GOTO 10
 900  FORMAT (A)
!-----------------------------------------------------------------------
!     Enter GDF 
!-----------------------------------------------------------------------
 11   WRITE (INTOP,903) FNAMES
 903  FORMAT (1X,/,1X,A20,'file cannot be opened or is incomplete; ',// &
     & ' Enter full name of geometric data file (GDF):          ',$)
      READ (INTIP,900) FILEN(1)
!-----------------------------------------------------------------------
!     Remove leading blanks from filename
!-----------------------------------------------------------------------
      FILEN(1)=ADJUSTL(FILEN(1))
  20  CLOSE (IFILEN(1))
  25  CONTINUE
!-----------------------------------------------------------------------
!     Derive filename for PRE input file from PREMOD to DEFMOD
!     (full name = GDF filename with PRE extension)
!-----------------------------------------------------------------------
      CALL NEWNAM(FILEN,EXTN,1,2)
!-----------------------------------------------------------------------
!     Derive filename for MOD output file from DEFMOD to POTEN
!     (full name = GDF filename with MOD extension)
!-----------------------------------------------------------------------
      CALL NEWNAM(FILEN,EXTN,2,3)
!-----------------------------------------------------------------------
!     Check to see if an old output file exists with the same name
!     (IFILEN(1) is used temporarily for this unit number)
!-----------------------------------------------------------------------
      OPEN (IFILEN(1), FILE=FILEN(3), STATUS='OLD', ERR=50)
 30   WRITE (INTOP,906) FILEN(3)
906   FORMAT(1X,/                                                       &
     & 4X,' WARNING!'//                                                 &
     & 4X,' An old MOD file exists with the specified filename ',A,/    &
     & 4X,' Options: overwrite old file or abort this run?:'//          &
     & 4X,' Overwrite old file?  (Y/N) :',$)
      READ (INTIP,900) YN
      IF(INDEX(YN,'Y')+INDEX(YN,'y').EQ.0) STOP
!-----------------------------------------------------------------------
!     Delete old output file
!-----------------------------------------------------------------------
      CLOSE (IFILEN(1), STATUS='DELETE')
907   FORMAT(1X,/,
     &' When NBODY=1 and IALTPOT=1, GDF filename can be read from FNAME
     &S.WAM.'//                                                         &
     &' Read GDF filename from FNAMES.WAM? (Y/N) :   ',$)
908   FORMAT (/,'Enter full name of geometric data file (GDF):    ',$)
50    RETURN
      END

      SUBROUTINE NEWNAM(FILEN,EXTN,IOLD,INEW)
!-----------------------------------------------------------------------
!
!     COPYRIGHT (C) 1988 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
!
!-----------------------------------------------------------------------
!
!     DESCRIPTION : Derives new file name with extension from old
!     filename
!
!-----------------------------------------------------------------------
!
!     Character arrays are as follows
!        EXT(N) = filename extensions defined in PREMODDATA statement
!      FILEN(N) = filenames defined in PREMODDATA statement and input:
!
!                       PREMODMODULE            FORCE MODULE
!                       ----------------------------------------------
!      FILEN(1)         GDF filename            OUT filename
!      FILEN(2)         POT control filename    FRC control filename
!      FILEN(3)         P2F intermediate file   P2F input filename
!
!     IOLD = index of old filename
!     INEW = index of new filename and extension
!
!     Maximum length of filename without extension is 16 characters
!     if input filename is blank, output filename=new extension only
!
!-----------------------------------------------------------------------
      CHARACTER*4 EXTN(*)
      CHARACTER*20 FILEN(*)
!-----------------------------------------------------------------------
!     Assign new filename and strip leading blanks
!-----------------------------------------------------------------------
      FILEN(INEW)=FILEN(IOLD)
      DO 10 N=1,20
        IF (INDEX(FILEN(INEW),' ').EQ.1) FILEN(INEW)=FILEN(INEW)(2:20)
  10  CONTINUE
      IENDP=INDEX(FILEN(INEW),'.')
      IF (IENDP.GT.0) THEN
        IEND=IENDP
      ELSE
        IEND=INDEX(FILEN(INEW),' ')
      ENDIF
      FILEN(INEW)(IEND:IEND+4)=EXTN(INEW)
      RETURN
      END

      SUBROUTINE DATTIM(ADATE,ATIME)
!-----------------------------------------------------------------------
!
!     COPYRIGHT (C) 2000 WAMIT Incorporated
!     COPYRIGHT (C) 1998 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
!
!-----------------------------------------------------------------------
!
!     Source code file: DATTIM.FOR
!
!-----------------------------------------------------------------------
!
!     Version : 6.0
!
!-----------------------------------------------------------------------
!
!     DESCRIPTION :  Sets date and time from system clock
!                    Returns each as a character string
!                      (dimensioned in main program)
!
!    The code in this subroutine is generic Fortran-90 and should not
!    require any system-dependant modification.
!
!-----------------------------------------------------------------------
      CHARACTER*(*) ADATE,ATIME
      CHARACTER*8 DATE
      CHARACTER*10 TIME
      CHARACTER*3 MON(12)
      DATA MON/ 'Jan','Feb','Mar','Apr','May','Jun',                    &
     &          'Jul','Aug','Sep','Oct','Nov','Dec'/
      CALL DATE_AND_TIME(DATE,TIME)
      M=10*(ICHAR(DATE(5:5))-48)+ICHAR(DATE(6:6))-48
      ADATE=DATE(7:8)//'-'//MON(M)//'-'//DATE(1:4)
      ATIME=TIME(1:2)//':'//TIME(3:4)//':'//TIME(5:6)
      RETURN
      END


      SUBROUTINE HEADRD(IOUT,VERHDR)
!-----------------------------------------------------------------------
!
!     COPYRIGHT (C) 2000 WAMIT Incorporated
!     COPYRIGHT (C) 1998 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
!
!-----------------------------------------------------------------------
!
!     Source-code file :  DEFMOD.FOR
!
!-----------------------------------------------------------------------
!
!     VERSION : 6.0
!
!-----------------------------------------------------------------------
!
!     DESCRIPTION : Writes copyright statement and other information
!
!-----------------------------------------------------------------------
!
!     Input Arguments (in alphabetical order):
!
!     Symbolic  Where value  Description
!       name    is assigned
!     ------------------------------------------------------------------
!     IOUT      POTEN        Device number screen
!     VERHDR    POTEN        Version number
!
!     MOD: ADD WAMIT COPYRIGHT
!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) VERHDR
      CHARACTER*16 WM
      CHARACTER*66 COP,COPOLD
      INTEGER IOUT
!-----------------------------------------------------------------------
!     Output WAMIT header and Copyright statement on the screen
!-----------------------------------------------------------------------
      WM  =' WAMIT  Version '
      COP =                                                             &
     &'Copyright (c) 2000 WAMIT Incorporated'
      COPOLD =                                                          &
     &'Copyright (c) 1998 Massachusetts Institute of Technology'
      WRITE(IOUT,200) WM,VERHDR,COP,COPOLD
200   FORMAT(1X,71(1H-),//,23X,2A,//,5X,A,/,5X,A,                       &
     &   //1X,71(1H-),//                                                &
     & 5X,'The WAMIT software performs computations of wave interactions&
     & with'/                                                           &
     & 5X,'floating or submerged vessels.  WAMIT is a registered tradema&
     &rk of'/                                                           &
     & 5X,'WAMIT Incorporated.  Further information is available at '// &
     & 5X,'                  www.wamit.com                          '// &
     & 5X,'DEFMOD is a user-modified program to define generalized modes&
     &.'/                                                               &
     & 5X,'DEFMOD is distributed in source-code to all licensed users of&
     & the'/                                                            &
     & 5X,'WAMIT software.'/                                            &
     & /,1X,71(1H-),//)
      RETURN
      END
