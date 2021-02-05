      PROGRAM SCHWDRIVER
C   THIS PROGRAM DRIVES DSCPACK. IT READS FROM THE FILE DSCDATA. ON
C   OUTPUT,IT GIVES THE COMPUTED ACCESSORY PARAMETERS.THE ACCURACY
C   OF THE PARAMETERS IS ALSO TESTED.
C
C
C     .. Scalars in Common ..
      DOUBLE PRECISION DLAM
      INTEGER IU
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION UARY(8),VARY(3)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX C,WW,ZZ,ZZ0,WP,W,E,ZZINV,INVPREF
      DOUBLE PRECISION EPS,PI,TOL,U
      INTEGER I,IGUESS,INV,IPOLY,ISHAPE,ISOLV,ITRY,K,LINEARC,M,N,NPTQ
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX W0(30),W1(30),Z0(30),Z1(30)
      DOUBLE PRECISION ALFA0(30),ALFA1(30),PHI0(30),PHI1(30),QWORK(1660)
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX WDSC,ZDSC,WTHETA,WPROD,ARGUM
      EXTERNAL WDSC,ZDSC,WTHETA,WPROD,ARGUM
C     ..
C     .. External Subroutines ..
      EXTERNAL ANGLES,CHECK,DSCDATA,DSCSOLV,DSCTEST,QINIT,THDATA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ACOS
C     ..
C     .. Common blocks ..
      COMMON /PARAM4/UARY,VARY,DLAM,IU
C     ..
      PI = ACOS(-1.D0)
C
C  INITIALISING ARRAYS:
      DO 10 K = 1,30
          ALFA0(K) = 0.D0
          ALFA1(K) = 0.D0
   10 CONTINUE

C  INITIALISING SHAPE
      M = 3
      N = 2
      Z0(1) = (2.0D0,0.D0)
      Z0(2) = (100.D0,00.D0)
      Z0(3) = (1.00001D0,0.D0)
      Z1(1) = (0.D0,0.D0)
      Z1(2) = CMPLX(1/1.00001D0,0.D0)
      ALFA0(1) = 1.D0
      ALFA0(2) = -2.D0
      ALFA0(3) = 2.D0
      CALL ANGLES(N,Z1,ALFA1,1)
   20 CONTINUE
C
C  GENERATE THE GAUSS-JACOBI WEIGHTS & NODES AND CHECK THE INPUT:
C     WRITE (6,FMT=*) 'INPUT NPTQ(THE NUMBER OF G-J POINTS)'
C     WRITE (6,FMT=*)
C    +  '(RECOMMENDED VALUES FOR NPTQ ARE: 2,3,4,5,6,7,OR 8)'
C     READ (5,FMT=*) NPTQ
      NPTQ = 8
      CALL QINIT(M,N,ALFA0,ALFA1,NPTQ,QWORK)
      ISHAPE = 1
C     IF (IPOLY.EQ.2 .OR. IPOLY.EQ.5 .OR. IPOLY.EQ.7) ISHAPE = 1
      CALL CHECK(ALFA0,ALFA1,M,N,ISHAPE)
C
C  SPECIFY SOME PARAMETERS OF THE CALLING SEQUENCE OF DSCSOLV:
      IGUESS = 1
      LINEARC = 1
      TOL = 1.D-10
C
C   SOLVE THE ACCESSORY PARAMETER PROBLEM:
      ISOLV = 1
      IF (ISOLV.EQ.1) CALL DSCSOLV(TOL,IGUESS,M,N,U,C,W0,W1,PHI0,PHI1,
     +                             Z0,Z1,ALFA0,ALFA1,NPTQ,QWORK,ISHAPE,
     +                             LINEARC)
C
C   OUTPUT WILL BE ARRANGED IN DSCSOLV.
C
C   COMPUTE DATA FOR THETA-FUNCTION AND TEST THE ACCURACY:
      CALL THDATA(U)
      CALL DSCTEST(M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,NPTQ,QWORK)

C     ABS OF PRE-IMAGE OF UNIT DISC
      ZZ = (1.D0,0.D0)
      EPS = 1.D-6
      WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +        NPTQ,QWORK,EPS,1)

      OPEN(11,FILE="output/sprime.txt")
      OPEN(12,FILE="output/wprods.txt")

      NPT = 2000
      WIDTH = 2
      DO 30 I = 0,NPT-1
        ZARG = -.5*WIDTH*PI + WIDTH*PI*I/NPT 
        ZZ = CMPLX(COS(ZARG),SIN(ZARG))
        EPS = 1.D-6
        WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +          NPTQ,QWORK,EPS,1)

        ARG = ATAN(IMAG(WW),REAL(WW))
        WRITE(14,FMT=*) arg

        WP = WPROD(WW,M,N,U,W0,W1,ALFA0,ALFA1)
        ZP0 = WP

        ZZINV = CMPLX(1.0D0,0.0D0)/(ZZ)
        INVPREF = 0.5*WW*ZZINV*(C*WP)

        WRITE(12,*) REAL(WP),"+I*",IMAG(WP)
        WRITE(11,*) REAL(INVPREF),"+I*",IMAG(INVPREF)

   30 CONTINUE
      ENDFILE(11)
      CLOSE(11)
      ENDFILE(12)
      CLOSE(12)

      END
