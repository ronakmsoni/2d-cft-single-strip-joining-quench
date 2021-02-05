C
C   THIS IS THE FINAL VERSION SUBMITTED TO ACM TOMS.
      PROGRAM TESTDRIVER
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
      DOUBLE COMPLEX C,WW,ZZ,ZZ0,WP,W
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
      Z0(2) = (100.D0,100.D0)
      Z0(3) = (1.00001D0,0.D0)
      Z1(1) = (0.D0,0.D0)
      Z1(2) = (.99999D0,0.D0)
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
C
C  FOLLOWING PARAGRAPH IS FOR TESTING INVERSE EVALUATIONS:
      WRITE (6,FMT=*)
     +  'WANT TO INVERSE THE MAP? (1 FOR "YES",2 FOR "NO")'
      READ (5,FMT=*) INV
      IF (INV.EQ.1) THEN
          DO 40 I = 1,10
              WRITE (6,FMT=*)
     +          'INPUT THE POINT INSIDE THE DOMAIN IN THE FORM OF'
              WRITE (6,FMT=*) '          ( X-COORD., Y-COORD. )'
              READ (5,FMT=*) ZZ
              EPS = 1.D-6
              WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +             NPTQ,QWORK,EPS,1)
              WRITE (6,FMT=*) 'THE PREIMAGE OF ZZ=',WW,
     +            " =",ABS(WW),"*E*I",ATAN(IMAG(WW),REAL(WW))
              IF (ABS(WW).LE.1.D-12) GO TO 30
              WRITE (6,FMT=*) 'CHECK BY MAPPING THE PREIMAGE BACK'
              ZZ0 = ZDSC(WW,0,2,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,
     +              PHI1,NPTQ,QWORK,1)
              WRITE (6,FMT=*) 'THE POINT ENTERED=',ZZ0
   30         WRITE (6,FMT=*)
     +          'WANT TO TRY ANOTHER PT ? (1 FOR "YES",2 FOR "NO")'
              READ (5,FMT=*) INV
              IF (INV.EQ.2) GO TO 50
   40     CONTINUE
      END IF
C
   50 WRITE (6,FMT=*)
     +  'WANT TO EVALUATE THE MAP? (1 FOR "YES",2 FOR "NO")'
      READ (5,FMT=*) EVAL
      IF (EVAL.EQ.1) THEN
          DO 70 I = 1,10
              WRITE (6,FMT=*)
     +          'INPUT THE POINT INSIDE THE ANNULUS IN THE FORM OF'
              WRITE (6,FMT=*) '          ( X-COORD., Y-COORD. )'
              READ (5,FMT=*) WW
              EPS = 1.D-6
              ZZ = ZDSC(WW,0,2,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,
     +              PHI1,NPTQ,QWORK,1)
              WRITE (6,FMT=*) 'THE IMAGE=',ZZ
   60         WRITE (6,FMT=*)
     +          'WANT TO TRY ANOTHER PT ? (1 FOR "YES",2 FOR "NO")'
              READ (5,FMT=*) EVAL
              IF (EVAL.EQ.2) GO TO 71
   70     CONTINUE
      END IF
C     PRINT '(/)'
C  60 WRITE (6,FMT=*)
C    +  'WANT TO TRY ANOTHER EXAMPLE? (1 FOR "YES",2 FOR "NO")'
C     READ (5,FMT=*) ITRY
C     IF (ITRY.EQ.1) GO TO 10

C     Checking inversion symmetry.
C     A = (1,0)
C     WRITE (6,FMT=*) REAL(A)

   71 NPT = 20
      ZZ = (1.D0,0.D0)
      EPS = 1.D-6
      WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +        NPTQ,QWORK,EPS,1)
      WA = ABS(WW)
      PHI = ATAN(IMAG(WW),REAL(WW))
      WRITE(6,FMT=*) "W1=",WW
      WRITE(6,FMT=*) "ABS1=",WA
      WRITE(6,FMT=*) "ARG1=",PHI

      DO 80 I = 0,NPT
        WW = WA*CMPLX(COS(PI*I/NPT),SIN(PI*I/NPT))
        ZZ = ZDSC(WW,0,2,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,
     +              PHI1,NPTQ,QWORK,1)
        ZA2 = ABS(ZZ)
        
        ZZ = CMPLX(COS(PI*I/NPT),SIN(PI*I/NPT))
        EPS = 1.D-6
        WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +          NPTQ,QWORK,EPS,1)
        WRITE (6,FMT=*) I," ",ZA2," ",ABS(WW)
        IF (I.EQ.N) GO TO 80
   80 CONTINUE



      OPEN(13,FILE="output/real-wprod.txt")
      OPEN(14,FILE="output/real-prefactor.txt")

      NPT = 100
      DO 90 I = 0,NPT
        ZZ = Z1(2) + I*(Z0(3) - Z1(2))/NPT
        EPS = 1.D-6
        WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +          NPTQ,QWORK,EPS,1)

        WP = WPROD(WW,M,N,U,W0,W1,ALFA0,ALFA1)
        WRITE(13,*) WP
        WRITE(14,*) 2.0*ZZ/(C*WP)
   90 CONTINUE
      ENDFILE(13)
      CLOSE(13)
      ENDFILE(14)
      CLOSE(14)

      DEL = 1.D-5
      WRITE(6,FMT=*) DEL
      WW = CMPLX(WA,0)
C     ZZ = ZDSC(WW,0,2,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,
C    +            PHI1,NPTQ,QWORK,1)
      WP = WPROD(WW,M,N,U,W0,W1,ALFA0,ALFA1)
      ZP0 = WP
      WW = (1-DEL)*WW
      WP = WPROD(WW,M,N,U,W0,W1,ALFA0,ALFA1)
      ZP1 = WP

      ZPP0 = (ZP0 - ZP1)/DEL
      PS0 = ZPP0/ZP0
      WRITE(6,FMT=*) "PS0 =",PS0

      
      WW = (1.0 - 2.0*DEL)*WW/(1-DEL)
      WP = WPROD(WW,M,N,U,W0,W1,ALFA0,ALFA1)
      ZP2 = WP

      ZPP1 = (ZP1 - ZP2)/DEL
      PS1 = ZPP1/ZP1
      WRITE(6,FMT=*) "PS1 =",PS1
      
      SCH0 = (PS0 - PS1)/DEL - .5*PS0*PS0
      WRITE(6,FMT=*) "{z,w}( z = 1 ) = ",SCH0

      NPT = 2000
      ZZ = (1.D0,0.D0)
      EPS = 1.D-6
      WIDTH = 2
      WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +        NPTQ,QWORK,EPS,1)
      WA = ABS(WW)
      PHI = ATAN(IMAG(WW),REAL(WW))

      OPEN(15,FILE="output/args.txt")
      DO 100 I = 0,NPT
        ZARG = WIDTH*PI*I/NPT
        ZZ = CMPLX(COS(ZARG),SIN(ZARG))
        EPS = 1.D-6
        WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +          NPTQ,QWORK,EPS,1)
        PHI = ATAN(IMAG(WW),REAL(WW))
        WRITE(15,*) PHI
        IF (I.EQ.N) GO TO 100
  100 CONTINUE
C   90 STOP

      END
