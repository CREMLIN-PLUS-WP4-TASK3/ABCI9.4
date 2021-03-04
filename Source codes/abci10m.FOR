      PROGRAM ABCI
C-----------------------------------------------------------------------
C
C              PROGRAM      A B C I
C
C    AZIMUTHAL BEAM CAVITY INTERACTION IN A CYLINDRICALLY
C                 SYMMETRIC STRUCTURE
C
C           YONG HO CHIN ,    CERN/SL-AP, FEBRUARY 1994
C
C    VERSION 1.0:  CERN LEP/TH , JAN.1988
C    VERSION 2.0:  ADDITION OF OPTION TO CALCULATE M=0 MODE
C    VERSION 2.1:  CORRECTION OF AZIMUTHAL WAKE CALCULATION
C    VERSION 3.0:  MOVING MESH INSTALLED,   MAY 1992
C    VERSION 3.1:  NEW MESH GENERATOR INSTALLED, AUG.1992
C    VERSION 3.2:  CIRCLULAR AND ELLIPTIC BOUNDARY INPUT, AUG.1992
C    VERSION 4.0:  JUMPING MESH, SEPT.1992
C    VERSION 4.1:  VARIABLE RADIAL MESH SIZES
C    VERSION 4.2:  INCREMENTAL MESH INPUT, 24 SEPT.1992
C    VERSION 4.3:  REPETITION COMMAND FOR MESH, 28 SEPT.1992
C    VERSION 5.0:  NAPOLY'S INTEGRATION LONGITUDINAL
C    VERSION 5.1:  NAPOLY'S INTEGRATION TRANSVERSE, 29 SEPT.1992
C    VERSION 6.0:  PLOTTING WITH TOPDRAWER, 1 OCT.1992
C    VERSION 6.1:  FFT OF WAKE POTENTIALS,  5 OCT.1992
C    VERSION 6.1.1:CORRECTION OF NAPOLY LONGI. W/ LOG TERM, 23 OCT.1992
C    VERSION 6.2:  NEW NAPOLY'S INTEGRATION TRANSVERSE, 23 OCT.1992
C    VERSION 6.2.1:CORRECTION OF ELLIPTIC CONNE. AND CEPS, 30 NOV.1992
C    VERSION 6.2.2:CORRECTION OF ESTIMATE OF NW, 2 DEC.1992
C    VERSION 6.2.3:CORRECTION OF VARIABLE RADIAL MESH SIZES, 6 JAN.1993
C    VERSION 6.2.4:CORRECTION OF TRANS. MAX/MIN IN PLOT, 21 FEB.1993
C    VERSION 6.3:  ARBITRARY BUNCH SHAPE, 15 MAR.1993
C    VERSION 6.4:  NEW FFT (IMPEDANCE) AND ENERGY SPECTRUM, 26 MAR.1993
C    VERSION 6.5:  CHANGE ON CHOICE OF NAPOLY-ZOTTER PATH, 31 MAR.1993
C    VERSION 6.6:  PLOT OF ELECTRIC FIELD (CURRENT) LINES, 14 APR.1993
C    VERSION 6.7:  SAVE WAKE POTENTIALS IN A FILE ON FT08, 12 JUL.1993
C    VERSION 6.8:  SAVE FFT RESULTS IN A FILE ON FT08, 13 JUL.1993
C    VERSION 7.0:  WINDOWING OF WAKE POTENTIALS FOR FFT, 1 AUG.1993
C    VERSION 7.0.1:WINDOW FUNCTION --> BLACKMAN-HARRIS, 5 AUG.1993
C    VERSION 7.1:  ISLANDS INPUT, 27 AUG.1993
C    VERSION 8.0:  OVERHAUL FOR SPPED-UP BY H. GROTE(CERN), 30 AUG.1993
C    VERSION 8.1:  OVERHAUL OF WAKE INTEGRATION ROUTINES, 3 SEP.1993
C    VERSION 8.1.1:SPLIT MAIN PROGRAM FOR OPTIMIZATION, 10 SEP.1993
C    VERSION 8.2:  THREE CHOICES IN WINDOW FUNCTION, 15 SEP.1993
C    VERSION 8.2.1:CORRECTION OF CHARACTER SET FOR CRAY, 11 NOV.1993
C    VERSION 8.3:  SIMPLER FIT EQS FOR MROT=1 NEAR AXIS, 21 DEC.1993
C    VERSION 8.3.1:CORRECTION OF AZIMUTHAL WAKE CALCULATION, 5 JAN.1994
C    VERSION 8.4:  MAXIMUM MT SPECIFIABLE BY ARRAY SIZE NPT, 5 JAN.1994
C    VERSION 8.4.1:CORRECTION IN SUBROUTINE EXPAND2, 20 JAN.1994
C    VERSION 8.4.2:CHANGE IN TRANSVERSE WAKE CALCULATION, 24 JAN.1994
C    VERSION 8.4.3:MORE CHOICES IN WAKE POTENTIAL PLOTS, 24 JAN.1994
C    VERSION 8.4.4:IMMEDIATE STOP AFTER AN OVERFLOW, 24 JAN.1994
C    VERSION 8.5  :CRITICAL MT FOR STABILITY IF DDR/DDZ<1, 24 JAN.1994
C    VERSION 8.6  :NEW LAYOUT IN WAKE AND IMPEDANCE PLOTS, 2 FEB.1994
C    VERSION 8.6.1:FURTHER CHANGE IN LAYOUT, 9 FEB.1994
C    VERSION 8.7  :CHANGE IN UNITS OF WAKE POTENTIALS, 10 FEB.1994
C    VERSION 8.7.1:CORRECTIONS IN BSHAPE, WPADD, & WPLOT, 24 FEB.1994
C    VERSION 8.8  :PLOT INTEGRATION OF IMPEDANCE OVER FREQ, 28 FEB.1994
C    VERSION 8.8.1:CORRECTION IN ISIG, 5 MAY.1994
C    VERSION 8.9  :CORRECTION IN NOMALIZATION OF IMPEDANCE, 20 MAY.1994
C    VERSION 9.0  :WAKE POTENTIALS FOR COUNTER-ROTATING BEAM 9 JUN.1994
C    VERSION 9.1  :ALL PLOTS' DIMENSIONS REDUCABLE, 10 JUN.1994
C    VERSION 9.2  :CORRECTION IN NAPOLY-ZOTTER CONTOUR, 17 NOV.1994
C    VERSION 9.2.1:PRINT LOSS FACTORS FOR FIRST/LAST HALF, 5 JUL.1995
C    VERSION 9.2.2:CORRECTION IN SAVE LONG. WAKE WITH LOG, 21 JAN.1997
C    VERSION 9.3  :MROT=1 COMP.FOR UNEQUAL BEAM PIPE RADII, 21 JAN.1997
C    VERSION 9.4  :OPTION TO SET EXPFAC AND CUTOFF FREQ. 21 JAN.1997
C    VERSION XP   :WINDOWS VERSION OF ABCI 9.4 4 JUN.2005   
C
C-----------------------------------------------------------------------
C    *** NOTE: How to increase the size of mesh points ***
C
C    VM/CMS COMMAND PROCEDURE...
C
C    XEDIT ABCI FORTRAN
C    CHA /NMESH=10000000/NMESH=VALUE YOU LIKE/ * *
C    CHA /NRN=20000/NRN=###/ * *
C    CHA /NWN=500000/NWN=###/ * *
C    CHA /NHIT=3000000/NHIT=###/ * *
C    CHA /NBP=1000/NBP=###/ * *
C    CHA /NF=100/NF=###/ * *
C    CHA /NSP=20000/NSP=###/ * *
C    CHA /NVR=100/NVR=###/ * *
C    CHA /NPF=624/NPF=###/ * *
C    CHA /NPT=20/NPT=###/ * *
C    SAVE
C
C    *** COMMENT ***
C
C    NRN=MAX(NR)
C    NWN=MAX(NW)+2
C    NMESH=MAX(NP)
C    NHIT=A FEW TIMES NMESH. AT LEAST NRN*3
C    NBP=MAXIMUM OF THE NUMBER OF INPUT LINES FOR CAVITY STRUCTURE
C    NSP=MAXIMUM OF THE NUMBER OF FITTING POINTS FOR CAVITY STRUCTURE
C    NVR=MAXIMUM NUMBER OF STRUCTURE VARIABLES
C    NF=MAXIMUM NUMBER OF DIFFERENT RADIAL MESH SIZES
C    NPF=MAXIMUM NUMBER OF POINT IN FREQUENCY IN FOURIER TRANSFORMS
C    NPT=MAXIMUM NUMBER OF MT+1
C   WHERE
C    NW=MAX(UBT/DDZ, 10*SIG/DDZ)=THE NUMBER OF WAKE POTENTIAL POINTS
C    NP=MIN(NR*(NW+4),NR*NZ)=THE NUMBER OF MESH POINTS IN THE WINDOW
C-----------------------------------------------------------------------
C    *** NOTE: How to install the ABCI source code into your computer
C
C    VAX users who don't use the DATLIB library must locate
C    "FOR VAX USERS" comments, and comment out "CHARACTER*8 NDATE,
C    NTIME" lines and include lines containing "REAL*8 NDATE, NTIME"
C    instead.
C
C    SUN and HP users must locate two "FOR SUN AND HP USERS" comments
C    in the subroutine PLTMSH and comment out lines which contain back
C    slash characters "\", and include lines which contain characters
C    "L" instead.
C
C    You must link ABCI with the CERN library to get the correct date,
C    time and the cpu time used. If you don't have one,
C    link with one of the DATLIB libraries provided by the author.
C
C    HP users must use the "-k" option during compilation.
C    They must also link ABCI with DATLIB and libU77.
C
C    Users of some computers may have to change the logical device
C    numbers.  This can be done by editing the routine BLOCKDATA
C    DEVNUM.
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NVR=100)
      PARAMETER (NF=100)
      PARAMETER (NPT=20)
      LOGICAL LSAV,LREC,LCFRON,SYM,LCBACK,LPLW,LPRW,LCHIN,LMATPR
      LOGICAL LMESHD,LNAPOLY,LBEAMD,LWAKED,LCAVIN,LCAVUS
      LOGICAL LADD,LPOT,LNTR,LNAPOLS,LSOS,LPLC,LELC
      LOGICAL LFFTL,LFFTA,LFFTT,LFFT,LCPUTM,LSPEC,LPLE
      LOGICAL LPALL,LFFTOR,LSVW,LSVF,LWNDW,LPATH,LNONAP
      LOGICAL LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,LINTZ,LCRBW
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
	CHARACTER*256 NMEA,NMEO,NMEP,NMET
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      CHARACTER*8 CAVT,CAVY,ADDT,ADDY,ADDH,ADDS
      CHARACTER*8 PRTNM,RDPRT,PRTNA
      CHARACTER*4 STOPC,CONTIN,BEAMD,BEAMA,DMESH,AMESH,WAKED,WAKEA
      CHARACTER*4 ENDFIL,BOUND,BOUNA,PLOTD,PLOTA,PRIND,PRINA
      CHARACTER*16 CAVHLF,CAVINP
      CHARACTER*1 RZIN(72,NBP)
      CHARACTER*2 ORD(3),ORDX
      CHARACTER*13 REDCHR
C-----------------------------------------------------------------------
      DIMENSION RMAK0(NF),DDRS(NF)
C-----------------------------------------------------------------------
C ASSIGN DEVICE NUMBER TO ND#. SEE BLOCKDATA DEVNUM
C-----------------------------------------------------------------------
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPL/ MROT
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /PLNZ/ LPATH
      COMMON /WPT/ NDATE,NTIME
      COMMON /WP2/ TT
      COMMON /OFFS/ RDRIVE,ZCF,ZCT,RWAK
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /CALS/ SUMQ,BSIG
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /SYMT/ SYM
      COMMON /ELINE/ NZF,NZL
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /VAR/ RZ(NVR)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /PIPE/ KOFF,KOFL
      COMMON /PINT/ KMIN
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
      COMMON /FFTS/ LSVF,LWNDW
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /ICON/ IZLFT,IZRGT,IRMIN
      COMMON /WFUN/ NWFUN,ALPHA
      COMMON /EXPF/ EXPFAC
      COMMON /INTR/ LINTZ
      COMMON /RED/ REDCHR
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      COMMON /WPS/ LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
C-----------------------------------------------------------------------
C INPUT VARIABLES
C-----------------------------------------------------------------------
      NAMELIST /FILE/ LSAV,ITEST,LREC,LCPUTM,TSOS,TMAX
      NAMELIST /BOUN/ IZL,IZR
      NAMELIST /MESH/ NR,NZ,DDR,DDZ,RAD,RZ,RMARK
      NAMELIST /BEAM/ MROT,SIG,ISIG,RDRIVE,NBUNCH,BSEP
      NAMELIST /TIME/ MT,TPS,NSHOT
      NAMELIST /WAKE/ UBT,LCFRON,LCBACK,LCHIN,LNAPOLY,ZCF,ZCT,RWAK,
     +LCRBW,ZSEP,LNONAP
      NAMELIST /PLOT/ LPLW,LCAVIN,LCAVUS,LFFTL,LFFTA,LFFTT,CUTOFF,LFFT,
     +LSPEC,LPLE,LPLC,NPLOT,LPALL,LWNDW,NWFUN,ALPHA,
     +LPLWL,LPLWA,LPLWT,LINTZ,REDUCE,LPATH,EXPFAC
      NAMELIST /PRINT/ LPRW,LMATPR,LSVW,LSVF,LSVWL,LSVWA,LSVWT
      NAMELIST /PRIN/ LPRW,LMATPR,LSVW,LSVF,LSVWL,LSVWA,LSVWT
C-----------------------------------------------------------------------
      DATA CAVT/'##CAVITY'/,CAVY/' ##CAVIT'/
      DATA ADDT/'@@CAVITY'/,ADDY/' @@CAVIT'/
      DATA ADDH/'@CAVITYS'/,ADDS/' @CAVITY'/
      DATA STOPC/'STOP'/,CONTIN/'CONT'/,BEAMD/' &BE'/,DMESH/' &ME'/
      DATA BEAMA/' $BE'/,AMESH/' $ME'/,BOUND/' &BO'/,BOUNA/' $BO'/
      DATA WAKED/' &WA'/,WAKEA/' $WA'/,PLOTD/' &PL'/,PLOTA/' $PL'/
      DATA PRIND/' &PR'/,PRINA/' $PR'/,PRTNM/' &PRINT '/
      DATA PRTNA/' $PRINT'/,ORD/'ND','RD','TH'/
      DATA CAVHLF/'HALF CELL INPUT '/
C-----------------------------------------------------------------------
C DEFAULT VALUES FOR PARAMETERS
C-----------------------------------------------------------------------
      CALL PRMINI(TMAX,TSOS,ITEST,UBT,CUTOFF,SYM,LADD,CAVINP,
     +LSAV,LREC,LCPUTM,LCFRON,LCBACK,LPRW,LPLW,LSVF,LSVW,LFFT,LFFTL,
     +LFFTA,LFFTT,LSPEC,LPLE,LPLC,LELC,LPALL,LWNDW,LCAVIN,LCAVUS,
     +LMATPR,LCHIN,LNAPOLY,LMESHD,LBEAMD,LWAKED,RDRIVE,ZCF,ZCT,RWAK,
     +TPS,NX1,NSHOT,ISIG,BSIG,NBUNCH,NPLOT,BSEP,NFST,RAD,ZL,
     +MT,MROT,NLOOP,NWFUN,ALPHA,LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,
     +LINTZ,LCRBW,ZSEP,REDUCE,LPATH,LNONAP,EXPFAC)
      CALL RZINIT
      DO 10 I=1,NF
   10 RMAK0(I)=0.0
      CALL IRDINI
C-----------------------------------------------------------------------
C OPEN STATEMENTS TO INITIALIZE INPUT AND OUTPUT FILES
C THE FOLLOWING IS AN EXAMPLE TO RUN ABCI FOR SAMPLE1 FILE ON WINDOWS
C-----------------------------------------------------------------------
C     OPEN(ND1, FORM='UNFORMATTED',FILE='sample1.sav',STATUS='UNKNOWN')
C     OPEN(ND2, FORM='UNFORMATTED',FILE='sample1.sa2',STATUS='UNKNOWN')
      OPEN(ND5, FILE=' ',STATUS='UNKNOWN')
	INQUIRE(UNIT=ND5, NAME=NMEA)
	CALL ASFILE(NMEA,NMEO,NMEP,NMET)
      OPEN(ND6, FILE=NMEO,STATUS='UNKNOWN')
      OPEN(ND8, FILE=NMEP,STATUS='UNKNOWN')
      OPEN(ND9, FILE=NMET,STATUS='UNKNOWN')
C-----------------------------------------------------------------------
C SET TIMER
C-----------------------------------------------------------------------
      CALL CPUTIM(TT,1)
C IBM/VM USER ONLY FOR IMMEDIATE TERMINATION AFTER ONE OVERFLOW
C     CALL ERRSET(207,1)
C-----------------------------------------------------------------------
C NOW READ INPUT DATA FROM FT05
C-----------------------------------------------------------------------
      READ(ND5,FILE)
      IF(.NOT.LREC) GO TO 30
      NBPN=(NBP-1)/400
      NMN=(NMESH-1)/8000
      READ(ND2) DDZ,DDR,RMARK,RZ,ZCAN,RCAN,SIG,ISIG,MT,
     +MROT,RDRIVE,NBUNCH,BSEP
      DO 20 I=1,NBPN+1
        ILAST=I*400
        IF(ILAST.GT.NBP) ILAST=NBP
   20 READ(ND2) ((RZIN(J,L),J=1,72),L=(I-1)*400+1,ILAST)
      READ(ND2) ZTXT,R,Z,NZP0,NLINE,NFST,UBT,LCFRON,LCBACK,
     +LCHIN,LNAPOLY
      READ(ND2) ZCF,ZCT,RWAK,LPLW,LCAVIN,LCAVUS,
     +LFFTL,LFFTA,LFFTT,LFFT,
     +LSPEC,LPLE,LPLC,NPLOT,TPS,NSHOT,NX1,TX1,TX2,LPALL,NVC0
      READ(ND2) CUTOFF,LPRW,LMATPR,SYM,LADD,NLOOP,LSVW,LSVF,LWNDW,
     +NWFUN,ALPHA,LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,LINTZ,LCRBW,ZSEP,
     +REDUCE,LPATH,LNONAP,EXPFAC
      LMESHD=.TRUE.
      NLOOP=NLOOP+1
      GO TO 40
C--TITLE----------------------------------------------------------------
   30 READ(ND5,10030) (ZTXT(JK),JK=1,9)
C-----------------------------------------------------------------------
C PRINT HEAD LINE
C-----------------------------------------------------------------------
      CALL DATIME(NDATE,NTIME)
      WRITE (ND6,10040) (ZTXT(II),II=1,9),NDATE,NTIME
C-----------------------------------------------------------------------
C READ THE REST OF INPUT DATA
C-----------------------------------------------------------------------
      READ(ND5,10000) RDPRT
      BACKSPACE(ND5)
      IF((RDPRT.EQ.BOUND).OR.(RDPRT.EQ.BOUNA)) THEN
        READ(ND5,BOUN)
      END IF
      READ(ND5,MESH)
      READ(ND5,10000) ZTXT(10)
      IF((ZTXT(10).EQ.CAVT).OR.(ZTXT(10).EQ.CAVY).OR.(ZTXT(10).EQ.ADDT)
     1.OR.(ZTXT(10).EQ.ADDY)) SYM=.TRUE.
      IF(SYM) CAVINP=CAVHLF
      IF((ZTXT(10).EQ.ADDT).OR.(ZTXT(10).EQ.ADDY).OR.(ZTXT(10).EQ.ADDH)
     1.OR.(ZTXT(10).EQ.ADDS)) LADD=.TRUE.
C--GEOMETRY OF STRUCTURE------------------------------------------------
      CALL INPUTI(LADD,NZP,NLINE,NVC,ZCAN,RCAN,RZIN)
      NZP0=NZP
      NVC0=NVC
      ZL00=Z(NVC)
      READ(ND5,BEAM)
      READ(ND5,TIME)
      READ(ND5,WAKE)
      READ(ND5,PLOT)
      READ(ND5,10000) RDPRT
      BACKSPACE(ND5)
      IF((RDPRT.EQ.PRTNM).OR.(RDPRT.EQ.PRTNA)) THEN
        READ(ND5,PRINT)
      ELSE
        READ(ND5,PRIN)
      END IF
   40 CONTINUE
C-----------------------------------------------------------------------
C SET UP PARAMETERS
C-----------------------------------------------------------------------
      DO 50 I=1,NF
        DDRS(I)=DDR(I)
   50 RMAK0(I)=RMARK(I)
      DDZS=DDZ
      NRS=NR
      NZS=NZ
      CUTOFFS=CUTOFF
      LNAPOLS=LNAPOLY
   60 CONTINUE
C TEXT
      IF((NLOOP.LE.2).AND.(NLOOP.GT.0)) ORDX=ORD(NLOOP)
      IF(NLOOP.GT.2) ORDX=ORD(3)
      IF(NLOOP.LE.0) GO TO 70
      IF(.NOT.LBEAMD) THEN
        CALL DATIME(NDATE,NTIME)
        WRITE (ND6,10050) NLOOP+1,ORDX,NDATE,NTIME
      END IF
      NZP=NZP0
      NVC=NVC0
      CALL INPUTD(LADD,NZP,NLINE,NVC,ZCAN,RCAN,RZIN)
   70 CONTINUE
C SET ARRAYS GR AND GDR
      CALL SETGRD(RVAC,ZCAN,RCAN,RMIN,ZLFT,ZRGT,RAD)
C
      LCFRON=.TRUE.
      M=MT
      IF(LCHIN.AND.(MROT.EQ.1)) M=4
      AM=FLOAT(M)
      UBTI=UBT
C
      CALL SETPRM(0,LELC,RAD,ZLFT,ZRGT,RMIN,RVAC)
C
C     IF(ZSEP.EQ.-9999.) ZSEP=ZL
      NZSEP=NINT(ZSEP/DDZ)
      ZSEP=FLOAT(NZSEP)*DDZ
C
      IF(REDUCE.EQ.0.0) WRITE(REDCHR,10100)
      IF(REDUCE.NE.0.0) WRITE(REDCHR,10110) REDUCE
C
      BSIG=FLOAT(ISIG)
      CDUM=BSHAPE(0)
      IF(UBTI.EQ.0.) UBTI=HLNGTH*2.
      IF(LCRBW) UBTI=AMAX1(UBTI,2.*HLNGTH+2.*FLOAT(NZ-NZSEP)*DDZ)
      IF(LSAV) LCPUTM=.TRUE.
      NW=NINT(UBTI/DDZ)+1
      NT=M*(NZ+NINT(UBTI/DDZ)+1)+1
      LELC=LPLE.OR.LPLC
      IF(TPS.NE.0.0) THEN
        NX1=NINT(TPS*FLOAT(M)*CLGHT/DDZ*1.E-12)
        TX2=FLOAT(NT-NX1)/FLOAT(NSHOT-1)
        TX1=TPS*FLOAT(M)*CLGHT/DDZ*1.E-12
      ELSE
        TX2=FLOAT(NT)/FLOAT(NSHOT)
        TX1=TX2
        NX1=NINT(TX2)
      END IF
C
      CALL SETPRM(2,LELC,RAD,ZLFT,ZRGT,RMIN,RVAC)
C
      CALL SETKMIN
C--TURN OFF NAPOLY INTEGRATION IF LNONAP=T
      IF(LNONAP) LNAPOLY=.FALSE.
C
      IF(CUTOFF.EQ.0.) CUTOFF=CUTOFFS
C CUTE IS UPPER LIMIT OF FREQUENCY SO THAT THE EXPONENTIAL CORRECT
C FACTOR DUE TO BUNCH LENGTH IS EXPFAC (DEFAULT=20).
      UPSM=0.5E-9*CLGHT/DDZ
      CUTE=SQRT(2.*LOG(EXPFAC))*0.30/(SIG*2.*PI)
      IF(UPSM.LT.CUTE) CUTE=UPSM
C     IF((CUTOFF.EQ.0.).OR.(CUTOFF.GT.CUTE)) CUTOFF=CUTE
      IF(CUTOFF.EQ.0.) CUTOFF=CUTE
      IF(LBEAMD) GO TO 80
C-----------------------------------------------------------------------
C PRINT OF PARAMETER
C-----------------------------------------------------------------------
C $FILE
      CALL PRFILE(LSAV,LREC,LCPUTM,TSOS,TMAX)
C-----------------------------------------------------------------------
C $BOUNDARY
      CALL PRBOUN
C-----------------------------------------------------------------------
C $MESH
   80 CALL PRMESH(NW,LBEAMD,RAD,MTSTAB,RDZDR)
      IF(LBEAMD) GO TO 90
C-----------------------------------------------------------------------
C #CAVITY / ##CAVITY
      IF(.NOT.LADD) WRITE (ND6,10010) CAVINP
      IF(LADD) WRITE (ND6,10020) CAVINP
      CALL PRNGEO(LADD,NLINE,ZL00,RZIN)
C-----------------------------------------------------------------------
C $BEAM
      CALL PRBEAM(ISIG)
C-----------------------------------------------------------------------
C $WAKE
      CALL PRWAKE(UBTI,NW)
C-----------------------------------------------------------------------
C PRINT TIME HERE IF $TIME USED WITHOUT NT ENABLING NT.GE.NTW
      CALL PRTIME(MTSTAB,RDZDR,LPLE,LPLC,TPS,NSHOT)
C-----------------------------------------------------------------------
C $PLOT
      CALL PRPLOT(LCAVIN,LCAVUS,LPLW,LFFT,LFFTA,LFFTT,LFFTL,LSPEC,
     +LPLE,LPLC,LPALL,NPLOT,CUTOFF,LFFTOR,REDUCE)
C-----------------------------------------------------------------------
C $PRINT
      CALL PRPRIN(LMATPR,LPRW,LSVW)
C-----------------------------------------------------------------------
C MESH GENERATING FOR...
C-----------------------------------------------------------------------
   90 CONTINUE
      IF(SYM) CALL EXPAND
      IF(SYM) CALL EXPAND2
  100 CONTINUE
C-----------------------------------------------------------------------
C CHECK VALIDITY OF GEOMETRY AND IF NAPOLY'S METHOD IS NEEDED
C-----------------------------------------------------------------------
      CALL EXAMIN (NLOOP,LBEAMD,LNTR,RMIN,RVAC,ZLFT,ZRGT,I1040)
      CALL SETKMIN
      AXY=GR(KMIN)
      IF(I1040.NE.0) GO TO 240
C-----------------------------------------------------------------------
C PLOT OF CAVITY SHAPE INPUT AND USED
C-----------------------------------------------------------------------
      LMATN=0
      IF (LMATPR) LMATN=1
      IF((NLOOP.EQ.0).AND.(LMATN.EQ.1)) CALL PLTMSH
      IF((NLOOP.GE.1).AND.LMESHD.AND.(LMATN.EQ.1)) CALL PLTMSH
      IF((NLOOP.EQ.0).AND.(LCAVIN))     CALL CPLOT(ZTXT)
      IF((NLOOP.GE.1).AND.LMESHD.AND.LCAVIN) CALL CPLOT(ZTXT)
      IF((NLOOP.EQ.0).AND.(LCAVUS))     CALL UPLOT(ZTXT)
      IF((NLOOP.GE.1).AND.LMESHD.AND.LCAVUS) CALL UPLOT(ZTXT)
C-----------------------------------------------------------------------
      CALL FLGCLR
      CALL INITJD
      CALL CALSUM
      CALL WAKECL
      IF(.NOT.LCHIN) CALL QKSET
      IF(LPLE.OR.LPLC) THEN
        QKMAX=0.
        DO 110 I=1,NW
          CALL CHAGE(1,(I-1)*M,QK)
          QKMAX=AMAX1(QK,QKMAX)
  110   CONTINUE
      END IF
C-----------------------------------------------------------------------
C RESTORE EM-FILEDS FROM A FILE IF LREC=.T. AND ADVANCE MESH
C-----------------------------------------------------------------------
      IF(ITEST.EQ.1) GO TO 170
      CALL ADVAN(LREC,IA,IB,NFST)
      IF(LNAPOLY.AND.(MROT.EQ.1).AND.LCRBW) THEN
        WRITE(ND6,10090)
        GO TO 240
      END IF
C-----------------------------------------------------------------------
C NOW START OF FIELD CALCULATION
C-----------------------------------------------------------------------
      DO 160 N=NFST,NT-1
        NSTK=N
C-MROT=0 CASE WITH NAPOLY'S INTEGRATION METHOD -------------------------
        IF(LNAPOLY.AND.(MROT.EQ.0)) THEN
          CALL MEFLD(N+1)
          CALL NAPL1(N,IA,NW)
          GO TO 120
        END IF
C-MROT=1 CASE WITH NAPOLY'S INTEGRATION METHOD -------------------------
        IF(LNAPOLY.AND.(MROT.EQ.1)) THEN
          CALL HEFLD(N+1)
          CALL NAPL2(N,IA,NW)
          GO TO 120
        END IF
C-MROT=0 CASE-----------------------------------------------------------
        IF(.NOT.LCHIN.OR.(MROT.EQ.0)) THEN
          IF (MROT.EQ.0) CALL MEFLD(N+1)
          IF (MROT.EQ.1) CALL HEFLD(N+1)
          CALL JMPBM(N,IA,NW)
        END IF
C-MROT=1 CASE-----------------------------------------------------------
        IF (LCHIN.AND.(MROT.EQ.1)) THEN
          CALL HEFLD(N+1)
          CALL CONBM(N,IA,IB,NW)
        END IF
C-----------------------------------------------------------------------
C CHECK CPU TIME LEFT
C-----------------------------------------------------------------------
  120   CONTINUE
        IF(.NOT.LCPUTM) GO TO 130
        CALL CPULFT(TSOS,TMAX,LSOS)
        IF(LSOS) GO TO 170
  130   CONTINUE
        IF((MROT.EQ.1).OR.(.NOT.LELC).OR.LCBACK) GO TO 160
        IF(LPLC) GO TO 150
        DO 140 I=1,NZL
          DO 140 J=1,NR
            K=(I-1)*NR+J
  140   HV(3,K)=HV(3,K)+HV(2,K)
  150   IF(N+1.LT.NX1) GO TO 160
        CALL EPLOT(ZTXT,NPLOT,N+1,LPLE,LPALL)
        TX1=TX1+TX2
        NX1=NINT(TX1)
  160 CONTINUE
  170 CONTINUE
C-----------------------------------------------------------------------
C SAVE RESULTS SO FAR INTO FT07
C-----------------------------------------------------------------------
      LSAVN=0
      IF (LSAV) LSAVN=1
      IF((NSTK.LT.NT-1).AND.(LSAVN.EQ.1)) THEN
        NSTK1=NSTK+1
        NBPN=(NBP-1)/400
        REWIND ND1
        WRITE(ND1) DDZ,DDR,RMARK,RZ,ZCAN,RCAN,SIG,ISIG,MT,MROT,RDRIVE,
     +  NBUNCH,BSEP
        DO 180 I=1,NBPN+1
          ILAST=I*400
          IF(ILAST.GT.NBP) ILAST=NBP
  180   WRITE(ND1) ((RZIN(J,L),J=1,72),L=(I-1)*400+1,ILAST)
        WRITE(ND1) ZTXT,R,Z,NZP0,NLINE,NSTK1,UBT,LCFRON,LCBACK,LCHIN,
     +  LNAPOLY
        WRITE(ND1) ZCF,ZCT,RWAK,LPLW,LCAVIN,LCAVUS,LFFTL,LFFTA,LFFTT,
     +  LFFT,LSPEC,LPLE,LPLC,NPLOT,TPS,NSHOT,NX1,TX1,TX2,LPALL,NVC0
        WRITE(ND1) CUTOFF,LPRW,LMATPR,SYM,LADD,NLOOP,LSVW,LSVF,LWNDW,
     +  NWFUN,ALPHA,LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,LINTZ,LCRBW,
     +  ZSEP,REDUCE,LPATH,LNONAP,EXPFAC
        WRITE(ND1) IA,IB
C RESTORE EM-FILEDS FROM A FILE IF LREC=.T. AND ADVANCE MESH
        CALL STORE
      END IF
C-----------------------------------------------------------------------
C CALCULATION OF THE TRANSVERSE WAKE POTENTIAL FROM THE LONGITUDINAL ONE
C-----------------------------------------------------------------------
      IF (MROT.EQ.0) GO TO 210
      IF (KOFF.NE.KOFL) THEN
        POTFAC=GR(JOFF)*GR(KMIN)*(1./GR(KOFF)**2-1./GR(KOFL)**2)
     +  /(2.*EPS*PI)
       DO 192 I=1,NW
          QK=0.
          IF(I.LE.IBT) CALL CHAGE(1,(I-1)*M,QK)
          WZ(I)=WZ(I)+QK*POTFAC
  192   CONTINUE
      END IF
      IF(.NOT.LCRBW) THEN
        WT(1)=-DDZ*WZ(1)/GR(KMIN)
        DO 190 I=2,NW
          WT(I)=WT(I-1)-DDZ*WZ(I)/GR(KMIN)
  190   CONTINUE
        NW1=NW-1
        DO 200 I=1,NW1
          WT(NW-I+1)=0.5*(WT(NW-I+1)+WT(NW-I))
  200   CONTINUE
        WT(1)=0.5*WT(1)
      ELSE
        DO 195 I=1,NW
          WT(I)=-WY(I)
  195   CONTINUE
      END IF
  210 CONTINUE
C-----------------------------------------------------------------------
C CALCULATION OF THE LONGITUDINAL WAKE POTENTIAL WITH POTENTIAL TERM DUE
C TO DIFFERENT BEAM PIPE RADII
C-----------------------------------------------------------------------
      LPOT=.FALSE.
      IF((MROT.EQ.0).AND.(LNAPOLY).AND.(KOFF.NE.KOFL)) LPOT=.TRUE.
      IF(LPOT) THEN
        DO 220 N=1,NW
          QK=0.
          IF(N.LE.IBT) CALL CHAGE(1,(N-1)*M,QK)
          WP(N)=WZ(N)-QK/(EPS*PI)*ALOG(GR(KOFL)/GR(KOFF))
          WZ(N)=WZ(N)-QK/(EPS*PI)*ALOG(GR(KOFL)/GR(KOFF))/2.
  220   CONTINUE
      END IF
C-----------------------------------------------------------------------
C PLOT WAKE POTENTIALS
C-----------------------------------------------------------------------
      CALL WKPLOT(NW,ZTXT,LPLW,LSVW,LPOT,LNTR,NSTK,AXY,ITEST)
C-----------------------------------------------------------------------
C PRINT WAKE POTENTIALS
C-----------------------------------------------------------------------
      IF(.NOT.LPRW) GO TO 230
      CALL PRINW(NW,LNTR,LPOT)
C-----------------------------------------------------------------------
C PLOT FOURIER TRANSFORM OF WAKE POTENTIALS
C-----------------------------------------------------------------------
  230 CONTINUE
      IF(ITEST.EQ.1) GO TO 240
      IF(LFFTOR.AND.(NBUNCH.GT.1)) GO TO 240
      IF(LFFT) LFFTL=.TRUE.
      IF(LFFT) LFFTA=.TRUE.
      IF(LFFT) LFFTT=.TRUE.
      IF(JOFF.GT.1) WFT=1./GR(JOFF)
      IF((MROT.EQ.1).AND.LFFTA.AND.(.NOT.LNTR)) CALL ZPLOT(ZTXT,WY,NW,1,
     +CUTOFF,WFT)
      IF((MROT.EQ.1).AND.LFFTT) CALL ZPLOT(ZTXT,WT,NW,2,CUTOFF,WFT)
      IF(LFFTL) CALL ZPLOT(ZTXT,WZ,NW,3,CUTOFF,WFT/AXY)
      IF(LSPEC.AND.(MROT.EQ.0).AND.(.NOT.LCRBW))
     + CALL ZPLOT(ZTXT,WZ,NW,5,CUTOFF,1.0)
      IF(LFFTL.AND.LPOT) CALL ZPLOT(ZTXT,WP,NW,4,CUTOFF,1.0)
      IF(LSPEC.AND.LPOT.AND.(.NOT.LCRBW))
     + CALL ZPLOT(ZTXT,WP,NW,6,CUTOFF,1.0)
  240 CONTINUE
C-----------------------------------------------------------------------
C REPEAT CALCULATION FOR ANOTHER SET OF BEAM PARAMETERS
C-----------------------------------------------------------------------
      LMESHD=.FALSE.
      LBEAMD=.FALSE.
      LWAKED=.FALSE.
      CUFOFF=0.0
      LNAPOLY=LNAPOLS
      CALL FLGINI(LFFTL,LFFTA,LFFTT,LFFTOR,LPLWOR,LSVWOR,LPLWLO,
     +LSVWLO,LFFT,LPLW,LSVW)
      CALL IRDINI
      DO 250 I=1,NF
  250 DDR(I)=DDRS(I)
      DDZ=DDZS
      NR=NRS
      NZ=NZS
      DO 260 I=1,NF
  260 RMARK(I)=RMAK0(I)
C
      READ(ND5,10060) ENDFIL
      IF(ENDFIL.NE.CONTIN) GO TO 310
      IF(ENDFIL.EQ.STOPC) GO TO 310
      IF(ENDFIL.EQ.CONTIN) READ(ND5,10060) ENDFIL
      IF((ENDFIL.NE.DMESH).AND.(ENDFIL.NE.AMESH)) GO TO 270
      BACKSPACE(ND5)
      READ(ND5,MESH)
      LMESHD=.TRUE.
      READ(ND5,10060) ENDFIL
      IF((ENDFIL.EQ.STOPC).OR.(ENDFIL.EQ.CONTIN)) THEN
        BACKSPACE(ND5)
        NLOOP=NLOOP+1
        GO TO 60
      END IF
  270 IF((ENDFIL.NE.BEAMD).AND.(ENDFIL.NE.BEAMA)) GO TO 280
      BACKSPACE(ND5)
      READ(ND5,BEAM)
      IF(.NOT.LMESHD) LBEAMD=.TRUE.
      READ(ND5,10060) ENDFIL
      IF((ENDFIL.EQ.STOPC).OR.(ENDFIL.EQ.CONTIN)) THEN
        BACKSPACE(ND5)
        NLOOP=NLOOP+1
        GO TO 60
      END IF
  280 IF((ENDFIL.NE.WAKED).AND.(ENDFIL.NE.WAKEA)) GO TO 290
      BACKSPACE(ND5)
      READ(ND5,WAKE)
      LBEAMD=.FALSE.
      READ(ND5,10060) ENDFIL
      IF((ENDFIL.EQ.STOPC).OR.(ENDFIL.EQ.CONTIN)) THEN
        BACKSPACE(ND5)
        NLOOP=NLOOP+1
        GO TO 60
      END IF
  290 IF((ENDFIL.NE.PLOTD).AND.(ENDFIL.NE.PLOTA)) GO TO 300
      BACKSPACE(ND5)
      READ(ND5,PLOT)
      READ(ND5,10060) ENDFIL
      IF((ENDFIL.EQ.STOPC).OR.(ENDFIL.EQ.CONTIN)) THEN
        BACKSPACE(ND5)
        NLOOP=NLOOP+1
        GO TO 60
      END IF
  300 IF((ENDFIL.NE.PRIND).AND.(ENDFIL.NE.PRINA)) GO TO 310
      BACKSPACE(ND5)
      READ(ND5,PRIN)
      NLOOP=NLOOP+1
      GO TO 60
  310 CONTINUE
      CALL CPUTIM(TT,3)
      WRITE(ND6,10070) TT
      STOP
10000 FORMAT(A8)
10010 FORMAT(/1X,'#CAVITYSHAPE (##CAVITYSHAPE) :  ',A16)
10020 FORMAT(/1X,'@CAVITYSHAPE (@@CAVITYSHAPE) :  ',A16)
10030 FORMAT(9A8)
10040 FORMAT(1H ,76(1H*)/1X,1H*,74X,1H* /1X,2H* ,
     +'                                 ABCI             ' ,22(1H ),2H *
     +/1X,2H* ,72(1H ),2H * /1X,2H* ,
     +'Azimuthal Beam Cavity Interaction in a cylindricall' ,
     +'y symmetric structure',2H * /1X,2H* ,72(1H ),2H * /1X,2H* ,9(A8),
     +2H * /1X,2H* ,72(1H ),2H * /1X,2H* ,'    DATE:',A8,'   TIME:',A8,
     +3X,'  VERSION 9.4,  JANUARY  1995  ',4X,1(1H ),2H * /1X,2H* ,72(
     +1H ),2H * /1X,76(1H*)/)
10050 FORMAT(1H1,76(1H*) /1X,2H* ,72(1H ),2H * /1X,2H* ,
     +'                           ABCI  ',I2,A2,' RUN    ' ,27(1H ),2H *
     +/1X,2H* ,72(1H ),2H * /1X,2H* ,'  DATE:',A8,'   TIME:',A8, 3X,
     +'  VERSION 9.4,  JANUARY  1995  ',6X,1(1H ),2H * /1X,2H* ,72(1H ),
C    2 /1X,2H* ,30(1H ),I2,'-TH RUN',33(1H ),2H *
     +2H * /1X,76(1H*)/)
10060 FORMAT(A4)
10070 FORMAT(1H1/1X,130(1H*) //1X,'ABCI ENDED WITHOUT ERROR',10X,
     +'TOTAL CPU TIME USED:',F10.1,' (S)' //1X,130(1H*))
10080 FORMAT(1X,'*** STOP ***  THE NUMBER OF MESH INPUT LINES ', 'NZP:',
     +I5,' EXCEEDS MAXIMUM NBP:',I5)
10090 FORMAT(/1X,'*** SKIP ***  LCRBW CANNOT BE TRUE WHEN LNAPOLY IS ',
     +'TRUE AND MROT=1 IN THIS VERSION OF ABCI.'/)
10100 FORMAT('             ')
10110 FORMAT('REDUCE ',F6.3)
      END
C
      BLOCKDATA DEVNUM
C-----------------------------------------------------------------------
C BLOCK DATA FOR LOGICAL NUMBER OF EACH DEVICE
C
C ND1 for saving intermediate results
C ND2 for recovering
C ND5 for input data
C ND6 for output on line printer
C ND8 for saving wake potentials and FFT results
C ND9 for TopDrawer input
C-----------------------------------------------------------------------
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
C-----------------------------------------------------------------------
      DATA ND1,ND2,ND5,ND6,ND8,ND9/1,2,5,6,8,9/
      END
C
      SUBROUTINE PRMINI(TMAX,TSOS,ITEST,UBT,CUTOFF,SYM,LADD,CAVINP,
     +LSAV,LREC,LCPUTM,LCFRON,LCBACK,LPRW,LPLW,LSVF,LSVW,LFFT,LFFTL,
     +LFFTA,LFFTT,LSPEC,LPLE,LPLC,LELC,LPALL,LWNDW,LCAVIN,LCAVUS,
     +LMATPR,LCHIN,LNAPOLY,LMESHD,LBEAMD,LWAKED,RDRIVE,ZCF,ZCT,RWAK,
     +TPS,NX1,NSHOT,ISIG,BSIG,NBUNCH,NPLOT,BSEP,NFST,RAD,ZL,
     +MT,MROT,NLOOP,NWFUN,ALPHA,LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,
     +LINTZ,LCRBW,ZSEP,REDUCE,LPATH,LNONAP,EXPFAC)
C-----------------------------------------------------------------------
C INITIALIZE ALL PARAMETERS
C-----------------------------------------------------------------------
      LOGICAL LSAV,LREC,LCFRON,SYM,LCBACK,LPLW,LPRW,LCHIN,LMATPR
      LOGICAL LMESHD,LNAPOLY,LBEAMD,LWAKED,LCAVIN,LCAVUS
      LOGICAL LADD,LPLC,LELC,LPALL,LSVW,LSVF,LWNDW,LPATH,LNONAP
      LOGICAL LFFTL,LFFTA,LFFTT,LFFT,LCPUTM,LSPEC,LPLE
      LOGICAL LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,LINTZ,LCRBW
      CHARACTER*16 CAVFUL,CAVINP
      DATA CAVFUL/'FULL CELL INPUT '/
      TMAX=3600.
      TSOS=5.
      ITEST=0
      UBT=0.
      CUTOFF=0.0
      SYM=.FALSE.
      LADD=.FALSE.
      CAVINP=CAVFUL
      LSAV=.FALSE.
      LREC=.FALSE.
      LCPUTM=.FALSE.
      LCFRON=.TRUE.
      LCBACK=.TRUE.
      LPRW=.FALSE.
      LPLW=.FALSE.
      LPLWL=.FALSE.
      LPLWA=.FALSE.
      LPLWT=.FALSE.
      LSVF=.FALSE.
      LSVW=.FALSE.
      LSVWL=.FALSE.
      LSVWA=.FALSE.
      LSVWT=.FALSE.
      LFFT=.FALSE.
      LFFTL=.FALSE.
      LFFTA=.FALSE.
      LFFTT=.FALSE.
      LSPEC=.FALSE.
      LINTZ=.FALSE.
      REDUCE=0.0
      LCRBW=.FALSE.
      LPATH=.FALSE.
      LNONAP=.FALSE.
C     ZSEP=-9999.
      ZSEP=0.0
      LPLE=.FALSE.
      LPLC=.FALSE.
      LELC=.FALSE.
      LPALL=.FALSE.
      LWNDW=.TRUE.
      LCAVIN=.FALSE.
      LCAVUS=.FALSE.
      LMATPR=.FALSE.
      LCHIN=.FALSE.
      LNAPOLY=.FALSE.
      LMESHD=.FALSE.
      LBEAMD=.FALSE.
      LWAKED=.FALSE.
      RDRIVE=-1.
      ZCF=-1.
      ZCT=-1.
      RWAK=-1.
      TPS=0.0
      NX1=0
      NSHOT=5
      ISIG=5
      BSIG=FLOAT(ISIG)
      NBUNCH=1
      NPLOT=10
      BSEP=0.0
      NFST=0
      RAD=0.
      ZL=0.
      MT=3
      MROT=1
      NLOOP=0
      NWFUN=1
      ALPHA=3.E0
      EXPFAC=20.
      RETURN
      END
C
      BLOCKDATA PHYPRM
C-----------------------------------------------------------------------
C BLOCK DATA FOR PHYSICAL PARAMETERS
C-----------------------------------------------------------------------
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
C-----------------------------------------------------------------------
      DATA EPS/8.854187818E-12/,CLGHT/2.99792458E8/,CEPS/1.E-6/,
     +PI/3.14159265359E0/,ZO/376.7303134E0/
      END
C
      SUBROUTINE FLGINI(LFFTL,LFFTA,LFFTT,LFFTOR,LPLWOR,LSVWOR,LPLWLO,
     +LSVWLO,LFFT,LPLW,LSVW)
C-----------------------------------------------------------------------
C INITIALIZE ALL FLAGS
C-----------------------------------------------------------------------
      LOGICAL LFFTL,LFFTA,LFFTT,LFFTOR,LPLWOR,LSVWOR,LPLWLO,LSVWLO
      LOGICAL LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
      LOGICAL LFFT,LPLW,LSVW
      COMMON /WPS/ LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
      IF(LFFT) THEN
        LFFTL=.FALSE.
        LFFTA=.FALSE.
        LFFTT=.FALSE.
        LFFTOR=.FALSE.
      END IF
      IF(LPLW) THEN
        LPLWL=.FALSE.
        LPLWA=.FALSE.
        LPLWT=.FALSE.
        LPLWOR=.FALSE.
        LPLWLO=.FALSE.
      END IF
      IF(LSVW) THEN
        LSVWL=.FALSE.
        LSVWA=.FALSE.
        LSVWT=.FALSE.
        LSVWOR=.FALSE.
        LSVWLO=.FALSE.
      END IF
      RETURN
      END
C
      SUBROUTINE FLGCLR
C-----------------------------------------------------------------------
C CLEAR FLAG
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      DO 10 K=1,NMESH
        IFLG(K)=5
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE INITHT
C-----------------------------------------------------------------------
C INITIRIZATION OF IHIT
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      LOGICAL IHIT
      COMMON /BC/ IHIT(NHIT,8)
      DO 10 K=1,NHIT
        DO 10 L=1,8
          IHIT(K,L)=.FALSE.
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE RZINIT
C-----------------------------------------------------------------------
C INITIRIZATION OF ARRAY RZ
C-----------------------------------------------------------------------
      PARAMETER (NVR=100)
      COMMON /VAR/ RZ(NVR)
      DO 10 K=1,NVR
        RZ(K)=-8888.
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE WAKECL
C-----------------------------------------------------------------------
C INITIRIZATION OF WZ TAND WY
C-----------------------------------------------------------------------
      PARAMETER (NWN=500000)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      DO 10 I=1,NWN
        WZ(I)=0.
        WY(I)=0.
        WT(I)=0.
        WP(I)=0.
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE INITJD
C-----------------------------------------------------------------------
C INITIRIZATION OF THE FLAGS FOR JUDGING OF BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NPT=20)
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      DO 20 L=1,3
        DO 10 K=1,NMESH
          SJDG(L,K)=0.
   10   CONTINUE
   20 CONTINUE
      DO 40 L=1,3
        DO 30 K=1,NMESH
          HV(L,K)=0.
          EV(L,K)=0.
   30   CONTINUE
   40 CONTINUE
      DO 70 I=1,NPT
        DO 60 L=1,2
          DO 50 K=1,NRN
            HVSTK(K,L,I)=0.
            HVSTL(K,L,I)=0.
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      RETURN
      END
C
      SUBROUTINE INPUTI(LADD,NZP,NLINE,NVC,ZCAN,RCAN,RZIN)
C-----------------------------------------------------------------------
C READ STRUCTURE GEOMETRY FROM INPUT FILE
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      CHARACTER*1 ABCD(72),RZIN(72,NBP),CMT
      LOGICAL LADD,LRZF,LREPS,LREPE,SYM
      DIMENSION XY(2),IXY(4)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /SYMT/ SYM
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      DATA CMT/'%'/
      N=0
      K=0
      NREP=0
      NPOL=0
      LRZF=.TRUE.
      LREPS=.FALSE.
      LREPE=.FALSE.
      ZCAN=0.
      RCAN=0.
      READ(ND5,*) RDUM
   10 READ(ND5,10000) ABCD
      IF(ABCD(1).EQ.CMT) GO TO 10
      IF(N+1.GT.NBP) GO TO 110
      DO 20 J=1,72
   20 RZIN(J,N+1)=ABCD(J)
   30 N=N+1
      CALL FLOT(ABCD,72,XY,IXY,IERROR)
      IF((IXY(4).GE.0).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).GE.0).AND.LADD) THEN
        IF(IERROR.LT.0) GO TO 50
        IF(LREPE) THEN
          WRITE(ND6,10010)
          GO TO 50
        ELSE
        END IF
        NREP=IXY(4)
        NSS=N
        LREPS=.TRUE.
        GO TO 10
      ELSE
      END IF
      IF((IXY(4).EQ.-9999).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).EQ.-9999).AND.LADD) THEN
        IF(IERROR.LT.0) GO TO 50
        IF(.NOT.LREPS) THEN
          WRITE(ND6,10020)
          GO TO 50
        ELSE
        END IF
        NREP=NREP-1
        LREPE=.TRUE.
        IF(NREP.LE.0) THEN
          LREPE=.FALSE.
          LREPS=.FALSE.
          GO TO 10
        ELSE
        END IF
        N=NSS
        DO 40 J=1,72
   40   ABCD(J)=RZIN(J,NSS+1)
        GO TO 30
      ELSE
      END IF
      IF(IERROR.GE.0) GO TO 80
   50 WRITE(ND6,10030)
      DO 70 L=1,N
        DO 60 J=1,72
   60   ABCD(J)=RZIN(J,L)
        WRITE(ND6,10040) L, ABCD
   70 CONTINUE
      STOP
   80 CONTINUE
      K=K+1
      IF(K.GT.NBP) GO TO 110
      IF((.NOT.LADD).OR.(IXY(3).LE.-1).OR.(XY(1).EQ.9999.).OR.LRZF)
     +THEN
        R(K)=XY(1)
        Z(K)=XY(2)
      ELSE
        IF(R(K-1).GT.-1.) THEN
          R(K)=R(K-1)*1.E2+XY(1)*1.E2
          Z(K)=Z(K-1)*1.E2+XY(2)*1.E2
          R(K)=R(K)*1.E-2
          Z(K)=Z(K)*1.E-2
        ELSE
          R(K)=R(K-2)*1.E2+XY(1)*1.E2
          Z(K)=Z(K-2)*1.E2+XY(2)*1.E2
          R(K)=R(K)*1.E-2
          Z(K)=Z(K)*1.E-2
        END IF
        IF(ABS(R(K)).LE.CEPS) R(K)=0.0
        IF(ABS(Z(K)).LE.CEPS) Z(K)=0.0
      END IF
      IF((NREP.GT.0).AND.LREPE) THEN
        DO 90 J=1,72
   90   ABCD(J)=RZIN(J,N+1)
        GO TO 30
      ELSE
        IF(R(K).EQ.9999.) GO TO 100
        IF(R(K).LT.0.) GO TO 10
        IF(K.GT.1) THEN
          IF(R(K-1).LE.-3.) GO TO 10
        ELSE
        END IF
        IF(RCAN.LT.R(K)) RCAN=R(K)
        IF(ZCAN.LT.Z(K)) ZCAN=Z(K)
        IF(LRZF) THEN
          RDUM=R(K)
          ZDUM=Z(K)
          LRZF=.FALSE.
        ELSE
          IF((ABS(R(K)-RDUM).LE.CEPS).AND.(ABS(Z(K)-ZDUM).LE.CEPS))
     +    THEN
            NPOL=NPOL+1
            IF(NPOL.EQ.1) NVC=K-2
          END IF
        END IF
        GO TO 10
      END IF
  100 NZP=K-1
      NLINE=N-1
      IF(NPOL.EQ.0) NVC=NZP-2
      IF(SYM.AND.(2*NZP-3.GT.NBP)) GO TO 110
      RETURN
  110 WRITE (ND6,10050) NBP
      STOP
10000 FORMAT(72A1)
10010 FORMAT(1X,'*** STOP *** YOU HAVE "<" BEFORE ">"')
10020 FORMAT(1X,'*** STOP *** YOU MAY HAVE FORGOTTEN ">" BEFORE "<"')
10030 FORMAT(1X,'??? DUMP OF THE MESH INPUT SO FAR ???'/)
10040 FORMAT(1X,I5,':',5X,72A1)
10050 FORMAT(1X,'*** STOP *** YOUR GEOMETRY INPUT IS TOO LONG.',
     +' INCLEASE THE PARAMETER NBP (CURRENTLY ONLY ',I7,' )')
      END
C
      SUBROUTINE INPUTD(LADD,NZP,NLINE,NVC,ZCAN,RCAN,RZIN)
C-----------------------------------------------------------------------
C READ STRUCTURE GEOMETRY FROM ARRAY RZIN
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      CHARACTER*1 ABCD(72),RZIN(72,NBP)
      LOGICAL LADD,LRZF
      DIMENSION XY(2),IXY(4)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      N=0
      K=0
      NPOL=0
      NREP=0
      LRZF=.TRUE.
      RCAN=0.0
      ZCAN=0.0
   10 N=N+1
      IF(N.GT.NLINE) GO TO 30
      DO 20 J=1,72
   20 ABCD(J)=RZIN(J,N)
      CALL FLOT(ABCD,72,XY,IXY,IERROR)
      IF((IXY(4).GE.0).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).GE.0).AND.LADD) THEN
        NREP=IXY(4)
        NSS=N
        GO TO 10
      ELSE
      END IF
      IF((IXY(4).EQ.-9999).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).EQ.-9999).AND.LADD) THEN
        NREP=NREP-1
        IF(NREP.LE.0) THEN
          GO TO 10
        ELSE
        END IF
        N=NSS
        GO TO 10
      ELSE
      END IF
      K=K+1
      IF((.NOT.LADD).OR.(IXY(3).LE.-1).OR.(XY(1).EQ.9999.).OR.LRZF)
     +THEN
        R(K)=XY(1)
        Z(K)=XY(2)
      ELSE
        IF(R(K-1).GT.-1.) THEN
          R(K)=R(K-1)*1.E2+XY(1)*1.E2
          Z(K)=Z(K-1)*1.E2+XY(2)*1.E2
          R(K)=R(K)*1.E-2
          Z(K)=Z(K)*1.E-2
        ELSE
          R(K)=R(K-2)*1.E2+XY(1)*1.E2
          Z(K)=Z(K-2)*1.E2+XY(2)*1.E2
          R(K)=R(K)*1.E-2
          Z(K)=Z(K)*1.E-2
        END IF
        IF(ABS(R(K)).LE.CEPS) R(K)=0.0
        IF(ABS(Z(K)).LE.CEPS) Z(K)=0.0
      END IF
      IF(R(K).LT.0.) GO TO 10
      IF(K.GT.1) THEN
        IF(R(K-1).LE.-3.) GO TO 10
      ELSE
      END IF
      IF(RCAN.LT.R(K)) RCAN=R(K)
      IF(ZCAN.LT.Z(K)) ZCAN=Z(K)
      IF(LRZF) THEN
        RDUM=R(K)
        ZDUM=Z(K)
        LRZF=.FALSE.
      ELSE
        IF((ABS(R(K)-RDUM).LE.CEPS).AND.(ABS(Z(K)-ZDUM).LE.CEPS))
     +  THEN
          NPOL=NPOL+1
          IF(NPOL.EQ.1) NVC=K-2
        END IF
      END IF
      GO TO 10
   30 CONTINUE
      IF(NPOL.EQ.0) NVC=NZP-2
      RETURN
      END
C
      SUBROUTINE PRNGEO(LADD,NLINE,ZL00,RZIN)
C-----------------------------------------------------------------------
C PRINT STRUCTURE GEOMETRY ON LINEPRINTER
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      CHARACTER*1 ABCD(72),RZIN(72,NBP),BLK,MNS
      CHARACTER*3 FG1,FG2
      LOGICAL LADD,LCAV,LVEX
      DIMENSION XY(2),IXY(4)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      DATA MNS/'-'/,BLK/' '/
      WRITE (ND6,10010)
      I=1
      NREP=0
      K=1
   10 IF(I.GT.NLINE) GO TO 30
      DO 20 J=1,72
   20 ABCD(J)=RZIN(J,I)
      CALL FLOT(ABCD,72,XY,IXY,IERROR)
      IF((IXY(4).GE.0).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).GE.0).AND.LADD) THEN
        NREP=IXY(4)
        NSS=I
        I=I+1
        GO TO 10
      ELSE
      END IF
      IF((IXY(4).EQ.-9999).AND.(.NOT.LADD)) GO TO 10
      IF((IXY(4).EQ.-9999).AND.LADD) THEN
        NREP=NREP-1
        IF(NREP.LE.0) I=I+1
        IF(NREP.LE.0) GO TO 10
        I=NSS+1
        GO TO 10
      ELSE
      END IF
      IF(IXY(1).NE.0) THEN
        FG1=BLK
        IF(IXY(1).LT.0) FG1=MNS
        IXY(1)=IABS(IXY(1))
      ELSE
      END IF
      IF(IXY(2).NE.0) THEN
        FG2=BLK
        IF(IXY(2).LT.0) FG2=MNS
        IXY(2)=IABS(IXY(2))
      ELSE
      END IF
C     K=K+1
      IF(R(K).GT.-1.) THEN
        FRR=GR(IRR(K))
        FZZ=FLOAT(IZZ(K)-1)*DDZ
        IF(K.EQ.NVC) FZZ=FLOAT(IZZ(K))*DDZ
        ZZZZ=Z(K)
        IF(K.EQ.NVC) ZZZZ=ZL00
        IF(K.EQ.NVC+1) ZZZZ=ZL00
        ABSR=ABS(R(K)-FRR)
        ABSZ=ABS(ZZZZ-FZZ)
        IF(.NOT.LADD) THEN
          IF((ABSR.LE.CEPS).AND.(ABSZ.LE.CEPS)) THEN
            IF((IXY(1).EQ.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10020) R (I),
     +      ZZZZ,IRR(I),IZZ(I)
            IF((IXY(1).NE.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10050) FG1,
     +      IXY(1),ZZZZ,IRR(I),IZZ(I),R(I),ZZZZ
            IF((IXY(1).EQ.0).AND.(IXY(2).NE.0)) WRITE (ND6,10060) R (I),
     +      FG2,IXY(2),IRR(I),IZZ(I),R(I),ZZZZ
            IF((IXY(1).NE.0).AND.(IXY(2).NE.0)) WRITE (ND6,10070) FG1,
     +      IXY(1),FG2,IXY(2),IRR(I),IZZ(I),R(I),ZZZZ
          ELSE
            IF((IXY(1).EQ.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10030) R (I),
     +      ZZZZ,IRR(I),IZZ(I),FRR,FZZ
            IF((IXY(1).NE.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10080) FG1,
     +      IXY(1),ZZZZ, IRR(I),IZZ(I),R(I),ZZZZ,FRR,FZZ
            IF((IXY(1).EQ.0).AND.(IXY(2).NE.0)) WRITE (ND6,10090) R (I),
     +      FG2,IXY(2), IRR(I),IZZ(I),R(I),ZZZZ,FRR,FZZ
            IF((IXY(1).NE.0).AND.(IXY(2).NE.0)) WRITE (ND6,10100) FG1,
     +      IXY(1),FG2,IXY(2), IRR(I),IZZ(I),R(I),ZZZZ,FRR,FZZ
          END IF
        ELSE
          IF(K.EQ.NVC+2) XY(2)=-FLOAT(IZZ(K-1)-1)*DDZ
          IF((ABSR.LE.CEPS).AND.(ABSZ.LE.CEPS)) THEN
            IF((IXY(1).EQ.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10030) XY
     +      (1),XY(2),IRR(K),IZZ(K),R(K),ZZZZ
            IF((IXY(1).NE.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10050) FG1,
     +      IXY(1),XY(2),IRR(K),IZZ(K),R(K),ZZZZ
            IF((IXY(1).EQ.0).AND.(IXY(2).NE.0)) WRITE (ND6,10060) XY
     +      (1),FG2,IXY(2),IRR(K),IZZ(K),R(K),ZZZZ
            IF((IXY(1).NE.0).AND.(IXY(2).NE.0)) WRITE (ND6,10070) FG1,
     +      IXY(1),FG2,IXY(2),IRR(K),IZZ(K),R(K),ZZZZ
          ELSE
            IF((IXY(1).EQ.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10040) XY
     +      (1),XY(2),IRR(K),IZZ(K),R(K),Z(K),FRR,FZZ
 
            IF((IXY(1).NE.0).AND.(IXY(2).EQ.0)) WRITE (ND6,10080) FG1,
     +      IXY(1),XY(2), IRR(K),IZZ(K),R(K),ZZZZ,FRR,FZZ
 
            IF((IXY(1).EQ.0).AND.(IXY(2).NE.0)) WRITE (ND6,10090) XY
     +      (1),FG2,IXY(2), IRR(K),IZZ(K),R(K),ZZZZ,FRR,FZZ
 
            IF((IXY(1).NE.0).AND.(IXY(2).NE.0)) WRITE (ND6,10100) FG1,
     +      IXY(1),FG2,IXY(2), IRR(K),IZZ(K),R(K),ZZZZ,FRR,FZZ
 
          END IF
        END IF
        I=I+1
        K=K+1
      ELSE
        IKIND=R(K)
        IF(IKIND.GE.-2) THEN
          LCAV=Z(K).GT.0.0
          LVEX=Z(K).LT.0.0
          IF(LCAV.AND.(IXY(2).EQ.0)) WRITE (ND6,10110) IKIND,Z(K)
          IF(LVEX.AND.(IXY(2).EQ.0)) WRITE (ND6,10120) IKIND,Z(K)
          IF(LCAV.AND.(IXY(2).NE.0)) WRITE (ND6,10130) IKIND,IXY(2),Z
     +    (K)
          IF(LVEX.AND.(IXY(2).NE.0)) WRITE (ND6,10140) IKIND,IXY(2),Z
     +    (K)
          I=I+1
          K=K+1
        ELSE
          IF(.NOT.LADD) WRITE (ND6,10150) IKIND,R(K+1),Z(K+1)
          IF(LADD) WRITE (ND6,10150) IKIND,R(K+1),Z(K+1)
          I=I+2
          K=K+2
        END IF
      END IF
      GO TO 10
   30 CONTINUE
      IF(LADD.AND.(ABS(R(NVC+2)).GE.0.5*GDR(1))) THEN
        N=NVC+2
        WRITE(ND6,10000)
        STOP
      ELSE
      END IF
      RETURN
10000 FORMAT(1X,'*** STOP *** THE FIRST POLYGON IS NOT CLOSED.'/1X,
     +'THERE IS A VERTICAL GAP MORE THAN DDR/2 AT RLC.')
10010 FORMAT(/1X,'INPUT:',7X,'(R, Z)',19X,'(IR, IZ)')
10020 FORMAT(1X,'( ',G13.5,G13.5,' )  ( ',I7,I7,' )')
10030 FORMAT(1X,'( ',G13.5,G13.5,' )  ( ',I7,I7,' ) ---> ', '( ',G13.5,
     +G13.5,' )')
10040 FORMAT(1X,'( ',G13.5,G13.5,' )  ( ',I7,I7,' ) ---> ', '( ',G13.5,
     +G13.5,' ) ---> ', '( ',G13.5,G13.5,' )')
10050 FORMAT(1X,'( ',4X,A1,'RZ',I3,3X,G13.5, ' )  ( ',I7,I7,' ) ---> ',
     +'( ',G13.5,G13.5,' )')
10060 FORMAT(1X,'( ',G13.5,4X,A1,'RZ',I3,3X, ' )  ( ',I7,I7,' ) ---> ',
     +'( ',G13.5,G13.5,' )')
10070 FORMAT(1X,'( ',4X,A1,'RZ',I3,3X,4X,A1,'RZ',I3,3X,' ) ', ' ( ',I7,
     +I7,' ) ---> ', '( ',G13.5,G13.5,' )')
10080 FORMAT(1X,'( ',4X,A1,'RZ',I3,3X,G13.5, ' )  ( ',I7,I7,' ) ---> ',
     +'( ',G13.5,G13.5,' ) ---> ', '( ',G13.5,G13.5,' )')
10090 FORMAT(1X,'( ',G13.5,4X,A1,'RZ',I3,3X, ' )  ( ',I7,I7,' ) ---> ',
     +'( ',G13.5,G13.5,' ) ---> ', '( ',G13.5,G13.5,' )')
10100 FORMAT(1X,'( ',4X,A1,'RZ',I3,3X,4X,A1,'RZ',I3,3X,' ) ', ' ( ',I7,
     +I7,' ) ---> ', '( ',G13.5,G13.5,' ) ---> ', '( ',G13.5,G13.5,' )')
10110 FORMAT(1X,'CONNECTED BY A CONCAVE CIRCLE, IKIND=',I2, ',  RADIUS='
     +,G13.5)
10120 FORMAT(1X,'CONNECTED BY A CONVEX  CIRCLE, IKIND=',I2, ',  RADIUS='
     +,G13.5)
10130 FORMAT(1X,'CONNECTED BY A CONCAVE CIRCLE, IKIND=',I2,
     +',  RADIUS= RZ',I3,' --->  ',G13.5)
10140 FORMAT(1X,'CONNECTED BY A CONVEX  CIRCLE, IKIND=',I2,
     +',  RADIUS= RZ',I3,' --->  ',G13.5)
10150 FORMAT(1X,'CONNECTED BY AN ELLIPSE,       IKIND=',I2,
     +',  (RM, ZM)=     ( ',G13.5,G13.5,' )')
      END
C
      SUBROUTINE SETGRD(RVAC,ZCAN,RCAN,RMIN,ZLFT,ZRGT,RAD)
C-----------------------------------------------------------------------
C SET GR AND GRD
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /SYMT/ SYM
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      IF(SYM) THEN
        RVAC=R(NVC)
      ELSE
        RVAC=AMAX1(R(2),R(NVC))
        IF(ABS(R(2)-R(NVC)).LE.CEPS) RVAC=AMIN1(R(2),R(NVC))
      END IF
      RAD2=0.0
      DO 10 I=1,NF
   10 IF(RAD2.LT.RMARK(I)) RAD2=RMARK(I)
      IF(RAD2.LT.RCAN) RAD2=RCAN
      NNX=0
   20 NNX=NNX+1
      IF(NNX.GE.3) GO TO 30
      CALL IRDSET(RAD2)
      CALL CENTER(RAD2,ZCAN,RMAX,ZMAX,RMIN,RVAC,ZLFT,ZRGT,NNX)
      IF((ABS(RAD2-RMAX).LE.CEPS).AND.(ABS(ZCAN-ZMAX).LE.CEPS)) GO TO 30
      RAD2=RMAX
      ZCAN=ZMAX
      GO TO 20
   30 RAD=RMAX
      ZL=ZMAX
      RETURN
      END
C
      SUBROUTINE SETPRM(KTYPE,LELC,RAD,ZLFT,ZRGT,RMIN,RVAC)
C-----------------------------------------------------------------------
C SET SOME IMPORTANT PARAMETERS
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM,LELC
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /SYMT/ SYM
      COMMON /WPL/ MROT
      COMMON /OFFS/ RDRIVE,ZCF,ZCT,RWAK
      COMMON /PIPE/ KOFF,KOFL
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /ICON/ IZLFT,IZRGT,IRMIN
      IF(KTYPE-1) 10,30,30
C DDR: STEP SIZE IN R
C DDZ: STEP SIZE IN Z
   10 IF(DDZ.EQ.0.) THEN
C FOR COMPATIBILITY WITH TBCI
        NZ=NZ-1
        DDZ=ZL/FLOAT(NZ)
        IF(SYM) DDZ=DDZ*2.
        DDR(1)=RAD/FLOAT(NR-1)
      ELSE
        NZ=NINT(ZL/DDZ)
        IF(SYM) NZ=NZ*2
        NR=IRMAX
      END IF
      ZL=FLOAT(NZ)*DDZ
      DO 20 N=1,NZP
        IF(R(N).LT.0.) GO TO 20
        IZZ(N)=NINT(Z(N)/DDZ)+1
        IRR(N)=IRD(R(N))
   20 CONTINUE
      IF(SYM) THEN
        IZZ(NVC+1)=NZ/2+1
        IZZ(NVC)=NZ/2
        Z(NVC+1)=ZL/2.
        Z(NVC)=ZL/2.
      ELSE
        IZZ(NVC+1)=NZ+1
        IZZ(NVC)=NZ
        Z(NVC+1)=ZL
        Z(NVC)=ZL
      END IF
      Z(NVC+2)=0.
      IZZ(NVC+2)=1
      IF(SYM) THEN
        ZLFT=ZL/2.-ZRGT
        ZRGT=ZL/2.+ZRGT
      ELSE
      END IF
      IZLFT=NINT(ZLFT/DDZ)+1
      IZRGT=NINT(ZRGT/DDZ)+1
      IRMIN=IRD(RMIN)
      RETURN
   30 IRVAC=IRD(RVAC)
C JOFF:RADIUS IN WHICH A BEAM GOES THROUGH A CAVITY
C KOFF:RADIUS IN WHICH A WAKE POTENTIAL IS INTEGRATED
      IF(SYM) THEN
        JOFF=IRR(NVC)
      ELSE
        JOFF=IRR(2)
      END IF
      KOFF=JOFF
      KOFL=KOFF
      IF(.NOT.SYM) KOFL=IRR(NVC)
      IF(RDRIVE.LT.0.AND.MROT.EQ.0.AND.LELC) JOFF=1
      IF(RDRIVE.GE.0.) JOFF=IRD(RDRIVE)
      IF((JOFF.EQ.1).AND.(MROT.EQ.1)) JOFF=KOFF
      IF((JOFF.EQ.2).AND.(MROT.EQ.1)) JOFF=KOFF
      IF(IRMIN.LT.JOFF) JOFF=IRMIN
C
      IF(RWAK.GE.0.) IRWAK=IRD(RWAK)
      IF(RWAK.LT.0.) IRWAK=MIN(IRMIN,IRVAC)
      RETURN
      END
C
      SUBROUTINE PRFILE(LSAV,LREC,LCPUTM,TSOS,TMAX)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &FILE
C-----------------------------------------------------------------------
      LOGICAL LSAV,LREC,LCPUTM
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      WRITE (ND6,10000)
      WRITE (ND6,10010) LSAV,LREC,LCPUTM
      IF(LCPUTM) WRITE (ND6,10020) TSOS,TMAX
      RETURN
10000 FORMAT(/1X,'$FILE:')
10010 FORMAT(1X,'SAVE FIELDS INTO FILE  (LSAV)  :  ',L1 /1X,
     +'RECOVER FROM FILE      (LREC)  :  ',L1 /1X,
     +'CPUTIME MONITOR ACTIVE (LCPUTM):  ',L1)
10020 FORMAT(1X,'CPUTIME LEFT FOR SAVING FIELDS : TSOS =',G13.5,'(S)'/ 1
     +X,'CPUTIME LIMIT OF BATCH JOB     : TMAX =',G13.5,'(S)')
      END
C
      SUBROUTINE PRBOUN
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &BOUN
C-----------------------------------------------------------------------
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      WRITE (ND6,10000)
      WRITE (ND6,10010)
      WRITE (ND6,10020)
      RETURN
10000 FORMAT(/1X,'$BOUN:')
10010 FORMAT(1X,'BOUNDARY CONDITIONS     LEFT   :' ,
     +' OPEN              ')
10020 FORMAT(1X,'BOUNDARY CONDITIONS     RIGHT  :' ,
     +' OPEN              ')
      END
C
      SUBROUTINE PRMESH(NW,LBEAMD,RAD,MTSTAB,RDZDR)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &MESH
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NVR=100)
      PARAMETER (NF=100)
      LOGICAL LCFRON,LCBACK,LCHIN,LBEAMD
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPL/ MROT
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /VAR/ RZ(NVR)
      IF(LBEAMD) GO TO 10
      WRITE (ND6,10000)
      WRITE (ND6,10010) NR
      WRITE (ND6,10020) NZ
      NFR=NZ/(NHIT/NR-4)+1
      IF(NFR.GE.10) WRITE(ND6,10090) NFR,NHIT
   10 NP=NR*(NW+4)
      IF(NP.GT.NR*NZ) NP=NR*NZ
      IF(.NOT.LCBACK) NP=NR*NZ
      RMP=FLOAT(NMESH)/FLOAT(NP)
      IF(.NOT.LBEAMD) THEN
        WRITE (ND6,10030) NP
        IF((RMP.LE.1.5).AND.LCBACK) WRITE(ND6,10040) RMP
      END IF
      IF(NP.GE.NMESH) THEN
        WRITE(ND6,10100) NP,NMESH
        STOP
      END IF
C     IF(LBEAMD) GO TO 300
      IF(LBEAMD) RETURN
      RDZDR=1.0
      IF(IMAX.EQ.1) THEN
        WRITE (ND6,10110) RAD,ZL
        WRITE (ND6,10050) DDR(1),DDZ
        RDZDR=AMAX1(RDZDR,DDZ/DDR(1))
      ELSE
        WRITE (ND6,10110) RAD,ZL
        RBOT=0.0
        DO 20 I=1,IMAX
          WRITE (ND6,10070) DDR(I),RBOT,RBDR(I)
          RDZDR=AMAX1(RDZDR,DDZ/DDR(I))
   20   RBOT=RBDR(I)
        WRITE (ND6,10060) DDZ
      END IF
      MTSTAB=INT(SQRT((1.2+0.4*FLOAT(MROT))*RDZDR**2+1.2))+1
      DO 30 I=1,NVR
        IF(RZ(I).EQ.-8888.) GO TO 30
        WRITE(ND6,10080) I,RZ(I)
   30 CONTINUE
      RETURN
10000 FORMAT(/1X,'$MESH:')
10010 FORMAT(1X,'NUMBER OF MESH LINES IN R      : NR   = ',I7)
10020 FORMAT(1X,'NUMBER OF MESH LINES IN Z      : NZ   = ',I7)
10030 FORMAT(1X,'NUMBER OF MESH POINTS          : NP   = ',I7)
10040 FORMAT(1X,'*** WARNING ***  THE RATIO OF NMESH TO NP:',G11.3,
     +' SHOULD BE MORE THAN 1.5 FOR A FASTER COMPUTATION')
10050 FORMAT(1X,'STEP SIZE IN R                 : DDR  =',G13.5,'(M)'/ 1
     +X,'STEP SIZE IN Z                 : DDZ  =',G13.5,'(M)')
10060 FORMAT(1X,'STEP SIZE IN Z                 : DDZ  =',G13.5,'(M)')
10070 FORMAT(1X,'STEP SIZE IN R                 : DDR  =',G13.5,'(M)',
     +'  FROM R=',G13.5,'  TO R=',G13.5)
10080 FORMAT(1X,'STRUCTURE INPUT VARIABLE       : RZ',I3, '=',G13.5,
     +'(M)')
10090 FORMAT(1X,'*** WARNING ***  BOUNDARY FRAMES ARE GENERATED:',I7,
     +' TIMES FOR NHIT:',I7)
10100 FORMAT(1X,'*** STOP ***  NP:',I7,' EXCEEDS MAXIMUM NMESH:',I7)
10110 FORMAT(1X,'TOTAL RADIUS OF THE STRUCTURE  : RAD  =',G13.5,'(M)' /1
     +X,'TOTAL LENGTH OF THE STRUCTURE  : ZL   =',G13.5,'(M)')
      END
C
      SUBROUTINE PRBEAM(ISIG)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &BEAM
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPL/ MROT
      COMMON /OFFS/ RDRIVE,ZCF,ZCT,RWAK
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      WRITE (ND6,10000)
      WRITE (ND6,10010) SIG
      IF((RDRIVE.GE.0.).AND.(ABS(RDRIVE-GR(JOFF)).GT.CEPS)) THEN
        WRITE (ND6,10030) ISIG,RDRIVE,GR(JOFF),MROT
      ELSE
        WRITE (ND6,10020) ISIG,GR(JOFF),MROT
      END IF
      IF(NBUNCH.GT.1) WRITE(ND6,10040) NBUNCH,BSEP
      IF(NBUNCH.LE.0) WRITE(ND6,10050)
      RETURN
10000 FORMAT(/1X,'$BEAM:')
10010 FORMAT(1X,'SIGMA OF THE GAUSSIAN BUNCH    : SIG  =',G13.5,'(M)')
10020 FORMAT(1X,'NUMBER OF STAND.DEV. USED      : ISIG = +/-',I3/ 1X,
     +'RADIAL BEAM OFFSET AT          :RDRIVE=',G13.5,'(M)'/ 1X,
     +'MODE NUMBER (MONOPOLE/DIPOLE)  : MROT =',I7)
10030 FORMAT(1X,'NUMBER OF STAND.DEV. USED      : ISIG = +/-',I3/ 1X,
     +'RADIAL BEAM OFFSET AT          :RDRIVE=',G13.5,'(M)', '  --->',
     +G13.5,'(M)'/ 1X,'ANGULAR MOMENTUM               : MROT =',I7)
10040 FORMAT(1X,'NUMBER OF GAUSSIAN BUNCHES     :NBUNCH=',I7/ 1X,
     +'SEPARATION BETWEEN BUNCHES     : BSEP =',G13.5,'(M)')
10050 FORMAT(1X,'*** ATTENTION *** USER NOW MUST PROVIDE THE CHARGE',
     +' DISTRIBUTION')
      END
C
      SUBROUTINE PRWAKE(UBTI,NW)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &WAKE
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      LOGICAL LCFRON,LCBACK,LCHIN,LNAPOLY,LCRBW,LNONAP
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /OFFS/ RDRIVE,ZCF,ZCT,RWAK
      COMMON /WPL/ MROT
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /PINT /KMIN
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      WRITE (ND6,10010)
      UBF=0.
      WRITE (ND6,10020) GR(KMIN)
      IF((LNAPOLY).AND.(MROT.EQ.0).AND.(KOFF.NE.KOFL)) WRITE
     +(ND6, 10030) GR(KOFL)
      WRITE (ND6,10040) UBF,UBTI
      WRITE (ND6,10050) LCFRON,LCBACK,LCHIN,LNAPOLY,LNONAP,LCRBW
      IF(LCRBW) WRITE(ND6,10055) ZSEP
      IF(LNAPOLY.AND.(ZCF.GE.0.)) WRITE (ND6,10060) ZCF
      IF(LNAPOLY.AND.(ZCT.GE.0.)) WRITE (ND6,10070) ZCT
      IF(LNAPOLY.AND.(RWAK.GE.0.)) WRITE (ND6,10080) RWAK
      WRITE (ND6,10090) NW
      IF(NW+2.GE.NWN) THEN
        WRITE(ND6,10000) NW+2,NWN
        STOP
      END IF
      RETURN
10000 FORMAT(1X,'*** STOP ***  NW+2:',I7,' EXCEEDS MAXIMUM NWN:',I7)
10010 FORMAT(/1X,'$WAKE:')
10020 FORMAT( 1X,'RADIAL WAKE OFFSET        (LEFT)   AT  R   =',G13.5,
     +'(M)')
10030 FORMAT( 1X,'                          (RIGHT)  AND R   =',G13.5,
     +'(M)')
10040 FORMAT( 1X,'WAKE    BETWEEN THE BUNCH  COORDINATES UBF =',G13.5,
     +'(M)'/ 1X,'                                   AND UBT =',G13.5,
     +'(M)')
10050 FORMAT(1X,'WINDOW FOR FRONT                 (LCFRON)  :  ',L1/ 1X,
     +'WINDOW FOR BACK                  (LCBACK)  :  ',L1/ 1X,
     +'CHIN WAKE INTEGRATION METHOD     (LCHIN )  :  ',L1/ 1X,
     +'NAPOLY WAKE INTEGRATION METHOD   (LNAPOLY) :  ',L1/ 1X,
     +'CONVENTIONAL INTEGRATION METHOD  (LNONAP)  :  ',L1/ 1X,
     +'WAKE FOR A COUNTER-ROTATING BEAM (LCRBW)   :  ',L1)
10055 FORMAT(1X,'INITIAL DISTANCE BTW. TWO BEAMS: ZSEP =',G13.5,'(M)')
10060 FORMAT(1X,'STARTING Z FOR INTEGRAL AT RWAK: ZCF  =',G13.5,'(M)')
10070 FORMAT(1X,'LAST Z FOR INTEGRAL AT RWAK    : ZCT  =',G13.5,'(M)')
10080 FORMAT(1X,'RADIAL OFFSET FOR WAKE INTEGRAL: RWAK =',G13.5,'(M)')
10090 FORMAT(1X,'NUMBER OF WAKE POTENTIAL POINTS: NW   =',I8)
      END
C
      SUBROUTINE PRTIME(MTSTAB,RDZDR,LPLE,LPLC,TPS,NSHOT)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &TIME
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NPT=20)
      LOGICAL LPLE,LPLC
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WIND/ NT
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      BETA=1.
      DT=DDZ/(CLGHT*FLOAT(M))*1.E12
      PT=FLOAT(NZ*M)*DT
      WRITE (ND6,10000) NT,M
      IF(M.GE.NPT) THEN
        WRITE(ND6,10010) M+1,NPT
        STOP
      END IF
      IF((M.LT.MTSTAB).AND.(RDZDR.GT.1.)) THEN
        WRITE(ND6,10020) RDZDR,MTSTAB
      END IF
      WRITE (ND6,10030) BETA,DT,PT
      IF(LPLE.OR.LPLC) WRITE(ND6,10040) TPS,NSHOT
      RETURN
10000 FORMAT(/1X,'$TIME:'/ 1X,'TIME STEPS TO BE PROCESSED     : NT   =',
     +I7/ 1X,'NUMBER OF TIME STEPS/MESH STEP : MT   =',I7)
10010 FORMAT(1X,'*** STOP ***  MT+1:',I7,' EXCEEDS MAXIMUM NPT:',I7)
10020 FORMAT(1X,'*** WARNING ***  NO GUARANTEE FOR NUMERICAL STABILITY',
     +' IN THE PRESENT CALCULATION WHERE MAX(DDZ/DDR) IS',G13.5/ 1X,
     +'*** A RECOMMENDED VALUE FOR MT : MT   =',I7,' ***')
10030 FORMAT(1X,'VELOCITY OF THE BUNCH / C      : BETA =',G13.5/ 1X,
     +'TIME-STEP VALUE                : DT   =',G13.5,'(PS)'/ 1X,
     +'TIME FOR A PARTICLE TO PASS    : PT   =',G13.5,'(PS)')
10040 FORMAT(1X,'START TIME FOR FILED LINE PLOT : TPS  =',G13.5,'(PS)'/
     +1X,'NUMBER OF FILED LINE PLOTS     : NSHOT=',I7)
      END
C
      SUBROUTINE PRPLOT(LCAVIN,LCAVUS,LPLW,LFFT,LFFTA,LFFTT,LFFTL,LSPEC,
     +LPLE,LPLC,LPALL,NPLOT,CUTOFF,LFFTOR,REDUCE)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &PLOT
C-----------------------------------------------------------------------
      LOGICAL LCAVIN,LCAVUS,LPLW,LFFT,LFFTA,LFFTT,LFFTL,LSPEC,
     +LPLE,LPLC,LPALL,LFFTOR,LWNDW,LPATH
      LOGICAL LPLWA,LPLWT,LPLWL,LSVWA,LSVWT,LSVWL,LINTZ,LCRBW
      LOGICAL LCBACK,LCFRON,LCHIN,LSVF
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /WPL/ MROT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /FFTS/ LSVF,LWNDW
      COMMON /WFUN/ NWFUN,ALPHA
      COMMON /EXPF/ EXPFAC
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      COMMON /WPS/ LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
      COMMON /INTR/ LINTZ
      COMMON /PLNZ/ LPATH
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      WRITE(ND6,10060)
      WRITE(ND6,10000) LCAVIN,LCAVUS,LPATH,LPLW,LPLWA,LPLWT,LPLWL,LFFT,
     +LFFTA,LFFTT,LFFTL,LINTZ,LSPEC
      IF(LSPEC.AND.LCRBW) WRITE(ND6,10002)
      WRITE(ND6,10004) LPLE,LPLC,LPALL,NPLOT,EXPFAC,CUTOFF,
     +LWNDW,NWFUN,ALPHA
      IF(MROT.EQ.0) LFFTOR=LFFT.OR.LFFTL.OR.LSPEC
      IF(MROT.EQ.1) LFFTOR=LFFT.OR.LFFTL.OR.LFFTA.OR.LFFTT.OR.LSPEC
      IF(LFFTOR.AND.(NBUNCH.GT.1)) WRITE(ND6,10050)
      IF(LPLE.AND.LCBACK) WRITE(ND6,10010)
      IF(LPLC.AND.LCBACK) WRITE(ND6,10020)
      IF(LPLE.AND.(MROT.NE.0)) WRITE(ND6,10030)
      IF(LPLC.AND.(MROT.NE.0)) WRITE(ND6,10040)
      IF(REDUCE.NE.0.0) WRITE(ND6,10001) REDUCE
      RETURN
10001 FORMAT(1X,'REDUCTION FACTOR OF ALL PLOTS  :REDUCE=',G13.5)
10000 FORMAT(1X,'PLOT OF CAVITY SHAPE INPUT       (LCAVIN)  :  ',L1/ 1X,
     +'PLOT OF CAVITY SHAPE USED        (LCAVUS)  :  ',L1/ 1X,
     +'PLOT OF NAPOLY-ZOTTER CONTOUR    (LPATH)   :  ',L1/ 1X,
     +'PLOT OF ALL WAKE POTENTIALS      (LPLW)    :  ',L1/ 1X,
     +'PLOT OF AZIMUTHAL WAKE POTENTIAL (LPLWA)   :  ',L1/ 1X,
     +'PLOT OF TRANSVERSE WAKE POTENT.  (LPLWT)   :  ',L1/ 1X,
     +'PLOT OF LONGITUDINAL WAKE POTENT.(LPLWL)   :  ',L1/ 1X,
     +'PLOT OF FFT OF WAKE POTENTIALS   (LFFT)    :  ',L1/ 1X,
     +'PLOT OF FFT OF AZIMUTHAL WAKE    (LFFTA)   :  ',L1/ 1X,
     +'PLOT OF FFT OF TRANSVERSE WAKE   (LFFTT)   :  ',L1/ 1X,
     +'PLOT OF FFT OF LONGITUDINAL WAKE (LFFTL)   :  ',L1/ 1X,
     +'PLOT OF INTEGRATION OF IMPEDANCE (LINTZ)   :  ',L1/ 1X,
     +'PLOT OF ENERGY LOSS SPECTRUM     (LSPEC)   :  ',L1)
10002 FORMAT(/1X,'*** SKIP ***  LSPEC CANNOT BE TRUE WHEN LCRBW IS ',
     +'TRUE IN THIS VERSION OF ABCI.'/)
10004 FORMAT(1X,'PLOT OF ELECTRIC FIELD LINES     (LPLE)    :  ',L1/ 1X,
     +'PLOT OF TOTAL CURRENT LINES      (LPLC)    :  ',L1/ 1X,
     +'SWITCH FROM ABSOLUTE TO RELATVE  (LPALL)   :  ',L1/ 1X,
     +'NUMBER OF CONTOUR LINES        : NPLOT=',I8/ 1X,
     +'ENHANCEMENT FACTOR FOR FFT PLOT:EXPFAC=',G13.5/ 1X,
     +'CUTOFF FREQUENCY FOR FFT PLOT  :CUTOFF=',G13.5, '(GHZ)'/ 1X,
     +'WINDOWING FOR FOURIER TRANSFORMS (LWNDW)   :  ',L1/ 1X,
     +'WINDOW FUNCTION                : NWFUN=',I8/ 1X,
     +'WINDOW PARAMETER FOR WIDTH     : ALPHA=',G13.5)
10010 FORMAT(1X,'*** WARNING ***  WHEN LPLE=.T., LCBACK MUST TO BE .F.',
     +'   NO PLOT OF ELECTRIC FIELD LINES WILL BE PRODUCED.')
10020 FORMAT(1X,'*** WARNING ***  WHEN LPLC=.T., LCBACK MUST TO BE .F.',
     +'   NO PLOT OF ELECTRIC FIELD LINES WILL BE PRODUCED.')
10030 FORMAT(1X,'*** WARNING ***  WHEN LPLE=.T., MROT MUST TO BE 0.',
     +'   NO PLOT OF ELECTRIC FIELD LINES WILL BE PRODUCED.')
10040 FORMAT(1X,'*** WARNING ***  WHEN LPLC=.T., MROT MUST TO BE 0.',
     +'   NO PLOT OF ELECTRIC FIELD LINES WILL BE PRODUCED.')
10050 FORMAT(1X,'*** WARNING ***  WHEN NBUNCH IS GREATER THAN 1,',
     +' NO FFT WILL BE PRODUCED.')
10060 FORMAT(/1X,'$PLOT:')
      END
C
      SUBROUTINE PRPRIN(LMATPR,LPRW,LSVW)
C-----------------------------------------------------------------------
C PRINT INPUT IN NAMELIST &PRIN
C-----------------------------------------------------------------------
      LOGICAL LPLWA,LPLWT,LPLWL,LSVWA,LSVWT,LSVWL
      LOGICAL LWNDW,LMATPR,LPRW,LSVW,LSVF
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /FFTS/ LSVF,LWNDW
      COMMON /WPS/ LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
      WRITE(ND6,10010)
      WRITE(ND6,10000) LMATPR,LPRW,LSVW,LSVWA,LSVWT,LSVWL,LSVF
      RETURN
10000 FORMAT(1X,'PRINTOUT OF CAVITY SHAPE USED    (LMATPR)  :  ',L1/ 1X,
     +'PRINTOUT OF WAKE POTENTIALS      (LPRW)    :  ',L1/ 1X,
     +'SAVE WAKE POTENTIALS IN A FILE   (LSVW)    :  ',L1/ 1X,
     +'SAVE AZIMUTHAL WAKE IN A FILE    (LSVWA)   :  ',L1/ 1X,
     +'SAVE TRANSVERSE WAKE IN A FILE   (LSVWT)   :  ',L1/ 1X,
     +'SAVE LONGITUDINAL WAKE IN A FILE (LSVWL)   :  ',L1/ 1X,
     +'SAVE FFT RESULTS IN A FILE       (LSVF)    :  ',L1)
10010 FORMAT(/1X,'$PRIN:')
      END
C
      SUBROUTINE WKPLOT(NW,ZTXT,LPLW,LSVW,LPOT,LNTR,NSTK,AXY,ITEST)
C-----------------------------------------------------------------------
C PLOT WAKE POTENTIALS ON LINEPRINTER
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NWN=500000)
      PARAMETER (NF=100)
      LOGICAL LPLW,LSVW,LPOT,LNTR,LCHIN,LNAPOLY,LCFRON,LCBACK,LNONAP
      LOGICAL LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT,LPLWOR,LSVWOR
      LOGICAL LPLWLO,LSVWLO
      CHARACTER*1 LZEIH,LZEIL,LZEIB,LZEIV,LZEIP
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      DIMENSION QXYMM(4),EFI(NWN)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPL/ MROT
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /WPT/ NDATE,NTIME
      COMMON /WP2/ TT
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /WPS/ LPLWL,LPLWA,LPLWT,LSVWL,LSVWA,LSVWT
      DATA LZEIH,LZEIL,LZEIB,LZEIV,LZEIP/'T','L','*','A','P'/
      DATA PC/1.E-12/
      Q1MAX=0.
      Q1MIN=0.
      Q2MAX=0.
      Q2MIN=0.
      Q3MAX=0.
      Q3MIN=0.
      Q4MAX=0.
      Q4MIN=0.
      DO 10 IS=0,NW
        SUMZ=WZ(IS+1)
        SUMT=WT(IS+1)
        SUMY=WY(IS+1)
        SUMP=WP(IS+1)
        Q1MAX=AMAX1(SUMZ,Q1MAX)
        Q1MIN=AMIN1(SUMZ,Q1MIN)
        Q2MAX=AMAX1(SUMT,Q2MAX)
        Q2MIN=AMIN1(SUMT,Q2MIN)
        Q3MAX=AMAX1(SUMY,Q3MAX)
        Q3MIN=AMIN1(SUMY,Q3MIN)
        Q4MAX=AMAX1(SUMP,Q4MAX)
        Q4MIN=AMIN1(SUMP,Q4MIN)
   10 CONTINUE
      ZSUM=AMAX1(ABS(Q1MAX),ABS(Q1MIN))
      TSUM=AMAX1(ABS(Q2MAX),ABS(Q2MIN))
      YSUM=AMAX1(ABS(Q3MAX),ABS(Q3MIN))
      PSUM=AMAX1(ABS(Q4MAX),ABS(Q4MIN))
C-----------------------------------------------------------------------
C LOSS PARAMETER PK
C-----------------------------------------------------------------------
      PKZ=0.
      PKT=0.
      PKY=0.
      PKP=0.
      DO 20 N=1,IBT
        CALL CHAGE(1,(N-1)*M,QK)
        PKZ=PKZ+WZ(N)*QK*DDZ
        PKT=PKT+WT(N)*QK*DDZ
        PKY=PKY+WY(N)*QK*DDZ
        PKP=PKP+WP(N)*QK*DDZ
   20 CONTINUE
      PKH=0.
      IBTH=IBT/2
      DO 25 N=1,IBTH-1
        CALL CHAGE(1,(N-1)*M,QK)
        PKH=PKH+WZ(N)*QK*DDZ
   25 CONTINUE
        CALL CHAGE(1,(IBTH-1)*M,QK)
        PKH=PKH+WZ(IBTH)*QK*DDZ/2.
        PKL=PKZ-PKH
C-----------------------------------------------------------------------
C PRINTS THE WAKE POTENTIAL AND GENERATES A CURVE VIA DRPLOT
C-----------------------------------------------------------------------
      IPRI=1
      QXYMM(1)=0.
      QXYMM(2)=FLOAT(NW-1)*DDZ
      QXYMM(3)=-1.00
      QXYMM(4)=+1.00
      QTOT=1.
C-----------------------------------------------------------------------
C GET CPU TIME USED FOR THIS RUN
C-----------------------------------------------------------------------
      CALL CPUTIM(TT,2)
C-----------------------------------------------------------------------
      WRITE (ND6,10000) (ZTXT(II),II=1,9),NDATE,NTIME,SIG,MROT,TT, DDR
     +(1),DDZ,LCHIN,LNAPOLY
C-----------------------------------------------------------------------
C WRITE HEAD PART OF TOPDRAW INPUT FILE
C-----------------------------------------------------------------------
      IF(LPLW) LPLWL=.TRUE.
      IF(LPLW) LPLWA=.TRUE.
      IF(LPLW) LPLWT=.TRUE.
      IF(LSVW) LSVWL=.TRUE.
      IF(LSVW) LSVWA=.TRUE.
      IF(LSVW) LSVWT=.TRUE.
      LPLWOR=LPLWL.OR.LPLWA.OR.LPLWT
      LSVWOR=LSVWL.OR.LSVWA.OR.LSVWT
      LPLWLO=(MROT.EQ.0).AND.(LPLWOR)
      LSVWLO=(MROT.EQ.0).AND.(LSVWOR)
C-----------------------------------------------------------------------
      IF(LPLWOR) CALL WPLOT(ZTXT,EFI,NW,0,QXYMM)
      IF(LPLWOR) CALL WPADD(QMIN,QMAX,PKY,MROT,0)
      IF(JOFF.GT.1) WFT=1./GR(JOFF)*PC
C-----------------------------------------------------------------------
C NOW PRINT LINE CHARGE DENSITY
      IF(ITEST.EQ.1) IPRI=0
      QKMAX=0.
      QKMIN=1.E15
      DO 30 I=1,NW
        CALL CHAGE(1,(I-1)*M,QK)
        QKMAX=AMAX1(QK*DDZ,QKMAX)
        QKMIN=AMIN1(QK*DDZ,QKMIN)
        EFI(I)=QK*DDZ
   30 CONTINUE
      DO 40 I=1,NW
   40 EFI(I)=EFI(I)/QKMAX
      CALL DRPLOT (ND6,DUMMY,EFI,NW,LZEIB,0,QXYMM,IPRI)
      IF(LPLWOR) CALL WPLOT(ZTXT,EFI,NW,1,QXYMM)
C-----------------------------------------------------------------------
C NOW WAKE AZIMUTHAL
      IF(ITEST.EQ.1) GO TO 130
   50 IF(YSUM.LT.ZSUM*1.E-10) GO TO 70
      DO 60 I=1,NW
   60 EFI(I)=WY(I)/YSUM
      CALL DRPLOT (ND6,DUMMY,EFI,NW,LZEIV,0,QXYMM,IPRI)
      IF(LPLWA) CALL WPLOT(ZTXT,EFI,NW,2,QXYMM)
      IF(LSVWA) CALL WPRIN(ZTXT,WY,NW,2,WFT)
C-----------------------------------------------------------------------
C NOW WAKE TRANSVERSE
   70 IF(ITEST.EQ.1) GO TO 130
   80 IF(TSUM.LT.ZSUM*1.E-10) GO TO 100
      DO 90 I=1,NW
   90 EFI(I)=WT(I)/TSUM
      CALL DRPLOT (ND6,DUMMY,EFI,NW,LZEIH,0,QXYMM,IPRI)
      IF(LPLWT) CALL WPLOT(ZTXT,EFI,NW,3,QXYMM)
      IF(LSVWT) CALL WPRIN(ZTXT,WT,NW,3,WFT)
C-----------------------------------------------------------------------
C NOW WAKE LONGITUDINAL
  100 CONTINUE
      IF(.NOT.LPOT) IPRI=0
      DO 110 I=1,NW
  110 EFI(I)=WZ(I)/ZSUM
      CALL DRPLOT (ND6,DUMMY,EFI,NW,LZEIL,0,QXYMM,IPRI)
      IF(LPLWL.OR.LPLWLO) CALL WPLOT(ZTXT,EFI,NW,4,QXYMM)
      IF((LSVWL.AND.MROT.EQ.0).OR.LSVWLO) CALL WPRIN(ZTXT,WZ,NW,1,PC)
      IF(LSVWL.AND.MROT.EQ.1) CALL WPRIN(ZTXT,WZ,NW,4,WFT/AXY)
      IF(LPOT) THEN
        IPRI=0
        DO 120 I=1,NW
  120   EFI(I)=WP(I)/PSUM
        CALL DRPLOT (ND6,DUMMY,EFI,NW,LZEIP,0,QXYMM,IPRI)
        IF(LPLWL.OR.LPLWLO) CALL WPLOT(ZTXT,EFI,NW,5,QXYMM)
        IF(LSVWL.OR.LSVWLO) CALL WPRIN(ZTXT,WP,NW,5,PC)
      ELSE
      END IF
C-->
  130 CONTINUE
      IF(JOFF.GT.1) YFT=1./GR(JOFF)*PC
      ZFT=YFT/AXY
      GRF=GR(JOFF)
      WRITE(ND6,10090) QKMIN/DDZ,QKMAX/DDZ,QKMAX/DDZ,GRF
      IF(.NOT.LNTR) THEN
        IF(MROT.EQ.1) WRITE(ND6,10100) Q3MIN*YFT,Q3MAX*YFT,YSUM*YFT,GRF
        IF(MROT.EQ.1) WRITE(ND6,10110) Q2MIN*YFT,Q2MAX*YFT,TSUM*YFT,GRF
        IF(MROT.EQ.1) WRITE(ND6,10140) Q1MIN*ZFT,Q1MAX*ZFT,ZSUM*ZFT,GRF
        IF(MROT.EQ.0) WRITE(ND6,10120) Q1MIN*PC,Q1MAX*PC,ZSUM*PC
        IF(LPOT) WRITE(ND6,10130) Q4MIN*PC,Q4MAX*PC,PSUM*PC
      ELSE
        IF(MROT.EQ.1) WRITE(ND6,10110) Q2MIN*YFT,Q2MAX*YFT,TSUM*YFT,GRF
        IF(MROT.EQ.1) WRITE(ND6,10140) Q1MIN*YFT,Q1MAX*YFT,ZSUM*ZFT,GRF
      END IF
      IF(((NSTK.LT.NT-1).AND.(ITEST.EQ.0)).OR.(SIG/DDZ.LT.4.0)) THEN
        WRITE(ND6,10160)
      ELSE
        WRITE(ND6,10150)
      END IF
      IF(.NOT.LNTR) THEN
        IF(MROT.EQ.1) WRITE(ND6,10040) PKY*YFT
        IF(MROT.EQ.1) WRITE(ND6,10050) PKT*YFT
        IF(MROT.EQ.1) WRITE(ND6,10010) PKZ*ZFT
        IF(MROT.EQ.0) WRITE(ND6,10020) PKZ*PC, PKH*PC, PKL*PC
        IF(LPLWOR) THEN
          IF(MROT.EQ.1) CALL WPADD(Q3MIN*YFT,Q3MAX*YFT,PKY*YFT,MROT,2)
          IF(MROT.EQ.1) CALL WPADD(Q2MIN*YFT,Q2MAX*YFT,PKT*YFT,MROT,3)
          IF(MROT.EQ.1) CALL WPADD(Q1MIN*ZFT,Q1MAX*ZFT,PKZ*ZFT,MROT,4)
          IF(MROT.EQ.0) CALL WPADD(Q1MIN*PC,Q1MAX*PC,PKZ*PC,MROT,4)
        ELSE
        END IF
        IF(LPOT) WRITE(ND6,10030) PKP*PC
        IF(LPOT.AND.LPLWOR) CALL WPADD(Q4MIN*PC,Q4MAX*PC,PKP*PC,MROT,5)
      ELSE
        IF(MROT.EQ.1) WRITE(ND6,10050) PKT*YFT
        IF(MROT.EQ.1) WRITE(ND6,10010) PKZ*ZFT
        IF((MROT.EQ.1).AND.LPLWOR) CALL WPADD(Q2MIN*YFT,Q2MAX*YFT,
     +  PKT*YFT,MROT,3)
        IF((MROT.EQ.1).AND.LPLWOR) CALL WPADD(Q1MIN*ZFT,Q1MAX*ZFT,
     +  PKZ*ZFT,MROT,4)
      END IF
      IF((NSTK.LT.NT-1).AND.(ITEST.EQ.0)) WRITE(ND6,10060) NSTK
      IF(SIG/DDZ+CEPS.LT.4.0) WRITE(ND6,10070) SIG/DDZ
      RETURN
10000 FORMAT(1H1,'PROBLEM :',9(A8),2X,A8,1X,A8,5X,'SIG/M=',1PE10.3, 2X,
     +'MROT=',I2/ 1X,71X,'CPUTIME USED:',1PE10.3,'(S)',3X,'DDR=',1PE10.
     +3, 2X,'DDZ=',1PE10.3/ 1X,
     +'(A,T,L,P) = WAKE POTENTIAL AS A FUNCTION OF',
     +' THE BUNCH COORDINATES S IN (M), (*)=BUNCH SHAPE,  FRONT LEFT', 5
     +X,'LCHIN= ',L1,3X,'LNAPOLY= ',L1)
10010 FORMAT(1X,'INTEGRATED  LONGITUDINAL WAKE * ', 'CHARGE DENSITY = ',
     +1PE10.3,' V/PC/M**2')
10020 FORMAT(1X,'INTEGRATED  LONGITUDINAL WAKE * ', 'CHARGE DENSITY = ',
     +1PE10.3,' V/PC/M**0',3X,
     +'(FIRST/LAST HALF = ',1PE10.3,'/',1PE10.3,' V/PC/M**0)')
10030 FORMAT(1X,'WITH LOG. TERM DUE TO UNEQUAL BE', 'AM PIPE RADII  = ',
     +1PE10.3,' V/PC/M**0')
10040 FORMAT(1X,'INTEGRATED  AZIMUTHAL    WAKE * ', 'CHARGE DENSITY = ',
     +1PE10.3,' V/PC/M**1')
10050 FORMAT(1X,'INTEGRATED  TRANSVERSE   WAKE * ', 'CHARGE DENSITY = ',
     +1PE10.3,' V/PC/M**1')
10060 FORMAT(/1X,'*** WARNING *** COMPUTATION IS NOT COMPLETED.',
     +'   THE LAST TIME STEP =',I7)
10070 FORMAT(/1X,'*** WARNING *** YOU SHOULD NOT BELIEVE THE ABOVE',
     +' RESULTS.  SIG/DDZ=',F7.2,' IS TOO SMALL (IT SHOULD BE AT LEAST',
     +' 4-5).')
10080 FORMAT(/' TOTAL CHARGE IN THE BUNCH =',1PE10.3,' AS')
10090 FORMAT(/' LINE CHARGE DENSITY MIN/MAX=',1PE10.3,'/',1PE10.3,
     +'    C/M**1,',
     +' SCALE=',1PE10.3,'    C/M**1,',
     +' PASSING    AT R=',1PE10.3,' M')
10100 FORMAT( ' AZIMUTHAL    WAKE   MIN/MAX=',1PE10.3,'/',
     +1PE10.3,' V/PC/M**1,',' SCALE=',1PE10.3,' V/PC/M**1,',
     +' NORMALIZED BY D=',1PE10.3,' M DIPOLE OFFSET')
10110 FORMAT( ' TRANSVERSE   WAKE   MIN/MAX=',1PE10.3,'/',
     +1PE10.3,' V/PC/M**1,',' SCALE=',1PE10.3,' V/PC/M**1,',
     +' NORMALIZED BY D=',1PE10.3,' M DIPOLE OFFSET')
10120 FORMAT( ' LONGITUDINAL WAKE   MIN/MAX=',1PE10.3,'/',
     +1PE10.3,' V/PC/M**0,',' SCALE=',1PE10.3,' V/PC/M**0')
10130 FORMAT( ' WAKE WITH LOG, TERM MIN/MAX=',1PE10.3,'/',
     +1PE10.3,' V/PC/M**0,',' SCALE=',1PE10.3,' V/PC/M**0')
10140 FORMAT( ' LONGITUDINAL WAKE   MIN/MAX=',1PE10.3,'/',
     +1PE10.3,' V/PC/M**2,',' SCALE=',1PE10.3,' V/PC/M**2,',
     +' NORMALIZED BY D=',1PE10.3,' M DIPOLE OFFSET')
10150 FORMAT(1H0)
10160 FORMAT(1X)
      END
C
      SUBROUTINE EXAMIN (NLOOP,LBEAMD,LNTR,RMIN,RVAC,ZLFT,ZRGT,I1040)
C-----------------------------------------------------------------------
C CHECK GEOMETRY AND WHETHER NAPOLY'S METHOD IS NEEDED
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NWN=500000)
      PARAMETER (NF=100)
      LOGICAL LNAPOLY,LNTR,LBEAMD,LPIPER,LPIPEL,LNONAP
      CHARACTER*2 ORD(3),ORDX
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPL/ MROT
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /OFFS/ RDRIVE,ZCF,ZCT,RWAK
      COMMON /PIPE/ KOFF,KOFL
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /ICON/ IZLFT,IZRGT,IRMIN
      DATA ORD/'ND','RD','TH'/
C-----------------------------------------------------------------------
C CHECK IF BEAM PIPES TOO SHORT
C-----------------------------------------------------------------------
      I1040=0
      LPIPEL=.FALSE.
      IZDEF=IZZ(3)-IZZ(2)
      IF((IRR(2).EQ.IRR(3)).AND.(IZDEF.GT.4)) LPIPEL=.TRUE.
      LPIPER=.FALSE.
      IZDEL=IZZ(NVC)-IZZ(NVC-1)+1
      IF((IRR(NVC-1).EQ.IRR(NVC)).AND.(IZDEL.GT.4))  LPIPER=.TRUE.
      IF((.NOT.LPIPEL).OR.(.NOT.LPIPER)) THEN
        IF((NLOOP.GT.0).AND.(LBEAMD)) THEN
          IF(NLOOP.LE.2) ORDX=ORD(NLOOP)
          IF(NLOOP.GT.2) ORDX=ORD(3)
          CALL DATIME(NDATE,NTIME)
          WRITE (ND6,10000) NLOOP+1,ORDX,NDATE,NTIME
        ELSE
        END IF
        WRITE(ND6,10010) IZDEF, IZDEL
        STOP
      ELSE
      END IF
C-----------------------------------------------------------------------
C CHECK IF BEAM PIPE RADII ARE EQUAL FOR MROT=1
C-----------------------------------------------------------------------
      IF((MROT.EQ.1).AND.(KOFF.NE.KOFL)) THEN
        IF((NLOOP.GT.0).AND.(LBEAMD)) THEN
          IF(NLOOP.LE.2) ORDX=ORD(NLOOP)
          IF(NLOOP.GT.2) ORDX=ORD(3)
          CALL DATIME(NDATE,NTIME)
          WRITE (ND6,10000) NLOOP+1,ORDX,NDATE,NTIME
        ELSE
        END IF
        IF(.NOT.LNONAP) THEN 
         WRITE(ND6,10020) GR(KOFF),GR(KOFL)
C     STOP
         I1040=1040
         RETURN
        ELSE
         WRITE(ND6,10025) GR(KOFF),GR(KOFL)
        END IF
      ELSE
      END IF
C-----------------------------------------------------------------------
C NAPOLY INTEGRATION FOR MROT=1 AND 0
C-----------------------------------------------------------------------
      LNTR=.FALSE.
C--
      IF((LNAPOLY).AND.(MROT.EQ.1)) LNTR=.TRUE.
      IF(LNAPOLY) THEN
        IF(ZCF.GE.0.) IZCF=NINT(ZCF/DDZ)+1
        IF(ZCT.GE.0.) IZCT=NINT(ZCT/DDZ)+1
        IF(RMIN.GE.RVAC) THEN
          IF(ZCF.LT.0.) IZCF=IZZ(3)-2
          IF(ZCT.LT.0.) IZCT=IZZ(NVC-1)+2
        ELSE
          IF(ZCF.LT.0.) IZCF=IZLFT-2
          IF(ZCT.LT.0.) IZCT=IZRGT+2
        END IF
        IF(IZCF.LT.2) IZCF=2
        IF(IZCT.GT.NZ) IZCT=NZ
        IF((ZCF.LT.0.).OR.(ZCT.LT.0.).OR.(RWAK.LT.0.)) THEN
          IF((NLOOP.GT.0).AND.(LBEAMD)) THEN
            IF(NLOOP.LE.2) ORDX=ORD(NLOOP)
            IF(NLOOP.GT.2) ORDX=ORD(3)
            CALL DATIME(NDATE,NTIME)
            WRITE (ND6,10000) NLOOP+1,ORDX,NDATE,NTIME
            LBEAMD=.FALSE.
          ELSE
          END IF
          WRITE(ND6,10030) FLOAT(IZCF-1)*DDZ,FLOAT(IZCT-1)*DDZ,GR(IRWAK)
        ELSE
        END IF
      ELSE
      END IF
C-----------------------------------------------------------------------
C CHECK IF NAPOLY'S INTEGRATION OPTION IS CHOSEN
C-----------------------------------------------------------------------
      IF((RMIN.LT.RVAC).OR.(KOFF.NE.KOFL)) THEN
        IF((NLOOP.GT.0).AND.(LBEAMD)) THEN
          IF(NLOOP.LE.2) ORDX=ORD(NLOOP)
          IF(NLOOP.GT.2) ORDX=ORD(3)
          CALL DATIME(NDATE,NTIME)
          WRITE (ND6,10000) NLOOP+1,ORDX,NDATE,NTIME
          LBEAMD=.FALSE.
        ELSE
        END IF
        IF(.NOT.LNONAP) WRITE(ND6,10040) LNAPOLY
        IF(LNAPOLY) THEN
          IF(RMIN.LT.RVAC) THEN
            IF((IZCF.GE.IZLFT).OR.(IZCT.LE.IZRGT).OR.(IRWAK.GT.IRMIN))
     +      THEN
              WRITE(ND6,10050) FLOAT(IZCF-1)*DDZ,ZLFT,FLOAT(IZCT-1)*DDZ,
     +        ZRGT,GR(IRWAK),RMIN
              IZCF=IZLFT-2
              IF(IZCF.LT.2) IZCF=2
              IZCT=IZRGT+2
              IF(IZCT.GT.NZ) IZCT=NZ
              IRWAK=IRMIN
              WRITE(ND6,10060) FLOAT(IZCF-1)*DDZ,FLOAT(IZCT-1)*DDZ,GR
     +        (IRWAK)
C
            ELSE
            END IF
          ELSE
          END IF
        ELSE
          IZCF=IZLFT-2
          IF(IZCF.LT.2) IZCF=2
          IZCT=IZRGT+2
          IF(IZCT.GT.NZ) IZCT=NZ
          IF(RMIN.LT.RVAC) THEN
            IRWAK=IRMIN
          ELSE
            IRWAK=MIN(KOFF,KOFL)
          END IF
          IF(.NOT.LNONAP) THEN 
            LNAPOLY=.TRUE.
            IF(MROT.EQ.1)  LNTR=.TRUE.
            WRITE(ND6,10070) FLOAT(IZCF-1)*DDZ,FLOAT(IZCT-1)*DDZ,
     +      GR(IRWAK)
          ELSE
C           LNTR=.TRUE.
            WRITE(ND6,10075) GR(IRWAK)
          END IF
        END IF
      ELSE
      END IF
C-----------------------------------------------------------------------
C IRWAK TOO CLOSE TO THE BEAM AZIS FOR MRTO=1 ?
C-----------------------------------------------------------------------
      IF(LNTR)  THEN
        IF(IRWAK.LE.4) THEN
          WRITE(ND6,10080) IRWAK
          I1040=1040
        ELSE
        END IF
      ELSE
      END IF
C-----------------------------------------------------------------------
C SET UP NAPOLY'S INTEGRAL PATH
C-----------------------------------------------------------------------
      IF(LNAPOLY) THEN
        SGF=1.0
        IF(IRWAK.GT.KOFF) SGF=-1.0
        KENF=KOFF-1
        IF(IRWAK.GT.KOFF) KENF=IRWAK-1
        KSTF=IRWAK
        IF(IRWAK.GT.KOFF) KSTF=KOFF
        SGL=1.0
        IF(IRWAK.GT.KOFL) SGL=-1.0
        KENL=KOFL-1
        IF(IRWAK.GT.KOFL) KENL=IRWAK-1
        KSTL=IRWAK
        IF(IRWAK.GT.KOFL) KSTL=KOFL
        ARM=GR(KOFF)/GR(IRWAK)
        RAM=1./ARM
        ARA=(ARM+RAM)
        ARD=(ARM-RAM)
      ELSE
      END IF
      RETURN
10000 FORMAT(1H1,76(1H*) /1X,2H* ,72(1H ),2H * /1X,2H* ,
     +'                           ABCI  ',I2,A2,' RUN    ' ,27(1H ),2H *
     +/1X,2H* ,72(1H ),2H * /1X,2H* ,'  DATE:',A8,'   TIME:',A8, 3X,
     +'  VERSION 9.4,  JANUARY  1995  ',6X,1(1H ),2H * /1X,2H* ,72(1H ),
C    2 /1X,2H* ,30(1H ),I2,'-TH RUN',33(1H ),2H *
     +2H * /1X,76(1H*)/)
10010 FORMAT(1H0,'*** STOP *** THE BEAM PIPES AT BOTH ENDS ARE',
     +' TOO SHORT.'/ 1X,
     +'THEY MUST HAVE AT LEAST 5 MESH LENGTH.  THEIR LENGTHES ARE:'/ 1X,
     +'L(LEFT)/DDZ=',I7,'   AND L(RIGHT)/DDZ=',I7/ 1X,
     +'INCREASE THE BEAM PIPE LENGTH AND RUN IT AGAIN.')
10020 FORMAT(1H0,'*** SKIP THIS CALCULATION *** UNEQUAL BEAM PIPE',
     +' RADII AT BOTH ENDS'/ 1X,'R=',G13.5,'(M) AND R=',G13.5,
     +'(M) ARE NOT ALLOWED FOR MROT=1.'/ 1X,'*** SET LNONAP=T ',
     +'IF YOU WANT TO COMPUTE WITH CONVENTIONAL INTEGRATION ***')
10025 FORMAT(1H0,'*** CAUTION *** UNEQUAL BEAM PIPE RADII AT BOTH ENDS'
     +/1X,'R=',G13.5,'(M) AND R=',G13.5,'(M) FOR MROT=1 CASE')
10030 FORMAT(1H ,'*** WARNING *** YOU DID NOT FULLY SPECIFY THE',
     +' INTEGRATION CONTOUR DESPITE OF LNAPOLY=T.'/ 1H ,
     +'*** ABCI CHOOSES THE FOLLOWING CONTOUR: ***'/ 1X,'ZCF=',G13.5,
     +'(M)    ZCT=',G13.5,'(M)    AND RWAK=', G13.5,'(M)'/ 1X,
     +'YOU SHOULD STILL CHECK IF THIS CONTOUR IS VALID OR NOT.')
10040 FORMAT(1H0,'*** CAUTION *** YOU MUST NOW USE NAPOLY INTEGRATION',
     +' FOR CALCULATION OF THE WAKE POTENTIALS.'/ 1X,
     +'YOUR CASE: LNAPOLY=',1X,L1/ 1X,
     +'IF LNAPOLY=F, ABCI NOW CHANGES IT TO LNAPOLY=T.'/ 1X,
     +'CHECK IF THE CONTOUR IS VALID.')
10050 FORMAT(1H ,'*** WARNING *** YOUR INTEGRATION CONTOUR MAY INTERS',
     +'ECT OR RUN TOO CLOSE TO THE MATRIAL.'/ 1X,
     +'CHECK THE CONTOUR: (IT MUST BE AT LEAST TWO MESHES AWAY FROM',
     +' THE VERTICAL BOUNDARY)'/ 1X,' ZCF=',G13.5,'(M)    AND ZLFT=',G13
     +.5,'(M)'/ 1X,' ZCT=',G13.5,'(M)    AND ZRGT=',G13.5,'(M)'/ 1X,
     +'RWAK=',G13.5,'(M)    AND RMIN=',G13.5,'(M)')
10060 FORMAT(1H ,'*** ABCI CHOOSES THE FOLLOWING CONTOUR INSTEAD: ***'/
     +1X,'ZCF=',G13.5,'(M)    ZCT=',G13.5,'(M)    AND RWAK=', G13.5,
     +'(M)'/ 1X,
     +'YOU SHOULD STILL CHECK IF THIS CONTOUR IS VALID OR NOT.')
10070 FORMAT(1H ,'*** WARNING *** NAPOLY INTEGRATION METHOD WAS NOT',
     +' CHOSEN DESPITE OF THE MATERIAL COMING DOWN '/ 1X,
     +'BELOW THE BEAM PIPE RADIUS',
     +' OR UNEQUAL BEAM PIPE RADII AT BOTH ENDS.'/ 1H ,
     +'*** ABCI CHOOSES THE FOLLOWING CONTOUR: ***'/ 1X,'ZCF=',G13.5,
     +'(M)    ZCT=',G13.5,'(M)    AND RWAK=', G13.5,'(M)'/ 1X,
     +'YOU SHOULD STILL CHECK IF THIS CONTOUR IS VALID OR NOT.')
10075 FORMAT(1H ,'*** WARNING *** THE CONVENTIONAL INTEGRATION METHOD'/
     +1X,'AT RWAK=',G13.5,'(M) WILL BE USED.')
10080 FORMAT(1H ,'*** SKIP THIS CALCULATION *** THE HORIZONATAL',
     +'INTEGRATION CONTOUR AT R=RWAK IS TOO CLOSE TO THE BEAM AXIS.'/ 1
     +X,'IT MUST BE AT LEAST 5 RADIAL MESHES AWAY FROM THE BEAM',
     +'AXIS FOR MROT=1.'/ 1X,'YOUR CASE: ONLY',I5,' MESHES AWAY'/ 1X,
     +'YOU SHOULD INCREASE RWAK OR USE A SMALLER RADIAL MESH SIZE')
      END
C
      SUBROUTINE SETKMIN
C-----------------------------------------------------------------------
C SET KMIN
C-----------------------------------------------------------------------
      LOGICAL LNAPOLY,LNONAP
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /PINT /KMIN
      COMMON /PIPE/ KOFF,KOFL
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
C
      IF(.NOT.LNONAP) THEN
        KMIN=KOFF
      ELSE
        KMIN=MIN(KOFF,KOFL)
        IF((IRWAK.GT.1).AND.(IRWAK.LT.KMIN)) KMIN=IRWAK
      END IF
      RETURN
      END 
C
      SUBROUTINE PLTMSH
C-----------------------------------------------------------------------
C PLOT MESH
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      CHARACTER*1 CHA(120*NRN),STAR,LSHP,SLSH,DOT,BEAM
      LOGICAL IHIT
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /BC/ IHIT(NHIT,8)
      COMMON /LFRBK/ LFR, LBK, MZF, MZL
C     DATA STAR,LSHP,SLSH,DOT,BEAM/'*','\','/','.','B'/
C FOR SUN AND HP USERS -- DELETE THE ABOVE LINE AND INCLUDE ONE BELOW
      DATA STAR,LSHP,SLSH,DOT,BEAM/'*','L','/','.','B'/
      WRITE(ND6,10000)
      LFR=1
      LBK=NHIT/NR-2
      MZF=1
      MZL=MZF+NMESH/NR-1
      CALL HIT1(0)
      NNZ=(NZ-1)/120
      IF(NNZ.EQ.0) THEN
        DO 30 I=1,NZ
          I1=I-1
          IF(I.EQ.1) I1=1
          CALL BCDET(I1,I,I1-1,I-1)
          IF(I.LE.MZL) GO TO 10
          MZF=I1
          MZL=MZF+NMESH/NR-1
   10     DO 20 J=1,NR
            K=(I-MZF)*NR+J
            KZ=(I-1)*NR+J
            IF(IFLG(K).EQ.5) CHA(KZ)=STAR
            IF(IFLG(K).EQ.4) CHA(KZ)=SLSH
            IF(IFLG(K).EQ.3) CHA(KZ)=SLSH
            IF(IFLG(K).EQ.2) CHA(KZ)=LSHP
            IF(IFLG(K).EQ.1) CHA(KZ)=LSHP
            IF(IFLG(K).EQ.0) CHA(KZ)=DOT
            IF((J.EQ.JOFF).AND.(I.EQ.1)) CHA(KZ)=BEAM
   20     CONTINUE
   30   CONTINUE
        DO 40 J=NR,1,-1
          WRITE(ND6,10020) (CHA((L-1)*NR+J),L=1,NZ)
   40   CONTINUE
      ELSE
        DO 90 N=1,NNZ
          N1=(N-1)*120
          N2=N*120
          DO 70 I=1+N1,N2
            I1=I-1
            IF(I.EQ.1) I1=1
            CALL BCDET(I1,I,I1-1,I-1)
            IF(I.LE.MZL) GO TO 50
            MZF=I1
            MZL=MZF+NMESH/NR-1
   50       DO 60 J=1,NR
              K=(I-MZF)*NR+J
              KZ=(I-N1-1)*NR+J
              IF(IFLG(K).EQ.5) CHA(KZ)=STAR
              IF(IFLG(K).EQ.4) CHA(KZ)=SLSH
              IF(IFLG(K).EQ.3) CHA(KZ)=SLSH
              IF(IFLG(K).EQ.2) CHA(KZ)=LSHP
              IF(IFLG(K).EQ.1) CHA(KZ)=LSHP
              IF(IFLG(K).EQ.0) CHA(KZ)=DOT
              IF((J.EQ.JOFF).AND.(I.EQ.1)) CHA(KZ)=BEAM
   60       CONTINUE
   70     CONTINUE
          DO 80 J=NR,1,-1
            WRITE(ND6,10020) (CHA((L-1)*NR+J),L=1,120)
   80     CONTINUE
          WRITE(ND6,10010)
   90   CONTINUE
        NNZ1=NNZ*120
        DO 120 I=1+NNZ1,NZ
          I1=I-1
          CALL BCDET(I1,I,I1-1,I-1)
          IF(I.LE.MZL) GO TO 100
          MZF=I1
          MZL=MZF+NMESH/NR-1
  100     DO 110 J=1,NR
            K=(I-MZF)*NR+J
            KZ=(I-NNZ1-1)*NR+J
            IF(IFLG(K).EQ.5) CHA(KZ)=STAR
            IF(IFLG(K).EQ.4) CHA(KZ)=SLSH
            IF(IFLG(K).EQ.3) CHA(KZ)=SLSH
            IF(IFLG(K).EQ.2) CHA(KZ)=LSHP
            IF(IFLG(K).EQ.1) CHA(KZ)=LSHP
            IF(IFLG(K).EQ.0) CHA(KZ)=DOT
            IF((J.EQ.JOFF).AND.(I.EQ.1)) CHA(KZ)=BEAM
  110     CONTINUE
  120   CONTINUE
        DO 130 J=NR,1,-1
          WRITE(ND6,10020) (CHA((L-1)*NR+J),L=1,NZ-NNZ*120)
  130   CONTINUE
      END IF
      RETURN
10000 FORMAT(1H1,'CAVITY SHAPE ' /1X,
C    +'.=VACUUM ; \,/=TRIANGULAR FILLED MESH CELL ; *=FULL',
C FOR SUN AND HP USERS -- DELETE THE ABOVE LINE AND INCLUDE ONE BELOW
     +'.=VACUUM ; L,/=TRIANGULAR FILLED MESH CELL ; *=FULL',
     +' METAL ; B=BEAM'/)
10010 FORMAT(1H1)
10020 FORMAT(1X,120A1)
      END
C
      SUBROUTINE EXPAND
C-----------------------------------------------------------------------
C EXPANDING THE MESH POINTS FOR SYMMETRY INPUT CASE
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      DIMENSION IRS(NBP),IZS(NBP)
      NZH=NZ/2
      DO 10 I=NVC+1,2,-1
        IZS(NVC+2-I)=NZH-IZZ(I)+2
        IF(I.EQ.NVC) IZS(NVC+2-I)=NZH-IZZ(I)+1
        IRS(NVC+2-I)=IRR(I)
        IF(R(I).GE.-2.) GO TO 10
        IRS(NVC+1-I)=IRR(I)
        IZS(NVC+1-I)=IZZ(I)
        IRS(NVC+2-I)=IRR(I+1)
        IZS(NVC+2-I)=NZH-IZZ(I+1)+2
   10 CONTINUE
      DO 20 I=2,NVC+1
        IZS(NVC+I-1)=IZZ(I)+NZH
        IRS(NVC+I-1)=IRR(I)
        IF(R(I).LE.-3.) IZS(NVC+I-1)=IZZ(I)
   20 CONTINUE
      IZS(2*NVC+1)=1
      IRS(2*NVC+1)=1
C
      IF(NZP.EQ.NVC+2) GO TO 50
      DO 30 I=NZP,NVC+3,-1
        IZS(NZP+2*NVC+2-I)=NZH-IZZ(I)+2
        IRS(NZP+2*NVC+2-I)=IRR(I)
        IF(R(I).GE.-2.) GO TO 30
        IRS(NZP+2*NVC+1-I)=IRR(I)
        IZS(NZP+2*NVC+1-I)=IZZ(I)
        IRS(NZP+2*NVC+2-I)=IRR(I+1)
        IZS(NZP+2*NVC+2-I)=NZH-IZZ(I+1)+2
   30 CONTINUE
      DO 40 I=NVC+3,NZP
        IZS(NZP+I-3)=IZZ(I)+NZH
        IRS(NZP+I-3)=IRR(I)
        IF(R(I).LE.-3.) IZS(NZP+I-3)=IZZ(I)
   40 CONTINUE
   50 DO 60 I=1,2*NZP-3
        IZZ(I)=IZS(I)
        IRR(I)=IRS(I)
   60 CONTINUE
      RETURN
      END
C
      SUBROUTINE EXPAND2
C-----------------------------------------------------------------------
C EXPANDING THE MESH POINTS FOR SYMMETRY INPUT CASE
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      DIMENSION RS(NBP),ZS(NBP)
      ZH=ZL/2.
      DO 10 I=NVC+1,2,-1
        ZS(NVC+2-I)=ZH-Z(I)
        RS(NVC+2-I)=R(I)
        IF(R(I).LE.-1.) ZS(NVC+2-I)=Z(I)
        IF(R(I).GE.-2.) GO TO 10
        RS(NVC+1-I)=R(I)
        ZS(NVC+1-I)=Z(I)
        RS(NVC+2-I)=R(I+1)
        ZS(NVC+2-I)=ZH-Z(I+1)
   10 CONTINUE
      DO 20 I=2,NVC+1
        ZS(NVC+I-1)=Z(I)+ZH
        RS(NVC+I-1)=R(I)
        IF(R(I).LE.-1.) ZS(NVC+I-1)=Z(I)
   20 CONTINUE
      ZS(2*NVC+1)=0.
      RS(2*NVC+1)=0.
C
      IF(NZP.EQ.NVC+2) GO TO 50
      DO 30 I=NZP,NVC+3,-1
        ZS(NZP+2*NVC+2-I)=ZH-Z(I)
        RS(NZP+2*NVC+2-I)=R(I)
        IF(R(I).LE.-1.) ZS(NZP+2*NVC+2-I)=Z(I)
        IF(R(I).GE.-2.) GO TO 30
        RS(NZP+2*NVC+1-I)=R(I)
        ZS(NZP+2*NVC+1-I)=Z(I)
        RS(NZP+2*NVC+2-I)=R(I+1)
        ZS(NZP+2*NVC+2-I)=ZH-Z(I+1)
   30 CONTINUE
      DO 40 I=NVC+3,NZP
        ZS(NZP+I-3)=Z(I)+ZH
        RS(NZP+I-3)=R(I)
        IF(R(I).LE.-1.) ZS(NZP+I-3)=Z(I)
   40 CONTINUE
   50 DO 60 I=1,2*NZP-3
        Z(I)=ZS(I)
        R(I)=RS(I)
   60 CONTINUE
      NVC=2*NVC-1
      NZP=2*NZP-3
      RETURN
      END
C
      SUBROUTINE CENTER(RCAN,ZCAN,RMAX,ZMAX,RMIN,RVAC,ZLFT,ZRGT,NNX)
C-----------------------------------------------------------------------
C FIND A CENTER OF A CIRCLE
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM
      DIMENSION AXS1(2),AXS2(2),AXS12(2),CENT(2),A1CENT(2)
      LOGICAL LLONG, LCAV,LVEX
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /SYMT/ SYM
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      ZLFT=ZCAN
      ZRGT=0.0
      RMIN=RVAC
      RMAX=RCAN
      ZMAX=ZCAN
      DZZ=DDZ
      IF(NR.NE.0) DRR=RCAN/FLOAT(NR-1)
      IF((NZ.NE.0).AND.(.NOT.SYM)) DZZ=ZCAN/FLOAT(NZ-1)
      IF((NZ.NE.0).AND.(SYM)) DZZ=ZCAN/FLOAT((NZ-1)/2)
      I=2
   10 IF(I.GT.NZP) GO TO 110
      IF((I.NE.NVC+1).AND.(I.NE.NVC+2)) GO TO 20
      I=I+1
      GO TO 10
   20 IF(R(I).GT.-1.) THEN
C-STRAIGHT LINE----
        IF(RMAX.LT.R(I)) RMAX=R(I)
        IF(ZMAX.LT.Z(I)) ZMAX=Z(I)
        IF(RMIN.GT.R(I)) RMIN=R(I)
        IF(SYM) THEN
          IF((RVAC.GT.R(I)).AND.(ZLFT.GT.Z(I)).AND.(I.NE.NVC)) ZLFT=Z(I)
          IF((RVAC.GT.R(I)).AND.(ZRGT.LT.Z(I)).AND.(I.NE.NVC)) ZRGT=Z(I)
        ELSE
          IF((RVAC.GT.R(I)).AND.(ZLFT.GT.Z(I)).AND.(I.NE.NVC).AND.
     +    (I.NE.2)) ZLFT=Z(I)
          IF((RVAC.GT.R(I)).AND.(ZRGT.LT.Z(I)).AND.(I.NE.NVC).AND.
     +    (I.NE.2)) ZRGT=Z(I)
        END IF
      IF((R(I+1).GT.-1.).AND.(I+1.NE.NVC)) THEN
        Z1=Z(I)
        Z2=Z(I+1)
        R1=R(I)
        R2=R(I+1)
        IF ((RVAC.LE.R1).AND.(RVAC.LE.R2)) GO TO 26
C       IF((.NOT.SYM).AND.(I.EQ.2)) GO TO 26
        DEL=DZZ/3.
        IZ1=NINT(Z1/DEL)+1
        IZ2=NINT(Z2/DEL)+1
        IR1=NINT(R1/DEL)+1
        IR2=NINT(R2/DEL)+1
        NSTPZ=1
        IF(IZ1.GT.IZ2) NSTPZ=-1
        NSTPR=1
        IF(IR1.GT.IR2) NSTPR=-1
C
        IF(IZ1.EQ.IZ2) GO TO 22
        DO 24 IZ=IZ1,IZ2,NSTPZ
          P=(R2-R1)/FLOAT(IZ2-IZ1)*FLOAT(IZ-IZ1)+R1
          Q=FLOAT(IZ-1)*DEL
          IF(SYM) THEN
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q).AND.(I.NE.NVC)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q).AND.(I.NE.NVC)) ZRGT=Q
          ELSE
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q).AND.(I.NE.NVC).AND.
     +      (I.NE.2)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q).AND.(I.NE.NVC).AND.
     +      (I.NE.2)) ZRGT=Q
          END IF
   24   CONTINUE
   22   IF(IR1.EQ.IR2) GO TO 26
        DO 28 IR=IR1,IR2,NSTPR
          Q=(Z2-Z1)/FLOAT(IR2-IR1)*FLOAT(IR-IR1)+Z1
          P=FLOAT(IR-1)*DEL
          IF(SYM) THEN
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q).AND.(I.NE.NVC)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q).AND.(I.NE.NVC)) ZRGT=Q
          ELSE
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q).AND.(I.NE.NVC).AND.
     +      (I.NE.2)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q).AND.(I.NE.NVC).AND.
     +      (I.NE.2)) ZRGT=Q
          END IF
   28   CONTINUE
   26   CONTINUE
        END IF
        I=I+1
      ELSE
        IKIND=R(I)
        IF(IKIND.GE.-2) THEN
C-CIRCLE----
          RAD=ABS(Z(I))
          LLONG=IKIND.EQ.-2
          LCAV=Z(I).GT.0.0
          LVEX=Z(I).LT.0.0
          AXS1(1)=R(I-1)
          AXS2(1)=R(I+1)
          AXS1(2)=Z(I-1)
          AXS2(2)=Z(I+1)
          AXS12(1)=AXS2(1)-AXS1(1)
          AXS12(2)=AXS2(2)-AXS1(2)
          DST12=SQRT(AXS12(1)**2+AXS12(2)**2)
          IF(DST12.LE.2.*RAD) GO TO 30
          DIF=ABS(DST12-2.*RAD)
          IF(DIF.GT.CEPS) WRITE(ND6,10000) AXS1(1),AXS1(2), RAD,0.50
     +    *DST12
          RAD=0.5*DST12
   30     ARG=DST12/(2.*RAD)
          PHIX=2.0*ASIN(ARG)
          PHX=ACOS(ARG)
          IF(LLONG) PHX=-PHX
          IF(LCAV) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)+SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LCAV) CENT(2)=AXS1(2)+(-SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2)
     +    ) *(RAD/DST12)
          IF(LVEX) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)-SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LVEX) CENT(2)=AXS1(2)+(SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2))
     +    *(RAD/DST12)
C
          IF(DDZ.NE.0.0) THEN
            IRR(I)=IRD(CENT(1))
            IZZ(I)=NINT(CENT(2)/DDZ)+1
            DRR=GDR(IRR(I))
          ELSE
            IRR(I)=NINT(CENT(1)/DRR)+1
            IZZ(I)=NINT(CENT(2)/DZZ)+1
          END IF
C
          IF((CENT(1)+RAD.LE.RCAN).AND.(CENT(2)+RAD.LE.ZCAN).AND. (CENT
     +    (1)-RAD.GE.RVAC)) GO TO 50
C
          A1CENT(1)=AXS1(1)-CENT(1)
          A1CENT(2)=AXS1(2)-CENT(2)
          IF(LLONG) PHIX=2.*PI-PHIX
          DDD=AMIN1(DZZ,DRR)
          NPOLY=(5.0*RAD/DDD)*(PHIX/PI)+1
          STEP=1.0/FLOAT(NPOLY)
          T=0.0
          DO 40 J=2,NPOLY
            T=T+STEP
            IF(LCAV) P=CENT(1)+COS(PHIX*T)*A1CENT(1)+SIN(PHIX*T)*A1CENT
     +      (2)
            IF(LCAV) Q=CENT(2)-SIN(PHIX*T)*A1CENT(1)+COS(PHIX*T)*A1CENT
     +      (2)
            IF(LVEX) P=CENT(1)+COS(PHIX*T)*A1CENT(1)-SIN(PHIX*T)*A1CENT
     +      (2)
            IF(LVEX) Q=CENT(2)+SIN(PHIX*T)*A1CENT(1)+COS(PHIX*T)*A1CENT
     +      (2)
            IF(RMAX.LT.P) RMAX=P
            IF(ZMAX.LT.Q) ZMAX=Q
            IF(RMIN.GT.P) RMIN=P
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q)) ZRGT=Q
   40     CONTINUE
   50     CONTINUE
          I=I+1
        ELSE
C-ELLIPSE
          R1=R(I-1)
          RM=R(I+1)
          R2=R(I+2)
          Z1=Z(I-1)
          ZM=Z(I+1)
          Z2=Z(I+2)
          IF(DDZ.NE.0.) DRR=GDR(IRD(RM))
          A=0.
          B=0.
          ELLP=(R1-RM)**2*(Z2-ZM)**2-(R2-RM)**2*(Z1-ZM)**2
          IF(SQRT(ABS(ELLP)).GE.CEPS) GO TO 60
          AB=Z(I)
          IF(AB.LE.0.) THEN
            WRITE(ND6,10010) RM,ZM,R1,Z1,R2,Z2,ELLP
            STOP
          ELSE
            IF(Z1.EQ.Z2) A=0.5*ABS(R1-R2)
            IF(Z1.EQ.Z2) B=AB
            IF(R1.EQ.R2) A=AB
            IF(R1.EQ.R2) B=0.5*ABS(Z1-Z2)
            IF((A.GT.0.).AND.(B.GT.0.).AND.(ABS(A).GE.DRR) .AND.(ABS(B)
     +      .GE.DZZ)) GO TO 80
            WRITE(ND6,10020) RM,ZM,A,B
            STOP
          END IF
   60     CONTINUE
          A=((Z2-ZM)**2-(Z1-ZM)**2)/ELLP
          B=((R1-RM)**2-(R2-RM)**2)/ELLP
          IF((A.GT.0.).AND.(B.GT.0.).AND.(ABS(A).GE.DRR) .AND.(ABS(B)
     +    .GE.DZZ)) GO TO 70
          WRITE(ND6,10020) RM,ZM,A,B
          STOP
   70     CONTINUE
          A=1.0/SQRT(A)
          B=1.0/SQRT(B)
   80     CONTINUE
C
          IRR(I)=NINT(A/DRR)+1
          IZZ(I)=NINT(B/DZZ)+1
          IF((RM+A.LE.RCAN).AND.(ZM+B.LE.ZCAN).AND.(RM-A.GE.RVAC))
     +    GO TO 100
C
          IF(R1.GT.RM) PHX=ACOS((Z1-ZM)/B)
          IF(R2.GT.RM) PHY=ACOS((Z2-ZM)/B)
          IF(R1.LT.RM) PHX=-ACOS((Z1-ZM)/B)
          IF(R2.LT.RM) PHY=-ACOS((Z2-ZM)/B)
          IF((R1.EQ.RM).AND.(Z1.GE.ZM)) PHX=0.
          IF((R1.EQ.RM).AND.(Z1.LT.ZM)) PHX=PI
          IF((R2.EQ.RM).AND.(Z2.GE.ZM)) PHY=0.
          IF((R2.EQ.RM).AND.(Z2.LT.ZM)) PHY=PI
          DPH1=PHY-PHX
          IF(DPH1.GE.2.0*PI) DPH1=DPH1-2.*PI
          IF(DPH1.LE.-2.0*PI) DPH1=DPH1+2.*PI
          IF(DPH1.GE.0.) DPH2=DPH1-2.*PI
          IF(DPH1.LT.0.) DPH2=DPH1+2.*PI
          IF((IKIND.EQ.-3).AND.(ABS(DPH1).LT.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-3).AND.(ABS(DPH2).LT.ABS(DPH1))) DELPH=DPH2
          IF((IKIND.EQ.-4).AND.(ABS(DPH1).GE.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-4).AND.(ABS(DPH2).GE.ABS(DPH1))) DELPH=DPH2
C
          NPOLY=5.*AMAX1(A/DRR,B/DZZ)+1
          STEP=1./FLOAT(NPOLY)
          T=0.0
          DO 90 J=2,NPOLY
            T=T+STEP
            P=RM+A*SIN(PHX+T*DELPH)
            Q=ZM+B*COS(PHX+T*DELPH)
            IF(RMAX.LT.P) RMAX=P
            IF(ZMAX.LT.Q) ZMAX=Q
            IF(RMIN.GT.P) RMIN=P
            IF((RVAC.GT.P).AND.(ZLFT.GT.Q)) ZLFT=Q
            IF((RVAC.GT.P).AND.(ZRGT.LT.Q)) ZRGT=Q
   90     CONTINUE
  100     CONTINUE
          I=I+2
        END IF
      END IF
      GO TO 10
  110 CONTINUE
C     IF((ZMAX.GT.Z(NZP)+CEPS).AND.(ZMAX.NE.ZCAN)) ZMAX=ZMAX+DZZ
      IF((ZMAX.GT.Z(NVC)+CEPS).AND.(NNX.EQ.1))     ZMAX=ZMAX+DZZ
      RETURN
10000 FORMAT(1X,'*** WARNING ***  AT (R, Z)=(',G13.5,G13.5,')',
     +' RADIUS: ',G13.5,' IS SHOTER THAN GAP/2: ',G13.5)
10010 FORMAT(1X,'*** ERROR ***', ' ELLIPSE GENERATION FAILED',
     +' AT (R,Z)=(',G13.5,G13.5,').'/ 1X,
     +'NO ELLIPSE IS UNIQUELY DEFINED IF R1=R2 OR Z1=Z2.'/ 1X,
     +'YOUR CASE: (R1,Z1)=(',G13.5,G13.5,')'/ 1X,'           (R2,Z2)=(',
     +G13.5,G13.5,')', 1X,' --> ELLP=',G13.5/ 1X,
     +'SPLIT IT INTO TWO IF THE ELLIPSE IS SYMMETRICAL')
10020 FORMAT(1X,'*** ERROR ***', ' ELLIPSE CONNECTION FAILED',
     +' AT (R,  Z)=(',G13.5,G13.5,').'/ 1X,'A=',G13.5,' B=',G13.5,
     +' ARE TOO SMALL.')
      END
C
      SUBROUTINE HIT1(KTYPE)
C-----------------------------------------------------------------------
C TRACE BOUNDARY POSITIONS ON MESHES
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NSP=20000)
      LOGICAL SYM,LRZF,LLONG,LCAV,LVEX
      DIMENSION IZS(NSP),IRS(NSP),ZR(NSP)
      DIMENSION AXS1(2),AXS2(2),AXS12(2),CENT(2),A1CENT(2)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /SYMT/ SYM
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      IF (KTYPE.EQ.0) CALL INITHT
      IR3=1
      IZ3=1
      N=2
      NTN=0
      LRZF=.FALSE.
      RDUM=0.0
      ZDUM=0.0
   10 IF(N.GT.NZP-1) GO TO 771
      IF(LRZF) THEN
        RDUM=R(N)
        ZDUM=Z(N)
        LRZF=.FALSE.
      ELSE
        IF(N.EQ.NVC) THEN
          N=N+2
          GO TO 10
        ELSE
          IF((ABS(RDUM-R(N)).GT.CEPS).OR.(ABS(ZDUM-Z(N)).GT.CEPS))
     +    GO TO 20
          IF(KTYPE.NE.0) WRITE(ND9,10000)
          LRZF=.TRUE.
          N=N+1
          GO TO 10
        END IF
      END IF
   20 ITOT=0
      DO 30 I=1,NSP
        IZS(I)=0
        IRS(I)=0
   30 ZR(I)=0.
      IF(R(N).LE.-1.) THEN
        IKIND=R(N)
        IF(IKIND.GE.-2) THEN
C-----------------------------------------------------------------------
C AUTO FITTING OF A CIRCLE
C-----------------------------------------------------------------------
          RAD=ABS(Z(N))
          IZF=IZZ(N)-NINT(RAD/DDZ)-1
          IF(IZF.LT.1) IZF=1
          IZL=IZZ(N)+NINT(RAD/DDZ)+1
          IF(IZL.GT.NZ) IZL=NZ
          IF((IZL.LT.LFR).OR.(IZF.GT.LBK)) N=N+1
          IF((IZL.LT.LFR).OR.(IZF.GT.LBK)) GO TO 10
          LLONG=IKIND.EQ.-2
          LCAV=Z(N).GT.0.0
          LVEX=Z(N).LT.0.0
          AXS1(1)=R(N-1)
          AXS2(1)=R(N+1)
          AXS1(2)=Z(N-1)
          AXS2(2)=Z(N+1)
          AXS12(1)=AXS2(1)-AXS1(1)
          AXS12(2)=AXS2(2)-AXS1(2)
          DST12=SQRT(AXS12(1)**2+AXS12(2)**2)
          IF(DST12.LE.2.*RAD) GO TO 40
          RAD=0.50*DST12
   40     ARG=DST12/(2.*RAD)
          PHIX=2.*ASIN(ARG)
          PHX=ACOS(ARG)
          IF(LLONG) PHX=-PHX
          IF(LCAV) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)+SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LCAV) CENT(2)=AXS1(2)+(-SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2)
     +    ) *(RAD/DST12)
          IF(LVEX) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)-SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LVEX) CENT(2)=AXS1(2)+(SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IRF=IRD(CENT(1)-RAD)-1
          IF(IRF.LT.1) IRF=1
          IRL=IRD(CENT(1)+RAD)+1
          IF(IRL.GT.NR) IRL=NR
          A1CENT(1)=AXS1(1)-CENT(1)
          A1CENT(2)=AXS1(2)-CENT(2)
          IF((IZZ(N-1).LT.LFR).OR.(IZZ(N-1).GT.LBK)) GO TO 50
          ITOT=ITOT+1
          IF(ITOT.GT.NSP) GO TO 260
          IRS(1)=IRR(N-1)
          IZS(1)=IZZ(N-1)-LFR
          ZR(1)=0.
   50     IF(LLONG) PHIX=2.*PI-PHIX
          DO 70 IZ=IZF,IZL
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 70
            ZZ=FLOAT(IZ-1)*DDZ
            SR=RAD**2-(ZZ-CENT(2))**2
            IF(SR.LT.0.) GO TO 70
            RP=CENT(1)+ABS(SQRT(SR))
            RM=CENT(1)-ABS(SQRT(SR))
            ATPS=(((RP-CENT(1))*A1CENT(2)-A1CENT(1)*(ZZ-CENT(2))) /RAD
     +      **2)
            ATMS=(((RM-CENT(1))*A1CENT(2)-A1CENT(1)*(ZZ-CENT(2))) /RAD
     +      **2)
            IF(LVEX) ATPS=-ATPS
            IF(LVEX) ATMS=-ATMS
            IF(ABS(ATPS).GT.1.0) ATPS=ATPS/ABS(ATPS)
            IF(ABS(ATMS).GT.1.0) ATMS=ATMS/ABS(ATMS)
            ATPC=(((RP-CENT(1))*A1CENT(1)+A1CENT(2)*(ZZ-CENT(2))) /RAD
     +      **2)
            ATMC=(((RM-CENT(1))*A1CENT(1)+A1CENT(2)*(ZZ-CENT(2))) /RAD
     +      **2)
            IF(ABS(ATPC).GT.1.0) ATPC=ATPC/ABS(ATPC)
            IF(ABS(ATMC).GT.1.0) ATMC=ATMC/ABS(ATMC)
            IF((ATPS.GE.0.).AND.(ATPC.GE.0.)) TPS=ASIN(ATPS)
            IF((ATPS.GE.0.).AND.(ATPC.LT.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.LE.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.GT.0.)) TPS=2.*PI+ASIN(ATPS)
            IF((ATMS.GE.0.).AND.(ATMC.GE.0.)) TMS=ASIN(ATMS)
            IF((ATMS.GE.0.).AND.(ATMC.LT.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.LE.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.GT.0.)) TMS=2.*PI+ASIN(ATMS)
            IF((TPS.LT.0.).OR.(TPS.GT.PHIX)) GO TO 60
            IR=IRD(RP)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=TPS
   60       IF((TMS.LT.0.).OR.(TMS.GT.PHIX)) GO TO 70
            IR=IRD(RM)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=TMS
   70     CONTINUE
          DO 90 IR=IRF,IRL
            RR=GR(IR)
            SZ=RAD**2-(RR-CENT(1))**2
            IF(SZ.LT.0.) GO TO 90
            ZP=CENT(2)+ABS(SQRT(SZ))
            ZM=CENT(2)-ABS(SQRT(SZ))
            ATPS=(((RR-CENT(1))*A1CENT(2)-A1CENT(1)*(ZP-CENT(2))) /RAD
     +      **2)
            ATMS=(((RR-CENT(1))*A1CENT(2)-A1CENT(1)*(ZM-CENT(2))) /RAD
     +      **2)
            IF(LVEX) ATPS=-ATPS
            IF(LVEX) ATMS=-ATMS
            IF(ABS(ATPS).GT.1.0) ATPS=ATPS/ABS(ATPS)
            IF(ABS(ATMS).GT.1.0) ATMS=ATMS/ABS(ATMS)
            ATPC=(((RR-CENT(1))*A1CENT(1)+A1CENT(2)*(ZP-CENT(2))) /RAD
     +      **2)
            ATMC=(((RR-CENT(1))*A1CENT(1)+A1CENT(2)*(ZM-CENT(2))) /RAD
     +      **2)
            IF(ABS(ATPC).GT.1.0) ATPC=ATPC/ABS(ATPC)
            IF(ABS(ATMC).GT.1.0) ATMC=ATMC/ABS(ATMC)
            IF((ATPS.GE.0.).AND.(ATPC.GE.0.)) TPS=ASIN(ATPS)
            IF((ATPS.GE.0.).AND.(ATPC.LT.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.LE.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.GT.0.)) TPS=2.*PI+ASIN(ATPS)
            IF((ATMS.GE.0.).AND.(ATMC.GE.0.)) TMS=ASIN(ATMS)
            IF((ATMS.GE.0.).AND.(ATMC.LT.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.LE.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.GT.0.)) TMS=2.*PI+ASIN(ATMS)
            IF((TPS.LT.0.).OR.(TPS.GT.PHIX)) GO TO 80
            IZ=NINT(ZP/DDZ-0.00000)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 80
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=TPS
   80       IF((TMS.LT.0.).OR.(TMS.GT.PHIX)) GO TO 90
            IZ=NINT(ZM/DDZ-0.00000)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 90
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=TMS
   90     CONTINUE
          IF((IZZ(N+1).LT.LFR).OR.(IZZ(N+1).GT.LBK)) GO TO 100
          ITOT=ITOT+1
          IF(ITOT.GT.NSP) GO TO 260
          IZS(ITOT)=IZZ(N+1)-LFR
          IRS(ITOT)=IRR(N+1)
          ZR(ITOT)=PHIX
  100     CONTINUE
          N=N+1
        ELSE
C-----------------------------------------------------------------------
C AUTO FITTING OF AN ELLIPSE
C-----------------------------------------------------------------------
          IA=IRR(N)
          IB=IZZ(N)
          IZF=IZZ(N+1)-IB-1
          IF(IZF.LT.1) IZF=1
          IZL=IZZ(N+1)+IB+1
          IF(IZL.GT.NZ) IZL=NZ
          IF((IZL.LT.LFR).OR.(IZF.GT.LBK)) N=N+2
          IF((IZL.LT.LFR).OR.(IZF.GT.LBK)) GO TO 10
          R1=R(N-1)
          RM=R(N+1)
          CENT(1)=RM
          R2=R(N+2)
          Z1=Z(N-1)
          ZM=Z(N+1)
          CENT(2)=ZM
          Z2=Z(N+2)
          ELLP=(R1-RM)**2*(Z2-ZM)**2-(R2-RM)**2*(Z1-ZM)**2
          IF(SQRT(ABS(ELLP)).GE.CEPS) GO TO 110
          AB=Z(N)
          IF(Z1.EQ.Z2) A=0.50*ABS(R1-R2)
          IF(Z1.EQ.Z2) B=AB
          IF(R1.EQ.R2) A=AB
          IF(R1.EQ.R2) B=0.50*ABS(Z1-Z2)
          GO TO 120
  110     A=((Z2-ZM)**2-(Z1-ZM)**2)/ELLP
          B=((R1-RM)**2-(R2-RM)**2)/ELLP
          A=1.0/SQRT(A)
          B=1.0/SQRT(B)
  120     CONTINUE
          IRF=IRD(CENT(1)-A)-1
          IF(IRF.LT.1) IRF=1
          IRL=IRD(CENT(1)+A)+1
          IF(IRL.GT.NR) IRL=NR
          ITOT=ITOT+1
          IF(ITOT.GT.NSP) GO TO 260
          IRS(1)=IRR(N-1)
          IZS(1)=IZZ(N-1)-LFR
          ZR(1)=-9999.
          IF(R1.GT.RM) PHX=ACOS((Z1-ZM)/B)
          IF(R2.GT.RM) PHY=ACOS((Z2-ZM)/B)
          IF(R1.LT.RM) PHX=-ACOS((Z1-ZM)/B)
          IF(R2.LT.RM) PHY=-ACOS((Z2-ZM)/B)
          IF((R1.EQ.RM).AND.(Z1.GE.ZM)) PHX=0.
          IF((R1.EQ.RM).AND.(Z1.LT.ZM)) PHX=PI
          IF((R2.EQ.RM).AND.(Z2.GE.ZM)) PHY=0.
          IF((R2.EQ.RM).AND.(Z2.LT.ZM)) PHY=PI
          DPH1=PHY-PHX
          IF(DPH1.GE.2.0*PI) DPH1=DPH1-2.*PI
          IF(DPH1.LE.-2.0*PI) DPH1=DPH1+2.*PI
          IF(DPH1.GE.0.) DPH2=DPH1-2.*PI
          IF(DPH1.LT.0.) DPH2=DPH1+2.*PI
          IF((IKIND.EQ.-3).AND.(ABS(DPH1).LT.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-3).AND.(ABS(DPH2).LT.ABS(DPH1))) DELPH=DPH2
          IF((IKIND.EQ.-4).AND.(ABS(DPH1).GE.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-4).AND.(ABS(DPH2).GE.ABS(DPH1))) DELPH=DPH2
          IF(DELPH.GE.0.) THEN
            AGF=PHX
            AGL=PHX+DELPH
          ELSE
            AGF=PHX+DELPH
            AGL=PHX
          END IF
          SGNR=DELPH/ABS(DELPH)
          DO 160 IZ=IZF,IZL
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 160
            ZZ=FLOAT(IZ-1)*DDZ
            SR=1.0-((ZZ-CENT(2))/B)**2
            IF(SR.LT.0.) GO TO 160
            RP=CENT(1)+A*ABS(SQRT(SR))
            RM=CENT(1)-A*ABS(SQRT(SR))
            ATPS=(RP-CENT(1))/A
            ATMS=(RM-CENT(1))/A
            IF(ABS(ATPS).GT.1.0) ATPS=ATPS/ABS(ATPS)
            IF(ABS(ATMS).GT.1.0) ATMS=ATMS/ABS(ATMS)
            ATPC=(ZZ-CENT(2))/B
            ATMC=(ZZ-CENT(2))/B
            IF(ABS(ATPC).GT.1.0) ATPC=ATPC/ABS(ATPC)
            IF(ABS(ATMC).GT.1.0) ATMC=ATMC/ABS(ATMC)
            IF((ATPS.GE.0.).AND.(ATPC.GE.0.)) TPS=ASIN(ATPS)
            IF((ATPS.GE.0.).AND.(ATPC.LT.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.LE.0.)) TPS=PI-ASIN(ATPS)
            IF((ATPS.LT.0.).AND.(ATPC.GT.0.)) TPS=2.*PI+ASIN(ATPS)
            IF((ATMS.GE.0.).AND.(ATMC.GE.0.)) TMS=ASIN(ATMS)
            IF((ATMS.GE.0.).AND.(ATMC.LT.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.LE.0.)) TMS=PI-ASIN(ATMS)
            IF((ATMS.LT.0.).AND.(ATMC.GT.0.)) TMS=2.*PI+ASIN(ATMS)
            TP2=TPS-2.*PI
            TM2=TMS-2.*PI
            IF((TPS.LT.AGF).OR.(TPS.GT.AGL)) GO TO 130
            IR=IRD(RP)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TPS-AGF)*SGNR
  130       IF((TMS.LT.AGF).OR.(TMS.GT.AGL)) GO TO 140
            IR=IRD(RM)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TMS-AGF)*SGNR
  140       CONTINUE
            IF((TP2.LT.AGF).OR.(TP2.GT.AGL)) GO TO 150
            IR=IRD(RP)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TP2-AGF)*SGNR
  150       IF((TM2.LT.AGF).OR.(TM2.GT.AGL)) GO TO 160
            IR=IRD(RM)
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TM2-AGF)*SGNR
  160     CONTINUE
          DO 200 IR=IRF,IRL
            RR=GR(IR)
            SZ=1.0-((RR-CENT(1))/A)**2
            IF(SZ.LT.0.) GO TO 200
            ZP=CENT(2)+B*ABS(SQRT(SZ))
            ZM=CENT(2)-B*ABS(SQRT(SZ))
            ATPS=(RR-CENT(1))/A
            ATMS=(RR-CENT(1))/A
            IF(ABS(ATPS).GT.1.0) ATPS=ATPS/ABS(ATPS)
            IF(ABS(ATMS).GT.1.0) ATMS=ATMS/ABS(ATMS)
            ATPC=(ZP-CENT(2))/B
            ATMC=(ZM-CENT(2))/B
            IF(ABS(ATPC).GT.1.0) ATPC=ATPC/ABS(ATPC)
            IF(ABS(ATMC).GT.1.0) ATMC=ATMC/ABS(ATMC)
            IF((ATPS.GE.0.).AND.(ATPC.GE.0.)) TPC=ACOS(ATPC)
            IF((ATPS.GE.0.).AND.(ATPC.LT.0.)) TPC=ACOS(ATPC)
            IF((ATPS.LT.0.).AND.(ATPC.LE.0.)) TPC=2.*PI-ACOS(ATPC)
            IF((ATPS.LT.0.).AND.(ATPC.GT.0.)) TPC=2.*PI-ACOS(ATPC)
            IF((ATMS.GE.0.).AND.(ATMC.GE.0.)) TMC=ACOS(ATMC)
            IF((ATMS.GE.0.).AND.(ATMC.LT.0.)) TMC=ACOS(ATMC)
            IF((ATMS.LT.0.).AND.(ATMC.LE.0.)) TMC=2.*PI-ACOS(ATMC)
            IF((ATMS.LT.0.).AND.(ATMC.GT.0.)) TMC=2.*PI-ACOS(ATMC)
            TP2=TPC-2.*PI
            TM2=TMC-2.*PI
            IF((TPC.LT.AGF).OR.(TPC.GT.AGL)) GO TO 170
            IZ=NINT(ZP/DDZ-0.00001)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 170
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TPC-AGF)*SGNR
  170       IF((TMC.LT.AGF).OR.(TMC.GT.AGL)) GO TO 180
            IZ=NINT(ZM/DDZ-0.00001)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 180
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TMC-AGF)*SGNR
  180       CONTINUE
            IF((TP2.LT.AGF).OR.(TP2.GT.AGL)) GO TO 190
            IZ=NINT(ZP/DDZ-0.00001)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 190
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TP2-AGF)*SGNR
  190       IF((TM2.LT.AGF).OR.(TM2.GT.AGL)) GO TO 200
            IZ=NINT(ZM/DDZ-0.00001)+1
            IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 200
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR
            IRS(ITOT)=IR
            ZR(ITOT)=(TM2-AGF)*SGNR
  200     CONTINUE
          IF((IZZ(N+2).LT.LFR).OR.(IZZ(N+2).GT.LBK)) GO TO 210
          ITOT=ITOT+1
          IZS(ITOT)=IZZ(N+2)-LFR
          IRS(ITOT)=IRR(N+2)
          ZR(ITOT)=9998.
  210     CONTINUE
          N=N+2
        END IF
      ELSE
C-----------------------------------------------------------------------
C AUTO FITTING OF A STRAIGHT LINE
C-----------------------------------------------------------------------
        IF(R(N+1).LE.-1.) N=N+1
        IF(R(N).LE.-1.) GO TO 10
        IZF=IZZ(N)
        IZL=IZZ(N+1)
        IRF=IRR(N)
        IRL=IRR(N+1)
        ZF=Z(N)
        ZE=Z(N+1)
        RF=R(N)
        RE=R(N+1)
        NSTPZ=1
        IF(IZF.GT.IZL) NSTPZ=-1
        NSTPR=1
        IF(IRF.GT.IRL) NSTPR=-1
        IF(NSTPZ.GE.0) THEN
          IF((IZL.LT.LFR).OR.(IZF.GT.LBK)) GO TO 250
        ELSE
          IF((IZF.LT.LFR).OR.(IZL.GT.LBK)) GO TO 250
        END IF
        ITOT=ITOT+1
        IZS(ITOT)=IZF-LFR
        IRS(ITOT)=IRF
        ZR(ITOT)=0.
        IF(IZL.EQ.IZF) GO TO 230
        DO 220 IZ=IZF,IZL,NSTPZ
          RR=(GR(IRL)-GR(IRF))/FLOAT(IZL-IZF)*FLOAT(IZ-IZF)+GR(IRF)
          IR=IRD(RR)
          IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 220
          ITOT=ITOT+1
          IF(ITOT.GT.NSP) GO TO 260
          IZS(ITOT)=IZ-LFR
          IRS(ITOT)=IR
          ZCON=(FLOAT(IZS(ITOT)+LFR-IZF)/FLOAT(IZL-IZF))**2
          RCON=0.
          IF(IRF.NE.IRL) RCON=((GR(IR)-RF)/(RE-RF))**2
          ZR(ITOT)=SQRT(ZCON+RCON)
          IF(N.NE.NVC-1) GO TO 220
          IF(IZ.EQ.NZ) THEN
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR+1
            IRS(ITOT)=IR
            ZCON=(FLOAT(IZS(ITOT)+LFR-IZF)/FLOAT(IZL-IZF))**2
            RCON=0.
            IF(IRF.NE.IRL) RCON=((GR(IR)-RF)/(RE-RF))**2
            ZR(ITOT)=SQRT(ZCON+RCON)
          END IF
  220   CONTINUE
  230   IF(IRL.EQ.IRF) GO TO 250
        DO 240 IR=IRF,IRL,NSTPR
          IZ=NINT(FLOAT(IZL-IZF)*(GR(IR)-GR(IRF))/(GR(IRL)-GR(IRF))
     +    -0.00001)+IZF
          IF((IZ.LT.LFR).OR.(IZ.GT.LBK)) GO TO 240
          ITOT=ITOT+1
          IF(ITOT.GT.NSP) GO TO 260
          IZS(ITOT)=IZ-LFR
          IRS(ITOT)=IR
          ZCON=0.
          IF(IZL.NE.IZF) ZCON=(FLOAT(IZS(ITOT)+LFR-IZF)/FLOAT(IZL-IZF))
     +    **2
          RCON=((GR(IR)-RF)/(RE-RF))**2
          ZR(ITOT)=SQRT(ZCON+RCON)
          IF(N.NE.NVC-1) GO TO 240
          IF(IZ.EQ.NZ) THEN
            ITOT=ITOT+1
            IF(ITOT.GT.NSP) GO TO 260
            IZS(ITOT)=IZ-LFR+1
            IRS(ITOT)=IR
            ZCON=0.
            IF(IZL.NE.IZF) ZCON=(FLOAT(IZS(ITOT)+LFR-IZF)/FLOAT(IZL-IZF)
     +      )**2
            RCON=((GR(IR)-RF)/(RE-RF))**2
            ZR(ITOT)=SQRT(ZCON+RCON)
          END IF
  240   CONTINUE
  250   CONTINUE
        N=N+1
      END IF
C-----------------------------------------------------------------------
C NOW SORT IZS AND IRS AND TRACE MOTION OF BOUNDARY
C-----------------------------------------------------------------------
      IF((IZL.EQ.IZF).AND.(IRL.EQ.IRF)) GO TO 10
      CALL TRACE(KTYPE,ITOT,IRS,IZS,ZR,NTN,IR3,IZ3)
      GO TO 10
  771 RETURN
  260 WRITE(ND6,10010)
      STOP
10000 FORMAT(1X,' JOIN 1 SOLID'/ 1X,' PLOT')
10010 FORMAT(1X,'*** STOP ***  THE NUMBER OF MESH FITTING POINTS', I5,
     +' EXCEEDS MAXIMUM NSP:',I5)
      END
C
      SUBROUTINE TRACE(KTYPE,ITOT,IRS,IZS,ZR,NTN,IR3,IZ3)
C-----------------------------------------------------------------------
C SORT IZS AND IRS AND TRACE MOTION OF BOUNDARY
C-----------------------------------------------------------------------
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NSP=20000)
      LOGICAL IHIT
      DIMENSION IZS(NSP),IRS(NSP),ZR(NSP)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /BC/ IHIT(NHIT,8)
C-----------------------------------------------------------------------
C NOW SORT IN ORDER
C-----------------------------------------------------------------------
      DO 20 I=1,ITOT-1
        ALOW=ZR(I)
        DO 10 J=I+1,ITOT
          IF(ZR(J).GE.ALOW) GO TO 10
          IZT=IZS(I)
          IRT=IRS(I)
          IZS(I)=IZS(J)
          IRS(I)=IRS(J)
          ZR(I)=ZR(J)
          IZS(J)=IZT
          IRS(J)=IRT
          ZR(J)=ALOW
          ALOW=ZR(I)
   10   CONTINUE
   20 CONTINUE
      IF(KTYPE.EQ.0) GO TO 40
      DO 30 I=1,ITOT
        IZ1=IZS(I)
        IR1=IRS(I)
        ZZ1=FLOAT(IZ1)*DDZ
        RR1=GR(IR1)
        IF((IZ1.NE.IZ3).OR.(IR1.NE.IR3)) THEN
          WRITE(ND9,10000) ZZ1,RR1
          NTN=NTN+1
          IF(NTN.GE.1021) THEN
            WRITE(ND9,10010)
            WRITE(ND9,10000) ZZ1,RR1
            NTN=1
          ELSE
          END IF
        ELSE
        END IF
        IZ3=IZ1
        IR3=IR1
   30 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C NOW DETERMINE HOW THE BOUNDARY MOVES ON MESHES
C-----------------------------------------------------------------------
   40 CONTINUE
      DO 50 I=1,ITOT-1
        IZ1=IZS(I)
        IZ2=IZS(I+1)
        IR1=IRS(I)
        IR2=IRS(I+1)
        IP=IZ1*NR+IR1
        IF((IR2.EQ.IR1).AND.(IZ2.EQ.IZ1)) GO TO 50
        IF((IR2.EQ.IR1+1).AND.(IZ2.EQ.IZ1))  IHIT(IP,4)=.TRUE.
        IF((IR2.EQ.IR1+1).AND.(IZ2.EQ.IZ1+1)) IHIT(IP,6)=.TRUE.
        IF((IR2.EQ.IR1).AND.(IZ2.EQ.IZ1+1))  IHIT(IP,2)=.TRUE.
        IF((IR2.EQ.IR1-1).AND.(IZ2.EQ.IZ1+1)) IHIT(IP-1,7)=.TRUE.
        IF((IR2.EQ.IR1-1).AND.(IZ2.EQ.IZ1))  IHIT(IP-1,3)=.TRUE.
        IF((IR2.EQ.IR1-1).AND.(IZ2.EQ.IZ1-1)) IHIT(IP-NR-1,5)=.TRUE.
        IF((IR2.EQ.IR1).AND.(IZ2.EQ.IZ1-1))  IHIT(IP-NR,1)=.TRUE.
        IF((IR2.EQ.IR1+1).AND.(IZ2.EQ.IZ1-1)) IHIT(IP-NR,8)=.TRUE.
   50 CONTINUE
      RETURN
10000 FORMAT(1X,2G13.5)
10010 FORMAT(1X,' JOIN 1 SOLID'/ 1X,' PLOT')
      END
C
      SUBROUTINE ADVAN(LREC,IA,IB,NFST)
C-----------------------------------------------------------------------
C RESTORE EM-FILEDS FROM A FILE IF LREC=.T. AND ADVANCE MESH
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN,LREC
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WPL/ MROT
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /ELINE/ NZF,NZL
      COMMON /NZIN/ NZF1,NZL1
      COMMON /PIPE/ KOFF,KOFL
      IF(.NOT.LCHIN.OR.(MROT.EQ.0)) IC=1
      IF(LCHIN.AND.(MROT.EQ.1)) IC=0
      ID=1
C-----------------------------------------------------------------------
C RESORE EM-FILEDS FROM A FILE IF LREC=.T.
C-----------------------------------------------------------------------
      IF(LREC) THEN
        NMN=(NMESH-1)/8000
        NWNN=(NWN-1)/8000
        NRNN=(NRN-1)/400
        READ(ND2) IC,ID
        DO 10 I=1,NWNN+1
          ILAST=I*8000
          IF(ILAST.GT.NWN) ILAST=NWN
   10   READ(ND2) (WZ(J),J=(I-1)*8000+1,ILAST)
        DO 20 I=1,NWNN+1
          ILAST=I*8000
          IF(ILAST.GT.NWN) ILAST=NWN
   20   READ(ND2) (WY(J),J=(I-1)*8000+1,ILAST)
        DO 30 K=1,NPT
          DO 30 L=1,2
            DO 30 I=1,NRNN+1
              ILAST=I*400
              IF(ILAST.GT.NRN) ILAST=NRN
   30   READ(ND2) (HVSTK(J,L,K),J=(I-1)*400+1,ILAST)
        DO 40 K=1,NPT
          DO 40 L=1,2
            DO 40 I=1,NRNN+1
              ILAST=I*400
              IF(ILAST.GT.NRN) ILAST=NRN
   40   READ(ND2) (HVSTL(J,L,K),J=(I-1)*400+1,ILAST)
        DO 50 K=1,3
          DO 50 I=1,NMN+1
            ILAST=I*8000
            IF(ILAST.GT.NMESH) ILAST=NMESH
   50   READ(ND2) (HV(K,J),J=(I-1)*8000+1,ILAST)
        DO 60 K=1,3
          DO 60 I=1,NMN+1
            ILAST=I*8000
            IF(ILAST.GT.NMESH) ILAST=NMESH
   60   READ(ND2) (EV(K,J),J=(I-1)*8000+1,ILAST)
      ELSE
      END IF
      IA=IC
      IB=ID
C IF .NOT.LREC, CALCULATION STARTS FROM THE POINT WHERE THE GEOMETRY IS
C CHANGED TO SAVE CPU TIME
      IF(JOFF.NE.KOFF) GO TO 70
      IF(IRR(2).NE.IRR(1)) GO TO 70
      IF(.NOT.LREC) THEN
        NFST=(IZZ(3)-IZZ(2))*M
        IA=IA+IZZ(3)-IZZ(2)
        IB=IB+IZZ(3)-IZZ(2)
      ELSE
      END IF
   70 CONTINUE
C-----------------------------------------------------------------------
C ADVANCE MESHES TO STARTING POINT
C-----------------------------------------------------------------------
      LFR=1
      LBK=NHIT/NR-2
      MZF=1
      MZL=MZF+NMESH/NR-1
      CALL HIT1(0)
      NZF1=1
      NZL1=1
      CALL BCDET(1,1,NZF1,NZL1)
      CALL BCDET(1,2,NZF1,NZL1)
      NZF=1
      NZL=2
      IF((KOFF.NE.JOFF).AND.(.NOT.LREC)) GO TO 100
      DO 90 N=0,NFST
        NZL=N/M+2
        IF(NZL.GE.NZ) NZL=NZ
        IF(LCBACK) THEN
          NZF=N/M-NT/M+NZ-2
          IF(NZF.LE.1) NZF=1
        ELSE
          NZF=1
        END IF
        IF((NZF.NE.NZF1).OR.(NZL.NE.NZL1)) CALL BCDET(NZF,NZL,NZF1,NZL1)
        IF (NZF.EQ.NZF1) GO TO 80
        IF((NZL.EQ.MZL+1).AND.(NZL.NE.NZL1)) THEN
          MZF=NZF
          MZL=MZF+NMESH/NR-1
        ELSE
        END IF
   80   CONTINUE
        NZF1=NZF
        NZL1=NZL
   90 CONTINUE
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE STORE
C-----------------------------------------------------------------------
C STORE EM-FILEDS IN A FILE IF LSAV=.T.
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      NMN=(NMESH-1)/8000
      NWNN=(NWN-1)/8000
      NRNN=(NRN-1)/400
      DO 10 I=1,NWNN+1
        ILAST=I*8000
        IF(ILAST.GT.NWN) ILAST=NWN
   10 WRITE(ND1) (WZ(J),J=(I-1)*8000+1,ILAST)
      DO 20 I=1,NWNN+1
        ILAST=I*8000
        IF(ILAST.GT.NWN) ILAST=NWN
   20 WRITE(ND1) (WY(J),J=(I-1)*8000+1,ILAST)
      DO 30 K=1,NPT
        DO 30 L=1,2
          DO 30 I=1,NRNN+1
            ILAST=I*400
            IF(ILAST.GT.NRN) ILAST=NRN
   30 WRITE(ND1) (HVSTK(J,L,K),J=(I-1)*400+1,ILAST)
      DO 40 K=1,NPT
        DO 40 L=1,2
          DO 40 I=1,NRNN+1
            ILAST=I*400
            IF(ILAST.GT.NRN) ILAST=NRN
   40 WRITE(ND1) (HVSTL(J,L,K),J=(I-1)*400+1,ILAST)
      DO 50 K=1,3
        DO 50 I=1,NMN+1
          ILAST=I*8000
          IF(ILAST.GT.NMESH) ILAST=NMESH
   50 WRITE(ND1) (HV(K,J),J=(I-1)*8000+1,ILAST)
      DO 60 K=1,3
        DO 60 I=1,NMN+1
          ILAST=I*8000
          IF(ILAST.GT.NMESH) ILAST=NMESH
   60 WRITE(ND1) (EV(K,J),J=(I-1)*8000+1,ILAST)
      RETURN
      END
C
      SUBROUTINE BCDET(NZF,NZL,NZF1,NZL1)
C-----------------------------------------------------------------------
C DETERMINATION OF BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NHIT=3000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL IHIT,OP,SYM
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /BC/ IHIT(NHIT,8)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /SYMT/ SYM
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
C
      LZF=MZF
      LZL=MZL
C
      IF((NZL.NE.NZL1).AND.(NZL.EQ.LBK)) THEN
        LFR=LBK-2
        LBK=LFR+NHIT/NR-2
        CALL HIT1(0)
      ELSE
      END IF
      IO=-(LFR-1)*NR
C
      IF(NZF.EQ.NZF1) GO TO 30
      IF((NZL.EQ.MZL+1).AND.(NZL.NE.NZL1)) THEN
        DO 10 L=NZF1-MZF+2,NZL1-MZF+1
          L1=(L-1)*NR
          LL=L-(NZF1-MZF+2)+1
          LL1=(LL-1)*NR
          DO 10 J=1,NR
            K=L1+J
            KR=LL1+J
            IFLG(KR)=IFLG(K)
            SJDG(1,KR)=SJDG(1,K)
            SJDG(2,KR)=SJDG(2,K)
            SJDG(3,KR)=SJDG(3,K)
   10   CONTINUE
        K1=(NZL1-NZF1)*NR+1
        DO 20 K=K1,NMESH
          IFLG(K)=5
          SJDG(1,K)=0.
          SJDG(2,K)=0.
          SJDG(3,K)=0.
   20   CONTINUE
        LZF=NZF
        LZL=LZF+NMESH/NR-1
      ELSE
      END IF
   30 CONTINUE
C
      IF(NZL1.EQ.NZ) GO TO 170
   40 L=NZL
      OP=.TRUE.
      J=1
      KEY=0
   50 CONTINUE
      IF(J.GT.NR) GO TO 150
C
      K=(NZL-LZF)*NR+J
      KL=(L-1)*NR+J+IO
C
      IFLG(K)=KEY
C
      IF(.NOT.IHIT(KL,5)) GO TO 60
      IF(.NOT.IHIT(KL,6)) THEN
        IFLG(K)=3
        KEY=0
      ELSE
        IFK=IFLG(K-1)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.3)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
   60 CONTINUE
      IF(.NOT.IHIT(KL,6)) GO TO 70
      IF(.NOT.IHIT(KL,5)) THEN
        IFLG(K)=4
        KEY=5
      ELSE
        IFK=IFLG(K-1)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.3)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
   70 CONTINUE
      IF(.NOT.IHIT(KL,7)) GO TO 80
      IF(.NOT.IHIT(KL,8)) THEN
        IFLG(K)=2
        KEY=5
      ELSE
        IFK=IFLG(K-1)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.3)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
   80 CONTINUE
      IF(.NOT.IHIT(KL,8)) GO TO 90
      IF(.NOT.IHIT(KL,7)) THEN
        IFLG(K)=1
        KEY=0
      ELSE
        IFK=IFLG(K-1)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.3)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
   90 CONTINUE
      IF(.NOT.IHIT(KL,1)) GO TO 100
      IF(.NOT.IHIT(KL,2)) THEN
        IFLG(K)=0
        KEY=0
      ELSE
        IFLG(K)=KEY
C     KEY=0
      END IF
      GO TO 140
  100 CONTINUE
      IF(.NOT.IHIT(KL,2)) GO TO 110
      IF(.NOT.IHIT(KL,1)) THEN
        IFLG(K)=5
        KEY=5
      ELSE
C     IF(KEY.EQ.5) GO TO 5
        IFLG(K)=KEY
C     KEY=0
      END IF
      GO TO 140
  110 CONTINUE
      IF(.NOT.IHIT(KL,3)) GO TO 120
      IF(.NOT.IHIT(KL,4)) THEN
        IFLG(K)=5
        KEY=5
      ELSE
        IFK=IFLG(K-NR)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.4)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
  120 CONTINUE
      IF(.NOT.IHIT(KL,4)) GO TO 130
      IF(.NOT.IHIT(KL,3)) THEN
        IFLG(K)=0
        KEY=0
      ELSE
        IFK=IFLG(K-NR)
        IF((IFK.EQ.0).OR.(IFK.EQ.1).OR.(IFK.EQ.4)) THEN
          IFLG(K)=0
          KEY=0
        ELSE
          IFLG(K)=5
          KEY=5
        END IF
      END IF
      GO TO 140
  130 CONTINUE
  140 CONTINUE
      J=J+1
      GO TO 50
  150 CONTINUE
C-----------------------------------------------------------------------
C JADGEMENT OF BOUNDARY CONDITIONS
      I=NZL
      DO 160 J=1,NR
        K=(I-LZF)*NR+J
        L=IFLG(K)
        IF(J.GT.1) THEN
          IF(I.GT.1) THEN
            LJ1=IFLG(K-1-NR)
          ELSE
            LJ1=IFLG(K-1)
          END IF
        ELSE
          LJ1=0
        END IF
C
        IF(I.EQ.1) THEN
          L1=IFLG(K)
        ELSE
          L1=IFLG(K-NR)
        END IF
C
        IF(J.EQ.1) THEN
          LJ=5
        ELSE
          LJ=IFLG(K-1)
        END IF
        IF(((L.EQ.0).OR.(L.EQ.2).OR.(L.EQ.3)).AND.((L1.EQ.0).OR.
     +  (L1.EQ.1) .OR.(L1.EQ.4))) SJDG(1,K)=1.
        IF(((L.EQ.0).OR.(L.EQ.2).OR.(L.EQ.4)).AND.((LJ.EQ.0).OR.
     +  (LJ.EQ.1) .OR.(LJ.EQ.3))) SJDG(3,K)=1.
        IF(((L.EQ.0).OR.(L.EQ.2)).AND.((L1.EQ.0).OR.(L1.EQ.4)).AND.
     +  ((LJ.EQ.0).OR.(LJ.EQ.3)).AND.((LJ1.EQ.0).OR.(LJ1.EQ.1))) SJDG
     +  (2,K)=1.
  160 CONTINUE
C
  170 CONTINUE
      RETURN
      END
C
      SUBROUTINE QKSET
C-----------------------------------------------------------------------
C SET VALUES OF CHARGE DENSITY ON MESH
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
      IBT=NINT(HLNGTH*2./DDZ)+2
      DO 10 N=1,NWN
   10 QKST(N)=0.0
      DO 20 N=1,IBT
        CALL CHAGE(1,(N-1)*M,QK)
        QKST(N)=QK
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE CHAGE(IZ,N,QK)
C-----------------------------------------------------------------------
C CALCULATION OF THE LINE CHARGE DENSITY
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /CALS/ SUMQ,BSIG
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      X=FLOAT(IZ-1)*DDZ-(FLOAT(N)*DDZ/FLOAT(M)-HLNGTH)
      IF(ABS(X).GT.HLNGTH)   THEN
        QK=0.
      ELSE
        QK=BSHAPE(X)/(SUMQ*DDZ)
      END IF
      RETURN
      END
C
      SUBROUTINE CALSUM
C-----------------------------------------------------------------------
C CALCULATION OF THE SUM OF LINE CHARGE DENSITY
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /CALS/ SUMQ,BSIG
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      I8SD=NINT(HLNGTH*2./DDZ)
      SUM=0.
      Z0=DDZ/FLOAT(M)-HLNGTH
      DO 10 I=-I8SD-NZ-1,I8SD+NZ+1
        X=FLOAT(I-1)*DDZ-Z0
        IF(ABS(X).GT.HLNGTH) GO TO 10
        SUM=SUM+BSHAPE(X)
   10 CONTINUE
      SUMQ=SUM
      RETURN
      END
C
      REAL FUNCTION BSHAPE(X)
C-----------------------------------------------------------------------
C BUNCH SHAPE (LINE CHARGE DISTRIBUTION) DEFAULT:GAUSSIAN
C USERS CAN CHANGE THE BUNCH SHAPE AS THEY WISH.
C RULES:
C  X=RELATIVE POSITION FROM THE CENTER OF THE CHARGE DISTRIBUTION.
C  HLNGTH=HALF LENGTH OF TOTAL DISTRIBUTION (DEFAULT:BSIG*SIG)
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /CALS/ SUMQ,BSIG
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      IF(NBUNCH-1) 40,10,20
C A SINGLE GAUSSIAN DISTRIBUTION:
   10 HLNGTH=BSIG*SIG
      BSHAPE=EXP(-X**2/(2.*SIG**2))
      RETURN
C GAUSSIAN BUNCHLETS
   20 HL=FLOAT(NBUNCH-1)/2.*BSEP
      HLNGTH=HL+BSIG*SIG
      XX=-HL-BSEP
      SUM=0.
      DO 30 N=1,NBUNCH
        XX=XX+BSEP
        ARG=(X-XX)/SIG
        IF(ABS(ARG).LE.5.) SUM=SUM+EXP(-0.5*ARG**2)
   30 CONTINUE
      BSHAPE=SUM
      RETURN
   40 CONTINUE
C *** USER MUST PROVIDE HIS CHARGE DISTRIBUTION BELOW ***
C HERE IS AN EXAMPLE FOR A TRIANGULAR DISTRIBUTION
      HLNGTH=SIG
      IF(ABS(X).LE.SIG) BSHAPE=SIG-ABS(X)
      IF(ABS(X).GT.SIG) BSHAPE=0.0
      RETURN
      END
C
      SUBROUTINE FLOT(IA,IK,X,IX,IERROR)
C-----------------------------------------------------------------------
C DECODE INPUT LINES OF MESHES
C-----------------------------------------------------------------------
      PARAMETER (NVR=100)
      DIMENSION X(2)
      DIMENSION IX(4)
      LOGICAL LL,LVAR,LCE,LREP
      CHARACTER*1 POS,MIN,DEC,E,D,IA(IK)
      CHARACTER*1 COM,RR,ZZ
      CHARACTER*1 NAM(11)
      CHARACTER*1 CL,LG,SH,EL
      CHARACTER*1 AIN,AOU
      DATA POS/'+'/,MIN/'-'/,DEC/'.'/,E/'E'/,D/'D'/
      DATA COM/','/,RR/'R'/,ZZ/'Z'/
      DATA NAM/'0','1','2','3','4','5','6','7','8','9',' '/
      DATA CL/'C'/,LG/'L'/,SH/'S'/,EL/'E'/
      DATA AIN/'>'/,AOU/'<'/
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /VAR/ RZ(NVR)
      X(1)=0.0
      X(2)=0.0
      IX(1)=0
      IX(2)=0
      IX(3)=0
      IX(4)=-1
      IKS=0
      LREP=.FALSE.
      IERROR=0
      DO 120 N=1,2
        LL=.FALSE.
        LCE=.FALSE.
        LVAR=.FALSE.
        R=0.0
        IE=0
        TEN=10.0
        IDEC=0
        SIG=1.0
        I=IKS
   10   I=I+1
        IF(I.GT.IK) GO TO 90
        IF(IA(I).EQ.NAM(11).AND.(.NOT.LL)) GO TO 10
        IF(IA(I).EQ.NAM(11).AND.(LL)) GO TO 90
        IF(IA(I).EQ.COM.AND.(.NOT.LL)) GO TO 10
        IF(IA(I).EQ.COM.AND.(LL)) GO TO 90
        IF(IA(I).EQ.POS) GO TO 10
        IF(IA(I).EQ.MIN) GO TO 20
        IF(IA(I).EQ.AOU) IX(4)=-9999
        IF(IA(I).EQ.AOU) RETURN
        IF((IA(I).EQ.E).AND.(LL)) GO TO 110
        IF((IA(I).EQ.D).AND.(LL)) GO TO 110
C     IF(IA(I).EQ.E.OR.IA(I).EQ.D) GO TO 25
        IF(IA(I).NE.DEC) GO TO 30
        IDEC=-1
        TEN=1.00
        GO TO 10
   20   SIG=-1.0
        GO TO 10
   30   DO 40 J=1,10
   40   IF(IA(I).EQ.NAM(J)) GO TO 80
        IF(IA(I).NE.RR) GO TO 50
        I=I+1
        IF(IA(I).NE.ZZ) GO TO 50
        LL=.TRUE.
        LVAR=.TRUE.
        GO TO 10
   50   CONTINUE
        IF(IA(I).NE.AIN) GO TO 60
C     LL=.TRUE.
        LREP=.TRUE.
        IKS=I
        GO TO 120
   60   CONTINUE
        IF(IA(I).EQ.CL) THEN
          I=I+1
          IF(IA(I).EQ.SH) THEN
            X(1)=-1.0
            IX(3)=-1
            LL=.TRUE.
            LCE=.TRUE.
            GO TO 10
          ELSE
          END IF
          IF(IA(I).EQ.LG) THEN
            X(1)=-2.0
            IX(3)=-2
            LL=.TRUE.
            LCE=.TRUE.
            GO TO 10
          ELSE
          END IF
          GO TO 70
        ELSE
        END IF
        IF(IA(I).EQ.EL) THEN
          I=I+1
          IF(IA(I).EQ.SH) THEN
            X(1)=-3.0
            IX(3)=-3
            LL=.TRUE.
            LCE=.TRUE.
            GO TO 10
          ELSE
          END IF
          IF(IA(I).EQ.LG) THEN
            X(1)=-4.0
            IX(3)=-4
            LL=.TRUE.
            LCE=.TRUE.
            GO TO 10
          ELSE
          END IF
          GO TO 70
        ELSE
        END IF
   70   WRITE(ND6,10000) IA(I),(IA(K),K=1,IK)
        IERROR=-1
        RETURN
   80   S=J-1
        IF(IDEC.NE.0) S=S*10D0**IDEC
        R=R*TEN+S
        IF(IDEC.NE.0) IDEC=IDEC-1
        LL=.TRUE.
        GO TO 10
   90   CONTINUE
        IF(LCE) GO TO 100
        X(N)=R*SIG*10D0**IE
        IF(.NOT.LVAR) GO TO 100
        IX(N)=X(N)
        IXN=ABS(IX(N))
        X(N)=SIG*RZ(IXN)
        IF(X(N).NE.-8888.) GO TO 100
        WRITE(ND6,10010) IX(N),(IA(K),K=1,IK)
  100   IKS=I
        IF(LREP) IX(4)=X(N)
C     WRITE(6,1010) X(N)
C1010 FORMAT(1X,'X=',G13.5)
        GO TO 120
  110   LK=IK-I
        IE=IPAK(IA(I+1),LK,IL)
        I=I+IL
        GOTO 90
  120 CONTINUE
      RETURN
10000 FORMAT(' *** ERROR *** ILLEGAL CHARACTER "',A1,
     +'" IN MESH INPUT LINE --> ',72A1)
10010 FORMAT(' *** ERROR *** UNGIVEN VARIABLE "RZ',I2,
     +'" IN MESH INPUT LINE --> ',40A1)
      END
C
      FUNCTION IPAK(IA,IK,IL)
C-----------------------------------------------------------------------
C DECODE EXPONENTIAL PART
C-----------------------------------------------------------------------
      CHARACTER*1 IA(IK)
      CHARACTER*1 NAM(10),M,P,BLK,COM
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      DATA NAM/'0','1','2','3','4','5','6','7','8','9'/, M/'-'/,P/'+'/,
     +BLK/' '/,COM/','/
      J=0
      IS=1
      DO 40 L=1,IK
        IF(IA(L).EQ.P) GOTO 40
        IF(IA(L).EQ.BLK) GOTO 50
        IF(IA(L).EQ.COM) GOTO 50
        IF(IA(L).NE.M) GO TO 10
        IS=-1
        GOTO 40
   10   DO 20 K=1,10
   20   IF(IA(L).EQ.NAM(K)) GOTO 30
        WRITE(ND6,10000) IA(L)
        STOP
   30   J=J*10+K-1
   40 CONTINUE
   50 IPAK=J*IS
      IL=L
      RETURN
10000 FORMAT(' *** ERROR *** ILLEGAL CHARACTER "',A1,'" ')
      END
C
      SUBROUTINE IRDSET(RAD)
C-----------------------------------------------------------------------
C SET ARRAYS FR AND FDR
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM
      DIMENSION DRR(NF)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /LAT3/ IR(NF)
      COMMON /SYMT/ SYM
      RMAX=0.0
      JMAX=1
      ID=0
      IX=0
      DO 10 I=1,NF
        IF(RMAX.LT.RMARK(I)) JMAX=I
        IF(RMAX.LT.RMARK(I)) RMAX=RMARK(I)
        IF(RMARK(I).NE.0.) IX=IX+1
        DRR(I)=DDR(I)
   10 IF(DDR(I).NE.0.) ID=ID+1
C
      IF(DDZ.EQ.0.) THEN
        DRR(1)=RAD/FLOAT(NR-1)
        RBDR(1)=RAD
      ELSE
        IF(RMAX.LT.RAD) RMARK(JMAX)=RAD
        IF(ID.EQ.1) THEN
          RBDR(1)=RAD
          IF((IX.NE.0).AND.(IX.NE.1)) GO TO 80
        ELSE
          RBDR(1)=RMARK(1)
          IF(ID.NE.IX) GO TO 80
        END IF
      END IF
C
      I=1
      IR(1)=NINT(RBDR(1)/DRR(1)-0.00001)+1
      RBDR(1)=FLOAT(IR(1)-1)*DRR(1)
      IF(RMARK(1).EQ.RAD) GO TO 30
      IF(NZ.NE.0) GO TO 30
      IF(I.EQ.ID) GO TO 30
   20 I=I+1
      IF(RMARK(I).EQ.0.) GO TO 30
      IF(DRR(I).LE.0.) GO TO 70
      IR(I)=NINT((RMARK(I)-RBDR(I-1))/DRR(I)-0.00001)+IR(I-1)
      RBDR(I)=FLOAT(IR(I)-IR(I-1))*DRR(I)+RBDR(I-1)
      IF(RMARK(I).EQ.RAD) GO TO 30
      IF(I.EQ.ID) GO TO 30
      GO TO 20
   30 CONTINUE
      IMAX=I
      IRMAX=IR(IMAX)
      IF(IRMAX.GE.NRN) THEN
        WRITE(ND6,10020) IRMAX,NRN
        STOP
      END IF
      DO 60 I=1,IRMAX
        J=1
        IF((I.GT.0).AND.(I.LE.IR(J))) GO TO 50
   40   J=J+1
        IF((I.GT.IR(J-1)).AND.(I.LE.IR(J))) GO TO 50
        GO TO 40
   50   CONTINUE
        GDR(I)=DRR(J)
        IF(J.EQ.1) GR(I)=FLOAT(I-1)*DRR(1)
        IF(J.GT.1) GR(I)=RBDR(J-1)+FLOAT(I-IR(J-1))*DRR(J)
   60 CONTINUE
      RETURN
   70 WRITE(ND6,10000) I
      STOP
   80 WRITE(ND6,10010) ID,IX
      STOP
10000 FORMAT(1X,'*** STOP *** THE VALUE OF THE', I3,
     +'TH DDR IS NOT GIVEN.')
10010 FORMAT(1X,'*** STOP *** THE NUMBER OF NON-ZERO DDR=',I2,
     +' DOES NOT MATCH TO THE NUMBER OF NON-ZERO RMARK=',I2)
10020 FORMAT(1X,'*** STOP ***  NR:',I7,' EXCEEDS MAXIMUM NRN:',I7)
      END
C
      INTEGER FUNCTION IRD(R)
C-----------------------------------------------------------------------
C RETURN J COORDINATE FOR R
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /LAT3/ IR(NF)
      IF(R.GT.RBDR(IMAX)) R=RBDR(IMAX)
      IF(R.LT.0.) R=0.0
      J=1
      IF((R.GE.0.).AND.(R.LE.RBDR(J))) GO TO 20
   10 J=J+1
      IF((R.GT.RBDR(J-1)).AND.(R.LE.RBDR(J))) GO TO 20
      GO TO 10
   20 CONTINUE
      IF(J.GT.1) IRD=NINT((R-RBDR(J-1))/DDR(J)-0.00001)+IR(J-1)
      IF(J.EQ.1) IRD=NINT(R/DDR(J)-0.00001)+1
      RETURN
      END
C
      SUBROUTINE IRDINI
C-----------------------------------------------------------------------
C SET ARRAYS FR AND FDR
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /LAT3/ IR(NF)
C
      COMMON / CNTCOM / NRCNT4, NRCNT7, NRCNTM
      SAVE / CNTCOM /
      COMMON / NAPCOM / NRNAP1, NRNAP2
      SAVE / NAPCOM /
C
      NRCNT4 = 0
      NRCNT7 = 0
      NRCNTM = 0
      NRNAP1 = 0
      NRNAP2 = 0
      NR=0
      NZ=0
      DDZ=0.0
      IMAX=0.0
      IRMAX=0.0
      DO 10 I=1,NF
        DDR(I)=0.0
        IR(I)=0.0
        RMARK(I)=0.0
   10 RBDR(I)=0.0
      DO 20 I=1,NRN
        GR(I)=0.0
   20 GDR(I)=0.0
      RETURN
      END
C
      SUBROUTINE DRPLOT (ND6,X,Y,N,KSYMB,KEQUAL,XYMM,KSTORE)
C-----------------------------------------------------------------------
C  PRINT-PLOT SIGULATION FOR ONE OR SEVERAL CURVES ON THE LINEPRINTER
C  ORIGINAL VERSION: H.BRANDIS, TH-DARMSTADT
C-----------------------------------------------------------------------
      REAL Y(N)
      DIMENSION X(1),XYMM(4),INFO(690), IPOT(10), TV(10),XB(10),IZEICH
     +(7)
      CHARACTER*1 LZM,LZP,LZI,LZB,IZIF(10),KSYMB,KPRI(121),LZEA
      CHARACTER*4 KFMT(7)
      CHARACTER*3 IFT(7)
C LSP LZE DEFINE SIZE OF PRINT IN COLUMNS/LINES
      DATA TV / 0.2,0.5,0.5,1.,1.,1.,1.,2.,2.,2.  /
      DATA LZM,LZP,LZI,LZB, KBL,IZL / '-','+','I',' ' ,  0,0  /
      DATA KFMT / '(1H0', ',6X,', '10(1', 'PE9.', '2,  ', ' ', '))' /
      DATA IFT / '21X', '15X', '11X', '8X', '6X', '4X', '3X'  /
      DATA IZIF / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
      DATA LSP, LZE, LST, NDIM / 121, 42,  9, 690 /
C
      LZEI=ICHAR(KSYMB)
      IF (LZEI.LT.10.AND.LZEI.GE.0) LZEI=ICHAR(IZIF(LZEI+1))
      IF (KBL) 470,10,160
   10 KBL=1
      SPZE=FLOAT(LSP+LZE)
      LRX=LSP-1
      LRY=LZE-1
      IPOT(1)=1
      DO 20 I=2,LST
        IPOT(I)=IPOT(I-1)*8
   20 CONTINUE
      DO 30 I=1,NDIM
        INFO(I)=0
   30 CONTINUE
      DF=XYMM(2)-XYMM(1)
      IF (DF-5.E-4*ABS(XYMM(2)+XYMM(1))) 40,40,60
   40 IF (KEQUAL.EQ.0) GO TO 460
      XYMM(1)=X(1)
      XYMM(2)=X(1)
      DO 50 I=1,N
        XYMM(1)=AMIN1(XYMM(1),X(I))
        XYMM(2)=AMAX1(XYMM(2),X(I))
   50 CONTINUE
      DF=XYMM(2)-XYMM(1)
      IF (DF-5.E-4*ABS(XYMM(2)+XYMM(1))) 460,460,60
   60 P=10.**(INT(ALOG10(DF)+199.999)-200)
      IST=DF/P
   70 XT=TV(IST)*P
      NTP=INT(XYMM(2)/XT+AMAX1(SIGN(0.99,XYMM(2)),0.01))
      NTM=INT(XYMM(1)/XT+AMIN1(SIGN(0.99,XYMM(1)),0.01))
      NTX=NTP-NTM
      IF (NTX-8) 90,90,80
   80 IST=IST+1
      GO TO 70
   90 ITX=LRX/NTX
      NT0=-NTM+1
      XMIN=XT*FLOAT(NTM)
      DX=XT/FLOAT(ITX)
      KPX=ITX*NTX+1
      KFMT(6)=IFT(NTX-3)
      DF=XYMM(4)-XYMM(3)
      IF (DF-5.E-4*ABS(XYMM(4)+XYMM(3))) 100,100,120
  100 XYMM(3)=Y(1)
      XYMM(4)=Y(1)
      DO 110 I=1,N
        XYMM(3)=AMIN1(XYMM(3),Y(I))
        XYMM(4)=AMAX1(XYMM(4),Y(I))
  110 CONTINUE
      DF=XYMM(4)-XYMM(3)
      IF (DF-5.E-4*ABS(XYMM(4)+XYMM(3))) 460,460,120
  120 P=10.**(INT(ALOG10(DF)+199.999)-200)
      IST=DF/P
  130 YT=TV(IST)*P
      NTP=INT(XYMM(4)/YT+AMAX1(SIGN(0.99,XYMM(4)),0.01))
      NTM=INT(XYMM(3)/YT+AMIN1(SIGN(0.99,XYMM(3)),0.01))
      NTY=NTP-NTM
      IF (NTY-8) 150,150,140
  140 IST=IST+1
      GO TO 130
  150 ITY=LRY/NTY
      NZE0=ITY*NTP+1
      YMAX=YT*FLOAT(NTP)
      DY=YT/FLOAT(ITY)
      KPY=ITY*NTY+1
      GO TO 190
  160 IF (ICHAR(KSYMB)+1) 170,200,170
  170 DO 180 KZL=1,KBL
        IF (LZEI.EQ.IZEICH(KZL))GO TO 200
  180 CONTINUE
      KBL=MOD(KBL,7)+1
  190 KZL=KBL
      LZEA = CHAR(LZEI)
      IF (LZEA.EQ.LZB)LZEI=ICHAR(IZIF(KBL+1))
      IZEICH(KZL)=LZEI
  200 IF (KEQUAL) 220,210,220
  210 DXX=(XYMM(2)-XYMM(1))/DX/AMAX0(N-1,1)
      XX=(XYMM(1)-XMIN)/DX+1.5-DXX
  220 IF (N) 310,310,230
  230 DO 300 I=1,N
        IF (KEQUAL) 250,240,250
  240   SP=DXX*FLOAT(I)+XX
        GO TO 260
  250   SP=(X(I)-XMIN)/DX+1.5
  260   ZE=(YMAX-Y(I))/DY+1.5
        IF (ABS(SP)+ABS(ZE)-SPZE) 270,300,300
  270   NSP=INT(SP)-1
        NZE=INT(ZE)-1
        IF (NSP*(LRX-NSP)) 300,280,280
  280   IF (NZE*(LRY-NZE)) 300,290,290
  290   NPKT=NZE*LSP+NSP
        NAD=NPKT/LST
        NEX=MOD(NPKT,LST)
        INFO(NAD+1)=INFO(NAD+1)+IPOT(NEX+1)*(KZL-MOD(INFO(NAD+1)/
     +  IPOT(NEX+1),8))
  300 CONTINUE
  310 IF (KSTORE) 771,320,771
  320 DO 440 IZE=1,KPY
        KR=MOD(IZE-1,ITY)
        IF (KR) 360,330,360
  330   DO 340 ISP=1,KPX
          KPRI(ISP)=LZM
  340   CONTINUE
        DO 350 ISP=1,KPX,ITX
          KPRI(ISP)=LZP
  350   CONTINUE
        GO TO 390
  360   DO 370 ISP=1,KPX
          KPRI(ISP)=LZB
  370   CONTINUE
        DO 380 ISP=1,KPX,ITX
          KPRI(ISP)=LZI
  380   CONTINUE
  390   NPA=(IZE-1)*LSP-1
        DO 410 ISP=1,KPX
          IPK=ISP+NPA
          NAD=IPK/LST
          NEX=MOD(IPK,LST)
          IZL=MOD(INFO(NAD+1)/IPOT(NEX+1),8)
          IF (IZL) 410,410,400
  400     KPRI(ISP)=CHAR(IZEICH(IZL))
  410   CONTINUE
        IF (KR) 430,420,430
  420   YB=DY*FLOAT(NZE0-IZE)
        WRITE (ND6,10000) YB,(KPRI(ISP),ISP=1,KPX)
        GO TO 440
  430   WRITE (ND6,10010) (KPRI(ISP),ISP=1,KPX)
  440 CONTINUE
      DO 450 IT=1,NTX
        XB(IT)=XT*FLOAT(IT-NT0)
  450 CONTINUE
      WRITE (ND6,KFMT) (XB(IT),IT=1,NTX)
      KBL=0
      GO TO 771
  460 WRITE (ND6,10020) XYMM
  470 KBL=-IABS(KSTORE)
  771 RETURN
C
10000 FORMAT (1PE10.2, 1X, 121A1)
10010 FORMAT (11X, 121A1)
10020 FORMAT(///1X,11H**WARNING**//1X,
     +'WRONG RANGE IN X OR Y GIVEN TO SUBROUTINE DRPLOT:'/ ' XMIN=',1
     +PE12.5/ ' XMAX=',1PE12.5/ ' YMIN=',1PE12.5/ ' YMAX=',1PE12.5//)
      END
C
      SUBROUTINE CPLOT(ZTXT)
C-----------------------------------------------------------------------
C TYPE A HEADLINE FOR A PLOT OF CAVITY SHAPE INPUT ON FT09
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LNAPOLY,LPATH,LNONAP
      CHARACTER*13 REDCHR
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /WPT/ NDATE,NTIME
      COMMON /RED/ REDCHR
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON /PLNZ/ LPATH
      COMMON /PINT/ KMIN
      COMMON /NZPATH/ ZFST,ZLST,RFST,RLST,XMAX,YMAX
      DATA SCM /1000./
      IMAX=1
   10 IF(DDR(IMAX+1).EQ.0.) GO TO 20
      IMAX=IMAX+1
      GO TO 10
   20 IF(IMAX.GT.8) IMAX=6
      WRITE(ND9,10000) REDCHR,(ZTXT(J),J=1,3),(ZTXT(L),L=4,9)
      WRITE(ND9,10010) NDATE,NTIME
      WRITE(ND9,10020) DDZ*SCM,DDR(1)*SCM
      IF(IMAX.GE.2) THEN
        DO 30 I=2,IMAX
   30   WRITE(ND9,10030) DDR(I)*SCM
      END IF
      CALL CAVPLT(YMIN)
C
      IF (LPATH) THEN
        ZFST=Z(2)
        RFST=R(2)
        ZLST=Z(NVC)
        RLST=R(NVC)
        YDEV=YMAX/200.
        XDEV=XMAX/200.
        IF(LNAPOLY) THEN
          ZCF=FLOAT(IZCF-1)*DDZ
          ZCT=FLOAT(IZCT-1)*DDZ
          RWAK=GR(IRWAK)
          WRITE(ND9,10040) ZFST,RFST-YDEV
          WRITE(ND9,10040) ZCF-XDEV,RFST-YDEV
          WRITE(ND9,10040) ZCF-XDEV,RWAK-YDEV
          WRITE(ND9,10040) ZCT+XDEV,RWAK-YDEV
          WRITE(ND9,10040) ZCT+XDEV,RLST-YDEV
          WRITE(ND9,10040) ZLST,RLST-YDEV
          WRITE(ND9,10050)
       ELSE
          RLINE=GR(KMIN)
          WRITE(ND9,10040) ZFST,RLINE-YDEV
          WRITE(ND9,10040) ZLST,RLINE-YDEV
          WRITE(ND9,10050)
        END IF
      END IF
C
      RETURN
10000 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/
     +1X,' SET WINDOW X FROM 2 12 Y FROM 1 TO 8'/
     +1X,' SET TICKS SIZE 0.05'/ 1X,' TITLE 7.0 9.0 CENTER SIZE 2.2  ',
     +1H', 'Cavity Shape Input',1H'/ 1X,' TITLE 2.0 8.6 SIZE 1.5  ',
     +1H','A B C I 9.4 : ',3A8,1H'/ 1X,' MORE ',1H',6A8,1H'/ 1X,
     +' TITLE BOTTOM ',1H','Z-axis (m)',1H'/
     +1X,' TITLE LEFT   ',1H','R-axis (m)',
     +1H')
10010 FORMAT(1X,' TITLE 10.9 8.9 SIZE 1.25 ',1H',A8,2X,A8,1H')
10020 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H', 'DDZ=',F6.3,' mm,',
     +2X,'DDR=',F6.3,' mm',1H')
10030 FORMAT(1X,' MORE',1X,1H',', ',F6.3,' mm',1H')
10040 FORMAT(1X,2G13.5)
10050 FORMAT(1X,' JOIN 1 DASH')
      END
C
      SUBROUTINE EPLOT(ZTXT,NPLOT,NN,LPLE,LPALL)
C-----------------------------------------------------------------------
C PLOT ELECTRIC FIELD LINES ON FT09
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NSP=20000)
      PARAMETER (NWN=500000)
      REAL Z(NMESH)
      CHARACTER*13 REDCHR
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      LOGICAL LPLE,LPALL,LNOPLT
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /CDIS/ HLNGTH,BSEP,NBUNCH
      COMMON /WPT/ NDATE,NTIME
      COMMON /RED/ REDCHR
      COMMON /WIND/ NT
      COMMON /WPL/ MROT
      COMMON /CNTR/ IX,IY,IDX,IDY,IS,ICL,VC,IMX(8),IMY(8),IQ(NSP),NP
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
      DATA SCM /100./
C SET UP PARAMETERS
      LNOPLT=.FALSE.
      ZMIN=1.E30
      ZMAX=-ZMIN
      FAC1=2.0*PI*DDZ/FLOAT(M)/CLGHT
      FAC2=2.0*PI
      IF(.NOT.LPLE) THEN
        DO 10 I=1,NZL
          DO 10 J=1,NR
            K=(I-1)*NR+J
            IF(IFLG(K).EQ.5) GO TO 10
            Z(K)=0.0
            IF(ABS(HV(2,K)).LE.CEPS) GO TO 10
            Z(K)=HV(2,K)*(GR(J)+GDR(J)/2.)*FAC2
            IF(Z(K).GT.ZMAX) ZMAX=Z(K)
            IF(Z(K).LT.ZMIN) ZMIN=Z(K)
   10   CONTINUE
      ELSE
        DO 20 I=1,NZL
          DO 20 J=1,NR
            K=(I-1)*NR+J
            IF(IFLG(K).EQ.5) GO TO 20
            Z(K)=0.0
            IF(ABS(HV(3,K)).LE.CEPS) GO TO 20
            Z(K)=HV(3,K)*(GR(J)+GDR(J)/2.)*FAC1
            IF(Z(K).GT.ZMAX) ZMAX=Z(K)
            IF(Z(K).LT.ZMIN) ZMIN=Z(K)
   20   CONTINUE
      END IF
      IF(LPALL) THEN
        IF(LPLE) NCON=NPLOT
        IF(.NOT.LPLE) NCON=NPLOT/2
        ZDEF=(ZMAX-ZMIN)/FLOAT(NCON)
        ZMIN=ZMIN+ZDEF*0.5
        ZMAX=ZMAX-ZDEF*0.5
      ELSE
        IF(ABS(ZMAX).LE.CEPS.AND.ABS(ZMIN).LE.CEPS) LNOPLT=.TRUE.
        IF(LPLE) ZDEF=1.0/FLOAT(NPLOT)/FLOAT(NBUNCH)
        IF(.NOT.LPLE) ZDEF=QKMAX*CLGHT/FLOAT(NPLOT/2)
        IF(ZMAX.GT.0.0) NCON1=ZMAX/ZDEF+0.5
        IF(ZMAX.LT.0.0) NCON1=ZMAX/ZDEF-0.5
        IF(ZMAX.EQ.0.0) NCON1=0
        IF(ZMIN.GT.0.0) NCON2=ZMIN/ZDEF+0.5
        IF(ZMIN.LT.0.0) NCON2=ZMIN/ZDEF-0.5
        IF(ZMIN.EQ.0.0) NCON2=0
        NCON=NCON1-NCON2
        ZMAX=FLOAT(NCON1)*ZDEF-ZDEF*0.5
        ZMIN=FLOAT(NCON2)*ZDEF+ZDEF*0.5
      END IF
      IF(NCON.LT.1) LNOPLT=.TRUE.
C PRINT TITLE
      IF(LPLE) THEN
        WRITE(ND9,10000) REDCHR,(ZTXT(J),J=1,3),(ZTXT(L),L=4,9)
      ELSE
        WRITE(ND9,10010) REDCHR,(ZTXT(J),J=1,3),(ZTXT(L),L=4,9)
      END IF
      TT=FLOAT(NN)/FLOAT(M)*DDZ/CLGHT*1.E12
      WRITE(ND9,10020) NDATE,NTIME
      WRITE(ND9,10030) TT
      IF(LPLE) THEN
        WRITE(ND9,10040) SIG*SCM,ZMIN,ZMAX,ZDEF
      ELSE
        WRITE(ND9,10050) SIG*SCM,ZMIN,ZMAX,ZDEF
      END IF
C PLOT A CAVITY SHAPE
      CALL CAVPLT(YMIN)
C PLOT A CHARGE DISTRIBUTION
      IF(NN.LT.NT) CALL CHGPLT(NN,YMIN)
C PLOT A CONTOUR MAP
      IF(LNOPLT) RETURN
      DO 110 N=1,NCON
        VC=ZMIN+ZDEF*FLOAT(N-1)
        NP=0
        ICL=0
        DO 40 IX1=2,NZL
          I=IX1-1
          DO 40 J=1,NR
            K=(I-1)*NR+J
            K1=(IX1-1)*NR+J
            K2=K1-1
            IF(J.EQ.1) K2=K1
            K4=K-1
            IF(J.EQ.1) K4=K
            K3=K+1
            IF(J.EQ.NR) K3=K
            K5=K1+1
            IF(J.EQ.NR) K5=K1
            IF(IFLG(K).EQ.5.OR.IFLG(K1).EQ.5) GO TO 40
            IF(Z(K).GE.VC.OR.Z(K1).LT.VC) GO TO 30
            IF(J.EQ.1.OR.IFLG(K2).EQ.5.OR.IFLG(K2).EQ.2.OR.IFLG(K2).EQ.
     +      4.OR. IFLG(K4).EQ.5.OR.IFLG(K4).EQ.2.OR.IFLG(K4).EQ.4.OR.
     +      IFLG(K).EQ.1.OR.IFLG(K).EQ.3.OR.IFLG(K1).EQ.1.OR.IFLG(K1)
     +      .EQ.3) THEN
              IX=IX1
              IY=J
              IDX=-1
              IDY=0
              IS=1
              CALL CONPNT(Z)
            ELSE
            END IF
   30       IF(Z(K1).GE.VC.OR.Z(K).LT.VC) GO TO 40
            IF(J.EQ.NR.OR.IFLG(K3).EQ.5.OR.IFLG(K3).EQ.1.OR.IFLG(K3).EQ.
     +      3.OR. IFLG(K5).EQ.5.OR.IFLG(K5).EQ.1.OR.IFLG(K5).EQ.3.OR.
     +      IFLG(K).EQ.2.OR.IFLG(K).EQ.4.OR.IFLG(K1).EQ.2.OR.IFLG(K1)
     +      .EQ.4) THEN
              IX=I
              IY=J
              IDX=1
              IDY=0
              IS=5
              CALL CONPNT(Z)
            ELSE
            END IF
   40   CONTINUE
        DO 60 JY1=2,NR
          J=JY1-1
          DO 60 I=1,NZL
            K=(I-1)*NR+J
            K1=(I-1)*NR+JY1
            K2=K1+NR
            IF(I.EQ.NZ) K2=K1
            K4=K+NR
            IF(I.EQ.NZ) K4=K
            K3=K-NR
            IF(I.EQ.1) K3=K
            K5=K1-NR
            IF(I.EQ.1) K5=K1
            IF(IFLG(K).EQ.5.OR.IFLG(K1).EQ.5) GO TO 60
            IF(Z(K).GE.VC.OR.Z(K1).LT.VC) GO TO 50
            IF(I.EQ.NZL.OR.IFLG(K2).EQ.5.OR.IFLG(K2).EQ.1.OR.IFLG(K2)
     +      .EQ.4.OR. IFLG(K4).EQ.5.OR.IFLG(K4).EQ.1.OR.IFLG(K4).EQ.4.OR
     +      . IFLG(K).EQ.2.OR.IFLG(K).EQ.3.OR.IFLG(K1).EQ.2.OR.IFLG(K1)
     +      .EQ.3) THEN
              IX=I
              IY=JY1
              IDX=0
              IDY=-1
              IS=7
              CALL CONPNT(Z)
            ELSE
            END IF
   50       IF(Z(K1).GE.VC.OR.Z(K).LT.VC) GO TO 60
            IF(I.EQ.1.OR.IFLG(K3).EQ.5.OR.IFLG(K3).EQ.2.OR.IFLG(K3).EQ.
     +      3.OR. IFLG(K5).EQ.5.OR.IFLG(K5).EQ.2.OR.IFLG(K5).EQ.3.OR.
     +      IFLG(K).EQ.1.OR.IFLG(K).EQ.4.OR.IFLG(K1).EQ.1.OR.IFLG(K1)
     +      .EQ.4) THEN
              IX=I
              IY=J
              IDX=0
              IDY=1
              IS=3
              CALL CONPNT(Z)
            ELSE
            END IF
   60   CONTINUE
C
        ICL=1
        DO 100 JY1=3,NR
          J=JY1-1
          DO 90 IX1=2,NZL
            I=IX1-1
            K=(I-1)*NR+J
            K1=(IX1-1)*NR+J
            IF(IFLG(K).EQ.5.OR.IFLG(K1).EQ.5) GO TO 90
            IF(Z(K).GE.VC.OR.Z(K1).LT.VC) GO TO 90
            IF(NP.EQ.0) GO TO 80
            DO 70 L=1,NP
              IF(IQ(L).EQ.K1) GO TO 90
   70       CONTINUE
   80       NP=NP+1
            IF(NP.GT.NSP) RETURN
            IQ(NP)=K1
            IX=IX1
            IY=J
            IDX=-1
            IDY=0
            IS=1
            CALL CONPNT(Z)
   90     CONTINUE
  100   CONTINUE
  110 CONTINUE
C
      RETURN
10000 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/
     +1X,' SET WINDOW X FROM 2 12 Y FROM 1 TO 8'/
     +1X,' SET TICKS SIZE 0.05'/ 1X,' TITLE 7.0 9.0 CENTER SIZE 2.2  ',
     +1H', 'Electric Field Lines',1H'/ 1X,
     +' TITLE 2.0 8.6 SIZE 1.5  ',1H','A B C I 9.4 : ',3A8,1H'/ 1X,
     +' MORE ' ,1H',6A8,1H'/
     +1X,' TITLE BOTTOM ',1H','Z-axis (m)',1H'/ 1X,
     +' TITLE LEFT   ',1H','R-axis (m)',1H')
10010 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/
     +1X,' SET WINDOW X FROM 2 12 Y FROM 1 TO 8'/
     +1X,' SET TICKS SIZE 0.05'/ 1X,' TITLE 7.0 9.0 CENTER SIZE 2.2  ',
     +1H', 'Total Current Lines',1H'/ 1X,' TITLE 2.0 8.6 SIZE 1.5  '
     +,1H','A B C I 9.4 : ',3A8,1H'/ 1X,' MORE ',1H',6A8,1H'/ 1X,
     +' TITLE BOTTOM ',1H','Z-axis (m)',1H'/
     +1X,' TITLE LEFT   ',1H','R-axis (m)',
     +1H')
10020 FORMAT(1X,' TITLE 10.9 8.9 SIZE 1.25 ',1H',A8,2X,A8,1H')
10030 FORMAT(1X,' TITLE 10.1 9.1 SIZE 1.25 ',1H','Time Passed:', F11.2,
     +' (ps)',1H')
10040 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H',
     +'SIG=',F6.3,' cm,',1H'/
     +1X,' MORE',1X,1H',2X, 'Min/Max=',1PE10.3,'/',1PE10.3,'(C),',1H'/
     +1X,' MORE',1X,1H',2X, 'Flux Between Lines=',1PE10.3,'(C)',1H')
10050 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H',
     +'SIG=',F6.3,' cm,',1H'/
     +1X,' MORE',1X,1H',2X, 'Min/Max=',1PE10.3,'/',1PE10.3,'(A),',1H'/
     +1X,' MORE',1X,1H',2X, 'Current Between Lines=',1PE10.3,'(A)',1H')
      END
C
      SUBROUTINE CONPNT(Z)
C-----------------------------------------------------------------------
C PLOT ELECTRIC FIELD LINES ON FT09
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NSP=20000)
      REAL Z(NMESH)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /CNTR/ IX,IY,IDX,IDY,IS,ICL,VC,IMX(8),IMY(8),IQ(NSP),NP
      COMMON /ELINE/ NZF,NZL
C
      KXY(IXX,IYY)=(IXX-1)*NR+IYY
      F(P,Q)=(P-VC)/(P-Q)
C
      IX0=IX
      IY0=IY
      IS0=IS
      IN=0
      ILINE=1
C
      IF(IDX.EQ.0) GO TO 10
      Y=GR(IY)+GDR(IY)/2.
      IDD=IX+IDX
      X=F(Z(KXY(IX,IY)),Z(KXY(IDD,IY)))*FLOAT(IDX)*DDZ +(FLOAT(IX-1)
     ++0.5)*DDZ
      GO TO 20
   10 X=(FLOAT(IX-1)+0.5)*DDZ
      IDD=IY+IDY
      Y=F(Z(KXY(IX,IY)),Z(KXY(IX,IDD)))*(GDR(IY)+GDR(IDD))/2. *FLOAT
     +(IDY)+GR(IY)+GDR(IY)/2.
   20 WRITE(ND9,10000) X,Y
      ILINE=ILINE+1
      IF(IS.NE.1) GO TO 30
      NP=NP+1
      IF(NP.GT.NSP) GO TO 100
      IQ(NP)=KXY(IX,IY)
   30 IS=IS+1
      IN=IN+1
      IF(IN.GT.8) GO TO 100
      IF(IS.GT.8) IS=IS-8
      IDX=IMX(IS)
      IDY=IMY(IS)
      IX2=IX+IDX
      IY2=IY+IDY
      IF(IX2.GT.NZL.OR.IY2.GT.NR.OR.IX2.LT.1.OR.IY2.LT.1) GO TO 100
      IF(IFLG(KXY(IX2,IY2)).EQ.5) GO TO 100
   40 IF(VC-Z(KXY(IX2,IY2))) 50, 50, 60
   50 IS=IS+4
      IX=IX2
      IY=IY2
      GO TO 30
   60 IF(IS/2*2.EQ.IS) GO TO 30
      IF(IDX.EQ.0) GO TO 70
      Y=GR(IY)+GDR(IY)/2.
      IDD=IX+IDX
      X=F(Z(KXY(IX,IY)),Z(KXY(IDD,IY)))*FLOAT(IDX)*DDZ +(FLOAT(IX-1)
     ++0.5)*DDZ
      JX=IDD
      JY=IY
      GO TO 80
   70 X=(FLOAT(IX-1)+0.5)*DDZ
      IDD=IY+IDY
      Y=F(Z(KXY(IX,IY)),Z(KXY(IX,IDD)))*(GDR(IY)+GDR(IDD))/2. *FLOAT
     +(IDY)+GR(IY)+GDR(IY)/2.
      JX=IX
      JY=IDD
   80 WRITE(ND9,10000) X,Y
      IF(ILINE.EQ.1021) THEN
        WRITE(ND9,10010)
        WRITE(ND9,10000) X,Y
        ILINE=1
      END IF
      ILINE=ILINE+1
      XOLD=X
      YOLD=Y
      IN=0
      IF(IS.NE.1) GO TO 90
      NP=NP+1
      IF(NP.GT.NSP) GO TO 100
      IQ(NP)=KXY(IX,IY)
   90 IF(IX.EQ.IX0.AND.IY.EQ.IY0.AND.IS.EQ.IS0) GO TO 100
      IF(IFLG(KXY(JX,JY)).EQ.5) GO TO 100
      IF(ICL.EQ.0) GO TO 30
      IF(IX.NE.IX0.OR.IY.NE.IY0.OR.IS.NE.IS0) GO TO 30
  100 WRITE(ND9,10010)
      RETURN
10000 FORMAT(1X,2G13.5)
10010 FORMAT(1X,' JOIN 1 SOLID')
      END
C
      BLOCKDATA EPLBD
C-----------------------------------------------------------------------
C BLOCK DATA FOR CONPNT
C-----------------------------------------------------------------------
      PARAMETER (NSP=20000)
      COMMON /CNTR/ IX,IY,IDX,IDY,IS,ICL,VC,IMX(8),IMY(8),IQ(NSP),NP
      DATA IMX/-1,-1,0,1,1,1,0,-1/
      DATA IMY/0,1,1,1,0,-1,-1,-1/
      END
C
      SUBROUTINE CAVPLT(YMIN)
C-----------------------------------------------------------------------
C TYPE X-Y COORDINATES OF INPUT CAVITY SHAPE ON FT09 FOR TOP DRAWER
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM,LRZF
      DIMENSION AXS1(2),AXS2(2),AXS12(2),CENT(2),A1CENT(2)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      LOGICAL LLONG, LCAV,LVEX
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /SYMT/ SYM
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /WPT/ NDATE,NTIME
      COMMON /WP2/ TT
      COMMON /NZPATH/ ZFST,ZLST,RFST,RLST,XMAX,YMAX
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      XMIN=-ZL/20.
      XMAX=ZL-XMIN
      YMIN=-GR(IRMAX)/20.
      YMAX=GR(IRMAX)-YMIN
      XZERO=0.0
      YZERO=0.0
      RDUM=0.0
      ZDUM=0.0
      LRZF=.FALSE.
      WRITE(ND9,10000) XMIN,XMAX,YMIN,YMAX
      WRITE(ND9,10020) Z(1),R(1)
      WRITE(ND9,10020) Z(2),R(2)
      I=2
      WRITE(ND9,10010)
      ZFST=Z(2)
      RFST=R(2)
      NTN=0
   10 IF(I.GT.NZP) GO TO 60
      IF(LRZF) THEN
        RDUM=R(I)
        ZDUM=Z(I)
        LRZF=.FALSE.
      ELSE
        IF(I.EQ.NVC) THEN
          WRITE(ND9,10020) Z(I),R(I)
          IF(NZP-2.NE.NVC) WRITE(ND9,10030)
          I=I+3
          LRZF=.TRUE.
          NTN=0
          GO TO 10
        ELSE
          IF((ABS(RDUM-R(I)).GT.CEPS).OR.(ABS(ZDUM-Z(I)).GT.CEPS))
     +    GO TO 20
          LRZF=.TRUE.
        END IF
      END IF
   20 IF(R(I).GT.-1.) THEN
C-STRAIGHT LINE----
        WRITE(ND9,10020) Z(I),R(I)
        NTN=NTN+1
        IF((NTN.GE.1021).AND.(I.NE.NZP)) THEN
          WRITE(ND9,10030)
          WRITE(ND9,10020) Z(I),R(I)
          NTN=1
        END IF
        IF((LRZF).AND.(I.NE.NZP)) THEN
          WRITE(ND9,10030)
        END IF
        I=I+1
      ELSE
        IKIND=R(I)
        IF(IKIND.GE.-2) THEN
C-CIRCLE----
          RAD=ABS(Z(I))
          LLONG=IKIND.EQ.-2
          LCAV=Z(I).GT.0.0
          LVEX=Z(I).LT.0.0
          AXS1(1)=R(I-1)
          AXS2(1)=R(I+1)
          AXS1(2)=Z(I-1)
          AXS2(2)=Z(I+1)
          AXS12(1)=AXS2(1)-AXS1(1)
          AXS12(2)=AXS2(2)-AXS1(2)
          DST12=SQRT(AXS12(1)**2+AXS12(2)**2)
          IF(DST12.LE.2.*RAD) GO TO 30
          RAD=0.50*DST12
   30     ARG=DST12/(2.*RAD)
          PHIX=2.0*ASIN(ARG)
          PHX=ACOS(ARG)
          IF(LLONG) PHX=-PHX
          IF(LCAV) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)+SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LCAV) CENT(2)=AXS1(2)+(-SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2)
     +    ) *(RAD/DST12)
          IF(LVEX) CENT(1)=AXS1(1)+(COS(PHX)*AXS12(1)-SIN(PHX)*AXS12(2))
     +    *(RAD/DST12)
          IF(LVEX) CENT(2)=AXS1(2)+(SIN(PHX)*AXS12(1)+COS(PHX)*AXS12(2))
     +    *(RAD/DST12)
C
          A1CENT(1)=AXS1(1)-CENT(1)
          A1CENT(2)=AXS1(2)-CENT(2)
          IF(LLONG) PHIX=2.*PI-PHIX
          DRR=GDR(IRR(I))
          DDD=AMIN1(DDZ,DRR)
          NPOLY=(5.0*RAD/DDD)*(PHIX/PI)+1
          STEP=1.0/FLOAT(NPOLY)
          T=0.0
          PL=R(I-1)
          QL=Z(I-1)
          DO 40 J=2,NPOLY
            T=T+STEP
            IF(LCAV) P=CENT(1)+COS(PHIX*T)*A1CENT(1)+SIN(PHIX*T)*A1CENT
     +      (2)
            IF(LCAV) Q=CENT(2)-SIN(PHIX*T)*A1CENT(1)+COS(PHIX*T)*A1CENT
     +      (2)
            IF(LVEX) P=CENT(1)+COS(PHIX*T)*A1CENT(1)-SIN(PHIX*T)*A1CENT
     +      (2)
            IF(LVEX) Q=CENT(2)+SIN(PHIX*T)*A1CENT(1)+COS(PHIX*T)*A1CENT
     +      (2)
            IF((ABS(P-PL).GT.CEPS).OR.(ABS(Q-QL).GT.CEPS)) THEN
              WRITE(ND9,10020) Q,P
              NTN=NTN+1
              IF(NTN.GE.1021) THEN
                WRITE(ND9,10030)
                WRITE(ND9,10020) Q,P
                NTN=1
              ELSE
              END IF
            ELSE
            END IF
            PL=P
            QL=Q
   40     CONTINUE
          I=I+1
        ELSE
C-ELLIPSE
          R1=R(I-1)
          RM=R(I+1)
          R2=R(I+2)
          Z1=Z(I-1)
          ZM=Z(I+1)
          Z2=Z(I+2)
          DRR=GDR(IRD(RM))
          ELLP=(R1-RM)**2*(Z2-ZM)**2-(R2-RM)**2*(Z1-ZM)**2
          A=((Z2-ZM)**2-(Z1-ZM)**2)/ELLP
          B=((R1-RM)**2-(R2-RM)**2)/ELLP
          A=1.0/SQRT(A)
          B=1.0/SQRT(B)
C
          IF(R1.GT.RM) PHX=ACOS((Z1-ZM)/B)
          IF(R2.GT.RM) PHY=ACOS((Z2-ZM)/B)
          IF(R1.LT.RM) PHX=-ACOS((Z1-ZM)/B)
          IF(R2.LT.RM) PHY=-ACOS((Z2-ZM)/B)
          IF((R1.EQ.RM).AND.(Z1.GE.ZM)) PHX=0.
          IF((R1.EQ.RM).AND.(Z1.LT.ZM)) PHX=PI
          IF((R2.EQ.RM).AND.(Z2.GE.ZM)) PHY=0.
          IF((R2.EQ.RM).AND.(Z2.LT.ZM)) PHY=PI
          DPH1=PHY-PHX
          IF(DPH1.GE.2.0*PI) DPH1=DPH1-2.*PI
          IF(DPH1.LE.-2.0*PI) DPH1=DPH1+2.*PI
          IF(DPH1.GE.0.) DPH2=DPH1-2.*PI
          IF(DPH1.LT.0.) DPH2=DPH1+2.*PI
          IF((IKIND.EQ.-3).AND.(ABS(DPH1).LT.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-3).AND.(ABS(DPH2).LT.ABS(DPH1))) DELPH=DPH2
          IF((IKIND.EQ.-4).AND.(ABS(DPH1).GE.ABS(DPH2))) DELPH=DPH1
          IF((IKIND.EQ.-4).AND.(ABS(DPH2).GE.ABS(DPH1))) DELPH=DPH2
C
          NPOLY=5.*AMAX1(A/DRR,B/DDZ)+1
          STEP=1./FLOAT(NPOLY)
          T=0.0
          PL=R(I-1)
          QL=Z(I-1)
          DO 50 J=2,NPOLY
            T=T+STEP
            P=RM+A*SIN(PHX+T*DELPH)
            Q=ZM+B*COS(PHX+T*DELPH)
            IF((ABS(P-PL).GT.CEPS).OR.(ABS(Q-QL).GT.CEPS)) THEN
              WRITE(ND9,10020) Q,P
              NTN=NTN+1
              IF(NTN.GE.1021) THEN
                WRITE(ND9,10030)
                WRITE(ND9,10020) Q,P
                NTN=1
              ELSE
              END IF
            ELSE
            END IF
            PL=P
            QL=Q
   50     CONTINUE
          I=I+2
        END IF
      END IF
      GO TO 10
   60 CONTINUE
      WRITE(ND9,10030)
      WRITE(ND9,10020) Z(NVC),R(NVC)
      WRITE(ND9,10020) Z(NVC+1),R(NVC+1)
      WRITE(ND9,10010)
      WRITE(ND9,10020) Z(NVC+1),R(NVC+1)
      WRITE(ND9,10020) XZERO,YZERO
      WRITE(ND9,10040)
      ZLST=Z(NVC)
      RLST=R(NVC)
      RETURN
10000 FORMAT(1X,' SET LIMITS X ',2G13.5,' Y ',2G13.5)
10010 FORMAT(1X,' JOIN 1 DOTS')
10020 FORMAT(1X,2G13.5)
10030 FORMAT(1X,' JOIN 1 SOLID')
10040 FORMAT(1X,' JOIN 1 DOTDASH')
      END
C
      SUBROUTINE CHGPLT(N,YMIN)
C-----------------------------------------------------------------------
C PLOT THE CHARGE DISTRIBUTION ON TOP OF E-FILED LINE PLOTS
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /ELINE/ NZF,NZL
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      OFFY=GR(JOFF)
      N3=N-1-(N-1)/M*M
      DO 30 J=5,0,-1
        ILINE=0
        FACT=FLOAT(J)/5.
        DO 20 IZ=1,NZL+1
          N2=(N-1)/M-(IZ-1)+1
          IF(N2.LT.1.OR.N2.GT.IBT) GO TO 20
          IF(IZ.NE.1.OR.J.LT.5) GO TO 10
          ZZ=FLOAT(IZ-1)*DDZ
          QK=OFFY
          WRITE(ND9,10000) ZZ,QK
          ILINE=ILINE+1
   10     ZZ=FLOAT(IZ-1)*DDZ
          QK=QKST(N2)*YMIN*0.75/QKMAX*FACT+OFFY
          WRITE(ND9,10000) ZZ,QK
          ILINE=ILINE+1
          IF(ILINE.EQ.1021) THEN
            WRITE(ND9,10010)
            WRITE(ND9,10000) ZZ,QK
            ILINE=1
          END IF
          IF(IZ.NE.NZ+1.OR.J.LT.5) GO TO 20
          ZZ=FLOAT(IZ-1)*DDZ
          QK=OFFY
          WRITE(ND9,10000) ZZ,QK
          ILINE=ILINE+1
          IF(ILINE.EQ.1021) THEN
            WRITE(ND9,10010)
            WRITE(ND9,10000) ZZ,QK
            ILINE=1
          END IF
   20   CONTINUE
        IF(ILINE.GT.1) WRITE(ND9,10010)
   30 CONTINUE
      RETURN
10000 FORMAT(1X,2G13.5)
10010 FORMAT(1X,' JOIN 1 SOLID')
      END
C
      SUBROUTINE UPLOT(ZTXT)
C-----------------------------------------------------------------------
C TYPE X-Y COORDINATES OF CAVITY SHAPE USED ON FT09 FOR TOP DRAWER
C-----------------------------------------------------------------------
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL SYM,LNAPOLY,LPATH,LNONAP
      CHARACTER*13 REDCHR
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /WPT/ NDATE,NTIME
      COMMON /RED/ REDCHR
      COMMON /WP2/ TT
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /SYMT/ SYM
      COMMON /GEO/ R(NBP),Z(NBP),ZL
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /WPL2/LNAPOLY,LNONAP
      COMMON/ PLNZ/ LPATH
      COMMON /PINT/ KMIN
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      DATA SCM /1000./
      IMAX=1
   10 IF(DDR(IMAX+1).EQ.0.) GO TO 20
      IMAX=IMAX+1
      GO TO 10
   20 IF(IMAX.GT.6) IMAX=6
      WRITE(ND9,10000) REDCHR,(ZTXT(J),J=1,3),(ZTXT(L),L=4,9)
      WRITE(ND9,10010) NDATE,NTIME
      WRITE(ND9,10020) DDZ*SCM,DDR(1)*SCM
      IF(IMAX.GE.2) THEN
        DO 30 I=2,IMAX
   30   WRITE(ND9,10030) DDR(I)*SCM
      END IF
      XMIN=-ZL/20.
      XMAX=ZL-XMIN
      XARS=XMAX-XMIN
      XARE=XMAX-XMIN/10.
      YMIN=-GR(IRMAX)/20.
      YMAX=GR(IRMAX)-YMIN
      XZERO=0.0
      YZERO=0.0
      WRITE(ND9,10040) XMIN,XMAX,YMIN,YMAX
C DRAW ARROWS TO POINT OUT WHERE DDR CHANGES
      IF(IMAX.GE.2) THEN
        WRITE(ND9,10055)
        DO 40 I=1,IMAX-1
          WRITE(ND9,10060) XARS,RMARK(I),XARE,RMARK(I)
   40   CONTINUE
      END IF
C START PLOTTING CAVITY SHAPE
      WRITE(ND9,10050) Z(1),R(1)
      WRITE(ND9,10050) Z(2),R(2)
      WRITE(ND9,10070)
C-----------------------------------------------------------------------
C TYPE MESH BOUNDARY
C-----------------------------------------------------------------------
      LFR=1
      LBK=NZ
      MZF=1
      MZL=NZ
      CALL HIT1(1)
      WRITE(ND9,10080)
      WRITE(ND9,10050) Z(NVC),R(NVC)
      WRITE(ND9,10050) Z(NVC+1),R(NVC+1)
      WRITE(ND9,10070)
      WRITE(ND9,10050) Z(NVC+1),R(NVC+1)
      WRITE(ND9,10050) XZERO,YZERO
      WRITE(ND9,10090)
C
      IF (LPATH) THEN
        ZFST=Z(2)
        RFST=R(2)
        ZLST=Z(NVC)
        RLST=R(NVC)
        YDEV=YMAX/200.
        XDEV=XMAX/200.
        IF (LNAPOLY) THEN
          ZCF=FLOAT(IZCF-1)*DDZ
          ZCT=FLOAT(IZCT-1)*DDZ
          RWAK=GR(IRWAK)
          WRITE(ND9,10050) ZFST,RFST-YDEV
          WRITE(ND9,10050) ZCF-XDEV,RFST-YDEV
          WRITE(ND9,10050) ZCF-XDEV,RWAK-YDEV
          WRITE(ND9,10050) ZCT+XDEV,RWAK-YDEV
          WRITE(ND9,10050) ZCT+XDEV,RLST-YDEV
          WRITE(ND9,10050) ZLST,RLST-YDEV
          WRITE(ND9,10100)
        ELSE
          RLINE=GR(KMIN)
          WRITE(ND9,10050) ZFST,RLINE-YDEV
          WRITE(ND9,10050) ZLST,RLINE-YDEV
          WRITE(ND9,10100)
        END IF
      END IF
      RETURN
10000 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/
     +1X,' SET WINDOW X FROM 2 12 Y FROM 1 TO 8'/
     +1X,' SET TICKS SIZE 0.05'/ 1X,' TITLE 7.0 9.0 CENTER SIZE 2.2  ',
     +1H', 'Cavity Shape Used',1H'/ 1X,' TITLE 2.0 8.6 SIZE 1.5  ',
     +1H','A B C I 9.4 : ',3A8,1H'/ 1X,' MORE ',1H',6A8,1H'/ 1X,
     +' TITLE BOTTOM ',1H','Z-axis (m)',1H'/
     +1X,' TITLE LEFT   ',1H','R-axis (m)',
     +1H')
10010 FORMAT(1X,' TITLE 10.9 8.9 SIZE 1.25 ',1H',A8,2X,A8,1H')
10020 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H', 'DDZ=',F6.3,' mm,',
     +2X,'DDR=',F6.3,' mm',1H')
10030 FORMAT(1X,' MORE',1X,1H',', ',F6.3,' mm',1H')
10040 FORMAT(1X,' SET LIMITS X ',2G13.5,' Y ',2G13.5)
10050 FORMAT(1X,2G13.5)
10055 FORMAT(1X,' SET ARROW FLARE 0.3')
10060 FORMAT(1X,' ARROW FROM ',2G11.4,' DATA TO ',2G11.4,' DATA')
10070 FORMAT(1X,' JOIN 1 DOTS')
10080 FORMAT(1X,' JOIN 1 SOLID'/ 1X,' PLOT')
10090 FORMAT(1X,' JOIN 1 DOTDASH')
10100 FORMAT(1X,' JOIN 1 DASH')
      END
C
      SUBROUTINE WPLOT(ZTXT,EFI,NW,KTYPE,QXYMM)
C-----------------------------------------------------------------------
C PRINT WAKE POTENTIALS ON FT09
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LNAPOLY,LCHIN,LCBACK,LCFRON,LNONAP
      CHARACTER*13 REDCHR
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      DIMENSION EFI(NW),QXYMM(4)
      CHARACTER*8 CURV(5,2)
      CHARACTER*14 ALEG(5)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WPL/ MROT
      COMMON /WPL2/ LNAPOLY,LNONAP
      COMMON /WPT/ NDATE,NTIME
      COMMON /RED/ REDCHR
      COMMON /WP2/ TT
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /WZERO/ XZERO,YZERO,YSTEP
      SAVE /WZERO/
      DATA CURV/'DOTS    ','DOTDASH ','DASH    ','SOLID   ','DASH    ',
     +'DOTS    ','DOTDASH ','SOLID   ','DASH    ','SOLID   '/
      DATA ALEG/'Charge Density','Azimuthal     ', 'Transverse    ',
     +'Longitudinal  ','With Log. Term'/
      DATA SCM /1000./
      IF(KTYPE.GT.0) GO TO 40
      IMAX=1
   10 IF(DDR(IMAX+1).EQ.0.) GO TO 20
      IMAX=IMAX+1
      GO TO 10
   20 IF(IMAX.GT.4) IMAX=4
      WRITE(ND9,10000) REDCHR,(ZTXT(J),J=1,3),(ZTXT(L),L=4,9)
      WRITE(ND9,10010) NDATE,NTIME
      WRITE(ND9,10020) TT
      WRITE(ND9,10030) MROT,SIG*100.,DDZ*SCM,DDR(1)*SCM
      IF(IMAX.GE.2) THEN
        DO 30 I=2,IMAX
   30   WRITE(ND9,10040) DDR(I)*SCM
      END IF
      XMIN=QXYMM(1)
      XMAX=QXYMM(2)
      YMIN=QXYMM(3)
      YMAX=QXYMM(4)
      XZERO=XMAX/20.
      YSTEP=YMAX/10.
      YZERO=YMAX-2.*YSTEP
      WRITE(ND9,10050) XMIN,XMAX,YMIN,YMAX
      RETURN
   40 CONTINUE
      NTN=0
      DO 60 I=1,NW
        S=FLOAT(I-1)*DDZ
        IF(KTYPE.GT.1) GO TO 50
        IF((I.EQ.1).OR.(I.EQ.NW)) GO TO 50
        IF((EFI(I).EQ.0.).AND.(EFI(I-1).EQ.0.)) GO TO 60
   50   WRITE(ND9,10060) S,EFI(I)
        NTN=NTN+1
        IF(NTN.GE.1021) THEN
          WRITE(ND9,10070) CURV(KTYPE,MROT+1)
          WRITE(ND9,10060) S,EFI(I)
          NTN=1
        ELSE
        END IF
   60 CONTINUE
      WRITE(ND9,10070) CURV(KTYPE,MROT+1)
      WRITE(ND9,10060) XZERO,YZERO
      WRITE(ND9,10060) 3.*XZERO,YZERO
      WRITE(ND9,10080) CURV(KTYPE,MROT+1)
      WRITE(ND9,10090) 3.5*XZERO,YZERO,ALEG(KTYPE)
      YZERO=YZERO-YSTEP
      RETURN
10000 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/ 1X,
     +' SET WINDOW X FROM 2 12 Y FROM 1.75 TO 8'/ 1X,
     +' SET TICKS SIZE 0.05'/ 1X,' TITLE 7.0 9.0 CENTER SIZE 2.2  ', 1H'
     +,'Wake Potentials',1H'/ 1X,' TITLE 2.0 8.6 SIZE 1.50 ',1H',
C    41X,' SET SYMBOL SIZE 0.0'/
     +'A B C I 9.4 : ',3A8,1H'/ 1X,' MORE ',1H',6A8,1H'/ 1X,
     +' TITLE BOTTOM ',1H','Distance from Bunch Head S (m)',1H'/
     +1X,' TITLE LEFT   ',1H','Scaled Wake Potential W (S)',1H')
10010 FORMAT(1X,' TITLE 10.9 8.9 SIZE 1.25 ',1H',A8,2X,A8,1H')
10020 FORMAT(1X,' TITLE 10.1 9.1 SIZE 1.25 ',1H','Cpu Time Used:', 1PE10
     +.3,'(s)',1H')
10030 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H', 'MROT=',I2,',',2X,
     + 'SIG=',F6.3,' cm,',1H'/
     +1X,' MORE',1X, 1H',2X,'DDZ=',F6.3,' mm,',2X,'DDR=',
     +F6.3,' mm',1H')
10040 FORMAT(1X,' MORE',1X,1H',', ',F6.3,' mm',1H')
10050 FORMAT(1X,' SET LIMITS X ',2G13.5,' Y ',2G13.5)
10060 FORMAT(1X,2G13.5)
10070 FORMAT(1X,' JOIN ',A8)
10080 FORMAT(1X,' JOIN 1 ',A8)
10090 FORMAT(1X,' TITLE ',2G13.5,' DATA SIZE 1.5 ',1H',A14,1H')
      END
C
      SUBROUTINE WPADD(QMIN,QMAX,PKY,MROT,KTYPE)
C-----------------------------------------------------------------------
C PRINT INFORMATION ABOUT WAKE POTENTIALS AND LOSS FACTORS ON FT09
C-----------------------------------------------------------------------
      CHARACTER*18 ALEG(5)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /AZERO/ XZERO,YZERO,YSTEP
      SAVE /AZERO/
      DATA ALEG/'Charge Density    ','Azimuthal Wake    ',
     +'Transverse Wake   ','Longitudinal Wake ','With Log. Term    '/
      IF(KTYPE.GT.0) GO TO 10
      XZERO=2.0
      YZERO=0.7
      YSTEP=0.23
      RETURN
   10 CONTINUE
      MM=1
      IF((KTYPE.EQ.4).AND.(MROT.EQ.0)) MM=0
      IF((KTYPE.EQ.4).AND.(MROT.EQ.1)) MM=2
      IF(KTYPE.EQ.5) MM=0
      WRITE(ND9,10000) XZERO,YZERO,ALEG(KTYPE)
      IF(MM-1) 20,30,40
   20   WRITE(ND9,10010) XZERO+2.0,YZERO,QMIN,QMAX
        WRITE(ND9,10020) XZERO+6.7,YZERO,PKY
        GO TO 50
   30   WRITE(ND9,10012) XZERO+2.0,YZERO,QMIN,QMAX
        WRITE(ND9,10022) XZERO+6.7,YZERO,PKY
        GO TO 50
   40   WRITE(ND9,10014) XZERO+2.0,YZERO,QMIN,QMAX,MM
        WRITE(ND9,10016)
        WRITE(ND9,10024) XZERO+6.7,YZERO,PKY,MM
        WRITE(ND9,10026)
   50 YZERO=YZERO-YSTEP
      RETURN
10000 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ',1H',A18,1H')
10010 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ', 1H','Min/Max= ',1PE10.3,
     +'/',1PE10.3,' V/pC,  ',3X,1H')
10012 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ', 1H','Min/Max= ',1PE10.3,
     +'/',1PE10.3,' V/pC/m,',3X,1H')
10014 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ', 1H','Min/Max= ',1PE10.3,
     +'/',1PE10.3,' V/pC/m2',I1,'3,',1H')
10016 FORMAT(1X,' CASE  ',1H',30X,'       X X ',1H')
10020 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ',1H',
     +'Loss Factor= ', 1PE10.3,' V/pC   ',2X,1H')
10022 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ',1H',
     +'Loss Factor= ', 1PE10.3,' V/pC/m ',2X,1H')
10024 FORMAT(1X,' TITLE ',2F5.2,' SIZE 1.40 ',1H',
     +'Loss Factor= ', 1PE10.3,' V/pC/m2',I1,'3',1H')
10026 FORMAT(1X,' CASE  ',1H',23X,'       X X',1H')
      END
C
      SUBROUTINE PRINW(NW,LNTR,LPOT)
C-----------------------------------------------------------------------
C PRINT WAKE POTENTIALS ON LINEPRINTER
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      DIMENSION QDUMMY(10)
      LOGICAL LNTR,LPOT
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WPL/ MROT
      IF(.NOT.LNTR) THEN
        IF(MROT.EQ.0) GO TO 50
        WRITE(ND6,10000)
        DO 20 I=1,NW/10+1
          DO 10 L=1,10
            LI=L+(I-1)*10-1
   10     QDUMMY(L)=FLOAT(LI)*DDZ
          WRITE(ND6,10070) (QDUMMY(L),L=1,10)
          WRITE(ND6,10010) (WY(L+(I-1)*10),L=1,10)
   20   CONTINUE
        WRITE(ND6,10020)
        DO 40 I=1,NW/10+1
          DO 30 L=1,10
            LI=L+(I-1)*10-1
   30     QDUMMY(L)=FLOAT(LI)*DDZ
          WRITE(ND6,10070) (QDUMMY(L),L=1,10)
          WRITE(ND6,10030) (WT(L+(I-1)*10),L=1,10)
   40   CONTINUE
   50   IF(MROT.EQ.1) WRITE(ND6,10040)
        IF(MROT.EQ.0) WRITE(ND6,10050)
        DO 70 I=1,NW/10+1
          DO 60 L=1,10
            LI=L+(I-1)*10-1
   60     QDUMMY(L)=FLOAT(LI)*DDZ
          WRITE(ND6,10070) (QDUMMY(L),L=1,10)
          WRITE(ND6,10060) (WZ(L+(I-1)*10),L=1,10)
   70   CONTINUE
        IF(LPOT) THEN
          DO 90 I=1,NW/10+1
            DO 80 L=1,10
              LI=L+(I-1)*10-1
   80       QDUMMY(L)=FLOAT(LI)*DDZ
            WRITE(ND6,10070) (QDUMMY(L),L=1,10)
            WRITE(ND6,10080) (WP(L+(I-1)*10),L=1,10)
   90     CONTINUE
        ELSE
        END IF
      ELSE
        WRITE(ND6,10020)
        DO 110 I=1,NW/10+1
          DO 100 L=1,10
            LI=L+(I-1)*10-1
  100     QDUMMY(L)=FLOAT(LI)*DDZ
          WRITE(ND6,10070) (QDUMMY(L),L=1,10)
          WRITE(ND6,10030) (WT(L+(I-1)*10),L=1,10)
  110   CONTINUE
        IF(MROT.EQ.1) WRITE(ND6,10040)
        DO 130 I=1,NW/10+1
          DO 120 L=1,10
            LI=L+(I-1)*10-1
  120     QDUMMY(L)=FLOAT(LI)*DDZ
          WRITE(ND6,10070) (QDUMMY(L),L=1,10)
          WRITE(ND6,10060) (WZ(L+(I-1)*10),L=1,10)
  130   CONTINUE
      END IF
      RETURN
10000 FORMAT(1H1/1X,'AZIMUTHAL WAKE FIELDS --->'/)
10010 FORMAT(1X,'AWAKE/V  =',1X,10(1X,1PE11.4)/)
10020 FORMAT(1H1/1X,'TRANSVERSE WAKE FIELDS --->'/)
10030 FORMAT(1X,'TWAKE/V  =',1X,10(1X,1PE11.4)/)
10040 FORMAT(1H0/1X,'LONGITUDINAL WAKE FIELDS --->'/)
10050 FORMAT(1H1/1X,'LONGITUDINAL WAKE FIELDS --->'/)
10060 FORMAT(1X,'LWAKE/V  =',1X,10(1X,1PE11.4)/)
10070 FORMAT(1X,'S/M      =',1X,10(1X,1PE11.4))
10080 FORMAT(1X,'PWAKE/V  =',1X,10(1X,1PE11.4)/)
      END
C
      SUBROUTINE WPRIN(ZTXT,EFI,NW,KTYPE,FACT)
C-----------------------------------------------------------------------
C PRINT WAKE POTENTIALS ON FT08
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      DIMENSION EFI(NW)
      CHARACTER*40 ALEG(5)
      CHARACTER*18 WKPT(5)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WPL/ MROT
      COMMON /WPT/ NDATE,NTIME
      DATA ALEG/
     1'LONGITUDINAL WAKE POTENTIAL             ',
     1'AZIMUTHAL WAKE POTENTIAL                ',
     1'TRANSVERSE WAKE POTENTIAL               ',
     1'LONGITUDINAL WAKE POTENTIAL             ',
     1'LONGITUDINAL WAKE POTENTIAL /W LOG. TERM'/
      DATA WKPT/'    WZ  (V/pC)    ','    WA (V/pC/m)   ',
     +'    WT (V/pC/m)   ','   WZ (V/pC/m**2) ', '    WZ  (V/pC)    '/
      WRITE(ND8,10000) ALEG(KTYPE), (ZTXT(II),II=1,9)
      WRITE(ND8,10010)
      WRITE(ND8,10020) NDATE,NTIME,MROT,NW
      WRITE(ND8,10030) WKPT(KTYPE)
      DO 10 I=1,NW
        S=FLOAT(I-1)*DDZ
        WRITE(ND8,10040) S,EFI(I)*FACT
   10 CONTINUE
      WRITE(ND8,10050)
      RETURN
10000 FORMAT(1X,1X,'TITLE: ',A40/1X,1X,9A8)
10010 FORMAT(1X,2X,'DATE:',5X,'TIME:',6X,'MROT:', 5X,'NO. OF POINTS:')
10020 FORMAT(1X,1X,A8,2X,A8,I7,7X,I8)
10030 FORMAT(1X,'     S  (m)',A18)
10040 FORMAT(1X,2G13.5)
10050 FORMAT(1X)
      END
C
      SUBROUTINE ZPLOT(ZTXT,W,NWW0,ITYPE,CUTOFF,WFT)
C-----------------------------------------------------------------------
C PRINT REAL AND IMAGINARY PARTS OF IMPEDANCE ON FT09
C THE ORIGINAL VERSION OF THIS ROUTINE WAS WRITTEN BY BILL NG
C AT FERMI NATIONAL LABORATORY IN 1992.
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NWN=500000)
      PARAMETER (NF=100)
      PARAMETER (NPF=624)
      LOGICAL LNAPOLY,LCHIN,LCBACK,LCFRON,LSVF,LWNDW,LINTZ,LCRBW
      LOGICAL LNONAP
      CHARACTER*13 REDCHR
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      DIMENSION W(NWN),XZ(NPF),BGAU(NPF)
      DIMENSION WD(2*NWN),FR(NPF),FI(NPF),FS(NPF)
      CHARACTER*40 ALEG(19)
      CHARACTER*34 TITLE(19)
      CHARACTER*34 CASET(19)
      CHARACTER*17 TITLB(19)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WPL/ MROT
      COMMON /WPL2/ LNAPOLY,LNONAP
      COMMON /WPT/ NDATE,NTIME
      COMMON /RED/ REDCHR
      COMMON /WP2/ TT
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /FFTS/ LSVF,LWNDW
      COMMON /WFUN/ NWFUN,ALPHA
      COMMON /EXPF/ EXPFAC
      COMMON /INTR/ LINTZ
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      DATA ZERO/0.0/
      DATA SCM /1000./
      DATA ALEG/
     1'    Real Part of Azimuthal Impedance    ',
     1' Imaginary Part of Azimuthal Impedance  ',
     1'   Real Part of Transverse Impedance    ',
     1' Imaginary Part of Transverse Impedance ',
     1'  Real Part of Longitudinal Impedance   ',
     1'Imaginary Part of Longitudinal Impedance',
     1'  Real Part of Longitudinal Impedance   ',
     1'Imaginary Part of Longitudinal Impedance',
     1'    Real Part of Long.+Log Impedance    ',
     1' Imaginary Part of Long.+Log Impedance  ',
     1'   Frequency Spectrum of Loss Factor    ',
     1' Loss Factor Spectrum Integrated upto F ',
     1'   Spectrum of Long.+Log Loss Factor    ',
     1'Long.+Log Loss Factor Integrated upto F ',
     1' Azimuthal Impedance Integrated upto F  ',
     1' Transverse Impedance Integrated upto F ',
     1'Longitudinal Impedance Integrated upto F',
     1'Longitudinal Impedance Integrated upto F',
     1' Long.+Log Impedance Integrated upto F  '/
      DATA TITLE/
     1'               Real Z2A1  (kW/m)  ',
     1'               Imag Z2A1  (kW/m)  ',
     1'               Real Z2T1  (kW/m)  ',
     1'               Imag Z2T1  (kW/m)  ',
     1'               Real Z2L1  (kW)    ',
     1'               Imag Z2L1  (kW)    ',
     1'            Real Z2L1  (kW/m**2)  ',
     1'            Imag Z2L1  (kW/m**2)  ',
     1'             Real Z2L+LOG1  (kW)  ',
     1'             Imag Z2L+LOG1  (kW)  ',
     1'             dk2L1/df  (V/pC/GHz) ',
     1'                  k2L1(F)  (V/pC) ',
     1'             dk2L1/df  (V/pC/GHz) ',
     1'                  k2L1(F)  (V/pC) ',
     1'   I02F10101  Re Z2A1df (kWXGHz/m)',
     1'   I02F10101  Re Z2T1df (kWXGHz/m)',
     1'     I02F10101  Re Z2L1df (kWXGHz)',
     1'I02F10101  Re Z2L1df (kWXGHz/m**2)',
     1'     I02F10101  Re Z2L1df (kWXGHz)'/
      DATA CASET/
     1'                     V V    F     ',
     1'                     V V    F     ',
     1'                     V V    F     ',
     1'                     V V    F     ',
     1'                     V V    F     ',
     1'                     V V    F     ',
     1'                  V V    F        ',
     1'                  V V    F        ',
     1'                   V     V    F   ',
     1'                   V     V    F   ',
     1'               V V                ',
     1'                   V V            ',
     1'               V V                ',
     1'                   V V            ',
     1'   MZX ZXY Z      V V     FM      ',
     1'   MZX ZXY Z      V V     FM      ',
     1'     MZX ZXY Z      V V     FM    ',
     1'MZX ZXY Z      V V     FM         ',
     1'     MZX ZXY Z      V V     FM    '/
      DATA TITLB/
     1'Frequency f (GHz)','Frequency f (GHz)','Frequency f (GHz)',
     1'Frequency f (GHz)','Frequency f (GHz)','Frequency f (GHz)',
     1'Frequency f (GHz)','Frequency f (GHz)','Frequency f (GHz)',
     1'Frequency f (GHz)','Frequency f (GHz)','Frequency F (GHz)',
     1'Frequency f (GHz)','Frequency F (GHz)','Frequency F (GHz)',
     1'Frequency F (GHz)','Frequency F (GHz)','Frequency F (GHz)',
     1'Frequency F (GHz)'/
C --- COMPUTE IMPEDANCE FROM WAKE --------------------------------------
C          :THE WAKE FIELD STARTS FROM
C          :-ISIG*SIG AND IS TRUNCATED AT A LENGTH OF (NWW-1)*DDZ.
C          :THE BUNCH CENTER IS AT 0.  WE WILL ADD IN THE WAKE WHICH
C          :IS ZERO AT THE NEGATIVE SIDE SO THAT THE TOTAL NUMBER OF
C          :POINTS BECEOMES NW=NNW+NWW WHERE NNW=NWW-2*ISIG*SIG/DDZ-1.
C
C     NN=624 IS THE NUMBER OF HORIZONTAL POINTS IN THE PLOTS, WHICH
C          IS THE RESOLUTION OF THE PRINTER
C
C-----------------------------------------------------------------------
      NWW=NWW0
      IF(LCRBW) THEN
      NNW=NWW-1-2*INT(BSIG*SIG/DDZ+.5)-2*(NZ-NZSEP)
      ELSE
      NNW=NWW-1-2*INT(BSIG*SIG/DDZ+.5)
      END IF
      NW=NNW+NWW
      NN=NPF
      CL9=CLGHT*1.E-9
C "FREQ" IS UPPER LIMIT OF FREQUENCY SO THAT THE EXPONENTIAL CORRECT    ion
C  FACTOR DUE TO BUNCH LENGTH IS EXPFAC (DEFAULT=20).
C "RES" IS THE FREQUENCY LIMITED BY THE MESH SIZE, OR THE UPPER LIMI    t
C  DUE TO SAMPLING.
      FREQ=SQRT(2.*LOG(EXPFAC))*0.30/(SIG*2.*PI)
      RES=0.5E-9*CLGHT/DDZ
      IF(RES.LT.FREQ) FREQ=RES
      IF(CUTOFF.GT.0.001.AND.CUTOFF.LT.RES) FREQ=CUTOFF
C     IF(CUTOFF.GT.0.001.AND.CUTOFF.LT.FREQ) FREQ=CUTOFF
      WFCT=WFT
      IF(MROT.EQ.0) WFCT=1.0
C
      GAP=FREQ/FLOAT(NN)
      DT=DDZ/CL9
      DTT=DT*1.E-9
      T2=DT*FLOAT(NW)*0.5E0
      T3=2.E0/(FLOAT(NW)-1.E0)
      INW=NNW+1
      SGNY=1.
      IF(ITYPE.EQ.1) SGNY=-1.
      T4=-1.E0
C
      IF(LWNDW) THEN
        ALPCO=ALPHA
        IF(NWFUN.EQ.1) THEN
          ALPCO=ALPHA*PI
          BESA=BESSI0(ALPCO)
        ENDIF
        DO 40 I=1,NW
          IF(NWFUN-2) 10,20,30
   10     ARGS=AMAX1(1.E0-T4**2,0.E0)
          TAPI=ALPCO*SQRT(ARGS)
          WD(I)=BESSI0(TAPI)/BESA
          GO TO 40
   20     TPI=T4*PI
          WD(I)=0.40217E0+0.49703E0*COS(TPI)+0.09392E0*COS(2.E0*TPI)
     +    +0.00183E0*COS(3.E0*TPI)
          GO TO 40
   30     WD(I)=EXP(-0.5E0*(ALPCO*T4)**2)
   40   T4=T4+T3
      ELSE
        DO 50 I=1,NW
          WD(I)=1.0
   50   T4=T4+T3
      END IF
C
      DO 70 K=1,NN
        FR(K)=0.
        FI(K)=0.
        OMEGA=2.*PI*FLOAT(K-1)*GAP
        ANG1=-T2+(FLOAT(INW)-0.5)*DT
        DO 60 I=1,NWW
          ANG=OMEGA*ANG1
          FR(K)=FR(K)-W(I)*COS(ANG)*SGNY*WD(I+NNW)
          FI(K)=FI(K)+W(I)*SIN(ANG)*SGNY*WD(I+NNW)
   60   ANG1=ANG1+DT
        FR(K)=FR(K)*DTT*WFCT*1.E-3
        FI(K)=FI(K)*DTT*WFCT*1.E-3
   70 CONTINUE
      DO 80 K=1,NN
        EXPON=2.*PI*GAP*FLOAT(K)*SIG/CL9
        ZCORRECT=EXP(0.5*EXPON*EXPON)
        IF(ITYPE.GE.5) ZCORRECT=2./ZCORRECT
        IF(ITYPE.GE.5) BGAU(K)=(ZCORRECT/2.)**2
        FR(K)=FR(K)*ZCORRECT
        FI(K)=FI(K)*ZCORRECT
        IF(MROT.NE.1) GOTO 80
        IF(ITYPE.EQ.3) GOTO 80
        A=FI(K)
        FI(K)=-FR(K)
        FR(K)=A
   80 CONTINUE
C-INTEGRATION OF FREQUENCY SPECTRUM OF DEPOSITED ENEYGY
      IF(ITYPE.GE.5) THEN
        SUM=0.0
        DO 90 K=1,NN
          SUM=SUM+FR(K)
          FI(K)=SUM*GAP
   90   CONTINUE
      ELSE
        IF(LINTZ) THEN
          SUM=0.0
          DO 100 K=1,NN
            SUM=SUM+FR(K)
            FS(K)=SUM*GAP
  100     CONTINUE
        END IF
      END IF
C----------------------------------------------------------------------
C PREPARE TOPDRAWER INPUT FILE AND WRITE IMPEDANCE ON FT09
C----------------------------------------------------------------------
      FRMAX=0.0
      DO 110 I=1,NN
        IF(FRMAX.LT.FR(I)) FRMAX=FR(I)
  110 XZ(I)=(FLOAT(I)-1.)*GAP
      KEYMAX=1
      IF(LINTZ) KEYMAX=2
      IF(ITYPE.GE.5) KEYMAX=1
      KEY=0
  120 IF(ITYPE.EQ.1.AND.KEY.EQ.0) KTYPE=1
      IF(ITYPE.EQ.1.AND.KEY.EQ.1) KTYPE=2
      IF(ITYPE.EQ.1.AND.KEY.EQ.2) KTYPE=15
      IF(ITYPE.EQ.2.AND.KEY.EQ.0) KTYPE=3
      IF(ITYPE.EQ.2.AND.KEY.EQ.1) KTYPE=4
      IF(ITYPE.EQ.2.AND.KEY.EQ.2) KTYPE=16
      IF(ITYPE.EQ.3.AND.KEY.EQ.0.AND.MROT.EQ.0) KTYPE=5
      IF(ITYPE.EQ.3.AND.KEY.EQ.1.AND.MROT.EQ.0) KTYPE=6
      IF(ITYPE.EQ.3.AND.KEY.EQ.2.AND.MROT.EQ.0) KTYPE=17
      IF(ITYPE.EQ.3.AND.KEY.EQ.0.AND.MROT.EQ.1) KTYPE=7
      IF(ITYPE.EQ.3.AND.KEY.EQ.1.AND.MROT.EQ.1) KTYPE=8
      IF(ITYPE.EQ.3.AND.KEY.EQ.2.AND.MROT.EQ.1) KTYPE=18
      IF(ITYPE.EQ.4.AND.KEY.EQ.0) KTYPE=9
      IF(ITYPE.EQ.4.AND.KEY.EQ.1) KTYPE=10
      IF(ITYPE.EQ.4.AND.KEY.EQ.2) KTYPE=19
      IF(ITYPE.EQ.5.AND.KEY.EQ.0) KTYPE=11
      IF(ITYPE.EQ.5.AND.KEY.EQ.1) KTYPE=12
      IF(ITYPE.EQ.6.AND.KEY.EQ.0) KTYPE=13
      IF(ITYPE.EQ.6.AND.KEY.EQ.1) KTYPE=14
      IMAX=1
  130 IF(DDR(IMAX+1).EQ.0.) GO TO 140
      IMAX=IMAX+1
      GO TO 130
  140 IF(IMAX.GT.4) IMAX=4
      WRITE(ND9,10010) REDCHR,ALEG(KTYPE),
     +(ZTXT(J),J=1,3),(ZTXT(L),L=4,9),
     +TITLB(KTYPE),TITLE(KTYPE),CASET(KTYPE)
      WRITE(ND9,10020) NDATE,NTIME
      WRITE(ND9,10030) MROT,SIG*100.,DDZ*SCM,DDR(1)*SCM
      IF(IMAX.GE.2) THEN
        DO 150 I=2,IMAX
  150   WRITE(ND9,10040) DDR(I)*SCM
      END IF
      WRITE(ND9,10050) ZERO,FREQ
      ILINE=1
      DO 160 I=1,NN
        IF(KEY.EQ.0) WRITE(ND9,10060) XZ(I),FR(I)
        IF(KEY.EQ.1) WRITE(ND9,10060) XZ(I),FI(I)
        IF(KEY.EQ.2) WRITE(ND9,10060) XZ(I),FS(I)
        IF(ILINE.EQ.1021) THEN
          IF(KTYPE.NE.12.AND.KTYPE.NE.14.AND.KTYPE.LT.15)
     +    WRITE(ND9,10070)
          IF(KTYPE.EQ.12.OR.KTYPE.EQ.14.OR.KTYPE.GE.15)
     +    WRITE(ND9,10080)
          ILINE=0
        ENDIF
        ILINE=ILINE+1
  160 CONTINUE
      IF(KTYPE.NE.12.AND.KTYPE.NE.14.AND.KTYPE.LT.15) WRITE(ND9,10070)
      IF(KTYPE.EQ.12.OR.KTYPE.EQ.14.OR.KTYPE.GE.15) WRITE(ND9,10080)
      IF(.NOT.LSVF) GO TO 170
      IF(KEY.EQ.0) CALL ZPRIN(ZTXT,XZ,FR,NN,KTYPE)
      IF(KEY.EQ.1) CALL ZPRIN(ZTXT,XZ,FI,NN,KTYPE)
      IF(KEY.EQ.2) CALL ZPRIN(ZTXT,XZ,FS,NN,KTYPE)
  170 IF((KTYPE.NE.11).AND.(KTYPE.NE.13)) GO TO 190
      ILINE=1
      DO 180 I=1,NN
        WRITE(ND9,10060) XZ(I),BGAU(I)*FRMAX
        IF(ILINE.EQ.1021) THEN
          WRITE(ND9,10090)
          ILINE=0
        ENDIF
        ILINE=ILINE+1
  180 CONTINUE
      WRITE(ND9,10090)
      XZERO=FREQ/20.
      XOFF=XZERO*12.
      YZERO=FRMAX*0.9
      WRITE(ND9,10060) XOFF,YZERO
      WRITE(ND9,10060) 3.*XZERO+XOFF,YZERO
      WRITE(ND9,10090)
      WRITE(ND9,10100) 3.5*XZERO+XOFF,YZERO
  190 CONTINUE
      WRITE(ND9,10110)
      WRITE(ND9,10060) ZERO,ZERO
      WRITE(ND9,10060) FREQ,ZERO
      WRITE(ND9,10120)
      IF(KEY.EQ.KEYMAX) GO TO 771
      KEY=KEY+1
      GO TO 120
  771 RETURN
10000 FORMAT(1X,'*** RETURN ***  NN:',I7,' EXCEEDS MAXIMUM NPF:',I7,
     +' IN SUBROUTINE ZPLOT FOR FOURIER TRANSFORM.')
10010 FORMAT(1X,' NEW FRAME'/ 1X,' SET FONT DUPLEX'/ 1X,
     +' SET SIZE 13 BY 10 ',A13/ 1X,
     +' SET WINDOW X FROM 2 12 Y FROM 1.00 TO 8'/ 1X,
     +' SET TICKS SIZE 0.05'/ 1X,' TITLE 3.4 9.0 SIZE 2.2  ', 1H' ,' ',
     +1H'/ 1X,' MORE ',1H',A40,1H'/ 1X,' TITLE 2.0 8.6 SIZE 1.5  ',1H',
     +'A B C I 9.4 : ',3A8,1H'/ 1X,' MORE ',1H',6A8,1H'/ 1X,
     +' TITLE BOTTOM ',1H',A17,1H'/ 1X,' TITLE LEFT   ',1H',A34,1H'/ 1X,
     +' CASE         ',1H',A34,1H')
10020 FORMAT(1X,' TITLE 10.9 8.9 SIZE 1.25 ',1H',A8,2X,A8,1H')
10030 FORMAT(1X,' TITLE 2.00 8.3 SIZE 1.40 ',1H', 'MROT=',I2,',',2X,
     + 'SIG=',F6.3,' cm,',1H'/
     +1X,' MORE',1X, 1H',2X,'DDZ=',F6.3,' mm,',2X,'DDR=',F6.3,' mm',1H')
10040 FORMAT(1X,' MORE',1X,1H',', ',F6.3,' mm',1H')
10050 FORMAT(1X,' SET LIMIT X ',2G13.5)
10060 FORMAT(1X,2G13.5)
10070 FORMAT(1X,' JOIN ')
10080 FORMAT(1X,' JOIN 1 ')
10090 FORMAT(1X,' JOIN DOTS ')
10100 FORMAT(1X,' TITLE ',2G13.5,' DATA SIZE 1.5 ',1H', 'Bunch Spectrum'
     +,1H')
10110 FORMAT(1X,' SET ORDER X Y ')
10120 FORMAT(1X,' JOIN DOTDASH ')
      END
C
      SUBROUTINE ZPRIN(ZTXT,FREQ,AMP,NN,KTYPE)
C-----------------------------------------------------------------------
C PRINT FFT RESULTS ON FT08
C-----------------------------------------------------------------------
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      CHARACTER*8 ZTXT(10)
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      DIMENSION FREQ(NN),AMP(NN)
      CHARACTER*57 ALEG(19)
      CHARACTER*17 TITLE(19)
      CHARACTER*13 TITLB(19)
      COMMON /NDEV/ ND1,ND2,ND5,ND6,ND8,ND9
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /WPL/ MROT
      COMMON /WPT/ NDATE,NTIME
      DATA ALEG/
     1'REAL PART OF AZIMUTHAL IMPEDANCE                         ',
     1'IMAGINARY PART OF AZIMUTHAL IMPEDANCE                    ',
     1'REAL PART OF TRANSVERSE IMPEDANCE                        ',
     1'IMAGINARY PART OF TRANSVERSE IMPEDANCE                   ',
     1'REAL PART OF LONGITUDINAL IMPEDANCE                      ',
     1'IMAGINARY PART OF LONGITUDINAL IMPEDANCE                 ',
     1'REAL PART OF LONGITUDINAL IMPEDANCE                      ',
     1'IMAGINARY PART OF LONGITUDINAL IMPEDANCE                 ',
     1'REAL PART OF LONGITUDINAL+LOG IMPEDANCE                  ',
     1'IMAGINARY PART OF LONGITUDINAL+LOG IMPEDANCE             ',
     1'FREQUENCY SPECTRUM OF LONGITUDINAL LOSS FACTOR           ',
     1'LONGITUDINAL LOSS FACTOR INTEGRATED OVER FREQUENCY UPTO F',
     1'FREQUENCY SPECTRUM OF LONGITUDINAL+LOG LOSS FACTOR       ',
     1'LONG.+LOG LOSS FACTOR INTEGRATED OVER FREQUENCY UPTO F   ',
     1'AZIMUTHAL IMPEDANCE INTEGRATED OVER FREQUENCY UPTO F     ',
     1'TRANSVERSE IMPEDANCE INTEGRATED OVER FREQUENCY UPTO F    ',
     1'LONGITUDINAL IMPEDANCE INTEGRATED OVER FREQUENCY UPTO F  ',
     1'LONGITUDINAL IMPEDANCE INTEGRATED OVER FREQUENCY UPTO F  ',
     1'LONG.+LOG IMPEDANCE INTEGRATED OVER FREQUENCY UPTO F     '/
      DATA TITLE/
     1'Real Z (kOhm/m)  ','Imag Z (kOhm/m)  ','Real Z (kOhm/m)  ',
     1'Imag Z (kOhm/m)  ','Real Z (kOhm)    ','Imag Z (kOhm)    ',
     1'Real Z(kOhm/m**2)','Imag Z(kOhm/m**2)','Real Z+LOG (kOhm)',
     1'Imag Z+LOG (kOhm)','dk/df (V/pC/GHz) ','k(F) (V/pC)      ',
     1'dk/df (V/pC/GHz) ','k(F) (V/pC)      ','ReZ*f(kOhm-GHz/m)',
     1'ReZ*f(kOhm-GHz/m)','ReZ*f (kOhm-GHz) ','ReZ*f(kOhmGHz/m2)',
     1'ReZ*f (kOhm-GHz) '/
      DATA TITLB/     'Freq. f (GHz)','Freq. f (GHz)','Freq. f (GHz)',
     1'Freq. f (GHz)','Freq. f (GHz)','Freq. f (GHz)','Freq. f (GHz)',
     1'Freq. f (GHz)','Freq. f (GHz)','Freq. f (GHz)','Freq. f (GHz)',
     1'Freq. F (GHz)','Freq. f (GHz)','Freq. F (GHz)','Freq. F (GHz)',
     1'Freq. F (GHz)','Freq. F (GHz)','Freq. F (GHz)','Freq. F (GHz)'/
      WRITE(ND8,10000) ALEG(KTYPE), (ZTXT(II),II=1,9)
      WRITE(ND8,10010)
      WRITE(ND8,10020) NDATE,NTIME,MROT,NN
      WRITE(ND8,10030) TITLB(KTYPE),TITLE(KTYPE)
      DO 10 I=1,NN
        WRITE(ND8,10040) FREQ(I),AMP(I)
   10 CONTINUE
      WRITE(ND8,10050)
      RETURN
10000 FORMAT(1X,1X,'TITLE: ',A57/1X,1X,9A8)
10010 FORMAT(1X,2X,'DATE:',5X,'TIME:',6X,'MROT:', 5X,'NO. OF POINTS:')
10020 FORMAT(1X,1X,A8,2X,A8,I7,7X,I8)
10030 FORMAT(1X,1X,A13,1X,A17)
10040 FORMAT(1X,2G13.5)
10050 FORMAT(1X)
      END
C
      FUNCTION BESSI0(X)
C-----------------------------------------------------------------------
C MODIFIED BESSEL FUNCTION OF 0-TH ORDER
C-----------------------------------------------------------------------
      REAL BESSI0,X
      REAL AX
      DOUBLE PRECISION P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y
      SAVE P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,0.2635537D-1,
     *-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     +  (Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
      SUBROUTINE NAPL1(N,IA,NW)
C-----------------------------------------------------------------------
C NAPOLY INTEGRATION FOR MROT=0
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LCFRON,LCBACK,LCHIN,LCRBW
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /NAPCOM/ NRNAP1,NRNAP2
      SAVE /NAPCOM/
      COMMON /NAPPRM/ DZM,SFG(NRN),SLG(NRN)
      SAVE /NAPPRM/
C
      IPM=1
      IF(LCRBW) IPM=-1
C
      IF(NRNAP1.LT.NR) THEN
        NRNAP1=NR
        DZM=DDZ/FLOAT(M)
        DO 10 J=1,NR
          GDRM=GDR(J)/FLOAT(M)
          SFG(J)=SGF*GDRM
          SLG(J)=SGL*GDRM
   10   CONTINUE
      END IF
C
      NJDGZ=(N+1)/M*M
      IF(LCBACK) THEN
        NZF=(N+1)/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
      NZL=(N+1)/M+2
      IF(NZL.GE.NZ) NZL=NZ
      IF(N+1.NE.NJDGZ) THEN
        DO 40 I=NZF,NZL
          L=I
          IF(LCRBW) L=NZ-I+1
          LL=I
          IF(LCRBW) LL=NZSEP-I+1
          IS=IA-LL+1
          IF(IS.LE.0) GO TO 40
          IF(IS.GT.NW+2) GO TO 40
          FCT=1.0
          IF((I.EQ.IZCF).OR.(I.EQ.IZCT)) FCT=0.5
          KK=(I-MZF)*NR
          IF(I.LE.IZCF) THEN
            K=KOFF+KK
            WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
          END IF
          IF(I.GE.IZCT) THEN
            K=KOFL+KK
            WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
          END IF
          IF(IZCT.GT.IZCF) THEN
            IF((I.GE.IZCF).AND.(I.LE.IZCT)) THEN
              K=IRWAK+KK
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
            END IF
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCF) THEN
              DO 20 J=KSTF,KENF
                K=J+KK
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)-(EV(1,K)+EV(1,K+IPM*NR))*SFG(J)/2.
                ELSE
                  WZ(IS)=WZ(IS)-EV(1,K)*SFG(J)
                END IF
                IF(N.EQ.NJDGZ) THEN
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)-ZO*(HV(2,K)+7.*HV
                    WZ(IS)=WZ(IS)-ZO*(7.*HV(2,K)+HV
     +              (2,K-IPM*NR))*SFG(J) /8.
                  ELSE
                    WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)
                  END IF
                ELSE
                  WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)
                END IF
   20         CONTINUE
            END IF
          END IF
          IF(IRWAK.NE.KOFL) THEN
            IF(I.EQ.IZCT) THEN
              DO 30 J=KSTL,KENL
                K=J+KK
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)+(EV(1,K)+EV(1,K+IPM*NR))*SLG(J)/2.
                ELSE
                  WZ(IS)=WZ(IS)+EV(1,K)*SLG(J)
                END IF
                IF(N.EQ.NJDGZ) THEN
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)+ZO*(HV(2,K)+7.*HV
                    WZ(IS)=WZ(IS)+ZO*(7.*HV(2,K)+HV
     +              (2,K-IPM*NR))*SLG(J) /8.
                  ELSE
                    WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)
                  END IF
                ELSE
                  WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)
                END IF
   30         CONTINUE
            END IF
          END IF
   40   CONTINUE
      ELSE
        DO 70 I=NZF,NZL
          L=I
          IF(LCRBW) L=NZ-I+1
          LL=I
          IF(LCRBW) LL=NZSEP-I+1
          IS=IA-LL+1
          IF(IS.LE.0) GO TO 70
          IF(IS.GT.NW+2) GO TO 70
          FCT=1.0
          IF((I.EQ.IZCF).OR.(I.EQ.IZCT)) FCT=0.5
          KK=(I-MZF)*NR
          IF(I.LE.IZCF) THEN
            K=KOFF+KK
            IF(L.NE.NZ) THEN
              WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,K+IPM*NR))*DZM/2.*FCT
            ELSE
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
            END IF
          END IF
          IF(I.GE.IZCT) THEN
            K=KOFL+KK
            IF(L.NE.NZ) THEN
              WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,K+IPM*NR))*DZM/2.*FCT
            ELSE
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
            END IF
          END IF
          IF(IZCT.GT.IZCF) THEN
            IF((I.GE.IZCF).AND.(I.LE.IZCT)) THEN
              K=IRWAK+KK
              IF(L.NE.NZ) THEN
                WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,K+IPM*NR))*DZM/2.*FCT
              ELSE
                WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
              END IF
            END IF
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCF) THEN
              DO 50 J=KSTF,KENF
                K=J+KK
                KP=K+IPM*NR
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)-(EV(1,K)+2.*EV
     +            (1,KP)+ EV(1,KP+IPM*NR))*SFG(J)/4.
                ELSE
                  WZ(IS)=WZ(IS)-EV(1,K)*SFG(J)
                END IF
                IF(L.NE.NZ) THEN
C                 WZ(IS)=WZ(IS)-ZO*(HV(2,K)+7.*HV(2,K+NR))*SFG(J)/8.
                  WZ(IS)=WZ(IS)-ZO*(7.*HV(2,K)+HV(2,KP))*SFG(J)/8.
                ELSE
                  WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)
                END IF
   50         CONTINUE
            END IF
          END IF
          IF(IRWAK.NE.KOFL) THEN
            IF(I.EQ.IZCT) THEN
              DO 60 J=KSTL,KENL
                K=J+KK
                KP=K+IPM*NR
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)+(EV(1,K)+2.*EV
     +            (1,KP)+ EV(1,KP+IPM*NR))*SLG(J)/4.
                ELSE
                  WZ(IS)=WZ(IS)+EV(1,K)*SLG(J)
                END IF
                IF(L.NE.NZ) THEN
C                 WZ(IS)=WZ(IS)+ZO*(HV(2,K)+7.*HV(2,K+NR))*SLG(J)/8.
                  WZ(IS)=WZ(IS)+ZO*(7.*HV(2,K)+HV(2,KP))*SLG(J)/8.
                ELSE
                  WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)
                END IF
   60         CONTINUE
            END IF
          END IF
   70   CONTINUE
        IA=IA+1
      END IF
      RETURN
      END
C
      SUBROUTINE NAPL2(N,IA,NW)
C-----------------------------------------------------------------------
C NAPOLY INTEGRATION FOR MROT=1
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LCFRON,LCBACK,LCHIN,LCRBW
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /WPRM/ IZCF,IZCT,IRWAK,KENF,KENL,KSTF,KSTL,SGF,SGL,ARA,ARD
      COMMON /NAPCOM/ NRNAP1,NRNAP2
      SAVE /NAPCOM/
      COMMON /NAPPRM/ DZM,SFG(NRN),SLG(NRN)
      SAVE /NAPPRM/
C
      IPM=1
      IF(LCRBW) IPM=-1
C
      IF(NRNAP2.LT.NR) THEN
        NRNAP2=NR
        DZM=DDZ/FLOAT(M)
        DO 10 J=1,NR
          GDRM=GDR(J)/FLOAT(M)
          SFG(J)=SGF*GDRM
          SLG(J)=SGL*GDRM
   10   CONTINUE
      END IF
C
      NJDGZ=(N+1)/M*M
      IF(LCBACK) THEN
        NZF=(N+1)/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
      NZL=(N+1)/M+2
      IF(NZL.GE.NZ) NZL=NZ
      IF(N+1.NE.NJDGZ) THEN
        DO 40 I=NZF,NZL
          L=I
          IF(LCRBW) L=NZ-I+1
          LL=I
          IF(LCRBW) LL=NZSEP-I+1
          IS=IA-LL+1
          IF(IS.LE.0) GO TO 40
          IF(IS.GT.NW+2) GO TO 40
          FCT=1.0
          IF((I.EQ.IZCF).OR.(I.EQ.IZCT)) FCT=0.5
          KK=(I-MZF)*NR
          IF(I.LE.IZCF) THEN
            K=KOFF+KK
            WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
          END IF
          IF(I.GE.IZCT) THEN
            K=KOFF+KK
            WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
          END IF
          IF(IZCT.GT.IZCF) THEN
            IF((I.GE.IZCF).AND.(I.LE.IZCT)) THEN
              K=IRWAK+KK
              KP=K+IPM*NR
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*ARA/2.*FCT
              IF(L.LT.NZ) THEN
                WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)-HV(3,K-2)
     +          +3.*HV(3,KP-1)-HV(3,KP-2))*DZM*ARD/8.*FCT
              ELSE
                WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)-HV(3,K-2))*DZM*ARD/4.*FCT
              END IF
            END IF
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCF) THEN
              DO 20 J=KSTF,KENF
                K=KK+J
                KP=K+IPM*NR
                KM=K-IPM*NR
                ARJ=GR(KOFF)/(GR(J)+0.5*GDR(J))
                ARJA=(ARJ+1./ARJ)
                ARJD=(ARJ-1./ARJ)
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)-(EV(1,K)+EV(1,KP))*SFG(J)*ARJA/4.
                ELSE
                  WZ(IS)=WZ(IS)-EV(1,K)*SFG(J)*ARJA/2.
                END IF
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)-(EV(2,K)+EV(2,KP)
     +            +EV(2,K+1)+EV(2,KP+1))*SFG(J)*ARJD/8.
                ELSE
                  WZ(IS)=WZ(IS)-(EV(2,K)+EV(2,K+1))*SFG(J)*ARJD/4.
                END IF
                IF(N.EQ.NJDGZ) THEN
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)-ZO*(HV(2,K)+7.*HV
                    WZ(IS)=WZ(IS)-ZO*(7.*HV(2,K)+HV
     +              (2,KM))*SFG(J)*ARJA/16.
                  ELSE
                    WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)*ARJA/2.
                  END IF
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)+ZO*(HV(1,K)+7.*HV(1,K-NR)
C    +              +HV(1,K+1)+7.*HV(1,K+1-NR))*SFG(J)*ARJD/32.
                    WZ(IS)=WZ(IS)+ZO*(7.*HV(1,K)+HV(1,KM)
     +              +7.*HV(1,K+1)+HV(1,KM+1))*SFG(J)*ARJD/32.
                  ELSE
                    WZ(IS)=WZ(IS)+ZO*(HV(1,K)+HV(1,K+1))*SFG(J)*ARJD/4.
                  END IF
                ELSE
                  WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)*ARJA/2.
                  WZ(IS)=WZ(IS)+ZO*(HV(1,K)+HV(1,K+1))*SFG(J)*ARJD/4.
                END IF
   20         CONTINUE
            END IF
          ELSE
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCT) THEN
              DO 30 J=KSTL,KENL
                K=KK+J
                KP=K+IPM*NR
                KM=K-IPM*NR
                ARJ=GR(KOFF)/(GR(J)+0.5*GDR(J))
                ARJA=(ARJ+1./ARJ)
                ARJD=(ARJ-1./ARJ)
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)+(EV(1,K)+EV(1,KP))*SLG(J)*ARJA/4.
                ELSE
                  WZ(IS)=WZ(IS)+EV(1,K)*SLG(J)*ARJA/2.
                END IF
                IF(L.LT.NZ) THEN
                  WZ(IS)=WZ(IS)+(EV(2,K)+EV(2,KP)
     +            +EV(2,K+1)+EV(2,KP+1))*SLG(J)*ARJD/8.
                ELSE
                  WZ(IS)=WZ(IS)+(EV(2,K)+EV(2,K+1))*SLG(J)*ARJD/4.
                END IF
                IF(N.EQ.NJDGZ) THEN
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)+ZO*(HV(2,K)+7.*HV
                    WZ(IS)=WZ(IS)+ZO*(7.*HV(2,K)+HV
     +              (2,KM))*SLG(J)*ARJA/16.
                  ELSE
                    WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)*ARJA/2.
                  END IF
                  IF(L.GT.1) THEN
C                   WZ(IS)=WZ(IS)-ZO*(HV(1,K)+7.*HV(1,K-NR)
C    +              +HV(1,K+1)+7.*HV(1,K+1-NR))*SLG(J)*ARJD/32.
                    WZ(IS)=WZ(IS)-ZO*(7.*HV(1,K)+HV(1,KM)
     +              +7.*HV(1,K+1)+HV(1,KM+1))*SLG(J)*ARJD/32.
                  ELSE
                    WZ(IS)=WZ(IS)-ZO*(HV(1,K)+HV(1,K+1))*SLG(J)*ARJD/4.
                  END IF
                ELSE
                  WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)*ARJA/2.
                  WZ(IS)=WZ(IS)-ZO*(HV(1,K)+HV(1,K+1))*SLG(J)*ARJD/4.
                END IF
   30         CONTINUE
            END IF
          ELSE
          END IF
   40   CONTINUE
      ELSE
        DO 70 I=NZF,NZL
          L=I
          IF(LCRBW) L=NZ-I+1
          LL=I
          IF(LCRBW) LL=NZSEP-I+1
          IS=IA-LL+1
          IF(IS.LE.0) GO TO 70
          IF(IS.GT.NW+2) GO TO 70
          FCT=1.0
          IF((I.EQ.IZCF).OR.(I.EQ.IZCT)) FCT=0.5
          KK=(I-MZF)*NR
          IF(I.LE.IZCF) THEN
            K=KOFF+KK
            IF(L.NE.NZ) THEN
              WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,K+IPM*NR))*DZM/2.*FCT
            ELSE
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
            END IF
          END IF
          IF(I.GE.IZCT) THEN
            K=KOFF+KK
            IF(L.NE.NZ) THEN
              WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,K+IPM*NR))*DZM/2.*FCT
            ELSE
              WZ(IS)=WZ(IS)+EV(3,K)*DZM*FCT
            END IF
          END IF
          IF(IZCT.GT.IZCF) THEN
            IF((I.GE.IZCF).AND.(I.LE.IZCT)) THEN
              K=IRWAK+KK
              KP=K+IPM*NR
              KP2=KP+IPM*NR
              IF((I.GT.IZCF).AND.(I.LT.IZCT)) THEN
                IF(L.NE.NZ) THEN
                  WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,KP))*DZM*ARA/4.
                  WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)+6.*HV(3,KP-1)+3.
     +            *HV(3,KP2-1) -HV(3,K-2)-2.*HV(3,KP-2)-HV(3,KP2
     +            -2))*DZM*ARD/16.
                ELSE
                  WZ(IS)=WZ(IS)+EV(3,K)*DZM*ARA/2.
                  WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)-HV(3,K-2))*DZM*ARD/4.
                END IF
              ELSE
                IF(L.NE.NZ) THEN
                  WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,KP))*DZM*ARA/8.
                  WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)+3.*HV(3,KP-1) -HV
     +            (3,K-2)-HV(3,KP-2))*DZM*ARD/16.
                ELSE
                  WZ(IS)=WZ(IS)+EV(3,K)*DZM*ARA/4.
                  WZ(IS)=WZ(IS)-ZO*(3.*HV(3,K-1)-HV(3,K-2))*DZM*ARD/8.
                END IF
              END IF
            END IF
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCF) THEN
              DO 50 J=KSTF,KENF
                K=KK+J
                KP=K+IPM*NR
                KP2=KP+IPM*NR
                ARJ=GR(KOFF)/(GR(J)+0.5*GDR(J))
                ARJA=(ARJ+1./ARJ)
                ARJD=(ARJ-1./ARJ)
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)-(EV(1,K)+2.*EV
     +            (1,KP)+ EV(1,KP2))*SFG(J)*ARJA/8.
                ELSE
                  WZ(IS)=WZ(IS)-EV(1,K)*SFG(J)*ARJA/2.
                END IF
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)-(EV(2,K)+2.*EV
     +            (2,KP)+ EV(2,KP2)+EV(2,K+1)+2.*EV(2,KP+1)+EV
     +            (2,KP2+1)) *SFG(J)*ARJD/16.
                ELSE
                  WZ(IS)=WZ(IS)-(EV(2,K)+EV(2,K+1))*SFG(J)*ARJD/4.
                END IF
                IF(L.NE.NZ) THEN
C                 WZ(IS)=WZ(IS)-ZO*(HV(2,K)+7.*HV(2,K+NR))*SFG(J)
                  WZ(IS)=WZ(IS)-ZO*(7.*HV(2,K)+HV(2,KP))*SFG(J)
     +            *ARJA/16.
C                 WZ(IS)=WZ(IS)+ZO*(HV(1,K)+7.*HV(1,K+NR)+HV
C    +            (1,K+1) +7.*HV(1,K+1+NR))*SFG(J)*ARJD/32.
                  WZ(IS)=WZ(IS)+ZO*(7.*HV(1,K)+HV(1,KP)+7.*HV
     +            (1,K+1) + HV(1,KP+1))*SFG(J)*ARJD/32.
                ELSE
                  WZ(IS)=WZ(IS)-ZO*HV(2,K)*SFG(J)*ARJA/2.
                  WZ(IS)=WZ(IS)+ZO*(HV(1,K)+HV(1,K+1))*SFG(J)*ARJD/4.
                END IF
   50         CONTINUE
            END IF
          END IF
          IF(IRWAK.NE.KOFF) THEN
            IF(I.EQ.IZCT) THEN
              DO 60 J=KSTL,KENL
                K=KK+J
                KP=K+IPM*NR
                KP2=KP+IPM*NR
                ARJ=GR(KOFF)/(GR(J)+0.5*GDR(J))
                ARJA=(ARJ+1./ARJ)
                ARJD=(ARJ-1./ARJ)
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)+(EV(1,K)+2.*EV
     +            (1,KP)+ EV(1,KP2))*SLG(J)*ARJA/8.
                ELSE
                  WZ(IS)=WZ(IS)+EV(1,K)*SLG(J)*ARJA/2.
                END IF
                IF (L.LT.NZ-1) THEN
                  WZ(IS)=WZ(IS)+(EV(2,K)+2.*EV
     +            (2,KP)+ EV(2,KP2)+EV(2,K+1)+2.*EV(2,KP+1)+EV
     +            (2,KP2+1)) *SLG(J)*ARJD/16.
                ELSE
                  WZ(IS)=WZ(IS)+(EV(2,K)+EV(2,K+1))*SLG(J)*ARJD/4.
                END IF
                IF(L.NE.NZ) THEN
C                 WZ(IS)=WZ(IS)+ZO*(HV(2,K)+7.*HV(2,K+NR))*SLG(J)
                  WZ(IS)=WZ(IS)+ZO*(7.*HV(2,K)+HV(2,KP))*SLG(J)
     +            *ARJA/16.
C                 WZ(IS)=WZ(IS)-ZO*(HV(1,K)+7.*HV(1,K+NR)+HV
C    +            (1,K+1) +7.*HV(1,K+1+NR))*SLG(J)*ARJD/32.
                  WZ(IS)=WZ(IS)-ZO*(7.*HV(1,K)+HV(1,KP)+7.*HV
     +            (1,K+1) + HV(1,KP+1))*SLG(J)*ARJD/32.
                ELSE
                  WZ(IS)=WZ(IS)+ZO*HV(2,K)*SLG(J)*ARJA/2.
                  WZ(IS)=WZ(IS)-ZO*(HV(1,K)+HV(1,K+1))*SLG(J)*ARJD/4.
                END IF
   60         CONTINUE
            END IF
          END IF
   70   CONTINUE
        IA=IA+1
      END IF
      RETURN
      END
C
      SUBROUTINE JMPBM(N,IA,NW)
C-----------------------------------------------------------------------
C INTEGRATION FOR WAKE POTENTIALS WHEN A BEAM JUMPS
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LCFRON,LCBACK,LCHIN,LCRBW
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WPL/ MROT
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PINT/ KMIN
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
C
      IPM=1
      IF(LCRBW) IPM=-1
      DZM=DDZ/FLOAT(M)
      DZMO=DZM*ZO
C
      NJDGZ=(N+1)/M*M
      IF(LCBACK) THEN
        NZF=(N+1)/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
      NZL=(N+1)/M+2
      IF(NZL.GE.NZ) NZL=NZ
      IF(N+1.NE.NJDGZ) THEN
        DO 10 I=NZF,NZL
          J=I
          IF(LCRBW) J=NZ-I+1
          JJ=I
          IF(LCRBW) JJ=NZSEP-I+1
          IS=IA-JJ+1
          IF(IS.LE.0) GO TO 10
          IF(IS.GT.NW+2) GO TO 10
          K=KMIN+(I-MZF)*NR
          WZ(IS)=WZ(IS)+EV(3,K)*DZM
          IF(J.LT.NZ) THEN
            WY(IS)=WY(IS)+(EV(2,K)+EV(2,K+IPM*NR))*DZM/2.
          ELSE
            WY(IS)=WY(IS)+EV(2,K)*DZM
          END IF
          IF(N.EQ.NJDGZ) THEN
            IF(J.GT.1) THEN
C             WY(IS)=WY(IS)+(HV(1,K)+7.*HV(1,K-NR))*DZMO/8.
              WY(IS)=WY(IS)+(7.*HV(1,K)+HV(1,K-IPM*NR))*DZMO/8.
            ELSE
              WY(IS)=WY(IS)+HV(1,K)*DZMO
            END IF
          ELSE
            WY(IS)=WY(IS)+HV(1,K)*DZMO
          END IF
   10   CONTINUE
      ELSE
        DO 20 I=NZF,NZL
          J=I
          IF(LCRBW) J=NZ-I+1
          JJ=I
          IF(LCRBW) JJ=NZSEP-I+1
          IS=IA-JJ+1
          IF(IS.LE.0) GO TO 20
          IF(IS.GT.NW+2) GO TO 20
          K=KMIN+(I-MZF)*NR
          KP=K+IPM*NR
          IF (J.LT.NZ-1) THEN
            WY(IS)=WY(IS)+(EV(2,K)+2.*EV(2,KP)+ EV
     +      (2,KP+IPM*NR))*DZM/4.
          ELSE
            WY(IS)=WY(IS)+EV(2,K)*DZM
          END IF
          IF(J.NE.NZ) THEN
            WZ(IS)=WZ(IS)+(EV(3,K)+EV(3,KP))*DZM/2.
C           WY(IS)=WY(IS)+(HV(1,K)+7.*HV(1,K+NR))*DZMO/8.
            WY(IS)=WY(IS)+(7.*HV(1,K)+HV(1,KP))*DZMO/8.
          ELSE
            WZ(IS)=WZ(IS)+EV(3,K)*DZM
            WY(IS)=WY(IS)+HV(1,K)*DZMO
          END IF
   20   CONTINUE
        IA=IA+1
      END IF
      RETURN
      END
C
      SUBROUTINE CONBM(N,IA,IB,NW)
C-----------------------------------------------------------------------
C INTEGRATION FOR WAKE POTENTIALS WHEN A BEAM MOVES SMOOTHLY
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NWN=500000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      LOGICAL LCFRON,LCBACK,LCHIN,LCRBW
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /WZWY/ WZ(NWN),WY(NWN),WT(NWN),WP(NWN)
      COMMON /WPL/ MROT
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PINT/ KMIN
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CRBW/ LCRBW,ZSEP,NZSEP
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
C
      DDZO=DDZ*ZO
C
      NJDGZ=((N+1)/M)*M+M/2
      NJDGX=((N+1)/M)*M
      IF(LCBACK) THEN
        NZF=(N+1)/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
      NZL=(N+1)/M+2
      IF(NZL.GE.NZ) NZL=NZ
      IF(N+1.NE.NJDGZ) GO TO 20
      IA=IA+1
      DO 10 I=NZF,NZL
        J=I
        IF(LCRBW) J=NZSEP-I+1
        IS=IA-J+1
        IF (IS.LE.0) GO TO 10
        IF (IS.GT.NW+2) GO TO 10
        K=KMIN+(I-MZF)*NR
        WZ(IS)=WZ(IS)+EV(3,K)*DDZ
        WY(IS)=WY(IS)+HV(1,K)*DDZO
   10 CONTINUE
C-----------------------------------------------------------------------
   20 IF(N+1.NE.NJDGX) GO TO 40
      IB=IB+1
      DO 30 I=NZF,NZL
        J=I
        IF(LCRBW) J=NZSEP-I+1
        IS=IB-J+1
        IF(IS.LE.0) GO TO 30
        IF(IS.GT.NW+2) GO TO 30
        K=KMIN+(I-MZF)*NR
        WY(IS)=WY(IS)+EV(2,K)*DDZ
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE HEFLD(N)
C-----------------------------------------------------------------------
C CALCULATION OF ELECTROMAGNETIC FIELDS FOR M=1
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
C
      ZOM=ZO/FLOAT(M)
      YOM=1./(ZO*FLOAT(M))
      IF(LCFRON) THEN
        NZL=N/M+2
        IF(NZL.GE.NZ) NZL=NZ
      ELSE
        NZL=NZ
      END IF
      IF(LCBACK) THEN
        NZF=N/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
C
      IF ((NZF.NE.NZF1).OR.(NZL.NE.NZL1)) CALL BCDET(NZF,NZL,NZF1,NZL1)
C
      IF(NZF1 .LE. 1) THEN
        CALL HEFF1(N)
      ENDIF
C
      IF(NZL1 .GE. NZ-1) THEN
        CALL HEFF2
      ENDIF
C
      IF (NZF .NE. NZF1) THEN
        IF((NZL.EQ.MZL+1).AND.(NZL.NE.NZL1)) THEN
          CALL HEFF3
        END IF
      ENDIF
      CALL HEFF4
      CALL HEFF5
      CALL HEFF6
      CALL HEFF7(N)
      CALL HEFF8(N)
      NZF1=NZF
      NZL1=NZL
      END
C
      SUBROUTINE HEFF1(N)
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
C
      IF(LCHIN) THEN
        N1=N-1-M/2
        N2=N-1+2*M-M/2
        CALL CHAGE(1,N1,QK1)
        CALL CHAGE(1,N2,QK2)
      ELSE
        N1=((N-2)/M)
        N2=((N-2+2*M)/M)
        QK1=QKST(N1+1)
        QK2=QKST(N2+1)
      END IF
C
      HPI=CLGHT/(PI*2.*GR(JOFF))
      AB=(GR(JOFF)/GR(KOFF))**2
C
      DO 30 JR=1,KOFF-1
        K=JR
        DO 20 L=1,2
          DO 10 J=NPT-M,NPT
            HVSTK(JR,L,J-1)=HVSTK(JR,L,J)
   10     CONTINUE
          ARP=(2.*GR(JOFF)/(GR(JR)+GR(JR+1)))**2
          IF(JR.GT.1) ARR=(GR(JOFF)/GR(JR))**2
          IF(JR.GE.JOFF) THEN
            IF(L.EQ.1) HVSTK(JR,L,NPT)=HV(L,K)-QK1*HPI*(AB-ARR) +QK2
     +      *HPI *(AB-ARR)
            IF(L.EQ.2) HVSTK(JR,L,NPT)=HV(L,K)-QK1*HPI*(AB+ARP) +QK2
     +      *HPI *(AB+ARP)
          ELSE
            HVSTK(JR,L,NPT)=HV(L,K)-QK1*HPI*(AB-1.)+QK2*HPI*(AB-1.)
          END IF
   20   CONTINUE
   30 CONTINUE
      END
C
      SUBROUTINE HEFF2
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
C
      IF(NZL1.EQ.NZ) THEN
        KK = (NZL1-MZF-1)*NR
      ELSE
        KK = (NZL1-MZF)*NR
      END IF
      DO 30 JR=1,KOFL-1
        K = KK + JR
        DO 20 L=1,2
          DO 10 J=NPT-M,NPT
            HVSTL(JR,L,J-1)=HVSTL(JR,L,J)
   10     CONTINUE
          HVSTL(JR,L,NPT)=HV(L,K)
   20   CONTINUE
   30 CONTINUE
      END
C
      SUBROUTINE HEFF3
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
C
      K3=(NZL1-NZF1)*NR+1
      DO 10 L=(NZF1-MZF)+2,NZL1-MZF+1
        LL=L-(NZF1-MZF+2)+1
        K1=(L-1)*NR
        K2=(LL-1)*NR
        DO 10 J=1,NR
          K=K1+J
          KR=K2+J
          DO 10 I=1,3
            HV(I,KR)=HV(I,K)
            EV(I,KR)=EV(I,K)
   10 CONTINUE
      DO 20 K=K3,NMESH
        DO 20 I=1,3
          HV(I,K)=0.
          EV(I,K)=0.
   20 CONTINUE
      MZF=NZF
      MZL=MZF+NMESH/NR-1
      END
C
      SUBROUTINE HEFF4
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
      COMMON / CNTCOM / NRCNT4, NRCNT7, NRCNTM
      SAVE / CNTCOM /
      COMMON / LOCCM1 / CX1(NRN), CX2(NRN), CX3(NRN), CX4(NRN), CX5
     +(NRN), CX6(NRN)
      SAVE /LOCCM1/
C
      IF (NR .GT. NRCNT4)  THEN
        NRCNT4 = NR
        DO 10 JR=1,NR
          CX1(JR) = DDZ/GDR(JR)
   10   CONTINUE
        DO 20 JR=2,NR-1
          CX2(JR) = DDZ/GR(JR)
          CX3(JR) = DDZ/(GR(JR)+GR(JR+1))
          BX7 = GR(JR+1)-GR(JR)
          CX5(JR) = GR(JR+1)/BX7
          CX6(JR) = GR(JR)/BX7
   20   CONTINUE
        CX3(NR) = CX3(NR-1)
        CX5(NR) = CX5(NR-1)
        CX6(NR) = CX6(NR-1)
      ENDIF
C
      DO 30 IZ=NZF,NZL-1
        K=(IZ-MZF)*NR+1
        HV(2,K)=HV(2,K)-YOM*(CX1(1)*(EV(3,K)-EV(3,K+1))
     +  +EV(1,K+NR)-EV(1,K))
        HV(1,K)=HV(2,K)
        HV(3,K)=HV(3,K)-YOM*CX1(1)*2.*(EV(2,K+1)+EV(1,K))
   30 CONTINUE
      END
C
      SUBROUTINE HEFF5
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
      COMMON / LOCCM1 / CX1(NRN), CX2(NRN), CX3(NRN), CX4(NRN), CX5
     +(NRN), CX6(NRN)
C
      DO 20 IZ=NZF,NZL-1
        KK = (IZ-MZF)*NR
        DO 10 JR=2,NR
          K = KK + JR
          IF(IFLG(K).EQ.5) GO TO 10
          BTP1 = HV(1,K)-SJDG(3,K)*YOM*(EV(2,K)-EV(2,K+NR) -CX2(JR)*EV
     +    (3,K))
          BTP3 = HV(3,K)-SJDG(1,K)*YOM*CX3(JR)*2.*(EV(2,K+1)*CX5(JR)-EV
     +    (2,K)*CX6(JR)+EV(1,K))
          FCT=1.
          IF(IFLG(K).NE.0) FCT=2.
          BTP2 = HV(2,K)-YOM*FCT*(CX1(JR)*(EV(3,K)-EV(3,K+1))+EV(1,K+NR)
     +    -EV(1,K))
          HV(1,K) = BTP1
          HV(3,K) = BTP3
          HV(2,K) = BTP2
   10   CONTINUE
   20 CONTINUE
      END
C
      SUBROUTINE HEFF6
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
      COMMON / LOCCM1 / CX1(NRN), CX2(NRN), CX3(NRN), CX4(NRN), CX5
     +(NRN), CX6(NRN)
C
      IF(NZL .EQ. NZ) THEN
        KK=(NZL-MZF)*NR
        K=KK+1
        NPTM=NPT+1-M
        HV(1,K)=HVSTL(1,1,NPTM)
        HV(2,K)=HVSTL(1,2,NPTM)
        HV(3,K)=HV(3,K)-YOM*CX1(1)*2.*(EV(2,K+1)+EV(1,K))
        DO 10 JR=2,KOFL-1
          K=KK+JR
          IF(SJDG(3,K) .NE. 0.) HV(1,K)=HVSTL(JR,1,NPTM)
          HV(3,K)=HV(3,K)-SJDG(1,K)*YOM*CX3(JR)*2.*(EV(2,K+1)*CX5(JR)
     +    -EV(2,K)*CX6(JR)+EV(1,K))
          HV(2,K)=HVSTL(JR,2,NPTM)
   10   CONTINUE
      ENDIF
      END
C
      SUBROUTINE HEFF7(N)
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
      COMMON / CNTCOM / NRCNT4, NRCNT7, NRCNTM
      SAVE / CNTCOM /
      COMMON / LOCCM2 / DX1(NRN), DX2(NRN), DX3(NRN), DX4(NRN), DX5
     +(NRN), DX6(NRN)
      SAVE /LOCCM2/
C
C     T1265 = 12. / 65.
      NPTM=NPT+1-M
      IF (NR .GT. NRCNT7)  THEN
        NRCNT7 = NR
        DO 10 JR=1,NR
          DX1(JR) = DDZ/GDR(JR)
          DX2(JR) = DDZ/(GR(JR)+GDR(JR)/2.)
   10   CONTINUE
        DO 20 JR=2,NR
          DX3(JR) = 2.*DDZ/(GDR(JR)+GDR(JR-1))
          DX4(JR) = DDZ/(GR(JR)+(GDR(JR)-GDR(JR-1))/4.)
          BX7 = DDZ/(GR(JR)-GDR(JR-1)/2.)
          DX5(JR) = DX3(JR) / DX2(JR)
          DX6(JR) = DX3(JR) / BX7
   20   CONTINUE
      ENDIF
C
      IF(NZF .LE. 1) THEN
        DO 30 JR=1,KOFF-1
          K=JR
          IF(JR.EQ.1) THEN
            EV(1,K)=EV(1,K)+ZOM*(HVSTK(JR,2,NPTM)
     +      -HV(2,K)+DX2(JR)*HV(3,K))
            EV(2,K)=EV(2,K)+ZOM*(HV(1,K)-HVSTK(JR,1,NPTM)
     +      -DX2(JR)*HV(3,K))
          ELSE
            EV(1,K)=EV(1,K)+SJDG(1,K)*ZOM*((HVSTK(JR,2,NPTM)
     +      -HV(2,K)) +DX2(JR)*HV(3,K))
            EV(2,K)=EV(2,K)+SJDG(2,K)*ZOM*((HV(1,K)-HVSTK
     +      (JR,1,NPTM)) +DX3(JR)*(HV(3,K-1)-HV(3,K)))
            EV(3,K)=EV(3,K)+SJDG(3,K)*ZOM*DX4(JR)*(HV(2,K)*DX5(JR)
     +      -HV(2,K-1) *DX6(JR) -HV(1,K))
          END IF
          IF(JR .EQ. JOFF .AND. SJDG(3,K) .NE. 0.) THEN
            IF (LCHIN) THEN
              N1=N-M/2
              CALL CHAGE(1,N1,QK)
            ELSE
C     N1=(N-1)/M*M
              N1=(N-1)/M
              QK=QKST(N1+1)
            END IF
C     EV(3,K)=EV(3,K)-QK*BX/(EPS*FLOAT(M)*DDR*PI*(FJR-1.))
            EV(3,K)=EV(3,K)-QK*DX1(JR)/(EPS*FLOAT(M)*PI*GR(JR))
          ENDIF
   30   CONTINUE
      ENDIF
      END
C
      SUBROUTINE HEFF8(N)
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /ELINE/ NZF,NZL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / LOCCOM /
     +AB, ARP, ARR, FCT, HPI, KR, LL, N1, N2, QK, QK1, QK2, YOM, ZOM
      COMMON / LOCCM2 / DX1(NRN), DX2(NRN), DX3(NRN), DX4(NRN), DX5
     +(NRN), DX6(NRN)
C
C     T1265 = 12. / 65.
      DO 20 IZ=NZF+1,NZL
        KK = (IZ-MZF)*NR
        K = KK + 1
        IF (LCHIN) THEN
          N1=N-M/2
          CALL CHAGE(IZ,N1,QK)
        ELSE
          N1=(N-1)/M
          N2=N1-(IZ-1)+1
          QK=0.0
          IF((N2.GE.1).AND.(N2.LE.IBT)) QK=QKST(N2)
        END IF
        IF(IFLG(K) .NE. 5) THEN
          JR = 1
          EV(1,K) = EV(1,K)+ZOM*(HV(2,K-NR)-HV(2,K)+DX2(JR)
     +    *HV(3,K))
          EV(2,K) = EV(2,K)+ZOM*(HV(1,K)-HV(1,K-NR)-DX2(JR)
     +    *HV(3,K))
        ENDIF
        DO 10 JR=2,NR
          K = KK + JR
          IF(IFLG(K).EQ.5) GO TO 10
          BTP1 = EV(1,K)+SJDG(1,K)*ZOM*((HV(2,K-NR) -HV(2,K))+DX2 (JR)
     +    *HV(3,K))
          BTP2 = EV(2,K)+SJDG(2,K)*ZOM*((HV(1,K)-HV(1,K-NR)) +DX3 (JR)
     +    *(HV(3,K-1)-HV(3,K)))
          BTP3 = EV(3,K)+SJDG(3,K)*ZOM*DX4(JR)*(HV(2,K) *DX5(JR) -HV
     +    (2,K-1)*DX6(JR)-HV(1,K))
          EV(1,K) = BTP1
          EV(2,K) = BTP2
          EV(3,K) = BTP3
          IF(JR .EQ. JOFF .AND. SJDG(3,K) .NE. 0.) THEN
            EV(3,K)=EV(3,K)-QK*DX1(JR)/(EPS*FLOAT(M)*PI*GR(JR))
          ENDIF
   10   CONTINUE
   20 CONTINUE
      END
C
      SUBROUTINE MEFLD(N)
C-----------------------------------------------------------------------
C CALCULATION OF ELECTROMAGNETIC FIELDS FOR MROT=0
C-----------------------------------------------------------------------
      PARAMETER (NMESH=10000000)
      PARAMETER (NRN=20000)
      PARAMETER (NBP=1000)
      PARAMETER (NF=100)
      PARAMETER (NWN=500000)
      PARAMETER (NPT=20)
      LOGICAL LCFRON,LCBACK,LCHIN
      COMMON /LATC/ NR,NZ,NZP,IRR(NBP),IZZ(NBP),JOFF,DDZ,DDR(NF),SIG,M
      COMMON /LAT2/ RBDR(NF),RMARK(NF),GR(NRN),GDR(NRN),IMAX,IRMAX,NVC
      COMMON /HE/ HV(3,NMESH),EV(3,NMESH)
      COMMON /FLG/ IFLG(NMESH),SJDG(3,NMESH)
      COMMON /STK/ HVSTK(NRN,2,NPT),HVSTL(NRN,2,NPT)
      COMMON /WIND/ NT
      COMMON /WIN2/ LCFRON,LCBACK,LCHIN
      COMMON /NZIN/ NZF1,NZL1
      COMMON /LFRBK/ LFR,LBK,MZF,MZL
      COMMON /ELINE/ NZF,NZL
      COMMON /PIPE/ KOFF,KOFL
      COMMON /CONS/ PI,EPS,CLGHT,CEPS,ZO
      COMMON /CALS/ SUMQ,BSIG
      COMMON /QKAR/ QKST(NWN),IBT,QKMAX
C
      COMMON / CNTCOM / NRCNT4, NRCNT7, NRCNTM
      SAVE / CNTCOM /
      COMMON / LOCCM1 / EX1(NRN), EX4(NRN), EX5(NRN), EX6(NRN),
     +EDUM1(NRN), EDUM2(NRN)
      SAVE /LOCCM1/
      COMMON / PARCOM / ZOM, YOM, CFCT
      SAVE /PARCOM/
C
      NPTM=NPT+1-M
      IF (NRCNTM .LT. NR)  THEN
        NRCNTM = NR
        ZOM=ZO/FLOAT(M)
        YOM=1./(ZO*FLOAT(M))
        CFCT=EPS*FLOAT(M)*PI
        DO 10 JR=1,NR
          EX1(JR)=DDZ/GDR(JR)
          BX2=DDZ/(GR(JR)+GDR(JR)/2.)
          IF(JR.EQ.1) GO TO 10
          BX7=DDZ/(GR(JR)-GDR(JR-1)/2.)
          BX3=2.*DDZ/(GDR(JR)+GDR(JR-1))
          EX4(JR)=DDZ/(GR(JR)+(GDR(JR)-GDR(JR-1))/4.)
          EX5(JR)=BX3/BX2
          EX6(JR)=BX3/BX7
   10   CONTINUE
      ENDIF
C
      IF(LCFRON) THEN
        NZL=N/M+2
        IF(NZL.GE.NZ) NZL=NZ
      ELSE
        NZL=NZ
      END IF
      IF(LCBACK) THEN
        NZF=N/M-NT/M+NZ-2
        IF(NZF.LE.1) NZF=1
      ELSE
        NZF=1
      END IF
C
      IF ((NZF.NE.NZF1).OR.(NZL.NE.NZL1)) CALL BCDET(NZF,NZL,NZF1,NZL1)
C
      IF(NZF1.GT.1) GO TO 50
      N1=(N-2)/M
      N2=(N-2+2*M)/M
      QK1=QKST(N1+1)
      QK2=QKST(N2+1)
      DO 40 JR=1,KOFF-1
        K=JR
        HPI=CLGHT/(PI*(GR(JR)+GR(JR+1)))
        DO 30 L=1,2
          DO 20 J=NPT-M,NPT
            HVSTK(JR,L,J-1)=HVSTK(JR,L,J)
   20     CONTINUE
          IF(JR.GE.JOFF) THEN
            HVSTK(JR,L,NPT)=HV(L,K) + (QK2 - QK1) * HPI
          ELSE
            HVSTK(JR,L,NPT)=HV(L,K)
          END IF
   30   CONTINUE
   40 CONTINUE
   50 CONTINUE
C
      IF(NZL1.LT.NZ-1) GO TO 90
      IF(NZL1.EQ.NZ) THEN
        KK = (NZL1-MZF-1)*NR
      ELSE
        KK = (NZL1-MZF)*NR
      ENDIF
      DO 80 JR=1,KOFL-1
        K = KK + JR
        DO 70 L=1,2
          DO 60 J=NPT-M,NPT
            HVSTL(JR,L,J-1)=HVSTL(JR,L,J)
   60     CONTINUE
          HVSTL(JR,L,NPT)=HV(L,K)
   70   CONTINUE
   80 CONTINUE
   90 CONTINUE
C
      IF (NZF.EQ.NZF1) GO TO 120
      IF((NZL.EQ.MZL+1).AND.(NZL.NE.NZL1)) THEN
C
        DO 100 L=(NZF1-MZF)+2,NZL1-MZF+1
          L1=(L-1)*NR
          LL=L-(NZF1-MZF+2)+1
          LL1=(LL-1)*NR
          DO 100 J=1,NR
            K=L1+J
            KR=LL1+J
            DO 100 I=1,3
              HV(I,KR)=HV(I,K)
              EV(I,KR)=EV(I,K)
  100   CONTINUE
        K1=(NZL1-NZF1)*NR+1
        DO 110 K=K1,NMESH
          DO 110 I=1,3
            HV(I,K)=0.
            EV(I,K)=0.
  110   CONTINUE
        MZF=NZF
        MZL=MZF+NMESH/NR-1
      END IF
  120 CONTINUE
C
      DO 140 JR=1,NR
        BX=DDZ/GDR(JR)
        DO 130 IZ=NZF,NZL-1
          K=(IZ-MZF)*NR+JR
          IF(IFLG(K).EQ.5) GO TO 130
          FCT=1.
          IF(IFLG(K).NE.0) FCT=2.
          HV(2,K)=HV(2,K)-YOM*FCT*(BX*(EV(3,K)-EV(3,K+1)) +EV(1,K+NR)-EV
     +    (1,K))
  130   CONTINUE
  140 CONTINUE
      IF(NZL.NE.NZ) GO TO 160
      KK=(NZL-MZF)*NR
      DO 150 JR=1,KOFL-1
        K=KK+JR
        HV(2,K)=HVSTL(JR,2,NPTM)
  150 CONTINUE
  160 CONTINUE
C
      N1=(N-1)/M
      QK=QKST(N1+1)
C
      IF(NZF .LE. 1) THEN
        DO 170 K=1,KOFF-1
          EV(1,K)=EV(1,K)+SJDG(1,K)*ZOM*((HVSTK(K,2,NPTM)
     +    -HV(2,K)))
          IF(K.GT.1) THEN
            EV(3,K)=EV(3,K)+SJDG(3,K)*ZOM*EX4(K)*(HV(2,K)*EX5(K) -HV
     +      (2,K-1)*EX6(K))
          ELSE
            EV(3,K)=EV(3,K)+ZOM*EX1(1)*4.*HV(2,K)
          END IF
          IF(K.NE.JOFF) GO TO 170
          IF(SJDG(3,K) .EQ. 0. .AND.(JOFF.GT.1)) GO TO 170
          IF(K.GT.1) THEN
            EV(3,K)=EV(3,K)-QK*EX1(K)/(CFCT*GR(K))/2.
          ELSE
            EV(3,K)=EV(3,K)-4.*QK*EX1(1)/(CFCT*GDR(K))
          END IF
  170   CONTINUE
      ENDIF
      DO 190 IZ=NZF+1,NZL
        N2=N1-(IZ-1)+1
        IF((N2.GE.1).AND.(N2.LE.IBT)) THEN
          QK=QKST(N2)
        ELSE
          QK=0.0
        ENDIF
        KK = (IZ - MZF) * NR
        K = KK + 1
        IF(IFLG(K) .NE. 5) THEN
          EV(1,K)=EV(1,K)+SJDG(1,K)*ZOM*((HV(2,K-NR) -HV(2,K)))
          EV(3,K)=EV(3,K)+ZOM*EX1(1)*4.*HV(2,K)
          IF(JOFF .EQ. 1) THEN
            EV(3,K)=EV(3,K)-4.*QK*EX1(1)/(CFCT*GDR(1))
          ENDIF
        ENDIF
        DO 180 JR = 2, NR
          K = KK + JR
          IF(IFLG(K).EQ.5) GO TO 180
          BTP1 = EV(1,K)+SJDG(1,K)*ZOM*((HV(2,K-NR) -HV(2,K)))
          BTP3 = EV(3,K)+SJDG(3,K)*ZOM*EX4(JR)*(HV(2,K)*EX5(JR) -HV
     +    (2,K-1)*EX6(JR))
          EV(1,K) = BTP1
          EV(3,K) = BTP3
  180   CONTINUE
        K = KK + JOFF
        IF(IFLG(K) .NE. 5 .AND. SJDG(3,K) .NE. 0.) THEN
          EV(3,K) = EV(3,K)
     +    - QK*EX1(JOFF)/(CFCT*GR(JOFF))/2.
        ENDIF
  190 CONTINUE
      NZF1=NZF
      NZL1=NZL
      END
C
      SUBROUTINE ASFILE(NME,NMEO,NMEP,NMET)
C-----------------------------------------------------------------------
C RETURN FILE NAMES FOR OUTPUT FILES
C-----------------------------------------------------------------------
	CHARACTER*256 NMEO,NMEP,NMET
	CHARACTER*1 NME(256),DOT,OCH,UCH,TCH,PCH
	DATA DOT/'.'/,OCH/'o'/,UCH/'u'/,TCH/'t'/,PCH/'p'/
C 	PRINT '(A256)', NMEA
	DO 1000 I=1,256
	IF(NME(I).NE.DOT) GO TO 1000
	N=I
 1000	CONTINUE
	WRITE(NMEO, 2000) (NME(I), I=1,N-1), DOT,OCH,UCH,TCH
	WRITE(NMEP, 2000) (NME(I), I=1,N-1), DOT,PCH,OCH,TCH
	WRITE(NMET, 2000) (NME(I), I=1,N-1), DOT,TCH,OCH,PCH
 2000 FORMAT(256A1)
C	PRINT *, '*** THE ASSIGNED OUTPUT FILES ARE AS FOLLOWS: ***'
C     PRINT '(A12,A256)', 'OUTPUT FILE=',NMEO
C 	PRINT '(A12,A256)', 'POTWAK FILE=',NMEP
C 	PRINT '(A12,A256)', 'TOPDRW FILE=',NMET
	PRINT *, '> ABCI IS NOW RUNNING. PLEASE WAIT...'
	RETURN
	END
C
      SUBROUTINE CPUTIM(TT,KTYPE)
C-----------------------------------------------------------------------
C RETURN THE CPU TIME SO FAR IN SECONDS
C-----------------------------------------------------------------------
      COMMON /CPUT/ T2
      SAVE /CPUT/
      IF(KTYPE-2) 10,20,30
C-SET TIMER TO START THE CLOCK------------------------------------------
   10 CONTINUE
      CALL TIMEST(TLIM)
      T2=0.
      RETURN
C-CPU TIME USED IN A RUN------------------------------------------------
   20 CONTINUE
      CALL TIMEX(T1)
      TT=T1-T2
      T2=T1
      RETURN
C-TOTAL CPU TIME USED---------------------------------------------------
   30 CONTINUE
      CALL TIMEX(TT)
      RETURN
      END
C
      SUBROUTINE CPULFT(TSOS,TMAX,LSOS)
C-----------------------------------------------------------------------
C CHECK IF THE CPU TIME LEFT (IN SECONDS) IS STILL LARGER THAN TSOS
C USE TIMEL IF THE CERN LIBRARY IS AVAIRABLE
C USE TIMEX IF YOUR CERN LIBRARY HAS NO TIMAL BUT TIMAX
C-----------------------------------------------------------------------
      LOGICAL LSOS
      LSOS=.FALSE.
      TMAX=TMAX
C     TT=TMAX-TIMEX(TT)
      CALL TIMEL(TT)
      IF((TT.GT.0.).AND.(TT.LT.TSOS)) LSOS=.TRUE.
      RETURN
      END
C
      SUBROUTINE DATIME(NDATE,NTIME)
C-----------------------------------------------------------------------
C RETURN THE DATE AND TIME OF THE DAY
C-----------------------------------------------------------------------
C FOR IBM USERS
      CHARACTER*8 NDATE,NTIME
C FOR VAX USERS
C     REAL*8 NDATE,NTIME
      CALL DATIMH(NDATE,NTIME)
      RETURN
      END
C------------------------------------------------------------------------
C File name:
C      DATLIB fortran library for HP users
C Created:
C      Yong Ho Chin, LBL 1/3/94
C Description:
C      Just copy it to the end of the ABCI source code or link its
C      object module with that of ABCI. You must also link with libU77.
C
C      If IDATE, TIME and/or SECNDS are not available, edit this library
C      according to the instructions written before each routine.
C------------------------------------------------------------------------
      SUBROUTINE DATIMH (NDATE,NTIME)
C
      CHARACTER*8 NDATE,NTIME
C     DATA MONTH,IDAY,IYEAR/1,1,94/
      NTIME='00.00.00'
C IF IDATE IS NOT AVAILABLE, COMMENT OUT THE NEXT LINE
      CALL IDATE(MONTH,IDAY,IYEAR)
C IF TIME IS NOT AVAILABLE, COMMENT OUT THE NEXT LINE
      CALL TIME(NTIME)
      WRITE(NDATE,1000) IDAY,MONTH,IYEAR
 1000 FORMAT(I2,'/',I2,'/',I2)
      RETURN
      END
C
      SUBROUTINE TIMEST(TIMLIM)
C
      TIMLIM=TIMLIM
      RETURN
C
      ENTRY TIMEX (TRES)
C
      TRES=ETIME(x)
      RETURN
      END
C
      SUBROUTINE TIMEL(TRES)
C
      TRES=3600.
      RETURN
      END
