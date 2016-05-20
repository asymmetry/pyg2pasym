C      Defined already in poltail.f
C----------------------------------
      BLOCK DATA FOR_INELASTIC
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)

       !DATA AMT/2.8094D0/,     TARA/3D0/, TARZ/2D0/
       DATA AMT/0.938272D0/, TARA/1D0/, TARZ/1D0/                 ! CPROT


      DATA AMM/2.7928456D0/
      DATA AMN/-1.913148D0/
      DATA CHBAR/.197328D0/
      DATA BARN/.389379D6/
      DATA AML/.511000D-3/
      DATA AML2/.261112D-6/
      DATA AL2/.522240D-6/
      DATA FERMOM/.164D0/     !  - ? FERMI MOMENTUM IN 3HE
      DATA ISF20/4/
      DATA PI/3.1415926D0/,PI2/9.869604D0/,ALFA/.729735D-2/
      DATA AMC2/1.151857D0/
      DATA AMP/.938272D0/,AMH/.938272D0/
      DATA I2/1,1,1,2,3,3,1,2/,I1/3,3,4,4,3,3,3,3/
      END





      SUBROUTINE POLSIG_2(E0,EP,TH_RAD,IELAS_IN,IPOL,SIGMA_BORN,SIGMA_RAD)
*----------------------------------------------------------------------
*     CALCULATION OF THE CROSS SECTION
*
*     INPUT
*       E0 = INCIDENT ELECTRON ENERGY IN MEV
*       EP = SCATTERED ELECTRON ENERGY IN MEV
*       TH_RAD = SCATTERING ANGLE IN RADIANS
*
*
*     KINEMATIC REGION
*
*       IELAS_IN = 0 : QUASI-ELASTIC AND INELASTIC
*       IELAS_IN = 1 : ELASTIC
*
*     POLARIZATION SPECIFICATION
*
*       IPOL = 0 : PARALLEL, UNPOLARIZED
*       IPOL = 1 : PARALLEL, POLARIZED
*       IPOL = 2 : PERPENDICULAR, UNPOLARIZED
*       IPOL = 3 : PERPENDICULAR, POLARIZED
*
*     THE OUTPUT SIGMA_BORN AND SIGMA_RAD ARE IN NB/MEV.SR
*
*---- FOR XJACOB:
*---- FACTOR 1.0D-3 CONVERTS NB/GEV.SR TO NB/MEV.SR FOR INELASTIC AND QUASI-
*---- ELASTIC CROSS SECTION AND ELASTIC TAIL. FOR ELASTIC CROSS SECTION, IT
*---- SIMPLY CONVERTS NB TO UB. AND ALSO REMEMBER THAT ELASTIC CROSS SECTION
*---- WILL LACK RECOIL FACTOR (E'/E)
*
*----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 E0,EP,TH_RAD,SIGMA_BORN,SIGMA_RAD

      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     &     SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
      COMMON/TAIL_2/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH

      COMMON/CHOI_2/W2_ELAS
      COMMON/SOFTPHOTON_2/EXTAI2
      COMMON/ELASTIC_2/IELAS
C
C----------------------------------------------------------------------
C

      IELAS = IELAS_IN  ! CMASS

      Q2    = 4.0*E0*EP*SIN(0.5*TH_RAD)**2*1.0E-6
      XNU   = (E0-EP)*1.0E-3
      
      IF (IELAS.EQ.0) THEN
         W2_N    = AMH*(AMH + 2.0*XNU) - Q2   ! JUST FOR REFERENCE
         W_N     = SQRT(W2_N)*1000.0
         W2      = AMT*(AMT + 2.0*XNU) - Q2
         XS      = Q2/(2.0*AMT*XNU)
         YS      = XNU/(E0*1.0E-3)

         
         XJACOB  = 2.0D0*PI*AMH*YS/(1.0-YS)   ! XJACOB = 1E-3 * (1-YS)/(2PI*MH*YS)
         XJACOB  = 1.0D-3/XJACOB              !        = (1-YS)/2PI*MH*YS) [1/MEV]

         Y_ELAS  = 1.0/(1.0 + AMT/(2.0*E0*1.0E-3*SIN(0.5*TH_RAD)**2))

         W2_ELAS = AMT**2                     ! CDEBUG ADD BREAKUP
         
      ELSE IF (IELAS.EQ.1) THEN
         STOP 'IELAS SHOULD BE 0'
         W2      = AMT*(AMT + 2.0*XNU) - Q2
         XS      = Q2/(2.0*AMT*XNU)
         YS      = XNU/(E0*1.0E-3)
         XJACOB  = 2.0D0*PI*AMH*YS/(1.0-YS)
         XJACOB  = 1.0D-3/XJACOB 
         Y_ELAS  = 1.0/(1.0 + AMT/(2.0*E0*1.0E-3*SIN(0.5*TH_RAD)**2))
         W2_ELAS = AMT**2       ! USE HE3 MASS FOR W2
      ELSE
         STOP 'IELAS'
      ENDIF

      IF ( (W2.LE.0.0).OR.(YS.LT.Y_ELAS-1.0E-02) ) THEN
         WRITE(6,*) "PROB IN POLSIG.  W<0 OR Y < Y_ELAS",SQRT(W2_N),W_N,W2
         SIGMA_BORN = 0.0
         SIGMA_RAD  = 0.0
         RETURN
      ENDIF

      SNUC = 2.0*AMH*SQRT((E0*1.0E-3)**2+AML2)       ! 2*MH*E0

      IF (IELAS.EQ.0) THEN
        Y   = SNUC*XS*YS*(AMT/AMH)
        CALL CONKIN_2(SNUC,AMT,IPOL)
      ELSE IF (IELAS.EQ.1) THEN
        Y   = SNUC*XS*YS*(AMT/AMH)
        CALL CONKIN_2(SNUC,AMT,IPOL)
      ENDIF
   

C--------------------------------------------------------------------------
C-----DELTA IS FACTORIZING PART OF VIRTUAL AND REAL LEPTONIC BREMSSTRAHLUNG
C--------------------------------------------------------------------------
      CALL DELTAS_2(TR,FACTOR1,DEL_SUB)

      IF (IPOL.EQ.0.OR.IPOL.EQ.2) THEN  ! CROSS SECTION FOR UNPOLARIZED HADRON TARGET   
         UN   = 1.0                     ! UN = 1 CALCULATES F1 AND F2 STRUCT FUNCT    
         PL   = 1.0                     ! LEPTON POLARIZATION                       
         PN   = 0.0                     ! UNPOLARIZED TARGET, G1 = G2 = 0        
         QN   = 0.0                     ! QN = 0,THE HADRONIC TENSOR FOR SPIN 1/2 PART.
         ISF1 = 1
         ISF2 = 2
         ISF3 = 1
      ELSEIF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN ! DIFF. BETWEEN 2 HADRON POLAR. DIRECTIONS 
         UN   = 0.0                        ! UN = 0 MEANS F1 = F2 = 0         
         PL   = 1.0                        ! LEPTON POLARIZATION                
         PN   = 1.0                        ! PN DEFINES HADRON POL. G1&G2 NON-ZERO
         QN   = 0.0
         ISF1 = 3
         ISF2 = 4
         ISF3 = 1
      ENDIF
           
      EXTAI2     = 1.0D0
      
      CALL BORNIN_2(SIB)                                ! SIB IS DSIG/(DXDY) IN NB

      EXTAI2     = ((SX-Y/TARA)**2/S/(S-Y/TARA))**TR    !USED ONLY IN ELASTIC CALC.
      
CDEBUG

      !TAIL       = AN*ALFA/PI*TAIL_INTEG_2(TAMIN,TAMAX)   
      !SIGMA_BORN = SIB*XJACOB                           ! NB/MEV-SR
      !SIGMA_RAD  = ( SIB*FACTOR1*(1.+ALFA/PI*DEL_SUB) + TAIL )*XJACOB
      !IF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN
      !  SIGMA_BORN = -SIGMA_BORN
      !  SIGMA_RAD  = -SIGMA_RAD
      !ENDIF





      TAIL        = AN*ALFA/PI*TAIL_INTEG_2(TAMIN,TAMAX) * XJACOB
      
      !SIGMA_BORN  = SIB*XJACOB                           ! NB/MEV-SR     ! HERE
      !WRITE(64,*) EP,SIB*XJACOB,SIGMA_BORN
      SIGMA_RAD   = SIGMA_BORN*FACTOR1*(1.+ALFA/PI*DEL_SUB) + TAIL
      CFACT_POL   = (1.+ALFA/PI*DEL_SUB) *FACTOR1
      IF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN
        SIGMA_BORN = -SIGMA_BORN
        SIGMA_RAD  = -SIGMA_RAD
        TAIL_POL   = -TAIL
      ENDIF

      RETURN
      END

C$$$  EXTAI1     = EXP(ALFA/PI*DELINF)
C$$$  EXTAI2     = 1.0                                ! TO CHECK PROTON ELASTIC TAIL
C$$$  SIGMA_RAD  = (SIB*EXTAI1*(1.+ALFA/PI*(DELTA-DELINF))+TAIL)*XJACOB
CSL   EXTAI3 = 1.0D0
CSL   EXTAI3 = ((SX-Y)**2/S/(S-Y))**TR                ! NOT USED



C
C----------------------------------------------------------------------
C
      SUBROUTINE CONKIN_2(SNUC,AMTAR,IPOL)
C     SET OF KINEMATICAL CONSTANTS
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/POL_2/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
*
*     POLARIZATION SPECIFICATION
*
*       IPOL = 0 : PARALLEL, UNPOLARIZED
*       IPOL = 1 : PARALLEL, POLARIZED
*       IPOL = 2 : PERPENDICULAR, UNPOLARIZED
*       IPOL = 3 : PERPENDICULAR, POLARIZED
*

      AMP  = AMTAR        !        M_TARGET
      AP   = 2.*AMP       !      2*M_TARGET
      AMP2 = AMP**2       !       (M_TARGET)**2
      AP2  = 2.*AMP**2    !  2 .0*(M_TARGET)**2
      S    = SNUC*AMP/AMH ! S = S*(M_TARGET/M_HADRON) = 2*M_TARGET*E0
      X    = S*(1.-YS)    ! X = (1-Y)S                = 2*M_TARGET*EP
      SX   = S-X          ! SX=                       = 2*M_TARGET*NU
      SXP  = S+X          ! SP=                       = 2*M_TARGET*(E0+EP)
      YM   = Y+AL2        ! Q^2 + 2M^2 = Q_M^2        ~ Q^2

      TPL  = S**2+X**2    ! 4*M_TARGET**2(E0**2 + EP**2)
      TMI  = S**2-X**2    ! 4*M_TARGET**2(E0**2 - EP**2)

                                              ! W2 IS ALREADY DEFINED IN MAIN PROGRAM
CSL   W2   = AMP2+S-Y-X                       ! USE PROTON MASS FOR W2
C$$$  W2   = AMT*(AMT + (S-X)/AMH) - Y        ! USE HE3 MASS FOR W2

      ALS  = S*S-AL2*AP2                      ! LAMBDA_S = S^2 - 2M^2*M^2 ~ S^2 
      ALX  = X*X-AL2*AP2                      ! LAMBDA_X = X^2 - 2M^2*M^2 ~ X^2
      ALM  = Y*Y+4.*AML2*Y                    ! LAMBDA_Y = Y^2 + 4M^2*Y   ~ Q2^2
      ALY  = SX**2+4.*AMP2*Y                  ! LAMBDA_Q = SX**2 + 4M_TARGET^2 * Q2
                                              !          = 4*M_TARGET^2 * |QVEC|^2

      SQLS = DSQRT(ALS)                       ! SQRT(LAMBDA_S) ~ S
      SQLX = DSQRT(ALX)                       ! SQRT(LAMBDA_X) ~ X
      SQLY = DSQRT(ALY)                       ! SQRT(LAMBDA_Q) = 2*M_TARGET * |QVEC| 
      SQLM = DSQRT(ALM)                       ! SQRT(LAMBDA_Y) ~ Q2
 
      ALLM = DLOG((SQLM+Y)/(SQLM-Y))/SQLM     ! ~ LOG( Q2/M^2 + 1) / Q2
      AXY  = PI*(S-X)                         ! 2*M_TARGET*NU
      AN   = 2.*ALFA**2/SQLS*AXY*BARN*AMH/AMP ! 2*ALPH^2*PI*(NU/E0)(M_P/M_T)*BARN

C     TAMIN= (SX-SQLY)/AP2                    ! EQUIVALENT TO EXPRESSION BELOW.
      TAMAX= (SX+SQLY)/AP2                    ! TAU_MAX
      TAMIN= -Y/AMP2/TAMAX                    ! TAU_MIN

      AS   = S/2./AML/SQLS                    ! ~1/M_E
      BS   = 0.
      CS   = -AML/SQLS                        ! -M_E/(2*M_T*E0)

      IF (IPOL/2.EQ.0) THEN                   ! PARALLEL CONFIGURATION
         AE  = AMP/SQLS                        ! M_T/(2*M_T*E0) = 1/(2*E0)
         BE  = 0.                              !
         CE  = -S/AP/SQLS                      ! ~ 1/(2*M_T)
      ELSE                                    ! PERPENDICULAR CONFIGURATION
         SQN = DSQRT(S*X*Y-ALY*AML2-AMP2*Y*Y)  ! ~SQRT(4M_T^2*E0*EP*Q2 - M_T^2*Q2*Q2)
         AE  = (-S*X+AP2*YM)/SQLS/SQN/2.       ! [1/GEV]
         BE  = SQLS/SQN/2.                     ! [1/GEV]
         CE  = -(S*Y+AL2*SX)/SQLS/SQN/2.       ! [1/GEV]
      ENDIF

      APQ   = -Y*(AE-BE)+CE*SX                                 ! Q*ETA        [GEV]
      APN   = (Y+4.*AML2)*(AE+BE)+CE*SXP                       ! (K1+K2)*ETA  [GEV]
      DK2KS = AS*YM+AL2*BS+CS*X                                ! K_2*KSI      [GEV]
      DKSP1 = AS*S+BS*X+CS*AP2                                 ! KSI*P        [GEV]
      DAPKS = 2.*(AL2*(AS*AE+BS*BE)+AP2*CS*CE+YM*(AS*BE+BS*AE)
     .        +S*(AS*CE+CS*AE)+X*(BS*CE+CS*BE))                ! KSI*ETA      [ 1 ]
      RETURN
      END

CDECK  ID>, BORNIN. 
****************** BORNIN *************************************

      SUBROUTINE BORNIN_2(SIBOR)
C
C     SIBOR IS BORN CROSS SECTION WITH POLARIZED INITIAL
C     LEPTON AND POLARIZED TARGET
C     SIAMM IS CONTRIBUTION OF ANOMALOUS MAGNETIC MOMENT.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/POL_2/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
      COMMON/TAIL_2/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH
      COMMON/PRINT_2/IPRI1
      DIMENSION SFM0(8),TM(8),SFM(8)
*
*-----FIRST DETERMINE REACTION REGION (ELASTIC, QUASI-ELASTIC OR INELASTIC)
*
      IPRI1 = 1
      CALL STRF_2(0D0,0D0,SFM0,SFM)
      IPRI1 = 0

C-----ONLY TM(3) AND TM(4) USED FOR INELASTIC.

      TM(1) = -(2.*AML2-Y)
      TM(2) = (-(AMP2*Y-S*X))/(2.*AMP2)
      TM(3) = (2.*(APQ*DK2KS-DAPKS*Y)*AML)/AMP           ! [GEV^2]
      TM(4) = APQ/AMP*(-(DK2KS*SX-2.*DKSP1*Y)*AML)/AMP2  ! [GEV^2]
      TM(7) = (-(4.*AML2+3.*APN**2-3.*APQ**2+Y))/2.
      TM(8) = APQ/AMP*(-3.*(APN*SXP-APQ*SX))/(2.*AP)
      EK    = (3.*APQ**2-Y)/AMP2
      TM(5) = -EK*TM(1)
      TM(6) = -EK*TM(2)
      SSUM  = 0.

      DO 1 ISF=ISF1,ISF2,ISF3
         PPOL = 1.

         IF(ISF.EQ.3.OR.ISF.EQ.4) PPOL = -PN
         IF(ISF.GE.5)             PPOL = QN/6

         SSUM = SSUM+TM(ISF)*SFM0(ISF)*PPOL               ! [GEV^2]
    1 CONTINUE

      SIBOR = SSUM*2.0*AN/Y**2.                           ! [GEV^2]*[GEV^2 NB]/[GEV^4]      
      RETURN                                              ! = [NB]
      END

CDECK  ID>, DELTAS. 
****************** DELTAS *************************************
      SUBROUTINE DELTAS_2(TR,FACTOR1,DEL_SUB)
C     SUBROUTINE DELTAS(DELTA,DELINF,TR,FACTOR1,DEL_SUB)
C
C-----DELTA IS FACTORIZING PART OF VIRTUAL AND REAL LEPTONIC BREMSSTRAHLUNG
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/CHOI_2/W2_ELAS
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)

      DEL1   = -YM*(ALM*ALLM**2/2.+2.*FSPEN_2(2D0*SQLM/(Y+SQLM))-PI2/2.)/SQLM
      DEL2   = (3.*Y/2.+4.*AML2)*ALLM-2.

      SUM    = VACPOL_2(Y)

      AJ0    = 2.*(YM*ALLM-1.)                               ! =2.*ln(Q2/me^2)
      DELTAI = AJ0*DLOG(DABS(W2-W2_ELAS)/AML/DSQRT(W2))

      SS     = X+Y
      XX     = S-Y
      ALSS   = SS**2-2.*W2*AL2
      ALXX   = XX**2-2.*W2*AL2
      SQLSS  = DSQRT(ALSS)
      SQLXX  = DSQRT(ALXX)
      ALLSS  = DLOG((SQLSS+SS)/(-SQLSS+SS))/SQLSS
      ALLXX  = DLOG((SQLXX+XX)/(-SQLXX+XX))/SQLXX
      DLM    = DLOG(Y/AML2)
      SFPR   = DLM**2/2.-DLM*DLOG(SS*XX/(AML2*W2))
     &          -(DLOG(SS/XX))**2/2.+FSPEN_2((S*X-Y*AMP2)/(SS*XX))-PI2/3.


      DELTA0 = (SS*ALLSS+XX*ALLXX)/2.+SFPR
      
      DELTA   = DELTAI+DELTA0+DEL1+DEL2+SUM
C$$$  DELINF  = (DLM-1.)*DLOG((W2-AMC2)**2/(SS*XX))
      DELINF  = (DLM-1.)*DLOG((W2-W2_ELAS)**2/(SS*XX))
      TR      = ALFA/PI*(DLM-1.)

C$$$  FACTOR1 = ((W2-W2_ELAS)**2/(SS*XX))**(ALFA/PI*0.5*AJ0)
      FACTOR1 = ((W2-W2_ELAS)**2/(SS*XX))**TR

      DEL_SUB = AJ0*DLOG(DSQRT(SS*XX)/AML/DSQRT(W2))
     &          + DELTA0 + DEL1 + DEL2 + SUM

      RETURN
      END
C
C----------------------------------------------------------------------
C
CDECK  ID>, VACPOL. 
      REAL*8 FUNCTION VACPOL_2(T)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
      DIMENSION AM2(3)
C
C     AM2 : SQUARED MASSES OF CHARGE LEPTONS
C
      DATA AM2/.26110D-6,.111637D-1,3.18301D0/

      SUML=0.
      DO 10 I=1,3
         A2=2.*AM2(I)
         SQLMI=DSQRT(T*T+2.*A2*T)
         ALLMI=DLOG((SQLMI+T)/(SQLMI-T))/SQLMI
  10  SUML=SUML+2.*(T+A2)*ALLMI/3.-10./9.+4.*A2*(1.-A2*ALLMI)/3./T
      IF(T.LT.1.D0)THEN
        AAA = -1.345D-9
        BBB = -2.302D-3
        CCC = 4.091
      ELSEIF(T.LT.64D0)THEN
        AAA = -1.512D-3
        BBB =  -2.822D-3
        CCC = 1.218
      ELSE
        AAA = -1.1344D-3
        BBB = -3.0680D-3
        CCC = 9.9992D-1
      ENDIF
      SUMH = -(AAA+BBB*LOG(1.+CCC*T)) *2*PI/ALFA

      VACPOL_2=SUML+SUMH

      END

CDECK  ID>, FSPENS. 
****************** FSPENS *************************************

      REAL*8 FUNCTION FSPENS_2(X)
C
C    SPENCE FUNCTION
C
      IMPLICIT REAL*8(A-H,O-Z)
      F=0.D0
      A=1.D0
      AN=0.D0
      TCH=1.D-16
  1   AN=AN+1.D0
      A=A*X
      B=A/AN**2
      F=F+B
      IF(B-TCH)2,2,1
  2   FSPENS_2=F
      RETURN
      END
CDECK  ID>, FSPEN.  
****************** FSPEN **************************************

      REAL*8 FUNCTION FSPEN_2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA F1/1.644934D0/
      IF(X)8,1,1
  1   IF(X-.5D0)2,2,3
    2 FSPEN_2=FSPENS_2(X)
      RETURN
    3 IF(X-1D0)4,4,5
    4 FSPEN_2=F1-DLOG(X)*DLOG(1D0-X+1D-10)-FSPENS_2(1D0-X)
      RETURN
    5 IF(X-2D0)6,6,7
    6 FSPEN_2=F1-.5*DLOG(X)*DLOG((X-1D0)**2/X)+FSPENS_2(1D0-1D0/X)
      RETURN
    7 FSPEN_2=2D0*F1-.5D0*DLOG(X)**2-FSPENS_2(1D0/X)
      RETURN
    8 IF(X+1D0)10,9,9
   9  FSPEN_2=-.5D0*DLOG(1D0-X)**2-FSPENS_2(X/(X-1D0))
      RETURN
  10  FSPEN_2=-.5*DLOG(1.-X)*DLOG(X**2/(1D0-X))-F1+FSPENS_2(1D0/(1D0-X))
      RETURN
      END


CDECK  ID>, TAILS.  
****************** TAILS **************************************

       SUBROUTINE TAILS_2(TA,TM)
       IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/POL_2/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
       COMMON/BSEO_2/OIS,OIR,OI12,EEIS,EEIR,EEI12,
     . EEI1I2,EB,EEB,TM3(6,4,3)
       DIMENSION TM(8,6),AJM2(2),AJM3(3),II(8)
      DATA II/1,2,3,4,1,2,5,6/

      B2=(-ALY*TA+SXP*SX*TA+2.*SXP*Y)/2. ! B_2
      B1=(-ALY*TA-SXP*SX*TA-2.*SXP*Y)/2. ! B_1
      C1=-(4.*(AMP2*TA**2-SX*TA-Y)*AML2-(S*TA+Y)**2) ! C_1
      C2=-(4.*(AMP2*TA**2-SX*TA-Y)*AML2-(TA*X-Y)**2) ! C_2
      BB=1./SQLY                ! F
      SC1=DSQRT(C1)
      SC2=DSQRT(C2)
      BI12=(SXP*(SX*TA+2.*Y))/(SC1*SC2*(SC1+SC2)) ! F_D
      BI1PI2=1./SC2+1./SC1      ! F_1+
      BIS=-B1/SC1/C1+B2/SC2/C2  ! F_2+
      BIR=B2/SC2/C2+B1/SC1/C1   ! F_2-
      B1I=-B1/ALY/SQLY          ! F_I
      B11I=(3.*B1**2-ALY*C1)/2./ALY**2/SQLY ! F_II
      SPS=AS+BS                 ! S_KSI
      SPE=AE+BE                 ! S_ETA
      CCPE=(AE-BE)*TA+2.*CE     ! R_ETA
      CCPS=(AS-BS)*TA+2.*CS     ! R_KSI

      SIS=(2.*BI1PI2*SPS+BIR*SPS*TA+BIS*CCPS)/2. ! F_{2+}^\KSI
      SIR=( (2.*BI12*SPS*TA+BIR*CCPS+BIS*SPS*TA))/2. ! F_{2-}^\KSI
      SI12=(BI12*CCPS+BI1PI2*SPS)/2. ! F_D^\KSI
      EIS=(2.*BI1PI2*SPE+BIR*SPE*TA+BIS*CCPE)/2. ! F_{2+}^\ETA
      EIR=( (2.*BI12*SPE*TA+BIR*CCPE+BIS*SPE*TA))/2. ! F_{2-}^\ETA
      EI12=(BI12*CCPE+BI1PI2*SPE)/2. ! F_D^\ETA

      OIS=((2.*BI1PI2+BIR*TA)*(CCPE*SPS+CCPS*SPE)+(CCPE*CCPS+
     . SPE*SPS*TA**2)*BIS+8.*BB*SPE*SPS+4.*BI12*SPE*SPS*TA**2)/
     .     4.                   ! F_{2+}^{\KSI\ETA}
      OIR=( ((2.*BI12+BIS)*(CCPE*SPS+CCPS*SPE)*TA+(CCPE*CCPS+
     . SPE*SPS*TA**2)*BIR+4.*BI1PI2*SPE*SPS*TA))/4. ! F_{2-}^{\KSI\ETA}
      OI12=((CCPE*CCPS+SPE*SPS*TA**2)*BI12+(CCPE*SPS+CCPS*SPE)*
     . BI1PI2+4.*BB*SPE*SPS)/4. ! F_D^{\KSI\ETA}
      EEIS=((CCPE**2+SPE**2*TA**2)*BIS+8.*BB*SPE**2+4.*BI12*SPE
     . **2*TA**2+4.*BI1PI2*CCPE*SPE+2.*BIR*CCPE*SPE*TA)/4. ! F_{1+}^{\ETA\ETA}
      EEIR=( ((CCPE**2+SPE**2*TA**2)*BIR+4.*BI12*CCPE*SPE*TA+4.
     . *BI1PI2*SPE**2*TA+2.*BIS*CCPE*SPE*TA))/4.
      EEI12=((CCPE**2+SPE**2*TA**2)*BI12+4.*BB*SPE**2+2.*BI1PI2
     . *CCPE*SPE)/4.
      EI1PI2=(4.*BB*SPE+BI12*SPE*TA**2+BI1PI2*CCPE)/2.
      EEI1I2=((CCPE**2+SPE**2*TA**2)*BI1PI2+4.*(2.*CCPE-SPE*TA)
     . *BB*SPE+8.*B1I*SPE**2+2.*BI12*CCPE*SPE*TA**2)/4.
      EB=((CCPE-SPE*TA)*BB+2.*B1I*SPE)/2.
      EEB=((CCPE-SPE*TA)**2*BB+4.*(CCPE-SPE*TA)*B1I*SPE+4.*B11I
     . *SPE**2)/4.
       CALL FFU_2(1,BB,BIS,BIR,BI12,BI1PI2,SIR,SIS,SI12
     .,EIS,EIR,EI12,EI1PI2,TA)
       CALL FFU_2(2,EB,EIS,EIR,EI12,EI1PI2,OIR,OIS,OI12
     .,EEIS,EEIR,EEI12,EEI1I2,TA)
       CALL FFU_2(3,EEB,EEIS,EEIR,EEI12,EEI1I2,0D0,0D0,0D0
     .,0D0,0D0,0D0,0D0,TA)
       AJM2(1)=APQ/AMP
       AJM2(2)=-1./AMP
       AJM3(1)=(Y-3.*APQ**2)/AMP2
       AJM3(2)=6.*APQ/AMP2
       AJM3(3)=-3./AMP2
       DO 15 I=1,8
       DO 13 L=1,6
   13  TM(I,L)=0
       DO 10 K=1,I2(I)
       AJK=1.
       IF(I.EQ.4.OR.I.EQ.8)AJK=AJM2(K)
       IF(I.EQ.5.OR.I.EQ.6)AJK=AJM3(K)
       DO 10 J=K,I1(I)+K-1
       TM(I,J)=TM(I,J)+TM3(II(I),J-K+1,K)*AJK
       IF((I.EQ.5.OR.I.EQ.6).AND.K.EQ.2)
     . TM(I,J)=TM(I,J)+TM3(II(I),J-K+1,1)*TA/AMP2
  10   CONTINUE
  15   CONTINUE
       RETURN
       END

CDECK  ID>, FFU.
****************** FFU ****************************************

       SUBROUTINE FFU_2(N,BB,BIS,BIR,BI12,BI1PI2,SIR,SIS,SI12
     .        ,EIS,EIR,EI12,EI1PI2,TA)
       IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/POL_2/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
       COMMON/BSEO_2/OIS,OIR,OI12,EEIS,EEIR,EEI12,
     . EEI1I2,EB,EEB,TM3(6,4,3)
      HI2=AML2*BIS-YM*BI12
      SHI2=AML2*SIS-YM*SI12
      EHI2=AML2*EIS-YM*EI12
      OHI2=AML2*OIS-YM*OI12
       GOTO(10,20,30)N
  10   CONTINUE
      TM3(3,1,N)=(8.*(APQ*DK2KS-DAPKS*Y)*AML*HI2)/AMP
      TM3(3,2,N)=(-2.*((2.*(BI12*DK2KS*TA-2.*SHI2)*APQ+(2.*SHI2-
     . SIR*Y+SIS*YM)*APN+4.*DAPKS*HI2*TA)-4.*((2.*EI12-EIS)*
     . DK2KS-(SI12-SIS)*APN)*AML2)*AML)/AMP
      TM3(3,3,N)=(2.*(((2.*SI12+SIR-SIS)*APN*TA-2.*DK2KS*EI12*TA
     . -6.*OHI2-OIR*Y+OIS*YM)-4.*AML2*OI12)*AML)/AMP
      TM3(3,4,N)=(2.*(2.*OI12-OIR+OIS)*AML*TA)/AMP
      TM3(5,1,N)=-2.*(4.*AML2+3.*APN**2-3.*APQ**2+Y)*HI2
      TM3(5,2,N)=-2.*(6.*AML2*APN*EIR-3.*APN**2*BI12*TA+3.*APN*
     . APQ*BI1PI2+6.*APQ*EHI2+HI2*TA)
      TM3(5,3,N)=-(24.*AML2*EEI12-6.*APN*EI1PI2-6.*APQ*EI12*TA-
     . 2.*BB-BI12*TA**2)
  20   CONTINUE
      TM3(4,1,N)=(-4.*(DK2KS*SX-2.*DKSP1*Y)*AML*HI2)/AMP2
      TM3(4,2,N)=(((2.*(SXP-2.*SX)*SHI2+2.*BI12*DK2KS*SX*TA+8.*
     . DKSP1*HI2*TA-SIR*SXP*Y+SIS*SXP*YM)-4.*(2.*BI12*DK2KS-BIS*
     . DK2KS-SI12*SXP+SIS*SXP)*AML2)*AML)/AMP2
      TM3(4,3,N)=((((SXP*TA-YM)*SIS-(SXP*TA-Y)*SIR+2.*BI12*DK2KS
     . *TA+6.*SHI2-2.*SI12*SXP*TA)+4.*AML2*SI12)*AML)/AMP2
      TM3(4,4,N)=(-(2.*SI12-SIR+SIS)*AML*TA)/AMP2
      TM3(6,1,N)=(-3.*(APN*SXP-APQ*SX)*HI2)/AMP
      TM3(6,2,N)=(-3.*(2.*(APN*BIR+EIR*SXP)*AML2-(2.*BI12*SXP*TA
     . -BI1PI2*SX)*APN+(BI1PI2*SXP+2.*HI2)*APQ+2.*EHI2*SX))/(2.*
     . AMP)
      TM3(6,3,N)=(-3.*(8.*AML2*EI12-APN*BI1PI2-APQ*BI12*TA-EI12*
     . SX*TA-EI1PI2*SXP))/(2.*AMP)
  30   CONTINUE
      TM3(1,1,N)=-4.*(2.*AML2-Y)*HI2
      TM3(1,2,N)=4.*HI2*TA
      TM3(1,3,N)=-2.*(2.*BB+BI12*TA**2)
      TM3(2,1,N)=(((SXP**2-SX**2)-4.*AMP2*Y)*HI2)/(2.*AMP2)
      TM3(2,2,N)=(2.*AML2*BIR*SXP-4.*AMP2*HI2*TA-BI12*SXP**2*TA+
     . BI1PI2*SXP*SX+2.*HI2*SX)/(2.*AMP2)
      TM3(2,3,N)=(2.*(2.*BB+BI12*TA**2)*AMP2+4.*AML2*BI12-BI12*
     . SX*TA-BI1PI2*SXP)/(2.*AMP2)
       RETURN
       END


CDECK  ID>, STRF.   
      SUBROUTINE STRF_2(TA,RR,SFM,SFM0)
C
C     THE PROGRAMM CALCULATES DEEP INELASTIC (ITA=1),
C     ELASTIC (ITA=2), QUASIELASTIC (ITA=3) STRUCTURE FUNCTIONS
C     IN KINEMATICAL POINT (TA,RR).
C          RR=SX-TT,
C          TA=(T-Y)/RR,
C     WHERE TT=T+AMF2-AMP2, AMF2 IS INVARINT MASS OF FINAL HADRONS
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/SOFTPHOTON_2/EXTAI2
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/TAIL_2/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH
      COMMON/PRINT_2/IPRI1
      DIMENSION SFM(8),SFM0(8) !,FIGI(4,2)
      COMMON/ELASTIC_2/IELAS
      !REAL*8 F1_FROM_DATA,F2_FROM_DATA,G1_FROM_DATA,G2_FROM_DATA
      COMMON/DEBUG_2/IDEBUG1

      !INTEGER ICONF,ITYPE,ASYMCHOICE
      !COMMON /CONFIG4/ICONF,ITYPE,ASYMCHOICE

*
*     INITIALIZE FORM FACTORS
*
      DO I = 1,8
         SFM(I)  = 0.0D0
         SFM0(I) = 0.0D0
      ENDDO
       
      T    = Y+RR*TA              ! SAME AS Q**2 WHEN THERE IS NO REAL PHOTON EMITTED
      TT   = SX-RR                ! SAME AS 2*M*NU WHEN THERE IS NO REAL PHOTON
      AMF2 = TT-T+AMP2            ! SAME AS W**2 WHEN THERE IS NO REAL PHOTON
      AKS  = T/TT                 ! SAME AS X_BJORKEN WHEN THERE IS NO REAL PHOTON
      ANU  = TT/AP                ! SAME AS NU WHEN THERE IS NO REAL PHOTON
      B1   = 0.D0
      B2   = 0.D0
      B3   = 0.D0
      B4   = 0.D0
*
*-----IF T<0, THEN INELASTIC (OR QUASI-ELASTIC) CHANNEL IS NOT POSSIBLE
*     
      IF (T.LE.0.0)    RETURN
      IF (AMF2.LT.0.0) RETURN
      IF (IELAS.GT.0)  GOTO 20

      EPSI = AP2/SX                  ! =2M^2/(2M NU) = M/NU
      AKS0 = Y/SX                    ! =Q2/2MNU      = XBJORKEN 
      F1 = 0.0 !F1_FROM_DATA(AKS0,Y)      ! SAME AS F1SFUN(X,Q^2)
      F2 = 0.0 !F2_FROM_DATA(AKS0,Y)      ! SAME AS F2SFUN(X,Q^2)
      !G1 = G1_FROM_DATA(AKS0,Y)
      !G2 = G2_FROM_DATA(AKS0,Y)

      CALL GET_G(AKS0,Y,G1,G2)
      !G1=-G1 !HERE
      !G2=-G2 
      ! AKS0,Y == X,Q2 AT THE KINEMATIC POINTS OF THE DATA.
      !WRITE(46,*) AKS0,Y,G1,G2


      SFM0(1) = UN*F1+QN/6.*B1
      SFM0(2) = EPSI    * (UN*F2+QN/6.*B2)

      SFM0(3) = EPSI    * (G1+G2)
      SFM0(4) = EPSI**2 * G2

      SFM0(5) = EPSI**2 * B1
      SFM0(6) = EPSI**3 * (B2/3.+B3+B4)
      SFM0(7) = EPSI    * (B2/3.-B3)
      SFM0(8) = EPSI**2 * (B2/3.-B4)

      IF (TA.EQ.0.0D0.AND.RR.EQ.0.0D0) THEN
         DO I = 1,8
            SFM(I) = SFM0(I)
         ENDDO
      ELSE ! CONTRIBUTION FROM THE INELASTIC CHANNEL
         EPSI = AP2/TT            ! SAME AS M/NU WHEN THERE IS NO REAL PHOTON

         F1 = 0.0 !F1_FROM_DATA(AKS,T) ! SAME AS F1SFUN(X,Q^2)
         F2 = 0.0 ! F2_FROM_DATA(AKS,T) ! SAME AS F2SFUN(X,Q^2)
         !G1 = G1_FROM_DATA(AKS,T)
         !G2 = G2_FROM_DATA(AKS,T)
         !WRITE(62,*) AKS,T*1.E6
         CALL GET_G(AKS,T,G1,G2)
         !G1=-G1 !HERE
         !G2=-G2
         !WRITE(46,*) AKS,G1,G2
 
         SFM(1) = UN*F1+QN/6.*B1
         SFM(2) = EPSI    * (UN*F2+QN/6.*B2)

         SFM(3) = EPSI    * (G1+G2)
         SFM(4) = EPSI**2 * G2

         SFM(5) = EPSI**2 * B1
         SFM(6) = EPSI**3 * (B2/3.+B3+B4)
         SFM(7) = EPSI    * (B2/3.-B3)
         SFM(8) = EPSI**2 * (B2/3.-B4)
      ENDIF
 10   CONTINUE
      IF (IELAS.LE.0) RETURN
 20   CONTINUE

*
*-----CONTRIBUTION FROM THE ELASTIC CHANNEL
*
      WRITE(6,*) "ELASTIC CONTRIBUTION"
      EPSI = AP2/T*(AMT/AMP)**2
      TAU  = T/4./AMT**2
      TAU1 = 1.+TAU

      CALL FFHE3_2(T,GE,GM)

C$$$  GE = (1.0D0 + T/0.71D0)**-2 ! PROTON FORM FACTOR USED BY MO AND TSAI
C$$$  GM = 2.793*GE

      XNU_EL = 0.5*T/AMT

      IF (ABS( ANU-XNU_EL ).GT.1.0D-3) THEN
         SPREAD = 0.0
      ELSE
         SPREAD = 1.0
      ENDIF
      
      F1 =     AMT * TAU    *            GM**2       * SPREAD
      F2 = 2.0*AMT * TAU    * (GE**2+TAU*GM**2)/TAU1 * SPREAD
      G1 =     AMT * TAU    *    GM*(GE+TAU*GM)/TAU1 * SPREAD
      G2 =     AMT * TAU**2 *    GM*(GE-GM)    /TAU1 * SPREAD

      FACTOR  = AMP/AMT
      FACTOR2 = FACTOR**2
      FACTOR3 = FACTOR2*FACTOR
C-----
      SFM(1) = SFM(1) + EXTAI2 * FACTOR            * (UN*F1+QN/6.*B1)
      SFM(2) = SFM(2) + EXTAI2 * FACTOR  * EPSI    * (UN*F2+QN/6.*B2)

      SFM(3) = SFM(3) + EXTAI2 * FACTOR2 * EPSI    * (G1+G2)
      SFM(4) = SFM(4) + EXTAI2 * FACTOR3 * EPSI**2 * G2

      SFM(5) = SFM(5) + EXTAI2 * FACTOR3 * EPSI**2 * B1
      SFM(6) = SFM(6) + EXTAI2 * FACTOR3 * EPSI**3 * (B2/3.+B3+B4)
      SFM(7) = SFM(7) + EXTAI2 * FACTOR  * EPSI    * (B2/3.-B3)
      SFM(8) = SFM(8) + EXTAI2 * FACTOR2 * EPSI**2 * (B2/3.-B4)
C-----
 30   CONTINUE
      RETURN
      END


CDECK  ID>, FFHE3.  
****************** FFHE3 **************************************

      SUBROUTINE FFHE3_2(T,GE,GM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      TF=T/CHBAR**2
      IF (TF.GE.0.0) THEN
         QF=SQRT(TF)
      ELSE
         QF=-SQRT(-TF)
      ENDIF
      A=.675
      B=.366
      C=.836
      AM=.654
      BM=.456
      CM=.821
      D=-6.78D-3
      P=.9
      Q0=3.98
      F0=DDEXP_2(-A**2*TF) - B**2*QF**2*DDEXP_2(-C**2*TF)
      FM=DDEXP_2(-AM**2*TF) - BM**2*QF**2*DDEXP_2(-CM**2*TF)
      DF=D*DDEXP_2(-((QF-Q0)/P)**2)
      GE=(F0+DF)*TARZ
      GM=FM*TARA      * (-2.13)
      END



CDECK  ID>, DDEXP.  
****************** DDEXP **************************************

      REAL*8 FUNCTION DDEXP_2(X)
      IMPLICIT REAL*8(A-H,O-Z)
        DDEXP_2=0.
        IF(X.GT.-50.)DDEXP_2=EXP(X)
      RETURN
      END





C
C----------------------------------------------------------------------
C
      REAL*8 FUNCTION DELTA_FTN_2(X,XM,SIG)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 X,XM,SIG,PI,SQ2PI
      !DATA PI/3.1415926D0/
      DATA SQ2PI/2.506628275D0/

      DELTA_FTN_2 = DEXP(-0.5D0*((X-XM)/SIG)**2)/(SIG*SQ2PI)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      
      REAL*8 FUNCTION TAIL_INTEGRAND_2(R)
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 SFM0(8),SFM(8),TM(8,6)
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/TAIL_2/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH
      COMMON/P_2/PI,PI2,ALFA,I1(8),I2(8)
      COMMON/TAIL_INTEGRAL_2/TAU,TM
      CALL STRF_2(TAU,R,SFM,SFM0)
      SUM = 0.0D0
      DO ISF = ISF1,ISF2,ISF3
         IF (ISF.EQ.3.OR.ISF.EQ.4) THEN
            PPOL = -PN
         ELSE IF (ISF.EQ.5) THEN
            PPOL = QN/6.0D0
         ELSE
            PPOL = 1.0D0
         ENDIF
         DO IRR = 1,I1(ISF)+I2(ISF)-1
            IF (IRR.EQ.1) THEN
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*(SFM(ISF)/(Y+
     &              R*TAU)**2-SFM0(ISF)/Y**2)
            ELSE
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*SFM(ISF)/(Y+
     &              R*TAU)**2
            ENDIF
            SUM = SUM + TEMP
         ENDDO
      ENDDO
      TAIL_INTEGRAND_2 = SUM

      !WRITE(68,'(5F15.3)') TAU,R,TAIL_INTEGRAND_2,TAMIN,TAMAX
      RETURN
      END
C
C----------------------------------------------------------------------
C

      REAL*8 FUNCTION TAIL_INT_OVER_R_2(TALN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 TALN,TAU,TAIL_INTEGRAND_2,TM(8,6)
      EXTERNAL TAIL_INTEGRAND_2
      REAL*8 R_EL,RCUT(100),SUM,TEMP
      REAL*8 TAU_PASS,DGQUAD
      COMMON/TAIL_INTEGRAL_2/TAU_PASS,TM
      REAL*8 AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,AMT,TARA,
     &     TARZ,FERMOM,AMM,AMN,CHBAR,BARN,S,X,SX,SXP,Y,YM,W2,ALS,
     &     ALX,ALM,ALY,SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,
     &     TMI
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     &     SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/ELASTIC_2/IELAS

      ! C030405
      COMMON/XGOOF_2/W2_START 

      DATA NENTRY/0/
      NENTRY = NENTRY+1

      TAU = DEXP(TALN)-XS
      !WRITE(68,*) W2,TALN
      CALL TAILS_2(TAU,TM)
      TAU_PASS = TAU       ! PASS TAU ARGUMENT TO TAIL_INTEGRAND_2

      IF (IELAS.EQ.0) THEN
         NINT = 1
         RCUT(1) = 1.0D-10

C$$$         W2      = AMP2+S-Y-X                        ! USE PROTON MASS FOR W2
C$$$         RCUT(2) = (W2-AMP2-1D-10)/(1.0D0+TAU)

C030405  RCUT(2) = ( (S-X)-Y -1.0D-10) /(1.0D0+TAU)  ! STARTING FROM ELASTIC

         W2      = AMP2+S-Y-X                        ! W IN TERMS OF TARGET MASS.
         RCUT(2) = (W2-W2_START)/(1.0D0+TAU)         ! STARTS INTEGRAL AT LOWEST INPUT W.
         !WRITE(68,*) W2,W2_START,TAU
         IF (RCUT(2).LT.RCUT(1)) THEN 
            RCUT(2) = RCUT(1)
         ENDIF
         IF (TAU.LT.0.0) THEN
            R_MAX = -Y/TAU
            IF (RCUT(2).GT.R_MAX) RCUT(2) = R_MAX
         ENDIF
         IF (RCUT(2).GT.SX) THEN
            RCUT(2) = SX
         ENDIF
C
C------- INNER INTEGRAL WRT R.
C------- FOR DGQUAD N CAN BE 2-16,20,24,32,40,48,64,80,OR 96
C
         SUM = 0.0D0
         DO I = 1,NINT
            IF (RCUT(I+1).GT.RCUT(I)) THEN
               !WRITE(68,*) SQRT(W2),RCUT(I),RCUT(I+1)
               TEMP = DGQUAD(TAIL_INTEGRAND_2,RCUT(I),RCUT(I+1),96)    !HERE
               !TEMP = DGAUSS(TAIL_INTEGRAND_2,RCUT(I),RCUT(I+1),1.0D-3)
               SUM  = SUM + TEMP
            ENDIF
         ENDDO
      ELSE IF (IELAS.EQ.1) THEN ! ADD CONTRIB. FROM THE QE PEAK ASSUMING DELTA FUNCTION
C$$$     TEMP = TAIL_INTEGRAND_2(R_QE)*2.0*AMP/(1.0+TAU)

*        ADD CONTRIBUTION FROM THE ELASTIC PEAK
*        FIRST IDENTIFY ELASTIC AND QUASIELASTIC R
*        ELASTIC R
*
         R_EL = (SX-Y)/(1.0 + TAU)
         SUM  = TAIL_INTEGRAND_2(R_EL)*2.0*AMT/(1.0+TAU)
      ELSE
         STOP 'IELAS'
      ENDIF
      TAIL_INT_OVER_R_2 = SUM*(XS+TAU)
      RETURN
      END
C
C----------------------------------------------------------------------
C

      REAL*8 FUNCTION TAIL_INTEG_2(TAU_MIN,TAU_MAX)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 TAU_MIN,TAU_MAX
      REAL*8 TAIL_INT_OVER_R_2,TEMP,DGQUAD
      EXTERNAL TAIL_INT_OVER_R_2
      REAL*8 TCUT(100)
      !DIMENSION NGAUSSPT(7)
      !DATA NGAUSSPT/8,16,16,8,16,16,8/
      COMMON/DEBUG_2/IDEBUG1
      COMMON/XGOOF2_2/IMARK,JMARK
      REAL*8 AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,AMT,TARA,
     &     TARZ,FERMOM,AMM,AMN,CHBAR,BARN,S,X,SX,SXP,Y,YM,W2,ALS,
     &     ALX,ALM,ALY,SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,
     &     TMI
      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG_2/AMT,TARA,TARZ
      COMMON/SXY_2/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     &     SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      REAL*8 TAU_PASS,TM(8,6),R,TAIL_INTEGRAND_2
      COMMON/TAIL_INTEGRAL_2/TAU_PASS,TM
      COMMON/ELASTIC_2/IELAS
      !REAL*8 TEMP1,TEMP2,TEMP3,WEIGHT
      IF (IELAS.EQ.0) THEN                  ! QUASI-ELASTIC AND INELASTIC
         TCUT(1) = DLOG(TAU_MIN+XS)
C$$$     TCUT(2) = DLOG(XS)
C$$$     TCUT(3) = DLOG(TAU_MAX+XS)
         TCUT(3) = DLOG(XS-Y/S)              
         TCUT(6) = DLOG(XS+Y/X)              
         DGUESS  = (TCUT(6)-TCUT(3))*1.0D-1
         TCUT(2) = TCUT(3) - DGUESS
         TCUT(4) = TCUT(3) + DGUESS
         TCUT(5) = TCUT(6) - DGUESS
         TCUT(7) = TCUT(6) + DGUESS
         TCUT(8) = DLOG(XS+TAU_MAX)
         ISTART  = 1                    
         IEND    = 7                    

         ! C030405
         ! THIS SHOULD COVER THE CANONICAL RC TRIANGLE FINE.
         ! NOT SURE WHAT THE PREVIOUS CHOICE OF LIMITS WAS FOR.

         TCUT(10) = DLOG(XS-Y/S)
         TCUT(11) = DLOG(XS+Y/X)
         ISTART   = 10                   
         IEND     = 10

      ELSE IF (IELAS.EQ.1) THEN
         ISTART  = 1             
         TCUT(1) = DLOG(TAU_MIN+XS)
         TCUT(2) = DLOG(XS)
         TCUT(3) = DLOG(TAU_MAX+XS)
         IEND = 2
      ELSE
         STOP 'IELAS'
      ENDIF

C      
C---- OUTER INTEGRAL WRT TAU.
C---- FOR DGQUAD N CAN BE 2-16,20,24,32,40,48,64,80,OR 96
C----
C---- DGAUSS AND DGQUAD GIVE ALMOST IDENTICAL RESULTS (FOR E94010 KINEMATICS)
C---- WITH 1.0D-3 PRECISION.  AT 1.0D-6 PRECISION IT IS JUST WAY TOO SLOW.

      TAIL_INTEG_2 = 0.0D0

      DO I = ISTART,IEND
         !WRITE(68,*) I,W2,TCUT(I),TCUT(I+1)
         !print *, I,W2,TCUT(I),TCUT(I+1)
         TEMP       = DGQUAD(TAIL_INT_OVER_R_2,TCUT(I),TCUT(I+1),96)    !HERE HERE
         !TEMP       = DGAUSS(TAIL_INT_OVER_R_2,TCUT(I),TCUT(I+1),1.0D-3)
         TAIL_INTEG_2 = TAIL_INTEG_2 + TEMP
      ENDDO

      RETURN
      END

C
C----------------------------------------------------------------------
C------------------------------NOTUSED---------------------------------
C----------------------------------------------------------------------
C
CNOTUSED
CNOTUSEDCDECK  ID>, F1SFUN. 
CNOTUSED********************** F1SFUN ***********************************
CNOTUSED      REAL*8 FUNCTION F1SFUN_2(AKS,T)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      F2=F2SFUN_2(AKS,T)
CNOTUSED      ANU=T/AP/AKS
CNOTUSED      F1SFUN_2=AMP*(1.+ANU**2/T)*F2/ANU/(1.+R1990_2(AKS,T))
CNOTUSED      END
CNOTUSED
CNOTUSEDCDECK  ID>, F2SFUN. 
CNOTUSED********************** F2SFUN ***********************************
CNOTUSED      REAL*8 FUNCTION F2SFUN_2(AKS,T)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      F2P=DF2H8(T,AKS)
CNOTUSED      F2D=DF2D8_2(T,AKS)
CNOTUSEDC$$$      F2SFUN=(F2P+2D0*F2D)/3D0
CNOTUSED*
CNOTUSED* MODIFIED BY CHOI
CNOTUSED*
CNOTUSED      F2SFUN_2=F2P+2D0*F2D
CNOTUSED      END
CNOTUSED
CNOTUSED
CNOTUSEDCDECK  ID>, R1990.  
CNOTUSED********************** R1990 ************************************
CNOTUSED
CNOTUSED      REAL*8 FUNCTION R1990_2(AKS,TC)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      T=TC
CNOTUSED      IF(TC.LT..35D0)T=0.35
CNOTUSED      TETA=1.+12.*T/(1.+T)*(.125**2/(AKS**2+.125**2))
CNOTUSED      ZN=TETA/LOG(T/.04)
CNOTUSED      RA=.0672*ZN+.4671/(T**4+1.8979**4)**(.25D0)
CNOTUSED      RB=.0635*ZN+.5747/T-.3534/(T**2+.09)
CNOTUSED      RC=.0599*ZN+.5088/SQRT((T-5.*(1.-AKS)**5)**2+2.1081**2)
CNOTUSED      RRR=(RA+RB+RC)/3.
CNOTUSED      R1990_2=RRR
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSED
CNOTUSED
CNOTUSEDCDECK  ID>, G1SFUN. 
CNOTUSED********************** G1SFUN ***********************************
CNOTUSED      REAL*8 FUNCTION G1SFUN_2(AKS,T)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSEDC        A1HE3=-0.01589*AKS**(-.663)*(1D0-DDEXP_2(-5.96*AKS))
CNOTUSED           DZ= 1./2.*LOG(1.+EXP(2.0-1000.*AKS))
CNOTUSED           DFN=1.
CNOTUSED           DF2NF2P=0.67225*(1.0-AKS)**1.6254-0.15436*AKS**0.048301
CNOTUSED     1            +(.41979+.047331*DZ-0.17816*DZ**2)
CNOTUSED           DF=DFN*(1./((2./DF2NF2P)+1))
CNOTUSED
CNOTUSED         A1NUE=0.00024-.00463*(AKS**0.1+AKS**0.5)
CNOTUSED     .          -3.48645*AKS+1.59218*AKS**1.5
CNOTUSED     .          +8.59393*AKS**2-5.74029*AKS**3
CNOTUSED         A1HE3=A1NUE*DF
CNOTUSED        G1SFUN_2=A1HE3*F1SFUN_2(AKS,T)
CNOTUSED
CNOTUSED      END
CNOTUSED
CNOTUSED
CNOTUSEDCDECK  ID>, G2SFUN. 
CNOTUSED********************** G2SFUN ***********************************
CNOTUSED      REAL*8 FUNCTION G2SFUN_2(AKS,T)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      G2SFUN_2 = 0.0D0
CNOTUSED      END
CNOTUSEDCDECK  ID>, FFPRO.  
CNOTUSED****************** FFPRO **************************************
CNOTUSED
CNOTUSED      SUBROUTINE FFPRO_2(T,GEP,GMP)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      GEP=1.2742/(1.+T/0.6394**2)-.2742/(1.+T/1.582**2)
CNOTUSED      GMP=(1.3262/(1.+T/0.6397**2)-.3262/(1.+T/1.3137**2))*AMM
CNOTUSEDC     GEP=1./((1.+.61*T)*(1.+2.31*T)*(1.+.04*T))
CNOTUSEDC     GMP=AMM*GEP
CNOTUSED      END
CNOTUSED
CNOTUSED
CNOTUSEDCDECK  ID>, FFQUAS.
CNOTUSED****************** FFQUAS **************************************
CNOTUSED
CNOTUSED      SUBROUTINE FFQUAS_2(T,GEUN,GMUN,GEPO,GMPO)
CNOTUSED      IMPLICIT REAL*8(A-H,O-Z)
CNOTUSED      COMMON/CMP_2/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
CNOTUSED     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
CNOTUSED      COMMON/POLTARG_2/AMT,TARA,TARZ
CNOTUSED      CALL FFPRO_2(T,GEP,GMP)
CNOTUSED      TF=T/CHBAR**2
CNOTUSED      TAU=T/4./AMP**2
CNOTUSED      TAU1=1.+TAU
CNOTUSEDC
CNOTUSEDC   T. DE FOREST AND D.J. VALECKA ADV. IN PHYS. 15(1966) NO.57
CNOTUSEDC
CNOTUSED        SUPELE=1.
CNOTUSED        SUPMAG=1.
CNOTUSED        IF(T.LE.(2.D0*FERMOM)**2)THEN
CNOTUSED           SQRAT=DSQRT(T)/FERMOM
CNOTUSED           SUPELE=0.75*SQRAT-SQRAT**3/16.
CNOTUSED           SUPMAG=SUPELE
CNOTUSED        ENDIF
CNOTUSED        GEUN=GEP*DSQRT(SUPELE*TARZ)
CNOTUSED        TARN=TARA-TARZ
CNOTUSED        GMUN=GEP*DSQRT(SUPMAG*(TARZ*AMM**2+TARN*AMN**2))
CNOTUSED        GEPO=0.
CNOTUSED        TARN=TARA-TARZ
CNOTUSED        GMPO=GEP*DSQRT(SUPMAG*(TARN*AMN**2))
CNOTUSED        END
CNOTUSED
CNOTUSEDCDECK  ID>, DF2D8.  
CNOTUSED      REAL*8 FUNCTION DF2D8_2(DQ2,DX)
CNOTUSED*:=====================================================================:
CNOTUSED*:                                                                     :
CNOTUSED*:      AUTHOR:    M.DUEREN        LAST UPDATE: 06.03.1991             :
CNOTUSED*:                                 TESTED: YES                         :
CNOTUSED*:                                                                     :
CNOTUSED*:      ARGUMENTS: DQ2,DX: DOUBLE PREC. INPUT XBJ,Q2                   :
CNOTUSED*:                 DF2H8* DOUBLE PREC F2  OUTPUT                       :
CNOTUSED*:                                                                     :
CNOTUSED*:      CALLED BY: MKF2                                                :
CNOTUSED*:                                                                     :
CNOTUSED*:      ACTION:    CALCULATE F2 STRUCTURE FUNCTION OF THE DEUTERON     :
CNOTUSED*:                 NMC FIT OF DIS-REGION WITH 8 PARAMETERS             :
CNOTUSED*:                 DATA OF NMC (1992), SLAC,BCDMS                      :
CNOTUSED*:                                                                     :
CNOTUSED*:                 PARAMETRIZED WITH A1,BI (I=1...4) AS                :
CNOTUSED*:                                                                     :
CNOTUSED*:                 F2_DIS(X,Q2) ~PROP.                                 :
CNOTUSED*:                   [1/B(N1,N2)*X**N1*(1-X)**N2 + 1/3*N3*(1-X)**N4 ]  :
CNOTUSED*:                   *S(X,Q2)                                          :
CNOTUSED*:                        WITH  X = (Q2+M_A)/(2M*NU + M_B**2)          :
CNOTUSED*:                              NI= AI+BI*S                            :
CNOTUSED*:                              S = LN(LN(Q2+M_A**2)/LAMBDA)/LN(..Q2_0):
CNOTUSED*:                 REFERENCE:                                          :
CNOTUSED*:                 THE NEW MUON COLLABORATION                          :
CNOTUSED*:                 NUCLEAR PHYSICS B 371 (1992) 3-31                   :
CNOTUSED*:=====================================================================:
CNOTUSEDC
CNOTUSED      IMPLICIT REAL*8 (D)
CNOTUSEDC
CNOTUSEDC
CNOTUSEDC *** D1,...,D8 = 8 PARAM OF NMC, SLAC, BCDMS (92)
CNOTUSEDC *** D9,...,D10 = 2 PARAMETERS: (1 FOR RESONANCE) + (1 FOR BACKGROUND)
CNOTUSEDC *** DAW,DBW =  WEIZMANN VARIABLE IN BODEK'S D2 FIT
CNOTUSEDC            VALUES: DAW=1.512(GEV2), DBW=0.351(GEV2)
CNOTUSEDC            REF:  BODEK ET AL., P.R.D20(1979)1427.
CNOTUSEDC            SEE P.1495, EQ(5.1) AND TABLE VIII
CNOTUSEDC
CNOTUSEDC *** DL2 = LAMDA**2 = 0.2**2 = 0.04 (GEV2)
CNOTUSEDC *** Q0**2 = 2 GEV2 ... (2+0.351)/0.04 = 58.771
CNOTUSEDC *** FIT BY Y.M.(25-NOV-88 19H43M14S)
CNOTUSEDC
CNOTUSED      DATA D1,D2,D3,D4,D5,D6,D7,D8
CNOTUSED     :     ,D9,D10
CNOTUSED     :     ,DAW,DBW
CNOTUSEDC
CNOTUSEDC     F2 FROM NMC, SLAC, BCDMS - DATA (92)
CNOTUSED     :     /.764053,-.171774,3.48979,.611064,.946086
CNOTUSED     :     ,1.06818,13.8472,-2.40967
CNOTUSEDC     RESONANCE-REGION
CNOTUSED     :     ,.89456,.16452
CNOTUSED     :     ,1.512,  .351 /
CNOTUSEDC
CNOTUSEDC
CNOTUSED      DF2D8_2=1.D-30
CNOTUSED      DW2 = .8803686078D0+DQ2*(1.D0/DX-1.D0)
CNOTUSED      DWMAS = DSQRT(DW2)
CNOTUSED      DDW = (DWMAS-1.03D0)
CNOTUSEDC
CNOTUSEDC *** DDW = W - (RESONANCE THRESHOLD - SMEARING 0F 0.05 GEV)
CNOTUSEDC *** LAMDA(QCD) = 0.2 GEV
CNOTUSEDC
CNOTUSED      IF(DDW.LE.0.D0) RETURN
CNOTUSEDC
CNOTUSEDC *** USE WEIZMANN VARIABLE FOR LOW Q2
CNOTUSEDC *** VALUES: DAWEIZ=1.512(GEV2), DBWEIZ=0.351(GEV2)
CNOTUSEDC *** REF:    BODEK ET AL., P.R.D20(1979)1427.
CNOTUSEDC ***         SEE P.1495, EQ(5.1) AND TABLE VIII
CNOTUSEDC
CNOTUSED      DQ2W=DQ2+DBW
CNOTUSED      DXW=DQ2W/(DQ2/DX+DAW)
CNOTUSEDC
CNOTUSED      DSBAR = DLOG(DLOG(DQ2W/.04D0)) - 1.404555751D0
CNOTUSEDC
CNOTUSED      DETAV1 = D1+D2*DSBAR
CNOTUSED      DETAV2 = D3+D4*DSBAR
CNOTUSED      DETAS1 = D5+D6*DSBAR
CNOTUSED      DETAS2 = D7+D8*DSBAR
CNOTUSEDC
CNOTUSED      DXW1=1.D0-DXW
CNOTUSED      DE1=DETAV1
CNOTUSED      DE2=DETAV2+1.D0
CNOTUSEDC
CNOTUSEDC *** SUPRESSION BUE TO "QUARK SUM RULE"
CNOTUSEDC ***       SUP.FACT.= 1 - GD**2, WITH GD=NUCLEON DIPOLE F.F.
CNOTUSEDC *** FURTHER SUPRESSION DUE TO LOW W PHASE SPACE VOLUME
CNOTUSEDC
CNOTUSED      DEN1=(1.D0+DQ2/.71D0)
CNOTUSED      DGD2=1.D0/DEN1**4
CNOTUSED      DSSUM = (1.D0-DGD2)
CNOTUSED      DS = DSSUM * (1.D0-DEXP(-4.177D0*DDW))
CNOTUSEDC
CNOTUSED      DF2D8_2 =
CNOTUSED     :     ( .83333333333333*DGAMMA(DE1+DE2)/DGAMMA(DE1)/DGAMMA(DE2)
CNOTUSED     :     *DXW**DETAV1*DXW1**DETAV2
CNOTUSED     :     + .33333333333333*DETAS1*DXW1**DETAS2  ) * DS
CNOTUSEDC$$$      SUPPRESS1 = 10.0*(1.7*DEXP(-6.0*(DQ2+0.05))+0.82 - 1.0)+1.0
CNOTUSEDC$$$      DF2H8 = DF2H8*SUPPRESS1
CNOTUSEDC
CNOTUSEDC *** RESONANCE CONTRIBUTION
CNOTUSEDC
CNOTUSED      DRES = 0.D0
CNOTUSEDC      DRES2=0.D0
CNOTUSEDC      DRES3=0.D0
CNOTUSEDC
CNOTUSEDC *** LORENTZIAN RESONANCE ( SMALL FERMI-SMEARING EFFECT)
CNOTUSEDC ***                        GAMMA(FERMI)=0.0883 GEV
CNOTUSEDC *** GAMMA(D2) = SQRT ( GAMMA(H2)**2 + GAMMA(FERMI)**2 )
CNOTUSEDC 1.232**2 = 1.518
CNOTUSEDC 1.520**2 = 2.310
CNOTUSEDC 1.681**2 = 2.826
CNOTUSEDC 1.232**2 * 0.15**2 = 0.0342
CNOTUSEDC 1.520**2 * 0.14**2 = 0.0453
CNOTUSEDC 1.681**2 * 0.14**2 = 0.0554
CNOTUSEDC
CNOTUSED      IF (DWMAS .LE. 2.3D0) THEN
CNOTUSED         DRES1 = D9**2*DEXP(-(DWMAS-1.232D0)**2/.0053D0)
CNOTUSED     :        /DEN1**3
CNOTUSEDC$$$         IF (DQ2.LT.1.0D0) THEN
CNOTUSEDC$$$            SUPPRESS1 = 0.5D0*(1.193-0.887*DQ2+0.445*DQ2**2)
CNOTUSEDC$$$         ELSE
CNOTUSEDC$$$            SUPPRESS1 = 0.5D0*0.751D0
CNOTUSEDC$$$         ENDIF
CNOTUSEDC$$$         DRES1 = DRES1*SUPPRESS1
CNOTUSEDC      DRES2 = D10**2/( (DW2-2.310D0)**2 + 0.0453D0 )
CNOTUSEDC    :        * DGD2
CNOTUSEDC      DRES3 = D11**2/( (DW2-2.826D0)**2 + 0.0554D0 )
CNOTUSEDC    :        * DGD2
CNOTUSEDC      ENDIF
CNOTUSEDC
CNOTUSEDC *** BACKGROUND UNDER RESONANCES
CNOTUSEDC
CNOTUSEDC    MP**2 = 0.8803686078
CNOTUSEDC    MP**2-M(PI)**2=0.8608892416
CNOTUSEDC *** DQS = MOMENTUM OF ONE PION, DECAYING FROM THE RESONANCE, IN CM
CNOTUSEDC FRAME
CNOTUSEDC
CNOTUSED         DW2M = (DWMAS+.05D0)**2
CNOTUSED         DQS=DSQRT((DW2M+0.8608892416D0)**2/4.D0/DW2M-0.8803686078D0)
CNOTUSED         DBG = (D10**2*DQS )
CNOTUSED     :        * DEXP(-0.5D0*DDW**2) / DEN1
CNOTUSED         DRES=(DRES1 + DBG) * DSSUM
CNOTUSED      ENDIF
CNOTUSEDC
CNOTUSEDC *** TOTAL F2 OF D2
CNOTUSEDC
CNOTUSED      DF2D8_2 = DF2D8_2 + DRES
CNOTUSEDC
CNOTUSED      IF (DF2D8_2.GT.0D0) RETURN
CNOTUSED      DF2D8_2 = 1.D-30
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSED
CNOTUSEDCDECK  ID>, DF2H8.  
CNOTUSED      REAL*8 FUNCTION DF2H8(DQ2,DX)
CNOTUSED*:=====================================================================:
CNOTUSED*:                                                                     :
CNOTUSED*:      AUTHOR:    M.DUEREN        LAST UPDATE: 06.03.1991             :
CNOTUSED*:                                 TESTED: YES                         :
CNOTUSED*:                                                                     :
CNOTUSED*:      ARGUMENTS: DQ2,DX: DOUBLE PREC. INPUT XBJ,Q2                   :
CNOTUSED*:                 DF2H8* DOUBLE PREC F2  OUTPUT                       :
CNOTUSED*:                                                                     :
CNOTUSED*:      CALLED BY: MKF2                                                :
CNOTUSED*:                                                                     :
CNOTUSED*:      ACTION:    CALCULATE F2 STRUCTURE FUNCTION OF THE PROTON       :
CNOTUSED*:                 NMC FIT OF DIS-REGION WITH 8 PARAMETERS             :
CNOTUSED*:                 DATA OF NMC (1992)SLAC,BCDMS                        :
CNOTUSED*:                                                                     :
CNOTUSED*:                 PARAMETRIZED WITH A1,BI (I=1...4) AS                :
CNOTUSED*:                                                                     :
CNOTUSED*:                 F2_DIS(X,Q2) ~PROP.                                 :
CNOTUSED*:                   [1/B(N1,N2)*X**N1*(1-X)**N2 + 1/3*N3*(1-X)**N4 ]  :
CNOTUSED*:                   *S(X,Q2)                                          :
CNOTUSED*:                        WITH  X = (Q2+M_A)/(2M*NU + M_B**2)          :
CNOTUSED*:                              NI= AI+BI*S                            :
CNOTUSED*:                              S = LN(LN(Q2+M_A**2)/LAMBDA)/LN(..Q2_0):
CNOTUSED*:                 REFERENCE:                                          :
CNOTUSED*:                 THE NEW MUON COLLABORATION                          :
CNOTUSED*:                 NUCLEAR PHYSICS B 371 (1992) 3-31                   :
CNOTUSED*:=====================================================================:
CNOTUSEDC
CNOTUSED      IMPLICIT REAL*8 (D)
CNOTUSED*
CNOTUSEDC *** D1,...,D8 = 8 PARAM OF NMC, SLAC, BCDMS (92)
CNOTUSEDC *** D9,...,D10 = 2 PARAMETERS: (1 FOR RESONANCE) + (1 FOR BACKGROUND)
CNOTUSEDC *** DAW,DBW =  WEIZMANN VARIABLE IN BODEK'S D2 FIT
CNOTUSEDC            VALUES: DAW=1.512(GEV2), DBW=0.351(GEV2)
CNOTUSEDC            REF:  BODEK ET AL., P.R.D20(1979)1427.
CNOTUSEDC            SEE P.1495, EQ(5.1) AND TABLE VIII
CNOTUSEDC
CNOTUSEDC *** DL2 = LAMDA**2 = 0.2**2 = 0.04 (GEV2)
CNOTUSEDC *** Q0**2 = 2 GEV2 ... (2+0.351)/0.04 = 58.771
CNOTUSEDC *** FIT BY Y.M.(25-NOV-88 19H43M14S)
CNOTUSED*
CNOTUSED      DATA D1,D2,D3,D4,D5,D6,D7,D8
CNOTUSED     :     ,D9,D10,D11,D12,D13,D14
CNOTUSED     :     ,D15,D16
CNOTUSED     :     ,DAW,DBW
CNOTUSEDC     F2 FROM NMC, SLAC, BCDMS  DATA '92 (FINAL)
CNOTUSED     :     /.886627,-.11191,3.3951,1.04064,1.02702,1.40335,12.4577,-
CNOTUSED     :     .100622
CNOTUSEDC     RESONANCE-REGION:
CNOTUSED     :     ,.1179, .044735, .038445, .27921, 8.8228D-5, 6.2099D-5
CNOTUSED     :     ,1.421,1.2582
CNOTUSED     :     ,1.642, .376/
CNOTUSEDC
CNOTUSED      DF2H8 =1.D-30
CNOTUSED      DW2 = .8803686078D0+DQ2*(1.D0/DX-1.D0)
CNOTUSED      DWMAS = DSQRT(DW2)
CNOTUSED      DDW = (DWMAS-1.08D0)
CNOTUSEDC
CNOTUSEDC *** DDW = W - (RESONANCE THRESHOLD)
CNOTUSEDC *** LAMDA(QCD) = D2: 0.20 GEV
CNOTUSEDC *** LAMDA(QCD) = H2: 0.15 GEV
CNOTUSEDC
CNOTUSED      IF(DDW.LE.0.D0) RETURN
CNOTUSEDC
CNOTUSEDC *** USE WEIZMANN VARIABLE FOR LOW Q2
CNOTUSEDC *** VALUES = D2 : DAWEIZ=1.512(GEV2), DBWEIZ=0.351(GEV2)
CNOTUSEDC *** VALUES = H2 : DAWEIZ=1.642(GEV2), DBWEIZ=0.376(GEV2)
CNOTUSEDC *** REF:    BODEK ET AL., P.R.D20(1979)1427.
CNOTUSEDC ***         SEE P.1495, EQ(5.1) AND TABLE VIII
CNOTUSEDC
CNOTUSED      DQ2W=DQ2+DBW
CNOTUSED      DXW=DQ2W/(DQ2/DX+DAW)
CNOTUSED      DSBAR=DLOG(DLOG(DQ2W/.0225D0)) - 1.538942135D0
CNOTUSEDC      DSBAR=DLOG(DLOG(DQ2W/.0225D0)/DLOG((2.D0+DBW)/.0225D0))
CNOTUSEDC      DSBAR = DLOG(DLOG(DQ2W/.04D0)) - 1.404555751D0
CNOTUSED      DETAV1 = D1+D2*DSBAR
CNOTUSED      DETAV2 = D3+D4*DSBAR
CNOTUSED      DETAS1 = D5+D6*DSBAR
CNOTUSED      DETAS2 = D7+D8*DSBAR
CNOTUSEDC
CNOTUSED      DXW1=1.D0-DXW
CNOTUSED      DE1=DETAV1
CNOTUSED      DE2=DETAV2+1.D0
CNOTUSEDC
CNOTUSEDC *** SUPRESSION DUE TO "QUARK SUM RULE"
CNOTUSEDC ***       SUP.FACT.= 1 - GD**2, WITH GD=NUCLEON DIPOLE F.F.
CNOTUSEDC *** FURTHER SUPRESSION DUE TO LOW W PHASE SPACE VOLUME
CNOTUSEDC
CNOTUSED      DEN1=(1.D0+DQ2/.71D0)
CNOTUSED      DGD2=1.D0/DEN1**4
CNOTUSED      DSSUM = (1.D0-DGD2)
CNOTUSED      DSTHR = 1.D0
CNOTUSED      IF( DDW .LE. 5.0D0) THEN
CNOTUSED         DTEMP = DEXP(DDW*3.98213222D0)
CNOTUSED         DSTHR = (DTEMP-1.D0)/(DTEMP+1.D0)
CNOTUSED      ENDIF
CNOTUSED      DS = DSSUM * DSTHR
CNOTUSEDC
CNOTUSED      DF2H8 =
CNOTUSED     :     ( .83333333333333*DGAMMA(DE1+DE2)/DGAMMA(DE1)/DGAMMA(DE2)
CNOTUSED     :     *DXW**DETAV1*DXW1**DETAV2
CNOTUSED     :     + .33333333333333*DETAS1*DXW1**DETAS2  ) * DS
CNOTUSEDC$$$      SUPPRESS1 = 10.0*(1.7*DEXP(-6.0*(DQ2+0.05))+0.82 - 1.0)+1.0
CNOTUSEDC$$$      DF2H8 = DF2H8*SUPPRESS1
CNOTUSEDC
CNOTUSEDC *** RESONANCE REGION
CNOTUSEDC
CNOTUSED      DRES = 0.D0
CNOTUSED
CNOTUSED      IF(DDW .LE. 5.0D0) THEN
CNOTUSEDC
CNOTUSEDC *** >>> + BACKGROUND UNDER THE RESONANCE
CNOTUSEDC
CNOTUSEDC *** APPROPRIATE KINEMATIC VARIABLES
CNOTUSEDC     ... DQS = MOMENTUM OF ONE PION IN PI-NUCLEON C.M.-SYSTEM
CNOTUSEDC                IN THE CASE OF SINGLE PION PRODUCTION
CNOTUSEDC *** N.B.  DQS = 0  AT  W = 1.08GEV (INELASTIC THRESHOLD)
CNOTUSEDC        MP**2 = 0.8803686078
CNOTUSEDC        MP**2-M(PI)**2=0.8608892416
CNOTUSEDC
CNOTUSED         DQS2 = (DW2+0.8608892416D0)**2/4.D0/DW2 - 0.8803686078D0
CNOTUSED         DQS = DSQRT(DQS2)
CNOTUSEDC
CNOTUSEDC *** >>> + RESONANCE SHAPE
CNOTUSEDC
CNOTUSEDC *** LORENTZIAN RESONANCE
CNOTUSEDC *** THIS INCLUDES THE W**2-DEPENDENCE OF THE RES.WIDTH.
CNOTUSEDC
CNOTUSEDC *** APPROPRIATE KINEMATIC VARIABLES
CNOTUSEDC  1) CORRECTION TO RES.WIDTH DUE TO RESONANCE THRESHOLD
CNOTUSEDC     ... DQS = MOMENTUM OF ONE PION IN PI-NUCLEON C.M.-SYSTEM
CNOTUSEDC                IN THE CASE OF SINGLE PION PRODUCTION
CNOTUSEDC  2) CORRECTION TO RES.WIDTH DUE TO THE Q2-DEPENDENCE
CNOTUSEDC     ... DKS = MOMENTUM OF VIRTUAL PHOTON IN PI-N C.M.-SYSTEM
CNOTUSEDC
CNOTUSED         IF(DDW .LE. 1.D0) THEN
CNOTUSEDC
CNOTUSED            DKS2 =
CNOTUSED     :           (DW2+0.8803686078D0+DQ2)**2/4.D0/DW2 - 0.8803686078D0
CNOTUSED            DKS = DSQRT(DKS2)
CNOTUSEDC
CNOTUSEDC *** RESONANCE FORM FACTOR (EFFECTIVE GUESS!)
CNOTUSEDC
CNOTUSED            DRESFF = 1. / DEN1**(D15**2)
CNOTUSEDC
CNOTUSEDC *** 1236
CNOTUSEDC     WRES**2         = 1.232**2            = 1.518
CNOTUSEDC     (WRES*GAMMA)**2 = 1.232**2 * 0.119**2 = 0.02149
CNOTUSEDC
CNOTUSED            DW2R = 1.518D0
CNOTUSED            DQSR2 =
CNOTUSED     :           (DW2R+0.8608892416D0)**2/4.D0/DW2R - 0.8803686078D0
CNOTUSED            DQSR = DSQRT(DQSR2)
CNOTUSED            DKSR2 =
CNOTUSED     :           (DW2R+0.8803686078+DQ2)**2/4./DW2R - 0.8803686078
CNOTUSED            DDCORR = (DQS/DQSR) * (1.+.16/DQSR2)/(1.+.16/DQS2)
CNOTUSED            DNCORR = DDCORR * (DKSR2+.16)/(DKS2+.16)
CNOTUSED            DDCORR = DDCORR**2
CNOTUSED            DRES1 = D9**2 * DNCORR
CNOTUSED     :           /( (DW2-1.518D0)**2 + 0.02149D0*DDCORR )
CNOTUSEDC$$$            IF (DQ2.LT.1.0D0) THEN
CNOTUSEDC$$$               SUPPRESS1 = 0.5D0*(1.193-0.887*DQ2+0.445*DQ2**2)
CNOTUSEDC$$$            ELSE
CNOTUSEDC$$$               SUPPRESS1 = 0.5D0*0.751D0
CNOTUSEDC$$$            ENDIF
CNOTUSEDC$$$            DRES1 = DRES1*SUPPRESS1
CNOTUSEDC
CNOTUSEDC *** 1520
CNOTUSEDC     WRES**2         = 1.520**2            = 2.310
CNOTUSEDC     (WRES*GAMMA)**2 = 1.520**2 * 0.097**2 = 0.02127
CNOTUSEDC     N.B. Q2-DEPENDENCE OF THE WIDTH IS NEGLECTED
CNOTUSEDC
CNOTUSED            DW2R = 2.310D0
CNOTUSED            DQSR2 =
CNOTUSED     :           (DW2R+0.8608892416D0)**2/4.D0/DW2R - 0.8803686078D0
CNOTUSED            DDCORR = DQS2/DQSR2
CNOTUSED            DRES2 = D10**2 * DDCORR
CNOTUSED     :           / ( (DW2-2.310D0)**2 + 0.02127D0*DDCORR )
CNOTUSEDC
CNOTUSEDC *** 1681
CNOTUSEDC     WRES**2         = 1.681**2            = 2.826
CNOTUSEDC     (WRES*GAMMA)**2 = 1.681**2 * 0.105**2 = 0.03115
CNOTUSEDC     N.B. Q2-DEPENDENCE OF THE WIDTH IS NEGLECTED
CNOTUSEDC
CNOTUSED            DW2R = 2.826D0
CNOTUSED            DQSR2 =
CNOTUSED     :           (DW2R+0.8608892416D0)**2/4.D0/DW2R - 0.8803686078D0
CNOTUSED            DQSR = DSQRT(DQSR2)
CNOTUSED            DDCORR = ( DQS/DQSR )**3
CNOTUSED            DRES3 = D11**2 * DDCORR
CNOTUSED     :           / ( (DW2-2.826D0)**2 + 0.03115D0*DDCORR )
CNOTUSEDC
CNOTUSEDC
CNOTUSEDC *** SUM OF ALL RESONANCES
CNOTUSEDC         * RESONANCE FORM FACTOR(Q2-DEPENDENCE)
CNOTUSEDC
CNOTUSED            DRES = (DRES1 +DRES2 +DRES3)*DRESFF
CNOTUSEDC
CNOTUSEDC *** END OF RESONANCE CALCULATION (ONLY IF   DDW < 1.0 GEV)
CNOTUSEDC
CNOTUSED         ENDIF
CNOTUSEDC
CNOTUSEDC *** BACKGROUND UNDER THE RESONANCES
CNOTUSEDC     N.B. EXP(-0.92**2*3.5393817) = 0.05
CNOTUSEDC        DBGFF = DQ2VAL/DXVAL /DEN1**(DP(16)**2)
CNOTUSEDC
CNOTUSED         DBGFF = 1. / DEN1**(D16**2)
CNOTUSED         DBG = (D12**2*DSQRT(DQS) +D13**2*DQS +D14**2*DQS2 )
CNOTUSED     :        * DBGFF * DEXP(-0.5D0*DDW**2)
CNOTUSEDC
CNOTUSEDC *** (RESONANCE REGION) = ( (RESONANCE) + (BACKGROUND) )
CNOTUSEDC                          * DSSUM(=SUPPRESSION)
CNOTUSEDC
CNOTUSED         DRES = (DRES + DBG) * DSSUM
CNOTUSEDC
CNOTUSEDC *** END OF RESONANCE CALCULATION
CNOTUSEDC
CNOTUSED      ENDIF
CNOTUSEDC
CNOTUSEDC *** (TOTAL) = (QCD PART) + (RESONANCE REGION)
CNOTUSEDC
CNOTUSED      DF2H8 = DF2H8 + DRES
CNOTUSEDC
CNOTUSED      IF(DF2H8 .GT. 0.D0) RETURN
CNOTUSEDC
CNOTUSED      DF2H8=1.D-30
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSED
CNOTUSEDC
CNOTUSEDC----------------------------------------------------------------------
CNOTUSEDC
CNOTUSED      REAL*8 FUNCTION F1_QE_SCOP_2(X,Q2)
CNOTUSED      IMPLICIT REAL*8 (A-H,O-Z)
CNOTUSED      DIMENSION C(10)
CNOTUSED
CNOTUSED      DATA XMP/0.938272/,XMP2/0.8803543459/,PI/3.1415926535/,
CNOTUSED     &     C/-10.121,28.493,-27.836,9.6378,-0.73175,
CNOTUSED     &       3.8071,-7.9896,5.8759,-3.7466,0.95000/
CNOTUSED
CNOTUSED      XNORM = POLYNOM_2(Q2,4,C(1))
CNOTUSED      IF (XNORM.GT.-80) THEN
CNOTUSED         XNORM = DEXP(XNORM)
CNOTUSED      ELSE
CNOTUSED         XNORM = 0.0
CNOTUSED      ENDIF
CNOTUSED      IF (XNORM.GT.0.0) THEN
CNOTUSED         WIDTH = POLYNOM_2(Q2,3,C(6))
CNOTUSED         WIDTH = 0.5*DEXP(WIDTH)
CNOTUSED      ELSE
CNOTUSED         WIDTH = 1.0D30
CNOTUSED      ENDIF
CNOTUSED      
CNOTUSED      W2 = XMP2 + (1.0/X-1.0)*Q2
CNOTUSED      IF (W2.GE.0.0) THEN
CNOTUSED         W = SQRT(W2)
CNOTUSED      ELSE
CNOTUSED         W = -SQRT(-W2)
CNOTUSED      ENDIF
CNOTUSED      F1_QE_SCOP_2 = XNORM*WIDTH/(PI*((W-C(10))**2+WIDTH**2))
CNOTUSEDC$$$  F1_QE_SCOP_2 = F1_QE_SCOP_2*0.2*(1.0-Q2*0.03)
CNOTUSEDC$$$  CALL HFILL(200,SNGL(W),SNGL(Q2),SNGL(F1_QE_SCOP_2))
CNOTUSEDC$$$  CALL HFILL(201,SNGL(W),SNGL(Q2),1.0)
CNOTUSEDC$$$  CALL HFILL(202,SNGL(X),SNGL(Q2),SNGL(F1_QE_SCOP_2))
CNOTUSEDC$$$  CALL HFILL(203,SNGL(X),SNGL(Q2),1.0)
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSEDC
CNOTUSEDC----------------------------------------------------------------------
CNOTUSEDC
CNOTUSED      REAL*8 FUNCTION F2_QE_SCOP_2(X,Q2)
CNOTUSED      IMPLICIT REAL*8 (A-H,O-Z)
CNOTUSED      DIMENSION C(10)
CNOTUSED
CNOTUSED      DATA XMP/0.938272/,XMP2/0.8803543459/,PI/3.1415926535/,
CNOTUSED     &     C/-11.312,29.037,-25.943,9.0082,-2.8383,
CNOTUSED     &       3.7698,-7.4814,5.2376,-3.5727,0.95526/
CNOTUSED
CNOTUSED      XNORM = POLYNOM_2(Q2,4,C(1))
CNOTUSED      IF (XNORM.GT.-80) THEN
CNOTUSED         XNORM = DEXP(XNORM)
CNOTUSED      ELSE
CNOTUSED         XNORM = 0.0
CNOTUSED      ENDIF
CNOTUSED      IF (XNORM.GT.0.0) THEN
CNOTUSED         WIDTH = POLYNOM_2(Q2,3,C(6))
CNOTUSED         WIDTH = 0.5*DEXP(WIDTH)
CNOTUSED      ELSE
CNOTUSED         WIDTH = 1.0D30
CNOTUSED      ENDIF
CNOTUSED
CNOTUSED      W2 = XMP2 + (1.0/X-1.0)*Q2
CNOTUSED      IF (W2.GE.0.0) THEN
CNOTUSED         W = SQRT(W2)
CNOTUSED      ELSE
CNOTUSED         W = -SQRT(-W2)
CNOTUSED      ENDIF
CNOTUSED      F2_QE_SCOP_2 = XNORM*WIDTH/(PI*((W-C(10))**2+WIDTH**2))
CNOTUSED
CNOTUSEDC$$$  F2_QE_SCOP_2 = F2_QE_SCOP_2*0.2*(1.0-Q2*0.03)
CNOTUSEDC$$$  CALL HFILL(210,SNGL(W),SNGL(Q2),SNGL(F2_QE_SCOP_2))
CNOTUSEDC$$$  CALL HFILL(211,SNGL(W),SNGL(Q2),1.0)
CNOTUSEDC$$$  CALL HFILL(212,SNGL(X),SNGL(Q2),SNGL(F2_QE_SCOP_2))
CNOTUSEDC$$$  CALL HFILL(213,SNGL(X),SNGL(Q2),1.0)
CNOTUSED      RETURN
CNOTUSED      END 
CNOTUSEDC
CNOTUSEDC----------------------------------------------------------------------
CNOTUSEDC
CNOTUSED
CNOTUSED      REAL*8 FUNCTION G1_QE_SCOP_2(X,Q2)
CNOTUSED      IMPLICIT REAL*8 (A-H,O-Z)
CNOTUSED      DIMENSION C(18)
CNOTUSED
CNOTUSED      DATA XMP/0.938272/,XMP2/0.8803543459/,PI/3.1415926535/,
CNOTUSED     &     C/-0.44581E-01, 0.77532E-01,-0.36036    , 0.61889E-01,-0.52041E-01,
CNOTUSED     &       -1.1221    , 0.94397    ,-0.73179E-02, 0.39891E-01, 0.24352E-01,
CNOTUSED     &       -15.881    ,-0.35483E-01, 0.68196E-01,  1.1320    ,  1.0018    ,
CNOTUSED     &       -0.62122E-01, -11.933    , 0.16062E-01/
CNOTUSED
CNOTUSED      XNORM1 = C(1)+C(2)*EXP(C(3)*Q2)
CNOTUSED      WIDTH1 = C(4)+C(5)*EXP(C(6)*Q2)
CNOTUSED      POSITION1 = C(7)+C(8)*Q2
CNOTUSED      XNORM2 = C(9)+C(10)*EXP(C(11)*Q2)
CNOTUSED      WIDTH2 = C(12)+C(13)*EXP(C(14)*Q2)
CNOTUSED      POSITION2 = C(15)+C(16)*EXP(C(17)*Q2)
CNOTUSED      OFFSET = C(18)
CNOTUSED      
CNOTUSED      W2 = XMP2 + (1.0/X-1.0)*Q2
CNOTUSED      IF (W2.GE.0.0) THEN
CNOTUSED         W = SQRT(W2)
CNOTUSED      ELSE
CNOTUSED         W = -SQRT(-W2)
CNOTUSED      ENDIF
CNOTUSED      
CNOTUSED      G1_QE_SCOP_2 = XNORM1*WIDTH1/(PI*((W-POSITION1)**2+WIDTH1**2))
CNOTUSED     &           - XNORM2*WIDTH2/(PI*((W-POSITION2)**2+WIDTH2**2))
CNOTUSED
CNOTUSEDC$$$  G1_QE_SCOP_2 = G1_QE_SCOP_2*0.1D0
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSEDC
CNOTUSEDC----------------------------------------------------------------------
CNOTUSEDC
CNOTUSED      REAL*8 FUNCTION G2_QE_SCOP_2(X,Q2)
CNOTUSED      IMPLICIT REAL*8 (A-H,O-Z)
CNOTUSED      DIMENSION C(9)
CNOTUSED
CNOTUSED      DATA XMP/0.938272/,XMP2/0.8803543459/,PI/3.1415926535/,
CNOTUSED     &     C/ 0.37890E-02, 0.13797    , -1.6161    ,-0.14017    , 0.16467    ,
CNOTUSED     &        0.32133    , 0.97506    , 0.89238E-02,-0.79519E-02/
CNOTUSED
CNOTUSED      XNORM1 = C(1)+C(2)*Q2*EXP(C(3)*Q2)
CNOTUSED      WIDTH1 = C(4)+C(5)*EXP(C(6)*Q2)
CNOTUSED      POSITION1 = C(7)+C(8)*Q2
CNOTUSED      OFFSET = C(9)
CNOTUSED      
CNOTUSED      W2 = XMP2 + (1.0/X-1.0)*Q2
CNOTUSED      IF (W2.GE.0.0) THEN
CNOTUSED         W = SQRT(W2)
CNOTUSED      ELSE
CNOTUSED         W = -SQRT(-W2)
CNOTUSED      ENDIF
CNOTUSED      
CNOTUSED      G2_QE_SCOP_2 = XNORM1*WIDTH1/(PI*((W-POSITION1)**2+WIDTH1**2))
CNOTUSEDC$$$  G2_QE_SCOP_2 = G2_QE_SCOP_2*0.1D0
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSED      
CNOTUSED
CNOTUSED      REAL*8 FUNCTION POLYNOM_2(X,N,C)
CNOTUSED      REAL*8 X,C
CNOTUSED      INTEGER N,I
CNOTUSED      DIMENSION C(*)
CNOTUSED
CNOTUSED      POLYNOM_2 = C(1)
CNOTUSED      DO I = 2,N+1
CNOTUSED         POLYNOM_2 = POLYNOM_2*X + C(I)
CNOTUSED      ENDDO
CNOTUSED      RETURN
CNOTUSED      END
CNOTUSED
