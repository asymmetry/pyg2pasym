******from Amroun et al data****************
C     AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)
C     A. Deur
      REAL*8 FUNCTION FORMC(Q2)
      IMPLICIT REAL*8 (A-H,O-Z)
c     Function to compute a SOG form factor
c     Q2 : momentum transfer squared (fm-2)
C     NR : number of Gaussians
C     GA : Gaussians rms (usually 0.8)
C     RR : Gaussians central positions
C     QQ : Gaussians amplitudes
      real*8 RR(12),QQ(12)
      real*8  Q2,GA,A,B
      DATA NR/12/,GA/0.65D0/
      DATA RR/0.1D0,0.5D0,0.9D0,1.3D0,1.6D0,2.0D0,2.4D0,
     &        2.9D0,3.4D0,4.0D0,4.6D0,5.2D0/
      DATA QQ/.027614D0,.170847D0,.219805D0,.170486D0,
     &        .134453D0,.100953D0,.074310D0,.053970D0,
     +        .023689D0,.017502D0,.002034D0,.004338D0/
      Q=DSQRT(Q2)
      G2=GA*GA
      A=DEXP(-Q2*G2/4.D0)
      S=0.D0
      DO I=1,NR
          B=2.D0*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.D0+B
          ELSE
              SS=DCOS(QR)+B*DSIN(QR)/QR
          END IF
          SS=QQ(I)/(1.D0+B)*SS
          S=S+SS
      END DO
      FORMC=A*S
      RETURN
      END


      REAL*8 FUNCTION FORMM(Q2)
C     AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)
C     A. Deur
      IMPLICIT REAL*8 (A-H,O-Z)
c     Function to compute a SOG form factor
c     Q2 : momentum transfer squared (fm-2)
c     NR : number of Gaussians
c     GA : Gaussians rms (usually 0.8)
C     RR : Gaussians central positions
C     QQ : Gaussians amplitudes
      real*8 RR(12),QQ(12)
      real*8  Q2,GA,A,B
      DATA NR/12/,GA/0.65D0/
      DATA RR/0.1D0,0.5D0,0.9D0,1.3D0,1.6D0,2.0D0,2.4D0,2.9D0,
     &        3.4D0,4.0D0,4.6D0,5.2D0/
      DATA QQ/.059785D0,.138368D0,.281326D0,.000037D0,
     &        .289808D0,.019056D0,.114825D0,.042296D0,
     +        .028345D0,.018312D0,.007843D0,.000000D0/
      Q=DSQRT(Q2)
      G2=GA*GA
      A=DEXP(-Q2*G2/4.D0)
      S=0.D0
      DO I=1,NR
          B=2.D0*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.D0+B
          ELSE
              SS=DCOS(QR)+B*DSIN(QR)/QR
          END IF
          SS=QQ(I)/(1.D0+B)*SS
          S=S+SS
      END DO
      FORMM=A*S
      RETURN
      END

