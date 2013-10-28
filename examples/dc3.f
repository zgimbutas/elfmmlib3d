c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains Okada's formulas for elastostatic potentials
c       in R^3 due to point sources and rectangles in lower half space.
c
c       Half space elastostatic Green's function (Mindlin's solution).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      SUBROUTINE  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,      00010000
     *               UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 00020002
      IMPLICIT REAL*8 (A-H,O-Z)                                         00030000
      REAL*8   ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,               00040000
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             00050000
C                                                                       00060000
C********************************************************************   00070000
C*****                                                          *****   00080000
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   00090000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   00100000
C*****                         CODED BY  Y.OKADA ... SEP.1991   *****   00110002
C*****                         REVISED   Y.OKADA ... NOV.1991   *****   00120002
C*****                                                          *****   00130000
C********************************************************************   00140000
C                                                                       00150000
C***** INPUT                                                            00160000
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           00170000
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          00180000
C*****   DEPTH : SOURCE DEPTH                                           00190000
C*****   DIP   : DIP-ANGLE (DEGREE)                                     00200000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        00210000
C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2   00220000
C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3     00230000
C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4     00240000
C                                                                       00250000
C***** OUTPUT                                                           00260000
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /          00270000
C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )  00280000
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /          00290000
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )  00300000
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     00310000
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  00320002
C                                                                       00330000
      COMMON /C1/DUMMY(8),R,DUMMY2(15)                                  00340000
      DIMENSION  U(12),DUA(12),DUB(12),DUC(12)                          00350000
      DATA  F0/0.D0/                                                    00360000
C-----                                                                  00370000
      IF(Z.GT.0.) WRITE(6,'(''0** POSITIVE Z WAS GIVEN IN SUB-DC3D0'')')00380000
      DO 111 I=1,12                                                     00390000
        U(I)=F0                                                         00400000
        DUA(I)=F0                                                       00410000
        DUB(I)=F0                                                       00420000
        DUC(I)=F0                                                       00430000
  111 CONTINUE                                                          00440000
      AALPHA=ALPHA                                                      00450000
      DDIP=DIP                                                          00460000
      CALL DCCON0(AALPHA,DDIP)                                          00470000
C======================================                                 00480000
C=====  REAL-SOURCE CONTRIBUTION  =====                                 00490000
C======================================                                 00500000
      XX=X                                                              00510000
      YY=Y                                                              00520000
      ZZ=Z                                                              00530000
      DD=DEPTH+Z                                                        00540000
      CALL DCCON1(XX,YY,DD)                                             00550000
      IF(R.EQ.F0) GO TO 99                                              00560000
      PP1=POT1                                                          00570000
      PP2=POT2                                                          00580000
      PP3=POT3                                                          00590000
      PP4=POT4                                                          00600000
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            00610000
C-----                                                                  00620000
      DO 222 I=1,12                                                     00630000
        IF(I.LT.10) U(I)=U(I)-DUA(I)                                    00640000
        IF(I.GE.10) U(I)=U(I)+DUA(I)                                    00650000
  222 CONTINUE                                                          00660000
C=======================================                                00670000
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                00680000
C=======================================                                00690000
      DD=DEPTH-Z                                                        00700000
      CALL DCCON1(XX,YY,DD)                                             00710000
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            00720000
      CALL UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUB)                         00730000
      CALL UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUC)                         00740000
C-----                                                                  00750000
      DO 333 I=1,12                                                     00760000
        DU=DUA(I)+DUB(I)+ZZ*DUC(I)                                      00770000
        IF(I.GE.10) DU=DU+DUC(I-9)                                      00780000
        U(I)=U(I)+DU                                                    00790000
  333 CONTINUE                                                          00800000
C=====                                                                  00810000
      UX=U(1)                                                           00820000
      UY=U(2)                                                           00830000
      UZ=U(3)                                                           00840000
      UXX=U(4)                                                          00850000
      UYX=U(5)                                                          00860000
      UZX=U(6)                                                          00870000
      UXY=U(7)                                                          00880000
      UYY=U(8)                                                          00890000
      UZY=U(9)                                                          00900000
      UXZ=U(10)                                                         00910000
      UYZ=U(11)                                                         00920000
      UZZ=U(12)                                                         00930000
      IRET=0                                                            00940002
      RETURN                                                            00950000
C=======================================                                00960000
C=====  IN CASE OF SINGULAR (R=0)  =====                                00970000
C=======================================                                00980000
   99 UX=F0                                                             00990000
      UY=F0                                                             01000000
      UZ=F0                                                             01010000
      UXX=F0                                                            01020000
      UYX=F0                                                            01030000
      UZX=F0                                                            01040000
      UXY=F0                                                            01050000
      UYY=F0                                                            01060000
      UZY=F0                                                            01070000
      UXZ=F0                                                            01080000
      UYZ=F0                                                            01090000
      UZZ=F0                                                            01100000
      IRET=1                                                            01110002
      RETURN                                                            01120000
      END                                                               01130000
      SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U)                      01140000
      IMPLICIT REAL*8 (A-H,O-Z)                                         01150000
      DIMENSION U(12),DU(12)                                            01160000
C                                                                       01170000
C********************************************************************   01180000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   01190000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   01200000
C********************************************************************   01210000
C                                                                       01220000
C***** INPUT                                                            01230000
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    01240000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        01250000
C***** OUTPUT                                                           01260000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     01270000
C                                                                       01280000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  01290000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     01300000
     *           UY,VY,WY,UZ,VZ,WZ                                      01310000
      DATA F0,F1,F3/0.D0,1.D0,3.D0/                                     01320000
      DATA PI2/6.283185307179586D0/                                     01330000
C-----                                                                  01340000
      DO 111  I=1,12                                                    01350000
  111 U(I)=F0                                                           01360000
C======================================                                 01370000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 01380000
C======================================                                 01390000
      IF(POT1.NE.F0) THEN                                               01400000
        DU( 1)= ALP1*Q/R3    +ALP2*X2*QR                                01410000
        DU( 2)= ALP1*X/R3*SD +ALP2*XY*QR                                01420000
        DU( 3)=-ALP1*X/R3*CD +ALP2*X*D*QR                               01430000
        DU( 4)= X*QR*(-ALP1 +ALP2*(F1+A5) )                             01440000
        DU( 5)= ALP1*A3/R3*SD +ALP2*Y*QR*A5                             01450000
        DU( 6)=-ALP1*A3/R3*CD +ALP2*D*QR*A5                             01460000
        DU( 7)= ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY                     01470000
        DU( 8)= F3*X/R5*(-ALP1*Y*SD +ALP2*(Y*UY+Q) )                    01480000
        DU( 9)= F3*X/R5*( ALP1*Y*CD +ALP2*D*UY )                        01490000
        DU(10)= ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ                     01500000
        DU(11)= F3*X/R5*( ALP1*D*SD +ALP2*Y*UZ )                        01510000
        DU(12)= F3*X/R5*(-ALP1*D*CD +ALP2*(D*UZ-Q) )                    01520000
        DO 222 I=1,12                                                   01530000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        01540000
      ENDIF                                                             01550000
C===================================                                    01560000
C=====  DIP-SLIP CONTRIBUTION  =====                                    01570000
C===================================                                    01580000
      IF(POT2.NE.F0) THEN                                               01590000
        DU( 1)=            ALP2*X*P*QR                                  01600000
        DU( 2)= ALP1*S/R3 +ALP2*Y*P*QR                                  01610000
        DU( 3)=-ALP1*T/R3 +ALP2*D*P*QR                                  01620000
        DU( 4)=                 ALP2*P*QR*A5                            01630000
        DU( 5)=-ALP1*F3*X*S/R5 -ALP2*Y*P*QRX                            01640000
        DU( 6)= ALP1*F3*X*T/R5 -ALP2*D*P*QRX                            01650000
        DU( 7)=                          ALP2*F3*X/R5*VY                01660000
        DU( 8)= ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR)         01670000
        DU( 9)=-ALP1*(C2D/R3-F3*Y*T/R5) +ALP2*F3*D/R5*VY                01680000
        DU(10)=                          ALP2*F3*X/R5*VZ                01690000
        DU(11)= ALP1*(C2D/R3+F3*D*S/R5) +ALP2*F3*Y/R5*VZ                01700000
        DU(12)= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR)         01710000
        DO 333 I=1,12                                                   01720000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        01730000
      ENDIF                                                             01740000
C========================================                               01750000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               01760000
C========================================                               01770000
      IF(POT3.NE.F0) THEN                                               01780000
        DU( 1)= ALP1*X/R3 -ALP2*X*Q*QR                                  01790000
        DU( 2)= ALP1*T/R3 -ALP2*Y*Q*QR                                  01800000
        DU( 3)= ALP1*S/R3 -ALP2*D*Q*QR                                  01810000
        DU( 4)= ALP1*A3/R3     -ALP2*Q*QR*A5                            01820000
        DU( 5)=-ALP1*F3*X*T/R5 +ALP2*Y*Q*QRX                            01830000
        DU( 6)=-ALP1*F3*X*S/R5 +ALP2*D*Q*QRX                            01840000
        DU( 7)=-ALP1*F3*XY/R5           -ALP2*X*QR*WY                   01850000
        DU( 8)= ALP1*(C2D/R3-F3*Y*T/R5) -ALP2*(Y*WY+Q)*QR               01860000
        DU( 9)= ALP1*(S2D/R3-F3*Y*S/R5) -ALP2*D*QR*WY                   01870000
        DU(10)= ALP1*F3*X*D/R5          -ALP2*X*QR*WZ                   01880000
        DU(11)=-ALP1*(S2D/R3-F3*D*T/R5) -ALP2*Y*QR*WZ                   01890000
        DU(12)= ALP1*(C2D/R3+F3*D*S/R5) -ALP2*(D*WZ-Q)*QR               01900000
        DO 444 I=1,12                                                   01910000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        01920000
      ENDIF                                                             01930000
C=========================================                              01940000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              01950000
C=========================================                              01960000
      IF(POT4.NE.F0) THEN                                               01970000
        DU( 1)=-ALP1*X/R3                                               01980000
        DU( 2)=-ALP1*Y/R3                                               01990000
        DU( 3)=-ALP1*D/R3                                               02000000
        DU( 4)=-ALP1*A3/R3                                              02010000
        DU( 5)= ALP1*F3*XY/R5                                           02020000
        DU( 6)= ALP1*F3*X*D/R5                                          02030000
        DU( 7)= DU(5)                                                   02040000
        DU( 8)=-ALP1*B3/R3                                              02050000
        DU( 9)= ALP1*F3*Y*D/R5                                          02060000
        DU(10)=-DU(6)                                                   02070000
        DU(11)=-DU(9)                                                   02080000
        DU(12)= ALP1*C3/R3                                              02090000
        DO 555 I=1,12                                                   02100000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        02110000
      ENDIF                                                             02120000
      RETURN                                                            02130000
      END                                                               02140000
      SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    02150000
      IMPLICIT REAL*8 (A-H,O-Z)                                         02160000
      DIMENSION U(12),DU(12)                                            02170000
C                                                                       02180000
C********************************************************************   02190000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   02200000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   02210000
C********************************************************************   02220000
C                                                                       02230000
C***** INPUT                                                            02240000
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  02250000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        02260000
C***** OUTPUT                                                           02270000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     02280000
C                                                                       02290000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  02300000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     02310000
     *           UY,VY,WY,UZ,VZ,WZ                                      02320000
      DATA F0,F1,F2,F3,F4,F5,F8,F9                                      02330000
     *        /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,8.D0,9.D0/                 02340000
      DATA PI2/6.283185307179586D0/                                     02350000
C-----                                                                  02360000
      C=D+Z                                                             02370000
      RD=R+D                                                            02380000
      D12=F1/(R*RD*RD)                                                  02390000
      D32=D12*(F2*R+D)/R2                                               02400000
      D33=D12*(F3*R+D)/(R2*RD)                                          02410000
      D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)                           02420000
      D54=D12*(F5*R2+F4*R*D+D2)/R3*D12                                  02430000
C-----                                                                  02440000
      FI1= Y*(D12-X2*D33)                                               02450000
      FI2= X*(D12-Y2*D33)                                               02460000
      FI3= X/R3-FI2                                                     02470000
      FI4=-XY*D32                                                       02480000
      FI5= F1/(R*RD)-X2*D32                                             02490000
      FJ1=-F3*XY*(D33-X2*D54)                                           02500000
      FJ2= F1/R3-F3*D12+F3*X2*Y2*D54                                    02510000
      FJ3= A3/R3-FJ2                                                    02520000
      FJ4=-F3*XY/R5-FJ1                                                 02530000
      FK1=-Y*(D32-X2*D53)                                               02540000
      FK2=-X*(D32-Y2*D53)                                               02550000
      FK3=-F3*X*D/R5-FK2                                                02560000
C-----                                                                  02570000
      DO 111  I=1,12                                                    02580000
  111 U(I)=F0                                                           02590000
C======================================                                 02600000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 02610000
C======================================                                 02620000
      IF(POT1.NE.F0) THEN                                               02630000
        DU( 1)=-X2*QR  -ALP3*FI1*SD                                     02640000
        DU( 2)=-XY*QR  -ALP3*FI2*SD                                     02650000
        DU( 3)=-C*X*QR -ALP3*FI4*SD                                     02660000
        DU( 4)=-X*QR*(F1+A5) -ALP3*FJ1*SD                               02670000
        DU( 5)=-Y*QR*A5      -ALP3*FJ2*SD                               02680000
        DU( 6)=-C*QR*A5      -ALP3*FK1*SD                               02690000
        DU( 7)=-F3*X2/R5*UY      -ALP3*FJ2*SD                           02700000
        DU( 8)=-F3*XY/R5*UY-X*QR -ALP3*FJ4*SD                           02710000
        DU( 9)=-F3*C*X/R5*UY     -ALP3*FK2*SD                           02720000
        DU(10)=-F3*X2/R5*UZ  +ALP3*FK1*SD                               02730000
        DU(11)=-F3*XY/R5*UZ  +ALP3*FK2*SD                               02740000
        DU(12)= F3*X/R5*(-C*UZ +ALP3*Y*SD)                              02750000
        DO 222 I=1,12                                                   02760000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        02770000
      ENDIF                                                             02780000
C===================================                                    02790000
C=====  DIP-SLIP CONTRIBUTION  =====                                    02800000
C===================================                                    02810000
      IF(POT2.NE.F0) THEN                                               02820000
        DU( 1)=-X*P*QR +ALP3*FI3*SDCD                                   02830000
        DU( 2)=-Y*P*QR +ALP3*FI1*SDCD                                   02840000
        DU( 3)=-C*P*QR +ALP3*FI5*SDCD                                   02850000
        DU( 4)=-P*QR*A5 +ALP3*FJ3*SDCD                                  02860000
        DU( 5)= Y*P*QRX +ALP3*FJ1*SDCD                                  02870000
        DU( 6)= C*P*QRX +ALP3*FK3*SDCD                                  02880000
        DU( 7)=-F3*X/R5*VY      +ALP3*FJ1*SDCD                          02890000
        DU( 8)=-F3*Y/R5*VY-P*QR +ALP3*FJ2*SDCD                          02900000
        DU( 9)=-F3*C/R5*VY      +ALP3*FK1*SDCD                          02910000
        DU(10)=-F3*X/R5*VZ -ALP3*FK3*SDCD                               02920000
        DU(11)=-F3*Y/R5*VZ -ALP3*FK1*SDCD                               02930000
        DU(12)=-F3*C/R5*VZ +ALP3*A3/R3*SDCD                             02940000
        DO 333 I=1,12                                                   02950000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        02960000
      ENDIF                                                             02970000
C========================================                               02980000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               02990000
C========================================                               03000000
      IF(POT3.NE.F0) THEN                                               03010000
        DU( 1)= X*Q*QR -ALP3*FI3*SDSD                                   03020000
        DU( 2)= Y*Q*QR -ALP3*FI1*SDSD                                   03030000
        DU( 3)= C*Q*QR -ALP3*FI5*SDSD                                   03040000
        DU( 4)= Q*QR*A5 -ALP3*FJ3*SDSD                                  03050000
        DU( 5)=-Y*Q*QRX -ALP3*FJ1*SDSD                                  03060000
        DU( 6)=-C*Q*QRX -ALP3*FK3*SDSD                                  03070000
        DU( 7)= X*QR*WY     -ALP3*FJ1*SDSD                              03080000
        DU( 8)= QR*(Y*WY+Q) -ALP3*FJ2*SDSD                              03090000
        DU( 9)= C*QR*WY     -ALP3*FK1*SDSD                              03100000
        DU(10)= X*QR*WZ +ALP3*FK3*SDSD                                  03110000
        DU(11)= Y*QR*WZ +ALP3*FK1*SDSD                                  03120000
        DU(12)= C*QR*WZ -ALP3*A3/R3*SDSD                                03130000
        DO 444 I=1,12                                                   03140000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        03150000
      ENDIF                                                             03160000
C=========================================                              03170000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              03180000
C=========================================                              03190000
      IF(POT4.NE.F0) THEN                                               03200000
        DU( 1)= ALP3*X/R3                                               03210000
        DU( 2)= ALP3*Y/R3                                               03220000
        DU( 3)= ALP3*D/R3                                               03230000
        DU( 4)= ALP3*A3/R3                                              03240000
        DU( 5)=-ALP3*F3*XY/R5                                           03250000
        DU( 6)=-ALP3*F3*X*D/R5                                          03260000
        DU( 7)= DU(5)                                                   03270000
        DU( 8)= ALP3*B3/R3                                              03280000
        DU( 9)=-ALP3*F3*Y*D/R5                                          03290000
        DU(10)=-DU(6)                                                   03300000
        DU(11)=-DU(9)                                                   03310000
        DU(12)=-ALP3*C3/R3                                              03320000
        DO 555 I=1,12                                                   03330000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        03340000
      ENDIF                                                             03350000
      RETURN                                                            03360000
      END                                                               03370000
      SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    03380000
      IMPLICIT REAL*8 (A-H,O-Z)                                         03390000
      DIMENSION U(12),DU(12)                                            03400000
C                                                                       03410000
C********************************************************************   03420000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   03430000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   03440000
C********************************************************************   03450000
C                                                                       03460000
C***** INPUT                                                            03470000
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  03480000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        03490000
C***** OUTPUT                                                           03500000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     03510000
C                                                                       03520000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  03530000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     02310000
     *           UY,VY,WY,UZ,VZ,WZ                                      02320000
c     COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3      03540000
      DATA F0,F1,F2,F3,F5,F7,F10,F15                                    03550000
     *        /0.D0,1.D0,2.D0,3.D0,5.D0,7.D0,10.D0,15.D0/               03560000
      DATA PI2/6.283185307179586D0/                                     03570000
C-----                                                                  03580000
      C=D+Z                                                             03590000
      Q2=Q*Q                                                            03600000
      R7=R5*R2                                                          03610000
      A7=F1-F7*X2/R2                                                    03620000
      B5=F1-F5*Y2/R2                                                    03630000
      B7=F1-F7*Y2/R2                                                    03640000
      C5=F1-F5*D2/R2                                                    03650000
      C7=F1-F7*D2/R2                                                    03660000
      D7=F2-F7*Q2/R2                                                    03670000
      QR5=F5*Q/R2                                                       03680000
      QR7=F7*Q/R2                                                       03690000
      DR5=F5*D/R2                                                       03700000
C-----                                                                  03710000
      DO 111  I=1,12                                                    03720000
  111 U(I)=F0                                                           03730000
C======================================                                 03740000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 03750000
C======================================                                 03760000
      IF(POT1.NE.F0) THEN                                               03770000
        DU( 1)=-ALP4*A3/R3*CD  +ALP5*C*QR*A5                            03780000
        DU( 2)= F3*X/R5*( ALP4*Y*CD +ALP5*C*(SD-Y*QR5) )                03790000
        DU( 3)= F3*X/R5*(-ALP4*Y*SD +ALP5*C*(CD+D*QR5) )                03800000
        DU( 4)= ALP4*F3*X/R5*(F2+A5)*CD   -ALP5*C*QRX*(F2+A7)           03810000
        DU( 5)= F3/R5*( ALP4*Y*A5*CD +ALP5*C*(A5*SD-Y*QR5*A7) )         03820000
        DU( 6)= F3/R5*(-ALP4*Y*A5*SD +ALP5*C*(A5*CD+D*QR5*A7) )         03830000
        DU( 7)= DU(5)                                                   03840000
        DU( 8)= F3*X/R5*( ALP4*B5*CD -ALP5*F5*C/R2*(F2*Y*SD+Q*B7) )     03850000
        DU( 9)= F3*X/R5*(-ALP4*B5*SD +ALP5*F5*C/R2*(D*B7*SD-Y*C7*CD) )  03860000
        DU(10)= F3/R5*   (-ALP4*D*A5*CD +ALP5*C*(A5*CD+D*QR5*A7) )      03870000
        DU(11)= F15*X/R7*( ALP4*Y*D*CD  +ALP5*C*(D*B7*SD-Y*C7*CD) )     03880000
        DU(12)= F15*X/R7*(-ALP4*Y*D*SD  +ALP5*C*(F2*D*CD-Q*C7) )        03890000
        DO 222 I=1,12                                                   03900000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        03910000
      ENDIF                                                             03920000
C===================================                                    03930000
C=====  DIP-SLIP CONTRIBUTION  =====                                    03940000
C===================================                                    03950000
      IF(POT2.NE.F0) THEN                                               03960000
        DU( 1)= ALP4*F3*X*T/R5          -ALP5*C*P*QRX                   03970000
        DU( 2)=-ALP4/R3*(C2D-F3*Y*T/R2) +ALP5*F3*C/R5*(S-Y*P*QR5)       03980000
        DU( 3)=-ALP4*A3/R3*SDCD         +ALP5*F3*C/R5*(T+D*P*QR5)       03990000
        DU( 4)= ALP4*F3*T/R5*A5              -ALP5*F5*C*P*QR/R2*A7      04000000
        DU( 5)= F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-ALP5*F5*C/R2*(S-Y*P*QR7)) 04010000
        DU( 6)= F3*X/R5*(ALP4*(F2+A5)*SDCD   -ALP5*F5*C/R2*(T+D*P*QR7)) 04020000
        DU( 7)= DU(5)                                                   04030000
        DU( 8)= F3/R5*(ALP4*(F2*Y*C2D+T*B5)                             04040000
     *                               +ALP5*C*(S2D-F10*Y*S/R2-P*QR5*B7)) 04050000
        DU( 9)= F3/R5*(ALP4*Y*A5*SDCD-ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7)) 04060000
        DU(10)= F3*X/R5*(-ALP4*(S2D-T*DR5) -ALP5*F5*C/R2*(T+D*P*QR7))   04070000
        DU(11)= F3/R5*(-ALP4*(D*B5*C2D+Y*C5*S2D)                        04080000
     *                                -ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7))04090000
        DU(12)= F3/R5*(-ALP4*D*A5*SDCD-ALP5*C*(S2D-F10*D*T/R2+P*QR5*C7))04100000
        DO 333 I=1,12                                                   04110000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        04120000
      ENDIF                                                             04130000
C========================================                               04140000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               04150000
C========================================                               04160000
      IF(POT3.NE.F0) THEN                                               04170000
        DU( 1)= F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))                     04180000
        DU( 2)= ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)04190000
        DU( 3)=-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)04200000
        DU( 4)=-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)          04210000
        DU( 5)= F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)                          04220000
     *                               -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z)) 04230000
        DU( 6)= F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)                        04240000
     *                               +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z)) 04250000
        DU( 7)= DU(5)                                                   04260000
        DU( 8)= F3/R5*(-ALP4*(F2*Y*S2D+S*B5)                            04270000
     *                -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5)) 04280000
        DU( 9)= F3/R5*( ALP4*Y*(F1-A5*SDSD)                             04290000
     *                +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))             04300000
        DU(10)= F3*X/R5*(-ALP4*(C2D+S*DR5)                              04310000
     *               +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))            04320000
        DU(11)= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)                        04330000
     *               +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))     04340000
        DU(12)= F3/R5*(-ALP4*D*(F1-A5*SDSD)                             04350000
     *               -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5))) 04360000
        DO 444 I=1,12                                                   04370000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        04380000
      ENDIF                                                             04390000
C=========================================                              04400000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              04410000
C=========================================                              04420000
      IF(POT4.NE.F0) THEN                                               04430000
        DU( 1)= ALP4*F3*X*D/R5                                          04440000
        DU( 2)= ALP4*F3*Y*D/R5                                          04450000
        DU( 3)= ALP4*C3/R3                                              04460000
        DU( 4)= ALP4*F3*D/R5*A5                                         04470000
        DU( 5)=-ALP4*F15*XY*D/R7                                        04480000
        DU( 6)=-ALP4*F3*X/R5*C5                                         04490000
        DU( 7)= DU(5)                                                   04500000
        DU( 8)= ALP4*F3*D/R5*B5                                         04510000
        DU( 9)=-ALP4*F3*Y/R5*C5                                         04520000
        DU(10)= DU(6)                                                   04530000
        DU(11)= DU(9)                                                   04540000
        DU(12)= ALP4*F3*D/R5*(F2+C5)                                    04550000
        DO 555 I=1,12                                                   04560000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        04570000
      ENDIF                                                             04580000
      RETURN                                                            04590000
      END                                                               04600000
      SUBROUTINE  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,                           04610005
     *              AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,                  04620005
     *              UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)  04630005
      IMPLICIT REAL*8 (A-H,O-Z)                                         04640005
      REAL*8   ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, 04650005
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             04660005
C                                                                       04670005
C********************************************************************   04680005
C*****                                                          *****   04690005
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   04700005
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   04710005
C*****              CODED BY  Y.OKADA ... SEP.1991              *****   04720005
C*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****   04730005
C*****                          JUL.1993                        *****   04740005
C********************************************************************   04750005
C                                                                       04760005
C***** INPUT                                                            04770005
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           04780005
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          04790005
C*****   DEPTH : DEPTH OF REFERENCE POINT                               04800005
C*****   DIP   : DIP-ANGLE (DEGREE)                                     04810005
C*****   AL1,AL2   : FAULT LENGTH RANGE                                 04820005
C*****   AW1,AW2   : FAULT WIDTH RANGE                                  04830005
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              04840005
C                                                                       04850005
C***** OUTPUT                                                           04860005
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               04870005
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             04880005
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )04890005
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     04900005
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  04910005
C                                                                       04920005
      COMMON /C0/DUMMY(5),SD,CD,DUMMY2(5)                               04930005
      DIMENSION  XI(2),ET(2),KXI(2),KET(2)                              04940005
      DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)                   04950005
      DATA  F0,EPS/ 0.D0, 1.D-6 /                                       04960005
C-----                                                                  04970005
      IF(Z.GT.0.) WRITE(6,'('' ** POSITIVE Z WAS GIVEN IN SUB-DC3D'')') 04980005
      DO 111 I=1,12                                                     04990005
        U  (I)=F0                                                       05000005
        DUA(I)=F0                                                       05010005
        DUB(I)=F0                                                       05020005
        DUC(I)=F0                                                       05030005
  111 CONTINUE                                                          05040005
      AALPHA=ALPHA                                                      05050005
      DDIP=DIP                                                          05060005
      CALL DCCON0(AALPHA,DDIP)                                          05070005
C-----                                                                  05080005
      ZZ=Z                                                              05090005
      DD1=DISL1                                                         05100005
      DD2=DISL2                                                         05110005
      DD3=DISL3                                                         05120005
      XI(1)=X-AL1                                                       05130005
      XI(2)=X-AL2                                                       05140005
      IF(DABS(XI(1)).LT.EPS) XI(1)=F0                                   05150005
      IF(DABS(XI(2)).LT.EPS) XI(2)=F0                                   05160005
C======================================                                 05170005
C=====  REAL-SOURCE CONTRIBUTION  =====                                 05180005
C======================================                                 05190005
      D=DEPTH+Z                                                         05200005
      P=Y*CD+D*SD                                                       05210005
      Q=Y*SD-D*CD                                                       05220005
      ET(1)=P-AW1                                                       05230005
      ET(2)=P-AW2                                                       05240005
      IF(DABS(Q).LT.EPS)  Q=F0                                          05250005
      IF(DABS(ET(1)).LT.EPS) ET(1)=F0                                   05260005
      IF(DABS(ET(2)).LT.EPS) ET(2)=F0                                   05270005
C--------------------------------                                       05280005
C----- REJECT SINGULAR CASE -----                                       05290005
C--------------------------------                                       05300005
C----- ON FAULT EDGE                                                    05310014
      IF(Q.EQ.F0                                                        05320014
     *   .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)          05330014
     *         .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))       05340014
     *   GO TO 99                                                       05350005
C----- ON NEGATIVE EXTENSION OF FAULT EDGE                              05360014
      KXI(1)=0                                                          05370005
      KXI(2)=0                                                          05380005
      KET(1)=0                                                          05390005
      KET(2)=0                                                          05400005
      R12=DSQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)                            05410005
      R21=DSQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)                            05420005
      R22=DSQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)                            05430005
      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   05440011
      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   05450011
      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   05460011
      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   05470011
C=====                                                                  05480015
      DO 223 K=1,2                                                      05490005
      DO 222 J=1,2                                                      05500005
        CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  05510014
        CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          05520005
C-----                                                                  05530005
        DO 220 I=1,10,3                                                 05540005
          DU(I)  =-DUA(I)                                               05550005
          DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD                              05560005
          DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD                              05570005
          IF(I.LT.10) GO TO 220                                         05580005
          DU(I)  =-DU(I)                                                05590005
          DU(I+1)=-DU(I+1)                                              05600005
          DU(I+2)=-DU(I+2)                                              05610005
  220   CONTINUE                                                        05620005
        DO 221 I=1,12                                                   05630005
          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  05640005
          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  05650005
  221   CONTINUE                                                        05660005
C-----                                                                  05670005
  222 CONTINUE                                                          05680005
  223 CONTINUE                                                          05690005
C=======================================                                05700005
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                05710005
C=======================================                                05720005
      D=DEPTH-Z                                                         05730005
      P=Y*CD+D*SD                                                       05740005
      Q=Y*SD-D*CD                                                       05750005
      ET(1)=P-AW1                                                       05760005
      ET(2)=P-AW2                                                       05770005
      IF(DABS(Q).LT.EPS)  Q=F0                                          05780005
      IF(DABS(ET(1)).LT.EPS) ET(1)=F0                                   05790005
      IF(DABS(ET(2)).LT.EPS) ET(2)=F0                                   05800005
C--------------------------------                                       05810005
C----- REJECT SINGULAR CASE -----                                       05820005
C--------------------------------                                       05830005
C----- ON FAULT EDGE                                                    05840015
      IF(Q.EQ.F0                                                        05850015
     *   .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)          05860015
     *         .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))       05870015
     *   GO TO 99                                                       05880015
C----- ON NEGATIVE EXTENSION OF FAULT EDGE                              05890015
      KXI(1)=0                                                          05900005
      KXI(2)=0                                                          05910005
      KET(1)=0                                                          05920005
      KET(2)=0                                                          05930005
      R12=DSQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)                            05940005
      R21=DSQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)                            05950005
      R22=DSQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)                            05960005
      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   05970011
      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   05980011
      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   05990011
      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   06000011
C=====                                                                  06010015
      DO 334 K=1,2                                                      06020005
      DO 333 J=1,2                                                      06030005
        CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  06040014
        CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          06050005
        CALL UB(XI(J),ET(K),Q,DD1,DD2,DD3,DUB)                          06060005
        CALL UC(XI(J),ET(K),Q,ZZ,DD1,DD2,DD3,DUC)                       06070005
C-----                                                                  06080005
        DO 330 I=1,10,3                                                 06090005
          DU(I)=DUA(I)+DUB(I)+Z*DUC(I)                                  06100005
          DU(I+1)=(DUA(I+1)+DUB(I+1)+Z*DUC(I+1))*CD                     06110005
     *           -(DUA(I+2)+DUB(I+2)+Z*DUC(I+2))*SD                     06120005
          DU(I+2)=(DUA(I+1)+DUB(I+1)-Z*DUC(I+1))*SD                     06130005
     *           +(DUA(I+2)+DUB(I+2)-Z*DUC(I+2))*CD                     06140005
          IF(I.LT.10) GO TO 330                                         06150005
          DU(10)=DU(10)+DUC(1)                                          06160005
          DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD                             06170005
          DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD                             06180005
  330   CONTINUE                                                        06190005
        DO 331 I=1,12                                                   06200005
          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  06210005
          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  06220005
  331   CONTINUE                                                        06230005
C-----                                                                  06240005
  333 CONTINUE                                                          06250005
  334 CONTINUE                                                          06260005
C=====                                                                  06270005
      UX=U(1)                                                           06280005
      UY=U(2)                                                           06290005
      UZ=U(3)                                                           06300005
      UXX=U(4)                                                          06310005
      UYX=U(5)                                                          06320005
      UZX=U(6)                                                          06330005
      UXY=U(7)                                                          06340005
      UYY=U(8)                                                          06350005
      UZY=U(9)                                                          06360005
      UXZ=U(10)                                                         06370005
      UYZ=U(11)                                                         06380005
      UZZ=U(12)                                                         06390005
      IRET=0                                                            06400005
      RETURN                                                            06410005
C===========================================                            06420005
C=====  IN CASE OF SINGULAR (ON EDGE)  =====                            06430005
C===========================================                            06440005
   99 UX=F0                                                             06450005
      UY=F0                                                             06460005
      UZ=F0                                                             06470005
      UXX=F0                                                            06480005
      UYX=F0                                                            06490005
      UZX=F0                                                            06500005
      UXY=F0                                                            06510005
      UYY=F0                                                            06520005
      UZY=F0                                                            06530005
      UXZ=F0                                                            06540005
      UYZ=F0                                                            06550005
      UZZ=F0                                                            06560005
      IRET=1                                                            06570005
      RETURN                                                            06580005
      END                                                               06590005
      SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U)                       06600005
      IMPLICIT REAL*8 (A-H,O-Z)                                         06610005
      DIMENSION U(12),DU(12)                                            06620005
C                                                                       06630005
C********************************************************************   06640005
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   06650005
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   06660005
C********************************************************************   06670005
C                                                                       06680005
C***** INPUT                                                            06690005
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  06700005
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              06710005
C***** OUTPUT                                                           06720005
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     06730005
C                                                                       06740005
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  06750005
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  06760005
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                06770005
      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/                     06780005
C-----                                                                  06790005
      DO 111  I=1,12                                                    06800005
  111 U(I)=F0                                                           06810005
      XY=XI*Y11                                                         06820005
      QX=Q *X11                                                         06830005
      QY=Q *Y11                                                         06840005
C======================================                                 06850005
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 06860005
C======================================                                 06870005
      IF(DISL1.NE.F0) THEN                                              06880005
        DU( 1)=    TT/F2 +ALP2*XI*QY                                    06890005
        DU( 2)=           ALP2*Q/R                                      06900005
        DU( 3)= ALP1*ALE -ALP2*Q*QY                                     06910005
        DU( 4)=-ALP1*QY  -ALP2*XI2*Q*Y32                                06920005
        DU( 5)=          -ALP2*XI*Q/R3                                  06930005
        DU( 6)= ALP1*XY  +ALP2*XI*Q2*Y32                                06940005
        DU( 7)= ALP1*XY*SD        +ALP2*XI*FY+D/F2*X11                  06950005
        DU( 8)=                    ALP2*EY                              06960005
        DU( 9)= ALP1*(CD/R+QY*SD) -ALP2*Q*FY                            06970005
        DU(10)= ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11                  06980005
        DU(11)=                    ALP2*EZ                              06990005
        DU(12)=-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ                            07000005
        DO 222 I=1,12                                                   07010005
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       07020005
      ENDIF                                                             07030005
C======================================                                 07040005
C=====    DIP-SLIP CONTRIBUTION   =====                                 07050005
C======================================                                 07060005
      IF(DISL2.NE.F0) THEN                                              07070005
        DU( 1)=           ALP2*Q/R                                      07080005
        DU( 2)=    TT/F2 +ALP2*ET*QX                                    07090005
        DU( 3)= ALP1*ALX -ALP2*Q*QX                                     07100005
        DU( 4)=        -ALP2*XI*Q/R3                                    07110005
        DU( 5)= -QY/F2 -ALP2*ET*Q/R3                                    07120005
        DU( 6)= ALP1/R +ALP2*Q2/R3                                      07130005
        DU( 7)=                      ALP2*EY                            07140005
        DU( 8)= ALP1*D*X11+XY/F2*SD +ALP2*ET*GY                         07150005
        DU( 9)= ALP1*Y*X11          -ALP2*Q*GY                          07160005
        DU(10)=                      ALP2*EZ                            07170005
        DU(11)= ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ                         07180005
        DU(12)=-ALP1*D*X11          -ALP2*Q*GZ                          07190005
        DO 333 I=1,12                                                   07200005
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       07210005
      ENDIF                                                             07220005
C========================================                               07230005
C=====  TENSILE-FAULT CONTRIBUTION  =====                               07240005
C========================================                               07250005
      IF(DISL3.NE.F0) THEN                                              07260005
        DU( 1)=-ALP1*ALE -ALP2*Q*QY                                     07270005
        DU( 2)=-ALP1*ALX -ALP2*Q*QX                                     07280005
        DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY)                            07290005
        DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32                                07300005
        DU( 5)=-ALP1/R   +ALP2*Q2/R3                                    07310005
        DU( 6)=-ALP1*QY  -ALP2*Q*Q2*Y32                                 07320005
        DU( 7)=-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY                           07330005
        DU( 8)=-ALP1*Y*X11         -ALP2*Q*GY                           07340005
        DU( 9)= ALP1*(D*X11+XY*SD) +ALP2*Q*HY                           07350005
        DU(10)= ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ                           07360005
        DU(11)= ALP1*D*X11         -ALP2*Q*GZ                           07370005
        DU(12)= ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ                           07380005
        DO 444 I=1,12                                                   07390005
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       07400005
      ENDIF                                                             07410005
      RETURN                                                            07420005
      END                                                               07430005
      SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U)                       07440005
      IMPLICIT REAL*8 (A-H,O-Z)                                         07450005
      DIMENSION U(12),DU(12)                                            07460005
C                                                                       07470005
C********************************************************************   07480005
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   07490005
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   07500005
C********************************************************************   07510005
C                                                                       07520005
C***** INPUT                                                            07530005
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  07540005
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              07550005
C***** OUTPUT                                                           07560005
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     07570005
C                                                                       07580005
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  07590005
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  07600005
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                07610005
      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            07620005
C-----                                                                  07630005
      RD=R+D                                                            07640005
      D11=F1/(R*RD)                                                     07650005
      AJ2=XI*Y/RD*D11                                                   07660005
      AJ5=-(D+Y*Y/RD)*D11                                               07670005
      IF(CD.NE.F0) THEN                                                 07680005
        IF(XI.EQ.F0) THEN                                               07690005
          AI4=F0                                                        07700005
        ELSE                                                            07710005
          X=DSQRT(XI2+Q2)                                               07720005
          AI4=F1/CDCD*( XI/RD*SDCD                                      07730005
     *       +F2*DATAN((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )        07740005
        ENDIF                                                           07750005
        AI3=(Y*CD/RD-ALE+SD*DLOG(RD))/CDCD                              07760005
        AK1=XI*(D11-Y11*SD)/CD                                          07770005
        AK3=(Q*Y11-Y*D11)/CD                                            07780005
        AJ3=(AK1-AJ2*SD)/CD                                             07790005
        AJ6=(AK3-AJ5*SD)/CD                                             07800005
      ELSE                                                              07810005
        RD2=RD*RD                                                       07820005
        AI3=(ET/RD+Y*Q/RD2-ALE)/F2                                      07830005
        AI4=XI*Y/RD2/F2                                                 07840005
        AK1=XI*Q/RD*D11                                                 07850005
        AK3=SD/RD*(XI2*D11-F1)                                          07860005
        AJ3=-XI/RD2*(Q2*D11-F1/F2)                                      07870005
        AJ6=-Y/RD2*(XI2*D11-F1/F2)                                      07880005
      ENDIF                                                             07890005
C-----                                                                  07900005
      XY=XI*Y11                                                         07910005
      AI1=-XI/RD*CD-AI4*SD                                              07920005
      AI2= DLOG(RD)+AI3*SD                                              07930005
      AK2= F1/R+AK3*SD                                                  07940005
      AK4= XY*CD-AK1*SD                                                 07950005
      AJ1= AJ5*CD-AJ6*SD                                                07960005
      AJ4=-XY-AJ2*CD+AJ3*SD                                             07970005
C=====                                                                  07980005
      DO 111  I=1,12                                                    07990005
  111 U(I)=F0                                                           08000005
      QX=Q*X11                                                          08010005
      QY=Q*Y11                                                          08020005
C======================================                                 08030005
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 08040005
C======================================                                 08050005
      IF(DISL1.NE.F0) THEN                                              08060005
        DU( 1)=-XI*QY-TT -ALP3*AI1*SD                                   08070005
        DU( 2)=-Q/R      +ALP3*Y/RD*SD                                  08080005
        DU( 3)= Q*QY     -ALP3*AI2*SD                                   08090005
        DU( 4)= XI2*Q*Y32 -ALP3*AJ1*SD                                  08100005
        DU( 5)= XI*Q/R3   -ALP3*AJ2*SD                                  08110005
        DU( 6)=-XI*Q2*Y32 -ALP3*AJ3*SD                                  08120005
        DU( 7)=-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD                           08130005
        DU( 8)=-EY          +ALP3*(F1/R+AJ5)*SD                         08140005
        DU( 9)= Q*FY        -ALP3*(QY-AJ6)*SD                           08150005
        DU(10)=-XI*FZ-Y*X11 +ALP3*AK1*SD                                08160005
        DU(11)=-EZ          +ALP3*Y*D11*SD                              08170005
        DU(12)= Q*FZ        +ALP3*AK2*SD                                08180005
        DO 222 I=1,12                                                   08190005
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       08200005
      ENDIF                                                             08210005
C======================================                                 08220005
C=====    DIP-SLIP CONTRIBUTION   =====                                 08230005
C======================================                                 08240005
      IF(DISL2.NE.F0) THEN                                              08250005
        DU( 1)=-Q/R      +ALP3*AI3*SDCD                                 08260005
        DU( 2)=-ET*QX-TT -ALP3*XI/RD*SDCD                               08270005
        DU( 3)= Q*QX     +ALP3*AI4*SDCD                                 08280005
        DU( 4)= XI*Q/R3     +ALP3*AJ4*SDCD                              08290005
        DU( 5)= ET*Q/R3+QY  +ALP3*AJ5*SDCD                              08300005
        DU( 6)=-Q2/R3       +ALP3*AJ6*SDCD                              08310005
        DU( 7)=-EY          +ALP3*AJ1*SDCD                              08320005
        DU( 8)=-ET*GY-XY*SD +ALP3*AJ2*SDCD                              08330005
        DU( 9)= Q*GY        +ALP3*AJ3*SDCD                              08340005
        DU(10)=-EZ          -ALP3*AK3*SDCD                              08350005
        DU(11)=-ET*GZ-XY*CD -ALP3*XI*D11*SDCD                           08360005
        DU(12)= Q*GZ        -ALP3*AK4*SDCD                              08370005
        DO 333 I=1,12                                                   08380005
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       08390005
      ENDIF                                                             08400005
C========================================                               08410005
C=====  TENSILE-FAULT CONTRIBUTION  =====                               08420005
C========================================                               08430005
      IF(DISL3.NE.F0) THEN                                              08440005
        DU( 1)= Q*QY           -ALP3*AI3*SDSD                           08450005
        DU( 2)= Q*QX           +ALP3*XI/RD*SDSD                         08460005
        DU( 3)= ET*QX+XI*QY-TT -ALP3*AI4*SDSD                           08470005
        DU( 4)=-XI*Q2*Y32 -ALP3*AJ4*SDSD                                08480005
        DU( 5)=-Q2/R3     -ALP3*AJ5*SDSD                                08490005
        DU( 6)= Q*Q2*Y32  -ALP3*AJ6*SDSD                                08500005
        DU( 7)= Q*FY -ALP3*AJ1*SDSD                                     08510005
        DU( 8)= Q*GY -ALP3*AJ2*SDSD                                     08520005
        DU( 9)=-Q*HY -ALP3*AJ3*SDSD                                     08530005
        DU(10)= Q*FZ +ALP3*AK3*SDSD                                     08540005
        DU(11)= Q*GZ +ALP3*XI*D11*SDSD                                  08550005
        DU(12)=-Q*HZ +ALP3*AK4*SDSD                                     08560005
        DO 444 I=1,12                                                   08570005
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       08580005
      ENDIF                                                             08590005
      RETURN                                                            08600005
      END                                                               08610005
      SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U)                     08620005
      IMPLICIT REAL*8 (A-H,O-Z)                                         08630005
      DIMENSION U(12),DU(12)                                            08640005
C                                                                       08650005
C********************************************************************   08660005
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   08670005
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   08680005
C********************************************************************   08690005
C                                                                       08700005
C***** INPUT                                                            08710005
C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              08720005
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              08730005
C***** OUTPUT                                                           08740005
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     08750005
C                                                                       08760005
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  08770005
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  08780005
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                08790005
      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/     08800005
C-----                                                                  08810005
      C=D+Z                                                             08820005
      X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2                     08830005
      Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2                     08840005
      H=Q*CD-Z                                                          08850005
      Z32=SD/R3-H*Y32                                                   08860005
      Z53=F3*SD/R5-H*Y53                                                08870005
      Y0=Y11-XI2*Y32                                                    08880005
      Z0=Z32-XI2*Z53                                                    08890005
      PPY=CD/R3+Q*Y32*SD                                                08900005
      PPZ=SD/R3-Q*Y32*CD                                                08910005
      QQ=Z*Y32+Z32+Z0                                                   08920005
      QQY=F3*C*D/R5-QQ*SD                                               08930005
      QQZ=F3*C*Y/R5-QQ*CD+Q*Y32                                         08940005
      XY=XI*Y11                                                         08950005
      QX=Q*X11                                                          08960005
      QY=Q*Y11                                                          08970005
      QR=F3*Q/R5                                                        08980005
      CQX=C*Q*X53                                                       08990005
      CDR=(C+D)/R3                                                      09000005
      YY0=Y/R3-Y0*CD                                                    09010005
C=====                                                                  09020005
      DO 111  I=1,12                                                    09030005
  111 U(I)=F0                                                           09040005
C======================================                                 09050005
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 09060005
C======================================                                 09070005
      IF(DISL1.NE.F0) THEN                                              09080005
        DU( 1)= ALP4*XY*CD           -ALP5*XI*Q*Z32                     09090005
        DU( 2)= ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3                       09100005
        DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)      09110005
        DU( 4)= ALP4*Y0*CD                  -ALP5*Q*Z0                  09120005
        DU( 5)=-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR               09130005
        DU( 6)=-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)    09140005
        DU( 7)=-ALP4*XI*PPY*CD    -ALP5*XI*QQY                          09150005
        DU( 8)= ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD                         09160005
     *                            -ALP5*(CDR*SD-ET/R3-C*Y*QR)           09170005
        DU( 9)=-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD) 09180005
        DU(10)= ALP4*XI*PPZ*CD    -ALP5*XI*QQZ                          09190005
        DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR)   09200005
        DU(12)=         YY0*CD    -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD) 09210005
        DO 222 I=1,12                                                   09220005
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       09230005
      ENDIF                                                             09240005
C======================================                                 09250005
C=====    DIP-SLIP CONTRIBUTION   =====                                 09260005
C======================================                                 09270005
      IF(DISL2.NE.F0) THEN                                              09280005
        DU( 1)= ALP4*CD/R -QY*SD -ALP5*C*Q/R3                           09290005
        DU( 2)= ALP4*Y*X11       -ALP5*C*ET*Q*X32                       09300005
        DU( 3)=     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32)                   09310005
        DU( 4)=-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD                09320005
        DU( 5)=-ALP4*Y/R3     +ALP5*C*ET*QR                             09330005
        DU( 6)=    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2)                  09340005
        DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR)                09350005
        DU( 8)= ALP4*(X11-Y*Y*X32) -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53) 09360005
        DU( 9)=  XI*PPY*SD+Y*D*X32 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)   09370005
        DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR)                09380005
        DU(11)= ALP4*Y*D*X32       -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53) 09390005
        DU(12)=-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53) 09400005
        DO 333 I=1,12                                                   09410005
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       09420005
      ENDIF                                                             09430005
C========================================                               09440005
C=====  TENSILE-FAULT CONTRIBUTION  =====                               09450005
C========================================                               09460005
      IF(DISL3.NE.F0) THEN                                              09470005
        DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32)                09480005
        DU( 2)= ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32)                09490005
        DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32)           09500005
        DU( 4)= ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)09510005
        DU( 5)= ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2)             09520005
        DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0)                 09530005
        DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)    09540005
        DU( 8)=-ALP4*F2*XI*PPY*SD-Y*D*X32                               09550005
     *                    +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)            09560005
        DU( 9)=-ALP4*(XI*PPY*CD-X11+Y*Y*X32)                            09570005
     *                    +ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY) 09580005
        DU(10)=  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)  09590005
        DU(11)= ALP4*F2*XI*PPZ*SD-X11+D*D*X32                           09600005
     *                    -ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)            09610005
        DU(12)= ALP4*(XI*PPZ*CD+Y*D*X32)                                09620005
     *                    +ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ) 09630005
        DO 444 I=1,12                                                   09640005
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       09650005
      ENDIF                                                             09660005
      RETURN                                                            09670005
      END                                                               09680005
      SUBROUTINE  DCCON0(ALPHA,DIP)                                     09690005
      IMPLICIT REAL*8 (A-H,O-Z)                                         09700005
C                                                                       09710005
C*******************************************************************    09720005
C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    09730005
C*******************************************************************    09740005
C                                                                       09750005
C***** INPUT                                                            09760005
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           09770005
C*****   DIP   : DIP-ANGLE (DEGREE)                                     09780005
C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   09790005
C                                                                       09800005
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  09810005
      DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             09820005
      DATA EPS/1.D-6/                                                   09830005
C-----                                                                  09840005
      ALP1=(F1-ALPHA)/F2                                                09850005
      ALP2= ALPHA/F2                                                    09860005
      ALP3=(F1-ALPHA)/ALPHA                                             09870005
      ALP4= F1-ALPHA                                                    09880005
      ALP5= ALPHA                                                       09890005
C-----                                                                  09900005
      P18=PI2/360.D0                                                    09910005
      SD=DSIN(DIP*P18)                                                  09920005
      CD=DCOS(DIP*P18)                                                  09930005
      IF(DABS(CD).LT.EPS) THEN                                          09940005
        CD=F0                                                           09950005
        IF(SD.GT.F0) SD= F1                                             09960005
        IF(SD.LT.F0) SD=-F1                                             09970005
      ENDIF                                                             09980005
      SDSD=SD*SD                                                        09990005
      CDCD=CD*CD                                                        10000005
      SDCD=SD*CD                                                        10010005
      S2D=F2*SDCD                                                       10020005
      C2D=CDCD-SDSD                                                     10030005
      RETURN                                                            10040005
      END                                                               10050005
      SUBROUTINE  DCCON1(X,Y,D)                                         10060005
      IMPLICIT REAL*8 (A-H,O-Z)                                         10070005
C                                                                       10080005
C********************************************************************** 10090005
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 10100005
C********************************************************************** 10110005
C                                                                       10120005
C***** INPUT                                                            10130005
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    10140005
C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  10150005
C                                                                       10160005
      COMMON /C0/DUMMY(5),SD,CD,DUMMY2(5)                               10170005
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     10180005
     *           UY,VY,WY,UZ,VZ,WZ                                      10190005
      DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/                  10200005
C-----                                                                  10210005
      IF(DABS(X).LT.EPS) X=F0                                           10220005
      IF(DABS(Y).LT.EPS) Y=F0                                           10230005
      IF(DABS(D).LT.EPS) D=F0                                           10240005
      P=Y*CD+D*SD                                                       10250005
      Q=Y*SD-D*CD                                                       10260005
      S=P*SD+Q*CD                                                       10270005
      T=P*CD-Q*SD                                                       10280005
      XY=X*Y                                                            10290005
      X2=X*X                                                            10300005
      Y2=Y*Y                                                            10310005
      D2=D*D                                                            10320005
      R2=X2+Y2+D2                                                       10330005
      R =DSQRT(R2)                                                      10340005
      IF(R.EQ.F0) RETURN                                                10350005
      R3=R *R2                                                          10360005
      R5=R3*R2                                                          10370005
      R7=R5*R2                                                          10380005
C-----                                                                  10390005
      A3=F1-F3*X2/R2                                                    10400005
      A5=F1-F5*X2/R2                                                    10410005
      B3=F1-F3*Y2/R2                                                    10420005
      C3=F1-F3*D2/R2                                                    10430005
C-----                                                                  10440005
      QR=F3*Q/R5                                                        10450005
      QRX=F5*QR*X/R2                                                    10460005
C-----                                                                  10470005
      UY=SD-F5*Y*Q/R2                                                   10480005
      UZ=CD+F5*D*Q/R2                                                   10490005
      VY=S -F5*Y*P*Q/R2                                                 10500005
      VZ=T +F5*D*P*Q/R2                                                 10510005
      WY=UY+SD                                                          10520005
      WZ=UZ+CD                                                          10530005
      RETURN                                                            10540005
      END                                                               10550005
      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD,KXI,KET)                         10560005
      IMPLICIT REAL*8 (A-H,O-Z)                                         10570005
C                                                                       10580005
C********************************************************************** 10590005
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 10600005
C********************************************************************** 10610005
C                                                                       10620005
C***** INPUT                                                            10630005
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  10640005
C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                10650005
C*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY  10660005
C                                                                       10670005
C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER010680005
C                                                                       10690005
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  10700005
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                10710005
      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/                          10720005
C-----                                                                  10730005
      IF(DABS(XI).LT.EPS) XI=F0                                         10740005
      IF(DABS(ET).LT.EPS) ET=F0                                         10750005
      IF(DABS( Q).LT.EPS)  Q=F0                                         10760005
      XI2=XI*XI                                                         10770005
      ET2=ET*ET                                                         10780005
      Q2=Q*Q                                                            10790005
      R2=XI2+ET2+Q2                                                     10800005
      R =DSQRT(R2)                                                      10810005
      IF(R.EQ.F0) RETURN                                                10820005
      R3=R *R2                                                          10830005
      R5=R3*R2                                                          10840005
      Y =ET*CD+Q*SD                                                     10850005
      D =ET*SD-Q*CD                                                     10860005
C-----                                                                  10870005
      IF(Q.EQ.F0) THEN                                                  10880005
        TT=F0                                                           10890005
      ELSE                                                              10900005
        TT=DATAN(XI*ET/(Q*R))                                           10910005
      ENDIF                                                             10920005
C-----                                                                  10930005
      IF(KXI.EQ.1) THEN                                                 10940005
        ALX=-DLOG(R-XI)                                                 10950005
        X11=F0                                                          10960005
        X32=F0                                                          10970005
      ELSE                                                              10980005
        RXI=R+XI                                                        10990005
        ALX=DLOG(RXI)                                                   11000005
        X11=F1/(R*RXI)                                                  11010005
        X32=(R+RXI)*X11*X11/R                                           11020005
      ENDIF                                                             11030005
C-----                                                                  11040005
      IF(KET.EQ.1) THEN                                                 11050005
        ALE=-DLOG(R-ET)                                                 11060005
        Y11=F0                                                          11070005
        Y32=F0                                                          11080005
      ELSE                                                              11090005
        RET=R+ET                                                        11100005
        ALE=DLOG(RET)                                                   11110005
        Y11=F1/(R*RET)                                                  11120005
        Y32=(R+RET)*Y11*Y11/R                                           11130005
      ENDIF                                                             11140005
C-----                                                                  11150005
      EY=SD/R-Y*Q/R3                                                    11160005
      EZ=CD/R+D*Q/R3                                                    11170005
      FY=D/R3+XI2*Y32*SD                                                11180005
      FZ=Y/R3+XI2*Y32*CD                                                11190005
      GY=F2*X11*SD-Y*Q*X32                                              11200005
      GZ=F2*X11*CD+D*Q*X32                                              11210005
      HY=D*Q*X32+XI*Q*Y32*SD                                            11220005
      HZ=Y*Q*X32+XI*Q*Y32*CD                                            11230005
      RETURN                                                            11240005
      END                                                               11250005
      SUBROUTINE  DCCON3(X,P,Q,AL1,AL2,AW1,AW2,KXI,KET)                 11260005
      IMPLICIT REAL*8 (A-H,O-Z)                                         11270005
      REAL*8  X,AL1,AL2,AW1,AW2                                         11280005
      DIMENSION KXI(2),KET(2)                                           11290005
C                                                                       11300005
C************************************************************           11310005
C*****   CHECK SINGULARITIES RELATED TO R+XI AND R+ET   *****           11320005
C************************************************************           11330005
C                                                                       11340005
C***** INPUT                                                            11350005
C*****   X,P,Q : STATION COORDINATE                                     11360005
C*****   AL1,AL2,AW1,AW2 : FAULT DIMENSIONS                             11370005
C***** OUTPUT                                                           11380005
C*****   KXI(2) : =1 IF R+XI<EPS                                        11390005
C*****   KET(2) : =1 IF R+ET<EPS                                        11400005
C                                                                       11410005
      DATA  EPS/1.D-6/                                                  11420005
C-----                                                                  11430005
      KXI(1)=0                                                          11440005
      KXI(2)=0                                                          11450005
      KET(1)=0                                                          11460005
      KET(2)=0                                                          11470005
      IF(X.GT.AL1 .AND. P.GT.AW2) RETURN                                11480005
C-----                                                                  11490005
      DO 222 K=1,2                                                      11500005
      IF(K.EQ.1) ET=P-AW1                                               11510005
      IF(K.EQ.2) ET=P-AW2                                               11520005
      DO 111 J=1,2                                                      11530005
        IF(J.EQ.1) XI=X-AL1                                             11540005
        IF(J.EQ.2) XI=X-AL2                                             11550005
        R=DSQRT(XI*XI+ET*ET+Q*Q)                                        11560005
        RXI=R+XI                                                        11570005
        RET=R+ET                                                        11580005
        IF(X.LT.AL1 .AND. RXI.LT.EPS) KXI(K)=1                          11590005
        IF(P.LT.AW1 .AND. RET.LT.EPS) KET(J)=1                          11600005
  111 CONTINUE                                                          11610005
  222 CONTINUE                                                          11620005
      RETURN                                                            11630005
      END                                                               11640005
