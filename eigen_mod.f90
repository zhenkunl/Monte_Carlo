!***********************************************************************
!  DESCRIPTION:
!       Module to compute all eigenvalues and corresponding eigenvector
!       numerically for symmetric matrix
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Shen Yu, SCC
!       Revised 12/2015 by Li Zhenkun, SCC
!
!***********************************************************************
    module eigen_mod
!
!      purpose:
!      find symmetric matrix all eigenvalue and its eigenvector numerically.
!      本程序的改编参阅国防科技大学出版的严宝勇等编著的《特征值特征向量库程序》
!      第二章求实对称阵的特征值和特征向量 P.96页
!
!****************************************************************************
       implicit none

       public  :: RS
       private :: TRED1, TRED2, IMTQL1, IMTQL2, TQL2, TQLRAT, SASUM, SDOT, ISMIN
!
       contains
!
!
       REAL FUNCTION SASUM(N,SX,INCX)
!         TAKES THE SUM OF THE ABSOLUTE VALUES.
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: N,INCX
          REAL(4),DIMENSION(*),INTENT(IN):: SX
          INTEGER:: I,IX,M,MP1
          REAL(4):: STEMP
!
          SASUM = 0.0E0
          STEMP = 0.0E0
          IF(N.LE.0) RETURN
          IF(INCX.EQ.1) GOTO 20
!         CODE FOR INCREMENT NOT EQUAL TO 1
          IX = 1
          IF(INCX.LT.0) IX = ( - N + 1) * INCX + 1
          DO 10 I = 1,N
              STEMP = STEMP + ABS(SX(IX))
              IX = IX + INCX
10        CONTINUE
          SASUM = STEMP
          RETURN
!         CODE FOR INCREMENT NOT EQUAL TO 1
!         CLEAN-UP LOOP
20        M = MOD(N,6)
          IF(M.EQ.0) GO TO 40
          DO 30 I = 1,M
              STEMP = STEMP + ABS(SX(I))
30        CONTINUE
          IF(N.LT.6) GO TO 60
40        MP1 = M + 1
          DO 50 I = MP1,N,6
              STEMP = STEMP + ABS(SX(I)) + ABS(SX(I + 1)) &
                      + ABS(SX(I + 2)) + ABS(SX(I + 3)) &
                      + ABS(SX(I + 4)) + ABS(SX(I + 5))
50        CONTINUE
60        SASUM = STEMP
          RETURN
       END FUNCTION SASUM
!
!
       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
!      FORMS THE DOT PRODUCT OF TWO VECTORS.
!      USES UNROLLED LOOPS FOR
!      INCREMENTS EQUAL TO ONE.
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: INCX,INCY,N
          REAL(4),DIMENSION(*),INTENT(IN):: SX,SY
          INTEGER:: I,IX,IY,M,MP1
          REAL(4):: STEMP
!
          STEMP = 0.0E0
          SDOT = 0.0E0
          IF(N.LE.0) RETURN
          IF(INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
!         CODE FOR UNEQUAL INCREMENTS OR
!         EQUAL INCREMENTS NOT EQUAL TO 1
          IX = 1
          IY = 1
          IF(INCX.LT.0) IX = ( - N + 1) * INCX + 1
          IF(INCY.LT.0) IX = ( - N + 1) * INCY + 1
          DO 10 I = 1,N
              STEMP = STEMP + SX(IX) * SY(IY)
              IX = IX + INCX
              IY = IY + INCY
10        CONTINUE
          SDOT=STEMP
          RETURN
!         CODE FOR BOTH INCREMENTS EQUAL TO 1
!         CLEAN-UP LOOP
20        M = MOD(N,5)
          IF(M.EQ.0) GOTO 40
          DO 30 I = 1,M
              STEMP = STEMP + SX(I) * SY(I)
30        CONTINUE
          IF(N.LT.5) GOTO 60
40        MP1 = M + 1
          DO 50 I = MP1,N,5
              STEMP = STEMP + SX(I) * SY(I) + SX(I + 1) &
                      * SY(I + 1) + SX(I + 2) * SY(I + 2) + SX(I + 3) &
                      * SY(I + 3) + SX(I + 4) * SY(I + 4)
50        CONTINUE
60        SDOT = STEMP
          RETURN
       END FUNCTION SDOT
!
!
       INTEGER FUNCTION ISMIN(N,SX,INCX)
!      FINDS THE INDEX OF ELEMENT HAVING MIN VALUE.
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: N,INCX
          REAL(4),DIMENSION(*),INTENT(IN):: SX
          INTEGER:: I,IX
          REAL(4):: SMIN
!
          ISMIN = -1
          IF(N.LT.1) RETURN
          ISMIN = 0
          IF(N.EQ.1) RETURN
          IF(INCX.EQ.1) GO TO 20
!         CODE FOR INCREMENT NOT EQUAL TO 1
          IX = 1
          IF(INCX.LT.0) IX = ( - N + 1) * INCX + 1
          SMIN = SX(IX)
          ISMIN = IX - 1
          IF(INCX.LT.0) GOTO 15
          DO 5 I = 2,N
              IX = IX + INCX
              IF(SX(IX).GE.SMIN) GOTO 5
              ISMIN = IX - 1
              SMIN = SX(IX)
5         CONTINUE
          RETURN
15        DO 10 I = 2,N
               IX = IX + INCX
               IF(SX(IX).GT.SMIN) GOTO 10
               ISMIN = IX - 1
               SMIN = SX(IX)
10        CONTINUE
          RETURN
!         CODE FOR INCREMENT EQUAL TO 1
20        SMIN = SX(1)
          DO 30 I = 2,N
              IF(SX(I).GE.SMIN) GOTO 30
                ISMIN = I - 1
                SMIN = SX(I)
30        CONTINUE
          RETURN
       END FUNCTION ISMIN
!
!
       SUBROUTINE TRED1(NM,N,AA,D,E,E2)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: NM,N
          REAL(4),DIMENSION(NM,N),INTENT(IN):: AA
          REAL(4),DIMENSION(*),INTENT(OUT):: D,E,E2
          REAL(4),DIMENSION(NM,N):: A
          INTEGER:: I,J,K,L,II,JP1
          REAL(4):: F,G,H,SCALE
!
          A = AA  ! avoiding the side effects of the sub-routine.
          DO 100 I = 1,N
100           D(I) = A(I,I)
          DO 300 II = 1,N
              I = N + 1 - II
              L = I - 1
              H = 0.
              SCALE = 0.
              IF(L.LT.1) GOTO 130
              SCALE = SASUM(L,A(I,1),NM)
              IF(SCALE.NE.0.0) GOTO 140
130           E(I) = 0.
              E2(I) = 0.
              GOTO 290
140           DO 150 K = 1,L
                  A(I,K) = A(I,K) / SCALE
150           CONTINUE
              H = SDOT(L,A(I,1),NM,A(I,1),NM)
              E2(I) = SCALE * SCALE * H
              F = A(I,L)
              G = -SIGN(SQRT(H),F)
              E(I) = SCALE*G
              H = H - F * G
              A(I,L) = F - G
              IF(L.EQ.1) GOTO 270
              F = 0.
              DO 240 J = 1,L
!                 * * FORM ELEMENT OF A*U * *
                  G = SDOT(J,A(J,1),NM,A(I,1),NM)
                  JP1 = J + 1
                  IF(L.LT.JP1) GOTO 220
                  G = G + SDOT(L-J,A(JP1,J),1,A(I,JP1),NM)
!                 * * FORM ELEMENT OF P * *
220               E(J) = G / H
                  F = F + E(J) * A(I,J)
240          CONTINUE
              H =F / (H + H)
!             * * FORM REDUCED A * *
              DO 260 J = 1,L
                  F = A(I,J)
                  G = E(J) - H*F
                  E(J) = G
                  DO 260 K = 1,J
                      A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
260           CONTINUE
270           DO 280 K = 1,L
280               A(I,K) = SCALE * A(I,K)
290           H = D(I)
              D(I) = A(I,I)
              A(I,I) = H
300       CONTINUE
          RETURN
       END SUBROUTINE TRED1
!
!
       SUBROUTINE TRED2(NM,N,AA,D,E,Z)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: N,NM
          REAL(4),DIMENSION(NM,N),INTENT(IN):: AA
          REAL(4),DIMENSION(NM,N),INTENT(OUT):: Z
          REAL(4),DIMENSION(*),INTENT(OUT):: D,E
          REAL(4),DIMENSION(NM,N):: A
          INTEGER:: I,J,K,L,II,JP1
          REAL(4):: F,G,H,HH,SCALE
!
          A = AA  ! avoiding the side effects of the sub-routine.
          DO 100 I = 1,N
              DO 100 J = 1,I
100               Z(I,J) = A(I,J)
          IF(N.EQ.1) GO TO 320
!           * * FOR I = N STEP -1 UNTIL 2 DO * *
          DO 300 II = 2,N
              I = N + 2 - II
              L = I - 1
              H = 0.0
              SCALE = 0.0
              IF(L.LT.2) GO TO 130
              SCALE = SASUM(L,Z(I,1),NM)
              IF(SCALE.NE.0.0) GO TO 140
130           E(I) = Z(I,L)
              GO TO 290
140           DO 141 J = 1,L
141               Z(I,J) = Z(I,J) / SCALE
              H = SDOT(L,Z(I,1),NM,Z(I,1),NM)
              F = Z(I,L)
              G = -SIGN(SQRT(H),F)
              E(I) = SCALE * G
              H = H - F * G
              Z(I,L) = F - G
              F = 0.0
              DO 240 J = 1,L
                  Z(J,I) = Z(I,J) / H
!               * * FORM ELEMENT OF A * U * *
                  G = SDOT(J,Z(J,1),NM,Z(I,1),NM)
                  JP1 = J + 1
                  IF(L.LT.JP1) GO TO 220
                  G = G + SDOT(L-J,Z(JP1,J),1,Z(I,JP1),NM)
!                 * * FORM ELEMENT OF P * *
220               E(J) = G / H
                  F = F + E(J) * Z(I,J)
240           CONTINUE
              HH = F / (H + H)
!             * * FORM REDUCED A * *
              DO 260 J = 1,L
                  F = Z(I,J)
                  G = E(J) - HH * F
                  E(J) = G
                  DO 260 K = 1,J
                      Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
260           CONTINUE
290           D(I) = H
300       CONTINUE
320       D(1) = 0.0
          E(1) = 0.0
!
!         * * ACCUMULATION OF TRANSFORMATION MATRICES * *
!
          DO 500 I = 1,N
              L = I - 1
              IF(D(I).EQ.0.0) GO TO 380
              IF(L.EQ.0) GO TO 380
              DO 360 J = 1,L
                  G = SDOT(L,Z(I,1),NM,Z(1,J),1)
                  DO 360 K = 1,L
                    Z(K,J) = Z(K,J) - G * Z(K,I)
360               CONTINUE
380           D(I) = Z(I,I)
              Z(I,I) = 1.0
              IF(L.LT.1) GO TO 500
              DO 400 J = 1,L
                  Z(I,J) = 0.0
                  Z(J,I) = 0.0
400           CONTINUE
500       CONTINUE
          RETURN
       END SUBROUTINE TRED2
!
!
       SUBROUTINE IMTQL1(N,D,EE,IERR)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: N
          REAL(4),DIMENSION(N),INTENT(INOUT):: D
          REAL(4),DIMENSION(N),INTENT(IN):: EE
          INTEGER,INTENT(OUT):: IERR
          REAL(4),DIMENSION(N):: E
          INTEGER:: I,J,L,M,II,MML
          REAL(4):: B,C,F,G,P,R,S,MACHEP
!         *** MACHEP IS A MACHINE DEPENDENP
!         PARAMETER SPECIFYING THE RELATIVE
!         PRECISION OF FLOATING POINT ARITHMETI! ***
!
          E = EE  ! avoiding the side effects of the sub-routine.
          MACHEP = 2. ** ( - 47)
          IERR = 0
          IF(N.EQ.1) GO TO 1001
          DO 100 I = 2,N
100           E(I - 1) = E(I)
          E(N) = 0.0
          DO 290 L = 1,N
              J = 0
!             *** LOOK FOR SMALL SUB-DIAGONAL ELEMENT ***
105           DO 110 M = L,N - 1
                  IF(ABS(E(M)).LE.MACHEP*(ABS(D(M)) &
                  + ABS(D(M + 1)))) GO TO 120
110           CONTINUE
              M = N
120           P = D(L)
              IF(M.EQ.L) GO TO 215
              IF(J.EQ.30) GO TO 1000
              J = J + 1
!           *** FORM SHIFT ***
              G = (D(L + 1) - P) / (2.0 * E(L))
              R = SQRT(G * G + 1.0)
              G = D(M) - P + E(L) / (G + SIGN(R,G))
              S = 1.0
              C = 1.0
              P = 0.0
              MML = M - L
!             *** FOR I=M-1 STEP -1 UNTIL L DO -- ***
              DO 200 II = 1,MML
                  I = M - II
                  F = S * E(I)
                  B = C * E(I)
                  IF(ABS(F).LT.ABS(G)) GO TO 150
                  C = G / F
                  R = SQRT(C * C + 1.0)
                  E(I + 1) = F * R
                  S = 1.0 / R
                  C = C * S
                  GO TO 160
150               S = F / G
                  R = SQRT(S * S + 1.0)
                  E(I + 1) = G * R
                  C = 1.0 / R
                  S = S * C
160               G = D(I + 1) - P
                  R = (D(I) - G) * S + 2.0 * C * B
                  P = S * R
                  D(I + 1) = G + P
                  G = C * R - B
200           CONTINUE
              D(L) = D(L) - P
              E(L) = G
              E(M) = 0.0
              GO TO 105
!             *** ORDER EIGENVALUES ***
215           IF(L.EQ.1) GO TO 250
!             *** FOR I=LSTEP-1 UNTIL 2 DO ***
              DO 230 II = 2,L
                  I = L + 2 - II
                  IF(P.GE.D(I - 1)) GO TO 270
                  D(I) = D(I - 1)
230           CONTINUE
250           I = 1
270           D(I) = P
290       CONTINUE
          GO TO 1001
!         *** SET ERROR--NO CONVERGENCE TO AN
!       EIGENVALUE AFTER 30 ITERATIONS ***
1000      IERR = L
1001      RETURN
       END SUBROUTINE IMTQL1
!
!
       SUBROUTINE IMTQL2(NM,N,D,EE,Z,IERR)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: NM,N
          REAL(4),DIMENSION(*),INTENT(INOUT):: D
          REAL(4),DIMENSION(N),INTENT(IN):: EE
          REAL(4),DIMENSION(NM,*),INTENT(INOUT):: Z
          INTEGER,INTENT(OUT):: IERR
          REAL(4),DIMENSION(N):: E
          INTEGER:: I,J,K,L,M,II,MML
          REAL(4):: B,C,F,G,P,R,S,MACHEP
!
          E = EE  ! avoiding the side effects of the sub-routine.
          MACHEP = 2. ** ( - 47)
          IERR = 0
          IF(N.EQ.1) GO TO 1001
          DO 100 I = 2,N
100           E(I - 1) = E(I)
          E(N) = 0.
          DO 240 L = 1,N
              J = 0
!            ** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **
105           DO 110 M = L,N - 1
                  IF(ABS(E(M)).LE.MACHEP*(ABS(D(M)) + &
                  ABS(D(M+1)))) GO TO 120
110           CONTINUE
              M = N
120           P = D(L)
              IF(M.EQ.L) GO TO 240
              IF(J.EQ.30) GO TO 1000
              J = J + 1
!           ************** FORM SHIFT ***************
              G = (D( L+ 1) - P) / (2. * E(L))
              R = SQRT(G * G + 1.0)
              G = D(M) - P + E(L) / (G + SIGN(R,G))
              S = 1.
              C = 1.
              P = 0.
              MML = M - L
!             **** FOR I = M - 1 STEP - 1 UNTIL L DO -- ****
              DO 200 II = 1,MML
                  I = M - II
                  F = S * E(I)
                  B = C * E(I)
                  IF(ABS(F).LT.ABS(G)) GO TO 150
                  C = G / F
                  R = SQRT(C * C + 1.)
                  E(I + 1) = F * R
                  S = 1. / R
                  C = C * S
                  GO TO 160
150               S = F / G
                  R = SQRT(S * S + 1.)
                  E(I + 1) = G * R
                  C = 1. / R
                  S = S * C
160               G = D(I + 1) - P
                  R = (D(I) - G) * S + 2. * C * B
                  P = S * R
                  D(I + 1) = G + P
                  G = C * R - B
!                 **********FPRM VECTOR***********
                  DO 180 K = 1,N
                      F = Z(K,I + 1)
                      Z(K,I + 1) = S * Z(K,I) + C*F
                      Z(K,I) = C * Z(K,I) - S * F
180               CONTINUE
200           CONTINUE
              D(L) = D(L) - P
              E(L) = G
              E(M) = 0.
              GO TO 105
240       CONTINUE
!         * * ORDER EIGENVALUES AND EIGENVECTORS * *
          DO 300 II = 2,N
              I = II - 1
              K = ISMIN(N - II + 2,D(I),1) + I
              P = D(K)
              IF(K.EQ.I) GO TO 300
              D(K) = D(I)
              D(I) = P
              DO 280 J = 1,N
                  P = Z(J,I)
                  Z(J,I) = Z(J,K)
                  Z(J,K) = P
280           CONTINUE
300       CONTINUE
          GO TO 1001
!         *** SET ERROR--NO CONVERGENCE TO AN
!         EIGENVALUE AFTER 30 ITERATIONS ***
1000      IERR = L
1001      RETURN
       END SUBROUTINE IMTQL2
!
!
       SUBROUTINE TQL2(NM,N,D,EE,Z,IERR)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: NM,N
          REAL(4),DIMENSION(N),INTENT(INOUT):: D
          REAL(4),DIMENSION(N),INTENT(IN):: EE
          REAL(4),DIMENSION(NM,N),INTENT(INOUT):: Z
          INTEGER,INTENT(OUT):: IERR
          REAL(4),DIMENSION(N):: E
          INTEGER:: I,J,K,L,M,II,L1,MML
          REAL(4):: C,F,B,H,HH,G,P,PP,R,S,MACHEP
!         * * * * * MACHEP IS A MACHINE DEPENDENT
!         PARAMETER SPECIFYING
!         THE RELATIVE PRECISION OF FLOATING
!         POINT ARITHMETIC. * * * * *
!
          E = EE  ! avoiding the side effects of the sub-routine.
          MACHEP = 2. ** ( - 47)
          IERR = 0
          IF(N.EQ.1) GOTO 1001
          DO 100 I = 2,N
100           E(I - 1) = E(I)
          F = 0.
          B = 0.
          E(N) = 0.
          DO 240 L = 1,N
              J = 0
              H = MACHEP * (ABS(D(L)) + ABS(E(L)))
              IF(B.LT.H) B = H
!             ** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **
              DO 110 M = L,N
                  IF(ABS(E(M)).LE.B) GOTO 120
!                 ** E(N) IS ALWAYS ZERO. SO THREE IS NO EXIT
!                 THROUGH THE BOTTOM OF THE LOOP **
110           CONTINUE
120           IF(M.EQ.L) GOTO 220
130           IF(J.EQ.30) GOTO 1000
              J = J + 1
!           * * * * FORM SHIFT * * * * *
              L1 = L + 1
              G = D(L)
              P = (D(L1) - G)/(2.0 * E(L))
              R = SQRT(P * P + 1.)
              D(L) = E(L) / (P + SIGN(R,P))
              H = G - D(L)
              DO 140 I = L1,N
140               D(I) = D(I) - H
              F = F + H
!             * * QL TRANSFORMATION * *
              P = D(M)
              C = 1.0
              S = 0.
              MML = M - L
!             * * FOR I = M - 1 STEP -1 UNTIL L DO- * *
              DO 200 II = 1,MML
                  I = M - II
                  G = C * E(I)
                  H = C * P
                  IF(ABS(P).LT.ABS(E(I))) GOTO 150
                  C = E(I) / P
                  R = SQRT(C * C + 1.0)
                  E(I + 1) = S * P * R
                  S = C / R
                  C = 1.0 / R
                  GOTO 160
150               C = P / E(I)
                  R = SQRT(C * C + 1.0)
                  E(I + 1) = S * E(I) * R
                  S = 1.0 / R
                  C = C * S
160               P = C * D(I) - S * G
                  D(I+1)=H+S*(C*G+S*D(I))
!                 * * * * * FORM VECTOR * * * * *
                  DO 180 K = 1,N
                      HH = Z(K,I + 1)
                      Z(K,I + 1) = S * Z(K,I) + C * HH
                      Z(K,I) = C * Z(K,I) - S * HH
180               CONTINUE
200           CONTINUE
              E(L) = S * P
              D(L) = C * P
              IF(ABS(E(L)).GT.B) GOTO 130
220           D(L) = D(L) + F
240       CONTINUE
!         * * ORDER EIGENVALUES AND EIGENVECTORS * *
          DO 300 II = 2,N
              I = II - 1
              K = ISMIN(N - II + 2,D(I),1) + I
              P = D(K)
              IF(K.EQ.I) GOTO 300
              D(K) = D(I)
              D(I) = P
              DO 280 J = 1,N
                  PP = Z(J,I)
                  Z(J,I) = Z(J,K)
                  Z(J,K) = PP
280           CONTINUE
300       CONTINUE
          GOTO 1001
!         * * SET ERROR--NO CONVERGENCE TO AN
!         EOGENVALUE AFTER 30 ITERATION * *
1000      IERR = L
1001      RETURN
       END SUBROUTINE TQL2
!
!
       SUBROUTINE TQLRAT(N,D,EE2,IERR)
!
          IMPLICIT NONE
          INTEGER,INTENT(IN):: N
          REAL(4),DIMENSION(N),INTENT(IN):: EE2
          REAL(4),DIMENSION(N),INTENT(INOUT):: D
          INTEGER,INTENT(OUT):: IERR
          REAL(4),DIMENSION(N):: E2
          INTEGER:: I,J,L,M,II,L1,MML
          REAL(4):: B,C,F,G,H,P,R,S,MACHEP
!         ***** MACHELP IS A MACHINE DEPLENDENT
!         PARAMETER SPECIFYING
!         THE RELATIVE PRECISION OF FLOATING
!         POINT ARITHMETI! *****
!
          E2 = EE2  ! avoiding the side effects of the sub-routine.
          MACHEP = 2. ** ( - 47)
          IERR = 0
          IF(N.EQ.1) GOTO 1001
          DO 100 I = 2,N
100           E2(I - 1) = E2(I)
          F = 0.
          B = 0.
          E2(N) = 0.
          DO 290 L = 1,N
              J = 0
              H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
              IF(B.GT.H) GOTO 105
              B = H
              C = B * B
!             ***** LOOK FOR SMALL SQUARED
!             SUB-DIAGONAL ELEMENT *****
105           DO 110 M = L,N
                  IF(E2(M).LE.C) GOTO 120
!                 ** E2(N) IS ALWAYS ZERO. SO THERE IS NOT
!                 EXIT THROUGH THE BOTTOM OF THE LOOP **
110           CONTINUE
120           IF(M.EQ.L) GOTO 210
130           IF(J.EQ.30) GOTO 1000
              J = J + J
!             ***** FORM SHIFT *****
              L1 = L + 1
              S = SQRT(E2(L))
              G = D(L)
              P = (D(L1) - G) / (2.0 * S)
              R = SQRT(P * P + 1.)
              D(L) = S / (P + SIGN(R,P))
              H = G - D(L)
              DO 140 I = L1,N
140               D(I) = D(I) - H
              F = F + H
!             ** RATIONAL QL TRANSFORMATION **
              G = D(M)
              IF(G.EQ.0.) G = B
              H = G
              S = 0.
              MML = M - L
!             ** FOR I=M-1 STEP -1 UNTIL L DO- **
              DO 200 II=1,MML
                  I = M - II
                  P = G * H
                  R = P + E2(I)
                  E2(I + 1) = S * R
                  S = E2(I) / R
                  D(I + 1) = H + S * (H + D(I))
                  G = D(I) - E2(I) / G
                  IF(G.EQ.0.) G = B
                  H = G * P / R
200           CONTINUE
              E2(L) = S * G
              D(L) = H
!             ** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **
              IF(H.EQ.0.) GOTO 210
              IF(ABS(E2(L)).LE.ABS(C/H)) GOTO 210
              E2(L) = H * E2(L)
              IF(E2(L).NE.0.) GOTO 130
210           P = D(L) + F
!             ***** ORDER EIGENVALUES *****
              IF(L.EQ.1) GOTO 250
!             ***** FOR I=L STEP -1 UNTIL 2 DO -*****
              DO 230 II = 2,L
                  I = L + 2 - II
                  IF(P.GE.D(I - 1)) GOTO 270
                  D(I) = D(I - 1)
230           CONTINUE
250           I = 1
270           D(I) = P
290       CONTINUE
          GOTO 1001
!         ** SET ERROR-NO CONVERGENCE TO AN
!         EIGENVALUE AFTER 30 ITERATIONS **
1000      IERR = L
1001      RETURN
       END SUBROUTINE TQLRAT
!
!
       SUBROUTINE RS(NM,N,AA,W,OPT,Z,IERR)
!
!         all eigenvalue & its eigenvector stored in array w and z are in
!         ascending order.
          IMPLICIT NONE
          INTEGER,INTENT(IN):: NM,N,OPT
!         NM is the number of row of matrix A.
!         N is the number of column of matrix A, and NM >= N.
!         parameter OPT is set to use different method for finding
!         eigenvalues its eigenvectors of symmetric matrix a.
!         OPT=1 for eigenvalue only by method 1.
!         OPT=2 for eigenvalue only by method 2.
!         OPT=3 for both eigenvalue and eigenvector by method 1.
!         OPT=4 for both eigenvalue and eigenvector by method 2.
          REAL(4),DIMENSION(NM,N),INTENT(IN):: AA
          INTEGER,INTENT(OUT):: IERR
!         All the eigenvalues of symmetric matrix a are found within 30
!         iterations, parameter IERR = 0, otherwise IERR will
!         give the eigenvalue subscript of evaluation failed.
!         if N > NM, sub RS will stop and return IERR = 10*N.
          REAL(4),DIMENSION(NM,N),INTENT(OUT):: Z
          REAL(4),DIMENSION(N),INTENT(OUT):: W
          REAL(4):: FV1(N),FV2(N)

          IF(NM >= N)THEN
              IF(OPT == 1)THEN
!                 FIND EIGENVALUES ONLY BY METHOD 1.
                  CALL TRED1(NM,N,AA,W,FV1,FV2)
                  CALL IMTQL1(N,W,FV1,IERR)
              ELSEIF(OPT == 2)THEN
!                 FIND EIGENVALUES ONLY BY METHOD 2.
                  CALL TRED1(NM,N,AA,W,FV1,FV2)
                  CALL TQLRAT(N,W,FV2,IERR)
              ELSEIF(OPT == 3)THEN
!                 FIND BOTH EIGENVALUES AND EIGENVECTORS BY METHOD 1.
                  CALL TRED2(NM,N,AA,W,FV1,Z)
                  CALL TQL2(NM,N,W,FV1,Z,IERR)
              ELSEIF(OPT == 4)THEN
!                 FIND BOTH EIGENVALUES AND EIGENVECTORS BY METHOD 2.
                  CALL TRED2(NM,N,AA,W,FV1,Z)
                  CALL IMTQL2(NM,N,W,FV1,Z,IERR)
              ELSE
                  WRITE(6, *)'OPT IS NOT SUPPORT, SUBROUTINE RS INPUT PARA ERROR (OPT) IS FOUND'
                  STOP
              ENDIF
          ELSE
              IERR = 10 * N
              WRITE(6, *)'IERR = 10 * N, SUBROUTINE RS INPUT PARA ERROR (N > NM) IS FOUND!'
              STOP
          ENDIF

          RETURN
       END SUBROUTINE RS

!
    end module eigen_mod
