*****************************************************
* 以下の行をコメントアウトしました。
* 必要に応じて復活させてください。
* (190行め前後にあります。)
*
*      WRITE (*,2003)
* 2003 FORMAT (' (SUBR.DEHINT) insufficient mesh',
*     $        ' refinement (warning only)')
******************************************************
*
      real*8 FUNCTION FUNC (x)
        implicit real*8(a-h,o-z)
        COMMON / COMETA / eta
        FUNC = x / (x + eta) * exp(- x ** 2)
     1           * CELI1(4 * eta * x / (eta + x)**2)
*       FUNC = exp(-x)
        RETURN
      END
*
      SUBROUTINE DEHINT (eta_, A, EPS, V)
*****************************************************
*       Integrate FUNC over (A,infinity) by         *
*     the double exponential formula (DE formula)   *
*               in double precision                 *
*      COPYRIGHT : M.Mori  JUNE 30 1989  V.1        *
*                                                   *
*     Before the first call of this routine         *
*     HIAB must be called                           *
*                                                   *
*   ---- input parameters ----                      *
*     FUNC...name of function subprogram for        *
*            integrand                              *
*     A...lower bound of integration                *
*     EPS...absolute error tolerance                *
*   ---- output parameter ----                      *
*     V...integral of FUNC                          *
*****************************************************
      real*8 FUNC,A,EPS,V,eta_,eta
      real*8 ONE,HALF
      real*8 H
      real*8 VOLD,VNEW,WM,WP
      real*8 EPS0,EPSQ,EPSM,EPSP
      real*8 AM,A0,AP,BM,B0,BP
*
      COMMON / COMDEF / AM(640,3),A0(3),AP(640,3),
     $                  BM(640,3),B0(3),BP(640,3)
      COMMON / COMDEN / NPOW,NEND
      COMMON / COMCNT / NEVAL,L
*
      PARAMETER (ONE = 1)
      PARAMETER (HALF = ONE / 2)
*
      DATA EPS0 / 1.0D-32 /
      DATA EPSM / 0 /, EPSP / 0 /
      COMMON / COMETA / eta; eta=eta_
*
      NEVAL = 0
*
      IF (ABS(EPS) .GE. EPS0) THEN
        EPSV = ABS(EPS)
      ELSE
        EPSV = EPS0
      END IF
*
      EPSQ = 0.2 * SQRT(EPSV)
*
      H = HALF
*
      IS = 2**NPOW
      IH = IS
*
      L = 1
*
  101 CONTINUE
*
      KM = 0
      KP = 0
      NM = 0
      NP = 0
*
      VNEW = 0
*
*     ---- initial step ----
*          integrate with mesh size = 0.5
*          and check decay of integrand
*
      DO 110 I = IS, NEND, IH
*
        IF (KP .LE. 1) THEN
          WP = FUNC(AP(I,L) + A) * BP(I,L)
          NEVAL = NEVAL + 1
          VNEW = VNEW + WP
          IF (ABS(WP) .LE. EPSV) THEN
            KP = KP + 1
            IF (KP .GE. 2) THEN
              NP = I - IH
              GO TO 111
            END IF
          ELSE
            KP = 0
          END IF
        END IF
  110 CONTINUE
*
  111 CONTINUE
*
      IF (L .LE. 2) THEN
        IF (NP .EQ. 0) THEN
          L = L + 1
          GO TO 101
        END IF
      END IF
*
      DO 120 I = IS, NEND, IH
*
        IF (KM .LE. 1) THEN
          WM = FUNC(AM(I,L) + A) * BM(I,L)
          NEVAL = NEVAL + 1
          VNEW = VNEW + WM
          IF (ABS(WM) .LE. EPSV) THEN
            KM = KM + 1
            IF (KM .GE. 2) THEN
              NM = I - IH
              GO TO 121
            END IF
          ELSE
            KM = 0
          END IF
        END IF
*
  120 CONTINUE
*
  121 CONTINUE
*
      VNEW = VNEW + FUNC(A0(L) + A) * B0(L)
      NEVAL = NEVAL + 1
*
      IF (NM .EQ. 0) THEN
        NM = NEND
        EPSM = 0.2 * SQRT(ABS(WM))
        WRITE (*,2001)
 2001   FORMAT (' (SUBR.DEHINT) slow decay on',
     $          ' negative side (warning only)')
      END IF
*
      IF (NP .EQ. 0) THEN
        NP = NEND
        EPSP = 0.2 * SQRT(ABS(WP))
        WRITE (*,2002)
 2002   FORMAT (' (SUBR.DEHINT) slow decay on',
     $          ' positive side (warning only)')
      END IF
*
      EPSQ = MAX(EPSQ, EPSM, EPSP)
*
*     ---- general step ----
*
      VOLD = H * VNEW
*
      DO 20 MSTEP = 1, NPOW
*
        VNEW = 0.0
*
        IH = IS
        IS = IS / 2
*
        DO 540 I = IS, NM, IH
          VNEW = VNEW
     $         + FUNC(AM(I,L) + A) * BM(I,L)
          NEVAL = NEVAL + 1
  540   CONTINUE
*
        DO 550 I = IS, NP, IH
          VNEW = VNEW
     $         + FUNC(AP(I,L) + A) * BP(I,L)
          NEVAL = NEVAL + 1
  550   CONTINUE
*
        VNEW = (VOLD + H * VNEW) * HALF
*
        IF (ABS(VNEW - VOLD) .LT. EPSQ) THEN
*
*        ---- converged and return ----
*
          V = VNEW
          RETURN
        END IF
*
        H = H * HALF
        VOLD = VNEW
*
   20 CONTINUE
*
*      WRITE (*,2003)
* 2003 FORMAT (' (SUBR.DEHINT) insufficient mesh',
*     $        ' refinement (warning only)')
      V = VNEW
      RETURN
*
      END
*
      SUBROUTINE HIAB
*******************************************************
*     Generate points and weights of the DE formula   *
*     over half infinite interval (A,infinity)        *
*                                                     *
*     Before the first call of DEHINT this routine    *
*     must be called                                  *
*                                                     *
*      L for DE half infinite transformation          *
*         L = 1  x = exp(0.5*t - exp(-t))             *
*         L = 2  x = exp(t - exp(-t))                 *
*         L = 3  x = exp(2 * sinh t)                  *
*******************************************************
      real*8 ONE,HALF
      real*8 H,EH,EN,ENI
      real*8 SH,CH
      real*8 AM,A0,AP,BM,B0,BP
*
      COMMON / COMDEF / AM(640,3),A0(3),AP(640,3),
     $                  BM(640,3),B0(3),BP(640,3)
      COMMON / COMDEN / NPOW,NEND
*
      PARAMETER (ONE = 1)
      PARAMETER (HALF = ONE / 2)
      PARAMETER (N6 = 6)
*
      NPOW = N6
      NEND = 5 * 2**(NPOW+1)
      H = ONE / 2**(NPOW+1)
      EH = EXP(H)
*
*     ---- DE transformation x = exp(0.5*t-exp(-t)) ----
*                        L = 1
*
      A0(1) = EXP(-ONE)
      B0(1) = 1.5D0 * A0(1)
      EN = 1.0
*
      DO 10 I = 1, NEND
        EN = EN * EH
        ENI = 1 / EN
        SH = HALF * H * I
        AP(I,1) = EXP(SH - ENI)
        BP(I,1) = (HALF + ENI) * AP(I,1)
        AM(I,1) = EXP(-SH - EN)
        BM(I,1) = (HALF + EN) * AM(I,1)
   10 CONTINUE
*
*     ---- DE transformation x = exp(t-exp(-t)) ----
*                        L = 2
*
      A0(2) = EXP(-ONE)
      B0(2) = 2.0 * A0(2)
      EN = 1.0
*
      DO 20 I = 1, NEND
        EN = EN * EH
        ENI = 1 / EN
        SH = H * I
        AP(I,2) = EXP(SH - ENI)
        BP(I,2) = (1 + ENI) * AP(I,2)
        AM(I,2) = EXP(-SH - EN)
        BM(I,2) = (1 + EN) * AM(I,2)
   20 CONTINUE
*
*     ---- DE transformation x = exp(2*sinh t) ----
*                        L = 3
*
      A0(3) = 1.0
      B0(3) = 2.0
      EN = 1.0
*
      DO 30 I = 1, NEND
        EN = EH * EN
        ENI = 1 / EN
        SH = EN - ENI
        CH = EN + ENI
        AP(I,3) = EXP(SH)
        BP(I,3) = CH * AP(I,3)
        AM(I,3) = 1 / AP(I,3)
        BM(I,3) = CH * AM(I,3)
   30 CONTINUE
*
      RETURN
      END
*
      FUNCTION  CELI1( KSQ )                                            CEK00100
C***********************************************************************CEK00200
C*                                                                     *CEK00300
C*  PURPOSE                                                            *CEK00400
C*      COMPUTES THE COMPLETE ELLIPITC INTEGRAL OF THE FIRST KIND      *CEK00500
C*      IN DOUBLE PRECISION.                                           *CEK00600
C*                                                                     *CEK00700
C*  FUNCTION SUBPROGRAM                                                *CEK00800
C*      CELI1( KSQ )                                                   *CEK00900
C*                                                                     *CEK01000
C*  DESCRIPTION OF PARAMETER                                           *CEK01100
C*      KSQ  - SQUARE OF THE MODULUS K.                                *CEK01200
C*                                                                     *CEK01300
C*  COPYRIGHT                                                          *CEK01400
C*      Y. ONODERA,    AUGUST, 1989.                                   *CEK01500
C*                                                                     *CEK01600
C***********************************************************************CEK01700
      IMPLICIT REAL*8 (A-Z)                                             CEK01800
      DATA  ONE / 1.0D0 /                                               CEK01900
C                                                                       CEK02000
      IF ( KSQ .LT. 0 ) GOTO 20                                         CEK02100
      X = ONE - KSQ                                                     CEK02200
C                                                                       CEK02300
      IF ( KSQ .LT. 0.9D0 ) THEN                                        CEK02400
C                             LANDEN TRANSFORMATION.                    CEK02500
        K1 = SQRT( X )                                                  CEK02600
        A = KSQ / ( K1 + ONE )**2                                       CEK02700
        Z = A * A                                                       CEK02800
        F = A + ONE                                                     CEK02900
        IF ( KSQ .GE. 0.27D0 ) THEN                                     CEK03000
C                             ONCE MORE LANDEN TRANSFORMATION.          CEK03100
          K1 = SQRT( ONE - Z )                                          CEK03200
          B = Z / ( K1 + ONE )**2                                       CEK03300
          Z = B * B                                                     CEK03400
          F = F * ( B + ONE )                                           CEK03500
        ENDIF                                                           CEK03600
C                             NOW  Z < 0.0062                           CEK03700
        CELI1 = F * ( ( ( ( ( (        0.0799 D0               *Z       CEK03800
     A  + 0.0951 308D0           )*Z + 0.1174 4540 4D0        )*Z       CEK03900
     B  + 0.1533 9807 879 D0     )*Z + 0.2208 9323 3455 5D0   )*Z       CEK04000
     C  + 0.3926 9908 1698 724D0 )*Z + 1.5707 9632 6794 89662D0 )       CEK04100
        RETURN                                                          CEK04200
C                                                                       CEK04300
      ELSE                                                              CEK04400
        IF ( X .EQ. 0 ) GOTO 10                                         CEK04500
        IF ( X .LT. 0 ) GOTO 20                                         CEK04600
        F = ((((( ((((( (((                                             CEK04700
     A      3.0D-2               * X  +  3.24D-2             ) * X      CEK04800
     B   +  3.495D-2           ) * X  +  3.7958D-2           ) * X      CEK04900
     C   +  4.1524 6D-2        ) * X  +  4.5829 54D-2        ) * X      CEK05000
     D   +  5.1127 765D-2      ) * X  +  5.7806 3778D-2      ) * X      CEK05100
     E   +  6.6482 4962 6D-2   ) * X  +  7.8202 0568 85D-2   ) * X      CEK05200
     F   +  9.4884 2366 536D-2 ) * X  +  1.2044 2708 3333D-1 ) * X      CEK05300
     G   +  0.1640 625D0       ) * X  +  0.25D0              ) * X      CEK05400
        G = ((((( ((((( (((                                             CEK05500
     A      1.1D-2               * X  +  1.20D-2             ) * X      CEK05600
     B   +  1.300D-2           ) * X  +  1.4144D-2           ) * X      CEK05700
     C   +  1.5522 7D-2        ) * X  +  1.7199 67D-2        ) * X      CEK05800
     D   +  1.9282 673D-2      ) * X  +  2.1939 3969D-2      ) * X      CEK05900
     E   +  2.5444 5076 0D-2   ) * X  +  3.0281 0668 95D-2   ) * X      CEK06000
     F   +  3.7384 0332 031D-2 ) * X  +  0.0488 2812 5D0     ) * X      CEK06100
     G   +  0.0703 125D0       ) * X  +  0.125D0             ) * X      CEK06200
     H   +  0.5D0                                                       CEK06300
        CELI1 = G * LOG( 16.0D0 / X ) - F                               CEK06400
        RETURN                                                          CEK06500
      ENDIF                                                             CEK06600
C                                                                       CEK06700
   10 CELI1 = 100.0D0                                                   CEK06800
      RETURN                                                            CEK06900
   20 WRITE(6,30) KSQ                                                  CEK07000
      CELI1 = 0                                                         CEK07100
      RETURN                                                            CEK07200
   30 FORMAT(' ARGUMENT OF CELI1 IS INVALID.  X=',1PD23.15)             CEK07300
C *****************************************************                 CEK07400
C              NUMERICAL TABLE OF CELI1                                 CEK07500
C       ARGUMENT                     CELI1                              CEK07600
C -----------------------------------------------------                 CEK07700
C   256/1024 ( 0.25       )   1.6857 5035 4812 5960                     CEK07800
C   276/1024 ( 0.269531...)   1.6964 8703 4351 1121                     CEK07900
C   277/1024 ( 0.270508...)   1.6970 3221 5452 5081                     CEK08000
C   512/1024 ( 0.5        )   1.8540 7467 7301 3719                     CEK08100
C   747/1024 ( 0.729492...)   2.1212 9467 8256 3108                     CEK08200
C   748/1024 ( 0.730469...)   2.1229 0614 6685 8874                     CEK08300
C   768/1024 ( 0.75       )   2.1565 1564 7499 6432                     CEK08400
C   921/1024 ( 0.899414...)   2.5753 4344 0273 8417                     CEK08500
C   922/1024 ( 0.900391...)   2.5799 3389 3410 7172                     CEK08600
C *****************************************************                 CEK08700
      END                                                               CEK08800
