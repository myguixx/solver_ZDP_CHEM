      SUBROUTINE SENKIN (LIN, LOUT, LINKCK, LSAVE, LIGN, LREST,
     1                   LR, R, LI, I, LC, C)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C     This subroutine computes the spontaneous ignition delay
C     for an adiabatic gas mixture in a closed system with either:
C     (1) pressure constant,
C     (2) volume a function of time,
C     (3) pressure constant and temperature a function of time.
C     The CHEMKIN package is used to handle chemical kinetics.
C     The time integration and kinetic sensitivity analysis
C     employ program DASAC.
C
C     Author:  Andy Lutz, Computational Mechanics Division, 8245,
C              Sandia National Laboratories, Livermore, CA, 94550
C              (415) 294-2761
C
C     VERSION 1.8
C     CHANGES FROM VERSION 1.1
C     1.  changed REAL*8 to DOUBLE PRECISION
C     2.  read linking file for array lengths
C     3.  remove need for COMMON/CKSTRT/
C     4.  allow upper/lower case input
C     CHANGES FROM VERSION 1.2
C     1.  read KERR and MAXTB from linking file
C     CHANGES FROM VERSION 1.3
C     1.  modified OPEN, 'tremain' statements for unicos
C     CHANGES FOR VERSION 1.5
C     1.  CALL CKLEN will check linking file version and precision
C     2.  Read statements for new linking file
C     CHANGES FOR VERSION 1.6
C     1.  Call CKRHOY to compute variable density in RCONT
C     CHANGES FOR VERSION 1.7
C     1.  Call CKLEN and CKINDX instead of reading the linking file,
C         and changed position of some pointers in I, R, and C
C         arrays
C     CHANGES FOR VERSION 1.8
C     1.  Place reaction constant (D) at front of RPAR (remove A,B,E).
C     2.  Call CKRDEX instead of CKA to perturb reaction constant.
C
      DIMENSION R(*), I(*)
      CHARACTER C(*)*(*)
C
      LOGICAL LSENS
C
      COMMON /POINT/ IPICK, IPRCK, IPWT, IPWDOT, IPU, IPRD
C
      WRITE (LOUT, 15)
      WRITE (LIGN, 15)
   15 FORMAT(
     1/' SENKIN:  Sensitity Analysis',
     2/'          Author: Andy Lutz',
     3/'          CHEMKIN-II Version 1.8, December 1991',
C*****precision > double
     4/'          DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     4/'          SINGLE PRECISION')
C*****END precision > single
C
C     DETERMINE MECHANISM SIZE
C
      IPICK = 20
C
      CALL CKLEN (LINKCK, LOUT, LENI, LENR, LENC)
      IF (LENI .LE. LI .AND. LENR.LE.LR .AND. LENC.LE.LC) THEN
         CALL CKINIT (LI, LR, LC, LINKCK, LOUT, I(IPICK), R, C)
         CALL CKINDX (I, R, MM, KK, II, NFIT)
      ENDIF
C
C        READ SENSITIVITY AND PROBLEM CHOICE KEYWORDS
C
      CALL REDSEN (ICASE, LIGN, LIN, LOUT, LSENS)
C
C         COMPUTE DASAC WORK SPACE
C
      IF (ICASE.LT.1 .OR. ICASE.GT.5) THEN
         WRITE (LOUT, '(/1X,A,/)') ' Stop, ICASE not found in SENKIN.'
         STOP
      ELSEIF (ICASE .GT. 3) THEN
         NSYS = KK
      ELSE
         NSYS  = KK + 1
      ENDIF
C
      IF (LSENS) THEN
         NEQ    = NSYS * (II + 1)
         LSDAS  = NSYS
      ELSE
         NEQ = NSYS
         LSDAS  = 1
      ENDIF
      LIDAS  = 20 + NEQ
      LRDAS  = 40 + NEQ*9 + NSYS*NSYS
C
C         SET RPAR POINTERS
C     NOTE: MUST HAVE IPRD 1ST IN RPAR!
C
      IPRD   = 1
      IPRCK  = IPRD  + II
      IPWT   = IPRCK + LENR
      IPWDOT = IPWT  + KK
      IPU    = IPWDOT+ KK
      IPTOT  = IPU   + KK
C
C         APPORTION REAL WORK SPACE
C
      NRPAR  = 1
      NRDAS  = NRPAR + IPTOT
      NSDAS  = NRDAS + LRDAS
      NATOL  = NSDAS + LSDAS
      NRTOL  = NATOL + NEQ
      NXMOL  = NRTOL + NEQ
      NZ     = NXMOL + KK
      NZP    = NZ    + NEQ
      LRTOT  = NZP   + NEQ
C
C        APPORTION INTEGER WORK SPACE
C
      NIPAR  = 1
      NIDAS  = NIPAR + LENI + IPICK
      LITOT  = NIDAS + LIDAS
C
C        APPORTION CHARACTER WORK SPACE
C
      IPCCK = 1
      NKSYM = IPCCK + LENC
      LCTOT = NKSYM + KK
C
C          CHECK FOR SUFFICIENT SPACE
C
      WRITE (LOUT, 7020) LI, LITOT, LR, LRTOT, LC, LCTOT
 7020 FORMAT (/, '                Working Space Requirements',
     1        /, '                 Provided        Required ',
     2        /, ' Integer  ', 2I15,
     3        /, ' Real     ', 2I15,
     4        /, ' Character', 2I15, /)
C
      IF (LRTOT.GT.LR .OR. LITOT.GT.LI .OR. LCTOT.GT.LC) THEN
         WRITE (LOUT, *) '  Stop, not enough work space provided.'
         STOP
      ENDIF
C
C          GO TO MAIN LEVEL
C
      CALL BEGIN (NSYS, NEQ, ICASE, II, KK, LENI, LENR, LENC,
     1            LINKCK, LIN, LOUT, LSAVE, LIGN, LREST, LSENS,
     2            LIDAS, LRDAS, LSDAS, I(NIDAS), R(NRDAS), R(NSDAS),
     3            R(NRPAR), I(NIPAR), R(NZ), R(NZP), R(NRTOL),
     4            R(NATOL), R(NXMOL), C(NKSYM), C(IPCCK))
C
      STOP
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE BEGIN (NSYS, NEQ, ICASE, II, KK, LENICK, LENRCK,
     1                  LENCCK, LINKCK, LIN, LOUT, LSAVE, LIGN, LREST,
     2                  LSENS, LIDAS, LRDAS, LSDAS, IDWORK, DWORK,
     3                  SDWORK, RPAR, IPAR, Z, ZP, RTOL, ATOL, XMOL,
     5                  KSYM, CCKWRK)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), XMOL(*), DWORK(*), IDWORK(*), SDWORK(*),
     1          RTOL(*), ATOL(*), RPAR(*), IPAR(*), TOLS(4)
      CHARACTER*(*) CCKWRK(*), KSYM(*)
C
      EXTERNAL RCONP, RCONV, RCONT, RVOLT, RTEMP
C
      LOGICAL LSENS, RESTRT, IERR
C
C        PHYSICAL VARIABLES IN COMMON
C
      COMMON /RES1/ P
      COMMON /RES2/ TOTMAS
      COMMON /RES3/ T
      COMMON /RES4/ RHO
C
C        POINTERS IN COMMON
C
      COMMON /POINT/ IPICK, IPRCK, IPWT, IPWDOT, IPU, IPRD
      IPAR(1) = KK
      IPAR(2) = IPRCK
      IPAR(3) = IPRD
      IPAR(6) = IPWT
      IPAR(7) = IPWDOT
      IPAR(8) = IPU
      IPAR(9) = IPICK
      IPAR(11) = II
C
C         INITIALIZE CHEMKIN
C         (RESET CHEMKIN WORKSPACE)
C
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     1             IPAR(IPICK), RPAR(IPRCK), CCKWRK )
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      CALL CKWT   (IPAR(IPICK), RPAR(IPRCK), RPAR(IPWT))
      CALL CKRP   (IPAR(IPICK), RPAR(IPRCK), RU, RUC, PATM)
C
C      LOAD THE REACTION PARAMETERS INTO RPAR
C
      DO 15 I = 1, II
        CALL CKRDEX (I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C        READ KEYWORD INPUT
C
      CALL REDKEY (DELT, KK, KSYM, LIN, LOUT, PA,
     1             RESTRT, T, TLIM, TOLS, TRES, TSTOP, XMOL)
      P = PA * PATM
C
C       PHYSICAL INITIAL CONDITIONS
C
C     CASES:
C
C     ICASE = 1 == CONSTANT PRESSURE
C     ICASE = 2 == CONSTANT VOLUME
C     ICASE = 3 == VOLUME GIVEN IN TIME
C     ICASE = 4 == CONSTANT TEMPERATURE, PRESSURE
C     ICASE = 5 == TEMPERATURE GIVEN IN TIME
C
      IF (RESTRT) THEN
C
C          READ RESTART DATA
C
         READ (LREST) IDUM
         READ (LREST) IDUM, KKR
C
         IF (KKR .NE. KK)
     1      WRITE (LOUT, '(/3X, A, /9X, A, /9X, A, I6, /9X, A, I6)')
     2      'Stop! Number of species in restart file does not match',
     3      'the number from the CHEMKIN link file.',
     4      'No. species in restart file = ', KKR,
     5      'No. species in linking file = ', KK
C
         READ (LREST) TFIL, P, T, (XMOL(K), K = 1, KK)
         IF (TRES .LT. 0.) THEN
            TIM = TFIL
         ELSE
            TIM = TRES
         ENDIF
C
         IF (ICASE .LE. 3) THEN
            Z(1) = T
            DO 50 K = 1, KK
               Z(K+1) = XMOL(K)
   50       CONTINUE
            CALL CKYTX  (Z(2), IPAR(IPICK), RPAR(IPRCK), XMOL)
            CALL CKRHOY (P, T, Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
            IF (ICASE .EQ. 3) THEN
               CALL VOLT (TIM, VZERO, DVDT)
               TOTMAS = RHO * VZERO
            ENDIF
C
         ELSEIF (ICASE .EQ. 4) THEN
            DO 60 K = 1, KK
               Z(K) = XMOL(K)
   60       CONTINUE
            CALL CKYTX  (Z, IPAR(IPICK), RPAR(IPRCK), XMOL)
            CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
C
         ELSEIF (ICASE .EQ. 5) THEN
            DO 70 K = 1, KK
               Z(K) = XMOL(K)
   70       CONTINUE
            CALL TEMPT  (TIM, T)
            CALL CKYTX  (Z, IPAR(IPICK), RPAR(IPRCK), XMOL)
            CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
C
         ELSE
            WRITE (LOUT, '(/1X,A,/)') ' Stop, ICASE not found in BEGIN.'
            STOP
         ENDIF
C
      ELSE
C
C      ORIGINAL JOB
C
         TIM = 0.E+0
         IF (ICASE .LE. 3) THEN
            Z(1) = T
            CALL CKXTY  (XMOL, IPAR(IPICK), RPAR(IPRCK), Z(2))
            CALL CKRHOY (P, T, Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
            IF (ICASE .EQ. 3) THEN
               CALL VOLT (TIM, VZERO, DVDT)
               TOTMAS = RHO * VZERO
            ENDIF
C
         ELSEIF (ICASE .EQ. 4) THEN
            CALL CKXTY  (XMOL, IPAR(IPICK), RPAR(IPRCK), Z)
            CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
C
         ELSEIF (ICASE .EQ. 5) THEN
            CALL TEMPT  (TIM, T)
            CALL CKXTY  (XMOL, IPAR(IPICK), RPAR(IPRCK), Z)
            CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
C
         ELSE
            WRITE (LOUT, '(/1X,A,/)') ' Stop, ICASE not found in BEGIN.'
            STOP
         ENDIF
      ENDIF
C
      IF (RHO .LE. 0.0) THEN
         WRITE (LOUT, '(/1X,A,/)') ' Stop, initial density < 0.'
         STOP
      ENDIF
C
C       PRINT INITIAL CONDITIONS
C
      IF (RESTRT) THEN
         WRITE (LOUT, 7000)
         WRITE (LIGN, 7000)
      ENDIF
C
      IF (ICASE .EQ. 1) THEN
         WRITE (LOUT, 7111)
         WRITE (LIGN, 7111)
      ELSEIF (ICASE .EQ. 2) THEN
         WRITE (LOUT, 7112)
         WRITE (LIGN, 7112)
      ELSEIF (ICASE .EQ. 3) THEN
         WRITE (LOUT, 7113)
         WRITE (LIGN, 7113)
      ELSEIF (ICASE .EQ. 4) THEN
         WRITE (LOUT, 7114)
         WRITE (LIGN, 7114)
      ELSEIF (ICASE .EQ. 5) THEN
         WRITE (LOUT, 7115)
         WRITE (LIGN, 7115)
      ENDIF
      WRITE (LOUT, 7103)
      WRITE (LIGN, 7103)
      WRITE (LOUT, 7100) PA, T, RHO
      WRITE (LIGN, 7100) PA, T, RHO
      IF (ICASE .EQ. 3) THEN
         WRITE (LOUT, 7104) TOTMAS
         WRITE (LIGN, 7104) TOTMAS
      ENDIF
      WRITE (LIGN, 7101)
      WRITE (LOUT, 7101)
      DO 130 K = 1, KK
         WRITE (LIGN, 7102) KSYM(K), XMOL(K)
         WRITE (LOUT, 7102) KSYM(K), XMOL(K)
130   CONTINUE
C
C        INTEGRATION ROUTINE TO RUN PROBLEM
C
      IF (ICASE .EQ. 1) THEN
         CALL SENS13 (RCONP, ICASE, NSYS, KK, II, DELT, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LIGN, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
      ELSEIF (ICASE .EQ. 2) THEN
         CALL SENS13 (RCONV, ICASE, NSYS, KK, II, DELT, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LIGN, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
      ELSEIF (ICASE .EQ. 3) THEN
         CALL SENS13 (RVOLT, ICASE, NSYS, KK, II, DELT, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LIGN, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
      ELSEIF (ICASE .EQ. 4) THEN
         CALL SENS45 (RCONT, ICASE, NSYS, II, DELT, TSTOP, TIM, PATM,
     1                TLIM, LOUT, LSAVE, LIGN, LSENS, LIDAS, LRDAS,
     2                LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR, IPAR,
     3                ATOL, RTOL, TOLS, XMOL, KSYM)
      ELSEIF (ICASE .EQ. 5) THEN
         CALL SENS45 (RTEMP, ICASE, NSYS, II, DELT, TSTOP, TIM, PATM,
     1                TLIM, LOUT, LSAVE, LIGN, LSENS, LIDAS, LRDAS,
     2                LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
      ENDIF
C
7000  FORMAT (/5X,'Restart calculation from previous solution.'/)
7100  FORMAT ('  Pressure (atm)  =', 1PE12.4,/,
     1        '  Temperature (K) =', 1PE12.4,/,
     2        '  Density (gm/cc) =', 1PE12.4)
7101  FORMAT (/,'  Mole Fractions:')
7102  FORMAT (1X,A10,2H =,1PE11.4)
7103  FORMAT (/,'  Initial Conditions:',/)
7104  FORMAT ('  Mass (gm)       =', 1PE12.4)
7111  FORMAT (/5X,'Pressure is held constant.'/)
7112  FORMAT (/5X,'Volume is held constant.'/)
7113  FORMAT (/5X,'Volume is a function of time.'/)
7114  FORMAT (/5X,'Temperature is held constant.'/)
7115  FORMAT (/5X,'Temperature is a function of time.'/)
      STOP
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SENS13 (RES, ICASE, NSYS, KK, II, DTOUT, TSTOP, TIM,
     1                   PATM, TLIM, LOUT, LSAVE, LIGN, LSENS, LIW,
     2                   LRW, LSW, Z, ZP, ELWRK, IELWRK, SWORK, RPAR,
     3                   IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
C
C  This module directs the integration for cases 1-3, where temperature
C  is not known and the energy equation is included.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(NSYS,*), ZP(NSYS,*), ELWRK(*), IELWRK(*), SWORK(*),
     1          ISEN(5), RPAR(*), IPAR(*), INFO(15), RTOL(NSYS,*),
     2          ATOL(NSYS,*), TOLS(*), XMOL(*)
C
      EXTERNAL RES, DRES
C
      LOGICAL LSENS
      CHARACTER KSYM(*)*(*)
C
      COMMON /RES1/ P
C
      DATA ISEN /5*0/, INFO /15*0/, DTLIM / 400. /, NOSAV /1/
C
C       SET OUTPUT FILE FOR USER COMFORMATION
C
      LDTS = 11
      OPEN (LDTS, FORM='FORMATTED', FILE = 'output/skout_datasheet')
C
C       SET PARAMETERS FOR DASAC
C
      INFO(3) = 1
      INFO(2) = 1
C
C        TOLERANCES ON DEPENDENT VARIABLES
C
      DO 10 I = 1, NSYS
         RTOL(I,1) = TOLS(1)
         ATOL(I,1) = TOLS(2)
10    CONTINUE
C
C        TOLERANCES ON SENSITIVITY COEF.
C
      IF (LSENS) THEN
         ISEN(1) = 1
         ISEN(5) = II
         DO 20 I = 1, NSYS
            DO 20 J = 2, II+1
               RTOL(I,J) = TOLS(3)
               ATOL(I,J) = TOLS(4)
20       CONTINUE
      ENDIF
C
C       COMPUTE INITIAL DERIVATIVES
C
C*****precision > single
CC
C      CALL SDINIT (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
C     1            IRES, ISEN(1), ISEN(4), RES, DRES)
CC
C*****END precision > single
C*****precision > double
C
      CALL DDINIT (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
     1             IRES, ISEN(1), ISEN(4), RES, DRES)
C
C*****END precision > double
C
C       PRINT    Initial condition.
C
      WRITE (LSAVE) LSENS
      WRITE (LSAVE) NSYS, KK, II
      WRITE (LSAVE) TIM, P, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) (( Z(I,J), I=1,NSYS), J = 2, II+1)
      WRITE (LOUT, '(//A/)') '  Time Integration:'
      WRITE (LIGN, '(//A/)') '  Time Integration:'
      CALL TEXT13 (IPAR, KK, KSYM, LIGN, LOUT,
     1             P, PATM, RPAR, TIM, XMOL, Z)
      WRITE (LDTS, 7720)
      WRITE (LDTS, *) (KSYM(I), I = 1, KK)
      CALL TEXT13DTS (IPAR, KK, KSYM, LDTS,
     1             P, PATM, RPAR, TIM, XMOL, Z)
 7720 FORMAT(/,' t(sec)     P(atm)     T(K)    ', $)
C
C       A variable used for sensitivity analysis. Z(ls,J) --> ls = 1: temperature, ls = 2...NSYS: species
C
      ls = 1
C
C       PRINT    Sensitivity coefficient.
C
      IF (LSENS) THEN
         WRITE (LIGN, '(A)') ' Sensitivity coefficient'
         WRITE (LIGN, '(I4, 1PE12.4)') (J-1, Z(ls,J), J = 2, II+1)
      ENDIF
C
C       INTEGRATION LOOP
C
      TPRINT = DTOUT + TIM
      IFLG = 0
250   CONTINUE
C
C         CALL THE O.D.E. SOLVER
C
      TIMOLD = TIM
      TOLD   = Z(1,1)
C
C*****precision > single
C      CALL SDASAC (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
C     1             RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
C     2             IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C*****END precision > single
C*****precision > double
      CALL DDASAC (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
     1             RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
     2             IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C*****END precision > double
C
      IF (IDID .LT. 0) THEN
         WRITE (LOUT, '(1X,A,I3)') 'IDID =', IDID
         STOP
      ENDIF
C
C           PRINT TO PLOT FILE
C
      WRITE (LSAVE) TIM, P, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) ((Z(I,J),I = 1,NSYS), J = 2,II+1)
      NOSAV = NOSAV + 1
C
C         CHECK FOR THERMAL RUNAWAY
C
      T = Z(1,1)
      IF (T.GE.TLIM .AND. IFLG.EQ.0) THEN
         DT = ELWRK(7)
         TIGN = TIMOLD + DT *(TLIM-TOLD)/(T-TOLD)
         IFLG = 1
      ENDIF
C
C           PRINT OUT SOLUTION
C
      IF (TIM .GE. TPRINT) THEN
         CALL TEXT13 (IPAR, KK, KSYM, LIGN, LOUT,
     1                P, PATM, RPAR, TIM, XMOL, Z)
         CALL TEXT13DTS (IPAR, KK, KSYM, LDTS,
     1                P, PATM, RPAR, TIM, XMOL, Z)
         TLASTP = TIM
         TPRINT = TPRINT + DTOUT
C
C       PRINT    Sensitivity coefficient.
C
         IF (LSENS) THEN
            WRITE (LIGN, '(A)') ' Sensitivity coefficient'
            WRITE (LIGN, '(I4, 1PE12.4)') (J-1, Z(ls,J), J = 2, II+1)
         ENDIF
      ENDIF
C
C*****unicos timelimit
C      CALL TREMAIN (RITIM)
C      IF (RITIM .LE. 30.0) THEN
C         WRITE (LOUT, '(/5X,A,/)')
C     1   'Stop!  Senkin job is nearing its timelimit.'
C         GO TO 1000
C      ENDIF
C*****END unicos timelimit
C
      IF (TIM .LT. TSTOP) GO TO 250
C
C        LAST TEXT PRINT
C
1000  CONTINUE
      IF (TLASTP .NE. TIM) THEN
         CALL TEXT13 (IPAR, KK, KSYM, LIGN,
     1                   LOUT, P, PATM, RPAR, TIM, XMOL, Z)
C
C       PRINT    Sensitivity coefficient.
C
         IF (LSENS) THEN
            WRITE (LIGN, '(A)') ' Sensitivity coefficient'
            WRITE (LIGN, '(I4, 1PE12.4)') (J-1, Z(ls,J), J = 2, II+1)
         ENDIF
      ENDIF
C
      WRITE (LIGN, 7040) TIGN
      WRITE (LIGN, 7045) TLIM
      WRITE (LOUT, 7040) TIGN
      WRITE (LOUT, 7045) TLIM
      WRITE (LIGN, 7050) NOSAV
      WRITE (LOUT, 7050) NOSAV
C
C         FORMATS
C
 7040 FORMAT(//1X,'  Integration completed:',/2X,
     1' Ignition Time (sec) = ',1PE12.4)
 7045 FORMAT(2X,' Temp criteria (K) = ',1PE12.4)
 7050 FORMAT(/1X,'  Binary file has ',I6,' time datasets.')
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SENS45 (RES, ICASE, NSYS, II, DTOUT, TSTOP, TIM,
     1                   PATM, TLIM, LOUT, LSAVE, LIGN, LSENS, LIW,
     2                   LRW, LSW, Z, ZP, ELWRK, IELWRK, SWORK, RPAR,
     3                   IPAR, ATOL, RTOL, TOLS, XMOL, KSYM)
C
C  This module directs the integration for cases 4-5, where temperature
C  is specified and so there is no energy equation. Here  NSYS = KK.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(NSYS,*), ZP(NSYS,*), ELWRK(*), IELWRK(*),
     1          SWORK(*), ISEN(5), RPAR(*), IPAR(*), INFO(15),
     2          RTOL(NSYS,*), ATOL(NSYS,*), TOLS(*), XMOL(*)
C
      EXTERNAL RES, DRES
C
      LOGICAL LSENS
      CHARACTER KSYM(*)*(*)
C
      COMMON /RES1/ P
      COMMON /RES3/ T
C
      DATA ISEN /5*0/, INFO /15*0/, NOSAV /1/
C
      KK = NSYS
C
C       SET PARAMETERS FOR DASAC
C
      INFO(3) = 1
      INFO(2) = 1
C
C        TOLERANCES ON DEPENDENT VARIABLES
C
      DO 10 I = 1, NSYS
         RTOL(I,1) = TOLS(1)
         ATOL(I,1) = TOLS(2)
10    CONTINUE
C
C        TOLERANCES ON SENSITIVITY COEF.
C
      IF (LSENS) THEN
         ISEN(1) = 1
         ISEN(5) = II
         DO 20 I = 1, NSYS
            DO 20 J = 2, II+1
               RTOL(I,J) = TOLS(3)
               ATOL(I,J) = TOLS(4)
20       CONTINUE
      ENDIF
C
C       COMPUTE INITIAL DERIVATIVES
C
C*****precision > single
C      CALL SDINIT (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
C     1             IRES, ISEN(1), ISEN(4), RES, DRES)
C*****END precision > single
C*****precision > double
      CALL DDINIT (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
     1             IRES, ISEN(1), ISEN(4), RES, DRES)
C*****END precision > double
C
C       PRINT
C
      WRITE (LSAVE) LSENS
      WRITE (LSAVE) NSYS, KK, II
      WRITE (LSAVE) TIM, P, T, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) (( Z(I,J), I=1,NSYS), J = 2, II+1)
      WRITE (LOUT, '(//A/)') '  Time Integration:'
      WRITE (LIGN, '(//A/)') '  Time Integration:'
      CALL TEXT45 (IPAR, KK, KSYM, LIGN, LOUT,
     1             P, PATM, RPAR, T, TIM, XMOL, Z)
C
C       INTEGRATION LOOP
C
      TPRINT = DTOUT + TIM
      IFLG = 0
250   CONTINUE
C
C         CALL THE O.D.E. SOLVER
C
C*****precision > single
C      CALL SDASAC (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
C     1             RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
C     2             IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C*****END precision > single
C*****precision > double
      CALL DDASAC (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
     1             RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
     2             IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C*****END precision > double
C
      IF (IDID .LT. 0) THEN
         WRITE (LOUT, '(1X,A,I3)') 'IDID =', IDID
         STOP
      ENDIF
C
C           PRINT TO PLOT FILE
C
      WRITE (LSAVE) TIM, P, T, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) ((Z(I,J), I = 1,NSYS), J = 2,II+1)
      NOSAV = NOSAV + 1
C
C           PRINT OUT SOLUTION
C
      IF (TIM .GE. TPRINT) THEN
         CALL TEXT45 (IPAR, KK, KSYM, LIGN, LOUT,
     1                P, PATM, RPAR, T, TIM, XMOL, Z)
         TLASTP = TIM
         TPRINT = TPRINT + DTOUT
      ENDIF
C
C*****unicos timelimit
C      CALL TREMAIN (RITIM)
C      IF (RITIM .LE. 30.) THEN
C         WRITE (LOUT, '(/5X,A,/)')
C     1   'Stop!  Senkin job is nearing its timelimit.'
C         GO TO 1000
C      ENDIF
C*****END unicos timelimit
C
      IF (TIM .LT. TSTOP) GO TO 250
C
C        LAST TEXT PRINT
C
1000  CONTINUE
      IF (TLASTP .NE. TIM) CALL TEXT45 (IPAR, KK, KSYM, LIGN,
     1                LOUT, P, PATM, RPAR, T, TIM, XMOL, Z)
C
      WRITE (LIGN, 7050) NOSAV
      WRITE (LOUT, 7050) NOSAV
C
C         FORMATS
C
 7040 FORMAT(//1X,'  Integration completed:',/2x,
     1' Ignition Time (sec) = ',1PE12.4)
 7050 FORMAT(/1X,'  Binary file has ',I6,' time datasets.')
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RTEMP (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
C
      COMMON /RES1/ P
      COMMON /RES3/ T
C
C  Residual of differential equations for case where pressure
C  is constant, and temperature is a user specified function of
C  time.
C
C  Variables:
C    Z(K) = species mass fractions
C    RPAR = array of reaction pre-exponential constants
C    T    = temperature (Kelvin)
C    P    = pressure (dyne/cm2) - constant in time
C    TIME = time (sec)
C
C  User supplys a subroutine for the temperature:
C    SUBROUTINE TEMPT (TIME, T)
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPH    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
C
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C       USER SUPPLIES TEMPERATURE
C
      CALL TEMPT (TIME, T)
C
C         CALL CHEMKIN SUBROUTINES
C
      CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
      CALL CKWYP  (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RPAR(IPWDOT))
      IF (RHO .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero density in RTEMP.'
         STOP
      ENDIF
      VOLSP = 1. / RHO
C
C         SPECIES EQUATIONS
C
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K) = ZP(K) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONT (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
C
      COMMON /RES1/ P
      COMMON /RES3/ T
C
C  Residual of differential equations for case where temperature,
C  and pressure are constant.  Density varies with the mean molecular
C  weight of the mixture.
C
C  Variables:
C    Z(K) = species mass fractions
C    T    = temperature (Kelvin) - constant in time
C    P    = pressure  (dyne/cm2) - constant in time
C    TIME = time (sec)
C    RHO  = density of mixture (gm/cm3)
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPH    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
C
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C         CALL CHEMKIN SUBROUTINES
C
      CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
      VOLSP = 1./ RHO
      CALL CKWYP (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RPAR(IPWDOT))
C
C         SPECIES EQUATIONS
C
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K) = ZP(K) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONP (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
C
      COMMON /RES1/ P
C
C  Residual of differential equations for constant pressure case
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2) - constant in time
C    RHO    = density (gm/cm3)
C    RPAR   = array of reaction pre-exponential constants
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPH    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
C
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C         CALL CHEMKIN SUBROUTINES
C
      CALL CKRHOY (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
      CALL CKCPBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CPB)
      CALL CKWYP  (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKHMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPH))
      IF (RHO .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero density in RCONP.'
         STOP
      ENDIF
      VOLSP = 1. / RHO
C
C         ENERGY EQUATION
C
      SUM = 0.
      DO 100 K = 1, KK
         K1 = K-1
         SUM = SUM + RPAR(IPH+K1) * RPAR(IPWDOT+K1) * RPAR(IPWT+K1)
 100  CONTINUE
      DELTA(1) = ZP(1) + VOLSP *SUM /CPB
C
C         SPECIES EQUATIONS
C
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RVOLT (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
C
      EXTERNAL VOLT
C
      COMMON /RES1/ P
      COMMON /RES2/ TOTMAS
C
C  Residual of differential equations for case where volume is
C  a user-specified function of time.
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2)
C    RHO    = density (gm/cm3)
C    TOTMAS = mass (gm) - constant for closed system
C    TIME   = time (sec)
C
C  User supplys a subroutine for the volume:
C    SUBROUTINE VOLT (TIME, VOL, DVDT)
C     VOL    = volume of system
C     DVDT   = time derivative of system volume
C     VOLSP  = specific volume
C     VDOT   = time derivative of specific volume
C  Note: Consistent units for volume are (cm3), but the volume can
C        be considered to be normalized and therefore unitless.
C        The problem is solved using intensive variables, so the
C        solution is independent of extensive variables such as
C        the volume and total mass.  Subroutine VOLT is called at
C        time zero and a total mass is computed, but the solution
C        only depends on the density, so the total mass and total
C        volume are not important.
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPU    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
C
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C        GET VOLUME AS FCN OF TIME
C
      CALL VOLT (TIME, VOL, DVDT)
      VOLSP = VOL  / TOTMAS
      VDOT  = DVDT / TOTMAS
      IF (VOLSP .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero volume in RVOLT.'
         STOP
      ENDIF
      RHO = 1./ VOLSP
C
C         CALL CHEMKIN SUBROUTINES
C
      CALL CKCVBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CVB)
      CALL CKWYR  (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKUMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPU))
      CALL CKPY   (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), P)
C
C         ENERGY EQUATION
C
      SUM = 0.
      DO 100 K = 1, KK
         K1 = K-1
         SUM = SUM + RPAR(IPU+K1) * RPAR(IPWDOT+K1) * RPAR(IPWT+K1)
 100  CONTINUE
      DELTA(1) = ZP(1) + (P *VDOT + VOLSP *SUM) /CVB
C
C         SPECIES EQUATIONS
C
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONV (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
C
      EXTERNAL VOLT
C
      COMMON /RES1/ P
      COMMON /RES4/ RHO
C
C  Residual of differential equations when volume is constant
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2)
C    RHO    = density (gm/cm3) - constant in time
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPU    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
C
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C        CONSTANT DENSITY
C
      VOLSP = 1. / RHO
C
C         CALL CHEMKIN SUBROUTINES
C
      CALL CKCVBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CVB)
      CALL CKWYR  (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKUMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPU))
      CALL CKPY   (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), P)
C
C         ENERGY EQUATION
C
      SUM = 0.
      DO 100 K = 1, KK
         K1 = K-1
         SUM = SUM + RPAR(IPU+K1) * RPAR(IPWDOT+K1) * RPAR(IPWT+K1)
 100  CONTINUE
      DELTA(1) = ZP(1) + VOLSP *SUM /CVB
C
C         SPECIES EQUATIONS
C
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE CKA (II, WORK, RA)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C        THIS SUBROUTINE PLACES THE PRE-EXPONENTIAL CONSTANT
C        FOR THE ARRHENIUS COEFFICIENT OF THE REACTIONS INTO THE
C        CHEMKIN WORK ARRAY.    A. LUTZ (2/87)
C
C      INPUT
C        WORK   - ARRAY OF REAL INTERNAL WORK SPACE.  THE WORK ARRAY IS
C                 INITIALIZED BY THE CALL TO SUBROUTINE CKINIT.
C                    DIMENSION WORK(*) AT LEAST LENWK. SEE CKINIT FOR
C                    DETAILS ON THE REQUIRED LENGTH OF WORK.
C
C      OUTPUT
C        RA     - ARRAY OF PRE-EXPONENTIAL CONSTANTS FOR THE II
C                 REACTIONS.
C                    CGS UNITS - DEPENDS ON REACTION.
C                    DIMENSION RA(*) AT LEAST II.
C**           END ARGUMENT DESCRIPTION FOR SUBROUTINE CKA
      DIMENSION RA(*), WORK(*)
C
      DO 50 I = 1, II
         CALL CKRAEX (-I, WORK, RA(I))
   50 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE REDKEY (DELT, KK, KSYM, LIN, LOUT, P,
     1                   RESTRT, T, TLIM, TOLS, TRES, TSTOP, REAC)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION REAC(KK), TOLS(*), VALUE(5)
C
      CHARACTER KEYWRD*4, KSYM(*)*(*), LINE*80
      LOGICAL KERR, IERR, RESTRT
      DATA IREAC /0/, KERR /.FALSE./
      ICASE = 0
C
C           ISSUE A PROMPT
C
      WRITE (LOUT, *) ' Enter keywords:'
      WRITE (LOUT, *) '   PRES, TEMP, REAC, TIME, DELT, TLIM, REST,'
      WRITE (LOUT, *) '   RTOL, ATOL, RTLS, ATLS, or END'
C
C        DEFAULTS
C
      TRES = -1.
      TLIM = 0.0
      RESTRT = .FALSE.
      DO 10 K = 1, KK
         REAC(K) = 0.
10    CONTINUE
      TOLS(1) = 1.E-8
      TOLS(2) = 1.E-20
      TOLS(3) = 1.E-5
      TOLS(4) = 1.E-5
C
C         READ NEXT INPUT LINE
C
      WRITE (LOUT, '(/A,/)') '     Keyword input:'
C
90    CONTINUE
      LINE = ' '
      READ  (LIN,  '(A4,A80)',END=100) KEYWRD, LINE
      WRITE (LOUT, '(A4,A)') KEYWRD, LINE
      CALL UPCASE (KEYWRD)
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/' .OR.
     1    KEYWRD(1:1) .EQ. '!') GO TO 90
C
C
C         TOLERANCES (OPTIONAL)
C
      IF (KEYWRD .EQ. 'RTOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(1), IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
      ELSEIF (KEYWRD .EQ. 'ATOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(2), IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
      ELSEIF (KEYWRD .EQ. 'RTLS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(3), IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
      ELSEIF (KEYWRD .EQ. 'ATLS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(4), IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         PRESSURE
C
      ELSEIF (KEYWRD .EQ. 'PRES') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         TEMPERATURE
C
      ELSEIF (KEYWRD .EQ. 'TEMP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         IGNITION TEMPERATURE
C
      ELSEIF (KEYWRD .EQ. 'TLIM') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TLIM, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         SPECIES
C
      ELSEIF (KEYWRD .EQ. 'REAC') THEN
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            CALL ERRPR (KERR, KEYWRD, LOUT)
         ELSE
            REAC(KSPEC) = VALUE(1)
         ENDIF
C
C         INTEGRATION TIME
C
      ELSEIF (KEYWRD .EQ. 'TIME') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TSTOP, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         TIME OUTPUT INTERVAL
C
      ELSEIF (KEYWRD .EQ. 'DELT') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DELT, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C         RESTART JOB
C
      ELSEIF (KEYWRD .EQ. 'REST') THEN
         RESTRT = .TRUE.
C
C         RESTART TIME
C
      ELSEIF (KEYWRD .EQ. 'TRES') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TRES, IERR)
         IF (IERR) CALL ERRPR (KERR, KEYWRD, LOUT)
C
C        LAST CARD
C
      ELSEIF (KEYWRD .EQ. 'END ') THEN
         GO TO 100
C
C        KEYWORD IS BOGUS
C
      ELSE
        WRITE (LOUT, '(/2X,A,A,A)')
     1  ' Warning, keyword ', KEYWRD, ' was not found.'
      ENDIF
      GO TO 90
C
  100 CONTINUE
C
C         DONE READING CARDS
C
C     SET LIMIT TEMPERATURE IF NOT INPUT
      IF (TLIM.EQ.0.) TLIM = T + 400.
C
      IF (RESTRT) RETURN
C
C         CHECK THE REACTANT SUM
C
      SUMR = 0.
      DO 6100 K = 1, KK
         SUMR = SUMR + REAC(K)
6100  CONTINUE
C
C          CHECK FOR NECESSARY INPUT
C
      IF (P .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "PRES" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (T .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "TEMP" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (TSTOP .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "TIME" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (DELT .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "DELT" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (SUMR .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "REAC" not properly specified.'
         KERR = .TRUE.
      ENDIF
C
C        STOP IF ERRORS ENCOUNTERED
C
      IF (KERR) STOP
C
C         NORMALIZE REACTANT MOLE FRACTIONS
C
      DO 6200 K = 1, KK
         REAC(K) = REAC(K) / SUMR
6200  CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE ERRPR (KERR, KEYWRD, LOUT)
      LOGICAL KERR
      CHARACTER *4 KEYWRD
      KERR =.TRUE.
      WRITE (LOUT, '(2X,A,A,A,/)')
     1 ' Error reading data for keyword ', KEYWRD,'.'
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE REDSEN (ICASE, LIGN, LIN, LOUT, LSENS)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*4 KEYWRD
      LOGICAL LSENS
      LSENS = .FALSE.
C
C           ISSUE A PROMPT
C
20    WRITE (LOUT, *) ' Enter keyword: '
      WRITE (LOUT, *) '   SENS, CONP, CONV, VTIM, CONT, TTIM'
C
C....
C         READ NEXT INPUT LINE
C
50    READ  (LIN,  '(A4)')     KEYWRD
      WRITE (LOUT, '(10X,A4)') KEYWRD
      CALL UPCASE (KEYWRD)
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/' .OR.
     1    KEYWRD(1:1) .EQ. '!') GO TO 50
C
C         SENSITIVITY KEYWORD
C
      IF (KEYWRD .EQ. 'SENS') THEN
         WRITE (LOUT, 9000)
         WRITE (LIGN, 9000)
         LSENS = .TRUE.
         GO TO 20
C
C         PROBLEM CHOICE KEYWORD
C
       ELSEIF (KEYWRD .EQ. 'CONP') THEN
         ICASE = 1
       ELSEIF (KEYWRD .EQ. 'CONV') THEN
         ICASE = 2
       ELSEIF (KEYWRD .EQ. 'VTIM') THEN
         ICASE = 3
       ELSEIF (KEYWRD .EQ. 'CONT') THEN
         ICASE = 4
       ELSEIF (KEYWRD .EQ. 'TTIM') THEN
         ICASE = 5
       ELSE
          WRITE (LOUT, 9002) 'Warning, keyword ',KEYWRD, ' not found.'
       ENDIF
C
C      WAS SENSITIVITY FOUND?
C
      IF (.NOT. LSENS) THEN
         WRITE (LOUT, 9001)
         WRITE (LIGN, 9001)
      ENDIF
C
C        CHECK FOR PROBLEM SELECTION
C
      IF (ICASE.GE.1 .AND. ICASE.LE.5) RETURN
      WRITE (LOUT, '(/10X, A, /5X, A, /5X, A)')
     1'Fatal error.',
     2'The problem choice keyword must appear on one of the first',
     3'two input lines (CONP, CONV, CONT, VTIM, or TTIM).'
      STOP
C
9000  FORMAT (/5X, 'Sensitivity analysis will be performed.')
9001  FORMAT (/5X, 'Sensitivity analysis will not be performed.')
9002  FORMAT (5X,A)
9003  FORMAT (5X,A,A,A)
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE TEXT13 (IPAR, KK, KSYM, LIGN, LOUT,
     1                   P, PATM, RPAR, TIM, XMOL, Z)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), RPAR(*), IPAR(*), XMOL(*)
      CHARACTER KSYM(*)*(*)
      DATA NPL /3/
C
C       COMPUTE MOLE FRACTIONS
C
      IPRCK  = IPAR(2)
      IPICK  = IPAR(9)
      CALL CKYTX (Z(2), IPAR(IPICK), RPAR(IPRCK), XMOL)
C
      NPRINT = KK / NPL
      R = MOD(KK,NPL)
      IF (R .NE. 0.) NPRINT = NPRINT + 1
C
C*****vax vms
C      WRITE (LOUT, '(1X,A,1PE11.4,3X,A,1PE11.4)')
C     1 'Time (sec) = ', TIM, 'T (K) = ', Z(1)
C*****END vax vms
C
      PA = P / PATM
      WRITE (LIGN, 7700) TIM, PA, Z(1)
      DO 30 J = 1, NPRINT
         IMIN = NPL*(J-1) + 1
         IMAX = IMIN + NPL-1
         IMAX = MIN (IMAX, KK)
         WRITE (LIGN, 7710)  (KSYM(I), XMOL(I), I = IMIN, IMAX)
30    CONTINUE
C
 7700 FORMAT(/,' t(sec) =',1PE11.4,'   P(atm) =',1PE11.4,
     1'   T(K) =',1PE11.4)
 7710 FORMAT(4(2X,A10,' =',1PE9.2))
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE TEXT13DTS (IPAR, KK, KSYM, LDTS,
     1                   P, PATM, RPAR, TIM, XMOL, Z)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), RPAR(*), IPAR(*), XMOL(*)
      CHARACTER KSYM(*)*(*)
C
C       COMPUTE MOLE FRACTIONS
C
      IPRCK  = IPAR(2)
      IPICK  = IPAR(9)
      CALL CKYTX (Z(2), IPAR(IPICK), RPAR(IPRCK), XMOL)
C
C*****vax vms
C      WRITE (LOUT, '(1X,A,1PE11.4,3X,A,1PE11.4)')
C     1 'Time (sec) = ', TIM, 'T (K) = ', Z(1)
C*****END vax vms
C
      PA = P / PATM
      WRITE (LDTS, *) TIM, PA, Z(1), (XMOL(I), I = 1, KK)
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE TEXT45 (IPAR, KK, KSYM, LIGN, LOUT,
     1                   P, PATM, RPAR, T, TIM, XMOL, Z)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), RPAR(*), IPAR(*), XMOL(*)
      CHARACTER KSYM(*)*(*)
      DATA NPL /4/
C
C       COMPUTE MOLE FRACTIONS
C
      IPRCK  = IPAR(2)
      IPICK  = IPAR(9)
      CALL CKYTX (Z, IPAR(IPICK), RPAR(IPRCK), XMOL)
C
      NPRINT = KK / NPL
      R = MOD(KK,NPL)
      IF (R .NE. 0.) NPRINT = NPRINT + 1
C
C*****vax vms
C      WRITE (LOUT, '(1X,A,1PE11.4,3X,A,1PE11.4)')
C     1 'Time (sec) = ', TIM, 'T (K) = ', T
C*****END vax vms
C
      PA = P / PATM
      WRITE (LIGN, 7700) TIM, PA, T
      DO 30 J = 1, NPRINT
         IMIN = NPL*(J-1) + 1
         IMAX = IMIN + NPL-1
         IMAX = MIN (IMAX, KK)
         WRITE (LIGN, 7710) (KSYM(I), XMOL(I), I = IMIN, IMAX)
30    CONTINUE
C
 7700 FORMAT(/,' t(sec)=',1PE11.4,'   P(atm)=',1PE11.4,
     1'   T(K)=',1PE11.4)
 7710 FORMAT(4(2X,A10,'=',1PE9.2))
      RETURN
      END
C
C----------------------------------------------------------
C
      SUBROUTINE DRES (T, Y, YP, D, IM, RP, IP)
      DIMENSION Y(*), YP(*), D(*), RP(*), IP(*)
C     DUMMY ROUTINE TO AVOID MISSING EXTERNAL
      RETURN
      END
C
      SUBROUTINE UPCASE(STR)
      CHARACTER STR*(*), LCASE(26)*1, UCASE(26)*1
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      DO 10 L = 1, LEN(STR)
         DO 10 N = 1, 26
            IF (STR(L:L) .EQ. LCASE(N)) STR(L:L) = UCASE(N)
   10 CONTINUE
      RETURN
      END
