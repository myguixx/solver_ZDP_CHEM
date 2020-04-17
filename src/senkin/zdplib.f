      SUBROUTINE ZDPINIT()
      use ZDPlasKin

      OPEN (17, FORM='FORMATTED', FILE = 'output/zdpout')

      call ZDPlasKin_init()

      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONP_ZDP (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
      use ZDPlasKin

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*),
     1          wdot_number(species_max), wdot_mole(species_max),
     2          X(species_max)
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

C ---------- ZDPlasKin ----------

      ! get mole fractions
      call CKYTX (Z(2), IPAR(IPICK), RPAR(IPRCK), X)

      ! convert to ZDPlasKin unit system
      R  = 8.31446d+0                        ! Gas constant m3*Pa/K/mole
      AN = 6.02214d+23                       ! Avogadro number 1/mole
      P_Pa = P*0.1d+0                        ! Pa
      density_mole   = 1.0d-6*P_Pa/(R*Z(1))  ! mole/cm3
      density_number = AN*density_mole       ! /cm3

      density(:) = density_number*X(:)       ! /cm3

      ! Set properties
      reduced_field_Td = 180.d0
      gas_temperature_K = Z(1)
      call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature_K,
     1                              REDUCED_FIELD=reduced_field_Td)

      ! get reaction rates of each species [cm-3*s-1]
      call ZDPlasKin_get_rates(SOURCE_TERMS=wdot_number)

      ! convert to CEMKIN unit system
      wdot_mole = wdot_number/AN

      ! check calcrated reaction rates
      ! write(17, *) R, AN, P, (Z(I), I = 1, KK+1)
!       write(17, *) P_Pa, density_mole, density_number, 
!      1             (density(I), I = 1, KK)
!       write(17, *) (wdot_number(I), I = 1, KK)
!       write(17, *) (wdot_mole(I), I = 1, KK)
!       write(17, *) (RPAR(IPWDOT+I), I = 1, KK)

      ! convert to CHEMKIN unit system


      RETURN
      END
C
C---------------------------------------------------------------
C
