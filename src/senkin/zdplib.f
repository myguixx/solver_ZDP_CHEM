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
      REAL(8), parameter :: eV_to_erg = 1.602d-12
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

C ---------- get production rates from CHEMKIN library ----------

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

C ---------- calcrate derivetives from CHEMKIN library ----------

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

C ---------- get production rates from ZDPlasKin module ----------

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
      ! reduced_field_Td = 180.d0
      call get_reduced_electricfield(time, reduced_field_Td)
      gas_temperature_K = Z(1)
      call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature_K,
     1                              REDUCED_FIELD=reduced_field_Td)

      ! get reaction rates of each species [cm-3*s-1]
      call ZDPlasKin_get_rates(SOURCE_TERMS=wdot_number)

      ! convert to CEMKIN unit system
      wdot_mole = wdot_number/AN

C ---------- calcrate derivetives from ZDPlasKin module ----------

C
C         ENERGY EQUATION
C

      call ZDPlasKin_get_density('E',density_electron)
      call ZDPlasKin_get_density_total(ALL_SPECIES=density_all)
      call ZDPlasKin_get_conditions(ELEC_POWER_ELASTIC_N=power_elastic)
      density_gas = density_all - density_electron
      DELTA(1) = DELTA(1) + eV_to_erg*power_elastic*density_electron
     1           *density_gas*VOLSP/CPB
C
C         SPECIES EQUATIONS
C
      DO 400 K = 1, KK
         K1 = K-1
         DELTA(K+1) = DELTA(K+1) - wdot_mole(K) *RPAR(IPWT+K1) *VOLSP
 400  CONTINUE

      RETURN
      END
C
C---------------------------------------------------------------
C

      subroutine get_reduced_electricfield(time, reduced_field)
      real(8), parameter :: reduced_field_high = 180.0d0     ! [Td]
      real(8), parameter :: frequency          = 30.0d3      ! [Hz]
      real(8), parameter :: num_pulse          = 300.0d0     ! total nuber of pulses [-]
      real(8), parameter :: duration_freq      = 1/frequency ! [s]
      real(8), parameter :: duration_pulse     = 10.0d-9      ! [s] 
      
      real(8), intent(in)  :: time
      real(8), intent(out) :: reduced_field

      if (time < duration_pulse) then
            reduced_field = reduced_field_high
      else
            reduced_field = 1.0d-10
      endif
      
      ! write(24, *) time, time_mod, ith_pulse, reduced_field

      end subroutine get_reduced_electricfield
