module zdp_chem
      implicit none

      integer, parameter :: unit_zdp_chem = 41
      integer, parameter :: unit_dts = 42

      integer,parameter :: num_spec = 83

      ! parameter for dishcarge
      real(8), parameter :: reduced_field_high = 180.0d0     ! [Td]
      real(8), parameter :: frequency          = 30.0d3      ! [Hz]
      real(8), parameter :: duration_freq      = 1/frequency ! [s]
      real(8), parameter :: duration_pulse     = 5.0d-9      ! [s] 
      integer, parameter :: num_pulse          = 10          ! total nuber of pulses [-]
      real(8) reduced_field
      integer ith_pulse

     ! parameter for ODE
      real(8), parameter :: p = 7999.342105d0            ! [Pa]
      real(8), parameter :: delta_t = duration_freq      ! [s]
      real(8), parameter :: tols(4) = (/1.d-8, 1.d-20, 1.d-5, 1.d-5/) ! Tolerances

      ! dependent variables for ODE
      real(8) :: t = 296d0       ! [K]
      real(8) :: y(num_spec) = & ! Mass fractions
            (/0.00d+00,0.00d+00,0.00d+00,5.52d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,3.10d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,1.38d-01, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
              0.00d+00,0.00d+00,5.68d-20/)


      contains 

      subroutine open_units()
            open(unit_zdp_chem, form='formatted', file='output/zdp_chem')
            open(unit_dts, form='formatted', file='output/datasheet')
      end subroutine open_units

      subroutine write_header_datasheet(ksym)
            character(6), intent(in) :: ksym(num_spec)*16

            if (ith_pulse == 1) then
                  write(unit_dts, *) 't(sec)  E/N(Td)  T(K)  ', ksym
            endif

      end subroutine write_header_datasheet

      subroutine write_datasheet(time, z)
            use chemkin, only: int_ckwk, real_ckwk, ipick, iprck

            real(8), intent(in)  :: time
            real(8), intent(in) :: z(num_spec+1)
            
            real(8) x(num_spec)
            real(8) time_ 
            time_ = (ith_pulse - 1)*duration_freq + time

            call ckytx(z(2), int_ckwk(ipick), real_ckwk(iprck), x)

            write(unit_dts, *) time_, reduced_field, z(1), x
      end subroutine write_datasheet

      subroutine calc_reduced_field(time)
            real(8), intent(in)  :: time
      
            if (time < duration_pulse) then
                  reduced_field = reduced_field_high
            else
                  reduced_field = 1.0d-10
            endif
      end subroutine calc_reduced_field

end module zdp_chem

program test_main
      use zdp_chem
      use chemkin, only: initialize_chemkin_workarray, get_next_TY
      use zdplaskin, only: zdplaskin_init
      
      !   ------- initialize section ---------

      call initialize_chemkin_workarray()

      call ZDPlasKin_init()
      
      call open_units()
      
      !   ------- chemistry section ---------

      do ith_pulse = 1, num_pulse

            call get_next_TY(p, t, y, delta_t, tols)

            write(unit_zdp_chem, *) ith_pulse, '-th pulse', t, y

      enddo

end program test_main

      SUBROUTINE TEMPT (TIME, TEMP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      SUBROUTINE VOLT (TIME, VOL, DVDT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      