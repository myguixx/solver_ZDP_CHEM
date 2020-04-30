program test_main
      use chemkin, only: initialize_chemkin_workarray, get_next_TY
      
      !   ------- user input section ---------

      integer,parameter :: num_spec = 83
      
      ! phisycal values from CFD calculation
      real(8) :: p_cfd = 7999.342105d0  ! Pa
      real(8) :: t_cfd = 296d0         ! K
      real(8) y_cfd(num_spec)           ! Mass fractions
      real(8) :: delta_t_cfd = 1.0d0    ! s
      real(8) :: tols_cfd(4)            ! Tolerances
            
      ! Assurme Y has a same secuence as species in ckout
      data y_cfd /0.00d+00,0.00d+00,0.00d+00,5.52d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,3.10d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,1.38d-01, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,5.68E-20/

      ! tolerance values for ODE solver
      data tols_cfd /1.d-8, 1.d-20, 1.d-5, 1.d-5/

      !   ------- initialize section ---------

      call initialize_chemkin_workarray()
      
      !   ------- chemistry section ---------

      call get_next_TY(p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

      write(6, *) 'temperature [K]'
      write(6, *) t_cfd
      write(6, *) 'mole fractions [-]'
      write(6, *) y_cfd

end program test_main

      SUBROUTINE TEMPT (TIME, TEMP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      SUBROUTINE VOLT (TIME, VOL, DVDT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      