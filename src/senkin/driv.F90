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
      
      ! output transport data
      ! mixture diffusion coefficient [CM**2/S]
      real(8) :: D_mix(num_spec) 
      ! mixture thermal conductivity [ERG/CM*K*S]
      real(8) :: Lambda_mix
      !  mean specific heat at constant pressure [ergs/(gm*K)]
      real(8) c_p
      
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
                  0.00d+00,0.00d+00,0.00d+00/

      ! tolerance values for ODE solver
      data tols_cfd /1.d-8, 1.d-20, 1.d-5, 1.d-5/

      ! flag to generate output
      make_output = .false.

      !   ------- initialize section ---------

      call initialize_chemkin_workarray()
      
      !   ------- chemistry section ---------

      call get_next_TY(p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

      write(6, *) 'temperature [K]'
      write(6, *) t_cfd
      write(6, *) 'mole fractions [-]'
      write(6, *) y_cfd

end program test_main


!       PROGRAM DRIVER
! C
! C*****DOUBLE PRECISION
!       IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
! C*****END DOUBLE PRECISION
! C*****SINGLE PRECISION
! C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
! C*****END SINGLE PRECISION
!       PARAMETER ( LENIWK = 1000000, LENRWK = 20000000, LENCWK = 10000,
!      1            LENSYM = 16)
!       DIMENSION IWORK (LENIWK), RWORK (LENRWK)
!       LOGICAL LEXIST
!       CHARACTER CWORK(LENCWK)*(LENSYM)
!       DATA LIN/15/, LOUT/16/, LINKCK/25/, LSAVE/7/, LIGN/9/, LREST/10/
! C
! C     LIN    = Unit number for Keyword input
! C     LOUT   = Unit number for text output to terminal
! C     LIGN   = Unit number for text output file
! C     LSAVE  = Unit number for binary output file
! C     LINKCK = Unit number for CHENKIN linking file
! C     LREST  = Unit number for binary restart file
! C     LENIWK = Length of integer work array
! C     LENRWK = Length of real work array
! C     LENCWK = Length of character work array
! C     LENSYM = Length of a character string in character work array
! C     IWORK  = Integer work array
! C     RWORK  = Real work array
! C     CWORK  = Character work array
! C
! C*****vms
! C      SET I/O UNITS AND OPEN FILES. OPERATING SYSTEM IS vms VMS.
! C      OPEN (LINKCK, STATUS='OLD', FORM='UNFORMATTED')
! C      OPEN (LSAVE, STATUS='NEW', FORM='UNFORMATTED')
! C      OPEN (LOUT, STATUS='NEW', FORM='FORMATTED')
! C      OPEN (LIGN, STATUS='NEW', FORM='FORMATTED')
! C      OPEN (LIN, STATUS='OLD', FORM='FORMATTED')
! C      INQUIRE (FILE='restart', EXIST=LEXIST)
! C      IF (LEXIST) OPEN (LREST,STATUS='OLD',FORM='UNFORMATTED')
! C*****END vms
! C
! C*****unix
!       OPEN (LINKCK, FORM='UNFORMATTED', FILE='data/cklink')
!       OPEN (LSAVE, FORM='UNFORMATTED', FILE ='output/save')
!       OPEN (LOUT, FORM='FORMATTED', FILE='output/terminalout')
!       OPEN (LIGN, FORM='FORMATTED', FILE = 'output/skout')
!       OPEN (LIN, FORM='FORMATTED', FILE='input/inp')
!       INQUIRE (FILE='output/restart', EXIST=LEXIST)
!       IF (LEXIST) THEN 
!         OPEN (LREST, FORM='UNFORMATTED',
!      1        FILE='output/restart')
!       ENDIF
! C*****END unix
! C
! C     PASS CONTROL TO SENKIN
! C
!       CALL SENKIN (LIN, LOUT, LINKCK, LSAVE, LIGN, LREST,
!      1             LENRWK, RWORK, LENIWK, IWORK, LENCWK, CWORK)
! C
!       STOP
!       END

      SUBROUTINE TEMPT (TIME, TEMP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      SUBROUTINE VOLT (TIME, VOL, DVDT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RETURN
      END

      