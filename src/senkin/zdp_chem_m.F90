module zdp_chem
    
    implicit none

    integer, parameter :: unit_dts = 41

    REAL(8), parameter :: eV_to_erg = 1.602d-12

    ! contains

    !     subroutine zdpinit()
    !         use zdplaskin

    !         open(unit_zdp_out, form='formatted', file='output/zdpout')

    !         call zdplaskin_init()

    !     end subroutine zdpinit

    !     subroutine rconp_zdp(time, z, zp, delta, ires, rpar, ipar)
    !         use zdplaskin
    !         use chemkin, only: kk, iprck, iprd, ipwt, ipwdot, iph, ipick, &
    !                            unit_stdout, ii

    !         real(8), intent(in) :: time
    !         real(8), intent(in) :: z(:)
    !         real(8), intent(in) :: zp(:)
    !         real(8), intent(inout) :: delta(:)
    !         real(8), intent(in) :: ires(:)
    !         real(8), intent(in) :: rpar(:)
    !         real(8), intent(in) :: ipar(:)

    !         real(8) :: wdot_number(species_max)
    !         real(8) :: wdot_mole(species_max)
    !         real(8) :: X(species_max)

    !         integer i

    !         common /res1/ p

    !         i_work_array = ipar(ipick)
    !         r_work_array = rpar(iprck)

    !         ! ---------- get production rates from CHEMKIN library ----------

    !         do i = 1, ii
    !             call ckrdex(-i, r_work_array, r_work_array(i))
    !         enddo

    !         call ckrhoy(p, z(1), z(2), ipar(ipick), rpar(iprck))

end module zdp_chem