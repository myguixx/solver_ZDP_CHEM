!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Fri Apr 10 08:09:53 2020
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS:  C  H  N  O AR HE KR  E 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 83, species_electrons = 83, species_length = 11, reactions_max = 103, reactions_length = 34
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 3, bolsig_species_length = 3, bolsig_rates_max = 15 
  character(*), parameter, private          :: bolsigfile = "data/bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1,-1/
  data species_name(1:species_max) &
  /"H          ","H2         ","O          ","O2         ","OH         ","H2O        ","N2         ","HO2        ","H2O2       ",&
   "AR         ","HE         ","CO         ","CO2        ","CH2O       ","HCO        ","HOCHO      ","OCHO       ","HOCH2O     ",&
   "CH3OH      ","CH2OH      ","CH3O       ","CH3O2H     ","CH3O2      ","CH4        ","CH3        ","CH2        ","CH2(S)     ",&
   "C          ","CH         ","HOCO       ","HCOH       ","C2H6       ","C2H5       ","C2H4       ","C2H3       ","C2H2       ",&
   "C2H        ","CH3CHO     ","CH3CO      ","CH2CHO     ","CH2CO      ","HCCO       ","H2CC       ","C2O        ","HCCOH      ",&
   "CH3CO3H    ","CH3CO3     ","CH3CO2     ","C2H5OH     ","C2H5O      ","CH2CH2OH   ","CH3CHOH    ","O2C2H4OH   ","C2H5O2H    ",&
   "C2H5O2     ","C2H4O2H    ","C2H4O(1)   ","C2H3O(2)   ","CH2CHOH    ","CH3OCHO    ","CH2OCHO    ","CH3OCO     ","C2H5OCHO   ",&
   "C2H5OCO    ","CH2CH2OCHO ","CH3CHOCHO  ","OH*        ","CH*        ","KR         ","CH4^+      ","CH3^+      ","CH2^+      ",&
   "O2(A1)     ","O2(B1)     ","O2(4.5EV)  ","O(1D)      ","O(1S)      ","O2^+       ","O^+        ","O3         ","HE(SINGLET)",&
   "HE^+       ","E          "/
  data reaction_sign(1:54) &
  /"bolsig:O2->O2(A1)                 ","bolsig:O2->O2(B1)                 ","bolsig:O2->O2(4.5EV)              ",&
   "bolsig:O2->O2(6.0EV)              ","bolsig:O2->O2(8.4EV)              ","bolsig:O2->O2(9.97EV)             ",&
   "bolsig:O2->O2^+                   ","bolsig:HE->HE(SINGLET)            ","bolsig:HE->HE^+                   ",&
   "bolsig:CH4->CH3+H                 ","bolsig:CH4->CH2+H2                ","bolsig:CH4->CH+H2+H               ",&
   "bolsig:CH4->C+2H2                 ","bolsig:CH4->CH4^+                 ","bolsig:CH4->CH3^+                 ",&
   "O2(A1)+CH4=>O2+CH4                ","O2(A1)+CH4=>CH3+HO2               ","O2(B1)+CH4=>O2(A1)+CH4            ",&
   "O2(4.5EV)+CH4=>O2(B1)+CH4         ","O(1D)+CH4=>CH3+OH                 ","O(1D)+CH4=>CH2OH+H                ",&
   "O(1D)+CH4=>CH2O+H2                ","HE(SINGLET)+CH4=>HE+CH+H2+H       ","HE(SINGLET)+CH4=>HE+CH4^++E       ",&
   "HE(SINGLET)+CH4=>HE+CH3^++H+E     ","HE(SINGLET)+CH4=>HE+CH2^++H2+E    ","O2(A1)+O+O2=>O2+O+O2              ",&
   "O2(A1)+O+HE=>O2+O+HE              ","O2(4.5EV)+HE=>O2(A1)+HE           ","O2(4.5EV)+HE=>O2(B1)+HE           ",&
   "O(1D)+HE=>O+HE                    ","O2(A1)+HE=>O2+HE                  ","O2(B1)+HE=>O2(A1)+HE              ",&
   "HE(SINGLET)+HE(SINGLET)=>HE+HE^++E","HE(SINGLET)+O2=>O2^++HE+E         ","HE(SINGLET)+O3=>O2^++O+HE+E       ",&
   "HE(SINGLET)+O2(A1)=>O2^++HE+E     ","HE(SINGLET)+O=>O^++HE+E           ","HE(SINGLET)+O(1D)=>O^++HE+E       ",&
   "HE(SINGLET)+O(1S)=>O^++HE+E       ","E+HE^+=>HE(SINGLET)               ","E+E+HE^+=>HE(SINGLET)+E           ",&
   "E+CH4^+=>CH2+H+H                  ","E+CH4^+=>CH3+H                    ","E+CH3^+=>CH2+H                    ",&
   "E+CH2^+=>CH+H                     ","CH4^++O2=>CH4+O2^+                ","HE^++O2=>O^++O+HE                 ",&
   "HE^++O3=>O^++O2+HE                ","HE^++O2=>O2^++HE                  ","HE^++O2(B1)=>O^++O+HE             ",&
   "HE^++O2(B1)=>O2^++HE              ","HE^++O=>O^++HE                    ","HE^++O(1D)=>O^++HE                "/
  data reaction_sign(55:103) &
  /"HE^++O(1S)=>O^++HE                ","O+O+HE=>O2+HE                     ","O+O+HE=>O2(A1)+HE                 ",&
   "O+O2+HE=>O3+HE                    ","E+O2^+=>O+O                       ","E+O2^+=>O+O(1D)                   ",&
   "E+O2^+=>O(1D)+O(1D)               ","E+O2^+=>O(1D)+O(1S)               ","E+O^++E=>O+E                      ",&
   "E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL  ","O2(A1)=>O2                        ","O2(B1)=>O2(A1)                    ",&
   "O2(B1)=>O2                        ","O2(4.5EV)=>O2                     ","O2(A1)+O=>O2+O                    ",&
   "O2(A1)+O2=>O2+O2                  ","O2(A1)+O3=>O2+O2+O(1D)            ","O2(A1)+O2(A1)=>O2+O2(B1)          ",&
   "O+O3=>O2+O2(A1)                   ","O2(B1)+O=>O2(A1)+O                ","O2(B1)+O=>O2+O(1D)                ",&
   "O2(B1)+O2=>O2(A1)+O2              ","O2(B1)+O3=>O2+O2+O                ","O2(4.5EV)+O=>O2+O(1S)             ",&
   "O2(4.5EV)+O2=>O2(B1)+O2(B1)       ","O(1D)+O=>O+O                      ","O(1D)+O2=>O+O2                    ",&
   "O(1D)+O2=>O+O2(A1)                ","O(1D)+O2=>O+O2(B1)                ","O(1D)+O3=>O2+O+O                  ",&
   "O(1D)+O3=>O2+O2                   ","O(1S)+O=>O(1D)+O                  ","O(1S)+O2=>O(1D)+O2                ",&
   "O(1S)+O2=>O+O+O                   ","O(1S)+O2(A1)=>O+O2(4.5EV)         ","O(1S)+O2(A1)=>O(1D)+O2(B1)        ",&
   "O(1S)+O2(A1)=>O+O+O               ","O(1S)+O3=>O2+O2                   ","O(1S)+O3=>O2+O+O(1D)              ",&
   "O2+O2=>O+O+O2                     ","O2+O=>O+O+O                       ","O3+O2=>O2+O+O2                    ",&
   "O3+O=>O2+O+O                      ","O+O+O2=>O2+O2                     ","O+O+O=>O2+O                       ",&
   "O+O2+O2=>O3+O2                    ","O+O2+O=>O3+O                      ","O^++O2=>O2^++O                    ",&
   "O^++O3=>O2^++O2                   "/
  data bolsig_species(1:bolsig_species_max) &
  /"O2 ","HE ","CH4"/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    EN =", ZDPlasKin_cfg(3)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A11)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A34)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(391x,83(1x,i11))
312 format(A3,1x,A34,1x,83(1x,A11))
313 format(i3,1x,A34,1x,83(1x,1pd11.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(103(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,83(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,103(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,83(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(104(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(04,001) = - reac_rate_local(001) 
  reac_source_local(73,001) = + reac_rate_local(001) 
  reac_source_local(04,002) = - reac_rate_local(002) 
  reac_source_local(74,002) = + reac_rate_local(002) 
  reac_source_local(04,003) = - reac_rate_local(003) 
  reac_source_local(75,003) = + reac_rate_local(003) 
  reac_source_local(03,004) = + reac_rate_local(004) * 2.d0
  reac_source_local(04,004) = - reac_rate_local(004) 
  reac_source_local(03,005) = + reac_rate_local(005) 
  reac_source_local(04,005) = - reac_rate_local(005) 
  reac_source_local(76,005) = + reac_rate_local(005) 
  reac_source_local(03,006) = + reac_rate_local(006) 
  reac_source_local(04,006) = - reac_rate_local(006) 
  reac_source_local(77,006) = + reac_rate_local(006) 
  reac_source_local(04,007) = - reac_rate_local(007) 
  reac_source_local(78,007) = + reac_rate_local(007) 
  reac_source_local(83,007) = + reac_rate_local(007) 
  reac_source_local(11,008) = - reac_rate_local(008) 
  reac_source_local(81,008) = + reac_rate_local(008) 
  reac_source_local(11,009) = - reac_rate_local(009) 
  reac_source_local(82,009) = + reac_rate_local(009) 
  reac_source_local(83,009) = + reac_rate_local(009) 
  reac_source_local(01,010) = + reac_rate_local(010) 
  reac_source_local(24,010) = - reac_rate_local(010) 
  reac_source_local(25,010) = + reac_rate_local(010) 
  reac_source_local(02,011) = + reac_rate_local(011) 
  reac_source_local(24,011) = - reac_rate_local(011) 
  reac_source_local(26,011) = + reac_rate_local(011) 
  reac_source_local(01,012) = + reac_rate_local(012) 
  reac_source_local(02,012) = + reac_rate_local(012) 
  reac_source_local(24,012) = - reac_rate_local(012) 
  reac_source_local(29,012) = + reac_rate_local(012) 
  reac_source_local(02,013) = + reac_rate_local(013) * 2.d0
  reac_source_local(24,013) = - reac_rate_local(013) 
  reac_source_local(28,013) = + reac_rate_local(013) 
  reac_source_local(24,014) = - reac_rate_local(014) 
  reac_source_local(70,014) = + reac_rate_local(014) 
  reac_source_local(83,014) = + reac_rate_local(014) 
  reac_source_local(01,015) = + reac_rate_local(015) 
  reac_source_local(24,015) = - reac_rate_local(015) 
  reac_source_local(71,015) = + reac_rate_local(015) 
  reac_source_local(83,015) = + reac_rate_local(015) 
  reac_source_local(04,016) = + reac_rate_local(016) 
  reac_source_local(73,016) = - reac_rate_local(016) 
  reac_source_local(08,017) = + reac_rate_local(017) 
  reac_source_local(24,017) = - reac_rate_local(017) 
  reac_source_local(25,017) = + reac_rate_local(017) 
  reac_source_local(73,017) = - reac_rate_local(017) 
  reac_source_local(73,018) = + reac_rate_local(018) 
  reac_source_local(74,018) = - reac_rate_local(018) 
  reac_source_local(74,019) = + reac_rate_local(019) 
  reac_source_local(75,019) = - reac_rate_local(019) 
  reac_source_local(05,020) = + reac_rate_local(020) 
  reac_source_local(24,020) = - reac_rate_local(020) 
  reac_source_local(25,020) = + reac_rate_local(020) 
  reac_source_local(76,020) = - reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(20,021) = + reac_rate_local(021) 
  reac_source_local(24,021) = - reac_rate_local(021) 
  reac_source_local(76,021) = - reac_rate_local(021) 
  reac_source_local(02,022) = + reac_rate_local(022) 
  reac_source_local(14,022) = + reac_rate_local(022) 
  reac_source_local(24,022) = - reac_rate_local(022) 
  reac_source_local(76,022) = - reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(02,023) = + reac_rate_local(023) 
  reac_source_local(11,023) = + reac_rate_local(023) 
  reac_source_local(24,023) = - reac_rate_local(023) 
  reac_source_local(29,023) = + reac_rate_local(023) 
  reac_source_local(81,023) = - reac_rate_local(023) 
  reac_source_local(11,024) = + reac_rate_local(024) 
  reac_source_local(24,024) = - reac_rate_local(024) 
  reac_source_local(70,024) = + reac_rate_local(024) 
  reac_source_local(81,024) = - reac_rate_local(024) 
  reac_source_local(83,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(11,025) = + reac_rate_local(025) 
  reac_source_local(24,025) = - reac_rate_local(025) 
  reac_source_local(71,025) = + reac_rate_local(025) 
  reac_source_local(81,025) = - reac_rate_local(025) 
  reac_source_local(83,025) = + reac_rate_local(025) 
  reac_source_local(02,026) = + reac_rate_local(026) 
  reac_source_local(11,026) = + reac_rate_local(026) 
  reac_source_local(24,026) = - reac_rate_local(026) 
  reac_source_local(72,026) = + reac_rate_local(026) 
  reac_source_local(81,026) = - reac_rate_local(026) 
  reac_source_local(83,026) = + reac_rate_local(026) 
  reac_source_local(04,027) = + reac_rate_local(027) 
  reac_source_local(73,027) = - reac_rate_local(027) 
  reac_source_local(04,028) = + reac_rate_local(028) 
  reac_source_local(73,028) = - reac_rate_local(028) 
  reac_source_local(73,029) = + reac_rate_local(029) 
  reac_source_local(75,029) = - reac_rate_local(029) 
  reac_source_local(74,030) = + reac_rate_local(030) 
  reac_source_local(75,030) = - reac_rate_local(030) 
  reac_source_local(03,031) = + reac_rate_local(031) 
  reac_source_local(76,031) = - reac_rate_local(031) 
  reac_source_local(04,032) = + reac_rate_local(032) 
  reac_source_local(73,032) = - reac_rate_local(032) 
  reac_source_local(73,033) = + reac_rate_local(033) 
  reac_source_local(74,033) = - reac_rate_local(033) 
  reac_source_local(11,034) = + reac_rate_local(034) 
  reac_source_local(81,034) = - reac_rate_local(034) * 2.d0
  reac_source_local(82,034) = + reac_rate_local(034) 
  reac_source_local(83,034) = + reac_rate_local(034) 
  reac_source_local(04,035) = - reac_rate_local(035) 
  reac_source_local(11,035) = + reac_rate_local(035) 
  reac_source_local(78,035) = + reac_rate_local(035) 
  reac_source_local(81,035) = - reac_rate_local(035) 
  reac_source_local(83,035) = + reac_rate_local(035) 
  reac_source_local(03,036) = + reac_rate_local(036) 
  reac_source_local(11,036) = + reac_rate_local(036) 
  reac_source_local(78,036) = + reac_rate_local(036) 
  reac_source_local(80,036) = - reac_rate_local(036) 
  reac_source_local(81,036) = - reac_rate_local(036) 
  reac_source_local(83,036) = + reac_rate_local(036) 
  reac_source_local(11,037) = + reac_rate_local(037) 
  reac_source_local(73,037) = - reac_rate_local(037) 
  reac_source_local(78,037) = + reac_rate_local(037) 
  reac_source_local(81,037) = - reac_rate_local(037) 
  reac_source_local(83,037) = + reac_rate_local(037) 
  reac_source_local(03,038) = - reac_rate_local(038) 
  reac_source_local(11,038) = + reac_rate_local(038) 
  reac_source_local(79,038) = + reac_rate_local(038) 
  reac_source_local(81,038) = - reac_rate_local(038) 
  reac_source_local(83,038) = + reac_rate_local(038) 
  reac_source_local(11,039) = + reac_rate_local(039) 
  reac_source_local(76,039) = - reac_rate_local(039) 
  reac_source_local(79,039) = + reac_rate_local(039) 
  reac_source_local(81,039) = - reac_rate_local(039) 
  reac_source_local(83,039) = + reac_rate_local(039) 
  reac_source_local(11,040) = + reac_rate_local(040) 
  reac_source_local(77,040) = - reac_rate_local(040) 
  reac_source_local(79,040) = + reac_rate_local(040) 
  reac_source_local(81,040) = - reac_rate_local(040) 
  reac_source_local(83,040) = + reac_rate_local(040) 
  reac_source_local(81,041) = + reac_rate_local(041) 
  reac_source_local(82,041) = - reac_rate_local(041) 
  reac_source_local(83,041) = - reac_rate_local(041) 
  reac_source_local(81,042) = + reac_rate_local(042) 
  reac_source_local(82,042) = - reac_rate_local(042) 
  reac_source_local(83,042) = - reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) * 2.d0
  reac_source_local(26,043) = + reac_rate_local(043) 
  reac_source_local(70,043) = - reac_rate_local(043) 
  reac_source_local(83,043) = - reac_rate_local(043) 
  reac_source_local(01,044) = + reac_rate_local(044) 
  reac_source_local(25,044) = + reac_rate_local(044) 
  reac_source_local(70,044) = - reac_rate_local(044) 
  reac_source_local(83,044) = - reac_rate_local(044) 
  reac_source_local(01,045) = + reac_rate_local(045) 
  reac_source_local(26,045) = + reac_rate_local(045) 
  reac_source_local(71,045) = - reac_rate_local(045) 
  reac_source_local(83,045) = - reac_rate_local(045) 
  reac_source_local(01,046) = + reac_rate_local(046) 
  reac_source_local(29,046) = + reac_rate_local(046) 
  reac_source_local(72,046) = - reac_rate_local(046) 
  reac_source_local(83,046) = - reac_rate_local(046) 
  reac_source_local(04,047) = - reac_rate_local(047) 
  reac_source_local(24,047) = + reac_rate_local(047) 
  reac_source_local(70,047) = - reac_rate_local(047) 
  reac_source_local(78,047) = + reac_rate_local(047) 
  reac_source_local(03,048) = + reac_rate_local(048) 
  reac_source_local(04,048) = - reac_rate_local(048) 
  reac_source_local(11,048) = + reac_rate_local(048) 
  reac_source_local(79,048) = + reac_rate_local(048) 
  reac_source_local(82,048) = - reac_rate_local(048) 
  reac_source_local(04,049) = + reac_rate_local(049) 
  reac_source_local(11,049) = + reac_rate_local(049) 
  reac_source_local(79,049) = + reac_rate_local(049) 
  reac_source_local(80,049) = - reac_rate_local(049) 
  reac_source_local(82,049) = - reac_rate_local(049) 
  reac_source_local(04,050) = - reac_rate_local(050) 
  reac_source_local(11,050) = + reac_rate_local(050) 
  reac_source_local(78,050) = + reac_rate_local(050) 
  reac_source_local(82,050) = - reac_rate_local(050) 
  reac_source_local(03,051) = + reac_rate_local(051) 
  reac_source_local(11,051) = + reac_rate_local(051) 
  reac_source_local(74,051) = - reac_rate_local(051) 
  reac_source_local(79,051) = + reac_rate_local(051) 
  reac_source_local(82,051) = - reac_rate_local(051) 
  reac_source_local(11,052) = + reac_rate_local(052) 
  reac_source_local(74,052) = - reac_rate_local(052) 
  reac_source_local(78,052) = + reac_rate_local(052) 
  reac_source_local(82,052) = - reac_rate_local(052) 
  reac_source_local(03,053) = - reac_rate_local(053) 
  reac_source_local(11,053) = + reac_rate_local(053) 
  reac_source_local(79,053) = + reac_rate_local(053) 
  reac_source_local(82,053) = - reac_rate_local(053) 
  reac_source_local(11,054) = + reac_rate_local(054) 
  reac_source_local(76,054) = - reac_rate_local(054) 
  reac_source_local(79,054) = + reac_rate_local(054) 
  reac_source_local(82,054) = - reac_rate_local(054) 
  reac_source_local(11,055) = + reac_rate_local(055) 
  reac_source_local(77,055) = - reac_rate_local(055) 
  reac_source_local(79,055) = + reac_rate_local(055) 
  reac_source_local(82,055) = - reac_rate_local(055) 
  reac_source_local(03,056) = - reac_rate_local(056) * 2.d0
  reac_source_local(04,056) = + reac_rate_local(056) 
  reac_source_local(03,057) = - reac_rate_local(057) * 2.d0
  reac_source_local(73,057) = + reac_rate_local(057) 
  reac_source_local(03,058) = - reac_rate_local(058) 
  reac_source_local(04,058) = - reac_rate_local(058) 
  reac_source_local(80,058) = + reac_rate_local(058) 
  reac_source_local(03,059) = + reac_rate_local(059) * 2.d0
  reac_source_local(78,059) = - reac_rate_local(059) 
  reac_source_local(83,059) = - reac_rate_local(059) 
  reac_source_local(03,060) = + reac_rate_local(060) 
  reac_source_local(76,060) = + reac_rate_local(060) 
  reac_source_local(78,060) = - reac_rate_local(060) 
  reac_source_local(83,060) = - reac_rate_local(060) 
  reac_source_local(76,061) = + reac_rate_local(061) * 2.d0
  reac_source_local(78,061) = - reac_rate_local(061) 
  reac_source_local(83,061) = - reac_rate_local(061) 
  reac_source_local(76,062) = + reac_rate_local(062) 
  reac_source_local(77,062) = + reac_rate_local(062) 
  reac_source_local(78,062) = - reac_rate_local(062) 
  reac_source_local(83,062) = - reac_rate_local(062) 
  reac_source_local(03,063) = + reac_rate_local(063) 
  reac_source_local(79,063) = - reac_rate_local(063) 
  reac_source_local(83,063) = - reac_rate_local(063) 
  reac_source_local(03,064) = + reac_rate_local(064) 
  reac_source_local(79,064) = - reac_rate_local(064) 
  reac_source_local(83,064) = - reac_rate_local(064) 
  reac_source_local(04,065) = + reac_rate_local(065) 
  reac_source_local(73,065) = - reac_rate_local(065) 
  reac_source_local(73,066) = + reac_rate_local(066) 
  reac_source_local(74,066) = - reac_rate_local(066) 
  reac_source_local(04,067) = + reac_rate_local(067) 
  reac_source_local(74,067) = - reac_rate_local(067) 
  reac_source_local(04,068) = + reac_rate_local(068) 
  reac_source_local(75,068) = - reac_rate_local(068) 
  reac_source_local(04,069) = + reac_rate_local(069) 
  reac_source_local(73,069) = - reac_rate_local(069) 
  reac_source_local(04,070) = + reac_rate_local(070) 
  reac_source_local(73,070) = - reac_rate_local(070) 
  reac_source_local(04,071) = + reac_rate_local(071) * 2.d0
  reac_source_local(73,071) = - reac_rate_local(071) 
  reac_source_local(76,071) = + reac_rate_local(071) 
  reac_source_local(80,071) = - reac_rate_local(071) 
  reac_source_local(04,072) = + reac_rate_local(072) 
  reac_source_local(73,072) = - reac_rate_local(072) * 2.d0
  reac_source_local(74,072) = + reac_rate_local(072) 
  reac_source_local(03,073) = - reac_rate_local(073) 
  reac_source_local(04,073) = + reac_rate_local(073) 
  reac_source_local(73,073) = + reac_rate_local(073) 
  reac_source_local(80,073) = - reac_rate_local(073) 
  reac_source_local(73,074) = + reac_rate_local(074) 
  reac_source_local(74,074) = - reac_rate_local(074) 
  reac_source_local(03,075) = - reac_rate_local(075) 
  reac_source_local(04,075) = + reac_rate_local(075) 
  reac_source_local(74,075) = - reac_rate_local(075) 
  reac_source_local(76,075) = + reac_rate_local(075) 
  reac_source_local(73,076) = + reac_rate_local(076) 
  reac_source_local(74,076) = - reac_rate_local(076) 
  reac_source_local(03,077) = + reac_rate_local(077) 
  reac_source_local(04,077) = + reac_rate_local(077) * 2.d0
  reac_source_local(74,077) = - reac_rate_local(077) 
  reac_source_local(80,077) = - reac_rate_local(077) 
  reac_source_local(03,078) = - reac_rate_local(078) 
  reac_source_local(04,078) = + reac_rate_local(078) 
  reac_source_local(75,078) = - reac_rate_local(078) 
  reac_source_local(77,078) = + reac_rate_local(078) 
  reac_source_local(04,079) = - reac_rate_local(079) 
  reac_source_local(74,079) = + reac_rate_local(079) * 2.d0
  reac_source_local(75,079) = - reac_rate_local(079) 
  reac_source_local(03,080) = + reac_rate_local(080) 
  reac_source_local(76,080) = - reac_rate_local(080) 
  reac_source_local(03,081) = + reac_rate_local(081) 
  reac_source_local(76,081) = - reac_rate_local(081) 
  reac_source_local(03,082) = + reac_rate_local(082) 
  reac_source_local(04,082) = - reac_rate_local(082) 
  reac_source_local(73,082) = + reac_rate_local(082) 
  reac_source_local(76,082) = - reac_rate_local(082) 
  reac_source_local(03,083) = + reac_rate_local(083) 
  reac_source_local(04,083) = - reac_rate_local(083) 
  reac_source_local(74,083) = + reac_rate_local(083) 
  reac_source_local(76,083) = - reac_rate_local(083) 
  reac_source_local(03,084) = + reac_rate_local(084) * 2.d0
  reac_source_local(04,084) = + reac_rate_local(084) 
  reac_source_local(76,084) = - reac_rate_local(084) 
  reac_source_local(80,084) = - reac_rate_local(084) 
  reac_source_local(04,085) = + reac_rate_local(085) * 2.d0
  reac_source_local(76,085) = - reac_rate_local(085) 
  reac_source_local(80,085) = - reac_rate_local(085) 
  reac_source_local(76,086) = + reac_rate_local(086) 
  reac_source_local(77,086) = - reac_rate_local(086) 
  reac_source_local(76,087) = + reac_rate_local(087) 
  reac_source_local(77,087) = - reac_rate_local(087) 
  reac_source_local(03,088) = + reac_rate_local(088) * 3.d0
  reac_source_local(04,088) = - reac_rate_local(088) 
  reac_source_local(77,088) = - reac_rate_local(088) 
  reac_source_local(03,089) = + reac_rate_local(089) 
  reac_source_local(73,089) = - reac_rate_local(089) 
  reac_source_local(75,089) = + reac_rate_local(089) 
  reac_source_local(77,089) = - reac_rate_local(089) 
  reac_source_local(73,090) = - reac_rate_local(090) 
  reac_source_local(74,090) = + reac_rate_local(090) 
  reac_source_local(76,090) = + reac_rate_local(090) 
  reac_source_local(77,090) = - reac_rate_local(090) 
  reac_source_local(03,091) = + reac_rate_local(091) * 3.d0
  reac_source_local(73,091) = - reac_rate_local(091) 
  reac_source_local(77,091) = - reac_rate_local(091) 
  reac_source_local(04,092) = + reac_rate_local(092) * 2.d0
  reac_source_local(77,092) = - reac_rate_local(092) 
  reac_source_local(80,092) = - reac_rate_local(092) 
  reac_source_local(03,093) = + reac_rate_local(093) 
  reac_source_local(04,093) = + reac_rate_local(093) 
  reac_source_local(76,093) = + reac_rate_local(093) 
  reac_source_local(77,093) = - reac_rate_local(093) 
  reac_source_local(80,093) = - reac_rate_local(093) 
  reac_source_local(03,094) = + reac_rate_local(094) * 2.d0
  reac_source_local(04,094) = - reac_rate_local(094) 
  reac_source_local(03,095) = + reac_rate_local(095) * 2.d0
  reac_source_local(04,095) = - reac_rate_local(095) 
  reac_source_local(03,096) = + reac_rate_local(096) 
  reac_source_local(04,096) = + reac_rate_local(096) 
  reac_source_local(80,096) = - reac_rate_local(096) 
  reac_source_local(03,097) = + reac_rate_local(097) 
  reac_source_local(04,097) = + reac_rate_local(097) 
  reac_source_local(80,097) = - reac_rate_local(097) 
  reac_source_local(03,098) = - reac_rate_local(098) * 2.d0
  reac_source_local(04,098) = + reac_rate_local(098) 
  reac_source_local(03,099) = - reac_rate_local(099) * 2.d0
  reac_source_local(04,099) = + reac_rate_local(099) 
  reac_source_local(03,100) = - reac_rate_local(100) 
  reac_source_local(04,100) = - reac_rate_local(100) 
  reac_source_local(80,100) = + reac_rate_local(100) 
  reac_source_local(03,101) = - reac_rate_local(101) 
  reac_source_local(04,101) = - reac_rate_local(101) 
  reac_source_local(80,101) = + reac_rate_local(101) 
  reac_source_local(03,102) = + reac_rate_local(102) 
  reac_source_local(04,102) = - reac_rate_local(102) 
  reac_source_local(78,102) = + reac_rate_local(102) 
  reac_source_local(79,102) = - reac_rate_local(102) 
  reac_source_local(04,103) = + reac_rate_local(103) 
  reac_source_local(78,103) = + reac_rate_local(103) 
  reac_source_local(79,103) = - reac_rate_local(103) 
  reac_source_local(80,103) = - reac_rate_local(103) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(84)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(04) * density(83) 
  rrt(002) = rrt(002) * density(04) * density(83) 
  rrt(003) = rrt(003) * density(04) * density(83) 
  rrt(004) = rrt(004) * density(04) * density(83) 
  rrt(005) = rrt(005) * density(04) * density(83) 
  rrt(006) = rrt(006) * density(04) * density(83) 
  rrt(007) = rrt(007) * density(04) * density(83) 
  rrt(008) = rrt(008) * density(11) * density(83) 
  rrt(009) = rrt(009) * density(11) * density(83) 
  rrt(010) = rrt(010) * density(24) * density(83) 
  rrt(011) = rrt(011) * density(24) * density(83) 
  rrt(012) = rrt(012) * density(24) * density(83) 
  rrt(013) = rrt(013) * density(24) * density(83) 
  rrt(014) = rrt(014) * density(24) * density(83) 
  rrt(015) = rrt(015) * density(24) * density(83) 
  rrt(016) = rrt(016) * density(24) * density(73) 
  rrt(017) = rrt(017) * density(24) * density(73) 
  rrt(018) = rrt(018) * density(24) * density(74) 
  rrt(019) = rrt(019) * density(24) * density(75) 
  rrt(020) = rrt(020) * density(24) * density(76) 
  rrt(021) = rrt(021) * density(24) * density(76) 
  rrt(022) = rrt(022) * density(24) * density(76) 
  rrt(023) = rrt(023) * density(24) * density(81) 
  rrt(024) = rrt(024) * density(24) * density(81) 
  rrt(025) = rrt(025) * density(24) * density(81) 
  rrt(026) = rrt(026) * density(24) * density(81) 
  rrt(027) = rrt(027) * density(03) * density(04) * density(73) 
  rrt(028) = rrt(028) * density(03) * density(11) * density(73) 
  rrt(029) = rrt(029) * density(11) * density(75) 
  rrt(030) = rrt(030) * density(11) * density(75) 
  rrt(031) = rrt(031) * density(11) * density(76) 
  rrt(032) = rrt(032) * density(11) * density(73) 
  rrt(033) = rrt(033) * density(11) * density(74) 
  rrt(034) = rrt(034) * density(81)**2 
  rrt(035) = rrt(035) * density(04) * density(81) 
  rrt(036) = rrt(036) * density(80) * density(81) 
  rrt(037) = rrt(037) * density(73) * density(81) 
  rrt(038) = rrt(038) * density(03) * density(81) 
  rrt(039) = rrt(039) * density(76) * density(81) 
  rrt(040) = rrt(040) * density(77) * density(81) 
  rrt(041) = rrt(041) * density(82) * density(83) 
  rrt(042) = rrt(042) * density(82) * density(83)**2 
  rrt(043) = rrt(043) * density(70) * density(83) 
  rrt(044) = rrt(044) * density(70) * density(83) 
  rrt(045) = rrt(045) * density(71) * density(83) 
  rrt(046) = rrt(046) * density(72) * density(83) 
  rrt(047) = rrt(047) * density(04) * density(70) 
  rrt(048) = rrt(048) * density(04) * density(82) 
  rrt(049) = rrt(049) * density(80) * density(82) 
  rrt(050) = rrt(050) * density(04) * density(82) 
  rrt(051) = rrt(051) * density(74) * density(82) 
  rrt(052) = rrt(052) * density(74) * density(82) 
  rrt(053) = rrt(053) * density(03) * density(82) 
  rrt(054) = rrt(054) * density(76) * density(82) 
  rrt(055) = rrt(055) * density(77) * density(82) 
  rrt(056) = rrt(056) * density(03)**2 * density(11) 
  rrt(057) = rrt(057) * density(03)**2 * density(11) 
  rrt(058) = rrt(058) * density(03) * density(04) * density(11) 
  rrt(059) = rrt(059) * density(78) * density(83) 
  rrt(060) = rrt(060) * density(78) * density(83) 
  rrt(061) = rrt(061) * density(78) * density(83) 
  rrt(062) = rrt(062) * density(78) * density(83) 
  rrt(063) = rrt(063) * density(79) * density(83)**2 
  rrt(064) = rrt(064) * density(79) * density(83) 
  rrt(065) = rrt(065) * density(73) 
  rrt(066) = rrt(066) * density(74) 
  rrt(067) = rrt(067) * density(74) 
  rrt(068) = rrt(068) * density(75) 
  rrt(069) = rrt(069) * density(03) * density(73) 
  rrt(070) = rrt(070) * density(04) * density(73) 
  rrt(071) = rrt(071) * density(73) * density(80) 
  rrt(072) = rrt(072) * density(73)**2 
  rrt(073) = rrt(073) * density(03) * density(80) 
  rrt(074) = rrt(074) * density(03) * density(74) 
  rrt(075) = rrt(075) * density(03) * density(74) 
  rrt(076) = rrt(076) * density(04) * density(74) 
  rrt(077) = rrt(077) * density(74) * density(80) 
  rrt(078) = rrt(078) * density(03) * density(75) 
  rrt(079) = rrt(079) * density(04) * density(75) 
  rrt(080) = rrt(080) * density(03) * density(76) 
  rrt(081) = rrt(081) * density(04) * density(76) 
  rrt(082) = rrt(082) * density(04) * density(76) 
  rrt(083) = rrt(083) * density(04) * density(76) 
  rrt(084) = rrt(084) * density(76) * density(80) 
  rrt(085) = rrt(085) * density(76) * density(80) 
  rrt(086) = rrt(086) * density(03) * density(77) 
  rrt(087) = rrt(087) * density(04) * density(77) 
  rrt(088) = rrt(088) * density(04) * density(77) 
  rrt(089) = rrt(089) * density(73) * density(77) 
  rrt(090) = rrt(090) * density(73) * density(77) 
  rrt(091) = rrt(091) * density(73) * density(77) 
  rrt(092) = rrt(092) * density(77) * density(80) 
  rrt(093) = rrt(093) * density(77) * density(80) 
  rrt(094) = rrt(094) * density(04)**2 
  rrt(095) = rrt(095) * density(03) * density(04) 
  rrt(096) = rrt(096) * density(04) * density(80) 
  rrt(097) = rrt(097) * density(03) * density(80) 
  rrt(098) = rrt(098) * density(03)**2 * density(04) 
  rrt(099) = rrt(099) * density(03)**3 
  rrt(100) = rrt(100) * density(03) * density(04)**2 
  rrt(101) = rrt(101) * density(03)**2 * density(04) 
  rrt(102) = rrt(102) * density(04) * density(79) 
  rrt(103) = rrt(103) * density(79) * density(80) 
  ydot(01) = +rrt(010)+rrt(012)+rrt(015)+rrt(021)+rrt(023)+rrt(025)+  2.d0 * rrt(043)+rrt(044)+rrt(045)+rrt(046) 
  ydot(02) = +rrt(011)+rrt(012)+  2.d0 * rrt(013)+rrt(022)+rrt(023)+rrt(026) 
  ydot(03) = +  2.d0 * rrt(004)+rrt(005)+rrt(006)+rrt(031)+rrt(036)-rrt(038)+rrt(048)+rrt(051)-rrt(053)-  2.d0 * rrt(056)&
             -  2.d0 * rrt(057)-rrt(058)+  2.d0 * rrt(059)+rrt(060)+rrt(063)+rrt(064)-rrt(073)-rrt(075)+rrt(077)-rrt(078)+rrt(080)&
             +rrt(081)+rrt(082)+rrt(083)+  2.d0 * rrt(084)+  3.d0 * rrt(088)+rrt(089)+  3.d0 * rrt(091)+rrt(093)+  2.d0 * rrt(094)&
             +  2.d0 * rrt(095)+rrt(096)+rrt(097)-  2.d0 * rrt(098)-  2.d0 * rrt(099)-rrt(100)-rrt(101)+rrt(102) 
  ydot(04) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)+rrt(016)+rrt(027)+rrt(028)+rrt(032)-rrt(035)-rrt(047)&
             -rrt(048)+rrt(049)-rrt(050)+rrt(056)-rrt(058)+rrt(065)+rrt(067)+rrt(068)+rrt(069)+rrt(070)+  2.d0 * rrt(071)+rrt(072)&
             +rrt(073)+rrt(075)+  2.d0 * rrt(077)+rrt(078)-rrt(079)-rrt(082)-rrt(083)+rrt(084)+  2.d0 * rrt(085)-rrt(088)&
             +  2.d0 * rrt(092)+rrt(093)-rrt(094)-rrt(095)+rrt(096)+rrt(097)+rrt(098)+rrt(099)-rrt(100)-rrt(101)-rrt(102)+rrt(103) 
  ydot(05) = +rrt(020) 
  ydot(06) = 0.0d0
  ydot(07) = 0.0d0
  ydot(08) = +rrt(017) 
  ydot(09) = 0.0d0
  ydot(10) = 0.0d0
  ydot(11) = -rrt(008)-rrt(009)+rrt(023)+rrt(024)+rrt(025)+rrt(026)+rrt(034)+rrt(035)+rrt(036)+rrt(037)+rrt(038)+rrt(039)+rrt(040)&
             +rrt(048)+rrt(049)+rrt(050)+rrt(051)+rrt(052)+rrt(053)+rrt(054)+rrt(055) 
  ydot(12) = 0.0d0
  ydot(13) = 0.0d0
  ydot(14) = +rrt(022) 
  ydot(15) = 0.0d0
  ydot(16) = 0.0d0
  ydot(17) = 0.0d0
  ydot(18) = 0.0d0
  ydot(19) = 0.0d0
  ydot(20) = +rrt(021) 
  ydot(21) = 0.0d0
  ydot(22) = 0.0d0
  ydot(23) = 0.0d0
  ydot(24) = -rrt(010)-rrt(011)-rrt(012)-rrt(013)-rrt(014)-rrt(015)-rrt(017)-rrt(020)-rrt(021)-rrt(022)-rrt(023)-rrt(024)-rrt(025)&
             -rrt(026)+rrt(047) 
  ydot(25) = +rrt(010)+rrt(017)+rrt(020)+rrt(044) 
  ydot(26) = +rrt(011)+rrt(043)+rrt(045) 
  ydot(27) = 0.0d0
  ydot(28) = +rrt(013) 
  ydot(29) = +rrt(012)+rrt(023)+rrt(046) 
  ydot(30) = 0.0d0
  ydot(31) = 0.0d0
  ydot(32) = 0.0d0
  ydot(33) = 0.0d0
  ydot(34) = 0.0d0
  ydot(35) = 0.0d0
  ydot(36) = 0.0d0
  ydot(37) = 0.0d0
  ydot(38) = 0.0d0
  ydot(39) = 0.0d0
  ydot(40) = 0.0d0
  ydot(41) = 0.0d0
  ydot(42) = 0.0d0
  ydot(43) = 0.0d0
  ydot(44) = 0.0d0
  ydot(45) = 0.0d0
  ydot(46) = 0.0d0
  ydot(47) = 0.0d0
  ydot(48) = 0.0d0
  ydot(49) = 0.0d0
  ydot(50) = 0.0d0
  ydot(51) = 0.0d0
  ydot(52) = 0.0d0
  ydot(53) = 0.0d0
  ydot(54) = 0.0d0
  ydot(55) = 0.0d0
  ydot(56) = 0.0d0
  ydot(57) = 0.0d0
  ydot(58) = 0.0d0
  ydot(59) = 0.0d0
  ydot(60) = 0.0d0
  ydot(61) = 0.0d0
  ydot(62) = 0.0d0
  ydot(63) = 0.0d0
  ydot(64) = 0.0d0
  ydot(65) = 0.0d0
  ydot(66) = 0.0d0
  ydot(67) = 0.0d0
  ydot(68) = 0.0d0
  ydot(69) = 0.0d0
  ydot(70) = +rrt(014)+rrt(024)-rrt(043)-rrt(044)-rrt(047) 
  ydot(71) = +rrt(015)+rrt(025)-rrt(045) 
  ydot(72) = +rrt(026)-rrt(046) 
  ydot(73) = +rrt(001)-rrt(016)-rrt(017)+rrt(018)-rrt(027)-rrt(028)+rrt(029)-rrt(032)+rrt(033)-rrt(037)+rrt(057)-rrt(065)+rrt(066)&
             -rrt(069)-rrt(070)-rrt(071)-  2.d0 * rrt(072)+rrt(073)+rrt(074)+rrt(076)+rrt(082)-rrt(089)-rrt(090)-rrt(091) 
  ydot(74) = +rrt(002)-rrt(018)+rrt(019)+rrt(030)-rrt(033)-rrt(051)-rrt(052)-rrt(066)-rrt(067)+rrt(072)-rrt(074)-rrt(075)-rrt(076)&
             -rrt(077)+  2.d0 * rrt(079)+rrt(083)+rrt(090) 
  ydot(75) = +rrt(003)-rrt(019)-rrt(029)-rrt(030)-rrt(068)-rrt(078)-rrt(079)+rrt(089) 
  ydot(76) = +rrt(005)-rrt(020)-rrt(021)-rrt(022)-rrt(031)-rrt(039)-rrt(054)+rrt(060)+  2.d0 * rrt(061)+rrt(062)+rrt(071)+rrt(075)&
             -rrt(080)-rrt(081)-rrt(082)-rrt(083)-rrt(084)-rrt(085)+rrt(086)+rrt(087)+rrt(090)+rrt(093) 
  ydot(77) = +rrt(006)-rrt(040)-rrt(055)+rrt(062)+rrt(078)-rrt(086)-rrt(087)-rrt(088)-rrt(089)-rrt(090)-rrt(091)-rrt(092)-rrt(093) 
  ydot(78) = +rrt(007)+rrt(035)+rrt(036)+rrt(037)+rrt(047)+rrt(050)+rrt(052)-rrt(059)-rrt(060)-rrt(061)-rrt(062)+rrt(102)+rrt(103) 
  ydot(79) = +rrt(038)+rrt(039)+rrt(040)+rrt(048)+rrt(049)+rrt(051)+rrt(053)+rrt(054)+rrt(055)-rrt(063)-rrt(064)-rrt(102)-rrt(103) 
  ydot(80) = -rrt(036)-rrt(049)+rrt(058)-rrt(071)-rrt(073)-rrt(077)-rrt(084)-rrt(085)-rrt(092)-rrt(093)-rrt(096)-rrt(097)+rrt(100)&
             +rrt(101)-rrt(103) 
  ydot(81) = +rrt(008)-rrt(023)-rrt(024)-rrt(025)-rrt(026)-  2.d0 * rrt(034)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)&
             +rrt(041)+rrt(042) 
  ydot(82) = +rrt(009)+rrt(034)-rrt(041)-rrt(042)-rrt(048)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)-rrt(054)-rrt(055) 
  ydot(83) = +rrt(007)+rrt(009)+rrt(014)+rrt(015)+rrt(024)+rrt(025)+rrt(026)+rrt(034)+rrt(035)+rrt(036)+rrt(037)+rrt(038)+rrt(039)&
             +rrt(040)-rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)-rrt(046)-rrt(059)-rrt(060)-rrt(061)-rrt(062)-rrt(063)-rrt(064) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(84) = 0.0d0
  if( lgas_heating ) then
    ydot(84) = ( ZDPlasKin_cfg(14)/k_B + ydot(84) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(84) = ydot(84) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(84)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(04,04) = pd(04,04) - rrt(001) * density(83) 
  pd(04,83) = pd(04,83) - rrt(001) * density(04) 
  pd(73,04) = pd(73,04) + rrt(001) * density(83) 
  pd(73,83) = pd(73,83) + rrt(001) * density(04) 
  pd(04,04) = pd(04,04) - rrt(002) * density(83) 
  pd(04,83) = pd(04,83) - rrt(002) * density(04) 
  pd(74,04) = pd(74,04) + rrt(002) * density(83) 
  pd(74,83) = pd(74,83) + rrt(002) * density(04) 
  pd(04,04) = pd(04,04) - rrt(003) * density(83) 
  pd(04,83) = pd(04,83) - rrt(003) * density(04) 
  pd(75,04) = pd(75,04) + rrt(003) * density(83) 
  pd(75,83) = pd(75,83) + rrt(003) * density(04) 
  pd(03,04) = pd(03,04) + rrt(004) * density(83) * 2.0d0
  pd(03,83) = pd(03,83) + rrt(004) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(004) * density(83) 
  pd(04,83) = pd(04,83) - rrt(004) * density(04) 
  pd(03,04) = pd(03,04) + rrt(005) * density(83) 
  pd(03,83) = pd(03,83) + rrt(005) * density(04) 
  pd(04,04) = pd(04,04) - rrt(005) * density(83) 
  pd(04,83) = pd(04,83) - rrt(005) * density(04) 
  pd(76,04) = pd(76,04) + rrt(005) * density(83) 
  pd(76,83) = pd(76,83) + rrt(005) * density(04) 
  pd(03,04) = pd(03,04) + rrt(006) * density(83) 
  pd(03,83) = pd(03,83) + rrt(006) * density(04) 
  pd(04,04) = pd(04,04) - rrt(006) * density(83) 
  pd(04,83) = pd(04,83) - rrt(006) * density(04) 
  pd(77,04) = pd(77,04) + rrt(006) * density(83) 
  pd(77,83) = pd(77,83) + rrt(006) * density(04) 
  pd(04,04) = pd(04,04) - rrt(007) * density(83) 
  pd(04,83) = pd(04,83) - rrt(007) * density(04) 
  pd(78,04) = pd(78,04) + rrt(007) * density(83) 
  pd(78,83) = pd(78,83) + rrt(007) * density(04) 
  pd(83,04) = pd(83,04) + rrt(007) * density(83) 
  pd(83,83) = pd(83,83) + rrt(007) * density(04) 
  pd(11,11) = pd(11,11) - rrt(008) * density(83) 
  pd(11,83) = pd(11,83) - rrt(008) * density(11) 
  pd(81,11) = pd(81,11) + rrt(008) * density(83) 
  pd(81,83) = pd(81,83) + rrt(008) * density(11) 
  pd(11,11) = pd(11,11) - rrt(009) * density(83) 
  pd(11,83) = pd(11,83) - rrt(009) * density(11) 
  pd(82,11) = pd(82,11) + rrt(009) * density(83) 
  pd(82,83) = pd(82,83) + rrt(009) * density(11) 
  pd(83,11) = pd(83,11) + rrt(009) * density(83) 
  pd(83,83) = pd(83,83) + rrt(009) * density(11) 
  pd(01,24) = pd(01,24) + rrt(010) * density(83) 
  pd(01,83) = pd(01,83) + rrt(010) * density(24) 
  pd(24,24) = pd(24,24) - rrt(010) * density(83) 
  pd(24,83) = pd(24,83) - rrt(010) * density(24) 
  pd(25,24) = pd(25,24) + rrt(010) * density(83) 
  pd(25,83) = pd(25,83) + rrt(010) * density(24) 
  pd(02,24) = pd(02,24) + rrt(011) * density(83) 
  pd(02,83) = pd(02,83) + rrt(011) * density(24) 
  pd(24,24) = pd(24,24) - rrt(011) * density(83) 
  pd(24,83) = pd(24,83) - rrt(011) * density(24) 
  pd(26,24) = pd(26,24) + rrt(011) * density(83) 
  pd(26,83) = pd(26,83) + rrt(011) * density(24) 
  pd(01,24) = pd(01,24) + rrt(012) * density(83) 
  pd(01,83) = pd(01,83) + rrt(012) * density(24) 
  pd(02,24) = pd(02,24) + rrt(012) * density(83) 
  pd(02,83) = pd(02,83) + rrt(012) * density(24) 
  pd(24,24) = pd(24,24) - rrt(012) * density(83) 
  pd(24,83) = pd(24,83) - rrt(012) * density(24) 
  pd(29,24) = pd(29,24) + rrt(012) * density(83) 
  pd(29,83) = pd(29,83) + rrt(012) * density(24) 
  pd(02,24) = pd(02,24) + rrt(013) * density(83) * 2.0d0
  pd(02,83) = pd(02,83) + rrt(013) * density(24) * 2.0d0
  pd(24,24) = pd(24,24) - rrt(013) * density(83) 
  pd(24,83) = pd(24,83) - rrt(013) * density(24) 
  pd(28,24) = pd(28,24) + rrt(013) * density(83) 
  pd(28,83) = pd(28,83) + rrt(013) * density(24) 
  pd(24,24) = pd(24,24) - rrt(014) * density(83) 
  pd(24,83) = pd(24,83) - rrt(014) * density(24) 
  pd(70,24) = pd(70,24) + rrt(014) * density(83) 
  pd(70,83) = pd(70,83) + rrt(014) * density(24) 
  pd(83,24) = pd(83,24) + rrt(014) * density(83) 
  pd(83,83) = pd(83,83) + rrt(014) * density(24) 
  pd(01,24) = pd(01,24) + rrt(015) * density(83) 
  pd(01,83) = pd(01,83) + rrt(015) * density(24) 
  pd(24,24) = pd(24,24) - rrt(015) * density(83) 
  pd(24,83) = pd(24,83) - rrt(015) * density(24) 
  pd(71,24) = pd(71,24) + rrt(015) * density(83) 
  pd(71,83) = pd(71,83) + rrt(015) * density(24) 
  pd(83,24) = pd(83,24) + rrt(015) * density(83) 
  pd(83,83) = pd(83,83) + rrt(015) * density(24) 
  pd(04,24) = pd(04,24) + rrt(016) * density(73) 
  pd(04,73) = pd(04,73) + rrt(016) * density(24) 
  pd(73,24) = pd(73,24) - rrt(016) * density(73) 
  pd(73,73) = pd(73,73) - rrt(016) * density(24) 
  pd(08,24) = pd(08,24) + rrt(017) * density(73) 
  pd(08,73) = pd(08,73) + rrt(017) * density(24) 
  pd(24,24) = pd(24,24) - rrt(017) * density(73) 
  pd(24,73) = pd(24,73) - rrt(017) * density(24) 
  pd(25,24) = pd(25,24) + rrt(017) * density(73) 
  pd(25,73) = pd(25,73) + rrt(017) * density(24) 
  pd(73,24) = pd(73,24) - rrt(017) * density(73) 
  pd(73,73) = pd(73,73) - rrt(017) * density(24) 
  pd(73,24) = pd(73,24) + rrt(018) * density(74) 
  pd(73,74) = pd(73,74) + rrt(018) * density(24) 
  pd(74,24) = pd(74,24) - rrt(018) * density(74) 
  pd(74,74) = pd(74,74) - rrt(018) * density(24) 
  pd(74,24) = pd(74,24) + rrt(019) * density(75) 
  pd(74,75) = pd(74,75) + rrt(019) * density(24) 
  pd(75,24) = pd(75,24) - rrt(019) * density(75) 
  pd(75,75) = pd(75,75) - rrt(019) * density(24) 
  pd(05,24) = pd(05,24) + rrt(020) * density(76) 
  pd(05,76) = pd(05,76) + rrt(020) * density(24) 
  pd(24,24) = pd(24,24) - rrt(020) * density(76) 
  pd(24,76) = pd(24,76) - rrt(020) * density(24) 
  pd(25,24) = pd(25,24) + rrt(020) * density(76) 
  pd(25,76) = pd(25,76) + rrt(020) * density(24) 
  pd(76,24) = pd(76,24) - rrt(020) * density(76) 
  pd(76,76) = pd(76,76) - rrt(020) * density(24) 
  pd(01,24) = pd(01,24) + rrt(021) * density(76) 
  pd(01,76) = pd(01,76) + rrt(021) * density(24) 
  pd(20,24) = pd(20,24) + rrt(021) * density(76) 
  pd(20,76) = pd(20,76) + rrt(021) * density(24) 
  pd(24,24) = pd(24,24) - rrt(021) * density(76) 
  pd(24,76) = pd(24,76) - rrt(021) * density(24) 
  pd(76,24) = pd(76,24) - rrt(021) * density(76) 
  pd(76,76) = pd(76,76) - rrt(021) * density(24) 
  pd(02,24) = pd(02,24) + rrt(022) * density(76) 
  pd(02,76) = pd(02,76) + rrt(022) * density(24) 
  pd(14,24) = pd(14,24) + rrt(022) * density(76) 
  pd(14,76) = pd(14,76) + rrt(022) * density(24) 
  pd(24,24) = pd(24,24) - rrt(022) * density(76) 
  pd(24,76) = pd(24,76) - rrt(022) * density(24) 
  pd(76,24) = pd(76,24) - rrt(022) * density(76) 
  pd(76,76) = pd(76,76) - rrt(022) * density(24) 
  pd(01,24) = pd(01,24) + rrt(023) * density(81) 
  pd(01,81) = pd(01,81) + rrt(023) * density(24) 
  pd(02,24) = pd(02,24) + rrt(023) * density(81) 
  pd(02,81) = pd(02,81) + rrt(023) * density(24) 
  pd(11,24) = pd(11,24) + rrt(023) * density(81) 
  pd(11,81) = pd(11,81) + rrt(023) * density(24) 
  pd(24,24) = pd(24,24) - rrt(023) * density(81) 
  pd(24,81) = pd(24,81) - rrt(023) * density(24) 
  pd(29,24) = pd(29,24) + rrt(023) * density(81) 
  pd(29,81) = pd(29,81) + rrt(023) * density(24) 
  pd(81,24) = pd(81,24) - rrt(023) * density(81) 
  pd(81,81) = pd(81,81) - rrt(023) * density(24) 
  pd(11,24) = pd(11,24) + rrt(024) * density(81) 
  pd(11,81) = pd(11,81) + rrt(024) * density(24) 
  pd(24,24) = pd(24,24) - rrt(024) * density(81) 
  pd(24,81) = pd(24,81) - rrt(024) * density(24) 
  pd(70,24) = pd(70,24) + rrt(024) * density(81) 
  pd(70,81) = pd(70,81) + rrt(024) * density(24) 
  pd(81,24) = pd(81,24) - rrt(024) * density(81) 
  pd(81,81) = pd(81,81) - rrt(024) * density(24) 
  pd(83,24) = pd(83,24) + rrt(024) * density(81) 
  pd(83,81) = pd(83,81) + rrt(024) * density(24) 
  pd(01,24) = pd(01,24) + rrt(025) * density(81) 
  pd(01,81) = pd(01,81) + rrt(025) * density(24) 
  pd(11,24) = pd(11,24) + rrt(025) * density(81) 
  pd(11,81) = pd(11,81) + rrt(025) * density(24) 
  pd(24,24) = pd(24,24) - rrt(025) * density(81) 
  pd(24,81) = pd(24,81) - rrt(025) * density(24) 
  pd(71,24) = pd(71,24) + rrt(025) * density(81) 
  pd(71,81) = pd(71,81) + rrt(025) * density(24) 
  pd(81,24) = pd(81,24) - rrt(025) * density(81) 
  pd(81,81) = pd(81,81) - rrt(025) * density(24) 
  pd(83,24) = pd(83,24) + rrt(025) * density(81) 
  pd(83,81) = pd(83,81) + rrt(025) * density(24) 
  pd(02,24) = pd(02,24) + rrt(026) * density(81) 
  pd(02,81) = pd(02,81) + rrt(026) * density(24) 
  pd(11,24) = pd(11,24) + rrt(026) * density(81) 
  pd(11,81) = pd(11,81) + rrt(026) * density(24) 
  pd(24,24) = pd(24,24) - rrt(026) * density(81) 
  pd(24,81) = pd(24,81) - rrt(026) * density(24) 
  pd(72,24) = pd(72,24) + rrt(026) * density(81) 
  pd(72,81) = pd(72,81) + rrt(026) * density(24) 
  pd(81,24) = pd(81,24) - rrt(026) * density(81) 
  pd(81,81) = pd(81,81) - rrt(026) * density(24) 
  pd(83,24) = pd(83,24) + rrt(026) * density(81) 
  pd(83,81) = pd(83,81) + rrt(026) * density(24) 
  pd(04,03) = pd(04,03) + rrt(027) * density(04) * density(73) 
  pd(04,04) = pd(04,04) + rrt(027) * density(03) * density(73) 
  pd(04,73) = pd(04,73) + rrt(027) * density(03) * density(04) 
  pd(73,03) = pd(73,03) - rrt(027) * density(04) * density(73) 
  pd(73,04) = pd(73,04) - rrt(027) * density(03) * density(73) 
  pd(73,73) = pd(73,73) - rrt(027) * density(03) * density(04) 
  pd(04,03) = pd(04,03) + rrt(028) * density(11) * density(73) 
  pd(04,11) = pd(04,11) + rrt(028) * density(03) * density(73) 
  pd(04,73) = pd(04,73) + rrt(028) * density(03) * density(11) 
  pd(73,03) = pd(73,03) - rrt(028) * density(11) * density(73) 
  pd(73,11) = pd(73,11) - rrt(028) * density(03) * density(73) 
  pd(73,73) = pd(73,73) - rrt(028) * density(03) * density(11) 
  pd(73,11) = pd(73,11) + rrt(029) * density(75) 
  pd(73,75) = pd(73,75) + rrt(029) * density(11) 
  pd(75,11) = pd(75,11) - rrt(029) * density(75) 
  pd(75,75) = pd(75,75) - rrt(029) * density(11) 
  pd(74,11) = pd(74,11) + rrt(030) * density(75) 
  pd(74,75) = pd(74,75) + rrt(030) * density(11) 
  pd(75,11) = pd(75,11) - rrt(030) * density(75) 
  pd(75,75) = pd(75,75) - rrt(030) * density(11) 
  pd(03,11) = pd(03,11) + rrt(031) * density(76) 
  pd(03,76) = pd(03,76) + rrt(031) * density(11) 
  pd(76,11) = pd(76,11) - rrt(031) * density(76) 
  pd(76,76) = pd(76,76) - rrt(031) * density(11) 
  pd(04,11) = pd(04,11) + rrt(032) * density(73) 
  pd(04,73) = pd(04,73) + rrt(032) * density(11) 
  pd(73,11) = pd(73,11) - rrt(032) * density(73) 
  pd(73,73) = pd(73,73) - rrt(032) * density(11) 
  pd(73,11) = pd(73,11) + rrt(033) * density(74) 
  pd(73,74) = pd(73,74) + rrt(033) * density(11) 
  pd(74,11) = pd(74,11) - rrt(033) * density(74) 
  pd(74,74) = pd(74,74) - rrt(033) * density(11) 
  pd(11,81) = pd(11,81) + rrt(034) * density(81) * 2.0d0
  pd(81,81) = pd(81,81) - rrt(034) * density(81) * 4.0d0
  pd(82,81) = pd(82,81) + rrt(034) * density(81) * 2.0d0
  pd(83,81) = pd(83,81) + rrt(034) * density(81) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(035) * density(81) 
  pd(04,81) = pd(04,81) - rrt(035) * density(04) 
  pd(11,04) = pd(11,04) + rrt(035) * density(81) 
  pd(11,81) = pd(11,81) + rrt(035) * density(04) 
  pd(78,04) = pd(78,04) + rrt(035) * density(81) 
  pd(78,81) = pd(78,81) + rrt(035) * density(04) 
  pd(81,04) = pd(81,04) - rrt(035) * density(81) 
  pd(81,81) = pd(81,81) - rrt(035) * density(04) 
  pd(83,04) = pd(83,04) + rrt(035) * density(81) 
  pd(83,81) = pd(83,81) + rrt(035) * density(04) 
  pd(03,80) = pd(03,80) + rrt(036) * density(81) 
  pd(03,81) = pd(03,81) + rrt(036) * density(80) 
  pd(11,80) = pd(11,80) + rrt(036) * density(81) 
  pd(11,81) = pd(11,81) + rrt(036) * density(80) 
  pd(78,80) = pd(78,80) + rrt(036) * density(81) 
  pd(78,81) = pd(78,81) + rrt(036) * density(80) 
  pd(80,80) = pd(80,80) - rrt(036) * density(81) 
  pd(80,81) = pd(80,81) - rrt(036) * density(80) 
  pd(81,80) = pd(81,80) - rrt(036) * density(81) 
  pd(81,81) = pd(81,81) - rrt(036) * density(80) 
  pd(83,80) = pd(83,80) + rrt(036) * density(81) 
  pd(83,81) = pd(83,81) + rrt(036) * density(80) 
  pd(11,73) = pd(11,73) + rrt(037) * density(81) 
  pd(11,81) = pd(11,81) + rrt(037) * density(73) 
  pd(73,73) = pd(73,73) - rrt(037) * density(81) 
  pd(73,81) = pd(73,81) - rrt(037) * density(73) 
  pd(78,73) = pd(78,73) + rrt(037) * density(81) 
  pd(78,81) = pd(78,81) + rrt(037) * density(73) 
  pd(81,73) = pd(81,73) - rrt(037) * density(81) 
  pd(81,81) = pd(81,81) - rrt(037) * density(73) 
  pd(83,73) = pd(83,73) + rrt(037) * density(81) 
  pd(83,81) = pd(83,81) + rrt(037) * density(73) 
  pd(03,03) = pd(03,03) - rrt(038) * density(81) 
  pd(03,81) = pd(03,81) - rrt(038) * density(03) 
  pd(11,03) = pd(11,03) + rrt(038) * density(81) 
  pd(11,81) = pd(11,81) + rrt(038) * density(03) 
  pd(79,03) = pd(79,03) + rrt(038) * density(81) 
  pd(79,81) = pd(79,81) + rrt(038) * density(03) 
  pd(81,03) = pd(81,03) - rrt(038) * density(81) 
  pd(81,81) = pd(81,81) - rrt(038) * density(03) 
  pd(83,03) = pd(83,03) + rrt(038) * density(81) 
  pd(83,81) = pd(83,81) + rrt(038) * density(03) 
  pd(11,76) = pd(11,76) + rrt(039) * density(81) 
  pd(11,81) = pd(11,81) + rrt(039) * density(76) 
  pd(76,76) = pd(76,76) - rrt(039) * density(81) 
  pd(76,81) = pd(76,81) - rrt(039) * density(76) 
  pd(79,76) = pd(79,76) + rrt(039) * density(81) 
  pd(79,81) = pd(79,81) + rrt(039) * density(76) 
  pd(81,76) = pd(81,76) - rrt(039) * density(81) 
  pd(81,81) = pd(81,81) - rrt(039) * density(76) 
  pd(83,76) = pd(83,76) + rrt(039) * density(81) 
  pd(83,81) = pd(83,81) + rrt(039) * density(76) 
  pd(11,77) = pd(11,77) + rrt(040) * density(81) 
  pd(11,81) = pd(11,81) + rrt(040) * density(77) 
  pd(77,77) = pd(77,77) - rrt(040) * density(81) 
  pd(77,81) = pd(77,81) - rrt(040) * density(77) 
  pd(79,77) = pd(79,77) + rrt(040) * density(81) 
  pd(79,81) = pd(79,81) + rrt(040) * density(77) 
  pd(81,77) = pd(81,77) - rrt(040) * density(81) 
  pd(81,81) = pd(81,81) - rrt(040) * density(77) 
  pd(83,77) = pd(83,77) + rrt(040) * density(81) 
  pd(83,81) = pd(83,81) + rrt(040) * density(77) 
  pd(81,82) = pd(81,82) + rrt(041) * density(83) 
  pd(81,83) = pd(81,83) + rrt(041) * density(82) 
  pd(82,82) = pd(82,82) - rrt(041) * density(83) 
  pd(82,83) = pd(82,83) - rrt(041) * density(82) 
  pd(83,82) = pd(83,82) - rrt(041) * density(83) 
  pd(83,83) = pd(83,83) - rrt(041) * density(82) 
  pd(81,82) = pd(81,82) + rrt(042) * density(83)**2 
  pd(81,83) = pd(81,83) + rrt(042) * density(82) * density(83) * 2.0d0
  pd(82,82) = pd(82,82) - rrt(042) * density(83)**2 
  pd(82,83) = pd(82,83) - rrt(042) * density(82) * density(83) * 2.0d0
  pd(83,82) = pd(83,82) - rrt(042) * density(83)**2 
  pd(83,83) = pd(83,83) - rrt(042) * density(82) * density(83) * 2.0d0
  pd(01,70) = pd(01,70) + rrt(043) * density(83) * 2.0d0
  pd(01,83) = pd(01,83) + rrt(043) * density(70) * 2.0d0
  pd(26,70) = pd(26,70) + rrt(043) * density(83) 
  pd(26,83) = pd(26,83) + rrt(043) * density(70) 
  pd(70,70) = pd(70,70) - rrt(043) * density(83) 
  pd(70,83) = pd(70,83) - rrt(043) * density(70) 
  pd(83,70) = pd(83,70) - rrt(043) * density(83) 
  pd(83,83) = pd(83,83) - rrt(043) * density(70) 
  pd(01,70) = pd(01,70) + rrt(044) * density(83) 
  pd(01,83) = pd(01,83) + rrt(044) * density(70) 
  pd(25,70) = pd(25,70) + rrt(044) * density(83) 
  pd(25,83) = pd(25,83) + rrt(044) * density(70) 
  pd(70,70) = pd(70,70) - rrt(044) * density(83) 
  pd(70,83) = pd(70,83) - rrt(044) * density(70) 
  pd(83,70) = pd(83,70) - rrt(044) * density(83) 
  pd(83,83) = pd(83,83) - rrt(044) * density(70) 
  pd(01,71) = pd(01,71) + rrt(045) * density(83) 
  pd(01,83) = pd(01,83) + rrt(045) * density(71) 
  pd(26,71) = pd(26,71) + rrt(045) * density(83) 
  pd(26,83) = pd(26,83) + rrt(045) * density(71) 
  pd(71,71) = pd(71,71) - rrt(045) * density(83) 
  pd(71,83) = pd(71,83) - rrt(045) * density(71) 
  pd(83,71) = pd(83,71) - rrt(045) * density(83) 
  pd(83,83) = pd(83,83) - rrt(045) * density(71) 
  pd(01,72) = pd(01,72) + rrt(046) * density(83) 
  pd(01,83) = pd(01,83) + rrt(046) * density(72) 
  pd(29,72) = pd(29,72) + rrt(046) * density(83) 
  pd(29,83) = pd(29,83) + rrt(046) * density(72) 
  pd(72,72) = pd(72,72) - rrt(046) * density(83) 
  pd(72,83) = pd(72,83) - rrt(046) * density(72) 
  pd(83,72) = pd(83,72) - rrt(046) * density(83) 
  pd(83,83) = pd(83,83) - rrt(046) * density(72) 
  pd(04,04) = pd(04,04) - rrt(047) * density(70) 
  pd(04,70) = pd(04,70) - rrt(047) * density(04) 
  pd(24,04) = pd(24,04) + rrt(047) * density(70) 
  pd(24,70) = pd(24,70) + rrt(047) * density(04) 
  pd(70,04) = pd(70,04) - rrt(047) * density(70) 
  pd(70,70) = pd(70,70) - rrt(047) * density(04) 
  pd(78,04) = pd(78,04) + rrt(047) * density(70) 
  pd(78,70) = pd(78,70) + rrt(047) * density(04) 
  pd(03,04) = pd(03,04) + rrt(048) * density(82) 
  pd(03,82) = pd(03,82) + rrt(048) * density(04) 
  pd(04,04) = pd(04,04) - rrt(048) * density(82) 
  pd(04,82) = pd(04,82) - rrt(048) * density(04) 
  pd(11,04) = pd(11,04) + rrt(048) * density(82) 
  pd(11,82) = pd(11,82) + rrt(048) * density(04) 
  pd(79,04) = pd(79,04) + rrt(048) * density(82) 
  pd(79,82) = pd(79,82) + rrt(048) * density(04) 
  pd(82,04) = pd(82,04) - rrt(048) * density(82) 
  pd(82,82) = pd(82,82) - rrt(048) * density(04) 
  pd(04,80) = pd(04,80) + rrt(049) * density(82) 
  pd(04,82) = pd(04,82) + rrt(049) * density(80) 
  pd(11,80) = pd(11,80) + rrt(049) * density(82) 
  pd(11,82) = pd(11,82) + rrt(049) * density(80) 
  pd(79,80) = pd(79,80) + rrt(049) * density(82) 
  pd(79,82) = pd(79,82) + rrt(049) * density(80) 
  pd(80,80) = pd(80,80) - rrt(049) * density(82) 
  pd(80,82) = pd(80,82) - rrt(049) * density(80) 
  pd(82,80) = pd(82,80) - rrt(049) * density(82) 
  pd(82,82) = pd(82,82) - rrt(049) * density(80) 
  pd(04,04) = pd(04,04) - rrt(050) * density(82) 
  pd(04,82) = pd(04,82) - rrt(050) * density(04) 
  pd(11,04) = pd(11,04) + rrt(050) * density(82) 
  pd(11,82) = pd(11,82) + rrt(050) * density(04) 
  pd(78,04) = pd(78,04) + rrt(050) * density(82) 
  pd(78,82) = pd(78,82) + rrt(050) * density(04) 
  pd(82,04) = pd(82,04) - rrt(050) * density(82) 
  pd(82,82) = pd(82,82) - rrt(050) * density(04) 
  pd(03,74) = pd(03,74) + rrt(051) * density(82) 
  pd(03,82) = pd(03,82) + rrt(051) * density(74) 
  pd(11,74) = pd(11,74) + rrt(051) * density(82) 
  pd(11,82) = pd(11,82) + rrt(051) * density(74) 
  pd(74,74) = pd(74,74) - rrt(051) * density(82) 
  pd(74,82) = pd(74,82) - rrt(051) * density(74) 
  pd(79,74) = pd(79,74) + rrt(051) * density(82) 
  pd(79,82) = pd(79,82) + rrt(051) * density(74) 
  pd(82,74) = pd(82,74) - rrt(051) * density(82) 
  pd(82,82) = pd(82,82) - rrt(051) * density(74) 
  pd(11,74) = pd(11,74) + rrt(052) * density(82) 
  pd(11,82) = pd(11,82) + rrt(052) * density(74) 
  pd(74,74) = pd(74,74) - rrt(052) * density(82) 
  pd(74,82) = pd(74,82) - rrt(052) * density(74) 
  pd(78,74) = pd(78,74) + rrt(052) * density(82) 
  pd(78,82) = pd(78,82) + rrt(052) * density(74) 
  pd(82,74) = pd(82,74) - rrt(052) * density(82) 
  pd(82,82) = pd(82,82) - rrt(052) * density(74) 
  pd(03,03) = pd(03,03) - rrt(053) * density(82) 
  pd(03,82) = pd(03,82) - rrt(053) * density(03) 
  pd(11,03) = pd(11,03) + rrt(053) * density(82) 
  pd(11,82) = pd(11,82) + rrt(053) * density(03) 
  pd(79,03) = pd(79,03) + rrt(053) * density(82) 
  pd(79,82) = pd(79,82) + rrt(053) * density(03) 
  pd(82,03) = pd(82,03) - rrt(053) * density(82) 
  pd(82,82) = pd(82,82) - rrt(053) * density(03) 
  pd(11,76) = pd(11,76) + rrt(054) * density(82) 
  pd(11,82) = pd(11,82) + rrt(054) * density(76) 
  pd(76,76) = pd(76,76) - rrt(054) * density(82) 
  pd(76,82) = pd(76,82) - rrt(054) * density(76) 
  pd(79,76) = pd(79,76) + rrt(054) * density(82) 
  pd(79,82) = pd(79,82) + rrt(054) * density(76) 
  pd(82,76) = pd(82,76) - rrt(054) * density(82) 
  pd(82,82) = pd(82,82) - rrt(054) * density(76) 
  pd(11,77) = pd(11,77) + rrt(055) * density(82) 
  pd(11,82) = pd(11,82) + rrt(055) * density(77) 
  pd(77,77) = pd(77,77) - rrt(055) * density(82) 
  pd(77,82) = pd(77,82) - rrt(055) * density(77) 
  pd(79,77) = pd(79,77) + rrt(055) * density(82) 
  pd(79,82) = pd(79,82) + rrt(055) * density(77) 
  pd(82,77) = pd(82,77) - rrt(055) * density(82) 
  pd(82,82) = pd(82,82) - rrt(055) * density(77) 
  pd(03,03) = pd(03,03) - rrt(056) * density(03) * density(11) * 4.0d0
  pd(03,11) = pd(03,11) - rrt(056) * density(03)**2 * 2.0d0
  pd(04,03) = pd(04,03) + rrt(056) * density(03) * density(11) * 2.0d0
  pd(04,11) = pd(04,11) + rrt(056) * density(03)**2 
  pd(03,03) = pd(03,03) - rrt(057) * density(03) * density(11) * 4.0d0
  pd(03,11) = pd(03,11) - rrt(057) * density(03)**2 * 2.0d0
  pd(73,03) = pd(73,03) + rrt(057) * density(03) * density(11) * 2.0d0
  pd(73,11) = pd(73,11) + rrt(057) * density(03)**2 
  pd(03,03) = pd(03,03) - rrt(058) * density(04) * density(11) 
  pd(03,04) = pd(03,04) - rrt(058) * density(03) * density(11) 
  pd(03,11) = pd(03,11) - rrt(058) * density(03) * density(04) 
  pd(04,03) = pd(04,03) - rrt(058) * density(04) * density(11) 
  pd(04,04) = pd(04,04) - rrt(058) * density(03) * density(11) 
  pd(04,11) = pd(04,11) - rrt(058) * density(03) * density(04) 
  pd(80,03) = pd(80,03) + rrt(058) * density(04) * density(11) 
  pd(80,04) = pd(80,04) + rrt(058) * density(03) * density(11) 
  pd(80,11) = pd(80,11) + rrt(058) * density(03) * density(04) 
  pd(03,78) = pd(03,78) + rrt(059) * density(83) * 2.0d0
  pd(03,83) = pd(03,83) + rrt(059) * density(78) * 2.0d0
  pd(78,78) = pd(78,78) - rrt(059) * density(83) 
  pd(78,83) = pd(78,83) - rrt(059) * density(78) 
  pd(83,78) = pd(83,78) - rrt(059) * density(83) 
  pd(83,83) = pd(83,83) - rrt(059) * density(78) 
  pd(03,78) = pd(03,78) + rrt(060) * density(83) 
  pd(03,83) = pd(03,83) + rrt(060) * density(78) 
  pd(76,78) = pd(76,78) + rrt(060) * density(83) 
  pd(76,83) = pd(76,83) + rrt(060) * density(78) 
  pd(78,78) = pd(78,78) - rrt(060) * density(83) 
  pd(78,83) = pd(78,83) - rrt(060) * density(78) 
  pd(83,78) = pd(83,78) - rrt(060) * density(83) 
  pd(83,83) = pd(83,83) - rrt(060) * density(78) 
  pd(76,78) = pd(76,78) + rrt(061) * density(83) * 2.0d0
  pd(76,83) = pd(76,83) + rrt(061) * density(78) * 2.0d0
  pd(78,78) = pd(78,78) - rrt(061) * density(83) 
  pd(78,83) = pd(78,83) - rrt(061) * density(78) 
  pd(83,78) = pd(83,78) - rrt(061) * density(83) 
  pd(83,83) = pd(83,83) - rrt(061) * density(78) 
  pd(76,78) = pd(76,78) + rrt(062) * density(83) 
  pd(76,83) = pd(76,83) + rrt(062) * density(78) 
  pd(77,78) = pd(77,78) + rrt(062) * density(83) 
  pd(77,83) = pd(77,83) + rrt(062) * density(78) 
  pd(78,78) = pd(78,78) - rrt(062) * density(83) 
  pd(78,83) = pd(78,83) - rrt(062) * density(78) 
  pd(83,78) = pd(83,78) - rrt(062) * density(83) 
  pd(83,83) = pd(83,83) - rrt(062) * density(78) 
  pd(03,79) = pd(03,79) + rrt(063) * density(83)**2 
  pd(03,83) = pd(03,83) + rrt(063) * density(79) * density(83) * 2.0d0
  pd(79,79) = pd(79,79) - rrt(063) * density(83)**2 
  pd(79,83) = pd(79,83) - rrt(063) * density(79) * density(83) * 2.0d0
  pd(83,79) = pd(83,79) - rrt(063) * density(83)**2 
  pd(83,83) = pd(83,83) - rrt(063) * density(79) * density(83) * 2.0d0
  pd(03,79) = pd(03,79) + rrt(064) * density(83) 
  pd(03,83) = pd(03,83) + rrt(064) * density(79) 
  pd(79,79) = pd(79,79) - rrt(064) * density(83) 
  pd(79,83) = pd(79,83) - rrt(064) * density(79) 
  pd(83,79) = pd(83,79) - rrt(064) * density(83) 
  pd(83,83) = pd(83,83) - rrt(064) * density(79) 
  pd(04,73) = pd(04,73) + rrt(065) 
  pd(73,73) = pd(73,73) - rrt(065) 
  pd(73,74) = pd(73,74) + rrt(066) 
  pd(74,74) = pd(74,74) - rrt(066) 
  pd(04,74) = pd(04,74) + rrt(067) 
  pd(74,74) = pd(74,74) - rrt(067) 
  pd(04,75) = pd(04,75) + rrt(068) 
  pd(75,75) = pd(75,75) - rrt(068) 
  pd(04,03) = pd(04,03) + rrt(069) * density(73) 
  pd(04,73) = pd(04,73) + rrt(069) * density(03) 
  pd(73,03) = pd(73,03) - rrt(069) * density(73) 
  pd(73,73) = pd(73,73) - rrt(069) * density(03) 
  pd(04,04) = pd(04,04) + rrt(070) * density(73) 
  pd(04,73) = pd(04,73) + rrt(070) * density(04) 
  pd(73,04) = pd(73,04) - rrt(070) * density(73) 
  pd(73,73) = pd(73,73) - rrt(070) * density(04) 
  pd(04,73) = pd(04,73) + rrt(071) * density(80) * 2.0d0
  pd(04,80) = pd(04,80) + rrt(071) * density(73) * 2.0d0
  pd(73,73) = pd(73,73) - rrt(071) * density(80) 
  pd(73,80) = pd(73,80) - rrt(071) * density(73) 
  pd(76,73) = pd(76,73) + rrt(071) * density(80) 
  pd(76,80) = pd(76,80) + rrt(071) * density(73) 
  pd(80,73) = pd(80,73) - rrt(071) * density(80) 
  pd(80,80) = pd(80,80) - rrt(071) * density(73) 
  pd(04,73) = pd(04,73) + rrt(072) * density(73) * 2.0d0
  pd(73,73) = pd(73,73) - rrt(072) * density(73) * 4.0d0
  pd(74,73) = pd(74,73) + rrt(072) * density(73) * 2.0d0
  pd(03,03) = pd(03,03) - rrt(073) * density(80) 
  pd(03,80) = pd(03,80) - rrt(073) * density(03) 
  pd(04,03) = pd(04,03) + rrt(073) * density(80) 
  pd(04,80) = pd(04,80) + rrt(073) * density(03) 
  pd(73,03) = pd(73,03) + rrt(073) * density(80) 
  pd(73,80) = pd(73,80) + rrt(073) * density(03) 
  pd(80,03) = pd(80,03) - rrt(073) * density(80) 
  pd(80,80) = pd(80,80) - rrt(073) * density(03) 
  pd(73,03) = pd(73,03) + rrt(074) * density(74) 
  pd(73,74) = pd(73,74) + rrt(074) * density(03) 
  pd(74,03) = pd(74,03) - rrt(074) * density(74) 
  pd(74,74) = pd(74,74) - rrt(074) * density(03) 
  pd(03,03) = pd(03,03) - rrt(075) * density(74) 
  pd(03,74) = pd(03,74) - rrt(075) * density(03) 
  pd(04,03) = pd(04,03) + rrt(075) * density(74) 
  pd(04,74) = pd(04,74) + rrt(075) * density(03) 
  pd(74,03) = pd(74,03) - rrt(075) * density(74) 
  pd(74,74) = pd(74,74) - rrt(075) * density(03) 
  pd(76,03) = pd(76,03) + rrt(075) * density(74) 
  pd(76,74) = pd(76,74) + rrt(075) * density(03) 
  pd(73,04) = pd(73,04) + rrt(076) * density(74) 
  pd(73,74) = pd(73,74) + rrt(076) * density(04) 
  pd(74,04) = pd(74,04) - rrt(076) * density(74) 
  pd(74,74) = pd(74,74) - rrt(076) * density(04) 
  pd(03,74) = pd(03,74) + rrt(077) * density(80) 
  pd(03,80) = pd(03,80) + rrt(077) * density(74) 
  pd(04,74) = pd(04,74) + rrt(077) * density(80) * 2.0d0
  pd(04,80) = pd(04,80) + rrt(077) * density(74) * 2.0d0
  pd(74,74) = pd(74,74) - rrt(077) * density(80) 
  pd(74,80) = pd(74,80) - rrt(077) * density(74) 
  pd(80,74) = pd(80,74) - rrt(077) * density(80) 
  pd(80,80) = pd(80,80) - rrt(077) * density(74) 
  pd(03,03) = pd(03,03) - rrt(078) * density(75) 
  pd(03,75) = pd(03,75) - rrt(078) * density(03) 
  pd(04,03) = pd(04,03) + rrt(078) * density(75) 
  pd(04,75) = pd(04,75) + rrt(078) * density(03) 
  pd(75,03) = pd(75,03) - rrt(078) * density(75) 
  pd(75,75) = pd(75,75) - rrt(078) * density(03) 
  pd(77,03) = pd(77,03) + rrt(078) * density(75) 
  pd(77,75) = pd(77,75) + rrt(078) * density(03) 
  pd(04,04) = pd(04,04) - rrt(079) * density(75) 
  pd(04,75) = pd(04,75) - rrt(079) * density(04) 
  pd(74,04) = pd(74,04) + rrt(079) * density(75) * 2.0d0
  pd(74,75) = pd(74,75) + rrt(079) * density(04) * 2.0d0
  pd(75,04) = pd(75,04) - rrt(079) * density(75) 
  pd(75,75) = pd(75,75) - rrt(079) * density(04) 
  pd(03,03) = pd(03,03) + rrt(080) * density(76) 
  pd(03,76) = pd(03,76) + rrt(080) * density(03) 
  pd(76,03) = pd(76,03) - rrt(080) * density(76) 
  pd(76,76) = pd(76,76) - rrt(080) * density(03) 
  pd(03,04) = pd(03,04) + rrt(081) * density(76) 
  pd(03,76) = pd(03,76) + rrt(081) * density(04) 
  pd(76,04) = pd(76,04) - rrt(081) * density(76) 
  pd(76,76) = pd(76,76) - rrt(081) * density(04) 
  pd(03,04) = pd(03,04) + rrt(082) * density(76) 
  pd(03,76) = pd(03,76) + rrt(082) * density(04) 
  pd(04,04) = pd(04,04) - rrt(082) * density(76) 
  pd(04,76) = pd(04,76) - rrt(082) * density(04) 
  pd(73,04) = pd(73,04) + rrt(082) * density(76) 
  pd(73,76) = pd(73,76) + rrt(082) * density(04) 
  pd(76,04) = pd(76,04) - rrt(082) * density(76) 
  pd(76,76) = pd(76,76) - rrt(082) * density(04) 
  pd(03,04) = pd(03,04) + rrt(083) * density(76) 
  pd(03,76) = pd(03,76) + rrt(083) * density(04) 
  pd(04,04) = pd(04,04) - rrt(083) * density(76) 
  pd(04,76) = pd(04,76) - rrt(083) * density(04) 
  pd(74,04) = pd(74,04) + rrt(083) * density(76) 
  pd(74,76) = pd(74,76) + rrt(083) * density(04) 
  pd(76,04) = pd(76,04) - rrt(083) * density(76) 
  pd(76,76) = pd(76,76) - rrt(083) * density(04) 
  pd(03,76) = pd(03,76) + rrt(084) * density(80) * 2.0d0
  pd(03,80) = pd(03,80) + rrt(084) * density(76) * 2.0d0
  pd(04,76) = pd(04,76) + rrt(084) * density(80) 
  pd(04,80) = pd(04,80) + rrt(084) * density(76) 
  pd(76,76) = pd(76,76) - rrt(084) * density(80) 
  pd(76,80) = pd(76,80) - rrt(084) * density(76) 
  pd(80,76) = pd(80,76) - rrt(084) * density(80) 
  pd(80,80) = pd(80,80) - rrt(084) * density(76) 
  pd(04,76) = pd(04,76) + rrt(085) * density(80) * 2.0d0
  pd(04,80) = pd(04,80) + rrt(085) * density(76) * 2.0d0
  pd(76,76) = pd(76,76) - rrt(085) * density(80) 
  pd(76,80) = pd(76,80) - rrt(085) * density(76) 
  pd(80,76) = pd(80,76) - rrt(085) * density(80) 
  pd(80,80) = pd(80,80) - rrt(085) * density(76) 
  pd(76,03) = pd(76,03) + rrt(086) * density(77) 
  pd(76,77) = pd(76,77) + rrt(086) * density(03) 
  pd(77,03) = pd(77,03) - rrt(086) * density(77) 
  pd(77,77) = pd(77,77) - rrt(086) * density(03) 
  pd(76,04) = pd(76,04) + rrt(087) * density(77) 
  pd(76,77) = pd(76,77) + rrt(087) * density(04) 
  pd(77,04) = pd(77,04) - rrt(087) * density(77) 
  pd(77,77) = pd(77,77) - rrt(087) * density(04) 
  pd(03,04) = pd(03,04) + rrt(088) * density(77) * 3.0d0
  pd(03,77) = pd(03,77) + rrt(088) * density(04) * 3.0d0
  pd(04,04) = pd(04,04) - rrt(088) * density(77) 
  pd(04,77) = pd(04,77) - rrt(088) * density(04) 
  pd(77,04) = pd(77,04) - rrt(088) * density(77) 
  pd(77,77) = pd(77,77) - rrt(088) * density(04) 
  pd(03,73) = pd(03,73) + rrt(089) * density(77) 
  pd(03,77) = pd(03,77) + rrt(089) * density(73) 
  pd(73,73) = pd(73,73) - rrt(089) * density(77) 
  pd(73,77) = pd(73,77) - rrt(089) * density(73) 
  pd(75,73) = pd(75,73) + rrt(089) * density(77) 
  pd(75,77) = pd(75,77) + rrt(089) * density(73) 
  pd(77,73) = pd(77,73) - rrt(089) * density(77) 
  pd(77,77) = pd(77,77) - rrt(089) * density(73) 
  pd(73,73) = pd(73,73) - rrt(090) * density(77) 
  pd(73,77) = pd(73,77) - rrt(090) * density(73) 
  pd(74,73) = pd(74,73) + rrt(090) * density(77) 
  pd(74,77) = pd(74,77) + rrt(090) * density(73) 
  pd(76,73) = pd(76,73) + rrt(090) * density(77) 
  pd(76,77) = pd(76,77) + rrt(090) * density(73) 
  pd(77,73) = pd(77,73) - rrt(090) * density(77) 
  pd(77,77) = pd(77,77) - rrt(090) * density(73) 
  pd(03,73) = pd(03,73) + rrt(091) * density(77) * 3.0d0
  pd(03,77) = pd(03,77) + rrt(091) * density(73) * 3.0d0
  pd(73,73) = pd(73,73) - rrt(091) * density(77) 
  pd(73,77) = pd(73,77) - rrt(091) * density(73) 
  pd(77,73) = pd(77,73) - rrt(091) * density(77) 
  pd(77,77) = pd(77,77) - rrt(091) * density(73) 
  pd(04,77) = pd(04,77) + rrt(092) * density(80) * 2.0d0
  pd(04,80) = pd(04,80) + rrt(092) * density(77) * 2.0d0
  pd(77,77) = pd(77,77) - rrt(092) * density(80) 
  pd(77,80) = pd(77,80) - rrt(092) * density(77) 
  pd(80,77) = pd(80,77) - rrt(092) * density(80) 
  pd(80,80) = pd(80,80) - rrt(092) * density(77) 
  pd(03,77) = pd(03,77) + rrt(093) * density(80) 
  pd(03,80) = pd(03,80) + rrt(093) * density(77) 
  pd(04,77) = pd(04,77) + rrt(093) * density(80) 
  pd(04,80) = pd(04,80) + rrt(093) * density(77) 
  pd(76,77) = pd(76,77) + rrt(093) * density(80) 
  pd(76,80) = pd(76,80) + rrt(093) * density(77) 
  pd(77,77) = pd(77,77) - rrt(093) * density(80) 
  pd(77,80) = pd(77,80) - rrt(093) * density(77) 
  pd(80,77) = pd(80,77) - rrt(093) * density(80) 
  pd(80,80) = pd(80,80) - rrt(093) * density(77) 
  pd(03,04) = pd(03,04) + rrt(094) * density(04) * 4.0d0
  pd(04,04) = pd(04,04) - rrt(094) * density(04) * 2.0d0
  pd(03,03) = pd(03,03) + rrt(095) * density(04) * 2.0d0
  pd(03,04) = pd(03,04) + rrt(095) * density(03) * 2.0d0
  pd(04,03) = pd(04,03) - rrt(095) * density(04) 
  pd(04,04) = pd(04,04) - rrt(095) * density(03) 
  pd(03,04) = pd(03,04) + rrt(096) * density(80) 
  pd(03,80) = pd(03,80) + rrt(096) * density(04) 
  pd(04,04) = pd(04,04) + rrt(096) * density(80) 
  pd(04,80) = pd(04,80) + rrt(096) * density(04) 
  pd(80,04) = pd(80,04) - rrt(096) * density(80) 
  pd(80,80) = pd(80,80) - rrt(096) * density(04) 
  pd(03,03) = pd(03,03) + rrt(097) * density(80) 
  pd(03,80) = pd(03,80) + rrt(097) * density(03) 
  pd(04,03) = pd(04,03) + rrt(097) * density(80) 
  pd(04,80) = pd(04,80) + rrt(097) * density(03) 
  pd(80,03) = pd(80,03) - rrt(097) * density(80) 
  pd(80,80) = pd(80,80) - rrt(097) * density(03) 
  pd(03,03) = pd(03,03) - rrt(098) * density(03) * density(04) * 4.0d0
  pd(03,04) = pd(03,04) - rrt(098) * density(03)**2 * 2.0d0
  pd(04,03) = pd(04,03) + rrt(098) * density(03) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) + rrt(098) * density(03)**2 
  pd(03,03) = pd(03,03) - rrt(099) * density(03)**2 * 6.0d0
  pd(04,03) = pd(04,03) + rrt(099) * density(03)**2 * 3.0d0
  pd(03,03) = pd(03,03) - rrt(100) * density(04)**2 
  pd(03,04) = pd(03,04) - rrt(100) * density(03) * density(04) * 2.0d0
  pd(04,03) = pd(04,03) - rrt(100) * density(04)**2 
  pd(04,04) = pd(04,04) - rrt(100) * density(03) * density(04) * 2.0d0
  pd(80,03) = pd(80,03) + rrt(100) * density(04)**2 
  pd(80,04) = pd(80,04) + rrt(100) * density(03) * density(04) * 2.0d0
  pd(03,03) = pd(03,03) - rrt(101) * density(03) * density(04) * 2.0d0
  pd(03,04) = pd(03,04) - rrt(101) * density(03)**2 
  pd(04,03) = pd(04,03) - rrt(101) * density(03) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(101) * density(03)**2 
  pd(80,03) = pd(80,03) + rrt(101) * density(03) * density(04) * 2.0d0
  pd(80,04) = pd(80,04) + rrt(101) * density(03)**2 
  pd(03,04) = pd(03,04) + rrt(102) * density(79) 
  pd(03,79) = pd(03,79) + rrt(102) * density(04) 
  pd(04,04) = pd(04,04) - rrt(102) * density(79) 
  pd(04,79) = pd(04,79) - rrt(102) * density(04) 
  pd(78,04) = pd(78,04) + rrt(102) * density(79) 
  pd(78,79) = pd(78,79) + rrt(102) * density(04) 
  pd(79,04) = pd(79,04) - rrt(102) * density(79) 
  pd(79,79) = pd(79,79) - rrt(102) * density(04) 
  pd(04,79) = pd(04,79) + rrt(103) * density(80) 
  pd(04,80) = pd(04,80) + rrt(103) * density(79) 
  pd(78,79) = pd(78,79) + rrt(103) * density(80) 
  pd(78,80) = pd(78,80) + rrt(103) * density(79) 
  pd(79,79) = pd(79,79) - rrt(103) * density(80) 
  pd(79,80) = pd(79,80) - rrt(103) * density(79) 
  pd(80,79) = pd(80,79) - rrt(103) * density(80) 
  pd(80,80) = pd(80,80) - rrt(103) * density(79) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(84,8) = eV_to_K * ZDPlasKin_cfg(11)
    pd(84,:) = pd(84,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: EN
  double precision :: Te
  double precision :: ANY_NEUTRAL
  DOUBLE PRECISION :: DTION, TIONN,TEFFN ! K
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  call ZDPlasKin_get_density_total(ALL_NEUTRAL=ANY_NEUTRAL)
  DTION = 2.0D0 / ( 3.0D0 * 1.3807D-16 ) * 1.6605D-24 * ( 1.0D-17 * EN )**2
  TIONN = TGAS + DTION * 14.0D0 * 8.0D19**2
  TEFFN = ( TIONN + 0.5D0 * TGAS ) / ( 1.0D0 + 0.5D0 )
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = bolsig_rates(bolsig_pointer(14))
  rrt(015) = bolsig_rates(bolsig_pointer(15))
  rrt(016) = 1.4D-18
  rrt(017) = 6.14D-12*EXP(-1.79D4/TGAS)
  rrt(018) = 1.08D-13
  rrt(019) = 1.08D-13
  rrt(020) = 1.13D-10
  rrt(021) = 3.00D-11
  rrt(022) = 7.50D-12
  rrt(023) = 5.6D-13
  rrt(024) = 7.9D-12
  rrt(025) = 8.3D-12
  rrt(026) = 0.71D-12
  rrt(027) = 1.0D-32
  rrt(028) = 6.3D-33
  rrt(029) = 1.0D-13
  rrt(030) = 9.0D-13
  rrt(031) = 1.0D-13
  rrt(032) = 8.0D-21*(TGAS/300)**0.5
  rrt(033) = 1.0D-17*(TGAS/300)**0.5
  rrt(034) = 1.0D-9*(TGAS/300)**0.5
  rrt(035) = 2.54D-10*(TGAS/300)**0.5
  rrt(036) = rrt(35)
  rrt(037) = rrt(35)
  rrt(038) = rrt(35)
  rrt(039) = rrt(35)
  rrt(040) = rrt(35)
  rrt(041) = 6.76D-13*TE**(-0.5)
  rrt(042) = 5.12D-27*TE**(-4.5)
  rrt(043) = 1.7D-7*(300.0D0/TE)**0.5
  rrt(044) = rrt(43)
  rrt(045) = 3.5D-7*(300.0D0/TE)**0.5
  rrt(046) = 2.5D-7*(300.0D0/TE)**0.5
  rrt(047) = 5.0D-10
  rrt(048) = 1.07D-9*(TGAS/300)**0.5
  rrt(049) = rrt(48)
  rrt(050) = 3.30D-11*(TGAS/300)**0.5
  rrt(051) = rrt(48)
  rrt(052) = rrt(50)
  rrt(053) = 5.00D-11*(TGAS/300)**0.5
  rrt(054) = rrt(53)
  rrt(055) = rrt(53)
  rrt(056) = 1.00D-33
  rrt(057) = 9.88D-35
  rrt(058) = 3.40D-34*(TGAS/300)**(-1.2)
  rrt(059) = 1.95D-7*(300.0D0/TE)**0.7*0.262D0
  rrt(060) = 1.95D-7*(300.0D0/TE)**0.7*0.435D0
  rrt(061) = 1.95D-7*(300.0D0/TE)**0.7*0.2575D0
  rrt(062) = 1.95D-7*(300.0D0/TE)**0.7*0.0405D0
  rrt(063) = 7.0D-20*(300.0D0/TE)**4.5
  rrt(064) = 6.0D-27*(300.0D0/TE)**1.5*ANY_NEUTRAL
  rrt(065) = 2.6D-4
  rrt(066) = 1.5D-3
  rrt(067) = 8.5D-2
  rrt(068) = 11.0D0
  rrt(069) = 7.0D-16
  rrt(070) = 3.8D-18*EXP(-205.0D0/TGAS)
  rrt(071) = 5.2D-11*EXP(-2840.0D0/TGAS)
  rrt(072) = 7.0D-28*TGAS**3.8*EXP(700.0D0/TGAS)
  rrt(073) = 1.0D-11*EXP(-2300.0D0/TGAS)
  rrt(074) = 8.1D-14
  rrt(075) = 3.4D-11*(300.0D0/TGAS)**0.1*EXP(-4200.0D0/TGAS)
  rrt(076) = 4.3D-22*TGAS**2.4*EXP(-281.0D0/TGAS)
  rrt(077) = 2.2D-11
  rrt(078) = 9.0D-12
  rrt(079) = 3.0D-13
  rrt(080) = 8.0D-12
  rrt(081) = 6.4D-12*EXP(67.0D0/TGAS)
  rrt(082) = 1.0D-12
  rrt(083) = 2.6D-11*EXP(67.0D0/TGAS)
  rrt(084) = 1.2D-10
  rrt(085) = 1.2D-10
  rrt(086) = 5.0D-11*EXP(-300.0D0/TGAS)
  rrt(087) = 1.3D-12*EXP(-850.0D0/TGAS)
  rrt(088) = 3.0D-12*EXP(-850.0D0/TGAS)
  rrt(089) = 1.1D-10
  rrt(090) = 2.9D-11
  rrt(091) = 3.2D-11
  rrt(092) = 2.9D-10
  rrt(093) = 2.9D-10
  rrt(094) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*5.9D0
  rrt(095) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*21.D0
  rrt(096) = 6.6D-10*EXP(-11600.0D0/TGAS)*0.38D0
  rrt(097) = 6.6D-10*EXP(-11600.0D0/TGAS)*6.3D0*EXP(170.0D0/TGAS)
  rrt(098) = 4.0D-33*(300.0D0/TGAS)**0.41*1.0D0
  rrt(099) = 4.0D-33*(300.0D0/TGAS)**0.41*3.6D0
  rrt(100) = 7.6D-34*(300.0D0/TGAS)**1.9
  rrt(101) = MIN(3.9D-33*(300.0D0/TGAS)**1.9,1.1D-34*EXP(1060.0D0/TGAS))
  rrt(102) = 2.0D-11*(300.0D0/TEFFN)**0.5
  rrt(103) = 1.0D-10
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
