!
! ZDPLASKIN version 2.0
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
! Wed Jan 24 12:14:56 2018
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E N O H 
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
  integer, parameter :: species_max = 79, species_electrons = 79, species_length = 9, reactions_max = 792, reactions_length = 51
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
  integer, parameter, private               :: bolsig_species_max = 20, bolsig_species_length = 7, bolsig_rates_max = 106, &
                                               bolsig_addsect_max = 10 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0, &
                                               bolsig_addsect(2,bolsig_addsect_max) 
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
    subroutine ZDPlasKin_bolsig_GetNumCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetNumCollisions
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
    subroutine ZDPlasKin_bolsig_GetEEDF(a, i)
      integer, intent(out) :: i
      double precision, dimension(:,:), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,&
    1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1, 0, 0, 1, 1, 1,-1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1/
  data species_name(1:species_max) &
  /"N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ","N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","N2(V9)   ",&
   "N2(V10)  ","N2(V11)  ","N2(V12)  ","N2(V13)  ","N2(V14)  ","N2(V15)  ","N2(A3)   ","N2(B3)   ","N2(A`1)  ","N2(C3)   ",&
   "N2(B1)   ","N2(B'1)  ","N        ","N(2D)    ","N(2P)    ","N^+      ","N2^+     ","N3^+     ","N4^+     ","O2       ",&
   "O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","O2(A1)   ","O2(B1)   ","O2(4.5EV)","O        ","O(1D)    ","O(1S)    ",&
   "O3       ","O^+      ","O2^+     ","O3^+     ","O4^+     ","O^-      ","O2^-     ","O3^-     ","O4^-     ","NO       ",&
   "N2O      ","NO2      ","NO3      ","N2O5     ","NO^+     ","N2O^+    ","NO2^+    ","NO^-     ","N2O^-    ","NO2^-    ",&
   "NO3^-    ","H2       ","H        ","H^+      ","H2^+     ","H3^+     ","H^-      ","H2O      ","H3O^+    ","NH3      ",&
   "NH2      ","NH       ","NH3^+    ","NH2^+    ","NH^+     ","N2H^+    ","NH4^+    ","O2^+N2   ","E        "/
  data reaction_sign(1:36) &
  /"bolsig:N2->N2(V1)                                  ","bolsig:N2->N2(V2)                                  ",&
   "bolsig:N2->N2(V3)                                  ","bolsig:N2->N2(V4)                                  ",&
   "bolsig:N2->N2(V5)                                  ","bolsig:N2->N2(V6)                                  ",&
   "bolsig:N2->N2(V7)                                  ","bolsig:N2->N2(V8)                                  ",&
   "bolsig:N2->N2(V9)                                  ","bolsig:N2->N2(V10)                                 ",&
   "bolsig:N2->N2(V11)                                 ","bolsig:N2->N2(V12)                                 ",&
   "bolsig:N2->N2(V13)                                 ","bolsig:N2->N2(V14)                                 ",&
   "bolsig:N2->N2(V15)                                 ","bolsig:N2(V1)->N2                                  ",&
   "bolsig:N2(V1)->N2(V2)                              ","bolsig:N2(V1)->N2(V3)                              ",&
   "bolsig:N2(V1)->N2(V4)                              ","bolsig:N2(V2)->N2                                  ",&
   "bolsig:N2(V2)->N2(V1)                              ","bolsig:N2(V2)->N2(V3)                              ",&
   "bolsig:N2(V2)->N2(V4)                              ","bolsig:N2(V2)->N2(V5)                              ",&
   "bolsig:N2(V3)->N2                                  ","bolsig:N2(V3)->N2(V1)                              ",&
   "bolsig:N2(V3)->N2(V2)                              ","bolsig:N2(V3)->N2(V4)                              ",&
   "bolsig:N2(V3)->N2(V5)                              ","bolsig:N2(V3)->N2(V6)                              ",&
   "bolsig:N2(V4)->N2                                  ","bolsig:N2(V4)->N2(V1)                              ",&
   "bolsig:N2(V4)->N2(V2)                              ","bolsig:N2(V4)->N2(V3)                              ",&
   "bolsig:N2(V5)->N2                                  ","bolsig:N2(V5)->N2(V2)                              "/
  data reaction_sign(37:72) &
  /"bolsig:N2(V5)->N2(V3)                              ","bolsig:N2(V6)->N2                                  ",&
   "bolsig:N2(V6)->N2(V3)                              ","bolsig:N2(V7)->N2                                  ",&
   "bolsig:N2(V8)->N2                                  ","bolsig:N2(V9)->N2                                  ",&
   "bolsig:N2(V10)->N2                                 ","bolsig:O2->O2(V1RES)                               ",&
   "bolsig:O2->O2(V1)                                  ","bolsig:O2->O2(V2RES)                               ",&
   "bolsig:O2->O2(V2)                                  ","bolsig:O2->O2(V3)                                  ",&
   "bolsig:O2->O2(V4)                                  ","N2(V1)+N2=>N2+N2+2.91D-1_EV                        ",&
   "N2(V2)+N2=>N2(V1)+N2+2.99D-1_EV                    ","N2(V3)+N2=>N2(V2)+N2+2.90D-1_EV                    ",&
   "N2(V4)+N2=>N2(V3)+N2+2.90D-1_EV                    ","N2(V5)+N2=>N2(V4)+N2+3.00D-1_EV                    ",&
   "N2(V6)+N2=>N2(V5)+N2+2.90D-1_EV                    ","N2(V7)+N2=>N2(V6)+N2+3.00D-1_EV                    ",&
   "N2(V8)+N2=>N2(V7)+N2+2.90D-1_EV                    ","N2+N2+2.91D-1_EV=>N2(V1)+N2                        ",&
   "N2(V1)+N2+2.99D-1_EV=>N2(V2)+N2                    ","N2(V2)+N2+2.90D-1_EV=>N2(V3)+N2                    ",&
   "N2(V3)+N2+2.90D-1_EV=>N2(V4)+N2                    ","N2(V4)+N2+3.00D-1_EV=>N2(V5)+N2                    ",&
   "N2(V5)+N2+2.90D-1_EV=>N2(V6)+N2                    ","N2(V6)+N2+3.00D-1_EV=>N2(V7)+N2                    ",&
   "N2(V7)+N2+2.90D-1_EV=>N2(V8)+N2                    ","N2(V1)+N=>N2+N+2.91D-1_EV                          ",&
   "N2(V2)+N=>N2(V1)+N+2.99D-1_EV                      ","N2(V3)+N=>N2(V2)+N+2.90D-1_EV                      ",&
   "N2(V4)+N=>N2(V3)+N+2.90D-1_EV                      ","N2(V5)+N=>N2(V4)+N+3.00D-1_EV                      ",&
   "N2(V6)+N=>N2(V5)+N+2.90D-1_EV                      ","N2(V7)+N=>N2(V6)+N+3.00D-1_EV                      "/
  data reaction_sign(73:108) &
  /"N2(V8)+N=>N2(V7)+N+2.90D-1_EV                      ","N2+N+2.91D-1_EV=>N2(V1)+N                          ",&
   "N2(V1)+N+2.99D-1_EV=>N2(V2)+N                      ","N2(V2)+N+2.90D-1_EV=>N2(V3)+N                      ",&
   "N2(V3)+N+2.90D-1_EV=>N2(V4)+N                      ","N2(V4)+N+3.00D-1_EV=>N2(V5)+N                      ",&
   "N2(V5)+N+2.90D-1_EV=>N2(V6)+N                      ","N2(V6)+N+3.00D-1_EV=>N2(V7)+N                      ",&
   "N2(V7)+N+2.90D-1_EV=>N2(V8)+N                      ","N2(V1)+O=>N2+O+2.91D-1_EV                          ",&
   "N2(V2)+O=>N2(V1)+O+2.99D-1_EV                      ","N2(V3)+O=>N2(V2)+O+2.90D-1_EV                      ",&
   "N2(V4)+O=>N2(V3)+O+2.90D-1_EV                      ","N2(V5)+O=>N2(V4)+O+3.00D-1_EV                      ",&
   "N2(V6)+O=>N2(V5)+O+2.90D-1_EV                      ","N2(V7)+O=>N2(V6)+O+3.00D-1_EV                      ",&
   "N2(V8)+O=>N2(V7)+O+2.90D-1_EV                      ","N2+O+2.91D-1_EV=>N2(V1)+O                          ",&
   "N2(V1)+O+2.99D-1_EV=>N2(V2)+O                      ","N2(V2)+O+2.90D-1_EV=>N2(V3)+O                      ",&
   "N2(V3)+O+2.90D-1_EV=>N2(V4)+O                      ","N2(V4)+O+3.00D-1_EV=>N2(V5)+O                      ",&
   "N2(V5)+O+2.90D-1_EV=>N2(V6)+O                      ","N2(V6)+O+3.00D-1_EV=>N2(V7)+O                      ",&
   "N2(V7)+O+2.90D-1_EV=>N2(V8)+O                      ","O2(V1)+O2=>O2+O2+1.90D-01_EV                       ",&
   "O2(V2)+O2=>O2(V1)+O2+1.90D-01_EV                   ","O2(V3)+O2=>O2(V2)+O2+1.90D-01_EV                   ",&
   "O2(V4)+O2=>O2(V3)+O2+1.80D-01_EV                   ","O2+O2+1.90D-01_EV=>O2(V1)+O2                       ",&
   "O2(V1)+O2+1.90D-01_EV=>O2(V2)+O2                   ","O2(V2)+O2+1.90D-01_EV=>O2(V3)+O2                   ",&
   "O2(V3)+O2+1.80D-01_EV=>O2(V4)+O2                   ","O2(V1)+O=>O2+O+1.90D-01_EV                         ",&
   "O2(V2)+O=>O2(V1)+O+1.90D-01_EV                     ","O2(V3)+O=>O2(V2)+O+1.90D-01_EV                     "/
  data reaction_sign(109:144) &
  /"O2(V4)+O=>O2(V3)+O+1.80D-01_EV                     ","O2+O+1.90D-01_EV=>O2(V1)+O                         ",&
   "O2(V1)+O+1.90D-01_EV=>O2(V2)+O                     ","O2(V2)+O+1.90D-01_EV=>O2(V3)+O                     ",&
   "O2(V3)+O+1.80D-01_EV=>O2(V4)+O                     ","bolsig:N2->N2(A3,V0-4)                             ",&
   "bolsig:N2->N2(A3,V5-9)                             ","bolsig:N2->N2(A3,V10-)                             ",&
   "bolsig:N2->N2(B3,V0-3)                             ","bolsig:N2->N2(B3,V4-16)                            ",&
   "bolsig:N2->N2(W3,V0-5)                             ","bolsig:N2->N2(W3,V6-10)                            ",&
   "bolsig:N2->N2(W3,V11-19)                           ","bolsig:N2->N2(B'3,V0-6)                            ",&
   "bolsig:N2->N2(B'3,V7-18)                           ","bolsig:N2->N2(A'1,V0-6)                            ",&
   "bolsig:N2->N2(A'1,V7-19)                           ","bolsig:N2->N2(A1,V0-3)                             ",&
   "bolsig:N2->N2(A1,V4-15)                            ","bolsig:N2->N2(W1,V0-5)                             ",&
   "bolsig:N2->N2(W1,V6-18)                            ","bolsig:N2->N2(C3,V0-4)                             ",&
   "bolsig:N2->N2(E3)                                  ","bolsig:N2->N2(A''1)                                ",&
   "bolsig:N2->N2(B1,V0-6)                             ","bolsig:N2->N2(B1,V7-14)                            ",&
   "bolsig:N2->N2(C'1)                                 ","bolsig:N2->N2(G3)                                  ",&
   "bolsig:N2->N2(C31)                                 ","bolsig:N2->N2(F3)                                  ",&
   "bolsig:N2->N2(B'1,V0-10)                           ","bolsig:N2->N2(B'1,V10-H)                           ",&
   "bolsig:N2->N2(O31)                                 ","bolsig:N2->N2(SUMSINGLETS)                         ",&
   "bolsig:O2->O2(A1)                                  ","bolsig:O2->O2(B1)                                  "/
  data reaction_sign(145:180) &
  /"bolsig:O2->O2(4.5EV)                               ","bolsig:O2->O2(6.0EV)                               ",&
   "bolsig:O2->O2(8.4EV)                               ","bolsig:O2->O2(9.97EV)                              ",&
   "bolsig:O2(A1)->O+O*(6EV)                           ","bolsig:O->O(1D)                                    ",&
   "bolsig:O->O(1S)                                    ","bolsig:NH3->NH3(E1)(5.72EV)                        ",&
   "bolsig:NH3->NH3(E2)(8.65EV)                        ","bolsig:N->N^+                                      ",&
   "bolsig:O->O^+                                      ","bolsig:N2->N2^+                                    ",&
   "bolsig:N2(A3)->N2^+                                ","bolsig:O2->O2^+                                    ",&
   "bolsig:O2(A1)->O2^+                                ","bolsig:NO->NO^+                                    ",&
   "bolsig:N2O->N2O^+                                  ","bolsig:O3->O3^+                                    ",&
   "bolsig:NH3->NH3^+                                  ","bolsig:NH3->NH2+H^+                                ",&
   "bolsig:NH3->NH+H2^+                                ","E+N2^+=>N+N                                        ",&
   "E+N2^+=>N+N(2D)                                    ","E+N2^+=>N+N(2P)                                    ",&
   "E+O2^+=>O+O+6.95D0_EV                              ","E+O2^+=>O+O(1D)+4.99D0_EV                          ",&
   "E+O2^+=>O+O(1S)                                    ","E+NO^+=>O+N                                        ",&
   "E+NO^+=>O+N(2D)                                    ","E+N3^+=>N2+N                                       ",&
   "E+N4^+=>N2+N2                                      ","E+N2O^+=>N2+O                                      ",&
   "E+NO2^+=>NO+O                                      ","E+O4^+=>O2+O2                                      ",&
   "E+O2^+N2=>O2+N2                                    ","E+N^++E=>N+E                                       "/
  data reaction_sign(181:216) &
  /"E+O^++E=>O+E                                       ","E+N^++ANY_NEUTRAL=>N+ANY_NEUTRAL                   ",&
   "E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL                   ","E+NH4^+=>H+H+NH2                                   ",&
   "bolsig:O2->O^-+O                                   ","bolsig:NO->N+O^-                                   ",&
   "bolsig:O3->O^-+O2                                  ","bolsig:O3->O2^-+O                                  ",&
   "E+NO2=>O^-+NO                                      ","E+O+O2=>O^-+O2                                     ",&
   "E+O+O2=>O2^-+O                                     ","E+O3+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL                 ",&
   "E+NO+ANY_NEUTRAL=>NO^-+ANY_NEUTRAL                 ","E+N2O+ANY_NEUTRAL=>N2O^-+ANY_NEUTRAL               ",&
   "E+O2+N2=>O2^-+N2                                   ","bolsig:NH3->NH2+H^-                                ",&
   "O^-+O=>O2+E                                        ","O^-+N=>NO+E                                        ",&
   "O^-+NO=>NO2+E                                      ","O^-+N2=>N2O+E                                      ",&
   "O^-+O2=>O3+E                                       ","O^-+O2(A1)=>O3+E                                   ",&
   "O^-+O2(B1)=>O+O2+E                                 ","O^-+N2(A3)=>O+N2+E                                 ",&
   "O^-+N2(B3)=>O+N2+E                                 ","O^-+O3=>O2+O2+E                                    ",&
   "O2^-+O=>O3+E                                       ","O2^-+N=>NO2+E                                      ",&
   "O2^-+O2=>O2+O2+E                                   ","O2^-+O2(A1)=>O2+O2+E                               ",&
   "O2^-+O2(B1)=>O2+O2+E                               ","O2^-+N2=>O2+N2+E                                   ",&
   "O2^-+N2(A3)=>O2+N2+E                               ","O2^-+N2(B3)=>O2+N2+E                               ",&
   "O3^-+O=>O2+O2+E                                    ","NO^-+N=>N2O+E                                      "/
  data reaction_sign(217:252) &
  /"O3^-+N=>NO+O2+E                                    ","N2O^-+N=>NO+N2+E                                   ",&
   "NO2^-+N=>NO+NO+E                                   ","NO3^-+N=>NO+NO2+E                                  ",&
   "NO^-+O=>NO2+E                                      ","N2O^-+O=>NO+NO+E                                   ",&
   "NO2^-+O=>NO+O2+E                                   ","NO3^-+O=>NO+O3+E                                   ",&
   "O3^-+N2(A3)=>O3+N2+E                               ","NO^-+N2(A3)=>NO+N2+E                               ",&
   "N2O^-+N2(A3)=>N2O+N2+E                             ","NO2^-+N2(A3)=>NO2+N2+E                             ",&
   "NO3^-+N2(A3)=>NO3+N2+E                             ","O3^-+N2(B3)=>O3+N2+E                               ",&
   "NO^-+N2(B3)=>NO+N2+E                               ","N2O^-+N2(B3)=>N2O+N2+E                             ",&
   "NO2^-+N2(B3)=>NO2+N2+E                             ","NO3^-+N2(B3)=>NO3+N2+E                             ",&
   "N2(A3)=>N2                                         ","N2(B3)=>N2(A3)                                     ",&
   "N2(A`1)=>N2                                        ","N2(C3)=>N2(B3)                                     ",&
   "O2(A1)=>O2                                         ","O2(B1)=>O2(A1)                                     ",&
   "O2(B1)=>O2                                         ","O2(4.5EV)=>O2                                      ",&
   "N2(A3)+O2=>N2+O+O+1.05D0_EV                        ","N2(A3)+O2=>N2+O2(B1)+4.54D0_EV                     ",&
   "N2(B3)+O2=>N2+O+O+2.23D0_EV                        ","N2(A`1)+O2=>N2+O+O(1D)+1.02D0_EV                   ",&
   "N2(C3)+O2=>N2+O+O+5.91D0_EV                        ","N2(C3)+O2=>N2+O+O(1D)+3.65D0_EV                    ",&
   "N2(C3)+O2=>N2+O+O(1S)+1.97D0_EV                    ","N2(A3)+O=>NO+N(2D)                                 ",&
   "N2(A3)+O=>N2+O(1S)+2.23D0_EV                       ","N2(A3)+N=>N2+N                                     "/
  data reaction_sign(253:288) &
  /"N2(A3)+N=>N2+N(2P)                                 ","N2(A3)+N2=>N2+N2                                   ",&
   "N2(A3)+NO=>N2+NO                                   ","N2(A3)+N2O=>N2+N+NO                                ",&
   "N2(A3)+NO2=>N2+O+NO                                ","N2(A3)+N2(A3)=>N2+N2(B3)+4.00D0_EV                 ",&
   "N2(A3)+N2(A3)=>N2+N2(C3)+1.31D0_EV                 ","N2(B3)+N2=>N2(A3)+N2                               ",&
   "N2(B3)+N2=>N2+N2                                   ","N2(B3)+NO=>N2(A3)+NO                               ",&
   "N2(C3)+N2=>N2(A`1)+N2                              ","N2(A`1)+N2=>N2(B3)+N2                              ",&
   "N2(A`1)+NO=>N2+N+O                                 ","N2(A`1)+N2(A3)=>N4^++E                             ",&
   "N2(A`1)+N2(A`1)=>N4^++E                            ","N2(A3)+H=>N2+H                                     ",&
   "N2(A3)+H2=>N2+2H                                   ","N2(A3)+NH3=>N2+NH3                                 ",&
   "N2(B3)+H2=>N2(A3)+H2                               ","N2(A`1)+H=>N2+H                                    ",&
   "N2(A`1)+H2=>N2+H2                                  ","N+N+N2=>N2(A3)+N2                                  ",&
   "N+N+H2=>N2(A3)+H2                                  ","N+N+NH3=>N2(A3)+NH3                                ",&
   "N+N+N=>N2(A3)+N                                    ","N+N+H=>N2(A3)+H                                    ",&
   "N+N+N2=>N2(B3)+N2                                  ","N+N+H2=>N2(B3)+H2                                  ",&
   "N+N+NH3=>N2(B3)+NH3                                ","N+N+N=>N2(B3)+N                                    ",&
   "N+N+H=>N2(B3)+H                                    ","N(2D)+O=>N+O(1D)                                   ",&
   "N(2D)+O2=>NO+O                                     ","N(2D)+NO=>N2+O                                     ",&
   "N(2D)+N2O=>NO+N2                                   ","N(2D)+N2=>N+N2                                     "/
  data reaction_sign(289:324) &
  /"N(2P)+N=>N+N                                       ","N(2P)+O=>N+O                                       ",&
   "N(2P)+N=>N(2D)+N                                   ","N(2P)+N2=>N+N2                                     ",&
   "N(2P)+N(2D)=>N2^++E                                ","N(2P)+O2=>NO+O                                     ",&
   "N(2P)+NO=>N2(A3)+O                                 ","N(2D)+H2=>H+NH                                     ",&
   "N(2D)+NH3=>NH+NH2                                  ","N(2P)+H2=>H+NH                                     ",&
   "O2(A1)+O=>O2+O                                     ","O2(A1)+N=>NO+O                                     ",&
   "O2(A1)+O2=>O2+O2                                   ","O2(A1)+N2=>O2+N2                                   ",&
   "O2(A1)+NO=>O2+NO                                   ","O2(A1)+O3=>O2+O2+O(1D)                             ",&
   "O2(A1)+O2(A1)=>O2+O2(B1)                           ","O+O3=>O2+O2(A1)                                    ",&
   "O2(B1)+O=>O2(A1)+O                                 ","O2(B1)+O=>O2+O(1D)                                 ",&
   "O2(B1)+O2=>O2(A1)+O2                               ","O2(B1)+N2=>O2(A1)+N2                               ",&
   "O2(B1)+NO=>O2(A1)+NO                               ","O2(B1)+O3=>O2+O2+O                                 ",&
   "O2(4.5EV)+O=>O2+O(1S)                              ","O2(4.5EV)+O2=>O2(B1)+O2(B1)                        ",&
   "O2(4.5EV)+N2=>O2(B1)+N2                            ","O(1D)+O=>O+O+2.26D0_EV                             ",&
   "O(1D)+O2=>O+O2+2.26D0_EV                           ","O(1D)+O2=>O+O2(A1)+1.28D0_EV                       ",&
   "O(1D)+O2=>O+O2(B1)+0.63D0_EV                       ","O(1D)+N2=>O+N2+2.26D0_EV                           ",&
   "O(1D)+O3=>O2+O+O                                   ","O(1D)+O3=>O2+O2                                    ",&
   "O(1D)+NO=>O2+N                                     ","O(1D)+N2O=>NO+NO                                   "/
  data reaction_sign(325:360) &
  /"O(1D)+N2O=>O2+N2                                   ","O(1S)+O=>O(1D)+O+1.68D0_EV                         ",&
   "O(1S)+N=>O+N+3.94D0_EV                             ","O(1S)+O2=>O(1D)+O2+1.68D0_EV                       ",&
   "O(1S)+O2+1.18D0_EV=>O+O+O                          ","O(1S)+N2=>O+N2+3.94D0_EV                           ",&
   "O(1S)+O2(A1)=>O+O2(4.5EV)+0.42D0_EV                ","O(1S)+O2(A1)=>O(1D)+O2(B1)+1.03D0_EV               ",&
   "O(1S)+O2(A1)+0.20D0_EV=>O+O+O                      ","O(1S)+NO=>O+NO+4.92D0_EV                           ",&
   "O(1S)+NO=>O(1D)+NO+1.68D0_EV                       ","O(1S)+O3=>O2+O2+3.94D0_EV                          ",&
   "O(1S)+O3=>O2+O+O(1D)+3.94D0_EV                     ","O(1S)+N2O=>O+N2O+3.94D0_EV                         ",&
   "O(1S)+N2O=>O(1D)+N2O+1.68D0_EV                     ","N+NO=>O+N2                                         ",&
   "N+O2=>O+NO                                         ","N+NO2=>O+O+N2                                      ",&
   "N+NO2=>O+N2O                                       ","N+NO2=>N2+O2                                       ",&
   "N+NO2=>NO+NO                                       ","O+N2=>N+NO                                         ",&
   "O+NO=>N+O2                                         ","O+NO=>NO2                                          ",&
   "O+N2O=>N2+O2                                       ","O+N2O=>NO+NO                                       ",&
   "O+NO2=>NO+O2                                       ","O+NO3=>O2+NO2                                      ",&
   "N2+O2=>O+N2O                                       ","NO+NO=>N+NO2                                       ",&
   "NO+NO=>O+N2O                                       ","NO+NO=>N2+O2                                       ",&
   "NO+O2=>O+NO2                                       ","NO+O3=>O2+NO2                                      ",&
   "NO+N2O=>N2+NO2                                     ","NO+NO3=>NO2+NO2                                    "/
  data reaction_sign(361:396) &
  /"O2+O2=>O+O3                                        ","O2+NO2=>NO+O3                                      ",&
   "NO2+NO2=>NO+NO+O2                                  ","NO2+NO2=>NO+NO3                                    ",&
   "NO2+O3=>O2+NO3                                     ","NO2+NO3=>NO+NO2+O2                                 ",&
   "NO3+O2=>NO2+O3                                     ","NO3+NO3=>O2+NO2+NO2                                ",&
   "N+N=>N2^++E                                        ","N+O=>NO^++E                                        ",&
   "H+NH=>N+H2                                         ","H+NH2=>NH+H2                                       ",&
   "H+NH3=>NH2+H2                                      ","H2+NH2=>NH3+H                                      ",&
   "N+NH=>N2+H                                         ","N+NH2=>N2+2H                                       ",&
   "N+NH2=>N2+H2                                       ","NH+NH=>N+NH2                                       ",&
   "NH+NH=>N2+2H                                       ","NH+NH=>N2+H2                                       ",&
   "NH+NH2=>NH3+N                                      ","H+NH2=>NH3                                         ",&
   "NH2+NH2=>NH+NH3                                    ","N2+N2=>N+N+N2                                      ",&
   "N2+O2=>N+N+O2                                      ","N2+NO=>N+N+NO                                      ",&
   "N2+O=>N+N+O                                        ","N2+N=>N+N+N                                        ",&
   "O2+N2=>O+O+N2                                      ","O2+O2=>O+O+O2                                      ",&
   "O2+O=>O+O+O                                        ","O2+N=>O+O+N                                        ",&
   "O2+NO=>O+O+NO                                      ","NO+N2=>N+O+N2                                      ",&
   "NO+O2=>N+O+O2                                      ","NO+O=>N+O+O                                        "/
  data reaction_sign(397:432) &
  /"NO+N=>N+O+N                                        ","NO+NO=>N+O+NO                                      ",&
   "O3+N2=>O2+O+N2                                     ","O3+O2=>O2+O+O2                                     ",&
   "O3+N=>O2+O+N                                       ","O3+O=>O2+O+O                                       ",&
   "N2O+N2=>N2+O+N2                                    ","N2O+O2=>N2+O+O2                                    ",&
   "N2O+NO=>N2+O+NO                                    ","N2O+N2O=>N2+O+N2O                                  ",&
   "NO2+N2=>NO+O+N2                                    ","NO2+O2=>NO+O+O2                                    ",&
   "NO2+NO=>NO+O+NO                                    ","NO2+NO2=>NO+O+NO2                                  ",&
   "NO3+N2=>NO2+O+N2                                   ","NO3+O2=>NO2+O+O2                                   ",&
   "NO3+NO=>NO2+O+NO                                   ","NO3+N=>NO2+O+N                                     ",&
   "NO3+O=>NO2+O+O                                     ","NO3+N2=>NO+O2+N2                                   ",&
   "NO3+O2=>NO+O2+O2                                   ","NO3+NO=>NO+O2+NO                                   ",&
   "NO3+N=>NO+O2+N                                     ","NO3+O=>NO+O2+O                                     ",&
   "N2O5+ANY_NEUTRAL=>NO2+NO3+ANY_NEUTRAL              ","N+N+N2=>N2+N2                                      ",&
   "N+N+O2=>N2+O2                                      ","N+N+NO=>N2+NO                                      ",&
   "N+N+N=>N2+N                                        ","N+N+O=>N2+O                                        ",&
   "O+O+N2=>O2+N2                                      ","O+O+O2=>O2+O2                                      ",&
   "O+O+N=>O2+N                                        ","O+O+O=>O2+O                                        ",&
   "O+O+NO=>O2+NO                                      ","N+O+N2=>NO+N2                                      "/
  data reaction_sign(433:468) &
  /"N+O+O2=>NO+O2                                      ","N+O+N=>NO+N                                        ",&
   "N+O+O=>NO+O                                        ","N+O+NO=>NO+NO                                      ",&
   "O+O2+N2=>O3+N2                                     ","O+O2+O2=>O3+O2                                     ",&
   "O+O2+NO=>O3+NO                                     ","O+O2+N=>O3+N                                       ",&
   "O+O2+O=>O3+O                                       ","O+N2+ANY_NEUTRAL=>N2O+ANY_NEUTRAL                  ",&
   "O+NO+N2=>NO2+N2                                    ","O+NO+O2=>NO2+O2                                    ",&
   "O+NO+NO=>NO2+NO                                    ","O+NO2+N2=>NO3+N2                                   ",&
   "O+NO2+O2=>NO3+O2                                   ","O+NO2+N=>NO3+N                                     ",&
   "O+NO2+O=>NO3+O                                     ","O+NO2+NO=>NO3+NO                                   ",&
   "NO2+NO3+ANY_NEUTRAL=>N2O5+ANY_NEUTRAL              ","N+N+N2=>N2+N2                                      ",&
   "N+N+H2=>N2+H2                                      ","N+N+NH3=>N2+NH3                                    ",&
   "H+H+H2=>H2+H2                                      ","H+H+N2=>H2+N2                                      ",&
   "N+H+N2=>NH+N2                                      ","N+H+H2=>NH+H2                                      ",&
   "N+H+NH3=>NH+NH3                                    ","N+H2+N2=>NH2+N2                                    ",&
   "N+H2+H2=>NH2+H2                                    ","N+H2+NH3=>NH2+NH3                                  ",&
   "H+NH+N2=>NH2+N2                                    ","H+NH+H2=>NH2+H2                                    ",&
   "H+NH+NH3=>NH2+NH3                                  ","H+NH2+N2=>NH3+N2                                   ",&
   "H+NH2+H2=>NH3+H2                                   ","H+NH2+NH3=>NH3+NH3                                 "/
  data reaction_sign(469:504) &
  /"NH+H2+N2=>NH3+N2                                   ","NH+H2+H2=>NH3+H2                                   ",&
   "NH+H2+NH3=>NH3+NH3                                 ","N^++O=>N+O^+                                       ",&
   "N^++O2=>O2^++N                                     ","N^++O2=>NO^++O                                     ",&
   "N^++O2=>O^++NO                                     ","N^++O3=>NO^++O2                                    ",&
   "N^++NO=>NO^++N                                     ","N^++NO=>N2^++O                                     ",&
   "N^++NO=>O^++N2                                     ","N^++N2O=>NO^++N2                                   ",&
   "O^++N2=>NO^++N                                     ","O^++O2=>O2^++O                                     ",&
   "O^++O3=>O2^++O2                                    ","O^++NO=>NO^++O                                     ",&
   "O^++NO=>O2^++N                                     ","O^++N(2D)=>N^++O                                   ",&
   "O^++N2O=>NO^++NO                                   ","O^++N2O=>N2O^++O                                   ",&
   "O^++N2O=>O2^++N2                                   ","O^++NO2=>NO2^++O                                   ",&
   "N2^++O2=>O2^++N2                                   ","N2^++O=>NO^++N                                     ",&
   "N2^++O3=>O2^++O+N2                                 ","N2^++N=>N^++N2                                     ",&
   "N2^++NO=>NO^++N2                                   ","N2^++N2O=>N2O^++N2                                 ",&
   "N2^++N2O=>NO^++N+N2                                ","O2^++N2=>NO^++NO                                   ",&
   "O2^++N=>NO^++O                                     ","O2^++NO=>NO^++O2                                   ",&
   "O2^++NO2=>NO^++O3                                  ","O2^++NO2=>NO2^++O2                                 ",&
   "N3^++O2=>O2^++N+N2                                 ","N3^++O2=>NO2^++N2                                  "/
  data reaction_sign(505:540) &
  /"N3^++N=>N2^++N2                                    ","N3^++NO=>NO^++N+N2                                 ",&
   "N3^++NO=>N2O^++N2                                  ","NO2^++NO=>NO^++NO2                                 ",&
   "N2O^++NO=>NO^++N2O                                 ","N4^++N2=>N2^++N2+N2                                ",&
   "N4^++O2=>O2^++N2+N2                                ","N4^++O=>O^++N2+N2                                  ",&
   "N4^++N=>N^++N2+N2                                  ","N4^++NO=>NO^++N2+N2                                ",&
   "O4^++N2=>O2^+N2+O2                                 ","O4^++O2=>O2^++O2+O2                                ",&
   "O4^++O2(A1)=>O2^++O2+O2                            ","O4^++O2(B1)=>O2^++O2+O2                            ",&
   "O4^++O=>O2^++O3                                    ","O4^++NO=>NO^++O2+O2                                ",&
   "O2^+N2+N2=>O2^++N2+N2                              ","O2^+N2+O2=>O4^++N2                                 ",&
   "N^++N2+N2=>N3^++N2                                 ","N^++O+ANY_NEUTRAL=>NO^++ANY_NEUTRAL                ",&
   "N^++N+ANY_NEUTRAL=>N2^++ANY_NEUTRAL                ","O^++N2+ANY_NEUTRAL=>NO^++N+ANY_NEUTRAL             ",&
   "O^++O+ANY_NEUTRAL=>O2^++ANY_NEUTRAL                ","O^++N+ANY_NEUTRAL=>NO^++ANY_NEUTRAL                ",&
   "N2^++N2+N2=>N4^++N2                                ","N2^++N+N2=>N3^++N2                                 ",&
   "O2^++O2+O2=>O4^++O2                                ","O2^++N2+N2=>O2^+N2+N2                              ",&
   "N^++H2=>NH^++H                                     ","N^++NH3=>NH2^++NH                                  ",&
   "N^++NH3=>N2H^++H2                                  ","H^++NH3=>NH3^++H                                   ",&
   "H2^++H=>H2+H^+                                     ","H2^++H2=>H3^++H                                    ",&
   "H2^++NH3=>NH3^++H2                                 ","H2^++N2=>N2H^++H                                   "/
  data reaction_sign(541:576) &
  /"NH^++H2=>H3^++N                                    ","NH^++H2=>NH2^++H                                   ",&
   "NH^++NH3=>NH^++NH3                                 ","NH^++NH3=>NH4^++N                                  ",&
   "NH^++N2=>N2H^++N                                   ","NH2^++H2=>NH3^++H                                  ",&
   "NH2^++NH3=>NH3^++NH2                               ","NH2^++NH3=>NH4^++NH                                ",&
   "N2H^++NH3=>NH4^++N2                                ","N2^++NH3=>N2+NH3^+                                 ",&
   "N^++NH3=>N+NH3^+                                   ","O2^++NH3=>O2+NH3^+                                 ",&
   "H3O^++NH3=>NH4^++H2O                               ","NH3^++NH3=>NH4^++NH2                               ",&
   "NH3^++NH2=>NH3+NH2^+                               ","O^-+O2(A1)=>O2^-+O                                 ",&
   "O^-+O3=>O3^-+O                                     ","O^-+NO2=>NO2^-+O                                   ",&
   "O^-+N2O=>NO^-+NO                                   ","O^-+N2O=>N2O^-+O                                   ",&
   "O2^-+O=>O^-+O2                                     ","O2^-+O3=>O3^-+O2                                   ",&
   "O2^-+NO2=>NO2^-+O2                                 ","O2^-+NO3=>NO3^-+O2                                 ",&
   "O3^-+O=>O2^-+O2                                    ","O3^-+NO=>NO3^-+O                                   ",&
   "O3^-+NO=>NO2^-+O2                                  ","O3^-+NO2=>NO2^-+O3                                 ",&
   "O3^-+NO2=>NO3^-+O2                                 ","O3^-+NO3=>NO3^-+O3                                 ",&
   "NO^-+O2=>O2^-+NO                                   ","NO^-+NO2=>NO2^-+NO                                 ",&
   "NO^-+N2O=>NO2^-+N2                                 ","NO2^-+O3=>NO3^-+O2                                 ",&
   "NO2^-+NO2=>NO3^-+NO                                ","NO2^-+NO3=>NO3^-+NO2                               "/
  data reaction_sign(577:612) &
  /"NO2^-+N2O5=>NO3^-+NO2+NO2                          ","NO3^-+NO=>NO2^-+NO2                                ",&
   "O4^-+N2=>O2^-+O2+N2                                ","O4^-+O2=>O2^-+O2+O2                                ",&
   "O4^-+O=>O3^-+O2                                    ","O4^-+O=>O^-+O2+O2                                  ",&
   "O4^-+O2(A1)=>O2^-+O2+O2                            ","O4^-+O2(B1)=>O2^-+O2+O2                            ",&
   "O4^-+NO=>NO3^-+O2                                  ","O^-+O2+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL               ",&
   "O^-+NO+ANY_NEUTRAL=>NO2^-+ANY_NEUTRAL              ","O2^-+O2+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL              ",&
   "O^-+N^+=>O+N+1.309D1_EV                            ","O^-+N2^+=>O+N2+1.414D1_EV                          ",&
   "O^-+O^+=>O+O+1.218D1_EV                            ","O^-+O2^+=>O+O2+1.063D1_EV                          ",&
   "O^-+NO^+=>O+NO+7.825D0_EV                          ","O^-+N2O^+=>O+N2O+1.145D1_EV                        ",&
   "O^-+NO2^+=>O+NO2+8.147D0_EV                        ","O2^-+N^+=>O2+N+1.409D1_EV                          ",&
   "O2^-+N2^+=>O2+N2+1.513D1_EV                        ","O2^-+O^+=>O2+O+1.317D1_EV                          ",&
   "O2^-+O2^+=>O2+O2+1.162D1_EV                        ","O2^-+NO^+=>O2+NO+8.816D0_EV                        ",&
   "O2^-+N2O^+=>O2+N2O+1.244D1_EV                      ","O2^-+NO2^+=>O2+NO2+9.138D0_EV                      ",&
   "O3^-+N^+=>O3+N+1.243D1_EV                          ","O3^-+N2^+=>O3+N2+1.348D1_EV                        ",&
   "O3^-+O^+=>O3+O+1.152D1_EV                          ","O3^-+O2^+=>O3+O2+9.967D0_EV                        ",&
   "O3^-+NO^+=>O3+NO+7.161D0_EV                        ","O3^-+N2O^+=>O3+N2O+1.079D1_EV                      ",&
   "O3^-+NO2^+=>O3+NO2+7.483D0_EV                      ","NO^-+N^+=>NO+N                                     ",&
   "NO^-+N2^+=>NO+N2                                   ","NO^-+O^+=>NO+O                                     "/
  data reaction_sign(613:648) &
  /"NO^-+O2^+=>NO+O2                                   ","NO^-+NO^+=>NO+NO                                   ",&
   "NO^-+N2O^+=>NO+N2O                                 ","NO^-+NO2^+=>NO+NO2                                 ",&
   "N2O^-+N^+=>N2O+N+1.431D1_EV                        ","N2O^-+N2^+=>N2O+N2+1.536D1_EV                      ",&
   "N2O^-+O^+=>N2O+O+1.340D1_EV                        ","N2O^-+O2^+=>N2O+O2+1.185D1_EV                      ",&
   "N2O^-+NO^+=>N2O+NO+9.044D0_EV                      ","N2O^-+N2O^+=>N2O+N2O+1.267D1_EV                    ",&
   "N2O^-+NO2^+=>N2O+NO2+9.366D0_EV                    ","NO2^-+N^+=>NO2+N+1.226D1_EV                        ",&
   "NO2^-+N2^+=>NO2+N2+1.331D1_EV                      ","NO2^-+O^+=>NO2+O+1.135D1_EV                        ",&
   "NO2^-+O2^+=>NO2+O2+9.797D0_EV                      ","NO2^-+NO^+=>NO2+NO+6.991D0_EV                      ",&
   "NO2^-+N2O^+=>NO2+N2O+1.062D1_EV                    ","NO2^-+NO2^+=>NO2+NO2+7.313D0_EV                    ",&
   "NO3^-+N^+=>NO3+N+1.060D1_EV                        ","NO3^-+N2^+=>NO3+N2+1.164D1_EV                      ",&
   "NO3^-+O^+=>NO3+O+9.681D0_EV                        ","NO3^-+O2^+=>NO3+O2+8.133D0_EV                      ",&
   "NO3^-+NO^+=>NO3+NO+5.327D0_EV                      ","NO3^-+N2O^+=>NO3+N2O+8.952D0_EV                    ",&
   "NO3^-+NO2^+=>NO3+NO2+5.649D0_EV                    ","O^-+N2^+=>O+N+N                                    ",&
   "O^-+N3^+=>O+N+N2                                   ","O^-+N4^+=>O+N2+N2                                  ",&
   "O^-+O2^+=>O+O+O                                    ","O^-+O4^+=>O+O2+O2                                  ",&
   "O^-+NO^+=>O+N+O                                    ","O^-+N2O^+=>O+N2+O                                  ",&
   "O^-+NO2^+=>O+N+O2                                  ","O^-+O2^+N2=>O+O2+N2                                ",&
   "O2^-+N2^+=>O2+N+N                                  ","O2^-+N3^+=>O2+N+N2                                 "/
  data reaction_sign(649:684) &
  /"O2^-+N4^+=>O2+N2+N2                                ","O2^-+O2^+=>O2+O+O                                  ",&
   "O2^-+O4^+=>O2+O2+O2                                ","O2^-+NO^+=>O2+N+O                                  ",&
   "O2^-+N2O^+=>O2+N2+O                                ","O2^-+NO2^+=>O2+N+O2                                ",&
   "O2^-+O2^+N2=>O2+O2+N2                              ","O3^-+N2^+=>O3+N+N                                  ",&
   "O3^-+N3^+=>O3+N+N2                                 ","O3^-+N4^+=>O3+N2+N2                                ",&
   "O3^-+O2^+=>O3+O+O                                  ","O3^-+O4^+=>O3+O2+O2                                ",&
   "O3^-+NO^+=>O3+N+O                                  ","O3^-+N2O^+=>O3+N2+O                                ",&
   "O3^-+NO2^+=>O3+N+O2                                ","O3^-+O2^+N2=>O3+O2+N2                              ",&
   "NO^-+N2^+=>NO+N+N                                  ","NO^-+N3^+=>NO+N+N2                                 ",&
   "NO^-+N4^+=>NO+N2+N2                                ","NO^-+O2^+=>NO+O+O                                  ",&
   "NO^-+O4^+=>NO+O2+O2                                ","NO^-+NO^+=>NO+N+O                                  ",&
   "NO^-+N2O^+=>NO+N2+O                                ","NO^-+NO2^+=>NO+N+O2                                ",&
   "NO^-+O2^+N2=>NO+O2+N2                              ","N2O^-+N2^+=>N2O+N+N                                ",&
   "N2O^-+N3^+=>N2O+N+N2                               ","N2O^-+N4^+=>N2O+N2+N2                              ",&
   "N2O^-+O2^+=>N2O+O+O                                ","N2O^-+O4^+=>N2O+O2+O2                              ",&
   "N2O^-+NO^+=>N2O+N+O                                ","N2O^-+N2O^+=>N2O+N2+O                              ",&
   "N2O^-+NO2^+=>N2O+N+O2                              ","N2O^-+O2^+N2=>N2O+O2+N2                            ",&
   "NO2^-+N2^+=>NO2+N+N                                ","NO2^-+N3^+=>NO2+N+N2                               "/
  data reaction_sign(685:720) &
  /"NO2^-+N4^+=>NO2+N2+N2                              ","NO2^-+O2^+=>NO2+O+O                                ",&
   "NO2^-+O4^+=>NO2+O2+O2                              ","NO2^-+NO^+=>NO2+N+O                                ",&
   "NO2^-+N2O^+=>NO2+N2+O                              ","NO2^-+NO2^+=>NO2+N+O2                              ",&
   "NO2^-+O2^+N2=>NO2+O2+N2                            ","NO3^-+N2^+=>NO3+N+N                                ",&
   "NO3^-+N3^+=>NO3+N+N2                               ","NO3^-+N4^+=>NO3+N2+N2                              ",&
   "NO3^-+O2^+=>NO3+O+O                                ","NO3^-+O4^+=>NO3+O2+O2                              ",&
   "NO3^-+NO^+=>NO3+N+O                                ","NO3^-+N2O^+=>NO3+N2+O                              ",&
   "NO3^-+NO2^+=>NO3+N+O2                              ","NO3^-+O2^+N2=>NO3+O2+N2                            ",&
   "O4^-+N^+=>O2+O2+N                                  ","O4^-+N2^+=>O2+O2+N2                                ",&
   "O4^-+O^+=>O2+O2+O                                  ","O4^-+O2^+=>O2+O2+O2                                ",&
   "O4^-+NO^+=>O2+O2+NO                                ","O4^-+N2O^+=>O2+O2+N2O                              ",&
   "O4^-+NO2^+=>O2+O2+NO2                              ","O4^-+N3^+=>O2+O2+N2+N                              ",&
   "O4^-+N4^+=>O2+O2+N2+N2                             ","O4^-+O4^+=>O2+O2+O2+O2                             ",&
   "O4^-+O2^+N2=>O2+O2+O2+N2                           ","O^-+N^++ANY_NEUTRAL=>O+N+ANY_NEUTRAL+1.309D1_EV    ",&
   "O^-+N2^++ANY_NEUTRAL=>O+N2+ANY_NEUTRAL+1.414D1_EV  ","O^-+O^++ANY_NEUTRAL=>O+O+ANY_NEUTRAL+1.218D1_EV    ",&
   "O^-+O2^++ANY_NEUTRAL=>O+O2+ANY_NEUTRAL+1.063D1_EV  ","O^-+NO^++ANY_NEUTRAL=>O+NO+ANY_NEUTRAL+7.825D0_EV  ",&
   "O2^-+N^++ANY_NEUTRAL=>O2+N+ANY_NEUTRAL+1.409D1_EV  ","O2^-+N2^++ANY_NEUTRAL=>O2+N2+ANY_NEUTRAL+1.513D1_EV",&
   "O2^-+O^++ANY_NEUTRAL=>O2+O+ANY_NEUTRAL+1.317D1_EV  ","O2^-+O2^++ANY_NEUTRAL=>O2+O2+ANY_NEUTRAL+1.162D1_EV"/
  data reaction_sign(721:756) &
  /"O2^-+NO^++ANY_NEUTRAL=>O2+NO+ANY_NEUTRAL+8.816D0_EV","O^-+N^++ANY_NEUTRAL=>NO+ANY_NEUTRAL                ",&
   "O^-+N2^++ANY_NEUTRAL=>N2O+ANY_NEUTRAL              ","O^-+O^++ANY_NEUTRAL=>O2+ANY_NEUTRAL                ",&
   "O^-+O2^++ANY_NEUTRAL=>O3+ANY_NEUTRAL               ","O^-+NO^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL              ",&
   "O2^-+N^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL              ","O2^-+O^++ANY_NEUTRAL=>O3+ANY_NEUTRAL               ",&
   "O2^-+NO^++ANY_NEUTRAL=>NO3+ANY_NEUTRAL             ","NH3^++O2^-=>NH3+O2+9.622D0_EV                      ",&
   "NH3^++O^-=>NH3+O+8.631D0_EV                        ","NH3^++O3^-=>NH3+O3+7.967D0_EV                      ",&
   "NH3^++NO^-=>NH3+NO                                 ","NH3^++NO2^-=>NH3+NO2+7.797D0_EV                    ",&
   "NH3^++NO3^-=>NH3+NO3+6.133D0_EV                    ","NH3^++N2O^-=>NH3+N2O+9.850D0_EV                    ",&
   "NO3^-+NH4^+=>NO3+NH3+H                             ","H^-+H2^+=>3H                                       ",&
   "H^-+H3^+=>H2+2H                                    ","H^-+N2^+=>N2+H+1.483D1_EV                          ",&
   "H^-+N4^+=>2N2+H                                    ","H^-+N2H^+=>H2+N2                                   ",&
   "H^-+H2^++N2=>H2+H+N2+1.467D1_EV                    ","H^-+H2^++H2=>H2+H+H2+1.467D1_EV                    ",&
   "H^-+H2^++NH3=>H2+H+NH3+1.467D1_EV                  ","H^-+H3^++N2=>2H2+N2                                ",&
   "H^-+H3^++H2=>2H2+H2                                ","H^-+H3^++NH3=>2H2+NH3                              ",&
   "H^-+N2^++N2=>N2+H+N2+1.483D1_EV                    ","H^-+N2^++H2=>N2+H+H2+1.483D1_EV                    ",&
   "H^-+N2^++NH3=>N2+H+NH3+1.483D1_EV                  ","H^-+N4^++N2=>2N2+H+N2                              ",&
   "H^-+N4^++H2=>2N2+H+H2                              ","H^-+N4^++NH3=>2N2+H+NH3                            ",&
   "H^-+N2H^++N2=>H2+N2+N2+7.045D0_EV                  ","H^-+N2H^++H2=>H2+N2+H2+7.045D0_EV                  "/
  data reaction_sign(757:792) &
  /"H^-+N2H^++NH3=>H2+N2+NH3+7.045D0_EV                ","O3^-+N^++ANY_NEUTRAL=>O3+N+ANY_NEUTRAL             ",&
   "O3^-+N2^++ANY_NEUTRAL=>O3+N2+ANY_NEUTRAL           ","O3^-+O^++ANY_NEUTRAL=>O3+O+ANY_NEUTRAL             ",&
   "O3^-+O2^++ANY_NEUTRAL=>O3+O2+ANY_NEUTRAL           ","O3^-+NO^++ANY_NEUTRAL=>O3+NO+ANY_NEUTRAL           ",&
   "O3^-+N2O^++ANY_NEUTRAL=>O3+N2O+ANY_NEUTRAL         ","O3^-+NO2^++ANY_NEUTRAL=>O3+NO2+ANY_NEUTRAL         ",&
   "NO^-+N^++ANY_NEUTRAL=>NO+N+ANY_NEUTRAL             ","NO^-+N2^++ANY_NEUTRAL=>NO+N2+ANY_NEUTRAL           ",&
   "NO^-+O^++ANY_NEUTRAL=>NO+O+ANY_NEUTRAL             ","NO^-+O2^++ANY_NEUTRAL=>NO+O2+ANY_NEUTRAL           ",&
   "NO^-+NO^++ANY_NEUTRAL=>NO+NO+ANY_NEUTRAL           ","NO^-+N2O^++ANY_NEUTRAL=>NO+N2O+ANY_NEUTRAL         ",&
   "NO^-+NO2^++ANY_NEUTRAL=>NO+NO2+ANY_NEUTRAL         ","N2O^-+N^++ANY_NEUTRAL=>N2O+N+ANY_NEUTRAL           ",&
   "N2O^-+N2^++ANY_NEUTRAL=>N2O+N2+ANY_NEUTRAL         ","N2O^-+O^++ANY_NEUTRAL=>N2O+O+ANY_NEUTRAL           ",&
   "N2O^-+O2^++ANY_NEUTRAL=>N2O+O2+ANY_NEUTRAL         ","N2O^-+NO^++ANY_NEUTRAL=>N2O+NO+ANY_NEUTRAL         ",&
   "N2O^-+N2O^++ANY_NEUTRAL=>N2O+N2O+ANY_NEUTRAL       ","N2O^-+NO2^++ANY_NEUTRAL=>N2O+NO2+ANY_NEUTRAL       ",&
   "NO2^-+N^++ANY_NEUTRAL=>NO2+N+ANY_NEUTRAL           ","NO2^-+N2^++ANY_NEUTRAL=>NO2+N2+ANY_NEUTRAL         ",&
   "NO2^-+O^++ANY_NEUTRAL=>NO2+O+ANY_NEUTRAL           ","NO2^-+O2^++ANY_NEUTRAL=>NO2+O2+ANY_NEUTRAL         ",&
   "NO2^-+NO^++ANY_NEUTRAL=>NO2+NO+ANY_NEUTRAL         ","NO2^-+N2O^++ANY_NEUTRAL=>NO2+N2O+ANY_NEUTRAL       ",&
   "NO2^-+NO2^++ANY_NEUTRAL=>NO2+NO2+ANY_NEUTRAL       ","NO3^-+N^++ANY_NEUTRAL=>NO3+N+ANY_NEUTRAL           ",&
   "NO3^-+N2^++ANY_NEUTRAL=>NO3+N2+ANY_NEUTRAL         ","NO3^-+O^++ANY_NEUTRAL=>NO3+O+ANY_NEUTRAL           ",&
   "NO3^-+O2^++ANY_NEUTRAL=>NO3+O2+ANY_NEUTRAL         ","NO3^-+NO^++ANY_NEUTRAL=>NO3+NO+ANY_NEUTRAL         ",&
   "NO3^-+N2O^++ANY_NEUTRAL=>NO3+N2O+ANY_NEUTRAL       ","NO3^-+NO2^++ANY_NEUTRAL=>NO3+NO2+ANY_NEUTRAL       "/
  data bolsig_species(1:bolsig_species_max) &
  /"N2     ","N2(V1) ","N2(V2) ","N2(V3) ","N2(V4) ","N2(V5) ","N2(V6) ","N2(V7) ","N2(V8) ","N2(V9) ","N2(V10)","N2(A3) ",&
   "O2     ","O2(A1) ","N      ","O      ","NO     ","O3     ","N2O    ","NH3    "/
  data bolsig_addsect(1:2,1:bolsig_addsect_max) &
  / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1,10, 1,11/
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
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0" // ") INIT:"
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
  integer :: i,j
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
    i = size(ELEC_EEDF(:,:),dim=2)
    call ZDPlasKin_bolsig_GetEEDF(ELEC_EEDF(1:2,:),j)
    if(lprint .and. i < j) write(*,"(A)") "ZDPlasKin WARNING: small size of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
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
111 format(i2,1x,A9)
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
211 format(i3,1x,A51)
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
311 format(561x,79(1x,i9))
312 format(A3,1x,A51,1x,79(1x,A9))
313 format(i3,1x,A51,1x,79(1x,1pd9.2))
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
      write(ifile_unit,"(792(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A12,79(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,792(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 4 )
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
      write(ifile_unit,"(80(1pe13.4))") densav(:,2)
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
      write(ifile_unit,"(793(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(01,001) = - reac_rate_local(001) 
  reac_source_local(02,001) = + reac_rate_local(001) 
  reac_source_local(01,002) = - reac_rate_local(002) 
  reac_source_local(03,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = - reac_rate_local(003) 
  reac_source_local(04,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = - reac_rate_local(004) 
  reac_source_local(05,004) = + reac_rate_local(004) 
  reac_source_local(01,005) = - reac_rate_local(005) 
  reac_source_local(06,005) = + reac_rate_local(005) 
  reac_source_local(01,006) = - reac_rate_local(006) 
  reac_source_local(07,006) = + reac_rate_local(006) 
  reac_source_local(01,007) = - reac_rate_local(007) 
  reac_source_local(08,007) = + reac_rate_local(007) 
  reac_source_local(01,008) = - reac_rate_local(008) 
  reac_source_local(09,008) = + reac_rate_local(008) 
  reac_source_local(01,009) = - reac_rate_local(009) 
  reac_source_local(10,009) = + reac_rate_local(009) 
  reac_source_local(01,010) = - reac_rate_local(010) 
  reac_source_local(11,010) = + reac_rate_local(010) 
  reac_source_local(01,011) = - reac_rate_local(011) 
  reac_source_local(12,011) = + reac_rate_local(011) 
  reac_source_local(01,012) = - reac_rate_local(012) 
  reac_source_local(13,012) = + reac_rate_local(012) 
  reac_source_local(01,013) = - reac_rate_local(013) 
  reac_source_local(14,013) = + reac_rate_local(013) 
  reac_source_local(01,014) = - reac_rate_local(014) 
  reac_source_local(15,014) = + reac_rate_local(014) 
  reac_source_local(01,015) = - reac_rate_local(015) 
  reac_source_local(16,015) = + reac_rate_local(015) 
  reac_source_local(01,016) = + reac_rate_local(016) 
  reac_source_local(02,016) = - reac_rate_local(016) 
  reac_source_local(02,017) = - reac_rate_local(017) 
  reac_source_local(03,017) = + reac_rate_local(017) 
  reac_source_local(02,018) = - reac_rate_local(018) 
  reac_source_local(04,018) = + reac_rate_local(018) 
  reac_source_local(02,019) = - reac_rate_local(019) 
  reac_source_local(05,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(03,020) = - reac_rate_local(020) 
  reac_source_local(02,021) = + reac_rate_local(021) 
  reac_source_local(03,021) = - reac_rate_local(021) 
  reac_source_local(03,022) = - reac_rate_local(022) 
  reac_source_local(04,022) = + reac_rate_local(022) 
  reac_source_local(03,023) = - reac_rate_local(023) 
  reac_source_local(05,023) = + reac_rate_local(023) 
  reac_source_local(03,024) = - reac_rate_local(024) 
  reac_source_local(06,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(04,025) = - reac_rate_local(025) 
  reac_source_local(02,026) = + reac_rate_local(026) 
  reac_source_local(04,026) = - reac_rate_local(026) 
  reac_source_local(03,027) = + reac_rate_local(027) 
  reac_source_local(04,027) = - reac_rate_local(027) 
  reac_source_local(04,028) = - reac_rate_local(028) 
  reac_source_local(05,028) = + reac_rate_local(028) 
  reac_source_local(04,029) = - reac_rate_local(029) 
  reac_source_local(06,029) = + reac_rate_local(029) 
  reac_source_local(04,030) = - reac_rate_local(030) 
  reac_source_local(07,030) = + reac_rate_local(030) 
  reac_source_local(01,031) = + reac_rate_local(031) 
  reac_source_local(05,031) = - reac_rate_local(031) 
  reac_source_local(02,032) = + reac_rate_local(032) 
  reac_source_local(05,032) = - reac_rate_local(032) 
  reac_source_local(03,033) = + reac_rate_local(033) 
  reac_source_local(05,033) = - reac_rate_local(033) 
  reac_source_local(04,034) = + reac_rate_local(034) 
  reac_source_local(05,034) = - reac_rate_local(034) 
  reac_source_local(01,035) = + reac_rate_local(035) 
  reac_source_local(06,035) = - reac_rate_local(035) 
  reac_source_local(03,036) = + reac_rate_local(036) 
  reac_source_local(06,036) = - reac_rate_local(036) 
  reac_source_local(04,037) = + reac_rate_local(037) 
  reac_source_local(06,037) = - reac_rate_local(037) 
  reac_source_local(01,038) = + reac_rate_local(038) 
  reac_source_local(07,038) = - reac_rate_local(038) 
  reac_source_local(04,039) = + reac_rate_local(039) 
  reac_source_local(07,039) = - reac_rate_local(039) 
  reac_source_local(01,040) = + reac_rate_local(040) 
  reac_source_local(08,040) = - reac_rate_local(040) 
  reac_source_local(01,041) = + reac_rate_local(041) 
  reac_source_local(09,041) = - reac_rate_local(041) 
  reac_source_local(01,042) = + reac_rate_local(042) 
  reac_source_local(10,042) = - reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) 
  reac_source_local(11,043) = - reac_rate_local(043) 
  reac_source_local(30,044) = - reac_rate_local(044) 
  reac_source_local(31,044) = + reac_rate_local(044) 
  reac_source_local(30,045) = - reac_rate_local(045) 
  reac_source_local(31,045) = + reac_rate_local(045) 
  reac_source_local(30,046) = - reac_rate_local(046) 
  reac_source_local(32,046) = + reac_rate_local(046) 
  reac_source_local(30,047) = - reac_rate_local(047) 
  reac_source_local(32,047) = + reac_rate_local(047) 
  reac_source_local(30,048) = - reac_rate_local(048) 
  reac_source_local(33,048) = + reac_rate_local(048) 
  reac_source_local(30,049) = - reac_rate_local(049) 
  reac_source_local(34,049) = + reac_rate_local(049) 
  reac_source_local(01,050) = + reac_rate_local(050) 
  reac_source_local(02,050) = - reac_rate_local(050) 
  reac_source_local(02,051) = + reac_rate_local(051) 
  reac_source_local(03,051) = - reac_rate_local(051) 
  reac_source_local(03,052) = + reac_rate_local(052) 
  reac_source_local(04,052) = - reac_rate_local(052) 
  reac_source_local(04,053) = + reac_rate_local(053) 
  reac_source_local(05,053) = - reac_rate_local(053) 
  reac_source_local(05,054) = + reac_rate_local(054) 
  reac_source_local(06,054) = - reac_rate_local(054) 
  reac_source_local(06,055) = + reac_rate_local(055) 
  reac_source_local(07,055) = - reac_rate_local(055) 
  reac_source_local(07,056) = + reac_rate_local(056) 
  reac_source_local(08,056) = - reac_rate_local(056) 
  reac_source_local(08,057) = + reac_rate_local(057) 
  reac_source_local(09,057) = - reac_rate_local(057) 
  reac_source_local(01,058) = - reac_rate_local(058) 
  reac_source_local(02,058) = + reac_rate_local(058) 
  reac_source_local(02,059) = - reac_rate_local(059) 
  reac_source_local(03,059) = + reac_rate_local(059) 
  reac_source_local(03,060) = - reac_rate_local(060) 
  reac_source_local(04,060) = + reac_rate_local(060) 
  reac_source_local(04,061) = - reac_rate_local(061) 
  reac_source_local(05,061) = + reac_rate_local(061) 
  reac_source_local(05,062) = - reac_rate_local(062) 
  reac_source_local(06,062) = + reac_rate_local(062) 
  reac_source_local(06,063) = - reac_rate_local(063) 
  reac_source_local(07,063) = + reac_rate_local(063) 
  reac_source_local(07,064) = - reac_rate_local(064) 
  reac_source_local(08,064) = + reac_rate_local(064) 
  reac_source_local(08,065) = - reac_rate_local(065) 
  reac_source_local(09,065) = + reac_rate_local(065) 
  reac_source_local(01,066) = + reac_rate_local(066) 
  reac_source_local(02,066) = - reac_rate_local(066) 
  reac_source_local(02,067) = + reac_rate_local(067) 
  reac_source_local(03,067) = - reac_rate_local(067) 
  reac_source_local(03,068) = + reac_rate_local(068) 
  reac_source_local(04,068) = - reac_rate_local(068) 
  reac_source_local(04,069) = + reac_rate_local(069) 
  reac_source_local(05,069) = - reac_rate_local(069) 
  reac_source_local(05,070) = + reac_rate_local(070) 
  reac_source_local(06,070) = - reac_rate_local(070) 
  reac_source_local(06,071) = + reac_rate_local(071) 
  reac_source_local(07,071) = - reac_rate_local(071) 
  reac_source_local(07,072) = + reac_rate_local(072) 
  reac_source_local(08,072) = - reac_rate_local(072) 
  reac_source_local(08,073) = + reac_rate_local(073) 
  reac_source_local(09,073) = - reac_rate_local(073) 
  reac_source_local(01,074) = - reac_rate_local(074) 
  reac_source_local(02,074) = + reac_rate_local(074) 
  reac_source_local(02,075) = - reac_rate_local(075) 
  reac_source_local(03,075) = + reac_rate_local(075) 
  reac_source_local(03,076) = - reac_rate_local(076) 
  reac_source_local(04,076) = + reac_rate_local(076) 
  reac_source_local(04,077) = - reac_rate_local(077) 
  reac_source_local(05,077) = + reac_rate_local(077) 
  reac_source_local(05,078) = - reac_rate_local(078) 
  reac_source_local(06,078) = + reac_rate_local(078) 
  reac_source_local(06,079) = - reac_rate_local(079) 
  reac_source_local(07,079) = + reac_rate_local(079) 
  reac_source_local(07,080) = - reac_rate_local(080) 
  reac_source_local(08,080) = + reac_rate_local(080) 
  reac_source_local(08,081) = - reac_rate_local(081) 
  reac_source_local(09,081) = + reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(02,082) = - reac_rate_local(082) 
  reac_source_local(02,083) = + reac_rate_local(083) 
  reac_source_local(03,083) = - reac_rate_local(083) 
  reac_source_local(03,084) = + reac_rate_local(084) 
  reac_source_local(04,084) = - reac_rate_local(084) 
  reac_source_local(04,085) = + reac_rate_local(085) 
  reac_source_local(05,085) = - reac_rate_local(085) 
  reac_source_local(05,086) = + reac_rate_local(086) 
  reac_source_local(06,086) = - reac_rate_local(086) 
  reac_source_local(06,087) = + reac_rate_local(087) 
  reac_source_local(07,087) = - reac_rate_local(087) 
  reac_source_local(07,088) = + reac_rate_local(088) 
  reac_source_local(08,088) = - reac_rate_local(088) 
  reac_source_local(08,089) = + reac_rate_local(089) 
  reac_source_local(09,089) = - reac_rate_local(089) 
  reac_source_local(01,090) = - reac_rate_local(090) 
  reac_source_local(02,090) = + reac_rate_local(090) 
  reac_source_local(02,091) = - reac_rate_local(091) 
  reac_source_local(03,091) = + reac_rate_local(091) 
  reac_source_local(03,092) = - reac_rate_local(092) 
  reac_source_local(04,092) = + reac_rate_local(092) 
  reac_source_local(04,093) = - reac_rate_local(093) 
  reac_source_local(05,093) = + reac_rate_local(093) 
  reac_source_local(05,094) = - reac_rate_local(094) 
  reac_source_local(06,094) = + reac_rate_local(094) 
  reac_source_local(06,095) = - reac_rate_local(095) 
  reac_source_local(07,095) = + reac_rate_local(095) 
  reac_source_local(07,096) = - reac_rate_local(096) 
  reac_source_local(08,096) = + reac_rate_local(096) 
  reac_source_local(08,097) = - reac_rate_local(097) 
  reac_source_local(09,097) = + reac_rate_local(097) 
  reac_source_local(30,098) = + reac_rate_local(098) 
  reac_source_local(31,098) = - reac_rate_local(098) 
  reac_source_local(31,099) = + reac_rate_local(099) 
  reac_source_local(32,099) = - reac_rate_local(099) 
  reac_source_local(32,100) = + reac_rate_local(100) 
  reac_source_local(33,100) = - reac_rate_local(100) 
  reac_source_local(33,101) = + reac_rate_local(101) 
  reac_source_local(34,101) = - reac_rate_local(101) 
  reac_source_local(30,102) = - reac_rate_local(102) 
  reac_source_local(31,102) = + reac_rate_local(102) 
  reac_source_local(31,103) = - reac_rate_local(103) 
  reac_source_local(32,103) = + reac_rate_local(103) 
  reac_source_local(32,104) = - reac_rate_local(104) 
  reac_source_local(33,104) = + reac_rate_local(104) 
  reac_source_local(33,105) = - reac_rate_local(105) 
  reac_source_local(34,105) = + reac_rate_local(105) 
  reac_source_local(30,106) = + reac_rate_local(106) 
  reac_source_local(31,106) = - reac_rate_local(106) 
  reac_source_local(31,107) = + reac_rate_local(107) 
  reac_source_local(32,107) = - reac_rate_local(107) 
  reac_source_local(32,108) = + reac_rate_local(108) 
  reac_source_local(33,108) = - reac_rate_local(108) 
  reac_source_local(33,109) = + reac_rate_local(109) 
  reac_source_local(34,109) = - reac_rate_local(109) 
  reac_source_local(30,110) = - reac_rate_local(110) 
  reac_source_local(31,110) = + reac_rate_local(110) 
  reac_source_local(31,111) = - reac_rate_local(111) 
  reac_source_local(32,111) = + reac_rate_local(111) 
  reac_source_local(32,112) = - reac_rate_local(112) 
  reac_source_local(33,112) = + reac_rate_local(112) 
  reac_source_local(33,113) = - reac_rate_local(113) 
  reac_source_local(34,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = - reac_rate_local(114) 
  reac_source_local(17,114) = + reac_rate_local(114) 
  reac_source_local(01,115) = - reac_rate_local(115) 
  reac_source_local(17,115) = + reac_rate_local(115) 
  reac_source_local(01,116) = - reac_rate_local(116) 
  reac_source_local(17,116) = + reac_rate_local(116) 
  reac_source_local(01,117) = - reac_rate_local(117) 
  reac_source_local(18,117) = + reac_rate_local(117) 
  reac_source_local(01,118) = - reac_rate_local(118) 
  reac_source_local(18,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = - reac_rate_local(119) 
  reac_source_local(18,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = - reac_rate_local(120) 
  reac_source_local(18,120) = + reac_rate_local(120) 
  reac_source_local(01,121) = - reac_rate_local(121) 
  reac_source_local(18,121) = + reac_rate_local(121) 
  reac_source_local(01,122) = - reac_rate_local(122) 
  reac_source_local(18,122) = + reac_rate_local(122) 
  reac_source_local(01,123) = - reac_rate_local(123) 
  reac_source_local(18,123) = + reac_rate_local(123) 
  reac_source_local(01,124) = - reac_rate_local(124) 
  reac_source_local(19,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = - reac_rate_local(125) 
  reac_source_local(19,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = - reac_rate_local(126) 
  reac_source_local(19,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = - reac_rate_local(127) 
  reac_source_local(19,127) = + reac_rate_local(127) 
  reac_source_local(01,128) = - reac_rate_local(128) 
  reac_source_local(19,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = - reac_rate_local(129) 
  reac_source_local(19,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = - reac_rate_local(130) 
  reac_source_local(20,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = - reac_rate_local(131) 
  reac_source_local(20,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = - reac_rate_local(132) 
  reac_source_local(20,132) = + reac_rate_local(132) 
  reac_source_local(01,133) = - reac_rate_local(133) 
  reac_source_local(21,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = - reac_rate_local(134) 
  reac_source_local(21,134) = + reac_rate_local(134) 
  reac_source_local(01,135) = - reac_rate_local(135) 
  reac_source_local(21,135) = + reac_rate_local(135) 
  reac_source_local(01,136) = - reac_rate_local(136) 
  reac_source_local(21,136) = + reac_rate_local(136) 
  reac_source_local(01,137) = - reac_rate_local(137) 
  reac_source_local(21,137) = + reac_rate_local(137) 
  reac_source_local(01,138) = - reac_rate_local(138) 
  reac_source_local(21,138) = + reac_rate_local(138) 
  reac_source_local(01,139) = - reac_rate_local(139) 
  reac_source_local(22,139) = + reac_rate_local(139) 
  reac_source_local(01,140) = - reac_rate_local(140) 
  reac_source_local(22,140) = + reac_rate_local(140) 
  reac_source_local(01,141) = - reac_rate_local(141) 
  reac_source_local(22,141) = + reac_rate_local(141) 
  reac_source_local(01,142) = - reac_rate_local(142) 
  reac_source_local(23,142) = + reac_rate_local(142) 
  reac_source_local(24,142) = + reac_rate_local(142) 
  reac_source_local(30,143) = - reac_rate_local(143) 
  reac_source_local(35,143) = + reac_rate_local(143) 
  reac_source_local(30,144) = - reac_rate_local(144) 
  reac_source_local(36,144) = + reac_rate_local(144) 
  reac_source_local(30,145) = - reac_rate_local(145) 
  reac_source_local(37,145) = + reac_rate_local(145) 
  reac_source_local(30,146) = - reac_rate_local(146) 
  reac_source_local(38,146) = + reac_rate_local(146) * 2.d0
  reac_source_local(30,147) = - reac_rate_local(147) 
  reac_source_local(38,147) = + reac_rate_local(147) 
  reac_source_local(39,147) = + reac_rate_local(147) 
  reac_source_local(30,148) = - reac_rate_local(148) 
  reac_source_local(38,148) = + reac_rate_local(148) 
  reac_source_local(40,148) = + reac_rate_local(148) 
  reac_source_local(35,149) = - reac_rate_local(149) 
  reac_source_local(38,149) = + reac_rate_local(149) * 2.d0
  reac_source_local(38,150) = - reac_rate_local(150) 
  reac_source_local(39,150) = + reac_rate_local(150) 
  reac_source_local(38,151) = - reac_rate_local(151) 
  reac_source_local(40,151) = + reac_rate_local(151) 
  reac_source_local(63,152) = + reac_rate_local(152) 
  reac_source_local(70,152) = - reac_rate_local(152) 
  reac_source_local(71,152) = + reac_rate_local(152) 
  reac_source_local(63,153) = + reac_rate_local(153) * 2.d0
  reac_source_local(70,153) = - reac_rate_local(153) 
  reac_source_local(72,153) = + reac_rate_local(153) 
  reac_source_local(23,154) = - reac_rate_local(154) 
  reac_source_local(26,154) = + reac_rate_local(154) 
  reac_source_local(79,154) = + reac_rate_local(154) 
  reac_source_local(38,155) = - reac_rate_local(155) 
  reac_source_local(42,155) = + reac_rate_local(155) 
  reac_source_local(79,155) = + reac_rate_local(155) 
  reac_source_local(01,156) = - reac_rate_local(156) 
  reac_source_local(27,156) = + reac_rate_local(156) 
  reac_source_local(79,156) = + reac_rate_local(156) 
  reac_source_local(17,157) = - reac_rate_local(157) 
  reac_source_local(27,157) = + reac_rate_local(157) 
  reac_source_local(79,157) = + reac_rate_local(157) 
  reac_source_local(30,158) = - reac_rate_local(158) 
  reac_source_local(43,158) = + reac_rate_local(158) 
  reac_source_local(79,158) = + reac_rate_local(158) 
  reac_source_local(35,159) = - reac_rate_local(159) 
  reac_source_local(43,159) = + reac_rate_local(159) 
  reac_source_local(79,159) = + reac_rate_local(159) 
  reac_source_local(50,160) = - reac_rate_local(160) 
  reac_source_local(55,160) = + reac_rate_local(160) 
  reac_source_local(79,160) = + reac_rate_local(160) 
  reac_source_local(51,161) = - reac_rate_local(161) 
  reac_source_local(56,161) = + reac_rate_local(161) 
  reac_source_local(79,161) = + reac_rate_local(161) 
  reac_source_local(41,162) = - reac_rate_local(162) 
  reac_source_local(44,162) = + reac_rate_local(162) 
  reac_source_local(79,162) = + reac_rate_local(162) 
  reac_source_local(70,163) = - reac_rate_local(163) 
  reac_source_local(73,163) = + reac_rate_local(163) 
  reac_source_local(79,163) = + reac_rate_local(163) 
  reac_source_local(64,164) = + reac_rate_local(164) 
  reac_source_local(70,164) = - reac_rate_local(164) 
  reac_source_local(71,164) = + reac_rate_local(164) 
  reac_source_local(79,164) = + reac_rate_local(164) 
  reac_source_local(65,165) = + reac_rate_local(165) 
  reac_source_local(70,165) = - reac_rate_local(165) 
  reac_source_local(72,165) = + reac_rate_local(165) 
  reac_source_local(79,165) = + reac_rate_local(165) 
  reac_source_local(23,166) = + reac_rate_local(166) * 2.d0
  reac_source_local(27,166) = - reac_rate_local(166) 
  reac_source_local(79,166) = - reac_rate_local(166) 
  reac_source_local(23,167) = + reac_rate_local(167) 
  reac_source_local(24,167) = + reac_rate_local(167) 
  reac_source_local(27,167) = - reac_rate_local(167) 
  reac_source_local(79,167) = - reac_rate_local(167) 
  reac_source_local(23,168) = + reac_rate_local(168) 
  reac_source_local(25,168) = + reac_rate_local(168) 
  reac_source_local(27,168) = - reac_rate_local(168) 
  reac_source_local(79,168) = - reac_rate_local(168) 
  reac_source_local(38,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(43,169) = - reac_rate_local(169) 
  reac_source_local(79,169) = - reac_rate_local(169) 
  reac_source_local(38,170) = + reac_rate_local(170) 
  reac_source_local(39,170) = + reac_rate_local(170) 
  reac_source_local(43,170) = - reac_rate_local(170) 
  reac_source_local(79,170) = - reac_rate_local(170) 
  reac_source_local(38,171) = + reac_rate_local(171) 
  reac_source_local(40,171) = + reac_rate_local(171) 
  reac_source_local(43,171) = - reac_rate_local(171) 
  reac_source_local(79,171) = - reac_rate_local(171) 
  reac_source_local(23,172) = + reac_rate_local(172) 
  reac_source_local(38,172) = + reac_rate_local(172) 
  reac_source_local(55,172) = - reac_rate_local(172) 
  reac_source_local(79,172) = - reac_rate_local(172) 
  reac_source_local(24,173) = + reac_rate_local(173) 
  reac_source_local(38,173) = + reac_rate_local(173) 
  reac_source_local(55,173) = - reac_rate_local(173) 
  reac_source_local(79,173) = - reac_rate_local(173) 
  reac_source_local(01,174) = + reac_rate_local(174) 
  reac_source_local(23,174) = + reac_rate_local(174) 
  reac_source_local(28,174) = - reac_rate_local(174) 
  reac_source_local(79,174) = - reac_rate_local(174) 
  reac_source_local(01,175) = + reac_rate_local(175) * 2.d0
  reac_source_local(29,175) = - reac_rate_local(175) 
  reac_source_local(79,175) = - reac_rate_local(175) 
  reac_source_local(01,176) = + reac_rate_local(176) 
  reac_source_local(38,176) = + reac_rate_local(176) 
  reac_source_local(56,176) = - reac_rate_local(176) 
  reac_source_local(79,176) = - reac_rate_local(176) 
  reac_source_local(38,177) = + reac_rate_local(177) 
  reac_source_local(50,177) = + reac_rate_local(177) 
  reac_source_local(57,177) = - reac_rate_local(177) 
  reac_source_local(79,177) = - reac_rate_local(177) 
  reac_source_local(30,178) = + reac_rate_local(178) * 2.d0
  reac_source_local(45,178) = - reac_rate_local(178) 
  reac_source_local(79,178) = - reac_rate_local(178) 
  reac_source_local(01,179) = + reac_rate_local(179) 
  reac_source_local(30,179) = + reac_rate_local(179) 
  reac_source_local(78,179) = - reac_rate_local(179) 
  reac_source_local(79,179) = - reac_rate_local(179) 
  reac_source_local(23,180) = + reac_rate_local(180) 
  reac_source_local(26,180) = - reac_rate_local(180) 
  reac_source_local(79,180) = - reac_rate_local(180) 
  reac_source_local(38,181) = + reac_rate_local(181) 
  reac_source_local(42,181) = - reac_rate_local(181) 
  reac_source_local(79,181) = - reac_rate_local(181) 
  reac_source_local(23,182) = + reac_rate_local(182) 
  reac_source_local(26,182) = - reac_rate_local(182) 
  reac_source_local(79,182) = - reac_rate_local(182) 
  reac_source_local(38,183) = + reac_rate_local(183) 
  reac_source_local(42,183) = - reac_rate_local(183) 
  reac_source_local(79,183) = - reac_rate_local(183) 
  reac_source_local(63,184) = + reac_rate_local(184) * 2.d0
  reac_source_local(71,184) = + reac_rate_local(184) 
  reac_source_local(77,184) = - reac_rate_local(184) 
  reac_source_local(79,184) = - reac_rate_local(184) 
  reac_source_local(30,185) = - reac_rate_local(185) 
  reac_source_local(38,185) = + reac_rate_local(185) 
  reac_source_local(46,185) = + reac_rate_local(185) 
  reac_source_local(79,185) = - reac_rate_local(185) 
  reac_source_local(23,186) = + reac_rate_local(186) 
  reac_source_local(46,186) = + reac_rate_local(186) 
  reac_source_local(50,186) = - reac_rate_local(186) 
  reac_source_local(79,186) = - reac_rate_local(186) 
  reac_source_local(30,187) = + reac_rate_local(187) 
  reac_source_local(41,187) = - reac_rate_local(187) 
  reac_source_local(46,187) = + reac_rate_local(187) 
  reac_source_local(79,187) = - reac_rate_local(187) 
  reac_source_local(38,188) = + reac_rate_local(188) 
  reac_source_local(41,188) = - reac_rate_local(188) 
  reac_source_local(47,188) = + reac_rate_local(188) 
  reac_source_local(79,188) = - reac_rate_local(188) 
  reac_source_local(46,189) = + reac_rate_local(189) 
  reac_source_local(50,189) = + reac_rate_local(189) 
  reac_source_local(52,189) = - reac_rate_local(189) 
  reac_source_local(79,189) = - reac_rate_local(189) 
  reac_source_local(38,190) = - reac_rate_local(190) 
  reac_source_local(46,190) = + reac_rate_local(190) 
  reac_source_local(79,190) = - reac_rate_local(190) 
  reac_source_local(30,191) = - reac_rate_local(191) 
  reac_source_local(47,191) = + reac_rate_local(191) 
  reac_source_local(79,191) = - reac_rate_local(191) 
  reac_source_local(41,192) = - reac_rate_local(192) 
  reac_source_local(48,192) = + reac_rate_local(192) 
  reac_source_local(79,192) = - reac_rate_local(192) 
  reac_source_local(50,193) = - reac_rate_local(193) 
  reac_source_local(58,193) = + reac_rate_local(193) 
  reac_source_local(79,193) = - reac_rate_local(193) 
  reac_source_local(51,194) = - reac_rate_local(194) 
  reac_source_local(59,194) = + reac_rate_local(194) 
  reac_source_local(79,194) = - reac_rate_local(194) 
  reac_source_local(30,195) = - reac_rate_local(195) 
  reac_source_local(47,195) = + reac_rate_local(195) 
  reac_source_local(79,195) = - reac_rate_local(195) 
  reac_source_local(67,196) = + reac_rate_local(196) 
  reac_source_local(70,196) = - reac_rate_local(196) 
  reac_source_local(71,196) = + reac_rate_local(196) 
  reac_source_local(79,196) = - reac_rate_local(196) 
  reac_source_local(30,197) = + reac_rate_local(197) 
  reac_source_local(38,197) = - reac_rate_local(197) 
  reac_source_local(46,197) = - reac_rate_local(197) 
  reac_source_local(79,197) = + reac_rate_local(197) 
  reac_source_local(23,198) = - reac_rate_local(198) 
  reac_source_local(46,198) = - reac_rate_local(198) 
  reac_source_local(50,198) = + reac_rate_local(198) 
  reac_source_local(79,198) = + reac_rate_local(198) 
  reac_source_local(46,199) = - reac_rate_local(199) 
  reac_source_local(50,199) = - reac_rate_local(199) 
  reac_source_local(52,199) = + reac_rate_local(199) 
  reac_source_local(79,199) = + reac_rate_local(199) 
  reac_source_local(01,200) = - reac_rate_local(200) 
  reac_source_local(46,200) = - reac_rate_local(200) 
  reac_source_local(51,200) = + reac_rate_local(200) 
  reac_source_local(79,200) = + reac_rate_local(200) 
  reac_source_local(30,201) = - reac_rate_local(201) 
  reac_source_local(41,201) = + reac_rate_local(201) 
  reac_source_local(46,201) = - reac_rate_local(201) 
  reac_source_local(79,201) = + reac_rate_local(201) 
  reac_source_local(35,202) = - reac_rate_local(202) 
  reac_source_local(41,202) = + reac_rate_local(202) 
  reac_source_local(46,202) = - reac_rate_local(202) 
  reac_source_local(79,202) = + reac_rate_local(202) 
  reac_source_local(30,203) = + reac_rate_local(203) 
  reac_source_local(36,203) = - reac_rate_local(203) 
  reac_source_local(38,203) = + reac_rate_local(203) 
  reac_source_local(46,203) = - reac_rate_local(203) 
  reac_source_local(79,203) = + reac_rate_local(203) 
  reac_source_local(01,204) = + reac_rate_local(204) 
  reac_source_local(17,204) = - reac_rate_local(204) 
  reac_source_local(38,204) = + reac_rate_local(204) 
  reac_source_local(46,204) = - reac_rate_local(204) 
  reac_source_local(79,204) = + reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(18,205) = - reac_rate_local(205) 
  reac_source_local(38,205) = + reac_rate_local(205) 
  reac_source_local(46,205) = - reac_rate_local(205) 
  reac_source_local(79,205) = + reac_rate_local(205) 
  reac_source_local(30,206) = + reac_rate_local(206) * 2.d0
  reac_source_local(41,206) = - reac_rate_local(206) 
  reac_source_local(46,206) = - reac_rate_local(206) 
  reac_source_local(79,206) = + reac_rate_local(206) 
  reac_source_local(38,207) = - reac_rate_local(207) 
  reac_source_local(41,207) = + reac_rate_local(207) 
  reac_source_local(47,207) = - reac_rate_local(207) 
  reac_source_local(79,207) = + reac_rate_local(207) 
  reac_source_local(23,208) = - reac_rate_local(208) 
  reac_source_local(47,208) = - reac_rate_local(208) 
  reac_source_local(52,208) = + reac_rate_local(208) 
  reac_source_local(79,208) = + reac_rate_local(208) 
  reac_source_local(30,209) = + reac_rate_local(209) 
  reac_source_local(47,209) = - reac_rate_local(209) 
  reac_source_local(79,209) = + reac_rate_local(209) 
  reac_source_local(30,210) = + reac_rate_local(210) * 2.d0
  reac_source_local(35,210) = - reac_rate_local(210) 
  reac_source_local(47,210) = - reac_rate_local(210) 
  reac_source_local(79,210) = + reac_rate_local(210) 
  reac_source_local(30,211) = + reac_rate_local(211) * 2.d0
  reac_source_local(36,211) = - reac_rate_local(211) 
  reac_source_local(47,211) = - reac_rate_local(211) 
  reac_source_local(79,211) = + reac_rate_local(211) 
  reac_source_local(30,212) = + reac_rate_local(212) 
  reac_source_local(47,212) = - reac_rate_local(212) 
  reac_source_local(79,212) = + reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) 
  reac_source_local(17,213) = - reac_rate_local(213) 
  reac_source_local(30,213) = + reac_rate_local(213) 
  reac_source_local(47,213) = - reac_rate_local(213) 
  reac_source_local(79,213) = + reac_rate_local(213) 
  reac_source_local(01,214) = + reac_rate_local(214) 
  reac_source_local(18,214) = - reac_rate_local(214) 
  reac_source_local(30,214) = + reac_rate_local(214) 
  reac_source_local(47,214) = - reac_rate_local(214) 
  reac_source_local(79,214) = + reac_rate_local(214) 
  reac_source_local(30,215) = + reac_rate_local(215) * 2.d0
  reac_source_local(38,215) = - reac_rate_local(215) 
  reac_source_local(48,215) = - reac_rate_local(215) 
  reac_source_local(79,215) = + reac_rate_local(215) 
  reac_source_local(23,216) = - reac_rate_local(216) 
  reac_source_local(51,216) = + reac_rate_local(216) 
  reac_source_local(58,216) = - reac_rate_local(216) 
  reac_source_local(79,216) = + reac_rate_local(216) 
  reac_source_local(23,217) = - reac_rate_local(217) 
  reac_source_local(30,217) = + reac_rate_local(217) 
  reac_source_local(48,217) = - reac_rate_local(217) 
  reac_source_local(50,217) = + reac_rate_local(217) 
  reac_source_local(79,217) = + reac_rate_local(217) 
  reac_source_local(01,218) = + reac_rate_local(218) 
  reac_source_local(23,218) = - reac_rate_local(218) 
  reac_source_local(50,218) = + reac_rate_local(218) 
  reac_source_local(59,218) = - reac_rate_local(218) 
  reac_source_local(79,218) = + reac_rate_local(218) 
  reac_source_local(23,219) = - reac_rate_local(219) 
  reac_source_local(50,219) = + reac_rate_local(219) * 2.d0
  reac_source_local(60,219) = - reac_rate_local(219) 
  reac_source_local(79,219) = + reac_rate_local(219) 
  reac_source_local(23,220) = - reac_rate_local(220) 
  reac_source_local(50,220) = + reac_rate_local(220) 
  reac_source_local(52,220) = + reac_rate_local(220) 
  reac_source_local(61,220) = - reac_rate_local(220) 
  reac_source_local(79,220) = + reac_rate_local(220) 
  reac_source_local(38,221) = - reac_rate_local(221) 
  reac_source_local(52,221) = + reac_rate_local(221) 
  reac_source_local(58,221) = - reac_rate_local(221) 
  reac_source_local(79,221) = + reac_rate_local(221) 
  reac_source_local(38,222) = - reac_rate_local(222) 
  reac_source_local(50,222) = + reac_rate_local(222) * 2.d0
  reac_source_local(59,222) = - reac_rate_local(222) 
  reac_source_local(79,222) = + reac_rate_local(222) 
  reac_source_local(30,223) = + reac_rate_local(223) 
  reac_source_local(38,223) = - reac_rate_local(223) 
  reac_source_local(50,223) = + reac_rate_local(223) 
  reac_source_local(60,223) = - reac_rate_local(223) 
  reac_source_local(79,223) = + reac_rate_local(223) 
  reac_source_local(38,224) = - reac_rate_local(224) 
  reac_source_local(41,224) = + reac_rate_local(224) 
  reac_source_local(50,224) = + reac_rate_local(224) 
  reac_source_local(61,224) = - reac_rate_local(224) 
  reac_source_local(79,224) = + reac_rate_local(224) 
  reac_source_local(01,225) = + reac_rate_local(225) 
  reac_source_local(17,225) = - reac_rate_local(225) 
  reac_source_local(41,225) = + reac_rate_local(225) 
  reac_source_local(48,225) = - reac_rate_local(225) 
  reac_source_local(79,225) = + reac_rate_local(225) 
  reac_source_local(01,226) = + reac_rate_local(226) 
  reac_source_local(17,226) = - reac_rate_local(226) 
  reac_source_local(50,226) = + reac_rate_local(226) 
  reac_source_local(58,226) = - reac_rate_local(226) 
  reac_source_local(79,226) = + reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(17,227) = - reac_rate_local(227) 
  reac_source_local(51,227) = + reac_rate_local(227) 
  reac_source_local(59,227) = - reac_rate_local(227) 
  reac_source_local(79,227) = + reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(17,228) = - reac_rate_local(228) 
  reac_source_local(52,228) = + reac_rate_local(228) 
  reac_source_local(60,228) = - reac_rate_local(228) 
  reac_source_local(79,228) = + reac_rate_local(228) 
  reac_source_local(01,229) = + reac_rate_local(229) 
  reac_source_local(17,229) = - reac_rate_local(229) 
  reac_source_local(53,229) = + reac_rate_local(229) 
  reac_source_local(61,229) = - reac_rate_local(229) 
  reac_source_local(79,229) = + reac_rate_local(229) 
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(18,230) = - reac_rate_local(230) 
  reac_source_local(41,230) = + reac_rate_local(230) 
  reac_source_local(48,230) = - reac_rate_local(230) 
  reac_source_local(79,230) = + reac_rate_local(230) 
  reac_source_local(01,231) = + reac_rate_local(231) 
  reac_source_local(18,231) = - reac_rate_local(231) 
  reac_source_local(50,231) = + reac_rate_local(231) 
  reac_source_local(58,231) = - reac_rate_local(231) 
  reac_source_local(79,231) = + reac_rate_local(231) 
  reac_source_local(01,232) = + reac_rate_local(232) 
  reac_source_local(18,232) = - reac_rate_local(232) 
  reac_source_local(51,232) = + reac_rate_local(232) 
  reac_source_local(59,232) = - reac_rate_local(232) 
  reac_source_local(79,232) = + reac_rate_local(232) 
  reac_source_local(01,233) = + reac_rate_local(233) 
  reac_source_local(18,233) = - reac_rate_local(233) 
  reac_source_local(52,233) = + reac_rate_local(233) 
  reac_source_local(60,233) = - reac_rate_local(233) 
  reac_source_local(79,233) = + reac_rate_local(233) 
  reac_source_local(01,234) = + reac_rate_local(234) 
  reac_source_local(18,234) = - reac_rate_local(234) 
  reac_source_local(53,234) = + reac_rate_local(234) 
  reac_source_local(61,234) = - reac_rate_local(234) 
  reac_source_local(79,234) = + reac_rate_local(234) 
  reac_source_local(01,235) = + reac_rate_local(235) 
  reac_source_local(17,235) = - reac_rate_local(235) 
  reac_source_local(17,236) = + reac_rate_local(236) 
  reac_source_local(18,236) = - reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(19,237) = - reac_rate_local(237) 
  reac_source_local(18,238) = + reac_rate_local(238) 
  reac_source_local(20,238) = - reac_rate_local(238) 
  reac_source_local(30,239) = + reac_rate_local(239) 
  reac_source_local(35,239) = - reac_rate_local(239) 
  reac_source_local(35,240) = + reac_rate_local(240) 
  reac_source_local(36,240) = - reac_rate_local(240) 
  reac_source_local(30,241) = + reac_rate_local(241) 
  reac_source_local(36,241) = - reac_rate_local(241) 
  reac_source_local(30,242) = + reac_rate_local(242) 
  reac_source_local(37,242) = - reac_rate_local(242) 
  reac_source_local(01,243) = + reac_rate_local(243) 
  reac_source_local(17,243) = - reac_rate_local(243) 
  reac_source_local(30,243) = - reac_rate_local(243) 
  reac_source_local(38,243) = + reac_rate_local(243) * 2.d0
  reac_source_local(01,244) = + reac_rate_local(244) 
  reac_source_local(17,244) = - reac_rate_local(244) 
  reac_source_local(30,244) = - reac_rate_local(244) 
  reac_source_local(36,244) = + reac_rate_local(244) 
  reac_source_local(01,245) = + reac_rate_local(245) 
  reac_source_local(18,245) = - reac_rate_local(245) 
  reac_source_local(30,245) = - reac_rate_local(245) 
  reac_source_local(38,245) = + reac_rate_local(245) * 2.d0
  reac_source_local(01,246) = + reac_rate_local(246) 
  reac_source_local(19,246) = - reac_rate_local(246) 
  reac_source_local(30,246) = - reac_rate_local(246) 
  reac_source_local(38,246) = + reac_rate_local(246) 
  reac_source_local(39,246) = + reac_rate_local(246) 
  reac_source_local(01,247) = + reac_rate_local(247) 
  reac_source_local(20,247) = - reac_rate_local(247) 
  reac_source_local(30,247) = - reac_rate_local(247) 
  reac_source_local(38,247) = + reac_rate_local(247) * 2.d0
  reac_source_local(01,248) = + reac_rate_local(248) 
  reac_source_local(20,248) = - reac_rate_local(248) 
  reac_source_local(30,248) = - reac_rate_local(248) 
  reac_source_local(38,248) = + reac_rate_local(248) 
  reac_source_local(39,248) = + reac_rate_local(248) 
  reac_source_local(01,249) = + reac_rate_local(249) 
  reac_source_local(20,249) = - reac_rate_local(249) 
  reac_source_local(30,249) = - reac_rate_local(249) 
  reac_source_local(38,249) = + reac_rate_local(249) 
  reac_source_local(40,249) = + reac_rate_local(249) 
  reac_source_local(17,250) = - reac_rate_local(250) 
  reac_source_local(24,250) = + reac_rate_local(250) 
  reac_source_local(38,250) = - reac_rate_local(250) 
  reac_source_local(50,250) = + reac_rate_local(250) 
  reac_source_local(01,251) = + reac_rate_local(251) 
  reac_source_local(17,251) = - reac_rate_local(251) 
  reac_source_local(38,251) = - reac_rate_local(251) 
  reac_source_local(40,251) = + reac_rate_local(251) 
  reac_source_local(01,252) = + reac_rate_local(252) 
  reac_source_local(17,252) = - reac_rate_local(252) 
  reac_source_local(01,253) = + reac_rate_local(253) 
  reac_source_local(17,253) = - reac_rate_local(253) 
  reac_source_local(23,253) = - reac_rate_local(253) 
  reac_source_local(25,253) = + reac_rate_local(253) 
  reac_source_local(01,254) = + reac_rate_local(254) 
  reac_source_local(17,254) = - reac_rate_local(254) 
  reac_source_local(01,255) = + reac_rate_local(255) 
  reac_source_local(17,255) = - reac_rate_local(255) 
  reac_source_local(01,256) = + reac_rate_local(256) 
  reac_source_local(17,256) = - reac_rate_local(256) 
  reac_source_local(23,256) = + reac_rate_local(256) 
  reac_source_local(50,256) = + reac_rate_local(256) 
  reac_source_local(51,256) = - reac_rate_local(256) 
  reac_source_local(01,257) = + reac_rate_local(257) 
  reac_source_local(17,257) = - reac_rate_local(257) 
  reac_source_local(38,257) = + reac_rate_local(257) 
  reac_source_local(50,257) = + reac_rate_local(257) 
  reac_source_local(52,257) = - reac_rate_local(257) 
  reac_source_local(01,258) = + reac_rate_local(258) 
  reac_source_local(17,258) = - reac_rate_local(258) * 2.d0
  reac_source_local(18,258) = + reac_rate_local(258) 
  reac_source_local(01,259) = + reac_rate_local(259) 
  reac_source_local(17,259) = - reac_rate_local(259) * 2.d0
  reac_source_local(20,259) = + reac_rate_local(259) 
  reac_source_local(17,260) = + reac_rate_local(260) 
  reac_source_local(18,260) = - reac_rate_local(260) 
  reac_source_local(01,261) = + reac_rate_local(261) 
  reac_source_local(18,261) = - reac_rate_local(261) 
  reac_source_local(17,262) = + reac_rate_local(262) 
  reac_source_local(18,262) = - reac_rate_local(262) 
  reac_source_local(19,263) = + reac_rate_local(263) 
  reac_source_local(20,263) = - reac_rate_local(263) 
  reac_source_local(18,264) = + reac_rate_local(264) 
  reac_source_local(19,264) = - reac_rate_local(264) 
  reac_source_local(01,265) = + reac_rate_local(265) 
  reac_source_local(19,265) = - reac_rate_local(265) 
  reac_source_local(23,265) = + reac_rate_local(265) 
  reac_source_local(38,265) = + reac_rate_local(265) 
  reac_source_local(50,265) = - reac_rate_local(265) 
  reac_source_local(17,266) = - reac_rate_local(266) 
  reac_source_local(19,266) = - reac_rate_local(266) 
  reac_source_local(29,266) = + reac_rate_local(266) 
  reac_source_local(79,266) = + reac_rate_local(266) 
  reac_source_local(19,267) = - reac_rate_local(267) * 2.d0
  reac_source_local(29,267) = + reac_rate_local(267) 
  reac_source_local(79,267) = + reac_rate_local(267) 
  reac_source_local(01,268) = + reac_rate_local(268) 
  reac_source_local(17,268) = - reac_rate_local(268) 
  reac_source_local(01,269) = + reac_rate_local(269) 
  reac_source_local(17,269) = - reac_rate_local(269) 
  reac_source_local(62,269) = - reac_rate_local(269) 
  reac_source_local(63,269) = + reac_rate_local(269) * 2.d0
  reac_source_local(01,270) = + reac_rate_local(270) 
  reac_source_local(17,270) = - reac_rate_local(270) 
  reac_source_local(17,271) = + reac_rate_local(271) 
  reac_source_local(18,271) = - reac_rate_local(271) 
  reac_source_local(01,272) = + reac_rate_local(272) 
  reac_source_local(19,272) = - reac_rate_local(272) 
  reac_source_local(01,273) = + reac_rate_local(273) 
  reac_source_local(19,273) = - reac_rate_local(273) 
  reac_source_local(17,274) = + reac_rate_local(274) 
  reac_source_local(23,274) = - reac_rate_local(274) * 2.d0
  reac_source_local(17,275) = + reac_rate_local(275) 
  reac_source_local(23,275) = - reac_rate_local(275) * 2.d0
  reac_source_local(17,276) = + reac_rate_local(276) 
  reac_source_local(23,276) = - reac_rate_local(276) * 2.d0
  reac_source_local(17,277) = + reac_rate_local(277) 
  reac_source_local(23,277) = - reac_rate_local(277) * 2.d0
  reac_source_local(17,278) = + reac_rate_local(278) 
  reac_source_local(23,278) = - reac_rate_local(278) * 2.d0
  reac_source_local(18,279) = + reac_rate_local(279) 
  reac_source_local(23,279) = - reac_rate_local(279) * 2.d0
  reac_source_local(18,280) = + reac_rate_local(280) 
  reac_source_local(23,280) = - reac_rate_local(280) * 2.d0
  reac_source_local(18,281) = + reac_rate_local(281) 
  reac_source_local(23,281) = - reac_rate_local(281) * 2.d0
  reac_source_local(18,282) = + reac_rate_local(282) 
  reac_source_local(23,282) = - reac_rate_local(282) * 2.d0
  reac_source_local(18,283) = + reac_rate_local(283) 
  reac_source_local(23,283) = - reac_rate_local(283) * 2.d0
  reac_source_local(23,284) = + reac_rate_local(284) 
  reac_source_local(24,284) = - reac_rate_local(284) 
  reac_source_local(38,284) = - reac_rate_local(284) 
  reac_source_local(39,284) = + reac_rate_local(284) 
  reac_source_local(24,285) = - reac_rate_local(285) 
  reac_source_local(30,285) = - reac_rate_local(285) 
  reac_source_local(38,285) = + reac_rate_local(285) 
  reac_source_local(50,285) = + reac_rate_local(285) 
  reac_source_local(01,286) = + reac_rate_local(286) 
  reac_source_local(24,286) = - reac_rate_local(286) 
  reac_source_local(38,286) = + reac_rate_local(286) 
  reac_source_local(50,286) = - reac_rate_local(286) 
  reac_source_local(01,287) = + reac_rate_local(287) 
  reac_source_local(24,287) = - reac_rate_local(287) 
  reac_source_local(50,287) = + reac_rate_local(287) 
  reac_source_local(51,287) = - reac_rate_local(287) 
  reac_source_local(23,288) = + reac_rate_local(288) 
  reac_source_local(24,288) = - reac_rate_local(288) 
  reac_source_local(23,289) = + reac_rate_local(289) 
  reac_source_local(25,289) = - reac_rate_local(289) 
  reac_source_local(23,290) = + reac_rate_local(290) 
  reac_source_local(25,290) = - reac_rate_local(290) 
  reac_source_local(24,291) = + reac_rate_local(291) 
  reac_source_local(25,291) = - reac_rate_local(291) 
  reac_source_local(23,292) = + reac_rate_local(292) 
  reac_source_local(25,292) = - reac_rate_local(292) 
  reac_source_local(24,293) = - reac_rate_local(293) 
  reac_source_local(25,293) = - reac_rate_local(293) 
  reac_source_local(27,293) = + reac_rate_local(293) 
  reac_source_local(79,293) = + reac_rate_local(293) 
  reac_source_local(25,294) = - reac_rate_local(294) 
  reac_source_local(30,294) = - reac_rate_local(294) 
  reac_source_local(38,294) = + reac_rate_local(294) 
  reac_source_local(50,294) = + reac_rate_local(294) 
  reac_source_local(17,295) = + reac_rate_local(295) 
  reac_source_local(25,295) = - reac_rate_local(295) 
  reac_source_local(38,295) = + reac_rate_local(295) 
  reac_source_local(50,295) = - reac_rate_local(295) 
  reac_source_local(24,296) = - reac_rate_local(296) 
  reac_source_local(62,296) = - reac_rate_local(296) 
  reac_source_local(63,296) = + reac_rate_local(296) 
  reac_source_local(72,296) = + reac_rate_local(296) 
  reac_source_local(24,297) = - reac_rate_local(297) 
  reac_source_local(70,297) = - reac_rate_local(297) 
  reac_source_local(71,297) = + reac_rate_local(297) 
  reac_source_local(72,297) = + reac_rate_local(297) 
  reac_source_local(25,298) = - reac_rate_local(298) 
  reac_source_local(62,298) = - reac_rate_local(298) 
  reac_source_local(63,298) = + reac_rate_local(298) 
  reac_source_local(72,298) = + reac_rate_local(298) 
  reac_source_local(30,299) = + reac_rate_local(299) 
  reac_source_local(35,299) = - reac_rate_local(299) 
  reac_source_local(23,300) = - reac_rate_local(300) 
  reac_source_local(35,300) = - reac_rate_local(300) 
  reac_source_local(38,300) = + reac_rate_local(300) 
  reac_source_local(50,300) = + reac_rate_local(300) 
  reac_source_local(30,301) = + reac_rate_local(301) 
  reac_source_local(35,301) = - reac_rate_local(301) 
  reac_source_local(30,302) = + reac_rate_local(302) 
  reac_source_local(35,302) = - reac_rate_local(302) 
  reac_source_local(30,303) = + reac_rate_local(303) 
  reac_source_local(35,303) = - reac_rate_local(303) 
  reac_source_local(30,304) = + reac_rate_local(304) * 2.d0
  reac_source_local(35,304) = - reac_rate_local(304) 
  reac_source_local(39,304) = + reac_rate_local(304) 
  reac_source_local(41,304) = - reac_rate_local(304) 
  reac_source_local(30,305) = + reac_rate_local(305) 
  reac_source_local(35,305) = - reac_rate_local(305) * 2.d0
  reac_source_local(36,305) = + reac_rate_local(305) 
  reac_source_local(30,306) = + reac_rate_local(306) 
  reac_source_local(35,306) = + reac_rate_local(306) 
  reac_source_local(38,306) = - reac_rate_local(306) 
  reac_source_local(41,306) = - reac_rate_local(306) 
  reac_source_local(35,307) = + reac_rate_local(307) 
  reac_source_local(36,307) = - reac_rate_local(307) 
  reac_source_local(30,308) = + reac_rate_local(308) 
  reac_source_local(36,308) = - reac_rate_local(308) 
  reac_source_local(38,308) = - reac_rate_local(308) 
  reac_source_local(39,308) = + reac_rate_local(308) 
  reac_source_local(35,309) = + reac_rate_local(309) 
  reac_source_local(36,309) = - reac_rate_local(309) 
  reac_source_local(35,310) = + reac_rate_local(310) 
  reac_source_local(36,310) = - reac_rate_local(310) 
  reac_source_local(35,311) = + reac_rate_local(311) 
  reac_source_local(36,311) = - reac_rate_local(311) 
  reac_source_local(30,312) = + reac_rate_local(312) * 2.d0
  reac_source_local(36,312) = - reac_rate_local(312) 
  reac_source_local(38,312) = + reac_rate_local(312) 
  reac_source_local(41,312) = - reac_rate_local(312) 
  reac_source_local(30,313) = + reac_rate_local(313) 
  reac_source_local(37,313) = - reac_rate_local(313) 
  reac_source_local(38,313) = - reac_rate_local(313) 
  reac_source_local(40,313) = + reac_rate_local(313) 
  reac_source_local(30,314) = - reac_rate_local(314) 
  reac_source_local(36,314) = + reac_rate_local(314) * 2.d0
  reac_source_local(37,314) = - reac_rate_local(314) 
  reac_source_local(36,315) = + reac_rate_local(315) 
  reac_source_local(37,315) = - reac_rate_local(315) 
  reac_source_local(38,316) = + reac_rate_local(316) 
  reac_source_local(39,316) = - reac_rate_local(316) 
  reac_source_local(38,317) = + reac_rate_local(317) 
  reac_source_local(39,317) = - reac_rate_local(317) 
  reac_source_local(30,318) = - reac_rate_local(318) 
  reac_source_local(35,318) = + reac_rate_local(318) 
  reac_source_local(38,318) = + reac_rate_local(318) 
  reac_source_local(39,318) = - reac_rate_local(318) 
  reac_source_local(30,319) = - reac_rate_local(319) 
  reac_source_local(36,319) = + reac_rate_local(319) 
  reac_source_local(38,319) = + reac_rate_local(319) 
  reac_source_local(39,319) = - reac_rate_local(319) 
  reac_source_local(38,320) = + reac_rate_local(320) 
  reac_source_local(39,320) = - reac_rate_local(320) 
  reac_source_local(30,321) = + reac_rate_local(321) 
  reac_source_local(38,321) = + reac_rate_local(321) * 2.d0
  reac_source_local(39,321) = - reac_rate_local(321) 
  reac_source_local(41,321) = - reac_rate_local(321) 
  reac_source_local(30,322) = + reac_rate_local(322) * 2.d0
  reac_source_local(39,322) = - reac_rate_local(322) 
  reac_source_local(41,322) = - reac_rate_local(322) 
  reac_source_local(23,323) = + reac_rate_local(323) 
  reac_source_local(30,323) = + reac_rate_local(323) 
  reac_source_local(39,323) = - reac_rate_local(323) 
  reac_source_local(50,323) = - reac_rate_local(323) 
  reac_source_local(39,324) = - reac_rate_local(324) 
  reac_source_local(50,324) = + reac_rate_local(324) * 2.d0
  reac_source_local(51,324) = - reac_rate_local(324) 
  reac_source_local(01,325) = + reac_rate_local(325) 
  reac_source_local(30,325) = + reac_rate_local(325) 
  reac_source_local(39,325) = - reac_rate_local(325) 
  reac_source_local(51,325) = - reac_rate_local(325) 
  reac_source_local(39,326) = + reac_rate_local(326) 
  reac_source_local(40,326) = - reac_rate_local(326) 
  reac_source_local(38,327) = + reac_rate_local(327) 
  reac_source_local(40,327) = - reac_rate_local(327) 
  reac_source_local(39,328) = + reac_rate_local(328) 
  reac_source_local(40,328) = - reac_rate_local(328) 
  reac_source_local(30,329) = - reac_rate_local(329) 
  reac_source_local(38,329) = + reac_rate_local(329) * 3.d0
  reac_source_local(40,329) = - reac_rate_local(329) 
  reac_source_local(38,330) = + reac_rate_local(330) 
  reac_source_local(40,330) = - reac_rate_local(330) 
  reac_source_local(35,331) = - reac_rate_local(331) 
  reac_source_local(37,331) = + reac_rate_local(331) 
  reac_source_local(38,331) = + reac_rate_local(331) 
  reac_source_local(40,331) = - reac_rate_local(331) 
  reac_source_local(35,332) = - reac_rate_local(332) 
  reac_source_local(36,332) = + reac_rate_local(332) 
  reac_source_local(39,332) = + reac_rate_local(332) 
  reac_source_local(40,332) = - reac_rate_local(332) 
  reac_source_local(35,333) = - reac_rate_local(333) 
  reac_source_local(38,333) = + reac_rate_local(333) * 3.d0
  reac_source_local(40,333) = - reac_rate_local(333) 
  reac_source_local(38,334) = + reac_rate_local(334) 
  reac_source_local(40,334) = - reac_rate_local(334) 
  reac_source_local(39,335) = + reac_rate_local(335) 
  reac_source_local(40,335) = - reac_rate_local(335) 
  reac_source_local(30,336) = + reac_rate_local(336) * 2.d0
  reac_source_local(40,336) = - reac_rate_local(336) 
  reac_source_local(41,336) = - reac_rate_local(336) 
  reac_source_local(30,337) = + reac_rate_local(337) 
  reac_source_local(38,337) = + reac_rate_local(337) 
  reac_source_local(39,337) = + reac_rate_local(337) 
  reac_source_local(40,337) = - reac_rate_local(337) 
  reac_source_local(41,337) = - reac_rate_local(337) 
  reac_source_local(38,338) = + reac_rate_local(338) 
  reac_source_local(40,338) = - reac_rate_local(338) 
  reac_source_local(39,339) = + reac_rate_local(339) 
  reac_source_local(40,339) = - reac_rate_local(339) 
  reac_source_local(01,340) = + reac_rate_local(340) 
  reac_source_local(23,340) = - reac_rate_local(340) 
  reac_source_local(38,340) = + reac_rate_local(340) 
  reac_source_local(50,340) = - reac_rate_local(340) 
  reac_source_local(23,341) = - reac_rate_local(341) 
  reac_source_local(30,341) = - reac_rate_local(341) 
  reac_source_local(38,341) = + reac_rate_local(341) 
  reac_source_local(50,341) = + reac_rate_local(341) 
  reac_source_local(01,342) = + reac_rate_local(342) 
  reac_source_local(23,342) = - reac_rate_local(342) 
  reac_source_local(38,342) = + reac_rate_local(342) * 2.d0
  reac_source_local(52,342) = - reac_rate_local(342) 
  reac_source_local(23,343) = - reac_rate_local(343) 
  reac_source_local(38,343) = + reac_rate_local(343) 
  reac_source_local(51,343) = + reac_rate_local(343) 
  reac_source_local(52,343) = - reac_rate_local(343) 
  reac_source_local(01,344) = + reac_rate_local(344) 
  reac_source_local(23,344) = - reac_rate_local(344) 
  reac_source_local(30,344) = + reac_rate_local(344) 
  reac_source_local(52,344) = - reac_rate_local(344) 
  reac_source_local(23,345) = - reac_rate_local(345) 
  reac_source_local(50,345) = + reac_rate_local(345) * 2.d0
  reac_source_local(52,345) = - reac_rate_local(345) 
  reac_source_local(01,346) = - reac_rate_local(346) 
  reac_source_local(23,346) = + reac_rate_local(346) 
  reac_source_local(38,346) = - reac_rate_local(346) 
  reac_source_local(50,346) = + reac_rate_local(346) 
  reac_source_local(23,347) = + reac_rate_local(347) 
  reac_source_local(30,347) = + reac_rate_local(347) 
  reac_source_local(38,347) = - reac_rate_local(347) 
  reac_source_local(50,347) = - reac_rate_local(347) 
  reac_source_local(38,348) = - reac_rate_local(348) 
  reac_source_local(50,348) = - reac_rate_local(348) 
  reac_source_local(52,348) = + reac_rate_local(348) 
  reac_source_local(01,349) = + reac_rate_local(349) 
  reac_source_local(30,349) = + reac_rate_local(349) 
  reac_source_local(38,349) = - reac_rate_local(349) 
  reac_source_local(51,349) = - reac_rate_local(349) 
  reac_source_local(38,350) = - reac_rate_local(350) 
  reac_source_local(50,350) = + reac_rate_local(350) * 2.d0
  reac_source_local(51,350) = - reac_rate_local(350) 
  reac_source_local(30,351) = + reac_rate_local(351) 
  reac_source_local(38,351) = - reac_rate_local(351) 
  reac_source_local(50,351) = + reac_rate_local(351) 
  reac_source_local(52,351) = - reac_rate_local(351) 
  reac_source_local(30,352) = + reac_rate_local(352) 
  reac_source_local(38,352) = - reac_rate_local(352) 
  reac_source_local(52,352) = + reac_rate_local(352) 
  reac_source_local(53,352) = - reac_rate_local(352) 
  reac_source_local(01,353) = - reac_rate_local(353) 
  reac_source_local(30,353) = - reac_rate_local(353) 
  reac_source_local(38,353) = + reac_rate_local(353) 
  reac_source_local(51,353) = + reac_rate_local(353) 
  reac_source_local(23,354) = + reac_rate_local(354) 
  reac_source_local(50,354) = - reac_rate_local(354) * 2.d0
  reac_source_local(52,354) = + reac_rate_local(354) 
  reac_source_local(38,355) = + reac_rate_local(355) 
  reac_source_local(50,355) = - reac_rate_local(355) * 2.d0
  reac_source_local(51,355) = + reac_rate_local(355) 
  reac_source_local(01,356) = + reac_rate_local(356) 
  reac_source_local(30,356) = + reac_rate_local(356) 
  reac_source_local(50,356) = - reac_rate_local(356) * 2.d0
  reac_source_local(30,357) = - reac_rate_local(357) 
  reac_source_local(38,357) = + reac_rate_local(357) 
  reac_source_local(50,357) = - reac_rate_local(357) 
  reac_source_local(52,357) = + reac_rate_local(357) 
  reac_source_local(30,358) = + reac_rate_local(358) 
  reac_source_local(41,358) = - reac_rate_local(358) 
  reac_source_local(50,358) = - reac_rate_local(358) 
  reac_source_local(52,358) = + reac_rate_local(358) 
  reac_source_local(01,359) = + reac_rate_local(359) 
  reac_source_local(50,359) = - reac_rate_local(359) 
  reac_source_local(51,359) = - reac_rate_local(359) 
  reac_source_local(52,359) = + reac_rate_local(359) 
  reac_source_local(50,360) = - reac_rate_local(360) 
  reac_source_local(52,360) = + reac_rate_local(360) * 2.d0
  reac_source_local(53,360) = - reac_rate_local(360) 
  reac_source_local(30,361) = - reac_rate_local(361) * 2.d0
  reac_source_local(38,361) = + reac_rate_local(361) 
  reac_source_local(41,361) = + reac_rate_local(361) 
  reac_source_local(30,362) = - reac_rate_local(362) 
  reac_source_local(41,362) = + reac_rate_local(362) 
  reac_source_local(50,362) = + reac_rate_local(362) 
  reac_source_local(52,362) = - reac_rate_local(362) 
  reac_source_local(30,363) = + reac_rate_local(363) 
  reac_source_local(50,363) = + reac_rate_local(363) * 2.d0
  reac_source_local(52,363) = - reac_rate_local(363) * 2.d0
  reac_source_local(50,364) = + reac_rate_local(364) 
  reac_source_local(52,364) = - reac_rate_local(364) * 2.d0
  reac_source_local(53,364) = + reac_rate_local(364) 
  reac_source_local(30,365) = + reac_rate_local(365) 
  reac_source_local(41,365) = - reac_rate_local(365) 
  reac_source_local(52,365) = - reac_rate_local(365) 
  reac_source_local(53,365) = + reac_rate_local(365) 
  reac_source_local(30,366) = + reac_rate_local(366) 
  reac_source_local(50,366) = + reac_rate_local(366) 
  reac_source_local(53,366) = - reac_rate_local(366) 
  reac_source_local(30,367) = - reac_rate_local(367) 
  reac_source_local(41,367) = + reac_rate_local(367) 
  reac_source_local(52,367) = + reac_rate_local(367) 
  reac_source_local(53,367) = - reac_rate_local(367) 
  reac_source_local(30,368) = + reac_rate_local(368) 
  reac_source_local(52,368) = + reac_rate_local(368) * 2.d0
  reac_source_local(53,368) = - reac_rate_local(368) * 2.d0
  reac_source_local(23,369) = - reac_rate_local(369) * 2.d0
  reac_source_local(27,369) = + reac_rate_local(369) 
  reac_source_local(79,369) = + reac_rate_local(369) 
  reac_source_local(23,370) = - reac_rate_local(370) 
  reac_source_local(38,370) = - reac_rate_local(370) 
  reac_source_local(55,370) = + reac_rate_local(370) 
  reac_source_local(79,370) = + reac_rate_local(370) 
  reac_source_local(23,371) = + reac_rate_local(371) 
  reac_source_local(62,371) = + reac_rate_local(371) 
  reac_source_local(63,371) = - reac_rate_local(371) 
  reac_source_local(72,371) = - reac_rate_local(371) 
  reac_source_local(62,372) = + reac_rate_local(372) 
  reac_source_local(63,372) = - reac_rate_local(372) 
  reac_source_local(71,372) = - reac_rate_local(372) 
  reac_source_local(72,372) = + reac_rate_local(372) 
  reac_source_local(62,373) = + reac_rate_local(373) 
  reac_source_local(63,373) = - reac_rate_local(373) 
  reac_source_local(70,373) = - reac_rate_local(373) 
  reac_source_local(71,373) = + reac_rate_local(373) 
  reac_source_local(62,374) = - reac_rate_local(374) 
  reac_source_local(63,374) = + reac_rate_local(374) 
  reac_source_local(70,374) = + reac_rate_local(374) 
  reac_source_local(71,374) = - reac_rate_local(374) 
  reac_source_local(01,375) = + reac_rate_local(375) 
  reac_source_local(23,375) = - reac_rate_local(375) 
  reac_source_local(63,375) = + reac_rate_local(375) 
  reac_source_local(72,375) = - reac_rate_local(375) 
  reac_source_local(01,376) = + reac_rate_local(376) 
  reac_source_local(23,376) = - reac_rate_local(376) 
  reac_source_local(63,376) = + reac_rate_local(376) * 2.d0
  reac_source_local(71,376) = - reac_rate_local(376) 
  reac_source_local(01,377) = + reac_rate_local(377) 
  reac_source_local(23,377) = - reac_rate_local(377) 
  reac_source_local(62,377) = + reac_rate_local(377) 
  reac_source_local(71,377) = - reac_rate_local(377) 
  reac_source_local(23,378) = + reac_rate_local(378) 
  reac_source_local(71,378) = + reac_rate_local(378) 
  reac_source_local(72,378) = - reac_rate_local(378) * 2.d0
  reac_source_local(01,379) = + reac_rate_local(379) 
  reac_source_local(63,379) = + reac_rate_local(379) * 2.d0
  reac_source_local(72,379) = - reac_rate_local(379) * 2.d0
  reac_source_local(01,380) = + reac_rate_local(380) 
  reac_source_local(62,380) = + reac_rate_local(380) 
  reac_source_local(72,380) = - reac_rate_local(380) * 2.d0
  reac_source_local(23,381) = + reac_rate_local(381) 
  reac_source_local(70,381) = + reac_rate_local(381) 
  reac_source_local(71,381) = - reac_rate_local(381) 
  reac_source_local(72,381) = - reac_rate_local(381) 
  reac_source_local(63,382) = - reac_rate_local(382) 
  reac_source_local(70,382) = + reac_rate_local(382) 
  reac_source_local(71,382) = - reac_rate_local(382) 
  reac_source_local(70,383) = + reac_rate_local(383) 
  reac_source_local(71,383) = - reac_rate_local(383) * 2.d0
  reac_source_local(72,383) = + reac_rate_local(383) 
  reac_source_local(01,384) = - reac_rate_local(384) 
  reac_source_local(23,384) = + reac_rate_local(384) * 2.d0
  reac_source_local(01,385) = - reac_rate_local(385) 
  reac_source_local(23,385) = + reac_rate_local(385) * 2.d0
  reac_source_local(01,386) = - reac_rate_local(386) 
  reac_source_local(23,386) = + reac_rate_local(386) * 2.d0
  reac_source_local(01,387) = - reac_rate_local(387) 
  reac_source_local(23,387) = + reac_rate_local(387) * 2.d0
  reac_source_local(01,388) = - reac_rate_local(388) 
  reac_source_local(23,388) = + reac_rate_local(388) * 2.d0
  reac_source_local(30,389) = - reac_rate_local(389) 
  reac_source_local(38,389) = + reac_rate_local(389) * 2.d0
  reac_source_local(30,390) = - reac_rate_local(390) 
  reac_source_local(38,390) = + reac_rate_local(390) * 2.d0
  reac_source_local(30,391) = - reac_rate_local(391) 
  reac_source_local(38,391) = + reac_rate_local(391) * 2.d0
  reac_source_local(30,392) = - reac_rate_local(392) 
  reac_source_local(38,392) = + reac_rate_local(392) * 2.d0
  reac_source_local(30,393) = - reac_rate_local(393) 
  reac_source_local(38,393) = + reac_rate_local(393) * 2.d0
  reac_source_local(23,394) = + reac_rate_local(394) 
  reac_source_local(38,394) = + reac_rate_local(394) 
  reac_source_local(50,394) = - reac_rate_local(394) 
  reac_source_local(23,395) = + reac_rate_local(395) 
  reac_source_local(38,395) = + reac_rate_local(395) 
  reac_source_local(50,395) = - reac_rate_local(395) 
  reac_source_local(23,396) = + reac_rate_local(396) 
  reac_source_local(38,396) = + reac_rate_local(396) 
  reac_source_local(50,396) = - reac_rate_local(396) 
  reac_source_local(23,397) = + reac_rate_local(397) 
  reac_source_local(38,397) = + reac_rate_local(397) 
  reac_source_local(50,397) = - reac_rate_local(397) 
  reac_source_local(23,398) = + reac_rate_local(398) 
  reac_source_local(38,398) = + reac_rate_local(398) 
  reac_source_local(50,398) = - reac_rate_local(398) 
  reac_source_local(30,399) = + reac_rate_local(399) 
  reac_source_local(38,399) = + reac_rate_local(399) 
  reac_source_local(41,399) = - reac_rate_local(399) 
  reac_source_local(30,400) = + reac_rate_local(400) 
  reac_source_local(38,400) = + reac_rate_local(400) 
  reac_source_local(41,400) = - reac_rate_local(400) 
  reac_source_local(30,401) = + reac_rate_local(401) 
  reac_source_local(38,401) = + reac_rate_local(401) 
  reac_source_local(41,401) = - reac_rate_local(401) 
  reac_source_local(30,402) = + reac_rate_local(402) 
  reac_source_local(38,402) = + reac_rate_local(402) 
  reac_source_local(41,402) = - reac_rate_local(402) 
  reac_source_local(01,403) = + reac_rate_local(403) 
  reac_source_local(38,403) = + reac_rate_local(403) 
  reac_source_local(51,403) = - reac_rate_local(403) 
  reac_source_local(01,404) = + reac_rate_local(404) 
  reac_source_local(38,404) = + reac_rate_local(404) 
  reac_source_local(51,404) = - reac_rate_local(404) 
  reac_source_local(01,405) = + reac_rate_local(405) 
  reac_source_local(38,405) = + reac_rate_local(405) 
  reac_source_local(51,405) = - reac_rate_local(405) 
  reac_source_local(01,406) = + reac_rate_local(406) 
  reac_source_local(38,406) = + reac_rate_local(406) 
  reac_source_local(51,406) = - reac_rate_local(406) 
  reac_source_local(38,407) = + reac_rate_local(407) 
  reac_source_local(50,407) = + reac_rate_local(407) 
  reac_source_local(52,407) = - reac_rate_local(407) 
  reac_source_local(38,408) = + reac_rate_local(408) 
  reac_source_local(50,408) = + reac_rate_local(408) 
  reac_source_local(52,408) = - reac_rate_local(408) 
  reac_source_local(38,409) = + reac_rate_local(409) 
  reac_source_local(50,409) = + reac_rate_local(409) 
  reac_source_local(52,409) = - reac_rate_local(409) 
  reac_source_local(38,410) = + reac_rate_local(410) 
  reac_source_local(50,410) = + reac_rate_local(410) 
  reac_source_local(52,410) = - reac_rate_local(410) 
  reac_source_local(38,411) = + reac_rate_local(411) 
  reac_source_local(52,411) = + reac_rate_local(411) 
  reac_source_local(53,411) = - reac_rate_local(411) 
  reac_source_local(38,412) = + reac_rate_local(412) 
  reac_source_local(52,412) = + reac_rate_local(412) 
  reac_source_local(53,412) = - reac_rate_local(412) 
  reac_source_local(38,413) = + reac_rate_local(413) 
  reac_source_local(52,413) = + reac_rate_local(413) 
  reac_source_local(53,413) = - reac_rate_local(413) 
  reac_source_local(38,414) = + reac_rate_local(414) 
  reac_source_local(52,414) = + reac_rate_local(414) 
  reac_source_local(53,414) = - reac_rate_local(414) 
  reac_source_local(38,415) = + reac_rate_local(415) 
  reac_source_local(52,415) = + reac_rate_local(415) 
  reac_source_local(53,415) = - reac_rate_local(415) 
  reac_source_local(30,416) = + reac_rate_local(416) 
  reac_source_local(50,416) = + reac_rate_local(416) 
  reac_source_local(53,416) = - reac_rate_local(416) 
  reac_source_local(30,417) = + reac_rate_local(417) 
  reac_source_local(50,417) = + reac_rate_local(417) 
  reac_source_local(53,417) = - reac_rate_local(417) 
  reac_source_local(30,418) = + reac_rate_local(418) 
  reac_source_local(50,418) = + reac_rate_local(418) 
  reac_source_local(53,418) = - reac_rate_local(418) 
  reac_source_local(30,419) = + reac_rate_local(419) 
  reac_source_local(50,419) = + reac_rate_local(419) 
  reac_source_local(53,419) = - reac_rate_local(419) 
  reac_source_local(30,420) = + reac_rate_local(420) 
  reac_source_local(50,420) = + reac_rate_local(420) 
  reac_source_local(53,420) = - reac_rate_local(420) 
  reac_source_local(52,421) = + reac_rate_local(421) 
  reac_source_local(53,421) = + reac_rate_local(421) 
  reac_source_local(54,421) = - reac_rate_local(421) 
  reac_source_local(01,422) = + reac_rate_local(422) 
  reac_source_local(23,422) = - reac_rate_local(422) * 2.d0
  reac_source_local(01,423) = + reac_rate_local(423) 
  reac_source_local(23,423) = - reac_rate_local(423) * 2.d0
  reac_source_local(01,424) = + reac_rate_local(424) 
  reac_source_local(23,424) = - reac_rate_local(424) * 2.d0
  reac_source_local(01,425) = + reac_rate_local(425) 
  reac_source_local(23,425) = - reac_rate_local(425) * 2.d0
  reac_source_local(01,426) = + reac_rate_local(426) 
  reac_source_local(23,426) = - reac_rate_local(426) * 2.d0
  reac_source_local(30,427) = + reac_rate_local(427) 
  reac_source_local(38,427) = - reac_rate_local(427) * 2.d0
  reac_source_local(30,428) = + reac_rate_local(428) 
  reac_source_local(38,428) = - reac_rate_local(428) * 2.d0
  reac_source_local(30,429) = + reac_rate_local(429) 
  reac_source_local(38,429) = - reac_rate_local(429) * 2.d0
  reac_source_local(30,430) = + reac_rate_local(430) 
  reac_source_local(38,430) = - reac_rate_local(430) * 2.d0
  reac_source_local(30,431) = + reac_rate_local(431) 
  reac_source_local(38,431) = - reac_rate_local(431) * 2.d0
  reac_source_local(23,432) = - reac_rate_local(432) 
  reac_source_local(38,432) = - reac_rate_local(432) 
  reac_source_local(50,432) = + reac_rate_local(432) 
  reac_source_local(23,433) = - reac_rate_local(433) 
  reac_source_local(38,433) = - reac_rate_local(433) 
  reac_source_local(50,433) = + reac_rate_local(433) 
  reac_source_local(23,434) = - reac_rate_local(434) 
  reac_source_local(38,434) = - reac_rate_local(434) 
  reac_source_local(50,434) = + reac_rate_local(434) 
  reac_source_local(23,435) = - reac_rate_local(435) 
  reac_source_local(38,435) = - reac_rate_local(435) 
  reac_source_local(50,435) = + reac_rate_local(435) 
  reac_source_local(23,436) = - reac_rate_local(436) 
  reac_source_local(38,436) = - reac_rate_local(436) 
  reac_source_local(50,436) = + reac_rate_local(436) 
  reac_source_local(30,437) = - reac_rate_local(437) 
  reac_source_local(38,437) = - reac_rate_local(437) 
  reac_source_local(41,437) = + reac_rate_local(437) 
  reac_source_local(30,438) = - reac_rate_local(438) 
  reac_source_local(38,438) = - reac_rate_local(438) 
  reac_source_local(41,438) = + reac_rate_local(438) 
  reac_source_local(30,439) = - reac_rate_local(439) 
  reac_source_local(38,439) = - reac_rate_local(439) 
  reac_source_local(41,439) = + reac_rate_local(439) 
  reac_source_local(30,440) = - reac_rate_local(440) 
  reac_source_local(38,440) = - reac_rate_local(440) 
  reac_source_local(41,440) = + reac_rate_local(440) 
  reac_source_local(30,441) = - reac_rate_local(441) 
  reac_source_local(38,441) = - reac_rate_local(441) 
  reac_source_local(41,441) = + reac_rate_local(441) 
  reac_source_local(01,442) = - reac_rate_local(442) 
  reac_source_local(38,442) = - reac_rate_local(442) 
  reac_source_local(51,442) = + reac_rate_local(442) 
  reac_source_local(38,443) = - reac_rate_local(443) 
  reac_source_local(50,443) = - reac_rate_local(443) 
  reac_source_local(52,443) = + reac_rate_local(443) 
  reac_source_local(38,444) = - reac_rate_local(444) 
  reac_source_local(50,444) = - reac_rate_local(444) 
  reac_source_local(52,444) = + reac_rate_local(444) 
  reac_source_local(38,445) = - reac_rate_local(445) 
  reac_source_local(50,445) = - reac_rate_local(445) 
  reac_source_local(52,445) = + reac_rate_local(445) 
  reac_source_local(38,446) = - reac_rate_local(446) 
  reac_source_local(52,446) = - reac_rate_local(446) 
  reac_source_local(53,446) = + reac_rate_local(446) 
  reac_source_local(38,447) = - reac_rate_local(447) 
  reac_source_local(52,447) = - reac_rate_local(447) 
  reac_source_local(53,447) = + reac_rate_local(447) 
  reac_source_local(38,448) = - reac_rate_local(448) 
  reac_source_local(52,448) = - reac_rate_local(448) 
  reac_source_local(53,448) = + reac_rate_local(448) 
  reac_source_local(38,449) = - reac_rate_local(449) 
  reac_source_local(52,449) = - reac_rate_local(449) 
  reac_source_local(53,449) = + reac_rate_local(449) 
  reac_source_local(38,450) = - reac_rate_local(450) 
  reac_source_local(52,450) = - reac_rate_local(450) 
  reac_source_local(53,450) = + reac_rate_local(450) 
  reac_source_local(52,451) = - reac_rate_local(451) 
  reac_source_local(53,451) = - reac_rate_local(451) 
  reac_source_local(54,451) = + reac_rate_local(451) 
  reac_source_local(01,452) = + reac_rate_local(452) 
  reac_source_local(23,452) = - reac_rate_local(452) * 2.d0
  reac_source_local(01,453) = + reac_rate_local(453) 
  reac_source_local(23,453) = - reac_rate_local(453) * 2.d0
  reac_source_local(01,454) = + reac_rate_local(454) 
  reac_source_local(23,454) = - reac_rate_local(454) * 2.d0
  reac_source_local(62,455) = + reac_rate_local(455) 
  reac_source_local(63,455) = - reac_rate_local(455) * 2.d0
  reac_source_local(62,456) = + reac_rate_local(456) 
  reac_source_local(63,456) = - reac_rate_local(456) * 2.d0
  reac_source_local(23,457) = - reac_rate_local(457) 
  reac_source_local(63,457) = - reac_rate_local(457) 
  reac_source_local(72,457) = + reac_rate_local(457) 
  reac_source_local(23,458) = - reac_rate_local(458) 
  reac_source_local(63,458) = - reac_rate_local(458) 
  reac_source_local(72,458) = + reac_rate_local(458) 
  reac_source_local(23,459) = - reac_rate_local(459) 
  reac_source_local(63,459) = - reac_rate_local(459) 
  reac_source_local(72,459) = + reac_rate_local(459) 
  reac_source_local(23,460) = - reac_rate_local(460) 
  reac_source_local(62,460) = - reac_rate_local(460) 
  reac_source_local(71,460) = + reac_rate_local(460) 
  reac_source_local(23,461) = - reac_rate_local(461) 
  reac_source_local(62,461) = - reac_rate_local(461) 
  reac_source_local(71,461) = + reac_rate_local(461) 
  reac_source_local(23,462) = - reac_rate_local(462) 
  reac_source_local(62,462) = - reac_rate_local(462) 
  reac_source_local(71,462) = + reac_rate_local(462) 
  reac_source_local(63,463) = - reac_rate_local(463) 
  reac_source_local(71,463) = + reac_rate_local(463) 
  reac_source_local(72,463) = - reac_rate_local(463) 
  reac_source_local(63,464) = - reac_rate_local(464) 
  reac_source_local(71,464) = + reac_rate_local(464) 
  reac_source_local(72,464) = - reac_rate_local(464) 
  reac_source_local(63,465) = - reac_rate_local(465) 
  reac_source_local(71,465) = + reac_rate_local(465) 
  reac_source_local(72,465) = - reac_rate_local(465) 
  reac_source_local(63,466) = - reac_rate_local(466) 
  reac_source_local(70,466) = + reac_rate_local(466) 
  reac_source_local(71,466) = - reac_rate_local(466) 
  reac_source_local(63,467) = - reac_rate_local(467) 
  reac_source_local(70,467) = + reac_rate_local(467) 
  reac_source_local(71,467) = - reac_rate_local(467) 
  reac_source_local(63,468) = - reac_rate_local(468) 
  reac_source_local(70,468) = + reac_rate_local(468) 
  reac_source_local(71,468) = - reac_rate_local(468) 
  reac_source_local(62,469) = - reac_rate_local(469) 
  reac_source_local(70,469) = + reac_rate_local(469) 
  reac_source_local(72,469) = - reac_rate_local(469) 
  reac_source_local(62,470) = - reac_rate_local(470) 
  reac_source_local(70,470) = + reac_rate_local(470) 
  reac_source_local(72,470) = - reac_rate_local(470) 
  reac_source_local(62,471) = - reac_rate_local(471) 
  reac_source_local(70,471) = + reac_rate_local(471) 
  reac_source_local(72,471) = - reac_rate_local(471) 
  reac_source_local(23,472) = + reac_rate_local(472) 
  reac_source_local(26,472) = - reac_rate_local(472) 
  reac_source_local(38,472) = - reac_rate_local(472) 
  reac_source_local(42,472) = + reac_rate_local(472) 
  reac_source_local(23,473) = + reac_rate_local(473) 
  reac_source_local(26,473) = - reac_rate_local(473) 
  reac_source_local(30,473) = - reac_rate_local(473) 
  reac_source_local(43,473) = + reac_rate_local(473) 
  reac_source_local(26,474) = - reac_rate_local(474) 
  reac_source_local(30,474) = - reac_rate_local(474) 
  reac_source_local(38,474) = + reac_rate_local(474) 
  reac_source_local(55,474) = + reac_rate_local(474) 
  reac_source_local(26,475) = - reac_rate_local(475) 
  reac_source_local(30,475) = - reac_rate_local(475) 
  reac_source_local(42,475) = + reac_rate_local(475) 
  reac_source_local(50,475) = + reac_rate_local(475) 
  reac_source_local(26,476) = - reac_rate_local(476) 
  reac_source_local(30,476) = + reac_rate_local(476) 
  reac_source_local(41,476) = - reac_rate_local(476) 
  reac_source_local(55,476) = + reac_rate_local(476) 
  reac_source_local(23,477) = + reac_rate_local(477) 
  reac_source_local(26,477) = - reac_rate_local(477) 
  reac_source_local(50,477) = - reac_rate_local(477) 
  reac_source_local(55,477) = + reac_rate_local(477) 
  reac_source_local(26,478) = - reac_rate_local(478) 
  reac_source_local(27,478) = + reac_rate_local(478) 
  reac_source_local(38,478) = + reac_rate_local(478) 
  reac_source_local(50,478) = - reac_rate_local(478) 
  reac_source_local(01,479) = + reac_rate_local(479) 
  reac_source_local(26,479) = - reac_rate_local(479) 
  reac_source_local(42,479) = + reac_rate_local(479) 
  reac_source_local(50,479) = - reac_rate_local(479) 
  reac_source_local(01,480) = + reac_rate_local(480) 
  reac_source_local(26,480) = - reac_rate_local(480) 
  reac_source_local(51,480) = - reac_rate_local(480) 
  reac_source_local(55,480) = + reac_rate_local(480) 
  reac_source_local(01,481) = - reac_rate_local(481) 
  reac_source_local(23,481) = + reac_rate_local(481) 
  reac_source_local(42,481) = - reac_rate_local(481) 
  reac_source_local(55,481) = + reac_rate_local(481) 
  reac_source_local(30,482) = - reac_rate_local(482) 
  reac_source_local(38,482) = + reac_rate_local(482) 
  reac_source_local(42,482) = - reac_rate_local(482) 
  reac_source_local(43,482) = + reac_rate_local(482) 
  reac_source_local(30,483) = + reac_rate_local(483) 
  reac_source_local(41,483) = - reac_rate_local(483) 
  reac_source_local(42,483) = - reac_rate_local(483) 
  reac_source_local(43,483) = + reac_rate_local(483) 
  reac_source_local(38,484) = + reac_rate_local(484) 
  reac_source_local(42,484) = - reac_rate_local(484) 
  reac_source_local(50,484) = - reac_rate_local(484) 
  reac_source_local(55,484) = + reac_rate_local(484) 
  reac_source_local(23,485) = + reac_rate_local(485) 
  reac_source_local(42,485) = - reac_rate_local(485) 
  reac_source_local(43,485) = + reac_rate_local(485) 
  reac_source_local(50,485) = - reac_rate_local(485) 
  reac_source_local(24,486) = - reac_rate_local(486) 
  reac_source_local(26,486) = + reac_rate_local(486) 
  reac_source_local(38,486) = + reac_rate_local(486) 
  reac_source_local(42,486) = - reac_rate_local(486) 
  reac_source_local(42,487) = - reac_rate_local(487) 
  reac_source_local(50,487) = + reac_rate_local(487) 
  reac_source_local(51,487) = - reac_rate_local(487) 
  reac_source_local(55,487) = + reac_rate_local(487) 
  reac_source_local(38,488) = + reac_rate_local(488) 
  reac_source_local(42,488) = - reac_rate_local(488) 
  reac_source_local(51,488) = - reac_rate_local(488) 
  reac_source_local(56,488) = + reac_rate_local(488) 
  reac_source_local(01,489) = + reac_rate_local(489) 
  reac_source_local(42,489) = - reac_rate_local(489) 
  reac_source_local(43,489) = + reac_rate_local(489) 
  reac_source_local(51,489) = - reac_rate_local(489) 
  reac_source_local(38,490) = + reac_rate_local(490) 
  reac_source_local(42,490) = - reac_rate_local(490) 
  reac_source_local(52,490) = - reac_rate_local(490) 
  reac_source_local(57,490) = + reac_rate_local(490) 
  reac_source_local(01,491) = + reac_rate_local(491) 
  reac_source_local(27,491) = - reac_rate_local(491) 
  reac_source_local(30,491) = - reac_rate_local(491) 
  reac_source_local(43,491) = + reac_rate_local(491) 
  reac_source_local(23,492) = + reac_rate_local(492) 
  reac_source_local(27,492) = - reac_rate_local(492) 
  reac_source_local(38,492) = - reac_rate_local(492) 
  reac_source_local(55,492) = + reac_rate_local(492) 
  reac_source_local(01,493) = + reac_rate_local(493) 
  reac_source_local(27,493) = - reac_rate_local(493) 
  reac_source_local(38,493) = + reac_rate_local(493) 
  reac_source_local(41,493) = - reac_rate_local(493) 
  reac_source_local(43,493) = + reac_rate_local(493) 
  reac_source_local(01,494) = + reac_rate_local(494) 
  reac_source_local(23,494) = - reac_rate_local(494) 
  reac_source_local(26,494) = + reac_rate_local(494) 
  reac_source_local(27,494) = - reac_rate_local(494) 
  reac_source_local(01,495) = + reac_rate_local(495) 
  reac_source_local(27,495) = - reac_rate_local(495) 
  reac_source_local(50,495) = - reac_rate_local(495) 
  reac_source_local(55,495) = + reac_rate_local(495) 
  reac_source_local(01,496) = + reac_rate_local(496) 
  reac_source_local(27,496) = - reac_rate_local(496) 
  reac_source_local(51,496) = - reac_rate_local(496) 
  reac_source_local(56,496) = + reac_rate_local(496) 
  reac_source_local(01,497) = + reac_rate_local(497) 
  reac_source_local(23,497) = + reac_rate_local(497) 
  reac_source_local(27,497) = - reac_rate_local(497) 
  reac_source_local(51,497) = - reac_rate_local(497) 
  reac_source_local(55,497) = + reac_rate_local(497) 
  reac_source_local(01,498) = - reac_rate_local(498) 
  reac_source_local(43,498) = - reac_rate_local(498) 
  reac_source_local(50,498) = + reac_rate_local(498) 
  reac_source_local(55,498) = + reac_rate_local(498) 
  reac_source_local(23,499) = - reac_rate_local(499) 
  reac_source_local(38,499) = + reac_rate_local(499) 
  reac_source_local(43,499) = - reac_rate_local(499) 
  reac_source_local(55,499) = + reac_rate_local(499) 
  reac_source_local(30,500) = + reac_rate_local(500) 
  reac_source_local(43,500) = - reac_rate_local(500) 
  reac_source_local(50,500) = - reac_rate_local(500) 
  reac_source_local(55,500) = + reac_rate_local(500) 
  reac_source_local(41,501) = + reac_rate_local(501) 
  reac_source_local(43,501) = - reac_rate_local(501) 
  reac_source_local(52,501) = - reac_rate_local(501) 
  reac_source_local(55,501) = + reac_rate_local(501) 
  reac_source_local(30,502) = + reac_rate_local(502) 
  reac_source_local(43,502) = - reac_rate_local(502) 
  reac_source_local(52,502) = - reac_rate_local(502) 
  reac_source_local(57,502) = + reac_rate_local(502) 
  reac_source_local(01,503) = + reac_rate_local(503) 
  reac_source_local(23,503) = + reac_rate_local(503) 
  reac_source_local(28,503) = - reac_rate_local(503) 
  reac_source_local(30,503) = - reac_rate_local(503) 
  reac_source_local(43,503) = + reac_rate_local(503) 
  reac_source_local(01,504) = + reac_rate_local(504) 
  reac_source_local(28,504) = - reac_rate_local(504) 
  reac_source_local(30,504) = - reac_rate_local(504) 
  reac_source_local(57,504) = + reac_rate_local(504) 
  reac_source_local(01,505) = + reac_rate_local(505) 
  reac_source_local(23,505) = - reac_rate_local(505) 
  reac_source_local(27,505) = + reac_rate_local(505) 
  reac_source_local(28,505) = - reac_rate_local(505) 
  reac_source_local(01,506) = + reac_rate_local(506) 
  reac_source_local(23,506) = + reac_rate_local(506) 
  reac_source_local(28,506) = - reac_rate_local(506) 
  reac_source_local(50,506) = - reac_rate_local(506) 
  reac_source_local(55,506) = + reac_rate_local(506) 
  reac_source_local(01,507) = + reac_rate_local(507) 
  reac_source_local(28,507) = - reac_rate_local(507) 
  reac_source_local(50,507) = - reac_rate_local(507) 
  reac_source_local(56,507) = + reac_rate_local(507) 
  reac_source_local(50,508) = - reac_rate_local(508) 
  reac_source_local(52,508) = + reac_rate_local(508) 
  reac_source_local(55,508) = + reac_rate_local(508) 
  reac_source_local(57,508) = - reac_rate_local(508) 
  reac_source_local(50,509) = - reac_rate_local(509) 
  reac_source_local(51,509) = + reac_rate_local(509) 
  reac_source_local(55,509) = + reac_rate_local(509) 
  reac_source_local(56,509) = - reac_rate_local(509) 
  reac_source_local(01,510) = + reac_rate_local(510) 
  reac_source_local(27,510) = + reac_rate_local(510) 
  reac_source_local(29,510) = - reac_rate_local(510) 
  reac_source_local(01,511) = + reac_rate_local(511) * 2.d0
  reac_source_local(29,511) = - reac_rate_local(511) 
  reac_source_local(30,511) = - reac_rate_local(511) 
  reac_source_local(43,511) = + reac_rate_local(511) 
  reac_source_local(01,512) = + reac_rate_local(512) * 2.d0
  reac_source_local(29,512) = - reac_rate_local(512) 
  reac_source_local(38,512) = - reac_rate_local(512) 
  reac_source_local(42,512) = + reac_rate_local(512) 
  reac_source_local(01,513) = + reac_rate_local(513) * 2.d0
  reac_source_local(23,513) = - reac_rate_local(513) 
  reac_source_local(26,513) = + reac_rate_local(513) 
  reac_source_local(29,513) = - reac_rate_local(513) 
  reac_source_local(01,514) = + reac_rate_local(514) * 2.d0
  reac_source_local(29,514) = - reac_rate_local(514) 
  reac_source_local(50,514) = - reac_rate_local(514) 
  reac_source_local(55,514) = + reac_rate_local(514) 
  reac_source_local(01,515) = - reac_rate_local(515) 
  reac_source_local(30,515) = + reac_rate_local(515) 
  reac_source_local(45,515) = - reac_rate_local(515) 
  reac_source_local(78,515) = + reac_rate_local(515) 
  reac_source_local(30,516) = + reac_rate_local(516) 
  reac_source_local(43,516) = + reac_rate_local(516) 
  reac_source_local(45,516) = - reac_rate_local(516) 
  reac_source_local(30,517) = + reac_rate_local(517) * 2.d0
  reac_source_local(35,517) = - reac_rate_local(517) 
  reac_source_local(43,517) = + reac_rate_local(517) 
  reac_source_local(45,517) = - reac_rate_local(517) 
  reac_source_local(30,518) = + reac_rate_local(518) * 2.d0
  reac_source_local(36,518) = - reac_rate_local(518) 
  reac_source_local(43,518) = + reac_rate_local(518) 
  reac_source_local(45,518) = - reac_rate_local(518) 
  reac_source_local(38,519) = - reac_rate_local(519) 
  reac_source_local(41,519) = + reac_rate_local(519) 
  reac_source_local(43,519) = + reac_rate_local(519) 
  reac_source_local(45,519) = - reac_rate_local(519) 
  reac_source_local(30,520) = + reac_rate_local(520) * 2.d0
  reac_source_local(45,520) = - reac_rate_local(520) 
  reac_source_local(50,520) = - reac_rate_local(520) 
  reac_source_local(55,520) = + reac_rate_local(520) 
  reac_source_local(01,521) = + reac_rate_local(521) 
  reac_source_local(43,521) = + reac_rate_local(521) 
  reac_source_local(78,521) = - reac_rate_local(521) 
  reac_source_local(01,522) = + reac_rate_local(522) 
  reac_source_local(30,522) = - reac_rate_local(522) 
  reac_source_local(45,522) = + reac_rate_local(522) 
  reac_source_local(78,522) = - reac_rate_local(522) 
  reac_source_local(01,523) = - reac_rate_local(523) 
  reac_source_local(26,523) = - reac_rate_local(523) 
  reac_source_local(28,523) = + reac_rate_local(523) 
  reac_source_local(26,524) = - reac_rate_local(524) 
  reac_source_local(38,524) = - reac_rate_local(524) 
  reac_source_local(55,524) = + reac_rate_local(524) 
  reac_source_local(23,525) = - reac_rate_local(525) 
  reac_source_local(26,525) = - reac_rate_local(525) 
  reac_source_local(27,525) = + reac_rate_local(525) 
  reac_source_local(01,526) = - reac_rate_local(526) 
  reac_source_local(23,526) = + reac_rate_local(526) 
  reac_source_local(42,526) = - reac_rate_local(526) 
  reac_source_local(55,526) = + reac_rate_local(526) 
  reac_source_local(38,527) = - reac_rate_local(527) 
  reac_source_local(42,527) = - reac_rate_local(527) 
  reac_source_local(43,527) = + reac_rate_local(527) 
  reac_source_local(23,528) = - reac_rate_local(528) 
  reac_source_local(42,528) = - reac_rate_local(528) 
  reac_source_local(55,528) = + reac_rate_local(528) 
  reac_source_local(01,529) = - reac_rate_local(529) 
  reac_source_local(27,529) = - reac_rate_local(529) 
  reac_source_local(29,529) = + reac_rate_local(529) 
  reac_source_local(23,530) = - reac_rate_local(530) 
  reac_source_local(27,530) = - reac_rate_local(530) 
  reac_source_local(28,530) = + reac_rate_local(530) 
  reac_source_local(30,531) = - reac_rate_local(531) 
  reac_source_local(43,531) = - reac_rate_local(531) 
  reac_source_local(45,531) = + reac_rate_local(531) 
  reac_source_local(01,532) = - reac_rate_local(532) 
  reac_source_local(43,532) = - reac_rate_local(532) 
  reac_source_local(78,532) = + reac_rate_local(532) 
  reac_source_local(26,533) = - reac_rate_local(533) 
  reac_source_local(62,533) = - reac_rate_local(533) 
  reac_source_local(63,533) = + reac_rate_local(533) 
  reac_source_local(75,533) = + reac_rate_local(533) 
  reac_source_local(26,534) = - reac_rate_local(534) 
  reac_source_local(70,534) = - reac_rate_local(534) 
  reac_source_local(72,534) = + reac_rate_local(534) 
  reac_source_local(74,534) = + reac_rate_local(534) 
  reac_source_local(26,535) = - reac_rate_local(535) 
  reac_source_local(62,535) = + reac_rate_local(535) 
  reac_source_local(70,535) = - reac_rate_local(535) 
  reac_source_local(76,535) = + reac_rate_local(535) 
  reac_source_local(63,536) = + reac_rate_local(536) 
  reac_source_local(64,536) = - reac_rate_local(536) 
  reac_source_local(70,536) = - reac_rate_local(536) 
  reac_source_local(73,536) = + reac_rate_local(536) 
  reac_source_local(62,537) = + reac_rate_local(537) 
  reac_source_local(63,537) = - reac_rate_local(537) 
  reac_source_local(64,537) = + reac_rate_local(537) 
  reac_source_local(65,537) = - reac_rate_local(537) 
  reac_source_local(62,538) = - reac_rate_local(538) 
  reac_source_local(63,538) = + reac_rate_local(538) 
  reac_source_local(65,538) = - reac_rate_local(538) 
  reac_source_local(66,538) = + reac_rate_local(538) 
  reac_source_local(62,539) = + reac_rate_local(539) 
  reac_source_local(65,539) = - reac_rate_local(539) 
  reac_source_local(70,539) = - reac_rate_local(539) 
  reac_source_local(73,539) = + reac_rate_local(539) 
  reac_source_local(01,540) = - reac_rate_local(540) 
  reac_source_local(63,540) = + reac_rate_local(540) 
  reac_source_local(65,540) = - reac_rate_local(540) 
  reac_source_local(76,540) = + reac_rate_local(540) 
  reac_source_local(23,541) = + reac_rate_local(541) 
  reac_source_local(62,541) = - reac_rate_local(541) 
  reac_source_local(66,541) = + reac_rate_local(541) 
  reac_source_local(75,541) = - reac_rate_local(541) 
  reac_source_local(62,542) = - reac_rate_local(542) 
  reac_source_local(63,542) = + reac_rate_local(542) 
  reac_source_local(74,542) = + reac_rate_local(542) 
  reac_source_local(75,542) = - reac_rate_local(542) 
  reac_source_local(23,544) = + reac_rate_local(544) 
  reac_source_local(70,544) = - reac_rate_local(544) 
  reac_source_local(75,544) = - reac_rate_local(544) 
  reac_source_local(77,544) = + reac_rate_local(544) 
  reac_source_local(01,545) = - reac_rate_local(545) 
  reac_source_local(23,545) = + reac_rate_local(545) 
  reac_source_local(75,545) = - reac_rate_local(545) 
  reac_source_local(76,545) = + reac_rate_local(545) 
  reac_source_local(62,546) = - reac_rate_local(546) 
  reac_source_local(63,546) = + reac_rate_local(546) 
  reac_source_local(73,546) = + reac_rate_local(546) 
  reac_source_local(74,546) = - reac_rate_local(546) 
  reac_source_local(70,547) = - reac_rate_local(547) 
  reac_source_local(71,547) = + reac_rate_local(547) 
  reac_source_local(73,547) = + reac_rate_local(547) 
  reac_source_local(74,547) = - reac_rate_local(547) 
  reac_source_local(70,548) = - reac_rate_local(548) 
  reac_source_local(72,548) = + reac_rate_local(548) 
  reac_source_local(74,548) = - reac_rate_local(548) 
  reac_source_local(77,548) = + reac_rate_local(548) 
  reac_source_local(01,549) = + reac_rate_local(549) 
  reac_source_local(70,549) = - reac_rate_local(549) 
  reac_source_local(76,549) = - reac_rate_local(549) 
  reac_source_local(77,549) = + reac_rate_local(549) 
  reac_source_local(01,550) = + reac_rate_local(550) 
  reac_source_local(27,550) = - reac_rate_local(550) 
  reac_source_local(70,550) = - reac_rate_local(550) 
  reac_source_local(73,550) = + reac_rate_local(550) 
  reac_source_local(23,551) = + reac_rate_local(551) 
  reac_source_local(26,551) = - reac_rate_local(551) 
  reac_source_local(70,551) = - reac_rate_local(551) 
  reac_source_local(73,551) = + reac_rate_local(551) 
  reac_source_local(30,552) = + reac_rate_local(552) 
  reac_source_local(43,552) = - reac_rate_local(552) 
  reac_source_local(70,552) = - reac_rate_local(552) 
  reac_source_local(73,552) = + reac_rate_local(552) 
  reac_source_local(68,553) = + reac_rate_local(553) 
  reac_source_local(69,553) = - reac_rate_local(553) 
  reac_source_local(70,553) = - reac_rate_local(553) 
  reac_source_local(77,553) = + reac_rate_local(553) 
  reac_source_local(70,554) = - reac_rate_local(554) 
  reac_source_local(71,554) = + reac_rate_local(554) 
  reac_source_local(73,554) = - reac_rate_local(554) 
  reac_source_local(77,554) = + reac_rate_local(554) 
  reac_source_local(70,555) = + reac_rate_local(555) 
  reac_source_local(71,555) = - reac_rate_local(555) 
  reac_source_local(73,555) = - reac_rate_local(555) 
  reac_source_local(74,555) = + reac_rate_local(555) 
  reac_source_local(35,556) = - reac_rate_local(556) 
  reac_source_local(38,556) = + reac_rate_local(556) 
  reac_source_local(46,556) = - reac_rate_local(556) 
  reac_source_local(47,556) = + reac_rate_local(556) 
  reac_source_local(38,557) = + reac_rate_local(557) 
  reac_source_local(41,557) = - reac_rate_local(557) 
  reac_source_local(46,557) = - reac_rate_local(557) 
  reac_source_local(48,557) = + reac_rate_local(557) 
  reac_source_local(38,558) = + reac_rate_local(558) 
  reac_source_local(46,558) = - reac_rate_local(558) 
  reac_source_local(52,558) = - reac_rate_local(558) 
  reac_source_local(60,558) = + reac_rate_local(558) 
  reac_source_local(46,559) = - reac_rate_local(559) 
  reac_source_local(50,559) = + reac_rate_local(559) 
  reac_source_local(51,559) = - reac_rate_local(559) 
  reac_source_local(58,559) = + reac_rate_local(559) 
  reac_source_local(38,560) = + reac_rate_local(560) 
  reac_source_local(46,560) = - reac_rate_local(560) 
  reac_source_local(51,560) = - reac_rate_local(560) 
  reac_source_local(59,560) = + reac_rate_local(560) 
  reac_source_local(30,561) = + reac_rate_local(561) 
  reac_source_local(38,561) = - reac_rate_local(561) 
  reac_source_local(46,561) = + reac_rate_local(561) 
  reac_source_local(47,561) = - reac_rate_local(561) 
  reac_source_local(30,562) = + reac_rate_local(562) 
  reac_source_local(41,562) = - reac_rate_local(562) 
  reac_source_local(47,562) = - reac_rate_local(562) 
  reac_source_local(48,562) = + reac_rate_local(562) 
  reac_source_local(30,563) = + reac_rate_local(563) 
  reac_source_local(47,563) = - reac_rate_local(563) 
  reac_source_local(52,563) = - reac_rate_local(563) 
  reac_source_local(60,563) = + reac_rate_local(563) 
  reac_source_local(30,564) = + reac_rate_local(564) 
  reac_source_local(47,564) = - reac_rate_local(564) 
  reac_source_local(53,564) = - reac_rate_local(564) 
  reac_source_local(61,564) = + reac_rate_local(564) 
  reac_source_local(30,565) = + reac_rate_local(565) 
  reac_source_local(38,565) = - reac_rate_local(565) 
  reac_source_local(47,565) = + reac_rate_local(565) 
  reac_source_local(48,565) = - reac_rate_local(565) 
  reac_source_local(38,566) = + reac_rate_local(566) 
  reac_source_local(48,566) = - reac_rate_local(566) 
  reac_source_local(50,566) = - reac_rate_local(566) 
  reac_source_local(61,566) = + reac_rate_local(566) 
  reac_source_local(30,567) = + reac_rate_local(567) 
  reac_source_local(48,567) = - reac_rate_local(567) 
  reac_source_local(50,567) = - reac_rate_local(567) 
  reac_source_local(60,567) = + reac_rate_local(567) 
  reac_source_local(41,568) = + reac_rate_local(568) 
  reac_source_local(48,568) = - reac_rate_local(568) 
  reac_source_local(52,568) = - reac_rate_local(568) 
  reac_source_local(60,568) = + reac_rate_local(568) 
  reac_source_local(30,569) = + reac_rate_local(569) 
  reac_source_local(48,569) = - reac_rate_local(569) 
  reac_source_local(52,569) = - reac_rate_local(569) 
  reac_source_local(61,569) = + reac_rate_local(569) 
  reac_source_local(41,570) = + reac_rate_local(570) 
  reac_source_local(48,570) = - reac_rate_local(570) 
  reac_source_local(53,570) = - reac_rate_local(570) 
  reac_source_local(61,570) = + reac_rate_local(570) 
  reac_source_local(30,571) = - reac_rate_local(571) 
  reac_source_local(47,571) = + reac_rate_local(571) 
  reac_source_local(50,571) = + reac_rate_local(571) 
  reac_source_local(58,571) = - reac_rate_local(571) 
  reac_source_local(50,572) = + reac_rate_local(572) 
  reac_source_local(52,572) = - reac_rate_local(572) 
  reac_source_local(58,572) = - reac_rate_local(572) 
  reac_source_local(60,572) = + reac_rate_local(572) 
  reac_source_local(01,573) = + reac_rate_local(573) 
  reac_source_local(51,573) = - reac_rate_local(573) 
  reac_source_local(58,573) = - reac_rate_local(573) 
  reac_source_local(60,573) = + reac_rate_local(573) 
  reac_source_local(30,574) = + reac_rate_local(574) 
  reac_source_local(41,574) = - reac_rate_local(574) 
  reac_source_local(60,574) = - reac_rate_local(574) 
  reac_source_local(61,574) = + reac_rate_local(574) 
  reac_source_local(50,575) = + reac_rate_local(575) 
  reac_source_local(52,575) = - reac_rate_local(575) 
  reac_source_local(60,575) = - reac_rate_local(575) 
  reac_source_local(61,575) = + reac_rate_local(575) 
  reac_source_local(52,576) = + reac_rate_local(576) 
  reac_source_local(53,576) = - reac_rate_local(576) 
  reac_source_local(60,576) = - reac_rate_local(576) 
  reac_source_local(61,576) = + reac_rate_local(576) 
  reac_source_local(52,577) = + reac_rate_local(577) * 2.d0
  reac_source_local(54,577) = - reac_rate_local(577) 
  reac_source_local(60,577) = - reac_rate_local(577) 
  reac_source_local(61,577) = + reac_rate_local(577) 
  reac_source_local(50,578) = - reac_rate_local(578) 
  reac_source_local(52,578) = + reac_rate_local(578) 
  reac_source_local(60,578) = + reac_rate_local(578) 
  reac_source_local(61,578) = - reac_rate_local(578) 
  reac_source_local(30,579) = + reac_rate_local(579) 
  reac_source_local(47,579) = + reac_rate_local(579) 
  reac_source_local(49,579) = - reac_rate_local(579) 
  reac_source_local(30,580) = + reac_rate_local(580) 
  reac_source_local(47,580) = + reac_rate_local(580) 
  reac_source_local(49,580) = - reac_rate_local(580) 
  reac_source_local(30,581) = + reac_rate_local(581) 
  reac_source_local(38,581) = - reac_rate_local(581) 
  reac_source_local(48,581) = + reac_rate_local(581) 
  reac_source_local(49,581) = - reac_rate_local(581) 
  reac_source_local(30,582) = + reac_rate_local(582) * 2.d0
  reac_source_local(38,582) = - reac_rate_local(582) 
  reac_source_local(46,582) = + reac_rate_local(582) 
  reac_source_local(49,582) = - reac_rate_local(582) 
  reac_source_local(30,583) = + reac_rate_local(583) * 2.d0
  reac_source_local(35,583) = - reac_rate_local(583) 
  reac_source_local(47,583) = + reac_rate_local(583) 
  reac_source_local(49,583) = - reac_rate_local(583) 
  reac_source_local(30,584) = + reac_rate_local(584) * 2.d0
  reac_source_local(36,584) = - reac_rate_local(584) 
  reac_source_local(47,584) = + reac_rate_local(584) 
  reac_source_local(49,584) = - reac_rate_local(584) 
  reac_source_local(30,585) = + reac_rate_local(585) 
  reac_source_local(49,585) = - reac_rate_local(585) 
  reac_source_local(50,585) = - reac_rate_local(585) 
  reac_source_local(61,585) = + reac_rate_local(585) 
  reac_source_local(30,586) = - reac_rate_local(586) 
  reac_source_local(46,586) = - reac_rate_local(586) 
  reac_source_local(48,586) = + reac_rate_local(586) 
  reac_source_local(46,587) = - reac_rate_local(587) 
  reac_source_local(50,587) = - reac_rate_local(587) 
  reac_source_local(60,587) = + reac_rate_local(587) 
  reac_source_local(30,588) = - reac_rate_local(588) 
  reac_source_local(47,588) = - reac_rate_local(588) 
  reac_source_local(49,588) = + reac_rate_local(588) 
  reac_source_local(23,589) = + reac_rate_local(589) 
  reac_source_local(26,589) = - reac_rate_local(589) 
  reac_source_local(38,589) = + reac_rate_local(589) 
  reac_source_local(46,589) = - reac_rate_local(589) 
  reac_source_local(01,590) = + reac_rate_local(590) 
  reac_source_local(27,590) = - reac_rate_local(590) 
  reac_source_local(38,590) = + reac_rate_local(590) 
  reac_source_local(46,590) = - reac_rate_local(590) 
  reac_source_local(38,591) = + reac_rate_local(591) * 2.d0
  reac_source_local(42,591) = - reac_rate_local(591) 
  reac_source_local(46,591) = - reac_rate_local(591) 
  reac_source_local(30,592) = + reac_rate_local(592) 
  reac_source_local(38,592) = + reac_rate_local(592) 
  reac_source_local(43,592) = - reac_rate_local(592) 
  reac_source_local(46,592) = - reac_rate_local(592) 
  reac_source_local(38,593) = + reac_rate_local(593) 
  reac_source_local(46,593) = - reac_rate_local(593) 
  reac_source_local(50,593) = + reac_rate_local(593) 
  reac_source_local(55,593) = - reac_rate_local(593) 
  reac_source_local(38,594) = + reac_rate_local(594) 
  reac_source_local(46,594) = - reac_rate_local(594) 
  reac_source_local(51,594) = + reac_rate_local(594) 
  reac_source_local(56,594) = - reac_rate_local(594) 
  reac_source_local(38,595) = + reac_rate_local(595) 
  reac_source_local(46,595) = - reac_rate_local(595) 
  reac_source_local(52,595) = + reac_rate_local(595) 
  reac_source_local(57,595) = - reac_rate_local(595) 
  reac_source_local(23,596) = + reac_rate_local(596) 
  reac_source_local(26,596) = - reac_rate_local(596) 
  reac_source_local(30,596) = + reac_rate_local(596) 
  reac_source_local(47,596) = - reac_rate_local(596) 
  reac_source_local(01,597) = + reac_rate_local(597) 
  reac_source_local(27,597) = - reac_rate_local(597) 
  reac_source_local(30,597) = + reac_rate_local(597) 
  reac_source_local(47,597) = - reac_rate_local(597) 
  reac_source_local(30,598) = + reac_rate_local(598) 
  reac_source_local(38,598) = + reac_rate_local(598) 
  reac_source_local(42,598) = - reac_rate_local(598) 
  reac_source_local(47,598) = - reac_rate_local(598) 
  reac_source_local(30,599) = + reac_rate_local(599) * 2.d0
  reac_source_local(43,599) = - reac_rate_local(599) 
  reac_source_local(47,599) = - reac_rate_local(599) 
  reac_source_local(30,600) = + reac_rate_local(600) 
  reac_source_local(47,600) = - reac_rate_local(600) 
  reac_source_local(50,600) = + reac_rate_local(600) 
  reac_source_local(55,600) = - reac_rate_local(600) 
  reac_source_local(30,601) = + reac_rate_local(601) 
  reac_source_local(47,601) = - reac_rate_local(601) 
  reac_source_local(51,601) = + reac_rate_local(601) 
  reac_source_local(56,601) = - reac_rate_local(601) 
  reac_source_local(30,602) = + reac_rate_local(602) 
  reac_source_local(47,602) = - reac_rate_local(602) 
  reac_source_local(52,602) = + reac_rate_local(602) 
  reac_source_local(57,602) = - reac_rate_local(602) 
  reac_source_local(23,603) = + reac_rate_local(603) 
  reac_source_local(26,603) = - reac_rate_local(603) 
  reac_source_local(41,603) = + reac_rate_local(603) 
  reac_source_local(48,603) = - reac_rate_local(603) 
  reac_source_local(01,604) = + reac_rate_local(604) 
  reac_source_local(27,604) = - reac_rate_local(604) 
  reac_source_local(41,604) = + reac_rate_local(604) 
  reac_source_local(48,604) = - reac_rate_local(604) 
  reac_source_local(38,605) = + reac_rate_local(605) 
  reac_source_local(41,605) = + reac_rate_local(605) 
  reac_source_local(42,605) = - reac_rate_local(605) 
  reac_source_local(48,605) = - reac_rate_local(605) 
  reac_source_local(30,606) = + reac_rate_local(606) 
  reac_source_local(41,606) = + reac_rate_local(606) 
  reac_source_local(43,606) = - reac_rate_local(606) 
  reac_source_local(48,606) = - reac_rate_local(606) 
  reac_source_local(41,607) = + reac_rate_local(607) 
  reac_source_local(48,607) = - reac_rate_local(607) 
  reac_source_local(50,607) = + reac_rate_local(607) 
  reac_source_local(55,607) = - reac_rate_local(607) 
  reac_source_local(41,608) = + reac_rate_local(608) 
  reac_source_local(48,608) = - reac_rate_local(608) 
  reac_source_local(51,608) = + reac_rate_local(608) 
  reac_source_local(56,608) = - reac_rate_local(608) 
  reac_source_local(41,609) = + reac_rate_local(609) 
  reac_source_local(48,609) = - reac_rate_local(609) 
  reac_source_local(52,609) = + reac_rate_local(609) 
  reac_source_local(57,609) = - reac_rate_local(609) 
  reac_source_local(23,610) = + reac_rate_local(610) 
  reac_source_local(26,610) = - reac_rate_local(610) 
  reac_source_local(50,610) = + reac_rate_local(610) 
  reac_source_local(58,610) = - reac_rate_local(610) 
  reac_source_local(01,611) = + reac_rate_local(611) 
  reac_source_local(27,611) = - reac_rate_local(611) 
  reac_source_local(50,611) = + reac_rate_local(611) 
  reac_source_local(58,611) = - reac_rate_local(611) 
  reac_source_local(38,612) = + reac_rate_local(612) 
  reac_source_local(42,612) = - reac_rate_local(612) 
  reac_source_local(50,612) = + reac_rate_local(612) 
  reac_source_local(58,612) = - reac_rate_local(612) 
  reac_source_local(30,613) = + reac_rate_local(613) 
  reac_source_local(43,613) = - reac_rate_local(613) 
  reac_source_local(50,613) = + reac_rate_local(613) 
  reac_source_local(58,613) = - reac_rate_local(613) 
  reac_source_local(50,614) = + reac_rate_local(614) * 2.d0
  reac_source_local(55,614) = - reac_rate_local(614) 
  reac_source_local(58,614) = - reac_rate_local(614) 
  reac_source_local(50,615) = + reac_rate_local(615) 
  reac_source_local(51,615) = + reac_rate_local(615) 
  reac_source_local(56,615) = - reac_rate_local(615) 
  reac_source_local(58,615) = - reac_rate_local(615) 
  reac_source_local(50,616) = + reac_rate_local(616) 
  reac_source_local(52,616) = + reac_rate_local(616) 
  reac_source_local(57,616) = - reac_rate_local(616) 
  reac_source_local(58,616) = - reac_rate_local(616) 
  reac_source_local(23,617) = + reac_rate_local(617) 
  reac_source_local(26,617) = - reac_rate_local(617) 
  reac_source_local(51,617) = + reac_rate_local(617) 
  reac_source_local(59,617) = - reac_rate_local(617) 
  reac_source_local(01,618) = + reac_rate_local(618) 
  reac_source_local(27,618) = - reac_rate_local(618) 
  reac_source_local(51,618) = + reac_rate_local(618) 
  reac_source_local(59,618) = - reac_rate_local(618) 
  reac_source_local(38,619) = + reac_rate_local(619) 
  reac_source_local(42,619) = - reac_rate_local(619) 
  reac_source_local(51,619) = + reac_rate_local(619) 
  reac_source_local(59,619) = - reac_rate_local(619) 
  reac_source_local(30,620) = + reac_rate_local(620) 
  reac_source_local(43,620) = - reac_rate_local(620) 
  reac_source_local(51,620) = + reac_rate_local(620) 
  reac_source_local(59,620) = - reac_rate_local(620) 
  reac_source_local(50,621) = + reac_rate_local(621) 
  reac_source_local(51,621) = + reac_rate_local(621) 
  reac_source_local(55,621) = - reac_rate_local(621) 
  reac_source_local(59,621) = - reac_rate_local(621) 
  reac_source_local(51,622) = + reac_rate_local(622) * 2.d0
  reac_source_local(56,622) = - reac_rate_local(622) 
  reac_source_local(59,622) = - reac_rate_local(622) 
  reac_source_local(51,623) = + reac_rate_local(623) 
  reac_source_local(52,623) = + reac_rate_local(623) 
  reac_source_local(57,623) = - reac_rate_local(623) 
  reac_source_local(59,623) = - reac_rate_local(623) 
  reac_source_local(23,624) = + reac_rate_local(624) 
  reac_source_local(26,624) = - reac_rate_local(624) 
  reac_source_local(52,624) = + reac_rate_local(624) 
  reac_source_local(60,624) = - reac_rate_local(624) 
  reac_source_local(01,625) = + reac_rate_local(625) 
  reac_source_local(27,625) = - reac_rate_local(625) 
  reac_source_local(52,625) = + reac_rate_local(625) 
  reac_source_local(60,625) = - reac_rate_local(625) 
  reac_source_local(38,626) = + reac_rate_local(626) 
  reac_source_local(42,626) = - reac_rate_local(626) 
  reac_source_local(52,626) = + reac_rate_local(626) 
  reac_source_local(60,626) = - reac_rate_local(626) 
  reac_source_local(30,627) = + reac_rate_local(627) 
  reac_source_local(43,627) = - reac_rate_local(627) 
  reac_source_local(52,627) = + reac_rate_local(627) 
  reac_source_local(60,627) = - reac_rate_local(627) 
  reac_source_local(50,628) = + reac_rate_local(628) 
  reac_source_local(52,628) = + reac_rate_local(628) 
  reac_source_local(55,628) = - reac_rate_local(628) 
  reac_source_local(60,628) = - reac_rate_local(628) 
  reac_source_local(51,629) = + reac_rate_local(629) 
  reac_source_local(52,629) = + reac_rate_local(629) 
  reac_source_local(56,629) = - reac_rate_local(629) 
  reac_source_local(60,629) = - reac_rate_local(629) 
  reac_source_local(52,630) = + reac_rate_local(630) * 2.d0
  reac_source_local(57,630) = - reac_rate_local(630) 
  reac_source_local(60,630) = - reac_rate_local(630) 
  reac_source_local(23,631) = + reac_rate_local(631) 
  reac_source_local(26,631) = - reac_rate_local(631) 
  reac_source_local(53,631) = + reac_rate_local(631) 
  reac_source_local(61,631) = - reac_rate_local(631) 
  reac_source_local(01,632) = + reac_rate_local(632) 
  reac_source_local(27,632) = - reac_rate_local(632) 
  reac_source_local(53,632) = + reac_rate_local(632) 
  reac_source_local(61,632) = - reac_rate_local(632) 
  reac_source_local(38,633) = + reac_rate_local(633) 
  reac_source_local(42,633) = - reac_rate_local(633) 
  reac_source_local(53,633) = + reac_rate_local(633) 
  reac_source_local(61,633) = - reac_rate_local(633) 
  reac_source_local(30,634) = + reac_rate_local(634) 
  reac_source_local(43,634) = - reac_rate_local(634) 
  reac_source_local(53,634) = + reac_rate_local(634) 
  reac_source_local(61,634) = - reac_rate_local(634) 
  reac_source_local(50,635) = + reac_rate_local(635) 
  reac_source_local(53,635) = + reac_rate_local(635) 
  reac_source_local(55,635) = - reac_rate_local(635) 
  reac_source_local(61,635) = - reac_rate_local(635) 
  reac_source_local(51,636) = + reac_rate_local(636) 
  reac_source_local(53,636) = + reac_rate_local(636) 
  reac_source_local(56,636) = - reac_rate_local(636) 
  reac_source_local(61,636) = - reac_rate_local(636) 
  reac_source_local(52,637) = + reac_rate_local(637) 
  reac_source_local(53,637) = + reac_rate_local(637) 
  reac_source_local(57,637) = - reac_rate_local(637) 
  reac_source_local(61,637) = - reac_rate_local(637) 
  reac_source_local(23,638) = + reac_rate_local(638) * 2.d0
  reac_source_local(27,638) = - reac_rate_local(638) 
  reac_source_local(38,638) = + reac_rate_local(638) 
  reac_source_local(46,638) = - reac_rate_local(638) 
  reac_source_local(01,639) = + reac_rate_local(639) 
  reac_source_local(23,639) = + reac_rate_local(639) 
  reac_source_local(28,639) = - reac_rate_local(639) 
  reac_source_local(38,639) = + reac_rate_local(639) 
  reac_source_local(46,639) = - reac_rate_local(639) 
  reac_source_local(01,640) = + reac_rate_local(640) * 2.d0
  reac_source_local(29,640) = - reac_rate_local(640) 
  reac_source_local(38,640) = + reac_rate_local(640) 
  reac_source_local(46,640) = - reac_rate_local(640) 
  reac_source_local(38,641) = + reac_rate_local(641) * 3.d0
  reac_source_local(43,641) = - reac_rate_local(641) 
  reac_source_local(46,641) = - reac_rate_local(641) 
  reac_source_local(30,642) = + reac_rate_local(642) * 2.d0
  reac_source_local(38,642) = + reac_rate_local(642) 
  reac_source_local(45,642) = - reac_rate_local(642) 
  reac_source_local(46,642) = - reac_rate_local(642) 
  reac_source_local(23,643) = + reac_rate_local(643) 
  reac_source_local(38,643) = + reac_rate_local(643) * 2.d0
  reac_source_local(46,643) = - reac_rate_local(643) 
  reac_source_local(55,643) = - reac_rate_local(643) 
  reac_source_local(01,644) = + reac_rate_local(644) 
  reac_source_local(38,644) = + reac_rate_local(644) * 2.d0
  reac_source_local(46,644) = - reac_rate_local(644) 
  reac_source_local(56,644) = - reac_rate_local(644) 
  reac_source_local(23,645) = + reac_rate_local(645) 
  reac_source_local(30,645) = + reac_rate_local(645) 
  reac_source_local(38,645) = + reac_rate_local(645) 
  reac_source_local(46,645) = - reac_rate_local(645) 
  reac_source_local(57,645) = - reac_rate_local(645) 
  reac_source_local(01,646) = + reac_rate_local(646) 
  reac_source_local(30,646) = + reac_rate_local(646) 
  reac_source_local(38,646) = + reac_rate_local(646) 
  reac_source_local(46,646) = - reac_rate_local(646) 
  reac_source_local(78,646) = - reac_rate_local(646) 
  reac_source_local(23,647) = + reac_rate_local(647) * 2.d0
  reac_source_local(27,647) = - reac_rate_local(647) 
  reac_source_local(30,647) = + reac_rate_local(647) 
  reac_source_local(47,647) = - reac_rate_local(647) 
  reac_source_local(01,648) = + reac_rate_local(648) 
  reac_source_local(23,648) = + reac_rate_local(648) 
  reac_source_local(28,648) = - reac_rate_local(648) 
  reac_source_local(30,648) = + reac_rate_local(648) 
  reac_source_local(47,648) = - reac_rate_local(648) 
  reac_source_local(01,649) = + reac_rate_local(649) * 2.d0
  reac_source_local(29,649) = - reac_rate_local(649) 
  reac_source_local(30,649) = + reac_rate_local(649) 
  reac_source_local(47,649) = - reac_rate_local(649) 
  reac_source_local(30,650) = + reac_rate_local(650) 
  reac_source_local(38,650) = + reac_rate_local(650) * 2.d0
  reac_source_local(43,650) = - reac_rate_local(650) 
  reac_source_local(47,650) = - reac_rate_local(650) 
  reac_source_local(30,651) = + reac_rate_local(651) * 3.d0
  reac_source_local(45,651) = - reac_rate_local(651) 
  reac_source_local(47,651) = - reac_rate_local(651) 
  reac_source_local(23,652) = + reac_rate_local(652) 
  reac_source_local(30,652) = + reac_rate_local(652) 
  reac_source_local(38,652) = + reac_rate_local(652) 
  reac_source_local(47,652) = - reac_rate_local(652) 
  reac_source_local(55,652) = - reac_rate_local(652) 
  reac_source_local(01,653) = + reac_rate_local(653) 
  reac_source_local(30,653) = + reac_rate_local(653) 
  reac_source_local(38,653) = + reac_rate_local(653) 
  reac_source_local(47,653) = - reac_rate_local(653) 
  reac_source_local(56,653) = - reac_rate_local(653) 
  reac_source_local(23,654) = + reac_rate_local(654) 
  reac_source_local(30,654) = + reac_rate_local(654) * 2.d0
  reac_source_local(47,654) = - reac_rate_local(654) 
  reac_source_local(57,654) = - reac_rate_local(654) 
  reac_source_local(01,655) = + reac_rate_local(655) 
  reac_source_local(30,655) = + reac_rate_local(655) * 2.d0
  reac_source_local(47,655) = - reac_rate_local(655) 
  reac_source_local(78,655) = - reac_rate_local(655) 
  reac_source_local(23,656) = + reac_rate_local(656) * 2.d0
  reac_source_local(27,656) = - reac_rate_local(656) 
  reac_source_local(41,656) = + reac_rate_local(656) 
  reac_source_local(48,656) = - reac_rate_local(656) 
  reac_source_local(01,657) = + reac_rate_local(657) 
  reac_source_local(23,657) = + reac_rate_local(657) 
  reac_source_local(28,657) = - reac_rate_local(657) 
  reac_source_local(41,657) = + reac_rate_local(657) 
  reac_source_local(48,657) = - reac_rate_local(657) 
  reac_source_local(01,658) = + reac_rate_local(658) * 2.d0
  reac_source_local(29,658) = - reac_rate_local(658) 
  reac_source_local(41,658) = + reac_rate_local(658) 
  reac_source_local(48,658) = - reac_rate_local(658) 
  reac_source_local(38,659) = + reac_rate_local(659) * 2.d0
  reac_source_local(41,659) = + reac_rate_local(659) 
  reac_source_local(43,659) = - reac_rate_local(659) 
  reac_source_local(48,659) = - reac_rate_local(659) 
  reac_source_local(30,660) = + reac_rate_local(660) * 2.d0
  reac_source_local(41,660) = + reac_rate_local(660) 
  reac_source_local(45,660) = - reac_rate_local(660) 
  reac_source_local(48,660) = - reac_rate_local(660) 
  reac_source_local(23,661) = + reac_rate_local(661) 
  reac_source_local(38,661) = + reac_rate_local(661) 
  reac_source_local(41,661) = + reac_rate_local(661) 
  reac_source_local(48,661) = - reac_rate_local(661) 
  reac_source_local(55,661) = - reac_rate_local(661) 
  reac_source_local(01,662) = + reac_rate_local(662) 
  reac_source_local(38,662) = + reac_rate_local(662) 
  reac_source_local(41,662) = + reac_rate_local(662) 
  reac_source_local(48,662) = - reac_rate_local(662) 
  reac_source_local(56,662) = - reac_rate_local(662) 
  reac_source_local(23,663) = + reac_rate_local(663) 
  reac_source_local(30,663) = + reac_rate_local(663) 
  reac_source_local(41,663) = + reac_rate_local(663) 
  reac_source_local(48,663) = - reac_rate_local(663) 
  reac_source_local(57,663) = - reac_rate_local(663) 
  reac_source_local(01,664) = + reac_rate_local(664) 
  reac_source_local(30,664) = + reac_rate_local(664) 
  reac_source_local(41,664) = + reac_rate_local(664) 
  reac_source_local(48,664) = - reac_rate_local(664) 
  reac_source_local(78,664) = - reac_rate_local(664) 
  reac_source_local(23,665) = + reac_rate_local(665) * 2.d0
  reac_source_local(27,665) = - reac_rate_local(665) 
  reac_source_local(50,665) = + reac_rate_local(665) 
  reac_source_local(58,665) = - reac_rate_local(665) 
  reac_source_local(01,666) = + reac_rate_local(666) 
  reac_source_local(23,666) = + reac_rate_local(666) 
  reac_source_local(28,666) = - reac_rate_local(666) 
  reac_source_local(50,666) = + reac_rate_local(666) 
  reac_source_local(58,666) = - reac_rate_local(666) 
  reac_source_local(01,667) = + reac_rate_local(667) * 2.d0
  reac_source_local(29,667) = - reac_rate_local(667) 
  reac_source_local(50,667) = + reac_rate_local(667) 
  reac_source_local(58,667) = - reac_rate_local(667) 
  reac_source_local(38,668) = + reac_rate_local(668) * 2.d0
  reac_source_local(43,668) = - reac_rate_local(668) 
  reac_source_local(50,668) = + reac_rate_local(668) 
  reac_source_local(58,668) = - reac_rate_local(668) 
  reac_source_local(30,669) = + reac_rate_local(669) * 2.d0
  reac_source_local(45,669) = - reac_rate_local(669) 
  reac_source_local(50,669) = + reac_rate_local(669) 
  reac_source_local(58,669) = - reac_rate_local(669) 
  reac_source_local(23,670) = + reac_rate_local(670) 
  reac_source_local(38,670) = + reac_rate_local(670) 
  reac_source_local(50,670) = + reac_rate_local(670) 
  reac_source_local(55,670) = - reac_rate_local(670) 
  reac_source_local(58,670) = - reac_rate_local(670) 
  reac_source_local(01,671) = + reac_rate_local(671) 
  reac_source_local(38,671) = + reac_rate_local(671) 
  reac_source_local(50,671) = + reac_rate_local(671) 
  reac_source_local(56,671) = - reac_rate_local(671) 
  reac_source_local(58,671) = - reac_rate_local(671) 
  reac_source_local(23,672) = + reac_rate_local(672) 
  reac_source_local(30,672) = + reac_rate_local(672) 
  reac_source_local(50,672) = + reac_rate_local(672) 
  reac_source_local(57,672) = - reac_rate_local(672) 
  reac_source_local(58,672) = - reac_rate_local(672) 
  reac_source_local(01,673) = + reac_rate_local(673) 
  reac_source_local(30,673) = + reac_rate_local(673) 
  reac_source_local(50,673) = + reac_rate_local(673) 
  reac_source_local(58,673) = - reac_rate_local(673) 
  reac_source_local(78,673) = - reac_rate_local(673) 
  reac_source_local(23,674) = + reac_rate_local(674) * 2.d0
  reac_source_local(27,674) = - reac_rate_local(674) 
  reac_source_local(51,674) = + reac_rate_local(674) 
  reac_source_local(59,674) = - reac_rate_local(674) 
  reac_source_local(01,675) = + reac_rate_local(675) 
  reac_source_local(23,675) = + reac_rate_local(675) 
  reac_source_local(28,675) = - reac_rate_local(675) 
  reac_source_local(51,675) = + reac_rate_local(675) 
  reac_source_local(59,675) = - reac_rate_local(675) 
  reac_source_local(01,676) = + reac_rate_local(676) * 2.d0
  reac_source_local(29,676) = - reac_rate_local(676) 
  reac_source_local(51,676) = + reac_rate_local(676) 
  reac_source_local(59,676) = - reac_rate_local(676) 
  reac_source_local(38,677) = + reac_rate_local(677) * 2.d0
  reac_source_local(43,677) = - reac_rate_local(677) 
  reac_source_local(51,677) = + reac_rate_local(677) 
  reac_source_local(59,677) = - reac_rate_local(677) 
  reac_source_local(30,678) = + reac_rate_local(678) * 2.d0
  reac_source_local(45,678) = - reac_rate_local(678) 
  reac_source_local(51,678) = + reac_rate_local(678) 
  reac_source_local(59,678) = - reac_rate_local(678) 
  reac_source_local(23,679) = + reac_rate_local(679) 
  reac_source_local(38,679) = + reac_rate_local(679) 
  reac_source_local(51,679) = + reac_rate_local(679) 
  reac_source_local(55,679) = - reac_rate_local(679) 
  reac_source_local(59,679) = - reac_rate_local(679) 
  reac_source_local(01,680) = + reac_rate_local(680) 
  reac_source_local(38,680) = + reac_rate_local(680) 
  reac_source_local(51,680) = + reac_rate_local(680) 
  reac_source_local(56,680) = - reac_rate_local(680) 
  reac_source_local(59,680) = - reac_rate_local(680) 
  reac_source_local(23,681) = + reac_rate_local(681) 
  reac_source_local(30,681) = + reac_rate_local(681) 
  reac_source_local(51,681) = + reac_rate_local(681) 
  reac_source_local(57,681) = - reac_rate_local(681) 
  reac_source_local(59,681) = - reac_rate_local(681) 
  reac_source_local(01,682) = + reac_rate_local(682) 
  reac_source_local(30,682) = + reac_rate_local(682) 
  reac_source_local(51,682) = + reac_rate_local(682) 
  reac_source_local(59,682) = - reac_rate_local(682) 
  reac_source_local(78,682) = - reac_rate_local(682) 
  reac_source_local(23,683) = + reac_rate_local(683) * 2.d0
  reac_source_local(27,683) = - reac_rate_local(683) 
  reac_source_local(52,683) = + reac_rate_local(683) 
  reac_source_local(60,683) = - reac_rate_local(683) 
  reac_source_local(01,684) = + reac_rate_local(684) 
  reac_source_local(23,684) = + reac_rate_local(684) 
  reac_source_local(28,684) = - reac_rate_local(684) 
  reac_source_local(52,684) = + reac_rate_local(684) 
  reac_source_local(60,684) = - reac_rate_local(684) 
  reac_source_local(01,685) = + reac_rate_local(685) * 2.d0
  reac_source_local(29,685) = - reac_rate_local(685) 
  reac_source_local(52,685) = + reac_rate_local(685) 
  reac_source_local(60,685) = - reac_rate_local(685) 
  reac_source_local(38,686) = + reac_rate_local(686) * 2.d0
  reac_source_local(43,686) = - reac_rate_local(686) 
  reac_source_local(52,686) = + reac_rate_local(686) 
  reac_source_local(60,686) = - reac_rate_local(686) 
  reac_source_local(30,687) = + reac_rate_local(687) * 2.d0
  reac_source_local(45,687) = - reac_rate_local(687) 
  reac_source_local(52,687) = + reac_rate_local(687) 
  reac_source_local(60,687) = - reac_rate_local(687) 
  reac_source_local(23,688) = + reac_rate_local(688) 
  reac_source_local(38,688) = + reac_rate_local(688) 
  reac_source_local(52,688) = + reac_rate_local(688) 
  reac_source_local(55,688) = - reac_rate_local(688) 
  reac_source_local(60,688) = - reac_rate_local(688) 
  reac_source_local(01,689) = + reac_rate_local(689) 
  reac_source_local(38,689) = + reac_rate_local(689) 
  reac_source_local(52,689) = + reac_rate_local(689) 
  reac_source_local(56,689) = - reac_rate_local(689) 
  reac_source_local(60,689) = - reac_rate_local(689) 
  reac_source_local(23,690) = + reac_rate_local(690) 
  reac_source_local(30,690) = + reac_rate_local(690) 
  reac_source_local(52,690) = + reac_rate_local(690) 
  reac_source_local(57,690) = - reac_rate_local(690) 
  reac_source_local(60,690) = - reac_rate_local(690) 
  reac_source_local(01,691) = + reac_rate_local(691) 
  reac_source_local(30,691) = + reac_rate_local(691) 
  reac_source_local(52,691) = + reac_rate_local(691) 
  reac_source_local(60,691) = - reac_rate_local(691) 
  reac_source_local(78,691) = - reac_rate_local(691) 
  reac_source_local(23,692) = + reac_rate_local(692) * 2.d0
  reac_source_local(27,692) = - reac_rate_local(692) 
  reac_source_local(53,692) = + reac_rate_local(692) 
  reac_source_local(61,692) = - reac_rate_local(692) 
  reac_source_local(01,693) = + reac_rate_local(693) 
  reac_source_local(23,693) = + reac_rate_local(693) 
  reac_source_local(28,693) = - reac_rate_local(693) 
  reac_source_local(53,693) = + reac_rate_local(693) 
  reac_source_local(61,693) = - reac_rate_local(693) 
  reac_source_local(01,694) = + reac_rate_local(694) * 2.d0
  reac_source_local(29,694) = - reac_rate_local(694) 
  reac_source_local(53,694) = + reac_rate_local(694) 
  reac_source_local(61,694) = - reac_rate_local(694) 
  reac_source_local(38,695) = + reac_rate_local(695) * 2.d0
  reac_source_local(43,695) = - reac_rate_local(695) 
  reac_source_local(53,695) = + reac_rate_local(695) 
  reac_source_local(61,695) = - reac_rate_local(695) 
  reac_source_local(30,696) = + reac_rate_local(696) * 2.d0
  reac_source_local(45,696) = - reac_rate_local(696) 
  reac_source_local(53,696) = + reac_rate_local(696) 
  reac_source_local(61,696) = - reac_rate_local(696) 
  reac_source_local(23,697) = + reac_rate_local(697) 
  reac_source_local(38,697) = + reac_rate_local(697) 
  reac_source_local(53,697) = + reac_rate_local(697) 
  reac_source_local(55,697) = - reac_rate_local(697) 
  reac_source_local(61,697) = - reac_rate_local(697) 
  reac_source_local(01,698) = + reac_rate_local(698) 
  reac_source_local(38,698) = + reac_rate_local(698) 
  reac_source_local(53,698) = + reac_rate_local(698) 
  reac_source_local(56,698) = - reac_rate_local(698) 
  reac_source_local(61,698) = - reac_rate_local(698) 
  reac_source_local(23,699) = + reac_rate_local(699) 
  reac_source_local(30,699) = + reac_rate_local(699) 
  reac_source_local(53,699) = + reac_rate_local(699) 
  reac_source_local(57,699) = - reac_rate_local(699) 
  reac_source_local(61,699) = - reac_rate_local(699) 
  reac_source_local(01,700) = + reac_rate_local(700) 
  reac_source_local(30,700) = + reac_rate_local(700) 
  reac_source_local(53,700) = + reac_rate_local(700) 
  reac_source_local(61,700) = - reac_rate_local(700) 
  reac_source_local(78,700) = - reac_rate_local(700) 
  reac_source_local(23,701) = + reac_rate_local(701) 
  reac_source_local(26,701) = - reac_rate_local(701) 
  reac_source_local(30,701) = + reac_rate_local(701) * 2.d0
  reac_source_local(49,701) = - reac_rate_local(701) 
  reac_source_local(01,702) = + reac_rate_local(702) 
  reac_source_local(27,702) = - reac_rate_local(702) 
  reac_source_local(30,702) = + reac_rate_local(702) * 2.d0
  reac_source_local(49,702) = - reac_rate_local(702) 
  reac_source_local(30,703) = + reac_rate_local(703) * 2.d0
  reac_source_local(38,703) = + reac_rate_local(703) 
  reac_source_local(42,703) = - reac_rate_local(703) 
  reac_source_local(49,703) = - reac_rate_local(703) 
  reac_source_local(30,704) = + reac_rate_local(704) * 3.d0
  reac_source_local(43,704) = - reac_rate_local(704) 
  reac_source_local(49,704) = - reac_rate_local(704) 
  reac_source_local(30,705) = + reac_rate_local(705) * 2.d0
  reac_source_local(49,705) = - reac_rate_local(705) 
  reac_source_local(50,705) = + reac_rate_local(705) 
  reac_source_local(55,705) = - reac_rate_local(705) 
  reac_source_local(30,706) = + reac_rate_local(706) * 2.d0
  reac_source_local(49,706) = - reac_rate_local(706) 
  reac_source_local(51,706) = + reac_rate_local(706) 
  reac_source_local(56,706) = - reac_rate_local(706) 
  reac_source_local(30,707) = + reac_rate_local(707) * 2.d0
  reac_source_local(49,707) = - reac_rate_local(707) 
  reac_source_local(52,707) = + reac_rate_local(707) 
  reac_source_local(57,707) = - reac_rate_local(707) 
  reac_source_local(01,708) = + reac_rate_local(708) 
  reac_source_local(23,708) = + reac_rate_local(708) 
  reac_source_local(28,708) = - reac_rate_local(708) 
  reac_source_local(30,708) = + reac_rate_local(708) * 2.d0
  reac_source_local(49,708) = - reac_rate_local(708) 
  reac_source_local(01,709) = + reac_rate_local(709) * 2.d0
  reac_source_local(29,709) = - reac_rate_local(709) 
  reac_source_local(30,709) = + reac_rate_local(709) * 2.d0
  reac_source_local(49,709) = - reac_rate_local(709) 
  reac_source_local(30,710) = + reac_rate_local(710) * 4.d0
  reac_source_local(45,710) = - reac_rate_local(710) 
  reac_source_local(49,710) = - reac_rate_local(710) 
  reac_source_local(01,711) = + reac_rate_local(711) 
  reac_source_local(30,711) = + reac_rate_local(711) * 3.d0
  reac_source_local(49,711) = - reac_rate_local(711) 
  reac_source_local(78,711) = - reac_rate_local(711) 
  reac_source_local(23,712) = + reac_rate_local(712) 
  reac_source_local(26,712) = - reac_rate_local(712) 
  reac_source_local(38,712) = + reac_rate_local(712) 
  reac_source_local(46,712) = - reac_rate_local(712) 
  reac_source_local(01,713) = + reac_rate_local(713) 
  reac_source_local(27,713) = - reac_rate_local(713) 
  reac_source_local(38,713) = + reac_rate_local(713) 
  reac_source_local(46,713) = - reac_rate_local(713) 
  reac_source_local(38,714) = + reac_rate_local(714) * 2.d0
  reac_source_local(42,714) = - reac_rate_local(714) 
  reac_source_local(46,714) = - reac_rate_local(714) 
  reac_source_local(30,715) = + reac_rate_local(715) 
  reac_source_local(38,715) = + reac_rate_local(715) 
  reac_source_local(43,715) = - reac_rate_local(715) 
  reac_source_local(46,715) = - reac_rate_local(715) 
  reac_source_local(38,716) = + reac_rate_local(716) 
  reac_source_local(46,716) = - reac_rate_local(716) 
  reac_source_local(50,716) = + reac_rate_local(716) 
  reac_source_local(55,716) = - reac_rate_local(716) 
  reac_source_local(23,717) = + reac_rate_local(717) 
  reac_source_local(26,717) = - reac_rate_local(717) 
  reac_source_local(30,717) = + reac_rate_local(717) 
  reac_source_local(47,717) = - reac_rate_local(717) 
  reac_source_local(01,718) = + reac_rate_local(718) 
  reac_source_local(27,718) = - reac_rate_local(718) 
  reac_source_local(30,718) = + reac_rate_local(718) 
  reac_source_local(47,718) = - reac_rate_local(718) 
  reac_source_local(30,719) = + reac_rate_local(719) 
  reac_source_local(38,719) = + reac_rate_local(719) 
  reac_source_local(42,719) = - reac_rate_local(719) 
  reac_source_local(47,719) = - reac_rate_local(719) 
  reac_source_local(30,720) = + reac_rate_local(720) * 2.d0
  reac_source_local(43,720) = - reac_rate_local(720) 
  reac_source_local(47,720) = - reac_rate_local(720) 
  reac_source_local(30,721) = + reac_rate_local(721) 
  reac_source_local(47,721) = - reac_rate_local(721) 
  reac_source_local(50,721) = + reac_rate_local(721) 
  reac_source_local(55,721) = - reac_rate_local(721) 
  reac_source_local(26,722) = - reac_rate_local(722) 
  reac_source_local(46,722) = - reac_rate_local(722) 
  reac_source_local(50,722) = + reac_rate_local(722) 
  reac_source_local(27,723) = - reac_rate_local(723) 
  reac_source_local(46,723) = - reac_rate_local(723) 
  reac_source_local(51,723) = + reac_rate_local(723) 
  reac_source_local(30,724) = + reac_rate_local(724) 
  reac_source_local(42,724) = - reac_rate_local(724) 
  reac_source_local(46,724) = - reac_rate_local(724) 
  reac_source_local(41,725) = + reac_rate_local(725) 
  reac_source_local(43,725) = - reac_rate_local(725) 
  reac_source_local(46,725) = - reac_rate_local(725) 
  reac_source_local(46,726) = - reac_rate_local(726) 
  reac_source_local(52,726) = + reac_rate_local(726) 
  reac_source_local(55,726) = - reac_rate_local(726) 
  reac_source_local(26,727) = - reac_rate_local(727) 
  reac_source_local(47,727) = - reac_rate_local(727) 
  reac_source_local(52,727) = + reac_rate_local(727) 
  reac_source_local(41,728) = + reac_rate_local(728) 
  reac_source_local(42,728) = - reac_rate_local(728) 
  reac_source_local(47,728) = - reac_rate_local(728) 
  reac_source_local(47,729) = - reac_rate_local(729) 
  reac_source_local(53,729) = + reac_rate_local(729) 
  reac_source_local(55,729) = - reac_rate_local(729) 
  reac_source_local(30,730) = + reac_rate_local(730) 
  reac_source_local(47,730) = - reac_rate_local(730) 
  reac_source_local(70,730) = + reac_rate_local(730) 
  reac_source_local(73,730) = - reac_rate_local(730) 
  reac_source_local(38,731) = + reac_rate_local(731) 
  reac_source_local(46,731) = - reac_rate_local(731) 
  reac_source_local(70,731) = + reac_rate_local(731) 
  reac_source_local(73,731) = - reac_rate_local(731) 
  reac_source_local(41,732) = + reac_rate_local(732) 
  reac_source_local(48,732) = - reac_rate_local(732) 
  reac_source_local(70,732) = + reac_rate_local(732) 
  reac_source_local(73,732) = - reac_rate_local(732) 
  reac_source_local(50,733) = + reac_rate_local(733) 
  reac_source_local(58,733) = - reac_rate_local(733) 
  reac_source_local(70,733) = + reac_rate_local(733) 
  reac_source_local(73,733) = - reac_rate_local(733) 
  reac_source_local(52,734) = + reac_rate_local(734) 
  reac_source_local(60,734) = - reac_rate_local(734) 
  reac_source_local(70,734) = + reac_rate_local(734) 
  reac_source_local(73,734) = - reac_rate_local(734) 
  reac_source_local(53,735) = + reac_rate_local(735) 
  reac_source_local(61,735) = - reac_rate_local(735) 
  reac_source_local(70,735) = + reac_rate_local(735) 
  reac_source_local(73,735) = - reac_rate_local(735) 
  reac_source_local(51,736) = + reac_rate_local(736) 
  reac_source_local(59,736) = - reac_rate_local(736) 
  reac_source_local(70,736) = + reac_rate_local(736) 
  reac_source_local(73,736) = - reac_rate_local(736) 
  reac_source_local(53,737) = + reac_rate_local(737) 
  reac_source_local(61,737) = - reac_rate_local(737) 
  reac_source_local(63,737) = + reac_rate_local(737) 
  reac_source_local(70,737) = + reac_rate_local(737) 
  reac_source_local(77,737) = - reac_rate_local(737) 
  reac_source_local(63,738) = + reac_rate_local(738) * 3.d0
  reac_source_local(65,738) = - reac_rate_local(738) 
  reac_source_local(67,738) = - reac_rate_local(738) 
  reac_source_local(62,739) = + reac_rate_local(739) 
  reac_source_local(63,739) = + reac_rate_local(739) * 2.d0
  reac_source_local(66,739) = - reac_rate_local(739) 
  reac_source_local(67,739) = - reac_rate_local(739) 
  reac_source_local(01,740) = + reac_rate_local(740) 
  reac_source_local(27,740) = - reac_rate_local(740) 
  reac_source_local(63,740) = + reac_rate_local(740) 
  reac_source_local(67,740) = - reac_rate_local(740) 
  reac_source_local(01,741) = + reac_rate_local(741) * 2.d0
  reac_source_local(29,741) = - reac_rate_local(741) 
  reac_source_local(63,741) = + reac_rate_local(741) 
  reac_source_local(67,741) = - reac_rate_local(741) 
  reac_source_local(01,742) = + reac_rate_local(742) 
  reac_source_local(62,742) = + reac_rate_local(742) 
  reac_source_local(67,742) = - reac_rate_local(742) 
  reac_source_local(76,742) = - reac_rate_local(742) 
  reac_source_local(62,743) = + reac_rate_local(743) 
  reac_source_local(63,743) = + reac_rate_local(743) 
  reac_source_local(65,743) = - reac_rate_local(743) 
  reac_source_local(67,743) = - reac_rate_local(743) 
  reac_source_local(62,744) = + reac_rate_local(744) 
  reac_source_local(63,744) = + reac_rate_local(744) 
  reac_source_local(65,744) = - reac_rate_local(744) 
  reac_source_local(67,744) = - reac_rate_local(744) 
  reac_source_local(62,745) = + reac_rate_local(745) 
  reac_source_local(63,745) = + reac_rate_local(745) 
  reac_source_local(65,745) = - reac_rate_local(745) 
  reac_source_local(67,745) = - reac_rate_local(745) 
  reac_source_local(62,746) = + reac_rate_local(746) * 2.d0
  reac_source_local(66,746) = - reac_rate_local(746) 
  reac_source_local(67,746) = - reac_rate_local(746) 
  reac_source_local(62,747) = + reac_rate_local(747) * 2.d0
  reac_source_local(66,747) = - reac_rate_local(747) 
  reac_source_local(67,747) = - reac_rate_local(747) 
  reac_source_local(62,748) = + reac_rate_local(748) * 2.d0
  reac_source_local(66,748) = - reac_rate_local(748) 
  reac_source_local(67,748) = - reac_rate_local(748) 
  reac_source_local(01,749) = + reac_rate_local(749) 
  reac_source_local(27,749) = - reac_rate_local(749) 
  reac_source_local(63,749) = + reac_rate_local(749) 
  reac_source_local(67,749) = - reac_rate_local(749) 
  reac_source_local(01,750) = + reac_rate_local(750) 
  reac_source_local(27,750) = - reac_rate_local(750) 
  reac_source_local(63,750) = + reac_rate_local(750) 
  reac_source_local(67,750) = - reac_rate_local(750) 
  reac_source_local(01,751) = + reac_rate_local(751) 
  reac_source_local(27,751) = - reac_rate_local(751) 
  reac_source_local(63,751) = + reac_rate_local(751) 
  reac_source_local(67,751) = - reac_rate_local(751) 
  reac_source_local(01,752) = + reac_rate_local(752) * 2.d0
  reac_source_local(29,752) = - reac_rate_local(752) 
  reac_source_local(63,752) = + reac_rate_local(752) 
  reac_source_local(67,752) = - reac_rate_local(752) 
  reac_source_local(01,753) = + reac_rate_local(753) * 2.d0
  reac_source_local(29,753) = - reac_rate_local(753) 
  reac_source_local(63,753) = + reac_rate_local(753) 
  reac_source_local(67,753) = - reac_rate_local(753) 
  reac_source_local(01,754) = + reac_rate_local(754) * 2.d0
  reac_source_local(29,754) = - reac_rate_local(754) 
  reac_source_local(63,754) = + reac_rate_local(754) 
  reac_source_local(67,754) = - reac_rate_local(754) 
  reac_source_local(01,755) = + reac_rate_local(755) 
  reac_source_local(62,755) = + reac_rate_local(755) 
  reac_source_local(67,755) = - reac_rate_local(755) 
  reac_source_local(76,755) = - reac_rate_local(755) 
  reac_source_local(01,756) = + reac_rate_local(756) 
  reac_source_local(62,756) = + reac_rate_local(756) 
  reac_source_local(67,756) = - reac_rate_local(756) 
  reac_source_local(76,756) = - reac_rate_local(756) 
  reac_source_local(01,757) = + reac_rate_local(757) 
  reac_source_local(62,757) = + reac_rate_local(757) 
  reac_source_local(67,757) = - reac_rate_local(757) 
  reac_source_local(76,757) = - reac_rate_local(757) 
  reac_source_local(23,758) = + reac_rate_local(758) 
  reac_source_local(26,758) = - reac_rate_local(758) 
  reac_source_local(41,758) = + reac_rate_local(758) 
  reac_source_local(48,758) = - reac_rate_local(758) 
  reac_source_local(01,759) = + reac_rate_local(759) 
  reac_source_local(27,759) = - reac_rate_local(759) 
  reac_source_local(41,759) = + reac_rate_local(759) 
  reac_source_local(48,759) = - reac_rate_local(759) 
  reac_source_local(38,760) = + reac_rate_local(760) 
  reac_source_local(41,760) = + reac_rate_local(760) 
  reac_source_local(42,760) = - reac_rate_local(760) 
  reac_source_local(48,760) = - reac_rate_local(760) 
  reac_source_local(30,761) = + reac_rate_local(761) 
  reac_source_local(41,761) = + reac_rate_local(761) 
  reac_source_local(43,761) = - reac_rate_local(761) 
  reac_source_local(48,761) = - reac_rate_local(761) 
  reac_source_local(41,762) = + reac_rate_local(762) 
  reac_source_local(48,762) = - reac_rate_local(762) 
  reac_source_local(50,762) = + reac_rate_local(762) 
  reac_source_local(55,762) = - reac_rate_local(762) 
  reac_source_local(41,763) = + reac_rate_local(763) 
  reac_source_local(48,763) = - reac_rate_local(763) 
  reac_source_local(51,763) = + reac_rate_local(763) 
  reac_source_local(56,763) = - reac_rate_local(763) 
  reac_source_local(41,764) = + reac_rate_local(764) 
  reac_source_local(48,764) = - reac_rate_local(764) 
  reac_source_local(52,764) = + reac_rate_local(764) 
  reac_source_local(57,764) = - reac_rate_local(764) 
  reac_source_local(23,765) = + reac_rate_local(765) 
  reac_source_local(26,765) = - reac_rate_local(765) 
  reac_source_local(50,765) = + reac_rate_local(765) 
  reac_source_local(58,765) = - reac_rate_local(765) 
  reac_source_local(01,766) = + reac_rate_local(766) 
  reac_source_local(27,766) = - reac_rate_local(766) 
  reac_source_local(50,766) = + reac_rate_local(766) 
  reac_source_local(58,766) = - reac_rate_local(766) 
  reac_source_local(38,767) = + reac_rate_local(767) 
  reac_source_local(42,767) = - reac_rate_local(767) 
  reac_source_local(50,767) = + reac_rate_local(767) 
  reac_source_local(58,767) = - reac_rate_local(767) 
  reac_source_local(30,768) = + reac_rate_local(768) 
  reac_source_local(43,768) = - reac_rate_local(768) 
  reac_source_local(50,768) = + reac_rate_local(768) 
  reac_source_local(58,768) = - reac_rate_local(768) 
  reac_source_local(50,769) = + reac_rate_local(769) * 2.d0
  reac_source_local(55,769) = - reac_rate_local(769) 
  reac_source_local(58,769) = - reac_rate_local(769) 
  reac_source_local(50,770) = + reac_rate_local(770) 
  reac_source_local(51,770) = + reac_rate_local(770) 
  reac_source_local(56,770) = - reac_rate_local(770) 
  reac_source_local(58,770) = - reac_rate_local(770) 
  reac_source_local(50,771) = + reac_rate_local(771) 
  reac_source_local(52,771) = + reac_rate_local(771) 
  reac_source_local(57,771) = - reac_rate_local(771) 
  reac_source_local(58,771) = - reac_rate_local(771) 
  reac_source_local(23,772) = + reac_rate_local(772) 
  reac_source_local(26,772) = - reac_rate_local(772) 
  reac_source_local(51,772) = + reac_rate_local(772) 
  reac_source_local(59,772) = - reac_rate_local(772) 
  reac_source_local(01,773) = + reac_rate_local(773) 
  reac_source_local(27,773) = - reac_rate_local(773) 
  reac_source_local(51,773) = + reac_rate_local(773) 
  reac_source_local(59,773) = - reac_rate_local(773) 
  reac_source_local(38,774) = + reac_rate_local(774) 
  reac_source_local(42,774) = - reac_rate_local(774) 
  reac_source_local(51,774) = + reac_rate_local(774) 
  reac_source_local(59,774) = - reac_rate_local(774) 
  reac_source_local(30,775) = + reac_rate_local(775) 
  reac_source_local(43,775) = - reac_rate_local(775) 
  reac_source_local(51,775) = + reac_rate_local(775) 
  reac_source_local(59,775) = - reac_rate_local(775) 
  reac_source_local(50,776) = + reac_rate_local(776) 
  reac_source_local(51,776) = + reac_rate_local(776) 
  reac_source_local(55,776) = - reac_rate_local(776) 
  reac_source_local(59,776) = - reac_rate_local(776) 
  reac_source_local(51,777) = + reac_rate_local(777) * 2.d0
  reac_source_local(56,777) = - reac_rate_local(777) 
  reac_source_local(59,777) = - reac_rate_local(777) 
  reac_source_local(51,778) = + reac_rate_local(778) 
  reac_source_local(52,778) = + reac_rate_local(778) 
  reac_source_local(57,778) = - reac_rate_local(778) 
  reac_source_local(59,778) = - reac_rate_local(778) 
  reac_source_local(23,779) = + reac_rate_local(779) 
  reac_source_local(26,779) = - reac_rate_local(779) 
  reac_source_local(52,779) = + reac_rate_local(779) 
  reac_source_local(60,779) = - reac_rate_local(779) 
  reac_source_local(01,780) = + reac_rate_local(780) 
  reac_source_local(27,780) = - reac_rate_local(780) 
  reac_source_local(52,780) = + reac_rate_local(780) 
  reac_source_local(60,780) = - reac_rate_local(780) 
  reac_source_local(38,781) = + reac_rate_local(781) 
  reac_source_local(42,781) = - reac_rate_local(781) 
  reac_source_local(52,781) = + reac_rate_local(781) 
  reac_source_local(60,781) = - reac_rate_local(781) 
  reac_source_local(30,782) = + reac_rate_local(782) 
  reac_source_local(43,782) = - reac_rate_local(782) 
  reac_source_local(52,782) = + reac_rate_local(782) 
  reac_source_local(60,782) = - reac_rate_local(782) 
  reac_source_local(50,783) = + reac_rate_local(783) 
  reac_source_local(52,783) = + reac_rate_local(783) 
  reac_source_local(55,783) = - reac_rate_local(783) 
  reac_source_local(60,783) = - reac_rate_local(783) 
  reac_source_local(51,784) = + reac_rate_local(784) 
  reac_source_local(52,784) = + reac_rate_local(784) 
  reac_source_local(56,784) = - reac_rate_local(784) 
  reac_source_local(60,784) = - reac_rate_local(784) 
  reac_source_local(52,785) = + reac_rate_local(785) * 2.d0
  reac_source_local(57,785) = - reac_rate_local(785) 
  reac_source_local(60,785) = - reac_rate_local(785) 
  reac_source_local(23,786) = + reac_rate_local(786) 
  reac_source_local(26,786) = - reac_rate_local(786) 
  reac_source_local(53,786) = + reac_rate_local(786) 
  reac_source_local(61,786) = - reac_rate_local(786) 
  reac_source_local(01,787) = + reac_rate_local(787) 
  reac_source_local(27,787) = - reac_rate_local(787) 
  reac_source_local(53,787) = + reac_rate_local(787) 
  reac_source_local(61,787) = - reac_rate_local(787) 
  reac_source_local(38,788) = + reac_rate_local(788) 
  reac_source_local(42,788) = - reac_rate_local(788) 
  reac_source_local(53,788) = + reac_rate_local(788) 
  reac_source_local(61,788) = - reac_rate_local(788) 
  reac_source_local(30,789) = + reac_rate_local(789) 
  reac_source_local(43,789) = - reac_rate_local(789) 
  reac_source_local(53,789) = + reac_rate_local(789) 
  reac_source_local(61,789) = - reac_rate_local(789) 
  reac_source_local(50,790) = + reac_rate_local(790) 
  reac_source_local(53,790) = + reac_rate_local(790) 
  reac_source_local(55,790) = - reac_rate_local(790) 
  reac_source_local(61,790) = - reac_rate_local(790) 
  reac_source_local(51,791) = + reac_rate_local(791) 
  reac_source_local(53,791) = + reac_rate_local(791) 
  reac_source_local(56,791) = - reac_rate_local(791) 
  reac_source_local(61,791) = - reac_rate_local(791) 
  reac_source_local(52,792) = + reac_rate_local(792) 
  reac_source_local(53,792) = + reac_rate_local(792) 
  reac_source_local(57,792) = - reac_rate_local(792) 
  reac_source_local(61,792) = - reac_rate_local(792) 
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(80)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(79) 
  rrt(002) = rrt(002) * density(01) * density(79) 
  rrt(003) = rrt(003) * density(01) * density(79) 
  rrt(004) = rrt(004) * density(01) * density(79) 
  rrt(005) = rrt(005) * density(01) * density(79) 
  rrt(006) = rrt(006) * density(01) * density(79) 
  rrt(007) = rrt(007) * density(01) * density(79) 
  rrt(008) = rrt(008) * density(01) * density(79) 
  rrt(009) = rrt(009) * density(01) * density(79) 
  rrt(010) = rrt(010) * density(01) * density(79) 
  rrt(011) = rrt(011) * density(01) * density(79) 
  rrt(012) = rrt(012) * density(01) * density(79) 
  rrt(013) = rrt(013) * density(01) * density(79) 
  rrt(014) = rrt(014) * density(01) * density(79) 
  rrt(015) = rrt(015) * density(01) * density(79) 
  rrt(016) = rrt(016) * density(02) * density(79) 
  rrt(017) = rrt(017) * density(02) * density(79) 
  rrt(018) = rrt(018) * density(02) * density(79) 
  rrt(019) = rrt(019) * density(02) * density(79) 
  rrt(020) = rrt(020) * density(03) * density(79) 
  rrt(021) = rrt(021) * density(03) * density(79) 
  rrt(022) = rrt(022) * density(03) * density(79) 
  rrt(023) = rrt(023) * density(03) * density(79) 
  rrt(024) = rrt(024) * density(03) * density(79) 
  rrt(025) = rrt(025) * density(04) * density(79) 
  rrt(026) = rrt(026) * density(04) * density(79) 
  rrt(027) = rrt(027) * density(04) * density(79) 
  rrt(028) = rrt(028) * density(04) * density(79) 
  rrt(029) = rrt(029) * density(04) * density(79) 
  rrt(030) = rrt(030) * density(04) * density(79) 
  rrt(031) = rrt(031) * density(05) * density(79) 
  rrt(032) = rrt(032) * density(05) * density(79) 
  rrt(033) = rrt(033) * density(05) * density(79) 
  rrt(034) = rrt(034) * density(05) * density(79) 
  rrt(035) = rrt(035) * density(06) * density(79) 
  rrt(036) = rrt(036) * density(06) * density(79) 
  rrt(037) = rrt(037) * density(06) * density(79) 
  rrt(038) = rrt(038) * density(07) * density(79) 
  rrt(039) = rrt(039) * density(07) * density(79) 
  rrt(040) = rrt(040) * density(08) * density(79) 
  rrt(041) = rrt(041) * density(09) * density(79) 
  rrt(042) = rrt(042) * density(10) * density(79) 
  rrt(043) = rrt(043) * density(11) * density(79) 
  rrt(044) = rrt(044) * density(30) * density(79) 
  rrt(045) = rrt(045) * density(30) * density(79) 
  rrt(046) = rrt(046) * density(30) * density(79) 
  rrt(047) = rrt(047) * density(30) * density(79) 
  rrt(048) = rrt(048) * density(30) * density(79) 
  rrt(049) = rrt(049) * density(30) * density(79) 
  rrt(050) = rrt(050) * density(01) * density(02) 
  rrt(051) = rrt(051) * density(01) * density(03) 
  rrt(052) = rrt(052) * density(01) * density(04) 
  rrt(053) = rrt(053) * density(01) * density(05) 
  rrt(054) = rrt(054) * density(01) * density(06) 
  rrt(055) = rrt(055) * density(01) * density(07) 
  rrt(056) = rrt(056) * density(01) * density(08) 
  rrt(057) = rrt(057) * density(01) * density(09) 
  rrt(058) = rrt(058) * density(01)**2 
  rrt(059) = rrt(059) * density(01) * density(02) 
  rrt(060) = rrt(060) * density(01) * density(03) 
  rrt(061) = rrt(061) * density(01) * density(04) 
  rrt(062) = rrt(062) * density(01) * density(05) 
  rrt(063) = rrt(063) * density(01) * density(06) 
  rrt(064) = rrt(064) * density(01) * density(07) 
  rrt(065) = rrt(065) * density(01) * density(08) 
  rrt(066) = rrt(066) * density(02) * density(23) 
  rrt(067) = rrt(067) * density(03) * density(23) 
  rrt(068) = rrt(068) * density(04) * density(23) 
  rrt(069) = rrt(069) * density(05) * density(23) 
  rrt(070) = rrt(070) * density(06) * density(23) 
  rrt(071) = rrt(071) * density(07) * density(23) 
  rrt(072) = rrt(072) * density(08) * density(23) 
  rrt(073) = rrt(073) * density(09) * density(23) 
  rrt(074) = rrt(074) * density(01) * density(23) 
  rrt(075) = rrt(075) * density(02) * density(23) 
  rrt(076) = rrt(076) * density(03) * density(23) 
  rrt(077) = rrt(077) * density(04) * density(23) 
  rrt(078) = rrt(078) * density(05) * density(23) 
  rrt(079) = rrt(079) * density(06) * density(23) 
  rrt(080) = rrt(080) * density(07) * density(23) 
  rrt(081) = rrt(081) * density(08) * density(23) 
  rrt(082) = rrt(082) * density(02) * density(38) 
  rrt(083) = rrt(083) * density(03) * density(38) 
  rrt(084) = rrt(084) * density(04) * density(38) 
  rrt(085) = rrt(085) * density(05) * density(38) 
  rrt(086) = rrt(086) * density(06) * density(38) 
  rrt(087) = rrt(087) * density(07) * density(38) 
  rrt(088) = rrt(088) * density(08) * density(38) 
  rrt(089) = rrt(089) * density(09) * density(38) 
  rrt(090) = rrt(090) * density(01) * density(38) 
  rrt(091) = rrt(091) * density(02) * density(38) 
  rrt(092) = rrt(092) * density(03) * density(38) 
  rrt(093) = rrt(093) * density(04) * density(38) 
  rrt(094) = rrt(094) * density(05) * density(38) 
  rrt(095) = rrt(095) * density(06) * density(38) 
  rrt(096) = rrt(096) * density(07) * density(38) 
  rrt(097) = rrt(097) * density(08) * density(38) 
  rrt(098) = rrt(098) * density(30) * density(31) 
  rrt(099) = rrt(099) * density(30) * density(32) 
  rrt(100) = rrt(100) * density(30) * density(33) 
  rrt(101) = rrt(101) * density(30) * density(34) 
  rrt(102) = rrt(102) * density(30)**2 
  rrt(103) = rrt(103) * density(30) * density(31) 
  rrt(104) = rrt(104) * density(30) * density(32) 
  rrt(105) = rrt(105) * density(30) * density(33) 
  rrt(106) = rrt(106) * density(31) * density(38) 
  rrt(107) = rrt(107) * density(32) * density(38) 
  rrt(108) = rrt(108) * density(33) * density(38) 
  rrt(109) = rrt(109) * density(34) * density(38) 
  rrt(110) = rrt(110) * density(30) * density(38) 
  rrt(111) = rrt(111) * density(31) * density(38) 
  rrt(112) = rrt(112) * density(32) * density(38) 
  rrt(113) = rrt(113) * density(33) * density(38) 
  rrt(114) = rrt(114) * density(01) * density(79) 
  rrt(115) = rrt(115) * density(01) * density(79) 
  rrt(116) = rrt(116) * density(01) * density(79) 
  rrt(117) = rrt(117) * density(01) * density(79) 
  rrt(118) = rrt(118) * density(01) * density(79) 
  rrt(119) = rrt(119) * density(01) * density(79) 
  rrt(120) = rrt(120) * density(01) * density(79) 
  rrt(121) = rrt(121) * density(01) * density(79) 
  rrt(122) = rrt(122) * density(01) * density(79) 
  rrt(123) = rrt(123) * density(01) * density(79) 
  rrt(124) = rrt(124) * density(01) * density(79) 
  rrt(125) = rrt(125) * density(01) * density(79) 
  rrt(126) = rrt(126) * density(01) * density(79) 
  rrt(127) = rrt(127) * density(01) * density(79) 
  rrt(128) = rrt(128) * density(01) * density(79) 
  rrt(129) = rrt(129) * density(01) * density(79) 
  rrt(130) = rrt(130) * density(01) * density(79) 
  rrt(131) = rrt(131) * density(01) * density(79) 
  rrt(132) = rrt(132) * density(01) * density(79) 
  rrt(133) = rrt(133) * density(01) * density(79) 
  rrt(134) = rrt(134) * density(01) * density(79) 
  rrt(135) = rrt(135) * density(01) * density(79) 
  rrt(136) = rrt(136) * density(01) * density(79) 
  rrt(137) = rrt(137) * density(01) * density(79) 
  rrt(138) = rrt(138) * density(01) * density(79) 
  rrt(139) = rrt(139) * density(01) * density(79) 
  rrt(140) = rrt(140) * density(01) * density(79) 
  rrt(141) = rrt(141) * density(01) * density(79) 
  rrt(142) = rrt(142) * density(01) * density(79) 
  rrt(143) = rrt(143) * density(30) * density(79) 
  rrt(144) = rrt(144) * density(30) * density(79) 
  rrt(145) = rrt(145) * density(30) * density(79) 
  rrt(146) = rrt(146) * density(30) * density(79) 
  rrt(147) = rrt(147) * density(30) * density(79) 
  rrt(148) = rrt(148) * density(30) * density(79) 
  rrt(149) = rrt(149) * density(35) * density(79) 
  rrt(150) = rrt(150) * density(38) * density(79) 
  rrt(151) = rrt(151) * density(38) * density(79) 
  rrt(152) = rrt(152) * density(70) * density(79) 
  rrt(153) = rrt(153) * density(70) * density(79) 
  rrt(154) = rrt(154) * density(23) * density(79) 
  rrt(155) = rrt(155) * density(38) * density(79) 
  rrt(156) = rrt(156) * density(01) * density(79) 
  rrt(157) = rrt(157) * density(17) * density(79) 
  rrt(158) = rrt(158) * density(30) * density(79) 
  rrt(159) = rrt(159) * density(35) * density(79) 
  rrt(160) = rrt(160) * density(50) * density(79) 
  rrt(161) = rrt(161) * density(51) * density(79) 
  rrt(162) = rrt(162) * density(41) * density(79) 
  rrt(163) = rrt(163) * density(70) * density(79) 
  rrt(164) = rrt(164) * density(70) * density(79) 
  rrt(165) = rrt(165) * density(70) * density(79) 
  rrt(166) = rrt(166) * density(27) * density(79) 
  rrt(167) = rrt(167) * density(27) * density(79) 
  rrt(168) = rrt(168) * density(27) * density(79) 
  rrt(169) = rrt(169) * density(43) * density(79) 
  rrt(170) = rrt(170) * density(43) * density(79) 
  rrt(171) = rrt(171) * density(43) * density(79) 
  rrt(172) = rrt(172) * density(55) * density(79) 
  rrt(173) = rrt(173) * density(55) * density(79) 
  rrt(174) = rrt(174) * density(28) * density(79) 
  rrt(175) = rrt(175) * density(29) * density(79) 
  rrt(176) = rrt(176) * density(56) * density(79) 
  rrt(177) = rrt(177) * density(57) * density(79) 
  rrt(178) = rrt(178) * density(45) * density(79) 
  rrt(179) = rrt(179) * density(78) * density(79) 
  rrt(180) = rrt(180) * density(26) * density(79)**2 
  rrt(181) = rrt(181) * density(42) * density(79)**2 
  rrt(182) = rrt(182) * density(26) * density(79) 
  rrt(183) = rrt(183) * density(42) * density(79) 
  rrt(184) = rrt(184) * density(77) * density(79) 
  rrt(185) = rrt(185) * density(30) * density(79) 
  rrt(186) = rrt(186) * density(50) * density(79) 
  rrt(187) = rrt(187) * density(41) * density(79) 
  rrt(188) = rrt(188) * density(41) * density(79) 
  rrt(189) = rrt(189) * density(52) * density(79) 
  rrt(190) = rrt(190) * density(30) * density(38) * density(79) 
  rrt(191) = rrt(191) * density(30) * density(38) * density(79) 
  rrt(192) = rrt(192) * density(41) * density(79) 
  rrt(193) = rrt(193) * density(50) * density(79) 
  rrt(194) = rrt(194) * density(51) * density(79) 
  rrt(195) = rrt(195) * density(01) * density(30) * density(79) 
  rrt(196) = rrt(196) * density(70) * density(79) 
  rrt(197) = rrt(197) * density(38) * density(46) 
  rrt(198) = rrt(198) * density(23) * density(46) 
  rrt(199) = rrt(199) * density(46) * density(50) 
  rrt(200) = rrt(200) * density(01) * density(46) 
  rrt(201) = rrt(201) * density(30) * density(46) 
  rrt(202) = rrt(202) * density(35) * density(46) 
  rrt(203) = rrt(203) * density(36) * density(46) 
  rrt(204) = rrt(204) * density(17) * density(46) 
  rrt(205) = rrt(205) * density(18) * density(46) 
  rrt(206) = rrt(206) * density(41) * density(46) 
  rrt(207) = rrt(207) * density(38) * density(47) 
  rrt(208) = rrt(208) * density(23) * density(47) 
  rrt(209) = rrt(209) * density(30) * density(47) 
  rrt(210) = rrt(210) * density(35) * density(47) 
  rrt(211) = rrt(211) * density(36) * density(47) 
  rrt(212) = rrt(212) * density(01) * density(47) 
  rrt(213) = rrt(213) * density(17) * density(47) 
  rrt(214) = rrt(214) * density(18) * density(47) 
  rrt(215) = rrt(215) * density(38) * density(48) 
  rrt(216) = rrt(216) * density(23) * density(58) 
  rrt(217) = rrt(217) * density(23) * density(48) 
  rrt(218) = rrt(218) * density(23) * density(59) 
  rrt(219) = rrt(219) * density(23) * density(60) 
  rrt(220) = rrt(220) * density(23) * density(61) 
  rrt(221) = rrt(221) * density(38) * density(58) 
  rrt(222) = rrt(222) * density(38) * density(59) 
  rrt(223) = rrt(223) * density(38) * density(60) 
  rrt(224) = rrt(224) * density(38) * density(61) 
  rrt(225) = rrt(225) * density(17) * density(48) 
  rrt(226) = rrt(226) * density(17) * density(58) 
  rrt(227) = rrt(227) * density(17) * density(59) 
  rrt(228) = rrt(228) * density(17) * density(60) 
  rrt(229) = rrt(229) * density(17) * density(61) 
  rrt(230) = rrt(230) * density(18) * density(48) 
  rrt(231) = rrt(231) * density(18) * density(58) 
  rrt(232) = rrt(232) * density(18) * density(59) 
  rrt(233) = rrt(233) * density(18) * density(60) 
  rrt(234) = rrt(234) * density(18) * density(61) 
  rrt(235) = rrt(235) * density(17) 
  rrt(236) = rrt(236) * density(18) 
  rrt(237) = rrt(237) * density(19) 
  rrt(238) = rrt(238) * density(20) 
  rrt(239) = rrt(239) * density(35) 
  rrt(240) = rrt(240) * density(36) 
  rrt(241) = rrt(241) * density(36) 
  rrt(242) = rrt(242) * density(37) 
  rrt(243) = rrt(243) * density(17) * density(30) 
  rrt(244) = rrt(244) * density(17) * density(30) 
  rrt(245) = rrt(245) * density(18) * density(30) 
  rrt(246) = rrt(246) * density(19) * density(30) 
  rrt(247) = rrt(247) * density(20) * density(30) 
  rrt(248) = rrt(248) * density(20) * density(30) 
  rrt(249) = rrt(249) * density(20) * density(30) 
  rrt(250) = rrt(250) * density(17) * density(38) 
  rrt(251) = rrt(251) * density(17) * density(38) 
  rrt(252) = rrt(252) * density(17) * density(23) 
  rrt(253) = rrt(253) * density(17) * density(23) 
  rrt(254) = rrt(254) * density(01) * density(17) 
  rrt(255) = rrt(255) * density(17) * density(50) 
  rrt(256) = rrt(256) * density(17) * density(51) 
  rrt(257) = rrt(257) * density(17) * density(52) 
  rrt(258) = rrt(258) * density(17)**2 
  rrt(259) = rrt(259) * density(17)**2 
  rrt(260) = rrt(260) * density(01) * density(18) 
  rrt(261) = rrt(261) * density(01) * density(18) 
  rrt(262) = rrt(262) * density(18) * density(50) 
  rrt(263) = rrt(263) * density(01) * density(20) 
  rrt(264) = rrt(264) * density(01) * density(19) 
  rrt(265) = rrt(265) * density(19) * density(50) 
  rrt(266) = rrt(266) * density(17) * density(19) 
  rrt(267) = rrt(267) * density(19)**2 
  rrt(268) = rrt(268) * density(17) * density(63) 
  rrt(269) = rrt(269) * density(17) * density(62) 
  rrt(270) = rrt(270) * density(17) * density(70) 
  rrt(271) = rrt(271) * density(18) * density(62) 
  rrt(272) = rrt(272) * density(19) * density(63) 
  rrt(273) = rrt(273) * density(19) * density(62) 
  rrt(274) = rrt(274) * density(01) * density(23)**2 
  rrt(275) = rrt(275) * density(23)**2 * density(62) 
  rrt(276) = rrt(276) * density(23)**2 * density(70) 
  rrt(277) = rrt(277) * density(23)**3 
  rrt(278) = rrt(278) * density(23)**2 * density(63) 
  rrt(279) = rrt(279) * density(01) * density(23)**2 
  rrt(280) = rrt(280) * density(23)**2 * density(62) 
  rrt(281) = rrt(281) * density(23)**2 * density(70) 
  rrt(282) = rrt(282) * density(23)**3 
  rrt(283) = rrt(283) * density(23)**2 * density(63) 
  rrt(284) = rrt(284) * density(24) * density(38) 
  rrt(285) = rrt(285) * density(24) * density(30) 
  rrt(286) = rrt(286) * density(24) * density(50) 
  rrt(287) = rrt(287) * density(24) * density(51) 
  rrt(288) = rrt(288) * density(01) * density(24) 
  rrt(289) = rrt(289) * density(23) * density(25) 
  rrt(290) = rrt(290) * density(25) * density(38) 
  rrt(291) = rrt(291) * density(23) * density(25) 
  rrt(292) = rrt(292) * density(01) * density(25) 
  rrt(293) = rrt(293) * density(24) * density(25) 
  rrt(294) = rrt(294) * density(25) * density(30) 
  rrt(295) = rrt(295) * density(25) * density(50) 
  rrt(296) = rrt(296) * density(24) * density(62) 
  rrt(297) = rrt(297) * density(24) * density(70) 
  rrt(298) = rrt(298) * density(25) * density(62) 
  rrt(299) = rrt(299) * density(35) * density(38) 
  rrt(300) = rrt(300) * density(23) * density(35) 
  rrt(301) = rrt(301) * density(30) * density(35) 
  rrt(302) = rrt(302) * density(01) * density(35) 
  rrt(303) = rrt(303) * density(35) * density(50) 
  rrt(304) = rrt(304) * density(35) * density(41) 
  rrt(305) = rrt(305) * density(35)**2 
  rrt(306) = rrt(306) * density(38) * density(41) 
  rrt(307) = rrt(307) * density(36) * density(38) 
  rrt(308) = rrt(308) * density(36) * density(38) 
  rrt(309) = rrt(309) * density(30) * density(36) 
  rrt(310) = rrt(310) * density(01) * density(36) 
  rrt(311) = rrt(311) * density(36) * density(50) 
  rrt(312) = rrt(312) * density(36) * density(41) 
  rrt(313) = rrt(313) * density(37) * density(38) 
  rrt(314) = rrt(314) * density(30) * density(37) 
  rrt(315) = rrt(315) * density(01) * density(37) 
  rrt(316) = rrt(316) * density(38) * density(39) 
  rrt(317) = rrt(317) * density(30) * density(39) 
  rrt(318) = rrt(318) * density(30) * density(39) 
  rrt(319) = rrt(319) * density(30) * density(39) 
  rrt(320) = rrt(320) * density(01) * density(39) 
  rrt(321) = rrt(321) * density(39) * density(41) 
  rrt(322) = rrt(322) * density(39) * density(41) 
  rrt(323) = rrt(323) * density(39) * density(50) 
  rrt(324) = rrt(324) * density(39) * density(51) 
  rrt(325) = rrt(325) * density(39) * density(51) 
  rrt(326) = rrt(326) * density(38) * density(40) 
  rrt(327) = rrt(327) * density(23) * density(40) 
  rrt(328) = rrt(328) * density(30) * density(40) 
  rrt(329) = rrt(329) * density(30) * density(40) 
  rrt(330) = rrt(330) * density(01) * density(40) 
  rrt(331) = rrt(331) * density(35) * density(40) 
  rrt(332) = rrt(332) * density(35) * density(40) 
  rrt(333) = rrt(333) * density(35) * density(40) 
  rrt(334) = rrt(334) * density(40) * density(50) 
  rrt(335) = rrt(335) * density(40) * density(50) 
  rrt(336) = rrt(336) * density(40) * density(41) 
  rrt(337) = rrt(337) * density(40) * density(41) 
  rrt(338) = rrt(338) * density(40) * density(51) 
  rrt(339) = rrt(339) * density(40) * density(51) 
  rrt(340) = rrt(340) * density(23) * density(50) 
  rrt(341) = rrt(341) * density(23) * density(30) 
  rrt(342) = rrt(342) * density(23) * density(52) 
  rrt(343) = rrt(343) * density(23) * density(52) 
  rrt(344) = rrt(344) * density(23) * density(52) 
  rrt(345) = rrt(345) * density(23) * density(52) 
  rrt(346) = rrt(346) * density(01) * density(38) 
  rrt(347) = rrt(347) * density(38) * density(50) 
  rrt(348) = rrt(348) * density(38) * density(50) 
  rrt(349) = rrt(349) * density(38) * density(51) 
  rrt(350) = rrt(350) * density(38) * density(51) 
  rrt(351) = rrt(351) * density(38) * density(52) 
  rrt(352) = rrt(352) * density(38) * density(53) 
  rrt(353) = rrt(353) * density(01) * density(30) 
  rrt(354) = rrt(354) * density(50)**2 
  rrt(355) = rrt(355) * density(50)**2 
  rrt(356) = rrt(356) * density(50)**2 
  rrt(357) = rrt(357) * density(30) * density(50) 
  rrt(358) = rrt(358) * density(41) * density(50) 
  rrt(359) = rrt(359) * density(50) * density(51) 
  rrt(360) = rrt(360) * density(50) * density(53) 
  rrt(361) = rrt(361) * density(30)**2 
  rrt(362) = rrt(362) * density(30) * density(52) 
  rrt(363) = rrt(363) * density(52)**2 
  rrt(364) = rrt(364) * density(52)**2 
  rrt(365) = rrt(365) * density(41) * density(52) 
  rrt(366) = rrt(366) * density(52) * density(53) 
  rrt(367) = rrt(367) * density(30) * density(53) 
  rrt(368) = rrt(368) * density(53)**2 
  rrt(369) = rrt(369) * density(23)**2 
  rrt(370) = rrt(370) * density(23) * density(38) 
  rrt(371) = rrt(371) * density(63) * density(72) 
  rrt(372) = rrt(372) * density(63) * density(71) 
  rrt(373) = rrt(373) * density(63) * density(70) 
  rrt(374) = rrt(374) * density(62) * density(71) 
  rrt(375) = rrt(375) * density(23) * density(72) 
  rrt(376) = rrt(376) * density(23) * density(71) 
  rrt(377) = rrt(377) * density(23) * density(71) 
  rrt(378) = rrt(378) * density(72)**2 
  rrt(379) = rrt(379) * density(72)**2 
  rrt(380) = rrt(380) * density(72)**2 
  rrt(381) = rrt(381) * density(71) * density(72) 
  rrt(382) = rrt(382) * density(63) * density(71) 
  rrt(383) = rrt(383) * density(71)**2 
  rrt(384) = rrt(384) * density(01)**2 
  rrt(385) = rrt(385) * density(01) * density(30) 
  rrt(386) = rrt(386) * density(01) * density(50) 
  rrt(387) = rrt(387) * density(01) * density(38) 
  rrt(388) = rrt(388) * density(01) * density(23) 
  rrt(389) = rrt(389) * density(01) * density(30) 
  rrt(390) = rrt(390) * density(30)**2 
  rrt(391) = rrt(391) * density(30) * density(38) 
  rrt(392) = rrt(392) * density(23) * density(30) 
  rrt(393) = rrt(393) * density(30) * density(50) 
  rrt(394) = rrt(394) * density(01) * density(50) 
  rrt(395) = rrt(395) * density(30) * density(50) 
  rrt(396) = rrt(396) * density(38) * density(50) 
  rrt(397) = rrt(397) * density(23) * density(50) 
  rrt(398) = rrt(398) * density(50)**2 
  rrt(399) = rrt(399) * density(01) * density(41) 
  rrt(400) = rrt(400) * density(30) * density(41) 
  rrt(401) = rrt(401) * density(23) * density(41) 
  rrt(402) = rrt(402) * density(38) * density(41) 
  rrt(403) = rrt(403) * density(01) * density(51) 
  rrt(404) = rrt(404) * density(30) * density(51) 
  rrt(405) = rrt(405) * density(50) * density(51) 
  rrt(406) = rrt(406) * density(51)**2 
  rrt(407) = rrt(407) * density(01) * density(52) 
  rrt(408) = rrt(408) * density(30) * density(52) 
  rrt(409) = rrt(409) * density(50) * density(52) 
  rrt(410) = rrt(410) * density(52)**2 
  rrt(411) = rrt(411) * density(01) * density(53) 
  rrt(412) = rrt(412) * density(30) * density(53) 
  rrt(413) = rrt(413) * density(50) * density(53) 
  rrt(414) = rrt(414) * density(23) * density(53) 
  rrt(415) = rrt(415) * density(38) * density(53) 
  rrt(416) = rrt(416) * density(01) * density(53) 
  rrt(417) = rrt(417) * density(30) * density(53) 
  rrt(418) = rrt(418) * density(50) * density(53) 
  rrt(419) = rrt(419) * density(23) * density(53) 
  rrt(420) = rrt(420) * density(38) * density(53) 
  rrt(421) = rrt(421) * density(54) 
  rrt(422) = rrt(422) * density(01) * density(23)**2 
  rrt(423) = rrt(423) * density(23)**2 * density(30) 
  rrt(424) = rrt(424) * density(23)**2 * density(50) 
  rrt(425) = rrt(425) * density(23)**3 
  rrt(426) = rrt(426) * density(23)**2 * density(38) 
  rrt(427) = rrt(427) * density(01) * density(38)**2 
  rrt(428) = rrt(428) * density(30) * density(38)**2 
  rrt(429) = rrt(429) * density(23) * density(38)**2 
  rrt(430) = rrt(430) * density(38)**3 
  rrt(431) = rrt(431) * density(38)**2 * density(50) 
  rrt(432) = rrt(432) * density(01) * density(23) * density(38) 
  rrt(433) = rrt(433) * density(23) * density(30) * density(38) 
  rrt(434) = rrt(434) * density(23)**2 * density(38) 
  rrt(435) = rrt(435) * density(23) * density(38)**2 
  rrt(436) = rrt(436) * density(23) * density(38) * density(50) 
  rrt(437) = rrt(437) * density(01) * density(30) * density(38) 
  rrt(438) = rrt(438) * density(30)**2 * density(38) 
  rrt(439) = rrt(439) * density(30) * density(38) * density(50) 
  rrt(440) = rrt(440) * density(23) * density(30) * density(38) 
  rrt(441) = rrt(441) * density(30) * density(38)**2 
  rrt(442) = rrt(442) * density(01) * density(38) 
  rrt(443) = rrt(443) * density(01) * density(38) * density(50) 
  rrt(444) = rrt(444) * density(30) * density(38) * density(50) 
  rrt(445) = rrt(445) * density(38) * density(50)**2 
  rrt(446) = rrt(446) * density(01) * density(38) * density(52) 
  rrt(447) = rrt(447) * density(30) * density(38) * density(52) 
  rrt(448) = rrt(448) * density(23) * density(38) * density(52) 
  rrt(449) = rrt(449) * density(38)**2 * density(52) 
  rrt(450) = rrt(450) * density(38) * density(50) * density(52) 
  rrt(451) = rrt(451) * density(52) * density(53) 
  rrt(452) = rrt(452) * density(01) * density(23)**2 
  rrt(453) = rrt(453) * density(23)**2 * density(62) 
  rrt(454) = rrt(454) * density(23)**2 * density(70) 
  rrt(455) = rrt(455) * density(62) * density(63)**2 
  rrt(456) = rrt(456) * density(01) * density(63)**2 
  rrt(457) = rrt(457) * density(01) * density(23) * density(63) 
  rrt(458) = rrt(458) * density(23) * density(62) * density(63) 
  rrt(459) = rrt(459) * density(23) * density(63) * density(70) 
  rrt(460) = rrt(460) * density(01) * density(23) * density(62) 
  rrt(461) = rrt(461) * density(23) * density(62)**2 
  rrt(462) = rrt(462) * density(23) * density(62) * density(70) 
  rrt(463) = rrt(463) * density(01) * density(63) * density(72) 
  rrt(464) = rrt(464) * density(62) * density(63) * density(72) 
  rrt(465) = rrt(465) * density(63) * density(70) * density(72) 
  rrt(466) = rrt(466) * density(01) * density(63) * density(71) 
  rrt(467) = rrt(467) * density(62) * density(63) * density(71) 
  rrt(468) = rrt(468) * density(63) * density(70) * density(71) 
  rrt(469) = rrt(469) * density(01) * density(62) * density(72) 
  rrt(470) = rrt(470) * density(62)**2 * density(72) 
  rrt(471) = rrt(471) * density(62) * density(70) * density(72) 
  rrt(472) = rrt(472) * density(26) * density(38) 
  rrt(473) = rrt(473) * density(26) * density(30) 
  rrt(474) = rrt(474) * density(26) * density(30) 
  rrt(475) = rrt(475) * density(26) * density(30) 
  rrt(476) = rrt(476) * density(26) * density(41) 
  rrt(477) = rrt(477) * density(26) * density(50) 
  rrt(478) = rrt(478) * density(26) * density(50) 
  rrt(479) = rrt(479) * density(26) * density(50) 
  rrt(480) = rrt(480) * density(26) * density(51) 
  rrt(481) = rrt(481) * density(01) * density(42) 
  rrt(482) = rrt(482) * density(30) * density(42) 
  rrt(483) = rrt(483) * density(41) * density(42) 
  rrt(484) = rrt(484) * density(42) * density(50) 
  rrt(485) = rrt(485) * density(42) * density(50) 
  rrt(486) = rrt(486) * density(24) * density(42) 
  rrt(487) = rrt(487) * density(42) * density(51) 
  rrt(488) = rrt(488) * density(42) * density(51) 
  rrt(489) = rrt(489) * density(42) * density(51) 
  rrt(490) = rrt(490) * density(42) * density(52) 
  rrt(491) = rrt(491) * density(27) * density(30) 
  rrt(492) = rrt(492) * density(27) * density(38) 
  rrt(493) = rrt(493) * density(27) * density(41) 
  rrt(494) = rrt(494) * density(23) * density(27) 
  rrt(495) = rrt(495) * density(27) * density(50) 
  rrt(496) = rrt(496) * density(27) * density(51) 
  rrt(497) = rrt(497) * density(27) * density(51) 
  rrt(498) = rrt(498) * density(01) * density(43) 
  rrt(499) = rrt(499) * density(23) * density(43) 
  rrt(500) = rrt(500) * density(43) * density(50) 
  rrt(501) = rrt(501) * density(43) * density(52) 
  rrt(502) = rrt(502) * density(43) * density(52) 
  rrt(503) = rrt(503) * density(28) * density(30) 
  rrt(504) = rrt(504) * density(28) * density(30) 
  rrt(505) = rrt(505) * density(23) * density(28) 
  rrt(506) = rrt(506) * density(28) * density(50) 
  rrt(507) = rrt(507) * density(28) * density(50) 
  rrt(508) = rrt(508) * density(50) * density(57) 
  rrt(509) = rrt(509) * density(50) * density(56) 
  rrt(510) = rrt(510) * density(01) * density(29) 
  rrt(511) = rrt(511) * density(29) * density(30) 
  rrt(512) = rrt(512) * density(29) * density(38) 
  rrt(513) = rrt(513) * density(23) * density(29) 
  rrt(514) = rrt(514) * density(29) * density(50) 
  rrt(515) = rrt(515) * density(01) * density(45) 
  rrt(516) = rrt(516) * density(30) * density(45) 
  rrt(517) = rrt(517) * density(35) * density(45) 
  rrt(518) = rrt(518) * density(36) * density(45) 
  rrt(519) = rrt(519) * density(38) * density(45) 
  rrt(520) = rrt(520) * density(45) * density(50) 
  rrt(521) = rrt(521) * density(01) * density(78) 
  rrt(522) = rrt(522) * density(30) * density(78) 
  rrt(523) = rrt(523) * density(01)**2 * density(26) 
  rrt(524) = rrt(524) * density(26) * density(38) 
  rrt(525) = rrt(525) * density(23) * density(26) 
  rrt(526) = rrt(526) * density(01) * density(42) 
  rrt(527) = rrt(527) * density(38) * density(42) 
  rrt(528) = rrt(528) * density(23) * density(42) 
  rrt(529) = rrt(529) * density(01)**2 * density(27) 
  rrt(530) = rrt(530) * density(01) * density(23) * density(27) 
  rrt(531) = rrt(531) * density(30)**2 * density(43) 
  rrt(532) = rrt(532) * density(01)**2 * density(43) 
  rrt(533) = rrt(533) * density(26) * density(62) 
  rrt(534) = rrt(534) * density(26) * density(70) 
  rrt(535) = rrt(535) * density(26) * density(70) 
  rrt(536) = rrt(536) * density(64) * density(70) 
  rrt(537) = rrt(537) * density(63) * density(65) 
  rrt(538) = rrt(538) * density(62) * density(65) 
  rrt(539) = rrt(539) * density(65) * density(70) 
  rrt(540) = rrt(540) * density(01) * density(65) 
  rrt(541) = rrt(541) * density(62) * density(75) 
  rrt(542) = rrt(542) * density(62) * density(75) 
  rrt(543) = rrt(543) * density(70) * density(75) 
  rrt(544) = rrt(544) * density(70) * density(75) 
  rrt(545) = rrt(545) * density(01) * density(75) 
  rrt(546) = rrt(546) * density(62) * density(74) 
  rrt(547) = rrt(547) * density(70) * density(74) 
  rrt(548) = rrt(548) * density(70) * density(74) 
  rrt(549) = rrt(549) * density(70) * density(76) 
  rrt(550) = rrt(550) * density(27) * density(70) 
  rrt(551) = rrt(551) * density(26) * density(70) 
  rrt(552) = rrt(552) * density(43) * density(70) 
  rrt(553) = rrt(553) * density(69) * density(70) 
  rrt(554) = rrt(554) * density(70) * density(73) 
  rrt(555) = rrt(555) * density(71) * density(73) 
  rrt(556) = rrt(556) * density(35) * density(46) 
  rrt(557) = rrt(557) * density(41) * density(46) 
  rrt(558) = rrt(558) * density(46) * density(52) 
  rrt(559) = rrt(559) * density(46) * density(51) 
  rrt(560) = rrt(560) * density(46) * density(51) 
  rrt(561) = rrt(561) * density(38) * density(47) 
  rrt(562) = rrt(562) * density(41) * density(47) 
  rrt(563) = rrt(563) * density(47) * density(52) 
  rrt(564) = rrt(564) * density(47) * density(53) 
  rrt(565) = rrt(565) * density(38) * density(48) 
  rrt(566) = rrt(566) * density(48) * density(50) 
  rrt(567) = rrt(567) * density(48) * density(50) 
  rrt(568) = rrt(568) * density(48) * density(52) 
  rrt(569) = rrt(569) * density(48) * density(52) 
  rrt(570) = rrt(570) * density(48) * density(53) 
  rrt(571) = rrt(571) * density(30) * density(58) 
  rrt(572) = rrt(572) * density(52) * density(58) 
  rrt(573) = rrt(573) * density(51) * density(58) 
  rrt(574) = rrt(574) * density(41) * density(60) 
  rrt(575) = rrt(575) * density(52) * density(60) 
  rrt(576) = rrt(576) * density(53) * density(60) 
  rrt(577) = rrt(577) * density(54) * density(60) 
  rrt(578) = rrt(578) * density(50) * density(61) 
  rrt(579) = rrt(579) * density(01) * density(49) 
  rrt(580) = rrt(580) * density(30) * density(49) 
  rrt(581) = rrt(581) * density(38) * density(49) 
  rrt(582) = rrt(582) * density(38) * density(49) 
  rrt(583) = rrt(583) * density(35) * density(49) 
  rrt(584) = rrt(584) * density(36) * density(49) 
  rrt(585) = rrt(585) * density(49) * density(50) 
  rrt(586) = rrt(586) * density(30) * density(46) 
  rrt(587) = rrt(587) * density(46) * density(50) 
  rrt(588) = rrt(588) * density(30) * density(47) 
  rrt(589) = rrt(589) * density(26) * density(46) 
  rrt(590) = rrt(590) * density(27) * density(46) 
  rrt(591) = rrt(591) * density(42) * density(46) 
  rrt(592) = rrt(592) * density(43) * density(46) 
  rrt(593) = rrt(593) * density(46) * density(55) 
  rrt(594) = rrt(594) * density(46) * density(56) 
  rrt(595) = rrt(595) * density(46) * density(57) 
  rrt(596) = rrt(596) * density(26) * density(47) 
  rrt(597) = rrt(597) * density(27) * density(47) 
  rrt(598) = rrt(598) * density(42) * density(47) 
  rrt(599) = rrt(599) * density(43) * density(47) 
  rrt(600) = rrt(600) * density(47) * density(55) 
  rrt(601) = rrt(601) * density(47) * density(56) 
  rrt(602) = rrt(602) * density(47) * density(57) 
  rrt(603) = rrt(603) * density(26) * density(48) 
  rrt(604) = rrt(604) * density(27) * density(48) 
  rrt(605) = rrt(605) * density(42) * density(48) 
  rrt(606) = rrt(606) * density(43) * density(48) 
  rrt(607) = rrt(607) * density(48) * density(55) 
  rrt(608) = rrt(608) * density(48) * density(56) 
  rrt(609) = rrt(609) * density(48) * density(57) 
  rrt(610) = rrt(610) * density(26) * density(58) 
  rrt(611) = rrt(611) * density(27) * density(58) 
  rrt(612) = rrt(612) * density(42) * density(58) 
  rrt(613) = rrt(613) * density(43) * density(58) 
  rrt(614) = rrt(614) * density(55) * density(58) 
  rrt(615) = rrt(615) * density(56) * density(58) 
  rrt(616) = rrt(616) * density(57) * density(58) 
  rrt(617) = rrt(617) * density(26) * density(59) 
  rrt(618) = rrt(618) * density(27) * density(59) 
  rrt(619) = rrt(619) * density(42) * density(59) 
  rrt(620) = rrt(620) * density(43) * density(59) 
  rrt(621) = rrt(621) * density(55) * density(59) 
  rrt(622) = rrt(622) * density(56) * density(59) 
  rrt(623) = rrt(623) * density(57) * density(59) 
  rrt(624) = rrt(624) * density(26) * density(60) 
  rrt(625) = rrt(625) * density(27) * density(60) 
  rrt(626) = rrt(626) * density(42) * density(60) 
  rrt(627) = rrt(627) * density(43) * density(60) 
  rrt(628) = rrt(628) * density(55) * density(60) 
  rrt(629) = rrt(629) * density(56) * density(60) 
  rrt(630) = rrt(630) * density(57) * density(60) 
  rrt(631) = rrt(631) * density(26) * density(61) 
  rrt(632) = rrt(632) * density(27) * density(61) 
  rrt(633) = rrt(633) * density(42) * density(61) 
  rrt(634) = rrt(634) * density(43) * density(61) 
  rrt(635) = rrt(635) * density(55) * density(61) 
  rrt(636) = rrt(636) * density(56) * density(61) 
  rrt(637) = rrt(637) * density(57) * density(61) 
  rrt(638) = rrt(638) * density(27) * density(46) 
  rrt(639) = rrt(639) * density(28) * density(46) 
  rrt(640) = rrt(640) * density(29) * density(46) 
  rrt(641) = rrt(641) * density(43) * density(46) 
  rrt(642) = rrt(642) * density(45) * density(46) 
  rrt(643) = rrt(643) * density(46) * density(55) 
  rrt(644) = rrt(644) * density(46) * density(56) 
  rrt(645) = rrt(645) * density(46) * density(57) 
  rrt(646) = rrt(646) * density(46) * density(78) 
  rrt(647) = rrt(647) * density(27) * density(47) 
  rrt(648) = rrt(648) * density(28) * density(47) 
  rrt(649) = rrt(649) * density(29) * density(47) 
  rrt(650) = rrt(650) * density(43) * density(47) 
  rrt(651) = rrt(651) * density(45) * density(47) 
  rrt(652) = rrt(652) * density(47) * density(55) 
  rrt(653) = rrt(653) * density(47) * density(56) 
  rrt(654) = rrt(654) * density(47) * density(57) 
  rrt(655) = rrt(655) * density(47) * density(78) 
  rrt(656) = rrt(656) * density(27) * density(48) 
  rrt(657) = rrt(657) * density(28) * density(48) 
  rrt(658) = rrt(658) * density(29) * density(48) 
  rrt(659) = rrt(659) * density(43) * density(48) 
  rrt(660) = rrt(660) * density(45) * density(48) 
  rrt(661) = rrt(661) * density(48) * density(55) 
  rrt(662) = rrt(662) * density(48) * density(56) 
  rrt(663) = rrt(663) * density(48) * density(57) 
  rrt(664) = rrt(664) * density(48) * density(78) 
  rrt(665) = rrt(665) * density(27) * density(58) 
  rrt(666) = rrt(666) * density(28) * density(58) 
  rrt(667) = rrt(667) * density(29) * density(58) 
  rrt(668) = rrt(668) * density(43) * density(58) 
  rrt(669) = rrt(669) * density(45) * density(58) 
  rrt(670) = rrt(670) * density(55) * density(58) 
  rrt(671) = rrt(671) * density(56) * density(58) 
  rrt(672) = rrt(672) * density(57) * density(58) 
  rrt(673) = rrt(673) * density(58) * density(78) 
  rrt(674) = rrt(674) * density(27) * density(59) 
  rrt(675) = rrt(675) * density(28) * density(59) 
  rrt(676) = rrt(676) * density(29) * density(59) 
  rrt(677) = rrt(677) * density(43) * density(59) 
  rrt(678) = rrt(678) * density(45) * density(59) 
  rrt(679) = rrt(679) * density(55) * density(59) 
  rrt(680) = rrt(680) * density(56) * density(59) 
  rrt(681) = rrt(681) * density(57) * density(59) 
  rrt(682) = rrt(682) * density(59) * density(78) 
  rrt(683) = rrt(683) * density(27) * density(60) 
  rrt(684) = rrt(684) * density(28) * density(60) 
  rrt(685) = rrt(685) * density(29) * density(60) 
  rrt(686) = rrt(686) * density(43) * density(60) 
  rrt(687) = rrt(687) * density(45) * density(60) 
  rrt(688) = rrt(688) * density(55) * density(60) 
  rrt(689) = rrt(689) * density(56) * density(60) 
  rrt(690) = rrt(690) * density(57) * density(60) 
  rrt(691) = rrt(691) * density(60) * density(78) 
  rrt(692) = rrt(692) * density(27) * density(61) 
  rrt(693) = rrt(693) * density(28) * density(61) 
  rrt(694) = rrt(694) * density(29) * density(61) 
  rrt(695) = rrt(695) * density(43) * density(61) 
  rrt(696) = rrt(696) * density(45) * density(61) 
  rrt(697) = rrt(697) * density(55) * density(61) 
  rrt(698) = rrt(698) * density(56) * density(61) 
  rrt(699) = rrt(699) * density(57) * density(61) 
  rrt(700) = rrt(700) * density(61) * density(78) 
  rrt(701) = rrt(701) * density(26) * density(49) 
  rrt(702) = rrt(702) * density(27) * density(49) 
  rrt(703) = rrt(703) * density(42) * density(49) 
  rrt(704) = rrt(704) * density(43) * density(49) 
  rrt(705) = rrt(705) * density(49) * density(55) 
  rrt(706) = rrt(706) * density(49) * density(56) 
  rrt(707) = rrt(707) * density(49) * density(57) 
  rrt(708) = rrt(708) * density(28) * density(49) 
  rrt(709) = rrt(709) * density(29) * density(49) 
  rrt(710) = rrt(710) * density(45) * density(49) 
  rrt(711) = rrt(711) * density(49) * density(78) 
  rrt(712) = rrt(712) * density(26) * density(46) 
  rrt(713) = rrt(713) * density(27) * density(46) 
  rrt(714) = rrt(714) * density(42) * density(46) 
  rrt(715) = rrt(715) * density(43) * density(46) 
  rrt(716) = rrt(716) * density(46) * density(55) 
  rrt(717) = rrt(717) * density(26) * density(47) 
  rrt(718) = rrt(718) * density(27) * density(47) 
  rrt(719) = rrt(719) * density(42) * density(47) 
  rrt(720) = rrt(720) * density(43) * density(47) 
  rrt(721) = rrt(721) * density(47) * density(55) 
  rrt(722) = rrt(722) * density(26) * density(46) 
  rrt(723) = rrt(723) * density(27) * density(46) 
  rrt(724) = rrt(724) * density(42) * density(46) 
  rrt(725) = rrt(725) * density(43) * density(46) 
  rrt(726) = rrt(726) * density(46) * density(55) 
  rrt(727) = rrt(727) * density(26) * density(47) 
  rrt(728) = rrt(728) * density(42) * density(47) 
  rrt(729) = rrt(729) * density(47) * density(55) 
  rrt(730) = rrt(730) * density(47) * density(73) 
  rrt(731) = rrt(731) * density(46) * density(73) 
  rrt(732) = rrt(732) * density(48) * density(73) 
  rrt(733) = rrt(733) * density(58) * density(73) 
  rrt(734) = rrt(734) * density(60) * density(73) 
  rrt(735) = rrt(735) * density(61) * density(73) 
  rrt(736) = rrt(736) * density(59) * density(73) 
  rrt(737) = rrt(737) * density(61) * density(77) 
  rrt(738) = rrt(738) * density(65) * density(67) 
  rrt(739) = rrt(739) * density(66) * density(67) 
  rrt(740) = rrt(740) * density(27) * density(67) 
  rrt(741) = rrt(741) * density(29) * density(67) 
  rrt(742) = rrt(742) * density(67) * density(76) 
  rrt(743) = rrt(743) * density(01) * density(65) * density(67) 
  rrt(744) = rrt(744) * density(62) * density(65) * density(67) 
  rrt(745) = rrt(745) * density(65) * density(67) * density(70) 
  rrt(746) = rrt(746) * density(01) * density(66) * density(67) 
  rrt(747) = rrt(747) * density(62) * density(66) * density(67) 
  rrt(748) = rrt(748) * density(66) * density(67) * density(70) 
  rrt(749) = rrt(749) * density(01) * density(27) * density(67) 
  rrt(750) = rrt(750) * density(27) * density(62) * density(67) 
  rrt(751) = rrt(751) * density(27) * density(67) * density(70) 
  rrt(752) = rrt(752) * density(01) * density(29) * density(67) 
  rrt(753) = rrt(753) * density(29) * density(62) * density(67) 
  rrt(754) = rrt(754) * density(29) * density(67) * density(70) 
  rrt(755) = rrt(755) * density(01) * density(67) * density(76) 
  rrt(756) = rrt(756) * density(62) * density(67) * density(76) 
  rrt(757) = rrt(757) * density(67) * density(70) * density(76) 
  rrt(758) = rrt(758) * density(26) * density(48) 
  rrt(759) = rrt(759) * density(27) * density(48) 
  rrt(760) = rrt(760) * density(42) * density(48) 
  rrt(761) = rrt(761) * density(43) * density(48) 
  rrt(762) = rrt(762) * density(48) * density(55) 
  rrt(763) = rrt(763) * density(48) * density(56) 
  rrt(764) = rrt(764) * density(48) * density(57) 
  rrt(765) = rrt(765) * density(26) * density(58) 
  rrt(766) = rrt(766) * density(27) * density(58) 
  rrt(767) = rrt(767) * density(42) * density(58) 
  rrt(768) = rrt(768) * density(43) * density(58) 
  rrt(769) = rrt(769) * density(55) * density(58) 
  rrt(770) = rrt(770) * density(56) * density(58) 
  rrt(771) = rrt(771) * density(57) * density(58) 
  rrt(772) = rrt(772) * density(26) * density(59) 
  rrt(773) = rrt(773) * density(27) * density(59) 
  rrt(774) = rrt(774) * density(42) * density(59) 
  rrt(775) = rrt(775) * density(43) * density(59) 
  rrt(776) = rrt(776) * density(55) * density(59) 
  rrt(777) = rrt(777) * density(56) * density(59) 
  rrt(778) = rrt(778) * density(57) * density(59) 
  rrt(779) = rrt(779) * density(26) * density(60) 
  rrt(780) = rrt(780) * density(27) * density(60) 
  rrt(781) = rrt(781) * density(42) * density(60) 
  rrt(782) = rrt(782) * density(43) * density(60) 
  rrt(783) = rrt(783) * density(55) * density(60) 
  rrt(784) = rrt(784) * density(56) * density(60) 
  rrt(785) = rrt(785) * density(57) * density(60) 
  rrt(786) = rrt(786) * density(26) * density(61) 
  rrt(787) = rrt(787) * density(27) * density(61) 
  rrt(788) = rrt(788) * density(42) * density(61) 
  rrt(789) = rrt(789) * density(43) * density(61) 
  rrt(790) = rrt(790) * density(55) * density(61) 
  rrt(791) = rrt(791) * density(56) * density(61) 
  rrt(792) = rrt(792) * density(57) * density(61) 
  ydot(01) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)-rrt(010)-rrt(011)-rrt(012)-rrt(013)&
             -rrt(014)-rrt(015)+rrt(016)+rrt(020)+rrt(025)+rrt(031)+rrt(035)+rrt(038)+rrt(040)+rrt(041)+rrt(042)+rrt(043)+rrt(050)&
             -rrt(058)+rrt(066)-rrt(074)+rrt(082)-rrt(090)-rrt(114)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)-rrt(120)-rrt(121)&
             -rrt(122)-rrt(123)-rrt(124)-rrt(125)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)-rrt(134)&
             -rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)-rrt(140)-rrt(141)-rrt(142)-rrt(156)+rrt(174)+  2.d0 * rrt(175)+rrt(176)&
             +rrt(179)-rrt(200)+rrt(204)+rrt(205)+rrt(213)+rrt(214)+rrt(218)+rrt(225)+rrt(226)+rrt(227)+rrt(228)+rrt(229)+rrt(230)&
             +rrt(231)+rrt(232)+rrt(233)+rrt(234)+rrt(235)+rrt(237)+rrt(243)+rrt(244)+rrt(245)+rrt(246)+rrt(247)+rrt(248)+rrt(249)&
             +rrt(251)+rrt(252)+rrt(253)+rrt(254)+rrt(255)+rrt(256)+rrt(257)+rrt(258)+rrt(259)+rrt(261)+rrt(265)+rrt(268)+rrt(269)&
             +rrt(270)+rrt(272)+rrt(273)+rrt(286)+rrt(287)+rrt(325)+rrt(340)+rrt(342)+rrt(344)-rrt(346)+rrt(349)-rrt(353)+rrt(356)&
             +rrt(359)+rrt(375)+rrt(376)+rrt(377)+rrt(379)+rrt(380)-rrt(384)-rrt(385)-rrt(386)-rrt(387)-rrt(388)+rrt(403)+rrt(404)&
             +rrt(405)+rrt(406)+rrt(422)+rrt(423)+rrt(424)+rrt(425)+rrt(426)-rrt(442)+rrt(452)+rrt(453)+rrt(454)+rrt(479)+rrt(480)&
             -rrt(481)+rrt(489)+rrt(491)+rrt(493)+rrt(494)+rrt(495)+rrt(496)+rrt(497)-rrt(498)+rrt(503)+rrt(504)+rrt(505)+rrt(506)&
             +rrt(507)+rrt(510)+  2.d0 * rrt(511)+  2.d0 * rrt(512)+  2.d0 * rrt(513)+  2.d0 * rrt(514)-rrt(515)+rrt(521)+rrt(522)&
             -rrt(523)-rrt(526)-rrt(529)-rrt(532)-rrt(540)-rrt(545)+rrt(549)+rrt(550)+rrt(573)+rrt(590)+rrt(597)+rrt(604)+rrt(611)&
             +rrt(618)+rrt(625)+rrt(632)+rrt(639)+  2.d0 * rrt(640)+rrt(644)+rrt(646)+rrt(648)+  2.d0 * rrt(649)+rrt(653)+rrt(655)&
             +rrt(657)+  2.d0 * rrt(658)+rrt(662)+rrt(664)+rrt(666)+  2.d0 * rrt(667)+rrt(671)+rrt(673)+rrt(675)+  2.d0 * rrt(676)&
             +rrt(680)+rrt(682)+rrt(684)+  2.d0 * rrt(685)+rrt(689)+rrt(691)+rrt(693)+  2.d0 * rrt(694)+rrt(698)+rrt(700)+rrt(702)&
             +rrt(708)+  2.d0 * rrt(709)+rrt(711)+rrt(713)+rrt(718)+rrt(740)+  2.d0 * rrt(741)+rrt(742)+rrt(749)+rrt(750)+rrt(751)&
             +  2.d0 * rrt(752)+  2.d0 * rrt(753)+  2.d0 * rrt(754)+rrt(755)+rrt(756)+rrt(757)+rrt(759)+rrt(766)+rrt(773)+rrt(780)
  ydot(01) = ydot(01) &
             +rrt(787) 
  ydot(02) = +rrt(001)-rrt(016)-rrt(017)-rrt(018)-rrt(019)+rrt(021)+rrt(026)+rrt(032)-rrt(050)+rrt(051)+rrt(058)-rrt(059)-rrt(066)&
             +rrt(067)+rrt(074)-rrt(075)-rrt(082)+rrt(083)+rrt(090)-rrt(091) 
  ydot(03) = +rrt(002)+rrt(017)-rrt(020)-rrt(021)-rrt(022)-rrt(023)-rrt(024)+rrt(027)+rrt(033)+rrt(036)-rrt(051)+rrt(052)+rrt(059)&
             -rrt(060)-rrt(067)+rrt(068)+rrt(075)-rrt(076)-rrt(083)+rrt(084)+rrt(091)-rrt(092) 
  ydot(04) = +rrt(003)+rrt(018)+rrt(022)-rrt(025)-rrt(026)-rrt(027)-rrt(028)-rrt(029)-rrt(030)+rrt(034)+rrt(037)+rrt(039)-rrt(052)&
             +rrt(053)+rrt(060)-rrt(061)-rrt(068)+rrt(069)+rrt(076)-rrt(077)-rrt(084)+rrt(085)+rrt(092)-rrt(093) 
  ydot(05) = +rrt(004)+rrt(019)+rrt(023)+rrt(028)-rrt(031)-rrt(032)-rrt(033)-rrt(034)-rrt(053)+rrt(054)+rrt(061)-rrt(062)-rrt(069)&
             +rrt(070)+rrt(077)-rrt(078)-rrt(085)+rrt(086)+rrt(093)-rrt(094) 
  ydot(06) = +rrt(005)+rrt(024)+rrt(029)-rrt(035)-rrt(036)-rrt(037)-rrt(054)+rrt(055)+rrt(062)-rrt(063)-rrt(070)+rrt(071)+rrt(078)&
             -rrt(079)-rrt(086)+rrt(087)+rrt(094)-rrt(095) 
  ydot(07) = +rrt(006)+rrt(030)-rrt(038)-rrt(039)-rrt(055)+rrt(056)+rrt(063)-rrt(064)-rrt(071)+rrt(072)+rrt(079)-rrt(080)-rrt(087)&
             +rrt(088)+rrt(095)-rrt(096) 
  ydot(08) = +rrt(007)-rrt(040)-rrt(056)+rrt(057)+rrt(064)-rrt(065)-rrt(072)+rrt(073)+rrt(080)-rrt(081)-rrt(088)+rrt(089)+rrt(096)&
             -rrt(097) 
  ydot(09) = +rrt(008)-rrt(041)-rrt(057)+rrt(065)-rrt(073)+rrt(081)-rrt(089)+rrt(097) 
  ydot(10) = +rrt(009)-rrt(042) 
  ydot(11) = +rrt(010)-rrt(043) 
  ydot(12) = +rrt(011) 
  ydot(13) = +rrt(012) 
  ydot(14) = +rrt(013) 
  ydot(15) = +rrt(014) 
  ydot(16) = +rrt(015) 
  ydot(17) = +rrt(114)+rrt(115)+rrt(116)-rrt(157)-rrt(204)-rrt(213)-rrt(225)-rrt(226)-rrt(227)-rrt(228)-rrt(229)-rrt(235)+rrt(236)&
             -rrt(243)-rrt(244)-rrt(250)-rrt(251)-rrt(252)-rrt(253)-rrt(254)-rrt(255)-rrt(256)-rrt(257)-  2.d0 * rrt(258)&
             -  2.d0 * rrt(259)+rrt(260)+rrt(262)-rrt(266)-rrt(268)-rrt(269)-rrt(270)+rrt(271)+rrt(274)+rrt(275)+rrt(276)+rrt(277)&
             +rrt(278)+rrt(295) 
  ydot(18) = +rrt(117)+rrt(118)+rrt(119)+rrt(120)+rrt(121)+rrt(122)+rrt(123)-rrt(205)-rrt(214)-rrt(230)-rrt(231)-rrt(232)-rrt(233)&
             -rrt(234)-rrt(236)+rrt(238)-rrt(245)+rrt(258)-rrt(260)-rrt(261)-rrt(262)+rrt(264)-rrt(271)+rrt(279)+rrt(280)+rrt(281)&
             +rrt(282)+rrt(283) 
  ydot(19) = +rrt(124)+rrt(125)+rrt(126)+rrt(127)+rrt(128)+rrt(129)-rrt(237)-rrt(246)+rrt(263)-rrt(264)-rrt(265)-rrt(266)&
             -  2.d0 * rrt(267)-rrt(272)-rrt(273) 
  ydot(20) = +rrt(130)+rrt(131)+rrt(132)-rrt(238)-rrt(247)-rrt(248)-rrt(249)+rrt(259)-rrt(263) 
  ydot(21) = +rrt(133)+rrt(134)+rrt(135)+rrt(136)+rrt(137)+rrt(138) 
  ydot(22) = +rrt(139)+rrt(140)+rrt(141) 
  ydot(23) = +rrt(142)-rrt(154)+  2.d0 * rrt(166)+rrt(167)+rrt(168)+rrt(172)+rrt(174)+rrt(180)+rrt(182)+rrt(186)-rrt(198)-rrt(208)&
             -rrt(216)-rrt(217)-rrt(218)-rrt(219)-rrt(220)-rrt(253)+rrt(256)+rrt(265)-  2.d0 * rrt(274)-  2.d0 * rrt(275)&
             -  2.d0 * rrt(276)-  2.d0 * rrt(277)-  2.d0 * rrt(278)-  2.d0 * rrt(279)-  2.d0 * rrt(280)-  2.d0 * rrt(281)&
             -  2.d0 * rrt(282)-  2.d0 * rrt(283)+rrt(284)+rrt(288)+rrt(289)+rrt(290)+rrt(292)-rrt(300)+rrt(323)-rrt(340)-rrt(341)&
             -rrt(342)-rrt(343)-rrt(344)-rrt(345)+rrt(346)+rrt(347)+rrt(354)-  2.d0 * rrt(369)-rrt(370)+rrt(371)-rrt(375)-rrt(376)&
             -rrt(377)+rrt(378)+rrt(381)+  2.d0 * rrt(384)+  2.d0 * rrt(385)+  2.d0 * rrt(386)+  2.d0 * rrt(387)+  2.d0 * rrt(388)&
             +rrt(394)+rrt(395)+rrt(396)+rrt(397)+rrt(398)-  2.d0 * rrt(422)-  2.d0 * rrt(423)-  2.d0 * rrt(424)-  2.d0 * rrt(425)&
             -  2.d0 * rrt(426)-rrt(432)-rrt(433)-rrt(434)-rrt(435)-rrt(436)-  2.d0 * rrt(452)-  2.d0 * rrt(453)-  2.d0 * rrt(454)&
             -rrt(457)-rrt(458)-rrt(459)-rrt(460)-rrt(461)-rrt(462)+rrt(472)+rrt(473)+rrt(477)+rrt(481)+rrt(485)+rrt(492)-rrt(494)&
             +rrt(497)-rrt(499)+rrt(503)-rrt(505)+rrt(506)-rrt(513)-rrt(525)+rrt(526)-rrt(528)-rrt(530)+rrt(541)+rrt(544)+rrt(545)&
             +rrt(551)+rrt(589)+rrt(596)+rrt(603)+rrt(610)+rrt(617)+rrt(624)+rrt(631)+  2.d0 * rrt(638)+rrt(639)+rrt(643)+rrt(645)&
             +  2.d0 * rrt(647)+rrt(648)+rrt(652)+rrt(654)+  2.d0 * rrt(656)+rrt(657)+rrt(661)+rrt(663)+  2.d0 * rrt(665)+rrt(666)&
             +rrt(670)+rrt(672)+  2.d0 * rrt(674)+rrt(675)+rrt(679)+rrt(681)+  2.d0 * rrt(683)+rrt(684)+rrt(688)+rrt(690)&
             +  2.d0 * rrt(692)+rrt(693)+rrt(697)+rrt(699)+rrt(701)+rrt(708)+rrt(712)+rrt(717)+rrt(758)+rrt(765)+rrt(772)+rrt(779)&
             +rrt(786) 
  ydot(24) = +rrt(142)+rrt(167)+rrt(173)+rrt(250)-rrt(284)-rrt(285)-rrt(286)-rrt(287)-rrt(288)+rrt(291)-rrt(293)-rrt(296)-rrt(297)&
             -rrt(486) 
  ydot(25) = +rrt(168)+rrt(253)-rrt(289)-rrt(290)-rrt(291)-rrt(292)-rrt(293)-rrt(294)-rrt(295)-rrt(298) 
  ydot(26) = +rrt(154)-rrt(180)-rrt(182)-rrt(472)-rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(478)-rrt(479)-rrt(480)+rrt(486)&
             +rrt(494)+rrt(513)-rrt(523)-rrt(524)-rrt(525)-rrt(533)-rrt(534)-rrt(535)-rrt(551)-rrt(589)-rrt(596)-rrt(603)-rrt(610)&
             -rrt(617)-rrt(624)-rrt(631)-rrt(701)-rrt(712)-rrt(717)-rrt(722)-rrt(727)-rrt(758)-rrt(765)-rrt(772)-rrt(779)-rrt(786) 
  ydot(27) = +rrt(156)+rrt(157)-rrt(166)-rrt(167)-rrt(168)+rrt(293)+rrt(369)+rrt(478)-rrt(491)-rrt(492)-rrt(493)-rrt(494)-rrt(495)&
             -rrt(496)-rrt(497)+rrt(505)+rrt(510)+rrt(525)-rrt(529)-rrt(530)-rrt(550)-rrt(590)-rrt(597)-rrt(604)-rrt(611)-rrt(618)&
             -rrt(625)-rrt(632)-rrt(638)-rrt(647)-rrt(656)-rrt(665)-rrt(674)-rrt(683)-rrt(692)-rrt(702)-rrt(713)-rrt(718)-rrt(723)&
             -rrt(740)-rrt(749)-rrt(750)-rrt(751)-rrt(759)-rrt(766)-rrt(773)-rrt(780)-rrt(787) 
  ydot(28) = -rrt(174)-rrt(503)-rrt(504)-rrt(505)-rrt(506)-rrt(507)+rrt(523)+rrt(530)-rrt(639)-rrt(648)-rrt(657)-rrt(666)-rrt(675)&
             -rrt(684)-rrt(693)-rrt(708) 
  ydot(29) = -rrt(175)+rrt(266)+rrt(267)-rrt(510)-rrt(511)-rrt(512)-rrt(513)-rrt(514)+rrt(529)-rrt(640)-rrt(649)-rrt(658)-rrt(667)&
             -rrt(676)-rrt(685)-rrt(694)-rrt(709)-rrt(741)-rrt(752)-rrt(753)-rrt(754) 
  ydot(30) = -rrt(044)-rrt(045)-rrt(046)-rrt(047)-rrt(048)-rrt(049)+rrt(098)-rrt(102)+rrt(106)-rrt(110)-rrt(143)-rrt(144)-rrt(145)&
             -rrt(146)-rrt(147)-rrt(148)-rrt(158)+  2.d0 * rrt(178)+rrt(179)-rrt(185)+rrt(187)-rrt(191)-rrt(195)+rrt(197)-rrt(201)&
             +rrt(203)+  2.d0 * rrt(206)+rrt(209)+  2.d0 * rrt(210)+  2.d0 * rrt(211)+rrt(212)+rrt(213)+rrt(214)+  2.d0 * rrt(215)&
             +rrt(217)+rrt(223)+rrt(239)+rrt(241)+rrt(242)-rrt(243)-rrt(244)-rrt(245)-rrt(246)-rrt(247)-rrt(248)-rrt(249)-rrt(285)&
             -rrt(294)+rrt(299)+rrt(301)+rrt(302)+rrt(303)+  2.d0 * rrt(304)+rrt(305)+rrt(306)+rrt(308)+  2.d0 * rrt(312)+rrt(313)&
             -rrt(314)-rrt(318)-rrt(319)+rrt(321)+  2.d0 * rrt(322)+rrt(323)+rrt(325)-rrt(329)+  2.d0 * rrt(336)+rrt(337)-rrt(341)&
             +rrt(344)+rrt(347)+rrt(349)+rrt(351)+rrt(352)-rrt(353)+rrt(356)-rrt(357)+rrt(358)-  2.d0 * rrt(361)-rrt(362)+rrt(363)&
             +rrt(365)+rrt(366)-rrt(367)+rrt(368)-rrt(389)-rrt(390)-rrt(391)-rrt(392)-rrt(393)+rrt(399)+rrt(400)+rrt(401)+rrt(402)&
             +rrt(416)+rrt(417)+rrt(418)+rrt(419)+rrt(420)+rrt(427)+rrt(428)+rrt(429)+rrt(430)+rrt(431)-rrt(437)-rrt(438)-rrt(439)&
             -rrt(440)-rrt(441)-rrt(473)-rrt(474)-rrt(475)+rrt(476)-rrt(482)+rrt(483)-rrt(491)+rrt(500)+rrt(502)-rrt(503)-rrt(504)&
             -rrt(511)+rrt(515)+rrt(516)+  2.d0 * rrt(517)+  2.d0 * rrt(518)+  2.d0 * rrt(520)-rrt(522)-rrt(531)+rrt(552)+rrt(561)&
             +rrt(562)+rrt(563)+rrt(564)+rrt(565)+rrt(567)+rrt(569)-rrt(571)+rrt(574)+rrt(579)+rrt(580)+rrt(581)+  2.d0 * rrt(582)&
             +  2.d0 * rrt(583)+  2.d0 * rrt(584)+rrt(585)-rrt(586)-rrt(588)+rrt(592)+rrt(596)+rrt(597)+rrt(598)+  2.d0 * rrt(599)&
             +rrt(600)+rrt(601)+rrt(602)+rrt(606)+rrt(613)+rrt(620)+rrt(627)+rrt(634)+  2.d0 * rrt(642)+rrt(645)+rrt(646)+rrt(647)&
             +rrt(648)+rrt(649)+rrt(650)+  3.d0 * rrt(651)+rrt(652)+rrt(653)+  2.d0 * rrt(654)+  2.d0 * rrt(655)+  2.d0 * rrt(660)&
             +rrt(663)+rrt(664)+  2.d0 * rrt(669)+rrt(672)+rrt(673)+  2.d0 * rrt(678)+rrt(681)+rrt(682)+  2.d0 * rrt(687)+rrt(690)&
             +rrt(691)+  2.d0 * rrt(696)+rrt(699)+rrt(700)+  2.d0 * rrt(701)+  2.d0 * rrt(702)+  2.d0 * rrt(703)+  3.d0 * rrt(704)&
             +  2.d0 * rrt(705)+  2.d0 * rrt(706)+  2.d0 * rrt(707)+  2.d0 * rrt(708)+  2.d0 * rrt(709)+  4.d0 * rrt(710)&
             +  3.d0 * rrt(711)+rrt(715)+rrt(717)+rrt(718)+rrt(719)+  2.d0 * rrt(720)+rrt(721)+rrt(724)+rrt(730)+rrt(761)+rrt(768)
  ydot(30) = ydot(30) &
             +rrt(775)+rrt(782)+rrt(789) 
  ydot(31) = +rrt(044)+rrt(045)-rrt(098)+rrt(099)+rrt(102)-rrt(103)-rrt(106)+rrt(107)+rrt(110)-rrt(111) 
  ydot(32) = +rrt(046)+rrt(047)-rrt(099)+rrt(100)+rrt(103)-rrt(104)-rrt(107)+rrt(108)+rrt(111)-rrt(112) 
  ydot(33) = +rrt(048)-rrt(100)+rrt(101)+rrt(104)-rrt(105)-rrt(108)+rrt(109)+rrt(112)-rrt(113) 
  ydot(34) = +rrt(049)-rrt(101)+rrt(105)-rrt(109)+rrt(113) 
  ydot(35) = +rrt(143)-rrt(149)-rrt(159)-rrt(202)-rrt(210)-rrt(239)+rrt(240)-rrt(299)-rrt(300)-rrt(301)-rrt(302)-rrt(303)-rrt(304)&
             -  2.d0 * rrt(305)+rrt(306)+rrt(307)+rrt(309)+rrt(310)+rrt(311)+rrt(318)-rrt(331)-rrt(332)-rrt(333)-rrt(517)-rrt(556)&
             -rrt(583) 
  ydot(36) = +rrt(144)-rrt(203)-rrt(211)-rrt(240)-rrt(241)+rrt(244)+rrt(305)-rrt(307)-rrt(308)-rrt(309)-rrt(310)-rrt(311)-rrt(312)&
             +  2.d0 * rrt(314)+rrt(315)+rrt(319)+rrt(332)-rrt(518)-rrt(584) 
  ydot(37) = +rrt(145)-rrt(242)-rrt(313)-rrt(314)-rrt(315)+rrt(331) 
  ydot(38) = +  2.d0 * rrt(146)+rrt(147)+rrt(148)+  2.d0 * rrt(149)-rrt(150)-rrt(151)-rrt(155)+  2.d0 * rrt(169)+rrt(170)+rrt(171)&
             +rrt(172)+rrt(173)+rrt(176)+rrt(177)+rrt(181)+rrt(183)+rrt(185)+rrt(188)-rrt(190)-rrt(197)+rrt(203)+rrt(204)+rrt(205)&
             -rrt(207)-rrt(215)-rrt(221)-rrt(222)-rrt(223)-rrt(224)+  2.d0 * rrt(243)+  2.d0 * rrt(245)+rrt(246)+  2.d0 * rrt(247)&
             +rrt(248)+rrt(249)-rrt(250)-rrt(251)+rrt(257)+rrt(265)-rrt(284)+rrt(285)+rrt(286)+rrt(294)+rrt(295)+rrt(300)-rrt(306)&
             -rrt(308)+rrt(312)-rrt(313)+rrt(316)+rrt(317)+rrt(318)+rrt(319)+rrt(320)+  2.d0 * rrt(321)+rrt(327)+  3.d0 * rrt(329)&
             +rrt(330)+rrt(331)+  3.d0 * rrt(333)+rrt(334)+rrt(337)+rrt(338)+rrt(340)+rrt(341)+  2.d0 * rrt(342)+rrt(343)-rrt(346)&
             -rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)+rrt(353)+rrt(355)+rrt(357)+rrt(361)-rrt(370)+  2.d0 * rrt(389)&
             +  2.d0 * rrt(390)+  2.d0 * rrt(391)+  2.d0 * rrt(392)+  2.d0 * rrt(393)+rrt(394)+rrt(395)+rrt(396)+rrt(397)+rrt(398)&
             +rrt(399)+rrt(400)+rrt(401)+rrt(402)+rrt(403)+rrt(404)+rrt(405)+rrt(406)+rrt(407)+rrt(408)+rrt(409)+rrt(410)+rrt(411)&
             +rrt(412)+rrt(413)+rrt(414)+rrt(415)-  2.d0 * rrt(427)-  2.d0 * rrt(428)-  2.d0 * rrt(429)-  2.d0 * rrt(430)&
             -  2.d0 * rrt(431)-rrt(432)-rrt(433)-rrt(434)-rrt(435)-rrt(436)-rrt(437)-rrt(438)-rrt(439)-rrt(440)-rrt(441)-rrt(442)&
             -rrt(443)-rrt(444)-rrt(445)-rrt(446)-rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(472)+rrt(474)+rrt(478)+rrt(482)+rrt(484)&
             +rrt(486)+rrt(488)+rrt(490)-rrt(492)+rrt(493)+rrt(499)-rrt(512)-rrt(519)-rrt(524)-rrt(527)+rrt(556)+rrt(557)+rrt(558)&
             +rrt(560)-rrt(561)-rrt(565)+rrt(566)-rrt(581)-rrt(582)+rrt(589)+rrt(590)+  2.d0 * rrt(591)+rrt(592)+rrt(593)+rrt(594)&
             +rrt(595)+rrt(598)+rrt(605)+rrt(612)+rrt(619)+rrt(626)+rrt(633)+rrt(638)+rrt(639)+rrt(640)+  3.d0 * rrt(641)+rrt(642)&
             +  2.d0 * rrt(643)+  2.d0 * rrt(644)+rrt(645)+rrt(646)+  2.d0 * rrt(650)+rrt(652)+rrt(653)+  2.d0 * rrt(659)+rrt(661)&
             +rrt(662)+  2.d0 * rrt(668)+rrt(670)+rrt(671)+  2.d0 * rrt(677)+rrt(679)+rrt(680)+  2.d0 * rrt(686)+rrt(688)+rrt(689)&
             +  2.d0 * rrt(695)+rrt(697)+rrt(698)+rrt(703)+rrt(712)+rrt(713)+  2.d0 * rrt(714)+rrt(715)+rrt(716)+rrt(719)+rrt(731)&
             +rrt(760)+rrt(767)+rrt(774)+rrt(781)+rrt(788) 
  ydot(39) = +rrt(147)+rrt(150)+rrt(170)+rrt(246)+rrt(248)+rrt(284)+rrt(304)+rrt(308)-rrt(316)-rrt(317)-rrt(318)-rrt(319)-rrt(320)&
             -rrt(321)-rrt(322)-rrt(323)-rrt(324)-rrt(325)+rrt(326)+rrt(328)+rrt(332)+rrt(335)+rrt(337)+rrt(339) 
  ydot(40) = +rrt(148)+rrt(151)+rrt(171)+rrt(249)+rrt(251)+rrt(313)-rrt(326)-rrt(327)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332)&
             -rrt(333)-rrt(334)-rrt(335)-rrt(336)-rrt(337)-rrt(338)-rrt(339) 
  ydot(41) = -rrt(162)-rrt(187)-rrt(188)-rrt(192)+rrt(201)+rrt(202)-rrt(206)+rrt(207)+rrt(224)+rrt(225)+rrt(230)-rrt(304)-rrt(306)&
             -rrt(312)-rrt(321)-rrt(322)-rrt(336)-rrt(337)-rrt(358)+rrt(361)+rrt(362)-rrt(365)+rrt(367)-rrt(399)-rrt(400)-rrt(401)&
             -rrt(402)+rrt(437)+rrt(438)+rrt(439)+rrt(440)+rrt(441)-rrt(476)-rrt(483)-rrt(493)+rrt(501)+rrt(519)-rrt(557)-rrt(562)&
             +rrt(568)+rrt(570)-rrt(574)+rrt(603)+rrt(604)+rrt(605)+rrt(606)+rrt(607)+rrt(608)+rrt(609)+rrt(656)+rrt(657)+rrt(658)&
             +rrt(659)+rrt(660)+rrt(661)+rrt(662)+rrt(663)+rrt(664)+rrt(725)+rrt(728)+rrt(732)+rrt(758)+rrt(759)+rrt(760)+rrt(761)&
             +rrt(762)+rrt(763)+rrt(764) 
  ydot(42) = +rrt(155)-rrt(181)-rrt(183)+rrt(472)+rrt(475)+rrt(479)-rrt(481)-rrt(482)-rrt(483)-rrt(484)-rrt(485)-rrt(486)-rrt(487)&
             -rrt(488)-rrt(489)-rrt(490)+rrt(512)-rrt(526)-rrt(527)-rrt(528)-rrt(591)-rrt(598)-rrt(605)-rrt(612)-rrt(619)-rrt(626)&
             -rrt(633)-rrt(703)-rrt(714)-rrt(719)-rrt(724)-rrt(728)-rrt(760)-rrt(767)-rrt(774)-rrt(781)-rrt(788) 
  ydot(43) = +rrt(158)+rrt(159)-rrt(169)-rrt(170)-rrt(171)+rrt(473)+rrt(482)+rrt(483)+rrt(485)+rrt(489)+rrt(491)+rrt(493)-rrt(498)&
             -rrt(499)-rrt(500)-rrt(501)-rrt(502)+rrt(503)+rrt(511)+rrt(516)+rrt(517)+rrt(518)+rrt(519)+rrt(521)+rrt(527)-rrt(531)&
             -rrt(532)-rrt(552)-rrt(592)-rrt(599)-rrt(606)-rrt(613)-rrt(620)-rrt(627)-rrt(634)-rrt(641)-rrt(650)-rrt(659)-rrt(668)&
             -rrt(677)-rrt(686)-rrt(695)-rrt(704)-rrt(715)-rrt(720)-rrt(725)-rrt(761)-rrt(768)-rrt(775)-rrt(782)-rrt(789) 
  ydot(44) = +rrt(162) 
  ydot(45) = -rrt(178)-rrt(515)-rrt(516)-rrt(517)-rrt(518)-rrt(519)-rrt(520)+rrt(522)+rrt(531)-rrt(642)-rrt(651)-rrt(660)-rrt(669)&
             -rrt(678)-rrt(687)-rrt(696)-rrt(710) 
  ydot(46) = +rrt(185)+rrt(186)+rrt(187)+rrt(189)+rrt(190)-rrt(197)-rrt(198)-rrt(199)-rrt(200)-rrt(201)-rrt(202)-rrt(203)-rrt(204)&
             -rrt(205)-rrt(206)-rrt(556)-rrt(557)-rrt(558)-rrt(559)-rrt(560)+rrt(561)+rrt(582)-rrt(586)-rrt(587)-rrt(589)-rrt(590)&
             -rrt(591)-rrt(592)-rrt(593)-rrt(594)-rrt(595)-rrt(638)-rrt(639)-rrt(640)-rrt(641)-rrt(642)-rrt(643)-rrt(644)-rrt(645)&
             -rrt(646)-rrt(712)-rrt(713)-rrt(714)-rrt(715)-rrt(716)-rrt(722)-rrt(723)-rrt(724)-rrt(725)-rrt(726)-rrt(731) 
  ydot(47) = +rrt(188)+rrt(191)+rrt(195)-rrt(207)-rrt(208)-rrt(209)-rrt(210)-rrt(211)-rrt(212)-rrt(213)-rrt(214)+rrt(556)-rrt(561)&
             -rrt(562)-rrt(563)-rrt(564)+rrt(565)+rrt(571)+rrt(579)+rrt(580)+rrt(583)+rrt(584)-rrt(588)-rrt(596)-rrt(597)-rrt(598)&
             -rrt(599)-rrt(600)-rrt(601)-rrt(602)-rrt(647)-rrt(648)-rrt(649)-rrt(650)-rrt(651)-rrt(652)-rrt(653)-rrt(654)-rrt(655)&
             -rrt(717)-rrt(718)-rrt(719)-rrt(720)-rrt(721)-rrt(727)-rrt(728)-rrt(729)-rrt(730) 
  ydot(48) = +rrt(192)-rrt(215)-rrt(217)-rrt(225)-rrt(230)+rrt(557)+rrt(562)-rrt(565)-rrt(566)-rrt(567)-rrt(568)-rrt(569)-rrt(570)&
             +rrt(581)+rrt(586)-rrt(603)-rrt(604)-rrt(605)-rrt(606)-rrt(607)-rrt(608)-rrt(609)-rrt(656)-rrt(657)-rrt(658)-rrt(659)&
             -rrt(660)-rrt(661)-rrt(662)-rrt(663)-rrt(664)-rrt(732)-rrt(758)-rrt(759)-rrt(760)-rrt(761)-rrt(762)-rrt(763)-rrt(764) 
  ydot(49) = -rrt(579)-rrt(580)-rrt(581)-rrt(582)-rrt(583)-rrt(584)-rrt(585)+rrt(588)-rrt(701)-rrt(702)-rrt(703)-rrt(704)-rrt(705)&
             -rrt(706)-rrt(707)-rrt(708)-rrt(709)-rrt(710)-rrt(711) 
  ydot(50) = -rrt(160)+rrt(177)-rrt(186)+rrt(189)-rrt(193)+rrt(198)-rrt(199)+rrt(217)+rrt(218)+  2.d0 * rrt(219)+rrt(220)&
             +  2.d0 * rrt(222)+rrt(223)+rrt(224)+rrt(226)+rrt(231)+rrt(250)+rrt(256)+rrt(257)-rrt(265)+rrt(285)-rrt(286)+rrt(287)&
             +rrt(294)-rrt(295)+rrt(300)-rrt(323)+  2.d0 * rrt(324)-rrt(340)+rrt(341)+  2.d0 * rrt(345)+rrt(346)-rrt(347)-rrt(348)&
             +  2.d0 * rrt(350)+rrt(351)-  2.d0 * rrt(354)-  2.d0 * rrt(355)-  2.d0 * rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(360)&
             +rrt(362)+  2.d0 * rrt(363)+rrt(364)+rrt(366)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-rrt(398)+rrt(407)+rrt(408)+rrt(409)&
             +rrt(410)+rrt(416)+rrt(417)+rrt(418)+rrt(419)+rrt(420)+rrt(432)+rrt(433)+rrt(434)+rrt(435)+rrt(436)-rrt(443)-rrt(444)&
             -rrt(445)+rrt(475)-rrt(477)-rrt(478)-rrt(479)-rrt(484)-rrt(485)+rrt(487)-rrt(495)+rrt(498)-rrt(500)-rrt(506)-rrt(507)&
             -rrt(508)-rrt(509)-rrt(514)-rrt(520)+rrt(559)-rrt(566)-rrt(567)+rrt(571)+rrt(572)+rrt(575)-rrt(578)-rrt(585)-rrt(587)&
             +rrt(593)+rrt(600)+rrt(607)+rrt(610)+rrt(611)+rrt(612)+rrt(613)+  2.d0 * rrt(614)+rrt(615)+rrt(616)+rrt(621)+rrt(628)&
             +rrt(635)+rrt(665)+rrt(666)+rrt(667)+rrt(668)+rrt(669)+rrt(670)+rrt(671)+rrt(672)+rrt(673)+rrt(705)+rrt(716)+rrt(721)&
             +rrt(722)+rrt(733)+rrt(762)+rrt(765)+rrt(766)+rrt(767)+rrt(768)+  2.d0 * rrt(769)+rrt(770)+rrt(771)+rrt(776)+rrt(783)&
             +rrt(790) 
  ydot(51) = -rrt(161)-rrt(194)+rrt(200)+rrt(216)+rrt(227)+rrt(232)-rrt(256)-rrt(287)-rrt(324)-rrt(325)+rrt(343)-rrt(349)-rrt(350)&
             +rrt(353)+rrt(355)-rrt(359)-rrt(403)-rrt(404)-rrt(405)-rrt(406)+rrt(442)-rrt(480)-rrt(487)-rrt(488)-rrt(489)-rrt(496)&
             -rrt(497)+rrt(509)-rrt(559)-rrt(560)-rrt(573)+rrt(594)+rrt(601)+rrt(608)+rrt(615)+rrt(617)+rrt(618)+rrt(619)+rrt(620)&
             +rrt(621)+  2.d0 * rrt(622)+rrt(623)+rrt(629)+rrt(636)+rrt(674)+rrt(675)+rrt(676)+rrt(677)+rrt(678)+rrt(679)+rrt(680)&
             +rrt(681)+rrt(682)+rrt(706)+rrt(723)+rrt(736)+rrt(763)+rrt(770)+rrt(772)+rrt(773)+rrt(774)+rrt(775)+rrt(776)&
             +  2.d0 * rrt(777)+rrt(778)+rrt(784)+rrt(791) 
  ydot(52) = -rrt(189)+rrt(199)+rrt(208)+rrt(220)+rrt(221)+rrt(228)+rrt(233)-rrt(257)-rrt(342)-rrt(343)-rrt(344)-rrt(345)+rrt(348)&
             -rrt(351)+rrt(352)+rrt(354)+rrt(357)+rrt(358)+rrt(359)+  2.d0 * rrt(360)-rrt(362)-  2.d0 * rrt(363)-  2.d0 * rrt(364)&
             -rrt(365)+rrt(367)+  2.d0 * rrt(368)-rrt(407)-rrt(408)-rrt(409)-rrt(410)+rrt(411)+rrt(412)+rrt(413)+rrt(414)+rrt(415)&
             +rrt(421)+rrt(443)+rrt(444)+rrt(445)-rrt(446)-rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(451)-rrt(490)-rrt(501)-rrt(502)&
             +rrt(508)-rrt(558)-rrt(563)-rrt(568)-rrt(569)-rrt(572)-rrt(575)+rrt(576)+  2.d0 * rrt(577)+rrt(578)+rrt(595)+rrt(602)&
             +rrt(609)+rrt(616)+rrt(623)+rrt(624)+rrt(625)+rrt(626)+rrt(627)+rrt(628)+rrt(629)+  2.d0 * rrt(630)+rrt(637)+rrt(683)&
             +rrt(684)+rrt(685)+rrt(686)+rrt(687)+rrt(688)+rrt(689)+rrt(690)+rrt(691)+rrt(707)+rrt(726)+rrt(727)+rrt(734)+rrt(764)&
             +rrt(771)+rrt(778)+rrt(779)+rrt(780)+rrt(781)+rrt(782)+rrt(783)+rrt(784)+  2.d0 * rrt(785)+rrt(792) 
  ydot(53) = +rrt(229)+rrt(234)-rrt(352)-rrt(360)+rrt(364)+rrt(365)-rrt(366)-rrt(367)-  2.d0 * rrt(368)-rrt(411)-rrt(412)-rrt(413)&
             -rrt(414)-rrt(415)-rrt(416)-rrt(417)-rrt(418)-rrt(419)-rrt(420)+rrt(421)+rrt(446)+rrt(447)+rrt(448)+rrt(449)+rrt(450)&
             -rrt(451)-rrt(564)-rrt(570)-rrt(576)+rrt(631)+rrt(632)+rrt(633)+rrt(634)+rrt(635)+rrt(636)+rrt(637)+rrt(692)+rrt(693)&
             +rrt(694)+rrt(695)+rrt(696)+rrt(697)+rrt(698)+rrt(699)+rrt(700)+rrt(729)+rrt(735)+rrt(737)+rrt(786)+rrt(787)+rrt(788)&
             +rrt(789)+rrt(790)+rrt(791)+rrt(792) 
  ydot(54) = -rrt(421)+rrt(451)-rrt(577) 
  ydot(55) = +rrt(160)-rrt(172)-rrt(173)+rrt(370)+rrt(474)+rrt(476)+rrt(477)+rrt(480)+rrt(481)+rrt(484)+rrt(487)+rrt(492)+rrt(495)&
             +rrt(497)+rrt(498)+rrt(499)+rrt(500)+rrt(501)+rrt(506)+rrt(508)+rrt(509)+rrt(514)+rrt(520)+rrt(524)+rrt(526)+rrt(528)&
             -rrt(593)-rrt(600)-rrt(607)-rrt(614)-rrt(621)-rrt(628)-rrt(635)-rrt(643)-rrt(652)-rrt(661)-rrt(670)-rrt(679)-rrt(688)&
             -rrt(697)-rrt(705)-rrt(716)-rrt(721)-rrt(726)-rrt(729)-rrt(762)-rrt(769)-rrt(776)-rrt(783)-rrt(790) 
  ydot(56) = +rrt(161)-rrt(176)+rrt(488)+rrt(496)+rrt(507)-rrt(509)-rrt(594)-rrt(601)-rrt(608)-rrt(615)-rrt(622)-rrt(629)-rrt(636)&
             -rrt(644)-rrt(653)-rrt(662)-rrt(671)-rrt(680)-rrt(689)-rrt(698)-rrt(706)-rrt(763)-rrt(770)-rrt(777)-rrt(784)-rrt(791) 
  ydot(57) = -rrt(177)+rrt(490)+rrt(502)+rrt(504)-rrt(508)-rrt(595)-rrt(602)-rrt(609)-rrt(616)-rrt(623)-rrt(630)-rrt(637)-rrt(645)&
             -rrt(654)-rrt(663)-rrt(672)-rrt(681)-rrt(690)-rrt(699)-rrt(707)-rrt(764)-rrt(771)-rrt(778)-rrt(785)-rrt(792) 
  ydot(58) = +rrt(193)-rrt(216)-rrt(221)-rrt(226)-rrt(231)+rrt(559)-rrt(571)-rrt(572)-rrt(573)-rrt(610)-rrt(611)-rrt(612)-rrt(613)&
             -rrt(614)-rrt(615)-rrt(616)-rrt(665)-rrt(666)-rrt(667)-rrt(668)-rrt(669)-rrt(670)-rrt(671)-rrt(672)-rrt(673)-rrt(733)&
             -rrt(765)-rrt(766)-rrt(767)-rrt(768)-rrt(769)-rrt(770)-rrt(771) 
  ydot(59) = +rrt(194)-rrt(218)-rrt(222)-rrt(227)-rrt(232)+rrt(560)-rrt(617)-rrt(618)-rrt(619)-rrt(620)-rrt(621)-rrt(622)-rrt(623)&
             -rrt(674)-rrt(675)-rrt(676)-rrt(677)-rrt(678)-rrt(679)-rrt(680)-rrt(681)-rrt(682)-rrt(736)-rrt(772)-rrt(773)-rrt(774)&
             -rrt(775)-rrt(776)-rrt(777)-rrt(778) 
  ydot(60) = -rrt(219)-rrt(223)-rrt(228)-rrt(233)+rrt(558)+rrt(563)+rrt(567)+rrt(568)+rrt(572)+rrt(573)-rrt(574)-rrt(575)-rrt(576)&
             -rrt(577)+rrt(578)+rrt(587)-rrt(624)-rrt(625)-rrt(626)-rrt(627)-rrt(628)-rrt(629)-rrt(630)-rrt(683)-rrt(684)-rrt(685)&
             -rrt(686)-rrt(687)-rrt(688)-rrt(689)-rrt(690)-rrt(691)-rrt(734)-rrt(779)-rrt(780)-rrt(781)-rrt(782)-rrt(783)-rrt(784)&
             -rrt(785) 
  ydot(61) = -rrt(220)-rrt(224)-rrt(229)-rrt(234)+rrt(564)+rrt(566)+rrt(569)+rrt(570)+rrt(574)+rrt(575)+rrt(576)+rrt(577)-rrt(578)&
             +rrt(585)-rrt(631)-rrt(632)-rrt(633)-rrt(634)-rrt(635)-rrt(636)-rrt(637)-rrt(692)-rrt(693)-rrt(694)-rrt(695)-rrt(696)&
             -rrt(697)-rrt(698)-rrt(699)-rrt(700)-rrt(735)-rrt(737)-rrt(786)-rrt(787)-rrt(788)-rrt(789)-rrt(790)-rrt(791)-rrt(792) 
  ydot(62) = -rrt(269)-rrt(296)-rrt(298)+rrt(371)+rrt(372)+rrt(373)-rrt(374)+rrt(377)+rrt(380)+rrt(455)+rrt(456)-rrt(460)-rrt(461)&
             -rrt(462)-rrt(469)-rrt(470)-rrt(471)-rrt(533)+rrt(535)+rrt(537)-rrt(538)+rrt(539)-rrt(541)-rrt(542)-rrt(546)+rrt(739)&
             +rrt(742)+rrt(743)+rrt(744)+rrt(745)+  2.d0 * rrt(746)+  2.d0 * rrt(747)+  2.d0 * rrt(748)+rrt(755)+rrt(756)+rrt(757) 
  ydot(63) = +rrt(152)+  2.d0 * rrt(153)+  2.d0 * rrt(184)+  2.d0 * rrt(269)+rrt(296)+rrt(298)-rrt(371)-rrt(372)-rrt(373)+rrt(374)&
             +rrt(375)+  2.d0 * rrt(376)+  2.d0 * rrt(379)-rrt(382)-  2.d0 * rrt(455)-  2.d0 * rrt(456)-rrt(457)-rrt(458)-rrt(459)&
             -rrt(463)-rrt(464)-rrt(465)-rrt(466)-rrt(467)-rrt(468)+rrt(533)+rrt(536)-rrt(537)+rrt(538)+rrt(540)+rrt(542)+rrt(546)&
             +rrt(737)+  3.d0 * rrt(738)+  2.d0 * rrt(739)+rrt(740)+rrt(741)+rrt(743)+rrt(744)+rrt(745)+rrt(749)+rrt(750)+rrt(751)&
             +rrt(752)+rrt(753)+rrt(754) 
  ydot(64) = +rrt(164)-rrt(536)+rrt(537) 
  ydot(65) = +rrt(165)-rrt(537)-rrt(538)-rrt(539)-rrt(540)-rrt(738)-rrt(743)-rrt(744)-rrt(745) 
  ydot(66) = +rrt(538)+rrt(541)-rrt(739)-rrt(746)-rrt(747)-rrt(748) 
  ydot(67) = +rrt(196)-rrt(738)-rrt(739)-rrt(740)-rrt(741)-rrt(742)-rrt(743)-rrt(744)-rrt(745)-rrt(746)-rrt(747)-rrt(748)-rrt(749)&
             -rrt(750)-rrt(751)-rrt(752)-rrt(753)-rrt(754)-rrt(755)-rrt(756)-rrt(757) 
  ydot(68) = +rrt(553) 
  ydot(69) = -rrt(553) 
  ydot(70) = -rrt(152)-rrt(153)-rrt(163)-rrt(164)-rrt(165)-rrt(196)-rrt(297)-rrt(373)+rrt(374)+rrt(381)+rrt(382)+rrt(383)+rrt(466)&
             +rrt(467)+rrt(468)+rrt(469)+rrt(470)+rrt(471)-rrt(534)-rrt(535)-rrt(536)-rrt(539)-rrt(544)-rrt(547)-rrt(548)-rrt(549)&
             -rrt(550)-rrt(551)-rrt(552)-rrt(553)-rrt(554)+rrt(555)+rrt(730)+rrt(731)+rrt(732)+rrt(733)+rrt(734)+rrt(735)+rrt(736)&
             +rrt(737) 
  ydot(71) = +rrt(152)+rrt(164)+rrt(184)+rrt(196)+rrt(297)-rrt(372)+rrt(373)-rrt(374)-rrt(376)-rrt(377)+rrt(378)-rrt(381)-rrt(382)&
             -  2.d0 * rrt(383)+rrt(460)+rrt(461)+rrt(462)+rrt(463)+rrt(464)+rrt(465)-rrt(466)-rrt(467)-rrt(468)+rrt(547)+rrt(554)&
             -rrt(555) 
  ydot(72) = +rrt(153)+rrt(165)+rrt(296)+rrt(297)+rrt(298)-rrt(371)+rrt(372)-rrt(375)-  2.d0 * rrt(378)-  2.d0 * rrt(379)&
             -  2.d0 * rrt(380)-rrt(381)+rrt(383)+rrt(457)+rrt(458)+rrt(459)-rrt(463)-rrt(464)-rrt(465)-rrt(469)-rrt(470)-rrt(471)&
             +rrt(534)+rrt(548) 
  ydot(73) = +rrt(163)+rrt(536)+rrt(539)+rrt(546)+rrt(547)+rrt(550)+rrt(551)+rrt(552)-rrt(554)-rrt(555)-rrt(730)-rrt(731)-rrt(732)&
             -rrt(733)-rrt(734)-rrt(735)-rrt(736) 
  ydot(74) = +rrt(534)+rrt(542)-rrt(546)-rrt(547)-rrt(548)+rrt(555) 
  ydot(75) = +rrt(533)-rrt(541)-rrt(542)-rrt(544)-rrt(545) 
  ydot(76) = +rrt(535)+rrt(540)+rrt(545)-rrt(549)-rrt(742)-rrt(755)-rrt(756)-rrt(757) 
  ydot(77) = -rrt(184)+rrt(544)+rrt(548)+rrt(549)+rrt(553)+rrt(554)-rrt(737) 
  ydot(78) = -rrt(179)+rrt(515)-rrt(521)-rrt(522)+rrt(532)-rrt(646)-rrt(655)-rrt(664)-rrt(673)-rrt(682)-rrt(691)-rrt(700)-rrt(711) 
  ydot(79) = +rrt(154)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+rrt(159)+rrt(160)+rrt(161)+rrt(162)+rrt(163)+rrt(164)+rrt(165)-rrt(166)&
             -rrt(167)-rrt(168)-rrt(169)-rrt(170)-rrt(171)-rrt(172)-rrt(173)-rrt(174)-rrt(175)-rrt(176)-rrt(177)-rrt(178)-rrt(179)&
             -rrt(180)-rrt(181)-rrt(182)-rrt(183)-rrt(184)-rrt(185)-rrt(186)-rrt(187)-rrt(188)-rrt(189)-rrt(190)-rrt(191)-rrt(192)&
             -rrt(193)-rrt(194)-rrt(195)-rrt(196)+rrt(197)+rrt(198)+rrt(199)+rrt(200)+rrt(201)+rrt(202)+rrt(203)+rrt(204)+rrt(205)&
             +rrt(206)+rrt(207)+rrt(208)+rrt(209)+rrt(210)+rrt(211)+rrt(212)+rrt(213)+rrt(214)+rrt(215)+rrt(216)+rrt(217)+rrt(218)&
             +rrt(219)+rrt(220)+rrt(221)+rrt(222)+rrt(223)+rrt(224)+rrt(225)+rrt(226)+rrt(227)+rrt(228)+rrt(229)+rrt(230)+rrt(231)&
             +rrt(232)+rrt(233)+rrt(234)+rrt(266)+rrt(267)+rrt(293)+rrt(369)+rrt(370) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(80) = 0.0d0
  if( lgas_heating ) then
    ydot(80) = +3.377D+03*rrt(050)+3.470D+03*rrt(051)+3.365D+03*rrt(052)+3.365D+03*rrt(053)+3.481D+03*rrt(054)+3.365D+03*rrt(055)&
               +3.481D+03*rrt(056)+3.365D+03*rrt(057)-3.376D+03*rrt(058)-3.468D+03*rrt(059)-3.364D+03*rrt(060)-3.364D+03*rrt(061)&
               -3.480D+03*rrt(062)-3.364D+03*rrt(063)-3.480D+03*rrt(064)-3.364D+03*rrt(065)+3.377D+03*rrt(066)+3.470D+03*rrt(067)&
               +3.365D+03*rrt(068)+3.365D+03*rrt(069)+3.481D+03*rrt(070)+3.365D+03*rrt(071)+3.481D+03*rrt(072)+3.365D+03*rrt(073)&
               -3.376D+03*rrt(074)-3.468D+03*rrt(075)-3.364D+03*rrt(076)-3.364D+03*rrt(077)-3.480D+03*rrt(078)-3.364D+03*rrt(079)&
               -3.480D+03*rrt(080)-3.364D+03*rrt(081)+3.377D+03*rrt(082)+3.470D+03*rrt(083)+3.365D+03*rrt(084)+3.365D+03*rrt(085)&
               +3.481D+03*rrt(086)+3.365D+03*rrt(087)+3.481D+03*rrt(088)+3.365D+03*rrt(089)-3.376D+03*rrt(090)-3.468D+03*rrt(091)&
               -3.364D+03*rrt(092)-3.364D+03*rrt(093)-3.480D+03*rrt(094)-3.364D+03*rrt(095)-3.480D+03*rrt(096)-3.364D+03*rrt(097)&
               +2.205D+03*rrt(098)+2.205D+03*rrt(099)+2.205D+03*rrt(100)+2.089D+03*rrt(101)-2.204D+03*rrt(102)-2.204D+03*rrt(103)&
               -2.204D+03*rrt(104)-2.088D+03*rrt(105)+2.205D+03*rrt(106)+2.205D+03*rrt(107)+2.205D+03*rrt(108)+2.089D+03*rrt(109)&
               -2.204D+03*rrt(110)-2.204D+03*rrt(111)-2.204D+03*rrt(112)-2.088D+03*rrt(113)+9.957D+03*rrt(142)+1.026D+04*rrt(146)&
               +1.188D+04*rrt(147)+1.061D+04*rrt(148)+8.065D+04*rrt(169)+5.791D+04*rrt(170)+1.218D+04*rrt(243)+5.268D+04*rrt(244)&
               +2.588D+04*rrt(245)+1.184D+04*rrt(246)+6.858D+04*rrt(247)+4.236D+04*rrt(248)+2.286D+04*rrt(249)+2.588D+04*rrt(251)&
               +4.642D+04*rrt(258)+1.520D+04*rrt(259)+2.623D+04*rrt(316)+2.623D+04*rrt(317)+1.485D+04*rrt(318)+7.311D+03*rrt(319)&
               +2.623D+04*rrt(320)+1.950D+04*rrt(326)+4.572D+04*rrt(327)+1.950D+04*rrt(328)-1.369D+04*rrt(329)+4.572D+04*rrt(330)&
               +4.874D+03*rrt(331)+1.195D+04*rrt(332)-2.320D+03*rrt(333)+5.709D+04*rrt(334)+1.950D+04*rrt(335)+4.572D+04*rrt(336)&
               +4.572D+04*rrt(337)+4.572D+04*rrt(338)+1.950D+04*rrt(339)+1.519D+05*rrt(589)+1.641D+05*rrt(590)+1.413D+05*rrt(591)&
               +1.234D+05*rrt(592)+9.081D+04*rrt(593)+1.329D+05*rrt(594)+9.454D+04*rrt(595)+1.635D+05*rrt(596)+1.756D+05*rrt(597)&
               +1.528D+05*rrt(598)+1.348D+05*rrt(599)+1.023D+05*rrt(600)+1.444D+05*rrt(601)+1.060D+05*rrt(602)+1.442D+05*rrt(603)
    ydot(80) = ydot(80) &
               +1.564D+05*rrt(604)+1.337D+05*rrt(605)+1.157D+05*rrt(606)+8.310D+04*rrt(607)+1.252D+05*rrt(608)+8.684D+04*rrt(609)&
               +1.661D+05*rrt(617)+1.782D+05*rrt(618)+1.555D+05*rrt(619)+1.375D+05*rrt(620)+1.050D+05*rrt(621)+1.470D+05*rrt(622)&
               +1.087D+05*rrt(623)+1.423D+05*rrt(624)+1.545D+05*rrt(625)+1.317D+05*rrt(626)+1.137D+05*rrt(627)+8.113D+04*rrt(628)&
               +1.232D+05*rrt(629)+8.486D+04*rrt(630)+1.230D+05*rrt(631)+1.351D+05*rrt(632)+1.123D+05*rrt(633)+9.438D+04*rrt(634)&
               +6.182D+04*rrt(635)+1.039D+05*rrt(636)+6.555D+04*rrt(637)+1.519D+05*rrt(712)+1.641D+05*rrt(713)+1.413D+05*rrt(714)&
               +1.234D+05*rrt(715)+9.081D+04*rrt(716)+1.635D+05*rrt(717)+1.756D+05*rrt(718)+1.528D+05*rrt(719)+1.348D+05*rrt(720)&
               +1.023D+05*rrt(721)+1.117D+05*rrt(730)+1.002D+05*rrt(731)+9.245D+04*rrt(732)+9.048D+04*rrt(734)+7.117D+04*rrt(735)&
               +1.143D+05*rrt(736)+1.721D+05*rrt(740)+1.702D+05*rrt(743)+1.702D+05*rrt(744)+1.702D+05*rrt(745)+1.721D+05*rrt(749)&
               +1.721D+05*rrt(750)+1.721D+05*rrt(751)+8.175D+04*rrt(755)+8.175D+04*rrt(756)+8.175D+04*rrt(757) 
    ydot(80) = ( ZDPlasKin_cfg(14)/k_B + ydot(80) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(80) = ydot(80) * ZDPlasKin_cfg(13)
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
  double precision              :: ysum
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(80)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(001) * density(79) 
  pd(01,79) = pd(01,79) - rrt(001) * density(01) 
  pd(02,01) = pd(02,01) + rrt(001) * density(79) 
  pd(02,79) = pd(02,79) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) - rrt(002) * density(79) 
  pd(01,79) = pd(01,79) - rrt(002) * density(01) 
  pd(03,01) = pd(03,01) + rrt(002) * density(79) 
  pd(03,79) = pd(03,79) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(79) 
  pd(01,79) = pd(01,79) - rrt(003) * density(01) 
  pd(04,01) = pd(04,01) + rrt(003) * density(79) 
  pd(04,79) = pd(04,79) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(79) 
  pd(01,79) = pd(01,79) - rrt(004) * density(01) 
  pd(05,01) = pd(05,01) + rrt(004) * density(79) 
  pd(05,79) = pd(05,79) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) - rrt(005) * density(79) 
  pd(01,79) = pd(01,79) - rrt(005) * density(01) 
  pd(06,01) = pd(06,01) + rrt(005) * density(79) 
  pd(06,79) = pd(06,79) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(79) 
  pd(01,79) = pd(01,79) - rrt(006) * density(01) 
  pd(07,01) = pd(07,01) + rrt(006) * density(79) 
  pd(07,79) = pd(07,79) + rrt(006) * density(01) 
  pd(01,01) = pd(01,01) - rrt(007) * density(79) 
  pd(01,79) = pd(01,79) - rrt(007) * density(01) 
  pd(08,01) = pd(08,01) + rrt(007) * density(79) 
  pd(08,79) = pd(08,79) + rrt(007) * density(01) 
  pd(01,01) = pd(01,01) - rrt(008) * density(79) 
  pd(01,79) = pd(01,79) - rrt(008) * density(01) 
  pd(09,01) = pd(09,01) + rrt(008) * density(79) 
  pd(09,79) = pd(09,79) + rrt(008) * density(01) 
  pd(01,01) = pd(01,01) - rrt(009) * density(79) 
  pd(01,79) = pd(01,79) - rrt(009) * density(01) 
  pd(10,01) = pd(10,01) + rrt(009) * density(79) 
  pd(10,79) = pd(10,79) + rrt(009) * density(01) 
  pd(01,01) = pd(01,01) - rrt(010) * density(79) 
  pd(01,79) = pd(01,79) - rrt(010) * density(01) 
  pd(11,01) = pd(11,01) + rrt(010) * density(79) 
  pd(11,79) = pd(11,79) + rrt(010) * density(01) 
  pd(01,01) = pd(01,01) - rrt(011) * density(79) 
  pd(01,79) = pd(01,79) - rrt(011) * density(01) 
  pd(12,01) = pd(12,01) + rrt(011) * density(79) 
  pd(12,79) = pd(12,79) + rrt(011) * density(01) 
  pd(01,01) = pd(01,01) - rrt(012) * density(79) 
  pd(01,79) = pd(01,79) - rrt(012) * density(01) 
  pd(13,01) = pd(13,01) + rrt(012) * density(79) 
  pd(13,79) = pd(13,79) + rrt(012) * density(01) 
  pd(01,01) = pd(01,01) - rrt(013) * density(79) 
  pd(01,79) = pd(01,79) - rrt(013) * density(01) 
  pd(14,01) = pd(14,01) + rrt(013) * density(79) 
  pd(14,79) = pd(14,79) + rrt(013) * density(01) 
  pd(01,01) = pd(01,01) - rrt(014) * density(79) 
  pd(01,79) = pd(01,79) - rrt(014) * density(01) 
  pd(15,01) = pd(15,01) + rrt(014) * density(79) 
  pd(15,79) = pd(15,79) + rrt(014) * density(01) 
  pd(01,01) = pd(01,01) - rrt(015) * density(79) 
  pd(01,79) = pd(01,79) - rrt(015) * density(01) 
  pd(16,01) = pd(16,01) + rrt(015) * density(79) 
  pd(16,79) = pd(16,79) + rrt(015) * density(01) 
  pd(01,02) = pd(01,02) + rrt(016) * density(79) 
  pd(01,79) = pd(01,79) + rrt(016) * density(02) 
  pd(02,02) = pd(02,02) - rrt(016) * density(79) 
  pd(02,79) = pd(02,79) - rrt(016) * density(02) 
  pd(02,02) = pd(02,02) - rrt(017) * density(79) 
  pd(02,79) = pd(02,79) - rrt(017) * density(02) 
  pd(03,02) = pd(03,02) + rrt(017) * density(79) 
  pd(03,79) = pd(03,79) + rrt(017) * density(02) 
  pd(02,02) = pd(02,02) - rrt(018) * density(79) 
  pd(02,79) = pd(02,79) - rrt(018) * density(02) 
  pd(04,02) = pd(04,02) + rrt(018) * density(79) 
  pd(04,79) = pd(04,79) + rrt(018) * density(02) 
  pd(02,02) = pd(02,02) - rrt(019) * density(79) 
  pd(02,79) = pd(02,79) - rrt(019) * density(02) 
  pd(05,02) = pd(05,02) + rrt(019) * density(79) 
  pd(05,79) = pd(05,79) + rrt(019) * density(02) 
  pd(01,03) = pd(01,03) + rrt(020) * density(79) 
  pd(01,79) = pd(01,79) + rrt(020) * density(03) 
  pd(03,03) = pd(03,03) - rrt(020) * density(79) 
  pd(03,79) = pd(03,79) - rrt(020) * density(03) 
  pd(02,03) = pd(02,03) + rrt(021) * density(79) 
  pd(02,79) = pd(02,79) + rrt(021) * density(03) 
  pd(03,03) = pd(03,03) - rrt(021) * density(79) 
  pd(03,79) = pd(03,79) - rrt(021) * density(03) 
  pd(03,03) = pd(03,03) - rrt(022) * density(79) 
  pd(03,79) = pd(03,79) - rrt(022) * density(03) 
  pd(04,03) = pd(04,03) + rrt(022) * density(79) 
  pd(04,79) = pd(04,79) + rrt(022) * density(03) 
  pd(03,03) = pd(03,03) - rrt(023) * density(79) 
  pd(03,79) = pd(03,79) - rrt(023) * density(03) 
  pd(05,03) = pd(05,03) + rrt(023) * density(79) 
  pd(05,79) = pd(05,79) + rrt(023) * density(03) 
  pd(03,03) = pd(03,03) - rrt(024) * density(79) 
  pd(03,79) = pd(03,79) - rrt(024) * density(03) 
  pd(06,03) = pd(06,03) + rrt(024) * density(79) 
  pd(06,79) = pd(06,79) + rrt(024) * density(03) 
  pd(01,04) = pd(01,04) + rrt(025) * density(79) 
  pd(01,79) = pd(01,79) + rrt(025) * density(04) 
  pd(04,04) = pd(04,04) - rrt(025) * density(79) 
  pd(04,79) = pd(04,79) - rrt(025) * density(04) 
  pd(02,04) = pd(02,04) + rrt(026) * density(79) 
  pd(02,79) = pd(02,79) + rrt(026) * density(04) 
  pd(04,04) = pd(04,04) - rrt(026) * density(79) 
  pd(04,79) = pd(04,79) - rrt(026) * density(04) 
  pd(03,04) = pd(03,04) + rrt(027) * density(79) 
  pd(03,79) = pd(03,79) + rrt(027) * density(04) 
  pd(04,04) = pd(04,04) - rrt(027) * density(79) 
  pd(04,79) = pd(04,79) - rrt(027) * density(04) 
  pd(04,04) = pd(04,04) - rrt(028) * density(79) 
  pd(04,79) = pd(04,79) - rrt(028) * density(04) 
  pd(05,04) = pd(05,04) + rrt(028) * density(79) 
  pd(05,79) = pd(05,79) + rrt(028) * density(04) 
  pd(04,04) = pd(04,04) - rrt(029) * density(79) 
  pd(04,79) = pd(04,79) - rrt(029) * density(04) 
  pd(06,04) = pd(06,04) + rrt(029) * density(79) 
  pd(06,79) = pd(06,79) + rrt(029) * density(04) 
  pd(04,04) = pd(04,04) - rrt(030) * density(79) 
  pd(04,79) = pd(04,79) - rrt(030) * density(04) 
  pd(07,04) = pd(07,04) + rrt(030) * density(79) 
  pd(07,79) = pd(07,79) + rrt(030) * density(04) 
  pd(01,05) = pd(01,05) + rrt(031) * density(79) 
  pd(01,79) = pd(01,79) + rrt(031) * density(05) 
  pd(05,05) = pd(05,05) - rrt(031) * density(79) 
  pd(05,79) = pd(05,79) - rrt(031) * density(05) 
  pd(02,05) = pd(02,05) + rrt(032) * density(79) 
  pd(02,79) = pd(02,79) + rrt(032) * density(05) 
  pd(05,05) = pd(05,05) - rrt(032) * density(79) 
  pd(05,79) = pd(05,79) - rrt(032) * density(05) 
  pd(03,05) = pd(03,05) + rrt(033) * density(79) 
  pd(03,79) = pd(03,79) + rrt(033) * density(05) 
  pd(05,05) = pd(05,05) - rrt(033) * density(79) 
  pd(05,79) = pd(05,79) - rrt(033) * density(05) 
  pd(04,05) = pd(04,05) + rrt(034) * density(79) 
  pd(04,79) = pd(04,79) + rrt(034) * density(05) 
  pd(05,05) = pd(05,05) - rrt(034) * density(79) 
  pd(05,79) = pd(05,79) - rrt(034) * density(05) 
  pd(01,06) = pd(01,06) + rrt(035) * density(79) 
  pd(01,79) = pd(01,79) + rrt(035) * density(06) 
  pd(06,06) = pd(06,06) - rrt(035) * density(79) 
  pd(06,79) = pd(06,79) - rrt(035) * density(06) 
  pd(03,06) = pd(03,06) + rrt(036) * density(79) 
  pd(03,79) = pd(03,79) + rrt(036) * density(06) 
  pd(06,06) = pd(06,06) - rrt(036) * density(79) 
  pd(06,79) = pd(06,79) - rrt(036) * density(06) 
  pd(04,06) = pd(04,06) + rrt(037) * density(79) 
  pd(04,79) = pd(04,79) + rrt(037) * density(06) 
  pd(06,06) = pd(06,06) - rrt(037) * density(79) 
  pd(06,79) = pd(06,79) - rrt(037) * density(06) 
  pd(01,07) = pd(01,07) + rrt(038) * density(79) 
  pd(01,79) = pd(01,79) + rrt(038) * density(07) 
  pd(07,07) = pd(07,07) - rrt(038) * density(79) 
  pd(07,79) = pd(07,79) - rrt(038) * density(07) 
  pd(04,07) = pd(04,07) + rrt(039) * density(79) 
  pd(04,79) = pd(04,79) + rrt(039) * density(07) 
  pd(07,07) = pd(07,07) - rrt(039) * density(79) 
  pd(07,79) = pd(07,79) - rrt(039) * density(07) 
  pd(01,08) = pd(01,08) + rrt(040) * density(79) 
  pd(01,79) = pd(01,79) + rrt(040) * density(08) 
  pd(08,08) = pd(08,08) - rrt(040) * density(79) 
  pd(08,79) = pd(08,79) - rrt(040) * density(08) 
  pd(01,09) = pd(01,09) + rrt(041) * density(79) 
  pd(01,79) = pd(01,79) + rrt(041) * density(09) 
  pd(09,09) = pd(09,09) - rrt(041) * density(79) 
  pd(09,79) = pd(09,79) - rrt(041) * density(09) 
  pd(01,10) = pd(01,10) + rrt(042) * density(79) 
  pd(01,79) = pd(01,79) + rrt(042) * density(10) 
  pd(10,10) = pd(10,10) - rrt(042) * density(79) 
  pd(10,79) = pd(10,79) - rrt(042) * density(10) 
  pd(01,11) = pd(01,11) + rrt(043) * density(79) 
  pd(01,79) = pd(01,79) + rrt(043) * density(11) 
  pd(11,11) = pd(11,11) - rrt(043) * density(79) 
  pd(11,79) = pd(11,79) - rrt(043) * density(11) 
  pd(30,30) = pd(30,30) - rrt(044) * density(79) 
  pd(30,79) = pd(30,79) - rrt(044) * density(30) 
  pd(31,30) = pd(31,30) + rrt(044) * density(79) 
  pd(31,79) = pd(31,79) + rrt(044) * density(30) 
  pd(30,30) = pd(30,30) - rrt(045) * density(79) 
  pd(30,79) = pd(30,79) - rrt(045) * density(30) 
  pd(31,30) = pd(31,30) + rrt(045) * density(79) 
  pd(31,79) = pd(31,79) + rrt(045) * density(30) 
  pd(30,30) = pd(30,30) - rrt(046) * density(79) 
  pd(30,79) = pd(30,79) - rrt(046) * density(30) 
  pd(32,30) = pd(32,30) + rrt(046) * density(79) 
  pd(32,79) = pd(32,79) + rrt(046) * density(30) 
  pd(30,30) = pd(30,30) - rrt(047) * density(79) 
  pd(30,79) = pd(30,79) - rrt(047) * density(30) 
  pd(32,30) = pd(32,30) + rrt(047) * density(79) 
  pd(32,79) = pd(32,79) + rrt(047) * density(30) 
  pd(30,30) = pd(30,30) - rrt(048) * density(79) 
  pd(30,79) = pd(30,79) - rrt(048) * density(30) 
  pd(33,30) = pd(33,30) + rrt(048) * density(79) 
  pd(33,79) = pd(33,79) + rrt(048) * density(30) 
  pd(30,30) = pd(30,30) - rrt(049) * density(79) 
  pd(30,79) = pd(30,79) - rrt(049) * density(30) 
  pd(34,30) = pd(34,30) + rrt(049) * density(79) 
  pd(34,79) = pd(34,79) + rrt(049) * density(30) 
  pd(01,01) = pd(01,01) + rrt(050) * density(02) 
  pd(01,02) = pd(01,02) + rrt(050) * density(01) 
  pd(02,01) = pd(02,01) - rrt(050) * density(02) 
  pd(02,02) = pd(02,02) - rrt(050) * density(01) 
  pd(02,01) = pd(02,01) + rrt(051) * density(03) 
  pd(02,03) = pd(02,03) + rrt(051) * density(01) 
  pd(03,01) = pd(03,01) - rrt(051) * density(03) 
  pd(03,03) = pd(03,03) - rrt(051) * density(01) 
  pd(03,01) = pd(03,01) + rrt(052) * density(04) 
  pd(03,04) = pd(03,04) + rrt(052) * density(01) 
  pd(04,01) = pd(04,01) - rrt(052) * density(04) 
  pd(04,04) = pd(04,04) - rrt(052) * density(01) 
  pd(04,01) = pd(04,01) + rrt(053) * density(05) 
  pd(04,05) = pd(04,05) + rrt(053) * density(01) 
  pd(05,01) = pd(05,01) - rrt(053) * density(05) 
  pd(05,05) = pd(05,05) - rrt(053) * density(01) 
  pd(05,01) = pd(05,01) + rrt(054) * density(06) 
  pd(05,06) = pd(05,06) + rrt(054) * density(01) 
  pd(06,01) = pd(06,01) - rrt(054) * density(06) 
  pd(06,06) = pd(06,06) - rrt(054) * density(01) 
  pd(06,01) = pd(06,01) + rrt(055) * density(07) 
  pd(06,07) = pd(06,07) + rrt(055) * density(01) 
  pd(07,01) = pd(07,01) - rrt(055) * density(07) 
  pd(07,07) = pd(07,07) - rrt(055) * density(01) 
  pd(07,01) = pd(07,01) + rrt(056) * density(08) 
  pd(07,08) = pd(07,08) + rrt(056) * density(01) 
  pd(08,01) = pd(08,01) - rrt(056) * density(08) 
  pd(08,08) = pd(08,08) - rrt(056) * density(01) 
  pd(08,01) = pd(08,01) + rrt(057) * density(09) 
  pd(08,09) = pd(08,09) + rrt(057) * density(01) 
  pd(09,01) = pd(09,01) - rrt(057) * density(09) 
  pd(09,09) = pd(09,09) - rrt(057) * density(01) 
  pd(01,01) = pd(01,01) - rrt(058) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(058) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(059) * density(02) 
  pd(02,02) = pd(02,02) - rrt(059) * density(01) 
  pd(03,01) = pd(03,01) + rrt(059) * density(02) 
  pd(03,02) = pd(03,02) + rrt(059) * density(01) 
  pd(03,01) = pd(03,01) - rrt(060) * density(03) 
  pd(03,03) = pd(03,03) - rrt(060) * density(01) 
  pd(04,01) = pd(04,01) + rrt(060) * density(03) 
  pd(04,03) = pd(04,03) + rrt(060) * density(01) 
  pd(04,01) = pd(04,01) - rrt(061) * density(04) 
  pd(04,04) = pd(04,04) - rrt(061) * density(01) 
  pd(05,01) = pd(05,01) + rrt(061) * density(04) 
  pd(05,04) = pd(05,04) + rrt(061) * density(01) 
  pd(05,01) = pd(05,01) - rrt(062) * density(05) 
  pd(05,05) = pd(05,05) - rrt(062) * density(01) 
  pd(06,01) = pd(06,01) + rrt(062) * density(05) 
  pd(06,05) = pd(06,05) + rrt(062) * density(01) 
  pd(06,01) = pd(06,01) - rrt(063) * density(06) 
  pd(06,06) = pd(06,06) - rrt(063) * density(01) 
  pd(07,01) = pd(07,01) + rrt(063) * density(06) 
  pd(07,06) = pd(07,06) + rrt(063) * density(01) 
  pd(07,01) = pd(07,01) - rrt(064) * density(07) 
  pd(07,07) = pd(07,07) - rrt(064) * density(01) 
  pd(08,01) = pd(08,01) + rrt(064) * density(07) 
  pd(08,07) = pd(08,07) + rrt(064) * density(01) 
  pd(08,01) = pd(08,01) - rrt(065) * density(08) 
  pd(08,08) = pd(08,08) - rrt(065) * density(01) 
  pd(09,01) = pd(09,01) + rrt(065) * density(08) 
  pd(09,08) = pd(09,08) + rrt(065) * density(01) 
  pd(01,02) = pd(01,02) + rrt(066) * density(23) 
  pd(01,23) = pd(01,23) + rrt(066) * density(02) 
  pd(02,02) = pd(02,02) - rrt(066) * density(23) 
  pd(02,23) = pd(02,23) - rrt(066) * density(02) 
  pd(02,03) = pd(02,03) + rrt(067) * density(23) 
  pd(02,23) = pd(02,23) + rrt(067) * density(03) 
  pd(03,03) = pd(03,03) - rrt(067) * density(23) 
  pd(03,23) = pd(03,23) - rrt(067) * density(03) 
  pd(03,04) = pd(03,04) + rrt(068) * density(23) 
  pd(03,23) = pd(03,23) + rrt(068) * density(04) 
  pd(04,04) = pd(04,04) - rrt(068) * density(23) 
  pd(04,23) = pd(04,23) - rrt(068) * density(04) 
  pd(04,05) = pd(04,05) + rrt(069) * density(23) 
  pd(04,23) = pd(04,23) + rrt(069) * density(05) 
  pd(05,05) = pd(05,05) - rrt(069) * density(23) 
  pd(05,23) = pd(05,23) - rrt(069) * density(05) 
  pd(05,06) = pd(05,06) + rrt(070) * density(23) 
  pd(05,23) = pd(05,23) + rrt(070) * density(06) 
  pd(06,06) = pd(06,06) - rrt(070) * density(23) 
  pd(06,23) = pd(06,23) - rrt(070) * density(06) 
  pd(06,07) = pd(06,07) + rrt(071) * density(23) 
  pd(06,23) = pd(06,23) + rrt(071) * density(07) 
  pd(07,07) = pd(07,07) - rrt(071) * density(23) 
  pd(07,23) = pd(07,23) - rrt(071) * density(07) 
  pd(07,08) = pd(07,08) + rrt(072) * density(23) 
  pd(07,23) = pd(07,23) + rrt(072) * density(08) 
  pd(08,08) = pd(08,08) - rrt(072) * density(23) 
  pd(08,23) = pd(08,23) - rrt(072) * density(08) 
  pd(08,09) = pd(08,09) + rrt(073) * density(23) 
  pd(08,23) = pd(08,23) + rrt(073) * density(09) 
  pd(09,09) = pd(09,09) - rrt(073) * density(23) 
  pd(09,23) = pd(09,23) - rrt(073) * density(09) 
  pd(01,01) = pd(01,01) - rrt(074) * density(23) 
  pd(01,23) = pd(01,23) - rrt(074) * density(01) 
  pd(02,01) = pd(02,01) + rrt(074) * density(23) 
  pd(02,23) = pd(02,23) + rrt(074) * density(01) 
  pd(02,02) = pd(02,02) - rrt(075) * density(23) 
  pd(02,23) = pd(02,23) - rrt(075) * density(02) 
  pd(03,02) = pd(03,02) + rrt(075) * density(23) 
  pd(03,23) = pd(03,23) + rrt(075) * density(02) 
  pd(03,03) = pd(03,03) - rrt(076) * density(23) 
  pd(03,23) = pd(03,23) - rrt(076) * density(03) 
  pd(04,03) = pd(04,03) + rrt(076) * density(23) 
  pd(04,23) = pd(04,23) + rrt(076) * density(03) 
  pd(04,04) = pd(04,04) - rrt(077) * density(23) 
  pd(04,23) = pd(04,23) - rrt(077) * density(04) 
  pd(05,04) = pd(05,04) + rrt(077) * density(23) 
  pd(05,23) = pd(05,23) + rrt(077) * density(04) 
  pd(05,05) = pd(05,05) - rrt(078) * density(23) 
  pd(05,23) = pd(05,23) - rrt(078) * density(05) 
  pd(06,05) = pd(06,05) + rrt(078) * density(23) 
  pd(06,23) = pd(06,23) + rrt(078) * density(05) 
  pd(06,06) = pd(06,06) - rrt(079) * density(23) 
  pd(06,23) = pd(06,23) - rrt(079) * density(06) 
  pd(07,06) = pd(07,06) + rrt(079) * density(23) 
  pd(07,23) = pd(07,23) + rrt(079) * density(06) 
  pd(07,07) = pd(07,07) - rrt(080) * density(23) 
  pd(07,23) = pd(07,23) - rrt(080) * density(07) 
  pd(08,07) = pd(08,07) + rrt(080) * density(23) 
  pd(08,23) = pd(08,23) + rrt(080) * density(07) 
  pd(08,08) = pd(08,08) - rrt(081) * density(23) 
  pd(08,23) = pd(08,23) - rrt(081) * density(08) 
  pd(09,08) = pd(09,08) + rrt(081) * density(23) 
  pd(09,23) = pd(09,23) + rrt(081) * density(08) 
  pd(01,02) = pd(01,02) + rrt(082) * density(38) 
  pd(01,38) = pd(01,38) + rrt(082) * density(02) 
  pd(02,02) = pd(02,02) - rrt(082) * density(38) 
  pd(02,38) = pd(02,38) - rrt(082) * density(02) 
  pd(02,03) = pd(02,03) + rrt(083) * density(38) 
  pd(02,38) = pd(02,38) + rrt(083) * density(03) 
  pd(03,03) = pd(03,03) - rrt(083) * density(38) 
  pd(03,38) = pd(03,38) - rrt(083) * density(03) 
  pd(03,04) = pd(03,04) + rrt(084) * density(38) 
  pd(03,38) = pd(03,38) + rrt(084) * density(04) 
  pd(04,04) = pd(04,04) - rrt(084) * density(38) 
  pd(04,38) = pd(04,38) - rrt(084) * density(04) 
  pd(04,05) = pd(04,05) + rrt(085) * density(38) 
  pd(04,38) = pd(04,38) + rrt(085) * density(05) 
  pd(05,05) = pd(05,05) - rrt(085) * density(38) 
  pd(05,38) = pd(05,38) - rrt(085) * density(05) 
  pd(05,06) = pd(05,06) + rrt(086) * density(38) 
  pd(05,38) = pd(05,38) + rrt(086) * density(06) 
  pd(06,06) = pd(06,06) - rrt(086) * density(38) 
  pd(06,38) = pd(06,38) - rrt(086) * density(06) 
  pd(06,07) = pd(06,07) + rrt(087) * density(38) 
  pd(06,38) = pd(06,38) + rrt(087) * density(07) 
  pd(07,07) = pd(07,07) - rrt(087) * density(38) 
  pd(07,38) = pd(07,38) - rrt(087) * density(07) 
  pd(07,08) = pd(07,08) + rrt(088) * density(38) 
  pd(07,38) = pd(07,38) + rrt(088) * density(08) 
  pd(08,08) = pd(08,08) - rrt(088) * density(38) 
  pd(08,38) = pd(08,38) - rrt(088) * density(08) 
  pd(08,09) = pd(08,09) + rrt(089) * density(38) 
  pd(08,38) = pd(08,38) + rrt(089) * density(09) 
  pd(09,09) = pd(09,09) - rrt(089) * density(38) 
  pd(09,38) = pd(09,38) - rrt(089) * density(09) 
  pd(01,01) = pd(01,01) - rrt(090) * density(38) 
  pd(01,38) = pd(01,38) - rrt(090) * density(01) 
  pd(02,01) = pd(02,01) + rrt(090) * density(38) 
  pd(02,38) = pd(02,38) + rrt(090) * density(01) 
  pd(02,02) = pd(02,02) - rrt(091) * density(38) 
  pd(02,38) = pd(02,38) - rrt(091) * density(02) 
  pd(03,02) = pd(03,02) + rrt(091) * density(38) 
  pd(03,38) = pd(03,38) + rrt(091) * density(02) 
  pd(03,03) = pd(03,03) - rrt(092) * density(38) 
  pd(03,38) = pd(03,38) - rrt(092) * density(03) 
  pd(04,03) = pd(04,03) + rrt(092) * density(38) 
  pd(04,38) = pd(04,38) + rrt(092) * density(03) 
  pd(04,04) = pd(04,04) - rrt(093) * density(38) 
  pd(04,38) = pd(04,38) - rrt(093) * density(04) 
  pd(05,04) = pd(05,04) + rrt(093) * density(38) 
  pd(05,38) = pd(05,38) + rrt(093) * density(04) 
  pd(05,05) = pd(05,05) - rrt(094) * density(38) 
  pd(05,38) = pd(05,38) - rrt(094) * density(05) 
  pd(06,05) = pd(06,05) + rrt(094) * density(38) 
  pd(06,38) = pd(06,38) + rrt(094) * density(05) 
  pd(06,06) = pd(06,06) - rrt(095) * density(38) 
  pd(06,38) = pd(06,38) - rrt(095) * density(06) 
  pd(07,06) = pd(07,06) + rrt(095) * density(38) 
  pd(07,38) = pd(07,38) + rrt(095) * density(06) 
  pd(07,07) = pd(07,07) - rrt(096) * density(38) 
  pd(07,38) = pd(07,38) - rrt(096) * density(07) 
  pd(08,07) = pd(08,07) + rrt(096) * density(38) 
  pd(08,38) = pd(08,38) + rrt(096) * density(07) 
  pd(08,08) = pd(08,08) - rrt(097) * density(38) 
  pd(08,38) = pd(08,38) - rrt(097) * density(08) 
  pd(09,08) = pd(09,08) + rrt(097) * density(38) 
  pd(09,38) = pd(09,38) + rrt(097) * density(08) 
  pd(30,30) = pd(30,30) + rrt(098) * density(31) 
  pd(30,31) = pd(30,31) + rrt(098) * density(30) 
  pd(31,30) = pd(31,30) - rrt(098) * density(31) 
  pd(31,31) = pd(31,31) - rrt(098) * density(30) 
  pd(31,30) = pd(31,30) + rrt(099) * density(32) 
  pd(31,32) = pd(31,32) + rrt(099) * density(30) 
  pd(32,30) = pd(32,30) - rrt(099) * density(32) 
  pd(32,32) = pd(32,32) - rrt(099) * density(30) 
  pd(32,30) = pd(32,30) + rrt(100) * density(33) 
  pd(32,33) = pd(32,33) + rrt(100) * density(30) 
  pd(33,30) = pd(33,30) - rrt(100) * density(33) 
  pd(33,33) = pd(33,33) - rrt(100) * density(30) 
  pd(33,30) = pd(33,30) + rrt(101) * density(34) 
  pd(33,34) = pd(33,34) + rrt(101) * density(30) 
  pd(34,30) = pd(34,30) - rrt(101) * density(34) 
  pd(34,34) = pd(34,34) - rrt(101) * density(30) 
  pd(30,30) = pd(30,30) - rrt(102) * density(30) * 2.0d0
  pd(31,30) = pd(31,30) + rrt(102) * density(30) * 2.0d0
  pd(31,30) = pd(31,30) - rrt(103) * density(31) 
  pd(31,31) = pd(31,31) - rrt(103) * density(30) 
  pd(32,30) = pd(32,30) + rrt(103) * density(31) 
  pd(32,31) = pd(32,31) + rrt(103) * density(30) 
  pd(32,30) = pd(32,30) - rrt(104) * density(32) 
  pd(32,32) = pd(32,32) - rrt(104) * density(30) 
  pd(33,30) = pd(33,30) + rrt(104) * density(32) 
  pd(33,32) = pd(33,32) + rrt(104) * density(30) 
  pd(33,30) = pd(33,30) - rrt(105) * density(33) 
  pd(33,33) = pd(33,33) - rrt(105) * density(30) 
  pd(34,30) = pd(34,30) + rrt(105) * density(33) 
  pd(34,33) = pd(34,33) + rrt(105) * density(30) 
  pd(30,31) = pd(30,31) + rrt(106) * density(38) 
  pd(30,38) = pd(30,38) + rrt(106) * density(31) 
  pd(31,31) = pd(31,31) - rrt(106) * density(38) 
  pd(31,38) = pd(31,38) - rrt(106) * density(31) 
  pd(31,32) = pd(31,32) + rrt(107) * density(38) 
  pd(31,38) = pd(31,38) + rrt(107) * density(32) 
  pd(32,32) = pd(32,32) - rrt(107) * density(38) 
  pd(32,38) = pd(32,38) - rrt(107) * density(32) 
  pd(32,33) = pd(32,33) + rrt(108) * density(38) 
  pd(32,38) = pd(32,38) + rrt(108) * density(33) 
  pd(33,33) = pd(33,33) - rrt(108) * density(38) 
  pd(33,38) = pd(33,38) - rrt(108) * density(33) 
  pd(33,34) = pd(33,34) + rrt(109) * density(38) 
  pd(33,38) = pd(33,38) + rrt(109) * density(34) 
  pd(34,34) = pd(34,34) - rrt(109) * density(38) 
  pd(34,38) = pd(34,38) - rrt(109) * density(34) 
  pd(30,30) = pd(30,30) - rrt(110) * density(38) 
  pd(30,38) = pd(30,38) - rrt(110) * density(30) 
  pd(31,30) = pd(31,30) + rrt(110) * density(38) 
  pd(31,38) = pd(31,38) + rrt(110) * density(30) 
  pd(31,31) = pd(31,31) - rrt(111) * density(38) 
  pd(31,38) = pd(31,38) - rrt(111) * density(31) 
  pd(32,31) = pd(32,31) + rrt(111) * density(38) 
  pd(32,38) = pd(32,38) + rrt(111) * density(31) 
  pd(32,32) = pd(32,32) - rrt(112) * density(38) 
  pd(32,38) = pd(32,38) - rrt(112) * density(32) 
  pd(33,32) = pd(33,32) + rrt(112) * density(38) 
  pd(33,38) = pd(33,38) + rrt(112) * density(32) 
  pd(33,33) = pd(33,33) - rrt(113) * density(38) 
  pd(33,38) = pd(33,38) - rrt(113) * density(33) 
  pd(34,33) = pd(34,33) + rrt(113) * density(38) 
  pd(34,38) = pd(34,38) + rrt(113) * density(33) 
  pd(01,01) = pd(01,01) - rrt(114) * density(79) 
  pd(01,79) = pd(01,79) - rrt(114) * density(01) 
  pd(17,01) = pd(17,01) + rrt(114) * density(79) 
  pd(17,79) = pd(17,79) + rrt(114) * density(01) 
  pd(01,01) = pd(01,01) - rrt(115) * density(79) 
  pd(01,79) = pd(01,79) - rrt(115) * density(01) 
  pd(17,01) = pd(17,01) + rrt(115) * density(79) 
  pd(17,79) = pd(17,79) + rrt(115) * density(01) 
  pd(01,01) = pd(01,01) - rrt(116) * density(79) 
  pd(01,79) = pd(01,79) - rrt(116) * density(01) 
  pd(17,01) = pd(17,01) + rrt(116) * density(79) 
  pd(17,79) = pd(17,79) + rrt(116) * density(01) 
  pd(01,01) = pd(01,01) - rrt(117) * density(79) 
  pd(01,79) = pd(01,79) - rrt(117) * density(01) 
  pd(18,01) = pd(18,01) + rrt(117) * density(79) 
  pd(18,79) = pd(18,79) + rrt(117) * density(01) 
  pd(01,01) = pd(01,01) - rrt(118) * density(79) 
  pd(01,79) = pd(01,79) - rrt(118) * density(01) 
  pd(18,01) = pd(18,01) + rrt(118) * density(79) 
  pd(18,79) = pd(18,79) + rrt(118) * density(01) 
  pd(01,01) = pd(01,01) - rrt(119) * density(79) 
  pd(01,79) = pd(01,79) - rrt(119) * density(01) 
  pd(18,01) = pd(18,01) + rrt(119) * density(79) 
  pd(18,79) = pd(18,79) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) - rrt(120) * density(79) 
  pd(01,79) = pd(01,79) - rrt(120) * density(01) 
  pd(18,01) = pd(18,01) + rrt(120) * density(79) 
  pd(18,79) = pd(18,79) + rrt(120) * density(01) 
  pd(01,01) = pd(01,01) - rrt(121) * density(79) 
  pd(01,79) = pd(01,79) - rrt(121) * density(01) 
  pd(18,01) = pd(18,01) + rrt(121) * density(79) 
  pd(18,79) = pd(18,79) + rrt(121) * density(01) 
  pd(01,01) = pd(01,01) - rrt(122) * density(79) 
  pd(01,79) = pd(01,79) - rrt(122) * density(01) 
  pd(18,01) = pd(18,01) + rrt(122) * density(79) 
  pd(18,79) = pd(18,79) + rrt(122) * density(01) 
  pd(01,01) = pd(01,01) - rrt(123) * density(79) 
  pd(01,79) = pd(01,79) - rrt(123) * density(01) 
  pd(18,01) = pd(18,01) + rrt(123) * density(79) 
  pd(18,79) = pd(18,79) + rrt(123) * density(01) 
  pd(01,01) = pd(01,01) - rrt(124) * density(79) 
  pd(01,79) = pd(01,79) - rrt(124) * density(01) 
  pd(19,01) = pd(19,01) + rrt(124) * density(79) 
  pd(19,79) = pd(19,79) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) - rrt(125) * density(79) 
  pd(01,79) = pd(01,79) - rrt(125) * density(01) 
  pd(19,01) = pd(19,01) + rrt(125) * density(79) 
  pd(19,79) = pd(19,79) + rrt(125) * density(01) 
  pd(01,01) = pd(01,01) - rrt(126) * density(79) 
  pd(01,79) = pd(01,79) - rrt(126) * density(01) 
  pd(19,01) = pd(19,01) + rrt(126) * density(79) 
  pd(19,79) = pd(19,79) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) - rrt(127) * density(79) 
  pd(01,79) = pd(01,79) - rrt(127) * density(01) 
  pd(19,01) = pd(19,01) + rrt(127) * density(79) 
  pd(19,79) = pd(19,79) + rrt(127) * density(01) 
  pd(01,01) = pd(01,01) - rrt(128) * density(79) 
  pd(01,79) = pd(01,79) - rrt(128) * density(01) 
  pd(19,01) = pd(19,01) + rrt(128) * density(79) 
  pd(19,79) = pd(19,79) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) - rrt(129) * density(79) 
  pd(01,79) = pd(01,79) - rrt(129) * density(01) 
  pd(19,01) = pd(19,01) + rrt(129) * density(79) 
  pd(19,79) = pd(19,79) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) - rrt(130) * density(79) 
  pd(01,79) = pd(01,79) - rrt(130) * density(01) 
  pd(20,01) = pd(20,01) + rrt(130) * density(79) 
  pd(20,79) = pd(20,79) + rrt(130) * density(01) 
  pd(01,01) = pd(01,01) - rrt(131) * density(79) 
  pd(01,79) = pd(01,79) - rrt(131) * density(01) 
  pd(20,01) = pd(20,01) + rrt(131) * density(79) 
  pd(20,79) = pd(20,79) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) - rrt(132) * density(79) 
  pd(01,79) = pd(01,79) - rrt(132) * density(01) 
  pd(20,01) = pd(20,01) + rrt(132) * density(79) 
  pd(20,79) = pd(20,79) + rrt(132) * density(01) 
  pd(01,01) = pd(01,01) - rrt(133) * density(79) 
  pd(01,79) = pd(01,79) - rrt(133) * density(01) 
  pd(21,01) = pd(21,01) + rrt(133) * density(79) 
  pd(21,79) = pd(21,79) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) - rrt(134) * density(79) 
  pd(01,79) = pd(01,79) - rrt(134) * density(01) 
  pd(21,01) = pd(21,01) + rrt(134) * density(79) 
  pd(21,79) = pd(21,79) + rrt(134) * density(01) 
  pd(01,01) = pd(01,01) - rrt(135) * density(79) 
  pd(01,79) = pd(01,79) - rrt(135) * density(01) 
  pd(21,01) = pd(21,01) + rrt(135) * density(79) 
  pd(21,79) = pd(21,79) + rrt(135) * density(01) 
  pd(01,01) = pd(01,01) - rrt(136) * density(79) 
  pd(01,79) = pd(01,79) - rrt(136) * density(01) 
  pd(21,01) = pd(21,01) + rrt(136) * density(79) 
  pd(21,79) = pd(21,79) + rrt(136) * density(01) 
  pd(01,01) = pd(01,01) - rrt(137) * density(79) 
  pd(01,79) = pd(01,79) - rrt(137) * density(01) 
  pd(21,01) = pd(21,01) + rrt(137) * density(79) 
  pd(21,79) = pd(21,79) + rrt(137) * density(01) 
  pd(01,01) = pd(01,01) - rrt(138) * density(79) 
  pd(01,79) = pd(01,79) - rrt(138) * density(01) 
  pd(21,01) = pd(21,01) + rrt(138) * density(79) 
  pd(21,79) = pd(21,79) + rrt(138) * density(01) 
  pd(01,01) = pd(01,01) - rrt(139) * density(79) 
  pd(01,79) = pd(01,79) - rrt(139) * density(01) 
  pd(22,01) = pd(22,01) + rrt(139) * density(79) 
  pd(22,79) = pd(22,79) + rrt(139) * density(01) 
  pd(01,01) = pd(01,01) - rrt(140) * density(79) 
  pd(01,79) = pd(01,79) - rrt(140) * density(01) 
  pd(22,01) = pd(22,01) + rrt(140) * density(79) 
  pd(22,79) = pd(22,79) + rrt(140) * density(01) 
  pd(01,01) = pd(01,01) - rrt(141) * density(79) 
  pd(01,79) = pd(01,79) - rrt(141) * density(01) 
  pd(22,01) = pd(22,01) + rrt(141) * density(79) 
  pd(22,79) = pd(22,79) + rrt(141) * density(01) 
  pd(01,01) = pd(01,01) - rrt(142) * density(79) 
  pd(01,79) = pd(01,79) - rrt(142) * density(01) 
  pd(23,01) = pd(23,01) + rrt(142) * density(79) 
  pd(23,79) = pd(23,79) + rrt(142) * density(01) 
  pd(24,01) = pd(24,01) + rrt(142) * density(79) 
  pd(24,79) = pd(24,79) + rrt(142) * density(01) 
  pd(30,30) = pd(30,30) - rrt(143) * density(79) 
  pd(30,79) = pd(30,79) - rrt(143) * density(30) 
  pd(35,30) = pd(35,30) + rrt(143) * density(79) 
  pd(35,79) = pd(35,79) + rrt(143) * density(30) 
  pd(30,30) = pd(30,30) - rrt(144) * density(79) 
  pd(30,79) = pd(30,79) - rrt(144) * density(30) 
  pd(36,30) = pd(36,30) + rrt(144) * density(79) 
  pd(36,79) = pd(36,79) + rrt(144) * density(30) 
  pd(30,30) = pd(30,30) - rrt(145) * density(79) 
  pd(30,79) = pd(30,79) - rrt(145) * density(30) 
  pd(37,30) = pd(37,30) + rrt(145) * density(79) 
  pd(37,79) = pd(37,79) + rrt(145) * density(30) 
  pd(30,30) = pd(30,30) - rrt(146) * density(79) 
  pd(30,79) = pd(30,79) - rrt(146) * density(30) 
  pd(38,30) = pd(38,30) + rrt(146) * density(79) * 2.0d0
  pd(38,79) = pd(38,79) + rrt(146) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(147) * density(79) 
  pd(30,79) = pd(30,79) - rrt(147) * density(30) 
  pd(38,30) = pd(38,30) + rrt(147) * density(79) 
  pd(38,79) = pd(38,79) + rrt(147) * density(30) 
  pd(39,30) = pd(39,30) + rrt(147) * density(79) 
  pd(39,79) = pd(39,79) + rrt(147) * density(30) 
  pd(30,30) = pd(30,30) - rrt(148) * density(79) 
  pd(30,79) = pd(30,79) - rrt(148) * density(30) 
  pd(38,30) = pd(38,30) + rrt(148) * density(79) 
  pd(38,79) = pd(38,79) + rrt(148) * density(30) 
  pd(40,30) = pd(40,30) + rrt(148) * density(79) 
  pd(40,79) = pd(40,79) + rrt(148) * density(30) 
  pd(35,35) = pd(35,35) - rrt(149) * density(79) 
  pd(35,79) = pd(35,79) - rrt(149) * density(35) 
  pd(38,35) = pd(38,35) + rrt(149) * density(79) * 2.0d0
  pd(38,79) = pd(38,79) + rrt(149) * density(35) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(150) * density(79) 
  pd(38,79) = pd(38,79) - rrt(150) * density(38) 
  pd(39,38) = pd(39,38) + rrt(150) * density(79) 
  pd(39,79) = pd(39,79) + rrt(150) * density(38) 
  pd(38,38) = pd(38,38) - rrt(151) * density(79) 
  pd(38,79) = pd(38,79) - rrt(151) * density(38) 
  pd(40,38) = pd(40,38) + rrt(151) * density(79) 
  pd(40,79) = pd(40,79) + rrt(151) * density(38) 
  pd(63,70) = pd(63,70) + rrt(152) * density(79) 
  pd(63,79) = pd(63,79) + rrt(152) * density(70) 
  pd(70,70) = pd(70,70) - rrt(152) * density(79) 
  pd(70,79) = pd(70,79) - rrt(152) * density(70) 
  pd(71,70) = pd(71,70) + rrt(152) * density(79) 
  pd(71,79) = pd(71,79) + rrt(152) * density(70) 
  pd(63,70) = pd(63,70) + rrt(153) * density(79) * 2.0d0
  pd(63,79) = pd(63,79) + rrt(153) * density(70) * 2.0d0
  pd(70,70) = pd(70,70) - rrt(153) * density(79) 
  pd(70,79) = pd(70,79) - rrt(153) * density(70) 
  pd(72,70) = pd(72,70) + rrt(153) * density(79) 
  pd(72,79) = pd(72,79) + rrt(153) * density(70) 
  pd(23,23) = pd(23,23) - rrt(154) * density(79) 
  pd(23,79) = pd(23,79) - rrt(154) * density(23) 
  pd(26,23) = pd(26,23) + rrt(154) * density(79) 
  pd(26,79) = pd(26,79) + rrt(154) * density(23) 
  pd(79,23) = pd(79,23) + rrt(154) * density(79) 
  pd(79,79) = pd(79,79) + rrt(154) * density(23) 
  pd(38,38) = pd(38,38) - rrt(155) * density(79) 
  pd(38,79) = pd(38,79) - rrt(155) * density(38) 
  pd(42,38) = pd(42,38) + rrt(155) * density(79) 
  pd(42,79) = pd(42,79) + rrt(155) * density(38) 
  pd(79,38) = pd(79,38) + rrt(155) * density(79) 
  pd(79,79) = pd(79,79) + rrt(155) * density(38) 
  pd(01,01) = pd(01,01) - rrt(156) * density(79) 
  pd(01,79) = pd(01,79) - rrt(156) * density(01) 
  pd(27,01) = pd(27,01) + rrt(156) * density(79) 
  pd(27,79) = pd(27,79) + rrt(156) * density(01) 
  pd(79,01) = pd(79,01) + rrt(156) * density(79) 
  pd(79,79) = pd(79,79) + rrt(156) * density(01) 
  pd(17,17) = pd(17,17) - rrt(157) * density(79) 
  pd(17,79) = pd(17,79) - rrt(157) * density(17) 
  pd(27,17) = pd(27,17) + rrt(157) * density(79) 
  pd(27,79) = pd(27,79) + rrt(157) * density(17) 
  pd(79,17) = pd(79,17) + rrt(157) * density(79) 
  pd(79,79) = pd(79,79) + rrt(157) * density(17) 
  pd(30,30) = pd(30,30) - rrt(158) * density(79) 
  pd(30,79) = pd(30,79) - rrt(158) * density(30) 
  pd(43,30) = pd(43,30) + rrt(158) * density(79) 
  pd(43,79) = pd(43,79) + rrt(158) * density(30) 
  pd(79,30) = pd(79,30) + rrt(158) * density(79) 
  pd(79,79) = pd(79,79) + rrt(158) * density(30) 
  pd(35,35) = pd(35,35) - rrt(159) * density(79) 
  pd(35,79) = pd(35,79) - rrt(159) * density(35) 
  pd(43,35) = pd(43,35) + rrt(159) * density(79) 
  pd(43,79) = pd(43,79) + rrt(159) * density(35) 
  pd(79,35) = pd(79,35) + rrt(159) * density(79) 
  pd(79,79) = pd(79,79) + rrt(159) * density(35) 
  pd(50,50) = pd(50,50) - rrt(160) * density(79) 
  pd(50,79) = pd(50,79) - rrt(160) * density(50) 
  pd(55,50) = pd(55,50) + rrt(160) * density(79) 
  pd(55,79) = pd(55,79) + rrt(160) * density(50) 
  pd(79,50) = pd(79,50) + rrt(160) * density(79) 
  pd(79,79) = pd(79,79) + rrt(160) * density(50) 
  pd(51,51) = pd(51,51) - rrt(161) * density(79) 
  pd(51,79) = pd(51,79) - rrt(161) * density(51) 
  pd(56,51) = pd(56,51) + rrt(161) * density(79) 
  pd(56,79) = pd(56,79) + rrt(161) * density(51) 
  pd(79,51) = pd(79,51) + rrt(161) * density(79) 
  pd(79,79) = pd(79,79) + rrt(161) * density(51) 
  pd(41,41) = pd(41,41) - rrt(162) * density(79) 
  pd(41,79) = pd(41,79) - rrt(162) * density(41) 
  pd(44,41) = pd(44,41) + rrt(162) * density(79) 
  pd(44,79) = pd(44,79) + rrt(162) * density(41) 
  pd(79,41) = pd(79,41) + rrt(162) * density(79) 
  pd(79,79) = pd(79,79) + rrt(162) * density(41) 
  pd(70,70) = pd(70,70) - rrt(163) * density(79) 
  pd(70,79) = pd(70,79) - rrt(163) * density(70) 
  pd(73,70) = pd(73,70) + rrt(163) * density(79) 
  pd(73,79) = pd(73,79) + rrt(163) * density(70) 
  pd(79,70) = pd(79,70) + rrt(163) * density(79) 
  pd(79,79) = pd(79,79) + rrt(163) * density(70) 
  pd(64,70) = pd(64,70) + rrt(164) * density(79) 
  pd(64,79) = pd(64,79) + rrt(164) * density(70) 
  pd(70,70) = pd(70,70) - rrt(164) * density(79) 
  pd(70,79) = pd(70,79) - rrt(164) * density(70) 
  pd(71,70) = pd(71,70) + rrt(164) * density(79) 
  pd(71,79) = pd(71,79) + rrt(164) * density(70) 
  pd(79,70) = pd(79,70) + rrt(164) * density(79) 
  pd(79,79) = pd(79,79) + rrt(164) * density(70) 
  pd(65,70) = pd(65,70) + rrt(165) * density(79) 
  pd(65,79) = pd(65,79) + rrt(165) * density(70) 
  pd(70,70) = pd(70,70) - rrt(165) * density(79) 
  pd(70,79) = pd(70,79) - rrt(165) * density(70) 
  pd(72,70) = pd(72,70) + rrt(165) * density(79) 
  pd(72,79) = pd(72,79) + rrt(165) * density(70) 
  pd(79,70) = pd(79,70) + rrt(165) * density(79) 
  pd(79,79) = pd(79,79) + rrt(165) * density(70) 
  pd(23,27) = pd(23,27) + rrt(166) * density(79) * 2.0d0
  pd(23,79) = pd(23,79) + rrt(166) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(166) * density(79) 
  pd(27,79) = pd(27,79) - rrt(166) * density(27) 
  pd(79,27) = pd(79,27) - rrt(166) * density(79) 
  pd(79,79) = pd(79,79) - rrt(166) * density(27) 
  pd(23,27) = pd(23,27) + rrt(167) * density(79) 
  pd(23,79) = pd(23,79) + rrt(167) * density(27) 
  pd(24,27) = pd(24,27) + rrt(167) * density(79) 
  pd(24,79) = pd(24,79) + rrt(167) * density(27) 
  pd(27,27) = pd(27,27) - rrt(167) * density(79) 
  pd(27,79) = pd(27,79) - rrt(167) * density(27) 
  pd(79,27) = pd(79,27) - rrt(167) * density(79) 
  pd(79,79) = pd(79,79) - rrt(167) * density(27) 
  pd(23,27) = pd(23,27) + rrt(168) * density(79) 
  pd(23,79) = pd(23,79) + rrt(168) * density(27) 
  pd(25,27) = pd(25,27) + rrt(168) * density(79) 
  pd(25,79) = pd(25,79) + rrt(168) * density(27) 
  pd(27,27) = pd(27,27) - rrt(168) * density(79) 
  pd(27,79) = pd(27,79) - rrt(168) * density(27) 
  pd(79,27) = pd(79,27) - rrt(168) * density(79) 
  pd(79,79) = pd(79,79) - rrt(168) * density(27) 
  pd(38,43) = pd(38,43) + rrt(169) * density(79) * 2.0d0
  pd(38,79) = pd(38,79) + rrt(169) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(169) * density(79) 
  pd(43,79) = pd(43,79) - rrt(169) * density(43) 
  pd(79,43) = pd(79,43) - rrt(169) * density(79) 
  pd(79,79) = pd(79,79) - rrt(169) * density(43) 
  pd(38,43) = pd(38,43) + rrt(170) * density(79) 
  pd(38,79) = pd(38,79) + rrt(170) * density(43) 
  pd(39,43) = pd(39,43) + rrt(170) * density(79) 
  pd(39,79) = pd(39,79) + rrt(170) * density(43) 
  pd(43,43) = pd(43,43) - rrt(170) * density(79) 
  pd(43,79) = pd(43,79) - rrt(170) * density(43) 
  pd(79,43) = pd(79,43) - rrt(170) * density(79) 
  pd(79,79) = pd(79,79) - rrt(170) * density(43) 
  pd(38,43) = pd(38,43) + rrt(171) * density(79) 
  pd(38,79) = pd(38,79) + rrt(171) * density(43) 
  pd(40,43) = pd(40,43) + rrt(171) * density(79) 
  pd(40,79) = pd(40,79) + rrt(171) * density(43) 
  pd(43,43) = pd(43,43) - rrt(171) * density(79) 
  pd(43,79) = pd(43,79) - rrt(171) * density(43) 
  pd(79,43) = pd(79,43) - rrt(171) * density(79) 
  pd(79,79) = pd(79,79) - rrt(171) * density(43) 
  pd(23,55) = pd(23,55) + rrt(172) * density(79) 
  pd(23,79) = pd(23,79) + rrt(172) * density(55) 
  pd(38,55) = pd(38,55) + rrt(172) * density(79) 
  pd(38,79) = pd(38,79) + rrt(172) * density(55) 
  pd(55,55) = pd(55,55) - rrt(172) * density(79) 
  pd(55,79) = pd(55,79) - rrt(172) * density(55) 
  pd(79,55) = pd(79,55) - rrt(172) * density(79) 
  pd(79,79) = pd(79,79) - rrt(172) * density(55) 
  pd(24,55) = pd(24,55) + rrt(173) * density(79) 
  pd(24,79) = pd(24,79) + rrt(173) * density(55) 
  pd(38,55) = pd(38,55) + rrt(173) * density(79) 
  pd(38,79) = pd(38,79) + rrt(173) * density(55) 
  pd(55,55) = pd(55,55) - rrt(173) * density(79) 
  pd(55,79) = pd(55,79) - rrt(173) * density(55) 
  pd(79,55) = pd(79,55) - rrt(173) * density(79) 
  pd(79,79) = pd(79,79) - rrt(173) * density(55) 
  pd(01,28) = pd(01,28) + rrt(174) * density(79) 
  pd(01,79) = pd(01,79) + rrt(174) * density(28) 
  pd(23,28) = pd(23,28) + rrt(174) * density(79) 
  pd(23,79) = pd(23,79) + rrt(174) * density(28) 
  pd(28,28) = pd(28,28) - rrt(174) * density(79) 
  pd(28,79) = pd(28,79) - rrt(174) * density(28) 
  pd(79,28) = pd(79,28) - rrt(174) * density(79) 
  pd(79,79) = pd(79,79) - rrt(174) * density(28) 
  pd(01,29) = pd(01,29) + rrt(175) * density(79) * 2.0d0
  pd(01,79) = pd(01,79) + rrt(175) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(175) * density(79) 
  pd(29,79) = pd(29,79) - rrt(175) * density(29) 
  pd(79,29) = pd(79,29) - rrt(175) * density(79) 
  pd(79,79) = pd(79,79) - rrt(175) * density(29) 
  pd(01,56) = pd(01,56) + rrt(176) * density(79) 
  pd(01,79) = pd(01,79) + rrt(176) * density(56) 
  pd(38,56) = pd(38,56) + rrt(176) * density(79) 
  pd(38,79) = pd(38,79) + rrt(176) * density(56) 
  pd(56,56) = pd(56,56) - rrt(176) * density(79) 
  pd(56,79) = pd(56,79) - rrt(176) * density(56) 
  pd(79,56) = pd(79,56) - rrt(176) * density(79) 
  pd(79,79) = pd(79,79) - rrt(176) * density(56) 
  pd(38,57) = pd(38,57) + rrt(177) * density(79) 
  pd(38,79) = pd(38,79) + rrt(177) * density(57) 
  pd(50,57) = pd(50,57) + rrt(177) * density(79) 
  pd(50,79) = pd(50,79) + rrt(177) * density(57) 
  pd(57,57) = pd(57,57) - rrt(177) * density(79) 
  pd(57,79) = pd(57,79) - rrt(177) * density(57) 
  pd(79,57) = pd(79,57) - rrt(177) * density(79) 
  pd(79,79) = pd(79,79) - rrt(177) * density(57) 
  pd(30,45) = pd(30,45) + rrt(178) * density(79) * 2.0d0
  pd(30,79) = pd(30,79) + rrt(178) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(178) * density(79) 
  pd(45,79) = pd(45,79) - rrt(178) * density(45) 
  pd(79,45) = pd(79,45) - rrt(178) * density(79) 
  pd(79,79) = pd(79,79) - rrt(178) * density(45) 
  pd(01,78) = pd(01,78) + rrt(179) * density(79) 
  pd(01,79) = pd(01,79) + rrt(179) * density(78) 
  pd(30,78) = pd(30,78) + rrt(179) * density(79) 
  pd(30,79) = pd(30,79) + rrt(179) * density(78) 
  pd(78,78) = pd(78,78) - rrt(179) * density(79) 
  pd(78,79) = pd(78,79) - rrt(179) * density(78) 
  pd(79,78) = pd(79,78) - rrt(179) * density(79) 
  pd(79,79) = pd(79,79) - rrt(179) * density(78) 
  pd(23,26) = pd(23,26) + rrt(180) * density(79)**2 
  pd(23,79) = pd(23,79) + rrt(180) * density(26) * density(79) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(180) * density(79)**2 
  pd(26,79) = pd(26,79) - rrt(180) * density(26) * density(79) * 2.0d0
  pd(79,26) = pd(79,26) - rrt(180) * density(79)**2 
  pd(79,79) = pd(79,79) - rrt(180) * density(26) * density(79) * 2.0d0
  pd(38,42) = pd(38,42) + rrt(181) * density(79)**2 
  pd(38,79) = pd(38,79) + rrt(181) * density(42) * density(79) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(181) * density(79)**2 
  pd(42,79) = pd(42,79) - rrt(181) * density(42) * density(79) * 2.0d0
  pd(79,42) = pd(79,42) - rrt(181) * density(79)**2 
  pd(79,79) = pd(79,79) - rrt(181) * density(42) * density(79) * 2.0d0
  pd(23,26) = pd(23,26) + rrt(182) * density(79) 
  pd(23,79) = pd(23,79) + rrt(182) * density(26) 
  pd(26,26) = pd(26,26) - rrt(182) * density(79) 
  pd(26,79) = pd(26,79) - rrt(182) * density(26) 
  pd(79,26) = pd(79,26) - rrt(182) * density(79) 
  pd(79,79) = pd(79,79) - rrt(182) * density(26) 
  pd(38,42) = pd(38,42) + rrt(183) * density(79) 
  pd(38,79) = pd(38,79) + rrt(183) * density(42) 
  pd(42,42) = pd(42,42) - rrt(183) * density(79) 
  pd(42,79) = pd(42,79) - rrt(183) * density(42) 
  pd(79,42) = pd(79,42) - rrt(183) * density(79) 
  pd(79,79) = pd(79,79) - rrt(183) * density(42) 
  pd(63,77) = pd(63,77) + rrt(184) * density(79) * 2.0d0
  pd(63,79) = pd(63,79) + rrt(184) * density(77) * 2.0d0
  pd(71,77) = pd(71,77) + rrt(184) * density(79) 
  pd(71,79) = pd(71,79) + rrt(184) * density(77) 
  pd(77,77) = pd(77,77) - rrt(184) * density(79) 
  pd(77,79) = pd(77,79) - rrt(184) * density(77) 
  pd(79,77) = pd(79,77) - rrt(184) * density(79) 
  pd(79,79) = pd(79,79) - rrt(184) * density(77) 
  pd(30,30) = pd(30,30) - rrt(185) * density(79) 
  pd(30,79) = pd(30,79) - rrt(185) * density(30) 
  pd(38,30) = pd(38,30) + rrt(185) * density(79) 
  pd(38,79) = pd(38,79) + rrt(185) * density(30) 
  pd(46,30) = pd(46,30) + rrt(185) * density(79) 
  pd(46,79) = pd(46,79) + rrt(185) * density(30) 
  pd(79,30) = pd(79,30) - rrt(185) * density(79) 
  pd(79,79) = pd(79,79) - rrt(185) * density(30) 
  pd(23,50) = pd(23,50) + rrt(186) * density(79) 
  pd(23,79) = pd(23,79) + rrt(186) * density(50) 
  pd(46,50) = pd(46,50) + rrt(186) * density(79) 
  pd(46,79) = pd(46,79) + rrt(186) * density(50) 
  pd(50,50) = pd(50,50) - rrt(186) * density(79) 
  pd(50,79) = pd(50,79) - rrt(186) * density(50) 
  pd(79,50) = pd(79,50) - rrt(186) * density(79) 
  pd(79,79) = pd(79,79) - rrt(186) * density(50) 
  pd(30,41) = pd(30,41) + rrt(187) * density(79) 
  pd(30,79) = pd(30,79) + rrt(187) * density(41) 
  pd(41,41) = pd(41,41) - rrt(187) * density(79) 
  pd(41,79) = pd(41,79) - rrt(187) * density(41) 
  pd(46,41) = pd(46,41) + rrt(187) * density(79) 
  pd(46,79) = pd(46,79) + rrt(187) * density(41) 
  pd(79,41) = pd(79,41) - rrt(187) * density(79) 
  pd(79,79) = pd(79,79) - rrt(187) * density(41) 
  pd(38,41) = pd(38,41) + rrt(188) * density(79) 
  pd(38,79) = pd(38,79) + rrt(188) * density(41) 
  pd(41,41) = pd(41,41) - rrt(188) * density(79) 
  pd(41,79) = pd(41,79) - rrt(188) * density(41) 
  pd(47,41) = pd(47,41) + rrt(188) * density(79) 
  pd(47,79) = pd(47,79) + rrt(188) * density(41) 
  pd(79,41) = pd(79,41) - rrt(188) * density(79) 
  pd(79,79) = pd(79,79) - rrt(188) * density(41) 
  pd(46,52) = pd(46,52) + rrt(189) * density(79) 
  pd(46,79) = pd(46,79) + rrt(189) * density(52) 
  pd(50,52) = pd(50,52) + rrt(189) * density(79) 
  pd(50,79) = pd(50,79) + rrt(189) * density(52) 
  pd(52,52) = pd(52,52) - rrt(189) * density(79) 
  pd(52,79) = pd(52,79) - rrt(189) * density(52) 
  pd(79,52) = pd(79,52) - rrt(189) * density(79) 
  pd(79,79) = pd(79,79) - rrt(189) * density(52) 
  pd(38,30) = pd(38,30) - rrt(190) * density(38) * density(79) 
  pd(38,38) = pd(38,38) - rrt(190) * density(30) * density(79) 
  pd(38,79) = pd(38,79) - rrt(190) * density(30) * density(38) 
  pd(46,30) = pd(46,30) + rrt(190) * density(38) * density(79) 
  pd(46,38) = pd(46,38) + rrt(190) * density(30) * density(79) 
  pd(46,79) = pd(46,79) + rrt(190) * density(30) * density(38) 
  pd(79,30) = pd(79,30) - rrt(190) * density(38) * density(79) 
  pd(79,38) = pd(79,38) - rrt(190) * density(30) * density(79) 
  pd(79,79) = pd(79,79) - rrt(190) * density(30) * density(38) 
  pd(30,30) = pd(30,30) - rrt(191) * density(38) * density(79) 
  pd(30,38) = pd(30,38) - rrt(191) * density(30) * density(79) 
  pd(30,79) = pd(30,79) - rrt(191) * density(30) * density(38) 
  pd(47,30) = pd(47,30) + rrt(191) * density(38) * density(79) 
  pd(47,38) = pd(47,38) + rrt(191) * density(30) * density(79) 
  pd(47,79) = pd(47,79) + rrt(191) * density(30) * density(38) 
  pd(79,30) = pd(79,30) - rrt(191) * density(38) * density(79) 
  pd(79,38) = pd(79,38) - rrt(191) * density(30) * density(79) 
  pd(79,79) = pd(79,79) - rrt(191) * density(30) * density(38) 
  pd(41,41) = pd(41,41) - rrt(192) * density(79) 
  pd(41,79) = pd(41,79) - rrt(192) * density(41) 
  pd(48,41) = pd(48,41) + rrt(192) * density(79) 
  pd(48,79) = pd(48,79) + rrt(192) * density(41) 
  pd(79,41) = pd(79,41) - rrt(192) * density(79) 
  pd(79,79) = pd(79,79) - rrt(192) * density(41) 
  pd(50,50) = pd(50,50) - rrt(193) * density(79) 
  pd(50,79) = pd(50,79) - rrt(193) * density(50) 
  pd(58,50) = pd(58,50) + rrt(193) * density(79) 
  pd(58,79) = pd(58,79) + rrt(193) * density(50) 
  pd(79,50) = pd(79,50) - rrt(193) * density(79) 
  pd(79,79) = pd(79,79) - rrt(193) * density(50) 
  pd(51,51) = pd(51,51) - rrt(194) * density(79) 
  pd(51,79) = pd(51,79) - rrt(194) * density(51) 
  pd(59,51) = pd(59,51) + rrt(194) * density(79) 
  pd(59,79) = pd(59,79) + rrt(194) * density(51) 
  pd(79,51) = pd(79,51) - rrt(194) * density(79) 
  pd(79,79) = pd(79,79) - rrt(194) * density(51) 
  pd(30,01) = pd(30,01) - rrt(195) * density(30) * density(79) 
  pd(30,30) = pd(30,30) - rrt(195) * density(01) * density(79) 
  pd(30,79) = pd(30,79) - rrt(195) * density(01) * density(30) 
  pd(47,01) = pd(47,01) + rrt(195) * density(30) * density(79) 
  pd(47,30) = pd(47,30) + rrt(195) * density(01) * density(79) 
  pd(47,79) = pd(47,79) + rrt(195) * density(01) * density(30) 
  pd(79,01) = pd(79,01) - rrt(195) * density(30) * density(79) 
  pd(79,30) = pd(79,30) - rrt(195) * density(01) * density(79) 
  pd(79,79) = pd(79,79) - rrt(195) * density(01) * density(30) 
  pd(67,70) = pd(67,70) + rrt(196) * density(79) 
  pd(67,79) = pd(67,79) + rrt(196) * density(70) 
  pd(70,70) = pd(70,70) - rrt(196) * density(79) 
  pd(70,79) = pd(70,79) - rrt(196) * density(70) 
  pd(71,70) = pd(71,70) + rrt(196) * density(79) 
  pd(71,79) = pd(71,79) + rrt(196) * density(70) 
  pd(79,70) = pd(79,70) - rrt(196) * density(79) 
  pd(79,79) = pd(79,79) - rrt(196) * density(70) 
  pd(30,38) = pd(30,38) + rrt(197) * density(46) 
  pd(30,46) = pd(30,46) + rrt(197) * density(38) 
  pd(38,38) = pd(38,38) - rrt(197) * density(46) 
  pd(38,46) = pd(38,46) - rrt(197) * density(38) 
  pd(46,38) = pd(46,38) - rrt(197) * density(46) 
  pd(46,46) = pd(46,46) - rrt(197) * density(38) 
  pd(79,38) = pd(79,38) + rrt(197) * density(46) 
  pd(79,46) = pd(79,46) + rrt(197) * density(38) 
  pd(23,23) = pd(23,23) - rrt(198) * density(46) 
  pd(23,46) = pd(23,46) - rrt(198) * density(23) 
  pd(46,23) = pd(46,23) - rrt(198) * density(46) 
  pd(46,46) = pd(46,46) - rrt(198) * density(23) 
  pd(50,23) = pd(50,23) + rrt(198) * density(46) 
  pd(50,46) = pd(50,46) + rrt(198) * density(23) 
  pd(79,23) = pd(79,23) + rrt(198) * density(46) 
  pd(79,46) = pd(79,46) + rrt(198) * density(23) 
  pd(46,46) = pd(46,46) - rrt(199) * density(50) 
  pd(46,50) = pd(46,50) - rrt(199) * density(46) 
  pd(50,46) = pd(50,46) - rrt(199) * density(50) 
  pd(50,50) = pd(50,50) - rrt(199) * density(46) 
  pd(52,46) = pd(52,46) + rrt(199) * density(50) 
  pd(52,50) = pd(52,50) + rrt(199) * density(46) 
  pd(79,46) = pd(79,46) + rrt(199) * density(50) 
  pd(79,50) = pd(79,50) + rrt(199) * density(46) 
  pd(01,01) = pd(01,01) - rrt(200) * density(46) 
  pd(01,46) = pd(01,46) - rrt(200) * density(01) 
  pd(46,01) = pd(46,01) - rrt(200) * density(46) 
  pd(46,46) = pd(46,46) - rrt(200) * density(01) 
  pd(51,01) = pd(51,01) + rrt(200) * density(46) 
  pd(51,46) = pd(51,46) + rrt(200) * density(01) 
  pd(79,01) = pd(79,01) + rrt(200) * density(46) 
  pd(79,46) = pd(79,46) + rrt(200) * density(01) 
  pd(30,30) = pd(30,30) - rrt(201) * density(46) 
  pd(30,46) = pd(30,46) - rrt(201) * density(30) 
  pd(41,30) = pd(41,30) + rrt(201) * density(46) 
  pd(41,46) = pd(41,46) + rrt(201) * density(30) 
  pd(46,30) = pd(46,30) - rrt(201) * density(46) 
  pd(46,46) = pd(46,46) - rrt(201) * density(30) 
  pd(79,30) = pd(79,30) + rrt(201) * density(46) 
  pd(79,46) = pd(79,46) + rrt(201) * density(30) 
  pd(35,35) = pd(35,35) - rrt(202) * density(46) 
  pd(35,46) = pd(35,46) - rrt(202) * density(35) 
  pd(41,35) = pd(41,35) + rrt(202) * density(46) 
  pd(41,46) = pd(41,46) + rrt(202) * density(35) 
  pd(46,35) = pd(46,35) - rrt(202) * density(46) 
  pd(46,46) = pd(46,46) - rrt(202) * density(35) 
  pd(79,35) = pd(79,35) + rrt(202) * density(46) 
  pd(79,46) = pd(79,46) + rrt(202) * density(35) 
  pd(30,36) = pd(30,36) + rrt(203) * density(46) 
  pd(30,46) = pd(30,46) + rrt(203) * density(36) 
  pd(36,36) = pd(36,36) - rrt(203) * density(46) 
  pd(36,46) = pd(36,46) - rrt(203) * density(36) 
  pd(38,36) = pd(38,36) + rrt(203) * density(46) 
  pd(38,46) = pd(38,46) + rrt(203) * density(36) 
  pd(46,36) = pd(46,36) - rrt(203) * density(46) 
  pd(46,46) = pd(46,46) - rrt(203) * density(36) 
  pd(79,36) = pd(79,36) + rrt(203) * density(46) 
  pd(79,46) = pd(79,46) + rrt(203) * density(36) 
  pd(01,17) = pd(01,17) + rrt(204) * density(46) 
  pd(01,46) = pd(01,46) + rrt(204) * density(17) 
  pd(17,17) = pd(17,17) - rrt(204) * density(46) 
  pd(17,46) = pd(17,46) - rrt(204) * density(17) 
  pd(38,17) = pd(38,17) + rrt(204) * density(46) 
  pd(38,46) = pd(38,46) + rrt(204) * density(17) 
  pd(46,17) = pd(46,17) - rrt(204) * density(46) 
  pd(46,46) = pd(46,46) - rrt(204) * density(17) 
  pd(79,17) = pd(79,17) + rrt(204) * density(46) 
  pd(79,46) = pd(79,46) + rrt(204) * density(17) 
  pd(01,18) = pd(01,18) + rrt(205) * density(46) 
  pd(01,46) = pd(01,46) + rrt(205) * density(18) 
  pd(18,18) = pd(18,18) - rrt(205) * density(46) 
  pd(18,46) = pd(18,46) - rrt(205) * density(18) 
  pd(38,18) = pd(38,18) + rrt(205) * density(46) 
  pd(38,46) = pd(38,46) + rrt(205) * density(18) 
  pd(46,18) = pd(46,18) - rrt(205) * density(46) 
  pd(46,46) = pd(46,46) - rrt(205) * density(18) 
  pd(79,18) = pd(79,18) + rrt(205) * density(46) 
  pd(79,46) = pd(79,46) + rrt(205) * density(18) 
  pd(30,41) = pd(30,41) + rrt(206) * density(46) * 2.0d0
  pd(30,46) = pd(30,46) + rrt(206) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(206) * density(46) 
  pd(41,46) = pd(41,46) - rrt(206) * density(41) 
  pd(46,41) = pd(46,41) - rrt(206) * density(46) 
  pd(46,46) = pd(46,46) - rrt(206) * density(41) 
  pd(79,41) = pd(79,41) + rrt(206) * density(46) 
  pd(79,46) = pd(79,46) + rrt(206) * density(41) 
  pd(38,38) = pd(38,38) - rrt(207) * density(47) 
  pd(38,47) = pd(38,47) - rrt(207) * density(38) 
  pd(41,38) = pd(41,38) + rrt(207) * density(47) 
  pd(41,47) = pd(41,47) + rrt(207) * density(38) 
  pd(47,38) = pd(47,38) - rrt(207) * density(47) 
  pd(47,47) = pd(47,47) - rrt(207) * density(38) 
  pd(79,38) = pd(79,38) + rrt(207) * density(47) 
  pd(79,47) = pd(79,47) + rrt(207) * density(38) 
  pd(23,23) = pd(23,23) - rrt(208) * density(47) 
  pd(23,47) = pd(23,47) - rrt(208) * density(23) 
  pd(47,23) = pd(47,23) - rrt(208) * density(47) 
  pd(47,47) = pd(47,47) - rrt(208) * density(23) 
  pd(52,23) = pd(52,23) + rrt(208) * density(47) 
  pd(52,47) = pd(52,47) + rrt(208) * density(23) 
  pd(79,23) = pd(79,23) + rrt(208) * density(47) 
  pd(79,47) = pd(79,47) + rrt(208) * density(23) 
  pd(30,30) = pd(30,30) + rrt(209) * density(47) 
  pd(30,47) = pd(30,47) + rrt(209) * density(30) 
  pd(47,30) = pd(47,30) - rrt(209) * density(47) 
  pd(47,47) = pd(47,47) - rrt(209) * density(30) 
  pd(79,30) = pd(79,30) + rrt(209) * density(47) 
  pd(79,47) = pd(79,47) + rrt(209) * density(30) 
  pd(30,35) = pd(30,35) + rrt(210) * density(47) * 2.0d0
  pd(30,47) = pd(30,47) + rrt(210) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(210) * density(47) 
  pd(35,47) = pd(35,47) - rrt(210) * density(35) 
  pd(47,35) = pd(47,35) - rrt(210) * density(47) 
  pd(47,47) = pd(47,47) - rrt(210) * density(35) 
  pd(79,35) = pd(79,35) + rrt(210) * density(47) 
  pd(79,47) = pd(79,47) + rrt(210) * density(35) 
  pd(30,36) = pd(30,36) + rrt(211) * density(47) * 2.0d0
  pd(30,47) = pd(30,47) + rrt(211) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(211) * density(47) 
  pd(36,47) = pd(36,47) - rrt(211) * density(36) 
  pd(47,36) = pd(47,36) - rrt(211) * density(47) 
  pd(47,47) = pd(47,47) - rrt(211) * density(36) 
  pd(79,36) = pd(79,36) + rrt(211) * density(47) 
  pd(79,47) = pd(79,47) + rrt(211) * density(36) 
  pd(30,01) = pd(30,01) + rrt(212) * density(47) 
  pd(30,47) = pd(30,47) + rrt(212) * density(01) 
  pd(47,01) = pd(47,01) - rrt(212) * density(47) 
  pd(47,47) = pd(47,47) - rrt(212) * density(01) 
  pd(79,01) = pd(79,01) + rrt(212) * density(47) 
  pd(79,47) = pd(79,47) + rrt(212) * density(01) 
  pd(01,17) = pd(01,17) + rrt(213) * density(47) 
  pd(01,47) = pd(01,47) + rrt(213) * density(17) 
  pd(17,17) = pd(17,17) - rrt(213) * density(47) 
  pd(17,47) = pd(17,47) - rrt(213) * density(17) 
  pd(30,17) = pd(30,17) + rrt(213) * density(47) 
  pd(30,47) = pd(30,47) + rrt(213) * density(17) 
  pd(47,17) = pd(47,17) - rrt(213) * density(47) 
  pd(47,47) = pd(47,47) - rrt(213) * density(17) 
  pd(79,17) = pd(79,17) + rrt(213) * density(47) 
  pd(79,47) = pd(79,47) + rrt(213) * density(17) 
  pd(01,18) = pd(01,18) + rrt(214) * density(47) 
  pd(01,47) = pd(01,47) + rrt(214) * density(18) 
  pd(18,18) = pd(18,18) - rrt(214) * density(47) 
  pd(18,47) = pd(18,47) - rrt(214) * density(18) 
  pd(30,18) = pd(30,18) + rrt(214) * density(47) 
  pd(30,47) = pd(30,47) + rrt(214) * density(18) 
  pd(47,18) = pd(47,18) - rrt(214) * density(47) 
  pd(47,47) = pd(47,47) - rrt(214) * density(18) 
  pd(79,18) = pd(79,18) + rrt(214) * density(47) 
  pd(79,47) = pd(79,47) + rrt(214) * density(18) 
  pd(30,38) = pd(30,38) + rrt(215) * density(48) * 2.0d0
  pd(30,48) = pd(30,48) + rrt(215) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(215) * density(48) 
  pd(38,48) = pd(38,48) - rrt(215) * density(38) 
  pd(48,38) = pd(48,38) - rrt(215) * density(48) 
  pd(48,48) = pd(48,48) - rrt(215) * density(38) 
  pd(79,38) = pd(79,38) + rrt(215) * density(48) 
  pd(79,48) = pd(79,48) + rrt(215) * density(38) 
  pd(23,23) = pd(23,23) - rrt(216) * density(58) 
  pd(23,58) = pd(23,58) - rrt(216) * density(23) 
  pd(51,23) = pd(51,23) + rrt(216) * density(58) 
  pd(51,58) = pd(51,58) + rrt(216) * density(23) 
  pd(58,23) = pd(58,23) - rrt(216) * density(58) 
  pd(58,58) = pd(58,58) - rrt(216) * density(23) 
  pd(79,23) = pd(79,23) + rrt(216) * density(58) 
  pd(79,58) = pd(79,58) + rrt(216) * density(23) 
  pd(23,23) = pd(23,23) - rrt(217) * density(48) 
  pd(23,48) = pd(23,48) - rrt(217) * density(23) 
  pd(30,23) = pd(30,23) + rrt(217) * density(48) 
  pd(30,48) = pd(30,48) + rrt(217) * density(23) 
  pd(48,23) = pd(48,23) - rrt(217) * density(48) 
  pd(48,48) = pd(48,48) - rrt(217) * density(23) 
  pd(50,23) = pd(50,23) + rrt(217) * density(48) 
  pd(50,48) = pd(50,48) + rrt(217) * density(23) 
  pd(79,23) = pd(79,23) + rrt(217) * density(48) 
  pd(79,48) = pd(79,48) + rrt(217) * density(23) 
  pd(01,23) = pd(01,23) + rrt(218) * density(59) 
  pd(01,59) = pd(01,59) + rrt(218) * density(23) 
  pd(23,23) = pd(23,23) - rrt(218) * density(59) 
  pd(23,59) = pd(23,59) - rrt(218) * density(23) 
  pd(50,23) = pd(50,23) + rrt(218) * density(59) 
  pd(50,59) = pd(50,59) + rrt(218) * density(23) 
  pd(59,23) = pd(59,23) - rrt(218) * density(59) 
  pd(59,59) = pd(59,59) - rrt(218) * density(23) 
  pd(79,23) = pd(79,23) + rrt(218) * density(59) 
  pd(79,59) = pd(79,59) + rrt(218) * density(23) 
  pd(23,23) = pd(23,23) - rrt(219) * density(60) 
  pd(23,60) = pd(23,60) - rrt(219) * density(23) 
  pd(50,23) = pd(50,23) + rrt(219) * density(60) * 2.0d0
  pd(50,60) = pd(50,60) + rrt(219) * density(23) * 2.0d0
  pd(60,23) = pd(60,23) - rrt(219) * density(60) 
  pd(60,60) = pd(60,60) - rrt(219) * density(23) 
  pd(79,23) = pd(79,23) + rrt(219) * density(60) 
  pd(79,60) = pd(79,60) + rrt(219) * density(23) 
  pd(23,23) = pd(23,23) - rrt(220) * density(61) 
  pd(23,61) = pd(23,61) - rrt(220) * density(23) 
  pd(50,23) = pd(50,23) + rrt(220) * density(61) 
  pd(50,61) = pd(50,61) + rrt(220) * density(23) 
  pd(52,23) = pd(52,23) + rrt(220) * density(61) 
  pd(52,61) = pd(52,61) + rrt(220) * density(23) 
  pd(61,23) = pd(61,23) - rrt(220) * density(61) 
  pd(61,61) = pd(61,61) - rrt(220) * density(23) 
  pd(79,23) = pd(79,23) + rrt(220) * density(61) 
  pd(79,61) = pd(79,61) + rrt(220) * density(23) 
  pd(38,38) = pd(38,38) - rrt(221) * density(58) 
  pd(38,58) = pd(38,58) - rrt(221) * density(38) 
  pd(52,38) = pd(52,38) + rrt(221) * density(58) 
  pd(52,58) = pd(52,58) + rrt(221) * density(38) 
  pd(58,38) = pd(58,38) - rrt(221) * density(58) 
  pd(58,58) = pd(58,58) - rrt(221) * density(38) 
  pd(79,38) = pd(79,38) + rrt(221) * density(58) 
  pd(79,58) = pd(79,58) + rrt(221) * density(38) 
  pd(38,38) = pd(38,38) - rrt(222) * density(59) 
  pd(38,59) = pd(38,59) - rrt(222) * density(38) 
  pd(50,38) = pd(50,38) + rrt(222) * density(59) * 2.0d0
  pd(50,59) = pd(50,59) + rrt(222) * density(38) * 2.0d0
  pd(59,38) = pd(59,38) - rrt(222) * density(59) 
  pd(59,59) = pd(59,59) - rrt(222) * density(38) 
  pd(79,38) = pd(79,38) + rrt(222) * density(59) 
  pd(79,59) = pd(79,59) + rrt(222) * density(38) 
  pd(30,38) = pd(30,38) + rrt(223) * density(60) 
  pd(30,60) = pd(30,60) + rrt(223) * density(38) 
  pd(38,38) = pd(38,38) - rrt(223) * density(60) 
  pd(38,60) = pd(38,60) - rrt(223) * density(38) 
  pd(50,38) = pd(50,38) + rrt(223) * density(60) 
  pd(50,60) = pd(50,60) + rrt(223) * density(38) 
  pd(60,38) = pd(60,38) - rrt(223) * density(60) 
  pd(60,60) = pd(60,60) - rrt(223) * density(38) 
  pd(79,38) = pd(79,38) + rrt(223) * density(60) 
  pd(79,60) = pd(79,60) + rrt(223) * density(38) 
  pd(38,38) = pd(38,38) - rrt(224) * density(61) 
  pd(38,61) = pd(38,61) - rrt(224) * density(38) 
  pd(41,38) = pd(41,38) + rrt(224) * density(61) 
  pd(41,61) = pd(41,61) + rrt(224) * density(38) 
  pd(50,38) = pd(50,38) + rrt(224) * density(61) 
  pd(50,61) = pd(50,61) + rrt(224) * density(38) 
  pd(61,38) = pd(61,38) - rrt(224) * density(61) 
  pd(61,61) = pd(61,61) - rrt(224) * density(38) 
  pd(79,38) = pd(79,38) + rrt(224) * density(61) 
  pd(79,61) = pd(79,61) + rrt(224) * density(38) 
  pd(01,17) = pd(01,17) + rrt(225) * density(48) 
  pd(01,48) = pd(01,48) + rrt(225) * density(17) 
  pd(17,17) = pd(17,17) - rrt(225) * density(48) 
  pd(17,48) = pd(17,48) - rrt(225) * density(17) 
  pd(41,17) = pd(41,17) + rrt(225) * density(48) 
  pd(41,48) = pd(41,48) + rrt(225) * density(17) 
  pd(48,17) = pd(48,17) - rrt(225) * density(48) 
  pd(48,48) = pd(48,48) - rrt(225) * density(17) 
  pd(79,17) = pd(79,17) + rrt(225) * density(48) 
  pd(79,48) = pd(79,48) + rrt(225) * density(17) 
  pd(01,17) = pd(01,17) + rrt(226) * density(58) 
  pd(01,58) = pd(01,58) + rrt(226) * density(17) 
  pd(17,17) = pd(17,17) - rrt(226) * density(58) 
  pd(17,58) = pd(17,58) - rrt(226) * density(17) 
  pd(50,17) = pd(50,17) + rrt(226) * density(58) 
  pd(50,58) = pd(50,58) + rrt(226) * density(17) 
  pd(58,17) = pd(58,17) - rrt(226) * density(58) 
  pd(58,58) = pd(58,58) - rrt(226) * density(17) 
  pd(79,17) = pd(79,17) + rrt(226) * density(58) 
  pd(79,58) = pd(79,58) + rrt(226) * density(17) 
  pd(01,17) = pd(01,17) + rrt(227) * density(59) 
  pd(01,59) = pd(01,59) + rrt(227) * density(17) 
  pd(17,17) = pd(17,17) - rrt(227) * density(59) 
  pd(17,59) = pd(17,59) - rrt(227) * density(17) 
  pd(51,17) = pd(51,17) + rrt(227) * density(59) 
  pd(51,59) = pd(51,59) + rrt(227) * density(17) 
  pd(59,17) = pd(59,17) - rrt(227) * density(59) 
  pd(59,59) = pd(59,59) - rrt(227) * density(17) 
  pd(79,17) = pd(79,17) + rrt(227) * density(59) 
  pd(79,59) = pd(79,59) + rrt(227) * density(17) 
  pd(01,17) = pd(01,17) + rrt(228) * density(60) 
  pd(01,60) = pd(01,60) + rrt(228) * density(17) 
  pd(17,17) = pd(17,17) - rrt(228) * density(60) 
  pd(17,60) = pd(17,60) - rrt(228) * density(17) 
  pd(52,17) = pd(52,17) + rrt(228) * density(60) 
  pd(52,60) = pd(52,60) + rrt(228) * density(17) 
  pd(60,17) = pd(60,17) - rrt(228) * density(60) 
  pd(60,60) = pd(60,60) - rrt(228) * density(17) 
  pd(79,17) = pd(79,17) + rrt(228) * density(60) 
  pd(79,60) = pd(79,60) + rrt(228) * density(17) 
  pd(01,17) = pd(01,17) + rrt(229) * density(61) 
  pd(01,61) = pd(01,61) + rrt(229) * density(17) 
  pd(17,17) = pd(17,17) - rrt(229) * density(61) 
  pd(17,61) = pd(17,61) - rrt(229) * density(17) 
  pd(53,17) = pd(53,17) + rrt(229) * density(61) 
  pd(53,61) = pd(53,61) + rrt(229) * density(17) 
  pd(61,17) = pd(61,17) - rrt(229) * density(61) 
  pd(61,61) = pd(61,61) - rrt(229) * density(17) 
  pd(79,17) = pd(79,17) + rrt(229) * density(61) 
  pd(79,61) = pd(79,61) + rrt(229) * density(17) 
  pd(01,18) = pd(01,18) + rrt(230) * density(48) 
  pd(01,48) = pd(01,48) + rrt(230) * density(18) 
  pd(18,18) = pd(18,18) - rrt(230) * density(48) 
  pd(18,48) = pd(18,48) - rrt(230) * density(18) 
  pd(41,18) = pd(41,18) + rrt(230) * density(48) 
  pd(41,48) = pd(41,48) + rrt(230) * density(18) 
  pd(48,18) = pd(48,18) - rrt(230) * density(48) 
  pd(48,48) = pd(48,48) - rrt(230) * density(18) 
  pd(79,18) = pd(79,18) + rrt(230) * density(48) 
  pd(79,48) = pd(79,48) + rrt(230) * density(18) 
  pd(01,18) = pd(01,18) + rrt(231) * density(58) 
  pd(01,58) = pd(01,58) + rrt(231) * density(18) 
  pd(18,18) = pd(18,18) - rrt(231) * density(58) 
  pd(18,58) = pd(18,58) - rrt(231) * density(18) 
  pd(50,18) = pd(50,18) + rrt(231) * density(58) 
  pd(50,58) = pd(50,58) + rrt(231) * density(18) 
  pd(58,18) = pd(58,18) - rrt(231) * density(58) 
  pd(58,58) = pd(58,58) - rrt(231) * density(18) 
  pd(79,18) = pd(79,18) + rrt(231) * density(58) 
  pd(79,58) = pd(79,58) + rrt(231) * density(18) 
  pd(01,18) = pd(01,18) + rrt(232) * density(59) 
  pd(01,59) = pd(01,59) + rrt(232) * density(18) 
  pd(18,18) = pd(18,18) - rrt(232) * density(59) 
  pd(18,59) = pd(18,59) - rrt(232) * density(18) 
  pd(51,18) = pd(51,18) + rrt(232) * density(59) 
  pd(51,59) = pd(51,59) + rrt(232) * density(18) 
  pd(59,18) = pd(59,18) - rrt(232) * density(59) 
  pd(59,59) = pd(59,59) - rrt(232) * density(18) 
  pd(79,18) = pd(79,18) + rrt(232) * density(59) 
  pd(79,59) = pd(79,59) + rrt(232) * density(18) 
  pd(01,18) = pd(01,18) + rrt(233) * density(60) 
  pd(01,60) = pd(01,60) + rrt(233) * density(18) 
  pd(18,18) = pd(18,18) - rrt(233) * density(60) 
  pd(18,60) = pd(18,60) - rrt(233) * density(18) 
  pd(52,18) = pd(52,18) + rrt(233) * density(60) 
  pd(52,60) = pd(52,60) + rrt(233) * density(18) 
  pd(60,18) = pd(60,18) - rrt(233) * density(60) 
  pd(60,60) = pd(60,60) - rrt(233) * density(18) 
  pd(79,18) = pd(79,18) + rrt(233) * density(60) 
  pd(79,60) = pd(79,60) + rrt(233) * density(18) 
  pd(01,18) = pd(01,18) + rrt(234) * density(61) 
  pd(01,61) = pd(01,61) + rrt(234) * density(18) 
  pd(18,18) = pd(18,18) - rrt(234) * density(61) 
  pd(18,61) = pd(18,61) - rrt(234) * density(18) 
  pd(53,18) = pd(53,18) + rrt(234) * density(61) 
  pd(53,61) = pd(53,61) + rrt(234) * density(18) 
  pd(61,18) = pd(61,18) - rrt(234) * density(61) 
  pd(61,61) = pd(61,61) - rrt(234) * density(18) 
  pd(79,18) = pd(79,18) + rrt(234) * density(61) 
  pd(79,61) = pd(79,61) + rrt(234) * density(18) 
  pd(01,17) = pd(01,17) + rrt(235) 
  pd(17,17) = pd(17,17) - rrt(235) 
  pd(17,18) = pd(17,18) + rrt(236) 
  pd(18,18) = pd(18,18) - rrt(236) 
  pd(01,19) = pd(01,19) + rrt(237) 
  pd(19,19) = pd(19,19) - rrt(237) 
  pd(18,20) = pd(18,20) + rrt(238) 
  pd(20,20) = pd(20,20) - rrt(238) 
  pd(30,35) = pd(30,35) + rrt(239) 
  pd(35,35) = pd(35,35) - rrt(239) 
  pd(35,36) = pd(35,36) + rrt(240) 
  pd(36,36) = pd(36,36) - rrt(240) 
  pd(30,36) = pd(30,36) + rrt(241) 
  pd(36,36) = pd(36,36) - rrt(241) 
  pd(30,37) = pd(30,37) + rrt(242) 
  pd(37,37) = pd(37,37) - rrt(242) 
  pd(01,17) = pd(01,17) + rrt(243) * density(30) 
  pd(01,30) = pd(01,30) + rrt(243) * density(17) 
  pd(17,17) = pd(17,17) - rrt(243) * density(30) 
  pd(17,30) = pd(17,30) - rrt(243) * density(17) 
  pd(30,17) = pd(30,17) - rrt(243) * density(30) 
  pd(30,30) = pd(30,30) - rrt(243) * density(17) 
  pd(38,17) = pd(38,17) + rrt(243) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(243) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) + rrt(244) * density(30) 
  pd(01,30) = pd(01,30) + rrt(244) * density(17) 
  pd(17,17) = pd(17,17) - rrt(244) * density(30) 
  pd(17,30) = pd(17,30) - rrt(244) * density(17) 
  pd(30,17) = pd(30,17) - rrt(244) * density(30) 
  pd(30,30) = pd(30,30) - rrt(244) * density(17) 
  pd(36,17) = pd(36,17) + rrt(244) * density(30) 
  pd(36,30) = pd(36,30) + rrt(244) * density(17) 
  pd(01,18) = pd(01,18) + rrt(245) * density(30) 
  pd(01,30) = pd(01,30) + rrt(245) * density(18) 
  pd(18,18) = pd(18,18) - rrt(245) * density(30) 
  pd(18,30) = pd(18,30) - rrt(245) * density(18) 
  pd(30,18) = pd(30,18) - rrt(245) * density(30) 
  pd(30,30) = pd(30,30) - rrt(245) * density(18) 
  pd(38,18) = pd(38,18) + rrt(245) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(245) * density(18) * 2.0d0
  pd(01,19) = pd(01,19) + rrt(246) * density(30) 
  pd(01,30) = pd(01,30) + rrt(246) * density(19) 
  pd(19,19) = pd(19,19) - rrt(246) * density(30) 
  pd(19,30) = pd(19,30) - rrt(246) * density(19) 
  pd(30,19) = pd(30,19) - rrt(246) * density(30) 
  pd(30,30) = pd(30,30) - rrt(246) * density(19) 
  pd(38,19) = pd(38,19) + rrt(246) * density(30) 
  pd(38,30) = pd(38,30) + rrt(246) * density(19) 
  pd(39,19) = pd(39,19) + rrt(246) * density(30) 
  pd(39,30) = pd(39,30) + rrt(246) * density(19) 
  pd(01,20) = pd(01,20) + rrt(247) * density(30) 
  pd(01,30) = pd(01,30) + rrt(247) * density(20) 
  pd(20,20) = pd(20,20) - rrt(247) * density(30) 
  pd(20,30) = pd(20,30) - rrt(247) * density(20) 
  pd(30,20) = pd(30,20) - rrt(247) * density(30) 
  pd(30,30) = pd(30,30) - rrt(247) * density(20) 
  pd(38,20) = pd(38,20) + rrt(247) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(247) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(248) * density(30) 
  pd(01,30) = pd(01,30) + rrt(248) * density(20) 
  pd(20,20) = pd(20,20) - rrt(248) * density(30) 
  pd(20,30) = pd(20,30) - rrt(248) * density(20) 
  pd(30,20) = pd(30,20) - rrt(248) * density(30) 
  pd(30,30) = pd(30,30) - rrt(248) * density(20) 
  pd(38,20) = pd(38,20) + rrt(248) * density(30) 
  pd(38,30) = pd(38,30) + rrt(248) * density(20) 
  pd(39,20) = pd(39,20) + rrt(248) * density(30) 
  pd(39,30) = pd(39,30) + rrt(248) * density(20) 
  pd(01,20) = pd(01,20) + rrt(249) * density(30) 
  pd(01,30) = pd(01,30) + rrt(249) * density(20) 
  pd(20,20) = pd(20,20) - rrt(249) * density(30) 
  pd(20,30) = pd(20,30) - rrt(249) * density(20) 
  pd(30,20) = pd(30,20) - rrt(249) * density(30) 
  pd(30,30) = pd(30,30) - rrt(249) * density(20) 
  pd(38,20) = pd(38,20) + rrt(249) * density(30) 
  pd(38,30) = pd(38,30) + rrt(249) * density(20) 
  pd(40,20) = pd(40,20) + rrt(249) * density(30) 
  pd(40,30) = pd(40,30) + rrt(249) * density(20) 
  pd(17,17) = pd(17,17) - rrt(250) * density(38) 
  pd(17,38) = pd(17,38) - rrt(250) * density(17) 
  pd(24,17) = pd(24,17) + rrt(250) * density(38) 
  pd(24,38) = pd(24,38) + rrt(250) * density(17) 
  pd(38,17) = pd(38,17) - rrt(250) * density(38) 
  pd(38,38) = pd(38,38) - rrt(250) * density(17) 
  pd(50,17) = pd(50,17) + rrt(250) * density(38) 
  pd(50,38) = pd(50,38) + rrt(250) * density(17) 
  pd(01,17) = pd(01,17) + rrt(251) * density(38) 
  pd(01,38) = pd(01,38) + rrt(251) * density(17) 
  pd(17,17) = pd(17,17) - rrt(251) * density(38) 
  pd(17,38) = pd(17,38) - rrt(251) * density(17) 
  pd(38,17) = pd(38,17) - rrt(251) * density(38) 
  pd(38,38) = pd(38,38) - rrt(251) * density(17) 
  pd(40,17) = pd(40,17) + rrt(251) * density(38) 
  pd(40,38) = pd(40,38) + rrt(251) * density(17) 
  pd(01,17) = pd(01,17) + rrt(252) * density(23) 
  pd(01,23) = pd(01,23) + rrt(252) * density(17) 
  pd(17,17) = pd(17,17) - rrt(252) * density(23) 
  pd(17,23) = pd(17,23) - rrt(252) * density(17) 
  pd(01,17) = pd(01,17) + rrt(253) * density(23) 
  pd(01,23) = pd(01,23) + rrt(253) * density(17) 
  pd(17,17) = pd(17,17) - rrt(253) * density(23) 
  pd(17,23) = pd(17,23) - rrt(253) * density(17) 
  pd(23,17) = pd(23,17) - rrt(253) * density(23) 
  pd(23,23) = pd(23,23) - rrt(253) * density(17) 
  pd(25,17) = pd(25,17) + rrt(253) * density(23) 
  pd(25,23) = pd(25,23) + rrt(253) * density(17) 
  pd(01,01) = pd(01,01) + rrt(254) * density(17) 
  pd(01,17) = pd(01,17) + rrt(254) * density(01) 
  pd(17,01) = pd(17,01) - rrt(254) * density(17) 
  pd(17,17) = pd(17,17) - rrt(254) * density(01) 
  pd(01,17) = pd(01,17) + rrt(255) * density(50) 
  pd(01,50) = pd(01,50) + rrt(255) * density(17) 
  pd(17,17) = pd(17,17) - rrt(255) * density(50) 
  pd(17,50) = pd(17,50) - rrt(255) * density(17) 
  pd(01,17) = pd(01,17) + rrt(256) * density(51) 
  pd(01,51) = pd(01,51) + rrt(256) * density(17) 
  pd(17,17) = pd(17,17) - rrt(256) * density(51) 
  pd(17,51) = pd(17,51) - rrt(256) * density(17) 
  pd(23,17) = pd(23,17) + rrt(256) * density(51) 
  pd(23,51) = pd(23,51) + rrt(256) * density(17) 
  pd(50,17) = pd(50,17) + rrt(256) * density(51) 
  pd(50,51) = pd(50,51) + rrt(256) * density(17) 
  pd(51,17) = pd(51,17) - rrt(256) * density(51) 
  pd(51,51) = pd(51,51) - rrt(256) * density(17) 
  pd(01,17) = pd(01,17) + rrt(257) * density(52) 
  pd(01,52) = pd(01,52) + rrt(257) * density(17) 
  pd(17,17) = pd(17,17) - rrt(257) * density(52) 
  pd(17,52) = pd(17,52) - rrt(257) * density(17) 
  pd(38,17) = pd(38,17) + rrt(257) * density(52) 
  pd(38,52) = pd(38,52) + rrt(257) * density(17) 
  pd(50,17) = pd(50,17) + rrt(257) * density(52) 
  pd(50,52) = pd(50,52) + rrt(257) * density(17) 
  pd(52,17) = pd(52,17) - rrt(257) * density(52) 
  pd(52,52) = pd(52,52) - rrt(257) * density(17) 
  pd(01,17) = pd(01,17) + rrt(258) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(258) * density(17) * 4.0d0
  pd(18,17) = pd(18,17) + rrt(258) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) + rrt(259) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(259) * density(17) * 4.0d0
  pd(20,17) = pd(20,17) + rrt(259) * density(17) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(260) * density(18) 
  pd(17,18) = pd(17,18) + rrt(260) * density(01) 
  pd(18,01) = pd(18,01) - rrt(260) * density(18) 
  pd(18,18) = pd(18,18) - rrt(260) * density(01) 
  pd(01,01) = pd(01,01) + rrt(261) * density(18) 
  pd(01,18) = pd(01,18) + rrt(261) * density(01) 
  pd(18,01) = pd(18,01) - rrt(261) * density(18) 
  pd(18,18) = pd(18,18) - rrt(261) * density(01) 
  pd(17,18) = pd(17,18) + rrt(262) * density(50) 
  pd(17,50) = pd(17,50) + rrt(262) * density(18) 
  pd(18,18) = pd(18,18) - rrt(262) * density(50) 
  pd(18,50) = pd(18,50) - rrt(262) * density(18) 
  pd(19,01) = pd(19,01) + rrt(263) * density(20) 
  pd(19,20) = pd(19,20) + rrt(263) * density(01) 
  pd(20,01) = pd(20,01) - rrt(263) * density(20) 
  pd(20,20) = pd(20,20) - rrt(263) * density(01) 
  pd(18,01) = pd(18,01) + rrt(264) * density(19) 
  pd(18,19) = pd(18,19) + rrt(264) * density(01) 
  pd(19,01) = pd(19,01) - rrt(264) * density(19) 
  pd(19,19) = pd(19,19) - rrt(264) * density(01) 
  pd(01,19) = pd(01,19) + rrt(265) * density(50) 
  pd(01,50) = pd(01,50) + rrt(265) * density(19) 
  pd(19,19) = pd(19,19) - rrt(265) * density(50) 
  pd(19,50) = pd(19,50) - rrt(265) * density(19) 
  pd(23,19) = pd(23,19) + rrt(265) * density(50) 
  pd(23,50) = pd(23,50) + rrt(265) * density(19) 
  pd(38,19) = pd(38,19) + rrt(265) * density(50) 
  pd(38,50) = pd(38,50) + rrt(265) * density(19) 
  pd(50,19) = pd(50,19) - rrt(265) * density(50) 
  pd(50,50) = pd(50,50) - rrt(265) * density(19) 
  pd(17,17) = pd(17,17) - rrt(266) * density(19) 
  pd(17,19) = pd(17,19) - rrt(266) * density(17) 
  pd(19,17) = pd(19,17) - rrt(266) * density(19) 
  pd(19,19) = pd(19,19) - rrt(266) * density(17) 
  pd(29,17) = pd(29,17) + rrt(266) * density(19) 
  pd(29,19) = pd(29,19) + rrt(266) * density(17) 
  pd(79,17) = pd(79,17) + rrt(266) * density(19) 
  pd(79,19) = pd(79,19) + rrt(266) * density(17) 
  pd(19,19) = pd(19,19) - rrt(267) * density(19) * 4.0d0
  pd(29,19) = pd(29,19) + rrt(267) * density(19) * 2.0d0
  pd(79,19) = pd(79,19) + rrt(267) * density(19) * 2.0d0
  pd(01,17) = pd(01,17) + rrt(268) * density(63) 
  pd(01,63) = pd(01,63) + rrt(268) * density(17) 
  pd(17,17) = pd(17,17) - rrt(268) * density(63) 
  pd(17,63) = pd(17,63) - rrt(268) * density(17) 
  pd(01,17) = pd(01,17) + rrt(269) * density(62) 
  pd(01,62) = pd(01,62) + rrt(269) * density(17) 
  pd(17,17) = pd(17,17) - rrt(269) * density(62) 
  pd(17,62) = pd(17,62) - rrt(269) * density(17) 
  pd(62,17) = pd(62,17) - rrt(269) * density(62) 
  pd(62,62) = pd(62,62) - rrt(269) * density(17) 
  pd(63,17) = pd(63,17) + rrt(269) * density(62) * 2.0d0
  pd(63,62) = pd(63,62) + rrt(269) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) + rrt(270) * density(70) 
  pd(01,70) = pd(01,70) + rrt(270) * density(17) 
  pd(17,17) = pd(17,17) - rrt(270) * density(70) 
  pd(17,70) = pd(17,70) - rrt(270) * density(17) 
  pd(17,18) = pd(17,18) + rrt(271) * density(62) 
  pd(17,62) = pd(17,62) + rrt(271) * density(18) 
  pd(18,18) = pd(18,18) - rrt(271) * density(62) 
  pd(18,62) = pd(18,62) - rrt(271) * density(18) 
  pd(01,19) = pd(01,19) + rrt(272) * density(63) 
  pd(01,63) = pd(01,63) + rrt(272) * density(19) 
  pd(19,19) = pd(19,19) - rrt(272) * density(63) 
  pd(19,63) = pd(19,63) - rrt(272) * density(19) 
  pd(01,19) = pd(01,19) + rrt(273) * density(62) 
  pd(01,62) = pd(01,62) + rrt(273) * density(19) 
  pd(19,19) = pd(19,19) - rrt(273) * density(62) 
  pd(19,62) = pd(19,62) - rrt(273) * density(19) 
  pd(17,01) = pd(17,01) + rrt(274) * density(23)**2 
  pd(17,23) = pd(17,23) + rrt(274) * density(01) * density(23) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(274) * density(23)**2 * 2.0d0
  pd(23,23) = pd(23,23) - rrt(274) * density(01) * density(23) * 4.0d0
  pd(17,23) = pd(17,23) + rrt(275) * density(23) * density(62) * 2.0d0
  pd(17,62) = pd(17,62) + rrt(275) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(275) * density(23) * density(62) * 4.0d0
  pd(23,62) = pd(23,62) - rrt(275) * density(23)**2 * 2.0d0
  pd(17,23) = pd(17,23) + rrt(276) * density(23) * density(70) * 2.0d0
  pd(17,70) = pd(17,70) + rrt(276) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(276) * density(23) * density(70) * 4.0d0
  pd(23,70) = pd(23,70) - rrt(276) * density(23)**2 * 2.0d0
  pd(17,23) = pd(17,23) + rrt(277) * density(23)**2 * 3.0d0
  pd(23,23) = pd(23,23) - rrt(277) * density(23)**2 * 6.0d0
  pd(17,23) = pd(17,23) + rrt(278) * density(23) * density(63) * 2.0d0
  pd(17,63) = pd(17,63) + rrt(278) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(278) * density(23) * density(63) * 4.0d0
  pd(23,63) = pd(23,63) - rrt(278) * density(23)**2 * 2.0d0
  pd(18,01) = pd(18,01) + rrt(279) * density(23)**2 
  pd(18,23) = pd(18,23) + rrt(279) * density(01) * density(23) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(279) * density(23)**2 * 2.0d0
  pd(23,23) = pd(23,23) - rrt(279) * density(01) * density(23) * 4.0d0
  pd(18,23) = pd(18,23) + rrt(280) * density(23) * density(62) * 2.0d0
  pd(18,62) = pd(18,62) + rrt(280) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(280) * density(23) * density(62) * 4.0d0
  pd(23,62) = pd(23,62) - rrt(280) * density(23)**2 * 2.0d0
  pd(18,23) = pd(18,23) + rrt(281) * density(23) * density(70) * 2.0d0
  pd(18,70) = pd(18,70) + rrt(281) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(281) * density(23) * density(70) * 4.0d0
  pd(23,70) = pd(23,70) - rrt(281) * density(23)**2 * 2.0d0
  pd(18,23) = pd(18,23) + rrt(282) * density(23)**2 * 3.0d0
  pd(23,23) = pd(23,23) - rrt(282) * density(23)**2 * 6.0d0
  pd(18,23) = pd(18,23) + rrt(283) * density(23) * density(63) * 2.0d0
  pd(18,63) = pd(18,63) + rrt(283) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(283) * density(23) * density(63) * 4.0d0
  pd(23,63) = pd(23,63) - rrt(283) * density(23)**2 * 2.0d0
  pd(23,24) = pd(23,24) + rrt(284) * density(38) 
  pd(23,38) = pd(23,38) + rrt(284) * density(24) 
  pd(24,24) = pd(24,24) - rrt(284) * density(38) 
  pd(24,38) = pd(24,38) - rrt(284) * density(24) 
  pd(38,24) = pd(38,24) - rrt(284) * density(38) 
  pd(38,38) = pd(38,38) - rrt(284) * density(24) 
  pd(39,24) = pd(39,24) + rrt(284) * density(38) 
  pd(39,38) = pd(39,38) + rrt(284) * density(24) 
  pd(24,24) = pd(24,24) - rrt(285) * density(30) 
  pd(24,30) = pd(24,30) - rrt(285) * density(24) 
  pd(30,24) = pd(30,24) - rrt(285) * density(30) 
  pd(30,30) = pd(30,30) - rrt(285) * density(24) 
  pd(38,24) = pd(38,24) + rrt(285) * density(30) 
  pd(38,30) = pd(38,30) + rrt(285) * density(24) 
  pd(50,24) = pd(50,24) + rrt(285) * density(30) 
  pd(50,30) = pd(50,30) + rrt(285) * density(24) 
  pd(01,24) = pd(01,24) + rrt(286) * density(50) 
  pd(01,50) = pd(01,50) + rrt(286) * density(24) 
  pd(24,24) = pd(24,24) - rrt(286) * density(50) 
  pd(24,50) = pd(24,50) - rrt(286) * density(24) 
  pd(38,24) = pd(38,24) + rrt(286) * density(50) 
  pd(38,50) = pd(38,50) + rrt(286) * density(24) 
  pd(50,24) = pd(50,24) - rrt(286) * density(50) 
  pd(50,50) = pd(50,50) - rrt(286) * density(24) 
  pd(01,24) = pd(01,24) + rrt(287) * density(51) 
  pd(01,51) = pd(01,51) + rrt(287) * density(24) 
  pd(24,24) = pd(24,24) - rrt(287) * density(51) 
  pd(24,51) = pd(24,51) - rrt(287) * density(24) 
  pd(50,24) = pd(50,24) + rrt(287) * density(51) 
  pd(50,51) = pd(50,51) + rrt(287) * density(24) 
  pd(51,24) = pd(51,24) - rrt(287) * density(51) 
  pd(51,51) = pd(51,51) - rrt(287) * density(24) 
  pd(23,01) = pd(23,01) + rrt(288) * density(24) 
  pd(23,24) = pd(23,24) + rrt(288) * density(01) 
  pd(24,01) = pd(24,01) - rrt(288) * density(24) 
  pd(24,24) = pd(24,24) - rrt(288) * density(01) 
  pd(23,23) = pd(23,23) + rrt(289) * density(25) 
  pd(23,25) = pd(23,25) + rrt(289) * density(23) 
  pd(25,23) = pd(25,23) - rrt(289) * density(25) 
  pd(25,25) = pd(25,25) - rrt(289) * density(23) 
  pd(23,25) = pd(23,25) + rrt(290) * density(38) 
  pd(23,38) = pd(23,38) + rrt(290) * density(25) 
  pd(25,25) = pd(25,25) - rrt(290) * density(38) 
  pd(25,38) = pd(25,38) - rrt(290) * density(25) 
  pd(24,23) = pd(24,23) + rrt(291) * density(25) 
  pd(24,25) = pd(24,25) + rrt(291) * density(23) 
  pd(25,23) = pd(25,23) - rrt(291) * density(25) 
  pd(25,25) = pd(25,25) - rrt(291) * density(23) 
  pd(23,01) = pd(23,01) + rrt(292) * density(25) 
  pd(23,25) = pd(23,25) + rrt(292) * density(01) 
  pd(25,01) = pd(25,01) - rrt(292) * density(25) 
  pd(25,25) = pd(25,25) - rrt(292) * density(01) 
  pd(24,24) = pd(24,24) - rrt(293) * density(25) 
  pd(24,25) = pd(24,25) - rrt(293) * density(24) 
  pd(25,24) = pd(25,24) - rrt(293) * density(25) 
  pd(25,25) = pd(25,25) - rrt(293) * density(24) 
  pd(27,24) = pd(27,24) + rrt(293) * density(25) 
  pd(27,25) = pd(27,25) + rrt(293) * density(24) 
  pd(79,24) = pd(79,24) + rrt(293) * density(25) 
  pd(79,25) = pd(79,25) + rrt(293) * density(24) 
  pd(25,25) = pd(25,25) - rrt(294) * density(30) 
  pd(25,30) = pd(25,30) - rrt(294) * density(25) 
  pd(30,25) = pd(30,25) - rrt(294) * density(30) 
  pd(30,30) = pd(30,30) - rrt(294) * density(25) 
  pd(38,25) = pd(38,25) + rrt(294) * density(30) 
  pd(38,30) = pd(38,30) + rrt(294) * density(25) 
  pd(50,25) = pd(50,25) + rrt(294) * density(30) 
  pd(50,30) = pd(50,30) + rrt(294) * density(25) 
  pd(17,25) = pd(17,25) + rrt(295) * density(50) 
  pd(17,50) = pd(17,50) + rrt(295) * density(25) 
  pd(25,25) = pd(25,25) - rrt(295) * density(50) 
  pd(25,50) = pd(25,50) - rrt(295) * density(25) 
  pd(38,25) = pd(38,25) + rrt(295) * density(50) 
  pd(38,50) = pd(38,50) + rrt(295) * density(25) 
  pd(50,25) = pd(50,25) - rrt(295) * density(50) 
  pd(50,50) = pd(50,50) - rrt(295) * density(25) 
  pd(24,24) = pd(24,24) - rrt(296) * density(62) 
  pd(24,62) = pd(24,62) - rrt(296) * density(24) 
  pd(62,24) = pd(62,24) - rrt(296) * density(62) 
  pd(62,62) = pd(62,62) - rrt(296) * density(24) 
  pd(63,24) = pd(63,24) + rrt(296) * density(62) 
  pd(63,62) = pd(63,62) + rrt(296) * density(24) 
  pd(72,24) = pd(72,24) + rrt(296) * density(62) 
  pd(72,62) = pd(72,62) + rrt(296) * density(24) 
  pd(24,24) = pd(24,24) - rrt(297) * density(70) 
  pd(24,70) = pd(24,70) - rrt(297) * density(24) 
  pd(70,24) = pd(70,24) - rrt(297) * density(70) 
  pd(70,70) = pd(70,70) - rrt(297) * density(24) 
  pd(71,24) = pd(71,24) + rrt(297) * density(70) 
  pd(71,70) = pd(71,70) + rrt(297) * density(24) 
  pd(72,24) = pd(72,24) + rrt(297) * density(70) 
  pd(72,70) = pd(72,70) + rrt(297) * density(24) 
  pd(25,25) = pd(25,25) - rrt(298) * density(62) 
  pd(25,62) = pd(25,62) - rrt(298) * density(25) 
  pd(62,25) = pd(62,25) - rrt(298) * density(62) 
  pd(62,62) = pd(62,62) - rrt(298) * density(25) 
  pd(63,25) = pd(63,25) + rrt(298) * density(62) 
  pd(63,62) = pd(63,62) + rrt(298) * density(25) 
  pd(72,25) = pd(72,25) + rrt(298) * density(62) 
  pd(72,62) = pd(72,62) + rrt(298) * density(25) 
  pd(30,35) = pd(30,35) + rrt(299) * density(38) 
  pd(30,38) = pd(30,38) + rrt(299) * density(35) 
  pd(35,35) = pd(35,35) - rrt(299) * density(38) 
  pd(35,38) = pd(35,38) - rrt(299) * density(35) 
  pd(23,23) = pd(23,23) - rrt(300) * density(35) 
  pd(23,35) = pd(23,35) - rrt(300) * density(23) 
  pd(35,23) = pd(35,23) - rrt(300) * density(35) 
  pd(35,35) = pd(35,35) - rrt(300) * density(23) 
  pd(38,23) = pd(38,23) + rrt(300) * density(35) 
  pd(38,35) = pd(38,35) + rrt(300) * density(23) 
  pd(50,23) = pd(50,23) + rrt(300) * density(35) 
  pd(50,35) = pd(50,35) + rrt(300) * density(23) 
  pd(30,30) = pd(30,30) + rrt(301) * density(35) 
  pd(30,35) = pd(30,35) + rrt(301) * density(30) 
  pd(35,30) = pd(35,30) - rrt(301) * density(35) 
  pd(35,35) = pd(35,35) - rrt(301) * density(30) 
  pd(30,01) = pd(30,01) + rrt(302) * density(35) 
  pd(30,35) = pd(30,35) + rrt(302) * density(01) 
  pd(35,01) = pd(35,01) - rrt(302) * density(35) 
  pd(35,35) = pd(35,35) - rrt(302) * density(01) 
  pd(30,35) = pd(30,35) + rrt(303) * density(50) 
  pd(30,50) = pd(30,50) + rrt(303) * density(35) 
  pd(35,35) = pd(35,35) - rrt(303) * density(50) 
  pd(35,50) = pd(35,50) - rrt(303) * density(35) 
  pd(30,35) = pd(30,35) + rrt(304) * density(41) * 2.0d0
  pd(30,41) = pd(30,41) + rrt(304) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(304) * density(41) 
  pd(35,41) = pd(35,41) - rrt(304) * density(35) 
  pd(39,35) = pd(39,35) + rrt(304) * density(41) 
  pd(39,41) = pd(39,41) + rrt(304) * density(35) 
  pd(41,35) = pd(41,35) - rrt(304) * density(41) 
  pd(41,41) = pd(41,41) - rrt(304) * density(35) 
  pd(30,35) = pd(30,35) + rrt(305) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(305) * density(35) * 4.0d0
  pd(36,35) = pd(36,35) + rrt(305) * density(35) * 2.0d0
  pd(30,38) = pd(30,38) + rrt(306) * density(41) 
  pd(30,41) = pd(30,41) + rrt(306) * density(38) 
  pd(35,38) = pd(35,38) + rrt(306) * density(41) 
  pd(35,41) = pd(35,41) + rrt(306) * density(38) 
  pd(38,38) = pd(38,38) - rrt(306) * density(41) 
  pd(38,41) = pd(38,41) - rrt(306) * density(38) 
  pd(41,38) = pd(41,38) - rrt(306) * density(41) 
  pd(41,41) = pd(41,41) - rrt(306) * density(38) 
  pd(35,36) = pd(35,36) + rrt(307) * density(38) 
  pd(35,38) = pd(35,38) + rrt(307) * density(36) 
  pd(36,36) = pd(36,36) - rrt(307) * density(38) 
  pd(36,38) = pd(36,38) - rrt(307) * density(36) 
  pd(30,36) = pd(30,36) + rrt(308) * density(38) 
  pd(30,38) = pd(30,38) + rrt(308) * density(36) 
  pd(36,36) = pd(36,36) - rrt(308) * density(38) 
  pd(36,38) = pd(36,38) - rrt(308) * density(36) 
  pd(38,36) = pd(38,36) - rrt(308) * density(38) 
  pd(38,38) = pd(38,38) - rrt(308) * density(36) 
  pd(39,36) = pd(39,36) + rrt(308) * density(38) 
  pd(39,38) = pd(39,38) + rrt(308) * density(36) 
  pd(35,30) = pd(35,30) + rrt(309) * density(36) 
  pd(35,36) = pd(35,36) + rrt(309) * density(30) 
  pd(36,30) = pd(36,30) - rrt(309) * density(36) 
  pd(36,36) = pd(36,36) - rrt(309) * density(30) 
  pd(35,01) = pd(35,01) + rrt(310) * density(36) 
  pd(35,36) = pd(35,36) + rrt(310) * density(01) 
  pd(36,01) = pd(36,01) - rrt(310) * density(36) 
  pd(36,36) = pd(36,36) - rrt(310) * density(01) 
  pd(35,36) = pd(35,36) + rrt(311) * density(50) 
  pd(35,50) = pd(35,50) + rrt(311) * density(36) 
  pd(36,36) = pd(36,36) - rrt(311) * density(50) 
  pd(36,50) = pd(36,50) - rrt(311) * density(36) 
  pd(30,36) = pd(30,36) + rrt(312) * density(41) * 2.0d0
  pd(30,41) = pd(30,41) + rrt(312) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(312) * density(41) 
  pd(36,41) = pd(36,41) - rrt(312) * density(36) 
  pd(38,36) = pd(38,36) + rrt(312) * density(41) 
  pd(38,41) = pd(38,41) + rrt(312) * density(36) 
  pd(41,36) = pd(41,36) - rrt(312) * density(41) 
  pd(41,41) = pd(41,41) - rrt(312) * density(36) 
  pd(30,37) = pd(30,37) + rrt(313) * density(38) 
  pd(30,38) = pd(30,38) + rrt(313) * density(37) 
  pd(37,37) = pd(37,37) - rrt(313) * density(38) 
  pd(37,38) = pd(37,38) - rrt(313) * density(37) 
  pd(38,37) = pd(38,37) - rrt(313) * density(38) 
  pd(38,38) = pd(38,38) - rrt(313) * density(37) 
  pd(40,37) = pd(40,37) + rrt(313) * density(38) 
  pd(40,38) = pd(40,38) + rrt(313) * density(37) 
  pd(30,30) = pd(30,30) - rrt(314) * density(37) 
  pd(30,37) = pd(30,37) - rrt(314) * density(30) 
  pd(36,30) = pd(36,30) + rrt(314) * density(37) * 2.0d0
  pd(36,37) = pd(36,37) + rrt(314) * density(30) * 2.0d0
  pd(37,30) = pd(37,30) - rrt(314) * density(37) 
  pd(37,37) = pd(37,37) - rrt(314) * density(30) 
  pd(36,01) = pd(36,01) + rrt(315) * density(37) 
  pd(36,37) = pd(36,37) + rrt(315) * density(01) 
  pd(37,01) = pd(37,01) - rrt(315) * density(37) 
  pd(37,37) = pd(37,37) - rrt(315) * density(01) 
  pd(38,38) = pd(38,38) + rrt(316) * density(39) 
  pd(38,39) = pd(38,39) + rrt(316) * density(38) 
  pd(39,38) = pd(39,38) - rrt(316) * density(39) 
  pd(39,39) = pd(39,39) - rrt(316) * density(38) 
  pd(38,30) = pd(38,30) + rrt(317) * density(39) 
  pd(38,39) = pd(38,39) + rrt(317) * density(30) 
  pd(39,30) = pd(39,30) - rrt(317) * density(39) 
  pd(39,39) = pd(39,39) - rrt(317) * density(30) 
  pd(30,30) = pd(30,30) - rrt(318) * density(39) 
  pd(30,39) = pd(30,39) - rrt(318) * density(30) 
  pd(35,30) = pd(35,30) + rrt(318) * density(39) 
  pd(35,39) = pd(35,39) + rrt(318) * density(30) 
  pd(38,30) = pd(38,30) + rrt(318) * density(39) 
  pd(38,39) = pd(38,39) + rrt(318) * density(30) 
  pd(39,30) = pd(39,30) - rrt(318) * density(39) 
  pd(39,39) = pd(39,39) - rrt(318) * density(30) 
  pd(30,30) = pd(30,30) - rrt(319) * density(39) 
  pd(30,39) = pd(30,39) - rrt(319) * density(30) 
  pd(36,30) = pd(36,30) + rrt(319) * density(39) 
  pd(36,39) = pd(36,39) + rrt(319) * density(30) 
  pd(38,30) = pd(38,30) + rrt(319) * density(39) 
  pd(38,39) = pd(38,39) + rrt(319) * density(30) 
  pd(39,30) = pd(39,30) - rrt(319) * density(39) 
  pd(39,39) = pd(39,39) - rrt(319) * density(30) 
  pd(38,01) = pd(38,01) + rrt(320) * density(39) 
  pd(38,39) = pd(38,39) + rrt(320) * density(01) 
  pd(39,01) = pd(39,01) - rrt(320) * density(39) 
  pd(39,39) = pd(39,39) - rrt(320) * density(01) 
  pd(30,39) = pd(30,39) + rrt(321) * density(41) 
  pd(30,41) = pd(30,41) + rrt(321) * density(39) 
  pd(38,39) = pd(38,39) + rrt(321) * density(41) * 2.0d0
  pd(38,41) = pd(38,41) + rrt(321) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(321) * density(41) 
  pd(39,41) = pd(39,41) - rrt(321) * density(39) 
  pd(41,39) = pd(41,39) - rrt(321) * density(41) 
  pd(41,41) = pd(41,41) - rrt(321) * density(39) 
  pd(30,39) = pd(30,39) + rrt(322) * density(41) * 2.0d0
  pd(30,41) = pd(30,41) + rrt(322) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(322) * density(41) 
  pd(39,41) = pd(39,41) - rrt(322) * density(39) 
  pd(41,39) = pd(41,39) - rrt(322) * density(41) 
  pd(41,41) = pd(41,41) - rrt(322) * density(39) 
  pd(23,39) = pd(23,39) + rrt(323) * density(50) 
  pd(23,50) = pd(23,50) + rrt(323) * density(39) 
  pd(30,39) = pd(30,39) + rrt(323) * density(50) 
  pd(30,50) = pd(30,50) + rrt(323) * density(39) 
  pd(39,39) = pd(39,39) - rrt(323) * density(50) 
  pd(39,50) = pd(39,50) - rrt(323) * density(39) 
  pd(50,39) = pd(50,39) - rrt(323) * density(50) 
  pd(50,50) = pd(50,50) - rrt(323) * density(39) 
  pd(39,39) = pd(39,39) - rrt(324) * density(51) 
  pd(39,51) = pd(39,51) - rrt(324) * density(39) 
  pd(50,39) = pd(50,39) + rrt(324) * density(51) * 2.0d0
  pd(50,51) = pd(50,51) + rrt(324) * density(39) * 2.0d0
  pd(51,39) = pd(51,39) - rrt(324) * density(51) 
  pd(51,51) = pd(51,51) - rrt(324) * density(39) 
  pd(01,39) = pd(01,39) + rrt(325) * density(51) 
  pd(01,51) = pd(01,51) + rrt(325) * density(39) 
  pd(30,39) = pd(30,39) + rrt(325) * density(51) 
  pd(30,51) = pd(30,51) + rrt(325) * density(39) 
  pd(39,39) = pd(39,39) - rrt(325) * density(51) 
  pd(39,51) = pd(39,51) - rrt(325) * density(39) 
  pd(51,39) = pd(51,39) - rrt(325) * density(51) 
  pd(51,51) = pd(51,51) - rrt(325) * density(39) 
  pd(39,38) = pd(39,38) + rrt(326) * density(40) 
  pd(39,40) = pd(39,40) + rrt(326) * density(38) 
  pd(40,38) = pd(40,38) - rrt(326) * density(40) 
  pd(40,40) = pd(40,40) - rrt(326) * density(38) 
  pd(38,23) = pd(38,23) + rrt(327) * density(40) 
  pd(38,40) = pd(38,40) + rrt(327) * density(23) 
  pd(40,23) = pd(40,23) - rrt(327) * density(40) 
  pd(40,40) = pd(40,40) - rrt(327) * density(23) 
  pd(39,30) = pd(39,30) + rrt(328) * density(40) 
  pd(39,40) = pd(39,40) + rrt(328) * density(30) 
  pd(40,30) = pd(40,30) - rrt(328) * density(40) 
  pd(40,40) = pd(40,40) - rrt(328) * density(30) 
  pd(30,30) = pd(30,30) - rrt(329) * density(40) 
  pd(30,40) = pd(30,40) - rrt(329) * density(30) 
  pd(38,30) = pd(38,30) + rrt(329) * density(40) * 3.0d0
  pd(38,40) = pd(38,40) + rrt(329) * density(30) * 3.0d0
  pd(40,30) = pd(40,30) - rrt(329) * density(40) 
  pd(40,40) = pd(40,40) - rrt(329) * density(30) 
  pd(38,01) = pd(38,01) + rrt(330) * density(40) 
  pd(38,40) = pd(38,40) + rrt(330) * density(01) 
  pd(40,01) = pd(40,01) - rrt(330) * density(40) 
  pd(40,40) = pd(40,40) - rrt(330) * density(01) 
  pd(35,35) = pd(35,35) - rrt(331) * density(40) 
  pd(35,40) = pd(35,40) - rrt(331) * density(35) 
  pd(37,35) = pd(37,35) + rrt(331) * density(40) 
  pd(37,40) = pd(37,40) + rrt(331) * density(35) 
  pd(38,35) = pd(38,35) + rrt(331) * density(40) 
  pd(38,40) = pd(38,40) + rrt(331) * density(35) 
  pd(40,35) = pd(40,35) - rrt(331) * density(40) 
  pd(40,40) = pd(40,40) - rrt(331) * density(35) 
  pd(35,35) = pd(35,35) - rrt(332) * density(40) 
  pd(35,40) = pd(35,40) - rrt(332) * density(35) 
  pd(36,35) = pd(36,35) + rrt(332) * density(40) 
  pd(36,40) = pd(36,40) + rrt(332) * density(35) 
  pd(39,35) = pd(39,35) + rrt(332) * density(40) 
  pd(39,40) = pd(39,40) + rrt(332) * density(35) 
  pd(40,35) = pd(40,35) - rrt(332) * density(40) 
  pd(40,40) = pd(40,40) - rrt(332) * density(35) 
  pd(35,35) = pd(35,35) - rrt(333) * density(40) 
  pd(35,40) = pd(35,40) - rrt(333) * density(35) 
  pd(38,35) = pd(38,35) + rrt(333) * density(40) * 3.0d0
  pd(38,40) = pd(38,40) + rrt(333) * density(35) * 3.0d0
  pd(40,35) = pd(40,35) - rrt(333) * density(40) 
  pd(40,40) = pd(40,40) - rrt(333) * density(35) 
  pd(38,40) = pd(38,40) + rrt(334) * density(50) 
  pd(38,50) = pd(38,50) + rrt(334) * density(40) 
  pd(40,40) = pd(40,40) - rrt(334) * density(50) 
  pd(40,50) = pd(40,50) - rrt(334) * density(40) 
  pd(39,40) = pd(39,40) + rrt(335) * density(50) 
  pd(39,50) = pd(39,50) + rrt(335) * density(40) 
  pd(40,40) = pd(40,40) - rrt(335) * density(50) 
  pd(40,50) = pd(40,50) - rrt(335) * density(40) 
  pd(30,40) = pd(30,40) + rrt(336) * density(41) * 2.0d0
  pd(30,41) = pd(30,41) + rrt(336) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(336) * density(41) 
  pd(40,41) = pd(40,41) - rrt(336) * density(40) 
  pd(41,40) = pd(41,40) - rrt(336) * density(41) 
  pd(41,41) = pd(41,41) - rrt(336) * density(40) 
  pd(30,40) = pd(30,40) + rrt(337) * density(41) 
  pd(30,41) = pd(30,41) + rrt(337) * density(40) 
  pd(38,40) = pd(38,40) + rrt(337) * density(41) 
  pd(38,41) = pd(38,41) + rrt(337) * density(40) 
  pd(39,40) = pd(39,40) + rrt(337) * density(41) 
  pd(39,41) = pd(39,41) + rrt(337) * density(40) 
  pd(40,40) = pd(40,40) - rrt(337) * density(41) 
  pd(40,41) = pd(40,41) - rrt(337) * density(40) 
  pd(41,40) = pd(41,40) - rrt(337) * density(41) 
  pd(41,41) = pd(41,41) - rrt(337) * density(40) 
  pd(38,40) = pd(38,40) + rrt(338) * density(51) 
  pd(38,51) = pd(38,51) + rrt(338) * density(40) 
  pd(40,40) = pd(40,40) - rrt(338) * density(51) 
  pd(40,51) = pd(40,51) - rrt(338) * density(40) 
  pd(39,40) = pd(39,40) + rrt(339) * density(51) 
  pd(39,51) = pd(39,51) + rrt(339) * density(40) 
  pd(40,40) = pd(40,40) - rrt(339) * density(51) 
  pd(40,51) = pd(40,51) - rrt(339) * density(40) 
  pd(01,23) = pd(01,23) + rrt(340) * density(50) 
  pd(01,50) = pd(01,50) + rrt(340) * density(23) 
  pd(23,23) = pd(23,23) - rrt(340) * density(50) 
  pd(23,50) = pd(23,50) - rrt(340) * density(23) 
  pd(38,23) = pd(38,23) + rrt(340) * density(50) 
  pd(38,50) = pd(38,50) + rrt(340) * density(23) 
  pd(50,23) = pd(50,23) - rrt(340) * density(50) 
  pd(50,50) = pd(50,50) - rrt(340) * density(23) 
  pd(23,23) = pd(23,23) - rrt(341) * density(30) 
  pd(23,30) = pd(23,30) - rrt(341) * density(23) 
  pd(30,23) = pd(30,23) - rrt(341) * density(30) 
  pd(30,30) = pd(30,30) - rrt(341) * density(23) 
  pd(38,23) = pd(38,23) + rrt(341) * density(30) 
  pd(38,30) = pd(38,30) + rrt(341) * density(23) 
  pd(50,23) = pd(50,23) + rrt(341) * density(30) 
  pd(50,30) = pd(50,30) + rrt(341) * density(23) 
  pd(01,23) = pd(01,23) + rrt(342) * density(52) 
  pd(01,52) = pd(01,52) + rrt(342) * density(23) 
  pd(23,23) = pd(23,23) - rrt(342) * density(52) 
  pd(23,52) = pd(23,52) - rrt(342) * density(23) 
  pd(38,23) = pd(38,23) + rrt(342) * density(52) * 2.0d0
  pd(38,52) = pd(38,52) + rrt(342) * density(23) * 2.0d0
  pd(52,23) = pd(52,23) - rrt(342) * density(52) 
  pd(52,52) = pd(52,52) - rrt(342) * density(23) 
  pd(23,23) = pd(23,23) - rrt(343) * density(52) 
  pd(23,52) = pd(23,52) - rrt(343) * density(23) 
  pd(38,23) = pd(38,23) + rrt(343) * density(52) 
  pd(38,52) = pd(38,52) + rrt(343) * density(23) 
  pd(51,23) = pd(51,23) + rrt(343) * density(52) 
  pd(51,52) = pd(51,52) + rrt(343) * density(23) 
  pd(52,23) = pd(52,23) - rrt(343) * density(52) 
  pd(52,52) = pd(52,52) - rrt(343) * density(23) 
  pd(01,23) = pd(01,23) + rrt(344) * density(52) 
  pd(01,52) = pd(01,52) + rrt(344) * density(23) 
  pd(23,23) = pd(23,23) - rrt(344) * density(52) 
  pd(23,52) = pd(23,52) - rrt(344) * density(23) 
  pd(30,23) = pd(30,23) + rrt(344) * density(52) 
  pd(30,52) = pd(30,52) + rrt(344) * density(23) 
  pd(52,23) = pd(52,23) - rrt(344) * density(52) 
  pd(52,52) = pd(52,52) - rrt(344) * density(23) 
  pd(23,23) = pd(23,23) - rrt(345) * density(52) 
  pd(23,52) = pd(23,52) - rrt(345) * density(23) 
  pd(50,23) = pd(50,23) + rrt(345) * density(52) * 2.0d0
  pd(50,52) = pd(50,52) + rrt(345) * density(23) * 2.0d0
  pd(52,23) = pd(52,23) - rrt(345) * density(52) 
  pd(52,52) = pd(52,52) - rrt(345) * density(23) 
  pd(01,01) = pd(01,01) - rrt(346) * density(38) 
  pd(01,38) = pd(01,38) - rrt(346) * density(01) 
  pd(23,01) = pd(23,01) + rrt(346) * density(38) 
  pd(23,38) = pd(23,38) + rrt(346) * density(01) 
  pd(38,01) = pd(38,01) - rrt(346) * density(38) 
  pd(38,38) = pd(38,38) - rrt(346) * density(01) 
  pd(50,01) = pd(50,01) + rrt(346) * density(38) 
  pd(50,38) = pd(50,38) + rrt(346) * density(01) 
  pd(23,38) = pd(23,38) + rrt(347) * density(50) 
  pd(23,50) = pd(23,50) + rrt(347) * density(38) 
  pd(30,38) = pd(30,38) + rrt(347) * density(50) 
  pd(30,50) = pd(30,50) + rrt(347) * density(38) 
  pd(38,38) = pd(38,38) - rrt(347) * density(50) 
  pd(38,50) = pd(38,50) - rrt(347) * density(38) 
  pd(50,38) = pd(50,38) - rrt(347) * density(50) 
  pd(50,50) = pd(50,50) - rrt(347) * density(38) 
  pd(38,38) = pd(38,38) - rrt(348) * density(50) 
  pd(38,50) = pd(38,50) - rrt(348) * density(38) 
  pd(50,38) = pd(50,38) - rrt(348) * density(50) 
  pd(50,50) = pd(50,50) - rrt(348) * density(38) 
  pd(52,38) = pd(52,38) + rrt(348) * density(50) 
  pd(52,50) = pd(52,50) + rrt(348) * density(38) 
  pd(01,38) = pd(01,38) + rrt(349) * density(51) 
  pd(01,51) = pd(01,51) + rrt(349) * density(38) 
  pd(30,38) = pd(30,38) + rrt(349) * density(51) 
  pd(30,51) = pd(30,51) + rrt(349) * density(38) 
  pd(38,38) = pd(38,38) - rrt(349) * density(51) 
  pd(38,51) = pd(38,51) - rrt(349) * density(38) 
  pd(51,38) = pd(51,38) - rrt(349) * density(51) 
  pd(51,51) = pd(51,51) - rrt(349) * density(38) 
  pd(38,38) = pd(38,38) - rrt(350) * density(51) 
  pd(38,51) = pd(38,51) - rrt(350) * density(38) 
  pd(50,38) = pd(50,38) + rrt(350) * density(51) * 2.0d0
  pd(50,51) = pd(50,51) + rrt(350) * density(38) * 2.0d0
  pd(51,38) = pd(51,38) - rrt(350) * density(51) 
  pd(51,51) = pd(51,51) - rrt(350) * density(38) 
  pd(30,38) = pd(30,38) + rrt(351) * density(52) 
  pd(30,52) = pd(30,52) + rrt(351) * density(38) 
  pd(38,38) = pd(38,38) - rrt(351) * density(52) 
  pd(38,52) = pd(38,52) - rrt(351) * density(38) 
  pd(50,38) = pd(50,38) + rrt(351) * density(52) 
  pd(50,52) = pd(50,52) + rrt(351) * density(38) 
  pd(52,38) = pd(52,38) - rrt(351) * density(52) 
  pd(52,52) = pd(52,52) - rrt(351) * density(38) 
  pd(30,38) = pd(30,38) + rrt(352) * density(53) 
  pd(30,53) = pd(30,53) + rrt(352) * density(38) 
  pd(38,38) = pd(38,38) - rrt(352) * density(53) 
  pd(38,53) = pd(38,53) - rrt(352) * density(38) 
  pd(52,38) = pd(52,38) + rrt(352) * density(53) 
  pd(52,53) = pd(52,53) + rrt(352) * density(38) 
  pd(53,38) = pd(53,38) - rrt(352) * density(53) 
  pd(53,53) = pd(53,53) - rrt(352) * density(38) 
  pd(01,01) = pd(01,01) - rrt(353) * density(30) 
  pd(01,30) = pd(01,30) - rrt(353) * density(01) 
  pd(30,01) = pd(30,01) - rrt(353) * density(30) 
  pd(30,30) = pd(30,30) - rrt(353) * density(01) 
  pd(38,01) = pd(38,01) + rrt(353) * density(30) 
  pd(38,30) = pd(38,30) + rrt(353) * density(01) 
  pd(51,01) = pd(51,01) + rrt(353) * density(30) 
  pd(51,30) = pd(51,30) + rrt(353) * density(01) 
  pd(23,50) = pd(23,50) + rrt(354) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(354) * density(50) * 4.0d0
  pd(52,50) = pd(52,50) + rrt(354) * density(50) * 2.0d0
  pd(38,50) = pd(38,50) + rrt(355) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(355) * density(50) * 4.0d0
  pd(51,50) = pd(51,50) + rrt(355) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(356) * density(50) * 2.0d0
  pd(30,50) = pd(30,50) + rrt(356) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(356) * density(50) * 4.0d0
  pd(30,30) = pd(30,30) - rrt(357) * density(50) 
  pd(30,50) = pd(30,50) - rrt(357) * density(30) 
  pd(38,30) = pd(38,30) + rrt(357) * density(50) 
  pd(38,50) = pd(38,50) + rrt(357) * density(30) 
  pd(50,30) = pd(50,30) - rrt(357) * density(50) 
  pd(50,50) = pd(50,50) - rrt(357) * density(30) 
  pd(52,30) = pd(52,30) + rrt(357) * density(50) 
  pd(52,50) = pd(52,50) + rrt(357) * density(30) 
  pd(30,41) = pd(30,41) + rrt(358) * density(50) 
  pd(30,50) = pd(30,50) + rrt(358) * density(41) 
  pd(41,41) = pd(41,41) - rrt(358) * density(50) 
  pd(41,50) = pd(41,50) - rrt(358) * density(41) 
  pd(50,41) = pd(50,41) - rrt(358) * density(50) 
  pd(50,50) = pd(50,50) - rrt(358) * density(41) 
  pd(52,41) = pd(52,41) + rrt(358) * density(50) 
  pd(52,50) = pd(52,50) + rrt(358) * density(41) 
  pd(01,50) = pd(01,50) + rrt(359) * density(51) 
  pd(01,51) = pd(01,51) + rrt(359) * density(50) 
  pd(50,50) = pd(50,50) - rrt(359) * density(51) 
  pd(50,51) = pd(50,51) - rrt(359) * density(50) 
  pd(51,50) = pd(51,50) - rrt(359) * density(51) 
  pd(51,51) = pd(51,51) - rrt(359) * density(50) 
  pd(52,50) = pd(52,50) + rrt(359) * density(51) 
  pd(52,51) = pd(52,51) + rrt(359) * density(50) 
  pd(50,50) = pd(50,50) - rrt(360) * density(53) 
  pd(50,53) = pd(50,53) - rrt(360) * density(50) 
  pd(52,50) = pd(52,50) + rrt(360) * density(53) * 2.0d0
  pd(52,53) = pd(52,53) + rrt(360) * density(50) * 2.0d0
  pd(53,50) = pd(53,50) - rrt(360) * density(53) 
  pd(53,53) = pd(53,53) - rrt(360) * density(50) 
  pd(30,30) = pd(30,30) - rrt(361) * density(30) * 4.0d0
  pd(38,30) = pd(38,30) + rrt(361) * density(30) * 2.0d0
  pd(41,30) = pd(41,30) + rrt(361) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(362) * density(52) 
  pd(30,52) = pd(30,52) - rrt(362) * density(30) 
  pd(41,30) = pd(41,30) + rrt(362) * density(52) 
  pd(41,52) = pd(41,52) + rrt(362) * density(30) 
  pd(50,30) = pd(50,30) + rrt(362) * density(52) 
  pd(50,52) = pd(50,52) + rrt(362) * density(30) 
  pd(52,30) = pd(52,30) - rrt(362) * density(52) 
  pd(52,52) = pd(52,52) - rrt(362) * density(30) 
  pd(30,52) = pd(30,52) + rrt(363) * density(52) * 2.0d0
  pd(50,52) = pd(50,52) + rrt(363) * density(52) * 4.0d0
  pd(52,52) = pd(52,52) - rrt(363) * density(52) * 4.0d0
  pd(50,52) = pd(50,52) + rrt(364) * density(52) * 2.0d0
  pd(52,52) = pd(52,52) - rrt(364) * density(52) * 4.0d0
  pd(53,52) = pd(53,52) + rrt(364) * density(52) * 2.0d0
  pd(30,41) = pd(30,41) + rrt(365) * density(52) 
  pd(30,52) = pd(30,52) + rrt(365) * density(41) 
  pd(41,41) = pd(41,41) - rrt(365) * density(52) 
  pd(41,52) = pd(41,52) - rrt(365) * density(41) 
  pd(52,41) = pd(52,41) - rrt(365) * density(52) 
  pd(52,52) = pd(52,52) - rrt(365) * density(41) 
  pd(53,41) = pd(53,41) + rrt(365) * density(52) 
  pd(53,52) = pd(53,52) + rrt(365) * density(41) 
  pd(30,52) = pd(30,52) + rrt(366) * density(53) 
  pd(30,53) = pd(30,53) + rrt(366) * density(52) 
  pd(50,52) = pd(50,52) + rrt(366) * density(53) 
  pd(50,53) = pd(50,53) + rrt(366) * density(52) 
  pd(53,52) = pd(53,52) - rrt(366) * density(53) 
  pd(53,53) = pd(53,53) - rrt(366) * density(52) 
  pd(30,30) = pd(30,30) - rrt(367) * density(53) 
  pd(30,53) = pd(30,53) - rrt(367) * density(30) 
  pd(41,30) = pd(41,30) + rrt(367) * density(53) 
  pd(41,53) = pd(41,53) + rrt(367) * density(30) 
  pd(52,30) = pd(52,30) + rrt(367) * density(53) 
  pd(52,53) = pd(52,53) + rrt(367) * density(30) 
  pd(53,30) = pd(53,30) - rrt(367) * density(53) 
  pd(53,53) = pd(53,53) - rrt(367) * density(30) 
  pd(30,53) = pd(30,53) + rrt(368) * density(53) * 2.0d0
  pd(52,53) = pd(52,53) + rrt(368) * density(53) * 4.0d0
  pd(53,53) = pd(53,53) - rrt(368) * density(53) * 4.0d0
  pd(23,23) = pd(23,23) - rrt(369) * density(23) * 4.0d0
  pd(27,23) = pd(27,23) + rrt(369) * density(23) * 2.0d0
  pd(79,23) = pd(79,23) + rrt(369) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(370) * density(38) 
  pd(23,38) = pd(23,38) - rrt(370) * density(23) 
  pd(38,23) = pd(38,23) - rrt(370) * density(38) 
  pd(38,38) = pd(38,38) - rrt(370) * density(23) 
  pd(55,23) = pd(55,23) + rrt(370) * density(38) 
  pd(55,38) = pd(55,38) + rrt(370) * density(23) 
  pd(79,23) = pd(79,23) + rrt(370) * density(38) 
  pd(79,38) = pd(79,38) + rrt(370) * density(23) 
  pd(23,63) = pd(23,63) + rrt(371) * density(72) 
  pd(23,72) = pd(23,72) + rrt(371) * density(63) 
  pd(62,63) = pd(62,63) + rrt(371) * density(72) 
  pd(62,72) = pd(62,72) + rrt(371) * density(63) 
  pd(63,63) = pd(63,63) - rrt(371) * density(72) 
  pd(63,72) = pd(63,72) - rrt(371) * density(63) 
  pd(72,63) = pd(72,63) - rrt(371) * density(72) 
  pd(72,72) = pd(72,72) - rrt(371) * density(63) 
  pd(62,63) = pd(62,63) + rrt(372) * density(71) 
  pd(62,71) = pd(62,71) + rrt(372) * density(63) 
  pd(63,63) = pd(63,63) - rrt(372) * density(71) 
  pd(63,71) = pd(63,71) - rrt(372) * density(63) 
  pd(71,63) = pd(71,63) - rrt(372) * density(71) 
  pd(71,71) = pd(71,71) - rrt(372) * density(63) 
  pd(72,63) = pd(72,63) + rrt(372) * density(71) 
  pd(72,71) = pd(72,71) + rrt(372) * density(63) 
  pd(62,63) = pd(62,63) + rrt(373) * density(70) 
  pd(62,70) = pd(62,70) + rrt(373) * density(63) 
  pd(63,63) = pd(63,63) - rrt(373) * density(70) 
  pd(63,70) = pd(63,70) - rrt(373) * density(63) 
  pd(70,63) = pd(70,63) - rrt(373) * density(70) 
  pd(70,70) = pd(70,70) - rrt(373) * density(63) 
  pd(71,63) = pd(71,63) + rrt(373) * density(70) 
  pd(71,70) = pd(71,70) + rrt(373) * density(63) 
  pd(62,62) = pd(62,62) - rrt(374) * density(71) 
  pd(62,71) = pd(62,71) - rrt(374) * density(62) 
  pd(63,62) = pd(63,62) + rrt(374) * density(71) 
  pd(63,71) = pd(63,71) + rrt(374) * density(62) 
  pd(70,62) = pd(70,62) + rrt(374) * density(71) 
  pd(70,71) = pd(70,71) + rrt(374) * density(62) 
  pd(71,62) = pd(71,62) - rrt(374) * density(71) 
  pd(71,71) = pd(71,71) - rrt(374) * density(62) 
  pd(01,23) = pd(01,23) + rrt(375) * density(72) 
  pd(01,72) = pd(01,72) + rrt(375) * density(23) 
  pd(23,23) = pd(23,23) - rrt(375) * density(72) 
  pd(23,72) = pd(23,72) - rrt(375) * density(23) 
  pd(63,23) = pd(63,23) + rrt(375) * density(72) 
  pd(63,72) = pd(63,72) + rrt(375) * density(23) 
  pd(72,23) = pd(72,23) - rrt(375) * density(72) 
  pd(72,72) = pd(72,72) - rrt(375) * density(23) 
  pd(01,23) = pd(01,23) + rrt(376) * density(71) 
  pd(01,71) = pd(01,71) + rrt(376) * density(23) 
  pd(23,23) = pd(23,23) - rrt(376) * density(71) 
  pd(23,71) = pd(23,71) - rrt(376) * density(23) 
  pd(63,23) = pd(63,23) + rrt(376) * density(71) * 2.0d0
  pd(63,71) = pd(63,71) + rrt(376) * density(23) * 2.0d0
  pd(71,23) = pd(71,23) - rrt(376) * density(71) 
  pd(71,71) = pd(71,71) - rrt(376) * density(23) 
  pd(01,23) = pd(01,23) + rrt(377) * density(71) 
  pd(01,71) = pd(01,71) + rrt(377) * density(23) 
  pd(23,23) = pd(23,23) - rrt(377) * density(71) 
  pd(23,71) = pd(23,71) - rrt(377) * density(23) 
  pd(62,23) = pd(62,23) + rrt(377) * density(71) 
  pd(62,71) = pd(62,71) + rrt(377) * density(23) 
  pd(71,23) = pd(71,23) - rrt(377) * density(71) 
  pd(71,71) = pd(71,71) - rrt(377) * density(23) 
  pd(23,72) = pd(23,72) + rrt(378) * density(72) * 2.0d0
  pd(71,72) = pd(71,72) + rrt(378) * density(72) * 2.0d0
  pd(72,72) = pd(72,72) - rrt(378) * density(72) * 4.0d0
  pd(01,72) = pd(01,72) + rrt(379) * density(72) * 2.0d0
  pd(63,72) = pd(63,72) + rrt(379) * density(72) * 4.0d0
  pd(72,72) = pd(72,72) - rrt(379) * density(72) * 4.0d0
  pd(01,72) = pd(01,72) + rrt(380) * density(72) * 2.0d0
  pd(62,72) = pd(62,72) + rrt(380) * density(72) * 2.0d0
  pd(72,72) = pd(72,72) - rrt(380) * density(72) * 4.0d0
  pd(23,71) = pd(23,71) + rrt(381) * density(72) 
  pd(23,72) = pd(23,72) + rrt(381) * density(71) 
  pd(70,71) = pd(70,71) + rrt(381) * density(72) 
  pd(70,72) = pd(70,72) + rrt(381) * density(71) 
  pd(71,71) = pd(71,71) - rrt(381) * density(72) 
  pd(71,72) = pd(71,72) - rrt(381) * density(71) 
  pd(72,71) = pd(72,71) - rrt(381) * density(72) 
  pd(72,72) = pd(72,72) - rrt(381) * density(71) 
  pd(63,63) = pd(63,63) - rrt(382) * density(71) 
  pd(63,71) = pd(63,71) - rrt(382) * density(63) 
  pd(70,63) = pd(70,63) + rrt(382) * density(71) 
  pd(70,71) = pd(70,71) + rrt(382) * density(63) 
  pd(71,63) = pd(71,63) - rrt(382) * density(71) 
  pd(71,71) = pd(71,71) - rrt(382) * density(63) 
  pd(70,71) = pd(70,71) + rrt(383) * density(71) * 2.0d0
  pd(71,71) = pd(71,71) - rrt(383) * density(71) * 4.0d0
  pd(72,71) = pd(72,71) + rrt(383) * density(71) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(384) * density(01) * 2.0d0
  pd(23,01) = pd(23,01) + rrt(384) * density(01) * 4.0d0
  pd(01,01) = pd(01,01) - rrt(385) * density(30) 
  pd(01,30) = pd(01,30) - rrt(385) * density(01) 
  pd(23,01) = pd(23,01) + rrt(385) * density(30) * 2.0d0
  pd(23,30) = pd(23,30) + rrt(385) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(386) * density(50) 
  pd(01,50) = pd(01,50) - rrt(386) * density(01) 
  pd(23,01) = pd(23,01) + rrt(386) * density(50) * 2.0d0
  pd(23,50) = pd(23,50) + rrt(386) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(387) * density(38) 
  pd(01,38) = pd(01,38) - rrt(387) * density(01) 
  pd(23,01) = pd(23,01) + rrt(387) * density(38) * 2.0d0
  pd(23,38) = pd(23,38) + rrt(387) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(388) * density(23) 
  pd(01,23) = pd(01,23) - rrt(388) * density(01) 
  pd(23,01) = pd(23,01) + rrt(388) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) + rrt(388) * density(01) * 2.0d0
  pd(30,01) = pd(30,01) - rrt(389) * density(30) 
  pd(30,30) = pd(30,30) - rrt(389) * density(01) 
  pd(38,01) = pd(38,01) + rrt(389) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(389) * density(01) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(390) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(390) * density(30) * 4.0d0
  pd(30,30) = pd(30,30) - rrt(391) * density(38) 
  pd(30,38) = pd(30,38) - rrt(391) * density(30) 
  pd(38,30) = pd(38,30) + rrt(391) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) + rrt(391) * density(30) * 2.0d0
  pd(30,23) = pd(30,23) - rrt(392) * density(30) 
  pd(30,30) = pd(30,30) - rrt(392) * density(23) 
  pd(38,23) = pd(38,23) + rrt(392) * density(30) * 2.0d0
  pd(38,30) = pd(38,30) + rrt(392) * density(23) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(393) * density(50) 
  pd(30,50) = pd(30,50) - rrt(393) * density(30) 
  pd(38,30) = pd(38,30) + rrt(393) * density(50) * 2.0d0
  pd(38,50) = pd(38,50) + rrt(393) * density(30) * 2.0d0
  pd(23,01) = pd(23,01) + rrt(394) * density(50) 
  pd(23,50) = pd(23,50) + rrt(394) * density(01) 
  pd(38,01) = pd(38,01) + rrt(394) * density(50) 
  pd(38,50) = pd(38,50) + rrt(394) * density(01) 
  pd(50,01) = pd(50,01) - rrt(394) * density(50) 
  pd(50,50) = pd(50,50) - rrt(394) * density(01) 
  pd(23,30) = pd(23,30) + rrt(395) * density(50) 
  pd(23,50) = pd(23,50) + rrt(395) * density(30) 
  pd(38,30) = pd(38,30) + rrt(395) * density(50) 
  pd(38,50) = pd(38,50) + rrt(395) * density(30) 
  pd(50,30) = pd(50,30) - rrt(395) * density(50) 
  pd(50,50) = pd(50,50) - rrt(395) * density(30) 
  pd(23,38) = pd(23,38) + rrt(396) * density(50) 
  pd(23,50) = pd(23,50) + rrt(396) * density(38) 
  pd(38,38) = pd(38,38) + rrt(396) * density(50) 
  pd(38,50) = pd(38,50) + rrt(396) * density(38) 
  pd(50,38) = pd(50,38) - rrt(396) * density(50) 
  pd(50,50) = pd(50,50) - rrt(396) * density(38) 
  pd(23,23) = pd(23,23) + rrt(397) * density(50) 
  pd(23,50) = pd(23,50) + rrt(397) * density(23) 
  pd(38,23) = pd(38,23) + rrt(397) * density(50) 
  pd(38,50) = pd(38,50) + rrt(397) * density(23) 
  pd(50,23) = pd(50,23) - rrt(397) * density(50) 
  pd(50,50) = pd(50,50) - rrt(397) * density(23) 
  pd(23,50) = pd(23,50) + rrt(398) * density(50) * 2.0d0
  pd(38,50) = pd(38,50) + rrt(398) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(398) * density(50) * 2.0d0
  pd(30,01) = pd(30,01) + rrt(399) * density(41) 
  pd(30,41) = pd(30,41) + rrt(399) * density(01) 
  pd(38,01) = pd(38,01) + rrt(399) * density(41) 
  pd(38,41) = pd(38,41) + rrt(399) * density(01) 
  pd(41,01) = pd(41,01) - rrt(399) * density(41) 
  pd(41,41) = pd(41,41) - rrt(399) * density(01) 
  pd(30,30) = pd(30,30) + rrt(400) * density(41) 
  pd(30,41) = pd(30,41) + rrt(400) * density(30) 
  pd(38,30) = pd(38,30) + rrt(400) * density(41) 
  pd(38,41) = pd(38,41) + rrt(400) * density(30) 
  pd(41,30) = pd(41,30) - rrt(400) * density(41) 
  pd(41,41) = pd(41,41) - rrt(400) * density(30) 
  pd(30,23) = pd(30,23) + rrt(401) * density(41) 
  pd(30,41) = pd(30,41) + rrt(401) * density(23) 
  pd(38,23) = pd(38,23) + rrt(401) * density(41) 
  pd(38,41) = pd(38,41) + rrt(401) * density(23) 
  pd(41,23) = pd(41,23) - rrt(401) * density(41) 
  pd(41,41) = pd(41,41) - rrt(401) * density(23) 
  pd(30,38) = pd(30,38) + rrt(402) * density(41) 
  pd(30,41) = pd(30,41) + rrt(402) * density(38) 
  pd(38,38) = pd(38,38) + rrt(402) * density(41) 
  pd(38,41) = pd(38,41) + rrt(402) * density(38) 
  pd(41,38) = pd(41,38) - rrt(402) * density(41) 
  pd(41,41) = pd(41,41) - rrt(402) * density(38) 
  pd(01,01) = pd(01,01) + rrt(403) * density(51) 
  pd(01,51) = pd(01,51) + rrt(403) * density(01) 
  pd(38,01) = pd(38,01) + rrt(403) * density(51) 
  pd(38,51) = pd(38,51) + rrt(403) * density(01) 
  pd(51,01) = pd(51,01) - rrt(403) * density(51) 
  pd(51,51) = pd(51,51) - rrt(403) * density(01) 
  pd(01,30) = pd(01,30) + rrt(404) * density(51) 
  pd(01,51) = pd(01,51) + rrt(404) * density(30) 
  pd(38,30) = pd(38,30) + rrt(404) * density(51) 
  pd(38,51) = pd(38,51) + rrt(404) * density(30) 
  pd(51,30) = pd(51,30) - rrt(404) * density(51) 
  pd(51,51) = pd(51,51) - rrt(404) * density(30) 
  pd(01,50) = pd(01,50) + rrt(405) * density(51) 
  pd(01,51) = pd(01,51) + rrt(405) * density(50) 
  pd(38,50) = pd(38,50) + rrt(405) * density(51) 
  pd(38,51) = pd(38,51) + rrt(405) * density(50) 
  pd(51,50) = pd(51,50) - rrt(405) * density(51) 
  pd(51,51) = pd(51,51) - rrt(405) * density(50) 
  pd(01,51) = pd(01,51) + rrt(406) * density(51) * 2.0d0
  pd(38,51) = pd(38,51) + rrt(406) * density(51) * 2.0d0
  pd(51,51) = pd(51,51) - rrt(406) * density(51) * 2.0d0
  pd(38,01) = pd(38,01) + rrt(407) * density(52) 
  pd(38,52) = pd(38,52) + rrt(407) * density(01) 
  pd(50,01) = pd(50,01) + rrt(407) * density(52) 
  pd(50,52) = pd(50,52) + rrt(407) * density(01) 
  pd(52,01) = pd(52,01) - rrt(407) * density(52) 
  pd(52,52) = pd(52,52) - rrt(407) * density(01) 
  pd(38,30) = pd(38,30) + rrt(408) * density(52) 
  pd(38,52) = pd(38,52) + rrt(408) * density(30) 
  pd(50,30) = pd(50,30) + rrt(408) * density(52) 
  pd(50,52) = pd(50,52) + rrt(408) * density(30) 
  pd(52,30) = pd(52,30) - rrt(408) * density(52) 
  pd(52,52) = pd(52,52) - rrt(408) * density(30) 
  pd(38,50) = pd(38,50) + rrt(409) * density(52) 
  pd(38,52) = pd(38,52) + rrt(409) * density(50) 
  pd(50,50) = pd(50,50) + rrt(409) * density(52) 
  pd(50,52) = pd(50,52) + rrt(409) * density(50) 
  pd(52,50) = pd(52,50) - rrt(409) * density(52) 
  pd(52,52) = pd(52,52) - rrt(409) * density(50) 
  pd(38,52) = pd(38,52) + rrt(410) * density(52) * 2.0d0
  pd(50,52) = pd(50,52) + rrt(410) * density(52) * 2.0d0
  pd(52,52) = pd(52,52) - rrt(410) * density(52) * 2.0d0
  pd(38,01) = pd(38,01) + rrt(411) * density(53) 
  pd(38,53) = pd(38,53) + rrt(411) * density(01) 
  pd(52,01) = pd(52,01) + rrt(411) * density(53) 
  pd(52,53) = pd(52,53) + rrt(411) * density(01) 
  pd(53,01) = pd(53,01) - rrt(411) * density(53) 
  pd(53,53) = pd(53,53) - rrt(411) * density(01) 
  pd(38,30) = pd(38,30) + rrt(412) * density(53) 
  pd(38,53) = pd(38,53) + rrt(412) * density(30) 
  pd(52,30) = pd(52,30) + rrt(412) * density(53) 
  pd(52,53) = pd(52,53) + rrt(412) * density(30) 
  pd(53,30) = pd(53,30) - rrt(412) * density(53) 
  pd(53,53) = pd(53,53) - rrt(412) * density(30) 
  pd(38,50) = pd(38,50) + rrt(413) * density(53) 
  pd(38,53) = pd(38,53) + rrt(413) * density(50) 
  pd(52,50) = pd(52,50) + rrt(413) * density(53) 
  pd(52,53) = pd(52,53) + rrt(413) * density(50) 
  pd(53,50) = pd(53,50) - rrt(413) * density(53) 
  pd(53,53) = pd(53,53) - rrt(413) * density(50) 
  pd(38,23) = pd(38,23) + rrt(414) * density(53) 
  pd(38,53) = pd(38,53) + rrt(414) * density(23) 
  pd(52,23) = pd(52,23) + rrt(414) * density(53) 
  pd(52,53) = pd(52,53) + rrt(414) * density(23) 
  pd(53,23) = pd(53,23) - rrt(414) * density(53) 
  pd(53,53) = pd(53,53) - rrt(414) * density(23) 
  pd(38,38) = pd(38,38) + rrt(415) * density(53) 
  pd(38,53) = pd(38,53) + rrt(415) * density(38) 
  pd(52,38) = pd(52,38) + rrt(415) * density(53) 
  pd(52,53) = pd(52,53) + rrt(415) * density(38) 
  pd(53,38) = pd(53,38) - rrt(415) * density(53) 
  pd(53,53) = pd(53,53) - rrt(415) * density(38) 
  pd(30,01) = pd(30,01) + rrt(416) * density(53) 
  pd(30,53) = pd(30,53) + rrt(416) * density(01) 
  pd(50,01) = pd(50,01) + rrt(416) * density(53) 
  pd(50,53) = pd(50,53) + rrt(416) * density(01) 
  pd(53,01) = pd(53,01) - rrt(416) * density(53) 
  pd(53,53) = pd(53,53) - rrt(416) * density(01) 
  pd(30,30) = pd(30,30) + rrt(417) * density(53) 
  pd(30,53) = pd(30,53) + rrt(417) * density(30) 
  pd(50,30) = pd(50,30) + rrt(417) * density(53) 
  pd(50,53) = pd(50,53) + rrt(417) * density(30) 
  pd(53,30) = pd(53,30) - rrt(417) * density(53) 
  pd(53,53) = pd(53,53) - rrt(417) * density(30) 
  pd(30,50) = pd(30,50) + rrt(418) * density(53) 
  pd(30,53) = pd(30,53) + rrt(418) * density(50) 
  pd(50,50) = pd(50,50) + rrt(418) * density(53) 
  pd(50,53) = pd(50,53) + rrt(418) * density(50) 
  pd(53,50) = pd(53,50) - rrt(418) * density(53) 
  pd(53,53) = pd(53,53) - rrt(418) * density(50) 
  pd(30,23) = pd(30,23) + rrt(419) * density(53) 
  pd(30,53) = pd(30,53) + rrt(419) * density(23) 
  pd(50,23) = pd(50,23) + rrt(419) * density(53) 
  pd(50,53) = pd(50,53) + rrt(419) * density(23) 
  pd(53,23) = pd(53,23) - rrt(419) * density(53) 
  pd(53,53) = pd(53,53) - rrt(419) * density(23) 
  pd(30,38) = pd(30,38) + rrt(420) * density(53) 
  pd(30,53) = pd(30,53) + rrt(420) * density(38) 
  pd(50,38) = pd(50,38) + rrt(420) * density(53) 
  pd(50,53) = pd(50,53) + rrt(420) * density(38) 
  pd(53,38) = pd(53,38) - rrt(420) * density(53) 
  pd(53,53) = pd(53,53) - rrt(420) * density(38) 
  pd(52,54) = pd(52,54) + rrt(421) 
  pd(53,54) = pd(53,54) + rrt(421) 
  pd(54,54) = pd(54,54) - rrt(421) 
  pd(01,01) = pd(01,01) + rrt(422) * density(23)**2 
  pd(01,23) = pd(01,23) + rrt(422) * density(01) * density(23) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(422) * density(23)**2 * 2.0d0
  pd(23,23) = pd(23,23) - rrt(422) * density(01) * density(23) * 4.0d0
  pd(01,23) = pd(01,23) + rrt(423) * density(23) * density(30) * 2.0d0
  pd(01,30) = pd(01,30) + rrt(423) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(423) * density(23) * density(30) * 4.0d0
  pd(23,30) = pd(23,30) - rrt(423) * density(23)**2 * 2.0d0
  pd(01,23) = pd(01,23) + rrt(424) * density(23) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(424) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(424) * density(23) * density(50) * 4.0d0
  pd(23,50) = pd(23,50) - rrt(424) * density(23)**2 * 2.0d0
  pd(01,23) = pd(01,23) + rrt(425) * density(23)**2 * 3.0d0
  pd(23,23) = pd(23,23) - rrt(425) * density(23)**2 * 6.0d0
  pd(01,23) = pd(01,23) + rrt(426) * density(23) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(426) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(426) * density(23) * density(38) * 4.0d0
  pd(23,38) = pd(23,38) - rrt(426) * density(23)**2 * 2.0d0
  pd(30,01) = pd(30,01) + rrt(427) * density(38)**2 
  pd(30,38) = pd(30,38) + rrt(427) * density(01) * density(38) * 2.0d0
  pd(38,01) = pd(38,01) - rrt(427) * density(38)**2 * 2.0d0
  pd(38,38) = pd(38,38) - rrt(427) * density(01) * density(38) * 4.0d0
  pd(30,30) = pd(30,30) + rrt(428) * density(38)**2 
  pd(30,38) = pd(30,38) + rrt(428) * density(30) * density(38) * 2.0d0
  pd(38,30) = pd(38,30) - rrt(428) * density(38)**2 * 2.0d0
  pd(38,38) = pd(38,38) - rrt(428) * density(30) * density(38) * 4.0d0
  pd(30,23) = pd(30,23) + rrt(429) * density(38)**2 
  pd(30,38) = pd(30,38) + rrt(429) * density(23) * density(38) * 2.0d0
  pd(38,23) = pd(38,23) - rrt(429) * density(38)**2 * 2.0d0
  pd(38,38) = pd(38,38) - rrt(429) * density(23) * density(38) * 4.0d0
  pd(30,38) = pd(30,38) + rrt(430) * density(38)**2 * 3.0d0
  pd(38,38) = pd(38,38) - rrt(430) * density(38)**2 * 6.0d0
  pd(30,38) = pd(30,38) + rrt(431) * density(38) * density(50) * 2.0d0
  pd(30,50) = pd(30,50) + rrt(431) * density(38)**2 
  pd(38,38) = pd(38,38) - rrt(431) * density(38) * density(50) * 4.0d0
  pd(38,50) = pd(38,50) - rrt(431) * density(38)**2 * 2.0d0
  pd(23,01) = pd(23,01) - rrt(432) * density(23) * density(38) 
  pd(23,23) = pd(23,23) - rrt(432) * density(01) * density(38) 
  pd(23,38) = pd(23,38) - rrt(432) * density(01) * density(23) 
  pd(38,01) = pd(38,01) - rrt(432) * density(23) * density(38) 
  pd(38,23) = pd(38,23) - rrt(432) * density(01) * density(38) 
  pd(38,38) = pd(38,38) - rrt(432) * density(01) * density(23) 
  pd(50,01) = pd(50,01) + rrt(432) * density(23) * density(38) 
  pd(50,23) = pd(50,23) + rrt(432) * density(01) * density(38) 
  pd(50,38) = pd(50,38) + rrt(432) * density(01) * density(23) 
  pd(23,23) = pd(23,23) - rrt(433) * density(30) * density(38) 
  pd(23,30) = pd(23,30) - rrt(433) * density(23) * density(38) 
  pd(23,38) = pd(23,38) - rrt(433) * density(23) * density(30) 
  pd(38,23) = pd(38,23) - rrt(433) * density(30) * density(38) 
  pd(38,30) = pd(38,30) - rrt(433) * density(23) * density(38) 
  pd(38,38) = pd(38,38) - rrt(433) * density(23) * density(30) 
  pd(50,23) = pd(50,23) + rrt(433) * density(30) * density(38) 
  pd(50,30) = pd(50,30) + rrt(433) * density(23) * density(38) 
  pd(50,38) = pd(50,38) + rrt(433) * density(23) * density(30) 
  pd(23,23) = pd(23,23) - rrt(434) * density(23) * density(38) * 2.0d0
  pd(23,38) = pd(23,38) - rrt(434) * density(23)**2 
  pd(38,23) = pd(38,23) - rrt(434) * density(23) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(434) * density(23)**2 
  pd(50,23) = pd(50,23) + rrt(434) * density(23) * density(38) * 2.0d0
  pd(50,38) = pd(50,38) + rrt(434) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(435) * density(38)**2 
  pd(23,38) = pd(23,38) - rrt(435) * density(23) * density(38) * 2.0d0
  pd(38,23) = pd(38,23) - rrt(435) * density(38)**2 
  pd(38,38) = pd(38,38) - rrt(435) * density(23) * density(38) * 2.0d0
  pd(50,23) = pd(50,23) + rrt(435) * density(38)**2 
  pd(50,38) = pd(50,38) + rrt(435) * density(23) * density(38) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(436) * density(38) * density(50) 
  pd(23,38) = pd(23,38) - rrt(436) * density(23) * density(50) 
  pd(23,50) = pd(23,50) - rrt(436) * density(23) * density(38) 
  pd(38,23) = pd(38,23) - rrt(436) * density(38) * density(50) 
  pd(38,38) = pd(38,38) - rrt(436) * density(23) * density(50) 
  pd(38,50) = pd(38,50) - rrt(436) * density(23) * density(38) 
  pd(50,23) = pd(50,23) + rrt(436) * density(38) * density(50) 
  pd(50,38) = pd(50,38) + rrt(436) * density(23) * density(50) 
  pd(50,50) = pd(50,50) + rrt(436) * density(23) * density(38) 
  pd(30,01) = pd(30,01) - rrt(437) * density(30) * density(38) 
  pd(30,30) = pd(30,30) - rrt(437) * density(01) * density(38) 
  pd(30,38) = pd(30,38) - rrt(437) * density(01) * density(30) 
  pd(38,01) = pd(38,01) - rrt(437) * density(30) * density(38) 
  pd(38,30) = pd(38,30) - rrt(437) * density(01) * density(38) 
  pd(38,38) = pd(38,38) - rrt(437) * density(01) * density(30) 
  pd(41,01) = pd(41,01) + rrt(437) * density(30) * density(38) 
  pd(41,30) = pd(41,30) + rrt(437) * density(01) * density(38) 
  pd(41,38) = pd(41,38) + rrt(437) * density(01) * density(30) 
  pd(30,30) = pd(30,30) - rrt(438) * density(30) * density(38) * 2.0d0
  pd(30,38) = pd(30,38) - rrt(438) * density(30)**2 
  pd(38,30) = pd(38,30) - rrt(438) * density(30) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(438) * density(30)**2 
  pd(41,30) = pd(41,30) + rrt(438) * density(30) * density(38) * 2.0d0
  pd(41,38) = pd(41,38) + rrt(438) * density(30)**2 
  pd(30,30) = pd(30,30) - rrt(439) * density(38) * density(50) 
  pd(30,38) = pd(30,38) - rrt(439) * density(30) * density(50) 
  pd(30,50) = pd(30,50) - rrt(439) * density(30) * density(38) 
  pd(38,30) = pd(38,30) - rrt(439) * density(38) * density(50) 
  pd(38,38) = pd(38,38) - rrt(439) * density(30) * density(50) 
  pd(38,50) = pd(38,50) - rrt(439) * density(30) * density(38) 
  pd(41,30) = pd(41,30) + rrt(439) * density(38) * density(50) 
  pd(41,38) = pd(41,38) + rrt(439) * density(30) * density(50) 
  pd(41,50) = pd(41,50) + rrt(439) * density(30) * density(38) 
  pd(30,23) = pd(30,23) - rrt(440) * density(30) * density(38) 
  pd(30,30) = pd(30,30) - rrt(440) * density(23) * density(38) 
  pd(30,38) = pd(30,38) - rrt(440) * density(23) * density(30) 
  pd(38,23) = pd(38,23) - rrt(440) * density(30) * density(38) 
  pd(38,30) = pd(38,30) - rrt(440) * density(23) * density(38) 
  pd(38,38) = pd(38,38) - rrt(440) * density(23) * density(30) 
  pd(41,23) = pd(41,23) + rrt(440) * density(30) * density(38) 
  pd(41,30) = pd(41,30) + rrt(440) * density(23) * density(38) 
  pd(41,38) = pd(41,38) + rrt(440) * density(23) * density(30) 
  pd(30,30) = pd(30,30) - rrt(441) * density(38)**2 
  pd(30,38) = pd(30,38) - rrt(441) * density(30) * density(38) * 2.0d0
  pd(38,30) = pd(38,30) - rrt(441) * density(38)**2 
  pd(38,38) = pd(38,38) - rrt(441) * density(30) * density(38) * 2.0d0
  pd(41,30) = pd(41,30) + rrt(441) * density(38)**2 
  pd(41,38) = pd(41,38) + rrt(441) * density(30) * density(38) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(442) * density(38) 
  pd(01,38) = pd(01,38) - rrt(442) * density(01) 
  pd(38,01) = pd(38,01) - rrt(442) * density(38) 
  pd(38,38) = pd(38,38) - rrt(442) * density(01) 
  pd(51,01) = pd(51,01) + rrt(442) * density(38) 
  pd(51,38) = pd(51,38) + rrt(442) * density(01) 
  pd(38,01) = pd(38,01) - rrt(443) * density(38) * density(50) 
  pd(38,38) = pd(38,38) - rrt(443) * density(01) * density(50) 
  pd(38,50) = pd(38,50) - rrt(443) * density(01) * density(38) 
  pd(50,01) = pd(50,01) - rrt(443) * density(38) * density(50) 
  pd(50,38) = pd(50,38) - rrt(443) * density(01) * density(50) 
  pd(50,50) = pd(50,50) - rrt(443) * density(01) * density(38) 
  pd(52,01) = pd(52,01) + rrt(443) * density(38) * density(50) 
  pd(52,38) = pd(52,38) + rrt(443) * density(01) * density(50) 
  pd(52,50) = pd(52,50) + rrt(443) * density(01) * density(38) 
  pd(38,30) = pd(38,30) - rrt(444) * density(38) * density(50) 
  pd(38,38) = pd(38,38) - rrt(444) * density(30) * density(50) 
  pd(38,50) = pd(38,50) - rrt(444) * density(30) * density(38) 
  pd(50,30) = pd(50,30) - rrt(444) * density(38) * density(50) 
  pd(50,38) = pd(50,38) - rrt(444) * density(30) * density(50) 
  pd(50,50) = pd(50,50) - rrt(444) * density(30) * density(38) 
  pd(52,30) = pd(52,30) + rrt(444) * density(38) * density(50) 
  pd(52,38) = pd(52,38) + rrt(444) * density(30) * density(50) 
  pd(52,50) = pd(52,50) + rrt(444) * density(30) * density(38) 
  pd(38,38) = pd(38,38) - rrt(445) * density(50)**2 
  pd(38,50) = pd(38,50) - rrt(445) * density(38) * density(50) * 2.0d0
  pd(50,38) = pd(50,38) - rrt(445) * density(50)**2 
  pd(50,50) = pd(50,50) - rrt(445) * density(38) * density(50) * 2.0d0
  pd(52,38) = pd(52,38) + rrt(445) * density(50)**2 
  pd(52,50) = pd(52,50) + rrt(445) * density(38) * density(50) * 2.0d0
  pd(38,01) = pd(38,01) - rrt(446) * density(38) * density(52) 
  pd(38,38) = pd(38,38) - rrt(446) * density(01) * density(52) 
  pd(38,52) = pd(38,52) - rrt(446) * density(01) * density(38) 
  pd(52,01) = pd(52,01) - rrt(446) * density(38) * density(52) 
  pd(52,38) = pd(52,38) - rrt(446) * density(01) * density(52) 
  pd(52,52) = pd(52,52) - rrt(446) * density(01) * density(38) 
  pd(53,01) = pd(53,01) + rrt(446) * density(38) * density(52) 
  pd(53,38) = pd(53,38) + rrt(446) * density(01) * density(52) 
  pd(53,52) = pd(53,52) + rrt(446) * density(01) * density(38) 
  pd(38,30) = pd(38,30) - rrt(447) * density(38) * density(52) 
  pd(38,38) = pd(38,38) - rrt(447) * density(30) * density(52) 
  pd(38,52) = pd(38,52) - rrt(447) * density(30) * density(38) 
  pd(52,30) = pd(52,30) - rrt(447) * density(38) * density(52) 
  pd(52,38) = pd(52,38) - rrt(447) * density(30) * density(52) 
  pd(52,52) = pd(52,52) - rrt(447) * density(30) * density(38) 
  pd(53,30) = pd(53,30) + rrt(447) * density(38) * density(52) 
  pd(53,38) = pd(53,38) + rrt(447) * density(30) * density(52) 
  pd(53,52) = pd(53,52) + rrt(447) * density(30) * density(38) 
  pd(38,23) = pd(38,23) - rrt(448) * density(38) * density(52) 
  pd(38,38) = pd(38,38) - rrt(448) * density(23) * density(52) 
  pd(38,52) = pd(38,52) - rrt(448) * density(23) * density(38) 
  pd(52,23) = pd(52,23) - rrt(448) * density(38) * density(52) 
  pd(52,38) = pd(52,38) - rrt(448) * density(23) * density(52) 
  pd(52,52) = pd(52,52) - rrt(448) * density(23) * density(38) 
  pd(53,23) = pd(53,23) + rrt(448) * density(38) * density(52) 
  pd(53,38) = pd(53,38) + rrt(448) * density(23) * density(52) 
  pd(53,52) = pd(53,52) + rrt(448) * density(23) * density(38) 
  pd(38,38) = pd(38,38) - rrt(449) * density(38) * density(52) * 2.0d0
  pd(38,52) = pd(38,52) - rrt(449) * density(38)**2 
  pd(52,38) = pd(52,38) - rrt(449) * density(38) * density(52) * 2.0d0
  pd(52,52) = pd(52,52) - rrt(449) * density(38)**2 
  pd(53,38) = pd(53,38) + rrt(449) * density(38) * density(52) * 2.0d0
  pd(53,52) = pd(53,52) + rrt(449) * density(38)**2 
  pd(38,38) = pd(38,38) - rrt(450) * density(50) * density(52) 
  pd(38,50) = pd(38,50) - rrt(450) * density(38) * density(52) 
  pd(38,52) = pd(38,52) - rrt(450) * density(38) * density(50) 
  pd(52,38) = pd(52,38) - rrt(450) * density(50) * density(52) 
  pd(52,50) = pd(52,50) - rrt(450) * density(38) * density(52) 
  pd(52,52) = pd(52,52) - rrt(450) * density(38) * density(50) 
  pd(53,38) = pd(53,38) + rrt(450) * density(50) * density(52) 
  pd(53,50) = pd(53,50) + rrt(450) * density(38) * density(52) 
  pd(53,52) = pd(53,52) + rrt(450) * density(38) * density(50) 
  pd(52,52) = pd(52,52) - rrt(451) * density(53) 
  pd(52,53) = pd(52,53) - rrt(451) * density(52) 
  pd(53,52) = pd(53,52) - rrt(451) * density(53) 
  pd(53,53) = pd(53,53) - rrt(451) * density(52) 
  pd(54,52) = pd(54,52) + rrt(451) * density(53) 
  pd(54,53) = pd(54,53) + rrt(451) * density(52) 
  pd(01,01) = pd(01,01) + rrt(452) * density(23)**2 
  pd(01,23) = pd(01,23) + rrt(452) * density(01) * density(23) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(452) * density(23)**2 * 2.0d0
  pd(23,23) = pd(23,23) - rrt(452) * density(01) * density(23) * 4.0d0
  pd(01,23) = pd(01,23) + rrt(453) * density(23) * density(62) * 2.0d0
  pd(01,62) = pd(01,62) + rrt(453) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(453) * density(23) * density(62) * 4.0d0
  pd(23,62) = pd(23,62) - rrt(453) * density(23)**2 * 2.0d0
  pd(01,23) = pd(01,23) + rrt(454) * density(23) * density(70) * 2.0d0
  pd(01,70) = pd(01,70) + rrt(454) * density(23)**2 
  pd(23,23) = pd(23,23) - rrt(454) * density(23) * density(70) * 4.0d0
  pd(23,70) = pd(23,70) - rrt(454) * density(23)**2 * 2.0d0
  pd(62,62) = pd(62,62) + rrt(455) * density(63)**2 
  pd(62,63) = pd(62,63) + rrt(455) * density(62) * density(63) * 2.0d0
  pd(63,62) = pd(63,62) - rrt(455) * density(63)**2 * 2.0d0
  pd(63,63) = pd(63,63) - rrt(455) * density(62) * density(63) * 4.0d0
  pd(62,01) = pd(62,01) + rrt(456) * density(63)**2 
  pd(62,63) = pd(62,63) + rrt(456) * density(01) * density(63) * 2.0d0
  pd(63,01) = pd(63,01) - rrt(456) * density(63)**2 * 2.0d0
  pd(63,63) = pd(63,63) - rrt(456) * density(01) * density(63) * 4.0d0
  pd(23,01) = pd(23,01) - rrt(457) * density(23) * density(63) 
  pd(23,23) = pd(23,23) - rrt(457) * density(01) * density(63) 
  pd(23,63) = pd(23,63) - rrt(457) * density(01) * density(23) 
  pd(63,01) = pd(63,01) - rrt(457) * density(23) * density(63) 
  pd(63,23) = pd(63,23) - rrt(457) * density(01) * density(63) 
  pd(63,63) = pd(63,63) - rrt(457) * density(01) * density(23) 
  pd(72,01) = pd(72,01) + rrt(457) * density(23) * density(63) 
  pd(72,23) = pd(72,23) + rrt(457) * density(01) * density(63) 
  pd(72,63) = pd(72,63) + rrt(457) * density(01) * density(23) 
  pd(23,23) = pd(23,23) - rrt(458) * density(62) * density(63) 
  pd(23,62) = pd(23,62) - rrt(458) * density(23) * density(63) 
  pd(23,63) = pd(23,63) - rrt(458) * density(23) * density(62) 
  pd(63,23) = pd(63,23) - rrt(458) * density(62) * density(63) 
  pd(63,62) = pd(63,62) - rrt(458) * density(23) * density(63) 
  pd(63,63) = pd(63,63) - rrt(458) * density(23) * density(62) 
  pd(72,23) = pd(72,23) + rrt(458) * density(62) * density(63) 
  pd(72,62) = pd(72,62) + rrt(458) * density(23) * density(63) 
  pd(72,63) = pd(72,63) + rrt(458) * density(23) * density(62) 
  pd(23,23) = pd(23,23) - rrt(459) * density(63) * density(70) 
  pd(23,63) = pd(23,63) - rrt(459) * density(23) * density(70) 
  pd(23,70) = pd(23,70) - rrt(459) * density(23) * density(63) 
  pd(63,23) = pd(63,23) - rrt(459) * density(63) * density(70) 
  pd(63,63) = pd(63,63) - rrt(459) * density(23) * density(70) 
  pd(63,70) = pd(63,70) - rrt(459) * density(23) * density(63) 
  pd(72,23) = pd(72,23) + rrt(459) * density(63) * density(70) 
  pd(72,63) = pd(72,63) + rrt(459) * density(23) * density(70) 
  pd(72,70) = pd(72,70) + rrt(459) * density(23) * density(63) 
  pd(23,01) = pd(23,01) - rrt(460) * density(23) * density(62) 
  pd(23,23) = pd(23,23) - rrt(460) * density(01) * density(62) 
  pd(23,62) = pd(23,62) - rrt(460) * density(01) * density(23) 
  pd(62,01) = pd(62,01) - rrt(460) * density(23) * density(62) 
  pd(62,23) = pd(62,23) - rrt(460) * density(01) * density(62) 
  pd(62,62) = pd(62,62) - rrt(460) * density(01) * density(23) 
  pd(71,01) = pd(71,01) + rrt(460) * density(23) * density(62) 
  pd(71,23) = pd(71,23) + rrt(460) * density(01) * density(62) 
  pd(71,62) = pd(71,62) + rrt(460) * density(01) * density(23) 
  pd(23,23) = pd(23,23) - rrt(461) * density(62)**2 
  pd(23,62) = pd(23,62) - rrt(461) * density(23) * density(62) * 2.0d0
  pd(62,23) = pd(62,23) - rrt(461) * density(62)**2 
  pd(62,62) = pd(62,62) - rrt(461) * density(23) * density(62) * 2.0d0
  pd(71,23) = pd(71,23) + rrt(461) * density(62)**2 
  pd(71,62) = pd(71,62) + rrt(461) * density(23) * density(62) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(462) * density(62) * density(70) 
  pd(23,62) = pd(23,62) - rrt(462) * density(23) * density(70) 
  pd(23,70) = pd(23,70) - rrt(462) * density(23) * density(62) 
  pd(62,23) = pd(62,23) - rrt(462) * density(62) * density(70) 
  pd(62,62) = pd(62,62) - rrt(462) * density(23) * density(70) 
  pd(62,70) = pd(62,70) - rrt(462) * density(23) * density(62) 
  pd(71,23) = pd(71,23) + rrt(462) * density(62) * density(70) 
  pd(71,62) = pd(71,62) + rrt(462) * density(23) * density(70) 
  pd(71,70) = pd(71,70) + rrt(462) * density(23) * density(62) 
  pd(63,01) = pd(63,01) - rrt(463) * density(63) * density(72) 
  pd(63,63) = pd(63,63) - rrt(463) * density(01) * density(72) 
  pd(63,72) = pd(63,72) - rrt(463) * density(01) * density(63) 
  pd(71,01) = pd(71,01) + rrt(463) * density(63) * density(72) 
  pd(71,63) = pd(71,63) + rrt(463) * density(01) * density(72) 
  pd(71,72) = pd(71,72) + rrt(463) * density(01) * density(63) 
  pd(72,01) = pd(72,01) - rrt(463) * density(63) * density(72) 
  pd(72,63) = pd(72,63) - rrt(463) * density(01) * density(72) 
  pd(72,72) = pd(72,72) - rrt(463) * density(01) * density(63) 
  pd(63,62) = pd(63,62) - rrt(464) * density(63) * density(72) 
  pd(63,63) = pd(63,63) - rrt(464) * density(62) * density(72) 
  pd(63,72) = pd(63,72) - rrt(464) * density(62) * density(63) 
  pd(71,62) = pd(71,62) + rrt(464) * density(63) * density(72) 
  pd(71,63) = pd(71,63) + rrt(464) * density(62) * density(72) 
  pd(71,72) = pd(71,72) + rrt(464) * density(62) * density(63) 
  pd(72,62) = pd(72,62) - rrt(464) * density(63) * density(72) 
  pd(72,63) = pd(72,63) - rrt(464) * density(62) * density(72) 
  pd(72,72) = pd(72,72) - rrt(464) * density(62) * density(63) 
  pd(63,63) = pd(63,63) - rrt(465) * density(70) * density(72) 
  pd(63,70) = pd(63,70) - rrt(465) * density(63) * density(72) 
  pd(63,72) = pd(63,72) - rrt(465) * density(63) * density(70) 
  pd(71,63) = pd(71,63) + rrt(465) * density(70) * density(72) 
  pd(71,70) = pd(71,70) + rrt(465) * density(63) * density(72) 
  pd(71,72) = pd(71,72) + rrt(465) * density(63) * density(70) 
  pd(72,63) = pd(72,63) - rrt(465) * density(70) * density(72) 
  pd(72,70) = pd(72,70) - rrt(465) * density(63) * density(72) 
  pd(72,72) = pd(72,72) - rrt(465) * density(63) * density(70) 
  pd(63,01) = pd(63,01) - rrt(466) * density(63) * density(71) 
  pd(63,63) = pd(63,63) - rrt(466) * density(01) * density(71) 
  pd(63,71) = pd(63,71) - rrt(466) * density(01) * density(63) 
  pd(70,01) = pd(70,01) + rrt(466) * density(63) * density(71) 
  pd(70,63) = pd(70,63) + rrt(466) * density(01) * density(71) 
  pd(70,71) = pd(70,71) + rrt(466) * density(01) * density(63) 
  pd(71,01) = pd(71,01) - rrt(466) * density(63) * density(71) 
  pd(71,63) = pd(71,63) - rrt(466) * density(01) * density(71) 
  pd(71,71) = pd(71,71) - rrt(466) * density(01) * density(63) 
  pd(63,62) = pd(63,62) - rrt(467) * density(63) * density(71) 
  pd(63,63) = pd(63,63) - rrt(467) * density(62) * density(71) 
  pd(63,71) = pd(63,71) - rrt(467) * density(62) * density(63) 
  pd(70,62) = pd(70,62) + rrt(467) * density(63) * density(71) 
  pd(70,63) = pd(70,63) + rrt(467) * density(62) * density(71) 
  pd(70,71) = pd(70,71) + rrt(467) * density(62) * density(63) 
  pd(71,62) = pd(71,62) - rrt(467) * density(63) * density(71) 
  pd(71,63) = pd(71,63) - rrt(467) * density(62) * density(71) 
  pd(71,71) = pd(71,71) - rrt(467) * density(62) * density(63) 
  pd(63,63) = pd(63,63) - rrt(468) * density(70) * density(71) 
  pd(63,70) = pd(63,70) - rrt(468) * density(63) * density(71) 
  pd(63,71) = pd(63,71) - rrt(468) * density(63) * density(70) 
  pd(70,63) = pd(70,63) + rrt(468) * density(70) * density(71) 
  pd(70,70) = pd(70,70) + rrt(468) * density(63) * density(71) 
  pd(70,71) = pd(70,71) + rrt(468) * density(63) * density(70) 
  pd(71,63) = pd(71,63) - rrt(468) * density(70) * density(71) 
  pd(71,70) = pd(71,70) - rrt(468) * density(63) * density(71) 
  pd(71,71) = pd(71,71) - rrt(468) * density(63) * density(70) 
  pd(62,01) = pd(62,01) - rrt(469) * density(62) * density(72) 
  pd(62,62) = pd(62,62) - rrt(469) * density(01) * density(72) 
  pd(62,72) = pd(62,72) - rrt(469) * density(01) * density(62) 
  pd(70,01) = pd(70,01) + rrt(469) * density(62) * density(72) 
  pd(70,62) = pd(70,62) + rrt(469) * density(01) * density(72) 
  pd(70,72) = pd(70,72) + rrt(469) * density(01) * density(62) 
  pd(72,01) = pd(72,01) - rrt(469) * density(62) * density(72) 
  pd(72,62) = pd(72,62) - rrt(469) * density(01) * density(72) 
  pd(72,72) = pd(72,72) - rrt(469) * density(01) * density(62) 
  pd(62,62) = pd(62,62) - rrt(470) * density(62) * density(72) * 2.0d0
  pd(62,72) = pd(62,72) - rrt(470) * density(62)**2 
  pd(70,62) = pd(70,62) + rrt(470) * density(62) * density(72) * 2.0d0
  pd(70,72) = pd(70,72) + rrt(470) * density(62)**2 
  pd(72,62) = pd(72,62) - rrt(470) * density(62) * density(72) * 2.0d0
  pd(72,72) = pd(72,72) - rrt(470) * density(62)**2 
  pd(62,62) = pd(62,62) - rrt(471) * density(70) * density(72) 
  pd(62,70) = pd(62,70) - rrt(471) * density(62) * density(72) 
  pd(62,72) = pd(62,72) - rrt(471) * density(62) * density(70) 
  pd(70,62) = pd(70,62) + rrt(471) * density(70) * density(72) 
  pd(70,70) = pd(70,70) + rrt(471) * density(62) * density(72) 
  pd(70,72) = pd(70,72) + rrt(471) * density(62) * density(70) 
  pd(72,62) = pd(72,62) - rrt(471) * density(70) * density(72) 
  pd(72,70) = pd(72,70) - rrt(471) * density(62) * density(72) 
  pd(72,72) = pd(72,72) - rrt(471) * density(62) * density(70) 
  pd(23,26) = pd(23,26) + rrt(472) * density(38) 
  pd(23,38) = pd(23,38) + rrt(472) * density(26) 
  pd(26,26) = pd(26,26) - rrt(472) * density(38) 
  pd(26,38) = pd(26,38) - rrt(472) * density(26) 
  pd(38,26) = pd(38,26) - rrt(472) * density(38) 
  pd(38,38) = pd(38,38) - rrt(472) * density(26) 
  pd(42,26) = pd(42,26) + rrt(472) * density(38) 
  pd(42,38) = pd(42,38) + rrt(472) * density(26) 
  pd(23,26) = pd(23,26) + rrt(473) * density(30) 
  pd(23,30) = pd(23,30) + rrt(473) * density(26) 
  pd(26,26) = pd(26,26) - rrt(473) * density(30) 
  pd(26,30) = pd(26,30) - rrt(473) * density(26) 
  pd(30,26) = pd(30,26) - rrt(473) * density(30) 
  pd(30,30) = pd(30,30) - rrt(473) * density(26) 
  pd(43,26) = pd(43,26) + rrt(473) * density(30) 
  pd(43,30) = pd(43,30) + rrt(473) * density(26) 
  pd(26,26) = pd(26,26) - rrt(474) * density(30) 
  pd(26,30) = pd(26,30) - rrt(474) * density(26) 
  pd(30,26) = pd(30,26) - rrt(474) * density(30) 
  pd(30,30) = pd(30,30) - rrt(474) * density(26) 
  pd(38,26) = pd(38,26) + rrt(474) * density(30) 
  pd(38,30) = pd(38,30) + rrt(474) * density(26) 
  pd(55,26) = pd(55,26) + rrt(474) * density(30) 
  pd(55,30) = pd(55,30) + rrt(474) * density(26) 
  pd(26,26) = pd(26,26) - rrt(475) * density(30) 
  pd(26,30) = pd(26,30) - rrt(475) * density(26) 
  pd(30,26) = pd(30,26) - rrt(475) * density(30) 
  pd(30,30) = pd(30,30) - rrt(475) * density(26) 
  pd(42,26) = pd(42,26) + rrt(475) * density(30) 
  pd(42,30) = pd(42,30) + rrt(475) * density(26) 
  pd(50,26) = pd(50,26) + rrt(475) * density(30) 
  pd(50,30) = pd(50,30) + rrt(475) * density(26) 
  pd(26,26) = pd(26,26) - rrt(476) * density(41) 
  pd(26,41) = pd(26,41) - rrt(476) * density(26) 
  pd(30,26) = pd(30,26) + rrt(476) * density(41) 
  pd(30,41) = pd(30,41) + rrt(476) * density(26) 
  pd(41,26) = pd(41,26) - rrt(476) * density(41) 
  pd(41,41) = pd(41,41) - rrt(476) * density(26) 
  pd(55,26) = pd(55,26) + rrt(476) * density(41) 
  pd(55,41) = pd(55,41) + rrt(476) * density(26) 
  pd(23,26) = pd(23,26) + rrt(477) * density(50) 
  pd(23,50) = pd(23,50) + rrt(477) * density(26) 
  pd(26,26) = pd(26,26) - rrt(477) * density(50) 
  pd(26,50) = pd(26,50) - rrt(477) * density(26) 
  pd(50,26) = pd(50,26) - rrt(477) * density(50) 
  pd(50,50) = pd(50,50) - rrt(477) * density(26) 
  pd(55,26) = pd(55,26) + rrt(477) * density(50) 
  pd(55,50) = pd(55,50) + rrt(477) * density(26) 
  pd(26,26) = pd(26,26) - rrt(478) * density(50) 
  pd(26,50) = pd(26,50) - rrt(478) * density(26) 
  pd(27,26) = pd(27,26) + rrt(478) * density(50) 
  pd(27,50) = pd(27,50) + rrt(478) * density(26) 
  pd(38,26) = pd(38,26) + rrt(478) * density(50) 
  pd(38,50) = pd(38,50) + rrt(478) * density(26) 
  pd(50,26) = pd(50,26) - rrt(478) * density(50) 
  pd(50,50) = pd(50,50) - rrt(478) * density(26) 
  pd(01,26) = pd(01,26) + rrt(479) * density(50) 
  pd(01,50) = pd(01,50) + rrt(479) * density(26) 
  pd(26,26) = pd(26,26) - rrt(479) * density(50) 
  pd(26,50) = pd(26,50) - rrt(479) * density(26) 
  pd(42,26) = pd(42,26) + rrt(479) * density(50) 
  pd(42,50) = pd(42,50) + rrt(479) * density(26) 
  pd(50,26) = pd(50,26) - rrt(479) * density(50) 
  pd(50,50) = pd(50,50) - rrt(479) * density(26) 
  pd(01,26) = pd(01,26) + rrt(480) * density(51) 
  pd(01,51) = pd(01,51) + rrt(480) * density(26) 
  pd(26,26) = pd(26,26) - rrt(480) * density(51) 
  pd(26,51) = pd(26,51) - rrt(480) * density(26) 
  pd(51,26) = pd(51,26) - rrt(480) * density(51) 
  pd(51,51) = pd(51,51) - rrt(480) * density(26) 
  pd(55,26) = pd(55,26) + rrt(480) * density(51) 
  pd(55,51) = pd(55,51) + rrt(480) * density(26) 
  pd(01,01) = pd(01,01) - rrt(481) * density(42) 
  pd(01,42) = pd(01,42) - rrt(481) * density(01) 
  pd(23,01) = pd(23,01) + rrt(481) * density(42) 
  pd(23,42) = pd(23,42) + rrt(481) * density(01) 
  pd(42,01) = pd(42,01) - rrt(481) * density(42) 
  pd(42,42) = pd(42,42) - rrt(481) * density(01) 
  pd(55,01) = pd(55,01) + rrt(481) * density(42) 
  pd(55,42) = pd(55,42) + rrt(481) * density(01) 
  pd(30,30) = pd(30,30) - rrt(482) * density(42) 
  pd(30,42) = pd(30,42) - rrt(482) * density(30) 
  pd(38,30) = pd(38,30) + rrt(482) * density(42) 
  pd(38,42) = pd(38,42) + rrt(482) * density(30) 
  pd(42,30) = pd(42,30) - rrt(482) * density(42) 
  pd(42,42) = pd(42,42) - rrt(482) * density(30) 
  pd(43,30) = pd(43,30) + rrt(482) * density(42) 
  pd(43,42) = pd(43,42) + rrt(482) * density(30) 
  pd(30,41) = pd(30,41) + rrt(483) * density(42) 
  pd(30,42) = pd(30,42) + rrt(483) * density(41) 
  pd(41,41) = pd(41,41) - rrt(483) * density(42) 
  pd(41,42) = pd(41,42) - rrt(483) * density(41) 
  pd(42,41) = pd(42,41) - rrt(483) * density(42) 
  pd(42,42) = pd(42,42) - rrt(483) * density(41) 
  pd(43,41) = pd(43,41) + rrt(483) * density(42) 
  pd(43,42) = pd(43,42) + rrt(483) * density(41) 
  pd(38,42) = pd(38,42) + rrt(484) * density(50) 
  pd(38,50) = pd(38,50) + rrt(484) * density(42) 
  pd(42,42) = pd(42,42) - rrt(484) * density(50) 
  pd(42,50) = pd(42,50) - rrt(484) * density(42) 
  pd(50,42) = pd(50,42) - rrt(484) * density(50) 
  pd(50,50) = pd(50,50) - rrt(484) * density(42) 
  pd(55,42) = pd(55,42) + rrt(484) * density(50) 
  pd(55,50) = pd(55,50) + rrt(484) * density(42) 
  pd(23,42) = pd(23,42) + rrt(485) * density(50) 
  pd(23,50) = pd(23,50) + rrt(485) * density(42) 
  pd(42,42) = pd(42,42) - rrt(485) * density(50) 
  pd(42,50) = pd(42,50) - rrt(485) * density(42) 
  pd(43,42) = pd(43,42) + rrt(485) * density(50) 
  pd(43,50) = pd(43,50) + rrt(485) * density(42) 
  pd(50,42) = pd(50,42) - rrt(485) * density(50) 
  pd(50,50) = pd(50,50) - rrt(485) * density(42) 
  pd(24,24) = pd(24,24) - rrt(486) * density(42) 
  pd(24,42) = pd(24,42) - rrt(486) * density(24) 
  pd(26,24) = pd(26,24) + rrt(486) * density(42) 
  pd(26,42) = pd(26,42) + rrt(486) * density(24) 
  pd(38,24) = pd(38,24) + rrt(486) * density(42) 
  pd(38,42) = pd(38,42) + rrt(486) * density(24) 
  pd(42,24) = pd(42,24) - rrt(486) * density(42) 
  pd(42,42) = pd(42,42) - rrt(486) * density(24) 
  pd(42,42) = pd(42,42) - rrt(487) * density(51) 
  pd(42,51) = pd(42,51) - rrt(487) * density(42) 
  pd(50,42) = pd(50,42) + rrt(487) * density(51) 
  pd(50,51) = pd(50,51) + rrt(487) * density(42) 
  pd(51,42) = pd(51,42) - rrt(487) * density(51) 
  pd(51,51) = pd(51,51) - rrt(487) * density(42) 
  pd(55,42) = pd(55,42) + rrt(487) * density(51) 
  pd(55,51) = pd(55,51) + rrt(487) * density(42) 
  pd(38,42) = pd(38,42) + rrt(488) * density(51) 
  pd(38,51) = pd(38,51) + rrt(488) * density(42) 
  pd(42,42) = pd(42,42) - rrt(488) * density(51) 
  pd(42,51) = pd(42,51) - rrt(488) * density(42) 
  pd(51,42) = pd(51,42) - rrt(488) * density(51) 
  pd(51,51) = pd(51,51) - rrt(488) * density(42) 
  pd(56,42) = pd(56,42) + rrt(488) * density(51) 
  pd(56,51) = pd(56,51) + rrt(488) * density(42) 
  pd(01,42) = pd(01,42) + rrt(489) * density(51) 
  pd(01,51) = pd(01,51) + rrt(489) * density(42) 
  pd(42,42) = pd(42,42) - rrt(489) * density(51) 
  pd(42,51) = pd(42,51) - rrt(489) * density(42) 
  pd(43,42) = pd(43,42) + rrt(489) * density(51) 
  pd(43,51) = pd(43,51) + rrt(489) * density(42) 
  pd(51,42) = pd(51,42) - rrt(489) * density(51) 
  pd(51,51) = pd(51,51) - rrt(489) * density(42) 
  pd(38,42) = pd(38,42) + rrt(490) * density(52) 
  pd(38,52) = pd(38,52) + rrt(490) * density(42) 
  pd(42,42) = pd(42,42) - rrt(490) * density(52) 
  pd(42,52) = pd(42,52) - rrt(490) * density(42) 
  pd(52,42) = pd(52,42) - rrt(490) * density(52) 
  pd(52,52) = pd(52,52) - rrt(490) * density(42) 
  pd(57,42) = pd(57,42) + rrt(490) * density(52) 
  pd(57,52) = pd(57,52) + rrt(490) * density(42) 
  pd(01,27) = pd(01,27) + rrt(491) * density(30) 
  pd(01,30) = pd(01,30) + rrt(491) * density(27) 
  pd(27,27) = pd(27,27) - rrt(491) * density(30) 
  pd(27,30) = pd(27,30) - rrt(491) * density(27) 
  pd(30,27) = pd(30,27) - rrt(491) * density(30) 
  pd(30,30) = pd(30,30) - rrt(491) * density(27) 
  pd(43,27) = pd(43,27) + rrt(491) * density(30) 
  pd(43,30) = pd(43,30) + rrt(491) * density(27) 
  pd(23,27) = pd(23,27) + rrt(492) * density(38) 
  pd(23,38) = pd(23,38) + rrt(492) * density(27) 
  pd(27,27) = pd(27,27) - rrt(492) * density(38) 
  pd(27,38) = pd(27,38) - rrt(492) * density(27) 
  pd(38,27) = pd(38,27) - rrt(492) * density(38) 
  pd(38,38) = pd(38,38) - rrt(492) * density(27) 
  pd(55,27) = pd(55,27) + rrt(492) * density(38) 
  pd(55,38) = pd(55,38) + rrt(492) * density(27) 
  pd(01,27) = pd(01,27) + rrt(493) * density(41) 
  pd(01,41) = pd(01,41) + rrt(493) * density(27) 
  pd(27,27) = pd(27,27) - rrt(493) * density(41) 
  pd(27,41) = pd(27,41) - rrt(493) * density(27) 
  pd(38,27) = pd(38,27) + rrt(493) * density(41) 
  pd(38,41) = pd(38,41) + rrt(493) * density(27) 
  pd(41,27) = pd(41,27) - rrt(493) * density(41) 
  pd(41,41) = pd(41,41) - rrt(493) * density(27) 
  pd(43,27) = pd(43,27) + rrt(493) * density(41) 
  pd(43,41) = pd(43,41) + rrt(493) * density(27) 
  pd(01,23) = pd(01,23) + rrt(494) * density(27) 
  pd(01,27) = pd(01,27) + rrt(494) * density(23) 
  pd(23,23) = pd(23,23) - rrt(494) * density(27) 
  pd(23,27) = pd(23,27) - rrt(494) * density(23) 
  pd(26,23) = pd(26,23) + rrt(494) * density(27) 
  pd(26,27) = pd(26,27) + rrt(494) * density(23) 
  pd(27,23) = pd(27,23) - rrt(494) * density(27) 
  pd(27,27) = pd(27,27) - rrt(494) * density(23) 
  pd(01,27) = pd(01,27) + rrt(495) * density(50) 
  pd(01,50) = pd(01,50) + rrt(495) * density(27) 
  pd(27,27) = pd(27,27) - rrt(495) * density(50) 
  pd(27,50) = pd(27,50) - rrt(495) * density(27) 
  pd(50,27) = pd(50,27) - rrt(495) * density(50) 
  pd(50,50) = pd(50,50) - rrt(495) * density(27) 
  pd(55,27) = pd(55,27) + rrt(495) * density(50) 
  pd(55,50) = pd(55,50) + rrt(495) * density(27) 
  pd(01,27) = pd(01,27) + rrt(496) * density(51) 
  pd(01,51) = pd(01,51) + rrt(496) * density(27) 
  pd(27,27) = pd(27,27) - rrt(496) * density(51) 
  pd(27,51) = pd(27,51) - rrt(496) * density(27) 
  pd(51,27) = pd(51,27) - rrt(496) * density(51) 
  pd(51,51) = pd(51,51) - rrt(496) * density(27) 
  pd(56,27) = pd(56,27) + rrt(496) * density(51) 
  pd(56,51) = pd(56,51) + rrt(496) * density(27) 
  pd(01,27) = pd(01,27) + rrt(497) * density(51) 
  pd(01,51) = pd(01,51) + rrt(497) * density(27) 
  pd(23,27) = pd(23,27) + rrt(497) * density(51) 
  pd(23,51) = pd(23,51) + rrt(497) * density(27) 
  pd(27,27) = pd(27,27) - rrt(497) * density(51) 
  pd(27,51) = pd(27,51) - rrt(497) * density(27) 
  pd(51,27) = pd(51,27) - rrt(497) * density(51) 
  pd(51,51) = pd(51,51) - rrt(497) * density(27) 
  pd(55,27) = pd(55,27) + rrt(497) * density(51) 
  pd(55,51) = pd(55,51) + rrt(497) * density(27) 
  pd(01,01) = pd(01,01) - rrt(498) * density(43) 
  pd(01,43) = pd(01,43) - rrt(498) * density(01) 
  pd(43,01) = pd(43,01) - rrt(498) * density(43) 
  pd(43,43) = pd(43,43) - rrt(498) * density(01) 
  pd(50,01) = pd(50,01) + rrt(498) * density(43) 
  pd(50,43) = pd(50,43) + rrt(498) * density(01) 
  pd(55,01) = pd(55,01) + rrt(498) * density(43) 
  pd(55,43) = pd(55,43) + rrt(498) * density(01) 
  pd(23,23) = pd(23,23) - rrt(499) * density(43) 
  pd(23,43) = pd(23,43) - rrt(499) * density(23) 
  pd(38,23) = pd(38,23) + rrt(499) * density(43) 
  pd(38,43) = pd(38,43) + rrt(499) * density(23) 
  pd(43,23) = pd(43,23) - rrt(499) * density(43) 
  pd(43,43) = pd(43,43) - rrt(499) * density(23) 
  pd(55,23) = pd(55,23) + rrt(499) * density(43) 
  pd(55,43) = pd(55,43) + rrt(499) * density(23) 
  pd(30,43) = pd(30,43) + rrt(500) * density(50) 
  pd(30,50) = pd(30,50) + rrt(500) * density(43) 
  pd(43,43) = pd(43,43) - rrt(500) * density(50) 
  pd(43,50) = pd(43,50) - rrt(500) * density(43) 
  pd(50,43) = pd(50,43) - rrt(500) * density(50) 
  pd(50,50) = pd(50,50) - rrt(500) * density(43) 
  pd(55,43) = pd(55,43) + rrt(500) * density(50) 
  pd(55,50) = pd(55,50) + rrt(500) * density(43) 
  pd(41,43) = pd(41,43) + rrt(501) * density(52) 
  pd(41,52) = pd(41,52) + rrt(501) * density(43) 
  pd(43,43) = pd(43,43) - rrt(501) * density(52) 
  pd(43,52) = pd(43,52) - rrt(501) * density(43) 
  pd(52,43) = pd(52,43) - rrt(501) * density(52) 
  pd(52,52) = pd(52,52) - rrt(501) * density(43) 
  pd(55,43) = pd(55,43) + rrt(501) * density(52) 
  pd(55,52) = pd(55,52) + rrt(501) * density(43) 
  pd(30,43) = pd(30,43) + rrt(502) * density(52) 
  pd(30,52) = pd(30,52) + rrt(502) * density(43) 
  pd(43,43) = pd(43,43) - rrt(502) * density(52) 
  pd(43,52) = pd(43,52) - rrt(502) * density(43) 
  pd(52,43) = pd(52,43) - rrt(502) * density(52) 
  pd(52,52) = pd(52,52) - rrt(502) * density(43) 
  pd(57,43) = pd(57,43) + rrt(502) * density(52) 
  pd(57,52) = pd(57,52) + rrt(502) * density(43) 
  pd(01,28) = pd(01,28) + rrt(503) * density(30) 
  pd(01,30) = pd(01,30) + rrt(503) * density(28) 
  pd(23,28) = pd(23,28) + rrt(503) * density(30) 
  pd(23,30) = pd(23,30) + rrt(503) * density(28) 
  pd(28,28) = pd(28,28) - rrt(503) * density(30) 
  pd(28,30) = pd(28,30) - rrt(503) * density(28) 
  pd(30,28) = pd(30,28) - rrt(503) * density(30) 
  pd(30,30) = pd(30,30) - rrt(503) * density(28) 
  pd(43,28) = pd(43,28) + rrt(503) * density(30) 
  pd(43,30) = pd(43,30) + rrt(503) * density(28) 
  pd(01,28) = pd(01,28) + rrt(504) * density(30) 
  pd(01,30) = pd(01,30) + rrt(504) * density(28) 
  pd(28,28) = pd(28,28) - rrt(504) * density(30) 
  pd(28,30) = pd(28,30) - rrt(504) * density(28) 
  pd(30,28) = pd(30,28) - rrt(504) * density(30) 
  pd(30,30) = pd(30,30) - rrt(504) * density(28) 
  pd(57,28) = pd(57,28) + rrt(504) * density(30) 
  pd(57,30) = pd(57,30) + rrt(504) * density(28) 
  pd(01,23) = pd(01,23) + rrt(505) * density(28) 
  pd(01,28) = pd(01,28) + rrt(505) * density(23) 
  pd(23,23) = pd(23,23) - rrt(505) * density(28) 
  pd(23,28) = pd(23,28) - rrt(505) * density(23) 
  pd(27,23) = pd(27,23) + rrt(505) * density(28) 
  pd(27,28) = pd(27,28) + rrt(505) * density(23) 
  pd(28,23) = pd(28,23) - rrt(505) * density(28) 
  pd(28,28) = pd(28,28) - rrt(505) * density(23) 
  pd(01,28) = pd(01,28) + rrt(506) * density(50) 
  pd(01,50) = pd(01,50) + rrt(506) * density(28) 
  pd(23,28) = pd(23,28) + rrt(506) * density(50) 
  pd(23,50) = pd(23,50) + rrt(506) * density(28) 
  pd(28,28) = pd(28,28) - rrt(506) * density(50) 
  pd(28,50) = pd(28,50) - rrt(506) * density(28) 
  pd(50,28) = pd(50,28) - rrt(506) * density(50) 
  pd(50,50) = pd(50,50) - rrt(506) * density(28) 
  pd(55,28) = pd(55,28) + rrt(506) * density(50) 
  pd(55,50) = pd(55,50) + rrt(506) * density(28) 
  pd(01,28) = pd(01,28) + rrt(507) * density(50) 
  pd(01,50) = pd(01,50) + rrt(507) * density(28) 
  pd(28,28) = pd(28,28) - rrt(507) * density(50) 
  pd(28,50) = pd(28,50) - rrt(507) * density(28) 
  pd(50,28) = pd(50,28) - rrt(507) * density(50) 
  pd(50,50) = pd(50,50) - rrt(507) * density(28) 
  pd(56,28) = pd(56,28) + rrt(507) * density(50) 
  pd(56,50) = pd(56,50) + rrt(507) * density(28) 
  pd(50,50) = pd(50,50) - rrt(508) * density(57) 
  pd(50,57) = pd(50,57) - rrt(508) * density(50) 
  pd(52,50) = pd(52,50) + rrt(508) * density(57) 
  pd(52,57) = pd(52,57) + rrt(508) * density(50) 
  pd(55,50) = pd(55,50) + rrt(508) * density(57) 
  pd(55,57) = pd(55,57) + rrt(508) * density(50) 
  pd(57,50) = pd(57,50) - rrt(508) * density(57) 
  pd(57,57) = pd(57,57) - rrt(508) * density(50) 
  pd(50,50) = pd(50,50) - rrt(509) * density(56) 
  pd(50,56) = pd(50,56) - rrt(509) * density(50) 
  pd(51,50) = pd(51,50) + rrt(509) * density(56) 
  pd(51,56) = pd(51,56) + rrt(509) * density(50) 
  pd(55,50) = pd(55,50) + rrt(509) * density(56) 
  pd(55,56) = pd(55,56) + rrt(509) * density(50) 
  pd(56,50) = pd(56,50) - rrt(509) * density(56) 
  pd(56,56) = pd(56,56) - rrt(509) * density(50) 
  pd(01,01) = pd(01,01) + rrt(510) * density(29) 
  pd(01,29) = pd(01,29) + rrt(510) * density(01) 
  pd(27,01) = pd(27,01) + rrt(510) * density(29) 
  pd(27,29) = pd(27,29) + rrt(510) * density(01) 
  pd(29,01) = pd(29,01) - rrt(510) * density(29) 
  pd(29,29) = pd(29,29) - rrt(510) * density(01) 
  pd(01,29) = pd(01,29) + rrt(511) * density(30) * 2.0d0
  pd(01,30) = pd(01,30) + rrt(511) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(511) * density(30) 
  pd(29,30) = pd(29,30) - rrt(511) * density(29) 
  pd(30,29) = pd(30,29) - rrt(511) * density(30) 
  pd(30,30) = pd(30,30) - rrt(511) * density(29) 
  pd(43,29) = pd(43,29) + rrt(511) * density(30) 
  pd(43,30) = pd(43,30) + rrt(511) * density(29) 
  pd(01,29) = pd(01,29) + rrt(512) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(512) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(512) * density(38) 
  pd(29,38) = pd(29,38) - rrt(512) * density(29) 
  pd(38,29) = pd(38,29) - rrt(512) * density(38) 
  pd(38,38) = pd(38,38) - rrt(512) * density(29) 
  pd(42,29) = pd(42,29) + rrt(512) * density(38) 
  pd(42,38) = pd(42,38) + rrt(512) * density(29) 
  pd(01,23) = pd(01,23) + rrt(513) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(513) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(513) * density(29) 
  pd(23,29) = pd(23,29) - rrt(513) * density(23) 
  pd(26,23) = pd(26,23) + rrt(513) * density(29) 
  pd(26,29) = pd(26,29) + rrt(513) * density(23) 
  pd(29,23) = pd(29,23) - rrt(513) * density(29) 
  pd(29,29) = pd(29,29) - rrt(513) * density(23) 
  pd(01,29) = pd(01,29) + rrt(514) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(514) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(514) * density(50) 
  pd(29,50) = pd(29,50) - rrt(514) * density(29) 
  pd(50,29) = pd(50,29) - rrt(514) * density(50) 
  pd(50,50) = pd(50,50) - rrt(514) * density(29) 
  pd(55,29) = pd(55,29) + rrt(514) * density(50) 
  pd(55,50) = pd(55,50) + rrt(514) * density(29) 
  pd(01,01) = pd(01,01) - rrt(515) * density(45) 
  pd(01,45) = pd(01,45) - rrt(515) * density(01) 
  pd(30,01) = pd(30,01) + rrt(515) * density(45) 
  pd(30,45) = pd(30,45) + rrt(515) * density(01) 
  pd(45,01) = pd(45,01) - rrt(515) * density(45) 
  pd(45,45) = pd(45,45) - rrt(515) * density(01) 
  pd(78,01) = pd(78,01) + rrt(515) * density(45) 
  pd(78,45) = pd(78,45) + rrt(515) * density(01) 
  pd(30,30) = pd(30,30) + rrt(516) * density(45) 
  pd(30,45) = pd(30,45) + rrt(516) * density(30) 
  pd(43,30) = pd(43,30) + rrt(516) * density(45) 
  pd(43,45) = pd(43,45) + rrt(516) * density(30) 
  pd(45,30) = pd(45,30) - rrt(516) * density(45) 
  pd(45,45) = pd(45,45) - rrt(516) * density(30) 
  pd(30,35) = pd(30,35) + rrt(517) * density(45) * 2.0d0
  pd(30,45) = pd(30,45) + rrt(517) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(517) * density(45) 
  pd(35,45) = pd(35,45) - rrt(517) * density(35) 
  pd(43,35) = pd(43,35) + rrt(517) * density(45) 
  pd(43,45) = pd(43,45) + rrt(517) * density(35) 
  pd(45,35) = pd(45,35) - rrt(517) * density(45) 
  pd(45,45) = pd(45,45) - rrt(517) * density(35) 
  pd(30,36) = pd(30,36) + rrt(518) * density(45) * 2.0d0
  pd(30,45) = pd(30,45) + rrt(518) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(518) * density(45) 
  pd(36,45) = pd(36,45) - rrt(518) * density(36) 
  pd(43,36) = pd(43,36) + rrt(518) * density(45) 
  pd(43,45) = pd(43,45) + rrt(518) * density(36) 
  pd(45,36) = pd(45,36) - rrt(518) * density(45) 
  pd(45,45) = pd(45,45) - rrt(518) * density(36) 
  pd(38,38) = pd(38,38) - rrt(519) * density(45) 
  pd(38,45) = pd(38,45) - rrt(519) * density(38) 
  pd(41,38) = pd(41,38) + rrt(519) * density(45) 
  pd(41,45) = pd(41,45) + rrt(519) * density(38) 
  pd(43,38) = pd(43,38) + rrt(519) * density(45) 
  pd(43,45) = pd(43,45) + rrt(519) * density(38) 
  pd(45,38) = pd(45,38) - rrt(519) * density(45) 
  pd(45,45) = pd(45,45) - rrt(519) * density(38) 
  pd(30,45) = pd(30,45) + rrt(520) * density(50) * 2.0d0
  pd(30,50) = pd(30,50) + rrt(520) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(520) * density(50) 
  pd(45,50) = pd(45,50) - rrt(520) * density(45) 
  pd(50,45) = pd(50,45) - rrt(520) * density(50) 
  pd(50,50) = pd(50,50) - rrt(520) * density(45) 
  pd(55,45) = pd(55,45) + rrt(520) * density(50) 
  pd(55,50) = pd(55,50) + rrt(520) * density(45) 
  pd(01,01) = pd(01,01) + rrt(521) * density(78) 
  pd(01,78) = pd(01,78) + rrt(521) * density(01) 
  pd(43,01) = pd(43,01) + rrt(521) * density(78) 
  pd(43,78) = pd(43,78) + rrt(521) * density(01) 
  pd(78,01) = pd(78,01) - rrt(521) * density(78) 
  pd(78,78) = pd(78,78) - rrt(521) * density(01) 
  pd(01,30) = pd(01,30) + rrt(522) * density(78) 
  pd(01,78) = pd(01,78) + rrt(522) * density(30) 
  pd(30,30) = pd(30,30) - rrt(522) * density(78) 
  pd(30,78) = pd(30,78) - rrt(522) * density(30) 
  pd(45,30) = pd(45,30) + rrt(522) * density(78) 
  pd(45,78) = pd(45,78) + rrt(522) * density(30) 
  pd(78,30) = pd(78,30) - rrt(522) * density(78) 
  pd(78,78) = pd(78,78) - rrt(522) * density(30) 
  pd(01,01) = pd(01,01) - rrt(523) * density(01) * density(26) * 2.0d0
  pd(01,26) = pd(01,26) - rrt(523) * density(01)**2 
  pd(26,01) = pd(26,01) - rrt(523) * density(01) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(523) * density(01)**2 
  pd(28,01) = pd(28,01) + rrt(523) * density(01) * density(26) * 2.0d0
  pd(28,26) = pd(28,26) + rrt(523) * density(01)**2 
  pd(26,26) = pd(26,26) - rrt(524) * density(38) 
  pd(26,38) = pd(26,38) - rrt(524) * density(26) 
  pd(38,26) = pd(38,26) - rrt(524) * density(38) 
  pd(38,38) = pd(38,38) - rrt(524) * density(26) 
  pd(55,26) = pd(55,26) + rrt(524) * density(38) 
  pd(55,38) = pd(55,38) + rrt(524) * density(26) 
  pd(23,23) = pd(23,23) - rrt(525) * density(26) 
  pd(23,26) = pd(23,26) - rrt(525) * density(23) 
  pd(26,23) = pd(26,23) - rrt(525) * density(26) 
  pd(26,26) = pd(26,26) - rrt(525) * density(23) 
  pd(27,23) = pd(27,23) + rrt(525) * density(26) 
  pd(27,26) = pd(27,26) + rrt(525) * density(23) 
  pd(01,01) = pd(01,01) - rrt(526) * density(42) 
  pd(01,42) = pd(01,42) - rrt(526) * density(01) 
  pd(23,01) = pd(23,01) + rrt(526) * density(42) 
  pd(23,42) = pd(23,42) + rrt(526) * density(01) 
  pd(42,01) = pd(42,01) - rrt(526) * density(42) 
  pd(42,42) = pd(42,42) - rrt(526) * density(01) 
  pd(55,01) = pd(55,01) + rrt(526) * density(42) 
  pd(55,42) = pd(55,42) + rrt(526) * density(01) 
  pd(38,38) = pd(38,38) - rrt(527) * density(42) 
  pd(38,42) = pd(38,42) - rrt(527) * density(38) 
  pd(42,38) = pd(42,38) - rrt(527) * density(42) 
  pd(42,42) = pd(42,42) - rrt(527) * density(38) 
  pd(43,38) = pd(43,38) + rrt(527) * density(42) 
  pd(43,42) = pd(43,42) + rrt(527) * density(38) 
  pd(23,23) = pd(23,23) - rrt(528) * density(42) 
  pd(23,42) = pd(23,42) - rrt(528) * density(23) 
  pd(42,23) = pd(42,23) - rrt(528) * density(42) 
  pd(42,42) = pd(42,42) - rrt(528) * density(23) 
  pd(55,23) = pd(55,23) + rrt(528) * density(42) 
  pd(55,42) = pd(55,42) + rrt(528) * density(23) 
  pd(01,01) = pd(01,01) - rrt(529) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(529) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(529) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(529) * density(01)**2 
  pd(29,01) = pd(29,01) + rrt(529) * density(01) * density(27) * 2.0d0
  pd(29,27) = pd(29,27) + rrt(529) * density(01)**2 
  pd(23,01) = pd(23,01) - rrt(530) * density(23) * density(27) 
  pd(23,23) = pd(23,23) - rrt(530) * density(01) * density(27) 
  pd(23,27) = pd(23,27) - rrt(530) * density(01) * density(23) 
  pd(27,01) = pd(27,01) - rrt(530) * density(23) * density(27) 
  pd(27,23) = pd(27,23) - rrt(530) * density(01) * density(27) 
  pd(27,27) = pd(27,27) - rrt(530) * density(01) * density(23) 
  pd(28,01) = pd(28,01) + rrt(530) * density(23) * density(27) 
  pd(28,23) = pd(28,23) + rrt(530) * density(01) * density(27) 
  pd(28,27) = pd(28,27) + rrt(530) * density(01) * density(23) 
  pd(30,30) = pd(30,30) - rrt(531) * density(30) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) - rrt(531) * density(30)**2 
  pd(43,30) = pd(43,30) - rrt(531) * density(30) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(531) * density(30)**2 
  pd(45,30) = pd(45,30) + rrt(531) * density(30) * density(43) * 2.0d0
  pd(45,43) = pd(45,43) + rrt(531) * density(30)**2 
  pd(01,01) = pd(01,01) - rrt(532) * density(01) * density(43) * 2.0d0
  pd(01,43) = pd(01,43) - rrt(532) * density(01)**2 
  pd(43,01) = pd(43,01) - rrt(532) * density(01) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(532) * density(01)**2 
  pd(78,01) = pd(78,01) + rrt(532) * density(01) * density(43) * 2.0d0
  pd(78,43) = pd(78,43) + rrt(532) * density(01)**2 
  pd(26,26) = pd(26,26) - rrt(533) * density(62) 
  pd(26,62) = pd(26,62) - rrt(533) * density(26) 
  pd(62,26) = pd(62,26) - rrt(533) * density(62) 
  pd(62,62) = pd(62,62) - rrt(533) * density(26) 
  pd(63,26) = pd(63,26) + rrt(533) * density(62) 
  pd(63,62) = pd(63,62) + rrt(533) * density(26) 
  pd(75,26) = pd(75,26) + rrt(533) * density(62) 
  pd(75,62) = pd(75,62) + rrt(533) * density(26) 
  pd(26,26) = pd(26,26) - rrt(534) * density(70) 
  pd(26,70) = pd(26,70) - rrt(534) * density(26) 
  pd(70,26) = pd(70,26) - rrt(534) * density(70) 
  pd(70,70) = pd(70,70) - rrt(534) * density(26) 
  pd(72,26) = pd(72,26) + rrt(534) * density(70) 
  pd(72,70) = pd(72,70) + rrt(534) * density(26) 
  pd(74,26) = pd(74,26) + rrt(534) * density(70) 
  pd(74,70) = pd(74,70) + rrt(534) * density(26) 
  pd(26,26) = pd(26,26) - rrt(535) * density(70) 
  pd(26,70) = pd(26,70) - rrt(535) * density(26) 
  pd(62,26) = pd(62,26) + rrt(535) * density(70) 
  pd(62,70) = pd(62,70) + rrt(535) * density(26) 
  pd(70,26) = pd(70,26) - rrt(535) * density(70) 
  pd(70,70) = pd(70,70) - rrt(535) * density(26) 
  pd(76,26) = pd(76,26) + rrt(535) * density(70) 
  pd(76,70) = pd(76,70) + rrt(535) * density(26) 
  pd(63,64) = pd(63,64) + rrt(536) * density(70) 
  pd(63,70) = pd(63,70) + rrt(536) * density(64) 
  pd(64,64) = pd(64,64) - rrt(536) * density(70) 
  pd(64,70) = pd(64,70) - rrt(536) * density(64) 
  pd(70,64) = pd(70,64) - rrt(536) * density(70) 
  pd(70,70) = pd(70,70) - rrt(536) * density(64) 
  pd(73,64) = pd(73,64) + rrt(536) * density(70) 
  pd(73,70) = pd(73,70) + rrt(536) * density(64) 
  pd(62,63) = pd(62,63) + rrt(537) * density(65) 
  pd(62,65) = pd(62,65) + rrt(537) * density(63) 
  pd(63,63) = pd(63,63) - rrt(537) * density(65) 
  pd(63,65) = pd(63,65) - rrt(537) * density(63) 
  pd(64,63) = pd(64,63) + rrt(537) * density(65) 
  pd(64,65) = pd(64,65) + rrt(537) * density(63) 
  pd(65,63) = pd(65,63) - rrt(537) * density(65) 
  pd(65,65) = pd(65,65) - rrt(537) * density(63) 
  pd(62,62) = pd(62,62) - rrt(538) * density(65) 
  pd(62,65) = pd(62,65) - rrt(538) * density(62) 
  pd(63,62) = pd(63,62) + rrt(538) * density(65) 
  pd(63,65) = pd(63,65) + rrt(538) * density(62) 
  pd(65,62) = pd(65,62) - rrt(538) * density(65) 
  pd(65,65) = pd(65,65) - rrt(538) * density(62) 
  pd(66,62) = pd(66,62) + rrt(538) * density(65) 
  pd(66,65) = pd(66,65) + rrt(538) * density(62) 
  pd(62,65) = pd(62,65) + rrt(539) * density(70) 
  pd(62,70) = pd(62,70) + rrt(539) * density(65) 
  pd(65,65) = pd(65,65) - rrt(539) * density(70) 
  pd(65,70) = pd(65,70) - rrt(539) * density(65) 
  pd(70,65) = pd(70,65) - rrt(539) * density(70) 
  pd(70,70) = pd(70,70) - rrt(539) * density(65) 
  pd(73,65) = pd(73,65) + rrt(539) * density(70) 
  pd(73,70) = pd(73,70) + rrt(539) * density(65) 
  pd(01,01) = pd(01,01) - rrt(540) * density(65) 
  pd(01,65) = pd(01,65) - rrt(540) * density(01) 
  pd(63,01) = pd(63,01) + rrt(540) * density(65) 
  pd(63,65) = pd(63,65) + rrt(540) * density(01) 
  pd(65,01) = pd(65,01) - rrt(540) * density(65) 
  pd(65,65) = pd(65,65) - rrt(540) * density(01) 
  pd(76,01) = pd(76,01) + rrt(540) * density(65) 
  pd(76,65) = pd(76,65) + rrt(540) * density(01) 
  pd(23,62) = pd(23,62) + rrt(541) * density(75) 
  pd(23,75) = pd(23,75) + rrt(541) * density(62) 
  pd(62,62) = pd(62,62) - rrt(541) * density(75) 
  pd(62,75) = pd(62,75) - rrt(541) * density(62) 
  pd(66,62) = pd(66,62) + rrt(541) * density(75) 
  pd(66,75) = pd(66,75) + rrt(541) * density(62) 
  pd(75,62) = pd(75,62) - rrt(541) * density(75) 
  pd(75,75) = pd(75,75) - rrt(541) * density(62) 
  pd(62,62) = pd(62,62) - rrt(542) * density(75) 
  pd(62,75) = pd(62,75) - rrt(542) * density(62) 
  pd(63,62) = pd(63,62) + rrt(542) * density(75) 
  pd(63,75) = pd(63,75) + rrt(542) * density(62) 
  pd(74,62) = pd(74,62) + rrt(542) * density(75) 
  pd(74,75) = pd(74,75) + rrt(542) * density(62) 
  pd(75,62) = pd(75,62) - rrt(542) * density(75) 
  pd(75,75) = pd(75,75) - rrt(542) * density(62) 
  pd(23,70) = pd(23,70) + rrt(544) * density(75) 
  pd(23,75) = pd(23,75) + rrt(544) * density(70) 
  pd(70,70) = pd(70,70) - rrt(544) * density(75) 
  pd(70,75) = pd(70,75) - rrt(544) * density(70) 
  pd(75,70) = pd(75,70) - rrt(544) * density(75) 
  pd(75,75) = pd(75,75) - rrt(544) * density(70) 
  pd(77,70) = pd(77,70) + rrt(544) * density(75) 
  pd(77,75) = pd(77,75) + rrt(544) * density(70) 
  pd(01,01) = pd(01,01) - rrt(545) * density(75) 
  pd(01,75) = pd(01,75) - rrt(545) * density(01) 
  pd(23,01) = pd(23,01) + rrt(545) * density(75) 
  pd(23,75) = pd(23,75) + rrt(545) * density(01) 
  pd(75,01) = pd(75,01) - rrt(545) * density(75) 
  pd(75,75) = pd(75,75) - rrt(545) * density(01) 
  pd(76,01) = pd(76,01) + rrt(545) * density(75) 
  pd(76,75) = pd(76,75) + rrt(545) * density(01) 
  pd(62,62) = pd(62,62) - rrt(546) * density(74) 
  pd(62,74) = pd(62,74) - rrt(546) * density(62) 
  pd(63,62) = pd(63,62) + rrt(546) * density(74) 
  pd(63,74) = pd(63,74) + rrt(546) * density(62) 
  pd(73,62) = pd(73,62) + rrt(546) * density(74) 
  pd(73,74) = pd(73,74) + rrt(546) * density(62) 
  pd(74,62) = pd(74,62) - rrt(546) * density(74) 
  pd(74,74) = pd(74,74) - rrt(546) * density(62) 
  pd(70,70) = pd(70,70) - rrt(547) * density(74) 
  pd(70,74) = pd(70,74) - rrt(547) * density(70) 
  pd(71,70) = pd(71,70) + rrt(547) * density(74) 
  pd(71,74) = pd(71,74) + rrt(547) * density(70) 
  pd(73,70) = pd(73,70) + rrt(547) * density(74) 
  pd(73,74) = pd(73,74) + rrt(547) * density(70) 
  pd(74,70) = pd(74,70) - rrt(547) * density(74) 
  pd(74,74) = pd(74,74) - rrt(547) * density(70) 
  pd(70,70) = pd(70,70) - rrt(548) * density(74) 
  pd(70,74) = pd(70,74) - rrt(548) * density(70) 
  pd(72,70) = pd(72,70) + rrt(548) * density(74) 
  pd(72,74) = pd(72,74) + rrt(548) * density(70) 
  pd(74,70) = pd(74,70) - rrt(548) * density(74) 
  pd(74,74) = pd(74,74) - rrt(548) * density(70) 
  pd(77,70) = pd(77,70) + rrt(548) * density(74) 
  pd(77,74) = pd(77,74) + rrt(548) * density(70) 
  pd(01,70) = pd(01,70) + rrt(549) * density(76) 
  pd(01,76) = pd(01,76) + rrt(549) * density(70) 
  pd(70,70) = pd(70,70) - rrt(549) * density(76) 
  pd(70,76) = pd(70,76) - rrt(549) * density(70) 
  pd(76,70) = pd(76,70) - rrt(549) * density(76) 
  pd(76,76) = pd(76,76) - rrt(549) * density(70) 
  pd(77,70) = pd(77,70) + rrt(549) * density(76) 
  pd(77,76) = pd(77,76) + rrt(549) * density(70) 
  pd(01,27) = pd(01,27) + rrt(550) * density(70) 
  pd(01,70) = pd(01,70) + rrt(550) * density(27) 
  pd(27,27) = pd(27,27) - rrt(550) * density(70) 
  pd(27,70) = pd(27,70) - rrt(550) * density(27) 
  pd(70,27) = pd(70,27) - rrt(550) * density(70) 
  pd(70,70) = pd(70,70) - rrt(550) * density(27) 
  pd(73,27) = pd(73,27) + rrt(550) * density(70) 
  pd(73,70) = pd(73,70) + rrt(550) * density(27) 
  pd(23,26) = pd(23,26) + rrt(551) * density(70) 
  pd(23,70) = pd(23,70) + rrt(551) * density(26) 
  pd(26,26) = pd(26,26) - rrt(551) * density(70) 
  pd(26,70) = pd(26,70) - rrt(551) * density(26) 
  pd(70,26) = pd(70,26) - rrt(551) * density(70) 
  pd(70,70) = pd(70,70) - rrt(551) * density(26) 
  pd(73,26) = pd(73,26) + rrt(551) * density(70) 
  pd(73,70) = pd(73,70) + rrt(551) * density(26) 
  pd(30,43) = pd(30,43) + rrt(552) * density(70) 
  pd(30,70) = pd(30,70) + rrt(552) * density(43) 
  pd(43,43) = pd(43,43) - rrt(552) * density(70) 
  pd(43,70) = pd(43,70) - rrt(552) * density(43) 
  pd(70,43) = pd(70,43) - rrt(552) * density(70) 
  pd(70,70) = pd(70,70) - rrt(552) * density(43) 
  pd(73,43) = pd(73,43) + rrt(552) * density(70) 
  pd(73,70) = pd(73,70) + rrt(552) * density(43) 
  pd(68,69) = pd(68,69) + rrt(553) * density(70) 
  pd(68,70) = pd(68,70) + rrt(553) * density(69) 
  pd(69,69) = pd(69,69) - rrt(553) * density(70) 
  pd(69,70) = pd(69,70) - rrt(553) * density(69) 
  pd(70,69) = pd(70,69) - rrt(553) * density(70) 
  pd(70,70) = pd(70,70) - rrt(553) * density(69) 
  pd(77,69) = pd(77,69) + rrt(553) * density(70) 
  pd(77,70) = pd(77,70) + rrt(553) * density(69) 
  pd(70,70) = pd(70,70) - rrt(554) * density(73) 
  pd(70,73) = pd(70,73) - rrt(554) * density(70) 
  pd(71,70) = pd(71,70) + rrt(554) * density(73) 
  pd(71,73) = pd(71,73) + rrt(554) * density(70) 
  pd(73,70) = pd(73,70) - rrt(554) * density(73) 
  pd(73,73) = pd(73,73) - rrt(554) * density(70) 
  pd(77,70) = pd(77,70) + rrt(554) * density(73) 
  pd(77,73) = pd(77,73) + rrt(554) * density(70) 
  pd(70,71) = pd(70,71) + rrt(555) * density(73) 
  pd(70,73) = pd(70,73) + rrt(555) * density(71) 
  pd(71,71) = pd(71,71) - rrt(555) * density(73) 
  pd(71,73) = pd(71,73) - rrt(555) * density(71) 
  pd(73,71) = pd(73,71) - rrt(555) * density(73) 
  pd(73,73) = pd(73,73) - rrt(555) * density(71) 
  pd(74,71) = pd(74,71) + rrt(555) * density(73) 
  pd(74,73) = pd(74,73) + rrt(555) * density(71) 
  pd(35,35) = pd(35,35) - rrt(556) * density(46) 
  pd(35,46) = pd(35,46) - rrt(556) * density(35) 
  pd(38,35) = pd(38,35) + rrt(556) * density(46) 
  pd(38,46) = pd(38,46) + rrt(556) * density(35) 
  pd(46,35) = pd(46,35) - rrt(556) * density(46) 
  pd(46,46) = pd(46,46) - rrt(556) * density(35) 
  pd(47,35) = pd(47,35) + rrt(556) * density(46) 
  pd(47,46) = pd(47,46) + rrt(556) * density(35) 
  pd(38,41) = pd(38,41) + rrt(557) * density(46) 
  pd(38,46) = pd(38,46) + rrt(557) * density(41) 
  pd(41,41) = pd(41,41) - rrt(557) * density(46) 
  pd(41,46) = pd(41,46) - rrt(557) * density(41) 
  pd(46,41) = pd(46,41) - rrt(557) * density(46) 
  pd(46,46) = pd(46,46) - rrt(557) * density(41) 
  pd(48,41) = pd(48,41) + rrt(557) * density(46) 
  pd(48,46) = pd(48,46) + rrt(557) * density(41) 
  pd(38,46) = pd(38,46) + rrt(558) * density(52) 
  pd(38,52) = pd(38,52) + rrt(558) * density(46) 
  pd(46,46) = pd(46,46) - rrt(558) * density(52) 
  pd(46,52) = pd(46,52) - rrt(558) * density(46) 
  pd(52,46) = pd(52,46) - rrt(558) * density(52) 
  pd(52,52) = pd(52,52) - rrt(558) * density(46) 
  pd(60,46) = pd(60,46) + rrt(558) * density(52) 
  pd(60,52) = pd(60,52) + rrt(558) * density(46) 
  pd(46,46) = pd(46,46) - rrt(559) * density(51) 
  pd(46,51) = pd(46,51) - rrt(559) * density(46) 
  pd(50,46) = pd(50,46) + rrt(559) * density(51) 
  pd(50,51) = pd(50,51) + rrt(559) * density(46) 
  pd(51,46) = pd(51,46) - rrt(559) * density(51) 
  pd(51,51) = pd(51,51) - rrt(559) * density(46) 
  pd(58,46) = pd(58,46) + rrt(559) * density(51) 
  pd(58,51) = pd(58,51) + rrt(559) * density(46) 
  pd(38,46) = pd(38,46) + rrt(560) * density(51) 
  pd(38,51) = pd(38,51) + rrt(560) * density(46) 
  pd(46,46) = pd(46,46) - rrt(560) * density(51) 
  pd(46,51) = pd(46,51) - rrt(560) * density(46) 
  pd(51,46) = pd(51,46) - rrt(560) * density(51) 
  pd(51,51) = pd(51,51) - rrt(560) * density(46) 
  pd(59,46) = pd(59,46) + rrt(560) * density(51) 
  pd(59,51) = pd(59,51) + rrt(560) * density(46) 
  pd(30,38) = pd(30,38) + rrt(561) * density(47) 
  pd(30,47) = pd(30,47) + rrt(561) * density(38) 
  pd(38,38) = pd(38,38) - rrt(561) * density(47) 
  pd(38,47) = pd(38,47) - rrt(561) * density(38) 
  pd(46,38) = pd(46,38) + rrt(561) * density(47) 
  pd(46,47) = pd(46,47) + rrt(561) * density(38) 
  pd(47,38) = pd(47,38) - rrt(561) * density(47) 
  pd(47,47) = pd(47,47) - rrt(561) * density(38) 
  pd(30,41) = pd(30,41) + rrt(562) * density(47) 
  pd(30,47) = pd(30,47) + rrt(562) * density(41) 
  pd(41,41) = pd(41,41) - rrt(562) * density(47) 
  pd(41,47) = pd(41,47) - rrt(562) * density(41) 
  pd(47,41) = pd(47,41) - rrt(562) * density(47) 
  pd(47,47) = pd(47,47) - rrt(562) * density(41) 
  pd(48,41) = pd(48,41) + rrt(562) * density(47) 
  pd(48,47) = pd(48,47) + rrt(562) * density(41) 
  pd(30,47) = pd(30,47) + rrt(563) * density(52) 
  pd(30,52) = pd(30,52) + rrt(563) * density(47) 
  pd(47,47) = pd(47,47) - rrt(563) * density(52) 
  pd(47,52) = pd(47,52) - rrt(563) * density(47) 
  pd(52,47) = pd(52,47) - rrt(563) * density(52) 
  pd(52,52) = pd(52,52) - rrt(563) * density(47) 
  pd(60,47) = pd(60,47) + rrt(563) * density(52) 
  pd(60,52) = pd(60,52) + rrt(563) * density(47) 
  pd(30,47) = pd(30,47) + rrt(564) * density(53) 
  pd(30,53) = pd(30,53) + rrt(564) * density(47) 
  pd(47,47) = pd(47,47) - rrt(564) * density(53) 
  pd(47,53) = pd(47,53) - rrt(564) * density(47) 
  pd(53,47) = pd(53,47) - rrt(564) * density(53) 
  pd(53,53) = pd(53,53) - rrt(564) * density(47) 
  pd(61,47) = pd(61,47) + rrt(564) * density(53) 
  pd(61,53) = pd(61,53) + rrt(564) * density(47) 
  pd(30,38) = pd(30,38) + rrt(565) * density(48) 
  pd(30,48) = pd(30,48) + rrt(565) * density(38) 
  pd(38,38) = pd(38,38) - rrt(565) * density(48) 
  pd(38,48) = pd(38,48) - rrt(565) * density(38) 
  pd(47,38) = pd(47,38) + rrt(565) * density(48) 
  pd(47,48) = pd(47,48) + rrt(565) * density(38) 
  pd(48,38) = pd(48,38) - rrt(565) * density(48) 
  pd(48,48) = pd(48,48) - rrt(565) * density(38) 
  pd(38,48) = pd(38,48) + rrt(566) * density(50) 
  pd(38,50) = pd(38,50) + rrt(566) * density(48) 
  pd(48,48) = pd(48,48) - rrt(566) * density(50) 
  pd(48,50) = pd(48,50) - rrt(566) * density(48) 
  pd(50,48) = pd(50,48) - rrt(566) * density(50) 
  pd(50,50) = pd(50,50) - rrt(566) * density(48) 
  pd(61,48) = pd(61,48) + rrt(566) * density(50) 
  pd(61,50) = pd(61,50) + rrt(566) * density(48) 
  pd(30,48) = pd(30,48) + rrt(567) * density(50) 
  pd(30,50) = pd(30,50) + rrt(567) * density(48) 
  pd(48,48) = pd(48,48) - rrt(567) * density(50) 
  pd(48,50) = pd(48,50) - rrt(567) * density(48) 
  pd(50,48) = pd(50,48) - rrt(567) * density(50) 
  pd(50,50) = pd(50,50) - rrt(567) * density(48) 
  pd(60,48) = pd(60,48) + rrt(567) * density(50) 
  pd(60,50) = pd(60,50) + rrt(567) * density(48) 
  pd(41,48) = pd(41,48) + rrt(568) * density(52) 
  pd(41,52) = pd(41,52) + rrt(568) * density(48) 
  pd(48,48) = pd(48,48) - rrt(568) * density(52) 
  pd(48,52) = pd(48,52) - rrt(568) * density(48) 
  pd(52,48) = pd(52,48) - rrt(568) * density(52) 
  pd(52,52) = pd(52,52) - rrt(568) * density(48) 
  pd(60,48) = pd(60,48) + rrt(568) * density(52) 
  pd(60,52) = pd(60,52) + rrt(568) * density(48) 
  pd(30,48) = pd(30,48) + rrt(569) * density(52) 
  pd(30,52) = pd(30,52) + rrt(569) * density(48) 
  pd(48,48) = pd(48,48) - rrt(569) * density(52) 
  pd(48,52) = pd(48,52) - rrt(569) * density(48) 
  pd(52,48) = pd(52,48) - rrt(569) * density(52) 
  pd(52,52) = pd(52,52) - rrt(569) * density(48) 
  pd(61,48) = pd(61,48) + rrt(569) * density(52) 
  pd(61,52) = pd(61,52) + rrt(569) * density(48) 
  pd(41,48) = pd(41,48) + rrt(570) * density(53) 
  pd(41,53) = pd(41,53) + rrt(570) * density(48) 
  pd(48,48) = pd(48,48) - rrt(570) * density(53) 
  pd(48,53) = pd(48,53) - rrt(570) * density(48) 
  pd(53,48) = pd(53,48) - rrt(570) * density(53) 
  pd(53,53) = pd(53,53) - rrt(570) * density(48) 
  pd(61,48) = pd(61,48) + rrt(570) * density(53) 
  pd(61,53) = pd(61,53) + rrt(570) * density(48) 
  pd(30,30) = pd(30,30) - rrt(571) * density(58) 
  pd(30,58) = pd(30,58) - rrt(571) * density(30) 
  pd(47,30) = pd(47,30) + rrt(571) * density(58) 
  pd(47,58) = pd(47,58) + rrt(571) * density(30) 
  pd(50,30) = pd(50,30) + rrt(571) * density(58) 
  pd(50,58) = pd(50,58) + rrt(571) * density(30) 
  pd(58,30) = pd(58,30) - rrt(571) * density(58) 
  pd(58,58) = pd(58,58) - rrt(571) * density(30) 
  pd(50,52) = pd(50,52) + rrt(572) * density(58) 
  pd(50,58) = pd(50,58) + rrt(572) * density(52) 
  pd(52,52) = pd(52,52) - rrt(572) * density(58) 
  pd(52,58) = pd(52,58) - rrt(572) * density(52) 
  pd(58,52) = pd(58,52) - rrt(572) * density(58) 
  pd(58,58) = pd(58,58) - rrt(572) * density(52) 
  pd(60,52) = pd(60,52) + rrt(572) * density(58) 
  pd(60,58) = pd(60,58) + rrt(572) * density(52) 
  pd(01,51) = pd(01,51) + rrt(573) * density(58) 
  pd(01,58) = pd(01,58) + rrt(573) * density(51) 
  pd(51,51) = pd(51,51) - rrt(573) * density(58) 
  pd(51,58) = pd(51,58) - rrt(573) * density(51) 
  pd(58,51) = pd(58,51) - rrt(573) * density(58) 
  pd(58,58) = pd(58,58) - rrt(573) * density(51) 
  pd(60,51) = pd(60,51) + rrt(573) * density(58) 
  pd(60,58) = pd(60,58) + rrt(573) * density(51) 
  pd(30,41) = pd(30,41) + rrt(574) * density(60) 
  pd(30,60) = pd(30,60) + rrt(574) * density(41) 
  pd(41,41) = pd(41,41) - rrt(574) * density(60) 
  pd(41,60) = pd(41,60) - rrt(574) * density(41) 
  pd(60,41) = pd(60,41) - rrt(574) * density(60) 
  pd(60,60) = pd(60,60) - rrt(574) * density(41) 
  pd(61,41) = pd(61,41) + rrt(574) * density(60) 
  pd(61,60) = pd(61,60) + rrt(574) * density(41) 
  pd(50,52) = pd(50,52) + rrt(575) * density(60) 
  pd(50,60) = pd(50,60) + rrt(575) * density(52) 
  pd(52,52) = pd(52,52) - rrt(575) * density(60) 
  pd(52,60) = pd(52,60) - rrt(575) * density(52) 
  pd(60,52) = pd(60,52) - rrt(575) * density(60) 
  pd(60,60) = pd(60,60) - rrt(575) * density(52) 
  pd(61,52) = pd(61,52) + rrt(575) * density(60) 
  pd(61,60) = pd(61,60) + rrt(575) * density(52) 
  pd(52,53) = pd(52,53) + rrt(576) * density(60) 
  pd(52,60) = pd(52,60) + rrt(576) * density(53) 
  pd(53,53) = pd(53,53) - rrt(576) * density(60) 
  pd(53,60) = pd(53,60) - rrt(576) * density(53) 
  pd(60,53) = pd(60,53) - rrt(576) * density(60) 
  pd(60,60) = pd(60,60) - rrt(576) * density(53) 
  pd(61,53) = pd(61,53) + rrt(576) * density(60) 
  pd(61,60) = pd(61,60) + rrt(576) * density(53) 
  pd(52,54) = pd(52,54) + rrt(577) * density(60) * 2.0d0
  pd(52,60) = pd(52,60) + rrt(577) * density(54) * 2.0d0
  pd(54,54) = pd(54,54) - rrt(577) * density(60) 
  pd(54,60) = pd(54,60) - rrt(577) * density(54) 
  pd(60,54) = pd(60,54) - rrt(577) * density(60) 
  pd(60,60) = pd(60,60) - rrt(577) * density(54) 
  pd(61,54) = pd(61,54) + rrt(577) * density(60) 
  pd(61,60) = pd(61,60) + rrt(577) * density(54) 
  pd(50,50) = pd(50,50) - rrt(578) * density(61) 
  pd(50,61) = pd(50,61) - rrt(578) * density(50) 
  pd(52,50) = pd(52,50) + rrt(578) * density(61) 
  pd(52,61) = pd(52,61) + rrt(578) * density(50) 
  pd(60,50) = pd(60,50) + rrt(578) * density(61) 
  pd(60,61) = pd(60,61) + rrt(578) * density(50) 
  pd(61,50) = pd(61,50) - rrt(578) * density(61) 
  pd(61,61) = pd(61,61) - rrt(578) * density(50) 
  pd(30,01) = pd(30,01) + rrt(579) * density(49) 
  pd(30,49) = pd(30,49) + rrt(579) * density(01) 
  pd(47,01) = pd(47,01) + rrt(579) * density(49) 
  pd(47,49) = pd(47,49) + rrt(579) * density(01) 
  pd(49,01) = pd(49,01) - rrt(579) * density(49) 
  pd(49,49) = pd(49,49) - rrt(579) * density(01) 
  pd(30,30) = pd(30,30) + rrt(580) * density(49) 
  pd(30,49) = pd(30,49) + rrt(580) * density(30) 
  pd(47,30) = pd(47,30) + rrt(580) * density(49) 
  pd(47,49) = pd(47,49) + rrt(580) * density(30) 
  pd(49,30) = pd(49,30) - rrt(580) * density(49) 
  pd(49,49) = pd(49,49) - rrt(580) * density(30) 
  pd(30,38) = pd(30,38) + rrt(581) * density(49) 
  pd(30,49) = pd(30,49) + rrt(581) * density(38) 
  pd(38,38) = pd(38,38) - rrt(581) * density(49) 
  pd(38,49) = pd(38,49) - rrt(581) * density(38) 
  pd(48,38) = pd(48,38) + rrt(581) * density(49) 
  pd(48,49) = pd(48,49) + rrt(581) * density(38) 
  pd(49,38) = pd(49,38) - rrt(581) * density(49) 
  pd(49,49) = pd(49,49) - rrt(581) * density(38) 
  pd(30,38) = pd(30,38) + rrt(582) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(582) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(582) * density(49) 
  pd(38,49) = pd(38,49) - rrt(582) * density(38) 
  pd(46,38) = pd(46,38) + rrt(582) * density(49) 
  pd(46,49) = pd(46,49) + rrt(582) * density(38) 
  pd(49,38) = pd(49,38) - rrt(582) * density(49) 
  pd(49,49) = pd(49,49) - rrt(582) * density(38) 
  pd(30,35) = pd(30,35) + rrt(583) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(583) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(583) * density(49) 
  pd(35,49) = pd(35,49) - rrt(583) * density(35) 
  pd(47,35) = pd(47,35) + rrt(583) * density(49) 
  pd(47,49) = pd(47,49) + rrt(583) * density(35) 
  pd(49,35) = pd(49,35) - rrt(583) * density(49) 
  pd(49,49) = pd(49,49) - rrt(583) * density(35) 
  pd(30,36) = pd(30,36) + rrt(584) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(584) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(584) * density(49) 
  pd(36,49) = pd(36,49) - rrt(584) * density(36) 
  pd(47,36) = pd(47,36) + rrt(584) * density(49) 
  pd(47,49) = pd(47,49) + rrt(584) * density(36) 
  pd(49,36) = pd(49,36) - rrt(584) * density(49) 
  pd(49,49) = pd(49,49) - rrt(584) * density(36) 
  pd(30,49) = pd(30,49) + rrt(585) * density(50) 
  pd(30,50) = pd(30,50) + rrt(585) * density(49) 
  pd(49,49) = pd(49,49) - rrt(585) * density(50) 
  pd(49,50) = pd(49,50) - rrt(585) * density(49) 
  pd(50,49) = pd(50,49) - rrt(585) * density(50) 
  pd(50,50) = pd(50,50) - rrt(585) * density(49) 
  pd(61,49) = pd(61,49) + rrt(585) * density(50) 
  pd(61,50) = pd(61,50) + rrt(585) * density(49) 
  pd(30,30) = pd(30,30) - rrt(586) * density(46) 
  pd(30,46) = pd(30,46) - rrt(586) * density(30) 
  pd(46,30) = pd(46,30) - rrt(586) * density(46) 
  pd(46,46) = pd(46,46) - rrt(586) * density(30) 
  pd(48,30) = pd(48,30) + rrt(586) * density(46) 
  pd(48,46) = pd(48,46) + rrt(586) * density(30) 
  pd(46,46) = pd(46,46) - rrt(587) * density(50) 
  pd(46,50) = pd(46,50) - rrt(587) * density(46) 
  pd(50,46) = pd(50,46) - rrt(587) * density(50) 
  pd(50,50) = pd(50,50) - rrt(587) * density(46) 
  pd(60,46) = pd(60,46) + rrt(587) * density(50) 
  pd(60,50) = pd(60,50) + rrt(587) * density(46) 
  pd(30,30) = pd(30,30) - rrt(588) * density(47) 
  pd(30,47) = pd(30,47) - rrt(588) * density(30) 
  pd(47,30) = pd(47,30) - rrt(588) * density(47) 
  pd(47,47) = pd(47,47) - rrt(588) * density(30) 
  pd(49,30) = pd(49,30) + rrt(588) * density(47) 
  pd(49,47) = pd(49,47) + rrt(588) * density(30) 
  pd(23,26) = pd(23,26) + rrt(589) * density(46) 
  pd(23,46) = pd(23,46) + rrt(589) * density(26) 
  pd(26,26) = pd(26,26) - rrt(589) * density(46) 
  pd(26,46) = pd(26,46) - rrt(589) * density(26) 
  pd(38,26) = pd(38,26) + rrt(589) * density(46) 
  pd(38,46) = pd(38,46) + rrt(589) * density(26) 
  pd(46,26) = pd(46,26) - rrt(589) * density(46) 
  pd(46,46) = pd(46,46) - rrt(589) * density(26) 
  pd(01,27) = pd(01,27) + rrt(590) * density(46) 
  pd(01,46) = pd(01,46) + rrt(590) * density(27) 
  pd(27,27) = pd(27,27) - rrt(590) * density(46) 
  pd(27,46) = pd(27,46) - rrt(590) * density(27) 
  pd(38,27) = pd(38,27) + rrt(590) * density(46) 
  pd(38,46) = pd(38,46) + rrt(590) * density(27) 
  pd(46,27) = pd(46,27) - rrt(590) * density(46) 
  pd(46,46) = pd(46,46) - rrt(590) * density(27) 
  pd(38,42) = pd(38,42) + rrt(591) * density(46) * 2.0d0
  pd(38,46) = pd(38,46) + rrt(591) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(591) * density(46) 
  pd(42,46) = pd(42,46) - rrt(591) * density(42) 
  pd(46,42) = pd(46,42) - rrt(591) * density(46) 
  pd(46,46) = pd(46,46) - rrt(591) * density(42) 
  pd(30,43) = pd(30,43) + rrt(592) * density(46) 
  pd(30,46) = pd(30,46) + rrt(592) * density(43) 
  pd(38,43) = pd(38,43) + rrt(592) * density(46) 
  pd(38,46) = pd(38,46) + rrt(592) * density(43) 
  pd(43,43) = pd(43,43) - rrt(592) * density(46) 
  pd(43,46) = pd(43,46) - rrt(592) * density(43) 
  pd(46,43) = pd(46,43) - rrt(592) * density(46) 
  pd(46,46) = pd(46,46) - rrt(592) * density(43) 
  pd(38,46) = pd(38,46) + rrt(593) * density(55) 
  pd(38,55) = pd(38,55) + rrt(593) * density(46) 
  pd(46,46) = pd(46,46) - rrt(593) * density(55) 
  pd(46,55) = pd(46,55) - rrt(593) * density(46) 
  pd(50,46) = pd(50,46) + rrt(593) * density(55) 
  pd(50,55) = pd(50,55) + rrt(593) * density(46) 
  pd(55,46) = pd(55,46) - rrt(593) * density(55) 
  pd(55,55) = pd(55,55) - rrt(593) * density(46) 
  pd(38,46) = pd(38,46) + rrt(594) * density(56) 
  pd(38,56) = pd(38,56) + rrt(594) * density(46) 
  pd(46,46) = pd(46,46) - rrt(594) * density(56) 
  pd(46,56) = pd(46,56) - rrt(594) * density(46) 
  pd(51,46) = pd(51,46) + rrt(594) * density(56) 
  pd(51,56) = pd(51,56) + rrt(594) * density(46) 
  pd(56,46) = pd(56,46) - rrt(594) * density(56) 
  pd(56,56) = pd(56,56) - rrt(594) * density(46) 
  pd(38,46) = pd(38,46) + rrt(595) * density(57) 
  pd(38,57) = pd(38,57) + rrt(595) * density(46) 
  pd(46,46) = pd(46,46) - rrt(595) * density(57) 
  pd(46,57) = pd(46,57) - rrt(595) * density(46) 
  pd(52,46) = pd(52,46) + rrt(595) * density(57) 
  pd(52,57) = pd(52,57) + rrt(595) * density(46) 
  pd(57,46) = pd(57,46) - rrt(595) * density(57) 
  pd(57,57) = pd(57,57) - rrt(595) * density(46) 
  pd(23,26) = pd(23,26) + rrt(596) * density(47) 
  pd(23,47) = pd(23,47) + rrt(596) * density(26) 
  pd(26,26) = pd(26,26) - rrt(596) * density(47) 
  pd(26,47) = pd(26,47) - rrt(596) * density(26) 
  pd(30,26) = pd(30,26) + rrt(596) * density(47) 
  pd(30,47) = pd(30,47) + rrt(596) * density(26) 
  pd(47,26) = pd(47,26) - rrt(596) * density(47) 
  pd(47,47) = pd(47,47) - rrt(596) * density(26) 
  pd(01,27) = pd(01,27) + rrt(597) * density(47) 
  pd(01,47) = pd(01,47) + rrt(597) * density(27) 
  pd(27,27) = pd(27,27) - rrt(597) * density(47) 
  pd(27,47) = pd(27,47) - rrt(597) * density(27) 
  pd(30,27) = pd(30,27) + rrt(597) * density(47) 
  pd(30,47) = pd(30,47) + rrt(597) * density(27) 
  pd(47,27) = pd(47,27) - rrt(597) * density(47) 
  pd(47,47) = pd(47,47) - rrt(597) * density(27) 
  pd(30,42) = pd(30,42) + rrt(598) * density(47) 
  pd(30,47) = pd(30,47) + rrt(598) * density(42) 
  pd(38,42) = pd(38,42) + rrt(598) * density(47) 
  pd(38,47) = pd(38,47) + rrt(598) * density(42) 
  pd(42,42) = pd(42,42) - rrt(598) * density(47) 
  pd(42,47) = pd(42,47) - rrt(598) * density(42) 
  pd(47,42) = pd(47,42) - rrt(598) * density(47) 
  pd(47,47) = pd(47,47) - rrt(598) * density(42) 
  pd(30,43) = pd(30,43) + rrt(599) * density(47) * 2.0d0
  pd(30,47) = pd(30,47) + rrt(599) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(599) * density(47) 
  pd(43,47) = pd(43,47) - rrt(599) * density(43) 
  pd(47,43) = pd(47,43) - rrt(599) * density(47) 
  pd(47,47) = pd(47,47) - rrt(599) * density(43) 
  pd(30,47) = pd(30,47) + rrt(600) * density(55) 
  pd(30,55) = pd(30,55) + rrt(600) * density(47) 
  pd(47,47) = pd(47,47) - rrt(600) * density(55) 
  pd(47,55) = pd(47,55) - rrt(600) * density(47) 
  pd(50,47) = pd(50,47) + rrt(600) * density(55) 
  pd(50,55) = pd(50,55) + rrt(600) * density(47) 
  pd(55,47) = pd(55,47) - rrt(600) * density(55) 
  pd(55,55) = pd(55,55) - rrt(600) * density(47) 
  pd(30,47) = pd(30,47) + rrt(601) * density(56) 
  pd(30,56) = pd(30,56) + rrt(601) * density(47) 
  pd(47,47) = pd(47,47) - rrt(601) * density(56) 
  pd(47,56) = pd(47,56) - rrt(601) * density(47) 
  pd(51,47) = pd(51,47) + rrt(601) * density(56) 
  pd(51,56) = pd(51,56) + rrt(601) * density(47) 
  pd(56,47) = pd(56,47) - rrt(601) * density(56) 
  pd(56,56) = pd(56,56) - rrt(601) * density(47) 
  pd(30,47) = pd(30,47) + rrt(602) * density(57) 
  pd(30,57) = pd(30,57) + rrt(602) * density(47) 
  pd(47,47) = pd(47,47) - rrt(602) * density(57) 
  pd(47,57) = pd(47,57) - rrt(602) * density(47) 
  pd(52,47) = pd(52,47) + rrt(602) * density(57) 
  pd(52,57) = pd(52,57) + rrt(602) * density(47) 
  pd(57,47) = pd(57,47) - rrt(602) * density(57) 
  pd(57,57) = pd(57,57) - rrt(602) * density(47) 
  pd(23,26) = pd(23,26) + rrt(603) * density(48) 
  pd(23,48) = pd(23,48) + rrt(603) * density(26) 
  pd(26,26) = pd(26,26) - rrt(603) * density(48) 
  pd(26,48) = pd(26,48) - rrt(603) * density(26) 
  pd(41,26) = pd(41,26) + rrt(603) * density(48) 
  pd(41,48) = pd(41,48) + rrt(603) * density(26) 
  pd(48,26) = pd(48,26) - rrt(603) * density(48) 
  pd(48,48) = pd(48,48) - rrt(603) * density(26) 
  pd(01,27) = pd(01,27) + rrt(604) * density(48) 
  pd(01,48) = pd(01,48) + rrt(604) * density(27) 
  pd(27,27) = pd(27,27) - rrt(604) * density(48) 
  pd(27,48) = pd(27,48) - rrt(604) * density(27) 
  pd(41,27) = pd(41,27) + rrt(604) * density(48) 
  pd(41,48) = pd(41,48) + rrt(604) * density(27) 
  pd(48,27) = pd(48,27) - rrt(604) * density(48) 
  pd(48,48) = pd(48,48) - rrt(604) * density(27) 
  pd(38,42) = pd(38,42) + rrt(605) * density(48) 
  pd(38,48) = pd(38,48) + rrt(605) * density(42) 
  pd(41,42) = pd(41,42) + rrt(605) * density(48) 
  pd(41,48) = pd(41,48) + rrt(605) * density(42) 
  pd(42,42) = pd(42,42) - rrt(605) * density(48) 
  pd(42,48) = pd(42,48) - rrt(605) * density(42) 
  pd(48,42) = pd(48,42) - rrt(605) * density(48) 
  pd(48,48) = pd(48,48) - rrt(605) * density(42) 
  pd(30,43) = pd(30,43) + rrt(606) * density(48) 
  pd(30,48) = pd(30,48) + rrt(606) * density(43) 
  pd(41,43) = pd(41,43) + rrt(606) * density(48) 
  pd(41,48) = pd(41,48) + rrt(606) * density(43) 
  pd(43,43) = pd(43,43) - rrt(606) * density(48) 
  pd(43,48) = pd(43,48) - rrt(606) * density(43) 
  pd(48,43) = pd(48,43) - rrt(606) * density(48) 
  pd(48,48) = pd(48,48) - rrt(606) * density(43) 
  pd(41,48) = pd(41,48) + rrt(607) * density(55) 
  pd(41,55) = pd(41,55) + rrt(607) * density(48) 
  pd(48,48) = pd(48,48) - rrt(607) * density(55) 
  pd(48,55) = pd(48,55) - rrt(607) * density(48) 
  pd(50,48) = pd(50,48) + rrt(607) * density(55) 
  pd(50,55) = pd(50,55) + rrt(607) * density(48) 
  pd(55,48) = pd(55,48) - rrt(607) * density(55) 
  pd(55,55) = pd(55,55) - rrt(607) * density(48) 
  pd(41,48) = pd(41,48) + rrt(608) * density(56) 
  pd(41,56) = pd(41,56) + rrt(608) * density(48) 
  pd(48,48) = pd(48,48) - rrt(608) * density(56) 
  pd(48,56) = pd(48,56) - rrt(608) * density(48) 
  pd(51,48) = pd(51,48) + rrt(608) * density(56) 
  pd(51,56) = pd(51,56) + rrt(608) * density(48) 
  pd(56,48) = pd(56,48) - rrt(608) * density(56) 
  pd(56,56) = pd(56,56) - rrt(608) * density(48) 
  pd(41,48) = pd(41,48) + rrt(609) * density(57) 
  pd(41,57) = pd(41,57) + rrt(609) * density(48) 
  pd(48,48) = pd(48,48) - rrt(609) * density(57) 
  pd(48,57) = pd(48,57) - rrt(609) * density(48) 
  pd(52,48) = pd(52,48) + rrt(609) * density(57) 
  pd(52,57) = pd(52,57) + rrt(609) * density(48) 
  pd(57,48) = pd(57,48) - rrt(609) * density(57) 
  pd(57,57) = pd(57,57) - rrt(609) * density(48) 
  pd(23,26) = pd(23,26) + rrt(610) * density(58) 
  pd(23,58) = pd(23,58) + rrt(610) * density(26) 
  pd(26,26) = pd(26,26) - rrt(610) * density(58) 
  pd(26,58) = pd(26,58) - rrt(610) * density(26) 
  pd(50,26) = pd(50,26) + rrt(610) * density(58) 
  pd(50,58) = pd(50,58) + rrt(610) * density(26) 
  pd(58,26) = pd(58,26) - rrt(610) * density(58) 
  pd(58,58) = pd(58,58) - rrt(610) * density(26) 
  pd(01,27) = pd(01,27) + rrt(611) * density(58) 
  pd(01,58) = pd(01,58) + rrt(611) * density(27) 
  pd(27,27) = pd(27,27) - rrt(611) * density(58) 
  pd(27,58) = pd(27,58) - rrt(611) * density(27) 
  pd(50,27) = pd(50,27) + rrt(611) * density(58) 
  pd(50,58) = pd(50,58) + rrt(611) * density(27) 
  pd(58,27) = pd(58,27) - rrt(611) * density(58) 
  pd(58,58) = pd(58,58) - rrt(611) * density(27) 
  pd(38,42) = pd(38,42) + rrt(612) * density(58) 
  pd(38,58) = pd(38,58) + rrt(612) * density(42) 
  pd(42,42) = pd(42,42) - rrt(612) * density(58) 
  pd(42,58) = pd(42,58) - rrt(612) * density(42) 
  pd(50,42) = pd(50,42) + rrt(612) * density(58) 
  pd(50,58) = pd(50,58) + rrt(612) * density(42) 
  pd(58,42) = pd(58,42) - rrt(612) * density(58) 
  pd(58,58) = pd(58,58) - rrt(612) * density(42) 
  pd(30,43) = pd(30,43) + rrt(613) * density(58) 
  pd(30,58) = pd(30,58) + rrt(613) * density(43) 
  pd(43,43) = pd(43,43) - rrt(613) * density(58) 
  pd(43,58) = pd(43,58) - rrt(613) * density(43) 
  pd(50,43) = pd(50,43) + rrt(613) * density(58) 
  pd(50,58) = pd(50,58) + rrt(613) * density(43) 
  pd(58,43) = pd(58,43) - rrt(613) * density(58) 
  pd(58,58) = pd(58,58) - rrt(613) * density(43) 
  pd(50,55) = pd(50,55) + rrt(614) * density(58) * 2.0d0
  pd(50,58) = pd(50,58) + rrt(614) * density(55) * 2.0d0
  pd(55,55) = pd(55,55) - rrt(614) * density(58) 
  pd(55,58) = pd(55,58) - rrt(614) * density(55) 
  pd(58,55) = pd(58,55) - rrt(614) * density(58) 
  pd(58,58) = pd(58,58) - rrt(614) * density(55) 
  pd(50,56) = pd(50,56) + rrt(615) * density(58) 
  pd(50,58) = pd(50,58) + rrt(615) * density(56) 
  pd(51,56) = pd(51,56) + rrt(615) * density(58) 
  pd(51,58) = pd(51,58) + rrt(615) * density(56) 
  pd(56,56) = pd(56,56) - rrt(615) * density(58) 
  pd(56,58) = pd(56,58) - rrt(615) * density(56) 
  pd(58,56) = pd(58,56) - rrt(615) * density(58) 
  pd(58,58) = pd(58,58) - rrt(615) * density(56) 
  pd(50,57) = pd(50,57) + rrt(616) * density(58) 
  pd(50,58) = pd(50,58) + rrt(616) * density(57) 
  pd(52,57) = pd(52,57) + rrt(616) * density(58) 
  pd(52,58) = pd(52,58) + rrt(616) * density(57) 
  pd(57,57) = pd(57,57) - rrt(616) * density(58) 
  pd(57,58) = pd(57,58) - rrt(616) * density(57) 
  pd(58,57) = pd(58,57) - rrt(616) * density(58) 
  pd(58,58) = pd(58,58) - rrt(616) * density(57) 
  pd(23,26) = pd(23,26) + rrt(617) * density(59) 
  pd(23,59) = pd(23,59) + rrt(617) * density(26) 
  pd(26,26) = pd(26,26) - rrt(617) * density(59) 
  pd(26,59) = pd(26,59) - rrt(617) * density(26) 
  pd(51,26) = pd(51,26) + rrt(617) * density(59) 
  pd(51,59) = pd(51,59) + rrt(617) * density(26) 
  pd(59,26) = pd(59,26) - rrt(617) * density(59) 
  pd(59,59) = pd(59,59) - rrt(617) * density(26) 
  pd(01,27) = pd(01,27) + rrt(618) * density(59) 
  pd(01,59) = pd(01,59) + rrt(618) * density(27) 
  pd(27,27) = pd(27,27) - rrt(618) * density(59) 
  pd(27,59) = pd(27,59) - rrt(618) * density(27) 
  pd(51,27) = pd(51,27) + rrt(618) * density(59) 
  pd(51,59) = pd(51,59) + rrt(618) * density(27) 
  pd(59,27) = pd(59,27) - rrt(618) * density(59) 
  pd(59,59) = pd(59,59) - rrt(618) * density(27) 
  pd(38,42) = pd(38,42) + rrt(619) * density(59) 
  pd(38,59) = pd(38,59) + rrt(619) * density(42) 
  pd(42,42) = pd(42,42) - rrt(619) * density(59) 
  pd(42,59) = pd(42,59) - rrt(619) * density(42) 
  pd(51,42) = pd(51,42) + rrt(619) * density(59) 
  pd(51,59) = pd(51,59) + rrt(619) * density(42) 
  pd(59,42) = pd(59,42) - rrt(619) * density(59) 
  pd(59,59) = pd(59,59) - rrt(619) * density(42) 
  pd(30,43) = pd(30,43) + rrt(620) * density(59) 
  pd(30,59) = pd(30,59) + rrt(620) * density(43) 
  pd(43,43) = pd(43,43) - rrt(620) * density(59) 
  pd(43,59) = pd(43,59) - rrt(620) * density(43) 
  pd(51,43) = pd(51,43) + rrt(620) * density(59) 
  pd(51,59) = pd(51,59) + rrt(620) * density(43) 
  pd(59,43) = pd(59,43) - rrt(620) * density(59) 
  pd(59,59) = pd(59,59) - rrt(620) * density(43) 
  pd(50,55) = pd(50,55) + rrt(621) * density(59) 
  pd(50,59) = pd(50,59) + rrt(621) * density(55) 
  pd(51,55) = pd(51,55) + rrt(621) * density(59) 
  pd(51,59) = pd(51,59) + rrt(621) * density(55) 
  pd(55,55) = pd(55,55) - rrt(621) * density(59) 
  pd(55,59) = pd(55,59) - rrt(621) * density(55) 
  pd(59,55) = pd(59,55) - rrt(621) * density(59) 
  pd(59,59) = pd(59,59) - rrt(621) * density(55) 
  pd(51,56) = pd(51,56) + rrt(622) * density(59) * 2.0d0
  pd(51,59) = pd(51,59) + rrt(622) * density(56) * 2.0d0
  pd(56,56) = pd(56,56) - rrt(622) * density(59) 
  pd(56,59) = pd(56,59) - rrt(622) * density(56) 
  pd(59,56) = pd(59,56) - rrt(622) * density(59) 
  pd(59,59) = pd(59,59) - rrt(622) * density(56) 
  pd(51,57) = pd(51,57) + rrt(623) * density(59) 
  pd(51,59) = pd(51,59) + rrt(623) * density(57) 
  pd(52,57) = pd(52,57) + rrt(623) * density(59) 
  pd(52,59) = pd(52,59) + rrt(623) * density(57) 
  pd(57,57) = pd(57,57) - rrt(623) * density(59) 
  pd(57,59) = pd(57,59) - rrt(623) * density(57) 
  pd(59,57) = pd(59,57) - rrt(623) * density(59) 
  pd(59,59) = pd(59,59) - rrt(623) * density(57) 
  pd(23,26) = pd(23,26) + rrt(624) * density(60) 
  pd(23,60) = pd(23,60) + rrt(624) * density(26) 
  pd(26,26) = pd(26,26) - rrt(624) * density(60) 
  pd(26,60) = pd(26,60) - rrt(624) * density(26) 
  pd(52,26) = pd(52,26) + rrt(624) * density(60) 
  pd(52,60) = pd(52,60) + rrt(624) * density(26) 
  pd(60,26) = pd(60,26) - rrt(624) * density(60) 
  pd(60,60) = pd(60,60) - rrt(624) * density(26) 
  pd(01,27) = pd(01,27) + rrt(625) * density(60) 
  pd(01,60) = pd(01,60) + rrt(625) * density(27) 
  pd(27,27) = pd(27,27) - rrt(625) * density(60) 
  pd(27,60) = pd(27,60) - rrt(625) * density(27) 
  pd(52,27) = pd(52,27) + rrt(625) * density(60) 
  pd(52,60) = pd(52,60) + rrt(625) * density(27) 
  pd(60,27) = pd(60,27) - rrt(625) * density(60) 
  pd(60,60) = pd(60,60) - rrt(625) * density(27) 
  pd(38,42) = pd(38,42) + rrt(626) * density(60) 
  pd(38,60) = pd(38,60) + rrt(626) * density(42) 
  pd(42,42) = pd(42,42) - rrt(626) * density(60) 
  pd(42,60) = pd(42,60) - rrt(626) * density(42) 
  pd(52,42) = pd(52,42) + rrt(626) * density(60) 
  pd(52,60) = pd(52,60) + rrt(626) * density(42) 
  pd(60,42) = pd(60,42) - rrt(626) * density(60) 
  pd(60,60) = pd(60,60) - rrt(626) * density(42) 
  pd(30,43) = pd(30,43) + rrt(627) * density(60) 
  pd(30,60) = pd(30,60) + rrt(627) * density(43) 
  pd(43,43) = pd(43,43) - rrt(627) * density(60) 
  pd(43,60) = pd(43,60) - rrt(627) * density(43) 
  pd(52,43) = pd(52,43) + rrt(627) * density(60) 
  pd(52,60) = pd(52,60) + rrt(627) * density(43) 
  pd(60,43) = pd(60,43) - rrt(627) * density(60) 
  pd(60,60) = pd(60,60) - rrt(627) * density(43) 
  pd(50,55) = pd(50,55) + rrt(628) * density(60) 
  pd(50,60) = pd(50,60) + rrt(628) * density(55) 
  pd(52,55) = pd(52,55) + rrt(628) * density(60) 
  pd(52,60) = pd(52,60) + rrt(628) * density(55) 
  pd(55,55) = pd(55,55) - rrt(628) * density(60) 
  pd(55,60) = pd(55,60) - rrt(628) * density(55) 
  pd(60,55) = pd(60,55) - rrt(628) * density(60) 
  pd(60,60) = pd(60,60) - rrt(628) * density(55) 
  pd(51,56) = pd(51,56) + rrt(629) * density(60) 
  pd(51,60) = pd(51,60) + rrt(629) * density(56) 
  pd(52,56) = pd(52,56) + rrt(629) * density(60) 
  pd(52,60) = pd(52,60) + rrt(629) * density(56) 
  pd(56,56) = pd(56,56) - rrt(629) * density(60) 
  pd(56,60) = pd(56,60) - rrt(629) * density(56) 
  pd(60,56) = pd(60,56) - rrt(629) * density(60) 
  pd(60,60) = pd(60,60) - rrt(629) * density(56) 
  pd(52,57) = pd(52,57) + rrt(630) * density(60) * 2.0d0
  pd(52,60) = pd(52,60) + rrt(630) * density(57) * 2.0d0
  pd(57,57) = pd(57,57) - rrt(630) * density(60) 
  pd(57,60) = pd(57,60) - rrt(630) * density(57) 
  pd(60,57) = pd(60,57) - rrt(630) * density(60) 
  pd(60,60) = pd(60,60) - rrt(630) * density(57) 
  pd(23,26) = pd(23,26) + rrt(631) * density(61) 
  pd(23,61) = pd(23,61) + rrt(631) * density(26) 
  pd(26,26) = pd(26,26) - rrt(631) * density(61) 
  pd(26,61) = pd(26,61) - rrt(631) * density(26) 
  pd(53,26) = pd(53,26) + rrt(631) * density(61) 
  pd(53,61) = pd(53,61) + rrt(631) * density(26) 
  pd(61,26) = pd(61,26) - rrt(631) * density(61) 
  pd(61,61) = pd(61,61) - rrt(631) * density(26) 
  pd(01,27) = pd(01,27) + rrt(632) * density(61) 
  pd(01,61) = pd(01,61) + rrt(632) * density(27) 
  pd(27,27) = pd(27,27) - rrt(632) * density(61) 
  pd(27,61) = pd(27,61) - rrt(632) * density(27) 
  pd(53,27) = pd(53,27) + rrt(632) * density(61) 
  pd(53,61) = pd(53,61) + rrt(632) * density(27) 
  pd(61,27) = pd(61,27) - rrt(632) * density(61) 
  pd(61,61) = pd(61,61) - rrt(632) * density(27) 
  pd(38,42) = pd(38,42) + rrt(633) * density(61) 
  pd(38,61) = pd(38,61) + rrt(633) * density(42) 
  pd(42,42) = pd(42,42) - rrt(633) * density(61) 
  pd(42,61) = pd(42,61) - rrt(633) * density(42) 
  pd(53,42) = pd(53,42) + rrt(633) * density(61) 
  pd(53,61) = pd(53,61) + rrt(633) * density(42) 
  pd(61,42) = pd(61,42) - rrt(633) * density(61) 
  pd(61,61) = pd(61,61) - rrt(633) * density(42) 
  pd(30,43) = pd(30,43) + rrt(634) * density(61) 
  pd(30,61) = pd(30,61) + rrt(634) * density(43) 
  pd(43,43) = pd(43,43) - rrt(634) * density(61) 
  pd(43,61) = pd(43,61) - rrt(634) * density(43) 
  pd(53,43) = pd(53,43) + rrt(634) * density(61) 
  pd(53,61) = pd(53,61) + rrt(634) * density(43) 
  pd(61,43) = pd(61,43) - rrt(634) * density(61) 
  pd(61,61) = pd(61,61) - rrt(634) * density(43) 
  pd(50,55) = pd(50,55) + rrt(635) * density(61) 
  pd(50,61) = pd(50,61) + rrt(635) * density(55) 
  pd(53,55) = pd(53,55) + rrt(635) * density(61) 
  pd(53,61) = pd(53,61) + rrt(635) * density(55) 
  pd(55,55) = pd(55,55) - rrt(635) * density(61) 
  pd(55,61) = pd(55,61) - rrt(635) * density(55) 
  pd(61,55) = pd(61,55) - rrt(635) * density(61) 
  pd(61,61) = pd(61,61) - rrt(635) * density(55) 
  pd(51,56) = pd(51,56) + rrt(636) * density(61) 
  pd(51,61) = pd(51,61) + rrt(636) * density(56) 
  pd(53,56) = pd(53,56) + rrt(636) * density(61) 
  pd(53,61) = pd(53,61) + rrt(636) * density(56) 
  pd(56,56) = pd(56,56) - rrt(636) * density(61) 
  pd(56,61) = pd(56,61) - rrt(636) * density(56) 
  pd(61,56) = pd(61,56) - rrt(636) * density(61) 
  pd(61,61) = pd(61,61) - rrt(636) * density(56) 
  pd(52,57) = pd(52,57) + rrt(637) * density(61) 
  pd(52,61) = pd(52,61) + rrt(637) * density(57) 
  pd(53,57) = pd(53,57) + rrt(637) * density(61) 
  pd(53,61) = pd(53,61) + rrt(637) * density(57) 
  pd(57,57) = pd(57,57) - rrt(637) * density(61) 
  pd(57,61) = pd(57,61) - rrt(637) * density(57) 
  pd(61,57) = pd(61,57) - rrt(637) * density(61) 
  pd(61,61) = pd(61,61) - rrt(637) * density(57) 
  pd(23,27) = pd(23,27) + rrt(638) * density(46) * 2.0d0
  pd(23,46) = pd(23,46) + rrt(638) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(638) * density(46) 
  pd(27,46) = pd(27,46) - rrt(638) * density(27) 
  pd(38,27) = pd(38,27) + rrt(638) * density(46) 
  pd(38,46) = pd(38,46) + rrt(638) * density(27) 
  pd(46,27) = pd(46,27) - rrt(638) * density(46) 
  pd(46,46) = pd(46,46) - rrt(638) * density(27) 
  pd(01,28) = pd(01,28) + rrt(639) * density(46) 
  pd(01,46) = pd(01,46) + rrt(639) * density(28) 
  pd(23,28) = pd(23,28) + rrt(639) * density(46) 
  pd(23,46) = pd(23,46) + rrt(639) * density(28) 
  pd(28,28) = pd(28,28) - rrt(639) * density(46) 
  pd(28,46) = pd(28,46) - rrt(639) * density(28) 
  pd(38,28) = pd(38,28) + rrt(639) * density(46) 
  pd(38,46) = pd(38,46) + rrt(639) * density(28) 
  pd(46,28) = pd(46,28) - rrt(639) * density(46) 
  pd(46,46) = pd(46,46) - rrt(639) * density(28) 
  pd(01,29) = pd(01,29) + rrt(640) * density(46) * 2.0d0
  pd(01,46) = pd(01,46) + rrt(640) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(640) * density(46) 
  pd(29,46) = pd(29,46) - rrt(640) * density(29) 
  pd(38,29) = pd(38,29) + rrt(640) * density(46) 
  pd(38,46) = pd(38,46) + rrt(640) * density(29) 
  pd(46,29) = pd(46,29) - rrt(640) * density(46) 
  pd(46,46) = pd(46,46) - rrt(640) * density(29) 
  pd(38,43) = pd(38,43) + rrt(641) * density(46) * 3.0d0
  pd(38,46) = pd(38,46) + rrt(641) * density(43) * 3.0d0
  pd(43,43) = pd(43,43) - rrt(641) * density(46) 
  pd(43,46) = pd(43,46) - rrt(641) * density(43) 
  pd(46,43) = pd(46,43) - rrt(641) * density(46) 
  pd(46,46) = pd(46,46) - rrt(641) * density(43) 
  pd(30,45) = pd(30,45) + rrt(642) * density(46) * 2.0d0
  pd(30,46) = pd(30,46) + rrt(642) * density(45) * 2.0d0
  pd(38,45) = pd(38,45) + rrt(642) * density(46) 
  pd(38,46) = pd(38,46) + rrt(642) * density(45) 
  pd(45,45) = pd(45,45) - rrt(642) * density(46) 
  pd(45,46) = pd(45,46) - rrt(642) * density(45) 
  pd(46,45) = pd(46,45) - rrt(642) * density(46) 
  pd(46,46) = pd(46,46) - rrt(642) * density(45) 
  pd(23,46) = pd(23,46) + rrt(643) * density(55) 
  pd(23,55) = pd(23,55) + rrt(643) * density(46) 
  pd(38,46) = pd(38,46) + rrt(643) * density(55) * 2.0d0
  pd(38,55) = pd(38,55) + rrt(643) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(643) * density(55) 
  pd(46,55) = pd(46,55) - rrt(643) * density(46) 
  pd(55,46) = pd(55,46) - rrt(643) * density(55) 
  pd(55,55) = pd(55,55) - rrt(643) * density(46) 
  pd(01,46) = pd(01,46) + rrt(644) * density(56) 
  pd(01,56) = pd(01,56) + rrt(644) * density(46) 
  pd(38,46) = pd(38,46) + rrt(644) * density(56) * 2.0d0
  pd(38,56) = pd(38,56) + rrt(644) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(644) * density(56) 
  pd(46,56) = pd(46,56) - rrt(644) * density(46) 
  pd(56,46) = pd(56,46) - rrt(644) * density(56) 
  pd(56,56) = pd(56,56) - rrt(644) * density(46) 
  pd(23,46) = pd(23,46) + rrt(645) * density(57) 
  pd(23,57) = pd(23,57) + rrt(645) * density(46) 
  pd(30,46) = pd(30,46) + rrt(645) * density(57) 
  pd(30,57) = pd(30,57) + rrt(645) * density(46) 
  pd(38,46) = pd(38,46) + rrt(645) * density(57) 
  pd(38,57) = pd(38,57) + rrt(645) * density(46) 
  pd(46,46) = pd(46,46) - rrt(645) * density(57) 
  pd(46,57) = pd(46,57) - rrt(645) * density(46) 
  pd(57,46) = pd(57,46) - rrt(645) * density(57) 
  pd(57,57) = pd(57,57) - rrt(645) * density(46) 
  pd(01,46) = pd(01,46) + rrt(646) * density(78) 
  pd(01,78) = pd(01,78) + rrt(646) * density(46) 
  pd(30,46) = pd(30,46) + rrt(646) * density(78) 
  pd(30,78) = pd(30,78) + rrt(646) * density(46) 
  pd(38,46) = pd(38,46) + rrt(646) * density(78) 
  pd(38,78) = pd(38,78) + rrt(646) * density(46) 
  pd(46,46) = pd(46,46) - rrt(646) * density(78) 
  pd(46,78) = pd(46,78) - rrt(646) * density(46) 
  pd(78,46) = pd(78,46) - rrt(646) * density(78) 
  pd(78,78) = pd(78,78) - rrt(646) * density(46) 
  pd(23,27) = pd(23,27) + rrt(647) * density(47) * 2.0d0
  pd(23,47) = pd(23,47) + rrt(647) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(647) * density(47) 
  pd(27,47) = pd(27,47) - rrt(647) * density(27) 
  pd(30,27) = pd(30,27) + rrt(647) * density(47) 
  pd(30,47) = pd(30,47) + rrt(647) * density(27) 
  pd(47,27) = pd(47,27) - rrt(647) * density(47) 
  pd(47,47) = pd(47,47) - rrt(647) * density(27) 
  pd(01,28) = pd(01,28) + rrt(648) * density(47) 
  pd(01,47) = pd(01,47) + rrt(648) * density(28) 
  pd(23,28) = pd(23,28) + rrt(648) * density(47) 
  pd(23,47) = pd(23,47) + rrt(648) * density(28) 
  pd(28,28) = pd(28,28) - rrt(648) * density(47) 
  pd(28,47) = pd(28,47) - rrt(648) * density(28) 
  pd(30,28) = pd(30,28) + rrt(648) * density(47) 
  pd(30,47) = pd(30,47) + rrt(648) * density(28) 
  pd(47,28) = pd(47,28) - rrt(648) * density(47) 
  pd(47,47) = pd(47,47) - rrt(648) * density(28) 
  pd(01,29) = pd(01,29) + rrt(649) * density(47) * 2.0d0
  pd(01,47) = pd(01,47) + rrt(649) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(649) * density(47) 
  pd(29,47) = pd(29,47) - rrt(649) * density(29) 
  pd(30,29) = pd(30,29) + rrt(649) * density(47) 
  pd(30,47) = pd(30,47) + rrt(649) * density(29) 
  pd(47,29) = pd(47,29) - rrt(649) * density(47) 
  pd(47,47) = pd(47,47) - rrt(649) * density(29) 
  pd(30,43) = pd(30,43) + rrt(650) * density(47) 
  pd(30,47) = pd(30,47) + rrt(650) * density(43) 
  pd(38,43) = pd(38,43) + rrt(650) * density(47) * 2.0d0
  pd(38,47) = pd(38,47) + rrt(650) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(650) * density(47) 
  pd(43,47) = pd(43,47) - rrt(650) * density(43) 
  pd(47,43) = pd(47,43) - rrt(650) * density(47) 
  pd(47,47) = pd(47,47) - rrt(650) * density(43) 
  pd(30,45) = pd(30,45) + rrt(651) * density(47) * 3.0d0
  pd(30,47) = pd(30,47) + rrt(651) * density(45) * 3.0d0
  pd(45,45) = pd(45,45) - rrt(651) * density(47) 
  pd(45,47) = pd(45,47) - rrt(651) * density(45) 
  pd(47,45) = pd(47,45) - rrt(651) * density(47) 
  pd(47,47) = pd(47,47) - rrt(651) * density(45) 
  pd(23,47) = pd(23,47) + rrt(652) * density(55) 
  pd(23,55) = pd(23,55) + rrt(652) * density(47) 
  pd(30,47) = pd(30,47) + rrt(652) * density(55) 
  pd(30,55) = pd(30,55) + rrt(652) * density(47) 
  pd(38,47) = pd(38,47) + rrt(652) * density(55) 
  pd(38,55) = pd(38,55) + rrt(652) * density(47) 
  pd(47,47) = pd(47,47) - rrt(652) * density(55) 
  pd(47,55) = pd(47,55) - rrt(652) * density(47) 
  pd(55,47) = pd(55,47) - rrt(652) * density(55) 
  pd(55,55) = pd(55,55) - rrt(652) * density(47) 
  pd(01,47) = pd(01,47) + rrt(653) * density(56) 
  pd(01,56) = pd(01,56) + rrt(653) * density(47) 
  pd(30,47) = pd(30,47) + rrt(653) * density(56) 
  pd(30,56) = pd(30,56) + rrt(653) * density(47) 
  pd(38,47) = pd(38,47) + rrt(653) * density(56) 
  pd(38,56) = pd(38,56) + rrt(653) * density(47) 
  pd(47,47) = pd(47,47) - rrt(653) * density(56) 
  pd(47,56) = pd(47,56) - rrt(653) * density(47) 
  pd(56,47) = pd(56,47) - rrt(653) * density(56) 
  pd(56,56) = pd(56,56) - rrt(653) * density(47) 
  pd(23,47) = pd(23,47) + rrt(654) * density(57) 
  pd(23,57) = pd(23,57) + rrt(654) * density(47) 
  pd(30,47) = pd(30,47) + rrt(654) * density(57) * 2.0d0
  pd(30,57) = pd(30,57) + rrt(654) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(654) * density(57) 
  pd(47,57) = pd(47,57) - rrt(654) * density(47) 
  pd(57,47) = pd(57,47) - rrt(654) * density(57) 
  pd(57,57) = pd(57,57) - rrt(654) * density(47) 
  pd(01,47) = pd(01,47) + rrt(655) * density(78) 
  pd(01,78) = pd(01,78) + rrt(655) * density(47) 
  pd(30,47) = pd(30,47) + rrt(655) * density(78) * 2.0d0
  pd(30,78) = pd(30,78) + rrt(655) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(655) * density(78) 
  pd(47,78) = pd(47,78) - rrt(655) * density(47) 
  pd(78,47) = pd(78,47) - rrt(655) * density(78) 
  pd(78,78) = pd(78,78) - rrt(655) * density(47) 
  pd(23,27) = pd(23,27) + rrt(656) * density(48) * 2.0d0
  pd(23,48) = pd(23,48) + rrt(656) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(656) * density(48) 
  pd(27,48) = pd(27,48) - rrt(656) * density(27) 
  pd(41,27) = pd(41,27) + rrt(656) * density(48) 
  pd(41,48) = pd(41,48) + rrt(656) * density(27) 
  pd(48,27) = pd(48,27) - rrt(656) * density(48) 
  pd(48,48) = pd(48,48) - rrt(656) * density(27) 
  pd(01,28) = pd(01,28) + rrt(657) * density(48) 
  pd(01,48) = pd(01,48) + rrt(657) * density(28) 
  pd(23,28) = pd(23,28) + rrt(657) * density(48) 
  pd(23,48) = pd(23,48) + rrt(657) * density(28) 
  pd(28,28) = pd(28,28) - rrt(657) * density(48) 
  pd(28,48) = pd(28,48) - rrt(657) * density(28) 
  pd(41,28) = pd(41,28) + rrt(657) * density(48) 
  pd(41,48) = pd(41,48) + rrt(657) * density(28) 
  pd(48,28) = pd(48,28) - rrt(657) * density(48) 
  pd(48,48) = pd(48,48) - rrt(657) * density(28) 
  pd(01,29) = pd(01,29) + rrt(658) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) + rrt(658) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(658) * density(48) 
  pd(29,48) = pd(29,48) - rrt(658) * density(29) 
  pd(41,29) = pd(41,29) + rrt(658) * density(48) 
  pd(41,48) = pd(41,48) + rrt(658) * density(29) 
  pd(48,29) = pd(48,29) - rrt(658) * density(48) 
  pd(48,48) = pd(48,48) - rrt(658) * density(29) 
  pd(38,43) = pd(38,43) + rrt(659) * density(48) * 2.0d0
  pd(38,48) = pd(38,48) + rrt(659) * density(43) * 2.0d0
  pd(41,43) = pd(41,43) + rrt(659) * density(48) 
  pd(41,48) = pd(41,48) + rrt(659) * density(43) 
  pd(43,43) = pd(43,43) - rrt(659) * density(48) 
  pd(43,48) = pd(43,48) - rrt(659) * density(43) 
  pd(48,43) = pd(48,43) - rrt(659) * density(48) 
  pd(48,48) = pd(48,48) - rrt(659) * density(43) 
  pd(30,45) = pd(30,45) + rrt(660) * density(48) * 2.0d0
  pd(30,48) = pd(30,48) + rrt(660) * density(45) * 2.0d0
  pd(41,45) = pd(41,45) + rrt(660) * density(48) 
  pd(41,48) = pd(41,48) + rrt(660) * density(45) 
  pd(45,45) = pd(45,45) - rrt(660) * density(48) 
  pd(45,48) = pd(45,48) - rrt(660) * density(45) 
  pd(48,45) = pd(48,45) - rrt(660) * density(48) 
  pd(48,48) = pd(48,48) - rrt(660) * density(45) 
  pd(23,48) = pd(23,48) + rrt(661) * density(55) 
  pd(23,55) = pd(23,55) + rrt(661) * density(48) 
  pd(38,48) = pd(38,48) + rrt(661) * density(55) 
  pd(38,55) = pd(38,55) + rrt(661) * density(48) 
  pd(41,48) = pd(41,48) + rrt(661) * density(55) 
  pd(41,55) = pd(41,55) + rrt(661) * density(48) 
  pd(48,48) = pd(48,48) - rrt(661) * density(55) 
  pd(48,55) = pd(48,55) - rrt(661) * density(48) 
  pd(55,48) = pd(55,48) - rrt(661) * density(55) 
  pd(55,55) = pd(55,55) - rrt(661) * density(48) 
  pd(01,48) = pd(01,48) + rrt(662) * density(56) 
  pd(01,56) = pd(01,56) + rrt(662) * density(48) 
  pd(38,48) = pd(38,48) + rrt(662) * density(56) 
  pd(38,56) = pd(38,56) + rrt(662) * density(48) 
  pd(41,48) = pd(41,48) + rrt(662) * density(56) 
  pd(41,56) = pd(41,56) + rrt(662) * density(48) 
  pd(48,48) = pd(48,48) - rrt(662) * density(56) 
  pd(48,56) = pd(48,56) - rrt(662) * density(48) 
  pd(56,48) = pd(56,48) - rrt(662) * density(56) 
  pd(56,56) = pd(56,56) - rrt(662) * density(48) 
  pd(23,48) = pd(23,48) + rrt(663) * density(57) 
  pd(23,57) = pd(23,57) + rrt(663) * density(48) 
  pd(30,48) = pd(30,48) + rrt(663) * density(57) 
  pd(30,57) = pd(30,57) + rrt(663) * density(48) 
  pd(41,48) = pd(41,48) + rrt(663) * density(57) 
  pd(41,57) = pd(41,57) + rrt(663) * density(48) 
  pd(48,48) = pd(48,48) - rrt(663) * density(57) 
  pd(48,57) = pd(48,57) - rrt(663) * density(48) 
  pd(57,48) = pd(57,48) - rrt(663) * density(57) 
  pd(57,57) = pd(57,57) - rrt(663) * density(48) 
  pd(01,48) = pd(01,48) + rrt(664) * density(78) 
  pd(01,78) = pd(01,78) + rrt(664) * density(48) 
  pd(30,48) = pd(30,48) + rrt(664) * density(78) 
  pd(30,78) = pd(30,78) + rrt(664) * density(48) 
  pd(41,48) = pd(41,48) + rrt(664) * density(78) 
  pd(41,78) = pd(41,78) + rrt(664) * density(48) 
  pd(48,48) = pd(48,48) - rrt(664) * density(78) 
  pd(48,78) = pd(48,78) - rrt(664) * density(48) 
  pd(78,48) = pd(78,48) - rrt(664) * density(78) 
  pd(78,78) = pd(78,78) - rrt(664) * density(48) 
  pd(23,27) = pd(23,27) + rrt(665) * density(58) * 2.0d0
  pd(23,58) = pd(23,58) + rrt(665) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(665) * density(58) 
  pd(27,58) = pd(27,58) - rrt(665) * density(27) 
  pd(50,27) = pd(50,27) + rrt(665) * density(58) 
  pd(50,58) = pd(50,58) + rrt(665) * density(27) 
  pd(58,27) = pd(58,27) - rrt(665) * density(58) 
  pd(58,58) = pd(58,58) - rrt(665) * density(27) 
  pd(01,28) = pd(01,28) + rrt(666) * density(58) 
  pd(01,58) = pd(01,58) + rrt(666) * density(28) 
  pd(23,28) = pd(23,28) + rrt(666) * density(58) 
  pd(23,58) = pd(23,58) + rrt(666) * density(28) 
  pd(28,28) = pd(28,28) - rrt(666) * density(58) 
  pd(28,58) = pd(28,58) - rrt(666) * density(28) 
  pd(50,28) = pd(50,28) + rrt(666) * density(58) 
  pd(50,58) = pd(50,58) + rrt(666) * density(28) 
  pd(58,28) = pd(58,28) - rrt(666) * density(58) 
  pd(58,58) = pd(58,58) - rrt(666) * density(28) 
  pd(01,29) = pd(01,29) + rrt(667) * density(58) * 2.0d0
  pd(01,58) = pd(01,58) + rrt(667) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(667) * density(58) 
  pd(29,58) = pd(29,58) - rrt(667) * density(29) 
  pd(50,29) = pd(50,29) + rrt(667) * density(58) 
  pd(50,58) = pd(50,58) + rrt(667) * density(29) 
  pd(58,29) = pd(58,29) - rrt(667) * density(58) 
  pd(58,58) = pd(58,58) - rrt(667) * density(29) 
  pd(38,43) = pd(38,43) + rrt(668) * density(58) * 2.0d0
  pd(38,58) = pd(38,58) + rrt(668) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(668) * density(58) 
  pd(43,58) = pd(43,58) - rrt(668) * density(43) 
  pd(50,43) = pd(50,43) + rrt(668) * density(58) 
  pd(50,58) = pd(50,58) + rrt(668) * density(43) 
  pd(58,43) = pd(58,43) - rrt(668) * density(58) 
  pd(58,58) = pd(58,58) - rrt(668) * density(43) 
  pd(30,45) = pd(30,45) + rrt(669) * density(58) * 2.0d0
  pd(30,58) = pd(30,58) + rrt(669) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(669) * density(58) 
  pd(45,58) = pd(45,58) - rrt(669) * density(45) 
  pd(50,45) = pd(50,45) + rrt(669) * density(58) 
  pd(50,58) = pd(50,58) + rrt(669) * density(45) 
  pd(58,45) = pd(58,45) - rrt(669) * density(58) 
  pd(58,58) = pd(58,58) - rrt(669) * density(45) 
  pd(23,55) = pd(23,55) + rrt(670) * density(58) 
  pd(23,58) = pd(23,58) + rrt(670) * density(55) 
  pd(38,55) = pd(38,55) + rrt(670) * density(58) 
  pd(38,58) = pd(38,58) + rrt(670) * density(55) 
  pd(50,55) = pd(50,55) + rrt(670) * density(58) 
  pd(50,58) = pd(50,58) + rrt(670) * density(55) 
  pd(55,55) = pd(55,55) - rrt(670) * density(58) 
  pd(55,58) = pd(55,58) - rrt(670) * density(55) 
  pd(58,55) = pd(58,55) - rrt(670) * density(58) 
  pd(58,58) = pd(58,58) - rrt(670) * density(55) 
  pd(01,56) = pd(01,56) + rrt(671) * density(58) 
  pd(01,58) = pd(01,58) + rrt(671) * density(56) 
  pd(38,56) = pd(38,56) + rrt(671) * density(58) 
  pd(38,58) = pd(38,58) + rrt(671) * density(56) 
  pd(50,56) = pd(50,56) + rrt(671) * density(58) 
  pd(50,58) = pd(50,58) + rrt(671) * density(56) 
  pd(56,56) = pd(56,56) - rrt(671) * density(58) 
  pd(56,58) = pd(56,58) - rrt(671) * density(56) 
  pd(58,56) = pd(58,56) - rrt(671) * density(58) 
  pd(58,58) = pd(58,58) - rrt(671) * density(56) 
  pd(23,57) = pd(23,57) + rrt(672) * density(58) 
  pd(23,58) = pd(23,58) + rrt(672) * density(57) 
  pd(30,57) = pd(30,57) + rrt(672) * density(58) 
  pd(30,58) = pd(30,58) + rrt(672) * density(57) 
  pd(50,57) = pd(50,57) + rrt(672) * density(58) 
  pd(50,58) = pd(50,58) + rrt(672) * density(57) 
  pd(57,57) = pd(57,57) - rrt(672) * density(58) 
  pd(57,58) = pd(57,58) - rrt(672) * density(57) 
  pd(58,57) = pd(58,57) - rrt(672) * density(58) 
  pd(58,58) = pd(58,58) - rrt(672) * density(57) 
  pd(01,58) = pd(01,58) + rrt(673) * density(78) 
  pd(01,78) = pd(01,78) + rrt(673) * density(58) 
  pd(30,58) = pd(30,58) + rrt(673) * density(78) 
  pd(30,78) = pd(30,78) + rrt(673) * density(58) 
  pd(50,58) = pd(50,58) + rrt(673) * density(78) 
  pd(50,78) = pd(50,78) + rrt(673) * density(58) 
  pd(58,58) = pd(58,58) - rrt(673) * density(78) 
  pd(58,78) = pd(58,78) - rrt(673) * density(58) 
  pd(78,58) = pd(78,58) - rrt(673) * density(78) 
  pd(78,78) = pd(78,78) - rrt(673) * density(58) 
  pd(23,27) = pd(23,27) + rrt(674) * density(59) * 2.0d0
  pd(23,59) = pd(23,59) + rrt(674) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(674) * density(59) 
  pd(27,59) = pd(27,59) - rrt(674) * density(27) 
  pd(51,27) = pd(51,27) + rrt(674) * density(59) 
  pd(51,59) = pd(51,59) + rrt(674) * density(27) 
  pd(59,27) = pd(59,27) - rrt(674) * density(59) 
  pd(59,59) = pd(59,59) - rrt(674) * density(27) 
  pd(01,28) = pd(01,28) + rrt(675) * density(59) 
  pd(01,59) = pd(01,59) + rrt(675) * density(28) 
  pd(23,28) = pd(23,28) + rrt(675) * density(59) 
  pd(23,59) = pd(23,59) + rrt(675) * density(28) 
  pd(28,28) = pd(28,28) - rrt(675) * density(59) 
  pd(28,59) = pd(28,59) - rrt(675) * density(28) 
  pd(51,28) = pd(51,28) + rrt(675) * density(59) 
  pd(51,59) = pd(51,59) + rrt(675) * density(28) 
  pd(59,28) = pd(59,28) - rrt(675) * density(59) 
  pd(59,59) = pd(59,59) - rrt(675) * density(28) 
  pd(01,29) = pd(01,29) + rrt(676) * density(59) * 2.0d0
  pd(01,59) = pd(01,59) + rrt(676) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(676) * density(59) 
  pd(29,59) = pd(29,59) - rrt(676) * density(29) 
  pd(51,29) = pd(51,29) + rrt(676) * density(59) 
  pd(51,59) = pd(51,59) + rrt(676) * density(29) 
  pd(59,29) = pd(59,29) - rrt(676) * density(59) 
  pd(59,59) = pd(59,59) - rrt(676) * density(29) 
  pd(38,43) = pd(38,43) + rrt(677) * density(59) * 2.0d0
  pd(38,59) = pd(38,59) + rrt(677) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(677) * density(59) 
  pd(43,59) = pd(43,59) - rrt(677) * density(43) 
  pd(51,43) = pd(51,43) + rrt(677) * density(59) 
  pd(51,59) = pd(51,59) + rrt(677) * density(43) 
  pd(59,43) = pd(59,43) - rrt(677) * density(59) 
  pd(59,59) = pd(59,59) - rrt(677) * density(43) 
  pd(30,45) = pd(30,45) + rrt(678) * density(59) * 2.0d0
  pd(30,59) = pd(30,59) + rrt(678) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(678) * density(59) 
  pd(45,59) = pd(45,59) - rrt(678) * density(45) 
  pd(51,45) = pd(51,45) + rrt(678) * density(59) 
  pd(51,59) = pd(51,59) + rrt(678) * density(45) 
  pd(59,45) = pd(59,45) - rrt(678) * density(59) 
  pd(59,59) = pd(59,59) - rrt(678) * density(45) 
  pd(23,55) = pd(23,55) + rrt(679) * density(59) 
  pd(23,59) = pd(23,59) + rrt(679) * density(55) 
  pd(38,55) = pd(38,55) + rrt(679) * density(59) 
  pd(38,59) = pd(38,59) + rrt(679) * density(55) 
  pd(51,55) = pd(51,55) + rrt(679) * density(59) 
  pd(51,59) = pd(51,59) + rrt(679) * density(55) 
  pd(55,55) = pd(55,55) - rrt(679) * density(59) 
  pd(55,59) = pd(55,59) - rrt(679) * density(55) 
  pd(59,55) = pd(59,55) - rrt(679) * density(59) 
  pd(59,59) = pd(59,59) - rrt(679) * density(55) 
  pd(01,56) = pd(01,56) + rrt(680) * density(59) 
  pd(01,59) = pd(01,59) + rrt(680) * density(56) 
  pd(38,56) = pd(38,56) + rrt(680) * density(59) 
  pd(38,59) = pd(38,59) + rrt(680) * density(56) 
  pd(51,56) = pd(51,56) + rrt(680) * density(59) 
  pd(51,59) = pd(51,59) + rrt(680) * density(56) 
  pd(56,56) = pd(56,56) - rrt(680) * density(59) 
  pd(56,59) = pd(56,59) - rrt(680) * density(56) 
  pd(59,56) = pd(59,56) - rrt(680) * density(59) 
  pd(59,59) = pd(59,59) - rrt(680) * density(56) 
  pd(23,57) = pd(23,57) + rrt(681) * density(59) 
  pd(23,59) = pd(23,59) + rrt(681) * density(57) 
  pd(30,57) = pd(30,57) + rrt(681) * density(59) 
  pd(30,59) = pd(30,59) + rrt(681) * density(57) 
  pd(51,57) = pd(51,57) + rrt(681) * density(59) 
  pd(51,59) = pd(51,59) + rrt(681) * density(57) 
  pd(57,57) = pd(57,57) - rrt(681) * density(59) 
  pd(57,59) = pd(57,59) - rrt(681) * density(57) 
  pd(59,57) = pd(59,57) - rrt(681) * density(59) 
  pd(59,59) = pd(59,59) - rrt(681) * density(57) 
  pd(01,59) = pd(01,59) + rrt(682) * density(78) 
  pd(01,78) = pd(01,78) + rrt(682) * density(59) 
  pd(30,59) = pd(30,59) + rrt(682) * density(78) 
  pd(30,78) = pd(30,78) + rrt(682) * density(59) 
  pd(51,59) = pd(51,59) + rrt(682) * density(78) 
  pd(51,78) = pd(51,78) + rrt(682) * density(59) 
  pd(59,59) = pd(59,59) - rrt(682) * density(78) 
  pd(59,78) = pd(59,78) - rrt(682) * density(59) 
  pd(78,59) = pd(78,59) - rrt(682) * density(78) 
  pd(78,78) = pd(78,78) - rrt(682) * density(59) 
  pd(23,27) = pd(23,27) + rrt(683) * density(60) * 2.0d0
  pd(23,60) = pd(23,60) + rrt(683) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(683) * density(60) 
  pd(27,60) = pd(27,60) - rrt(683) * density(27) 
  pd(52,27) = pd(52,27) + rrt(683) * density(60) 
  pd(52,60) = pd(52,60) + rrt(683) * density(27) 
  pd(60,27) = pd(60,27) - rrt(683) * density(60) 
  pd(60,60) = pd(60,60) - rrt(683) * density(27) 
  pd(01,28) = pd(01,28) + rrt(684) * density(60) 
  pd(01,60) = pd(01,60) + rrt(684) * density(28) 
  pd(23,28) = pd(23,28) + rrt(684) * density(60) 
  pd(23,60) = pd(23,60) + rrt(684) * density(28) 
  pd(28,28) = pd(28,28) - rrt(684) * density(60) 
  pd(28,60) = pd(28,60) - rrt(684) * density(28) 
  pd(52,28) = pd(52,28) + rrt(684) * density(60) 
  pd(52,60) = pd(52,60) + rrt(684) * density(28) 
  pd(60,28) = pd(60,28) - rrt(684) * density(60) 
  pd(60,60) = pd(60,60) - rrt(684) * density(28) 
  pd(01,29) = pd(01,29) + rrt(685) * density(60) * 2.0d0
  pd(01,60) = pd(01,60) + rrt(685) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(685) * density(60) 
  pd(29,60) = pd(29,60) - rrt(685) * density(29) 
  pd(52,29) = pd(52,29) + rrt(685) * density(60) 
  pd(52,60) = pd(52,60) + rrt(685) * density(29) 
  pd(60,29) = pd(60,29) - rrt(685) * density(60) 
  pd(60,60) = pd(60,60) - rrt(685) * density(29) 
  pd(38,43) = pd(38,43) + rrt(686) * density(60) * 2.0d0
  pd(38,60) = pd(38,60) + rrt(686) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(686) * density(60) 
  pd(43,60) = pd(43,60) - rrt(686) * density(43) 
  pd(52,43) = pd(52,43) + rrt(686) * density(60) 
  pd(52,60) = pd(52,60) + rrt(686) * density(43) 
  pd(60,43) = pd(60,43) - rrt(686) * density(60) 
  pd(60,60) = pd(60,60) - rrt(686) * density(43) 
  pd(30,45) = pd(30,45) + rrt(687) * density(60) * 2.0d0
  pd(30,60) = pd(30,60) + rrt(687) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(687) * density(60) 
  pd(45,60) = pd(45,60) - rrt(687) * density(45) 
  pd(52,45) = pd(52,45) + rrt(687) * density(60) 
  pd(52,60) = pd(52,60) + rrt(687) * density(45) 
  pd(60,45) = pd(60,45) - rrt(687) * density(60) 
  pd(60,60) = pd(60,60) - rrt(687) * density(45) 
  pd(23,55) = pd(23,55) + rrt(688) * density(60) 
  pd(23,60) = pd(23,60) + rrt(688) * density(55) 
  pd(38,55) = pd(38,55) + rrt(688) * density(60) 
  pd(38,60) = pd(38,60) + rrt(688) * density(55) 
  pd(52,55) = pd(52,55) + rrt(688) * density(60) 
  pd(52,60) = pd(52,60) + rrt(688) * density(55) 
  pd(55,55) = pd(55,55) - rrt(688) * density(60) 
  pd(55,60) = pd(55,60) - rrt(688) * density(55) 
  pd(60,55) = pd(60,55) - rrt(688) * density(60) 
  pd(60,60) = pd(60,60) - rrt(688) * density(55) 
  pd(01,56) = pd(01,56) + rrt(689) * density(60) 
  pd(01,60) = pd(01,60) + rrt(689) * density(56) 
  pd(38,56) = pd(38,56) + rrt(689) * density(60) 
  pd(38,60) = pd(38,60) + rrt(689) * density(56) 
  pd(52,56) = pd(52,56) + rrt(689) * density(60) 
  pd(52,60) = pd(52,60) + rrt(689) * density(56) 
  pd(56,56) = pd(56,56) - rrt(689) * density(60) 
  pd(56,60) = pd(56,60) - rrt(689) * density(56) 
  pd(60,56) = pd(60,56) - rrt(689) * density(60) 
  pd(60,60) = pd(60,60) - rrt(689) * density(56) 
  pd(23,57) = pd(23,57) + rrt(690) * density(60) 
  pd(23,60) = pd(23,60) + rrt(690) * density(57) 
  pd(30,57) = pd(30,57) + rrt(690) * density(60) 
  pd(30,60) = pd(30,60) + rrt(690) * density(57) 
  pd(52,57) = pd(52,57) + rrt(690) * density(60) 
  pd(52,60) = pd(52,60) + rrt(690) * density(57) 
  pd(57,57) = pd(57,57) - rrt(690) * density(60) 
  pd(57,60) = pd(57,60) - rrt(690) * density(57) 
  pd(60,57) = pd(60,57) - rrt(690) * density(60) 
  pd(60,60) = pd(60,60) - rrt(690) * density(57) 
  pd(01,60) = pd(01,60) + rrt(691) * density(78) 
  pd(01,78) = pd(01,78) + rrt(691) * density(60) 
  pd(30,60) = pd(30,60) + rrt(691) * density(78) 
  pd(30,78) = pd(30,78) + rrt(691) * density(60) 
  pd(52,60) = pd(52,60) + rrt(691) * density(78) 
  pd(52,78) = pd(52,78) + rrt(691) * density(60) 
  pd(60,60) = pd(60,60) - rrt(691) * density(78) 
  pd(60,78) = pd(60,78) - rrt(691) * density(60) 
  pd(78,60) = pd(78,60) - rrt(691) * density(78) 
  pd(78,78) = pd(78,78) - rrt(691) * density(60) 
  pd(23,27) = pd(23,27) + rrt(692) * density(61) * 2.0d0
  pd(23,61) = pd(23,61) + rrt(692) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(692) * density(61) 
  pd(27,61) = pd(27,61) - rrt(692) * density(27) 
  pd(53,27) = pd(53,27) + rrt(692) * density(61) 
  pd(53,61) = pd(53,61) + rrt(692) * density(27) 
  pd(61,27) = pd(61,27) - rrt(692) * density(61) 
  pd(61,61) = pd(61,61) - rrt(692) * density(27) 
  pd(01,28) = pd(01,28) + rrt(693) * density(61) 
  pd(01,61) = pd(01,61) + rrt(693) * density(28) 
  pd(23,28) = pd(23,28) + rrt(693) * density(61) 
  pd(23,61) = pd(23,61) + rrt(693) * density(28) 
  pd(28,28) = pd(28,28) - rrt(693) * density(61) 
  pd(28,61) = pd(28,61) - rrt(693) * density(28) 
  pd(53,28) = pd(53,28) + rrt(693) * density(61) 
  pd(53,61) = pd(53,61) + rrt(693) * density(28) 
  pd(61,28) = pd(61,28) - rrt(693) * density(61) 
  pd(61,61) = pd(61,61) - rrt(693) * density(28) 
  pd(01,29) = pd(01,29) + rrt(694) * density(61) * 2.0d0
  pd(01,61) = pd(01,61) + rrt(694) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(694) * density(61) 
  pd(29,61) = pd(29,61) - rrt(694) * density(29) 
  pd(53,29) = pd(53,29) + rrt(694) * density(61) 
  pd(53,61) = pd(53,61) + rrt(694) * density(29) 
  pd(61,29) = pd(61,29) - rrt(694) * density(61) 
  pd(61,61) = pd(61,61) - rrt(694) * density(29) 
  pd(38,43) = pd(38,43) + rrt(695) * density(61) * 2.0d0
  pd(38,61) = pd(38,61) + rrt(695) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(695) * density(61) 
  pd(43,61) = pd(43,61) - rrt(695) * density(43) 
  pd(53,43) = pd(53,43) + rrt(695) * density(61) 
  pd(53,61) = pd(53,61) + rrt(695) * density(43) 
  pd(61,43) = pd(61,43) - rrt(695) * density(61) 
  pd(61,61) = pd(61,61) - rrt(695) * density(43) 
  pd(30,45) = pd(30,45) + rrt(696) * density(61) * 2.0d0
  pd(30,61) = pd(30,61) + rrt(696) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(696) * density(61) 
  pd(45,61) = pd(45,61) - rrt(696) * density(45) 
  pd(53,45) = pd(53,45) + rrt(696) * density(61) 
  pd(53,61) = pd(53,61) + rrt(696) * density(45) 
  pd(61,45) = pd(61,45) - rrt(696) * density(61) 
  pd(61,61) = pd(61,61) - rrt(696) * density(45) 
  pd(23,55) = pd(23,55) + rrt(697) * density(61) 
  pd(23,61) = pd(23,61) + rrt(697) * density(55) 
  pd(38,55) = pd(38,55) + rrt(697) * density(61) 
  pd(38,61) = pd(38,61) + rrt(697) * density(55) 
  pd(53,55) = pd(53,55) + rrt(697) * density(61) 
  pd(53,61) = pd(53,61) + rrt(697) * density(55) 
  pd(55,55) = pd(55,55) - rrt(697) * density(61) 
  pd(55,61) = pd(55,61) - rrt(697) * density(55) 
  pd(61,55) = pd(61,55) - rrt(697) * density(61) 
  pd(61,61) = pd(61,61) - rrt(697) * density(55) 
  pd(01,56) = pd(01,56) + rrt(698) * density(61) 
  pd(01,61) = pd(01,61) + rrt(698) * density(56) 
  pd(38,56) = pd(38,56) + rrt(698) * density(61) 
  pd(38,61) = pd(38,61) + rrt(698) * density(56) 
  pd(53,56) = pd(53,56) + rrt(698) * density(61) 
  pd(53,61) = pd(53,61) + rrt(698) * density(56) 
  pd(56,56) = pd(56,56) - rrt(698) * density(61) 
  pd(56,61) = pd(56,61) - rrt(698) * density(56) 
  pd(61,56) = pd(61,56) - rrt(698) * density(61) 
  pd(61,61) = pd(61,61) - rrt(698) * density(56) 
  pd(23,57) = pd(23,57) + rrt(699) * density(61) 
  pd(23,61) = pd(23,61) + rrt(699) * density(57) 
  pd(30,57) = pd(30,57) + rrt(699) * density(61) 
  pd(30,61) = pd(30,61) + rrt(699) * density(57) 
  pd(53,57) = pd(53,57) + rrt(699) * density(61) 
  pd(53,61) = pd(53,61) + rrt(699) * density(57) 
  pd(57,57) = pd(57,57) - rrt(699) * density(61) 
  pd(57,61) = pd(57,61) - rrt(699) * density(57) 
  pd(61,57) = pd(61,57) - rrt(699) * density(61) 
  pd(61,61) = pd(61,61) - rrt(699) * density(57) 
  pd(01,61) = pd(01,61) + rrt(700) * density(78) 
  pd(01,78) = pd(01,78) + rrt(700) * density(61) 
  pd(30,61) = pd(30,61) + rrt(700) * density(78) 
  pd(30,78) = pd(30,78) + rrt(700) * density(61) 
  pd(53,61) = pd(53,61) + rrt(700) * density(78) 
  pd(53,78) = pd(53,78) + rrt(700) * density(61) 
  pd(61,61) = pd(61,61) - rrt(700) * density(78) 
  pd(61,78) = pd(61,78) - rrt(700) * density(61) 
  pd(78,61) = pd(78,61) - rrt(700) * density(78) 
  pd(78,78) = pd(78,78) - rrt(700) * density(61) 
  pd(23,26) = pd(23,26) + rrt(701) * density(49) 
  pd(23,49) = pd(23,49) + rrt(701) * density(26) 
  pd(26,26) = pd(26,26) - rrt(701) * density(49) 
  pd(26,49) = pd(26,49) - rrt(701) * density(26) 
  pd(30,26) = pd(30,26) + rrt(701) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(701) * density(26) * 2.0d0
  pd(49,26) = pd(49,26) - rrt(701) * density(49) 
  pd(49,49) = pd(49,49) - rrt(701) * density(26) 
  pd(01,27) = pd(01,27) + rrt(702) * density(49) 
  pd(01,49) = pd(01,49) + rrt(702) * density(27) 
  pd(27,27) = pd(27,27) - rrt(702) * density(49) 
  pd(27,49) = pd(27,49) - rrt(702) * density(27) 
  pd(30,27) = pd(30,27) + rrt(702) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(702) * density(27) * 2.0d0
  pd(49,27) = pd(49,27) - rrt(702) * density(49) 
  pd(49,49) = pd(49,49) - rrt(702) * density(27) 
  pd(30,42) = pd(30,42) + rrt(703) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(703) * density(42) * 2.0d0
  pd(38,42) = pd(38,42) + rrt(703) * density(49) 
  pd(38,49) = pd(38,49) + rrt(703) * density(42) 
  pd(42,42) = pd(42,42) - rrt(703) * density(49) 
  pd(42,49) = pd(42,49) - rrt(703) * density(42) 
  pd(49,42) = pd(49,42) - rrt(703) * density(49) 
  pd(49,49) = pd(49,49) - rrt(703) * density(42) 
  pd(30,43) = pd(30,43) + rrt(704) * density(49) * 3.0d0
  pd(30,49) = pd(30,49) + rrt(704) * density(43) * 3.0d0
  pd(43,43) = pd(43,43) - rrt(704) * density(49) 
  pd(43,49) = pd(43,49) - rrt(704) * density(43) 
  pd(49,43) = pd(49,43) - rrt(704) * density(49) 
  pd(49,49) = pd(49,49) - rrt(704) * density(43) 
  pd(30,49) = pd(30,49) + rrt(705) * density(55) * 2.0d0
  pd(30,55) = pd(30,55) + rrt(705) * density(49) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(705) * density(55) 
  pd(49,55) = pd(49,55) - rrt(705) * density(49) 
  pd(50,49) = pd(50,49) + rrt(705) * density(55) 
  pd(50,55) = pd(50,55) + rrt(705) * density(49) 
  pd(55,49) = pd(55,49) - rrt(705) * density(55) 
  pd(55,55) = pd(55,55) - rrt(705) * density(49) 
  pd(30,49) = pd(30,49) + rrt(706) * density(56) * 2.0d0
  pd(30,56) = pd(30,56) + rrt(706) * density(49) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(706) * density(56) 
  pd(49,56) = pd(49,56) - rrt(706) * density(49) 
  pd(51,49) = pd(51,49) + rrt(706) * density(56) 
  pd(51,56) = pd(51,56) + rrt(706) * density(49) 
  pd(56,49) = pd(56,49) - rrt(706) * density(56) 
  pd(56,56) = pd(56,56) - rrt(706) * density(49) 
  pd(30,49) = pd(30,49) + rrt(707) * density(57) * 2.0d0
  pd(30,57) = pd(30,57) + rrt(707) * density(49) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(707) * density(57) 
  pd(49,57) = pd(49,57) - rrt(707) * density(49) 
  pd(52,49) = pd(52,49) + rrt(707) * density(57) 
  pd(52,57) = pd(52,57) + rrt(707) * density(49) 
  pd(57,49) = pd(57,49) - rrt(707) * density(57) 
  pd(57,57) = pd(57,57) - rrt(707) * density(49) 
  pd(01,28) = pd(01,28) + rrt(708) * density(49) 
  pd(01,49) = pd(01,49) + rrt(708) * density(28) 
  pd(23,28) = pd(23,28) + rrt(708) * density(49) 
  pd(23,49) = pd(23,49) + rrt(708) * density(28) 
  pd(28,28) = pd(28,28) - rrt(708) * density(49) 
  pd(28,49) = pd(28,49) - rrt(708) * density(28) 
  pd(30,28) = pd(30,28) + rrt(708) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(708) * density(28) * 2.0d0
  pd(49,28) = pd(49,28) - rrt(708) * density(49) 
  pd(49,49) = pd(49,49) - rrt(708) * density(28) 
  pd(01,29) = pd(01,29) + rrt(709) * density(49) * 2.0d0
  pd(01,49) = pd(01,49) + rrt(709) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(709) * density(49) 
  pd(29,49) = pd(29,49) - rrt(709) * density(29) 
  pd(30,29) = pd(30,29) + rrt(709) * density(49) * 2.0d0
  pd(30,49) = pd(30,49) + rrt(709) * density(29) * 2.0d0
  pd(49,29) = pd(49,29) - rrt(709) * density(49) 
  pd(49,49) = pd(49,49) - rrt(709) * density(29) 
  pd(30,45) = pd(30,45) + rrt(710) * density(49) * 4.0d0
  pd(30,49) = pd(30,49) + rrt(710) * density(45) * 4.0d0
  pd(45,45) = pd(45,45) - rrt(710) * density(49) 
  pd(45,49) = pd(45,49) - rrt(710) * density(45) 
  pd(49,45) = pd(49,45) - rrt(710) * density(49) 
  pd(49,49) = pd(49,49) - rrt(710) * density(45) 
  pd(01,49) = pd(01,49) + rrt(711) * density(78) 
  pd(01,78) = pd(01,78) + rrt(711) * density(49) 
  pd(30,49) = pd(30,49) + rrt(711) * density(78) * 3.0d0
  pd(30,78) = pd(30,78) + rrt(711) * density(49) * 3.0d0
  pd(49,49) = pd(49,49) - rrt(711) * density(78) 
  pd(49,78) = pd(49,78) - rrt(711) * density(49) 
  pd(78,49) = pd(78,49) - rrt(711) * density(78) 
  pd(78,78) = pd(78,78) - rrt(711) * density(49) 
  pd(23,26) = pd(23,26) + rrt(712) * density(46) 
  pd(23,46) = pd(23,46) + rrt(712) * density(26) 
  pd(26,26) = pd(26,26) - rrt(712) * density(46) 
  pd(26,46) = pd(26,46) - rrt(712) * density(26) 
  pd(38,26) = pd(38,26) + rrt(712) * density(46) 
  pd(38,46) = pd(38,46) + rrt(712) * density(26) 
  pd(46,26) = pd(46,26) - rrt(712) * density(46) 
  pd(46,46) = pd(46,46) - rrt(712) * density(26) 
  pd(01,27) = pd(01,27) + rrt(713) * density(46) 
  pd(01,46) = pd(01,46) + rrt(713) * density(27) 
  pd(27,27) = pd(27,27) - rrt(713) * density(46) 
  pd(27,46) = pd(27,46) - rrt(713) * density(27) 
  pd(38,27) = pd(38,27) + rrt(713) * density(46) 
  pd(38,46) = pd(38,46) + rrt(713) * density(27) 
  pd(46,27) = pd(46,27) - rrt(713) * density(46) 
  pd(46,46) = pd(46,46) - rrt(713) * density(27) 
  pd(38,42) = pd(38,42) + rrt(714) * density(46) * 2.0d0
  pd(38,46) = pd(38,46) + rrt(714) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(714) * density(46) 
  pd(42,46) = pd(42,46) - rrt(714) * density(42) 
  pd(46,42) = pd(46,42) - rrt(714) * density(46) 
  pd(46,46) = pd(46,46) - rrt(714) * density(42) 
  pd(30,43) = pd(30,43) + rrt(715) * density(46) 
  pd(30,46) = pd(30,46) + rrt(715) * density(43) 
  pd(38,43) = pd(38,43) + rrt(715) * density(46) 
  pd(38,46) = pd(38,46) + rrt(715) * density(43) 
  pd(43,43) = pd(43,43) - rrt(715) * density(46) 
  pd(43,46) = pd(43,46) - rrt(715) * density(43) 
  pd(46,43) = pd(46,43) - rrt(715) * density(46) 
  pd(46,46) = pd(46,46) - rrt(715) * density(43) 
  pd(38,46) = pd(38,46) + rrt(716) * density(55) 
  pd(38,55) = pd(38,55) + rrt(716) * density(46) 
  pd(46,46) = pd(46,46) - rrt(716) * density(55) 
  pd(46,55) = pd(46,55) - rrt(716) * density(46) 
  pd(50,46) = pd(50,46) + rrt(716) * density(55) 
  pd(50,55) = pd(50,55) + rrt(716) * density(46) 
  pd(55,46) = pd(55,46) - rrt(716) * density(55) 
  pd(55,55) = pd(55,55) - rrt(716) * density(46) 
  pd(23,26) = pd(23,26) + rrt(717) * density(47) 
  pd(23,47) = pd(23,47) + rrt(717) * density(26) 
  pd(26,26) = pd(26,26) - rrt(717) * density(47) 
  pd(26,47) = pd(26,47) - rrt(717) * density(26) 
  pd(30,26) = pd(30,26) + rrt(717) * density(47) 
  pd(30,47) = pd(30,47) + rrt(717) * density(26) 
  pd(47,26) = pd(47,26) - rrt(717) * density(47) 
  pd(47,47) = pd(47,47) - rrt(717) * density(26) 
  pd(01,27) = pd(01,27) + rrt(718) * density(47) 
  pd(01,47) = pd(01,47) + rrt(718) * density(27) 
  pd(27,27) = pd(27,27) - rrt(718) * density(47) 
  pd(27,47) = pd(27,47) - rrt(718) * density(27) 
  pd(30,27) = pd(30,27) + rrt(718) * density(47) 
  pd(30,47) = pd(30,47) + rrt(718) * density(27) 
  pd(47,27) = pd(47,27) - rrt(718) * density(47) 
  pd(47,47) = pd(47,47) - rrt(718) * density(27) 
  pd(30,42) = pd(30,42) + rrt(719) * density(47) 
  pd(30,47) = pd(30,47) + rrt(719) * density(42) 
  pd(38,42) = pd(38,42) + rrt(719) * density(47) 
  pd(38,47) = pd(38,47) + rrt(719) * density(42) 
  pd(42,42) = pd(42,42) - rrt(719) * density(47) 
  pd(42,47) = pd(42,47) - rrt(719) * density(42) 
  pd(47,42) = pd(47,42) - rrt(719) * density(47) 
  pd(47,47) = pd(47,47) - rrt(719) * density(42) 
  pd(30,43) = pd(30,43) + rrt(720) * density(47) * 2.0d0
  pd(30,47) = pd(30,47) + rrt(720) * density(43) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(720) * density(47) 
  pd(43,47) = pd(43,47) - rrt(720) * density(43) 
  pd(47,43) = pd(47,43) - rrt(720) * density(47) 
  pd(47,47) = pd(47,47) - rrt(720) * density(43) 
  pd(30,47) = pd(30,47) + rrt(721) * density(55) 
  pd(30,55) = pd(30,55) + rrt(721) * density(47) 
  pd(47,47) = pd(47,47) - rrt(721) * density(55) 
  pd(47,55) = pd(47,55) - rrt(721) * density(47) 
  pd(50,47) = pd(50,47) + rrt(721) * density(55) 
  pd(50,55) = pd(50,55) + rrt(721) * density(47) 
  pd(55,47) = pd(55,47) - rrt(721) * density(55) 
  pd(55,55) = pd(55,55) - rrt(721) * density(47) 
  pd(26,26) = pd(26,26) - rrt(722) * density(46) 
  pd(26,46) = pd(26,46) - rrt(722) * density(26) 
  pd(46,26) = pd(46,26) - rrt(722) * density(46) 
  pd(46,46) = pd(46,46) - rrt(722) * density(26) 
  pd(50,26) = pd(50,26) + rrt(722) * density(46) 
  pd(50,46) = pd(50,46) + rrt(722) * density(26) 
  pd(27,27) = pd(27,27) - rrt(723) * density(46) 
  pd(27,46) = pd(27,46) - rrt(723) * density(27) 
  pd(46,27) = pd(46,27) - rrt(723) * density(46) 
  pd(46,46) = pd(46,46) - rrt(723) * density(27) 
  pd(51,27) = pd(51,27) + rrt(723) * density(46) 
  pd(51,46) = pd(51,46) + rrt(723) * density(27) 
  pd(30,42) = pd(30,42) + rrt(724) * density(46) 
  pd(30,46) = pd(30,46) + rrt(724) * density(42) 
  pd(42,42) = pd(42,42) - rrt(724) * density(46) 
  pd(42,46) = pd(42,46) - rrt(724) * density(42) 
  pd(46,42) = pd(46,42) - rrt(724) * density(46) 
  pd(46,46) = pd(46,46) - rrt(724) * density(42) 
  pd(41,43) = pd(41,43) + rrt(725) * density(46) 
  pd(41,46) = pd(41,46) + rrt(725) * density(43) 
  pd(43,43) = pd(43,43) - rrt(725) * density(46) 
  pd(43,46) = pd(43,46) - rrt(725) * density(43) 
  pd(46,43) = pd(46,43) - rrt(725) * density(46) 
  pd(46,46) = pd(46,46) - rrt(725) * density(43) 
  pd(46,46) = pd(46,46) - rrt(726) * density(55) 
  pd(46,55) = pd(46,55) - rrt(726) * density(46) 
  pd(52,46) = pd(52,46) + rrt(726) * density(55) 
  pd(52,55) = pd(52,55) + rrt(726) * density(46) 
  pd(55,46) = pd(55,46) - rrt(726) * density(55) 
  pd(55,55) = pd(55,55) - rrt(726) * density(46) 
  pd(26,26) = pd(26,26) - rrt(727) * density(47) 
  pd(26,47) = pd(26,47) - rrt(727) * density(26) 
  pd(47,26) = pd(47,26) - rrt(727) * density(47) 
  pd(47,47) = pd(47,47) - rrt(727) * density(26) 
  pd(52,26) = pd(52,26) + rrt(727) * density(47) 
  pd(52,47) = pd(52,47) + rrt(727) * density(26) 
  pd(41,42) = pd(41,42) + rrt(728) * density(47) 
  pd(41,47) = pd(41,47) + rrt(728) * density(42) 
  pd(42,42) = pd(42,42) - rrt(728) * density(47) 
  pd(42,47) = pd(42,47) - rrt(728) * density(42) 
  pd(47,42) = pd(47,42) - rrt(728) * density(47) 
  pd(47,47) = pd(47,47) - rrt(728) * density(42) 
  pd(47,47) = pd(47,47) - rrt(729) * density(55) 
  pd(47,55) = pd(47,55) - rrt(729) * density(47) 
  pd(53,47) = pd(53,47) + rrt(729) * density(55) 
  pd(53,55) = pd(53,55) + rrt(729) * density(47) 
  pd(55,47) = pd(55,47) - rrt(729) * density(55) 
  pd(55,55) = pd(55,55) - rrt(729) * density(47) 
  pd(30,47) = pd(30,47) + rrt(730) * density(73) 
  pd(30,73) = pd(30,73) + rrt(730) * density(47) 
  pd(47,47) = pd(47,47) - rrt(730) * density(73) 
  pd(47,73) = pd(47,73) - rrt(730) * density(47) 
  pd(70,47) = pd(70,47) + rrt(730) * density(73) 
  pd(70,73) = pd(70,73) + rrt(730) * density(47) 
  pd(73,47) = pd(73,47) - rrt(730) * density(73) 
  pd(73,73) = pd(73,73) - rrt(730) * density(47) 
  pd(38,46) = pd(38,46) + rrt(731) * density(73) 
  pd(38,73) = pd(38,73) + rrt(731) * density(46) 
  pd(46,46) = pd(46,46) - rrt(731) * density(73) 
  pd(46,73) = pd(46,73) - rrt(731) * density(46) 
  pd(70,46) = pd(70,46) + rrt(731) * density(73) 
  pd(70,73) = pd(70,73) + rrt(731) * density(46) 
  pd(73,46) = pd(73,46) - rrt(731) * density(73) 
  pd(73,73) = pd(73,73) - rrt(731) * density(46) 
  pd(41,48) = pd(41,48) + rrt(732) * density(73) 
  pd(41,73) = pd(41,73) + rrt(732) * density(48) 
  pd(48,48) = pd(48,48) - rrt(732) * density(73) 
  pd(48,73) = pd(48,73) - rrt(732) * density(48) 
  pd(70,48) = pd(70,48) + rrt(732) * density(73) 
  pd(70,73) = pd(70,73) + rrt(732) * density(48) 
  pd(73,48) = pd(73,48) - rrt(732) * density(73) 
  pd(73,73) = pd(73,73) - rrt(732) * density(48) 
  pd(50,58) = pd(50,58) + rrt(733) * density(73) 
  pd(50,73) = pd(50,73) + rrt(733) * density(58) 
  pd(58,58) = pd(58,58) - rrt(733) * density(73) 
  pd(58,73) = pd(58,73) - rrt(733) * density(58) 
  pd(70,58) = pd(70,58) + rrt(733) * density(73) 
  pd(70,73) = pd(70,73) + rrt(733) * density(58) 
  pd(73,58) = pd(73,58) - rrt(733) * density(73) 
  pd(73,73) = pd(73,73) - rrt(733) * density(58) 
  pd(52,60) = pd(52,60) + rrt(734) * density(73) 
  pd(52,73) = pd(52,73) + rrt(734) * density(60) 
  pd(60,60) = pd(60,60) - rrt(734) * density(73) 
  pd(60,73) = pd(60,73) - rrt(734) * density(60) 
  pd(70,60) = pd(70,60) + rrt(734) * density(73) 
  pd(70,73) = pd(70,73) + rrt(734) * density(60) 
  pd(73,60) = pd(73,60) - rrt(734) * density(73) 
  pd(73,73) = pd(73,73) - rrt(734) * density(60) 
  pd(53,61) = pd(53,61) + rrt(735) * density(73) 
  pd(53,73) = pd(53,73) + rrt(735) * density(61) 
  pd(61,61) = pd(61,61) - rrt(735) * density(73) 
  pd(61,73) = pd(61,73) - rrt(735) * density(61) 
  pd(70,61) = pd(70,61) + rrt(735) * density(73) 
  pd(70,73) = pd(70,73) + rrt(735) * density(61) 
  pd(73,61) = pd(73,61) - rrt(735) * density(73) 
  pd(73,73) = pd(73,73) - rrt(735) * density(61) 
  pd(51,59) = pd(51,59) + rrt(736) * density(73) 
  pd(51,73) = pd(51,73) + rrt(736) * density(59) 
  pd(59,59) = pd(59,59) - rrt(736) * density(73) 
  pd(59,73) = pd(59,73) - rrt(736) * density(59) 
  pd(70,59) = pd(70,59) + rrt(736) * density(73) 
  pd(70,73) = pd(70,73) + rrt(736) * density(59) 
  pd(73,59) = pd(73,59) - rrt(736) * density(73) 
  pd(73,73) = pd(73,73) - rrt(736) * density(59) 
  pd(53,61) = pd(53,61) + rrt(737) * density(77) 
  pd(53,77) = pd(53,77) + rrt(737) * density(61) 
  pd(61,61) = pd(61,61) - rrt(737) * density(77) 
  pd(61,77) = pd(61,77) - rrt(737) * density(61) 
  pd(63,61) = pd(63,61) + rrt(737) * density(77) 
  pd(63,77) = pd(63,77) + rrt(737) * density(61) 
  pd(70,61) = pd(70,61) + rrt(737) * density(77) 
  pd(70,77) = pd(70,77) + rrt(737) * density(61) 
  pd(77,61) = pd(77,61) - rrt(737) * density(77) 
  pd(77,77) = pd(77,77) - rrt(737) * density(61) 
  pd(63,65) = pd(63,65) + rrt(738) * density(67) * 3.0d0
  pd(63,67) = pd(63,67) + rrt(738) * density(65) * 3.0d0
  pd(65,65) = pd(65,65) - rrt(738) * density(67) 
  pd(65,67) = pd(65,67) - rrt(738) * density(65) 
  pd(67,65) = pd(67,65) - rrt(738) * density(67) 
  pd(67,67) = pd(67,67) - rrt(738) * density(65) 
  pd(62,66) = pd(62,66) + rrt(739) * density(67) 
  pd(62,67) = pd(62,67) + rrt(739) * density(66) 
  pd(63,66) = pd(63,66) + rrt(739) * density(67) * 2.0d0
  pd(63,67) = pd(63,67) + rrt(739) * density(66) * 2.0d0
  pd(66,66) = pd(66,66) - rrt(739) * density(67) 
  pd(66,67) = pd(66,67) - rrt(739) * density(66) 
  pd(67,66) = pd(67,66) - rrt(739) * density(67) 
  pd(67,67) = pd(67,67) - rrt(739) * density(66) 
  pd(01,27) = pd(01,27) + rrt(740) * density(67) 
  pd(01,67) = pd(01,67) + rrt(740) * density(27) 
  pd(27,27) = pd(27,27) - rrt(740) * density(67) 
  pd(27,67) = pd(27,67) - rrt(740) * density(27) 
  pd(63,27) = pd(63,27) + rrt(740) * density(67) 
  pd(63,67) = pd(63,67) + rrt(740) * density(27) 
  pd(67,27) = pd(67,27) - rrt(740) * density(67) 
  pd(67,67) = pd(67,67) - rrt(740) * density(27) 
  pd(01,29) = pd(01,29) + rrt(741) * density(67) * 2.0d0
  pd(01,67) = pd(01,67) + rrt(741) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(741) * density(67) 
  pd(29,67) = pd(29,67) - rrt(741) * density(29) 
  pd(63,29) = pd(63,29) + rrt(741) * density(67) 
  pd(63,67) = pd(63,67) + rrt(741) * density(29) 
  pd(67,29) = pd(67,29) - rrt(741) * density(67) 
  pd(67,67) = pd(67,67) - rrt(741) * density(29) 
  pd(01,67) = pd(01,67) + rrt(742) * density(76) 
  pd(01,76) = pd(01,76) + rrt(742) * density(67) 
  pd(62,67) = pd(62,67) + rrt(742) * density(76) 
  pd(62,76) = pd(62,76) + rrt(742) * density(67) 
  pd(67,67) = pd(67,67) - rrt(742) * density(76) 
  pd(67,76) = pd(67,76) - rrt(742) * density(67) 
  pd(76,67) = pd(76,67) - rrt(742) * density(76) 
  pd(76,76) = pd(76,76) - rrt(742) * density(67) 
  pd(62,01) = pd(62,01) + rrt(743) * density(65) * density(67) 
  pd(62,65) = pd(62,65) + rrt(743) * density(01) * density(67) 
  pd(62,67) = pd(62,67) + rrt(743) * density(01) * density(65) 
  pd(63,01) = pd(63,01) + rrt(743) * density(65) * density(67) 
  pd(63,65) = pd(63,65) + rrt(743) * density(01) * density(67) 
  pd(63,67) = pd(63,67) + rrt(743) * density(01) * density(65) 
  pd(65,01) = pd(65,01) - rrt(743) * density(65) * density(67) 
  pd(65,65) = pd(65,65) - rrt(743) * density(01) * density(67) 
  pd(65,67) = pd(65,67) - rrt(743) * density(01) * density(65) 
  pd(67,01) = pd(67,01) - rrt(743) * density(65) * density(67) 
  pd(67,65) = pd(67,65) - rrt(743) * density(01) * density(67) 
  pd(67,67) = pd(67,67) - rrt(743) * density(01) * density(65) 
  pd(62,62) = pd(62,62) + rrt(744) * density(65) * density(67) 
  pd(62,65) = pd(62,65) + rrt(744) * density(62) * density(67) 
  pd(62,67) = pd(62,67) + rrt(744) * density(62) * density(65) 
  pd(63,62) = pd(63,62) + rrt(744) * density(65) * density(67) 
  pd(63,65) = pd(63,65) + rrt(744) * density(62) * density(67) 
  pd(63,67) = pd(63,67) + rrt(744) * density(62) * density(65) 
  pd(65,62) = pd(65,62) - rrt(744) * density(65) * density(67) 
  pd(65,65) = pd(65,65) - rrt(744) * density(62) * density(67) 
  pd(65,67) = pd(65,67) - rrt(744) * density(62) * density(65) 
  pd(67,62) = pd(67,62) - rrt(744) * density(65) * density(67) 
  pd(67,65) = pd(67,65) - rrt(744) * density(62) * density(67) 
  pd(67,67) = pd(67,67) - rrt(744) * density(62) * density(65) 
  pd(62,65) = pd(62,65) + rrt(745) * density(67) * density(70) 
  pd(62,67) = pd(62,67) + rrt(745) * density(65) * density(70) 
  pd(62,70) = pd(62,70) + rrt(745) * density(65) * density(67) 
  pd(63,65) = pd(63,65) + rrt(745) * density(67) * density(70) 
  pd(63,67) = pd(63,67) + rrt(745) * density(65) * density(70) 
  pd(63,70) = pd(63,70) + rrt(745) * density(65) * density(67) 
  pd(65,65) = pd(65,65) - rrt(745) * density(67) * density(70) 
  pd(65,67) = pd(65,67) - rrt(745) * density(65) * density(70) 
  pd(65,70) = pd(65,70) - rrt(745) * density(65) * density(67) 
  pd(67,65) = pd(67,65) - rrt(745) * density(67) * density(70) 
  pd(67,67) = pd(67,67) - rrt(745) * density(65) * density(70) 
  pd(67,70) = pd(67,70) - rrt(745) * density(65) * density(67) 
  pd(62,01) = pd(62,01) + rrt(746) * density(66) * density(67) * 2.0d0
  pd(62,66) = pd(62,66) + rrt(746) * density(01) * density(67) * 2.0d0
  pd(62,67) = pd(62,67) + rrt(746) * density(01) * density(66) * 2.0d0
  pd(66,01) = pd(66,01) - rrt(746) * density(66) * density(67) 
  pd(66,66) = pd(66,66) - rrt(746) * density(01) * density(67) 
  pd(66,67) = pd(66,67) - rrt(746) * density(01) * density(66) 
  pd(67,01) = pd(67,01) - rrt(746) * density(66) * density(67) 
  pd(67,66) = pd(67,66) - rrt(746) * density(01) * density(67) 
  pd(67,67) = pd(67,67) - rrt(746) * density(01) * density(66) 
  pd(62,62) = pd(62,62) + rrt(747) * density(66) * density(67) * 2.0d0
  pd(62,66) = pd(62,66) + rrt(747) * density(62) * density(67) * 2.0d0
  pd(62,67) = pd(62,67) + rrt(747) * density(62) * density(66) * 2.0d0
  pd(66,62) = pd(66,62) - rrt(747) * density(66) * density(67) 
  pd(66,66) = pd(66,66) - rrt(747) * density(62) * density(67) 
  pd(66,67) = pd(66,67) - rrt(747) * density(62) * density(66) 
  pd(67,62) = pd(67,62) - rrt(747) * density(66) * density(67) 
  pd(67,66) = pd(67,66) - rrt(747) * density(62) * density(67) 
  pd(67,67) = pd(67,67) - rrt(747) * density(62) * density(66) 
  pd(62,66) = pd(62,66) + rrt(748) * density(67) * density(70) * 2.0d0
  pd(62,67) = pd(62,67) + rrt(748) * density(66) * density(70) * 2.0d0
  pd(62,70) = pd(62,70) + rrt(748) * density(66) * density(67) * 2.0d0
  pd(66,66) = pd(66,66) - rrt(748) * density(67) * density(70) 
  pd(66,67) = pd(66,67) - rrt(748) * density(66) * density(70) 
  pd(66,70) = pd(66,70) - rrt(748) * density(66) * density(67) 
  pd(67,66) = pd(67,66) - rrt(748) * density(67) * density(70) 
  pd(67,67) = pd(67,67) - rrt(748) * density(66) * density(70) 
  pd(67,70) = pd(67,70) - rrt(748) * density(66) * density(67) 
  pd(01,01) = pd(01,01) + rrt(749) * density(27) * density(67) 
  pd(01,27) = pd(01,27) + rrt(749) * density(01) * density(67) 
  pd(01,67) = pd(01,67) + rrt(749) * density(01) * density(27) 
  pd(27,01) = pd(27,01) - rrt(749) * density(27) * density(67) 
  pd(27,27) = pd(27,27) - rrt(749) * density(01) * density(67) 
  pd(27,67) = pd(27,67) - rrt(749) * density(01) * density(27) 
  pd(63,01) = pd(63,01) + rrt(749) * density(27) * density(67) 
  pd(63,27) = pd(63,27) + rrt(749) * density(01) * density(67) 
  pd(63,67) = pd(63,67) + rrt(749) * density(01) * density(27) 
  pd(67,01) = pd(67,01) - rrt(749) * density(27) * density(67) 
  pd(67,27) = pd(67,27) - rrt(749) * density(01) * density(67) 
  pd(67,67) = pd(67,67) - rrt(749) * density(01) * density(27) 
  pd(01,27) = pd(01,27) + rrt(750) * density(62) * density(67) 
  pd(01,62) = pd(01,62) + rrt(750) * density(27) * density(67) 
  pd(01,67) = pd(01,67) + rrt(750) * density(27) * density(62) 
  pd(27,27) = pd(27,27) - rrt(750) * density(62) * density(67) 
  pd(27,62) = pd(27,62) - rrt(750) * density(27) * density(67) 
  pd(27,67) = pd(27,67) - rrt(750) * density(27) * density(62) 
  pd(63,27) = pd(63,27) + rrt(750) * density(62) * density(67) 
  pd(63,62) = pd(63,62) + rrt(750) * density(27) * density(67) 
  pd(63,67) = pd(63,67) + rrt(750) * density(27) * density(62) 
  pd(67,27) = pd(67,27) - rrt(750) * density(62) * density(67) 
  pd(67,62) = pd(67,62) - rrt(750) * density(27) * density(67) 
  pd(67,67) = pd(67,67) - rrt(750) * density(27) * density(62) 
  pd(01,27) = pd(01,27) + rrt(751) * density(67) * density(70) 
  pd(01,67) = pd(01,67) + rrt(751) * density(27) * density(70) 
  pd(01,70) = pd(01,70) + rrt(751) * density(27) * density(67) 
  pd(27,27) = pd(27,27) - rrt(751) * density(67) * density(70) 
  pd(27,67) = pd(27,67) - rrt(751) * density(27) * density(70) 
  pd(27,70) = pd(27,70) - rrt(751) * density(27) * density(67) 
  pd(63,27) = pd(63,27) + rrt(751) * density(67) * density(70) 
  pd(63,67) = pd(63,67) + rrt(751) * density(27) * density(70) 
  pd(63,70) = pd(63,70) + rrt(751) * density(27) * density(67) 
  pd(67,27) = pd(67,27) - rrt(751) * density(67) * density(70) 
  pd(67,67) = pd(67,67) - rrt(751) * density(27) * density(70) 
  pd(67,70) = pd(67,70) - rrt(751) * density(27) * density(67) 
  pd(01,01) = pd(01,01) + rrt(752) * density(29) * density(67) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(752) * density(01) * density(67) * 2.0d0
  pd(01,67) = pd(01,67) + rrt(752) * density(01) * density(29) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(752) * density(29) * density(67) 
  pd(29,29) = pd(29,29) - rrt(752) * density(01) * density(67) 
  pd(29,67) = pd(29,67) - rrt(752) * density(01) * density(29) 
  pd(63,01) = pd(63,01) + rrt(752) * density(29) * density(67) 
  pd(63,29) = pd(63,29) + rrt(752) * density(01) * density(67) 
  pd(63,67) = pd(63,67) + rrt(752) * density(01) * density(29) 
  pd(67,01) = pd(67,01) - rrt(752) * density(29) * density(67) 
  pd(67,29) = pd(67,29) - rrt(752) * density(01) * density(67) 
  pd(67,67) = pd(67,67) - rrt(752) * density(01) * density(29) 
  pd(01,29) = pd(01,29) + rrt(753) * density(62) * density(67) * 2.0d0
  pd(01,62) = pd(01,62) + rrt(753) * density(29) * density(67) * 2.0d0
  pd(01,67) = pd(01,67) + rrt(753) * density(29) * density(62) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(753) * density(62) * density(67) 
  pd(29,62) = pd(29,62) - rrt(753) * density(29) * density(67) 
  pd(29,67) = pd(29,67) - rrt(753) * density(29) * density(62) 
  pd(63,29) = pd(63,29) + rrt(753) * density(62) * density(67) 
  pd(63,62) = pd(63,62) + rrt(753) * density(29) * density(67) 
  pd(63,67) = pd(63,67) + rrt(753) * density(29) * density(62) 
  pd(67,29) = pd(67,29) - rrt(753) * density(62) * density(67) 
  pd(67,62) = pd(67,62) - rrt(753) * density(29) * density(67) 
  pd(67,67) = pd(67,67) - rrt(753) * density(29) * density(62) 
  pd(01,29) = pd(01,29) + rrt(754) * density(67) * density(70) * 2.0d0
  pd(01,67) = pd(01,67) + rrt(754) * density(29) * density(70) * 2.0d0
  pd(01,70) = pd(01,70) + rrt(754) * density(29) * density(67) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(754) * density(67) * density(70) 
  pd(29,67) = pd(29,67) - rrt(754) * density(29) * density(70) 
  pd(29,70) = pd(29,70) - rrt(754) * density(29) * density(67) 
  pd(63,29) = pd(63,29) + rrt(754) * density(67) * density(70) 
  pd(63,67) = pd(63,67) + rrt(754) * density(29) * density(70) 
  pd(63,70) = pd(63,70) + rrt(754) * density(29) * density(67) 
  pd(67,29) = pd(67,29) - rrt(754) * density(67) * density(70) 
  pd(67,67) = pd(67,67) - rrt(754) * density(29) * density(70) 
  pd(67,70) = pd(67,70) - rrt(754) * density(29) * density(67) 
  pd(01,01) = pd(01,01) + rrt(755) * density(67) * density(76) 
  pd(01,67) = pd(01,67) + rrt(755) * density(01) * density(76) 
  pd(01,76) = pd(01,76) + rrt(755) * density(01) * density(67) 
  pd(62,01) = pd(62,01) + rrt(755) * density(67) * density(76) 
  pd(62,67) = pd(62,67) + rrt(755) * density(01) * density(76) 
  pd(62,76) = pd(62,76) + rrt(755) * density(01) * density(67) 
  pd(67,01) = pd(67,01) - rrt(755) * density(67) * density(76) 
  pd(67,67) = pd(67,67) - rrt(755) * density(01) * density(76) 
  pd(67,76) = pd(67,76) - rrt(755) * density(01) * density(67) 
  pd(76,01) = pd(76,01) - rrt(755) * density(67) * density(76) 
  pd(76,67) = pd(76,67) - rrt(755) * density(01) * density(76) 
  pd(76,76) = pd(76,76) - rrt(755) * density(01) * density(67) 
  pd(01,62) = pd(01,62) + rrt(756) * density(67) * density(76) 
  pd(01,67) = pd(01,67) + rrt(756) * density(62) * density(76) 
  pd(01,76) = pd(01,76) + rrt(756) * density(62) * density(67) 
  pd(62,62) = pd(62,62) + rrt(756) * density(67) * density(76) 
  pd(62,67) = pd(62,67) + rrt(756) * density(62) * density(76) 
  pd(62,76) = pd(62,76) + rrt(756) * density(62) * density(67) 
  pd(67,62) = pd(67,62) - rrt(756) * density(67) * density(76) 
  pd(67,67) = pd(67,67) - rrt(756) * density(62) * density(76) 
  pd(67,76) = pd(67,76) - rrt(756) * density(62) * density(67) 
  pd(76,62) = pd(76,62) - rrt(756) * density(67) * density(76) 
  pd(76,67) = pd(76,67) - rrt(756) * density(62) * density(76) 
  pd(76,76) = pd(76,76) - rrt(756) * density(62) * density(67) 
  pd(01,67) = pd(01,67) + rrt(757) * density(70) * density(76) 
  pd(01,70) = pd(01,70) + rrt(757) * density(67) * density(76) 
  pd(01,76) = pd(01,76) + rrt(757) * density(67) * density(70) 
  pd(62,67) = pd(62,67) + rrt(757) * density(70) * density(76) 
  pd(62,70) = pd(62,70) + rrt(757) * density(67) * density(76) 
  pd(62,76) = pd(62,76) + rrt(757) * density(67) * density(70) 
  pd(67,67) = pd(67,67) - rrt(757) * density(70) * density(76) 
  pd(67,70) = pd(67,70) - rrt(757) * density(67) * density(76) 
  pd(67,76) = pd(67,76) - rrt(757) * density(67) * density(70) 
  pd(76,67) = pd(76,67) - rrt(757) * density(70) * density(76) 
  pd(76,70) = pd(76,70) - rrt(757) * density(67) * density(76) 
  pd(76,76) = pd(76,76) - rrt(757) * density(67) * density(70) 
  pd(23,26) = pd(23,26) + rrt(758) * density(48) 
  pd(23,48) = pd(23,48) + rrt(758) * density(26) 
  pd(26,26) = pd(26,26) - rrt(758) * density(48) 
  pd(26,48) = pd(26,48) - rrt(758) * density(26) 
  pd(41,26) = pd(41,26) + rrt(758) * density(48) 
  pd(41,48) = pd(41,48) + rrt(758) * density(26) 
  pd(48,26) = pd(48,26) - rrt(758) * density(48) 
  pd(48,48) = pd(48,48) - rrt(758) * density(26) 
  pd(01,27) = pd(01,27) + rrt(759) * density(48) 
  pd(01,48) = pd(01,48) + rrt(759) * density(27) 
  pd(27,27) = pd(27,27) - rrt(759) * density(48) 
  pd(27,48) = pd(27,48) - rrt(759) * density(27) 
  pd(41,27) = pd(41,27) + rrt(759) * density(48) 
  pd(41,48) = pd(41,48) + rrt(759) * density(27) 
  pd(48,27) = pd(48,27) - rrt(759) * density(48) 
  pd(48,48) = pd(48,48) - rrt(759) * density(27) 
  pd(38,42) = pd(38,42) + rrt(760) * density(48) 
  pd(38,48) = pd(38,48) + rrt(760) * density(42) 
  pd(41,42) = pd(41,42) + rrt(760) * density(48) 
  pd(41,48) = pd(41,48) + rrt(760) * density(42) 
  pd(42,42) = pd(42,42) - rrt(760) * density(48) 
  pd(42,48) = pd(42,48) - rrt(760) * density(42) 
  pd(48,42) = pd(48,42) - rrt(760) * density(48) 
  pd(48,48) = pd(48,48) - rrt(760) * density(42) 
  pd(30,43) = pd(30,43) + rrt(761) * density(48) 
  pd(30,48) = pd(30,48) + rrt(761) * density(43) 
  pd(41,43) = pd(41,43) + rrt(761) * density(48) 
  pd(41,48) = pd(41,48) + rrt(761) * density(43) 
  pd(43,43) = pd(43,43) - rrt(761) * density(48) 
  pd(43,48) = pd(43,48) - rrt(761) * density(43) 
  pd(48,43) = pd(48,43) - rrt(761) * density(48) 
  pd(48,48) = pd(48,48) - rrt(761) * density(43) 
  pd(41,48) = pd(41,48) + rrt(762) * density(55) 
  pd(41,55) = pd(41,55) + rrt(762) * density(48) 
  pd(48,48) = pd(48,48) - rrt(762) * density(55) 
  pd(48,55) = pd(48,55) - rrt(762) * density(48) 
  pd(50,48) = pd(50,48) + rrt(762) * density(55) 
  pd(50,55) = pd(50,55) + rrt(762) * density(48) 
  pd(55,48) = pd(55,48) - rrt(762) * density(55) 
  pd(55,55) = pd(55,55) - rrt(762) * density(48) 
  pd(41,48) = pd(41,48) + rrt(763) * density(56) 
  pd(41,56) = pd(41,56) + rrt(763) * density(48) 
  pd(48,48) = pd(48,48) - rrt(763) * density(56) 
  pd(48,56) = pd(48,56) - rrt(763) * density(48) 
  pd(51,48) = pd(51,48) + rrt(763) * density(56) 
  pd(51,56) = pd(51,56) + rrt(763) * density(48) 
  pd(56,48) = pd(56,48) - rrt(763) * density(56) 
  pd(56,56) = pd(56,56) - rrt(763) * density(48) 
  pd(41,48) = pd(41,48) + rrt(764) * density(57) 
  pd(41,57) = pd(41,57) + rrt(764) * density(48) 
  pd(48,48) = pd(48,48) - rrt(764) * density(57) 
  pd(48,57) = pd(48,57) - rrt(764) * density(48) 
  pd(52,48) = pd(52,48) + rrt(764) * density(57) 
  pd(52,57) = pd(52,57) + rrt(764) * density(48) 
  pd(57,48) = pd(57,48) - rrt(764) * density(57) 
  pd(57,57) = pd(57,57) - rrt(764) * density(48) 
  pd(23,26) = pd(23,26) + rrt(765) * density(58) 
  pd(23,58) = pd(23,58) + rrt(765) * density(26) 
  pd(26,26) = pd(26,26) - rrt(765) * density(58) 
  pd(26,58) = pd(26,58) - rrt(765) * density(26) 
  pd(50,26) = pd(50,26) + rrt(765) * density(58) 
  pd(50,58) = pd(50,58) + rrt(765) * density(26) 
  pd(58,26) = pd(58,26) - rrt(765) * density(58) 
  pd(58,58) = pd(58,58) - rrt(765) * density(26) 
  pd(01,27) = pd(01,27) + rrt(766) * density(58) 
  pd(01,58) = pd(01,58) + rrt(766) * density(27) 
  pd(27,27) = pd(27,27) - rrt(766) * density(58) 
  pd(27,58) = pd(27,58) - rrt(766) * density(27) 
  pd(50,27) = pd(50,27) + rrt(766) * density(58) 
  pd(50,58) = pd(50,58) + rrt(766) * density(27) 
  pd(58,27) = pd(58,27) - rrt(766) * density(58) 
  pd(58,58) = pd(58,58) - rrt(766) * density(27) 
  pd(38,42) = pd(38,42) + rrt(767) * density(58) 
  pd(38,58) = pd(38,58) + rrt(767) * density(42) 
  pd(42,42) = pd(42,42) - rrt(767) * density(58) 
  pd(42,58) = pd(42,58) - rrt(767) * density(42) 
  pd(50,42) = pd(50,42) + rrt(767) * density(58) 
  pd(50,58) = pd(50,58) + rrt(767) * density(42) 
  pd(58,42) = pd(58,42) - rrt(767) * density(58) 
  pd(58,58) = pd(58,58) - rrt(767) * density(42) 
  pd(30,43) = pd(30,43) + rrt(768) * density(58) 
  pd(30,58) = pd(30,58) + rrt(768) * density(43) 
  pd(43,43) = pd(43,43) - rrt(768) * density(58) 
  pd(43,58) = pd(43,58) - rrt(768) * density(43) 
  pd(50,43) = pd(50,43) + rrt(768) * density(58) 
  pd(50,58) = pd(50,58) + rrt(768) * density(43) 
  pd(58,43) = pd(58,43) - rrt(768) * density(58) 
  pd(58,58) = pd(58,58) - rrt(768) * density(43) 
  pd(50,55) = pd(50,55) + rrt(769) * density(58) * 2.0d0
  pd(50,58) = pd(50,58) + rrt(769) * density(55) * 2.0d0
  pd(55,55) = pd(55,55) - rrt(769) * density(58) 
  pd(55,58) = pd(55,58) - rrt(769) * density(55) 
  pd(58,55) = pd(58,55) - rrt(769) * density(58) 
  pd(58,58) = pd(58,58) - rrt(769) * density(55) 
  pd(50,56) = pd(50,56) + rrt(770) * density(58) 
  pd(50,58) = pd(50,58) + rrt(770) * density(56) 
  pd(51,56) = pd(51,56) + rrt(770) * density(58) 
  pd(51,58) = pd(51,58) + rrt(770) * density(56) 
  pd(56,56) = pd(56,56) - rrt(770) * density(58) 
  pd(56,58) = pd(56,58) - rrt(770) * density(56) 
  pd(58,56) = pd(58,56) - rrt(770) * density(58) 
  pd(58,58) = pd(58,58) - rrt(770) * density(56) 
  pd(50,57) = pd(50,57) + rrt(771) * density(58) 
  pd(50,58) = pd(50,58) + rrt(771) * density(57) 
  pd(52,57) = pd(52,57) + rrt(771) * density(58) 
  pd(52,58) = pd(52,58) + rrt(771) * density(57) 
  pd(57,57) = pd(57,57) - rrt(771) * density(58) 
  pd(57,58) = pd(57,58) - rrt(771) * density(57) 
  pd(58,57) = pd(58,57) - rrt(771) * density(58) 
  pd(58,58) = pd(58,58) - rrt(771) * density(57) 
  pd(23,26) = pd(23,26) + rrt(772) * density(59) 
  pd(23,59) = pd(23,59) + rrt(772) * density(26) 
  pd(26,26) = pd(26,26) - rrt(772) * density(59) 
  pd(26,59) = pd(26,59) - rrt(772) * density(26) 
  pd(51,26) = pd(51,26) + rrt(772) * density(59) 
  pd(51,59) = pd(51,59) + rrt(772) * density(26) 
  pd(59,26) = pd(59,26) - rrt(772) * density(59) 
  pd(59,59) = pd(59,59) - rrt(772) * density(26) 
  pd(01,27) = pd(01,27) + rrt(773) * density(59) 
  pd(01,59) = pd(01,59) + rrt(773) * density(27) 
  pd(27,27) = pd(27,27) - rrt(773) * density(59) 
  pd(27,59) = pd(27,59) - rrt(773) * density(27) 
  pd(51,27) = pd(51,27) + rrt(773) * density(59) 
  pd(51,59) = pd(51,59) + rrt(773) * density(27) 
  pd(59,27) = pd(59,27) - rrt(773) * density(59) 
  pd(59,59) = pd(59,59) - rrt(773) * density(27) 
  pd(38,42) = pd(38,42) + rrt(774) * density(59) 
  pd(38,59) = pd(38,59) + rrt(774) * density(42) 
  pd(42,42) = pd(42,42) - rrt(774) * density(59) 
  pd(42,59) = pd(42,59) - rrt(774) * density(42) 
  pd(51,42) = pd(51,42) + rrt(774) * density(59) 
  pd(51,59) = pd(51,59) + rrt(774) * density(42) 
  pd(59,42) = pd(59,42) - rrt(774) * density(59) 
  pd(59,59) = pd(59,59) - rrt(774) * density(42) 
  pd(30,43) = pd(30,43) + rrt(775) * density(59) 
  pd(30,59) = pd(30,59) + rrt(775) * density(43) 
  pd(43,43) = pd(43,43) - rrt(775) * density(59) 
  pd(43,59) = pd(43,59) - rrt(775) * density(43) 
  pd(51,43) = pd(51,43) + rrt(775) * density(59) 
  pd(51,59) = pd(51,59) + rrt(775) * density(43) 
  pd(59,43) = pd(59,43) - rrt(775) * density(59) 
  pd(59,59) = pd(59,59) - rrt(775) * density(43) 
  pd(50,55) = pd(50,55) + rrt(776) * density(59) 
  pd(50,59) = pd(50,59) + rrt(776) * density(55) 
  pd(51,55) = pd(51,55) + rrt(776) * density(59) 
  pd(51,59) = pd(51,59) + rrt(776) * density(55) 
  pd(55,55) = pd(55,55) - rrt(776) * density(59) 
  pd(55,59) = pd(55,59) - rrt(776) * density(55) 
  pd(59,55) = pd(59,55) - rrt(776) * density(59) 
  pd(59,59) = pd(59,59) - rrt(776) * density(55) 
  pd(51,56) = pd(51,56) + rrt(777) * density(59) * 2.0d0
  pd(51,59) = pd(51,59) + rrt(777) * density(56) * 2.0d0
  pd(56,56) = pd(56,56) - rrt(777) * density(59) 
  pd(56,59) = pd(56,59) - rrt(777) * density(56) 
  pd(59,56) = pd(59,56) - rrt(777) * density(59) 
  pd(59,59) = pd(59,59) - rrt(777) * density(56) 
  pd(51,57) = pd(51,57) + rrt(778) * density(59) 
  pd(51,59) = pd(51,59) + rrt(778) * density(57) 
  pd(52,57) = pd(52,57) + rrt(778) * density(59) 
  pd(52,59) = pd(52,59) + rrt(778) * density(57) 
  pd(57,57) = pd(57,57) - rrt(778) * density(59) 
  pd(57,59) = pd(57,59) - rrt(778) * density(57) 
  pd(59,57) = pd(59,57) - rrt(778) * density(59) 
  pd(59,59) = pd(59,59) - rrt(778) * density(57) 
  pd(23,26) = pd(23,26) + rrt(779) * density(60) 
  pd(23,60) = pd(23,60) + rrt(779) * density(26) 
  pd(26,26) = pd(26,26) - rrt(779) * density(60) 
  pd(26,60) = pd(26,60) - rrt(779) * density(26) 
  pd(52,26) = pd(52,26) + rrt(779) * density(60) 
  pd(52,60) = pd(52,60) + rrt(779) * density(26) 
  pd(60,26) = pd(60,26) - rrt(779) * density(60) 
  pd(60,60) = pd(60,60) - rrt(779) * density(26) 
  pd(01,27) = pd(01,27) + rrt(780) * density(60) 
  pd(01,60) = pd(01,60) + rrt(780) * density(27) 
  pd(27,27) = pd(27,27) - rrt(780) * density(60) 
  pd(27,60) = pd(27,60) - rrt(780) * density(27) 
  pd(52,27) = pd(52,27) + rrt(780) * density(60) 
  pd(52,60) = pd(52,60) + rrt(780) * density(27) 
  pd(60,27) = pd(60,27) - rrt(780) * density(60) 
  pd(60,60) = pd(60,60) - rrt(780) * density(27) 
  pd(38,42) = pd(38,42) + rrt(781) * density(60) 
  pd(38,60) = pd(38,60) + rrt(781) * density(42) 
  pd(42,42) = pd(42,42) - rrt(781) * density(60) 
  pd(42,60) = pd(42,60) - rrt(781) * density(42) 
  pd(52,42) = pd(52,42) + rrt(781) * density(60) 
  pd(52,60) = pd(52,60) + rrt(781) * density(42) 
  pd(60,42) = pd(60,42) - rrt(781) * density(60) 
  pd(60,60) = pd(60,60) - rrt(781) * density(42) 
  pd(30,43) = pd(30,43) + rrt(782) * density(60) 
  pd(30,60) = pd(30,60) + rrt(782) * density(43) 
  pd(43,43) = pd(43,43) - rrt(782) * density(60) 
  pd(43,60) = pd(43,60) - rrt(782) * density(43) 
  pd(52,43) = pd(52,43) + rrt(782) * density(60) 
  pd(52,60) = pd(52,60) + rrt(782) * density(43) 
  pd(60,43) = pd(60,43) - rrt(782) * density(60) 
  pd(60,60) = pd(60,60) - rrt(782) * density(43) 
  pd(50,55) = pd(50,55) + rrt(783) * density(60) 
  pd(50,60) = pd(50,60) + rrt(783) * density(55) 
  pd(52,55) = pd(52,55) + rrt(783) * density(60) 
  pd(52,60) = pd(52,60) + rrt(783) * density(55) 
  pd(55,55) = pd(55,55) - rrt(783) * density(60) 
  pd(55,60) = pd(55,60) - rrt(783) * density(55) 
  pd(60,55) = pd(60,55) - rrt(783) * density(60) 
  pd(60,60) = pd(60,60) - rrt(783) * density(55) 
  pd(51,56) = pd(51,56) + rrt(784) * density(60) 
  pd(51,60) = pd(51,60) + rrt(784) * density(56) 
  pd(52,56) = pd(52,56) + rrt(784) * density(60) 
  pd(52,60) = pd(52,60) + rrt(784) * density(56) 
  pd(56,56) = pd(56,56) - rrt(784) * density(60) 
  pd(56,60) = pd(56,60) - rrt(784) * density(56) 
  pd(60,56) = pd(60,56) - rrt(784) * density(60) 
  pd(60,60) = pd(60,60) - rrt(784) * density(56) 
  pd(52,57) = pd(52,57) + rrt(785) * density(60) * 2.0d0
  pd(52,60) = pd(52,60) + rrt(785) * density(57) * 2.0d0
  pd(57,57) = pd(57,57) - rrt(785) * density(60) 
  pd(57,60) = pd(57,60) - rrt(785) * density(57) 
  pd(60,57) = pd(60,57) - rrt(785) * density(60) 
  pd(60,60) = pd(60,60) - rrt(785) * density(57) 
  pd(23,26) = pd(23,26) + rrt(786) * density(61) 
  pd(23,61) = pd(23,61) + rrt(786) * density(26) 
  pd(26,26) = pd(26,26) - rrt(786) * density(61) 
  pd(26,61) = pd(26,61) - rrt(786) * density(26) 
  pd(53,26) = pd(53,26) + rrt(786) * density(61) 
  pd(53,61) = pd(53,61) + rrt(786) * density(26) 
  pd(61,26) = pd(61,26) - rrt(786) * density(61) 
  pd(61,61) = pd(61,61) - rrt(786) * density(26) 
  pd(01,27) = pd(01,27) + rrt(787) * density(61) 
  pd(01,61) = pd(01,61) + rrt(787) * density(27) 
  pd(27,27) = pd(27,27) - rrt(787) * density(61) 
  pd(27,61) = pd(27,61) - rrt(787) * density(27) 
  pd(53,27) = pd(53,27) + rrt(787) * density(61) 
  pd(53,61) = pd(53,61) + rrt(787) * density(27) 
  pd(61,27) = pd(61,27) - rrt(787) * density(61) 
  pd(61,61) = pd(61,61) - rrt(787) * density(27) 
  pd(38,42) = pd(38,42) + rrt(788) * density(61) 
  pd(38,61) = pd(38,61) + rrt(788) * density(42) 
  pd(42,42) = pd(42,42) - rrt(788) * density(61) 
  pd(42,61) = pd(42,61) - rrt(788) * density(42) 
  pd(53,42) = pd(53,42) + rrt(788) * density(61) 
  pd(53,61) = pd(53,61) + rrt(788) * density(42) 
  pd(61,42) = pd(61,42) - rrt(788) * density(61) 
  pd(61,61) = pd(61,61) - rrt(788) * density(42) 
  pd(30,43) = pd(30,43) + rrt(789) * density(61) 
  pd(30,61) = pd(30,61) + rrt(789) * density(43) 
  pd(43,43) = pd(43,43) - rrt(789) * density(61) 
  pd(43,61) = pd(43,61) - rrt(789) * density(43) 
  pd(53,43) = pd(53,43) + rrt(789) * density(61) 
  pd(53,61) = pd(53,61) + rrt(789) * density(43) 
  pd(61,43) = pd(61,43) - rrt(789) * density(61) 
  pd(61,61) = pd(61,61) - rrt(789) * density(43) 
  pd(50,55) = pd(50,55) + rrt(790) * density(61) 
  pd(50,61) = pd(50,61) + rrt(790) * density(55) 
  pd(53,55) = pd(53,55) + rrt(790) * density(61) 
  pd(53,61) = pd(53,61) + rrt(790) * density(55) 
  pd(55,55) = pd(55,55) - rrt(790) * density(61) 
  pd(55,61) = pd(55,61) - rrt(790) * density(55) 
  pd(61,55) = pd(61,55) - rrt(790) * density(61) 
  pd(61,61) = pd(61,61) - rrt(790) * density(55) 
  pd(51,56) = pd(51,56) + rrt(791) * density(61) 
  pd(51,61) = pd(51,61) + rrt(791) * density(56) 
  pd(53,56) = pd(53,56) + rrt(791) * density(61) 
  pd(53,61) = pd(53,61) + rrt(791) * density(56) 
  pd(56,56) = pd(56,56) - rrt(791) * density(61) 
  pd(56,61) = pd(56,61) - rrt(791) * density(56) 
  pd(61,56) = pd(61,56) - rrt(791) * density(61) 
  pd(61,61) = pd(61,61) - rrt(791) * density(56) 
  pd(52,57) = pd(52,57) + rrt(792) * density(61) 
  pd(52,61) = pd(52,61) + rrt(792) * density(57) 
  pd(53,57) = pd(53,57) + rrt(792) * density(61) 
  pd(53,61) = pd(53,61) + rrt(792) * density(57) 
  pd(57,57) = pd(57,57) - rrt(792) * density(61) 
  pd(57,61) = pd(57,61) - rrt(792) * density(57) 
  pd(61,57) = pd(61,57) - rrt(792) * density(61) 
  pd(61,61) = pd(61,61) - rrt(792) * density(57) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  pd(80,:) = 0.0d0
  if( lgas_heating ) then
    ysum = sum(density(1:species_max)) - density(species_electrons)
    rrt(050) = rrt(050) * 3.377D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(050) * density(02) 
    pd(80,02) = pd(80,02) + rrt(050) * density(01) 
    rrt(051) = rrt(051) * 3.470D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(051) * density(03) 
    pd(80,03) = pd(80,03) + rrt(051) * density(01) 
    rrt(052) = rrt(052) * 3.365D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(052) * density(04) 
    pd(80,04) = pd(80,04) + rrt(052) * density(01) 
    rrt(053) = rrt(053) * 3.365D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(053) * density(05) 
    pd(80,05) = pd(80,05) + rrt(053) * density(01) 
    rrt(054) = rrt(054) * 3.481D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(054) * density(06) 
    pd(80,06) = pd(80,06) + rrt(054) * density(01) 
    rrt(055) = rrt(055) * 3.365D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(055) * density(07) 
    pd(80,07) = pd(80,07) + rrt(055) * density(01) 
    rrt(056) = rrt(056) * 3.481D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(056) * density(08) 
    pd(80,08) = pd(80,08) + rrt(056) * density(01) 
    rrt(057) = rrt(057) * 3.365D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(057) * density(09) 
    pd(80,09) = pd(80,09) + rrt(057) * density(01) 
    rrt(058) = rrt(058) * 3.376D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(058) * density(01) 
    rrt(059) = rrt(059) * 3.468D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(059) * density(02) 
    pd(80,02) = pd(80,02) - rrt(059) * density(01) 
    rrt(060) = rrt(060) * 3.364D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(060) * density(03) 
    pd(80,03) = pd(80,03) - rrt(060) * density(01) 
    rrt(061) = rrt(061) * 3.364D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(061) * density(04) 
    pd(80,04) = pd(80,04) - rrt(061) * density(01) 
    rrt(062) = rrt(062) * 3.480D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(062) * density(05) 
    pd(80,05) = pd(80,05) - rrt(062) * density(01) 
    rrt(063) = rrt(063) * 3.364D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(063) * density(06) 
    pd(80,06) = pd(80,06) - rrt(063) * density(01) 
    rrt(064) = rrt(064) * 3.480D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(064) * density(07) 
    pd(80,07) = pd(80,07) - rrt(064) * density(01) 
    rrt(065) = rrt(065) * 3.364D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(065) * density(08) 
    pd(80,08) = pd(80,08) - rrt(065) * density(01) 
    rrt(066) = rrt(066) * 3.377D+03 / ysum
    pd(80,02) = pd(80,02) + rrt(066) * density(23) 
    pd(80,23) = pd(80,23) + rrt(066) * density(02) 
    rrt(067) = rrt(067) * 3.470D+03 / ysum
    pd(80,03) = pd(80,03) + rrt(067) * density(23) 
    pd(80,23) = pd(80,23) + rrt(067) * density(03) 
    rrt(068) = rrt(068) * 3.365D+03 / ysum
    pd(80,04) = pd(80,04) + rrt(068) * density(23) 
    pd(80,23) = pd(80,23) + rrt(068) * density(04) 
    rrt(069) = rrt(069) * 3.365D+03 / ysum
    pd(80,05) = pd(80,05) + rrt(069) * density(23) 
    pd(80,23) = pd(80,23) + rrt(069) * density(05) 
    rrt(070) = rrt(070) * 3.481D+03 / ysum
    pd(80,06) = pd(80,06) + rrt(070) * density(23) 
    pd(80,23) = pd(80,23) + rrt(070) * density(06) 
    rrt(071) = rrt(071) * 3.365D+03 / ysum
    pd(80,07) = pd(80,07) + rrt(071) * density(23) 
    pd(80,23) = pd(80,23) + rrt(071) * density(07) 
    rrt(072) = rrt(072) * 3.481D+03 / ysum
    pd(80,08) = pd(80,08) + rrt(072) * density(23) 
    pd(80,23) = pd(80,23) + rrt(072) * density(08) 
    rrt(073) = rrt(073) * 3.365D+03 / ysum
    pd(80,09) = pd(80,09) + rrt(073) * density(23) 
    pd(80,23) = pd(80,23) + rrt(073) * density(09) 
    rrt(074) = rrt(074) * 3.376D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(074) * density(23) 
    pd(80,23) = pd(80,23) - rrt(074) * density(01) 
    rrt(075) = rrt(075) * 3.468D+03 / ysum
    pd(80,02) = pd(80,02) - rrt(075) * density(23) 
    pd(80,23) = pd(80,23) - rrt(075) * density(02) 
    rrt(076) = rrt(076) * 3.364D+03 / ysum
    pd(80,03) = pd(80,03) - rrt(076) * density(23) 
    pd(80,23) = pd(80,23) - rrt(076) * density(03) 
    rrt(077) = rrt(077) * 3.364D+03 / ysum
    pd(80,04) = pd(80,04) - rrt(077) * density(23) 
    pd(80,23) = pd(80,23) - rrt(077) * density(04) 
    rrt(078) = rrt(078) * 3.480D+03 / ysum
    pd(80,05) = pd(80,05) - rrt(078) * density(23) 
    pd(80,23) = pd(80,23) - rrt(078) * density(05) 
    rrt(079) = rrt(079) * 3.364D+03 / ysum
    pd(80,06) = pd(80,06) - rrt(079) * density(23) 
    pd(80,23) = pd(80,23) - rrt(079) * density(06) 
    rrt(080) = rrt(080) * 3.480D+03 / ysum
    pd(80,07) = pd(80,07) - rrt(080) * density(23) 
    pd(80,23) = pd(80,23) - rrt(080) * density(07) 
    rrt(081) = rrt(081) * 3.364D+03 / ysum
    pd(80,08) = pd(80,08) - rrt(081) * density(23) 
    pd(80,23) = pd(80,23) - rrt(081) * density(08) 
    rrt(082) = rrt(082) * 3.377D+03 / ysum
    pd(80,02) = pd(80,02) + rrt(082) * density(38) 
    pd(80,38) = pd(80,38) + rrt(082) * density(02) 
    rrt(083) = rrt(083) * 3.470D+03 / ysum
    pd(80,03) = pd(80,03) + rrt(083) * density(38) 
    pd(80,38) = pd(80,38) + rrt(083) * density(03) 
    rrt(084) = rrt(084) * 3.365D+03 / ysum
    pd(80,04) = pd(80,04) + rrt(084) * density(38) 
    pd(80,38) = pd(80,38) + rrt(084) * density(04) 
    rrt(085) = rrt(085) * 3.365D+03 / ysum
    pd(80,05) = pd(80,05) + rrt(085) * density(38) 
    pd(80,38) = pd(80,38) + rrt(085) * density(05) 
    rrt(086) = rrt(086) * 3.481D+03 / ysum
    pd(80,06) = pd(80,06) + rrt(086) * density(38) 
    pd(80,38) = pd(80,38) + rrt(086) * density(06) 
    rrt(087) = rrt(087) * 3.365D+03 / ysum
    pd(80,07) = pd(80,07) + rrt(087) * density(38) 
    pd(80,38) = pd(80,38) + rrt(087) * density(07) 
    rrt(088) = rrt(088) * 3.481D+03 / ysum
    pd(80,08) = pd(80,08) + rrt(088) * density(38) 
    pd(80,38) = pd(80,38) + rrt(088) * density(08) 
    rrt(089) = rrt(089) * 3.365D+03 / ysum
    pd(80,09) = pd(80,09) + rrt(089) * density(38) 
    pd(80,38) = pd(80,38) + rrt(089) * density(09) 
    rrt(090) = rrt(090) * 3.376D+03 / ysum
    pd(80,01) = pd(80,01) - rrt(090) * density(38) 
    pd(80,38) = pd(80,38) - rrt(090) * density(01) 
    rrt(091) = rrt(091) * 3.468D+03 / ysum
    pd(80,02) = pd(80,02) - rrt(091) * density(38) 
    pd(80,38) = pd(80,38) - rrt(091) * density(02) 
    rrt(092) = rrt(092) * 3.364D+03 / ysum
    pd(80,03) = pd(80,03) - rrt(092) * density(38) 
    pd(80,38) = pd(80,38) - rrt(092) * density(03) 
    rrt(093) = rrt(093) * 3.364D+03 / ysum
    pd(80,04) = pd(80,04) - rrt(093) * density(38) 
    pd(80,38) = pd(80,38) - rrt(093) * density(04) 
    rrt(094) = rrt(094) * 3.480D+03 / ysum
    pd(80,05) = pd(80,05) - rrt(094) * density(38) 
    pd(80,38) = pd(80,38) - rrt(094) * density(05) 
    rrt(095) = rrt(095) * 3.364D+03 / ysum
    pd(80,06) = pd(80,06) - rrt(095) * density(38) 
    pd(80,38) = pd(80,38) - rrt(095) * density(06) 
    rrt(096) = rrt(096) * 3.480D+03 / ysum
    pd(80,07) = pd(80,07) - rrt(096) * density(38) 
    pd(80,38) = pd(80,38) - rrt(096) * density(07) 
    rrt(097) = rrt(097) * 3.364D+03 / ysum
    pd(80,08) = pd(80,08) - rrt(097) * density(38) 
    pd(80,38) = pd(80,38) - rrt(097) * density(08) 
    rrt(098) = rrt(098) * 2.205D+03 / ysum
    pd(80,30) = pd(80,30) + rrt(098) * density(31) 
    pd(80,31) = pd(80,31) + rrt(098) * density(30) 
    rrt(099) = rrt(099) * 2.205D+03 / ysum
    pd(80,30) = pd(80,30) + rrt(099) * density(32) 
    pd(80,32) = pd(80,32) + rrt(099) * density(30) 
    rrt(100) = rrt(100) * 2.205D+03 / ysum
    pd(80,30) = pd(80,30) + rrt(100) * density(33) 
    pd(80,33) = pd(80,33) + rrt(100) * density(30) 
    rrt(101) = rrt(101) * 2.089D+03 / ysum
    pd(80,30) = pd(80,30) + rrt(101) * density(34) 
    pd(80,34) = pd(80,34) + rrt(101) * density(30) 
    rrt(102) = rrt(102) * 2.204D+03 / ysum
    pd(80,30) = pd(80,30) - rrt(102) * density(30) 
    rrt(103) = rrt(103) * 2.204D+03 / ysum
    pd(80,30) = pd(80,30) - rrt(103) * density(31) 
    pd(80,31) = pd(80,31) - rrt(103) * density(30) 
    rrt(104) = rrt(104) * 2.204D+03 / ysum
    pd(80,30) = pd(80,30) - rrt(104) * density(32) 
    pd(80,32) = pd(80,32) - rrt(104) * density(30) 
    rrt(105) = rrt(105) * 2.088D+03 / ysum
    pd(80,30) = pd(80,30) - rrt(105) * density(33) 
    pd(80,33) = pd(80,33) - rrt(105) * density(30) 
    rrt(106) = rrt(106) * 2.205D+03 / ysum
    pd(80,31) = pd(80,31) + rrt(106) * density(38) 
    pd(80,38) = pd(80,38) + rrt(106) * density(31) 
    rrt(107) = rrt(107) * 2.205D+03 / ysum
    pd(80,32) = pd(80,32) + rrt(107) * density(38) 
    pd(80,38) = pd(80,38) + rrt(107) * density(32) 
    rrt(108) = rrt(108) * 2.205D+03 / ysum
    pd(80,33) = pd(80,33) + rrt(108) * density(38) 
    pd(80,38) = pd(80,38) + rrt(108) * density(33) 
    rrt(109) = rrt(109) * 2.089D+03 / ysum
    pd(80,34) = pd(80,34) + rrt(109) * density(38) 
    pd(80,38) = pd(80,38) + rrt(109) * density(34) 
    rrt(110) = rrt(110) * 2.204D+03 / ysum
    pd(80,30) = pd(80,30) - rrt(110) * density(38) 
    pd(80,38) = pd(80,38) - rrt(110) * density(30) 
    rrt(111) = rrt(111) * 2.204D+03 / ysum
    pd(80,31) = pd(80,31) - rrt(111) * density(38) 
    pd(80,38) = pd(80,38) - rrt(111) * density(31) 
    rrt(112) = rrt(112) * 2.204D+03 / ysum
    pd(80,32) = pd(80,32) - rrt(112) * density(38) 
    pd(80,38) = pd(80,38) - rrt(112) * density(32) 
    rrt(113) = rrt(113) * 2.088D+03 / ysum
    pd(80,33) = pd(80,33) - rrt(113) * density(38) 
    pd(80,38) = pd(80,38) - rrt(113) * density(33) 
    rrt(142) = rrt(142) * 9.957D+03 / ysum
    pd(80,01) = pd(80,01) + rrt(142) * density(79) 
    pd(80,79) = pd(80,79) + rrt(142) * density(01) 
    rrt(146) = rrt(146) * 1.026D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(146) * density(79) 
    pd(80,79) = pd(80,79) + rrt(146) * density(30) 
    rrt(147) = rrt(147) * 1.188D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(147) * density(79) 
    pd(80,79) = pd(80,79) + rrt(147) * density(30) 
    rrt(148) = rrt(148) * 1.061D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(148) * density(79) 
    pd(80,79) = pd(80,79) + rrt(148) * density(30) 
    rrt(169) = rrt(169) * 8.065D+04 / ysum
    pd(80,43) = pd(80,43) + rrt(169) * density(79) 
    pd(80,79) = pd(80,79) + rrt(169) * density(43) 
    rrt(170) = rrt(170) * 5.791D+04 / ysum
    pd(80,43) = pd(80,43) + rrt(170) * density(79) 
    pd(80,79) = pd(80,79) + rrt(170) * density(43) 
    rrt(243) = rrt(243) * 1.218D+04 / ysum
    pd(80,17) = pd(80,17) + rrt(243) * density(30) 
    pd(80,30) = pd(80,30) + rrt(243) * density(17) 
    rrt(244) = rrt(244) * 5.268D+04 / ysum
    pd(80,17) = pd(80,17) + rrt(244) * density(30) 
    pd(80,30) = pd(80,30) + rrt(244) * density(17) 
    rrt(245) = rrt(245) * 2.588D+04 / ysum
    pd(80,18) = pd(80,18) + rrt(245) * density(30) 
    pd(80,30) = pd(80,30) + rrt(245) * density(18) 
    rrt(246) = rrt(246) * 1.184D+04 / ysum
    pd(80,19) = pd(80,19) + rrt(246) * density(30) 
    pd(80,30) = pd(80,30) + rrt(246) * density(19) 
    rrt(247) = rrt(247) * 6.858D+04 / ysum
    pd(80,20) = pd(80,20) + rrt(247) * density(30) 
    pd(80,30) = pd(80,30) + rrt(247) * density(20) 
    rrt(248) = rrt(248) * 4.236D+04 / ysum
    pd(80,20) = pd(80,20) + rrt(248) * density(30) 
    pd(80,30) = pd(80,30) + rrt(248) * density(20) 
    rrt(249) = rrt(249) * 2.286D+04 / ysum
    pd(80,20) = pd(80,20) + rrt(249) * density(30) 
    pd(80,30) = pd(80,30) + rrt(249) * density(20) 
    rrt(251) = rrt(251) * 2.588D+04 / ysum
    pd(80,17) = pd(80,17) + rrt(251) * density(38) 
    pd(80,38) = pd(80,38) + rrt(251) * density(17) 
    rrt(258) = rrt(258) * 4.642D+04 / ysum
    pd(80,17) = pd(80,17) + rrt(258) * density(17) 
    rrt(259) = rrt(259) * 1.520D+04 / ysum
    pd(80,17) = pd(80,17) + rrt(259) * density(17) 
    rrt(316) = rrt(316) * 2.623D+04 / ysum
    pd(80,38) = pd(80,38) + rrt(316) * density(39) 
    pd(80,39) = pd(80,39) + rrt(316) * density(38) 
    rrt(317) = rrt(317) * 2.623D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(317) * density(39) 
    pd(80,39) = pd(80,39) + rrt(317) * density(30) 
    rrt(318) = rrt(318) * 1.485D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(318) * density(39) 
    pd(80,39) = pd(80,39) + rrt(318) * density(30) 
    rrt(319) = rrt(319) * 7.311D+03 / ysum
    pd(80,30) = pd(80,30) + rrt(319) * density(39) 
    pd(80,39) = pd(80,39) + rrt(319) * density(30) 
    rrt(320) = rrt(320) * 2.623D+04 / ysum
    pd(80,01) = pd(80,01) + rrt(320) * density(39) 
    pd(80,39) = pd(80,39) + rrt(320) * density(01) 
    rrt(326) = rrt(326) * 1.950D+04 / ysum
    pd(80,38) = pd(80,38) + rrt(326) * density(40) 
    pd(80,40) = pd(80,40) + rrt(326) * density(38) 
    rrt(327) = rrt(327) * 4.572D+04 / ysum
    pd(80,23) = pd(80,23) + rrt(327) * density(40) 
    pd(80,40) = pd(80,40) + rrt(327) * density(23) 
    rrt(328) = rrt(328) * 1.950D+04 / ysum
    pd(80,30) = pd(80,30) + rrt(328) * density(40) 
    pd(80,40) = pd(80,40) + rrt(328) * density(30) 
    rrt(329) = rrt(329) * 1.369D+04 / ysum
    pd(80,30) = pd(80,30) - rrt(329) * density(40) 
    pd(80,40) = pd(80,40) - rrt(329) * density(30) 
    rrt(330) = rrt(330) * 4.572D+04 / ysum
    pd(80,01) = pd(80,01) + rrt(330) * density(40) 
    pd(80,40) = pd(80,40) + rrt(330) * density(01) 
    rrt(331) = rrt(331) * 4.874D+03 / ysum
    pd(80,35) = pd(80,35) + rrt(331) * density(40) 
    pd(80,40) = pd(80,40) + rrt(331) * density(35) 
    rrt(332) = rrt(332) * 1.195D+04 / ysum
    pd(80,35) = pd(80,35) + rrt(332) * density(40) 
    pd(80,40) = pd(80,40) + rrt(332) * density(35) 
    rrt(333) = rrt(333) * 2.320D+03 / ysum
    pd(80,35) = pd(80,35) - rrt(333) * density(40) 
    pd(80,40) = pd(80,40) - rrt(333) * density(35) 
    rrt(334) = rrt(334) * 5.709D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(334) * density(50) 
    pd(80,50) = pd(80,50) + rrt(334) * density(40) 
    rrt(335) = rrt(335) * 1.950D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(335) * density(50) 
    pd(80,50) = pd(80,50) + rrt(335) * density(40) 
    rrt(336) = rrt(336) * 4.572D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(336) * density(41) 
    pd(80,41) = pd(80,41) + rrt(336) * density(40) 
    rrt(337) = rrt(337) * 4.572D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(337) * density(41) 
    pd(80,41) = pd(80,41) + rrt(337) * density(40) 
    rrt(338) = rrt(338) * 4.572D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(338) * density(51) 
    pd(80,51) = pd(80,51) + rrt(338) * density(40) 
    rrt(339) = rrt(339) * 1.950D+04 / ysum
    pd(80,40) = pd(80,40) + rrt(339) * density(51) 
    pd(80,51) = pd(80,51) + rrt(339) * density(40) 
    rrt(589) = rrt(589) * 1.519D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(589) * density(46) 
    pd(80,46) = pd(80,46) + rrt(589) * density(26) 
    rrt(590) = rrt(590) * 1.641D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(590) * density(46) 
    pd(80,46) = pd(80,46) + rrt(590) * density(27) 
    rrt(591) = rrt(591) * 1.413D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(591) * density(46) 
    pd(80,46) = pd(80,46) + rrt(591) * density(42) 
    rrt(592) = rrt(592) * 1.234D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(592) * density(46) 
    pd(80,46) = pd(80,46) + rrt(592) * density(43) 
    rrt(593) = rrt(593) * 9.081D+04 / ysum
    pd(80,46) = pd(80,46) + rrt(593) * density(55) 
    pd(80,55) = pd(80,55) + rrt(593) * density(46) 
    rrt(594) = rrt(594) * 1.329D+05 / ysum
    pd(80,46) = pd(80,46) + rrt(594) * density(56) 
    pd(80,56) = pd(80,56) + rrt(594) * density(46) 
    rrt(595) = rrt(595) * 9.454D+04 / ysum
    pd(80,46) = pd(80,46) + rrt(595) * density(57) 
    pd(80,57) = pd(80,57) + rrt(595) * density(46) 
    rrt(596) = rrt(596) * 1.635D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(596) * density(47) 
    pd(80,47) = pd(80,47) + rrt(596) * density(26) 
    rrt(597) = rrt(597) * 1.756D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(597) * density(47) 
    pd(80,47) = pd(80,47) + rrt(597) * density(27) 
    rrt(598) = rrt(598) * 1.528D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(598) * density(47) 
    pd(80,47) = pd(80,47) + rrt(598) * density(42) 
    rrt(599) = rrt(599) * 1.348D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(599) * density(47) 
    pd(80,47) = pd(80,47) + rrt(599) * density(43) 
    rrt(600) = rrt(600) * 1.023D+05 / ysum
    pd(80,47) = pd(80,47) + rrt(600) * density(55) 
    pd(80,55) = pd(80,55) + rrt(600) * density(47) 
    rrt(601) = rrt(601) * 1.444D+05 / ysum
    pd(80,47) = pd(80,47) + rrt(601) * density(56) 
    pd(80,56) = pd(80,56) + rrt(601) * density(47) 
    rrt(602) = rrt(602) * 1.060D+05 / ysum
    pd(80,47) = pd(80,47) + rrt(602) * density(57) 
    pd(80,57) = pd(80,57) + rrt(602) * density(47) 
    rrt(603) = rrt(603) * 1.442D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(603) * density(48) 
    pd(80,48) = pd(80,48) + rrt(603) * density(26) 
    rrt(604) = rrt(604) * 1.564D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(604) * density(48) 
    pd(80,48) = pd(80,48) + rrt(604) * density(27) 
    rrt(605) = rrt(605) * 1.337D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(605) * density(48) 
    pd(80,48) = pd(80,48) + rrt(605) * density(42) 
    rrt(606) = rrt(606) * 1.157D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(606) * density(48) 
    pd(80,48) = pd(80,48) + rrt(606) * density(43) 
    rrt(607) = rrt(607) * 8.310D+04 / ysum
    pd(80,48) = pd(80,48) + rrt(607) * density(55) 
    pd(80,55) = pd(80,55) + rrt(607) * density(48) 
    rrt(608) = rrt(608) * 1.252D+05 / ysum
    pd(80,48) = pd(80,48) + rrt(608) * density(56) 
    pd(80,56) = pd(80,56) + rrt(608) * density(48) 
    rrt(609) = rrt(609) * 8.684D+04 / ysum
    pd(80,48) = pd(80,48) + rrt(609) * density(57) 
    pd(80,57) = pd(80,57) + rrt(609) * density(48) 
    rrt(617) = rrt(617) * 1.661D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(617) * density(59) 
    pd(80,59) = pd(80,59) + rrt(617) * density(26) 
    rrt(618) = rrt(618) * 1.782D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(618) * density(59) 
    pd(80,59) = pd(80,59) + rrt(618) * density(27) 
    rrt(619) = rrt(619) * 1.555D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(619) * density(59) 
    pd(80,59) = pd(80,59) + rrt(619) * density(42) 
    rrt(620) = rrt(620) * 1.375D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(620) * density(59) 
    pd(80,59) = pd(80,59) + rrt(620) * density(43) 
    rrt(621) = rrt(621) * 1.050D+05 / ysum
    pd(80,55) = pd(80,55) + rrt(621) * density(59) 
    pd(80,59) = pd(80,59) + rrt(621) * density(55) 
    rrt(622) = rrt(622) * 1.470D+05 / ysum
    pd(80,56) = pd(80,56) + rrt(622) * density(59) 
    pd(80,59) = pd(80,59) + rrt(622) * density(56) 
    rrt(623) = rrt(623) * 1.087D+05 / ysum
    pd(80,57) = pd(80,57) + rrt(623) * density(59) 
    pd(80,59) = pd(80,59) + rrt(623) * density(57) 
    rrt(624) = rrt(624) * 1.423D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(624) * density(60) 
    pd(80,60) = pd(80,60) + rrt(624) * density(26) 
    rrt(625) = rrt(625) * 1.545D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(625) * density(60) 
    pd(80,60) = pd(80,60) + rrt(625) * density(27) 
    rrt(626) = rrt(626) * 1.317D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(626) * density(60) 
    pd(80,60) = pd(80,60) + rrt(626) * density(42) 
    rrt(627) = rrt(627) * 1.137D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(627) * density(60) 
    pd(80,60) = pd(80,60) + rrt(627) * density(43) 
    rrt(628) = rrt(628) * 8.113D+04 / ysum
    pd(80,55) = pd(80,55) + rrt(628) * density(60) 
    pd(80,60) = pd(80,60) + rrt(628) * density(55) 
    rrt(629) = rrt(629) * 1.232D+05 / ysum
    pd(80,56) = pd(80,56) + rrt(629) * density(60) 
    pd(80,60) = pd(80,60) + rrt(629) * density(56) 
    rrt(630) = rrt(630) * 8.486D+04 / ysum
    pd(80,57) = pd(80,57) + rrt(630) * density(60) 
    pd(80,60) = pd(80,60) + rrt(630) * density(57) 
    rrt(631) = rrt(631) * 1.230D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(631) * density(61) 
    pd(80,61) = pd(80,61) + rrt(631) * density(26) 
    rrt(632) = rrt(632) * 1.351D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(632) * density(61) 
    pd(80,61) = pd(80,61) + rrt(632) * density(27) 
    rrt(633) = rrt(633) * 1.123D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(633) * density(61) 
    pd(80,61) = pd(80,61) + rrt(633) * density(42) 
    rrt(634) = rrt(634) * 9.438D+04 / ysum
    pd(80,43) = pd(80,43) + rrt(634) * density(61) 
    pd(80,61) = pd(80,61) + rrt(634) * density(43) 
    rrt(635) = rrt(635) * 6.182D+04 / ysum
    pd(80,55) = pd(80,55) + rrt(635) * density(61) 
    pd(80,61) = pd(80,61) + rrt(635) * density(55) 
    rrt(636) = rrt(636) * 1.039D+05 / ysum
    pd(80,56) = pd(80,56) + rrt(636) * density(61) 
    pd(80,61) = pd(80,61) + rrt(636) * density(56) 
    rrt(637) = rrt(637) * 6.555D+04 / ysum
    pd(80,57) = pd(80,57) + rrt(637) * density(61) 
    pd(80,61) = pd(80,61) + rrt(637) * density(57) 
    rrt(712) = rrt(712) * 1.519D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(712) * density(46) 
    pd(80,46) = pd(80,46) + rrt(712) * density(26) 
    rrt(713) = rrt(713) * 1.641D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(713) * density(46) 
    pd(80,46) = pd(80,46) + rrt(713) * density(27) 
    rrt(714) = rrt(714) * 1.413D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(714) * density(46) 
    pd(80,46) = pd(80,46) + rrt(714) * density(42) 
    rrt(715) = rrt(715) * 1.234D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(715) * density(46) 
    pd(80,46) = pd(80,46) + rrt(715) * density(43) 
    rrt(716) = rrt(716) * 9.081D+04 / ysum
    pd(80,46) = pd(80,46) + rrt(716) * density(55) 
    pd(80,55) = pd(80,55) + rrt(716) * density(46) 
    rrt(717) = rrt(717) * 1.635D+05 / ysum
    pd(80,26) = pd(80,26) + rrt(717) * density(47) 
    pd(80,47) = pd(80,47) + rrt(717) * density(26) 
    rrt(718) = rrt(718) * 1.756D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(718) * density(47) 
    pd(80,47) = pd(80,47) + rrt(718) * density(27) 
    rrt(719) = rrt(719) * 1.528D+05 / ysum
    pd(80,42) = pd(80,42) + rrt(719) * density(47) 
    pd(80,47) = pd(80,47) + rrt(719) * density(42) 
    rrt(720) = rrt(720) * 1.348D+05 / ysum
    pd(80,43) = pd(80,43) + rrt(720) * density(47) 
    pd(80,47) = pd(80,47) + rrt(720) * density(43) 
    rrt(721) = rrt(721) * 1.023D+05 / ysum
    pd(80,47) = pd(80,47) + rrt(721) * density(55) 
    pd(80,55) = pd(80,55) + rrt(721) * density(47) 
    rrt(730) = rrt(730) * 1.117D+05 / ysum
    pd(80,47) = pd(80,47) + rrt(730) * density(73) 
    pd(80,73) = pd(80,73) + rrt(730) * density(47) 
    rrt(731) = rrt(731) * 1.002D+05 / ysum
    pd(80,46) = pd(80,46) + rrt(731) * density(73) 
    pd(80,73) = pd(80,73) + rrt(731) * density(46) 
    rrt(732) = rrt(732) * 9.245D+04 / ysum
    pd(80,48) = pd(80,48) + rrt(732) * density(73) 
    pd(80,73) = pd(80,73) + rrt(732) * density(48) 
    rrt(734) = rrt(734) * 9.048D+04 / ysum
    pd(80,60) = pd(80,60) + rrt(734) * density(73) 
    pd(80,73) = pd(80,73) + rrt(734) * density(60) 
    rrt(735) = rrt(735) * 7.117D+04 / ysum
    pd(80,61) = pd(80,61) + rrt(735) * density(73) 
    pd(80,73) = pd(80,73) + rrt(735) * density(61) 
    rrt(736) = rrt(736) * 1.143D+05 / ysum
    pd(80,59) = pd(80,59) + rrt(736) * density(73) 
    pd(80,73) = pd(80,73) + rrt(736) * density(59) 
    rrt(740) = rrt(740) * 1.721D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(740) * density(67) 
    pd(80,67) = pd(80,67) + rrt(740) * density(27) 
    rrt(743) = rrt(743) * 1.702D+05 / ysum
    pd(80,01) = pd(80,01) + rrt(743) * density(65) * density(67) 
    pd(80,65) = pd(80,65) + rrt(743) * density(01) * density(67) 
    pd(80,67) = pd(80,67) + rrt(743) * density(01) * density(65) 
    rrt(744) = rrt(744) * 1.702D+05 / ysum
    pd(80,62) = pd(80,62) + rrt(744) * density(65) * density(67) 
    pd(80,65) = pd(80,65) + rrt(744) * density(62) * density(67) 
    pd(80,67) = pd(80,67) + rrt(744) * density(62) * density(65) 
    rrt(745) = rrt(745) * 1.702D+05 / ysum
    pd(80,65) = pd(80,65) + rrt(745) * density(67) * density(70) 
    pd(80,67) = pd(80,67) + rrt(745) * density(65) * density(70) 
    pd(80,70) = pd(80,70) + rrt(745) * density(65) * density(67) 
    rrt(749) = rrt(749) * 1.721D+05 / ysum
    pd(80,01) = pd(80,01) + rrt(749) * density(27) * density(67) 
    pd(80,27) = pd(80,27) + rrt(749) * density(01) * density(67) 
    pd(80,67) = pd(80,67) + rrt(749) * density(01) * density(27) 
    rrt(750) = rrt(750) * 1.721D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(750) * density(62) * density(67) 
    pd(80,62) = pd(80,62) + rrt(750) * density(27) * density(67) 
    pd(80,67) = pd(80,67) + rrt(750) * density(27) * density(62) 
    rrt(751) = rrt(751) * 1.721D+05 / ysum
    pd(80,27) = pd(80,27) + rrt(751) * density(67) * density(70) 
    pd(80,67) = pd(80,67) + rrt(751) * density(27) * density(70) 
    pd(80,70) = pd(80,70) + rrt(751) * density(27) * density(67) 
    rrt(755) = rrt(755) * 8.175D+04 / ysum
    pd(80,01) = pd(80,01) + rrt(755) * density(67) * density(76) 
    pd(80,67) = pd(80,67) + rrt(755) * density(01) * density(76) 
    pd(80,76) = pd(80,76) + rrt(755) * density(01) * density(67) 
    rrt(756) = rrt(756) * 8.175D+04 / ysum
    pd(80,62) = pd(80,62) + rrt(756) * density(67) * density(76) 
    pd(80,67) = pd(80,67) + rrt(756) * density(62) * density(76) 
    pd(80,76) = pd(80,76) + rrt(756) * density(62) * density(67) 
    rrt(757) = rrt(757) * 8.175D+04 / ysum
    pd(80,67) = pd(80,67) + rrt(757) * density(70) * density(76) 
    pd(80,70) = pd(80,70) + rrt(757) * density(67) * density(76) 
    pd(80,76) = pd(80,76) + rrt(757) * density(67) * density(70) 
    pd(80,1) = pd(80,1) + eV_to_K * ZDPlasKin_cfg(11)
    pd(80,:) = pd(80,:) * ZDPlasKin_cfg(13)
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
  DOUBLE PRECISION :: DTION, TIONN, TIONN2, TIONN3, TIONN4, TEFFN, TEFFN2, TEFFN3, TEFFN4 ! K
  DOUBLE PRECISION, PARAMETER :: ENERGY_VIBN2 = 0.290D0*11605.0D0, ENERGY_VIBO2 = 0.190D0*11605.0D0 ! K
  DOUBLE PRECISION :: QVIBN2, KVT10_N2N2, KVT01_N2N2, KVT10_N2N, KVT01_N2N, KVT10_N2O, KVT01_N2O ! CM3.S-1
  DOUBLE PRECISION :: QVIBO2, KVT10_O2O2, KVT01_O2O2, KVT10_O2O, KVT01_O2O ! CM3.S-1
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  call ZDPlasKin_get_density_total(ALL_NEUTRAL=ANY_NEUTRAL)
  DTION = 2.0D0 / ( 3.0D0 * 1.3807D-16 ) * 1.6605D-24 * ( 1.0D-17 * EN )**2
  TIONN = TGAS + DTION * 14.0D0 * 8.0D19**2
  TIONN2 = TGAS + DTION * 28.0D0 * 4.1D19**2
  TIONN3 = TGAS + DTION * 42.0D0 * 6.1D19**2
  TIONN4 = TGAS + DTION * 56.0D0 * 7.0D19**2
  TEFFN = ( TIONN + 0.5D0 * TGAS ) / ( 1.0D0 + 0.5D0 )
  TEFFN2 = ( TIONN2 + 1.0D0 * TGAS ) / ( 1.0D0 + 1.0D0 )
  TEFFN3 = ( TIONN3 + 1.5D0 * TGAS ) / ( 1.0D0 + 1.5D0 )
  TEFFN4 = ( TIONN4 + 2.0D0 * TGAS ) / ( 1.0D0 + 2.0D0 )
  QVIBN2 = EXP( - ENERGY_VIBN2 / TGAS )
  KVT10_N2N2 = 7.80D-12 * TGAS * EXP( - 218.0 / TGAS**(1.0/3.0) + 690.0 / TGAS ) / ( 1.0 - QVIBN2 )
  KVT10_N2N = 4.00D-16 * ( TGAS / 300.0D0 )**0.5
  KVT10_N2O = 1.20D-13 * EXP( - 27.6 / TGAS**(1.0/3.0) )
  KVT01_N2N2 = KVT10_N2N2 * QVIBN2
  KVT01_N2N = KVT10_N2N * QVIBN2
  KVT01_N2O = KVT10_N2O * QVIBN2
  QVIBO2 = EXP( - ENERGY_VIBO2 / TGAS )
  KVT10_O2O2 = 1.35D-12 * TGAS * EXP( - 137.9 / TGAS**(1.0/3.0) ) / ( 1.0 - QVIBO2 )
  KVT10_O2O = 4.50D-15 * TGAS
  KVT01_O2O2 = KVT10_O2O2 * QVIBO2
  KVT01_O2O = KVT10_O2O * QVIBO2
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
  rrt(016) = bolsig_rates(bolsig_pointer(16))
  rrt(017) = bolsig_rates(bolsig_pointer(17))
  rrt(018) = bolsig_rates(bolsig_pointer(18))
  rrt(019) = bolsig_rates(bolsig_pointer(19))
  rrt(020) = bolsig_rates(bolsig_pointer(20))
  rrt(021) = bolsig_rates(bolsig_pointer(21))
  rrt(022) = bolsig_rates(bolsig_pointer(22))
  rrt(023) = bolsig_rates(bolsig_pointer(23))
  rrt(024) = bolsig_rates(bolsig_pointer(24))
  rrt(025) = bolsig_rates(bolsig_pointer(25))
  rrt(026) = bolsig_rates(bolsig_pointer(26))
  rrt(027) = bolsig_rates(bolsig_pointer(27))
  rrt(028) = bolsig_rates(bolsig_pointer(28))
  rrt(029) = bolsig_rates(bolsig_pointer(29))
  rrt(030) = bolsig_rates(bolsig_pointer(30))
  rrt(031) = bolsig_rates(bolsig_pointer(31))
  rrt(032) = bolsig_rates(bolsig_pointer(32))
  rrt(033) = bolsig_rates(bolsig_pointer(33))
  rrt(034) = bolsig_rates(bolsig_pointer(34))
  rrt(035) = bolsig_rates(bolsig_pointer(35))
  rrt(036) = bolsig_rates(bolsig_pointer(36))
  rrt(037) = bolsig_rates(bolsig_pointer(37))
  rrt(038) = bolsig_rates(bolsig_pointer(38))
  rrt(039) = bolsig_rates(bolsig_pointer(39))
  rrt(040) = bolsig_rates(bolsig_pointer(40))
  rrt(041) = bolsig_rates(bolsig_pointer(41))
  rrt(042) = bolsig_rates(bolsig_pointer(42))
  rrt(043) = bolsig_rates(bolsig_pointer(43))
  rrt(044) = bolsig_rates(bolsig_pointer(44))
  rrt(045) = bolsig_rates(bolsig_pointer(45))
  rrt(046) = bolsig_rates(bolsig_pointer(46))
  rrt(047) = bolsig_rates(bolsig_pointer(47))
  rrt(048) = bolsig_rates(bolsig_pointer(48))
  rrt(049) = bolsig_rates(bolsig_pointer(49))
  rrt(050) = KVT10_N2N2*1.0D0
  rrt(051) = KVT10_N2N2*2.0D0
  rrt(052) = KVT10_N2N2*3.0D0
  rrt(053) = KVT10_N2N2*4.0D0
  rrt(054) = KVT10_N2N2*5.0D0
  rrt(055) = KVT10_N2N2*6.0D0
  rrt(056) = KVT10_N2N2*7.0D0
  rrt(057) = KVT10_N2N2*8.0D0
  rrt(058) = KVT01_N2N2*1.0D0
  rrt(059) = KVT01_N2N2*2.0D0
  rrt(060) = KVT01_N2N2*3.0D0
  rrt(061) = KVT01_N2N2*4.0D0
  rrt(062) = KVT01_N2N2*5.0D0
  rrt(063) = KVT01_N2N2*6.0D0
  rrt(064) = KVT01_N2N2*7.0D0
  rrt(065) = KVT01_N2N2*8.0D0
  rrt(066) = KVT10_N2N*1.0D0
  rrt(067) = KVT10_N2N*2.0D0
  rrt(068) = KVT10_N2N*3.0D0
  rrt(069) = KVT10_N2N*4.0D0
  rrt(070) = KVT10_N2N*5.0D0
  rrt(071) = KVT10_N2N*6.0D0
  rrt(072) = KVT10_N2N*7.0D0
  rrt(073) = KVT10_N2N*8.0D0
  rrt(074) = KVT01_N2N*1.0D0
  rrt(075) = KVT01_N2N*2.0D0
  rrt(076) = KVT01_N2N*3.0D0
  rrt(077) = KVT01_N2N*4.0D0
  rrt(078) = KVT01_N2N*5.0D0
  rrt(079) = KVT01_N2N*6.0D0
  rrt(080) = KVT01_N2N*7.0D0
  rrt(081) = KVT01_N2N*8.0D0
  rrt(082) = KVT10_N2O*1.0D0
  rrt(083) = KVT10_N2O*2.0D0
  rrt(084) = KVT10_N2O*3.0D0
  rrt(085) = KVT10_N2O*4.0D0
  rrt(086) = KVT10_N2O*5.0D0
  rrt(087) = KVT10_N2O*6.0D0
  rrt(088) = KVT10_N2O*7.0D0
  rrt(089) = KVT10_N2O*8.0D0
  rrt(090) = KVT01_N2O*1.0D0
  rrt(091) = KVT01_N2O*2.0D0
  rrt(092) = KVT01_N2O*3.0D0
  rrt(093) = KVT01_N2O*4.0D0
  rrt(094) = KVT01_N2O*5.0D0
  rrt(095) = KVT01_N2O*6.0D0
  rrt(096) = KVT01_N2O*7.0D0
  rrt(097) = KVT01_N2O*8.0D0
  rrt(098) = KVT10_O2O2*1.0D0
  rrt(099) = KVT10_O2O2*2.0D0
  rrt(100) = KVT10_O2O2*3.0D0
  rrt(101) = KVT10_O2O2*4.0D0
  rrt(102) = KVT01_O2O2*1.0D0
  rrt(103) = KVT01_O2O2*2.0D0
  rrt(104) = KVT01_O2O2*3.0D0
  rrt(105) = KVT01_O2O2*4.0D0
  rrt(106) = KVT10_O2O*1.0D0
  rrt(107) = KVT10_O2O*2.0D0
  rrt(108) = KVT10_O2O*3.0D0
  rrt(109) = KVT10_O2O*4.0D0
  rrt(110) = KVT01_O2O*1.0D0
  rrt(111) = KVT01_O2O*2.0D0
  rrt(112) = KVT01_O2O*3.0D0
  rrt(113) = KVT01_O2O*4.0D0
  rrt(114) = bolsig_rates(bolsig_pointer(50))
  rrt(115) = bolsig_rates(bolsig_pointer(51))
  rrt(116) = bolsig_rates(bolsig_pointer(52))
  rrt(117) = bolsig_rates(bolsig_pointer(53))
  rrt(118) = bolsig_rates(bolsig_pointer(54))
  rrt(119) = bolsig_rates(bolsig_pointer(55))
  rrt(120) = bolsig_rates(bolsig_pointer(56))
  rrt(121) = bolsig_rates(bolsig_pointer(57))
  rrt(122) = bolsig_rates(bolsig_pointer(58))
  rrt(123) = bolsig_rates(bolsig_pointer(59))
  rrt(124) = bolsig_rates(bolsig_pointer(60))
  rrt(125) = bolsig_rates(bolsig_pointer(61))
  rrt(126) = bolsig_rates(bolsig_pointer(62))
  rrt(127) = bolsig_rates(bolsig_pointer(63))
  rrt(128) = bolsig_rates(bolsig_pointer(64))
  rrt(129) = bolsig_rates(bolsig_pointer(65))
  rrt(130) = bolsig_rates(bolsig_pointer(66))
  rrt(131) = bolsig_rates(bolsig_pointer(67))
  rrt(132) = bolsig_rates(bolsig_pointer(68))
  rrt(133) = bolsig_rates(bolsig_pointer(69))
  rrt(134) = bolsig_rates(bolsig_pointer(70))
  rrt(135) = bolsig_rates(bolsig_pointer(71))
  rrt(136) = bolsig_rates(bolsig_pointer(72))
  rrt(137) = bolsig_rates(bolsig_pointer(73))
  rrt(138) = bolsig_rates(bolsig_pointer(74))
  rrt(139) = bolsig_rates(bolsig_pointer(75))
  rrt(140) = bolsig_rates(bolsig_pointer(76))
  rrt(141) = bolsig_rates(bolsig_pointer(77))
  rrt(142) = bolsig_rates(bolsig_pointer(78))
  rrt(143) = bolsig_rates(bolsig_pointer(79))
  rrt(144) = bolsig_rates(bolsig_pointer(80))
  rrt(145) = bolsig_rates(bolsig_pointer(81))
  rrt(146) = bolsig_rates(bolsig_pointer(82))
  rrt(147) = bolsig_rates(bolsig_pointer(83))
  rrt(148) = bolsig_rates(bolsig_pointer(84))
  rrt(149) = bolsig_rates(bolsig_pointer(85))
  rrt(150) = bolsig_rates(bolsig_pointer(86))
  rrt(151) = bolsig_rates(bolsig_pointer(87))
  rrt(152) = bolsig_rates(bolsig_pointer(88))
  rrt(153) = bolsig_rates(bolsig_pointer(89))
  rrt(154) = bolsig_rates(bolsig_pointer(90))
  rrt(155) = bolsig_rates(bolsig_pointer(91))
  rrt(156) = bolsig_rates(bolsig_pointer(92))
  rrt(157) = bolsig_rates(bolsig_pointer(93))
  rrt(158) = bolsig_rates(bolsig_pointer(94))
  rrt(159) = bolsig_rates(bolsig_pointer(95))
  rrt(160) = bolsig_rates(bolsig_pointer(96))
  rrt(161) = bolsig_rates(bolsig_pointer(97))
  rrt(162) = bolsig_rates(bolsig_pointer(98))
  rrt(163) = bolsig_rates(bolsig_pointer(99))
  rrt(164) = bolsig_rates(bolsig_pointer(100))
  rrt(165) = bolsig_rates(bolsig_pointer(101))
  rrt(166) = 1.8D-7*(300.0D0/TE)**0.39*0.50D0
  rrt(167) = 1.8D-7*(300.0D0/TE)**0.39*0.45D0
  rrt(168) = 1.8D-7*(300.0D0/TE)**0.39*0.05D0
  rrt(169) = 2.7D-7*(300.0D0/TE)**0.7*0.55D0
  rrt(170) = 2.7D-7*(300.0D0/TE)**0.7*0.40D0
  rrt(171) = 2.7D-7*(300.0D0/TE)**0.7*0.05D0
  rrt(172) = 4.2D-7*(300.0D0/TE)**0.85*0.20D0
  rrt(173) = 4.2D-7*(300.0D0/TE)**0.85*0.80D0
  rrt(174) = 2.0D-7*(300.0D0/TE)**0.5
  rrt(175) = 2.3D-6*(300.0D0/TE)**0.53
  rrt(176) = rrt(174)
  rrt(177) = rrt(174)
  rrt(178) = 1.4D-6*(300.0D0/TE)**0.5
  rrt(179) = 1.3D-6*(300.0D0/TE)**0.5
  rrt(180) = 7.0D-20*(300.0D0/TE)**4.5
  rrt(181) = rrt(180)
  rrt(182) = 6.0D-27*(300.0D0/TE)**1.5*ANY_NEUTRAL
  rrt(183) = rrt(182)
  rrt(184) = 3.00D-7
  rrt(185) = bolsig_rates(bolsig_pointer(102))
  rrt(186) = bolsig_rates(bolsig_pointer(103))
  rrt(187) = bolsig_rates(bolsig_pointer(104))
  rrt(188) = bolsig_rates(bolsig_pointer(105))
  rrt(189) = 1.0D-11
  rrt(190) = 1.0D-31
  rrt(191) = 1.0D-31
  rrt(192) = 1.0D-31*ANY_NEUTRAL
  rrt(193) = 8.0D-31*ANY_NEUTRAL
  rrt(194) = 6.0D-33*ANY_NEUTRAL
  rrt(195) = 1.1D-31*(300.0D0/TE)**2*EXP(-70.0D0/TGAS)*EXP(1500.0D0*(TE-TGAS)/(TE*TGAS))
  rrt(196) = bolsig_rates(bolsig_pointer(106))
  rrt(197) = 1.4D-10
  rrt(198) = 2.6D-10
  rrt(199) = 2.6D-10
  rrt(200) = 5.0D-13
  rrt(201) = 5.0D-15
  rrt(202) = 3.0D-10
  rrt(203) = 6.9D-10
  rrt(204) = 2.2D-9
  rrt(205) = 1.9D-9
  rrt(206) = 3.0D-10
  rrt(207) = 1.5D-10
  rrt(208) = 5.0D-10
  rrt(209) = 2.7D-10*(TEFFN2/300.0D0)**0.5*EXP(-5590.0D0/TEFFN2)
  rrt(210) = 2.0D-10
  rrt(211) = 3.6D-10
  rrt(212) = 1.9D-12*(TEFFN2/300.0D0)**0.5*EXP(-4990.0D0/TEFFN2)
  rrt(213) = 2.1D-9
  rrt(214) = 2.5D-9
  rrt(215) = 3.0D-10
  rrt(216) = 5.0D-10
  rrt(217) = 5.0D-10
  rrt(218) = 5.0D-10
  rrt(219) = 5.0D-10
  rrt(220) = 5.0D-10
  rrt(221) = 1.5D-10
  rrt(222) = 1.5D-10
  rrt(223) = 1.5D-10
  rrt(224) = 1.5D-10
  rrt(225) = 2.1D-9
  rrt(226) = 2.1D-9
  rrt(227) = 2.1D-9
  rrt(228) = 2.1D-9
  rrt(229) = 2.1D-9
  rrt(230) = 2.5D-9
  rrt(231) = 2.5D-9
  rrt(232) = 2.5D-9
  rrt(233) = 2.5D-9
  rrt(234) = 2.5D-9
  rrt(235) = 0.50D0
  rrt(236) = 1.34D5
  rrt(237) = 1.0D2
  rrt(238) = 2.45D7
  rrt(239) = 2.6D-4
  rrt(240) = 1.5D-3
  rrt(241) = 8.5D-2
  rrt(242) = 11.0D0
  rrt(243) = 1.7D-12
  rrt(244) = 7.5D-13
  rrt(245) = 3.0D-10
  rrt(246) = 2.8D-11
  rrt(247) = 1.15D-10
  rrt(248) = 1.30D-10
  rrt(249) = 5.00D-11
  rrt(250) = 7.0D-12
  rrt(251) = 2.1D-11
  rrt(252) = 2.0D-12
  rrt(253) = 4.0D-11*(300.0D0/TGAS)**0.667
  rrt(254) = 3.0D-16
  rrt(255) = 6.9D-11
  rrt(256) = 1.0D-11
  rrt(257) = 1.0D-12
  rrt(258) = 3.0D-10
  rrt(259) = 1.5D-10
  rrt(260) = 3.0D-11
  rrt(261) = 2.0D-12
  rrt(262) = 2.4D-10
  rrt(263) = 1.0D-11
  rrt(264) = 1.9D-13
  rrt(265) = 3.6D-10
  rrt(266) = 4.0D-12
  rrt(267) = 1.0D-11
  rrt(268) = 5.0D-11
  rrt(269) = 2.0D-10
  rrt(270) = 1.6D-10
  rrt(271) = 2.5D-11
  rrt(272) = 1.5D-11
  rrt(273) = 2.6D-11
  rrt(274) = 4.4D-36
  rrt(275) = 4.4D-36
  rrt(276) = 4.4D-36
  rrt(277) = 2.6D-35
  rrt(278) = 2.6D-35
  rrt(279) = 6.2D-36
  rrt(280) = 6.2D-36
  rrt(281) = 6.2D-36
  rrt(282) = 3.6D-35
  rrt(283) = 3.6D-35
  rrt(284) = 4.0D-13
  rrt(285) = 5.2D-12
  rrt(286) = 1.8D-10
  rrt(287) = 3.5D-12
  rrt(288) = 1.0D-13*EXP(-510.0D0/TGAS)
  rrt(289) = 1.8D-12
  rrt(290) = 1.0D-12
  rrt(291) = 6.0D-13
  rrt(292) = 6.0D-14
  rrt(293) = 1.0D-13
  rrt(294) = 2.6D-12
  rrt(295) = 3.0D-11
  rrt(296) = 2.3D-12
  rrt(297) = 1.1D-10
  rrt(298) = 2.5D-14
  rrt(299) = 7.0D-16
  rrt(300) = 2.0D-14*EXP(-600.0D0/TGAS)
  rrt(301) = 3.8D-18*EXP(-205.0D0/TGAS)
  rrt(302) = 3.0D-21
  rrt(303) = 2.5D-11
  rrt(304) = 5.2D-11*EXP(-2840.0D0/TGAS)
  rrt(305) = 7.0D-28*TGAS**3.8*EXP(700.0D0/TGAS)
  rrt(306) = 1.0D-11*EXP(-2300.0D0/TGAS)
  rrt(307) = 8.1D-14
  rrt(308) = 3.4D-11*(300.0D0/TGAS)**0.1*EXP(-4200.0D0/TGAS)
  rrt(309) = 4.3D-22*TGAS**2.4*EXP(-281.0D0/TGAS)
  rrt(310) = 1.7D-15*(TGAS/300.0D0)
  rrt(311) = 6.0D-14
  rrt(312) = 2.2D-11
  rrt(313) = 9.0D-12
  rrt(314) = 3.0D-13
  rrt(315) = 9.0D-15
  rrt(316) = 8.0D-12
  rrt(317) = 6.4D-12*EXP(67.0D0/TGAS)
  rrt(318) = 1.0D-12
  rrt(319) = 2.6D-11*EXP(67.0D0/TGAS)
  rrt(320) = 2.3D-11
  rrt(321) = 1.2D-10
  rrt(322) = 1.2D-10
  rrt(323) = 1.7D-10
  rrt(324) = 7.2D-11
  rrt(325) = 4.4D-11
  rrt(326) = 5.0D-11*EXP(-300.0D0/TGAS)
  rrt(327) = 1.0D-12
  rrt(328) = 1.3D-12*EXP(-850.0D0/TGAS)
  rrt(329) = 3.0D-12*EXP(-850.0D0/TGAS)
  rrt(330) = 1.0D-17
  rrt(331) = 1.1D-10
  rrt(332) = 2.9D-11
  rrt(333) = 3.2D-11
  rrt(334) = 2.9D-10
  rrt(335) = 5.1D-10
  rrt(336) = 2.9D-10
  rrt(337) = 2.9D-10
  rrt(338) = 6.3D-12
  rrt(339) = 3.1D-12
  rrt(340) = 1.8D-11*(TGAS/300.0)**0.5
  rrt(341) = 3.2D-12*(TGAS/300.0)*EXP(-3150.0D0/TGAS)
  rrt(342) = 9.1D-13
  rrt(343) = 3.0D-12
  rrt(344) = 7.0D-13
  rrt(345) = 2.3D-12
  rrt(346) = 3.0D-10*EXP(-38370.0D0/TGAS)
  rrt(347) = 7.5D-12*(TGAS/300.0)*EXP(-19500.0D0/TGAS)
  rrt(348) = 4.2D-18
  rrt(349) = 8.3D-12*EXP(-14000.0D0/TGAS)
  rrt(350) = 1.5D-10*EXP(-14090.0D0/TGAS)
  rrt(351) = 9.1D-12*(TGAS/300.0D0)**0.18
  rrt(352) = 1.0D-11
  rrt(353) = 2.5D-10*EXP(-50390.0D0/TGAS)
  rrt(354) = 3.3D-16*(300.0D0/TGAS)**0.5*EXP(-39200.0D0/TGAS)
  rrt(355) = 2.2D-12*EXP(-32100.0D0/TGAS)
  rrt(356) = 5.1D-13*EXP(-33660.0D0/TGAS)
  rrt(357) = 2.8D-12*EXP(-23400.0D0/TGAS)
  rrt(358) = 2.5D-13*EXP(-765.0D0/TGAS)
  rrt(359) = 4.6D-10*EXP(-25170.0D0/TGAS)
  rrt(360) = 1.7D-11
  rrt(361) = 2.0D-11*EXP(-49800.0D0/TGAS)
  rrt(362) = 2.8D-12*EXP(-25400.0D0/TGAS)
  rrt(363) = 3.3D-12*EXP(-13500.0D0/TGAS)
  rrt(364) = 4.5D-10*EXP(-18500.0D0/TGAS)
  rrt(365) = 1.2D-13*EXP(-2450.0D0/TGAS)
  rrt(366) = 2.3D-13*EXP(-1600.0D0/TGAS)
  rrt(367) = 1.5D-12*EXP(-15020.0D0/TGAS)
  rrt(368) = 4.3D-12*EXP(-3850.0D0/TGAS)
  rrt(369) = 2.7D-11*EXP(-6.74D4/TGAS)
  rrt(370) = 1.6D-12*(TGAS/300.0D0)**0.5*(0.19D0+8.6D0*TGAS)*EXP(-32000.0D0/TGAS)
  rrt(371) = 5.4D-11*EXP(-165/TGAS)
  rrt(372) = 6.6D-11*EXP(-1840/TGAS)
  rrt(373) = 8.4D-14*(TGAS/300)**4.1*EXP(-4760/TGAS)
  rrt(374) = 5.4D-11*EXP(-6492/TGAS)
  rrt(375) = 5.0D-11
  rrt(376) = 1.2D-10
  rrt(377) = 1.2D-10
  rrt(378) = 1.7D-12*(TGAS/300)**1.5
  rrt(379) = 8.5D-11
  rrt(380) = 5.0D-14*(TGAS/300)
  rrt(381) = 1.66D-12
  rrt(382) = 1.0D-13
  rrt(383) = 1.4D-12
  rrt(384) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*1.0D0
  rrt(385) = rrt(384)
  rrt(386) = rrt(384)
  rrt(387) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*6.6D0
  rrt(388) = rrt(387)
  rrt(389) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*1.0D0
  rrt(390) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*5.9D0
  rrt(391) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*21.D0
  rrt(392) = rrt(389)
  rrt(393) = rrt(389)
  rrt(394) = 8.7D-9*EXP(-75994.0D0/TGAS)*1.0D0
  rrt(395) = rrt(394)
  rrt(396) = 8.7D-9*EXP(-75994.0D0/TGAS)*20.D0
  rrt(397) = rrt(396)
  rrt(398) = rrt(396)
  rrt(399) = 6.6D-10*EXP(-11600.0D0/TGAS)*1.0D0
  rrt(400) = 6.6D-10*EXP(-11600.0D0/TGAS)*0.38D0
  rrt(401) = 6.6D-10*EXP(-11600.0D0/TGAS)*6.3D0*EXP(170.0D0/TGAS)
  rrt(402) = rrt(401)
  rrt(403) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*1.0D0
  rrt(404) = rrt(403)
  rrt(405) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*2.0D0
  rrt(406) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*4.0D0
  rrt(407) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*1.0D0
  rrt(408) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*0.78D0
  rrt(409) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*7.8D0
  rrt(410) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*5.9D0
  rrt(411) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(412) = rrt(411)
  rrt(413) = rrt(411)
  rrt(414) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*10.D0
  rrt(415) = rrt(414)
  rrt(416) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(417) = rrt(416)
  rrt(418) = rrt(416)
  rrt(419) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*12.D0
  rrt(420) = rrt(419)
  rrt(421) = 2.1D-11*(300.0D0/TGAS)**4.4*EXP(-11080.0D0/TGAS)*ANY_NEUTRAL
  rrt(422) = MAX(8.3D-34*EXP(500.0D0/TGAS),1.91D-33)
  rrt(423) = 1.8D-33*EXP(435.0D0/TGAS)*1.0D0
  rrt(424) = rrt(423)
  rrt(425) = 1.8D-33*EXP(435.0D0/TGAS)*3.0D0
  rrt(426) = rrt(425)
  rrt(427) = MAX(2.8D-34*EXP(720.0D0/TGAS),1.0D-33*(300.0D0/TGAS)**0.41)
  rrt(428) = 4.0D-33*(300.0D0/TGAS)**0.41*1.0D0
  rrt(429) = 4.0D-33*(300.0D0/TGAS)**0.41*0.8D0
  rrt(430) = 4.0D-33*(300.0D0/TGAS)**0.41*3.6D0
  rrt(431) = 4.0D-33*(300.0D0/TGAS)**0.41*0.17D0
  rrt(432) = 1.0D-32*(300.0D0/TGAS)**0.5
  rrt(433) = rrt(432)
  rrt(434) = 1.8D-31*(300.0D0/TGAS)
  rrt(435) = rrt(434)
  rrt(436) = rrt(434)
  rrt(437) = MAX(5.8D-34*(300.0D0/TGAS)**2.8,5.4D-34*(300.0D0/TGAS)**1.9)
  rrt(438) = 7.6D-34*(300.0D0/TGAS)**1.9
  rrt(439) = rrt(438)
  rrt(440) = MIN(3.9D-33*(300.0D0/TGAS)**1.9,1.1D-34*EXP(1060.0D0/TGAS))
  rrt(441) = rrt(440)
  rrt(442) = 3.9D-35*EXP(-10400.0D0/TGAS)*ANY_NEUTRAL
  rrt(443) = 1.2D-31*(300.0D0/TGAS)**1.8*1.0D0
  rrt(444) = 1.2D-31*(300.0D0/TGAS)**1.8*0.78D0
  rrt(445) = rrt(444)
  rrt(446) = 8.9D-32*(300.0D0/TGAS)**2*1.0D0
  rrt(447) = rrt(446)
  rrt(448) = 8.9D-32*(300.0D0/TGAS)**2*13.D0
  rrt(449) = rrt(448)
  rrt(450) = 8.9D-32*(300.0D0/TGAS)**2*2.4D0
  rrt(451) = 3.7D-30*(300.0D0/TGAS)**4.1*ANY_NEUTRAL
  rrt(452) = 2.4D-36*EXP(500/TGAS)
  rrt(453) = rrt(452)
  rrt(454) = rrt(452)
  rrt(455) = 2.3D-35*(300/TGAS)**0.6
  rrt(456) = 2.2D-35*(300/TGAS)
  rrt(457) = 2.6D-36
  rrt(458) = 2.6D-36
  rrt(459) = 2.6D-36
  rrt(460) = 2.6D-37
  rrt(461) = 2.6D-37
  rrt(462) = 2.6D-37
  rrt(463) = 2.6D-35
  rrt(464) = 2.6D-35
  rrt(465) = 2.6D-35
  rrt(466) = 1.4D-32
  rrt(467) = 1.4D-32
  rrt(468) = 1.4D-32
  rrt(469) = 6.5D-38*(TGAS/300)*EXP(1700/TGAS)
  rrt(470) = rrt(469)
  rrt(471) = rrt(469)
  rrt(472) = 1.0D-12
  rrt(473) = 2.8D-10
  rrt(474) = 2.5D-10
  rrt(475) = 2.8D-11
  rrt(476) = 5.0D-10
  rrt(477) = 8.0D-10
  rrt(478) = 3.0D-12
  rrt(479) = 1.0D-12
  rrt(480) = 5.5D-10
  rrt(481) = (1.5D0-2.0D-3*TEFFN+9.6D-7*TEFFN**2)*1.0D-12
  rrt(482) = 2.0D-11*(300.0D0/TEFFN)**0.5
  rrt(483) = 1.0D-10
  rrt(484) = 2.4D-11
  rrt(485) = 3.0D-12
  rrt(486) = 1.3D-10
  rrt(487) = 2.3D-10
  rrt(488) = 2.2D-10
  rrt(489) = 2.0D-11
  rrt(490) = 1.6D-9
  rrt(491) = 6.0D-11*(300.0D0/TEFFN2)**0.5
  rrt(492) = 1.3D-10*(300.0D0/TEFFN2)**0.5
  rrt(493) = 1.0D-10
  rrt(494) = 7.2D-13*(TEFFN2/300.0D0)
  rrt(495) = 3.3D-10
  rrt(496) = 5.0D-10
  rrt(497) = 4.0D-10
  rrt(498) = 1.0D-17
  rrt(499) = 1.2D-10
  rrt(500) = 6.3D-10
  rrt(501) = 1.0D-11
  rrt(502) = 6.6D-10
  rrt(503) = 2.3D-11
  rrt(504) = 4.4D-11
  rrt(505) = 6.6D-11
  rrt(506) = 7.0D-11
  rrt(507) = 7.0D-11
  rrt(508) = 2.9D-10
  rrt(509) = 2.9D-10
  rrt(510) = MIN(2.1D-16*EXP(TEFFN4/121.0D0),1.0D-10)
  rrt(511) = 2.5D-10
  rrt(512) = 2.5D-10
  rrt(513) = 1.0D-11
  rrt(514) = 4.0D-10
  rrt(515) = 4.6D-12*(TEFFN4/300.0D0)**2.5*EXP(-2650.0D0/TEFFN4)
  rrt(516) = 3.3D-6*(300.0D0/TEFFN4)**4*EXP(-5030.0D0/TEFFN4)
  rrt(517) = 1.0D-10
  rrt(518) = 1.0D-10
  rrt(519) = 3.0D-10
  rrt(520) = 1.0D-10
  rrt(521) = 1.1D-6*(300.0D0/TEFFN4)**5.3*EXP(-2360.0D0/TEFFN4)
  rrt(522) = 1.0D-9
  rrt(523) = 1.7D-29*(300.0D0/TEFFN)**2.1
  rrt(524) = 1.0D-29*ANY_NEUTRAL
  rrt(525) = rrt(524)
  rrt(526) = 6.0D-29*(300.0D0/TEFFN)**2*ANY_NEUTRAL
  rrt(527) = rrt(524)
  rrt(528) = rrt(524)
  rrt(529) = 5.2D-29*(300.0D0/TEFFN2)**2.2
  rrt(530) = 9.0D-30*EXP(400.0D0/TEFFN2)
  rrt(531) = 2.4D-30*(300.0D0/TEFFN2)**3.2
  rrt(532) = 9.0D-31*(300.0D0/TEFFN2)**2
  rrt(533) = 5.0D-10
  rrt(534) = 4.7D-10
  rrt(535) = 2.12D-10
  rrt(536) = 5.20D-9
  rrt(537) = 6.40D-10
  rrt(538) = 2.00D-9
  rrt(539) = 5.70D-9
  rrt(540) = 2.00D-9
  rrt(541) = 1.85D-10
  rrt(542) = 1.05D-9
  rrt(543) = 1.80D-9
  rrt(544) = 6.00D-10
  rrt(545) = 6.50D-10
  rrt(546) = 1.95D-10
  rrt(547) = 1.15D-9
  rrt(548) = 1.15D-9
  rrt(549) = 2.3D-9
  rrt(550) = 1.14D-9
  rrt(551) = 2.39D-9
  rrt(552) = 9.96D-10
  rrt(553) = 2.49D-9
  rrt(554) = 2.20D-9
  rrt(555) = 1.00D-9
  rrt(556) = 1.0D-10
  rrt(557) = 8.0D-10
  rrt(558) = 1.2D-9
  rrt(559) = 2.0D-10
  rrt(560) = 2.0D-12
  rrt(561) = 3.3D-10
  rrt(562) = 3.5D-10
  rrt(563) = 7.0D-10
  rrt(564) = 5.0D-10
  rrt(565) = 1.0D-11
  rrt(566) = 1.0D-11
  rrt(567) = 2.6D-12
  rrt(568) = 7.0D-11
  rrt(569) = 2.0D-11
  rrt(570) = 5.0D-10
  rrt(571) = 5.0D-10
  rrt(572) = 7.4D-10
  rrt(573) = 2.8D-14
  rrt(574) = 1.8D-11
  rrt(575) = 4.0D-12
  rrt(576) = 5.0D-10
  rrt(577) = 7.0D-10
  rrt(578) = 3.0D-15
  rrt(579) = 1.0D-10*EXP(-1044.0D0/TEFFN4)
  rrt(580) = rrt(579)
  rrt(581) = 4.0D-10
  rrt(582) = 3.0D-10
  rrt(583) = 1.0D-10
  rrt(584) = 1.0D-10
  rrt(585) = 2.5D-10
  rrt(586) = 1.1D-30*(300.0D0/TEFFN)*ANY_NEUTRAL
  rrt(587) = rrt(524)
  rrt(588) = 3.5D-31*(300.0D0/TEFFN2)*ANY_NEUTRAL
  rrt(589) = 2.0D-7*(300.0D0/TIONN)**0.5
  rrt(590) = rrt(589)
  rrt(591) = rrt(589)
  rrt(592) = rrt(589)
  rrt(593) = rrt(589)
  rrt(594) = rrt(589)
  rrt(595) = rrt(589)
  rrt(596) = rrt(589)
  rrt(597) = rrt(589)
  rrt(598) = rrt(589)
  rrt(599) = rrt(589)
  rrt(600) = rrt(589)
  rrt(601) = rrt(589)
  rrt(602) = rrt(589)
  rrt(603) = rrt(589)
  rrt(604) = rrt(589)
  rrt(605) = rrt(589)
  rrt(606) = rrt(589)
  rrt(607) = rrt(589)
  rrt(608) = rrt(589)
  rrt(609) = rrt(589)
  rrt(610) = rrt(589)
  rrt(611) = rrt(589)
  rrt(612) = rrt(589)
  rrt(613) = rrt(589)
  rrt(614) = rrt(589)
  rrt(615) = rrt(589)
  rrt(616) = rrt(589)
  rrt(617) = rrt(589)
  rrt(618) = rrt(589)
  rrt(619) = rrt(589)
  rrt(620) = rrt(589)
  rrt(621) = rrt(589)
  rrt(622) = rrt(589)
  rrt(623) = rrt(589)
  rrt(624) = rrt(589)
  rrt(625) = rrt(589)
  rrt(626) = rrt(589)
  rrt(627) = rrt(589)
  rrt(628) = rrt(589)
  rrt(629) = rrt(589)
  rrt(630) = rrt(589)
  rrt(631) = rrt(589)
  rrt(632) = rrt(589)
  rrt(633) = rrt(589)
  rrt(634) = rrt(589)
  rrt(635) = rrt(589)
  rrt(636) = rrt(589)
  rrt(637) = rrt(589)
  rrt(638) = 1.0D-7
  rrt(639) = 1.0D-7
  rrt(640) = 1.0D-7
  rrt(641) = 1.0D-7
  rrt(642) = 1.0D-7
  rrt(643) = 1.0D-7
  rrt(644) = 1.0D-7
  rrt(645) = 1.0D-7
  rrt(646) = 1.0D-7
  rrt(647) = 1.0D-7
  rrt(648) = 1.0D-7
  rrt(649) = 1.0D-7
  rrt(650) = 1.0D-7
  rrt(651) = 1.0D-7
  rrt(652) = 1.0D-7
  rrt(653) = 1.0D-7
  rrt(654) = 1.0D-7
  rrt(655) = 1.0D-7
  rrt(656) = 1.0D-7
  rrt(657) = 1.0D-7
  rrt(658) = 1.0D-7
  rrt(659) = 1.0D-7
  rrt(660) = 1.0D-7
  rrt(661) = 1.0D-7
  rrt(662) = 1.0D-7
  rrt(663) = 1.0D-7
  rrt(664) = 1.0D-7
  rrt(665) = 1.0D-7
  rrt(666) = 1.0D-7
  rrt(667) = 1.0D-7
  rrt(668) = 1.0D-7
  rrt(669) = 1.0D-7
  rrt(670) = 1.0D-7
  rrt(671) = 1.0D-7
  rrt(672) = 1.0D-7
  rrt(673) = 1.0D-7
  rrt(674) = 1.0D-7
  rrt(675) = 1.0D-7
  rrt(676) = 1.0D-7
  rrt(677) = 1.0D-7
  rrt(678) = 1.0D-7
  rrt(679) = 1.0D-7
  rrt(680) = 1.0D-7
  rrt(681) = 1.0D-7
  rrt(682) = 1.0D-7
  rrt(683) = 1.0D-7
  rrt(684) = 1.0D-7
  rrt(685) = 1.0D-7
  rrt(686) = 1.0D-7
  rrt(687) = 1.0D-7
  rrt(688) = 1.0D-7
  rrt(689) = 1.0D-7
  rrt(690) = 1.0D-7
  rrt(691) = 1.0D-7
  rrt(692) = 1.0D-7
  rrt(693) = 1.0D-7
  rrt(694) = 1.0D-7
  rrt(695) = 1.0D-7
  rrt(696) = 1.0D-7
  rrt(697) = 1.0D-7
  rrt(698) = 1.0D-7
  rrt(699) = 1.0D-7
  rrt(700) = 1.0D-7
  rrt(701) = 1.0D-7
  rrt(702) = 1.0D-7
  rrt(703) = 1.0D-7
  rrt(704) = 1.0D-7
  rrt(705) = 1.0D-7
  rrt(706) = 1.0D-7
  rrt(707) = 1.0D-7
  rrt(708) = 1.0D-7
  rrt(709) = 1.0D-7
  rrt(710) = 1.0D-7
  rrt(711) = 1.0D-7
  rrt(712) = 2.0D-25*(300.0D0/TIONN)**2.5*ANY_NEUTRAL
  rrt(713) = rrt(712)
  rrt(714) = rrt(712)
  rrt(715) = rrt(712)
  rrt(716) = rrt(712)
  rrt(717) = rrt(712)
  rrt(718) = rrt(712)
  rrt(719) = rrt(712)
  rrt(720) = rrt(712)
  rrt(721) = rrt(712)
  rrt(722) = rrt(712)
  rrt(723) = rrt(712)
  rrt(724) = rrt(712)
  rrt(725) = rrt(712)
  rrt(726) = rrt(712)
  rrt(727) = rrt(712)
  rrt(728) = rrt(712)
  rrt(729) = rrt(712)
  rrt(730) = 3.47D-6*TE**(-0.5)
  rrt(731) = rrt(730)
  rrt(732) = rrt(730)
  rrt(733) = rrt(730)
  rrt(734) = rrt(730)
  rrt(735) = rrt(730)
  rrt(736) = rrt(730)
  rrt(737) = 2.99D-6
  rrt(738) = 2D-7*(300/TGAS)
  rrt(739) = rrt(738)
  rrt(740) = rrt(738)
  rrt(741) = rrt(738)
  rrt(742) = rrt(738)
  rrt(743) = 5D-28*(300/TGAS)**2.5
  rrt(744) = rrt(743)
  rrt(745) = rrt(743)
  rrt(746) = rrt(743)
  rrt(747) = rrt(743)
  rrt(748) = rrt(743)
  rrt(749) = rrt(743)
  rrt(750) = rrt(743)
  rrt(751) = rrt(743)
  rrt(752) = rrt(743)
  rrt(753) = rrt(743)
  rrt(754) = rrt(743)
  rrt(755) = rrt(743)
  rrt(756) = rrt(743)
  rrt(757) = rrt(743)
  rrt(758) = 2.0D-25*(300.0D0/TIONN2)**2.5*ANY_NEUTRAL
  rrt(759) = rrt(758)
  rrt(760) = rrt(758)
  rrt(761) = rrt(758)
  rrt(762) = rrt(758)
  rrt(763) = rrt(758)
  rrt(764) = rrt(758)
  rrt(765) = rrt(758)
  rrt(766) = rrt(758)
  rrt(767) = rrt(758)
  rrt(768) = rrt(758)
  rrt(769) = rrt(758)
  rrt(770) = rrt(758)
  rrt(771) = rrt(758)
  rrt(772) = rrt(758)
  rrt(773) = rrt(758)
  rrt(774) = rrt(758)
  rrt(775) = rrt(758)
  rrt(776) = rrt(758)
  rrt(777) = rrt(758)
  rrt(778) = rrt(758)
  rrt(779) = rrt(758)
  rrt(780) = rrt(758)
  rrt(781) = rrt(758)
  rrt(782) = rrt(758)
  rrt(783) = rrt(758)
  rrt(784) = rrt(758)
  rrt(785) = rrt(758)
  rrt(786) = rrt(758)
  rrt(787) = rrt(758)
  rrt(788) = rrt(758)
  rrt(789) = rrt(758)
  rrt(790) = rrt(758)
  rrt(791) = rrt(758)
  rrt(792) = rrt(758)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
