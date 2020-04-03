program calc_rop

  ! declare variables and modules
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_temperature_K  = 800.0d0,  & 
                                 density_ini_nh3   = 2.01d+18,  &
                                 density_ini_o2    = 1.51d+18,  &
                                 density_ini_n2    = 5.66d+18,  &
                                 density_ini_elec  = 9.17d-02,  &
                                 reduced_field_Td  = 100d0,     &
                                 spec_heat_ratio   = 7.d0/5.d0
  ! double precision            :: time  = 0.0d0, dtime = 1.0d-10
  integer                     :: i, j, unit_stm = 10
  ! logical                     :: gas_heating      = .true.

  double precision,dimension(species_max, reactions_max)  :: source_terms_matrix
  open(unit_stm, file='output/source_terms_matrix', form='formatted')

  call ZDPlasKin_init()

  ! Set intial conditions
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature_K, REDUCED_FIELD=reduced_field_Td )
  ! call ZDPlasKin_set_conditions(SPEC_HEAT_RATIO=spec_heat_ratio, GAS_HEATING=gas_heating )
  call ZDPlasKin_set_density( 'NH3', density_ini_nh3)
  call ZDPlasKin_set_density(  'O2', density_ini_o2)
  call ZDPlasKin_set_density(  'N2', density_ini_n2)
  call ZDPlasKin_set_density(   'e', density_ini_elec)
  call ZDPlasKin_set_density('N2^+', density_ini_elec) 
	

  ! Flag for qt-style output
  ! call ZDPlaskin_set_config(QTPLASKIN_SAVE=.true.)

  ! print column headers and initial values
  write(*,'(4(A12))') ( trim(species_name(i)), i = 1, species_max )
  write(*,'(4(1pe12.4))') density(:)

  call ZDPlasKin_get_rates(SOURCE_TERMS_MATRIX=source_terms_matrix)

  ! do i = 1, species_max
  !   write(unit_stm, *) (source_terms_matrix(i, j), j = 1, reactions_max)
  ! enddo

  print*, 'aaa'
  write(unit_stm, *) sum(source_terms_matrix, dim=2)

  ! call ZDPlasKin_timestep(time,dtime)

end program calc_rop
