program particle_species_test

use particle_species_mod

implicit none

integer :: i_dim, species
integer :: namelist_file, n_char

character(100) :: lat_name, lat_path, base_name, in_file
character(30), parameter :: r_name = 'particle_species_test'


character(20) :: example_names(1:100)  = ''

logical :: verbose

namelist / particle_species_test_params / &
    example_names, verbose

!------------------------------------------
! Defaults for namelist
example_names(1:11) =  [character(20) :: 'Ag+76', 'NH3--', 'CH3++', 'CH3+2', 'NH3@M37.5-', 'C+', '#12C+3', &
                        '#293Lv-2', '#291Og', '#295Ts', '#400Mc']
verbose = .false.

! Read namelist

if (command_argument_count() > 0) then
  call get_command_argument(1, in_file)
  namelist_file = lunget()
  print *, 'Opening: ', trim(in_file)
  open (namelist_file, file = in_file, status = "old")
  read (namelist_file, nml = particle_species_test_params)
  close (namelist_file)
endif

! 
print *, 'atomic mass unit (eV): ', atomic_mass_unit

! Decode
! call print_species(1)

! Encode
! species = species_id('#2H--')
! call print_species(species)

open (1, file = 'output.now')

! species = species_id('H2O ')
! call print_species(species)
! call write_mass_and_charge(species)

call print_all(example_names, .true.)
call print_all(atomic_name, .true.)
call print_all(molecular_name, .true.)

call check_fundamental('ref_particle', ref_particle$, 0.0_rp, 0, 0.0_rp, anti_ref_particle$)
call check_fundamental('anti_ref_particle', anti_ref_particle$, 0.0_rp, 0, 0.0_rp, ref_particle$)
call check_fundamental('deuteron', deuteron$, m_deuteron, +1, anomalous_mag_moment_deuteron, anti_deuteron$)
call check_fundamental('pion0', pion_0$, m_pion_0, 0, 0.0_rp, pion_0$)
call check_fundamental('pion+', pion_plus$, m_pion_charged, +1, 0.0_rp, pion_minus$)
call check_fundamental('antimuon', antimuon$, m_muon, +1, anomalous_mag_moment_muon, muon$)
call check_fundamental('proton', proton$, m_proton, +1, anomalous_mag_moment_proton, antiproton$)
call check_fundamental('positron', positron$, m_electron, +1, anomalous_mag_moment_electron, electron$)
call check_fundamental('photon', photon$, 0.0_rp, 0, 0.0_rp, photon$)
call check_fundamental('electron', electron$, m_electron, -1, anomalous_mag_moment_electron, positron$)
call check_fundamental('antiproton', antiproton$, m_proton, -1, anomalous_mag_moment_proton, proton$)
call check_fundamental('muon', muon$, m_muon, -1, anomalous_mag_moment_muon, antimuon$)
call check_fundamental('pion-', pion_minus$, m_pion_charged, -1, 0.0_rp, pion_plus$)
call check_fundamental('anti_deuteron', anti_deuteron$, m_deuteron, -1, anomalous_mag_moment_deuteron, deuteron$)

call magnetic_moment(deuteron$, 1.0_rp, 4.330735094d-27)
call magnetic_moment(proton$, 0.5_rp, 1.41060679736d-26)
call magnetic_moment(electron$, 0.5_rp, -9.2847647043d-24)
call magnetic_moment(muon$, 0.5_rp, -4.49044830d-26)

!---------------------------------------------------------------------------
contains

subroutine magnetic_moment (species, spin, mu_meas)

real(rp) spin, g, mu_meas, mu_calc, factor, a_meas
integer species

!

g = 2 * anomalous_moment_of(species) + 2
factor = charge_of(species) * spin * h_bar_planck * e_charge * c_light**2 / (2 * mass_of(species))
mu_calc = g * factor
a_meas = (mu_meas / factor - 2.0_rp) / 2.0_rp
!print '(a12, 5es20.12)', species_name(species), g, mu_calc, mu_meas, mu_calc - mu_meas, a_meas

write (1, '(3a, i0, 3x, es14.6)') '"dMag-', trim(species_name(species)), '"   ABS 1E', nint(log10(abs(mu_meas)))-10, mu_calc - mu_meas

end subroutine magnetic_moment

!---------------------------------------------------------------------------
! contains

subroutine check_fundamental (name, id, m_part, charge, anom_mag_moment, anti_part)
character(*) name
integer id, charge, anti_part
real(rp) m_part, anom_mag_moment

write (1, '(3a, t30, a, 5l1, a)') '"', name, '"', 'STR   "', upcase(name) == upcase(species_name(id)), m_part == mass_of(id), &
               charge == charge_of(id), anom_mag_moment == anomalous_moment_of(id), anti_part == antiparticle(id), '"' 


end subroutine

!---------------------------------------------------------------------------
! contains

subroutine print_all(p_array, convert_to_amu)
integer :: i, species
character(*)  :: p_array(:)
logical convert_to_amu

!
write(1, *) ''

do i=lbound(p_array, 1), ubound(p_array, 1)
  if (trim(p_array(i)) =='') cycle ! Skip empty
  species = species_id(p_array(i))
  if (verbose) call print_species(species, convert_to_amu)
  call write_mass_and_charge(species, convert_to_amu)
enddo

end subroutine

!---------------------------------------------------------------------------
! contains

subroutine print_species(species, convert_to_amu)

integer :: species
logical convert_to_amu

write (*, '(a, i0)')    'species: ', species
write (*, '(a, a)')     'name:    ', species_name(species)
if (convert_to_amu) then
  write (*, '(a, f20.5)') 'mass:    ', mass_of(species)/atomic_mass_unit
else
  write (*, '(a, es20.10)') 'mass:    ', mass_of(species)
endif
write (*, '(a, i0)')    'charge   ', charge_of(species)
write (*,*) ''

end subroutine 

!---------------------------------------------------------------------------
! contains

subroutine write_mass_and_charge(species, convert_to_amu)

integer :: species
character(20) :: name
logical convert_to_amu

name = species_name(species)
if (convert_to_amu) then
  write (1, '(a, a, f16.9)') '"'//trim(name)//':mass"', ' ABS  1E-10 ', mass_of(species) / atomic_mass_unit
else
  write (1, '(a, a, es20.10)') '"'//trim(name)//':mass"', ' ABS  1E-10 ', mass_of(species)
endif
write (1, '(a, a, i6)') '"'//trim(name)//':charge"'   , ' ABS  0 ', charge_of(species)

end subroutine

end program
