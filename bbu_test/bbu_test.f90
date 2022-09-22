program bbu_test

use bbu_track_mod

implicit none

type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in
type (beam_init_struct) beam_init

integer irep
real(rp) hom_voltage_gain, growth_rate
logical lost

integer status

type (coord_struct), allocatable :: orb(:) 
namelist / bbu_params / bbu_param, beam_init


! Defaults for namelist
beam_init%n_particle = 1

open (1, file = 'bbu.init', status = 'old')
read (1, nml = bbu_params)
close (1)

! Define distance between bunches
beam_init%dt_bunch = 1 / bbu_param%bunch_freq

! Seed random number generator
call ran_seed_put (bbu_param%ran_seed)
if (bbu_param%ran_gauss_sigma_cut > 0) then
  call ran_gauss_converter (set_sigma_cut = bbu_param%ran_gauss_sigma_cut)
endif

! Init and parse
call bmad_parser (bbu_param%lat_filename, lat_in) !! lat_in is the parsed lattice

!! Closed orbit computed
call twiss_and_track(lat_in,orb,status,0,.true.)

lat = lat_in  

call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)

open (2, file = 'output.now', recl = 200)

beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
write (2, '(a, l1, a)')     '"lost_boolean_1A"      STR "', lost, '"'
write (2, '(a, es22.12E3)') '"hom_voltage_gain_1A"  ABS 1E-8', hom_voltage_gain
write (2, '(a, es22.12E3)') '"growth_rate_1A"       ABS 1E-8', growth_rate
 
bbu_param%current = 0.001
beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
write (2, '(a, l1, a)')     '"lost_boolean_1mA"      STR "', lost, '"' 
write (2, '(a, es22.12E3)') '"hom_voltage_gain_1mA"  ABS 1E-8', hom_voltage_gain
write (2, '(a, es22.12E3)') '"growth_rate_1mA"       ABS 1E-8', growth_rate

bbu_param%current = 100
beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
write (2, '(a, l1, a)')     '"lost_boolean_100A"      STR "', lost, '"'
write (2, '(a, es22.12E3)') '"hom_voltage_gain_100A"  ABS 1E-8', hom_voltage_gain
write (2, '(a, es22.12E3)') '"growth_rate_100A"       ABS 4E-8', growth_rate

close(2)

end program
