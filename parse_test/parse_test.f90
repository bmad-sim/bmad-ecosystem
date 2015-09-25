!+
! Program parse_test
!
! This program is part of the Bmad regression testing suite.
!-

program parse_test

use bmad
use bmad_parser_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orbit(:)
real(rp) value
integer i
character(1) delim
logical err, delim_found

! 

open (1, file = 'output.now')

call bmad_parser ('parse_test.bmad', lat)

write (1, '(a, f10.2)') '"bmad_com[max_aperture_limit]"          ABS 0', bmad_com%max_aperture_limit
write (1, '(a, i4)')    '"bmad_com[ptc_max_fringe_order]"        ABS 0', bmad_com%ptc_max_fringe_order
write (1, '(a, l1, a)') '"bmad_com[convert_to_kinetic_momentum]" STR   "', bmad_com%convert_to_kinetic_momentum, '"'

bp_com%input_from_file = .false.
bp_com%parse_line = '-2*7)'
call evaluate_value ('ERR', value, lat, delim, delim_found, err, ',)')
write (1, '(a, f10.4)') '"EVAL 1"  ABS 0', value

write (1, '(a, f10.4)') '"1 REL"  ABS 0', lat%ele(1)%value(k1$)
write (1, '(a, f10.4)') '"2 REL"  ABS 0', lat%ele(2)%value(k1$)
write (1, '(a, f10.4)') '"3 REL"  ABS 0', lat%ele(3)%value(k1$)
write (1, '(a, f10.4)') '"4 REL"  ABS 0', lat%ele(4)%value(k1$)
write (1, '(a, f10.4)') '"5 REL"  ABS 0', lat%ele(5)%value(k2$)
write (1, '(4a)')       '"TM1"     STR ', '"', trim(tracking_method_name(lat%ele(1)%tracking_method)), '"'
write (1, '(4a)')       '"TM5"     STR ', '"', trim(tracking_method_name(lat%ele(5)%tracking_method)), '"'
write (1, '(a, i3)')       '"N3"     ABS 0', lat%branch(2)%n_ele_track
write (1, '(4a)')       '"Custom"  STR ', '"', trim(attribute_name(lat%ele(1), custom_attribute1$)), '"'
do i = lbound(mass_of, 1), ubound(mass_of, 1)
  write (1, '(3a, es20.12, i6)')         '"', trim(particle_name(i)), '"  REL 1e-12', mass_of(i), charge_of(i)
enddo

!

bp_com%always_parse = .true.
call bmad_and_xsif_parser ('DCO4.xsif', lat)
lat%ele(156)%value(hkick$) = 0.00001
lat%ele(156)%value(vkick$) = 0.00002

call twiss_at_start (lat)
call twiss_propagate_all (lat)
call reallocate_coord (orbit, lat)
call closed_orbit_calc (lat, orbit, 4)

write (1, '(a, 2f12.8)')  '"XSIF:Twiss"  REL 1e-8', lat%ele(0)%a%beta, lat%ele(0)%b%beta
write (1, '(a, 2es15.8)') '"XSIF:Orbit"  REL 1e-8', orbit(0)%vec(1), orbit(0)%vec(3)

!

call bmad_parser ('parse_test.bmad', lat, use_line = 'l2')

write (1, '(4a)')       '"L2-1"     STR ', '"', trim(lat%ele(1)%name), '"'
write (1, '(2a, i0, a)')       '"L2-N"     STR ', '"', lat%n_ele_max, '"'


close(1)

end program
