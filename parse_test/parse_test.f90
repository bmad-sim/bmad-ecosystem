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

real(rp) value
integer i
character(1) delim
logical err, delim_found

! 

open (1, file = 'output.now')

call bmad_parser ('parse_test.bmad', lat)

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
write (1, '(4a)')       '"Custom"  STR ', '"', trim(attribute_name(lat%ele(1), custom_attribute1$)), '"'
do i = lbound(particle_name, 1), ubound(particle_name, 1)
  write (1, '(3a, es20.12, i6)')         '"', trim(particle_name(i)), '"  REL 1e-12', mass_of(i), charge_of(i)
enddo

close(1)

end program
