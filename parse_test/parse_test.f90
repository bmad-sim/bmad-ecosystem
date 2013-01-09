!+
! Program parse_test
!
! This program is part of the Bmad regression testing suite.
!-

program parse_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

! 

call bmad_parser ('parse_test.bmad', lat)


open (1, file = 'output.now')

write (1, '(a, f10.4)') '"1 REL"  1E-10', lat%ele(1)%value(k1$)
write (1, '(a, f10.4)') '"2 REL"  1E-10', lat%ele(2)%value(k1$)
write (1, '(a, f10.4)') '"3 REL"  1E-10', lat%ele(3)%value(k1$)
write (1, '(a, f10.4)') '"4 REL"  1E-10', lat%ele(4)%value(k1$)
write (1, '(a, f10.4)') '"5 REL"  1E-10', lat%ele(5)%value(k2$)
write (1, '(2a)')       '"TM1"     STR ', tracking_method_name(lat%ele(1)%tracking_method)
write (1, '(2a)')       '"TM5"     STR ', tracking_method_name(lat%ele(5)%tracking_method)

close(1)

end program
