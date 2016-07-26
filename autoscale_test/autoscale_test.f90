program test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

!

call bmad_parser ('autoscale.bmad', lat)
open (1, file = 'output.now')

ele => lat%ele(1)
write (1, '(a, f16.12)') '"Auto_Phase" ABS 2e-4', ele%value(phi0_autoscale$)
write (1, '(a, es18.10)') '"Auto_Amp"   REL 1e-6', ele%value(field_autoscale$)

end program
