! Function tao_constraint_type_name (datum) result (datum_name)
!
! Function to return the constraint type in the form:
!   data_type <merit_type>
! For example:
!   eta.x <target>
!
! Input:
!   datum      -- Tao_data_struct: Datum
!
! Output:
!   datum_name -- Character(400): Appropriate name.
!-

function tao_constraint_type_name(datum) result (datum_name)

use tao_struct

implicit none

type (tao_data_struct) datum
character(400) datum_name
integer ix

! Expressions strings can be long. So shorten the name a wee bit.

datum_name = trim(datum%data_type) // ' <' // trim(datum%merit_type) // '>'
if (datum_name(1:11) == 'expression:') datum_name = adjustl(datum_name(12:))

end function tao_constraint_type_name
