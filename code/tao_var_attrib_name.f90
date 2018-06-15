!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_var_attrib_name(var) result (var_attrib_name)
!
! Function to return a string encoding the attributes controlled by a 
! variable in the form:
!   {universe@}element[attribute]
! For example:
!   Q03W[k1]
!   [2,3]@V11E[HKICK]
!
! Input:
!   var -- Tao_var_struct: Variable
!
! Output:
!   var_attrib_name -- Character(60): Attribute list.
!-

function tao_var_attrib_name(var) result (var_attrib_name)

use location_encode_mod
use tao_struct

implicit none

type (tao_var_struct) var

character(60) var_attrib_name
integer i

!

if (size(s%u) > 1 .and. size(var%slave) > 0) then
  call location_encode (var_attrib_name, var%slave%ix_uni, ',')
  if (index(var_attrib_name, ',') /= 0 .or. index(var_attrib_name, ':') /= 0) then  
    var_attrib_name = '[' // trim(var_attrib_name) // ']' 
  endif
  var_attrib_name = trim(var_attrib_name) // '@' // &
                    trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
else
  var_attrib_name = trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
endif

end function tao_var_attrib_name
