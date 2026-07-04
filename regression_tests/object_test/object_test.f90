!+
! Module object_test_priv
!
! Private module holding the get_more_text_func callback used by the object_test program. Made a
! module procedure (not nested) so that gfortran does not generate a stack trampoline (which would
! force an executable stack). See issue #1960.
!-

module object_test_priv

use object_model_mod

implicit none

contains

function get_more_text_func (line, end_of_doc, why_invalid) result (valid)

character(:), allocatable :: line, why_invalid
logical end_of_doc
logical valid

call str_set (line, '{a = {b = 7}}')
valid = .true.

end function

end module object_test_priv

!-------------------------------------------------------------------------------------

program test

use object_model_mod
use object_test_priv

implicit none

type (object_struct) obj
character(:), allocatable :: line, why_invalid
logical valid

valid = object_document_parse (obj, 'head', line, why_invalid, get_more_text_func)

print *, object_tree_name(obj%child(1)%child(1))

end program
