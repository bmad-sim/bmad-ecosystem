subroutine parse_converter_distribution(ele, delim, delim_found)

use bmad_parser_mod
use object_model_mod

implicit none

type (ele_struct) ele
type (object_struct) obj

integer i, ix_word

character(*) delim
character(200) str
character(100) why_invalid
character(:), allocatable :: line

logical delim_found, valid

!

valid = object_document_parse (obj, line, why_invalid, get_more_text_func, parse_one = .true.)
bp_com%parse_line = trim(line) // ' ' // bp_com%parse_line
call get_next_word(str, ix_word, '}],', delim, delim_found)

if (.not. valid) then
  call parser_error (why_invalid)
  return
endif

do i = 1, obj%n_obj
  print *, obj%obj(i)%name
enddo

!------------------------
contains

function get_more_text_func (line, why_invalid) result (valid)

integer ix_word
character(:), allocatable :: line
character(*) why_invalid
character(200) str
logical valid

!

call get_next_word(str, ix_word, '}],', delim, delim_found)
valid = (len(str) /= 0 .and. .not. delim_found)
if (.not. valid) why_invalid = 'END OF ELEMENT DEFINITION BEFORE PARSING COMPLETE.'
line = trim(line) // ' ' // trim(str) // delim


end function get_more_text_func

end subroutine parse_converter_distribution

