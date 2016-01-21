!+
! Function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)
!
! Function to convert a character string  (eg: "drift") to an index (eg: drift$).
! Wildcard "*" is translated to key_index = 0 (match all)
!
! Modules needed:
!   use bmad
!
! Input:
!   key_str        -- Character(*): Name of the key. Result is case insensitive.
!   abbrev_allowed -- Logical, optional: Abbreviations (eg: "quad") allowed?
!                       Default is False. At least 3 characters are needed 
!                       (except for rfcavity elements) if True.
!
! Output:
!   key_index -- Integer: Index of the key. Set to -1 if key_name not recognized.
!-

function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)

use basic_attribute_mod, dummy => key_name_to_key_index

implicit none

character(*) key_str
character(16) name

logical, optional :: abbrev_allowed
logical abbrev

integer key_index
integer i, n_name, n_match

!

if (key_str == '*') then
  key_index = 0
  return
endif

n_match = 0
key_index = -1
if (key_str == '') return

!

call str_upcase (name, key_str)
call string_trim (name, name, n_name)

abbrev = logic_option(.false., abbrev_allowed)

do i = 1, n_key$
  if (abbrev .and. (n_name > 2 .or. name(1:2) == "RF")) then
    if (name(:n_name) == key_name(i)(1:n_name)) then
      key_index = i
      n_match = n_match + 1
    endif
  else
    if (name == key_name(i)) then
      key_index = i
      return
    endif
  endif
enddo

if (abbrev .and. n_match > 1) key_index = -1  ! Multiple matches are not valid

end function key_name_to_key_index 

