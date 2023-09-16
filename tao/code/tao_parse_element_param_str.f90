!+
! Subroutine tao_parse_element_param_str (err, in_str, uni, element, parameter, where, component)
!
! Routine to parse apart lattice element parameter strings like: 
!   "4@ele_mid::q01[k1]|model".
!   ele_begin::quadrupole::q*[x_offset]
!
! Input:
!   in_str      -- character(*): String specifying a parameter of an element or elements.
!
! Output:
!   err         -- logical: Set True if there is a parse error. False otherwise.
!   uni         -- character(*): Universe substring.
!   element     -- character(*): Element name.
!   parameter   -- character(*): Element parameter name.
!   where       -- integer: One of not_set$, anchor_beginning$, anchor_center$, or anchor_end$.
!   component   -- character(*): One of "model", "design", or "base".
!-

subroutine tao_parse_element_param_str (err, in_str, uni, element, parameter, where, component)

use tao_interface, dummy => tao_parse_element_param_str

implicit none

character(*) in_str, uni, element, parameter, component
character(100) str

integer where, ix, ix1, ix2
logical err

!

uni = ''
element = ''
parameter = ''
where = not_set$
component = ''

str = in_str
err = .true.

ix = index(str, '@')
if (ix /= 0) then
  uni = str(:ix-1)
  str = str(ix+1:)
endif

ix = index(str, '::')
if (ix == 0) return
select case (str(:ix-1))
case ('ele');         where = anchor_end$
case ('ele_mid');     where = anchor_center$
case ('ele_begin');   where = anchor_beginning$
end select

if (where /= not_set$) str = str(ix+2:)

ix = index(str, '[')
ix1 = index(str, '|')

if (ix /= 0) then
  ix2 = index(str, ']')
  if (ix2 == 0) return
  element = str(:ix-1)
  parameter = str(ix+1: ix2-1)
  str = str(ix2+1:)
elseif (ix1 /= 0) then
  element = str(:ix1-1)
  str = str(ix1:)
endif

if (str == '') then
  err = .false.
  return
endif

if (str(1:1) /= '|') return
component = str
err = .false.

end subroutine tao_parse_element_param_str
