!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)
!
! Routine to split a string in the form str = "xxx::yyy" into an element class
! and an element name. Example: 
!   str = "quad::q*".
! gives
!   ix_class = quadrupole$
!   ele_name = "Q*"
!
! If "::" is not found then ix_class is set to 0 (all classes).
! If str is of the form: "*::yyy" then ix_class is set to 0 (all classes).
! Class abbreviations and lower case names accepted 
! ele_name will be converted to upper case
!
! Input:
!   str       -- Character(*): Character string to parse.
!   print_err -- Logical, optional: If True then print an error message if 
!                 there is a problem. Default is True.
!
! Output:
!   ix_class  -- Integer: Element class. 0 => all classes. -1 => Not an element [EG: str = "var::..."].
!   ele_name  -- Character(*): Element name.
!   err       -- Set true if there is a problem translating the element class.
!-

subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)

use bmad_routine_interface

implicit none

integer ix, ix_class

character(*) str, ele_name
character(40) :: r_name = 'tao_string_to_element_id'
character(20) class

logical, optional :: print_err
logical err

! 

if (index(str, 'dat::') /= 0) then
  call out_io (s_error$, r_name, 'NAME USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // str)
  call err_exit
endif

!

err = .false.
ix_class = -1

if (str(1:6) == 'data::') return
if (str(1:5) == 'var::') return
if (str(1:5) == 'lat::') return
if (str(1:6) == 'wall::') return

ix = index(str, '::')

if (ix == 0) then
  ix_class = 0
  ele_name = str
  call str_upcase (ele_name, ele_name)
  return
endif

class = str(:ix-1)
ele_name = str(ix+2:)
call str_upcase (ele_name, ele_name)

if (class == '*') then
  ix_class = 0
  return
endif

ix_class = key_name_to_key_index (class, .true.)
if (ix_class < 1) then
  if (logic_option (.true., print_err)) call out_io (s_error$, r_name, 'BAD CLASS NAME: ' // class)
  err = .true.
endif

end subroutine tao_string_to_element_id
