!+
! Subroutine element_locator (ele_name, lat, ix_ele)
!
! Subroutine to locate an element in a lattice. If multipole elements
! match then ix_ele will point to the first one. Also see:
!   elements_locator
! The element name may be of the form: "S:<number>". This will match
! to the element whose longitudinal position at the exit end is
! closest to <number>.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele_name -- Character(40): Name of the element to find.
!                 The name is case insensitive
!   lat      -- lat_struct: Lattice to search through.
!
! Output:
!   ix_ele -- Integer: Index of element in lat%ele(:) array. 
!               ix_ele set to -1 if not found.
!-

#include "CESR_platform.inc"

subroutine element_locator (ele_name, lat, ix_ele)

use bmad_struct
use bmad_interface, except_dummy => element_locator
  
implicit none

type (lat_struct) lat
integer ix_ele, ios
real(rp) s
character(*) ele_name
character(40) e_name

! Init

call str_upcase (e_name, ele_name)

! Longitudinal position search

if (e_name(1:2) == 'S:') then
  read (e_name(3:), *, iostat = ios) s
  if (ios /= 0) return
  call ele_at_s (lat, s, ix_ele)
  ! See if previous element is closer
  if (ix_ele > 0) then
    if (s - lat%ele(ix_ele-1)%s < lat%ele(ix_ele)%s - s) ix_ele = ix_ele - 1
  endif
  return
endif

! Name search

do ix_ele = 0, lat%n_ele_max
  if (lat%ele(ix_ele)%name == e_name) return
enddo

ix_ele = -1

end subroutine
