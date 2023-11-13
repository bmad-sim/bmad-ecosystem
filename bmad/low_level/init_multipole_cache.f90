!+
! Subroutine init_multipole_cache(ele)
!
! Routine to init ele%multipole_cache component.
!
! Input:
!   ele       -- ele_struct: Element to init
!
! Output:
!   ele       -- ele_struct: Initalized element.
!-

subroutine init_multipole_cache(ele)

use attribute_mod, dummy => init_multipole_cache

implicit none

type (ele_struct) ele

!

if (has_attribute(ele, 'A0') .or. has_attribute(ele, 'A0_ELEC') .or. ele%key == multipole$) then
  if (allocated(ele%multipole_cache)) then
    ele%multipole_cache%mag_valid = .false.
    ele%multipole_cache%elec_valid = .false.
  else
    allocate(ele%multipole_cache)
  endif
endif


end subroutine
