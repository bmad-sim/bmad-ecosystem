!+
! Function congruent_lattice_elements (ele1, ele2) result (is_congruent)
!
! Routine to determine if two lattice elements are congruent. Two lattice elements are congruent
! If all the parameters of the two elements that are independent variables have the same value.
!
! For example, if two quadrupole with %field_master = True only have differing reference energy
! and differing k1 then the two quadrupoles are congruent.
!
! Input:
!   ele1  -- ele_struct: A Lattice element
!   ele2  -- ele_struct: A second lattice element to compare to.
!
! Output:
!   is_congruent  -- logical: Set True if ele1 is congruent to ele2. False otherwise.
!-

function congruent_lattice_elements (ele1, ele2) result (is_congruent)

use equality_mod
use attribute_mod, dummy2 => congruent_lattice_elements

implicit none

type (ele_struct) ele1, ele2
type (ele_attribute_struct) info

integer i
logical is_congruent

!

is_congruent = .true.

do i = 1, num_ele_attrib$
  if (ele1%value(i) == ele2%value(i)) cycle
  info = attribute_info(ele1, i)
  select case (info%state)
  case (is_free$)
    is_congruent = .false.
    return
  case (quasi_free$)
    if (attribute_free(ele1, info%name)) then
      is_congruent = .false.
      return
    endif
  end select
  return
enddo

!

if (attribute_name(ele1, tracking_method$) == 'TRACKING_METHOD') then
  if (ele1%tracking_method /= ele2%tracking_method) is_congruent = .false.
endif

if (attribute_name(ele1, mat6_calc_method$) == 'MAT6_CALC_METHOD') then
  if (ele1%mat6_calc_method /= ele2%mat6_calc_method) is_congruent = .false.
endif

if (attribute_name(ele1, spin_tracking_method$) == 'SPIN_TRACKING_METHOD') then
  if (ele1%spin_tracking_method /= ele2%spin_tracking_method) is_congruent = .false.
endif

if (attribute_name(ele1, aperture_at$) == 'APERTURE_AT') then
  if (ele1%aperture_at /= ele2%aperture_at) is_congruent = .false.
endif

if (attribute_name(ele1, aperture_type$) == 'APERTURE_TYPE') then
  if (ele1%aperture_type /= ele2%aperture_type) is_congruent = .false.
endif

if (attribute_name(ele1, field_calc$) == 'FIELD_CALC') then
  if (ele1%field_calc /= ele2%field_calc) is_congruent = .false.
endif

if (attribute_name(ele1, is_on$) == 'IS_ON') then
  if (ele1%is_on .neqv. ele2%is_on) is_congruent = .false.
endif

if (attribute_name(ele1, csr_method$) == 'CSR_METHOD') then
  if (ele1%csr_method /= ele2%csr_method) is_congruent = .false.
endif

if (attribute_name(ele1, space_charge_method$) == 'SPACE_CHARGE_METHOD') then
  if (ele1%space_charge_method /= ele2%space_charge_method) is_congruent = .false.
endif

!

if (associated(ele1%grid_field) .eqv. associated(ele2%grid_field)) then
  if (associated(ele1%grid_field)) then
    if (size(ele1%grid_field) == size(ele2%grid_field)) then
      do i = 1, size(ele1%grid_field)
        if (.not. ele1%grid_field(i) == ele2%grid_field(i)) is_congruent = .false.
      enddo
    else
      is_congruent = .false.
    endif
  endif
else
  is_congruent = .false.
endif

!

if (associated(ele1%cartesian_map) .eqv. associated(ele2%cartesian_map)) then
  if (associated(ele1%cartesian_map)) then
    if (size(ele1%cartesian_map) == size(ele2%cartesian_map)) then
      do i = 1, size(ele1%cartesian_map)
        if (.not. ele1%cartesian_map(i) == ele2%cartesian_map(i)) is_congruent = .false.
      enddo
    else
      is_congruent = .false.
    endif
  endif
else
  is_congruent = .false.
endif

!

if (associated(ele1%cylindrical_map) .eqv. associated(ele2%cylindrical_map)) then
  if (associated(ele1%cylindrical_map)) then
    if (size(ele1%cylindrical_map) == size(ele2%cylindrical_map)) then
      do i = 1, size(ele1%cylindrical_map)
        if (.not. ele1%cylindrical_map(i) == ele2%cylindrical_map(i)) is_congruent = .false.
      enddo
    else
      is_congruent = .false.
    endif
  endif
else
  is_congruent = .false.
endif

!

if (associated(ele1%taylor_field) .eqv. associated(ele2%taylor_field)) then
  if (associated(ele1%taylor_field)) then
    if (size(ele1%taylor_field) == size(ele2%taylor_field)) then
      do i = 1, size(ele1%taylor_field)
        if (.not. ele1%taylor_field(i) == ele2%taylor_field(i)) is_congruent = .false.
      enddo
    else
      is_congruent = .false.
    endif
  endif
else
  is_congruent = .false.
endif



end function congruent_lattice_elements
