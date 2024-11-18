!+
! Function tao_curve_ele_ref (curve, point_to_track) result (ele_ref)
!
! Routine to point to the the corresponding curve%ref_ele.
! If point_to_track_ele is True and this element is a super_lord, ele_ref will point to the last super_slave.
!
! Input:
!   curve               -- tao_curve_struct: Curve with ref ele.
!   point_to_track_ele  -- logical: If True, point to tracking element. If False, 
!
! Output:
!   ele_ref     -- ele_struct: Corresponding element in the tracking 
!                         part of the lattice. Set to null.
!-

function tao_curve_ele_ref (curve, point_to_track) result (ele_ref)

use tao_interface, dummy => tao_curve_ele_ref

implicit none

type (tao_curve_struct) curve
type (ele_struct), pointer :: ele_ref
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele, slave
type (ele_pointer_struct), allocatable, target :: eles(:)

integer ix_branch, ix_ele
integer i_uni, ix_c

logical point_to_track, err

!

i_uni = tao_universe_index(curve%ix_universe)
lat => s%u(i_uni)%model%lat
call tao_locate_elements (curve%ele_ref_name, i_uni, eles, err, ignore_blank = .true., ix_branch = curve%ix_branch)

ele_ref => null()
if (size(eles) == 0) return
ele => eles(1)%ele

if (point_to_track) then
  if (ele%lord_status == super_lord$) then
    ele_ref => pointer_to_slave (ele, ele%n_slave)

  elseif (ele%lord_status == not_a_lord$) then
    ele_ref => ele
  endif

else
  ele_ref => ele
endif

end function
