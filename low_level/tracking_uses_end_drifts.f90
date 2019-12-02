!+
! Function tracking_uses_end_drifts (ele, use_hard_edge_model) result (has_drifts)
!
! Function to determine if the tracking for an element uses a "hard edge model"
! where the tracking looks like (drift, model, drift). For example,
! RF cavity fields with ele%field_calc = bmad_standard$ use a hard edge model where
! the length of the cavity is c_light / (2 * freq).
!
! If the use_hard_edge model argument is set to False, this function will always return False.
!
! Input:
!   ele                   -- ele_struct: Element.
!   use_hard_edge_drifts  -- logical, optional: Use bmad hard edge model for cavities, etc?
!                              Default is set by bmad_com%use_hard_edge_drifts.
!                              Default bmad_com%use_hard_edge_drifts is True.
!
! Output:
!   has_drifts -- Logical: True if tracking uses end drifts.
!-

function tracking_uses_end_drifts (ele, use_hard_edge_model) result (has_drifts)

use bmad_struct
implicit none

type (ele_struct) ele
logical has_drifts
logical, optional :: use_hard_edge_model

!

has_drifts = .false.
if (.not. logic_option(bmad_com%use_hard_edge_drifts, use_hard_edge_model)) return

select case (ele%key)
case (lcavity$, rfcavity$)
  if (ele%field_calc == bmad_standard$) has_drifts = .true.
  if (ele%value(l_hard_edge$) >= ele%value(l$)) has_drifts = .false. ! Avoid creating negative length end drifts 
end select

end function tracking_uses_end_drifts
