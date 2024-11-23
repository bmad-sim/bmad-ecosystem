!+
! Subroutine ele_to_spin_taylor(ele, param, orb0)
!
! Routine to create a spin Taylor map.
!
! Input:
!   ele       -- ele_struct: Lattice element.
!   param     -- lat_param_sruct: Branch parameters.
!   orb0      -- coord_struct: Starting ref coords.
!
! Output:
!   ele       -- ele_struct: Element with spin map.
!-

subroutine ele_to_spin_taylor(ele, param, orb0)

use bmad

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb0, end_orb
integer i
logical st_on, err_flag

character(*), parameter :: r_name = 'ele_to_spin_taylor'
!

if (ele%spin_tracking_method == sprint$) then
  call sprint_spin_taylor_map(ele, orb0%vec)

elseif (valid_tracking_method(ele, orb0%species, symp_lie_ptc$)) then
  st_on = bmad_com%spin_tracking_on
  bmad_com%spin_tracking_on = .true.
  call ele_to_taylor(ele, orb0, include_damping = bmad_com%radiation_damping_on)
  bmad_com%spin_tracking_on = st_on

else
  call make_mat6_tracking(ele, param, orb0, end_orb, err_flag, spin_only = .true.)
endif

end subroutine
