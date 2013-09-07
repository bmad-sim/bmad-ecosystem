!+
! Subroutine track1_backup (end_orb, ele, param, start_orb, err_flag)
!
! "Backup" tracking through a single element. 
! That is, the particle is going backward in time from the downstream end
! of the element back to the upstream end. 
!
! This is different from tracking a particle going forward in time in the 
! -s direction.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   end_orb   -- Coord_struct: Position at downstream end.
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!     %particle     -- Particle type
!
! Output:
!   start_orb  -- Coord_struct: Position at upstream end.
!   err_flag  -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine track1_backup (end_orb, ele, param, start_orb, err_flag)

use bmad, dummy => track1_backup

implicit none

type (ele_struct) ele
type (coord_struct) :: start_orb, end_orb, e_orb
type (lat_param_struct) param

real(rp) beta_ref, dt_ref, z_end, t_end
logical, optional :: err_flag

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

e_orb = end_orb 
z_end = end_orb%vec(5)
t_end = end_orb%t

! flip to reversed coords

e_orb%vec(2) = -e_orb%vec(2)
e_orb%vec(4) = -e_orb%vec(4)

! backup tracking

select case (ele%key)
case (elseparator$)
  ! Nothing to do
case (lcavity$, rfcavity$)
  if (ele%key == lcavity$) then
    ele%value(gradient$) = -ele%value(gradient$)
    ele%value(gradient_err$) = -ele%value(gradient_err$)
  else
    ele%value(voltage$) = -ele%value(voltage$)
  endif
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  dt_ref = ele%value(l$) / (c_light * beta_ref)
  e_orb%vec(5) = e_orb%vec(5) - dt_ref * e_orb%beta
  e_orb%t = e_orb%t + dt_ref
case default
  param%rel_tracking_charge = -param%rel_tracking_charge
end select



ele%orientation = -ele%orientation
param%reverse_time_tracking = .true.

call track1 (e_orb, ele, param, start_orb, err_flag = err_flag)

ele%orientation = -ele%orientation
param%reverse_time_tracking = .false.

select case (ele%key)
case (elseparator$)
  ! Nothing to do
case (lcavity$)
  ele%value(gradient$) = -ele%value(gradient$)
  ele%value(gradient_err$) = -ele%value(gradient_err$)
case (rfcavity$)
  ele%value(voltage$) = -ele%value(voltage$)
case default
  param%rel_tracking_charge = -param%rel_tracking_charge
end select

! flip back to normal coords

start_orb%vec(2) = -start_orb%vec(2)
start_orb%vec(4) = -start_orb%vec(4)
start_orb%vec(5) = z_end - (start_orb%vec(5) - e_orb%vec(5))
start_orb%t = t_end - (start_orb%t - e_orb%t)

end subroutine track1_backup
