! NOTE: THIS ROUTINE IS NOT CURRENTLY OPERATIONAL!
!       DO NOT USE!
!       IF YOU ARE INTERESTED IN USING THIS PLEASE CONTACT DAVID SAGAN.
!+
! Subroutine track1_backwards_time (start_orb, ele, param, end_orb, err_flag)
!
! Backwards tracking in time through a single element. 
! That is, the particle is going backward in time from the downstream end
! of the element back to the upstream end. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Position at downstream end.
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!     %particle     -- Particle type
!
! Output:
!   end_orb   -- Coord_struct: Position at upstream end.
!   err_flag  -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine track1_backwards_time (start_orb, ele, param, end_orb, err_flag)

use bmad, dummy => track1_backwards_time

implicit none

type (ele_struct) ele
type (coord_struct) :: end_orb, start_orb, e_orb
type (lat_param_struct) param

real(rp) beta_ref, dt_ref, z_end, t_end
logical, optional :: err_flag
character(*), parameter :: r_name = 'track1_backwards_time'

!

call out_io (s_abort$, r_name, 'ROUTINE IS NOT CURRENTLY OPERATIONAL!')
call err_exit

! initially set start_orb = end_orb

if (present(err_flag)) err_flag = .false.

e_orb = start_orb 
z_end = start_orb%vec(5)
t_end = start_orb%t

! flip to reversed coords

e_orb%vec(2) = -e_orb%vec(2)
e_orb%vec(4) = -e_orb%vec(4)

! backwards_time tracking

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
  e_orb%species = -e_orb%species
end select



ele%orientation = -ele%orientation
param%backwards_time_tracking = .true.

call track1 (e_orb, ele, param, end_orb, err_flag = err_flag)

ele%orientation = -ele%orientation
param%backwards_time_tracking = .false.

select case (ele%key)
case (elseparator$)
  ! Nothing to do
case (lcavity$)
  ele%value(gradient$) = -ele%value(gradient$)
  ele%value(gradient_err$) = -ele%value(gradient_err$)
case (rfcavity$)
  ele%value(voltage$) = -ele%value(voltage$)
end select

! flip back to normal coords

end_orb%vec(2) = -end_orb%vec(2)
end_orb%vec(4) = -end_orb%vec(4)
end_orb%vec(5) = z_end - (end_orb%vec(5) - e_orb%vec(5))
end_orb%t = t_end - (end_orb%t - e_orb%t)
end_orb%species = start_orb%species

end subroutine track1_backwards_time
