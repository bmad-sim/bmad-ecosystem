!+
! Subroutine track1_bmad_photon (start_orb, ele, param, end_orb, err_flag)
!
! Particle tracking through a single element BMAD_standard style.
! This routine is NOT meant for long term tracking since it does not get 
! all the 2nd order terms for the longitudinal motion.
!
! Note: track1_bmad_photon *never* relies on ele%mat6 for tracking excect for 
! hybrid elements.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   end_orb   -- Coord_struct: End position
!   err_flag  -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine track1_bmad_photon (start_orb, ele, param, end_orb, err_flag)

use capillary_mod, dummy => track1_bmad_photon
use track1_photon_mod, dummy2 => track1_bmad_photon
use lat_geometry_mod, dummy4 => track1_bmad_photon

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb, temp_orb
type (ele_struct) :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

real(rp) length, w_mat_inv(3,3), vec0(6), mat6(6,6), r_vec(3), p_vec(3)
real(rp) vel_vec(3), hit_point(3), cos_g, sin_g

integer i, n, n_slice, key

logical, optional :: err_flag
logical err

character(16) :: r_name = 'track1_bmad_photon'

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

start2_orb = start_orb ! In case start_orb and end_orb share the same memory.

end_orb = start_orb     ! transfer start to end
length = ele%value(l$)

!-----------------------------------------------
! Select
! If element is off looks like a drift. LCavities will still do wakefields.

key = ele%key
if (.not. ele%is_on .and. key /= lcavity$) key = drift$  

select case (key)

!-----------------------------------------------
! capillary

case (capillary$) 

  call offset_photon (ele, end_orb, set$)
  call track_a_capillary (end_orb, ele)
  call offset_photon (ele, end_orb, unset$)  

!-----------------------------------------------
! crystal

case (crystal$) 

  call offset_photon (ele, end_orb, set$)
  call track1_crystal (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$)

!-----------------------------------------------
! drift
 
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

  call offset_photon (ele, end_orb, set$)
  call track_a_drift_photon (end_orb, ele, length)
  call offset_photon (ele, end_orb, unset$)

!-----------------------------------------------
! marker, etc.

case (marker$, branch$, photon_branch$, floor_shift$, fiducial$)

  return

!-----------------------------------------------
! match

case (match$)

  if (ele%value(match_end_orbit$) /= 0) then
    ele%value(x0$)  = start2_orb%vec(1)
    ele%value(px0$) = start2_orb%vec(2)
    ele%value(y0$)  = start2_orb%vec(3)
    ele%value(py0$) = start2_orb%vec(4)
    ele%value(z0$)  = start2_orb%vec(5)
    ele%value(pz0$) = start2_orb%vec(6)
    end_orb%vec = [ ele%value(x1$), ele%value(px1$), &
                ele%value(y1$), ele%value(py1$), &
                ele%value(z1$), ele%value(pz1$) ]
    return
  endif

  call match_ele_to_mat6 (ele, vec0, mat6, err)
  if (err) then
    ! Since there are cases where this error may be raised many 
    ! times, do not print an error message.
    if (present(err_flag)) err_flag = .true.
    end_orb%state = lost$
    return
  endif

  end_orb%vec = matmul (mat6, end_orb%vec) + vec0
  end_orb%t = start2_orb%t + (ele%value(l$) + start2_orb%vec(5) - end_orb%vec(5)) / (c_light)
  end_orb%s = ele%s

!-----------------------------------------------
! mirror

case (mirror$)

  call offset_photon (ele, end_orb, set$)

  call to_crystal_surface_coords (ele, ele%value(graze_angle$), end_orb, vel_vec, hit_point, cos_g, sin_g)

  ! Check aperture

  if (ele%aperture_at == surface$) then
    temp_orb%vec(1:5:2) = hit_point 
    call check_aperture_limit (temp_orb, ele, surface$, param)
    if (end_orb%state /= alive$) return
  endif

  ! Reflect

  if (has_curved_surface(ele)) then
    call err_exit
  else
    end_orb%vec(1:4) = [-end_orb%vec(1), -end_orb%vec(2), end_orb%vec(3), end_orb%vec(4)]
  endif

  call offset_photon (ele, end_orb, unset$)

!-----------------------------------------------
! multilayer_mirror

case (multilayer_mirror$) 

  call offset_photon (ele, end_orb, set$)
  call track1_multilayer_mirror (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$)

!-----------------------------------------------
! patch

case (patch$)

  end_orb%vec(1) = end_orb%vec(1) - ele%value(x_offset$)
  end_orb%vec(3) = end_orb%vec(3) - ele%value(y_offset$)
  r_vec = [end_orb%vec(1), end_orb%vec(3), -ele%value(z_offset$)]
  
  p_vec = end_orb%vec(2:6:2)

  if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
    call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = w_mat_inv)
    p_vec = matmul(w_mat_inv, p_vec)
    r_vec = matmul(w_mat_inv, r_vec)
    end_orb%vec(2) = p_vec(1)
    end_orb%vec(4) = p_vec(2)
    end_orb%vec(6) = p_vec(3)
  endif

  end_orb%vec(1) = r_vec(1) - r_vec(3) * p_vec(1) / p_vec(3)
  end_orb%vec(3) = r_vec(2) - r_vec(3) * p_vec(2) / p_vec(3)

  end_orb%vec(5) = end_orb%vec(5) + r_vec(3) / p_vec(3) + &
                    end_orb%beta * (ele%value(z_offset$) + c_light * ele%value(t_offset$))

!-----------------------------------------------
! Taylor

case (taylor$)

  call track1_taylor (start_orb, ele, param, end_orb)
  end_orb%t = start2_orb%t + (ele%value(l$) + start2_orb%vec(5) - end_orb%vec(5)) / (c_light)
  end_orb%s = ele%s

end select

end subroutine track1_bmad_photon
