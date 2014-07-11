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
use geometry_mod, dummy4 => track1_bmad_photon

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb, temp_orb
type (ele_struct) :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

real(rp) length, w(3,3), vec0(6), mat6(6,6)
real(rp) vel_vec(3), hit_point(3), cos_g, sin_g

integer i, n, n_slice, key

logical, optional :: err_flag
logical err

character(*), parameter :: r_name = 'track1_bmad_photon'

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

start2_orb = start_orb ! In case start_orb and end_orb share the same memory.

end_orb = start_orb     ! transfer start to end
length = ele%value(l$)

!-----------------------------------------------
! Select
! If element is off looks like a drift. LCavities will still do wakefields.

key = ele%key
if (.not. ele%is_on) key = drift$  

select case (key)

!-----------------------------------------------
! Capillary

case (capillary$) 

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track_a_capillary (end_orb, ele)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Crystal

case (crystal$) 

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_crystal (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Diffraction_plate
 
case (diffraction_plate$)

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_diffraction_plate (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Drift
 
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

  if (end_orb%vec(6) * end_orb%direction < 0) then  ! Heading backwards
    end_orb%state = lost_z_aperture$
    return
  endif

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track_a_drift_photon (end_orb, length, .true.); if (end_orb%state /= alive$) return
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

  end_orb%s = ele%s

!-----------------------------------------------
! Marker, etc.

case (marker$, x_ray_init$, detector$, fork$, photon_fork$, floor_shift$, fiducial$)

  end_orb%vec(5) = 0
  return

!-----------------------------------------------
! Match

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
! Mirror

case (mirror$)

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_mirror (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Multilayer_Mirror

case (multilayer_mirror$) 

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_multilayer_mirror (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Patch

case (patch$)

  if (end_orb%direction == 1) then
    ! Translate (x, y, z) to coordinate system with respect to downstream origin.
    end_orb%vec(1) = end_orb%vec(1) - ele%value(x_offset$)
    end_orb%vec(3) = end_orb%vec(3) - ele%value(y_offset$)
    end_orb%vec(5) = -ele%value(z_offset$)   
    
    if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = w)
      end_orb%vec(2:6:2) = matmul(w, end_orb%vec(2:6:2))
      end_orb%vec(1:5:2) = matmul(w, end_orb%vec(1:5:2))
    endif

    call track_a_drift_photon (end_orb, -end_orb%vec(5), .false.)
    end_orb%s = ele%s

  else
    end_orb%s = ele%s - ele%value(l$)
    end_orb%vec(5) = 0   ! Assume particle starts at downstream face

    if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat = w)
      end_orb%vec(2:6:2) = matmul(w, end_orb%vec(2:6:2))
      end_orb%vec(1:5:2) = matmul(w, end_orb%vec(1:5:2))
    endif

    end_orb%vec(1) = end_orb%vec(1) + ele%value(x_offset$)
    end_orb%vec(3) = end_orb%vec(3) + ele%value(y_offset$)

    call track_a_drift_photon (end_orb, ele%value(z_offset$), .false.)

  endif

!-----------------------------------------------
! Sample

case (sample$)

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_sample (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Taylor

case (taylor$)

  call track1_taylor (start_orb, ele, param, end_orb)
  end_orb%t = start2_orb%t + (ele%value(l$) + start2_orb%vec(5) - end_orb%vec(5)) / (c_light)
  end_orb%s = ele%s

!-----------------------------------------------
! Not recognized

case default

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'BMAD_STANDARD TRACKING_METHOD NOT IMPLMENTED FOR PHOTONS FOR: ' // key_name(ele%key), &
          'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

end select

end subroutine track1_bmad_photon
