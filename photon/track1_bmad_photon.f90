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
use bmad_interface, dummy4 => track1_bmad_photon

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb, temp_orb
type (ele_struct) :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

real(rp) w(3,3), vec0(6), mat6(6,6)
real(rp) vel_vec(3), hit_point(3), cos_g, sin_g

integer i, n, n_slice, key

logical, optional :: err_flag
logical err

character(*), parameter :: r_name = 'track1_bmad_photon'

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .true.

start2_orb = start_orb ! In case start_orb and end_orb share the same memory.

end_orb = start_orb     ! transfer start to end
end_orb%vec(5) = 0      ! start at beginning of element

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
  if (is_true(ele%value(is_mosaic$))) then
    call track1_mosaic_crystal (ele, param, end_orb)
  else
    call track1_crystal (ele, param, end_orb)
  endif
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Diffraction_plate
 
case (diffraction_plate$)

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_diffraction_plate_or_mask (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Drift
 
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

  if (end_orb%vec(6) * end_orb%direction < 0) then  ! Heading backwards
    end_orb%state = lost_pz_aperture$
    return
  endif

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track_a_drift_photon (end_orb, ele%value(l$), .true.)
  if (end_orb%state /= alive$) return
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

  end_orb%s = ele%s

!-----------------------------------------------
! Lens

case (lens$) 

  call offset_photon (ele, end_orb, set$); if (end_orb%state /= alive$) return
  call track1_lens (ele, param, end_orb)
  call offset_photon (ele, end_orb, unset$); if (end_orb%state /= alive$) return

!-----------------------------------------------
! Marker, etc.

case (marker$, detector$, fork$, photon_fork$, floor_shift$, fiducial$)

  end_orb%vec(5) = 0

!-----------------------------------------------
! Match

case (match$)

  if (is_true(ele%value(match_end_orbit$))) then
    ele%value(x0$)  = start2_orb%vec(1)
    ele%value(px0$) = start2_orb%vec(2)
    ele%value(y0$)  = start2_orb%vec(3)
    ele%value(py0$) = start2_orb%vec(4)
    ele%value(z0$)  = start2_orb%vec(5)
    ele%value(pz0$) = start2_orb%vec(6)
    end_orb%vec = [ele%value(x1$), ele%value(px1$), &
                   ele%value(y1$), ele%value(py1$), &
                   ele%value(z1$), ele%value(pz1$)]

  else
    call track1_bmad (start2_orb, ele, param, end_orb, err_flag)
  endif

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

  call track_a_patch_photon (ele, end_orb)

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
! Photon_Init

case (photon_init$)
  ! For tracking purposes, a phton_init element is like a marker

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

!---------------------------------------------------------------------------------------------------

end_orb%location = downstream_end$
if (present(err_flag)) err_flag = .false.

end subroutine track1_bmad_photon
