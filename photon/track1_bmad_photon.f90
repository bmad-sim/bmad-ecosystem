!+
! Subroutine track1_bmad_photon (orbit, ele, param, err_flag)
!
! Particle tracking through a single element BMAD_standard style.
! This routine is NOT meant for long term tracking since it does not get 
! all the 2nd order terms for the longitudinal motion.
!
! Note: track1_bmad_photon *never* relies on ele%mat6 for tracking excect for 
! hybrid elements.
! 
! Input:
!   orbit      -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   orbit     -- Coord_struct: End position
!   err_flag  -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine track1_bmad_photon (orbit, ele, param, err_flag)

use capillary_mod, dummy => track1_bmad_photon
use track1_photon_mod, dummy2 => track1_bmad_photon
use bmad_interface, dummy4 => track1_bmad_photon

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct) :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

real(rp) w(3,3), vec0(6), mat6(6,6)
real(rp) vel_vec(3), hit_point(3), cos_g, sin_g

integer i, n, n_slice, key

logical, optional :: err_flag
logical err

character(*), parameter :: r_name = 'track1_bmad_photon'

!-----------------------------------------------
! Select
! If element is off looks like a drift. LCavities will still do wakefields.

if (present(err_flag)) err_flag = .true.

start_orb = orbit
orbit%vec(5) = 0      ! start at beginning of element

key = ele%key
if (.not. ele%is_on) key = drift$  

select case (key)

!-----------------------------------------------
! Capillary

case (capillary$) 

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track_a_capillary (orbit, ele)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Crystal

case (crystal$) 

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  if (is_true(ele%value(is_mosaic$))) then
    call track1_mosaic_crystal (ele, param, orbit)
  else
    call track1_crystal (ele, param, orbit)
  endif
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Diffraction_plate
 
case (diffraction_plate$)

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_diffraction_plate_or_mask (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Drift
 
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

  if (orbit%vec(6) * orbit%direction < 0) then  ! Heading backwards
    orbit%state = lost_pz_aperture$
    return
  endif

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track_a_drift_photon (orbit, ele%value(l$), .true.)
  if (orbit%state /= alive$) return
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

  orbit%s = ele%s

!-----------------------------------------------
! Lens

case (lens$) 

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_lens (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Marker, etc.

case (marker$, detector$, fork$, photon_fork$, floor_shift$, fiducial$)

  orbit%vec(5) = 0

!-----------------------------------------------
! Mask                                                                                                                         

case (mask$)

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_diffraction_plate_or_mask (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return
   
!-----------------------------------------------
! Match

case (match$)

  if (is_true(ele%value(match_end_orbit$))) then
    ele%value(x0$)  = start_orb%vec(1)
    ele%value(px0$) = start_orb%vec(2)
    ele%value(y0$)  = start_orb%vec(3)
    ele%value(py0$) = start_orb%vec(4)
    ele%value(z0$)  = start_orb%vec(5)
    ele%value(pz0$) = start_orb%vec(6)
    orbit%vec = [ele%value(x1$), ele%value(px1$), &
                   ele%value(y1$), ele%value(py1$), &
                   ele%value(z1$), ele%value(pz1$)]

  else
    call track1_bmad (orbit, ele, param, err_flag)
  endif

!-----------------------------------------------
! Mirror

case (mirror$)

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_mirror (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Multilayer_Mirror

case (multilayer_mirror$) 

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_multilayer_mirror (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Patch

case (patch$)

  call track_a_patch_photon (ele, orbit)

!-----------------------------------------------
! Sample

case (sample$)

  call offset_photon (ele, orbit, set$); if (orbit%state /= alive$) return
  call track1_sample (ele, param, orbit)
  call offset_photon (ele, orbit, unset$); if (orbit%state /= alive$) return

!-----------------------------------------------
! Taylor

case (taylor$)

  call track1_taylor (orbit, ele, param)
  orbit%t = start_orb%t + (ele%value(l$) + start_orb%vec(5) - orbit%vec(5)) / (c_light)
  orbit%s = ele%s

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

orbit%location = downstream_end$
if (present(err_flag)) err_flag = .false.

end subroutine track1_bmad_photon
