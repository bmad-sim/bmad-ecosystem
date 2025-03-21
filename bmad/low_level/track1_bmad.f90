!+
! Subroutine track1_bmad (start_orb, ele, param, err_flag, track, mat6, make_matrix)
!
! Particle tracking through a single element BMAD_standard style.
!
! Input:
!   orbit       -- coord_struct: Starting position
!   ele         -- ele_struct: Element
!   param       -- lat_param_struct:
!     %particle     -- Particle type
!   mat6(6,6)   -- real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit         -- coord_struct: End position.
!   err_flag      -- logical, optional: Set true if there is an error. False otherwise.
!   track         -- track_struct, optional: Structure holding the track information if the 
!                      lattice element does tracking step-by-step. See track1 for more details.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix propagated through the element.
!-

subroutine track1_bmad (orbit, ele, param, err_flag, track, mat6, make_matrix)

use bmad_interface, dummy4 => track1_bmad

implicit none

type (coord_struct) :: orbit    
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp), optional :: mat6(6,6)
real(rp)  knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer, parameter :: do_nothing$ = 9999
integer key, ix_pole_max

logical, optional :: err_flag
logical, optional :: make_matrix
logical err

character(*), parameter :: r_name = 'track1_bmad'

! If element is off... 

if (present(err_flag)) err_flag = .false.

key = ele%key

if (.not. ele%is_on) then
  select case (key)
  case (ab_multipole$, multipole$, taylor$, match$, fiducial$, floor_shift$)
    key = do_nothing$
  case (lcavity$, sbend$, patch$)
    ! Note: LCavities will still do wakefields.
  case default
    key = drift$
  end select
endif

! Select.

select case (key)

!-----------------------------------------------
! beambeam
                        
case (beambeam$)
  call track_a_beambeam (orbit, ele, param, track, mat6, make_matrix)

!-----------------------------------------------
! crab_cavity
                        
case (converter$)
  call track_a_converter (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! crab_cavity
                        
case (crab_cavity$)
  call track_a_crab_cavity(orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Thick multipoles

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$, ac_kicker$, kicker$, hkicker$, vkicker$) 
  call track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! drift
 
case (drift$) 
  call track_a_drift (orbit, ele%value(l$)*orbit%time_dir, mat6, make_matrix)

!-----------------------------------------------
! elseparator

case (elseparator$)
  call track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! GKicker

case (gkicker$)
  call track_a_gkicker(orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! LCavity: Linac rf cavity.

case (lcavity$)
  call track_a_lcavity (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! marker, etc.
! Note: floor_shift elements can have finite length in the case where it is a slice_slave of a taylor 
! element (the first slice is a taylor element and all other slices are floor_shifts).

case (marker$, fork$, photon_fork$, floor_shift$, fiducial$, detector$, beginning_ele$)

  orbit%t = orbit%t + ele%value(delta_ref_time$)

!-----------------------------------------------
! mask

case (mask$)
  call track_a_mask (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! match

case (match$)
  call track_a_match (orbit, ele, param, err_flag, mat6, make_matrix)

!-----------------------------------------------
! multipole, ab_multipole

case (multipole$, ab_multipole$) 

  if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, set$, orbit, set_tilt = .false.)

  call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilt)

  if (ix_pole_max > -1) then
    call multipole_kicks (knl, tilt, ele, orbit, ref_orb_offset = (ele%key == multipole$))

    if (logic_option(.false., make_matrix)) then
      call multipole_kick_mat (knl, tilt, param%particle, ele, orbit, 1.0_rp, ele%mat6)

      ! if knl(0) is non-zero then the reference orbit itself is bent
      ! and we need to account for this.

      if (knl(0) /= 0 .and. ele%key == multipole$) then
        ele%mat6(2,6) = knl(0) * cos(tilt(0))
        ele%mat6(4,6) = knl(0) * sin(tilt(0))
        ele%mat6(5,1) = -ele%mat6(2,6)
        ele%mat6(5,3) = -ele%mat6(4,6)
      endif
    endif
  endif

 if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, unset$, orbit, set_tilt = .false.)

!-----------------------------------------------
! octupole
! The octupole is modeled using kick-drift.

case (octupole$, thick_multipole$)
  call track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! patch

case (patch$)
  call track_a_patch (ele, orbit, mat6 = mat6, make_matrix = make_matrix)

!-----------------------------------------------
! pickup

case (pickup$)
  call track_a_pickup(orbit, ele, param, err_flag, mat6, make_matrix)

!-----------------------------------------------
! quadrupole

case (quadrupole$)
  call track_a_quadrupole (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! rfcavity

case (rfcavity$)
  call track_a_rfcavity (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sad_multipole

case (sad_mult$)
  call track_a_sad_mult (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sbend

case (sbend$)
  call track_a_bend (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sextupole
! The sextupole is modeled using kick-drift.

case (sextupole$)
  call track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Solenoid

case (sol_quad$, solenoid$)
  call track_a_sol_quad (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! foil

case (foil$)
  call track_a_foil (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Taylor

case (taylor$)
  call track1_taylor (orbit, ele, mat6 = mat6, make_matrix = make_matrix)

!-----------------------------------------------
! wiggler:

case (wiggler$, undulator$)
  call track_a_wiggler (orbit, ele, param, mat6, make_matrix)

!-----------------------------------------------
! do nothing case

case (do_nothing$)

!-----------------------------------------------
! Confused case

case (capillary$, crystal$, diffraction_plate$, lens$, mirror$, multilayer_mirror$, sample$, photon_init$)

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'ELEMENT: ' // ele%name, &
          'WHICH IS OF TYPE: ' // key_name(ele%key), &
          'CANNOT BE USED TO TRACK A NON-PHOTON PARTICLE.')
  if (global_com%exit_on_error) call err_exit
  return

!-----------------------------------------------
! unknown

case default

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'BMAD_STANDARD TRACKING_METHOD NOT IMPLMENTED FOR: ' // key_name(ele%key), &
          'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

end select

!-----------------------------------------------------------------------------------
! Set s-position

if (orbit%direction*orbit%time_dir == 1) then
  orbit%s = ele%s
else
  orbit%s = ele%s_start
endif

end subroutine track1_bmad
