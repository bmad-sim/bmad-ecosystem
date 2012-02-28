!+
! Subroutine track1 (start_orb, ele, param, end_orb, track, err_flag)
!
! Particle tracking through a single element. 
! Optionally synchrotron radiation and space charge kicks can included.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb -- Coord_struct: Starting position.
!   ele       -- Ele_struct: Element to track through.
!   param     -- lat_param_struct:
!     %aperture_limit_on -- If True check if particle is lost by going outside
!                of the element aperture. 
!
! Output:
!   end_orb   -- Coord_struct: End position.
!   param
!     %lost          -- Set True If the particle cannot make it through an element.
!                         Set False otherwise.
!     %end_lost_at   -- entrance_end$ or exit_end$.
!     %plane_lost_at -- x_plane$, y_plane$ (for apertures), or 
!                         z_plane$ (turned around in an lcavity).
!   track     -- track_struct, optional: Structure holding the track information if the 
!                  tracking method does tracking step-by-step.
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!
! Notes:
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag)

use bmad, except_dummy1 => track1
use mad_mod, only: track1_mad
use boris_mod, only: track1_boris
use space_charge_mod, except_dummy2 => track1
use spin_mod, except_dummy3 => track1

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
type (ele_struct)   :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) p0c_start
integer tracking_method

character(8), parameter :: r_name = 'track1'

logical, optional :: err_flag
logical err

!

if (present(err_flag)) err_flag = .true.
start2_orb = start_orb

! Correct start_orb %beta and %p0c.
! Doing this here to be compatible with programs that do not set this.

if (ele%tracking_method /= time_runge_kutta$ .and. start_orb%beta < 0) then
  p0c_start = ele%value(p0c_start$)

  call convert_pc_to (p0c_start * (1 + start2_orb%vec(6)), param%particle, beta = start2_orb%beta)
  start2_orb%p0c = p0c_start
  start2_orb%status = outside$
endif

! custom

if (ele%tracking_method == custom$) then
  call track1_custom (start2_orb, ele, param, end_orb, err, track)
  if (present(err_flag)) err_flag = err
  return
endif

! Init

param%lost = .false.  ! assume everything will be OK

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele, param)

! check for particles outside aperture

if (ele%aperture_at == entrance_end$ .or. ele%aperture_at == both_ends$ .or. ele%aperture_at == continuous$) &
                              call check_aperture_limit (start2_orb, ele, entrance_end$, param)
if (start2_orb%status == dead$) then
  param%end_lost_at = entrance_end$
  end_orb = start2_orb
  if (present(err_flag)) err_flag = .false.
  return
endif

! Radiation damping and/or fluctuations for the 1st half of the element.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) .and. ele%is_on) then
  call track1_radiation (start2_orb, ele, param, start2_orb, start_edge$) 
endif

! bmad_standard handles the case when the element is turned off.

tracking_method = ele%tracking_method
if (.not. ele%is_on) tracking_method = bmad_standard$

select case (tracking_method)

case (bmad_standard$)
  call track1_bmad (start2_orb, ele, param, end_orb)

case (custom$)
  call track1_custom (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case (runge_kutta$) 
  call track1_runge_kutta (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case (linear$) 
  call track1_linear (start2_orb, ele, param, end_orb)

case (taylor$) 
  call track1_taylor (start2_orb, ele, param, end_orb)

case (symp_map$) 
  call track1_symp_map (start2_orb, ele, param, end_orb)

case (symp_lie_bmad$) 
  call symp_lie_bmad (ele, param, start2_orb, end_orb, .false., track)

case (symp_lie_ptc$) 
  call track1_symp_lie_ptc (start2_orb, ele, param, end_orb)

case (boris$) 
  call track1_boris (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case (mad$)
  call track1_mad (start2_orb, ele, param, end_orb)

case (custom2$)
  call track1_custom2 (start2_orb, ele, param, end_orb, err)
  if (err) return

case (time_runge_kutta$)
  call track1_time_runge_kutta (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING_METHOD: \i0\ ', ele%tracking_method)
  if (bmad_status%exit_on_error) call err_exit
  return

end select

! Radiation damping and/or fluctuations for the last half of the element

if ((bmad_com%radiation_damping_on .or. &
                  bmad_com%radiation_fluctuations_on) .and. ele%is_on) then
  call track1_radiation (end_orb, ele, param, end_orb, end_edge$) 
endif

! space charge

if (bmad_com%space_charge_on) &
      call track1_ultra_rel_space_charge (end_orb, ele, param, end_orb)

! spin tracking
 
if (bmad_com%spin_tracking_on) call track1_spin (start2_orb, ele, param, end_orb)

! check for particles outside aperture

if (.not. param%lost) then
  if (ele%aperture_at == exit_end$ .or. ele%aperture_at == both_ends$ .or. ele%aperture_at == continuous$) then
    call check_aperture_limit (end_orb, ele, exit_end$, param)
    if (param%lost) param%end_lost_at = exit_end$
  endif
endif

if (end_orb%p0c < 0 .and. end_orb%status /= dead$) then
  param%lost = .false. ! Temp
  if (ele%aperture_at == entrance_end$ .or. ele%aperture_at == both_ends$ .or. ele%aperture_at == continuous$) &
                  call check_aperture_limit (start2_orb, ele, entrance_end$, param)
  if (param%lost) then
    end_orb%status = dead$
    param%end_lost_at = entrance_end$
  endif
  param%lost = .true.
endif

if (present(err_flag)) err_flag = .false.

end subroutine
