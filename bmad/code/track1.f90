!+
! Subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation, make_map1, init_to_edge)
!
! Particle tracking through a single element. 
! This includes synchrotron radiation and space charge kicks if enabled by the appropriate
! bmad_com%radiation_damping_on, etc. components. See the code for more details. 
!
! Note: The ignore_radiation argument is meant as a "temporary" override to turn off
! all radiation and space-charge effects independent of the settings of the bmad_com%xxx
! switches. This is used by routines that are tracking the reference particle.
!
! Note: The init_to_edge argument is needed for historical reasons since calling routines may
! not call track1 with start_orb properly initialized to be at the element's edge. 
! Initializing start_orb to the element's edge is fine for the vast majority of cases.
! Setting init_to_edge = False will handle the other cases (EG: Dark current tracking).
! 
!
! Input:
!   start_orb     -- Coord_struct: Starting position.
!   ele           -- Ele_struct: Element to track through.
!   param         -- lat_param_struct: Reference particle info.
!   track         -- track_struct, optional: Structure holding existing track.
!   ignore_radiation
!                 -- Logical, optional: If present and True then do not include radiation
!                    effects along with space charge effects. 
!   make_map1     -- logical, optional: Make ele%mat6 and ele%spin_q components? Default is false.
!   init_to_edge  -- logical, optional: Default is True. If True then force the tracked particle to
!                    begin at the element's edge. See above. 
!                    Do not use this argument unless you know what you are doing.
!
! Output:
!   ele           -- ele_struct: Modified if make_map1 is True.
!   end_orb       -- coord_struct: End position.
!   track         -- track_struct, optional: Structure holding the track information if the 
!                      tracking method does tracking step-by-step.
!                      When tracking through multiple elements, the trajectory in an element
!                      is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!   err_flag      -- logical, optional: Set true if there is an error. False otherwise.
!                      Note: The particle getting lost (EG hitting an aperture) is *not* an error.
!                      An error is something like start_orb not being properly initialized.
!-

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation, make_map1, init_to_edge)

use bmad_routine_interface, except_dummy1 => track1
use mad_mod, only: track1_mad
use high_energy_space_charge_mod, only: track1_high_energy_space_charge
use radiation_mod, only: track1_radiation

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
type (ele_struct)   :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) p0c_start
integer n, tracking_method, stm

character(*), parameter :: r_name = 'track1'

logical, optional :: make_map1
logical, optional :: err_flag, ignore_radiation, init_to_edge
logical err, do_extra, finished, radiation_included, do_spin_tracking, time_RK_tracking, do_save

! Set end_orb%location appropriate for tracking through the element.
! For time runge-kutta the particle may be starting inside the element.
! In this case, do not set the location.

if (present(err_flag)) err_flag = .true.

start2_orb = start_orb    ! Don't want to modify start_orb

if (logic_option(.true., init_to_edge)) then
  if (start2_orb%direction*start2_orb%time_dir == 1) then
    start2_orb%location = upstream_end$
    start2_orb%s = ele%s_start
  else
    start2_orb%location = downstream_end$
    start2_orb%s = ele%s
  endif
endif

! Set ix_ele. If the element is a slice_slave then the appropriate ix_ele is given by the lord.

if (ele%slave_status == slice_slave$) then
  start2_orb%ix_ele    = ele%lord%ix_ele
  start2_orb%ix_branch = ele%lord%ix_branch
else
  start2_orb%ix_ele    = ele%ix_ele
  start2_orb%ix_branch = ele%ix_branch
endif

! For historical reasons, the calling routine may not have correctly 
! initialized the starting orbit. If so, we do an init here.
! Time Runge-Kutta is tricky so do not attempt to do a set.

time_RK_tracking = (ele%tracking_method == time_runge_kutta$ .or. ele%tracking_method == fixed_step_time_runge_kutta$)
if (start2_orb%state == not_set$) then
  if (time_RK_tracking) then
    call out_io (s_error$, r_name, 'STARTING ORBIT NOT PROPERLY INITIALIZED FOR TIME_RUNGE_KUTTA!', &
                                   'PLEASE REPORT THIS!')
    return
  endif
  call init_coord(start2_orb, start2_orb, ele, upstream_end$, particle = default_tracking_species(param)) 
endif

! Other tracking methods save their own trajectories

do_save = present(track)
select case (ele%tracking_method)
case (bmad_standard$)
  if (ele%key == beambeam$) do_save = .false.
case (linear$, taylor$, mad$)
case default
  do_save = .false.
end select

if (do_save) then
  if (start2_orb%direction*start2_orb%time_dir == 1) then
    call save_a_step(track, ele, param, .false., start2_orb, 0.0_rp)
  else
    call save_a_step(track, ele, param, .false., start2_orb, ele%value(l$))
  endif
endif

! "Preprocess" with custom code.
! Note: Energy ramping can be done here so track1_preprocess is called before the energy sanity checks below.

radiation_included = .false.

if (associated(track1_preprocess_ptr))then
  finished = .false.
  call track1_preprocess_ptr (start2_orb, ele, param, err, finished, radiation_included, track)
  if (finished .or. err .or. start2_orb%state /= alive$) then
    if (present(err_flag)) err_flag = err
    end_orb = start2_orb
    return
  endif
endif

! symp_lie_bmad tracking does include radiation effects

do_extra = .not. logic_option(.false., ignore_radiation)
radiation_included = (radiation_included .or. ele%tracking_method == symp_lie_bmad$) 

! Energy sanity checks.

if (start2_orb%species /= photon$ .and. start2_orb%state == alive$) then
  err = .false.

  if (ele%key == beginning_ele$) then
    if (significant_difference(ele%value(p0c$), start2_orb%p0c, rel_tol = small_rel_change$)) err = .true. ! For e_gun case

  else if (start2_orb%location == upstream_end$) then
    if (significant_difference(ele%value(p0c_start$), start2_orb%p0c, rel_tol = small_rel_change$)) err = .true.

  else if (start2_orb%location == downstream_end$) then
    if (significant_difference(ele%value(p0c$), start2_orb%p0c, rel_tol = small_rel_change$)) err = .true. ! For e_gun case
  endif

  if (err) then
    call out_io (s_error$, r_name, 'STARTING ORBIT NOT PROPERLY INITIALIZED! [NEEDS AN INIT_COORD CALL?]', &
                                   'TRACKING THROUGH ELEMENT: ' // ele_full_name(ele, '@N (&#)'))
    return
  endif
endif

! Photons get the z-coordinate reset to zero.

if (start2_orb%species == photon$ .and. start2_orb%location == downstream_end$) then
  start2_orb%vec(5) = 0
endif

! Check for particle lost even before we begin

if (start2_orb%state /= alive$) then
  end_orb = start2_orb
  start2_orb%vec = 0
  if (present(err_flag)) err_flag = .false.
  return
endif

! If a particle is inside the element then only time_runge_kutta
! can handle this situation.

if (start2_orb%state == inside$ .and. .not. time_RK_tracking) then
  call out_io (s_error$, r_name, 'PARTICLE''S STARTING POSITION IS INSIDE THE ELEMENT! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Init

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele)

! check for particles outside aperture.

call check_aperture_limit (start2_orb, ele, first_track_edge$, param)

! Radiation damping and/or fluctuations for the 1st half of the element.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) .and. &
                                           .not. radiation_included .and. ele%is_on .and. do_extra) then
  call track1_radiation (start2_orb, ele, start_edge$) 
endif

!

end_orb = start2_orb

if (end_orb%state /= alive$) then
  if (present(err_flag)) err_flag = .false.
  return
endif

! bmad_standard handles the case when the element is turned off.

do_spin_tracking = bmad_com%spin_tracking_on

tracking_method = ele%tracking_method
if (.not. ele%is_on) tracking_method = bmad_standard$

select case (tracking_method)

case (bmad_standard$)
  if (end_orb%species == photon$) then
    call track1_bmad_photon (end_orb, ele, param, err)
  else
    call track1_bmad (end_orb, ele, param, err, track, mat6 = ele%mat6, make_matrix = make_map1)

    select case (ele%key)
    case (beambeam$);  do_spin_tracking = .false.
    end select
  endif
  if (err) return

case (runge_kutta$, fixed_step_runge_kutta$) 
  call track1_runge_kutta (end_orb, ele, param, err, track, mat6 = ele%mat6, make_matrix = make_map1)
  if (err) return

case (linear$) 
  call track1_linear (end_orb, ele, param)

case (taylor$) 
  call track1_taylor (end_orb, ele)

case (symp_lie_bmad$) 
  call symp_lie_bmad (ele, param, end_orb, track, mat6 = ele%mat6, make_matrix = make_map1)

case (symp_lie_ptc$)
  call track1_symp_lie_ptc (end_orb, ele, param, track)
  stm = ele%spin_tracking_method
  if (stm == tracking$ .or. stm == symp_lie_ptc$) do_spin_tracking = .false.

case (mad$)
  call track1_mad (end_orb, ele, param)

case (custom$)
  if (.not. associated(track1_custom_ptr)) then
    call out_io (s_error$, r_name, 'TRACK1_CUSTOM_PTR HAS NOT BEEN SET IN THIS PROGRAM!', &
                                   'NEEDED FOR CUSTOM ELEMENT OR TRACKING_METHOD = CUSTOM: ' // ele%name)
    end_orb%state = lost$
    return
  endif
  call track1_custom_ptr (end_orb, ele, param, err, finished, track)
  if (err) return
  if (finished) then
    if (present(err_flag)) err_flag = err
    return
  endif

case (time_runge_kutta$, fixed_step_time_runge_kutta$)
  call track1_time_runge_kutta (end_orb, ele, param, err, track)
  if (err) return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING_METHOD: \i0\ ', ele%tracking_method)
  if (global_com%exit_on_error) call err_exit
  return
end select

!-----------------------------------------------------------------------------------
! Check

if (orbit_too_large (end_orb, param)) then
  if (present(err_flag)) err_flag = .false.
  return
endif

if (.not. time_RK_tracking) then
  if (end_orb%state /= alive$) then
    end_orb%location = inside$
  elseif (end_orb%direction*end_orb%time_dir == 1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
endif

if (end_orb%state /= alive$) then
  if (present(err_flag)) err_flag = .false.
  return
endif

! spin tracking. Must do after regular tracking in the case of spin_tracking_method = bmad_standard
 
if (do_spin_tracking) call track1_spin (start2_orb, ele, param, end_orb, make_map1)


! Radiation damping and/or fluctuations for the last half of the element

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) .and. &
                                            .not. radiation_included .and. ele%is_on .and. do_extra) then
  call track1_radiation (end_orb, ele, end_edge$)
endif

! space charge

if (bmad_com%high_energy_space_charge_on .and. do_extra) call track1_high_energy_space_charge (ele, param, end_orb)

! check for particles outside aperture

call check_aperture_limit (end_orb, ele, second_track_edge$, param)

if (associated(track1_postprocess_ptr)) call track1_postprocess_ptr (start_orb, ele, param, end_orb)

if (do_save) then
  if (start_orb%direction*start_orb%time_dir == 1) then
    call save_a_step(track, ele, param, .false., end_orb, ele%value(l$))
  else
    call save_a_step(track, ele, param, .false., end_orb, 0.0_rp)
  endif
endif

if (present(err_flag)) err_flag = .false.

end subroutine
