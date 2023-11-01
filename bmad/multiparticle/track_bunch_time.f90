!+
! Subroutine track_bunch_time (bunch, branch, t_end, s_end, dt_step, extra_field)
!
! Routine to track a particle bunch for a given time step (or if the particle position exceeds lattice end).
!
! Note: If a particle passes the end of the lattice, the particle will be tracked as if it was in a drift.
!
! Input:
!   bunch          -- bunch_struct: Coordinates must be time-coords in element body frame.
!   branch         -- branch_struct: Lattice branch being tracked through.
!   t_end          -- real(rp): Ending time.
!   s_end          -- real(rp): Ending s-position.
!   dt_step(:)     -- real(rp), optional: Initial step to take for each particle. 
!                      Overrides bmad_com%init_ds_adaptive_tracking.
!   extra_field(:) -- em_field_struct, optional: Per particle static field to be added to the lattice element field.
!                        Eg used with space charge.
!
! Output:
!   bunch          -- bunch_struct: Coordinates will be time-coords in element body frame.
!   dt_step(:)     -- real(rp), optional: Next RK time step that this tracker would take based on the error tolerance.
!-

subroutine track_bunch_time (bunch, branch, t_end, s_end, dt_step, extra_field)

use bmad_interface, dummy => track_bunch_time
!$ use omp_lib

implicit none

type (ele_struct), target :: drift_ele
type (ele_struct), pointer :: ele_here
type (bunch_struct), target :: bunch
type (branch_struct), target :: branch
type (em_field_struct), optional :: extra_field(:)
type (coord_struct), pointer :: orbit

real(rp) t_end, s_end
real(rp), optional :: dt_step(:)
real(rp) dt, significant_time

integer i, j, n
logical err

character(*), parameter :: r_name = 'track_bunch_time'

!

significant_time = bmad_com%significant_length / c_light

n = branch%n_ele_track
call init_ele (drift_ele, drift$, ix_ele = n + 1, branch = branch)
drift_ele%value(l$) = 100
drift_ele%name = 'end_drift'
call ele_compute_ref_energy_and_time(branch%ele(n), drift_ele, branch%param, err)

!$OMP parallel do default(shared) private(dt, orbit, ele_here, err)

do i = 1, size(bunch%particle)
  if (present(dt_step)) then;  dt = dt_step(i)
  else;                        dt = bmad_com%init_ds_adaptive_tracking
  endif

  orbit => bunch%particle(i)

  do
    if (orbit%t >= t_end - significant_time) exit
    if (orbit%s >= s_end - 0.1_rp * bmad_com%significant_length) exit
    if (orbit%state /= alive$) exit

    ele_here => pointer_to_next_track_ele(orbit, branch, drift_ele)
    if (orbit%state /= alive$) exit

    if (present(extra_field)) then
      call track1_time_RK (orbit, ele_here, branch%param, err, t_end, dt, extra_field(i))
    else
      call track1_time_RK (orbit, ele_here, branch%param, err, t_end, dt)
    endif
  enddo

  if (present(dt_step)) dt_step(i) = dt
enddo

!$OMP end parallel do

!-------------------------------------------------------------------------
contains

function pointer_to_next_track_ele(orbit, branch, drift_ele) result (ele_here)

type (coord_struct) orbit
type (branch_struct) branch
type (ele_struct), pointer :: ele_here
type (ele_struct), target :: drift_ele

integer ie
logical err_flag

!

ele_here => branch%ele(orbit%ix_ele)

select case (orbit%location)
case (upstream_end$)
  if (orbit%direction*orbit%time_dir /= -1) return

  ! If needed, setup a drift to take care of particles past the beginning of the branch.
  ! Except for branches that have e_guns and therefore have a cathode.
  if (ele_here%ix_ele == 1) then
    do ie = 1, branch%n_ele_track
      if (branch%ele(ie)%key == e_gun$) then
        orbit%state = lost_z$
        return
      endif
      if (branch%ele(ie)%key == drift$) exit
    enddo
    ele_here => drift_ele
    return
  endif

  ele_here => pointer_to_next_ele(ele_here, -1)
  orbit%location = downstream_end$
  orbit%ix_ele = ele_here%ix_ele
  orbit%vec(5) = ele_here%value(l$)

case (downstream_end$)
  if (orbit%direction*orbit%time_dir /= 1) return

  ! If needed, setup a drift to take care of particles past the end of the branch.
  if (ele_here%ix_ele == branch%n_ele_track) then
    ele_here => drift_ele
    drift_ele%ix_ele = -1
    return
  endif

  ele_here => pointer_to_next_ele(ele_here)
  orbit%location = upstream_end$
  orbit%ix_ele = ele_here%ix_ele
  orbit%vec(5) = 0

! Can happen that a particle that is before the beginning of the branch is inside$
case (inside$)
  if (orbit%s < branch%ele(0)%s) then
    ele_here => drift_ele
    drift_ele%ix_ele = -1
  endif
end select

end function pointer_to_next_track_ele

!-------------------------------------------------------------------------
! contains
! Similar to track1_time_runge_kutta except that input and output are time coords

subroutine track1_time_RK (orbit, ele, param, err, t_end, dt_step, extra_field)

use time_tracker_mod, only: odeint_bmad_time

type (coord_struct) orbit
type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct), optional :: extra_field
real(rp) t_end, dt_step, rf_time
logical err, set_spin

! If at edge of element 

if (orbit%location /= inside$) then
  call convert_particle_coordinates_t_to_s(orbit, ele, s_body_calc(orbit, ele))

  if (ele%key /= patch$ .and. ele%value(l$) == 0) then
    call track_a_zero_length_element (orbit, ele, param, err)
    call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation)
    return
  endif

  set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., set_spin = set_spin)
  call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation)
endif

! Track

rf_time = particle_rf_time (orbit, ele, .true., orbit%s - ele%s_start, time_coords = .true.)
call odeint_bmad_time(orbit, ele, param, +1, rf_time, err, t_end = t_end, dt_step = dt_step, extra_field = extra_field)

! If at edge of element.

if (orbit%location /= inside$) then
  call convert_particle_coordinates_t_to_s(orbit, ele, s_body_calc(orbit, ele))
  set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., set_spin = set_spin)
  call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation)
endif

end subroutine track1_time_RK

end subroutine track_bunch_time

