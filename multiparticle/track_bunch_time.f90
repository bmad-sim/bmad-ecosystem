!+
! Subroutine track_bunch_time (lat, bunch, t_end, s_end, dt_step, extra_field)
!
! Routine to track a particle bunch for a given time step (or if the 
! particle position exceeds s_end).
!
! Input:
!   lat            -- lat_struct: Lattice to track through.
!   bunch          -- bunch_struct: Coordinates must be time-coords in element body frame.
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

subroutine track_bunch_time (lat, bunch, t_end, s_end, dt_step, extra_field)

use bmad_interface, dummy => track_bunch_time
!$ use omp_lib

implicit none

type (lat_struct), target :: lat
type (bunch_struct), target :: bunch
type (em_field_struct), optional :: extra_field(:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct), pointer :: orbit

real(rp) t_end, s_end
real(rp), optional :: dt_step(:)
real(rp) dt, significant_time

integer i, j
logical err

character(*), parameter :: r_name = 'track_bunch_time'

!

significant_time = bmad_com%significant_length / (10 * c_light)
branch => lat%branch(bunch%particle(1)%ix_branch)

!$OMP parallel do default(shared) private(dt, orbit, ele)

do i = 1, size(bunch%particle)
  if (present(dt_step)) then;  dt = dt_step(i)
  else;                        dt = bmad_com%init_ds_adaptive_tracking
  endif

  orbit => bunch%particle(i)

  do
    if (orbit%t >= t_end - significant_time) exit
    if (orbit%s >= s_end - 0.1_rp * bmad_com%significant_length) exit
    if (orbit%state /= alive$) exit

    ele => pointer_to_next_track_ele(orbit, branch)
    if (orbit%state /= alive$) exit

    if (present(extra_field)) then
      call track1_time_RK (orbit, ele, branch%param, err, t_end, dt, extra_field(i))
    else
      call track1_time_RK (orbit, ele, branch%param, err, t_end, dt)
    endif
  enddo

  if (present(dt_step)) dt_step(i) = dt
enddo

!$OMP end parallel do

!-------------------------------------------------------------------------
contains

function pointer_to_next_track_ele(orbit, branch) result (ele)

type (coord_struct) orbit
type (branch_struct) branch
type (ele_struct), pointer :: ele

!

ele => branch%ele(orbit%ix_ele)

select case (orbit%location)
case (upstream_end$)
  if (orbit%direction /= -1) return

  if (ele%ix_ele == 1) then
    orbit%state = lost_z_aperture$
    return
  endif

  ele => pointer_to_next_ele(ele, -1)
  orbit%location = downstream_end$
  orbit%ix_ele = ele%ix_ele
  orbit%vec(5) = ele%value(l$)

case (downstream_end$)
  if (orbit%direction /= 1) return

  if (ele%ix_ele == branch%n_ele_track) then
    call out_io (s_error$, r_name, 'PARTICLE AT END OF LATTICE. WILL BE MARKED AS DEAD')
    orbit%state = lost_z_aperture$
    return
  endif

  ele => pointer_to_next_ele(ele)
  orbit%location = upstream_end$
  orbit%ix_ele = ele%ix_ele
  orbit%vec(5) = 0
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
real(rp) t_end, dt_step, dt_ref, rf_time
logical err, set_spin

! If at edge of element 

if (orbit%location /= inside$) then
  if (orbit%direction == 1) then
    dt_ref = orbit%t - (ele%ref_time - ele%value(delta_ref_time$))
  else
    dt_ref = orbit%t - ele%ref_time
  endif

  call convert_particle_coordinates_t_to_s(orbit, dt_ref, ele, s_body_calc(orbit, ele))

  if (ele%key /= patch$ .and. ele%value(l$) == 0) then
    call track_a_zero_length_element (orbit, ele, param, orbit, err)
    call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation, dt_ref)
    return
  endif

  set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., set_spin = set_spin)
  call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation, dt_ref)
endif

! Track

dt_ref = 0  ! 
rf_time = particle_rf_time (orbit, ele, .true., orbit%s - ele%s_start)
call odeint_bmad_time(orbit, dt_ref, ele, param, +1, rf_time, err, &
                                        t_end = t_end, dt_step = dt_step, extra_field = extra_field)

! If at edge of element.

if (orbit%location /= inside$) then
  if (orbit%direction == 1) then
    dt_ref = orbit%t - (ele%ref_time - ele%value(delta_ref_time$))
  else
    dt_ref = orbit%t - ele%ref_time
  endif

  call convert_particle_coordinates_t_to_s(orbit, dt_ref, ele, s_body_calc(orbit, ele))
  set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., set_spin = set_spin)
  call convert_particle_coordinates_s_to_t(orbit, s_body_calc(orbit, ele), ele%orientation, dt_ref)
endif

end subroutine track1_time_RK

end subroutine track_bunch_time

