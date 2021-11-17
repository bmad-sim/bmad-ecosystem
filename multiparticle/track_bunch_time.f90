!+
! Subroutine track_bunch_time (lat, bunch, t_end, s_end, dt_step)
!
! Routine to track a particle bunch for a given time step (or if the 
! particle position exceeds s_end).
!
! Input:
!   lat         -- lat_struct: Lattice to track through.
!   bunch       -- bunch_struct: Bunch at end of element ix1.
!   t_end       -- real(rp): Ending time.
!   s_end       -- real(rp): Ending s-position.
!   dt_step(:)  -- real(rp), optional: Initial step to take for each particle. 
!                   Overrides bmad_com%init_ds_adaptive_tracking.
!
! Output:
!   bunch       -- bunch_struct: Bunch at end of element ix2.
!   dt_step(:)  -- real(rp), optional: Next RK time step that this tracker would take based on the error tolerance.
!-

subroutine track_bunch_time (lat, bunch, t_end, s_end, dt_step)

use bmad_interface, dummy => track_bunch_time

implicit none

type (lat_struct), target :: lat
type (bunch_struct), target :: bunch
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

    call track1_time_runge_kutta (orbit, ele, branch%param, orbit, err, t_end = t_end, dt_step = dt)
  enddo

  if (present(dt_step)) dt_step(i) = dt
enddo

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

case (downstream_end$)
  if (orbit%direction /= 1) return

  if (ele%ix_ele == branch%n_ele_track) then
    call out_io (s_error$, r_name, 'PARTICLE AT END OF LATTICE. WILL BE MARKED AS DEAD')
    orbit%state = lost_z_aperture$
    return
  endif

  ele => pointer_to_next_ele(ele)
  orbit%location = upstream_end$
end select

end function pointer_to_next_track_ele

end subroutine track_bunch_time

