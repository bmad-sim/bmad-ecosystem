!+
! Subroutine track1_gun_space_charge (bunch, ele, err, to_s_coords)
!
! Subroutine to track a bunch of particles through an e_gun.
!
! Input:
!   bunch       -- bunch_struct: Starting bunch position.
!   ele         -- ele_struct: E_gun element to track through. Must be part of a lattice.
!   to_s_coords -- logical, optional: Default is True. If False lease bunch in time coords at end of tracking.
!
! Output:
!   bunch     -- bunch_struct: Ending bunch position.
!   err       -- logical: Set true if there is an error. EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch_space_charge (bunch, ele, err, to_s_coords)

use space_charge_mod, dummy => track1_bunch_space_charge

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (coord_struct), pointer :: p

integer i
logical err, finished
logical, optional :: to_s_coords
real(rp) :: dt_step, dt_next, t_now, t_end

integer, parameter :: fixed_time_step$ = 1, adaptive_step$ = 2 ! Need this in bmad_struct

character(*), parameter :: r_name = 'track1_bunch_space_charge'

! Initialize variables

t_now = 1e30_rp
branch => pointer_to_branch(ele)
dt_step = bmad_com%init_ds_adaptive_tracking / c_light  ! Init time step.
dt_next = dt_step
t_now = minval(bunch%particle(:)%t)

! Don't track zero-length elements
if (ele%value(l$) == 0) then
  if (logic_option(.true., to_s_coords)) call drift_to_s(bunch, ele%s, branch%lat, bunch)
  err = .false.
  return
endif

! Track
do
  if (ele%tracking_method==fixed_step_time_runge_kutta$) then
    call sc_step(bunch, ele, t_now+dt_step, dt_step)
  else
    call sc_adaptive_step(bunch, ele, t_now, dt_step, dt_next)
  end if

  t_now = t_now + dt_step
  dt_step = dt_next

  ! Check if all particles are finished
  finished = .true.
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    if (p%state == pre_born$ .or. (p%s < ele%s - 0.1_rp * bmad_com%significant_length .and. p%state == alive$)) then
      finished = .false.
      exit
    endif
  enddo
  if (finished) exit
enddo

if (logic_option(.true., to_s_coords)) call drift_to_s(bunch, ele%s, branch%lat, bunch)

err = .false.

end subroutine track1_bunch_space_charge

!+
! Drift a bunch of particles to the same s coordinate
!
! Input:
!   bunch_in  -- bunch_struct: input bunch position in t-based coordinate
!   s         -- real(rp): target s coordinate
!   lat       -- lat_struct: lattice particles is tracking through
!
! Output:
!   bunch_out -- bunch_struct: output bunch position in t-based coordinate. Particles will be at the same s coordinate
!-

subroutine drift_to_s (bunch_in, s, lat, bunch_out)
  
  use bmad
  
  implicit none
  
  type (bunch_struct), target :: bunch_in, bunch_out
  type (coord_struct), pointer :: p
  type (lat_struct) :: lat

  integer i
  real(rp) s, pz0, E_tot, dt

  bunch_out = bunch_in

  ! Convert bunch to s-based coordinates
  do i = 1, size(bunch_out%particle) 
    call convert_particle_coordinates_t_to_s(bunch_out%particle(i), 0.0_rp, lat%ele(bunch_out%particle(i)%ix_ele))
  enddo

  do i = 1, size(bunch_out%particle)
    p => bunch_out%particle(i)
    pz0 = sqrt( (1.0_rp + p%vec(6))**2 - p%vec(2)**2 - p%vec(4)**2 ) ! * p0 
    E_tot = sqrt((1.0_rp + p%vec(6))**2 + (mass_of(p%species)/p%p0c)**2) ! * p0
    dt = (s-p%s)/(c_light*pz0/E_tot)
  
    p%vec(1) = p%vec(1) + dt*c_light*p%vec(2)/E_tot
    p%vec(3) = p%vec(3) + dt*c_light*p%vec(4)/E_tot
    p%s = s
    p%t = p%t + dt
  enddo

end subroutine drift_to_s

