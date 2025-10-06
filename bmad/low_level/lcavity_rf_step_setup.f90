!+
! Subroutine lcavity_rf_step_setup(ele)
!
! Routine to construct the RF step parameters in ele%rf.
!
! Input:
!   ele                    -- ele_struct: Lcavity element.
!
! Output:
!   ele                    -- ele_struct: Element with ele%rf properly setup.
!-

recursive subroutine lcavity_rf_step_setup(ele)

use bmad_routine_interface, dummy => lcavity_rf_step_setup
use super_recipes_mod, only: super_zbrent, super_bracket_root

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
real(rp) field_scale, E_tot_final, x_range(2), t_final, delta_ref_time, dt_err
integer nn, status, i
logical err_flag

character(*), parameter :: r_name = 'lcavity_rf_step_setup'

!

nn = nint(ele%value(n_rf_steps$))
if (nn < 1) return
if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) return

!

if (.not. associated(ele%rf)) allocate(ele%rf)
if (allocated(ele%rf%steps)) then
  if (ubound(ele%rf%steps, 1) /= nn+1) deallocate(ele%rf%steps)
endif

if (.not. allocated(ele%rf%steps)) then
  allocate(ele%rf%steps(0:nn+1))
  do i = 0, nn+1
    ele%rf%steps(i)%ix_step = i
  enddo
endif

if (ele%slave_status /= multipass_slave$) then
  ele%value(phi0_autoscale$) = 0
  ele%value(field_autoscale$) = 1
endif

!

if (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  call this_rf_multipass_slave_setup(ele, lord)

elseif (ele%lord_status == multipass_lord$) then
  ! The idea is that we want the kick reference times to have mirror symmetry to handle ERLs 
  ! where particles are both accelerated and decellarated.

  ! Initial guess. If the reference energy gain is small compared to the voltage, accept
  ! the initial guess as a reasonable approximation.

  call this_rf_free_ele_setup(ele)
  if (abs(ele%value(E_tot$) - ele%value(E_tot_start$)) < 0.01*ele%value(voltage$)) return

  ! Converge to a solution

  x_range = super_bracket_root(this_dt_track, 0.99_rp*ele%value(delta_ref_time$), 1.01_rp*ele%value(delta_ref_time$), status)
  delta_ref_time = super_zbrent(this_dt_track, x_range(1), x_range(2), 1e-12_rp, 0.0_rp, status)
  dt_err = this_dt_track(delta_ref_time, status)
  !! print *, 'dE_rel_err: ', this_dE_track(ele%value(field_autoscale$), status)

else ! Not multipass
  call this_rf_free_ele_setup(ele)
endif

!----------------------------------------------------------------------------------------------------
contains

function this_dt_track(delta_ref_time, status) result (dt_err)

type (rf_stair_step_struct), pointer :: step, step2

real(rp), intent(in) :: delta_ref_time
real(rp) dt_err, x_range(2), field_scale, dt, dE_rel_err

integer status, i, j, nn

! Symetrize reference times

nn = nint(ele%value(n_rf_steps$))
ele%rf%steps(nn+1)%time = delta_ref_time
ele%value(delta_ref_time$) = delta_ref_time
ele%ref_time = ele%value(ref_time_start$) + delta_ref_time

do i = 0, (nn+1)/2
  j = nn - i
  step => ele%rf%steps(i)
  step2 => ele%rf%steps(j)
  dt = 0.5_rp * (delta_ref_time - step%time - step2%time)
  step%time  = step%time  + dt
  if (j /= i) step2%time = step2%time + dt
enddo

!

x_range = super_bracket_root(this_dE_track, 0.99_rp*ele%value(field_autoscale$), 1.01_rp*ele%value(field_autoscale$), status)
field_scale = super_zbrent(this_dE_track, x_range(1), x_range(2), 1e-12_rp, 0.0_rp, status)
dE_rel_err = this_dE_track(field_scale, status)
dt_err = t_final - delta_ref_time

end function this_dt_track

!----------------------------------------------------------------------------------------------------
! contains

function this_dE_track(field_scale, status) result (dE_rel_err)

type (rf_stair_step_struct), pointer :: step

real(rp), intent(in) :: field_scale
real(rp) dE_rel_err, dE_track, E_tot, pc, beta, t, phase, dE, dt_err, dt
integer i, nn, status

! Track through lcavity to compute the difference between "actual" energy gain and wanted (reference) energy gain.

ele%value(field_autoscale$) = field_scale
nn = nint(ele%value(n_rf_steps$))
E_tot = ele%value(E_tot_start$)
pc = ele%value(p0c_start$)
t = 0

do i = 0, nn+1
  step => ele%rf%steps(i)
  beta = pc / E_tot
  t = t + (step%s - step%s0) / (c_light * beta)
  dt = t - step%time
  phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(rf_frequency$) * dt)
  dE = step%scale * ele%value(voltage$) * ele%value(field_autoscale$) * cos(phase)
  E_tot = E_tot + dE
  call convert_total_energy_to(E_tot, ele%ref_species, pc = pc)
enddo

t_final = t
dE_rel_err = (E_tot - ele%value(E_tot$)) / ele%value(E_tot$)

end function this_dE_track

!----------------------------------------------------------------------------------------------------
! contains

subroutine this_rf_free_ele_setup(ele)

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch

real(rp) t, ds, beta, E_tot0, E_tot1, fac, scale, dE_ref, p0c, p1c, mass, phi, s0
integer i, nn

! Set e_tot$ and p0c$

phi = twopi * ele%value(phi0$)
if (.not. bmad_com%absolute_time_tracking) phi = phi + twopi * ele%value(phi0_multipass$)

!

nn = nint(ele%value(n_rf_steps$))
scale = 1.0_rp / nn
ds = ele%value(l_active$) * scale
ele%rf%ds_step = ds
mass = mass_of(ele%ref_species)

dE_ref = (ele%value(E_tot$) - ele%value(E_tot_start$)) * scale
E_tot0 = ele%value(E_tot_start$)
E_tot1 = E_tot0 + 0.5_rp * dE_ref
p0c = ele%value(p0c_start$)
s0 = 0.5_rp * (ele%value(l$) - ele%value(l_active$))
beta = p0c / E_tot0
t = s0 / (c_light * beta)
call convert_total_energy_to(E_tot1, ele%ref_species, pc = p1c)
ele%rf%steps(0) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, 0.5_rp * scale, t, 0.0_rp, s0, 0)

fac = 1.0_rp
do i = 1, nn
  E_tot0 = E_tot1
  if (i == nn) fac = 0.5_rp
  E_tot1 = E_tot0 + fac * dE_ref
  p0c = p1c
  call convert_total_energy_to(E_tot1, ele%ref_species, pc = p1c)
  beta = p0c / E_tot0
  t = t + ds / (c_light * beta)
  ele%rf%steps(i) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, fac*scale, t, (i-1)*ds+s0, i*ds+s0, i)
enddo

ele%rf%steps(nn)%E_tot1 = ele%value(E_tot$)
ele%rf%steps(nn)%p1c = ele%value(p0c$)
t = t + s0 / (c_light * ele%value(p0c$) / ele%value(E_tot$))
ele%rf%steps(nn+1) = rf_stair_step_struct(ele%value(E_tot$), ele%value(E_tot$), ele%value(p0c$), &
                                         ele%value(p0c$), 0.0_rp, t, ele%value(l$)-s0, ele%value(l$), nn+1)
ele%ref_time = ele%value(ref_time_start$) + t
ele%value(delta_ref_time$) = t
 
end subroutine this_rf_free_ele_setup

!----------------------------------------------------------------------------------------------------
! contains

subroutine this_rf_multipass_slave_setup(ele, lord)

type (ele_struct), target :: ele, lord
type (ele_struct), pointer :: slave1
type (rf_stair_step_struct), pointer :: steps(:)
type (rf_stair_step_struct), pointer :: step

real(rp) time_ref, p0c, phase, dE, beta
integer i, nn

! It can happen that if the slave tracking_method is switched to bmad_standard, the lord has
! not yet been setup.

if (.not. associated(lord%rf)) call lcavity_rf_step_setup(lord)

! Now need to correct for each rf%step: %E_tot0, %E_tot1, %p0c, %p1c. Also element delta_ref_time.
! Everything else is the same as the lord.

ele%rf = lord%rf
steps => ele%rf%steps
nn = nint(lord%value(n_rf_steps$))

step => steps(0)
step%E_tot0 = ele%value(E_tot_start$)
step%p0c = ele%value(p0c_start$)
beta = step%p0c / step%E_tot0    ! Ref beta in drift after this step
time_ref = 0

do i = 0, nn+1
  step => steps(i)

  if (i > 0) then
    step%E_tot0 = steps(i-1)%E_tot1
    step%p0c = steps(i-1)%p1c
  endif

  beta = step%p0c / step%E_tot0
  time_ref = time_ref + (step%s - step%s0) / (c_light * beta)

  phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + (time_ref - step%time) * ele%value(rf_frequency$))
  dE = step%scale * lord%value(voltage$) * lord%value(field_autoscale$) * cos(phase)

  if (i == nn) then   ! At the end must match desired ref energy
    steps(nn)%p1c      = ele%value(p0c$)
    steps(nn)%E_tot1   = ele%value(E_tot$)
  else
    step%E_tot1 = step%E_tot0 + dE
    call convert_total_energy_to(step%E_tot1, ele%ref_species, pc = step%p1c)
  endif
enddo

ele%value(delta_ref_time$) = time_ref
ele%ref_time = ele%value(ref_time_start$) + time_ref

end subroutine this_rf_multipass_slave_setup

end subroutine lcavity_rf_step_setup
