!+
! Subroutine lcavity_rf_step_setup(ele)
!
! Routine to construct the RF step parameters in ele%rf.
!
! Input:
!   ele         -- ele_struct: Lcavity element.
!
! Output:
!   ele         -- ele_struct: Element with ele%rf properly setup.
!-

recursive subroutine lcavity_rf_step_setup(ele)

use bmad_routine_interface, dummy => lcavity_rf_step_setup

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer n

! Note: The calling routine guarantees that ele is not a slice_slave nor a super_slave.

n = nint(ele%value(n_rf_steps$))
if (.not. associated(ele%rf)) allocate(ele%rf)
if (allocated(ele%rf%steps)) then
  if (ubound(ele%rf%steps, 1) /= n+1) deallocate(ele%rf%steps)
endif
if (.not. allocated(ele%rf%steps)) allocate(ele%rf%steps(0:n+1))

!

if (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  call this_rf_multipass_slave_setup(ele, lord)
else
  call this_rf_setup(ele)
endif

!----------------------------------------------------------------------------------------------------
contains

subroutine this_rf_setup(ele)

type (ele_struct) ele

real(rp) t, ds, beta, E_tot0, E_tot1, fac, scale, dE_ref, dE_amp, p0c, p1c, mass
integer i, n

!

n = nint(ele%value(n_rf_steps$))
t = 0.0_rp
scale = 1.0_rp / n
ds = ele%value(l$) * scale
ele%rf%ds_step = ds
mass = mass_of(ele%ref_species)

dE_ref = (ele%value(E_tot$) - ele%value(E_tot_start$)) * scale
dE_amp = (ele%value(voltage$) + ele%value(voltage_err$)) * scale ! Amplitude of RF field
E_tot0 = ele%value(E_tot_start$)
E_tot1 = E_tot0 + 0.5_rp * dE_ref
p0c = ele%value(p0c_start$)
call convert_total_energy_to(E_tot1, ele%ref_species, pc = p1c)
ele%rf%steps(0) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, 0.5_rp * dE_amp, 0.5_rp * scale, t, 0.0_rp, 0.0_rp, 0)

fac = 1.0_rp
do i = 1, n
  E_tot0 = E_tot1
  if (i == n) fac = 0.5_rp
  E_tot1 = E_tot0 + fac * dE_ref
  p0c = p1c
  call convert_total_energy_to(E_tot1, ele%ref_species, pc = p1c)
  beta = p0c / E_tot0
  t = t + ds / (c_light * beta)
  ele%rf%steps(i) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, fac*dE_amp, fac*scale, t, 0.0_rp, i * ds, i)
enddo

ele%rf%steps(n)%E_tot1 = ele%value(E_tot$)
ele%rf%steps(n)%p1c = ele%value(p0c$)
ele%rf%steps(n+1) = rf_stair_step_struct(ele%value(E_tot$), ele%value(E_tot$), ele%value(p0c$), &
                                         ele%value(p0c$), 0.0_rp, 0.0_rp, t, 0.0_rp, n * ds, n+1)
ele%ref_time = ele%value(ref_time_start$) + t
 
end subroutine this_rf_setup

!----------------------------------------------------------------------------------------------------
! contains

subroutine this_rf_multipass_slave_setup(ele, lord)

type (ele_struct), target :: ele, lord
type (ele_struct), pointer :: slave1
type (rf_stair_step_struct), pointer :: steps(:)
type (rf_stair_step_struct), pointer :: step, step1

real(rp) dt_rf, time_ref, dt_ref, p0c, mass, phase, dE, beta, p1c, dt
integer i, n

! It can happen that if the slave tracking_method is switched to bmad_standard, the lord has
! not yet been setup.

if (.not. associated(lord%rf)) call lcavity_rf_step_setup(lord)

! Correct for the fact that reference particle transit time through the slave will be different than the 
! reference transit time through the multipass_lord and this will give RF phase shifts and will shift
! the reference energy of the slave.

if (bmad_com%absolute_time_tracking) then
  slave1 => pointer_to_slave(lord, 1)
  dt = ele%value(ref_time_start$) - slave1%value(ref_time_start$)
endif

ele%rf = lord%rf
steps => ele%rf%steps
mass = mass_of(ele%ref_species)
n = nint(lord%value(n_rf_steps$))

! Now need to correct for each rf%step: %E_tot0, %E_tot1, %p0c, %p1c, %time, %dt_rf
! Everything else is the same as the lord.

time_ref = 0.0_rp
dt_rf = 0.0_rp    ! ref time of slave relative to the lord
steps(0)%E_tot0 = ele%value(E_tot_start$)
steps(0)%p0c = ele%value(p0c_start$)

do i = 0, n
  step => steps(i)
  phase = ele%value(rf_frequency$) * dt_rf + ele%value(phi0$)
  if (bmad_com%absolute_time_tracking) then
    phase = phase + dt * ele%value(rf_frequency$)
  else
    phase = phase + ele%value(phi0_multipass$)
  endif

  dE = step%scale * lord%value(voltage$) * cos(twopi * phase)
  step%E_tot1 = step%E_tot0 + dE
  call convert_total_energy_to(step%E_tot1, ele%ref_species, pc = p1c)
  step%p1c = p1c
  step%dt_rf = dt_rf
  step%time = time_ref

  ! Calc at next step kick point
  step1 => steps(i+1)
  step1%E_tot0 = step%E_tot1
  step1%p0c = p1c
  if (i /= n) then
    beta = step1%p0c / step1%E_tot0    ! Ref beta in drift after this step
    dt_ref = (step1%s - step%s) / (c_light * beta)
    time_ref = time_ref + dt_ref
    dt_rf = time_ref - step1%time
  endif
enddo

step => steps(n+1)
step%E_tot1 = step%E_tot0
step%p1c = step%p0c
step%time = steps(n)%time

ele%value(E_tot$) = step%E_tot1
ele%value(p0c$)   = step%p1c
ele%ref_time = ele%value(ref_time_start$) + time_ref

end subroutine this_rf_multipass_slave_setup

end subroutine lcavity_rf_step_setup
