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

subroutine lcavity_rf_step_setup(ele)

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

real(rp) t, ds, beta, E_tot0, E_tot1, fac, dp0c, scale, dE_ref, dE_amp, p0c, p1c, mass
integer i, n

!

n = nint(ele%value(n_rf_steps$))
t = 0.0_rp
scale = 1.0_rp / n
ds = ele%value(l$) * scale
ele%rf%ds_step = ds
mass = mass_of(ele%ref_species)

dE_ref = (ele%value(E_tot$) - ele%value(E_tot_start$)) * scale
dE_amp = (ele%value(voltage$) + ele%value(voltage_err$)) * scale
E_tot0 = ele%value(E_tot_start$)
E_tot1 = E_tot0 + 0.5_rp * dE_ref
p0c = ele%value(p0c_start$)
dp0c = dpc_given_dE(p0c, mass, 0.5_rp * dE_ref)
p1c = p0c + dp0c
ele%rf%steps(0) = rf_stair_step_struct(E_tot0, E_tot1, p0c, dp0c, 0.5_rp * dE_amp, 0.5_rp * scale, t, 0.0_rp)

fac = 1.0_rp
do i = 1, n
  E_tot0 = E_tot1
  if (i == n) fac = 0.5_rp
  E_tot1 = E_tot0 + fac * dE_ref
  p0c = p1c
  dp0c = dpc_given_dE(p0c, mass, fac * dE_ref)
  p1c = p0c + dp0c
  beta = p0c / E_tot0
  t = t + ds / (c_light * beta)
  ele%rf%steps(i) = rf_stair_step_struct(E_tot0, E_tot1, p0c, dp0c, fac*dE_amp, fac*scale, t, i * ds)
enddo

ele%rf%steps(n+1) = rf_stair_step_struct(ele%value(E_tot$), ele%value(E_tot$), ele%value(p0c$), &
                                                    0.0_rp, 0.0_rp, real_garbage$, real_garbage$, real_garbage$)
ele%ref_time = ele%value(ref_time_start$) + t
 
end subroutine this_rf_setup

!----------------------------------------------------------------------------------------------------
! contains

subroutine this_rf_multipass_slave_setup(ele, lord)

type (ele_struct), target :: ele
type (ele_struct) lord
type (rf_stair_step_struct), pointer :: steps(:)
type (rf_stair_step_struct), pointer :: step, step1

real(rp) t, p0c, mass, phase, dE, beta
integer i, n

! Correct for the fact that reference particle transit time through the slave will be different than the 
! reference transit time through the multipass_lord and this will give RF phase shifts and will shift
! the reference energy of the slave.

ele%rf = lord%rf
steps => ele%rf%steps
mass = mass_of(ele%ref_species)
n = nint(lord%value(n_rf_steps$))

! Now need to correct for each rf%step: %E_tot0, %E_tot1, %p0c, %dp0c.

t = 0.0_rp    ! Ref time of slave relative to the lord
steps(0)%E_tot0 = ele%value(E_tot_start$)
steps(0)%p0c = ele%value(p0c_start$)


do i = 0, n
  step => steps(i)
  phase = ele%value(rf_frequency$) * t + ele%value(phi0$)
  dE = step%scale * lord%value(voltage$) * cos(twopi * phase)
  step%E_tot1 = step%E_tot0 + dE
  step%dp0c = dpc_given_dE(step%p0c, mass, dE)

  ! Calc at next step kick point
  step1 => steps(i+1)
  beta = (step%p0c + step%dp0c) / step%E_tot1    ! Beta in drift after this step
  t = t + (step1%s - step%s) / (c_light * beta) - (step1%dtime - step%dtime)
  step1%E_tot0 = step%E_tot1
  step1%p0c = step%p0c + step%dp0c
enddo

step => steps(n+1)
step%E_tot1 = step%E_tot0

ele%value(E_tot$) = step%E_tot1
ele%value(p0c$)   = step%p0c

end subroutine this_rf_multipass_slave_setup

end subroutine lcavity_rf_step_setup
