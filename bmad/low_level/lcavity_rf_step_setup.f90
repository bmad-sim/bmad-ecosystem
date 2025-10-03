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

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
integer nn
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
if (.not. allocated(ele%rf%steps)) allocate(ele%rf%steps(0:nn+1))

!

if (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  call this_rf_multipass_slave_setup(ele, lord)
else
  call this_rf_free_ele_setup(ele)
endif

!----------------------------------------------------------------------------------------------------
contains

subroutine this_rf_free_ele_setup(ele)

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch

real(rp) t, ds, beta, E_tot0, E_tot1, fac, scale, dE_ref, dE_amp, p0c, p1c, mass, phi
integer i, n

! Set e_tot$ and p0c$

phi = twopi * ele%value(phi0$)
if (.not. bmad_com%absolute_time_tracking) phi = phi + twopi * ele%value(phi0_multipass$)

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
ele%rf%steps(0) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, 0.5_rp * dE_amp, 0.5_rp * scale, t, 0.0_rp, 0)

fac = 1.0_rp
do i = 1, n
  E_tot0 = E_tot1
  if (i == n) fac = 0.5_rp
  E_tot1 = E_tot0 + fac * dE_ref
  p0c = p1c
  call convert_total_energy_to(E_tot1, ele%ref_species, pc = p1c)
  beta = p0c / E_tot0
  t = t + ds / (c_light * beta)
  ele%rf%steps(i) = rf_stair_step_struct(E_tot0, E_tot1, p0c, p1c, fac*dE_amp, fac*scale, t, i * ds, i)
enddo

ele%rf%steps(n)%E_tot1 = ele%value(E_tot$)
ele%rf%steps(n)%p1c = ele%value(p0c$)
ele%rf%steps(n+1) = rf_stair_step_struct(ele%value(E_tot$), ele%value(E_tot$), ele%value(p0c$), &
                                         ele%value(p0c$), 0.0_rp, 0.0_rp, t, n * ds, n+1)
ele%ref_time = ele%value(ref_time_start$) + t
 
end subroutine this_rf_free_ele_setup

!----------------------------------------------------------------------------------------------------
! contains

subroutine this_rf_multipass_slave_setup(ele, lord)

type (ele_struct), target :: ele, lord
type (ele_struct), pointer :: slave1
type (rf_stair_step_struct), pointer :: steps(:), lord_steps(:)
type (rf_stair_step_struct), pointer :: step, step1

real(rp) time_ref, p0c, mass, phase, dE, beta, p1c, s0
integer i, nn

! It can happen that if the slave tracking_method is switched to bmad_standard, the lord has
! not yet been setup.

if (.not. associated(lord%rf)) call lcavity_rf_step_setup(lord)

! Correct for the fact that reference particle transit time through the slave will be different than the 
! reference transit time through the multipass_lord and this will give RF phase shifts which will shift
! the reference energy of the slave.

time_ref = 0
s0 = 0.5_rp * (ele%value(l$) - ele%value(l_active$))
time_ref = time_ref + (ele%value(p0c_start$) / ele%value(E_tot_start$)) * c_light * s0

ele%rf = lord%rf
lord_steps => lord%rf%steps
steps => ele%rf%steps
mass = mass_of(ele%ref_species)
nn = nint(lord%value(n_rf_steps$))

! Now need to correct for each rf%step: %E_tot0, %E_tot1, %p0c, %p1c, %time
! Everything else is the same as the lord.

beta = ele%value(p0c_start$) / ele%value(E_tot_start$)
time_ref = s0 / (beta * c_light)
steps(0)%E_tot0 = ele%value(E_tot_start$)
steps(0)%p0c = ele%value(p0c_start$)

do i = 0, nn
  step => steps(i)

  phase = ele%value(phi0$) + ele%value(phi0_multipass_ref$) + ele%value(rf_frequency$) * (time_ref - lord_steps(i)%time)
  if (.not. bmad_com%absolute_time_tracking) phase = phase + ele%value(phi0_multipass$)

  phase = modulo2(twopi * phase, pi)
  dE = step%scale * lord%value(voltage$) * cos(phase)

  step%E_tot1 = step%E_tot0 + dE
  call convert_total_energy_to(step%E_tot1, ele%ref_species, pc = p1c)
  step%p1c = p1c
  step%time = time_ref

  ! Calc at next step kick point
  step1 => steps(i+1)
  step1%E_tot0 = step%E_tot1
  step1%p0c = p1c
  beta = step1%p0c / step1%E_tot0    ! Ref beta in drift after this step
  time_ref = time_ref + (step1%s - step%s) / (c_light * beta)
enddo

step => steps(nn+1)
step%E_tot1 = step%E_tot0
step%p1c = step%p0c
step%time = time_ref

ele%ref_time = ele%value(ref_time_start$) + time_ref

end subroutine this_rf_multipass_slave_setup

end subroutine lcavity_rf_step_setup
