!+
! Subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an crab_cavity element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: crab_cavity element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_crab_cavity (orbit, ele, param, mat6, make_matrix)

use fringe_mod, except_dummy => track_a_crab_cavity

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) voltage, phase0, phase, t0, z, length, charge_dir, dt_ref, beta_ref

integer i, n_slice, orientation

logical, optional :: make_matrix

!

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

length = ele%value(l$)
n_slice = max(1, nint(length / ele%value(ds_step$))) 
charge_dir = rel_tracking_charge_to_mass(orbit, param) * ele%orientation
voltage = e_accel_field(ele, voltage$) * charge_dir
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
dt_ref = length / (c_light * beta_ref)

phase0 = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
        (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))
if (ele%orientation == -1) phase0 = phase0 + twopi * ele%value(rf_frequency$) * dt_ref
phase = phase0

!call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

! Track through slices.

t0 = orbit%t

do i = 0, n_slice

  if (logic_option(.false., make_matrix)) then
  endif

  orbit%vec(2) = orbit%vec(2)

  if (i /= n_slice) then
    z = orbit%vec(5)
    call track_a_drift (orbit, length/n_slice, mat6, make_matrix)
    phase = phase + twopi * ele%value(rf_frequency$) * (orbit%vec(5) - z) / (c_light * orbit%beta)
  endif

enddo

! coupler kick, multipoles, back to lab coords.

!call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

end subroutine
