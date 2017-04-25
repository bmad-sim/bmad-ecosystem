!+
! Subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an rfcavity element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: rfcavity element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_rfcavity

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) beta_ref, dt_ref, voltage, phase, dE, z, charge_dir, length

integer i, n_step, orientation

logical, optional :: make_matrix
logical drifting

!

length = ele%value(l$)
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
n_step = max(1, nint(length / ele%value(ds_step$))) 
dt_ref = length / (c_light * beta_ref)
charge_dir = rel_tracking_charge_to_mass(orbit, param) * ele%orientation

call offset_particle (ele, param, set$, orbit)

! The cavity field is modeled as a standing wave antisymmetric wrt the center.
! Thus if the cavity is flipped (orientation = -1), the wave of interest, which is 
! always the accelerating wave, is the "backward" wave.

voltage = e_accel_field(ele, voltage$) * charge_dir

phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
        (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit)

! Track through slices.
! The phase of the accelerating wave traveling in the same direction as the particle is
! assumed to be traveling with a phase velocity the same speed as the reference velocity.

do i = 0, n_step

  dE = voltage * sin(phase) / n_step
  if (i == 0 .or. i == n_step) dE = dE / 2

  call apply_energy_kick (dE, orbit)
  
  if (orbit%vec(6) == -1) then
    orbit%state = lost_z_aperture$
    return
  endif

  if (i /= n_step) then
    z = orbit%vec(5)
    call track_a_drift (orbit, length/n_step)
    phase = phase + twopi * ele%value(rf_frequency$) * (orbit%vec(5) - z) / (c_light * orbit%beta)
  endif

enddo

! coupler kick

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit)

call offset_particle (ele, param, unset$, orbit)

!

end subroutine
