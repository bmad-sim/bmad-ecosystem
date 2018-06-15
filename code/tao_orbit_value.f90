!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_orbit_value (component, orbit, value, err)
!
! Routine to return the orbit component indicated by component
!
! Input:
!   component -- character(*): 'orbit.x', 'orbit.px', 'intensity.x', 'phase.y', 'energy', 'pc', etc.
!   orbit     -- coord_struct: Particle orbit.
!
! Output:
!   value -- real(rp): orbit component.
!   err   -- logical: Set True if component is not recognized. False otherwise.
!-

subroutine tao_orbit_value (component, orbit, value, err)

use bmad_routine_interface

implicit none

type (coord_struct) orbit

real(rp) value
character(*) component
logical err

!

err = .true.

if (component == 'state') then
  value = orbit%state
  err = .false.
  return
endif

if (orbit%state /= alive$) return

select case (component)
case ('orbit_x', 'orbit.x')
  value = orbit%vec(1)
case ('orbit_px', 'orbit.px')
  value = orbit%vec(2)
case ('orbit_y', 'orbit.y')
  value = orbit%vec(3)
case ('orbit_py', 'orbit.py')
  value = orbit%vec(4)
case ('orbit_z', 'orbit.z')
  value = orbit%vec(5)
case ('orbit_pz', 'orbit.pz')
  value = orbit%vec(6)
case ('spin.x', 'spin_x')
  value = orbit%spin(1)
case ('spin.y', 'spin_y')
  value = orbit%spin(2)
case ('spin.z', 'spin_z')
  value = orbit%spin(3)
case ('intensity')
  value = orbit%field(1)**2 + orbit%field(2)**2
case ('intensity_x', 'intensity.x')
  value = orbit%field(1)**2
case ('intensity_y', 'intensity.y')
  value = orbit%field(2)**2
case ('phase_x', 'phase.x')
  value = orbit%phase(1)
case ('phase_y', 'phase.y')
  value = orbit%phase(2)
case ('t', 'time')
  value = orbit%t
case ('beta')
  value = orbit%beta
case ('energy')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, e_tot = value)
  endif
case ('pc')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    value = orbit%p0c * (1 + orbit%vec(6))
  endif
case default
  return
end select

err = .false.

end subroutine tao_orbit_value
