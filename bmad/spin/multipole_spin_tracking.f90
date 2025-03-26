!+
! Subroutine multipole_spin_tracking (ele, param, orbit)
!
! Subroutine to track the spins in a multipole field.
!
! Input:
!   ele         -- Ele_struct: Element
!   param       -- Lat_param_struct
!   orbit       -- coord_struct: Particle coordinates.
!
! Output:
!   orbit       -- coord_struct: Particle coordinates.
!-

subroutine multipole_spin_tracking (ele, param, orbit)

use bmad_routine_interface, dummy => multipole_spin_tracking

implicit none

type (ele_struct) :: ele
type (lat_param_struct) param
type (coord_struct) orbit

complex(rp) kick, pos

real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), knl, Ex, Ey

integer n, sign_z_vel, ix_pole_max

! spin tracking for magnetic multipoles

call multipole_ele_to_ab(ele, .false., ix_pole_max, an, bn, magnetic$, include_kicks$)

! calculate kick_angle (for particle) and unit vector (Bx, By) parallel to B-field
! according to bmad manual, chapter "physics", section "Magnetic Fields"
! kick = qL/P_0*(B_y+i*Bx) = \sum_n (b_n+i*a_n)*(x+i*y)^n

if (ix_pole_max > -1) then
  kick = bn(0) + i_imag * an(0)
  pos = orbit%vec(1) + i_imag * orbit%vec(3)
  if (pos /= 0) then
    kick = kick + (bn(1) + i_imag * an(1)) * pos
    do n = 2, ix_pole_max
      pos = pos * (orbit%vec(1) + i_imag * orbit%vec(3))
      kick = kick + (bn(n) + i_imag * an(n)) * pos
    enddo
  endif

  ! Rotate spin

  sign_z_vel = orbit%direction * ele%orientation

  if (kick /= 0) then
    call rotate_spin_given_field (orbit, sign_z_vel, &
              [aimag(kick), real(kick), 0.0_rp] * (ele%value(p0c$) / (charge_of(param%particle) * c_light)))
  endif

  ! calculate rotation of local coordinate system due to dipole component

  if (ele%key == multipole$ .and. (bn(0) /= 0 .or. an(0) /= 0)) then
    kick = bn(0) + i_imag * an(0)
    call rotate_spin (orbit%time_dir*[-aimag(kick), -real(kick), 0.0_rp], orbit%spin)
  endif
endif

! Spin tracking for electric multipoles

call multipole_ele_to_ab(ele, .false., ix_pole_max, an, bn, electric$)
do n = 0, ix_pole_max
  if (an(n) == 0 .and. bn(n) == 0) cycle
  call elec_multipole_field(an(n), bn(n), n, orbit, Ex, Ey)
  call rotate_spin_given_field (orbit, sign_z_vel, EL = [Ex, Ey, 0.0_rp] * ele%value(l$))
enddo

end subroutine multipole_spin_tracking


