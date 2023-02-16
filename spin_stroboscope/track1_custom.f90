!+
! Subroutine track1_custom (orbit, ele, param, err_flag, finished, track)
!
! Dummy routine for custom tracking. 
! This routine needs to be replaced for a custom calculation.
! If not replaced and this routine is called, this routine will generate an error message.
!
! Also see:
!   track1_preprocess
!   track1_postprocess
!
! If this routine takes into account radiation damping and/or excitation when bmad_com%radiation_damping_on 
! and/or bmad_com%radiation_fluctuations_on is True, a custom version of track1_preprocess should be 
! constructed to set its radiation_included argument to True.
! If not, the track1 routine will use track1_radiation to include the radiation effects.
! Note: If this routine calles symp_lie_bmad, the symp_lie_bmad routine does take into account radiation effects.
!
! Note: If ele%spin_tracking_method = tracking$, then it is expected that this routine will also handle
! spin tracking. The alternative is when ele%spin_tracking_method = custom$ in which case track1_spin_custom will
! be called after this routine. If doing spin tracking here, bmad_com%spin_tracking_on should be checked
! to see if spin tracking is actually wanted.
!
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Input:
!   orbit      -- coord_struct: Starting position.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   orbit       -- coord_struct: End position.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_custom (orbit, ele, param, err_flag, finished, track)

use bmad, except_dummy => track1_custom

implicit none

type (coord_struct) :: orbit
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) cx, sx, cy, sy, cz, sz, dphi_x, dphi_y, dphi_z, angle
real(rp) phi_x, phi_y, phi_z
real(rp) mat_x(2,2), mat_y(2,2), mat_z(2,2), om(3)

real(rp) Qx, Qy, Qz
real(rp) Qx0, Qy0, Qz0
real(rp) eps_x, eps_y, eps_z
real(rp) nu_0

integer i, n_step

logical err_flag, finished

character(*), parameter :: r_name = 'track1_custom'

!

Qx     = get_this_attribute('QX')
Qy     = get_this_attribute('QY')
Qz     = get_this_attribute('QZ')
Qx0    = get_this_attribute('QX0')
Qy0    = get_this_attribute('QY0')
Qz0    = get_this_attribute('QZ0')
eps_x  = get_this_attribute('EPS_X')
eps_y  = get_this_attribute('EPS_Y')
eps_z  = get_this_attribute('EPS_Z')
nu_0   = get_this_attribute('NU_0')
n_step = nint(get_this_attribute('N_STEP'))

! 3D resonance model

err_flag = .false.
finished = .false.

dphi_x = Qx * twopi / n_step
cx = cos(dphi_x); sx = sin(dphi_x)
mat_x(1,:) = [ cx, sx]
mat_x(2,:) = [-sx, cx]

dphi_y = Qy * twopi / n_step
cy = cos(dphi_y); sy = sin(dphi_y)
mat_y(1,:) = [ cy, sy]
mat_y(2,:) = [-sy, cy]

dphi_z = Qz * twopi / n_step
cz = cos(dphi_z); sz = sin(dphi_z)
mat_z(1,:) = [ cz, sz]
mat_z(2,:) = [-sz, cz]

!

do i = 1, n_step

  ! Rotate orbital coords

  phi_x = -atan2(orbit%vec(2), orbit%vec(1)) + dphi_x / 2
  orbit%vec(1:2) = matmul(mat_x, orbit%vec(1:2))

  phi_y = -atan2(orbit%vec(4), orbit%vec(3)) + dphi_y / 2
  orbit%vec(3:4) = matmul(mat_y, orbit%vec(3:4))

  phi_z = -atan2(orbit%vec(6), orbit%vec(5)) + dphi_z / 2
  orbit%vec(5:6) = matmul(mat_z, orbit%vec(5:6))

  ! Spin transformation.

  if (orbit%vec(1) == 0 .and. orbit%vec(2) == 0) eps_x = 0
  if (orbit%vec(3) == 0 .and. orbit%vec(4) == 0) eps_y = 0
  if (orbit%vec(5) == 0 .and. orbit%vec(6) == 0) eps_z = 0

  om = [eps_x * cos(phi_x + twopi * Qx0) + &
        eps_y * cos(phi_y + twopi * Qy0) + &
        eps_z * cos(phi_z + twopi * Qz0), &
        eps_x * sin(phi_x + twopi * Qx0) + &
        eps_y * sin(phi_y + twopi * Qy0) + &
        eps_z * sin(phi_z + twopi * Qz0), &
        nu_0] * (twopi / n_step)
  angle = norm2(om)
  orbit%spin = rotate_vec_given_axis_angle (orbit%spin, om/angle, angle)

enddo

!-----------------------------------------------------
contains

function get_this_attribute(name) result (value)

type (all_pointer_struct) a_ptr
real(rp) value
logical err
character(*) name

!

call pointer_to_attribute (ele, name, .true., a_ptr, err)
value = a_ptr%r

end function

end subroutine
