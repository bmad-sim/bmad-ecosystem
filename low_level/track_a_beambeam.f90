!+
! Subroutine track_a_beambeam (orbit, ele, param, track, mat6, make_matrix)
!
! Bmad_standard tracking through a beambeam element.
!
! If the track arg is present the track is storage order is:
!   0) starting orbit at element center
!   1) orbit before first slice
!   2) orbit after first slice
!   .... similarly for all other slices
!   N) ending orbit at element center    ! N = track%n_pt
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Beambeam element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit         -- coord_struct: End position.
!   track         -- track_struct, optional: Structure holding the track information if the 
!                      lattice element does tracking step-by-step. See track1 for more details.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_beambeam (orbit, ele, param, track, mat6, make_matrix)

use fringe_mod, except_dummy => track_a_beambeam

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track
type (em_field_struct) field
type (strong_beam_struct) sbb

real(rp), optional :: mat6(6,6)
real(rp) sig_x, sig_y, x_center, y_center, s_lab, s_body, ff
real(rp) s_slice, s_slice_old, k0_x, k0_y, k_xx1, k_xy1, k_yx1, k_yy1, k_xx2, k_xy2, k_yx2, k_yy2, coef, del
real(rp) mat21, mat23, mat41, mat43, del_s, x_pos, y_pos, ratio, bbi_const, z, dx, dy, dcoef
real(rp), allocatable :: z_slice(:)
real(rp) om(3), quat(0:3)

integer i, n_slice

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_beambeam'

!

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

del = 0.001
if (ele%value(sig_x$) == 0 .or. ele%value(sig_x$) == 0) then
  if (strong_beam_strength(ele) == 0) return
  call out_io (s_error$, r_name, 'STRONG BEAM SIGMAS NOT SET FOR BEAMBEAM ELEMENT: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

if (present(track)) call save_a_step(track, ele, param, .false., orbit, 0.0_rp, strong_beam = strong_beam_struct())

s_slice = 0    ! Begin at the IP
call offset_particle (ele, set$, orbit, s_pos = s_slice, s_out = s_slice, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)

n_slice = max(1, nint(ele%value(n_slice$)))
allocate(z_slice(n_slice))
call bbi_slice_calc (ele, n_slice, z_slice)

do i = 1, n_slice
  z = z_slice(i)      ! Positive z_slice is the tail of the strong beam.
  s_slice_old = s_slice   ! Where particle is (body frame).
  s_slice = (orbit%vec(5) + z) / 2  ! Where particle needs to drift to.

  call offset_particle(ele, unset$, orbit, s_pos = s_slice_old, s_out = s_lab, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
  del_s = s_slice - s_slice_old
  call solenoid_track_and_mat (ele, del_s, param, orbit, orbit, mat6, make_matrix)
  call offset_particle(ele, set$, orbit, s_pos = s_lab+del_s, s_out = s_body, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
  ! The drifting above will be slightly inaccurate if there is a x_pitch or y_pitch so try again. 
  ! Since drifting is linear, this should be fairly accurate. 
  ! There will still be a small inaccuracy due trajectories not being straight in a solenoid field.
  if (s_body /= s_slice) then
    call offset_particle(ele, unset$, orbit, s_pos = s_body, s_out = s_lab, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
    del_s = del_s * (s_slice - s_body) / (s_body - s_slice_old)
    call solenoid_track_and_mat (ele, del_s, param, orbit, orbit, mat6, make_matrix)
    call offset_particle(ele, set$, orbit, s_pos = s_lab+del_s, s_out = s_slice, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
  endif

  call strong_beam_sigma_calc (ele, s_slice, -z, sig_x, sig_y, bbi_const, x_center, y_center)

  dx = orbit%vec(1) - x_center
  dy = orbit%vec(3) - y_center
  x_pos = dx / sig_x
  y_pos = dy / sig_y

  if (present(track)) then
    sbb = strong_beam_struct(i, x_center, y_center, sig_x, sig_y, dx, dy)
    call save_a_step(track, ele, param, .true., orbit, s_slice, strong_beam = sbb)
  endif

  ratio = sig_y / sig_x
  call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)

  coef = bbi_const / n_slice
  orbit%vec(2) = orbit%vec(2) + k0_x * coef
  orbit%vec(4) = orbit%vec(4) + k0_y * coef

  if (present(track)) then
    call save_a_step(track, ele, param, .true., orbit, s_slice, strong_beam = sbb)
  endif

  if (logic_option(.false., make_matrix)) then
    call bbi_kick (x_pos-del, y_pos, ratio, k_xx1, k_yx1)
    call bbi_kick (x_pos, y_pos-del, ratio, k_xy1, k_yy1)
    call bbi_kick (x_pos+del, y_pos, ratio, k_xx2, k_yx2)
    call bbi_kick (x_pos, y_pos+del, ratio, k_xy2, k_yy2)

    dcoef = bbi_const / (ele%value(n_slice$) * del)
    mat21 = dcoef * (k_xx2 - k_xx1) / (2 * sig_x)
    mat23 = dcoef * (k_xy2 - k_xy1) / (2 * sig_y)
    mat41 = dcoef * (k_yx2 - k_yx1) / (2 * sig_x)
    mat43 = dcoef * (k_yy2 - k_yy1) / (2 * sig_y)

    mat6(2,:) = mat6(2,:) + mat21 * mat6(1,:) + mat23 * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat41 * mat6(1,:) + mat43 * mat6(3,:)
  endif

  if (bmad_com%spin_tracking_on) then
    ff = 1.0_rp + orbit%beta**2
    field%E = [ k0_x, k0_y, 0.0_rp] * (coef * orbit%p0c / (ff * charge_of(orbit%species)))
    field%B = [k0_y, -k0_x, 0.0_rp] * (coef * orbit%p0c * orbit%beta / (ff * c_light * charge_of(orbit%species)))
    om = spin_omega (field, orbit, +1)
    quat = omega_to_quat(om)
    orbit%spin = quat_rotate(quat, orbit%spin)
  endif
enddo

call offset_particle(ele, unset$, orbit, s_pos = s_slice, s_out = s_lab, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
call solenoid_track_and_mat (ele, -s_lab, param, orbit, orbit, mat6, make_matrix)

if (present(track)) call save_a_step(track, ele, param, .false., orbit, 0.0_rp, strong_beam = strong_beam_struct())

end subroutine

