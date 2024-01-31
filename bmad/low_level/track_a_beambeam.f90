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

use bmad_interface, except_dummy => track_a_beambeam
use super_recipes_mod, only: super_zbrent
implicit none

type (coord_struct), target :: orbit, orb_save
type (coord_struct), pointer :: orb_ptr
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track
type (em_field_struct) field
type (strong_beam_struct) sbb

real(rp), optional :: mat6(6,6)
real(rp) sigma(2), dsigma_ds(2), ff, ds_slice
real(rp) k0_x, k0_y, k_xx1, k_xy1, k_yx1, k_yy1, k_xx2, k_xy2, k_yx2, k_yy2, coef, del, s0
real(rp) xmat(6,6), del_s, x_rel, y_rel, ratio, bbi_const, dx, dy, dcoef, z, d
real(rp) om(3), quat(0:3), beta_strong, s_body_save, p_rel
real(rp) f_factor, new_beta, px_old, py_old, e_factor, dpx_dx, dpx_dy, dpy_dx, dpy_dy, ds_dz
real(rp), allocatable :: z_slice(:)
real(rp), target :: slice_center(3), s_body, s_lab
real(rp), pointer :: center_ptr(:), s_body_ptr, s_lab_ptr ! To get around ifort bug.

integer i, n_slice, status

logical, optional :: make_matrix
logical final_calc, make_mat

character(*), parameter :: r_name = 'track_a_beambeam'

!

center_ptr => slice_center
s_body_ptr => s_body
s_lab_ptr => s_lab
orb_ptr => orbit

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

del = 0.0001
if (ele%value(sig_x$) == 0 .or. ele%value(sig_x$) == 0) then
  if (strong_beam_strength(ele) == 0) return
  call out_io (s_error$, r_name, 'STRONG BEAM SIGMAS NOT SET FOR BEAMBEAM ELEMENT: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

s_lab = 0    ! Begin at the IP
if (present(track)) call save_a_step(track, ele, param, .false., orbit, s_lab, strong_beam = strong_beam_struct())

if (ele%value(species_strong$) /= real_garbage$ .and. ele%value(e_tot_strong$) > 0) then
  call convert_total_energy_to(ele%value(e_tot_strong$), nint(ele%value(species_strong$)), beta = beta_strong)
else
  beta_strong = ele%value(p0c$) / ele%value(E_tot$)
endif
s0 = 0.1_rp * ele%value(sig_z$) + &
        abs(particle_rf_time(orbit, ele, rf_freq = ele%value(repetition_frequency$)) * c_light * beta_strong)

call offset_particle (ele, set$, orbit, s_pos = s_lab, s_out = s_body, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)

n_slice = max(1, nint(ele%value(n_slice$)))
allocate(z_slice(n_slice))
call bbi_slice_calc (ele, n_slice, z_slice)
if (orbit%time_dir == -1) z_slice = z_slice(n_slice:1:-1)

do i = 1, n_slice
  z = z_slice(i)        ! Distance along strong beam axis. Positive z_slice is the tail of the strong beam.
  slice_center = strong_beam_center(ele, z) ! with respect to 

  final_calc = .false.; make_mat = .false.
  s_body_save = s_body
  orb_save = orbit
  s_lab = super_zbrent(at_slice_func, -abs(z)-s0, abs(z)+s0, 1e-12_rp, 1e-12_rp, status)

  final_calc = .true.; make_mat = logic_option(.false., make_matrix)
  s_body = s_body_save
  orbit = orb_save
  ds_slice = at_slice_func(s_lab, status)

  !

  call strong_beam_sigma_calc (ele, s_lab, sigma, bbi_const, dsigma_ds)

  dx = orbit%vec(1) - slice_center(1)
  dy = orbit%vec(3) - slice_center(2)
  x_rel = dx / sigma(1)
  y_rel = dy / sigma(2)
  px_old = orbit%vec(2)
  py_old = orbit%vec(4)
  p_rel = 1.0_rp + orbit%vec(6)

  if (present(track)) then
    sbb = strong_beam_struct(i, slice_center(1), slice_center(2), sigma(1), sigma(2), dx, dy)
    call save_a_step(track, ele, param, .true., orbit, s_body, strong_beam = sbb)
  endif

  ratio = sigma(2) / sigma(1)
  call bbi_kick (x_rel, y_rel, ratio, k0_x, k0_y)

  coef = orbit%time_dir * bbi_const / n_slice
  k0_x = k0_x * coef; k0_y = k0_y * coef

  call bbi_kick (x_rel-del, y_rel, ratio, k_xx1, k_yx1)
  call bbi_kick (x_rel, y_rel-del, ratio, k_xy1, k_yy1)
  call bbi_kick (x_rel+del, y_rel, ratio, k_xx2, k_yx2)
  call bbi_kick (x_rel, y_rel+del, ratio, k_xy2, k_yy2)

  dcoef = orbit%time_dir * bbi_const / (ele%value(n_slice$) * del)
  dpx_dx = dcoef * (k_xx2 - k_xx1) / (2 * sigma(1))
  dpx_dy = dcoef * (k_xy2 - k_xy1) / (2 * sigma(2))
  dpy_dx = dcoef * (k_yx2 - k_yx1) / (2 * sigma(1))
  dpy_dy = dcoef * (k_yy2 - k_yy1) / (2 * sigma(2))

  e_factor = 0.25_rp / p_rel
  orbit%vec(6) = orbit%vec(6) + e_factor * (k0_x * (k0_x + 2 * px_old) + k0_y * (k0_y + 2 * py_old)) + &
                        0.5_rp * (dpx_dx * dsigma_ds(1) * sigma(1) + dpy_dy * dsigma_ds(2) * sigma(2))
  call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, beta = new_beta)
  orbit%vec(5) = orbit%vec(5) * new_beta / orbit%beta
  orbit%beta = new_beta

  orbit%vec(2) = px_old + k0_x
  orbit%vec(4) = py_old + k0_y

  if (present(track)) then
    call save_a_step(track, ele, param, .true., orbit, s_body, strong_beam = sbb)
  endif

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit(xmat)
    ds_dz = 0.5_rp    ! Change in collision point with change in z

    xmat(2,1) = dpx_dx
    xmat(2,3) = dpx_dy
    xmat(4,1) = dpy_dx
    xmat(4,3) = dpy_dy

    xmat(1,5) = -ds_dz * k0_x / p_rel
    xmat(2,5) =  ds_dz * (px_old * dpx_dx + py_old * dpx_dy) / p_rel
    xmat(3,5) = -ds_dz * k0_y / p_rel
    xmat(4,5) =  ds_dz * (px_old * dpy_dx + py_old * dpy_dy) / p_rel
    xmat(6,5) =  e_factor * 2.0_rp * (xmat(2,5) * (k0_x + px_old) + xmat(4,5) * (k0_y + py_old))

    xmat(6,1) = e_factor * dpx_dx * 2.0_rp * (k0_x + px_old)
    xmat(6,2) = e_factor * 2.0_rp * k0_x
    xmat(6,3) = e_factor * dpy_dy * 2.0_rp * (k0_y + py_old)
    xmat(6,4) = e_factor * 2.0_rp * k0_y
    xmat(6,6) = 1.0_rp - e_factor * (k0_x * (k0_x + 2 * px_old) + k0_y * (k0_y + 2 * py_old)) / p_rel

    mat6 = matmul(xmat, mat6)
  endif

  if (bmad_com%spin_tracking_on) then
    ff = 1.0_rp + orbit%beta**2
    field%E = [k0_x,  k0_y, 0.0_rp] * (orbit%p0c / (ff * charge_of(orbit%species)))
    field%B = [k0_y, -k0_x, 0.0_rp] * (orbit%p0c * orbit%beta / (ff * c_light * charge_of(orbit%species)))
    om = spin_omega (field, orbit, +1)
    quat = omega_to_quat(om)
    orbit%spin = quat_rotate(quat, orbit%spin)
  endif
enddo

call offset_particle(ele, unset$, orbit, s_pos = s_body, s_out = s_lab, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
call solenoid_track_and_mat (ele, -s_lab, param, orbit, orbit, mat6, make_matrix)

if (present(track)) call save_a_step(track, ele, param, .false., orbit, 0.0_rp, strong_beam = strong_beam_struct())

!-------------------------------------------------------
contains

function at_slice_func(s_lab_target, status) result (ds_slice)

real(rp), intent(in) :: s_lab_target
real(rp) ds_slice
real(rp) s_body_target, s_lab_slice, del_s, s_lab, s_beam_center_strong, s_weak
integer status

!

call offset_particle(ele, unset$, orbit, s_pos = s_body, s_out = s_lab, set_spin = final_calc, mat6 = mat6, make_matrix = make_mat)
del_s = s_lab_target - s_lab
call solenoid_track_and_mat (ele, del_s, param, orbit, orbit, mat6, make_mat)
call offset_particle(ele, set$, orbit, s_pos = s_lab_target, s_out = s_body, set_spin = final_calc, mat6 = mat6, make_matrix = make_mat)

s_weak = -c_light * orbit%beta * particle_rf_time(orbit, ele, rf_freq = ele%value(repetition_frequency$)) 
s_beam_center_strong = (s_weak - ele%value(z_crossing$) - s_lab_target) * beta_strong / orbit%beta
s_body_target = 0.5_rp * (s_beam_center_strong + slice_center(3) + s_body)
ds_slice = s_body - s_body_target

end function at_slice_func

!-------------------------------------------------------
! contains

!+
! Function strong_beam_center(ele, z) result (center)
!
!   z_strong      -- real(rp): Position of slice within strong beam. Positive z_strong is at the tail of the bunch.
!   center(3)     -- real(rp): (x,y) position of slice in body coordinates.
!-

function strong_beam_center(ele, z) result (center)

type (ele_struct) ele
real(rp) r, z, z_strong, center(3)

!

z_strong = -z
r = (((ele%value(crab_x5$) * z_strong + ele%value(crab_x4$)) * z_strong + ele%value(crab_x3$)) * z_strong + & 
                                                              ele%value(crab_x2$)) * z_strong + ele%value(crab_x1$)

center(1:2) = z_strong * sin(r) * [cos(ele%value(crab_tilt$)), sin(ele%value(crab_tilt$))]
center(3) = -z_strong * cos(r)

end function strong_beam_center

end subroutine

