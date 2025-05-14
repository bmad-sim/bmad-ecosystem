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
real(rp) xmat(6,6), del_s, bbi_const, dx, dy, dcoef, z, d, coef, ds
real(rp) om(3), quat(0:3), beta_strong, s_body_save, p_rel, nk(2), dnk(2,2)
real(rp) f_factor, new_beta, px_old, py_old, e_factor, ds_dz, s0, s00, s0_factor
real(rp), allocatable :: z_slice(:)
real(rp), target :: slice_center(3), s_body, s_lab, part_time0, part_time1
real(rp), pointer :: center_ptr(:), s_body_ptr, s_lab_ptr, part_time0_ptr, part_time1_ptr ! To get around ifort bug.

integer i, n_slice, status

logical, optional :: make_matrix
logical final_calc, make_mat

character(*), parameter :: r_name = 'track_a_beambeam'

!

center_ptr => slice_center
s_body_ptr => s_body
s_lab_ptr => s_lab
orb_ptr => orbit
part_time0_ptr => part_time0
part_time1_ptr => part_time1

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
  if (strong_beam_strength(ele) == 0) return
  call out_io (s_error$, r_name, 'STRONG BEAM SIGMAS NOT SET FOR BEAMBEAM ELEMENT: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

s_lab = 0    ! Begin at the IP

if (ele%value(species_strong$) /= real_garbage$ .and. ele%value(e_tot_strong$) > 0) then
  call convert_total_energy_to(ele%value(e_tot_strong$), nint(ele%value(species_strong$)), beta = beta_strong)
else
  beta_strong = ele%value(p0c$) / ele%value(E_tot$)
endif

part_time0 = particle_rf_time(orbit, ele, rf_freq = ele%value(repetition_frequency$))
part_time1 = particle_rf_time(orbit, ele, rf_freq = ele%value(repetition_frequency$), abs_time = .true.)

s0_factor = orbit%beta / (orbit%beta + beta_strong)
s00 = (ele%value(z_offset_tot$) + (ele%value(crossing_time$) - part_time0) * (c_light * orbit%beta)) * s0_factor

call offset_particle (ele, set$, orbit, s_pos = s_lab, s_out = s_body, set_spin = .true., mat6 = mat6, make_matrix = make_matrix)
if (present(track)) call save_a_step(track, ele, param, .true., orbit, s_body, strong_beam = strong_beam_struct())

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
  s0 = s00 + slice_center(3) * s0_factor
  ds = 1.0_rp
  s_lab = super_zbrent(at_slice_func, s0-ds, s0+ds, 1e-12_rp, 1e-12_rp, status)

  final_calc = .true.; make_mat = logic_option(.false., make_matrix)
  s_body = s_body_save
  orbit = orb_save
  ds_slice = at_slice_func(s_lab, status)

  !

  call strong_beam_sigma_calc (ele, s_lab, sigma, bbi_const, dsigma_ds)

  dx = orbit%vec(1) - slice_center(1)
  dy = orbit%vec(3) - slice_center(2)
  px_old = orbit%vec(2)
  py_old = orbit%vec(4)
  p_rel = 1.0_rp + orbit%vec(6)

  if (present(track)) then
    sbb = strong_beam_struct(i, slice_center(1), slice_center(2), sigma(1), sigma(2), dx, dy)
    call save_a_step(track, ele, param, .true., orbit, s_body, strong_beam = sbb)
  endif

  call bbi_kick (dx, dy, sigma, nk, dnk)

  coef = orbit%time_dir * bbi_const / n_slice
  nk = nk * coef

  dcoef = orbit%time_dir * bbi_const / n_slice
  dnk = dnk * dcoef

  e_factor = 0.25_rp / p_rel
  orbit%vec(6) = orbit%vec(6) + e_factor * (nk(1) * (nk(1) + 2 * px_old) + nk(2) * (nk(2) + 2 * py_old)) + &
                        0.5_rp * (dnk(1,1) * dsigma_ds(1) * sigma(1) + dnk(2,2) * dsigma_ds(2) * sigma(2))
  call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, beta = new_beta)
  orbit%vec(5) = orbit%vec(5) * new_beta / orbit%beta
  orbit%beta = new_beta

  orbit%vec(2) = px_old + nk(1)
  orbit%vec(4) = py_old + nk(2)

  if (present(track)) then
    call save_a_step(track, ele, param, .true., orbit, s_body, strong_beam = sbb)
  endif

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit(xmat)
    ds_dz = 0.5_rp    ! Change in collision point with change in z

    xmat(2,1) = dnk(1,1)
    xmat(2,3) = dnk(1,2)
    xmat(4,1) = dnk(2,1)
    xmat(4,3) = dnk(2,2)

    xmat(1,5) = -ds_dz * nk(1) / p_rel
    xmat(2,5) =  ds_dz * (px_old * dnk(1,1) + py_old * dnk(1,2)) / p_rel
    xmat(3,5) = -ds_dz * nk(2) / p_rel
    xmat(4,5) =  ds_dz * (px_old * dnk(2,1) + py_old * dnk(2,2)) / p_rel
    xmat(6,5) =  e_factor * 2.0_rp * (xmat(2,5) * (nk(1) + px_old) + xmat(4,5) * (nk(2) + py_old))

    xmat(6,1) = e_factor * dnk(1,1) * 2.0_rp * (nk(1) + px_old)
    xmat(6,2) = e_factor * 2.0_rp * nk(1)
    xmat(6,3) = e_factor * dnk(2,2) * 2.0_rp * (nk(2) + py_old)
    xmat(6,4) = e_factor * 2.0_rp * nk(2)
    xmat(6,6) = 1.0_rp - e_factor * (nk(1) * (nk(1) + 2 * px_old) + nk(2) * (nk(2) + 2 * py_old)) / p_rel

    mat6 = matmul(xmat, mat6)
  endif

  if (bmad_com%spin_tracking_on) then
    ff = 1.0_rp + orbit%beta**2
    field%E = [nk(1),  nk(2), 0.0_rp] * (orbit%p0c / (ff * charge_of(orbit%species)))
    field%B = [nk(2), -nk(1), 0.0_rp] * (orbit%p0c * orbit%beta / (ff * c_light * charge_of(orbit%species)))
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

!+
! function at_slice_func(s_lab_target, status) result (ds)
!
! Routine to propagate the weak particle to s-position s_lab_target and return the difference between the particle and the strong slice 
! 
! Input:
!   s_lab_target    -- real(rp): Particle s-position in lab frame
! 
! Output:
!   ds              -- real(rp): difference in s-position between particle and slice in the body frame.
!-

function at_slice_func(s_lab_target, status) result (ds)

real(rp), intent(in) :: s_lab_target
real(rp) ds
real(rp) del_s, s_lab, s_slice, dpart_time
integer status

! Calc particle s-position...

! Convert from body to lab coords
call offset_particle(ele, unset$, orbit, s_pos = s_body, s_out = s_lab, set_spin = final_calc, mat6 = mat6, make_matrix = make_mat)
! del_s = distance to propagate in lab coords. Note: Solenoid field is always defined in lab coords.
del_s = s_lab_target - s_lab
! Propagate in lab coords
call solenoid_track_and_mat (ele, del_s, param, orbit, orbit, mat6, make_mat)
! Convert back from lab coords to body coords
call offset_particle(ele, set$, orbit, s_pos = s_lab_target, s_out = s_body, set_spin = final_calc, mat6 = mat6, make_matrix = make_mat)

! Calc slice s-position...
! part_time is the time the weak particle is at s_body.

dpart_time = particle_rf_time(orbit, ele, rf_freq = ele%value(repetition_frequency$), abs_time = .true.) - part_time1
s_slice = slice_center(3) + (ele%value(crossing_time$) - (part_time0 + dpart_time)) * beta_strong * c_light

! Difference between particle and slice s-positions

ds = s_body - s_slice

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

