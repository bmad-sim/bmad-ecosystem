!+
! Subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag, mat6, make_matrix)
!
! Particle tracking through a single element BMAD_standard style.
!
! Input:
!   start_orb   -- Coord_struct: Starting position
!   ele         -- Ele_struct: Element
!   param       -- lat_param_struct:
!     %particle     -- Particle type
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element. Only used with bmad_standard tracking.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   end_orb     -- Coord_struct: End position
!   err_flag    -- Logical, optional: Set true if there is an error. False otherwise.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!-

subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag, mat6, make_matrix)

use track1_mod, dummy2 => track1_bmad
use geometry_mod, dummy4 => track1_bmad

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb, temp_orb
type (ele_struct) :: ele, temp_ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) k1, k2, k2l, k3l, length, hkick, vkick, kick
real(rp) coef, knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp) ks, kss, ksr, beta, xmat(6,6), mat4(4,4), vec0(6)
real(rp) rel_p, k_z, x_pos, y_pos, x, y, z, px, py, pz, k, dE0, L, E
real(rp) xp_start, yp_start, dz4_coef(4,4), kk, cos_a, sin_a, f, z_start, t_start
real(rp) dcos_phi, dpz, step_len, r_step, g, e_tot, pc, ps, charge_dir
real(rp) cosh1_k, sinh_k, dt, e_rel, p_factor, sin_g, angle, rel_tracking_charge

integer i, n, key, orientation, ix_sec, n_step, ix_pole_max
integer, parameter :: do_nothing$ = 9999

logical, optional :: err_flag
logical, optional :: make_matrix
logical err

character(*), parameter :: r_name = 'track1_bmad'

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

end_orb = start_orb     ! transfer start to end
if (end_orb%species /= photon$) then
  end_orb%p0c = ele%value(p0c$)
endif
length = ele%value(l$)
rel_p = 1 + start_orb%vec(6)
orientation = ele%orientation * start_orb%direction
rel_tracking_charge = rel_tracking_charge_to_mass(start_orb, param)
charge_dir = rel_tracking_charge * orientation

!-----------------------------------------------
! If element is off... 

key = ele%key
if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

if (.not. ele%is_on) then
  select case (key)
  case (ab_multipole$, multipole$, taylor$, match$, fiducial$, floor_shift$)
    key = do_nothing$
  case (lcavity$, sbend$, patch$)
    ! Note: LCavities will still do wakefields.
  case default
    key = drift$  
  end select
endif

! Select.

select case (key)

!-----------------------------------------------
! beambeam
                        
case (beambeam$)

  call track_a_beambeam(end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Thick multipoles

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$) 

  call track_a_thick_multipole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! drift
 
case (drift$) 

  call track_a_drift (end_orb, length, mat6, make_matrix)

!-----------------------------------------------
! elseparator

case (elseparator$)

  call track_a_elseparator (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! LCavity: Linac rf cavity.

case (lcavity$)

  call track_a_lcavity (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! marker, etc.
! Note: floor_shift elements can have finite length in the case where it is a slice_slave of a taylor 
! element (the first slice is a taylor element and all other slices are floor_shifts).

case (marker$, fork$, photon_fork$, floor_shift$, fiducial$, detector$)

  end_orb%t = end_orb%t + ele%value(delta_ref_time$)

!-----------------------------------------------
! mask

case (mask$)

  call track_a_mask (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! match

case (match$)

  call track_a_match (end_orb, ele, param, err_flag, mat6, make_matrix)

!-----------------------------------------------
! multipole, ab_multipole

case (multipole$, ab_multipole$) 

  call offset_particle (ele, param, set$, end_orb, set_multipoles = .false., set_tilt = .false.)

  call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilt)
  if (ix_pole_max > -1) call multipole_kicks (knl, tilt, param%particle, ele, end_orb, ref_orb_offset = (ele%key == multipole$))

  call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false., set_tilt = .false.)

!-----------------------------------------------
! octupole
! The octupole is modeled using kick-drift.

case (octupole$)

  n_step = max(1, nint(length / ele%value(ds_step$)))

  k3l = charge_dir * ele%value(k3$) * length / n_step

  call offset_particle (ele, param, set$, end_orb)

  end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 12
  end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 12

  do i = 1, n_step

    call track_a_drift (end_orb, length / n_step)

    if (i == n_step) then
      end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 12
      end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 12
    else
      end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 6
      end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 6
    endif

  enddo

  call offset_particle (ele, param, unset$, end_orb)

!-----------------------------------------------
! patch

case (patch$)

  call track_a_patch (ele, end_orb, mat6 = mat6, make_matrix = make_matrix)

!-----------------------------------------------
! quadrupole

case (quadrupole$)

  call track_a_quadrupole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! rfcavity

case (rfcavity$)

  call track_a_rfcavity (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sad_multipole

case (sad_mult$)

  call track_a_sad_mult (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sbend

case (sbend$)

  call track_a_bend (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sextupole
! The sextupole is modeled using kick-drift.

case (sextupole$)

  n_step = max(1, nint(length / ele%value(ds_step$)))

  call offset_particle (ele, param, set$, end_orb)

  do i = 0, n_step
    k2l = charge_dir * ele%value(k2$) * length / n_step
    if (i == 0 .or. i == n_step) k2l = k2l / 2
    end_orb%vec(2) = end_orb%vec(2) + k2l * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
    end_orb%vec(4) = end_orb%vec(4) + k2l * end_orb%vec(1) * end_orb%vec(3)
    if (i /= n_step) call track_a_drift (end_orb, length/n_step)
  enddo

  call offset_particle (ele, param, unset$, end_orb)

!-----------------------------------------------
! Solenoid
! Notice that ks is independent of the ele orientation

case (solenoid$)

  z_start = end_orb%vec(5)
  t_start = end_orb%t

  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.)
  call solenoid_track_and_mat (ele, param, end_orb, end_orb)
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.)

  if (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0) then
    ks = rel_tracking_charge * ele%value(ks$)
    ksr = ks / rel_p
    kss = ksr * length
    if (abs(kss) < 1d-2) then
      cos_a = length * (1 - kss**2 / 12 + kss**4 / 360) / 2
      sin_a = (1 - kss**2 / 12 + kss**4 / 240)
      f = length * kss * (1 - kss**2 / 20 + kss**4 / 840) / 6
    else
      cos_a = (1 - cos(kss)) / (kss * ksr)
      sin_a = (kss + sin(kss)) / (2 * kss)
      f = (kss - sin(kss)) / (kss * ksr)
    endif

    end_orb%vec(1:4) = end_orb%vec(1:4) + ele%value(hkick$) * [cos_a, sin_a, -f, -cos_a * ksr / 2] + &
                                          ele%value(vkick$) * [f, cos_a * ksr / 2, cos_a, sin_a]
  endif

  end_orb%t = t_start + ele%value(delta_ref_time$) + (z_start - end_orb%vec(5)) / (end_orb%beta * c_light)

!-----------------------------------------------
! sol_quad

case (sol_quad$)

  z_start = end_orb%vec(5)
  t_start = end_orb%t

  call offset_particle (ele, param, set$, end_orb)

  ks = rel_tracking_charge * ele%value(ks$)
  k1 = charge_dir * ele%value(k1$)
  vec0 = 0
  vec0(6) = end_orb%vec(6)
  call sol_quad_mat6_calc (ks, k1, length, vec0, xmat, dz4_coef)
  end_orb%vec(5) = end_orb%vec(5) + low_energy_z_correction (end_orb, ele, length) + &
                                      sum(end_orb%vec(1:4) * matmul(dz4_coef, end_orb%vec(1:4))) 

  end_orb%vec(1:4) = matmul (xmat(1:4,1:4), end_orb%vec(1:4))

  call offset_particle (ele, param, unset$, end_orb)

  end_orb%t = t_start + ele%value(delta_ref_time$) + (z_start - end_orb%vec(5)) / (end_orb%beta * c_light)

!-----------------------------------------------
! Taylor

case (taylor$)

  call track_a_taylor (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! wiggler:
! Only periodic type wigglers are handled here.
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

case (wiggler$, undulator$)

  call track_a_wiggler (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! do nothing case

case (do_nothing$)

!-----------------------------------------------
! unknown

case default

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'BMAD_STANDARD TRACKING_METHOD NOT IMPLMENTED FOR: ' // key_name(ele%key), &
          'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

end select

!-----------------------------------------------------------------------------------
! Set s-position

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine track1_bmad
