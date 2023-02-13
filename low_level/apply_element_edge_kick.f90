!+
! Subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix, rf_time, apply_sol_fringe)
!
! Subroutine, used with runge_kutta tracking, to track through the edge fringe field of an element.
! This routine is used and with the bmad_standard field_calc where the field can have an abrubt, 
! unphysical termination of the field at the edges of the element.
!
! This routine is used in conjunction with calc_next_fringe_edge. See, for example, odeint_bmad.
!
! Elements that have kicks due to unphysical edge field termination include:
!   sbend
!   solenoid
!   sol_quad
!   lcavity
!   rfcavity 
!   e_gun
! Additionally, Any element that has an electric multipole has an edge kick.
!
! Input:
!   orb               -- Coord_struct: Starting coords in element reference frame.
!   fringe_info       -- fringe_field_info_struct: Fringe information.
!   track_ele         -- ele_struct: Element being tracked through. Is different from fringe_info%hard_ele
!                          when there are superpositions and track_ele can be a super_slave of fringe_info%hard_ele.
!   param             -- lat_param_struct: lattice parameters.
!   track_spin        -- logical: Track the spin?
!   mat6(6,6)         -- Real(rp), optional: Transfer matrix before fringe.
!   make_matrix       -- logical, optional: Propagate the transfer matrix? Default is false.
!   rf_time           -- real(rp), optional: RF clock time. If not present then the time will be calculated using the standard algorithm.
!   apply_sol_fringe  -- logical, optional: Apply the solenoid fringe kick? Default is True.
!
! Output:
!   fringe_info   -- fringe_field_info_struct: Fringe information.
!   orb           -- Coord_struct: Coords after application of the edge fringe field.
!   mat6(6,6)     -- Real(rp), optional: Transfer matrix transfer matrix including fringe.
!-

subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix, rf_time, apply_sol_fringe)

use fringe_mod, except_dummy => apply_element_edge_kick

implicit none

type (ele_struct), target :: track_ele
type (coord_struct) orb
type (lat_param_struct) param
type (em_field_struct) field
type (fringe_field_info_struct) fringe_info
type (ele_struct), pointer :: hard_ele, lord

real(rp), optional :: mat6(6,6), rf_time
real(rp) ff, fac, l_drift, s_edge, s, phi, omega(3), pc, z_saved, beta_ref, ds
real(rp) a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)
complex(rp) xiy, c_vec

integer physical_end, dir, i, fringe_at, at_sign, sign_z_vel, particle_at, ix_elec_max
integer hard_ele_field_calc

logical, optional :: make_matrix, apply_sol_fringe
logical finished, track_spin, track_spn, err_flag

character(*), parameter :: r_name = 'apply_element_edge_kick'

! The setting of fringe_info%hard_location is used by calc_next_fringe_edge to calculate the next fringe location.

particle_at = fringe_info%particle_at
physical_end = physical_ele_end (particle_at, orb, track_ele%orientation)
hard_ele => track_ele

if (associated(fringe_info%hard_ele)) then
  if (associated(fringe_info%hard_ele)) hard_ele => fringe_info%hard_ele
  s_edge = fringe_info%s_edge_hard

  if (particle_at == first_track_edge$) then
    fringe_info%hard_location = inside$
  elseif (orb%direction * orb%time_dir * track_ele%orientation == 1) then
    fringe_info%hard_location = exit_end$
  else
    fringe_info%hard_location = entrance_end$
  endif

else
  if (physical_end == entrance_end$) then
    s_edge = 0
  else
    s_edge = track_ele%value(l$)
  endif
endif

hard_ele_field_calc = hard_ele%field_calc
if (hard_ele_field_calc == refer_to_lords$) then
  lord => pointer_to_lord(hard_ele, 1)
  hard_ele_field_calc = lord%field_calc
endif

! Custom edge kick?

call apply_element_edge_kick_hook (orb, fringe_info, track_ele, param, finished, mat6, make_matrix, rf_time)
if (finished) return

! With a solenoid must always apply the fringe kick due to the longitudinal field. 
! If not done the matrix calc will not be symplectic.
! For other elements, especially quadrupoles, this is problematic due to the soft edge kick not being being exactly the reverse
! going from inside to outside and vice versa (it is confusing if a superimposed marker shifts the tracking).

fringe_at = nint(track_ele%value(fringe_at$))
if (hard_ele%key /= solenoid$ .and. hard_ele%key /= sol_quad$ .and. hard_ele%key /= sad_mult$) then
  if (.not. at_this_ele_end(physical_end, fringe_at)) return
endif
track_spn = (track_spin .and. bmad_com%spin_tracking_on .and. is_true(hard_ele%value(spin_fringe_on$)))

if (particle_at == first_track_edge$) then
  at_sign = 1
else
  at_sign = -1
endif

sign_z_vel = orb%direction * track_ele%orientation

if (orb%beta == 0) then
  call out_io(s_error$, r_name, 'FRINGE OF ELEMENT: ' // hard_ele%name, &
                                'AT SPOT WHERE VELOCITY IS ZERO. WILL MARK AS DEAD.', &
                                'POSSIBLE SOLUTION: SET "FRINGE_AT = EXIT_END" FOR ELEMENT.')
  orb%state = lost$
endif

! Edge field when %field_calc = fieldmap$.

if (hard_ele_field_calc /= bmad_standard$ .and. hard_ele%tracking_method /= bmad_standard$) then
  call em_field_calc(hard_ele, param, s_edge, orb, .true., field, .false., err_flag, .true.)
  ff = at_sign * charge_of(orb%species) 
  fac = ff * c_light / orb%p0c

  if (at_sign == 1) then
    call apply_energy_kick (-ff * field%phi, orb, ff * field%E(1:2), mat6, make_matrix)
    if (track_spn) call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, -orb%time_dir*ff * field%phi])
    orb%vec(2) = orb%vec(2) - fac * field%A(1)
    orb%vec(4) = orb%vec(4) - fac * field%A(2)
  else
    orb%vec(2) = orb%vec(2) - fac * field%A(1)
    orb%vec(4) = orb%vec(4) - fac * field%A(2)
    if (track_spn) call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, -orb%time_dir*ff * field%phi])
    call apply_energy_kick (-ff * field%phi, orb, ff * field%E(1:2), mat6, make_matrix)
  endif

  return
endif

!------------------------------------------------------------------------------------
! Only need this section when the field_calc or tracking is bmad_standard.
! Note: track_ele, if a slave, will have track_ele%field_calc = refer_to_lords$.

if (hard_ele%key == e_gun$ .and. physical_end == entrance_end$) return ! E_gun does not have an entrance kick

! Static electric longitudinal field
! Note: magnetic fringe effects, which may shift the particle's (x, y) position, need to
! be applied before the electric fringe on entrance and after the electric fringe on exit in
! order to conserve the particle's energy.

call multipole_ele_to_ab (hard_ele, .false., ix_elec_max, a_pole_elec, b_pole_elec, electric$)
if (particle_at == second_track_edge$) call electric_longitudinal_fringe(orb, hard_ele, &
                                              a_pole_elec, b_pole_elec, ix_elec_max, sign_z_vel, at_sign, hard_ele_field_calc)

! Static magnetic and electromagnetic fringes

select case (hard_ele%key)
case (quadrupole$)
  if (particle_at == first_track_edge$) then
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
  else
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
  endif

case (sextupole$)
  if (particle_at == first_track_edge$) then
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
  else
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
  endif

case (sbend$)
  call bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix, track_spn)

! Sad_mult edge fields are 

case (sad_mult$)
  if (hard_ele%value(l$) == 0) return
  if (logic_option(.true., apply_sol_fringe)) call apply_this_sol_fringe(orb, hard_ele, at_sign, sign_z_vel, track_spn)

  if (particle_at == first_track_edge$) then
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    if (orbit_too_large (orb, param)) return
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    call sad_mult_hard_bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    if (orb%state /= alive$) return
    call sad_soft_bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
  else
    call sad_soft_bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    call sad_mult_hard_bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    if (orb%state /= alive$) return
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix)
    if (orbit_too_large (orb, param)) return
  endif

case (solenoid$, sol_quad$)
  if (logic_option(.true., apply_sol_fringe)) call apply_this_sol_fringe(orb, hard_ele, at_sign, sign_z_vel, track_spn)

case (lcavity$, rfcavity$, e_gun$)

  ! Add on bmad_com%significant_length to make sure we are just inside the cavity.

  if (at_this_ele_end(physical_end, nint(hard_ele%value(fringe_at$))) .and. hard_ele%value(l$) /= 0) then
    s = s_edge
    z_saved = orb%vec(5)
    beta_ref = hard_ele%value(p0c$) / hard_ele%value(e_tot$)
    ds = track_ele%s_start - hard_ele%s_start
    orb%vec(5) = orb%vec(5) - c_light * orb%beta * &
          ((track_ele%value(ref_time_start$) - hard_ele%value(ref_time_start$)) - ds / (beta_ref * c_light))
    if (physical_end == entrance_end$) then
      s = s + bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, orb, .true., field, rf_time = rf_time)
    else
      s = s - bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, orb, .true., field, rf_time = rf_time)
    endif
    orb%vec(5) = z_saved

    ff = at_sign * charge_of(orb%species) / (2 * orb%p0c)
    orb%vec(2) = orb%vec(2) - field%e(3) * orb%vec(1) * ff + c_light * field%b(3) * orb%vec(3) * ff
    orb%vec(4) = orb%vec(4) - field%e(3) * orb%vec(3) * ff - c_light * field%b(3) * orb%vec(1) * ff

    if (track_spn) then
    select case (hard_ele%key)
      case (lcavity$, rfcavity$)
        ff = orb%time_dir * at_sign * charge_of(orb%species) / 4.0_rp  ! Notice factor of 4 here
        call rotate_spin_given_field (orb, sign_z_vel, [-orb%vec(3), orb%vec(1), 0.0_rp] * (ff * field%e(3) / c_light), &
                                                       -[orb%vec(1), orb%vec(3), 0.0_rp] * (ff * field%e(3)))
      case default
        ff = orb%time_dir * at_sign * charge_of(orb%species) / 2.0_rp
        call rotate_spin_given_field (orb, sign_z_vel, -[orb%vec(1), orb%vec(3), 0.0_rp] * (ff * field%b(3)), &
                                                       -[orb%vec(1), orb%vec(3), 0.0_rp] * (ff * field%e(3)))
      end select
    endif

    ! orb%phase(1) is set by em_field_calc.

    call rf_coupler_kick (hard_ele, param, particle_at, orb%phase(1), orb)
  endif

case (elseparator$)
  ! Longitudinal fringe field
  if (hard_ele%value(l$) /= 0) then
    ff = at_sign * charge_of(orb%species) * (hard_ele%value(p0c$) / hard_ele%value(l$))
    phi = ff * (hard_ele%value(hkick$) * orb%vec(1) + hard_ele%value(vkick$) * orb%vec(3))
    call apply_energy_kick (phi, orb, [ff * hard_ele%value(hkick$), ff * hard_ele%value(vkick$)], mat6, make_matrix)
    if (track_spn) then
      call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, orb%time_dir*phi])
    endif
  endif
end select

! Entrance Static electric longitudinal field

if (particle_at == first_track_edge$) call electric_longitudinal_fringe(orb, hard_ele, &
                                                   a_pole_elec, b_pole_elec, ix_elec_max, sign_z_vel, at_sign, hard_ele_field_calc)

!--------------------------------------------------------------------------------
contains

subroutine electric_longitudinal_fringe(orb, hard_ele, a_pole_elec, b_pole_elec, ix_elec_max, sign_z_vel, at_sign, hard_ele_field_calc)

type (ele_struct) hard_ele
type (coord_struct) orb
type (em_field_struct) field
type (cartesian_map_struct), pointer :: ct
type (cylindrical_map_struct), pointer :: cy

real(rp) a_pole_elec(0:), b_pole_elec(0:), ff, E_r(2)
complex(rp) ab_elec, xiy_old
integer sign_z_vel, at_sign, hard_ele_field_calc
integer i, ix_elec_max
logical err_flag

! Multipole fringe

if (hard_ele%value(l$) == 0) return  ! Can get divide by zero problems with zero length elements.

if (hard_ele_field_calc == bmad_standard$) then
  if (ix_elec_max > -1) then
    ff = at_sign * charge_of(orb%species) 

    if (hard_ele%key == sbend$ .and. nint(hard_ele%value(exact_multipoles$)) /= off$ .and. hard_ele%value(g$) /= 0) then
      call bend_exact_multipole_field (hard_ele, param, orb, .true., field, .false., .true.)
      call apply_energy_kick (-ff * field%phi, orb, [ff * field%E(1), ff * field%E(2)], mat6, make_matrix)
      if (track_spn) call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, -orb%time_dir * ff * field%phi])

    else
      xiy = 1
      c_vec = cmplx(orb%vec(1), orb%vec(3), rp)
      phi = 0
      E_r = 0
      do i = 0, ix_elec_max
        xiy_old = xiy
        xiy = xiy * c_vec
        if (a_pole_elec(i) == 0 .and. b_pole_elec(i) == 0) cycle
        ab_elec = cmplx(b_pole_elec(i), -a_pole_elec(i), rp)
        phi = phi - real(ab_elec * xiy) / (i + 1)
        E_r = E_r + [real(ab_elec * xiy_old), -imag(ab_elec * xiy_old)]
      enddo
      call apply_energy_kick (-ff * phi, orb, ff * E_r, mat6, make_matrix)
      if (track_spn) call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, -orb%time_dir * ff * phi])
    endif
  endif
endif

end subroutine electric_longitudinal_fringe

!--------------------------------------------------------------------------------
! contains

subroutine apply_this_sol_fringe(orb, hard_ele, at_sign, sign_z_vel, track_spn)

type (coord_struct) orb
type (ele_struct) hard_ele

real(rp) ks4, ff, xy_orb(2)
integer at_sign, sign_z_vel
logical track_spn

! To make reverse tracking the same as forward tracking, use a symmetrical orbital-spin-orbital kick scheme.
! Note: Cannot trust hard_ele%value(ks$) here since element may be superimposed with an lcavity with changing
! ref energy. So use hard_ele%value(bs_field$).

ks4 = at_sign * charge_of(orb%species) * hard_ele%value(bs_field$) * c_light / (4.0_rp * orb%p0c)
xy_orb = [orb%vec(1), orb%vec(3)]
if (hard_ele%key == sad_mult$) then
  xy_orb = xy_orb + rot_2d ([hard_ele%value(x_offset_mult$), hard_ele%value(y_offset_mult$)], -hard_ele%value(tilt$))
endif

orb%vec(2) = orb%vec(2) + ks4 * xy_orb(2)
orb%vec(4) = orb%vec(4) - ks4 * xy_orb(1)
if (track_spn) then
  ff = orb%time_dir * at_sign * sign_z_vel * hard_ele%value(bs_field$) / 2
  call rotate_spin_given_field (orb, sign_z_vel, -[xy_orb(1), xy_orb(2), 0.0_rp] * ff)
endif
orb%vec(2) = orb%vec(2) + ks4 * xy_orb(2)
orb%vec(4) = orb%vec(4) - ks4 * xy_orb(1)

end subroutine apply_this_sol_fringe

end subroutine apply_element_edge_kick
