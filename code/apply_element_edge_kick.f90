!+
! Subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix)
!
! Subroutine, used with runge_kutta and boris tracking, to track through the edge fringe field of an element.
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
!   orb         -- Coord_struct: Starting coords in element reference frame.
!   fringe_info -- fringe_edge_info_struct: Fringe information.
!   track_ele   -- ele_struct: Element being tracked through. Is different from fringe_info%hard_ele
!                    when there are superpositions and track_ele can be a super_slave of fringe_info%hard_ele.
!   param       -- lat_param_struct: lattice parameters.
!   track_spin  -- logical: Track the spin?
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before fringe.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orb        -- Coord_struct: Coords after application of the edge fringe field.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including fringe.
!-

subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix)

use track1_mod, except_dummy => apply_element_edge_kick

implicit none

type (ele_struct), target :: ele, track_ele
type (coord_struct) orb
type (lat_param_struct) param
type (em_field_struct) field
type (fringe_edge_info_struct) fringe_info
type (ele_struct), pointer :: hard_ele

real(rp), optional :: mat6(6,6)
real(rp) f, l_drift, ks, s_edge, s, phi, omega(3), pc
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
complex(rp) xiy, c_vec

integer physical_end, dir, i, fringe_at, at_sign, sign_z_vel, particle_at

logical, optional :: make_matrix
logical finished, track_spin, track_spn, has_nonzero_elec

! The setting of fringe_info%hard_location is used by calc_next_fringe_edge to calculate the next fringe location.


particle_at = fringe_info%particle_at

if (associated(fringe_info%hard_ele)) then
  hard_ele => fringe_info%hard_ele
  s_edge = fringe_info%s_edge_hard

  if (particle_at == first_track_edge$) then
    fringe_info%hard_location = inside$
  else
    dir = orb%direction
    if (hard_ele%value(l$) < 0) dir = -dir
    if (dir == 1) then
      fringe_info%hard_location = downstream_end$
    else
      fringe_info%hard_location = upstream_end$
    endif
  endif

else
  hard_ele => track_ele
  if (particle_at == first_track_edge$) then
    s_edge = 0
  else
    s_edge = track_ele%value(l$)
  endif
endif

! Custom edge kick?

call apply_element_edge_kick_hook (orb, fringe_info, track_ele, param, finished, mat6, make_matrix)
if (finished) return

!------------------------------------------------------------------------------------
! Only need this routine when the field_calc or tracking is bmad_standard.
! Note: track_ele, if a slave, will have track_ele%field_calc = refer_to_lords$.

if (hard_ele%field_calc /= bmad_standard$ .and. hard_ele%tracking_method /= bmad_standard$) return

physical_end = physical_ele_end (particle_at, orb%direction, track_ele%orientation)
fringe_at = nint(track_ele%value(fringe_at$))
if (.not. at_this_ele_end(physical_end, fringe_at)) return
track_spn = (track_spin .and. bmad_com%spin_tracking_on .and. is_true(hard_ele%value(spin_fringe_on$)))

if (hard_ele%key == e_gun$ .and. physical_end == entrance_end$) return ! E_gun does not have an entrance kick

if (particle_at == first_track_edge$) then
  at_sign = 1
else
  at_sign = -1
endif

sign_z_vel = orb%direction * track_ele%orientation
call multipole_ele_to_ab (hard_ele, .false., has_nonzero_elec, a_pole, b_pole, electric$)

! Static electric longitudinal field
! Note: magnetic fringe effects, which may shift the particle's (x, y) position, need to
! be applied before the electric fringe on entrance and after the electric fringe on exit in
! order to conserve the particle's energy.

if (particle_at == second_track_edge$) call electric_longitudinal_fringe()

! Static magnetic and electromagnetic fringes

select case (hard_ele%key)
case (quadrupole$)
  if (particle_at == first_track_edge$) then
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb)
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb)
  else
    call soft_quadrupole_edge_kick (hard_ele, param, particle_at, orb)
    call hard_multipole_edge_kick (hard_ele, param, particle_at, orb)
  endif

case (sbend$)
  call bend_edge_kick (hard_ele, param, particle_at, orb, mat6, make_matrix, track_spn)

! Note: Cannot trust hard_ele%value(ks$) here since element may be superimposed with an lcavity.
! So use hard_ele%value(bs_field$).

case (solenoid$, sol_quad$, bend_sol_quad$)
  ks = at_sign * relative_tracking_charge(orb, param) * hard_ele%value(bs_field$) * c_light / orb%p0c
  orb%vec(2) = orb%vec(2) + ks * orb%vec(3) / 2
  orb%vec(4) = orb%vec(4) - ks * orb%vec(1) / 2
  if (track_spn) then
    f = at_sign * sign_z_vel * hard_ele%value(bs_field$) / 2
    call rotate_spin_given_field (orb, sign_z_vel, -[orb%vec(1), orb%vec(3), 0.0_rp] * f)
  endif

case (lcavity$, rfcavity$, e_gun$)

  ! Add on bmad_com%significant_length to make sure we are just inside the cavity.
  f = at_sign * charge_of(orb%species) / (2 * orb%p0c)
  s = s_edge

  if (at_this_ele_end(physical_end, nint(hard_ele%value(fringe_at$)))) then
    if (particle_at == first_track_edge$) then
      s = s + bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, orb, .true., field)
    else
      s = s - bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, orb, .true., field)
    endif

    orb%vec(2) = orb%vec(2) - field%e(3) * orb%vec(1) * f + c_light * field%b(3) * orb%vec(3) * f
    orb%vec(4) = orb%vec(4) - field%e(3) * orb%vec(3) * f - c_light * field%b(3) * orb%vec(1) * f

    if (track_spn) then
      f = at_sign * charge_of(orb%species) / 2.0_rp
      call rotate_spin_given_field (orb, sign_z_vel, &
                                           -[orb%vec(1), orb%vec(3), 0.0_rp] * (f * field%b(3)), &
                                           -[orb%vec(1), orb%vec(3), 0.0_rp] * (f * field%e(3)))
    endif

    ! orb%phase(1) is set by em_field_calc.

    call rf_coupler_kick (hard_ele, param, particle_at, orb%phase(1), orb)
  endif

case (elseparator$)
  ! Longitudinal fringe field
  if (hard_ele%value(l$) /= 0) then
    f = at_sign * charge_of(orb%species) * (hard_ele%value(p0c$) / hard_ele%value(l$))
    phi = f * (hard_ele%value(hkick$) * orb%vec(1) + hard_ele%value(vkick$) * orb%vec(3))
    call apply_energy_kick (phi, orb)
    if (track_spn) then
      call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, phi])
    endif
  endif
end select

! Entrance Static electric longitudinal field

if (particle_at == first_track_edge$) call electric_longitudinal_fringe()


!--------------------------------------------------------------------------------
contains

subroutine electric_longitudinal_fringe()

if (has_nonzero_elec) then
  xiy = 1
  c_vec = cmplx(orb%vec(1), orb%vec(3), rp)
  do i = 0, max_nonzero(0, a_pole, b_pole)
    xiy = xiy * c_vec
    if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
    phi = at_sign * charge_of(orb%species) * real(cmplx(b_pole(i), -a_pole(i), rp) * xiy) / (i + 1)
    call apply_energy_kick (phi, orb)

    if (track_spn) then
      call rotate_spin_given_field (orb, sign_z_vel, EL = [0.0_rp, 0.0_rp, phi])
    endif
  enddo
endif

end subroutine electric_longitudinal_fringe

end subroutine apply_element_edge_kick
