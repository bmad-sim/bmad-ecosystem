!+
! Subroutine update_fibre_from_ele (ele)
!
! Routine to update a fibre when the associated Bmad ele has been modified.
!
! Input:
!   ele           -- ele_struct: Element with corresponding PTC fibre.
!
! Output:
!   ele%ptc_fibre -- PTC fibre.
!-

subroutine update_fibre_from_ele (ele)

use ptc_interface_mod, dum1 => update_fibre_from_ele
use pointer_lattice, dum2 => twopi, dum3 => pi, dum4 => sqrt

type (ele_struct), target :: ele, m_ele
type (branch_struct), pointer :: branch
type (keywords) ptc_key
type (element), pointer :: mag
type (elementp), pointer :: magp
type (magnet_chart), pointer :: p, pp

real(rp) value, hk, vk, phi_tot, fh
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp), pointer :: val(:)

integer i, ix, ix_pole_max

character(*), parameter :: r_name = 'update_fibre_from_ele'

! "0" argument in add routine means set k/ks to value given.
! As opposed to "1" which means add to existing value.

branch => pointer_to_branch(ele)
val => ele%value

! Must set all poles even if zero since they might have been non-zero beforehand.
! Note: On ptc side bn(1) is error field when creating a fibre but is total field when fibre is being modified.	 

! Magnetic

call ele_to_ptc_magnetic_an_bn (ele, branch%param, b_pole, a_pole) ! Yes arg order is b_pole, a_pole.
if (ele%key == sbend$) b_pole(1) = b_pole(1) + ele%value(g$)	 

do i = n_pole_maxx, 0, -1
  if (b_pole(i) /= 0) call add (ele%ptc_fibre,  (i+1), 0, b_pole(i))
  if (a_pole(i) /= 0) call add (ele%ptc_fibre, -(i+1), 0, a_pole(i))
enddo

! Electric. Notice that PTC assumes horizontally_pure bend multipoles

fh = 1d-9 * sign_of(charge_of(branch%param%particle)) / VOLT_C
if (ele%key == sbend$ .and. nint(val(exact_multipoles$)) == vertically_pure$) then
  call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$, include_kicks$)
  ! Notice that a_pole and b_pole are reversed for electric multipoles.
  call convert_bend_exact_multipole(ele%value(g$), horizontally_pure$, b_pole, a_pole)
else
  call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$)
endif

if (ele%key == elseparator$) then
  a_pole(0) = a_pole(0) + val(vkick$) * val(p0c$) / val(l$)
  b_pole(0) = b_pole(0) + val(hkick$) * val(p0c$) / val(l$)
endif

do i = 0, n_pole_maxx
  if (b_pole(i) /= 0) call add (ele%ptc_fibre,  (i+1), 0, fh*b_pole(i), electric = .true.)
  if (a_pole(i) /= 0) call add (ele%ptc_fibre, -(i+1), 0, fh*a_pole(i), electric = .true.)
enddo

! Note: ele_to_an_bn takes care of such stuff as sextupole strength conversion so
! only have to worry about non-multipole components here.

mag  => ele%ptc_fibre%mag
magp => ele%ptc_fibre%magp

p => mag%p
pp => magp%p

! call set_integer (p%method, pp%method, nint(val(integrator_order$)))
! call set_integer (p%nst, pp%nst, nint(val(num_steps$)))

select case (ele%key)

case (solenoid$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (sol_quad$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (rfcavity$, lcavity$)
  phi_tot = twopi * (val(phi0$) + val(phi0_multipass$) + val(phi0_err$) + val(phi0_autoscale$))
  if (ele%key == lcavity$) phi_tot = pi / 2 - twopi * phi_tot
  call set_real (mag%phas, magp%phas, -mag%lag)

  call set_real (mag%volt, magp%volt, 2d-6 * e_accel_field(ele, voltage$))

case (sad_mult$)

  if (val(l$) /= 0) then
    call set_real (mag%b_sol, magp%b_sol, val(ks$))
    call set_real (mag%va, magp%va, -sign(sqrt(24 * abs(val(fq1$))), val(fq1$)))
    call set_real (mag%vs, magp%vs, val(fq2$))
  endif

  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (sbend$)
  call set_real (mag%hgap(1), magp%hgap(1), val(hgap$))
  call set_real (mag%fint(1), magp%fint(1), val(fint$))
  call set_real (mag%hgap(2), magp%hgap(2), val(hgapx$))
  call set_real (mag%fint(2), magp%fint(2), val(fintx$))
  call set_real_all (mag%p%edge(1), magp%p%edge(1), val(e1$))
  call set_real_all (mag%p%edge(2), magp%p%edge(2), val(e2$))

end select

! Fringe

if (ele%key == sbend$) then
  ix = nint(val(ptc_fringe_geometry$))
  call set_logic (p%bend_fringe, pp%bend_fringe, (ix == x_invariant$))

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (basic_bend$, linear_edge$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (full$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (hard_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 2)
  case (sad_full$)
    call set_integer (p%permfringe, pp%permfringe, 3)
  end select

elseif (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute
  call set_logic (p%bend_fringe, pp%bend_fringe, .false.)

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (hard_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 2)
  case (full$)
    call set_integer (p%permfringe, pp%permfringe, 3)
  end select

  if (ele%key == sad_mult$ .and. val(l$) == 0) call set_integer (p%permfringe, pp%permfringe, 0)
endif

if (attribute_index(ele, 'FRINGE_AT') > 0) then  ! If fringe_at is a valid attribute
  ix = nint(val(fringe_at$))
  call set_logic (p%kill_ent_fringe, pp%kill_ent_fringe, (ix == downstream_end$ .or. ix == no_end$))
  call set_logic (p%kill_exi_fringe, pp%kill_exi_fringe, (ix == upstream_end$ .or. ix == no_end$))
endif

! misalign

call misalign_ele_to_fibre (ele, .true., ele%ptc_fibre)

ele%bookkeeping_state%ptc = ok$

!-------------------------------------------------------------------------
contains

subroutine set_real (to1, to2, value)
real(rp) value, to1
type(real_8) to2

to1 = value
to2 = value

end subroutine set_real

!-------------------------------------------------------------------------
! contains

subroutine set_real_all (to1, to2, value)
real(rp) value, to1, to2

to1 = value
to2 = value

end subroutine set_real_all

!-------------------------------------------------------------------------
! contains

subroutine set_integer (to1, to2, value)
integer to1, to2, value

to1 = value
to2 = value

end subroutine set_integer

!-------------------------------------------------------------------------
! contains

subroutine set_logic (to1, to2, value)
logical to1, to2, value

to1 = value
to2 = value

end subroutine set_logic

end subroutine update_fibre_from_ele

