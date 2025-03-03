!+
! Subroutine update_fibre_from_ele (ele, survey_needed)
!
! Routine to update a fibre when the associated Bmad ele has been modified.
!
! Input:
!   ele           -- ele_struct: Element with corresponding PTC fibre.
!
! Output:
!   ele%ptc_fibre -- PTC fibre.
!   survey_needed -- logical: Set True if a call to survey will be needed.
!                      Calling survey is avoided in this routine to save time if multiple elements
!                      are being updated.
!-

subroutine update_fibre_from_ele (ele, survey_needed)

use ptc_interface_mod, dum1 => update_fibre_from_ele
use pointer_lattice, dum2 => twopi, dum3 => pi, dum4 => sqrt

implicit none

type (ele_struct), target :: ele, m_ele
type (fibre), pointer :: fib
type (branch_struct), pointer :: branch
type (keywords) ptc_key
type (element), pointer :: mag
type (elementp), pointer :: magp
type (magnet_chart), pointer :: mp, mpp
type (work) ptc_work

real(rp) value, hk, vk, phi_tot, fh, volt, delta_p, e1, e2
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) a_ptc(0:n_pole_maxx), b_ptc(0:n_pole_maxx)
real(rp), pointer :: val(:)

integer i, n, ns, ix, ix_pole_max, cavity_type
logical survey_needed, kill_spin_fringe

character(*), parameter :: r_name = 'update_fibre_from_ele'

! Note: change_settings_fibre (in PTC) is an alternative way of setting some element parameters.

survey_needed = .false.
branch => pointer_to_branch(ele)
fib => ele%ptc_fibre
val => ele%value
FEED_P0C = .true.

mag  => fib%mag
magp => fib%magp

mp => mag%p
mpp => magp%p

ptc_work = fib
ptc_work = 1e-9_rp * val(p0c$) - mp%p0c 
fib = ptc_work

cavity_type = not_set$
if (ele%key == rfcavity$ .or. ele%key == lcavity$) then
  cavity_type = nint(ele%value(cavity_type$))
  if (ele%tracking_method == bmad_standard$) cavity_type = ptc_standard$
endif

!

if (ele%key == sbend$) then
  call set_real_all (mp%ld, mpp%ld, val(l$))
  call set_real_all (mp%lc, mpp%lc, val(l_chord$))
elseif (cavity_type /= standing_wave$) then
  call set_real_all (mp%ld, mpp%ld, val(l$))
  call set_real_all (mp%lc, mpp%lc, val(l$))
endif

!

if (ele%key == marker$) return

! Must set all poles even if zero since they might have been non-zero beforehand.
! Note: On ptc side bn(1) is error field when creating a fibre but is total field when fibre is being modified.	 

! "0" argument in add routine means set k/ks to value given.
! As opposed to "1" which means add to existing value.

! Magnetic

a_ptc = 0
b_ptc = 0

if (associated(fib%mag%an)) then
  n = size(fib%mag%an)
  a_ptc(0:n-1) = fib%mag%an
  b_ptc(0:n-1) = fib%mag%bn
endif

call ele_to_ptc_magnetic_bn_an (ele, b_pole, a_pole)
if (ele%key == sbend$) b_pole(0) = b_pole(0) + ele%value(g$)	 

do i = n_pole_maxx, 0, -1
  if (b_pole(i) /= b_ptc(i)) call add (fib,  (i+1), 0, b_pole(i))
  if (a_pole(i) /= a_ptc(i)) call add (fib, -(i+1), 0, a_pole(i))
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
  if (b_pole(i) /= 0) call add (fib,  (i+1), 0, fh*b_pole(i), electric = .true.)
  if (a_pole(i) /= 0) call add (fib, -(i+1), 0, fh*a_pole(i), electric = .true.)
enddo

! Note: ele_to_an_bn takes care of such stuff as sextupole strength conversion so
! only have to worry about non-multipole components here.

! PTC does not like %nst for a marker to be anything other than 1.
! And don't set drifts.

if (nint(val(num_steps$)) /= mp%nst .and. (mag%kind /= kind0 .and. mag%kind /= kind1)) then
  ns = nint(val(num_steps$))
  if (cavity_type == standing_wave$) ns = min(5, ns) ! To avoid bug
  call set_integer (mp%nst, mpp%nst, ns)
  if (mag%kind /= kind4 .and. mag%kind /= kind21) call add (fib, 1, 1, 0.0_rp)  ! Triggers recompute of matrices.
  survey_needed = .true.
endif

n = nint(val(integrator_order$))
if (n /= 0) call set_integer (mp%method, mpp%method, n)

select case (ele%key)

case (solenoid$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (sol_quad$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (rfcavity$, lcavity$, crab_cavity$)
  phi_tot = twopi * (val(phi0$) + val(phi0_multipass$) + val(phi0_err$) + val(phi0_autoscale$))
  if (ele%key == lcavity$) phi_tot = pi / 2 - twopi * phi_tot

  select case (cavity_type)
  case (standing_wave$)
    call set_real_all (mp%ld, mpp%ld, val(l_active$))
    call set_real_all (mp%lc, mpp%lc, val(l_active$))
    call set_real (mag%h1, magp%h1, (val(l$) - val(l_active$)) / 2)
    call set_real (mag%h2, magp%h2, (val(l$) - val(l_active$)) / 2)
    volt = 2d-6 * e_accel_field(ele, voltage$)
  case default
    volt = 1d-6 * e_accel_field(ele, voltage$) 
  end select

  mag%lag = phi_tot  ! There is no magp%lat
  call set_real (mag%phas, magp%phas, -phi_tot)
  call set_real (mag%volt, magp%volt, volt)
  call set_real (mag%freq, magp%freq, val(rf_frequency$))

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

  ix = both_ends$
  if (attribute_index(ele, 'FRINGE_AT') > 0) ix = nint(ele%value(fringe_at$))
  kill_spin_fringe = is_false(ele%value(spin_fringe_on$))

  call set_logic(mag%p%kill_ent_fringe, magp%p%kill_ent_fringe, (ix == exit_end$ .or. ix == no_end$))
  call set_logic(mag%p%kill_exi_fringe, magp%p%kill_exi_fringe, (ix == entrance_end$ .or. ix == no_end$))

  call set_logic(mag%p%kill_ent_spin, magp%p%kill_ent_spin, (ix == exit_end$ .or. ix == no_end$ .or. kill_spin_fringe))
  call set_logic(mag%p%kill_exi_spin, magp%p%kill_exi_spin, (ix == entrance_end$ .or. ix == no_end$ .or. kill_spin_fringe))

  ix = nint(ele%value(fringe_type$))

  e1 = ele%value(e1$)
  if (ptc_key%list%kill_ent_fringe .or. ix == none$) e1 = 0

  e2 = ele%value(e2$)
  if (ptc_key%list%kill_exi_fringe .or. ix == none$) e2 = 0

  call set_real_all (mag%p%edge(1), magp%p%edge(1), e1)
  call set_real_all (mag%p%edge(2), magp%p%edge(2), e2)

  !!  if (nint(ele%value(ptc_field_geometry$)) == straight$) then
  !!    ptc_key%list%t1   = e1 - ele%value(angle$)/2
  !!    ptc_key%list%t2   = e2 - ele%value(angle$)/2
  !!  endif

end select

! Fringe

if (ele%key == sbend$) then
  ix = nint(val(ptc_fringe_geometry$))
  call set_logic (mp%bend_fringe, mpp%bend_fringe, (ix == x_invariant$))

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (mp%permfringe, mpp%permfringe, 0)
  case (basic_bend$, linear_edge$)
    call set_integer (mp%permfringe, mpp%permfringe, 0)
  case (full$)
    call set_integer (mp%permfringe, mpp%permfringe, 1)
  case (hard_edge_only$)
    call set_integer (mp%permfringe, mpp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (mp%permfringe, mpp%permfringe, 2)
  case (sad_full$)
    call set_integer (mp%permfringe, mpp%permfringe, 3)
  end select

elseif (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute
  call set_logic (mp%bend_fringe, mpp%bend_fringe, .false.)

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (mp%permfringe, mpp%permfringe, 0)
  case (hard_edge_only$)
    call set_integer (mp%permfringe, mpp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (mp%permfringe, mpp%permfringe, 2)
  case (full$)
    call set_integer (mp%permfringe, mpp%permfringe, 3)
  end select

  if (ele%key == sad_mult$ .and. val(l$) == 0) call set_integer (mp%permfringe, mpp%permfringe, 0)
endif

if (attribute_index(ele, 'FRINGE_AT') > 0) then  ! If fringe_at is a valid attribute
  ix = nint(val(fringe_at$))
  call set_logic (mp%kill_ent_fringe, mpp%kill_ent_fringe, (ix == downstream_end$ .or. ix == no_end$))
  call set_logic (mp%kill_exi_fringe, mpp%kill_exi_fringe, (ix == upstream_end$ .or. ix == no_end$))
endif

! misalign

call misalign_ptc_fibre (ele, .true., fib, .true.)

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

