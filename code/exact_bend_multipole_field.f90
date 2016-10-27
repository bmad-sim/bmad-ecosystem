!+
! Subroutine exact_bend_multipole_field (ele, param, orbit, local_ref_frame, field, potential, calc_dfield)
!
! Routine to calculate the electric and magnetic field at a given point in a bend element due to any multipoles.
! The field due to a multipole in a bend is different from a straight element since Maxwell's equations
! are modified due to the curvature of the reference orbit.
!
! Note: The returned field does not include the bend field (g) itself or the bend field error (g_err).
! Note: init_exact_bend_multipole_coefs must be called before this routine is called.
!
! This routine is a port of the getelectricr routine from PTC.
!
! Input:
!   ele             -- ele_stuct: Bend element.
!   param           -- lat_param_struct: Lattice branch parameters.
!   orbit           -- coord_struct: particle position.
!   local_ref_frame -- logical: Is the particle position in the local element ref 
!                         frame (as opposed to the lab frame)?
!   calc_dfield     -- logical, optional: If present and True then calculate the field derivatives.
!
! Output:
!   field           -- em_field_struct: Field
!   potential       -- em_potential_struct: Potential. Note that potential%A (mag vector potential) is not computed.
!-

subroutine exact_bend_multipole_field (ele, param, orbit, local_ref_frame, field, potential, calc_dfield)

use bmad_interface, except_dummy => exact_bend_multipole_field

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) orbit
type (em_field_struct) field
type (em_potential_struct) potential
type (exact_bend_multipole_struct), pointer :: eb

real(rp) x1, x3, bx, by, btx, bty, btyt, phit, phi_mag
real(rp) ex, etx, ey, ety, b0
real(rp) dbx_dx, dbx_dy, dby_dx, dby_dy, dex_dx, dex_dy, dey_dx, dey_dy
real(rp) dbtx_dy, dbty_dy, detx_dy, dety_dy

integer j, m, a, k

logical local_ref_frame
logical, optional :: calc_dfield
logical has_nonzero_pole, do_dfield_calc

! potential%phi_b is the integrated magnetic fringe field (the scalar potential)

x1 = orbit%vec(1)
x3 = orbit%vec(3)
b0 = ele%value(g$)

field = em_field_struct()
potential = em_potential_struct()

if (.not. associated(ele%exact_bend_multipole)) return

do_dfield_calc = logic_option(.false., calc_dfield)
eb => ele%exact_bend_multipole

bx = 0.0_rp
by = 0.0_rp
ex = 0.0_rp
ey = 0.0_rp

if (do_dfield_calc) then
  dbx_dx = 0;    dbx_dy = 0;    dby_dx = 0;    dby_dy = 0
  dex_dx = 0;    dex_dy = 0;    dey_dx = 0;    dey_dy = 0
  dbtx_dy = 0;   dbty_dy = 0;   detx_dy = 0;   dety_dy = 0
endif

m = n_pole_maxx + 1  ! PTC multipole index off by 1 compaired to Bmad.
k = 0

do a = m, 1, -1
  etx = 0.0_rp
  ety = 0.0_rp
  btx = 0.0_rp
  bty = 0.0_rp
  phit = 0.0_rp
  phi_mag = 0.0_rp
  do j = m - a, 1, -1
    k = k + 1

    if (do_dfield_calc) then
      dbtx_dy = dbtx_dy * x3 + (btx + eb%bf_x(k))
      dbty_dy = dbty_dy * x3 + (bty + eb%bf_y(k))
      detx_dy = detx_dy * x3 + (etx + eb%e_x(k))
      dety_dy = dety_dy * x3 + (ety + eb%e_y(k))
    endif

    btx = (btx + eb%bf_x(k))*x3 
    bty = (bty + eb%bf_y(k))*x3
    etx = (etx + eb%e_x(k))*x3 
    ety = (ety + eb%e_y(k))*x3
    phit = (phit + eb%phi(k))*x3
    phi_mag = (phi_mag + eb%vm(k))*x3
  enddo

  k = k + 1
  btx = (btx + eb%bf_x(k))
  bty = (bty + eb%bf_y(k))
  etx = (etx + eb%e_x(k))
  ety = (ety + eb%e_y(k))
  phit = (phit + eb%phi(k))
  phi_mag = (phi_mag + eb%vm(k))

  if (do_dfield_calc) then
    dbx_dx = dbx_dx * x1 + (bx + btx)
    dbx_dy = (dbx_dy + dbtx_dy) * x1
    dby_dx = dby_dx * x1 + (by + bty)
    dby_dy = (dby_dy + dbty_dy) * x1
    dex_dx = dex_dx * x1 + (ex + etx)
    dex_dy = (dex_dy + detx_dy) * x1
    dey_dx = dey_dx * x1 + (ey + ety)
    dey_dy = (dey_dy + dety_dy) * x1
  endif

  bx = (bx + btx)*x1
  by = (by + bty)*x1
  ex = (ex + etx)*x1
  ey = (ey + ety)*x1
  potential%phi = (potential%phi + phit)*x1
  potential%phi_b = (potential%phi_b + phi_mag)*x1
enddo

btx = 0.0_rp
bty = 0.0_rp
etx = 0.0_rp
ety = 0.0_rp
phit = 0.0_rp
phi_mag = 0.0_rp

do j = m, 1, -1
  k = k + 1

  if (do_dfield_calc) then
    dbtx_dy = dbtx_dy * x3 + (btx + eb%bf_x(k))
    dbty_dy = dbty_dy * x3 + (bty + eb%bf_y(k))
    detx_dy = detx_dy * x3 + (etx + eb%e_x(k))
    dety_dy = dety_dy * x3 + (ety + eb%e_y(k))
  endif

  btx = (btx + eb%bf_x(k))*x3
  bty = (bty + eb%bf_y(k))*x3
  etx = (etx + eb%e_x(k))*x3
  ety = (ety + eb%e_y(k))*x3
  phit = (phit + eb%phi(k))*x3
  phi_mag = (phi_mag + eb%vm(k))*x3
enddo

k = k + 1

if (do_dfield_calc) then
  dbx_dy = dbx_dy + dbtx_dy
  dby_dy = dby_dy + dbty_dy
  dex_dy = dex_dy + detx_dy
  dey_dy = dey_dy + dety_dy
endif

bx = bx + btx + eb%bf_x(k)
by = by + bty + eb%bf_y(k)
ex = ex + etx + eb%e_x(k) 
ey = ey + ety + eb%e_y(k)

potential%phi = potential%phi + phit + eb%phi(k)
potential%phi_b = potential%phi_b + phi_mag + eb%vm(k)

! Kick is:
!  kick_x = -by*(1.0_rp + b0*x1)
!  kick_y =  bx*(1.0_rp + b0*x1)

if (do_dfield_calc) then
  field%db(1,1:2) = [dbx_dx, dbx_dy]
  field%db(2,1:2) = [dby_dx, dby_dy]
  field%de(1,1:2) = [dex_dx, dex_dy]
  field%de(2,1:2) = [dey_dx, dey_dy]
endif

field%b(1) = bx
field%b(2) = by

field%e(1) = ex
field%e(2) = ey

potential%phi   = potential%phi   * c_light / ele%value(p0c$)
potential%phi_B = potential%phi_B * ele%value(p0c$) / c_light

end subroutine exact_bend_multipole_field

!------------------------------------------------------------------------------------------------------
! This is a port of the getanbnr routine from the forest library
! Note: k1 and k2 multipoles will be folded in to the multipole array.

subroutine init_exact_bend_multipole_coefs (ele, param, local_ref_frame, has_nonzero_pole)

use multipole_mod, except_dummy => init_exact_bend_multipole_coefs
use s_status, only: s_b_from_v, s_e, sector_nmul_max

implicit none

type (ele_struct), target :: ele
type (exact_bend_multipole_struct), pointer :: eb
type (lat_param_struct) param

real(rp) b0, f, a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)

integer i, j, k, pow, n_mono, n_max

logical :: local_ref_frame, has_nonzero_pole, has_nonzero_mag, has_nonzero_elec

!

call multipole_ele_to_ab (ele, .not. local_ref_frame, has_nonzero_mag, a_pole, b_pole)
call multipole_ele_to_ab (ele, .not. local_ref_frame, has_nonzero_elec, a_pole_elec, b_pole_elec, electric$)

has_nonzero_pole = (has_nonzero_mag .or. has_nonzero_elec .or. ele%value(k1$) /= 0 .or. ele%value(k2$) /= 0) 
if (.not. has_nonzero_pole) return

if (.not. associated(ele%exact_bend_multipole)) allocate(ele%exact_bend_multipole)
eb => ele%exact_bend_multipole

if (ele%value(k1$) /= 0) eb%bn(1) = eb%bn(1) + ele%value(k1$) * ele%value(l$)
if (ele%value(k2$) /= 0) eb%bn(2) = eb%bn(2) + ele%value(k2$) * ele%value(l$) / 2

n_mono = s_b_from_v%n_mono
n_max = sector_nmul_max

call re_allocate (eb%bf_x, n_mono)
call re_allocate (eb%bf_y, n_mono)
call re_allocate (eb%vm,   n_mono)

call re_allocate (eb%e_x, n_mono)
call re_allocate (eb%e_y, n_mono)
call re_allocate (eb%phi, n_mono)

b0 = ele%value(g$)
a_pole = a_pole / ele%value(l$)
b_pole = b_pole / ele%value(l$)

if (eb%b0 == b0 .and. all(a_pole == eb%an) .and. all(b_pole == eb%bn) .and. &
                          all(a_pole_elec == eb%ae) .and. all(b_pole_elec == eb%be)) return  ! Nothing to do if already set.

eb%b0 = b0
eb%an = a_pole
eb%bn = b_pole
eb%ae = a_pole_elec
eb%be = b_pole_elec

do i = n_pole_maxx, 0, -1
  if (eb%an(i) == 0 .and. eb%bn(i) == 0 .and. eb%ae(i) == 0 .and. eb%be(i) == 0) cycle
  eb%n_pole_max = i
  exit
enddo

eb%bf_x = 0
eb%bf_y = 0

do i = 1, n_max
  do j = 1, n_mono
    k = s_b_from_v%i(j) + s_b_from_v%j(j)
    pow = k + 1 - i
    if(k + 1 >= i) then
       eb%e_x(j)=eb%e_x(j) + (eb%ae(i-1)*s_e%a_x(i,j) + eb%be(i-1)*s_e%b_x(i,j))*b0**pow
       eb%e_y(j)=eb%e_y(j) + (eb%ae(i-1)*s_e%a_y(i,j) + eb%be(i-1)*s_e%b_y(i,j))*b0**pow
       eb%bf_x(j) = eb%bf_x(j) + (eb%an(i-1)*s_b_from_v%a_x(i,j) + eb%bn(i-1)*s_b_from_v%b_x(i,j))*b0**pow
       eb%bf_y(j) = eb%bf_y(j) + (eb%an(i-1)*s_b_from_v%a_y(i,j) + eb%bn(i-1)*s_b_from_v%b_y(i,j))*b0**pow
    endif
  enddo
enddo

eb%vm = 0
eb%phi = 0
do i = 1, n_max
  do j = 1, n_mono
    k = s_b_from_v%i(j) + s_b_from_v%j(j)
    pow = k - i
    if(k >= i) then
       eb%phi(j)= eb%phi(j) + (eb%ae(i-1)*s_e%va(i,j) + eb%be(i-1)*s_e%vb(i,j))*b0**pow
       eb%vm(j) = eb%vm(j)  + (eb%an(i-1)*s_b_from_v%va(i,j) + eb%bn(i-1)*s_b_from_v%vb(i,j))*b0**pow
    endif
  enddo
enddo

f = charge_of(param%particle) * ele%value(p0c$) / c_light
eb%bf_x = eb%bf_x * f
eb%bf_y = eb%bf_y * f
eb%phi  = eb%phi * f

end subroutine init_exact_bend_multipole_coefs

