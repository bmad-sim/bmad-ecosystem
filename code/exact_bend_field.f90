!+
! Subroutine exact_bend_field (ele, param, orbit, local_ref_frame, field, potential, return_kick)
!
! Routine to calculate the electric and magnetic field at a given point in a bend element with multipoles.
! The field due to a multipole in a bend is different from a straight element due Maxwell's equations
! being modified due to the curvature of the reference orbit.
!
! This routine is a port of the getelectricr routine from the forest library
!
! Input:
!   ele             -- ele_stuct: Bend element.
!   param           -- lat_param_struct: Lattice branch parameters.
!   orbit           -- coord_struct: particle position.
!   local_ref_frame -- logical: Is the particle position in the local element ref 
!                         frame (as opposed to the lab frame)?
!   return_kick     -- logical, optional: Return the kick instead of the field? Default is False.
!
! Output:
!   field           -- em_field_struct: Field
!   potential       -- em_potential_struct: Potential
!-

subroutine exact_bend_field (ele, param, orbit, local_ref_frame, field, potential, return_kick)

use bmad_interface, except_dummy => exact_bend_field
use s_status, only: s_b_from_v, s_e, sector_nmul_max

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) orbit
type (em_field_struct) field
type (em_potential_struct) potential
type (exact_bend_struct), pointer :: eb

real(rp) x1, x3, bx, by, btx, bty, btyt, phit, phi_mag
real(rp) ex, etx, ey, ety, b0

integer j, m, a, k

logical local_ref_frame
logical, optional :: return_kick
logical has_nonzero_pole

! potential%phi_b is the integrated magnetic fringe field (the scalar potential)

x1 = orbit%vec(1)
x3 = orbit%vec(3)
b0 = ele%value(g$)

field = em_field_struct()
potential = em_potential_struct()

call init_exact_bend_coefs(ele, param, has_nonzero_pole)
if (.not. has_nonzero_pole) return

eb => ele%exact_bend

bx = 0.0_rp
by = 0.0_rp
ex = 0.0_rp
ey = 0.0_rp

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
    btx = (btx + eb%bf_x(k))*x3  !x1
    bty = (bty + eb%bf_y(k))*x3
    etx = (etx + eb%e_x(k))*x3  !x1
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
  btx = (btx + eb%bf_x(k))*x3
  bty = (bty + eb%bf_y(k))*x3
  etx = (etx + eb%e_x(k))*x3
  ety = (ety + eb%e_y(k))*x3
  phit = (phit + eb%phi(k))*x3
  phi_mag = (phi_mag + eb%vm(k))*x3
enddo

k = k + 1
bx = bx + btx + eb%bf_x(k)  ! + x3
by = by + bty + eb%bf_y(k)  ! + x3
ex = ex + etx + eb%e_x(k)  ! + x3
ey = ey + ety + eb%e_y(k)  ! + x3
potential%phi = potential%phi + phit + eb%phi(k)  ! + x3
potential%phi_b = potential%phi_b + phi_mag + eb%vm(k)  ! + x3

if (logic_option(.false., return_kick)) then
  field%b(1) = -by*(1.0_rp + b0*x1)
  field%b(2) =  bx*(1.0_rp + b0*x1)
  field%b(3) = 0.0_rp

else
  field%b(1) = bx
  field%b(2) = by
  field%b(3) = 0.0_rp
endif

field%e(1) = ex/ele%value(p0c$)
field%e(2) = ey/ele%value(p0c$)
field%e(3) = 0.0_rp
potential%phi = potential%phi/ele%value(p0c$)

!--------------------------------------------------------------------------
contains

! This is a port of the getanbnr routine from the forest library

subroutine init_exact_bend_coefs (ele, param, has_nonzero_pole)

implicit none

type (ele_struct), target :: ele
type (exact_bend_struct), pointer :: eb
type (lat_param_struct) param

real(rp) b0, f
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer i, j, k, pow, n_mono, n_max
logical has_nonzero_pole

!

call multipole_ele_to_ab (ele, .not. local_ref_frame, has_nonzero_pole, a_pole, b_pole)
if (.not. has_nonzero_pole) return

n_mono = s_b_from_v%n_mono
n_max = sector_nmul_max

if (.not. associated(ele%exact_bend)) allocate(ele%exact_bend)
eb => ele%exact_bend

call re_allocate (eb%bf_x, n_mono)
call re_allocate (eb%bf_y, n_mono)
call re_allocate (eb%vm,   n_mono)

call re_allocate (eb%e_x, n_mono)
call re_allocate (eb%e_y, n_mono)
call re_allocate (eb%phi, n_mono)

b0 = ele%value(g$)

if (eb%b0 == b0 .and. all(a_pole == eb%an) .and. all(b_pole == eb%bn)) return  ! Nothing to do if already set.

eb%b0 = b0
eb%an = a_pole
eb%bn = b_pole

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

end subroutine init_exact_bend_coefs

end subroutine exact_bend_field

