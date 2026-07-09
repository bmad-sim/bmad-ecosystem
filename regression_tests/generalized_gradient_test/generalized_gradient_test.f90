!+
! Program generalized_gradient_test
!
! Tests the generalized gradient (GG) field machinery:
!   1) That the GG field and vector potential are self-consistent (field vs the gg_taylor
!      monomial map, and B = curl(A) for the vector potential).
!   2) That tracking a GG element through PTC (the "pancake" translation) agrees with
!      Bmad's runge_kutta integrator, and that the three PTC pancake settings
!      (pancake_symplectic, pancake_canonical, and plain) give the same result.
!-

program generalized_gradient_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele2, curl_ele
type (gg_taylor_struct) emt(3)
type (coord_struct) orb, curl_orb, start_orb, end_rk, end_symp, end_canon, end_plain
type (em_field_struct) field

real(rp) em_field(3), err, bcurl(3), b0(3), a0(3), s0, curl_s
integer i, j

!

open (1, file = 'output.now')

call bmad_parser ('generalized_gradient_test.bmad', lat)

!--------------------------------------------------------------------
! Part 1: General GG field / potential consistency (no PTC).

! Straight-frame element (ele 1): em_field_calc field and the gg_taylor monomial map are
! both built from the coefficient table and must agree.

ele => lat%ele(1)
orb%vec = 0
call init_coord (orb, orb, ele, downstream_end$)

call gen_grad_at_s_to_gg_taylor (ele, ele%gen_gradients(1), 0.05_rp, emt)
err = 0

do i = -1, 2
do j = -1, 2
  orb%vec(1) = 0.01_rp * i
  orb%vec(3) = 0.01_rp * j
  call em_field_calc (ele, lat%param, 0.05_rp, orb, .true., field)   ! same s as the gg_taylor above
  call evaluate_gg_taylor ([orb%vec(1), orb%vec(3)], emt, em_field)
  err = err + sum(abs(field%b - em_field))
enddo
enddo

write (1, '(a, es20.12)') '"GG-Field-vs-Taylor" ABS 1e-12', err

! Curved-frame element (ele 2): field + vector potential, and a B = curl(A) check.

ele2 => lat%ele(2)
orb%vec = 0
call init_coord (orb, orb, ele2, downstream_end$)
orb%vec(1) = 0.011_rp
orb%vec(3) = -0.007_rp
s0 = 0.05_rp

call em_field_calc (ele2, lat%param, s0, orb, .true., field, calc_potential = .true.)
b0 = field%b
a0 = field%a

write (1, '(a, 3es20.12)') '"Curved-B" ABS 1e-10', b0
write (1, '(a, 3es20.12)') '"Curved-A" ABS 1e-10', a0

! B = curl(A) check. Done on the straight-frame element (ele 1), where the flat-space curl
! of A equals B. (In the curved frame the identity carries the 1 + g_ref*x metric factor.)

curl_ele => ele
curl_orb%vec = 0
call init_coord (curl_orb, curl_orb, curl_ele, downstream_end$)
curl_orb%vec(1) = 0.013_rp
curl_orb%vec(3) = -0.009_rp
curl_s = 0.05_rp

call em_field_calc (curl_ele, lat%param, curl_s, curl_orb, .true., field, calc_potential = .true.)
b0 = field%b

bcurl(1) = avec(3, 1) - avec(2, 2)   ! dAs/dy - dAy/ds
bcurl(2) = avec(1, 2) - avec(3, 0)   ! dAx/ds - dAs/dx
bcurl(3) = avec(2, 0) - avec(1, 1)   ! dAy/dx - dAx/dy

write (1, '(a, es20.12)') '"Curl-Err" ABS 1e-6', maxval(abs(bcurl - b0))

!--------------------------------------------------------------------
! Part 2: PTC pancake tracking.
! Track the straight-frame GG element with runge_kutta (reference) and with symp_lie_ptc
! for each of the three PTC pancake settings. The settings are read by ele_to_fibre when the
! fibre is (re)built, which track1_symp_lie_ptc does on every call, so changing ptc_com
! between track1 calls takes effect.

ele => lat%ele(1)
call init_coord (start_orb, lat%particle_start, ele, upstream_end$)

call track1 (start_orb, ele, lat%param, end_rk)

ele%tracking_method  = symp_lie_ptc$
ele%mat6_calc_method = symp_lie_ptc$

ptc_com%pancake_symplectic = .true.        ! symplectic (canonical is forced True)
ptc_com%pancake_canonical  = .false.
call track1 (start_orb, ele, lat%param, end_symp)

ptc_com%pancake_symplectic = .false.       ! non-symplectic, canonical coordinates
ptc_com%pancake_canonical  = .true.
call track1 (start_orb, ele, lat%param, end_canon)

ptc_com%pancake_symplectic = .false.       ! non-symplectic, non-canonical
ptc_com%pancake_canonical  = .false.
call track1 (start_orb, ele, lat%param, end_plain)

! Each pancake setting must reproduce the runge_kutta result, and hence each other.

write (1, '(a, es16.6)') '"PTC-symp-vs-RK"     ABS 1e-5', maxval(abs(end_symp%vec  - end_rk%vec))
write (1, '(a, es16.6)') '"PTC-canon-vs-RK"    ABS 1e-5', maxval(abs(end_canon%vec - end_rk%vec))
write (1, '(a, es16.6)') '"PTC-plain-vs-RK"    ABS 1e-5', maxval(abs(end_plain%vec - end_rk%vec))
write (1, '(a, es16.6)') '"PTC-symp-vs-canon"  ABS 1e-5', maxval(abs(end_symp%vec  - end_canon%vec))
write (1, '(a, es16.6)') '"PTC-canon-vs-plain" ABS 1e-8', maxval(abs(end_canon%vec - end_plain%vec))

close (1)

contains

! Central-difference partial dA_icomp/du_idir at the curl test point.
! icomp = 1..3 -> Ax,Ay,As ; idir = 0,1,2 -> perturb x,y,s.

function avec (icomp, idir) result (deriv)
real(rp) deriv, dd, ap(3), am(3)
integer icomp, idir
type (coord_struct) o
type (em_field_struct) f
real(rp) s

dd = 1e-6_rp
o = curl_orb
s = curl_s
select case (idir)
case (0); o%vec(1) = curl_orb%vec(1) + dd
case (1); o%vec(3) = curl_orb%vec(3) + dd
case (2); s = curl_s + dd
end select
call em_field_calc (curl_ele, lat%param, s, o, .true., f, calc_potential = .true.)
ap = f%a

o = curl_orb
s = curl_s
select case (idir)
case (0); o%vec(1) = curl_orb%vec(1) - dd
case (1); o%vec(3) = curl_orb%vec(3) - dd
case (2); s = curl_s - dd
end select
call em_field_calc (curl_ele, lat%param, s, o, .true., f, calc_potential = .true.)
am = f%a

deriv = (ap(icomp) - am(icomp)) / (2 * dd)
end function avec

end program
