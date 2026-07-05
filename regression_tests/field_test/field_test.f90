program field_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele2, curl_ele
type (gg_taylor_struct) emt(3)
type (coord_struct) orb, curl_orb
type (em_field_struct) field

real(rp) em_field(3), err, bcurl(3), b0(3), a0(3), s0, curl_s
integer i, j

!

open (1, file = 'output.now')

call bmad_parser ('field_test.bmad', lat)

! Straight-frame gen_gradients (ele 1): the em_field_calc field and the gg_taylor
! monomial map are both built from the coefficient table and must agree.

ele => lat%ele(1)
orb%vec = 0
call init_coord (orb, orb, ele, downstream_end$)

call gen_grad1_to_gg_taylor(ele, ele%gen_gradients(1), 0, emt)
err = 0

do i = -1, 2
do j = -1, 2
  orb%vec(1) = 0.01_rp * i
  orb%vec(3) = 0.01_rp * j
  call em_field_calc (ele, lat%param, 0.0_rp, orb, .true., field)
  call evaluate_gg_taylor ([orb%vec(1), orb%vec(3)], emt, em_field)
  err = err + sum(abs(field%b-em_field))
enddo
enddo

write (1, '(a, es20.12)') '"GG-Taylor-vs-Field" ABS 1e-12', err

! Curved-frame gen_gradients (ele 2): field + vector potential, and a B = curl(A) check.

ele2 => lat%ele(2)
orb%vec = 0
call init_coord (orb, orb, ele2, downstream_end$)
orb%vec(1) = 0.011_rp
orb%vec(3) = -0.007_rp
s0 = 0.0_rp

call em_field_calc (ele2, lat%param, s0, orb, .true., field, calc_potential = .true.)
b0 = field%b
a0 = field%a

write (1, '(a, 3es20.12)') '"Curved-B"  ABS 1e-10', b0
write (1, '(a, 3es20.12)') '"Curved-A"  ABS 1e-10', a0

! B = curl(A) check. Done on the straight-frame element (ele 1), where the flat-space
! curl of A equals B. (In the curved frame the identity carries the 1 + g_ref*x metric
! factor, so a Cartesian finite-difference curl would not match B there.)

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

write (1, '(a, es20.12)') '"Curl-Err"  ABS 1e-6', maxval(abs(bcurl - b0))

!

write (1, '(a, es20.12)') '"Symp_Err" ABS 1e-10', mat_symp_error(lat%branch(1)%ele(1)%mat6)

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
