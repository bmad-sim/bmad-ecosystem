!+
! Subroutine ele_to_sprint_spin_taylor_map (ele)
!
! Routine to calculate the spin Taylor map for a lattice element using the sprint formalism.
!
! Input:
!   ele       -- ele_struct: Element to form map for.
!
! Output:
!   ele       -- ele_struct: Element with map.
!     %spin_taylor(:)   -- Taylor map.
!-

subroutine ele_to_sprint_spin_taylor_map (ele)

use bmad, dummy => ele_to_sprint_spin_taylor_map

implicit none

type (ele_struct) ele
type (coord_struct) orb_start, orb_end, orb_ele
type (taylor_struct) spin_taylor(0:3)
type (fringe_field_info_struct) fringe_info
type (branch_struct), pointer :: branch
type (spin_orbit_map1_struct) map_start, map_ele, map_end
type (track_struct) track

real(rp) gma, l, g, k1, k0, ks, kx, m, a, q
real(rp) cx, sx, cy, sy, omega, omegax, omegay, taux, tauy
real(rp) chi, zeta, psi, alpha, beta, sigma, xi
real(rp) d, cd, cdh, sd, sdh, e, ce, ceh, se, seh
real(rp) s, cs, csh, ss, ssh, t, ct, cth, st, sth
logical err_flag, spin_fringe
integer i, j, fringe_at

! Constants

m = mass_of(ele%ref_species)
e = ele%value(e_tot$)
gma = e/m
a = anomalous_moment_of(ele%ref_species)
q = charge_of(ele%ref_species)
l = ele%value(l$)
k1 = ele%value(k1$)
ks = ele%value(ks$)
g = ele%value(g$)
k0 = g + ele%value(dg$)
kx = k1+g*k0

chi = 1 + a*gma
zeta = gma - 1
psi = gma**2 - 1

spin_fringe = is_true(ele%value(spin_fringe_on$))
fringe_at = nint(ele%value(fringe_at$))

ele%spin_taylor_ref_orb_in = 0  ! Sprint ref is always the zero orbit
do i = 0, 3
  call init_taylor_series(ele%spin_taylor(i), 0)
enddo

call mat_make_unit(map_start%orb_mat)
call mat_make_unit(map_end%orb_mat)
call mat_make_unit(map_ele%orb_mat)

map_start%vec0 = ele%spin_taylor_ref_orb_in
map_end%vec0 = ele%spin_taylor_ref_orb_in
map_ele%vec0 = ele%spin_taylor_ref_orb_in

call init_coord(orb_start, ele%spin_taylor_ref_orb_in, ele, upstream_end$, ele%ref_species)
call init_coord(orb_end, ele%spin_taylor_ref_orb_in, ele, downstream_end$, ele%ref_species)
call init_coord(orb_ele, ele%spin_taylor_ref_orb_in, ele, upstream_end$, ele%ref_species)


! call init_coord (orb, ele%spin_taylor_ref_orb_in, ele, upstream_end$, ele%ref_species)
branch => pointer_to_branch(ele)

!update transfer matrices

err_flag = .false.
call init_fringe_info (fringe_info, ele)
i = ele%value(fringe_at$)
fringe_info%particle_at = first_track_edge$
call apply_element_edge_kick(orb_start, fringe_info, ele, branch%param, .false., map_start%orb_mat, .true.)
ele%value(fringe_at$) = no_end$
call track1_bmad(orb_ele, ele, branch%param, orb_end, err_flag, track, map_ele%orb_mat, .true.)
ele%value(fringe_at$) = i
fringe_info%particle_at = second_track_edge$
call apply_element_edge_kick(orb_end, fringe_info, ele, branch%param, .false., map_end%orb_mat, .true.)

map_start%spin_q(0,0) = 1
map_end%spin_q(0,0) = 1

select case (ele%key)

! Drift

case (drift$)
  call add_taylor_term(ele%spin_taylor(0), 1.0_rp, [0,0,0,0,0,0])

! Kicker

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$)
  call add_taylor_term(ele%spin_taylor(0), 1.0_rp, [0,0,0,0,0,0])

! Quadrupole

case (quadrupole$)

  if (k1 > 0) then
    omega = sqrt(k1)
    sx = sin(l*omega)/omega
    cx = (1-cos(l*omega))/omega**2
    sy = sinh(l*omega)/omega
    cy = (-1+cosh(l*omega))/omega**2
  else
    omega = sqrt(-k1)
    sx = sinh(l*omega)/omega
    cx = (-1+cosh(l*omega))/omega**2
    sy = sin(l*omega)/omega
    cy = (1-cos(l*omega))/omega**2
  endif
  !Scalar terms
  map_ele%spin_q(0,0) = 1.
  !q_x terms
  map_ele%spin_q(1,3) = 0.5_rp*chi*k1*sy
  map_ele%spin_q(1,4) = 0.5_rp*chi*k1*cy
  !q_y terms
  map_ele%spin_q(2,1) = 0.5_rp*chi*k1*sx
  map_ele%spin_q(2,2) = 0.5_rp*chi*k1*cx

! SBend

case (sbend$)

  d = k0*l
  cd = cos(d)
  cdh = cos(0.5_rp*d)
  sd = sin(d)
  sdh = sin(0.5_rp*d)

  e = a*k0*l*gma
  ce = cos(e)
  ceh = cos(0.5_rp*e)
  se = sin(e)
  seh = sin(0.5_rp*e)

  map_start%spin_q(3,3) = -(1+a)*k0/2.
  map_end%spin_q(3,3) = (1+a)*k0/2.
  !constants
  if (k1 == 0) then
    map_ele%spin_q(0,0) = ceh
    map_ele%spin_q(0,1) = -0.5_rp*g*chi*sd*seh
    map_ele%spin_q(0,2) = 0.5_rp*chi*(cd-1)*seh
    map_ele%spin_q(0,6) = (1/(2*gma))*(gma*chi*sd-a*psi*d)*seh
    map_ele%spin_q(2,0) = -seh
    map_ele%spin_q(2,1) = -0.5_rp*g*chi*sd*ceh
    map_ele%spin_q(2,2) = 0.5_rp*chi*(cd-1)*ceh
    map_ele%spin_q(2,6) = (1/(2*gma))*(gma*chi*sd-a*psi*d)*ceh
    map_ele%spin_q(3,4) = (1/gma)*zeta*seh

  else
    if (kx > 0) then
      omegax = sqrt(kx)
      cx = cos(l*omegax)
      sx = sin(l*omegax)
      taux = -1
    else
      omegax = sqrt(-kx)
      cx = cosh(l*omegax)
      sx = sinh(l*omegax)
      taux = 1
    endif

    if (k1 > 0) then
      omegay = sqrt(k1)
      cy = cosh(l*omegay)
      sy = sinh(l*omegay)
      tauy = 1
    else
      omegay = sqrt(-k1)
      cy = cos(l*omegay)
      sy = sin(l*omegay)
      tauy = -1
    endif

    alpha = 2*(a**2*g**2*gma**2+k1)
    beta = a*g*k1*(gma*chi-zeta)
    sigma = (k1+a*k1*gma+a**2*g**2*zeta*gma)*omegay
    xi = (k1*chi+a**2*g**2*zeta*gma)*omegay
    map_ele%spin_q(0,0) = ceh
    map_ele%spin_q(0,1) = -(1/(2*omegax))*kx*chi*sx*seh
    map_ele%spin_q(0,2) = (1/(2*omegax**2))*kx*chi*taux*(1-cx)*seh
    map_ele%spin_q(0,6) = -0.5_rp*g*((a*l*psi/gma)-(chi*sx/omegax))*seh
    map_ele%spin_q(1,3) = -(1/alpha)*(beta*(1+cy)*seh + tauy*sigma*sy*ceh)
    map_ele%spin_q(1,4) = -(1/(omegay*alpha))*(xi*(-1+cy)*ceh + beta*sy*seh)
    map_ele%spin_q(2,0) = -seh
    map_ele%spin_q(2,1) = -(1/(2*omegax))*kx*chi*sx*ceh
    map_ele%spin_q(2,2) = (1/(2*omegax**2))*kx*chi*taux*(1-cx)*ceh
    map_ele%spin_q(2,6) = -0.5_rp*g*((a*l*psi/gma)-(chi*sx/omegax))*ceh
    map_ele%spin_q(3,3) = -(1/alpha)*(beta*(-1+cy)*ceh - tauy*sigma*sy*seh)
    map_ele%spin_q(3,4) = (1/(omegay*alpha))*(xi*(1 + cy)*seh - beta*ceh*sy)
  endif

! Solenoid

case (solenoid$)

  s = a*ks*l
  cs = cos(s)
  csh = cos(s/2.)
  ss = sin(s)
  ssh = sin(s/2.)

  t = (1+a)*ks*l
  ct = cos(t)
  cth = cos(t/2.)
  st = sin(t)
  sth = sin(t/2.)

  map_start%spin_q(1,1) = ks*chi/4.
  map_start%spin_q(2,3) = ks*chi/4.

  map_end%spin_q(1,1) = -ks*chi/4.
  map_end%spin_q(2,3) = -ks*chi/4.

  map_ele%spin_q(0, 0) =  cth
  map_ele%spin_q(0,6) = 0.5_rp*t*sth

  map_ele%spin_q(1, 1) = 0.25_rp*ks*zeta*((1-cs)*cth - ss * sth)
  map_ele%spin_q(1, 2) = 0.5_rp*zeta*((1-cs)*sth + ss*cth)
  map_ele%spin_q(1, 3) = 0.25_rp*ks*zeta*((1-cs)*sth + ss * cth)
  map_ele%spin_q(1, 4) = 0.5_rp*zeta*((-1+cs)*cth + ss*sth)
  map_ele%spin_q(2, 1) = 0.25_rp*ks*zeta*((-1+cs)*sth - ss * cth)
  map_ele%spin_q(2, 2) = 0.5_rp*zeta*((1-cs)*cth - ss*sth)
  map_ele%spin_q(2, 3) = 0.25_rp*ks*zeta*((1-cs)*cth - ss * sth)
  map_ele%spin_q(2, 4) = 0.5_rp*zeta*((1-cs)*sth + ss*cth)

  map_ele%spin_q(3,0) = -sth
  map_ele%spin_q(3,6) = 0.5_rp*t*cth

case default
  print *, 'HELP! I SHOULD NOT BE HERE!!!!'
  return

end select

!

if (spin_fringe) then
  if (fringe_at == both_ends$ .or. fringe_at == entrance_end$) then
    map_ele = map_ele * map_start
  endif
  if (fringe_at == both_ends$ .or. fringe_at == exit_end$) then
    map_ele = map_end * map_ele
  endif
endif

do i = 0,3
  do j = 0,6
    if (j == 0) then
      call add_taylor_term(ele%spin_taylor(i), map_ele%spin_q(i, j), [0,0,0,0,0,0])
  else
    call add_taylor_term(ele%spin_taylor(i), map_ele%spin_q(i, j), taylor_expn([j]))
  endif
enddo
enddo

end subroutine
