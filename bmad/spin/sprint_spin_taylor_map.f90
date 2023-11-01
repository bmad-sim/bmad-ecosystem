!+
! Subroutine sprint_spin_taylor_map (ele, start_orbit)
!
! Routine to calculate the spin Taylor map for a lattice element using the sprint formalism.
!
! Input:
!   ele             -- ele_struct: Element to form map for.
!   start_orbit(6)  -- real(rp), optional: Reference orbit for the map. Default is zero orbit. 
!
! Output:
!   ele             -- ele_struct: Element with map.
!     %spin_taylor(0:3)   -- Taylor map.
!-

subroutine sprint_spin_taylor_map (ele, start_orbit)

use bmad, dummy => sprint_spin_taylor_map

implicit none

type (ele_struct) ele
type (coord_struct) orb1, orb2, orb3, orb4
type (taylor_struct) spin_taylor(0:3)
type (fringe_field_info_struct) fringe_info
type (branch_struct), pointer :: branch
type (spin_orbit_map1_struct) map_start, map_ele, map_end, map_mis
type (track_struct) track

real(rp), optional :: start_orbit(6)
real(rp) gma, l, g, k1, k0, ks, kx, m, a, q, e1, e2
real(rp) cx, sx, cy, sy, omega, omegax, omegay, taux, tauy
real(rp) chi, zeta, psi, alpha, beta, sigma, xi, hkick, vkick
real(rp) d, c_d, s_d, e, c_e2, s_e2
real(rp) s, c_s, s_s, t, c_t2, s_t2

integer i, j, fringe_at, key
logical err_flag, spin_fringe

! Constants

m = mass_of(ele%ref_species)
e = ele%value(e_tot$)
gma = e/m
a = anomalous_moment_of(ele%ref_species)
q = charge_of(ele%ref_species)
l = ele%value(l$)

chi = 1 + a*gma
zeta = gma - 1
psi = gma**2 - 1

spin_fringe = is_true(ele%value(spin_fringe_on$))
fringe_at = nint(ele%value(fringe_at$))
branch => pointer_to_branch(ele)

do i = 0, 3
  call init_taylor_series(ele%spin_taylor(i), 0)
enddo

call mat_make_unit(map_start%orb_mat)
call mat_make_unit(map_end%orb_mat)
call mat_make_unit(map_ele%orb_mat)

! Transfer matrices.
! Currently, the sprint spin transport is with respect to the zero orbit

call init_coord(orb1, vec0$, ele, upstream_end$, ele%ref_species)

err_flag = .false.
call init_fringe_info (fringe_info, ele)
fringe_info%particle_at = first_track_edge$
orb2 = orb1
call apply_element_edge_kick(orb2, fringe_info, ele, branch%param, .false., map_start%orb_mat, .true.)
map_start%vec0 = orb2%vec - matmul(map_start%orb_mat, orb1%vec)

i = ele%value(fringe_at$)
ele%value(fringe_at$) = no_end$
orb3 = orb2
call track1_bmad(orb3, ele, branch%param, err_flag, track, map_ele%orb_mat, .true.)
ele%value(fringe_at$) = i
map_ele%vec0 = orb3%vec - matmul(map_ele%orb_mat, orb2%vec)

fringe_info%particle_at = second_track_edge$
orb4 = orb3
call apply_element_edge_kick(orb4, fringe_info, ele, branch%param, .false., map_end%orb_mat, .true.)
map_end%vec0 = orb4%vec - matmul(map_end%orb_mat, orb3%vec)

map_start%spin_q(0,0) = 1
map_end%spin_q(0,0) = 1

key = ele%key
if (key == quadrupole$ .and. abs(ele%value(k1$)) < 1e-20) key = pipe$

select case (key)

! Kicker
! Note: Kicks are put in at the end by offset_particle.

case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$, &
      kicker$, hkicker$, vkicker$, sextupole$, octupole$, thick_multipole$, rfcavity$)
  map_ele%spin_q(0,0) = 1

! Quadrupole

case (quadrupole$)
  k1 = ele%value(k1$)

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
  map_ele%spin_q(0,0) = 1.

  map_ele%spin_q(1,3) = -0.5_rp*chi*k1*sy
  map_ele%spin_q(1,4) = -0.5_rp*chi*k1*cy

  map_ele%spin_q(2,1) = -0.5_rp*chi*k1*sx
  map_ele%spin_q(2,2) = -0.5_rp*chi*k1*cx

! SBend

case (sbend$)
  e1 = ele%value(e1$)
  e2 = ele%value(e2$)
  k1 = ele%value(k1$)
  g = ele%value(g$)
  k0 = g + ele%value(dg$)
  kx = k1+g*k0

  d = k0*l
  c_d = cos(d)
  s_d = sin(d)

  e = a*k0*l*gma
  c_e2 = cos(0.5_rp*e)
  s_e2 = sin(0.5_rp*e)

  if (spin_fringe) then
    map_start%spin_q(1,3) = 0.5_rp*(1+a)*k0*sin(e1)
    map_start%spin_q(2,1) = 0.5_rp*chi*k0*tan(e1)
    map_start%spin_q(3,3) = -0.5_rp*(1+a)*k0*cos(e1)
    map_end%spin_q(1,3) = 0.5_rp*(1+a)*k0*sin(e2)
    map_end%spin_q(2,1) = 0.5_rp*chi*k0*tan(e2)
    map_end%spin_q(3,3) = 0.5_rp*(1+a)*k0*cos(e2)
  endif

  if (k1 == 0) then
    map_ele%spin_q(0,0) = c_e2
    map_ele%spin_q(0,1) = -0.5_rp*g*chi*s_d*s_e2
    map_ele%spin_q(0,2) = 0.5_rp*chi*(c_d-1)*s_e2
    map_ele%spin_q(0,6) = (1/(2*gma))*(gma*chi*s_d-a*psi*d)*s_e2

    map_ele%spin_q(2,0) = -s_e2
    map_ele%spin_q(2,1) = -0.5_rp*g*chi*s_d*c_e2
    map_ele%spin_q(2,2) = 0.5_rp*chi*(c_d-1)*c_e2
    map_ele%spin_q(2,6) = (1/(2*gma))*(gma*chi*s_d-a*psi*d)*c_e2

    map_ele%spin_q(3,4) = (1/gma)*zeta*s_e2

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

    map_ele%spin_q(0,0) = c_e2
    map_ele%spin_q(0,1) = -(1/(2*omegax))*kx*chi*sx*s_e2
    map_ele%spin_q(0,2) = (1/(2*omegax**2))*kx*chi*taux*(1-cx)*s_e2
    map_ele%spin_q(0,6) = -0.5_rp*g*((a*l*psi/gma)-(chi*sx/omegax))*s_e2

    map_ele%spin_q(1,3) = -(1/alpha)*(beta*(1+cy)*s_e2 + tauy*sigma*sy*c_e2)
    map_ele%spin_q(1,4) = -(1/(omegay*alpha))*(xi*(-1+cy)*c_e2 + beta*sy*s_e2)

    map_ele%spin_q(2,0) = -s_e2
    map_ele%spin_q(2,1) = -(1/(2*omegax))*kx*chi*sx*c_e2
    map_ele%spin_q(2,2) = (1/(2*omegax**2))*kx*chi*taux*(1-cx)*c_e2
    map_ele%spin_q(2,6) = -0.5_rp*g*((a*l*psi/gma)-(chi*sx/omegax))*c_e2

    map_ele%spin_q(3,3) = -(1/alpha)*(beta*(-1+cy)*c_e2 - tauy*sigma*sy*s_e2)
    map_ele%spin_q(3,4) = (1/(omegay*alpha))*(xi*(1 + cy)*s_e2 - beta*c_e2*sy)
  endif

! Solenoid

case (solenoid$)
  ks = ele%value(ks$)

  s = a*ks*l
  c_s = cos(s)
  s_s = sin(s)

  t = (1+a)*ks*l
  c_t2 = cos(t/2.)
  s_t2 = sin(t/2.)

  if (spin_fringe) then
    map_start%spin_q(1,1) = ks*chi/4.
    map_start%spin_q(2,3) = ks*chi/4.

    map_end%spin_q(1,1) = -ks*chi/4.
    map_end%spin_q(2,3) = -ks*chi/4.
  endif

  map_ele%spin_q(0, 0) = c_t2
  map_ele%spin_q(0,6) = 0.5_rp*t*s_t2

  map_ele%spin_q(1, 1) = 0.25_rp*ks*zeta*((1-c_s)*c_t2 - s_s * s_t2)
  map_ele%spin_q(1, 2) = 0.5_rp*zeta*((1-c_s)*s_t2 + s_s*c_t2)
  map_ele%spin_q(1, 3) = 0.25_rp*ks*zeta*((1-c_s)*s_t2 + s_s * c_t2)
  map_ele%spin_q(1, 4) = 0.5_rp*zeta*((-1+c_s)*c_t2 + s_s*s_t2)
  map_ele%spin_q(2, 1) = 0.25_rp*ks*zeta*((-1+c_s)*s_t2 - s_s * c_t2)
  map_ele%spin_q(2, 2) = 0.5_rp*zeta*((1-c_s)*c_t2 - s_s*s_t2)
  map_ele%spin_q(2, 3) = 0.25_rp*ks*zeta*((1-c_s)*c_t2 - s_s * s_t2)
  map_ele%spin_q(2, 4) = 0.5_rp*zeta*((1-c_s)*s_t2 + s_s*c_t2)

  map_ele%spin_q(3,0) = -s_t2
  map_ele%spin_q(3,6) = 0.5_rp*t*c_t2

case default
  print *, 'HELP! I SHOULD NOT BE HERE!!!!'
  return

end select

! Concatenate quaternions and maps.

if (fringe_at == both_ends$ .or. fringe_at == entrance_end$) then
  map_ele = map_ele * map_start
endif
if (fringe_at == both_ends$ .or. fringe_at == exit_end$) then
  map_ele = map_end * map_ele
endif

! Add in misalignments and h/vkicks.
! Note: 1st order quaternion for h/vkicks is absent.

if (present(start_orbit)) then
  ele%spin_taylor_ref_orb_in = start_orbit
else
  ele%spin_taylor_ref_orb_in = 0
endif

if (ele_has_nonzero_offset(ele) .or. ele_has_nonzero_kick(ele)) then
  call init_coord (orb1, ele%spin_taylor_ref_orb_in, ele, upstream_end$)
  orb2 = orb1
  call mat_make_unit(map_mis%orb_mat)
  call offset_particle (ele, set$, orb2, mat6 = map_mis%orb_mat, make_matrix = .true., spin_qrot = map_mis%spin_q)
  map_mis%vec0 = orb2%vec - matmul(map_mis%orb_mat, orb1%vec)

  call shift_to_reference_orbit(map_ele, orb2%vec)
  map_ele = map_ele * map_mis

  orb2%vec = matmul(map_ele%orb_mat, orb1%vec) + map_ele%vec0
  orb3 = orb2
  call mat_make_unit(map_mis%orb_mat)
  call offset_particle (ele, unset$, orb3, mat6 = map_mis%orb_mat, make_matrix = .true., spin_qrot = map_mis%spin_q)
  map_mis%vec0 = orb3%vec - matmul(map_mis%orb_mat, orb2%vec)

  map_ele = map_mis * map_ele

else
  call shift_to_reference_orbit(map_ele, start_orbit)
endif

! Convert map%spin_q to ele%spin_taylor

do i = 0, 3
do j = 0, 6
  if (j == 0) then
    call add_taylor_term(ele%spin_taylor(i), map_ele%spin_q(i,j), [0,0,0,0,0,0])
  else
    call add_taylor_term(ele%spin_taylor(i), map_ele%spin_q(i,j), taylor_expn([j]))
  endif
enddo
enddo

!-----------------------------------------------------------
contains

subroutine shift_to_reference_orbit (map_ele, ref_orbit)

type (spin_orbit_map1_struct) map_ele
real(rp), optional :: ref_orbit(6)
real(rp) f_norm
integer j

!

if (.not. present(ref_orbit)) return
if (all(ref_orbit == 0)) return

do j = 1, 6
  map_ele%spin_q(:,0) = map_ele%spin_q(:,0) + ref_orbit(j) * map_ele%spin_q(:,j)
enddo

! Renormalize to make 0th order qaternion have unit length
f_norm = 1.0_rp / norm2(map_ele%spin_q(:,0))
map_ele%spin_q(:,:) = map_ele%spin_q(:,:)  * f_norm

! Now make first order quaternions perpendicular to the 0th order quaternions.
do j = 1, 6
  f_norm = dot_product(map_ele%spin_q(:,0), map_ele%spin_q(:,j))
  map_ele%spin_q(:,j) = map_ele%spin_q(:,j) - f_norm * map_ele%spin_q(:,0)
enddo

end subroutine shift_to_reference_orbit

end subroutine
