!+
! Module fringe_mod
!
! Module of helper routines for track1 routines
!-

module fringe_mod

use bmad_routine_interface

! Private routines for exact_bend_edge_kick
private ptc_rot_xz, ptc_wedger, ptc_fringe_dipoler

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, track_spin)
!
! Subroutine to track through the edge field of an sbend.
! This routine is called by apply_element_edge_kick only.
!
! Input:
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   orb         -- Coord_struct: Starting coords.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before fringe.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!   track_spin  -- logical, optional: If True then track the spin through the edge fields. Default: False.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including fringe.
!-

subroutine bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, track_spin)

implicit none

type (ele_struct) ele
type (coord_struct) orb, ave_orb
type (lat_param_struct) param
type (em_field_struct) field

real(rp), optional :: mat6(6,6)
real(rp) vec(6), e_ang, omega(3), x, y, tan_e_x

integer fringe_type, physical_end, particle_at

logical, optional :: make_matrix, track_spin
character(*), parameter :: r_name = 'bend_edge_kick'

!

if (.not. fringe_here(ele, orb, particle_at)) return
physical_end = physical_ele_end (particle_at, orb, ele%orientation)

! Higher order fringes. 

if (particle_at == first_track_edge$) then
  call hard_multipole_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
endif

! Bend Fringe

fringe_type = nint(ele%value(fringe_type$))

select case (fringe_type)
case (full$)
  call exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (basic_bend$, hard_edge_only$)
  call hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (sad_full$)
  if (particle_at == first_track_edge$) then
    call sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
    call hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
  else
    call hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
    call sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
  endif

case (soft_edge_only$)
  call sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (linear_edge$)
  call linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (none$)

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN FRINGE_TYPE: \i0\ ', i_array = [fringe_type])
  if (global_com%exit_on_error) call err_exit
end select

! Higher order fringe

if (particle_at == second_track_edge$) then
  call hard_multipole_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
endif

! spin

if (logic_option(.false., track_spin)) then
  ave_orb = orb
  ave_orb%vec = (ave_orb%vec + orb%vec) / 2   ! Use average position
  field%E = 0
  if (physical_end == entrance_end$) then
    e_ang = ele%value(e1$)
  else
    e_ang = ele%value(e2$)
  endif
  x = ave_orb%vec(1);  y = ave_orb%vec(3)
  tan_e_x = tan(e_ang) * x
  field%B = (ele%value(b_field$) + ele%value(db_field$)) * [-sin(e_ang)*y, -tan_e_x, cos(e_ang)*y]
  if (ele%value(b1_gradient$) /= 0) field%B = field%B - ele%value(b1_gradient$) * tan_e_x * [x*y, x*x - y*y, 0.0_rp]
  if (ele%value(b2_gradient$) /= 0) field%B = field%B - ele%value(b2_gradient$) * tan_e_x * [3*x*x*y - y**3, x**3 - 3*x*y*y, 0.0_rp]
  if (physical_ele_end(particle_at, orb, ele%orientation) == downstream_end$) field%B(3) = -field%B(3)
  omega = spin_omega (field, ave_orb, ave_orb%direction * ele%orientation) * ave_orb%time_dir
  call rotate_spin (omega, orb%spin)
endif

end subroutine bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the edge field of an sbend.
! Apply only the first order kick, which is edge focusing.
!
! Input:
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   particle_at -- Integer: first_track_edge$, or second_track_edge$, 
!   orb         -- Coord_struct: Starting coords.
!   mat6        -- Real(rp), optional: Transfer matrix up to the edge.
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix including the edge.
!-

subroutine linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) e, g_tot, fint, hgap, ht_x, ht_y, cos_e, sin_e, tan_e, v0(6)
real(rp) c_dir
real(rp) ht2, hs2, sec_e
integer particle_at, element_end

logical, optional :: make_matrix

character(*), parameter :: r_name = 'linear_bend_edge_kick'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a dg.
! Apply only the first order kick. That is, only edge focusing.

c_dir = rel_tracking_charge_to_mass(orb, param%particle) * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb, ele%orientation)

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(dg$)) * c_dir
else
  g_tot = 0
endif

if (element_end == entrance_end$) then
  e = ele%value(e1$); fint = ele%value(fint$); hgap = ele%value(hgap$)
else
  e = ele%value(e2$); fint = ele%value(fintx$); hgap = ele%value(hgapx$)
endif

cos_e = cos(e); sin_e = sin(e); tan_e = sin_e / cos_e; sec_e = 1 / cos_e
ht_x = g_tot * tan_e
ht2 = g_tot * tan_e**2
hs2 = g_tot * sec_e**2

if (fint == 0) then
  ht_y = -ht_x
else
  ht_y = -g_tot * tan(e - 2 * fint * g_tot * hgap * (1 + sin_e**2) / cos_e)
endif

v0 = orb%vec
ht_x =  orb%time_dir * ht_x
ht_y =  orb%time_dir * ht_y

orb%vec(1) = v0(1)
orb%vec(2) = v0(2) + ht_x * v0(1)
orb%vec(3) = v0(3)
orb%vec(4) = v0(4) + ht_y * v0(3)
if (logic_option(.false., make_matrix)) then
  mat6(2,:) = mat6(2,:) + ht_x * mat6(1,:)
  mat6(4,:) = mat6(4,:) + ht_y * mat6(3,:)
end if

end subroutine linear_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the edge field of an sbend using a 2nd order map.
! Adapted from:
!   Hwang and S. Y. Lee, 
!   "Dipole Fringe Field Thin Map for Compact Synchrotrons",
!   Phys. Rev. ST Accel. Beams, 12, 122401, (2015).
! See the Bmad manual for details.
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   particle_at -- Integer: first_track_edge$, or second_track_edge$
!   mat6        -- Real(rp), optional: Transfer matrix up to the edge.
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix including the edge.
!-

subroutine hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) e, g_tot, fint_gap, gt, cos_e, sin_e, tan_e, sec_e, v(6), k1_tane
real(rp) gt2, gs2, c_dir, k1, kmat(6,6), e_factor, fg_factor
real(rp) dx, dpx, dy, dpy, dz, td
integer particle_at, element_end, fringe_type

logical, optional :: make_matrix
character(*), parameter :: r_name = 'hwang_bend_edge_kick'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a dg.

c_dir = rel_tracking_charge_to_mass(orb, param%particle) * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb, ele%orientation)
fringe_type = nint(ele%value(fringe_type$))
e_factor = 1 / (1 + orb%vec(6))
td = rp8(orb%time_dir)

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(dg$)) * c_dir
  k1 = ele%value(k1$)
else
  g_tot = 0
  k1 = 0
endif

if (element_end == entrance_end$) then
  e = ele%value(e1$); fint_gap = ele%value(fint$) * ele%value(hgap$)
else
  e = ele%value(e2$); fint_gap = ele%value(fintx$) * ele%value(hgapx$)
endif

if (fringe_type == hard_edge_only$ .or. fringe_type == sad_full$) fint_gap = 0

cos_e = cos(e); sin_e = sin(e); tan_e = sin_e / cos_e; sec_e = 1 / cos_e
gt = g_tot * tan_e
gt2 = g_tot * tan_e**2
gs2 = g_tot * sec_e**2
k1_tane = k1 * tan_e * c_dir
fg_factor = 2 * fint_gap * gs2 * g_tot * sec_e * (1 + sin_e**2)

v = orb%vec

if (entering_element(orb, particle_at)) then
  dx  = (-gt2 * v(1)**2 + gs2 * v(3)**2) * e_factor / 2
  dpx = (gt * g_tot * (1 + 2 * tan_e**2) * v(3)**2 / 2 + gt2 * (v(1) * v(2) - v(3) * v(4)) + k1_tane * (v(1)**2 - v(3)**2)) * e_factor
  dy  = gt2 * v(1) * v(3) * e_factor
  dpy = (fg_factor * v(3) - gt2 * v(1) * v(4) - (g_tot + gt2) * v(2) * v(3) - 2 * k1_tane * v(1) * v(3)) * e_factor
  dz = e_factor**2 * 0.5_rp * (v(3)**2 * fg_factor &
            + v(1)**3 * (4.0_rp * k1_tane - gt * gt2) / 6.0_rp + 0.5_rp * v(1)*v(3)**2 * (-4.0_rp * k1_tane + gt * gs2) &
            + (v(1)**2*v(2) - 2.0_rp * v(1)*v(3)*v(4)) * gt2 - v(2)*v(3)**2 * gs2)

  orb%vec(1) = v(1) + td * dx
  orb%vec(2) = v(2) + td * (dpx + gt * v(1))
  orb%vec(3) = v(3) + td * dy
  orb%vec(4) = v(4) + td * (dpy - gt * v(3))
  orb%vec(5) = v(5) + td * dz

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit (kmat)
    kmat(1,1) = 1 - gt2 * v(1) * e_factor
    kmat(1,3) = gs2 * v(3) * e_factor
    kmat(2,1) = gt + (gt2 * v(2) + 2 * k1_tane * v(1)) * e_factor
    kmat(2,2) = 1 + (gt2 * v(1)) * e_factor
    kmat(2,3) = (-gt2 * v(4) - 2 * k1_tane * v(3) + gt * g_tot * (1 + 2 * tan_e**2) * v(3)) * e_factor
    kmat(2,4) = (-gt2 * v(3)) * e_factor
    kmat(3,1) = gt2 * v(3) * e_factor
    kmat(3,3) = 1 +  gt2 * v(1) * e_factor
    kmat(4,1) = (-gt2 * v(4) - 2 * k1_tane * v(3)) * e_factor
    kmat(4,2) = -(g_tot + gt2) * v(3) * e_factor
    kmat(4,3) = -gt + (fg_factor - (g_tot + gt2) * v(2) - 2 * k1_tane * v(1)) * e_factor
    kmat(4,4) = 1 - gt2 * v(1) * e_factor
    kmat(5,6) = (fg_factor * v(3)**2 / 2 + (4*k1_tane - gt*gt2) * v(1)**3 / 12 + (-4*k1_tane + gt*gs2) * v(1) * v(3)**2 /4 + &
                  gt2 * (v(1)**2 * v(2) - 2 * v(1) * v(3) * v(4)) / 2 - (g_tot + gt2) * v(2) * v(3)**2 / 2) * e_factor**2
  end if

else
  dx  = (gt2 * v(1)**2 - gs2 * v(3)**2) * e_factor / 2
  dpx = (gt2 * (v(3) * v(4) - v(1) * v(2)) + k1_tane * (v(1)**2 - v(3)**2) - gt * gt2 * (v(1)**2 + v(3)**2) / 2) * e_factor
  dy  = -gt2 * v(1) * v(3) * e_factor
  dpy = (fg_factor * v(3) + gt2 * v(1) * v(4) + (g_tot + gt2) * v(2) * v(3) + (gt * gs2 - 2 * k1_tane) * v(1) * v(3)) * e_factor
  dz = e_factor**2 * 0.5_rp * (v(3)**2 * fg_factor &
            + v(1)**3 * (4.0_rp * k1_tane - gt * gt2) / 6.0_rp + 0.5_rp * v(1)*v(3)**2 * (-4.0_rp * k1_tane + gt * gs2) &
            - (v(1)**2*v(2) - 2.0_rp * v(1)*v(3)*v(4)) * gt2 + v(2)*v(3)**2 * gs2)

  orb%vec(1) = v(1) + td * dx
  orb%vec(2) = v(2) + td * (dpx + gt * v(1))
  orb%vec(3) = v(3) + td * dy
  orb%vec(4) = v(4) + td * (dpy - gt * v(3))
  orb%vec(5) = v(5) + td * dz

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit (kmat)
    kmat(1,1) = 1 + gt2 * v(1) * e_factor
    kmat(1,3) = -gs2 * v(3) * e_factor
    kmat(2,1) = gt + (-gt2 * v(2) + 2 * k1_tane * v(1) - gt * gt2 * v(1)) * e_factor
    kmat(2,2) = 1 - gt2 * v(1) * e_factor
    kmat(2,3) = (gt2 * v(4) - 2 * k1_tane * v(3) - gt * gt2 * v(3)) * e_factor
    kmat(2,4) = gt2 * v(3) * e_factor
    kmat(3,1) = -gt2 * v(3) * e_factor
    kmat(3,3) = 1 -  gt2 * v(1) * e_factor
    kmat(4,1) = (gt2 * v(4) + (gt * gs2 - 2 * k1_tane) * v(3)) * e_factor
    kmat(4,2) = (g_tot + gt2) * v(3) * e_factor
    kmat(4,3) = -gt + (fg_factor + (g_tot + gt2) * v(2) + (gt * gs2 - 2 * k1_tane) * v(1)) * e_factor
    kmat(4,4) = 1 + gt2 * v(1) * e_factor
    kmat(5,6) = (fg_factor * v(3)**2 / 2 + (4*k1_tane - gt*gt2) * v(1)**3 / 12 + (-4*k1_tane + gt*gs2) * v(1) * v(3)**2 /4 + &
                  gt2 * (-v(1)**2 * v(2) + 2 * v(1) * v(3) * v(4)) / 2 + (g_tot + gt2) * v(2) * v(3)**2 / 2) * e_factor**2
  end if
endif

!

if (logic_option(.false., make_matrix)) then
  kmat(1,6) = -dx  * e_factor
  kmat(2,6) = -dpx * e_factor
  kmat(3,6) = -dy  * e_factor
  kmat(4,6) = -dpy * e_factor
  ! The m(5,x) terms follow from the symplectic condition.
  kmat(5,1) = -kmat(2,6)*kmat(1,1) + kmat(1,6)*kmat(2,1) - kmat(4,6)*kmat(3,1) + kmat(3,6)*kmat(4,1)
  kmat(5,2) = -kmat(2,6)*kmat(1,2) + kmat(1,6)*kmat(2,2) - kmat(4,6)*kmat(3,2) + kmat(3,6)*kmat(4,2)
  kmat(5,3) = -kmat(2,6)*kmat(1,3) + kmat(1,6)*kmat(2,3) - kmat(4,6)*kmat(3,3) + kmat(3,6)*kmat(4,3)
  kmat(5,4) = -kmat(2,6)*kmat(1,4) + kmat(1,6)*kmat(2,4) - kmat(4,6)*kmat(3,4) + kmat(3,6)*kmat(4,4)
  if (td == -1) call mat_inverse(kmat, kmat)
  mat6 = matmul (kmat, mat6)
endif

end subroutine hwang_bend_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine sad_mult_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the hard edge bend fringe field for a bend or sad_mult element.
! Only the bend field is taken into account here. Higher order multipolse must be handled elsewhere.
!
! This routine assumes that the particle coordinates are with respect to the actual magnet face.
! Thus finite e1/e2 must be taken into account by other routines.
!
! SAD calls this the "linear" fringe even though it is nonlinear.
!
! Input:
!   ele         -- ele_struct: Element with fringe.
!   param       -- lat_param_struct: Tracking parameters.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!   orbit       -- coord_struct: Starting coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix up to the fringe.
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine sad_mult_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) g, px, y, y2, rel_p, p_zy, yg, kmat(6,6), c_dir, t0, ppx2, td

integer fringe_type, particle_at
integer i, i_max, element_end

logical, optional :: make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= full$ .and. fringe_type /= hard_edge_only$) return
if (.not. fringe_here(ele, orbit, particle_at)) return

! Extract params from ele.
! Currently this routine is never called with sbend elements.

if (ele%key == sbend$) then
  g = ele%value(g$) + ele%value(dg$)
  t0 = 0
else  ! sad_mult
  if (.not. associated(ele%a_pole)) return
  call multipole1_ab_to_kt(ele%a_pole(0), ele%b_pole(0), 0, g, t0)
  g = g / ele%value(l$)
endif

!

if (g == 0) return
td = rp8(orbit%time_dir)
c_dir = rel_tracking_charge_to_mass(orbit, param%particle) * ele%orientation * orbit%direction
g = g * c_dir
if (particle_at == second_track_edge$) then
  g = -g
endif

! Rotate

if (t0 /= 0) call tilt_coords (t0, orbit%vec, mat6, make_matrix)

px = orbit%vec(2)
y = orbit%vec(3)
y2 = y**2
rel_p = 1 + orbit%vec(6)
ppx2 = rel_p**2 - px**2
if (ppx2 < 0) then
  orbit%state = lost$
  return
endif
p_zy = sqrt(ppx2)
yg = y2 * g**2 / 12

orbit%vec(1) = orbit%vec(1) + td * g * y2 * (1 - yg) * rel_p**2 / (2 * p_zy**3)
orbit%vec(4) = orbit%vec(4) - td * g * px * y * (1 - 2 * yg) / p_zy
orbit%vec(5) = orbit%vec(5) - td * g * y2 * px * (1 - yg) * rel_p / (2 * p_zy**3)

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(kmat)
  kmat(1,2) = 3 * g * px * y2 * (1 - yg) * rel_p**2 / (2 * p_zy**5)
  kmat(1,3) =  g * (y - 2 * y * yg) * rel_p**2 / p_zy**3
  kmat(1,6) = -g * y2 * (1 - yg) * rel_p * (rel_p**2 + 2 * px**2) / (2 * p_zy**5)
  kmat(4,2) = -g * y * (1 - 2 * yg) * rel_p**2 / p_zy**3
  kmat(4,3) = -g * px * (1 - 6 * yg) / p_zy
  kmat(4,6) =  g * px * y * (1 - 2 * yg) * rel_p / p_zy**3
  kmat(5,2) = -g * y2 * (1 - yg) * (rel_p**2 + 2 * px**2) * rel_p / (2 * p_zy**5)
  kmat(5,3) = -g * px * (y - 2 * y * yg) * rel_p / p_zy**3
  kmat(5,6) =  g * px * y2 * (1 - 2 * yg) * (2 * rel_p**2 + px**2) / (2 * p_zy**5)
  if (orbit%time_dir == -1) call mat_inverse(kmat, kmat)
  mat6 = matmul(kmat, mat6)
endif

! Rotate

if (t0 /= 0) call tilt_coords (-t0, orbit%vec, mat6, make_matrix)

end subroutine sad_mult_hard_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine soft_quadrupole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
!
! Routine to add the SAD "linear" soft edge (for finite f1 or f2).
! This routine assumes that the particle orbit has been rotated to the element reference frame.
! This routine is called with sad_mult and quadrupole elements.
!
! Input:
!   ele           -- ele_struct: Element being tracked through
!   param         -- lat_param_struct: Tracking parameters.
!   particle_at   -- integer: first_track_edge$, or second_track_edge$.
!   orbit         -- coord_struct: Position before kick.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix up to the edge.
!   make_matrix   -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit         -- coord_struct: Position after kick.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix with edge kick added on.
!-

subroutine soft_quadrupole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: m_ele
type (coord_struct) orbit
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) k1_rel, x, y, px, py, charge_dir, kmat(6,6)
real(rp) f1, f2, ef1, vec(4), rel_p, vx, vy, tilt1, td

integer particle_at
integer fringe_type

logical, optional :: make_matrix

!

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= soft_edge_only$ .and. fringe_type /= full$) return
if (.not. fringe_here(ele, orbit, particle_at)) return

charge_dir = ele%orientation * orbit%direction
if (associated(ele%branch)) charge_dir = charge_dir * rel_tracking_charge_to_mass(orbit, param%particle)

rel_p = 1 + orbit%vec(6)
td = rp8(orbit%time_dir)

!

select case (ele%key)
case (quadrupole$)
  k1_rel = charge_dir * ele%value(k1$) / rel_p
  tilt1 = 0

case (sad_mult$)
  ! Slice slaves and super slaves have their associated multipoles stored in the lord
  if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
    m_ele => pointer_to_super_lord(ele)
  else
    m_ele => ele
  endif
  if (.not. associated(m_ele%a_pole)) return
  call multipole1_ab_to_kt (m_ele%a_pole(1), m_ele%b_pole(1), 1, k1_rel, tilt1)
  k1_rel = charge_dir * k1_rel / (m_ele%value(l$) * rel_p)

case default
  return
end select

!

f1 = k1_rel * ele%value(fq1$)
f2 = k1_rel * ele%value(fq2$)
if (f1 == 0 .and. f2 == 0) return
if (particle_at == second_track_edge$) f1 = -f1

if (orbit%time_dir == -1) then
  f2 = -f2
endif

ef1 = exp(f1)

! 

call tilt_coords(tilt1, orbit%vec, mat6, make_matrix)

vec = orbit%vec(1:4)
vx = vec(2) / rel_p;  vy = vec(4) / rel_p

orbit%vec(5) = orbit%vec(5) - (f1 * vec(1) + f2 * (1 + f1/2) * vx / ef1) * vx + &
                              (f1 * vec(3) + f2 * (1 - f1/2) * vy * ef1) * vy

orbit%vec(1:2) = [vec(1) * ef1 + vx * f2,  vec(2) / ef1]
orbit%vec(3:4) = [vec(3) / ef1 - vy * f2,  vec(4) * ef1]

!

if (logic_option(.false., make_matrix)) then
  kmat = 0
  kmat(1,1) = ef1
  kmat(1,2) = f2 / rel_p 
  kmat(1,6) = -vec(1) * f1 * ef1 / rel_p - 2 * vx * f2 / rel_p
  kmat(2,2) = 1 / ef1
  kmat(2,6) = vx * f1 / ef1 
  kmat(3,3) = 1 / ef1
  kmat(3,4) = -f2 / rel_p
  kmat(3,6) =  vec(3) * f1 / ef1 / rel_p + 2 * vy * f2 / rel_p
  kmat(4,4) = ef1
  kmat(4,6) = -vy * f1 * ef1
  kmat(5,1) = -f1 * vx
  kmat(5,2) = -(f1 * vec(1) + f2 * (2 + f1) * vx / ef1) / rel_p
  kmat(5,3) =  f1 * vy
  kmat(5,4) =  (f1 * vec(3) + f2 * (2 - f1) * vy * ef1) / rel_p
  kmat(5,5) = 1
  kmat(5,6) = 2 * f1 * vec(1) * vx / rel_p + &
              f2 * (1 + f1/2) * vx**2 * (3 - f1) / (ef1 * rel_p) + & 
              f2 * f1 * vx**2 / (2 * ef1 * rel_p) - &
              2 * f1 * vec(3) * vy / rel_p - &
              f2 * (1 - f1/2) * vy**2 * ef1 * (3 + f1) / rel_p + & 
              f2 * f1 * vy**2 * ef1 / (2 * rel_p)
  kmat(6,6) = 1
  mat6 = matmul(kmat, mat6)
endif

call tilt_coords(-tilt1, orbit%vec, mat6, make_matrix)

end subroutine soft_quadrupole_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine hard_multipole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the hard edge field of a multipole.
! The dipole component is ignored and only quadrupole and higher multipoles are included.
!
! This routine handles elements of type:
!   sad_mult, sbend, quadrupole, sextupole
!
! For sad_mult elements, ele%a_pole and ele%b_pole ae used for the multipole values.
! For the other elements, k1 or k2 is used and it is assumed that we are in the element
! frame of reference so tilt = 0.
! 
! Input:
!   ele         -- ele_struct: Element with fringe.
!   param       -- lat_param_struct: Tracking parameters.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!   orbit       -- coord_struct: Starting coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix up to the fringe.
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine hard_multipole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit

real(rp), optional :: mat6(6,6)
real(rp) rel_p, cn, x, y, px, py, denom, ddenom_dx, ddenom_dy, ddenom_dpz, kmat(6,6), ap, bp, charge_dir
real(rp) fx, dfx_dx, dfx_dy, d2fx_dxx, d2fx_dxy, d2fx_dyy
real(rp) fy, dfy_dx, dfy_dy, d2fy_dxx, d2fy_dxy, d2fy_dyy

complex(rp) poly, poly_n1, poly_n2, dpoly_dx, dpoly_dy, d2poly_dxx, d2poly_dxy, d2poly_dyy
complex(rp) xy, xny, dxny_dx, dxny_dy, cab

integer fringe_type, particle_at
integer n, n_max

logical, optional :: make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
select case (fringe_type)
case (hard_edge_only$)
  if (ele%key == sbend$) return  ! Uses hwang_bend_edge which has k1 hard edge.
case (full$)
case default
  return
end select

if (.not. fringe_here(ele, orbit, particle_at)) return

charge_dir = ele%orientation * orbit%direction
if (associated(ele%branch)) charge_dir = charge_dir * rel_tracking_charge_to_mass(orbit, param%particle)

!

select case (ele%key)
case (sad_mult$)
  if (.not. associated(ele%a_pole)) return
  do n = ubound(ele%a_pole, 1), 0, -1
    if (ele%a_pole(n) /= 0 .or. ele%b_pole(n) /= 0) exit
  enddo
  n_max = n
case (quadrupole$, sbend$)
  n_max = 1
case (sextupole$)
  n_max = 2
case default
  if (global_com%exit_on_error) call err_exit
end select

!

x = orbit%vec(1)
y = orbit%vec(3)
xy = cmplx(x, y, rp)

poly_n1 = 1
poly = xy

rel_p = 1 + orbit%vec(6)

fx = 0
dfx_dx = 0
dfx_dy = 0
d2fx_dxx = 0
d2fx_dxy = 0
d2fx_dyy = 0

fy = 0
dfy_dx = 0
dfy_dy = 0
d2fy_dxx = 0
d2fy_dxy = 0
d2fy_dyy = 0

do n = 1, n_max

  poly_n2 = poly_n1
  poly_n1 = poly
  poly = poly * xy

  select case (ele%key)
  case (sad_mult$)
    if (ele%a_pole(n) == 0 .and. ele%b_pole(n) == 0) cycle
    ap = ele%a_pole(n) / ele%value(l$)
    bp = ele%b_pole(n) / ele%value(l$)
  case (quadrupole$, sbend$)
    ap = 0
    bp = ele%value(k1$)
  case (sextupole$)
    ap = 0
    bp = ele%value(k2$)
    if (n /= 2) cycle
  end select

  dpoly_dx = (n+1) * poly_n1
  dpoly_dy = i_imag * dpoly_dx

  d2poly_dxx = n * (n+1) * poly_n2
  d2poly_dxy = i_imag * d2poly_dxx
  d2poly_dyy = -d2poly_dxx

  cab = charge_dir * cmplx(bp, ap, rp) / (4 * (n + 2) * rel_p)
  if (particle_at == first_track_edge$) cab = -cab

  cn = real(n+3, rp) / (n+1) 
  xny = cmplx(x, -cn * y, rp)
  dxny_dy = cmplx(0.0_rp, -cn, rp)

  fx = fx + real(cab * poly * xny)
  dfx_dx = dfx_dx + real(cab * (dpoly_dx * xny + poly))
  dfx_dy = dfx_dy + real(cab * (dpoly_dy * xny + poly * dxny_dy))
  d2fx_dxx = d2fx_dxx + real(cab * (d2poly_dxx * xny + 2 * dpoly_dx))
  d2fx_dxy = d2fx_dxy + real(cab * (d2poly_dxy * xny + dpoly_dx * dxny_dy + dpoly_dy))
  d2fx_dyy = d2fx_dyy + real(cab * (d2poly_dyy * xny + 2 * dpoly_dy * dxny_dy))

  xny = cmplx(y, cn * x, rp)
  dxny_dx = cmplx(0.0_rp, cn, rp)

  fy = fy + real(cab * poly * xny)
  dfy_dx = dfy_dx + real(cab * (dpoly_dx * xny + poly * dxny_dx))
  dfy_dy = dfy_dy + real(cab * (dpoly_dy * xny + poly))
  d2fy_dxx = d2fy_dxx + real(cab * (d2poly_dxx * xny + 2 * dpoly_dx * dxny_dx))
  d2fy_dxy = d2fy_dxy + real(cab * (d2poly_dxy * xny + dpoly_dx + dpoly_dy * dxny_dx))
  d2fy_dyy = d2fy_dyy + real(cab * (d2poly_dyy * xny + 2 * dpoly_dy))
enddo

px = orbit%vec(2)
py = orbit%vec(4)
denom = (1 - dfx_dx) * (1 - dfy_dy) - dfx_dy * dfy_dx
ddenom_dx = -d2fx_dxx - d2fy_dxy + d2fx_dxx * dfy_dy + dfx_dx * d2fy_dxy - d2fx_dxy * dfy_dx - dfx_dy * d2fy_dxx 
ddenom_dy = -d2fx_dxy - d2fy_dyy + d2fx_dxy * dfy_dy + dfx_dx * d2fy_dyy - d2fx_dyy * dfy_dx - dfx_dy * d2fy_dxy 
ddenom_dpz = (dfx_dx + dfy_dy - 2 * dfx_dx * dfy_dy + 2 * dfx_dy * dfy_dx) / rel_p

orbit%vec(1) = orbit%vec(1) - fx
orbit%vec(2) = px           + ((1.0_rp - dfy_dy - denom) * px + dfy_dx * py) / denom
orbit%vec(3) = orbit%vec(3) - fy
orbit%vec(4) = py           + (dfx_dy * px + (1.0_rp - dfx_dx - denom) * py) / denom
orbit%vec(5) = orbit%vec(5) + (orbit%vec(2) * fx + orbit%vec(4) * fy ) / rel_p

if (logic_option(.false., make_matrix)) then
  kmat = 0
  kmat(1,1) = 1 - dfx_dx
  kmat(1,3) = -dfx_dy
  kmat(1,6) = fx / rel_p
  kmat(2,1) = (-d2fy_dxy * px + d2fy_dxx * py) / denom - orbit%vec(2) * ddenom_dx / denom
  kmat(2,2) = (1 - dfy_dy) / denom
  kmat(2,3) = (-d2fy_dyy * px + d2fy_dxy * py) / denom - orbit%vec(2) * ddenom_dy / denom
  kmat(2,4) = dfy_dx / denom
  kmat(2,6) = (dfy_dy * px - dfy_dx * py) / (denom * rel_p) - orbit%vec(2) * ddenom_dpz / denom
  kmat(3,1) = -dfy_dx
  kmat(3,3) = 1 - dfy_dy
  kmat(3,6) = fy / rel_p
  kmat(4,1) = (d2fx_dxy * px - d2fx_dxx * py) / denom - orbit%vec(4) * ddenom_dx / denom
  kmat(4,2) = dfx_dy / denom
  kmat(4,3) = (d2fx_dyy * px - d2fx_dxy * py) / denom - orbit%vec(4) * ddenom_dy / denom
  kmat(4,4) = (1 - dfx_dx) / denom
  kmat(4,6) = (-dfx_dy * px + dfx_dx * py) / (denom * rel_p) - orbit%vec(4) * ddenom_dpz / denom
  kmat(5,1) = (kmat(2,1) * fx + orbit%vec(2) * dfx_dx + kmat(4,1) * fy + orbit%vec(4) *dfy_dx) / rel_p
  kmat(5,2) = (kmat(2,2) * fx + kmat(4,2) * fy) / rel_p
  kmat(5,3) = (kmat(2,3) * fx + orbit%vec(2) * dfx_dy + kmat(4,3) * fy + orbit%vec(4) *dfy_dy) / rel_p
  kmat(5,4) = (kmat(2,4) * fx + kmat(4,4) * fy) / rel_p
  kmat(5,5) = 1
  kmat(5,6) = (kmat(2,6) * fx + kmat(4,6) * fy) / rel_p - 2 * (orbit%vec(2) * fx + orbit%vec(4) * fy) / rel_p**2
  kmat(6,6) = 1
  !!! if (orbit%time_dir == -1) call mat_inverse(kmat, kmat)
  mat6 = matmul (kmat, mat6)
endif

end subroutine hard_multipole_edge_kick 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the ("linear") bend soft edge field of an sbend or sad_mult.
!
! Input:
!   ele         -- ele_struct: SBend or sad_mult element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   orb         -- Coord_struct: Starting coords.
!   mat6(6,6)   -- real(rp), optional: Starting matrix 
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!   k0l         -- real(rp), optional: Used with sad_mult. 
!                     If present, use this instead of ele%a_pole/%b_pole.
!   t0          -- real(rp), optional: Used with sad_mult. 
!                     If present, use this instead of ele%a_pole/%b_pole.
!                     Must be present if k0l is.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix after fringe field
!-

subroutine sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) :: kmat(6,6), tilt
real(rp) :: f1, el_p, g, ct, c1, c2, c3, y, px, rel_p, sin_e
real(rp) k0, sk0, c1_k, c2_k, c3_k, c1_sk, c2_sk, c3_sk 

integer :: particle_at, c_dir, element_end, fringe_type

logical, optional :: make_matrix

character(*), parameter :: r_name = 'sad_soft_bend_edge_kick'

!

fringe_type = nint(ele%value(fringe_type$))

if (ele%key == sbend$) then
  if (fringe_type /= soft_edge_only$ .and. fringe_type /= sad_full$) return
else  ! Must be sad_mult
  if (fringe_type /= soft_edge_only$ .and. fringe_type /=  full$) return
endif

! 

element_end = physical_ele_end(particle_at, orb, ele%orientation)

select case (ele%key)
case (sbend$)
  if (element_end == entrance_end$) then
    sin_e = sin(ele%value(e1$))
    f1 = 12 * ele%value(fint$) * ele%value(hgap$) 
  else
    sin_e = sin(ele%value(e2$))
    f1 = 12 * ele%value(fintx$) * ele%value(hgapx$)
  endif

  if (f1 == 0) return

  g = ele%value(g$) + ele%value(dg$)
  tilt = 0

case (sad_mult$)
  if (element_end == entrance_end$) then
    f1 = ele%value(fb1$)
  else
    f1 = ele%value(fb2$)
  endif

  if (f1 == 0) return

  call multipole1_ab_to_kt(ele%a_pole(0), ele%b_pole(0), 0, g, tilt)
  g = g / ele%value(l$)
end select

!

if (g == 0) return
c_dir = rel_tracking_charge_to_mass(orb, param%particle) * ele%orientation * orb%direction * orb%time_dir
g = g * c_dir

if (particle_at == second_track_edge$) then
  sin_e = -sin_e
  g = -g
endif

if (tilt /= 0) call tilt_coords(tilt, orb%vec, mat6, make_matrix)

px = orb%vec(2) ! + sin_e ???
y  = orb%vec(3)
rel_p = 1 + orb%vec(6)

c1 = f1**2 * g / (24 * rel_p)  ! * px
c2 = f1 * g**2 / (6 * rel_p)  ! * y^2
c3 = 2 * g**2 / (3 * f1 * rel_p)   ! * y^4

if (logic_option(.false., make_matrix)) then
  call mat_make_unit (kmat)
  kmat(1,6) =  c1 / rel_p
  kmat(4,3) =  c2 - 3 * c3 * y**2
  kmat(4,6) = (-c2 * y + c3 * y**3) / rel_p
  kmat(5,2) = c1 / rel_p
  kmat(5,3) = (c2 * y - c3 * y**3) / rel_p
  kmat(5,6) = -2*(c1 * px + c2 * y**2/2 - c3 * y**4/4) / (rel_p**2) ! Derivative of the orb%vec(5) line below

  mat6(1,:) = mat6(1,:) + kmat(1,6) * mat6(6,:)
  mat6(4,:) = mat6(4,:) + kmat(4,3) * mat6(3,:) + kmat(4,6) * mat6(6,:)
  mat6(5,:) = mat6(5,:) + kmat(5,2) * mat6(2,:) + kmat(5,3) * mat6(3,:) + kmat(5,6) * mat6(6,:)
endif

orb%vec(1) = orb%vec(1) + c1 * orb%vec(6)
orb%vec(4) = orb%vec(4) + c2 * y - c3 * y**3
orb%vec(5) = orb%vec(5) + (c1 * px + c2 * y**2/2 - c3 * y**4/4) / rel_p

if (tilt /= 0) call tilt_coords (-tilt, orb%vec, mat6, make_matrix)

end subroutine sad_soft_bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_wedger (a, g_tot, beta0, X, err_flag, mat6, make_matrix)
!
! Subroutine to track PTC coordinates through a wedge
!
! Adapted from forest/code/Sh_def_kind.f90 : WEDGER
!
! Input:
!   a           -- real(rp): wedge angle (rad)
!   g_tot       -- real(rp): reference bending radius
!   beta0       -- real(rp): reference relativistic beta
!   X(6)        -- real(rp): PTC phase space coordinates
!   make_matrix -- logical, optional: Propagate transfer matrix? Default is False.
!
! Output:
!   X(6)        -- real(rp): PTC phase space coordinates
!   err_flag    -- real(rp): Set true if there is a problem.
!   mat6        -- real(rp), optional: Transfer matrix.
!-
subroutine ptc_wedger (a, g_tot, beta0, X, err_flag, mat6, make_matrix)

implicit none

real(rp) :: X(6)
real(rp) :: a, beta0, g_tot
real(rp) :: Xn(6),pz,pzs,pt,b1
character(*), parameter :: r_name = 'ptc_wedger'
real(rp), optional :: mat6(6,6)

real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx4, dpt_dx5, fac, radix
real(rp) dpzs_dx1, dpzs_dx2, dpzs_dx4, dpzs_dx5, factor1, factor2

logical err_flag
logical, optional :: make_matrix

! No net field case...

err_flag = .true.
b1 = g_tot

if (b1==0) then
  call ptc_rot_xz(a, X, beta0, err_flag, mat6, make_matrix)
  return
endif

! Normal case

fac = 1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(2)**2-X(4)**2
if (fac < 0) return
radix = 1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(4)**2
if (radix < 1e-10) return
pz=sqrt(fac)
pt=sqrt(radix)

Xn(2)=X(2)*cos(a)+(pz-b1*X(1))*sin(a)

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(mat6)
  dpz_dx2 = -X(2)/pz
  dpz_dx4 = -X(4)/pz
  dpz_dx5 = (1/beta0+X(5))/pz
  dpt_dx4 = -X(4)/pt
  dpt_dx5 = (1/beta0+X(5))/pt
  mat6(2,1) = -b1*sin(a)
  mat6(2,2) = cos(a)+sin(a)*dpz_dx2
  mat6(2,4) = sin(a)*dpz_dx4
  mat6(2,5) = sin(a)*dpz_dx5
end if

radix = 1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-Xn(2)**2-X(4)**2
if (radix < 1e-10) return
pzs=sqrt(radix)

if (logic_option(.false., make_matrix)) then
  dpzs_dx1 = -Xn(2)*mat6(2,1)/pzs
  dpzs_dx2 = -Xn(2)*mat6(2,2)/pzs
  dpzs_dx4 = (-X(4)-Xn(2)*mat6(2,4))/pzs
  dpzs_dx5 = (1/beta0+X(5)-Xn(2)*mat6(2,5))/pzs
end if

Xn(1)=X(1)*cos(a)+(X(1)*X(2)*sin(2.0_rp*a)+sin(a)**2*(2.0_rp*X(1)*pz-b1*X(1)**2)) / (pzs+pz*cos(a)-X(2)*sin(a))

if (logic_option(.false., make_matrix)) then
  factor1 = cos(a)*pz+pzs-X(2)*sin(a)
  factor2 = (-b1*X(1)**2+2*X(1)*pz)*sin(a)**2+X(1)*X(2)*sin(2*a)
  mat6(1,1) = cos(a)+((-2*b1*X(1)+2*pz)*sin(a)**2+X(2)*sin(2*a))/factor1 - (factor2*dpzs_dx1)/factor1**2
  mat6(1,2) = (X(1)*sin(2*a)+2*X(1)*sin(a)**2*dpz_dx2)/factor1 &
              -(factor2*(-sin(a)+cos(a)*dpz_dx2+dpzs_dx2))/factor1**2
  mat6(1,4) = (2*X(1)*sin(a)**2*dpz_dx4)/factor1 - (factor2*(cos(a)*dpz_dx4+dpzs_dx4))/factor1**2
  mat6(1,5) = (2*X(1)*sin(a)**2*dpz_dx5)/factor1 - (factor2*(cos(a)*dpz_dx5+dpzs_dx5))/factor1**2
end if

Xn(3)=(a+asin(X(2)/pt)-asin(Xn(2)/pt))/b1

if (logic_option(.false., make_matrix)) then
  factor1 = sqrt(1-(Xn(2)/pt)**2)
  factor2 = sqrt(1-(X(2)/pt)**2)
  mat6(3,1) = - mat6(2,1)/(b1*factor1*pt)
  mat6(3,2) = (1/factor2-mat6(2,2)/factor1)/(b1*pt)
  mat6(3,4) = (-X(2)*dpt_dx4/(factor2*pt**2)-(-Xn(2)*dpt_dx4/pt**2+mat6(2,4)/pt)/factor1)/b1
  mat6(3,5) = (-X(2)*dpt_dx5/(factor2*pt**2)-(-Xn(2)*dpt_dx5/pt**2+mat6(2,5)/pt)/factor1)/b1
  mat6(6,1) = (1/beta0+X(5))*mat6(3,1)
  mat6(6,2) = (1/beta0+X(5))*mat6(3,2)
  mat6(6,4) = (1/beta0+X(5))*mat6(3,4)
  mat6(6,5) = Xn(3)+(1/beta0+X(5))*mat6(3,5)
  mat6(3,1) = X(4)*mat6(3,1)
  mat6(3,2) = X(4)*mat6(3,2)
  mat6(3,4) = Xn(3)+X(4)*mat6(3,4)
  mat6(3,5) = X(4)*mat6(3,5)
end if

Xn(6)=X(6)+Xn(3)*(1.0_rp/beta0+X(5))

Xn(3)=X(3)+X(4)*Xn(3)

X(1)=Xn(1)
X(2)=Xn(2)
X(3)=Xn(3)
X(6)=Xn(6)

err_flag = .false.

end subroutine ptc_wedger

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, err_flag, mat6, make_matrix)
!
! Subroutine to compute the exact hard edge fringe field of a bend.
! 
! Adapted from forest/code/Sh_def_kind.f90 : FRINGE_DIPOLER
!
!
! Input:
!   X(6)        -- real(rp): PTC phase space coordinates
!   beta0       -- real(rp): reference relativistic beta
!   g_tot       -- real(rp): reference bending radius
!   fint        -- real(rp): field integral for pole face
!   hgap        -- real(rp): gap height at pole face in meters. Only the product fint*hgap is used.
!   particle_at -- Integer: first_track_edge$, or second_track_edge$
!   make_matrix -- logical, optional: Make the matrix? Default is False. 
!
! Output:
!   X(6)        -- real(rp): PTC phase space coordinates
!   err_flag    -- real(rp): Set true if ther is a problem.
!   mat6        -- Real(rp), optional: Transfer matrix.
!-

subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, err_flag, mat6, make_matrix)

implicit none

integer  :: i
integer  :: particle_at

real(rp), optional :: mat6(6,6)
real(rp) :: X(6) !PTC phase space coordinates
real(rp) :: FINT, HGAP
real(rp) :: beta0, g_tot, fac
real(rp) :: pz, xp, yp, time_fac
real(rp) :: d(3,3), fi(3), fi0, fi0t, b, bb, co1, co2
real(rp) :: dpz_dx2, dpz_dx4, dpz_dx5, dtime_fac_dx5
real(rp) :: dxp_dx2, dxp_dx4, dxp_dx5, dyp_dx2, dyp_dx4, dyp_dx5
real(rp) :: d11_dx2, d11_dx4, d11_dx5, d21_dx2, d21_dx4, d21_dx5, d31_dx2, d31_dx4, d31_dx5
real(rp) :: d12_dx2, d12_dx4, d12_dx5, d22_dx2, d22_dx4, d22_dx5, d32_dx2, d32_dx4, d32_dx5
real(rp) :: d13_dx2, d13_dx4, d13_dx5, d23_dx2, d23_dx4, d23_dx5, d33_dx2, d33_dx4, d33_dx5
real(rp) :: dco1_dx2, dco1_dx4, dco1_dx5, dco2_dx2, dco2_dx4, dco2_dx5
real(rp) :: dfi0_dx2, dfi0_dx4, dfi0_dx5, dfi1_dx2, dfi1_dx4, dfi1_dx5
real(rp) :: dfi2_dx2, dfi2_dx4, dfi2_dx5, dfi3_dx2, dfi3_dx4, dfi3_dx5
real(rp) :: factor1, factor2 

character(*), parameter :: r_name = 'ptc_fringe_dipoler'

logical err_flag
logical, optional :: make_matrix

!

err_flag = .true.

if (particle_at == second_track_edge$) then
   B = -g_tot  !EL%CHARGE*BN(1)
else if (particle_at == first_track_edge$) then
   B = g_tot       
else
  call out_io (s_fatal$, r_name, 'INVALID PARTICLE_AT')
  if (global_com%exit_on_error) call err_exit
endif

fac = 1.0_rp+2.0_rp*x(5)/beta0+x(5)**2-x(2)**2-x(4)**2
if (fac == 0) return
pz=sqrt(fac)
time_fac=1.0_rp/beta0+x(5)
dtime_fac_dx5 = 1

xp=x(2)/pz
yp=x(4)/pz

d(1,1)=(1.0_rp+xp**2)/pz
d(2,1)=xp*yp/pz
d(3,1)=-xp
d(1,2)=xp*yp/pz
d(2,2)=(1.0_rp+yp**2)/pz
d(3,2)=-yp
d(1,3)=-time_fac*xp/pz**2
d(2,3)=-time_fac*yp/pz**2
d(3,3)= time_fac/pz

fi0= atan((xp/(1.0_rp+yp**2)))-b*fint*hgap*2.0_rp*( 1.0_rp + xp**2*(2.0_rp+yp**2) )*pz

co2=b/cos(fi0)**2
co1=co2/(1.0_rp+(xp/(1.0_rp+yp**2))**2 )

fi(1)=co1/(1.0_rp+yp**2)-co2*b*fint*hgap*2.0_rp*( 2.0_rp*xp*(2.0_rp+yp**2)*pz )
fi(2)=-co1*2.0_rp*xp*yp/(1.0_rp+yp**2)**2-co2*b*fint*hgap*2.0_rp*( 2.0_rp*xp**2*yp)*pz
fi(3)=-co2*b*fint*hgap*2.0_rp*( 1.0_rp + xp**2*(2.0_rp+yp**2) )

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(mat6)

  dpz_dx2 = -X(2)/pz
  dpz_dx4 = -X(4)/pz
  dpz_dx5 = (1/beta0+X(5))/pz

  dxp_dx2 = 1/pz-x(2)*dpz_dx2/pz**2
  dxp_dx4 = -x(2)*dpz_dx4/pz**2
  dxp_dx5 = -x(2)*dpz_dx5/pz**2

  dyp_dx2 = -x(4)*dpz_dx2/pz**2
  dyp_dx4 = 1/pz-x(4)*dpz_dx4/pz**2
  dyp_dx5 = -x(4)*dpz_dx5/pz**2

  d11_dx2 = 2*xp*dxp_dx2/pz-(1+xp**2)*dpz_dx2/pz**2
  d11_dx4 = 2*xp*dxp_dx4/pz-(1+xp**2)*dpz_dx4/pz**2
  d11_dx5 = 2*xp*dxp_dx5/pz-(1+xp**2)*dpz_dx5/pz**2

  d21_dx2 = xp*dyp_dx2/pz+yp*dxp_dx2/pz-xp*yp*dpz_dx2/pz**2
  d21_dx4 = xp*dyp_dx4/pz+yp*dxp_dx4/pz-xp*yp*dpz_dx4/pz**2
  d21_dx5 = xp*dyp_dx5/pz+yp*dxp_dx5/pz-xp*yp*dpz_dx5/pz**2

  d31_dx2 = -dxp_dx2
  d31_dx4 = -dxp_dx4
  d31_dx5 = -dxp_dx5

  d12_dx2 = d21_dx2
  d12_dx4 = d21_dx4
  d12_dx5 = d21_dx5

  d22_dx2 = 2*yp*dyp_dx2/pz-(1+yp**2)*dpz_dx2/pz**2
  d22_dx4 = 2*yp*dyp_dx4/pz-(1+yp**2)*dpz_dx4/pz**2
  d22_dx5 = 2*yp*dyp_dx5/pz-(1+yp**2)*dpz_dx5/pz**2

  d32_dx2 = -dyp_dx2
  d32_dx4 = -dyp_dx4
  d32_dx5 = -dyp_dx5

  d13_dx2 = 2*time_fac*xp*dpz_dx2/pz**3-time_fac*dxp_dx2/pz**2
  d13_dx4 = 2*time_fac*xp*dpz_dx4/pz**3-time_fac*dxp_dx4/pz**2
  d13_dx5 = 2*time_fac*xp*dpz_dx5/pz**3-time_fac*dxp_dx5/pz**2-xp*dtime_fac_dx5/pz**2

  d23_dx2 = 2*time_fac*yp*dpz_dx2/pz**3-time_fac*dyp_dx2/pz**2
  d23_dx4 = 2*time_fac*yp*dpz_dx4/pz**3-time_fac*dyp_dx4/pz**2
  d23_dx5 = 2*time_fac*yp*dpz_dx5/pz**3-time_fac*dyp_dx5/pz**2-yp*dtime_fac_dx5/pz**2

  d33_dx2 = -time_fac*dpz_dx2/pz**2
  d33_dx4 = -time_fac*dpz_dx4/pz**2
  d33_dx5 = -time_fac*dpz_dx5/pz**2-dtime_fac_dx5/pz

  factor1 = b*fint*hgap
  factor2 = 1+yp**2
  dfi0_dx2 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx2-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx2+2*xp**2*yp*dyp_dx2) &
             +(dxp_dx2/factor2-2*xp*yp*dyp_dx2/factor2**2)/(1+(xp/factor2)**2) 
  dfi0_dx4 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx4-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx4+2*xp**2*yp*dyp_dx4) &
             +(dxp_dx4/factor2-2*xp*yp*dyp_dx4/factor2**2)/(1+(xp/factor2)**2)
  dfi0_dx5 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx5-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx5+2*xp**2*yp*dyp_dx5) &
             +(dxp_dx5/factor2-2*xp*yp*dyp_dx5/factor2**2)/(1+(xp/factor2)**2)

  dco2_dx2 = 2*b*tan(fi0)*dfi0_dx2/cos(fi0)**2
  dco2_dx4 = 2*b*tan(fi0)*dfi0_dx4/cos(fi0)**2
  dco2_dx5 = 2*b*tan(fi0)*dfi0_dx5/cos(fi0)**2

  dco1_dx2 = dco2_dx2/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx2/factor2**2-4*xp**2*yp*dyp_dx2/factor2**3)/(1+(xp/factor2)**2)**2
  dco1_dx4 = dco2_dx4/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx4/factor2**2-4*xp**2*yp*dyp_dx4/factor2**3)/(1+(xp/factor2)**2)**2
  dco1_dx5 = dco2_dx5/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx5/factor2**2-4*xp**2*yp*dyp_dx5/factor2**3)/(1+(xp/factor2)**2)**2

  dfi1_dx2 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx2-4*factor1*co2*xp*(2+yp**2)*dpz_dx2-4*factor1*co2*pz*(2+yp**2)*dxp_dx2 &
             -8*factor1*co2*pz*xp*yp*dyp_dx2+dco1_dx2/factor2-2*co1*yp*dyp_dx2/factor2**2
  dfi1_dx4 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx4-4*factor1*co2*xp*(2+yp**2)*dpz_dx4-4*factor1*co2*pz*(2+yp**2)*dxp_dx4 &
             -8*factor1*co2*pz*xp*yp*dyp_dx4+dco1_dx4/factor2-2*co1*yp*dyp_dx4/factor2**2
  dfi1_dx5 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx5-4*factor1*co2*xp*(2+yp**2)*dpz_dx5-4*factor1*co2*pz*(2+yp**2)*dxp_dx5 &
             -8*factor1*co2*pz*xp*yp*dyp_dx5+dco1_dx5/factor2-2*co1*yp*dyp_dx5/factor2**2

  dfi2_dx2 = -4*factor1*pz*xp**2*yp*dco2_dx2-4*factor1*co2*xp**2*yp*dpz_dx2-8*factor1*co2*pz*xp*yp*dxp_dx2 &
             -4*factor1*co2*pz*xp**2*dyp_dx2-2*xp*yp*dco1_dx2/factor2**2-2*co1*yp*dxp_dx2/factor2**2 &
             +8*co1*xp*yp**2*dyp_dx2/factor2**3-2*co1*xp*dyp_dx2/factor2**2
  dfi2_dx4 = -4*factor1*pz*xp**2*yp*dco2_dx4-4*factor1*co2*xp**2*yp*dpz_dx4-8*factor1*co2*pz*xp*yp*dxp_dx4 &
             -4*factor1*co2*pz*xp**2*dyp_dx4-2*xp*yp*dco1_dx4/factor2**2-2*co1*yp*dxp_dx4/factor2**2 &
             +8*co1*xp*yp**2*dyp_dx4/factor2**3-2*co1*xp*dyp_dx4/factor2**2
  dfi2_dx5 = -4*factor1*pz*xp**2*yp*dco2_dx5-4*factor1*co2*xp**2*yp*dpz_dx5-8*factor1*co2*pz*xp*yp*dxp_dx5 &
             -4*factor1*co2*pz*xp**2*dyp_dx5-2*xp*yp*dco1_dx5/factor2**2-2*co1*yp*dxp_dx5/factor2**2 &
             +8*co1*xp*yp**2*dyp_dx5/factor2**3-2*co1*xp*dyp_dx5/factor2**2

  dfi3_dx2 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx2-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx2+2*xp**2*yp*dyp_dx2)
  dfi3_dx4 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx4-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx4+2*xp**2*yp*dyp_dx4)
  dfi3_dx5 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx5-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx5+2*xp**2*yp*dyp_dx5)

  dfi0_dx2 = b*dfi0_dx2/cos(fi0)**2
  dfi0_dx4 = b*dfi0_dx4/cos(fi0)**2
  dfi0_dx5 = b*dfi0_dx5/cos(fi0)**2
endif

fi0t=b*tan(fi0)
bb = sum(fi * d(:,2))

if (logic_option(.false., make_matrix)) then
  factor1 = sqrt(1-2*bb*x(3))
  factor2 = (1+sqrt(1-2*bb*x(3)))**2
  mat6(3,2) = 2*x(3)**2*(fi(1)*d12_dx2+fi(2)*d22_dx2+fi(3)*d32_dx2+d(1,2)*dfi1_dx2+d(2,2)*dfi2_dx2+d(3,2)*dfi3_dx2)/(factor1*factor2)
  mat6(3,3) = 2/(1+factor1)+2*x(3)*bb/(factor1*factor2)
  mat6(3,4) = 2*x(3)**2*(fi(1)*d12_dx4+fi(2)*d22_dx4+fi(3)*d32_dx4+d(1,2)*dfi1_dx4+d(2,2)*dfi2_dx4+d(3,2)*dfi3_dx4)/(factor1*factor2)
  mat6(3,5) = 2*x(3)**2*(fi(1)*d12_dx5+fi(2)*d22_dx5+fi(3)*d32_dx5+d(1,2)*dfi1_dx5+d(2,2)*dfi2_dx5+d(3,2)*dfi3_dx5)/(factor1*factor2)
end if

fac = 1.0_rp-2.0_rp*bb*x(3)
if (fac < 0) return
x(3)=2.0_rp*x(3)/(1.0_rp+ sqrt(fac))
x(4)=x(4)-fi0t*x(3)

if (logic_option(.false., make_matrix)) then
  mat6(4,2) = -fi0t*mat6(3,2)-x(3)*dfi0_dx2
  mat6(4,3) = -fi0t*mat6(3,3)
  mat6(4,4) = -fi0t*mat6(3,4)-x(3)*dfi0_dx4+1
  mat6(4,5) = -fi0t*mat6(3,5)-x(3)*dfi0_dx5
end if

bb = sum(fi * d(:,1))
x(1)=x(1)+0.5_rp*bb*x(3)**2

if (logic_option(.false., make_matrix)) then
  mat6(1,2) = x(3)*bb*mat6(3,2)+0.5*x(3)**2*(fi(1)*d11_dx2+fi(2)*d21_dx2+fi(3)*d31_dx2+d(1,1)*dfi1_dx2+d(2,1)*dfi2_dx2+d(3,1)*dfi3_dx2)
  mat6(1,3) = x(3)*bb*mat6(3,3)
  mat6(1,4) = x(3)*bb*mat6(3,4)+0.5*x(3)**2*(fi(1)*d11_dx4+fi(2)*d21_dx4+fi(3)*d31_dx4+d(1,1)*dfi1_dx4+d(2,1)*dfi2_dx4+d(3,1)*dfi3_dx4)
  mat6(1,5) = x(3)*bb*mat6(3,5)+0.5*x(3)**2*(fi(1)*d11_dx5+fi(2)*d21_dx5+fi(3)*d31_dx5+d(1,1)*dfi1_dx5+d(2,1)*dfi2_dx5+d(3,1)*dfi3_dx5)
end if

bb = sum(fi * d(:,3))
x(6)=x(6)-0.5_rp*bb*x(3)**2

if (logic_option(.false., make_matrix)) then
  mat6(6,2) = -x(3)*bb*mat6(3,2)-0.5*x(3)**2*(fi(1)*d13_dx2+fi(2)*d23_dx2+fi(3)*d33_dx2+d(1,3)*dfi1_dx2+d(2,3)*dfi2_dx2+d(3,3)*dfi3_dx2)
  mat6(6,3) = -x(3)*bb*mat6(3,3)
  mat6(6,4) = -x(3)*bb*mat6(3,4)-0.5*x(3)**2*(fi(1)*d13_dx4+fi(2)*d23_dx4+fi(3)*d33_dx4+d(1,3)*dfi1_dx4+d(2,3)*dfi2_dx4+d(3,3)*dfi3_dx4)
  mat6(6,5) = -x(3)*bb*mat6(3,5)-0.5*x(3)**2*(fi(1)*d13_dx5+fi(2)*d23_dx5+fi(3)*d33_dx5+d(1,3)*dfi1_dx5+d(2,3)*dfi2_dx5+d(3,3)*dfi3_dx5)
end if

err_flag = .false.
    
end subroutine ptc_fringe_dipoler

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_rot_xz(a, X, beta0, err_flag, mat6, make_matrix)
!
! Subroutine to rotate the local reference frame about the Y axis in PTC coordinates.  
! Adapted from forest/code/Sc_euclidean.f90 : ROT_XZ
!
! Input: 
!   a           -- real(rp): rotation angle (rad)
!   X(6)        -- real(rp): PTC phase space coordinates
!   beta0       -- real(rp): reference relativistic beta
!   make_matrix -- logical, optional: Make the matrix? Default is false.
!
! Output:
!   X(6)        -- real(rp): PTC phase space coordinates
!   err_flag    -- real(rp): Set true if ther is a problem.
!   mat6        -- real(rp), optional: Transfer matrix.
!-

subroutine ptc_rot_xz(a, X, beta0, err_flag, mat6, make_matrix)

implicit none

real(rp) :: x(6)
real(rp) :: x1,pz,pt
real(rp) :: a, beta0
character(*), parameter :: r_name = 'ptc_rot_xz'
real(rp), optional :: mat6(6,6)
real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx2, dpt_dx4, dpt_dx5, arg

logical err_flag
logical, optional :: make_matrix

!

err_flag = .true.

arg = 1.0_rp+2.0_rp*x(5)/ beta0+x(5)**2-x(2)**2-x(4)**2
if (arg < 0) return
pz=sqrt(arg)
pt=1.0_rp-x(2)*tan(a)/pz

if (logic_option(.false., make_matrix)) then 
  call mat_make_unit(mat6)

  dpz_dx2 = -X(2)/pz
  dpz_dx4 = -X(4)/pz
  dpz_dx5 = (1/beta0+X(5))/pz
  dpt_dx2 = -tan(a)/pz+x(2)*tan(a)*dpz_dx2/pz**2
  dpt_dx4 = x(2)*tan(a)*dpz_dx4/pz**2
  dpt_dx5 = x(2)*tan(a)*dpz_dx5/pz**2

  mat6(1,1) = 1/(cos(a)*pt)
  mat6(1,2) = -x(1)*dpt_dx2/(cos(a)*pt**2)
  mat6(1,4) = -x(1)*dpt_dx4/(cos(a)*pt**2)
  mat6(1,5) = -x(1)*dpt_dx5/(cos(a)*pt**2)
  mat6(2,2) = cos(a)+sin(a)*dpz_dx2
  mat6(2,4) = sin(a)*dpz_dx4
  mat6(2,5) = sin(a)*dpz_dx5
  mat6(3,1) = x(4)*tan(a)/(pt*pz)
  mat6(3,2) = -x(1)*x(4)*tan(a)*(dpt_dx2/(pt**2*pz)+dpz_dx2/(pt*pz**2))
  mat6(3,4) = -x(1)*x(4)*tan(a)*(dpt_dx4/(pt**2*pz)+dpz_dx4/(pt*pz**2))+x(1)*tan(a)/(pt*pz)
  mat6(3,5) = -x(1)*x(4)*tan(a)*(dpt_dx5/(pt**2*pz)+dpz_dx5/(pt*pz**2))
  mat6(6,1) = (1/beta0+x(5))*tan(a)/(pt*pz)
  mat6(6,2) = -x(1)*(1/beta0+x(5))*tan(a)*(dpt_dx2/(pt**2*pz)+dpz_dx2/(pt*pz**2))
  mat6(6,4) = -x(1)*(1/beta0+x(5))*tan(a)*(dpt_dx4/(pt**2*pz)+dpz_dx4/(pt*pz**2))
  mat6(6,5) = -x(1)*(1/beta0+x(5))*tan(a)*(dpt_dx5/(pt**2*pz)+dpz_dx5/(pt*pz**2))+x(1)*tan(a)/(pt*pz)
end if

x1=x(1)
x(1)=x1/cos(a)/pt
x(2)=x(2)*cos(a)+sin(a)*pz
x(3)=x(3)+x(4)*x1*tan(a)/pz/pt
x(6)=x(6)+x1*tan(a)/pz/pt*(1.0_rp/ beta0+x(5))

err_flag = .false.

end subroutine ptc_rot_xz

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the edge field of an sbend.
! Uses routines adapted from PTC
!
! Input:
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix up to the edge.
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix through the edge.
!-

subroutine exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

use ptc_interface_mod

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param
real(rp), optional :: mat6(6,6)
real(rp) :: mat6_int(6,6)
real(rp) :: X(6), ct
real(rp) :: beta0, g_tot, edge_angle, hgap, fint
integer :: particle_at

logical err_flag
logical, optional :: make_matrix
character(*), parameter :: r_name = 'exact_bend_edge_kick'

! Get reference beta0

beta0 = ele%value(p0c$) / ele%value(e_tot$)
if (ele%is_on) then
  g_tot = ele%value(g$) + ele%value(dg$)
else
  g_tot = 0
endif

! Convert to PTC coordinates
if (logic_option(.false., make_matrix)) then 
  call vec_bmad_to_ptc(orb%vec, beta0, X, mat6_int)
  mat6 = matmul(mat6_int, mat6)
else
  call vec_bmad_to_ptc(orb%vec, beta0, X)
end if

! get edge parameters
 
if (physical_ele_end(particle_at, orb, ele%orientation) == entrance_end$) then
  edge_angle = orb%time_dir * ele%value(e1$)
  fint = ele%value(FINT$)
  hgap = ele%value(HGAP$)
else
  edge_angle = orb%time_dir * ele%value(e2$)
  fint = ele%value(FINTX$)
  hgap = ele%value(HGAPX$)
endif

! Save time

ct = X(6)

if (particle_at == first_track_edge$) then
  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, orb%time_dir*fint, hgap, particle_at, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int,mat6)

  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

else if (particle_at == second_track_edge$) then
  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, orb%time_dir*fint, hgap, particle_at, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_pz$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

else
  !error!
  if (global_com%exit_on_error) call err_exit
endif

! Convert back to bmad coordinates
call vec_ptc_to_bmad (X, beta0, orb%vec, mat6_int)
if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

! Correct time
orb%t = orb%t + (X(6) - ct)/c_light

!--------------------------------------------
contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, conversion_mat)
!
! Routine to convert a Bmad vector map to PTC vector,
!
! Input:
!   vec_bmad(6) -- real(rp): Bmad coordinates.
!   beta0       -- real(rp): Reference particle velocity
!
! Output:
!   vec_ptc(6)     -- real(rp): PTC coordinates.
!   conversion_mat -- real(rp), optional: Jacobian matrix of Bmad -> PTC conversion map.
!-

subroutine vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, conversion_mat)

implicit none

real(rp) vec_bmad(:), vec_ptc(:)
real(rp) beta0, vec_temp(6)
real(rp), optional :: conversion_mat(6,6)
real(rp) factor1, factor2

! vec_ptc(5) = (E - E0) / P0c
! vec_ptc(6) = c (t - t0)
! 1/beta0 + vec_ptc(5) == E / P0c

vec_temp = vec_bmad
vec_temp(5) = (vec_bmad(6)**2 + 2*vec_bmad(6)) / (1/beta0 + sqrt(1/beta0**2+vec_bmad(6)**2+2*vec_bmad(6)) )
vec_temp(6) = -vec_bmad(5) * (1/beta0 + vec_temp(5)) / (1 + vec_bmad(6))

if (present(conversion_mat)) then
  call mat_make_unit(conversion_mat)
  factor1 = sqrt(1/beta0**2+vec_bmad(6)**2+2*vec_bmad(6))
  factor2 = 1+beta0**2*vec_bmad(6)*(2+vec_bmad(6))
  conversion_mat(5,5) = 0
  conversion_mat(5,6) = beta0**2*(1+vec_bmad(6))*factor1/factor2
  conversion_mat(6,5) = -(1/beta0+beta0*vec_bmad(6)*(2+vec_bmad(6))/(1+beta0*factor1))/(1+vec_bmad(6))
  conversion_mat(6,6) = -((beta0**2-1)*vec_bmad(5)*factor1)/((1+vec_bmad(6))**2*factor2)
end if

vec_ptc = vec_temp

end subroutine vec_bmad_to_ptc 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad, conversion_mat, state)
!
! Routine to convert a PTC orbit vector to a Bmad orbit vector.
!
! Input:
!   vec_ptc(6)  -- real(rp): PTC coordinates.
!   beta0       -- real(rp): Reference particle velocity
!
! Output:
!   vec_bmad(6)    -- real(rp): Bmad coordinates.
!   conversion_mat -- real(rp), optional: Jacobian matrix of PTC -> Bmad conversion map.
!   state          -- integer, optional: Set to lost_pz$ if energy is too low. Set to alive$ otherwise.
!-

subroutine vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad, conversion_mat, state)

implicit none

real(rp) vec_bmad(:), vec_ptc(:)
real(rp) beta0, vec_temp(6)
real(rp), optional :: conversion_mat(6,6)
real(rp) factor1, factor2
integer, optional :: state

!

if (present(state)) state = alive$

factor1 = 1+2*vec_ptc(5)/beta0+vec_ptc(5)**2
if (factor1 <= 0) then
  if (present(state)) state = lost_pz$
  return
endif

vec_temp = vec_ptc
vec_temp(6) = (2*vec_ptc(5)/beta0+vec_ptc(5)**2)/(sqrt(factor1)+1)
vec_temp(5) = -vec_ptc(6) * (1 + vec_temp(6)) / (1/beta0 + vec_ptc(5))

if (present(conversion_mat)) then
  call mat_make_unit(conversion_mat)
  factor1 = sqrt(1+2*vec_ptc(5)/beta0+vec_ptc(5)**2)
  factor2 = beta0+2*vec_ptc(5)+beta0*vec_ptc(5)**2 
  conversion_mat(5,5) = beta0*(beta0**2-1)*factor1*vec_ptc(6)/((1+beta0*vec_ptc(5))**2*factor2)
  conversion_mat(5,6) = -(1+vec_ptc(5)*(2+beta0*vec_ptc(5))/(beta0*(1+factor1)))/(1/beta0+vec_ptc(5))
  conversion_mat(6,5) = (1+beta0*vec_ptc(5))*factor1/factor2
  conversion_mat(6,6) = 0
end if

vec_bmad = vec_temp

end subroutine vec_ptc_to_bmad 

end subroutine exact_bend_edge_kick

end module
