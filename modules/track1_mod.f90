!+
! Module track1_mod
!
! Module of helper routines for track1 routines
!-

module track1_mod

use bmad_interface
use make_mat6_mod
use em_field_mod
use geometry_mod
use spin_mod

! Private routines for exact_bend_edge_kick
private ptc_rot_xz, ptc_wedger, ptc_fringe_dipoler

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, particle_at, param, check_momentum)
!
! Subroutine to check if an orbit is outside an element's aperture.
! Note: A particle will also be considered to have hit an aperture
! if |p_x| or |p_y| > 1 
!
! Also see:
!   orbit_too_large
!
! Modules needed:
!   use track1_mod
!
! Input:
!   orb            -- Coord_struct: coordinates of a particle.
!   ele            -- Ele_struct: Element holding the aperture
!   particle_at    -- Integer: first_track_edge$, second_track_edge$, surface$, in_between$
!   param          -- lat_param_struct: Parameter structure
!   check_momentum -- Logical, optional: If present and false then checking of
!                       p_x and p_y will be disabled.
!
! Output:
!   orb   -- Coord_struct: coordinates of a particle.
!     %state -- State of the particle
!-

recursive subroutine check_aperture_limit (orb, ele, particle_at, param, check_momentum)

implicit none

type (coord_struct) :: orb
type (coord_struct) orb2 
type (ele_struct) :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct), intent(inout) :: param

real(rp) x_lim, y_lim, x_particle, y_particle, s_here, r, rel_p, x0, y0
real(rp) x1_lim, x2_lim, y1_lim, y2_lim, dx1, dx2, dy1, dy2, d_max

integer i, particle_at, physical_end
logical do_tilt, err
logical, optional :: check_momentum
character(*), parameter :: r_name = 'check_aperture_limit'

! Super_slave elements have the aperture info stored in the lord

physical_end = physical_ele_end (particle_at, orb%direction, ele%orientation)

if (ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    if (.not. lord_edge_aligned (ele, physical_end, lord) .and. lord%aperture_at /= continuous$) cycle
    call check_aperture_limit (orb, lord, particle_at, param, check_momentum)
    if (orb%state /= alive$) return
  enddo
  return
endif

! Custom

if (ele%aperture_type == custom_aperture$) then
  call check_aperture_limit_custom (orb, ele, particle_at, param, err)
  return
endif

! Check p_x and p_y

if (logic_option(.true., check_momentum)) then
  if (orbit_too_large (orb, param)) return
endif

! Check if there is an aperture here. If not, simply return.

if (.not. at_this_ele_end (physical_end, ele%aperture_at)) return

! A photon at the surface will have the appropriate coords already so do not need to offset.

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)

  select case (ele%aperture_at)
  case (surface$)
    if (.not. ele%offset_moves_aperture) then 
      call out_io (s_error$, r_name, 'Surface aperture must have offset_moves_aperture = True for element: ' // ele%name)
    endif
  case default
    if (ele%offset_moves_aperture) then 
      call out_io (s_error$, r_name, 'Non-Surface aperture must have offset_moves_aperture = False for element: ' // ele%name)
    endif
  end select

  x_particle = orb%vec(1)
  y_particle = orb%vec(3)

case default
  if (ele%offset_moves_aperture .and. physical_end /= surface$) then
    do_tilt = .false.
    if (ele%key == ecollimator$ .or. ele%key == rcollimator$) do_tilt = .true.
    orb2 = orb
    s_here = orb2%s - (ele%s - ele%value(l$))
    if (orb2%species == photon$) then
      call offset_photon (ele, orb2, set$)
    else
      call offset_particle (ele, param, set$, orb2, set_tilt = do_tilt, &
                                 set_multipoles = .false., set_hvkicks = .false., ds_pos = s_here)
    endif
    x_particle = orb2%vec(1)
    y_particle = orb2%vec(3)

  else
    x_particle = orb%vec(1)
    y_particle = orb%vec(3)
  endif
end select

!

x_lim = (ele%value(x1_limit$) + ele%value(x2_limit$)) / 2
y_lim = (ele%value(y1_limit$) + ele%value(y2_limit$)) / 2
x0 = (ele%value(x2_limit$) - ele%value(x1_limit$)) / 2
y0 = (ele%value(y2_limit$) - ele%value(y1_limit$)) / 2
x_particle = x_particle - x0
y_particle = y_particle - y0

if (.not. bmad_com%aperture_limit_on .or. x_lim == 0 .or. y_lim == 0) then
  x_lim = bmad_com%max_aperture_limit
  y_lim = bmad_com%max_aperture_limit
endif

!

select case (ele%aperture_type)

case (elliptical$)

  r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
  if (r < 1) return

  if (abs(x_particle / x_lim) > abs(y_particle / y_lim)) then
    if (x_particle > 0) then; orb%state = lost_pos_x_aperture$
    else;                     orb%state = lost_neg_x_aperture$
    endif
  else
    if (y_particle > 0) then; orb%state = lost_pos_y_aperture$
    else;                     orb%state = lost_neg_y_aperture$
    endif
  endif
  param%unstable_factor = sqrt(r) - 1

case (rectangular$, auto_aperture$)

  if (abs(x_particle) < x_lim .and. abs(y_particle) < y_lim) return

  if (abs(x_particle/x_lim) > abs(y_particle/y_lim)) then
    if (x_particle > 0) then
      orb%state = lost_pos_x_aperture$
    else
      orb%state = lost_neg_x_aperture$
    endif
    param%unstable_factor = abs(x_particle/x_lim) - 1

  else
    if (y_particle > 0) then
      orb%state = lost_pos_y_aperture$
    else
      orb%state = lost_neg_y_aperture$
    endif
    param%unstable_factor = abs(y_particle/y_lim) - 1
  endif

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
end select

end subroutine check_aperture_limit

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function orbit_too_large (orbit, param) result (is_too_large)
!
! Routine to check if an orbit is too large.
! This routine is used to prevent floating point overflow.
! Too large is defined by:
!   |x|, |y| > bmad_com%max_aperture_limit or 
!   |pz| > 1 (photons only) or 
!   px^2 + py^2 > (1 + pz)^2 (non-photons only)
!
! Also see:
!   check_aperture_limit
!
! Input:
!   orbit         -- coord_struct: Particle orbit.
!
! Output:
!   orbit         -- coord_struct: Particle orbit.
!     %state          -- Particle status.
!   is_too_large  -- logical: True if orbit is too large. False otherwise.
!   param         -- lat_param_struct, optional: 
!     %unstable_factor  -- Set if orbit is too large. Otherwise not set
!-

function orbit_too_large (orbit, param) result (is_too_large)

implicit none

type (coord_struct) orbit
type (lat_param_struct), optional :: param

logical is_too_large
real(rp) rel_p

! Assume the worst

is_too_large = .true.

! Test aperture

if (orbit%vec(1) > bmad_com%max_aperture_limit) then
  orbit%state = lost_pos_x_aperture$
  if (present(param)) param%unstable_factor = abs(orbit%vec(1)) - bmad_com%max_aperture_limit
  return
elseif (-orbit%vec(1) > bmad_com%max_aperture_limit) then
  orbit%state = lost_neg_x_aperture$
  if (present(param)) param%unstable_factor = abs(orbit%vec(1)) - bmad_com%max_aperture_limit
  return
endif

if (orbit%vec(3) > bmad_com%max_aperture_limit) then
  orbit%state = lost_pos_y_aperture$
  if (present(param)) param%unstable_factor = abs(orbit%vec(3)) - bmad_com%max_aperture_limit
  return
elseif (-orbit%vec(3) > bmad_com%max_aperture_limit) then
  orbit%state = lost_neg_y_aperture$
  if (present(param)) param%unstable_factor = abs(orbit%vec(3)) - bmad_com%max_aperture_limit
  return
endif

! Test photons

if (orbit%species == photon$) then
  if (abs(orbit%vec(6)) > 1) then
    orbit%state = lost_z_aperture$
    if (present(param)) param%unstable_factor = abs(orbit%vec(6)) - 1
    return
  endif

! charged particle test

else
  rel_p = 1 + orbit%vec(6)

  if (orbit%vec(2)**2 + orbit%vec(4)**2 > rel_p**2) then
    if (present(param)) param%unstable_factor = sqrt(orbit%vec(2)**2 + orbit%vec(4)**2 - rel_p**2)

    if (abs(orbit%vec(2)) > abs(orbit%vec(4))) then
      if (orbit%vec(2) > 0) then
        orbit%state = lost_pos_x_aperture$
      else
        orbit%state = lost_neg_x_aperture$
      endif

    else
      if (orbit%vec(4) > 0) then
        orbit%state = lost_pos_y_aperture$
      else
        orbit%state = lost_neg_y_aperture$
      endif
    endif

    return
  endif

endif

! Passed tests.

is_too_large = .false.

end function orbit_too_large

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_bend (start_orb, ele, param, end_orb)
!
! Particle tracking through a bend element. 
! For e1 or e2 non-zero this subroutine treats the edges as thin quads.
! This routine is not meant for general use.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position.
!   ele        -- Ele_struct: Bend element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb    -- Coord_struct: End position.
!-

subroutine track_a_bend (start_orb, ele, param, end_orb)

use multipole_mod

implicit none

type (coord_struct), intent(in)  :: start_orb
type (coord_struct), intent(out) :: end_orb
type (ele_struct),   intent(inout)  :: ele
type (lat_param_struct), intent(inout) :: param

real(rp) angle, ct, st, x, px, y, py, z, pz, dpx_t, p_long
real(rp) rel_p, rel_p2, Dy, px_t, factor, g, g_err, c_dir, stg
real(rp) length, g_tot, eps, pxy2, ff, fg, k_2, alpha, beta, one_ct
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, z_2, s_y, c_y, r(6)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer n, n_step

logical has_nonzero_pole, drifting

!-----------------------------------------------------------------------

end_orb = start_orb
call offset_particle (ele, param, set$, end_orb, set_multipoles = .false.)

! Entrance edge kick

c_dir = ele%orientation * end_orb%direction * relative_tracking_charge(start_orb, param)
call bend_edge_kick (ele, param, first_track_edge$, end_orb)

! If we have a sextupole component then step through in steps of length ds_step

n_step = 1

call multipole_ele_to_kt(ele, .false., has_nonzero_pole, knl, tilt)
if (ele%value(k2$) /= 0 .or. has_nonzero_pole) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)
if (has_nonzero_pole) knl = knl * c_dir / n_step

! Set some parameters

length = ele%value(l$) / n_step
g = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * c_dir
g_err = g_tot - g
angle = g * length
pz = start_orb%vec(6)
rel_p  = 1 + pz
rel_p2 = rel_p**2
k_1 = ele%value(k1$) * c_dir
k_2 = ele%value(k2$) * c_dir
drifting = .false.

if (.not. ele%is_on) then
  g_err = -g
  g_tot = 0
  k_1 = 0
  k_2 = 0
endif

! 1/2 sextupole kick at the beginning.

if (k_2 /= 0) then
  end_orb%vec(2) = end_orb%vec(2) + k_2/2 * length * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
  end_orb%vec(4) = end_orb%vec(4) + k_2/2 * length * end_orb%vec(1) * end_orb%vec(3)
end if
if (has_nonzero_pole) call multipole_kicks (knl/2, tilt, end_orb)

! And track with n_step steps

do n = 1, n_step

  ! with k1 /= 0 use small angle approximation

  if (k_1 /= 0) then

    call sbend_body_with_k1_map (ele, param, n_step, end_orb, end_orb = end_orb)

  elseif (g == 0 .and. g_err == 0) then
    call track_a_drift (end_orb, length)
    drifting = .true.

  !-----------------------------------------------------------------------
  ! Track through main body...
  ! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  else

    ct = cos(angle)
    st = sin(angle)
    if (abs(angle) < 1d-7) then
      stg = length * (1 - angle**2 / 6)
      one_ct = length * angle / 2
    else
      stg = sin(angle) / g
      one_ct = (1 - ct) / g
    endif

    x  = end_orb%vec(1)
    px = end_orb%vec(2)
    y  = end_orb%vec(3)
    py = end_orb%vec(4)
    z  = end_orb%vec(5)
   
    pxy2 = px**2 + py**2
    if (rel_p2 - pxy2 < 0.01) then  ! somewhat arbitrary cutoff
      end_orb%state = lost$
      end_orb%vec(1) = 2 * bmad_com%max_aperture_limit
      end_orb%vec(3) = 2 * bmad_com%max_aperture_limit
      return
    endif 

    p_long = sqrt(rel_p2 - pxy2)

    ! The following is to make sure that a beam entering on-axis remains 
    ! *exactly* on-axis.

    if (pxy2 < 1d-5) then  
      ff = pxy2 / (2 * rel_p)
      fg = g * (pz - ff - ff*ff/2 - g_tot*x) - g_err
    else
      fg = g * p_long - g_tot * (1 + x * g)
    endif

    Dy  = sqrt(rel_p2 - py**2)
    px_t = px*ct + fg*stg
    dpx_t = -px*st*g + fg*ct

    if (abs(px) > Dy .or. abs(px_t) > Dy) then
      if (max(px, px_t) > 0) then; end_orb%state = lost_pos_x_aperture$
      else;                        end_orb%state = lost_neg_x_aperture$
      endif
      return
    endif

    if (abs(angle) < 1d-5 .and. abs(g_tot * length) < 1d-5) then
      end_orb%vec(1) = end_orb%vec(1) + length * px / p_long - &
                       g_tot * (length * Dy)**2 / (2 * p_long**3) + &
                       g * length * (length * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) / (2 * p_long**2)
    elseif (abs(g_tot) < 1d-5 * abs(g)) then
      alpha = p_long * ct - px * st
      end_orb%vec(1) = (p_long * (1 + g * x) - alpha) / (g * alpha) - &
                   g_tot * (Dy * (1 + g * x) * st)**2 / (2 * alpha**3 * g**2) + &
                   g_tot**2 * Dy**2 * ((1 + g * x) * st)**3 * (px * ct + p_long * st) / (2 * alpha**5 * g**3)
    else
      eps = px_t**2 + py**2
      if (eps < 1d-5 * rel_p2) then  ! use small angle approximation
        eps = eps / (2 * rel_p)
        end_orb%vec(1) = (pz + px * st - ct * p_long + g_tot * x * ct + ct - g_err * one_ct + &
                                                                 eps * (eps / (2 * rel_p) - 1)) / g_tot
      else
        end_orb%vec(1) = (sqrt(rel_p2 - eps) + px*st + g_tot*x*ct - p_long*ct) / g_tot - one_ct
      endif
    endif

    end_orb%vec(2) = px_t
    end_orb%vec(4) = py

    if (abs(g_tot) < 1d-5 * abs(g)) then
      beta = (1 + g * x) * st / (g * alpha) - &
             g_tot * (px * ct + p_long * st) * (st * (1 + g * x))**2 / (2 * g**2 * alpha**3)
      end_orb%vec(3) = y + py * beta
      end_orb%vec(5) = z + length  - rel_p * beta 
    else
      factor = (asin(px/Dy) - asin(px_t/Dy)) / g_tot
      end_orb%vec(3) = y + py * (angle/g_tot + factor)
      end_orb%vec(5) = z + length * (g_err - g*pz) / g_tot - rel_p * factor
    endif

  endif

  ! sextupole kick

  if (n == n_step) then
    if (k_2 /= 0) then
      end_orb%vec(2) = end_orb%vec(2) + k_2/2 * length * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
      end_orb%vec(4) = end_orb%vec(4) + k_2/2 * length * end_orb%vec(1) * end_orb%vec(3)
    end if
    if (has_nonzero_pole) call multipole_kicks (knl/2, tilt, end_orb)
  else
    if (k_2 /= 0) then
       end_orb%vec(2) = end_orb%vec(2) + k_2 * length * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
       end_orb%vec(4) = end_orb%vec(4) + k_2 * length * end_orb%vec(1) * end_orb%vec(3)
    end if
    if (has_nonzero_pole) call multipole_kicks (knl, tilt, end_orb)
  endif

enddo

! Track through the exit face. Treat as thin lens.
! Need low energy z correction except when using track_a_drift.

if (orbit_too_large(end_orb, param)) return
call bend_edge_kick (ele, param, second_track_edge$, end_orb)

call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false.)

if (.not. drifting) call track1_low_energy_z_correction (end_orb, ele, param)

end subroutine track_a_bend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the edge field of an sbend.
! Apply only the first order kick, which is edge focusing.
! This routine is the same as mad_hard_bend_edge_kick, except that only
! first order terms are applied.
!
! Module needed:
!   use track1_mod
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
! See MAD physics guide for writeup. Note that MAD does not have a g_err.
! Apply only the first order kick.  i.e. only edge focusing.

c_dir = relative_tracking_charge(orb, param) * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(g_err$)) * c_dir
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

if (particle_at == first_track_edge$) then
  orb%vec(1) = v0(1)
  orb%vec(2) = v0(2) + ht_x * v0(1)
  orb%vec(3) = v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3)
  if (logic_option(.false., make_matrix)) then
    mat6(2,:) = mat6(2,:) + ht_x * mat6(1,:)
    mat6(4,:) = mat6(4,:) + ht_y * mat6(3,:)
  end if
else
  orb%vec(1) = v0(1)
  orb%vec(2) = v0(2) + ht_x * v0(1)
  orb%vec(3) = v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3)
  if (logic_option(.false., make_matrix)) then
    mat6(2,:) = mat6(2,:) + ht_x * mat6(1,:)
    mat6(4,:) = mat6(4,:) + ht_y * mat6(3,:)
  end if
endif

end subroutine linear_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
!
! Subroutine to track through the edge field of an sbend.
!
! Module needed:
!   use track1_mod
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

subroutine mad_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) e, g_tot, fint_gap, ht_x, ht_y, cos_e, sin_e, tan_e, sec_e, v0(6), k1_eff
real(rp) ht2, hs2, c_dir, k1, kmat(6,6)
integer particle_at, element_end, fringe_type

logical, optional :: make_matrix
character(24), parameter :: r_name = 'mad_bend_edge_kick'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a g_err.

c_dir = relative_tracking_charge(orb, param) * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)
fringe_type = nint(ele%value(fringe_type$))

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(g_err$)) * c_dir
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
ht_x = g_tot * tan_e
ht2 = g_tot * tan_e**2
hs2 = g_tot * sec_e**2
k1_eff = k1 * c_dir

if (fint_gap == 0) then
  ht_y = -ht_x
else
  ht_y = -g_tot * tan(e - 2 * fint_gap * g_tot * (1 + sin_e**2) / cos_e)
endif

v0 = orb%vec

if (particle_at == first_track_edge$) then
  orb%vec(1) = v0(1) - ht2 * v0(1)**2 / 2 + hs2 * v0(3)**2 / 2
  orb%vec(2) = v0(2) + ht_x * v0(1) + ht2 * (v0(1) * v0(2) - v0(3) * v0(4)) + &
                       k1_eff * tan_e * (v0(1)**2 - v0(3)**2) + &
                       ht_x * (ht2 + hs2) * v0(3)**2 / 2
  orb%vec(3) = v0(3) + ht2 * v0(1) * v0(3) 
  orb%vec(4) = v0(4) + ht_y * v0(3) - ht2 * v0(1) * v0(4) - hs2 * v0(2) * v0(3) - &
                       2 * k1_eff * tan_e * v0(1) * v0(3)
  if (logic_option(.false., make_matrix)) then
    call mat_make_unit (kmat)
    kmat(1,1) = 1 - ht2 * v0(1)
    kmat(1,3) = hs2 * v0(3)
    kmat(2,1) = ht_x + ht2 * v0(2) + 2 * k1_eff * tan_e * v0(1)
    kmat(2,2) = 1 + ht2 * v0(1)
    kmat(2,3) = - ht2 * v0(4) - 2 * k1_eff * tan_e * v0(3) + ht_x * (ht2 + hs2) * v0(3)
    kmat(2,4) = - ht2 * v0(3)
    kmat(3,1) = ht2 * v0(3)
    kmat(3,3) = 1 +  ht2 * v0(1)
    kmat(4,1) = - ht2 * v0(4) - 2 * tan_e * k1_eff * v0(3)
    kmat(4,2) = - hs2 * v0(3)
    kmat(4,3) = ht_y - hs2 * v0(2) - 2 * tan_e * k1_eff * v0(1)
    kmat(4,4) = 1 - ht2 * v0(1)
    mat6 = matmul (kmat, mat6)
  end if
else
  orb%vec(1) = v0(1)
  orb%vec(2) = v0(2) + ht_x * v0(1)
  orb%vec(3) = v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3)
  orb%vec(1) = v0(1) + ht2 * v0(1)**2 / 2 - hs2 * v0(3)**2 / 2
  orb%vec(2) = v0(2) + ht_x * v0(1) + ht2 * (v0(3) * v0(4) - v0(1) * v0(2)) + &
                       k1_eff * tan_e * (v0(1)**2 - v0(3)**2) - &
                       ht_x * ht2 * (v0(1)**2 + v0(3)**2) / 2
  orb%vec(3) = v0(3) - ht2 * v0(1) * v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3) + ht2 * v0(1) * v0(4) + hs2 * v0(2) * v0(3) + &
                       (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(1) * v0(3)
  if (logic_option(.false., make_matrix)) then
    call mat_make_unit (kmat)
    kmat(1,1) = 1 + ht2 * v0(1)
    kmat(1,3) = - hs2 * v0(3)
    kmat(2,1) = ht_x - ht2 * v0(2) + 2 * k1_eff * tan_e * v0(1) - ht_x * ht2 * v0(1)
    kmat(2,2) = 1 - ht2 * v0(1)
    kmat(2,3) = ht2 * v0(4) - 2 * k1_eff * tan_e * v0(3) - ht_x * ht2 * v0(3)
    kmat(2,4) = ht2 * v0(3)
    kmat(3,1) = - ht2 * v0(3)
    kmat(3,3) = 1 -  ht2 * v0(1)
    kmat(4,1) = ht2 * v0(4) + (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(3)
    kmat(4,2) = hs2 * v0(3)
    kmat(4,3) = ht_y + hs2 * v0(2) + (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(1)
    kmat(4,4) = 1 + ht2 * v0(1)
    mat6 = matmul (kmat, mat6)
  end if
endif

end subroutine mad_bend_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine bmad_standard_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the hard edge bend fringe field for a bend.
! Only the bend field is taken into account here. Higher order multipolse must be handled elsewhere.
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

subroutine bmad_standard_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) p_vec(3), r_vec(3), sin_e, cos_e, vec_in(6), e_angle
real(rp) rel_pc

integer fringe_type, fringe_at, physical_end, particle_at
integer i, i_max, element_end

logical, optional :: make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= full$ .and. fringe_type /= hard_edge_only$) return

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

!

element_end = physical_ele_end(particle_at, orbit%direction, ele%orientation)
if (element_end == entrance_end$) then
  e_angle = ele%value(e1$)
else
  e_angle = ele%value(e2$)
endif

sin_e = sin(e_angle)
cos_e = cos(e_angle)

!

rel_pc = 1 + orbit%vec(6)
vec_in = orbit%vec

p_vec = [orbit%vec(2), orbit%vec(4), orbit%direction * sqrt(rel_pc**2 - orbit%vec(2)**2 - orbit%vec(4)**2)]
r_vec = [orbit%vec(1), orbit%vec(3), 0.0_rp]

if (particle_at == first_track_edge$) then
  call track_between_physical_and_nominal_edges (.false.)
  call no_edge_angle_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
  call track_between_physical_and_nominal_edges (.true.)
  
else
  call track_between_physical_and_nominal_edges (.true.)
  call no_edge_angle_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
  call track_between_physical_and_nominal_edges (.false.)
endif


!----------------------------------------------------------------------------------------------
contains

! Transport between the nominal edge and the physical bend edge.
! The nominal edge is the same as the physical edge when e1/e2 are zero.

subroutine track_between_physical_and_nominal_edges (in_field_region)

real(rp) g, ll
logical in_field_region

! If in the bend field region then we are transporting "back" from the physical edge to the nominal edge.

if (in_field_region) sin_e = -sin_e

! Field free transport

p_vec = [p_vec(1) * cos_e + p_vec(3) * sin_e, p_vec(2), -p_vec(1) * sin_e + p_vec(3) * cos_e]
r_vec = [r_vec(1) * cos_e, r_vec(2), -r_vec(1) * sin_e]
ll = -r_vec(3) * rel_pc / p_vec(3)

orbit%vec(1) = orbit%vec(1) - r_vec(3)  * p_vec(1) / p_vec(3)
orbit%vec(2) = p_vec(1)
orbit%vec(3) = orbit%vec(3) - r_vec(3)  * p_vec(2) / p_vec(3)
orbit%vec(4) = p_vec(2)
orbit%vec(5) = orbit%vec(5) - ll

! Treat the bend field as a perturbation to the field free transport.
! That is, the path length is small compared to the bending radius.

if (in_field_region) then
  g = (ele%value(g$) + ele%value(g_err$) + vec_in(1) + ele%value(k1$) * ll * vec_in(2) / (2 * rel_pc)) / rel_pc
  orbit%vec(1) = orbit%vec(1) - ll**2 * g / 2
  orbit%vec(2) = orbit%vec(2) - ll * g
endif

!

if (make_matrix) then

endif

end subroutine track_between_physical_and_nominal_edges

end subroutine bmad_standard_hard_bend_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine no_edge_angle_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix, g_bend, tilt)
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
!   g_bend      -- real(rp), optional: Bend bend strength. Used (and must be present) only with sad_mult.
!   tilt        -- real(rp), optional: field rotation. Used (and must be present) only with sad_mult.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine no_edge_angle_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix, g_bend, tilt)

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (lat_param_struct) param

real(rp), optional :: mat6(6,6), g_bend, tilt
real(rp) g, px, y, y2, rel_p, p_zy, yg, kmat(6,6), c_dir, t0, ppx2

integer fringe_type, fringe_at, physical_end, particle_at
integer i, i_max, element_end

logical, optional :: make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= full$ .and. fringe_type /= hard_edge_only$) return

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

! extract params from ele

if (ele%key == sbend$) then
  g = ele%value(g$) + ele%value(g_err$)
  t0 = 0
else
  g = g_bend
  t0 = tilt
endif

!

if (g == 0) return
c_dir = relative_tracking_charge(orbit, param) * ele%orientation * orbit%direction
g = g * c_dir

if (particle_at == second_track_edge$) then
  g = -g
endif

! Rotate

if (t0 /= 0) call tilt_coords (t0, orbit%vec)

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

orbit%vec(1) = orbit%vec(1) + g * y2 * (1 - yg) * rel_p**2 / (2 * p_zy**3)
orbit%vec(4) = orbit%vec(4) - g * px * y * (1 - 2 * yg) / p_zy
orbit%vec(5) = orbit%vec(5) - g * y2 * px * (1 - yg) * rel_p / (2 * p_zy**3)

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
  call tilt_mat6(kmat, t0)
  mat6 = matmul(kmat, mat6)
endif

! Rotate

if (t0 /= 0) call tilt_coords (-t0, orbit%vec)

end subroutine no_edge_angle_hard_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine soft_quadrupole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)
!
! Routine to add the SAD "linear" soft edge (for finite f1 or f2).
! This routine assumes that the particle orbit has been rotated to the element reference frame.
!
! Moudle needed:
!   use track1_mod
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
real(rp) f1, f2, ef1, vec(4), rel_p, vx, vy

integer particle_at
integer fringe_at, physical_end, fringe_type

logical, optional :: make_matrix

!

if (ele%key == sbend$) then
 fringe_type = nint(ele%value(higher_order_fringe_type$))
else
 fringe_type = nint(ele%value(fringe_type$))
endif

if (fringe_type /= soft_edge_only$ .and. fringe_type /= full$) return


fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

charge_dir = ele%orientation * orbit%direction
if (associated(ele%branch)) charge_dir = charge_dir * relative_tracking_charge(orbit, param)

rel_p = 1 + orbit%vec(6)

select case (ele%key)
case (quadrupole$)
  k1_rel = charge_dir * ele%value(k1$) / rel_p
case (sad_mult$)
  ! Slice slaves and super slaves have their associated multipoles stored in the lord
  if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
    m_ele => pointer_to_lord(ele, 1)
  else
    m_ele => ele
  endif
  k1_rel = charge_dir * sqrt(m_ele%a_pole(1)**2 + m_ele%b_pole(1)**2)  / m_ele%value(l$) / rel_p
end select

f1 = k1_rel * ele%value(fq1$)
f2 = k1_rel * ele%value(fq2$)
if (f1 == 0 .and. f2 == 0) return
if (particle_at == second_track_edge$) f1 = -f1

! 

ef1 = exp(f1)

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

end subroutine soft_quadrupole_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine hard_multipole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix, a_pole, b_pole)
!
! Routine to track through the hard edge field of a multipole.
! The dipole component is ignored and only quadrupole and higher multipoles are included.
!
! This routine handles elements of type:
!   sad_mult, sbend, quadrupole, sextupole
!
! For sad_mult elements, a_pole and b_pole ae used for the multipole values.
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
!   a_pole(0:)  -- real(rp), optional: sad_mult skew components. Must be present for sad_mult.
!   b_pole(0:)  -- real(rp), optional: sad_mult normal components. Must be present for sad_mult.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine hard_multipole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix, a_pole, b_pole)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit

real(rp), optional :: mat6(6,6)
real(rp), optional ::  a_pole(0:), b_pole(0:)
real(rp) rel_p, cn, x, y, px, py, denom, ddenom_dx, ddenom_dy, ddenom_dpz, kmat(6,6), ap, bp, charge_dir
real(rp) fx, dfx_dx, dfx_dy, d2fx_dxx, d2fx_dxy, d2fx_dyy
real(rp) fy, dfy_dx, dfy_dy, d2fy_dxx, d2fy_dxy, d2fy_dyy

complex(rp) poly, poly_n1, poly_n2, dpoly_dx, dpoly_dy, d2poly_dxx, d2poly_dxy, d2poly_dyy
complex(rp) xy, xny, dxny_dx, dxny_dy, cab

integer fringe_type, fringe_at, physical_end, particle_at
integer n, n_max

logical, optional :: make_matrix

! Fringe here?

if (ele%key == sbend$) then
  fringe_type = nint(ele%value(higher_order_fringe_type$))
else
  fringe_type = nint(ele%value(fringe_type$))
endif

select case (fringe_type)
case (hard_edge_only$, full$)
case default
  return
end select

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

charge_dir = ele%orientation * orbit%direction
if (associated(ele%branch)) charge_dir = charge_dir * relative_tracking_charge(orbit, param)

!

select case (ele%key)
case (sad_mult$)
  do n = ubound(a_pole, 1), 0, -1
    if (a_pole(n) /= 0 .or. b_pole(n) /= 0) exit
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
    if (a_pole(n) == 0 .and. b_pole(n) == 0) cycle
    ap = a_pole(n) / ele%value(l$)
    bp = b_pole(n) / ele%value(l$)
  case (quadrupole$, sbend$)
    ap = 0
    bp = ele%value(k1$)
  case (sextupole$)
    ap = 0
    bp = ele%value(k2$)
    if (n /= 2) cycle
  end select

  dpoly_dx = (n+1) * poly_n1
  dpoly_dy = i_imaginary * dpoly_dx

  d2poly_dxx = n * (n+1) * poly_n2
  d2poly_dxy = i_imaginary * d2poly_dxx
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
orbit%vec(2) = ((1 - dfy_dy) * px + dfy_dx * py) / denom
orbit%vec(3) = orbit%vec(3) - fy
orbit%vec(4) = (dfx_dy * px + (1 - dfx_dx) * py) / denom
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
  mat6 = matmul (kmat, mat6)
endif

end subroutine hard_multipole_edge_kick 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, track_spin)
!
! Subroutine to track through the edge field of an sbend.
!
! Module needed:
!   use track1_mod
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

integer particle_at, fringe_type, physical_end

logical, optional :: make_matrix, track_spin
character(*), parameter :: r_name = 'bend_edge_kick'

!

physical_end = physical_ele_end (particle_at, orb%direction, ele%orientation)
if (.not. at_this_ele_end (physical_end, nint(ele%value(fringe_at$)))) return

! Higher order fringes. 
! Remember: In a bend, these fringes are turned on/off by ele%value(higher_order_fringe_type$)

if (logic_option(.false., track_spin)) ave_orb = orb

if (particle_at == first_track_edge$) then
  call hard_multipole_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
endif

! Bend Fringe

fringe_type = nint(ele%value(fringe_type$))

select case (fringe_type)
case (full$)
  call exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (basic_bend$, hard_edge_only$)
  call mad_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (sad_full$)
  if (particle_at == first_track_edge$) then
    call soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
    call mad_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
  else
    call mad_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
    call soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)
  endif

case (soft_edge_only$)
  call soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (linear_edge$)
  call linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

case (test_edge$)
  call bmad_standard_hard_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

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
  ave_orb%vec = (ave_orb%vec + orb%vec) / 2   ! Use average position
  field%E = 0
  if (physical_end == entrance_end$) then
    e_ang = ele%value(e1$)
  else
    e_ang = ele%value(e2$)
  endif
  x = ave_orb%vec(1);  y = ave_orb%vec(3)
  tan_e_x = tan(e_ang) * x
  field%B = (ele%value(b_field$) + ele%value(b_field_err$)) * [-sin(e_ang)*y, -tan_e_x, cos(e_ang)*y]
  if (ele%value(b1_gradient$) /= 0) field%B = field%B - ele%value(b1_gradient$) * tan_e_x * [x*y, x*x - y*y, 0.0_rp]
  if (ele%value(b2_gradient$) /= 0) field%B = field%B - ele%value(b2_gradient$) * tan_e_x * [3*x*x*y - y**3, x**3 - 3*x*y*y, 0.0_rp]
  if (particle_at == second_track_edge$) field%B(3) = -field%B(3)
  omega = spin_omega (field, ave_orb)
  call rotate_spinor (omega, orb%spin)
endif

end subroutine bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, k0l, t0)
!
! Subroutine to track through the ("linear") bend soft edge field of an sbend or sad_mult.
!
! Module needed:
!   use track1_mod
!
! Input:
!   ele         -- ele_struct: SBend or sad_mult element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   orb         -- Coord_struct: Starting coords.
!   mat6(6,6)   -- real(rp), optional: Starting matrix 
!   make_matrix -- real(rp), optional: Propagate the transfer matrix? Default is False.
!   k0l         -- real(rp), optional: Used with sad_mult. 
!                     If present, use this instead of ele%a/b_pole.
!   t0          -- real(rp), optional: Used with sad_mult. 
!                     If present, use this instead of ele%a/b_pole.
!                     Must be present if k0l is.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix after fringe field
!-

subroutine soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, k0l, t0)

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6), k0l, t0
real(rp) :: kmat(6,6), tilt
real(rp) :: f1, el_p, g, ct, c1, c2, c3, y, px, rel_p, sin_e
real(rp) k0, sk0, c1_k, c2_k, c3_k, c1_sk, c2_sk, c3_sk 

integer :: particle_at, c_dir, element_end, fringe_type

logical, optional :: make_matrix

character(*), parameter :: r_name = 'soft_bend_edge_kick'

!

fringe_type = nint(ele%value(fringe_type$))

if (ele%key == sbend$) then
  if (fringe_type /= soft_edge_only$ .and. fringe_type /= sad_full$) return
else  ! Must be sad_mult
  if (fringe_type /= soft_edge_only$ .and. fringe_type /=  full$) return
endif

! 

element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)

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

  g = ele%value(g$) + ele%value(g_err$)
  tilt = 0

case (sad_mult$)

  if (element_end == entrance_end$) then
    f1 = ele%value(fb1$)
  else
    f1 = ele%value(fb2$)
  endif

  if (f1 == 0) return

  if (present(k0l)) then
    g = k0l / ele%value(l$)
    tilt = t0
  else
    g = sqrt(ele%b_pole(0)**2 + ele%a_pole(0)**2) / ele%value(l$)
    tilt = -atan2(ele%a_pole(0), ele%b_pole(0))
  endif

end select

!

if (g == 0) return
c_dir = relative_tracking_charge(orb, param) * ele%orientation * orb%direction
g = g * c_dir

if (particle_at == second_track_edge$) then
  sin_e = -sin_e
  g = -g
endif

if (tilt /= 0) call tilt_coords(tilt, orb%vec)

px = orb%vec(2) ! + sin_e ???
y  = orb%vec(3)
rel_p = 1 + orb%vec(6)

c1 = f1**2 * g / (24 * rel_p)  ! * px
c2 = f1 * g**2 / (6 * rel_p)  ! * y^2
if (ele%key == sad_mult$ .or. (ele%key == sbend$ .and. &
                (fringe_type == soft_edge_only$ .or. fringe_type == sad_full$))) then
  c3 = 2 * g**2 / (3 * f1 * rel_p)   ! * y^4
else
  c3 = 0
endif

if (logic_option(.false., make_matrix)) then
  call mat_make_unit (kmat)
  kmat(1,6) =  c1 / rel_p
  kmat(4,3) =  c2 - 3 * c3 * y**2
  kmat(4,6) = (-c2 * y + c3 * y**3) / rel_p
  kmat(5,2) = c1 / rel_p
  kmat(5,3) = (c2 * y - c3 * y**3) / rel_p
  kmat(5,6) = (-c1 * px - 2 * c2 * y**2 + c3 * y**4/2) / rel_p**2 + c1

  if (tilt /= 0) call tilt_mat6 (kmat, tilt)

  mat6(1,:) = mat6(1,:) + kmat(1,6) * mat6(6,:)
  mat6(4,:) = mat6(4,:) + kmat(4,3) * mat6(3,:) + kmat(4,6) * mat6(6,:)
  mat6(5,:) = mat6(5,:) + kmat(5,2) * mat6(2,:) + kmat(5,3) * mat6(3,:) + kmat(5,6) * mat6(6,:)
endif

orb%vec(1) = orb%vec(1) + c1 * orb%vec(6)
orb%vec(4) = orb%vec(4) + c2 * y - c3 * y**3
orb%vec(5) = orb%vec(5) + (c1 * px + c2 * y**2/2 - c3 * y**4/4) / rel_p

if (tilt /= 0) call tilt_coords (-tilt, orb%vec)

end subroutine soft_bend_edge_kick

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
character(20) :: r_name = 'ptc_wedger'
real(rp), optional :: mat6(6,6)

real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx4, dpt_dx5, fac
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
pz=sqrt(fac)
dpz_dx2 = -X(2)/pz
dpz_dx4 = -X(4)/pz
dpz_dx5 = (1/beta0+X(5))/pz

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

Xn(2)=X(2)*cos(a)+(pz-b1*X(1))*sin(a)
if (logic_option(.false., make_matrix)) then
  mat6(2,1) = -b1*sin(a)
  mat6(2,2) = cos(a)+sin(a)*dpz_dx2
  mat6(2,4) = sin(a)*dpz_dx4
  mat6(2,5) = sin(a)*dpz_dx5
end if

pt=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(4)**2)
dpt_dx4 = -X(4)/pt
dpt_dx5 = (1/beta0+X(5))/pt

pzs=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-Xn(2)**2-X(4)**2)
if (logic_option(.false., make_matrix)) then
  dpzs_dx1 = -Xn(2)*mat6(2,1)/pzs
  dpzs_dx2 = -Xn(2)*mat6(2,2)/pzs
  dpzs_dx4 = (-X(4)-Xn(2)*mat6(2,4))/pzs
  dpzs_dx5 = (1/beta0+X(5)-Xn(2)*mat6(2,5))/pzs
end if

Xn(1)=X(1)*cos(a)+(X(1)*X(2)*sin(2.0_rp*a)+sin(a)**2*(2.0_rp*X(1)*pz-b1*X(1)**2) )&
    / (pzs+pz*cos(a)-X(2)*sin(a))
if (logic_option(.false., make_matrix)) then
  factor1 = cos(a)*pz+pzs-X(2)*sin(a)
  factor2 = (-b1*X(1)**2+2*X(1)*pz)*sin(a)**2+X(1)*X(2)*sin(2*a)
  mat6(1,1) = cos(a)+((-2*b1*X(1)+2*pz)*sin(a)**2+X(2)*sin(2*a))/factor1 &
              -(factor2*dpzs_dx1)/factor1**2
  mat6(1,2) = (X(1)*sin(2*a)+2*X(1)*sin(a)**2*dpz_dx2)/factor1 &
              -(factor2*(-sin(a)+cos(a)*dpz_dx2+dpzs_dx2))/factor1**2
  mat6(1,4) = (2*X(1)*sin(a)**2*dpz_dx4)/factor1 &
              -(factor2*(cos(a)*dpz_dx4+dpzs_dx4))/factor1**2
  mat6(1,5) = (2*X(1)*sin(a)**2*dpz_dx5)/factor1 &
              -(factor2*(cos(a)*dpz_dx5+dpzs_dx5))/factor1**2
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
real(rp) :: PZ,XP,YP,TIME_FAC
real(rp) :: D(3,3),FI(3),FI0,B,co1,co2
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
dpz_dx2 = -X(2)/pz
dpz_dx4 = -X(4)/pz
dpz_dx5 = (1/beta0+X(5))/pz

time_fac=1.0_rp/beta0+x(5)
dtime_fac_dx5 = 1

xp=x(2)/pz
dxp_dx2 = 1/pz-x(2)*dpz_dx2/pz**2
dxp_dx4 = -x(2)*dpz_dx4/pz**2
dxp_dx5 = -x(2)*dpz_dx5/pz**2

yp=x(4)/pz
dyp_dx2 = -x(4)*dpz_dx2/pz**2
dyp_dx4 = 1/pz-x(4)*dpz_dx4/pz**2
dyp_dx5 = -x(4)*dpz_dx5/pz**2

d(1,1)=(1.0_rp+xp**2)/pz
d11_dx2 = 2*xp*dxp_dx2/pz-(1+xp**2)*dpz_dx2/pz**2
d11_dx4 = 2*xp*dxp_dx4/pz-(1+xp**2)*dpz_dx4/pz**2
d11_dx5 = 2*xp*dxp_dx5/pz-(1+xp**2)*dpz_dx5/pz**2

d(2,1)=xp*yp/pz
d21_dx2 = xp*dyp_dx2/pz+yp*dxp_dx2/pz-xp*yp*dpz_dx2/pz**2
d21_dx4 = xp*dyp_dx4/pz+yp*dxp_dx4/pz-xp*yp*dpz_dx4/pz**2
d21_dx5 = xp*dyp_dx5/pz+yp*dxp_dx5/pz-xp*yp*dpz_dx5/pz**2

d(3,1)=-xp
d31_dx2 = -dxp_dx2
d31_dx4 = -dxp_dx4
d31_dx5 = -dxp_dx5

d(1,2)=xp*yp/pz
d12_dx2 = d21_dx2
d12_dx4 = d21_dx4
d12_dx5 = d21_dx5

d(2,2)=(1.0_rp+yp**2)/pz
d22_dx2 = 2*yp*dyp_dx2/pz-(1+yp**2)*dpz_dx2/pz**2
d22_dx4 = 2*yp*dyp_dx4/pz-(1+yp**2)*dpz_dx4/pz**2
d22_dx5 = 2*yp*dyp_dx5/pz-(1+yp**2)*dpz_dx5/pz**2

d(3,2)=-yp
d32_dx2 = -dyp_dx2
d32_dx4 = -dyp_dx4
d32_dx5 = -dyp_dx5

d(1,3)=-time_fac*xp/pz**2
d13_dx2 = 2*time_fac*xp*dpz_dx2/pz**3-time_fac*dxp_dx2/pz**2
d13_dx4 = 2*time_fac*xp*dpz_dx4/pz**3-time_fac*dxp_dx4/pz**2
d13_dx5 = 2*time_fac*xp*dpz_dx5/pz**3-time_fac*dxp_dx5/pz**2-xp*dtime_fac_dx5/pz**2

d(2,3)=-time_fac*yp/pz**2
d23_dx2 = 2*time_fac*yp*dpz_dx2/pz**3-time_fac*dyp_dx2/pz**2
d23_dx4 = 2*time_fac*yp*dpz_dx4/pz**3-time_fac*dyp_dx4/pz**2
d23_dx5 = 2*time_fac*yp*dpz_dx5/pz**3-time_fac*dyp_dx5/pz**2-yp*dtime_fac_dx5/pz**2

d(3,3)= time_fac/pz
d33_dx2 = -time_fac*dpz_dx2/pz**2
d33_dx4 = -time_fac*dpz_dx4/pz**2
d33_dx5 = -time_fac*dpz_dx5/pz**2-dtime_fac_dx5/pz

fi0= atan((xp/(1.0_rp+yp**2)))-b*fint*hgap*2.0_rp*( 1.0_rp + xp**2*(2.0_rp+yp**2) )*pz
factor1 = b*fint*hgap
factor2 = 1+yp**2
dfi0_dx2 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx2-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx2+2*xp**2*yp*dyp_dx2) &
           +(dxp_dx2/factor2-2*xp*yp*dyp_dx2/factor2**2)/(1+(xp/factor2)**2) 
dfi0_dx4 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx4-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx4+2*xp**2*yp*dyp_dx4) &
           +(dxp_dx4/factor2-2*xp*yp*dyp_dx4/factor2**2)/(1+(xp/factor2)**2)
dfi0_dx5 = -2*factor1*(1+xp**2*(2+yp**2))*dpz_dx5-2*factor1*pz*(2*xp*(2+yp**2)*dxp_dx5+2*xp**2*yp*dyp_dx5) &
           +(dxp_dx5/factor2-2*xp*yp*dyp_dx5/factor2**2)/(1+(xp/factor2)**2)

co2=b/cos(fi0)**2
dco2_dx2 = 2*b*tan(fi0)*dfi0_dx2/cos(fi0)**2
dco2_dx4 = 2*b*tan(fi0)*dfi0_dx4/cos(fi0)**2
dco2_dx5 = 2*b*tan(fi0)*dfi0_dx5/cos(fi0)**2

co1=co2/(1.0_rp+(xp/(1.0_rp+yp**2))**2 )
dco1_dx2 = dco2_dx2/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx2/factor2**2-4*xp**2*yp*dyp_dx2/factor2**3)/(1+(xp/factor2)**2)**2
dco1_dx4 = dco2_dx4/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx4/factor2**2-4*xp**2*yp*dyp_dx4/factor2**3)/(1+(xp/factor2)**2)**2
dco1_dx5 = dco2_dx5/(1+(xp/factor2)**2)-co2*(2*xp*dxp_dx5/factor2**2-4*xp**2*yp*dyp_dx5/factor2**3)/(1+(xp/factor2)**2)**2

fi(1)=co1/(1.0_rp+yp**2)-co2*b*fint*hgap*2.0_rp*( 2.0_rp*xp*(2.0_rp+yp**2)*pz )
dfi1_dx2 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx2-4*factor1*co2*xp*(2+yp**2)*dpz_dx2-4*factor1*co2*pz*(2+yp**2)*dxp_dx2 &
           -8*factor1*co2*pz*xp*yp*dyp_dx2+dco1_dx2/factor2-2*co1*yp*dyp_dx2/factor2**2
dfi1_dx4 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx4-4*factor1*co2*xp*(2+yp**2)*dpz_dx4-4*factor1*co2*pz*(2+yp**2)*dxp_dx4 &
           -8*factor1*co2*pz*xp*yp*dyp_dx4+dco1_dx4/factor2-2*co1*yp*dyp_dx4/factor2**2
dfi1_dx5 = -4*factor1*pz*xp*(2+yp**2)*dco2_dx5-4*factor1*co2*xp*(2+yp**2)*dpz_dx5-4*factor1*co2*pz*(2+yp**2)*dxp_dx5 &
           -8*factor1*co2*pz*xp*yp*dyp_dx5+dco1_dx5/factor2-2*co1*yp*dyp_dx5/factor2**2

fi(2)=-co1*2.0_rp*xp*yp/(1.0_rp+yp**2)**2-co2*b*fint*hgap*2.0_rp*( 2.0_rp*xp**2*yp)*pz
dfi2_dx2 = -4*factor1*pz*xp**2*yp*dco2_dx2-4*factor1*co2*xp**2*yp*dpz_dx2-8*factor1*co2*pz*xp*yp*dxp_dx2 &
           -4*factor1*co2*pz*xp**2*dyp_dx2-2*xp*yp*dco1_dx2/factor2**2-2*co1*yp*dxp_dx2/factor2**2 &
           +8*co1*xp*yp**2*dyp_dx2/factor2**3-2*co1*xp*dyp_dx2/factor2**2
dfi2_dx4 = -4*factor1*pz*xp**2*yp*dco2_dx4-4*factor1*co2*xp**2*yp*dpz_dx4-8*factor1*co2*pz*xp*yp*dxp_dx4 &
           -4*factor1*co2*pz*xp**2*dyp_dx4-2*xp*yp*dco1_dx4/factor2**2-2*co1*yp*dxp_dx4/factor2**2 &
           +8*co1*xp*yp**2*dyp_dx4/factor2**3-2*co1*xp*dyp_dx4/factor2**2
dfi2_dx5 = -4*factor1*pz*xp**2*yp*dco2_dx5-4*factor1*co2*xp**2*yp*dpz_dx5-8*factor1*co2*pz*xp*yp*dxp_dx5 &
           -4*factor1*co2*pz*xp**2*dyp_dx5-2*xp*yp*dco1_dx5/factor2**2-2*co1*yp*dxp_dx5/factor2**2 &
           +8*co1*xp*yp**2*dyp_dx5/factor2**3-2*co1*xp*dyp_dx5/factor2**2

fi(3)=-co2*b*fint*hgap*2.0_rp*( 1.0_rp + xp**2*(2.0_rp+yp**2) )
dfi3_dx2 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx2-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx2+2*xp**2*yp*dyp_dx2)
dfi3_dx4 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx4-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx4+2*xp**2*yp*dyp_dx4)
dfi3_dx5 = -2*factor1*(1+xp**2*(2+yp**2))*dco2_dx5-2*co2*factor1*(2*xp*(2+yp**2)*dxp_dx5+2*xp**2*yp*dyp_dx5)

dfi0_dx2 = b*dfi0_dx2/cos(fi0)**2
dfi0_dx4 = b*dfi0_dx4/cos(fi0)**2
dfi0_dx5 = b*dfi0_dx5/cos(fi0)**2
fi0=b*tan(fi0)

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

b=0
do i=1,3
   b=fi(i)*d(i,2)+b
enddo

if (logic_option(.false., make_matrix)) then
  factor1 = sqrt(1-2*b*x(3))
  factor2 = (1+sqrt(1-2*b*x(3)))**2
  mat6(3,2) = 2*x(3)**2*(fi(1)*d12_dx2+fi(2)*d22_dx2+fi(3)*d32_dx2+d(1,2)*dfi1_dx2+d(2,2)*dfi2_dx2+d(3,2)*dfi3_dx2)/(factor1*factor2)
  mat6(3,3) = 2/(1+factor1)+2*x(3)*b/(factor1*factor2)
  mat6(3,4) = 2*x(3)**2*(fi(1)*d12_dx4+fi(2)*d22_dx4+fi(3)*d32_dx4+d(1,2)*dfi1_dx4+d(2,2)*dfi2_dx4+d(3,2)*dfi3_dx4)/(factor1*factor2)
  mat6(3,5) = 2*x(3)**2*(fi(1)*d12_dx5+fi(2)*d22_dx5+fi(3)*d32_dx5+d(1,2)*dfi1_dx5+d(2,2)*dfi2_dx5+d(3,2)*dfi3_dx5)/(factor1*factor2)
end if
fac = 1.0_rp-2.0_rp*b*x(3)
if (fac < 0) return
x(3)=2.0_rp*x(3)/(1.0_rp+ sqrt(fac))

x(4)=x(4)-fi0*x(3)
if (logic_option(.false., make_matrix)) then
  mat6(4,2) = -fi0*mat6(3,2)-x(3)*dfi0_dx2
  mat6(4,3) = -fi0*mat6(3,3)
  mat6(4,4) = -fi0*mat6(3,4)-x(3)*dfi0_dx4+1
  mat6(4,5) = -fi0*mat6(3,5)-x(3)*dfi0_dx5
end if

b=0
do i=1,3
   b=fi(i)*d(i,1)+b
enddo

x(1)=x(1)+0.5_rp*b*x(3)**2
if (logic_option(.false., make_matrix)) then
  mat6(1,2) = x(3)*b*mat6(3,2)+0.5*x(3)**2*(fi(1)*d11_dx2+fi(2)*d21_dx2+fi(3)*d31_dx2+d(1,1)*dfi1_dx2+d(2,1)*dfi2_dx2+d(3,1)*dfi3_dx2)
  mat6(1,3) = x(3)*b*mat6(3,3)
  mat6(1,4) = x(3)*b*mat6(3,4)+0.5*x(3)**2*(fi(1)*d11_dx4+fi(2)*d21_dx4+fi(3)*d31_dx4+d(1,1)*dfi1_dx4+d(2,1)*dfi2_dx4+d(3,1)*dfi3_dx4)
  mat6(1,5) = x(3)*b*mat6(3,5)+0.5*x(3)**2*(fi(1)*d11_dx5+fi(2)*d21_dx5+fi(3)*d31_dx5+d(1,1)*dfi1_dx5+d(2,1)*dfi2_dx5+d(3,1)*dfi3_dx5)
end if

b=0
do i=1,3
   b=fi(i)*d(i,3)+b
enddo

x(6)=x(6)-0.5_rp*b*x(3)**2
if (logic_option(.false., make_matrix)) then
  mat6(6,2) = -x(3)*b*mat6(3,2)-0.5*x(3)**2*(fi(1)*d13_dx2+fi(2)*d23_dx2+fi(3)*d33_dx2+d(1,3)*dfi1_dx2+d(2,3)*dfi2_dx2+d(3,3)*dfi3_dx2)
  mat6(6,3) = -x(3)*b*mat6(3,3)
  mat6(6,4) = -x(3)*b*mat6(3,4)-0.5*x(3)**2*(fi(1)*d13_dx4+fi(2)*d23_dx4+fi(3)*d33_dx4+d(1,3)*dfi1_dx4+d(2,3)*dfi2_dx4+d(3,3)*dfi3_dx4)
  mat6(6,5) = -x(3)*b*mat6(3,5)-0.5*x(3)**2*(fi(1)*d13_dx5+fi(2)*d23_dx5+fi(3)*d33_dx5+d(1,3)*dfi1_dx5+d(2,3)*dfi2_dx5+d(3,3)*dfi3_dx5)
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
real(rp) :: xn(6),pz,pt
real(rp) :: a, beta0
character(20) :: r_name = 'ptc_rot_xz'
real(rp), optional :: mat6(6,6)
real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx2, dpt_dx4, dpt_dx5, arg

logical err_flag
logical, optional :: make_matrix

!

err_flag = .true.

arg = 1.0_rp+2.0_rp*x(5)/ beta0+x(5)**2-x(2)**2-x(4)**2
if (arg < 0) return
pz=sqrt(arg)
dpz_dx2 = -X(2)/pz
dpz_dx4 = -X(4)/pz
dpz_dx5 = (1/beta0+X(5))/pz

pt=1.0_rp-x(2)*tan(a)/pz
dpt_dx2 = -tan(a)/pz+x(2)*tan(a)*dpz_dx2/pz**2
dpt_dx4 = x(2)*tan(a)*dpz_dx4/pz**2
dpt_dx5 = x(2)*tan(a)*dpz_dx5/pz**2

xn(1)=x(1)/cos(a)/pt
xn(2)=x(2)*cos(a)+sin(a)*pz
xn(3)=x(3)+x(4)*x(1)*tan(a)/pz/pt
xn(6)=x(6)+x(1)*tan(a)/pz/pt*(1.0_rp/ beta0+x(5))

if (logic_option(.false., make_matrix)) then 
  call mat_make_unit(mat6)
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

x(1)=xn(1)
x(2)=xn(2)
x(3)=xn(3)
x(6)=xn(6)

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
! Module needed:
!   use track1_mod
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
character(20) :: r_name = 'exact_bend_edge_kick'

! Get reference beta0

beta0 = ele%value(p0c$) / ele%value(e_tot$)
if (ele%is_on) then
  g_tot = ele%value(g$) + ele%value(g_err$)
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
 
if (physical_ele_end(particle_at, orb%direction, ele%orientation) == entrance_end$) then
  edge_angle = ele%value(e1$)
  fint = ele%value(FINT$)
  hgap = ele%value(HGAP$)
else
  edge_angle = ele%value(e2$)
  fint = ele%value(FINTX$)
  hgap = ele%value(HGAPX$)
endif

! Save time

ct = X(6)

if (particle_at == first_track_edge$) then
  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int,mat6)

  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

else if (particle_at == second_track_edge$) then
  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
    return
  endif
  if (logic_option(.false., make_matrix)) mat6 = matmul(mat6_int, mat6)

  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X, err_flag, mat6_int, make_matrix)
  if (err_flag) then
    orb%state = lost_z_aperture$
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

end subroutine exact_bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine track1_low_energy_z_correction (orbit, ele, param)
! 
! Routine to add a correction to z due to speed < c corrections when tracking through an element.
! This routine assumes a constant velocity.
!
! Moudle needed:
!   use track1_mod
!
! Input:
!   orbit   -- coord_struct: Position before correction
!   ele     -- ele_struct: Element being tracked through
!   param   -- lat_param_struct: Species info.
!
! Output:
!   orbit   -- coord_struct: Position after correction.
!-

subroutine track1_low_energy_z_correction (orbit, ele, param)

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (lat_param_struct) param

real(rp) p0c, pc, beta0, mass, e_tot

!

mass = mass_of(orbit%species)
e_tot = ele%value(e_tot$)
p0c = ele%value(p0c$)

if (abs(orbit%vec(6)) < 1d-6 * mass**2 * p0c / e_tot**3) then
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * orbit%vec(6) * (1 - 3 * orbit%vec(6) / 2) * (mass / e_tot)**2
else
  pc = (1 + orbit%vec(6)) * ele%value(p0c$)
  beta0 = ele%value(p0c$) / ele%value(e_tot$)
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * (orbit%beta - beta0) / beta0
endif

end subroutine track1_low_energy_z_correction

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit)
! 
! Routine to add a RF cavity coupler kicks
!
! Moudle needed:
!   use track1_mod
!
! Input:
!   ele         -- ele_struct: Element being tracked through
!   param       -- lat_param_struct: branch parameters.
!   particle_at -- integer: first_track_edge$, or second_track_edge$.
!   phase       -- real(rp): phase of cavity
!   orbit       -- coord_struct: Position before kick.
!
! Output:
!   orbit       -- coord_struct: Position after kick.
!-

subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit

real(rp) dp, dp_x, dp_y, phase, ph, dE

integer particle_at, physical_end

!

physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end (physical_end, nint(ele%value(coupler_at$)))) return
if (ele%value(coupler_strength$) == 0) return

ph = phase
if (ele%key == rfcavity$) ph = pi/2 - ph
ph = ph + twopi * ele%value(coupler_phase$)

dp = e_accel_field(ele, gradient$) * ele%value(coupler_strength$)

dp_x = dp * cos(twopi * ele%value(coupler_angle$))
dp_y = dp * sin(twopi * ele%value(coupler_angle$))

if (nint(ele%value(coupler_at$)) == both_ends$) then
  dp_x = dp_x / 2
  dp_y = dp_y / 2
endif

! Since the transverse kick is time dependent there must be an associated energy kick
! to make things symplectic.

orbit%vec(2) = orbit%vec(2) + dp_x * cos(ph) / orbit%p0c
orbit%vec(4) = orbit%vec(4) + dp_y * cos(ph) / orbit%p0c

dE = (dp_x * orbit%vec(1) + dp_y * orbit%vec(3)) * sin(ph) * twopi * ele%value(rf_frequency$) / c_light
call apply_energy_kick (dE, orbit)

end subroutine rf_coupler_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, w_mat, p_s, track_spin)
! 
! Routine to track through a patch element.
! The steps for tracking are:
!   1) Transform from entrance to exit coordinates.
!   2) Drift particle from the entrance to the exit coordinants.
!
! Moudle needed:
!   use track1_mod
!
! Input:
!   ele           -- ele_struct: patch element.
!   orbit         -- coord_struct: Starting phase space coords
!   drift_to_exit -- logical, optional: If False then do not drift the particle from
!                      beginning to end face. Also do not correct for a reference energy shift.
!                      Default is True. 
!   track_spin    -- logical, optional: If True rotate the spin vector appropriately. 
!                       If ele%spin_tracking_method = symp_lie_ptc -> default = True. Else default = False.
!
! Output:
!   orbit      -- coord_struct: Coords after applying a patch transformation.
!   s_ent      -- real(rp), optional: Longitudinal coordinate of initial particle position at
!                   the entrance face in the frame of reference of the exit face.
!                   For a patch with positive z_offset and all other attributes zero, s_ent = -z_offset.
!   ds_ref     -- real(rp), optional: Distance the reference particle travels from the entrance
!                   to the exit face.
!   w_mat(3,3) -- real(rp), optional: Rotation matrix used in tracking.
!   p_s        -- real(rp), optional: Longitudinal momentum.
!-

subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, w_mat, p_s, track_spin)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit

real(rp), pointer :: v(:)
real(rp) p_vec(3), r_vec(3), rel_pc, ww(3,3), beta0, ds0, s_vec(3)
real(rp), optional :: s_ent, ds_ref, w_mat(3,3), p_s

logical, optional :: drift_to_exit, track_spin


! Transform to exit face coords.

v => ele%value

rel_pc = 1 + orbit%vec(6)
p_vec = [orbit%vec(2), orbit%vec(4), orbit%direction * sqrt(rel_pc**2 - orbit%vec(2)**2 - orbit%vec(4)**2)]

! For other types of elements, the ele%orientation is the same as the upstream and downstream 
! elements orientation. For a patch this is not necessarily true which is why a patch element
! needs to store the upstream and downstream orientations.

! orbit%direction * ele%orientation = -1 means we are going from the downstream end to the upstream end.
! In this case, must reverse the order of application of the offsets and pitches.

if (orbit%direction * ele%orientation == 1) then
  p_vec(3) = p_vec(3) * ele%value(upstream_ele_dir$)
  r_vec = [orbit%vec(1) - v(x_offset$), orbit%vec(3) - v(y_offset$), -v(z_offset$)]
  if (v(x_pitch$) /= 0 .or. v(y_pitch$) /= 0 .or. v(tilt$) /= 0) then
    call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat_inv = ww)
    p_vec = matmul(ww, p_vec)
    r_vec = matmul(ww, r_vec)
    orbit%vec(2) = p_vec(1)
    orbit%vec(4) = p_vec(2)
    if (logic_option((ele%spin_tracking_method /= symp_lie_ptc$), track_spin)) then
      s_vec = spinor_to_vec(orbit%spin)
      orbit%spin = vec_to_spinor(matmul(ww, s_vec))
    endif
  else
    call mat_make_unit (ww)
  endif
  ds0 = ww(3,1) * v(x_offset$) + ww(3,2) * v(y_offset$) + ww(3,3) * v(z_offset$)
  orbit%vec(5) = orbit%vec(5) + orbit%beta * c_light * v(t_offset$)

else
  p_vec(3) = p_vec(3) * ele%value(downstream_ele_dir$)
  r_vec = [orbit%vec(1), orbit%vec(3), 0.0_rp]
  if (v(x_pitch$) /= 0 .or. v(y_pitch$) /= 0 .or. v(tilt$) /= 0) then
    call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat = ww)
    p_vec = matmul(ww, p_vec)
    r_vec = matmul(ww, r_vec)
    orbit%vec(2) = p_vec(1)
    orbit%vec(4) = p_vec(2)
    if (logic_option((ele%spin_tracking_method /= symp_lie_ptc$), track_spin)) then
      s_vec = spinor_to_vec(orbit%spin)
      orbit%spin = vec_to_spinor(matmul(ww, s_vec))
    endif
  else
    call mat_make_unit (ww)
  endif
  r_vec = r_vec + [v(x_offset$), v(y_offset$), v(z_offset$)]
  ds0 = ww(1,3) * v(x_offset$) + ww(2,3) * v(y_offset$) + ww(3,3) * v(z_offset$)
  orbit%vec(5) = orbit%vec(5) - orbit%beta * c_light * v(t_offset$)
endif

!

orbit%vec(1) = r_vec(1)
orbit%vec(3) = r_vec(2)

if (present(s_ent))  s_ent = r_vec(3)
if (present(ds_ref)) ds_ref = ds0
if (present(w_mat))  w_mat = ww
if (present(p_s))    p_s = p_vec(3)

! Drift to exit face.
! Notice that the drift distance is -r_vec(3). 

if (logic_option(.true., drift_to_exit)) then
  ! Set track edge so that energy correction does not ignore an energy shift.
  if (orbit%direction == 1) then
    call reference_energy_correction (ele, orbit, first_track_edge$)
  else
    call reference_energy_correction (ele, orbit, second_track_edge$)
  endif
  beta0 = v(p0c$) / v(e_tot$)
  orbit%vec(1) = orbit%vec(1) - r_vec(3) * p_vec(1) / p_vec(3)
  orbit%vec(3) = orbit%vec(3) - r_vec(3) * p_vec(2) / p_vec(3)
  orbit%vec(5) = orbit%vec(5) + r_vec(3) * rel_pc / p_vec(3) + ds0 * orbit%beta / beta0
  orbit%t = orbit%t - r_vec(3) * rel_pc / (p_vec(3) * orbit%beta * c_light)
  orbit%s = orbit%s + ds0
endif

end subroutine track_a_patch

end module
