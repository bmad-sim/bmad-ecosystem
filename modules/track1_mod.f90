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
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. 
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

real(rp) x_lim, y_lim, x_particle, y_particle, s_here, r, rel_p
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
    if (.not. (lord_edge_aligned (ele, physical_end, lord) .or. lord%aperture_at == continuous$)) cycle
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

if (x_particle < 0) then
  x_lim = ele%value(x1_limit$)
else
  x_lim = ele%value(x2_limit$)
endif

if (x_lim <= 0 .or. .not. param%aperture_limit_on) x_lim = bmad_com%max_aperture_limit

if (y_particle < 0) then
  y_lim = ele%value(y1_limit$)
else
  y_lim = ele%value(y2_limit$)
endif

if (y_lim <= 0 .or. .not. param%aperture_limit_on) y_lim = bmad_com%max_aperture_limit

if (x_lim == 0 .and. y_lim == 0) return

select case (ele%aperture_type)

case (elliptical$)
  if (x_lim == 0 .or. y_lim == 0) then
    call out_io (s_fatal$, r_name, &
              'ECOLLIMATOR HAS ONE LIMIT ZERO AND THE OTHER NOT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  endif

  r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
  if (r > 1) then
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
  endif

case (rectangular$, auto_aperture$)

  if (abs(x_particle) > x_lim .or. abs(y_particle) > y_lim) then
    if (abs(x_particle)/x_lim > abs(y_particle)/y_lim) then
      if (x_particle > 0) then; orb%state = lost_pos_x_aperture$
      else;                     orb%state = lost_neg_x_aperture$
      endif
      param%unstable_factor = abs(x_particle) / x_lim - 1
    else
      if (y_particle > 0) then; orb%state = lost_pos_y_aperture$
      else;                     orb%state = lost_neg_y_aperture$
      endif
      param%unstable_factor = abs(y_particle) / y_lim - 1
    endif
  endif

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
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
!   px^2 + py^2 > 1 + pz (non-photons only)
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
! Subroutine track_a_drift (orb, ele, length)
!
! Subroutine to track a particle as through a drift.
!
! Modules needed:
!   use track1_mod
!
! Input:
!   orb      -- coord_struct: Orbit at start of the drift.
!   ele      -- Ele_struct: Element tracked through.
!   length   -- Real(rp): Length to drift through.
!
! Output:
!   orb      -- coord_struct: Orbit at end of the drift
!-

subroutine track_a_drift (orb, ele, length)

implicit none

type (coord_struct) orb
type (ele_struct) ele
type (lat_param_struct) param
real(rp) length, rel_pc, dz, px, py, pz, pxy2

! Everything but photons

rel_pc = 1 + orb%vec(6)
px = orb%vec(2) / rel_pc
py = orb%vec(4) / rel_pc
pxy2 = px**2 + py**2
if (pxy2 >= 1) then
  orb%state = lost_z_aperture$
  return
endif
pz = sqrt(1 - pxy2)

orb%vec(1) = orb%vec(1) + length * px / pz
orb%vec(3) = orb%vec(3) + length * py / pz

if (orb%beta > 0) then
  dz = length * (orb%beta * ele%value(e_tot$) / ele%value(p0c$) - 1/pz)
  orb%t = orb%t + length / (orb%beta * pz * c_light)
else
  dz = length * (1 - 1/pz)
endif

orb%vec(5) = orb%vec(5) + dz
orb%s = orb%s + length

end subroutine track_a_drift

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

c_dir = ele%orientation * end_orb%direction * param%rel_tracking_charge
call bend_edge_kick (end_orb, ele, param, first_track_edge$)

! If we have a sextupole component then step through in steps of length ds_step

n_step = 1

call multipole_ele_to_kt(ele, param, .false., has_nonzero_pole, knl, tilt)
if (ele%value(k2$) /= 0 .or. has_nonzero_pole) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)
if (has_nonzero_pole) knl = knl / n_step

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
    call track_a_drift (end_orb, ele, length)
    drifting = .true.

  !-----------------------------------------------------------------------
  ! Track through main body...
  ! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  else

    ct = cos(angle)
    st = sin(angle)
    if (angle < 1d-7) then
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
    if (rel_p2 - pxy2 < 0.1) then  ! somewhat arbitrary cutoff
      end_orb%state = lost$
      end_orb%vec(1) = 2 * bmad_com%max_aperture_limit
      end_orb%vec(3) = 2 * bmad_com%max_aperture_limit
      return
    endif 

    p_long = sqrt(rel_p2 - pxy2)

    ! The following is to make sure that a beam entering on-axis remains 
    ! *exactly* on-axis.

    if (pxy2 < 1e-5) then  
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

    if (abs(angle) < 1e-5 .and. abs(g_tot * length) < 1e-5) then
      end_orb%vec(1) = end_orb%vec(1) + length * px / p_long - &
                       g_tot * (length * Dy)**2 / (2 * p_long**3) + &
                       g * length * (length * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) / (2 * p_long**2)
    elseif (abs(g_tot) < 1e-5 * abs(g)) then
      alpha = p_long * ct - px * st
      end_orb%vec(1) = (p_long * (1 + g * x) - alpha) / (g * alpha) - &
                   g_tot * (Dy * (1 + g * x) * st)**2 / (2 * alpha**3 * g**2) + &
                   g_tot**2 * Dy**2 * ((1 + g * x) * st)**3 * (px * ct + p_long * st) / (2 * alpha**5 * g**3)
    else
      eps = px_t**2 + py**2
      if (eps < 1e-5 * rel_p2) then  ! use small angle approximation
        eps = eps / (2 * rel_p)
        end_orb%vec(1) = (pz + px * st - ct * p_long + g_tot * x * ct + ct - g_err * one_ct + &
                                                                 eps * (eps / (2 * rel_p) - 1)) / g_tot
      else
        end_orb%vec(1) = (sqrt(rel_p2 - eps) + px*st + g_tot*x*ct - p_long*ct) / g_tot - one_ct
      endif
    endif

    end_orb%vec(2) = px_t
    end_orb%vec(4) = py

    if (abs(g_tot) < 1e-5 * abs(g)) then
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
call bend_edge_kick (end_orb, ele, param, second_track_edge$)

call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false.)

if (.not. drifting) call track1_low_energy_z_correction (end_orb, ele, param)

end subroutine track_a_bend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine edge_focus_only (orb, ele, param, particle_at, mat6)
!
! Subroutine to track through the edge field of an sbend.
! Apply only the first order kick, which is edge focusing.
! This routine is the same as approx_bend_edge_kick, except that only
! first order terms are applied.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   particle_at -- Integer: first_track_edge$, or second_track_edge$, 
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix.
!-

subroutine edge_focus_only (orb, ele, param, particle_at, mat6)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) e, g_tot, fint, hgap, ht_x, ht_y, cos_e, sin_e, tan_e, v0(6)
real(rp) c_dir
real(rp) ht2, hs2, sec_e, k1_eff, k1
integer particle_at, element_end

character(24), parameter :: r_name = 'edge_focus_only'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a g_err.
! Apply only the first order kick.  i.e. only edge focusing.

c_dir = param%rel_tracking_charge * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(g_err$)) * c_dir
  k1 = ele%value(k1$)
else
  g_tot = 0
  k1 = 0
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
k1_eff = k1 * c_dir

if (fint == 0) then
  ht_y = -ht_x
else
  ht_y = -g_tot * tan(e - 2 * fint * g_tot * hgap * (1 + sin_e**2) / cos_e)
endif

v0 = orb%vec

if (present(mat6)) call mat_make_unit(mat6)

if (particle_at == first_track_edge$) then
  orb%vec(1) = v0(1)
  orb%vec(2) = v0(2) + ht_x * v0(1)
  orb%vec(3) = v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3)
  if (present(mat6)) then
    mat6(1,1) = 1
    mat6(2,1) = ht_x
    mat6(2,2) = 1
    mat6(3,3) = 1
    mat6(4,3) = ht_y
    mat6(4,4) = 1
  end if
else
  orb%vec(1) = v0(1)
  orb%vec(2) = v0(2) + ht_x * v0(1)
  orb%vec(3) = v0(3)
  orb%vec(4) = v0(4) + ht_y * v0(3)
  if (present(mat6)) then
    mat6(1,1) = 1
    mat6(2,1) = ht_x
    mat6(2,2) = 1
    mat6(3,3) = 1
    mat6(4,3) = ht_y
    mat6(4,4) = 1
  end if
endif

end subroutine edge_focus_only

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine approx_bend_edge_kick (orb, ele, param, particle_at, mat6)
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
!   particle_at -- Integer: first_track_edge$, or second_track_edge$, 
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix.
!-

subroutine approx_bend_edge_kick (orb, ele, param, particle_at, mat6)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) e, g_tot, fint, hgap, ht_x, ht_y, cos_e, sin_e, tan_e, sec_e, v0(6), k1_eff
real(rp) ht2, hs2, c_dir, k1
integer particle_at, element_end

character(24), parameter :: r_name = 'approx_bend_edge_kick'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a g_err.

c_dir = param%rel_tracking_charge * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)

if (ele%is_on) then
  g_tot = (ele%value(g$) + ele%value(g_err$)) * c_dir
  k1 = ele%value(k1$)
else
  g_tot = 0
  k1 = 0
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
k1_eff = k1 * c_dir

if (fint == 0) then
  ht_y = -ht_x
else
  ht_y = -g_tot * tan(e - 2 * fint * g_tot * hgap * (1 + sin_e**2) / cos_e)
endif

v0 = orb%vec

if (present(mat6)) call mat_make_unit(mat6)

if (particle_at == first_track_edge$) then
  orb%vec(1) = v0(1) - ht2 * v0(1)**2 / 2 + hs2 * v0(3)**2 / 2
  orb%vec(2) = v0(2) + ht_x * v0(1) + ht2 * (v0(1) * v0(2) - v0(3) * v0(4)) + &
                       k1_eff * tan_e * (v0(1)**2 - v0(3)**2) + &
                       ht_x * (ht2 + hs2) * v0(3)**2 / 2
  orb%vec(3) = v0(3) + ht2 * v0(1) * v0(3) 
  orb%vec(4) = v0(4) + ht_y * v0(3) - ht2 * v0(1) * v0(4) - hs2 * v0(2) * v0(3) - &
                       2 * k1_eff * tan_e * v0(1) * v0(3)
  if (present(mat6)) then
    mat6(1,1) = 1 - ht2 * v0(1)
    mat6(1,3) = hs2 * v0(3)
    mat6(2,1) = ht_x + ht2 * v0(2) + 2 * k1_eff * tan_e * v0(1)
    mat6(2,2) = 1 + ht2 * v0(1)
    mat6(2,3) = - ht2 * v0(4) - 2 * k1_eff * tan_e * v0(3) + ht_x * (ht2 + hs2) * v0(3)
    mat6(2,4) = - ht2 * v0(3)
    mat6(3,1) = ht2 * v0(3)
    mat6(3,3) = 1 +  ht2 * v0(1)
    mat6(4,1) = - ht2 * v0(4) - 2 * tan_e * k1_eff * v0(3)
    mat6(4,2) = - hs2 * v0(3)
    mat6(4,3) = ht_y - hs2 * v0(2) - 2 * tan_e * k1_eff * v0(1)
    mat6(4,4) = 1 - ht2 * v0(1)
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
  if (present(mat6)) then
    mat6(1,1) = 1 + ht2 * v0(1)
    mat6(1,3) = - hs2 * v0(3)
    mat6(2,1) = ht_x - ht2 * v0(2) + 2 * k1_eff * tan_e * v0(1) - ht_x * ht2 * v0(1)
    mat6(2,2) = 1 - ht2 * v0(1)
    mat6(2,3) = ht2 * v0(4) - 2 * k1_eff * tan_e * v0(3) - ht_x * ht2 * v0(3)
    mat6(2,4) = ht2 * v0(3)
    mat6(3,1) = - ht2 * v0(3)
    mat6(3,3) = 1 -  ht2 * v0(1)
    mat6(4,1) = ht2 * v0(4) + (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(3)
    mat6(4,2) = hs2 * v0(3)
    mat6(4,3) = ht_y + hs2 * v0(2) + (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(1)
    mat6(4,4) = 1 + ht2 * v0(1)
  end if
endif

end subroutine approx_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine apply_hard_edge_kick (orb, s_edge, t_rel, hard_ele, track_ele, param, particle_at)
!
! Subroutine to track through the edge field of an element.
! This routine is used with bmad_standard field_calc where the field
! can have an abrubt, unphysical termination of the longitudinal field at the
! edges of the element. 
!
! Elements that have kicks due to unphysical edge field terminations:
!   sbend
!   solenoid
!   sol_quad
!   sad_mult
!   bend_sol_quad
!   lcavity
!   rfcavity 
!   e_gun
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords in element reference frame.
!   s_edge      -- real(rp): Hard edge relative to start of hard_ele.
!   t_rel       -- real(rp): Time relative to track_ele entrance edge
!   hard_ele    -- ele_struct: Element with hard edges.
!   track_ele   -- ele_struct: Element being tracked through. 
!                    Is different from hard_ele when there are superpositions.
!   param       -- lat_param_struct: lattice parameters.
!   particle_at -- Integer: first_track_edge$ or second_track_edge$
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!-

subroutine apply_hard_edge_kick (orb, s_edge, t_rel, hard_ele, track_ele, param, particle_at)

implicit none

type (ele_struct) hard_ele, track_ele
type (coord_struct) orb
type (lat_param_struct) param
type (em_field_struct) field

real(rp) t, f, l_drift, ks, t_rel, s_edge, s

integer particle_at, physical_end

! 

if (hard_ele%field_calc /= bmad_standard$) return
physical_end = physical_ele_end (particle_at, orb%direction, track_ele%orientation)

select case (hard_ele%key)
case (quadrupole$, sad_mult$)
  call quadrupole_edge_kick (hard_ele, particle_at, orb)

case (sbend$)
  call bend_edge_kick (orb, hard_ele, param, particle_at)

! Note: Cannot trust hard_ele%value(ks$) here since element may be superimposed with an lcavity.
! So use hard_ele%value(bs_field$).

case (solenoid$, sol_quad$, bend_sol_quad$)
  ks = param%rel_tracking_charge * hard_ele%value(bs_field$) * c_light / orb%p0c
  if (particle_at == first_track_edge$) then
    orb%vec(2) = orb%vec(2) + ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) - ks * orb%vec(1) / 2
  else
    orb%vec(2) = orb%vec(2) - ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) + ks * orb%vec(1) / 2
  endif

case (lcavity$, rfcavity$, e_gun$)

  ! Add on bmad_com%significant_length to make sure we are just inside the cavity.
  f = charge_of(orb%species) / (2 * orb%p0c)
  t = t_rel + track_ele%value(ref_time_start$) - hard_ele%value(ref_time_start$) 
  s = s_edge

  if (at_this_ele_end(physical_end, nint(track_ele%value(fringe_at$)))) then

    if (particle_at == first_track_edge$) then
      ! E_gun does not have an entrance kick
      s = s + bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, t, orb, .true., field)

      orb%vec(2) = orb%vec(2) - field%e(3) * orb%vec(1) * f + c_light * field%b(3) * orb%vec(3) * f
      orb%vec(4) = orb%vec(4) - field%e(3) * orb%vec(3) * f - c_light * field%b(3) * orb%vec(1) * f

    else
      s = s - bmad_com%significant_length / 10 ! Make sure inside field region
      call em_field_calc (hard_ele, param, s, t, orb, .true., field)

      orb%vec(2) = orb%vec(2) + field%e(3) * orb%vec(1) * f - c_light * field%b(3) * orb%vec(3) * f
      orb%vec(4) = orb%vec(4) + field%e(3) * orb%vec(3) * f + c_light * field%b(3) * orb%vec(1) * f
    endif

    ! orb%phase(1) is set by em_field_calc.

    call rf_coupler_kick (hard_ele, param, particle_at, orb%phase(1), orb)

  endif

end select

end subroutine apply_hard_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine quadrupole_edge_kick (ele, particle_at, orbit)
!
! Routine to add the 3rd order quadrupolar edge kick.
! This routine assumes that the particle orbit has been rotated to the element reference frame.
!
! Moudle needed:
!   use track1_mod
!
! Input:
!   ele         -- ele_struct: Element being tracked through
!   particle_at -- integer: first_track_edge$, or second_track_edge$.
!   orbit       -- coord_struct: Position before kick.
!
! Output:
!   orbit   -- coord_struct: Position after kick.
!-

subroutine quadrupole_edge_kick (ele, particle_at, orbit)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: m_ele
type (coord_struct) orbit

real(rp) k1_rel, x, y, px, py, charge_dir
real(rp) f1, f2, ef1, vec(4), rel_pc, vx, vy

integer particle_at
integer fringe_at, physical_end, fringe_type

!

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type == none$) return

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

charge_dir = ele%orientation * orbit%direction
if (associated(ele%branch)) charge_dir = charge_dir * ele%branch%param%rel_tracking_charge

rel_pc = 1 + orbit%vec(6)

select case (ele%key)
case (quadrupole$)
  k1_rel = charge_dir * ele%value(k1$) / rel_pc
case (sad_mult$)
  ! Slice slaves and super slaves have their associated multipoles stored in the lord
  if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
    m_ele => pointer_to_lord(ele, 1)
  else
    m_ele => ele
  endif
  k1_rel = charge_dir * sqrt(m_ele%a_pole(1)**2 + m_ele%b_pole(1)**2)  / m_ele%value(l$) / rel_pc
end select

! Everything but SAD nonlinear

select case (fringe_type)
case (full_straight$, full_bend$)
  if (particle_at == second_track_edge$) k1_rel = -k1_rel

  x = orbit%vec(1); px = orbit%vec(2); y = orbit%vec(3); py = orbit%vec(4)
  orbit%vec(1) = x  + k1_rel * (x**3/12 + x*y**2/4)
  orbit%vec(2) = px + k1_rel * (x*y*py/2 - px*(x**2 + y**2)/4)
  orbit%vec(3) = y  - k1_rel * (y**3/12 + y*x**2/4)
  orbit%vec(4) = py - k1_rel * (y*x*px/2 - py*(y**2 + x**2)/4)
  orbit%vec(5) = orbit%vec(5) + k1_rel * (y**3*py/12 - x**3*px/12 + x**2*y*py/4 - x*y**2*px/4) / (1 + orbit%vec(6))

case (sad_linear$, sad_full$)
  f1 = -k1_rel * ele%value(f1$) * abs(ele%value(f1$)) / 24
  f2 =  k1_rel * ele%value(f2$)
  if (f1 /= 0 .or. f2 /= 0) then
    if (particle_at == second_track_edge$) f1 = -f1

    ef1 = exp(f1)

    vec = orbit%vec(1:4)
    vx = vec(2) / rel_pc;  vy = vec(4) / rel_pc

    orbit%vec(5) = orbit%vec(5) - (f1 * vec(1) + f2 * (1 + f1/2) * vx / ef1) * vx + &
                                  (f1 * vec(3) + f2 * (1 - f1/2) * vy * ef1) * vy

    orbit%vec(1:2) = [vec(1) * ef1 + vx * f2,  vec(2) / ef1]
    orbit%vec(3:4) = [vec(3) / ef1 - vy * f2,  vec(4) * ef1]
  endif
end select

! SAD nonlinear
! Not yet implemented

select case (fringe_type)
case (sad_nonlin_only$, sad_full$)
  
end select

end subroutine quadrupole_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bend_edge_kick (orb, ele, param, particle_at, mat6)
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
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix.
!-

subroutine bend_edge_kick (orb, ele, param, particle_at, mat6)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
integer particle_at, fringe_type, physical_end

character(*), parameter :: r_name = 'bend_edge_kick'

!

physical_end = physical_ele_end (particle_at, orb%direction, ele%orientation)
if (.not. at_this_ele_end (physical_end, nint(ele%value(fringe_at$)))) then
  if (present(mat6)) call mat_make_unit(mat6)
  return
endif

! Fringe due to finite e1 or e2

fringe_type = nint(ele%value(fringe_type$))
select case (fringe_type)
case (full_straight$, full_bend$)
  call exact_bend_edge_kick (orb, ele, param, particle_at, mat6)
case (basic_bend$, sad_full$, sad_linear$, sad_nonlin_only$)
  call approx_bend_edge_kick (orb, ele, param, particle_at, mat6)
case (edge_focus_only$)
  call edge_focus_only (orb, ele, param, particle_at, mat6)
case (none$)
  if (present(mat6)) call mat_make_unit (mat6)
case default
  call out_io (s_fatal$, r_name, 'UNKNOWN FRINGE_TYPE: \i0\ ', i_array = [fringe_type])
  if (global_com%exit_on_error) call err_exit
end select

! Sad fringe

select case (fringe_type)
case (sad_full$, sad_linear$) 
  call add_sad_bend_edge_kick (orb, ele, param, particle_at, mat6)
end select

end subroutine bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine add_sad_bend_edge_kick (orb, ele, param, particle_at, mat6)
!
! Subroutine to track through the edge field of an sbend.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   mat6        -- real(rp), optional: Starting matrix 
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix after fringe field
!-

subroutine add_sad_bend_edge_kick (orb, ele, param, particle_at, mat6)

use ptc_interface_mod

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) :: sad_mat(6,6)
real(rp) :: f1, el_p, g, ct, c1, c2, c3, y, px, rel_p, sin_e, e

integer :: particle_at, c_dir, element_end

character(*), parameter :: r_name = 'add_sad_bend_edge_kick'

! Finite f1 fringe

f1 = ele%value(f1$)
if (f1 == 0) return

c_dir = param%rel_tracking_charge * ele%orientation * orb%direction
element_end = physical_ele_end(particle_at, orb%direction, ele%orientation)

if (element_end == entrance_end$) then
  e = sin(ele%value(e1$))
else
  e = sin(ele%value(e2$))
endif

g = ele%value(g$) + ele%value(g_err$)
sin_e = sin(e)
if (particle_at == second_track_edge$) then
  sin_e = -sin_e
  g = -g
endif

px = orb%vec(2) + sin_e
y  = orb%vec(3)
rel_p = 1 + orb%vec(6)

c1 = f1**2 * g / (24 * rel_p)  ! * px
c2 = f1 * g**2 / (12 * rel_p)  ! * y^2
c3 = g**2 / (6 * f1 * rel_p)   ! * y^4

if (present(mat6)) then
  call mat_make_unit (sad_mat)
  sad_mat(1,6) =  c1 / rel_p
  sad_mat(4,3) =  2 * c2 - 12 * c3 * y**2
  sad_mat(4,6) = (-2 * c2 * y + 4 * c3 * y**3) / rel_p
  sad_mat(5,2) = c1 / rel_p
  sad_mat(5,3) = (2 * c2 * y - 4 * c3 * y**3) / rel_p
  sad_mat(5,6) = -2 * (c1 * px + c2 * y**2 - c3 * y**4) / rel_p**2 + c1

  mat6(1,:) = mat6(1,:) + sad_mat(1,6) * mat6(6,:)
  mat6(4,:) = mat6(4,:) + sad_mat(4,3) * mat6(3,:) + sad_mat(4,6) * mat6(6,:)
  mat6(5,:) = mat6(5,:) + sad_mat(5,2) * mat6(2,:) + sad_mat(5,3) * mat6(3,:) + sad_mat(5,6) * mat6(6,:)
endif

orb%vec(1) = orb%vec(1) + c1 * orb%vec(6)
orb%vec(4) = orb%vec(4) + 2 * c2 * y - 4 * c3 * y**3
orb%vec(5) = orb%vec(5) + (c1 * px + c2 * y**2 - c3 * y**4) / rel_p

end subroutine add_sad_bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! THIS ROUTINE IS SLATED FOR REMOVAL !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+
! Subroutine add_sad_bend_hard_edge_edge_kick (orb, ele, param, particle_at, mat6)
!
! Subroutine to track through the edge field of an sbend.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!   mat6        -- real(rp), optional: Starting matrix 
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix after fringe field
!-

subroutine add_sad_bend_hard_edge_edge_kick (orb, ele, param, particle_at, mat6)

use ptc_interface_mod

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) :: sad_mat(6,6), dx, dpy, dz, sqt
real(rp) :: f1, el_p, g, ct, c1, c2, c3, y, px, rel_p, p_long, yyy

integer :: particle_at

character(*), parameter :: r_name = 'add_sad_bend_hard_edge_edge_kick'

!

g = ele%value(g$)
if (particle_at == second_track_edge$) g = -g

px = orb%vec(2)
y  = orb%vec(3)
rel_p = 1 + orb%vec(6)

c1 = g**2 / 12
yyy = y**2 * (1 - c1 * y**2)

p_long = rel_p**2 - px**2
sqt = sqrt(p_long)
dx = rel_p**2 * g * yyy / (2 * sqt**3)
dpy = -px * g * y * (1 - 2 * c1 * y**2) / sqt
dz = -rel_p * g * px * yyy / (2 * sqt**3)

if (present(mat6)) then
  call mat_make_unit (sad_mat)
  sad_mat(1,2) = dx * 3 * px / p_long
  sad_mat(1,3) = rel_p**2 * g * (y - 2 * c1 * y**3) / sqt**3
  sad_mat(1,6) = -dx * (rel_p**2 + 2 * px**2) / (rel_p * p_long)
  sad_mat(4,2) = -rel_p**2 * g * y * (1 - 2 * c1 * y**2) / sqt**3
  sad_mat(4,3) = -px * g * (1 - 6 * c1 * y**2) / sqt
  sad_mat(4,6) = -dpy * rel_p / p_long 
  sad_mat(5,2) = -rel_p * g * (rel_p**2 + 2 * px**2) * yyy / (2 * sqt**5)
  sad_mat(5,3) = -rel_p * g * px * (y - 2 * c1 * y**3) / sqt**3
  sad_mat(5,6) = g * px * (2 * rel_p**2 + px**2) * yyy / (2 * sqt**5)

  mat6(1,:) = mat6(1,:) + sad_mat(1,2) * mat6(2,:) + sad_mat(1,3) * mat6(3,:) + sad_mat(1,6) * mat6(6,:)
  mat6(4,:) = mat6(4,:) + sad_mat(4,2) * mat6(2,:) + sad_mat(4,3) * mat6(3,:) + sad_mat(4,6) * mat6(6,:)
  mat6(5,:) = mat6(5,:) + sad_mat(5,2) * mat6(2,:) + sad_mat(5,3) * mat6(3,:) + sad_mat(5,6) * mat6(6,:)
endif

orb%vec(1) = orb%vec(1) + dx
orb%vec(4) = orb%vec(4) + dpy
orb%vec(5) = orb%vec(5) + dz

end subroutine add_sad_bend_hard_edge_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_wedger (a, g_tot, beta0, X, mat6)
!
! Subroutine to track PTC coordinates through a wedge
!
! Adapted from forest/code/Sh_def_kind.f90 : WEDGER
!
! Input:
!   a      -- real(rp): wedge angle (rad)
!   g_tot  -- real(rp): reference bending radius
!   beta0  -- real(rp): reference relativistic beta
!   X(6)   -- real(rp): PTC phase space coordinates
!
! Output:
!   X(6)   -- real(rp): PTC phase space coordinates
!   mat6   -- Real(rp), optional: Transfer matrix.
!-
subroutine ptc_wedger (a, g_tot, beta0, X, mat6)

implicit none

real(rp) :: X(6)
real(rp) :: a, beta0, g_tot
real(rp) :: Xn(6),pz,pzs,pt,b1
character(20) :: r_name = 'ptc_wedger'
real(rp), optional :: mat6(6,6)

real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx4, dpt_dx5
real(rp) dpzs_dx1, dpzs_dx2, dpzs_dx4, dpzs_dx5, factor1, factor2

! No net field case...

b1 = g_tot

if(b1==0) then
   if (present(mat6)) then
      call ptc_rot_xz(a, X, beta0, mat6)
   else
      call ptc_rot_xz(a, X, beta0)
   end if
   return
endif

! Normal case

pz=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(2)**2-X(4)**2)
dpz_dx2 = -X(2)/pz
dpz_dx4 = -X(4)/pz
dpz_dx5 = (1/beta0+X(5))/pz

if (present(mat6)) call mat_make_unit(mat6)

Xn(2)=X(2)*cos(a)+(pz-b1*X(1))*sin(a)
if (present(mat6)) then
  mat6(2,1) = -b1*sin(a)
  mat6(2,2) = cos(a)+sin(a)*dpz_dx2
  mat6(2,4) = sin(a)*dpz_dx4
  mat6(2,5) = sin(a)*dpz_dx5
end if

pt=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(4)**2)
dpt_dx4 = -X(4)/pt
dpt_dx5 = (1/beta0+X(5))/pt

pzs=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-Xn(2)**2-X(4)**2)
if (present(mat6)) then
  dpzs_dx1 = -Xn(2)*mat6(2,1)/pzs
  dpzs_dx2 = -Xn(2)*mat6(2,2)/pzs
  dpzs_dx4 = (-X(4)-Xn(2)*mat6(2,4))/pzs
  dpzs_dx5 = (1/beta0+X(5)-Xn(2)*mat6(2,5))/pzs
end if

Xn(1)=X(1)*cos(a)+(X(1)*X(2)*sin(2.0_rp*a)+sin(a)**2*(2.0_rp*X(1)*pz-b1*X(1)**2) )&
    / (pzs+pz*cos(a)-X(2)*sin(a))
if (present(mat6)) then
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
if (present(mat6)) then
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

end subroutine ptc_wedger

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, mat6)
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
!
! Output:
!   X(6)   -- real(rp): PTC phase space coordinates
!   mat6   -- Real(rp), optional: Transfer matrix.
!-

subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, mat6)

implicit none

real(rp) :: X(6) !PTC phase space coordinates
real(rp) :: FINT, HGAP
real(rp) :: beta0, g_tot
real(rp) :: PZ,XP,YP,TIME_FAC
real(rp) :: D(3,3),FI(3),FI0,B,co1,co2
integer  :: i
integer  :: particle_at
character(20) :: r_name = 'ptc_fringe_dipoler'
real(rp), optional :: mat6(6,6)
real(rp) :: dpz_dx2, dpz_dx4, dpz_dx5, dtime_fac_dx5
real(rp) :: dxp_dx2, dxp_dx4, dxp_dx5, dyp_dx2, dyp_dx4, dyp_dx5
real(rp) :: d11_dx2, d11_dx4, d11_dx5, d21_dx2, d21_dx4, d21_dx5, d31_dx2, d31_dx4, d31_dx5
real(rp) :: d12_dx2, d12_dx4, d12_dx5, d22_dx2, d22_dx4, d22_dx5, d32_dx2, d32_dx4, d32_dx5
real(rp) :: d13_dx2, d13_dx4, d13_dx5, d23_dx2, d23_dx4, d23_dx5, d33_dx2, d33_dx4, d33_dx5
real(rp) :: dco1_dx2, dco1_dx4, dco1_dx5, dco2_dx2, dco2_dx4, dco2_dx5
real(rp) :: dfi0_dx2, dfi0_dx4, dfi0_dx5, dfi1_dx2, dfi1_dx4, dfi1_dx5
real(rp) :: dfi2_dx2, dfi2_dx4, dfi2_dx5, dfi3_dx2, dfi3_dx4, dfi3_dx5
real(rp) :: factor1, factor2 

!

if (particle_at == second_track_edge$) then
   B = -g_tot  !EL%CHARGE*BN(1)
else if (particle_at == first_track_edge$) then
   B = g_tot       
else
  call out_io (s_fatal$, r_name, 'INVALID PARTICLE_AT')
  call err_exit
endif

pz=sqrt(1.0_rp+2.0_rp*x(5)/beta0+x(5)**2-x(2)**2-x(4)**2)
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

if (present(mat6)) call mat_make_unit(mat6)

b=0
do i=1,3
   b=fi(i)*d(i,2)+b
enddo

if (present(mat6)) then
  factor1 = sqrt(1-2*b*x(3))
  factor2 = (1+sqrt(1-2*b*x(3)))**2
  mat6(3,2) = 2*x(3)**2*(fi(1)*d12_dx2+fi(2)*d22_dx2+fi(3)*d32_dx2+d(1,2)*dfi1_dx2+d(2,2)*dfi2_dx2+d(3,2)*dfi3_dx2)/(factor1*factor2)
  mat6(3,3) = 2/(1+factor1)+2*x(3)*b/(factor1*factor2)
  mat6(3,4) = 2*x(3)**2*(fi(1)*d12_dx4+fi(2)*d22_dx4+fi(3)*d32_dx4+d(1,2)*dfi1_dx4+d(2,2)*dfi2_dx4+d(3,2)*dfi3_dx4)/(factor1*factor2)
  mat6(3,5) = 2*x(3)**2*(fi(1)*d12_dx5+fi(2)*d22_dx5+fi(3)*d32_dx5+d(1,2)*dfi1_dx5+d(2,2)*dfi2_dx5+d(3,2)*dfi3_dx5)/(factor1*factor2)
end if
x(3)=2.0_rp*x(3)/(1.0_rp+ sqrt(1.0_rp-2.0_rp*b*x(3)) )

x(4)=x(4)-fi0*x(3)
if (present(mat6)) then
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
if (present(mat6)) then
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
if (present(mat6)) then
  mat6(6,2) = -x(3)*b*mat6(3,2)-0.5*x(3)**2*(fi(1)*d13_dx2+fi(2)*d23_dx2+fi(3)*d33_dx2+d(1,3)*dfi1_dx2+d(2,3)*dfi2_dx2+d(3,3)*dfi3_dx2)
  mat6(6,3) = -x(3)*b*mat6(3,3)
  mat6(6,4) = -x(3)*b*mat6(3,4)-0.5*x(3)**2*(fi(1)*d13_dx4+fi(2)*d23_dx4+fi(3)*d33_dx4+d(1,3)*dfi1_dx4+d(2,3)*dfi2_dx4+d(3,3)*dfi3_dx4)
  mat6(6,5) = -x(3)*b*mat6(3,5)-0.5*x(3)**2*(fi(1)*d13_dx5+fi(2)*d23_dx5+fi(3)*d33_dx5+d(1,3)*dfi1_dx5+d(2,3)*dfi2_dx5+d(3,3)*dfi3_dx5)
end if
    
end subroutine ptc_fringe_dipoler

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_rot_xz(a, X, beta0, mat6)
!
! Subroutine to rotate the local reference frame about the Y axis in PTC coordinates.  
! Adapted from forest/code/Sc_euclidean.f90 : ROT_XZ
!
!
! Input:
!   a      -- real(rp): rotation angle (rad)
!   X(6)   -- real(rp): PTC phase space coordinates
!   beta0  -- real(rp): reference relativistic beta
!
! Output:
!   X(6)   -- real(rp): PTC phase space coordinates
!   mat6   -- Real(rp), optional: Transfer matrix.
!-

subroutine ptc_rot_xz(a, X, beta0, mat6)

implicit none

real(rp) :: x(6)
real(rp) :: xn(6),pz,pt
real(rp) :: a, beta0
character(20) :: r_name = 'ptc_rot_xz'
real(rp), optional :: mat6(6,6)
real(rp) dpz_dx2, dpz_dx4, dpz_dx5, dpt_dx2, dpt_dx4, dpt_dx5

!

pz=sqrt(1.0_rp+2.0_rp*x(5)/ beta0+x(5)**2-x(2)**2-x(4)**2)
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

if (present(mat6)) then 
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

end subroutine ptc_rot_xz

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine exact_bend_edge_kick (orb, ele, param, particle_at, mat6)
!
! Subroutine to track through the edge field of an sbend.
! Uses routines adapted from PTC
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: 
!   particle_at -- Integer: first_track_edge$, or second_track_edge$.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   mat6       -- Real(rp), optional: Transfer matrix.
!-

subroutine exact_bend_edge_kick (orb, ele, param, particle_at, mat6)

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

character(20) :: r_name = 'exact_bend_edge_kick'

! Get reference beta0

beta0 = ele%value(p0c$) / ele%value(e_tot$)
if (ele%is_on) then
  g_tot = ele%value(g$) + ele%value(g_err$)
else
  g_tot = 0
endif

! Convert to PTC coordinates
if (present(mat6)) then 
  call vec_bmad_to_ptc(orb%vec, beta0, X, mat6_int)
  mat6 = mat6_int
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
  if (present(mat6)) then
    call ptc_wedger(edge_angle, 0.0_rp, beta0, X, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_wedger(edge_angle, 0.0_rp, beta0, X)
  end if
  ! Edge kick
  if (present(mat6)) then
    call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at)
  end if
  ! Backtrack
  if (present(mat6)) then
    call ptc_wedger(-edge_angle, g_tot, beta0, X, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_wedger(-edge_angle, g_tot, beta0, X)
  end if
else if (particle_at == second_track_edge$) then
  ! Backtrack
  if (present(mat6)) then
    call ptc_wedger(-edge_angle, g_tot, beta0, X, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_wedger(-edge_angle, g_tot, beta0, X)
  end if
  ! Edge kick
  if (present(mat6)) then
    call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, particle_at)
  end if
  ! Drift forward
  if (present(mat6)) then
    call ptc_wedger(edge_angle, 0.0_rp, beta0, X, mat6_int)
    mat6 = matmul(mat6_int,mat6)
  else
    call ptc_wedger(edge_angle, 0.0_rp, beta0, X)
  end if
else
  !error!
  call err_exit
endif

! Convert back to bmad coordinates
if (present(mat6)) then
  call vec_ptc_to_bmad (X, beta0, orb%vec, mat6_int)
  mat6 = matmul(mat6_int,mat6)
else
  call vec_ptc_to_bmad (X, beta0, orb%vec)
end if   

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

real(rp) p0c, pc, beta, beta0, mass, e_tot

!

mass = mass_of(orbit%species)
e_tot = ele%value(e_tot$)
p0c = ele%value(p0c$)

if (abs(orbit%vec(6)) < 1e-6 * mass**2 * p0c / e_tot**3) then
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * orbit%vec(6) * (1 - 3 * orbit%vec(6) / 2) * (mass / e_tot)**2
else
  pc = (1 + orbit%vec(6)) * ele%value(p0c$)
  call convert_pc_to (pc, orbit%species, beta = beta)
  beta0 = ele%value(p0c$) / ele%value(e_tot$)
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * (beta - beta0) / beta0
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
! Subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, w_inv, p_s)
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
!   ele           -- ele_struct: patch, fork, or photon_fork element.
!   orbit         -- coord_struct: Starting phase space coords
!   drift_to_exit -- Logical, optional: If False then do not drift the particle from
!                      Entrance to exit faces. Also do not correct for a reference energy shift.
!                      Default is True. 
!
! Output:
!   orbit      -- coord_struct: Coords after applying a patch transformation.
!   s_ent      -- real(rp), optional: Longitudinal coordinate of initial particle position at
!                   the entrance face in the frame of reference of the exit face.
!                   For a patch with positive z_offset and all other attributes zero, s_ent = -z_offset.
!   ds_ref     -- real(rp), optional: Distance the reference particle travels from the entrance
!                   to the exit face.
!   w_inv(3,3) -- real(rp), optional: Rotation matrix used in tracking.
!   p_s        -- real(rp), optional: Longitudinal momentum.
!-

subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, w_inv, p_s)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit

real(rp), pointer :: v(:)
real(rp) p_vec(3), r_vec(3), rel_pc, winv(3,3), beta0, ds0
real(rp), optional :: s_ent, ds_ref, w_inv(3,3), p_s

integer rel_dir1

logical, optional :: drift_to_exit

! Transform to exit face coords.

v => ele%value
r_vec = [orbit%vec(1) - v(x_offset$), orbit%vec(3) - v(y_offset$), -v(z_offset$)]

rel_pc = 1 + orbit%vec(6)
p_vec = [orbit%vec(2), orbit%vec(4), orbit%direction * sqrt(rel_pc**2 - orbit%vec(2)**2 - orbit%vec(4)**2)]
p_vec(3) = p_vec(3) * ele%value(upstream_ele_dir$)

if (v(x_pitch$) /= 0 .or. v(y_pitch$) /= 0 .or. v(tilt$) /= 0) then
  call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat_inv = winv)
  if (present(w_inv)) w_inv = winv
  p_vec = matmul(winv, p_vec)
  r_vec = matmul(winv, r_vec)
  orbit%vec(2) = p_vec(1)
  orbit%vec(4) = p_vec(2)
else
  call mat_make_unit (winv)
endif

orbit%vec(1) = r_vec(1)
orbit%vec(3) = r_vec(2)
orbit%vec(5) = orbit%vec(5) + orbit%beta * c_light * v(t_offset$)

ds0 = winv(3,1) * v(x_offset$) + winv(3,2) * v(y_offset$) + winv(3,3) * v(z_offset$)

if (present(s_ent))  s_ent = r_vec(3)
if (present(ds_ref)) ds_ref = ds0
if (present(w_inv))  w_inv = winv
if (present(p_s))    p_s = p_vec(3)

! Drift to exit face.
! Notice that the drift distance is -r_vec(3). 

if (logic_option(.true., drift_to_exit)) then
  call reference_energy_correction (ele, orbit)
  beta0 = v(p0c$) / v(e_tot$)
  orbit%vec(1) = orbit%vec(1) - r_vec(3) * p_vec(1) / p_vec(3)
  orbit%vec(3) = orbit%vec(3) - r_vec(3) * p_vec(2) / p_vec(3)
  orbit%vec(5) = orbit%vec(5) + r_vec(3) * rel_pc / p_vec(3) + ds0 * orbit%beta / beta0
  orbit%t = orbit%t - r_vec(3) * rel_pc / (p_vec(3) * orbit%beta * c_light)
  orbit%s = orbit%s + ds0
endif

end subroutine track_a_patch

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine reference_energy_correction (ele, orbit)
!
! For elements where the reference energy is changing the reference energy in the body is 
! taken by convention to be the reference energy at the exit end.
! Elements where the reference energy can change:
!   lcavity
!   patch
!   custom
!
! This routine should be called at the start of any tracking integration.
!
! Input:
!   ele       -- Ele_struct: Element being tracked through.
!   orbit     -- Coord_struct: Coordinates to correct.
!
! Output:
!   orbit     -- Coord_struct: Coordinates to correct.
!-

subroutine reference_energy_correction (ele, orbit)

implicit none

type (ele_struct) :: ele
type (coord_struct) :: orbit

real(rp) p0, p1, e_start, p_rel
character(*), parameter :: r_name = 'reference_energy_correction'

!

if (ele%value(p0c$) == ele%value(p0c_start$)) return

if (orbit%direction == 1) then
  p_rel = ele%value(p0c_start$) / ele%value(p0c$)
  orbit%p0c = ele%value(p0c$)
else
  p_rel = ele%value(p0c$) / ele%value(p0c_start$)
  orbit%p0c = ele%value(p0c_start$)
endif

  orbit%vec(2) = orbit%vec(2) * p_rel
  orbit%vec(4) = orbit%vec(4) * p_rel
  orbit%vec(6) = (1 + orbit%vec(6)) * p_rel - 1

end subroutine reference_energy_correction

end module
