!+
! Module track1_mod
!
! Module of helper routines for track1 routines
!-

module track1_mod

use bmad_struct
use bmad_interface
use make_mat6_mod
use em_field_mod

! Private routines for exact_bend_edge_kick
private ptc_rot_xz, ptc_wedger, ptc_fringe_dipoler

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, at, param, check_momentum)
!
! Subroutine to check if an orbit is outside an element's aperture.
! Note: A particle will also be considered to have hit an aperture
! if |p_x| or |p_y| > 1.
!
! Modules needed:
!   use bmad
!
! Input:
!   orb   -- Coord_struct: coordinates of a particle.
!   ele   -- Ele_struct: Element holding the aperture
!     %value(x1_limit$) -- Horizontal negative side aperture.
!     %value(x2_limit$) -- Horizontal positive side aparture.
!     %value(y1_limit$) -- Vertical negative side aperture.
!     %value(y2_limit$) -- Vertical positive side aparture.
!     %offset_moves_aperture -- If True then aperture moves with the element.
!   at    -- Integer: upstream_end$ or downstream_end$
!   param -- lat_param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. 
!   check_momentum
!         -- Logical, optional: If present and false then checking of p_x and 
!               p_y will be disabled.
!
! Output:
!   orb   -- Coord_struct: coordinates of a particle.
!     %state -- State of the particle
!-

subroutine check_aperture_limit (orb, ele, at, param, check_momentum)

implicit none

type (coord_struct) :: orb
type (coord_struct) orb2 
type (ele_struct) :: ele
type (lat_param_struct), intent(inout) :: param

real(rp) x_lim, y_lim, x_particle, y_particle, s_here, r
integer at, at2
logical do_tilt, err
logical, optional :: check_momentum
character(20) :: r_name = 'check_aperture_limit'

! Check if there is a limit here. If not, simply return.

at2 = physical_ele_end (at, ele%orientation)

if (at2 == entrance_end$) then
  if (ele%aperture_at /= entrance_end$ .and. ele%aperture_at /= both_ends$ .and. ele%aperture_at /= continuous$) return
elseif (at2 == exit_end$) then
  if (ele%aperture_at /= exit_end$ .and. ele%aperture_at /= both_ends$ .and. ele%aperture_at /= continuous$) return
endif

! Custom

if (ele%aperture_type == custom$) then
  call check_aperture_limit_custom (orb, ele, at, param, err)
  return
endif

! Check p_x and p_y

if (logic_option(.true., check_momentum)) then
  if (abs(orb%vec(2)) > 1+orb%vec(6) .and. abs(orb%vec(2)) > abs(orb%vec(4))) then
    if (orb%vec(2) > 0) then; orb%state = lost_pos_x_aperture$
    else;                     orb%state = lost_neg_x_aperture$
    endif
    param%unstable_factor = 100 * abs(orb%vec(2)) 
    return
  elseif (abs(orb%vec(4)) > 1 + orb%vec(6)) then
    if (orb%vec(4) > 0) then; orb%state = lost_pos_y_aperture$
    else;                     orb%state = lost_neg_y_aperture$
    endif
    param%unstable_factor = 100 * abs(orb%vec(4)) 
    return
  endif
endif

!

if (ele%offset_moves_aperture .and. (at2 == entrance_end$ .or. at2 == exit_end$)) then
  do_tilt = .false.
  if (ele%key == ecollimator$ .or. ele%key == rcollimator$) do_tilt = .true.
  orb2 = orb
  s_here = orb2%s - (ele%s - ele%value(l$))
  call offset_particle (ele, orb2, param, set$, set_canonical = .false., &
               set_tilt = do_tilt, set_multipoles = .false., set_hvkicks = .false., &
               ds_pos = s_here)
  x_particle = orb2%vec(1)
  y_particle = orb2%vec(3)
elseif (ele%aperture_at == surface$) then
  x_particle = orb%vec(5)
  y_particle = orb%vec(3)
else
  x_particle = orb%vec(1)
  y_particle = orb%vec(3)
endif

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

case (rectangular$)

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
! Subroutine track_a_drift (orb, ele, length)
!
! Subroutine to track a particle as through a drift.
!
! Modules needed:
!   use precision_def
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
! Subroutine track_a_drift_photon (orb, ele, length)
!
! Subroutine to track a particle as through a drift.
!
! Modules needed:
!   use precision_def
!
! Input:
!   orb      -- coord_struct: Orbit at start of the drift.
!   ele      -- Ele_struct: Element tracked through.
!   length   -- Real(rp): Length to drift through.
!
! Output:
!   orb      -- coord_struct: Orbit at end of the drift
!-

subroutine track_a_drift_photon (orb, ele, length)

implicit none

type (coord_struct) orb
type (ele_struct) ele
type (lat_param_struct) param
real(rp) length, rel_pc, dz, px, py, pz, pxy2

! Photon tracking uses a different coordinate system. 
! Notice that if orb%vec(6) is negative then the photon will be going back in time.

orb%vec(1) = orb%vec(1) + length * orb%vec(2) / orb%vec(6)
orb%vec(3) = orb%vec(3) + length * orb%vec(4) / orb%vec(6)
orb%vec(5) = orb%vec(5) + length
orb%s      = orb%s      + length
orb%t = orb%t + length * orb%vec(6) / c_light

end subroutine track_a_drift_photon

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
!   ele    -- Ele_struct: Bend element.
!   param  -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb     -- Coord_struct: End position.
!-

subroutine track_a_bend (start_orb, ele, param, end_orb)

use multipole_mod

implicit none

type (coord_struct), intent(in)  :: start_orb
type (coord_struct), intent(out) :: end_orb
type (ele_struct),   intent(inout)  :: ele
type (lat_param_struct), intent(inout) :: param

real(rp) angle, ct, st, x, px, y, py, z, pz, dpx_t, p_long
real(rp) rel_p, rel_p2, Dy, px_t, factor, rho, g, g_err, c_dir
real(rp) length, g_tot, del_p, eps, pxy2, f, k_2, alpha, beta
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, z_2, s_y, c_y, r(6)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer n, n_step, ix_fringe

logical has_nonzero_pole

!-----------------------------------------------------------------------

end_orb = start_orb
call offset_particle (ele, end_orb, param, set$, set_canonical = .false., set_multipoles = .false.)

! Entrance edge kick

call bend_edge_kick (end_orb, ele, param, upstream_end$, .false.)

! If we have a sextupole component then step through in steps of length ds_step

n_step = 1

call multipole_ele_to_kt(ele, param, .true., has_nonzero_pole, knl, tilt)
if (ele%value(k2$) /= 0 .or. has_nonzero_pole) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)
knl = knl / n_step

! Set some parameters

c_dir = ele%orientation * param%rel_tracking_charge
length = ele%value(l$) / n_step
g = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * c_dir
g_err = g_tot - g
angle = g * length
rho = 1 / g
del_p = start_orb%vec(6)
rel_p  = 1 + del_p
rel_p2 = rel_p**2
k_1 = ele%value(k1$) * c_dir
k_2 = ele%value(k2$) * c_dir

! 1/2 sextupole kick at the beginning.

if (k_2 /= 0) then
  end_orb%vec(2) = end_orb%vec(2) + k_2/2 * length * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
  end_orb%vec(4) = end_orb%vec(4) + k_2/2 * length * end_orb%vec(1) * end_orb%vec(3)
end if
if (has_nonzero_pole) call multipole_kicks (knl/2, tilt, end_orb)

! And track with n_step steps

do n = 1, n_step

  ! with k1 /= 0 use small angle approximation

  if (g == 0 .or. k_1 /= 0) then

    call sbend_body_with_k1_map (ele, param, n_step, end_orb, end_orb = end_orb)

  !-----------------------------------------------------------------------
  ! Track through main body...
  ! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  else

    ct = cos(angle)
    st = sin(angle)

    x  = end_orb%vec(1)
    px = end_orb%vec(2)
    y  = end_orb%vec(3)
    py = end_orb%vec(4)
    z  = end_orb%vec(5)
    pz = end_orb%vec(6)
   
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
      f = pxy2 / (2 * rel_p)
      f = del_p - f - f*f/2 - g_err*rho - g_tot*x
    else
      f = p_long - g_tot * (1 + x * g) / g
    endif

    Dy  = sqrt(rel_p2 - py**2)
    px_t = px*ct + f*st
    dpx_t = -px*st*g + f*ct*g

    if (abs(px) > Dy .or. abs(px_t) > Dy) then
      if (max(px, px_t) > 0) then; end_orb%state = lost_pos_x_aperture$
      else;                        end_orb%state = lost_neg_x_aperture$
      endif
      return
    endif    

    if (abs(g_tot) < 1e-5 * abs(g)) then
      alpha = p_long * ct - px * st
      end_orb%vec(1) = (p_long * (1 + g * x) - alpha) / (g * alpha) - &
                   g_tot * (Dy * (1 + g * x) * st)**2 / (2 * alpha**3 * g**2) + &
                   g_tot**2 * Dy**2 * ((1 + g * x) * st)**3 * (px * ct + p_long * st) / (2 * alpha**5 * g**3)
    else
      eps = px_t**2 + py**2
      if (eps < 1e-5 * rel_p2 ) then  ! use small angle approximation
        eps = eps / (2 * rel_p)
        end_orb%vec(1) = (del_p - g_err / g - rho*dpx_t + eps * (eps / (2 * rel_p) - 1)) / g_tot
      else
        end_orb%vec(1) = (sqrt(rel_p2 - eps) - rho*dpx_t - rho*g_tot) / g_tot
      endif
    endif

    end_orb%vec(2) = px_t
    end_orb%vec(4) = py
    end_orb%vec(6) = pz

    if (abs(g_tot) < 1e-5 * abs(g)) then
      beta = (1 + g * x) * st / (g * alpha) - &
             g_tot * (px * ct + p_long * st) * (st * (1 + g * x))**2 / (2 * g**2 * alpha**3)
      end_orb%vec(3) = y + py * beta
      end_orb%vec(5) = z + length  - rel_p * beta 
    else
      factor = (asin(px/Dy) - asin(px_t/Dy)) / g_tot
      end_orb%vec(3) = y + py * (angle/g_tot + factor)
      end_orb%vec(5) = z + length * (g_err - g*del_p) / g_tot - rel_p * factor
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

call bend_edge_kick (end_orb, ele, param, downstream_end$, .false.)

call offset_particle (ele, end_orb, param, unset$, set_canonical = .false., set_multipoles = .false.)

call track1_low_energy_z_correction (end_orb, ele, param)

end subroutine track_a_bend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine approx_bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)
!
! Subroutine to track through the edge field of an sbend.
! In_to_out tracking starts with the particle just outside the bend and
! returns the orbit that the particle had just inside the bend.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   stream_end  -- Integer: upstream_end$ or downstream_end$
!   in_to_out   -- Logical: If True then make the inverse transformation.
!                    That is, for the entrance end take the input orb as the coordinates
!                    just inside the entrance end of the bend and return the coordinates 
!                    just oustide the entrance end.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   kx, ky     -- Real(rp), optional: Horizontal and vertical edge focusing strengths.
!                  Useful for constructing the edge transfer matrix.
!                  The values of kx and ky are not affected by the in_to_out argument.
!-

subroutine approx_bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: kx, ky
real(rp) e, g, g_tot, fint, hgap, ht_x, ht_y, cos_e, sin_e, tan_e, sec_e, v0(6), k1_eff
real(rp) ht2, hs2, c_dir
integer stream_end, element_end
logical in_to_out

character(24), parameter :: r_name = 'approx_bend_edge_kick'

! Track through the entrence face. 
! See MAD physics guide for writeup. Note that MAD does not have a g_err.

c_dir = param%rel_tracking_charge * ele%orientation
element_end = physical_ele_end(stream_end, ele%orientation)

g     = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * c_dir

if (element_end == entrance_end$) then
  e = ele%value(e1$); fint = ele%value(fint$); hgap = ele%value(hgap$)
else
  e = ele%value(e2$); fint = ele%value(fintx$); hgap = ele%value(hgapx$)
endif

cos_e = cos(e); sin_e = sin(e); tan_e = sin_e / cos_e; sec_e = 1 / cos_e
ht_x = g_tot * tan_e
ht2 = g * tan_e**2
hs2 = g * sec_e**2
k1_eff = ele%value(k1$) * c_dir

if (fint == 0) then
  ht_y = -ht_x
else
  ht_y = -g_tot * tan(e - 2 * fint * g_tot * hgap * (1 + sin_e**2) / cos_e)
endif

v0 = orb%vec

if (in_to_out) then
  if (stream_end == upstream_end$) then
    orb%vec(1) = v0(1) + ht2 * v0(1)**2 / 2 - hs2 * v0(3)**2 / 2
    orb%vec(2) = v0(2) - ht_x * v0(1) + ht2 * (v0(3) * v0(4) - v0(1) * v0(2)) - &
                         k1_eff * tan_e * (v0(1)**2 - v0(3)**2) + &
                         ht_x * ht2 * (v0(1)**2 + v0(3)**2) / 2
    orb%vec(3) = v0(3) - ht2 * v0(1) * v0(3)
    orb%vec(4) = v0(4) - ht_y * v0(3) + ht2 * v0(1) * v0(4) + hs2 * v0(2) * v0(3) - &
                         (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(1) * v0(3)
  else
    orb%vec(1) = v0(1) - ht2 * v0(1)**2 / 2 + hs2 * v0(3)**2 / 2
    orb%vec(2) = v0(2) - ht_x * v0(1) + ht2 * (v0(1) * v0(2) - v0(3) * v0(4)) - &
                         k1_eff * tan_e * (v0(1)**2 - v0(3)**2) - &
                         ht_x * (ht2 + hs2) * v0(3)**2 / 2
    orb%vec(3) = v0(3) + ht2 * v0(1) * v0(3) 
    orb%vec(4) = v0(4) - ht_y * v0(3) - ht2 * v0(1) * v0(4) - hs2 * v0(2) * v0(3) + &
                         2 * k1_eff * tan_e * v0(1) * v0(3) 
  endif

else
  if (stream_end == upstream_end$) then
    orb%vec(1) = v0(1) - ht2 * v0(1)**2 / 2 + hs2 * v0(3)**2 / 2
    orb%vec(2) = v0(2) + ht_x * v0(1) + ht2 * (v0(1) * v0(2) - v0(3) * v0(4)) + &
                         k1_eff * tan_e * (v0(1)**2 - v0(3)**2) + &
                         ht_x * (ht2 + hs2) * v0(3)**2 / 2
    orb%vec(3) = v0(3) + ht2 * v0(1) * v0(3) 
    orb%vec(4) = v0(4) + ht_y * v0(3) - ht2 * v0(1) * v0(4) - hs2 * v0(2) * v0(3) - &
                         2 * k1_eff * tan_e * v0(1) * v0(3) 
  else
    orb%vec(1) = v0(1) + ht2 * v0(1)**2 / 2 - hs2 * v0(3)**2 / 2
    orb%vec(2) = v0(2) + ht_x * v0(1) + ht2 * (v0(3) * v0(4) - v0(1) * v0(2)) + &
                         k1_eff * tan_e * (v0(1)**2 - v0(3)**2) - &
                         ht_x * ht2 * (v0(1)**2 + v0(3)**2) / 2
    orb%vec(3) = v0(3) - ht2 * v0(1) * v0(3)
    orb%vec(4) = v0(4) + ht_y * v0(3) + ht2 * v0(1) * v0(4) + hs2 * v0(2) * v0(3) + &
                         (ht_x * hs2 - 2 * tan_e * k1_eff) * v0(1) * v0(3)
  endif
endif

if (present(kx)) kx = ht_x 
if (present(ky)) ky = ht_y

end subroutine approx_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine apply_hard_edge_kick (orb, s_edge, t_rel, hard_ele, track_ele, param, element_end)
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
!   element_end -- Integer: upstream_end$ or downstream_end$.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!-

subroutine apply_hard_edge_kick (orb, s_edge, t_rel, hard_ele, track_ele, param, element_end)

implicit none

type (ele_struct) hard_ele, track_ele
type (coord_struct) orb
type (lat_param_struct) param
type (em_field_struct) field

real(rp) t, f, l_drift, ks, t_rel, s_edge, s

integer element_end, ix_fringe

! 

if (hard_ele%field_calc /= bmad_standard$) return

select case (hard_ele%key)
case (sbend$)
  ix_fringe = nint(hard_ele%value(fringe_type$))
  if (ix_fringe == full_straight$ .or. ix_fringe == full_bend$) then
    call exact_bend_edge_kick (orb, hard_ele, param, element_end, .false.)
  elseif (ix_fringe == basic_bend$) then
    call approx_bend_edge_kick (orb, hard_ele, param, element_end, .false.)
  endif

! Note: Cannot trust hard_ele%value(ks$) here since element may be superimposed with an lcavity.
! So use hard_ele%value(bs_field$).

case (solenoid$, sol_quad$, bend_sol_quad$)
  ks = param%rel_tracking_charge * hard_ele%value(bs_field$) * c_light / orb%p0c
  if (element_end == upstream_end$) then
    orb%vec(2) = orb%vec(2) + ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) - ks * orb%vec(1) / 2
  else
    orb%vec(2) = orb%vec(2) - ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) + ks * orb%vec(1) / 2
  endif

case (lcavity$, rfcavity$, e_gun$)

  ! Add on bmad_com%significant_length to make sure we are just inside the cavity.
  f = charge_of(param%particle) / (2 * orb%p0c)
  t = t_rel + track_ele%value(ref_time_start$) - hard_ele%value(ref_time_start$) 
  s = s_edge

  if (element_end == upstream_end$) then

    if (hard_ele%key == e_gun$) return  ! E_gun does not have an entrance kick
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
end select

end subroutine apply_hard_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)
!
! Subroutine to track through the edge field of an sbend.
! In_to_out tracking starts with the particle just outside the bend and
! returns the orbit that the particle had just inside the bend.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: Rel charge.
!   stream_end  -- Integer: upstream_end$ or downstream_end$
!   in_to_out   -- Logical: If True then make the inverse transformation.
!                    That is, for the entrance end take the input orb as the coordinates
!                    just inside the entrance end of the bend and return the coordinates 
!                    just oustide the entrance end.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   kx, ky     -- Real(rp), optional: Horizontal and vertical edge focusing strengths.
!                  Useful for constructing the edge transfer matrix.
!                  The values of kx and ky are not affected by the in_to_out argument.
!-

subroutine bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp), optional :: kx, ky
integer stream_end, element_end, ix_fringe
logical in_to_out

ix_fringe = nint(ele%value(fringe_type$))
if (present(kx) .and. present(ky)) then 
   if (ix_fringe == full_straight$ .or. ix_fringe == full_bend$) then
      call exact_bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)
   elseif (ix_fringe == basic_bend$) then
      call approx_bend_edge_kick (orb, ele, param, stream_end, in_to_out, kx, ky)
   endif
else
   if (ix_fringe == full_straight$ .or. ix_fringe == full_bend$) then
      call exact_bend_edge_kick (orb, ele, param, stream_end, in_to_out)
   elseif (ix_fringe == basic_bend$) then
      call approx_bend_edge_kick (orb, ele, param, stream_end, in_to_out)
   endif
end if

end subroutine bend_edge_kick

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine pitches_to_rotation_matrix (x_pitch, y_pitch, set, rot_mat)
!
! Routine to create a rotation matrix form the pitches and tilt
!
! Input:
!   x_pitch -- Real(rp): X-pitch
!   y_pitch -- Real(rp): Y-pitch
!   set     -- Logical: set$ (True)    -> rot translates from lab to element coords.
!                       unset$ (False) -> rot translates from element to lab coords.
! Output:
!   rot_mat(3,3) -- Real(rp): Rotation matrix.
!-

subroutine pitches_to_rotation_matrix (x_pitch, y_pitch, set, rot_mat)

implicit none

real(rp) x_pitch, y_pitch, rot_mat(3,3)
real(rp) sx, sy, nx, ny, cos_t, norm

logical set

! Degenerate case

if (x_pitch == 0 .and. y_pitch == 0) then
  call mat_make_unit(rot_mat)
  return
endif

!

sx = sin(x_pitch)
sy = sin(y_pitch)

if (set) then
  sx = -sx
  sy = -sy
endif

norm = sqrt(sx**2 + sy**2)
nx = -sy / norm
ny =  sx / norm
cos_t = sqrt(1 - norm**2)

rot_mat(1,:) = [nx**2 + ny**2 * cos_t, nx * ny * (1 - cos_t), sx]
rot_mat(2,:) = [nx * ny * (1 - cos_t), ny**2 + nx**2 * cos_t, sy]
rot_mat(3,:) = [-sx,                   -sy,                   cos_t]

end subroutine pitches_to_rotation_matrix

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_wedger (a, g_tot, beta0, X)
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
!-
subroutine ptc_wedger (a, g_tot, beta0, X)

implicit none

real(rp) :: X(6)
real(rp) :: a, beta0, g_tot
real(rp) :: Xn(6),pz,pzs,pt,b1
character(20) :: r_name = 'ptc_wedger'

! No net field case...

b1 = g_tot

if(b1==0) then
   call ptc_rot_xz(a, X, beta0)
   return
endif

! Normal case

pz=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(2)**2-X(4)**2)

Xn(2)=X(2)*cos(a)+(pz-b1*X(1))*sin(a)

pt=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-X(4)**2)
pzs=sqrt(1.0_rp+2.0_rp*X(5)/beta0+X(5)**2-Xn(2)**2-X(4)**2)

Xn(1)=X(1)*cos(a)+(X(1)*X(2)*sin(2.0_rp*a)+sin(a)**2*(2.0_rp*X(1)*pz-b1*X(1)**2) )&
    / (pzs+pz*cos(a)-X(2)*sin(a))
Xn(3)=(a+asin(X(2)/pt)-asin(Xn(2)/pt))/b1

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
! private subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, element_end)
!
! Subroutine to compute the exact hard edge fringe field of a bend.
! 
! Adapted from forest/code/Sh_def_kind.f90 : FRINGE_DIPOLER
!
!
! Input:
!   X(6)   -- real(rp): PTC phase space coordinates
!   beta0  -- real(rp): reference relativistic beta
!   g_tot  -- real(rp): reference bending radius
!   fint   -- real(rp): field integral for pole face
!   hgap   -- real(rp): gap height at pole face in meters
!                       Only the product fint*hgap is used
!   element_end -- Integer: upstream_end$ or downstream_end$
!
! Output:
!   X(6)   -- real(rp): PTC phase space coordinates
!-

subroutine ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, element_end)

implicit none

real(rp) :: X(6) !PTC phase space coordinates
real(rp) :: FINT, HGAP
real(rp) :: beta0, g_tot
real(rp) :: PZ,XP,YP,TIME_FAC
real(rp) :: D(3,3),FI(3),FI0,B,co1,co2
integer  :: i
integer  :: element_end
character(20) :: r_name = 'ptc_fringe_dipoler'

!

if (element_end == downstream_end$) then
   B = -g_tot  !EL%CHARGE*BN(1)
else if (element_end == upstream_end$) then
   B = g_tot       
else
  call out_io (s_fatal$, r_name, 'INVALID ELEMENT_END')
  call err_exit
endif


pz=sqrt(1.0_rp+2.0_rp*x(5)/beta0+x(5)**2-x(2)**2-x(4)**2)
time_fac=1.0_rp/beta0+x(5)

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

fi0=b*tan(fi0)

b=0
do i=1,3
   b=fi(i)*d(i,2)+b
enddo
x(3)=2.0_rp*x(3)/(1.0_rp+ sqrt(1.0_rp-2.0_rp*b*x(3)) )
x(4)=x(4)-fi0*x(3)

b=0
do i=1,3
   b=fi(i)*d(i,1)+b
enddo
x(1)=x(1)+0.5_rp*b*x(3)**2

b=0
do i=1,3
   b=fi(i)*d(i,3)+b
enddo
x(6)=x(6)-0.5_rp*b*x(3)**2
    
end subroutine ptc_fringe_dipoler

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! private subroutine ptc_rot_xz(a, X, beta0)
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
!-

subroutine ptc_rot_xz(a, X, beta0)

implicit none

real(rp) :: x(6)
real(rp) :: xn(6),pz,pt
real(rp) :: a, beta0
character(20) :: r_name = 'ptc_rot_xz'

!

pz=sqrt(1.0_rp+2.0_rp*x(5)/ beta0+x(5)**2-x(2)**2-x(4)**2)
pt=1.0_rp-x(2)*tan(a)/pz
xn(1)=x(1)/cos(a)/pt
xn(2)=x(2)*cos(a)+sin(a)*pz
xn(3)=x(3)+x(4)*x(1)*tan(a)/pz/pt
xn(6)=x(6)+x(1)*tan(a)/pz/pt*(1.0_rp/ beta0+x(5))

x(1)=xn(1)
x(2)=xn(2)
x(3)=xn(3)
x(6)=xn(6)

end subroutine ptc_rot_xz

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine exact_bend_edge_kick (orb, ele, param, element_end, in_to_out, kx, ky)
!
! Subroutine to track through the edge field of an sbend.
! In_to_out tracking starts with the particle just outside the bend and
! returns the orbit that the particle had just inside the bend.
!
! Uses routines adapted from PTC
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   param       -- lat_param_struct: 
!   element_end -- Integer: upstream_end$ or downstream_end$
!   in_to_out   -- Logical: If True then make the inverse transformation.
!                     That is, for the entrance end take the input orb as the coordinates
!                     just inside the entrance end of the bend and return the coordinates 
!                     just oustide the entrance end.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   kx, ky     -- Real(rp), optional: Horizontal and vertical edge focusing strengths.
!                  Useful for constructing the edge transfer matrix.
!                  The values of kx and ky are not affected by the in_to_out argument.
!-

subroutine exact_bend_edge_kick (orb, ele, param, element_end, in_to_out, kx, ky)

use ptc_interface_mod

implicit none

type(coord_struct) :: orb
type(ele_struct) :: ele
type (lat_param_struct) param
real(rp), optional :: kx, ky
real(rp) :: X(6), ct
real(rp) :: beta0, g_tot, edge_angle, hgap, fint
integer :: element_end
logical :: in_to_out

character(20) :: r_name = 'exact_bend_edge_kick'

!
if (in_to_out) then
  call out_io (s_fatal$, r_name, 'IN_TO_OUT NOT IMPLEMENTED')
  call err_exit
endif

!Get reference beta0

beta0 = ele%value(e_tot$) / ele%value(p0c$)
g_tot = ele%value(g$) + ele%value(g_err$)


! Convert to PTC coordinates
call vec_bmad_to_ptc(orb%vec, beta0, X)
!Save time
ct = X(6)


if (element_end == upstream_end$) then
  edge_angle = ele%value(e1$)
  fint = ele%value(FINT$)
  hgap = ele%value(HGAP$)
  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X)
  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, element_end)
  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X)

else if (element_end == downstream_end$) then
  edge_angle = ele%value(e2$)
  fint = ele%value(FINTX$)
  hgap = ele%value(HGAPX$)
  ! Backtrack
  call ptc_wedger(-edge_angle, g_tot, beta0, X)
  ! Edge kick
  call ptc_fringe_dipoler(X, g_tot, beta0, fint, hgap, element_end)
  ! Drift forward
  call ptc_wedger(edge_angle, 0.0_rp, beta0, X)

else
  !error!
  call err_exit
endif

! Convert back to bmad coordinates
call vec_ptc_to_bmad (X, beta0, orb%vec)

! Correct time
orb%t = orb%t + (X(6) - ct)/c_light

! Focusing terms 
if (present(kx)) kx = g_tot * tan(edge_angle)
if (present(ky)) then
  if (fint == 0) then
    ky = -g_tot * tan(edge_angle)
  else
    ky = -g_tot * tan(edge_angle - 2 * fint * g_tot * hgap * (1 + sin(edge_angle)**2) / cos(edge_angle))
  endif
endif

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

mass = mass_of(param%particle)
e_tot = ele%value(e_tot$)
p0c = ele%value(p0c$)

if (abs(orbit%vec(6)) < 1e-6 * mass**2 * p0c / e_tot**3) then
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * orbit%vec(6) * (1 - 3 * orbit%vec(6) / 2) * (mass / e_tot)**2
else
  pc = (1 + orbit%vec(6)) * ele%value(p0c$)
  call convert_pc_to (pc, param%particle, beta = beta)
  beta0 = ele%value(p0c$) / ele%value(e_tot$)
  orbit%vec(5) = orbit%vec(5) + ele%value(l$) * (beta - beta0) / beta0
endif

end subroutine track1_low_energy_z_correction

end module
