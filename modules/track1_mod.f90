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
!   at    -- Integer: entrance_end$ or exit_end$
!   param -- lat_param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. In this case param%lost will
!               be set to True.
!   check_momentum
!         -- Logical, optional -- If present and false then checking of p_x and 
!               p_y will be disabled.
!
! Output:
!   param -- lat_param_struct: Parameter structure:
!     %lost -- Set True if the orbit is outside the aperture. 
!              Note: %lost is NOT set False if the orbit is inside 
!                the aperture.
!     %plane_lost_at   -- Integer: Plane where particle is lost:
!                           x_plane$ or y_plane$
!     %unstable_factor -- Real(rp): |orbit_amp/limit|
!-

subroutine check_aperture_limit (orb, ele, at, param, check_momentum)

implicit none

type (coord_struct) :: orb
type (coord_struct) orb2 
type (ele_struct) :: ele
type (lat_param_struct), intent(inout) :: param

real(rp) x_lim, y_lim, x_beam, y_beam, s_here, r
integer at
logical do_tilt, err
logical, optional :: check_momentum
character(20) :: r_name = 'check_aperture_limit'

! Custom

if (ele%aperture_type == custom$) then
  call check_aperture_limit_custom (orb, ele, at, param, err)
  return
endif

! Check p_x and p_y

if (logic_option(.true., check_momentum)) then
  if (abs(orb%vec(2)) > 1 .or. abs(orb%vec(4)) > 1) then
    param%lost = .true.
    orb%status = dead$
    if (abs(orb%vec(2)) > abs(orb%vec(4))) then
      param%plane_lost_at = x_plane$
      param%unstable_factor = 100 * abs(orb%vec(2)) 
    else
      param%plane_lost_at = y_plane$
      param%unstable_factor = 100 * abs(orb%vec(4)) 
    endif
    return
  endif
endif

!

if (ele%offset_moves_aperture) then
  do_tilt = .false.
  if (ele%key == ecollimator$ .or. ele%key == rcollimator$) do_tilt = .true.
  orb2 = orb
  s_here = 0
  if (at == exit_end$) s_here = ele%value(l$)
  call offset_particle (ele, param, orb2, set$, set_canonical = .false., &
               set_tilt = do_tilt, set_multipoles = .false., set_hvkicks = .false., &
               ds_pos = s_here)
  x_beam = orb2%vec(1)
  y_beam = orb2%vec(3)
else
  x_beam = orb%vec(1)
  y_beam = orb%vec(3)
endif

!

if (x_beam < 0) then
  x_lim = ele%value(x1_limit$)
else
  x_lim = ele%value(x2_limit$)
endif

if (x_lim <= 0 .or. .not. param%aperture_limit_on) x_lim = bmad_com%max_aperture_limit

if (y_beam < 0) then
  y_lim = ele%value(y1_limit$)
else
  y_lim = ele%value(y2_limit$)
endif

if (y_lim <= 0 .or. .not. param%aperture_limit_on) &
                                  y_lim = bmad_com%max_aperture_limit

if (x_lim == 0 .and. y_lim == 0) return

select case (ele%aperture_type)

case (elliptical$)
  if (x_lim == 0 .or. y_lim == 0) then
    call out_io (s_fatal$, r_name, &
              'ECOLLIMATOR HAS ONE LIMIT ZERO AND THE OTHER NOT: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
  endif

  r = (x_beam / x_lim)**2 + (y_beam / y_lim)**2
  if (r > 1) then
    param%lost = .true.
    orb%status = dead$
    if (abs(x_beam / x_lim) > abs(y_beam / y_lim)) then
      param%plane_lost_at = x_plane$
    else
      param%plane_lost_at = y_plane$
    endif
    param%unstable_factor = sqrt(r) - 1
  endif

case (rectangular$)

  if (abs(x_beam) > x_lim .or. abs(y_beam) > y_lim) then
    param%lost = .true.
    orb%status = dead$
    if (abs(x_beam)/x_lim > abs(y_beam)/y_lim) then
      param%plane_lost_at = x_plane$
      param%unstable_factor = abs(x_beam) / x_lim - 1
    else
      param%plane_lost_at = y_plane$
      param%unstable_factor = abs(y_beam) / y_lim - 1
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
! Subroutine track_a_drift (orb, length)
!
! Subroutine to track a particle as through a drift.
!
! Modules needed:
!   use precision_def
!
! Input:
!   orb    -- coord_struct: Orbit at start of the drift
!   length -- Real(rp): Length of drift.
!
! Output:
!   orb%vec(:) -- coord_struct: Orbit at end of the drift
!-

subroutine track_a_drift (orb, length)

implicit none

type (coord_struct) orb
type (lat_param_struct) param
real(rp) length, rel_pc, dz

!
rel_pc = 1 + orb%vec(6)

orb%vec(1) = orb%vec(1) + length * orb%vec(2) / rel_pc
orb%vec(3) = orb%vec(3) + length * orb%vec(4) / rel_pc
dz = -length * (orb%vec(2)**2 + orb%vec(4)**2) / (2 * rel_pc**2)
orb%vec(5) = orb%vec(5) + dz
if (orb%beta /= 0) orb%t = orb%t + (length - dz) / (orb%beta * c_light)
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
real(rp) rel_p, rel_p2, Dy, px_t, factor, rho, g, g_err
real(rp) length, g_tot, del_p, eps, pxy2, f, k_2, alpha, beta
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, z_2, s_y, c_y, r(6)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer n, n_step

!-----------------------------------------------------------------------
! simple case

if (ele%value(g$) == 0) then
  length = ele%value(l$)
  end_orb = start_orb
  end_orb%vec(2) = end_orb%vec(2) - length * ele%value(g_err$) / 2
  call track_a_drift (end_orb, length)
  end_orb%vec(2) = end_orb%vec(2) - length * ele%value(g_err$) / 2
  return
endif

!-----------------------------------------------------------------------

end_orb = start_orb
call offset_particle (ele, param, end_orb, set$, set_canonical = .false., set_multipoles = .false.)
call apply_bend_edge_kick (end_orb, ele, entrance_end$, .false.)

! If we have a sextupole component then step through in steps of length ds_step

n_step = 1
if (ele%value(k2$) /= 0 .or. associated(ele%a_pole)) &
                  n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)

if (associated(ele%a_pole)) then
  call multipole_ele_to_kt(ele, param%particle, knl, tilt, .true.)
  knl = knl / n_step
endif

! Set some parameters

length = ele%value(l$) / n_step
g = ele%value(g$)
g_err = ele%value(g_err$)
g_tot = g + g_err
angle = ele%value(g$) * length
rho = 1 / g
del_p = start_orb%vec(6)
rel_p  = 1 + del_p
rel_p2 = rel_p**2
k_1 = ele%value(k1$)
k_2 = ele%value(k2$) * length

! 1/2 sextupole kick at the beginning.

if (k_2 /= 0) call multipole_kick (k_2/2, 0.0_rp, 2, end_orb)
if (associated(ele%a_pole)) call multipole_kicks (knl/2, tilt, end_orb)

! And track with n_step steps

do n = 1, n_step

  ! with k1 /= 0 use small angle approximation

  if (k_1 /= 0) then

    call sbend_body_with_k1_map (g, g_err, length, k_1, end_orb%vec, end = end_orb%vec)

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
      param%lost = .true.
      end_orb%status = dead$
      param%plane_lost_at = x_plane$
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
      param%lost = .true.
      param%plane_lost_at = x_plane$
      end_orb%status = dead$
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
    if (k_2 /= 0) call multipole_kick (k_2/2, 0.0_rp, 2, end_orb)
    if (associated(ele%a_pole)) call multipole_kicks (knl/2, tilt, end_orb)
  else
    if (k_2 /= 0) call multipole_kick (k_2, 0.0_rp, 2, end_orb)
    if (associated(ele%a_pole)) call multipole_kicks (knl, tilt, end_orb)
  endif

enddo

! Track through the exit face. Treat as thin lens.

call apply_bend_edge_kick (end_orb, ele, exit_end$, .false.)
call offset_particle (ele, param, end_orb, unset$, set_canonical = .false., set_multipoles = .false.)

end subroutine track_a_bend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine apply_bend_edge_kick (orb, ele, element_end, reverse, kx, ky)
!
! Subroutine to track through the edge field of an sbend.
! Reverse tracking starts with the particle just outside the bend and
! returns the orbit that the particle had just inside the bend.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   ele         -- ele_struct: SBend element.
!   element_end -- Integer: entrance_end$ or exit_end$
!   reverse     -- Logical: If True then take the input orb as the position
!                    just outside the bend and output the position just inside the bend. 
!                    This does not affect the values of kx and ky
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!   kx, ky     -- Real(rp), optional: Horizontal and vertical edge focusing strengths.
!                  Useful for constructing the edge transfer matrix.
!                  The values of kx and ky are not affected by the reverse argument.
!-

subroutine apply_bend_edge_kick (orb, ele, element_end, reverse, kx, ky)

type (ele_struct) ele
type (coord_struct) orb
real(rp), optional :: kx, ky
real(rp) e, g_tot, fint, hgap, k1x, k1y
integer element_end
logical reverse

! Track through the entrence face. Treat as thin lens.

g_tot = ele%value(g$) + ele%value(g_err$)

if (element_end == entrance_end$) then
  e = ele%value(e1$); fint = ele%value(fint$); hgap = ele%value(hgap$)
else
  e = ele%value(e2$); fint = ele%value(fintx$); hgap = ele%value(hgapx$)
endif

k1x = g_tot * tan(e)
if (fint == 0) then
  k1y = -k1x
else
  k1y = -g_tot * tan(e - 2 * fint * g_tot * hgap * (1 + sin(e)**2) / cos(e))
endif

if (reverse) then
  orb%vec(2) = orb%vec(2) - k1x * orb%vec(1)
  orb%vec(4) = orb%vec(4) - k1y * orb%vec(3)
else
  orb%vec(2) = orb%vec(2) + k1x * orb%vec(1)
  orb%vec(4) = orb%vec(4) + k1y * orb%vec(3)
endif

if (present(kx)) kx = k1x 
if (present(ky)) ky = k1y

end subroutine apply_bend_edge_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine apply_hard_edge_kick (orb, t_rel, hard_ele, track_ele, param, element_end)
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
!   lcavity and rfcavity 
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords in element reference frame.
!   t_rel       -- real(rp): Time relative to track_ele entrance edge
!   hard_ele    -- ele_struct: Element with hard edges.
!   track_ele   -- ele_struct: Element being tracked through. 
!                    Is different from hard_ele when there are superpositions.
!   param       -- lat_param_struct: lattice parameters.
!   element_end -- Integer: entrance_end$ or exit_end$.
!
! Output:
!   orb        -- Coord_struct: Coords after tracking.
!-

subroutine apply_hard_edge_kick (orb, t_rel, hard_ele, track_ele, param, element_end)

implicit none

type (ele_struct) hard_ele, track_ele
type (coord_struct) orb
type (lat_param_struct) param
type (em_field_struct) field

real(rp) t, f, l_drift, dref_time, ks, t_rel

integer element_end

! 

if (hard_ele%field_calc /= bmad_standard$) return

select case (hard_ele%key)
case (sbend$)
  call apply_bend_edge_kick (orb, hard_ele, element_end, .false.)

! Note: Cannot trust hard_ele%value(ks$) here since element may be superimposed with an lcavity.
! So use hard_ele%value(bs_field$).

case (solenoid$, sol_quad$)
  ks = hard_ele%value(bs_field$) * c_light / orb%p0c
  if (element_end == entrance_end$) then
    orb%vec(2) = orb%vec(2) + ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) - ks * orb%vec(1) / 2
  else
    orb%vec(2) = orb%vec(2) - ks * orb%vec(3) / 2
    orb%vec(4) = orb%vec(4) + ks * orb%vec(1) / 2
  endif

case (lcavity$, rfcavity$)

  ! Add on bmad_com%significant_length to make sure we are just inside the cavity.
  l_drift = (hard_ele%value(l$) - hard_ele%value(l_hard_edge$)) / 2 + bmad_com%significant_length
  f = charge_of(param%particle) / (2 * orb%p0c)
  t = t_rel + hard_ele%value(ref_time_start$) - track_ele%value(ref_time_start$)

  if (element_end == entrance_end$) then
    call em_field_calc (hard_ele, param, l_drift, t, orb, .true., field)

    orb%vec(2) = orb%vec(2) - field%e(3) * orb%vec(1) * f + c_light * field%b(3) * orb%vec(3) * f
    orb%vec(4) = orb%vec(4) - field%e(3) * orb%vec(3) * f - c_light * field%b(3) * orb%vec(1) * f

  else
    ! dref_time is the reference time to cross just the cavity. That is, this subtracts off the drift times.
    call em_field_calc (hard_ele, param, hard_ele%value(l$)-l_drift, t, orb, .true., field)

    orb%vec(2) = orb%vec(2) + field%e(3) * orb%vec(1) * f - c_light * field%b(3) * orb%vec(3) * f
    orb%vec(4) = orb%vec(4) + field%e(3) * orb%vec(3) * f + c_light * field%b(3) * orb%vec(1) * f

  endif
end select

end subroutine apply_hard_edge_kick

end module
