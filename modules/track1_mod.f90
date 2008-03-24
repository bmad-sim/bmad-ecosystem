!+
! Module track1_mod
!
! Module of helper routines for track1 routines
!-

#include "CESR_platform.inc"

module track1_mod

  use bmad_struct
  use bmad_interface
  use make_mat6_mod
  use mad_mod

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine orbit_to_bpm_reading (orb, ele, axis, reading)
!
! Find the reading on a bpm given the orbit at the BPM and the 
! BPM offsets, noise, etc.
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Input: 
!  orb        -- Coord_struct: Orbit position at BPM
!  ele        -- Ele_struct: the BPM
!    %value(noise$)         -- relative bpm resolution RMS
!    %value(tilt_tot$)      -- angle error in radians rms.
!    %value(x_gain_calib$)  -- Horizontal gain correction.
!    %value(y_gain_err$)    -- Horizontal gain error.
!    ... etc ...
!  axis       -- Integer: x_plane$ or y_plane$
!
! Output:
!  reading  -- Real(rp): BPM reading
!-

subroutine orbit_to_bpm_reading (orb, ele, axis, reading)

use random_mod

implicit none

type (coord_struct) orb
type (ele_struct) ele

real(rp) reading, x_noise_factor, y_noise_factor
real(rp) ran_num(2), angle, x, y, x_gain, y_gain

integer axis

character(20) :: r_name = "orbit_to_bpm_reading"

logical err

  if (.not. ele%is_on) then
    reading = 0.0
    return
  endif

  if (ele%key /= monitor$ .and. ele%key /= instrument$ .and. ele%key /= marker$) then
    reading = 0.0
    return
  endif

  if (ele%value(noise$) /= 0) then
    call ran_gauss (ran_num)
    x_noise_factor = ele%value(noise$) * ran_num(1) 
    y_noise_factor = ele%value(noise$) * ran_num(2)
  else
    x_noise_factor = 0
    y_noise_factor = 0
  endif

  x_gain = 1 + ele%value(x_gain_err$) - ele%value(x_gain_calib$)
  y_gain = 1 + ele%value(y_gain_err$) - ele%value(y_gain_calib$)
  x = orb%vec(1) - ele%value(x_offset_tot$) + ele%value(x_offset_calib$)
  y = orb%vec(3) - ele%value(y_offset_tot$) + ele%value(y_offset_calib$)

  if (axis == x_plane$) then
    angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) + &
                  (ele%value(crunch$) - ele%value(crunch_calib$))
    reading = x_noise_factor + x_gain * (x * cos(angle) + y * sin(angle))
  elseif (axis == y_plane$) then
    angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) - &
                  (ele%value(crunch$) - ele%value(crunch_calib$))
    reading = y_noise_factor + y_gain * (-x * sin(angle) + y * cos(angle))
  else
    reading = 0.0
    call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
    call err_exit
  endif
    
end subroutine orbit_to_bpm_reading

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, at, param)
!
! Subroutine to check if an orbit is outside the aperture.
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
!
! Output:
!   param -- lat_param_struct: Parameter structure:
!     %lost -- Set True if the orbit is outside the aperture. 
!              Note: %lost is NOT set False if the orbit is inside 
!                the aperture.
!     %plane_lost_at -- Integer: Plane where particle is lost:
!                     x_plane$ or y_plane$
!-

subroutine check_aperture_limit (orb, ele, at, param)

  implicit none

  type (coord_struct) :: orb
  type (coord_struct) orb2 
  type (ele_struct) :: ele
  type (lat_param_struct), intent(inout) :: param

  real(rp) x_lim, y_lim, x_beam, y_beam, s_here
  integer at

!

  if (ele%offset_moves_aperture) then
    orb2 = orb
    s_here = 0
    if (at == exit_end$) s_here = ele%value(l$)
    call offset_particle (ele, param, orb2, set$, set_canonical = .false., &
                 set_tilt = .false., set_multipoles = .false., set_hvkicks = .false., &
                 s_pos = s_here)
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

  if (x_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    x_lim = bmad_com%max_aperture_limit

  if (y_beam < 0) then
    y_lim = ele%value(y1_limit$)
  else
    y_lim = ele%value(y2_limit$)
  endif

  if (y_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    y_lim = bmad_com%max_aperture_limit

  if (x_lim == 0 .and. y_lim == 0) return

  if (ele%key == ecollimator$) then
    if (x_lim == 0 .or. y_lim == 0) then
      print *, 'ERROR IN CHECK_APERTURE_LIMIT: ECOLLIMATOR HAS ONE LIMIT ZERO'
      print *, '      AND THE OTHER NOT: ', ele%name
      call err_exit
    endif
    if ((x_beam / x_lim)**2 + (y_beam / y_lim)**2 > 1) then
      param%lost = .true.
      if (abs(x_beam / x_lim) > abs(y_beam / y_lim)) then
        param%plane_lost_at = x_plane$
      else
        param%plane_lost_at = y_plane$
      endif
    endif
  else
    if (abs(x_beam) > x_lim) then
      param%lost = .true.
      param%plane_lost_at = x_plane$
    endif
    if (abs(y_beam) > y_lim) then
      param%lost = .true.
      param%plane_lost_at = y_plane$
    endif
  endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_drift (orb, length)
!
! Subroutine to track through a drift.
!
! Modules needed:
!   use precision_def
!
! Input:
!   orb(6) -- Real(rp): Orbit at start of the drift
!   length -- Real(rp): Length of drift.
!
! Output:
!   orb(6) -- Real(rp): Orbit at end of the drift
!-

subroutine track_a_drift (orb, length)

  implicit none

  real(rp) orb(6), length, rel_E

  rel_E = 1 + orb(6)

  orb(1) = orb(1) + length * orb(2) / rel_E
  orb(3) = orb(3) + length * orb(4) / rel_E
  orb(5) = orb(5) - length * (orb(2)**2 + orb(4)**2) / (2 * rel_E**2)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_bend (start, ele, param, end)
!
! Particle tracking through a bend element. 
! For e1 or e2 non-zero this subroutine treats the edges as thin quads.
! This routine is not meant for general use.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position.
!   ele    -- Ele_struct: Bend element.
!   param  -- lat_param_struct: Lattice parameters.
!
! Output:
!   end     -- Coord_struct: End position.
!-

subroutine track_a_bend (start, ele, param, end)

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (lat_param_struct), intent(inout) :: param

  real(rp) b1, angle, ct, st, x, px, y, py, z, pz, dpx_t
  real(rp) rel_p, rel_p2, Dy, px_t, factor, rho, g, g_err
  real(rp) length, g_tot, dP, eps, pxy2, f
  real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, z_2, s_y, c_y, r(6)

! simple case

  if (ele%value(g$) == 0) then
    end = start
    end%vec(2) = end%vec(2) - length * ele%value(g_err$) / 2
    call track_a_drift (end%vec, ele%value(l$))
    end%vec(2) = end%vec(2) - length * ele%value(g_err$) / 2
    return
  endif

!-----------------------------------------------------------------------

  end = start
  call offset_particle (ele, param, end, set$, set_canonical = .false.)
  call track_bend_edge (end, ele, .true., .false.)

  length = ele%value(l$)
  g = ele%value(g$)
  g_err = ele%value(g_err$)
  g_tot = g + g_err
  b1 = g_tot
  angle = ele%value(g$) * length
  rho = 1 / g
  dP = start%vec(6)
  rel_p  = 1 + dP
  rel_p2 = rel_p**2
  k_1 = ele%value(k1$)

! with k1 /= 0 use small angle approximation

  if (k_1 /= 0) then

    call sbend_body_with_k1_map (g, g_err, length, k_1, end%vec, end = end%vec)

!-----------------------------------------------------------------------
! Track through main body...
! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  else

    ct = cos(angle)
    st = sin(angle)

    x  = end%vec(1)
    px = end%vec(2)
    y  = end%vec(3)
    py = end%vec(4)
    z  = end%vec(5)
    pz = end%vec(6)
 
    pxy2 = px**2 + py**2
    if (rel_p2 - pxy2 < 0.1) then  ! somewhat arbitrary cutoff
      param%lost = .true.
      param%plane_lost_at = x_plane$
      end%vec(1) = 2 * bmad_com%max_aperture_limit
      end%vec(3) = 2 * bmad_com%max_aperture_limit
      return
    endif 

! The following is to make sure that a beam entering on-axis remains *exactly* on-axis.

    if (pxy2 < 1e-5) then  
      f = pxy2 / (2 * rel_p)
      f = dP - f - f*f/2 - g_err*rho - b1*x
    else
      f = sqrt(rel_p2 - pxy2) - 1 - g_err*rho - b1*x
    endif

    Dy  = sqrt(rel_p2 - py**2)
    px_t = px*ct + f*st
    dpx_t = -px*st/rho + f*ct/rho

    if (abs(px) > Dy .or. abs(px_t) > Dy) then
      param%lost = .true.
      param%plane_lost_at = x_plane$
      return
    endif    

    factor = (asin(px/Dy) - asin(px_t/Dy)) / b1

    eps = px_t**2 + py**2
    if (eps < 1e-5 * rel_p2 ) then  ! use small angle approximation
      eps = eps / (2 * rel_p)
      end%vec(1) = (dP - g_err / g - rho*dpx_t + eps * (eps / (2 * rel_p) - 1)) / b1
    else
      end%vec(1) = (sqrt(rel_p2 - eps) - rho*dpx_t - rho*b1) / b1
    endif

    end%vec(2) = px_t
    end%vec(3) = y + py * (angle/b1 + factor)
    end%vec(4) = py
    end%vec(5) = end%vec(5) + length * (g_err - g*dP) / g_tot - rel_p * factor
    end%vec(6) = pz

  endif

! Track through the exit face. Treat as thin lens.

  call track_bend_edge (end, ele, .false., .false.)
  call offset_particle (ele, param, end, unset$, set_canonical = .false.)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_bend_edge (orb, e, g, start_edge, reverse, kx, ky)
!
! Subroutine to track through the edge field of a bend.
! This routine is not meant for general use.
!-

subroutine track_bend_edge (orb, ele, start_edge, reverse, kx, ky)

  type (ele_struct) ele
  type (coord_struct) orb
  real(rp), optional :: kx, ky
  real(rp) e, g_tot, del
  logical start_edge, reverse

! Track through the entrence face. Treat as thin lens.

  g_tot = ele%value(g$) + ele%value(g_err$)
  if (reverse) g_tot = -g_tot

  if (start_edge) then
    e = ele%value(e1$)
    del = tan(e) * g_tot
    if (present(kx)) kx = del 
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    if (ele%value(fint$) /= 0) del = g_tot * tan(e - 2 * ele%value(fint$) * &
                      abs(g_tot) * ele%value(hgap$) *  (1 + sin(e)**2) / cos(e))
    if (present(ky)) ky = -del
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

  else
    e = ele%value(e2$)
    del = tan(e) * g_tot
    if (present(ky)) kx = del
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    if (ele%value(fintx$) /= 0) del = g_tot * tan(e - 2 * ele%value(fintx$) * &
                      abs(g_tot) * ele%value(hgapx$) *  (1 + sin(e)**2) / cos(e))
    if (present(ky)) ky = -del
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

  endif

  if (reverse) then
    if (present(kx)) kx = -kx
    if (present(ky)) ky = -ky
  endif

end subroutine

end module
