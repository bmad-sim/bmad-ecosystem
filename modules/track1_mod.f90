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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, param, plane_lost)
!
! Subroutine to check if an orbit is outside the aperture.
!
! Modules needed:
!   use bmad
!
! Input:
!   orb -- Coord_struct: coordinates of a particle.
!   ele -- Ele_struct: Element holding the aperture
!     %value(x_limit$) -- Horizontal aperture.
!     %value(y_limit$) -- Vertical aparture.
!   param -- Param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. In this case param%lost will
!               be set to True.
!
! Output:
!   param -- Param_struct: Parameter structure:
!     %lost -- Set True if the orbit is outside the aperture. 
!              Note: %lost is NOT set False if the orbit is inside 
!                the aperture.
!   plane_lost -- Integer: Plane where particle is lost:
!                     x_plane$ or y_plane$
!-

subroutine check_aperture_limit (orb, ele, param, plane_lost)

  implicit none

  type (coord_struct), intent(in) :: orb
  type (ele_struct),   intent(in)  :: ele
  type (param_struct), intent(inout) :: param

  real(rp) x_lim, y_lim, x_beam, y_beam, l2
  integer plane_lost

!

  x_lim = ele%value(x_limit$)
  if (x_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    x_lim = bmad_com%max_aperture_limit

  y_lim = ele%value(y_limit$)
  if (y_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    y_lim = bmad_com%max_aperture_limit

  if (x_lim == 0 .and. y_lim == 0) return

  l2 = ele%value(l$) / 2
  x_beam = orb%vec(1)
  y_beam = orb%vec(3)

  if (ele%key == ecollimator$) then
    if (x_lim == 0 .or. y_lim == 0) then
      print *, 'ERROR IN CHECK_APERTURE_LIMIT: ECOLLIMATOR HAS ONE LIMIT ZERO'
      print *, '      AND THE OTHER NOT: ', ele%name
      call err_exit
    endif
    if ((x_beam / x_lim)**2 + (y_beam / y_lim)**2 > 1) then
      param%lost = .true.
      if (abs(x_beam / x_lim) > abs(y_beam / y_lim)) then
        plane_lost = x_plane$
      else
        plane_lost = y_plane$
      endif
    endif
  else
    if (abs(x_beam) > x_lim) then
      param%lost = .true.
      plane_lost = x_plane$
    endif
    if (abs(y_beam) > y_lim) then
      param%lost = .true.
      plane_lost = y_plane$
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
!   param  -- Param_struct: Lattice parameters.
!
! Output:
!   end     -- Coord_struct: End position.
!-

subroutine track_a_bend (start, ele, param, end)

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  real(rp) b1, angle, ct, st, x, px, y, py, z, pz, dpx_t
  real(rp) rel_E, rel_E2, Dxy, Dy, px_t, factor, rho, g0, dg
  real(rp) length, g_tot, dE, re_xy

! simple case

  if (ele%value(g$) == 0) then
    end = start
    end%vec(2) = end%vec(2) + length * ele%value(delta_g$) / 2
    call track_a_drift (end%vec, ele%value(l$))
    end%vec(2) = end%vec(2) + length * ele%value(delta_g$) / 2
    return
  endif

! with k1 == 0 use the MAD 2nd order map.

  if (ele%value(k1$) /= 0) then
    call track1_mad (start, ele, param, end)
    return
  endif

!-----------------------------------------------------------------------

  end = start
  call offset_particle (ele, param, end, set$, set_canonical = .false.)

  length = ele%value(l$)
  g0 = ele%value(g$)
  dg = ele%value(delta_g$)
  g_tot = g0 + dg
  b1 = g_tot
  angle = ele%value(g$) * length
  rho = 1 / g0
  dE = start%vec(6)
  rel_E  = 1 + dE
  rel_E2 = rel_E**2

  call track_bend_edge (end, ele, .true., .false.)

! Track through main body...
! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  ct = cos(angle)
  st = sin(angle)

  x  = end%vec(1)
  px = end%vec(2)
  y  = end%vec(3)
  py = end%vec(4)
  z  = end%vec(5)
  pz = end%vec(6)
 
  re_xy = rel_E2 - px**2 - py**2
  if (re_xy < 0.1) then  ! somewhat arbitrary cutoff
    param%lost = .true.
    end%vec(1) = 2 * bmad_com%max_aperture_limit
    end%vec(3) = 2 * bmad_com%max_aperture_limit
    return
  endif 

  Dxy = sqrt(re_xy)
  Dy  = sqrt(rel_E2 - py**2)

  px_t = px*ct + (Dxy - b1*(rho+x))*st
  dpx_t = -px*st/rho + (Dxy - b1*(rho+x))*ct/rho
  factor = (asin(px/Dy) - asin(px_t/Dy)) / b1

  end%vec(1) = (sqrt(rel_E2 - px_t**2 -py**2) - rho*dpx_t - rho*b1) / b1
  end%vec(2) = px_t
  end%vec(3) = y + py * (angle/b1 + factor)
  end%vec(4) = py
  end%vec(5) = end%vec(5) + length * (dg - g0*dE) / g_tot - rel_E * factor
  end%vec(6) = pz

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
  real(rp) e, g, del
  logical start_edge, reverse

! Track through the entrence face. Treat as thin lens.

  g = ele%value(g$) + ele%value(delta_g$)
  if (reverse) g = -g

  if (start_edge) then
    e = ele%value(e1$)
    del = tan(e) * g
    if (present(kx)) kx = del 
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    if (ele%value(fint$) /= 0) del = g * tan(e - 2 * ele%value(fint$) * &
                      abs(g) * ele%value(hgap$) *  (1 + sin(e)**2) / cos(e))
    if (present(ky)) ky = -del
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

  else
    e = ele%value(e2$)
    del = tan(e) * g
    if (present(ky)) kx = del
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    if (ele%value(fintx$) /= 0) del = g * tan(e - 2 * ele%value(fintx$) * &
                      abs(g) * ele%value(hgapx$) *  (1 + sin(e)**2) / cos(e))
    if (present(ky)) ky = -del
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

  endif

  if (reverse) then
    if (present(kx)) kx = -kx
    if (present(ky)) ky = -ky
  endif

end subroutine

end module
