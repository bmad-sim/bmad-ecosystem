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

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, param)
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
!-

subroutine check_aperture_limit (orb, ele, param)

  implicit none

  type (coord_struct), intent(in) :: orb
  type (ele_struct),   intent(in)  :: ele
  type (param_struct), intent(inout) :: param

  real(rp) x_lim, y_lim, x_beam, y_beam, l2

!

  x_lim = ele%value(x_limit$)
  if (x_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    x_lim = bmad_com%max_aperture_limit

  y_lim = ele%value(y_limit$)
  if (y_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    y_lim = bmad_com%max_aperture_limit

  if (x_lim == 0 .and. y_lim == 0) return

  l2 = ele%value(l$) / 2
  x_beam = orb%vec(1) - ele%value(x_offset_tot$)  ! - ele%value(x_pitch$) * l2
  y_beam = orb%vec(3) - ele%value(y_offset_tot$)  ! - ele%value(y_pitch$) * l2

  if (ele%key == ecollimator$) then
    if (x_lim == 0 .or. y_lim == 0) then
      print *, 'ERROR IN CHECK_APERTURE_LIMIT: ECOLLIMATOR HAS ONE LIMIT ZERO'
      print *, '      AND THE OTHER NOT: ', ele%name
      call err_exit
    endif
    if ((x_beam / x_lim)**2 + (y_beam / y_lim)**2 > 1) param%lost = .true.
  else
    if (abs(x_beam) > x_lim) param%lost = .true.
    if (abs(y_beam) > y_lim) param%lost = .true.
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
  type (coord_struct) start2

  real(rp) b1, angle, ct, st, x, px, y, py, z, pz, dpx_t, phi
  real(rp) rel_E, rel_E2, Dxy, Dy, px_t, factor, rho, fact, g0, dg
  real(rp) length, g_tot, dE, del, f, k1, kc, mat2(2,2), mat_i6(6)

! simple case

  if (ele%value(g$) == 0) then
    end = start
    end%vec(2) = end%vec(2) + length * ele%value(delta_g$) / 2
    call track_a_drift (end%vec, ele%value(l$))
    end%vec(2) = end%vec(2) + length * ele%value(delta_g$) / 2
    return
  endif

!

  end = start
  call offset_particle (ele, param, end, set$, set_canonical = .false.)

  length = ele%value(l$)
  g0 = ele%value(g$)
  dg = ele%value(delta_g$)
  g_tot = g0 + dg
  b1 = g_tot
  angle = ele%value(angle$)
  rho = 1 / g0
  dE = start%vec(6)
  rel_E  = 1 + dE
  rel_E2 = rel_E**2

  call track_bend_edge (end, ele%value(e1$), g_tot, .true.)

!-----------------------------------------------------------------------
! Track through main body...

  k1 = ele%value(k1$)
  select case (k1 == 0)

! For k1 == 0 then use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  case (.true.)  ! k1 == 0

    ct = cos(angle)
    st = sin(angle)

    x  = end%vec(1)
    px = end%vec(2)
    y  = end%vec(3)
    py = end%vec(4)
    z  = end%vec(5)
    pz = end%vec(6)
 
    Dxy = sqrt(rel_E2 - px**2 - py**2)
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

!  If k1 /= 0 then just use the 1st order transfer matrix.

  case (.false.)

    end%vec(2) = end%vec(2) / rel_E  ! convert px to x'
    end%vec(4) = end%vec(4) / rel_E  ! convert py to y'

    k1 = k1 / rel_E
    kc = g_tot**2 + k1
    start2 = end  ! Save coords after entrence face

    call quad_mat2_calc (-kc, length, mat2)
    end%vec(1:2) = matmul (mat2, end%vec(1:2))

    call quad_mat2_calc (k1, length, mat2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))

    phi = sqrt(abs(kc)) * length
    if (kc < 0) then
      mat_i6(1) = (1 - cosh(phi)) * g_tot / kc
      mat_i6(2) = sinh(phi) * g_tot / sqrt(-kc)
      mat_i6(5) = (phi - sinh(phi)) * g_tot**2 / abs(kc)**1.5
    else
      mat_i6(1) = (1 - cos(phi)) * g_tot / kc
      mat_i6(2) = sin(phi) * g_tot / sqrt(kc)
      mat_i6(5) = (sin(phi) - phi) * g_tot**2 / kc**1.5
    endif

    end%vec(1) = end%vec(1) + mat_i6(1) * dE
    end%vec(2) = end%vec(2) + mat_i6(2) * dE
    end%vec(5) = end%vec(5) + mat_i6(5) * dE

    fact = start2%vec(4) * (mat_i6(5) + length)

    end%vec(5) = end%vec(5) - &
                    mat_i6(2) * start2%vec(1) - mat_i6(1) * start2%vec(2) + &
                    fact * start2%vec(4) / 2
    end%vec(3) = end%vec(3) + fact * dE

    end%vec(2) = end%vec(2) * rel_E
    end%vec(4) = end%vec(4) * rel_E

  end select
  
!------------------------------------------------------------------------
! Track through the exit face. Treat as thin lens.

  call track_bend_edge (end, ele%value(e2$), g_tot, .false.)

  call offset_particle (ele, param, end, unset$, set_canonical = .false.)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_bend_edge (orb, e, g, start_edge)
!
! Subroutine to track through the edge field of a bend.
! This routine is not meant for general use.
!-

subroutine track_bend_edge (orb, e, g, start_edge)

  type (coord_struct) orb
  real(rp) e, g, del, f
  logical start_edge

! Track through the entrence face. Treat as thin lens.
! The second order terms come from the Hamiltonian term:
!       H = (g * sec^2(e1$) / 2) * p_x * y^2
! (See the MAD8 Physics writeup)

  if (start_edge) then
    del = tan(e) * g
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

!    f = g / (2 * (1 - e**2))
!    orb%vec(1) = orb%vec(1) + f * orb%vec(3)**2
!    orb%vec(4) = orb%vec(4) - 2 * f * orb%vec(2) * orb%vec(3)

  else
!    f = g / (2 * (1 - e**2))
!    orb%vec(1) = orb%vec(1) - f * orb%vec(3)**2
!    orb%vec(4) = orb%vec(4) + 2 * f * orb%vec(2) * orb%vec(3)

    del = tan(e) * g
    orb%vec(2) = orb%vec(2) + del * orb%vec(1)
    orb%vec(4) = orb%vec(4) - del * orb%vec(3)

  endif

end subroutine

end module
