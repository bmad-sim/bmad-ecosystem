!+
! Subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)
!
!   Subroutine to convert accelerator coordinates (x, x', y, y', z, z') to
! cartesian coordinates and time derivatives (x, x_dot, y, y_dot, z, z_dot) or
! vice-versa.  The conversions assume no small-angle nor other approximations
! whatsoever EXCEPT that there is no (local) bend--the design orbit is locally
! straight.
!   Note that if the Lorentz factor (gamma) is around 100 or higher, converting
! _from_ accelerator coordinates and then back _to_ accelerator coordinates will
! produce a slight change in z'.
!   Also, due to rounding errors in iterations, particles may acquire a speed
! higher than c.  Such speeds will be reset (if TO_CART is false) so that the
! Lorentz factor is about 1000.
!  These two problems are corrected if coil_track is the calling routine.
! -- Created by Daniel Fromowitz, November 1998.
!
! Modules Needed:
!   use bmad
!
! Input:
!   coord      -- Coord_struct: Coordinates of particle
!   ref_energy -- Real(rp): Reference energy of beam
!   ref_z      -- Real(rp): Reference longitudinal position of beam
!   to_cart    -- Logical: True if converting to cartesian coordinates
!                          False if converting to accelerator coordinates
!   If to_cart == .false.:
!     time_disp -- Real(rp): Time displacement of particle
!
! Output:
!   coord  -- Coord_struct: Converted coordinates
!   If to_cart == .true.:
!     time_disp -- Real(rp): Time displacement of particle
!-


#include "CESR_platform.inc"

subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct)  coord, saved
  real(rp) beta2, ref_energy, ref_z, time_disp, rvar
  logical to_cart

!

  saved%vec(1) = coord%vec(1)
  saved%vec(3) = coord%vec(3)

  if (to_cart) then
    rvar = (m_electron / (ref_energy * (coord%vec(6)+1)))**2
    if (rvar <= 0.001) then
      saved%vec(6) = c_light / sqrt(1 + coord%vec(2)**2 + coord%vec(4)**2)  &
        * (1 - rvar/2)
    else
      saved%vec(6) = c_light * sqrt((1 - rvar)  &
        / (1 + coord%vec(2)**2 + coord%vec(4)**2))
    endif
    saved%vec(2) = coord%vec(2) * saved%vec(6)
    saved%vec(4) = coord%vec(4) * saved%vec(6)
!         Convert z from a relative position to an absolute position:
    saved%vec(5) = ref_z
    time_disp =  - coord%vec(5) / saved%vec(6)
  else
    saved%vec(2) = coord%vec(2) / coord%vec(6)
    saved%vec(4) = coord%vec(4) / coord%vec(6)
    beta2 = (coord%vec(2)**2 + coord%vec(4)**2 + coord%vec(6)**2) / c_light**2
!         Prevent sqrt(negative number) due to rounding errors.  Set such
!         particles to have gamma = 1000 (approximately).
    if (beta2 >= 0.999999) beta2 = 0.999999
    saved%vec(6) = m_electron / (ref_energy * sqrt(1 - beta2)) - 1
!         Convert z from an absolute position to a relative position:
    saved%vec(5) = - coord%vec(6) * time_disp
  endif
  coord = saved

end subroutine
