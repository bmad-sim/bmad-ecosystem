!+
! Subroutine CHANGE_BASIS (COORD, REF_ENERGY, REF_Z, TO_CART, TIME_DISP)
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
!     COORD      -- Coord_struct: Coordinates of particle
!     REF_ENERGY -- Real(rdef): Reference energy of beam
!     REF_Z      -- Real(rdef): Reference longitudinal position of beam
!     TO_CART    -- Logical: True if converting to cartesian coordinates
!                            False if converting to accelerator coordinates
!     If to_cart == .false.:
!       TIME_DISP -- Real(rdef): Time displacement of particle
!
! Output:
!     COORD  -- Coord_struct: Converted coordinates
!     If to_cart == .true.:
!       TIME_DISP -- Real(rdef): Time displacement of particle
!-

!$Id$
!$Log$
!Revision 1.6  2003/01/27 14:40:31  dcs
!bmad_version = 56
!
!Revision 1.5  2002/06/13 14:54:23  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:12  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:37  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct)  coord, saved
  real(rdef) beta2, ref_energy, ref_z, time_disp, rvar
  logical to_cart

!

  saved%x%pos = coord%x%pos
  saved%y%pos = coord%y%pos

  if (to_cart) then
    rvar = (e_mass / (ref_energy * (coord%z%vel+1)))**2
    if (rvar <= 0.001) then
      saved%z%vel = c_light / sqrt(1 + coord%x%vel**2 + coord%y%vel**2)  &
        * (1 - rvar/2)
    else
      saved%z%vel = c_light * sqrt((1 - rvar)  &
        / (1 + coord%x%vel**2 + coord%y%vel**2))
    endif
    saved%x%vel = coord%x%vel * saved%z%vel
    saved%y%vel = coord%y%vel * saved%z%vel
!         Convert z from a relative position to an absolute position:
    saved%z%pos = ref_z
    time_disp =  - coord%z%pos / saved%z%vel
  else
    saved%x%vel = coord%x%vel / coord%z%vel
    saved%y%vel = coord%y%vel / coord%z%vel
    beta2 = (coord%x%vel**2 + coord%y%vel**2 + coord%z%vel**2) / c_light**2
!         Prevent sqrt(negative number) due to rounding errors.  Set such
!         particles to have gamma = 1000 (approximately).
    if (beta2 >= 0.999999) beta2 = 0.999999
    saved%z%vel = e_mass / (ref_energy * sqrt(1 - beta2)) - 1
!         Convert z from an absolute position to a relative position:
    saved%z%pos = - coord%z%vel * time_disp
  endif
  coord = saved

  return
  end
