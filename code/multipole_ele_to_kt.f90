!+
! Subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)
!
! Subroutine to put the multipole components (strength and tilt)
! into 2 vectors along with the appropriate scaling.
! Note: tilt(:) does includes ele%value(tilt$).
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele          -- Ele_struct: Multipole element.
!   particle     -- Integer: Particle species (+1 or -1).
!   use_ele_tilt -- Logical: If True then include ele%value(tilt$) 
!                     in calculations.
!
! Output:
!   knl(0:n_pole_maxx)  -- Real: Vector of strengths, MAD units.
!   tilt(0:n_pole_maxx) -- Real: Vector of tilts.
!-

!$Id$
!$Log$
!Revision 1.1  2002/01/08 21:48:14  dcs
!Align with VMS version
!

#include "CESR_platform.inc"

subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct)  ele

  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx), signn, a_n, b_n
  real value(n_attrib_maxx), a(0:n_pole_maxx), b(0:n_pole_maxx)

  integer n, particle, n_fact

  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on) then
    knl = 0
    tilt = 0
    return
  endif

! multipole
                    
  if (ele%key == multipole$) then
    knl  = ele%value(ix1_m$:ix2_m$-1:2)
    tilt = ele%value(ix1_m$+1:ix2_m$:2) + ele%value(tilt$)
    return
  endif

! ab_multiple, etc...

  call multipole_ele_to_ab (ele, particle, a, b, .false.)
  call multipole_ab_to_kt (a, b, knl, tilt)
  if (use_ele_tilt) tilt = tilt + ele%value(tilt$)
                       
end subroutine
