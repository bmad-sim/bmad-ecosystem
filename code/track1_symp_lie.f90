!+
! Subroutine track1_symp_lie (start, ele, param, end)
!
! Particle tracking through a single element using a hamiltonian
! and a symplectic integrator.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!     %aperture_limit_on -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   end   -- Coord_struct: End position
!   param
!     %lost -- Set .true. If the particle is outside the aperture and
!                %aperture_limit_on is set. Also: %lost is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %aperture_limit_on.
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

#include "CESR_platform.inc"

subroutine track1_symp_lie (start, ele, param, end)

  use accelerator

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param
  type (fibre), pointer :: fibre_ele

  real*8 re(6)
  integer charge

!

  call alloc_fibre (fibre_ele)
  call ele_to_fibre (ele, fibre_ele, param)

  if (param%particle > 0) then
    charge = +1
  else
    charge = -1
  endif
  
  call vec_bmad_to_ptc (start%vec, re)
  call track (fibre_ele, re, DEFAULT, charge)
  call vec_ptc_to_bmad (re, end%vec)

  call kill(fibre_ele)

end subroutine
