!+
! Subroutine track1_symp_lie_ptc (start, ele, param, end)
!
! Particle tracking through a single element using a hamiltonian
! and a symplectic integrator. This uses Etienne's PTC code. For a 
! "native" BMAD version see track1_symnp_lie_bmad.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

#include "CESR_platform.inc"

subroutine track1_symp_lie_ptc (start, ele, param, end)

  use ptc_interface_mod
  use s_tracking, only: DEFAULT, alloc_fibre

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param
  type (fibre), pointer :: fibre_ele

  real*8 re(6)
  integer charge

! Construct a PTC fibre out of the ele element.
! A fibre is PTC's structure analogous to BMAD's ele_struct.  

  call alloc_fibre (fibre_ele)
  call ele_to_fibre (ele, fibre_ele, param)

  if (param%particle > 0) then
    charge = +1
  else
    charge = -1
  endif

! call the PTC routines to track through the fibre.

  call vec_bmad_to_ptc (start%vec, re)  ! convert BMAD coords to PTC coords
  call ptc_track (fibre_ele, re, DEFAULT, charge)  ! "track" in PTC
  call vec_ptc_to_bmad (re, end%vec)

  call kill(fibre_ele)  ! clean up allocated memory.

end subroutine
