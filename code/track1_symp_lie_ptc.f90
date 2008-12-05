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
!   param  -- lat_param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

#include "CESR_platform.inc"

subroutine track1_symp_lie_ptc (start, ele, param, end)

use ptc_interface_mod, except_dummy => track1_symp_lie_ptc
use s_tracking, only: DEFAULT, alloc_fibre
use mad_like, only: kill, ptc_track => track

implicit none

type (coord_struct) :: start
type (coord_struct) :: end
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (fibre), pointer :: fibre_ele

real(dp) re(6)

character(20) :: r_name = 'track1_symp_lie_ptc'

! Construct a PTC fibre out of the ele element.
! A fibre is PTC's structure analogous to BMAD's ele_struct.  

call alloc_fibre (fibre_ele)
call ele_to_fibre (ele, fibre_ele, param, .true.)

! call the PTC routines to track through the fibre.

if (ele%key == wiggler$ .and. ele%value(z_patch$) == 0) then
  call out_io (s_fatal$, r_name, 'WIGGLER Z_PATCH VALUE HAS NOT BEEN COMPUTED!')
  call err_exit 
endif

call vec_bmad_to_ptc (start%vec, re)  ! convert BMAD coords to PTC coords
call ptc_track (fibre_ele, re, DEFAULT, +1)  ! "track" in PTC
call vec_ptc_to_bmad (re, end%vec)

if (ele%key == wiggler$) then
  end%vec(5) = end%vec(5) - ele%value(z_patch$)
endif

call kill(fibre_ele)  ! clean up allocated memory.

end subroutine
