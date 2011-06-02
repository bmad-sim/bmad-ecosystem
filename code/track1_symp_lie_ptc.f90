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

real(dp) re(6), beta0, m2_rel

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
call ptc_track (fibre_ele, re, DEFAULT)  ! "track" in PTC
call vec_ptc_to_bmad (re, end%vec)

if (ele%key == wiggler$) then
  end%vec(5) = end%vec(5) - ele%value(z_patch$)
endif

call kill(fibre_ele)  ! clean up allocated memory.

! Correct map since in ptc z_ptc (vec(6)) is: 
!     z_ptc = v * t - v_ref * t_ref 
! and what is wanted to do a simple conversion to Bmad is to transform:
!     z_ptc -> v * (t - t_ref) = -z_bmad

beta0 = ele%value(p0c$) / ele%value(e_tot$)
m2_rel = (mass_of(param%particle) / ele%value(p0c$))**2

end%vec(5) = end%vec(5) + ele%value(l$) * beta0 * m2_rel * (2.0_rp * end%vec(6) + end%vec(6)**2) / &
      (((1.0_rp + end%vec(6)) / sqrt((1.0_rp + end%vec(6))**2 + m2_rel) + beta0) * ((1.0_rp + end%vec(6))**2 + m2_rel))

end subroutine
