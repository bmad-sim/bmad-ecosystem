!+
! Subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb)
!
! Particle tracking through a single element using a hamiltonian
! and a symplectic integrator. This uses Etienne's PTC code. For a 
! "native" BMAD version see track1_symnp_lie_bmad.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct:
!
! Output:
!   end_orb   -- Coord_struct: End position
!-

#include "CESR_platform.inc"

subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb)

use ptc_interface_mod, except_dummy => track1_symp_lie_ptc
use s_tracking, only: DEFAULT, alloc_fibre
use mad_like, only: fibre, kill, ptc_track => track

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (fibre), pointer :: fibre_ele

real(dp) re(6), beta0

character(20) :: r_name = 'track1_symp_lie_ptc'

! Error test

if (ele%key == wiggler$ .and. ele%value(z_patch$) == 0) then
  call out_io (s_fatal$, r_name, 'WIGGLER Z_PATCH VALUE HAS NOT BEEN COMPUTED!')
  call err_exit 
endif

! Construct a PTC fibre out of the ele element.
! A fibre is PTC's structure analogous to BMAD's ele_struct.  

call ele_to_fibre (ele, fibre_ele, param, .true.)

! call the PTC routines to track through the fibre.

if (ele_has_constant_reference_energy(ele)) then
  beta0 = ele%value(p0c$) / ele%value(e_tot$)
else
  beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
endif

call vec_bmad_to_ptc (start_orb%vec, beta0, re)
call ptc_track (fibre_ele, re, DEFAULT)  ! "track" in PTC
call vec_ptc_to_bmad (re, beta0, end_orb%vec)
if (.not. ele_has_constant_reference_energy(ele)) &
      call vec_bmad_ref_energy_correct(end_orb%vec, ele%value(p0c$) / ele%value(p0c_start$))

if (ele%key == wiggler$) then
  end_orb%vec(5) = end_orb%vec(5) - ele%value(z_patch$)
endif

end subroutine
