!+
! Subroutine track1_symp_lie_ptc (orbit, ele, param, track)
!
! Particle tracking through a single element using a hamiltonian
! and a symplectic integrator. This uses Etienne's PTC code. For a 
! "native" BMAD version see track1_symp_lie_bmad.
!
! Input:
!   orbit      -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   orbit      -- Coord_struct: End position
!   track      -- Track_struct, optional: Structure holding the track information.
!-

subroutine track1_symp_lie_ptc (orbit, ele, param, track)

use ptc_interface_mod, except_dummy => track1_symp_lie_ptc
use ptc_spin, only: probe, assignment(=), operator(+), SPIN0, TOTALPATH0, &
                                              track_probe, track_probe_x, CONVERSION_XPRIME_IN_ABELL
use s_tracking, only: alloc_fibre, integration_node, check_stable
use mad_like, only: fibre, kill

implicit none

type (coord_struct) :: orbit, start_orb
type (track_struct), optional :: track
type (ele_struct) :: ele, drift_ele
type (lat_param_struct) :: param
type (fibre), pointer :: ptc_fibre
type (probe) ptc_probe
type (integration_node), pointer :: ptc_track
type (internal_state) state

real(dp) re(6)
integer i, stm
logical err_flag

character(20) :: r_name = 'track1_symp_lie_ptc'

! call the PTC routines to track through the fibre.

CONVERSION_XPRIME_IN_ABELL = (.not. bmad_com%convert_to_kinetic_momentum) ! Only affects cylindrical map eles
state = ptc_private%base_state

re = orbit%vec

!-----------------------------
! track element

start_orb = orbit   ! Save initial state

call ele_to_fibre (ele, ptc_fibre, .true., err_flag, ref_in = orbit)
if (err_flag) then
  orbit%state = lost$
  return
endif

!

stm = ele%spin_tracking_method
if (bmad_com%spin_tracking_on .and. (stm == tracking$ .or. stm == symp_lie_ptc$) .or. present(track)) then
  if (bmad_com%spin_tracking_on) state = state + spin0

  ptc_probe = re
  ptc_probe%q%x = [1, 0, 0, 0]

  if (present(track)) then
    ptc_track => ptc_fibre%t1
    call save_this_step(track, ptc_probe, ele)

    do while (.not. associated(ptc_track, ptc_fibre%t2))
      call track_probe (ptc_probe, state, node1 = ptc_track, node2 = ptc_track%next)
      call save_this_step(track, ptc_probe, ele)
      if (.not. check_stable) exit
      ptc_track => ptc_track%next
    enddo

    call track_probe (ptc_probe, state, node1 = ptc_track, node2 = ptc_track%next)
    call save_this_step(track, ptc_probe, ele)

  else
    call track_probe (ptc_probe, state, fibre1 = ptc_fibre)
  endif

  orbit%spin = quat_rotate(ptc_probe%q%x, start_orb%spin)
  re = ptc_probe%x

else
  ! Orignally used track (ptc_fibre, re, state) but this will not track taylor elements correctly.
  call track_probe_x (re, state, fibre1 = ptc_fibre)
endif

!-----------------------------

orbit%vec = re

! 

if (ele%value(p0c$) /= ele%value(p0c_start$) .or. start_orb%vec(6) /= orbit%vec(6)) then
  call convert_pc_to (ele%value(p0c$) * (1 + orbit%vec(6)), orbit%species, beta = orbit%beta)
endif

orbit%s = ele%s
orbit%p0c = ele%value(p0c$)

if (state%totalpath == 1) then
  orbit%t = start_orb%t + start_orb%vec(5) / (start_orb%beta * c_light) - orbit%vec(5) / (orbit%beta * c_light)
else
  orbit%t = start_orb%t + ele%value(delta_ref_time$) + &
                          start_orb%vec(5) / (start_orb%beta * c_light) - orbit%vec(5) / (orbit%beta * c_light)
endif

if (.not. check_stable) orbit%state = lost$

CONVERSION_XPRIME_IN_ABELL = .true. ! Reset to normal.

!---------------------------------------------------------------------
contains

subroutine save_this_step(track, ptc_probe, ele)

type (track_struct) track
type (probe) ptc_probe
type (ele_struct) ele
type (coord_struct) orbit
real(dp) re(6)

! The complication is that PTC pz is the true canonical momentum which includes the electrostatic potential.
! But Bmad pz does not include the electrostatic potential.
! Note: ptc_track%s(2) is the same as %s(1) except in a true rbend where %s(2) is the chord distance.

orbit = start_orb
orbit%vec = ptc_probe%x
orbit%s = ptc_track%s(1) + ele%s_start
orbit%spin = quat_rotate(ptc_probe%q%x, start_orb%spin)
if (.not. check_stable) orbit%state = lost$
call save_a_step (track, ele, param, .false., orbit, ptc_track%s(1))

end subroutine save_this_step

end subroutine track1_symp_lie_ptc
