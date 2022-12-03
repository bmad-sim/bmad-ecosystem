!+
! Subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb, track)
!
! Particle tracking through a single element using a hamiltonian
! and a symplectic integrator. This uses Etienne's PTC code. For a 
! "native" BMAD version see track1_symp_lie_bmad.
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   end_orb    -- Coord_struct: End position
!   track      -- Track_struct, optional: Structure holding the track information.
!-

subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb, track)

use ptc_interface_mod, except_dummy => track1_symp_lie_ptc
use ptc_spin, only: probe, assignment(=), operator(+), SPIN0, TOTALPATH0, &
                                              track_probe, track_probe_x, CONVERSION_XPRIME_IN_ABELL
use s_tracking, only: alloc_fibre, integration_node, check_stable
use mad_like, only: fibre, kill

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
type (track_struct), optional :: track
type (ele_struct) :: ele, drift_ele
type (lat_param_struct) :: param
type (fibre), pointer :: ptc_fibre
type (probe) ptc_probe
type (integration_node), pointer :: ptc_track
type (internal_state) state, state0

real(dp) re(6)
integer i, stm
logical err_flag

character(20) :: r_name = 'track1_symp_lie_ptc'

! call the PTC routines to track through the fibre.

CONVERSION_XPRIME_IN_ABELL = (.not. bmad_com%convert_to_kinetic_momentum) ! Only affects cylindrical map eles

STATE0 = ptc_private%base_state
if (ptc_private%use_totalpath) STATE0 = STATE0 + TOTALPATH0

STATE = STATE0
if (bmad_com%spin_tracking_on) STATE = STATE0 + SPIN0

re = start_orb%vec

!-----------------------------
! track element

start2_orb = start_orb   ! Save initial state
end_orb = start_orb

call ele_to_fibre (ele, ptc_fibre, param, .true., err_flag, ref_in = start_orb)
if (err_flag) then
  end_orb%state = lost$
  return
endif

!

stm = ele%spin_tracking_method
if (bmad_com%spin_tracking_on .and. (stm == tracking$ .or. stm == symp_lie_ptc$) .or. present(track)) then
  ptc_probe = re
  ptc_probe%q%x = [1, 0, 0, 0]

  if (present(track)) then
    ptc_track => ptc_fibre%t1
    call save_this_step(track, ptc_probe, ele)

    do while (.not. associated(ptc_track, ptc_fibre%t2))
      call track_probe (ptc_probe, STATE, node1 = ptc_track, node2 = ptc_track%next)
      call save_this_step(track, ptc_probe, ele)
      if (.not. check_stable) exit
      ptc_track => ptc_track%next
    enddo

    call track_probe (ptc_probe, STATE, node1 = ptc_track, node2 = ptc_track%next)
    call save_this_step(track, ptc_probe, ele)

  else
    call track_probe (ptc_probe, STATE, fibre1 = ptc_fibre)
  endif

  end_orb%spin = quat_rotate(ptc_probe%q%x, start2_orb%spin)
  re = ptc_probe%x

else
  ! Orignally used track (ptc_fibre, re, STATE) but this will not track taylor elements correctly.
  call track_probe_x (re, STATE0, fibre1 = ptc_fibre)
endif

!-----------------------------

end_orb%vec = re

! 

if (ele%value(p0c$) /= ele%value(p0c_start$) .or. start2_orb%vec(6) /= end_orb%vec(6)) then
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
endif

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

if (ptc_private%use_totalpath) then
  end_orb%t = start2_orb%t + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
else
  end_orb%t = start2_orb%t + ele%value(delta_ref_time$) + &
                          start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

if (.not. check_stable) end_orb%state = lost$

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

! print '(i4, f10.6, 4x, 6f10.6, 4x, es16.8)', ptc_track%cas, ptc_track%s(1), ptc_probe%x, ptc_probe%E

orbit = start2_orb
orbit%vec = ptc_probe%x
orbit%s = ptc_track%s(1) + ele%s_start
orbit%spin = quat_rotate(ptc_probe%q%x, start2_orb%spin)
if (.not. check_stable) orbit%state = lost$
call save_a_step (track, ele, param, .false., orbit, ptc_track%s(1))

end subroutine save_this_step

end subroutine track1_symp_lie_ptc
