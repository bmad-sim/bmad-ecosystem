!+
! Subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb, track)
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
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   end_orb    -- Coord_struct: End position
!   track      -- Track_struct, optional: Structure holding the track information.
!-

subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb, track)

use ptc_interface_mod, except_dummy => track1_symp_lie_ptc
use ptc_spin, only: probe, assignment(=), operator(+), SPIN0, DEFAULT, track_probe, track_probe_x
use s_tracking, only: DEFAULT, alloc_fibre, integration_node
use mad_like, only: fibre, kill

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb, orbit
type (track_struct), optional :: track
type (ele_struct) :: ele, drift_ele
type (lat_param_struct) :: param
type (fibre), pointer :: fibre_ele
type (probe) ptc_probe
type (integration_node), pointer :: ptc_track

real(dp) re(6), beta0, beta1
integer i, stm

character(20) :: r_name = 'track1_symp_lie_ptc'

! call the PTC routines to track through the fibre.

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
beta1 = ele%value(p0c$) / ele%value(e_tot$)

call vec_bmad_to_ptc (start_orb%vec, beta0, re)

! Track a drift if using hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, upstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, fibre_ele, param, .true.)
  call track_probe_x (re, DEFAULT, fibre1 = fibre_ele)
endif  

!-----------------------------
! track element

start2_orb = start_orb
end_orb = start_orb
orbit = end_orb

call ele_to_fibre (ele, fibre_ele, param, .true., track_particle = start_orb)

stm = ele%spin_tracking_method
if (bmad_com%spin_tracking_on .and. (stm == tracking$ .or. stm == symp_lie_ptc$) .or. present(track)) then
  ptc_probe = re
  ptc_probe%s(1)%x = real(start_orb%spin, dp)

  if (present(track)) then
    ptc_track => fibre_ele%t1
    call save_this_step()

    do while (.not. associated(ptc_track, fibre_ele%t2))
      call track_probe (ptc_probe, DEFAULT+SPIN0, node1 = ptc_track, node2 = ptc_track%next)
      call save_this_step()
      ptc_track => ptc_track%next
    enddo

    call track_probe (ptc_probe, DEFAULT+SPIN0, node1 = ptc_track, node2 = ptc_track%next)
    call save_this_step()

  else
    call track_probe (ptc_probe, DEFAULT+SPIN0, fibre1 = fibre_ele)
  endif

  end_orb%spin = ptc_probe%s(1)%x
  re = ptc_probe%x

else
  ! Orignally used track (fibre_ele, re, DEFAULT) but this will not track taylor elements correctly.
  call track_probe_x (re, DEFAULT, fibre1 = fibre_ele)
endif

!-----------------------------

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, downstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, fibre_ele, param, .true.)
  call track_probe_x (re, DEFAULT, fibre1 = fibre_ele)
endif  

call vec_ptc_to_bmad (re, beta1, end_orb%vec)

! 

if (ele%value(p0c$) /= ele%value(p0c_start$) .or. start2_orb%vec(6) /= end_orb%vec(6)) then
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
endif

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

end_orb%t = start2_orb%t + ele%value(delta_ref_time$) + &
                          start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)


!---------------------------------------------------------------------
contains

subroutine save_this_step()

! The complication is that PTC pz is the true canonical momentum which includes the electrostatic potential.
! But Bmad pz does not include the electrostatic potential.

! Note: ptc_track%s(2) is the same as %s(1) except in a true rbend where %s(2) is the chord distance.

! print '(i4, f10.6, 4x, 6f10.6, 4x, es16.8)', ptc_track%cas, ptc_track%s(1), ptc_probe%x, ptc_probe%E

re = ptc_probe%x
re(5) = 1e9 * ptc_probe%E / end_orb%p0c       ! ptc_probe%E = Delta E in Gev
call vec_ptc_to_bmad (re, beta1, orbit%vec)
orbit%s = ptc_track%s(1) + ele%s_start
orbit%spin = ptc_probe%s(1)%x
call save_a_step (track, ele, param, .false., orbit)

end subroutine save_this_step

end subroutine track1_symp_lie_ptc
