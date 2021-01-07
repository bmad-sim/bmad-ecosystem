!+
! Subroutine track1_bunch_hook (bunch_start, lat, ele, bunch_end, err, centroid, direction, finished)
!
! Routine that can be customized for tracking a bunch through a single element.
!
! Input:
!   bunch_start   -- Bunch_struct: Starting bunch position.
!   lat          -- lat_struct: Lattice containing element to be tracked through.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   bunch_end    -- bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_bunch code will not be called.
!-

subroutine track1_bunch_hook (bunch_start, lat, ele, bunch_end, err, centroid, direction, finished)

use lt_tracking_mod, dummy => track1_bunch_hook

implicit none

type (bunch_struct), target :: bunch_start
type (bunch_struct) :: bunch_end
type (lat_struct) :: lat
type (ele_struct) ele
type (coord_struct), optional :: centroid(0:)
type (coord_struct), pointer :: orb

real(rp) t, r

integer, optional :: direction
integer ip, ir, ie

logical err, finished

! To save time, rampers are only applied to the element once per bunch. 
! That is, it is assumed that the ramper control function variation is negligible over the
! time scale of a bunch passage.

finished = .false.
if (size(ltt_com_global%ix_ramper) == 0) return

t = sum(bunch_start%particle%t, bunch_start%particle%state == alive$) / &
                     count(bunch_start%particle%state == alive$) + 0.5_rp * ele%value(delta_ref_time$)

do ir = 1, size(ltt_com_global%ix_ramper)
  ie = ltt_com_global%ix_ramper(ir)
  if (ie == -1) exit
  ltt_com_global%lat%ele(ie)%control%var(1)%value = t
  call apply_ramper (ele, ltt_com_global%lat%ele(ie), err)
enddo

! Adjust particle reference energy if needed.

if (bunch_start%particle(1)%p0c == ele%value(p0c_start$)) return

do ip = 1, size(bunch_start%particle)
  orb => bunch_start%particle(ip)
  if (orb%state /= alive$) cycle
  r = orb%p0c / ele%value(p0c_start$)
  orb%vec(2) = r * orb%vec(2)
  orb%vec(4) = r * orb%vec(4)
  orb%vec(6) = r * orb%vec(6) + (orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
  orb%p0c = ele%value(p0c_start$)
enddo

end subroutine
