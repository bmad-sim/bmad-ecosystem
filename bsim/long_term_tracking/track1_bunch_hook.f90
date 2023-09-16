!+
! Subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
!
! Routine that can be customized for tracking a bunch through a single element.
!
! Input:
!   bunch_start   -- Bunch_struct: Starting bunch position.
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

subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

use lt_tracking_mod, dummy => track1_bunch_hook

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (bunch_track_struct), optional :: bunch_track
type (coord_struct), optional :: centroid(0:)
type (coord_struct), pointer :: orb

real(rp) t, r

integer, optional :: direction
integer ip, ir, ie, n, iv

logical err, finished

! Rampers are only applied to the element once per bunch. That is, it is assumed 
! that the ramper control function variation is negligible over the time scale of a bunch passage. 
! To evaluate multiple times in a bunch passage would, in general, be wrong if using ran() or ran_gauss().

err = .false.
finished = .false.
if (.not. ltt_params_global%ramping_on) return 
if (ltt_params_global%ramp_update_each_particle) return 

n = ltt_com_global%n_ramper_loc
if (n == 0) return

t = sum(bunch%particle%t, bunch%particle%state == alive$) / &
           count(bunch%particle%state == alive$) + 0.5_rp * ele%value(delta_ref_time$) + &
           ltt_params_global%ramping_start_time

do ir = 1, ltt_com_global%n_ramper_loc
  do iv = 1, size(ltt_com_global%ramper(ir)%ele%control%var)
    if (ltt_com_global%ramper(ir)%ele%control%var(iv)%name /= 'TIME') cycle
    ltt_com_global%ramper(ir)%ele%control%var(iv)%value = t
  enddo
enddo

n = ltt_com_global%n_ramper_loc
call ltt_apply_rampers_to_slave (ele, ltt_com_global%ramper(1:n), err)

! The beginning element is never tracked through. If there is energy ramping and the user is writing out 
! p0c or E_tot from the beginning element, the user may be confused since these values will not change. 
! So adjust the beginning element's p0c and E_tot to keep users happy.

if (ele%ix_ele == 1) then
  ele0 => pointer_to_next_ele(ele, -1)
  ele0%value(p0c$) = ele%value(p0c_start$)
  ele0%value(E_tot$) = ele%value(E_tot_start$)
endif

! Adjust particle reference energy if needed.

if (bunch%particle(1)%p0c == ele%value(p0c_start$)) return

do ip = 1, size(bunch%particle)
  orb => bunch%particle(ip)
  if (orb%state /= alive$) cycle
  r = orb%p0c / ele%value(p0c_start$)
  orb%vec(2) = r * orb%vec(2)
  orb%vec(4) = r * orb%vec(4)
  orb%vec(6) = r * orb%vec(6) + (orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
  orb%p0c = ele%value(p0c_start$)
enddo

end subroutine
