!+
! Subroutine track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)
!
! Dummy routine for pre-processing at the start of the track1 routine.
!
! Also see:
!   track1_postprocess
!   track1_custom
!
! The radiation_included argument should be set to True if this routine (or a modified version of track1_custom)
! takes into account radiation damping and/or excitation. This will prevent track1 from calling track1_radiation.
! Note: If symp_lie_bmad is being called by this routine, symp_lie_bmad does take into account radiation effects.
! 
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Input:
!   start_orb  -- coord_struct: Starting position at the beginning of ele.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   start_orb   -- coord_struct: Modified starting position for track1 to use.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt further processing and use the modified 
!                     start_orb as the final end position of the particle.
!   radiation_included
!               -- logical: Should be set True if radiation damping/excitation is included in the tracking.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)

use da_program_mod, except_dummy => track1_preprocess

implicit none

type (coord_struct) :: start_orb
type (ele_struct) :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param
type (track_struct), optional :: track
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) r, t
integer ir, n, iu
logical err_flag, finished, radiation_included, is_there

character(*), parameter :: r_name = 'track1_preprocess'

!

err_flag = .false.
if (.not. da_com%ramping_on) return 
t = start_orb%t + 0.5_rp * ele%value(delta_ref_time$) + da_com%ramping_start_time

n = da_com%n_ramper_loc
allocate(eles(n))
do ir = 1, n
  ele0 => pointer_to_ele(ele%branch%lat, da_com%ramper(ir))
  eles(ir)%ele => ele0
  if (ele0%control%var(1)%name /= 'TIME') cycle
  ele0%control%var(1)%value = t
enddo

call apply_ramper (ele, eles, err_flag)

! The beginning element is never tracked through. If there is energy ramping and the user is writing out 
! p0c or E_tot from the beginning element, the user may be confused since these values will not change. 
! So adjust the beginning element's p0c and E_tot to keep users happy.

if (ele%ix_ele == 1) then
  ele0 => pointer_to_next_ele(ele, -1)
  ele0%value(p0c$) = ele%value(p0c_start$)
  ele0%value(E_tot$) = ele%value(E_tot_start$)
endif

! Adjust particle reference energy if needed.

if (start_orb%p0c == ele%value(p0c_start$)) return

r = start_orb%p0c / ele%value(p0c_start$)
start_orb%vec(2) = r * start_orb%vec(2)
start_orb%vec(4) = r * start_orb%vec(4)
start_orb%vec(6) = r * start_orb%vec(6) + (start_orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
start_orb%p0c = ele%value(p0c_start$)

end subroutine
