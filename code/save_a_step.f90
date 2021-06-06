!+
! Subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, mat6, make_matrix, rf_time)
!
! Routine used by the Runge-Kutta tracking routine to save the trajectory through an element.
!
! Notes: 
!   * It is assumed by this routine that here(:) is the orbit in local 
!     element coordinates. The actual track saved will be in laboratory coordinates.
!   * It is up to the calling routine to keep track of track%n_ok and track%n_bad if desired.
!     These numbers are used by adaptive trackers to record how many times the step size needed 
!     to be shortened.
!     
!
! Input:
!   track           -- track_struct: Track up to now. If track%n_pt < 0, the structure will be reinitialized.
!   ele             -- ele_struct: Element being tracked through.
!   param           -- lat_param_struct: Lattice parameters.
!   local_ref_frame -- Logical: If True then coordinates are wrt the frame of ref of the element.
!   orb             -- coord_struct: trajectory at s with respect to element coordinates.
!   s_rel           -- real(rp): Longitudinal position wrt the element. Lab coords if local_ref_frame = F
!                       and body coords if local_ref_fram = T
!   save_field      -- logical, optional: Save electric and magnetic field values? Default is False.
!   mat6(6,6)       -- real(rp), optional: Matrix to store.
!   make_matrix     -- logical, optional: Is mat6 a valid matrix? Default is False.
!   rf_time         -- real(rp), optional: RF clock time used for calculating the field.. 
!                       If not present then the time will be calculated using the standard algorithm.
!                       This is only needed if save_field = True.
!
! Ouput:
!   track           -- track_struct: Track with current trajectory info appended on.
!-

subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, mat6, make_matrix, rf_time)

use em_field_mod, dummy => save_a_step

implicit none

type (track_struct), target :: track, track2
type (ele_struct), target :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) orb
type (track_point_struct), pointer :: tp

integer n_pt, n, n_old
real(rp) s_rel
real(rp), optional :: mat6(6,6), rf_time
real(rp) s_lab
logical local_ref_frame
logical, optional :: save_field, make_matrix

! Init

if (.not. allocated (track%pt)) then
  allocate(track%pt(0:100))
  track%n_pt = -1
endif

if (track%n_pt < 0) then
  track%n_ok = 0
  track%n_bad = 0
  track%n_pt = -1
endif

!

track%n_pt = track%n_pt + 1
n_pt = track%n_pt
n_old = ubound(track%pt, 1)

if (n_pt > n_old) then
  n = 1.5 * n_pt
  call move_alloc (track%pt, track2%pt)
  allocate(track%pt(0:n))
  track%pt(:n_old) = track2%pt
end if

!

tp => track%pt(n_pt)
tp%orb = orb
if (logic_option(.false., make_matrix)) then
  tp%mat6 = mat6
else
  tp%mat6 = 0
endif

if (local_ref_frame) then
  call offset_particle (ele, param, unset$, tp%orb, drift_to_edge = .false., set_hvkicks = .false., &
                                                  s_pos = s_rel, s_out = s_lab, mat6 = tp%mat6, make_matrix = make_matrix)
else
  s_lab = s_rel
endif

tp%s_body = s_rel
tp%orb%ix_ele = ele%ix_ele

if (logic_option(.false., save_field)) then
  call em_field_calc (ele, param, s_lab, orb, local_ref_frame, tp%field, .false., rf_time = rf_time)
endif

track%n_ok = track%n_ok + 1

end subroutine save_a_step
