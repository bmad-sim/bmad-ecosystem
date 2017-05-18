!+
! Subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, rf_time)
!
! Routine used by Runge-Kutta and Boris tracking routines to save
! the trajectory through an element.
!
! Note: It is assumed by this routine that here(:) is the orbit in local 
! element coordinates. The actual track saved will be in laboratory coordinates.
!
! Input:
!   ele        -- ele_struct: Element being tracked through.
!   param      -- lat_param_struct: Lattice parameters.
!   local_ref_frame -- Logical: If True then coordinates are wrt the frame of ref of the element.
!   orb        -- Coord_struct: trajectory at s with respect to element coordinates.
!   s_rel      -- Real(rp), optional: Position with respect to start of element. Only used for calculating the field.
!   save_field -- logical, optional: Save electric and magnetic field values? Default is False.
!   rf_time    -- real(rp), optional: RF clock time used for calculating the field.. 
!                   If not present then the time will be calculated using the standard algorithm.
!                   This is only needed if save_field = True.
!
! Ouput:
!   track    -- track_struct: Trajectory structure to save to.
!   s_sav    -- Real(rp): Set equal to s.
!-

subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, rf_time)

use em_field_mod, dummy => save_a_step

implicit none

type (track_struct) track, track2
type (ele_struct), target :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) orb, orb2
integer n_pt, n, n_old
real(rp), optional :: rf_time, s_rel
logical local_ref_frame
logical, optional :: save_field

! Not allocated

if (.not. allocated (track%orb)) then
  allocate(track%orb(0:100))
  allocate(track%field(0:100))
  allocate(track%map(0:100))
  track%n_ok = 0
  track%n_bad = 0
  track%n_pt = -1
endif

!

track%n_pt = track%n_pt + 1
n_pt = track%n_pt
n_old = ubound(track%orb, 1)

if (n_pt > n_old) then
  n = 1.5 * n_pt
  call move_alloc (track%orb, track2%orb)
  call move_alloc (track%field, track2%field)
  call move_alloc (track%map, track2%map)
  allocate(track%orb(0:n), track%field(0:n), track%map(0:n))
  track%orb(:n_old) = track2%orb; track%field(:n_old) = track2%field; track%map(:n_old) = track2%map
  deallocate(track2%orb, track2%field, track2%map)
end if

! Notice that a translation due to a finite ele%value(z_offset$) is not wanted here.

orb2 = orb
if (local_ref_frame) call offset_particle (ele, param, unset$, orb2, &
          set_z_offset = .false., set_multipoles = .false., set_hvkicks = .false., ds_pos = s_rel)

track%orb(n_pt) = orb2
track%orb(n_pt)%ix_ele = ele%ix_ele
track%map(n_pt)%mat6 = 0

if (logic_option(.false., save_field)) then
  call em_field_calc (ele, param, s_rel, orb, local_ref_frame, track%field(n_pt), .false., rf_time = rf_time)
endif

end subroutine save_a_step
