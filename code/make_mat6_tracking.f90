!+
! Subroutine make_mat6_tracking (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element using the
! Present tracking method.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!   bmad_com
!     %d_orb(6) -- Real(rp): Vector of offsets to use. 
!     %mat6_track_symmetric 
!               -- Logical: If True then track with +/- d_orb offsets so
!                  the tracking routine is called 12 times. 
!                  If False use only +d_orb offsets using only 6 tracks.
!                  Default is True.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %vec0  -- 0th order map component
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

subroutine make_mat6_tracking (ele, param, c0, c1)

use bmad_interface, except_dummy => make_mat6_tracking

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: c0, c00, c1, start, end1, end2
type (lat_param_struct)  param

real(rp) del_orb(6), dorb6, abs_p
integer i

! This computation is singular when start%vec(6) = -1 (zero starting velocity).
! In this case, shift start%vec(6) slightly to avoid the singularity.

del_orb = bmad_com%d_orb

if (bmad_com%mat6_track_symmetric) then
  abs_p = max(abs(c0%vec(2)) + abs(del_orb(2)), abs(c0%vec(4)) + abs(del_orb(4)), abs(del_orb(6)))
else
  abs_p = max(abs(c0%vec(2) + del_orb(2)), abs(c0%vec(4) + del_orb(4)), -del_orb(6)) 
endif

! The factor of 1.01 is used to avoid roundoff problems.

dorb6 = max(0.0_rp, 1.01 * (abs_p - (1 + c0%vec(6))))   ! Shift in start%vec(6) to apply.

call init_coord(c00, c0, ele, upstream_end$, c0%species, shift_vec6 = .false.)
call track1 (c00, ele, param, c1)

! Symmetric tracking uses more tracks but is more accurate

if (bmad_com%mat6_track_symmetric) then
  do i = 1, 6
    start = c00
    start%vec(6) = start%vec(6) + dorb6
    start%vec(i) = start%vec(i) + del_orb(i)
    call init_coord(start, start, ele, upstream_end$, c0%species, shift_vec6 = .false.)
    call track1 (start, ele, param, end2)

    start = c00
    start%vec(6) = start%vec(6) + dorb6
    start%vec(i) = start%vec(i) - del_orb(i)
    call init_coord(start, start, ele, upstream_end$, c0%species, shift_vec6 = .false.)
    call track1 (start, ele, param, end1)

    ele%mat6(1:6, i) = (end2%vec - end1%vec) / (2 * del_orb(i))
  enddo

! Else non-symmetric tracking only uses one offset for each phase space coordinate.

else  

  do i = 1, 6
    start = c00
    start%vec(6) = start%vec(6) + dorb6
    start%vec(i) = start%vec(i) + del_orb(i)
    call init_coord(start, start, ele, upstream_end$, c0%species, shift_vec6 = .false.)
    call track1 (start, ele, param, end1)
    ele%mat6(1:6, i) = (end1%vec - c1%vec) / del_orb(i)
  enddo
endif

! vec0 calc

ele%vec0 = c1%vec - matmul(ele%mat6, c0%vec)

end subroutine

