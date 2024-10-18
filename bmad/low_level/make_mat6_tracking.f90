!+
! Subroutine make_mat6_tracking (ele, param, start_orb, end_orb, err_flag)
!
! Subroutine to make the 6x6 transfer matrix for an element using the
! Present tracking method.
!
! bmad_com common block settings:
!   bmad_com
!     %d_orb(6)  -- Real(rp): Vector of offsets to use. 
!
! Input:
!   ele       -- Ele_struct: Element with transfer matrix
!   param     -- lat_param_struct: Parameters are needed for some elements.
!   start_orb -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0      -- 0th order map component
!     %mat6      -- 6x6 transfer matrix.
!   end_orb   -- Coord_struct: Coordinates at the end of element.
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_tracking (ele, param, start_orb, end_orb, err_flag)

use bmad_interface, except_dummy => make_mat6_tracking

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, start_orb0, end_orb, start, end1, end2
type (lat_param_struct)  param

real(rp) del_orb(6), dorb6, abs_p, mat6(6,6)
integer i
logical err_flag
character(*), parameter :: r_name = 'make_mat6_tracking'

! This computation is singular when start%vec(6) = -1 (zero starting velocity).
! In this case, shift start%vec(6) slightly to avoid the singularity.

err_flag = .true.
del_orb = bmad_com%d_orb
abs_p = max(abs(start_orb%vec(2)) + abs(del_orb(2)), abs(start_orb%vec(4)) + abs(del_orb(4)), abs(del_orb(6)))
bmad_private%random_on = .false.

! The factor of 1.01 is used to avoid roundoff problems.
! Note: init_coord is avoided since init_coord will make z and t consistent with the element's t_ref.
! However, the reference time of the particle is not necessarily the same as the element's ref time.

dorb6 = max(0.0_rp, 1.01 * (abs_p - (1 + start_orb%vec(6))))   ! Shift in start%vec(6) to apply.
start_orb0 = start_orb

call track1 (start_orb0, ele, param, end_orb)
if (end_orb%state /= alive$) then
  call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING CENTRAL PARTICLE. MATRIX NOT CALCULATED FOR ELEMENT: ' // ele%name)
  bmad_private%random_on = .true.
  return
endif

! Tracking

do i = 1, 6
  start = start_orb0
  start%vec(6) = start%vec(6) + dorb6
  start%vec(i) = start%vec(i) + del_orb(i)
  call adjust_this
  call track1 (start, ele, param, end2)
  if (end2%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING (+). MATRIX NOT CALCULATED FOR ELEMENT: ' // ele%name)
    bmad_private%random_on = .true.
    return
  endif

  start = start_orb0
  start%vec(6) = start%vec(6) + dorb6
  start%vec(i) = start%vec(i) - del_orb(i)
  call adjust_this
  call track1 (start, ele, param, end1)
  if (end1%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING (-). MATRIX NOT CALCULATED FOR ELEMENT: ' // ele%name)
    bmad_private%random_on = .true.
    return
  endif

  mat6(1:6, i) = (end2%vec - end1%vec) / (2 * del_orb(i))
enddo

! vestart_orb calc

ele%vec0 = end_orb%vec - matmul(mat6, start_orb%vec)
ele%mat6 = mat6

bmad_private%random_on = .true.
err_flag = .false.

!------------------------------------------------------
contains

subroutine adjust_this

if (start_orb%species == photon$) then
 call init_coord(start, start, ele, start_end$, start_orb%species)
 return
endif

!

call convert_pc_to (start%p0c * (1 + start%vec(6)), start%species, beta = start%beta)

if (start_orb%beta == 0) then
  start%t = start_orb%t - start%vec(5) / (c_light * start%beta)
else
  start%t = start_orb%t - start%vec(5) / (c_light * start%beta) + start_orb%vec(5) / (c_light * start_orb%beta)
endif

end subroutine adjust_this

end subroutine

