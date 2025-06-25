!+
! Subroutine make_mat6_tracking (ele, param, start_orb, end_orb, err_flag, spin_only)
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
!   spin_only -- logical, optional: Default False. If True, only calculate ele%spin_taylor. 
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0         -- 0th order map component
!     %mat6         -- 6x6 transfer matrix.
!     %spin_taylor  -- Spin taylor map (if spin tracking is on).
!   end_orb   -- Coord_struct: Coordinates at the end of element.
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_tracking (ele, param, start_orb, end_orb, err_flag, spin_only)

use bmad_interface, except_dummy => make_mat6_tracking

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, start_orb0, end_orb, start, end1, end2
type (lat_param_struct)  param

real(rp) del_orb(6), dorb6, abs_p, mat6(6,6)
real(rp) q1(0:3), q2(0:3), q_map(0:3,0:6)
integer i
logical err_flag, do_spin
logical, optional :: spin_only

character(*), parameter :: r_name = 'make_mat6_tracking'

! This computation is singular when start%vec(6) = -1 (zero starting velocity).
! In this case, shift start%vec(6) slightly to avoid the singularity.

err_flag = .true.
del_orb = bmad_com%d_orb
abs_p = max(abs(start_orb%vec(2)) + abs(del_orb(2)), abs(start_orb%vec(4)) + abs(del_orb(4)), abs(del_orb(6)))
bmad_private%random_on = .false.
do_spin = (logic_option(.false., spin_only) .or. bmad_com%spin_tracking_on)

! The factor of 1.01 is used to avoid roundoff problems.
! Note: init_coord is avoided since init_coord will make z and t consistent with the element's t_ref.
! However, the reference time of the particle is not necessarily the same as the element's ref time.

dorb6 = max(0.0_rp, 1.01 * (abs_p - (1 + start_orb%vec(6))))   ! Shift in start%vec(6) to apply.
start_orb0 = start_orb

call track_this (start_orb0, ele, param, end_orb, q_map(:,0), .false.)
if (end_orb%state /= alive$) then
  call out_io (s_error$, r_name, 'CENTRAL PARTICLE LOST IN TRACKING TO CONSTRUCT TRANSFER MATRIX FOR ELEMENT: ' // ele%name)
  bmad_private%random_on = .true.
  return
endif

! Tracking

do i = 1, 6
  start = start_orb0
  start%vec(6) = start%vec(6) + dorb6
  start%vec(i) = start%vec(i) + del_orb(i)
  call track_this (start, ele, param, end2, q1, .true.)
  if (end2%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING TO CONSTRUCT TRANSFER MATRIX FOR ELEMENT: ' // ele%name, &
                                   'WITH POSITIVE OFFSET IN PHASE SPACE COORDINATE: ' // int_str(i))
    bmad_private%random_on = .true.
    return
  endif

  start = start_orb0
  start%vec(6) = start%vec(6) + dorb6
  start%vec(i) = start%vec(i) - del_orb(i)
  call track_this (start, ele, param, end1, q2, .true.)
  if (end1%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING TO CONSTRUCT TRANSFER MATRIX FOR ELEMENT: ' // ele%name, &
                                   'WITH NEGATIVE OFFSET IN PHASE SPACE COORDINATE: ' // int_str(i))
    bmad_private%random_on = .true.
    return
  endif

  mat6(1:6, i) = (end2%vec - end1%vec) / (2 * del_orb(i))
  if (do_spin) q_map(:, i) = (q2 - q1) / (2 * del_orb(i))
enddo

! vestart_orb calc

if (.not. logic_option(.false., spin_only)) then
  ele%vec0 = end_orb%vec - matmul(mat6, start_orb%vec)
  ele%mat6 = mat6
endif

if (do_spin) then
  ele%spin_q = q_map
  call linear_to_spin_taylor(q_map, ele%spin_taylor)
endif

bmad_private%random_on = .true.
err_flag = .false.

!------------------------------------------------------
contains

subroutine track_this(starting, ele, param, ending, q, adjust)

type (coord_struct) starting, ending
type (ele_struct) ele
type (lat_param_struct) param

real(rp) q(0:3), m(3,3)
integer i
logical adjust

!

if (adjust) then
  if (start_orb%species == photon$) then
   call init_coord(starting, starting, ele, start_end$, start_orb%species)
   return
  endif

  !

  call convert_pc_to (starting%p0c * (1 + starting%vec(6)), starting%species, beta = starting%beta)

  if (start_orb%beta == 0) then
    starting%t = start_orb%t - starting%vec(5) / (c_light * starting%beta)
  else
    starting%t = start_orb%t - starting%vec(5) / (c_light * starting%beta) + start_orb%vec(5) / (c_light * start_orb%beta)
  endif
endif

if (do_spin) then
  do i = 1, 3
    starting%spin = 0
    starting%spin(i) = 1
    call track1(starting, ele, param, ending)
    m(:,i) = ending%spin
  enddo
  q = w_mat_to_quat(m)
  ele%spin_taylor_ref_orb_in = start_orb%vec
else
  call track1(starting, ele, param, ending)
endif

end subroutine track_this

end subroutine

