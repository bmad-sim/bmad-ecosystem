!+
! Subroutine ptc_transfer_map_with_spin (branch, t_map, s_map, orb0, err_flag, ix1, ix2, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between two elements.
! To calculate just the first order transfer matrices see the routine:
!   transfer_matrix_calc
!
! The transfer map is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn map is calculated.
!
! If ix2 < ix1 and branch%param%geometry is closed$ then the
! calculation will "wrap around" the lattice end.
! For example, if ix1 = 900 and ix2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If ix2 < ix1 and branch%param%geometry is open$ then the backwards
! transfer matrix is computed.
!
! If ix2 = ix1 then you get the unit map except if one_turn = True and the lattice is circular.
!
! Note: If integrate = False and if a taylor map does not exist for an
! element this routine will make one and store it in the element.
!
! Input:
!   branch     -- branch_struct: Lattice branch used in the calculation.
!   t_map(6)   -- taylor_struct: Initial orbital map (used when unit_start = False)
!   s_map(3,3) -- taylor_struct: Initial spin map (used when unit_start = False)
!   orb0       -- coord_struct: Initial orbit around which the map is made.
!   ix1        -- integer, optional: Element start index for the calculation.
!                   Default is 0.
!   ix2        -- integer, optional: Element end index for the calculation.
!                   Default is branch%n_ele_track.
!   one_turn   -- logical, optional: If present and True, and if ix1 = ix2,
!                   and the lattice is circular, then construct the one-turn
!                   map from ix1 back to ix1. Default = False.
!   unit_start -- logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!   t_map(6)   -- Taylor_struct: Orbital transfer map.
!   s_map(4)   -- Taylor_struct: Quaternion spin transfer map.
!   err_flag   -- logical: Set True if problem like number overflow, etc.
!-

subroutine ptc_transfer_map_with_spin (branch, t_map, s_map, orb0, err_flag, ix1, ix2, one_turn, unit_start)

use ptc_layout_mod, dummy => ptc_transfer_map_with_spin 
use pointer_lattice

implicit none

type (branch_struct) :: branch
type (taylor_struct) :: t_map(6), s_map(4)
type (ele_struct), pointer :: ele
type (coord_struct) orb0
type (internal_state) ptc_state
type (probe) ptc_probe
type (probe_8) ptc_probe8
type (c_damap) ptc_c_map
type (real_8) y0(6), y2(6)
type (fibre), pointer :: fib1, fib2

real(rp) x(6)

integer, optional :: ix1, ix2
integer i, i1, i2

logical err_flag
logical, optional :: one_turn, unit_start

!

if (.not. associated (branch%ptc%m_t_layout)) call lat_to_ptc_layout (branch%lat)

i1 = integer_option(0, ix1) 
i2 = integer_option(branch%n_ele_track, ix2)

!

call alloc (ptc_c_map)
call alloc (ptc_probe8)
call alloc (y0)

t_map(:)%ref = orb0%vec
x = orb0%vec

ptc_state = ptc_private%base_state + SPIN0
ptc_c_map = 1
ptc_probe = orb0%vec
ptc_probe8 = ptc_probe + ptc_c_map

! Track

fib1 => pointer_to_fibre(branch%ele(i1))
fib2 => pointer_to_fibre(branch%ele(i2))

if (logic_option(.false., one_turn) .and. i1 == i2 .and. branch%param%geometry == closed$) then
  call propagate (branch%ptc%m_t_layout, ptc_probe8, +ptc_state, fib1%pos)
else
  call propagate (branch%ptc%m_t_layout, ptc_probe8, +ptc_state, fib1%pos, fib2%pos)
endif

! take out the offset

y0 = ptc_probe8%x

if (any(x /= 0)) then
  call alloc(y2)
  y2 = -x  ! y2 = IdentityMap - x
  call concat_real_8 (y2, y0, y0, keep_y1_const_terms = .true.)
  call kill(y2)
endif

t_map = y0

do i = 1, 4
  s_map(i) = ptc_probe8%q%x(i-1)%t
enddo

call kill (ptc_c_map)
call kill (y0)
call kill (ptc_probe8)

end subroutine ptc_transfer_map_with_spin
