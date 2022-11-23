!+         
! Subroutine transfer_map_calc (lat, t_map, err_flag, ix1, ix2, ref_orb, ix_branch, one_turn, 
!                                                                          unit_start, concat_if_possible)
!
! Subroutine to calculate the transfer map between two elements.
!
! To calculate just the first order transfer matrices see the routine: 
!   transfer_matrix_calc
!
! The transfer map is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn map is calculated.
!
! If ix2 < ix1 and lat%param%geometry is closed$ then the
! calculation will "wrap around" the lattice end.
! For example, if ix1 = 900 and ix2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If ix2 < ix1 and lat%param%geometry is open$ then the backwards
! transfer matrix is computed.
!
! If ix2 = ix1 then you get the unit map except if one_turn = True and the lattice is circular.
!
! The default calculation uses tracking to propagate the map through the lattice elements.
! This, while accurate, can be slow. If the concat_if_possible argument is set to True, and if a
! lattice element has a Taylor map, map concatenation rather than tracking will be used. Care must
! be taken with concatination since if the reference orbits about which the two maps to be concatenated
! are different, the resulting concatenated map may be off.
!
! Input:
!   lat        -- lat_struct: Lattice used in the calculation.
!   t_map(6)   -- taylor_struct: Initial map (used when unit_start = False)
!   ix1        -- integer, optional: Element start index for the calculation.
!                   Default is 0.
!   ix2        -- integer, optional: Element end index for the calculation.
!                   Default is lat%n_ele_track.
!   ref_orb    -- coord_struct, optional: Reference orbit/particle at s1 around which the map is made.
!                   This arg is needed if: unit_start = True or particle is not the same as the reference 
!                   particle of the lattice.
!   ix_branch  -- integer, optional: Lattice branch index. Default is 0.
!   one_turn   -- logical, optional: If present and True, and if ix1 = ix2,
!                   and the lattice is circular, then construct the one-turn 
!                   map from ix1 back to ix1. Default = False.
!   unit_start -- logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!   concat_if_possible
!              -- logical, optional: If present and True then use map concatenation rather than tracking 
!                   if a map is present for a given lattice element. See above. Default is False.
!
! Output:
!   t_map(6)   -- Taylor_struct: Transfer map.
!   err_flag   -- logical: Set True if problem like number overflow, etc.
!-

subroutine transfer_map_calc (lat, t_map, err_flag, ix1, ix2, ref_orb, ix_branch, one_turn, &
                                                                        unit_start, concat_if_possible)

use bmad_interface, except_dummy => transfer_map_calc
use ptc_interface_mod, only: concat_ele_taylor, taylor_propagate1, taylor_inverse, concat_taylor

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (taylor_struct) :: t_map(:), a_map(6)
type (coord_struct), optional :: ref_orb
type (coord_struct) orb0

real(rp) :: ex_factor(0:ptc_private%taylor_order_ptc)

integer, intent(in), optional :: ix1, ix2, ix_branch
integer i, i1, i2

logical, optional :: one_turn, unit_start, concat_if_possible
logical unit_start_this, err_flag

character(*), parameter :: r_name = 'transfer_map_calc'

!

unit_start_this = logic_option(.true., unit_start)

branch => lat%branch(integer_option(0, ix_branch))

i1 = integer_option(0, ix1) 
i2 = integer_option(branch%n_ele_track, ix2)

err_flag = .false.

if (unit_start_this) then
  call taylor_make_unit (t_map, ref_orb%vec)
endif

if (i1 == i2 .and. (lat%param%geometry == open$ .or. .not. logic_option (.false., one_turn))) return

if (present(ref_orb)) then
  call init_coord (orb0, ref_orb, branch%ele(i1), downstream_end$)
else
  orb0 = coord_struct()
  orb0%species = branch%param%particle
endif

! Map term overflow defined by |term| > 10^(20*(n+1)) where n = sum(term_exponents)

ex_factor(0) = 1d20
do i = 1, ubound(ex_factor, 1)
  ex_factor(i) = ex_factor(i-1) * 1d20
enddo

! Normal: i1 < i2.

if (i1 < i2) then 
  do i = i1+1, i2
    call add_on_to_map (t_map, branch%ele(i))
    if (err_flag) return
  enddo

! Circular lattice with i1 > i2: Track through origin.

elseif (branch%param%geometry == closed$) then
  do i = i1+1, branch%n_ele_track
    call add_on_to_map (t_map, branch%ele(i))
    if (err_flag) return
  enddo
  do i = 1, i2
    call add_on_to_map (t_map, branch%ele(i))
    if (err_flag) return
  enddo

! Open lattice with i1 > i2: Track backwards.

else

  if (unit_start_this) then
    do i = i2+1, i1
      call add_on_to_map (t_map, branch%ele(i))
      if (err_flag) return
    enddo
    call taylor_inverse (t_map, t_map)

  else
    call taylor_make_unit (a_map)
    do i = i2+1, i1
      call add_on_to_map (a_map, branch%ele(i))
      if (err_flag) return
    enddo
    call taylor_inverse (a_map, a_map)
    call concat_taylor (t_map, a_map, t_map)
    call kill_taylor (a_map)
  endif

endif

!--------------------------------------------------------
contains

subroutine add_on_to_map (map, ele)

type (taylor_struct) :: map(:)
type (ele_struct) ele
integer i, k

!

if (logic_option(.false., concat_if_possible) .and. associated(ele%taylor(1)%term)) then
  call concat_ele_taylor (map, ele, map)
else
  call taylor_propagate1 (map, ele, branch%param, orb0)
endif

! Check for overflow
! Map term overflow defined by |term| > 10^(20*n) where n = sum(term_exponents)

do i = 1, 6
  do k = 1, size(map(i)%term)
    if (abs(map(i)%term(k)%coef) < ex_factor(sum(map(i)%term(k)%expn))) cycle
    err_flag = .true.
    !! call out_io (s_error$, r_name, 'TAYLOR MAP TERM FLOATING OVERFLOW.')
    return
  enddo
enddo

end subroutine add_on_to_map

end subroutine transfer_map_calc
