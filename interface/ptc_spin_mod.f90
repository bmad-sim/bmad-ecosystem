module ptc_spin_mod

use ptc_layout_mod

implicit none

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine transfer_map_calc_with_spin (lat, t_map, s_map, orb0, err_flag, &
!																							ix1, ix2, ix_branch, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between two elements.
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
! Note: If integrate = False and if a taylor map does not exist for an
! element this routine will make one and store it in the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- lat_struct: Lattice used in the calculation.
!   t_map(6)   -- taylor_struct: Initial orbital map (used when unit_start = False)
!   s_map(3,3) -- taylor_struct: Initial spin map (used when unit_start = False)
!		orb0			 -- coord_struct: Initial orbit around which the map is made.
!   ix1        -- integer, optional: Element start index for the calculation.
!                   Default is 0.
!   ix2        -- integer, optional: Element end index for the calculation.
!                   Default is lat%n_ele_track.
!   ix_branch  -- integer, optional: Lattice branch index. Default is 0.
!   one_turn   -- logical, optional: If present and True, and if ix1 = ix2,
!                   and the lattice is circular, then construct the one-turn
!                   map from ix1 back to ix1. Default = False.
!   unit_start -- logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!   t_map(6)   -- Taylor_struct: Orbital transfer map.
!   s_map(3,3) -- Taylor_struct: Spin transfer map.
!   err_flag   -- logical: Set True if problem like number overflow, etc.
!-

subroutine transfer_map_calc_with_spin (lat, t_map, s_map, orb0, err_flag, &
																								ix1, ix2, ix_branch, one_turn, unit_start)

use pointer_lattice, only: assignment(=), operator(+), operator(**), operator(-), operator(/), operator(*), &
                      kill, internal_state, probe, probe_8, sqrt, c_damap, real_8, default, spin0, alloc, propagate

type (lat_struct), target :: lat
type (taylor_struct) :: t_map(6), s_map(3,3)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) orb0
type (internal_state) ptc_state
type (probe) ptc_probe
type (probe_8) ptc_probe8
type (c_damap) ptc_c_map
type (real_8) y0(6), y2(6), bet

real(rp) x(6), beta

integer, optional :: ix1, ix2, ix_branch
integer i, i1, i2

logical, optional :: one_turn, unit_start
logical err_flag


!

call set_ptc(init_complex = .true.)

branch => lat%branch(integer_option(0, ix_branch))
if (.not. associated (lat%branch(0)%ptc%m_t_layout)) call lat_to_ptc_layout (lat)

i1 = integer_option(0, ix1) 
i2 = integer_option(branch%n_ele_track, ix2)

!

call alloc (ptc_c_map)
call alloc (bet)

t_map(:)%ref = orb0%vec
x = orb0%vec

ptc_state = DEFAULT + SPIN0
ptc_c_map = 1
ptc_probe = orb0%vec
ptc_probe8 = ptc_probe + ptc_c_map

! Bmad to PTC

ele => branch%ele(i1)
beta = ele%value(p0c$) / ele%value(e_tot$)
call real_8_init(y0)
y0 = x ! = IdentityMap + const
ptc_probe8%x = y0
ptc_probe8%x(5) = (y0(6)**2+2.d0*y0(6))/(1.d0/beta+sqrt( 1.d0/beta**2+y0(6)**2+2.d0*y0(6)) )
bet = (1.d0+y0(6))/(1.d0/beta+ptc_probe8%x(5))
ptc_probe8%x(6) = -y0(5)/bet

! Track

call propagate (branch%ptc%m_t_layout, ptc_probe8, +ptc_state, &
														branch%ele(i1)%ptc_fibre%pos+1, branch%ele(i2)%ptc_fibre%pos+1)

! PTC to Bmad

ele => branch%ele(i2)
beta = ele%value(p0c$) / ele%value(e_tot$)
y0 = ptc_probe8%x
y0(6) = (2.d0*ptc_probe8%x(5)/beta+ptc_probe8%x(5)**2)/(sqrt(1.d0+2.d0*ptc_probe8%x(5)/beta+ptc_probe8%x(5)**2)+1.d0)
bet = (1.d0+y0(6))/(1.d0/beta+ptc_probe8%x(5))
y0(5) = -bet*ptc_probe8%x(6)

! take out the offset

if (any(x /= 0)) then
  call real_8_init(y2)
  y2 = -x  ! y2 = IdentityMap - x

  call concat_real_8 (y2, y0, y0)

!  do i = 1, 3
!    call concat_real_8 (y2, ptc_probe8%s(i)%x, ptc_probe8%s(i)%x)
!  enddo

  call kill(y2)
endif

!

t_map = y0

do i = 1, 3
  s_map(:,i) = ptc_probe8%s(i)%x
enddo

call kill (ptc_c_map)
call kill (y0)
call kill (bet)

end subroutine transfer_map_calc_with_spin

end module
