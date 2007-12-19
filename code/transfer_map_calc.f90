!+         
! Subroutine transfer_map_calc (lat, t_map, ix1, ix2, &
!                                         integrate, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between two elements.
! To calculate just the first order transfer matrices see the routine: 
!   transfer_matrix_calc
!
! The transfer map is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn map is calculated.
!
! If ix2 < ix1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will "wrap around" the lattice end.
! For example, if ix1 = 900 and ix2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If ix2 < ix1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
!
! If ix2 = ix1 then you get the unit map except if one_turn = True.
!
! Note: If integrate = False and if a taylor map does not exist for an 
! element this routine will make one and store it in the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Lat_struct: Lattice used in the calculation.
!   t_map(6)   -- Taylor_struct: Initial map (used when unit_start = False)
!   ix1        -- Integer, optional: Element start index for the calculation.
!                   Default is 0.
!   ix2        -- Integer, optional: Element end index for the calculation.
!                   Default is lat%n_ele_track.
!   integrate  -- Logical, optional: If present and True then do symplectic
!                   integration instead of concatenation. 
!                   Default = False.
!   one_turn   -- Logical, optional: If present and True, and if ix1 = ix2,
!                   then construct the one-turn map from ix1 back to ix1.
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

subroutine transfer_map_calc (lat, t_map, ix1, ix2, &
                                           integrate, one_turn, unit_start)

  use bmad_struct
  use bmad_interface, except_dummy => transfer_map_calc
  use ptc_interface_mod, only: concat_ele_taylor, ele_to_taylor, &
                               taylor_propagate1, taylor_inverse

  implicit none

  type (lat_struct) lat

  type (taylor_struct) :: t_map(:)

  integer, intent(in), optional :: ix1, ix2
  integer i, i1, i2

  logical, optional :: integrate, one_turn, unit_start
  logical integrate_this, one_turn_this, unit_start_this

  character(20) :: r_name = "transfer_map_calc"

!

  integrate_this  = logic_option (.false., integrate)
  one_turn_this   = logic_option (.false., one_turn)
  unit_start_this = logic_option(.true., unit_start)

  i1 = integer_option(0, ix1) 
  i2 = integer_option(lat%n_ele_track, ix2)
 
  if (unit_start_this) call taylor_make_unit (t_map)

  if (i1 == i2 .and. .not. one_turn_this) return

! Normal: i1 < i2.

  if (i1 < i2) then 
    do i = i1+1, i2
      call add_on_to_t_map
    enddo

! Circular lattice with i1 > i2: Track through origin.

  elseif (lat%param%lattice_type == circular_lattice$) then
    do i = i1+1, lat%n_ele_track
      call add_on_to_t_map
    enddo
    do i = 1, i2
      call add_on_to_t_map
    enddo

! Linear lattice with i1 > i2: Track backwards.

  else
    if (.not. unit_start_this) then
      call out_io (s_fatal$, r_name, 'Backwards propagation with a non-unit starting map!')
      call err_exit
    endif

    do i = i2+1, i1
      call add_on_to_t_map
    enddo
    call taylor_inverse (t_map, t_map)

  endif

!--------------------------------------------------------
contains

subroutine add_on_to_t_map

  integer k

  ! match, lcavity and taylor elements do not have corresponding ptc fibre elements.
  ! In this case we must concat.

  k = lat%ele(i)%key
  if (integrate_this .and. k /= lcavity$ .and. k /= match$ .and. k /= taylor$) then
    call taylor_propagate1 (t_map, lat%ele(i), lat%param)
  else
    if (.not. associated(lat%ele(i)%taylor(1)%term)) then
      call ele_to_taylor (lat%ele(i), lat%param)
    endif

    call concat_ele_taylor (t_map, lat%ele(i), t_map)
  endif

end subroutine

end subroutine
