!+
! Subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr)
!
! Routine to set flags in the model lattice indicating that a parameter value has changed.
! Call this routine *after* setting the variable.
!
! Input:
!   u        -- tao_universe_sturct: Universe containing the lattice.
!   ele_name -- character(*): Associated "element" of the changed parameter.
!   ele_ptr  -- ele_struct, pointer, optional: Pointer to the element. 
!                 May be null, for example, if ele_name = "PARTICLE_START".
!   val_ptr  -- real(rp):, pointer, optional: Pointer to the attribute that was changed.
!-

subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr)

use tao_interface, dummy => tao_set_flags_for_changed_attribute
use bookkeeper_mod, only: set_flags_for_changed_attribute

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer, optional :: ele_ptr
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:) 

real(rp), pointer, optional :: val_ptr
integer ib, ie, n_loc
logical err

character(*) ele_name

! If the beginning element is modified, need to reinit any beam distribution.

u%calc%lattice = .true.
lat => u%model%lat

if (ele_name == 'PARTICLE_START') then
  if (lat%branch(0)%param%geometry == closed$) then
    u%model%tao_branch(0)%orb0%vec(6) = lat%particle_start%vec(6)
  endif
  return
endif

if (present(ele_ptr)) then
  if (associated(ele_ptr)) then

    if (ele_ptr%key == ramper$) then
      call lat_ele_locator ('RAMPER::*', lat, eles, n_loc, err)
      do ib = 0, ubound(lat%branch, 1)
        branch => lat%branch(ib)
        do ie = 0, branch%n_ele_max
          call apply_ramper(branch%ele(ie), eles(1:n_loc), err)
        enddo
      enddo
      call set_flags_for_changed_attribute (lat)

    else
      if (ele_ptr%ix_ele == 0) u%model_branch(0)%beam%init_starting_distribution = .true.
      if (present(val_ptr)) call set_flags_for_changed_attribute (ele_ptr, val_ptr)
    endif

  endif
endif

end subroutine tao_set_flags_for_changed_attribute
