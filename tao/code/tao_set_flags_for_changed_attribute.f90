!+
! Subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr, who)
!
! Routine to set flags in the model lattice indicating that a parameter value has changed.
! Call this routine *after* setting the variable.
!
! Input:
!   u        -- tao_universe_sturct: Universe containing the lattice.
!   ele_name -- character(*): Associated "element" of the changed parameter.
!   ele_ptr  -- ele_struct, pointer, optional: Pointer to the element. 
!                 May be null, for example, if ele_name = "PARTICLE_START".
!   val_ptr  -- all_pointer_struct: optional: Pointer to the attribute that was changed.
!                 Must be present if ele_ptr is present.
!   who      -- character(*), optional: Name of changed attribute. Only used with PARTICLE_START.
!-

subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr, who)

use tao_interface, dummy => tao_set_flags_for_changed_attribute
use bookkeeper_mod, only: set_flags_for_changed_attribute

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer, optional :: ele_ptr
type (lat_struct), pointer :: lat
type (all_pointer_struct), optional :: val_ptr

integer ib, ie, n_loc, ix_branch
logical err

character(*) ele_name
character(*), optional :: who

! If the beginning element is modified, need to reinit any beam distribution.

u%calc%lattice = .true.
lat => u%model%lat

if (ele_name == 'PARTICLE_START') then
  if (lat%branch(0)%param%geometry == closed$) then
    u%model%tao_branch(0)%orb0%vec(6) = lat%particle_start%vec(6)
  endif
  ! For photon bookkeeping
  select case (who)
  case ('PZ');        lat%particle_start%p0c = 0
  case ('E_PHOTON');  lat%particle_start%vec(6) = 0
  case ('T');         lat%particle_start%vec(5) = real_garbage$   ! So time is used instead of z.
  end select
  return
endif

if (present(ele_ptr)) then
  if (associated(ele_ptr)) then
    ix_branch = ele_ptr%ix_branch

    if (ele_ptr%key == ramper$) then
      call apply_all_rampers(lat, err)
    else
      if (ele_ptr%ix_ele <= u%model_branch(ix_branch)%beam%ix_track_start) &
                               u%model_branch(ix_branch)%beam%init_starting_distribution = .true.
      call set_flags_for_changed_attribute (ele_ptr, val_ptr)
    endif
  endif
endif

end subroutine tao_set_flags_for_changed_attribute
