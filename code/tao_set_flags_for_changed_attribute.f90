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
!                 May be null, for example, if ele_name = "BEAM_START".
!   val_ptr  -- real(rp):, pointer, optional: Pointer to the attribute that was changed.
!-

subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr)

use tao_mod, dummy => tao_set_flags_for_changed_attribute

implicit none

type (tao_universe_struct) u
type (ele_struct), pointer, optional :: ele_ptr

real(rp), pointer, optional :: val_ptr

character(*) ele_name

! Beam_start. 

if (ele_name == 'BEAM_START') then
  u%beam%beam_init%center            = u%model%lat%beam_start%vec

  u%model%tao_branch(0)%orbit(0)%vec = u%model%lat%beam_start%vec
  u%model%tao_branch(0)%orb0%vec     = u%model%lat%beam_start%vec

  u%model%tao_branch(0)%orbit(0)%t   = u%model%lat%beam_start%t
  u%model%tao_branch(0)%orb0%t       = u%model%lat%beam_start%t

  u%model%tao_branch(0)%orbit(0)%p0c   = u%model%lat%beam_start%p0c
  u%model%tao_branch(0)%orb0%p0c       = u%model%lat%beam_start%p0c

  u%beam%init_beam0 = .true.
endif

!

if (present(ele_ptr)) then
  if (associated(ele_ptr)) call set_flags_for_changed_attribute (ele_ptr, val_ptr)
endif

u%calc%lattice = .true.

end subroutine tao_set_flags_for_changed_attribute
