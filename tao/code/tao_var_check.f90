!+
! Subroutine tao_var_check (eles, attribute, silent)
!
! Routine to check that setting element attributes has caused a problem.
!
! A problem arises if a Tao variable has multiple slaves and only some of them got set. 
! For example, if a variable controls multiple elements named "Q" and if only "Q##2" 
! has been changed there is a problem.
! 
! This routine will try to fix the problem if possible.
!
! Input:
!   eles(:)     -- ele_pointer_struct: Array of elements which have a changed attribute.
!   attribute   -- character(*): Name of attribute changed.
!   silent      -- logical: If True and the problem can be fixed, do not issue an error message.
!-

subroutine tao_var_check (eles, attribute, silent)

use tao_interface, dummy => tao_var_check

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (tao_universe_struct), pointer :: u
type (tao_var_slave_struct), pointer :: vs
type (all_pointer_struct) a_ptr
type (tao_var_struct), pointer :: var

integer i, iv, is, ix_found

logical silent, err, mismatch

character(*) attribute
character(*), parameter :: r_name = 'tao_var_check'

!

do i = 1, size(eles)
  if (.not. associated(eles(i)%ele)) cycle ! Can happen, for example with particle_start var.
  u => s%u(eles(i)%id)
  call pointer_to_attribute(eles(i)%ele, attribute, .true., a_ptr, err, .false.)
  if (err) cycle
  if (.not. associated(a_ptr%r)) cycle

  do iv = 1, s%n_var_used
    var => s%var(iv)
    if (size(var%slave) == 1) cycle
    mismatch = .false.
    ix_found = 1    ! Index of parameter that got changed. Default to 1.
    do is = 1, size(var%slave)
      vs => var%slave(is)
      if (vs%model_value /= var%slave(1)%model_value) mismatch = .true.
      if (associated(a_ptr%r, vs%model_value)) ix_found = is
    enddo

    if (mismatch) then
      if (.not. silent) then
        call out_io (s_warn$, r_name, 'TAO VARIABLE: ' // tao_var1_name(var), &
              'WHICH HAS MULTIPLE SLAVE PARAMETERS NOW DOES NOT HAVE ALL SLAVE PARAMETER VALUES THE SAME.', &
              'WILL FIX...')
      endif

      call tao_set_var_model_value(var, var%slave(ix_found)%model_value)
    endif

  enddo
enddo

end subroutine tao_var_check

