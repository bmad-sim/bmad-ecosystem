!+
! Subroutine kill_gen_field (gen_fieled)
!
! Subroutine to kill a gen_field.\
!
! Modules needed:
!   use bmad
!
! Input:
!   gen_field -- Genfield, pointer: gen_field to kill.
!
! Output:
!   gen_field -- Genfield, pointer: Killed gen_field.
!-

subroutine kill_gen_field (gen_field)

  use accelerator

  type (genfield), pointer :: gen_field

!

  if (associated(gen_field)) then
    call kill (gen_field)
    deallocate (gen_field)
  endif

end subroutine
