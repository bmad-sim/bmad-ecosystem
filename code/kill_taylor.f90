!+
! Subroutine kill_taylor (bmad_taylor)
!
! Subroutine to deallocate a BMAD taylor_struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: New structure.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: deallocated structure.
!-

subroutine kill_taylor (bmad_taylor)

  use bmad

  implicit none

  type (taylor_struct) bmad_taylor(:)

  integer i, ix

!

  do i = 1, size(bmad_taylor)
    deallocate (bmad_taylor(i)%term, stat = ix)
  enddo

end subroutine
