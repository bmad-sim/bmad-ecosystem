!+
! Subroutine init_taylor (bmad_taylor)
!
! Subroutine to initialize a BMAD taylor_struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: New structure.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Initalized structure.
!-

subroutine init_taylor (bmad_taylor)

  use bmad

  implicit none

  type (taylor_struct) bmad_taylor(:)

  integer i

!

  do i = 1, size(bmad_taylor)
    nullify (bmad_taylor(i)%term)
  enddo

end subroutine
