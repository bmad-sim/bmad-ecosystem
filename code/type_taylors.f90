!+
! Subroutine type_taylors (bmad_taylor)
!
! Subroutine to type out the taylor series from a BMAD taylor array
!
! Modules needed:
!   use accelerator
!
! Input
!   bmad_taylor(6) -- Taylor_struct: 6 taylor series: (x, P_x, y, P_y, z, P_z) 
!-

subroutine type_taylors (bmad_taylor)

  use bmad

  implicit none

  type (taylor_struct), target :: bmad_taylor(:)
  integer i, n_lines
  character*80, pointer :: lines(:)

!

  nullify (lines)

  call type2_taylors (bmad_taylor, lines, n_lines)

  do i = 1, n_lines
    print *, trim(lines(i))
  enddo

  deallocate(lines)

end subroutine
