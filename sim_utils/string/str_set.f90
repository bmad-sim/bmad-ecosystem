!+
! Subroutine str_set(str_out, str_in)
!
! Routine to set a variable length string.
! Trailing blanks will be trimmed.
! 
! Input:
!   str_in      -- character(*): String to set to
!
! Output:
!   str_out     -- character(:), allocatable: Variable length string.
!-

subroutine str_set(str_out, str_in)

implicit none

character(:), allocatable :: str_out
character(*) str_in
integer n

!

n = len_trim(str_in)
if (.not. allocated(str_out)) allocate (character(n):: str_out)
str_out = str_in(1:n)

end subroutine
