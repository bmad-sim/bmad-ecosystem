!+
! Subroutine tao_remove_blank_characters (str)
!
! Routine to remove blank characters from a string.
! This is used for graph and curve %component parameters to get them into a standard form.
! EG: "model - design" => "model-design"
!
! Input:
!   str   -- character(*): Input string.
!
! Output:
!   str   -- character(*): String with blank characters removed.
!-

subroutine tao_remove_blank_characters (str)

implicit none

character(*) str
integer i, j

!

j = 1
do i = 1, len_trim(str)
  if (str(j:j) == ' ') then
    str = str(:j-1) // str(j+1:)
  else
    j = j + 1
  endif
enddo

end subroutine tao_remove_blank_characters
