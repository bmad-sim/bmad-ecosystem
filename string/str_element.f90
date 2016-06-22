!+
! Subroutine str_element(dest_string, element_number, delim_string, source_string)
!
! Subroutine to trim a string of leading blanks and/or tabs and also to 
! return the length of the first word.
!
! Input:
!  element_number -  integer: Element # to return (0 is first element, 1 is second)
!  source_string  - Character(*): String constant to be parsed.
!  delim_string   - Character(*): String constant to be parsed.
!
! Output:
!   dest_string - Character(*): String constant: Left shifted to trim leading blanks and/or tabs.
!-

subroutine str_element(dest_string, element_number, delim_string, source_string)

implicit none

character(*) :: source_string, delim_string, dest_string
character :: delim*1, temp_string*120
integer :: i, element_number, ix_word
logical :: delim_found

ix_word = 0
dest_string = ''
temp_string = source_string
do i = 0, element_number
  call word_read(temp_string, delim_string, dest_string, &
       ix_word, delim, delim_found, temp_string)
enddo

end subroutine str_element




