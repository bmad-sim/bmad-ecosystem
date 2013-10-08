integer function lenstr(string)
!  LENSTR      SUBROUTINE  STRING      C.TYPSCN    TAP         95.08.02
!tpl  length=lenstr(string)  !return len excl trailing blanks

implicit none
character(*) string
integer last_position, length
logical repeat

repeat = .TRUE.
length = len(string)
last_position = length

do while(repeat)
  if(ichar(string(last_position:last_position)) <= 32) then
     last_position = last_position - 1
     if(last_position < 1)  repeat = .FALSE.
  else
     repeat = .FALSE.
  endif
enddo

lenstr = last_position

return
end
