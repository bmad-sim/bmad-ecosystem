!+
! Subroutine read_string
! 
! Routine to read a string at the command line 
! On a carriage return, with no string supplied, string = default_string
!
! Input:
!  character(*) -- query string
!  character(*) -- default_string
!  character(*) -- string
!
! 
subroutine  read_string(query_str,default_string, string)

implicit none

 character(*) query_str, string, default_string
 character(100) temp
 integer ix

  print '(a,$)', query_str
  read(*,'(a)') string
  call string_trim(string, string, ix)
  if(ix == 0 .or. ix == len(string))then
     string = default_string
  endif
 return
 end
