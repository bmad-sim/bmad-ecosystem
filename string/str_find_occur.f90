!+
! Subroutine str_find_occur (string, str_match, locs,num)
! Routine to find and count all occurances of a char in a string
! Input:
!   string      -- Character(*): Character string.
!   str_match   -- Character(1): character to find. 
! Output:
!   locs     integer*(*) locations in string where str_match found
!   num      integer number of occurances found ;
!   note if str_match=' ', trailing blanks are not counted 
!-

     subroutine str_find_occur(string, str_match, locs,num)
     implicit none
     character(*) string,str_match*1
     integer     i,locs(*), num,n,len
     logical not_notif
     not_notif=.true.
     len=len_trim(string) ; num=0 
     do i=1,len
      if(string(i:i).eq.str_match) then
       num=num+1 ; locs(num)=i
      endif
     enddo
     end subroutine
