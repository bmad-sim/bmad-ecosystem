!+
! Function str_match_wild(str, pat) result (a_match)
!
! Function to match a character string against a regular expression pattern.
!
! Note: trailing blanks are included in the match!
! Use the routine match_wild if trailing blanks are to be ignored.
!
! Wild card characters are:
!  "*"  -- Match to any number of characters.
!  "%"  -- Match to any single character.
!
! Input:
!   str -- Character(*): String to match against.
!   pat -- Character(*): Regular expression pattern.
!
! Output:
!   a_match -- Integer:
!                 0 -> No match.
!                >0 -> Index of where match begins.
!-

recursive function str_match_wild(str, pat) result (a_match)

implicit none
  
character(*) pat, str
integer p, s, s2, plen, slen
logical a_match
  
!

plen = len(pat)
slen = len(str)
s = 1
p = 1
a_match = .true.
  
do while ((p <= plen) .and. (s <= slen) .and. (a_match))

  select case (pat(p:p))

  case('%')
    s = s + 1

  case('*')
    if (p == plen) then
       a_match = .true.
       return
    else
      s2 = s
      a_match = .false.
      do s2 = s,slen
        if (str_match_wild(str(s2:), pat(p + 1:))) then
           a_match = .true.
        endif
      enddo
      return
    endif

  case default
    if (pat(p:p) /= str(s:s)) then
      a_match = .false.
    endif
    s = s + 1
  end select

  p = p + 1

enddo
    
do while (a_match .and. (p <= plen) )
 if (pat(p:p) /= '*') exit
  p = p + 1
end do
    
a_match = ((s > slen) .and. (p > plen) .and. a_match)

end function
