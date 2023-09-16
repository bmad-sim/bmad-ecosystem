!+
! Function match_wild (string, template) result (this_match)
!
! Function to do wild card matches.
!
! Note: trailing blanks will be discarded before any matching is done. 
! If this is not desired, use str_match_wild to include trailing blanks.
!
! Wild card characters are:
!  "*"  -- Match to any number of characters.
!  "%"  -- Match to any single character.
!
! Input:
!   string   -- Character(*): String to be matched to.
!   template -- Character(*): String with wild cards.
!
! Output:
!   this_match -- Logical: Set true if matched.
!-

function match_wild (string, template) result (this_match)

use sim_utils, except => match_wild

implicit none

character(*) string, template
integer is, it
logical this_match

!

is = max(1, len_trim(string))
it = max(1, len_trim(template))

this_match = str_match_wild (string(:is), template(:it))

end function match_wild
