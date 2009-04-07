!+
! Function match_wild (string, template) result (this_match)
!
! Function to do wild card matches. Note: trailing blanks will be discarded
! before any matching is done. See str_match_wild for more details.
!
! Input:
!   string   -- Character(*): String to be matched do
!   template -- Character(*): String with wild cards
!
! Output:
!   this_match -- Logical: Set true if matched
!-

#include "CESR_platform.inc"

function match_wild (string, template) result (this_match)

  use sim_utils

  implicit none

  character(*) string, template
  integer is, it
  logical this_match

!

  is = max(1, len_trim(string))
  it = max(1, len_trim(template))

  this_match = str_match_wild (string(:is), template(:it))

end function match_wild
