!+
! Subroutine bracket_index2 (s_arr, i_min, i_max, s, ix0, ix)
!
! Subroutine to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! If s <  s_arr(i_min) then ix = i_min - 1
! If s >= s_arr(i_max) then ix = i_max
! If i_max < i_min (array size zero) then ix = i_min - 1
!
! This routine assumes that s_arr is in assending order.
!
! This routine uses ix0 as an initial guess.
! See also bracket_index which does not need an initial guess.
!
! The reason for using this routine over bracket_index is speed when 
! dealing with large arrays.
! However, if ix0 is not near ix, this routine can be up to a 
! factor of 2 slower than bracket_index.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   s_arr(i_min:) -- Real(rp): Sequence of real numbers.
!   i_min         -- Integer: lower bound of s_arr
!   i_max         -- Integer: upper bound of s_arr
!   s             -- Real(rp): Number to bracket.
!   ix0           -- Integer: Initial guess for ix. 
!
! Output:
!   ix    -- Integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!-

subroutine bracket_index2 (s_arr, i_min, i_max, s, ix0, ix)

use precision_def

implicit none

real(rp) s_arr(i_min:), s

integer i_min, i_max
integer ix0, ix, n1, n2, n3, n_del

! Easy case

if (i_max < i_min) then
  ix = i_min - 1
  return
endif

! Start from ix0

n2 = min(max(ix0, i_min), i_max)
n_del = 1

! Bracket solution

if (s < s_arr(n2)) then
  n3 = n2
  do
    n1 = max(n3 - n_del, i_min)
    if (s >= s_arr(n1)) exit
    if (n1 == i_min) then
      ix = i_min - 1
      return
    endif
    n3 = n1
    n_del = 2 * n_del
  enddo

else
  n1 = n2
  do 
    n3 = min(n1 + n_del, i_max)
    if (s < s_arr(n3)) exit
    if (n3 == i_max) then
      ix = i_max
      return
    endif
    n1 = n3
    n_del = 2 * n_del
  enddo

endif

! Now find solution by successive divitions.

do

  if (n3 == n1 + 1) then
    ix = n1
    return
  endif

  n2 = (n1 + n3) / 2

  if (s < s_arr(n2)) then
    n3 = n2
  else
    n1 = n2
  endif

enddo

end subroutine

