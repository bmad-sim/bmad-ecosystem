!+
! Function bracket_index2 (s, ix0, s_arr, i_min, dr) result (ix)
!
! Function to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! If s <  s_arr(i_min) then ix = i_min - 1
! If s >= s_arr(i_max) then ix = i_max
! If array size zero, ix = i_min - 1
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
! Input:
!   s             -- real(rp): Number to bracket.
!   ix0           -- integer: Initial guess for ix. 
!   s_arr(i_min:) -- real(rp): Sequence of real numbers.
!   i_min         -- integer: lower bound of s_arr.
!
! Output:
!   ix         -- integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!   dr        -- real(rp), optional: Fraction (s-s_arr(ix)) / (s_arr(ix+1)-s_arr(ix))
!                   Undefined if s is outside of the region.
!-

function bracket_index2 (s, ix0, s_arr, i_min, dr) result (ix)

use precision_def

implicit none

real(rp) s_arr(i_min:), s
real(rp), optional :: dr

integer i_min, i_max
integer ix0, ix, n1, n2, n3, n_del

! Easy case

if (present(dr)) dr = 0
i_max = ubound(s_arr, 1)

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
    if (present(dr)) then
      if (s_arr(ix) /= s_arr(ix+1)) dr = (s - s_arr(ix)) / (s_arr(ix+1) - s_arr(ix))
    endif
    return
  endif

  n2 = (n1 + n3) / 2

  if (s < s_arr(n2)) then
    n3 = n2
  else
    n1 = n2
  endif
enddo

end function

