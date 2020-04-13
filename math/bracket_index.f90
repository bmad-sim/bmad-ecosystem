!+
! Function bracket_index (s, s_arr, i_min, dr) result (ix)
!
! Function to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! If s <  s_arr(i_min) then ix = i_min - 1
! If s >= s_arr(i_max) then ix = i_max
! If array size zero then ix = i_min - 1
!
! This routine assumes that s_arr is in assending order.
! Also see bracket_index2
!
! Input:
!   s             -- real(rp): Number to bracket.
!   s_arr(i_min:) -- real(rp): Sequence of real numbers.
!   i_min         -- integer: Lower bound of s_arr
!
! Output:
!   ix        -- integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!   dr        -- real(rp), optional: Fraction (s-s_arr(ix)) / (s_arr(ix+1)-s_arr(ix))
!                   Undefined if s is outside of the region.
!-

function bracket_index (s, s_arr, i_min, dr) result (ix)

use precision_def

implicit none

integer i_min, i_max
real(rp) s_arr(i_min:), s
real(rp), optional :: dr

integer ix, n1, n2, n3

! Easy cases

if (present(dr)) dr = 0
i_max = ubound(s_arr, 1)

if (size(s_arr) == 0) then
  ix = i_min - 1
  return
endif

if (s < s_arr(i_min)) then
  ix = i_min - 1
  return
endif

if (s >= s_arr(i_max)) then
  ix = i_max
  return
endif

!

n1 = i_min
n3 = i_max

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

