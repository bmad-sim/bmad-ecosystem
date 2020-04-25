!+
! Function bracket_index (s, s_arr, i_min, dr, restrict) result (ix)
!
! Function to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! Boundary cases:
!   If s <  s_arr(i_min) then ix = i_min - 1
!   If s >= s_arr(i_max) then ix = i_max
!   If array size zero then ix = i_min - 1
! Exception: If restrict_index = True, and i_max > i_min, ix is restricted to the domain [i_min, i_max-1].
!   If s <  s_arr(i_min+1) then ix = i_min 
!   If s >= s_arr(i_max-1) then ix = i_max-1
!
! This routine assumes that s_arr is in assending order.
! Also see bracket_index2
!
! Input:
!   s             -- real(rp): Number to bracket.
!   s_arr(i_min:) -- real(rp): Sequence of real numbers.
!   i_min         -- integer: Lower bound of s_arr.
!   restrict      -- logical, optional: Restrict output ix to be in the range [i_min, i_max-1]?
!                     Default = False.
!
! Output:
!   ix        -- integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!   dr        -- real(rp), optional: Fraction (s-s_arr(ix)) / (s_arr(ix+1)-s_arr(ix))
!                   Undefined if ix is outside the domain [i_min, i_max-1].
!                   Set to zero if s_arr(ix) == s_arr(ix+1).
!-

function bracket_index (s, s_arr, i_min, dr, restrict) result (ix)

use sim_utils, dummy => bracket_index

implicit none

integer i_min, i_max
real(rp) s_arr(i_min:), s
real(rp), optional :: dr

integer ix, n1, n2, n3
logical, optional :: restrict

! Easy cases

if (present(dr)) dr = 0
i_max = ubound(s_arr, 1)

if (size(s_arr) == 0) then
  ix = i_min - 1
  return
endif

if (i_max > i_min .and. logic_option(.false., restrict)) then
  if (s < s_arr(i_min)) then
    ix = i_min
    if (present(dr)) then
      if (s_arr(ix) /= s_arr(ix+1)) dr = (s - s_arr(ix)) / (s_arr(ix+1) - s_arr(ix))
    endif
    return

  elseif (s >= s_arr(i_max)) then
    ix = i_max-1
    if (present(dr)) then
      if (s_arr(ix) /= s_arr(ix+1)) dr = (s - s_arr(ix)) / (s_arr(ix+1) - s_arr(ix))
    endif
    return
  endif

else
  if (s < s_arr(i_min)) then
    ix = i_min - 1
    return

  elseif (s >= s_arr(i_max)) then
    ix = i_max
    return
  endif
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

