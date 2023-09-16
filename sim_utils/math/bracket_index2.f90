!+
! Function bracket_index2 (s, ix0, s_arr, i_min, dr, restrict) result (ix)
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
!   restrict      -- logical, optional: Restrict output ix to be in the range [i_min, i_max-1]?
!                     Default = False.
!
! Output:
!   ix         -- integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!   dr         -- real(rp), optional: Fraction (s-s_arr(ix)) / (s_arr(ix+1)-s_arr(ix))
!                   Undefined if ix is outside the domain [i_min, i_max-1].
!                   Set to zero if s_arr(ix) == s_arr(ix+1).
!-

function bracket_index2 (s, ix0, s_arr, i_min, dr, restrict) result (ix)

use sim_utils, dummy => bracket_index2

implicit none

real(rp) s_arr(i_min:), s
real(rp), optional :: dr

integer i_min, i_max
integer ix0, ix, n1, n2, n3, n_del
logical, optional :: restrict

! Easy cases

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
      if (i_max > i_min .and. logic_option(.false., restrict)) then
        ix = i_min
        if (present(dr)) then
          if (s_arr(ix) /= s_arr(ix+1)) dr = (s - s_arr(ix)) / (s_arr(ix+1) - s_arr(ix))
        endif
      else
        ix = i_min - 1
      endif
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
      if (i_max > i_min .and. logic_option(.false., restrict)) then
        ix = i_max - 1
        if (present(dr)) then
          if (s_arr(ix) /= s_arr(ix+1)) dr = (s - s_arr(ix)) / (s_arr(ix+1) - s_arr(ix))
        endif
      else
        ix = i_max
      endif
      return
    endif
    n1 = n3
    n_del = 2 * n_del
  enddo

endif

! Now find solution by successive divisions.

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

