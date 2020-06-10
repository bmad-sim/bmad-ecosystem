!+
! Function bracket_index_int (s, int_arr, i_min, restrict) result (ix)
!
! Function to find the index ix of int_arr(i_min:i_max) so that int_arr(ix) <= s < int_arr(ix+1).
! Boundary cases:
!   If s <  int_arr(i_min) then ix = i_min - 1
!   If s >= int_arr(i_max) then ix = i_max
!   If array size zero then ix = i_min - 1
! Exception: If restrict_index = True, and i_max > i_min, ix is restricted to the domain [i_min, i_max-1].
!   If s <  int_arr(i_min+1) then ix = i_min 
!   If s >= int_arr(i_max-1) then ix = i_max-1
!
! This routine assumes that int_arr is in assending order.
!
! Input:
!   s               -- integer: Number to bracket.
!   int_arr(i_min:) -- integer: Sequence of integers.
!   i_min           -- integer: Lower bound of int_arr.
!   restrict        -- logical, optional: Restrict output ix to be in the range [i_min, i_max-1]?
!                       Default = False.
!
! Output:
!   ix        -- integer: Index so that int_arr(ix) <= s < int_arr(ix+1).
!-

function bracket_index_int (s, int_arr, i_min, restrict) result (ix)

use sim_utils, dummy => bracket_index_int

implicit none

integer i_min, i_max
integer int_arr(i_min:), s

integer ix, n1, n2, n3
logical, optional :: restrict

! Easy cases

i_max = ubound(int_arr, 1)

if (size(int_arr) == 0) then
  ix = i_min - 1
  return
endif

if (i_max > i_min .and. logic_option(.false., restrict)) then
  if (s < int_arr(i_min)) then
    ix = i_min
    return

  elseif (s >= int_arr(i_max)) then
    ix = i_max-1
    return
  endif

else
  if (s < int_arr(i_min)) then
    ix = i_min - 1
    return

  elseif (s >= int_arr(i_max)) then
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
    return
  endif

  n2 = (n1 + n3) / 2

  if (s < int_arr(n2)) then
    n3 = n2
  else
    n1 = n2
  endif
enddo

end function

