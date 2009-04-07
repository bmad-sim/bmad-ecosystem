!+
! Subroutine bracket_index (s_arr, i_min, i_max, s, ix)
!
! Subroutine to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! If s <  s_arr(i_min) then ix = i_min - 1
! If s >= s_arr(i_max) then ix = i_max  [where i_max = ubound(s_arr)]
!
! This routine assumes that s_arr is in assending order.
!
! Modules needed:
!   use dcslib
!
! Input:
!   s_arr(i_min:) -- Real(rp): Sequence of real numbers.
!   i_min         -- Integer: lower bound of s_arr
!   i_max         -- Integer: upper bound of s_arr
!   s             -- Real(rp): Number to bracket.
!
! Output:
!   ix    -- Integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!-

subroutine bracket_index (s_arr, i_min, i_max, s, ix)

  use precision_def

  implicit none

  integer i_min, i_max
  real(rp) s_arr(i_min:), s

  integer ix, n1, n2, n3

!

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

