!  ABS_SORT    SUBROUTINE  MATH        C.DCS.LIB   DCS         96.7
!+
! Subroutine ABS_SORT (ARRAY, INDEX, N)
!
! Input:
!     Real(rp)  ARRAY(*)    -- Array to be sorted by absolute value.
!     Int   N           -- Sort elements ARRAY(1) through ARRAY(N).
!
! Output
!     Int   INDEX(*)    -- INDEX(1) points to smallest abs(ARRAY(I))
!                       -- INDEX(N) points to largest abs(ARRAY(I))
!
! This subroutine indexes an array using the absolute value of the
! elements of the array. this subroutine is essentually identical to
! the subroutine INDEXX in "Numerical Recipes" by Press et. al
!-

!$Id$
!$Log$
!Revision 1.4  2003/07/09 01:29:29  dcs
!new bmad
!
!Revision 1.3  2002/02/23 20:34:38  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine abs_sort (array, index, n)


  use precision_def

  implicit none
  real(rp) array(*), q
  integer index(*), l, ir, n, i, j, indxt

  do j = 1, n
    index(j) = j
  enddo

  l = n/2 + 1
  ir = n

10    continue

  if (l > 1) then
    l = l - 1
    indxt = index(l)
    q = abs(array(indxt))
  else
    indxt = index(ir)
    q = abs(array(indxt))
    index(ir) = index(1)
    ir = ir - 1
    if (ir == 1) then
      index(1) = indxt
      return
    endif
  endif
  i = l
  j = l + l

20    if (j <= ir) then
    if (j < ir) then
      if (abs(array(index(j))) < abs(array(index(j + 1)))) &
                                                        j = j + 1
    endif
    if (q < abs(array(index(j)))) then
      index(i) = index(j)
      i = j
      j = j + j
    else
      j = ir + 1
    endif
    go to 20
  endif
  index(i) = indxt
  go to 10
  end
