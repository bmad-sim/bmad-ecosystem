!+
! Subroutine ele_at_s (lat, s, ix_ele)
!
! Subroutine to return the index of the element at position s.
! That is, ix_ele is choisen such that:
!     lat%ele_(ix_ele-1)%s < s <= lat%ele_(ix_ele)%s
!
! Note: s is evaluated modulo the lattice length:
!     s -> s - lat_length * floor(s/lat_length)
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- Ring_struct: Lattice of elements.
!   s   -- Real(rp): Longitudinal position.
!
! Output:
!   ix_ele -- Integer: Index of element at s.
!-

subroutine ele_at_s (lat, s, ix_ele)

  use bmad

  implicit none

  type (ring_struct) lat
  real(rp) s, ss, ll
  integer ix_ele, n1, n2, n3

!

  ll = lat%param%total_length
  ss = s - ll * floor(s/ll)

  n1 = 0
  n3 = lat%n_ele_use

  do

    if (n3 == n1 + 1) then
      ix_ele = n1
      return
    endif

    n2 = (n1 + n3) / 2

    if (ss < lat%ele_(n2)%s) then
      n3 = n2
    else
      n1 = n2
    endif

  enddo

end subroutine
