!+
! Subroutine ele_at_s (lat, s, ix_ele)
!
! Subroutine to return the index of the element at position s.
! That is, ix_ele is choisen such that:
!     lat%ele(ix_ele-1)%s < s <= lat%ele(ix_ele)%s
!
! Note: For a circular lattice s is evaluated modulo the 
! lattice length:
!     s -> s - lat_length * floor(s/lat_length)
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lattice of elements.
!   s   -- Real(rp): Longitudinal position.
!
! Output:
!   ix_ele -- Integer: Index of element at s.
!-

subroutine ele_at_s (lat, s, ix_ele)

  use bmad

  implicit none

  type (lat_struct) lat
  real(rp) s, ss, ll
  integer ix_ele, n1, n2, n3
  character(16) :: r_name = 'ele_at_s'

!

  ll = lat%param%total_length

  if (lat%param%lattice_type == circular_lattice$) then
    ss = s - ll * floor((s-lat%ele(0)%s)/ll)
  else
    ss = s
    if (s < lat%ele(0)%s .or. s > lat%ele(lat%n_ele_track)%s) then
      call out_io (s_fatal$, r_name, 'S POSITION OUT OF BOUNDS \f10.2\ ' , s)
      call err_exit
    endif
  endif

  n1 = 0
  n3 = lat%n_ele_track

  do

    if (n3 == n1 + 1) then
      ix_ele = n3
      return
    endif

    n2 = (n1 + n3) / 2

    if (ss < lat%ele(n2)%s) then
      n3 = n2
    else
      n1 = n2
    endif

  enddo

end subroutine
