!+
! Subroutine taylor_to_mat6 (a_taylor, c0, mat6, c1)
!
! Subroutine to calculate the linear (Jacobian) matrix about 
! some trajectory from a taylor series.
!
! Modules needed:
!   use bmad
!
! Input:
!   a_taylor(6) -- Taylor_struct: Taylor series.
!   c0          -- Coord_struct: Coordinates at the input. 
!
! Output:
!   mat6(6,6) -- Real(rdef): Jacobian.
!   c1        -- Coord_struct: Oth order transport vector.
!-

subroutine taylor_to_mat6 (a_taylor, c0, mat6, c1)

  use bmad

  implicit none

  type (taylor_struct), target, intent(in) :: a_taylor(6)
  type (coord_struct), intent(in) :: c0
  type (coord_struct), intent(out) :: c1
  type (taylor_term_struct), pointer :: term

  real(rdef), intent(out) :: mat6(6,6)
  real(rdef) prod

  integer i, j, k, l

! mat6 calc

  mat6 = 0
  
  do i = 1, 6
    do j = 1, 6
      do k = 1, size(a_taylor(i)%term)

        term => a_taylor(i)%term(k)

        if (term%exp(j) == 0) cycle
        if (term%exp(j) > 1 .and. c0%vec(j) == 0) cycle

        if (term%exp(j) > 1)then
          prod = term%coef * term%exp(j) * c0%vec(j) ** (term%exp(j)-1)
        else
          prod = term%coef
        endif

        do l = 1, 6
          if (term%exp(l) == 0) cycle
          if (l == j) cycle
          prod = prod * c0%vec(l) ** term%exp(l)
        enddo

        mat6(i,j) = mat6(i,j) + prod

      enddo
    enddo
  enddo

! c1 calc

  c1%vec = 0

  do i = 1, 6
    do k = 1, size(a_taylor(i)%term)
      if (all(a_taylor(i)%term(k)%exp == 0)) then
        c1%vec(i) = a_taylor(i)%term(k)%coef
        exit
      endif
    enddo
  enddo

end subroutine
