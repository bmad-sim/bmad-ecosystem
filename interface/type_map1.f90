!+
! Subroutine type_map1 (y, type0, n_dim, style)
!
! Subroutine to type the transfer map up to first order.
!
! Modules needed:
!   use accelerator
!
! Input:
!   y(:)  -- Real_8: 
!   type0 -- Logical: Type zeroth order map
!   n_dim -- Integer: Number of dimensions to type: 4 or 6
!   style -- Integer, optional: bmad_std$ or ptc_std$.
!             bmad_std$ --> (x, p_x, y, p_y, z, dE/E) for phase space vars.
!             ptc_std$  --> (x, p_x, y, p_y, dE/E, ct ~ -z) 
!             default = ptc_std$.
!-

subroutine type_map1 (y, type0, n_dim, style)

  use accelerator

  implicit none

  type (real_8), intent(in) :: y(:)

  integer, intent(in) :: n_dim
  integer, optional, intent(in) :: style
  integer :: i, j

  logical, intent(in) :: type0

!

  if (type0) then
    type *, '0th Order Map:'
    type '(6f11.5)', (map_coef(y(:), i, style=style), i = 1, n_dim)
    type *
  endif

  type *, '1st Order Map:'
  do i = 1, n_dim
    type '(6f11.5)', (map_coef(y(:), i, j, style=style), j = 1, n_dim)
  enddo

end subroutine
