!+
! Function map_coef(y, i, j, k, l, style)
!
! Function to return the coefficient of the map y(:) up to 3rd order.
! Example: 
!   In 4-dimensional space with
!      X = (x_1, x_2, x_3, x_4)
!   And if y(:) is the map of X(in) to X(out) then
!      map_coef(y, 2, 1, 4)
!   Gives the coefficient for
!     x_2(out) = ... + coef * x_1(in) * x_4(in) + ...
!
! Notice that map_coef(y, i, j) just gives the linear (matrix) part of the map.
!
! Modules Needed:                    
!   use accelerator
!
! Input:
!   y(:)  -- Real_8: Taylor Map.
!   i     -- Integer: output index.
!   j     -- Integer, optional: 1st input index, needed for 1st order and above.
!   k     -- Integer, optional: 2nd input index, needed for 2nd order and above.
!   l     -- Integer, optional: 3rd input index, needed for 3rd order.
!   style -- Integer, optional: bmad_std$ or ptc_std$.
!             bmad_std$ --> (x, p_x, y, p_y, z, dE/E) for phase space vars.
!             ptc_std$  --> (x, p_x, y, p_y, dE/E, ct ~ -z)
!             default = ptc_std$.
!
! Output:
!   map_coef -- Real*8: Coefficient.
!-

function map_coef (y, i, j, k, l, style) result (m_coef)

  use accelerator

  implicit none

  type (real_8) y(:)


  real*8 m_coef

  integer i
  integer, optional :: j, k, l, style
  integer arr(40), n_max, sgn, ii

  character str*40
  character, parameter :: str1(4) = (/ '1', '2', '3', '4' /)

  logical use_bmad

!

  use_bmad = .false.
  if (present(style)) then
    if (style == bmad_std$) use_bmad = .true.
  endif

  arr = 0
  sgn = 1
  str = '0000000000000000000000000000000000000000'
  n_max = 1

  call map_index(j)
  call map_index(k)
  call map_index(l)

  call map_index(i, ii)
  m_coef = sgn * y(ii)%t.sub.str(1:n_max)

!--------------------------------------------------------------
contains

subroutine map_index (iz, i_in)
  
  integer, optional :: iz, i_in
  integer n0

!

  if (.not. present(iz)) return

    n0 = iz

    if (use_bmad) then
      if (iz == 5) then
        n0 = 6
        sgn = -sgn
      elseif (iz == 6) then
        n0 = 5
      endif
    endif

! i_in is present only with the input index.

    if (present(i_in)) then
      i_in = n0
      return
    endif

    arr(n0) = arr(n0) + 1
    str(n0:n0) = str1(arr(n0))

    n_max = max(n_max, n0)

end subroutine

end function
