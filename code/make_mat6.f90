!+
! Subroutine make_mat6 (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- Param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct, optional: Coordinates at the beginning of element. 
!               If not present then effectively c0 = 0.
!
! Output:
!   ele    -- Ele_struct: Element
!     %mat6  -- Real(rdef): 6x6 transfer matrix.
!   c1     -- Coord_struct, optional: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6 (ele, param, c0, c1)

  use bmad

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct), optional :: c0, c1
  type (param_struct)  param
  type (coord_struct) c00, c11

  integer mat6_calc_method

!--------------------------------------------------------
! init

  call attribute_bookkeeper (ele, param)

  mat6_calc_method = ele%mat6_calc_method
  if (.not. ele%is_on) mat6_calc_method = bmad_standard$

!

  if (present(c0)) then
    c00 = c0
  else
    c00%vec = 0
  endif

  select case (mat6_calc_method)

  case (taylor$)
    call make_mat6_taylor (ele, param, c00, c11)

  case (runge_kutta$)
    call make_mat6_runge_kutta (ele, param, c00, c11)

  case (custom$) 
    call make_mat6_custom (ele, param, c00, c11)

  case (bmad_standard$)
    call make_mat6_bmad (ele, param, c00, c11)

  case (symp_lie$)
    call make_mat6_symp_lie (ele, param, c00, c11)

  case (tracking$)
    call make_mat6_tracking (ele, param, c00, c11)

  case (none$)
    return

  case default
    print *, 'ERROR IN MAKE_MAT6: UNKNOWN MAT6_CALC_METHOD: ', &
                                    calc_method_name(ele%mat6_calc_method)
    call err_exit
  end select

! symplectify if wanted

  if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6)

! make the 0th order transfer vector

  ele%vec0 = c11%vec - matmul(ele%mat6, c00%vec)
  if (present (c1)) c1 = c11

end subroutine

