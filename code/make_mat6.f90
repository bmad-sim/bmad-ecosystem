!+
! Subroutine make_mat6 (ele, param, start, end)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- Param_struct: Parameters are needed for some elements.
!   start  -- Coord_struct, optional: Coordinates at the beginning of element. 
!               If not present then effectively start = 0.
!
! Output:
!   ele    -- Ele_struct: Element
!     %mat6  -- Real(rp): 6x6 transfer matrix.
!   end    -- Coord_struct, optional: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6 (ele, param, start, end)

  use bmad_struct
  use bmad_interface
  use make_mat6_mod
  use symp_lie_mod
  use bookkeeper_mod

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct), optional :: start, end
  type (param_struct)  param
  type (coord_struct) a_start, a_end

  integer mat6_calc_method

!--------------------------------------------------------
! init

  call attribute_bookkeeper (ele, param)

  mat6_calc_method = ele%mat6_calc_method
  if (.not. ele%is_on) mat6_calc_method = bmad_standard$

!

  if (present(start)) then
    a_start = start
  else
    a_start%vec = 0
  endif

  select case (mat6_calc_method)

  case (taylor$)
    call make_mat6_taylor (ele, param, a_start, a_end)

  case (runge_kutta$)
    call make_mat6_runge_kutta (ele, param, a_start, a_end)

  case (custom$) 
    call make_mat6_custom (ele, param, a_start, a_end)

  case (bmad_standard$)
    call make_mat6_bmad (ele, param, a_start, a_end)

  case (symp_lie_ptc$)
    call make_mat6_symp_lie_ptc (ele, param, a_start, a_end)

  case (symp_lie_bmad$)
    call symp_lie_bmad (ele, param, a_start, a_end, .true.)

  case (tracking$)
    call make_mat6_tracking (ele, param, a_start, a_end)

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

  ele%vec0 = a_end%vec - matmul(ele%mat6, a_start%vec)
  if (present (end)) end = a_end

end subroutine

