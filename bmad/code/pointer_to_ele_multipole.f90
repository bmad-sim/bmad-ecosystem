!+
! Subroutine pointer_to_ele_multipole (ele, a_pole, b_pole, ksl_pole, pole_type)
!
! Routine to point to the appropriate magnetic or electric poles in an element.
!
! Input:
!   ele          -- Ele_struct: Lattice element.
!   pole_type    -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!
! Output:
!   a_pole(:)    -- real(rp), pointer: Pointer to skew electric or magnetic poles. KL for multipole elements.
!   b_pole(:)    -- real(rp), pointer: Pointer to normal electric or magnetic poles. Tilt for multipole elements.
!   ksl_pole(:)  -- real(rp), pointer: For multipole elements only.
!-

subroutine pointer_to_ele_multipole (ele, a_pole, b_pole, ksl_pole, pole_type)

use bmad_struct

implicit none

type (ele_struct), target :: ele

real(rp), pointer :: a_pole(:), b_pole(:), ksl_pole(:)
integer, optional :: pole_type

!

ksl_pole => null()

select case (integer_option(magnetic$, pole_type))
case (magnetic$)
  a_pole   => ele%a_pole
  b_pole   => ele%b_pole
  if (ele%key == multipole$) ksl_pole => ele%a_pole_elec

case (electric$)
  if (ele%key == multipole$) then
    a_pole => null()   ! Multipoles do not have electric fields.
    b_pole => null()
  else
    a_pole => ele%a_pole_elec
    b_pole => ele%b_pole_elec
  endif
end select

end subroutine pointer_to_ele_multipole

