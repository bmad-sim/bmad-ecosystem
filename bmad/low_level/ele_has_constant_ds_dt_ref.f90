!+
! Function ele_has_constant_ds_dt_ref (ele) result (is_const)
!
! Function to determine if an element has a constant longitudinal reference velocity.
! When in doubt, the assumption is that the longitudinal velocity is not constant.
! EG: Cavity with zero field is marked as non-constant ds_dt(ref).
!
! Input:
!   ele -- ele_struct: Element.
!
! Output:
!   is_const -- Logical: True if reference velocity must be a constant.
!-

recursive function ele_has_constant_ds_dt_ref (ele) result (is_const)

use bmad_routine_interface, except_dummy => ele_has_constant_ds_dt_ref

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer ie
logical is_const

! Anything with longitudinal electric fields or anything
! where the "zero-orbit" is not a straight line down the middle
! has a varying ds/dt(ref).

select case (ele%key)
case (lcavity$, custom$, hybrid$, wiggler$, undulator$, rfcavity$, e_gun$)
  is_const = .false.

case (em_field$)
  if (ele%slave_status == slice_slave$) then
    is_const = ele_has_constant_ds_dt_ref(ele%lord)

  elseif (ele%slave_status == super_slave$) then
    do ie = 1, ele%n_lord
      lord => pointer_to_lord(ele, ie)
      is_const = ele_has_constant_ds_dt_ref(lord)
      if (.not. is_const) return
    enddo

  else
    is_const = is_true(ele%value(constant_ref_energy$))
  endif

case default
  is_const = .true.
end select

end function ele_has_constant_ds_dt_ref

