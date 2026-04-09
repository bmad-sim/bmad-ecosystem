!+
! Function ele_has_constant_ds_dt_ref (ele) result (is_const)
!
! Function to determine if an element has a constant longitudinal reference velocity.
! When in doubt, the assumption is that the longitudinal velocity is not constant.
! For example, a cavity with zero field is marked as non-constant ref ds/dt.
!
! Notable examples:
!  *) A wiggler has non-constant ref ds/dt.
!  *) A super_lord element that normally has a const ds/dt (say, a solenoid), if superimposed with 
!       a non-const ds/dt element (say, a lcavity), will have a non-cost ds/dt.
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
type (ele_struct), pointer :: lord, slave
integer ie
logical is_const

! For a super_lord or multipass_lord check the slaves.

if (ele%lord_status == multipass_lord$) then
  slave => pointer_to_slave(ele, 1)
  is_const = ele_has_constant_ds_dt_ref(slave)
  return
endif

! A super_lord element that normally has a const ds/dt (say, a solenoid), if superimposed with 
! a non-const ds/dt element (say, a lcavity), will have a non-cost ds/dt.

if (ele%lord_status == super_lord$) then

  is_const = .true.
  do ie = 1, ele%n_slave
    slave => pointer_to_slave(ele, 1)

    select case (slave%key)
    case (lcavity$, custom$, hybrid$, wiggler$, undulator$, rfcavity$, e_gun$)
      is_const = .false.

    case (em_field$)
      is_const = (is_const .and. is_true(slave%value(constant_ref_energy$)))
    end select
  enddo

  return
endif

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

