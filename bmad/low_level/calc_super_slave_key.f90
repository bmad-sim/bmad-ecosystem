!+
! Subroutine calc_super_slave_key (lord1, lord2, slave, create_jumbo_slave)
!
! Function to decide what the type of element a super_slave should be based
! upon the types of its lords.
!
! Input:
!   lord1     -- Ele_struct: First slave.
!     %key
!     %sub_key
!   lord2     -- Ele_struct: Second slave.
!     %key
!     %sub_key
!   create_jumbo_slave
!              -- Logical, optional: If True then slave%key will be set to em_field.
!                   Default is False. 
!
! Output:
!   slave -- Ele_struct: Super_slave element.
!     %key        -- Set to -1 if there is an error.
!     %sub_key
!-

subroutine calc_super_slave_key (lord1, lord2, slave, create_jumbo_slave)

use bmad_struct

implicit none

type (ele_struct), target :: lord1, lord2, slave
integer key1, key2
logical, optional :: create_jumbo_slave

!

key1 = lord1%key
key2 = lord2%key

slave%key = -1  ! Default if no superimpose possible
slave%sub_key = 0

! If one element is a drift then slave%key = key of other element.
! Control elements, naturally zero length elements, etc. cannot be superimposed upon

select case (key1)
case (overlay$, group$, ramper$, girder$, taylor$, match$, patch$, fiducial$, floor_shift$, multipole$, ab_multipole$); return
end select

select case (key2)
case (overlay$, group$, ramper$, girder$, taylor$, match$, patch$, fiducial$, floor_shift$, multipole$, ab_multipole$); return
end select

select case (key1)
case (drift$, pipe$)
  slave%key = key2
  if (key1 == pipe$ .and. key2 == drift$) slave%key = pipe$
  slave%sub_key = lord2%sub_key
  return
end select

select case (key2)
case (drift$, pipe$)
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
end select

if (key1 == sbend$ .or. key1 == rf_bend$ .or. key2 == sbend$ .or. key2 == rf_bend$) return

! If there are misalignments then no superposition is possible

if (lord1%value(x_offset$) /= 0 .or. lord1%value(y_offset$) /= 0 .or. lord1%value(z_offset$) /= 0 .or. &
    lord1%value(x_pitch$) /= 0 .or. lord1%value(y_pitch$) /= 0 .or. &
    lord2%value(x_offset$) /= 0 .or. lord2%value(y_offset$) /= 0 .or. lord2%value(z_offset$) /= 0 .or. &
    lord2%value(x_pitch$) /= 0 .or. lord2%value(y_pitch$) /= 0) then
  slave%key = em_field$
  if (key1 == lcavity$ .or. key2 == lcavity$) then
    slave%value(constant_ref_energy$) = false$
  elseif (key1 == em_field$) then
    slave%value(constant_ref_energy$) = lord1%value(constant_ref_energy$)
  elseif (key2 == em_field$) then
    slave%value(constant_ref_energy$) = lord2%value(constant_ref_energy$)
  else
    slave%value(constant_ref_energy$) = true$
  endif
  return
endif

! Superimposing two of like kind...

if (key1 == key2) then
  select case (key1)
  case (rfcavity$, wiggler$, undulator$)
    slave%key = em_field$
    slave%value(constant_ref_energy$) = true$
  case (lcavity$)
    slave%key = em_field$
    slave%value(constant_ref_energy$) = false$
  case (em_field$)
    slave%key = em_field$
    if (is_false(lord1%value(constant_ref_energy$)) .or. is_false(lord2%value(constant_ref_energy$))) then
      slave%value(constant_ref_energy$) = false$
    else
      slave%value(constant_ref_energy$) = true$
    endif
  case default
    slave%key = key1
  end select
  return
endif

! If one element is a rcollimator, monitor, or instrument then slave%key = key of other element.

if (lord1%aperture_type == elliptical$ .or. lord2%aperture_type == elliptical$ ) slave%aperture_type = elliptical$

if (any(key1 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  slave%key = key2
  slave%sub_key = lord2%sub_key
  return
endif

if (any(key2 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! If one element is a kicker then slave%key = key of other element.

if (any(key1 == [kicker$, hkicker$, vkicker$])) then
  if (any(key2 == [kicker$, hkicker$, vkicker$])) then
    slave%key = kicker$
  else
    slave%key = key2
  endif
  return
endif

if (any(key2 == [kicker$, hkicker$, vkicker$])) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! General case...

! em_field wanted

if (logic_option(.false., create_jumbo_slave)) then
  slave%key = em_field$
  if (key1 == lcavity$ .or. key2 == lcavity$) then
    slave%value(constant_ref_energy$) = false$
  elseif (key1 == em_field$) then
    slave%value(constant_ref_energy$) = lord1%value(constant_ref_energy$)
  elseif (key2 == em_field$) then
    slave%value(constant_ref_energy$) = lord2%value(constant_ref_energy$)
  else
    slave%value(constant_ref_energy$) = true$
  endif
  return
endif

!

select case (key1)
case (quadrupole$,  solenoid$, sol_quad$) 
  select case (key2)
  case (quadrupole$);    slave%key = sol_quad$
  case (solenoid$);      slave%key = sol_quad$
  case (sol_quad$);      slave%key = sol_quad$
  end select
end select

if (slave%key /= -1) return  ! Have found something

! Only thing left is to use em_field type element.

slave%key = em_field$
if (key1 == lcavity$ .or. key2 == lcavity$) then
  slave%value(constant_ref_energy$) = false$
elseif (key1 == em_field$) then
  slave%value(constant_ref_energy$) = lord1%value(constant_ref_energy$)
elseif (key2 == em_field$) then
  slave%value(constant_ref_energy$) = lord2%value(constant_ref_energy$)
else
  slave%value(constant_ref_energy$) = true$
endif

end subroutine calc_super_slave_key

