!+
! Function attribute_free (ele, ix_attrib, ring, err_print_flag) result (free)
!
! Subroutine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly are super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay_lord.
!
! Also dependent variables such as the angle of a bend cannot be 
!   freely variable.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: Element
!   ix_attrib       -- Integer: Index to the attribute in ele%value() array.
!   ring            -- Ring_struct: Ring structure.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message if attribute is not free.
!
! Output:
!   free   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"


function attribute_free (ele, ix_attrib, ring, err_print_flag) result (free)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele

  integer ix_attrib, i, ir, ix

  logical free, do_print
  logical, optional :: err_print_flag

! init

  free = .false.
  do_print = logic_option (err_print_flag, .true.)

! if the attribute is controled by an overlay lord then it cannot be varied

  do i = ele%ic1_lord, ele%ic2_lord
    ix = ring%ic_(i)
    ir = ring%control_(ix)%ix_lord
    if (ring%ele_(ir)%control_type == overlay_lord$) then
      if (ring%control_(ix)%ix_attrib == ix_attrib) then
        call print_error (ele, &
            'THE ATTRIBUTE: ' // attribute_name(ele, ix_attrib), &
            'OF ELEMENT: ' // ele%name, &
            'IS CONTROLLED BY OVERLAY_LORD: ' // ring%ele_(ir)%name)
        return
      endif
    endif
  enddo

! further checking

  call check_this_attribute_free (ele, ix_attrib)

!-------------------------------------------------------------------------
contains

recursive subroutine check_this_attribute_free (ele2, ix_attrib)

  type (ele_struct) ele2
  integer ix_attrib, i

! super_slaves attributes cannot be varied

  free = .false.

  if (ele2%control_type == super_slave$) then
    call print_error (ele2, &
            'TRYING TO VARY AN ATTRIBUTE OF: ' // ele2%name, &
            'WHICH IS A SUPER_SLAVE WILL NOT WORK.', &
            'VARY THE ATTRIBUTE OF ONE OF ITS SUPER_LORDS INSTEAD.')
    return
  endif

! only one particular attribute of an overlay lord is allowed to be adjusted

  if (ele2%control_type == overlay_lord$) then
    if (ix_attrib /= ele2%ix_value) then
      call print_error (ele2, &
                     'OVERLAYS HAVE ONLY ONE ATTRIBUTE TO VARY.', &
                     'FOR THE OVERLAY: ' // trim(ele2%name), &
                     'THAT ATTRIBUTE IS: ' // ele2%attribute_name)
      return
    endif
  endif

  if (ele2%control_type == group_lord$) then
    if (ix_attrib /= command$ .and. ix_attrib /= old_command$) then
      call print_error (ele2, &
                     'GROUPS CAN ONLY VARY THE COMMAND OR OLD_COMMAND ATTRIBUTE.', &
                     'FOR THE GROUP: ' // trim(ele2%name))
      return
    endif
  endif

! Everything OK so far

  free = .true.

! check if it is a dependent variable.

  select case (ele2%key)
  case (sbend$)
    if (any(ix_attrib == (/ angle$, l_chord$, rho$ /))) free = .false.
  case (rfcavity$)
    if (ix_attrib == rf_frequency$ .and. ele2%value(harmon$) /= 0) free = .false.
  case (beambeam$)
    if (ix_attrib == bbi_const$) free = .false.
  case (wiggler$)
    if (ix_attrib == k1$ .or. ix_attrib == rho$) free = .false. 
  case (lcavity$)
    if (ix_attrib == delta_e$) free = .false.
  case (elseparator$)
    if (ix_attrib == e_field$ .or. ix_attrib == voltage$) free = .false.
  end select

  if (.not.free) then
    call print_error (ele2, &
                   'THE ATTRIBUTE: ' // attribute_name(ele2, ix_attrib), &
                   'OF ELEMENT: ' // ele2%name, &
                   'IS A DEPENDENT VARIABLE.')
    return
  endif

! field_master on means that the b_field and b_gradient values control
! the strength.

  if (ele2%field_master) then
    select case (ele2%key)
    case (quadrupole$)
      if (ix_attrib == k1$) free = .false.
    case (sextupole$)
      if (ix_attrib == k2$) free = .false.
    case (octupole$)
      if (ix_attrib == k3$) free = .false.
    case (solenoid$)
      if (ix_attrib == ks$) free = .false.
    case (sol_quad$)
      if (ix_attrib == ks$) free = .false.
      if (ix_attrib == k1$) free = .false.
    case (sbend$)
      if (ix_attrib == g$) free = .false.
    case (hkicker$, vkicker$)
      if (ix_attrib == kick$) free = .false.
    end select

    if (ix_attrib == hkick$) free = .false.
    if (ix_attrib == vkick$) free = .false.

  else
    select case (ele2%key)
    case (quadrupole$)
      if (ix_attrib == b_gradient$) free = .false.
    case (sextupole$)
      if (ix_attrib == b_gradient$) free = .false.
    case (octupole$)
      if (ix_attrib == b_gradient$) free = .false.
    case (solenoid$)
      if (ix_attrib == b_field$) free = .false.
    case (sol_quad$)
      if (ix_attrib == b_field$) free = .false.
      if (ix_attrib == b_gradient$) free = .false.
    case (sbend$)
      if (ix_attrib == b_field$) free = .false.
    case (hkicker$, vkicker$)
      if (ix_attrib == bl_kick$) free = .false.
    end select

    if (ix_attrib == bl_hkick$) free = .false.
    if (ix_attrib == bl_vkick$) free = .false.

  endif

  if (.not. free) then
    call print_error (ele2, &
             'THE ATTRIBUTE: ' // attribute_name(ele2, ix_attrib), &
             'OF ELEMENT: ' // ele2%name, &
             'IS A DEPENDENT VARIABLE SINCE FIELD_MASTER IS ' // &
                                             on_off_logic (ele2%field_master))
    return
  endif

! check slaves

  if (ele2%control_type == group_lord$ .or. ele2%control_type == overlay_lord$) then
    do i = ele2%ix1_slave, ele2%ix2_slave
      call check_this_attribute_free (ring%ele_(ring%control_(i)%ix_slave), &
                                                       ring%control_(i)%ix_attrib)
      if (.not. free) return
    enddo
  endif

end subroutine

!-------------------------------------------------------
! contains

subroutine print_error (ele2, l1, l2, l3)

  type (ele_struct) ele2
  character(*) l1, l2
  character(*), optional :: l3

!

  if (.not. do_print) return

  print *, 'ERROR IN ATTRIBUTE_FREE:'
  if (ele2%name /= ele%name) then
    print *, '      THE GROUP/OVERLAY LORD: ', ele%name    
  endif
  print '(5x, a)', trim(l1)
  print '(5x, a)', trim(l2)
  if (present(l3)) print '(5x, a)', trim(l3)

end subroutine

!-------------------------------------------------------
! contains

function on_off_logic (logic) result (name)

  logical logic
  character(3) name

!

  if (logic) then
    name = 'ON'
  else
    name = 'OFF'
  endif

end function

end function
