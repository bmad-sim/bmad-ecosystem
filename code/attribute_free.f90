!+
! Function attribute_free (ix_ele, ix_attrib, lat, err_print_flag) result (free)
!
! Subroutine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly include super_slave attributes (since
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
!   ix_ele          -- Integer: Index of element in lat%ele_(:) array.
!   ix_attrib       -- Integer: Index to the attribute in ele%value(:) array.
!   lat             -- Ring_struct: Lattice structure.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message if attribute is not free.
!
! Output:
!   free   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"


function attribute_free (ix_ele, ix_attrib, lat, err_print_flag) result (free)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target :: lat

  integer ix_attrib, ix_ele, ix_ele0, ix_attrib0

  logical free, do_print
  logical, optional :: err_print_flag
  character(16) :: r_name = 'attribute_free'

! init & check

  do_print = logic_option (err_print_flag, .true.)
  ix_ele0 = ix_ele
  ix_attrib0 = ix_attrib
  call check_this_attribute_free (ix_ele, ix_attrib)

!-------------------------------------------------------------------------
contains

recursive subroutine check_this_attribute_free (ix_ele, ix_attrib, ix_lord)

  type (ele_struct), pointer :: ele
  integer ix_ele, ix_attrib, i, ir, ix
  integer, optional :: ix_lord

! super_slaves attributes cannot be varied

  free = .false.
  ele => lat%ele_(ix_ele)

! if the attribute is controled by an overlay lord then it cannot be varied

  do i = ele%ic1_lord, ele%ic2_lord
    ix = lat%ic_(i)
    ir = lat%control_(ix)%ix_lord
    if (present(ix_lord)) then
      if (ix_lord == ir) cycle
    endif
    if (lat%ele_(ir)%control_type == overlay_lord$) then
      if (lat%control_(ix)%ix_attrib == ix_attrib) then
        call print_error (ix_ele, ix_attrib, &
            'IT IS CONTROLLED BY THE OVERLAY_LORD: ' // lat%ele_(ir)%name)
        return
      endif
    endif
  enddo

  if (ele%control_type == super_slave$) then
    call print_error (ix_ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.')
    return
  endif

! only one particular attribute of an overlay lord is allowed to be adjusted

  if (ele%control_type == overlay_lord$) then
    if (ix_attrib /= ele%ix_value) then
      call print_error (ix_ele, ix_attrib, &
             'FOR THIS OVERLAY ELEMENT THE ATTRIBUTE TO VARY IS: ' // ele%attribute_name)
      return
    endif
  endif

  if (ele%control_type == group_lord$) then
    if (ix_attrib /= command$ .and. ix_attrib /= old_command$) then
      call print_error (ix_ele, ix_attrib, &
            'FOR THIS GROUP ELEMENT THE ATTRIBUTE TO VARY IS: "COMMAND" OR "OLD_COMMAND"')
      return
    endif
  endif

! Everything OK so far

  free = .true.

! check if it is a dependent variable.

  select case (ele%key)
  case (sbend$)
    if (any(ix_attrib == (/ angle$, l_chord$, rho$ /))) free = .false.
  case (rfcavity$)
    if (ix_attrib == rf_frequency$ .and. ele%value(harmon$) /= 0) free = .false.
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
    call print_error (ix_ele, ix_attrib, 'THE ATTRIBUTE IS A DEPENDENT VARIABLE.')
    return
  endif

! field_master on means that the b_field and b_gradient values control
! the strength.

  if (ele%field_master) then
    select case (ele%key)
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
    select case (ele%key)
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
    call print_error (ix_ele, ix_attrib, &
         'THE ATTRIBUTE IS A DEPENDENT VARIABLE SINCE THE FIELD_MASTER ATTRIBUTE IS ' // &
                                             on_off_logic (ele%field_master))
    return
  endif

! check slaves

  if (ele%control_type == group_lord$ .or. ele%control_type == overlay_lord$) then
    do i = ele%ix1_slave, ele%ix2_slave
      call check_this_attribute_free (lat%control_(i)%ix_slave, &
                                               lat%control_(i)%ix_attrib, ix_ele)
      if (.not. free) return
    enddo
  endif

end subroutine

!-------------------------------------------------------
! contains

subroutine print_error (ix_ele, ix_attrib, l1)

  integer ix_ele, ix_attrib
  character(*) l1
  character(80) li(8)

!

  if (.not. do_print) return

  li(1) = 'ERROR IN ATTRIBUTE_FREE:'
  li(2) = '     THE ATTRIBUTE: ' // attribute_name(lat%ele_(ix_ele0), ix_attrib0)
  li(3) = '     OF THE ELEMENT: ' // lat%ele_(ix_ele0)%name

  if (ix_ele == ix_ele0) then
    li(4) = '   IS NOT FREE TO VARY SINCE:'
    li(5) = '     ' // l1
    call out_io (s_error$, r_name, li(1:5))   
  else 
    li(4) = '   IS NOT FREE TO VARY SINCE IT IS TRYING TO CONTROL:'
    li(5) = '     THE ATTRIBUTE: ' // attribute_name(lat%ele_(ix_ele), ix_attrib)
    li(6) = '     OF THE ELEMENT: ' // lat%ele_(ix_ele)%name
    li(7) = '   AND THIS IS NOT FREE TO VARY SINCE:'
    li(8) = '     ' // l1
    call out_io (s_error$, r_name, li)        
  endif

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
