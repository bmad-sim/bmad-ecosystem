!+
! Subroutine check_attrib_free (ele, ix_attrib, ring, free, err_print_flag)
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

subroutine check_attrib_free (ele, ix_attrib, ring, free, err_print_flag)

  use bmad_struct
  use bmad_interface, except => check_attrib_free

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele

  integer ix_attrib, i, ir, ix

  logical free, do_print
  logical, optional :: err_print_flag

! init

  free = .true.

  do_print = .true.
  if (present(err_print_flag)) do_print = err_print_flag

! super_slaves attributes cannot be varied

  if (ele%control_type == super_slave$) then
    if (do_print) then
      print *, 'ERROR IN CHECK_ATTRIB_FREE:'
      print *, '      TRYING TO VARY AN ATTRIBUTE OF: ', ele%name
      print *, '      WHICH IS A SUPER_SLAVE WILL NOT WORK.'
      print *, '      VARY THE ATTRIBUTE OF ONE OF ITS SUPER_LORDS INSTEAD.'
    endif
    return
  endif

! if the attribute is controled by an overlay lord then it cannot be varied

  do i = ele%ic1_lord, ele%ic2_lord
    ix = ring%ic_(i)
    ir = ring%control_(ix)%ix_lord
    if (ring%ele_(ir)%control_type == overlay_lord$) then
      if (ring%control_(ix)%ix_attrib == ix_attrib) then
        if (do_print) print '((1x, a))', &
            'ERROR IN CHECK_ATTRIB_FREE. THE ATTRIBUTE: ' // &
                                             attribute_name(ele, ix_attrib), &
            '      OF ELEMENT: ' // ele%name, &
            '      IS CONTROLLED BY OVERLAY_LORD: ' // ring%ele_(ir)%name, &
            '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
        return
      endif
    endif
  enddo

! only one particular attribute of an overlay lord is allowed to be adjusted

  if (ele%control_type == overlay_lord$) then
    if (ix_attrib /= ele%ix_value) then
      if (do_print) print '((1x, a))', &
              'ERROR IN CHECK_ATTRIB_FREE:' // &
                        ' OVERLAYS HAVE ONLY ONE ATTRIBUTE TO VARY.', &
              '      FOR THE OVERLAY: ' // trim(ele%name), &
              '      THAT ATTRIBUTE IS: ' // ele%attribute_name
      return
    endif
  endif

! Everything OK so far

  free = .false.

! check if it is a dependent variable.

  select case (ele%key)
  case (sbend$)
    if (any(ix_attrib == (/ angle$, l_chord$, rho$ /))) free = .true.
  case (rfcavity$)
    if (ix_attrib == rf_frequency$ .and. ele%value(harmon$) /= 0) free = .true.
  case (beambeam$)
    if (ix_attrib == bbi_const$) free = .true.
  case (wiggler$)
    if (ix_attrib == k1$ .or. ix_attrib == rho$) free = .true. 
  case (lcavity$)
    if (ix_attrib == delta_e$) free = .true.
  case (elseparator$)
    if (ix_attrib == e_field$ .or. ix_attrib == voltage$) free = .true.
  end select

  if (free .and. do_print) then
    print '((1x, a))', &
            'ERROR IN CHECK_ATTRIB_FREE. THE ATTRIBUTE: ' // &
                                             attribute_name(ele, ix_attrib), &
            '      OF ELEMENT: ' // ele%name, &
            '      IS A DEPENDENT VARIABLE.', &
            '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
  endif

! field_master on means that the b_field and b_gradient values control
! the strength.

  if (ele%field_master) then
    select case (ele%key)
    case (quadrupole$)
      if (ix_attrib == k1$) free = .true.
    case (sextupole$)
      if (ix_attrib == k2$) free = .true.
    case (octupole$)
      if (ix_attrib == k3$) free = .true.
    case (solenoid$)
      if (ix_attrib == ks$) free = .true.
    case (sbend$)
      if (ix_attrib == g$) free = .true.
    end select
  endif

  if (free .and. do_print) then
    print '((1x, a))', &
            'ERROR IN CHECK_ATTRIB_FREE. THE ATTRIBUTE: ' // &
                                             attribute_name(ele, ix_attrib), &
            '      OF ELEMENT: ' // ele%name, &
            '      IS A DEPENDENT VARIABLE SINCE FIELD_MASTER IS ON.', &
            '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
  endif


end subroutine
