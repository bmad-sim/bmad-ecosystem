!+
! Subroutine type2_ele (ele, type_zero_attrib, type_mat6, twiss_out, 
!                                           type_control, lines, n_lines, ring)
!
! Subroutine to put information on an element in a string array. 
! See also the subroutine: type_ele.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele          -- Ele_struct: Element
!   type_zero_attrib  
!                -- Logical: If true then type all attributes even if the
!                   attribute value is 0 (except for multipoles).
!   type_mat6    -- Integer:
!                      type_mat6 = 0   => Do not type ELE%MAT6
!                      type_mat6 = 4   => Type 4X4 XY submatrix
!                      type_mat6 = 6   => Type full 6x6 matrix
!   twiss_out    -- Integer: Type the Twiss parameters at
!                   the end of the element?
!                       = 0  => Do not type the Twiss parameters
!                       = radians$  => Type Twiss, use radians for phi.
!                       = degrees$  => Type Twiss, use degrees for phi.
!                       = cycles$   => Type Twiss, use cycles (1 = 2pi) units.
!   type_control -- Logical: If true then type control status.
!   ring         -- Ring_struct, optional: Needed for control typeout.
!              
! Output       
!   lines(:)     -- Character*(*): Character array to hold the output.
!   n_lines      -- Number of lines used
!-

!$Id$
!$Log$
!Revision 1.6  2002/01/23 20:08:29  dcs
!Bug fix
!
!Revision 1.5  2002/01/08 21:44:44  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.4  2001/11/29 19:39:54  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.3  2001/10/12 20:53:35  rwh24
!DCS changes and two files added
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine type2_ele (ele, type_zero_attrib, type_mat6, twiss_out,  &
                                          type_control, lines, n_lines, ring)

  use bmad_struct
  use bmad_interface         

  implicit none

  type (ele_struct)  ele
  type (ring_struct), optional :: ring

  integer i, j, n, twiss_out, type_mat6, ix, iv, ic, ct, n_lines
  integer nl, nl2, i_max, particle

  real coef, value(n_attrib_maxx), value2(n_attrib_maxx)
  real a(0:n_pole_maxx), b(0:n_pole_maxx)
  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

  character*16 a_name, name, null_type / ' ' /
  character*(*) lines(:)
  character val_str*12, str*80

  logical type_control, type_zero_attrib

! Encode element name and type

  ct = ele%control_type
  
  write (lines(1), *) 'Element Name: ', ele%name
  nl = 1

  if (ele%type /= null_type) then
    nl = nl + 1
    write (lines(nl), *) 'Element Type: "', ele%type, '"'
  endif

! Encode element key and attributes

  if (ele%key <= 0) then
    nl = nl + 1
    write (lines(nl), *) 'Key: UNKNOWN!', ele%key

  else

    nl = nl + 1
    write (lines(nl), *) 'Key: ', key_name(ele%key)

    if (.not. type_zero_attrib) then
      nl = nl + 1
      write (lines(nl), *) 'Attribute values [Only non-zero values shown]:'
    endif

    if (ct == overlay_lord$) then
      i = ele%ix_value
      name = ele%attribute_name
      nl = nl + 1
      write (lines(nl), '(i6, 3x, 2a, 1pe12.4)') i, name, ' =', ele%value(i)
    else
      if (ele%key == multipole$) then
        i_max = n_attrib_maxx
      else
        i_max = ix1_m$ - 1 
      endif  
      do i = 1, i_max
        if (i == energy$) cycle   ! energy not yet implemented.
        if (i == aperture$) cycle ! aperture not yet implemented
        if (attribute_name(ele, i) /= null_name) then
          if (ele%value(i) /= 0 .or. type_zero_attrib) then
            nl = nl + 1
            write (lines(nl), '(i6, 3x, 2a, 1pe12.4)')  i, &
                         attribute_name(ele, i), ' =', ele%value(i)
          endif
        endif
      enddo

      if (ele%key /= multipole$) then
        particle = +1
        if (present(ring)) particle = ring%param%particle
        call multipole_ele_to_ab (ele, particle, a, b, .false.)
        call multipole_ab_to_value (a, b, value)
        call multipole_ele_to_kt (ele, particle, knl, tilt, .true.)
        call multipole_kt_to_ab (knl, tilt, a, b)
        call multipole_ab_to_value (a, b, value2)
        do i = ix1_m$, ix2_m$
          if (ele%value(i) /= 0 .or. value(i) /= 0 .or. value2(i) /= 0) then
            nl = nl + 1
            if (ele%key == ab_multipole$) then
              write (lines(nl), '(5x, a, 3(a, 1pe11.3))') &
                   attribute_name(ele, i), ' =', ele%value(i), &
                   '   W/Tilt:', value2(i)
            else
              write (lines(nl), '(5x, a, 3(a, 1pe11.3))') &
                   attribute_name(ele, i), ' =', ele%value(i), &
                   '   Scaled:', value(i), '   W/Tilt:', value2(i)
            endif
          endif
        enddo
      endif

    endif
  endif

! Encode on/off status and s_position

  if (.not. ele%is_on) then
    nl = nl + 1
    write (lines(nl), *) '*** Note: Element is turned OFF ***'
  endif

  if (.not. ele%multipoles_on) then
    nl = nl + 1
    write (lines(nl), *) '*** Note: Element Multipoles are turned OFF ***'
  endif

  nl = nl + 1
  write (lines(nl), '(1x, a, f13.4)') 'S:', ele%s

! Encode methods

  write (lines(nl+1), *)
  write (lines(nl+2), '(2a)') ' Tracking_method:  ', &
                                    calc_method_name(ele%tracking_method)
  write (lines(nl+3), '(2a)') ' Mat6_calc_method: ', &
                                    calc_method_name(ele%mat6_calc_method)
  write (lines(nl+4), '(2a)') ' Num_steps: ', ele%num_steps 
  nl = nl + 4

! Encode lord info

  if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
    write (lines(nl+1), *)
    write (lines(nl+2), *) &
                   'ERROR: N_SLAVE DOES NOT MATCH IX2_SLAVE-IX1_SLAVE+1 !!'
    write (lines(nl+3), *) ele%n_slave, ele%ix2_slave, ele%ix1_slave
    nl = nl + 3
  endif          

  if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
    write (lines(nl+1), *)
    write (lines(nl+2), *) 'ERROR: N_LORD DOES NOT MATCH Ic1_LORD-Ic2_LORD+1 !!'
    write (lines(nl+3), *) ele%n_lord, ele%ic1_lord, ele%ic2_lord
    nl = nl + 3
  endif

! Encode slave info.
! For super_lords there is no attribute_name associated with a slave.
! For slaves who are overlay_lords then the attribute_name is obtained by
!   looking at the overlay_lord's 1st slave (slave of slave of the input ele).

  if (type_control) then

    nl = nl + 1
    write (lines(nl), *)

    if (ct <= 0) then
      nl = nl + 1
      write (lines(nl), *) 'Control_type: UNKNOWN!', ct
    else
      nl = nl + 1
      write (lines(nl), *) 'Control_type: ', control_name(ct)

      if (ele%n_slave /= 0) then
        write (lines(nl+1), '(1x, a, i4)') 'Slaves: Number:', ele%n_slave
        write (lines(nl+2), *) &
          '    Name            Ring_index  Attribute       Coefficient'
        nl = nl + 2
        do i = ele%ix1_slave, ele%ix2_slave
          j = ring%control_(i)%ix_slave
          iv = ring%control_(i)%ix_attrib
          coef = ring%control_(i)%coef
          if (ct == super_lord$) then
            a_name = '--------'
          elseif (ring%ele_(j)%control_type == overlay_lord$) then
            if (iv == ring%ele_(j)%ix_value) then
              ix = ring%control_(ring%ele_(j)%ix1_slave)%ix_slave
              a_name = attribute_name(ring%ele_(ix), iv)
            else
              a_name = '** BAD POINTER! **'
            endif            
          else
            a_name = attribute_name(ring%ele_(j), iv)
          endif
          nl = nl + 1
          write (lines(nl), '(5x, a, i10, 2x, a16, 1pe11.3, 1pe12.3)') &
                                ring%ele_(j)%name, j, a_name, coef
        enddo
      endif

      if (ele%n_lord /= 0) then
        write (lines(nl+1), '(1x, a, i4)') 'Lords: Number:', ele%n_lord
        write (lines(nl+2), *) &
  '    Name            Ring_index  Attribute       Coefficient       Value'
        nl = nl + 2
        do i = ele%ic1_lord, ele%ic2_lord
          ic = ring%ic_(i)
          j = ring%control_(ic)%ix_lord
          coef = ring%control_(ic)%coef
          ct = ele%control_type
          if (ct == super_slave$) then
            a_name = '--------'
            val_str = '    --------'
          else
            iv = ring%control_(ic)%ix_attrib
            a_name = attribute_name(ele, iv)
            ix = ring%ele_(j)%ix_value
            write (val_str, '(1pe12.3)') ring%ele_(j)%value(ix)
          endif
          nl = nl + 1  
          write (lines(nl), '(5x, a, i10, 2x, a16, 1pe11.3, a12)') &
                             ring%ele_(j)%name, j, a_name, coef, val_str
        enddo
      endif

    endif
  endif

! Encode Twiss info

  if (twiss_out > 0) then
    call type2_twiss (ele, twiss_out, lines(nl+1:), nl2)
    nl = nl + nl2
  endif

! Encode mat6 info

  n = type_mat6

  if (n /= 0) then
    nl = nl + 1
    write (lines(nl), *)
  endif

  if (any(abs(ele%mat6(1:n,1:n)) >= 1000)) then
    do i = 1, n
      nl = nl + 1
      write (lines(nl), '(1p6e11.3)') (ele%mat6(i, j), j = 1, n)
    enddo
  else
    do i = 1, n
      nl = nl + 1
      write (lines(nl), '(6f10.5)') (ele%mat6(i, j), j = 1, n)
    enddo
  endif

  n_lines = nl

end subroutine
