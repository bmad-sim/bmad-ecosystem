!+
! Subroutine type2_ele (ele, type_zero_attrib, type_mat6, type_taylor, 
!                               twiss_out, type_control, lines, n_lines, ring)
!
! Subroutine to put information on an element in a string array. 
! See also the subroutine: type_ele.
!
! Modules needed:
!   use bmad
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
!   type_taylor  -- Logical: Type out taylor series?
!                   if ele%taylor is not allocated then this is ignored.
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
!   lines(:)     -- Character*80, pointer: Character array to hold the 
!                     output. The array size of lines(:) will be set by
!                     this subroutine. Note: You need to deallocate lines
!                     between each call to type2_ele.
!   n_lines      -- Number of lines in lines(:).
!-

!$Id$
!$Log$
!Revision 1.13  2003/01/27 14:40:47  dcs
!bmad_version = 56
!
!Revision 1.11  2002/08/23 20:19:02  dcs
!Update Documentation.
!
!Revision 1.10  2002/07/31 14:33:51  dcs
!Extended lines(:) array.
!
!Revision 1.9  2002/06/13 14:54:31  dcs
!Interfaced with FPP/PTC
!
!Revision 1.8  2002/02/23 20:32:29  dcs
!Double/Single Real toggle added
!
!Revision 1.7  2002/01/24 16:37:35  dcs
!Bug fix
!
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

subroutine type2_ele (ele, type_zero_attrib, type_mat6, type_taylor, &
                               twiss_out, type_control, lines, n_lines, ring)

  use bmad_struct
  use bmad_interface
  use multipole_mod

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (ring_struct), optional, intent(in) :: ring
  type (wig_term_struct), pointer :: term

  integer, intent(in) :: type_mat6
  integer, intent(out) :: n_lines
  integer i, j, k, n, twiss_out, ix, iv, ic, ct, nl2
  integer nl, nt, n_max, i_max, particle

  real(rdef) coef
  real(rdef) a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rdef) a2(0:n_pole_maxx), b2(0:n_pole_maxx)

  character(80), pointer :: lines(:)
  character(80), pointer :: li(:), li2(:)
  character(16) a_name, name
  character(12) val_str
  character(80) str
  character(2) str_i

  logical, intent(in) :: type_taylor
  logical, intent(in) :: type_control, type_zero_attrib

! init

  allocate (li(300))
  
! Encode element name and type

  ct = ele%control_type
  
  write (li(1), *) 'Element Name: ', ele%name
  nl = 1

  if (ele%type /= blank_name) then
    nl = nl + 1
    write (li(nl), *) 'Element Type: "', ele%type, '"'
  endif

  if (associated(ele%descrip)) then
    nl = nl + 1
    write (li(nl), *), 'Descrip: ', trim(ele%descrip)
  endif

! Encode element key and attributes

  if (ele%key <= 0) then
    nl = nl + 1
    write (li(nl), *) 'Key: UNKNOWN!', ele%key

  else

    nl = nl + 1
    write (li(nl), *) 'Key: ', key_name(ele%key)

    if (ele%sub_key /= 0) then
      nl = nl + 1
      write (li(nl), *) 'Sub Key: ', sub_key_name(ele%sub_key)
    endif

    if (.not. type_zero_attrib) then
      nl = nl + 1
      write (li(nl), *) 'Attribute values [Only non-zero values shown]:'
    endif

    if (ct == overlay_lord$) then
      i = ele%ix_value
      name = ele%attribute_name
      nl = nl + 1
      write (li(nl), '(i6, 3x, 2a, 1pe12.4)') i, name, ' =', ele%value(i)

    else
      do i = 1, n_attrib_maxx
        if (i == aperture$) cycle ! aperture not yet implemented
        if (attribute_name(ele, i) /= null_name) then
          if (ele%value(i) /= 0 .or. type_zero_attrib) then
            nl = nl + 1
            write (li(nl), '(i6, 3x, 2a, 1pe12.4)')  i, &
                         attribute_name(ele, i), ' =', ele%value(i)
          endif
        endif
      enddo

      if (associated(ele%a)) then
        particle = +1
        if (present(ring)) particle = ring%param%particle

        if (ele%key == multipole$) then
          call multipole_ele_to_kt (ele, particle, a,  b,  .false.)
          call multipole_ele_to_kt (ele, particle, a2, b2, .true.)
        else
          call multipole_ele_to_ab (ele, particle, a,  b,  .false.)
          call multipole_ele_to_ab (ele, particle, a2, b2, .true.)
        endif

        do i = 0, n_pole_maxx
          if (ele%a(i) == 0 .and. ele%b(i) == 0) cycle
          write (str_i, '(i2)') i
          call string_trim (str_i, str_i, ix)
          if (ele%key == ab_multipole$) then
            write (li(nl+1), '(5x, 2a, 2(a, 1pe11.3))') &
                   'A', str_i, ' =', ele%a(i), '   W/Tilt:', a2(i)
            write (li(nl+2), '(5x, 2a, 2(a, 1pe11.3))') &
                   'B', str_i, ' =', ele%b(i), '   W/Tilt:', b2(i)
          elseif (ele%key == multipole$) then
            write (li(nl+1), '(5x, 2a, 2(a, 1pe11.3))') &
                   'K', trim(str_i), 'L =', ele%a(i), '   W/Tilt:', a2(i)
            write (li(nl+2), '(5x, 2a, 2(a, 1pe11.3))') &
                   'T', trim(str_i), '  =', ele%b(i), '   W/Tilt:', b2(i)
          else
            write (li(nl+1), '(5x, 2a, 3(a, 1pe11.3))') 'A', str_i, ' =', &
                   ele%a(i), '   Scaled:', a(i), '   W/Tilt:', a2(i)
            write (li(nl+2), '(5x, 2a, 3(a, 1pe11.3))') 'B', str_i, ' =', &
                   ele%b(i), '   Scaled:', b(i), '   W/Tilt:', b2(i)
          endif

          nl = nl + 2
        enddo

      endif

    endif
  endif

! wiggler terms

  if (ele%key == wiggler$ .and. ele%sub_key == map_type$) then
    nl = nl + 1
    write (li(nl), '(a, 6x, a, 3(9x, a), 7x, a)') ' Term#', &
                              'Coef', 'K_x', 'K_y', 'K_z', 'phi_z   Type'
    do i = 1, size(ele%wig_term)
      term => ele%wig_term(i)
      write (li(nl+i), '(i4, 5f12.6, 3x, a)') i, term%coef, &
          term%kx, term%ky, term%kz, term%phi_z, wig_term_type_name(term%type)
    enddo
    nl = nl + size(ele%wig_term)
  endif

! Encode on/off status and s_position

  if (.not. ele%is_on) then
    nl = nl + 1
    write (li(nl), *) '*** Note: Element is turned OFF ***'
  endif

  if (.not. ele%multipoles_on) then
    nl = nl + 1
    write (li(nl), *) '*** Note: Element Multipoles are turned OFF ***'
  endif

  nl = nl + 1
  write (li(nl), '(1x, a, f13.4)') 'S:', ele%s

! Encode methods

  write (li(nl+1), *)
  write (li(nl+2), '(2a)')      ' Tracking_method:  ', &
                                    calc_method_name(ele%tracking_method)
  write (li(nl+3), '(2a)')      ' Mat6_calc_method: ', &
                                    calc_method_name(ele%mat6_calc_method)
  write (li(nl+4), '(a, i4)')   ' Integration_order:', ele%integration_order 
  write (li(nl+5), '(a, i4)')   ' Num_steps:        ', ele%num_steps 
  write (li(nl+6), '(a, l)')    ' Symplectify:      ', ele%symplectify
  nl = nl + 6

  if (ele%ptc_kind /= 0) then
    write (li(nl+1), '(a, i4)') ' PTC_kind:        ', ele%ptc_kind
    nl = nl + 1
  endif
  
! Encode lord info

  if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
    write (li(nl+1), *)
    write (li(nl+2), *) &
                   'ERROR: N_SLAVE DOES NOT MATCH IX2_SLAVE-IX1_SLAVE+1 !!'
    write (li(nl+3), *) ele%n_slave, ele%ix2_slave, ele%ix1_slave
    nl = nl + 3
  endif          

  if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
    write (li(nl+1), *)
    write (li(nl+2), *) 'ERROR: N_LORD DOES NOT MATCH Ic1_LORD-Ic2_LORD+1 !!'
    write (li(nl+3), *) ele%n_lord, ele%ic1_lord, ele%ic2_lord
    nl = nl + 3
  endif

! Encode slave info.
! For super_lords there is no attribute_name associated with a slave.
! For slaves who are overlay_lords then the attribute_name is obtained by
!   looking at the overlay_lord's 1st slave (slave of slave of the input ele).

  if (type_control) then

    nl = nl + 1
    write (li(nl), *)

    if (ct <= 0) then
      nl = nl + 1
      write (li(nl), *) 'Control_type: UNKNOWN!', ct
    else
      nl = nl + 1
      write (li(nl), *) 'Control_type: ', control_name(ct)

      if (ele%n_slave /= 0) then
        write (li(nl+1), '(1x, a, i4)') 'Slaves: Number:', ele%n_slave
        write (li(nl+2), *) &
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
          write (li(nl), '(5x, a, i10, 2x, a16, 1pe11.3, 1pe12.3)') &
                                ring%ele_(j)%name, j, a_name, coef
        enddo
      endif

      if (ele%n_lord /= 0) then
        write (li(nl+1), '(1x, a, i4)') 'Lords: Number:', ele%n_lord
        write (li(nl+2), *) &
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
          write (li(nl), '(5x, a, i10, 2x, a16, 1pe11.3, a12)') &
                             ring%ele_(j)%name, j, a_name, coef, val_str
        enddo
      endif

    endif
  endif

! Encode Twiss info

  if (twiss_out > 0) then
    call type2_twiss (ele, twiss_out, li(nl+1:), nl2)
    nl = nl + nl2
  endif

! Encode mat6 info

  n = type_mat6

  if (n /= 0) then
    nl = nl + 1
    write (li(nl), *)
  endif

  if (any(abs(ele%mat6(1:n,1:n)) >= 1000)) then
    do i = 1, n
      nl = nl + 1
      write (li(nl), '(1p6e11.3)') (ele%mat6(i, j), j = 1, n)
    enddo
  else
    do i = 1, n
      nl = nl + 1
      write (li(nl), '(6f10.5)') (ele%mat6(i, j), j = 1, n)
    enddo
  endif

! Encode taylor series

  if (type_taylor .and. ele%pointer_init == has_been_inited$ .and. &
                                    associated(ele%taylor(1)%term)) then
    call type2_taylors (ele%taylor, li2, nt)
    n_max = nl + nt + 1
    n_lines = n_max
    allocate (lines(n_max+1))
    write (li(nl+1), *)
    do i = 1, nt
      lines(i+nl+1) = li2(i)
    enddo
    deallocate (li2)
  else
    allocate(lines(nl))
    n_lines = nl
  endif

! finish

  do i = 1, nl
    lines(i) = li(i)
  enddo
  
  deallocate (li)
  
end subroutine
