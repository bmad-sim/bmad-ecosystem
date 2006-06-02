!+
! Subroutine type2_ele (ele, lines, n_lines, type_zero_attrib, type_mat6, 
!                      type_taylor, twiss_out, type_control, lattice, type_wake)
!
! Subroutine to put information on an element in a string array. 
! See also the subroutine: type_ele.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele              -- Ele_struct: Element
!   type_zero_attrib -- Logical, optional: If False then surpress printing of
!                          attributes whose value is 0. Default is False.
!   type_mat6      -- Integer, optional:
!                         = 0   => Do not type ele%mat6
!                         = 4   => Type 4X4 xy submatrix
!                         = 6   => Type full 6x6 matrix (Default)
!   type_taylor    -- Logical, optional: Print out taylor map terms?
!                       If ele%taylor is not allocated then this is ignored.
!                       Default is False.
!   twiss_out      -- Integer, optional: Print the Twiss parameters at the 
!                         element end?
!                       = 0         => Do not print the Twiss parameters
!                       = radians$  => Print Twiss, phi in radians (Default).
!                       = degrees$  => Print Twiss, phi in degrees.
!                       = cycles$   => Print Twiss, phi in radians/2pi.
!   type_control   -- Logical, optional: If True then print control status.
!                       Default: False if lattice is not present. Otherwise True.
!   lattice        -- Ring_struct, optional: Needed for control typeout.
!   type_wake      -- Logical, optional: If True then print the long-range and 
!                       short-range wakes information. If False then just print
!                       how many terms the wake has. Default is True.
!                       If ele%wake is not allocated then this is ignored.
!
! Output       
!   lines(:)     -- Character(80), pointer: Character array to hold the 
!                     output. The array size of lines(:) will be set by
!                     this subroutine. Note: You need to deallocate lines
!                     between each call to type2_ele.
!   n_lines      -- Number of lines in lines(:).
!-

#include "CESR_platform.inc"

subroutine type2_ele (ele, lines, n_lines, type_zero_attrib, type_mat6, &
                        type_taylor, twiss_out, type_control, lattice, type_wake)

  use bmad_struct
  use bmad_interface, except => type2_ele
  use multipole_mod

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (ring_struct), optional, intent(in) :: lattice
  type (wig_term_struct), pointer :: term
  type (lr_wake_struct), pointer :: lr
  type (sr1_wake_struct), pointer :: sr1
  type (sr2_wake_struct), pointer :: sr2

  integer, optional, intent(in) :: type_mat6, twiss_out
  integer, intent(out) :: n_lines
  integer i, j, n, ix, iv, ic, ct, nl2
  integer nl, nt, n_max, particle
  integer pos_tot(n_attrib_maxx)

  real(rp) coef
  real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rp) a2(0:n_pole_maxx), b2(0:n_pole_maxx)

  character(100), pointer :: lines(:)
  character(100), pointer :: li(:), li2(:)
  character(40) a_name, name
  character(12) val_str
  character(9) angle
  character(2) str_i

  logical, optional, intent(in) :: type_taylor, type_wake
  logical, optional, intent(in) :: type_control, type_zero_attrib
  logical type_zero

! init

  allocate (li(300))

  pos_tot = 0
  if (ele%control_type /= group_lord$  .and. ele%control_type /= overlay_lord$) then
    pos_tot(x_offset$) = x_offset_tot$
    pos_tot(y_offset$) = y_offset_tot$
    pos_tot(s_offset$) = s_offset_tot$
    pos_tot(tilt$)     = tilt_tot$
    pos_tot(x_pitch$)  = x_pitch_tot$
    pos_tot(y_pitch$)  = y_pitch_tot$
  endif

  type_zero = logic_option(.false., type_zero_attrib)

! Encode element name and type

  ct = ele%control_type
  
  write (li(1), *) 'Element Name: ', ele%name
  nl = 1

  if (ele%type /= blank_name) then
    nl=nl+1; write (li(nl), *) 'Element Type: "', ele%type, '"'
  endif

  if (ele%alias /= blank_name) then
    nl=nl+1; write (li(nl), *) 'Element Alias: "', ele%alias, '"'
  endif

  if (associated(ele%descrip)) then
    nl=nl+1; write (li(nl), *) 'Descrip: ', trim(ele%descrip)
  endif

! Encode element key and attributes

  if (ele%key <= 0) then
    nl=nl+1; write (li(nl), *) 'Key: UNKNOWN!', ele%key

  else

    nl=nl+1; write (li(nl), *) 'Key: ', key_name(ele%key)

    if (ele%sub_key /= 0) then
      nl=nl+1; write (li(nl), *) 'Sub Key: ', sub_key_name(ele%sub_key)
    endif

    if (.not. type_zero) then
      nl=nl+1; write (li(nl), *) 'Attribute values [Only non-zero values shown]:'
    endif

    if (ct == overlay_lord$) then
      i = ele%ix_value
      name = ele%attribute_name
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, 1pe15.7)') i, name, ' =', ele%value(i)

    else
      do i = 1, n_attrib_maxx
        if (attribute_name(ele, i) == null_name) cycle
        ix = pos_tot(i)
        if (ix == 0) then
          if (ele%value(i) == 0 .and. .not. type_zero) cycle
          nl=nl+1; write (li(nl), '(i6, 3x, 2a, 1pe15.7)')  i, &
                        attribute_name(ele, i), ' =', ele%value(i)
        else
          if (ele%value(i) == 0 .and. ele%value(ix) == 0 .and. &
                                                 .not. type_zero) cycle
          nl=nl+1; write (li(nl), '(i6, 3x, 2a, 1pe15.7, 4x, a, e15.7)')  i, &
                        attribute_name(ele, i), ' =', ele%value(i), &
                        'Total:', ele%value(ix)
        endif
      enddo

      nl=nl+1; write (li(nl), '(9x, 2a)')   'APERTURE_AT      = ', &
                                  element_end_name(ele%aperture_at)
      if (ele%key == lcavity$) then
        nl=nl+1; write (li(nl), '(9x, 2a)') 'COUPLER_AT       = ', &
                                  element_end_name(ele%coupler_at)
      endif

      if (associated(ele%a)) then
        particle = +1
        if (present(lattice)) particle = lattice%param%particle

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
    nl=nl+1; write (li(nl), '(a, 6x, a, 3(9x, a), 7x, a)') ' Term#', &
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
    nl=nl+1; write (li(nl), *) '*** Note: Element is turned OFF ***'
  endif

  if (.not. ele%multipoles_on) then
    nl=nl+1; write (li(nl), *) '*** Note: Element Multipoles are turned OFF ***'
  endif

  nl=nl+1; write (li(nl), '(1x, a, f13.4)') 'S:', ele%s

! Encode methods, etc.

  nl=nl+1; write (li(nl), *) ' '

  if (attribute_index(ele, 'TRACKING_METHOD') /= 0) then
    nl=nl+1; write (li(nl), '(2a)') ' Tracking_method:  ', &
                                    calc_method_name(ele%tracking_method)
  endif

  if (attribute_index(ele, 'MAT6_CALC_METHOD') /= 0) then
    nl=nl+1; write (li(nl), '(2a)') ' Mat6_calc_method: ', &
                                    calc_method_name(ele%mat6_calc_method)
  endif

  if (attribute_index(ele, 'FIELD_CALC') /= 0) then
    nl=nl+1; write (li(nl), '(2a)') ' Field_calc:       ', &
                                    calc_method_name(ele%field_calc)
  endif

  if (attribute_index(ele, 'INTEGRATOR_ORDER') /= 0) then
    nl=nl+1; write (li(nl), '(a, i4)') ' Integrator_order:', ele%integrator_order 
  endif

  if (attribute_index(ele, 'NUM_STEPS') /= 0) then
    write (li(nl), '(a, i4)') ' Num_steps:       ', ele%num_steps 
  endif

  if (attribute_index(ele, 'SYMPLECTIFY') /= 0) then
    nl=nl+1; write (li(nl), '(a, l1)') ' Symplectify:       ', ele%symplectify
  endif
  
! Encode lord info

  if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
    write (li(nl+1), *) ' '
    write (li(nl+2), *) &
                   'ERROR: N_SLAVE DOES NOT MATCH IX2_SLAVE-IX1_SLAVE+1 !!'
    write (li(nl+3), *) ele%n_slave, ele%ix2_slave, ele%ix1_slave
    nl = nl + 3
  endif          

  if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
    write (li(nl+1), *) ' '
    write (li(nl+2), *) 'ERROR: N_LORD DOES NOT MATCH Ic1_LORD-Ic2_LORD+1 !!'
    write (li(nl+3), *) ele%n_lord, ele%ic1_lord, ele%ic2_lord
    nl = nl + 3
  endif

! Encode slave info.
! For super_lords there is no attribute_name associated with a slave.
! For slaves who are overlay_lords then the attribute_name is obtained by
!   looking at the overlay_lord's 1st slave (slave of slave of the input ele).

  if (logic_option(present(lattice), type_control)) then

    if (.not. present (lattice)) then
      print *, 'ERROR IN TYPE2_ELE: TYPE_CONTROL IS TRUE BUT NO LATTICE PRESENT.'
      call err_exit
    endif

    nl=nl+1; write (li(nl), *) ' '

    if (ct <= 0) then
      nl=nl+1; write (li(nl), *) 'Control_type: UNKNOWN!', ct
    else
      nl=nl+1; write (li(nl), *) 'Control_type: ', control_name(ct)

      if (ele%n_slave /= 0) then
        write (li(nl+1), '(1x, a, i4)') 'Slaves: Number:', ele%n_slave
        write (li(nl+2), *) &
          '    Name             Lat_index  Attribute       Coefficient'
        nl = nl + 2
        do i = ele%ix1_slave, ele%ix2_slave
          j = lattice%control_(i)%ix_slave
          iv = lattice%control_(i)%ix_attrib
          coef = lattice%control_(i)%coef
          select case (ct)
          case (super_lord$, i_beam_lord$, multipass_lord$) 
            a_name = '--------'
          case default
            if (lattice%ele_(j)%control_type == overlay_lord$) then
              if (iv == lattice%ele_(j)%ix_value) then
                ix = lattice%control_(lattice%ele_(j)%ix1_slave)%ix_slave
                a_name = attribute_name(lattice%ele_(ix), iv)
              else
                a_name = '** BAD POINTER! **'
              endif            
            else
              a_name = attribute_name(lattice%ele_(j), iv)
            endif
          end select
          nl=nl+1; write (li(nl), '(5x, a, i10, 2x, a20, es11.3, es12.3)') &
                                lattice%ele_(j)%name, j, a_name, coef
        enddo
      endif

      if (ele%n_lord /= 0) then
        write (li(nl+1), '(1x, a, i4)') 'Lords: Number:', ele%n_lord
        write (li(nl+2), *) &
  '    Name             Lat_index  Attribute       Coefficient       Value'
        nl = nl + 2
        do i = ele%ic1_lord, ele%ic2_lord
          ic = lattice%ic_(i)
          j = lattice%control_(ic)%ix_lord
          coef = lattice%control_(ic)%coef
          ct = ele%control_type
          if (ct == super_slave$ .or. ct == multipass_slave$ .or. &
                            lattice%ele_(j)%control_type == i_beam_lord$) then
            a_name = '--------'
            val_str = '    --------'
          else
            iv = lattice%control_(ic)%ix_attrib
            a_name = attribute_name(ele, iv)
            ix = lattice%ele_(j)%ix_value
            if (ix == 0) then
              val_str = '  GARBAGE!'
            else
              write (val_str, '(1p, e12.3)') lattice%ele_(j)%value(ix)
            endif
          endif
          nl=nl+1; write (li(nl), '(5x, a, i10, 2x, a20, es11.3, a12)') &
                             lattice%ele_(j)%name, j, a_name, coef, val_str
        enddo
      endif

    endif
  endif

! Encode Twiss info

  if (integer_option(radians$, twiss_out) /= 0) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; li(nl) = 'Twiss at end of element:'
    call type2_twiss (ele, li(nl+1:), nl2, twiss_out)
    nl = nl + nl2
  endif

! Encode mat6 info

  n = integer_option (6, type_mat6)

  if (n /= 0) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; li(nl) = 'Transfer Matrix:'
  endif

  if (any(abs(ele%mat6(1:n,1:n)) >= 1000)) then
    do i = 1, n
      nl=nl+1; write (li(nl), '(1p, 6e11.3)') (ele%mat6(i, j), j = 1, n)
    enddo
  else
    do i = 1, n
      nl=nl+1; write (li(nl), '(6f10.5)') (ele%mat6(i, j), j = 1, n)
    enddo
  endif

! Encode taylor series

  if (associated(ele%taylor(1)%term)) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; write (li(nl), '(a, l1)') &
                        'map_with_offsets: ', ele%map_with_offsets
    if (logic_option(.false., type_taylor)) then
      call type2_taylors (ele%taylor, li2, nt)
      call re_associate (li,  len(li(1)), nl+nt+100)
      li(1+nl:nt+nl) = li2(1:nt)
      deallocate (li2)
      nl = nl + nt
    else
      ct = 0
      do i = 1, 6
        ct = ct + size(ele%taylor(i)%term)
      enddo
      nl=nl+1; write (li(nl), *) 'Taylor map total number of terms:', ct
    endif
  endif

! Encode HOM info

  if (associated(ele%wake)) then

    if (size(ele%wake%sr1) /= 0) then
      nl=nl+1; write (li(nl), *)
      if (logic_option (.true., type_wake)) then
        call re_associate (li,  len(li(1)), nl+size(ele%wake%sr1)+100)
        nl=nl+1; li(nl) = 'Short-range wake table:'
        nl=nl+1; li(nl) = &
            '   #           Z   Longitudinal     Transverse'
        do i = 0, ubound(ele%wake%sr1,1)
          sr1 => ele%wake%sr1(i)
          nl=nl+1; write (li(nl), '(i4, es12.4, 2es15.4)') i, sr1%z, sr1%long, sr1%trans
        enddo
      else
        nl=nl+1; write (li(nl), *) 'Number of short-range wake table rows:', size(ele%wake%sr1)
      endif
    endif

    if (size(ele%wake%sr2_long) /= 0) then
      nl=nl+1; write (li(nl), *)
      if (logic_option (.true., type_wake)) then
        nl=nl+1; li(nl) = 'Short-range pseudo modes:'
        nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi'
        do i = 1, size(ele%wake%sr2_long)
          sr2 => ele%wake%sr2_long(i)
          nl=nl+1; write (li(nl), '(i4, 4es12.4)') i, sr2%amp, sr2%damp, sr2%k, sr2%phi
        enddo
      else
        nl=nl+1; write (li(nl), *) &
                  'Number of short-range longitudinal pseudo modes:', size(ele%wake%sr2_long)
      endif
    endif

    if (size(ele%wake%sr2_trans) /= 0) then
      nl=nl+1; write (li(nl), *)
      if (logic_option (.true., type_wake)) then
        nl=nl+1; li(nl) = 'Short-range pseudo modes:'
        nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi'
        do i = 1, size(ele%wake%sr2_trans)
          sr2 => ele%wake%sr2_trans(i)
          nl=nl+1; write (li(nl), '(i4, 4es12.4)') i, sr2%amp, sr2%damp, sr2%k, sr2%phi
        enddo
      else
        nl=nl+1; write (li(nl), *) &
                  'Number of short-range transitudinal pseudo-modes:', size(ele%wake%sr2_trans)
      endif
    endif

    if (size(ele%wake%lr) /= 0) then
      nl=nl+1; write (li(nl), *)
      if (logic_option (.true., type_wake)) then
        nl=nl+1; li(nl) = 'Long-range HOM modes:'
        nl=nl+1; li(nl) = &
            '   #        Freq         R/Q           Q   m  Polarization_Angle'
        do i = 1, size(ele%wake%lr)
          lr => ele%wake%lr(i)
          angle = '-'
          if (lr%polarized) write (angle, '(f9.4)') lr%angle
          nl=nl+1; write (li(nl), '(i4, 3es12.4, i4, a)') i, &
                  lr%freq, lr%R_over_Q, lr%Q, lr%m, angle
        enddo
      else
        nl=nl+1; write (li(nl), *) 'Number of long-range HOM modes:', size(ele%wake%lr)
      endif
    endif

  endif

! finish

  allocate(lines(nl))
  n_lines = nl
  lines(1:nl) = li(1:nl)
  deallocate (li)
  
end subroutine
