!+
! Subroutine tao_show_cmd (s, word1, word2, word3, word4)
!
! Show information on variable, parameters, elements, etc...
!
!
! Input:
!   s        -- tao_super_universe_struct
!   word1 -- Character(*):
!   word2 -- Character(*):
!   word3 -- Character(*):
!   word4 -- Character(*):
!-

subroutine tao_show_cmd (s, word1, word2, word3, word4)

use tao_mod
use tao_common

implicit none

type (tao_super_universe_struct), target :: s
type (tao_universe_struct), pointer :: u
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (coord_struct) orb

type (ele_struct), pointer :: ele
type (ele_struct) ele3

real(rp) f_phi

integer, parameter :: max_lines = 300   ! maximum lines of output
character(16) :: max_lines_char = '300'   

character(24) :: data_class, var_class
character(24)  :: plane, fmt, imt, lmt, amt
character(*) :: word1, word2, word3, word4
character(40) :: show_word1, show_word2, show_word3, show_word4
character(8) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name
character(80), pointer :: ptr_lines(:)
character(16) ele_name

character(24) :: show_names(8) = (/ &
        'data     ', 'var      ', 'global   ',  'alias    ', 'top10    ', &
        'optimizer', 'ele      ', 'lattice  ' /)

character(200) lines(max_lines)

integer :: data_number, ix_plane
integer nl, loc
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index

logical err, found, at_ends
logical show_all, name_found

! The Show string data is stored in the 'lines' array
! and the number of entries is stored in nl

err = .false.

lines = " "
nl = 0

fmt = '(a, 1pe16.8)'
imt = '(a, i8)'
lmt = '(a, l)'
amt = '(2a)'

show_word1 = word1
show_word2 = word2
show_word3 = word3
show_word4 = word4

u => s%u(s%global%u_view)

if (s%global%phase_units == radians$) f_phi = 1
if (s%global%phase_units == degrees$) f_phi = 180 / pi
if (s%global%phase_units == cycles$)  f_phi = 1 / twopi

! find what to show

ix1 = index(show_word1, '/')
if (ix1 /= 0) show_word1 = show_word1(:ix1-1) // ' ' // show_word1(ix1:)

call string_trim (show_word1, show_word1, ix)
if (ix == 0) then
  call out_io (s_error$, r_name, 'SHOW WHAT?')
  return
else
  call match_word (show_word1, show_names, ix)
  if (ix == 0) then
    call out_io (s_error$, r_name, 'SHOW WHAT? WORD NOT RECOGNIZED: ' // show_word1)
    return
  endif
  show_name = show_names(ix)
endif


select case (show_name)

!----------------------------------------------------------------------
! alias

case ('alias')

  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, tao_com%n_alias
    nl=nl+1; lines(nl) = trim(tao_com%alias(i)%name) // ' = "' // &
                                    trim(tao_com%alias(i)%string) // '"'
    call out_io (s_blank$, r_name, lines(1:nl))
  enddo

!----------------------------------------------------------------------
! data

case ('data')

  if (.not. associated (u%d2_data)) return 

! If just "show data" then show all classes

  if (show_word2 == ' ') then
    do i = 1, size(u%d2_data)
      if (u%d2_data(i)%class == ' ') cycle
      if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many d2_datas! Listing first "// &
                                      		max_lines_char // " matches.")
	      exit
      endif
      nl=nl+1; write (lines(nl), '(i4, 2x, a)') i, u%d2_data(i)%class
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif


! get pointers to the data

  if (show_word3 == ' ') show_word3 = 'null'
  call tao_find_data (err, u, show_word2, d2_ptr, d1_ptr, show_word3, d_ptr)
  if (err) return

! If d_ptr points to something then show the datum info.

  if (associated(d_ptr)) then

    nl=nl+1; write(lines(nl), amt)  'Name:              ', d_ptr%name
    nl=nl+1; write(lines(nl), amt)  'Alias:             ', d_ptr%alias
    nl=nl+1; write(lines(nl), amt)  'Ele_name:          ', d_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  'Ele_name2:         ', d_ptr%ele_name2
    nl=nl+1; write(lines(nl), imt)  'Ix_ele:            ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)  'Ix_ele2:           ', d_ptr%ix_ele2
    nl=nl+1; write(lines(nl), imt)  'Ix_dModel:         ', d_ptr%ix_dModel
    nl=nl+1; write(lines(nl), imt)  'Ix_d:              ', d_ptr%ix_d
    nl=nl+1; write(lines(nl), imt)  'Ix_data:           ', d_ptr%ix_data
    nl=nl+1; write(lines(nl), fmt)  'Data_value:        ', d_ptr%data_value
    nl=nl+1; write(lines(nl), fmt)  'Ref_value:         ', d_ptr%ref_value
    nl=nl+1; write(lines(nl), fmt)  'Model_value:       ', d_ptr%model_value
    nl=nl+1; write(lines(nl), fmt)  'base_value:        ', d_ptr%base_value
    nl=nl+1; write(lines(nl), fmt)  'Delta:             ', d_ptr%delta
    nl=nl+1; write(lines(nl), fmt)  'Design_value:      ', d_ptr%design_value
    nl=nl+1; write(lines(nl), fmt)  'Old_value:         ', d_ptr%old_value
    nl=nl+1; write(lines(nl), fmt)  'Fit_value:         ', d_ptr%fit_value
    nl=nl+1; write(lines(nl), fmt)  'Merit:             ', d_ptr%merit
    nl=nl+1; write(lines(nl), fmt)  'Conversion_factor: ', d_ptr%conversion_factor
    nl=nl+1; write(lines(nl), fmt)  'Weight:            ', d_ptr%weight
    nl=nl+1; write(lines(nl), amt)  'Merit_type:        ', d_ptr%merit_type
    nl=nl+1; write(lines(nl), lmt)  'Exists:            ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  'Good_data:         ', d_ptr%good_data
    nl=nl+1; write(lines(nl), lmt)  'Good_ref:          ', d_ptr%good_ref
    nl=nl+1; write(lines(nl), lmt)  'Good_user:         ', d_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  'Useit_plot:        ', d_ptr%useit_plot
    nl=nl+1; write(lines(nl), lmt)  'Useit_opt:         ', d_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '[d2_data%good_opt:]', d_ptr%d1%d2%good_opt

! Else show the d1_data info.

  elseif (associated(d1_ptr)) then

    write(lines(1), '(2a)') 'Data type:    ', d2_ptr%class
    write(lines(2), '(2a)') 'Data sub_class:', d1_ptr%sub_class
    lines(3) = ' '
    write (lines(3), *) '        Name              Data         Model        Design'
    nl = 4
    do i = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
      if ((nl+1) .gt. max_lines) then
      call out_io (s_abort$, r_name, "Found too many datams! Listing first "//&
                    max_lines_char// " matches.")
      exit
      endif
      if (d1_ptr%d(i)%exists) then
        nl=nl+1
        write(lines(nl), '(i, 2x, a7, 3es14.4)') i, &
               d1_ptr%d(i)%name, d1_ptr%d(i)%data_value, &
               d1_ptr%d(i)%model_value, d1_ptr%d(i)%design_value
      endif
    enddo
    nl=nl+1
    write (lines(nl), *) '        Name              Data         Model        Design'

! else we must have a valid d2_ptr.

  else 

    write(lines(1), '(2a)') 'Data type:    ', d2_ptr%class
    do i = 1, size(d2_ptr%d1)
      if ((nl+1) .gt. max_lines) then
      call out_io (s_abort$, r_name, "Found too many d1_data! Listing first "//&
                    max_lines_char// " matches.")
      exit
      endif
      nl=nl+1; write(lines(nl), '(5x, a, i3, a, i3)') d2_ptr%d1(i)%sub_class, &
                lbound(d2_ptr%d1(i)%d, 1), ':', ubound(d2_ptr%d1(i)%d, 1)
    enddo
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! ele

case ('ele')


  call str_upcase (ele_name, show_word2)

  if (index(ele_name, '*') /= 0 .or. index(ele_name, '%') /= 0) then
    write (lines(1), *) 'Matches to name:'
    nl = 1
    do loc = 1, u%model%n_ele_max
      if (match_wild(u%model%ele_(loc)%name, ele_name)) then
      if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many elements! Listing first "//&
                    max_lines_char// " matches.")
        exit
      endif
        nl = nl + 1
        write (lines(nl), '(i8, 2x, a)') loc, u%model%ele_(loc)%name
        name_found = .true.
      endif
    enddo
    if (.not. name_found) then
      nl = nl + 1
      write (lines(nl), *) '   *** No Matches to Name Found ***'
    endif

! else no wild cards

  else  

    call tao_locate_element (ele_name, u%model, loc)
    if (loc < 0) return

    write (lines(nl+1), *) 'Element #', loc
    nl = nl + 1

    call type2_ele (u%model%ele_(loc), .true., 6, .false., &
                  s%global%phase_units, .true., ptr_lines, n, u%model)
    lines(nl+1:nl+n) = ptr_lines(1:n)
    nl = nl + n
    deallocate (ptr_lines)

    orb = u%model_orb(loc)
    fmt = '(2x, a, 3p2f11.4)'
    write (lines(nl+1), *) ' '
    write (lines(nl+2), *)   'Orbit: [mm, mrad]'
    write (lines(nl+3), fmt) "X  X':", orb%vec(1:2)
    write (lines(nl+4), fmt) "Y  Y':", orb%vec(3:4)
    write (lines(nl+5), fmt) "Z  Z':", orb%vec(5:6)
    nl = nl + 5

    found = .false.
    do i = loc + 1, u%model%n_ele_max
      if (u%model%ele_(i)%name == ele_name) then
        if (found) then
          nl = nl + 1
          write (lines(nl), *)
          found = .true.
        endif
        nl = nl + 1
        write (lines(nl), *) &
                        'Note: Found another element with same name at:', i
      endif
    enddo

  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! global

case ('global')

  nl=nl+1; write (lines(nl), imt) 'n_universes: ', size(s%u)
  nl=nl+1; write (lines(nl), imt) 'u_view:      ', s%global%u_view
  nl=nl+1; write (lines(nl), amt) 'optimizer:     ', s%global%optimizer
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on: ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), lmt) 'opt_with_ref:  ', s%global%opt_with_ref 
  nl=nl+1; write (lines(nl), lmt) 'opt_with_base: ', s%global%opt_with_base

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! lattice

case ('lattice')

  call tao_locate_element (show_word3, u%model, ix1); if (ix1 < 0) return
  call tao_locate_element (show_word4, u%model, ix2); if (ix2 < 0) return
  
  if (ix1 > u%model%n_ele_ring) then
    ix = u%model%ele_(ix1)%ix1_slave
    ix1 = u%model%control_(ix)%ix_slave
  endif

  if (ix2 > u%model%n_ele_ring) then
    ix = u%model%ele_(ix2)%ix2_slave
    ix2 = u%model%control_(ix)%ix_slave
  endif

  if ((ix2 - ix1 + 1) .gt. max_lines) then
    call out_io (s_fatal$, r_name, "Too many elements!")
    call out_io (s_fatal$, r_name, "Listing first " //&
        max_lines_char // " elements after element" // show_word3 )
    ix2 = ix1 + max_lines - 1
  endif

  if (.true.) then
    at_ends = .true.
    write (lines(nl+1), '(36x, a)') 'Model values at End of Element:'
  else
    at_ends = .false.
    write (lines(nl+1), '(36x, a)') 'Model values at Center of Element:'
  endif


  write (lines(nl+2), '(28x, a)') &
                     '|              X           |             Y        '
  write (lines(nl+3), '(5x, a, 16x, a)') 'Name', &
                  'S  |  Beta     Phi  Eta   Orb | Beta     Phi  Eta   Orb'

  nl=nl+3
  do ix = ix1, ix2
    ele => u%model%ele_(ix)
    if (ix == 0 .or. at_ends) then
      ele3 = ele
      orb = u%model_orb(ix)
    else
      call twiss_and_track_partial (u%model%ele_(ix-1), ele, &
                u%model%param, ele%value(l$)/2, ele3, u%model_orb(ix-1), orb)
    endif
    nl=nl+1
    write (lines(nl), '(i4, 1x, a16, f8.3, 2(f7.2, f8.3, f5.1, f8.3))') &
          ix, ele%name, ele%s-ele%value(l$)/2, &
          ele3%x%beta, f_phi*ele3%x%phi, ele3%x%eta, 1000*orb%vec(1), &
          ele3%y%beta, f_phi*ele3%y%phi, ele3%y%eta, 1000*orb%vec(3)
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! optimizer

case ('optimizer')

  do i = 1, size(s%u)
    u => s%u(i)
    call out_io (s_blank$, r_name, ' ', 'Data Used:')
    write (lines(1), '(a, i4)') 'Universe: ', i
    if (size(s%u) > 1) call out_io (s_blank$, r_name, lines(1))
    do j = 1, size(u%d2_data)
      if (u%d2_data(j)%class == ' ') cycle
      call tao_data_show_use (u%d2_data(j))
    enddo
  enddo

  call out_io (s_blank$, r_name, ' ', 'Variables Used:')
  if (s%global%parallel_vars) then
    u => s%u(1)
    do j = 1, size(u%v1_var)
      if (u%v1_var(j)%class == ' ') cycle
      call tao_var_show_use (u%v1_var(j))
    enddo
  else
    do i = 1, size(s%u)
      u => s%u(i)
      write (lines(1), '(a, i4)') 'Universe: ', i
      call out_io (s_blank$, r_name, lines(1))
      do j = 1, size(u%v1_var)
        if (u%v1_var(j)%class == ' ') cycle
        call tao_var_show_use (u%v1_var(j))
      enddo
    enddo
  endif


!----------------------------------------------------------------------
! top10

case ('top10')

  call tao_top10_print (s)

!----------------------------------------------------------------------
! variable
    
case ('var')

  if (.not. associated (u%v1_var)) return 

! If just "show var" then show all classes

  if (show_word2 == ' ') then
    do i = 1, size(u%v1_var)
      v1_ptr => u%v1_var(i)
      if (v1_ptr%class == ' ') cycle
      if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many v1_vars! Listing first "//&
                    max_lines_char// " matches")
        exit
      endif
      nl = nl+1
      write (lines(nl), '(i4, 2x, 2a)') i, v1_ptr%class, v1_ptr%class
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the variables

  if (show_word3 == ' ') show_word3 = 'null'
  call tao_find_var(err, u, show_word2, v1_ptr, show_word3, v_ptr) 
  if (err) return

! v_ptr is valid then show the variable info.

  if (associated(v_ptr)) then

    nl=nl+1; write(lines(nl), amt)  'Name:          ', v_ptr%name        
    nl=nl+1; write(lines(nl), amt)  'Alias:         ', v_ptr%alias       
    nl=nl+1; write(lines(nl), amt)  'Ele_name:      ', v_ptr%ele_name    
    nl=nl+1; write(lines(nl), amt)  'Attrib_name:   ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  'Ix_ele:        ', v_ptr%ix_ele            
    nl=nl+1; write(lines(nl), imt)  'Ix_v:          ', v_ptr%ix_v
    nl=nl+1; write(lines(nl), imt)  'Ix_var:        ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  'Ix_dvar:       ', v_ptr%ix_dvar           
    if (associated (v_ptr%model_value)) then
      nl=nl+1; write(lines(nl), fmt)  'Model_value:   ', v_ptr%model_value
    else
      nl=nl+1; write(lines(nl), fmt)  'Model_value:   <not associated>'
    endif
    if (associated (v_ptr%base_value)) then
      nl=nl+1; write(lines(nl), fmt)  'base_value:    ', v_ptr%base_value
    else
      nl=nl+1; write(lines(nl), fmt)  'base_value:    <not associated>'
    endif
    nl=nl+1; write(lines(nl), fmt)  'Design_value:  ', v_ptr%design_value     
    nl=nl+1; write(lines(nl), fmt)  'Old_value:     ', v_ptr%old_value        
    nl=nl+1; write(lines(nl), fmt)  'Data_value:    ', v_ptr%data_value       
    nl=nl+1; write(lines(nl), fmt)  'Ref_value:     ', v_ptr%ref_value        
    nl=nl+1; write(lines(nl), fmt)  'Target_value:  ', v_ptr%target_value     
    nl=nl+1; write(lines(nl), fmt)  'High_lim_value:', v_ptr%high_lim_value   
    nl=nl+1; write(lines(nl), fmt)  'Low_lim_value: ', v_ptr%low_lim_value    
    nl=nl+1; write(lines(nl), fmt)  'Step:          ', v_ptr%step             
    nl=nl+1; write(lines(nl), fmt)  'Weight:        ', v_ptr%weight           
    nl=nl+1; write(lines(nl), fmt)  'Merit:         ', v_ptr%merit            
    nl=nl+1; write(lines(nl), fmt)  'DMerit_dVar:   ', v_ptr%dMerit_dVar      
    nl=nl+1; write(lines(nl), lmt)  'Exists:        ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  'Good_var:      ', v_ptr%good_var  
    nl=nl+1; write(lines(nl), lmt)  'Good_user:     ', v_ptr%good_user 
    nl=nl+1; write(lines(nl), lmt)  'Useit_opt:     ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  'Useit_plot:    ', v_ptr%useit_plot   
    nl=nl+1; write(lines(nl), lmt)  '[v1_var%Good_opt:]', v_ptr%v1%good_opt 

! check if there is a variable number
! if no variable number requested, show it all

  else

    write(lines(1), '(2a)') 'Variable class:   ', v1_ptr%class
    lines(2) = ' '
    write (lines(3), *) '        Name              Data         Model        Design'
    nl = 3
    do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
      if (v1_ptr%v(i)%exists) then
        if ((nl+1) .gt. max_lines) then
          call out_io (s_abort$, r_name, "Found too many variable datams! Listing first "// &
                    max_lines_char// " matches.")
          exit
        endif
        nl = nl + 1
        write(lines(nl), '(i6, 2x, a7, 3es14.4)') i, &
                 v1_ptr%v(i)%name, v1_ptr%v(i)%data_value, &
                 v1_ptr%v(i)%model_value, v1_ptr%v(i)%design_value
      endif
    enddo
    nl=nl+1
    write (lines(nl), *) '        Name              Data         Model        Design'
  endif

! print out results

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "INTERNAL ERROR, SHOULDN'T BE HERE!")
  return

end select

end subroutine tao_show_cmd
