!+
! Subroutine tao_show_cmd (word1, word2, word3, word4)
!
! Show information on variable, parameters, elements, etc...
!
!
! Input:
!   word1 -- Character(*):
!   word2 -- Character(*):
!   word3 -- Character(*):
!   word4 -- Character(*):
!-

subroutine tao_show_cmd (word1, word2, word3, word4)

use tao_mod
use tao_top10_mod

implicit none

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

character(24) :: var_name
character(24)  :: plane, fmt, imt, lmt, amt, ffmt, iimt
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
character(100) line1, line2

integer :: data_number, ix_plane
integer nl, loc
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju

logical err, found, at_ends
logical show_all, name_found
logical, automatic :: picked(size(s%u))

! The Show string data is stored in the 'lines' array
! and the number of entries is stored in nl

err = .false.

lines = " "
nl = 0

fmt  = '(a, es16.8)'
ffmt = '(a, i1, a, es16.8)'
imt  = '(a, i8)'
iimt = '(a, i1, a, i8)'
lmt  = '(a, l)'
amt  = '(2a)'

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

  call tao_pick_universe (show_word2, show_word2, picked, err)
  if (err) return

  u_loop: do ju = 1, size(s%u)

    if (.not. picked(ju)) cycle
    if (.not. associated (s%u(ju)%d2_data)) return 
    u => s%u(ju)

! If just "show data" then show all namees

    if (show_word2 == ' ') then
      do i = 1, size(u%d2_data)
        if (u%d2_data(i)%name == ' ') cycle
        if ((nl+1) .gt. max_lines) then
          call out_io (s_abort$, r_name, "Found too many d2_datas! Listing first "// &
                                        		max_lines_char // " matches.")
	        exit
        endif
        nl=nl+1; write (lines(nl), '(i4, 2x, a)') i, u%d2_data(i)%name
      enddo
      call out_io (s_blank$, r_name, lines(1:nl))
      cycle u_loop
    endif


! get pointers to the data

    if (show_word3 == ' ') show_word3 = 'null'
    call tao_find_data (err, u, show_word2, d2_ptr, d1_ptr, show_word3, d_ptr)
    if (err) return

! If d_ptr points to something then show the datum info.

    if (associated(d_ptr)) then

      nl=nl+1; write(lines(nl), amt)  '%Name:              ', d_ptr%name
      nl=nl+1; write(lines(nl), amt)  '%Ele_name:          ', d_ptr%ele_name
      nl=nl+1; write(lines(nl), amt)  '%Ele2_name:         ', d_ptr%ele2_name
      nl=nl+1; write(lines(nl), amt)  '%Data_type:         ', d_ptr%data_type
      nl=nl+1; write(lines(nl), imt)  '%Ix_ele:            ', d_ptr%ix_ele
      nl=nl+1; write(lines(nl), imt)  '%Ix_ele2:           ', d_ptr%ix_ele2
      nl=nl+1; write(lines(nl), imt)  '%Ix_ele_merit:      ', d_ptr%ix_ele_merit
      nl=nl+1; write(lines(nl), imt)  '%Ix_dModel:         ', d_ptr%ix_dModel
      nl=nl+1; write(lines(nl), imt)  '%Ix_d1:             ', d_ptr%ix_d1
      nl=nl+1; write(lines(nl), imt)  '%Ix_data:           ', d_ptr%ix_data
      nl=nl+1; write(lines(nl), fmt)  '%Data_value:        ', d_ptr%data_value
      nl=nl+1; write(lines(nl), fmt)  '%Ref_value:         ', d_ptr%ref_value
      nl=nl+1; write(lines(nl), fmt)  '%Model_value:       ', d_ptr%model_value
      nl=nl+1; write(lines(nl), fmt)  '%base_value:        ', d_ptr%base_value
      nl=nl+1; write(lines(nl), fmt)  '%Delta:             ', d_ptr%delta
      nl=nl+1; write(lines(nl), fmt)  '%Design_value:      ', d_ptr%design_value
      nl=nl+1; write(lines(nl), fmt)  '%Old_value:         ', d_ptr%old_value
      nl=nl+1; write(lines(nl), fmt)  '%Fit_value:         ', d_ptr%fit_value
      nl=nl+1; write(lines(nl), fmt)  '%Merit:             ', d_ptr%merit
      nl=nl+1; write(lines(nl), fmt)  '%Conversion_factor: ', d_ptr%conversion_factor
      nl=nl+1; write(lines(nl), fmt)  '%S:                 ', d_ptr%s
      nl=nl+1; write(lines(nl), fmt)  '%Weight:            ', d_ptr%weight
      nl=nl+1; write(lines(nl), amt)  '%Merit_type:        ', d_ptr%merit_type
      nl=nl+1; write(lines(nl), lmt)  '%Exists:            ', d_ptr%exists
      nl=nl+1; write(lines(nl), lmt)  '%Good_data:         ', d_ptr%good_data
      nl=nl+1; write(lines(nl), lmt)  '%Good_ref:          ', d_ptr%good_ref
      nl=nl+1; write(lines(nl), lmt)  '%Good_user:         ', d_ptr%good_user
      nl=nl+1; write(lines(nl), lmt)  '%Good_opt:          ', d_ptr%good_opt
      nl=nl+1; write(lines(nl), lmt)  '%Good_plot:         ', d_ptr%good_plot
      nl=nl+1; write(lines(nl), lmt)  '%Useit_plot:        ', d_ptr%useit_plot
      nl=nl+1; write(lines(nl), lmt)  '%Useit_opt:         ', d_ptr%useit_opt

! Else show the d1_data info.

    elseif (associated(d1_ptr)) then

      write(lines(1), '(2a)') 'Data name: ', trim(d2_ptr%name) // ':' // d1_ptr%name
      lines(2) = ' '
      line1 = '        Name              Data         Model        Design'
      write (lines(3), *) line1
      nl = 3
      do i = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
        if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many datams! Listing first " // &
                                   max_lines_char // " matches.")
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
      write (lines(nl), *) line1

! else we must have a valid d2_ptr.

    else 

      write(lines(1), '(2a)') 'Data type:    ', d2_ptr%name
      do i = 1, size(d2_ptr%d1)
        if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many d1_data! Listing first " // &
                                       max_lines_char // " matches.")
        exit
        endif
        nl=nl+1; write(lines(nl), '(5x, a, i3, a, i3)') d2_ptr%d1(i)%name, &
                  lbound(d2_ptr%d1(i)%d, 1), ':', ubound(d2_ptr%d1(i)%d, 1)
      enddo
    endif

    call out_io (s_blank$, r_name, lines(1:nl))

  enddo u_loop

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
      if (u%d2_data(j)%name == ' ') cycle
      call tao_data_show_use (u%d2_data(j))
    enddo
  enddo

  call out_io (s_blank$, r_name, ' ', 'Variables Used:')
  do j = 1, size(s%v1_var)
    if (s%v1_var(j)%name == ' ') cycle
    call tao_var_show_use (s%v1_var(j))
  enddo


!----------------------------------------------------------------------
! top10

case ('top10')

  call tao_top10_print ()

!----------------------------------------------------------------------
! variable
    
case ('var')

  if (.not. associated (s%v1_var)) return 

! If just "show var" then show all namees

  if (show_word2 == ' ') then
    do i = 1, size(s%v1_var)
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      if ((nl+1) .gt. max_lines) then
        call out_io (s_abort$, r_name, "Found too many v1_vars! Listing first "//&
                    max_lines_char// " matches")
        exit
      endif
      nl = nl+1
      write (lines(nl), '(i4, 2x, 2a)') i, v1_ptr%name, v1_ptr%name
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the variables

  if (show_word3 == ' ') show_word3 = 'null'
  call tao_find_var(err, show_word2, v1_ptr, show_word3, v_ptr) 
  if (err) return

! v_ptr is valid then show the variable info.

  if (associated(v_ptr)) then

    nl=nl+1; write(lines(nl), amt)  'Name:          ', v_ptr%name        
    nl=nl+1; write(lines(nl), amt)  'Alias:         ', v_ptr%alias       
    nl=nl+1; write(lines(nl), amt)  'Ele_name:      ', v_ptr%ele_name    
    nl=nl+1; write(lines(nl), amt)  'Attrib_name:   ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  'Ix_var:        ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  'Ix_dvar:       ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  'Ix_v1:         ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), fmt)  'Model_value:   ', v_ptr%model_value
    nl=nl+1; write(lines(nl), fmt)  'Base_value:    ', v_ptr%base_value

    if (.not. associated (v_ptr%this)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      do i = 1, size(v_ptr%this)
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%Ix_uni:        ', &
                                                            v_ptr%this(i)%ix_uni
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%Ix_ele:        ', v_ptr%this(i)%ix_ele
        if (associated (v_ptr%this(i)%model_ptr)) then
          nl=nl+1; write(lines(nl), ffmt)  '%this(', i, ')%Model_ptr:   ', &
                                                            v_ptr%this(i)%model_ptr
        else
          nl=nl+1; write(lines(nl), ffmt)  '%this(', i, ')%Model_ptr:   <not associated>'
        endif
        if (associated (v_ptr%this(i)%base_ptr)) then
          nl=nl+1; write(lines(nl), ffmt)  '%this(', i, ')%Base_ptr:    ', &
                                                            v_ptr%this(i)%base_ptr
        else
          nl=nl+1; write(lines(nl), ffmt)  '%this(', i, ')%Base_ptr:    <not associated>'
        endif
      enddo
    endif

    nl=nl+1; write(lines(nl), fmt)  '%Design_value:   ', v_ptr%design_value
    nl=nl+1; write(lines(nl), fmt)  '%Old_value:      ', v_ptr%old_value        
    nl=nl+1; write(lines(nl), fmt)  '%Data_value:     ', v_ptr%data_value       
    nl=nl+1; write(lines(nl), fmt)  '%Ref_value:      ', v_ptr%ref_value        
    nl=nl+1; write(lines(nl), fmt)  '%correction_alue:', v_ptr%correction_value     
    nl=nl+1; write(lines(nl), fmt)  '%High_lim:       ', v_ptr%high_lim   
    nl=nl+1; write(lines(nl), fmt)  '%Low_lim:        ', v_ptr%low_lim    
    nl=nl+1; write(lines(nl), fmt)  '%Step:           ', v_ptr%step             
    nl=nl+1; write(lines(nl), fmt)  '%Weight:         ', v_ptr%weight           
    nl=nl+1; write(lines(nl), fmt)  '%Ddelta:         ', v_ptr%delta
    nl=nl+1; write(lines(nl), fmt)  '%Merit:          ', v_ptr%merit            
    nl=nl+1; write(lines(nl), fmt)  '%DMerit_dVar:    ', v_ptr%dMerit_dVar      
    nl=nl+1; write(lines(nl), lmt)  '%Exists:         ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%Good_var:       ', v_ptr%good_var  
    nl=nl+1; write(lines(nl), lmt)  '%Good_user:      ', v_ptr%good_user 
    nl=nl+1; write(lines(nl), lmt)  '%Good_opt:       ', v_ptr%good_opt 
    nl=nl+1; write(lines(nl), lmt)  '%Useit_opt:      ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_plot:     ', v_ptr%useit_plot   

! check if there is a variable number
! if no variable number requested, show it all

  else

    write(lines(1), '(2a)') 'Variable name:   ', v1_ptr%name
    lines(2) = ' '
    line1 = '       Name            Data         Model        Design'
    write (lines(3), *) line1
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
    write (lines(nl), *) line1
  endif

! print out results

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "INTERNAL ERROR, SHOULDN'T BE HERE!")
  return

end select

end subroutine tao_show_cmd
