!+
! Subroutine tao_show_cmd (word1, word2, word3)
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

subroutine tao_show_cmd (word1, word2, word3)

use tao_mod
use tao_top10_mod
use tao_single_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_plot_struct), pointer :: plot
type (tao_plot_region_struct), pointer :: region
type (coord_struct) orb

type (ele_struct), pointer :: ele
type (ele_struct) ele3

real(rp) f_phi, s_pos

character(24) :: var_name
character(24)  :: plane, fmt, imt, lmt, amt, ffmt, iimt
character(*) :: word1, word2, word3
character(40) :: show_word1, show_word2, show_word3
character(8) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name
character(80), pointer :: ptr_lines(:)
character(16) ele_name, name

character(16) :: show_names(10) = (/ &
   'data       ', 'var        ', 'global     ', 'alias      ', 'top10      ', &
   'optimizer  ', 'ele        ', 'lattice    ', 'constraints', 'plots      ' /)

character(200), allocatable :: lines(:)
character(100) line1, line2
character(4) null_word

integer :: data_number, ix_plane
integer nl, loc
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju
integer num_locations, max_lines

logical err, found, at_ends
logical show_all, name_found
logical, automatic :: picked(size(s%u))
logical, allocatable :: show_here(:)

max_lines = n_output_lines_maxx 
if (.not. allocated (lines)) allocate (lines(max_lines))
if (size(lines) .ne. max_lines) allocate (lines(max_lines))

null_word = 'null'

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
  enddo
  
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, 'ALL')
  call tao_show_constraints (0, 'TOP10')

!----------------------------------------------------------------------
! data

case ('data')

  call tao_pick_universe (show_word2, show_word2, picked, err)
  if (err) return

  u_loop: do ju = 1, size(s%u)

    if (.not. picked(ju)) cycle
    if (.not. associated (s%u(ju)%d2_data)) return 
    u => s%u(ju)

    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), *) ' '
      nl=nl+1; write(lines(nl), *) 'Universe:', ju
    endif

! If just "show data" then show all names

    if (show_word2 == ' ') then
      nl=nl+1; write (lines(nl), '(62x, a)') 'Bounds' 
      nl=nl+1; write (lines(nl), '(5x, a, 23x, a, 14x, a)') &
                                             'd2_Name', 'Ix  d1_Name', 'Lower  Upper'
      do i = 1, size(u%d2_data)
        d2_ptr => u%d2_data(i)
        if (d2_ptr%name == ' ') cycle
        nl=nl+1; write (lines(nl), '(5x, a)') d2_ptr%name
        do j = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
          d1_ptr => d2_ptr%d1(j)
          name = d1_ptr%name
          if (name == ' ') name = '<blank>'
          nl=nl+1; write (lines(nl), '(32x, i5, 2x, a, 2x, 2i6)') &
                        j, name, lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
        enddo
      enddo
      call out_io (s_blank$, r_name, lines(1:nl))
      cycle u_loop
    endif

! get pointers to the data

    call string_trim (show_word3, show_word3, ix)
    ! are we looking at a range of locations?
    if ((show_word3 .eq. ' ') .or. (index(trim(show_word3), ' ') .ne. 0) &
                                           .or. index(show_word3, ':') .ne. 0) then
      call tao_find_data (err, u, show_word2, d2_ptr, d1_ptr, null_word, d_ptr)
    else
      call tao_find_data (err, u, show_word2, d2_ptr, d1_ptr, show_word3, d_ptr)
    endif
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
      nl=nl+1; write(lines(nl), fmt)  '%meas_value:        ', d_ptr%meas_value
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
      nl=nl+1; write(lines(nl), lmt)  '%Good_meas:         ', d_ptr%good_meas
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
      nl = 2

      line1 = '                                                                  |   Useit'
      line2 = '     Name                     Data         Model        Design    | Opt  Plot'
      nl=nl+1; lines(nl) = line1
      nl=nl+1; lines(nl) = line2

! if a range is specified, show the data range   

      allocate (show_here(lbound(d1_ptr%d,1):ubound(d1_ptr%d,1)))
      if (show_word3 == ' ') then
        show_here = .true.
      else
        call location_decode (show_word3, show_here, lbound(d1_ptr%d,1), num_locations)
        if (num_locations .eq. -1) then
          call out_io (s_error$, r_name, "Syntax error in range list!")
          deallocate(show_here)
          return
        endif
      endif

      do i = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
        if (.not. (show_here(i) .and. d1_ptr%d(i)%exists)) cycle
        if (nl+2 .gt. max_lines) then
          call out_io (s_blank$, r_name, "Too many elements!")
          call out_io (s_blank$, r_name, "Listing first \i5\ selected elements", max_lines-4)
          exit
        endif
        nl=nl+1; write(lines(nl), '(i5, 2x, a16, 3es14.4, 2l6)') i, &
                     d1_ptr%d(i)%name, d1_ptr%d(i)%meas_value, &
                     d1_ptr%d(i)%model_value, d1_ptr%d(i)%design_value, &
                     d1_ptr%d(i)%useit_opt, d1_ptr%d(i)%useit_plot
      enddo

      deallocate(show_here)

      nl=nl+1; lines(nl) = line2
      nl=nl+1; lines(nl) = line1

! else we must have a valid d2_ptr.

    else 

      nl=nl+1; write(lines(nl), '(2a)') 'D2_Data type:    ', d2_ptr%name
      nl=nl+1; write(lines(nl), '(5x, a)') '                   Bounds'
      nl=nl+1; write(lines(nl), '(5x, a)') 'D1_Data name    lower: Upper' 
      do i = 1, size(d2_ptr%d1)
        if (nl+10 .gt. max_lines) then
          call out_io (s_blank$, r_name, &
              "Found too many d1_data! Listing first \i5\ matches", max_lines-1)
          exit
        endif
        nl=nl+1; write(lines(nl), '(5x, a, i5, a, i5)') d2_ptr%d1(i)%name, &
                  lbound(d2_ptr%d1(i)%d, 1), ':', ubound(d2_ptr%d1(i)%d, 1)
      enddo

      if (any(d2_ptr%descrip /= ' ')) then
        nl=nl+1; write (lines(nl), *)
        nl=nl+1; write (lines(nl), '(a)') 'Descrip:'
        do i = 1, size(d2_ptr%descrip)
          if (d2_ptr%descrip(i) /= ' ') then
            nl=nl+1; write (lines(nl), '(i4, 2a)') i, ': ', d2_ptr%descrip(i)
          endif
        enddo
      endif

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
        if (nl+1 .gt. max_lines) then
          call out_io (s_blank$, r_name, "Too many elements!")
          call out_io (s_blank$, r_name, "Listing first \i5\ selected elements", max_lines-1)
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

    ! Show the element info
    call type2_ele (u%model%ele_(loc), ptr_lines, n, .true., 6, .false., &
                                        s%global%phase_units, .true., u%model)
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

    ! Show data associated with this element
    call show_ele_data (u, loc, lines, nl)

    found = .false.
    do i = loc + 1, u%model%n_ele_max
      if (u%model%ele_(i)%name == ele_name) then
        if (nl+2 .gt. max_lines) then
          nl = nl + 1
          write (lines(nl), *) "Found too many elements!"
          exit
        endif
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

  nl=nl+1; write (lines(nl), imt) 'n_universes:       ', size(s%u)
  nl=nl+1; write (lines(nl), imt) 'u_view:            ', s%global%u_view
  nl=nl+1; write (lines(nl), imt) 'phase_units:       ', s%global%phase_units
  nl=nl+1; write (lines(nl), imt) 'n_opti_cycles:     ', s%global%n_opti_cycles
  nl=nl+1; write (lines(nl), amt) 'track_type:        ', s%global%track_type
  if (s%global%track_type .eq. 'macro') &
  nl=nl+1; write (lines(nl), imt) 'bunch_to_plot::    ', s%global%bunch_to_plot
  nl=nl+1; write (lines(nl), amt) 'optimizer:         ', s%global%optimizer
  nl=nl+1; write (lines(nl), amt) 'prompt_string:     ', s%global%prompt_string
  nl=nl+1; write (lines(nl), amt) 'var_out_file:      ', s%global%var_out_file
  nl=nl+1; write (lines(nl), amt) 'opt_var_out_file:  ', s%global%opt_var_out_file
  nl=nl+1; write (lines(nl), amt) 'print_command:     ', s%global%print_command
  nl=nl+1; write (lines(nl), amt) 'current_init_file: ',s%global%current_init_file
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on:     ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), lmt) 'opt_with_ref:      ', s%global%opt_with_ref 
  nl=nl+1; write (lines(nl), lmt) 'opt_with_base:     ', s%global%opt_with_base
  nl=nl+1; write (lines(nl), lmt) 'plot_on:           ', s%global%plot_on
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on:     ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), amt) 'curren_init_file:  ', s%global%current_init_file

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! lattice

case ('lattice')
  
  if (show_word2 .eq. ' ') then
    nl=nl+1
    write (lines(nl), '(a, i2, a, i6, a)') "Universe ", s%global%u_view, &
               " has ", u%model%n_ele_use, " regular elements."
    nl=nl+1
    if (u%is_on) then
      write (lines(nl), '(a)') "This universe is turned ON"
    else
      write (lines(nl), '(a)') "This universe is turned OFF"
    endif
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif
  
  allocate (show_here(0:u%model%n_ele_use))
  show_word3 = trim(show_word2) // trim(show_word3)
  call location_decode (show_word3, show_here, 0, num_locations)
  if (num_locations .eq. -1) then
    call out_io (s_error$, r_name, "Syntax error in range list!")
    deallocate(show_here)
    return
  endif

  if (.true.) then
    at_ends = .true.
    write (lines(nl+1), '(37x, a)') 'Model values at End of Element:'
  else
    at_ends = .false.
    write (lines(nl+1), '(37x, a)') 'Model values at Center of Element:'
  endif


  write (lines(nl+2), '(29x, 22x, a)') &
                     '|              X           |             Y        '
  write (lines(nl+3), '(6x, a, 16x, a)') ' Name             key', &
                  '   S    |  Beta   Phi   Eta  Orb   | Beta    Phi    Eta   Orb'

  nl=nl+3
  do ix = lbound(show_here,1), ubound(show_here,1)
    if (.not. show_here(ix)) cycle
    if (nl+2 .gt. max_lines) then
      call out_io (s_blank$, r_name, "Too many elements!")
      call out_io (s_blank$, r_name, "Listing first \i5\ selected elements", max_lines-4)
      exit
    endif
    ele => u%model%ele_(ix)
    if (ix == 0 .or. at_ends) then
      ele3 = ele
      orb = u%model_orb(ix)
      s_pos = ele3%s
    else
      call twiss_and_track_partial (u%model%ele_(ix-1), ele, &
                u%model%param, ele%value(l$)/2, ele3, u%model_orb(ix-1), orb)
      s_pos = ele%s-ele%value(l$)/2
    endif
    nl=nl+1
    write (lines(nl), '(i6, 1x, a16, 1x, a16, f10.3, 2(f7.2, f8.3, f5.1, f8.3))') &
          ix, ele%name, key_name(ele%key), s_pos, &
          ele3%x%beta, f_phi*ele3%x%phi, ele3%x%eta, 1000*orb%vec(1), &
          ele3%y%beta, f_phi*ele3%y%phi, ele3%y%eta, 1000*orb%vec(3)
  enddo

  write (lines(nl+1), '(6x, a, 16x, a)') ' Name             key', &
                  '   S    |  Beta   Phi   Eta  Orb   | Beta    Phi    Eta   Orb'
  write (lines(nl+2), '(29x, 22x, a)') &
                     '|              X           |             Y        '
  nl=nl+2
  
  call out_io (s_blank$, r_name, lines(1:nl))

  deallocate(show_here)

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

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; write (lines(nl), amt) 'optimizer:        ', s%global%optimizer
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! plots

case ('plots')

  nl=nl+1; lines(nl) = '   '
  nl=nl+1; lines(nl) = 'Template Plots: Graphs'
  do i = 1, size(s%template_plot)
    plot => s%template_plot(i)
    if (plot%name == ' ') cycle
    nl=nl+1; write (lines(nl), '(4x, 2a)') trim(plot%name), ' :'
    if (associated(plot%graph)) then
      do j = 1, size(plot%graph)
        nl=nl+1; write (lines(nl), '(20x, a)') plot%graph(j)%name
      enddo
    endif
  enddo

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; lines(nl) = '[Visible]     Plot Region     <-->  Template' 
  nl=nl+1; lines(nl) = '---------     -----------           ------------'
  do i = 1, size(s%plot_page%region)
    region => s%plot_page%region(i)
    nl=nl+1; write (lines(nl), '(3x l1, 10x, 3a)') region%visible, &
                                    region%name, '<-->  ', region%plot%name
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))

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
    write (lines(1), '(5x, a)') '                  Bounds'
    write (lines(2), '(5x, a)') 'Name            Lower  Upper'
    nl = 2
    do i = 1, size(s%v1_var)
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      if (nl+1 .gt. max_lines) then
        call out_io (s_blank$, r_name, &
          & "Found too many v1_vars! Listing first \i5\ matches", max_lines)
        exit
      endif
      nl=nl+1
      write(lines(nl), '(5x, a, i5, i7)') v1_ptr%name, &
                  lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the variables

  call string_trim (show_word3, show_word3, ix)
! are we looking at a range of locations?
  if ((show_word3 .eq. ' ') .or. (index(trim(show_word3), ' ') .ne. 0) &
                                       .or. index(show_word3, ':') .ne. 0) then
    call tao_find_var(err, show_word2, v1_ptr, null_word, v_ptr) 
  else
    call tao_find_var(err, show_word2, v1_ptr, show_word3, v_ptr) 
  endif
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

    nl=nl+1; write(lines(nl), fmt)  '%Design_value:    ', v_ptr%design_value
    nl=nl+1; write(lines(nl), fmt)  '%Old_value:       ', v_ptr%old_value
    nl=nl+1; write(lines(nl), fmt)  '%Meas_value:      ', v_ptr%meas_value
    nl=nl+1; write(lines(nl), fmt)  '%Ref_value:       ', v_ptr%ref_value
    nl=nl+1; write(lines(nl), fmt)  '%Correction_value:', v_ptr%correction_value
    nl=nl+1; write(lines(nl), fmt)  '%High_lim:        ', v_ptr%high_lim
    nl=nl+1; write(lines(nl), fmt)  '%Low_lim:         ', v_ptr%low_lim
    nl=nl+1; write(lines(nl), fmt)  '%Step:            ', v_ptr%step
    nl=nl+1; write(lines(nl), fmt)  '%Weight:          ', v_ptr%weight
    nl=nl+1; write(lines(nl), fmt)  '%Delta:           ', v_ptr%delta
    nl=nl+1; write(lines(nl), amt)  '%Merit_type:      ', v_ptr%merit_type
    nl=nl+1; write(lines(nl), fmt)  '%Merit:           ', v_ptr%merit
    nl=nl+1; write(lines(nl), fmt)  '%dMerit_dVar:     ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), lmt)  '%Exists:          ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%Good_var:        ', v_ptr%good_var
    nl=nl+1; write(lines(nl), lmt)  '%Good_user:       ', v_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%Good_opt:        ', v_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_opt:       ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_plot:      ', v_ptr%useit_plot

! check if there is a variable number
! if no variable number requested, show a range

  else

    write(lines(1), '(2a)') 'Variable name:   ', v1_ptr%name
    lines(2) = ' '
    line1 = '       Name                     Data         Model        Design  Useit_opt'
    write (lines(3), *) line1
    nl = 3
    ! if a range is specified, show the variable range   
    if (show_word3 .ne. ' ') then
      allocate (show_here(lbound(v1_ptr%v,1):ubound(v1_ptr%v,1)))
      call location_decode (show_word3, show_here, lbound(v1_ptr%v,1), num_locations)
      if (num_locations .eq. -1) then
        call out_io (s_error$, r_name, "Syntax error in range list!")
        deallocate(show_here)
        return
      endif
      do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
        if (.not. (show_here(i) .and. v1_ptr%v(i)%exists)) cycle
        if (nl+2 .gt. max_lines) then
          call out_io (s_blank$, r_name, "Too many elements!")
          call out_io (s_blank$, r_name, "Listing first \i5\ selected elements", max_lines-4)
          exit
        endif
        nl=nl+1
        write(lines(nl), '(i6, 2x, a16, 3es14.4, 7x, l)') i, &
                 v1_ptr%v(i)%name, v1_ptr%v(i)%meas_value, &
                 v1_ptr%v(i)%model_value, v1_ptr%v(i)%design_value, v1_ptr%v(i)%useit_opt
      enddo 
      nl=nl+1
      write (lines(nl), *) line1
      deallocate(show_here)
    else
      do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
        if (.not. v1_ptr%v(i)%exists) cycle
        if (nl+2 .gt. max_lines) then
          call out_io (s_blank$, r_name, "Too many elements!")
          call out_io (s_blank$, r_name, "Listing first \i5\ elements", max_lines-4)
          exit
        endif
        nl=nl+1
        write(lines(nl), '(i6, 2x, a16, 3es14.4, 7x, l)') i, &
                 v1_ptr%v(i)%name, v1_ptr%v(i)%meas_value, &
                 v1_ptr%v(i)%model_value, v1_ptr%v(i)%design_value, v1_ptr%v(i)%useit_opt
      enddo
      nl=nl+1
      write (lines(nl), *) line1
    endif
  endif

! print out results

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "INTERNAL ERROR, SHOULDN'T BE HERE!")
  return

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine show_ele_data (u, i_ele, lines, nl)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum
character(*) :: lines(:)
integer i_ele, nl, i

character(30) :: dmt = "(a, 3(1x, es15.5)) "

logical :: found_one = .false.

  nl=nl+1; write (lines(nl), '(a)') "  "
  nl=nl+1; write (lines(nl), '(a)') &
        "   Data Type      |  Model Value  |  Design Value |  Base Value"

  do i = 1, size(u%data)
    if (u%data(i)%ix_ele .eq. i_ele) then
      found_one = .true.
      datum => u%data(i)
      nl = nl + 1
      write (lines(nl), dmt) datum%data_type, datum%model_value, &
                             datum%design_value, datum%base_value 
    endif
  enddo

  if (.not. found_one) then
    nl = nl +1 
    write (lines(nl), '(a)') "No data types associated with this element."
  endif

  nl=nl+1; write (lines(nl), '(a)') "  "
  nl=nl+1; write (lines(nl), '(a)') &
        "   Data Type      |  Model Value  |  Design Value |  Base Value"


end subroutine show_ele_data

end subroutine tao_show_cmd
