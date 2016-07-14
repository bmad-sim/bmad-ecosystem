!+
! Subroutine tao_python_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
! All lines are saved in linex
!
! Note: The syntax for "serial output form" is:
!   <component_name>;<type>;<variable>;<component_value>;
! <type> is one of:
!   STRING
!   INTEGER
!   REAL
!   LOGICAL
! <variable> indicates if the component can be varied. It is one of:
!   T
!   F
!
! Input:
!   input_str  -- Character(*): What to show.
!-


subroutine tao_python_cmd (input_str)

use tao_mod
use tao_command_mod
use location_encode_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, save, target :: v_array(:)
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: c1
type (tao_plot_region_struct), pointer :: region
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (wake_lr_struct), pointer :: lr
type (ele_struct), pointer :: ele
type (coord_struct), target :: orb
type (ele_struct), target :: ele3
type (bunch_params_struct) bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch
type (random_state_struct) ran_state
type (tao_scratch_space_struct), pointer :: ss

character(*) input_str
character(24) imt, lmt, amt, iamt, ramt, f3mt, rmt, irmt, iimt
character(40) max_loc, loc_ele, name1(40), name2(40)
character(200) line, file_name
character(20) cmd, command
character(20) :: r_name = 'tao_python_cmd'
character(20) :: cmd_names(24)= &
          ['plot_visible   ', 'plot_template  ', 'graph          ', 'curve          ', &
           'plot1          ', 'var_all        ', 'var_v1         ', 'var1           ', &
           'help           ', 'curve1         ', 'curve_sym      ', 'curve_line     ', &
           'lat_global     ', 'data_d2        ', 'data_d1        ', 'data1          ', &
           'ele_all        ', 'ele1_all       ', 'beam_all       ', 'ele1_attrib    ', &
           'constraint_data', 'constraint_vars', 'global         ', 'lat_ele_list   ']

real(rp) target_value, angle

integer :: i, j, ie, iu, md, nl, ct, n1, nl2, n, ix, iu_write
integer :: ix_ele, ix_ele1, ix_ele2, ix_branch, ix_universe
integer :: ios
logical :: err, print_flag, opened, doprint

character(20) switch

!

line = input_str
doprint = .true.
opened = .false.

do
  call tao_next_switch (line, ['-append ', '-write  ', '-noprint'], .false., switch, err, ix)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-noprint')
    doprint = .false.

  case ('-append', '-write')
    call string_trim(line, line, ix)
    file_name = line(:ix)
    call string_trim(line(ix+1:), line, ix)

    iu_write = lunget()
    if (switch == '-append') then
      open (iu_write, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = 200)
    else
      open (iu_write, file = file_name, status = 'REPLACE', recl = 200)
    endif

    opened = .true.
  end select
enddo

call string_trim(line, line, ix)
cmd = line(1:ix)
call string_trim(line(ix+1:), line, ix)

call match_word (cmd, cmd_names, ix, matched_name = command)
if (ix == 0) then
  call out_io (s_error$, r_name, '***PYTHON WHAT? WORD NOT RECOGNIZED: ' // command)
  return
endif

if (ix < 0) then
  call out_io (s_error$, r_name, '***PYTHON COMMAND? AMBIGUOUS: ' // command)
  return
endif

amt = '(3a)'
imt = '(a,i0,a)'
rmt = '(a,es15.8,a)'
lmt = '(a,l1,a)'

nl = 0
ss => scratch
if (.not. allocated(scratch%lines)) allocate (scratch%lines(200))

select case (command)

!----------------------------------------------------------------------
! help
! returns list of "help xxx" topics

case ('help')

  call tao_help ('help-list', '', ss%lines, n)

  nl2 = 0
  do i = 1, n
    if (ss%lines(i) == '') cycle
    call string_trim(ss%lines(i), line, ix)
    nl=nl+1; name1(nl) = line(1:ix)
    call string_trim(line(ix+1:), line, ix)
    if (ix == 0) cycle
    nl2=nl2+1; name2(nl2) = line
  enddo

  ss%lines(1:nl) = name1(1:nl)
  ss%lines(nl+1:nl+nl2) = name2(1:nl2)
  nl = nl + nl2

!----------------------------------------------------------------------
! data used in optimizations (with datum%useit_opt = T). 
! Input syntax: 
!   python constrint_data
! Output syntax:
!   <datum_name>;<constraint_type>;<ref_ele>;<start_ele>;<ele>;<meas_value>;<model_value>;<merit_value>;<max_merit_location>
! Example output line:
!   orbit.x[0];orbit.x <target>;;;DET_00W;  0.0000E+00;  1.4693E-03;  2.1589E+00;;


case ('constraint_data')

do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, size(s%u(i)%data)
    d_ptr => s%u(i)%data(j)
    if (.not. d_ptr%useit_opt) cycle
    
    branch => s%u(i)%model%lat%branch(d_ptr%ix_branch)
    ie = d_ptr%ix_ele_merit

    if (ie < 0) then
      max_loc = ''
    else
      max_loc = branch%ele(ie)%name
    endif

    if (nl == size(ss%lines)) call re_allocate(ss%lines, 2*nl)
    nl=nl+1; write (ss%lines(nl), '(10a, 3(es12.4, a), 2a)') &
        trim(tao_datum_name(d_ptr)), ';', trim(tao_constraint_type_name(d_ptr)), ';', &
        trim(d_ptr%ele_ref_name), ';', trim(d_ptr%ele_start_name), ';', trim(d_ptr%ele_name), ';', &
        d_ptr%meas_value, ';', d_ptr%model_value, ';', d_ptr%merit, ';', trim(max_loc), ';'

  enddo
enddo

!----------------------------------------------------------------------
! Variables used in optimizations (with var%useit_opt = T). 
! Input syntax:  
!   python constrint_vars
! Output syntax:
!   <var_name>;<attribute_name>;<lattice_element>;<target_value>;<model_value>;<merit_value>

case ('constraint_vars')

do i = 1, s%n_var_used
  v_ptr => s%var(i)
  if (.not. v_ptr%useit_opt) cycle

  u => s%u(v_ptr%this(1)%ix_uni)

  loc_ele = ''

  if (v_ptr%this(1)%ix_ele >= 0) then
    branch => u%model%lat%branch(v_ptr%this(1)%ix_branch)
    ct = branch%ele(v_ptr%this(1)%ix_ele)%lord_status

    if (ct /= group_lord$ .and. ct /= overlay_lord$ .and. ct /= multipass_lord$) then
      write (loc_ele, '(f8.2)') branch%ele(v_ptr%this(1)%ix_ele)%s
      call string_trim (loc_ele, loc_ele, ix)
    endif
  endif

  if (v_ptr%merit_type == 'target') then
    target_value = v_ptr%meas_value
  elseif (v_ptr%merit_type == 'limit') then
    if (abs(v_ptr%model_value - v_ptr%high_lim) < abs(v_ptr%model_value - v_ptr%low_lim)) then
      target_value = v_ptr%high_lim
    else
      target_value = v_ptr%low_lim
    endif
  endif

  if (nl == size(ss%lines)) call re_allocate (ss%lines, 2*nl)
  nl=nl+1; write (ss%lines(nl), '(6a, 3(es12.4, a))') &
        trim(tao_var1_name(v_ptr)), ';', trim(tao_var_attrib_name(v_ptr)), ';', &
        trim(loc_ele), ';', target_value, ';', v_ptr%model_value, ';', v_ptr%merit, ';'

enddo

!----------------------------------------------------------------------
! Curve information for a plot
! Input syntax:
!   pyton curve <curve_name>
! Output syntax is serial output form. See documentation at beginning of this file.

case ('curve')

  call tao_find_plots (err, line, 'BOTH', curve = curve, always_allocate = .true.)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=nl+1; ss%lines(nl) = 'INVALID;;'
  else
    c1 => curve(1)%c
    nl=nl+1; write (ss%lines(nl), amt)  'data_source;STRING;F;',               trim(c1%data_source), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'data_index;STRING;F;',                trim(c1%data_index), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'data_type_x;STRING;F;',               trim(c1%data_type_x), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'data_type_z;STRING;F;',               trim(c1%data_type_z), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'data_type;STRING;F;',                 trim(c1%data_type), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'legend_text;STRING;F;',               trim(c1%legend_text), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'ele_ref_name;STRING;F;',              trim(c1%ele_ref_name), ';'
    nl=nl+1; write (ss%lines(nl), imt)  'ix_branch;INTEGER;F;',                c1%ix_branch, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'ix_ele_ref;INTEGER;F;',               c1%ix_ele_ref, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'ix_ele_ref_track;INTEGER;F;',         c1%ix_ele_ref_track, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'ix_bunch;INTEGER;F;',                 c1%ix_bunch, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'ix_universe;INTEGER;F;',              c1%ix_universe, ';'
    nl=nl+1; write (ss%lines(nl), rmt)  'z_color0;REAL;F;',                    c1%z_color0, ';'
    nl=nl+1; write (ss%lines(nl), rmt)  'z_color1;REAL;F;',                    c1%z_color1, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'symbol_every;INTEGER;F;',             c1%symbol_every, ';'
    nl=nl+1; write (ss%lines(nl), rmt)  'y_axis_scale_factor;REAL;F;',         c1%y_axis_scale_factor, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'use_y2;LOGICAL;F;',                   c1%use_y2, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'use_z_color;LOGICAL;F;',              c1%use_z_color, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'draw_line;LOGICAL;F;',                c1%draw_line, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'draw_symbols;LOGICAL;F;',             c1%draw_symbols, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'draw_symbol_index;LOGICAL;F;',        c1%draw_symbol_index, ';'
    nl=nl+1; write (ss%lines(nl), lmt)  'smooth_line_calc;LOGICAL;F;',         c1%smooth_line_calc, ';'
    nl=nl+1; write (ss%lines(nl), imt)  'line%width;INTEGER;F;',               c1%line%width, ';'
    nl=nl+1; write (ss%lines(nl), amt)  'line%color;STRING;F;',                trim(qp_color_name(c1%line%color)), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'line%pattern;STRING;F;',              trim(qp_line_pattern_name(c1%line%pattern)), ';'
    nl=nl+1; write (ss%lines(nl), amt)  'symbol%type;STRING;F;',               trim(qp_symbol_type_name(c1%symbol%type)), ';'
    nl=nl+1; write (ss%lines(nl), rmt)  'symbol%height;REAL;F;',               c1%symbol%height, ';'
    nl=nl+1; write (ss%lines(nl), amt)  'symbol%fill_pattern;STRING;F;',       trim(qp_fill_name(c1%symbol%fill_pattern)), ';'
    nl=nl+1; write (ss%lines(nl), imt)  'symbol%line_width;INTEGER;F;',        c1%symbol%line_width, ';'
  endif

!----------------------------------------------------------------------
! Points used to construct a smooth line for a plot curve.
! Input syntax:
!   python curve <curve>
! Output syntax: 
!   <index>;<x>;<y>;

case ('curve_line')

  call tao_find_plots (err, line, 'BOTH', curve = curve, always_allocate = .true.)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=nl+1; ss%lines(nl) = 'INVALID;;;'
  else
    c1 => curve(1)%c
    call re_allocate (ss%lines, nl+size(c1%x_line)+100, .false.)
    do i = 1, size(c1%x_line)
      nl=nl+1; write (ss%lines(nl), '(i0, a, 2(es15.8, a))') i, ';', c1%x_line(i), ';', c1%y_line(i), ';'
    enddo
  endif

!----------------------------------------------------------------------
! Locations to draw symbols for a plot curve.
! Input syntax:
!   python curve <curve>
! Output syntax: 
!   <index>;<symbol_index>;<x>;<y>;

case ('curve_sym')

  call tao_find_plots (err, line, 'BOTH', curve = curve, always_allocate = .true.)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=nl+1; ss%lines(nl) = 'INVALID;;;'
  else
    c1 => curve(1)%c
    call re_allocate (ss%lines, nl+size(c1%x_symb)+100, .false.)
    do i = 1, size(c1%x_symb)
      nl=nl+1; write (ss%lines(nl), '(2(i0, a), 2(es15.8, a))') i, ';', c1%ix_symb(i), ';', c1%x_symb(i), ';', c1%y_symb(i)
    enddo
  endif

!----------------------------------------------------------------------
! List of datums in a given data d1 array.
! Use the "python data_d2" command to get a list of d2.d1 arrays. 
! Use the "python data1" command to get detailed information on a particular datum.
! Input syntax:
!   python data_d1 <d1_datum>
! Output syntax: 
!   <index>;<data_type>;<merit_type>;<ele_ref_name>;<ele_start_name>;<ele_name>;<meas_value>;<model_value>;<design_value>;<good_user>;<useit_opt>;<useit_plot>;

case ('data_d1')

  call tao_find_data (err, line, d1_array = d1_array)

  if (.not. allocated(d1_array) .or. size(d1_array) /= 1) then
    nl=nl+1; ss%lines(nl) = '0;INVALID;;;;0;0;0;F;F;'
  else
    d1_ptr => d1_array(1)%d1
    do i = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
      d_ptr => d1_ptr%d(i)
      if (.not. d_ptr%exists) cycle
      if (nl == size(ss%lines)) call re_allocate (ss%lines, nl+200, .false.)
!      nl=nl+1; write(ss%lines(nl), '(i0, 11a, 3(es15.8, a), 3(l1, a))') &
!                       d_ptr%ix_d1, ';', trim(d_ptr%data_type), ';', trim(d_ptr%merit_type), ';', &
!                       trim(d_ptr%ele_ref_name), ';', trim(d_ptr%ele_start_name), ';', &
!                       trim(d_ptr%ele_name), ';', d_ptr%meas_value, ';', d_ptr%model_value, ';', &
!                       d_ptr%design_value, ';', d_ptr%good_user, ';', d_ptr%useit_opt, ';', d_ptr%useit_plot, ';'
!      nl=nl+1; write(ss%lines(nl), '( (a, i0,a), 15a, 3(a, es15.8, a), 9a)') &
!                       "ix_d1=",           d_ptr%ix_d1, ";", &
!                       "data_type='",      trim(d_ptr%data_type), "';", &
!                       "merit_type='",     trim(d_ptr%merit_type), "';", &
!                       "ele_ref_name='",   trim(d_ptr%ele_ref_name), "';", &
!                       "ele_start_name='", trim(d_ptr%ele_start_name), "';", &
!                       "ele_name='",       trim(d_ptr%ele_name), "';", &
!                       "meas_value=",      d_ptr%meas_value, ";", &
!                       "model_value=",     d_ptr%model_value, ";", &
!                       "design_value=",    d_ptr%design_value, ";", &
!                       "good_user=",       py_bool(d_ptr%good_user), ";", &
!                       "useit_opt=",       py_bool(d_ptr%useit_opt), ";", &
!                       "useit_plot=",      py_bool(d_ptr%useit_plot), ";"                         
      nl=nl+1; write(ss%lines(nl), '(i0, 11a, 3(es15.8, a), 6a)') &
                       d_ptr%ix_d1, ';', trim(d_ptr%data_type), ';', &
                       py_string(d_ptr%merit_type),     ';', py_string(d_ptr%ele_ref_name), ';', &
                       py_string(d_ptr%ele_start_name), ';', py_string(d_ptr%ele_name), ';', &
                       d_ptr%meas_value, ';', d_ptr%model_value, ';', &
                       d_ptr%design_value, ';',  py_bool(d_ptr%good_user), ';',  &
                       py_bool(d_ptr%useit_opt), ';',  py_bool(d_ptr%useit_plot)
    enddo
  endif

!----------------------------------------------------------------------
! D2 level data. 
! Use the "python data-d1" command to get detailed info on a specific d1 array.
! Input syntax:
!   python data_d2
! Output syntax: 
!   <ix_universe>;<d2_name>;<d1_name>;<lbound_d_array>;<ubound_d_array>;<useit_list>;

case ('data_d2')

  do iu = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(iu)

    do i = 1, u%n_d2_data_used
      d2_ptr => u%d2_data(i)
      if (d2_ptr%name == ' ') cycle

      do j = 1, size(d2_ptr%d1)
        d1_ptr => d2_ptr%d1(j)

        call location_encode (line, d1_ptr%d%useit_opt, d1_ptr%d%exists, lbound(d1_ptr%d, 1))

        nl=nl+1; write (ss%lines(nl), '(i0, 5a, i0, a, i0, 3a)') iu, ';', trim(d2_ptr%name), ';', &
              trim(d1_ptr%name), ';', lbound(d1_ptr%d, 1), ';', ubound(d1_ptr%d, 1), ';', trim(line), ';'

      enddo

    enddo

  enddo

!----------------------------------------------------------------------
! Individual datum info.
! Use the "python data-d1" command to get detailed info on a specific d1 array.
! Input syntax:
!   python data1 <datum_name>
! Output syntax is serial output form. See documentation at beginning of this file.

case ('data1')

  call tao_find_data (err, line, d_array = d_array)

  if (.not. allocated(d_array) .or. size(d_array) /= 1) then
    nl=nl+1; ss%lines(nl) = 'INVALID;STRING;F;0;'
  else
    d_ptr => d_array(1)%d
    nl=nl+1; write(ss%lines(nl), amt)  'ele_name;STRING;F;',           trim(d_ptr%ele_name), ';'
    nl=nl+1; write(ss%lines(nl), amt)  'ele_start_name;STRING;F;',     trim(d_ptr%ele_start_name), ';'
    nl=nl+1; write(ss%lines(nl), amt)  'ele_ref_name;STRING;F;',       trim(d_ptr%ele_ref_name), ';'
    nl=nl+1; write(ss%lines(nl), amt)  'data_type;STRING;F;',          trim(d_ptr%data_type), ';'
    nl=nl+1; write(ss%lines(nl), amt)  'data_source;STRING;F;',        trim(d_ptr%data_source), ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_branch;INTEGER;F;',         d_ptr%ix_branch, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_ele;INTEGER;F;',            d_ptr%ix_ele, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_ele_start;INTEGER;F;',      d_ptr%ix_ele_start, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_ele_ref;INTEGER;F;',        d_ptr%ix_ele_ref, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_ele_merit;INTEGER;F;',      d_ptr%ix_ele_merit, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_dmodel;INTEGER;F;',         d_ptr%ix_dModel, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_d1;INTEGER;F;',             d_ptr%ix_d1, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_data;INTEGER;F;',           d_ptr%ix_data, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_bunch;INTEGER;F;',          d_ptr%ix_bunch, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'model;REAL;F;',                d_ptr%model_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'design;REAL;F;',               d_ptr%design_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'meas;REAL;T;',                 d_ptr%meas_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'ref;REAL;T;',                  d_ptr%ref_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'base;REAL;F;',                 d_ptr%base_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'old;REAL;F;',                  d_ptr%old_value   , ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'invalid;REAL;F;',              d_ptr%invalid_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  's;REAL;F;',                    d_ptr%s, ';'
    nl=nl+1; write(ss%lines(nl), amt)  'merit_type;STRING;F;',         trim(d_ptr%merit_type), ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'merit;REAL;F;',                d_ptr%merit, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'delta_merit;REAL;F;',          d_ptr%delta_merit, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'weight;REAL;F;',               d_ptr%weight, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'exists;LOGICAL;F;',            d_ptr%exists, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_model;LOGICAL;F;',        d_ptr%good_model, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_design;LOGICAL;F;',       d_ptr%good_design, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_base;LOGICAL;F;',         d_ptr%good_base , ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_meas;LOGICAL;F;',         d_ptr%good_meas, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_ref;LOGICAL;F;',          d_ptr%good_ref, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_user;LOGICAL;T;',         d_ptr%good_user, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_opt;LOGICAL;F;',          d_ptr%good_opt, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_plot;LOGICAL;F;',         d_ptr%good_plot, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'useit_plot;LOGICAL;F;',        d_ptr%useit_plot, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'useit_opt;LOGICAL;F;',         d_ptr%useit_opt, ';'
  endif


!----------------------------------------------------------------------
! All parameters associated with given element. 
! Input syntax: 
!   python ele1_all <ix_ele> <ix_branch> <ix_universe>

case ('ele1_all')
  read (line, *,  iostat = ios) ix_ele, ix_branch, ix_universe
  
  if (ios /= 0) then
    call out_io (s_error$, r_name, "Correct input: ele1_all <ix_ele> <ix_branch> <ix_universe>")
    return
  endif
  
  ele => s%u(ix_universe)%model%lat%branch(ix_branch)%ele(ix_ele)
  
  if (ele%key == sbend$ .or. ele%key == rbend$) then
    angle = ele%value(angle$)
  else 
    angle = 0.0_rp
  endif
  
  if (nl == size(ss%lines)) call re_allocate (ss%lines, nl+200, .false.)
  nl=nl+1; write(ss%lines(nl), '(i0, 4a, 2(a, es15.8))') &
    ele%ix_ele, ';', trim(ele%name), ';', trim(key_name(ele%key)),  ';', ele%value(L$), ';', angle
  nl=nl+1; write(ss%lines(nl), '( 5(es15.8,a), es15.8)') &
    ele%floor%r(1), ';', &
    ele%floor%r(2), ';', &
    ele%floor%r(3), ';', &
    ele%floor%theta, ';', &
    ele%floor%phi, ';', &
    ele%floor%psi


!----------------------------------------------------------------------
! Attribute list. Shows which attributes are adjustable.
! Input syntax: 
!   python ele1_all <ix_ele> <ix_branch> <ix_universe>

case ('ele1_attrib')


!----------------------------------------------------------------------
! Global parameters
! Input syntax: 
!   python global

case ('global')

!----------------------------------------------------------------------
! Lattice element list.
! Input syntax:
!   python lat_ele <ix_ele_start> <ix_ele_end> <ix_branch>
! Output syntax is serial output form. See documentation at beginning of this file.

case ('graph')

    call tao_find_plots (err, line, 'BOTH', graph = graph, always_allocate = .true.)

    if (.not. allocated(graph) .and. size(graph) /= 1) then
      nl=nl+1; ss%lines(nl) = 'INVALID;STRING;F;0;'
    else
      g => graph(1)%g
      nl=nl+1; write (ss%lines(nl), amt) 'type;STRING;F;',                       trim(g%type), ';'
      nl=nl+1; write (ss%lines(nl), amt) 'title;STRING;F;',                      trim(g%title), ';'
      nl=nl+1; write (ss%lines(nl), amt) 'title_suffix;STRING;F;',               trim(g%title_suffix), ';'
      nl=nl+1; write (ss%lines(nl), amt) 'component;STRING;F;',                  trim(g%component), ';'
      nl=nl+1; write (ss%lines(nl), imt) 'box1;INTEGER;F;',                      g%box(1), ';'
      nl=nl+1; write (ss%lines(nl), imt) 'box2;INTEGER;F;',                      g%box(2), ';'
      nl=nl+1; write (ss%lines(nl), imt) 'box3;INTEGER;F;',                      g%box(3), ';'
      nl=nl+1; write (ss%lines(nl), imt) 'box4;INTEGER;F;',                      g%box(4), ';'
      nl=nl+1; write (ss%lines(nl), imt) 'ix_universe;INTEGER;F;',               g%ix_universe, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'valid;LOGICAL;F;',                     g%valid, ';'

      nl=nl+1; write (ss%lines(nl), rmt) 'x_axis_scale_factor;REAL;F;',          g%x_axis_scale_factor, ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'symbol_size_scale;REAL;F;',            g%symbol_size_scale, ';'
      nl=nl+1; write (ss%lines(nl), amt) 'x%label;STRING;F;',                    trim(g%x%label), ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'x%max;REAL;F;',                        g%x%max, ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'x%min;REAL;F;',                        g%x%min, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'x%major_div;INTEGER;F;',               g%x%major_div, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'x%major_div_nominal;INTEGER;F;',       g%x%major_div_nominal, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'x%places;INTEGER;F;',                  g%x%places, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'x%draw_label;LOGICAL;F;',              g%x%draw_label, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'x%draw_numbers;LOGICAL;F;',            g%x%draw_numbers, ';'

      nl=nl+1; write (ss%lines(nl), lmt) 'y2_mirrors_y;LOGICAL;F;',              g%y2_mirrors_y, ';'
      nl=nl+1; write (ss%lines(nl), amt) 'y%label;STRING;F;',                    trim(g%y%label), ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'y%max;REAL;F;',                        g%y%max, ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'y%min;REAL;F;',                        g%y%min, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y%major_div;INTEGER;F;',               g%y%major_div, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y%major_div_nominal;INTEGER;F;',       g%y%major_div_nominal, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y%places;INTEGER;F;',                  g%y%places, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'y%draw_label;LOGICAL;F;',              g%y%draw_label, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'y%draw_numbers;LOGICAL;F;',            g%y%draw_numbers, ';'

      nl=nl+1; write (ss%lines(nl), amt) 'y2%label;STRING;F;',                   trim(g%y2%label), ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'y2%max;REAL;F;',                       g%y2%max, ';'
      nl=nl+1; write (ss%lines(nl), rmt) 'y2%min;REAL;F;',                       g%y2%min, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y2%major_div;INTEGER;F;',              g%y2%major_div, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y2%major_div_nominal;INTEGER;F;',      g%y2%major_div_nominal, ';'
      nl=nl+1; write (ss%lines(nl), imt) 'y2%places;INTEGER;F;',                 g%y2%places, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'y2%draw_label;LOGICAL;F;',             g%y2%draw_label, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'y2%draw_numbers;LOGICAL;F;',           g%y2%draw_numbers, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'limited;LOGICAL;F;',                   g%limited, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'clip;LOGICAL;F;',                      g%clip, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'draw_axes;LOGICAL;F;',                 g%draw_axes, ';'
      nl=nl+1; write (ss%lines(nl), lmt) 'correct_xy_distortion;LOGICAL;F;',     g%correct_xy_distortion, ';'
      do i = 1, size(g%curve)
        nl=nl+1; write (ss%lines(nl), '(a, i0, 3a)') 'curve(', i, ');STRING;F;',   trim(g%curve(i)%name), ';'
      enddo
    endif

!----------------------------------------------------------------------
! Lattice element list.
! Input syntax:
!   python lat_ele <ix_ele_start> <ix_ele_end> <ix_branch>

case ('lat_ele_list')



!----------------------------------------------------------------------
! Lattice globals.
! Input syntax:
!   python lat_global <ix_universe>

case ('lat_global')
  read (line, *,  iostat = ios)  ix_universe
  
  if (ios /= 0) then
    call out_io (s_error$, r_name, "Correct input: lat_global <ix_universe>")
    return
  endif
  
  lat => s%u(ix_universe)%model%lat

  if (nl == size(ss%lines)) call re_allocate (ss%lines, nl+200, .false.)
  nl=nl+1; write(ss%lines(nl), '(5a)') &
    trim(lat%use_name), ';', trim(lat%lattice), ';', trim(lat%input_file_name) 
  nl=nl+1; write(ss%lines(nl), '(i0, a, i0)') &
    lat%n_ele_track, ';', lat%n_ele_max
  nl=nl+1; write(ss%lines(nl), '(3a)') &
    trim(species_name(lat%param%particle)), ';', trim(geometry_name(lat%param%geometry))




!----------------------------------------------------------------------
! Info on a given plot.
! Input syntax:
!   python plot1 <plot_name>
! Output syntax is serial output form. See documentation at beginning of this file.

case ('plot1')

  call tao_find_plots (err, line, 'BOTH', plot, print_flag = .false.)
  if (err) return

  if (allocated(plot)) then
    p => plot(1)%p
    nl=nl+1; write (ss%lines(nl), amt) 'x_axis_type;STRING;F;',          trim(p%x_axis_type), ';'
    nl=nl+1; write (ss%lines(nl), amt) 'x%label;STRING;F;',              trim(p%x%label), ';'
    nl=nl+1; write (ss%lines(nl), rmt) 'x%max;REAL;F;',                  p%x%max, ';'
    nl=nl+1; write (ss%lines(nl), rmt) 'x%min;REAL;F;',                  p%x%min, ';'
    nl=nl+1; write (ss%lines(nl), imt) 'x%major_div;INTEGER;F;',         p%x%major_div, ';'
    nl=nl+1; write (ss%lines(nl), imt) 'x%major_div_nominal;INTEGER;F;', p%x%major_div_nominal, ';'
    nl=nl+1; write (ss%lines(nl), imt) 'x%places;INTEGER;F;',            p%x%places, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'x%draw_label;LOGICAL;F;',        p%x%draw_label, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'x%draw_numbers;LOGICAL;F;',      p%x%draw_numbers, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'autoscale_x;LOGICAL;F;',         p%autoscale_x, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'autoscale_y;LOGICAL;F;',         p%autoscale_y, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'autoscale_gang_x;LOGICAL;F;',    p%autoscale_gang_x, ';'
    nl=nl+1; write (ss%lines(nl), lmt) 'autoscale_gang_y;LOGICAL;F;',    p%autoscale_gang_y, ';'
    do i = 1, size(p%graph)
      nl=nl+1; write (ss%lines(nl), amt) 'graph;STRING;F;',  trim(p%graph(i)%name), ';'
    enddo
  endif

!----------------------------------------------------------------------
! List of visible plot names.
! Input syntax: 
!   python plot_visible
! Output syntax:
!   <plot_name>;

case ('plot_visible')

  do i = 1, size(s%plot_page%region)
    region => s%plot_page%region(i)
    if (region%name == '') cycle
    if (.not. region%visible) cycle
    nl=nl+1; write (ss%lines(nl), '(2a)') trim(region%plot%name), ';'
  enddo

!----------------------------------------------------------------------
! List of plot templates.
! Input syntax:  
!   python plot_template
! Output syntax:

case ('plot_template')
  do i = 1, size(s%plot_page%template)
    p => s%plot_page%template(i)
    if (p%name == '') cycle
    if (p%name == 'scratch') cycle
    if (allocated(p%graph)) then
        nl=nl+1; write (ss%lines(nl), '(100a)') &
                          trim(p%name), (';.', trim(p%graph(j)%name), j = 1, size(p%graph)), ';'
    else
      nl=nl+1; write (ss%lines(nl), '(3x, a)') p%name 
    endif
  enddo

!----------------------------------------------------------------------
! List of all variable v1 arrays
! Input syntax: 
!   python var_all
! Output syntax:
!   <v1_var name>;<v1_var%v lower bound>;<v1_var%v upper bound>;<vars used in optimization list>;

case ('var_all')

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == '') cycle
    if (nl == size(ss%lines)) call re_allocate (ss%lines, nl+200, .false.)
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
    nl=nl+1; write(ss%lines(nl), '(2a, 2(i0, a),  a)') &
                    trim(v1_ptr%name), ';', lbound(v1_ptr%v, 1), ';', &
                    ubound(v1_ptr%v, 1), ';', trim(line)
    
  enddo
  
!----------------------------------------------------------------------
! List of variables in a given variable v1 array
! Input syntax: 
!   python var_v1 <v1_var>
! Output syntax:
!   <index>;<lat_ele_name>;<attribute_name>;<meas_value>;<model_value>;<design_value>;<good_user>;<useit_opt>;

case ('var_v1')

  call tao_find_var (err, line, v1_array = v1_array)

  if (.not. allocated(v1_array) .or. size(v1_array) /= 1) then
    nl=nl+1; ss%lines(nl) = '0;INVALID;;0;0;0;F;F'
  else
    v1_ptr => v1_array(1)%v1
    do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
      v_ptr => v1_ptr%v(i)
      if (.not. v_ptr%exists) cycle
      if (nl == size(ss%lines)) call re_allocate (ss%lines, nl+200, .false.)
!      nl=nl+1; write(ss%lines(nl), '(i0, 5a, 3(es15.8, a), 2(l1, a))') &
!                       v_ptr%ix_v1, ';', trim(v_ptr%ele_name), ';', trim(v_ptr%attrib_name), ';', &
!                       v_ptr%meas_value, ';', v_ptr%model_value, ';', &
!                       v_ptr%design_value, ';', v_ptr%good_user, ';', v_ptr%useit_opt, ';'
      nl=nl+1; write(ss%lines(nl), '(i0, 5a, 3(es15.8, a), 4a)') &
                       v_ptr%ix_v1, ';', trim(v_ptr%ele_name), ';', trim(v_ptr%attrib_name), ';', &
                       v_ptr%meas_value, ';', v_ptr%model_value, ';', &
                       v_ptr%design_value, ';', py_bool(v_ptr%good_user), ';', py_bool(v_ptr%useit_opt)
    enddo
  endif

!----------------------------------------------------------------------
! Info on an individual variable
! Input syntax: 
!   python var1 <var>
! Output syntax is serial output form. See documentation at beginning of this file.

case ('var1')

  call tao_find_var (err, line, v_array = v_array)

  if (.not. allocated(v_array) .or. size(v_array) /= 1) then
    nl=nl+1; ss%lines(nl) = '0;INVALID;;0;0;0;F;F;'
  else
    v_ptr => v_array(1)%v
    nl=nl+1; write(ss%lines(nl), amt)  'ele_name;STRING;F;',           trim(v_ptr%ele_name), ';'
    nl=nl+1; write(ss%lines(nl), amt)  'attrib_name;STRING;F;',        trim(v_ptr%attrib_name), ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_v1;INTEGER;F;',             v_ptr%ix_v1, ';'
    nl=nl+1; write(ss%lines(nl), imt)  'ix_var;INTEGER;F;',            v_ptr%ix_var, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'model_value;REAL;F;',          v_ptr%model_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'base_value;REAL;F;',           v_ptr%base_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'design_value;REAL;F;',         v_ptr%design_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'meas_value;REAL;T;',           v_ptr%meas_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'ref_value;REAL;T;',            v_ptr%ref_value, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'high_lim;REAL;F;',             v_ptr%high_lim, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'low_lim;REAL;F;',              v_ptr%low_lim, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'step;REAL;F;',                 v_ptr%step, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'weight;REAL;F;',               v_ptr%weight, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'delta_merit;REAL;F;',          v_ptr%delta_merit, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  'merit;REAL;F;',                v_ptr%merit, ';'
    nl=nl+1; write(ss%lines(nl), rmt)  's;REAL;F;',                    v_ptr%s, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'exists;LOGICAL;F;',            v_ptr%exists, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_var;LOGICAL;F;',          v_ptr%good_var, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_user;LOGICAL;T;',         v_ptr%good_user, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'good_opt;LOGICAL;F;',          v_ptr%good_opt, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'goov_plot;LOGICAL;F;',         v_ptr%good_plot, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'useit_plot;LOGICAL;F;',        v_ptr%useit_plot, ';'
    nl=nl+1; write(ss%lines(nl), lmt)  'useit_opt;LOGICAL;F;',         v_ptr%useit_opt, ';'

  endif

!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "***INTERNAL ERROR, SHOULDN'T BE HERE!")

end select

! return through scratch

scratch%n_lines = nl

if (doprint) call out_io (s_blank$, r_name, ss%lines(1:nl))

if (opened) then
  do i = 1, nl
    write (iu_write, '(a)') trim(ss%lines(i))
  enddo
  close (iu_write)
endif

!----------------------------------------------------------------------
! Helper function to write 'True' or 'False' strings from a logical

contains

function py_bool(bool) result(boolstring)
logical :: bool
character(5) :: boolstring 
if (bool) then
  boolstring = trim('True')
else
  boolstring = trim('False')
endif
end function

!----------------------------------------------------------------------
! contains

function py_string(chars) result(pystring)
character(*)::  chars
character(len_trim(chars)+2) :: pystring
pystring="'"//trim(chars)//"'"
end function

end subroutine




