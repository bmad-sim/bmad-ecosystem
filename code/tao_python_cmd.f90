!+
! Subroutine tao_python_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_python_cmd (input_str)

use tao_mod
use location_encode_mod

implicit none

character(n_char_show), allocatable, save :: lines(:)
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
type (tao_ele_shape_struct), pointer :: shape
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (rf_wake_lr_struct), pointer :: lr
type (ele_struct), pointer :: ele
type (coord_struct), target :: orb
type (ele_struct), target :: ele3
type (bunch_params_struct) bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch
type (random_state_struct) ran_state


character(*) input_str
character(24) imt, lmt, amt, iamt, ramt, f3mt, rmt, irmt, iimt
character(40) max_loc, loc_ele
character(200) line
character(20) cmd, command
character(20) :: r_name = 'tao_python_cmd'
character(20) :: cmd_names(32)= &
          ['plot_visible   ', 'plot_template  ', 'graph          ', 'curve          ', &
           'plot1          ', 'var_all        ', 'var_v1         ', 'var1           ', &
           'graph1         ', 'curve1         ', 'curve_sym      ', 'curve_line     ', &
           'data_all       ', 'data_d2        ', 'data_d1        ', 'data1          ', &
           'ele_all        ', 'ele1_all       ', 'beam_all       ', 'ele1_attrib    ', &
           'constraint_data', 'constraint_var ', 'global         ', 'lat_ele_list   ', &
           'lat_global     ', '               ', '               ', '               ', &
           '               ', '               ', '               ', '               ' ]

real(rp) target_value

integer i, j, ie, iu, ix, md, nl, ct
integer ix_ele, ix_ele1, ix_ele2, ix_branch, ix_universe

logical err, print_flag

!

call string_trim(input_str, line, ix)
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

nl = 0
if (.not. allocated(lines)) allocate (lines(200))

select case (command)

!----------------------------------------------------------------------
! datums used in optimizations (%useit_opt = T). 
! Syntax:  constrint_data

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

    if (nl == size(lines)) call re_allocate(lines, 2*nl)
    nl=nl+1; write (lines(nl), '(10a, 3(es12.4, a), 2a)') &
        trim(tao_datum_name(d_ptr)), ';', trim(tao_datum_type_name(d_ptr)), ';', &
        trim(d_ptr%ele_ref_name), ';', trim(d_ptr%ele_start_name), ';', trim(d_ptr%ele_name), ';', &
        d_ptr%meas_value, ';', d_ptr%model_value, ';', d_ptr%merit, ';', trim(max_loc), ';'

  enddo
enddo

!----------------------------------------------------------------------
! variables used in optimizations (%useit_opt = T). 
! Syntax:  constrint_data

case ('constraint_vars')

do i = 1, s%n_var_used
  v_ptr => s%var(i)
  if (.not. v_ptr%useit_opt) cycle

  u => s%u(v_ptr%this(1)%ix_uni)
  branch => u%model%lat%branch(v_ptr%this(1)%ix_branch)
  ct = branch%ele(v_ptr%this(1)%ix_ele)%lord_status

  loc_ele = ''
  if (ct /= group_lord$ .and. ct /= overlay_lord$ .and. ct /= multipass_lord$) then
    write (loc_ele, '(f8.2)') branch%ele(v_ptr%this(1)%ix_ele)%s
    call string_trim (loc_ele, loc_ele, ix)
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

  if (nl == size(lines)) call re_allocate(lines, 2*nl)
  nl=nl+1; write (lines(nl), '(6a, 3(es12.4, a))') &
        trim(tao_var1_name(v_ptr)), ';', trim(tao_var_attrib_name(v_ptr)), ';', &
        trim(loc_ele), ';', target_value, ';', v_ptr%model_value, ';', v_ptr%merit, ';'

enddo

!----------------------------------------------------------------------
! D1 level data. 
! Syntax: 
!   <index>;<data_type>;<merit_type>;<ele_ref_name>;<ele_start_name>;<ele_name>;<meas_value>;<model_value>;<design_value>;<good_user>;<useit_opt>;<useit_plot>;

case ('data_d1')

  call tao_find_data (err, line, d1_array = d1_array)

  if (.not. allocated(d1_array) .or. size(d1_array) /= 1) then
    nl=nl+1; lines(nl) = '0;INVALID;;;;0;0;0;F;F;'
  else
    d1_ptr => d1_array(1)%d1
    do i = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
      d_ptr => d1_ptr%d(i)
      if (.not. d_ptr%exists) cycle
      if (nl == size(lines)) call re_allocate (lines, nl+200, .false.)
      nl=nl+1; write(lines(nl), '(i0, 11a, 3(es16.8, a), 3(l1, a))') &
                       d_ptr%ix_d1, ';', trim(d_ptr%data_type), ';', trim(d_ptr%merit_type), ';', &
                       trim(d_ptr%ele_ref_name), ';', trim(d_ptr%ele_start_name), ';', &
                       trim(d_ptr%ele_name), ';', d_ptr%meas_value, ';', d_ptr%model_value, ';', &
                       d_ptr%design_value, ';', d_ptr%good_user, ';', d_ptr%useit_opt, ';', d_ptr%useit_plot, ';'
    enddo
  endif

!----------------------------------------------------------------------
! D2 level data. 
! Syntax: 
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

        nl=nl+1; write (lines(nl), '(i0, 5a, i0, a, i0, 3a)') iu, ';', trim(d2_ptr%name), ';', &
              trim(d1_ptr%name), ';', lbound(d1_ptr%d, 1), ';', ubound(d1_ptr%d, 1), ';', trim(line), ';'

      enddo

    enddo

  enddo

!----------------------------------------------------------------------
! Individual datum
! Syntax:
!   <component_name>;<type>;<variable>;<component_value>;
! <type> is one of:
!   STRING
!   INTEGER
!   REAL
!   LOGICAL
! <variable> indicates if the component can be varied. It is one of:
!   T
!   F

case ('data1')

  call tao_find_data (err, line, d_array = d_array)

  if (.not. allocated(d_array) .or. size(d_array) /= 1) then
    nl=nl+1; lines(nl) = 'INVALID;STRING;F;0;'
  else
    d_ptr => d_array(1)%d
    amt = '(3a)'
    imt = '(a,i0,a)'
    rmt = '(a,es16.8,a)'
    lmt = '(a,l1,a)'
    nl=nl+1; write(lines(nl), amt)  'ele_name;STRING;F;',           trim(d_ptr%ele_name), ';'
    nl=nl+1; write(lines(nl), amt)  'ele_start_name;STRING;F;',     trim(d_ptr%ele_start_name), ';'
    nl=nl+1; write(lines(nl), amt)  'ele_ref_name;STRING;F;',       trim(d_ptr%ele_ref_name), ';'
    nl=nl+1; write(lines(nl), amt)  'data_type;STRING;F;',          trim(d_ptr%data_type), ';'
    nl=nl+1; write(lines(nl), amt)  'data_source;STRING;F;',        trim(d_ptr%data_source), ';'
    nl=nl+1; write(lines(nl), imt)  'ix_branch;INTEGER;F;',         d_ptr%ix_branch, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_ele;INTEGER;F;',            d_ptr%ix_ele, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_ele_start;INTEGER;F;',      d_ptr%ix_ele_start, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_ele_ref;INTEGER;F;',        d_ptr%ix_ele_ref, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_ele_merit;INTEGER;F;',      d_ptr%ix_ele_merit, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_dmodel;INTEGER;F;',         d_ptr%ix_dModel, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_d1;INTEGER;F;',             d_ptr%ix_d1, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_data;INTEGER;F;',           d_ptr%ix_data, ';'
    nl=nl+1; write(lines(nl), imt)  'ix_bunch;INTEGER;F;',          d_ptr%ix_bunch, ';'
    nl=nl+1; write(lines(nl), rmt)  'model;REAL;F;',                d_ptr%model_value, ';'
    nl=nl+1; write(lines(nl), rmt)  'design;REAL;F;',               d_ptr%design_value, ';'
    nl=nl+1; write(lines(nl), rmt)  'meas;REAL;T;',                 d_ptr%meas_value, ';'
    nl=nl+1; write(lines(nl), rmt)  'ref;REAL;T;',                  d_ptr%ref_value, ';'
    nl=nl+1; write(lines(nl), rmt)  'base;REAL;F;',                 d_ptr%base_value, ';'
    nl=nl+1; write(lines(nl), rmt)  'old;REAL;F;',                  d_ptr%old_value   , ';'
    nl=nl+1; write(lines(nl), rmt)  'invalid;REAL;F;',              d_ptr%invalid_value, ';'
    nl=nl+1; write(lines(nl), rmt)  's;REAL;F;',                    d_ptr%s, ';'
    nl=nl+1; write(lines(nl), amt)  'merit_type;STRING;F;',         trim(d_ptr%merit_type), ';'
    nl=nl+1; write(lines(nl), rmt)  'merit;REAL;F;',                d_ptr%merit, ';'
    nl=nl+1; write(lines(nl), rmt)  'delta_merit;REAL;F;',          d_ptr%delta_merit, ';'
    nl=nl+1; write(lines(nl), rmt)  'weight;REAL;F;',               d_ptr%weight, ';'
    nl=nl+1; write(lines(nl), lmt)  'exists;LOGICAL;F;',            d_ptr%exists, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_model;LOGICAL;F;',        d_ptr%good_model, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_design;LOGICAL;F;',       d_ptr%good_design, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_base;LOGICAL;F;',         d_ptr%good_base , ';'
    nl=nl+1; write(lines(nl), lmt)  'good_meas;LOGICAL;F;',         d_ptr%good_meas, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_ref;LOGICAL;F;',          d_ptr%good_ref, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_user;LOGICAL;T;',         d_ptr%good_user, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_opt;LOGICAL;F;',          d_ptr%good_opt, ';'
    nl=nl+1; write(lines(nl), lmt)  'good_plot;LOGICAL;F;',         d_ptr%good_plot, ';'
    nl=nl+1; write(lines(nl), lmt)  'useit_plot;LOGICAL;F;',        d_ptr%useit_plot, ';'
    nl=nl+1; write(lines(nl), lmt)  'useit_opt;LOGICAL;F;',         d_ptr%useit_opt, ';'
  endif


!----------------------------------------------------------------------
! All parameters associated with given element. 
! Synrax: ele1_all <ix_ele> <ix_branch> <ix_universe>

case ('ele1_all')


!----------------------------------------------------------------------
! Attribute list. Shows which attributes are adjustable.
! Synrax: ele1_all <ix_ele> <ix_branch> <ix_universe>

case ('ele1_attrib')


!----------------------------------------------------------------------
! global parameters
! Syntax: global

case ('global')

!----------------------------------------------------------------------
! Lattice element list.
! lat_ele <ix_ele_start> <ix_ele_end> <ix_branch>

case ('lat_ele_list')



!----------------------------------------------------------------------
! Lattice globals.
! Syntax: lat_global <ix_universe>

case ('lat_global')

!----------------------------------------------------------------------
! Output info on a given plot.
! Syntax: plot1 <plot_name>

case ('plot1')

  call tao_find_plots (err, line, 'BOTH', plot, print_flag = .false.)
  if (err) return

  if (allocated(plot)) then
    nl=nl+1; write (lines(nl), amt) 'x_axis_type          = ', p%x_axis_type
    nl=nl+1; write (lines(nl), amt) 'x%label              = ', p%x%label
    nl=nl+1; write (lines(nl), rmt) 'x%max                = ', p%x%max
    nl=nl+1; write (lines(nl), rmt) 'x%min                = ', p%x%min
    nl=nl+1; write (lines(nl), imt) 'x%major_div          = ', p%x%major_div
    nl=nl+1; write (lines(nl), imt) 'x%major_div_nominal  = ', p%x%major_div_nominal
    nl=nl+1; write (lines(nl), imt) 'x%places             = ', p%x%places
    nl=nl+1; write (lines(nl), lmt) 'x%draw_label         = ', p%x%draw_label
    nl=nl+1; write (lines(nl), lmt) 'x%draw_numbers       = ', p%x%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'autoscale_x          = ', p%autoscale_x
    nl=nl+1; write (lines(nl), lmt) 'autoscale_y          = ', p%autoscale_y
    nl=nl+1; write (lines(nl), lmt) 'autoscale_gang_x     = ', p%autoscale_gang_x
    nl=nl+1; write (lines(nl), lmt) 'autoscale_gang_y     = ', p%autoscale_gang_y
    do i = 1, size(p%graph)
      nl=nl+1; write (lines(nl), amt) 'graph = ', p%graph(i)%name
    enddo
  endif

!----------------------------------------------------------------------
! output list of visible plot names.
! Syntax: plot_visible

case ('plot_visible')

  do i = 1, size(s%plot_region)
    region => s%plot_region(i)
    if (region%name == '') cycle
    if (.not. region%visible) cycle
    nl=nl+1; write (lines(nl), '(a)') region%plot%name
  enddo

!----------------------------------------------------------------------
! output list of plot templates.
! Syntax:  plot_template

case ('plot_template')
  do i = 1, size(s%template_plot)
    p => s%template_plot(i)
    if (p%name == '') cycle
    if (p%name == 'scratch') cycle
    if (allocated(p%graph)) then
        nl=nl+1; write (lines(nl), '(100a)') &
                          trim(p%name), (';.', trim(p%graph(j)%name), j = 1, size(p%graph))
    else
      nl=nl+1; write (lines(nl), '(3x, a)') p%name 
    endif
  enddo

!----------------------------------------------------------------------
! Output list of variable v1 arrays
! Syntax: var_all

case ('var_all')

    do i = 1, s%n_v1_var_used
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == '') cycle
      if (nl == size(lines)) call re_allocate (lines, nl+200, .false.)
      call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
      nl=nl+1; write(lines(nl), '(i5, 2x, 2a, i0, a, i0, a, t50, a)') v1_ptr%ix_v1, &
                      trim(v1_ptr%name), '[', lbound(v1_ptr%v, 1), ':', &
                      ubound(v1_ptr%v, 1), ']', trim(line)
      
    enddo
  
!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "***INTERNAL ERROR, SHOULDN'T BE HERE!")

end select

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine
