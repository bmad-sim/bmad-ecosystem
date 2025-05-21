module tao_graph_setup_mod

use tao_interface
use tao_command_mod
use expression_mod, only: expression_stack_value

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_setup (plot, graph)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve

integer i, iu
logical found

character(*), parameter :: r_name = 'tao_graph_setup'

!

if (.not. graph%is_valid) return
graph%text_legend_out = graph%text_legend

if (graph%type == 'floor_plan') return  ! Nothing to do.
if (graph%type == 'key_table') return   ! Nothing to do

if (graph%type == 'lat_layout') then
  u => tao_pointer_to_universe(graph%ix_universe, .true.)
  if (.not. associated(u)) then
    graph%is_valid = .false.
    write (graph%why_invalid, '(a, i0, a)') 'BAD UNIVERSE INDEX', graph%ix_universe
  elseif (.not. u%is_on) then
    graph%is_valid = .false.
    write (graph%why_invalid, '(a, i0, a)') 'UNIVERSE ', u%ix_uni, ' IS OFF!'
  endif
  return
endif

if (.not. allocated (graph%curve)) then
  call out_io (s_warn$, r_name, 'NO CURVES ASSOCIATED WITH: ' // tao_graph_name(graph))
  graph%is_valid = .false.
  graph%why_invalid = 'NO ASSOCIATED CURVES'
  return
endif

!

do i = 1, size(graph%curve)
  graph%curve(i)%valid = .true.  ! Assume no problem
  graph%curve(i)%message_text = ''
  call tao_remove_blank_characters(graph%curve(i)%component)
enddo

if (associated(tao_hook_graph_setup_ptr)) then
  call tao_hook_graph_setup_ptr (plot, graph, found)
  if (found) return
endif

!

select case (graph%type)
case ('phase_space')
  call tao_graph_phase_space_setup (plot, graph)

case ('data', 'lat_layout')
  if (plot%x_axis_type == 'data') then
    call tao_graph_data_slice_setup(plot, graph)
  else
    call tao_graph_data_setup(plot, graph)
  endif

case ('histogram')
  call tao_graph_histogram_setup (plot, graph)

case ('dynamic_aperture')
  call tao_graph_dynamic_aperture_setup (plot, graph)

case ('control_curve')  ! Not yet implemented
  call tao_graph_controller_setup (graph)

end select

! Renormalize

if (allocated (graph%curve)) then
  do i = 1, size(graph%curve)
    curve => graph%curve(i)
    if (.not. curve%valid) cycle
    ! Unlike other curves, multi_turn_orbit curves gets calculated in tao_lattice_calc which takes care of the scaling.
    if (curve%data_source == 'multi_turn_orbit' .or. curve%data_source == 'rel_multi_turn_orbit') cycle 

    if (allocated(curve%x_symb)) then
        curve%x_symb = curve%x_symb * graph%x_axis_scale_factor
        curve%y_symb = curve%y_symb * curve%y_axis_scale_factor
    endif

    if (allocated(curve%err_symb)) curve%err_symb = curve%err_symb * curve%y_axis_scale_factor

    if (allocated(curve%x_line)) then
      curve%x_line = curve%x_line * graph%x_axis_scale_factor
      curve%y_line = curve%y_line * curve%y_axis_scale_factor
    endif

    if (graph%type == 'histogram') then
      curve%hist%minimum = curve%hist%minimum * graph%x_axis_scale_factor
      curve%hist%maximum = curve%hist%maximum * graph%x_axis_scale_factor
      curve%hist%width   = curve%hist%width * graph%x_axis_scale_factor
    endif
  enddo
endif

if (associated(tao_hook_graph_postsetup_ptr)) call tao_hook_graph_postsetup_ptr (plot, graph)

end subroutine tao_graph_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_controller_setup (graph)

implicit none

type (tao_graph_struct), target :: graph
type (tao_plot_struct), pointer :: plot
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele, slave
type (control_struct), pointer :: ctl
type (expression_atom_struct), pointer :: stack(:) 

real(rp) var, var0, value
real(rp), allocatable :: x(:), y(:)
real(rp), pointer :: y_knot(:)

integer i, j, ix, ix_slave, n_curve_pts, n_loc

logical err, ok, allocated_stack

character(100) err_str
character(40) name
character(*), parameter :: r_name = 'tao_graph_controller_setup'

!

plot => graph%p

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

allocate (x(0:n_curve_pts), y(0:n_curve_pts))

do i = 1, size(graph%curve)
  curve => graph%curve(i)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))

  ix = index(curve%data_type, ':')
  if (ix == 0) then
    ix_slave = 1
    name = curve%data_type
  else
    name = curve%data_type(:ix-1)
    if (.not. is_integer(curve%data_type(ix+1:), ix_slave)) then
      curve%valid = .false.
      call out_io (s_warn$, r_name, 'CURVE DATA_TYPE HAS NON-INTEGER SLAVE INDEX: ' // curve%data_type, &
                                     'FOR CURVE: ' // tao_curve_name(curve))
      cycle
    endif
  endif

  call lat_ele_locator (name, u%model%lat, eles, n_loc, err)
  if (n_loc == 0) then
    call out_io (s_warn$, r_name, 'CANNOT FIND CONTROLLER ELEMENT: ' // name, &
                                   'FOR CURVE: ' // tao_curve_name(curve))
    curve%valid = .false.
    cycle
  endif

  if (n_loc > 0) then
    call out_io (s_warn$, r_name, 'MULTIPLE ELEMENTS FOUND THAT MATCH: ' // name, &
                                  'USING THE FIRST ONE.', &
                                  'FOR CURVE: ' // tao_curve_name(curve))
  endif

  ele => eles(1)%ele
  if (ele%key == group$ .or. ele%key == overlay$) then
    if (ix_slave < 1 .or. ix_slave > ele%n_slave) then
      call out_io (s_warn$, r_name, 'SLAVE INDEX OF CONTROLLER ELEMENT OUT OF RANGE: ' // curve%data_type, &
                                     'FOR CURVE: ' // tao_curve_name(curve))
      curve%valid = .false.
      cycle
    endif
    slave => pointer_to_slave (ele, ix_slave, ctl)
    stack => ctl%stack
    allocated_stack = allocated(ctl%stack)
    y_knot => ctl%y_knot

  elseif (ele%key == ramper$) then
    if (ix_slave < 1 .or. ix_slave > size(ele%control%ramp)) then
      call out_io (s_warn$, r_name, 'SLAVE INDEX OF CONTROLLER ELEMENT OUT OF RANGE: ' // curve%data_type, &
                                     'FOR CURVE: ' // tao_curve_name(curve))
      curve%valid = .false.
      cycle
    endif
    stack => ele%control%ramp(ix_slave)%stack
    allocated_stack = allocated(ele%control%ramp(ix_slave)%stack)
    y_knot => ele%control%ramp(ix_slave)%y_knot

  else
    call out_io (s_warn$, r_name, 'ELEMENT IS NOT A GROUP, RAMPER, OR OVERLAY: ' // name, &
                                   'FOR CURVE: ' // tao_curve_name(curve))
    curve%valid = .false.
    cycle
  endif

  var0 = ele%control%var(1)%value

  do j = 1, n_curve_pts
    var = graph%x%eval_min + (j - 1) * (graph%x%eval_max - graph%x%eval_min) / n_curve_pts
    ele%control%var(1)%value = var
    if (allocated_stack) then
      value = expression_stack_value(stack, err, err_str, ele%control%var, .false.)
    else
      call spline_akima_interpolate (ele%control%x_knot, y_knot, value, ok, value)
    endif

    x(j) = var
    y(j) = value
  enddo

  call re_allocate (curve%x_symb, n_curve_pts)
  call re_allocate (curve%y_symb, n_curve_pts)
  curve%x_symb = x
  curve%y_symb = y

  call re_allocate (curve%x_line, n_curve_pts)
  call re_allocate (curve%y_line, n_curve_pts)
  curve%x_line = x
  curve%y_line = y
enddo

end subroutine tao_graph_controller_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_data_slice_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

real(rp), pointer :: symb(:)
real(rp) value

integer i, j, k, m, n_symb, ix

character(200) name
character(40) :: r_name = 'tao_graph_data_slice_setup'

logical err

!

if (size(graph%curve) == 0) then
  graph%is_valid = .false.
  graph%why_invalid = 'No curves asscoiated with graph.'
  return
endif

!

graph%title_suffix = ''

do k = 1, size(graph%curve)
  curve => graph%curve(k)

  if (index(curve%data_type, '#ref') /= 0) graph%title_suffix = &
                  trim(graph%title_suffix) // '[At: ' // trim(curve%ele_ref_name) // ']'
enddo

!

if (all(graph%curve%component == graph%curve(1)%component)) graph%title_suffix = trim(graph%title_suffix) // &
                                                                     ' [' // trim(graph%curve(1)%component) // ']'

! loop over all curves

curve_loop: do k = 1, size(graph%curve)

  curve => graph%curve(k)

  ! Find data points

  do i = 1, 2  ! x-axis, y-axis
    if (i == 1) name = curve%data_type_x
    if (i == 2) name = curve%data_type
    call tao_data_type_substitute (name, name, curve, graph)
    if (i == 1) call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, &
                                                 scratch%x, err, info = scratch%info_x, dflt_component = curve%component)
    if (i == 2) call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, &
                                                 scratch%y, err, info = scratch%info_y, dflt_component = curve%component)
    if (err) then
      call tao_set_curve_invalid (curve, 'CANNOT FIND DATA.')
      cycle curve_loop
    endif
  enddo

  ! How many good points?

  if (size(scratch%x) /= size(scratch%y)) then
    call tao_set_curve_invalid (curve, 'ARRAY SIZES ARE NOT THE SAME FOR BOTH AXES.')
    cycle
  endif

  n_symb = count(scratch%info_x%good .and. scratch%info_y%good)

  call re_allocate (curve%x_symb, n_symb)
  call re_allocate (curve%y_symb, n_symb)
  call re_allocate (curve%ix_symb, n_symb)

  ! Transfer the values

  curve%x_symb = pack (scratch%x, mask = scratch%info_x%good .and. scratch%info_y%good)
  curve%y_symb = pack (scratch%y, mask = scratch%info_x%good .and. scratch%info_y%good)

  ! Calc symbol index

  if (curve%data_index == '') then
    curve%ix_symb = [(i, i = 1, n_symb)]
  else
    call tao_data_type_substitute (curve%data_index, name, curve, graph)
    call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, &
                                          scratch%x, err, info = scratch%info_ix, dflt_component = curve%component)
    if (size(scratch%info_x) == size(scratch%info_y)) then
      curve%ix_symb = pack (nint(scratch%x), mask = scratch%info_x%good .and. scratch%info_y%good)
    else
      call out_io (s_warn$, r_name, &
          'SIZE OF SYMBOL INDEX ARRAY IS WRONG IN CURVE: ' // tao_curve_name(curve), &
          'CURVE%DATA_INDEX: ' // curve%data_index)
    endif
  endif

  ! The line data just goes through the symbols

  call re_allocate (curve%x_line, n_symb)
  call re_allocate (curve%y_line, n_symb)
  curve%x_line = curve%x_symb
  curve%y_line = curve%y_symb

enddo curve_loop

end subroutine tao_graph_data_slice_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_type_substitute (template, str_out, curve, graph)
!
! Routine substitute the appropriate data type string for instances of "#ref" and 
! "#comp" in template. 
!
! Additionally, if template does not have a "|" character,
! the string "|" + component will be added at the end of str_out.
!
! Input:
!   template    -- character(*): String template.
!   curve       -- tao_curve_struct: curve%ele_ref_name is substituted for all instances of "#ref".
!   graph       -- tao_graph_struct: 
!
! Output:
!   str_out     -- character(*): String with substitutions.
!-

subroutine tao_data_type_substitute (template, str_out, curve, graph)

implicit none

type (tao_curve_struct) curve
type (tao_graph_struct) graph
character(*) template, str_out
integer ix

!

call string_trim(template, str_out, ix)
if (str_out(1:11) == 'expression:') str_out = str_out(12:)

!

do
  ix = index(str_out, '#ref')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(curve%ele_ref_name) // trim(str_out(ix+4:))
enddo

do
  ix = index(str_out, '#comp')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(curve%component) // trim(str_out(ix+5:))
enddo

if (index(str_out, '|') == 0) str_out = trim(str_out) // '|' // curve%component

end subroutine tao_data_type_substitute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_phase_space_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (beam_struct), pointer :: beam
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_x, d1_y
type (coord_struct), pointer :: p(:)
type (tao_beam_branch_struct), pointer :: bb

real(rp) v_mat(4,4), v_inv_mat(4,4), g_mat(4,4), g_inv_mat(4,4)
real(rp) mat4(4,4), sigma_mat(4,4), theta, theta_xy, rx, ry, phi
real(rp) emit_a, emit_b

integer k, n, m, ib, ix1_ax, ix2_ax, ix3_ax, ix, i
integer n_curve_pts

logical err, same_uni

character(40) name
character(*), parameter :: r_name = 'tao_graph_phase_space_setup'

! Set up the graph suffix

if (size(graph%curve) == 0) then
  graph%is_valid = .false.
  graph%why_invalid = 'NO CURVES ASSCOIATED WITH GRAPH.'
  return
endif

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

same_uni = .true.
ix = tao_universe_index(tao_curve_ix_uni(graph%curve(1)))
do k = 2, size(graph%curve)
  curve => graph%curve(k)
  if (tao_universe_index(tao_curve_ix_uni(curve)) /= ix) same_uni = .false.
enddo

graph%title_suffix = ''
do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
  if (.not. tao_curve_check_universe(curve, u)) cycle
  ele => tao_curve_ele_ref(curve, .true.)
  if (.not. associated(ele)) then
    call tao_set_curve_invalid (curve, 'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name)
    cycle
  endif

  name = curve%ele_ref_name
  if (name == '') name = ele%name

  if (len_trim(graph%title_suffix) + len_trim(name) + 10 > len(graph%title_suffix)) then
    graph%title_suffix = trim(graph%title_suffix) // ' etc...'  ! prevent overflow
    exit
  endif

  if (same_uni) then
    write (graph%title_suffix, '(2a, i0, 3a)') trim(graph%title_suffix), '[', ele%ix_ele, ': ', trim(name), ']'
  else
    write (graph%title_suffix, '(2a, i0, a, i0, 3a)') trim(graph%title_suffix), &
                                                          '[', u%ix_uni, '@', ele%ix_ele, ': ', trim(name), ']'
  endif
enddo

! loop over all curves

do k = 1, size(graph%curve)

  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
  ele => tao_curve_ele_ref(curve, .true.)
  if (.not. associated(ele)) then
    call out_io (s_error$, r_name, 'CURVE REFERENCE ELEMENT IS NOT AN ELEMENT THAT IS TRACKED THROUGH: ' // curve%ele_ref_name, &
                                   'FOR CURVE: ' // tao_curve_name(curve))
    cycle
  endif

  select case (curve%component)
  case('model', '')
    ele => u%model%lat%branch(ele%ix_branch)%ele(ele%ix_ele)
  case('design')
    ele => u%design%lat%branch(ele%ix_branch)%ele(ele%ix_ele)
  case('base')
    ele => u%base%lat%branch(ele%ix_branch)%ele(ele%ix_ele)
  case default
    call out_io (s_error$, r_name, 'IF THE CURVE COMPONENT IS SET, IT MUST BE SET TO ONE OF: "model", "design", OR "base".', &
                                   'FOR CURVE: ' // tao_curve_name(curve))
  end select

  ! fill the curve data arrays

  if (allocated (curve%ix_symb) .and. curve%data_source /= 'multi_turn_orbit' .and. &
            curve%data_source /= 'rel_multi_turn_orbit') deallocate (curve%ix_symb, curve%x_symb, curve%y_symb)
  if (allocated (curve%x_line))  deallocate (curve%x_line, curve%y_line)

  !----------------------------

  if (curve%data_source == 'beam') then
    if (.not. associated(ele)) then
      call out_io (s_warn$, r_name, 'REFERENCE ELEMENT DOES NOT EXIST: ' // trim(curve%ele_ref_name), &
                    'CANNOT DO PHASE_SPACE PLOTTING FOR CURVE: ' // tao_curve_name(curve))
      curve%g%why_invalid = 'NO BEAM AT ELEMENT'
      return
    endif

    beam => u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
    if (.not. allocated(beam%bunch)) then
      call out_io (s_warn$, r_name, 'NO BEAM AT ELEMENT: ' // ele_full_name(ele), &
                    'CANNOT DO PHASE_SPACE PLOTTING FOR CURVE: ' // tao_curve_name(curve))
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      curve%g%why_invalid = 'NO BEAM AT ELEMENT'
      return
    endif

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1, size(beam%bunch)
        n = n + count(beam%bunch(ib)%particle%state == alive$)
      enddo
    else
      n = count(beam%bunch(curve%ix_bunch)%particle%state == alive$)
    endif

    if (n == 0) then
      call out_io (s_warn$, r_name, 'NO LIVE BEAM PARTICLES PRESENT AT ELEMENT: ' // trim(ele%name), &
                    'CANNOT DO PHASE_SPACE PLOTTING FOR CURVE: ' // tao_curve_name(curve))
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      curve%g%why_invalid = 'NO LIVE BEAM PARTICLES PRESENT'
      return
    endif

    call re_allocate (curve%ix_symb, n)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)
    if (graph%symbol_size_scale > 0) call re_allocate (curve%symb_size, n)
    if (curve%z_color%is_on) call re_allocate (curve%z_symb, n)

    n = 0
    do ib = 1, size(beam%bunch)
      if (curve%ix_bunch /= 0 .and. curve%ix_bunch /= ib) cycle
      p => beam%bunch(ib)%particle
      m = count(beam%bunch(ib)%particle%state == alive$)
      call tao_particle_data_value (curve%data_type_x, p, scratch%axis1, err, ele, ib); if (err) return
      call tao_particle_data_value (curve%data_type,   p, scratch%axis2, err, ele, ib); if (err) return
      curve%x_symb(n+1:n+m) = pack(scratch%axis1, mask = (p%state == alive$))
      curve%y_symb(n+1:n+m) = pack(scratch%axis2, mask = (p%state == alive$))
      if (curve%z_color%is_on) then
        call tao_particle_data_value (curve%z_color%data_type, p, scratch%axis3, err, ele, ib); if (err) return
        curve%z_symb(n+1:n+m) = pack(scratch%axis3, mask = (p%state == alive$))
      endif
      if (graph%symbol_size_scale > 0) curve%symb_size(n+1:n+m) = pack(graph%symbol_size_scale * &
                           sqrt(p(:)%field(1)**2 + p(:)%field(2)**2), mask = (p%state == alive$))
      curve%ix_symb(n+1:n+m) = pack([(i, i = 1,m)], mask = (p%state == alive$))
      n = n + count(p%state == alive$)
    enddo

  !----------------------------

  elseif (curve%data_source == 'multi_turn_orbit' .or. curve%data_source == 'rel_multi_turn_orbit') then
    ! Everything is handled in tao_lattice_calc

  elseif (curve%data_source == 'twiss') then

    n = 2 * n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    call make_v_mats (ele, v_mat, v_inv_mat)
    call make_g_mats (ele, g_mat, g_inv_mat)

    bb => u%model_branch(0)%beam
    mat4 = matmul(v_mat, g_inv_mat)
    emit_a = bb%beam_init%a_emit
    if (emit_a == 0) emit_a = 1e-6  ! default value
    emit_b = bb%beam_init%b_emit
    if (emit_b == 0) emit_b = 1e-6  ! default value

    sigma_mat =  0
    sigma_mat(1,1) = emit_a
    sigma_mat(2,2) = emit_a
    sigma_mat(3,3) = emit_b
    sigma_mat(4,4) = emit_b
    sigma_mat = matmul (matmul (mat4, sigma_mat), transpose(mat4))

    ! find phase space axes to plot

    ix1_ax = tao_phase_space_axis_index (curve%data_type_x, err); if (err) return
    ix2_ax = tao_phase_space_axis_index (curve%data_type, err); if (err) return

    if (ix1_ax > 4 .or. ix2_ax > 4) then
      call out_io (s_warn$, r_name, &
        'Z OR PZ PHASE SPACE PLOTTING NOT YET IMPLEMENTED FOR "twiss" DATA_SOURCE.')
      return
    endif

    rx = sqrt(sigma_mat(ix1_ax, ix1_ax))
    ry = sqrt(sigma_mat(ix2_ax, ix2_ax))

    do n = size(graph%text_legend_out), 1, -1
      if (graph%text_legend_out(n) /= '') exit
    enddo

    write (graph%text_legend_out(n+1), '(a, es9.2)') 'emit_a:', emit_a
    write (graph%text_legend_out(n+2), '(a, es9.2)') 'emit_b:', emit_b

    if(rx == 0 .or. ry == 0) then
      theta_xy = 0
      write (graph%text_legend_out(n+3), '(a, f10.4)') 'Theta_tilt (rad):', 0
    else
      theta_xy =  asin(sigma_mat(ix1_ax, ix2_ax) / (rx * ry))
      phi = 0.5 *atan2((rx**2+ry**2) * sin(2*theta_xy), &
                              (rx**2-ry**2) * cos(2*theta_xy)) - theta_xy
      write (graph%text_legend_out(n+3), '(a, f10.4)') 'Theta_tilt (rad):', phi
    endif

    n = 2 * n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    do i = 1, n
      theta = (i-1) * twopi / (n-1)
      curve%x_line(i) = rx * cos(theta)
      curve%y_line(i) = ry * sin(theta + theta_xy)
    enddo

  else
    call out_io (s_warn$, r_name, &
        'INVALID CURVE%DATA_SOURCE: ' // curve%data_source, &
        'FOR CURVE: '// tao_curve_name(curve))
    return
  endif

enddo

end subroutine tao_graph_phase_space_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_dynamic_aperture_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_curve_struct), allocatable :: temp_curve(:)
type (tao_universe_struct), pointer :: u
type (tao_dynamic_aperture_struct), pointer :: da
type (ele_struct), pointer :: ref_ele
type (aperture_scan_struct), pointer :: scan
type (coord_struct), allocatable :: orbit(:)
integer :: i, j, k, ib, n_da, n_da_curve, n

logical err

character(40) name
character(*), parameter :: r_name = 'tao_graph_dynamic_aperture_setup'

! Only the graph's universe for now
u => tao_pointer_to_universe (graph%ix_universe)

da => u%dynamic_aperture
n_da = size(da%scan)

! Count DA curves 

n_da_curve = 0
do i = 1, size(graph%curve)
  curve => graph%curve(i)
  if (substr(curve%data_type,1,12) == 'beam_ellipse') cycle
  n_da_curve = n_da_curve + 1
  ! Assign a default type 
  if (curve%data_type == '') curve%data_type = 'dynamic_aperture'
  select case (curve%data_type)
  case ('dynamic_aperture_centered')
    call out_io (s_warn$, r_name, 'curve%data_type = "dynamic_aperture_centered" is now "dynamic_aperture".', &
                                  'Note: "dynamic_aperture_ref0" gives curves centered on the zero orbit.')
    curve%data_type = 'dynamic_aperture'
  case ('dynamic_aperture', 'dynamic_aperture_ref0')
  case default
    call tao_set_curve_invalid (curve, 'BAD DATA_TYPE FOR DYNAMIC_APERTURE PLOT CURVE: ' // quote(curve%data_type))
  end select
enddo 

graph%is_valid = (n_da_curve /= 0)
if (.not. graph%is_valid) then
  write (graph%why_invalid, '(a, i0)') 'NO DYNAMIC APERTURES CURVES DEFINED FOR UNIVERSE ', u%ix_uni
  return
endif

err = (.not. allocated(da%scan))
if (allocated(da%scan)) then
  if (.not. allocated(da%scan(1)%point)) err = .true.
endif

if (err) then
  graph%is_valid = .false.
  write (graph%why_invalid, '(a, i0)') 'DYNAMIC APERTURE NOT CALCULATED FOR UNIVERSE ', u%ix_uni
  return
endif

! Loop over curves.
! It is allowed to have more curves than aperture scans. In this case, ignore unused curves.

do i = 1, size(graph%curve)
  curve => graph%curve(i)
  if (substr(curve%data_type,1,12) == 'beam_ellipse') then
    call tao_curve_beam_ellipse_setup(curve)
    cycle
  endif

  if (.not. curve%valid) cycle

  if (.not. is_integer(curve%data_index, k)) then
    call tao_set_curve_invalid (curve, 'DATA_INDEX IS NOT AN INTEGER: ' // quote(curve%data_index))
    cycle
  endif

  if (k > size(da%scan)) then
    curve%valid = .false.
    curve%why_invalid = 'IGNORE'
    cycle
  endif

  scan => da%scan(k)
  n = size(scan%point)
  
  call re_allocate (curve%x_line, n)
  call re_allocate (curve%y_line, n)
  call re_allocate (curve%x_symb, n)
  call re_allocate (curve%y_symb, n)

  write(curve%legend_text, '(a, f6.2, a)') '\gd:', 100*da%pz(k), ' %'

  if (curve%data_type == 'dynamic_aperture') then
    curve%x_line(:) = scan%point(:)%x
    curve%y_line(:) = scan%point(:)%y
  else  ! "dynamic_aperture_ref0"
    curve%x_line(:) = scan%point(:)%x + scan%ref_orb%vec(1)
    curve%y_line(:) = scan%point(:)%y + scan%ref_orb%vec(3)
  endif

  curve%x_symb = curve%x_line
  curve%y_symb = curve%y_line
enddo

end subroutine tao_graph_dynamic_aperture_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! routine for drawing the physical aperture for dynamic aperture plots.

subroutine tao_curve_beam_ellipse_setup (curve)

use wall3d_mod

type (tao_curve_struct) :: curve
type (tao_universe_struct), pointer :: u
type (tao_dynamic_aperture_struct), pointer :: da
type (ele_pointer_struct), allocatable :: eles(:)
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele

real(rp) :: phi, min_angle, max_angle, angle, xx, yy, dr_dtheta
integer :: i, n_loc
integer, parameter :: n = 25 ! points per elliptical quadrant
character(*), parameter :: r_name = 'tao_curve_beam_ellipse_setup'
logical err

! Try to find a lattice element at s = 0 that has apertures defined

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
lat => u%model%lat
da => u%dynamic_aperture

if (da%a_emit <= 0 .or. da%b_emit <= 0) then
  curve%valid = .false.
  curve%why_invalid = ''  ! Don't want any error messages printed since this is not a real error.
  call out_io (s_warn$, r_name, 'a_emit and/or b_emit not set in tao_dynamic_aperture namelist. Will not draw beam ellipse.')
  return
endif

if (da%param%start_ele == '') then
  ele => lat%ele(0)
else
  call lat_ele_locator(da%param%start_ele, lat, eles, n_loc)
  ele => eles(1)%ele
endif

!

if (curve%data_type == 'beam_ellipse_full') then
  min_angle = 0
  max_angle = twopi
else
  min_angle = u%dynamic_aperture%param%min_angle
  max_angle = u%dynamic_aperture%param%max_angle
endif

xx = da%ellipse_scale * sqrt(ele%a%beta * da%a_emit)
yy = da%ellipse_scale * sqrt(ele%b%beta * da%b_emit)

call re_allocate (curve%x_line, 101)
call re_allocate (curve%y_line, 101)
if (allocated(curve%x_symb)) deallocate (curve%x_symb, curve%y_symb)

do i = 0, 100
  angle = (max_angle-min_angle) * i / 100.0_rp + min_angle
  curve%x_line(i+1) = xx * cos(angle) 
  curve%y_line(i+1) = yy * sin(angle) 
enddo

curve%valid = .true.

end subroutine tao_curve_beam_ellipse_setup


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_histogram_setup (plot, graph)

use bin_mod

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ref_ele, track_ele
type (beam_struct), pointer :: beam
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1
type (coord_struct), pointer :: p(:)
type (bin_struct) :: bins

real(rp) v_mat(4,4), v_inv_mat(4,4), g_mat(4,4), g_inv_mat(4,4)
real(rp) mat4(4,4), sigma_mat(4,4), theta, theta_xy, rx, ry, phi
real(rp) emit_a, emit_b, bin_shift
real(rp), allocatable :: data(:), weight(:)

integer k, n, m, ib, ix1_ax, ix, i
integer, allocatable :: number_in_bin(:)

logical err, same_uni

character(40) name
character(40) :: r_name = 'tao_graph_histogram_setup'

! Valid?

if (size(graph%curve) == 0) then
  graph%is_valid = .false.
  return
endif

! Set up the graph suffix

same_uni = .true.
ix = tao_universe_index(tao_curve_ix_uni(graph%curve(1)))
do k = 2, size(graph%curve)
  curve => graph%curve(k)
  if (tao_universe_index(tao_curve_ix_uni(curve)) /= ix) same_uni = .false.
enddo

graph%title_suffix = ''
do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
  if (.not. associated(u)) then
    call tao_set_curve_invalid (curve, 'NO ASSOCIATED UNIVERSE')
    cycle
  endif

  if (.not. u%is_on) then
    call tao_set_curve_invalid (curve, 'NO UNIVERSE IS TURNED OFF')
    cycle
  endif

  ref_ele => tao_curve_ele_ref(curve, .false.)
  track_ele => tao_curve_ele_ref(curve, .true.)

  if (.not. associated(track_ele)) then
    call tao_set_curve_invalid (curve, 'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name)
    cycle
  endif

  name = curve%ele_ref_name
  if (name == '') name = ref_ele%name
  if (same_uni) then
    write (graph%title_suffix, '(2a, i0, 3a)') trim(graph%title_suffix), &
                                '[', ref_ele%ix_ele, ': ', trim(name), ']'
  else
    write (graph%title_suffix, '(2a, i0, a, i0, 3a)') trim(graph%title_suffix), &
            '[', u%ix_uni, '@', ref_ele%ix_ele, ': ', trim(name), ']'
  endif
enddo

! loop over all curves

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))

  ! fill the data array

  if (allocated (curve%ix_symb)) deallocate (curve%ix_symb, curve%x_symb, curve%y_symb)
  if (allocated (curve%x_line))  deallocate (curve%x_line, curve%y_line)

  if (curve%data_source == 'beam') then
    track_ele => tao_curve_ele_ref(curve, .true.)
    beam => u%model_branch(track_ele%ix_branch)%ele(track_ele%ix_ele)%beam
    if (.not. allocated(beam%bunch)) then
      call out_io (s_warn$, r_name, 'NO ALLOCATED BEAM WITH PHASE_SPACE PLOTTING.')
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      return
    endif

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1, size(beam%bunch)
        n = n + count(beam%bunch(ib)%particle%state == alive$)
      enddo
    else
      n = count(beam%bunch(curve%ix_bunch)%particle%state == alive$)
    endif

    allocate (data(n))
    if (curve%hist%weight_by_charge) allocate (weight(n))
    
    n = 0
    do ib = 1, size(beam%bunch)
      p => beam%bunch(ib)%particle
      m = count(p%state == alive$)
      call tao_particle_data_value (curve%data_type, p, scratch%axis1, err, ref_ele, ib)
      data(n+1:n+m) = pack(scratch%axis1, mask = (p%state == alive$))
      if (curve%hist%weight_by_charge) weight(n+1:n+m) = pack(p%charge, mask = (p%state == alive$))
      n = n + m
    enddo

  ! Unrecognized

  else
    call tao_set_curve_invalid (curve, 'INVALID CURVE%DATA_SOURCE: ' // curve%data_source)
    cycle
  endif

  ! Bin the data

  curve%hist%minimum = minval(data)
  curve%hist%maximum = maxval(data)

  ! Select the number of bins if needed
  if (curve%hist%number == 0) then
    if (curve%hist%width == 0) then
      curve%hist%number = n_bins_automatic(size((data)))
    else
      ! Set the number according to this width, and set a new maximum
      curve%hist%number = nint((curve%hist%maximum - curve%hist%minimum)/curve%hist%width)
      curve%hist%maximum = curve%hist%number*curve%hist%width + curve%hist%minimum
    endif
  else
    ! Number is set, now set the width
    curve%hist%width = (curve%hist%maximum - curve%hist%minimum)/curve%hist%number
  endif
   
  ! Shift for the center bin
  bin_shift = curve%hist%center - bin_x_center(bin_index(curve%hist%center, curve%hist%minimum, curve%hist%width), &
                                                                                 curve%hist%minimum, curve%hist%width)
  if (bin_shift > 0) then
    curve%hist%maximum =  curve%hist%maximum + 2 * bin_shift
  else
    curve%hist%minimum =  curve%hist%minimum + 2 * bin_shift
  endif
                            
  curve%hist%width = (curve%hist%maximum - curve%hist%minimum)/(curve%hist%number - 1)
                            
  if (curve%hist%weight_by_charge) then
    bins = bin_data (data, weight = weight, min = curve%hist%minimum, max = curve%hist%maximum, n_bins = curve%hist%number)
  else
    bins = bin_data (data, min = curve%hist%minimum, max = curve%hist%maximum, n_bins = curve%hist%number)
           
  endif
  ! Set width actually used
 
  curve%hist%width = bins%delta
  
  allocate(curve%x_line(bins%n))
  allocate(curve%y_line(bins%n))
  do i=1, bins%n
    curve%x_line(i) = bin_x_center(i, bins%min, bins%delta)
    curve%y_line(i) = bins%count(i)
  enddo
  if (curve%hist%density_normalized) curve%y_line = curve%y_line/bins%delta
  
enddo

end subroutine tao_graph_histogram_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function tao_phase_space_axis_index (data_type, err) result (ix_axis)
!
! Routine to calculate the phase space axis index for a given data type.
!
! Input:
!   data_type   -- character(*): Type of data.
!   err         -- logical: Set True if there is an error.
!
! Output:
!   ix_axis     -- integer: Axis index.
!-

function tao_phase_space_axis_index (data_type, err) result (ix_axis)

implicit none

integer ix_axis

logical :: err

character(*) data_type
character(*), parameter :: r_name = 'tao_phase_space_axis_index'

!

err = .false.

select case (data_type)
case ('x');   ix_axis = 1
case ('px');  ix_axis = 2
case ('y');   ix_axis = 3
case ('py');  ix_axis = 4
case ('z');   ix_axis = 5
case ('pz');  ix_axis = 6

case default
  call out_io (s_warn$, r_name, 'BAD PHASE_SPACE DATA_TYPE: ' // data_type, &
                                'SHOULD BE ONE OF "x", "px", "y", "py", "z", OR "pz".')
  err = .true.
end select

end function tao_phase_space_axis_index

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_particle_data_value (data_type, p, value, err, ele, ix_bunch)
!
! Routine to calculate the value array of a data_type for an array of particles.
!
! Input:
!   data_type   -- character(*): Type of data.
!   p(:)        -- coord_struct, Array of particles containing the data.
!   ele         -- ele_struct: Needed for "Ja" evaluation.
!   ix_bunch    -- integer: Bunch index.
!
! Output:
!   value(:)    -- real(rp): Array of values.
!   err         -- logical: Set True if there is an error. False otherwise.
!-

subroutine tao_particle_data_value (data_type, p, value, err, ele, ix_bunch)

implicit none

type (coord_struct), target :: p(:)
type (coord_struct) :: p1
type (ele_struct) :: ele

real(rp), allocatable :: value(:)

integer ix_bunch
integer i

logical :: err

character(*) data_type
character(*), parameter :: r_name = 'tao_particle_data_value'

!

err = .false.
call re_allocate (value, size(p))

select case (data_type)
case ('x');   value = p%vec(1)
case ('px');  value = p%vec(2)
case ('y');   value = p%vec(3)
case ('py');  value = p%vec(4)
case ('z');   value = p%vec(5)
case ('pz');  value = p%vec(6)
case ('intensity_x'); value = p%field(1)**2
case ('intensity_y'); value = p%field(2)**2
case ('phase_x');     value = p%phase(1)
case ('phase_y');     value = p%phase(2)
case ('t', 'time');   value = p%t
case ('bunch_index');    value = ix_bunch

case ('intensity')
  value = p%field(1)**2 + p%field(2)**2
  
case ('Ja')
  do i=1, size(p)
    call convert_coords('LAB', p(i), ele, 'ACTION-ANGLE', p1)
    value(i) = p1%vec(1)
  enddo

case ('Jb')
  do i=1, size(p)
    call convert_coords('LAB', p(i), ele, 'ACTION-ANGLE', p1)
    value(i) = p1%vec(3)
  enddo

case ('energy')
  do i=1, size(p)
    call convert_pc_to((1 + p(i)%vec(6)) * p(i)%p0c, p(i)%species, e_tot =  value(i))
  enddo
  
case default
  call out_io (s_warn$, r_name, 'BAD PHASE_SPACE CURVE DATA_TYPE: ' // data_type)
  err = .true.
end select

end subroutine tao_particle_data_value

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_data_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), target :: branch_curve
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

integer n, ic, ib, n0_line, n0_symb
logical err

!

if (allocated(graph%curve)) then
  if (all(graph%curve%component == graph%curve(1)%component) .and. graph%curve(1)%component /= '') then
    graph%title_suffix = graph%curve(1)%component
  else
    graph%title_suffix = ''
  endif

  if (size(s%u) > 1) then
    if (all(graph%curve%ix_universe == graph%curve(1)%ix_universe) .and. graph%curve(1)%ix_universe /= -1) then
      graph%title_suffix = trim(graph%title_suffix) // ' uni:' // int_str(graph%curve(1)%ix_universe)
    elseif (graph%ix_universe /= -1) then
      graph%title_suffix = trim(graph%title_suffix) // ' uni:' // int_str(graph%ix_universe)
    else
      graph%title_suffix = trim(graph%title_suffix) // ' uni:' // int_str(s%global%default_universe)
    endif
  endif

else
  graph%title_suffix = ''
endif

if (graph%title_suffix /= '') graph%title_suffix = '[' // trim(adjustl(graph%title_suffix)) // ']'

! Attach x-axis type to title suffix if needed.
! Needed %label is blank and %draw_label = F.
! Note: if %label is blank and %draw_label = T then the x-axis_type is printed elsewhere.
 
if (graph%x%label == '' .and. .not. graph%x%draw_label) then
  if (plot%x_axis_type == "lat" .or. plot%x_axis_type == "var") then
    graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: ' // &
              trim(plot%x_axis_type) // '::' // graph%curve(1)%data_type_x
  else
    graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: ' // plot%x_axis_type
  endif
endif

! Loop over all curves in the graph

if (allocated(graph%curve)) then
  do ic = 1, size(graph%curve)
    curve => graph%curve(ic)
    u => tao_pointer_to_universe (tao_curve_ix_uni(graph%curve(ic)))
    if (.not. tao_curve_check_universe (curve, u)) cycle
    call tao_curve_data_setup (plot, graph, curve)
  enddo
endif

end subroutine tao_graph_data_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_curve_data_setup (plot, graph, curve)

use swap_mod

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), target :: curve
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: model_lat, base_lat
type (tao_ele_shape_struct), pointer :: ele_shape
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (ele_struct), pointer :: ele, ele1, ele2, slave
type (branch_struct), pointer :: branch
type (tao_curve_array_struct), allocatable :: curves(:)
type (all_pointer_struct) var_ptr

real(rp) f, eps, gs, l_tot, s0, s1, x_max, x_min, val, val0, dx, limit, len_branch
real(rp), allocatable :: value_arr(:), x_arr(:), y_arr(:)


integer ii, k, m, n, n_dat, n2_dat, ib, ie, jj, iv, ic
integer ix, ir, jg, i, j, ix_this, ix_uni, ix1, ix2, n_curve_pts, ix_slave
integer, allocatable :: xx_arr(:)

logical err, err_flag, smooth_curve, found, zero_average_phase, ok
logical straight_line_between_syms, valid, in_graph
logical, allocatable :: good(:), this_u(:)

character(200) data_type, name
character(100) str
character(16) data_source, dflt_index
character(*), parameter :: r_name = 'tao_curve_data_setup'

!

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

call re_allocate_eles (scratch%eles, 1, exact = .true.)

if (allocated(curve%x_line)) deallocate (curve%x_line, curve%y_line)
if (allocated(curve%y2_line)) deallocate (curve%y2_line)
if (allocated(curve%ix_line)) deallocate (curve%ix_line)
if (allocated(curve%x_symb)) deallocate (curve%x_symb, curve%y_symb)
if (allocated(curve%z_symb)) deallocate (curve%z_symb)
if (allocated(curve%ix_symb)) deallocate (curve%ix_symb)
if (allocated(curve%err_symb)) deallocate (curve%err_symb)
if (allocated(curve%symb_size)) deallocate (curve%symb_size)

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
if (.not. tao_curve_check_universe(curve, u)) return

model_lat => u%model%lat
base_lat => u%base%lat

ib = tao_branch_index(curve%ix_branch)
branch => model_lat%branch(ib)

if (curve%ele_ref_name == '') then
  zero_average_phase = .true.
else
  zero_average_phase = .false.
  ele => tao_curve_ele_ref(curve, .true.)
  if (.not. associated(ele)) then
    call tao_set_curve_invalid (curve, 'CANNOT LOCATE ELEMENT: ' // trim(curve%ele_ref_name))
    return
  endif
endif

!----------------------------------------------------------------------------
! Calculate where the symbols are to be drawn on the graph.

data_source = curve%data_source

if (plot%x_axis_type == 'lat' .or. plot%x_axis_type == 'var') data_source = 'plot_x_axis_var'
if (substr(curve%data_type,1,11) == 'expression:' .and. (data_source == 'dat' .or. data_source == 'var')) data_source = 'expression'

select case (data_source)

!----------------------------------------------------------------------------
! Case: expression

case ('expression')

  if (plot%x_axis_type == 'index') then
    n_dat = nint(graph%x%eval_max) - nint(graph%x%eval_min) + 1

    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
    call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

    n_dat = 0
    do i = nint(graph%x%eval_min), nint(graph%x%eval_max)
      write (dflt_index, '(i0)') i
      call tao_evaluate_expression  (curve%data_type(12:), 0, graph%draw_only_good_user_data_or_vars, value_arr, &
                   err, .false., scratch%info, scratch%stack, curve%component, curve%data_source, dflt_dat_or_var_index = dflt_index)
      if (err .or. .not. scratch%info(1)%good) cycle
      n_dat = n_dat + 1

      curve%x_symb(n_dat) = i
      curve%y_symb(n_dat) = value_arr(1)
      curve%ix_symb(n_dat) = 0
    enddo

  else
    call tao_evaluate_expression  (curve%data_type(12:), 0, graph%draw_only_good_user_data_or_vars, value_arr, &
                                            err, .true., scratch%info, scratch%stack, curve%component, curve%data_source)
    if (err) return
    n_dat = count(scratch%info%good)
    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
    call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

    n_dat = 0
    do i = 1, size(value_arr)
      if (.not. scratch%info(i)%good) cycle
      n_dat = n_dat + 1
      curve%y_symb(n_dat) = value_arr(i)

      dx = (graph%x%eval_max - graph%x%eval_min) / 100.0_rp

      select case (plot%x_axis_type)
      case ('s')
        if (scratch%info(i)%s < graph%x%eval_min - dx .or. scratch%info(i)%s > graph%x%eval_max + dx) cycle
        curve%x_symb(n_dat) = scratch%info(i)%s
      case ('ele_index')
        if (.not. associated(scratch%info(i)%ele)) cycle
        if (scratch%info(i)%ele%ix_ele < graph%x%eval_min - dx .or. scratch%info(i)%ele%ix_ele > graph%x%eval_max + dx) cycle
        curve%x_symb(n_dat) = scratch%info(i)%ele%ix_ele
      end select
    enddo
  endif

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
  call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

  !

  if (curve%draw_line) then
    call re_allocate (curve%x_line,  n_dat) ! allocate space for the data
    call re_allocate (curve%y_line,  n_dat) ! allocate space for the data
    curve%x_line = curve%x_symb
    curve%y_line = curve%y_symb
  else
    if (allocated (curve%x_line)) deallocate (curve%x_line, curve%y_line)
  endif

!----------------------------------------------------------------------------
! Case: x-axis uses a variable.

case ('plot_x_axis_var')

  call re_allocate (curve%ix_symb, n_curve_pts)
  call re_allocate (curve%x_symb, n_curve_pts)
  call re_allocate (curve%y_symb, n_curve_pts)
  call re_allocate (curve%x_line, n_curve_pts)
  call re_allocate (curve%y_line, n_curve_pts)
  var_ptr = all_pointer_struct()

  if (plot%x_axis_type == 'lat') then

    call tao_pick_universe (curve%data_type_x, name, this_u, err, ix_uni)
    if (err .or. count(this_u) /= 1) then
      call tao_set_curve_invalid (curve, 'BAD UNIVERSE CONSTRUCT IN CURVE%DATA_TYPE_X: ' //curve%data_type_x)
      return
    endif

    call upcase_string(name)
    ix1 = index(name, '[')
    ix2 = index(name, ']')
    if (ix1 == 0 .or. ix2 == 0 .or. ix2 /= len_trim(name)) then
      call tao_set_curve_invalid (curve, 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // curve%data_type_x)
      return
    endif

    u => tao_pointer_to_universe(ix_uni)
    call pointers_to_attribute (u%model%lat, name(1:ix1-1), name(ix1+1:ix2-1), .true., scratch%attribs, &
                                                                                 err, eles = scratch%eles)
    if (err .or. size(scratch%attribs) /= 1) then
      call tao_set_curve_invalid (curve, 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // curve%data_type_x)
      return
    endif
    var_ptr%r => scratch%attribs(1)%r

  else  ! x_axis_type == 'var'
    call tao_find_var (err, curve%data_type_x, v_array = scratch%var_array)
    if (err .or. size(scratch%var_array) /= 1) then
      call tao_set_curve_invalid (curve, 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // curve%data_type_x)
      return
    endif
    var_ptr%r => scratch%var_array(1)%v%model_value
  endif

  ! Get datum values as a function of the variable

  val0 = var_ptr%r

  j = 0
  do i = 1, n_curve_pts 
    val = graph%x%eval_min + (graph%x%eval_max - graph%x%eval_min) * (i - 1.0_rp) / (n_curve_pts - 1)
    if (plot%x_axis_type == 'lat')then
      var_ptr%r = val
      call tao_set_flags_for_changed_attribute (u, name(1:ix1-1), scratch%eles(1)%ele, var_ptr)
      s%u(ix_uni)%calc%lattice = .true.
    else
      call tao_set_var_model_value (scratch%var_array(1)%v, val)
    endif
    call tao_lattice_calc (valid, .false.)

    call tao_evaluate_expression (curve%data_type, 0, .false., value_arr, err, &
                          dflt_component = curve%component, dflt_source = curve%data_source)

    if (.not. valid .or. err .or. size(value_arr) /= 1) cycle
    j = j + 1
    curve%x_symb(j) = val      
    curve%y_symb(j) = value_arr(1)
    curve%x_line(j) = val      
    curve%y_line(j) = value_arr(1)
  enddo

  if (j == 0) then
    if (.not. valid) then
      call tao_set_curve_invalid (curve, 'PROBLEM CALCULATING LATTICE PARAMETERS (INSTABILITY?) FOR: ' // curve%data_type, .true.)
    elseif(err .or. size(value_arr) /= 1) then
      call tao_set_curve_invalid (curve, 'PROBLEM EVALUATING DATA FOR: ' // curve%data_type, .true.)
    endif
  endif

  call re_allocate(curve%x_symb, j)
  call re_allocate(curve%y_symb, j)
  call re_allocate(curve%x_line, j)
  call re_allocate(curve%y_line, j)

  ! Reset

  if (plot%x_axis_type == 'lat')then
    var_ptr%r = val0
    call tao_set_flags_for_changed_attribute (u, name(1:ix1-1), scratch%eles(1)%ele, var_ptr)
    s%u(ix_uni)%calc%lattice = .true.
  else
    call tao_set_var_model_value (scratch%var_array(1)%v, val0)
  endif
  call tao_lattice_calc (valid)

!----------------------------------------------------------------------------
! Case: data_source is a data_array

case ('data')

  ! Calculate values

  call tao_data_type_substitute (curve%data_type, data_type, curve, graph)
  call tao_evaluate_expression  (data_type, 0, graph%draw_only_good_user_data_or_vars, value_arr, err, &
                          stack = scratch%stack, dflt_component = curve%component, dflt_source = 'data', dflt_uni = u%ix_uni)
  if (err) then
    call tao_set_curve_invalid (curve, 'CANNOT FIND DATA CORRESPONDING: ' // data_type)
    return
  end if

  ! point d1_array to the data to be plotted

  do i = 1, size(scratch%stack)
    if (scratch%stack(i)%type == data_num$) exit
    if (i == size(scratch%stack)) then
      call tao_set_curve_invalid (curve, 'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type)
      return
    endif
  enddo

  call tao_find_data (err, scratch%stack(i)%name, d2_array, scratch%d1_array, ix_uni = tao_curve_ix_uni(curve))
  if (err .or. size(scratch%d1_array) /= 1) then
    call tao_set_curve_invalid (curve, 'CANNOT FIND VALID DATA ARRAY TO PLOT CURVE: ' // curve%data_type)
    return
  endif

  d1_ptr => scratch%d1_array(1)%d1
  if (index(curve%data_type, 'phase') /= 0) then
    if (all(d1_ptr%d(:)%ele_ref_name == '')) then
      zero_average_phase = .true.
    else
      zero_average_phase = .false.
    endif
  endif

  ! Set %good_plot True for all data that is within the x-axis limits.
  ! For a circular lattice "wrap around" at s = 0 may mean some data points show up twice.

  d1_ptr => scratch%d1_array(1)%d1
  d1_ptr%d%good_plot = .false.
  if (graph%x%min /= graph%x%max) then
    eps = 1e-6 * (graph%x%max - graph%x%min)
    if (plot%x_axis_type == 'index') then
      where (d1_ptr%d%ix_d1 > graph%x%min-eps .and. &
             d1_ptr%d%ix_d1 < graph%x%max+eps) d1_ptr%d%good_plot = .true.
    elseif (plot%x_axis_type == 'ele_index') then
      where (d1_ptr%d%ix_ele > graph%x%min-eps .and. &
             d1_ptr%d%ix_ele < graph%x%max+eps) d1_ptr%d%good_plot = .true.
    else ! s
      where (d1_ptr%d%s > graph%x%min-eps .and. &
             d1_ptr%d%s < graph%x%max+eps) d1_ptr%d%good_plot = .true.
      ! Wrap around can be problematical in some cases. For example, with beam tracking data.
      ! So here use min+eps and max-eps for the cutoff instead of min-eps and max+eps
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then 
        l_tot = branch%param%total_length
        where (d1_ptr%d%s-l_tot > graph%x%min+eps .and. &
               d1_ptr%d%s-l_tot < graph%x%max-eps) d1_ptr%d%good_plot = .true.
        where (d1_ptr%d%s+l_tot > graph%x%min+eps .and. &
               d1_ptr%d%s+l_tot < graph%x%max-eps) d1_ptr%d%good_plot = .true.
      endif

    endif
  endif

  ! make sure %useit_plot up-to-date & count the number of data points

  call tao_data_useit_plot_calc (curve, graph, d1_ptr%d, plot%x_axis_type == 's', curve%why_invalid) 
  n_dat = count (d1_ptr%d%useit_plot)       

  ! resize the curve data arrays

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%x_symb, n_dat) ! allocate space for the data

  ! 

  curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)
  curve%y_symb  = pack(value_arr, mask = d1_ptr%d%useit_plot)

  if (curve%draw_error_bars) then
    call re_allocate (curve%err_symb, n_dat)
    curve%err_symb = pack(d1_ptr%d%error_rms, mask = d1_ptr%d%useit_plot)
  endif

  if (plot%x_axis_type == 'index') then
    curve%x_symb = curve%ix_symb
  elseif (plot%x_axis_type == 'ele_index') then
    curve%x_symb = d1_ptr%d(curve%ix_symb)%ix_ele
  elseif (plot%x_axis_type == 's') then
    curve%x_symb = d1_ptr%d(curve%ix_symb)%s
    ! If there is a wrap-around then reorder data
    if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
      do i = 1, n_dat
        if (curve%x_symb(i) > graph%x%max+eps) curve%x_symb(i) = curve%x_symb(i)-l_tot
        if (curve%x_symb(i) < graph%x%min-eps) curve%x_symb(i) = curve%x_symb(i)+l_tot
      enddo
    endif
    ! Super lords will be out of order so reorder in increasing s.
    do i = 2, n_dat
      do j = i, 2, -1
        if (curve%x_symb(j-1) > curve%x_symb(j)) then
          call swap(curve%x_symb(j-1), curve%x_symb(j))
          call swap(curve%y_symb(j-1), curve%y_symb(j))
          call swap(curve%ix_symb(j-1), curve%ix_symb(j))
        else
          exit
        endif
      enddo
    enddo

  else
    call tao_set_curve_invalid (curve, 'UNKNOWN X_AXIS_TYPE: ' // quote(plot%x_axis_type))
    return
  endif

!----------------------------------------------------------------------------
! Case: data_source is another curve.

case ('curve')

  if (plot%x_axis_type /= 'curve') then
    call tao_set_curve_invalid (curve, 'IF data_source IS SET TO "CURVE" THEN THE PLOT''S x_axis_type MUST BE SET TO "CURVE".', .true.)
    return
  endif

  call tao_find_plots (err, curve%data_type, 'REGION', curve = curves, blank_means_all = .true., only_visible = .false.)
  if (err .or. size(curves) < 2) then
    call tao_set_curve_invalid (curve, 'CANNOT LOCATE CURVES CORRESPONDING TO: ' // curve%data_type)
    return
  endif

  if (size(curves) > 2) then
    call tao_set_curve_invalid (curve, 'TOO MANY CURVES ASSOCIATED WITH: ' // curve%data_type)
    return
  endif

  ! Line setup

  if (curve%draw_line) then
    if (.not. allocated(curves(1)%c%y_line) .or. .not. allocated (curves(2)%c%y_line)) then
      call tao_set_curve_invalid (curve, 'DRAW_LINE = T BUT ASSOCIATED CURVES HAVE NO LINE DATA.')
      return
    endif

    n = size(curves(1)%c%y_line)
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)
    curve%x_line = curves(1)%c%y_line
    curve%y_line = curves(2)%c%y_line
  else
    if (allocated (curve%x_line)) deallocate (curve%x_line, curve%y_line)
  endif

  ! Curve setup

  if (curve%draw_symbols) then
    if (.not. allocated(curves(1)%c%y_symb) .or. .not. allocated (curves(2)%c%y_symb)) then
      call tao_set_curve_invalid (curve, 'DRAW_SYMBOLS = T BUT ASSOCIATED CURVES HAVE NO SYMBOL DATA.')
      return
    endif

    n = size(curves(1)%c%y_symb)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)
    call re_allocate (curve%ix_symb, n)
    curve%x_symb  = curves(1)%c%y_symb
    curve%y_symb  = curves(2)%c%y_symb
    curve%ix_symb = curves(1)%c%ix_symb
  else
    if (allocated (curve%x_symb)) deallocate (curve%x_symb, curve%y_symb, curve%ix_symb)
  endif

!----------------------------------------------------------------------------
! Case: data_source is a var_array

case ('var')

  call tao_find_var (err, curve%data_type, scratch%v1_array)
  if (err .or. size(scratch%v1_array) /= 1) return
  v1_ptr => scratch%v1_array(1)%v1

  ! find which universe we're viewing
  ix_this = -1
  v_loop: do iv = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
    v_ptr => v1_ptr%v(iv)
    if (.not. v_ptr%exists) cycle
    do jj = 1, size(v_ptr%slave)
      if (v_ptr%slave(jj)%ix_uni .eq. s%global%default_universe) then
        ix_this = jj
        exit v_loop
      endif
    enddo
  enddo v_loop
  if (ix_this .eq. -1) then
    call out_io (s_warn$, r_name, "This variable doesn't point to the currently displayed  universe.")
    return
  endif

  v1_ptr%v%good_plot = .true.
  if (graph%x%min /= graph%x%max) then
    eps = 1e-4 * (graph%x%max - graph%x%min)
    if (plot%x_axis_type == 'index') then
      where (v1_ptr%v%ix_v1 < graph%x%min-eps) v1_ptr%v%good_plot = .false.
      where (v1_ptr%v%ix_v1 > graph%x%max+eps) v1_ptr%v%good_plot = .false.
    elseif (plot%x_axis_type == 'ele_index') then
      do jj = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
        if (v1_ptr%v(jj)%slave(ix_this)%ix_ele < graph%x%min-eps) v1_ptr%v%good_plot = .false.
        if (v1_ptr%v(jj)%slave(ix_this)%ix_ele > graph%x%max+eps) v1_ptr%v%good_plot = .false.
      enddo
    else
      where (v1_ptr%v%s < graph%x%min-eps) v1_ptr%v%good_plot = .false.
      where (v1_ptr%v%s > graph%x%max+eps) v1_ptr%v%good_plot = .false.
    endif
  endif

  call tao_var_useit_plot_calc (graph, v1_ptr%v) ! make sure %useit_plot up-to-date
  n_dat = count(v1_ptr%v%useit_plot)             ! count the number of data points

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data
  call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data

  curve%ix_symb = pack(v1_ptr%v%ix_v1, mask = v1_ptr%v%useit_plot)

  graph%x%label = plot%x_axis_type

  if (plot%x_axis_type == 'index') then
    curve%x_symb = curve%ix_symb
  elseif (plot%x_axis_type == 'ele_index') then
    do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
      curve%x_symb(jj) = v1_ptr%v(curve%ix_symb(jj))%slave(ix_this)%ix_ele
    enddo
  elseif (plot%x_axis_type == 's') then
    do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
      ele => branch%ele(v1_ptr%v(curve%ix_symb(jj))%slave(ix_this)%ix_ele)
      if (ele%lord_status == multipass_lord$) ele => pointer_to_slave(ele, 1)
      curve%x_symb(jj) = ele%s
    enddo
  endif

  ! calculate the y-axis data point values.

  data_type = trim(curve%data_type) // '|' // trim(curve%component)
  call tao_evaluate_expression (data_type, 0, graph%draw_only_good_user_data_or_vars, value_arr, err)
  if (err) then
    call out_io (s_warn$, r_name, 'BAD CURVE DATA_TYPE: ' // data_type)
    return
  end if

  curve%y_symb = pack(value_arr, mask = v1_ptr%v%useit_plot)

!----------------------------------------------------------------------------
! Case: data_source is from lattice, or beam

case ('lat', 'beam')

  ! Find how many symbol points there are...
  ! Here 'index' and 'ele_index' mean the same thing.

  select case (plot%x_axis_type)
  case ('index', 'ele_index')
    x_min = 1
    x_max = branch%n_ele_track
    if (graph%x%min /= graph%x%max) then
      x_min = max(x_min, graph%x%min)
      x_max = min(x_max, graph%x%max)
    endif 
    n_dat = max(0, nint(x_max+1-x_min))
    call re_allocate_eles (scratch%eles, n_dat, exact = .true.)
    do i = 1, n_dat
      scratch%eles(i)%ele => pointer_to_ele (model_lat, nint(i+x_min-1), ib)
    enddo

  ! x_axis_type == 's':

  case ('s')
    ! Symbols are to be put at the ends of displayed elements in the lat_layout
    eps = 1e-4 * (graph%x%max - graph%x%min)       ! a small number
    branch%ele%logic = .false.                     ! Mark if ele is in the graph

    ! Mark all eles in branch if they match a shape.
    do i = 0, branch%n_ele_track
      ele => branch%ele(i)
      ele_shape => tao_pointer_to_ele_shape (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape)
      if (.not. associated(ele_shape)) cycle
      if (.not. ele_shape%draw) cycle
      call find_element_ends (ele, ele1, ele2)
      ele1%logic = .true.
      ele2%logic = .true.
    enddo

    ! Mark slaves of lord elements that match a shape.
    do i = model_lat%n_ele_track+1, model_lat%n_ele_max
      ele => model_lat%ele(i)
      ele_shape => tao_pointer_to_ele_shape (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape)
      if (.not. associated(ele_shape)) cycle
      if (.not. ele_shape%draw) cycle
      if (ele%lord_status == multipass_lord$) then
        do j = 1, ele%n_slave
          slave => pointer_to_slave (ele, j)
          call find_element_ends (slave, ele1, ele2)
          ele1%logic = .true.
          ele2%logic = .true.
        enddo
      else
        call find_element_ends (ele, ele1, ele2)
        ele1%logic = .true.
        ele2%logic = .true.
      endif
    enddo

    ! Now unmark all elements in the branch that are not within the graph boundries.
    do i = 0, branch%n_ele_track
      ele => branch%ele(i)
      if (.not. ele%logic) cycle
      if (graph%x%min == graph%x%max) cycle
      s0 = ele%s_start
      s1 = ele%s
      in_graph = (s0 >= graph%x%min-eps) .and. (s1 <= graph%x%max+eps)
      l_tot = branch%param%total_length
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) in_graph = in_graph .or. &
                      (s0-l_tot >= graph%x%min-eps) .and. (s1-l_tot <= graph%x%max+eps)
      ele%logic = ele%logic .and. in_graph                                 
    enddo

    ! Allocate eles(:) array and set the eles(:)%ele pointers to point to the marked elements.
    n_dat = count (branch%ele(:)%logic)
    call re_allocate_eles (scratch%eles, n_dat, exact = .true.)

    n = 0
    do i = 0, branch%n_ele_max
      ele => branch%ele(i)
      if (.not. ele%logic) cycle
      n = n + 1
      scratch%eles(n)%ele => ele 
    enddo      

  ! Error for x_axis_type unrecognized.

  case default
    call out_io (s_error$, r_name, 'BAD PLOT%X_AXIS_TYPE: ' // plot%x_axis_type)
    return
  end select

  if (curve%draw_symbols .and. substr(curve%data_type,1,9) /= 'aperture.') then
    call tao_curve_datum_calc (scratch%eles, plot, curve, 'SYMBOL')
    if (.not. curve%valid) return
  endif

!----------------------------------------------------------------------------
! Case: Aperture

case ('aperture')

  branch => u%model%lat%branch(ib)
  call re_allocate (xx_arr, 100)
  call re_allocate (x_arr, 100)
  call re_allocate (y_arr, 100)
  ir = 0

  do i = 0, branch%n_ele_track
    ele => branch%ele(i)

    select case (curve%data_type)
    case ('+x'); limit = -ele%value(x1_limit$)
    case ('-x'); limit =  ele%value(x2_limit$)
    case ('+y'); limit = -ele%value(y1_limit$)
    case ('-y'); limit =  ele%value(y2_limit$)
    case default
      call out_io (s_error$, r_name, &
              'BAD CURVE%DATA_TYPE VALUE FOR CURVE WITH %DATA_SOURCE SET TO "APERTURE": ' // quote(curve%data_type), &
              'SHOULD BE ONE OF: "+x", "-x", "+y", or "-y".')
    end select

    if (limit == 0) cycle

    if (ir + 2 > size(x_arr)) then
      call re_allocate (xx_arr, 2*size(xx_arr))
      call re_allocate (x_arr, 2*size(x_arr))
      call re_allocate (y_arr, 2*size(y_arr))
    endif    

    if (at_this_ele_end(physical_ele_end(first_track_edge$, coord_struct(), ele%orientation), ele%aperture_at)) then
      ir = ir + 1
      y_arr(ir) = limit
      xx_arr(ir) = i - 1 
      x_arr(ir) = ele%s_start
    endif

    if (at_this_ele_end(physical_ele_end(second_track_edge$, coord_struct(), ele%orientation), ele%aperture_at)) then
      ir = ir + 1
      y_arr(ir) = limit
      xx_arr(ir) = i
      x_arr(ir) = ele%s
    endif
  enddo

  select case (plot%x_axis_type)
  case ('s')
    len_branch = branch%param%total_length
    call tao_graph_s_min_max_calc(graph, branch, x_min, x_max)
    if (.not. graph%is_valid) return

    if (x_min < branch%ele(0)%s) then  ! Wrap case
      ix1 = bracket_index(x_min+len_branch-10*bmad_com%significant_length, x_arr(1:ir), 1, restrict = .true.) 
      if (x_max < branch%ele(0)%s) then
        ix2 = bracket_index(x_max+len_branch+10*bmad_com%significant_length, x_arr(1:ir), 1, restrict = .true.) + 1
      else
        ix2 = ir
      endif

      n_dat = ix2 - ix1 + 1
      call re_allocate(curve%x_symb, n_dat)
      call re_allocate(curve%y_symb, n_dat)
      curve%x_symb = x_arr(ix1:ix2)
      curve%y_symb = y_arr(ix1:ix2)

      if (x_max > branch%ele(0)%s) then
        ix2 = bracket_index(x_max+10*bmad_com%significant_length, x_arr(1:ir), 1, restrict = .true.) + 1
        n2_dat = ix2 + n_dat
        call re_allocate(curve%x_symb, n2_dat)
        call re_allocate(curve%y_symb, n2_dat)
        curve%x_symb(n_dat+1:n2_dat) = x_arr(1:ix2)
        curve%y_symb(n_dat+1:n2_dat) = y_arr(1:ix2)
      endif

    else  ! No wrap case.
      ix1 = bracket_index(x_min-10*bmad_com%significant_length, x_arr(1:ir), 1, restrict = .true.) 
      ix2 = bracket_index(x_max+10*bmad_com%significant_length, x_arr(1:ir), 1, restrict = .true.) + 1
      n_dat = ix2 - ix1 + 1
      call re_allocate(curve%x_symb, n_dat)
      call re_allocate(curve%y_symb, n_dat)
      curve%x_symb = x_arr(ix1:ix2)
      curve%y_symb = y_arr(ix1:ix2)
    endif


  ! Index x-axis.
  case default
    ix1 = bracket_index_int(nint(x_min), xx_arr(1:ir), 1, restrict = .true.)
    ix2 = bracket_index_int(nint(x_max), xx_arr(1:ir), 1, restrict = .true.) + 1
    n_dat = ix2 - ix1 + 1
    call re_allocate(curve%x_symb, n_dat)
    call re_allocate(curve%y_symb, n_dat)
    curve%x_symb = xx_arr(ix1:ix2)
    curve%y_symb = y_arr(ix1:ix2)
  end select

  call re_allocate(curve%x_line, size(curve%x_symb))
  call re_allocate(curve%y_line, size(curve%y_symb))
  curve%x_line = curve%x_symb
  curve%y_line = curve%y_symb
    
  return

!----------------------------------------------------------------------------
! Case: Bad data_source

case default
  call out_io (s_error$, r_name, 'UNKNOWN DATA_SOURCE: ' // curve%data_source)
  return
end select

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Line calc
! If the x-axis is by index or ele_index then these points are the same as the symbol points.
! That is, for x-axis = 'index' or 'ele_index' the line is piece-wise linear between the symbols.

if (curve%draw_line) then
  select case (plot%x_axis_type)
  case ('index', 'ele_index')
    call re_allocate (curve%y_line, size(curve%x_symb)) ! allocate space for the data
    call re_allocate (curve%x_line, size(curve%y_symb)) ! allocate space for the data
    curve%x_line = curve%x_symb
    curve%y_line = curve%y_symb

  ! If the axis is by s-value then, if possible, the line is a "smooth" curve with n_curve_pts points.

  case ('s')

    ! beam data_source is not interpolated.
    smooth_curve = (curve%data_source == 'lat' .and. curve%smooth_line_calc .and. .not. s%global%disable_smooth_line_calc)
    if (index(curve%data_type, 'emit.') /= 0) smooth_curve = .false.
    if (substr(curve%data_type,1,4) == 'bpm_') smooth_curve = .false.
    if (substr(curve%data_type,1,6) == 'bunch_') smooth_curve = .false.
    if (substr(curve%data_type,1,7) == 'smooth.') smooth_curve = .false.
    if (substr(curve%data_type,1,7) == 'spin_dn') smooth_curve = .false.
    if (substr(curve%data_type,1,15) == 'element_attrib.') smooth_curve = .false.
    if (substr(curve%data_type,1,6) == 'chrom.' .and. substr(curve%data_type,1,7) /= 'chrom.w') smooth_curve = .false.

    if (index(curve%component, 'meas') /= 0 .or. index(curve%component, 'ref') /= 0 .or. &
        curve%data_source == 'data' .or. curve%data_source == 'var') then
      straight_line_between_syms = .true.
      smooth_curve = .false.
    else
      straight_line_between_syms = .false.
    endif

    ! Smooth curves using expressions...

    if (smooth_curve) then
      ! Allocate data space. 
      ! Tracking is problematical if the step size is less than significant_length so adjust if needed.

      n_curve_pts = nint(min(1.0_rp*n_curve_pts, &
                          2+0.1_rp*u%model%lat%branch(ib)%param%total_length/bmad_com%significant_length))
      call re_allocate (curve%y_line, n_curve_pts) 
      call re_allocate (curve%x_line, n_curve_pts)
      call re_allocate (good, n_curve_pts) 
      good = .true.
      curve%y_line = 0

      call tao_split_component(curve%component, scratch%comp, err)
      if (err) then
        call tao_set_curve_invalid (curve, 'BAD CURVE COMPONENT EXPRESSION: ' // curve%component)
        return
      endif
      if (curve%component == '') then
        call tao_set_curve_invalid (curve, 'BLANK CURVE COMPONENT STRING.')
        return
      endif

      do m = 1, size(scratch%comp)
        select case (scratch%comp(m)%name)
        case ('') 
          cycle
        case ('model')
          call tao_calc_data_at_s_pts (u%model, curve, scratch%comp(m)%sign, good)
        case ('base')  
          call tao_calc_data_at_s_pts (u%base, curve, scratch%comp(m)%sign, good)
        case ('design')  
          call tao_calc_data_at_s_pts (u%design, curve, scratch%comp(m)%sign, good)
        case default
          call tao_set_curve_invalid (curve, 'BAD CURVE COMPONENT: ' // curve%component)
          return
        end select
      enddo

      if (substr(curve%data_type, 1, 3) == 'b0_' .or. substr(curve%data_type, 1, 3) == 'e0_') then
        ix = index(curve%legend_text, '[')
        if (ix == 0) ix = len_trim(curve%legend_text) + 2
        str = ''
        if (curve%orbit%x /= 0) str = trim(str) // ', x=' // real_str(curve%orbit%x, 5)
        if (curve%orbit%y /= 0) str = trim(str) // ', y=' // real_str(curve%orbit%y, 5)
        if (curve%orbit%t /= 0) str = trim(str) // ', t=' // real_str(curve%orbit%t, 5)
        if (len_trim(str) == 0) then
          curve%legend_text = curve%legend_text(:ix-1) // '[x=y=t=0]'
        else
          curve%legend_text = curve%legend_text(:ix-1) // '[' // trim(str(3:)) // ']'
        endif
      endif

      !! if (all(.not. good)) exit
      n_dat = count(good)
      curve%x_line(1:n_dat) = pack(curve%x_line, mask = good)
      curve%y_line(1:n_dat) = pack(curve%y_line, mask = good)
      call re_allocate (curve%y_line, n_dat) ! allocate space for the data
      call re_allocate (curve%x_line, n_dat) ! allocate space for the data

    ! For non-smooth curves: Draw straight lines through the symbols if
    ! the data uses "ref" or "meas" values. Else evaluate at the element ends.

    else if (straight_line_between_syms) then
      if (allocated (curve%x_symb)) then
        n_dat = size(curve%x_symb)
        call re_allocate (curve%y_line, n_dat) 
        call re_allocate (curve%x_line, n_dat) 
        curve%x_line = curve%x_symb 
        curve%y_line = curve%y_symb 
      else
        call tao_curve_datum_calc (scratch%eles, plot, curve, 'LINE')
        if (.not. curve%valid) return
      endif

    ! Evaluate at element ends

    else

      eps = 1e-4 * (graph%x%max - graph%x%min)             ! a small number
      l_tot = branch%param%total_length
      branch%ele%logic = .false.
      do i = 0, branch%n_ele_track
        ele => branch%ele(i)
        if (graph%x%min == graph%x%max) cycle
        s1 = ele%s
        ele%logic = (s1 >= graph%x%min-eps) .and. (s1 <= graph%x%max+eps)
        if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
          ele%logic = ele%logic .or. ((s1-l_tot >= graph%x%min-eps) .and. (s1-l_tot <= graph%x%max+eps))
        endif
      enddo
      n_dat = count (branch%ele(:)%logic)
      call re_allocate_eles (scratch%eles, n_dat, exact = .true.)
      i = 0
      do j = 0, ubound(branch%ele, 1)
        if (.not. branch%ele(j)%logic) cycle
        i = i + 1
        scratch%eles(i)%ele => branch%ele(j)
      enddo

      ! If there is a wrap-around then reorder the data

      ix1 = 0
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
        do i = 1, n_dat
          if (ix1 == 0 .and. branch%ele(scratch%eles(i)%ele%ix_ele)%s - l_tot > graph%x%min) ix1 = i
          if (branch%ele(scratch%eles(i)%ele%ix_ele)%s < graph%x%max+eps) ix2 = i
        enddo
        if (ix1 /= 0) then
          call re_allocate_eles(scratch%eles, n_dat + ix2 + 1 - ix1, .true., .true.)
          scratch%eles = [scratch%eles(ix1:n_dat), scratch%eles(1:ix2)]
        endif
      endif

      call tao_curve_datum_calc (scratch%eles, plot, curve, 'LINE')
      if (.not. curve%valid) return

      if (substr(curve%data_type, 1, 15) == 'element_attrib.') then
        n = size(curve%x_line)
        call re_allocate(curve%x_line, 4*n)
        call re_allocate(curve%y_line, 4*n)
        do i = n, 1, -1
          ele => branch%ele(scratch%eles(i)%ele%ix_ele)
          j = 4*i - 4

          if (ele%slave_status == super_slave$) then
            ele => pointer_to_lord(ele, 1, ix_slave_back = ix_slave)
            if (ix_slave /= ele%n_slave) then ! Will be duplicate
              curve%x_line(j+1:j+4) = real_garbage$
              cycle
            endif
          endif

          curve%y_line(j+2:j+3) = curve%y_line(i)
          curve%y_line(j+1) = 0
          curve%y_line(j+4) = 0
          curve%x_line(j+1:j+2) = ele%s_start
          curve%x_line(j+3:j+4) = ele%s
          if (ix1 /= 0 .and. i <= n_dat - ix1 + 1) curve%x_line(j+1:j+4) = curve%x_line(j+1:j+4) - l_tot
        enddo

        k = 0
        do j = 0, 4*(n-1), 4
          if (curve%x_line(j+1) == real_garbage$) cycle
          curve%x_line(k+1:k+4) = curve%x_line(j+1:j+4)
          curve%y_line(k+1:k+4) = curve%y_line(j+1:j+4)
          k = k + 4
        enddo
        call re_allocate(curve%x_line, k)
        call re_allocate(curve%y_line, k)

      else
        do i = 1, size(curve%x_line)
          curve%x_line(i) = branch%ele(scratch%eles(i)%ele%ix_ele)%s
          if (ix1 /= 0 .and. i <= n_dat - ix1 + 1) curve%x_line(i) = curve%x_line(i) - l_tot
        enddo
      endif

    endif
  end select
endif

!----------------------------------------------------------------------------
! Note: Since there is an arbitrary overall phase, the phase data 
! gets renormalized so that the average value is zero.

if (index(curve%data_type, 'phase') /= 0 .and. n_dat /= 0 .and. zero_average_phase) then
  if (allocated(curve%y_symb)) then
    f = sum(curve%y_symb) / n_dat
    curve%y_symb = curve%y_symb - f
    if (allocated(curve%y_line)) curve%y_line = curve%y_line - f 
  elseif (allocated(curve%y_line)) then
    f = sum(curve%y_line) / n_dat
    curve%y_line = curve%y_line - f
  endif
endif 

err_flag = .false.

end subroutine tao_curve_data_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_calc_data_at_s_pts (tao_lat, curve, comp_sign, good)

use transfer_map_mod
use twiss_and_track_mod, only: twiss_and_track_at_s

implicit none

type (tao_lattice_struct), target :: tao_lat
type (tao_curve_struct) curve
type (tao_lattice_branch_struct), pointer :: tao_branch
type (bunch_track_struct), pointer :: bunch_track
type (bunch_params_struct), pointer :: bunch_params0, bunch_params1
type (bunch_params_struct) :: bunch_params
type (coord_struct), pointer :: orb(:), orb_ref
type (coord_struct) orbit_end, orbit_last, orbit
type (lat_struct), pointer :: lat
type (ele_struct), target :: ele_to_s
type (ele_struct), pointer :: ele_last, ele_here, ele_ref, this_ele, ele0
type (taylor_struct) t_map(6)
type (branch_struct), pointer :: branch
type (all_pointer_struct) a_ptr

real(rp) x1, x2, cbar(2,2), s_last, s_now, value, mat6(6,6), vec0(6), mat0(6,6)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), one_pz, gamma, len_tot
real(rp) comp_sign, vec3(3), r_bunch, ds, dt, time
real(rp), allocatable :: val_arr(:)

integer i, ii, ix, j, k, np, expnt(6), ix_ele_here, ix_ref, ix_branch, idum, n_ele_track
integer cache_status
integer, parameter :: loading_cache$ = 1, using_cache$ = 2

character(40) name, sub_data_type, data_type_select, data_source
character(100) why_invalid
character(200) data_type
character(*), parameter :: r_name = 'tao_calc_data_at_s_pts'
logical err_flag, good(:), first_time, radiation_fluctuations_on, ok

! Some init

data_type = curve%data_type

ix_branch = tao_branch_index(curve%ix_branch)
lat => tao_lat%lat
tao_branch => tao_lat%tao_branch(ix_branch)
orb => tao_branch%orbit
branch => lat%branch(ix_branch)
first_time = .true.
n_ele_track = branch%n_ele_track

ele_ref => tao_curve_ele_ref(curve, .true.)
if (associated(ele_ref)) then
  ix_ref = ele_ref%ix_ele
else
  ix_ref = 0   ! If beginning element is needed as ref.
endif

if (lat%param%geometry == closed$ .and. .not. lat%param%stable) then
  curve%g%why_invalid = 'Unstable Lattice'
  good = .false.
  return
endif

if (curve%data_source == 'lat') then
  select case (data_type(1:5))
  case ('emitt', 'norm_')
    curve%g%why_invalid = 'curve%data_source = "lat" is not compatable with data_type: ' // data_type
    call out_io (s_warn$, r_name, curve%g%why_invalid)
    call out_io (s_blank$, r_name, "Will not perform any plot smoothing")
    good = .false.
    return
  end select 
endif

! x1 and x2 are the longitudinal end points of the plot

call tao_graph_s_min_max_calc(curve%g, branch, x1, x2)
if (.not. curve%g%is_valid) return

radiation_fluctuations_on = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

ele_ref => branch%ele(ix_ref)
orb_ref => orb(ix_ref)
s_last = ele_ref%s
orbit_last = orbit

ix = index(data_type, '.')
if (data_type(1:11) == 'expression:') then
  data_type_select = 'expression'
elseif (ix == 0) then
  data_type_select = data_type
else
  data_type_select = data_type(:ix-1)
  sub_data_type = data_type(ix+1:)
endif

! Only cache plot data if the number of points is equal to s%plot_page%n_curve_pts

if (curve%data_source == 'lat') then
  np = size(curve%x_line)
  if (tao_branch%plot_cache_valid .and. tao_branch%cache_x_min == x1 .and. &
            tao_branch%cache_x_max == x2 .and. tao_branch%cache_n_pts == np) then
    cache_status = using_cache$

  else
    cache_status = loading_cache$
    if (allocated(tao_branch%plot_cache)) then
      if (size(tao_branch%plot_cache) /= np) call tao_deallocate_plot_cache(tao_branch%plot_cache)
    endif
    if (.not. allocated(tao_branch%plot_cache)) allocate (tao_branch%plot_cache(np))
    tao_branch%cache_x_min = x1
    tao_branch%cache_x_max = x2 
    tao_branch%cache_n_pts = size(curve%x_line)
    tao_branch%plot_cache_valid = .true.
  endif
endif

!

select case (data_type_select)
case ('momentum_compaction', 'r56_compaction', 'r')
  call mat6_from_s_to_s (lat, mat0, vec0, branch%ele(0)%s, ele_ref%s, orb(0), ix_branch = ix_branch)
  call mat_inverse(mat0, mat0)
end select

!

ele_last => branch%ele(0)

do ii = 1, size(curve%x_line)

  ! Good(ii) may be false if this is not first time tao_calc_data_at_s_pts is called from tao_curve_data_setup.
  ! For example, tao_calc_data_at_s_pts is called twice when plotting "meas - design".

  if (.not. good(ii)) then
    first_time = .true.
    cycle
  endif

  s_now = x1 + (ii-1) * (x2-x1) / (size(curve%x_line)-1)
  if (associated(tao_hook_curve_s_pt_ptr))s_now = &
                        tao_hook_curve_s_pt_ptr (s_now, ii, x1, x2, size(curve%x_line), tao_lat, curve)

  if (s_now > branch%ele(n_ele_track)%s) s_now = branch%ele(n_ele_track)%s
  value = 0

  ! Check if in a hybrid or taylor element within which interpolation cannot be done.

  ix_ele_here = element_at_s (lat, s_now, .true., ix_branch, err_flag)
  ele_here => branch%ele(ix_ele_here)

  if (ele_here%key == hybrid$ .or. ele_here%key == taylor$ .or. err_flag .or. &
      ele_last%key == hybrid$ .or. ele_last%key == taylor$) then
    if (err_flag .or. s_last == ele_here%s) then
      good(ii) = .false.
      first_time = .true.
      cycle
    elseif (s_now < (ele_here%s_start + ele_here%s)/2 .or. s_last >= ele_here%s_start) then
      s_now = ele_here%s
    else
      s_now = ele_here%s_start
    endif

    if (s_now == s_last) then
      good(ii) = .false.
      first_time = .true.
      cycle
    endif
  endif

  curve%x_line(ii) = s_now
  ele_last => ele_here

  if (ii > 1 .and. s_now <= s_last) then
    good(ii) = .false.
    first_time = .true.
    cycle
  endif

  !-----------------------------

  select case (curve%data_source)
  case ('beam')
    if (.not. allocated(tao_branch%bunch_params)) then
      call out_io (s_error$, r_name, 'BUNCH_PARAMS NOT ALLOCATED.')
      return
    endif
 
    if (allocated(tao_branch%bunch_params_comb)) then
      ix = max(1, curve%ix_bunch)
      bunch_track => tao_branch%bunch_params_comb(ix)
      np = bunch_track%n_pt
      ix = bracket_index (s_now, bunch_track%pt(0:np)%s, 0)
      bunch_params0 => bunch_track%pt(ix)
      bunch_params1 => bunch_track%pt(min(ix,np))
    else
      np = n_ele_track
      ix = bracket_index (s_now, tao_branch%bunch_params(0:np)%s, 0)
      bunch_params0 => tao_branch%bunch_params(ix)
      bunch_params1 => tao_branch%bunch_params(min(ix,np))
    endif


    if (bunch_params0%s == bunch_params1%s) then
      r_bunch = 0
    else
      r_bunch = (s_now - bunch_params0%s) / (bunch_params1%s - bunch_params0%s)
    endif

    ! Cannot return if no particles since bunch may be injected not at the start of the branch.
    if (bunch_params0%n_particle_live == 0) then
      good(ii) = .false.
      cycle
    endif

    orbit%vec = (1-r_bunch) * bunch_params0%centroid%vec + r_bunch * bunch_params1%centroid%vec

    bunch_params%sigma = (1-r_bunch) * bunch_params0%sigma + (1-r_bunch) * bunch_params1%sigma

    bunch_params%a%emit = (1-r_bunch) * bunch_params0%a%emit + (1-r_bunch) * bunch_params1%a%emit
    bunch_params%b%emit = (1-r_bunch) * bunch_params0%b%emit + (1-r_bunch) * bunch_params1%b%emit

    bunch_params%x%emit = (1-r_bunch) * bunch_params0%x%emit + (1-r_bunch) * bunch_params1%x%emit
    bunch_params%y%emit = (1-r_bunch) * bunch_params0%y%emit + (1-r_bunch) * bunch_params1%y%emit

    bunch_params%a%norm_emit = (1-r_bunch) * bunch_params0%a%norm_emit + (1-r_bunch) * bunch_params1%a%norm_emit
    bunch_params%b%norm_emit = (1-r_bunch) * bunch_params0%b%norm_emit + (1-r_bunch) * bunch_params1%b%norm_emit
    bunch_params%z%norm_emit = (1-r_bunch) * bunch_params0%z%norm_emit + (1-r_bunch) * bunch_params1%z%norm_emit

  case ('lat')
    if (cache_status == using_cache$) then
      ele_to_s  = tao_branch%plot_cache(ii)%ele_to_s
      orbit     = tao_branch%plot_cache(ii)%orbit
      mat6      = tao_branch%plot_cache(ii)%ele_to_s%mat6
      vec0      = tao_branch%plot_cache(ii)%ele_to_s%vec0
      err_flag  = tao_branch%plot_cache(ii)%err

    else
      ! Note: first_time may be set True when a Taylor or Hybrid element is encountered.
      if (first_time) then
        call twiss_and_track_at_s (lat, s_now, ele_to_s, orb, orbit, ix_branch, err_flag, compute_floor_coords = .true.)
        call mat6_from_s_to_s (lat, mat6, vec0, branch%ele(0)%s, s_now, orb(0), ix_branch = ix_branch)
        ele_to_s%vec0 = vec0
        ele_to_s%mat6 = mat6
        orbit_end = orbit
        first_time = .false.

      else
        call twiss_and_track_from_s_to_s (branch, orbit, s_now, orbit_end, ele_to_s, ele_to_s, err_flag, &
                                                               compute_floor_coords = .true., compute_twiss = .false.)
        orbit = orbit_end
        vec0 = matmul(ele_to_s%mat6, vec0) + ele_to_s%vec0
        mat6 = matmul(ele_to_s%mat6, mat6)  ! Matrix from beginning of branch.
        ele_to_s%vec0 = vec0
        ele_to_s%mat6 = mat6
        ele_to_s%key = hybrid$                  ! So twiss_propagate1 does not get confused.
        call twiss_propagate1(branch%ele(0), ele_to_s, err_flag)
        ! Correct phase that can be off by factors of twopi
        ele0 => branch%ele(ele_here%ix_ele-1)
        ele_to_s%a%phi = ele_to_s%a%phi + twopi * nint((0.5_rp *(ele0%a%phi + ele_here%a%phi) - ele_to_s%a%phi) / twopi)
        ele_to_s%b%phi = ele_to_s%b%phi + twopi * nint((0.5_rp *(ele0%b%phi + ele_here%b%phi) - ele_to_s%b%phi) / twopi)
      endif


      if (cache_status == loading_cache$) then
        tao_branch%plot_cache(ii)%ele_to_s = ele_to_s
        tao_branch%plot_cache(ii)%orbit    = orbit
        tao_branch%plot_cache(ii)%err      = err_flag
      endif

      if (err_flag) then
        tao_branch%plot_cache(ii:)%err = .true.
        good(ii:) = .false.
        bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
        return
      endif
    endif

    if (orbit%state /= alive$) then
      good(ii:) = .false.
      tao_branch%plot_cache(ii:)%err = .true.
      write (curve%message_text, '(f10.3)') s_now
      curve%message_text = trim(curve%data_type) // ': Particle lost at s = ' // &
                           trim(adjustl(curve%message_text))
      bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
      return
    endif

  case default
    call out_io (s_fatal$, r_name, 'I DO NOT KNOW HOW TO HANDLE THIS curve%data_source: ' // curve%data_source)
    return
  end select

  call this_value_at_s (data_type_select, sub_data_type, value, ii, s_last, s_now, &
                          tao_branch, orbit, lat, branch, ele_to_s, ele_ref, ele_here, mat6, err_flag);  if (err_flag) return


  curve%y_line(ii) = curve%y_line(ii) + comp_sign * value
  s_last = s_now
enddo

! Subtract reference?

if (curve%ele_ref_name /= '') then
  select case (data_type_select)
  case ('r', 't', 'tt')
  case default
    if (curve%data_source /= 'lat') then
      call tao_set_curve_invalid(curve, 'SETTING ele_ref_name FOR CURVE WHERE curve%data_source IS NOT "lat" IS NOT VALID.')
      good = .false.
      bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
      return
    endif

    ele_ref => tao_curve_ele_ref(curve, .false.)
    s_now = ele_ref%s
    call twiss_and_track_at_s (lat, s_now, ele_to_s, orb, orbit, ix_branch, err_flag, compute_floor_coords = .true.)

    if (err_flag) then
      tao_branch%plot_cache(ii:)%err = .true.
      good(ii:) = .false.
      bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
      return
    endif

    if (orbit%state /= alive$) then
      good(ii:) = .false.
      tao_branch%plot_cache(ii:)%err = .true.
      write (curve%message_text, '(f10.3)') s_now
      curve%message_text = trim(curve%data_type) // ': Particle lost at s = ' // &
                           trim(adjustl(curve%message_text))
      bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
      return
    endif

    call this_value_at_s (data_type_select, sub_data_type, value, ii, s_last, s_now, &
                  tao_branch, orbit, lat, branch, ele_to_s, ele_ref, ele_here, mat6, err_flag);  if (err_flag) return

    curve%y_line = curve%y_line - comp_sign * value
  end select
endif
!

bmad_com%radiation_fluctuations_on = radiation_fluctuations_on

!--------------------------------------------------------
contains

subroutine this_value_at_s (data_type_select, sub_data_type, value, ii, s_last, s_now, &
                                   tao_branch, orbit, lat, branch, ele_to_s, ele_ref, ele_here, mat6, err_flag)

type (coord_struct) orbit, orb_end
type (tao_lattice_branch_struct) :: tao_branch
type (ele_struct), target :: ele_to_s, ele_ref, ele_dum, high_ele, low_ele, ele_here
type (lat_struct) lat
type (lat_struct), pointer :: this_lat
type (branch_struct) branch
type (branch_struct), pointer :: this_branch
type (twiss_struct), pointer :: z0, z1, z2

real(rp) value, s_last, s_now, ds, m6(6,6), dalpha, dbeta, aa, bb, dE, mat6(6,6)
integer status, ii, i, j
logical is_ok, err_flag
character(*) data_type_select, sub_data_type
character(40) name

!

select case (data_type_select)

case ('apparent_emit', 'norm_apparent_emit')
  select case (data_type)
  case ('apparent_emit.x', 'norm_apparent_emit.x')
    if (curve%data_source == 'beam') then
      value = tao_beam_emit_calc (x_plane$, apparent_emit$, ele_to_s, bunch_params)
    else
      value = tao_lat_emit_calc (x_plane$, apparent_emit$, ele_to_s, tao_branch%modes_6d)
    endif
    if (data_type_select(1:4) == 'norm') value = value * ele_to_s%value(E_tot$) / mass_of(branch%param%particle)
  case ('apparent_emit.y', 'norm_apparent_emit.y')
    if (curve%data_source == 'beam') then
      value = tao_beam_emit_calc (y_plane$, apparent_emit$, ele_to_s, bunch_params)
    else
      value = tao_lat_emit_calc (y_plane$, apparent_emit$, ele_to_s, tao_branch%modes_6d)
    endif
    if (data_type_select(1:4) == 'norm') value = value * ele_to_s%value(E_tot$) / mass_of(branch%param%particle)
  case default
    goto 9000  ! Error message & Return
  end select

case ('chrom')
  dE = 2 * s%global%delta_e_chrom  ! Actually this is the change in pz
  ds = s_now - ele_here%s_start
  i = ele_here%ix_ele

  this_lat => tao_lat%high_E_lat
  this_branch => this_lat%branch(ix_branch)
  call twiss_and_track_intra_ele (this_branch%ele(i), this_branch%param, 0.0_rp, ds, .true., .true., &
              this_branch%ele(i)%map_ref_orb_in, orb_end, this_branch%ele(i-1), high_ele)
  
  this_lat => tao_lat%low_E_lat
  this_branch => this_lat%branch(ix_branch)
  call twiss_and_track_intra_ele (this_branch%ele(i), this_branch%param, 0.0_rp, ds, .true., .true., &
              this_branch%ele(i)%map_ref_orb_in, orb_end, this_branch%ele(i-1), low_ele)

  if (data_type == 'chrom.w.a') then
    z2 => high_ele%a
    z1 => low_ele%a
    Z0 => ele_to_s%a
  else
    z2 => high_ele%b
    z1 => low_ele%b
    z0 => ele_to_s%b
  endif

  dalpha = (z2%alpha - z1%alpha) / dE
  dbeta  = (z2%beta - z1%beta) / dE
  aa = dalpha - z0%alpha * dbeta / z0%beta
  bb = dbeta / z0%beta
  value = sqrt(aa**2 + bb**2)

case ('element_attrib')
  name = upcase(curve%data_source(16:))
  ele_dum%key = overlay$  ! so entire attribute name table will be searched
  i = attribute_index(ele_dum, name)
  if (i < 1) goto 9000
  call pointer_to_attribute (ele_ref, name, .false., a_ptr, err_flag, .false.)
  if (associated (a_ptr%r)) value = a_ptr%r

case ('emit')
  select case (data_type)
  case ('emit.a')
    value = bunch_params%a%emit
  case ('emit.b')
    value = bunch_params%b%emit
  case ('emit.x', 'norm_emit.x')
    if (curve%data_source == 'beam') then
      value = bunch_params%x%emit
    else
      value = tao_lat_emit_calc (x_plane$, projected_emit$, ele_to_s, tao_branch%modes_6d)
    endif
    if (data_type_select(1:4) == 'norm') value = value * ele_to_s%value(E_tot$) / mass_of(branch%param%particle)
  case ('emit.y', 'norm_emit.y')
    if (curve%data_source == 'beam') then
      value = bunch_params%y%emit
    else
      value = tao_lat_emit_calc (y_plane$, projected_emit$, ele_to_s, tao_branch%modes_6d)
    endif
    if (data_type_select(1:4) == 'norm') value = value * ele_to_s%value(E_tot$) / mass_of(branch%param%particle)
  case default
    goto 9000  ! Error message & Return
  end select

case ('expression')
  this_ele => ele_to_s  ! Need a pointer to an ele_to_s
  call tao_evaluate_expression (data_type(12:), 1, .false., val_arr, err_flag, .true., &
            dflt_source = 'at_ele', dflt_ele = this_ele, dflt_uni = tao_lat%u%ix_uni, &
            dflt_orbit = orbit);  if (err_flag) goto 9000  ! Error message & Return
  value = val_arr(1)

case ('momentum_compaction')
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb_ref%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb_ref%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb_ref%vec(4) / one_pz
  ds = ele_to_s%s - branch%ele(0)%s
  if (ds == 0) then
    value = 0
  else
    m6 = matmul(mat6, mat0)
    value = -(sum(m6(5,1:4) * eta_vec) + m6(5,6)) / ds
  endif

case ('norm_emit')
  select case (data_type)
  case ('norm_emit.a')
    value = bunch_params%a%norm_emit
  case ('norm_emit.b')
    value = bunch_params%b%norm_emit
  case ('norm_emit.z')
    value = bunch_params%z%norm_emit
  case default
    goto 9000  ! Error message & Return
  end select

case ('r')
  i = tao_read_phase_space_index (data_type, 3); if (i == 0) goto 9000
  j = tao_read_phase_space_index (data_type, 4); if (j == 0) goto 9000
  m6 = matmul(mat6, mat0)
  value = m6(i, j)

case ('r56_compaction')
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb_ref%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb_ref%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb_ref%vec(4) / one_pz
  m6 = matmul(mat6, mat0)
  value = sum(m6(5,1:4) * eta_vec) + m6(5,6)

case ('sigma')
  if (curve%data_source == 'beam') then
    select case (data_type)
    case ('sigma.x');    value = sqrt(bunch_params%sigma(1,1))
    case ('sigma.px');   value = sqrt(bunch_params%sigma(2,2))
    case ('sigma.y');    value = sqrt(bunch_params%sigma(3,3))
    case ('sigma.py');   value = sqrt(bunch_params%sigma(4,4))
    case ('sigma.z');    value = sqrt(bunch_params%sigma(5,5))
    case ('sigma.pz');   value = sqrt(bunch_params%sigma(6,6))
    case default;        goto 9000  ! Error message & Return
    end select
  else
    m6 = matmul(matmul(mat6, tao_branch%lat_sigma(0)%mat), transpose(mat6))
    select case (data_type)
    case ('sigma.x');    value = sqrt(m6(1,1))
    case ('sigma.px');   value = sqrt(m6(2,2))
    case ('sigma.y');    value = sqrt(m6(3,3))
    case ('sigma.py');   value = sqrt(m6(4,4))
    case ('sigma.z');    value = sqrt(m6(5,5))
    case ('sigma.pz');   value = sqrt(m6(6,6))
    case default;        goto 9000  ! Error message & Return
    end select
  endif

case ('t', 'tt')
  if (ii == 1) then
    call twiss_and_track_at_s (lat, s_last, orb = orb, orb_at_s = orbit, ix_branch = ix_branch)
    call taylor_make_unit (t_map, orbit%vec)
  endif
  if (s_now < s_last) return
  expnt = 0
  i = tao_read_phase_space_index (sub_data_type, 1); if (i == 0) goto 9000
  do j = 2, 20
    if (sub_data_type(j:j) == '') exit
    k = tao_read_phase_space_index (sub_data_type, j); if (k == 0) goto 9000
    expnt(k) = expnt(k) + 1
  enddo
  call transfer_map_from_s_to_s (lat, t_map, s_last, s_now, ix_branch = ix_branch, &
                          unit_start = .false., err_flag = err_flag, concat_if_possible = s%global%concatenate_maps)
  if (err_flag) return
  value = taylor_coef (t_map(i), expnt)

case default
  value = tao_param_value_at_s (data_type, ele_to_s, ele_here, orbit, err_flag, why_invalid)
  if (err_flag) then
    call tao_set_curve_invalid(curve, why_invalid, .true.)
    goto 9100  ! Cleanup & Return
  endif
end select

return

! Error message

9000 continue
call tao_set_curve_invalid(curve, 'CANNOT EVALUATE FOR SMOOTH LINE CALC: ' // data_type, .true.)
9100 continue
good = .false.
bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
if (cache_status == loading_cache$) tao_branch%plot_cache_valid = .false.

end subroutine this_value_at_s

end subroutine tao_calc_data_at_s_pts

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_useit_plot_calc (curve, graph, data, check_s_position, most_invalid)
!
! Routine to set the data for plotting.
!
! Input:
!   graph             -- tao_graph_struct
!   curve             -- tao_curve_struct
!   check_s_position  -- logical: If present and True then 
!                           veto data that does not have an s-position.
! Output:
!   data(:)           -- Tao_data_struct:
!     %useit_plot         -- True if good for plotting.
!   most_invalid      -- character(*): String documenting biggest invalid data problem. 
!-

subroutine tao_data_useit_plot_calc (curve, graph, data, check_s_position, most_invalid)

implicit none

type (tao_curve_struct) curve
type (tao_graph_struct) graph
type (tao_data_struct) data(:)
integer id, n_bad_exists, n_bad_plot, n_bad_user, n_bad_meas, n_bad_ref, n_bad_model, n_bad_s
logical check_s_position
character(*) most_invalid

!

n_bad_exists = 0;  n_bad_plot = 0;  n_bad_user = 0;  n_bad_meas = 0;  
n_bad_ref = 0;  n_bad_model = 0;  n_bad_s = 0

do id = 1, size(data)
  data(id)%useit_plot = data(id)%exists
  if (is_bad(data(id)%useit_plot, n_bad_exists)) cycle

  data(id)%useit_plot = data(id)%useit_plot .and. data(id)%good_plot
  if (is_bad(data(id)%useit_plot, n_bad_plot)) cycle

  if (graph%draw_only_good_user_data_or_vars) then
    data(id)%useit_plot = data(id)%useit_plot .and. data(id)%good_user
    if (is_bad(data(id)%useit_plot, n_bad_user)) cycle
  endif

  if (index(curve%component, 'meas') /= 0) then
    data(id)%useit_plot = data(id)%useit_plot .and. data(id)%good_meas
    if (is_bad(data(id)%useit_plot, n_bad_user)) cycle
  endif

  if (index(curve%component, 'ref') /= 0) then
    data(id)%useit_plot = data(id)%useit_plot .and. data(id)%good_ref
    if (is_bad(data(id)%useit_plot, n_bad_user)) cycle
  endif

  if (index(curve%component, 'model') /= 0) then
    data(id)%useit_plot = data(id)%useit_plot .and. data(id)%good_model
    if (is_bad(data(id)%useit_plot, n_bad_user)) cycle
  endif

  if (check_s_position) then
    data(id)%useit_plot = data(id)%useit_plot .and. (data(id)%s /= real_garbage$)
    if (is_bad(data(id)%useit_plot, n_bad_user)) cycle
  endif
enddo

if (max(n_bad_exists,n_bad_plot, n_bad_user, n_bad_meas, n_bad_ref, n_bad_model, n_bad_s) < size(data)/2) then
  most_invalid = ''
elseif (n_bad_exists >= max(n_bad_plot, n_bad_user, n_bad_meas, n_bad_ref, n_bad_model, n_bad_s)) then
  most_invalid = 'Biggest problem: Data does not exist.'
elseif (n_bad_plot >= max(n_bad_user, n_bad_meas, n_bad_ref, n_bad_model, n_bad_s)) then
  most_invalid = 'Biggest problem: Data ouside of plot range.'
elseif (n_bad_user >= max(n_bad_meas, n_bad_ref, n_bad_model, n_bad_s)) then
  most_invalid = 'Biggest problem: Data vetoed by user.'
elseif (n_bad_meas >= max(n_bad_ref, n_bad_model, n_bad_s)) then
  most_invalid = 'Biggest problem: Lack of valid meas data values'
elseif (n_bad_ref >= max(n_bad_model, n_bad_s)) then
  most_invalid = 'Biggest problem: Lack of valid ref data values.'
elseif (n_bad_model >= n_bad_s) then
  most_invalid = 'Biggest problem: Lack of valid model data values.'
else
  most_invalid = 'Data not having a well defined s-position. (plot by datum index instead?)'
endif
!-------------------------------------
contains

function is_bad(useit_plot, n_bad) result (bad)
integer n_bad
logical useit_plot, bad
!
bad = (.not. useit_plot)
if (bad) n_bad = n_bad + 1
end function

end subroutine tao_data_useit_plot_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_curve_datum_calc (eles, plot, curve, who)
!
! Routine to calculate datum values. 
! The values are calculated at the end of each eles(:)%ele element.
!
! Input:
!   eles(:)   -- ele_pointer_struct: Array of elements.
!   plot      -- Tao_plot_struct:
!   curve     -- Tao_curve_struct:
!   who       -- Character(*): Where to put the data. 
!                  Either: "SYMBOL" or "LINE".
!
! Output:
!   curve -- Tao_curve_struct: Structure holding the datum values
!-

subroutine tao_curve_datum_calc (eles, plot, curve, who)

implicit none

type (tao_plot_struct) plot
type (tao_curve_struct) curve
type (tao_universe_struct), pointer :: u
type (tao_data_struct) datum
type (taylor_struct) t_map(6)
type (ele_pointer_struct), allocatable, target :: eles(:)
type (ele_struct), pointer :: ele, ele_ref

real(rp) y_val

integer i, j, m, ie, n_dat

logical err_flag, valid
logical, allocatable :: good(:)

character(*) who
character(20) :: r_name = 'tao_curve_datum_calc'
character(80) why_invalid

! calculate the y-axis data point values.

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
n_dat = size(eles)

call re_allocate (good, n_dat)   ! allocate space for the data
call re_allocate (scratch%y_value, n_dat) ! allocate space for the data

scratch%y_value = 0
good = .true.
ele_ref => tao_curve_ele_ref(curve, .true.)

datum%exists         = .true.
if (associated(ele_ref)) then
  datum%ix_ele_ref     = ele_ref%ix_ele
else
  datum%ix_ele_ref     = -1
endif
datum%ix_ele_start   = -1
datum%ele_start_name = ''
datum%merit_type     = 'target'
datum%data_type      = curve%data_type
datum%ele_ref_name   = curve%ele_ref_name
datum%data_source    = curve%data_source
datum%ix_branch      = tao_branch_index(curve%ix_branch)

call tao_split_component (curve%component, scratch%comp, err_flag)
if (err_flag) then
  call tao_set_curve_invalid (curve, 'BAD CURVE COMPONENT EXPRESSION: ' // curve%component)
  return
endif
if (curve%component == '') then
  call tao_set_curve_invalid (curve, 'BLANK CURVE COMPONENT STRING.')
  return
endif

do m = 1, size(scratch%comp)

  do ie = 1, n_dat

    datum%ix_ele = eles(ie)%ele%ix_ele
    datum%ix_branch = eles(ie)%ele%ix_branch
    datum%exists = .true.

    select case (scratch%comp(m)%name)
    case ('') 
      cycle
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, y_val, valid, why_invalid)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, y_val, valid, why_invalid)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, y_val, valid, why_invalid)
    case ('ref', 'meas')
      call tao_set_curve_invalid(curve, 'CURVE COMPONENT: ' // trim(scratch%comp(m)%name) // &
                                                       ' NOT ALLOWED WITH DATA_SOURCE IS: ' // curve%data_source)
      return
    case default
      call tao_set_curve_invalid(curve, 'BAD PLOT COMPONENT: ' // scratch%comp(m)%name)
      return
    end select
    scratch%y_value(ie) = scratch%y_value(ie) + scratch%comp(m)%sign * y_val
    if (.not. valid) good(ie) = .false.
    if (substr(datum%data_type,1,3) == 'tt.' .or. substr(datum%data_type,1,2) == 't.') then
      if (datum%ix_ele < datum%ix_ele_ref) datum%ix_ele_ref = datum%ix_ele
    endif

  enddo
enddo

if (n_dat > 0 .and. all(.not. good)) then
  call tao_set_curve_invalid(curve, why_invalid)
  return
endif

n_dat = count(good)

if (who == 'SYMBOL') then
  call re_allocate (curve%x_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%ix_symb, n_dat)
  j = 0
  do i = 1, size(eles)
    if (.not. good(i)) cycle
    j = j + 1
    if (plot%x_axis_type == 's') then
      if (substr(datum%data_type, 1, 15) == 'element_attrib.') then
        ! Element attributes reference ele center.
        ele => eles(i)%ele
        if (ele%slave_status == super_slave$) ele => pointer_to_lord(ele, 1)
        curve%x_symb(j)  = 0.5_rp * (ele%s_start + ele%s)  
      else        
        curve%x_symb(j)  = eles(i)%ele%s
      endif
    else  ! 'index' or 'ele_index'
      curve%x_symb(j)  = eles(i)%ele%ix_ele
    endif
    curve%ix_symb(j) = eles(i)%ele%ix_ele
    curve%y_symb(j)  = scratch%y_value(i)
  enddo

else
  call re_allocate (curve%x_line, n_dat) ! allocate space for the data
  call re_allocate (curve%y_line, n_dat) ! allocate space for the data
  j = 0
  do i = 1, size(eles)
    if (.not. good(i)) cycle
    j = j + 1
    curve%x_line(j)  = eles(i)%ele%ix_ele
    curve%y_line(j)  = scratch%y_value(i)
    eles(j) = eles(i)
  enddo
  call re_allocate_eles (eles, n_dat, .true., .true.)
endif

valid = .true.

end subroutine tao_curve_datum_calc 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_curve_invalid (curve, why_invalid, print_err)
!
! Routine to set curve%valid to False.
!
! Input:
!   curve       -- tao_curve_struct: Curve to set.
!   why_invalid -- character(*): Invalid information.
!   print_err   -- logical, optional: If present and True then also print an error message.
!
! Output:
!   curve       -- tao_curve_struct: Curve properly set.
!-

subroutine tao_set_curve_invalid (curve, why_invalid, print_err)

type (tao_curve_struct) curve
character(*) why_invalid
logical, optional :: print_err
character(*), parameter :: r_name = 'tao_set_curve_invalid'

!

curve%valid = .false.
curve%why_invalid = trim(curve%name) // ': ' // why_invalid
if (logic_option(.false., print_err)) call out_io (s_error$, r_name, why_invalid, 'FOR CURVE: ' // tao_curve_name(curve))

end subroutine tao_set_curve_invalid

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function tao_curve_check_universe (curve, uni) result (is_ok)
!
! Routine to check if the universe associated with a curve exists and is on.
!
! Input:
!   curve       -- tao_curve_struct: Curve to check.
!   uni         -- tao_universe_struct, pointer: Associated universe
!
! Output:
!   curve       -- tao_curve_struct: Curve%valid set to False if needed.
!   is_ok       -- logical: Set True if associated universe exists and is on.
!-

function tao_curve_check_universe (curve, uni) result (is_ok)

type (tao_curve_struct) curve
type (tao_universe_struct), pointer :: uni
logical is_ok

!

if (.not. associated(uni)) then
  call tao_set_curve_invalid (curve, 'NO ASSOCIATED UNIVERSE.')

else
  curve%valid = uni%is_on
  if (.not. curve%valid) then
    write (curve%why_invalid, '(a, i0, a)') trim(curve%name) // ': UNIVERSE ', uni%ix_uni, ' IS OFF!'
  endif
endif

is_ok = curve%valid

end function tao_curve_check_universe

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_graph_s_min_max_calc(graph, branch, s_min, s_max)
!
! Routine to calculate min and max for a graph when plot%x_axis_type is set to "s".
!
! Input:
!   graph   -- tao_graph_struct: Graph to calculate for.
!   branch  -- branch_struct: Associated lattice branch.
!
! Output:
!   s_min   -- real(rp): Graph min. May be negative with graph%allow_wrap_around = T.
!   s_max   -- real(rp): Graph max.
!-

subroutine tao_graph_s_min_max_calc(graph, branch, s_min, s_max)

implicit none

type (tao_graph_struct) graph
type (branch_struct) branch

real(rp) s_min, s_max, len_tot
integer n

!

n = branch%n_ele_track

s_min = branch%ele(0)%s
s_max = branch%ele(n)%s

len_tot = s_max - s_min

if (graph%x%eval_min /= graph%x%eval_max) then
  if (branch%param%geometry == closed$) then
    s_min = min(branch%ele(n)%s, max(graph%x%eval_min, s_min-len_tot))
    s_max = min(s_max, max(graph%x%eval_max, branch%ele(0)%s-len_tot))
  else
    s_min = min(branch%ele(n)%s, max(graph%x%eval_min, s_min))
    s_max = min(s_max, max(graph%x%eval_max, branch%ele(0)%s))
  endif
endif

if (s_min == s_max) then
  graph%is_valid = .false.
  graph%why_invalid = 'Graph S-position range does not overlap lattice branch range.'
endif

end subroutine tao_graph_s_min_max_calc

end module
