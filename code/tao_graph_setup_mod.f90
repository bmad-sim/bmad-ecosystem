module tao_graph_setup_mod

use tao_mod
use tao_lattice_calc_mod
use tao_plot_mod
use tao_command_mod
use tao_evaluate_mod

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

!

graph%valid = .true.   ! assume everything OK
graph%why_invalid = ''
graph%text_legend = ''

call tao_hook_graph_setup (plot, graph, found)
if (found) return

u => tao_pointer_to_universe(graph%ix_universe)
if (.not. u%is_on) then
  graph%valid = .false.
  write (graph%why_invalid, '(a, i0, a)') 'UNIVERSE ', u%ix_uni, ' IS OFF!'
  return
endif

select case (graph%type)
case ('phase_space')
  call tao_phase_space_graph_setup (plot, graph)
case ('data')
  if (plot%x_axis_type == 'data') then
    call tao_data_slice_graph_setup(plot, graph)
  else
    call tao_data_graph_setup(plot, graph)
  endif
end select

! Renormalize

if (allocated (graph%curve)) then
  do i = 1, size(graph%curve)
    curve => graph%curve(i)
    if (allocated(curve%x_symb)) then
        curve%x_symb = curve%x_symb * curve%x_axis_scale_factor
        curve%y_symb = curve%y_symb * curve%y_axis_scale_factor
    endif
    if (allocated(curve%x_line)) then
      curve%x_line = curve%x_line * curve%x_axis_scale_factor
      curve%y_line = curve%y_line * curve%y_axis_scale_factor
    endif
  enddo
endif

call tao_hook_graph_postsetup (plot, graph)

end subroutine tao_graph_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_data_slice_graph_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve

real(rp), allocatable, save :: x(:), y(:)
real(rp), pointer :: symb(:)
real(rp) value

integer i, j, k, m, n_symb, ix

character(160) name
character(40) :: r_name = 'tao_data_slice_graph_setup'

logical err
logical, allocatable, save :: gx(:), gy(:), gix(:)

! setup the graph suffix

graph%valid = .false.
if (size(graph%curve) == 0) return

graph%title_suffix = ''

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  if (index(curve%data_type, '#ref') /= 0) graph%title_suffix = &
                  trim(graph%title_suffix) // '[At: ' // trim(curve%ele_ref_name) // ']'
enddo

if (graph%component /= '') graph%title_suffix = &
              trim(graph%title_suffix) // ' [' // trim(graph%component) // ']'

! loop over all curves

curve_loop: do k = 1, size(graph%curve)

  curve => graph%curve(k)

  ! Find data points

  do i = 1, 2  ! x-axis, y-axis

    if (i == 1) name = curve%data_type_x
    if (i == 2) name = curve%data_type
    call tao_data_type_substitute (name, name, curve%ele_ref_name, graph%component)
    if (i == 1) call tao_to_real_vector (name, 'BOTH', 0, .true., x, gx, err)
    if (i == 2) call tao_to_real_vector (name, 'BOTH', 0, .true., y, gy, err)
    if (err) then
      call out_io (s_error$, r_name, &
                'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // tao_curve_name(curve))   
      return
    endif

  enddo

  ! How many good points?

  if (size(x) /= size(y)) then
    call out_io (s_error$, r_name, &
                  'ARRAY SIZES ARE NOT THE SAME FOR BOTH AXES.', &
                  'FOR: ' // tao_curve_name(curve))
    return
  endif

  n_symb = count(gx .and. gy)

  call re_allocate (curve%x_symb, n_symb)
  call re_allocate (curve%y_symb, n_symb)
  call re_allocate (curve%ix_symb, n_symb)

  ! Transfer the values

  curve%x_symb = pack (x, mask = gx .and. gy)
  curve%y_symb = pack (y, mask = gx .and. gy)

  ! Calc symbol index

  if (curve%data_index == '') then
    curve%ix_symb = (/ (i, i = 1, n_symb) /)
  else
    call tao_data_type_substitute (curve%data_index, name, curve%ele_ref_name, graph%component)
    call tao_to_real_vector (name, 'BOTH', 0, .true., x, gix, err)
    if (size(gx) == size(gy)) then
      curve%ix_symb = pack (nint(x), mask = gx .and. gy)
    else
      call out_io (s_error$, r_name, &
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

graph%valid = .true.

end subroutine tao_data_slice_graph_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_data_type_substitute (template, str_out, ref_name, component)

implicit none

character(*) template, str_out, ref_name, component
integer ix

!

str_out = template

do
  ix = index(str_out, '#ref')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(ref_name) // trim(str_out(ix+4:))
enddo

do
  ix = index(str_out, '#comp')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(component) // trim(str_out(ix+5:))
enddo

if (index(str_out, '|') == 0) str_out = trim(str_out) // '|' // component

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_phase_space_graph_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (beam_struct), pointer :: beam
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_x, d1_y

real(rp) v_mat(4,4), v_inv_mat(4,4), g_mat(4,4), g_inv_mat(4,4)
real(rp) mat4(4,4), sigma_mat(4,4), theta, theta_xy, rx, ry, phi
real(rp) emit_a, emit_b

integer k, n, m, ib, ix1_ax, ix2_ax, ix, i

logical err, same_uni

character(40) name
character(40) :: r_name = 'tao_phase_space_graph_setup'

! Set up the graph suffix

graph%valid = .false.
if (size(graph%curve) == 0) return

same_uni = .true.
ix = tao_universe_number(graph%curve(1)%ix_universe)
do k = 2, size(graph%curve)
  curve => graph%curve(k)
  if (tao_universe_number(curve%ix_universe) /= ix) same_uni = .false.
enddo

graph%title_suffix = ''
do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (curve%ix_universe)
  if (.not. associated(u)) return
  if (curve%ix_ele_ref_track < 0) then
    call out_io (s_error$, r_name, &
                'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name, &
                'CANNOT PLOT PHASE SPACE FOR: ' // tao_curve_name(curve))
    return
  endif
  ele => u%model%lat%ele(curve%ix_ele_ref_track)
  name = curve%ele_ref_name
  if (name == ' ') name = ele%name
  if (same_uni) then
    write (graph%title_suffix, '(2a, i0, 3a)') trim(graph%title_suffix), &
                                '[', curve%ix_ele_ref, ': ', trim(name), ']'
  else
    write (graph%title_suffix, '(2a, i0, a, i0, 3a)') trim(graph%title_suffix), &
            '[', u%ix_uni, '@', curve%ix_ele_ref, ': ', trim(name), ']'
  endif
enddo

! loop over all curves

do k = 1, size(graph%curve)

  curve => graph%curve(k)
  u => tao_pointer_to_universe (curve%ix_universe)

  ! find phase space axes to plot

  err = .false.
  call tao_phase_space_axis (curve%data_type_x, ix1_ax, err); if (err) return
  call tao_phase_space_axis (curve%data_type,   ix2_ax, err); if (err) return

  ! fill the curve data arrays

  if (allocated (curve%ix_symb)) deallocate (curve%ix_symb, curve%x_symb, curve%y_symb)
  if (allocated (curve%x_line))  deallocate (curve%x_line, curve%y_line)

  if (curve%data_source == 'beam') then
    beam => u%ele(curve%ix_ele_ref_track)%beam
    if (.not. allocated(beam%bunch)) then
      call out_io (s_abort$, r_name, 'NO ALLOCATED BEAM WITH PHASE_SPACE PLOTTING.')
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      return
    endif

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1,  size(beam%bunch)
        n = n + size(beam%bunch(ib)%particle)
      enddo
    else
      n = size(beam%bunch(curve%ix_bunch)%particle)
    endif

    call re_allocate (curve%ix_symb, n)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1, size(beam%bunch)
        m = size(beam%bunch(ib)%particle)
        curve%x_symb(n+1:n+m) = beam%bunch(ib)%particle(:)%r%vec(ix1_ax)
        curve%y_symb(n+1:n+m) = beam%bunch(ib)%particle(:)%r%vec(ix2_ax)
        forall (i = 1:m) curve%ix_symb(n+i) = i
        n = n + m
      enddo
    else
      curve%x_symb = beam%bunch(curve%ix_bunch)%particle(:)%r%vec(ix1_ax)
      curve%y_symb = beam%bunch(curve%ix_bunch)%particle(:)%r%vec(ix2_ax)
      forall (i = 1:m) curve%ix_symb(i) = i
    endif

  !----------------------------

  elseif (curve%data_source == 'multi_turn_orbit') then
    
    call tao_find_data (err, curve%data_source, d2_ptr, ix_uni = curve%ix_universe)
    if (err) then
      call out_io (s_error$, r_name, &
                'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type)
      graph%valid = .false.
      return
    endif

    nullify (d1_x, d1_y)
    do i = 1, size(d2_ptr%d1)
      if (curve%data_type_x == d2_ptr%d1(i)%name) d1_x => d2_ptr%d1(i)
      if (curve%data_type   == d2_ptr%d1(i)%name) d1_y => d2_ptr%d1(i)
    enddo
    if (.not. associated(d1_x)) then
      call out_io (s_error$, r_name, &
              'CANNOT FIND DATA FOR PHASE SPACE COORDINATE: ' // curve%data_type_x, &
              'FOR CURVE: ' // curve%name)
      call err_exit
    endif
    if (.not. associated(d1_y)) then
      call out_io (s_error$, r_name, &
              'CANNOT FIND DATA FOR PHASE SPACE COORDINATE: ' // curve%data_type, &
              'FOR CURVE: ' // curve%name)
      call err_exit
    endif

    if (lbound(d1_x%d, 1) /= lbound(d1_y%d, 1) .or. &
                                        ubound(d1_x%d, 1) /= ubound(d1_y%d, 1)) then 
      call out_io (s_error$, r_name, &
              'BOUNDS FOR X-AXIS AND Y-AXIS DATA OF PHASE SPACE PLOTTING MISMATCHED.', &
              'FOR CURVE: ' // curve%name)
      call err_exit
    endif

    n = size(d1_x%d)
    call re_allocate (curve%ix_symb, n)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)

    do ib = 1, n
      i = ib + lbound(d1_x%d, 1) - 1
      curve%x_symb(ib) = d1_x%d(i)%model_value
      curve%y_symb(ib) = d1_y%d(i)%model_value
    enddo


  elseif (curve%data_source == 'twiss') then

    n = 2 * s%plot_page%n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    call make_v_mats (ele, v_mat, v_inv_mat)
    call make_g_mats (ele, g_mat, g_inv_mat)

    mat4 = matmul(v_mat, g_inv_mat)
    emit_a = u%model%lat%a%emit
    if (emit_a == 0) emit_a = 1e-6  ! default value
    emit_b = u%model%lat%b%emit
    if (emit_b == 0) emit_b = 1e-6  ! default value

    sigma_mat =  0
    sigma_mat(1,1) = emit_a
    sigma_mat(2,2) = emit_a
    sigma_mat(3,3) = emit_b
    sigma_mat(4,4) = emit_b
    sigma_mat = matmul (matmul (mat4, sigma_mat), transpose(mat4))

    if (ix1_ax > 4 .or. ix2_ax > 4) then
      call out_io (s_error$, r_name, &
        'Z OR P_Z PHASE SPACE PLOTTING NOT YET IMPLEMENTED FOR "twiss" DATA_SOURCE.')
      return
    endif

    rx = sqrt(sigma_mat(ix1_ax, ix1_ax))
    ry = sqrt(sigma_mat(ix2_ax, ix2_ax))
    write (graph%text_legend(1), '(a, es9.2)') 'emit_a:', emit_a
    write (graph%text_legend(2), '(a, es9.2)') 'emit_b:', emit_b

    if(rx == 0 .or. ry == 0) then
      theta_xy = 0
      write (graph%text_legend(3), '(a, f10.4)') 'Theta_tilt (rad):', 0
    else
      theta_xy =  asin(sigma_mat(ix1_ax, ix2_ax) / (rx * ry))
      phi = 0.5 *atan2((rx**2+ry**2) * sin(2*theta_xy), &
                              (rx**2-ry**2) * cos(2*theta_xy)) - theta_xy
      write (graph%text_legend(3), '(a, f10.4)') 'Theta_tilt (rad):', phi
  endif

    n = 2 * s%plot_page%n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    do i = 1, n
      theta = (i-1) * twopi / (n-1)
      curve%x_line(i) = rx * cos(theta)
      curve%y_line(i) = ry * sin(theta + theta_xy)
    enddo

  else
    call out_io (s_abort$, r_name, &
        'INVALID CURVE%DATA_SOURCE: ' // curve%data_source, &
        'FOR CURVE: '// curve%name)
    call err_exit
  endif

enddo

graph%valid = .true.

end subroutine tao_phase_space_graph_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_phase_space_axis (data_type, ix_axis, err)

implicit none

character(*) data_type
integer ix_axis
logical err
character(16) :: r_name = 'phase_space_axis'

!

select case (data_type)
case ('x');   ix_axis = 1
case ('p_x'); ix_axis = 2
case ('y');   ix_axis = 3
case ('p_y'); ix_axis = 4
case ('z');   ix_axis = 5
case ('p_z'); ix_axis = 6
case default
  call out_io (s_abort$, r_name, 'BAD PHASE_SPACE CURVE DATA_TYPE: ' // data_type)
  call err_exit
  err = .true.
end select

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_data_graph_setup (plot, graph)

use tao_data_mod
use nrutil, only: swap

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: model_lat, base_lat
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_var_struct), pointer :: v_ptr
type (ele_struct), pointer :: ele, ele1, ele2
type (tao_data_var_component_struct), allocatable, save :: comp(:)

real(rp) f, eps, gs, l_tot, s0, s1, x_max, x_min
real(rp), allocatable :: value(:)

integer ii, k, m, n, n_dat, ie, jj, iv, ic
integer ix, ir, jg, i, j, ix_this
integer, allocatable, save :: ix_ele(:)

logical err, smooth_curve, found, zero_average_phase, ok, ref_or_meas, valid
logical, allocatable, save :: good(:)

character(60) data_type
character(12) :: u_view_char
character(30) :: r_name = 'tao_data_graph_setup'

!

graph%valid = .false.
call re_allocate (ix_ele, 1)

graph%title_suffix = '[' // trim(graph%component) // ']'

! attach x-axis type to title suffix 
 
 if (plot%x_axis_type == 'index') then
   graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: index'
 elseif (plot%x_axis_type == 'ele_index') then
   graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: ele_index'
 elseif (plot%x_axis_type == 's') then
   graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: s'
 endif  

!-------------------------------------------------------------------------------
! Loop over all curves in the graph

do k = 1, size(graph%curve)

  curve => graph%curve(k)

  u => tao_pointer_to_universe (curve%ix_universe)
  if (.not. associated(u)) then
    graph%why_invalid = 'NO ASSOCIATED UNIVERSE!'
    return
  endif

  if (tao_com%common_lattice) then
    call tao_lattice_calc (ok, u%ix_uni, model$)
  endif

  model_lat => u%model%lat
  base_lat => u%base%lat

  if (curve%ele_ref_name == ' ') then
    zero_average_phase = .true.
    curve%ix_ele_ref = 0
    curve%ix_ele_ref_track = 0
  else
    zero_average_phase = .false.
    call tao_locate_elements (curve%ele_ref_name, curve%ix_universe, ix_ele, .true.)
    if (ix_ele(1) < 0) then
      graph%why_invalid = 'CANNOT LOCATE ELEMENT: ' // trim(curve%ele_ref_name)
      return
    endif
    curve%ix_ele_ref = ix_ele(1)
    call tao_ele_ref_to_ele_ref_track(curve%ix_universe, ix_ele(1), curve%ix_ele_ref_track)
  endif



  !----------------------------------------------------------------------------
  ! Calculate where the symbols are to be drawn on the graph.

  select case (curve%data_source)

  !----------------------------------------------------------------------------
  ! Case: data_source is a data_array

  case ('data_array')

    ! point d1_array to the data to be plotted

    call tao_find_data (err, curve%data_type, d2_ptr, d1_array, ix_uni = curve%ix_universe)
    if (err .or. size(d1_array) /= 1) then
      graph%why_invalid = 'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type
      return
    endif

    if (d2_ptr%name == 'phase' .or. d2_ptr%name == 'bpm_phase') then
      if (all(d1_array(1)%d1%d(:)%ele0_name == '')) then
        zero_average_phase = .true.
      else
        zero_average_phase = .false.
      endif
    endif

    ! Set %good_plot True for all data that is within the x-axis limits.
    ! For a circular lattice "wrap around" at s = 0 may mean 
    !   some data points show up twice.

    d1_ptr => d1_array(1)%d1
    d1_ptr%d%good_plot = .false.
    if (graph%x%min /= graph%x%max) then
      eps = 1e-4 * (graph%x%max - graph%x%min)
      if (plot%x_axis_type == 'index') then
        where (d1_ptr%d%ix_d1 > graph%x%min-eps .and. &
               d1_ptr%d%ix_d1 < graph%x%max+eps) d1_ptr%d%good_plot = .true.
      elseif (plot%x_axis_type == 'ele_index') then
        where (d1_ptr%d%ix_ele > graph%x%min-eps .and. &
               d1_ptr%d%ix_ele < graph%x%max+eps) d1_ptr%d%good_plot = .true.
      else ! s
        where (d1_ptr%d%s > graph%x%min-eps .and. &
               d1_ptr%d%s < graph%x%max+eps) d1_ptr%d%good_plot = .true.
        if (model_lat%param%lattice_type == circular_lattice$) then 
          l_tot = model_lat%param%total_length
          where (d1_ptr%d%s-l_tot > graph%x%min-eps .and. &
                 d1_ptr%d%s-l_tot < graph%x%max+eps) d1_ptr%d%good_plot = .true.
        endif

      endif
    endif

    ! make sure %useit_plot up-to-date & count the number of data points

    call tao_data_useit_plot_calc (graph, d1_ptr%d) 
    if (plot%x_axis_type == 's') then
      ! veto non-regular elements when plotting s
      forall (m = lbound(d1_ptr%d,1):ubound(d1_ptr%d,1), &
                     d1_ptr%d(m)%ix_ele > model_lat%n_ele_track)
        d1_ptr%d(m)%useit_plot = .false.
      end forall
    endif
    n_dat = count (d1_ptr%d%useit_plot)       

    ! resize the curve data arrays

    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
    call re_allocate (curve%x_symb, n_dat) ! allocate space for the data

    ! 

    curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)

    if (plot%x_axis_type == 'index') then
      curve%x_symb = curve%ix_symb
    elseif (plot%x_axis_type == 'ele_index') then
      curve%x_symb = d1_ptr%d(curve%ix_symb)%ix_ele
    elseif (plot%x_axis_type == 's') then
      curve%x_symb = model_lat%ele(d1_ptr%d(curve%ix_symb)%ix_ele)%s
      ! If there is a wrap-around then reorder data
      if (model_lat%param%lattice_type == circular_lattice$) then
        do i = 1, n_dat
          if (curve%x_symb(i) > graph%x%max+eps .and. & 
                        curve%x_symb(i)-l_tot > graph%x%min-eps) then
            curve%ix_symb = (/ curve%ix_symb(i:), curve%ix_symb(:i-1) /)
            curve%x_symb = (/ curve%x_symb(i:)-l_tot, curve%x_symb(:i-1) /)
            exit
          endif
        enddo
      endif
      ! Super lords will be out of order so reorder in increasing s.
      do i = 2, n_dat
        do j = i, 2, -1
          if (curve%x_symb(j-1) > curve%x_symb(j)) then
            call swap(curve%x_symb(j-1), curve%x_symb(j))
            call swap(curve%ix_symb(j-1), curve%ix_symb(j))
          else
            exit
          endif
        enddo
      enddo

    else
      graph%why_invalid = 'UNKNOWN AXIS TYPE!'
      return
    endif

    ! calculate the y-axis data point values.

    data_type = trim(curve%data_type) // '|' // graph%component
    call tao_to_real_vector (data_type, 'DATA', 0, .true., value, good, err)
    if (err) then
      graph%why_invalid = 'BAD PLOT COMPONENT: ' // data_type
      return
    end if

    curve%y_symb = pack(value, mask = d1_ptr%d%useit_plot)

  !----------------------------------------------------------------------------
  ! Case: data_source is a var_array

  case ('var_array')
    call tao_find_var (err, curve%data_type, v1_array)
    if (err .or. size(v1_array) /= 1) return
    v1_ptr => v1_array(1)%v1

    ! find which universe we're viewing
    ix_this = -1
    v_loop: do iv = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
      v_ptr => v1_ptr%v(iv)
      if (.not. v_ptr%exists) cycle
      do jj = 1, size(v_ptr%this)
        if (v_ptr%this(jj)%ix_uni .eq. s%global%u_view) then
          ix_this = jj
          exit v_loop
        endif
      enddo
    enddo v_loop
    if (ix_this .eq. -1) then
      call out_io (s_error$, r_name, &
                     "This variable doesn't point to the currently displayed  universe.")
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
          if (v1_ptr%v(jj)%this(ix_this)%ix_ele < graph%x%min-eps) v1_ptr%v%good_plot = .false.
          if (v1_ptr%v(jj)%this(ix_this)%ix_ele > graph%x%max+eps) v1_ptr%v%good_plot = .false.
        enddo
      else
        where (v1_ptr%v%s < graph%x%min-eps) v1_ptr%v%good_plot = .false.
        where (v1_ptr%v%s > graph%x%max+eps) v1_ptr%v%good_plot = .false.
      endif
    endif

    call tao_var_useit_plot_calc (graph, v1_ptr%v) ! make sure %useit_plot up-to-date
    n_dat = count (v1_ptr%v%useit_plot)       ! count the number of data points

    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
    call re_allocate (curve%x_symb, n_dat) ! allocate space for the data

    curve%ix_symb = pack(v1_ptr%v%ix_v1, mask = v1_ptr%v%useit_plot)

    graph%x%label = plot%x_axis_type

    if (plot%x_axis_type == 'index') then
      curve%x_symb = curve%ix_symb
    elseif (plot%x_axis_type == 'ele_index') then
      do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
        curve%x_symb(jj) = v1_ptr%v(curve%ix_symb(jj))%this(ix_this)%ix_ele
      enddo
    elseif (plot%x_axis_type == 's') then
      do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
        curve%x_symb(jj) = model_lat%ele(v1_ptr%v(curve%ix_symb(jj))%this(ix_this)%ix_ele)%s
      enddo
    endif

    ! calculate the y-axis data point values.

    data_type = trim(curve%data_type) // '|' // graph%component
    call tao_to_real_vector (data_type, 'VAR', 0, .true., value, good, err)
    if (err) then
      call out_io (s_error$, r_name, 'BAD PLOT COMPONENT: ' // data_type)
      return
    end if

    curve%y_symb = pack(value, mask = d1_ptr%d%useit_plot)

  !----------------------------------------------------------------------------
  ! Case: data_source is from lattice, or beam

  case ('lattice', 'beam')

    ! Find how many symbol points there are

    if (plot%x_axis_type == 'index' .or. plot%x_axis_type == 'ele_index') then
      x_min = 1
      x_max = model_lat%n_ele_track
      if (graph%x%min /= graph%x%max) then
        x_min = max(x_min, graph%x%min)
        x_max = min(x_max, graph%x%max)
      endif 
      n_dat = max(0, nint(x_max+1-x_min))
      call re_allocate (ix_ele, n_dat)
      if (n_dat > 0) ix_ele = (/ (i, i = nint(x_min), nint(x_max)) /)

    elseif (plot%x_axis_type == 's') then
      call re_allocate (ix_ele, 0)

    else
      call out_io (s_error$, r_name, 'BAD PLOT%X_AXIS_TYPE: ' // plot%x_axis_type)
      call err_exit
    endif

    call tao_curve_datum_calc (ix_ele, curve, 'SYMBOL', valid)
    if (.not. valid) return

  !----------------------------------------------------------------------------
  ! Case: Bad data_source

  case default
    call out_io (s_error$, r_name, 'UNKNOWN DATA_SOURCE: ' // curve%data_source)
    return
  end select

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Now calculate the points for drawing the curve through the symbols...

  ! If the x-axis is by index or ele_index then these points are the same as the symbol points.
  ! That is, for x-axis = 'index' or 'ele_index' the line is piece-wise linear between the symbols.

  select case (plot%x_axis_type)
  case ('index', 'ele_index')
    call re_allocate (curve%y_line, n_dat) ! allocate space for the data
    call re_allocate (curve%x_line, n_dat) ! allocate space for the data
    curve%x_line = curve%x_symb
    curve%y_line = curve%y_symb

  ! If the axis is by s-value then, if possible, the line is a "smooth" curve with n_curve_pts points.

  case ('s')

    smooth_curve = (curve%data_source == 'lattice') .or. &
                   (curve%data_source == 'beam' .and. allocated(u%model%bunch_params2))
    smooth_curve = smooth_curve .and. curve%smooth_line_calc

    if (curve%data_source == 'lattice' .and. index(curve%data_type, 'emit.') /= 0) &
                                                                       smooth_curve = .false.

    if (index(graph%component, 'meas') /= 0 .or. index(graph%component, 'ref') /= 0) then
      ref_or_meas = .true.
      smooth_curve = .false.
    else
      ref_or_meas = .false.
    endif

    if (smooth_curve) then

      ! allocate data space

      call re_allocate (curve%y_line, s%plot_page%n_curve_pts) 
      call re_allocate (curve%x_line, s%plot_page%n_curve_pts) 
      call re_allocate (good,         s%plot_page%n_curve_pts) 
      curve%y_line = 0
      good = .true.

      call tao_split_component(graph%component, comp, err)
      if (err) return
      do m = 1, size(comp)
        select case (comp(m)%name)
        case (' ') 
          cycle
        case ('model')
          call calc_data_at_s (u%model, curve, comp(m)%sign, good)
        case ('base')  
          call calc_data_at_s (u%base, curve, comp(m)%sign, good)
        case ('design')  
          call calc_data_at_s (u%design, curve, comp(m)%sign, good)
        case default
          call out_io (s_error$, r_name, &
                       'BAD PLOT COMPONENT WITH "S" X-AXIS: ' // comp(m)%name)
          return
        end select
      enddo

      if (all(.not. good)) return
      n_dat = count(good)
      curve%x_line(1:n_dat) = pack(curve%x_line, mask = good)
      curve%y_line(1:n_dat) = pack(curve%y_line, mask = good)
      call re_allocate (curve%y_line, n_dat) ! allocate space for the data
      call re_allocate (curve%x_line, n_dat) ! allocate space for the data

    ! For x_axis_type = 's' and non smooth curves: Just evaluate at the element ends.

    else if (plot%x_axis_type == 's' .and. .not. ref_or_meas) then

      eps = 1e-4 * (graph%x%max - graph%x%min)             ! a small number
      l_tot = model_lat%param%total_length
      model_lat%ele%logic = .false.
      do i = 1, model_lat%n_ele_track
        ele => model_lat%ele(i)
        if (graph%x%min == graph%x%max) cycle
        s0 = ele%s - ele%value(l$)
        s1 = ele%s
        ele%logic = (s0 >= graph%x%min-eps) .and. (s1 <= graph%x%max+eps)
        if (model_lat%param%lattice_type == circular_lattice$) then
          ele%logic = ele%logic .or. &
                     ((s0-l_tot >= graph%x%min-eps) .and. (s1-l_tot <= graph%x%max+eps))
        endif
      enddo
      n_dat = count (model_lat%ele(:)%logic)
      call re_allocate (ix_ele, n_dat)
      ix_ele = pack(model_lat%ele(:)%ix_ele, mask = model_lat%ele(:)%logic)
      ! If there is a wrap-around then reorder the data
      do i = 1, n_dat
        if (model_lat%ele(ix_ele(i))%s - l_tot > graph%x%min) then
          ix_ele = (/ ix_ele(i:), ix_ele(:i-1) /)
          exit
        endif
      enddo

      call tao_curve_datum_calc (ix_ele, curve, 'LINE', graph%valid)
      if (.not. graph%valid) return

      do i = 1, size(curve%x_line)
        curve%x_line(i) = model_lat%ele(ix_ele(i))%s
      enddo

    ! For all else just draw straight lines through the symbols.

    else
      ! allocate space for the data
      call re_allocate (curve%y_line, n_dat) 
      call re_allocate (curve%x_line, n_dat) 
      curve%x_line = curve%x_symb 
      curve%y_line = curve%y_symb 
    endif

  end select

  !----------------------------------------------------------------------------
  ! Note: Since there is an arbitrary overall phase, the phase data 
  ! gets renormalized so that the average value is zero.

  if ((curve%data_type(1:6) == 'phase.' .or. curve%data_type(1:10) == 'bpm_phase.') &
                                      .and. n_dat /= 0 .and. zero_average_phase) then
    f = sum(curve%y_symb) / n_dat
    curve%y_symb = curve%y_symb - f
    curve%y_line = curve%y_line - f 
  endif 

  graph%valid = .true.

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_data_at_s (tao_lat, curve, comp_sign, good)

use transfer_map_mod

implicit none

type (tao_lattice_struct), target :: tao_lat
type (tao_curve_struct) curve
type (bunch_params_struct), pointer :: bunch_params
type (coord_struct), pointer :: orb(:)
type (coord_struct), pointer :: orb0
type (lat_struct), pointer :: lat
type (ele_struct) ele
type (ele_struct), pointer :: ele0
type (coord_struct) here
type (taylor_struct) t_map(6)

real(rp) x1, x2, cbar(2,2), s_last, s_now, value, mat6(6,6), vec0(6)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), one_pz, gamma, len_tot
real(rp) comp_sign

integer i, ii, ix, j, k, expnt(6), ix_ele, ix0

character(40) data_type
character(40) data_type_select, data_source
character(20) :: r_name = 'calc_data_at_s'
logical err, good(:)

!

data_type = curve%data_type

lat => tao_lat%lat
orb => tao_lat%orb
ix0 = curve%ix_ele_ref_track
if (ix0 < 0) ix0 = 0

if (lat%param%lattice_type == circular_lattice$ .and. .not. lat%param%stable) then
  good = .false.
  return
endif

if (curve%data_source == 'lattice') then
  select case (data_type(1:5))
  case ('sigma', 'emitt', 'norm_')
    call out_io (s_fatal$, r_name, &
              'CURVE%DATA_SOURCE = "lattice" IS NOT COMPATABLE WITH DATA_TYPE: ' &
              // data_type)
    call out_io (s_blank$, r_name, "Will not perfrom any plot smoothing")
    good = .false.
    return
  end select 
endif

! x1 and x2 are the longitudinal end points of the plot

x1 = lat%ele(0)%s
x2 = lat%ele(lat%n_ele_track)%s
len_tot = lat%param%total_length
if (curve%g%x%min /= curve%g%x%max) then
  if (lat%param%lattice_type == circular_lattice$) then
    x1 = min(lat%ele(lat%n_ele_track)%s, max(curve%g%x%min, x1-len_tot))
    x2 = min(x2, max(curve%g%x%max, lat%ele(0)%s-len_tot))
  else
    x1 = min(lat%ele(lat%n_ele_track)%s, max(curve%g%x%min, x1))
    x2 = min(x2, max(curve%g%x%max, lat%ele(0)%s))
  endif
endif
ele0 => lat%ele(ix0)
orb0 => orb(ix0)
s_last = ele0%s

data_type_select = data_type
if (data_type_select(1:2) == 'r.') data_type_select = 'r.'
if (data_type_select(1:2) == 't.') data_type_select = 't.'
if (data_type_select(1:3) == 'tt.') data_type_select = 'tt.'

!

do ii = 1, size(curve%x_line)

  if (.not. good(ii)) cycle

  s_now = x1 + (ii-1) * (x2-x1) / (size(curve%x_line)-1)
  if (s_now .ge. lat%ele(lat%n_ele_track)%s) s_now = lat%ele(lat%n_ele_track)%s - 1e-9
  curve%x_line(ii) = s_now
  value = 0

  select case (curve%data_source)
  case ('lattice')   
    call twiss_and_track_at_s (lat, s_now, ele, orb, here, err)
    if (err) then
      good(ii:) = .false.
      return
    endif

  case ('beam')
    if (.not. allocated(tao_lat%bunch_params2)) then
      call out_io (s_fatal$, r_name, 'BUNCH_PARAMS2 NOT ALLOCATED.')
      call err_exit
    endif
 
    call bracket_index (tao_lat%bunch_params2(:)%s, 1, tao_lat%n_bunch_params2, s_now, ix)
    if (abs(tao_lat%bunch_params2(ix)%s - s_now) < abs(tao_lat%bunch_params2(ix+1)%s - s_now)) then
      bunch_params => tao_lat%bunch_params2(ix)
    else
      bunch_params => tao_lat%bunch_params2(ix+1)
    endif

    call ele_at_s (lat, s_now, ix_ele)
    ele = lat%ele(ix_ele)
    here = bunch_params%centroid

  case default
    call out_io (s_fatal$, r_name, &
            'I DO NOT KNOW HOW TO HANDLE THIS curve%data_source: ' // curve%data_source)
    call err_exit
  end select

!

  select case (data_type_select)
  case ('orbit.x')
    value = here%vec(1)
  case ('orbit.y')
    value = here%vec(3)
  case ('orbit.z')
    value = here%vec(5)
  case ('orbit.p_x')
    value = here%vec(2)
  case ('orbit.p_y')
    value = here%vec(4)
  case ('orbit.p_z')
    value = here%vec(6)
  case ('orbit.amp_a')
    call orbit_amplitude_calc (ele, here, amp_a = value)
  case ('orbit.amp_b')
    call orbit_amplitude_calc (ele, here, amp_b = value)
  case ('orbit.norm_amp_a')
    call orbit_amplitude_calc (ele, here, amp_na = value)
  case ('orbit.norm_amp_b')
    call orbit_amplitude_calc (ele, here, amp_nb = value)
  case ('phase.a')
    value = ele%a%phi
  case ('phase.b')
    value = ele%b%phi
  case ('beta.a')
    value = ele%a%beta
  case ('beta.b')
    value = ele%b%beta
  case ('alpha.a')
    value = ele%a%alpha
  case ('alpha.b')
    value = ele%b%alpha
  case ('eta.x')
    value = ele%x%eta
  case ('eta.y')
    value = ele%y%eta
  case ('etap.x')
    value = ele%x%etap
  case ('etap.y')
    value = ele%y%etap
  case ('eta.a')
    value = ele%a%eta
  case ('eta.b')
    value = ele%b%eta
  case ('etap.a')
    value = ele%a%etap
  case ('etap.b')
    value = ele%b%etap
  case ('floor.x')
    value = ele%floor%x
  case ('floor.y')
    value = ele%floor%y
  case ('floor.z')
    value = ele%floor%z
  case ('e_tot')
    value = ele%value(e_tot$) * (1 + here%vec(6))
  case ('%e_tot')
    value = here%vec(6)
  case ('cbar.11')
    call c_to_cbar (ele, cbar)
    value = cbar(1,1)
  case ('cbar.12')
    call c_to_cbar (ele, cbar)
    value = cbar(1,2)
  case ('cbar.21')
    call c_to_cbar (ele, cbar)
    value = cbar(2,1)
  case ('cbar.22')
    call c_to_cbar (ele, cbar)
    value = cbar(2,2)
  case ('coupling.11b')
    call c_to_cbar (ele, cbar)
    value = cbar(1,1) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
  case ('coupling.12a')
    call c_to_cbar (ele, cbar)
    value = cbar(1,2) * sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
  case ('coupling.12b')
    call c_to_cbar (ele, cbar)
    value = cbar(1,2) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
  case ('coupling.22a')
    call c_to_cbar (ele, cbar)
    value = cbar(2,2)* sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
  case ('sigma.x')
    value = sqrt(bunch_params%sigma(s11$))
  case ('sigma.p_x')
    value = sqrt(bunch_params%sigma(s22$))
  case ('sigma.y')
    value = sqrt(bunch_params%sigma(s33$))
  case ('sigma.p_y')
    value = sqrt(bunch_params%sigma(s44$))
  case ('sigma.z')
    value = sqrt(bunch_params%sigma(s55$))
  case ('sigma.p_z')
    value = sqrt(bunch_params%sigma(s66$))
  case ('norm_emit.a')
    value = bunch_params%a%norm_emitt
  case ('norm_emit.b')
    value = bunch_params%b%norm_emitt
  case ('norm_emit.x')
    value = bunch_params%x%norm_emitt
  case ('norm_emit.y')
    value = bunch_params%y%norm_emitt
  case ('norm_emit.z')
    value = bunch_params%z%norm_emitt
  case ('emit.a')
    value = bunch_params%a%norm_emitt
    call convert_total_energy_to (ele%value(E_tot$), &
                                              lat%param%particle, gamma)
    value = value / gamma
  case ('emit.b')
    value = bunch_params%b%norm_emitt
    call convert_total_energy_to (ele%value(E_tot$), &
                                              lat%param%particle, gamma)
    value = value / gamma
  case ('r.')
    if (ii == 1) call mat_make_unit (mat6)
    if (s_now < s_last) cycle
    i = tao_read_this_index (data_type, 3); if (i == 0) return
    j = tao_read_this_index (data_type, 4); if (j == 0) return
    call mat6_calc_at_s (lat, mat6, vec0, s_last, s_now, unit_start = .false.)
    value = mat6(i, j)
  case ('t.')
    if (ii == 1) call taylor_make_unit (t_map)
    if (s_now < s_last) cycle
    i = tao_read_this_index (data_type, 3); if (i == 0) return
    j = tao_read_this_index (data_type, 4); if (j == 0) return
    k = tao_read_this_index (data_type, 5); if (k == 0) return
    call transfer_map_calc_at_s (lat, t_map, s_last, s_now, unit_start = .false.)
    value = taylor_coef (t_map(i), j, k)
  case ('tt.')
    if (ii == 1) call taylor_make_unit (t_map)
    if (s_now < s_last) cycle
    expnt = 0
    i = tao_read_this_index (data_type, 4); if (i == 0) return
    do j = 5, 15
      if (data_type(j:j) == ' ') exit
      k = tao_read_this_index (data_type, j); if (k == 0) return
      expnt(k) = expnt(k) + 1
    enddo
    call transfer_map_calc_at_s (lat, t_map, s_last, s_now, unit_start = .false.)
    value = taylor_coef (t_map(i), expnt)
  case ('momentum_compaction')
    if (ii == 1) call mat_make_unit (mat6)
    call mat6_calc_at_s (lat, mat6, vec0, s_last, s_now, unit_start = .false.)
    call make_v_mats (ele0, v_mat, v_inv_mat)
    eta_vec = (/ ele0%a%eta, ele0%a%etap, ele0%b%eta, ele0%b%etap /)
    eta_vec = matmul (v_mat, eta_vec)
    one_pz = 1 + orb0%vec(6)
    eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
    eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz
    value = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)

  case default
    call out_io (s_fatal$, r_name, &
                  'DO NOT KNOW ABOUT THIS DATA_TYPE: ' // data_type)
    call out_io (s_blank$, r_name, "Will not perfrom any plot smoothing")
    good = .false.
    return
  end select

  curve%y_line(ii) = curve%y_line(ii) + comp_sign * value
  s_last = s_now

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_useit_plot_calc (graph, data)
!
! Subroutine to set the data for plotting.
!
! Input:
!
! Output:
!   data     -- Tao_data_struct:
!     %useit_plot -- True if good for plotting.
!                  = %exists & %good_plot (w/o measured & reference data)
!                  = %exists & %good_plot & %good_user & %good_meas (w/ meas data)
!                  = %exists & %good_plot & %good_user & %good_ref (w/ reference data)
!                  = %exists & %good_plot & %good_user & %good_meas & %good_ref 
!                                                        (w/ measured & reference data)
!-

subroutine tao_data_useit_plot_calc (graph, data)

implicit none

type (tao_graph_struct) graph
type (tao_data_struct) data(:)

!

data%useit_plot = data%exists .and. data%good_plot .and. data%good_user
if (index(graph%component, 'meas') /= 0) &
         data%useit_plot = data%useit_plot .and. data%good_meas
if (index(graph%component, 'ref') /= 0)  &
         data%useit_plot = data%useit_plot .and. data%good_ref
if (index(graph%component, 'model') /= 0)  &
         data%useit_plot = data%useit_plot .and. data%good_model

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_curve_datum_calc (ix_ele, curve, who, valid)
!
! Input:
!   curve -- Tao_curve_struct:
!   who   -- Character(*): Where to put the data. 
!               Either: "SYMBOL" or "LINE".
!
! Output:
!   curve -- Tao_curve_struct: 
!   valid -- Logical: 
!-

subroutine tao_curve_datum_calc (ix_ele, curve, who, valid)

implicit none

type (tao_curve_struct) curve
type (tao_data_var_component_struct), allocatable, save :: comp(:)
type (tao_universe_struct), pointer :: u
type (tao_data_struct) datum
type (taylor_struct) t_map(6)

integer, allocatable :: ix_ele(:)
real(rp), allocatable, save :: y_value(:)
real(rp) y_val

logical, allocatable, save :: good(:)
logical valid, err

character(*) who
character(20) :: r_name = 'tao_curve_datum_calc'
character(40) why_invalid

integer m, ie, n_dat

! calculate the y-axis data point values.

u => tao_pointer_to_universe (curve%ix_universe)
n_dat = size(ix_ele)

call re_allocate (good, n_dat)   ! allocate space for the data
call re_allocate (y_value, n_dat) ! allocate space for the data

y_value = 0
good = .true.

datum%ix_ele0     = curve%ix_ele_ref_track
datum%merit_type  = 'target'
datum%data_type   = curve%data_type
datum%ele0_name   = curve%ele_ref_name
datum%data_source = curve%data_source

call tao_split_component (curve%g%component, comp, err)
if (err) return
do m = 1, size(comp)

  do ie = 1, n_dat

    if (datum%data_type(1:3) == 'tt.' .or. datum%data_type(1:2) == 't.') then
      if (ie == 1) then
        call taylor_make_unit (t_map)
      else
        datum%ix_ele0 = datum%ix_ele
      endif
    endif

    datum%ix_ele = ix_ele(ie)

    select case (comp(m)%name)
    case (' ') 
      cycle
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, y_val, valid, why_invalid, t_map)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, y_val, valid, why_invalid, t_map)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, y_val, valid, why_invalid, t_map)
    case ('ref', 'meas')
      call out_io (s_error$, r_name, &
              'PLOT COMPONENT WHICH IS: ' // comp(m)%name, &
              '    FOR DATA_TYPE: ' // curve%data_type, &
              '    NOT ALLOWED SINCE DATA_SOURCE IS SET TO: ' // curve%data_source)
      return
    case default
      call out_io (s_error$, r_name, &
              'BAD PLOT COMPONENT: ' // comp(m)%name, &
              '    FOR DATA_TYPE: ' // curve%data_type)
      return
    end select
    y_value(ie) = y_value(ie) + comp(m)%sign * y_val
    if (.not. valid) good(ie) = .false.
    if (datum%data_type(1:3) == 'tt.' .or. datum%data_type(1:2) == 't.') then
      if (datum%ix_ele < datum%ix_ele0) datum%ix_ele0 = datum%ix_ele
    endif

  enddo
enddo

if (n_dat > 0 .and. all(.not. good)) then
  valid = .false.
  curve%g%why_invalid = why_invalid
  return
endif

n_dat = count(good)

if (who == 'SYMBOL') then
  call re_allocate (curve%x_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%ix_symb, n_dat)
  curve%x_symb  = pack(ix_ele, mask = good)
  curve%y_symb  = pack(y_value, mask = good)
  curve%ix_symb = curve%x_symb
else
  call re_allocate (curve%x_line, n_dat) ! allocate space for the data
  call re_allocate (curve%y_line, n_dat) ! allocate space for the data
  curve%x_line  = pack(ix_ele, mask = good)
  curve%y_line  = pack(y_value, mask = good)
  ix_ele(1:n_dat) = pack(ix_ele, mask = good)
  call re_allocate (ix_ele, n_dat)
endif

valid = .true.

end subroutine

end module
