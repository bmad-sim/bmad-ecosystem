!+
! Subroutine tao_plot_data_setup ()
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_page%region(:)%plot%graph(:)%curve(:) arrays.
!
! Input/Output:
!-

subroutine tao_plot_data_setup ()

use tao_mod
use tao_data_mod

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_universe_struct), pointer :: u
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_v1_var_struct) , pointer :: v1_ptr
type (tao_data_struct) datum
type (ele_struct), pointer :: ele
type (taylor_struct) t_map(6)

real(rp) f, x1, x2, y_val, eps
real(rp), pointer :: value(:)

integer i, ii, j, k, m, n_dat, i_uni, ie, jj
integer ix_this, ix

logical err, smooth_curve

character(20) :: r_name = 'tao_plot_data_setup'
character(12)  :: u_view_char

! Find which elements are to be drawn for a lattice layout.

if (any(s%plot_page%ele_shape(:)%key /= 0)) then
  do i = 1, size(s%u)
    s%u(i)%model%ele_(:)%ix_pointer = 0
    s%u(i)%base%ele_(:)%ix_pointer = 0
    do j = 1, s%u(i)%model%n_ele_max
      ele => s%u(i)%model%ele_(j)
      if (ele%control_type == group_lord$) cycle
      if (ele%control_type == overlay_lord$) cycle
      if (ele%control_type == super_slave$) cycle
      if (ele%control_type == multipass_lord$) cycle
      do k = 1, size(s%plot_page%ele_shape(:))
        if (s%plot_page%ele_shape(k)%key == 0) cycle
        if (ele%key == s%plot_page%ele_shape(k)%key .and. &
                 match_wild(ele%name, s%plot_page%ele_shape(k)%ele_name)) then
          ele%ix_pointer = k
          s%u(i)%base%ele_(j)%ix_pointer = j
          s%u(i)%base%ele_(j-1)%ix_pointer = j-1
          exit
        endif
      enddo
    enddo
  enddo
endif

! setup the plots

plot_loop: do i = 1, size(s%plot_page%region)

  plot => s%plot_page%region(i)%plot

  if (.not. s%plot_page%region(i)%visible) cycle  ! Don't worry about invisable graphs

  if (plot%x_axis_type /= 'index' .and. plot%x_axis_type /= 's' .and. &
      plot%x_axis_type /= 'ele_index') then
    call out_io (s_abort$, r_name, 'BAD X_AXIS_TYPE: ' // plot%x_axis_type)
    plot%graph%valid = .false.
    cycle
  endif

! loop over all graphs and curves

  graph_loop: do j = 1, size(plot%graph)

    graph => plot%graph(j)
    graph%valid = .true.   ! assume everything OK
    graph%legend = ' '

    if (graph%type /= 'data') cycle ! Don't worry about other types

    do k = 1, size(graph%curve)

      curve => graph%curve(k)

      if (curve%data_source == 'lat_layout' .and. &
            (plot%x_axis_type == 'index' .or. plot%x_axis_type == 'ele_index')) then
        call out_io (s_error$, r_name, 'CURVE%DATA_SOURCE = "LAT_LAYOUT" ' // &
                        'AND PLOT%X_AXIS_TYPE = "INDEX" DO NOT GO TOGETHER.')
        graph%valid = .false.
        cycle graph_loop
      endif

      i_uni = s%global%u_view  ! universe where the data comes from
      if (curve%ix_universe /= 0) i_uni = curve%ix_universe 
      u => s%u(i_uni)
      
      if (curve%ele2_name /= ' ') then
        call element_locator (curve%ele2_name, u%design, curve%ix_ele2)
        if (curve%ix_ele2 < 0) then
          curve%ix_ele2 = 0
          call out_io (s_error$, r_name, &
                  'Curve%ele2_name cannot be found in lattice: ' // curve%ele2_name)
        endif
      endif

!----------------------------------------------------------------------------
! data_source is a data array

      select case (curve%data_source)
      case ('data_array')
        call tao_find_data (err, u, curve%data_type, d2_ptr, d1_ptr)
        if (err) then
          graph%valid = .false.
          cycle graph_loop
        endif
	
        d1_ptr%d%good_plot = .true.
        eps = 1e-4 * (plot%x%max - plot%x%min)
        if (plot%x_axis_type == 'index') then
          where (d1_ptr%d%ix_d1 < plot%x%min-eps) d1_ptr%d%good_plot = .false.
          where (d1_ptr%d%ix_d1 > plot%x%max+eps) d1_ptr%d%good_plot = .false.
        elseif (plot%x_axis_type == 'ele_index') then
          where (d1_ptr%d%ix_ele < plot%x%min-eps) d1_ptr%d%good_plot = .false.
          where (d1_ptr%d%ix_ele > plot%x%max+eps) d1_ptr%d%good_plot = .false.
        else ! s
          where (d1_ptr%d%s < plot%x%min-eps) d1_ptr%d%good_plot = .false.
          where (d1_ptr%d%s > plot%x%max+eps) d1_ptr%d%good_plot = .false.
        endif

        ! make sure %useit_plot up-to-date & count the number of data points
        call tao_data_useit_plot_calc (graph, d1_ptr%d) 
        n_dat = count (d1_ptr%d%useit_plot)       

        call reassociate_integer (curve%ix_symb, n_dat)
        call reassociate_real (curve%y_symb, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_symb, n_dat) ! allocate space for the data

        curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)

        if (plot%x_axis_type == 'index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 'ele_index') then
          curve%x_symb = d1_ptr%d(curve%ix_symb)%ix_ele
        elseif (plot%x_axis_type == 's') then
          curve%x_symb = u%model%ele_(d1_ptr%d(curve%ix_symb)%ix_ele)%s
        else
          call out_io (s_error$, r_name, "Unknown axis type!")
          graph%valid = .false.
          cycle graph_loop
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0

        do m = 1, size(graph%who)
          select case (graph%who(m)%name)
          case (' ') 
            cycle
          case ('model') 
            value => d1_ptr%d%model_value
          case ('base')  
            value => d1_ptr%d%base_value
          case ('design')  
            value => d1_ptr%d%design_value
          case ('ref')     
            value => d1_ptr%d%ref_value
          case ('meas')    
            value => d1_ptr%d%meas_value
          case default
            call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // graph%who(m)%name)
            graph%valid = .false.
            cycle graph_loop
          end select
          curve%y_symb = curve%y_symb + &
                   graph%who(m)%sign * pack(value, mask = d1_ptr%d%useit_plot)
        enddo

        if (curve%convert) curve%y_symb = curve%y_symb * &
                           pack(d1_ptr%d%conversion_factor, d1_ptr%d%useit_plot)


!----------------------------------------------------------------------------
! data_source is a var array

      case ('var_array')
        call tao_find_var (err, curve%data_type, v1_ptr)
        if (err) then
          graph%valid = .false.
          cycle graph_loop
        endif

        ! find which universe we're viewing

        ix_this = -1
        do jj = 1, size(v1_ptr%v(1)%this)
          if (v1_ptr%v(1)%this(jj)%ix_uni .eq. s%global%u_view) ix_this = jj
        enddo
        if (ix_this .eq. -1) then
          call out_io (s_error$, r_name, &
                     "This variable doesn't point to the currently displayed  universe.")
          graph%valid = .false.
          cycle graph_loop
        endif
      
        v1_ptr%v%good_plot = .true.
        eps = 1e-4 * (plot%x%max - plot%x%min)
        if (plot%x_axis_type == 'index') then
          where (v1_ptr%v%ix_v1 < plot%x%min-eps) v1_ptr%v%good_plot = .false.
          where (v1_ptr%v%ix_v1 > plot%x%max+eps) v1_ptr%v%good_plot = .false.
        elseif (plot%x_axis_type == 'ele_index') then
          do jj = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
            if (v1_ptr%v(jj)%this(ix_this)%ix_ele < plot%x%min-eps) v1_ptr%v%good_plot = .false.
            if (v1_ptr%v(jj)%this(ix_this)%ix_ele > plot%x%max+eps) v1_ptr%v%good_plot = .false.
          enddo
        else
          where (v1_ptr%v%s < plot%x%min-eps) v1_ptr%v%good_plot = .false.
          where (v1_ptr%v%s > plot%x%max+eps) v1_ptr%v%good_plot = .false.
        endif

        call tao_var_useit_plot_calc (graph, v1_ptr%v) ! make sure %useit_plot up-to-date
        n_dat = count (v1_ptr%v%useit_plot)       ! count the number of data points

        call reassociate_integer (curve%ix_symb, n_dat)
        call reassociate_real (curve%y_symb, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_symb, n_dat) ! allocate space for the data

        curve%ix_symb = pack(v1_ptr%v%ix_v1, mask = v1_ptr%v%useit_plot)

        plot%x%label = plot%x_axis_type

        if (plot%x_axis_type == 'index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 'ele_index') then
          do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
            curve%x_symb(jj) = v1_ptr%v(curve%ix_symb(jj))%this(ix_this)%ix_ele
          enddo
        elseif (plot%x_axis_type == 's') then
          do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
            curve%x_symb(jj) = u%model%ele_(v1_ptr%v(curve%ix_symb(jj))%this(ix_this)%ix_ele)%s
          enddo
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0

        do m = 1, size(graph%who)
          select case (graph%who(m)%name)
          case (' ') 
            cycle

          ! set value to whatever it is in currently viewed universe
          case ('model') 
            do jj = lbound(v1_ptr%v,1), ubound(v1_ptr%v,1)
              if (associated(v1_ptr%v(jj)%this(ix_this)%model_ptr)) &
                      v1_ptr%v(jj)%plot_model_value = v1_ptr%v(jj)%this(ix_this)%model_ptr
            enddo
            value => v1_ptr%v%plot_model_value

          ! set value to whatever it is in currently viewed universe
          case ('base')  
            do jj = lbound(v1_ptr%v,1), ubound(v1_ptr%v,1)
              if (associated(v1_ptr%v(jj)%this(ix_this)%base_ptr)) &
                      v1_ptr%v(jj)%plot_base_value = v1_ptr%v(jj)%this(ix_this)%base_ptr
            enddo
            value => v1_ptr%v%plot_base_value

          case ('design')  
            value => v1_ptr%v%design_value
          case ('ref')     
            value => v1_ptr%v%ref_value
          case ('meas')    
            value => v1_ptr%v%meas_value
          case default
            call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // graph%who(m)%name)
            graph%valid = .false.
            cycle graph_loop
          end select

          curve%y_symb = curve%y_symb + &
                   graph%who(m)%sign * pack(value, mask = v1_ptr%v%useit_plot)
        enddo

        if (curve%convert) curve%y_symb = curve%y_symb * &
                           pack(v1_ptr%v%conversion_factor, v1_ptr%v%useit_plot)


!----------------------------------------------------------------------------
! data source is from the lattice_layout

      case ('lat_layout')
 
        eps = 1e-4 * (plot%x%max - plot%x%min)
        u%base%ele_(:)%logic = (u%base%ele_(:)%ix_pointer > 0) .and. &
            (u%base%ele_(:)%s >= plot%x%min-eps) .and. (u%base%ele_(:)%s <= plot%x%max+eps)
        n_dat = count (u%base%ele_(:)%logic)

        call reassociate_integer (curve%ix_symb, n_dat)
        call reassociate_real (curve%y_symb, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_symb, n_dat) ! allocate space for the data

        curve%ix_symb = pack(u%base%ele_(:)%ix_pointer, mask = u%base%ele_(:)%logic)

        if (plot%x_axis_type == 'index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 'ele_index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 's') then
          curve%x_symb = u%model%ele_(curve%ix_symb)%s
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0
        datum%ix_ele2 = curve%ix_ele2
        datum%merit_type = 'target'
        datum%data_type = curve%data_type
        datum%ele2_name = curve%ele2_name

        do m = 1, size(graph%who)
          do ie = 1, n_dat

            datum%ix_ele = curve%ix_symb(ie)

            if (datum%data_type(1:3) == 'tt:' .or. datum%data_type(1:2) == 't:') then
              if (ie == 1) call taylor_make_unit (t_map)
            endif

            select case (graph%who(m)%name)
            case (' ') 
              cycle
            case ('model')   
              call tao_evaluate_a_datum (datum, u, u%model, u%model_orb, y_val, t_map)
            case ('base')  
              call tao_evaluate_a_datum (datum, u, u%base, u%base_orb, y_val, t_map)
            case ('design')  
              call tao_evaluate_a_datum (datum, u, u%design, u%design_orb, y_val, t_map)
            case default
              call out_io (s_error$, r_name, &
                          'BAD PLOT "WHO" FOR LAT_LAYOUT DATA_SOURCE: ' // graph%who(m)%name, &
                          '    FOR DATA_TYPE: ' // curve%data_type)
              graph%valid = .false.
              cycle graph_loop
            end select
            curve%y_symb(ie) = curve%y_symb(ie) + graph%who(m)%sign * y_val

            if (datum%data_type(1:3) == 'tt:' .or. datum%data_type(1:2) == 't:') then
              if (datum%ix_ele > datum%ix_ele2) datum%ix_ele2 = datum%ix_ele
            endif

          enddo
        enddo

!----------------------------------------------------------------------------
! Bad data_source

      case default
        call out_io (s_error$, r_name, 'UNKNOWN DATA_SOURCE: ' // curve%data_source)
        graph%valid = .false.
        cycle graph_loop
      end select

!----------------------------------------------------------------------------
! Calculate the points for drawing the curve through the symbols.
! If the x-axis is by index or ele_index then these points are the same as the symbol points.
! That is, for x-axis = index or ele_index the line is piece-wise linear between the symbols.
! If the axis is by s-value then the line is a "smooth" curve with n_curve_pts points if
! plotting model, base or design data. It's the same as the symbol points otherwise.
! Smoothing will only be performed if performing single particle tracking.

      if (plot%x_axis_type == 'index') then
        call reassociate_real (curve%y_line, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_line, n_dat) ! allocate space for the data
        curve%x_line = curve%x_symb
        curve%y_line = curve%y_symb
      elseif (plot%x_axis_type == 'ele_index') then
        call reassociate_real (curve%y_line, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_line, n_dat) ! allocate space for the data
        curve%x_line = curve%x_symb
        curve%y_line = curve%y_symb
      elseif (plot%x_axis_type == 's') then
        smooth_curve = .true.
        do m = 1, size(graph%who)
          if (graph%who(m)%name .eq. 'meas' .or. graph%who(m)%name .eq. 'ref' .or. &
	      s%global%track_type .ne. 'single') smooth_curve = .false.
        enddo
        if (smooth_curve) then
          ! allocate data space
          call reassociate_real (curve%y_line, s%global%n_curve_pts) 
          call reassociate_real (curve%x_line, s%global%n_curve_pts) 
          curve%y_line = 0
          do m = 1, size(graph%who)
            select case (graph%who(m)%name)
            case (' ') 
              cycle
            case ('model')
              call s_data_to_plot (u%model, u%model_orb, curve%data_type, &
                                                       graph%who(m), curve, err)
              if (err) cycle graph_loop
            case ('base')  
              call s_data_to_plot (u%base, u%base_orb, curve%data_type, &
                                                       graph%who(m), curve, err)
              if (err) cycle graph_loop
            case ('design')  
              call s_data_to_plot (u%design, u%design_orb, curve%data_type, &
                                                       graph%who(m), curve, err)
              if (err) cycle graph_loop
            case default
              call out_io (s_error$, r_name, &
                       'BAD PLOT "WHO" WITH "S" X-AXIS: ' // graph%who(m)%name)
              graph%valid = .false.
              cycle graph_loop
            end select
          enddo
        else
          ! allocate space for the data
          call reassociate_real (curve%y_line, n_dat) 
          call reassociate_real (curve%x_line, n_dat) 
          curve%x_line = curve%x_symb 
          curve%y_line = curve%y_symb 
        endif
        
      endif

!----------------------------------------------------------------------------
! Renormalize and check for limited graph
! Note: Since there is an arbitrary overall phase, the phase data 
! gets renormalized so that the average value is zero.

      curve%y_symb = curve%y_symb * curve%units_factor
      curve%y_line = curve%y_line * curve%units_factor

      if (curve%data_type(1:6) == 'phase:' .and. n_dat /= 0 .and. curve%ele2_name == ' ') then
        f = sum(curve%y_symb) / n_dat
        curve%y_symb = curve%y_symb - f
        curve%y_line = curve%y_line - f 
      endif 

      curve%limited = &
              any(curve%y_symb .lt. graph%y%min .or. curve%y_symb .gt. graph%y%max)

! For the title_suffix: strip off leading "+" and enclose in "[ ]".

      graph%title_suffix = ' '
      do m = 1, size(graph%who)
        if (graph%who(m)%name == ' ') cycle
        if (graph%who(m)%sign == 1) then
          graph%title_suffix = trim(graph%title_suffix) // ' + ' // graph%who(m)%name
        elseif (graph%who(m)%sign == -1) then
          graph%title_suffix = trim(graph%title_suffix) // ' - ' // graph%who(m)%name
        endif
      enddo

      if (graph%title_suffix(2:2) == '+') graph%title_suffix = graph%title_suffix(4:)

      graph%title_suffix = '[' // trim(graph%title_suffix) // ']'

! attach x-axis type to title suffix 	 
  	 
       if (plot%x_axis_type .eq. 'index') then 	 
         graph%title_suffix = trim(graph%title_suffix) // ', X-axis: index, ' 	 
       elseif (plot%x_axis_type .eq. 'ele_index') then 	 
         graph%title_suffix = trim(graph%title_suffix) // ', X-axis: ele_index, ' 	 
       elseif (plot%x_axis_type .eq. 's') then 	 
         graph%title_suffix = trim(graph%title_suffix) // ', X-axis: s, ' 	 
       endif 	 

    enddo

  enddo graph_loop
enddo plot_loop


!----------------------------------------------------------------------------
contains

subroutine s_data_to_plot (lat, orb, data_type, who, curve, err)

type (tao_plot_who_struct) who
type (coord_struct) orb(0:)
type (ring_struct) lat
type (ele_struct) ele
type (coord_struct) here
type (tao_curve_struct) curve
type (taylor_struct) t_map(6)

real(rp) cbar(2,2), s_last, s_now, value, mat6(6,6)
integer i, j, k, expnt(6)
character(16) data_type
logical err

!

err = .false.

x1 = min (u%model%ele_(u%model%n_ele_use)%s, max (plot%x%min, u%model%ele_(0)%s))
x2 = min (u%model%ele_(u%model%n_ele_use)%s, max (plot%x%max, u%model%ele_(0)%s))
s_last = lat%ele_(curve%ix_ele2)%s

if (data_type(1:2) == 'r:') data_type = 'r:'
if (data_type(1:2) == 't:') data_type = 't:'
if (data_type(1:3) == 'tt:') data_type = 'tt:'

!

do ii = 1, size(curve%x_line)

  s_now = x1 + (ii-1) * (x2-x1) / (size(curve%x_line)-1)
  curve%x_line(ii) = s_now
  value = 0

  call twiss_and_track_at_s (lat, curve%x_line(ii), ele, orb, here)

  select case (data_type)
  case ('orbit:x')
    value = here%vec(1)
  case ('orbit:y')
    value = here%vec(3)
  case ('orbit:z')
    value = here%vec(5)
  case ('orbit:p_x')
    value = here%vec(2)
  case ('orbit:p_y')
    value = here%vec(4)
  case ('orbit:p_z')
    value = here%vec(6)
  case ('phase:x')
    value = ele%x%phi
  case ('phase:y')
    value = ele%y%phi
  case ('beta:x')
    value = ele%x%beta
  case ('beta:y')
    value = ele%y%beta
  case ('alpha:x')
    value = ele%x%alpha
  case ('alpha:y')
    value = ele%y%alpha
  case ('eta:x')
    value = ele%x%eta
  case ('eta:y')
    value = ele%y%eta
  case ('etap:x')
    value = ele%x%etap
  case ('etap:y')
    value = ele%y%etap
  case ('beam_energy')
    value = ele%value(beam_energy$)
  case ('cbar:11')
    call c_to_cbar (ele, cbar)
    value = cbar(1,1)
  case ('cbar:12')
    call c_to_cbar (ele, cbar)
    value = cbar(1,2)
  case ('cbar:21')
    call c_to_cbar (ele, cbar)
    value = cbar(2,1)
  case ('cbar:22')
    call c_to_cbar (ele, cbar)
    value = cbar(2,2)
  case ('r:')
    if (ii == 1) call mat_make_unit (mat6)
    if (s_now < s_last) cycle
    i = read_this_index (data_type, 3); if (i == 0) return
    j = read_this_index (data_type, 4); if (j == 0) return
    call tao_mat6_calc_at_s (lat, mat6, s_last, s_now, unit_start = .false.)
    value = mat6(i, j)
  case ('t:')
    if (ii == 1) call taylor_make_unit (t_map)
    if (s_now < s_last) cycle
    i = read_this_index (data_type, 3); if (i == 0) return
    j = read_this_index (data_type, 4); if (j == 0) return
    k = read_this_index (data_type, 5); if (k == 0) return
    call tao_transfer_map_calc_at_s (lat, t_map, s_last, s_now, unit_start = .false.)
    value = taylor_coef (t_map(i), j, k)
  case ('tt:')
    if (ii == 1) call taylor_make_unit (t_map)
    if (s_now < s_last) cycle
    expnt = 0
    i = read_this_index (data_type, 4); if (i == 0) return
    do j = 5, 15
      if (data_type(j:j) == ' ') exit
      k = read_this_index (data_type, j); if (k == 0) return
      expnt(k) = expnt(k) + 1
    enddo
    call tao_transfer_map_calc_at_s (lat, t_map, s_last, s_now, unit_start = .false.)
    value = taylor_coef (t_map(i), expnt)
  case default
   call out_io (s_fatal$, r_name, &
                  'DO NOT KNOW ABOUT THIS DATA_TYPE: ' // curve%data_type)
   call out_io (s_blank$, r_name, "Will not perfrom any plot smoothing")
    err = .true.
    return
  end select

  curve%y_line(ii) = curve%y_line(ii) + who%sign * value
  s_last = s_now

enddo

end subroutine

end subroutine
