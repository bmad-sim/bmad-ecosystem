!+
! Subroutine tao_plot_data_setup ()
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_page%plot(:)%graph(:)%curve(:) arrays.
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

real(rp) f, x1, x2, y_val, eps
real(rp), pointer :: value(:)

integer i, ii, j, k, m, n_dat, i_uni, ie
logical err, smoothing

character(20) :: r_name = 'tao_plot_data_setup'

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
      do k = 1, size(s%plot_page%ele_shape(:))
        if (s%plot_page%ele_shape(k)%key == 0) exit
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

plot_loop: do i = 1, size(s%plot_page%plot)

  plot => s%plot_page%plot(i)
  plot%valid = .true.   ! assume everything OK

  if (.not. plot%visible) cycle  ! Don't worry about invisable graphs

  if (plot%x_axis_type /= 'index' .and. plot%x_axis_type /= 's') then
    call out_io (s_abort$, r_name, 'BAD X_AXIS_TYPE: ' // plot%x_axis_type)
    plot%valid = .false.
    cycle
  endif

! loop over all graphs and curves

  do j = 1, size(plot%graph)

    graph => plot%graph(j)
    if (graph%type /= 'data') cycle ! Don't worry about other types

    do k = 1, size(graph%curve)

      curve => graph%curve(k)

      if (curve%data_source == 'lat_layout' .and. &
                                       plot%x_axis_type == 'index') then
        call out_io (s_error$, r_name, 'CURVE%DATA_SOURCE = "LAT_LAYOUT" ' // &
                        'AND PLOT%X_AXIS_TYPE = "INDEX" DO NOT GO TOGETHER.')
        plot%valid = .false.
        cycle plot_loop
      endif

      i_uni = s%global%u_view  ! universe where the data comes from
      if (curve%ix_universe /= 0) i_uni = curve%ix_universe 
      u => s%u(i_uni)

!----------------------------------------------------------------------------
! data_source is a data array

      select case (curve%data_source)
      case ('data_array')
        call tao_find_data (err, u, curve%data_type, d2_ptr, d1_ptr)
        if (err) then
          plot%valid = .false.
          cycle plot_loop
        endif

        d1_ptr%d%good_plot = .true.
        eps = 1e-4 * (plot%x%max - plot%x%min)
        if (plot%x_axis_type == 'index') then
          where (d1_ptr%d%ix_d1 < plot%x%min-eps) d1_ptr%d%good_plot = .false.
          where (d1_ptr%d%ix_d1 > plot%x%max+eps) d1_ptr%d%good_plot = .false.
        else
          where (d1_ptr%d%s < plot%x%min-eps) d1_ptr%d%good_plot = .false.
          where (d1_ptr%d%s > plot%x%max+eps) d1_ptr%d%good_plot = .false.
        endif

        call tao_useit_data_plot_calc (plot, d1_ptr%d) ! make sure %useit_plot up-to-date
        n_dat = count (d1_ptr%d%useit_plot)       ! count the number of data points

        call reassociate_integer (curve%ix_symb, n_dat)
        call reassociate_real (curve%y_symb, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_symb, n_dat) ! allocate space for the data

        curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)

        if (plot%x_axis_type == 'index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 's') then
          curve%x_symb = u%model%ele_(d1_ptr%d(curve%ix_symb)%ix_ele)%s
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0

        do m = 1, size(plot%who)
          select case (plot%who(m)%name)
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
            call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // plot%who(m)%name)
            plot%valid = .false.
            cycle plot_loop
          end select
          curve%y_symb = curve%y_symb + &
                   plot%who(m)%sign * pack(value, mask = d1_ptr%d%useit_plot)
        enddo

        if (curve%convert) curve%y_symb = curve%y_symb * &
                           pack(d1_ptr%d%conversion_factor, d1_ptr%d%useit_plot)


!----------------------------------------------------------------------------
! data_source is a var array

      case ('var_array')
        call tao_find_var (err, curve%data_type, v1_ptr)
        if (err) then
          plot%valid = .false.
          cycle plot_loop
        endif

        v1_ptr%v%good_plot = .true.
        eps = 1e-4 * (plot%x%max - plot%x%min)
        if (plot%x_axis_type == 'index') then
          where (v1_ptr%v%ix_v1 < plot%x%min-eps) v1_ptr%v%good_plot = .false.
          where (v1_ptr%v%ix_v1 > plot%x%max+eps) v1_ptr%v%good_plot = .false.
        else
          where (v1_ptr%v%s < plot%x%min-eps) v1_ptr%v%good_plot = .false.
          where (v1_ptr%v%s > plot%x%max+eps) v1_ptr%v%good_plot = .false.
        endif

        call tao_useit_var_plot_calc (plot, v1_ptr%v) ! make sure %useit_plot up-to-date
        n_dat = count (v1_ptr%v%useit_plot)       ! count the number of data points

        call reassociate_integer (curve%ix_symb, n_dat)
        call reassociate_real (curve%y_symb, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_symb, n_dat) ! allocate space for the data

        curve%ix_symb = pack(v1_ptr%v%ix_v1, mask = v1_ptr%v%useit_plot)

        if (plot%x_axis_type == 'index') then
          curve%x_symb = curve%ix_symb
        elseif (plot%x_axis_type == 's') then
	  ! FIX ME!!! need to set up for different lattices in diff universes
          !curve%x_symb = u%model%ele_(v1_ptr%v(curve%ix_symb)%this(1)%ix_ele)%s
          plot%valid = .false.
          cycle plot_loop
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0

        do m = 1, size(plot%who)
          select case (plot%who(m)%name)
          case (' ') 
            cycle
          case ('model') 
            value => v1_ptr%v%model_value
          case ('base')  
            value => v1_ptr%v%base_value
          case ('design')  
            value => v1_ptr%v%design_value
          case ('ref')     
            value => v1_ptr%v%ref_value
          case ('meas')    
            value => v1_ptr%v%meas_value
          case default
            call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // plot%who(m)%name)
            plot%valid = .false.
            cycle plot_loop
          end select
          curve%y_symb = curve%y_symb + &
                   plot%who(m)%sign * pack(value, mask = v1_ptr%v%useit_plot)
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
        elseif (plot%x_axis_type == 's') then
          curve%x_symb = u%model%ele_(curve%ix_symb)%s
        endif

! calculate the y-axis data point values.

        curve%y_symb = 0
        datum%ix_ele2 = -1
        datum%merit_type = 'target'
        datum%data_type = curve%data_type

        do ie = 1, n_dat
          datum%ix_ele = curve%ix_symb(ie)
          do m = 1, size(plot%who)
            select case (plot%who(m)%name)
            case (' ') 
              cycle
            case ('model')   
              call tao_evaluate_a_datum (datum, u%model, u%model_orb, y_val)
            case ('base')  
              call tao_evaluate_a_datum (datum, u%base, u%base_orb, y_val)
            case ('design')  
              call tao_evaluate_a_datum (datum, u%design, u%design_orb, y_val)
            case default
              call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // plot%who(m)%name)
              plot%valid = .false.
              cycle plot_loop
            end select
            curve%y_symb(ie) = curve%y_symb(ie) + plot%who(m)%sign * y_val
          enddo
        enddo

!----------------------------------------------------------------------------
! Bad data_source

      case default
        call out_io (s_error$, r_name, 'UNKNOWN DATA_SOURCE: ' // curve%data_source)
        plot%valid = .false.
        cycle plot_loop
      end select


!----------------------------------------------------------------------------
! Note: Since there is an arbitrary overall phase, the phase data 
! gets renormalized so that the average value is zero.

      curve%y_symb = curve%y_symb * curve%units_factor

      if (curve%data_type(1:6) == 'phase:' .and. n_dat /= 0) then
         f = sum(curve%y_symb) / n_dat
         curve%y_symb = curve%y_symb - f
      endif 

! Calculate the points for drawing the curve through the symbols.
! If the x-axis is by index then these points are the same as the symbol points.
!  That is, for x-axis = index the line is piece-wise linear between the symbols.
! If the axis is by s-value then the line is a "smooth" curve with 400 points if
! plotting model, base or design data. It's the same as the symbol points
! otherwise.

      if (plot%x_axis_type == 'index') then
        call reassociate_real (curve%y_line, n_dat) ! allocate space for the data
        call reassociate_real (curve%x_line, n_dat) ! allocate space for the data
        curve%x_line = curve%x_symb
        curve%y_line = curve%y_symb
      elseif (plot%x_axis_type == 's') then
        smoothing = .true.
	do m = 1, size(plot%who)
	  if(plot%who(m)%name .eq. 'meas' .or. plot%who(m)%name .eq. 'ref') &
                     smoothing = .false.
        enddo
	if (smoothing) then
          call reassociate_real (curve%y_line, 400) ! allocate space for the data
          call reassociate_real (curve%x_line, 400) ! allocate space for the data
          curve%y_line = 0
          x1 = max (plot%x%min, u%model%ele_(0)%s)
          x2 = min (plot%x%max, u%model%ele_(u%model%n_ele_use)%s)
          do ii = 1, size(curve%x_line)
            curve%x_line(ii) = x1 + (ii-1) * (x2-x1) / (size(curve%x_line)-1)
            do m = 1, size(plot%who)
              select case (plot%who(m)%name)
              case (' ') 
                cycle
              case ('model')
                call s_data_to_plot (u%model, u%model_orb,  & 
                            curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
                if (err) cycle plot_loop
              case ('base')  
                call s_data_to_plot (u%base, u%base_orb, &
                            curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
                if (err) cycle plot_loop
              case ('design')  
                call s_data_to_plot (u%design, u%design_orb, &
                            curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
                if (err) cycle plot_loop
              case default
                call out_io (s_error$, r_name, &
                                'BAD PLOT "WHO" WITH "S" X-AXIS: ' // plot%who(m)%name)
                plot%valid = .false.
                cycle plot_loop
              end select
            enddo
          enddo
	else
          call reassociate_real (curve%y_line, n_dat) ! allocate space for the data
          call reassociate_real (curve%x_line, n_dat) ! allocate space for the data
          curve%x_line = curve%x_symb
          curve%y_line = curve%y_symb
	endif
        
      endif

! For the title_suffix: strip off leading "+" and enclose in "[ ]".

      graph%title_suffix = ' '
      do m = 1, size(plot%who)
        if (plot%who(m)%name == ' ') cycle
        if (plot%who(m)%sign == 1) then
          graph%title_suffix = trim(graph%title_suffix) // ' + ' // plot%who(m)%name
        elseif (plot%who(m)%sign == -1) then
          graph%title_suffix = trim(graph%title_suffix) // ' - ' // plot%who(m)%name
        endif
      enddo

      if (graph%title_suffix(2:2) == '+') graph%title_suffix = graph%title_suffix(4:)

      graph%title_suffix = '[' // trim(graph%title_suffix) // ']'

! attach x-axis type to title suffix

      if (plot%x_axis_type .eq. 'index') then
	graph%title_suffix = trim(graph%title_suffix) // ' index '
      elseif (plot%x_axis_type .eq. 's') then
	graph%title_suffix = trim(graph%title_suffix) // ' s '
      endif

    enddo
  enddo
enddo plot_loop


!----------------------------------------------------------------------------
contains

subroutine s_data_to_plot (lat, orb, s_pos, who, y, err)

type (tao_plot_who_struct) who
type (coord_struct) orb(0:)
type (ring_struct) lat
type (ele_struct) ele
type (coord_struct) here
real(rp) s_pos, y, cbar(2,2)
logical err

!

call twiss_and_track_at_s (lat, s_pos, ele, orb, here)

err = .false.

select case (curve%data_type)
case ('orbit:x')
  y = y + who%sign * here%vec(1)
case ('orbit:y')
  y = y + who%sign * here%vec(3)
case ('phase:x')
  y = y + who%sign * ele%x%phi
case ('phase:y')
  y = y + who%sign * ele%y%phi
case ('beta:x')
  y = y + who%sign * ele%x%beta
case ('beta:y')
  y = y + who%sign * ele%y%beta
case ('alpha:x')
  y = y + who%sign * ele%x%alpha
case ('alpha:y')
  y = y + who%sign * ele%y%alpha
case ('eta:x')
  y = y + who%sign * ele%x%eta
case ('eta:y')
  y = y + who%sign * ele%y%eta
case ('etap:x')
  y = y + who%sign * ele%x%etap
case ('etap:y')
  y = y + who%sign * ele%y%etap
case ('cbar:11')
  call c_to_cbar (ele, cbar)
  y = y + who%sign * cbar(1,1)
case ('cbar:12')
  call c_to_cbar (ele, cbar)
  y = y + who%sign * cbar(1,2)
case ('cbar:21')
  call c_to_cbar (ele, cbar)
  y = y + who%sign * cbar(2,1)
case ('cbar:22')
  call c_to_cbar (ele, cbar)
  y = y + who%sign * cbar(2,2)
case default
  call out_io (s_fatal$, r_name, 'DO NOT KNOW ABOUT THIS DATA_TYPE: ' // curve%data_type)
  call out_io (s_blank$, r_name, "Will not perfrom any plot smoothing")
  err = .true.
end select

end subroutine

end subroutine
