!+
! Subroutine tao_plot_data_setup (s)
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_page%plot(:)%graph(:)%curve(:) arrays.
!
! Input/Output:
!   s     -- Tao_super_universe_struct
!-

subroutine tao_plot_data_setup (s)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr

real(rp) f
integer i, ii, j, k, m, n_dat, i_uni
logical err

character(20) :: r_name = 'tao_plot_data_setup'

!

plot_loop: do i = 1, size(s%plot_page%plot)

  plot => s%plot_page%plot(i)
  plot%valid = .true.   ! assume everything OK

  if (.not. plot%visible) cycle  ! Don't worry about invisable graphs
  if (plot%type /= 'data') cycle ! Don't worry about other types

  do j = 1, size(plot%graph)

    graph => plot%graph(j)

    do k = 1, size(graph%curve)

      curve => graph%curve(k)

      i_uni = s%global%u_view  ! universe where the data comes from
      if (curve%ix_universe /= 0) i_uni = curve%ix_universe 
      call tao_find_data (err, s%u(i_uni), curve%data_class, d2_ptr, d1_ptr)
      if (err) cycle
      call tao_useit_plot_calc (plot, d1_ptr%d) ! make sure %useit_plot up-to-date
      n_dat = count (d1_ptr%d%useit_plot)       ! count the number of data points
      call reallocate_real (curve%y_symb, n_dat) ! allocate space for the data
      call reallocate_real (curve%x_symb, n_dat) ! allocate space for the data
      call reallocate_integer (curve%ix_symb, n_dat)
      curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)

      graph%title_suffix = ' '

! calculate the y-axis data point values.

      curve%y_symb = 0
      do m = 1, size(plot%who)
        select case (plot%who(m)%name)
        case (' ') 
          cycle
        case ('model')   
          call data_to_plot (d1_ptr%d%model_value, plot%who(m))
        case ('base')  
          call data_to_plot (d1_ptr%d%base_value, plot%who(m))
        case ('design')  
          call data_to_plot (d1_ptr%d%design_value, plot%who(m))
        case ('ref')     
          call data_to_plot (d1_ptr%d%ref_value, plot%who(m))
        case ('data')    
          call data_to_plot (d1_ptr%d%data_value, plot%who(m))
        case default
          call out_io (s_error$, r_name, 'BAD PLOT "WHO": ' // plot%who(m)%name)
          plot%valid = .false.
          cycle plot_loop
        end select
      enddo

! Note: Since there is an arbitrary overall phase, the phase data 
! gets renormalized so that the average value is zero.

      curve%y_symb = curve%y_symb * curve%units_factor
      if (plot%convert) curve%y_symb = curve%y_symb * &
                         pack(d1_ptr%d%conversion_factor, d1_ptr%d%useit_plot)
      if (curve%data_class(1:6) == 'phase:' .and. n_dat /= 0) then
         f = sum(curve%y_symb) / n_dat
         curve%y_symb = curve%y_symb - f
      endif 

! Calculate the line points.
! If the x-axis is by index then the line points are the same as the symbol points.
!  That is the line is piece-wise linear through the data points.
! If the axis is by s-value then the line is a smooth curve

      if (plot%x_axis_type == 'index') then
        curve%x_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)
        call reallocate_real (curve%y_line, n_dat) ! allocate space for the data
        call reallocate_real (curve%x_line, n_dat) ! allocate space for the data
        curve%x_line = curve%x_symb
        curve%y_line = curve%y_symb
      elseif (plot%x_axis_type == 's') then
        curve%x_symb = s%u(i_uni)%model%ele_(d1_ptr%d(curve%ix_symb)%ix_ele)%s
        call reallocate_real (curve%y_line, 400) ! allocate space for the data
        call reallocate_real (curve%x_line, 400) ! allocate space for the data
        curve%y_line = 0
        do ii = 1, size(curve%x_line)
          curve%x_line(ii) = plot%x%min + &
                    (ii-1) * (plot%x%max-plot%x%min) / (size(curve%x_line)-1)
          select case (plot%who(m)%name)
          case (' ') 
            cycle
          case ('model')
            call s_data_to_plot (s%u(i_uni)%model, s%u(i_uni)%model_orb,  & 
                          curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
            if (err) cycle plot_loop
          case ('base')  
            call s_data_to_plot (s%u(i_uni)%base, s%u(i_uni)%base_orb, &
                          curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
            if (err) cycle plot_loop
          case ('design')  
            call s_data_to_plot (s%u(i_uni)%design, s%u(i_uni)%design_orb, &
                          curve%x_line(ii), plot%who(m), curve%y_line(ii), err)
            if (err) cycle plot_loop
          case default
              call out_io (s_error$, r_name, &
                              'BAD PLOT "WHO" WITH "S" X-AXIS: ' // plot%who(m)%name)
              plot%valid = .false.
              cycle plot_loop
            end select
        enddo
        
      else
        call out_io (s_abort$, r_name, 'BAD X_AXIS_TYPE: ' // plot%x_axis_type)
        call err_exit
      endif

! For the title_suffix: strip off leading "+" and enclose in "[ ]".

      if (graph%title_suffix(2:2) == '+') graph%title_suffix = graph%title_suffix(4:)
      graph%title_suffix = '[' // trim(graph%title_suffix) // ']'

    enddo
  enddo
enddo plot_loop


!----------------------------------------------------------------------------
contains

subroutine data_to_plot (value, who)

type (tao_plot_who_struct) who
real(rp) value(:)
character(16) name
character(80) line

!

curve%y_symb = curve%y_symb + who%sign * pack(value, mask = d1_ptr%d%useit_plot)  

name = who%name
call str_upcase (name(1:1), name(1:1))

if (who%sign == 1) then
  graph%title_suffix = trim(graph%title_suffix) // ' + ' // name
elseif (who%sign == -1) then
  graph%title_suffix = trim(graph%title_suffix) // ' - ' // name
else
  write (line, '(a, i4)') 'BAD PLOT WHO SIGN: ', who%sign
  call out_io (s_error$, r_name, line)
  call err_exit
endif

end subroutine

!----------------------------------------------------------------------------
! contains

subroutine s_data_to_plot (lat, orb, s_pos, who, y, err)

type (tao_plot_who_struct) who
type (coord_struct) orb(:)
type (ring_struct) lat
type (ele_struct) ele
type (coord_struct) here
real(rp) s_pos, y, cbar(2,2)
logical err

!

call twiss_and_track_at_s (lat, s_pos, ele, orb, here)

err = .false.

select case (curve%data_class)
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
case ('eta:x')
  y = y + who%sign * ele%x%eta
case ('eta:y')
  y = y + who%sign * ele%y%eta
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
  call out_io (s_fatal$, r_name, 'DO NOT KNOW ABOUT THIS DATA_CLASS: ' // curve%data_class)
  err = .true.
end select

end subroutine

end subroutine