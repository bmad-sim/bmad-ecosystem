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

integer i, j, k, m, n_dat, i_uni
logical err

character(20) :: r_name = 'tao_plot_data_setup'
character(40) class

!

plot_loop: do i = 1, size(s%plot_page%plot)

  plot => s%plot_page%plot(i)
  plot%valid = .true.   ! assume everything OK

  if (.not. plot%visible) cycle  ! Don't worry about invisable graphs

  do j = 1, size(plot%graph)

    graph => plot%graph(j)

    do k = 1, size(graph%curve)

      curve => graph%curve(k)

      i_uni = s%global%u_view  ! universe where the data comes from
      if (curve%ix_universe /= 0) i_uni = curve%ix_universe
      class = trim(plot%class) // ':' // curve%sub_class
      call tao_find_data (err, s%u(i_uni), class, d2_ptr, d1_ptr)
      if (err) cycle
      call tao_useit_plot_calc (plot, d1_ptr%d) ! make sure %useit_plot up-to-date
      n_dat = count (d1_ptr%d%useit_plot)       ! count the number of data points
      call reallocate_real (curve%y_dat, n_dat) ! allocate space for the data
      call reallocate_real (curve%x_dat, n_dat) ! allocate space for the data
      call reallocate_integer (curve%ix_data, n_dat)
      curve%x_dat = pack(d1_ptr%d%ix_d, mask = d1_ptr%d%useit_plot)
      curve%ix_data = pack(d1_ptr%d%ix_d, mask = d1_ptr%d%useit_plot)
      graph%title_suffix = ' '

! calculate the y-axis data point values.

      curve%y_dat = 0
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

      curve%y_dat = curve%y_dat * curve%units_factor
      if (plot%convert) curve%y_dat = curve%y_dat * &
                         pack(d1_ptr%d%conversion_factor, d1_ptr%d%useit_plot)
      if (plot%class == 'phase' .and. n_dat /= 0) curve%y_dat = &
                                    curve%y_dat - sum(curve%y_dat) / n_dat

! For the title_suffix: strip off leading "+" and enclose in "[ ]".

      if (graph%title_suffix(2:2) == '+') graph%title_suffix = graph%title_suffix(4:)
      graph%title_suffix = '[' // trim(graph%title_suffix) // ']'

    enddo
  enddo
enddo plot_loop


!----------------------------------------------------------------------------
contains

subroutine data_to_plot (value, who)

real(rp) value(:)
type (tao_plot_who_struct) who
character(16) name
character(80) line

!

curve%y_dat = curve%y_dat + who%sign * pack(value, mask = d1_ptr%d%useit_plot)  

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

end subroutine
