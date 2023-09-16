program compare_tracking_methods_plot

use bmad
use quick_plot

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (all_pointer_struct), allocatable :: ptr_array(:)

character(40) :: input_file  = 'compare_tracking_methods_plot.bmad'
character(40) :: output_file = 'compare_tracking_methods_plot'
character(20) :: base_method, element_to_vary, attrib_to_vary
character(20) :: veto_methods(10)

real(rp) :: scan_start, scan_end, step_size
real(rp), dimension(:), allocatable   :: scan_var
integer :: nsteps, nmethods, i, j, k, l, ix_base, id, xbox, ybox
logical :: is_valid(n_methods$+1)
logical :: err_flag

type mystep           ! stores the phase space coordinates for a particular iteration of the scan
   real(rp) :: vec(6) ! (x, px, y, py, z, pz)
end type mystep

type mymethod                            ! stores the relevant information for a particular method
   type (mystep), allocatable :: step(:)
   character(18) :: method_name
   character(4)  :: short_method_name
end type mymethod

type (mymethod), allocatable :: method(:), dmethod(:) 
type (coord_struct) end_orb
real(rp) :: ymin(6), ymax(6)
character(200) :: start_orb_desc, x_axis_desc, y_axis_desc
character(1) ans
type (qp_line_struct), dimension(:), allocatable :: lines
character(4) :: short_method_name(16)
character(10), parameter :: y_axis_label(6) = [ "\gDx      ", "\gDp\dx\u ", "\gDy      ", &
                                                "\gDp\dy\u ", "\gDz      ", "\gDp\dz\u " ]

namelist / scan_params / output_file, base_method, veto_methods, element_to_vary, attrib_to_vary, scan_start, scan_end, nsteps

!

short_method_name(bmad_standard$)     = 'STD'
short_method_name(symp_lie_ptc$)      = 'SLP'
short_method_name(runge_kutta$)       = 'RK'
short_method_name(linear$)            = 'LIN'
short_method_name(taylor$)            = 'TAY'
short_method_name(symp_lie_bmad$)     = 'SLB'
short_method_name(mad$)               = 'MAD'
short_method_name(time_runge_kutta$)  = 'TRK'

if (command_argument_count() == 1) call get_command_argument (1, input_file)

veto_methods = ''

call bmad_parser (input_file, lat)
ele => lat%ele(1)

open (1, file = input_file)
read (1, nml = scan_params)
close (1)

start_orb_desc = "Starting orbit: (" // trim(adjustl(convert_to_string(lat%particle_start%vec(1)))) // "," // trim(adjustl(convert_to_string(lat%particle_start%vec(2)))) // "," // trim(adjustl(convert_to_string(lat%particle_start%vec(3)))) // "," // trim(adjustl(convert_to_string(lat%particle_start%vec(4)))) // "," // trim(adjustl(convert_to_string(lat%particle_start%vec(5)))) // "," // trim(adjustl(convert_to_string(lat%particle_start%vec(6)))) // ")"

! Make a list of methods to use

do i = 1, size(is_valid)
  if (i /= n_methods$+1) then 
     is_valid(i) = valid_tracking_method(ele, electron$, i)
  else
     is_valid(i) = valid_tracking_method(ele, electron$, symp_lie_ptc$)
  end if
enddo

DO i = 1, size(veto_methods)
   call match_word (veto_methods(i), tracking_method_name(1:), j)
   if(j > 0) is_valid(j) = .false.
   if(j == symp_lie_ptc$) is_valid(n_methods$+1) = .false.
END DO

is_valid(custom$) = .false.

!

nmethods = count(is_valid)

allocate (scan_var(nsteps+1), method(nmethods), dmethod(nmethods-1), lines(nmethods-1))
do i = 1, nmethods
   allocate (method(i)%step(nsteps+1))
   if(.not. (i == nmethods)) allocate (dmethod(i)%step(nsteps+1))
end do

open (1, file = trim(adjustl(output_file)) // '.abs')
write (1,'(a)',advance='no') "Compare tracking methods as "; write (1,'(a)',advance='no') trim(adjustl(element_to_vary)); write (1,'(a)',advance='no') "%"; write (1,'(a)',advance='no') trim(adjustl(attrib_to_vary))
write (1,'(a)',advance='no') " is varied. Look at phase space coordinates at exit for lattice with only a "; write (1,'(a)',advance='no') trim(adjustl(key_name(ele%key))); write (1,'(a)',advance='no') "."
write (1,*); write (1,*)

open (2, file = trim(adjustl(output_file)) // '.diff')
write (2,'(a)',advance='no') "Compare tracking methods as "; write (2,'(a)',advance='no') trim(adjustl(element_to_vary)); write (2,'(a)',advance='no') "%"; write (2,'(a)',advance='no') trim(adjustl(attrib_to_vary))
write (2,'(a)',advance='no') " is varied. Look at difference in phase space coordinates at exit for lattice with only a "; write (2,'(a)',advance='no') trim(adjustl(key_name(ele%key))); write (2,'(a)',advance='no') "."
write (2,*); write (2,*)

step_size = (scan_end - scan_start) / nsteps

! Figure out which method is the base method

k = 0
DO i = 1, size(is_valid)
   if (.not. is_valid(i)) cycle
   k = k+1
   if (i /= n_methods$+1) then
      method(k)%method_name = tracking_method_name(i)
      method(k)%short_method_name = short_method_name(i)
   else   
      method(k)%method_name = 'Symp_Lie_PTC_Exact'
      method(k)%short_method_name = 'SLPE'
   end if
enddo

call match_word (base_method, method%method_name, ix_base, can_abbreviate = .false.)
if (ix_base == 0) then 
   print *, "Base method is not valid"; stop
end if 

! For each method, perform scan and store phase space coordinates at exit

k = 0
DO i = 1, size(is_valid)
   if (.not. is_valid(i)) cycle
   k = k+1
   if (i /= n_methods$+1) then
      ele%tracking_method = i  
      write (1,*) "Tracking Method = ", tracking_method_name(i)
   else 
      ptc_com%exact_model = .true.
      ele%tracking_method = symp_lie_ptc$
      write (1,*) "Tracking Method = Symp_Lie_PTC_Exact"
   end if
   DO j = 1, nsteps+1
      if (k == 1) scan_var(j) = scan_start + (j-1) * step_size 
      call pointers_to_attribute (lat, element_to_vary, attrib_to_vary, .false., ptr_array, err_flag)
      ptr_array(1)%r = scan_start + (j-1) * step_size    
      if (j == 1 .and. i == linear$) then
        ele%tracking_method = ix_base
        call make_mat6 (ele, lat%param, lat%particle_start)
        ele%tracking_method = i
      endif
      call init_coord (lat%particle_start, lat%particle_start, ele = ele, element_end = upstream_end$)
      call track1 (lat%particle_start, ele, lat%param, end_orb)
      method(k)%step(j)%vec = end_orb%vec
      write (1,'(es24.15,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)',advance='no') ptr_array(1)%r, end_orb%vec
      write (1,*)
   END DO
   write (1,*)
END DO

close(1)

! For each method, calculate difference in phase space coordinates w.r.t. base method   
k = 0
DO i = 1, nmethods
   if (i == ix_base) cycle
   k = k+1
   write (2,*) "Tracking Method = ", method(i)%method_name
   dmethod(k)%method_name = method(i)%method_name
   dmethod(k)%short_method_name = method(i)%short_method_name
   DO j = 1, nsteps+1
      dmethod(k)%step(j)%vec  = method(i)%step(j)%vec - method(ix_base)%step(j)%vec
      write (2,'(es24.15,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)',advance='no') scan_start + (j-1) * step_size, dmethod(k)%step(j)%vec
      write (2,*)
   END DO
   write (2,*)
END DO

close(2)

ymin = 0.; ymax = 0.

! Figure out the y axis range for various plots
DO i = 1, nmethods -1
   DO l = 1, 6
      ymin(l) = min(ymin(l),minval(dmethod(i)%step(:)%vec(l)))
      ymax(l) = max(ymax(l),maxval(dmethod(i)%step(:)%vec(l)))
   END DO
END DO

call qp_open_page ("PS-L", id, 1000.0_rp, 770.0_rp, "POINTS")
call create_plot
call qp_close_page

call qp_open_page ("X", id, 1000.0_rp, 770.0_rp, "POINTS")
call create_plot
write (*, "(a)", advance = "NO") " Hit any class to end program: "
read (*, "(a)") ans

deallocate (scan_var, method, dmethod)

contains
character(8) function convert_to_string(a)
  real(rp) :: a
  write(convert_to_string, "(F8.3)") a
end function convert_to_string

subroutine which_box(num,xbox,ybox)
  integer :: num, xbox, ybox
  SELECT CASE (num)
   CASE (1)
      xbox = 1; ybox = 3
   CASE (2)
      xbox = 2; ybox = 3
   CASE (3)
      xbox = 1; ybox = 2
   CASE (4)
      xbox = 2; ybox = 2
   CASE (5)
      xbox = 1; ybox = 1
   CASE (6)
      xbox = 2; ybox = 1
   END SELECT
end subroutine which_box

subroutine my_draw_axes(xbox,ybox,y_min,y_max,ylabel)
  integer :: xbox, ybox
  real(rp) :: y_min,y_max
  character(*) :: ylabel
  call qp_set_box (xbox, ybox, 2, 3)
  call qp_calc_and_set_axis ("Y", y_min, y_max, 4, 8, "GENERAL")
  call qp_draw_axes (y_lab=ylabel)
end subroutine my_draw_axes

subroutine plot_it(xbox,ybox,y_min,y_max,xvec,yvec)
  integer :: xbox, ybox
  real(rp) :: y_min,y_max
  real(rp), dimension(:) :: xvec, yvec
  call qp_set_box (xbox, ybox, 2, 3)
  call qp_calc_and_set_axis ("Y", y_min, y_max, 4, 8, "GENERAL")
  call qp_draw_data (xvec, yvec, symbol_every = 0)
end subroutine

subroutine create_plot
  call qp_set_page_border (0.02_rp, 0.01_rp, 0.01_rp, 0.1_rp, "%PAGE")
  call qp_set_margin (0.07_rp, 0.05_rp, 0.02_rp, 0.01_rp, "%PAGE")
  call qp_calc_and_set_axis ("X", scan_start, scan_end, 4, 8, "GENERAL")

  ! Draw axes
  DO i = 1, 6
     call which_box(i, xbox, ybox)
     call my_draw_axes (xbox, ybox, ymin(i), ymax(i), y_axis_label(i))
  END DO

  ! Draw curves
  DO i = 1, nmethods-1 
     call qp_set_line_attrib ("PLOT", color = qp_color_name(i), pattern = 'solid')
     lines(i)%color = qp_color_name(i)
     DO j = 1, 6
        call which_box(j, xbox, ybox)
        call plot_it (xbox, ybox, ymin(j), ymax(j), scan_var, dmethod(i)%step(:)%vec(j))
        k = max(1,int(nsteps*(0.1+i*0.09)))
        call qp_draw_text (dmethod(i)%short_method_name, scan_var(k), dmethod(i)%step(k)%vec(j), height = 10.0_rp, color = qp_color_name(i)) 
     END DO
  END DO

  ! Draw legend
  if (nmethods > 6) then
     if (mod(nmethods+1,2) == 1) then
        call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines(1:(nmethods+1)/2), text = dmethod(1:(nmethods+1)/2)%method_name, draw_symbol = .false.)
        call qp_draw_curve_legend (0.28_rp, 0.99_rp, "%PAGE", line = lines((nmethods+3)/2:(nmethods-1)), text = dmethod((nmethods+3)/2:(nmethods-1))%method_name, draw_symbol = .false.)
     else
        call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines(1:nmethods/2), text = dmethod(1:nmethods/2)%method_name, draw_symbol = .false.)
        call qp_draw_curve_legend (0.28_rp, 0.99_rp, "%PAGE", line = lines(nmethods/2+1:(nmethods-1)), text = dmethod(nmethods/2+1:(nmethods-1))%method_name, draw_symbol = .false.)
     end if
  else
     call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines, text = dmethod(:)%method_name, draw_symbol = .false.)   
  end if

  ! Write description on plot
  call qp_draw_text (start_orb_desc, 0.53_rp, 0.98_rp, "%PAGE", height = 15.0_rp) 
  x_axis_desc = "X axis: parameter being varied (" // trim(adjustl(element_to_vary)) // "%" // trim(adjustl(attrib_to_vary)) // ")."
  call qp_draw_text (x_axis_desc, 0.53_rp, 0.95_rp, "%PAGE", height = 15.0_rp)
  y_axis_desc = "Y axis: absolute difference w.r.t. " // trim(adjustl(method(ix_base)%method_name)) // " tracking."
  call qp_draw_text (y_axis_desc, 0.53_rp, 0.92_rp, "%PAGE", height = 15.0_rp)
end subroutine

end program
