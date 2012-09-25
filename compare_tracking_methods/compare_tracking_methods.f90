program compare_tracking_methods

use bmad
use quick_plot

implicit none

type (lat_struct), target :: lat
character(20) :: base_method, element_to_vary, attrib_to_vary
character(20) :: veto_methods(10)
real(rp) :: scan_start, scan_end, step_size
integer :: nsteps, nmethods, i, j, k, l, m, id, xbox, ybox
real(rp), dimension(:), allocatable   :: scan_var
type (real_pointer_struct), allocatable :: ptr_array(:)
logical :: err_flag
logical :: valid_tracking_method(49,16) ! x coordinate indexes the element, y coordinate indexes the tracking method

type mystep           ! stores the phase space coordinates for a particular iteration of the scan
   real(rp) :: vec(6) ! (x, px, y, py, z, pz)
end type mystep

type mymethod                            ! stores the relevant information for a particular method
   type (mystep), allocatable :: step(:)
   character(16) :: method_name
   character(3)  :: short_method_name
end type mymethod

type (mymethod), allocatable :: method(:), dmethod(:) 
type (coord_struct) end_orb
real(rp) :: ymin(6), ymax(6)
character(200) :: start_orb_desc, x_axis_desc, y_axis_desc
character(1) ans
type (qp_line_struct), dimension(:), allocatable :: lines
character(3), parameter :: short_method_name(16) = [ 'STD', 'SLP', 'RK', 'LIN', '', 'SM', '', 'TAY', '', 'SLB', '', 'BOR', '', 'MAD', 'TRK', '']
character(10), parameter :: y_axis_label(6) = [ "\gDx", "\gDp\dx\u", "\gDy", "\gDp\dy\u", "\gDz", "\gDp\dz\u" ]

veto_methods = ''

namelist / scan_params / base_method, veto_methods, element_to_vary, attrib_to_vary, scan_start, scan_end, nsteps

call bmad_parser ('compare_tracking_methods/compare_tracking_methods.bmad', lat)

open (1, file = 'compare_tracking_methods/compare_tracking_methods.bmad')
read (1, nml = scan_params)
close (1)

start_orb_desc = "Starting orbit: (" // trim(adjustl(convert_to_string(lat%beam_start%vec(1)))) // "," // trim(adjustl(convert_to_string(lat%beam_start%vec(2)))) // "," // trim(adjustl(convert_to_string(lat%beam_start%vec(3)))) // "," // trim(adjustl(convert_to_string(lat%beam_start%vec(4)))) // "," // trim(adjustl(convert_to_string(lat%beam_start%vec(5)))) // "," // trim(adjustl(convert_to_string(lat%beam_start%vec(6)))) // ")"

! Define valid tracking methods for different lattice elements 
valid_tracking_method = .false.
valid_tracking_method(ab_multipole$, [1, 2, 4, 6, 8]) = .true.
valid_tracking_method(beambeam$, [1, 4]) = .true.
valid_tracking_method(drift$, [1, 2, 3, 4, 6, 8, 12, 14, 15]) = .true.
valid_tracking_method(e_gun$, [3, 12, 15]) = .true.
valid_tracking_method(ecollimator$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(elseparator$, [1, 2, 3, 4, 6, 8, 12, 14, 15]) = .true.
valid_tracking_method(em_field$, [3, 12, 15]) = .true.
valid_tracking_method(hkicker$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(instrument$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(kicker$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(lcavity$, [1, 2, 3, 4, 6, 8, 10, 12, 15]) = .true.
valid_tracking_method(marker$, [1, 2, 4, 6, 8]) = .true.
valid_tracking_method(match$, [1, 8]) = .true.
valid_tracking_method(monitor$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(multipole$, [1, 2, 4, 6, 8]) = .true.
valid_tracking_method(octupole$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(patch$, [1, 8]) = .true.
valid_tracking_method(quadrupole$, [1, 2, 3, 4, 6, 8, 10, 12, 14, 15]) = .true.
valid_tracking_method(rbend$, [1, 2, 4, 6, 8, 14]) = .true.
valid_tracking_method(rcollimator$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
valid_tracking_method(rfcavity$, [1, 2, 3, 4, 6, 8, 10, 12, 14, 15]) = .true.
valid_tracking_method(sbend$, [1, 2, 4, 6, 8, 14]) = .true.
valid_tracking_method(sextupole$, [1, 2, 3, 4, 6, 8, 12, 14, 15]) = .true.
valid_tracking_method(solenoid$, [1, 2, 3, 4, 6, 8, 10, 12, 14, 15]) = .true.
valid_tracking_method(sol_quad$, [1, 2, 3, 4, 6, 8, 10, 12, 14, 15]) = .true.
valid_tracking_method(taylor$, [1, 4]) = .true.
valid_tracking_method(vkicker$, [1, 2, 3, 4, 6, 8, 12, 15]) = .true.
if(lat%ele(1)%key == wiggler$ .and. lat%ele(1)%sub_key == map_type$) valid_tracking_method(wiggler$, [1, 2, 3, 4, 6, 8, 10, 12, 15]) = .true.
if(lat%ele(1)%key == wiggler$ .and. lat%ele(1)%sub_key == periodic_type$) valid_tracking_method(wiggler$, [1, 2, 3, 4, 6, 8, 10]) = .true.

! Remove vetoed methods from list of valid methods 
DO i = 1, size(veto_methods)
   call match_word (veto_methods(i), calc_method_name(1:), j)
   if(.not. j == 0) valid_tracking_method(lat%ele(1)%key,j) = .false.
END DO

nmethods = count(valid_tracking_method(lat%ele(1)%key,:))

allocate (scan_var(nsteps+1), method(nmethods), dmethod(nmethods-1), lines(nmethods-1))
do i = 1, nmethods
   allocate (method(i)%step(nsteps+1))
   if(.not. (i == nmethods)) allocate (dmethod(i)%step(nsteps+1))
end do

open (1, file = 'compare_tracking_methods/compare_tracking_methods.out')
write (1,'(a)',advance='no') "Compare tracking methods as "; write (1,'(a)',advance='no') trim(adjustl(element_to_vary)); write (1,'(a)',advance='no') "%"; write (1,'(a)',advance='no') trim(adjustl(attrib_to_vary))
write (1,'(a)',advance='no') " is varied. Look at phase space coordinates at exit for lattice with only a "; write (1,'(a)',advance='no') trim(adjustl(key_name(lat%ele(1)%key))); write (1,'(a)',advance='no') "."
write (1,*); write (1,*)

step_size = (scan_end - scan_start) / nsteps

! For each method, perform scan and store phase space coordinates at exit
k = 0
DO i = 1, size(calc_method_name)-1
   if(.not. valid_tracking_method(lat%ele(1)%key,i)) cycle
   lat%ele(1)%tracking_method = i
   k = k+1
   method(k)%method_name = calc_method_name(i)
   method(k)%short_method_name = short_method_name(i)
   write (1,*) "Tracking Method = ", calc_method_name(i)
   DO j = 1, nsteps+1
      if (k == 1) scan_var(j) = scan_start + (j-1) * step_size 
      call pointers_to_attribute (lat, element_to_vary, attrib_to_vary, .false., ptr_array, err_flag)
      ptr_array(1)%r = scan_start + (j-1) * step_size    
      call init_coord (lat%beam_start, lat%beam_start, ele = lat%ele(1), at_exit_end = .false.)
      call track1 (lat%beam_start, lat%ele(1), lat%param, end_orb)
      DO l = 1, 6
         method(k)%step(j)%vec(l) = end_orb%vec(l)
      END DO
      write (1,'(es24.15,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)',advance='no') ptr_array(1)%r, end_orb%vec(1), end_orb%vec(2), end_orb%vec(3), end_orb%vec(4), end_orb%vec(5), end_orb%vec(6)
      write (1,*)
   END DO
   write (1,*)
END DO

close(1)

! Figure out which method is the base method
call match_word (base_method, method%method_name, m, can_abbreviate = .false.)
if (m == 0) then 
   print *, "Base method is not valid"; stop
end if 

! For each method, calculate difference in phase space coordinates w.r.t. base method   
k = 0
DO i = 1, nmethods
   if (i == m) cycle
   k = k+1
   dmethod(k)%method_name = method(i)%method_name
   dmethod(k)%short_method_name = method(i)%short_method_name
   DO j = 1, nsteps+1
      DO l = 1, 6
         dmethod(k)%step(j)%vec(l)  = method(i)%step(j)%vec(l) - method(m)%step(j)%vec(l)
      END DO
   END DO
END DO

ymin = 0.; ymax = 0.

! Figure out the y axis range for various plots
DO i = 1, nmethods -1
   DO l = 1, 6
      ymin(l) = min(ymin(l),minval(dmethod(i)%step(:)%vec(l)))
      ymax(l) = max(ymax(l),maxval(dmethod(i)%step(:)%vec(l)))
   END DO
END DO

call qp_open_page ("X", id, 1000.0_rp, 770.0_rp, "POINTS")

call qp_set_page_border (0.02_rp, 0.01_rp, 0.01_rp, 0.1_rp, "%PAGE")
call qp_set_margin (0.07_rp, 0.05_rp, 0.02_rp, 0.01_rp, "%PAGE")
call qp_calc_and_set_axis ("X", scan_start, scan_end, 4, 8, "ZERO_AT_END")

! Draw axes
DO i = 1, 6
   call which_box(i, xbox, ybox)
   call my_draw_axes (xbox, ybox, ymin(i), ymax(i), y_axis_label(i))
END DO  

! Draw curves
DO i = 1, nmethods-1 
   call qp_set_line_attrib ("PLOT", color = i, pattern = 1)
   lines(i)%color = i
   DO j = 1, 6
      call which_box(j, xbox, ybox)
      call plot_it (xbox, ybox, ymin(j), ymax(j), scan_var, dmethod(i)%step(:)%vec(j))
      call qp_draw_text (dmethod(i)%short_method_name, scan_var(int(nsteps*(0.1+i*0.09))), dmethod(i)%step(int(nsteps*(0.1+i*0.09)))%vec(j), height = 10.0_rp, color = i) 
   END DO
END DO

! Draw legend
if (nmethods > 6) then
   if (mod(nmethods+1,2)) then
      call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines(1:(nmethods+1)/2), text = dmethod(1:(nmethods+1)/2)%method_name, draw_symbol = .false.)
      call qp_draw_curve_legend (0.3_rp, 0.99_rp, "%PAGE", line = lines((nmethods+3)/2:(nmethods-1)), text = dmethod((nmethods+3)/2:(nmethods-1))%method_name, draw_symbol = .false.)
   else
      call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines(1:nmethods/2), text = dmethod(1:nmethods/2)%method_name, draw_symbol = .false.)
      call qp_draw_curve_legend (0.3_rp, 0.99_rp, "%PAGE", line = lines(nmethods/2+1:(nmethods-1)), text = dmethod(nmethods/2+1:(nmethods-1))%method_name, draw_symbol = .false.)
   end if
else
   call qp_draw_curve_legend (0.09_rp, 0.99_rp, "%PAGE", line = lines, text = dmethod(:)%method_name, draw_symbol = .false.)   
end if

! Write description on plot
call qp_draw_text (start_orb_desc, 0.53_rp, 0.98_rp, "%PAGE", height = 15.0_rp) 
x_axis_desc = "X axis: parameter being varied (" // trim(adjustl(element_to_vary)) // "%" // trim(adjustl(attrib_to_vary)) // ")."
call qp_draw_text (x_axis_desc, 0.53_rp, 0.95_rp, "%PAGE", height = 15.0_rp)
y_axis_desc = "Y axis: absolute difference w.r.t. " // trim(adjustl(method(m)%method_name)) // " tracking."
call qp_draw_text (y_axis_desc, 0.53_rp, 0.92_rp, "%PAGE", height = 15.0_rp)

write (*, "(a)", advance = "NO") " Hit any class to end program: "
accept "(a)", ans

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
  real(rp), dimension(100) :: xvec, yvec
  call qp_set_box (xbox, ybox, 2, 3)
  call qp_calc_and_set_axis ("Y", y_min, y_max, 4, 8, "GENERAL")
  call qp_draw_data (xvec, yvec, symbol_every = 0)
end subroutine

end program
