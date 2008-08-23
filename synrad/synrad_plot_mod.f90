module synrad_plot_mod

use cesrv_struct
use cesrv_interface
use sr_mod
use quick_plot

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! Integrates the power hit info for a window or target location

subroutine burn_plot (window, iw, ring, gen, logic)

implicit none

type (crotch_window_struct), target :: window(:)
type (logic_struct) :: logic
type (lat_struct) :: ring
type (synrad_param_struct) :: gen

integer  j, iw, m, n, lun, lunget
real(rp) :: x(window(iw)%n_ray_hit), y(window(iw)%n_ray_hit)
real(rp) :: sigma
type (coord_struct), pointer :: coord
type (ray_hit_struct), pointer :: ray_hit
character*16 line
character*60 title
logical target
real(rp) xmax, xmin, ymax, ymin
integer, parameter :: n_grid$ = 100
real(rp) :: grid(1:n_grid$,1:n_grid$), dx, dy, gen_factor
real(rp) :: power_factor, r, dpower, track_len, level(10)

print *, 'Do you want the burn plot at the projected target '
print '(a,$)', ' or at the window (default) ?  (Enter t or w) '
accept '(a)', line
call str_upcase(line, line)
line = adjustl(line)
target = .false.
if (line(1:1) == 'T') target = .true.
 
grid = 0

do j=1,window(iw)%n_ray_hit

  if (target) then
    coord => window(iw)%ray_hit_(j)%target_coord
  else
    coord => window(iw)%ray_hit_(j)%hit_coord
  endif
  x(j) = coord%vec(1) 
  y(j) = coord%vec(3) 

enddo

print '(a,$)', ' x min?  (Enter in meters, or enter 999 for auto) '
accept '(f)', xmin

type *, xmin
if (xmin == 999.0) then
  if (target) then
    xmax =  real((maxval(x)))
    xmin =  real((minval(x)))
  else
    xmax =  window(iw)%length
    xmin =  0.0
  endif
else
  print '(a,$)', ' x max?  (Enter in meters) '
  accept '(f)', xmax
endif

print '(a,$)', ' y min?  (Enter in meters, or enter 999 for auto) '
accept '(f)', ymin

if (ymin == 999.0) then
  ymax =  real((maxval(y)))
  ymin =  real((minval(y)))
else
  print '(a,$)', ' y max?  (Enter in meters) '
  accept '(f)', ymax
endif

dx = (xmax - xmin)/n_grid$
dy = (ymax - ymin)/n_grid$
gen_factor = 14.1e3 * gen%i_beam * ring%ele(0)%value(e_tot$)**4 
!  write(*,*) 'ymin, ymax, dy, xmin, xmax, dx'
!  write(*,*) ymin, ymax, dy, xmin, xmax, dx

do j=1,window(iw)%n_ray_hit

  ray_hit => window(iw)%ray_hit_(j)
  if (target) then
    sigma = ray_hit%sig_y_eff 
    track_len = ray_hit%dist + ray_hit%ray%track_len
    coord => window(iw)%ray_hit_(j)%target_coord
  else
    sigma = ray_hit%window_sig_y
    track_len = ray_hit%ray%track_len
    coord => window(iw)%ray_hit_(j)%hit_coord
  endif
  power_factor = gen_factor * ray_hit%ray%g_bend * dx * dy / &
                     (sqrt(twopi) * sigma * track_len)

  do m=1,n_grid$
    do n=1,n_grid$

      r = sqrt(((xmin+m*dx) - coord%vec(1))**2 + &
                ((ymin+n*dy) - coord%vec(3))**2)
      dpower = exp(-r**2/2/sigma**2) * power_factor
      grid(m,n) = grid(m,n) + dpower

    enddo
  enddo

enddo

lun = lunget()
open (lun, file = 'burn.out', status = 'new')
write(lun,*) grid

do m=1,n_grid$
  write(lun,*) m*dx+xmin
enddo
do m=1,n_grid$
  write(lun,*) m*dy+ymin
enddo

close(lun)

do m= 1,10

  level(m) =  (minval(grid) + maxval(grid)) / 11 * m

enddo


call qp_clear_page
call qp_set_box (1, 1, 1, 1)
call qp_set_margin (120.0_rp, 40.0_rp, 30.0_rp, 20.0_rp, 'POINTS')
call qp_set_page_border (0.0_rp, 0.0_rp, 0.0_rp, 20.0_rp, 'POINTS')
call qp_set_axis ('X', xmin, xmax, 10, 3) 
call qp_set_axis ('Y', ymin, ymax, 25, 4) 
level(1) =  (minval(grid) + maxval(grid)) / 2
title = "Burn Plot for " // logic%lattice
call plot_contour (grid, level, 'x', 'y', title)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! subroutine plot_contour (z_dat, levels, x_lab, y_lab, title)
!
! Subroutine to plot contour data, axes, and a title.   9/2001 mjf
!
! Input:
!     Z_DAT(:,:)         -- Real(Rp) matrix 
!     LEVELS(:)          -- Real(Rp) array of values for contours
!     X_LAB, Y_LAB       -- Character*(*) x and y axes labels.
!                             If = 'NULL' then label will not be printed.
!     TITLE              -- Character*(*) graph Title.
!-

subroutine plot_contour (z_dat, levels, x_lab, y_lab, title)

implicit none

real(rp) :: x_dat(2), y_dat(2)

real(rp) z_dat(:,:), levels(:), pt(1:4), currlevel
real(rp) pttop, ptright, ptbot, ptleft, ptmid
real(rp) rm, rn

character(*) x_lab, y_lab, title

integer :: ilevel, m, n, count, ipt, index(4), iptnext, i

! init

call qp_draw_axes (x_lab, y_lab, title)


rm = 1.0 / (size(z_dat,1) - 1)
rn = 1.0 / (size(z_dat,2) - 1)


do ilevel=1,size(levels)       ! Find the grid walls crossed for one contour
                           ! at a time
  currlevel = levels(ilevel)

  do m = 1, size(z_dat,1) - 1             ! loop through the vertical walls
                         
    do n = 1, size(z_dat,2) - 1           ! loop through the horizontal walls 

      count = 0
      pt(1) = z_dat(m,n)       !  Upper left point
      pt(2) = z_dat(m+1,n)     !  Upper right point
      pt(3) = z_dat(m+1,n+1)   !  Lower right point
      pt(4) = z_dat(m,n+1)     !  Lower left point

      do ipt=1,4

        if (pt(ipt) == currlevel) then

          count = count + 1
          index(count) = ipt

        endif

      enddo

      if (count == 4) cycle
      if (count == 3) then

        ! ------ Connect the points diagonally across from each other ----
        if ((pt(1) == currlevel) .and. (pt(3) == currlevel)) then

          ! ------ plot line between pt1 and pt3
          x_dat(1) = (m - 1) * rm

          x_dat(1) = (m - 1) * rm 
          y_dat(1) = (n - 1) * rn 
          x_dat(2) = (m) * rm 
          y_dat(2) = (n) * rn 
          call qp_draw_line (x_dat(1), x_dat(2), y_dat(1), y_dat(2), '%GRAPH')

        else

          ! ------ plot line between pt2 and pt4
          x_dat(1) = (m) * rm 
          y_dat(1) = (n - 1) * rn 
          x_dat(2) = (m - 1) * rm 
          y_dat(2) = (n) * rn 
          call qp_draw_line (x_dat(1), x_dat(2), y_dat(1), y_dat(2), '%GRAPH')

        endif
        cycle

      elseif (count == 2) then

        ! ---- plot line between points pt(index(1)) and pt(index(2))
        do i=1,2

          if (index(i) == 1) then

            x_dat(i) = (m - 1) * rm 
            y_dat(i) = (n - 1) * rn 

          elseif (index(i) == 2) then

            x_dat(i) = (m) * rm 
            y_dat(i) = (n - 1) * rn 

          elseif (index(i) == 3) then

            x_dat(i) = (m) * rm 
            y_dat(i) = (n) * rn 

          elseif (index(i) == 4) then

            x_dat(i) = (m - 1) * rm 
            y_dat(i) = (n) * rn

          endif

        enddo

        call qp_draw_line (x_dat(1), x_dat(2), y_dat(1), y_dat(2), '%GRAPH')
        cycle

      endif

      ! ------------ Now consider contour through edges ----------------
      count = 0
      index = 0

      do ipt=1,4

        if (ipt == 4) then
          iptnext = 1
        else
          iptnext = ipt +1
        endif

        if (pt(ipt) < pt(iptnext)) then

          if ((pt(ipt) < currlevel) .and. (currlevel <= pt(iptnext))) then
            count = count + 1
            index(count) = ipt
          endif

        elseif (pt(ipt) > pt(iptnext)) then

          if ((pt(ipt) > currlevel) .and. (currlevel >= pt(iptnext))) then
            count = count + 1
            index(count) = ipt
          endif

        endif

      enddo

      if (count == 1) cycle
      if (count == 2) then

        ! ---- plot line between points
        do i=1,2

          if (index(i) == 1) then
            x_dat(i) = (m + (currlevel - pt(1)) / (pt(2) - pt(1)) - 1) * rm 
            y_dat(i) = (n - 1) * rn 

          elseif (index(i) == 2) then

            x_dat(i) = m * rm 
            y_dat(i) = (n  + (currlevel - pt(2)) / (pt(3) - pt(2)) - 1) * rn

          elseif (index(i) == 3) then

            x_dat(i) = (m + (currlevel - pt(4)) / (pt(3) - pt(4)) - 1) * rm 
            y_dat(i) = n * rn 

          elseif (index(i) == 4) then

            x_dat(i) = (m - 1) * rm 
            y_dat(i) = (n  + (currlevel - pt(1)) / (pt(4) - pt(1)) - 1) * rn 

          endif

        enddo

        call qp_draw_line (x_dat(1), x_dat(2), y_dat(1), y_dat(2), '%GRAPH')
        cycle

      endif

      ! -- Divide the grid square into four squares and use mean values --
      ptmid   = (pt(1) + pt(2) + pt(3) + pt(4)) / 4
      pttop   = (pt(1) + pt(2)) / 2
      ptright = (pt(2) + pt(3)) / 2
      ptbot   = (pt(3) + pt(4)) / 2
      ptleft  = (pt(4) + pt(1)) / 2

    enddo

  enddo

enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Types out the crotch window hit information

subroutine ray_plot (window, iw)

implicit none

type (crotch_window_struct), target :: window(:)
type (coord_struct), pointer :: coord

integer  j, iw
real(rp) :: x(window(iw)%n_ray_hit), y(window(iw)%n_ray_hit)
character*16 line
logical target

!

print *, 'Do you want ray information at the projected target '
print '(a,$)', ' or at the window (default) ?  (Enter t or w) '
accept '(a)', line
call str_upcase(line, line)
line = adjustl(line)
target = .false.
if (line(1:1) == 'T') target = .true.
call qp_clear_page

do j=1,window(iw)%n_ray_hit

  if (target) then
    coord => window(iw)%ray_hit_(j)%target_coord
  else
    coord => window(iw)%ray_hit_(j)%hit_coord
  endif
  x(j) = coord%vec(1) * 1e3
  y(j) = coord%vec(3) * 1e6

enddo

call qp_set_box (1, 1, 1, 1)
call qp_set_margin (120.0_rp, 40.0_rp, 30.0_rp, 20.0_rp, 'POINTS')
call qp_set_page_border (0.0_rp, 0.0_rp, 0.0_rp, 20.0_rp, 'POINTS')

if (target) then
  call qp_set_axis ('X', real(floor(minval(x)), rp), real(ceiling(maxval(x)), rp), 5, 5)
else
  call qp_set_axis ('X', 0.0_rp, window(iw)%length*1e3, 5, 5)
endif

call qp_set_axis ('Y', real(floor(minval(y)), rp), real(ceiling(maxval(y)), rp), 20, 5)
call qp_draw_graph (x, y, 'X (in millimeters)', 'Y (in microns)', window(iw)%name)

end subroutine

end module
