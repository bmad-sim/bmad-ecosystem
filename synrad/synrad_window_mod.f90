module synrad_window_mod

use sr_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine find_windows (outside, window)

  implicit none

  type (wall_struct) outside
  type (crotch_window_struct) window(:)

  integer i, ix
  real(rp) ds_window, dx_window

  !

  i = 0

  do ix=1,outside%n_pt_tot          ! search through outside wall points

    if (outside%pt(ix)%phantom) then   ! checking for phantoms,
      i = i + 1
      window(i)%ix_pt = ix            ! adding them as windows
      window(i)%n_ray_hit = 0         ! initial ray hit counter

! This next part should be automated at some point...

      if ((outside%pt(ix)%s >18) .and. (outside%pt(ix)%s <22)) then
        window(i)%name = "CROTCH 4W"    ! set window name
        window(i)%layout = reverse$
        window(i)%side = w_west$
      elseif ((outside%pt(ix)%s >23) .and. (outside%pt(ix)%s <27)) then
        window(i)%name = "CROTCH 5W"    ! set window name
        window(i)%layout = reverse$
        window(i)%side = w_west$
      elseif ((outside%pt(ix)%s >28) .and. (outside%pt(ix)%s <32)) then
        window(i)%name = "CROTCH 6W"    ! set window name
        window(i)%layout = reverse$
        window(i)%side = w_west$
      elseif ((outside%pt(ix)%s >60) .and. (outside%pt(ix)%s <70)) then
        window(i)%name = "CROTCH 11W"    ! set window name
        window(i)%layout = forward$
        window(i)%side = w_west$
      elseif ((outside%pt(ix)%s >736) .and. (outside%pt(ix)%s <740)) then
        window(i)%name = "CROTCH 6E"    ! set window name
        window(i)%layout = forward$
        window(i)%side = w_east$
      elseif ((outside%pt(ix)%s >741) .and. (outside%pt(ix)%s <745)) then
        window(i)%name = "CROTCH 5E"    ! set window name
        window(i)%layout = forward$
        window(i)%side = w_east$
      elseif ((outside%pt(ix)%s >746) .and. (outside%pt(ix)%s <750)) then
        window(i)%name = "CROTCH 4E"    ! set window name        
        window(i)%layout = forward$
        window(i)%side = w_east$
      endif
      if (window(i)%layout == reverse$) then
        ds_window = outside%pt(ix)%s - outside%pt(ix-1)%s
        dx_window = outside%pt(ix)%x - outside%pt(ix-1)%x
        window(i)%angle = atan2 (ds_window, dx_window)
      else
        ds_window = outside%pt(ix-1)%s - outside%pt(ix)%s
        dx_window = outside%pt(ix-1)%x - outside%pt(ix)%x
        window(i)%angle = atan2 (ds_window, dx_window)
      endif
      window(i)%length = sqrt(dx_window**2 + ds_window**2)
      outside%pt(ix)%phantom = .false.
    endif

  enddo

  print *, 'Found ',i,' windows.'

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Finds vertical sigma effective at the window
! 
! input:
!           window   :  array of windows with window ray hit info
!           gen       :  contains vertical emittance of beam
! output:
!           window   :  array of windows with sigma information


subroutine sigma_at_windows (window, gen)

implicit none

type (crotch_window_struct),target :: window(:)
type (synrad_param_struct) gen

integer i, j, dist, count, ix
type (coord_struct) wind,targ
type (ray_hit_struct), pointer :: ray_hit

!

do i=1,n_windows$

  do j=1,window(i)%n_ray_hit

    ray_hit => window(i)%ray_hit_(j)
    ray_hit%sig_y = sqrt(gen%epsilon_y * ray_hit%ray%y_twiss%beta)
    ray_hit%sig_yp = sqrt(gen%epsilon_y / ray_hit%ray%y_twiss%beta)
    ray_hit%window_sig_y = sqrt(ray_hit%sig_y**2 + &
              (ray_hit%ray%track_len * ray_hit%sig_yp)**2)
      
  enddo
enddo

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Projects rays from windows to a target distance.
! Gives target coordinates and sigmas
! input:
!           window   :  array of windows with window hit info
! output:
!           window   :  array of windows with target information

subroutine project_from_windows (window)

  implicit none

  type (crotch_window_struct),target :: window(n_windows$)

  integer iw_(n_windows$)
  integer i, j, dist
  type (coord_struct) wind,targ
  type (ray_hit_struct), pointer :: ray_hit

  !

  call get_window_numbers(window, iw_)
  print *, 'Enter distance from window to project rays:'
  print '(a,$)', ' (in meters) > '
  read *, dist

  if ((dist > 0) .and. (dist < 1000)) then
    do i=1,n_windows$

      if (iw_(i) == 0) exit
      do j=1,window(iw_(i))%n_ray_hit

        ray_hit => window(iw_(i))%ray_hit_(j)
        wind = ray_hit%hit_coord
        targ = wind
        targ%vec(1) = wind%vec(1) + (wind%vec(2) * dist)
        targ%vec(3) = wind%vec(3) + (wind%vec(4) * dist)
        ray_hit%target_coord = targ
        ray_hit%dist = dist
        ray_hit%sig_y_eff = sqrt(ray_hit%sig_y**2 + &
                  ((dist + ray_hit%ray%track_len) * ray_hit%sig_yp)**2)
        
      enddo
    enddo
  else
    print *, 'Distance must be between 0 and 1000 meters.'
  endif

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Asks for which crotches and
! gives corresponding window index numbers
! input:
!           window   :  array of windows
! output:
!           iw       :  array of integers for requested windows

subroutine get_window_numbers (window, iw)

  implicit none

  type (crotch_window_struct) window(:)
  integer iw(n_windows$)

  integer i, j, length, count
  character*80 line

  iw(1:n_windows$) = 0
  print *, 'Enter crotch window(s) of interest, separated by a space,'
  print '(a,$)', ' (4E 5E 6E 4W 5W 6W 11W WEST EAST or ALL) > '
  accept '(a)', line
  call str_upcase(line, line)
  count = 0
  do

    line = adjustl(line)
    line = trim(line)
    length = len(line)  
    j = index(line, ' ')
    if ((length == 0) .or. (j == 1)) exit
    if (j == 0) j = length

    if (line(1:1) == 'A') then            ! for command 'ALL' select all
      do i=1,n_windows$

        count = count + 1
        iw(count) = i

      enddo

    elseif (line(1:1) == 'E') then        ! select all EAST windows
      do i=1,n_windows$

        if (window(i)%side == w_east$) then
          count = count + 1
          iw(count) = i
        endif

      enddo
    elseif (line(1:1) == 'W') then        ! select all WEST windows
      do i=1,n_windows$

        if (window(i)%side == w_west$) then
          count = count + 1
          iw(count) = i
        endif

      enddo
    else                                  ! select individual windows
      do i=1,n_windows$

        if (index(window(i)%name, line(1:j-1)) /= 0) then
          count = count + 1
          iw(count) = i
        endif

      enddo
    endif

    if (j == length) exit
    line = line(j+1:length)
    if (count == n_windows$) exit
  enddo

  if (count == 0) print *,'No valid windows were found.'
end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine window_ray_reset (window)

  implicit none

  type (crotch_window_struct) window(:)

  integer i

  !

  do i=1,n_windows$
      window(i)%n_ray_hit = 0         ! initial ray hit counter
  enddo

  type *, 'Window ray counter reset'

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine calculate_wind_power (ring, orb, direction, power, walls, gen, window)

  implicit none

  type (walls_struct), target :: walls
  type (lat_struct) ring
  type (coord_struct) orb(:)
  type (synrad_param_struct) gen
  type (ele_power_struct) power(:)
  type (crotch_window_struct) window(:)

  integer direction, ie

! init

  power%at_wall = 0
  power%radiated = 0

! loop over all elements

  do ie = 1, ring%n_ele_track
    call ele_wind_power (ring, ie, orb, direction, power, walls, gen, window)
  enddo

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine ele_wind_power (ring, ie, orb, direction, power, walls, gen, window)

  implicit none

  type (walls_struct), target :: walls
  type (lat_struct), target :: ring
  type (coord_struct) orb(:)
  type (ele_struct), pointer :: ele
  type (ray_struct) rays(3000), ray, ray_temp
  type (synrad_param_struct) gen
  type (ele_power_struct) power(:)
  type (crotch_window_struct) window(:)

  integer direction, ie, i_ray, n_slice, ns

  real(rp) l_off, del_l, l0, l1, l_try

! power calculation is only for bends, quads and wigglers.

  ele => ring%ele(ie)

  if (ele%key /= sbend$ .and. ele%key /= quadrupole$ .and. &
                     ele%key /= wiggler$ .and. ele%key /= sol_quad$) return
  if (ele%s > 190.0 .and. ele%s < 570) return   ! Ignore north half of ring

! partition the quad/bend/wiggler into sections and track the synch
! radiation comming from the ends of the sections to see where they hit
! the wall

! If different rays hit both inside and outside then we need to find the ray
! that just misses the inside wall and break up the fan accordingly into
! two pieces

!!      if (abs(ele%value(tilt$)) > 20*twopi/360) cycle   ! ignore skew quads

  i_ray = 0
  n_slice = 50

  if (ele%key == wiggler$) then
    if (ele%value(n_pole$) == 0) then
      type *, 'WARNING IN ELE_WIND_POWER: "N_POLE" FOR WIGGLER = 0.'
      type *, '      CALCULATED RADIATION FROM ', trim(ele%name),' WILL BE 0!'
      return
    endif
    n_slice = 50 * ele%value(n_pole$)
  endif

  del_l = ele%value(l$) / n_slice

  do ns = 0, n_slice

    l_off = ns * del_l
    i_ray = i_ray + 1
    call init_ray (rays(i_ray), ring, ie, l_off, orb, direction)
    call track_ray_to_window (rays(i_ray), ring, walls, window)

    if (i_ray > 1 .and. rays(i_ray)%wall%side /= rays(i_ray-1)%wall%side) then
      ray_temp = rays(i_ray)
      rays(i_ray) = rays(i_ray-1)
      l0 = l_off - del_l
      l1 = l_off
      do
        l_try = (l0 + l1) / 2
        call init_ray (ray, ring, ie, l_try, orb, direction)
        call track_ray_to_window (ray, ring, walls, window)
        if (ray%wall%side == rays(i_ray)%wall%side) then
          rays(i_ray) = ray
          l0 = l_try
        else
          ray_temp = ray
          l1 = l_try
        endif
        if ((l1 - l0) .le. 5e-4) then
          if (abs(rays(i_ray)%start%vec(5) - &
                          rays(i_ray-1)%start%vec(5)) < 5e-4) i_ray = i_ray - 1
          call seg_power_calc (rays, i_ray, walls, ring, gen, power(ie))
          rays(1) = ray_temp
          i_ray = 1
          exit
        endif
      enddo
    endif

  enddo

! now compute the sr power hitting the wall using the tracking results
! and interpolating inbetween.

  call seg_power_calc (rays, i_ray, walls, ring, gen, power(ie))

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine track_ray_to_window (ray, ring, walls, window, hit_flag, track_max)

  implicit none

  type (lat_struct), target :: ring
  type (ray_struct), target :: ray
  type (walls_struct), target :: walls
  type (wall_struct), pointer :: inside, outside
  type (crotch_window_struct) window(:)

  logical, optional :: hit_flag
  real(rp), optional :: track_max

  integer ix_in, ix_out

  real(rp) s_next

  logical is_hit, passed_end

! init

inside => walls%negative_x_wall
outside => walls%positive_x_wall
passed_end = .false.

!  type *, hit_flag
!  if (present(hit_flag)) hit_flag = .true.  ! assume that we will hit
  if (present(hit_flag)) type *,"present"    ! assume that we will hit
 
! ix_in and ix_out are the next inside and outside wall points that
! are at or just "downstream" of the ray.

  call get_initial_pt (ray, inside, ix_in, ring)
  call get_initial_pt (ray, outside, ix_out, ring)

! propagation loop:
! Propagate the ray. Figure out how far to advance in s.
! Do not advance past the next wall point (either inside or outside).

  do

    if (ray%direction == 1) then
      s_next = min(inside%pt(ix_in)%s, outside%pt(ix_out)%s, ray%now%vec(5) + 1.0)
    else
      s_next = max(inside%pt(ix_in)%s, outside%pt(ix_out)%s, ray%now%vec(5) - 1.0)
    endif

    call propagate_ray (ray, s_next, ring, .true., walls%circular)

! See if we are outside the beam pipe.
! If so we calculate the exact hit spot where the ray crossed the
! wall boundry and return

    call hit_spot_calc_wind (ray, inside, ix_in, is_hit, ring, window, walls%circular)
    if (is_hit) return

    call hit_spot_calc_wind (ray, outside, ix_out, is_hit, ring, window, walls%circular)
    if (is_hit) return

    if (present(track_max)) then
      if (ray%track_len .ge. track_max) then
        hit_flag = .false.
        return
      endif
    endif

    if (ray%now%vec(5) == inside%pt(ix_in)%s) &
              call next_pt (ray, inside, ix_in, walls%circular, passed_end)
    if (ray%now%vec(5) == outside%pt(ix_out)%s) &
              call next_pt (ray, outside, ix_out, walls%circular, passed_end)

  enddo

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine hit_spot_calc_wind (ray, wall, ix_wall, is_hit, ring, window, circular)

  implicit none

  type (lat_struct) ring
  type (ray_struct) :: ray, ray0, ray1, ray2
  type (wall_struct), target :: wall
  type (wall_pt_struct), pointer :: pt(:)
  type (wall_pt_struct), pointer :: pt0
  type (crotch_window_struct) window(:)

  integer ix_wall, ix0, ix1, ix2, i

  real(rp) dx_wall, ds_wall, del_s, s1, denom
  real(rp) del0, del1, del2, x, r_wall

  type (coord_struct) temp

  logical is_hit, circular

!  type *, "hit_spot_calc_wind starting"

! init

  pt => wall%pt
  pt0 => pt(ix_wall)

  is_hit = .false.  ! assume no hit

! if the ray is on the opposite side of the beam centerline from the wall there
! is no hit

  if (sign(1.0, ray%now%vec(1)) /= sign(1.0, pt0%x)) return

! figure out if there is a hit.
! if there is no alley here then we have not hit if x_pos < x_wall

  if (pt0%type == no_alley$) then
    ray%alley_status = no_local_alley$
    ix1 = ix_wall - ray%direction
    if (ray%now%vec(5) == pt0%s) then
      x = pt0%x
    else
      x = pt0%x + (pt(ix1)%x - pt0%x) * &             ! interpolate to get x
                    (ray%now%vec(5) - pt0%s) / (pt(ix1)%s - pt0%s)
    endif
    if (abs(ray%now%vec(1)) < abs(x)) return

! Here for an alley...
! Note. To simplify the program we do not allow a calculation of whether
! we have hit if we are in an ally but between wall points (that is, in
! an alley the wall points cannot be further apart than 1 meter).

  elseif (ray%now%vec(5) == pt0%s) then   ! in alley at a wall pt

    if (pt0%type == open_end$) then
      if (ray%direction == pt0%closed_end_direct) then
        if (abs(ray%now%vec(1)) .gt. abs(pt0%x)) then
          ray%alley_status = in_alley$
        else
          ray%alley_status = out_of_alley$
        endif
        ix1 = ix_wall - 1
        if (pt(ix1)%type == open_end$) then
          if (abs(ray%now%vec(1)) > max(abs(pt(ix1)%x), abs(pt0%x))) return
          if (abs(ray%now%vec(1)) < min(abs(pt(ix1)%x), abs(pt0%x))) return
        else
          return
        endif
      else
        if (ray%alley_status == in_alley$) then
          ix1 = ix_wall + ray%direction
          if (pt(ix1)%type == open_end$) return
          if (abs(ray%now%vec(1)) > abs(pt0%x)) return
        else
          ix1 = ix_wall - ray%direction
          if (pt(ix1)%type == open_end$) return
          if (abs(ray%now%vec(1)) < abs(pt0%x)) return
        endif
      endif

    elseif (pt0%type == closed_end$) then
      if (ray%direction /= pt0%closed_end_direct) return
      if (ray%alley_status /= in_alley$) return
      if (abs(ray%now%vec(1)) < abs(pt0%x)) then
        ix1 = ix_wall + ray%direction
        if (pt(ix1)%type == closed_end$ .and. &
                                     abs(ray%now%vec(1)) < abs(pt(ix1)%x)) return
      else
        ix1 = ix_wall - ray%direction
        if (pt(ix1)%type == closed_end$ .and. &
                                    abs(ray%now%vec(1)) > abs(pt(ix1)%x)) return
      endif

    elseif (pt0%type == inner_wall$) then
      if (ray%alley_status == in_alley$) return
      ix1 = ix_wall - ray%direction
      if (abs(ray%now%vec(1)) < abs(pt0%x)) return

    elseif (pt0%type == middle_wall$) then
      if (ray%alley_status /= in_alley$) return
      ix1 = ix_wall + ray%direction
      if (abs(ray%now%vec(1)) > abs(pt0%x)) return

    elseif (pt0%type == outer_wall$) then
      if (ray%alley_status /= in_alley$) return
      ix1 = ix_wall - ray%direction
      if (abs(ray%now%vec(1)) < abs(pt0%x)) return

    else
      type *, 'ERROR IN HIT_SPOT_CALC: ALLEY ERROR'
      call err_exit
    endif

  else
    if (pt0%type /= no_alley$) then
      type *, 'ERROR IN HIT_SPOT_CALC: CALCULATION IN ALLEY BUT BETWEEN POINTS.'
      call err_exit
    endif
  endif

  ix0 = min(ix_wall, ix1)
  ix2 = max(ix_wall, ix1)
  is_hit = .true.
  ray%ix_wall_pt = ix2

! Here if we have a hit.
! We need to find where exactly the ray hit the wall.
! ray%r_wall is the percentage distance along the wall piece from pt(ix0)
! where the ray hits.
! ray%r_wall = 0.0 => the hit is at pt(ix0)
! ray%r_wall = 1.0 => the hit is at pt(ix2)
! we need to iterate in a bend since the wall is actually curved.


  ray0 = ray
  ray1 = ray
  ray2 = ray

  if (ray%now%vec(5) < ray%old%vec(5)) then
    ray2%now = ray%old
  elseif (ray%now%vec(5) > ray%old%vec(5)) then
    ray0%now = ray%old
  endif

  dx_wall = pt(ix2)%x - pt(ix0)%x
  ds_wall = pt(ix2)%s - pt(ix0)%s
  denom = sqrt (dx_wall**2 + ds_wall**2)

  del0 = (dx_wall*(ray0%now%vec(5) - pt(ix0)%s) - &
                             ds_wall*(ray0%now%vec(1) - pt(ix0)%x)) / denom
  del2 = (dx_wall*(ray2%now%vec(5) - pt(ix0)%s) - &
                             ds_wall*(ray2%now%vec(1) - pt(ix0)%x)) / denom
                                    
  do i = 1, 20

    if (abs(del0) < 1.0e-4) then
      ray1 = ray0
      exit
    elseif (abs(del2) < 1.0e-4) then
      ray1 = ray2
      exit
    endif

    if (i == 20) then
      type *, 'ERROR IN HIT_SPOT_CALC: CALCULATION IS NOT CONVERGING'
      call err_exit
    endif  

    s1 = (del2 * ray0%now%vec(5) - del0 * ray2%now%vec(5)) / (del2 - del0)


    if (s1 < min(pt(ix0)%s, pt(ix2)%s) .or. s1 > max(pt(ix0)%s, pt(ix2)%s)) then
      type *, 'ERROR IN HIT_SPOT_CALC: INTERPOLATION ERROR'
      call err_exit
    endif

    if (s1 < ray1%now%vec(5)) then
      ray1%direction = -1
    else
      ray1%direction = +1
    endif
    call propagate_ray (ray1, s1, ring, .false., circular)

    del1 = (dx_wall*(ray1%now%vec(5) - pt(ix0)%s) - &
                        ds_wall*(ray1%now%vec(1) - pt(ix0)%x)) / denom
                   
    if (s1 < ray0%now%vec(5)) then
      ray2 = ray0; del2 = del0
      ray0 = ray1; del0 = del1
    elseif (s1 > ray2%now%vec(5)) then
      ray0 = ray2; del0 = del2
      ray2 = ray1; del2 = del1
    elseif (sign(1.0, del0) == sign(1.0, del1)) then
      ray0 = ray1; del0 = del1
    else
      ray2 = ray1; del2 = del1
    endif

  enddo

! cleanup

  r_wall = (dx_wall*(ray1%now%vec(1) - pt(ix0)%x) + &
                ds_wall*(ray1%now%vec(5) - pt(ix0)%s)) / denom**2
  if (r_wall > 1 .and. r_wall <  1.0001) r_wall = 1
  if (r_wall < 0 .and. r_wall > -0.0001) r_wall = 0

  if (r_wall > 1 .or. r_wall < 0) then
    type *, 'ERROR IN HIT_SPOT_CALC: R_WALL OUT OF BOUNDS.', r_wall
    call err_exit
  endif

  ray%now = ray1%now
  ray%r_wall = r_wall

  del_s = ray%now%vec(5) - ray%start%vec(5)
  if (ray%crossed_end) then
    ray%track_len = ring%param%total_length - abs(del_s)
  else      
    ray%track_len = abs(del_s)
  endif

  ray%wall => wall

! check for hit crotch window

  do i = 1, n_windows$

    if ( window(i)%ix_pt == ray%ix_wall_pt ) then
      window(i)%n_ray_hit = window(i)%n_ray_hit + 1
      window(i)%ray_hit_(window(i)%n_ray_hit)%ray = ray
      temp = ray%now
      if (window(i)%side == forward$) then
        temp%vec(1) = (1 - ray%r_wall) * window(i)%length
      else
        temp%vec(1) = ray%r_wall * window(i)%length
      endif
      temp%vec(2) = tan(atan(ray%now%vec(2)) + window(i)%angle)
      window(i)%ray_hit_(window(i)%n_ray_hit)%hit_coord = temp
    endif

  enddo
end subroutine

end module
