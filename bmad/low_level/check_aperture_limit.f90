!+
! Subroutine check_aperture_limit (orb, ele, particle_at, param, old_orb, check_momentum)
!
! Subroutine to check if an orbit is outside an element's aperture.
! Note: A particle will also be considered to have hit an aperture
! if |p_x| or |p_y| > 1 
!
! Also see:
!   orbit_too_large
!
! Input:
!   orb            -- Coord_struct: coordinates of a particle.
!   ele            -- Ele_struct: Element holding the aperture
!   particle_at    -- Integer: first_track_edge$, second_track_edge$, surface$, in_between$
!   param          -- lat_param_struct: Lattice global parameter structure.
!   old_orb        -- Coord_struct, optional: Old coordinates at last check. 
!                       Needed if ele%aperture_at = wall_transition$. 
!                       If not present then wall transitions will be ignored.
!   check_momentum -- Logical, optional: If present and false then checking of
!                       p_x and p_y will be disabled.
!
! Output:
!   orb      -- Coord_struct: coordinates of a particle.
!     %state           -- State of the particle. Not set if particle is not past an aperture.
!   param    -- lat_param_struct: Lattice global parameter structure.
!     %unstable_factor -- Value indicates how much the parcticle is outside the aperture.
!                           This can be used in optimizations to keep particles within the aperture.
!-

recursive subroutine check_aperture_limit (orb, ele, particle_at, param, old_orb, check_momentum)

use wall3d_mod, dummy => check_aperture_limit

implicit none

type (coord_struct) :: orb
type (coord_struct), optional :: old_orb
type (coord_struct) orb2 
type (ele_struct) :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct), intent(inout) :: param

real(rp) x_width2, y_width2, x_particle, y_particle, s_here, r, rel_p, x0, y0
real(rp) x1_lim, x2_lim, y1_lim, y2_lim, dx1, dx2, dy1, dy2, d_max, r_wall
real(rp) d_radius, d_old, position(6), x_old, y_old, r_old, f, unstable_factor

integer i, particle_at, physical_end
logical do_tilt, err, no_wall, inside, old_inside
logical, optional :: check_momentum
character(*), parameter :: r_name = 'check_aperture_limit'

! Care must be taken with OpenMP since param%unstable_factor is common to all threads.

param%unstable_factor = 0
unstable_factor = 0

! Check p_x and p_y

if (logic_option(.true., check_momentum)) then
  if (orbit_too_large (orb, param)) return
endif

if (.not. bmad_com%aperture_limit_on) return

! Super_slave elements have the aperture info stored in the lord(s).

physical_end = physical_ele_end (particle_at, orb, ele%orientation)

if (ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    if (.not. lord_edge_aligned (ele, physical_end, lord) .and. lord%aperture_at /= continuous$) cycle
    call check_aperture_limit (orb, lord, particle_at, param, old_orb, check_momentum = .false.)
    if (orb%state /= alive$) return
  enddo
  return
endif

! Custom

if (ele%aperture_type == custom_aperture$) then
  call check_aperture_limit_custom (orb, ele, particle_at, param, err)
  return
endif

! Check if there is an aperture here. If not, simply return.

if (.not. at_this_ele_end (physical_end, ele%aperture_at)) return

! Check momentum limits

select case (ele%key)
case (ecollimator$)
  x_width2 = ele%value(px_aperture_width2$)
  y_width2 = ele%value(py_aperture_width2$)
  x0 = ele%value(px_aperture_center$)
  y0 = ele%value(py_aperture_center$)
  if (x_width2 == 0 .or. y_width2 == 0) then
    call check_this_limit (orb, param, orb%vec(2), x_width2, x0, unstable_factor)
    call check_this_limit (orb, param, orb%vec(4), y_width2, y0, unstable_factor)
  else
    r = ((orb%vec(2) - x0)/x_width2)**2 + ((orb%vec(4) - y0)/y_width2)**2
    if (r > 1) then
      orb%state = lost$
      unstable_factor = sqrt(r) - 1
    endif
  endif

  x_width2 = ele%value(z_aperture_width2$)
  y_width2 = ele%value(pz_aperture_width2$)
  x0 = ele%value(z_aperture_center$)
  y0 = ele%value(pz_aperture_center$)
  if (x_width2 == 0 .or. y_width2 == 0) then
    call check_this_limit (orb, param, orb%vec(5), x_width2, x0, unstable_factor)
    call check_this_limit (orb, param, orb%vec(6), y_width2, y0, unstable_factor)
  else
    r = ((orb%vec(5) - x0)/x_width2)**2 + ((orb%vec(6) - y0)/y_width2)**2
    if (r > 1) then
      orb%state = lost$
      unstable_factor = sqrt(r) - 1
    endif
  endif

case (rcollimator$)
  call check_this_limit (orb, param, orb%vec(2), ele%value(px_aperture_width2$), ele%value(px_aperture_center$), unstable_factor)
  call check_this_limit (orb, param, orb%vec(4), ele%value(py_aperture_width2$), ele%value(py_aperture_center$), unstable_factor)
  call check_this_limit (orb, param, orb%vec(5), ele%value(z_aperture_width2$), ele%value(z_aperture_center$), unstable_factor)
  call check_this_limit (orb, param, orb%vec(6), ele%value(pz_aperture_width2$), ele%value(pz_aperture_center$), unstable_factor)
end select

! A photon at the surface will have the appropriate coords already so do not need to offset.

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)

  select case (ele%aperture_at)
  case (surface$)
    if (.not. ele%offset_moves_aperture) then 
      call out_io (s_error$, r_name, 'Surface aperture must have offset_moves_aperture = True for element: ' // ele%name)
    endif
  case default
    if (ele%offset_moves_aperture) then 
      call out_io (s_error$, r_name, 'Non-Surface aperture must have offset_moves_aperture = False for element: ' // ele%name)
    endif
  end select

  x_particle = orb%vec(1)
  y_particle = orb%vec(3)

case default
  if (ele%offset_moves_aperture .and. physical_end /= surface$) then
    do_tilt = .false.
    if (ele%key == ecollimator$ .or. ele%key == rcollimator$) do_tilt = .true.
    orb2 = orb
    s_here = orb2%s - ele%s_start
    if (orb2%species == photon$) then
      call offset_photon (ele, orb2, set$)
    else
      call offset_particle (ele, set$, orb2, set_tilt = do_tilt, set_hvkicks = .false., s_pos = s_here)
    endif
    x_particle = orb2%vec(1)
    y_particle = orb2%vec(3)

  else
    x_particle = orb%vec(1)
    y_particle = orb%vec(3)
  endif
end select

! Wall transition calc. Lost when particle crosses wall from inside to out or vice versa.

x1_lim = -ele%value(x1_limit$)
x2_lim =  ele%value(x2_limit$)
y1_lim = -ele%value(y1_limit$)
y2_lim =  ele%value(y2_limit$)

if (ele%aperture_at == wall_transition$) then
  if (.not. present(old_orb)) then
    param%unstable_factor = unstable_factor
    return
  endif

  select case (ele%aperture_type)
  case (elliptical$)
    call elliptical_params_setup(ele, x1_lim, x2_lim, x_particle, x_width2, x0, err); if (err) return
    call elliptical_params_setup(ele, y1_lim, y2_lim, y_particle, y_width2, y0, err); if (err) return
    x_old = old_orb%vec(1) - x0
    y_old = old_orb%vec(3) - y0

    r_old = (x_old / x_width2)**2 + (y_old / y_width2)**2
    r = (x_particle / x_width2)**2 + (y_particle / y_width2)**2
    if (r_old > 1 .neqv. r > 1) then
      orb%state = lost$
      unstable_factor = abs(sqrt(r) - 1)
    endif

  case (rectangular$, auto_aperture$)
    x_old = old_orb%vec(1)
    y_old = old_orb%vec(3)

    inside = (x1_lim == 0 .or. x_particle > x1_lim) .and. (x2_lim == 0 .or. x_particle < x2_lim) .and. &
             (y1_lim == 0 .or. y_particle > y1_lim) .and. (y2_lim == 0 .or. y_particle < y2_lim)
    old_inside = (x1_lim == 0 .or. x_old > x1_lim) .and. (x2_lim == 0 .or. x_old < x2_lim) .and. &
                 (y1_lim == 0 .or. y_old > y1_lim) .and. (y2_lim == 0 .or. y_old < y2_lim)

    if (inside .neqv. old_inside) then
      orb%state = lost$

      if (x1_lim /= 0 .and. ((x_particle <= x1_lim .and. x_old >= x1_lim) .or. (x_particle >= x1_lim .and. x_old <= x1_lim))) then
        unstable_factor = max(unstable_factor, abs((x_particle - x1_lim) / (x2_lim - x1_lim)))
      endif

      if (x2_lim /= 0 .and. ((x_particle <= x2_lim .and. x_old >= x2_lim) .or. (x_particle >= x2_lim .and. x_old <= x2_lim))) then
        unstable_factor = max(unstable_factor, abs((x_particle - x2_lim) / (x2_lim - x1_lim)))
      endif

      if (y1_lim /= 0 .and. ((y_particle <= y1_lim .and. y_old >= y1_lim) .or. (y_particle >= y1_lim .and. y_old <= y1_lim))) then
        unstable_factor = max(unstable_factor, abs((y_particle - y1_lim) / (y2_lim - y1_lim)))
      endif

      if (y2_lim /= 0 .and. ((y_particle <= y2_lim .and. y_old >= y2_lim) .or. (y_particle >= y2_lim .and. y_old <= y2_lim))) then
        unstable_factor = max(unstable_factor, abs((y_particle - y2_lim) / (y2_lim - y1_lim)))
      endif
    endif

  case (wall3d$)
    position = wall3d_to_position (old_orb, ele)
    d_old = wall3d_d_radius (position, ele, no_wall_here = no_wall, radius_wall = r_wall)
    position = wall3d_to_position (orb, ele)
    d_radius = wall3d_d_radius (position, ele)
    if (.not. no_wall .and. (d_radius > 0 .neqv. d_old > 0)) then
      orb%state = lost$
      unstable_factor = abs(d_radius) / r_wall
    endif

  case default
    call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  end select

  param%unstable_factor = unstable_factor
  return
endif

! Non-wall transition aperture.

select case (ele%aperture_type)

case (elliptical$)
  call elliptical_params_setup(ele, x1_lim, x2_lim, x_particle, x_width2, x0, err); if (err) return
  call elliptical_params_setup(ele, y1_lim, y2_lim, y_particle, y_width2, y0, err); if (err) return
  r = (x_particle / x_width2)**2 + (y_particle / y_width2)**2
  if (r < 1) then
    param%unstable_factor = unstable_factor
    return
  endif

  if (abs(x_particle / x_width2) > abs(y_particle / y_width2)) then
    if (x_particle > 0) then; orb%state = lost_pos_x_aperture$
    else;                     orb%state = lost_neg_x_aperture$
    endif
  else
    if (y_particle > 0) then; orb%state = lost_pos_y_aperture$
    else;                     orb%state = lost_neg_y_aperture$
    endif
  endif

  unstable_factor = sqrt(r) - 1

case (rectangular$, auto_aperture$)
  if (x1_lim /= 0 .and. x_particle < x1_lim) then
    f = abs((x_particle - x1_lim) / (x2_lim - x1_lim))
    if (f > unstable_factor) then
      orb%state = lost_neg_x_aperture$
      unstable_factor = f
    endif
  endif
  
  if (x2_lim /= 0 .and. x_particle > x2_lim) then
    f = abs((x_particle - x2_lim) / (x2_lim - x1_lim))
    if (f > unstable_factor) then
      orb%state = lost_pos_x_aperture$
      unstable_factor = f
    endif
  endif

  if (y1_lim /= 0 .and. y_particle < y1_lim) then
    f = abs((y_particle - y1_lim) / (y2_lim - y1_lim))
    if (f > unstable_factor) then
      orb%state = lost_neg_y_aperture$
      unstable_factor = f
    endif
  endif
  
  if (y2_lim /= 0 .and. y_particle > y2_lim) then
    f = abs((y_particle - y2_lim) / (y2_lim - y1_lim))
    if (f > unstable_factor) then
      orb%state = lost_pos_y_aperture$
      unstable_factor = f
    endif
  endif

case (wall3d$)
  position = wall3d_to_position (orb, ele)
  d_radius = wall3d_d_radius (position, ele, no_wall_here = no_wall, radius_wall = r_wall)
  if (.not. no_wall .and. d_radius > 0) then
    orb%state = lost$
    unstable_factor = d_radius / r_wall
  endif

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
end select

param%unstable_factor = unstable_factor

!----------------------------------------------------------------------------------
contains

subroutine check_this_limit (orbit, param, coord, width2, center, unstable_factor)

type (coord_struct) orbit
type (lat_param_struct) param
real(rp) coord, width2, center, unstable_factor

if (width2 == 0) return
if (coord > center - width2 .and. coord < center + width2) return

orbit%state = lost$
unstable_factor = max((center - width2) - coord, coord - (center + width2))

end subroutine check_this_limit

!----------------------------------------------------------------------------------
! contains

subroutine elliptical_params_setup (ele, lim1, lim2, pos, width2, center, err)

type (ele_struct) ele
real(rp) lim1, lim2, pos, width2, center
logical err

!

if ((lim1 == 0 .and. lim2 /= 0) .or. (lim1 /= 0 .and. lim2 == 0)) then
  call out_io (s_error$, r_name, 'FOR AN ELLIPTICAL APERTURE ALL FOUR X1_LIMIT, X2_LIMIT, Y1_LIMIT, AND Y2_LIMIT MUST BE SET!', &
                                 '  THIS IS NOT SO FOR ELEMENT: ' // ele%name)
  err = .true.
  return
endif

width2 = (lim2 - lim1) / 2
center = (lim2 + lim1) / 2
pos = pos - center

if (width2 == 0) width2 = bmad_com%max_aperture_limit
width2 = max(width2, 0.0_rp) ! Make sure not negative

err = .false.

end subroutine elliptical_params_setup

end subroutine check_aperture_limit
