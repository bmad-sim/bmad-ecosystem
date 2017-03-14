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
! Modules needed:
!   use track1_mod
!
! Input:
!   orb            -- Coord_struct: coordinates of a particle.
!   ele            -- Ele_struct: Element holding the aperture
!   particle_at    -- Integer: first_track_edge$, second_track_edge$, surface$, in_between$
!   param          -- lat_param_struct: Parameter structure
!   old_orb        -- Coord_struct, optional: Old coordinates at last check. 
!                       Needed if ele%aperture_at = wall_transition$. 
!                       If not present then wall transitions will be ignored.
!   check_momentum -- Logical, optional: If present and false then checking of
!                       p_x and p_y will be disabled.
!
! Output:
!   orb   -- Coord_struct: coordinates of a particle.
!     %state -- State of the particle
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

real(rp) x_lim, y_lim, x_particle, y_particle, s_here, r, rel_p, x0, y0
real(rp) x1_lim, x2_lim, y1_lim, y2_lim, dx1, dx2, dy1, dy2, d_max
real(rp) d_radius, d_old, position(6), x_old, y_old, r_old

integer i, particle_at, physical_end
logical do_tilt, err
logical, optional :: check_momentum
character(*), parameter :: r_name = 'check_aperture_limit'

! Super_slave elements have the aperture info stored in the lord

param%unstable_factor = -1

physical_end = physical_ele_end (particle_at, orb%direction, ele%orientation)

if (ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    if (.not. lord_edge_aligned (ele, physical_end, lord) .and. lord%aperture_at /= continuous$) cycle
    call check_aperture_limit (orb, lord, particle_at, param, old_orb, check_momentum)
    if (orb%state /= alive$) return
  enddo
  return
endif

! Custom

if (ele%aperture_type == custom_aperture$) then
  call check_aperture_limit_custom (orb, ele, particle_at, param, err)
  return
endif

! Check p_x and p_y

if (logic_option(.true., check_momentum)) then
  if (orbit_too_large (orb, param)) return
endif

! Check if there is an aperture here. If not, simply return.

if (.not. at_this_ele_end (physical_end, ele%aperture_at)) return

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
      call offset_particle (ele, param, set$, orb2, set_tilt = do_tilt, &
                                 set_multipoles = .false., set_hvkicks = .false., ds_pos = s_here)
    endif
    x_particle = orb2%vec(1)
    y_particle = orb2%vec(3)

  else
    x_particle = orb%vec(1)
    y_particle = orb%vec(3)
  endif
end select

!

x_lim = (ele%value(x1_limit$) + ele%value(x2_limit$)) / 2
y_lim = (ele%value(y1_limit$) + ele%value(y2_limit$)) / 2
x0 = (ele%value(x2_limit$) - ele%value(x1_limit$)) / 2
y0 = (ele%value(y2_limit$) - ele%value(y1_limit$)) / 2
x_particle = x_particle - x0
y_particle = y_particle - y0

if (.not. bmad_com%aperture_limit_on .or. x_lim == 0) x_lim = bmad_com%max_aperture_limit
if (.not. bmad_com%aperture_limit_on .or. y_lim == 0) y_lim = bmad_com%max_aperture_limit

! 

if (ele%aperture_at == wall_transition$) then
  if (.not. present(old_orb)) return

  x_old = old_orb%vec(1) - x0
  y_old = old_orb%vec(3) - y0

  select case (ele%aperture_type)
  case (elliptical$)
    r_old = (x_old / x_lim)**2 + (y_old / y_lim)**2
    r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
    param%unstable_factor = r - 1
    if (r_old > 1 .xor. r > 1) orb%state = lost$

  case (rectangular$, auto_aperture$)
    if (abs(x_particle) < x_lim .and. abs(y_particle) < y_lim) then
      param%unstable_factor = 1 - max(abs(x_particle)/x_lim, abs(y_particle)/y_lim)
    else
      param%unstable_factor = min(abs(x_particle)/x_lim, abs(y_particle)/y_lim) - 1
    endif

    if ((abs(x_particle) < x_lim .and. abs(y_particle) < y_lim) .xor. &
        (abs(x_old) < x_lim .and. abs(y_old) < y_lim)) then
      orb%state = lost$
    endif

  case (wall3d$)
    position = [old_orb%vec(1:4), old_orb%s-ele%s_start, 1.0_rp]
    d_old = wall3d_d_radius (position, ele)
    position = [orb%vec(1:4), orb%s-ele%s_start, 1.0_rp]
    d_radius = wall3d_d_radius (position, ele)
    param%unstable_factor = d_radius - 1
    if (d_radius > 1 .xor. d_old > 1) orb%state = lost$

  case default
    call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  end select

  return
endif

!

select case (ele%aperture_type)

case (elliptical$)
  r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
  param%unstable_factor = sqrt(r) - 1
  if (r < 1) return

  if (abs(x_particle / x_lim) > abs(y_particle / y_lim)) then
    if (x_particle > 0) then; orb%state = lost_pos_x_aperture$
    else;                     orb%state = lost_neg_x_aperture$
    endif
  else
    if (y_particle > 0) then; orb%state = lost_pos_y_aperture$
    else;                     orb%state = lost_neg_y_aperture$
    endif
  endif

case (rectangular$, auto_aperture$)
  if (abs(x_particle/x_lim) > abs(y_particle/y_lim)) then
    param%unstable_factor = abs(x_particle/x_lim) - 1
  else
    param%unstable_factor = abs(y_particle/y_lim) - 1
  endif

  if (abs(x_particle) < x_lim .and. abs(y_particle) < y_lim) return

  if (abs(x_particle/x_lim) > abs(y_particle/y_lim)) then
    if (x_particle > 0) then
      orb%state = lost_pos_x_aperture$
    else
      orb%state = lost_neg_x_aperture$
    endif

  else
    if (y_particle > 0) then
      orb%state = lost_pos_y_aperture$
    else
      orb%state = lost_neg_y_aperture$
    endif
  endif

case (wall3d$)
  position = [orb%vec(1:4), orb%s-ele%s_start, 1.0_rp]
  d_radius = wall3d_d_radius (position, ele)
  param%unstable_factor = d_radius - 1
  if (d_radius > 1) then
    orb%state = lost$
  endif

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN APERTURE_TYPE FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
end select

end subroutine check_aperture_limit
