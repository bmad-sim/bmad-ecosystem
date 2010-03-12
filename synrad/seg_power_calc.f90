!+
! subroutine seg_power_calc (fan, i_ray, walls, wall_side, lat, gen, power)
!
! Subroutine to calculate the synch radiation power hitting a wall.
! It is assumed that the radiation is hitting only one wall.
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   fan(:) -- ray_struct: array of rays from one lat element
!   i_ray  -- integer: index of highest ray used
!   lat    -- lat_struct: with twiss propagated and mat6s made
!   walls  -- walls_struct: Both walls and ends
!   wall_side -- Integer: +x or -x side where rays are hitting.
!   gen    -- synrad_param_struct: Contains lat name, vert emittance, and beam current.
!
! Output:
!   power(*)  -- ele_power_struct: power radiated from a lat ele
!   walls   -- wall_struct: both walls and ends with power information
!-

subroutine seg_power_calc (fan, i_ray, walls, wall_side, lat, gen, power)

use synrad_struct
use synrad_interface, except => seg_power_calc

implicit none

type (ray_struct), target :: fan(:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: wall
type (wall_seg_struct), pointer :: seg
type (ray_struct) :: ray
type (ray_struct), pointer :: ray1, ray2
type (synrad_param_struct) gen
type (lat_struct) lat
type (ele_power_struct) power

real(rp) energy, dx_seg, x1, x2
real(rp) sig_y, sig_yp, sig_y_eff, dx, ds, theta_rel
real(rp) rr, dr2, dr1, theta, power_factor
real(rp) dx_wall, ds_wall, r1, r2, track_len
real(rp) r_i1, r_i2, frac_illum, gamma
real(rp) theta_base, theta_floor1, theta_floor2, theta_floor_seg

integer i_ray, iw, i, ip, is, iss, wall_side
integer type0, type1, i_pt, n_pos, n_neg

logical hit_flag

! get energy from source element

energy = lat%ele(fan(1)%ix_source)%value(e_tot$)

! 14.1e3 [m (eV)^-3] (used below) is  Cgamma * 1e9 / (2 pi) 
!    Cgamma is from Sands (p98) which is 8.85e-5 m (GeV)^-3 
!    for electrons.  Not valid for protons, etc!!!

power_factor = 14.1e3 * gen%i_beam * (energy/1e9)**4
do i = 2, i_ray

  ray1 => fan(i-1)
  ray2 => fan(i)

  ! Sum up the power radiated by the element:
  ! above power_factor * delta S * avg(gbend^2)
  ! which is an integrated segment of Sands eq. 4.5
  ! * the current to get the power

  power%radiated = power%radiated + power_factor * &
         (ray2%start%vec(5) - ray1%start%vec(5)) * (ray2%g_bend**2 + ray1%g_bend**2) / 2
enddo
 
! Ignore element if the radiated power is 0
if (power%radiated == 0) return


! wall is the side that is being hit by the radiation.
! If wall_side = 0 (only hitting ends) then does not matter what wall is set to.

if (wall_side == positive_x$) then
  wall => walls%positive_x_wall
else
  wall => walls%negative_x_wall
endif

! calculate the power factors for each ray 
do i = 1, i_ray

  ray2 => fan(i)

  iw = ray2%ix_wall_pt
  

  ! should update to include energy and dispersion component!!!
  ! sig_y is the sqrt of the vert emitt * vert beta
  sig_y = sqrt(gen%epsilon_y * ray2%y_twiss%beta)

  ! sig_yprime is the vertical opening angle

  call convert_total_energy_to (energy, lat%param%particle, gamma)
  sig_yp = sqrt(gen%epsilon_y * ray2%y_twiss%gamma + 1 / gamma**2)

  ! the effective sig_y includes the increase from sig_yp
  sig_y_eff = sqrt(sig_y**2 + (ray2%track_len * sig_yp)**2)

  ! the p1 factor is used in power/length
  ray2%p1_factor = power_factor * ray2%g_bend

  ! the p2 factor is used in power/area
  ray2%p2_factor = ray2%p1_factor / (sqrt(twopi) * sig_y_eff)

enddo

! fill the segment array by interpolating between the fan rays.
! theta_rel is the relative angle between the ray and the wall
! Note: the power at a wall segment could be due to more than one source.
! The code tries to set %ix_ele_source and %s_source to reflect the
! largest source

do i = 2, i_ray

  ray1 => fan(i-1)
  ray2 => fan(i)

  ! loop over all wall segments whose mid-point falls between the hit points
  ! of rays ray1 and ray2. We must take into account rays that straddle
  ! the lattice ends for circular lattices.

  if (ray1%wall_side == start_side$) then
    r1 = 0
  elseif (ray1%wall_side == exit_side$) then
    r1 = wall%n_seg_tot
  else
    i_pt = ray1%ix_wall_pt
    r1 = wall%pt(i_pt)%ix_seg + ray1%r_wall * wall%pt(i_pt)%n_seg + 1
  endif

  if (ray2%wall_side == start_side$) then
    r2 = 0
  elseif (ray2%wall_side == exit_side$) then
    r2 = wall%n_seg_tot
  else
    i_pt = ray2%ix_wall_pt
    r2 = wall%pt(i_pt)%ix_seg + ray2%r_wall * wall%pt(i_pt)%n_seg + 1
  endif

  ! There is a complication if, for a circular lattice, the fan is 
  ! illuminating the lattice ends.

  if (ray1%crossed_end .neqv. ray2%crossed_end) then
    r_i1 = max(r1, r2)
    r_i2 = min(r1, r2) + wall%n_seg_tot   ! This addition will be "mod"ed out below
  else
    r_i1 = min(r1, r2)
    r_i2 = max(r1, r2)
  endif

  ! theta_floor1 and theta_floor2 are the ray angles in the floor coordinate system.

  theta_base = theta_floor(ray1%now%vec(5), lat) 
  theta_floor1 = ray1%now%vec(2) + theta_base
  theta_floor2 = ray2%now%vec(2) + theta_floor(ray2%now%vec(5), lat, theta_base)

  ! Loop over all segments that have light hitting it.

  do iss = int(r_i1), int(r_i2)

    if (lat%param%lattice_type == linear_lattice$ .and. iss > wall%n_seg_tot) exit

    is = iss
    if (is > wall%n_seg_tot) is = is - wall%n_seg_tot
    seg => wall%seg(is)

    ip = seg%ix_pt              ! wall point associated with this segment

    ! dr1 is the perpendicular distance between the ray1 ray and the segment mid-point.
    ! dr2 is the perpendicular distance between the ray2 ray and the segment mid-point.
    ! Notice that typically dr1 and dr2 have opposite signs.

    dr1 = ray_to_seg_distance (ray1, seg)
    dr2 = ray_to_seg_distance (ray2, seg)

    if ((dr1 - dr2) == 0) then ! If both rays hit the same spot...
      rr = 1  ! then just use theta_floor2
    else
      rr = dr1 / (dr1 - dr2)
    endif

    ! If both rays are on the same side of the segment midpoint then just 
    !   use the closest ray.

    if (rr > 1) rr = 1
    if (rr < 0) rr = 0

    ! In a bend, the angle w.r.t. the center line can change rapidly with s.
    ! For the interpolation it is therefore better to use the floor angle.
    ! Note: Since increasing floor angle is opposite the bend angle, there is 
    ! a plus sign.

    ! theta_floor_seg is the ray angle averaged over theta_floor1 and theta_floor2.
    ! The averaging factor is rr.

    theta_floor_seg = (1 - rr) * theta_floor1 + rr * theta_floor2
    theta = theta_floor_seg - theta_floor(seg%s_mid, lat, theta_base)

    ! theta_rel gives the graze angle of the ray with respect to the wall surface.
    ! theta_rel = 0 means the ray and wall surface are parallel.

    dx_wall = wall%pt(ip)%x - wall%pt(ip-1)%x
    ds_wall = wall%pt(ip)%s - wall%pt(ip-1)%s
    theta_rel = modulo2(theta - atan2(dx_wall, ds_wall), pi) ! -pi <= theta_rel < pi

    ! A given wall segment may be shadowed by a part of the wall upstream so we need
    !   to check if this is the case. 
    ! First, if the ray is stricking the wall from the back side then this segment 
    !   is shadowed.

    if (wall%side == positive_x$ .and. theta_rel*ray2%direction < 0) cycle
    if (wall%side == negative_x$ .and. theta_rel*ray2%direction > 0) cycle

    ! Normally wall%pt(i)%s is an increasing function of i. 
    ! When there are x-ray lines with crothes things are more complicated.
    ! An x-ray line branching off from the beam pipe is called an "alley". 
    ! Where there is an alley, at a given s position, there can be three wall pieces:
    !  The outer and middle wall pieces define the x-ray line.
    !  The inner wall piece defines the beam pipe.

    type0 = wall%pt(ip-1)%type
    type1 = wall%pt(ip)%type
    if (type0 == outer_wall$ .or. type1 == outer_wall$) then
      ray%alley_status = in_alley$
    elseif (type0 == middle_wall$ .or. type1 == middle_wall$) then
      ray%alley_status = in_alley$
    elseif (type0 == inner_wall$ .or. type1 == inner_wall$) then
      ray%alley_status = out_of_alley$
    else
      ray%alley_status = no_local_alley$
    endif

    ! Track backwards from the segment midpoint to see if this segment is shadowed.
    ! If a ray can travel backwards a distance equal to the distance the rays 
    ! traveled forward without hitting a wall then there is no shadow.

    ray%start%vec(1) = seg%x_mid
    ray%start%vec(5) = seg%s_mid
    ray%start%vec(2) = theta
    ray%ix_ele = ray2%ix_ele
    ray%ix_source = ray2%ix_ele
    ray%direction = -ray2%direction  ! track backwards
    ray%now = ray%start
    ray%now%vec(1) = ray%start%vec(1) + 1.0e-7 * ray%direction * tan(theta)
    ray%now%vec(5) = ray%start%vec(5) + 1.0e-7 * ray%direction
    ray%track_len = 0
    ray%crossed_end = ray2%crossed_end
    track_len = abs((1 - rr)*ray1%start%vec(5) + &
                           rr*ray2%start%vec(5) - seg%s_mid)
    ! Correct for segment being on the opposite side of the ip.
    if ((seg%s_mid - ray2%start%vec(5)) * ray2%direction < 0) &
                                      track_len = lat%param%total_length - track_len 

    call track_ray_to_wall (ray, lat, walls, hit_flag, track_len)

    ! Do not count if shadowed. 
    ! Because of inaccuracies, test if hit is on opposite side. 
    ! If so then this is not a real shaddow

    if (hit_flag .and. ray%start%vec(1) * ray%now%vec(1) > 0) then
      track_len = 0
      cycle
    endif

    ! frac_illum is the fraction of the segment illuminated by the fan

    frac_illum = 1.0
    if (r_i1 > iss) frac_illum = 1.0 - (r_i1 - iss)
    if (r_i2 < iss + 1) frac_illum = frac_illum - (iss + 1 - r_i2)

    ! calculate the power absorbed by the segment

    call this_seg_power_calc (seg)

  enddo  ! iss

  ! End wall case.
  ! It is assumed that it is impossible for one ray to point to one end wall and 
  ! the other ray to point to the other end wall.

  if (ray1%wall_side == start_side$ .or. ray1%wall_side == exit_side$ .or. &
      ray2%wall_side == start_side$ .or. ray2%wall_side == exit_side$) then

    if (ray1%wall_side == start_side$ .or. ray2%wall_side == start_side$) then
      seg => walls%start_end
      wall_side = start_side$
      n_pos = 0
      n_neg = 0
      dx_seg = walls%positive_x_wall%pt(0)%x - walls%negative_x_wall%pt(0)%x
    else
      seg => walls%exit_end
      wall_side = exit_side$
      n_pos = walls%positive_x_wall%n_pt_tot
      n_neg = walls%negative_x_wall%n_pt_tot
      dx_seg = walls%positive_x_wall%pt(n_pos)%x - walls%negative_x_wall%pt(n_neg)%x
    endif

    if (ray1%wall_side == start_side$ .or. ray1%wall_side == exit_side$) then
      theta_rel = modulo2(theta_floor1 - pi/2, pi)
      track_len = ray1%track_len
      rr = 0
      x1 = ray1%now%vec(1)
    elseif (ray1%wall_side == positive_x$) then
      x1 = walls%positive_x_wall%pt(n_pos)%x
    elseif (ray1%wall_side == negative_x$) then
      x1 = walls%negative_x_wall%pt(n_neg)%x
    endif

    if (ray2%wall_side == start_side$ .or. ray2%wall_side == exit_side$) then
      theta_rel = modulo2(theta_floor1 - pi/2, pi)
      track_len = ray2%track_len
      rr = 1
      x2 = ray2%now%vec(1)
    elseif (ray2%wall_side == positive_x$) then
      x2 = walls%positive_x_wall%pt(n_pos)%x
    elseif (ray2%wall_side == negative_x$) then
      x2 = walls%negative_x_wall%pt(n_neg)%x
    endif

    frac_illum = abs(x2 - x1) / dx_seg    

  endif

enddo  ! fan

!-----------------------------------------------------------------------------
contains

!+
! Perpendicular distance between the ray and the segment.
! Notice that distances can be negative.
!-

function ray_to_seg_distance (ray, seg) result (dist)

type (ray_struct) ray
type (wall_seg_struct) seg
real(rp) dist, now(6)

!

now = ray%now%vec
if (ray%crossed_end) now(5) = now(5) + lat%param%total_length * ray%direction
dx = seg%x_mid - now(1)
ds = seg%s_mid - now(5)
dist = dx * cos(now(2)) - ds * sin(now(2))

end function ray_to_seg_distance

!-----------------------------------------------------------------------------
! contains

subroutine this_seg_power_calc (seg)

type (wall_seg_struct), target :: seg
type (seg_power_struct), pointer :: ep

real(rp) illum_factor, power_per_len

! Only change ep%ix_ele_source and ep%s_source if the contribution to the
! power is larger than what has so far been accumulated.

illum_factor = abs(sin(theta_rel)) * frac_illum / track_len
power_per_len = ((1 - rr)*ray1%p1_factor + rr*ray2%p1_factor) * illum_factor
power%at_wall = power%at_wall + power_per_len * seg%len

ep => seg%power
ep%n_source = ep%n_source + 1

ep%power_per_len  = ep%power_per_len + power_per_len
ep%power_per_area = ep%power_per_area + abs(sin(theta_rel)) * frac_illum * &
        ((1 - rr)*ray1%p2_factor + rr*ray2%p2_factor) / track_len
ep%power_tot = ep%power_tot + power_per_len * seg%len

if (power_per_len > ep%main_source%power_per_len) then
  ep%main_source%power_per_len = power_per_len
  ep%main_source%s             = ray1%start%vec(5) * (1 - rr) + ray2%start%vec(5) * rr
  ep%main_source%ix_ele        = ray2%ix_source
endif

! from Sands and Chao/Tigner p188
! Flux is 3.248 * Power / critical photon energy 
! 3.248 is 15*sqrt(3)/8
! Also: Convert keV to Joules 

ep%photons_per_sec = ep%photons_per_sec + 3.248 * illum_factor * &
                       power_factor * seg%len / (2.218 * (energy/1e9)**3 * 1.602e-16) 

end subroutine this_seg_power_calc

end subroutine seg_power_calc
