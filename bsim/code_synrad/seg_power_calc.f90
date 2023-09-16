!+
! subroutine seg_power_calc (fan, i_ray, walls, wall_side, branch, gen, ele_power)
!
! Subroutine to calculate the synch radiation power hitting a wall.
! It is assumed that the radiation is hitting only one wall.
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   fan(:)    -- ray_struct: array of rays from one lattice element
!   i_ray     -- integer: index of highest ray used
!   branch    -- branch_struct: with twiss propagated and mat6s made
!   walls     -- walls_struct: Both walls and ends
!   wall_side -- Integer: +x or -x side where rays are hitting.
!   sr_param  -- synrad_param_struct: Contains lattice name, vert emittance, and beam current.
!
! Output:
!   power(*)  -- ele_power_struct: power radiated from a lattice ele
!   walls     -- wall_struct: both walls and ends with power information
!-

subroutine seg_power_calc (fan, i_ray, walls, wall_side, branch, sr_param, ele_power)

use synrad_struct
use synrad_interface, except => seg_power_calc

implicit none

type (ray_struct), target :: fan(:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: wall
type (wall_seg_struct), pointer :: seg
type (ray_struct) :: ray
type (ray_struct), pointer :: ray1, ray2
type (synrad_param_struct) sr_param
type (branch_struct) branch
type (ele_power_struct) ele_power

real(rp) energy, dx_seg, x1, x2
real(rp) sig_y, sig_yp, sig_y_eff, dx, ds, theta_rel
real(rp) rr, dr2, dr1, theta, power_factor
real(rp) dx_wall, ds_wall, track_len
real(rp) frac_illum, gamma, len1, len2
real(rp) theta_ray1, theta_ray2, theta_ray, dr_end_pt(4)

integer i_ray, iw, i, ip, is, iss, wall_side, ix_seg1, ix_seg2
integer type0, type1, i_pt, n_pos, n_neg, ix_seg_start, ix_seg_end

! get energy from source element

energy = branch%ele(fan(1)%start%ix_ele)%value(e_tot$)

! 14.1e3 [m (eV)^-3] (used below) is  Cgamma * 1e9 / (2 pi) 
!    Cgamma is from Sands (p98) which is 8.85e-5 m (GeV)^-3 for electrons.  

power_factor = 14.1e3 * sr_param%i_beam * (energy/1e9)**4 * (mass_of(electron$) / mass_of(branch%param%particle))**3

do i = 2, i_ray

  ray1 => fan(i-1)
  ray2 => fan(i)

  ! Sum up the power radiated by the element:
  ! above power_factor * delta S * avg(gbend^2)
  ! which is an integrated segment of Sands eq. 4.5
  ! * the current to get the power

  ele_power%radiated = ele_power%radiated + power_factor * &
         (ray2%start%s - ray1%start%s) * (ray2%g_bend**2 + ray1%g_bend**2) / 2
enddo
 
! Ignore element if the radiated power is 0
if (ele_power%radiated == 0) return


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
  sig_y = sqrt(sr_param%epsilon_y * ray2%y_twiss%beta)

  ! sig_yprime is the vertical opening angle

  call convert_total_energy_to (energy, branch%param%particle, gamma)
  sig_yp = sqrt(sr_param%epsilon_y * ray2%y_twiss%gamma + 1 / gamma**2)

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
    ix_seg1 = 0
  elseif (ray1%wall_side == exit_side$) then
    ix_seg1 = wall%n_seg_max
  else
    ix_seg1 = ray1%ix_seg_pt
  endif

  if (ray2%wall_side == start_side$) then
    ix_seg2 = 0
  elseif (ray2%wall_side == exit_side$) then
    ix_seg2 = wall%n_seg_max
  else
    ix_seg2 = ray2%ix_seg_pt
  endif

  ! There is a complication if, for a circular lattice, the fan is 
  ! illuminating the lattice ends.

  if (walls%lat_geometry == closed$ .and. abs(ix_seg1 - ix_seg2) > wall%n_seg_max/2) then
    ix_seg_start = max(ix_seg1, ix_seg2)
    ix_seg_end = min(ix_seg1, ix_seg2) + wall%n_seg_max   ! This addition will be "mod"ed out below
  else
    ix_seg_start = min(ix_seg1, ix_seg2)
    ix_seg_end = max(ix_seg1, ix_seg2)
  endif

  ! theta_ray1 and theta_ray2 are the ray angles in the floor coordinate system.
  ! Make sure |theta_ray1 - theta_ray2| < pi

  theta_ray1 = ray1%now_floor%theta
  theta_ray2 = modulo2(ray2%now_floor%theta - theta_ray1, pi) + theta_ray1

  ! Loop over all segments that have light hitting it.

  do iss = ix_seg_start, ix_seg_end

    if (branch%param%geometry == open$ .and. iss > wall%n_seg_max) exit

    is = iss
    if (is > wall%n_seg_max) is = is - wall%n_seg_max
    seg => wall%seg(is)

    ip = seg%ix_pt              ! wall point associated with this segment

    ! dr1 is the perpendicular distance between the ray1 ray and the segment mid-point.
    ! dr2 is the perpendicular distance between the ray2 ray and the segment mid-point.
    ! Notice that typically dr1 and dr2 have opposite signs.

    dr1 = ray_to_seg_perp_distance (ray1, seg, dr_end_pt(1:2))
    dr2 = ray_to_seg_perp_distance (ray2, seg, dr_end_pt(3:4))

    ! If all the dr_end_pt values are positive or negative then the segment is outside the fan.

    if (minval(dr_end_pt) >= 0 .or. maxval(dr_end_pt) <= 0) cycle

    if ((dr1 - dr2) == 0) then ! If both rays hit the same spot...
      rr = 1  ! then just use theta_ray2
    else
      rr = dr1 / (dr1 - dr2)
    endif

    ! theta_ray is the ray angle averaged over theta_ray1 and theta_ray2.

    theta_ray = (1 - rr) * theta_ray1 + rr * theta_ray2

    len1 = sqrt((ray1%start_floor%r(1) - seg%r_floor_mid(1))**2 + (ray1%start_floor%r(3) - seg%r_floor_mid(3))**2)
    len2 = sqrt((ray2%start_floor%r(1) - seg%r_floor_mid(1))**2 + (ray2%start_floor%r(3) - seg%r_floor_mid(3))**2)
    track_len = (1 - rr) * len1 + rr * len2

    ! theta_rel gives the graze angle of the ray with respect to the wall surface.
    ! theta_rel = 0 means the ray and wall surface are parallel.

    theta_rel = modulo2(theta_ray - seg%theta, pi)

    ! A given wall segment may be shadowed by a part of the wall upstream so we need
    !   to check if this is the case. 
    ! First, if the ray is stricking the wall from the back side then this segment 
    !   is shadowed.

    if (wall%side == positive_x$ .and. theta_rel < 0) cycle
    if (wall%side == negative_x$ .and. theta_rel > 0) cycle

    ! Track backwards from the segment midpoint to see if this segment is shadowed.
    ! If a ray can travel backwards a distance equal to the distance the rays 
    ! traveled forward without hitting a wall then there is no shadow.

    ray%direction         = -ray2%direction  ! track backwards
    !! ray%start%vec(1)      = seg%x_mid
    !! ray%start%s           = seg%s_mid
    ray%start_floor%r     = seg%r_floor_mid
    ray%start_floor%theta = modulo2(theta_ray + pi, pi)
    ray%now = ray%start

    call track_ray_to_wall (ray, walls)

    ! Do not count if shadowed. 
    ! Because of inaccuracies, test if hit is on opposite side. 
    ! If so then this is not a real shaddow

    if (ray%track_len < track_len .and. ray%wall_side == ray1%wall_side) cycle

    ! frac_illum is the fraction of the segment illuminated by the fan

    frac_illum = 1.0

    if (is == ray1%ix_seg_pt .and. is == ray2%ix_seg_pt) then
      frac_illum = abs(ray2%r_seg - ray1%r_seg)

    elseif (is == ray1%ix_seg_pt) then
      if (ix_seg1 == ix_seg_start) then
        frac_illum = 1 - ray1%r_seg
      else
        frac_illum = ray1%r_seg
      endif

    elseif (is == ray2%ix_seg_pt) then
      if (ix_seg2 == ix_seg_start) then
        frac_illum = 1 - ray2%r_seg
      else
        frac_illum = ray2%r_seg
      endif
    endif

    ! Calculate the power absorbed by the segment

    call this_seg_power_calc (seg, ele_power, frac_illum)

    ! For custom calculations.

    ray%start     = ray1%start
    ray%start%s   = (1 - rr) * ray1%start%s + rr * ray2%start%s
    ray%start%vec = (1 - rr) * ray1%start%vec + rr * ray2%start%vec

    ray%now     = ray1%now
    ray%now%s   = (1 - rr) * ray1%now%s + rr * ray2%now%s
    ray%now%vec = (1 - rr) * ray1%now%vec + rr * ray2%now%vec

    ray%direction = ray1%direction
    ray%track_len = track_len

    ray%x_twiss = average_twiss(1-rr, ray1%x_twiss, ray2%x_twiss)
    ray%y_twiss = average_twiss(1-rr, ray1%y_twiss, ray2%y_twiss)

    ray%g_bend = (1 - rr) * ray1%g_bend + rr * ray2%g_bend

    call synrad_custom_seg_calc (wall, ray, seg, frac_illum)

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
      n_pos = walls%positive_x_wall%n_pt_max
      n_neg = walls%negative_x_wall%n_pt_max
      dx_seg = walls%positive_x_wall%pt(n_pos)%x - walls%negative_x_wall%pt(n_neg)%x
    endif

    if (ray1%wall_side == start_side$ .or. ray1%wall_side == exit_side$) then
      theta_rel = modulo2(theta_ray1 - pi/2, pi)
      track_len = ray1%track_len
      rr = 0
      x1 = ray1%now%vec(1)
    elseif (ray1%wall_side == positive_x$) then
      x1 = walls%positive_x_wall%pt(n_pos)%x
    elseif (ray1%wall_side == negative_x$) then
      x1 = walls%negative_x_wall%pt(n_neg)%x
    endif

    if (ray2%wall_side == start_side$ .or. ray2%wall_side == exit_side$) then
      theta_rel = modulo2(theta_ray1 - pi/2, pi)
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

function ray_to_seg_perp_distance (ray, seg, dist_at_ends) result (dist)

type (ray_struct) ray
type (wall_seg_struct) seg
real(rp) dist, t, dist_at_ends(2), r_floor(3)

! Calc at segment center

t = ray%now_floor%theta
dist = cos(t) * (seg%r_floor_mid(1) - ray%now_floor%r(1)) - sin(t) * (seg%r_floor_mid(3) - ray%now_floor%r(3))

! calc distance at segment ends

dist_at_ends(2) = cos(t) * (seg%r_floor(1) - ray%now_floor%r(1)) - sin(t) * (seg%r_floor(3) - ray%now_floor%r(3))
r_floor = 2*seg%r_floor_mid - seg%r_floor  ! r at beginning of segment
dist_at_ends(1) = cos(t) * (r_floor(1) - ray%now_floor%r(1)) - sin(t) * (r_floor(3) - ray%now_floor%r(3))

end function ray_to_seg_perp_distance

!-----------------------------------------------------------------------------
! contains

subroutine this_seg_power_calc (seg, ele_power, frac_illum)

type (wall_seg_struct), target :: seg
type (seg_power_struct), pointer :: ep
type (ele_power_struct) ele_power

real(rp) illum_factor, power_per_len, frac_illum

if (sr_param%filter_phantom_photons .and. wall%pt(seg%ix_pt)%phantom) return

! Only change ep%ix_ele_source and ep%s_source if the contribution to the
! power is larger than what has so far been accumulated.

illum_factor = abs(sin(theta_rel)) * frac_illum / track_len
power_per_len = ((1 - rr)*ray1%p1_factor + rr*ray2%p1_factor) * illum_factor

ele_power%at_wall = ele_power%at_wall + power_per_len * seg%len

ep => seg%power
ep%n_source = ep%n_source + 1

ep%power_per_len  = ep%power_per_len + power_per_len
ep%power_per_area = ep%power_per_area + abs(sin(theta_rel)) * frac_illum * &
        ((1 - rr)*ray1%p2_factor + rr*ray2%p2_factor) / track_len
ep%power_tot = ep%power_tot + power_per_len * seg%len

if (power_per_len > ep%main_source%power_per_len) then
  ep%main_source%power_per_len = power_per_len
  ep%main_source%s             = ray1%start%s * (1 - rr) + ray2%start%s * rr
  ep%main_source%ix_ele        = ray2%start%ix_ele
endif

! from Sands and Chao/Tigner p188
! Flux is 3.248 * Power / critical photon energy 
! 3.248 is 15*sqrt(3)/8
! Also: Convert keV to Joules 

ep%photons_per_sec = ep%photons_per_sec + 3.248 * illum_factor * &
                       power_factor * seg%len / (2.218 * (energy/1e9)**3 * 1.602e-16) 

end subroutine this_seg_power_calc

end subroutine seg_power_calc
