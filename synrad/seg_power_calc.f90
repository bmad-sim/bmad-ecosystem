!+
! subroutine seg_power_calc (fan, i_ray, walls, lat, gen, power)
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
!   gen    -- synrad_param_struct: Contains lat name, vert emittance, and beam current.
!
! Output:
!   power(*)  -- ele_power_struct: power radiated from a lat ele
!   walls   -- wall_struct: both walls and ends with power information
!-

subroutine seg_power_calc (fan, i_ray, walls, lat, gen, power)

use synrad_struct
use synrad_interface, except => seg_power_calc

implicit none

type (ray_struct) :: fan(:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: negative_x_wall, positive_x_wall
type (wall_struct), pointer :: wall
type (ray_struct) :: ray
type (synrad_param_struct) gen
type (lat_struct) lat
type (wall_seg_struct), pointer :: ws
type (seg_power_struct), pointer :: ep
type (ele_power_struct) power

real(rp) energy
real(rp) sig_y, sig_yp, sig_y_eff, dx, ds, theta_rel
real(rp) rr, dr2, dr1, theta, factor
real(rp) dx_wall, ds_wall, r1, r2, dpower, track_len
real(rp) r_i1, r_i2, frac_illum, gamma, end_s_scale, end_x_scale
real(rp) theta_base, theta_floor1, theta_floor2, theta_floor_seg
integer i_ray, iw, i, ip, is, iss
integer type0, type1, i1_pt, i2_pt

logical hit_flag

type (source_struct), pointer :: temp_sources(:)
integer temp_size

! set pointers
positive_x_wall => walls%positive_x_wall
negative_x_wall => walls%negative_x_wall


! calculate power radiated

! get energy from source element
energy = lat%ele(fan(1)%ix_source)%value(E_TOT$)


! 14.1e3 [m (eV)^-3] (used below) is  Cgamma * 1e9 / (2 pi) 
!    Cgamma is from Sands (p98) which is 8.85e-5 m (GeV)^-3 
!    for electrons.  Not valid for protons, etc!!!

factor = 14.1e3 * gen%i_beam * (energy/1e9)**4
do i = 2, i_ray

  ! Sum up the power radiated by the element:
  ! above factor * delta S * avg(gbend^2)
  ! which is an integrated segment of Sands eq. 4.5
  ! * the current to get the power
  power%radiated = power%radiated + factor * &
                 (fan(i)%start%vec(5) - fan(i-1)%start%vec(5)) * &
                 (fan(i)%g_bend**2 + fan(i-1)%g_bend**2) / 2
enddo
 
! Ignore element if the radiated power is 0
if (power%radiated == 0) return


! wall is the side that is being hit by the radiation.

wall => fan(1)%wall

! calculate the power factors for each ray 
do i = 1, i_ray
  iw = fan(i)%ix_wall_pt
  

  ! should update to include energy and dispersion component!!!
  ! sig_y is the sqrt of the vert emitt * vert beta
  sig_y = sqrt(gen%epsilon_y * fan(i)%y_twiss%beta)

  ! sig_yprime is the vertical opening angle

  call convert_total_energy_to (energy, lat%param%particle, gamma)
  sig_yp = sqrt(gen%epsilon_y * fan(i)%y_twiss%gamma + 1 / gamma**2)

  ! the effective sig_y includes the increase from sig_yp
  sig_y_eff = sqrt(sig_y**2 + (fan(i)%track_len * sig_yp)**2)

  ! the p1 factor is used in power/length
  fan(i)%p1_factor = factor * fan(i)%g_bend

  ! the p2 factor is used in power/area
  fan(i)%p2_factor = fan(i)%p1_factor / (sqrt(twopi) * sig_y_eff)

enddo

! fill the segment array by interpolating between the fan rays.
! theta_rel is the relative angle between the ray and the wall
! Note: the power at a wall segment could be due to more than one source.
! The code tries to set %ix_ele_source and %s_source to reflect the
! largest source

do i = 2, i_ray

  ! loop over all wall segments whose mid-point falls between the hit points
  ! of rays fan(i-1) and fan(i). We must take into account rays that straddle
  ! the lattice ends for circular lattices.

  i1_pt = fan(i-1)%ix_wall_pt
  r1 = wall%pt(i1_pt)%ix_seg + fan(i-1)%r_wall * wall%pt(i1_pt)%n_seg + 1

  i2_pt = fan(i)%ix_wall_pt
  r2 = wall%pt(i2_pt)%ix_seg + fan(i)%r_wall * wall%pt(i2_pt)%n_seg + 1

  ! There is a complication if, for a circular lattice, the fan is 
  ! illuminating the lattice ends.

  if (fan(i-1)%crossed_end .neqv. fan(i)%crossed_end) then
    r_i1 = max(r1, r2)
    r_i2 = min(r1, r2) + wall%n_seg_tot   ! This addition will be "mod"ed out below
  else
    r_i1 = min(r1, r2)
    r_i2 = max(r1, r2)
  endif

  ! Loop over all segments that have light hitting it.

  do iss = int(r_i1), int(r_i2)

    is = mod(iss-1, wall%n_seg_tot) + 1  ! segment index
    ip = wall%seg(is)%ix_pt              ! wall point associated with this segment

    ! theta_floor1 and theta_floor2 are the ray angles in the floor coordinate system.

    theta_base = theta_floor(fan(i-1)%now%vec(5), lat) 
    theta_floor1 = fan(i-1)%now%vec(2) + theta_base
    theta_floor2 = fan(i)%now%vec(2) + &
                              theta_floor(fan(i)%now%vec(5), lat, theta_base)

    ! fan(j)%now is the point where the j^th ray hits the wall.
    ! dx, ds are the coordinates of the segment midpoint w.r.t the fan hit point
    ! dr1 is the perpendicular distance between the fan(i-1) ray and a parallel ray
    !  that is aimed to strike the segment midpoint.

    dx = wall%seg(is)%x_mid - fan(i-1)%now%vec(1)
    ds = wall%seg(is)%s_mid - fan(i-1)%now%vec(5)
    if (ds >  lat%param%total_length/2) ds = ds - lat%param%total_length
    if (ds < -lat%param%total_length/2) ds = ds + lat%param%total_length
    dr1 = dx * cos(fan(i-1)%now%vec(2)) - ds * sin(fan(i-1)%now%vec(2))

    ! dr2 is the perpendicular distance between the fan(i) ray and a parallel ray
    !  that is aimed to strike the segment midpoint.
    ! Notice that dr1 and/or dr2 can be negative.

    dx = wall%seg(is)%x_mid - fan(i)%now%vec(1)
    ds = wall%seg(is)%s_mid - fan(i)%now%vec(5)
    if (ds >  lat%param%total_length/2) ds = ds - lat%param%total_length
    if (ds < -lat%param%total_length/2) ds = ds + lat%param%total_length
    dr2 = dx * cos(fan(i)%now%vec(2)) - ds * sin(fan(i)%now%vec(2))

    ! theta_floor_seg is the ray angle averaged over theta_floor1 and theta_floor2.
    ! The averaging factor is rr.

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

    theta_floor_seg = (1 - rr) * theta_floor1 + rr * theta_floor2
    theta = theta_floor_seg - theta_floor(wall%seg(is)%s_mid, lat, theta_base)

    ! theta_rel gives the ray angle with respect to the wall surface.
    ! theta_rel = 0 means the ray and wall surface are parallel.

    dx_wall = wall%pt(ip)%x - wall%pt(ip-1)%x
    ds_wall = wall%pt(ip)%s - wall%pt(ip-1)%s
    theta_rel = theta - atan2(dx_wall, ds_wall)
    if (theta_rel > pi) theta_rel = theta_rel - twopi
    if (theta_rel < -pi) theta_rel = theta_rel + twopi

    ! A given wall segment may be shadowed by a part of the wall upstream so we need
    !   to check if this is the case. 
    ! First, if the ray is stricking the wall from the back side then this segment 
    !   is shadowed.

    if (wall%side == positive_x$ .and. theta_rel*fan(i)%direction < 0) cycle
    if (wall%side == negative_x$ .and. theta_rel*fan(i)%direction > 0) cycle

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

    ray%start%vec(1) = wall%seg(is)%x_mid
    ray%start%vec(5) = wall%seg(is)%s_mid
    ray%start%vec(2) = theta
    ray%ix_ele = fan(i)%ix_ele
    ray%ix_source = fan(i)%ix_ele
    ray%direction = -fan(i)%direction  ! track backwards
    ray%now = ray%start
    ray%now%vec(1) = ray%start%vec(1) + 1.0e-7 * ray%direction * tan(theta)
    ray%now%vec(5) = ray%start%vec(5) + 1.0e-7 * ray%direction
    ray%track_len = 0
    ray%crossed_end = fan(i)%crossed_end
    track_len = abs((1 - rr)*fan(i-1)%start%vec(5) + &
                           rr*fan(i)%start%vec(5) - wall%seg(is)%s_mid)
    ! We assume that the travel length cannot be greater then half the circumference.
    if (abs(track_len) > lat%param%total_length / 2) &
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

    ! if the lattice is not circular, then rays that go past the end of the 
    ! lattice should not wrap around

    ws => wall%seg(is)
    if (lat%param%lattice_type == linear_lattice$) then
      if ( (fan(i)%direction == 1) .and. (iss > wall%n_seg_tot) ) then
        ws =>walls%final_end
      else if ( (fan(i)%direction == -1) .and. (iss <= wall%n_seg_tot) ) then
        ws =>walls%initial_end
      endif

      if (  ((fan(i)%direction == 1) .and. (iss > wall%n_seg_tot)) .or.  &
             ((fan(i)%direction == -1) .and. (iss <= wall%n_seg_tot))  ) then
        end_s_scale = (ws%s_mid - fan(i)%start%vec(5)) / &
                (ws%s_mid - fan(i)%start%vec(5) + fan(i)%now%vec(5))
        end_x_scale = end_s_scale * (fan(i)%now%vec(1) - fan(i)%start%vec(1)) + fan(i)%start%vec(1)
        theta_rel = pi/2.0 + &
                atan2( end_x_scale, ws%s_mid - fan(i)%start%vec(5) + fan(i)%now%vec(5) )
        if (end_x_scale>0) then 
          end_x_scale = walls%positive_x_wall%pt(walls%positive_x_wall%n_pt_tot)%x - end_x_scale
        else
          end_x_scale = ws%x - end_x_scale
        endif 

        frac_illum = abs(end_x_scale) / ws%len
      endif
    endif

    ep => ws%power

    ! only change ep%ix_ele_source and ep%s_source if the contribution to the
    ! power is larger than what has so far been accumulated.

    dpower = ((1 - rr)*fan(i-1)%p1_factor + rr*fan(i)%p1_factor) * &
                        abs(sin(theta_rel)) * frac_illum / track_len

    ep%power_per_len  = ep%power_per_len + dpower
    ep%power_per_area = ep%power_per_area + abs(sin(theta_rel)) * frac_illum * &
            ((1 - rr)*fan(i-1)%p2_factor + rr*fan(i)%p2_factor) / track_len
    ep%power_tot = ep%power_tot + dpower * ws%len

    ! from Sands and Chao/Tigner p188
    ! Flux is 3.248 * Power / critical photon energy 
    ! 3.248 is 15*sqrt(3)/8

    ep%photons_per_sec = ep%photons_per_sec + 3.248 * (dpower * ws%len) / &
            (2.218 * (energy/1e9)**3 * ( (1 - rr) * fan(i)%g_bend + rr * fan(i-1)%g_bend) &
            * 1.602e-16 ) ! Convert keV to Joules 

    ep%n_source = ep%n_source + 1

    if (associated(ep%sources)) then
      temp_size = size(ep%sources)
      if (temp_size < ep%n_source) then
        allocate(temp_sources(temp_size))
        temp_sources = ep%sources
        deallocate(ep%sources)
        allocate (ep%sources(ep%n_source))
        ep%sources(1:temp_size) = temp_sources(1:temp_size)
        deallocate(temp_sources)
      endif
    else
      allocate (ep%sources(5))
    endif

    power%at_wall = power%at_wall + dpower * ws%len
    ep%sources(ep%n_source)%start%vec = fan(i-1)%start%vec * &
                                         (1 - rr) + fan(i)%start%vec * rr
    ep%sources(ep%n_source)%now%vec = fan(i-1)%now%vec * &
                                         (1 - rr) + fan(i)%now%vec * rr
    ep%sources(ep%n_source)%ix_ele_source = fan(i)%ix_source
    if (dpower > ep%power_per_len) then
      ep%ix_ele_source = fan(i)%ix_source
      ep%s_source = ep%sources(ep%n_source)%start%vec(5)
    endif
  enddo

enddo

end subroutine seg_power_calc
