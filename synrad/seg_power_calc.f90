!+
! subroutine seg_power_calc (rays, i_ray, inside, outside, 
!                             lat, gen, power)
!
! subroutine to calculate the synch radiation power for
!      segments of the wall from one lat element
!
! Modules needed:
!   use sr_mod
!
! Input:
!   rays(*) -- ray_struct: array of rays from one lat element
!   i_ray   -- integer: index of highest ray used
!   lat    -- lat_struct: with twiss propagated and mat6s made
!   inside  -- wall_struct: inside wall with outline ready
!   outside -- wall_struct: outside wall with outline ready
!   gen    -- synrad_param_struct: Contains lat name,
!                     vert emittance, and beam current
!
! Output:
!   power(*)  -- ele_power_struct: power radiated from a lat ele
!   inside  -- wall_struct: inside wall with power information
!   outside -- wall_struct: outside wall with power information
!-

subroutine seg_power_calc (rays, i_ray, inside, outside, lat, gen, power)

  use sr_struct
  use sr_interface

  implicit none

  type (ray_struct) :: rays(:)
  type (wall_struct) inside, outside
  type (wall_struct), pointer :: wall
  type (ray_struct) :: ray
  type (synrad_param_struct) gen
  type (lat_struct) lat
  type (sr_power_struct), pointer :: ep
  type (ele_power_struct) power

  real(rp) energy
  real(rp) sig_y, sig_yp, sig_y_eff, dx, ds, theta_rel
  real(rp) rr, dr2, dr1, theta, factor
  real(rp) dx_wall, ds_wall, r1, r2, dpower, track_len
  real(rp) r_i1, r_i2, frac_illum, gamma
  real(rp) theta_base, theta_floor1, theta_floor2, theta_floor_seg
  integer i_ray, iw, i, ip, is, iss
  integer type0, type1, i1_pt, i2_pt

  logical hit_flag

  type (source_struct), pointer :: temp_rays(:)
  integer temp_size

! calculate power radiated


  ! get energy from source element
  energy = lat%ele(rays(1)%ix_source)%value(E_TOT$)


  ! 14.1e3 [m (eV)^-3] (used below) is  Cgamma * 1e9 / (2 pi) 
  !    Cgamma is from Sands (p98) which is 8.85e-5 m (GeV)^-3 
  !    for electrons.  Not valid for protons, etc!!!
  factor = 14.1e3 * gen%i_beam * (energy/1e9)**4
  do i = 2, i_ray

    ! sum up the power radiated by the element:
    ! above factor * delta S * avg(gbend^2)
    ! which is an integrated segment of Sands eq. 4.5
    ! * the current to get the power
    power%radiated = power%radiated + factor * &
                   (rays(i)%start%vec(5) - rays(i-1)%start%vec(5)) * &
                   (rays(i)%g_bend**2 + rays(i-1)%g_bend**2) / 2
  enddo
 
! Let user determine cutoff from all data instead of 1 W cutoff:
!  if (power%radiated < 1) return

! calculate the power factors for each ray 

  wall => rays(1)%wall

  do i = 1, i_ray
    iw = rays(i)%ix_wall_pt
    

    ! should update to include energy and dispersion component!!!
    ! sig_y is the sqrt of the vert emitt * vert beta
    sig_y = sqrt(gen%epsilon_y * rays(i)%y_twiss%beta)

    ! sig_yprime is the vertical opening angle

    call convert_total_energy_to (energy, lat%param%particle, gamma)
    sig_yp = sqrt(gen%epsilon_y * rays(i)%y_twiss%gamma + 1 / gamma**2)

    ! the effective sig_y includes the increase from sig_yp
    sig_y_eff = sqrt(sig_y**2 + (rays(i)%track_len * sig_yp)**2)

    ! the p1 factor is used in power/length
    rays(i)%p1_factor = factor * rays(i)%g_bend

    ! the p2 factor is used in power/area
    rays(i)%p2_factor = rays(i)%p1_factor / (sqrt(twopi) * sig_y_eff)

  enddo

! fill the segment array by interpolating between fan rays.
! theta_rel is the relative angle between the ray and the wall
! Note: the power at a wall segment could be due to more than one source.
! The code tries to set %ix_ele_source and %s_source to reflect the
! largest source

  do i = 2, i_ray

! loop over all segments whose mid-point falls between the hit points
! of rays(i-1) and rays(i). We must take into account rays that straddle
! the IP

    i1_pt = rays(i-1)%ix_wall_pt
    r1 = wall%pt(i1_pt)%ix_seg + rays(i-1)%r_wall * wall%pt(i1_pt)%n_seg + 1

    i2_pt = rays(i)%ix_wall_pt
    r2 = wall%pt(i2_pt)%ix_seg + rays(i)%r_wall * wall%pt(i2_pt)%n_seg + 1

    if (rays(i-1)%crossed_end .xor. rays(i)%crossed_end) then
      r_i1 = max(r1, r2)
      r_i2 = min(r1, r2) + wall%n_seg_tot
    else
      r_i1 = min(r1, r2)
      r_i2 = max(r1, r2)
    endif

    do iss = int(r_i1), int(r_i2)

      is = mod(iss-1, wall%n_seg_tot) + 1
      ip = wall%seg(is)%ix_pt

      theta_base = theta_floor(rays(i-1)%now%vec(5), lat)
      theta_floor1 = rays(i-1)%now%vec(2) + theta_base
      theta_floor2 = rays(i)%now%vec(2) + &
                                theta_floor(rays(i)%now%vec(5), lat, theta_base)

      dx = wall%seg(is)%x_mid - rays(i-1)%now%vec(1)
      ds = wall%seg(is)%s_mid - rays(i-1)%now%vec(5)
      if (ds >  lat%param%total_length/2) ds = ds - lat%param%total_length
      if (ds < -lat%param%total_length/2) ds = ds + lat%param%total_length
      dr1 = dx * cos(rays(i-1)%now%vec(2)) - ds * sin(rays(i-1)%now%vec(2))

      dx = wall%seg(is)%x_mid - rays(i)%now%vec(1)
      ds = wall%seg(is)%s_mid - rays(i)%now%vec(5)
      if (ds >  lat%param%total_length/2) ds = ds - lat%param%total_length
      if (ds < -lat%param%total_length/2) ds = ds + lat%param%total_length
      dr2 = dx * cos(rays(i)%now%vec(2)) - ds * sin(rays(i)%now%vec(2))

      if ((dr1 - dr2) == 0) then
        rr = 1
      else
        rr = dr1 / (dr1 - dr2)
      endif

      if (rr > 1) rr = 1
      if (rr < 0) rr = 0

      ! In a bend, the angle w.r.t. the center line can change rapidly with s.
      ! For the interpolation it is therefore better to use the floor angle.
      ! Since increasing floor angle is opposite the bend angle, there is 
      ! a plus sign.

      theta_floor_seg = (1 - rr) * theta_floor1 + rr * theta_floor2
      theta = theta_floor_seg - theta_floor(wall%seg(is)%s_mid, lat, theta_base)

      dx_wall = wall%pt(ip)%x - wall%pt(ip-1)%x
      ds_wall = wall%pt(ip)%s - wall%pt(ip-1)%s
      theta_rel = theta - atan2(dx_wall, ds_wall)
      if (theta_rel > pi) theta_rel = theta_rel - twopi
      if (theta_rel < -pi) theta_rel = theta_rel + twopi

! Track backwards to see if there is a shadow.
! If so cycle before we calculate the power

      if (wall%side == outside$ .and. theta_rel*rays(i)%direction < 0) cycle
      if (wall%side == inside$ .and. theta_rel*rays(i)%direction > 0) cycle
!      if (rr < 0 .or. rr > 1) cycle  ! not in fan

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

      ray%start%vec(1) = wall%seg(is)%x_mid
      ray%start%vec(5) = wall%seg(is)%s_mid
      ray%start%vec(2) = theta
      ray%ix_ele = rays(i)%ix_ele
      ray%direction = -rays(i)%direction
      ray%now = ray%start
      ray%now%vec(1) = ray%start%vec(1) + 2.0e-4 * ray%direction * tan(theta)
      ray%now%vec(5) = ray%start%vec(5) + 2.0e-4 * ray%direction
      ray%track_len = 0
      ray%crossed_end = rays(i)%crossed_end
      track_len = abs((1 - rr)*rays(i-1)%start%vec(5) + &
                             rr*rays(i)%start%vec(5) - wall%seg(is)%s_mid)
      ! We assume that the travel length cannot be greater then half the circumference.
      if (abs(track_len) > lat%param%total_length / 2) track_len = lat%param%total_length - track_len 
      call track_ray_to_wall (ray, lat, inside, outside, hit_flag, track_len)
      ! Do not count if shadowed. 
      ! Because of inaccuracies test if hit is on opposite side. 
      ! If so then this is not a real shaddow
      if (hit_flag .and. ray%start%vec(1) * ray%now%vec(1) > 0) then
        track_len = 0
        cycle
      endif

      ep => wall%seg(is)%sr_power

! frac_illum is the fraction of the segment illuminated by the fan

      frac_illum = 1.0
      if (r_i1 > iss) frac_illum = 1.0 - (r_i1 - iss)
      if (r_i2 < iss + 1) frac_illum = frac_illum - (iss + 1 - r_i2)

! only change ep%ix_ele_source and ep%s_source if the contribution to the
! power is larger than what has so far been accumulated.

      dpower = ((1 - rr)*rays(i-1)%p1_factor + rr*rays(i)%p1_factor) * &
                        abs(sin(theta_rel)) * frac_illum / track_len

      ep%power_per_len  = ep%power_per_len + dpower
      ep%power_per_area = ep%power_per_area + &
            abs(sin(theta_rel)) * frac_illum * &
            ((1 - rr)*rays(i-1)%p2_factor + rr*rays(i)%p2_factor) / track_len
      ep%power = ep%power + dpower * wall%seg(is)%len

      ! from Sands and Chao/Tigner p188
      ! Flux is 3.248 * Power / critical photon energy 
      ! 3.248 is 15*sqrt(3)/8
      ep%photons_per_sec = ep%photons_per_sec + &
            3.248 * (dpower * wall%seg(is)%len) / &
            (2.218 * (energy/1e9)**3 * &
            ( (1 - rr) * rays(i)%g_bend + rr * rays(i-1)%g_bend) &
            * 1.602e-16 ) ! Convert keV to Joules 

      ep%n_source = ep%n_source + 1

      if (associated(ep%rays)) then
        temp_size = size(ep%rays)
        if (temp_size < ep%n_source) then
          allocate(temp_rays(temp_size))
          temp_rays = ep%rays
          deallocate(ep%rays)
          allocate (ep%rays(ep%n_source))
          ep%rays(1:temp_size) = temp_rays(1:temp_size)
          deallocate(temp_rays)
        endif
      else
        allocate (ep%rays(5))
      endif

      power%at_wall = power%at_wall + dpower * wall%seg(is)%len
      ep%rays(ep%n_source)%start%vec = rays(i-1)%start%vec * &
                                         (1 - rr) + rays(i)%start%vec * rr
      ep%rays(ep%n_source)%now%vec = rays(i-1)%now%vec * &
                                         (1 - rr) + rays(i)%now%vec * rr
      ep%rays(ep%n_source)%ix_ele_source = rays(i)%ix_source
      if (dpower > ep%power_per_len) then
        ep%ix_ele_source = rays(i)%ix_source
        ep%s_source = ep%rays(ep%n_source)%start%vec(5)
      endif
    enddo

  enddo

end subroutine
