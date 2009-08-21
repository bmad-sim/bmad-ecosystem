!+
! subroutine ele_synrad_power (lat, ie, orb, direction, power, walls, gen)
!
! subroutine to calculate the synch radiation power from
!      one element of the lat
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   lat     -- lat_struct: with twiss propagated and mat6s made
!   ie      -- integer: index of the current element in lat
!   orb(0:) -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 for in direction of s
!                         -1 for against s
!   walls  -- walls_struct: both walls with outline ready
!   gen    -- synrad_param_struct: Contains lat name,
!                     vert emittance, and beam current
!
! Output:
!   power(:) -- ele_power_struct: power radiated from a lat ele
!   wall   -- walls_struct: both walls with power information
!-

subroutine ele_synrad_power (lat, ie, orb, direction, power, walls, gen)

  use synrad_struct
  use synrad_interface, except => ele_synrad_power

  implicit none

  type (lat_struct), target :: lat
  type (coord_struct) orb(0:)
  type (ele_struct), pointer :: ele
  type (walls_struct), target :: walls
  type (wall_struct), pointer :: negative_x_wall, positive_x_wall
  type (ray_struct), allocatable, save :: fan(:)
  type (ray_struct) ray, ray_temp
  type (synrad_param_struct) gen
  type (ele_power_struct) power(:)

  integer direction, ie, i_ray, n_slice, ns, old_wall_side

  real(rp) l_off, del_l, l0, l1, l_try

  ! set pointers
  positive_x_wall => walls%positive_x_wall
  negative_x_wall => walls%negative_x_wall

  ! power calculation is only for bends, quads and wigglers.

  ele => lat%ele(ie)

  if (ele%key /= sbend$ .and. ele%key /= quadrupole$ .and. &
                     ele%key /= wiggler$ .and. ele%key /= sol_quad$) return

  ! check if ele is on
  if (.not. ele%is_on) return

  ! Partition the quad/bend/wiggler into slices and track the synch
  ! radiation comming from each slice to see where it hits the wall.
  ! this is called a "ray".

  ! seg_power_calc takes a set of rays (called a "fan") from a set of consecutive slices.
  ! seg_power_calc demands that all the rays of a fan are all hitting the same wall.
  ! Thus, if one ray hits one wall and the ray from the next ray hits the other wall, 
  ! then we need to divide the slace to find the "transition" ray.

  i_ray = 0
  n_slice = gen%n_slice 

  if (ele%key == wiggler$) then

    ! If periodic wiggler, do n_slices per pole
    if (ele%sub_key == periodic_type$) then
      if (ele%value(b_max$) == 0) return
      if (ele%value(n_pole$) == 0) then
        print *, 'WARNING IN ELE_SYNRAD_POWER: "N_POLE" FOR WIGGLER = 0.'
        print *, '      CALCULATED RADIATION FROM THIS WIGGLER WILL BE 0!'
      endif
      n_slice = max(1, nint(n_slice * ele%value(n_pole$)))
    else
      ! Rather arbitrary choice of 10 * n_slice for non-periodic wigglers
      n_slice = max(10*n_slice, ele%num_steps)
    endif      
  endif

  if (.not. allocated (fan)) allocate(fan(n_slice+1))
  if (n_slice+1 > size(fan)) then
    deallocate(fan)
    allocate(fan(n_slice+1))
  endif

  ! each change in position is the element length / n_slice

  old_wall_side = 0  ! Are we hitting +x or -x side? (ignore ends here)
  del_l = ele%value(l$) / n_slice

  do ns = 0, n_slice

    l_off = ns * del_l
    i_ray = i_ray + 1

    ! start a ray from each slice location
    call init_ray (fan(i_ray), lat, ie, l_off, orb, direction)
    ! track the ray until it hits something
    call track_ray_to_wall (fan(i_ray), lat, walls)

    ! check if this ray hit a different wall than the previous rays

    if (fan(i_ray)%wall_side /= 0 .and. old_wall_side /= 0 .and. &
        fan(i_ray)%wall_side /= old_wall_side) then

      ray_temp = fan(i_ray)
      fan(i_ray) = fan(i_ray-1)
      l0 = l_off - del_l
      l1 = l_off

      ! binary search for transition point
      do
        l_try = (l0 + l1) / 2
        call init_ray (ray, lat, ie, l_try, orb, direction)
        call track_ray_to_wall (ray, lat, walls)
        if (ray%wall_side /= 0 .and. ray%wall_side == old_wall_side) then
          fan(i_ray) = ray
          l0 = l_try
        else
          ray_temp = ray
          l1 = l_try
        endif

        ! if the transition point has been found to within 0.5 mm,
        ! calc the first wall's power then reset the ray array 
        ! for the second wall
        if ((l1 - l0) .le. 5e-4) then
          if (abs(fan(i_ray)%start%vec(5) - fan(i_ray-1)%start%vec(5)) < 5e-4) i_ray = i_ray - 1
          call seg_power_calc (fan, i_ray, walls, old_wall_side, lat, gen, power(ie))
          fan(1) = ray_temp  ! Reset fan. First ray in fan is the transition ray.
          i_ray = 1
          old_wall_side = 0
          exit
        endif
      enddo
    endif

    if (fan(i_ray)%wall_side /= 0) old_wall_side = fan(i_ray)%wall_side

  enddo

  ! Final call to seg_power_calc to compute the power hitting the wall 
  ! using the tracking results and interpolating in between.

  call seg_power_calc (fan, i_ray, walls, old_wall_side, lat, gen, power(ie))

end subroutine ele_synrad_power
