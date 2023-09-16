!+
! subroutine ele_synrad_power (branch, ie, orb, direction, power, walls, gen)
!
! subroutine to calculate the synch radiation power from one element of the lattice
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   branch  -- branch_struct: with twiss propagated and mat6s made
!   ie      -- integer: index of the current element in branch
!   orb(0:) -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 for in direction of s
!                         -1 for against s
!   walls  -- walls_struct: both walls with outline ready
!   gen    -- synrad_param_struct: Contains lattice name,
!                     vert emittance, and beam current
!
! Output:
!   power(:) -- ele_power_struct: power radiated from a lattice ele
!   wall   -- walls_struct: both walls with power information
!-

subroutine ele_synrad_power (branch, ie, orb, direction, power, walls, gen)

  use synrad_struct
  use synrad_interface, except => ele_synrad_power

  implicit none

  type (branch_struct), target :: branch
  type (coord_struct) orb(0:)
  type (ele_struct), pointer :: ele, field_ele
  type (walls_struct), target :: walls
  type (wall_struct), pointer :: negative_x_wall, positive_x_wall
  type (ray_struct), allocatable, target :: fan(:)
  type (ray_struct) this_ray, ray_temp
  type (ray_struct), pointer :: ray
  type (synrad_param_struct) gen
  type (ele_power_struct) power(:)

  integer direction, ie, i_ray, n_slice, ns, old_wall_side
  integer, save :: n_warn = 0

  real(rp) l_off, del_l, l0, l1, l_try

  character(*), parameter :: r_name = 'ele_synrad_power'

  ! set pointers
  positive_x_wall => walls%positive_x_wall
  negative_x_wall => walls%negative_x_wall

  ! power calculation is only for bends, quads and wigglers.

  ele => branch%ele(ie)
  select case (ele%key)
  case (elseparator$, sbend$, quadrupole$, sad_mult$, wiggler$, undulator$, sol_quad$)
  case default
    return
  end select

  ! check if ele is on
  if (.not. ele%is_on) return

  ! Zero length elements cannot be handled

  if (ele%value(l$) == 0) then
    n_warn = n_warn + 1
    if (n_warn <= 5) then
      call out_io (s_warn$, r_name, 'Element has zero length: ' // ele%name, &
                                    'No radiation will be generated from this element.')
    endif
    if (n_warn == 5) then
      call out_io (s_info$, r_name, 'Enough! Zero length warnings will now be suppressed...')
    endif
    return
  endif

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
    field_ele => pointer_to_field_ele(ele, 1)
    if (field_ele%field_calc == planar_model$ .or. field_ele%field_calc == helical_model$) then
      if (ele%value(b_max$) == 0) return
      if (ele%value(n_period$) == 0) then
        call out_io (s_warn$, r_name, '"N_Period" for Wiggler = 0.', &
                                      'Calculated radiation from this wiggler will be 0!')
      endif
      n_slice = max(1, nint(2 * n_slice * ele%value(n_period$)))
    else
      ! Rather arbitrary choice of 10 * n_slice for non-periodic wigglers
      n_slice = max(10*n_slice, nint(ele%value(num_steps$)))
    endif      
  endif

  allocate(fan(n_slice+1))

  ! each change in position is the element length / n_slice

  old_wall_side = 0  ! Are we hitting +x or -x side? (ignore ends here)
  del_l = ele%value(l$) / n_slice

  do ns = 0, n_slice

    l_off = ns * del_l
    i_ray = i_ray + 1
    ray => fan(i_ray)

    ! start a ray from each slice location
    call init_ray (ray, branch, ie, l_off, orb, direction)
    if (ray_is_outside_wall(ray, walls)) then
      print *, 'PHOTON GENERATED OUTSIDE WALL AT (X, S):', ray%now%vec(1), ray%now%s
      call err_exit
    endif

    ! track the ray until it hits something
    call track_ray_to_wall (ray, walls)

    ! check if this ray hit a different wall than the previous rays

    if (ray%wall_side /= 0 .and. old_wall_side /= 0 .and. &
        ray%wall_side /= old_wall_side) then

      ray_temp = ray
      ray = fan(i_ray-1)
      l0 = l_off - del_l
      l1 = l_off

      ! binary search for transition point
      do
        l_try = (l0 + l1) / 2
        call init_ray (this_ray, branch, ie, l_try, orb, direction)
        if (ray_is_outside_wall(this_ray, walls)) then
          print *, 'PHOTON GENERATED OUTSIDE WALL AT (X, S):', this_ray%now%vec(1), this_ray%now%s
          call err_exit
        endif
        call track_ray_to_wall (this_ray, walls)
        if (this_ray%wall_side /= 0 .and. this_ray%wall_side == old_wall_side) then
          ray = this_ray
          l0 = l_try
        else
          ray_temp = this_ray
          l1 = l_try
        endif

        ! If the transition point has been found to within 0.5 mm,
        ! calc the first wall's power then reset the ray array 
        ! for the second wall
        if ((l1 - l0) .le. 5e-4) then
          if (abs(fan(i_ray)%start%s - fan(i_ray-1)%start%s) < 5e-4) i_ray = i_ray - 1
          call seg_power_calc (fan, i_ray, walls, old_wall_side, branch, gen, power(ie))
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

  call seg_power_calc (fan, i_ray, walls, old_wall_side, branch, gen, power(ie))

end subroutine ele_synrad_power
