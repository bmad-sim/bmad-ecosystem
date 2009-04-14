!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine init_ray (ray, lat, ix_ele, l_offset, orb, direction)
!
! subroutine to start a new synch radiation ray
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   lat   -- lat_struct with twiss propagated and mat6s made
!   ix_ele -- integer: index of lat element to start ray from
!   l_offset -- real(rp): offset along the length of the element 
!                         to use as a starting point for ray
!   orb(0:*) -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 for in direction of s
!                         -1 for against s
!
! Output:
!   ray    -- ray_struct: synch radiation ray with starting
!                         parameters set
!-

subroutine init_ray (ray, lat, ix_ele, l_offset, orb, direction)

  use synrad_struct
  use synrad_interface, except => init_ray
  use boris_mod

  implicit none

  type (lat_struct), target :: lat
  type (coord_struct) :: orb(0:*), orb0, orb1
  type (ele_struct), pointer :: ele0, ele
  type (ele_struct) :: runt_ele
  type (ray_struct) :: ray
  type (em_field_struct) :: field

  real(rp) l_offset, k_wig, g_max, l_small

  integer direction, ix_ele

! initialize the ray

  ele0 => lat%ele(ix_ele-1)
  ele => lat%ele(ix_ele)

  runt_ele = ele                             ! make a runt element

  ! set the ray's initial twiss values

  call twiss_and_track_partial (ele0, ele, lat%param, l_offset, runt_ele, & 
                                                      orb(ix_ele-1), orb0)
  ray%x_twiss = runt_ele%a
  ray%y_twiss = runt_ele%b

  ! set the ray's g_bend value (inverse bending radius at src pt) 

  if (ele%key == sbend$) then  

    ! sbends are easy
    ray%g_bend = abs(1 / ele%value(rho$))

  elseif (ele%key == quadrupole$ .or. ele%key == sol_quad$) then

    ! for quads or sol_quads, get the bending radius
    ! from the change in x' and y' over a small 
    ! distance in the element
    l_small = 1e-2      ! something small
    runt_ele%value(l$) = l_small
    call make_mat6 (runt_ele, lat%param, orb0, orb0)
    call track1 (orb0, runt_ele, lat%param, orb1)
    orb1%vec = orb1%vec - orb0%vec
    ray%g_bend = sqrt(orb1%vec(2)**2 + orb1%vec(4)**2) / l_small

  elseif (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) then

    ! for periodic wigglers, get the max g_bend from 
    ! the max B field of the wiggler, then scale it 
    ! by the cos of the position along the poles
    ! Note: assumes particles are relativistic!!
    k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))

    ! changed to use the E_TOT at the ele 2006.11.24 mjf
    g_max = c_light * ele%value(b_max$) / (ele%value(E_TOT$))
    ray%g_bend = abs(g_max * cos (k_wig * l_offset))
    orb0%vec(2) = orb0%vec(2) + (g_max / k_wig) * sin (k_wig * l_offset)

  elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then

    ! for mapped wigglers, find the B field at the source point
    ! and extract the g_bend
    ! Note: assumes particles are relativistic!!
    call em_field_calc (runt_ele, lat%param, l_offset, orb0, .false., field)

    ! changed to use the E_TOT at the ele 2006.11.24 mjf
    ray%g_bend = sqrt(sum(field%b(1:2)**2)) * c_light / ele%value(E_TOT$)

  else

    print *, 'ERROR: UNKNOWN ELEMENT HERE ', ele%name

  endif

  ray%ix_source = ix_ele
  ray%start = orb0
  ray%start%vec(5) = ele0%s + l_offset                    ! s position
  ray%now = ray%start
  ray%ix_ele = ix_ele
  ray%crossed_end = .false.
  ray%direction = direction
  ray%track_len = 0
  ray%alley_status = no_local_alley$

end subroutine
