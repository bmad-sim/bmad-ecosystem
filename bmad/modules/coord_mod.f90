module coord_mod

use equal_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine reallocate_coord (...)
!
! Routine to allocate or reallocate at allocatable coord_struct array.
! reallocate_coord is an overloaded name for:
!   reallocate_coord_n (coord, n_coord)
!   reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable coord_struct array to at least:
!     coord(0:n_coord)                            if n_coord arg is used.
!     coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.
!
! The old coordinates are saved
! If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
! In any case, coord(n)%vec for n > 0 is set to zero.
!
! Input:
!   coord(:)  -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!   lat       -- lat_struct: Lattice 
!   ix_branch -- Integer, optional: Branch to use. Default is 0 (main branch).
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

interface reallocate_coord
  module procedure reallocate_coord_n
  module procedure reallocate_coord_lat
end interface

private reallocate_coord_n, reallocate_coord_lat

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord (...)
!
! Routine to initialize a coord_struct. 
!
! This routine is an overloaded name for:
!   Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)
!   Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)
!   Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)
!
! Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice), 
! or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6) 
! corresponding to the beginning_ele's value of ele%value(p0c_start$).
!
! Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead 
! of the standard which is to set orb%t from orb%vec(5).
!
! For photons:
!   orb%vec(5) is set depending upon where the photon is relative to the element.
!   If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle 
!       except if direction is set.
!
! Input:
!   orb_in       -- coord_struct: Input orbit.
!   vec(6)       -- real(rp), optional: Coordinate vector. If not present then taken to be zero.
!   ele          -- ele_struct, optional: Particle is initialized to start at element_end of this ele.
!   element_end  -- integer, optional: upstream_end$, downstream_end$, inside$, or start_end$.
!                     Must be present if ele argument is present.
!                     start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1.
!                     Default is upstream_end$. Note: If ele is the beginning element (index zero), the
!                     setting of element_end will not matter.
!   particle     -- integer, optional: Particle type (electron$, etc.). 
!                     If particle = not_set$ and orb_in is present, use orb_in%species instead.
!   dirction     -- integer, optional: +1 -> moving downstream +s direciton, -1 -> moving upstream.
!                     0 -> Ignore. Default is to not change orb%direction except for photons which get set
!                     according to orb%vec(6).
!   E_photon     -- real(rp), optional: Photon energy if particle is a photon. Ignored otherwise.
!   t_offset     -- real(rp), optional: Offset of the reference time. This is non-zero when
!                     there are multiple bunches and the reference time for a particular particle
!                     is pegged to the time of the center of the bunch.
!   shift_vec6   -- logical, optional: If present and False, prevent the shift of orb%vec(6).
!   spin(3)      -- real(rp), optional: Particle spin. Taken to be zero if not present.
!   s_pos        -- real(rp), optional: Particle s-position. Only relavent if element_end = inside$.
!
! Output:
!   orb -- Coord_struct: Initialized coordinate.
!                 Note: For photons, orb%vec(6) is computed as sqrt(1 - vec(2)^2 - vec(4)^2) if needed.
!-

interface init_coord
  module procedure init_coord1
  module procedure init_coord2
  module procedure init_coord3
end interface

private init_coord1, init_coord2

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_n (coord, n_coord)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_n (coord, n_coord)

type (coord_struct), allocatable :: coord(:)
type (coord_struct), allocatable :: old(:)

integer, intent(in) :: n_coord
integer i, n_old

character(*), parameter :: r_name = 'reallocate_coord_n'

!

if (allocated (coord)) then

  if (lbound(coord, 1) /= 0) then
    call out_io (s_fatal$, r_name, 'ORBIT ARRAY LOWER BOUND NOT EQUAL TO ZERO!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

  n_old = ubound(coord, 1)
  if (n_old >= n_coord) return
  allocate(old(0:n_old))

  do i = 0, n_old
    old(i) = coord(i)
  enddo

  deallocate (coord)
  allocate (coord(0:n_coord))

  do i = 0, n_old
    coord(i) = old(i)
  enddo

  deallocate(old)

else
  allocate (coord(0:n_coord))
endif

end subroutine reallocate_coord_n

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_lat (coord, lat, ix_branch)

type (coord_struct), allocatable :: coord(:)
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer, optional :: ix_branch

!

branch => lat%branch(integer_option(0, ix_branch))

if (allocated(coord)) then
  call reallocate_coord_n (coord, branch%n_ele_max)
else
  allocate (coord(0:branch%n_ele_max))
endif

end subroutine reallocate_coord_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_array (coord_array, lat)
!
! Subroutine to allocate an allocatable coord_array_struct array to
! the proper size for a lattice.
!
! Note: Any old coordinates are not saved except for coord_array(:)%orbit(0).
! If, at input, coord_array is not allocated, coord_array(:)%orbit(0)%vec is set to zero.
! In any case, all other %vec components are set to zero.
!
! Input:
!   coord(:) -- Coord_struct, allocatable: Allocatable array.
!   lat      -- lat_struct: 
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord_array (coord_array, lat)

implicit none

type (coord_array_struct), allocatable :: coord_array(:)
type (lat_struct) lat
type (coord_struct), allocatable :: start(:)

integer i, j, nb

!

if (.not. allocated(lat%branch)) return
nb = ubound(lat%branch, 1)

if (allocated (coord_array)) then
  if (size(coord_array) /= nb + 1) then
    call reallocate_coord(start, nb)
    do i = 0, nb
      start(i) = coord_array(i)%orbit(0)
    enddo
    deallocate (coord_array)
    allocate (coord_array(0:nb))
    do i = 0, nb
      call reallocate_coord (coord_array(i)%orbit, lat%branch(i)%n_ele_max)
      coord_array(i)%orbit(0) = start(i)
    enddo
  endif
else
  allocate (coord_array(0:nb))
  do i = 0, nb
    call reallocate_coord (coord_array(i)%orbit, lat%branch(i)%n_ele_max)
  enddo
endif

end subroutine reallocate_coord_array

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)

implicit none

type (coord_struct) orb, orb_temp
type (ele_struct), optional :: ele
real(rp) :: vec(6)
real(rp), optional :: t_offset, E_photon, spin(3), s_pos
integer, optional :: element_end, particle, direction
logical, optional :: shift_vec6

!

orb_temp = coord_struct()
orb_temp%vec = vec

call init_coord2 (orb, orb_temp, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)

end subroutine init_coord1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord2 (orb_out, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord2 (orb_out, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos)

implicit none

type (coord_struct) orb_out, orb_in, orb
type (ele_struct), optional, target :: ele
type (branch_struct), pointer :: branch

real(rp), optional :: E_photon, t_offset, spin(3), s_pos
real(rp) p0c, e_tot, ref_time

integer, optional :: element_end, particle, direction
integer species, dir

logical, optional :: shift_vec6

character(*), parameter :: r_name = 'init_coord2'

!

orb = orb_in
if (present(ele)) then
  branch => pointer_to_branch(ele)
else
  branch => null()
endif

species = orb_in%species
if (present(particle)) then
  if (particle /= not_set$) species = particle
endif

if (present(particle)) orb%species = particle

if (orb%species == not_set$ .and. present(ele) .and. associated (branch)) then
  orb%species = default_tracking_species(branch%param)
endif

if (orb%species == not_set$) then
  call out_io (s_warn$, r_name, 'NO PARTICLE SPECIES GIVEN. USING POSITRONS AS DEFAULT!')
  orb%species = positron$
endif

orb%state = alive$

dir = integer_option(0, direction)
if (dir == 0) then
  if (orb%species == photon$ .and. orb%vec(6) > 0) orb%direction = 1
  if (orb%species == photon$ .and. orb%vec(6) < 0) orb%direction = -1
else
  orb%direction = dir
  if (orb%species == photon$) orb%vec(6) = orb%direction * abs(orb%vec(6))
endif

! Set location and species

if (present(element_end)) orb%location = element_end

if (orb%location == start_end$) then
  if (orb%direction == 1) then
    orb%location = upstream_end$
  else
    orb%location = downstream_end$
  endif
endif

! spin

if (present(spin)) orb%spin = spin

! Energy values

if (present(ele)) then
  if (.not. present(element_end)) then
    call out_io (s_fatal$, r_name, 'RULE: "ELEMENT_END" ARGUMENT MUST BE PRESENT IF "ELE" ARGUMENT IS.')
    call err_exit
  endif

  if (orb%location == downstream_end$ .or. ele%key == beginning_ele$) then
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    ref_time = ele%ref_time
    orb%s = ele%s
  elseif (orb%location == upstream_end$) then
    p0c = ele%value(p0c_start$)
    e_tot = ele%value(e_tot_start$)
    ref_time = ele%value(ref_time_start$)
    orb%s = ele%s_start
  else
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    orb%s = real_option(ele%s, s_pos)
    ref_time = (ele%value(ref_time_start$) * (ele%s - orb%s) + ele%ref_time * (orb%s - ele%s_start)) / ele%value(l$)
  endif
endif

! Photon

if (orb%species == photon$) then
  if (present(ele)) then
    if (present(ele)) orb%p0c = p0c
    if (ele%key == photon_init$) then
      call init_a_photon_from_a_photon_init_ele (ele, branch%param, orb)
    endif
  endif

  orb%dt_ref = 0
  orb%beta = 1

  if (present(E_photon)) then
    if (E_photon /= 0) orb%p0c = E_photon
  endif

  ! If original particle is not a photon, photon is launched in same direciton as particle 
  if (orb_in%species /= photon$ .and. orb_in%species /= not_set$) then
    orb%vec(2:4:2) = orb%vec(2:4:2) / (1 + orb%vec(6))
    if (dir == 0) orb%direction = orb_in%direction
  endif
  orb%vec(6) = orb%direction * sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  if (orb%location == downstream_end$) then
    orb%vec(5) = ele%value(l$)
  else
    orb%vec(5) = 0
  endif
endif

! If ele is present...

orb%ix_ele = -1
orb%ix_branch = -1

if (present(ele)) then

  if (ele%slave_status == slice_slave$) then
    orb%ix_ele = ele%lord%ix_ele
    orb%ix_branch = ele%lord%ix_branch
  else
    orb%ix_ele = ele%ix_ele
    orb%ix_branch = ele%ix_branch
  endif

  if (ele%key == beginning_ele$) orb%location = downstream_end$

  if (orb%species /= photon$) then

    orb%p0c = p0c

    ! E_gun shift. Only time p0c_start /= p0c for an init_ele is when there is an e_gun present in the branch.
    if (logic_option(.true., shift_vec6)) then
      if (ele%key == beginning_ele$) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_start$) - ele%value(p0c$)) / ele%value(p0c$)
      elseif ((ele%key == e_gun$ .and. orb%location == upstream_end$) .or. (ele%key == marker$ .and. ele%value(e_tot_ref_init$) /= 0)) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_ref_init$) - ele%value(p0c$)) / ele%value(p0c$)
      endif
    endif

    if (orb%vec(6) == 0) then
      orb%beta = p0c / e_tot
    else
      call convert_pc_to (p0c * (1 + orb%vec(6)), orb%species, beta = orb%beta)
    endif

    ! Do not set %t if %beta = 0 since %t may be a good value.

    if (orb%vec(5) == real_garbage$) then
      orb%vec(5) = orb%beta * c_light * (ref_time - orb%t)

    elseif (orb%beta == 0) then
      if (orb%vec(5) /= 0) then
        call out_io (s_error$, r_name, 'Z-POSITION IS NONZERO WITH BETA = 0.', &
                                       'THIS IS NONSENSE SO SETTING Z TO ZERO.')
        orb%vec(5) = 0
      endif
    else
      orb%t = ref_time - orb%vec(5) / (orb%beta * c_light)
      if (present(t_offset)) orb%t = orb%t + t_offset
      if (orb%ix_turn /= 0 .and. associated(ele%branch)) then
        orb%t = orb%t + orb%ix_turn * (ele%branch%ele(ele%branch%n_ele_track)%ref_time - ele%branch%ele(0)%ref_time)
      endif
    endif
  endif

endif

orb_out = orb

end subroutine init_coord2

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)

implicit none

type (coord_struct) orb
type (ele_struct), optional :: ele
real(rp), optional :: t_offset, E_photon, spin(3)
integer, optional :: element_end, particle, direction
logical, optional :: shift_vec6

!

call init_coord2 (orb, coord_struct(), ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)

end subroutine init_coord3

end module
