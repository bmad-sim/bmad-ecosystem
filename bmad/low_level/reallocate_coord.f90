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

use bmad_routine_interface, dummy => reallocate_coord_n
implicit none

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

use bmad_routine_interface, dummy => reallocate_coord_lat
implicit none

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
!   coord(:) -- Coord_array_struct, allocatable: Allocatable array.
!   lat      -- lat_struct: 
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord_array (coord_array, lat)

use bmad_routine_interface
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

