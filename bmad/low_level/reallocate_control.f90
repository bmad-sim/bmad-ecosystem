!+
! Subroutine reallocate_control(lat, n) 
!
! Routine to reallocate the lat%control(:) and lat%ic(:) arrays.
! The old data in the arrays will be saved.
! 
! Input:
!   lat  -- Lat_struct: Lattice.
!   n    -- Integer: Array size for lat%control(:) and lat%ic(:).
!
! Output:
!   lat  -- Lat_struct: Lattice.
!     %control(:) -- Control Array with size at least n.
!     %ic(:)      -- Control Array.
!-

subroutine reallocate_control (lat, n)

use bmad_struct

implicit none

type (lat_struct) lat
type (control_struct), allocatable :: control(:)
integer, intent(in) :: n
integer i, n_old

!

if (.not. allocated(lat%control)) then
  allocate (lat%control(n), lat%ic(n))
  lat%ic = 0
  return
endif

n_old = size(lat%control)
if (n_old >= n) return

call move_alloc (lat%control, control)

allocate (lat%control(n))
do i = 1, n_old
  call move_alloc(control(i)%stack, lat%control(i)%stack)
  call move_alloc(control(i)%y_knot, lat%control(i)%y_knot)
  lat%control(i)%lord      = control(i)%lord
  lat%control(i)%slave     = control(i)%slave
  lat%control(i)%ix_attrib = control(i)%ix_attrib
  lat%control(i)%attribute = control(i)%attribute
enddo

call re_allocate(lat%ic, max(n, size(lat%ic) + n - n_old))
lat%ic(n_old+1:) = 0

end subroutine reallocate_control

