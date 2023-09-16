!+
! Subroutine create_field_overlap (lat, lord_name, slave_name, err_flag)
!
! Subroutine to add the bookkeeping information to a lattice for an element's field
! overlapping another element.
!
! The number of lord/slave overlaps is equal to the number of lord elements that match lord_name. 
! Each lord is matched with the nearest slave element.
! "Nearest" is decided by distance in global coordinates.
!
! Input:
!   lat         -- lat_struct: Lattice
!   lord_name   -- character(*): Name of the element with a field extending beyound it's bounds.
!   slave_name  -- character(*): Name of the element the lord's field overlaps.
!
! Output:
!   lat         -- lat_struct: Lattice
!   err_flag    -- logical: Set true if there is a problem (like no elements found).
!-

subroutine create_field_overlap (lat, lord_name, slave_name, err_flag)

use bmad_interface, except_dummy => create_field_overlap

implicit none

type (lat_struct) :: lat
type (ele_pointer_struct), allocatable :: lords(:), slaves(:)
type (ele_struct), pointer :: lord, slave
type (floor_position_struct) lord_floor, slave_floor

real(rp) dr(3), min_dist2, dist2

integer i, j, n_lord, n_slave, j_min, n, n2

logical err_flag, err

character(*) lord_name, slave_name

!

err_flag = .true.

call lat_ele_locator (lord_name, lat, lords, n_lord, err)
if (err .or. n_lord == 0) return

call lat_ele_locator (slave_name, lat, slaves, n_slave, err)
if (err .or. n_slave == 0) return

do i = 1, n_lord
  lord => lords(i)%ele
  call ele_geometry(lord%floor, lord, lord_floor, -0.5_rp)

  min_dist2 = 1e30
  do j = 1, n_slave
    slave => slaves(j)%ele
    call ele_geometry(slave%floor, slave, slave_floor, -0.5_rp)
    dr = lord_floor%r - slave_floor%r
    dist2 = dot_product(dr, dr)
    if (dist2 < min_dist2) then
      j_min = j
      min_dist2 = dist2
    endif
  enddo
  slave => slaves(j_min)%ele

  call add_lattice_control_structs (lord, n_add_slave_field = 1)
  n = lord%ix1_slave + lord%n_slave + lord%n_slave_field - 1
  lat%control(n)%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
  lat%control(n)%ix_attrib = field_overlaps$
  lat%control(n)%attribute = 'FIELD_OVERLAPS'

  call add_lattice_control_structs (slave, n_add_lord_field = 1)
  n2 = slave%ic1_lord + slave%n_lord + slave%n_lord_field - 1
  lat%ic(n2) = n
enddo

err_flag = .false.

end subroutine
