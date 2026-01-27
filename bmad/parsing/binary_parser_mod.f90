module binary_parser_mod

use attribute_mod    ! Uses attribute_name and attribute_index

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_cartesian_map (file_name, ele, cart_map, err_flag)
!
! Routine to write a binary cartesian_map structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   cart_map      -- cartesian_map_struct: Cartesian map.
!
! Output:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_cartesian_map (file_name, ele, cart_map, err_flag)

implicit none

type (cartesian_map_struct), target :: cart_map
type (ele_struct) ele
integer i, j, iu, ios
logical err_flag

character(*) file_name
character(400) f_name
character(*), parameter :: r_name = 'write_binary_cartesian_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) cart_map%field_scale, cart_map%r0, attribute_name(ele, cart_map%master_parameter), &
           cart_map%ele_anchor_pt, cart_map%field_type

f_name = file_name
write (iu, err = 9000) size(cart_map%ptr%term), f_name

do j = 1, size(cart_map%ptr%term)
  write (iu, err = 9000) cart_map%ptr%term(j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR WRITING FIELDMAP STRUCTURE. FILE: ' // file_name)
close (iu)
return

end subroutine write_binary_cartesian_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_cartesian_map (file_name, ele, cart_map, err_flag)
!
! Routine to read a binary cartesian_map structure.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!
! Output:
!   cart_map      -- cartesian_map_struct, cartesian map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_cartesian_map (file_name, ele, cart_map, err_flag)

implicit none

type (cartesian_map_struct), target :: cart_map
type (ele_struct) ele

integer i, j, iu, nt, iver
logical err_flag

character(*) file_name
character(40) master_name
character(*), parameter :: r_name = 'read_binary_cartesian_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) cart_map%field_scale, cart_map%r0, master_name, &
                      cart_map%ele_anchor_pt, cart_map%field_type
cart_map%master_parameter = attribute_index(ele, master_name)

allocate (cart_map%ptr)
read (iu, err = 9000) nt, cart_map%ptr%file
allocate (cart_map%ptr%term(nt))

do j = 1, nt
  read (iu, err = 9000) cart_map%ptr%term(j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR READING BINARY FIELDMAP FILE. FILE: ' // file_name)
close (iu)
return

end subroutine read_binary_cartesian_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_cylindrical_map (file_name, ele, cl_map, err_flag)
!
! Routine to write a binary cylindrical_map structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   cl_map        -- cylindrical_map_struct: Cylindrical map.
!
! Output:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_cylindrical_map (file_name, ele, cl_map, err_flag)

implicit none

type (cylindrical_map_struct), target :: cl_map
type (ele_struct) ele

integer i, j, iu, ios
logical err_flag

character(*) file_name
character(400) f_name
character(*), parameter :: r_name = 'write_binary_cylindrical_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) cl_map%field_scale, attribute_name(ele, cl_map%master_parameter), cl_map%harmonic, &
                cl_map%phi0_fieldmap, cl_map%theta0_azimuth, cl_map%ele_anchor_pt, cl_map%m, cl_map%dz, cl_map%r0


f_name = file_name
write (iu, err = 9000) size(cl_map%ptr%term), f_name

do j = 1, size(cl_map%ptr%term)
  write (iu, err = 9000) cl_map%ptr%term(j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR READING BINARY FIELDMAP FILE. FILE: ' // file_name)
close (iu)
return

end subroutine write_binary_cylindrical_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_cylindrical_map (file_name, ele, cl_map, err_flag)
!
! Routine to read a binary cylindrical_map structure.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!
! Output:
!   cl_map        -- cylindrical_map_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_cylindrical_map (file_name, ele, cl_map, err_flag)

implicit none

type (cylindrical_map_struct), target :: cl_map
type (ele_struct) ele

integer i, j, iu, nt, iver
logical err_flag

character(*) file_name
character(40) master_name
character(*), parameter :: r_name = 'read_binary_cylindrical_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) cl_map%field_scale, master_name, cl_map%harmonic, &
                cl_map%phi0_fieldmap, cl_map%theta0_azimuth, cl_map%ele_anchor_pt, cl_map%m, cl_map%dz, cl_map%r0
cl_map%master_parameter = attribute_index(ele, master_name)

allocate (cl_map%ptr)
read (iu, err = 9000) nt, cl_map%ptr%file
allocate (cl_map%ptr%term(nt))

do j = 1, nt
  read (iu, err = 9000) cl_map%ptr%term(j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR READING BINARY FIELDMAP FILE. FILE: ' // file_name)
close (iu)
return

end subroutine read_binary_cylindrical_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to write a binary grid_field structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   g_field       -- grid_field_struct: Cylindrical map.
!
! Output:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_grid_field (file_name, ele, g_field, err_flag)

implicit none

type (grid_field_struct), target :: g_field
type (ele_struct) ele

integer i, j, k, n, iu, ios
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'write_binary_grid_field'
character(400) f_name

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) g_field%field_scale, attribute_name(ele, g_field%master_parameter), &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type, g_field%interpolation_order

f_name = file_name
write (iu, err = 9000) lbound(g_field%ptr%pt), ubound(g_field%ptr%pt), f_name

do j = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)
  write (iu, err = 9000) g_field%ptr%pt(:, :, j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR WRITING FIELDMAP STRUCTURE. FILE: ' // file_name)
close (iu)
return

end subroutine write_binary_grid_field

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to read a binary grid_field structure.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!
! Output:
!   g_field       -- grid_field_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_grid_field (file_name, ele, g_field, err_flag)

implicit none

type (grid_field_struct), target :: g_field
type (ele_struct) ele

integer i, j, k, iu, n0(3), n1(3), iver
logical err_flag

character(*) file_name
character(40) master_name
character(*), parameter :: r_name = 'read_binary_grid_field'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) g_field%field_scale, master_name, &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type, g_field%interpolation_order
g_field%master_parameter = attribute_index(ele, master_name)

allocate (g_field%ptr)
read (iu, err = 9000) n0, n1, g_field%ptr%file
allocate (g_field%ptr%pt(n0(1):n1(1), n0(2):n1(2), n0(3):n1(3)))

do j = n0(3), n1(3)
  read (iu, err = 9000) g_field%ptr%pt(:, :, j)
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR READING BINARY FIELDMAP FILE. FILE: ' // file_name)
close (iu)
return

end subroutine read_binary_grid_field

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function open_binary_file (file_name, action, iu, r_name, iver) result (is_ok)
!
! Routine to open a binary file for reading or writing.
!
! Input:
!   file_name     -- character(*): File to create.
!   action        -- character(*): 'READ' or 'WRITE'
!   r_name        -- character(*): Calling routine name for error messages.
!
! Output:
!   iu            -- integer: Unit number of opened file.
!   iver          -- integer: Version number if action = 'READ'
!   is_ok         -- logical: Open OK?
!-


function open_binary_file (file_name, action, iu, r_name, iver) result (is_ok)

implicit none

integer iu, ios, iostat, iver
character(*) file_name, action, r_name
logical is_ok

! Open file

is_ok = .false.
iu = lunget()

select case (action)
case ('READ')
  open (iu, file = file_name, form = 'unformatted', action = 'READ', status = 'old', iostat = ios)

case ('WRITE')
  open (iu, file = file_name, form = 'unformatted', iostat = ios)

case default
  call err_exit
end select

if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name)
  return
endif

! Read or Write the version number

select case (action)
case ('READ')
  read (iu, err = 9000) iver     ! Version number
  if (iver /= 1) then
    call out_io (s_error$, r_name, 'VERSION NUMBER NOT RECOGNIZED. FILE: ' // file_name)
    close (iu)
    return
  endif

case ('WRITE')
  write (iu, err = 9100) iver     ! Version number
end select

!

is_ok = .true.
return

!
 
9000 continue
call out_io (s_error$, r_name, 'ERROR READING VERSION NUMBER FROM BINARY FILE. FILE: ' // file_name)
close (iu)
return

9100 continue
call out_io (s_error$, r_name, 'ERROR WRITING VERSION NUMBER TO FIELDMAP BINARY FILE. FILE: ' // file_name)
close (iu)
return

end function open_binary_file

end module
