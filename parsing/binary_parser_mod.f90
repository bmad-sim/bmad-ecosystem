module binary_parser_mod

use bmad_interface

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_cartesian_map (file_name, cart_map, err_flag)
!
! Routine to write a binary cartesian_map structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   cart_map      -- cartesian_map_struct: Cartesian map.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_cartesian_map (file_name, cart_map, err_flag)

implicit none

type (cartesian_map_struct), target :: cart_map

integer i, j, iu, ios
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'write_binary_cartesian_map'
character(200) f_name

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) cart_map%field_scale, cart_map%r0, cart_map%master_parameter, &
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
call out_io (s_error$, r_name, 'WRITE ERROR!')
close (iu)
return

end subroutine write_binary_cartesian_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_cartesian_map (file_name, cart_map, err_flag)
!
! Routine to read a binary cartesian_map structure.
!
! Input:
!   file_name     -- character(*): File to create.
!
! Ouput:
!   cart_map      -- cartesian_map_struct, cartesian map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_cartesian_map (file_name, cart_map, err_flag)

implicit none

type (cartesian_map_struct), target :: cart_map

integer i, j, iu, nt, iver
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'read_binary_cartesian_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) cart_map%field_scale, cart_map%r0, cart_map%master_parameter, &
                      cart_map%ele_anchor_pt, cart_map%field_type

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
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

end subroutine read_binary_cartesian_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_cylindrical_map (file_name, cl_map, err_flag)
!
! Routine to write a binary cylindrical_map structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   cl_map        -- cylindrical_map_struct: Cylindrical map.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_cylindrical_map (file_name, cl_map, err_flag)

implicit none

type (cylindrical_map_struct), target :: cl_map

integer i, j, iu, ios
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'write_binary_cylindrical_map'
character(200) f_name

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) cl_map%field_scale, cl_map%master_parameter, cl_map%harmonic, &
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
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

end subroutine write_binary_cylindrical_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_cylindrical_map (file_name, cl_map, err_flag)
!
! Routine to read a binary cylindrical_map structure.
!
! Input:
!   file_name     -- character(*): File to create.
!
! Ouput:
!   cl_map        -- cylindrical_map_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_cylindrical_map (file_name, cl_map, err_flag)

implicit none

type (cylindrical_map_struct), target :: cl_map

integer i, j, iu, nt, iver
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'read_binary_cylindrical_map'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) cl_map%field_scale, cl_map%master_parameter, cl_map%harmonic, &
                cl_map%phi0_fieldmap, cl_map%theta0_azimuth, cl_map%ele_anchor_pt, cl_map%m, cl_map%dz, cl_map%r0

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
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

end subroutine read_binary_cylindrical_map

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_grid_field (file_name, g_field, err_flag)
!
! Routine to write a binary grid_field structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   g_field       -- grid_field_struct: Cylindrical map.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_grid_field (file_name, g_field, err_flag)

implicit none

type (grid_field_struct), target :: g_field

integer i, j, k, n, iu, ios
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'write_binary_grid_field'
character(200) f_name

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) g_field%field_scale, g_field%master_parameter, &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type

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
call out_io (s_error$, r_name, 'WRITE ERROR!')
close (iu)
return

end subroutine write_binary_grid_field

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_grid_field (file_name, g_field, err_flag)
!
! Routine to read a binary grid_field structure.
!
! Input:
!   file_name     -- character(*): File to create.
!
! Ouput:
!   g_field       -- grid_field_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_grid_field (file_name, g_field, err_flag)

implicit none

type (grid_field_struct), target :: g_field

integer i, j, k, iu, n0(3), n1(3), iver
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'read_binary_grid_field'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) g_field%field_scale, g_field%master_parameter, &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type

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
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

end subroutine read_binary_grid_field

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine write_binary_taylor_field (file_name, t_field, err_flag)
!
! Routine to write a binary taylor_field structure.
! Note: The file name should have a ".bin" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   t_field       -- taylor_field_struct: Cylindrical map.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine write_binary_taylor_field (file_name, t_field, err_flag)

implicit none

type (taylor_field_struct), target :: t_field

integer i, j, k, n, iu, ios
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'write_binary_taylor_field'
character(200) f_name

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'WRITE', iu, r_name, 1)) return

!

write (iu, err = 9000) t_field%field_scale, t_field%master_parameter, t_field%curved_ref_frame, &
          t_field%ele_anchor_pt, t_field%field_type, t_field%dz, t_field%r0, t_field%canonical_tracking

f_name = file_name
write (iu, err = 9000) lbound(t_field%ptr%plane, 1), ubound(t_field%ptr%plane, 1), f_name

do j = lbound(t_field%ptr%plane, 1), ubound(t_field%ptr%plane, 1)
  do k = 1, 3
    write (iu, err = 9000) size(t_field%ptr%plane(j)%field(k)%term)
    do n = 1, size(t_field%ptr%plane(j)%field(k)%term)
      write (iu, err = 9000) t_field%ptr%plane(j)%field(k)%term(n)
    enddo
  enddo
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'WRITE ERROR!')
close (iu)
return

end subroutine write_binary_taylor_field

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine read_binary_taylor_field (file_name, t_field, err_flag)
!
! Routine to read a binary taylor_field structure.
!
! Input:
!   file_name     -- character(*): File to create.
!
! Ouput:
!   t_field       -- taylor_field_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine read_binary_taylor_field (file_name, t_field, err_flag)

implicit none

type (taylor_field_struct), target :: t_field

integer i, j, k, iu, n0, n1, n, nn, iver
logical err_flag

character(*) file_name
character(*), parameter :: r_name = 'read_binary_taylor_field'

!

err_flag = .true.
if (.not. open_binary_file(file_name, 'READ', iu, r_name, iver)) return

!

read (iu, err = 9000) t_field%field_scale, t_field%master_parameter, t_field%curved_ref_frame, &
          t_field%ele_anchor_pt, t_field%field_type, t_field%dz, t_field%r0, t_field%canonical_tracking

allocate (t_field%ptr)
read (iu, err = 9000) n0, n1, t_field%ptr%file
allocate (t_field%ptr%plane(n0:n1))

do j = n0, n1
  do k = 1, 3
    read (iu, err = 9000) nn
    allocate (t_field%ptr%plane(j)%field(k)%term(nn))
    do n = 1, nn
       read (iu, err = 9000) t_field%ptr%plane(j)%field(k)%term(n)
    enddo
  enddo
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

end subroutine read_binary_taylor_field

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
    call out_io (s_error$, r_name, 'VERSION NUMBER NOT RECOGNIZED!')
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
call out_io (s_error$, r_name, 'READ ERROR!')
close (iu)
return

9100 continue
call out_io (s_error$, r_name, 'WRITE ERROR!')
close (iu)
return

end function open_binary_file

end module
