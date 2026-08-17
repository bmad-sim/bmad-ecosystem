!+
! Subroutine tao_init_building_wall (wall_file)
!
! Subroutine to initialize the building wall.
!
! Input:
!   wall_file  -- character(*): file name containing the wall definition.
!                   If blank then there is no wall.
!-

subroutine tao_init_building_wall (wall_file)

use tao_interface
use tao_input_struct
use tao_nml_mod
use tao_attrib_resolve_mod

implicit none

type (tao_building_wall_point_struct), target :: point(100)
type (tao_building_wall_point_struct), pointer :: pt(:)
type (tao_nml_group_struct) nml_group
type (tao_nml_ref_struct) ref
type (tao_ptr_struct) ptr

real(rp) x_mid, z_mid, dx, dz, a, a2, theta, x_offset, z_offset
integer i, j, k, iu, ios, n_wall, n_val

character(*) wall_file
character(n_file_max_len) complete_file_name
character(40) name
character(16) constraint
character(200) why
character(:), allocatable :: nml_val(:)
character(*), parameter :: r_name = 'tao_init_building_wall'

logical err, nml_eof

! Open file

if (wall_file == '') then
  allocate (s%building_wall%section(0))
  return
endif

call out_io (s_blank$, r_name, '*Init: Opening Building Wall File: ' // wall_file)
call tao_open_file (wall_file, iu, complete_file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_error$, r_name, 'ERROR OPENING BUILDING WALL FILE...')
  return
endif

! Wall orientation

theta = 0
x_offset = 0
z_offset = 0

call tao_nml_group_read (iu, complete_file_name, 'building_wall_orientation', nml_group, nml_eof, err, why)
if (err) then
  call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_ORINETATION NAMELIST', why)
  return
elseif (.not. nml_eof) then
  do i = 1, nml_group%n_item
    call tao_nml_item_ref (nml_group%item(i), ref, err)
    if (err) cycle
    select case (ref%head)
    case ('theta');     call tao_nml_value_set (nml_group%item(i), theta, err)
    case ('x_offset');  call tao_nml_value_set (nml_group%item(i), x_offset, err)
    case ('z_offset');  call tao_nml_value_set (nml_group%item(i), z_offset, err)
    case default;       call tao_nml_unknown (nml_group%item(i), 'building_wall_orientation', err)
    end select
  enddo
  s%building_wall%orientation%theta = theta
  s%building_wall%orientation%x_offset = x_offset
  s%building_wall%orientation%z_offset = z_offset
endif

call tao_nml_rewind (iu)
! Count the number of walls

n_wall = 0
do
  call tao_nml_group_read (iu, complete_file_name, 'building_wall_section', nml_group, nml_eof, err, why)
  if (err) then
    call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_SECTION NAMELIST', why)
    return
  endif
  if (nml_eof) exit
  n_wall = n_wall + 1
enddo

call out_io (s_blank$, r_name, '  Number of building walls: \I2\ ', n_wall)

allocate (s%building_wall%section(n_wall))
if (n_wall == 0) return ! no walls

! Now transfer the information

call tao_nml_rewind (iu)
do i = 1, n_wall

  name = ''
  constraint = 'none'
  point%radius = 0
  point%x = real_garbage$

  call tao_nml_group_read (iu, complete_file_name, 'building_wall_section', nml_group, nml_eof, err, why)
  if (err) then
    call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_SECTION NAMELIST', why)
    return
  endif

  do j = 1, nml_group%n_item
    call tao_nml_item_ref (nml_group%item(j), ref, err)
    if (err) cycle
    select case (ref%head)
    case ('name');        call tao_nml_value_set (nml_group%item(j), name, err)
    case ('constraint');  call tao_nml_value_set (nml_group%item(j), constraint, err)

    case ('point')
      if (.not. ref%has_sub) then
        call tao_nml_err (nml_group%item(j), 'POINT IS AN ARRAY SO A SUBSCRIPT IS NEEDED')
        cycle
      endif
      if (ref%isub < 1 .or. ref%isub > size(point)) then
        call tao_nml_err (nml_group%item(j), 'POINT SUBSCRIPT OUT OF RANGE')
        cycle
      endif

      if (ref%rest == '') then    ! Whole structure assignment. Eg: point(3) = 1.0, 2.0
        call tao_nml_split_values (nml_group%item(j)%value, nml_val, n_val, err, why)
        if (err) then
          call tao_nml_err (nml_group%item(j), why)
          cycle
        endif
        do k = 1, n_val
          if (nml_val(k) == '') cycle
          call tao_res_tao_building_wall_point_struct_slot (point(ref%isub), '', k, ptr, err, why)
          if (.not. err) call tao_set_ptr_value (ptr, nml_val(k), err, why)
          if (err) then
            call tao_nml_err (nml_group%item(j), why)
            exit
          endif
        enddo
      else
        call tao_res_tao_building_wall_point_struct (point(ref%isub), ref%rest, ptr, err, why)
        if (err) then
          call tao_nml_err (nml_group%item(j), why)
        else
          call tao_nml_value_set (nml_group%item(j), ptr, err)
        endif
      endif

    case default
      call tao_nml_unknown (nml_group%item(j), 'building_wall_section', err)
    end select
  enddo

  s%building_wall%section(i)%name = name

  select case (constraint)
  case ('left_side', 'right_side', 'none')
    s%building_wall%section(i)%constraint = constraint
  case default
    call out_io (s_error$, r_name, 'BAD "CONSTRAINT" FOR BUILDING WALL: ' // constraint)
    return
  end select

  do j = size(point), 1, -1
    if (point(j)%x == real_garbage$) cycle
    allocate (s%building_wall%section(i)%point(j))
    exit
  enddo

  if (.not. allocated(s%building_wall%section(i)%point)) then
    call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_SECTION NAMELIST POINT ARRAY', &
                                   'IN BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', i_array = [i])
    return
  endif

  if (point(1)%radius /= 0) then
    call out_io (s_fatal$, r_name, 'ERROR IN POINT ARRAY OF BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', &
                                   'FIRST POINT HAS NON-ZERO RADIUS', i_array = [i])
    return
  endif

  pt => s%building_wall%section(i)%point

  do j = 1, size(pt)
    
    pt(j)%x = point(j)%x
    pt(j)%z = point(j)%z
    pt(j)%radius = point(j)%radius
    if (pt(j)%radius /= 0) then
      x_mid = (pt(j)%x + pt(j-1)%x) / 2; z_mid = (pt(j)%z + pt(j-1)%z) / 2 
      dx    = (pt(j)%x - pt(j-1)%x) / 2; dz    = (pt(j)%z - pt(j-1)%z) / 2 
      a2 = (pt(j)%radius**2 - dx**2 - dz**2) / (dx**2 + dz**2)
      if (a2 < 0) then
        call out_io (s_fatal$, r_name, 'ERROR IN POINT ARRAY OF BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', &
                                       'WALL POINTS TOO FAR APART FOR CIRCLE AT POINT \i0\ ', i_array = [i, j])
        return
      endif
      a = sqrt(a2)
      if (pt(j)%radius < 0) a = -a
      pt(j)%x_center = x_mid - a * dz
      pt(j)%z_center = z_mid + a * dx
    endif
  enddo

enddo

close (iu)

end subroutine tao_init_building_wall
