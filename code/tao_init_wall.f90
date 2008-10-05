!+
! Subroutine tao_init_wall (wall_file)
!
! Subroutine to initialize the tunnel wall.
!
! Input:
!   wall_file  -- character(*): file name containing the wall definition.
!                   If blank then there is no wall.
!-

subroutine tao_init_wall (wall_file)

  use tao_mod
  use tao_input_struct

  implicit none

  type tao_wall_point_input
    character(8) type
    real(rp) x, z, r, theta1, theta2
  end type

  type (tao_wall_point_input) point(100)
  type (tao_wall_point_struct), pointer :: pt(:)

  integer i, j, iu, ios, n_wall

  character(*) wall_file
  character(200) complete_file_name
  character(20) :: r_name = 'tao_init_wall'
  character(8) side

  namelist / tunnel_wall / side, point

! Open file

  if (wall_file == '') return

  call tao_open_file ('TAO_INIT_DIR', wall_file, iu, complete_file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // complete_file_name)
  if (iu == 0) then
    call out_io (s_fatal$, r_name, 'ERROR OPENING TUNNEL WALL FILE. WILL EXIT HERE...')
    call err_exit
  endif

! Count the number of walls

  n_wall = 0
  do
    read (iu, nml = tunnel_wall, iostat = ios)
    if (ios > 0) then
      call out_io (s_fatal$, r_name, 'ERROR READING TUNNEL_WALL NAMELIST')
      do   ! Generate informational message
        read (iu, nml = tunnel_wall)
      enddo
    endif
    if (ios < 0) exit
    n_wall = n_wall + 1
  enddo

  call out_io (s_blank$, r_name, '  Number of tunnel walls: \I2\ ', n_wall)

  allocate (s%wall(n_wall))
  if (n_wall == 0) return ! no walls

! Now transfer the information

  rewind (iu)
  do i = 1, n_wall

    side = ''
    point%type = ''
    read (iu, nml = tunnel_wall, iostat = ios)

    if (side == 'left') then
      s%wall(i)%side = left_side$
    elseif (side == 'right') then
      s%wall(i)%side = right_side$
    else
      call out_io (s_error$, r_name, 'BAD "SIDE" FOR TUNNEL WALL: ' // side)
      call err_exit
    endif

    do j = size(point), 1, -1
      if (point(j)%type == '') cycle
      if (.not. allocated(s%wall(i)%point)) then
        allocate (s%wall(i)%point(j))
        pt => s%wall(i)%point
      endif
      if (point(j)%type == 'point') then
        pt(j)%type = point$
      elseif (point(j)%type == 'arc') then
        pt(j)%type = arc$
      else
        call out_io (s_error$, r_name, 'BAD POINT "TYPE" FOR TUNNEL WALL: ' // point%type)
        call err_exit
      endif
      pt(j)%x = point(j)%x
      pt(j)%z = point(j)%z
      pt(j)%r = point(j)%r
      pt(j)%theta1 = point(j)%theta1
      pt(j)%theta2 = point(j)%theta2
      if (pt(j)%theta2 < pt(j)%theta1) pt(j)%theta2 = pt(j)%theta2 + &
                twopi * int((pt(j)%theta1 - pt(j)%theta2) / twopi)
    enddo

  enddo

  close (iu)

end subroutine tao_init_wall
