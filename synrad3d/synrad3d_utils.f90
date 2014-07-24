module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod
use capillary_mod

type sr3d_wall_section_input
  real(rp) s                      ! Longitudinal position.
  character(40) name              ! Name of setion
  character(60) basic_shape       ! "elliptical", "rectangular", "triangular:xxx", or "gen_shape:xxx"
  real(rp) width2                 ! Half width ignoring antechamber.
  real(rp) height2                ! Half height ignoring antechamber.
  real(rp) width2_plus            ! Distance from pipe center to +x side edge.
  real(rp) ante_height2_plus      ! Antechamber half height on +x side of the wall
  real(rp) width2_minus           ! Distance from pipe center -x side edge.
  real(rp) ante_height2_minus     ! Antechamber half height on -x side of the wall
end type

type surface_input
  character(40) name
  logical is_local
end type

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_read_wall_file (wall_file, s_lat, geometry, wall, err_flag)
!
! Routine to check the vacuum chamber wall for problematic values.
! Also compute some wall parameters
!
! Input:
!   wall_file   -- character(*): Name of the wall file.
!   s_lat       -- Real(rp): Lattice length
!   geometry    -- Integer: Type of lattice. open$ or closed$
!
! Output:
!   wall      -- sr3d_wall_struct: Wall structure with computed parameters.
!   err_flag  -- logical, optional: Set true if there is a problem
!-

subroutine sr3d_read_wall_file (wall_file, s_lat, geometry, wall, err_flag)

implicit none

! Needed since Fortran does not allow pointers to be part of a namelist

type (sr3d_wall_struct), target :: wall
type (lat_struct) lat
type (sr3d_wall_section_struct), pointer :: sec, sec0, sec2
type (sr3d_wall_section_struct), allocatable :: temp_section(:)
type (sr3d_wall_section_struct) ref_section
type (sr3d_wall_section_input) section
type (wall3d_vertex_struct) v(100)
type (wall3d_section_struct), pointer :: wall3d_section
type (sr3d_multi_section_struct), pointer :: m_sec
type (surface_input) surface
type (photon_reflect_surface_struct), pointer :: surface_ptr  => null()

real(rp) ix_vertex_ante(2), ix_vertex_ante2(2), s_lat
real(rp) rad, radius(4), area, max_area, cos_a, sin_a, angle, dr_dtheta

integer i, j, k, im, n, ix, iu, n_wall_section_max, ios, n_shape, n_surface, n_repeat
integer m_max, n_add, geometry

character(28), parameter :: r_name = 'sr3d_read_wall_file'
character(40) name
character(200) reflectivity_file, file
character(*) wall_file

logical, optional :: err_flag
logical err, multi_local, in_ante
logical, allocatable ::  is_local(:) 


namelist / section_def / section, surface
namelist / gen_shape_def / name, v, ix_vertex_ante, ix_vertex_ante2
namelist / surface_def / name, reflectivity_file

! Open file

wall%has_triangular_sections = .false.

iu = lunget()
call fullfilename (wall_file, file)
open (iu, file = file, status = 'old')

! Read in the surface info

rewind(iu)
n_surface = 0
do
  read (iu, nml = surface_def, iostat = ios)
  if (ios > 0) then ! error
    print *, 'ERROR READING SURFACE_DEF NAMELIST'
    rewind (iu)
    do
      read (iu, nml = surface_def) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  n_surface = n_surface + 1
enddo

if (allocated(wall%surface)) deallocate(wall%surface)
allocate (wall%surface(n_surface+1))
wall%surface(1)%reflectivity_file = ''
wall%surface(1)%descrip = 'default'
call photon_reflection_std_surface_init(wall%surface(1))

rewind(iu)
do i = 2, n_surface+1
  read (iu, nml = surface_def, iostat = ios)
  call read_surface_reflection_file (reflectivity_file, wall%surface(i))
  wall%surface(i)%descrip = name
  wall%surface(i)%reflectivity_file = reflectivity_file
enddo

! Read multi_section

call sr3d_read_wall_multi_section (iu, wall)

! Get wall info
! First count the cross-section number

rewind(iu)
n_wall_section_max = -1
do
  read (iu, nml = section_def, iostat = ios)
  if (ios > 0) then ! error
    rewind (iu)
    do
      read (iu, nml = section_def) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  n_wall_section_max = n_wall_section_max + 1
enddo

print *, 'number of wall cross-sections read:', n_wall_section_max + 1
if (n_wall_section_max < 1) then
  print *, 'NO WALL SPECIFIED. WILL STOP HERE.'
  call err_exit
endif

if (allocated(wall%section)) deallocate(wall%section)
allocate (wall%section(0:n_wall_section_max), is_local(0:n_wall_section_max))
wall%n_section_max = n_wall_section_max

! Now transfer info from the file to the wall%section array

rewind (iu)
do i = 0, n_wall_section_max
  section%basic_shape = ''
  section%width2 = -1
  section%height2 = -1
  section%ante_height2_plus = -1
  section%ante_height2_minus = -1
  section%width2_plus = -1
  section%width2_minus = -1
  section%name = ''
  surface%name = ''
  surface%is_local = .false.

  read (iu, nml = section_def)

  sec => wall%section(i)
  if (section%basic_shape(1:1) == '+') then
    section%s = wall%section(i-1)%s + section%s
    section%basic_shape = section%basic_shape(2:)
  endif

  sec = sr3d_wall_section_struct(section%name, &
          section%basic_shape, '', section%s, section%width2, section%height2, &
          section%width2_plus, section%ante_height2_plus, &
          section%width2_minus, section%ante_height2_minus, &
          -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, null(), null())

  ix = index(section%basic_shape, ':')
  if (ix /= 0) then
    sec%basic_shape = section%basic_shape(1:ix-1)
    sec%shape_name  = section%basic_shape(ix+1:)
  endif

  call sr3d_associate_surface (sec%surface, surface%name, wall%surface)
  is_local(i) = surface%is_local

enddo

! Get the gen_shape info

n_shape = 0
rewind(iu)
do
  read (iu, nml = gen_shape_def, iostat = ios)
  if (ios > 0) then ! If error
    print *, 'ERROR READING GEN_SHAPE_DEF NAMELIST.'
    rewind (iu)
    do
      read (iu, nml = gen_shape_def) ! Generate error message
    enddo
  endif
  if (ios < 0) exit  ! End of file
  n_shape = n_shape + 1
enddo

if (allocated(wall%gen_shape)) deallocate(wall%gen_shape)
allocate (wall%gen_shape(n_shape))

rewind(iu)
do i = 1, n_shape
  ix_vertex_ante = 0
  ix_vertex_ante2 = 0
  v = wall3d_vertex_struct()
  name = ''
  read (iu, nml = gen_shape_def, iostat = ios)

  if (name == '') then
    print *, 'GEN_SHAPE_DEF DOES NOT HAVE A NAME! GEN_SHAPE NUMBER:', i
    call err_exit
  endif

  do j = 1, i-1
    if (wall%gen_shape(j)%name == name) then
      print *, 'TWO GEN_SHAPE_DEFS HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo
  wall%gen_shape(i)%name = name

  ! Count number of vertices and calc angles.

  wall3d_section => wall%gen_shape(i)%wall3d_section
  do n = 1, size(v)
    if (v(n)%x == 0 .and. v(n)%y == 0 .and. v(n)%radius_x == 0) exit
  enddo

  if (any(v(n:)%x /= 0) .or. any(v(n:)%y /= 0) .or. &
      any(v(n:)%radius_x /= 0) .or. any(v(n:)%radius_y /= 0)) then
    print *, 'MALFORMED GEN_SHAPE:', name
    call err_exit
  endif

  allocate(wall3d_section%v(n-1))
  wall3d_section%v = v(1:n-1)
  wall3d_section%n_vertex_input = n-1    

  call wall3d_section_initializer (wall3d_section, err)
  if (err) then
    print *, 'ERROR AT GEN_SHAPE: ', trim(name)
    call err_exit
  endif

  wall%gen_shape(i)%ix_vertex_ante = ix_vertex_ante
  if (ix_vertex_ante(1) > 0 .or. ix_vertex_ante(2) > 0) then
    if (ix_vertex_ante(1) < 1 .or. ix_vertex_ante(1) > size(wall3d_section%v) .or. &
        ix_vertex_ante(2) < 1 .or. ix_vertex_ante(2) > size(wall3d_section%v)) then
      print *, 'ERROR IN IX_VERTEX_ANTE:', ix_vertex_ante
      print *, '      FOR GEN_SHAPE: ', trim(name)
      call err_exit
    endif
  endif

  wall%gen_shape(i)%ix_vertex_ante2 = ix_vertex_ante2
  if (ix_vertex_ante2(1) > 0 .or. ix_vertex_ante2(2) > 0) then
    if (ix_vertex_ante2(1) < 1 .or. ix_vertex_ante2(1) > size(wall3d_section%v) .or. &
        ix_vertex_ante2(2) < 1 .or. ix_vertex_ante2(2) > size(wall3d_section%v)) then
      print *, 'ERROR IN IX_VERTEX_ANTE2:', ix_vertex_ante2
      print *, '      FOR GEN_SHAPE: ', trim(name)
      call err_exit
    endif
  endif

enddo

close (iu)

! Expand multi_sections

i = -1
outer: do 
  i = i + 1
  if (i > wall%n_section_max) exit
  ref_section = wall%section(i)
  if (ref_section%basic_shape /= 'multi_section') cycle
  n_repeat = nint(ref_section%width2) 

  if (n_repeat < 0) then
    print *, 'ERROR: MULTI_SECTION DOES NOT HAVE THE REPEAT COUNT SET.'
    call err_exit
  endif

  do j = 1, size(wall%multi_section)
    m_sec => wall%multi_section(j)
    if (ref_section%shape_name /= m_sec%name) cycle
    m_max = ubound(m_sec%section, 1)

    if (m_sec%section(m_max)%name == 'closed') then
      n_add = n_repeat * m_max + 1
    elseif (m_sec%section(m_max)%name == 'open') then
      n_add = n_repeat * m_max
    else
      print *, 'ERROR: LAST SECTION IN MULTI_SECTION IS NOT "closed" NOR "open"'
      call err_exit
    endif

    call move_alloc (wall%section, temp_section)
    n = wall%n_section_max
    allocate (wall%section(0:n+n_add-1))
    wall%section(0:i-1) = temp_section(0:i-1)
    wall%section(i+n_add:n+n_add-1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall%n_section_max = n + n_add - 1

    call re_allocate2 (is_local, 0, n+n_add-1)
    is_local(i+n_add:n+n_add-1) = is_local(i+1:n)
    is_local(i:i+n_add-1) = is_local(i)

    do k = 1, n_repeat
      do im = 0, m_max - 1
        sec2 => wall%section(i+(k-1)*m_max+im)
        sec2 = m_sec%section(im)
        sec2%s = ref_section%s + (k-1) * m_sec%section(m_max)%s + m_sec%section(im)%s
      enddo
    enddo
    if (m_sec%section(m_max)%name == 'closed') then
      sec2 => wall%section(i+n_repeat*m_max)
      sec2 = m_sec%section(1)
      sec2%s = ref_section%s + n_repeat * m_sec%section(m_max)%s
    endif      

    cycle outer

  enddo

  print *, 'CANNOT FIND MATCHING MULTI_SECTION NAME: ', trim(ref_section%shape_name)
  call err_exit

enddo outer

! point to gen_shapes

section_loop: do i = 0, wall%n_section_max
  sec => wall%section(i)
  if (sec%basic_shape == 'triangular') wall%has_triangular_sections = .true.
  if (sec%basic_shape /= 'gen_shape' .and. sec%basic_shape /= 'triangular') cycle
  do j = 1, size(wall%gen_shape)
    if (sec%shape_name /= wall%gen_shape(j)%name) cycle
    sec%gen_shape => wall%gen_shape(j)
    cycle section_loop
  enddo

  print *, 'CANNOT FIND MATCHING SHAPE FOR: ', trim(sec%shape_name)
  call err_exit

enddo section_loop

! Last section has s adjusted to match the lattice length.

if (abs(wall%section(wall%n_section_max)%s - s_lat) > 0.01) then
  call out_io (s_info$, r_name, &
        'Wall ends at: \f12.4\ ', &
        'And not at lattice end of: \f12.4\ ', &
        '[But last point is always adjusted to have s = s_lat]', &
        r_array = [wall%section(wall%n_section_max)%s, s_lat])
endif

wall%section(wall%n_section_max)%s = s_lat
wall%geometry = geometry

! Checks

do i = 0, wall%n_section_max
  sec => wall%section(i)

  ! Check s ordering

  if (i > 0) then
    if (sec%s == wall%section(i-1)%s) sec%s = sec%s + 1000*sr3d_params%significant_length
    if (sec%s < wall%section(i-1)%s) then
      call out_io (s_fatal$, r_name, &
                'WALL%SECTION(i)%S: \f0.4\ ', &
                '    IS LESS THAN SECTION(i-1)%S: \f0.4\ ', &
                '    FOR I = \i0\ ', &
                r_array = [sec%s, wall%section(i-1)%s], i_array = [i])
      call err_exit
    endif
  endif

  ! Check %basic_shape

  if (.not. any(sec%basic_shape == ['elliptical    ', 'rectangular   ', &
                                    'gen_shape     ', 'triangular    '])) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%SECTION(i)%BASIC_SHAPE: ' // sec%basic_shape, &
              '    FOR I = \i0\ ', i_array = [i])
    call err_exit
  endif

  ! Gen_shape and triangular checks

  if (sec%basic_shape == 'gen_shape' .or. sec%basic_shape == 'triangular') then
    if (.not. associated (sec%gen_shape)) then
      call out_io (s_fatal$, r_name, 'BAD SHAPE ASSOCIATION')
      call err_exit
    endif
    if (sec%basic_shape == 'gen_shape') cycle
    if (i == 0) cycle

    sec0 => wall%section(i-1)
    if (sec0%basic_shape /= 'gen_shape' .and. sec0%basic_shape /= 'triangular') then
      call out_io (s_fatal$, r_name, &
              'BASIC_SHAPE FOR SECTION PRECEEDING "triangular" SECTION MUST BE ', &
              '"gen_shape" OR "triangular" SECTION NUMBER \i0\ ', i_array = [i])
      call err_exit
    endif

    if (size(sec0%gen_shape%wall3d_section%v) /= size(sec%gen_shape%wall3d_section%v)) then
      call out_io (s_fatal$, r_name, &
              '"triangular" CONSTRUCT MUST HAVE THE SAME NUMBER OF VERTEX POINTS ON', &
              'SUCCESIVE CROSS-SECTIONS  \2i0\ ', i_array = [i-1, i])
      call err_exit
    endif

    cycle
  endif

  ! Checks for everything else

  if (sec%width2 <= 0) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%SECTION(i)%WIDTH2: \f0.4\ ', &
              '    FOR I = \i0\ ', r_array = [sec%width2], i_array = [i])
    call err_exit
  endif

  if (sec%width2 <= 0) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%SECTION(i)%HEIGHT2: \f0.4\ ', &
              '    FOR I = \i0\ ', r_array = [sec%height2], i_array = [i])
    call err_exit
  endif

  ! +x side check

  if (sec%ante_height2_plus < 0 .and.sec%width2_plus > 0) then
    if (sec%width2_plus > sec%width2) then
      call out_io (s_fatal$, r_name, &
              'WITHOUT AN ANTECHAMBER: WALL%SECTION(i)%WIDTH2_PLUS \f0.4\ ', &
              '    MUST BE LESS THEN WIDTH2 \f0.4\ ', &
              '    FOR I = \i0\ ', &
              r_array = [sec%width2_plus, sec%width2], i_array = [i])
      call err_exit
    endif
  endif

  ! -x side check

  if (sec%ante_height2_minus < 0 .and. sec%width2_minus > 0) then
    if (sec%width2_minus > sec%width2) then
      call out_io (s_fatal$, r_name, &
              'WITHOUT AN ANTECHAMBER: WALL%SECTION(i)%WIDTH2_MINUS \f0.4\ ', &
              '    MUST BE LESS THEN WIDTH2 \f0.4\ ', &
              '    FOR I = \i0\ ', &
              r_array = [sec%width2_minus, sec%width2], i_array = [i])
      call err_exit
    endif
  endif

enddo

! If circular lattice then start and end shapes must match

if (wall%geometry == closed$) then
  sec0 => wall%section(0)
  sec  => wall%section(wall%n_section_max)
  if (sec0%basic_shape /= sec%basic_shape .or. sec0%width2 /= sec%width2 .or. sec0%height2 /= sec%height2 .or. &
        sec0%ante_height2_plus /= sec%ante_height2_plus .or. sec0%width2_plus /= sec%width2_plus .or. &
        sec0%ante_height2_minus /= sec%ante_height2_minus .or. sec0%width2_minus /= sec%width2_minus) then
      call out_io (s_fatal$, r_name, &
              'FOR A "CLOSED" LATTICE THE LAST WALL CROSS-SECTION MUST BE THE SAME AS THE FIRST.')
      call err_exit
  endif
endif

! computations

do i = 0, wall%n_section_max
  sec => wall%section(i)

  ! +x side computation...
  ! If ante_height2_plus > 0 --> Has +x antechamber

  if (sec%ante_height2_plus > 0) then
    if (sec%basic_shape == 'elliptical') then
      sec%ante_x0_plus = sec%width2 * sqrt (1 - (sec%ante_height2_plus / sec%height2)**2)
    else
      sec%ante_x0_plus = sec%width2
    endif

    if (sec%width2_plus <= sec%ante_x0_plus) then
      call out_io (s_fatal$, r_name, &
              'WITH AN ANTECHAMBER: WALL%SECTION(i)%WIDTH2_PLUS \f0.4\ ', &
              '    MUST BE GREATER THEN: \f0.4\ ', &
              '    FOR I = \i0\ ', &
              r_array = [sec%width2_plus, sec%ante_x0_plus], i_array = [i])
      call err_exit
    endif

  ! if width2_plus > 0 (and ante_height2_plus < 0) --> beam stop

  elseif (sec%width2_plus > 0) then
    if (sec%basic_shape == 'elliptical') then
      sec%y0_plus = sec%height2 * sqrt (1 - (sec%width2_plus / sec%width2)**2)
    else
      sec%y0_plus = sec%height2
    endif
  endif

  ! -x side computation

  if (sec%ante_height2_minus > 0) then
    if (sec%basic_shape == 'elliptical') then
      sec%ante_x0_minus = sec%width2 * sqrt (1 - (sec%ante_height2_minus / sec%height2)**2)
    else
      sec%ante_x0_minus = sec%width2
    endif

    if (sec%width2_minus <= sec%ante_x0_minus) then
      call out_io (s_fatal$, r_name, &
              'WITH AN ANTECHAMBER: WALL%SECTION(i)%WIDTH2_MINUS \f0.4\ ', &
              '    MUST BE GREATER THEN: \f0.4\ ', &
              '    FOR I = \i0\ ', &
              r_array = [sec%width2_minus, sec%ante_x0_minus], i_array = [i])

      call err_exit
    endif

  elseif (sec%width2_minus > 0) then
    if (sec%basic_shape == 'elliptical') then
      sec%y0_minus = sec%height2 * sqrt (1 - (sec%width2_minus / sec%width2)**2)
    else
      sec%y0_minus = sec%height2
    endif
  endif

enddo

! Calculate largest "safe" box for each section.
! For triangular mesh, this is a complicated calc so, for now, the safe box calc is not 
! implemented in this instance. 

wall%gen_shape(:)%x_safe = 0
wall%gen_shape(:)%y_safe = 0

do i = 0, wall%n_section_max
  sec => wall%section(i)

  if (sec%basic_shape == 'triangular') cycle

  if (associated(sec%gen_shape)) then
    if (sec%gen_shape%x_safe > 0) then
      sec%x_safe = sec%gen_shape%x_safe
      sec%y_safe = sec%gen_shape%y_safe
      cycle
    endif
  endif

  max_area = 0
  do j = 1, 39
    angle = j * pi / 39
    cos_a = cos(angle); sin_a = sin(angle)
    call sr3d_wall_section_params (sec,  cos_a,  sin_a, radius(1), dr_dtheta, in_ante)
    call sr3d_wall_section_params (sec,  cos_a, -sin_a, radius(2), dr_dtheta, in_ante)
    call sr3d_wall_section_params (sec, -cos_a,  sin_a, radius(3), dr_dtheta, in_ante)
    call sr3d_wall_section_params (sec, -cos_a, -sin_a, radius(4), dr_dtheta, in_ante)
    rad = minval(radius)
    area = rad**2 * cos_a * sin_a
    if (area > max_area) then
      sec%x_safe = 0.999 * rad * cos_a
      sec%y_safe = 0.999 * rad * sin_a
      max_area = area
    endif
  enddo

  if (associated(sec%gen_shape)) then
    sec%gen_shape%x_safe = sec%x_safe
    sec%gen_shape%y_safe = sec%y_safe
  endif

enddo

! Surface info

surface_ptr => wall%surface(1)  ! Default surface

do i = wall%n_section_max, 0, -1
  sec => wall%section(i)

  if (associated(sec%surface)) then
    if (.not. is_local(i)) surface_ptr => sec%surface
  else
    sec%surface => surface_ptr
  endif

enddo

end subroutine sr3d_read_wall_file 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_associate_surface (surface_ptr, surface_name, surfaces)
!
! Routine to assocaite a surface pointer with one of an array os surfaces.
!
! Input:
!   surface_name  -- Character(*): Name of surface.
!   surfaces(:)   -- photon_reflect_surface_struct: Array of surfaces.
!
! Output:
!   surface_ptr   -- photon_reflect_surface_struct, pointer: pointer to a surface.
!-

subroutine sr3d_associate_surface (surface_ptr, surface_name, surfaces)

implicit none

type (photon_reflect_surface_struct), pointer :: surface_ptr
type (photon_reflect_surface_struct), target :: surfaces(:)

character(*) surface_name

integer i

!

nullify(surface_ptr)
if (surface_name == '') return

do i = 1, size(surfaces)
  if (surfaces(i)%descrip == surface_name) then
    surface_ptr => surfaces(i)
    return
  endif
enddo

print *, 'NO SURFACE CORRESPONDING TO: ', trim(surface_name)
call err_exit

end subroutine sr3d_associate_surface

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_read_wall_multi_section (iu, wall)
!
! Routine to read in the multi_section part of the wall definition.
!
! Input:
!   iu   -- integer: File unit number.
!
! Output:
!   wall -- sr3d_wall_struct: Wall structure with computed parameters.
!-

subroutine sr3d_read_wall_multi_section (iu, wall)

type (sr3d_wall_struct), target :: wall
type (sr3d_wall_section_input) section(0:100)

integer i, j, iu, ix, n_multi, n_section, ios

character(40) surface, name

namelist / multi_section_def / name, section

! Count multi_sections

rewind(iu)
n_multi = 0
do
  read (iu, nml = multi_section_def, iostat = ios)
  if (ios > 0) then ! error
    rewind (iu)
    do
      read (iu, nml = multi_section_def) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  n_multi = n_multi + 1
enddo

if (allocated(wall%multi_section)) deallocate (wall%multi_section)
allocate (wall%multi_section(n_multi))

! Read each multi_section

rewind(iu)
do i = 1, n_multi
  section%basic_shape = ''
  section%ante_height2_plus = -1
  section%ante_height2_minus = -1
  section%width2_plus = -1
  section%width2_minus = -1
  section%name = ''
  name = ''
  read (iu, nml = multi_section_def)

  if (name == '') then
    print *, 'MULTI_SECTION_DEF DOES NOT HAVE A NAME! MULTI_SECTION_DEF NUMBER:', i
    call err_exit
  endif

  do j = 1, i-1
    if (wall%multi_section(j)%name == name) then
      print *, 'TWO GEN_SHAPE_DEFS HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo

  wall%multi_section(i)%name = name

  n_section = count(section%basic_shape /= '')
  if (any(section(n_section:)%basic_shape /= '')) then
    print *, 'CONFUSED MULTI_SECTION_DEF: ', trim(name)
    call err_exit
  endif
  allocate (wall%multi_section(i)%section(0:n_section-1))

  do j = 0, n_section-1
    ix = index(section(j)%basic_shape, ':')
    if (ix == 0) ix = len_trim(section(j)%basic_shape) + 1
    wall%multi_section(i)%section(j) = sr3d_wall_section_struct(section(j)%name, &
          section(j)%basic_shape(:ix-1), section(j)%basic_shape(ix+1:), &
          section(j)%s, section(j)%width2, section(j)%height2, section(j)%width2_plus, &
          section(j)%ante_height2_plus, section(j)%width2_minus, section(j)%ante_height2_minus, &
          -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, null(), null())
  enddo
enddo

end subroutine sr3d_read_wall_multi_section

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
!
! Routine to get the parameters at a photon emission point.
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   lat       -- lat_struct with twiss propagated and mat6s made.
!   orb(0:*)  -- coord_struct: orbit of particles to use as source of ray.
!   ix_ele    -- integer: index of lat element to start ray from.
!   s_offset  -- real(rp): Distance from beginning of element to the point where the photon is emitted.
!
! Output:
!   ele_here  -- ele_struct: Twiss parameters at emission point.
!   orb_here  -- coord_struct: Beam center coords at emission point.
!   gx        -- Real(rp): Horizontal 1/bending_radius.
!   gy        -- Real(rp): Vertical 1/bending_radius.
!-

subroutine sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

use em_field_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele
type (ele_struct) ele_here
type (sr3d_photon_coord_struct) :: photon
type (em_field_struct) :: field

real(rp) s_offset, k_wig, g_max, l_small, gx, gy
real(rp), save :: s_old_offset = 0

integer ix_ele

logical err
logical, save :: init_needed = .true.

! Init

if (init_needed) then
  call init_ele (ele_here)
  init_needed = .false.
endif

ele  => lat%ele(ix_ele)

! Calc the photon's initial twiss values.
! Tracking through a wiggler can take time so use twiss_and_track_intra_ele to
!   minimize the length over which we track.

if (ele_here%ix_ele /= ele%ix_ele .or. ele_here%ix_branch /= ele%ix_branch .or. s_old_offset > s_offset) then
  ele_here = lat%ele(ix_ele-1)
  ele_here%ix_ele = ele%ix_ele
  ele_here%ix_branch = ele%ix_branch
  orb_here = orb(ix_ele-1)
  s_old_offset = 0
endif

call twiss_and_track_intra_ele (ele, lat%param, s_old_offset, s_offset, .true., .true., &
                                            orb_here, orb_here, ele_here, ele_here, err)
if (err) call err_exit
s_old_offset = s_offset

! Calc the photon's g_bend value (inverse bending radius at src pt) 

select case (ele%key)
case (sbend$)  

  ! sbends are easy
  gx = ele%value(g$) + ele%value(g_err$)
  gy = 0
  if (ele%value(roll$) /= 0) then
    gy = gx * sin(ele%value(roll$))
    gx = gx * cos(ele%value(roll$))
  endif

case (quadrupole$, sol_quad$, elseparator$, sad_mult$)

  ! for quads or sol_quads, get the bending radius
  ! from the change in x' and y' over a small 
  ! distance in the element

  l_small = 1e-2      ! something small
  ele_here%value(l$) = l_small
  call make_mat6 (ele_here, lat%param, orb_here, orb_here, .true.)
  call track1 (orb_here, ele_here, lat%param, orb1)
  orb1%vec = orb1%vec - orb_here%vec
  gx = orb1%vec(2) / l_small
  gy = orb1%vec(4) / l_small

case (wiggler$)

  if (ele%sub_key == periodic_type$) then

    ! for periodic wigglers, get the max g_bend from 
    ! the max B field of the wiggler, then scale it 
    ! by the cos of the position along the poles

    k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))
    g_max = c_light * ele%value(b_max$) / (ele%value(p0c$) * (1 + orb_here%vec(6)))
    gx = g_max * sin (k_wig * s_offset)
    gy = 0
    orb_here%vec(1) = (g_max / k_wig) * sin (k_wig * s_offset)
    orb_here%vec(2) = (g_max / k_wig) * cos (k_wig * s_offset)

  else

    ! for mapped wigglers, find the B field at the source point
    ! Note: assumes particles are relativistic!!

    call em_field_calc (ele_here, lat%param, ele_here%value(l$), 0.0_rp, orb_here, .false., field)
    gx = field%b(2) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))
    gy = field%b(1) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))

  endif

case default

  ! for mapped wigglers, find the B field at the source point
  ! Note: assumes particles are relativistic!!

  call em_field_calc (ele_here, lat%param, ele_here%value(l$), 0.0_rp, orb_here, .false., field)
  gx = field%b(2) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))
  gy = field%b(1) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))

end select

end subroutine sr3d_get_emission_pt_params 


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon)
!
! subroutine sr3d_to initialize a new photon
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   ele_here  -- Ele_struct: Element emitting the photon. Emission is at the exit end of the element.
!   orb_here  -- coord_struct: orbit of particles emitting the photon.
!   gx, gy    -- Real(rp): Horizontal and vertical bending strengths.
!   emit_a    -- Real(rp): Emittance of the a-mode.
!   emit_b    -- Real(rp): Emittance of the b-mode.
!   photon_direction 
!             -- Integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- photon_coord_struct: Generated photon.
!-

subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, p_orb)

implicit none

type (ele_struct), target :: ele_here
type (coord_struct) :: orb_here, orb_init
type (sr3d_photon_coord_struct) :: p_orb
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_e, gx, gy, g_tot, gamma, v2
real(rp) r(3), vec(4), v_mat(4,4)

integer photon_direction

! Get photon energy and "vertical angle".

call convert_total_energy_to (ele_here%value(E_tot$), ele_here%branch%param%particle, gamma) 
call bend_photon_init (gx, gy, gamma, orb_init)
p_orb%energy = orb_init%p0c
p_orb%vec = orb_init%vec
p_orb%s   = orb_init%s

! Offset due to finite beam size

call ran_gauss(r)
t => ele_here%a
vec(1:2) = (/ sqrt(t%beta*emit_a) * r(1)                    + t%eta  * sig_e * r(3), &
              sqrt(emit_a/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_e * r(3) /)

call ran_gauss(r)
t => ele_here%b
vec(3:4) = (/ sqrt(t%beta*emit_b) * r(1)                    + t%eta  * sig_e * r(3), &
              sqrt(emit_b/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_e * r(3) /)

call make_v_mats (ele_here, v_mat)

p_orb%vec(1:4) = p_orb%vec(1:4) + matmul(v_mat, vec)

! Offset due to non-zero orbit.

p_orb%vec(1:4) = p_orb%vec(1:4) + orb_here%vec(1:4)

! Longitudinal position

p_orb%s = ele_here%s

! Above equations are valid in the small angle limit.
! Sometimes a large-angle photon is generated so make sure
! there is no problem with the sqrt() evaluation.

v2 = p_orb%vec(2)**2 + p_orb%vec(4)**2
if (v2 >= 0.99) then
  p_orb%vec(2) = p_orb%vec(2) * 0.99 / v2
  p_orb%vec(4) = p_orb%vec(4) * 0.99 / v2
  v2 = p_orb%vec(2)**2 + p_orb%vec(4)**2
endif

p_orb%vec(6) = photon_direction * sqrt(1 - v2)

end subroutine sr3d_emit_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_d_radius (p_orb, wall, d_radius, lat, dw_perp, in_antechamber, check_safe)
!
! Routine to calculate the (transverse) radius of the photon  relative to the wall.
! Optionally can also caluclate the outwrd normal vector perpendicular to the wall.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall       -- sr3d_wall_struct: Wall
!   s          -- Real(rp): Longitudinal position.
!   lat        -- Lat_struct, optional: Lattice. Only needed when dw_perp is calculated.
!   check_safe -- logical, optional: If True, check if photon is safely in the "safe" box far from
!                   the wall. This is used to speed up computations when d_radius value is not needed.
!                   Default is False.
!
! Output:
!   d_radius       -- real(rp): r_photon - r_wall
!   dw_perp(3)     -- real(rp), optional: Outward normal vector perpendicular to the wall.
!   in_antechamber -- Logical, optional: At antechamber wall?
!-

Subroutine sr3d_photon_d_radius (p_orb, wall, d_radius, lat, dw_perp, in_antechamber, check_safe)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_photon_coord_struct), target :: p_orb
type (lat_struct), optional, target :: lat
type (ele_struct), pointer :: ele

real(rp) d_radius
real(rp), optional :: dw_perp(3)
real(rp) radius0, radius1, f, cos_ang, sin_ang, r_photon, disp
real(rp) dr0_dtheta, dr1_dtheta, pt0(3), pt1(3), pt2(3), dp1(3), dp2(3)

integer ix, ix_ele

logical, optional :: in_antechamber, check_safe
logical in_ante0, in_ante1

! If the photon's position (abs(x), abs(y)) is within the box (x_safe, y_safe) then the photon
! is guaranteed to be within the vacuum chamber.

call sr3d_get_wall_index (p_orb, wall, ix)

if (logic_option(.false., check_safe)) then
  if (abs(p_orb%vec(1)) < min(wall%section(ix)%x_safe, wall%section(ix+1)%x_safe) .and. &
      abs(p_orb%vec(3)) < min(wall%section(ix)%y_safe, wall%section(ix+1)%y_safe)) then
    d_radius = -1
    return
  endif
endif

! triangular calc.
! The wall outward normal is just given by the cross product: (pt1-pt0) x (pt2-pt2)

if (wall%section(ix+1)%basic_shape == 'triangular') then
  if (present(in_antechamber)) in_antechamber = .false.
  if (.not. present(dw_perp)) return
  call sr3d_get_mesh_wall_triangle_pts (wall%section(ix), wall%section(ix+1), p_orb%ix_triangle, pt0, pt1, pt2)
  dp1 = pt1 - pt0
  dp2 = pt2 - pt0
  dw_perp = [dp1(2)*dp2(3) - dp1(3)*dp2(2), dp1(3)*dp2(1) - dp1(1)*dp2(3), dp1(1)*dp2(2) - dp1(2)*dp2(1)]

! Not triangular calc.

else
  ! Get the parameters at the defined cross-sections to either side of the photon position.

  if (p_orb%vec(1) == 0 .and. p_orb%vec(3) == 0) then
    r_photon = 0
    cos_ang = 1
    sin_ang = 0
  else
    r_photon = sqrt(p_orb%vec(1)**2 + p_orb%vec(3)**2)
    cos_ang = p_orb%vec(1) / r_photon
    sin_ang = p_orb%vec(3) / r_photon
  endif

  f = (p_orb%s - wall%section(ix)%s) / (wall%section(ix+1)%s - wall%section(ix)%s)

  ! If f is close to 0 or 1 and dw_perp is not to be calculated then can simplify calc and save time.

  if (.not. present(dw_perp)) then
    if (abs(f) < sr3d_params%significant_length) then
      call sr3d_wall_section_params (wall%section(ix),   cos_ang, sin_ang, radius0, dr0_dtheta, in_ante0)
      d_radius = r_photon - radius0
      if (present(in_antechamber)) in_antechamber = in_ante0
      return

    elseif (abs(f - 1) < sr3d_params%significant_length) then
      call sr3d_wall_section_params (wall%section(ix+1), cos_ang, sin_ang, radius1, dr1_dtheta, in_ante1)
      d_radius = r_photon - radius1
      if (present(in_antechamber)) in_antechamber = in_ante1
      return
    endif
  endif

  ! 

  call sr3d_wall_section_params (wall%section(ix),   cos_ang, sin_ang, radius0, dr0_dtheta, in_ante0)
  call sr3d_wall_section_params (wall%section(ix+1), cos_ang, sin_ang, radius1, dr1_dtheta, in_ante1)

  d_radius = r_photon - ((1 - f) * radius0 + f * radius1)

  if (present(in_antechamber)) in_antechamber = (in_ante0 .and. in_ante1)

  if (present (dw_perp)) then
    dw_perp(1:2) = [cos_ang, sin_ang] - [-sin_ang, cos_ang] * &
                              ((1 - f) * dr0_dtheta + f * dr1_dtheta) / r_photon
    dw_perp(3) = (radius0 - radius1) / (wall%section(ix+1)%s - wall%section(ix)%s)
  endif

endif

! In a bend dw_perp must be corrected since the true longitudinal "length" at the particle
! is, for a horizontal bend, ds * (1 + x/rho) where ds is the length along the reference 
! trajectory, x is the transverse displacement, and rho is the bend radius.

! Also dw_perp needs to be normalized to 1.

if (present(dw_perp)) then
  ix_ele = element_at_s (lat, p_orb%s, .true.)
  ele => lat%ele(ix_ele)
  if (ele%key == sbend$) then
    if (ele%value(ref_tilt_tot$) == 0) then
      disp = p_orb%vec(1) 
    else
      disp = p_orb%vec(1) * cos(ele%value(ref_tilt_tot$)) + p_orb%vec(3) * sin(ele%value(ref_tilt_tot$))
    endif
    dw_perp(3) = dw_perp(3) / (1 + disp * ele%value(g$))
  endif

  dw_perp = dw_perp / sqrt(sum(dw_perp**2))  ! Normalize
endif

end subroutine sr3d_photon_d_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_wall_index (p_orb, wall, ix_wall)
!
! Routine to get the wall index such that 
! For p_orb%vec(6) > 0 (forward motion):
!   wall%section(ix_wall)%s < p_orb%s <= wall%section(ix_wall+1)%s
! For p_orb%vec(6) < 0 (backward motion):
!   wall%section(ix_wall)%s <= p_orb%s < wall%section(ix_wall+1)%s
! Exceptions:
!   If p_orb%s == wall%section(0)%s (= 0)       -> ix_wall = 0
!   If p_orb%s == wall%section(wall%n_section_max)%s -> ix_wall = wall%n_section_max - 1
!
! Input:
!   p_orb  -- sr3d_photon_coord_struct: Photon position.
!   wall   -- sr3d_wall_struct: Wall structure
!
! Output:
!   ix_wall -- Integer: Wall index
!-

subroutine sr3d_get_wall_index (p_orb, wall, ix_wall)

implicit none

type (sr3d_photon_coord_struct) :: p_orb
type (sr3d_wall_struct), target :: wall

integer ix_wall

! 

ix_wall = p_orb%ix_wall
if (p_orb%s < wall%section(ix_wall)%s .or. p_orb%s > wall%section(ix_wall+1)%s) then
  call bracket_index2 (wall%section%s, 0, wall%n_section_max, p_orb%s, p_orb%ix_wall, ix_wall)
  p_orb%ix_wall = ix_wall
  if (ix_wall == wall%n_section_max) ix_wall = wall%n_section_max - 1
endif

! vec(5) at boundary cases

if (p_orb%s == wall%section(ix_wall)%s   .and. p_orb%vec(6) > 0 .and. ix_wall /= 0)               ix_wall = ix_wall - 1
if (p_orb%s == wall%section(ix_wall+1)%s .and. p_orb%vec(6) < 0 .and. ix_wall /= wall%n_section_max-1) ix_wall = ix_wall + 1

end subroutine sr3d_get_wall_index

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_mesh_wall_triangle_pts (pt1, pt2, ix_tri, tri_vert0, tri_vert1, tri_vert2)
!
! Routine to return the three vertex points for a triangular wall surface element between
! two cross-sections.
!
! Input:
!   pt1 -- sr3d_wall_section_struct: A gen_shape or triangular cross-section.
!   pt2 -- sr3d_wall_section_struct: Second cross-section. Should be triangular.
!   ix_tr  -- Integer: Triangle index. Must be between 1 and 2*size(pt1%gen_shape%wall3d_section%v).
!               [Note: size(pt1%gen_shape%wall3d_section%v) = size(pt2%gen_shape%wall3d_section%v)]
!
! Output:
!   tri_vert0(3), tri_vert1(3), tri_vert2(3)
!         -- Real(rp): (x, y, s) vertex points for the triangle.
!             Looking from the outside, the points are in counter-clockwise order.
!             This is important in determining the outward normal vector
!-

subroutine sr3d_get_mesh_wall_triangle_pts (pt1, pt2, ix_tri, tri_vert0, tri_vert1, tri_vert2)

implicit none

type (sr3d_wall_section_struct) pt1, pt2

integer ix_tri
integer ix1, ix2

real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)

! 

ix1 = (ix_tri + 1) / 2
ix2 = ix1 + 1
if (ix2 > size(pt1%gen_shape%wall3d_section%v)) ix2 = 1

if (odd(ix_tri)) then
  tri_vert0 = [pt1%gen_shape%wall3d_section%v(ix1)%x, pt1%gen_shape%wall3d_section%v(ix1)%y, pt1%s]
  tri_vert1 = [pt1%gen_shape%wall3d_section%v(ix2)%x, pt1%gen_shape%wall3d_section%v(ix2)%y, pt1%s]
  tri_vert2 = [pt2%gen_shape%wall3d_section%v(ix1)%x, pt2%gen_shape%wall3d_section%v(ix1)%y, pt2%s]
else
  tri_vert0 = [pt1%gen_shape%wall3d_section%v(ix2)%x, pt1%gen_shape%wall3d_section%v(ix2)%y, pt1%s]
  tri_vert1 = [pt2%gen_shape%wall3d_section%v(ix2)%x, pt2%gen_shape%wall3d_section%v(ix2)%y, pt2%s]
  tri_vert2 = [pt2%gen_shape%wall3d_section%v(ix1)%x, pt2%gen_shape%wall3d_section%v(ix1)%y, pt2%s]
endif

end subroutine sr3d_get_mesh_wall_triangle_pts

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_wall_section_params (wall_section, cos_photon, sin_photon, r_wall, dr_dtheta, in_antechamber)
!
! Routine to compute parameters needed by sr3d_photon_d_radius routine.
!
! Input:
!   wall_section -- sr3d_wall_section_struct: Wall outline at a particular longitudinal location.
!   cos_photon -- Real(rp): Cosine of the photon transverse position.
!   sin_photon -- Real(rp): Sine of the photon transverse position.
!
! Output:
!   r_wall         -- Real(rp): Radius of the wall
!   dr_dtheta      -- Real(rp): Transverse directional derivatives: d(r_wall)/d(theta)
!   in_antechamber -- Logical: Set true of particle is in antechamber
!-

subroutine sr3d_wall_section_params (wall_section, cos_photon, sin_photon, r_wall, dr_dtheta, in_antechamber)

implicit none

type (sr3d_wall_section_struct) wall_section, section
type (wall3d_vertex_struct), pointer :: v(:)

real(rp) dr_dtheta, cos_photon, sin_photon 
real(rp) r_wall

integer ix, ix_vertex, ixv(2)

logical in_antechamber

! Init

in_antechamber = .false.

! general shape

if (wall_section%basic_shape == 'gen_shape') then
  call calc_wall_radius (wall_section%gen_shape%wall3d_section%v, cos_photon, sin_photon, r_wall, dr_dtheta, ix_vertex)

  ixv = wall_section%gen_shape%ix_vertex_ante
  if (ixv(1) > 0) then
    if (ixv(2) > ixv(1)) then
      if (ix_vertex >= ixv(1) .and. ix_vertex < ixv(2)) in_antechamber = .true.
    else
      if (ix_vertex >= ixv(1) .or. ix_vertex < ixv(2)) in_antechamber = .true.
    endif
  endif

  ixv = wall_section%gen_shape%ix_vertex_ante2
  if (ixv(1) > 0) then
    if (ixv(2) > ixv(1)) then
      if (ix_vertex >= ixv(1) .and. ix_vertex < ixv(2)) in_antechamber = .true.
    else
      if (ix_vertex >= ixv(1) .or. ix_vertex < ixv(2)) in_antechamber = .true.
    endif
  endif

  return
endif


! general shape: Should not be here

if (wall_section%basic_shape == 'triangular') then
  call err_exit
endif

! Check for antechamber or beam stop...
! If the line extending from the origin through the photon intersects the
! antechamber or beam stop then pretend the chamber is rectangular with the 
! antechamber or beam stop dimensions.

! Positive x side check.

section = wall_section

if (cos_photon > 0) then

  ! If there is an antechamber...
  if (section%ante_height2_plus > 0) then

    if (abs(sin_photon/cos_photon) < section%ante_height2_plus/section%ante_x0_plus) then  
      section%basic_shape = 'rectangular'
      section%width2 = section%width2_plus
      section%height2 = section%ante_height2_plus
      if (cos_photon >= section%ante_x0_plus) in_antechamber = .true.
    endif

  ! If there is a beam stop...
  elseif (section%width2_plus > 0) then
    if (abs(sin_photon/cos_photon) < section%y0_plus/section%width2_plus) then 
      section%basic_shape = 'rectangular'
      section%width2 = section%width2_plus
    endif

  endif

! Negative x side check

elseif (cos_photon < 0) then

  ! If there is an antechamber...
  if (section%ante_height2_minus > 0) then

    if (abs(sin_photon/cos_photon) < section%ante_height2_minus/section%ante_x0_minus) then  
      section%basic_shape = 'rectangular'
      section%width2 = section%width2_minus
      section%height2 = section%ante_height2_minus
      if (cos_photon >= section%ante_x0_minus) in_antechamber = .true.
    endif

  ! If there is a beam stop...
  elseif (section%width2_minus > 0) then
    if (abs(sin_photon / cos_photon) < section%y0_minus/section%width2_minus) then 
      section%basic_shape = 'rectangular'
      section%width2 = section%width2_minus
    endif

  endif

endif

! Compute parameters

if (section%basic_shape == 'rectangular') then
  if (abs(cos_photon/section%width2) > abs(sin_photon/section%height2)) then
    r_wall = section%width2 / abs(cos_photon)
    dr_dtheta = r_wall * sin_photon / cos_photon
  else
    r_wall = section%height2 / abs(sin_photon)
    dr_dtheta = -r_wall * cos_photon / sin_photon
  endif

elseif (section%basic_shape == 'elliptical') then
  r_wall = 1 / sqrt((cos_photon/section%width2)**2 + (sin_photon/section%height2)**2)
  dr_dtheta = r_wall**3 * cos_photon * sin_photon * (1/section%width2**2 - 1/section%height2**2)
endif

end subroutine sr3d_wall_section_params

end module
