!+
! Module synrad3d_parse_wall
!
! Module for routines needed to parse the vacuum chamber wall file.
!-

module synrad3d_parse_wall

use synrad3d_struct

implicit none

type sr3d_section_struct
  character(60) :: shape_name = '' 
  character(40) :: name = ''
  character(40) :: sub_chamber_name = ''
  character(40) :: surface_name = ''
  real(rp) :: s = 0
  integer :: type = normal$
  integer :: repeat_count = 0
  logical :: surface_is_local = .false.
  type (wall3d_section_struct) :: section 
  type (sr3d_section_struct), pointer :: m_sec => null()    ! Multi-section pointer
  integer :: ix_branch = 0
end type

! multi_section structure

type sr3d_multi_section_struct
  character(40) name
  type (sr3d_section_struct), allocatable :: section(:)
end type

! Vacuum chamber description structure.

type sr3d_wall_struct
  type (sr3d_section_struct), allocatable :: section(:)  ! indexed from 1
  type (sr3d_multi_section_struct), allocatable :: multi_section(:)
  integer n_place
end type

!------------------

type sr3d_section_input
  real(rp) s                      ! Longitudinal position.
  character(40) name              ! Name of setion
  character(60) shape_name       
  integer repeat_count
end type

type surface_input
  character(40) name
  logical surface_is_local
end type

type old_wall_section_input
  real(rp) s                      ! Longitudinal position.
  character(40) name              ! Name of setion
  character(60) shape_name        ! "elliptical", "rectangular", or "gen_shape:xxx"
  real(rp) width2                 ! Half width ignoring antechamber.
  real(rp) height2                ! Half height ignoring antechamber.
  real(rp) width2_plus            ! Distance from pipe center to +x side edge.
  real(rp) ante_height2_plus      ! Antechamber half height on +x side of the wall
  real(rp) width2_minus           ! Distance from pipe center -x side edge.
  real(rp) ante_height2_minus     ! Antechamber half height on -x side of the wall
end type

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_read_wall_file (wall_file, lat, err_flag)
!
! Routine to parse a vacuum chamber wall file.
!
! Input:
!   wall_file   -- character(*): Name of the wall file.
!   lat         -- lat_struct: lattice.
!
! Output:
!   branch%wall3d -- wall3d_struct: Wall structure with computed parameters.
!   err_flag      -- logical, optional: Set true if there is a problem
!-

subroutine sr3d_read_wall_file (wall_file, lat, err_flag)

type this_surface_ptr_struct
  type (photon_reflect_surface_struct), pointer :: ptr  => null()
end type

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch, branch2
type (surface_input) surface
type (this_surface_ptr_struct), allocatable :: surface_ptr(:)
type (wall3d_struct), pointer :: wall3d
type (wall3d_vertex_struct) v(100)
type (wall3d_section_struct), allocatable, target :: shape(:)
type (wall3d_section_struct), pointer :: sec3d, s0, s1
type (sr3d_wall_struct), target :: wall_in
type (sr3d_section_input) section
type (sr3d_multi_section_struct), pointer :: m_sec
type (sr3d_section_struct), allocatable :: temp_section(:)
type (wall3d_section_struct), allocatable :: temp_section3d(:)
type (sr3d_section_struct), pointer :: sec, sec0, sec2
type (sr3d_section_struct) ref_section
type (ele_struct), pointer :: ele, ele2, ele0
type (sr3d_branch_overlap_struct), allocatable :: branch_overlap(:)

real(rp) s_branch, r0(2), max_r1, max_r2, max_r

integer i, j, k, n, ix, iu, im, iw, it, iss, ns, ios, ib, ie, ix0, ix2, i2, nv
integer, allocatable :: n_sub_sec(:), ix_sort(:)
integer n_shape, n_repeat, n_surface, last_type, n_overlap
integer m_max, n_add, n_sub, ix_ele0, ix_ele1, ix_bend, ix_patch, ix_b, ix_p
integer ix_slow, ix_fast, chamber_end_geometry, n_sec
integer ix_ele1_start, ix_ele1_end, ix_ele2_start, ix_ele2_end, ixe

logical, optional :: err_flag
logical err, absolute_vertices

character(*) wall_file
character(40) name, slow, fast, branch_name, subchamber_name
character(40), allocatable :: sub_name(:)
character(200) reflectivity_file, file
character(*), parameter :: r_name = 'sr3d_read_wall_file'

namelist / place / section, surface
namelist / shape_def / name, v, r0, absolute_vertices
namelist / surface_def / reflectivity_file
namelist / slow_fast / slow, fast
namelist / subchamber_branch / subchamber_name, branch_name 

! Open file

if (present(err_flag)) err_flag = .true.
iu = lunget()
call fullfilename (wall_file, file)
open (iu, file = file, status = 'old')

! Old file format?

if (sr3d_old_wall_file_format(iu)) then
  call out_io (s_fatal$, r_name, 'WALL FILE IS USING AN OLD FORMAT!! ' // trim(wall_file), &
               'PLEASE SEE THE SYNRAD3D MANUAL FOR INSTRUCTIONS ON HOW TO CONVERT TO THE NEW FORMAT.', &
               'STOPPING NOW.')
  stop
endif

! Read in the surface info

rewind(iu)
n_surface = 0
do
  read (iu, nml = surface_def, iostat = ios)
  if (ios > 0) then ! error
    call out_io (s_fatal$, r_name, 'ERROR READING SURFACE_DEF NAMELIST IN WALL FILE: ' // trim(wall_file))
    rewind (iu)
    do
      read (iu, nml = surface_def) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  n_surface = n_surface + 1
enddo

if (allocated(sr3d_com%surface)) deallocate(sr3d_com%surface)
allocate (sr3d_com%surface(n_surface+3))

call photon_reflection_std_surface_init (sr3d_com%surface(1))

sr3d_com%surface(2)%reflectivity_file = '<none>'
sr3d_com%surface(2)%name = 'ABSORBER'
sr3d_com%surface(2)%description = 'Perfect Absorber'

sr3d_com%surface(3)%reflectivity_file = '<none>'
sr3d_com%surface(3)%name = 'PHANTOM'
sr3d_com%surface(3)%description = 'Virtual Wall'

rewind(iu)
do i = 4, n_surface+3
  read (iu, nml = surface_def, iostat = ios)
  call read_surface_reflection_file (reflectivity_file, sr3d_com%surface(i))
enddo

! Read multi_section

call sr3d_read_wall_multi_section (iu, wall_in)

! Get wall info
! First count the cross-section number

rewind(iu)
wall_in%n_place = 0
do
  read (iu, nml = place, iostat = ios)
  if (ios > 0) then ! error
    rewind (iu)
    do
      read (iu, nml = place) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  wall_in%n_place = wall_in%n_place + 1
enddo


call out_io (s_info$, r_name, 'Number of wall cross-sections read: \i0\ ', wall_in%n_place)
if (wall_in%n_place < 1) then
  call out_io (s_fatal$, r_name, 'NO WALL SPECIFIED IN WALL FILE: ' // trim(wall_file), &
                                 'WILL STOP HERE.')
  call err_exit
endif

if (allocated(wall_in%section)) deallocate(wall_in%section)
allocate (wall_in%section(1:wall_in%n_place))

! Now transfer info from the file to the wall_in%section array

rewind (iu)
do i = 1, wall_in%n_place
  section%shape_name = ''
  section%repeat_count = -1
  section%name = ''

  surface%name = ''
  surface%surface_is_local = .false.

  read (iu, nml = place)

  sec => wall_in%section(i)
  sec%shape_name        = section%shape_name
  sec%repeat_count      = section%repeat_count
  sec%name              = section%name
  sec%s                 = section%s

  ix = index(section%shape_name, ':')
  if (ix /= 0) then
    sec%sub_chamber_name  = section%shape_name(1:ix-1)
    sec%shape_name = section%shape_name(ix+1:)
  endif

  ix = index(sec%shape_name, '@')
  if (ix /= 0) then
    select case (sec%shape_name(ix+1:))
    case ('')
      sec%type = normal$
    case ('START')
      sec%type = wall_start$
    case ('END')
      sec%type = wall_end$
    case default
      call out_io (s_fatal$, r_name, 'BAD SECTION STOP POINT SPECIFICATION: ' // trim(section%shape_name))
    end select
    sec%shape_name = sec%shape_name(1:ix-1)
  endif

  sec%surface_name = surface%name
  sec%surface_is_local = surface%surface_is_local
enddo

! Get the info as to what branch the subchambers are associated with

rewind(iu)
do
  
  read (iu, nml = subchamber_branch, iostat = ios)
  if (ios > 0) then ! If error
    call out_io (s_fatal$, r_name, 'ERROR READING SUBCHAMBER_BRANCH NAMELIST IN WALL FILE:' // trim(wall_file))
    rewind (iu)
    do
      read (iu, nml = subchamber_branch) ! generate error message
    enddo
  endif
  if (ios < 0) exit  ! End of file

  branch => pointer_to_branch (branch_name, lat)
  if (.not. associated(branch)) then
    call out_io (s_fatal$, r_name, 'ERROR READING SUBCHAMBER_BRANCH NAMELIST', &
                    'CANNOT FIND LATTICE BRANCH ASSOCIATED WITH BRANCH NAME: ' // branch_name)
    stop
  endif

  do i = 1, wall_in%n_place
    sec => wall_in%section(i)
    if (sec%sub_chamber_name /= subchamber_name) cycle
    sec%ix_branch = branch%ix_branch
  enddo

enddo

! Get the shape info

n_shape = 0
rewind(iu)
do
  read (iu, nml = shape_def, iostat = ios)
  if (ios > 0) then ! If error
    call out_io (s_fatal$, r_name, 'ERROR READING SHAPE_DEF NAMELIST IN WALL FILE:' //  trim(wall_file))
    rewind (iu)
    do
      read (iu, nml = shape_def) ! Generate error message
    enddo
  endif
  if (ios < 0) exit  ! End of file
  n_shape = n_shape + 1
enddo

if (allocated(shape)) deallocate(shape)
allocate (shape(n_shape))

rewind(iu)
do i = 1, n_shape
  v = wall3d_vertex_struct()
  name = ''
  absolute_vertices = .false.
  r0 = 0

  read (iu, nml = shape_def, iostat = ios)

  if (name == '') then
    call out_io (s_fatal$, r_name, 'SHAPE_DEF DOES NOT HAVE A NAME! SHAPE NUMBER: \i0\ ', i)
    call err_exit
  endif

  do j = 1, i-1
    if (shape(j)%name == name) then
      call out_io (s_fatal$, r_name, 'TWO SHAPE_DEFS HAVE THE SAME NAME: ' // trim(name))
      call err_exit
    endif
  enddo
  shape(i)%name = name

  ! Count number of vertices and calc angles.

  sec3d => shape(i)
  do nv = 1, size(v)
    if (v(nv)%x == 0 .and. v(nv)%y == 0 .and. v(nv)%radius_x == 0) exit
  enddo

  if (any(v(nv:)%x /= 0) .or. any(v(nv:)%y /= 0) .or. &
      any(v(nv:)%radius_x /= 0) .or. any(v(nv:)%radius_y /= 0)) then
    call out_io (s_fatal$, r_name, 'MALFORMED SHAPE:' // name)
    call err_exit
  endif

  allocate(sec3d%v(nv-1))
  sec3d%v = v(1:nv-1)
  sec3d%r0 = r0
  sec3d%n_vertex_input = nv-1

  if (absolute_vertices) then
    sec3d%vertices_state = absolute$
  else
    sec3d%vertices_state = relative$
  endif

  call wall3d_section_initializer (sec3d, err)
  if (err) then
    call out_io (s_fatal$, r_name, 'ERROR AT SHAPE: ' // trim(name))
    call err_exit
  endif

enddo

! Expand multi_sections

i = 0
outer: do 
  i = i + 1
  if (i > wall_in%n_place) exit
  ref_section = wall_in%section(i)
  if (all(ref_section%shape_name /= wall_in%multi_section%name)) cycle
  n_repeat = ref_section%repeat_count

  if (n_repeat < 0) then
    call out_io (s_fatal$, r_name, 'ERROR: MULTI_SECTION DOES NOT HAVE THE REPEAT COUNT SET.')
    call err_exit
  endif

  do j = 1, size(wall_in%multi_section)
    m_sec => wall_in%multi_section(j)
    if (ref_section%shape_name /= m_sec%name) cycle
    m_max = ubound(m_sec%section, 1)

    if (m_sec%section(m_max)%name == 'closed') then
      n_add = n_repeat * (m_max - 1) + 1
    elseif (m_sec%section(m_max)%name == 'open') then
      n_add = n_repeat * (m_max - 1)
    else
      call out_io (s_fatal$, r_name, 'ERROR: LAST SECTION IN MULTI_SECTION IS NOT "closed" NOR "open"')
      call err_exit
    endif

    call move_alloc (wall_in%section, temp_section)
    n = wall_in%n_place
    allocate (wall_in%section(1:n+n_add-1))
    wall_in%section(1:i-1) = temp_section(1:i-1)
    wall_in%section(i+n_add:n+n_add-1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall_in%n_place = n + n_add - 1

    do k = 1, n_repeat
      do im = 1, m_max - 1
        sec2   => wall_in%section(i+(k-1)*(m_max-1)+(im-1))
        sec2   = m_sec%section(im)
        sec2%s = ref_section%s + (k-1) * m_sec%section(m_max)%s + m_sec%section(im)%s
        sec2%sub_chamber_name = ref_section%sub_chamber_name
        sec2%surface_name     = ref_section%surface_name
        sec2%m_sec => m_sec%section(im)
      enddo
    enddo

    if (m_sec%section(m_max)%name == 'closed') then
      sec2   => wall_in%section(i+n_repeat*(m_max-1))
      sec2   = m_sec%section(1)
      sec2%s = ref_section%s + n_repeat * m_sec%section(m_max)%s
      sec2%sub_chamber_name = ref_section%sub_chamber_name
      sec2%surface_name     = ref_section%surface_name
      sec2%m_sec => m_sec%section(1)
    endif

    cycle outer

  enddo

  call out_io (s_fatal$, r_name, 'CANNOT FIND MATCHING MULTI_SECTION NAME: ' // trim(ref_section%shape_name))
  call err_exit

enddo outer

! point to shapes

section_loop: do i = 1, wall_in%n_place
  sec => wall_in%section(i)
  do j = 1, size(shape)
    if (sec%shape_name /= shape(j)%name) cycle
    sec%section = shape(j)
    cycle section_loop
  enddo

  call out_io (s_fatal$, r_name, 'CANNOT FIND MATCHING SHAPE FOR: ' // trim(sec%shape_name))
  call err_exit

enddo section_loop

! Checks

do i = 1, wall_in%n_place
  sec => wall_in%section(i)

  ! Shape checks

  if (sec%shape_name == 'shape') then
    if (sec%section%name == '') then
      call out_io (s_fatal$, r_name, 'BAD SHAPE ASSOCIATION')
      call err_exit
    endif
    cycle
  endif

enddo

! Count number of subchambers

allocate (sub_name(10), ix_sort(10))
n_sub = 0

do i = 1, wall_in%n_place
  if (n_sub == size(sub_name)) then
    call re_allocate(sub_name, 2*n_sub)
    call re_allocate(ix_sort, 2*n_sub)
  endif
  sec => wall_in%section(i)
  write(name, '(2a, i0)') trim(sec%sub_chamber_name), ':', sec%ix_branch
  call find_index(name, sub_name, ix_sort, n_sub, ix, add_to_list = .true.)
enddo

! Associate surface

do i = 1, wall_in%n_place
  sec => wall_in%section(i)
  call sr3d_associate_surface (sec%section%surface, sec%surface_name, sr3d_com%surface)
enddo

allocate (surface_ptr(n_sub))
do i = 1, n_sub
  surface_ptr(i)%ptr => sr3d_com%surface(1)  ! Default surface
enddo

do i = wall_in%n_place, 1, -1
  sec => wall_in%section(i)
  write(name, '(2a, i0)') trim(sec%sub_chamber_name), ':', sec%ix_branch
  call find_index(name, sub_name, ix_sort, n_sub, ix, add_to_list = .false.)

  if (associated(sec%section%surface)) then
    if (.not. sec%surface_is_local) surface_ptr(ix)%ptr => sec%section%surface
  else
    sec%section%surface => surface_ptr(ix)%ptr
  endif
enddo

!--------------------------------------------------------------------------------------------
! Transfer info to branch%wall3d

! Count how many subchambers there are and allocate branch%wall3d(:) array.

if (allocated(sr3d_com%branch)) deallocate (sr3d_com%branch)
allocate(sr3d_com%branch(0:ubound(lat%branch,1)))
allocate (n_sub_sec(n_sub))

do ib = 0, ubound(lat%branch, 1)

  n_sub = 0
  n_sub_sec = 0
  branch => lat%branch(ib)
  s_branch = branch%ele(branch%n_ele_track)%s

  chamber_end_geometry = sr3d_params%chamber_end_geometry
  if (chamber_end_geometry < 0 .or. ib > 0) chamber_end_geometry = branch%param%geometry

  do i = 1, wall_in%n_place
    sec => wall_in%section(i)
    if (sec%ix_branch /= ib) cycle
    name = sec%sub_chamber_name
    call find_index(name, sub_name, ix_sort, n_sub, ix, add_to_list = .true.)
    n_sub_sec(ix) = n_sub_sec(ix) + 1
  enddo

  ! Allocate space for walls

  if (associated(branch%wall3d)) deallocate (branch%wall3d)
  if (n_sub == 0) cycle

  if (allocated(sr3d_com%branch(ib)%fast)) deallocate (sr3d_com%branch(ib)%fast)
  allocate (branch%wall3d(n_sub), sr3d_com%branch(ib)%fast(n_sub))

  do iw = 1, n_sub
    wall3d => branch%wall3d(iw)
    wall3d%name = sub_name(iw)
    wall3d%ix_wall3d = iw
    allocate (wall3d%section(n_sub_sec(iw)))
    allocate(sr3d_com%branch(ib)%fast(iw)%ix_wall3d(0))
  enddo

  ! Sort sections

  n_sub_sec = 0

  do i = 1, wall_in%n_place
    sec => wall_in%section(i)
    if (sec%ix_branch /= ib) cycle
    call find_index(sec%sub_chamber_name, sub_name, ix_sort, n_sub, iss)
    wall3d => branch%wall3d(iss)
    n_sub_sec(iss) = n_sub_sec(iss) + 1
    ns = n_sub_sec(iss)  
    sec3d           => wall3d%section(ns)
    sec3d           = sec%section
    sec3d%s         = sec%s
    sec3d%type      = sec%type
    sec3d%ix_branch = branch%ix_branch
    sec3d%name      = sec%name
    if (sec3d%name == '') sec3d%name = sec%shape_name

    ! Check s ordering...
    ! The problem with two sections at the same s-value is that the sr3d_photon_hit_spot_calc routine
    ! cannot handle a photon hitting a wall perpendicular to the longitudinal azis. To avoid this 
    ! situation, move a section by a tiny amount longitudinally if two sections have the same s.

    if (ns > 1) then
      if (sec3d%s == wall3d%section(ns-1)%s) sec3d%s = sec3d%s + 1000*sr3d_params%significant_length
      if (sec3d%s < wall3d%section(ns-1)%s) then
        call out_io (s_fatal$, r_name, &
                  'SECTION(i)%S: \f0.4\ ', &
                  'IS LESS THAN SECTION(i-1)%S: \f0.4\ ', &
                  'FOR I = \i0\ ', &
                  r_array = [sec3d%s, wall3d%section(ns-1)%s], i_array = [i])
        call err_exit
      endif
    endif

    ! This check is currently not needed since the above code shifts sections longitudanally.

    if (ns > 2) then
      if (sec3d%s < wall3d%section(ns-2)%s) then
        call out_io (s_fatal$, r_name, &
                  'THREE WALL SECTIONS WITH THE SAME LONGITUDINAL S VALUE NOT ALLOWED.', &
                  'AT S = \f0.4\ ', &
                  'FOR I = \i0\ ', &
                  r_array = [sec3d%s], i_array = [i])
        call err_exit
      endif
    endif

  enddo

  ! Start and End sections must come in non-overlapping pairs

  do iw = 1, n_sub
    wall3d => branch%wall3d(iw)
    last_type = normal$
    n_sec = 0
    do iss = 1, size(wall3d%section)
      it = wall3d%section(iss)%type
      if (it == normal$) cycle

      if (it == wall_start$ .and. last_type == wall_start$) then
        call out_io (s_fatal$, r_name, 'TWO STARTING WALL SECTIONS WITHOUT AN INTROVENING END SECTION IN SUBCHAMBER: ' // wall3d%name)
        call err_exit
      elseif (it == wall_end$ .and. last_type == wall_end$) then
        call out_io (s_fatal$, r_name, 'TWO ENDING WALL SECTIONS WITHOUT AN INTROVENING START SECTION IN SUBCHAMBER: ' // wall3d%name)
        call err_exit
      endif

      last_type = it
      n_sec = n_sec + 1
    enddo

    If (mod(n_sec, 2) == 1) then
      call out_io (s_fatal$, r_name, 'NUMBER OF START SECTIONS NOT EQUAL TO NUMBER OF END SECTIONS IN SUBCHAMBER: ' // wall3d%name)
      call err_exit
    endif
  enddo

  ! For an open lattice geometry if the first section of a subchamber is 
  ! not marked as a starting or ending section it must have s = 0

  if (chamber_end_geometry == open$) then
    do i = 1, size(branch%wall3d)
      wall3d => branch%wall3d(i)
      if (wall3d%section(1)%type /= normal$) cycle  
      if (wall3d%section(1)%s /= 0) then
        call out_io (s_info$, r_name, &
              'subchamber named "' // trim(wall3d%name) // '" begins at: \f12.4\ ', &
              'And not at lattice beginning.', &
              r_array = [wall3d%section(1)%s])
      endif
    enddo
  endif

  ! Some wall checks

  do i = 1, size(branch%wall3d)
    wall3d => branch%wall3d(i)
    n_sec = ubound(wall3d%section, 1)

    ! Last section gets shifted to lattice end if it is beyound.
    if (wall3d%section(n_sec)%s > s_branch) wall3d%section(n_sec)%s = s_branch

    ! For an open lattice geometry, subchambers must either stop or have cross section at end of lattice.
    if (chamber_end_geometry == open$) then
      if (n_sec == 1) then
        call out_io (s_fatal$, r_name, 'subchamber named "' // trim(wall3d%name) // '" has only one section.')
        call err_exit
      endif

      if (wall3d%section(1)%type /= wall_start$ .and. wall3d%section(1)%s /= 0) then
        call out_io (s_info$, r_name, &
              'subchamber named "' // trim(wall3d%name) // '" has first cross-section at s =: \f12.4\ ', &
              'And not at s = 0.', &
              r_array = [wall3d%section(1)%s])
        call err_exit
      endif

      if (wall3d%section(n_sec)%type /= wall_end$ .and. abs(wall3d%section(n_sec)%s - s_branch) > 0.01) then
        call out_io (s_info$, r_name, &
              'subchamber named "' // trim(wall3d%name) // '" ends at: \f12.4\ ', &
              'And not at lattice end of: \f12.4\ ', &
              '[But last point is always adjusted to have s = s_branch]', &
              r_array = [wall3d%section(n_sec)%s, s_branch])
      endif

      if (wall3d%section(n_sec)%type /= wall_end$) wall3d%section(n_sec)%s = s_branch

    endif
  enddo

  ! Regions between wall sections are not allowed to contain both bends and patch elements.
  ! If this is the case then add additional sections to avoid this situation

  do iw = 1, size(branch%wall3d)
    wall3d => branch%wall3d(iw)

    i = 0
    do
      i = i + 1
      n_sec = size(wall3d%section)

      if (i > n_sec) exit
      if (i == n_sec .and. wall3d%section(n_sec)%type == wall_end$) exit

      ix_ele0 = element_at_s(lat, wall3d%section(i)%s,   .true., branch%ix_branch)
      i2 = mod(i, n_sec) + 1
      ix_ele1 = element_at_s(lat, wall3d%section(i2)%s, .false., branch%ix_branch)

      ix_bend = -1    ! Index of last bend before patch or first bend after patch
      ix_patch = -1   ! Index of last patch before bend or first patch after bend
      ele => branch%ele(ix_ele0)

      ! Must be able to handle case where the wall has only one wall section

      if (ix_ele0 == ix_ele1 .and. wall3d%section(i)%s < wall3d%section(i2)%s) cycle

      do 
        if (ele%key == sbend$) ix_bend = ele%ix_ele
        if (ele%key == patch$) ix_patch = ele%ix_ele
        if (ix_bend /= -1 .and. ix_patch /= -1) exit

        if (ele%ix_ele == ix_ele1 .and. ix_ele1 /= ix_ele0) exit

        if (ele%key == patch$ .and. ele%orientation == -1) then
          call out_io (s_fatal$, r_name, 'PATCH ELEMENTS WITH ORIENTATION = -1 NOT YET IMPLEMENTED!')
          call err_exit
        endif

        ele => pointer_to_next_ele(ele)
        if (ele%ix_ele == ix_ele1 .and. ix_ele1 == ix_ele0) exit
      enddo

      if (ix_bend == -1 .or. ix_patch == -1) cycle

      ! Need to add an additional section

      call move_alloc (wall3d%section, temp_section3d)
      allocate (wall3d%section(1:n_sec+1))

      ! Wall wrapping around case.
      if (ix_bend < ix_ele0 .and. ix_patch < ix_ele0) then
        wall3d%section(2:n_sec+1) = temp_section3d(1:n_sec)
        sec3d => wall3d%section(1)
        ix0 = n_sec + 1
        ix2 = 2
        i = 0  ! Reset
      ! Normal case.
      else
        wall3d%section(1:i) = temp_section3d(1:i)
        if (i /= n_sec) wall3d%section(i+2:n_sec+1) = temp_section3d(i+1:n_sec)
        sec3d => wall3d%section(i+1)
        ix0 = i
        ix2 = mod(i+1, n_sec+1) + 1
      endif

      deallocate (temp_section3d)

      ! Remember: Patch can have negative length

      ix_b = ix_bend
      if (ix_b < ix_ele0) ix_b = ix_b + branch%n_ele_track

      ix_p = ix_patch
      if (ix_p < ix_ele0) ix_p = ix_p + branch%n_ele_track

      if (ix_b < ix_p) then  ! Add section at end of bend, before patch
        sec3d = wall3d%section(ix0)
        sec3d%s = min(branch%ele(ix_bend)%s, branch%ele(ix_patch)%s_start)
      else  ! Add section at beginning of bend, after patch
        sec3d = wall3d%section(ix2)
        sec3d%s = max(branch%ele(ix_bend)%s_start, branch%ele(ix_patch)%s)
      endif

      if (sec3d%name(1:6) /= 'ADDED:') then
        j = len(sec3d%name)
        sec3d%name = 'ADDED:' // sec3d%name(1:j-6)
      endif

      call out_io (s_info$, r_name, 'Extra section added to separate bend and patch at s = \f10.2\ ', r_array = [sec3d%s])

    enddo
  enddo

  ! Associate lattice element with section
  ! Error if the section is inside a patch element. 
  ! OK if on the edge but associated element (%ix_ele) may not point to the patch since there may
  ! be a ambiguity of the coordinate system associated with the section.

  do iw = 1, size(branch%wall3d)
    wall3d => branch%wall3d(iw)

    do iss = 1, size(wall3d%section)
      sec3d => wall3d%section(iss)
      if (iss == 1) then
        sec3d%ix_ele    = element_at_s(lat, sec3d%s, .false., branch%ix_branch)
      else
        sec3d%ix_ele    = element_at_s(lat, sec3d%s, .true., branch%ix_branch)
      endif

      ix = sec3d%ix_ele
      if (branch%ele(ix)%key == patch$) then
        if (sec3d%s == branch%ele(ix-1)%s) then
          sec3d%ix_ele = ix - 1
        elseif (sec3d%s == branch%ele(ix)%s) then
          sec3d%ix_ele = ix + 1
        else
          call out_io (s_fatal$, r_name, 'WALL CROSS-SECTION AT S = \f10.3\ ', &
                                         'IS WITHIN A PATCH ELEMENT. THIS IS NOT ALLOWED.', r_array = [sec3d%s])
          call err_exit
        endif
      endif

    enddo
  enddo 

enddo ! Branch loop

! Associate slow/fast walls

rewind(iu)
do
  read (iu, nml = slow_fast, iostat = ios)
  if (ios > 0) then ! error
    call out_io (s_fatal$, r_name, 'ERROR READING SLOW_FAST NAMELIST IN WALL FILE: ' // trim(wall_file))
    rewind (iu)
    do
      read (iu, nml = slow_fast) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  
  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    ix_slow = 0
    ix_fast = 0

    do iw = 1, size(branch%wall3d)
      wall3d => branch%wall3d(iw)
      if (wall3d%name == slow) ix_slow = iw
      if (wall3d%name == fast) ix_fast = iw
    enddo

    if (ix_slow /= 0 .and. ix_fast /= 0) exit
  enddo

  if (ix_slow == 0 .or. ix_fast == 0) then
    call out_io (s_fatal$, r_name, 'SLOW AND/OR FAST SUBCHAMBER NOT FOUND: ' // trim(slow), '  ' // trim(fast))
    call err_exit
  endif

  n = size(sr3d_com%branch(ib)%fast(ix_slow)%ix_wall3d) + 1
  call re_allocate (sr3d_com%branch(ib)%fast(ix_slow)%ix_wall3d, n)
  sr3d_com%branch(ib)%fast(ix_slow)%ix_wall3d(n) = ix_fast
enddo

close (iu)

!--------------------------------------------------------------------------------------------
! Find branch overlap regions

n_overlap = 0
if (allocated(sr3d_com%branch_overlap)) deallocate (sr3d_com%branch_overlap)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (.not. associated(branch%wall3d)) cycle

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= fork$ .and. ele%key /= photon_fork$) cycle

    branch2 => lat%branch(nint(ele%value(ix_to_branch$)))
    if (.not. associated(branch2%wall3d)) cycle
    ele2 => branch2%ele(nint(ele%value(ix_to_element$)))

    ! In branch:  Overlap region are elements with indexes in the range [ix_ele1_start, ix_ele1_end]
    !             max_r1 is the maximum chamber radius in this range
    ! In branch2: Overlap region are elements with indexes in the range [ix_ele2_start, ix_ele2_end]
    !             max_r2 is the maximum chamber radius in this range

    call overlap_region_rough_guess (ele, branch, ix_ele1_start, ix_ele1_end, max_r1) 
    call overlap_region_rough_guess (ele2, branch2, ix_ele2_start, ix_ele2_end, max_r2)
    max_r = 2 * (max_r1 + max_r2) ! Factor of two is a safety factor.

    n_overlap = n_overlap + 1
    call move_alloc (sr3d_com%branch_overlap, branch_overlap)
    allocate (sr3d_com%branch_overlap(n_overlap))
    if (n_overlap > 1) sr3d_com%branch_overlap(1:n_overlap-1) = branch_overlap(1:n_overlap-1)

    sr3d_com%branch_overlap(n_overlap)%ix_branch1 = ib
    sr3d_com%branch_overlap(n_overlap)%ix_branch2 = branch2%ix_branch

    ! Refine ix_ele1_start

    ixe = (ix_ele2_start + ix_ele2_end) / 2
    ele0 => branch%ele(ix_ele1_start)
    do
      if (overlaps(0.0_rp, ele0, branch2%ele(ixe)) .or. overlaps(ele0%value(l$), ele0, branch2%ele(ixe))) exit
      if (ele0%ix_ele == ix_ele1_end) then
        call out_io (s_fatal$, r_name, 'CANNOT FIND OVERLAP REGION FOR THE WALLS OF TWO LATTICE BRANCHES.')
        stop
      endif
      ele0 => pointer_to_next_ele (ele0, +1)
    enddo
    sr3d_com%branch_overlap(n_overlap)%ix_ele1_start = ele0%ix_ele

    ! Refine ix_ele1_end

    ele0 => branch%ele(ix_ele1_end)
    do
      if (overlaps(0.0_rp, ele0, branch2%ele(ixe)) .or. overlaps(ele0%value(l$), ele0, branch2%ele(ixe))) exit
      if (ele0%ix_ele == ix_ele1_start) then
        call out_io (s_fatal$, r_name, 'CANNOT FIND OVERLAP REGION FOR THE WALLS OF TWO LATTICE BRANCHES.')
        stop
      endif
      ele0 => pointer_to_next_ele (ele0, -1)
    enddo
    sr3d_com%branch_overlap(n_overlap)%ix_ele1_end = ele0%ix_ele


    ! Refine ix_ele2_start

    ixe = (ix_ele1_start + ix_ele1_end) / 2
    ele0 => branch2%ele(ix_ele2_start)
    do
      if (overlaps(0.0_rp, ele0, branch2%ele(ixe)) .or. overlaps(ele0%value(l$), ele0, branch2%ele(ixe))) exit
      if (ele0%ix_ele == ix_ele2_end) then
        call out_io (s_fatal$, r_name, 'CANNOT FIND OVERLAP REGION FOR THE WALLS OF TWO LATTICE BRANCHES.')
        stop
      endif
      ele0 => pointer_to_next_ele (ele0, +1)
    enddo
    sr3d_com%branch_overlap(n_overlap)%ix_ele2_start = ele0%ix_ele

    ! Refine ix_ele2_end

    ele0 => branch2%ele(ix_ele2_end)
    do
      if (overlaps(0.0_rp, ele0, branch2%ele(ixe)) .or. overlaps(ele0%value(l$), ele0, branch2%ele(ixe))) exit
      if (ele0%ix_ele == ix_ele2_start) then
        call out_io (s_fatal$, r_name, 'CANNOT FIND OVERLAP REGION FOR THE WALLS OF TWO LATTICE BRANCHES.')
        stop
      endif
      ele0 => pointer_to_next_ele (ele0, -1)
    enddo
    sr3d_com%branch_overlap(n_overlap)%ix_ele2_end = ele0%ix_ele

  enddo
enddo

! Cleanup

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (.not. associated(branch%wall3d)) cycle
  call mark_patch_regions (branch)
enddo

deallocate(shape)
if (present(err_flag)) err_flag = .false.

!-------------------------------------------------------------------------
contains

! Init to 10 meters from branch point

subroutine overlap_region_rough_guess (ele, branch, ix_ele_start, ix_ele_end, max_r)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele2
type (branch_struct) branch
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: sec

integer ix_ele_start, ix_ele_end
integer iw, is, j
real(rp) max_r, len_region, ll, angle, r_wall, dr_dtheta

!

len_region = min(10.0_rp, branch%param%total_length/2)

ele2 => ele
ll = 0

do
  if (chamber_end_geometry == open$ .and. ele2%ix_ele == 1) exit
  ele2 => pointer_to_next_ele(ele2, -1)
  ll = ll + ele2%value(l$)
  if (ll > len_region) exit
enddo

ix_ele_start = ele2%ix_ele 

!

ele2 => ele
ll = 0

do
  if (chamber_end_geometry == open$ .and. ele2%ix_ele == branch%n_ele_track) exit
  ele2 => pointer_to_next_ele(ele2, 1)
  ll = ll + ele2%value(l$)
  if (ll > len_region) exit
enddo

ix_ele_end = ele2%ix_ele 

!

max_r = 0

do iw = 1, size(branch%wall3d)
  wall3d => branch%wall3d(iw)
  do is = 1, size(wall3d%section)
    sec => wall3d%section(is)
    if (ix_ele_start <= ix_ele_end) then
      if (sec%ix_ele < ix_ele_start .or. sec%ix_ele > ix_ele_end) cycle
    else
      if (sec%ix_ele > ix_ele_start .or. sec%ix_ele < ix_ele_end) cycle
    endif

    do j = 0, 7
      angle = j * twopi / 8
      call calc_wall_radius(sec%v, cos(angle), sin(angle), r_wall, dr_dtheta)
      r_wall = sqrt((sec%r0(1) + r_wall*cos(angle))**2 + (sec%r0(2) + r_wall*sin(angle))**2)
      max_r = max(max_r, r_wall)
    enddo
  enddo
enddo

end subroutine overlap_region_rough_guess

!-------------------------------------------------------------------------
! contains

function overlaps (s, ele, ele2) result (has_overlap)

type (branch_struct) branch2
type (ele_struct) ele, ele2
type (ele_struct), pointer :: ele2_p
type (floor_position_struct) floor

real(rp) s
integer ix_ele_start, ix_ele2_start, ix_ele2_end, status
logical has_overlap

!

has_overlap = .false.

floor%r = [0.0_rp, 0.0_rp, s]
floor = coords_local_curvilinear_to_floor(floor, ele, calculate_angles = .false.)
floor = coords_floor_to_curvilinear (floor, ele2, ele2_p, status)

if (status == outside$ .or. status == patch_problem$) return
if (floor%r(1)**2 + floor%r(2)**2 > max_r) return

has_overlap = .true.

end function overlaps

end subroutine sr3d_read_wall_file 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function sr3d_old_wall_file_format (iunit) result (is_old)
!
! Routine to see if a wall file has the old format.
!
! Input:

function sr3d_old_wall_file_format (iunit) result (is_old)

type (old_wall_section_input) section
type (surface_input) surface

integer iunit, ios
logical is_old

namelist / section_def / section, surface

!

read (iunit, nml = section_def, iostat = ios)
rewind (iunit)

is_old = (ios == 0)

end function sr3d_old_wall_file_format

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
character(*), parameter :: r_name = 'sr3d_associate_surface'
integer i

!

nullify(surface_ptr)
if (surface_name == '') return

do i = 1, size(surfaces)
  if (surfaces(i)%name == surface_name) then
    surface_ptr => surfaces(i)
    return
  endif
enddo

call out_io (s_fatal$, r_name, 'NO SURFACE CORRESPONDING TO: ' // trim(surface_name))
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
type (sr3d_section_input), target :: section(1:100)
type (sr3d_section_input), pointer :: sec_in
type (sr3d_section_struct), pointer :: sec

integer i, j, iu, ix, n_multi, n_section, ios

character(40) surface, name
character(*), parameter :: r_name = 'sr3d_read_wall_multi_section'

namelist / multi_place / name, section

! Count multi_sections

rewind(iu)
n_multi = 0
do
  read (iu, nml = multi_place, iostat = ios)
  if (ios > 0) then ! error
    rewind (iu)
    do
      read (iu, nml = multi_place) ! will bomb program with error message
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
  section%shape_name = ''
  section%name = ''
  name = ''
  read (iu, nml = multi_place)

  if (name == '') then
    call out_io (s_fatal$, r_name, 'MULTI_PLACE DOES NOT HAVE A NAME! MULTI_PLACE NUMBER:', i)
    call err_exit
  endif

  do j = 1, i-1
    if (wall%multi_section(j)%name == name) then
      call out_io (s_fatal$, r_name, 'TWO MULTI_PLACES HAVE THE SAME NAME: ' // trim(name))
      call err_exit
    endif
  enddo

  wall%multi_section(i)%name = name

  n_section = count(section%shape_name /= '')
  if (any(section(n_section+1:)%shape_name /= '')) then
    call out_io (s_fatal$, r_name, 'CONFUSED MULTI_PLACE: ' // trim(name))
    call err_exit
  endif
  allocate (wall%multi_section(i)%section(1:n_section))

  do j = 1, n_section
    sec => wall%multi_section(i)%section(j)
    sec_in => section(j)
    sec%name             = sec_in%name
    sec%shape_name       = sec_in%shape_name
    sec%s                = sec_in%s
  enddo
enddo

end subroutine sr3d_read_wall_multi_section


end module
