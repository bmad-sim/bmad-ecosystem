!+
! Module synrad3d_parse_wall
!
! Module for routines needed to parse the vacuum chamber wall file.
!-

module synrad3d_parse_wall

use synrad3d_struct

implicit none

type sr3d_section_struct
  character(60) :: section_id = '' 
  character(40) :: name = ''
  character(40) :: sub_chamber_name = ''
  character(40) :: surface_name = ''
  real(rp) :: s = 0
  integer :: type = normal$
  integer :: repeat_count = 0
  logical :: surface_is_local = .false.
  type (wall3d_section_struct) :: section 
  type (sr3d_section_struct), pointer :: m_sec => null()    ! Multi-section pointer
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
  integer n_section_max
end type

!------------------

type sr3d_section_input
  real(rp) s                      ! Longitudinal position.
  character(40) name              ! Name of setion
  character(60) section_id       
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
! Subroutine sr3d_read_wall_file (wall_file, branch, err_flag)
!
! Routine to parse a vacuum chamber wall file.
!
! Input:
!   wall_file   -- character(*): Name of the wall file.
!   branch      -- branch_struct: Lattice branch associated with the wall..
!
! Output:
!   branch%wall3d -- wall3d_struct: Wall structure with computed parameters.
!   err_flag      -- logical, optional: Set true if there is a problem
!-

subroutine sr3d_read_wall_file (wall_file, branch, err_flag)

type (branch_struct), target :: branch
type (surface_input) surface
type (photon_reflect_surface_struct), pointer :: surface_ptr  => null()
type (wall3d_struct), pointer :: wall3d
type (wall3d_vertex_struct) v(100)
type (wall3d_section_struct), allocatable, target :: shape(:)
type (wall3d_section_struct), pointer :: sec3d
type (sr3d_wall_struct), target :: wall_in
type (sr3d_section_input) section
type (sr3d_multi_section_struct), pointer :: m_sec
type (sr3d_section_struct), allocatable :: temp_section(:)
type (sr3d_section_struct), pointer :: sec, sec0, sec2
type (sr3d_section_struct) ref_section

real(rp) s_lat, r0(2)

integer i, j, k, n, ix, iu, im, iw, it, iss, ns, ios
integer, allocatable :: n_sub_sec(:), ix_sort(:)
integer n_place, n_shape, n_repeat, n_surface, last_type
integer m_max, n_add, n_sub, ix_ele0, ix_ele1, ix_bend, ix_patch

logical, optional :: err_flag
logical err, absolute_vertices

character(*) wall_file
character(40) name
character(40), allocatable :: sub_name(:)
character(200) reflectivity_file, file
character(28), parameter :: r_name = 'sr3d_read_wall_file'

namelist / place / section, surface
namelist / shape_def / name, v, r0, absolute_vertices
namelist / surface_def / name, reflectivity_file

! Open file

iu = lunget()
call fullfilename (wall_file, file)
open (iu, file = file, status = 'old')

! Old file format?

if (sr3d_old_wall_file_format(iu)) then
  print '(a)', 'WALL FILE IS USING AN OLD FORMAT!!'
  print '(a)', 'PLEASE SEE THE SYNRAD3D MANUAL FOR INSTRUCTIONS ON HOW TO CONVERT TO THE NEW FORMAT.'
  print '(a)', 'STOPPING NOW.'
  stop
endif

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

if (associated(branch%lat%surface)) deallocate(branch%lat%surface)
allocate (branch%lat%surface(n_surface+2))

branch%lat%surface(1)%reflectivity_file = ''
branch%lat%surface(1)%descrip = 'default'
call photon_reflection_std_surface_init (branch%lat%surface(1))

branch%lat%surface(2)%reflectivity_file = ''
branch%lat%surface(2)%descrip = 'ABSORBER'

rewind(iu)
do i = 3, n_surface+2
  read (iu, nml = surface_def, iostat = ios)
  call read_surface_reflection_file (reflectivity_file, branch%lat%surface(i))
  branch%lat%surface(i)%descrip = name
  branch%lat%surface(i)%reflectivity_file = reflectivity_file
enddo

! Read multi_section

call sr3d_read_wall_multi_section (iu, wall_in)

! Get wall info
! First count the cross-section number

rewind(iu)
n_place = 0
do
  read (iu, nml = place, iostat = ios)
  if (ios > 0) then ! error
    rewind (iu)
    do
      read (iu, nml = place) ! will bomb program with error message
    enddo  
  endif
  if (ios < 0) exit   ! End of file reached
  n_place = n_place + 1
enddo


print *, 'number of wall cross-sections read:', n_place
if (n_place < 1) then
  print *, 'NO WALL SPECIFIED. WILL STOP HERE.'
  call err_exit
endif

if (allocated(wall_in%section)) deallocate(wall_in%section)
allocate (wall_in%section(1:n_place))
wall_in%n_section_max = n_place


! Now transfer info from the file to the wall_in%section array

rewind (iu)
do i = 1, n_place
  section%section_id = ''
  section%repeat_count = -1
  section%name = ''

  surface%name = ''
  surface%surface_is_local = .false.

  read (iu, nml = place)

  sec => wall_in%section(i)
  sec%section_id        = section%section_id
  sec%repeat_count      = section%repeat_count
  sec%name              = section%name
  sec%s                 = section%s

  ix = index(section%section_id, ':')
  if (ix /= 0) then
    sec%sub_chamber_name  = section%section_id(1:ix-1)
    sec%section_id = section%section_id(ix+1:)
  endif

  ix = index(sec%section_id, '@')
  if (ix /= 0) then
    select case (sec%section_id(ix+1:))
    case ('')
      sec%type = normal$
    case ('START')
      sec%type = wall_start$
    case ('END')
      sec%type = wall_end$
    case default
      print *, 'BAD SECTION STOP POINT SPECIFICATION: ', trim(section%section_id)
    end select
    sec%section_id = sec%section_id(1:ix-1)
  endif

  sec%surface_name = surface%name
  sec%surface_is_local = surface%surface_is_local
enddo

! Get the shape info

n_shape = 0
rewind(iu)
do
  read (iu, nml = shape_def, iostat = ios)
  if (ios > 0) then ! If error
    print *, 'ERROR READING SHAPE_DEF NAMELIST.'
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

  read (iu, nml = shape_def, iostat = ios)

  if (name == '') then
    print *, 'SHAPE_DEF DOES NOT HAVE A NAME! SHAPE NUMBER:', i
    call err_exit
  endif

  do j = 1, i-1
    if (shape(j)%name == name) then
      print *, 'TWO SHAPE_DEFS HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo
  shape(i)%name = name

  ! Count number of vertices and calc angles.

  sec3d => shape(i)
  do n = 1, size(v)
    if (v(n)%x == 0 .and. v(n)%y == 0 .and. v(n)%radius_x == 0) exit
  enddo

  if (any(v(n:)%x /= 0) .or. any(v(n:)%y /= 0) .or. &
      any(v(n:)%radius_x /= 0) .or. any(v(n:)%radius_y /= 0)) then
    print *, 'MALFORMED SHAPE:', name
    call err_exit
  endif

  allocate(sec3d%v(n-1))
  sec3d%v = v(1:n-1)
  sec3d%r0 = r0
  sec3d%n_vertex_input = n-1    

  sec3d%absolute_vertices_input = absolute_vertices

  call wall3d_section_initializer (sec3d, err)
  if (err) then
    print *, 'ERROR AT SHAPE: ', trim(name)
    call err_exit
  endif

enddo

close (iu)

! Expand multi_sections

i = 0
outer: do 
  i = i + 1
  if (i > wall_in%n_section_max) exit
  ref_section = wall_in%section(i)
  if (ref_section%section_id /= 'multi_section') cycle
  n_repeat = ref_section%repeat_count

  if (n_repeat < 0) then
    print *, 'ERROR: MULTI_SECTION DOES NOT HAVE THE REPEAT COUNT SET.'
    call err_exit
  endif

  do j = 1, size(wall_in%multi_section)
    m_sec => wall_in%multi_section(j)
    if (ref_section%section_id /= m_sec%name) cycle
    m_max = ubound(m_sec%section, 1)

    if (m_sec%section(m_max)%name == 'closed') then
      n_add = n_repeat * (m_max - 1) + 1
    elseif (m_sec%section(m_max)%name == 'open') then
      n_add = n_repeat * (m_max - 1)
    else
      print *, 'ERROR: LAST SECTION IN MULTI_SECTION IS NOT "closed" NOR "open"'
      call err_exit
    endif

    call move_alloc (wall_in%section, temp_section)
    n = wall_in%n_section_max
    allocate (wall_in%section(1:n+n_add-1))
    wall_in%section(1:i-1) = temp_section(1:i-1)
    wall_in%section(i+n_add:n+n_add-1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall_in%n_section_max = n + n_add - 1

    do k = 1, n_repeat
      do im = 1, m_max - 1
        sec2 => wall_in%section(i+(k-1)*(m_max-1)+(im-1))
        sec2 = m_sec%section(im)
        sec2%s = ref_section%s + (k-1) * m_sec%section(m_max)%s + m_sec%section(im)%s
        sec2%m_sec => m_sec%section(im)
      enddo
    enddo

    if (m_sec%section(m_max)%name == 'closed') then
      sec2 => wall_in%section(i+n_repeat*(m_max-1))
      sec2 = m_sec%section(1)
      sec2%s = ref_section%s + n_repeat * m_sec%section(m_max)%s
      sec2%m_sec => m_sec%section(1)
    endif

    cycle outer

  enddo

  print *, 'CANNOT FIND MATCHING MULTI_SECTION NAME: ', trim(ref_section%section_id)
  call err_exit

enddo outer

! point to shapes

section_loop: do i = 1, wall_in%n_section_max
  sec => wall_in%section(i)
  do j = 1, size(shape)
    if (sec%section_id /= shape(j)%name) cycle
    sec%section = shape(j)
    cycle section_loop
  enddo

  print *, 'CANNOT FIND MATCHING SHAPE FOR: ', trim(sec%section_id)
  call err_exit

enddo section_loop

! Checks

do i = 1, wall_in%n_section_max
  sec => wall_in%section(i)

  ! Shape checks

  if (sec%section_id == 'shape') then
    if (sec%section%name == '') then
      call out_io (s_fatal$, r_name, 'BAD SHAPE ASSOCIATION')
      call err_exit
    endif
    cycle
  endif

enddo

! If circular lattice then begin and end shapes must match

if (branch%param%geometry == closed$) then
  sec0 => wall_in%section(1)
  sec  => wall_in%section(wall_in%n_section_max)
  if (sec0%section_id /= sec%section_id) then
    call out_io (s_fatal$, r_name, 'FOR A "CLOSED" LATTICE THE LAST WALL CROSS-SECTION MUST BE THE SAME AS THE FIRST.')
    call err_exit
  endif
endif

! Regions between wall sections are not allowed to contain both bends and patch elements.
! If this is the case then add additional sections to avoid this situation

i = 0
do
  i = i + 1
  if (i == wall_in%n_section_max) exit

  ix_ele0 = element_at_s(branch%lat, wall_in%section(i)%s, .true., branch%ix_branch)
  ix_ele1 = element_at_s(branch%lat, wall_in%section(i+1)%s, .false., branch%ix_branch)

  ix_bend = -1    ! Index of last bend before patch or first bend after patch
  ix_patch = -1   ! Index of last patch before bend or first patch after bend
  do j = ix_ele0, ix_ele1
    if (branch%ele(j)%key == sbend$ .and. (ix_bend == -1 .or. ix_patch == -1)) ix_bend = j
    if (branch%ele(j)%key == patch$ .and. (ix_bend == -1 .or. ix_patch == -1)) ix_patch = j
  enddo

  if (branch%ele(j)%key == patch$ .and. branch%ele(j)%orientation == -1) then
    call out_io (s_fatal$, r_name, 'PATCH ELEMENTS WITH ORIENTATION = -1 NOT YET IMPLEMENTED!')
    call err_exit
  endif

  if (ix_bend == -1 .or. ix_patch == -1) cycle

  ! Need to add an additional section

  if (size(wall_in%section) == wall_in%n_section_max) then
    call move_alloc (wall_in%section, temp_section)
    n = wall_in%n_section_max
    allocate (wall_in%section(1:n+10))
    wall_in%section(1:i) = temp_section(1:i)
    wall_in%section(i+2:n+1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall_in%n_section_max = n + 1
  endif

  sec => wall_in%section(i+1)

  if (ix_bend < ix_patch) then  ! Add section at end of bend, before patch
    sec = wall_in%section(i)
    ! Remember: Patch can have negative length
    sec%s = min(branch%ele(ix_bend)%s, branch%ele(ix_patch)%s - branch%ele(ix_patch)%value(l$))
  else  ! Add section at beginning of bend, after patch
    sec = wall_in%section(i+2)
    ! Remember: Patch can have negative length
    sec%s = max(branch%ele(ix_bend)%s - branch%ele(ix_bend)%value(l$), branch%ele(ix_patch)%s)
  endif

  if (sec%name(1:6) /= 'ADDED:') then
    n = len(sec%name)
    sec%name = 'ADDED:' // sec%name(1:n-6)
  endif

  call out_io (s_info$, r_name, 'Extra section added to separate bend and patch at s = \f10.2\ ', &
                                r_array = [sec%s])

enddo

! Associate surface

do i = 1, wall_in%n_section_max
  sec => wall_in%section(i)
  call sr3d_associate_surface (sec%section%surface, sec%surface_name, branch%lat%surface)
enddo

surface_ptr => branch%lat%surface(1)  ! Default surface

do i = wall_in%n_section_max, 1, -1
  sec => wall_in%section(i)

  if (associated(sec%section%surface)) then
    if (.not. sec%surface_is_local) surface_ptr => sec%section%surface
  else
    sec%section%surface => surface_ptr
  endif

enddo

!---------------------------------------------------------------
! Transfer info to branch%wall3d

! Count how many sub-chambers there are and allocate branch%wall3d(:) array.

allocate (sub_name(10), ix_sort(10), n_sub_sec(10))

n_sub = 0
n_sub_sec = 0

do i = 1, n_place
  if (n_sub == size(sub_name)) then
    call re_allocate(sub_name, n_sub+10)
    call re_allocate(ix_sort, n_sub+10, init_val = 0)
    call re_allocate(n_sub_sec, n_sub+10, init_val = 0)
  endif

  name = wall_in%section(i)%sub_chamber_name
  call find_indexx(name, sub_name, ix_sort, n_sub, ix, add_to_list = .true.)
  n_sub_sec(ix) = n_sub_sec(ix) + 1
enddo

if (associated(branch%wall3d)) deallocate (branch%wall3d, sr3d_com%fast)
allocate (branch%wall3d(n_sub), sr3d_com%fast(n_sub))
do iw = 1, n_sub
  wall3d => branch%wall3d(iw)
  wall3d%name = sub_name(iw)
  wall3d%ix_wall3d = iw
  allocate (wall3d%section(n_sub_sec(iw)))
  sr3d_com%fast(iw)%wall3d => null()
enddo

outer_loop: do iw = 1, n_sub
  wall3d => branch%wall3d(iw)
  if (wall3d%name /= 'FAST' .and. wall3d%name(1:5) /= 'FAST_') cycle
  name = wall3d%name(6:)
  if (wall3d%name == 'FAST') name = wall3d%name(5:)
  do i = 1, n_sub
    if (branch%wall3d(i)%name /= wall3d%name(6:)) cycle
    sr3d_com%fast(i)%wall3d => wall3d
    cycle outer_loop
  enddo
  call out_io (s_fatal$, r_name, 'CANNOT FIND CORRESPONDING SUB-CHAMBER TO FAST SUB-CHAMBER: ' // wall3d%name)
  call err_exit
enddo outer_loop

n_sub_sec = 0

do i = 1, n_place
  call find_indexx(wall_in%section(i)%sub_chamber_name, sub_name, ix_sort, n_sub, iss)
  wall3d => branch%wall3d(iss)
  n_sub_sec(iss) = n_sub_sec(iss) + 1
  ns = n_sub_sec(iss)  
  sec3d           => wall3d%section(ns)
  sec3d           = wall_in%section(i)%section
  sec3d%s         = wall_in%section(i)%s
  sec3d%name      = wall_in%section(i)%name
  sec3d%type      = wall_in%section(i)%type
  sec3d%ix_branch = branch%ix_branch
  sec3d%ix_ele    = element_at_s(branch%lat, sec3d%s, .true., branch%ix_branch)

  ! Check s ordering

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



  ! Error if the section is inside a patch element. 
  ! OK if on the edge but associated element (%ix_ele) may not point to the patch since there may
  ! be a ambiguity of the coordinate system associated with the section.

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

! Begin and end sections must come in non-overlapping pairs

do iw = 1, n_sub
  wall3d => branch%wall3d(iw)
  last_type = normal$
  n = 0
  do iss = 1, size(wall3d%section)
    it = wall3d%section(iss)%type
    if (it == normal$) cycle

    if (it == wall_start$ .and. last_type == wall_start$) then
      call out_io (s_fatal$, r_name, 'TWO STARTING WALL SECTIONS WITHOUT AN INTROVENING END SECTION IN SUB-CHAMBER: ' // wall3d%name)
      call err_exit
    elseif (it == wall_end$ .and. last_type == wall_end$) then
      call out_io (s_fatal$, r_name, 'TWO ENDING WALL SECTIONS WITHOUT AN INTROVENING START SECTION IN SUB-CHAMBER: ' // wall3d%name)
      call err_exit
    endif

    last_type = it
    n = n + 1
  enddo

  if (mod(n, 2) == 1) then
    call out_io (s_fatal$, r_name, 'NUMBER OF START SECTIONS NOT EQUAL TO NUMBER OF END SECTIONS IN SUB-CHAMBER: ' // wall3d%name)
    call err_exit
  endif
enddo

! First section of a sub-chamber must have s = 0
! Last section of a sub-chamber has s adjusted to match the lattice length.
! Except: if first/last section represents a sub-section beginning or ending.

s_lat = branch%ele(branch%n_ele_track)%s

do i = 1, size(branch%wall3d)
  wall3d => branch%wall3d(i)
  if (wall3d%section(1)%type /= normal$) cycle  
  if (wall3d%section(1)%s /= 0) then
    call out_io (s_info$, r_name, &
          'Sub-chamber named "' // trim(wall3d%name) // '" begins at: \f12.4\ ', &
          'And not at lattice beginning.', &
          r_array = [wall3d%section(1)%s])
  endif
enddo

do i = 1, size(branch%wall3d)
  wall3d => branch%wall3d(i)
  n = ubound(wall3d%section, 1)
  if (wall3d%section(n)%type /= normal$) cycle  
  if (abs(wall3d%section(n)%s - s_lat) > 0.01) then
    call out_io (s_info$, r_name, &
          'Sub-chamber named "' // trim(wall3d%name) // '" ends at: \f12.4\ ', &
          'And not at lattice end of: \f12.4\ ', &
          '[But last point is always adjusted to have s = s_lat]', &
          r_array = [wall3d%section(n)%s, s_lat])
  endif
  wall3d%section(n)%s = s_lat
enddo

! Cleanup

call mark_patch_regions (branch)

deallocate(shape)

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
type (sr3d_section_input), target :: section(1:100)
type (sr3d_section_input), pointer :: sec_in
type (sr3d_section_struct), pointer :: sec

integer i, j, iu, ix, n_multi, n_section, ios

character(40) surface, name

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
  section%section_id = ''
  section%name = ''
  name = ''
  read (iu, nml = multi_place)

  if (name == '') then
    print *, 'MULTI_PLACE DOES NOT HAVE A NAME! MULTI_PLACE NUMBER:', i
    call err_exit
  endif

  do j = 1, i-1
    if (wall%multi_section(j)%name == name) then
      print *, 'TWO MULTI_PLACES HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo

  wall%multi_section(i)%name = name

  n_section = count(section%section_id /= '')
  if (any(section(n_section+1:)%section_id /= '')) then
    print *, 'CONFUSED MULTI_PLACE: ', trim(name)
    call err_exit
  endif
  allocate (wall%multi_section(i)%section(1:n_section))

  do j = 1, n_section
    sec => wall%multi_section(i)%section(j)
    sec_in => section(j)
    sec%section%name     = sec_in%name
    sec%section_id       = sec_in%section_id
    sec%section%s        = sec_in%s
  enddo
enddo

end subroutine sr3d_read_wall_multi_section


end module
