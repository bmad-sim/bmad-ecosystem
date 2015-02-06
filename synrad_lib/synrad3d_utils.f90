module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod
use capillary_mod
use track1_mod

implicit none

type sr3d_wall_section_input
  real(rp) s                      ! Longitudinal position.
  character(40) name              ! Name of setion
  character(60) basic_shape       ! "elliptical", "rectangular", or "gen_shape:xxx"
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
! Subroutine sr3d_read_wall_file (wall_file, branch, err_flag)
!
! Routine to check the vacuum chamber wall for problematic values.
! Also compute some wall parameters
!
! Input:
!   wall_file   -- character(*): Name of the wall file.
!   branch      -- branch_struct: Lattice branch.
!
! Output:
!   branch%wall3d -- wall3d_struct: Wall structure with computed parameters.
!   err_flag      -- logical, optional: Set true if there is a problem
!-

subroutine sr3d_read_wall_file (wall_file, branch, err_flag)

implicit none

type (branch_struct), target :: branch
type (sr3d_wall_struct), target :: wall
type (sr3d_wall_section_struct), pointer :: sec, sec0, sec1, sec2
type (sr3d_wall_section_struct), allocatable :: temp_section(:)
type (sr3d_wall_section_struct) ref_section
type (sr3d_wall_section_input) section
type (wall3d_vertex_struct) v(100)
type (wall3d_section_struct), pointer :: sec3d
type (sr3d_multi_section_struct), pointer :: m_sec
type (surface_input) surface
type (photon_reflect_surface_struct), pointer :: surface_ptr  => null()
type (wall3d_struct), pointer :: wall3d
type (sr3d_gen_shape_struct), allocatable, target :: gen_shape(:)

real(rp) ix_vertex_ante(2), ix_vertex_ante2(2), s_lat
real(rp) rad, radius(4), area, max_area, cos_a, sin_a, angle, dr_dtheta

integer i, j, k, im, n, ig, ix, iu, iv, n_wall_section_max, ios, n_shape, n_surface, n_repeat
integer m_max, n_add, ix1, ix2, ix_ele0, ix_ele1, ix_bend, ix_patch

character(28), parameter :: r_name = 'sr3d_read_wall_file'
character(40) name
character(200) reflectivity_file, file
character(*) wall_file

logical, optional :: err_flag
logical err, multi_local

namelist / section_def / section, surface
namelist / gen_shape_def / name, v, ix_vertex_ante, ix_vertex_ante2
namelist / surface_def / name, reflectivity_file

! Open file

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

if (associated(branch%lat%surface)) deallocate(branch%lat%surface)
allocate (branch%lat%surface(n_surface+1))
branch%lat%surface(1)%reflectivity_file = ''
branch%lat%surface(1)%descrip = 'default'
call photon_reflection_std_surface_init(branch%lat%surface(1))

rewind(iu)
do i = 2, n_surface+1
  read (iu, nml = surface_def, iostat = ios)
  call read_surface_reflection_file (reflectivity_file, branch%lat%surface(i))
  branch%lat%surface(i)%descrip = name
  branch%lat%surface(i)%reflectivity_file = reflectivity_file
enddo

! Read multi_section

call sr3d_read_wall_multi_section (iu, wall)

! Get wall info
! First count the cross-section number

rewind(iu)
n_wall_section_max = 0
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

print *, 'number of wall cross-sections read:', n_wall_section_max
if (n_wall_section_max < 1) then
  print *, 'NO WALL SPECIFIED. WILL STOP HERE.'
  call err_exit
endif

if (allocated(wall%section)) deallocate(wall%section)
allocate (wall%section(1:n_wall_section_max))
wall%n_section_max = n_wall_section_max

! Now transfer info from the file to the wall%section array

rewind (iu)
do i = 1, n_wall_section_max
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

  sec%name                = section%name
  sec%basic_shape         = section%basic_shape
  sec%shape_name          = section%basic_shape
  sec%s                   = section%s
  sec%width2              = section%width2
  sec%height2             = section%height2
  sec%width2_plus         = section%width2_plus
  sec%ante_height2_plus   = section%ante_height2_plus
  sec%width2_minus        = section%width2_minus
  sec%ante_height2_minus  = section%ante_height2_minus

  ix = index(section%basic_shape, ':')
  if (ix /= 0) then
    sec%basic_shape = section%basic_shape(1:ix-1)
    sec%shape_name  = section%basic_shape(ix+1:)
  endif

  sec%surface_name = surface%name
  sec%is_local = surface%is_local

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

if (allocated(gen_shape)) deallocate(gen_shape)
allocate (gen_shape(n_shape))

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
    if (gen_shape(j)%name == name) then
      print *, 'TWO GEN_SHAPE_DEFS HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo
  gen_shape(i)%name = name

  ! Count number of vertices and calc angles.

  sec3d => gen_shape(i)%wall3d_section
  do n = 1, size(v)
    if (v(n)%x == 0 .and. v(n)%y == 0 .and. v(n)%radius_x == 0) exit
  enddo

  if (any(v(n:)%x /= 0) .or. any(v(n:)%y /= 0) .or. &
      any(v(n:)%radius_x /= 0) .or. any(v(n:)%radius_y /= 0)) then
    print *, 'MALFORMED GEN_SHAPE:', name
    call err_exit
  endif

  allocate(sec3d%v(n-1))
  sec3d%v = v(1:n-1)
  sec3d%n_vertex_input = n-1    

  call wall3d_section_initializer (sec3d, err)
  if (err) then
    print *, 'ERROR AT GEN_SHAPE: ', trim(name)
    call err_exit
  endif

  gen_shape(i)%ix_vertex_ante = ix_vertex_ante
  if (ix_vertex_ante(1) > 0 .or. ix_vertex_ante(2) > 0) then
    if (ix_vertex_ante(1) < 1 .or. ix_vertex_ante(1) > size(sec3d%v) .or. &
        ix_vertex_ante(2) < 1 .or. ix_vertex_ante(2) > size(sec3d%v)) then
      print *, 'ERROR IN IX_VERTEX_ANTE:', ix_vertex_ante
      print *, '      FOR GEN_SHAPE: ', trim(name)
      call err_exit
    endif
  endif

  gen_shape(i)%ix_vertex_ante2 = ix_vertex_ante2
  if (ix_vertex_ante2(1) > 0 .or. ix_vertex_ante2(2) > 0) then
    if (ix_vertex_ante2(1) < 1 .or. ix_vertex_ante2(1) > size(sec3d%v) .or. &
        ix_vertex_ante2(2) < 1 .or. ix_vertex_ante2(2) > size(sec3d%v)) then
      print *, 'ERROR IN IX_VERTEX_ANTE2:', ix_vertex_ante2
      print *, '      FOR GEN_SHAPE: ', trim(name)
      call err_exit
    endif
  endif

enddo

close (iu)

! Expand multi_sections

i = 0
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
    allocate (wall%section(1:n+n_add-1))
    wall%section(1:i-1) = temp_section(1:i-1)
    wall%section(i+n_add:n+n_add-1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall%n_section_max = n + n_add - 1

    do k = 1, n_repeat
      do im = 1, m_max
        sec2 => wall%section(i+(k-1)*m_max+im)
        sec2 = m_sec%section(im)
        sec2%s = ref_section%s + (k-1) * m_sec%section(m_max)%s + m_sec%section(im)%s
        sec2%m_sec => m_sec%section(im)
      enddo
    enddo
    if (m_sec%section(m_max)%name == 'closed') then
      sec2 => wall%section(i+n_repeat*m_max)
      sec2 = m_sec%section(1)
      sec2%s = ref_section%s + n_repeat * m_sec%section(m_max)%s
      sec2%m_sec => m_sec%section(1)
    endif

    cycle outer

  enddo

  print *, 'CANNOT FIND MATCHING MULTI_SECTION NAME: ', trim(ref_section%shape_name)
  call err_exit

enddo outer

! point to gen_shapes

section_loop: do i = 1, wall%n_section_max
  sec => wall%section(i)
  if (sec%basic_shape /= 'gen_shape') cycle
  do j = 1, size(gen_shape)
    if (sec%shape_name /= gen_shape(j)%name) cycle
    sec%gen_shape = gen_shape(j)
    cycle section_loop
  enddo

  print *, 'CANNOT FIND MATCHING SHAPE FOR: ', trim(sec%shape_name)
  call err_exit

enddo section_loop

! Last section has s adjusted to match the lattice length.

s_lat = branch%ele(branch%n_ele_track)%s

if (abs(wall%section(wall%n_section_max)%s - s_lat) > 0.01) then
  call out_io (s_info$, r_name, &
        'Wall ends at: \f12.4\ ', &
        'And not at lattice end of: \f12.4\ ', &
        '[But last point is always adjusted to have s = s_lat]', &
        r_array = [wall%section(wall%n_section_max)%s, s_lat])
endif

wall%section(wall%n_section_max)%s = s_lat

! Checks

do i = 1, wall%n_section_max
  sec => wall%section(i)

  ! Check s ordering

  if (i > 1) then
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
                                    'gen_shape     '])) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%SECTION(i)%BASIC_SHAPE: ' // sec%basic_shape, &
              '    FOR I = \i0\ ', i_array = [i])
    call err_exit
  endif

  ! Gen_shape checks

  if (sec%basic_shape == 'gen_shape') then
    if (sec%gen_shape%name == '') then
      call out_io (s_fatal$, r_name, 'BAD SHAPE ASSOCIATION')
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

if (branch%param%geometry == closed$) then
  sec0 => wall%section(1)
  sec  => wall%section(wall%n_section_max)
  if (sec0%basic_shape /= sec%basic_shape .or. sec0%width2 /= sec%width2 .or. sec0%height2 /= sec%height2 .or. &
        sec0%ante_height2_plus /= sec%ante_height2_plus .or. sec0%width2_plus /= sec%width2_plus .or. &
        sec0%ante_height2_minus /= sec%ante_height2_minus .or. sec0%width2_minus /= sec%width2_minus) then
      call out_io (s_fatal$, r_name, &
              'FOR A "CLOSED" LATTICE THE LAST WALL CROSS-SECTION MUST BE THE SAME AS THE FIRST.')
      call err_exit
  endif
endif

! Regions between wall sections are not allowed to contain both bends and patch elements.
! If this is the case then add additional sections to avoid this situation

i = 0
do
  i = i + 1
  if (i == wall%n_section_max) exit

  ix_ele0 = element_at_s(branch%lat, wall%section(i)%s, .true., branch%ix_branch)
  ix_ele1 = element_at_s(branch%lat, wall%section(i+1)%s, .false., branch%ix_branch)

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

  if (size(wall%section) == wall%n_section_max) then
    call move_alloc (wall%section, temp_section)
    n = wall%n_section_max
    allocate (wall%section(1:n+10))
    wall%section(1:i) = temp_section(1:i)
    wall%section(i+2:n+1) = temp_section(i+1:n)
    deallocate (temp_section)
    wall%n_section_max = n + 1
  endif

  sec => wall%section(i+1)

  if (ix_bend < ix_patch) then  ! Add section at end of bend, before patch
    sec = wall%section(i)
    ! Remember: Patch can have negative length
    sec%s = min(branch%ele(ix_bend)%s, branch%ele(ix_patch)%s - branch%ele(ix_patch)%value(l$))
  else  ! Add section at beginning of bend, after patch
    sec = wall%section(i+2)
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

! Computations

do i = 1, wall%n_section_max
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

! Create a gen_shape for each section that does not have one.
! Also set antechamber segments.

ig = n_shape

do i = 1, wall%n_section_max
  sec => wall%section(i)

  ! If there is an associated shape then mark where antichamber is.

  if (sec%gen_shape%name /= '') then
    sec3d => sec%gen_shape%wall3d_section

    ix1 = sec%gen_shape%ix_vertex_ante(1)
    ix2 = sec%gen_shape%ix_vertex_ante(2)
    if (ix2 > ix1) then
      do iv = ix1+1, ix2
        sec3d%v(iv)%type = antechamber$
      enddo
    elseif (ix2 < ix1) then
      do iv = ix1+1, size(sec3d%v)
        sec3d%v(iv)%type = antechamber$
      enddo
      do iv = 1, ix2
        sec3d%v(iv)%type = antechamber$
      enddo
    endif

    ix1 = sec%gen_shape%ix_vertex_ante2(1)
    ix2 = sec%gen_shape%ix_vertex_ante2(2)
    if (ix2 > ix1) then
      do iv = ix1+1, ix2
        sec3d%v(iv)%type = antechamber$
      enddo
    elseif (ix2 < ix1) then
      do iv = ix1+1, size(sec3d%v)
        sec3d%v(iv)%type = antechamber$
      enddo
      do iv = 1, ix2
        sec3d%v(iv)%type = antechamber$
      enddo
    endif

    cycle
  endif

  ! Here if there is no associated gen_shape
  ! All equivalent multi-section sections point to the same gen_shape.

  if (associated(sec%m_sec)) then  
    if (sec%m_sec%gen_shape%name /= '') then
      sec%gen_shape = sec%m_sec%gen_shape
      cycle
    endif
  endif

  ! Simple rectangle or ellipse

  sec3d => sec%gen_shape%wall3d_section

  if (sec%width2_plus <= 0 .and. sec%width2_minus <= 0) then
    allocate (sec3d%v(1))
    sec3d%n_vertex_input = 1
    if (sec%basic_shape == 'elliptical') then
      sec3d%v(1)%radius_x = sec%width2
      sec3d%v(1)%radius_y = sec%height2
    else
      sec3d%v(1)%x = sec%width2
      sec3d%v(1)%y = sec%height2
    endif  

  ! Has antichamber(s) or beam stop(s)

  else
    do n = 1, 10
      v(n) = wall3d_vertex_struct()
    enddo

    if (sec%width2_plus <= 0) then   ! No antechamber or beam stop
      iv = 1
      v(1)%x = sec%width2
      if (sec%basic_shape == 'rectangular') v(1)%y = sec%height2
    elseif (sec%width2_plus < sec%width2) then  ! Beam stop
      iv = 1
      v(1)%x = sec%width2_plus
      v(1)%y = sec%y0_plus
    else                                        ! Antechamber
      v(1)%x    = sec%width2_plus
      v(1)%y    = sec%ante_height2_plus
      v(1)%type = antechamber$
      v(2)%x    = sec%ante_x0_plus
      v(2)%y    = sec%ante_height2_plus
      v(2)%type = antechamber$
      iv = 2
      if (sec%basic_shape == 'rectangular') then
        v(3)%x    = sec%width2
        v(3)%y    = sec%height2
        iv = 3
      endif
    endif

    if (sec%basic_shape == 'elliptical') then
      v(iv+1)%radius_x = sec%width2
      v(iv+1)%radius_y = sec%height2
    endif

    if (sec%width2_minus <= 0) then   ! No antechamber or beam stop
      iv = iv + 1
      v(iv)%x = -sec%width2
      if (sec%basic_shape == 'rectangular') v(iv)%y =  sec%height2
    elseif (sec%width2_minus < sec%width2) then  ! Beam stop
      iv = iv + 1
      v(iv)%x = -sec%width2_minus
      v(iv)%y =  sec%y0_minus
    else                                        ! Antechamber
      if (sec%basic_shape == 'rectangular') then
        iv = iv + 1
        v(iv)%x    = -sec%width2
        v(iv)%y    = sec%height2
      endif
      v(iv+1)%x    = -sec%ante_x0_minus
      v(iv+1)%y    =  sec%ante_height2_minus
      v(iv+2)%x    = -sec%width2_minus
      v(iv+2)%y    =  sec%ante_height2_minus
      v(iv+2)%type = antechamber$
      v(iv+3)%x    = -sec%width2_minus
      v(iv+3)%y    =  0
      v(iv+3)%type = antechamber$
      iv = iv + 3
    endif

    allocate (sec3d%v(iv))
    sec3d%v = v(1:iv)
    sec3d%n_vertex_input = iv

    if (sec%width2_plus > sec%width2) sec%gen_shape%ix_vertex_ante = [size(sec3d%v)-1, 2]
    if (sec%width2_minus > sec%width2) sec%gen_shape%ix_vertex_ante2 = [iv-2, iv+2]

  endif

  call wall3d_section_initializer (sec3d, err)
  if (err) call err_exit

enddo

! Associate surface

do i = 1, wall%n_section_max
  sec => wall%section(i)
  call sr3d_associate_surface (sec%gen_shape%wall3d_section%surface, sec%surface_name, branch%lat%surface)
enddo

surface_ptr => branch%lat%surface(1)  ! Default surface

do i = wall%n_section_max, 1, -1
  sec => wall%section(i)

  if (associated(sec%gen_shape%wall3d_section%surface)) then
    if (.not. sec%is_local) surface_ptr => sec%gen_shape%wall3d_section%surface
  else
    sec%gen_shape%wall3d_section%surface => surface_ptr
  endif

enddo

! Transfer info to branch%wall3d

if (associated(branch%wall3d)) deallocate (branch%wall3d)
allocate(branch%wall3d)
wall3d => branch%wall3d

allocate (wall3d%section(wall%n_section_max))
do i = 1, wall%n_section_max
  sec3d           => wall3d%section(i)
  sec3d           = wall%section(i)%gen_shape%wall3d_section
  sec3d%s         = wall%section(i)%s
  sec3d%name      = wall%section(i)%name
  sec3d%ix_branch = branch%ix_branch
  sec3d%ix_ele    = element_at_s(branch%lat, sec3d%s, .true., branch%ix_branch)

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

deallocate(gen_shape)

call mark_patch_regions (branch)

! Calculate largest "safe" box for each section.
! implemented in this instance. 

wall3d%section%x_safe = -1
wall3d%section%y_safe = -1

do i = 1, wall%n_section_max
  sec3d => wall3d%section(i)
  if (sec3d%patch_in_region) cycle  ! Cannot compute safe box.

  max_area = 0
  do j = 1, 39
    angle = j * pi / 39
    cos_a = cos(angle); sin_a = sin(angle)
    call calc_wall_radius (sec3d%v,  cos_a,  sin_a, radius(1), dr_dtheta)
    call calc_wall_radius (sec3d%v,  cos_a, -sin_a, radius(2), dr_dtheta)
    call calc_wall_radius (sec3d%v, -cos_a,  sin_a, radius(3), dr_dtheta)
    call calc_wall_radius (sec3d%v, -cos_a, -sin_a, radius(4), dr_dtheta)
    rad = minval(radius)
    area = rad**2 * cos_a * sin_a
    if (area > max_area) then
      sec3d%x_safe = 0.999 * rad * cos_a
      sec3d%y_safe = 0.999 * rad * sin_a
      max_area = area
    endif
  enddo

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
type (sr3d_wall_section_input), target :: section(1:100)
type (sr3d_wall_section_struct), pointer :: sec
type (sr3d_wall_section_input), pointer :: sec_in

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
      print *, 'TWO MULTI_SECTION_DEFS HAVE THE SAME NAME: ', trim(name)
      call err_exit
    endif
  enddo

  wall%multi_section(i)%name = name

  n_section = count(section%basic_shape /= '')
  if (any(section(n_section:)%basic_shape /= '')) then
    print *, 'CONFUSED MULTI_SECTION_DEF: ', trim(name)
    call err_exit
  endif
  allocate (wall%multi_section(i)%section(1:n_section))

  do j = 1, n_section
    ix = index(section(j)%basic_shape, ':')
    if (ix == 0) ix = len_trim(section(j)%basic_shape) + 1
    sec => wall%multi_section(i)%section(j)
    sec_in => section(j)
    sec%name                = sec_in%name
    sec%basic_shape(:ix-1)  = sec_in%basic_shape
    sec%shape_name(ix+1:)   = sec_in%basic_shape
    sec%s                   = sec_in%s
    sec%width2              = sec_in%width2
    sec%height2             = sec_in%height2
    sec%width2_plus         = sec_in%width2_plus
    sec%ante_height2_plus   = sec_in%ante_height2_plus
    sec%width2_minus        = sec_in%width2_minus
    sec%ante_height2_minus  = sec_in%ante_height2_minus
  enddo
enddo

end subroutine sr3d_read_wall_multi_section

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_get_emission_pt_params (branch, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
!
! Routine to get the parameters at a photon emission point.
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   branch    -- branch_struct with twiss propagated and mat6s made.
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

subroutine sr3d_get_emission_pt_params (branch, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

use em_field_mod

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele
type (ele_struct) ele_here
type (coord_struct) :: photon
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

ele  => branch%ele(ix_ele)

! Calc the photon's initial twiss values.
! Tracking through a wiggler can take time so use twiss_and_track_intra_ele to
!   minimize the length over which we track.

if (ele_here%ix_ele /= ele%ix_ele .or. ele_here%ix_branch /= ele%ix_branch .or. s_old_offset > s_offset) then
  ele_here = branch%ele(ix_ele-1)
  ele_here%ix_ele = ele%ix_ele
  ele_here%ix_branch = ele%ix_branch
  orb_here = orb(ix_ele-1)
  s_old_offset = 0
endif

call twiss_and_track_intra_ele (ele, branch%param, s_old_offset, s_offset, .true., .true., &
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
  call make_mat6 (ele_here, branch%param, orb_here, orb_here, .true.)
  call track1 (orb_here, ele_here, branch%param, orb1)
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

    call em_field_calc (ele_here, branch%param, ele_here%value(l$), 0.0_rp, orb_here, .false., field)
    gx = field%b(2) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))
    gy = field%b(1) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))

  endif

case default

  ! for mapped wigglers, find the B field at the source point
  ! Note: assumes particles are relativistic!!

  call em_field_calc (ele_here, branch%param, ele_here%value(l$), 0.0_rp, orb_here, .false., field)
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
!   ele_here  -- Ele_struct: Emission is at the exit end of this element (which is a slice_slave).
!   orb_here  -- coord_struct: orbit of particle emitting the photon.
!   gx, gy    -- Real(rp): Horizontal and vertical bending strengths.
!   emit_a    -- Real(rp): Emittance of the a-mode.
!   emit_b    -- Real(rp): Emittance of the b-mode.
!   photon_direction 
!             -- Integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- coord_struct: Generated photon.
!-

subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, p_orb)

implicit none

type (ele_struct), target :: ele_here
type (ele_struct), pointer :: ele
type (coord_struct) :: orb_here
type (coord_struct) :: p_orb
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_e, gx, gy, g_tot, gamma, v2
real(rp) r(3), vec(4), v_mat(4,4)

integer photon_direction

! Get photon energy and "vertical angle".

call convert_total_energy_to (ele_here%value(E_tot$), ele_here%branch%param%particle, gamma) 
call bend_photon_init (gx, gy, gamma, p_orb)

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
p_orb%ix_ele = ele_here%ix_ele

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
p_orb%direction = photon_direction
p_orb%ix_ele = element_at_s(ele_here%branch, ele_here%s, .false.)
p_orb%location = inside$

end subroutine sr3d_emit_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_d_radius (p_orb, branch, d_radius, dw_perp, in_antechamber, check_safe)
!
! Routine to calculate the (transverse) radius of the photon  relative to the wall.
! Optionally can also caluclate the outwrd normal vector perpendicular to the wall.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   p_orb          -- coord_struct): Position.
!   branch         -- branch_struct: Lattice branch containing the wall.
!   check_safe     -- logical, optional: If True, check if photon is safely in the "safe" box far from
!                       the wall. This is used to speed up computations when d_radius value is not needed.
!                       Default is False.
!
! Output:
!   d_radius       -- real(rp): r_photon - r_wall
!   dw_perp(3)     -- real(rp), optional: Outward normal vector perpendicular to the wall.
!   in_antechamber -- Logical, optional: At antechamber wall?
!-

subroutine sr3d_photon_d_radius (p_orb, branch, d_radius, dw_perp, in_antechamber, check_safe)

implicit none

type (sr3d_coord_struct), target :: p_orb
type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d
type (ele_struct), pointer :: ele
type (floor_position_struct) here

real(rp) d_radius, position(6), origin(3)
real(rp), optional :: dw_perp(3)

integer ix

logical, optional :: in_antechamber, check_safe
logical err_flag

! If the photon's position (abs(x), abs(y)) is within the box (x_safe, y_safe) then the photon
! is guaranteed to be within the vacuum chamber.

call sr3d_get_section_index (p_orb, branch)
ix = p_orb%ix_wall_section
wall3d => branch%wall3d

if (logic_option(.false., check_safe)) then
  if (abs(p_orb%orb%vec(1)) < min(wall3d%section(ix)%x_safe, wall3d%section(ix+1)%x_safe) .and. &
      abs(p_orb%orb%vec(3)) < min(wall3d%section(ix)%y_safe, wall3d%section(ix+1)%y_safe)) then
    d_radius = -1
    return
  endif
endif

! 

ele => branch%ele(p_orb%orb%ix_ele)

position = p_orb%orb%vec
position(5) = p_orb%orb%s - branch%ele(p_orb%orb%ix_ele-1)%s

d_radius = wall3d_d_radius (position, ele, dw_perp, ix, in_antechamber, origin, err_flag)

end subroutine sr3d_photon_d_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_section_index (p_orb, branch)
!
! Routine to get the wall index such that 
! For p_orb%orb%vec(6) > 0 (forward motion):
!   wall%section(ix_section)%s < p_orb%s <= wall%section(ix_section+1)%s
! For p_orb%orb%vec(6) < 0 (backward motion):
!   wall%section(ix_section)%s <= p_orb%s < wall%section(ix_section+1)%s
! Exceptions:
!   If p_orb%s == wall%section(1)%s (= 0)       -> ix_section = 0
!   If p_orb%s == wall%section(wall%n_section_max)%s -> ix_section = wall%n_section_max - 1
!
! Input:
!   p_orb  -- sr3d_coord_struct: Photon position.
!   branch -- branch_struct: Lattice branch containing the wall.
!
! Output:
!   p_orb%ix_section -- Integer: Wall index
!-

subroutine sr3d_get_section_index (p_orb, branch)

implicit none

type (sr3d_coord_struct), target :: p_orb
type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d

integer, pointer :: ix_sec
integer n_max

! 

wall3d => branch%wall3d
n_max = ubound(wall3d%section, 1)

ix_sec => p_orb%ix_wall_section
if (ix_sec == not_set$) ix_sec = 1
if (ix_sec == n_max) ix_sec = n_max - 1

if (p_orb%orb%s < wall3d%section(ix_sec)%s .or. p_orb%orb%s > wall3d%section(ix_sec+1)%s) then
  call bracket_index2 (wall3d%section%s, 1, n_max, p_orb%orb%s, ix_sec, ix_sec)
  if (ix_sec == n_max) ix_sec = n_max - 1
endif

! vec(5) at boundary cases

if (p_orb%orb%s == wall3d%section(ix_sec)%s   .and. p_orb%orb%vec(6) > 0 .and. ix_sec /= 1)       ix_sec = ix_sec - 1
if (p_orb%orb%s == wall3d%section(ix_sec+1)%s .and. p_orb%orb%vec(6) < 0 .and. ix_sec /= n_max-1) ix_sec = ix_sec + 1

end subroutine sr3d_get_section_index

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_print_photon_info (photon)
!
! Routine to print information on the photon being tracked.
!
! Input:
!   photon -- sr3d_photon_track_struct: Photon being tracked.
!-

subroutine sr3d_print_photon_info (photon)

type (sr3d_photon_track_struct) photon

print '(8x, a, 3i8, f12.4)', 'Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%orb%p0c
print '(8x, a, 6es13.5, f13.6)', 'Start:  ', photon%start%orb%vec, photon%start%orb%s
print '(8x, a, 6es13.5, f13.6)', 'Now:    ', photon%now%orb%vec, photon%now%orb%s

end subroutine sr3d_print_photon_info

end module
