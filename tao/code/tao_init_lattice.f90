!+
! Subroutine tao_init_lattice (namelist_file, err_flag)
!
! Subroutine to initialize the design lattices.
!
! Input:
!   namelist_file  -- character(*): Name of file containing lattice file info.
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_init_lattice (namelist_file, err_flag)

use tao_interface, except => tao_init_lattice
use tao_input_struct
use ptc_interface_mod

implicit none

type (tao_design_lat_input), target :: design_lattice(0:200)
type (tao_design_lat_input), pointer :: design_lat, dl0
type (lat_struct), pointer :: lat
type (lat_struct) parse_lat
type (ele_struct), pointer :: ele1, ele2
type (tao_universe_struct), pointer :: u, u_work
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orbit(:)

character(*) namelist_file
character(200) full_input_name
character(100) unique_name_suffix, suffix
character(40) name
character(*), parameter :: r_name = 'tao_init_lattice'

integer i_uni, j, k, n, iu, ios, version, ix, key, n_universes, ib, ie, status

logical err_flag, custom_init, combine_consecutive_elements_of_like_name
logical alternative_lat_file_exists
logical err, err1, err2

namelist / tao_design_lattice / design_lattice, &
       combine_consecutive_elements_of_like_name, unique_name_suffix, &
       n_universes

! Defaults

err_flag = .true.
design_lattice = tao_design_lat_input()
design_lattice(0)%file = 'Garbage!!!'

alternative_lat_file_exists = (s%init%hook_lat_file /= '' .or. s%init%lattice_file_arg /= '')

! Read lattice info

if (associated(tao_hook_init_read_lattice_info_ptr)) call tao_hook_init_read_lattice_info_ptr (namelist_file)

if (s%com%init_read_lat_info) then
  ! namelist_file == '' means there is no lattice file so just use the defaults.

  if (namelist_file /= '') then
    call tao_open_file (namelist_file, iu, full_input_name, s_fatal$)
    call out_io (s_blank$, r_name, '*Init: Opening Lattice Info File: ' // namelist_file)
    if (iu == 0) then
      call out_io (s_fatal$, r_name, 'ERROR OPENING TAO LATTICE INFO FILE. WILL EXIT HERE...')
      return
    endif
  endif

  ! Defaults

  n_universes = s%com%n_universes
  combine_consecutive_elements_of_like_name = .false.
  unique_name_suffix = ''

  if (namelist_file /= '') then
    read (iu, nml = tao_design_lattice, iostat = ios)
    if (ios > 0 .or. (ios < 0 .and. .not. alternative_lat_file_exists)) then
      call out_io (s_abort$, r_name, 'TAO_DESIGN_LATTICE NAMELIST READ ERROR.')
      rewind (iu)
      do
        read (iu, nml = tao_design_lattice)  ! force printing of error message
      enddo
    endif
    close (iu)
  endif

  s%com%combine_consecutive_elements_of_like_name = combine_consecutive_elements_of_like_name
  s%com%n_universes = n_universes
  s%init%unique_name_suffix = trim(unique_name_suffix)
endif

!

allocate (s%u(s%com%n_universes))
nullify (s%com%u_working)

! Read in the lattices

do i_uni = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i_uni)
  u%is_on = .true.          ! turn universe on
  u%ix_uni = i_uni
  u%calc = tao_universe_calc_struct()

  ! Get the name of the lattice file
  ! The presedence is: 
  !   From the command line (preferred).
  !   From the tao_design_lattice namelist.
  !   Specified by hook code.

  design_lat => design_lattice(i_uni)

  ! Element range?

  if (design_lat%use_element_range(1) /= '') then
    design_lat%slice_lattice = trim(design_lat%use_element_range(1)) // ':' // trim(design_lat%use_element_range(2))
    call out_io (s_warn$, r_name, 'In the tao_design_lattice namelist in the init file: ' // namelist_file, &
                '"design_lattice(i_uni)%use_element_range" is now "design_lattice(i_uni)%slice_lattice".', &
                'Please modify your file.')
  endif

  !

  if (s%init%lattice_file_arg /= '') then
    design_lat%file = s%init%lattice_file_arg
    design_lat%language = ''
    design_lat%file2 = ''
  elseif (s%init%hook_lat_file /= '' .and. design_lat%file == '') then
    design_lat%file = s%init%hook_lat_file
    design_lat%language = ''
    design_lat%file2 = ''
  endif

  if (design_lat%file == '' .and. i_uni > 1) design_lat = design_lattice(i_uni-1)

  ix = index (design_lat%file, '|') ! Indicates multiple lattices
  if (ix /= 0) then
    design_lat%file2 = design_lat%file(ix+1:)
    design_lat%file  = design_lat%file(1:ix-1)
  endif

  ! Split off language name and/or use_line if needed

  ix = index(design_lat%file(1:10), '::')

  if (ix == 0) then
    design_lat%language = 'bmad'
  else
    design_lat%language = design_lat%file(1:ix-1)
    design_lat%file = design_lat%file(ix+2:)
  endif

  !

  ix = index(design_lat%file, '@')
  if (ix /= 0) then
    design_lat%use_line = design_lat%file(ix+1:)
    design_lat%file = design_lat%file(:ix-1)
  else
    design_lat%use_line = ''
  endif

  ! If the lattice file was obtained from the tao init file and if the lattice name 
  ! is relative, then it is relative to the directory where the tao init file is.

  if (.not. alternative_lat_file_exists .and. file_name_is_relative(design_lat%file)) &
                design_lat%file = trim(s%init%init_file_arg_path) // design_lat%file 

  ! Read in the design lattice. 
  ! A blank means use the lattice form universe 1.

  allocate (u%design, u%base, u%model)
  dl0 => design_lattice(i_uni-1)
  if (design_lat%file == dl0%file .and. design_lat%slice_lattice == dl0%slice_lattice .and. &
            design_lat%file2 == dl0%file2 .and. design_lat%use_line == dl0%use_line .and. &
           (design_lat%reverse_lattice .eqv. dl0%reverse_lattice) .and. &
            design_lat%start_branch_at == dl0%start_branch_at) then
    u%design_same_as_previous = .true.
    u%design%lat = s%u(i_uni-1)%design%lat
  else
    u%design_same_as_previous = .false.
    select case (design_lat%language)
    case ('bmad')
      call out_io (s_blank$, r_name, 'Reading Bmad file: ' // design_lat%file)
      if (design_lat%use_line /= '') call out_io (s_blank$, r_name, '  Line used is: ' // design_lat%use_line)
      call bmad_parser (design_lat%file, u%design%lat, use_line = design_lat%use_line, err_flag = err, parse_lat = parse_lat)
    case ('digested')
      call out_io (s_blank$, r_name, "Reading digested Bmad file: " // trim(design_lat%file))
      call read_digested_bmad_file (design_lat%file, u%design%lat, version, err)
    case default
      call out_io (s_abort$, r_name, 'LANGUAGE NOT RECOGNIZED: ' // design_lat%language)
      return
    end select

    if (err .and. .not. s%global%debug_on) then
      call out_io (s_fatal$, r_name, &
              'PARSER ERROR DETECTED FOR UNIVERSE: \i0\ ', &
              'EXITING...', i_array = [i_uni])
      return
    endif

    ! Call bmad_parser2 if wanted

    if (design_lat%file2 /= '') then
      call bmad_parser2 (design_lat%file2, u%design%lat, parse_lat = parse_lat)
    endif

    if (design_lat%language == 'bmad') call deallocate_lat_pointers(parse_lat)

    !

    if (s%com%combine_consecutive_elements_of_like_name) call combine_consecutive_elements(u%design%lat, err)
    if (s%global%symbol_import) call tao_symbol_import_from_lat(u%design%lat)

    unique_name_suffix = s%init%unique_name_suffix
    do
      if (unique_name_suffix == '') exit
      call string_trim(unique_name_suffix, unique_name_suffix, ix)
      name = unique_name_suffix(1:ix)
      unique_name_suffix = unique_name_suffix(ix+1:)
      if (index(name, '##') /= 0) then
        call out_io (s_error$, r_name, 'USE OF "##" IN UNIQUE_NAME_SUFFIX CONFLICTS WITH BMAD "##" CONSTRUCT TO', &
                                       'IDENTIFY THE Nth ELEMENT OF A GIVEN NAME. SUFFIX NOT APPLIED.')
        cycle
      endif

      ix = index(name, '::')
      if (ix == 0) then
        call create_unique_ele_names (u%design%lat, 0, name)
      else
        j = key_name_to_key_index(name(1:ix-1), .true.)
        call create_unique_ele_names (u%design%lat, j, name(ix+2:))
      endif
    enddo

    if (s%init%start_branch_at_arg /= '') design_lat%start_branch_at = s%init%start_branch_at_arg
    if (s%init%slice_lattice_arg /= '') design_lat%slice_lattice = s%init%slice_lattice_arg

    if (design_lat%start_branch_at /= '') then
      call start_branch_at (u%design%lat, design_lat%start_branch_at, .false., err)
    endif

    if (design_lat%slice_lattice /= '') then
      call slice_lattice (u%design%lat, design_lat%slice_lattice, err)
      if (err) call out_io (s_error$, r_name, 'ERROR SLICING LATTICE FOR UNIVERSE: ' // int_str(i_uni))
    endif

    if (design_lat%reverse_lattice .neqv. (s%init%reverse_arg == present_str)) then
      lat => u%design%lat
      call reallocate_coord (orbit, lat)
      call init_coord (orbit(0), lat%particle_start, lat%ele(0), downstream_end$)
      call twiss_and_track(lat, orbit, status)
      if (status /= ok$) then
        call out_io (s_error$, r_name, 'PROBLEM CALCULATING TWISS/ORBIT FOR UNIVERSE: '// int_str(i_uni)) 
      endif
      call reverse_lat(lat, lat)
      lat%particle_start = orbit(lat%n_ele_track)
      lat%particle_start%vec(2) = -lat%particle_start%vec(2)
      lat%particle_start%vec(4) = -lat%particle_start%vec(4)
      lat%particle_start%vec(5) = -lat%particle_start%vec(5)
      if (lat%particle_start%beta == 0) then
        lat%particle_start%t = lat%ele(0)%ref_time
      else
        lat%particle_start%t = lat%ele(0)%ref_time - lat%particle_start%vec(5) / (lat%particle_start%beta * c_light)
      endif
    endif
  endif

  !

  u%calc%one_turn_map     = design_lat%one_turn_map_calc
  u%calc%dynamic_aperture = design_lat%dynamic_aperture_calc

  ! Custom stuff

  if (associated(tao_hook_init_lattice_post_parse_ptr)) call tao_hook_init_lattice_post_parse_ptr (u)

  ! In case there is a match element with recalc = T, propagate the twiss parameters which
  ! makes the beginning Twiss of the match element match the end Twiss of the previous element

  if (u%design%lat%param%geometry == open$) then
    call twiss_propagate_all (u%design%lat)
  endif

  ! Init model, base, and u%ele

  n = ubound(u%design%lat%branch, 1)
  allocate (u%model%tao_branch(0:n))
  allocate (u%design%tao_branch(0:n))
  allocate (u%base%tao_branch(0:n))
  allocate (u%model_branch(0:n))

  do ib = 0, ubound(u%design%lat%branch, 1)
    u%design%tao_branch(ib)%modes_ri%a%emittance = u%design%lat%branch(ib)%a%emit
    u%design%tao_branch(ib)%modes_ri%b%emittance = u%design%lat%branch(ib)%b%emit
    n = u%design%lat%branch(ib)%n_ele_max

    allocate (u%model%tao_branch(ib)%orbit(0:n), u%model%tao_branch(ib)%bunch_params(0:n))
    allocate (u%design%tao_branch(ib)%orbit(0:n), u%design%tao_branch(ib)%bunch_params(0:n))
    allocate (u%base%tao_branch(ib)%orbit(0:n), u%base%tao_branch(ib)%bunch_params(0:n))

    allocate (u%model%tao_branch(ib)%lat_sigma(0:n))
    allocate (u%design%tao_branch(ib)%lat_sigma(0:n))
    allocate (u%base%tao_branch(ib)%lat_sigma(0:n))

    allocate (u%model_branch(ib)%ele(-1:n))
  enddo

  u%model = u%design
  u%base  = u%design

  u%design%name = 'design'
  u%model%name  = 'model'
  u%base%name   = 'base'

  u%design%u => u
  u%model%u => u
  u%base%u => u

  u%model%tao_branch(0)%orb0  = u%model%lat%particle_start
  u%design%tao_branch(0)%orb0 = u%design%lat%particle_start
  u%base%tao_branch(0)%orb0   = u%base%lat%particle_start

  ! Check for match element with recalc = True

  do ib = 0, ubound(u%design%lat%branch, 1)
    branch => u%design%lat%branch(ib)
    u%model%tao_branch(ib)%has_open_match_element = .false.
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (ele%key /= match$) cycle
      if (.not. is_true(ele%value(recalc$))) cycle
      if (nint(ele%value(matrix$)) /= match_twiss$ .and. nint(ele%value(matrix$)) /= phase_trombone$ .and. &
                                        nint(ele%value(kick0$)) /= match_orbit$) cycle  
      u%model%tao_branch(ib)%has_open_match_element = .true.
      exit
    enddo
    u%design%tao_branch(ib)%has_open_match_element = u%model%tao_branch(ib)%has_open_match_element
    u%base%tao_branch(ib)%has_open_match_element = u%model%tao_branch(ib)%has_open_match_element

    u%design%tao_branch(ib)%tao_lat => u%design
    u%model%tao_branch(ib)%tao_lat  => u%model
    u%base%tao_branch(ib)%tao_lat   => u%base
  enddo

  call ele_order_calc(u%design%lat, u%ele_order)

enddo

err_flag = .false.

end subroutine tao_init_lattice
