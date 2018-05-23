!+
! Subroutine tao_init_lattice (input_file_name)
!
! Subroutine to initialize the design lattices.
!
! Input:
!   input_file_name  -- character(*): file name containing lattice file
!                                             namestructs
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_init_lattice (input_file_name)

use tao_interface, except => tao_init_lattice
use tao_input_struct
use ptc_interface_mod

implicit none

type (tao_design_lat_input) design_lattice(0:200), design_lat
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele1, ele2
type (tao_universe_struct), pointer :: u, u_work
type (branch_struct), pointer :: branch

character(*) input_file_name
character(200) full_input_name, init_lat_file
character(40) unique_name_suffix, suffix
character(20) :: r_name = 'tao_init_lattice'

integer i, j, k, n, iu, ios, version, ix, key, n_universes, ib, ie

logical custom_init, combine_consecutive_elements_of_like_name
logical common_lattice
logical err, err1, err2, err3

namelist / tao_design_lattice / design_lattice, &
       combine_consecutive_elements_of_like_name, unique_name_suffix, &
       common_lattice, n_universes

! Defaults

design_lattice = tao_design_lat_input()
design_lat = tao_design_lat_input()

! Read lattice info

call tao_hook_init_read_lattice_info (input_file_name)

if (s%com%init_read_lat_info) then

  ! input_file_name == '' means there is no lattice file so just use the defaults.

  if (input_file_name /= '') then
    call tao_open_file (input_file_name, iu, full_input_name, s_fatal$)
    call out_io (s_blank$, r_name, '*Init: Opening Lattice Info File: ' // input_file_name)
    if (iu == 0) then
      call out_io (s_fatal$, r_name, 'ERROR OPENING TAO LATTICE INFO FILE. WILL EXIT HERE...')
      call err_exit
    endif
  endif

  ! Defaults

  common_lattice = .false.
  n_universes = s%com%n_universes
  combine_consecutive_elements_of_like_name = .false.
  unique_name_suffix = ''

  if (input_file_name /= '') then
    read (iu, nml = tao_design_lattice, iostat = ios)
    if (ios > 0 .or. (ios < 0 .and. s%com%lat_file == '')) then
      call out_io (s_abort$, r_name, 'TAO_DESIGN_LATTICE NAMELIST READ ERROR.')
      rewind (iu)
      do
        read (iu, nml = tao_design_lattice)  ! force printing of error message
      enddo
    endif
    close (iu)
  endif

  s%com%combine_consecutive_elements_of_like_name = combine_consecutive_elements_of_like_name
  s%com%common_lattice = common_lattice
  s%com%n_universes = n_universes
  s%com%unique_name_suffix = unique_name_suffix
endif

!

if (s%com%common_lattice) then
  allocate (s%u(0:s%com%n_universes))
  allocate (s%com%u_working)

  if (any(design_lattice(2:)%file /= '')) then
    call out_io (s_fatal$, r_name, 'ONLY ONE LATTICE MAY BE SPECIFIED WHEN USING COMMON_LATTICE')
    call err_exit
  endif

else
  allocate (s%u(s%com%n_universes))
  nullify (s%com%u_working)
endif

! Read in the lattices

init_lat_file = s%com%lat_file

do i = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i)
  u%is_on = .true.          ! turn universe on
  u%ix_uni = i
  u%calc%rad_int_for_data          = .false.
  u%calc%srdt_for_data             = 0
  u%calc%rad_int_for_plotting      = .false.
  u%calc%chrom_for_data            = .false.
  u%calc%chrom_for_plotting        = .false.
  u%calc%beam_sigma_for_data       = .false.
  u%calc%beam_sigma_for_plotting   = .false.

  ! If unified then only read in a lattice for the common universe.

  if (s%com%common_lattice .and. i /= ix_common_uni$) cycle

  ! Get the name of the lattice file

  if (design_lattice(i)%file /= '' .and. (init_lat_file(1:6) == 'HOOK::' .or. init_lat_file == '')) then
    design_lat = design_lattice(i)

  elseif (init_lat_file /= '') then
    if (init_lat_file(1:6) == 'HOOK::') init_lat_file = init_lat_file(7:)
    ix = index (init_lat_file, '|') ! Indicates multiple lattices
    if (ix == 0) then
      design_lat%file = init_lat_file
    else
      design_lat%file = init_lat_file(1:ix-1)
      init_lat_file = init_lat_file(ix+1:)
    endif
    design_lat%language = ''
    design_lat%file2 = ''
  endif

  ! Can happen when design_lattice(1) is set and not design_lattice(0)

  if (s%com%common_lattice .and. design_lat%file == '') then
    design_lat = design_lattice(i+1)
  endif

  ! Split off language name and/or use_line if needed

  ix = index(design_lat%file(1:10), '::')

  if (ix == 0) then
    if (index(design_lat%file, '.xsif ') /= 0) then
      design_lat%language = 'xsif'
    else
      design_lat%language = 'bmad'
    endif
  else
    design_lat%language = design_lat%file(1:ix-1)
    design_lat%file = design_lat%file(ix+2:)
  endif

  ! "#reverse" construct is old deprecated syntax.

  ix = index(design_lat%file, '#reverse')
  if (ix /= 0) then
    design_lat%reverse_tracking = .true.
    design_lat%file = design_lat%file(1:ix-1) // design_lat%file(ix+8:)
  endif

  u%reverse_tracking = design_lat%reverse_tracking

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

  if (s%com%lat_file == '' .and. file_name_is_relative(design_lat%file)) &
                design_lat%file = trim(s%com%init_tao_file_path) // design_lat%file 

  ! Read in the design lattice. 
  ! A blank means use the lattice form universe 1.

  allocate (u%design, u%base, u%model)
  select case (design_lat%language)
  case ('bmad')
    call bmad_parser (design_lat%file, u%design%lat, use_line = design_lat%use_line, err_flag = err)
  case ('xsif')
    call xsif_parser (design_lat%file, u%design%lat, use_line = design_lat%use_line, err_flag = err)
  case ('aml')
    call aml_parser (design_lat%file, u%design%lat, err_flag = err)
  case ('digested')
    call out_io (s_blank$, r_name, "Reading digested BMAD file " // trim(design_lat%file))
    call read_digested_bmad_file (design_lat%file, u%design%lat, version, err)
  case default
    call out_io (s_abort$, r_name, 'LANGUAGE NOT RECOGNIZED: ' // design_lat%language)
    call err_exit
  end select

  if (err .and. .not. s%global%debug_on) then
    call out_io (s_fatal$, r_name, &
            'PARSER ERROR DETECTED FOR UNIVERSE: \i0\ ', &
            'EXITING...', i_array = (/ i /))
    if (s%global%stop_on_error) stop
  endif

  if (s%com%combine_consecutive_elements_of_like_name) call combine_consecutive_elements(u%design%lat)

  if (s%com%unique_name_suffix /= '') then
    call tao_string_to_element_id (s%com%unique_name_suffix, key, suffix, err, .true.)
    if (err) call err_exit
    call create_unique_ele_names (u%design%lat, key, suffix)
  endif

  ! Element range?

  if (design_lat%use_element_range(1) /= '') then
    lat => u%design%lat
    call pointer_to_this_ele (1, design_lat%use_element_range(1), err1, ele1)
    call pointer_to_this_ele (2, design_lat%use_element_range(2), err2, ele2)

    err3 = (ubound(u%design%lat%branch, 1) > 0)
    if (err3) then
      call out_io (s_fatal$, r_name, 'USE_ELEMENT_RANGE ERROR. CURRENTLY RANGES CAN ONLY BE USED WITH LATTICES WITH A SINGLE BRANCH')
      if (s%global%stop_on_error) stop
    endif

    if (.not. err1 .and. .not. err2 .and. .not. err3) then
      if (ele1%ix_ele > ele2%ix_ele) then
        call out_io (s_fatal$, r_name, 'USE_ELEMENET_RANGE_ERROR. END ELEMENT1 IS AFTER END ELEMENT2: ' // &
                          trim(design_lat%use_element_range(1)) // ', ' // trim(design_lat%use_element_range(2)))
        if (s%global%stop_on_error) stop

      else
        lat%ele(1:ele1%ix_ele-1)%key = -1
        lat%ele(ele2%ix_ele+1:lat%n_ele_track)%key = -1
        call remove_eles_from_lat (lat)
      endif
    endif
  endif

  ! Call bmad_parser2 if wanted

  if (design_lat%file2 /= '') then
    call bmad_parser2 (design_lat%file2, u%design%lat)
  endif

  !

  u%calc%one_turn_map     = design_lat%one_turn_map_calc
  u%calc%dynamic_aperture = design_lat%dynamic_aperture_calc

  ! Custom stuff

  call tao_hook_init_lattice_post_parse (u)

  ! In case there is a match element with match_end = T, propagate the twiss parameters which
  ! makes the beginning Twiss of the match element match the end Twiss of the previous element

  if (u%design%lat%param%geometry == open$) then
    call twiss_propagate_all (u%design%lat)
  endif

  ! Init model, base, and u%ele

  n = ubound(u%design%lat%branch, 1)
  allocate (u%model%tao_branch(0:n))
  allocate (u%design%tao_branch(0:n))
  allocate (u%base%tao_branch(0:n))
  allocate (u%uni_branch(0:n))

  do ib = 0, ubound(u%design%lat%branch, 1)
    u%design%tao_branch(ib)%modes%a%emittance = u%design%lat%branch(ib)%a%emit
    u%design%tao_branch(ib)%modes%b%emittance = u%design%lat%branch(ib)%b%emit
    n = u%design%lat%branch(ib)%n_ele_max
    allocate (u%model%tao_branch(ib)%orbit(0:n), u%model%tao_branch(ib)%bunch_params(0:n))
    allocate (u%design%tao_branch(ib)%orbit(0:n), u%design%tao_branch(ib)%bunch_params(0:n))
    allocate (u%base%tao_branch(ib)%orbit(0:n), u%base%tao_branch(ib)%bunch_params(0:n))
    allocate (u%uni_branch(ib)%ele(-1:n))
  enddo

  u%model = u%design
  u%base  = u%design

  u%design%name = 'design'
  u%model%name  = 'model'
  u%base%name   = 'base'

  u%model%tao_branch(0)%orb0  = u%model%lat%beam_start
  u%design%tao_branch(0)%orb0 = u%design%lat%beam_start
  u%base%tao_branch(0)%orb0   = u%base%lat%beam_start

  ! Check for match element with match_end = True

  do ib = 0, ubound(u%design%lat%branch, 1)
    branch => u%design%lat%branch(ib)
    u%model%tao_branch(ib)%has_open_match_element = .false.
    do ie = 1, branch%n_ele_track
      if (branch%ele(ie)%key /= match$) cycle
      if (.not. is_true(branch%ele(ie)%value(match_end$)) .and. &
          .not. is_true(branch%ele(ie)%value(match_end_orbit$))) cycle
      u%model%tao_branch(ib)%has_open_match_element = .true.
      exit
    enddo
    u%design%tao_branch(ib)%has_open_match_element = u%model%tao_branch(ib)%has_open_match_element
    u%base%tao_branch(ib)%has_open_match_element = u%model%tao_branch(ib)%has_open_match_element
  enddo

enddo

! Working lattice setup

if (s%com%common_lattice) then

  u_work => s%com%u_working
  u_work%common     => s%u(ix_common_uni$)
  allocate (u_work%design, u_work%base, u_work%model)
  u_work%design%lat = u_work%common%design%lat
  u_work%base%lat   = u_work%common%base%lat
  u_work%model%lat  = u_work%common%model%lat

  n = ubound(u_work%design%lat%branch, 1)
  allocate (u_work%model%tao_branch(0:n))
  allocate (u_work%design%tao_branch(0:n))
  allocate (u_work%base%tao_branch(0:n))
  allocate (u_work%uni_branch(0:n))

  do k = 0, ubound(u_work%design%lat%branch, 1)
    n = u_work%design%lat%branch(k)%n_ele_max
    allocate (u_work%model%tao_branch(k)%orbit(0:n), u_work%model%tao_branch(k)%bunch_params(0:n))
    allocate (u_work%design%tao_branch(k)%orbit(0:n), u_work%design%tao_branch(k)%bunch_params(0:n))
    allocate (u_work%base%tao_branch(k)%orbit(0:n), u_work%base%tao_branch(k)%bunch_params(0:n))
    allocate (u_work%uni_branch(k)%ele(-1:n))
  enddo

  ! If unified then point back to the common universe (#1) and the working universe (#2)

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (i == ix_common_uni$) cycle
    u => s%u(i)
    u%common     => s%u(ix_common_uni$)
    u%uni_branch => s%u(ix_common_uni$)%uni_branch
    u%design => s%u(ix_common_uni$)%design
    u%base   => s%u(ix_common_uni$)%base
    u%model  => s%com%u_working%model
  enddo

endif

!------------------------------------------------------------------
contains

subroutine pointer_to_this_ele (ix_range, ele_name, err, ele)

type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
character(*) ele_name
integer n_loc, ix_range
logical err

!

call lat_ele_locator (ele_name, u%design%lat, eles, n_loc, err)

if (err .or. n_loc == 0) then
  call out_io (s_fatal$, r_name, 'USE_ELEMENT_RANGE ERROR. CANNOT FIND ELEMENT IN LATTICE: ' // ele_name)
  if (s%global%stop_on_error) stop
  err = .true.
  return
endif

if (n_loc > 1) then
  call out_io (s_fatal$, r_name, 'USE_ELEMENT_RANGE ERROR. MULTIPLE ELEMENTS IN LATTICE: ' // ele_name)
  if (s%global%stop_on_error) stop
  err = .true.
  return
endif

ele => eles(1)%ele
if (ele%lord_status == super_lord$) then
  if (ix_range == 1) ele => pointer_to_slave(ele, 1)
  if (ix_range == 2) ele => pointer_to_slave(ele, ele%n_slave)
endif

if (ele%ix_ele > u%design%lat%n_ele_track) then
  call out_io (s_fatal$, r_name, 'USE_ELEMENT_RANGE ERROR. ELEMENT IN: ' // ele_name)
  if (s%global%stop_on_error) stop
  err = .true.
  return
endif

end subroutine pointer_to_this_ele

end subroutine tao_init_lattice
