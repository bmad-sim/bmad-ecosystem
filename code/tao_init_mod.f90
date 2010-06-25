module tao_init_mod

use tao_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_global (init_file)
!
! Subroutine to initialize the tao global structures.
! If init_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   init_file  -- Character(*): Tao initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_global (init_file)

!use tao_lattice_calc_mod
!use bmad_parser_mod
use random_mod
use csr_mod, only: csr_param
use opti_de_mod, only: opti_de_param

implicit none

type (tao_global_struct), save :: global, default_global

integer ios, iu, i, j, k, ix, num
integer n_data_max, n_var_max
integer n_d2_data_max, n_v1_var_max ! Deprecated variables
integer n, iostat

character(*) init_file
character(40) :: r_name = 'tao_init_global'
character(200) file_name
character(40) name, universe

character(100) line

logical err
logical, save :: init_needed = .true.

namelist / tao_params / global, bmad_com, csr_param, opti_de_param, &
          n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  
!-----------------------------------------------------------------------
! First time through capture the default global (could have been set via command line arg.)

if (init_needed) then
  default_global = s%global
  init_needed = .false.
endif

global = default_global         ! establish defaults

call tao_hook_init_global (init_file, global)

! read global structure from tao_params namelist
! init_file == '' means there is no lattice file so just use the defaults.

if (init_file /= '') then
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening Init File: ' // file_name)
  if (iu == 0) then
    call out_io (s_blank$, r_name, "Note: Cannot open init file for tao_params namelist read")
  else
    call out_io (s_blank$, r_name, 'Init: Reading tao_params namelist')
    bmad_com%rel_tolerance = 1e-8   ! Need tighter tol for calculating derivatives
    bmad_com%abs_tolerance = 1e-11  ! Need tighter tol for calculating derivatives
    read (iu, nml = tao_params, iostat = ios)
    if (ios > 0) then
      call out_io (s_error$, r_name, 'ERROR READING TAO_PARAMS NAMELIST.')
      rewind (iu)
      read (iu, nml = tao_params)  ! To give error message
    endif
    if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_params namelist found')
    close (iu)
  endif
endif

! transfer global to s%global

s%global = global
if (tao_com%noplot_arg_found) s%global%plot_on = .false.

! This for backwards compatibility.

if (global%n_curve_pts > 0) s%plot_page%n_curve_pts = global%n_curve_pts  

if (s%global%track_type == "macro") then
  call out_io (s_error$, r_name, &
             'MACROPARTICLE TRACKING IS NOT ACTIVELY SUPPORTED!', &
             'PLEASE USE "beam" TRACKING INSTEAD.', &
             'IF YOU NEED MACROPARTICLE TRACKING PLEASE SEE DAVID SAGAN.')
endif

!-----------------------------------------------------------------------
! Tao does its own bookkeeping

bmad_com%auto_bookkeeper = .false.
tao_com%valid_plot_who(1:5) = (/ 'model ', 'base  ', 'ref   ', 'design', 'meas  ' /)

! Seed random number generator

call ran_seed_put (s%global%random_seed)
call ran_engine (s%global%random_engine)
call ran_gauss_converter (s%global%random_gauss_converter, s%global%random_sigma_cutoff)

end subroutine tao_init_global

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_connected_universes (init_file)
!
! Subroutine to initialize universe connections for
! continuing tracking from one universe to another.
!
! If init_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   init_file  -- Character(*): Tao initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_connected_universes (init_file)

use tao_input_struct

implicit none

type (tao_universe_struct), pointer :: u
type (tao_connected_uni_input) connect

integer i, k, iu, ios, ib, n_uni
integer n, iostat, ix_universe, to_universe
integer ix_track_start, ix_track_end

character(*) init_file
character(40) :: r_name = 'tao_init_connected_universes'
character(200) file_name

logical err

namelist / tao_connected_uni_init / to_universe, connect
  
!-----------------------------------------------------------------------
! Init connected universes

! defaults

do i = lbound(s%u, 1), ubound(s%u, 1)
  s%u(i)%connect%connected = .false.
  s%u(i)%connect%match_to_design = .false.
  s%u(i)%connect%from_uni = -1
  s%u(i)%connect%to_uni = -1
  s%u(i)%connect%from_uni_s = -1
  s%u(i)%connect%from_uni_ix_ele = -1
enddo

call tao_hook_init_connected_uni ()

if (.not. tao_com%init_connected_uni .or. init_file == '') return

!

call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)

do
  to_universe = -1
  connect%match_to_design = .false.
  connect%from_universe = -1
  connect%at_element = ''
  connect%at_ele_index = -1
  connect%at_s = -1

  read (iu, nml = tao_connected_uni_init, iostat = ios)
  if (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_CONNECTED_UNI_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_connected_uni_init)  ! generate an error message
    enddo
  endif
  if (ios < 0) exit

  if (to_universe == -1) then
    call out_io (s_abort$, r_name, &
          'INIT: READ TAO_CONNECTED_UNI_INIT NAMELIST HAS NOT SET TO_UNIVERSE!')
    call err_exit
  endif
  call out_io (s_blank$, r_name, &
         'Init: Read tao_connected_uni_init namelist for universe \i3\ ', to_universe)
  call init_connected_uni (connect, to_universe)

enddo

close (iu)

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

subroutine init_connected_uni (connect, this_uni_index)

implicit none

type (tao_universe_struct), pointer :: u_to
type (tao_universe_struct), pointer ::  u_from
type (tao_connected_uni_input) connect 
integer this_uni_index

character(40) class, ele_name
character(40) :: r_name = 'init_connected_uni'

integer j, ix

!

u_from => s%u(connect%from_universe)
u_to   => s%u(to_universe)

if (connect%from_universe .eq. 0 .or. connect%at_element .eq. "none") then
  u_to%connect%connected = .false.
  return
endif

if (connect%from_universe .ge. this_uni_index) then
  call out_io (s_abort$, r_name, &
        "A UNIVERSE CAN ONLY INJECT INTO A UNIVERSE WITH A GREATER UNIVERSE INDEX!")
  call err_exit
endif
  
u_to%connect%connected = .true.
u_to%connect%from_uni = connect%from_universe
u_from%connect%to_uni = to_universe

call init_ele (u_to%connect%match_ele)
u_to%connect%match_ele%key = match$
u_to%connect%match_to_design = connect%match_to_design
          
! find extraction element
call string_trim (connect%at_element, ele_name, ix)
call str_upcase (ele_name, ele_name)
if (ix /= 0) then
  if (connect%at_s /= -1 .or. connect%at_ele_index /= -1) then
    call out_io (s_error$, r_name, &
        "CANNOT SPECIFY AN ELEMENT, IT'S INDEX OR POSITION AT SAME TIME!", &
        "WILL USE ELEMENT NAME.")
  endif
  if (ele_name == "END") then
    u_to%connect%from_uni_s  = u_from%design%lat%ele(u_from%design%lat%n_ele_track)%s
    u_to%connect%from_uni_ix_ele = u_from%design%lat%n_ele_track
  else
    ! using element name 
    ! find last element with name
    do j = u_from%design%lat%n_ele_track, 0, -1
      if (ele_name == trim(u_from%design%lat%ele(j)%name)) then
        u_to%connect%from_uni_s = u_from%design%lat%ele(j)%s
        u_to%connect%from_uni_ix_ele = j
        return
      endif
      if (j == 0) then
        call out_io (s_abort$, r_name, &
                    "COULDN'T FIND CONNECTION ELEMENT IN UNIVERSE \I\ ", &
                     connect%from_universe)
        call err_exit
      endif
    enddo
  endif
elseif (connect%at_ele_index /= -1) then
  if (connect%at_s /= -1) then
    call out_io (s_error$, r_name, &
        "YOU CANNOT SPECIFY AN ELEMENT, IT'S INDEX OR POSITION AT SAME TIME!", &
        "WILL USE ELEMENT INDEX.")
  endif
    u_to%connect%from_uni_s = u_from%design%lat%ele(connect%at_ele_index)%s
    u_to%connect%from_uni_ix_ele = connect%at_ele_index
else
  ! using s position
  if (s%global%track_type /= 'single' ) then
    call out_io (s_abort$, r_name, &
     "AN ARBITRARY S POSITION FOR COUPLING IF NOT TRACKING A SINGLE PARTICLE IS ILLEGAL")
    call err_exit
  endif
  !FIX_ME: get ix_ele for element right before this s position
  u_to%connect%from_uni_s = connect%at_s
endif

end subroutine init_connected_uni

end subroutine tao_init_connected_universes

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_beams (init_file)
!
! Subroutine to initialize beam stuff.
!
! If init_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   init_file  -- Character(*): Tao initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_beams (init_file)

use spin_mod
use tao_input_struct

implicit none

type (tao_universe_struct), pointer :: u
type (beam_init_struct) beam_init
type (spin_polar_struct) spin

integer i, k, iu, ios, ib, n_uni
integer n, iostat, ix_universe, to_universe
integer ix_track_start, ix_track_end

character(*) init_file
character(40) :: r_name = 'tao_init_beams'
character(160) beam_saved_at
character(200) file_name, beam0_file, beam_all_file
character(60), target :: save_beam_at(100)   ! old style syntax

logical err

namelist / tao_beam_init / ix_universe, beam0_file, &
  ix_track_start, ix_track_end, beam_all_file, beam_init, beam_saved_at, &
  beam_saved_at
         
!-----------------------------------------------------------------------
! Init Beams

call tao_hook_init_beam ()

if (.not. tao_com%init_beam .or. init_file == '') return

!

call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

do i = lbound(s%u, 1), ubound(s%u, 1)
  do ib = 0, ubound(s%u(i)%uni_branch, 1)
    s%u(i)%beam_info%beam0_file = ''
  enddo
enddo

do 

  ! defaults
  ix_universe = -1
  beam_init%distribution_type = ''
  beam_init%a_norm_emitt  = 0.0
  beam_init%b_norm_emitt  = 0.0
  beam_init%dPz_dz = 0.0
  beam_init%center(:) = 0.0
  beam_init%bunch_charge = 0.0
  beam_init%dt_bunch = 0
  beam_init%sig_z   = 0.0
  beam_init%sig_e   = 0.0
  beam_init%renorm_center = .true.
  beam_init%renorm_sigma = .true.
  beam_init%n_bunch = -1
  beam_init%n_particle  = -1
  beam0_file = tao_com%beam0_file        ! From the command line
  beam_all_file = tao_com%beam_all_file  ! From the command line
  beam_saved_at = ''
  save_beam_at  = ''
  ix_track_start = 0
  ix_track_end = -1

  ! Read beam parameters

  read (iu, nml = tao_beam_init, iostat = ios)
  if (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_BEAM_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_beam_init)  ! generate an error message
    enddo
  endif
  ! transfer info from old style save_beam_at(:) to beam_saved_at
  do i = 1, size(save_beam_at)
    if (save_beam_at(i) == '') cycle
    beam_saved_at = trim(beam_saved_at) // ', ' // trim(save_beam_at(i))
  enddo

  if (ios /= 0) exit

  ! Error checking

  if (beam_init%n_bunch < 1) then
    call out_io (s_fatal$, r_name, &
      'BEAM_INIT%N_BUNCH NOT PROPERLY SET.')
    call err_exit
  endif

  ! init

  call out_io (s_blank$, r_name, &
        'Init: Read tao_beam_init namelist for universe \i3\ ', ix_universe)
  if (ix_universe == -1) then
    do i = lbound(s%u, 1), ubound(s%u, 1)
      call init_beam(s%u(i))
    enddo
  else
    if (ix_universe < lbound(s%u, 1) .or. ix_universe > ubound(s%u, 1)) then
      call out_io (s_error$, r_name, &
            'BAD IX_UNIVERSE IN TAO_BEAM_INIT NAMELIST: \i0\ ', ix_universe)
      call err_exit
    endif
    call init_beam(s%u(ix_universe))
  endif

enddo

close (iu)

!----------------------------------------------------------------
!----------------------------------------------------------------
contains
!
! Initialize the beams. Determine which element to track beam to
!

subroutine init_beam (u)

use tao_read_beam_mod

implicit none

type (tao_universe_struct), target :: u
type (ele_pointer_struct), allocatable, save, target :: eles(:)
type (ele_struct), pointer :: ele
type (tao_universe_branch_struct), pointer :: uni_branch0
type (branch_struct), pointer :: branch

real(rp) v(6), bunch_charge, gamma
integer i, j, ix, iu, n_part, ix_class, n_bunch, n_particle
character(60) at, class, ele_name, line

! Set tracking start/stop

uni_branch0 => u%uni_branch(0)

uni_branch0%ix_track_start = ix_track_start
uni_branch0%ix_track_end   = ix_track_end

u%beam_info%beam_init = beam_init
u%beam_info%beam0_file = beam0_file
u%beam_info%beam_all_file = beam_all_file
u%design%lat_branch(0)%orbit(0)%vec = beam_init%center

! No initialization for a circular lattice

if (u%model%lat%param%lattice_type == circular_lattice$) return

! Find where to save the beam at.
! Always save at branch points.

do i = 0, ubound(u%model%lat%branch, 1)
  branch => u%design%lat%branch(i)
  u%uni_branch(i)%ele%save_beam = .false.
  u%uni_branch(i)%ele(0)%save_beam = .true.
  u%uni_branch(i)%ele(branch%n_ele_track)%save_beam = .true.
  do j = 1, ubound(branch%ele, 1)
    if (branch%ele(j)%key == branch$) u%uni_branch(i)%ele(j)%save_beam = .true.
  enddo
enddo

if (beam_saved_at /= '') then
  call tao_locate_elements (beam_saved_at, u%ix_uni, eles, err, .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD BEAM_SAVED_AT ELEMENT: ' // beam_saved_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      u%uni_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam = .true.
    enddo
  endif
endif

u%beam_saved_at = beam_saved_at

! If beam_all_file is set, read in the beam distributions.

if (u%beam_info%beam_all_file /= '') then
  tao_com%use_saved_beam_in_tracking = .true.
  call tao_open_beam_file (beam_all_file)
  call tao_read_beam_file_header (j, u%beam_info%beam_init%n_bunch, &
                                             u%beam_info%beam_init%n_particle, err)  
  if (err) call err_exit
  do
    if (j == -1) exit
    call tao_read_beam (uni_branch0%ele(j)%beam, err)
    if (err) call err_exit
  enddo  
  call out_io (s_info$, r_name, 'Read beam_all file: ' // u%beam_info%beam_all_file)
  call tao_close_beam_file ()
endif

if (u%connect%connected) u%connect%injecting_beam = u%current_beam
if (allocated(eles)) deallocate (eles)

end subroutine init_beam

end subroutine tao_init_beams

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!+
! Subroutine tao_init_data (data_file)
!
! Subroutine to initialize the tao data structures.
! If data_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   data_file -- Character(*): Tao data initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_data (data_file)

use tao_data_and_eval_mod
use tao_lattice_calc_mod
use tao_input_struct
use bmad_parser_mod
use random_mod
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data
type (tao_data_input) data(n_data_minn:n_data_maxx)
type (tao_datum_input) datum(n_data_minn:n_data_maxx) 

real(rp) default_weight        ! default merit function weight

integer ios, iu, i, j, j1, k, ix, n_uni, num
integer n, iostat
integer n_d1_data, ix_ele, ix_min_data, ix_max_data, ix_d1_data

integer :: n_d2_data(lbound(s%u, 1) : ubound(s%u, 1))

character(*) data_file
character(40) :: r_name = 'tao_init_data'
character(200) file_name
character(40) name,  universe, default_universe, d_typ
character(40) default_merit_type, default_attribute, data_type, default_data_source
character(40) use_same_lat_eles_as, source
character(100) line, default_data_type, search_for_lat_eles

logical err, free, gang
logical :: good_unis(lbound(s%u, 1) : ubound(s%u, 1))
logical :: mask(lbound(s%u, 1) : ubound(s%u, 1))

namelist / tao_d2_data / d2_data, n_d1_data, default_merit_type, universe

namelist / tao_d1_data / d1_data, data, datum, ix_d1_data, &
               default_weight, default_data_type, default_data_source, &
               use_same_lat_eles_as, search_for_lat_eles, ix_min_data, ix_max_data

!-----------------------------------------------------------------------
! Find out how many d2_data structures we need for each universe

call tao_hook_init_data() 
if (.not. tao_com%init_data .or. data_file == '') then
  call tao_init_data_end_stuff ()
  return
endif

!---

call tao_open_file ('TAO_INIT_DIR', data_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Data File: ' // file_name)

n_d2_data = 0

do 
  universe = '*'
  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'TAO_D2_DATA NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_d2_data)  ! force printing of error message
    enddo
  endif
  if (ios /= 0) exit

  if (universe == '*') then
    good_unis = .true.
  else
    call location_decode (universe, good_unis, lbound(s%u, 1), num)
    if (num < 0) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
  endif

  where (good_unis) n_d2_data = n_d2_data + 1

enddo

rewind (iu)

! Allocate space for the data

do i = lbound(s%u, 1), ubound(s%u, 1)
  call tao_init_data_in_universe (s%u(i), n_d2_data(i))
enddo

! Init data

do 
  mask(:) = .true.      ! set defaults
  d2_data%name       = ''
  universe           = '*'
  default_merit_type = 'target'

  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, 'Init: Read tao_d2_data namelist: ' // d2_data%name)

  if (universe == '*') then
    good_unis = .true.
  else
    call location_decode (universe, good_unis, lbound(s%u, 1), num)
  endif

  uni_loop: do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. good_unis(i)) cycle  
    ! check if this data type has already been defined for this universe
    do k = 1, size(s%u(i)%d2_data)
      if (trim(s%u(i)%d2_data(k)%name) == trim(d2_data%name)) then
        mask(i) = .false.
        cycle uni_loop
      endif
    enddo
      
    call tao_d2_data_stuffit (s%u(i), d2_data%name, n_d1_data)
  enddo uni_loop

  do k = 1, n_d1_data
    use_same_lat_eles_as   = ''
    search_for_lat_eles    = ''
    d1_data%name           = ''
    default_weight         = 0      ! set default
    default_data_type      = ''
    default_data_source    = 'lat'
    ix_min_data            = int_garbage$
    ix_max_data            = int_garbage$

    data(:)%data_type      = ''
    data(:)%merit_type     = default_merit_type 
    data(:)%merit_type     = ''
    data(:)%ele_name       = ''
    data(:)%ele0_name      = ''
    data(:)%meas           = real_garbage$  ! used to tag when %meas_value is set in file
    data(:)%weight         = 0.0
    data(:)%invalid_value  = 0.0
    data(:)%ix_bunch       = 0
    data(:)%data_source    = ''
    data(:)%good_user      = .true.

    datum(:)%data_type      = ''
    datum(:)%merit_type     = default_merit_type 
    datum(:)%merit_type     = ''
    datum(:)%ele_name       = ''
    datum(:)%ele_start_name = ''
    datum(:)%ele_ref_name   = ''
    datum(:)%meas           = real_garbage$  ! used to tag when %meas_value is set in file
    datum(:)%weight         = 0.0
    datum(:)%invalid_value  = 0.0
    datum(:)%ix_bunch       = 0
    datum(:)%data_source    = ''
    datum(:)%good_user      = .true.

    read (iu, nml = tao_d1_data, iostat = ios)
    if (ios > 0) then
      call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST READ ERROR.')
      rewind (iu)
      do
        read (iu, nml = tao_d1_data)  ! force printing of error message
      enddo
    endif

    ! Transfer data(:) to datum(:) if needed

    if (any(data%data_type /= '') .or. any(data%ele_name /= '')) then
      datum(:)%data_type      = data(:)%data_type
      datum(:)%merit_type     = data(:)%merit_type
      datum(:)%merit_type     = data(:)%merit_type
      datum(:)%ele_name       = data(:)%ele_name
      datum(:)%meas           = data(:)%meas
      datum(:)%weight         = data(:)%weight
      datum(:)%invalid_value  = data(:)%invalid_value
      datum(:)%ix_bunch       = data(:)%ix_bunch
      datum(:)%data_source    = data(:)%data_source
      datum(:)%good_user      = data(:)%good_user
      do i = lbound(datum, 1), ubound(datum, 1)
        if (datum(i)%data_type == '') cycle
        d_typ = datum(i)%data_type
        if (d_typ(1:2) == 'i5') datum(i)%data_type = 'rad_int.' // trim(d_typ) ! Convert old style
        if (d_typ(1:6) == 'floor.' .or. d_typ == 'momentum_compaction' .or. &
            d_typ(1:12) == 'periodic.tt.' .or. d_typ(1:5) == 'phase' .or. &
            d_typ(1:2) == 'r.' .or. d_typ(1:10) == 'rel_floor.' .or. &
            d_typ == 's_position' .or. d_typ(1:2) == 't.' .or. &
            d_typ(1:3) == 'tt.') then
          datum(:)%ele_ref_name   = data(:)%ele0_name          
          if (datum(i)%ele_ref_name == '' .and. datum(i)%ele_name /= '') &
                                                       datum(i)%ele_ref_name = 'BEGINNING'
        else
          datum(i)%ele_start_name = data(i)%ele0_name
        endif
        ! convert old style to new
        if (d_typ == 'lattice')       datum(i)%data_type = 'lat'
        if (d_typ == 'beam_tracking') datum(i)%data_source = 'beam'
      enddo
    endif

    ! Convert old format to new

    if (datum(0)%ele_name(1:7) == 'SEARCH:') then
      call string_trim(datum(0)%ele_name(8:), search_for_lat_eles, ix)
    elseif (datum(0)%ele_name(1:5) == 'SAME:') then
      call string_trim (datum(0)%ele_name(6:), use_same_lat_eles_as, ix)
    endif

    ! Check that we read the correct namelist

    if (ix_d1_data /= k) then
      write (line, '(a, 2i4)') ', k, ix_d1_data'
      call out_io (s_abort$, r_name, &
                'ERROR: IX_D1_DATA MISMATCH FOR D2_DATA: ' // d2_data%name, &
                '       THE D1_DATA HAD THE NAME: ' // d1_data%name, &
                '       I EXPECTED IX_D1_DATA TO BE: \i3\ ', &
                '       I READ IX_D1_DATA TO BE: \i3\ ', &
                i_array = (/ k, ix_d1_data /) )  
      call err_exit
    endif
    do i = lbound(datum, 1), ubound(datum, 1)
      if ((datum(i)%ele_ref_name /= '' .or. datum(i)%ele_start_name /= '') .and. datum(i)%ele_name == '') then
        write (line, '(4a, i0, a)') trim(d2_data%name), '.', trim(d1_data%name), '[', i, ']'
        call out_io (s_abort$, r_name, &
              'ERROR: ELE_NAME IS BLANK BUT ELE_REF_NAME OR ELE_START_NAME IS NOT FOR: ' // line)
        call err_exit
      endif
    enddo
    call out_io (s_blank$, r_name, 'Init: Read tao_d1_data namelist: ' // d1_data%name)

    do i = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. good_unis(i)) cycle
      if (.not. mask(i)) cycle
      call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used)
    enddo

  enddo

enddo

close (iu)

! Init ix_data array

call tao_init_data_end_stuff ()

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine d1_data_stuffit (i_d1, u, n_d2)

type (tao_universe_struct), target :: u
type (tao_d1_data_struct), pointer :: d1_this
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (ele_pointer_struct), allocatable, save :: eles(:)
type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: dat

integer i, n1, n2, ix, k, ix1, ix2, j, jj, n_d2

integer i_d1

character(20) fmt

!

d1_this => u%d2_data(n_d2)%d1(i_d1)  
if (d1_data%name == '') then
  write (d1_this%name, '(i0)') i_d1
else
  d1_this%name = d1_data%name    ! stuff in the data
endif

! Check if we are searching for elements or repeating elements
! and record the element names in the data structs.
    
if (search_for_lat_eles /= '') then
  call tao_init_find_elements (u, search_for_lat_eles, eles)
  if (size(eles) == 0) then
    call out_io (s_warn$, r_name, &
      'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_for_lat_eles, &
      'WHILE SETTING UP DATA ARRAY: ' // tao_d2_d1_name(d1_this))
    return
  endif
  ! finish finding data array limits
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(eles)
  call tao_allocate_data_array (u, n2)

  if (ix_min_data == int_garbage$) ix_min_data = 1
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)

  ! get element names
  jj = n1
  do k = lbound(eles, 1), ubound(eles, 1)
    if (jj .gt. n2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    u%data(jj)%ele_name  = eles(k)%ele%name
    u%data(jj)%ix_ele    = eles(k)%ele%ix_ele
    u%data(jj)%ix_branch = eles(k)%ele%ix_branch
    u%data(jj)%exists    = .true.
    jj = jj + 1
  enddo

  u%data(n1:n2)%good_user      = datum(ix1:ix2)%good_user
  u%data(n1:n2)%weight         = datum(ix1:ix2)%weight
  u%data(n1:n2)%invalid_value  = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%ele_start_name = datum(ix1:ix2)%ele_start_name
  u%data(n1:n2)%ele_ref_name   = datum(ix1:ix2)%ele_ref_name
  u%data(n1:n2)%ix_bunch       = datum(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source    = datum(ix1:ix2)%data_source

! use_same_lat_eles_as

elseif (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_data (err, name, d1_array = d1_array, ix_uni = u%ix_uni)
  if (err .or. size(d1_array) /= 1) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(d1_array(1)%d1%d)
  call tao_allocate_data_array (u, n2)

  ix_min_data = lbound(d1_array(1)%d1%d, 1)
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)

  u%data(n1:n2)%ele_name        = d1_array(1)%d1%d%ele_name
  u%data(n1:n2)%ix_ele          = d1_array(1)%d1%d%ix_ele
  u%data(n1:n2)%ele_ref_name    = d1_array(1)%d1%d%ele_ref_name
  u%data(n1:n2)%ix_ele_ref      = d1_array(1)%d1%d%ix_ele_ref
  u%data(n1:n2)%ele_start_name  = d1_array(1)%d1%d%ele_start_name
  u%data(n1:n2)%ix_ele_start    = d1_array(1)%d1%d%ix_ele_start
  u%data(n1:n2)%exists          = d1_array(1)%d1%d%exists
  u%data(n1:n2)%data_source     = d1_array(1)%d1%d%data_source
  u%data(n1:n2)%invalid_value   = d1_array(1)%d1%d%invalid_value

! Not SEARCH or SAME:

else

  if (ix_min_data == int_garbage$) ix_min_data = 1
  if (ix_max_data == int_garbage$) then
    do i = ubound(datum, 1), lbound(datum, 1), -1
      if (datum(i)%ele_name /= '' .or. datum(i)%data_type /= '') then
        ix_max_data = i
        exit
      endif
    enddo
  endif

  if (ix_max_data == int_garbage$) then
    call out_io (s_error$, r_name, 'NO DATA FOUND FOR: ' // tao_d2_d1_name(d1_this))
    return
  endif

  n1 = u%n_data_used + 1
  n2 = u%n_data_used + ix_max_data - ix_min_data + 1
  ix1 = ix_min_data
  ix2 = ix_max_data
  call tao_allocate_data_array (u, n2)

  ! Transfer info from the input structure

  u%data(n1:n2)%good_user      = datum(ix1:ix2)%good_user
  u%data(n1:n2)%weight         = datum(ix1:ix2)%weight
  u%data(n1:n2)%ele_name       = datum(ix1:ix2)%ele_name
  u%data(n1:n2)%invalid_value  = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%ele_ref_name   = datum(ix1:ix2)%ele_ref_name
  u%data(n1:n2)%ele_start_name = datum(ix1:ix2)%ele_start_name
  u%data(n1:n2)%ix_bunch       = datum(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source    = datum(ix1:ix2)%data_source

  ! Find elements associated with the data

  do j = n1, n2

    if (u%data(j)%exists) cycle

    if (u%data(j)%ele_name == '') cycle
    call str_upcase (u%data(j)%ele_name, u%data(j)%ele_name)
    call element_locator (u%data(j)%ele_name, u%design%lat, ix)
    if (ix < 0) then
      call out_io (s_error$, r_name, 'ELEMENT NOT LOCATED: ' // u%data(j)%ele_name)
      u%data(j)%exists = .false.
      cycle
    endif

    u%data(j)%ix_ele = ix
    u%data(j)%exists = .true.

    if (u%data(j)%ele_ref_name /= '') then
      call str_upcase (u%data(j)%ele_ref_name, u%data(j)%ele_ref_name)
      call element_locator (u%data(j)%ele_ref_name, u%design%lat, ix)
      if (ix < 0) then
        call out_io (s_error$, r_name, 'ELE_REF NOT LOCATED: ' // u%data(j)%ele_ref_name)
        u%data(j)%exists = .false.
        cycle
      endif
      u%data(j)%ix_ele_ref = ix
    endif

    if (u%data(j)%ele_start_name /= '') then
      call str_upcase (u%data(j)%ele_start_name, u%data(j)%ele_start_name)
      call element_locator (u%data(j)%ele_start_name, u%design%lat, ix)
      if (ix < 0) then
        call out_io (s_error$, r_name, 'ELE_START NOT LOCATED: ' // u%data(j)%ele_start_name)
        u%data(j)%exists = .false.
        cycle
      endif
      u%data(j)%ix_ele_start = ix
    endif

  enddo

endif

!-----------------------------------------------------------
! If %meas_value was set then %good_meas is set to True

u%data(n1:n2)%data_type     = datum(ix1:ix2)%data_type
u%data(n1:n2)%merit_type    = datum(ix1:ix2)%merit_type
u%data(n1:n2)%weight        = datum(ix1:ix2)%weight
u%data(n1:n2)%invalid_value = datum(ix1:ix2)%invalid_value

u%data(n1:n2)%meas_value = datum(ix1:ix2)%meas
where (u%data(n1:n2)%meas_value == real_garbage$)  ! where %meas_value was set
  u%data(n1:n2)%meas_value = 0  
elsewhere
  u%data(n1:n2)%good_meas = .true.
end where

!--------------------------------------------
! use default_data_type if given, if not, auto-generate the data_type

if (default_data_type == '') then
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = &
                            trim(d2_data%name) // '.' // d1_data%name
else
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = default_data_type
endif

! Point the %data back to the d1_data_struct

call tao_point_d1_to_data (d1_this, u%data(n1:n2), ix_min_data)

! In a d1_data array, not all the datums need to exist. 
! If a datum is not associated with an element, that generally means that
! it does not exist. There are, however, a few exceptions. EG: unstable.ring, etc.
! Here we mark data%exists for such datums.
! Also determine if we need to do the radiation integrals. This can save a lot of time.

do j = n1, n2
  dat => u%data(j)
  if (dat%weight == 0) dat%weight = default_weight
  if (dat%merit_type == '') dat%merit_type = default_merit_type
  if (dat%data_source == '') dat%data_source = default_data_source

  ! Convert old style to new style

  ix = index(dat%data_type, 'emittance.')
  if (ix /= 0) dat%data_type = dat%data_type(1:ix-1) // 'emit.' // dat%data_type(ix+10:)
  if (dat%data_type(1:9) == 'unstable_') dat%data_type(9:9) = '.'

  !
  data_type = dat%data_type
  source = dat%data_source

  if (tao_rad_int_calc_needed(data_type, source)) then
    u%do_rad_int_calc = .true. 
    if (dat%ix_branch /= 0) then
      call out_io (s_fatal$, r_name, 'EVALUATING A DATUM OF TYPE: ' // data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
      call err_exit
    endif
  endif

  if (data_type(1:6) == 'chrom.') u%do_chrom_calc = .true.

  if (data_type == 'unstable.orbit' .or. &
              (u%design%lat%param%lattice_type == circular_lattice$ .and. &
              (data_type(1:6)  == 'chrom.' .or. data_type(1:17) == 'multi_turn_orbit.' .or. &
               data_type(1:13) == 'unstable.ring' .or. index(data_type, 'emit.') /= 0))) then
    dat%exists = .true.
    if (dat%ele_name /= '') then
      call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // data_type, &
                        'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // dat%ele_name)
      call err_exit
    endif
  endif

enddo

if (.not. any(u%data(n1:n2)%exists)) then
  call out_io (s_warn$, r_name, &
            'Note: All datums in: ' // tao_d2_d1_name(d1_this), &
            'are marked as non-existent')
endif

end subroutine d1_data_stuffit

end subroutine tao_init_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Defines what datums to evaluate at each element in specified universe

subroutine tao_init_data_end_stuff ()

implicit none

type (tao_universe_struct), pointer :: u
type (tao_data_struct), pointer :: data
type (tao_element_struct), pointer :: uni_ele(:)
integer i, ib, j, k, ix_ele, n_max

logical err

!

do i = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i)
  call tao_allocate_data_array (u, u%n_data_used, .true.) ! Trim u%data size

  do ib = 0, ubound(u%model%lat%branch, 1)
    uni_ele => u%uni_branch(ib)%ele
    uni_ele(:)%n_datum = 0
  end do

  ! Since some information gets lost during tracking (like beam distributions),
  !   find where each datum gets evaluated when tao_load_data_array is called.
  ! ix_ele = -1  -->  Gets evaluated after all tracking

  do j = 1, size(u%data)
    data => u%data(j)
    if (.not. data%exists) cycle
    if (data%data_type(1:17) == 'multi_turn_orbit.') then
      cycle ! Does not get evaluated by tao_lattice_calc_mod
    elseif (data%data_type(1:8) == 'rad_int.') then
      ix_ele = -1
    elseif (data%ix_ele > s%u(data%d1%d2%ix_uni)%model%lat%n_ele_track) then
      ix_ele = -1
    elseif (data%ix_ele == -1) then
      ix_ele = -1
    elseif (index(data%data_type, 'emit.') /= 0 .and. data%data_source == 'lat') then
      ix_ele = -1
    elseif (data%ix_ele_ref > data%ix_ele) then
      ix_ele = u%model%lat%n_ele_track
    else
      ix_ele = data%ix_ele
    endif

    uni_ele => u%uni_branch(data%ix_branch)%ele
    uni_ele(ix_ele)%n_datum = uni_ele(ix_ele)%n_datum + 1 
  enddo
    
  ! allocate ix_datum array for each element

  do ib = 0, ubound(u%model%lat%branch, 1)
    uni_ele => u%uni_branch(ib)%ele
    do j = -1, ubound(uni_ele, 1)
      if (uni_ele(j)%n_datum == 0) cycle
      allocate (uni_ele(j)%ix_datum(uni_ele(j)%n_datum))
    enddo
    uni_ele(:)%n_datum = 0
  end do

  ! setup ix_ele array for each element
  ! This is the point where the datum is evaluated
  ! if ix_ele_ref > ix_ele then there is "wrap around"

  do j = 1, size(u%data)
    data => u%data(j)
    if (.not. data%exists) cycle
    if (data%data_type(1:17) == 'multi_turn_orbit.') then
      cycle ! Does not get evaluated by tao_lattice_calc_mod
    elseif (data%data_type(1:8) == 'rad_int.') then
      ix_ele = -1
    elseif (data%ix_ele > s%u(data%d1%d2%ix_uni)%model%lat%n_ele_track) then
      ix_ele = -1
    elseif (data%ix_ele == -1) then
      ix_ele = -1
    elseif (index(data%data_type, 'emit.') /= 0 .and. data%data_source == 'lat') then
      ix_ele = -1
    elseif (data%ix_ele_ref > data%ix_ele) then
      ix_ele = u%model%lat%n_ele_track
    else
      ix_ele = data%ix_ele
    endif

    uni_ele => u%uni_branch(data%ix_branch)%ele
    k = uni_ele(ix_ele)%n_datum + 1
    uni_ele(ix_ele)%ix_datum(k) = j
    uni_ele(ix_ele)%n_datum = k
  enddo

enddo

call tao_data_check (err)
if (err) stop

end subroutine tao_init_data_end_stuff

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_allocate_data_array (u, n_data, exact)

type (tao_universe_struct) :: u
type (tao_data_struct), allocatable :: data(:)
type (tao_d1_data_struct), pointer :: d1

integer i, j1, j2, n0, n_data
logical, optional :: exact  ! Default = False

! Exact means that size(u%data) must end up to be n_data.
! Not exact means that size(u%data) must be at least n_data.

u%n_data_used = n_data
  
if (n_data == size(u%data)) return 
if (.not. logic_option(.false., exact) .and. n_data < size(u%data)) return 

! Reallocate the data array. 
! If not exact then allocate more space than needed to reduce the number
! of times we need to reallocate stuff.

if (allocated(u%data)) then
  n0 = min(n_data, size(u%data))
  allocate (data(n0))
  data = u%data(1:n0)
  deallocate (u%data)
  if (logic_option(.false., exact)) then
    allocate (u%data(n_data))
  else
    allocate (u%data(2*n_data))
  endif
  u%data(1:n0) = data
  deallocate (data)
else
  allocate(u%data(n_data))
endif

! Since the data array gets reallocated the pointer from d1 to the datums must 
! be reestablished.

j2 = 0
do
  j1 = j2 + 1
  if (j1 > n0) exit
  d1 => u%data(j1)%d1
  if (.not. associated(d1)) exit
  do 
    if (j2 == n0) exit
    if (.not. associated(u%data(j2+1)%d1, d1)) exit
    j2 = j2 + 1
  enddo
  call tao_point_d1_to_data (d1, u%data(j1:j2), u%data(j1)%ix_d1)
enddo

! set defaults

do i = n0+1, size(u%data)
  u%data(i)%exists         = .false.       ! set default
  u%data(i)%good_meas      = .false.   ! set default
  u%data(i)%good_ref       = .false.   ! set default
  u%data(i)%good_user      = .true.    ! set default
  u%data(i)%good_opt       = .true.
  u%data(i)%merit_type     = 'target'  ! set default
  u%data(i)%ele_name       = ''
  u%data(i)%ix_ele         = -1
  u%data(i)%ele_ref_name   = ''
  u%data(i)%ix_ele_ref     = -1 
  u%data(i)%ele_start_name = ''
  u%data(i)%ix_ele_start   = -1 
  u%data(i)%ix_data        = i
  u%data(i)%ix_branch      = 0
enddo

end subroutine tao_allocate_data_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_d2_data_stuffit (u, d2_name, n_d1_data)

type (tao_universe_struct), target :: u
type (tao_d2_data_struct), pointer :: d2

integer i, nn, n_d1_data
character(*) d2_name
character(40) :: r_name = 'tao_d2_data_stuffit'

! Setup another d2_data structure.

u%n_d2_data_used = u%n_d2_data_used + 1
nn = u%n_d2_data_used

if (size(u%d2_data) < nn) then
  call out_io (s_error$, r_name, 'D2_DATA ARRAY OVERFLOW!')
  call err_exit
endif

d2 => u%d2_data(nn)

d2%name = d2_name
d2%ix_uni = u%ix_uni

! allocate memory for the u%d1_data structures

if (allocated(d2%d1)) deallocate (d2%d1)
allocate(d2%d1(n_d1_data))

do i = 1, n_d1_data
  d2%d1(i)%d2 => d2
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine tao_init_data_in_universe (u, n_d2_data)

implicit none

type (tao_universe_struct) u
integer j, n_d2_data

!

allocate (u%data(0))
u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
u%n_data_used = 0         ! size of s%u(i)%data(:) array
u%ix_rad_int_cache = 0

if (n_d2_data == 0) return
if (allocated(u%d2_data)) deallocate (u%d2_data)

allocate (u%d2_data(n_d2_data))

do j = 1, n_d2_data
  u%d2_data(j)%descrip = ''
enddo

u%d2_data%name = ''  ! blank name means it doesn't exist

! This is needed to keep the totalview debugger happy.

if (allocated(u%dmodel_dvar)) deallocate (u%dmodel_dvar)
allocate (u%dmodel_dvar(1,1))

end subroutine tao_init_data_in_universe

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_d1_to_data (d1, data, n_min)
!
! Routine used for arbitrary data pointer indexing
!
! d1     -- tao_data_struct: the pointer
! data   -- tao_data_struct: the data
! n_min  -- integer: starting index for the pointer
!-

subroutine tao_point_d1_to_data (d1, data, n_min)

implicit none

integer n, n_min, i, n0, n1

type (tao_d1_data_struct), target :: d1
type (tao_data_struct), target :: data(n_min:)

d1%d => data

do n = lbound(data, 1), ubound(data, 1)
  data(n)%d1 => d1
  data(n)%ix_d1 = n
enddo

end subroutine tao_point_d1_to_data

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!+
! Subroutine tao_init_variables (var_file)
!
! Subroutine to initialize the tao variable structures.
! If var_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   var_file  -- Character(*): Tao variable initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_variables (var_file)

use tao_lattice_calc_mod
use tao_input_struct
use bmad_parser_mod
use random_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_v1_var_input) v1_var
type (tao_var_struct), pointer :: var_ptr
type (tao_v1_var_struct), pointer :: v1_var_ptr
type (tao_this_var_struct), pointer :: this
type (tao_var_input) var(n_var_minn:n_var_maxx)

real(rp) default_weight        ! default merit function weight
real(rp) default_step          ! default "small" step size
real(rp) default_low_lim, default_high_lim, default_key_delta
real(rp), allocatable, save :: default_key_d(:)

integer ios, iu, i, j, j1, j2, k, ix, num
integer n, iostat, n_list
integer ix_min_var, ix_max_var, ix_ele, n_v1, n_v1_var_max

character(*) var_file
character(40) :: r_name = 'tao_init_variables'
character(40) name, universe, default_universe
character(40) default_merit_type, default_attribute
character(40) use_same_lat_eles_as
character(8) default_key_bound
character(200) file_name
character(8), allocatable, save :: default_key_b(:)
character(100) line, search_for_lat_eles

logical err, free, gang
logical searching
logical, allocatable, save :: dflt_good_unis(:), good_unis(:)

namelist / tao_var / v1_var, var, default_weight, default_step, default_key_delta, &
                    ix_min_var, ix_max_var, default_universe, default_attribute, &
                    default_low_lim, default_high_lim, default_merit_type, &
                    use_same_lat_eles_as, search_for_lat_eles, default_key_bound

!-----------------------------------------------------------------------
! Init

allocate (s%var(1))
s%n_var_used = 0

! Call the hook routine and then return if the standard init is not wanted.

call tao_hook_init_var() 
if (.not. tao_com%init_var .or. var_file == '') then
  call tao_init_var_end_stuff ()
  return
endif

! Standard var init:
! First open the var init file.

call tao_open_file ('TAO_INIT_DIR', var_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Variable File: ' // file_name)

! Count how many v1_var definitions there are

allocate (default_key_b(100), default_key_d(100))
n = 0
do
  default_key_bound = ''
  default_key_delta = 0
  read (iu, nml = tao_var, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'TAO_VAR NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_var)  ! force printing of error message
    enddo
  endif
  if (ios < 0) exit
  n = n + 1
  if (n >= size(default_key_b)) then
    call re_allocate (default_key_b, 2*size(default_key_b))
    call re_allocate (default_key_d, 2*size(default_key_d))
  endif
  default_key_b(n) = default_key_bound
  default_key_d(n) = default_key_delta
enddo

call tao_allocate_v1_var (n)
n_list = n

! Now fill in all the information

rewind (iu)

allocate (dflt_good_unis(lbound(s%u,1):ubound(s%u, 1)), good_unis(lbound(s%u,1):ubound(s%u,1)))

n_v1 = 0
do
  n_v1 = n_v1 + 1
  if (n_v1 > n_list) exit

  use_same_lat_eles_as = ''
  search_for_lat_eles  = ''
  v1_var%name        = " "         ! set default
  default_merit_type = 'limit'
  default_weight     = 0     ! set default
  default_step       = 0       ! set default
  default_attribute  = ''
  default_universe   = ''
  default_low_lim    = -1e30
  default_high_lim   = 1e30
  ix_min_var         = 1
  ix_max_var         = 0
  var%ele_name       = ''
  var%merit_type     = ''
  var%weight         = 0         ! set default
  var%step           = 0         ! set default
  var%attribute      = ''
  var%universe       = ''
  var%low_lim        = default_low_lim
  var%high_lim       = default_high_lim
  var%good_user      = ''
  var%key_bound      = default_key_b(n_v1)
  var%key_delta      = default_key_d(n_v1)

  read (iu, nml = tao_var, iostat = ios)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, 'Init: Read tao_var namelist: ' // v1_var%name)

  ! Convert old format to new

  call str_upcase (default_attribute, default_attribute)
  if (var(0)%ele_name(1:7) == 'SEARCH:') then
    call string_trim(var(0)%ele_name(8:), search_for_lat_eles, ix)
  elseif (var(0)%ele_name(1:5) == 'SAME:') then
    call string_trim (var(0)%ele_name(6:), use_same_lat_eles_as, ix)
  endif

  ! Convert to upper case.

  do i = lbound(var, 1), ubound(var, 1)
    call str_upcase (var(i)%attribute, var(i)%attribute)
    call str_upcase (var(i)%ele_name, var(i)%ele_name)
  enddo

  if (v1_var%name == '') then
    call out_io (s_error$, r_name, 'FOUND TAO_VAR NAMELIST WITHOUT V1_VAR%NAME PARAMETER!')
    cycle
  endif

  ! Gang or clone?

  gang = .true.
  if (default_universe(1:5) == 'clone') then
    gang = .false.
    call string_trim (default_universe(6:), default_universe, ix)
  endif
  if (default_universe(1:4) == 'gang') then
    call string_trim (default_universe(5:), default_universe, ix)
  endif

  ! Read universe numbers

  if (default_universe == '*' .or. default_universe == '') then
    dflt_good_unis = .true.
    if (tao_com%common_lattice .and. gang) then
      dflt_good_unis = .false.
      dflt_good_unis(ix_common_uni$) = .true.
    endif

  else
    call location_decode (default_universe, dflt_good_unis, 1, num)
    if (num == 0) dflt_good_unis = .true.  ! blank => all
    if (num < 0) then
      call out_io (s_error$, r_name, 'ERROR READING DEFAULT_UNIVERSE FOR: ' // v1_var%name)
      cycle
    endif
  endif

  ! Set up variable lists

  if (gang) then
    call tao_var_stuffit1 (var, v1_var_ptr, v1_var, -1, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

    if (.not. searching) then
      do j = lbound(v1_var_ptr%v, 1), ubound(v1_var_ptr%v, 1)

        ! Find which universes
        good_unis = dflt_good_unis

        if (.not. searching) then
          if (var(j)%universe /= '') then
            if (var(j)%universe == '*') then
              good_unis = .true.
            else
              call location_decode (var(j)%universe, good_unis, 1, num)
              if (num < 0) then
                call out_io (s_error$, r_name, 'ERROR READING UNIVERSE FOR: ' // v1_var%name)
                cycle
              endif
            endif
          endif
        endif

        if (count(good_unis) == 0) then
          call out_io (s_error$, r_name, 'ERROR: NO UNIVERSE FOR: ' // v1_var%name)
          call err_exit
        endif

        call tao_var_stuffit2 (good_unis, v1_var_ptr%v(j))
      enddo
    endif

  else   ! If clone...
    do i = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. dflt_good_unis(i)) cycle
      call tao_var_stuffit1 (var, v1_var_ptr, v1_var, i, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

      write (v1_var_ptr%name, '(2a, i0)') trim(v1_var_ptr%name), '_u', i
      if (.not. searching) then
        good_unis = .false.
        good_unis(i) = .true.
        do j = lbound(v1_var_ptr%v, 1), ubound(v1_var_ptr%v, 1)
          call tao_var_stuffit2 (good_unis, v1_var_ptr%v(j))
        enddo
      endif
    enddo
  endif

enddo

close (iu)
deallocate (dflt_good_unis, good_unis)
deallocate (default_key_b, default_key_d)

call tao_init_var_end_stuff ()

end subroutine tao_init_variables

!-----------------------------------------------------------------------
!------------------------------------------------------------------------
! stuff common to all universes

subroutine tao_var_stuffit1 (var, v1_var_ptr, v1_var, ix_uni, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

use tao_input_struct

implicit none

type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_var_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: var_ptr
type (ele_pointer_struct), allocatable, save :: eles(:)
type (tao_v1_var_input) v1_var
type (tao_var_input) var(n_var_minn:n_var_maxx)

real(rp) default_weight, default_step, default_low_lim, default_high_lim

character(*) use_same_lat_eles_as, search_for_lat_eles, default_attribute
character(*) default_merit_type
character(20) fmt
character(40) name
character(40), allocatable :: ele_names(:)
character(60) search_string
character(20) :: r_name = 'tao_var_stuffit1'

integer i, iu, ip, j, jj, k, kk, n, nn, n1, n2, ix1, ix2, ix
integer num_ele, ios, ix_uni, ixm, ix2m
integer, allocatable, save :: an_indexx(:)
integer ix_min_var, ix_max_var

logical :: dflt_good_unis(lbound(s%u,1):)
logical searching, grouping, found_one, bound, err

! count number of v1 entries

s%n_v1_var_used = s%n_v1_var_used + 1
nn = s%n_v1_var_used
v1_var_ptr => s%v1_var(nn)
v1_var_ptr%name = v1_var%name

! If reusing a previous element list...

if (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_var (err, name, v1_array = v1_array)
  if (err .or. size(v1_array) /= 1) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "USE_SAME_LAT_ELES_AS": ' // name)
    call err_exit
  endif
  v1_ptr => v1_array(1)%v1
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + size(v1_ptr%v)
  call tao_allocate_var_array (n2)

  ix_min_var = lbound(v1_ptr%v, 1)
  ix1 = ix_min_var
  ix2 = ix1 + (n2 - n1)

  s%var(n1:n2)%ele_name    = v1_ptr%v%ele_name
  s%var(n1:n2)%s           = v1_ptr%v%s

  do n = n1, n2
    ix = ix1 + (n - n1)
    ip = 1 + (n - n1) 

    
    s%var(n)%good_user = v1_ptr%v(ip)%good_user
    if (var(ix)%good_user /= '') read (var(ix)%good_user, *, iostat = ios) s%var(n)%good_user
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "GOOD_USER" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif

    s%var(n)%ix_key_table = v1_ptr%v(ip)%ix_key_table
    if (var(ix)%key_bound /= '') read (var(ix)%key_bound, *, iostat = ios) bound
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "KEY_BOUND" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif
    if (bound) s%var(n)%ix_key_table = 1

    s%var(n)%key_delta = v1_ptr%v(ip)%key_delta
    if (var(ix)%key_delta /= 0) s%var(n)%key_delta = var(ix)%key_delta

    s%var(n)%attrib_name = v1_ptr%v(ip)%attrib_name
    if (default_attribute /= '') s%var(n)%attrib_name = default_attribute
    if (var(ix)%attribute /= '') s%var(n)%attrib_name = var(ix)%attribute

    s%var(n)%weight = v1_ptr%v(ip)%weight
    if (default_weight /= 0) s%var(n)%weight = default_weight
    if (var(ix)%weight /= 0) s%var(n)%weight = var(ix)%weight

    s%var(n)%step = v1_ptr%v(ip)%step
    if (default_step /= 0) s%var(n)%step = default_step
    if (var(ix)%step /= 0) s%var(n)%step = var(ix)%step

    s%var(n)%merit_type = v1_ptr%v(ip)%merit_type
    if (default_merit_type /= '') s%var(n)%merit_type = default_merit_type
    if (var(ix)%merit_type /= '') s%var(n)%merit_type = var(ix)%merit_type

    s%var(n)%low_lim = v1_ptr%v(ip)%low_lim
    if (default_low_lim /= -1e30) s%var(n)%low_lim = default_low_lim
    if (var(ix)%low_lim /= -1e30) s%var(n)%low_lim = var(ix)%low_lim

    s%var(n)%high_lim = v1_ptr%v(ip)%high_lim
    if (default_high_lim /= 1e30) s%var(n)%high_lim = default_high_lim
    if (var(ix)%high_lim /= 1e30) s%var(n)%high_lim = var(ix)%high_lim

    s%var(n)%ix_key_table = v1_ptr%v(ip)%ix_key_table

    s%var(n)%ix_v1 = ix_min_var + n - n1
    s%var(n)%v1 => v1_var_ptr
  enddo

  call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var)

  return

endif

!------------------------------
! Are we searching for elements?

if (search_for_lat_eles /= '') then

  searching = .true.
  call string_trim (search_for_lat_eles, search_string, ix)
  grouping = .true.  ! Default
  if (search_string(1:12) == '-no_grouping') then
    grouping = .false.
    search_string = search_string(13:)
  endif
  search_string = '-no_slaves ' // trim(search_string)
  if (any(var%universe /= '')) then
    call out_io (s_abort$, r_name, &
           "CANNOT SPECIFY INDIVIDUAL UNIVERSES WHEN SEARCHING FOR VARIABLES")
    call err_exit
  endif

  ! find matching elements...
  ! first count how many variables we need

  num_ele = 0
  allocate(ele_names(100), an_indexx(100))
  do iu = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. dflt_good_unis(iu)) cycle
    if (ix_uni > -1 .and. iu /= ix_uni) cycle
    call tao_init_find_elements (s%u(iu), search_string, eles, default_attribute)
    if (grouping) then
      do kk = 1, size(eles)
        call find_indexx(eles(kk)%ele%name, ele_names, an_indexx, num_ele, ixm, ix2m)
        if (ixm == 0) then
          if (num_ele+1 > size(ele_names)) then
            call re_allocate(ele_names, size(ele_names) + 100)
            call re_allocate(an_indexx, size(an_indexx) + 100)
          endif
          an_indexx(ix2m+1:num_ele+1) = an_indexx(ix2m:num_ele)
          an_indexx(ix2m) = num_ele + 1
          ele_names(num_ele+1) = eles(kk)%ele%name
          num_ele = num_ele + 1
        endif
      enddo
    else
      num_ele = num_ele + size(eles)
    endif
  enddo
  deallocate(ele_names, an_indexx)

  if (num_ele == 0) then
    call out_io (s_error$, r_name, &
              'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_for_lat_eles, &
              'WHILE SETTING UP VARIABLE ARRAY: ' // v1_var%name)
    return
  endif

  ! Now load the information into the variable structures

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + num_ele
  call tao_allocate_var_array (n2)

  nn = n1 - 1
  do iu = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. dflt_good_unis(iu)) cycle
    if (ix_uni > -1 .and. iu /= ix_uni) cycle
    call tao_init_find_elements (s%u(iu), search_string, eles, default_attribute)
    kk_loop: do kk = 1, size(eles)
      ! If the name matches an existing variable then use that variable
      if (grouping) then
        do j = n1, nn
          if (s%var(j)%ele_name == eles(kk)%ele%name) then
            call tao_pointer_to_var_in_lattice (s%var(j), iu, eles(kk)%ele, err)
            cycle kk_loop
          endif
        enddo
      endif
      ! Here if there is no match...
      ! create a new variable and stuff the info into it.
      nn = nn + 1
      var_ptr => s%var(nn)
      var_ptr%v1 => v1_var_ptr  ! Used for error messages in pointer_to_var
      var_ptr%ix_v1 = ix_min_var + nn - n1
      var_ptr%exists = .true.
      var_ptr%ele_name = eles(kk)%ele%name
      var_ptr%s = eles(kk)%ele%s
      var_ptr%attrib_name = default_attribute
      call tao_pointer_to_var_in_lattice (var_ptr, iu, eles(kk)%ele, err)
    enddo kk_loop
  enddo 

  ix1 = ix_min_var
  ix2 = (n2 - n1) + ix_min_var

!------------------------------
! If not searching or reusing...

else  
  searching = .false.

  ! if ix_max/min_var has not been set then just search for elements that have been named
  ! and use this info to set ix_max/min

  if (ix_max_var < ix_min_var) then
    found_one = .false.
    do i = lbound(var, 1), ubound(var, 1)
      if (var(i)%ele_name /= '') then
        ix_max_var = i
        if (.not. found_one) ix_min_var = i
        found_one = .true.
      endif
    enddo
    if (ix_min_var > 1) ix_min_var = 1 
  endif

  if (ix_max_var < ix_min_var) then
    call out_io (s_abort$, r_name, 'NO ELEMENTS FOR: ' // v1_var%name)
    call err_exit
  endif

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + ix_max_var - ix_min_var + 1
  ix1 = ix_min_var
  ix2 = ix_max_var
 
  call tao_allocate_var_array (n2)

  s%var(n1:n2)%ele_name = var(ix1:ix2)%ele_name

endif

!---------------------------------------------------------------
! Common init stuff

do n = n1, n2
  i = ix1 + n - n1

  s%var(n)%ix_key_table = -1
  if (var(i)%key_bound /= '') then
    read (var(i)%key_bound, *, iostat = ios) bound
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "KEY_BOUND" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif
    if (bound) s%var(n)%ix_key_table = 1
  endif

  if (var(i)%good_user == '') then
    s%var(n)%good_user = .true.
  else
    read (var(i)%good_user, *, iostat = ios) s%var(n)%good_user
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "GOOD_USER" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif
  endif

  s%var(n)%ix_v1 = ix_min_var + n - n1
  s%var(n)%v1 => v1_var_ptr

enddo

call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var)

s%var(n1:n2)%key_delta = var(ix1:ix2)%key_delta

s%var(n1:n2)%attrib_name = var(ix1:ix2)%attribute
where (s%var(n1:n2)%attrib_name == '') s%var(n1:n2)%attrib_name = default_attribute

s%var(n1:n2)%weight = var(ix1:ix2)%weight
where (s%var(n1:n2)%weight == 0) s%var(n1:n2)%weight = default_weight
 
s%var(n1:n2)%step = var(ix1:ix2)%step
where (s%var(n1:n2)%step == 0) s%var(n1:n2)%step = default_step
 
s%var(n1:n2)%merit_type = var(ix1:ix2)%merit_type
where (s%var(n1:n2)%merit_type == '') s%var(n1:n2)%merit_type = default_merit_type
 
s%var(n1:n2)%low_lim = var(ix1:ix2)%low_lim
where (s%var(n1:n2)%low_lim == -1e30) s%var(n1:n2)%low_lim = default_low_lim
 
s%var(n1:n2)%high_lim = var(ix1:ix2)%high_lim
where (s%var(n1:n2)%high_lim == 1e30) s%var(n1:n2)%high_lim = default_high_lim

end subroutine tao_var_stuffit1 
  
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine tao_allocate_v1_var (n_v1)

implicit none

integer i, n_v1

!

if (allocated(s%v1_var)) deallocate (s%v1_var)
allocate (s%v1_var(n_v1))

do i = 1, n_v1
  s%v1_var(i)%name = ''  ! blank name means it doesn't (yet) exist
  s%v1_var(i)%ix_v1 = i
enddo

s%n_v1_var_used = 0
s%n_var_used = 0

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine tao_var_stuffit2 (good_unis, var)

implicit none

type (tao_var_struct), target :: var
type (tao_this_var_struct), pointer :: this
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat

integer i, j, n, n1, n2, ie, iu, ib

character(20) :: r_name = 'tao_var_stuffit2'
logical err, good_unis(lbound(s%u, 1):), found

! 

if (var%ele_name == '') then
  var%exists = .false.
  var%ix_key_table = -1
  return
endif

found = .false.
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. good_unis(iu)) cycle
  lat => s%u(iu)%model%lat
  do ib = 0, ubound(lat%branch, 1)
    do ie = 0, lat%branch(ib)%n_ele_max
      if (var%ele_name /= lat%branch(ib)%ele(ie)%name) cycle
      call tao_pointer_to_var_in_lattice (var, iu, lat%branch(ib)%ele(ie), err)
      if (err) return
      found = .true.
    enddo
  enddo
enddo

if (.not. found) then
  call out_io (s_error$, r_name, &
            'CANNOT FIND LATTICE ELEMENT WITH NAME: "' // trim(var%ele_name) // '"', &
            'FOR VARIABLE: ' // var%v1%name)
  return
endif

if (size(var%this) > 0) then
  var%exists = .true.
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Put the variables marked by key_bound in the key table.

subroutine tao_init_var_end_stuff ()

implicit none

integer i, j

! Key table setup

allocate (s%key(count(s%var%ix_key_table > 0)))

j = 0
do i = 1, s%n_var_used
  if (s%var(i)%ix_key_table < 1) cycle
  j = j + 1
  s%key(j) = i
  s%var(i)%key_val0 = s%var(i)%model_value
  s%var(i)%ix_key_table = j
enddo

end subroutine

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
! This searches the lattice for the specified element and flags eles(:)

subroutine tao_init_find_elements (u, search_string, eles, attribute)

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)

character(*) search_string
character(*), optional :: attribute
character(80) string
character(40) ele_name, key_name_in
character(20) :: r_name = 'tao_init_find_elements'

integer key, found_key, ix_attrib, t
integer i, k, ix, ii, j, ix0, ix1, ix2, n_ele

logical no_slaves, no_lords, err, warn_given

! Sort switches

no_slaves = .false.
no_lords = .false.

call string_trim(search_string, string, ix)

do

  if (string(1:1) /= '-') exit

  select case (string(1:ix))
  case ('-no_lords') 
    no_lords = .true.
  case ('-no_slaves') 
    no_slaves = .true.
  case default
    call out_io (s_abort$, r_name, "BAD SEARCH SWITCH: " // search_string)
    call err_exit
  end select

  call string_trim (string(ix+1:), string, ix)

enddo

! Find elements

call tao_locate_elements (string, u%ix_uni, eles, err)

warn_given = .false.
n_ele = 0
do j = 1, size(eles)
  ele => eles(j)%ele
  t = ele%slave_status
  if ((t == multipass_slave$ .or. t == super_slave$) .and. no_slaves) cycle
  t = ele%lord_status 
  if ((t == girder_lord$ .or. t == overlay_lord$ .or. t == group_lord$ .or. &
       t == super_lord$ .or. t == multipass_lord$) .and. no_lords) cycle
  ! If attribute is not free then don't count it
  if (present(attribute)) then
    ix_attrib = attribute_index(ele, attribute)
    if (ix_attrib < 1) then
      call out_io (s_error$, r_name, &
                      'BAD ATTRIBUTE: ' // attribute, &
                      'FOR VARIABLE SEARCH: ' // search_string, &
                      'FOR ELEMENT: ' // ele%name)
      return
    endif
    if (.not. attribute_free (eles(j)%ele, attribute, u%model%lat, .false.)) then
      if (.not. warn_given) then
        call out_io (s_info$, r_name, &
                  'Non-free attribute ' // attribute, &
                  'For variable search: ' // search_string, &
                  'For element: ' // ele%name)
        warn_given = .true.
      endif
      cycle
    endif
  endif
  ! Passes test so add it to the list
  n_ele = n_ele + 1
  eles(n_ele)%ele => eles(j)%ele
enddo

call re_allocate_eles (eles, n_ele, .true., .true.)

end subroutine tao_init_find_elements

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pointer_to_var_in_lattice (var, ix_uni, ele, err)
! 
! Routine to set a pointer to the appropriate variable in a lattice
!
! Input:
!   var       -- Tao_var_struct: Structure has the info of where to point.
!   ix_uni    -- Integer: the universe to use
!   ix_ele    -- Integer: Index of element.
!
! Output:
!   var%this(ix_this) -- Tao_this_var_struct: 
!     %model_ptr
!     %base_ptr
!     %ix_ele
!     %ix_uni
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice (var, ix_uni, ele, err)

implicit none

type (tao_var_struct), target :: var
type (ele_struct) ele
type (ele_struct), pointer :: ele2
type (tao_universe_struct), pointer :: u
type (tao_this_var_struct), pointer :: this
type (tao_this_var_struct) :: this_saved(size(var%this))

integer ix, ix_uni, ix_this
logical :: err
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

err = .true.

u => s%u(ix_uni)
if (tao_com%common_lattice) u => tao_com%u_working

! allocate space for var%this.

ix_this = size(var%this) + 1
this_saved = var%this
deallocate (var%this)  
allocate (var%this(ix_this))
var%this(1:ix_this-1) = this_saved
this => var%this(ix_this)

! locate attribute

ele2 => pointer_to_ele (u%model%lat, ele%ix_ele, ele%ix_branch)
call pointer_to_attribute (ele2, var%attrib_name, .true., &
                          this%model_value, err, .false., ix_attrib = var%ix_attrib)
if (err) then
  call out_io (s_error$, r_name, &
            'IN VARIBALE: ' // tao_var1_name(var), &
            '  INVALID ATTRIBUTE: ' // var%attrib_name, &
            '  FOR ELEMENT: ' // ele%name)
  var%model_value => var%old_value
  var%base_value  => var%old_value
  var%exists = .false.
  var%ix_key_table = -1
  return
endif

ele2 => pointer_to_ele (u%base%lat, ele%ix_ele, ele%ix_branch)
call pointer_to_attribute (ele2,  var%attrib_name, .true., this%base_value,  err)

var%model_value => var%this(1)%model_value
var%base_value  => var%this(1)%base_value
var%design_value = var%this(1)%model_value

this%ix_ele    = ele%ix_ele
this%ix_branch = ele%ix_branch
this%ix_uni    = ix_uni

! Common pointer

if (associated(u%common)) then
  ele2 => pointer_to_ele (u%common%model%lat, ele%ix_ele, ele%ix_branch)
  call pointer_to_attribute (ele,  var%attrib_name, .true., var%common%model_value,  err)
  ele2 => pointer_to_ele (u%common%base%lat, ele%ix_ele, ele%ix_branch)
  call pointer_to_attribute (ele, var%attrib_name, .true., var%common%base_value,  err)
endif

! With unified lattices: model_value and base_value get their own storage
! instead of pointing to var%this(1). 
! Exception: If variable controls a common parameter

if (tao_com%common_lattice) then
  if (this%ix_uni == ix_common_uni$) then
    var%model_value => var%common%model_value
    var%base_value => var%common%base_value
  else
    allocate (var%model_value, var%base_value)
    var%model_value = this%model_value
    var%base_value = this%base_value
  endif
endif

end subroutine tao_pointer_to_var_in_lattice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_allocate_var_array (n_var)
!
! Routine to increase the s%var(:) array size.
!
! Input:
!   n_var -- Integer: Size of s%var(:) wanted.
!-

subroutine tao_allocate_var_array (n_var)

type (tao_var_struct) :: var(size(s%var))
type (tao_v1_var_struct), pointer :: v1

integer i, j1, j2, n0, n_var

! First save the information presently in s%var(:) in the var(:) array.
! s%n_var_used gives the present upper bound to s%var(:).

if (allocated(s%var)) then
  n0 = s%n_var_used
  do i = 1, n0
    allocate(var(i)%this(size(s%var(i)%this)))
  enddo
  var(1:n0) = s%var(1:n0)
  deallocate (s%var)
  allocate (s%var(n_var))
  do i = 1, n0
    allocate(s%var(i)%this(size(var(i)%this)))
  enddo
  s%var(1:n0) = var(1:n0)
else
  n0 = 0
  allocate (s%var(n_var))
endif

! Since the s%var(:) array has been reallocated, the pointer from s%var(:)%v1 
! to the datums must be reestablished.

j2 = 0
do
  j1 = j2 + 1
  if (j1 > n0) exit
  v1 => s%var(j1)%v1
  do 
    if (j2 == n0) exit
    if (.not. associated(s%var(j2+1)%v1, v1)) exit
    j2 = j2 + 1
  enddo
  call tao_point_v1_to_var (v1, s%var(j1:j2), s%var(j1)%ix_v1)
enddo

! Set the defaults for the newly created slots in the s%var(:) array

do i = n0+1, size(s%var)
  s%var(i)%ix_var    = i
  s%var(i)%good_opt  = .true.
  s%var(i)%exists    = .false.
  s%var(i)%good_var  = .true.
  s%var(i)%good_user = .true.
  s%var(i)%model_value => tao_com%dummy_target  ! Just to point to somewhere
  s%var(i)%base_value  => tao_com%dummy_target  ! Just to point to somewhere
  s%var(i)%low_lim = -1e30
  s%var(i)%high_lim = 1e30
  allocate(s%var(i)%this(0))
enddo
  
s%n_var_used = n_var


end subroutine tao_allocate_var_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_v1_to_var (v1, var, n)
!
! used for arbitrary variable pointer indexing
!
! v1       -- tao_v1_var_struct: Contains the pointer.
! var(n:)  -- tao_var_struct: the variable
! n        -- integer: starting index for the var array.
!-

subroutine tao_point_v1_to_var (v1, var, n)

implicit none

integer n, i, n_var

type (tao_v1_var_struct), target :: v1
type (tao_var_struct), target :: var(n:)

v1%v => var

do i = lbound(var, 1), ubound(var, 1)
  var(i)%v1 => v1
  var(i)%ix_v1 = i
enddo

end subroutine 

end module
