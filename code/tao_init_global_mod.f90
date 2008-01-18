module tao_init_global_mod

use tao_mod

contains

!+
! Subroutine tao_init_global (init_file)
!
! Subroutine to initialize the tao global structures.
! If init_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   init_file      -- Character(*): Tao initialization file.
!-

subroutine tao_init_global (init_file)

use tao_lattice_calc_mod
use tao_input_struct
use macroparticle_mod
use bmad_parser_mod
use random_mod
use csr_mod, only: csr_param
use spin_mod
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_global_struct), save :: global, default_global
type (beam_init_struct) beam_init
type (macro_init_struct) macro_init
type (tao_connected_uni_input) connect
type (spin_polar_struct) spin

integer ios, iu, i, j, i2, j2, k, ix, n_uni, num
integer n_data_max, n_var_max, ix_track_start, ix_track_end
integer n_d2_data_max, n_v1_var_max ! Deprecated variables
integer n, n_universes, iostat, ix_universe, n_max

character(*) init_file
character(40) :: r_name = 'tao_init_global'
character(200) file_name, beam0_file, beam_all_file
character(40) name,  universe

character(60) save_beam_at(100)
character(100) line

logical err

namelist / tao_params / global, bmad_com, csr_param, &
          n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  
namelist / tao_connected_uni_init / ix_universe, connect
  
namelist / tao_beam_init / ix_universe, beam0_file, &
  ix_track_start, ix_track_end, beam_all_file, beam_init, save_beam_at
         
namelist / tao_macro_init / ix_universe, macro_init

!-----------------------------------------------------------------------
! Init lattaces
! read global structure from tao_params namelist

call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Init File: ' // file_name)
if (iu == 0) then
  call out_io (s_abort$, r_name, "Error opening init file")
  call err_exit
endif

global = default_global         ! establish defaults
global%valid_plot_who(1:5) = (/ 'model ', 'base  ', 'ref   ', 'design', 'meas  ' /)
global%default_key_merit_type = 'limit'
n_var_max = 0
n_data_max = 0

call out_io (s_blank$, r_name, 'Init: Reading tao_params namelist')
read (iu, nml = tao_params, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_PARAMS NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_params)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_params namelist found')

s%global = global  ! transfer global to s%global
  
if (s%global%track_type == "macro") then
  call out_io (s_error$, r_name, &
             'MACROPARTICLE TRACKING IS NOT ACTIVELY SUPPORTED!', &
             'PLEASE USE "beam" TRACKING INSTEAD.', &
             'IF YOU NEED MACROPARTICLE TRACKING PLEASE SEE DAVID SAGAN.')
endif

n = size(s%u)
n_max = 0
do i = 1, size(s%u)
  s%u(i)%ix_uni = i
  s%u(i)%do_synch_rad_int_calc = .false.
  s%u(i)%do_chrom_calc         = .false.
  n_max = max(n_max, s%u(i)%design%lat%n_ele_max)
enddo

tao_com%n_data_max = n_data_max
tao_com%n_var_max = n_var_max

!-----------------------------------------------------------------------
! Seed random number generator

call ran_seed_put (s%global%random_seed)

select case (s%global%random_engine)
case ('pseudo')
  call ran_engine (pseudo_random$)
case ('quasi')
  call ran_engine (quasi_random$)
case default
  call out_io (s_error$, r_name, &
        'BAD GLOBAL%RANDOM_ENGINE SWITCH: ' // s%global%random_engine)
end select

select case (s%global%random_gauss_converter)
case ('exact')
  call ran_gauss_converter (exact_gaussian$)
case ('limited')
  call ran_gauss_converter (limited_gaussian$, sigma_cut = s%global%random_gauss_cutoff)
case default
  call out_io (s_error$, r_name, &
        'BAD GLOBAL%RANDOM_GAUSS_CONVERTER SWITCH: ' // s%global%random_gauss_converter)
end select

!-----------------------------------------------------------------------
! allocate lattice coord_structs and equate model and base to design

call init_orbits ()
  
! set model/base = design

do i = 1, size(s%u)
  s%u(i)%model%lat = s%u(i)%design%lat
  s%u(i)%base%lat  = s%u(i)%design%lat
enddo

!-----------------------------------------------------------------------
! Init connected universes

! defaults
do i = 1, size(s%u)
  s%u(i)%connect%connected = .false.
  s%u(i)%connect%match_to_design = .false.
  s%u(i)%connect%use_connect_ele = .false.
  s%u(i)%connect%from_uni = -1
  s%u(i)%connect%from_uni_s = -1
  s%u(i)%connect%from_uni_ix_ele = -1
enddo

do
  ix_universe = -1
  connect%from_universe = -1
  connect%at_element = ''
  connect%at_ele_index = -1
  connect%at_s = -1
  connect%match_to_design = .false.

  read (iu, nml = tao_connected_uni_init, iostat = ios)

  if (ios == 0) then
    if (ix_universe == -1) then
      call out_io (s_abort$, r_name, &
            'INIT: READ TAO_CONNECTED_UNI_INIT NAMELIST HAS NOT SET IX_UNIVERSE!')
      call err_exit
    endif
    call out_io (s_blank$, r_name, &
        'Init: Read tao_connected_uni_init namelist for universe \i3\ ', ix_universe)
    i = ix_universe
    call init_connected_uni (s%u(i), connect, i)
    cycle
  elseif (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_CONNECTED_UNI_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_connected_uni_init)  ! generate an error message
    enddo
  endif

  close (iu)
  exit

enddo


!-----------------------------------------------------------------------
! Init Beam

! Do not initialize both beam and macro
if (s%global%track_type /= 'macro') then

  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

  do i = 1, size(s%u)
    s%u(i)%beam0_file = ''
  enddo

  ! defaults
  do 
    ix_universe = -1
    beam_init%a_norm_emitt  = 0.0
    beam_init%b_norm_emitt  = 0.0
    beam_init%dPz_dz = 0.0
    beam_init%center(:) = 0.0
    beam_init%bunch_charge = 0.0
    beam_init%ds_bunch = 1
    beam_init%sig_z   = 0.0
    beam_init%sig_e   = 0.0
    beam_init%renorm_center = .true.
    beam_init%renorm_sigma = .true.
    beam_init%n_bunch = 1
    beam_init%n_particle  = -1
    beam0_file = tao_com%beam0_file        ! From the command line
    beam_all_file = tao_com%beam_all_file  ! From the command line
    save_beam_at = ''
    ix_track_start = -1
    ix_track_end = -1

    read (iu, nml = tao_beam_init, iostat = ios)

    if (ios == 0) then
      call out_io (s_blank$, r_name, &
              'Init: Read tao_beam_init namelist for universe \i3\ ', ix_universe)
      if (ix_universe == -1) then
        do i = 1, size(s%u)
          call init_beam(s%u(i))
        enddo
      else
        call init_beam(s%u(ix_universe))
      endif
      cycle
    elseif (ios > 0 .and. s%global%track_type == "beam") then
      call out_io (s_abort$, r_name, 'INIT: TAO_BEAM_INIT NAMELIST READ ERROR!')
      rewind (iu)
      do
        read (iu, nml = tao_beam_init)  ! generate an error message
      enddo
    endif

    close (iu)
    exit

  enddo


!-----------------------------------------------------------------------
! Init macroparticles
 
elseif(s%global%track_type == 'macro') then
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

  ! defaults
  do
    ix_universe = -1
    macro_init%a%norm_emit  = 0.0
    macro_init%b%norm_emit  = 0.0
    macro_init%dPz_dz = 0.0
    macro_init%center(:) = 0.0
    macro_init%ds_bunch = 1
    macro_init%sig_z   = 10e-6
    macro_init%sig_e   = 10e-3
    macro_init%sig_e_cut = 3
    macro_init%sig_z_cut = 3
    macro_init%n_bunch = 1
    macro_init%n_slice = 1
    macro_init%n_macro = 1
    macro_init%n_part  = 1e10
    read (iu, nml = tao_macro_init, iostat = ios)
    if (ios == 0) then
      if (ix_universe == -1) then
        call out_io (s_abort$, r_name, &
                'INIT: READ TAO_MACRO_INIT NAMELIST HAS NOT SET IX_UNIVERSE!')
        call err_exit
      endif
      call out_io (s_blank$, r_name, &
              'Init: Read tao_macro_init namelist for universe \i3\ ', ix_universe)
      i = ix_universe
      call init_macro(s%u(i), macro_init)  ! generate an error message
      cycle
    elseif (ios > 0) then
      call out_io (s_abort$, r_name, 'INIT: TAO_MACRO_INIT NAMELIST READ ERROR!')
      rewind (iu)
      do
        read (iu, nml = tao_macro_init)  ! generate an error message
      enddo
    endif

    close (iu)
    exit

  enddo

endif

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

subroutine init_orbits ()

implicit none

integer i

do i = 1, size(s%u)

  n = s%u(i)%design%lat%n_ele_max
  if (allocated(s%u(i)%model%orb)) then
    deallocate (s%u(i)%model%orb, s%u(i)%model%bunch_params)
    deallocate (s%u(i)%design%orb, s%u(i)%design%bunch_params)
    deallocate (s%u(i)%base%orb, s%u(i)%base%bunch_params)
  endif
  allocate (s%u(i)%model%orb(0:n), s%u(i)%model%bunch_params(0:n))
  allocate (s%u(i)%design%orb(0:n), s%u(i)%design%bunch_params(0:n))
  allocate (s%u(i)%base%orb(0:n), s%u(i)%base%bunch_params(0:n))
  ! Specify initial conditions
  s%u(i)%design%orb(0)%vec = 0.0
  call polar_to_spinor (spin, s%u(i)%design%orb(0))
  call init_lat (s%u(i)%model%lat, s%u(i)%design%lat%n_ele_max)
  call init_lat (s%u(i)%base%lat, s%u(i)%design%lat%n_ele_max)

enddo
  
end subroutine init_orbits

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize universe connections
!-

subroutine init_connected_uni (u, connect, this_uni_index)

implicit none

type (tao_universe_struct) u
type (tao_universe_struct), pointer ::  from_uni
type (tao_connected_uni_input) connect 
integer this_uni_index

character(40) class, ele_name

integer j, ix

!
if (connect%from_universe .eq. 0 .or. connect%at_element .eq. "none") then
  u%connect%connected = .false.
  return
endif

if (connect%from_universe .ge. this_uni_index) then
  call out_io (s_abort$, r_name, &
        "A universe can only inject into a universe with a greater universe index")
  call err_exit
endif
  
u%connect%connected = .true.
u%connect%from_uni = connect%from_universe
from_uni => s%u(connect%from_universe)

call init_ele (u%connect%connect_ele)
u%connect%connect_ele%key = match$
u%connect%match_to_design = connect%match_to_design
if (u%connect%match_to_design) u%connect%use_connect_ele = .true.
          
! find extraction element
call string_trim (connect%at_element, ele_name, ix)
if (ix /= 0) then
  if (connect%at_s /= -1 .or. connect%at_ele_index /= -1) then
    call out_io (s_error$, r_name, &
        "INIT UNIVERSE CONNECTIONS: CANNOT SPECIFY AN ELEMENT, IT'S INDEX OR POSITION AT SAME TIME!", &
        "Will use element name.")
  endif
  if (ele_name == "end") then
    u%connect%from_uni_s  = from_uni%design%lat%ele(from_uni%design%lat%n_ele_track)%s
    u%connect%from_uni_ix_ele = from_uni%design%lat%n_ele_track
  else
    ! using element name 
    ! find last element with name
    do j = from_uni%design%lat%n_ele_track, 0, -1
      if (ele_name(1:ix) == trim(from_uni%design%lat%ele(j)%name)) then
        u%connect%from_uni_s = from_uni%design%lat%ele(j)%s
        u%connect%from_uni_ix_ele = j
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
        "INIT UNIVERSE CONNECTION: CANNOT SPECIFY AN ELEMENT, IT'S INDEX OR POSITION AT SAME TIME!", &
        "Will use element index.")
  endif
    u%connect%from_uni_s = from_uni%design%lat%ele(connect%at_ele_index)%s
    u%connect%from_uni_ix_ele = connect%at_ele_index
else
  ! using s position
  if (s%global%track_type /= 'single' ) then
    call out_io (s_abort$, r_name, &
     "CANNOT SPECIFY ARBITRARY S POSITION FOR COUPLING IF NOT TRACKING A SINGLE PARTICLE")
    call err_exit
  endif
  !FIX_ME: get ix_ele for element right before this s position
  u%connect%from_uni_s = connect%at_s
endif

end subroutine init_connected_uni

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize the beams. Determine which element to track beam to
!

subroutine init_beam (u)

use tao_read_beam_mod

implicit none

type (tao_universe_struct), target :: u
type (lat_struct), pointer :: lat

real(rp) v(6), bunch_charge, gamma
integer i, ix, iu, n_part, ix_class, n_bunch, n_particle
character(60) at, class, ele_name, line
integer, allocatable, save :: ix_eles(:)

! Set tracking start/stop

u%ix_track_start = ix_track_start
u%ix_track_end   = ix_track_end

! The emittance set in the tao init file takes priority over the emittance set
! in the lattice file.

lat => u%design%lat
call convert_total_energy_to (lat%E_tot, lat%param%particle, gamma)

if (beam_init%a_norm_emitt /= 0) then
  lat%a%emit = beam_init%a_norm_emitt / gamma
else
  beam_init%a_norm_emitt = lat%a%emit * gamma
endif

if (beam_init%b_norm_emitt /= 0) then
  lat%b%emit = beam_init%b_norm_emitt / gamma
else
  beam_init%b_norm_emitt = lat%b%emit * gamma
endif

u%beam_init = beam_init
u%beam0_file = beam0_file
u%beam_all_file = beam_all_file
u%design%orb(0)%vec = beam_init%center

! No initialization for a circular lattice
if (u%design%lat%param%lattice_type == circular_lattice$) return

! Find where to save the beam at

u%ele%save_beam = .false.
u%ele(0)%save_beam = .true.
u%ele(u%design%lat%n_ele_track)%save_beam = .true.

do i = 1, size(save_beam_at)
  if (save_beam_at(i) == '') exit
  call tao_ele_locations_given_name(u, save_beam_at(i), ix_eles, err, .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD SAVE_BEAM_AT ELEMENT: ' // save_beam_at(i))
    cycle
  endif
  do k = 1, size(ix_eles)
    j = ix_eles(k)
    u%ele(j)%save_beam = .true.
  enddo
enddo
  
if (allocated (u%save_beam_at)) deallocate (u%save_beam_at)
allocate (u%save_beam_at(i-1))
if (i > 1) u%save_beam_at(1:i-1) = save_beam_at(1:i-1)

! If beam_all_file is set then read in the beam distributions.

if (u%beam_all_file /= '') then
  tao_com%use_saved_beam_in_tracking = .true.
  call open_beam_file (beam_all_file)
  
  do
    call read_beam_params (j, n_bunch, n_particle, bunch_charge)
    if (j == -1) exit
    u%beam_init%n_bunch = n_bunch
    u%beam_init%n_particle = n_particle
    u%beam_init%bunch_charge = bunch_charge
    call read_beam (u%ele(j)%beam)
  enddo  

  call out_io (s_info$, r_name, 'Read beam_all file: ' // u%beam_all_file)
  call close_beam_file ()

endif

if (u%connect%connected) u%connect%injecting_beam = u%current_beam
if (allocated(ix_eles)) deallocate (ix_eles)

end subroutine init_beam

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize the macroparticles. Determine which element to track beam to
!

subroutine init_macro(u, macro_init)

implicit none

type (tao_universe_struct) u
type (macro_init_struct) macro_init

!

if (u%design%lat%param%lattice_type == circular_lattice$) then
  call out_io (s_blank$, r_name, "***")
  call out_io (s_blank$, r_name, &
                 "Macroparticle tracking through a circular lattice.")
  call out_io (s_blank$, r_name, &
         "Twiss parameters and initial orbit will be found from the closed orbit.")
  call out_io (s_blank$, r_name, "***")
endif

u%macro_beam%macro_init = macro_init
u%design%orb(0)%vec = macro_init%center

! Don't initialize beams in circular lattice
if (u%design%lat%param%lattice_type == circular_lattice$) return
    
! This is just to get things allocated
call init_macro_distribution (u%macro_beam%beam, macro_init, u%design%lat%ele(0), .true.)
if (u%connect%connected) &
  call init_macro_distribution (u%connect%injecting_macro_beam, &
                                             macro_init, u%design%lat%ele(0), .true.)

! keep track of where macros are lost
if (associated (u%macro_beam%ix_lost)) deallocate (u%macro_beam%ix_lost)
allocate (u%macro_beam%ix_lost(macro_init%n_bunch, macro_init%n_slice, macro_init%n_macro))
u%macro_beam%ix_lost(:,:,:) = -1

end subroutine init_macro

end subroutine tao_init_global

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
!   data_file      -- Character(*): Tao data initialization file.
!-

subroutine tao_init_data (data_file)

use tao_data_mod
use tao_lattice_calc_mod
use tao_input_struct
use bmad_parser_mod
use random_mod
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data
type (tao_data_input) data(n_data_minn:n_data_maxx) ! individual weight 
type (tao_d1_data_struct), pointer :: d1_ptr

real(rp) default_weight        ! default merit function weight
real(rp) default_step          ! default "small" step size
real(rp) default_data_noise         ! default noise for data type

integer ios, iu, i, j, i2, j2, k, ix, n_uni, num
integer n, n_universes, iostat, ix_universe, n_max
integer n_d1_data, ix_ele, ix_min_data, ix_max_data, ix_d1_data

integer, automatic :: n_d2_data(size(s%u))

character(*) data_file
character(40) :: r_name = 'tao_init_data'
character(200) file_name, beam0_file, beam_all_file
character(40) name,  universe, default_universe, default_data_type
character(40) default_merit_type, default_attribute, data_type, default_data_source
character(40) use_same_lat_eles_as, search_for_lat_eles

character(60) save_beam_at(100)
character(100) line

logical err, free, gang
logical, automatic :: mask(size(s%u))

namelist / tao_d2_data / d2_data, n_d1_data, default_merit_type, universe, &
                  default_data_noise
  
namelist / tao_d1_data / d1_data, data, ix_d1_data, ix_min_data, ix_max_data, &
                   default_weight, default_data_type, default_data_source, &
                   use_same_lat_eles_as, search_for_lat_eles
                     

!-----------------------------------------------------------------------
! Find out how many d2_data structures we need for each universe

call tao_open_file ('TAO_INIT_DIR', data_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Data File: ' // file_name)

n_d2_data = 0

do 
  universe = '*'
  read (iu, nml = tao_d2_data, iostat = ios, err = 9100)
  if (ios /= 0) exit
  if (universe == '*') then
    n_d2_data = n_d2_data + 1
  else
    read (universe, *, iostat = ios) n_uni
    if (ios /= 0 .or. n_uni > size(s%u)) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
    n_d2_data(n_uni) = n_d2_data(n_uni) + 1
  endif
enddo

do i = 1, size(s%u)
  call init_data_in_universe (s%u(i), n_d2_data(i))
enddo

! Init data

rewind (iu)

do 
  mask(:) = .true.      ! set defaults
  d2_data%name       = ''
  universe           = '*'
  default_merit_type = 'target'
  default_data_noise = 0.0

  read (iu, nml = tao_d2_data, iostat = ios, err = 9100)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, &
                      'Init: Read tao_d2_data namelist: ' // d2_data%name)
    
  if (universe == '*') then
    uni_loop1: do i = 1, size(s%u)

    ! check if this data type has already been defined for this universe
    do k = 1, size(s%u(i)%d2_data)
      if (trim(s%u(i)%d2_data(k)%name) == trim(d2_data%name)) then
        mask(i) = .false.
        cycle uni_loop1
      endif
    enddo
      
    call d2_data_stuffit (s%u(i), i)
    enddo uni_loop1
  else
    read (universe, *, iostat = ios) n_uni
    if (ios /= 0 .or. n_uni > size(s%u)) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
    if (n_uni == 0) then
      call out_io (s_error$, r_name, &
              '"UNIVERSE == 0" MUST BE REPLACED BY "UNIVERSE = *"', &
              ' IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
    call d2_data_stuffit (s%u(n_uni), n_uni)
  endif

  do k = 1, n_d1_data
    use_same_lat_eles_as = ''
    search_for_lat_eles  = ''
    d1_data%name        = ''
    default_weight      = 0      ! set default
    default_data_type   = ''
    default_data_source = 'lattice'
    data(:)%data_type   = ''
    data(:)%merit_type  = default_merit_type 
    data(:)%name        = ''
    data(:)%merit_type  = ''
    data(:)%ele_name    = ''
    data(:)%ele0_name   = ''
    data(:)%meas        = real_garbage$  ! used to tag when %meas_value is set in file
    data(:)%weight      = 0.0
    data(:)%ix_bunch    = 0
    data(:)%data_noise  = real_garbage$
    data(:)%data_source = ''
    data(:)%good_user   = .true.
    data(:)%relative    = ''
    ix_min_data         = int_garbage$
    ix_max_data         = int_garbage$

    read (iu, nml = tao_d1_data, err = 9150)

    ! Convert old format to new

    if (data(0)%ele_name(1:7) == 'SEARCH:') then
      call string_trim(data(0)%ele_name(8:), search_for_lat_eles, ix)
    elseif (data(0)%ele_name(1:5) == 'SAME:') then
      call string_trim (data(0)%ele_name(6:), use_same_lat_eles_as, ix)
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
    do i = lbound(data, 1), ubound(data, 1)
      ! 'beam_tracking' is old syntax.
      if (data(i)%data_source == 'beam_tracking') data(i)%data_source = 'beam'
      if (data(i)%ele0_name /= '' .and. data(i)%ele_name == '') then
        write (line, '(4a, i0, a)') trim(d2_data%name), '.', trim(d1_data%name), '[', i, ']'
        call out_io (s_abort$, r_name, &
              'ERROR: ELE_NAME IS BLANK BUT ELE0_NAME IS NOT FOR: ' // line)
        call err_exit
      endif
    enddo
    call out_io (s_blank$, r_name, &
                      'Init: Read tao_d1_data namelist: ' // d1_data%name)
    if (universe == '*') then          ! * => use all universes
      uni_loop2: do i = 1, size(s%u)

      ! check if this data type has already been defined for this universe
      if (.not. mask(i)) cycle uni_loop2
      
        call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used)
      enddo uni_loop2
    else
      call d1_data_stuffit (k, s%u(n_uni), s%u(n_uni)%n_d2_data_used)
    endif
  enddo

enddo

close (iu)

!-----------------------------------------------------------------------
! Init ix_data array

do i = 1, size(s%u)
  call init_ix_data (s%u(i))
enddo

return

!-----------------------------------------------------------------------
! namelist read error.

9100 continue
call out_io (s_error$, r_name, 'TAO_D2_DATA NAMELIST READ ERROR.')
rewind (iu)
do
  read (iu, nml = tao_d2_data)  ! force printing of error message
enddo

9150 continue
call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST READ ERROR.')
rewind (iu)
do
  read (iu, nml = tao_d1_data)  ! force printing of error message
enddo

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
contains

subroutine init_data_in_universe (u, n_d2_data)

type (tao_universe_struct) :: u
integer i, n_d2_data

!

u%is_on = .true.          ! turn universe on
u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
u%n_data_used = 0         ! size of s%u(i)%data(:) array
u%ix_rad_int_cache = 0

! allocate and set defaults

if (n_d2_data /= 0) then
  if (associated(u%d2_data)) deallocate (u%d2_data)
  allocate (u%d2_data(n_d2_data))
  do i = 1, n_d2_data
    u%d2_data(i)%descrip = ''
  enddo
  u%d2_data%name = ''  ! blank name means it doesn't exist
endif

if (tao_com%n_data_max /= 0) then
  if (associated(u%data)) deallocate (u%data)
  allocate (u%data(tao_com%n_data_max))
  u%data(:)%exists = .false.       ! set default
  u%data(:)%good_meas  = .false.   ! set default
  u%data(:)%good_ref   = .false.   ! set default
  u%data(:)%good_user  = .true.    ! set default
  u%data(:)%good_opt   = .true.
  u%data(:)%merit_type = 'target'  ! set default
  u%data(:)%ele_name   = ''
  u%data(:)%ix_ele     = -1
  u%data(:)%ele0_name  = ''
  u%data(:)%ix_ele0    = 0 ! by default, data relative to beginning of lattice
endif

! This is needed to keep the totalview debugger happy.

if (associated(u%dmodel_dvar)) deallocate (u%dmodel_dvar)
allocate (u%dmodel_dvar(1,1))
  
end subroutine init_data_in_universe

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine d2_data_stuffit (u, ix_uni)

type (tao_universe_struct), target :: u

integer nn, ix_uni

! Setup another d2_data structure.

u%n_d2_data_used = u%n_d2_data_used + 1
nn = u%n_d2_data_used

if (size(u%d2_data) < nn) then
  call out_io (s_error$, r_name, &
              'N_D2_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
  call err_exit
endif

u%d2_data(nn)%name = d2_data%name 
u%d2_data(nn)%ix_uni = ix_uni

! allocate memory for the u%d1_data structures

if (associated(u%d2_data(nn)%d1)) deallocate (u%d2_data(nn)%d1)
allocate(u%d2_data(nn)%d1(n_d1_data))

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine d1_data_stuffit (i_d1, u, n_d2)

type (tao_universe_struct), target :: u
integer i, n1, n2, ix, k, ix1, ix2, j, jj, n_d2

integer i_d1, num_hashes

character(40) search_string, d2_d1_name
character(20) fmt

integer, allocatable, save :: ix_eles(:)
logical emit_here

!

u%d2_data(n_d2)%d1(i_d1)%d2 => u%d2_data(n_d2)  ! point back to the parent
if (d1_data%name == '') then
  write (u%d2_data(n_d2)%d1(i_d1)%name, '(i0)') i_d1
else
  u%d2_data(n_d2)%d1(i_d1)%name = d1_data%name    ! stuff in the data
endif

d2_d1_name = u%d2_data(n_d2)%name // '.' // u%d2_data(n_d2)%d1(i_d1)%name

! Check if we are searching for elements or repeating elements
! and record the element names in the data structs.
    
if (search_for_lat_eles /= '') then
  call tao_find_elements (u, search_for_lat_eles, ix_eles)
  if (size(ix_eles) == 0) then
    call out_io (s_warn$, r_name, &
      'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_string, &
      'WHILE SETTING UP DATA ARRAY: ' // d1_data%name)
    return
  endif
  ! finish finding data array limits
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(ix_eles)
  u%n_data_used = n2
  if (ix_min_data == int_garbage$) ix_min_data = 1
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                  'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif
  ! get element names
  jj = n1
  do k = lbound(ix_eles, 1), ubound(ix_eles, 1)
    j = ix_eles(k)
    if (jj .gt. n2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    u%data(jj)%ele_name = u%design%lat%ele(j)%name
    u%data(jj)%ix_ele   = j
    u%data(jj)%exists   = .true.
    jj = jj + 1
  enddo

  u%data(n1:n2)%good_user   = data(ix1:ix2)%good_user
  u%data(n1:n2)%weight      = data(ix1:ix2)%weight
  u%data(n1:n2)%ele0_name   = data(ix1:ix2)%ele0_name
  u%data(n1:n2)%ix_bunch    = data(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source = data(ix1:ix2)%data_source

! use_same_lat_eles_as

elseif (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_data (err, name, d1_ptr = d1_ptr, ix_uni = u%ix_uni)
  if (err .or. .not. associated(d1_ptr)) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(d1_ptr%d)
  u%n_data_used = n2
  ix_min_data = lbound(d1_ptr%d, 1)
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  u%data(n1:n2)%ele_name    = d1_ptr%d%ele_name
  u%data(n1:n2)%ix_ele      = d1_ptr%d%ix_ele
  u%data(n1:n2)%ele0_name   = d1_ptr%d%ele0_name
  u%data(n1:n2)%ix_ele0     = d1_ptr%d%ix_ele0
  u%data(n1:n2)%exists      = d1_ptr%d%exists
  u%data(n1:n2)%data_source = d1_ptr%d%data_source
  u%data(n1:n2)%relative    = d1_ptr%d%relative

! Not SEARCH or SAME:

else

  if (ix_min_data == int_garbage$) ix_min_data = 1
  if (ix_max_data == int_garbage$) then
    do i = ubound(data, 1), lbound(data, 1), -1
      if (data(i)%ele_name /= '' .or. data(i)%data_type /= '') then
        ix_max_data = i
        exit
      endif
    enddo
  endif

  if (ix_max_data == int_garbage$) then
    call out_io (s_error$, r_name, 'NO DATA FOUND FOR: ' // d2_d1_name)
    return
  endif

  n1 = u%n_data_used + 1
  n2 = u%n_data_used + ix_max_data - ix_min_data + 1
  ix1 = ix_min_data
  ix2 = ix_max_data
  u%n_data_used = n2
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  ! Transfer info from the input structure

  u%data(n1:n2)%good_user   = data(ix1:ix2)%good_user
  u%data(n1:n2)%weight      = data(ix1:ix2)%weight
  u%data(n1:n2)%ele_name    = data(ix1:ix2)%ele_name
  u%data(n1:n2)%ele0_name   = data(ix1:ix2)%ele0_name
  u%data(n1:n2)%ix_bunch    = data(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source = data(ix1:ix2)%data_source

  ! Find elements associated with the data

  do j = n1, n2

    call tao_hook_does_data_exist (u%data(j))
    if (u%data(j)%exists) cycle

    if (u%data(j)%ele_name == '') cycle
    call str_upcase (u%data(j)%ele_name, u%data(j)%ele_name)
    call element_locator (u%data(j)%ele_name, u%design%lat, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT NOT LOCATED: ' // &
                                                       u%data(j)%ele_name)
      u%data(j)%exists = .false.
      cycle
    endif

    u%data(j)%ix_ele = ix
    u%data(j)%exists = .true.

    if (u%data(j)%ele0_name == '') cycle
    call str_upcase (u%data(j)%ele0_name, u%data(j)%ele0_name)
    call element_locator (u%data(j)%ele0_name, u%design%lat, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT2 NOT LOCATED: ' // &
                                                       u%data(j)%ele0_name)
      u%data(j)%exists = .false.
      cycle
    endif
    u%data(j)%ix_ele0 = ix
  enddo

endif

!-----------------------------------------------------------
! If %meas_value was set then %good_meas is set to True

u%data(n1:n2)%data_type   = data(ix1:ix2)%data_type
u%data(n1:n2)%merit_type  = data(ix1:ix2)%merit_type
u%data(n1:n2)%weight      = data(ix1:ix2)%weight

u%data(n1:n2)%meas_value = data(ix1:ix2)%meas
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
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = &
                                                    default_data_type
endif


! set data noise (not applicable to all data types)

if (d2_data%name .eq. "bpm" .or. d2_data%name .eq. "wire") then
  do j = n1, n2
    u%ele(u%data(j)%ix_ele)%data_noise = default_data_noise
  enddo
  do j = lbound(data,1), ubound(data,1)
    if (data(j)%data_noise /= real_garbage$) &
      u%ele(u%data(n1+j-ix1)%ix_ele)%data_noise = data(j)%data_noise
  enddo
endif                   

! now for some family guidance...
! point the children to the grandchildren in the big data array

call tao_point_d1_to_data (u%d2_data(n_d2)%d1(i_d1)%d, &
                                      u%data(n1:n2), ix_min_data, n1)

! point the %data back to the d1_data_struct

do j = n1, n2
  u%data(j)%d1 => u%d2_data(n_d2)%d1(i_d1)
  if (u%data(j)%weight == 0) u%data(j)%weight = default_weight
  if (u%data(j)%merit_type == '') u%data(j)%merit_type = default_merit_type
  if (u%data(j)%data_source == '') u%data(j)%data_source = default_data_source
  ! old style is to use "emittance." instead of "emit."
  ix = index(u%data(j)%data_type, 'emittance.')
  if (ix /= 0) u%data(j)%data_type = u%data(j)%data_type(1:ix-1) // &
                                         'emit.' // u%data(j)%data_type(ix+10:)
enddo

! point the children back to the mother    

u%d2_data(n_d2)%d1(i_d1)%d2 => u%d2_data(n_d2)

! do we need to do the radiation integrals?

do j = n1, n2
  data_type = u%data(j)%data_type
  emit_here = (index(data_type, 'emit.') /= 0)
  if (emit_here .and. u%data(j)%data_source == 'lattice') u%do_synch_rad_int_calc = .true. 
  if (data_type(1:2) == 'i5') u%do_synch_rad_int_calc = .true. 
  if (data_type(1:6) == 'chrom.') u%do_chrom_calc = .true.

  if (u%design%lat%param%lattice_type == circular_lattice$ .and. &
              (data_type(1:6)  == 'chrom.' .or. data_type(1:2) == 'i5' .or. &
               data_type(1:13) == 'unstable_ring' .or. emit_here .or. &
               data_type(1:17) == 'multi_turn_orbit.')) then
    u%data(j)%exists = .true.
    if (u%data(j)%ele_name /= '') then
      call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // data_type, &
                        'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // u%data(j)%ele_name)
      call err_exit
    endif
    cycle
  endif

enddo

end subroutine d1_data_stuffit

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Defines what datums to evaluate at each element in specified universe

subroutine init_ix_data (u)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: data

integer, automatic :: n_data(-2:u%design%lat%n_ele_max)
integer, automatic :: ix_next(-2:u%design%lat%n_ele_max)

integer j, k, ix_ele

n_data(:) = 0

! allocate the ix_data array
if (associated(u%ix_data)) deallocate(u%ix_data)
allocate(u%ix_data(-2:u%design%lat%n_ele_max))

! Since some information gets lost during tracking (like beam distributions),
!   find where each datum gets evaluated when tao_load_data_array is called.
! ix_ele = -1  ! Gets evaluated after all tracking
! ix_ele = -2  ! Does not get evaluated by tao_lattice_calc_mod

do j = 1, size(u%data)
  data => u%data(j)
  if (.not. data%exists) cycle
  if (data%data_type(1:17) == 'multi_turn_orbit.') then
    ix_ele = -2
  elseif (data%data_type(1:2) == 'i5') then
    ix_ele = -1
  elseif (data%ix_ele == -1) then
    ix_ele = -1
  elseif (index(data%data_type, 'emit.') /= 0 .and. data%data_source == 'lattice') then
    ix_ele = -1
  elseif (data%ix_ele0 > data%ix_ele) then
    ix_ele = u%model%lat%n_ele_track
  else
    ix_ele = data%ix_ele
  endif
  n_data(ix_ele) = n_data(ix_ele) + 1
enddo
  
! allocate ix_ele array for each element
do j = lbound(u%ix_data, 1), ubound(u%ix_data, 1)
  if (associated(u%ix_data(j)%ix_datum)) deallocate (u%ix_data(j)%ix_datum)
  if (n_data(j) == 0) cycle
  allocate (u%ix_data(j)%ix_datum(n_data(j)))
enddo

! used for keeping track of current datum index in each ix_ele element
ix_next(:) = 1
  
! setup ix_ele array for each element
! This is the point where the datum is evaluated
! if ix_ele0 > ix_ele then there is "wrap around"
do j = 1, size(u%data)
  data => u%data(j)
  if (.not. data%exists) cycle
  if (data%data_type(1:17) == 'multi_turn_orbit.') then
    ix_ele = -2
  elseif (data%data_type(1:2) == 'i5') then
    ix_ele = -1
  elseif (data%ix_ele == -1) then
    ix_ele = -1
  elseif (index(data%data_type, 'emit.') /= 0 .and. data%data_source == 'lattice') then
    ix_ele = -1
  elseif (data%ix_ele0 > data%ix_ele) then
    ix_ele = u%model%lat%n_ele_track
  else
    ix_ele = data%ix_ele
  endif
  u%ix_data(ix_ele)%ix_datum(ix_next(ix_ele)) = j
  ix_next(ix_ele) = ix_next(ix_ele) + 1
enddo

end subroutine init_ix_data

end subroutine tao_init_data

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
!   var_file       -- Character(*): Tao variable initialization file.
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

integer ios, iu, i, j, i2, j2, k, ix, n_uni, num
integer n, n_universes, iostat, ix_universe, n_max
integer ix_min_var, ix_max_var, ix_ele, n_v1, n_v1_var_max

character(*) var_file
character(40) :: r_name = 'tao_init_variables'
character(40) name, universe, default_universe, default_data_type
character(40) default_merit_type, default_attribute
character(40) use_same_lat_eles_as, search_for_lat_eles
character(8) default_key_bound
character(200) file_name
character(8), allocatable, save :: default_key_b(:)

character(100) line

logical err, free, gang
logical searching
logical, allocatable, save :: dflt_good_unis(:), good_unis(:)
                     
namelist / tao_var / v1_var, var, default_weight, default_step, default_key_delta, &
                    ix_min_var, ix_max_var, default_universe, default_attribute, &
                    default_low_lim, default_high_lim, default_merit_type, &
                    use_same_lat_eles_as, search_for_lat_eles, default_key_bound

!-----------------------------------------------------------------------
! Init

if (associated(s%var)) deallocate (s%var)
allocate (s%var(tao_com%n_var_max))
s%var(:)%good_opt  = .true.
s%var(:)%exists    = .false.
s%var(:)%good_var  = .true.
s%var(:)%good_user = .true.

s%n_var_used = 0

! First count how many v1_var definitions there are

if (associated(s%v1_var)) deallocate (s%v1_var)

call tao_open_file ('TAO_INIT_DIR', var_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Variable File: ' // file_name)

n_v1_var_max = 0

do
  read (iu, nml = tao_var, iostat = ios, err = 9200)
  if (ios < 0) exit
  n_v1_var_max = n_v1_var_max + 1
enddo

n = n_v1_var_max + 1
allocate (s%v1_var(n), default_key_b(n), default_key_d(n))
s%v1_var%name = ''  ! blank name means it doesn't (yet) exist
s%n_v1_var_used = 0       ! size of s%v1_var(:) array

! Read some defaults

rewind (iu)

n_v1 = 0
do
  default_key_bound = ''
  default_key_delta = 0
  read (iu, nml = tao_var, iostat = ios, err = 9200)
  if (ios < 0) exit
  n_v1 = n_v1 + 1
  default_key_b(n_v1) = default_key_bound
  default_key_d(n_v1) = default_key_delta
enddo

! Now fill in all the information

rewind (iu)

allocate (dflt_good_unis(size(s%u)), good_unis(size(s%u)))

n_v1 = 0
do
  n_v1 = n_v1 + 1
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
  var%name           = ''
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

  read (iu, nml = tao_var, iostat = ios, err = 9200)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, &
                        'Init: Read tao_var namelist: ' // v1_var%name)
  call str_upcase (default_attribute, default_attribute)

  ! Convert old format to new

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
    call var_stuffit1 (v1_var_ptr)
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

      call var_stuffit2 (good_unis, v1_var_ptr%v(j), searching)

    enddo

  else   ! If clone...
    do i = 1, size(s%u)
      if (.not. dflt_good_unis(i)) cycle
      call var_stuffit1 (v1_var_ptr)
      write (v1_var_ptr%name, '(2a, i0)') trim(v1_var_ptr%name), '_u', i
      good_unis = .false.
      good_unis(i) = .true.
      do j = lbound(v1_var_ptr%v, 1), ubound(v1_var_ptr%v, 1)
        call var_stuffit2 (good_unis, v1_var_ptr%v(j), searching)
      enddo
    enddo
  endif

enddo

close (iu)
deallocate (dflt_good_unis, good_unis)
deallocate (default_key_b, default_key_d)

return

!-----------------------------------------------------------------------
! namelist read error.

9200 continue
call out_io (s_error$, r_name, 'TAO_VAR NAMELIST READ ERROR.')
rewind (iu)
do
  read (iu, nml = tao_var)  ! force printing of error message
enddo

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

! stuff common to all universes

subroutine var_stuffit1 (v1_var_ptr)

type (tao_v1_var_struct), pointer :: v1_var_ptr
type (tao_v1_var_struct), pointer :: v1_ptr

character(20) fmt
character(60) search_string

integer i, iu, ip, j, jj, k, kk, nn, n1, n2, ix1, ix2, num_hashes, ix
integer num_ele, ios, ixx1, ixx2

integer, allocatable, save :: ix_eles(:)

! count number of v1 entries

s%n_v1_var_used = s%n_v1_var_used + 1
nn = s%n_v1_var_used
v1_var_ptr => s%v1_var(nn)
v1_var_ptr%name = v1_var%name

! If reusing a previous element list...

if (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_var (err, name, v1_ptr = v1_ptr)
  if (err .or. .not. associated(v1_ptr)) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "USE_SAME_LAT_ELES_AS": ' // name)
    call err_exit
  endif
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + size(v1_ptr%v)
  s%n_var_used = n2
  ix_min_var = lbound(v1_ptr%v, 1)
  ix1 = ix_min_var
  ix2 = ix1 + (n2 - n1)
  if (n2 > size(s%var)) then
    call out_io (s_abort$, r_name, &
                'N_VAR_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  s%var(n1:n2)%ele_name    = v1_ptr%v%ele_name
  s%var(n1:n2)%s           = v1_ptr%v%s
  s%var(n1:n2)%ix          = v1_ptr%v%ix

  do n = n1, n2
    ix = ix1 + (n - n1)
    ip = 1 + (n - n1) 

    
    s%var(n)%good_user = v1_ptr%v(ip)%good_user
    if (var(ix)%good_user /= '') read (var(ix)%good_user, *, iostat = ios) s%var(n)%good_user
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "GOOD_USER" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif

    s%var(n)%key_bound = v1_ptr%v(ip)%key_bound
    if (var(ix)%key_bound /= '') read (var(ix)%key_bound, *, iostat = ios) s%var(n)%key_bound
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "KEY_BOUND" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif

    s%var(n)%key_delta = v1_ptr%v(ip)%key_delta
    if (var(ix)%key_delta /= '') s%var(n)%key_delta = var(ix)%key_delta

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

    s%var(n)%key_bound = v1_ptr%v(ip)%key_bound

  enddo

  call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var, n1)
  return

endif

!------------------------------
! are we searching for elements?

if (search_for_lat_eles /= '') then
  searching = .true.
  search_string = '-no_slaves ' // trim(search_for_lat_eles)
  if (any(var%universe /= '')) then
    call out_io (s_abort$, r_name, &
           "CANNOT SPECIFY INDIVIDUAL UNIVERSES WHEN SEARCHING FOR VARIABLES")
    call err_exit
  endif
  ! search through all universes specified
  num_ele = 0
  do iu = 1, size(s%u)
    if (.not. dflt_good_unis(iu)) cycle
    call tao_find_elements (s%u(iu), search_string, ix_eles)
    num_ele = num_ele + size(ix_eles)
  enddo
  if (size(ix_eles) == 0) then
    call out_io (s_warn$, r_name, &
                'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_string, &
                'WHILE SETTING UP VARIABLE ARRAY: ' // v1_var%name)
    return
  endif

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + num_ele
  ix1 = ix_min_var
  ix2 = num_ele - (1-ix_min_var)
  s%n_var_used = n2
  jj = n1

  do iu = 1, size(s%u)
    if (.not. dflt_good_unis(iu)) cycle
    call tao_find_elements (s%u(iu), search_string, ix_eles)
    do kk = 1, size(ix_eles)
      k = ix_eles(kk)
      if (jj .gt. n2) then
        call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT SEARCHING")
        call err_exit
      endif
      s%var(jj)%ele_name = s%u(iu)%design%lat%ele(k)%name
      s%var(jj)%s = s%u(iu)%design%lat%ele(k)%s
      s%var(jj)%ix = k
      jj = jj + 1
    enddo
  enddo

!------------------------------
! If not searching or reusing...

else  
  searching = .false.
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + ix_max_var - ix_min_var + 1
  ix1 = ix_min_var
  ix2 = ix_max_var
 
  s%n_var_used = n2

  s%var(n1:n2)%ele_name    = var(ix1:ix2)%ele_name

endif

!---------------------------

do n = n1, n2
  i = ix1 + n - n1

  if (var(i)%key_bound == '') then
    s%var(n)%key_bound = .false.
  else
    read (var(i)%key_bound, *, iostat = ios) s%var(n)%key_bound
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD "KEY_BOUND" COMPONENT FOR VAR: ' // v1_var%name)
      call err_exit
    endif
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

enddo

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
 
! Point the v1_var mother to the appropriate children in the big var array

call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var, n1)

end subroutine
  
end subroutine tao_init_variables

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine var_stuffit2 (good_unis, var, searching)

implicit none

type (tao_var_struct), target :: var
type (tao_this_var_struct), pointer :: this
type (tao_universe_struct), pointer :: u

integer i, j, n, n1, n2, iv, iu, n_tot, n_ele, ie
integer, automatic :: n_ele_in_uni(size(s%u))

character(20) :: r_name = 'var_stuffit2'
logical err, searching, good_unis(:)

! point the children back to the mother

if (allocated(var%this)) deallocate (var%this)
if (var%ele_name == '') then
  allocate (var%this(0))
  var%exists = .false.
  var%key_bound = .false.
  return
endif

! Problem: var%ix gives element index but this might be different in different universes.
! Solution: Move this to var_stuffit1

if (searching) then
  allocate (var%this(count(good_unis)))
  j = 0
  do iu = 1, size(s%u)
    if (.not. good_unis(iu)) cycle
    j = j + 1
    call tao_pointer_to_var_in_lattice (var, var%this(j), iu, &
                                                  err, .false., ix_ele = var%ix)
    if (err) then
      var%exists = .false.
      var%key_bound = .false.
      return
    endif
  enddo

  ! Veto any variable that is not free

  if (var%ix_attrib > 0) then
    do j = 1, size(var%this)
      this => var%this(j)
      u => s%u(this%ix_uni)
      if (.not. attribute_free (this%ix_ele, var%ix_attrib, u%model%lat, .false.)) then
        call out_io (s_info$, r_name, &
              'Warning: Variable: ' // tao_var1_name(var), &
            '         is trying to control a non-free attribute. %exists will be set to False.')
        var%exists = .false.
        var%key_bound = .false.
        return
      endif
    enddo
  endif

else
  n_tot = 0
  do iu = 1, size(s%u)
    if (.not. good_unis(iu)) cycle

    n_ele = 0
    do iv = 0, s%u(iu)%model%lat%n_ele_max
      if (var%ele_name /= s%u(iu)%model%lat%ele(iv)%name) cycle
      n_ele = n_ele + 1
      s%u(iu)%model%lat%ele(n_ele)%ixx = iv
    enddo

    n_ele_in_uni(iu) = n_ele
    n_tot = n_tot + n_ele
  enddo

  if (n_tot == 0) then
    call out_io (s_error$, r_name, 'ELEMENT DOES NOT EXIST: ' // var%ele_name)
    var%exists = .false.
    return
  endif

  allocate (var%this(n_tot))

  n = 0
  do iu = 1, size(s%u)
    if (.not. good_unis(iu)) cycle
    do ie = 1, n_ele_in_uni(iu)
      call tao_pointer_to_var_in_lattice (var, var%this(ie+n), iu, &
                                    err, .true., ix_ele = s%u(iu)%model%lat%ele(ie)%ixx)
      if (err) then
        var%exists = .false.
        var%key_bound = .false.
        return
      endif
    enddo
    n = n + n_ele_in_uni(iu)
  enddo

endif

var%model_value = var%this(1)%model_ptr
var%design_value = var%this(1)%model_ptr
var%base_value = var%this(1)%base_ptr
var%exists = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! This searches the lattice for the specified element and flags ix_eles(:)

subroutine tao_find_elements (u, search_string, ix_eles)

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele

character(*) search_string
character(80) string
character(40) ele_name, key_name_in
character(20) :: r_name = 'tao_find_elements'

integer key, found_key
integer i, k, ix, ii, i2, j

integer, allocatable :: ix_eles(:)
logical no_slaves, no_lords, err

!

no_slaves = .false.
no_lords = .false.

call string_trim(search_string, string, ix)

select case (string(1:ix))
case ('-no_lords') 
  no_lords = .true.
  call string_trim (string(ix+1:), string, ix)
case ('-no_slaves') 
  no_slaves = .true.
  call string_trim (string(ix+1:), string, ix)
case default
  if (string(1:1) == '-') then
    call out_io (s_abort$, r_name, "BAD SEARCH RESTRICTION: " // search_string)
    call err_exit
  endif
end select

!

call tao_ele_locations_given_name (u, string, ix_eles, err)

k = 0
do j = 1, size(ix_eles)
  ele => u%design%lat%ele(ix_eles(j))
  select case (ele%control_type)
  case (multipass_slave$, super_slave$)
    if (no_slaves) cycle
  case (girder_lord$, overlay_lord$, super_lord$)
    if (no_lords) cycle
  end select
  k = k + 1
  ix_eles(k) = ix_eles(j)
enddo

call re_allocate (ix_eles, k)

end subroutine tao_find_elements

end module
