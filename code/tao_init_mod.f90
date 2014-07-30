module tao_init_mod

use tao_mod
 
implicit none

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

logical err, xxx
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

if (init_file == '') then
  call end_bookkeeping()
  return
endif

call tao_open_file (init_file, iu, file_name, s_blank$)
call out_io (s_blank$, r_name, '*Init: Opening Init File: ' // file_name)
if (iu == 0) then
  call out_io (s_blank$, r_name, "Note: Cannot open init file for tao_params namelist read")
  call end_bookkeeping()
  return
endif

! Read tao_params

call out_io (s_blank$, r_name, 'Init: Reading tao_params namelist')
bmad_com%rel_tol_tracking = 1e-8   ! Need tighter tol for calculating derivatives
bmad_com%abs_tol_tracking = 1e-11  ! Need tighter tol for calculating derivatives
read (iu, nml = tao_params, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_PARAMS NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_params)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_params namelist found')

! transfer global to s%global
s%global = global

close (iu)

call end_bookkeeping()

!-----------------------------------------------------------------------
contains

subroutine end_bookkeeping ()

! Tao does its own bookkeeping

bmad_com%auto_bookkeeper = .false.
s%com%valid_plot_who(1:5) = (/ 'model ', 'base  ', 'ref   ', 'design', 'meas  ' /)

! Seed random number generator

call ran_seed_put (s%global%random_seed)
call ran_engine (s%global%random_engine)
call ran_gauss_converter (s%global%random_gauss_converter, s%global%random_sigma_cutoff)

if (s%com%noplot_arg_set) s%global%plot_on = .false.

end subroutine end_bookkeeping

end subroutine tao_init_global

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
character(40) track_start, track_end
character(160) beam_saved_at
character(200) file_name, beam0_file, beam_all_file
character(60), target :: save_beam_at(100)   ! old style syntax

logical err

namelist / tao_beam_init / ix_universe, beam0_file, &
  ix_track_start, ix_track_end, beam_all_file, beam_init, beam_saved_at, &
  beam_saved_at, track_start, track_end
         
!-----------------------------------------------------------------------
! Init Beams

call tao_hook_init_beam ()

if (.not. s%com%init_beam .or. init_file == '') return

!

call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)
call tao_open_file (init_file, iu, file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'CANNOT OPEN BEAM INIT FILE: ' // init_file)
  call err_exit
endif

do i = lbound(s%u, 1), ubound(s%u, 1)
  s%u(i)%beam%beam0_file = s%com%beam0_file
  s%u(i)%beam%beam_all_file = s%com%beam_all_file
  do ib = 0, ubound(s%u(i)%uni_branch, 1)
    s%u(i)%uni_branch(ib)%track_start    = ''
    s%u(i)%uni_branch(ib)%track_end      = ''
    s%u(i)%uni_branch(ib)%ix_track_start = 0
    s%u(i)%uni_branch(ib)%ix_track_end   = -1
  enddo
enddo

do 

  ! defaults
  ix_universe = -1
  beam_init%distribution_type = ''
  beam_init%a_norm_emit   = 0.0
  beam_init%b_norm_emit   = 0.0
  beam_init%a_emit        = 0.0
  beam_init%b_emit        = 0.0
  beam_init%dPz_dz        = 0.0
  beam_init%center(:)     = 0.0
  beam_init%bunch_charge  = 0.0
  beam_init%dt_bunch      = 0
  beam_init%sig_z         = 0.0
  beam_init%sig_e         = 0.0
  beam_init%renorm_center = .true.
  beam_init%renorm_sigma  = .true.
  beam_init%n_bunch       = 1
  beam_init%n_particle    = -1
  beam0_file = s%com%beam0_file        ! From the command line
  beam_all_file = s%com%beam_all_file  ! From the command line
  beam_saved_at = ''
  save_beam_at  = ''
  track_start = ''
  track_end = ''
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

  if (ios < 0 .and. ix_universe == -1) exit  ! Exit on end-of-file and no namelist read

  ! Error checking

  if (beam_init%n_bunch < 1) then
    call out_io (s_fatal$, r_name, 'BEAM_INIT%N_BUNCH NOT PROPERLY SET.')
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
integer i, j, ix, iu, n_part, ix_class, n_bunch, n_particle, n_loc
character(60) at, class, ele_name, line

! Set tracking start/stop

uni_branch0 => u%uni_branch(0)

uni_branch0%track_start = track_start
uni_branch0%track_end   = track_end

if (track_start /= '') then
  call lat_ele_locator (track_start, u%design%lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    call out_io (s_fatal$, r_name, 'TRACK_START ELEMENT NOT FOUND: ' // track_start)
    call err_exit
  endif
  if (n_loc > 1) then
    call out_io (s_fatal$, r_name, 'MULTIPLE TRACK_START ELEMENTS FOUND: ' // track_start)
    call err_exit
  endif
  uni_branch0%ix_track_start = eles(1)%ele%ix_ele
else
  uni_branch0%ix_track_start = ix_track_start
endif

if (track_end /= '') then
  call lat_ele_locator (track_end, u%design%lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    call out_io (s_fatal$, r_name, 'TRACK_END ELEMENT NOT FOUND: ' // track_end)
    call err_exit
  endif
  if (n_loc > 1) then
    call out_io (s_fatal$, r_name, 'MULTIPLE TRACK_END ELEMENTS FOUND: ' // track_end)
    call err_exit
  endif
  uni_branch0%ix_track_end = eles(1)%ele%ix_ele
else
  uni_branch0%ix_track_end = ix_track_end
endif

u%beam%beam_init = beam_init
u%beam%beam0_file = beam0_file
u%beam%beam_all_file = beam_all_file
call init_coord(u%design%lat_branch(0)%orbit(0), beam_init%center, &
                    u%design%lat%ele(0), downstream_end$, u%design%lat%param%particle)

! No initialization for a circular lattice

if (u%model%lat%param%geometry == closed$) return

! Find where to save the beam at.
! Always save at branch points.

do i = 0, ubound(u%model%lat%branch, 1)
  branch => u%design%lat%branch(i)
  u%uni_branch(i)%ele%save_beam = .false.
  u%uni_branch(i)%ele(0)%save_beam = .true.
  u%uni_branch(i)%ele(branch%n_ele_track)%save_beam = .true.
  do j = 1, ubound(branch%ele, 1)
    if (branch%ele(j)%key == fork$) u%uni_branch(i)%ele(j)%save_beam = .true.
  enddo
enddo

if (beam_saved_at /= '') then
  call tao_locate_elements (beam_saved_at, u%ix_uni, eles, err, ignore_blank = .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD BEAM_SAVED_AT ELEMENT: ' // beam_saved_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      u%uni_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam = .true.
    enddo
  endif
endif

u%beam%saved_at = beam_saved_at

! If beam_all_file is set, read in the beam distributions.

if (u%beam%beam_all_file /= '') then
  s%com%use_saved_beam_in_tracking = .true.
  call tao_open_beam_file (beam_all_file, err)
  if (err) call err_exit
  call tao_read_beam_file_header (j, u%beam%beam_init%n_bunch, &
                                             u%beam%beam_init%n_particle, err)  
  if (err) call err_exit
  do
    if (j == -1) exit
    call tao_read_beam (uni_branch0%ele(j)%beam, err)
    if (err) call err_exit
  enddo  
  call out_io (s_info$, r_name, 'Read beam_all file: ' // u%beam%beam_all_file)
  call tao_close_beam_file ()
endif

if (allocated(eles)) deallocate (eles)

end subroutine init_beam

end subroutine tao_init_beams


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_dynamic_aperture (init_file)
!
!-
subroutine tao_init_dynamic_aperture(init_file)
implicit none
type (tao_dynamic_aperture_init_struct)  :: da_init(200)
type (tao_universe_struct), pointer :: u

integer :: ios, iu, i, j, n_pz

character(*) init_file
character(200) file_name
character(40) :: r_name = 'tao_init_dynamic_aperture'

namelist / tao_dynamic_aperture / da_init

!

if (init_file == '') return

!call out_io (s_blank$, r_name, '*Init: Opening Init File: ' // file_name)
call tao_open_file (init_file, iu, file_name, s_blank$)
if (iu == 0) then
  call out_io (s_blank$, r_name, "Note: Cannot open init file for tao_dynamic_aperture namelist read")
  return
endif

! Read tao_dynamic_aperture
call out_io (s_blank$, r_name, 'Init: Reading tao_dynamic_aperture namelist')
read (iu, nml = tao_dynamic_aperture, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_DYNAMIC_APERTURE NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_dynamic_aperture)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_dynamic_aperture namelist found')

close(iu)

do i = lbound(s%u, 1), ubound(s%u, 1)
 ! Count the list of pz
  do n_pz=1, 200
    if (da_init(i)%pz(n_pz) == real_garbage$) exit
  enddo
  n_pz = n_pz - 1
  if (n_pz == 0 ) cycle
  
  ! Set 
  u => s%u(i)
  allocate(u%dynamic_aperture%scan(n_pz))
  allocate(u%dynamic_aperture%pz(n_pz))
  call out_io (s_blank$, r_name, 'Found n_pz: ', n_pz)
  u%dynamic_aperture%scan(:)%param = da_init(i)%param
  u%dynamic_aperture%scan(:)%min_angle = da_init(i)%min_angle
  u%dynamic_aperture%scan(:)%max_angle = da_init(i)%max_angle
  u%dynamic_aperture%scan(:)%n_angle = da_init(i)%n_angle
  u%dynamic_aperture%pz(1:n_pz) = da_init(i)%pz(1:n_pz)
  
enddo


end subroutine tao_init_dynamic_aperture


end module
