module lux_module

use random_mod
use photon_init_mod
use photon_target_mod
use track1_photon_mod
use em_field_mod
use bmad
use quick_plot
use input_mod

implicit none

! Break a bend or a wiggler into longitudinal slices. To save time, if the photons from a given slice will not 
! make it to the detector, do not emit photons from that slice.

type lux_bend_slice_struct
  type (ele_struct) ele
  type (em_field_struct) field
  type (coord_struct) orbit
  real(rp) integrated_emit_prob    ! Probability of photon emission from element start to slice end
  real(rp) emit_prob               ! Probability of photon emission withing slice
  logical good_emit                ! Only generate photons in slices with %good_emit = True. Photons generated 
                                   !   in slices with good_emit = False will not make it through the aperture
end type

! %bend_slice(i) are the emission parameters at the end of the i^th slice

type lux_param_struct
  character(100) :: photon1_out_file = ''
  character(100) :: det_pix_out_file = ''
  character(100) :: track_out_file = 'track.dat'
  character(100) :: param_file = 'lux.init'
  character(100) :: lattice_file = ''
  character(200) :: bmad_parameters_out = ''
  character(40) :: photon_init_element = ''       ! element name
  character(40) :: detector_element = ''          ! element name
  character(40) :: photon1_element = ''           ! element name
  character(20) :: plotting = ''
  character(16) :: random_engine = 'pseudo'
  character(16) :: histogram_variable = ''
  real(rp) :: histogram_bin_width = 0
  real(rp) :: intensity_min_det_pixel_cutoff = 1e-6
  real(rp) :: intensity_min_photon1_cutoff = 1e-6
  real(rp) :: stop_total_intensity = 10           ! stop intensity per energy
  real(rp) :: window_width = 800.0_rp, window_height = 400.0_rp  ! For plotting
  real(rp) :: intensity_normalization_coef = 1e6
  real(rp) :: stop_num_photons = 0               ! stop number. Use real for clearer input
  real(rp) :: mpi_run_size = 0.1                 ! Normalized number of photons to track
  real(rp) :: timer_print_dtime = 300
  integer :: random_seed = 0
  integer :: n_photon1_out_max = 100     ! Max number of photons allowed in photon1_out_file.
  logical :: debug = .false.             ! For debugging
  logical :: reject_dead_at_det_photon1 = .false.
  logical :: scale_initial_field_to_1 = .true.
  logical :: normalization_includes_pixel_area = .true.
end type

type lux_photon_struct
  integer n_photon_generated
  type (coord_struct), allocatable :: orb(:)
end type

type lux_common_struct
  type (lat_struct) :: lat
  type (branch_struct), pointer :: physical_source_branch, tracking_branch
  type (ele_struct), pointer :: physical_source_ele, photon_init_ele, fork_ele, detec_ele, photon1_ele
  type (pixel_detec_struct), pointer :: pixel
  type (bunch_params_struct), allocatable :: stat(:)
  integer n_bend_slice             ! Number of slices
  integer(8) :: n_photon_stop1 = 0 ! Number of photons to track per processor.
                                   ! This is equal to lux_param%stop_num_photons when not using mpi.
  integer :: mpi_rank  = -1
  integer :: mpi_n_proc = 1        ! Number of processeses including master
  integer :: n_photon1_out = 0     ! Count of photons in photon1_out_file
  integer :: iu_photon1_out        ! File I/O unit number.
  type (lux_bend_slice_struct), allocatable :: bend_slice(:) ! Size: (0:n_bend_slice)
  type (pixel_pt_struct), allocatable :: histogram_bin(:)
  real(rp) E_min, E_max                                      ! Photon energy range for bends and wigglers
  logical :: verbose = .false.
  logical :: using_mpi = .false.
end type

type lux_output_data_struct
  integer :: nx_min = 0, nx_max = 0, ny_min = 0, ny_max = 0
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_run_serial()
!
! Top level routine to run lux in serial (non-mpi) mode.
!-

subroutine lux_run_serial()

implicit none

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lat_struct), pointer :: lat

integer nx, ny

!------------------------------------------

call lux_init (lux_param, lux_com)
call lux_init_data (lux_param, lux_com)
call lux_track_photons (lux_param, lux_com)
if (lux_param%photon1_out_file /= '') close (lux_com%iu_photon1_out)
call lux_write_data (lux_param, lux_com)

end subroutine lux_run_serial

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_init (lux_param, lux_com)
!
! Routine to init Lux.
!
! Output:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!-

subroutine lux_init (lux_param, lux_com)

type (lux_param_struct) lux_param
type (lux_common_struct), target :: lux_com
type (lat_struct), pointer :: lat
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: photon_init_ele

real(rp) cut, r

integer i, j, n, ix, ie, ir
integer n_throw, n_loc, iu, rank_id

character(3) num_str
character(40) arg, plotting, number_file
character(100) photon1_out_file

logical ok, is_there, err

character(*), parameter :: r_name = 'lux_init'

namelist / params / lux_param

! Get inputs

lux_com%verbose = .true.

i = 0
do while (i < command_argument_count())
  i = i + 1
  call get_command_argument(i, arg)
  select case (arg)
  case ('-silent')
    lux_com%verbose = .false.
  case ('-plot')
    i = i + 1
    call get_command_argument (i, plotting)
  case default
    if (arg(1:1) == '-') then
      call out_io (s_fatal$, r_name, 'I DO NOT UNDERSTAND: ' // trim(arg))
      ok = .false.
    endif
    lux_param%param_file = arg
  end select
enddo

lux_param%intensity_normalization_coef = 1e6

open (1, file = lux_param%param_file, status = 'old', action = 'read')
read (1, nml = params)
close (1) 
 
if (lux_param%debug) then
  call out_io (s_info$, r_name, 'Note: lux_param%debug = True')
endif

if (.not. lux_com%verbose) call output_direct (-1, .false., max_level = s_success$) ! Suppress bmad_parser info output.

! Negative radom_seed is used for testing and turns off the dithering with MPI

call ran_seed_put (abs(lux_param%random_seed), mpi_offset = lux_com%mpi_rank)

if (lux_param%random_engine == 'quasi') then
  if (lux_param%random_seed == 0) then
    call ran_uniform(r)
    n_throw = 10000 * r
  else
    n_throw = modulo(lux_param%random_seed-1, 10000)
  endif
  call ran_engine(set = lux_param%random_engine)
  do i = 1, n_throw
    call ran_uniform(r)
  enddo
endif

! Add number to file name

if (lux_com%mpi_rank < 1 .and. (index(lux_param%photon1_out_file, '#') /= 0 .or. index(lux_param%det_pix_out_file, '#') /= 0)) then
  number_file = 'lux_out_file.number'
  inquire(file = number_file, exist = is_there)
  if (.not. is_there) then
    open (1, file = number_file)
    write (1, *) 0
    close (1)
  endif

  call increment_file_number (number_file, 3, n, num_str)
  call sub_in (lux_param%photon1_out_file)
  call sub_in (lux_param%det_pix_out_file)
  call sub_in (lux_param%track_out_file)

  open (1, file = 'lux_out_file.names')
  write (1, '(a)') '# Lux output file names for use with python post-processing scripts'
  write (1, '(3a)') 'photon1_out_file = "', trim(lux_param%photon1_out_file), '"'
  write (1, '(3a)') 'det_pix_out_file = "', trim(lux_param%det_pix_out_file), '"'
  write (1, '(3a)') 'track_out_file   = "', trim(lux_param%track_out_file), '"'
  close (1)
endif

! Read in lattice

bmad_com%auto_bookkeeper = .false.
call bmad_parser (lux_param%lattice_file, lux_com%lat)
lat => lux_com%lat

! locate photon_init element

call lat_ele_locator (lux_param%photon_init_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  call out_io (s_fatal$, r_name, 'NO PHOTON_INIT ELEMENT FOUND MATCHING NAME: "' // trim(lux_param%photon_init_element) // '"')
  stop
elseif (n_loc > 1) then
  call out_io (s_fatal$, r_name, 'MULTIPLE PHOTON_INIT ELEMENTS FOUND MATCHING NAME: "' // trim(lux_param%photon_init_element) // '"')
  stop
endif

lux_com%photon_init_ele => eles(1)%ele
photon_init_ele => lux_com%photon_init_ele

if (photon_init_ele%key /= photon_init$) then
  call out_io (s_fatal$, r_name, 'CANNOT SIMULATE PHOTONS GENERATED IN ELEMENT OF TYPE: ' // key_name(photon_init_ele%key))
  stop   
endif

! Is there an assocaited physical source element?

if (photon_init_ele%component_name == '') then
  nullify (lux_com%physical_source_ele, lux_com%physical_source_branch)
else
  call lat_ele_locator (photon_init_ele%component_name, lat, eles, n_loc, err)
  if (n_loc == 0) then
    call out_io (s_fatal$, r_name, 'PHYSICAL SOURCE ELEMENT ASSOCIATED WITH PHOTON_INIT ELEMENT NO FOUND: ' // photon_init_ele%component_name)
    stop
  elseif (n_loc > 1) then
    call out_io (s_fatal$, r_name, 'MULTIPLE PHYSICAL SOURCE ELEMENTS ASSOCIATED WITH PHOTON_INIT ELEMENT FOUND: ' // photon_init_ele%component_name)
    stop
  endif
  lux_com%physical_source_ele => eles(1)%ele
  lux_com%physical_source_branch => lux_com%physical_source_ele%branch
  lux_com%fork_ele => pointer_to_ele (lux_com%lat, photon_init_ele%branch%ix_from_ele, photon_init_ele%branch%ix_from_branch)
  if (lux_com%fork_ele%ix_branch /= lux_com%physical_source_ele%ix_branch) then
    call out_io (s_fatal$, r_name, 'FORK ELEMENT TO BRANCH CONTAINING THE PHOTON_INIT ELEMEMENT NOT THE SAME AS THE BRANCH OF THE PHYSICAL SOURCE ELEMENT!')
    stop
  endif
endif

! Locate detector element

call lat_ele_locator (lux_param%detector_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  call out_io (s_fatal$, r_name, 'NO DETECTOR ELEMENT FOUND MATCHING NAME: "' // trim(lux_param%detector_element) // '"')
  stop
elseif (n_loc > 1) then
  call out_io (s_fatal$, r_name, 'MULTIPLE DETECTOR ELEMENTS FOUND MATCHING NAME: "' // trim(lux_param%detector_element) // '"')
  stop
endif

lux_com%detec_ele => eles(1)%ele
lux_com%pixel => lux_com%detec_ele%photon%pixel
lux_com%tracking_branch => lux_com%detec_ele%branch
allocate(lux_com%stat(0:lux_com%tracking_branch%n_ele_track))

if (lux_com%detec_ele%ix_branch /= photon_init_ele%ix_branch) then
  call out_io (s_fatal$, r_name, 'PHOTON_INIT ELEMENT AND DETECTOR ELEMENT NOT IN THE SAME BRANCH!')
  stop
endif

! Locate photon1 element

if (lux_param%photon1_element == '') then
  lux_com%photon1_ele => lux_com%detec_ele
else
  call lat_ele_locator (lux_param%photon1_element, lat, eles, n_loc, err)
  if (n_loc == 0) then
    call out_io (s_fatal$, r_name, 'NO PHOTON1 ELEMENT FOUND MATCHING NAME: "' // trim(lux_param%photon1_element) // '"')
    stop
  elseif (n_loc > 1) then
    call out_io (s_fatal$, r_name, 'MULTIPLE PHOTON1 ELEMENTS FOUND MATCHING NAME: "' // trim(lux_param%photon1_element) // '"')
    stop
  endif
  lux_com%photon1_ele => eles(1)%ele 
endif

! Photon1 file init

if (lux_param%photon1_out_file /= '') then
  lux_com%n_photon1_out = 0
  iu = lunget()
  lux_com%iu_photon1_out = iu

  photon1_out_file = lux_param%photon1_out_file 
  ix = index(photon1_out_file, '@')
  if (ix /= 0) then
    rank_id = 0
    if (lux_com%using_mpi) rank_id = lux_com%mpi_rank
    write (photon1_out_file, '(a, i0, a)') photon1_out_file(1:ix-1), rank_id, photon1_out_file(ix+1:)
  endif
  open (iu, file = photon1_out_file, recl = 200)
  write (iu, '(a)') '#      |                               Start                            |                                    End                             |                 End'
  write (iu, '(a)') '#   Ix |     x (mm)           vx       y (mm)           vy            z |        x (mm)           vx       y (mm)           vy            z  |     Energy     Intens_x     Intens_y'
endif

! Histogram setup

if (lux_param%det_pix_out_file /= '' .and. lux_param%histogram_variable /= '') then
  if (lux_param%histogram_bin_width == 0) then
    call out_io (s_fatal$, r_name, 'LUX_PARAM%HISTOGRAM_BIN_WIDTH IS ZERO!')
    stop
  endif

  allocate(lux_com%histogram_bin(-2:2))
endif

! Tracking init

call run_timer('START')

if (lux_com%using_mpi) then
  lux_com%n_photon_stop1 = 0.1 + lux_param%stop_num_photons * lux_param%mpi_run_size / (lux_com%mpi_n_proc - 1)
else
  lux_com%n_photon_stop1 = lux_param%stop_num_photons
endif

call lux_tracking_setup (lux_param, lux_com)

if (lat%photon_type == coherent$) then
  if (photon_init_ele%value(e_field_x$) == 0 .and. photon_init_ele%value(e_field_y$) == 0) then
    call out_io (s_fatal$, r_name, 'WARNING: INPUT E_FIELD IS ZERO SO RANDOM FIELD WILL BE GENERATED WITH COHERENT PHOTONS!')
  endif
endif

!--------------------------------------------------------------------
contains

subroutine sub_in (file_name)

character(*) file_name
integer ix

ix = index(file_name, '#')
if (ix == 0) return
file_name = file_name(:ix-1) // num_str // trim(file_name(ix+1:))

end subroutine sub_in

end subroutine lux_init

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_init_data (lux_param, lux_com)
!
! Routine to initialize the output data.
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!-

subroutine lux_init_data (lux_param, lux_com)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param

character(*), parameter :: r_name = 'lux_init_data'

!

lux_com%pixel => lux_com%pixel

if (.not. allocated(lux_com%pixel%pt)) then
  call out_io (s_fatal$, r_name, 'DETECTOR GRID NOT SET!')
  stop
endif

lux_com%pixel%n_track_tot = 0
lux_com%pixel%n_hit_detec = 0
lux_com%pixel%n_hit_pixel = 0

lux_com%pixel%pt = pixel_pt_struct()

end subroutine lux_init_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_generate_photon (orb, lux_param, lux_com)
!
! Routine to generate the starting photon coordinates
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!
! Output:
!   orb         -- coord_struct: Initialized starting coords.
!-

subroutine lux_generate_photon (orb, lux_param, lux_com)

type (lat_struct), pointer :: lat
type (coord_struct) charged_orb
type (coord_struct) :: orb
type (lux_param_struct) lux_param
type (lux_common_struct), target :: lux_com
type (lux_bend_slice_struct), pointer :: sl(:)
type (ele_struct) ele_ave
type (ele_struct), pointer :: photon_init_ele, physical_source_ele

real(rp) x, y, phi, r(3), dir(2), ds, rr, r_emit(5), prob, f
real(rp) v_mat(4,4), v_inv_mat(4,4), vec(4), dE, e(3), b(3), sig_vec(6)
real(rp) g_bend(3), gamma_electron

integer ix, n_slice

! Init

photon_init_ele => lux_com%photon_init_ele
physical_source_ele => lux_com%physical_source_ele
lat => lux_com%lat

!-----------------------------------------------------
! Everything but the field handled by Bmad

if (.not. associated(lux_com%physical_source_ele)) then
  x = photon_init_ele%value(e_field_x$); y = photon_init_ele%value(e_field_y$)
  if (x == 0 .and. y == 0) then
    call ran_uniform(rr)
    orb%field(1) = cos(twopi * rr)
    orb%field(2) = sin(twopi * rr)
  else
    if (lux_param%scale_initial_field_to_1) then
      orb%field(1) = x / sqrt(x**2 + y**2)
      orb%field(2) = y / sqrt(x**2 + y**2)
    else
      orb%field(1) = x
      orb%field(2) = y
    endif
  endif
  call init_coord (orb, orb, photon_init_ele, upstream_end$, random_on = .true.)
  return
endif

!-----------------------------------------------------
! bend, wiggler, undulator source

sl => lux_com%bend_slice
n_slice = ubound(sl, 1)

! Find where photon emitted

call ran_uniform(rr)  ! longitudinal position
ix = bracket_index (rr, sl%integrated_emit_prob, 0)
ix = ix + 1
if (ix == n_slice) ix = n_slice - 1
f = (rr - sl(ix-1)%integrated_emit_prob) / (sl(ix)%integrated_emit_prob - sl(ix-1)%integrated_emit_prob)

! Calculate electron average position

charged_orb = sl(ix)%orbit
charged_orb%vec = (1-f) * sl(ix)%orbit%vec + f * sl(ix+1)%orbit%vec

ele_ave%a = average_twiss(1-f, sl(ix)%ele%a, sl(ix+1)%ele%a)
ele_ave%b = average_twiss(1-f, sl(ix)%ele%b, sl(ix+1)%ele%b)
ele_ave%c_mat   = (1-f) * sl(ix)%ele%c_mat   + f * sl(ix+1)%ele%c_mat
ele_ave%gamma_c = (1-f) * sl(ix)%ele%gamma_c + f * sl(ix+1)%ele%gamma_c
call make_v_mats (ele_ave, v_mat, v_inv_mat)

! Add offsets due to finite bunch size to the electron position.
! To do this must transform to the normal mode coords

call ran_gauss(r_emit)  ! electron momentum offset.

charged_orb%vec(6) = charged_orb%vec(6) + r_emit(5) * lat%z%sigmap

dE = charged_orb%vec(6)
vec = matmul (v_inv_mat, charged_orb%vec(1:4))
vec(1:2) = vec(1:2) + charged_orb%vec(1:2) + [ele_ave%a%eta, ele_ave%a%etap] * dE
vec(3:4) = vec(3:4) + charged_orb%vec(1:2) + [ele_ave%b%eta, ele_ave%b%etap] * dE

vec(1) = vec(1) + sqrt(lat%a%emit * ele_ave%a%beta) * r_emit(1)
vec(2) = vec(2) + sqrt(lat%a%emit / ele_ave%a%beta) * (r_emit(2) - ele_ave%a%alpha * r_emit(1))

vec(3) = vec(3) + sqrt(lat%b%emit * ele_ave%b%beta) * r_emit(3)
vec(4) = vec(4) + sqrt(lat%b%emit / ele_ave%b%beta) * (r_emit(4) - ele_ave%b%alpha * r_emit(3))

charged_orb%vec(1:4) = matmul(v_mat, vec)

! Calculate bending strength

B = (1-f) * sl(ix)%field%b + f * sl(ix+1)%field%b
E = 0
g_bend = g_bend_from_em_field (B, E, charged_orb)

! Init photon

gamma_electron = physical_source_ele%value(p0c$) * &
                                (1 + sl(ix)%orbit%vec(6)) / sl(ix)%orbit%beta / mass_of(sl(ix)%orbit%species)
!! Note: Energy slices for coherent photons is not yet implemented.
!! if (lux_com%lat%photon_type == coherent$ .and. ix_energy > 0) then
!!   rr = (ix_energy - 0.5_rp) / lux_param%n_energy_pts
!!   call bend_photon_init (g_bend(1), g_bend(2), gamma_electron, orb, lux_com%E_min, lux_com%E_max, rr)
!! endif

call bend_photon_init (g_bend(1), g_bend(2), gamma_electron, orb, lux_com%E_min, lux_com%E_max)
call absolute_photon_position (charged_orb, orb)

orb%s = sl(ix-1)%ele%s + f * sl(ix)%ele%value(l$)

! Track to fork element.

ds = lux_com%fork_ele%s - orb%s  

if (physical_source_ele%key == sbend$) then
  call track_a_bend_photon (orb, physical_source_ele, ds)
else
  call track_a_drift_photon (orb, ds, .true.)
endif

end subroutine lux_generate_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_tracking_setup (lux_param, lux_com)
!
! Routine 
!
! Input:
!   lux_param     -- lux_param_struct: Lux input parameters.
!   lux_com       -- lux_common_struct: Common parameters.
!
! Output:
!   lux_com       -- Lux common block.
!-

subroutine lux_tracking_setup (lux_param, lux_com)

type (lat_struct), pointer :: lat
type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (floor_position_struct) floor
type (ele_struct) twiss_ele
type (ele_struct), pointer :: ele, fork_ele, photon_init_ele, physical_source_ele
type (lux_bend_slice_struct), pointer :: sl(:)
type (coord_struct) orb
type (coord_struct), pointer :: orbit
type (branch_struct), pointer :: branch

type coord4_struct
  real(rp) vec(4)
end type
type (coord4_struct) p_coord(4)

real(rp) vz, rho, x, f
real(rp) phi, y, ds, s_now, prob1, prob2, g_bend(3), g_abs
real(rp) gamma, v_mat(4,4)

integer i, j, k, n, n_phi, n_y, ip, iy, iz, ie, iz2
integer track_state, iy0, iy1, ip0, ip1, iz0, iz1
integer n_slice, n_z, ix

logical err, hit_below_top, hit_above_bottom, hit_right_of_left_edge, hit_left_of_right_edge
logical old_hit_below_top, old_hit_above_bottom, old_hit_right_of_left_edge, old_hit_left_of_right_edge

character(*), parameter :: r_name = 'lux_tracking_setup'

!-------------------------------------------------------------
! photon_init source

lat => lux_com%lat
branch => lux_com%tracking_branch
photon_init_ele => lux_com%photon_init_ele
physical_source_ele => lux_com%physical_source_ele

if (.not. associated(lux_com%physical_source_ele)) then
  call photon_target_setup (photon_init_ele)

!-------------------------------------------------------------
! Sbend or wiggler source

else
  call photon_target_setup (lux_com%fork_ele)

  if (photon_init_ele%value(ds_slice$) == 0) then
    call out_io (s_fatal$, r_name, 'DS_SLICE IS ZERO FOR ELEMENT: ' // photon_init_ele%name)
    stop
  endif

  n_slice = max(1, nint(physical_source_ele%value(l$) / photon_init_ele%value(ds_slice$)))
  lux_com%n_bend_slice = n_slice
  allocate (lux_com%bend_slice(0:n_slice))

  lux_com%E_min = photon_init_ele%value(E_center$) - photon_init_ele%value(sig_E$)
  if (is_true(photon_init_ele%value(E_center_relative_to_ref$))) lux_com%E_min = lux_com%E_min + lux_com%detec_ele%value(p0c$) 
  if (photon_init_ele%value(sig_E$) == 0) then
    lux_com%E_max = lux_com%E_min + 1d-10  ! Need some small offset for the calculation
  else
    lux_com%E_max = lux_com%E_min + 2 * photon_init_ele%value(sig_E$)
  endif

  ! Track through physical source ele and gather data

  call init_coord (orb, lat%particle_start, physical_source_ele, upstream_end$)
  twiss_ele = pointer_to_next_ele (physical_source_ele, -1)
  ds = physical_source_ele%value(l$) / n_slice
  s_now = 0

  gamma = (orb%p0c / orb%beta) / mass_of(orb%species)
  sl => lux_com%bend_slice

  old_hit_below_top          = .false.
  old_hit_above_bottom       = .false.
  old_hit_left_of_right_edge = .false.
  old_hit_right_of_left_edge = .false.

  do i = 0, n_slice

    call transfer_ele (twiss_ele, sl(i)%ele)
    sl(i)%orbit   = orb
    call em_field_calc (physical_source_ele, lat%param, s_now, orb, .false., sl(i)%field)

    g_bend = g_bend_from_em_field (sl(i)%field%b, sl(i)%field%e, orb)
    g_abs = norm2(g_bend)
    prob1 = bend_photon_energy_integ_prob (lux_com%E_min, g_abs, gamma)  ! Probability per radian of bend
    prob2 = bend_photon_energy_integ_prob (lux_com%E_max, g_abs, gamma)
    sl(i)%emit_prob = g_abs * ds * (prob2 - prob1)      ! Probability per slice

    ! See if any photons will make it hit inside the aperture. If so, set %good_emit = True.
    ! Orientation: +x = left, +y = up

    hit_below_top = .false.
    hit_above_bottom = .false.
    hit_right_of_left_edge = .false.
    hit_left_of_right_edge = .false.

    call make_v_mats (twiss_ele, v_mat)
    p_coord(1)%vec = matmul(v_mat, [1.0_rp, -twiss_ele%a%alpha, 0.0_rp, 0.0_rp]) * &
                              sqrt(lux_com%physical_source_branch%a%emit * twiss_ele%a%beta) * photon_init_ele%value(transverse_sigma_cut$)
    p_coord(2)%vec = matmul(v_mat, [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]) * &
                              sqrt(lux_com%physical_source_branch%a%emit / twiss_ele%a%beta) * photon_init_ele%value(transverse_sigma_cut$)
    p_coord(3)%vec = matmul(v_mat, [0.0_rp, 0.0_rp, 1.0_rp, -twiss_ele%b%alpha]) * &
                              sqrt(lux_com%physical_source_branch%b%emit * twiss_ele%b%beta) * photon_init_ele%value(transverse_sigma_cut$)
    p_coord(4)%vec = matmul(v_mat, [0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) * &
                              sqrt(lux_com%physical_source_branch%b%emit / twiss_ele%b%beta) * photon_init_ele%value(transverse_sigma_cut$)

    fork_ele => lux_com%fork_ele
    do k = 1, fork_ele%photon%target%n_corner
      floor = coords_relative_to_floor (fork_ele%floor, fork_ele%photon%target%corner(k)%r)
      floor = coords_floor_to_relative (twiss_ele%floor, floor, .false.)
      do j = 1, 4
        x = p_coord(j)%vec(1) + p_coord(j)%vec(2) * floor%r(3)
        if ( x > floor%r(1)) hit_left_of_right_edge = .true.
        if (-x < floor%r(1)) hit_right_of_left_edge = .true.
        y = p_coord(j)%vec(3) + p_coord(j)%vec(4) * floor%r(3)
        if ( y > floor%r(2)) hit_above_bottom = .true.
        if (-y < floor%r(2)) hit_below_top    = .true.
      enddo
    enddo 

    sl(i)%good_emit = ((old_hit_below_top .or. hit_below_top) .and. (old_hit_above_bottom .or. hit_above_bottom) .and. &
                       (old_hit_left_of_right_edge .or. hit_left_of_right_edge) .and. (old_hit_right_of_left_edge .or. hit_right_of_left_edge))

    old_hit_below_top          = hit_below_top
    old_hit_above_bottom       = hit_above_bottom
    old_hit_left_of_right_edge = hit_left_of_right_edge
    old_hit_right_of_left_edge = hit_right_of_left_edge

    !

    call twiss_and_track_intra_ele (physical_source_ele, lat%param, s_now, s_now+ds, &
                                      .true., .true., orb, orb, twiss_ele, twiss_ele, err, .true.)
    if (err) call err_exit
    s_now = s_now + ds

  enddo

  sl(0)%integrated_emit_prob = 0
  do i = 1, n_slice
    if (sl(i)%good_emit) then
      sl(i)%integrated_emit_prob = sl(i)%emit_prob + sl(i-1)%integrated_emit_prob
    else
      sl(i)%integrated_emit_prob = sl(i-1)%integrated_emit_prob
    endif
  enddo

  if (sl(n_slice)%integrated_emit_prob == 0) then
    call out_io (s_fatal$, r_name, 'INTEGRATED PHOTON EMISSION PROBABILITY IS ZERO! NO BENDING FIELD?')
    stop
  endif

  sl%integrated_emit_prob = sl%integrated_emit_prob / sl(n_slice)%integrated_emit_prob

  if (lux_com%verbose) print '(a, i4)', &
            'Number of slices of physical_source element to be used for photon generation:', count(sl%good_emit) 

endif

end subroutine lux_tracking_setup 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_track_photons (lux_param, lux_com)
!
! Routine to track photons.
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!-

subroutine lux_track_photons (lux_param, lux_com)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lux_photon_struct), target :: photon
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: detec_ele, photon_init_ele, photon1_ele
type (branch_struct), pointer :: s_branch, t_branch
type (bunch_struct) bunch, bunch_stop1, bunch_start
type (pixel_detec_struct), pointer :: detec_pixel
type (pixel_pt_struct), pointer :: pix
type (pixel_pt_struct) :: pixel

real(rp) phase, intensity_tot, intens, time0, time

integer i, ix, n, nt, nx, ny, track_state, ip, ie

logical err_flag

character(*), parameter :: r_name = 'lux_track_photons'

!

lat => lux_com%lat
detec_ele => lux_com%detec_ele
detec_pixel => lux_com%pixel
photon_init_ele => lux_com%photon_init_ele
photon1_ele => lux_com%photon1_ele
t_branch => lux_com%tracking_branch
s_branch => lux_com%physical_source_branch

nt = detec_ele%ix_ele

call reallocate_coord (photon%orb, lat, t_branch%ix_branch)

!

intensity_tot = 0
time0 = 0

do
  if (lux_com%n_photon_stop1 > 0 .and. detec_pixel%n_track_tot >= lux_com%n_photon_stop1) exit
  if (.not. lux_com%using_mpi .and. lux_param%stop_total_intensity > 0 .and. &
                                             intensity_tot >= lux_param%stop_total_intensity) exit
  detec_pixel%n_track_tot = detec_pixel%n_track_tot + 1
  photon%n_photon_generated = detec_pixel%n_track_tot

  call lux_generate_photon (photon%orb(0), lux_param, lux_com)
  !if (lux_param%debug) then
  !  call init_coord (photon%orb(0), lat%particle_start, t_branch%ele(0), downstream_end$, photon$, &
  !                                      1, t_branch%ele(0)%value(E_tot$) * (1 + lat%particle_start%vec(6)))
  !endif

  call track_all (lat, photon%orb, t_branch%ix_branch, track_state)

  if (.not. lux_com%using_mpi) then
    call run_timer('READ', time)
    if (time > time0 + lux_param%timer_print_dtime) then
      time0 = time
      print '(a, f10.2, a, i0, a, es16.5)', 'Time (min): ', time/60, &
                    ': N_track_tot: ', detec_pixel%n_track_tot, '  Intensity_tot:', intensity_tot
    endif
  endif

  call add_this_to_detector_statistics (photon%orb(0), photon%orb(nt), intens, intensity_tot)
  do i = 1, nt
    if (photon%orb(i)%state == alive$) lux_com%stat(i)%n_particle_live = lux_com%stat(i)%n_particle_live + 1
  enddo

  ! Write to photon1_out_file

  if (lux_param%photon1_out_file /= '') then
    call photon1_out (photon%orb(1), photon%orb(photon1_ele%ix_ele), photon%orb(detec_ele%ix_ele), int(detec_pixel%n_track_tot))
  endif
enddo

!--------------------------------------------------------------------
contains

subroutine photon1_out (orb_start, orb1, orb_det, ix_photon)

type (coord_struct) orb_start, orb1, orb_det
real(rp) v_start(5), v1(5)
integer ix_photon

!

if (lux_param%reject_dead_at_det_photon1 .and. orb_det%state /= alive$) return
if (orb1%state /= alive$) return
if (intens < lux_param%intensity_min_photon1_cutoff) return
lux_com%n_photon1_out = lux_com%n_photon1_out + 1
if (lux_com%n_photon1_out > lux_param%n_photon1_out_max) return

v_start = orb_start%vec(1:5); v_start(1) = 1d3* v_start(1); v_start(3) = 1d3* v_start(3)
v1 = orb1%vec(1:5); v1(1) = 1d3* v1(1); v1(3) = 1d3* v1(3)

write (lux_com%iu_photon1_out, '(i6, 2(5f13.6, 3x), f11.3, 2es13.4)') ix_photon, &
                        v_start, v1, orb1%p0c, orb_start%field(1)**2, orb1%field(2)**2

end subroutine photon1_out

!--------------------------------------------------------------------
! contains

subroutine add_this_to_detector_statistics (start_orb, det_orb, intens, intensity_tot)

type (coord_struct) start_orb, det_orb, orb0, orb, surf_orb
type (pixel_pt_struct), allocatable :: temp_bin(:)
real(rp) intens, intensity_tot
integer ix, i1, i2

!

if (det_orb%state /= alive$) return

call to_surface_coords (det_orb, detec_ele, surf_orb)

if (surf_orb%state /= alive$) then
  lux_com%pixel%n_hit_detec = lux_com%pixel%n_hit_detec + 1
  return
endif

intens = surf_orb%field(1)**2 + surf_orb%field(2)**2
intensity_tot = intensity_tot + intens

call photon_add_to_detector_statistics (start_orb, surf_orb, detec_ele)

! Histogram

if (lux_param%det_pix_out_file /= '' .and. lux_param%histogram_variable /= '') then
  orb  = to_photon_angle_coords (surf_orb, detec_ele)
  orb0 = to_photon_angle_coords (start_orb, detec_ele)

  select case (lux_param%histogram_variable)
  case ('init_x')
    ix = nint(orb0%vec(1) / lux_param%histogram_bin_width)
  case ('init_y')
    ix = nint(orb0%vec(3) / lux_param%histogram_bin_width)
  case ('init_x_angle')
    ix = nint(orb0%vec(2) / lux_param%histogram_bin_width)
  case ('init_y_angle')
    ix = nint(orb0%vec(4) / lux_param%histogram_bin_width)
  case ('x_angle')
    ix = nint(orb%vec(2) / lux_param%histogram_bin_width)
  case ('y_angle')
    ix = nint(orb%vec(4) / lux_param%histogram_bin_width)
  case ('energy')
    ix = nint(orb%vec(6) / lux_param%histogram_bin_width)
  case default
    call out_io (s_fatal$, r_name, 'BAD HISTOGRAM_VARIABLE SETTING: ' // lux_param%histogram_variable, &
                                   'SHOULD BE ONE OF: "init_x", "y_angle", "energy", etc.')
    stop
  end select

  if (ix < lbound(lux_com%histogram_bin,1) .or. ix > ubound(lux_com%histogram_bin,1)) then
    i1 = min(ix, lbound(lux_com%histogram_bin,1)); i2 = max(ix, ubound(lux_com%histogram_bin,1))
    call move_alloc(lux_com%histogram_bin, temp_bin)
    allocate (lux_com%histogram_bin(i1:i2))
    i1 = lbound(temp_bin,1); i2 = ubound(temp_bin,1)
    lux_com%histogram_bin(i1:i2) = temp_bin
  endif

  call photon_add_to_detector_statistics (start_orb, surf_orb, detec_ele, pixel_pt = lux_com%histogram_bin(ix))
endif

end subroutine add_this_to_detector_statistics

end subroutine lux_track_photons


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_add_in_slave_data (slave_pt, lux_param, lux_com)
!
! Routine to combine data from an MPI slave to the master data structure.
!
! Input:
!   slave_pixel   -- pixel_detec_struct: Slave data.
!   lux_param     -- lux_param_struct: Lux input parameters.
!   lux_com       -- lux_common_struct: Common parameters.
!
! Output:
!   lux_com%lat%pixel%pt
!-

subroutine lux_add_in_mpi_slave_data (slave_pixel, lux_param, lux_com)

type (lux_param_struct) lux_param
type (lux_common_struct), target :: lux_com
type (pixel_detec_struct) slave_pixel
integer i, j

!

do i = lbound(lux_com%pixel%pt, 1), ubound(lux_com%pixel%pt, 1)
do j = lbound(lux_com%pixel%pt, 2), ubound(lux_com%pixel%pt, 2)
  lux_com%pixel%pt(i,j)%n_photon       = lux_com%pixel%pt(i,j)%n_photon       + slave_pixel%pt(i,j)%n_photon
  lux_com%pixel%pt(i,j)%E_x            = lux_com%pixel%pt(i,j)%E_x            + slave_pixel%pt(i,j)%E_x
  lux_com%pixel%pt(i,j)%E_y            = lux_com%pixel%pt(i,j)%E_y            + slave_pixel%pt(i,j)%E_y
  lux_com%pixel%pt(i,j)%intensity_x    = lux_com%pixel%pt(i,j)%intensity_x    + slave_pixel%pt(i,j)%intensity_x
  lux_com%pixel%pt(i,j)%intensity_y    = lux_com%pixel%pt(i,j)%intensity_y    + slave_pixel%pt(i,j)%intensity_y
  lux_com%pixel%pt(i,j)%intensity      = lux_com%pixel%pt(i,j)%intensity      + slave_pixel%pt(i,j)%intensity
  lux_com%pixel%pt(i,j)%orbit          = lux_com%pixel%pt(i,j)%orbit          + slave_pixel%pt(i,j)%orbit
  lux_com%pixel%pt(i,j)%orbit_rms      = lux_com%pixel%pt(i,j)%orbit_rms      + slave_pixel%pt(i,j)%orbit_rms
  lux_com%pixel%pt(i,j)%init_orbit     = lux_com%pixel%pt(i,j)%init_orbit     + slave_pixel%pt(i,j)%init_orbit
  lux_com%pixel%pt(i,j)%init_orbit_rms = lux_com%pixel%pt(i,j)%init_orbit_rms + slave_pixel%pt(i,j)%init_orbit_rms
enddo
enddo

lux_com%pixel%n_track_tot   = lux_com%pixel%n_track_tot   + slave_pixel%n_track_tot
lux_com%pixel%n_hit_detec   = lux_com%pixel%n_hit_detec   + slave_pixel%n_hit_detec
lux_com%pixel%n_hit_pixel   = lux_com%pixel%n_hit_pixel   + slave_pixel%n_hit_pixel

end subroutine lux_add_in_mpi_slave_data 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_write_data (lux_param, lux_com)
!
! Routine to write the photon tracking data.
!
! Output:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!-

subroutine lux_write_data (lux_param, lux_com)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (pixel_detec_struct), pointer :: pixel
type (pixel_pt_struct), pointer :: pix
type (pixel_pt_struct) :: p
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (all_pointer_struct), allocatable :: ptr_array(:)

real(rp) normalization, cut, dtime, nrm2
real(rp) total_dead_intens, pix_in_file_intens
real(rp) :: intens_tot, intens_tot_x, intens_tot_y, intens_max
real(rp) x_sum, y_sum, x2_sum, y2_sum, x_ave, y_ave, e_rms, e_ave
real(rp) phase_x, phase_y, x, y

integer i, j, ix, ix2, nx, ny, n_min, n_max, n
integer nx_active_min, nx_active_max, ny_active_min, ny_active_max

logical err_flag

character(200) params_out, param

!------------------------------------------

pixel => lux_com%pixel
lat => lux_com%lat

do nx = lbound(pixel%pt, 1), ubound(pixel%pt, 1)
do ny = lbound(pixel%pt, 2), ubound(pixel%pt, 2)
  pix => pixel%pt(nx,ny)

  if (lat%photon_type == coherent$) then
    pix%intensity_x = abs(pix%E_x)
    pix%intensity_y = abs(pix%E_y)
    pix%intensity   = pix%intensity_x + pix%intensity_y 
  endif

  if (pix%intensity == 0) cycle
  pix%orbit          = pix%orbit / pix%intensity
  pix%orbit_rms      = pix%orbit_rms / pix%intensity
  pix%init_orbit     = pix%init_orbit / pix%intensity
  pix%init_orbit_rms = pix%init_orbit_rms / pix%intensity
enddo
enddo

!

do n = 1, 6
  intens_tot = sum(pixel%pt%intensity)
  p%orbit(n) = sum(pixel%pt%orbit(n)) / intens_tot
  p%orbit_rms(n) = sqrt(max(0.0_rp, sum(pixel%pt%orbit_rms(n)) / intens_tot - p%orbit(n)**2))
enddo

!------------------------------------------
! lux_param%det_pix_out_file

pix_in_file_intens = 0
intens_tot_x = 0
intens_tot_y = 0
intens_tot = 0

normalization = lux_param%intensity_normalization_coef / pixel%n_track_tot 
if (lux_param%normalization_includes_pixel_area) then
  normalization = normalization / (pixel%dr(1) * pixel%dr(2))
endif

if (lux_param%det_pix_out_file /= '') then
  open (3, file = lux_param%det_pix_out_file, recl = 300)
  open (4, file = trim(lux_param%det_pix_out_file) // '.x', recl = 240)
  open (5, file = trim(lux_param%det_pix_out_file) // '.y', recl = 240)

  intens_max = maxval(pixel%pt%intensity)
  cut = max(0.0_rp, intens_max * lux_param%intensity_min_det_pixel_cutoff)

  nx_active_min = ubound(pixel%pt, 1);  nx_active_max = lbound(pixel%pt, 1)
  ny_active_min = ubound(pixel%pt, 2);  ny_active_max = lbound(pixel%pt, 2)

  do i = lbound(pixel%pt, 1), ubound(pixel%pt, 1)
  do j = lbound(pixel%pt, 2), ubound(pixel%pt, 2)
    if (pixel%pt(i,j)%intensity <= cut) cycle
    nx_active_min = min(nx_active_min, i);  nx_active_max = max(nx_active_max, i)
    ny_active_min = min(ny_active_min, j);  ny_active_max = max(ny_active_max, j)
  enddo
  enddo

  call lux_write_header ('master_parameter_file             = ' // quote(lux_param%param_file))
  call lux_write_header ('lattice_file                      = ' // quote(lux_param%lattice_file))
  call lux_write_header ('intensity_normalization_coef      = ', 'es12.4', re = lux_param%intensity_normalization_coef)
  call lux_write_header ('normalization_includes_pixel_area = ', 'l2', logic = lux_param%normalization_includes_pixel_area)
  call lux_write_header ('normalization       = ', 'es14.6', re = normalization)
  call lux_write_header ('intensity_x_unnorm  = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity_x))
  call lux_write_header ('intensity_x_norm    = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity_x) * normalization)
  call lux_write_header ('intensity_y_unnorm  = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity_y))
  call lux_write_header ('intensity_y_norm    = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity_y) * normalization)
  call lux_write_header ('intensity_unnorm    = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity))
  call lux_write_header ('intensity_norm      = ', 'es16.5', re = sum(pixel%pt(:,:)%intensity) * normalization)
  call lux_write_header ('n_track_tot         = ', 'i0',    int8 = pixel%n_track_tot)
  call lux_write_header ('n_at_detec          = ', 'i0',    int8 = pixel%n_hit_detec)
  call lux_write_header ('n_hit_pixel         = ', 'i0',    int8 = pixel%n_hit_pixel)
  call lux_write_header ('dx_pixel            = ', 'f10.6',  re = pixel%dr(1))
  call lux_write_header ('dy_pixel            = ', 'f10.6',  re = pixel%dr(2))
  call lux_write_header ('nx_active_min       = ', 'i8',    int = nx_active_min)
  call lux_write_header ('nx_active_max       = ', 'i8',    int = nx_active_max)
  call lux_write_header ('ny_active_min       = ', 'i8',    int = ny_active_min)
  call lux_write_header ('ny_active_max       = ', 'i8',    int = ny_active_max)

  params_out = lux_param%bmad_parameters_out
  param_loop: do 
    if (params_out == '') exit
    ix = index(params_out, ';')
    if (ix == 0) then
      param = upcase(params_out)
      params_out = ''
    else
      param = upcase(params_out(1:ix-1))
      params_out = params_out(ix+1:)
    endif

    ix = index(param, '[')
    if (ix == 0) then  ! Must be a constant
      if (allocated(lat%constant)) then
        do i = 1, size(lat%constant)
          if (param /= lat%constant(i)%name) cycle
          write (3, '(2a, es16.8)') trim(downcase(param)), ' = ', lat%constant(i)%value
          cycle param_loop
        enddo
      endif

    else
      ix2 = index(param, ']')
      if (ix2 == 0 .or. param(ix2+1:) /= '') then
        write (3, '(2a, es16.8)') trim(downcase(param)), ' = "UNKNOWN"'
        cycle
      endif
      call pointers_to_attribute (lat, param(1:ix-1), param(ix+1:ix2-1), .true., ptr_array, err_flag)
      if (size(ptr_array) /= 0) then
        param = param(1:ix-1) // '__' // param(ix+1:ix2-1)
        write (3, '(2a, es16.8)') trim(downcase(param)), ' = ', ptr_array(1)%r
        cycle
      endif
    endif

    write (3, '(2a, es16.8)') trim(downcase(param)), ' = "UNKNOWN"'

  enddo param_loop

  write (3, '(a, es16.5, 3a)') '# x_center_det        =', p%orbit(1),     '  # ', to_str(p%orbit(1) / pixel%dr(1), 4), ' pixels'
  write (3, '(a, es16.5, 3a)') '# y_center_det        =', p%orbit(3),     '  # ', to_str(p%orbit(3) / pixel%dr(2), 4), ' pixels'
  write (3, '(a, es16.5, 3a)') '# x_rms_det           =', p%orbit_rms(1), '  # ', to_str(p%orbit_rms(1) / pixel%dr(1), 4), ' pixels'
  write (3, '(a, es16.5, 3a)') '# y_rms_det           =', p%orbit_rms(3), '  # ', to_str(p%orbit_rms(3) / pixel%dr(2), 4), ' pixels'
  write (3, '(a)') '#-----------------------------------------------------'
  write (3, '(a)') '#                                                                                                                                                                    |                                          Init'
  write (3, '(a)') '#   ix    iy      x_pix      y_pix   Intens_x    Phase_x   Intens_y    Phase_y  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave  Ang_y_rms|     X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms  Ang_y_ave  Ang_y_rms'

  do i = lbound(pixel%pt, 1), ubound(pixel%pt, 1)
  do j = lbound(pixel%pt, 2), ubound(pixel%pt, 2)
    pix => pixel%pt(i,j)
    intens_tot = intens_tot + pix%intensity 
    if (pix%intensity <= cut .or. pix%n_photon == 0) cycle
    phase_x = 0; phase_y = 0

    if (lat%photon_type == coherent$) then
      phase_x = atan2(aimag(pix%e_x), real(pix%e_x))
      phase_y = atan2(aimag(pix%e_y), real(pix%e_y))
    endif

    p = pix
    do n = 1, 6
      p%orbit_rms(n) = sqrt(max(0.0_rp, pix%orbit_rms(n) - pix%orbit(n)**2))
      p%init_orbit_rms(n) = sqrt(max(0.0_rp, pix%init_orbit_rms(n) - pix%init_orbit(n)**2))
    enddo

    write (3, '(2i6, 2es11.3, 5es11.3, i12, 2f10.3, 12es11.3)') i, j, [i,j]*pixel%dr+pixel%r0, &
           p%intensity_x * normalization, phase_x, p%intensity_y * normalization, phase_y, &
           p%intensity * normalization, p%n_photon, p%orbit(6), p%orbit_rms(6), &
           p%orbit(2), p%orbit_rms(2), p%orbit(4), p%orbit_rms(4), &
           p%init_orbit(1), p%init_orbit_rms(1), p%init_orbit(2), p%init_orbit_rms(2), p%init_orbit(3), p%init_orbit_rms(3), p%init_orbit(4), p%init_orbit_rms(4)
    pix_in_file_intens = pix_in_file_intens + p%intensity 
  enddo
  enddo
  close(3)

  ! det_pix_out_file.x

  write (4, '(a)') '#-----------------------------------------------------'
  write (4, '(a)') '#                                                                                                                             |                    Init'
  write (4, '(a)') '#   ix      x_pix   Intens_x   Intens_y  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave  Ang_y_rms|     X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms  Ang_y_ave  Ang_y_rms'

  do i = lbound(pixel%pt, 1), ubound(pixel%pt, 1)
    p = pixel_pt_struct()
    p%intensity_x = sum(pixel%pt(i,:)%intensity_x)
    p%intensity_y = sum(pixel%pt(i,:)%intensity_y)
    p%intensity   = sum(pixel%pt(i,:)%intensity)
    p%n_photon    = sum(pixel%pt(i,:)%n_photon)

    if (p%intensity /= 0) then
      do n = 1, 6
        p%orbit(n) = sum(pixel%pt(i,:)%orbit(n) * pixel%pt(i,:)%intensity) / p%intensity
        p%orbit_rms(n) = sqrt(max(0.0_rp, sum(pixel%pt(i,:)%orbit_rms(n) * pixel%pt(i,:)%intensity) / p%intensity - p%orbit(n)**2))
        p%init_orbit(n) = sum(pixel%pt(i,:)%init_orbit(n) * pixel%pt(i,:)%intensity) / p%intensity
        p%init_orbit_rms(n) = sqrt(max(0.0_rp, sum(pixel%pt(i,:)%init_orbit_rms(n) * pixel%pt(i,:)%intensity) / p%intensity - p%init_orbit(n)**2))
      enddo
    endif

    write (4, '(i6, es11.3, 3es11.3, i12, 2f10.3, 12es11.3)') i, i*pixel%dr(1)+pixel%r0(1), &
                       p%intensity_x * normalization, p%intensity_y * normalization, p%intensity * normalization, &
                       p%n_photon, p%orbit(6), p%orbit_rms(6), p%orbit(2), p%orbit_rms(2), p%orbit(4), p%orbit_rms(4), &
                       p%init_orbit(1), p%init_orbit_rms(1), p%init_orbit(2), p%init_orbit_rms(2), p%init_orbit(3), p%init_orbit_rms(3), p%init_orbit(4), p%init_orbit_rms(4)
  enddo
  close(4)

  ! det_pix_out_file.y

  write (5, '(a)') '#-----------------------------------------------------'
  write (5, '(a)') '#                                                                                                                             |                    Init'
  write (5, '(a)') '#   iy      y_pix   Intens_x   Intens_y  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave  Ang_y_rms|     X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms  Ang_y_ave  Ang_y_rms'
  do j = lbound(pixel%pt, 2), ubound(pixel%pt, 2)
    p = pixel_pt_struct()
    p%intensity_x = sum(pixel%pt(:,j)%intensity_x)
    p%intensity_y = sum(pixel%pt(:,j)%intensity_y)
    p%intensity   = sum(pixel%pt(:,j)%intensity)
    p%n_photon    = sum(pixel%pt(:,j)%n_photon)

    if (p%intensity /= 0) then
      do n = 1, 6
        p%orbit(n) = sum(pixel%pt(:,j)%orbit(n) * pixel%pt(:,j)%intensity) / p%intensity
        p%orbit_rms(n) = sqrt(max(0.0_rp, sum(pixel%pt(:,j)%orbit_rms(n) * pixel%pt(:,j)%intensity) / p%intensity - p%orbit(n)**2)) 
        p%init_orbit(n) = sum(pixel%pt(:,j)%init_orbit(n) * pixel%pt(:,j)%intensity) / p%intensity
        p%init_orbit_rms(n) = sqrt(max(0.0_rp, sum(pixel%pt(:,j)%init_orbit_rms(n) * pixel%pt(:,j)%intensity) / p%intensity - p%init_orbit(n)**2)) 
      enddo
    endif

    write (5, '(i6, es11.3, 3es11.3, i12, 2f10.3, 12es11.3)') j, j*pixel%dr(2)+pixel%r0(2), &
                       p%intensity_x * normalization, p%intensity_y * normalization, p%intensity * normalization, &
                       p%n_photon, p%orbit(6), p%orbit_rms(6), p%orbit(2), p%orbit_rms(2), p%orbit(4), p%orbit_rms(4), &
                       p%init_orbit(1), p%init_orbit_rms(1), p%init_orbit(2), p%init_orbit_rms(2), p%init_orbit(3), p%init_orbit_rms(3), p%init_orbit(4), p%init_orbit_rms(4)
  enddo
  close(5)

  !------------------------------------------
  ! det_pix_out_file.histogram

  if (lux_param%histogram_variable /= '') then
    nrm2 = normalization / lux_param%histogram_bin_width
    open (6, file = trim(lux_param%det_pix_out_file) // '.histogram', recl = 240)
    write (6, '(a)')  '#-----------------------------------------------------'
    write (6, '(a, t130, a)')  '#', '|                    Init'
    write (6, '(3a)') '# indx', adjustr(lux_param%histogram_variable(1:14)), '   Intens_x   Intens_y  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave  Ang_y_rms |    X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms  Ang_y_ave  Ang_y_rms'

    do j = ubound(lux_com%histogram_bin,1), lbound(lux_com%histogram_bin,1), -1
      pix => lux_com%histogram_bin(j)
      if (pix%intensity == 0) cycle

      pix%orbit          = pix%orbit / pix%intensity
      pix%orbit_rms      = pix%orbit_rms / pix%intensity
      pix%init_orbit     = pix%init_orbit / pix%intensity
      pix%init_orbit_rms = pix%init_orbit_rms / pix%intensity

      p = pix
      do n = 1, 6
        p%orbit_rms(n) = sqrt(max(0.0_rp, p%orbit_rms(n) - p%orbit(n)**2))
        p%init_orbit_rms(n) = sqrt(max(0.0_rp, p%init_orbit_rms(n) - p%init_orbit(n)**2))
      enddo

      write (6, '(i6, es14.3, 3es11.3, i12, 2f10.3, 12es11.3)') j, j * lux_param%histogram_bin_width, p%intensity_x * nrm2, &
                         p%intensity_y * nrm2, p%intensity * nrm2, p%n_photon, p%orbit(6), p%orbit_rms(6), &
                         p%orbit(2), p%orbit_rms(2), p%orbit(4), p%orbit_rms(4), &
                         p%init_orbit(1), p%init_orbit_rms(1), p%init_orbit(2), p%init_orbit_rms(2), p%init_orbit(3), p%init_orbit_rms(3), p%init_orbit(4), p%init_orbit_rms(4)
    enddo

    close(6)
  endif
endif

!------------------------------------------
! Track out file

if (lux_param%track_out_file /= '') then
  open (7, file = lux_param%track_out_file)
  write (7, '(3x, a, 2x, a, 16x, a, 20x, a, 4x, a)') 'Ix', 'Name', 'Class', 'S', '#Lost'
  do n = 0, lux_com%tracking_branch%n_ele_track
    ele => lux_com%tracking_branch%ele(n)
    write (7, '(i4, 2x, a20, a16, f10.3, i8)') n, ele%name, key_name(ele%key), ele%s, ele%ix_pointer
  enddo
  close(7)
endif

!------------------------------------------

call run_timer ('READ', dtime)

if (lux_com%verbose) then
  if (pixel%n_hit_pixel /= sum(pixel%pt%n_photon)) then
    print '(a, i0)', 'PARITY CHECK FAILED! PLEASE CONTACT DAVID SAGAN! ', pixel%n_hit_pixel - sum(pixel%pt%n_photon)
  endif
  print '(a, t35, i0, 5x, es10.2)', 'Photons Tracked:', pixel%n_track_tot, float(pixel%n_track_tot)
  print '(a, t35, i0)',      'Photons at detector:', pixel%n_hit_detec
  print '(a, t35, i0)',      'Photons hitting pix:', pixel%n_hit_pixel
  print '(a, t45, es12.4)', 'Normalization factor:', normalization
  print '(a, t45, es12.4)', 'Total intensity on pix (unnormalized):', intens_tot
  print '(a, t45, es12.4)', 'Total intensity on pix (normalized):', intens_tot * normalization

  if (intens_tot /= 0) then
    x_sum = 0; x2_sum = 0; y_sum = 0; y2_sum = 0
    e_ave = 0; e_rms = 0
    do nx = lbound(pixel%pt, 1), ubound(pixel%pt, 1)
    do ny = lbound(pixel%pt, 2), ubound(pixel%pt, 2)
      pix => pixel%pt(nx,ny)
      x = nx*pixel%dr(1) + pixel%r0(1)
      y = ny*pixel%dr(2) + pixel%r0(2)
      x_sum  = x_sum  + pix%intensity * x
      x2_sum = x2_sum + pix%intensity * x**2
      y_sum  = y_sum  + pix%intensity * y
      y2_sum = y2_sum + pix%intensity * y**2
      e_ave = e_ave   + pix%intensity * pix%orbit(6)
      e_rms = e_rms   + pix%intensity * pix%orbit_rms(6)
    enddo
    enddo

    x_ave = x_sum / intens_tot; y_ave = y_sum / intens_tot
    e_ave = e_ave / intens_tot; e_rms = e_rms / intens_tot
    print '(a, f10.3)',  '%Intensity not in det_pix file (due to %intensity_min_det_pixel_cutoff):', &
                                                               100 * (intens_tot - pix_in_file_intens) / intens_tot
    print '(a, 2f12.5)', 'Average position at det (x, y) (mm):     ', 1000 * x_ave, 1000 * y_ave
    print '(a, 2f12.5)', 'RMS at det (x, y) (mm):                  ', 1000 * sqrt(x2_sum / intens_tot - x_ave**2), &
                                                                      1000 * sqrt(y2_sum / intens_tot - y_ave**2)
    print '(a, 2f12.5)', 'Average energy deviation at det (eV):    ', e_ave
    print '(a, 2f12.5)', 'RMS energy at det (eV):                  ', sqrt(max(0.0_rp, e_rms - e_ave**2))
  endif

  print '(a, f10.2)', &
          ' Simulation time (min):               ', dtime/60
  print *
  print *, 'Tracking data file:       ', trim(lux_param%track_out_file)
  print *, 'Detector pixel data file: ', trim(lux_param%det_pix_out_file)
  print *, 'Photon1 data file:        ', trim(lux_param%photon1_out_file)
endif

! End plotting

if (lux_param%plotting /= '') then
  if (index('detector', trim(lux_param%plotting)) == 1) then
    call lux_plot_detector (lux_param, lux_com)
  endif
endif

end subroutine lux_write_data

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine lux_write_header (who, fmt, re, int, logic)
!
! Routine to write header info to data files.
!-

subroutine lux_write_header (who, fmt, re, int, int8, logic)

character(*) who
character(*), optional :: fmt
character(200) str

real(rp), optional :: re
integer, optional :: int
integer(8), optional :: int8
logical, optional :: logic

!

if (present(re)) then
  write (str, '('//fmt//')') re
  str = who // str
elseif (present(int)) then
  write (str, '('//fmt//')') int
  str = who // str
elseif (present(int8)) then
  write (str, '('//fmt//')') int8
  str = who // str
elseif (present(logic)) then
  write (str, '('//fmt//')') logic
  str = who // str
else
  str = who
endif

write (3, '(a)') trim('# ' // str)
write (4, '(a)') trim('# ' // str)
write (5, '(a)') trim('# ' // str)

end subroutine lux_write_header

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine lux_plot_detector (lux_param, lux_com)
!
! Routine to plot the intensity at the detector pixels
!-

subroutine lux_plot_detector (lux_param, lux_com)

type (lat_struct), pointer:: lat
type (ele_struct), pointer :: ele
type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param

real(rp) x1, x2, y1, y2
integer ix, iy, i_chan
character(16) ans

!

call qp_open_page ('X', i_chan, lux_param%window_width, lux_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.10_rp, 0.05_rp, '%PAGE')

lat => lux_com%lat
ele => lat%ele(lat%n_ele_track)

call qp_calc_and_set_axis ('X', -ele%value(x1_limit$), ele%value(x2_limit$), 6, 10, 'GENERAL')
call qp_calc_and_set_axis ('Y', -ele%value(y1_limit$), ele%value(y2_limit$), 6, 10, 'GENERAL')

call qp_draw_axes ('X', 'Y', 'Detector Pixels')

do ix = 1, lbound(lux_com%pixel%pt, 1), ubound(lux_com%pixel%pt, 1)
do iy = 1, lbound(lux_com%pixel%pt, 2), ubound(lux_com%pixel%pt, 2)
  call qp_paint_rectangle (x1, x2, y1, y2, color = 'red')
enddo
enddo

call read_a_line ('CR to exit: ', ans)

end subroutine lux_plot_detector

end module
