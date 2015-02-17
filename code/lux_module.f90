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

!

type lux_bend_slice_struct
  type (ele_struct) ele
  type (em_field_struct) field
  type (coord_struct) orbit
  real(rp) integrated_emit_prob    ! Probability of photon emission from element start to slice end
  real(rp) emit_prob               ! Probability of photon emission withing slice
  real(rp) g_x, g_y
  logical good_emit                ! Only generate photons in slices with %good_emit = True. Photons generated 
                                   !   in slices with good_emit = False will not make it through the aperture
end type

! %bend_slice(i) are the emission parameters at the end of the i^th slice

type lux_param_struct
  character(100) :: photon1_out_file = ''
  character(100) :: det_pix_out_file = ''
  character(100) :: param_file = 'lux.init'
  character(100) :: lattice_file = ''
  character(40) :: source_element = ''            ! element name
  character(40) :: detector_element = ''          ! element name
  character(40) :: photon1_element = ''           ! element name
  character(20) :: plotting = ''
  character(16) :: random_engine = 'pseudo'
  real(rp) :: intensity_min_det_pixel_cutoff = 1e-6
  real(rp) :: intensity_min_photon1_cutoff = 1e-6
  real(rp) :: stop_total_intensity = 10           ! stop intensity per energy
  real(rp) :: window_width = 800.0_rp, window_height = 400.0_rp  ! For plotting
  real(rp) :: intensity_normalization_coef = 1e6
  real(rp) :: stop_num_photons = 0                  ! stop number. Use real for clearer input
  real(rp) :: mpi_run_size = 0.1                 ! Normalized number of photons to track
  integer :: n_energy_bin_pts = 40
  integer :: random_seed = 0
  logical :: debug = .false.     ! For debugging
  logical :: reject_dead_at_det_photon1 = .false.
  logical :: scale_initial_field_to_1 = .true.
end type

type lux_photon_struct
  integer n_photon_generated
  type (coord_struct), allocatable :: orb(:)
end type

type lux_common_struct
  type (lat_struct) :: lat
  type (branch_struct), pointer :: source_branch, detec_branch
  type (ele_struct), pointer :: source_ele, fork_ele, detec_ele, photon1_ele
  integer n_bend_slice                                       ! Number of slices
  integer :: mpi_rank  = -1
  integer :: mpi_n_proc = 1                      ! Number of processeses including master
  type (lux_bend_slice_struct), allocatable :: bend_slice(:) ! Size: (0:n_bend_slice)
  type (surface_grid_pt_struct), allocatable :: energy_bin(:)
  real(rp) dE_bin
  real(rp) E_min, E_max                                      ! Photon energy range 
  logical :: verbose = .false.
  logical :: using_mpi = .false.
end type

type lux_output_data_struct
  integer :: n_track_tot = 0
  integer :: nx_min = 0, nx_max = 0, ny_min = 0, ny_max = 0
  integer :: n_live = 0
  integer :: n_lost = 0
end type

contains

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
type (ele_struct), pointer :: source_ele

real(rp) cut, r

integer i, j, n, ix, ie, ir
integer n_throw, n_loc

character(3) num_str
character(40) arg, plotting, number_file

logical ok, is_there, err


namelist / params / lux_param

! Get inputs

lux_com%verbose = .true.

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  select case (arg)
  case ('-silent')
    lux_com%verbose = .false.
  case ('-plot')
    i = i + 1
    call cesr_getarg (i, plotting)
  case default
    if (arg(1:1) == '-') then
      print *, 'I DO NOT UNDERSTAND: ', trim(arg)
      print *
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
  print *, 'Note: lux_param%debug = True'
endif

if (.not. lux_com%verbose) call output_direct (0, .false., max_level = s_success$) ! Suppress bmad_parser info output.

call ran_seed_put (lux_param%random_seed)
if (lux_com%using_mpi) then
  call ran_seed_get (ir)
  call ran_sed_put (ir + 100 * lux_com%mpi_rank)
endif

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

if (index(lux_param%photon1_out_file, '#') /= 0 .or. index(lux_param%det_pix_out_file, '#') /= 0) then

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

endif

! Read in lattice

call bmad_parser (lux_param%lattice_file, lux_com%lat)
lat => lux_com%lat

! locate source element

call lat_ele_locator (lux_param%source_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  print *, 'NO SOURCE ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
elseif (n_loc > 1) then
  print *, 'MULTIPLE SOURCE ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
endif

lux_com%source_ele => eles(1)%ele
lux_com%source_branch => lux_com%source_ele%branch

source_ele => lux_com%source_ele

! Locate detector element

call lat_ele_locator (lux_param%detector_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  print *, 'NO DETECTOR ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
elseif (n_loc > 1) then
  print *, 'MULTIPLE DETECTOR ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
endif

lux_com%detec_ele => eles(1)%ele
lux_com%detec_branch => lux_com%detec_ele%branch

! Locate photon1 element

if (lux_param%photon1_element == '') then
  lux_com%photon1_ele => lux_com%detec_ele
else
  call lat_ele_locator (lux_param%photon1_element, lat, eles, n_loc, err)
  if (n_loc == 0) then
    print *, 'NO PHOTON1 ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  elseif (n_loc > 1) then
    print *, 'MULTIPLE PHOTON1 ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  endif
  lux_com%photon1_ele => eles(1)%ele 
endif

!

select case (source_ele%key)
case (x_ray_source$)

case (sbend$, wiggler$, undulator$)
  nullify (lux_com%fork_ele)
  do i = 1, lux_com%source_branch%n_ele_track
    if (lux_com%source_branch%ele(i)%key /= photon_fork$) cycle
    lux_com%fork_ele => lux_com%source_branch%ele(i)
    exit
  enddo
  if (.not. associated(lux_com%fork_ele)) then
    print *, 'NO APPROPRIATE FORK/PHOTON_FORK ELEMENT FOUND IN LATTICE!'
    stop
  endif

case default
  print *, 'CANNOT SIMULATE PHOTONS GENERATED IN ELEMENT OF TYPE: ', trim(key_name(source_ele%key))
  stop
end select

! Some init

call run_timer('START')

!

call lux_tracking_setup (lux_param, lux_com)

if (lat%photon_type == coherent$) then
  if (source_ele%value(e_field_x$) == 0 .and. source_ele%value(e_field_y$) == 0) then
    print *, 'WARNING: INPUT E_FIELD IS ZERO SO RANDOM FILED WILL BE GENERATED WITH COHERENT PHOTONS!'
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

end subroutine

end subroutine lux_init

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_init_data (lux_param, lux_com, lux_data)
!
! Routine to initialize the output data.
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!
! Ouput:
!   lux_data    -- lux_output_data_struct: Initialized data.
!-

subroutine lux_init_data (lux_param, lux_com, lux_data)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lux_output_data_struct) lux_data
type (surface_grid_struct), pointer :: detec_grid


!

detec_grid => lux_com%detec_ele%photon%surface%grid

if (.not. allocated(detec_grid%pt)) then
  print *, 'DETECTOR GRID NOT SET!'
  stop
endif

lux_data%nx_min = lbound(detec_grid%pt, 1)
lux_data%nx_max = ubound(detec_grid%pt, 1)
lux_data%ny_min = lbound(detec_grid%pt, 2)
lux_data%ny_max = ubound(detec_grid%pt, 2)

lux_data%n_track_tot = 0



end subroutine lux_init_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_generate_photon (photon, lux_param, lux_com)
!
! Routine to generate the starting photon coordinates
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!
! Ouput:
!   photon     -- lux_photon_struct: Initialized starting coords.
!-

subroutine lux_generate_photon (photon, lux_param, lux_com)

use nr

type (lat_struct), pointer :: lat
type (lux_photon_struct), target :: photon
type (coord_struct) charged_orb
type (coord_struct), pointer :: orb
type (lux_param_struct) lux_param
type (lux_common_struct), target :: lux_com
type (lux_bend_slice_struct), pointer :: sl(:)
type (ele_struct) ele
type (ele_struct), pointer :: source_ele

real(rp) x, y, phi, r(3), dir(2), ds, rr, r_emit(5), prob, f
real(rp) v_mat(4,4), v_inv_mat(4,4), vec(4), dE, e(3), b(3), sig_vec(6)
real(rp) g_bend(3), gamma_electron

integer ix, n_slice

! Init

orb => photon%orb(0)
source_ele => lux_com%source_ele
lat => lux_com%lat

! Field

x = source_ele%value(e_field_x$); y = source_ele%value(e_field_y$)
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

!-----------------------------------------------------
! x_ray_source source

select case (source_ele%key)
case (x_ray_source$)

  ! Set position

  if (nint(source_ele%value(spatial_distribution$)) == uniform$) then
    call ran_uniform(r)
    r = 2 * r - 1
  elseif (nint(source_ele%value(spatial_distribution$)) == gaussian$) then
    call ran_gauss(r)
  else
    print *, 'BAD SOURCE_ELE SPATIAL_DISTRIBUTION: ', distribution_name(nint(source_ele%value(spatial_distribution$)))
    call err_exit
  endif

  orb%vec(1:5:2) = lat%beam_start%vec(1:5:2) + &
                    [source_ele%value(sig_x$) * r(1), source_ele%value(sig_y$) * r(2), source_ele%value(sig_z$) * r(3)]

  ! Set direction

  select case (nint(source_ele%value(velocity_distribution$)))
  case (spherical$)
    call isotropic_photon_emission (source_ele, lat%param, orb, +1, twopi)

  case (uniform$)
    call ran_uniform(dir)
    dir = 2 * dir - 1
    orb%vec(2:4:2) = lat%beam_start%vec(2:4:2) + dir * [source_ele%value(sig_vx$), source_ele%value(sig_vy$)]
    orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  case (gaussian$)
    call ran_gauss(dir)
    orb%vec(2:4:2) = lat%beam_start%vec(2:4:2) + dir * [source_ele%value(sig_vx$), source_ele%value(sig_vy$)]
    orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  end select

  ! Set energy

  if (nint(source_ele%value(energy_distribution$)) == uniform$) then
    call ran_uniform(rr)
    rr = 2 * rr - 1
  else if (nint(source_ele%value(energy_distribution$)) == gaussian$) then
    call ran_gauss(rr)
  else
    print *, 'BAD SOURCE_ELE%VALUE(ENERGY_DISTRIBUTION SETTING: ', distribution_name(nint(source_ele%value(energy_distribution$)))
  endif

  orb%p0c = source_ele%value(sig_E$) * rr + source_ele%value(dE_center$)
  if (is_true(source_ele%value(dE_relative_to_ref$))) orb%p0c = orb%p0c + source_ele%value(p0c$) 

  call init_coord (orb, orb%vec, source_ele, upstream_end$, photon$, 1, orb%p0c) 
  orb%s = orb%vec(5) + orb%s + source_ele%value(z_offset_tot$)
  orb%t = 0

  ! Translate from element to lab coordinates
  ! and track to entrance end of source_ele

  call offset_photon (source_ele, orb, unset$)

  call track_a_drift_photon (orb, -orb%s, .true.)

  return

!-----------------------------------------------------
! bend, wiggler, undulator source

case default
  sl => lux_com%bend_slice
  n_slice = ubound(sl, 1)

  ! Find where photon emitted

  call ran_uniform(rr)  ! longitudinal position
  call bracket_index (sl%integrated_emit_prob, 0, n_slice, rr, ix)
  ix = ix + 1
  if (ix == n_slice) ix = n_slice - 1
  f = (rr - sl(ix-1)%integrated_emit_prob) / (sl(ix)%integrated_emit_prob - sl(ix-1)%integrated_emit_prob)

  ! Calculate electron average position

  charged_orb = sl(ix)%orbit
  charged_orb%vec = (1-f) * sl(ix)%orbit%vec + f * sl(ix+1)%orbit%vec

  ele%a = average_twiss(1-f, sl(ix)%ele%a, sl(ix+1)%ele%a)
  ele%b = average_twiss(1-f, sl(ix)%ele%b, sl(ix+1)%ele%b)
  ele%c_mat   = (1-f) * sl(ix)%ele%c_mat   + f * sl(ix+1)%ele%c_mat
  ele%gamma_c = (1-f) * sl(ix)%ele%gamma_c + f * sl(ix+1)%ele%gamma_c
  call make_v_mats (ele, v_mat, v_inv_mat)

  ! Add offsets due to finite bunch size to the electron position.
  ! To do this must transform to the normal mode coords

  call ran_gauss(r_emit)  ! electron momentum offset.
  dE = r_emit(5) * source_ele%value(sig_E$) / source_ele%value(p0c$)
  charged_orb%vec(6) = charged_orb%vec(6) + dE

  vec = matmul (v_inv_mat, charged_orb%vec(1:4))
  vec(1:2) = vec(1:2) + charged_orb%vec(1:2) + [ele%a%eta, ele%a%etap] * dE
  vec(3:4) = vec(3:4) + charged_orb%vec(1:2) + [ele%b%eta, ele%b%etap] * dE

  vec(1) = vec(1) + sqrt(lux_com%source_branch%a%emit * ele%a%beta) * r_emit(1)
  vec(2) = vec(2) + sqrt(lux_com%source_branch%a%emit / ele%a%beta) * (r_emit(2) - ele%a%alpha * r_emit(1))

  vec(3) = vec(3) + sqrt(lux_com%source_branch%b%emit * ele%b%beta) * r_emit(3)
  vec(4) = vec(4) + sqrt(lux_com%source_branch%b%emit / ele%b%beta) * (r_emit(4) - ele%b%alpha * r_emit(3))

  charged_orb%vec(1:4) = matmul(v_mat, vec)

  ! Calculate bending strength

  B = (1-f) * sl(ix)%field%b + f * sl(ix+1)%field%b
  E = 0
  g_bend = g_bend_from_em_field (B, E, charged_orb)
  
  ! Init photon

  gamma_electron = source_ele%value(p0c$) * (1 + sl(ix)%orbit%vec(6)) / sl(ix)%orbit%beta / mass_of(sl(ix)%orbit%species)
  !! Note: Energy slices not yet implemented.
  !! if (lux_com%lat%photon_type == coherent$ .and. ix_energy > 0) then
  !!   rr = (ix_energy - 0.5_rp) / lux_param%n_energy_pts
  !!   call bend_photon_init (g_bend(1), g_bend(2), gamma_electron, orb, lux_com%E_min, lux_com%E_max, rr)
  !! endif

  call bend_photon_init (g_bend(1), g_bend(2), gamma_electron, orb, lux_com%E_min, lux_com%E_max)
  call absolute_photon_position (charged_orb, orb)

  orb%s = sl(ix-1)%ele%s + f * sl(ix)%ele%value(l$)

  ! Track to fork element.

  ds = lux_com%fork_ele%s - orb%s  

  if (source_ele%key == sbend$) then
    call track_a_bend_photon (orb, ele, ds)
  else
    call track_a_drift_photon (orb, ds, .true.)
  endif

  return
end select

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
type (ele_struct), pointer :: ele, fork_ele, source_ele
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

!-------------------------------------------------------------
! Lattice has sample element.

lat => lux_com%lat
branch => lux_com%detec_branch
source_ele => lux_com%source_ele

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  select case (ele%key)
  case (sample$, diffraction_plate$) 
    call photon_target_setup (ele)
  end select
enddo

!-------------------------------------------------------------
! x_ray_source source

select case (source_ele%key)
case (x_ray_source$)
  call photon_target_setup (source_ele)

  if (nint(source_ele%value(energy_distribution$)) == uniform$) then
    f = 1
  else
    f = 3
  endif
  lux_com%E_min = source_ele%value(dE_center$) - f * source_ele%value(sig_E$)
  lux_com%E_max = source_ele%value(dE_center$) + f * source_ele%value(sig_E$)
  if (is_false(source_ele%value(dE_relative_to_ref$))) then
    lux_com%E_min = lux_com%E_min - source_ele%value(p0c$) 
    lux_com%E_max = lux_com%E_max - source_ele%value(p0c$) 
  endif

!-------------------------------------------------------------
! Sbend or wiggler source

case default
  call photon_target_setup (lux_com%fork_ele)
  
  n_slice = max(1, nint(source_ele%value(l$) / source_ele%value(ds_slice$)))
  lux_com%n_bend_slice = n_slice
  allocate (lux_com%bend_slice(0:n_slice))

  lux_com%E_min = source_ele%value(dE_center$) - source_ele%value(sig_E$)
  if (is_true(source_ele%value(dE_relative_to_ref$))) lux_com%E_min = lux_com%E_min + lux_com%detec_ele%value(p0c$) 
  if (source_ele%value(sig_E$) == 0) then
    lux_com%E_max = lux_com%E_min + 1d-10  ! Need some small offset for the calculation
  else
    lux_com%E_max = lux_com%E_min + 2 * source_ele%value(sig_E$)
  endif

  ! Track through source ele and gather data

  call init_coord (orb, lat%beam_start, source_ele, upstream_end$)
  twiss_ele = pointer_to_next_ele (source_ele, -1)
  ds = source_ele%value(l$) / n_slice
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
    call em_field_calc (source_ele, lat%param, s_now, 0.0_rp, orb, .false., sl(i)%field)

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
                              sqrt(lux_com%source_branch%a%emit * twiss_ele%a%beta) * source_ele%value(transverse_sigma_cut$)
    p_coord(2)%vec = matmul(v_mat, [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]) * &
                              sqrt(lux_com%source_branch%a%emit / twiss_ele%a%beta) * source_ele%value(transverse_sigma_cut$)
    p_coord(3)%vec = matmul(v_mat, [0.0_rp, 0.0_rp, 1.0_rp, -twiss_ele%b%alpha]) * &
                              sqrt(lux_com%source_branch%b%emit * twiss_ele%b%beta) * source_ele%value(transverse_sigma_cut$)
    p_coord(4)%vec = matmul(v_mat, [0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) * &
                              sqrt(lux_com%source_branch%b%emit / twiss_ele%b%beta) * source_ele%value(transverse_sigma_cut$)

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

    call twiss_and_track_intra_ele (source_ele, lat%param, s_now, s_now+ds, &
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
  sl%integrated_emit_prob = sl%integrated_emit_prob / sl(n_slice)%integrated_emit_prob

  if (lux_com%verbose) print '(a, i4)', &
            'Number of slices of source element to be used for photon generation:', count(sl%good_emit) 

end select

! Setup energy binning

if (lux_param%n_energy_bin_pts == 0) lux_param%n_energy_bin_pts = 1
allocate (lux_com%energy_bin(lux_param%n_energy_bin_pts))

lux_com%dE_bin = (lux_com%E_max - lux_com%E_min) / lux_param%n_energy_bin_pts
if (lux_com%dE_bin == 0) lux_com%de_bin = 1e-5

do i = 1, lux_param%n_energy_bin_pts
  lux_com%energy_bin(i)%energy_ave = lux_com%E_min + lux_com%dE_bin * (i - 0.5)
enddo
if (is_true(source_ele%value(dE_relative_to_ref$))) lux_com%energy_bin%energy_ave = &
                                                  lux_com%energy_bin%energy_ave - lux_com%detec_ele%value(p0c$) 

end subroutine lux_tracking_setup 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_tracking (lux_param, lux_com, lux_data)
!
! Routine to init Lux.
!
! Input:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!
! Output:
!   lux_data    -- lux_output_data_struct: 
!-

subroutine lux_tracking (lux_param, lux_com, lux_data)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lux_output_data_struct) lux_data
type (lux_photon_struct), target :: photon
type (lat_struct), pointer :: lat
type (coord_struct), pointer :: end_orb, this_orb
type (ele_struct), pointer :: detec_ele, source_ele
type (branch_struct), pointer :: s_branch, d_branch
type (surface_grid_struct), pointer :: detec_grid
type (surface_grid_pt_struct), pointer :: pix
type (surface_grid_pt_struct) :: pixel

real(rp) intens, intens_x, intens_y
real(rp) e_ref, phase

integer ix, nt, nx, ny, track_state, stop_num_photons

logical accept

!

lat => lux_com%lat
detec_ele => lux_com%detec_ele
detec_grid => lux_com%detec_ele%photon%surface%grid
source_ele => lux_com%source_ele
d_branch => lux_com%detec_branch
s_branch => lux_com%source_branch

nt = d_branch%n_ele_track
stop_number_of_photons = lux_param%stop_number_of_photons
if (lux_com%using_mpi) stop_number_of_photons = 1 + &
                              stop_number_of_photons * lux_param%mpi_run_size / (lux_com%mpi_n_proc - 1)

call reallocate_coord (photon%orb, lat, d_branch%ix_branch)

if (lux_param%photon1_out_file /= '') then
  open (1, file = lux_param%photon1_out_file, recl = 240)
  write (1, '(a)') '#      |                               Start (x 1e3)                    |                                    End (x 1e3)                     |                 End'
  write (1, '(a)') '#   Ix |          x           vx            y           vy            z |             x           vx            y           vy            z  |     Energy     Intens_x     Intens_y'
endif

!

do 
  if (stop_num_photons > 0 .and. lux_data%n_track_tot >= stop_num_photons) exit
  lux_data%n_track_tot = lux_data%n_track_tot + 1
  photon%n_photon_generated = lux_data%n_track_tot

  call lux_generate_photon (photon, lux_param, lux_com)
  if (lux_param%debug) then
    call init_coord (photon%orb(0), lat%beam_start, d_branch%ele(0), downstream_end$, photon$, &
                                        1, d_branch%ele(0)%value(E_tot$) * (1 + lat%beam_start%vec(6)))
  endif

  call track_all (lat, photon%orb, d_branch%ix_branch, track_state)

  end_orb => photon%orb(nt)
  call offset_photon (detec_ele, end_orb, set$)

  intens_x = end_orb%field(1)**2
  intens_y = end_orb%field(2)**2
  intens = intens_x + intens_y

  ! Write results

  if (lux_param%photon1_out_file /= '' .and. intens >= lux_param%intensity_min_photon1_cutoff) then
    this_orb => photon%orb(lux_com%photon1_ele%ix_ele)
    accept = .true.
    if (lux_param%reject_dead_at_det_photon1 .and. track_state /= moving_forward$) accept = .false.
    if (this_orb%state /= alive$) accept = .false.
    if (accept) then
      write (1, '(i6, 5f13.6, 3x, 5f13.6, 3x, f11.3, 2es13.4)') lux_data%n_track_tot, 1d3*photon%orb(1)%vec(1:5), &
                      1d3*this_orb%vec(1:5), this_orb%p0c, end_orb%field(1)**2, end_orb%field(2)**2
    endif
  endif

  !

  if (track_state /= moving_forward$) then
    lux_data%n_lost = lux_data%n_lost + 1
    cycle
  endif

  ! Go to coordinates of the detector

  lux_data%n_live = lux_data%n_live + 1

  nx = nint((end_orb%vec(1) - detec_grid%r0(1)) / detec_grid%dr(1))
  ny = nint((end_orb%vec(3) - detec_grid%r0(2)) / detec_grid%dr(2))

  if (lux_data%nx_min <= nx .and. nx <= lux_data%nx_max .and. lux_data%ny_min <= ny .and. ny <= lux_data%ny_max) then
    pix => detec_grid%pt(nx,ny)
    pix%n_photon  = pix%n_photon + 1
    if (lat%photon_type == coherent$) then
      phase = end_orb%phase(1) 
      pix%E_x = pix%E_x + end_orb%field(1) * [cos(phase), sin(phase)]
      phase = end_orb%phase(2) 
      pix%E_y = pix%E_y + end_orb%field(2) * [cos(phase), sin(phase)]
    else
      pix%intensity_x = pix%intensity_x + intens_x
      pix%intensity_y = pix%intensity_y + intens_y
      pix%intensity   = pix%intensity   + intens
      pix%energy_ave  = pix%energy_ave  + intens * (end_orb%p0c - e_ref)
      pix%energy_rms  = pix%energy_rms  + intens * (end_orb%p0c - e_ref)**2
    endif

    if (is_true(source_ele%value(dE_relative_to_ref$))) then
      ix = nint((end_orb%p0c - lux_com%detec_ele%value(p0c$) - lux_com%energy_bin(1)%energy_ave) / lux_com%dE_bin) + 1 
    else
      ix = nint((end_orb%p0c - lux_com%energy_bin(1)%energy_ave) / lux_com%dE_bin) + 1 
    endif
    if (ix < 1) ix = 1
    if (ix > ubound(lux_com%energy_bin, 1)) ix = ubound(lux_com%energy_bin, 1)
    lux_com%energy_bin(ix)%intensity   = intens
    lux_com%energy_bin(ix)%intensity_x = intens_x
    lux_com%energy_bin(ix)%intensity_y = intens_y
    lux_com%energy_bin(ix)%n_photon    = lux_com%energy_bin(ix)%n_photon + 1

  endif

enddo

!

close(1)

if (lat%photon_type == coherent$) then
  do nx = lux_data%nx_min, lux_data%nx_max; do ny = lux_data%ny_min, lux_data%ny_max
    pix => detec_grid%pt(nx,ny)
    intens_x = pix%E_x(1)**2 + pix%E_x(2)**2
    intens_y = pix%E_y(1)**2 + pix%E_y(2)**2
    intens = intens_x + intens_y
    pix%intensity_x = intens_x
    pix%intensity_y = intens_y
    pix%intensity   = pix%intensity + intens
  enddo; enddo
endif

  do nx = lux_data%nx_min, lux_data%nx_max; do ny = lux_data%ny_min, lux_data%ny_max
    pix => detec_grid%pt(nx,ny)
    if (pix%intensity == 0) cycle
    pix%energy_ave = pix%energy_ave / pix%intensity
    pix%energy_rms = pix%energy_rms / pix%intensity
  enddo; enddo



end subroutine lux_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_write_data (lux_param, lux_com, lux_data)
!
! Routine to write the photon tracking data .
!
! Output:
!   lux_param   -- lux_param_struct: Lux input parameters.
!   lux_com     -- lux_common_struct: Common parameters.
!   lux_data    -- lux_output_data_struct: Tracking data.
!-

subroutine lux_write_data (lux_param, lux_com, lux_data)

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lux_output_data_struct) lux_data
type (surface_grid_struct), pointer :: detec_grid
type (surface_grid_pt_struct), pointer :: pix
type (surface_grid_pt_struct) :: pixel
type (lat_struct), pointer :: lat

real(rp) normalization, area, cut, dtime
real(rp) total_dead_intens, pix_in_file_intens
real(rp) :: intens_tot, intens_tot_x, intens_tot_y, intens_max
real(rp) x_sum, y_sum, x2_sum, y2_sum, x_ave, y_ave, e_rms, e_ave, e_ref
real(rp) phase_x, phase_y, x, y

integer i, j, nx, ny
integer nx_active_min, nx_active_max, ny_active_min, ny_active_max

!------------------------------------------

detec_grid => lux_com%detec_ele%photon%surface%grid
lat => lux_com%lat

!------------------------------------------
! lux_param%det_pix_out_file

pix_in_file_intens = 0
intens_tot_x = 0
intens_tot_y = 0
intens_tot = 0

area = fourpi
normalization = lux_param%intensity_normalization_coef * area / (lux_data%n_track_tot * fourpi)

if (lux_param%det_pix_out_file /= '') then
  open (3, file = lux_param%det_pix_out_file, recl = 160)
  intens_max = maxval(detec_grid%pt%intensity)
  cut = intens_max * lux_param%intensity_min_det_pixel_cutoff

  nx_active_min = lux_data%nx_max;    nx_active_max = lux_data%nx_min
  ny_active_min = lux_data%ny_max;    ny_active_max = lux_data%ny_min

  do i = lux_data%nx_min, lux_data%nx_max
  do j = lux_data%ny_min, lux_data%ny_max
    if (detec_grid%pt(i,j)%intensity <= cut) cycle
    nx_active_min = min(nx_active_min, i);  nx_active_max = max(nx_active_max, i)
    ny_active_min = min(ny_active_min, j);  ny_active_max = max(ny_active_max, j)
  enddo
  enddo

  write (3, '(3a)')        'master_input_file   = "', trim(lux_param%param_file), '"'
  write (3, '(3a)')        'lattice_file        = "', trim(lux_param%lattice_file), '"'
  write (3, '(a, es14.6)') 'normalization       =', normalization
  write (3, '(a, es16.5)') 'intensity_x_unnorm  =', sum(detec_grid%pt(:,:)%intensity_x)
  write (3, '(a, es16.5)') 'intensity_x_norm    =', sum(detec_grid%pt(:,:)%intensity_x) * normalization
  write (3, '(a, es16.5)') 'intensity_y_unnorm  =', sum(detec_grid%pt(:,:)%intensity_y)
  write (3, '(a, es16.5)') 'intensity_y_norm    =', sum(detec_grid%pt(:,:)%intensity_y) * normalization
  write (3, '(a, es16.5)') 'intensity_unnorm    =', sum(detec_grid%pt(:,:)%intensity)
  write (3, '(a, es16.5)') 'intensity_norm      =', sum(detec_grid%pt(:,:)%intensity) * normalization
  write (3, '(a, f10.6)')  'dx_pixel            =', detec_grid%dr(1)
  write (3, '(a, f10.6)')  'dy_pixel            =', detec_grid%dr(2)
  write (3, '(a, i8)')     'nx_active_min       =', nx_active_min
  write (3, '(a, i8)')     'nx_active_max       =', nx_active_max
  write (3, '(a, i8)')     'ny_active_min       =', ny_active_min
  write (3, '(a, i8)')     'ny_active_max       =', ny_active_max
  write (3, '(a)')         '#-----------------------------------------------------'
  write (3, '(a)')         '#     ix      iy        x_pix        y_pix   Intensity_x       Phase_x   Intensity_y       Phase_y     Intensity  N_photn     E_ave     E_rms'

  do i = lux_data%nx_min, lux_data%nx_max
  do j = lux_data%ny_min, lux_data%ny_max
    pix => detec_grid%pt(i,j)
    intens_tot = intens_tot + pix%intensity 
    if (pix%intensity <= cut .or. pix%n_photon == 0) cycle
    pix%energy_ave = pix%energy_ave / pix%intensity
    pix%energy_rms = sqrt(max(0.0_rp, pix%energy_rms / pix%intensity - pix%energy_ave**2)) 
    phase_x = 0; phase_y = 0
    if (lat%photon_type == coherent$) then
      phase_x = atan2(pix%e_x(1), pix%e_x(2))
      phase_y = atan2(pix%e_y(1), pix%e_y(2))
    endif
    write (3, '(2i8, 2es13.5, 5es14.5, i8, 2f10.3)') i, j, [i,j]*detec_grid%dr+detec_grid%r0, &
           pix%intensity_x * normalization, phase_x, pix%intensity_y * normalization, phase_y, &
           pix%intensity * normalization, pix%n_photon, pix%energy_ave, pix%energy_rms
    pix_in_file_intens = pix_in_file_intens + pix%intensity 
  enddo
  enddo
  close(3)

  ! det_pix_out_file.x

  open (3, file = trim(lux_param%det_pix_out_file) // '.x')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     ix        x_pix      ntensity_x     Intensity_y       Intensity  N_photn     E_ave     E_rms'
  do i = lux_data%nx_min, lux_data%nx_max
    pixel = surface_grid_pt_struct()
    pixel%intensity_x = sum(detec_grid%pt(i,:)%intensity_x)
    pixel%intensity_y = sum(detec_grid%pt(i,:)%intensity_y)
    pixel%intensity   = sum(detec_grid%pt(i,:)%intensity)
    pixel%n_photon    = sum(detec_grid%pt(i,:)%n_photon)
    if (pixel%intensity == 0) then
      pixel%energy_ave = 0
      pixel%energy_rms = 0
    else
      pixel%energy_ave = sum(detec_grid%pt(i,:)%energy_ave) / pixel%intensity
      pixel%energy_rms = sqrt(max(0.0_rp, sum(detec_grid%pt(i,:)%energy_rms) / pixel%intensity - pixel%energy_ave**2))
    endif
    write (3, '(i8, es13.5, 3es16.5, i9, 2f10.3)') i, i*detec_grid%dr(1)+detec_grid%r0(1), &
                       pixel%intensity_x * normalization, pixel%intensity_y * normalization, pixel%intensity * normalization, &
                       pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

  ! det_pix_out_file.y

  open (3, file = trim(lux_param%det_pix_out_file) // '.y')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     iy        y_pix      Intensity_x     Intensity_y       Intensity  N_photn     E_ave     E_rms'
  do j = lux_data%ny_min, lux_data%ny_max
    pixel = surface_grid_pt_struct()
    pixel%intensity_x = sum(detec_grid%pt(:,j)%intensity_x)
    pixel%intensity_y = sum(detec_grid%pt(:,j)%intensity_y)
    pixel%intensity   = sum(detec_grid%pt(:,j)%intensity)
    pixel%n_photon    = sum(detec_grid%pt(:,j)%n_photon)
    if (pixel%intensity == 0) then
      pixel%energy_ave = 0
      pixel%energy_rms = 0
    else
      pixel%energy_ave = sum(detec_grid%pt(:,j)%energy_ave) / pixel%intensity
      pixel%energy_rms = sqrt(max(0.0_rp, sum(detec_grid%pt(:,j)%energy_rms) / pixel%intensity - pixel%energy_ave**2)) 
    endif
    write (3, '(i8, es13.5, 3es16.5, i8, 2f10.3)') j, j*detec_grid%dr(2)+detec_grid%r0(2), &
                       pixel%intensity_x * normalization, pixel%intensity_y * normalization, pixel%intensity * normalization, &
                       pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

  ! det_pix_out_file.energy

  open (3, file = trim(lux_param%det_pix_out_file) // '.energy')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#      Energy     Intensity_x     Intensity_y       Intensity  N_photn'
  do j = 1, ubound(lux_com%energy_bin, 1)
    pix => lux_com%energy_bin(j)
    write (3, '(f13.5, 3es16.5, i9)') pix%energy_ave, pix%intensity_x * normalization, &
                       pix%intensity_y * normalization, pix%intensity * normalization, pix%n_photon
  enddo
  close(3)

endif

!------------------------------------------

call run_timer ('READ', dtime)

if (lux_com%verbose) then
  print *, 'Photons Tracked:                       ', lux_data%n_track_tot
  print *, 'Photons at detector:                   ', lux_data%n_live
  print *, 'Normalization factor:                  ', normalization
  print *, 'Total intensity (unnormalized):        ', intens_tot
  print *, 'Total intensity (normalized):          ', intens_tot * normalization

  if (intens_tot /= 0) then
    x_sum = 0; x2_sum = 0; y_sum = 0; y2_sum = 0
    e_ave = 0; e_rms = 0
    do nx = lux_data%nx_min, lux_data%nx_max
    do ny = lux_data%ny_min, lux_data%ny_max
      pix => detec_grid%pt(nx,ny)
      x = nx*detec_grid%dr(1) + detec_grid%r0(1)
      y = ny*detec_grid%dr(2) + detec_grid%r0(2)
      x_sum  = x_sum  + pix%intensity * x
      x2_sum = x2_sum + pix%intensity * x**2
      y_sum  = y_sum  + pix%intensity * y
      y2_sum = y2_sum + pix%intensity * y**2
      e_ave = e_ave   + pix%intensity * pix%energy_ave
      e_rms = e_rms   + pix%intensity * pix%energy_rms
    enddo
    enddo

    x_ave = x_sum / intens_tot; y_ave = y_sum / intens_tot
    e_ave = e_ave / intens_tot; e_rms = e_rms / intens_tot
    print '(a, f10.3)', &
            ' %Intensity at pixels not in det_pix file:', 100 * (intens_tot - pix_in_file_intens) / intens_tot
    print '(a, 2f12.5)', 'Average position at det (x, y) (mm):     ', 1000 * x_ave, 1000 * y_ave
    print '(a, 2f12.5)', 'RMS at det (x, y) (mm):                  ', 1000 * sqrt(x2_sum / intens_tot - x_ave**2), &
                                                                      1000 * sqrt(y2_sum / intens_tot - y_ave**2)
    print '(a, 2f12.5)', 'Average energy deviation at det (eV):    ', e_ave
    print '(a, 2f12.5)', 'RMS energy at det (eV):                  ', sqrt(max(0.0_rp, e_rms - e_ave**2))
  endif

  print '(a, f10.2)', &
          ' Simulation time (min):               ', dtime/60
  print *
  print *, 'Photon1 data file:        ', trim(lux_param%photon1_out_file)
  print *, 'Detector pixel data file: ', trim(lux_param%det_pix_out_file)
endif

! End plotting

if (lux_param%plotting /= '') then
  if (index('detector', trim(lux_param%plotting)) == 1) then
    call lux_plot_detector (lux_param, lux_com, detec_grid)
  endif
endif

end subroutine lux_write_data

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine lux_plot_detector (lux_param, lux_com, detector)
!
! Routine to plot the intensity at the detector pixels
!-

subroutine lux_plot_detector (lux_param, lux_com, detector)

type (lat_struct), pointer:: lat
type (ele_struct), pointer :: ele
type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (surface_grid_struct) detector

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

do ix = 1, lbound(detector%pt, 1), ubound(detector%pt, 1)
do iy = 1, lbound(detector%pt, 2), ubound(detector%pt, 2)
  call qp_paint_rectangle (x1, x2, y1, y2, color = red$)
enddo
enddo

call read_a_line ('CR to exit: ', ans)

end subroutine lux_plot_detector

end module
