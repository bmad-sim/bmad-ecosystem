module beam_utils

use beam_file_io
use wake_mod

private init_random_distribution, init_grid_distribution
private init_ellipse_distribution, init_kv_distribution
private combine_bunch_distributions, bunch_init_end_calc

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch_hom (bunch, ele, direction, bunch_track)
!
! Subroutine to track a bunch of particles through an element including wakefields.
!
! Input:
!   bunch         -- bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: The element to track through.
!   direction     -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- Bunch_struct: Ending bunch position.
!   bunch_track -- bunch_track_struct, optional: Track information appended to track.
!-

subroutine track1_bunch_hom (bunch, ele, direction, bunch_track)

use ptc_interface_mod, only: ele_to_taylor
use radiation_mod, only: radiation_map_setup
implicit none

type (bunch_struct) bunch
type (ele_struct) ele, half_ele
type (bunch_track_struct), optional :: bunch_track
type (ele_struct), pointer :: wake_ele

type (branch_struct), pointer :: branch

real(rp) charge, ds_wake

integer, optional :: direction
integer i, j, n, jj
logical err_flag, finished, thread_safe

character(*), parameter :: r_name = 'track1_bunch_hom'

! Note: PTC tracking is not not thread safe.

branch => pointer_to_branch(ele)
bunch%particle%direction = integer_option(1, direction)
if (ele%tracking_method == taylor$ .and. .not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele)
thread_safe = (ele%tracking_method /= symp_lie_ptc$ .and. global_com%mp_threading_is_safe)

call save_a_bunch_step (ele, bunch, bunch_track, 0.0_rp)

!------------------------------------------------
! Without wakefields just track through.

wake_ele => pointer_to_wake_ele(ele, ds_wake)
if (.not. associated (wake_ele) .or. (.not. bmad_com%sr_wakes_on .and. .not. bmad_com%lr_wakes_on)) then

  if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) call radiation_map_setup(ele, err_flag)
  !$OMP parallel do if (thread_safe)
  do j = 1, size(bunch%particle)
    if (bunch%particle(j)%state /= alive$) cycle
    call track1 (bunch%particle(j), ele, branch%param, bunch%particle(j))
  enddo
  !$OMP end parallel do

  bunch%charge_live = sum (bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$))
  return
endif

!------------------------------------------------
! This calculation is for an element with wakefields.
! Put the wakefield kicks at the half way point.
! For zero length elements just track the element.

if (ele%value(l$) == 0) then
  if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) call radiation_map_setup(ele, err_flag)
  !$OMP parallel do if (thread_safe)
  do j = 1, size(bunch%particle)
    if (bunch%particle(j)%state /= alive$) cycle
    call track1 (bunch%particle(j), ele, branch%param, bunch%particle(j))
  enddo
  !$OMP end parallel do

else
  call transfer_ele (ele, half_ele, .true.)

  if (integer_option(1, direction) == 1) then
    call create_element_slice (half_ele, ele, ds_wake, 0.0_rp, branch%param, .true., .false., err_flag)
  else
    call create_element_slice (half_ele, ele, ele%value(l$)-ds_wake, ds_wake, branch%param, .false., .true., err_flag, half_ele)
  endif

  if (err_flag) then
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (half_ele%tracking_method == taylor$ .and. .not. associated(half_ele%taylor(1)%term)) call ele_to_taylor(half_ele)

  if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) call radiation_map_setup(half_ele, err_flag)
  !$OMP parallel do if (thread_safe)
  do j = 1, size(bunch%particle)
    if (bunch%particle(j)%state /= alive$) cycle
    call track1 (bunch%particle(j), half_ele, branch%param, bunch%particle(j))
  enddo
  !$OMP end parallel do
endif

! Wakefields

call order_particles_in_z (bunch)  

finished = .false.
if (associated(track1_wake_hook_ptr)) call track1_wake_hook_ptr (bunch, ele, finished)

if (.not. finished) then
  call track1_sr_wake (bunch, wake_ele)
  call track1_lr_wake (bunch, wake_ele)
endif

! Track the last half of the cavity.

if (ele%value(l$) == 0) then
  bunch%charge_live = sum (bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$))
  return
endif

if (integer_option(1, direction) == 1) then
  call create_element_slice (half_ele, ele, ele%value(l$)-ds_wake, ds_wake, branch%param, .false., .true., err_flag, half_ele)
else
  call create_element_slice (half_ele, ele, ds_wake, 0.0_rp, branch%param, .true., .false., err_flag)
endif

if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) call radiation_map_setup(ele, err_flag)
!$OMP parallel do if (thread_safe)
do j = 1, size(bunch%particle)
  if (bunch%particle(j)%state /= alive$) cycle
  call track1 (bunch%particle(j), half_ele, branch%param, bunch%particle(j))
enddo
!$OMP end parallel do

bunch%charge_live = sum (bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$))
if (bunch%charge_live == 0) then
  call out_io (s_info$, r_name, 'Note: Wakes are on but bunch charge is zero!')
endif

call save_a_bunch_step (ele, bunch, bunch_track, ele%value(l$))

end subroutine track1_bunch_hom

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_beam_distribution (ele, param, beam_init, beam, err_flag, modes, beam_init_set, 
!                                                                     print_p0c_shift_warning, conserve_momentum)
!
! Subroutine to initialize a beam of particles. 
! Initialization uses the downstream parameters of ele.
! 
! Note: This routine sets the random number generator according to the settings
! in beam_int and at the end resets things to their initial state.
!
! For more information on individual bunch initialization, see the 
! init_bunch_distribution routine.
! 
! Note: The optional "modes" argument generally is used to pass in normal mode parameters as 
! calculated from the lattice. If present, and if a parameter like beam_init%a_emit are 
! set negative, then the corresponding parameter in the modes structure is used.
! If not present, a warning message is issued and the parameter is set to zero.
! This is only used for parameters that cannot be negative.
!
! Input:
!   ele                 -- Ele_struct: element to initialize distribution at (downstream end).
!   param               -- Lat_param_struct: Lattice parameters
!     %particle              -- Type of particle.
!   beam_init           -- beam_init_struct: Use "getf beam_init_struct" for more details.
!   modes               -- normal_modes_struct, optional: Normal mode parameters. See above.
!   print_p0c_shift_warning -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!   shift_momentum          -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!
! Output:
!   beam                -- Beam_struct: Structure with initialized particles.
!   err_flag            -- logical, optional: Set true if there is an error, false otherwise.
!   beam_init_set       -- beam_init_struct: Set to input beam_init with components like %a_emit set what is used in
!                           constructing the beam (which is different from beam_init%a_emit if this is set negative).
!-

subroutine init_beam_distribution (ele, param, beam_init, beam, err_flag, modes, beam_init_set, &
                                                                      print_p0c_shift_warning, conserve_momentum)
 
use random_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct), target :: beam_init
type (beam_init_struct), optional :: beam_init_set
type (beam_struct), target :: beam
type (normal_modes_struct), optional :: modes
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p

integer i_bunch, i, n, n_kv
logical, optional :: err_flag
logical, optional :: print_p0c_shift_warning, conserve_momentum
logical err_here

character(*), parameter :: r_name = "init_beam_distribution"

!

if (present(err_flag)) err_flag = .true.

if (beam_init%file_name /= '') then
  call out_io(s_abort$, r_name, '"BEAM_INIT%FILE_NAME" SHOULD BE "BEAM_INIT%POSITION_FILE". PLEASE CHANGE.')
  return
endif

! Init from file

if (beam_init%position_file /= '') then
  call read_beam_file (beam_init%position_file, beam, beam_init, err_here, ele, print_p0c_shift_warning, conserve_momentum)
  if (err_here) then
    call out_io (s_abort$, r_name, "PROBLEM READING BEAM POSITION FILE: "// quote(beam_init%position_file))
    return
  endif

  do i_bunch = 1, size(beam%bunch)
    call bunch_init_end_calc(beam%bunch(i_bunch), beam_init, i_bunch-1, ele)
  enddo

  if (present(err_flag)) err_flag = .false.
  return
endif

! Non-file init

call reallocate_beam (beam, max(beam_init%n_bunch, 1), 0)

do i_bunch = 1, size(beam%bunch)
  bunch => beam%bunch(i_bunch)
  call init_bunch_distribution (ele, param, beam_init, i_bunch-1, bunch, err_here, modes, beam_init_set)
  if (err_here) return
enddo
  
if (present(err_flag)) err_flag = .false.

end subroutine init_beam_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_bunch_distribution (ele, param, beam_init, ix_bunch, bunch, err_flag, modes, beam_init_used,
!                                                                          print_p0c_shift_warning, conserve_momentum)
!
! Subroutine to initialize a distribution of particles of a bunch.
! Initialization uses the downstream parameters of ele.
!
! There are four distributions available: 
!   '', or 'ran_gauss' -- Random gaussian distribution.
!   'ellipse'  -- concentric ellipses representing a Gaussian distribution
!   'grid'     -- uniform rectangular grid
!   'KV'       -- Kapchinsky-Vladimirsky distribution
! See the Bmad manual for more information.
!
! The distribution is matched to the Twiss parameters, centroid position, and Energy - z 
! correlation as specified. Coupling in the element ele is incorporated into the distribution.
!
! Note: Except for the random number seed, the random number generator 
! parameters used for this routine are set from the beam_init argument.
! That is, these parameters are independent of what is used everywhere else.
!
! Note: Make sure: |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z
!
! Note: The optional "modes" argument generally is used to pass in normal mode parameters as 
! calculated from the lattice. If present, and if a parameter like beam_init%a_emit are 
! set negative, then the corresponding parameter in the modes structure is used.
! If not present, a warning message is issued and the parameter is set to zero.
! This is only used for parameters that cannot be negative.
!
! Note: To get good results, It is important to make sure that for 
! circular rings that beam_init%center is the correct closed orbit. 
! The closed orbit will shift if, for example, radiation damping is turned on.
!
! Input:
!   ele                 -- Ele_struct: element to initialize distribution at (downstream end).
!   param               -- Lat_param_struct: Lattice parameters
!   beam_init           -- beam_init_struct: Use "getf beam_init_struct" for more details.
!   ix_bunch            -- integer: Bunch index. 0 = bunch generated at time = 0.
!   modes               -- normal_modes_struct, optional: Normal mode parameters. See above.
!   print_p0c_shift_warning -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!   shift_momentum          -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!
! Output:
!   bunch               -- bunch_struct: Structure with initialized particles.
!   err_flag            -- logical, optional: Set True if there is an error. False otherwise.
!   beam_init_used      -- beam_init_struct: Set to input beam_init with components like %a_emit set what is used in
!                           constructing the beam (which can be different from beam_init%a_emit if this is set negative).
!                           If reading from a file, beam_init_used will equal beam_init.
!-

subroutine init_bunch_distribution (ele, param, beam_init, ix_bunch, bunch, err_flag, modes, beam_init_used, &
                                                                                print_p0c_shift_warning, conserve_momentum)

use mode3_mod
use random_mod

implicit none

type (ele_struct) ele
type (ele_struct) twiss_ele
type (lat_param_struct) param
type (beam_struct) beam
type (beam_init_struct), target :: beam_init, beam_init_temp
type (beam_init_struct), optional, target :: beam_init_used
type (beam_init_struct), pointer :: b_init
type (bunch_struct), target :: bunch
type (coord_struct) p_temp
type (coord_struct), pointer :: p
type (kv_beam_init_struct), pointer :: kv
type (normal_modes_struct), optional :: modes

real(rp) beta(3), alpha(3), emit(3), covar, ran(6), center(6), abz_tunes(3)
real(rp) v_mat(4,4), v_inv(4,4), beta_vel
real(rp) old_cutoff
real(rp) tunes(1:3), g6mat(6,6), g6inv(6,6), v6mat(6,6), t6(6,6)

integer ix_bunch
integer i, j, k, n, species
integer :: n_kv     ! counts how many phase planes are of KV type
integer :: ix_kv(3) ! indices (1,2,3) of the two KV planes or 0 if uninitialized

character(16) old_engine, old_converter  
character(*), parameter :: r_name = "init_bunch_distribution"

logical, optional :: err_flag
logical, optional :: print_p0c_shift_warning, conserve_momentum
logical ok, correct_for_coupling(6)
logical ran_gauss_here, err

!

if (present(err_flag)) err_flag = .true.

if (present(beam_init_used)) then
  b_init => beam_init_used
else
  b_init => beam_init_temp
endif
b_init = beam_init

! Read from file?

if (beam_init%position_file /= '') then
  call read_beam_file (beam_init%position_file, beam, beam_init, err, ele, print_p0c_shift_warning, conserve_momentum)
  if (err) then
    call out_io (s_error$, r_name, "PROBLEM READING BEAM POSITION FILE: " // beam_init%position_file)
    return
  endif
  bunch = beam%bunch(1)
  bunch%n_good = 0
  bunch%n_bad = 0

  call bunch_init_end_calc (bunch, beam_init, ix_bunch, ele)

  if (present(err_flag)) err_flag = .false.
  return
endif

! Use modes info if present and parameter is set negative.

species = species_id(beam_init%species)
if (species == not_set$) species = default_tracking_species(param)

b_init = set_emit_from_beam_init (beam_init, ele, species, modes, err); if (err) return

! Save and set the random number generator parameters.

call ran_engine (b_init%random_engine, old_engine)
call ran_gauss_converter (b_init%random_gauss_converter, b_init%random_sigma_cutoff, old_converter, old_cutoff)

! Compute the Twiss parameters beta and alpha, and the emittance for each plane
! 1 = (x,px), 2 = (y,py), 3 = (z,pz)

if (b_init%full_6D_coupling_calc) then
  twiss_ele%a%beta = 1.0d0
  twiss_ele%a%alpha = 0.0d0
  twiss_ele%b%beta = 1.0d0
  twiss_ele%b%alpha = 0.0d0
  twiss_ele%value(e_tot$) = ele%value(e_tot$)
else
  call transfer_ele (ele, twiss_ele, .true.)
endif

covar = b_init%dPz_dz * b_init%sig_z**2
twiss_ele%z%emit = sqrt((b_init%sig_z*b_init%sig_pz)**2 - covar**2)
if (twiss_ele%z%emit == 0 .or. b_init%full_6D_coupling_calc) then
  twiss_ele%z%beta = 1
  twiss_ele%z%alpha = 0
else
  twiss_ele%z%beta = b_init%sig_z**2 / twiss_ele%z%emit
  twiss_ele%z%alpha = - covar / twiss_ele%z%emit
endif

twiss_ele%a%emit = b_init%a_emit
twiss_ele%b%emit = b_init%b_emit

! Init

bunch%n_good = 0
bunch%n_bad = 0

n_kv = 0
ix_kv = 0
ran_gauss_here = .false.
correct_for_coupling = .true.

! Fill the corresponding struct and generate the distribution for each phase plane.
! init_random_distribution must be called last since the other distributions "multiply"
! the number of particles (see the combine_bunch_distribution routine).

call reallocate_bunch (bunch, 0)

do i = 1, 3
  call str_upcase (b_init%distribution_type(i), b_init%distribution_type(i))
  select case (b_init%distribution_type(i))
  case ('', 'RAN_GAUSS')
    ran_gauss_here = .true.
  case ('ELLIPSE')
    call init_ellipse_distribution (i, b_init%ellipse(i), twiss_ele, bunch)
  case ('GRID')
    call init_grid_distribution (i, b_init%grid(i), bunch)
    ! The grid distribution ignores the local twiss parameters.
    correct_for_coupling(2*i-1:2*i) = .false.
  case ('KV') 
    n_kv = n_kv + 1
    ix_kv(n_kv) = i
  case default
    call out_io (s_abort$, r_name, 'PHASE SPACE DISTRIBUTION TYPE NOT RECOGNIZED: '//trim(b_init%distribution_type(i)))
    if (global_com%exit_on_error) call err_exit
    return
  end select
enddo

if (n_kv == 2) call init_KV_distribution (ix_kv(1), ix_kv(2), b_init%kv, twiss_ele, bunch)

if (ran_gauss_here) then
  call init_random_distribution (twiss_ele, species, b_init, bunch, err)
  if (err) return
endif

if (b_init%full_6D_coupling_calc) then
  call transfer_matrix_calc(ele%branch%lat, t6, ix1=ele%ix_ele, one_turn=.true.)
  call calc_z_tune(ele%branch)
  abz_tunes(1) = mod(ele%branch%ele(ele%branch%n_ele_track)%a%phi,twopi)
  abz_tunes(2) = mod(ele%branch%ele(ele%branch%n_ele_track)%b%phi,twopi)
  abz_tunes(3) = ele%branch%z%tune
  call normal_mode3_calc(t6, tunes, G6mat, V6mat, .true., abz_tunes)
  call mat_inverse(G6mat, G6inv, ok = ok)
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    p%vec = matmul(G6inv, p%vec)
  enddo
else
  call make_v_mats(ele, v_mat, v_inv)
endif

bunch%ix_ele = ele%ix_ele
bunch%charge_tot = b_init%bunch_charge
bunch%charge_live = b_init%bunch_charge
bunch%n_live = size(bunch%particle)

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  p%charge = bunch%charge_tot * p%charge
  ! Include Dispersion and coupling
  if (b_init%full_6D_coupling_calc) then
    p%vec(1:6) = matmul(V6mat, p%vec(1:6))

  else
    p_temp%vec(1:4) = p%vec(1:4) + p%vec(6) * [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap]
    p_temp%vec(1:4) = matmul(v_mat, p_temp%vec(1:4))
    where (correct_for_coupling(1:4)) p%vec(1:4) = p_temp%vec(1:4)
  endif

enddo

! Photons:
! For now just give one half e_field_x = 1 and one half e_field_y = 1

if (species == photon$) then
  n = size(bunch%particle)
  bunch%particle(1:n:2)%field(1) = 1/sqrt(2.0_rp)
  bunch%particle(1:n:2)%field(1) = 1/sqrt(2.0_rp)
endif

! End stuff

bunch%particle%species = species
call bunch_init_end_calc(bunch, b_init, ix_bunch, ele)
call init_spin_distribution (b_init, bunch, ele)

! Reset the random number generator parameters.

call ran_engine (old_engine)
call ran_gauss_converter (old_converter, old_cutoff)
if (present(err_flag)) err_flag = .false.
  
end subroutine init_bunch_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_random_distribution (ele, species, beam_init, bunch)
!
! Subroutine to initialize a random bunch of particles matched to
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified. Coupling in the element ele is incorporated into the
! distribution.
!
! Note: This routine is private. Use init_bunch_distribution instead.
!-

subroutine init_random_distribution (ele, species, beam_init, bunch, err_flag, modes)
 
use random_mod

implicit none

type (ele_struct) ele
type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch
type (normal_modes_struct), optional :: modes
type (coord_struct), allocatable :: p(:)
  
real(rp) dpz_dz, denom
real(rp) a_emit, b_emit, y, a, b
real(rp) ave(6), sigma(6), alpha(6), sig_mat(6,6), r(6)
real(rp) center(6), ran_g(2), old_cutoff

integer i, j, j2, n, i_bunch, n_particle, species

logical is_ran_plane(3), err_flag

character(16) old_engine, old_converter  
character(28) :: r_name = "init_random_distribution"

! If random is to be combined with other distributions, the number
! of particles is set by the other distributions.

err_flag = .true.
is_ran_plane = (beam_init%distribution_type == '' .or. beam_init%distribution_type == 'RAN_GAUSS')

n_particle = beam_init%n_particle
if (any(.not. is_ran_plane)) n_particle = size(bunch%particle)

if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
  return
endif

allocate(p(n_particle))

sig_mat = 0
ave = 0
do n = 1, n_particle
  call ran_gauss(r)
  p(n)%vec = r
  ave = ave + r
  forall (i=1:6, j=1:6) sig_mat(i,j) = sig_mat(i,j) + r(i) * r(j)
enddo  

ave = ave / n_particle
sig_mat = sig_mat / n_particle

! Since we are dealing with a finite number of particles, the sigmas of the distributions in 
! each dimension will not be exactly 1 and there will be some correlation between different dimensions.
! If beam_init%renorm_sigma = True then take this out. That is, make sig_mat = the unit matrix.
! Exception: Ignore if n_particle = 1.
! Also: If n_particle < 7 then cannot remove correlations between 
! dimensions so, in this case, only normalize sig_mat(i,i) to 1.

if (beam_init%renorm_sigma) then

  ! Make sure average is zero.

  do n = 1, n_particle
    p(n)%vec = p(n)%vec - ave
  enddo

  forall (i = 1:6, j = 1:6) sig_mat(i,j) = sig_mat(i,j) - ave(i) * ave(j)

  ! The first step is to zero the off-diagonal elements.
  ! We have to do this in the correct order otherwise zeroing one element
  ! might unzero others that have already been zeroed.

  if (n_particle > 6) then
    do i = 5, 1, -1
      do j = i+1, 6

        b = -sig_mat(i,j) / sig_mat(j,j)
        ! Transform the distribution
        do n = 1, n_particle
          p(n)%vec(i) = p(n)%vec(i) + b * p(n)%vec(j)
        enddo
        ! Since we have transformed the distribution we need to transform
        ! sig_mat to keep things consistant.
        sig_mat(i,i) = sig_mat(i,i) + 2 * b * sig_mat(i,j) + b**2 * sig_mat(j,j)
        do j2 = 1, 6
          if (j2 == i) cycle
          sig_mat(i,j2) = sig_mat(i,j2) + b * sig_mat(j ,j2)
          sig_mat(j2,i) = sig_mat(i,j2)
        enddo

      enddo
    enddo
  endif

  ! Now we make the diagonal elements unity

  if (n_particle > 1) then
    do i = 1, 6
      alpha(i) = sqrt(1/sig_mat(i,i))
    enddo

    do n = 1, n_particle
      p(n)%vec = p(n)%vec * alpha
    enddo
  endif

endif

! Compute sigmas

dpz_dz = beam_init%dpz_dz
  
call ran_gauss(ran_g) 
sigma(1) = sqrt(beam_init%a_emit * ele%a%beta)
sigma(2) = sqrt(beam_init%a_emit / ele%a%beta)
sigma(3) = sqrt(beam_init%b_emit * ele%b%beta)
sigma(4) = sqrt(beam_init%b_emit / ele%b%beta)
if (beam_init%full_6D_coupling_calc) then
  sigma(5) = sqrt(beam_init%sig_z*beam_init%sig_pz)
  sigma(6) = sqrt(beam_init%sig_z*beam_init%sig_pz)
else
  sigma(5) = beam_init%sig_z * (1 + beam_init%sig_z_jitter*ran_g(1))
  sigma(6) = beam_init%sig_pz * (1 + beam_init%sig_pz_jitter*ran_g(2))
endif

if (sigma(6) == 0 .or. dpz_dz == 0) then
  a = 0
else if (abs(dpz_dz * sigma(5)) > sigma(6)) then
  call out_io (s_abort$, r_name, "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  if (global_com%exit_on_error) call err_exit
  return
else
  a = dpz_dz * sigma(5) / sigma(6)
endif

b = sqrt(1-a**2)
     
! Put everything together to distribute the particles.

do i = 1, n_particle
  r = p(i)%vec
  p(i)%vec(1) =  sigma(1) *  r(1)
  p(i)%vec(2) = -sigma(2) * (r(2) + r(1) * ele%a%alpha)
  p(i)%vec(3) =  sigma(3) *  r(3)
  p(i)%vec(4) = -sigma(4) * (r(4) + r(3) * ele%b%alpha)
  p(i)%vec(5) =  sigma(5) *  r(5)
  p(i)%vec(6) =  sigma(6) * (r(6) * b + r(5) * a)
end do

! Renormalize the beam centroid

if (beam_init%renorm_center) then
  ave = 0
  do n = 1, n_particle
    ave = ave + p(n)%vec 
  enddo
  ave = ave / n_particle
  do n = 1, n_particle
    p(n)%vec = p(n)%vec - ave
  enddo
endif

! Set particle charge and transfer info the the bunch

p(:)%charge = 1.0_rp / n_particle
call combine_bunch_distributions (bunch, p, is_ran_plane, .false.)
err_flag = .false.

end subroutine init_random_distribution

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_grid_distribution (ix_plane, grid, bunch)
!
! Subroutine to initialize a uniform rectangular grid as the phase space
! distribution of a bunch.
!
! Input:
!    ix_plane       -- Integer: Index of plane of this distribution: 1, 2, or 3
!   grid            -- grid_beam_init_struct: Grid info.
!     %n_x              -- number of columns
!     %n_px             -- number of rows
!     %x_min, %x_max    -- upper and lower limits in beam size
!     %px_min, %px_max  -- upper and lower limits in divergence
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!-

subroutine init_grid_distribution (ix_plane, grid, bunch)

implicit none

type (grid_beam_init_struct) grid
type (coord_struct), allocatable :: p(:)
type (bunch_struct) bunch

integer i, j, k, ix_plane, n_particle

real(rp) x, px

logical where(3)

character(28) :: r_name = 'init_grid_distribution'

!

n_particle = grid%n_x * grid%n_px       ! total number of particles
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
endif

allocate (p(n_particle))

k = 1

do i = 1, grid%n_x
   if (grid%n_x == 1) then
      x = grid%x_min
   else
      x = grid%x_min + real(i - 1)/(grid%n_x - 1) * (grid%x_max - grid%x_min)
   endif

   do j = 1, grid%n_px
      if (grid%n_px == 1) then
         px = grid%px_min
      else
         px = grid%px_min + real(j - 1)/(grid%n_px - 1) * (grid%px_max - grid%px_min)
      endif

      p(k)%vec(2*ix_plane-1) = x
      p(k)%vec(2*ix_plane)   = px
      p(k)%charge = 1.0_rp / n_particle     ! total charge = 1

      k = k + 1
   enddo
enddo

! Combine with bunch distribution

where = .false.;  where(ix_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_grid_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_ellipse_distribution (ix_plane, ellipse, ix_plane, twiss_ele, bunch)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Gaussian distribution.
!
! Input:
!   ellipse           -- ellipse_distribution_struct: Init info.
!     %n_ellipse         -- number of ellipses (>= 1)
!     %part_per_ellipse  -- number of particles per ellipse
!     %sigma_cutoff      -- sigma cutoff of the representation
!   ix_plane          -- Integer: Plane of distribution. 1, 2, or 3.
!   twiss_ele%a       -- Twiss parameters
!   twiss_ele%b       -- Twiss parameters
!   twiss_ele%z       -- Twiss parameters
!   twiss_ele%a%emit  -- emittance
!   twiss_ele%b%emit  -- emittance
!   twiss_ele%z%emit  -- emittance
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!
! See manual for more details.
!-

subroutine init_ellipse_distribution (ix_plane, ellipse, twiss_ele, bunch)

implicit none

type (bunch_struct) bunch
type (coord_struct), allocatable :: p(:)
type (ellipse_beam_init_struct), target :: ellipse
type (ele_struct) twiss_ele
type (ellipse_beam_init_struct), pointer :: e

real(rp) beta, alpha, emit

integer ix_plane, n_particle
integer n, m, k

real(rp) b_inner, b_outer                  ! B_{n-1}/epsilon  and B_{n}/epsilon in the bmad manual

real(rp) J, phi
real(rp) x, px, charge

logical where(3)

character(28) :: r_name = 'init_ellipse_distribution'

!

if (ix_plane == 1) then
  beta = twiss_ele%a%beta
  alpha = twiss_ele%a%alpha
  emit = twiss_ele%a%emit
elseif (ix_plane == 2) then
  beta = twiss_ele%b%beta
  alpha = twiss_ele%b%alpha
  emit = twiss_ele%b%emit
elseif (ix_plane == 3) then
  beta = twiss_ele%z%beta
  alpha = twiss_ele%z%alpha
  emit = twiss_ele%z%emit
endif

e => ellipse
n_particle = e%n_ellipse * e%part_per_ellipse
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
endif

allocate (p(n_particle))

k = 0
b_outer = 0

do n = 1, e%n_ellipse

  b_inner = b_outer
  b_outer = e%sigma_cutoff**2/2.0 * (real(n)/e%n_ellipse)**2

  if (n == e%n_ellipse) then
    ! This is the ellipse that represents the distribution out to infinity
    charge = exp(-b_inner)       ! q_n
    J = emit * (b_inner + 1.0) * exp(-b_inner) / charge   ! J_n
  else
    charge = exp(-b_inner) - exp(-b_outer)   ! q_n
    J = emit * ((b_inner + 1.0) * exp(-b_inner) - (b_outer + 1.0) * exp(-b_outer)) / charge   ! J_n
  endif

  do m = 1, e%part_per_ellipse
    phi = (twopi * m) / e%part_per_ellipse
    k = k + 1
    p(k)%vec(2*ix_plane-1) =  sqrt(2 * J * beta) * cos(phi)
    p(k)%vec(2*ix_plane)   = -sqrt(2 * J / beta) * (alpha * cos(phi) + sin(phi))
    p(k)%charge = charge / e%part_per_ellipse
  enddo

enddo

! Combine with bunch distribution

where = .false.;  where(ix_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_ellipse_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_KV_distribution (ix1_plane, ix2_plane, kv, twiss_ele, bunch)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Kapchinsky-Vladimirsky distribution.
!
! See manual for more details.
!
! Input:
!   ix1_plane        -- Integer: Index of first plane.
!   ix2_plane        -- Integer: Index of second plane.
!   kv               -- kv_beam_init_struct: KV info.
!   twiss_ele%a,b,z  -- Twiss parameters of each phase plane
!   twiss_ele%a%emit -- emittance of a phase plane
!   twiss_ele%b%emit -- emittance of b phase plane
!   twiss_ele%z%emit -- emittance of z phase plane
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!-

subroutine init_KV_distribution (ix1_plane, ix2_plane, kv, twiss_ele, bunch)

implicit none

type (bunch_struct) bunch
type (kv_beam_init_struct) kv
type (ele_struct) twiss_ele
type (coord_struct), allocatable :: p(:)

real(rp) beta1, beta2, alpha1, alpha2, emit1, emit2

integer i_I2, i_phi1, i_phi2, k, n_particle, ix1_plane, ix2_plane, n_p1, n_p2

real(rp) emit_tot
real(rp) I1, I2
real(rp) J1, J2, phi1, phi2
real(rp) x1, x2, px1, px2, charge

logical where(3)

character(28) :: r_name = 'init_kv_distribution'

!

if (ix1_plane == 1) then
  beta1 = twiss_ele%a%beta
  alpha1 = twiss_ele%a%alpha
  emit1 = twiss_ele%a%emit
elseif (ix1_plane == 2) then
  beta1 = twiss_ele%b%beta
  alpha1 = twiss_ele%b%alpha
  emit1 = twiss_ele%b%emit
elseif (ix1_plane == 3) then
  beta1 = twiss_ele%z%beta
  alpha1 = twiss_ele%z%alpha
  emit1 = twiss_ele%z%emit
endif

if (ix2_plane == 1) then
  beta2 = twiss_ele%a%beta
  alpha2 = twiss_ele%a%alpha
  emit2 = twiss_ele%a%emit
elseif (ix2_plane == 2) then
  beta2 = twiss_ele%b%beta
  alpha2 = twiss_ele%b%alpha
  emit2 = twiss_ele%b%emit
elseif (ix2_plane == 3) then
  beta2 = twiss_ele%z%beta
  alpha2 = twiss_ele%z%alpha
  emit2 = twiss_ele%z%emit
endif

n_p1 = kv%part_per_phi(1)
n_p2 = kv%part_per_phi(2)

n_particle = kv%n_i2 * n_p1 * n_p2
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
endif

allocate (p(n_particle))

emit_tot = 1.0 / sqrt(1.0 / emit1**2 + 1.0 / emit2**2)
I1 = kv%A * emit_tot

k = 1

do i_I2 = 1, kv%n_i2
  I2 = -emit1/emit2 * I1 + emit1*emit2/emit_tot**2 * I1 * real(i_I2 - 0.5)/kv%n_i2

  J1 = (I1/emit1 - I2/emit2) * emit_tot
  J2 = (I1/emit2 + I2/emit1) * emit_tot
   
  do i_phi1 = 1, n_p1
    phi1 = 2.0 * pi * real(i_phi1 - 1)/n_p1
 
    do i_phi2 = 1, n_p2
      phi2 = 2.0 * pi * real(i_phi2 - 1)/n_p2

      x1 = sqrt(2.0 * J1 * beta1) * cos(phi1)
      px1 = -sqrt(2.0 * J1 / beta1) * (alpha1 * cos(phi1) + sin(phi1))
      x2 = sqrt(2.0 * J2 * beta2) * cos(phi2)
      px2 = -sqrt(2.0 * J2 / beta2) * (alpha2 * cos(phi2) + sin(phi2))
     
      p(k)%vec(2*ix1_plane-1) = x1
      p(k)%vec(2*ix1_plane)   = px1
      p(k)%vec(2*ix2_plane-1) = x2
      p(k)%vec(2*ix2_plane)   = px2
      p(k)%charge = 1.0_rp / n_particle

      k = k + 1
    enddo
  enddo
enddo

! Combine with bunch distribution

where = .false.;  where(ix1_plane) = .true.;  where(ix2_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_KV_distribution

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine combine_bunch_distributions (bunch, particle, where, do_multiply)
!
! This subroutine combines two bunch distributions together.
!
! If do_multiply is True:
!   size(combined distribution) = size(bunch%particle) * size(particle)
! If do_multiply is False:
!   size(combined distribution) = size(bunch%particle) = size(particle)
!
! Input:
!   bunch       -- bunch_struct: Structure holding the old distribution
!   particle(:) -- coord_struct, allocatable: A new distribution.
!                   This array will be deallocated.
!   where(3)    -- logical: Which planes of particle have the new distribution.
!   do_multiply -- logical: Determines type of combination.
!
! Output:
!   bunch       -- bunch_struct: Structure holding the combined distribution.
!-

subroutine combine_bunch_distributions (bunch, particle, where, do_multiply)

implicit none

type (bunch_struct) bunch
type (coord_struct), allocatable :: particle(:), p(:)

integer i, j, k, m

logical where(:), do_multiply

! If bunch%particle do not contain a distribution, just transfer particle to it.

if (size(bunch%particle) == 0) then
  call reallocate_bunch (bunch, size(particle))
  bunch%particle = particle
  deallocate(particle)
  return
endif

!------------------------
! Multiply combination

if (do_multiply) then
  allocate (p(size(bunch%particle) * size(particle)))
  m = 0
  do i = 1, size(bunch%particle)
    do j = 1, size(particle)
      m = m + 1
      p(m)%charge = bunch%particle(i)%charge * particle(j)%charge
      p(m)%vec = bunch%particle(i)%vec
      do k = 1, 3
        if (.not. where(k)) cycle
        p(m)%vec(2*k-1:2*k) = particle(j)%vec(2*k-1:2*k)
      enddo
    enddo
  enddo

  ! Transfer to bunch

  call reallocate_bunch (bunch, size(p))
  bunch%particle = p
  deallocate (p, particle)

! Overlap combination

else
  do i = 1, size(bunch%particle)
    do k = 1, 3
      if (.not. where(k)) cycle
      bunch%particle(i)%vec(2*k-1:2*k) = particle(i)%vec(2*k-1:2*k)
    enddo
  enddo
  deallocate(particle)
endif

end subroutine combine_bunch_distributions

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_spin_distribution (beam_init, bunch, ele)
!
! Initializes a spin distribution according to beam_init%spin.
!
! Input:
!   beam_init -- beam_init_struct: Initialization parameters 
!     %spin(3)  -- (x, y, z) spin coordinates
!   ele 
!
! Output:
!  bunch    -- bunch_struct: Bunch of particles.
!   %particle(:)%spin
!-

subroutine init_spin_distribution (beam_init, bunch, ele)

implicit none

type (beam_init_struct) beam_init
type (bunch_struct) bunch
type (ele_struct) ele
real(rp) spin(3)
integer i
character(*), parameter :: r_name = 'init_spin_distribution'

!

if (beam_init%use_particle_start) then
  if (.not. associated (ele%branch)) then
    call out_io (s_error$, r_name, 'NO ASSOCIATED LATTICE WITH BEAM_INIT%USE_PARTICLE_START = T.')
    return
  endif
  if (.not. associated (ele%branch%lat)) then
    call out_io (s_error$, r_name, 'NO ASSOCIATED LATTICE WITH BEAM_INIT%USE_PARTICLE_START = T.')
    return
  endif
  spin = ele%branch%lat%particle_start%spin
else
  spin = beam_init%spin
endif

do i = 1, size(bunch%particle)
  bunch%particle(i)%spin = spin
enddo

end subroutine init_spin_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine calc_bunch_params_slice (bunch, bunch_params, plane, slice_center, slice_spread, err, print_err, is_time_coords, ele)
!
! Finds bunch parameters for a slice of the beam.
!
! Input:
!   bunch           -- bunch_struct
!   plane           -- Integer: plane to slice through (x$, px$, & etc...)
!   slice_center    -- real(rp): Center to take slice about
!   slice_spread    -- real(rp): +/- spread in slice about center.
!   print_err       -- logical, optional: If present and False then suppress 
!                       "no eigen-system found" messages.
!   is_time_coords  -- logical, optional: Default is False. If True, input bunch is using time coordinates in which
!                       case there will be a conversion to s-coords before bunch_params are computed.
!   ele             -- ele_struct, optional: Element being tracked through. Must be present if is_time_coords = True.
!
! Output     
!   params -- bunch_params_struct:
!   err    -- Logical: Set True if there is an error.
!-

subroutine calc_bunch_params_slice (bunch, bunch_params, plane, slice_center, slice_spread, err, print_err, is_time_coords, ele)

implicit none

type (bunch_struct) :: bunch
type (bunch_struct) :: sliced_bunch
type (bunch_params_struct) bunch_params
type (ele_struct), optional :: ele

real(rp) slice_center, slice_spread

integer plane
integer i, n_part

logical, optional :: print_err, is_time_coords
logical err

!

n_part = 0

do i = 1, size(bunch%particle)
  if (bunch%particle(i)%vec(plane) > slice_center + abs(slice_spread) .or. &
      bunch%particle(i)%vec(plane) < slice_center - abs(slice_spread)) cycle
  n_part = n_part + 1
enddo

call reallocate_bunch (sliced_bunch, n_part)

sliced_bunch%z_center   = bunch%z_center
sliced_bunch%t_center   = bunch%t_center

n_part = 1
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%vec(plane) > slice_center + abs(slice_spread) .or. &
      bunch%particle(i)%vec(plane) < slice_center - abs(slice_spread)) cycle
  sliced_bunch%particle(n_part) = bunch%particle(i)
  n_part = n_part + 1
enddo

call calc_bunch_params (sliced_bunch, bunch_params, err, print_err, is_time_coords = is_time_coords, ele = ele)

end subroutine calc_bunch_params_slice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine calc_bunch_params_z_slice (bunch, bunch_params, slice_bounds, err, print_err, is_time_coords, ele)
!
! Finds bunch parameters for a slice of the beam.
!
! The slice is specified in terms of percentage of particles ordered by z-position.
! For example, slice_bounds = [0.0, 0.5] specifies the trailing half of the bunch
!
! Input:
!   bunch           -- bunch_struct
!   slice_bounds(2) -- real(rp): Slice bounds in percentage of particles ordered by z-position.
!                         0.0 is the back of the bunch and 1.0 is the front of the bunch.
!   print_err       -- logical, optional: If present and False then suppress 
!                       "no eigen-system found" messages.
!   is_time_coords  -- logical, optional: Default is False. If True, input bunch is using time coordinates in which
!                       case there will be a conversion to s-coords before bunch_params are computed.
!   ele             -- ele_struct, optional: Element being tracked through. Must be present if is_time_coords = True.
!
! Output     
!   params -- bunch_params_struct:
!   err    -- Logical: Set True if there is an error.
!-

subroutine calc_bunch_params_z_slice (bunch, bunch_params, slice_bounds, err, print_err, is_time_coords, ele)

type (bunch_struct) :: bunch
type (bunch_struct) :: sliced_bunch
type (bunch_params_struct) bunch_params
type (ele_struct), optional :: ele

real(rp) slice_bounds(2)
integer i, j, n_part, n0
logical err
logical, optional :: print_err, is_time_coords

!

n0 = nint(slice_bounds(1)*size(bunch%particle))
n_part = nint((slice_bounds(2)-slice_bounds(1))*size(bunch%particle))
call reallocate_bunch(sliced_bunch, n_part)

call indexer (bunch%particle%vec(5), bunch%ix_z)
do i = 1, n_part
  j = bunch%ix_z(i+n0)
  sliced_bunch%particle(i) = bunch%particle(j)
enddo

call calc_bunch_params (sliced_bunch, bunch_params, err, print_err, is_time_coords = is_time_coords, ele = ele)

end subroutine calc_bunch_params_z_slice 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_params (bunch, bunch_params, error, print_err, n_mat, is_time_coords, ele)
!
! Finds all bunch parameters defined in bunch_params_struct, both normal-mode
! and projected. Projected parameters are found purely from the geometrical
! distribution of the beam. Normal-Mode parameters are found using the method
! developed in:
!   "Alternate approach to general coupled linear optics" 
!    A. Wolski, PRST AB 9, 024001 (2006)
!
! Note: If less than two particle remain then the various parameters will be
! set to zero.
! 
! Input:
!   bunch           -- Bunch_struct
!   print_err       -- Logical, optional: If present and False then suppress 
!                       "no eigen-system found" messages.
!   is_time_coords  -- logical, optional: Are particle coords using time coords. Default is False.
!   ele             -- ele_struct, optional: Element being tracked through. Must be present if is_time_coords = True.
!
! Output     
!   bunch_params -- bunch_params_struct:
!   error        -- Logical: Set True if there is an error.
!   n_mat(6,6)   -- real(rp), optional: N matrix defined in Wolski Eq 44 and used to convert 
!                     from action-angle coords to lab coords (Wolski Eq 51.).
!-

subroutine calc_bunch_params (bunch, bunch_params, error, print_err, n_mat, is_time_coords, ele)

implicit none

type (bunch_struct) :: bunch
type (bunch_params_struct) bunch_params
type (ele_struct), optional :: ele

real(rp), optional :: n_mat(6,6)
real(rp) eta, etap
real(rp) :: sigma(6,6) = 0.0
real(rp) :: charge_live, average_pc
real(rp), allocatable :: charge(:)

integer i, j, species

logical, optional :: print_err, is_time_coords
logical error

character(*), parameter :: r_name = "calc_bunch_params"

! Init

bunch_params = bunch_params_struct()
bunch_params%twiss_valid = .false.  ! Assume the worst
error = .true.

call re_allocate (charge, size(bunch%particle))

species = bunch%particle(1)%species

! Get s-position from first live particle

do i = 1, size(bunch%particle)
  if (bunch%particle(i)%state /= alive$) cycle
  bunch_params%ix_ele = bunch%particle(i)%ix_ele
  bunch_params%location = bunch%particle(i)%location
  bunch_params%centroid = bunch%particle(i)  ! Note: centroid%vec gets set by calc_bunch_sigma_matrix
enddo

! n_particle and centroid

bunch_params%n_good_steps = bunch%n_good
bunch_params%n_bad_steps = bunch%n_bad

bunch_params%n_particle_tot = size(bunch%particle)
bunch_params%n_particle_live = count(bunch%particle%state == alive$)
bunch_params%charge_live = sum(bunch%particle%charge, mask = (bunch%particle%state == alive$))
bunch_params%charge_tot = sum(bunch%particle%charge)
bunch%charge_tot  = bunch_params%charge_tot
bunch%charge_live = bunch_params%charge_live

bunch_params%centroid%field(1) = sum(bunch%particle%field(1), mask = (bunch%particle%state == alive$))
bunch_params%centroid%field(2) = sum(bunch%particle%field(2), mask = (bunch%particle%state == alive$))

if (bunch%charge_tot == 0) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'CHARGE OF PARTICLES IN BUNCH NOT SET. CALCULATION CANNOT BE DONE.')
  return
endif

if (species == photon$) then
  charge = bunch%particle%field(1)**2 + bunch%particle%field(2)**2
else
  charge = bunch%particle%charge
endif

charge_live = sum(charge, mask = (bunch%particle%state == alive$))
bunch_params%centroid%charge = charge_live

!

if (bunch_params%n_particle_live == 0) return

if (charge_live == 0) then
  charge = 1
  charge_live = bunch_params%n_particle_live
endif

bunch_params%centroid%s = sum(bunch%particle%s * charge, mask = (bunch%particle%state == alive$)) / charge_live
bunch_params%centroid%t = sum(bunch%particle%t * charge, mask = (bunch%particle%state == alive$)) / charge_live
bunch_params%s = bunch_params%centroid%s
bunch_params%t = bunch_params%centroid%t
bunch_params%sigma_t = sum(bunch%particle%t**2 * charge, mask = (bunch%particle%state == alive$)) / charge_live
bunch_params%sigma_t = sqrt(max(0.0_rp, bunch_params%sigma_t - bunch_params%t**2))

if (bmad_com%spin_tracking_on) call calc_spin_params (bunch, bunch_params)
  
! Sigma matrix calc

call calc_bunch_sigma_matrix_etc (bunch%particle, charge, bunch_params, is_time_coords, ele)

if (species == photon$) return

average_pc = (1+bunch_params%centroid%vec(6)) * bunch_params%centroid%p0c
call convert_pc_to (average_pc, species, beta = bunch_params%centroid%beta)

! Rather arbitrary cutoff: If less than 6 particles, calculation of sigma matrix, etc is declared invalid

if (bunch_params%n_particle_live < 6) return

call calc_emittances_and_twiss_from_sigma_matrix (bunch_params%sigma, bunch_params, error, print_err, n_mat)
if (error .and. logic_option(.true., print_err)) then
  if (present(ele)) call out_io(s_blank$, r_name, 'This at element: ' // ele_full_name(ele))
endif

end subroutine calc_bunch_params

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine calc_emittances_and_twiss_from_sigma_matrix(sigma_mat, bunch_params, error, print_err, n_mat)
!
! Routine to calc emittances and Twiss function from a beam sigma matrix.
! See: Andy Wolski "Alternative approach to general coupled linear optics".
!
! Input:
!   sigma_mat(6,6)  -- real(rp): Sigma matrix.
!   print_err       -- logical, optional: If present and False then suppress 
!                        "no eigen-system found" messages.
!
! Output:
!   bunch_params    -- bunch_params_struct: Holds Twiss and emittance info.
!   error           -- logical: Set True if there is an error. Can happen if the emittance of a mode is zero.
!   n_mat(6,6)      -- real(rp), optional: N matrix defined in Wolski Eq 44 and used to convert 
!                       from action-angle coords to lab coords (Wolski Eq 51.).
!-

subroutine calc_emittances_and_twiss_from_sigma_matrix (sigma_mat, bunch_params, error, print_err, n_mat)

implicit none

type (bunch_params_struct) bunch_params

real(rp) sigma_mat(6,6), sigma_s(6,6), avg_energy, n_real(6,6), f_emit, cut
real(rp), optional :: n_mat(6,6)

complex(rp) :: eigen_val(6) = 0.0, eigen_vec(6,6)
complex(rp) :: n_cmplx(6,6), q(6,6)

integer dim

logical, optional :: print_err
logical error, good(3)

character(*), parameter :: r_name = 'calc_emittances_and_twiss_from_sigma_matrix'

! Init

bunch_params%x = twiss_struct()
bunch_params%y = twiss_struct()
bunch_params%z = twiss_struct()
bunch_params%a = twiss_struct()
bunch_params%b = twiss_struct()
bunch_params%c = twiss_struct()

! X, Y, & Z Projected Parameters

f_emit = bunch_params%centroid%p0c * (1.0_rp + bunch_params%centroid%vec(6)) / mass_of(bunch_params%centroid%species)

call projected_twiss_calc ('X', bunch_params%x, sigma_mat(1,1), sigma_mat(2,2), sigma_mat(1,2), sigma_mat(1,6), sigma_mat(2,6))
call projected_twiss_calc ('Y', bunch_params%y, sigma_mat(3,3), sigma_mat(4,4), sigma_mat(3,4), sigma_mat(3,6), sigma_mat(4,6))
call projected_twiss_calc ('Z', bunch_params%z, sigma_mat(5,5), sigma_mat(6,6), sigma_mat(5,6), 0.0_rp, 0.0_rp)
     
! Normal-Mode Parameters.
! find eigensystem of sigma.S 

sigma_s(:,1) = -sigma_mat(:,2)
sigma_s(:,2) =  sigma_mat(:,1)
sigma_s(:,3) = -sigma_mat(:,4)
sigma_s(:,4) =  sigma_mat(:,3)
sigma_s(:,5) = -sigma_mat(:,6)
sigma_s(:,6) =  sigma_mat(:,5)

! If z or pz has zero sigmas then just do the transverse 

dim = 6
cut = 1e-20_rp * maxval(abs(sigma_mat))
if (abs(sigma_mat(5,5)) < cut .or. abs(sigma_mat(6,6)) < cut) dim = 4  ! No energy oscillations.

call mat_eigen (sigma_s(1:dim,1:dim), eigen_val(1:dim), eigen_vec(1:dim,1:dim), error, print_err)
if (error) return

call normalize_e (eigen_vec(1:dim,1:dim), dim, good, error)
if (error) then
  if (logic_option(.true., print_err)) call out_io (s_warn$, r_name, 'Cannot normalize some eigenvectors.', &
                         'Note: This can happen if the emittance of a normal mode is very small or zero.', &
                         'This will throw off the emittance and other calculations.')
  return
endif

! The eigen-values of Sigma.S are the normal-mode emittances (Wolski Eq. 32)

bunch_params%a%emit = aimag(eigen_val(1))
bunch_params%b%emit = aimag(eigen_val(3))
if (dim == 6) bunch_params%c%emit = aimag(eigen_val(5))

bunch_params%a%norm_emit = bunch_params%a%emit * f_emit
bunch_params%b%norm_emit = bunch_params%b%emit * f_emit
if (dim == 6) bunch_params%c%norm_emit = bunch_params%c%emit * f_emit

! Now find normal-mode sigma matrix and twiss parameters
! Wolski: N = E.Q from Eq. 44 and Eq. 14
! mat_eigen finds row vectors, so switch to column vectors

eigen_vec(1:dim,1:dim) = transpose(eigen_vec(1:dim,1:dim))

q = 0.0
q(1,1) = 1.0/sqrt(2.0)
q(2,1) = 1.0/sqrt(2.0)
q(3,3) = 1.0/sqrt(2.0)
q(4,3) = 1.0/sqrt(2.0)
q(5,5) = 1.0/sqrt(2.0)
q(6,5) = 1.0/sqrt(2.0)
q(1,2) =  i_imag / sqrt(2.0) 
q(2,2) = -i_imag / sqrt(2.0)
q(3,4) =  i_imag / sqrt(2.0) 
q(4,4) = -i_imag / sqrt(2.0)
q(5,6) =  i_imag / sqrt(2.0) 
q(6,6) = -i_imag / sqrt(2.0)

! Compute N in Wolski Eq. 44
n_cmplx(1:dim,1:dim) = matmul(eigen_vec(1:dim,1:dim), q(1:dim,1:dim))  ! Imaginary part is zero.
n_real(1:dim,1:dim) = real(n_cmplx(1:dim,1:dim))

if (present(n_mat)) n_mat = n_real

! Twiss parameters come from Wolski equations 59, 63 and 64

if (good(1)) then
  bunch_params%a%beta = n_real(1,1)**2 + n_real(1,2)**2
  bunch_params%a%alpha = -(n_real(1,1)*n_real(2,1) + n_real(1,2)*n_real(2,2))
  bunch_params%a%gamma = n_real(2,1)**2 + n_real(2,2)**2
endif

if (good(2)) then
  bunch_params%b%beta = n_real(3,3)**2 + n_real(3,4)**2
  bunch_params%b%alpha = -(n_real(3,3)*n_real(4,3) + n_real(3,4)*n_real(4,4))
  bunch_params%b%gamma = n_real(4,3)**2 + n_real(4,4)**2
endif

! Dispersion comes from Wolski equations 69 and 70

if (dim == 6 .and. good(3)) then
  bunch_params%c%beta = n_real(5,5)**2 + n_real(5,6)**2
  bunch_params%c%alpha = -(n_real(5,5)*n_real(6,5) + n_real(5,6)*n_real(6,6))
  bunch_params%c%gamma = n_real(6,5)**2 + n_real(6,6)**2

  bunch_params%a%eta  = n_real(1,5)*n_real(6,5) + n_real(1,6)*n_real(6,6)
  bunch_params%a%etap = n_real(2,5)*n_real(6,5) + n_real(2,6)*n_real(6,6)

  bunch_params%b%eta  = n_real(3,5)*n_real(6,5) + n_real(3,6)*n_real(6,6)
  bunch_params%b%etap = n_real(4,5)*n_real(6,5) + n_real(4,6)*n_real(6,6)

  bunch_params%c%eta  = n_real(5,5)*n_real(6,5) + n_real(5,6)*n_real(6,6)
  bunch_params%c%etap = n_real(6,5)*n_real(6,5) + n_real(6,6)*n_real(6,6)
endif

bunch_params%twiss_valid = .true.

!----------------------------------------------------------------------
contains

subroutine projected_twiss_calc (plane, twiss, exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d)

implicit none

type (twiss_struct) :: twiss

real(rp), intent(in) :: exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d
real(rp) emit, x2, x_px, px2

logical err

character(*) plane

!

if (sigma_mat(6,6) /= 0) then
  twiss%eta   = exp_x_d / sigma_mat(6,6)
  twiss%etap  = exp_px_d / sigma_mat(6,6)
endif

if (sigma_mat(6,6) == 0) then
  x2   = exp_x2 
  x_px = exp_x_px 
  px2  = exp_px2  
else
  x2   = exp_x2 - exp_x_d**2 / sigma_mat(6,6)
  x_px = exp_x_px - exp_x_d * exp_px_d / sigma_mat(6,6)
  px2  = exp_px2  - exp_px_d**2 / sigma_mat(6,6)
endif

twiss%sigma = sqrt(max(0.0_rp, x2))          ! Roundoff may give negative argument.
twiss%sigma_p = sqrt(max(0.0_rp, px2))       ! Roundoff may give negative argument.

emit = sqrt(max(0.0_rp, x2*px2 - x_px**2))   ! Roundoff may give negative argument.

twiss%emit      = emit
twiss%norm_emit = f_emit * emit

if (emit /= 0) then
  twiss%alpha = -x_px / emit
  twiss%beta  = x2 / emit
  twiss%gamma = px2 / emit
endif

end subroutine projected_twiss_calc

!----------------------------------------------------------------------
! contains
! Wolski Eq. 14 But using row vectors to conform to BMAD's mat_eigen

subroutine normalize_e (e, dim, good, err)

implicit none

integer dim
complex(rp) :: s(dim,dim)

complex(rp) :: e(dim,dim), temp(dim)
complex(rp) :: wronsk, factor

integer i2, i, j, k

logical good(3), err

!

err = .false.
good = .false.

s = 0.0
s(1,2) = ( 1.0,0.0) 
s(2,1) = (-1.0,0.0)
s(3,4) = ( 1.0,0.0) 
s(4,3) = (-1.0,0.0)
if (dim == 6) then
  s(5,6) = ( 1.0,0.0) 
  s(6,5) = (-1.0,0.0)
endif

do i2 = 1, dim/2
  i = 2*i2 - 1
  e(i,:) = conjg(e(i+1,:)) ! Eq. 14b
  ! set up the normaization factor
  temp = matmul(s,e(i+1,:))
  wronsk = 0.0
  do j = 1, dim
    wronsk = e(i,j)*temp(j) + wronsk
  enddo
  if (wronsk == 0) then
    err = .true.
    cycle
  endif
  factor = sqrt(i_imag) / sqrt(wronsk)
  ! this next step is the actual normalization (Eq. 14a)
  e(i+1,:) = e(i+1,:) * factor 
  e(i,:) = conjg(e(i+1,:))
  good(i2) = .true.
enddo

end subroutine normalize_e
  
end subroutine calc_emittances_and_twiss_from_sigma_matrix

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine calc_spin_params (bunch, bunch_params)
!
! Rotine to calculate spin averages
!
! Input:
!   bunch -- bunch_struct: Bunch of spins
!
! Output:
!   bunch_params -- bunch_param_struct: Structure holding average
!     centroid%spin(3) -- (x,y,z) polarization.
!-

subroutine calc_spin_params (bunch, bunch_params)

implicit none

type (bunch_params_struct) bunch_params
type (bunch_struct) bunch

real(rp) ave_vec(3), charge_live

integer i

! polarization vector

bunch_params%centroid%spin = 0.0
charge_live = 0

ave_vec = 0.0
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%state /= alive$) cycle
  ave_vec = ave_vec + bunch%particle(i)%spin * bunch%particle(i)%charge
  charge_live = charge_live + bunch%particle(i)%charge
enddo

bunch_params%centroid%spin = ave_vec / charge_live

end subroutine calc_spin_params

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine calc_bunch_sigma_matrix_etc (particle, charge, bunch_params, is_time_coords, ele)
!
! Routine to find the sigma matrix elements of a particle distribution.
! 
! Input:
!   particle(:) -- Coord_struct: Array of particles.
!   charge(:)   -- real(rp): Particle charge or photon intensity.
!
! Output:
!   bunch_params -- bunch_params_struct: Bunch parameters.
!     %sigma(6,6)   
!     %centroid%vec(6)
!     %centroid%p0c
!     %centroid
!     %rel_max(6)
!     %rel_min(6)
!-

subroutine calc_bunch_sigma_matrix_etc (particle, charge, bunch_params, is_time_coords, ele)

implicit none

type (coord_struct) :: particle(:)
type (bunch_params_struct), target :: bunch_params
type (ele_struct), optional :: ele

real(rp) charge_live, p0c_avg, vec(7)
real(rp) charge(:)
real(rp) :: avg(7)

integer i, j2
logical, optional :: is_time_coords
logical is_time

character(*), parameter :: r_name = 'calc_bunch_sigma_matrix'

!

charge_live = sum(charge, mask = (particle%state == alive$))

if (charge_live == 0) then
  call out_io (s_error$, r_name, 'Charge of live particle in bunch is zero! Aborting calculation.')
  return
endif

is_time = logic_option(.false., is_time_coords)
p0c_avg = sum(particle%p0c*charge, mask = (particle%state == alive$)) / charge_live
if (is_time) bunch_params%centroid%p0c = p0c_avg  ! Round-off error is problematic with track1 code if not time coords.

avg = 0
bunch_params%sigma = 0
bunch_params%rel_max = -1e30_rp
bunch_params%rel_min =  1e30_rp

! It is important to use sigma = <r-r_ave>^2 rather than sigma = <r>^2 - r_ave^2 since the latter
! has round-off issues when r_ave is large.

do i = 1, size(particle)
  if (particle(i)%state /= alive$) cycle
  vec = to_basis_coords(particle(i), ele, is_time)
  avg = avg + vec * charge(i)
  bunch_params%rel_max = max(bunch_params%rel_max, vec)
  bunch_params%rel_min = min(bunch_params%rel_min, vec)
enddo

avg = avg / charge_live

do i = 1, size(particle)
  if (particle(i)%state /= alive$) cycle
  vec = to_basis_coords(particle(i), ele, is_time) - avg
  
  forall (j2 = 1:6)
    bunch_params%sigma(:,j2) = bunch_params%sigma(:,j2) + vec(1:6)*vec(j2)*charge(i)
  end forall
enddo

bunch_params%sigma        = bunch_params%sigma / charge_live
bunch_params%rel_max      = bunch_params%rel_max - avg
bunch_params%rel_min      = bunch_params%rel_min - avg
bunch_params%centroid%vec = avg(1:6)

!------------------------
contains

function to_basis_coords (particle, ele, is_time) result (vec)

type (coord_struct) particle, p
type (ele_struct) ele
real(rp) vec(7)
logical is_time

! Time coords uses a different basis for the sigma matrix.

if (.not. is_time) then
  vec = [real(rp):: particle%vec, particle%t]
  return
endif

p = particle
call convert_particle_coordinates_t_to_s(p, ele)
vec = [real(rp):: p%vec(1), p%vec(2) * p%p0c / p0c_avg, p%vec(3), p%vec(4) * p%p0c / p0c_avg, &
                                                p%s, (p%vec(6)*p%p0c + p%p0c - p0c_avg) / p0c_avg, p%t]

end function to_basis_coords

end subroutine calc_bunch_sigma_matrix_etc

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine bunch_init_end_calc (bunch, beam_init, ix_bunch, ele)
!
! Private routine to do the dependent parameter bookkeeping after either reading in particle 
! positions from a beam file or generating positions with init_bunch_distribution.
!
! Input:
!   bunch     -- bunch_struct: Structure with info from the beam file.
!   beam_init -- beam_init_struct: Displace the centroid?
!   ix_bunch  -- integer: Bunch index.
!   ele       -- ele_struct: Lattice element to initalize at.
!
! Output:
!   bunch     -- bunch_struct: Bunch after dependent parameter bookkeeping.
!-

subroutine bunch_init_end_calc (bunch, beam_init, ix_bunch, ele)

implicit none

type (bunch_struct), target :: bunch
type (beam_init_struct) beam_init
type (ele_struct) ele
type (coord_struct), pointer :: p
type (branch_struct), pointer :: branch

real(rp) center(6), ran_vec(6), old_charge, pz_min, t_offset, dz, dz_tot
integer ix_bunch, i, n
character(*), parameter :: r_name = 'bunch_init_end_calc'
logical from_file, h5_file

! Adjust center
! Note: %use_particle_start_for_center is old deprecated name for %use_particle_start.

call ran_gauss(ran_vec)
center = beam_init%center_jitter * ran_vec

if (beam_init%use_particle_start) then
  if (.not. associated (ele%branch)) then
    call out_io (s_error$, r_name, 'NO ASSOCIATED LATTICE WITH BEAM_INIT%USE_PARTICLE_START = T.')
    return
  endif
  if (.not. associated (ele%branch%lat)) then
    call out_io (s_error$, r_name, 'NO ASSOCIATED LATTICE WITH BEAM_INIT%USE_PARTICLE_START = T.')
    return
  endif
  center = center + ele%branch%lat%particle_start%vec
else
  center = center + beam_init%center
endif

!

from_file = (beam_init%position_file /= '')
n = len_trim(beam_init%position_file)
h5_file = (beam_init%position_file(max(1,n-4):n) == '.hdf5' .or. beam_init%position_file(max(1,n-2):n) == '.h5')
t_offset = beam_init%t_offset

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state == alive$) p%ix_turn = beam_init%ix_turn
enddo

if (.not. h5_file) then
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    p%vec = p%vec + center
    p%s = ele%s

    ! Time coordinates. HDF5 files have full particle information so time conversion is ignored.

    if (beam_init%use_t_coords) then
      if (beam_init%use_z_as_t) then
        ! Fixed s, particles distributed in time using vec(5)
        p%t = p%vec(5)
        p%location = downstream_end$
        p%vec(5) = 0 ! init_coord will not complain when beta == 0 and vec(5) == 0
        
      else
        ! Fixed time, particles distributed in space using vec(5)
        p%s = p%vec(5)
        p%t = ele%ref_time + t_offset
        p%location = inside$
      endif

      ! Convert to s coordinates
      p%p0c = ele%value(p0c$)
      call convert_pc_to (sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2), p%species, beta = p%beta)
      p%dt_ref = p%t - ele%ref_time
      call convert_particle_coordinates_t_to_s (p, ele)
      if (.not. from_file) p%state = alive$
      p%ix_ele    = ele%ix_ele
      p%ix_branch = ele%ix_branch

    ! Usual s-coordinates
    else
      call convert_pc_to (ele%value(p0c$) * (1 + p%vec(6)), p%species, beta = p%beta)
      if (from_file .and. p%state /= alive$) then
        p%t = ele%ref_time + t_offset - p%vec(5) / (p%beta * c_light)
        cycle  ! Don't want init_coord to raise the dead.
      endif

      ! If from a file then no vec6 shift needed.
      call init_coord (p, p, ele, downstream_end$, p%species, t_offset = t_offset, shift_vec6 = .not. from_file)

      ! With an e_gun, the particles will have nearly zero momentum (pz ~ -1).
      ! In this case, we need to take care that (1 + pz)^2 >= px^2 + py^2 otherwise, with
      ! an unphysical pz, the particle will be considered to be dead.
      pz_min = 1.000001_rp * sqrt(p%vec(2)**2 + p%vec(4)**2) - 1 ! 1.000001 factor to avoid roundoff problems.
      p%vec(6) = max(p%vec(6), pz_min)
    endif
  enddo
endif

!

bunch%t_center           = ix_bunch * beam_init%dt_bunch
bunch%particle(:)%t      = bunch%particle(:)%t + bunch%t_center
bunch%n_live             = size(bunch%particle)
bunch%charge_live        = sum(bunch%particle%charge)
bunch%ix_bunch           = ix_bunch

dz_tot = 0
do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  dz = -bunch%t_center * c_light * p%beta
  p%vec(5) = p%vec(5) + dz
  dz_tot = dz_tot + dz
enddo
bunch%z_center = dz_tot / bunch%n_live

! If from a file, scale the particle charge.

if (from_file .and. beam_init%bunch_charge /= 0) then
  old_charge = sum(bunch%particle%charge)
  if (old_charge == 0) then
    bunch%particle%charge = beam_init%bunch_charge / size(bunch%particle)
  else
    bunch%particle%charge = bunch%particle%charge * (beam_init%bunch_charge / old_charge)
  endif
  bunch%charge_tot = beam_init%bunch_charge
endif

end subroutine bunch_init_end_calc

end module
