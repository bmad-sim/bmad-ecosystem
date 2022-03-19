!+
! Module: tao_lattice_calc_mod
!
! Lattice calculation routines are here. It's a module so that custom lattice
! calculations has access to the subroutines.
!-

module tao_lattice_calc_mod

use tao_data_and_eval_mod
use beam_mod
use random_mod
use rad_int_common

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_single_track (tao_lat, calc_ok, ix_branch, print_err)
!
! Routine to track a single particle and calculate lattice functions through a lattice.
!
! Input:
!   tao_lat   -- tao_lattice_struct: Structure containing the lattice.
!   ix_branch -- integer: Branch index to track through.
!   print_err -- logical, optional: Default False. Print error messages if, eg, lattice is unstable?
!
! Output:
!   calc_ok   -- logical: Set True if there were no problems, False otherwise.
!-

subroutine tao_single_track (tao_lat, calc_ok, ix_branch, print_err)

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct), pointer :: u
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ele_struct), pointer :: ele

real(rp), parameter :: vec0(6) = 0

integer i, ii, n, nn, ix_branch, status, ix_lost, i_dim

character(80) :: lines(10)
character(*), parameter :: r_name = "tao_single_track"

logical, optional :: print_err
logical calc_ok, err, radiation_fluctuations_on

!

u => tao_lat%u
lat => tao_lat%lat
branch => tao_lat%lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
orbit => tao_branch%orbit

calc_ok = .true.
tao_branch%track_state = moving_forward$

! Track.
! By design, Tao turns off radiation fluctuations (but not damping) for single particle tracking.

ix_lost = branch%n_ele_track + 1

if (u%calc%track) then

  radiation_fluctuations_on = bmad_com%radiation_fluctuations_on
  bmad_com%radiation_fluctuations_on = .false.

  if (branch%param%geometry == closed$) then
    if (rf_is_on(branch)) then
      i_dim = 6
    else
      i_dim = 4
    endif

    call closed_orbit_calc (lat, tao_branch%orbit, i_dim, 1, ix_branch, err_flag = err, print_err = .false.)

    if (err) then
      ! In desperation try a different starting point
      call init_coord (orbit(0), lat%particle_start, branch%ele(0), downstream_end$, orbit(0)%species)
      call closed_orbit_calc (lat, tao_branch%orbit, i_dim, 1, ix_branch, err_flag = err, print_err = .false.)
    endif

    if (err) then
      ! In utmost desperation try a 3rd starting point
      if (i_dim == 4) then
        orbit(0)%vec(1:5) = 0
      else
        orbit(0)%vec = 0
      endif
      call closed_orbit_calc (lat, tao_branch%orbit, i_dim, 1, ix_branch, err_flag = err)
    endif

    if (err) then
      tao_branch%track_state = lost$
      calc_ok = .false.
      if (i_dim == 4) then
        orbit(0)%vec(1:5) = 0
      else
        orbit(0)%vec = 0
      endif
      do i = 1, ubound(orbit, 1)
        orbit(i)%vec = orbit(0)%vec
      enddo
    else
      tao_branch%track_state = moving_forward$
      tao_branch%orb0 = orbit(0)   ! Save beginning orbit as initial guess next time.
    endif

  else
    if (tao_branch%has_open_match_element) then
      do n = 1, branch%n_ele_track
        ele => branch%ele(n)
        call make_mat6(ele, branch%param, orbit(n-1), orbit(n), err_flag = err)
        if (err .or. .not. particle_is_moving_forward(orbit(n))) then
          tao_branch%track_state = n
          exit
        endif
        call twiss_propagate1(branch%ele(n-1), ele)
        if (ele%key == match$) then
          ele%value(match_end_orbit$) = false$
          ele%value(match_end$) = false$
        endif
      enddo

    else
      call track_all (lat, tao_branch%orbit, ix_branch, tao_branch%track_state, &
                                                            orbit0 = tao_lat%tao_branch(0)%orbit)
    endif

    if (tao_branch%track_state /= moving_forward$) then
      calc_ok = .false.
      ix_lost = tao_branch%track_state
      orbit(ix_lost+1:branch%n_ele_track) = coord_struct()
      call out_io (s_blank$, r_name, &
              "tao_single_track: particle lost in single particle tracking at branch>>element \I0\>>\I0\: " // &
              trim(branch%ele(ix_lost)%name) // '  [s =\F9.2\]', &
              r_array = [branch%ele(ix_lost)%s], i_array = [ix_branch, ix_lost])
    endif
  endif

  bmad_com%radiation_fluctuations_on = radiation_fluctuations_on

endif

! Twiss

if (u%calc%twiss .and. branch%param%particle /= photon$) then
  if (tao_branch%track_state == moving_forward$) then
    call lat_make_mat6 (lat, -1, orbit, ix_branch)

    if (branch%param%geometry == closed$) then
      call twiss_at_start (lat, status, branch%ix_branch, print_err)
      if (status /= ok$) then
        calc_ok = .false.
        return
      endif
    endif

    call twiss_propagate_all (lat, ix_branch, err, 0, ix_lost - 1)

  else
    branch%param%stable = .false.
    branch%param%unstable_factor = 0  ! Unknown
  endif
endif

! Calc radiation integrals when track_type == "beam" since init_beam_distribution may need this.

if (branch%param%particle /= photon$ .and. s%global%rad_int_calc_on .and. &
            (u%calc%rad_int_for_data .or. u%calc%rad_int_for_plotting .or. s%global%track_type == 'beam')) then
  call radiation_integrals (tao_lat%lat, tao_branch%orbit, &
                      tao_branch%modes, tao_branch%ix_rad_int_cache, ix_branch, tao_lat%rad_int)
endif

end subroutine tao_single_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lat_sigma_track (tao_lat, calc_ok, ix_branch, print_err)
!
! Routine to track the 6x6 sigma matrix through the lattice using the lattice linear transfer matrices.
!
! Input:
!   tao_lat     -- tao_lattice_struct: Structure containing the lattice.
!   ix_branch   -- integer: Branch index to track through.
!   print_err   -- logical, optional: Default False. Print error messages if, eg, lattice is unstable?
!
! Output:
!   calc_ok     -- logical: Set True if there were no problems, False otherwise.
!-

subroutine tao_lat_sigma_track (tao_lat, calc_ok, ix_branch, print_err)

use mode3_mod

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct), pointer :: u
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (beam_init_struct), pointer :: beam_init
type (tao_model_branch_struct), pointer :: model_branch
type (bunch_params_struct) :: bunch_params

real(rp) covar, radix, tune3(3), N_mat(6,6), D_mat(6,6), G_inv(6,6)

integer ix_branch
integer i, n, ie0, ibf, ief

logical, optional :: print_err
logical calc_ok, err

character(80) :: lines(10)
character(*), parameter :: r_name = "tao_lat_sigma_track"

!

u => tao_lat%u
lat => tao_lat%lat
branch => tao_lat%lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
model_branch => u%model_branch(ix_branch)
ie0 = model_branch%beam%ix_track_start
if (ie0 == not_set$) return

if ((.not. u%calc%lat_sigma_for_data .and. s%com%optimizer_running) .or. branch%param%particle == photon$) return

beam_init => u%model_branch(0)%beam%beam_init
ele => branch%ele(ie0)
if (ele%lord_status == super_lord$) ele => pointer_to_slave(ele, ele%n_slave)
if (.not. associated(ele%mode3)) allocate (ele%mode3)

! Sigma mat at beginning

if (branch%ix_from_branch >= 0) then  ! Propagate through fork
  ibf = branch%ix_from_branch
  ief = branch%ix_from_ele
  tao_branch%lat_sigma(ie0)%mat = tao_lat%tao_branch(ibf)%lat_sigma(ief)%mat

elseif (s%global%init_lat_sigma_from_beam) then
  call calc_bunch_params (model_branch%ele(ie0)%beam%bunch(s%global%bunch_to_plot), bunch_params, err, print_err)
  tao_branch%lat_sigma(ie0)%mat = bunch_params%sigma

else
  call calc_emit_from_beam_init(beam_init, ele, ele%ref_species)
  D_mat = 0
  D_mat(1,1) = ele%a%emit   ! Set by calc_this_emit
  D_mat(2,2) = ele%a%emit
  D_mat(3,3) = ele%b%emit
  D_mat(4,4) = ele%b%emit
  D_mat(5,5) = beam_init%sig_z * beam_init%sig_pz
  D_mat(6,6) = beam_init%sig_z * beam_init%sig_pz

  if (branch%param%geometry == closed$) then
    call transfer_matrix_calc (lat, branch%param%t1_with_RF, ix_branch = ix_branch, one_turn=.true.)
    call make_N (branch%param%t1_with_RF, N_mat, err, tunes_out = tune3)
    if (err) then
      call mat_type (branch%param%t1_with_RF, &
            header = 'SINGULAR ONE-TURN MATRIX WITH RF. WILL NOT BE ABLE TO COMPUTE SIGMAS.', &
            lines = lines, n_lines = n)
      call out_io (s_error$, r_name, lines(1:n))
      return
    endif

    tao_branch%lat_sigma(ie0)%mat = matmul(matmul(N_mat, D_mat), transpose(N_mat))

  else
    covar = beam_init%dPz_dz * beam_init%sig_z**2
    radix = (beam_init%sig_z * beam_init%sig_pz)**2 - covar**2
    if (radix < 0) then
      call out_io (s_error$, r_name, 'Beam init: |dPz_dz| must be less than Sig_Pz/Sig_z')
      calc_ok = .false.
      return
    endif

    ele%z%emit  = sqrt(radix)

    if (ele%z%emit == 0) then
      ele%z%beta  = 1
      ele%z%alpha = 0
      ele%z%gamma = 1
    else
      ele%z%beta  = beam_init%sig_z**2 / ele%z%emit
      ele%z%alpha = -covar / ele%z%emit
      ele%z%gamma = (1 + ele%z%alpha**2) / ele%z%beta
    endif

    call twiss3_from_twiss2 (branch%ele(0))

    G_inv = 0
    G_inv(1,1) = sqrt(ele%a%beta)
    G_inv(2,1:2) = [-ele%a%alpha / sqrt(ele%a%beta), 1/sqrt(ele%a%beta)]
    G_inv(3,3) = sqrt(ele%b%beta)
    G_inv(4,3:4) = [-ele%b%alpha / sqrt(ele%b%beta), 1/sqrt(ele%b%beta)]
    G_inv(5,5) = sqrt(ele%z%beta)
    G_inv(6,5:6) = [-ele%z%alpha / sqrt(ele%z%beta), 1/sqrt(ele%z%beta)]

    N_mat = matmul(ele%mode3%v, G_inv)

    tao_branch%lat_sigma(ie0)%mat = matmul(matmul(N_mat, D_mat), transpose(N_mat))
  endif
endif

! And propagate

i = ie0
do 
  i = i + 1
  if (i > branch%n_ele_track) then
    if (ie0 == 0) exit
    tao_branch%lat_sigma(0)%mat = tao_branch%lat_sigma(branch%n_ele_track)%mat
    i = 1
  endif
  if (i == ie0) exit
  ele => branch%ele(i)
  tao_branch%lat_sigma(i)%mat = matmul(matmul(ele%mat6, tao_branch%lat_sigma(i-1)%mat), transpose(ele%mat6))
enddo

end subroutine tao_lat_sigma_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_beam_track (u, tao_lat, ix_branch, beam, calc_ok)
!
! Routine to track a a beam of particles.
!
! Input:
!   u         -- tao_universe_struct: Universe to track through.
!   tao_lat   -- tao_lattice_struct: Structure containing the lattice.
!   ix_branch -- integer: Branch index to track through.
!   beam      -- beam_struct: Initial beam distribution
!   print_err -- logical, optional: Default False. Print error messages if, eg, lattice is unstable?
!
! Output:
!   beam      -- beam_struct: Final beam distribution.
!   calc_ok   -- logical: Set True if there were no problems, False otherwise.
!-

subroutine tao_beam_track (u, tao_lat, ix_branch, beam, calc_ok)

use wake_mod, only: zero_lr_wakes_in_lat
use beam_utils, only: calc_bunch_params

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (tao_universe_struct), target :: u
type (beam_struct) :: init_beam, beam
type (normal_modes_struct) :: modes
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_model_element_struct), pointer :: tao_model_ele(:)
type (tao_model_branch_struct), pointer :: model_branch
type (bunch_params_struct), pointer :: bunch_params

real(rp) sig(6,6)

integer what_lat, n_lost_old, ie_start, ie_end
integer i, n, i_uni, ip, ig, ic, ie
integer n_bunch, n_part, i_uni_to, ix_track
integer n_lost, ix_branch
integer, allocatable :: ix_ele(:)

character(*), parameter :: r_name = "tao_beam_track"

real(rp) :: value1, value2, f, time, old_time

logical calc_ok, print_err, err, lost, new_beam_file, can_save

! Initialize 

call re_allocate (ix_ele, 1)

s%com%have_tracked_beam = .true.
branch => tao_lat%lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
model_branch => u%model_branch(ix_branch)
tao_model_ele => model_branch%ele

lat => tao_lat%lat

tao_branch%track_state = moving_forward$  ! Needed by tao_evaluate_a_datum
ix_track = moving_forward$
lost = .false.
calc_ok = .true.
new_beam_file = .true.

ie_start = model_branch%beam%ix_track_start
ie_end   = model_branch%beam%ix_track_end

if (ie_start == not_set$ .or. ie_end == not_set$) then
  call out_io (s_error$, r_name, 'BEAM END POSITION NOT PROPERLY SET. NO TRACKING DONE.')
  return
endif

tao_branch%bunch_params(:)%n_particle_lost_in_ele = 0
tao_branch%bunch_params(:)%n_particle_live = 0
tao_branch%bunch_params(:)%n_particle_live = 0
tao_branch%bunch_params(:)%twiss_valid = .false.

! Transfer wakes from  design

do i = 1, branch%n_ele_max
  if (associated(branch%ele(i)%wake)) branch%ele(i)%wake%lr%mode = u%design%lat%branch(ix_branch)%ele(i)%wake%lr%mode
enddo

call zero_lr_wakes_in_lat (lat)

! track through the lattice elements.
! The reference orbit for the transfer matrix calculation is taken to be the nominal 
! bunch center rather then the computed center to prevent jitter from changing things.

print_err = .true.

if (s%global%beam_timer_on) then
  call run_timer ('START')
  old_time = 0
endif

if (ie_end < ie_start .and. branch%param%geometry == open$) then
  call out_io (s_abort$, r_name, 'BEAM TRACKING WITH STARTING POINT AFTER ENDING POINT IN AN OPEN LATTICE DOES NOT MAKE SENSE.')
  return
endif

n_lost_old = 0
ie = ie_start - 1

do 
  ie = ie + 1
  if (ie > branch%n_ele_track) ie = 1

  bunch_params => tao_branch%bunch_params(ie)
  ele => branch%ele(ie)

  ! track to the element and save for phase space plot

  if (s%com%use_saved_beam_in_tracking) then
    beam = tao_model_ele(ie)%beam

  else
    if (ie /= ie_start) then
      call track_beam (lat, beam, branch%ele(ie-1), ele, err, centroid = tao_branch%orbit)
      if (err) then
        calc_ok = .false.
        return
      endif
    endif

    can_save = (ie == ie_start .or. ie == ie_end .or. ele%key == fork$ .or. ele%key == photon_fork$)
    if (tao_model_ele(ie)%save_beam_internally .or. can_save) tao_model_ele(ie)%beam = beam

    if (u%beam%dump_file /= '' .and. tao_model_ele(ie)%save_beam_to_file) then
      if (index(u%beam%dump_file, '.h5') == 0 .and. index(u%beam%dump_file, '.hdf5') == 0) then
        call write_beam_file (u%beam%dump_file, beam, new_beam_file, ascii$, lat)
      else
        call write_beam_file (u%beam%dump_file, beam, new_beam_file, hdf5$, lat)
      endif
      new_beam_file = .false.
    endif
  endif
 
  ! Lost particles

  n_bunch = s%global%bunch_to_plot
  n_lost = count(beam%bunch(n_bunch)%particle(:)%state /= alive$)
  if (n_lost /= n_lost_old) then
    n = size(beam%bunch(n_bunch)%particle(:))
    if (size(s%u) == 1) then
      call out_io (s_blank$, r_name, &
            '\i0\ particle(s) lost at element ' // trim(ele_loc_name(ele)) // ': ' // trim(ele%name) // &
            '  Total lost: \i0\  of \i0\ ', &
            i_array = [n_lost-n_lost_old, n_lost, n])
    else
      call out_io (s_blank$, r_name, &
            '\i0\ particle(s) lost in universe \i0\ at element ' // trim(ele_loc_name(ele)) // &
                                                                           ': ' // trim(ele%name) // &
            '  Total lost: \i0\  of \i0\ ', &
            i_array = [n_lost-n_lost_old, u%ix_uni, n_lost, n])
    endif
    n_lost_old = n_lost
  endif

  if (tao_no_beam_left(beam, branch%param%particle) .and. .not. lost) then
    ix_track = ie
    lost = .true.
    call out_io (s_warn$, r_name, &
            'TOO MANY PARTICLES HAVE BEEN LOST AT ELEMENT ' // trim(ele_loc_name(ele)) // &
                                          ': ' // trim(ele%name))
  endif

  ! calc bunch params

  call calc_bunch_params (beam%bunch(s%global%bunch_to_plot), bunch_params, err, print_err)
  ! Only generate error message once per tracking
  if (err .and. print_err) then
    sig = bunch_params%sigma
    if (any(sig /= 0)) then
      call out_io (s_warn$, r_name, &
            'Beam parameters not computed at: ' // trim(ele%name) // '  ' // ele_loc_name(ele, parens = '()') , &
            '[This will happen with round beams. Ignore this problem if the beam parameters at problem locations are not needed.]', &
            'The singular sigma matrix is:', &
            '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', &
            'Will not print any more singular sigma matrices for this track...', &
            r_array = [sig(1,:), sig(2,:), sig(3,:), sig(4,:), sig(5,:), sig(6,:)])
      print_err = .false.  
    endif
  endif

  if (ie == ie_start) then
    bunch_params%n_particle_lost_in_ele = 0
  else
    bunch_params%n_particle_lost_in_ele = tao_branch%bunch_params(ie-1)%n_particle_live - &
                                                     bunch_params%n_particle_live
  endif

  ! Timer 

  if (s%global%beam_timer_on) then
    call run_timer ('READ', time)
    if (time - old_time > 60) then
      call out_io (s_blank$, r_name, 'Beam at Element: \i0\. Time: \i0\ min', &
                          i_array = [ie, nint(time/60)])
      old_time = time
    endif
  endif

  if (ie == ie_end) exit
enddo

! only post total lost if no extraction or extracting to a turned off lattice

n_lost = 0
do n_bunch = 1, size(beam%bunch)
  n_lost = n_lost + count(beam%bunch(n_bunch)%particle%state /= alive$)
enddo
if (n_lost /= 0) &
  call out_io (s_blank$, r_name, &
      "Total number of lost particles by the end of universe \I2\: \I5\.", &
                                  i_array = [u%ix_uni, n_lost])

tao_branch%track_state = ix_track
 
end subroutine tao_beam_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

function tao_no_beam_left (beam, particle) result (no_beam)

implicit none

type (beam_struct), target :: beam
type (coord_struct), pointer :: p(:)

real(rp) charge_tot
integer n_bunch, particle
logical no_beam, all_lost
character(*), parameter :: r_name = "tao_no_beam_left"

!

no_beam = .false.

n_bunch = s%global%bunch_to_plot
p =>beam%bunch(n_bunch)%particle(:)
all_lost = all(p%state /= alive$)

if (particle == photon$) then
  if (sum(p%field(1)**2) + sum(p%field(2)**2) == 0 .or. all_lost) no_beam = .true.
else
  if (sum(p%charge) == 0 .or. all_lost) no_beam = .true.
endif

end function tao_no_beam_left

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a particle from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated. If no injection then will set beginning orbit to
! whatever the user has specified. As always, tracking only occures in the model
! lattice.

subroutine tao_inject_particle (u, model, ix_branch)

implicit none

type (tao_universe_struct) u
type (tao_lattice_struct), target :: model
type (ele_struct) :: extract_ele
type (ele_struct), pointer :: from_ele, ele0
type (coord_struct) pos
type (coord_struct), pointer :: orb_out, orb_in
type (branch_struct), pointer :: branch, branch_from

integer ix_branch, i_ele_from, i_br_from

character(*), parameter :: r_name = "tao_inject_particle"

! If injecting a particle from another branch.

branch => model%lat%branch(ix_branch)
i_br_from  = branch%ix_from_branch

if (i_br_from > -1) then
  i_ele_from = branch%ix_from_ele
  branch_from => model%lat%branch(i_br_from)
  from_ele => branch_from%ele(i_ele_from)
  ele0 => branch%ele(0)
  ! Note: reference energy set by lat_compute_ref_energy_and_time
  if (is_true(ele0%value(inherit_from_fork$)) .and. branch%ix_to_ele == 0) call transfer_twiss(from_ele, ele0)
  orb_out => model%tao_branch(ix_branch)%orbit(0)
  call init_coord (orb_out, model%tao_branch(i_br_from)%orbit(i_ele_from), ele0, &
                  downstream_end$, default_tracking_species(branch%param), 1, ele0%value(E_tot$))
  if (nint(from_ele%value(direction$)) == -1) then
    if (orb_out%species == photon$) then
      orb_out%vec(2:6:2) = -orb_out%vec(2:6:2)
    else
      orb_out%vec(2:4:2) = -orb_out%vec(2:4:2)
    endif
  endif
  return
endif

! Not injecting from another branch.
! In model%tao_branch()%orb0 is saved the last computed orbit. 
! This is important when using a common_lattice since tao_lat%tao_branch()%orbit(0) has been overwritten.

orb_out => model%tao_branch(ix_branch)%orbit(0)

if (branch%param%geometry == closed$) then
  orb_in => model%tao_branch(ix_branch)%orb0
  if (.not. rf_is_on(branch)) orb_in%vec(6) = model%lat%particle_start%vec(6)
else
  orb_in => branch%lat%particle_start
endif

call init_coord (orb_out, orb_in, branch%ele(0), downstream_end$, default_tracking_species(branch%param), 1, orb_in%p0c)

end subroutine tao_inject_particle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subroutine tao_inject_beam (u, model, ix_branch, beam, init_ok)
!
! This will initialize the beam for a given lattice branch.
!
! Trying to inject a beam of one species into a branch with a different ref species
! (example: electron bunch into photon branch) is problematical. To avoid problems, Tao
! will set not inject (init_ok = False) if there is a mismatch.
!
! Input:
!   u         -- tao_universe_struct: Universe containing the lattice.
!   model     -- tao_lattice_struct: Universe parameters.
!   ix_branch -- integer: Lattice branch index to inject into.
!
! Output:
!   beam          -- beam_struct: Initial beam.
!   init_ok       -- logical: Set False if there are problems. True otherwise.
!-

subroutine tao_inject_beam (u, model, ix_branch, beam, init_ok)

use beam_utils
use beam_file_io

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), target :: model
type (tao_model_branch_struct), pointer :: model_branch
type (branch_struct), pointer :: branch
type (beam_struct) :: beam
type (ele_struct), pointer :: ele0
type (coord_struct), pointer :: orbit
type (tao_beam_branch_struct), pointer :: bb

real(rp) v(6)
integer i, j, n, iu, ios, n_in_file, n_in, ix_branch, ib0, ie0, ie_start

character(*), parameter :: r_name = "tao_inject_beam"
character(100) line

logical err, init_ok

!

init_ok = .false.

model_branch => u%model_branch(ix_branch)
bb => model_branch%beam
branch => model%lat%branch(ix_branch)
ie_start = bb%ix_track_start
if (ie_start == not_set$) then
  call out_io (s_error$, r_name, 'BEAM STARTING POSITION NOT PROPERLY SET. NO TRACKING DONE.')
  return
endif

! If using beam info from a file then no init necessary.

if (s%com%use_saved_beam_in_tracking) then
  init_ok = .true.
  beam = model_branch%ele(ie_start)%beam
  u%model_branch(ix_branch)%ele(ie_start)%beam = beam
  return
endif

! if injecting into a branch then use the branch point as the starting distribution.

if (branch%ix_from_branch >= 0) then  ! Injecting from other branch
  ib0 = branch%ix_from_branch
  ie0 = branch%ix_from_ele
  ele0 => u%model%lat%branch(ib0)%ele(ie0)

  if (.not. allocated (u%model_branch(ib0)%ele(ie0)%beam%bunch)) then
    call out_io (s_error$, r_name, 'CANNOT INJECT INTO BRANCH: ' // int_str(ix_branch) // &
                                   ' FROM: ' // ele0%name // ' IN UNIVERSE: ' // int_str(u%ix_uni))
    return
  endif

  ! How to handle injection from one branch with one species into another branch with a different
  ! species is not well defined. In this case, don't inject.

  beam = u%model_branch(ib0)%ele(ie0)%beam
  u%model_branch(ix_branch)%ele(ie_start)%beam = beam

  if (beam%bunch(1)%particle(1)%species /= u%model%lat%branch(ix_branch)%param%particle) return
  init_ok = .true.

  return
endif

! Init for a root branch...
! If init is via a call to init_beam_distribution: If the distribution is generated with the help of a random 
! number generator, a different distribution is generated each time init_beam_distribution is called. 
! This is a problem if, say, we are looking at changes to the beam transport due to changes in lattice parameters. 
! Of course if, for example, the beam_init structure is modified then we do want a distribution recalc.

bb => u%model_branch(ix_branch)%beam

if (bb%beam_init%use_particle_start .and. any(bb%beam_init%center /= u%model%lat%particle_start%vec)) then
  bb%beam_init%center = u%model%lat%particle_start%vec
  bb%beam_init%spin   = u%model%lat%particle_start%spin
  bb%init_starting_distribution = .true.
endif

ele0 => branch%ele(bb%ix_track_start)

if (bb%init_starting_distribution .or. u%beam%always_reinit) then
  call init_beam_distribution (ele0, branch%param, bb%beam_init, beam, err, u%model%tao_branch(ix_branch)%modes)
  if (err) then
    call out_io (s_error$, r_name, 'BEAM_INIT INITIAL BEAM PROPERTIES NOT PROPERLY SET FOR UNIVERSE: ' // int_str(u%ix_uni))
    return
  endif

  if (beam%bunch(1)%charge_tot == 0) then
    call out_io (s_warn$, r_name, 'Total beam charge is zero. Set beam_init%bunch_charge to change this.')
  endif

  if (tao_no_beam_left(beam, branch%param%particle)) then
    call out_io (s_warn$, r_name, "NOT ENOUGH PARTICLES OR NO CHARGE/INTENSITY FOR BEAM INITIALIZATION IN BRANCH " // int_str(ix_branch))
    return
  endif

  bb%init_starting_distribution = .false.
  bb%beam_at_start = beam

else
  beam = bb%beam_at_start
endif

u%model_branch(ix_branch)%ele(ie_start)%beam = beam
init_ok = .true.

end subroutine tao_inject_beam
 
end module tao_lattice_calc_mod
