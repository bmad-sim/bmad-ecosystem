!+
! Module: tao_lattice_calc_mod
!
! Lattice calculation routines are here. It's a module so that custom lattice
! calculations has access to the subroutines.
!-

module tao_lattice_calc_mod

use tao_interface
use tao_data_and_eval_mod
use beam_mod
use random_mod
use rad_int_common

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lattice_calc (calc_ok)
!
! For all universes track and calculate the lattice functions for the model lattice.
! Also compute model values for all data.
!
! Output:
!   calc_ok       -- logical: Set False if there was an error in the 
!                     calculation like a particle was lost or a lat is unstable.
!-

subroutine tao_lattice_calc (calc_ok)

use srdt_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_dat
type (tao_d1_data_struct), pointer :: d1_dat
type (coord_struct), allocatable, save :: orb(:)
type (tao_lattice_struct), pointer :: tao_lat
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (aperture_scan_struct), pointer :: scan

real(rp) tt

integer iuni, j, ib, ix, n_max, iu, it, id, ix_ele0

character(20) :: r_name = "tao_lattice_calc"
character(20) track_type, name

logical calc_ok, this_calc_ok, err, mat_changed

! Lattice bookkeeping

calc_ok = .true.

do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)

  if (.not. u%is_on .or. .not. u%calc%lattice) cycle
  tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.
  call tao_lat_bookkeeper (u, tao_lat, err)
  if (err) then
    do id = 1, size(u%data)
      if (u%data(id)%data_type /= 'unstable.lattice') cycle
      call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
    enddo
    calc_ok = .false.
    return
  endif
enddo

! do a custom lattice calculation if desired

if (.not. s%global%lattice_calc_on) return

s%com%ix_ref_taylor = -999   ! Reset taylor map
if (allocated(scratch%spin_map)) deallocate(scratch%spin_map)

call tao_hook_lattice_calc (calc_ok)

if (s%global%track_type /= 'single' .and. s%global%track_type /= 'beam') then
  call out_io (s_error$, r_name, 'UNKNOWN TRACK_TYPE: ' // quote(s%global%track_type), 'DEFAULTING TO "single"')
  s%global%track_type = 'single'
endif

! To save time, s%u(:)%calc%lattice are used to determine what gets calculated. 

uni_loop: do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)
  if (.not. u%is_on) cycle

  if (u%calc%lattice) then
    ! Pointer to appropriate lattice and zero data array

    s%com%lattice_calc_done = .true.
    tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.

    ! Loop over all branches

    branch_loop: do ib = 0, ubound(tao_lat%lat%branch, 1)
   
      branch => tao_lat%lat%branch(ib)
      tao_branch => tao_lat%tao_branch(ib)
      tao_branch%spin_valid = .false.

      u%model%tao_branch(:)%plot_cache_valid = .false.
      u%design%tao_branch(:)%plot_cache_valid = .false.
      u%base%tao_branch(:)%plot_cache_valid = .false.

      call tao_data_coupling_init(branch)

      if (.not. branch%param%live_branch) cycle
      
      do j = 1, 6
        tao_lat%tao_branch(ib)%orbit%vec(j) = 0.0
      enddo

      track_type = s%global%track_type
      if (ib > 0 .and. branch%param%particle == photon$) track_type = 'single'

      ! Even when beam tracking we need to calculate the lattice parameters with single tracking.

      call tao_inject_particle (u, tao_lat, ib)
      call tao_single_track (u, tao_lat, this_calc_ok, ib)
      if (.not. this_calc_ok) calc_ok = .false.
      if (.not. this_calc_ok) exit

      ! Calc radiation integrals when track_type == "beam" since init_beam_distribution may need this.

      if (s%global%rad_int_calc_on .and. &
            (u%calc%rad_int_for_data .or. u%calc%rad_int_for_plotting .or. track_type == 'beam')) then
        call radiation_integrals (tao_lat%lat, tao_branch%orbit, &
                              tao_branch%modes, tao_branch%ix_rad_int_cache, ib, tao_lat%rad_int)
      endif

      !

      if (track_type == 'beam') then
        call tao_inject_beam (u, tao_lat, ib, this_calc_ok, ix_ele0)
        if (.not. this_calc_ok) calc_ok = .false.
        if (.not. this_calc_ok) then
          if (ib == 0) then
            call out_io (s_error$, r_name, 'CANNOT INJECT BEAM. WILL NOT TRACK BEAM...')
          else
            call out_io (s_error$, r_name, &
                'CANNOT INJECT BEAM. WILL NOT TRACK BEAM FOR BRANCH \i0\ ', i_array = [ib])
          endif
          tao_branch%bunch_params(:)%n_particle_lost_in_ele = 0
          tao_branch%bunch_params(:)%n_particle_live = 0
          exit
        endif
        call tao_beam_track (u, tao_lat, this_calc_ok, ib, ix_ele0)
        if (.not. this_calc_ok) calc_ok = .false.
        if (.not. this_calc_ok) exit
      endif

      !

      if (tao_lat%lat%param%geometry == closed$ .and. (u%calc%chrom_for_data .or. u%calc%chrom_for_plotting)) then
        call chrom_calc (tao_lat%lat, s%global%delta_e_chrom, tao_branch%a%chrom, tao_branch%b%chrom, err, &
                    tao_branch%orbit(0)%vec(6), low_E_lat=tao_branch%low_E_lat, high_E_lat=tao_branch%high_E_lat)
      endif

      ! do multi-turn tracking if needed. This is always the main lattice. 

      if (ib == 0) then
        write (name, '(i0, a)') iuni, '@multi_turn_orbit'
        call tao_find_data (err, name, d2_array, print_err = .false.)

        if (size(d2_array) > 0) then
          if (size(d2_array) /= 1) then
            call out_io (s_error$, r_name, 'MULTIPLE D2 DATA ARRAYS ASSOCIATED WITH: ' // name)
          endif
          d2_dat => d2_array(1)%d2
          n_max = 0
          do id = 1, size(d2_dat%d1)
            n_max = max(n_max, ubound(d2_dat%d1(id)%d, 1))
          enddo
          call reallocate_coord (orb, branch%n_ele_max)
          orb(0) = tao_lat%lat%particle_start
          do it = 0, n_max
            do id = 1, size(d2_dat%d1)
              d1_dat => d2_dat%d1(id)
              if (it >= lbound(d1_dat%d, 1) .and. it <= ubound(d1_dat%d, 1)) then
                select case (d1_dat%name)
                case ('x');   d1_dat%d(it)%model_value = orb(0)%vec(1)
                case ('px');  d1_dat%d(it)%model_value = orb(0)%vec(2)
                case ('y');   d1_dat%d(it)%model_value = orb(0)%vec(3)
                case ('py');  d1_dat%d(it)%model_value = orb(0)%vec(4)
                case ('z');   d1_dat%d(it)%model_value = orb(0)%vec(5)
                case ('pz');  d1_dat%d(it)%model_value = orb(0)%vec(6)
                case default
                  call out_io (s_fatal$, r_name, &
                              'BAD MULTI_TURN_ORBIT D1_DATA%NAME: ' // d1_dat%name)
                  call err_exit
                end select
              endif
            enddo
            call track_all (tao_lat%lat, orb, ib)
            orb(0) = orb(branch%n_ele_track)
          enddo
        endif
      endif

      ! Dynamic aperture calc. Only for rings

      if (u%calc%dynamic_aperture .and. u%dynamic_aperture%ix_branch == ib) then  
        if (.not. rf_is_on(branch)) call reallocate_coord (orb, branch%n_ele_track)
        do j = 1, size(u%dynamic_aperture%pz)
          scan => u%dynamic_aperture%scan(j)
          scan%param = u%dynamic_aperture%param
          ! Check for open lattice. Only 1 turn is allowed
          if (branch%param%geometry == open$ .and. (scan%param%n_turn > 1))then
            call out_io (s_fatal$, r_name, 'DYNAMIC APERTURE CALC n_turn > 1 FOR OPEN LATTICE')
            call err_exit       
          endif
            
          if (rf_is_on(branch)) then
            ! pz surrounds the closed orbit
            scan%ref_orb = tao_branch%orbit(0)
            scan%ref_orb%vec(6) = scan%ref_orb%vec(6) + u%dynamic_aperture%pz(j)
          else
            ! If the RF is off, new fixed points will be calculated for various pz
            orb(0)%vec(6) = u%dynamic_aperture%pz(j)
            call  closed_orbit_calc (tao_lat%lat, orb, 4, ix_branch = ib)
            scan%ref_orb = orb(0)
          endif
      
          call out_io (s_info$, r_name, 'Starting Dynamic aperture scan for pz: \f10.6\ ... ', r_array = [u%dynamic_aperture%pz(j)])
          call run_timer ('START')
          call dynamic_aperture_scan(tao_lat%lat%branch(ib), scan, parallel = .true.)
          call run_timer ('READ', tt)
          call out_io (s_info$, r_name, 'Computation time for aperture scan at this energy (min): \f8.2\ ', tt/60)
        enddo
      endif

      ! Note: The SRDT calc does not involve PTC.

      if (u%calc%srdt_for_data > 0) then
        if (u%calc%srdt_for_data >= 2) then
          if (s%global%srdt_use_cache) then
            call srdt_calc_with_cache(tao_lat%lat, tao_branch%srdt, u%calc%srdt_for_data, s%global%srdt_gen_n_slices, s%global%srdt_sxt_n_slices, scratch%srdt_cache)
            if(.not. allocated(scratch%srdt_cache)) s%global%srdt_use_cache = .false. ! there was insufficient memory available.
          else
            call srdt_calc (tao_lat%lat, tao_branch%srdt, u%calc%srdt_for_data, s%global%srdt_gen_n_slices, s%global%srdt_sxt_n_slices)
          endif

        else
          call srdt_calc (tao_lat%lat, tao_branch%srdt, u%calc%srdt_for_data, s%global%srdt_gen_n_slices, s%global%srdt_sxt_n_slices)
        endif
      endif
      
      ! PTC one-turn-map and normal form calc. Only for rings. 

      call tao_ptc_normal_form (u%calc%one_turn_map, tao_lat, ib)

      ! Custom calc.

      call tao_hook_branch_calc (u, tao_lat, branch)

    enddo branch_loop

    ! If calc is on common model then transfer data to base of all other universes

    if (s%com%common_lattice .and. iuni == ix_common_uni$) then
      do j = lbound(s%u, 1), ubound(s%u, 1)
        if (j == ix_common_uni$) cycle
        s%u(j)%data(:)%base_value = u%data(:)%model_value
      enddo
    endif
  endif  ! if (u%calc%lattice)

  ! Calculate non-expression data 

  do id = 1, size(u%data)
    if (u%data(id)%data_type(1:11) == 'expression:') cycle
    call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
  enddo

  ! Mark this lattice as done 

  u%calc%lattice = .false.

  call tao_scale_ping_data(u)

enddo uni_loop

! do any post-processing

do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)
  if (.not. u%is_on) cycle

  do id = 1, size(u%data)
    if (u%data(id)%data_type(1:11) /= 'expression:') cycle
    call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
  enddo
enddo

call tao_hook_post_process_data ()

call tao_set_var_useit_opt
call tao_set_data_useit_opt

end subroutine tao_lattice_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_single_track (u, tao_lat, calc_ok, ix_branch)

use mode3_mod

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct), target :: u
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ele_struct), pointer :: ele
type (beam_init_struct), pointer :: beam_init

real(rp), parameter :: vec0(6) = 0
real(rp) covar, radix, tune3(3), N_mat(6,6), D_mat(6,6), G_inv(6,6)

integer i, ii, n, nn, ix_branch, status, ix_lost, i_dim

character(80) :: lines(10)
character(20) :: r_name = "tao_single_track"

logical calc_ok, err, radiation_fluctuations_on

!

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
              "particle lost in single particle tracking at branch>>element \I0\>>\I0\: " // &
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
      call twiss_at_start (lat, status, branch%ix_branch)
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

! Sigma matric calc.

if (u%calc%beam_sigma_for_data .or. u%calc%beam_sigma_for_plotting) then
  beam_init => u%beam%beam_init
  ele => branch%ele(0)
  if (.not. associated(ele%mode3)) allocate (ele%mode3)
  if (.not. allocated(tao_branch%linear)) allocate(tao_branch%linear(0:branch%n_ele_max))

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
  
    D_mat = 0
    D_mat(1,1) = beam_init%a_emit
    D_mat(2,2) = beam_init%a_emit
    D_mat(3,3) = beam_init%b_emit
    D_mat(4,4) = beam_init%b_emit
    D_mat(5,5) = beam_init%sig_z * beam_init%sig_pz
    D_mat(6,6) = beam_init%sig_z * beam_init%sig_pz

    tao_branch%linear(0)%sigma = matmul(matmul(N_mat, D_mat), transpose(N_mat))

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
      ele%z%alpha = - covar / ele%z%emit
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

    D_mat = 0
    D_mat(1,1) = beam_init%a_emit
    D_mat(2,2) = beam_init%a_emit
    D_mat(3,3) = beam_init%b_emit
    D_mat(4,4) = beam_init%b_emit
    D_mat(5,5) = ele%z%emit
    D_mat(6,6) = ele%z%emit

    tao_branch%linear(0)%sigma = matmul(matmul(N_mat, D_mat), transpose(N_mat))
  endif

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    tao_branch%linear(i)%sigma = matmul(matmul(ele%mat6, tao_branch%linear(i-1)%sigma), transpose(ele%mat6))
  enddo
endif

end subroutine tao_single_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Right now, there is no beam tracking in circular lattices. 
! If extracting from a lat then the beam is generated at the extraction point.

subroutine tao_beam_track (u, tao_lat, calc_ok, ix_branch, ix_ele0)

use wake_mod, only: zero_lr_wakes_in_lat
use beam_utils, only: calc_bunch_params

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (tao_universe_struct), target :: u
type (beam_struct) :: beam
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

integer ix_ele0, what_lat, n_lost_old
integer i, n, i_uni, ip, ig, ic, ie1, ie2, ie
integer n_bunch, n_part, i_uni_to, ix_track
integer n_lost, ix_branch
integer, allocatable, save :: ix_ele(:)

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

ie1 = ix_ele0
ie2 = branch%n_ele_track
if (u%beam%ix_track_end > -1 .and. ix_branch == 0) ie2 = u%beam%ix_track_end

print_err = .true.

if (s%global%beam_timer_on) then
  call run_timer ('START')
  old_time = 0
endif

beam = u%beam%beam_at_start

if (ie2 < ie1 .and. branch%param%geometry == open$) then
  call out_io (s_abort$, r_name, 'BEAM TRACKING WITH STARTING POINT AFTER ENDING POINT IN AN OPEN LATTICE DOES NOT MAKE SENSE.')
  return
endif

n_lost_old = 0
ie = ie1 - 1

do 
  ie = ie + 1
  if (ie == ie2 + 1) exit
  if (ie > branch%n_ele_track) ie = 1

  bunch_params => tao_branch%bunch_params(ie)
  ele => branch%ele(ie)

  ! track to the element and save for phase space plot

  if (s%com%use_saved_beam_in_tracking) then
    beam = tao_model_ele(ie)%beam

  else
    if (ie /= ie1) then
      call track_beam (lat, beam, branch%ele(ie-1), ele, err, centroid = tao_branch%orbit)
      if (err) then
        calc_ok = .false.
        return
      endif
    endif

    can_save = (ie == ie1 .or. ie == ie2 .or. ele%key == fork$ .or. ele%key == photon_fork$)
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
    call out_io (s_error$, r_name, [character(80):: 'Singular sigma matrix is:', &
            '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', '  \6es15.7\', &
            'Will not print any more singular sigma matrices for this track...'], &
            r_array = [sig(1,:), sig(2,:), sig(3,:), sig(4,:), sig(5,:), sig(6,:)])
    print_err = .false.  
  endif

  if (ie == ie1) then
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
character(24) :: r_name = "tao_no_beam_left"

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
type (ele_struct), save :: extract_ele
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
! Subroutine tao_inject_beam (u, model, ix_branch, init_ok, ix_ele0)
!
! This will initialize the beam for a given lattice branch.
!
! Input:
!   u         -- tao_universe_struct: Universe containing the lattice.
!   model     -- tao_lattic_struct: Universe parameters.
!   ix_branch -- integer: Lattice branch index to inject into.
!
! Output:
!   init_ok   -- logical: Set False if there are problems. True otherwise.
!   ix_ele0   -- integer: Index of element injected into in branch. Typically = 0.
!-
subroutine tao_inject_beam (u, model, ix_branch, init_ok, ix_ele0)

use beam_utils
use beam_file_io

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), target :: model
type (tao_model_branch_struct), pointer :: model_branch
type (ele_struct), save :: extract_ele
type (lat_param_struct), pointer :: param
type (branch_struct), pointer :: branch
type (beam_init_struct), pointer :: beam_init
type (beam_struct), pointer :: beam
type (coord_struct), pointer :: orbit

real(rp) v(6)
integer ix_ele0
integer i, j, n, iu, ios, n_in_file, n_in, ix_branch, ib, ie

character(20) :: r_name = "tao_inject_beam"
character(100) line

logical err, init_ok

! If using beam info from a file then no init necessary.

init_ok = .false.

if (s%com%use_saved_beam_in_tracking) then
  init_ok = .true.
  return
endif

! if injecting into a branch then use the branch point as the starting distribution.

model_branch => u%model_branch(ix_branch)
branch => model%lat%branch(ix_branch)

if (ix_branch > 0) then
  ib = branch%ix_from_branch
  ie = branch%ix_from_ele

  if (.not. allocated (u%model_branch(ib)%ele(ie)%beam%bunch)) then
    call out_io (s_error$, r_name, 'CANNOT INJECT INTO BRANCH FROM: ' // u%model%lat%branch(ib)%ele(ie)%name)
    return
  endif

  u%beam%beam_at_start = u%model_branch(ib)%ele(ie)%beam
  init_ok = .true.
  return
endif

! Init for main branch...
! If init is via a call to init_beam_distribution: If the distribution is generated with the help of a random 
! number generator a different distribution is generated on each time init_beam_distribution is called. 
! This is a problem if, say, we are looking at changes to the beam transport due to changes in lattice parameters. 
! Of course if, for example, the beam_init structure is modified then we do want a distribution recalc.
! So only reinit the distribution if the distribution has not already been initialized or if commanded via %init_starting_distribution.

if (u%beam%beam_init%use_particle_start_for_center .and. any(u%beam%beam_init%center /= u%model%lat%particle_start%vec)) then
  u%beam%beam_init%center = u%model%lat%particle_start%vec
  u%beam%init_starting_distribution = .true.
endif

ix_ele0 = u%beam%ix_track_start
if (ix_ele0 == -999) then
  call out_io(s_error$, r_name, 'NOTE: Bad start or stop beam track locations. No beam tracking done.')
  return
endif

beam_init => u%beam%beam_init
beam => u%beam%beam_at_start

if (u%beam%init_starting_distribution .or. .not. allocated(beam%bunch) .or. u%beam%beam_init%position_file /= "") then
  call init_beam_distribution (branch%ele(ix_ele0), branch%param, beam_init, beam, err, &
                                                                 u%model%tao_branch(ix_branch)%modes)
  if (err) then
    call out_io (s_error$, r_name, 'BEAM_INIT INITIAL BEAM PROPERTIES NOT SET FOR UNIVERSE: \i4\ ', u%ix_uni)
    return
  endif

  if (beam%bunch(1)%charge_tot == 0) then
    call out_io (s_warn$, r_name, 'Total beam charge is zero. Set beam_init%bunch_charge to change this.')
  endif

  if (tao_no_beam_left(beam, branch%param%particle)) then
    call out_io (s_warn$, r_name, "Not enough particles or no charge/intensity for beam init!")
    return
  endif

  u%beam%init_starting_distribution = .false.
endif

init_ok = .true.

end subroutine tao_inject_beam
 
end module tao_lattice_calc_mod
