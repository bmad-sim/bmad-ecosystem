!+
! Module: tao_lattice_calc_mod
!
! Lattice calculation routines are here. It's a module so that custom lattice
! calculations has access to the subroutines.
!-

module tao_lattice_calc_mod

use tao_mod
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
! Routine to calculate the lattice functions and TAO data. 
! Always tracks through the model lattice. 
! If initing the design lattices then the tracking will still be done
! through the model lattices. tao_init then transfers this into the
! design lattices.
!
! Output:
!   calc_ok -- Logical: Set False if there was an error in the 
!                calculation like a particle was lost or a lat is unstable.
!-

subroutine tao_lattice_calc (calc_ok)

use ptc_layout_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_dat
type (tao_d1_data_struct), pointer :: d1_dat
type (coord_struct), allocatable, save :: orb(:)
type (tao_lattice_struct), pointer :: tao_lat
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: lat_branch
type (normal_form_struct), pointer :: normal_form
type (aperture_scan_struct), pointer :: scan

real(rp) tt
integer iuni, j, ib, ix, n_max, iu, it, id, ix_ele0

character(20) :: r_name = "tao_lattice_calc"
character(20) track_type, name

logical calc_ok, this_calc_ok, err, rf_on

!

calc_ok = .true.

! do a custom lattice calculation if desired

if (.not. s%global%lattice_calc_on) return

s%com%ix_ref_taylor = -999   ! Reset taylor map
call tao_hook_lattice_calc (calc_ok)
    
! To save time, s%u(:)%calc%lattice are used to determine what gets calculated. 

uni_loop: do iuni = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(iuni)
  if (.not. u%is_on .or. .not. u%calc%lattice) cycle

  ! Pointer to appropriate lattice and zero data array

  tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.
  u%data(:)%good_model = .false. ! reset
  u%data%model_value = tiny(1.0_rp)

  ! Lattice bookkeeping

  call tao_lat_bookkeeper (u, tao_lat)

  ! Loop over all branches

  u%info%lat_len_tot = 0

  branch_loop: do ib = 0, ubound(tao_lat%lat%branch, 1)
 
    branch => tao_lat%lat%branch(ib)
    lat_branch => tao_lat%lat_branch(ib)
    
    u%info%lat_len_tot = u%info%lat_len_tot + branch%param%total_length

    call tao_data_coupling_init(branch)

    do j = 1, 6
      tao_lat%lat_branch(ib)%orbit%vec(j) = 0.0
    enddo

    track_type = s%global%track_type
    if (ib > 0 .and. branch%param%particle == photon$) track_type = 'single'

    select case (track_type)
    case ('single') 
      call tao_inject_particle (u, tao_lat, ib)
      call tao_single_track (u, tao_lat, this_calc_ok, ib)
      if (.not. this_calc_ok) calc_ok = .false.
      if (.not. this_calc_ok) exit

    case ('beam')  ! Even when beam tracking we need to calculate the lattice parameters.
      call tao_inject_particle (u, tao_lat, ib)
      call tao_single_track (u, tao_lat, this_calc_ok, ib)
      call tao_inject_beam (u, tao_lat, ib, this_calc_ok, ix_ele0)
      if (.not. this_calc_ok) calc_ok = .false.
      if (.not. this_calc_ok) then
        if (ib == 0) then
          call out_io (s_error$, r_name, 'CANNOT INJECT BEAM. WILL NOT TRACK BEAM...')
        else
          call out_io (s_error$, r_name, &
              'CANNOT INJECT BEAM. WILL NOT TRACK BEAM FOR BRANCH \i0\ ', i_array = [ib])
        endif      
        tao_lat%lat_branch(ib)%bunch_params(:)%n_particle_lost_in_ele = 0
        tao_lat%lat_branch(ib)%bunch_params(:)%n_particle_live = 0
        exit
      endif
      call tao_beam_track (u, tao_lat, this_calc_ok, ib, ix_ele0)
      if (.not. this_calc_ok) calc_ok = .false.
      if (.not. this_calc_ok) exit

    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING TYPE: ' // track_type)
      call err_exit
    end select

    ! Radiation integrals. At some point should extend this to non-main branches.

    if (ib == 0) then
      if (u%calc%rad_int_for_data .or. u%calc%rad_int_for_plotting) then
        call radiation_integrals (tao_lat%lat, tao_lat%lat_branch(ib)%orbit, &
                              tao_lat%modes, tao_lat%ix_rad_int_cache, ib, tao_lat%rad_int)
      endif
    endif

    if (u%calc%chrom_for_data .or. u%calc%chrom_for_plotting) then
      call chrom_calc (tao_lat%lat, s%global%delta_e_chrom, tao_lat%a%chrom, &
                           tao_lat%b%chrom, err, low_E_lat=tao_lat%low_E_lat, high_E_lat=tao_lat%high_E_lat)
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
        call reallocate_coord (orb, tao_lat%lat%n_ele_max)
        orb(0) = tao_lat%lat%beam_start
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
          call track_all (tao_lat%lat, orb)
          orb(0) = orb(tao_lat%lat%n_ele_track)
        enddo
      endif
    endif

    ! Dynamic aperture calc. Only for rings

    if (u%calc%dynamic_aperture) then  
      if (.not. rf_is_on(tao_lat%lat%branch(0))) call reallocate_coord (orb, tao_lat%lat%n_ele_track)
      do j=1, size(u%dynamic_aperture%pz)
        scan => u%dynamic_aperture%scan(j)
        ! Check for open lattice. Only 1 turn is allowed
        if (branch%param%geometry == open$ .and. (scan%param%n_turn > 1 ))then
          call out_io (s_fatal$, r_name, 'DYNAMIC APERTURE CALC n_turn > 1 FOR OPEN LATTICE')
          call err_exit       
        endif
          
        if (rf_is_on(tao_lat%lat%branch(0))) then
          ! pz surrounds the closed orbit
          scan%ref_orb = lat_branch%orb0
          scan%ref_orb%vec(6) = scan%ref_orb%vec(6) + u%dynamic_aperture%pz(j)
        else
          ! If the RF is off, new fixed points will be calculated for various pz
          orb(0)%vec(6) = u%dynamic_aperture%pz(j)
          call  closed_orbit_calc (tao_lat%lat, orb, 4)
          scan%ref_orb = orb(0)
        endif
    
        call out_io (s_info$, r_name, 'Starting Dynamic aperture scan for pz: \f10.6\ ... ', r_array = [u%dynamic_aperture%pz(j)])
        call run_timer ('START')
        call dynamic_aperture_scan(tao_lat%lat, scan, parallel = .true.)
        call run_timer ('READ', tt)
        call out_io (s_info$, r_name, 'Computation time for aperture scan at this energy (min): \f10.2\ ', tt/60)
      enddo
    endif
    
    ! PTC one-turn-map and normal form calc. Only for rings. 
    if (u%calc%one_turn_map .and. branch%param%geometry == closed$) then
      call set_ptc_verbose(.false.)
      rf_on = rf_is_on(branch)
      if (rf_on) then
        normal_form => branch%normal_form_with_rf
      else
        normal_form => branch%normal_form_no_rf
      endif
      
      if (.not. associated(normal_form%ele_origin)) normal_form%ele_origin => branch%ele(0)
      if (.not. associated(branch%ptc%m_t_layout)) call lat_to_ptc_layout (tao_lat%lat)
      
      ! Get one-turn-map
      call ptc_one_turn_map_at_ele (normal_form%ele_origin, normal_form%m, rf_on, pz = 0.0_rp )

      ! Get A, A_inv, dhdj
      call normal_form_taylors(normal_form%m, rf_on, dhdj = normal_form%dhdj, &
                                       A = normal_form%A, A_inverse = normal_form%A_inv)
      ! Get complex L and F
      call normal_form_complex_taylors(normal_form%m, rf_on, F = normal_form%F, L = normal_form%L)
    endif

    !

    !if (branch%param%geometry == closed$) then


    call tao_hook_branch_calc (u, tao_lat, branch)

  enddo branch_loop

  call tao_load_data_array (u, -1, 0, model$)

  ! If calc is on common model then transfer data to base of all other universes

  if (s%com%common_lattice .and. iuni == ix_common_uni$) then
    do j = lbound(s%u, 1), ubound(s%u, 1)
      if (j == ix_common_uni$) cycle
      s%u(j)%data(:)%base_value = u%data(:)%model_value
    enddo
  endif

  ! Mark this lattice as done 

  u%calc%lattice = .false.

enddo uni_loop

! do any post-processing

call tao_hook_post_process_data ()

call tao_set_var_useit_opt
call tao_set_data_useit_opt

end subroutine tao_lattice_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_single_track (u, tao_lat, calc_ok, ix_branch)

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct) u
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: lat_branch
type (ele_struct), pointer :: ele

real(rp), parameter :: vec0(6) = 0
integer i, ii, n, nn, ix_branch, status, ix_lost, i_dim

character(20) :: r_name = "tao_single_track"

logical calc_ok, err, radiation_fluctuations_on

!

lat => tao_lat%lat
branch => tao_lat%lat%branch(ix_branch)
lat_branch => tao_lat%lat_branch(ix_branch)
orbit => lat_branch%orbit

calc_ok = .true.
lat_branch%track_state = moving_forward$

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

    call closed_orbit_calc (lat, lat_branch%orbit, i_dim, 1, ix_branch, err_flag = err)
    if (err) then
      ! In desperation try a different starting point
      call init_coord (orbit(0), lat%beam_start, branch%ele(0), downstream_end$, orbit(0)%species)
      call closed_orbit_calc (lat, lat_branch%orbit, i_dim, 1, ix_branch, err_flag = err)
    endif
    if (err) then
      ! In desperation try a different starting point
      if (i_dim == 4) then
        orbit(0)%vec(1:5) = 0
      else
        orbit(0)%vec = 0
      endif
      call closed_orbit_calc (lat, lat_branch%orbit, i_dim, 1, ix_branch, err_flag = err)
    endif

    if (err) then
      calc_ok = .false.
      do i = 0, ubound(orbit, 1)
        orbit(i)%vec = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      enddo
    else
      lat_branch%orb0 = orbit(0)   ! Save beginning orbit as initial guess next time.
    endif

  else
    if (lat_branch%has_open_match_element) then
      do n = 1, branch%n_ele_track
        ele => branch%ele(n)
        call make_mat6(ele, branch%param, orbit(n-1), orbit(n), err_flag = err)
        if (err .or. .not. particle_is_moving_forward(orbit(n))) then
          lat_branch%track_state = n
          exit
        endif
        call twiss_propagate1(branch%ele(n-1), ele)
        branch%ele(n-1)%value(match_end_orbit$) = false$
      enddo

    else
      call track_all (lat, lat_branch%orbit, ix_branch, lat_branch%track_state, &
                                                            orbit0 = tao_lat%lat_branch(0)%orbit)
    endif

    if (lat_branch%track_state /= moving_forward$) then
      calc_ok = .false.
      ix_lost = lat_branch%track_state
      orbit(ix_lost+1:branch%n_ele_track) = coord_struct()
      call out_io (s_blank$, r_name, &
              "particle lost in single particle tracking at branch>>element \I0\>>\I0\: " // &
              trim(branch%ele(ix_lost)%name) // '  [s =\F9.2\]', &
              r_array = (/ branch%ele(ix_lost)%s /), i_array = (/ ix_branch, ix_lost /))
    endif
  endif

  bmad_com%radiation_fluctuations_on = radiation_fluctuations_on

endif

! Twiss

if (u%calc%mat6 .and. branch%param%particle /= photon$) then

  do i = 1, ix_lost - 1
    if (branch%ele(i)%tracking_method == linear$) then
      call lat_make_mat6 (lat, i, ix_branch = ix_branch)
    else
      call lat_make_mat6 (lat, i, orbit, ix_branch)
    endif
  enddo

  if (branch%param%geometry == closed$) then
    call twiss_at_start (lat, status, branch%ix_branch)
    if (status /= ok$) then
      calc_ok = .false.
      return
    endif
  endif

  call twiss_propagate_all (lat, ix_branch, err, 0, ix_lost - 1)

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
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (normal_modes_struct) :: modes
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: lat_branch
type (tao_element_struct), pointer :: uni_ele(:)
type (tao_universe_branch_struct), pointer :: uni_branch
type (bunch_params_struct), pointer :: bunch_params

integer ix_ele0, what_lat, n_lost_old
integer i, j, n, i_uni, ip, ig, ic, ie1, ie2
integer n_bunch, n_part, i_uni_to, ix_track
integer n_lost, ix_branch
integer, allocatable, save :: ix_ele(:)

character(20) :: r_name = "tao_beam_track"

real(rp) :: value1, value2, f, time, old_time

logical calc_ok, too_many_lost, print_err, err, lost

! Initialize 

call re_allocate (ix_ele, 1)

branch => tao_lat%lat%branch(ix_branch)
lat_branch => tao_lat%lat_branch(ix_branch)
uni_branch => u%uni_branch(ix_branch)
uni_ele => uni_branch%ele

beam_init => u%beam%beam_init

beam => u%beam%current

lat => tao_lat%lat

lat_branch%track_state = moving_forward$  ! Needed by tao_evaluate_a_datum
ix_track = moving_forward$
lost = .false.
calc_ok = .true.
too_many_lost = .false.

lat_branch%bunch_params(:)%n_particle_lost_in_ele = 0
lat_branch%bunch_params(:)%n_particle_live = 0

if (branch%param%geometry == closed$) then
  call out_io (s_fatal$, r_name, &
               'BEAM TRACKING WITH CIRCULAR LATTICE NOT YET IMPLEMENTED.', &
               'PLEASE SEE DCS IF YOU NEED THIS.')
  call err_exit
endif

! Transfer wakes from  design

do i = 1, branch%n_ele_max
  if (associated(branch%ele(i)%wake)) branch%ele(i)%wake%lr_mode = u%design%lat%branch(ix_branch)%ele(i)%wake%lr_mode
enddo

call zero_lr_wakes_in_lat (lat)

! Don't know what to do if the lattice is circular.

if (branch%param%geometry == closed$) then
  call out_io (s_abort$, r_name, 'BEAM TRACKING IN CIRCULAR LATTICE NOT YET IMPLEMENTED!')
  call err_exit
endif

! track through the lattice elements.
! The reference orbit for the transfer matrix calculation is taken to be the nominal 
! bunch center rather then the computed center to prevent jitter from changing things.

ie1 = ix_ele0
ie2 = branch%n_ele_track
if (uni_branch%ix_track_end > -1 .and. ix_branch == 0) ie2 = uni_branch%ix_track_end

print_err = .true.

if (s%global%beam_timer_on) then
  call run_timer ('START')
  old_time = 0
endif

beam = uni_branch%ele(ie1)%beam

n_lost_old = 0
do j = ie1, ie2

  bunch_params => lat_branch%bunch_params(j)
  ele => branch%ele(j)

  ! track to the element and save for phase space plot

  if (s%com%use_saved_beam_in_tracking) then
    beam = uni_ele(j)%beam

  else
    if (j /= ie1) then 
      call track_beam (lat, beam, branch%ele(j-1), ele, too_many_lost, centroid = lat_branch%orbit)
    endif

    if (uni_ele(j)%save_beam .or. ele%key == fork$ .or. ele%key == photon_fork$) uni_ele(j)%beam = beam
  endif
 
  ! Lost particles

  n_bunch = s%global%bunch_to_plot
  n_lost = count(beam%bunch(n_bunch)%particle(:)%state /= alive$)
  if (n_lost /= n_lost_old) then
    n = size(beam%bunch(n_bunch)%particle(:))
    if (size(s%u) == 1) then
      call out_io (s_blank$, r_name, &
            '\i0\ particle(s) lost at element ' // trim(ele_loc_to_string(ele)) // ': ' // trim(ele%name) // &
            '  Total lost: \i0\  of \i0\ ', &
            i_array = [n_lost-n_lost_old, n_lost, n])
    else
      call out_io (s_blank$, r_name, &
            '\i0\ particle(s) lost in universe \i0\ at element ' // trim(ele_loc_to_string(ele)) // &
                                                                           ': ' // trim(ele%name) // &
            '  Total lost: \i0\  of \i0\ ', &
            i_array = [n_lost-n_lost_old, u%ix_uni, n_lost, n])
    endif
    n_lost_old = n_lost
  endif

  if (tao_no_beam_left(beam, branch%param%particle) .and. .not. lost) then
    ix_track = j
    lost = .true.
    call out_io (s_warn$, r_name, &
            'TOO MANY PARTICLES HAVE BEEN LOST AT ELEMENT ' // trim(ele_loc_to_string(ele)) // &
                                          ': ' // trim(ele%name))
  endif

  ! calc bunch params

  call calc_bunch_params (u%beam%current%bunch(s%global%bunch_to_plot), &
                                                  bunch_params, err, print_err)
  if (err) print_err = .false.  ! Only generate one message.
  call tao_load_data_array (u, j, ix_branch, model$) 

  if (j == ie1) then
    bunch_params%n_particle_lost_in_ele = 0
  else
    bunch_params%n_particle_lost_in_ele = lat_branch%bunch_params(j-1)%n_particle_live - &
                                                     bunch_params%n_particle_live
  endif

  ! Timer 

  if (s%global%beam_timer_on) then
    call run_timer ('READ', time)
    if (time - old_time > 60) then
      call out_io (s_blank$, r_name, 'Beam at Element: \i0\. Time: \i0\ min', &
                          i_array = (/ j, nint(time/60) /) )
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
                                  i_array = (/u%ix_uni, n_lost /))

lat_branch%track_state = ix_track
 
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
type (spin_polar_struct) :: polar
type (branch_struct), pointer :: branch

integer ix_branch, i_ele_from, i_br_from

character(24) :: r_name = "tao_inject_particle"

! Not main branch case

branch => model%lat%branch(ix_branch)
i_br_from  = branch%ix_from_branch

if (i_br_from > -1) then
  i_ele_from = branch%ix_from_ele
  from_ele => model%lat%branch(i_br_from)%ele(i_ele_from)
  ele0 => branch%ele(0)
  call transfer_twiss (from_ele, ele0)
  orb_out => model%lat_branch(ix_branch)%orbit(0)
  call init_coord (orb_out, model%lat_branch(i_br_from)%orbit(i_ele_from), ele0, &
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

! In model%lat_branch()%orb0 is saved the last computed orbit. 
! This is important with common_lattice since tao_lat%lat_branch()%orbit(0) has been overwritten.

orb_in => model%lat_branch(ix_branch)%orb0
orb_out => model%lat_branch(ix_branch)%orbit(0)

call init_coord (orb_out, orb_in, branch%ele(0), downstream_end$, default_tracking_species(branch%param), 1, orb_in%p0c)

polar%theta = u%beam%beam_init%spin%theta
polar%phi = u%beam%beam_init%spin%phi
orb_out%spin = polar_to_vec (polar)

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

use tao_read_beam_mod
use beam_utils

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), target :: model
type (tao_universe_branch_struct), pointer :: uni_branch
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
ix_ele0 = 0

if (s%com%use_saved_beam_in_tracking) return

! if injecting into a branch then use the branch point as the starting distribution.

uni_branch => u%uni_branch(ix_branch)
branch => model%lat%branch(ix_branch)

if (ix_branch > 0) then
  ib = branch%ix_from_branch
  ie = branch%ix_from_ele

  if (.not. allocated (u%uni_branch(ib)%ele(ie)%beam%bunch)) then
    call out_io (s_error$, r_name, 'CANNOT INJECT INTO BRANCH FROM: ' // u%model%lat%branch(ib)%ele(ie)%name)
    return
  endif

  ix_ele0 = nint(model%lat%branch(ib)%ele(ie)%value(ix_to_element$))
  uni_branch%ele(ix_ele0)%beam = u%uni_branch(ib)%ele(ie)%beam
  init_ok = .true.
  return
endif

! Init for main branch...
! If there is an init file then read from the file

ix_ele0 = uni_branch%ix_track_start
beam_init => u%beam%beam_init
beam => uni_branch%ele(ix_ele0)%beam

if (u%beam%beam0_file /= "") then
  if (u%beam%init_beam0 .or. .not. allocated(beam%bunch)) then
    call tao_open_beam_file (u%beam%beam0_file, err)
    if (err) return
    call tao_set_beam_params (beam_init%n_bunch, beam_init%n_particle, beam_init%bunch_charge)
    call tao_read_beam (beam, err)
    if (err) return
    call tao_close_beam_file()
    n = 0
    do i = 1, size(beam%bunch)
      n = n + size(beam%bunch(i)%particle)
      do j = 1, size(beam%bunch(i)%particle)
        orbit => beam%bunch(i)%particle(j)
        orbit%vec = orbit%vec + model%lat%beam_start%vec
        if (orbit%state /= alive$) cycle  ! Don't want init_coord to raise the dead.
        call init_coord (orbit, orbit, branch%ele(ix_ele0), downstream_end$, branch%param%particle, +1, orbit%p0c, beam%bunch(i)%t_center)
      enddo
    enddo
    call out_io (s_info$, r_name, &
                  'Read initial beam distribution from: ' // u%beam%beam0_file, &
                  'Centroid Offset: \6es12.3\ ', &
                  'Number of particles: \i0\ ', &
                  r_array = model%lat%beam_start%vec, i_array = [n])
  endif

  if (tao_no_beam_left (beam, branch%param%particle)) then
    call out_io (s_warn$, r_name, "Not enough particles or no charge/intensity for beam init!")
    return
  endif

  init_ok = .true.
  return
endif

! Only reinit beam has not already been initialized or if commanded via %init_beam0.

if (u%beam%init_beam0 .or. .not. allocated(beam%bunch)) then
  if (beam_init%n_bunch < 1) beam_init%n_bunch = 1   ! Default if not set.
  call init_beam_distribution (branch%ele(ix_ele0), branch%param, beam_init, beam, err)
  if (err) then
    call out_io (s_fatal$, r_name, 'BEAM_INIT INITIAL BEAM PROPERTIES NOT SET FOR UNIVERSE: \i4\ ', u%ix_uni)
    return
  endif

  if (tao_no_beam_left(beam, branch%param%particle)) then
    call out_io (s_warn$, r_name, "Not enough particles or no charge/intensity for beam init!")
    return
  endif

  u%beam%init_beam0 = .false.
endif

init_ok = .true.

end subroutine tao_inject_beam
 
end module tao_lattice_calc_mod
