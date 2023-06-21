!+
! Subroutine tao_lattice_calc (calc_ok, print_err)
!
! For all universes track and calculate the lattice functions for the model lattice.
! Also compute model values for all data.
!
! Output:
!   calc_ok       -- logical: Set False if there was an error in the 
!                     calculation like a particle was lost or a lat is unstable.
!   print_err     -- logical, optional: Default True. If False, do not print error messages 
!                     if, for example, the lattice is unstable.
!-

subroutine tao_lattice_calc (calc_ok, print_err)

use tao_lattice_calc_mod, dummy => tao_lattice_calc
use srdt_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_dat
type (tao_d1_data_struct), pointer :: d1_dat
type (coord_struct), allocatable, target :: orb(:)
type (coord_struct), pointer :: orbit
type (tao_lattice_struct), pointer :: tao_lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_dynamic_aperture_struct), pointer :: da
type (beam_struct) beam
type (tao_plot_struct), pointer :: plot
type (tao_curve_struct), pointer :: curve
type (tao_curve_array_struct), allocatable :: crv(:), crv_temp(:)
type (lat_struct) lat_save

real(rp) pz0, dvec(6)
integer iuni, i, j, k, nc, ib, ix, iy, n_max, iu, it, id, n_turn

logical, optional :: print_err
logical calc_ok, this_calc_ok, err, mat_changed, track_this_beam

character(*), parameter :: r_name = "tao_lattice_calc"
character(20) name

! Lattice bookkeeping

calc_ok = .true.
s%com%is_err_message_printed = .false.  ! Reset for this round of computations
s%com%n_err_messages_printed = 0

do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)

  if (.not. u%is_on .or. .not. u%calc%lattice) cycle
  tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.
  call tao_lat_bookkeeper (u, tao_lat, err)
  if (err) then
    do id = 1, size(u%data)
      if (u%data(id)%data_type /= 'unstable.lattice') cycle
      if (.not. u%data(id)%exists) cycle
      call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
    enddo
    calc_ok = .false.
    return
  endif
enddo

! do a custom lattice calculation if desired

if (.not. s%global%lattice_calc_on) return

s%com%ix_ref_taylor = -999   ! Reset taylor map

call tao_hook_lattice_calc (calc_ok)

if (s%global%track_type /= 'single' .and. s%global%track_type /= 'beam') then
  call out_io (s_error$, r_name, 'UNKNOWN TRACK_TYPE: ' // quote(s%global%track_type), 'DEFAULTING TO "single"')
  s%global%track_type = 'single'
endif

if (s%global%track_type == 'beam' .and. count(s%u(:)%beam%track_beam_in_universe) == 0) then
  call out_io (s_error$, r_name, 'BEAM TRACKING CANNOT BE DONE UNLESS A TAO_BEAM_INIT NAMELIST HAS', &
                                 'BEEN DEFINED IN THE APPROPRIATE INIT FILE.')
endif

! To save time, s%u(:)%calc%lattice are used to determine what gets calculated. 

uni_loop: do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)
  if (.not. u%is_on) cycle

  if (u%calc%lattice) then
    ! Pointer to appropriate lattice and zero data array

    s%com%lattice_calc_done = .true.
    tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.
    u%spin_map%valid = .false.

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

      ! Even when beam tracking we need to calculate the lattice parameters with single tracking.

      call tao_inject_particle (u, tao_lat, ib)
      call tao_single_track (tao_lat, this_calc_ok, ib, print_err)
      if (.not. this_calc_ok) calc_ok = .false.

      ! Need to beam track even if single tracking is not OK since the merit function may depend
      ! upon the beam tracking.

      track_this_beam = (s%global%track_type == 'beam' .and. branch%param%particle /= photon$ .and. u%beam%track_beam_in_universe) 

      if (track_this_beam .or. s%global%init_lat_sigma_from_beam) call tao_inject_beam (u, tao_lat, ib, beam, this_calc_ok)

      if (track_this_beam) then
        if (.not. this_calc_ok) calc_ok = .false.
        if (this_calc_ok) then
          call tao_beam_track (u, tao_lat, ib, beam, this_calc_ok)
          if (.not. this_calc_ok) calc_ok = .false.
        else
          tao_branch%bunch_params(:)%n_particle_lost_in_ele = 0
          tao_branch%bunch_params(:)%n_particle_live = 0
        endif
      endif

      if (this_calc_ok) call tao_lat_sigma_track (tao_lat, this_calc_ok, ib, print_err)

      !

      if (s%com%force_chrom_calc .or. u%calc%chrom_for_data .or. u%calc%chrom_for_plotting) then
        call chrom_calc (tao_lat%lat, s%global%delta_e_chrom, tao_branch%a%chrom, tao_branch%b%chrom, err, &
                tao_branch%orbit(0)%vec(6), low_E_lat=tao_lat%low_E_lat, high_E_lat=tao_lat%high_E_lat, &
                ix_branch = ib, orb0 = tao_branch%orbit(0))
      endif

      ! do multi-turn tracking if needed.

      nc = 0
      n_turn = 0
      do i = 1, size(s%plot_page%region)
        if (.not. s%plot_page%region(i)%visible) cycle
        plot => s%plot_page%region(i)%plot
        if (.not. allocated(plot%graph)) cycle
        do j = 1, size(plot%graph)
          if (plot%graph(j)%type /= 'phase_space') cycle
          if (.not. allocated(plot%graph(j)%curve)) cycle
          do k = 1, size(plot%graph(j)%curve)
            curve => plot%graph(j)%curve(k)
            if (curve%data_source /= 'multi_turn_orbit' .and. curve%data_source /= 'rel_multi_turn_orbit') cycle
            if ((curve%ix_universe == -1 .and. iuni /= s%global%default_universe) .or. &
                                  (curve%ix_universe /= -1 .and. iuni /= curve%ix_universe)) cycle
            if (tao_branch_index(curve%ix_branch) /= ib) cycle
            nc = nc + 1
            s%com%multi_turn_orbit_is_plotted = .true.
            if (nc == 1) then
              if (allocated(crv)) deallocate(crv)
              allocate(crv(1))
            else
              call move_alloc(crv, crv_temp)
              allocate(crv(nc))
              crv(1:nc-1) = crv_temp
              deallocate(crv_temp)
            endif

            if (curve%n_turn < 1) then
              call out_io (s_warn$, r_name, 'Multi-turn orbit curve has %n_turn component non-positive! ' // &
                                                tao_curve_name(curve))
              curve%why_invalid = '%n_turn component is not positive.'
              curve%valid = .false.
              cycle
            endif

            crv(nc)%c => curve
            n_turn = max(n_turn, curve%n_turn)
            curve%valid = .true.
            call re_allocate(curve%x_symb, curve%n_turn)
            call re_allocate(curve%y_symb, curve%n_turn)
            call re_allocate(curve%ix_symb, curve%n_turn)
          enddo
        enddo
      enddo

      if (nc > 0) then
        call reallocate_coord (orb, branch%n_ele_max)
        orb(0) = tao_lat%lat%particle_start
        pz0 = orb(0)%vec(6)
        dvec = 0
        do it = 1, n_turn
          call track_all (tao_lat%lat, orb, ib)
          do k = 1, nc
            curve => crv(k)%c
            if (.not. curve%valid) cycle
            if (it > curve%n_turn) cycle

            orbit => orb(curve%ix_ele_ref)
            if (orbit%state /= alive$) then
              call re_allocate(curve%x_symb, it-1)
              call re_allocate(curve%y_symb, it-1)
              call re_allocate(curve%ix_symb, it-1)
              cycle
            endif

            call match_word (curve%data_type_x, coord_name, ix)
            call match_word (curve%data_type, coord_name, iy)
            if (ix == 0) then
              call out_io(s_error$, r_name, 'MALFORMED CURVE%DATA_TYPE_X: ' // quote(curve%data_type_x), &
                                            'FOR: ' //tao_curve_name(curve))
              curve%why_invalid = 'MALFORMED CURVE%DATA_TYPE_X'
              curve%valid = .false.
              cycle
            endif
            if (iy == 0) then
              call out_io(s_error$, r_name, 'MALFORMED CURVE%DATA_TYPE: ' // quote(curve%data_type), &
                                            'FOR: ' //tao_curve_name(curve))
              curve%why_invalid = 'MALFORMED CURVE%DATA_TYPE'
              curve%valid = .false.
              cycle
            endif
            if (curve%data_source == 'rel_multi_turn_orbit') then
              ele => branch%ele(curve%ix_ele_ref)
              dvec = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, 0.0_rp, 0.0_rp] * pz0
            endif
            curve%x_symb(it) = (orbit%vec(ix) - dvec(ix)) * curve%g%x_axis_scale_factor
            curve%y_symb(it) = (orbit%vec(iy) - dvec(iy)) * curve%y_axis_scale_factor
            curve%ix_symb(it) = it
          enddo
          orb(0) = orb(branch%n_ele_track)
          if (orb(0)%state /= alive$) exit
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

    ! Dynamic aperture calc. Only for rings

    if (u%calc%dynamic_aperture) then
      lat_save = tao_lat%lat  ! Save since DA calc varies pz
      da => u%dynamic_aperture
      call dynamic_aperture_scan(da%scan, da%param, da%pz, tao_lat%lat, print_timing = (.not. s%com%optimizer_running))
      tao_lat%lat = lat_save  ! Save since DA calc varies pz
    endif
  endif  ! if (u%calc%lattice)

  ! Calculate non-expression data 

  do id = 1, size(u%data)
    if (substr(u%data(id)%data_type,1,11) == 'expression:') cycle
    if (.not. u%data(id)%exists) cycle
    call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
  enddo

  ! Mark this lattice as done 

  u%calc%lattice = .false.

  call tao_scale_ping_data(u)

enddo uni_loop

s%com%force_chrom_calc   = .false.
s%com%force_rad_int_calc = .false.

! do any post-processing

do iuni = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(iuni)
  if (.not. u%is_on) cycle

  do id = 1, size(u%data)
    if (substr(u%data(id)%data_type,1,11) /= 'expression:') cycle
    if (.not. u%data(id)%exists) cycle
    call tao_evaluate_a_datum (u%data(id), u, u%model, u%data(id)%model_value, u%data(id)%good_model)
  enddo
enddo

call tao_hook_post_process_data ()

call tao_set_var_useit_opt
call tao_set_data_useit_opt

end subroutine tao_lattice_calc
