!+
! Module: tao_lattice_calc_mod
!
! Lattice calculation routines are here. It's a module so that custom lattice
! calculations has access to the subroutines.

module tao_lattice_calc_mod

use tao_mod
use tao_data_mod
use tao_calc_params_mod
use macroparticle_mod
use beam_mod

integer, parameter :: design$ = 1
integer, parameter :: model$ = 2

logical all_lost_already

contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lattice_calc (init_design)
!
! Routine to calculate the lattice functions and TAO data. Always tracks through the model
! lattice. If initing the design lattices then the tracking will still be done
! through the model lattices. tao_init then transfers this into the
! design lattices.
! 
! Input:
!  init_design      -- Logical: (optional) To initialize design lattices
!
! Output:
!-

subroutine tao_lattice_calc (init_design)

implicit none

logical, optional :: init_design

type (ring_struct), pointer :: lattice
type (coord_struct), pointer :: orbit(:)

integer i, j

character(20) :: r_name = "tao_lattice_calc"

logical initing_design
! used refers to if a custom lattice calculation is performed
logical, automatic :: used(size(s%u))

  used(:) = .false.

  ! have to do this because the logical init_design is optional
  initing_design = .false.
  if (present(init_design)) then
    if (init_design) initing_design = .true.
  endif

  ! make sure useit is up-to-date
  if (.not. initing_design) then
    call tao_set_var_useit_opt
    call tao_set_data_useit_opt
  endif
  
  ! do a custom lattice calculation if desired
  if (s%global%lattice_recalc) then
    do i = 1, size(s%u)
      lattice => s%u(i)%model
      orbit => s%u(i)%model_orb
      call tao_lat_bookkeeper (s%u(i), lattice)
      call tao_hook_lattice_calc (s%u(i), lattice, orbit, used(i))
    enddo
    if (.not. any(.not. used)) s%global%lattice_recalc = .false.
  endif
    
  ! Closed orbit and Twiss calculation.
  ! This can be slow for large lattices so only do it if the lattice changed.
  if (s%global%lattice_recalc) then
    if (s%global%track_type == 'single') then
      do i = 1, size(s%u)
        if (.not. s%u(i)%is_on .or. used(i)) cycle
        lattice => s%u(i)%model
        orbit => s%u(i)%model_orb
        ! set up matching element
        if (initing_design) call tao_match_lats_init (s%u(i))
        call tao_inject_particle (s%u(i), lattice, orbit)
        call tao_lat_bookkeeper (s%u(i), lattice)
	call tao_single_track (i, lattice, orbit)
      enddo
    elseif (s%global%track_type == 'beam') then
      do i = 1, size(s%u)
        if (.not. s%u(i)%is_on .or. used(i)) cycle
        lattice => s%u(i)%model
        orbit => s%u(i)%model_orb
        ! set up matching element
        if (initing_design) call tao_match_lats_init (s%u(i))
        call tao_inject_beam (s%u(i), lattice, orbit)
        call tao_lat_bookkeeper (s%u(i), lattice)
        call tao_beam_track (i, lattice, orbit)
      enddo
    elseif (s%global%track_type == 'macro') then
      do i = 1, size(s%u)
        if (.not. s%u(i)%is_on .or. used(i)) cycle
        lattice => s%u(i)%model
        orbit => s%u(i)%model_orb
        ! set up matching element
        if (initing_design) call tao_match_lats_init (s%u(i))
        call tao_inject_macro_beam (s%u(i), lattice, orbit)
        call tao_lat_bookkeeper (s%u(i), lattice)
        call tao_macro_track (i, lattice, orbit)
      enddo
    else
      call out_io (s_fatal$, r_name, &
                     "This tracking type has yet to be implemented!")
      call out_io (s_blank$, r_name, &
                     "No tracking or twiss calculations will be perfomred.")
    endif
    ! do any post-processing
    call tao_hook_post_process_data ()
    s%global%lattice_recalc = .false.
  endif

end subroutine tao_lattice_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_single_track (uni, lat, orb)

implicit none

type (ring_struct) lat
type (coord_struct) :: orb(0:)

integer uni, i, ii

character(20) :: r_name = "tao_single_track"

  if (lat%param%lattice_type == circular_lattice$) then
    call ring_make_mat6 (lat, -1)
    call twiss_at_start (lat)
    call closed_orbit_at_start (lat, orb(0), 4, .true.)
  endif

  lat%param%ix_lost = not_lost$

  i = 0
  if (lat%param%lattice_type == circular_lattice$) call twiss_at_start (lat)
  call tao_load_data_array (s%u(uni), i)

  do i = 1, lat%n_ele_use
    call track1 (orb(i-1), lat%ele_(i), lat%param, orb(i))

    if (lat%param%lost) then
      lat%param%ix_lost = i
      do ii = i+1, lat%n_ele_use
        orb(ii)%vec = 0
      enddo
      return
    endif

    call tao_calc_params (s%u(uni), i)
    call tao_load_data_array (s%u(uni), i)
    
  enddo
  
  if (lat%param%lost) &
    call out_io (s_blank$, r_name, "particle lost at element \I\.", &
                                            lat%param%ix_lost)

end subroutine tao_single_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Right now, there is no beam tracking in rings. If extracting from a
! ring then the beam is generated at the extraction point.

subroutine tao_beam_track (uni, lat, orb)

implicit none

integer uni, what_lat
type (ring_struct) :: lat
type (coord_struct) :: orb(0:)
type (tao_universe_struct), pointer :: u
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (modes_struct) :: modes

integer j, i_uni
integer n_bunch, n_part
integer extract_at_ix_ele, n_lost

character(20) :: r_name = "tao_beam_track"

real(rp) :: value1, value2, m_particle

logical post

  u => s%u(uni)
  beam => u%beam%beam
  beam_init => u%beam%beam_init

  all_lost_already = .false.

  ! Find if injecting into another lattice
  extract_at_ix_ele = -1
  inject_loop: do i_uni = uni+1, size(s%u)
    if (s%u(i_uni)%coupling%coupled) then
      if (s%u(i_uni)%coupling%from_uni == uni) then
	if (s%u(i_uni)%coupling%from_uni_ix_ele .ne. -1) then
	  extract_at_ix_ele = s%u(i_uni)%coupling%from_uni_ix_ele
	  exit inject_loop ! save i_uni for coupled universe
	else
	  call out_io (s_abort$, r_name, &
	       "Must specify an element when coupling lattices with a beam.")
	  call err_exit
	endif
      endif
    endif
  enddo inject_loop
  
  ! If beam is injected into this lattice then no initialization wanted.
  if (.not. u%coupling%coupled) then
    ! no beam tracking in rings
    if (lat%param%lattice_type == circular_lattice$ ) then
      call tao_single_track (uni, lat, orb) 
      if (u%beam%calc_emittance) then
        call radiation_integrals (lat, orb, modes)
        if (extract_at_ix_ele .ne. -1) then
	  if (lat%param%particle == electron$ .or. &
	      lat%param%particle == positron$       ) then
	    m_particle = m_electron
	  elseif (lat%param%particle == proton$ .or. &
	          lat%param%particle == antiproton$   ) then
	    m_particle = m_proton
	  endif
          beam_init%a_norm_emitt  = modes%a%emittance * &
	                  ((lat%ele_(extract_at_ix_ele)%value(beam_energy$)*&
			         (1+orb(extract_at_ix_ele)%vec(6))) / m_particle)
          beam_init%b_norm_emitt  = modes%b%emittance * &
	                  ((lat%ele_(extract_at_ix_ele)%value(beam_energy$)*&
			         (1+orb(extract_at_ix_ele)%vec(6))) / m_particle)
        endif
      endif
      !transfer extracted particle info into macro_init
      if (extract_at_ix_ele .ne. -1) then
        beam_init%center  = orb(extract_at_ix_ele)%vec
        ! other beam_init parameters will be as in init.tao, or as above
        call init_beam_distribution (lat%ele_(extract_at_ix_ele), &
                                     beam_init, s%u(i_uni)%coupling%injecting_beam, &
	                             .true., .true.)
      endif
      return
    elseif (lat%param%lattice_type == linear_lattice$) then
      ! set initial beam centroid
      beam_init%center = orb(0)%vec
    else
      call out_io (s_error$, r_name, &
                   "This lattice type not yet implemented for beam tracking!")
      call err_exit
    endif
    
    call init_beam_distribution (lat%ele_(0), beam_init, beam, .true., .true.)
  endif

  ! beginning element calculations
  call tao_load_data_array (u, 0) 

  ! track through every element
  do j = 1, lat%n_ele_use
    ! track to the element
    call track_beam (lat, beam, j-1, j)
 
    ! Save beam at location if injecting into another lattice
    if (extract_at_ix_ele == j) then
      call beam_equal_beam (s%u(i_uni)%coupling%injecting_beam, beam)
    endif
    
    ! compute centroid orbit
    call tao_find_beam_centroid (beam, orb(j), uni, j)

    ! Find lattice and beam parameters
    call tao_calc_params (u, j)
    
    ! load data
    call tao_load_data_array (u, j) 
  enddo

  ! only post total lost if no extraction or extracting to a turned off lattice
  post = .false.
  if (extract_at_ix_ele == -1) post = .true.
  if (uni .lt. size(s%u)) then
    if (.not. s%u(uni+1)%is_on) post = .true.
  endif
  if (post) then
    n_lost = 0 
    do n_bunch = 1, size(beam%bunch)
      n_lost = n_lost + count(beam%bunch(n_bunch)%particle%ix_lost /= not_lost$)
    enddo
    call out_io (s_blank$, r_name, &
      "Total number of lost particles by the end of universe \I2\: \I5\.", &
		                              i_array = (/uni, n_lost /))
  endif
  
end subroutine tao_beam_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Right now, there is no macroparticle tracking in rings. If extracting from a
! ring then the macroparticle beam is generated at the extraction point.

subroutine tao_macro_track (uni, lat, orb)

use macro_utils_mod

implicit none

integer uni, what_lat
type (ring_struct) :: lat
type (coord_struct) :: orb(0:)
type (tao_universe_struct), pointer :: u
type (macro_beam_struct), pointer :: beam
type (macro_init_struct), pointer :: macro_init
type (modes_struct) :: modes

integer j, i_uni
integer n_bunch, n_slice, n_macro
integer extract_at_ix_ele, n_lost

character(20) :: r_name = "macro_track"

real(rp) :: value1, value2, m_particle

logical post

  u => s%u(uni)
  beam => u%macro_beam%beam
  macro_init => u%macro_beam%macro_init

  all_lost_already = .false.

  ! Find if injecting into another lattice
  extract_at_ix_ele = -1
  inject_loop: do i_uni = uni+1, size(s%u)
    if (s%u(i_uni)%coupling%coupled) then
      if (s%u(i_uni)%coupling%from_uni == uni) then
	if (s%u(i_uni)%coupling%from_uni_ix_ele .ne. -1) then
	  extract_at_ix_ele = s%u(i_uni)%coupling%from_uni_ix_ele
	  exit inject_loop ! save i_uni for coupled universe
	else
	  call out_io (s_abort$, r_name, &
	       "Must specify an element when coupling lattices with macroparticles")
	  call err_exit
	endif
      endif
    endif
  enddo inject_loop
  
  ! If beam is injected to here then no initialization wanted.
  if (.not. u%coupling%coupled) then
    ! no macroparticle tracking in rings
    if (lat%param%lattice_type == circular_lattice$ ) then
      call tao_single_track (uni, lat, orb) 
      if (u%macro_beam%calc_emittance) then
        call radiation_integrals (lat, orb, modes)
        if (extract_at_ix_ele .ne. -1) then
	  if (lat%param%particle == electron$ .or. &
	      lat%param%particle == positron$       ) then
	    m_particle = m_electron
	  elseif (lat%param%particle == proton$ .or. &
	          lat%param%particle == antiproton$   ) then
	    m_particle = m_proton
	  endif
          macro_init%x%norm_emit  = modes%a%emittance * &
	                  ((lat%ele_(extract_at_ix_ele)%value(beam_energy$)*&
			         (1+orb(extract_at_ix_ele)%vec(6))) / m_particle)
          macro_init%y%norm_emit  = modes%b%emittance * &
	                  ((lat%ele_(extract_at_ix_ele)%value(beam_energy$)*&
			         (1+orb(extract_at_ix_ele)%vec(6))) / m_particle)
        endif
      endif
      !transfer extracted particle info into macro_init
      if (extract_at_ix_ele .ne. -1) then
        macro_init%center  = orb(extract_at_ix_ele)%vec
        ! other macro_init parameters will be as in init.tao, or as above
        call init_macro_distribution (s%u(i_uni)%coupling%injecting_macro_beam, &
	                              macro_init, lat%ele_(extract_at_ix_ele), .true.)
      endif
      return
    elseif (lat%param%lattice_type == linear_lattice$) then
      ! set initial beam centroid
      macro_init%center = orb(0)%vec
    else
      call out_io (s_error$, r_name, &
                   "This lattice type not yet implemented for macroparticles!")
      call err_exit
    endif
    
    call init_macro_distribution (beam, macro_init, lat%ele_(0), .true.)
    
    u%macro_beam%ix_lost(:,:,:) = not_lost$
  endif

  ! beginning element calculations
  call tao_load_data_array (u, 0) 

  ! track through every element
  do j = 1, lat%n_ele_use
    ! track to the element
    call track_macro_beam (lat, beam, j-1, j)
 
    ! Save beam at location if injecting into another lattice
    if (extract_at_ix_ele == j) then
      s%u(i_uni)%coupling%injecting_macro_beam = beam
    endif
    
    ! compute centroid orbit
    call tao_find_macro_beam_centroid (beam, orb(j), uni, j, u%macro_beam)

    ! Find lattice and beam parameters
    call tao_calc_params (u, j)
    
    ! load data
    call tao_load_data_array (u, j) 
  enddo

  ! only post total lost if no extraction or extracting to a turned off lattice
  post = .false.
  if (extract_at_ix_ele == -1) post = .true.
  if (uni .lt. size(s%u)) then
    if (.not. s%u(uni+1)%is_on) post = .true.
  endif
  if (post) then
    n_lost = 0 
    do n_bunch = 1, size(beam%bunch)
      do n_slice = 1, size(beam%bunch(n_bunch)%slice)
        do n_macro = 1, size(beam%bunch(n_bunch)%slice(n_slice)%macro)
          if (beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)%lost) &
	    n_lost = n_lost + 1
	enddo
      enddo
    enddo
    call out_io (s_blank$, r_name, &
      "Total number of lost macroparticles by the end of universe \I2\: \I5\.", &
		                              i_array = (/uni, n_lost /))
  endif
  
end subroutine tao_macro_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Find the centroid of all particles in viewed bunches
! Also keep track of lost macroparticles if the optional arguments are passed

subroutine tao_find_beam_centroid (beam, orb, uni, ix_ele)

implicit none

type (beam_struct), target :: beam
type (coord_struct) :: orb, coord
integer, optional :: uni, ix_ele

integer n_bunch, n_part, n_lost, tot_part, i_ele

logical record_lost

character(20) :: r_name = "tao_find_beam_centroid"

  coord%vec = 0.0
  tot_part = 0
  n_lost = 0

  n_bunch = s%global%bunch_to_plot
  
  ! If optional arguments not present no verbose and
  !  just check if particles lost
  if (present(uni) .or. present(ix_ele)) then
    record_lost = .true.
    i_ele = ix_ele
  else
    record_lost = .false.
    i_ele = 0
  endif

  
  do n_part = 1, size(beam%bunch(n_bunch)%particle)
    ! only average over particles that haven't been lost
    if (record_lost .and. beam%bunch(n_bunch)%particle(n_part)%ix_lost == i_ele) then
      n_lost = n_lost + 1
      cycle
    elseif (beam%bunch(n_bunch)%particle(n_part)%ix_lost .ne. not_lost$) then
      cycle
    endif
    tot_part = tot_part + 1
    coord%vec = coord%vec + beam%bunch(n_bunch)%particle(n_part)%r%vec
  enddo
      
  ! Post lost particles
  if (record_lost) then
    if (n_lost == 1) then
      call out_io (s_blank$, r_name, &
            "   1 particle  lost at element \I\ in universe \I\.", &
                                           i_array = (/ ix_ele, uni /) )
    elseif (n_lost .gt. 1) then
      call out_io (s_blank$, r_name, &
            "\I4\ particles lost at element \I\ in universe \I\.", &
                                           i_array = (/ n_lost, ix_ele, uni /) )
    endif
  endif
  
  ! average
  if (tot_part .ne. 0) then
    orb%vec = coord%vec / tot_part
  else 
    ! lost all particles
    if (record_lost .and. .not. all_lost_already) &
      call out_io (s_warn$, r_name, "All particles have been lost!!!!!!!!!!!!!")
      all_lost_already = .true.
    orb%vec = 0.0
  endif
 
end subroutine tao_find_beam_centroid

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Find the centroid of all macroparticles in all slices in all bunches
! Also keep track of lost macroparticles if the optional arguments are passed

subroutine tao_find_macro_beam_centroid (beam, orb, uni, ix_ele, u_beam)

implicit none

type (macro_beam_struct), target :: beam
type (coord_struct) :: orb, coord
integer, optional :: uni, ix_ele
type (tao_macro_beam_struct), optional :: u_beam
type (macro_struct), pointer :: macro

real charge, tot_charge

integer n_bunch, n_slice, n_macro

logical record_lost

character(20) :: r_name = "tao_find_macro_beam_centroid"

  ! are we checking keep track of lost particles?
  if (present(ix_ele) .and. present(u_beam) .and. present(uni)) then
    record_lost = .true.
  else
    record_lost = .false.
  endif

  tot_charge = 0.0
  coord%vec = 0.0

  n_bunch = s%global%bunch_to_plot
  slice_loop: do n_slice = 1, size(beam%bunch(n_bunch)%slice)
    macro_loop: do n_macro = 1, size(beam%bunch(n_bunch)%slice(n_slice)%macro)
      macro => beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)
      ! only average over particles that haven't been lost
      if (record_lost) then
        if (u_beam%ix_lost(n_bunch,n_slice,n_macro) /= not_lost$) cycle macro_loop
      else
        if (macro%lost) cycle macro_loop
      endif
      ! check for lost macro through this element
      if (macro%lost .and. record_lost) then
        u_beam%ix_lost(n_bunch,n_slice,n_macro) = ix_ele
        call out_io (s_blank$, r_name, &
              "Macroparticle lost at element \I\ in universe \I\.", &
                                           i_array = (/ ix_ele, uni /) )
        cycle macro_loop
      endif
      charge = macro%charge
      tot_charge = tot_charge + charge
      coord%vec = coord%vec + charge * macro%r%vec
    enddo macro_loop
  enddo slice_loop
      
  ! average

  if (tot_charge .ne. 0) then
    orb%vec = coord%vec / tot_charge
  else 
    ! lost all macros
    if (.not. all_lost_already) then
      call out_io (s_warn$, r_name, &
            "All macroparticles have been lost!!!!!!!!!!!!!")
      all_lost_already = .true.
      orb%vec = 0.0
    endif
  endif
 
end subroutine tao_find_macro_beam_centroid

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a particle from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated. If no injection then will set beginning orbit to
! whatever the user has specified. As always, tracking only occure in the model
! lattice.

subroutine tao_inject_particle (u, lat, orbit)

implicit none

type (tao_universe_struct) u
type (ring_struct) lat
type (coord_struct) orbit(0:)

type (ele_struct), save :: extract_ele
type (coord_struct) pos

character(20) :: r_name = "inject_particle"

  if (.not. u%coupling%coupled) return
    
  if (.not. s%u(u%coupling%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed")
    return
  endif
  
  call init_ele (extract_ele)
  
  ! get particle perameters from previous universe at position s
  call twiss_and_track_at_s (s%u(u%coupling%from_uni)%model, &
                             u%coupling%from_uni_s, extract_ele, &
		             s%u(u%coupling%from_uni)%model_orb, pos)

  !track through coupling element
  if (u%coupling%use_coupling_ele) then
    call make_mat6 (u%coupling%coupling_ele, s%u(u%coupling%from_uni)%design%param)
    call twiss_propagate1 (extract_ele, u%coupling%coupling_ele)
    call track1 (pos, u%coupling%coupling_ele, &
                    s%u(u%coupling%from_uni)%design%param, pos)
    u%coupling%coupling_ele%value(beam_energy$) = extract_ele%value(beam_energy$)
    u%coupling%coupling_ele%floor = extract_ele%floor
    extract_ele = u%coupling%coupling_ele
  endif
  
  ! transfer to this lattice
  lat%ele_(0)%x = extract_ele%x
  lat%ele_(0)%y = extract_ele%y
  lat%ele_(0)%z = extract_ele%z
  lat%ele_(0)%value(beam_energy$) = extract_ele%value(beam_energy$)
  lat%ele_(0)%c_mat   = extract_ele%c_mat
  lat%ele_(0)%gamma_c = extract_ele%gamma_c
  lat%ele_(0)%floor   = extract_ele%floor
  orbit(0)      = pos
		    
end subroutine tao_inject_particle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a beam from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated.

subroutine tao_inject_beam (u, lat, orbit)

implicit none
 
type (tao_universe_struct) u
type (ring_struct) lat
type (coord_struct) orbit(0:)
type (ele_struct), save :: extract_ele

type (param_struct), pointer :: param

integer n_bunch, n_part

character(20) :: r_name = "tao_inject_beam"

  if (.not. u%coupling%coupled) return
   
  if (.not. s%u(u%coupling%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed.")
    return
  endif
  
  ! beam from previous universe at end of extracting element should already be set
  ! but we still need the twiss parameters and everything else
  extract_ele = s%u(u%coupling%from_uni)%model%ele_(u%coupling%from_uni_ix_ele)

  !track through coupling element
  if (u%coupling%use_coupling_ele) then
    param => s%u(u%coupling%from_uni)%model%param
    call make_mat6 (u%coupling%coupling_ele, param)
    call twiss_propagate1 (extract_ele, u%coupling%coupling_ele)
    call track1_beam (u%coupling%injecting_beam, u%coupling%coupling_ele, &
                      param, u%coupling%injecting_beam)
    u%coupling%coupling_ele%value(beam_energy$) = extract_ele%value(beam_energy$)
    u%coupling%coupling_ele%floor = extract_ele%floor
    extract_ele = u%coupling%coupling_ele
  endif
  
  ! transfer to this lattice
  lat%ele_(0)%x = extract_ele%x
  lat%ele_(0)%y = extract_ele%y
  lat%ele_(0)%z = extract_ele%z
  lat%ele_(0)%value(beam_energy$) = extract_ele%value(beam_energy$)
  lat%ele_(0)%c_mat   = extract_ele%c_mat
  lat%ele_(0)%gamma_c = extract_ele%gamma_c
  lat%ele_(0)%floor   = extract_ele%floor
  call beam_equal_beam (u%beam%beam, u%coupling%injecting_beam)
  call tao_find_beam_centroid (u%coupling%injecting_beam, orbit(0))

end subroutine tao_inject_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a beam from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated.

subroutine tao_inject_macro_beam (u, lat, orbit)

implicit none
 
type (tao_universe_struct) u
type (ring_struct) lat
type (coord_struct) orbit(0:)
type (ele_struct), save :: extract_ele

type (param_struct), pointer :: param

integer n_bunch, n_slice, n_macro

character(20) :: r_name = "tao_inject_macro_beam"

  if (.not. u%coupling%coupled) return
   
  if (.not. s%u(u%coupling%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed")
    return
  endif
  
  ! beam from previous universe at end of element should already be set
  ! but we still need the twiss parameters and everything else
  extract_ele = s%u(u%coupling%from_uni)%model%ele_(u%coupling%from_uni_ix_ele)

  !track through coupling element
  if (u%coupling%use_coupling_ele) then
    param => s%u(u%coupling%from_uni)%model%param
    call make_mat6 (u%coupling%coupling_ele, param)
    call twiss_propagate1 (extract_ele, u%coupling%coupling_ele)
    call track1_macro_beam (u%coupling%injecting_macro_beam, u%coupling%coupling_ele, &
                                            param, u%coupling%injecting_macro_beam)
    u%coupling%coupling_ele%value(beam_energy$) = extract_ele%value(beam_energy$)
    u%coupling%coupling_ele%floor = extract_ele%floor
    extract_ele = u%coupling%coupling_ele
  endif
  
  ! transfer to this lattice
  lat%ele_(0)%x = extract_ele%x
  lat%ele_(0)%y = extract_ele%y
  lat%ele_(0)%z = extract_ele%z
  lat%ele_(0)%value(beam_energy$) = extract_ele%value(beam_energy$)
  lat%ele_(0)%c_mat   = extract_ele%c_mat
  lat%ele_(0)%gamma_c = extract_ele%gamma_c
  lat%ele_(0)%floor   = extract_ele%floor
  u%macro_beam%beam = u%coupling%injecting_macro_beam
  call tao_find_macro_beam_centroid (u%coupling%injecting_macro_beam, orbit(0))

  ! deterine if macroparticle already lost
  do n_bunch = 1, size(u%coupling%injecting_macro_beam%bunch)
    do n_slice = 1, size(u%coupling%injecting_macro_beam%bunch(n_bunch)%slice)
      do n_macro = 1, size(u%coupling%injecting_macro_beam%bunch(n_bunch)%slice(n_slice)%macro)
	if (u%coupling%injecting_macro_beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)%lost) &
	  u%macro_beam%ix_lost(n_bunch, n_slice, n_macro) = 0
      enddo
    enddo
  enddo
  
end subroutine tao_inject_macro_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will set up the matching element if needed. The lattice where the
! extraction occurs will have already been calculated but the injected lattice
! will not have been calculated yet.

subroutine tao_match_lats_init (u)

implicit none

type (tao_universe_struct), target :: u
type (coord_struct) extract_pos
type (ele_struct), save :: extract_ele, inject_ele
type (ele_struct), pointer :: coupling_ele
type (beam_struct), pointer :: injecting_beam
type (macro_beam_struct), pointer :: injecting_macro_beam

character(20) :: r_name = "match_lats_init"

  if (.not. (u%coupling%coupled .and. u%coupling%use_coupling_ele)) return

  coupling_ele => u%coupling%coupling_ele
  call init_ele (extract_ele)
  call init_ele (inject_ele)
  
  ! set up coupling%injecting_beam or macro_beam 
  if (s%global%track_type == 'beam') then
    injecting_beam => s%u(u%coupling%from_uni)%beam%beam
    call reallocate_beam (u%coupling%injecting_beam, size(injecting_beam%bunch), &
            size(injecting_beam%bunch(1)%particle))
  elseif (s%global%track_type == 'macro') then
    injecting_macro_beam => s%u(u%coupling%from_uni)%macro_beam%beam
    call reallocate_macro_beam (u%coupling%injecting_macro_beam, &
                                size(injecting_macro_beam%bunch), &
                                size(injecting_macro_beam%bunch(1)%slice), &
                                size(injecting_macro_beam%bunch(1)%slice(1)%macro))
  endif

  ! match design lattices
  if (u%coupling%match_to_design) then
    ! get twiss parameters from extracted lattice
    if (s%global%track_type == 'single') then
      call twiss_and_track_at_s (s%u(u%coupling%from_uni)%design, &
                     u%coupling%from_uni_s, extract_ele, &
          	   s%u(u%coupling%from_uni)%design_orb, extract_pos)
    elseif (s%global%track_type == 'beam') then
      extract_ele = s%u(u%coupling%from_uni)%design%ele_(u%coupling%from_uni_ix_ele)
    elseif (s%global%track_type == 'macro') then
      extract_ele = s%u(u%coupling%from_uni)%design%ele_(u%coupling%from_uni_ix_ele)
    endif
    
    ! get twiss parameters for injected lattice
    ! This is performed before the standard lattice calculation so the design
    ! twiss parameters in ele_(0) will still be as set in the lattice file.
    inject_ele = u%design%ele_(0)
  else
    call out_io (s_warn$, r_name, &
                "The coupling element will only match the design lattices")
    return
  endif

  ! set up matching element
  coupling_ele%value( beta_x0$) = extract_ele%x%beta
  coupling_ele%value(alpha_x0$) = extract_ele%x%alpha
  coupling_ele%value(  eta_x0$) = extract_ele%x%eta
  coupling_ele%value( etap_x0$) = extract_ele%x%etap
  coupling_ele%value( beta_y0$) = extract_ele%y%beta
  coupling_ele%value(alpha_y0$) = extract_ele%y%alpha
  coupling_ele%value(  eta_y0$) = extract_ele%y%eta
  coupling_ele%value( etap_y0$) = extract_ele%y%etap
  
  coupling_ele%value( beta_x1$) = inject_ele%x%beta
  coupling_ele%value(alpha_x1$) = inject_ele%x%alpha
  coupling_ele%value(  eta_x1$) = inject_ele%x%eta
  coupling_ele%value( etap_x1$) = inject_ele%x%etap
  coupling_ele%value( beta_y1$) = inject_ele%y%beta
  coupling_ele%value(alpha_y1$) = inject_ele%y%alpha
  coupling_ele%value(  eta_y1$) = inject_ele%y%eta
  coupling_ele%value( etap_y1$) = inject_ele%y%etap
  
  coupling_ele%value(dphi_x$)   = mod(inject_ele%x%phi - extract_ele%x%phi,twopi)
  coupling_ele%value(dphi_y$)   = mod(inject_ele%y%phi - extract_ele%y%phi,twopi)
  
  ! it's a linear element so no orbit need be passed
  call make_mat6 (coupling_ele, s%u(u%coupling%from_uni)%design%param)
  
end subroutine  tao_match_lats_init
 
end module tao_lattice_calc_mod
