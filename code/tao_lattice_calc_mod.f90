!+
! Module: tao_lattice_calc_mod
!
! Lattice calculation routines are here. It's a module so that custom lattice
! calculations has access to the subroutines.

module tao_lattice_calc_mod

use tao_mod
use tao_data_mod
use macroparticle_mod

integer, parameter :: from_design$ = 1
integer, parameter :: from_model$ = 2

contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lattice_calc (init_design)
!
! Routine to calculate the lattice functions.
! 
! Input:
!  init_design      -- Logical: To initialize design lattices
!
! Output:
!-

subroutine tao_lattice_calc (init_design)

implicit none

logical, optional :: init_design

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
! only applicable for standard calculation
if (.not. initing_design) then
  call tao_set_var_useit_opt
  call tao_set_data_useit_opt
endif

! do a custom lattice calculation if desired
if (s%global%lattice_recalc) then
  if (initing_design) then
    do i = 1, size(s%u)
      call tao_lat_bookkeeper (s%u(i)%design)
      call tao_hook_lattice_calc (s%u(i), s%u(i)%design, s%u(i)%design_orb, &
                                                                   used(i))
    enddo
    if (.not. any(.not. used)) s%global%lattice_recalc = .false.
  else
    do i = 1, size(s%u)
      call tao_lat_bookkeeper (s%u(i)%model)
      call tao_hook_lattice_calc (s%u(i), s%u(i)%model, s%u(i)%model_orb, &
                                                                   used(i))
    enddo
    if (.not. any(.not. used)) s%global%lattice_recalc = .false.
  endif
endif
  
! Closed orbit and Twiss calculation.
! This can be slow for large lattices so only do it if the lattice changed.
if (s%global%lattice_recalc) then
  if (s%global%track_type .eq. 'single') then
    if (initing_design) then
      do i = 1, size(s%u)
        if (.not. used(i)) then
	  call tao_lat_bookkeeper (s%u(i)%design)
          ! set up matching element
	  call tao_match_lats_init (s%u(i))
	  call compute_element_energy (s%u(i)%design)
          call tao_inject_particle (s%u(i), s%u(i)%design, s%u(i)%design_orb, from_design$)
          call twiss_and_track (s%u(i)%design, s%u(i)%design_orb)
          if (s%u(i)%design%param%lost) &
            call out_io (s_blank$, r_name, "particle lost at element \I\.", &
                                                s%u(i)%design%param%ix_lost)
        endif
      enddo
    else
      do i = 1, size(s%u)
        if (.not. used(i)) then
	  call tao_lat_bookkeeper (s%u(i)%model)
	  call compute_element_energy (s%u(i)%model)
          call tao_inject_particle (s%u(i), s%u(i)%model, s%u(i)%model_orb, from_model$)
          call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
          if (s%u(i)%model%param%lost) &
            call out_io (s_blank$, r_name, "particle lost at element \I\.", &
                                                s%u(i)%model%param%ix_lost)
        endif
      enddo
    endif
    s%global%lattice_recalc = .false.
  elseif (s%global%track_type .eq. 'macro') then
    if (initing_design) then
      do i = 1, size(s%u)
        if (.not. used(i)) then
	  call tao_lat_bookkeeper (s%u(i)%design)
          ! set up matching element
	  call tao_match_lats_init (s%u(i))
	  call compute_element_energy (s%u(i)%design)
          call tao_inject_beam (s%u(i), s%u(i)%design, s%u(i)%design_orb, from_design$)
          call tao_macro_track (i, s%u(i)%design, s%u(i)%design_orb)
	endif
      enddo
    else
      do i = 1, size(s%u)
        if (.not. used(i)) then
	  call tao_lat_bookkeeper (s%u(i)%model)
	  call compute_element_energy (s%u(i)%model)
          call tao_inject_beam (s%u(i), s%u(i)%model, s%u(i)%model_orb, from_model$)
          call tao_macro_track (i, s%u(i)%model, s%u(i)%model_orb)
	endif
      enddo
    endif
    s%global%lattice_recalc = .false.
  else
    call out_io (s_fatal$, r_name, &
                   "This tracking type has yet to be implemented!")
    call out_io (s_blank$, r_name, &
                   "No tracking or twiss calculations will be perfomred.")
  endif
endif

! Transfer info from %model to %data arrays.
! only applicable for standard calculation
if (.not. initing_design) &
  call tao_load_data_array ()

end subroutine tao_lattice_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_macro_track (uni, lat, orb)

implicit none

integer uni
type (ring_struct) :: lat
type (coord_struct) :: orb(0:)
type (tao_universe_struct), pointer :: u
type (beam_struct), pointer :: beam
type (macro_init_struct), pointer :: macro_init
type (modes_struct) :: modes
type (coord_struct), allocatable :: temp_orb(:)

integer j, jj, ix1, ix2
integer, save :: ix_cache(n_universe_maxx) = 0
integer extract_at_ix_ele 

character(20) :: r_name = "macro_track"

  extract_at_ix_ele = -1

  u => s%u(uni)
  beam => u%beam%beam
  macro_init => u%beam%macro_init

  ! If beam is injected to here then no initialization wanted.
  if (.not. u%coupling%coupled) then
    if (lat%param%lattice_type .eq. circular_lattice$ .and. &
                                    .not. u%coupling%coupled) then
      ! set up radiation integral cache
      if (ix_cache(uni) .eq. 0) then
        call out_io (s_blank$, r_name, &
                     "Initializing radiation integeral cache for universe \I\ ", uni)
        call twiss_and_track (lat, temp_orb)
        orb = temp_orb
        call radiation_integrals (lat, orb, modes, ix_cache(uni))
      endif
      ! find closed orbit
      call twiss_and_track (lat, temp_orb)
      orb = temp_orb
      ! use this to find emittance
      call radiation_integrals (lat, orb, modes, ix_cache(uni))
      !transfer closed orbit info into macro_init
      macro_init%x%beta  = lat%ele_(0)%x%beta
      macro_init%x%alpha = lat%ele_(0)%x%alpha
      macro_init%x%emit  = modes%a%emittance
      macro_init%y%beta  = lat%ele_(0)%y%beta
      macro_init%y%alpha = lat%ele_(0)%y%alpha
      macro_init%y%emit  = modes%b%emittance
      macro_init%center  = orb(0)%vec
      macro_init%E_0 = lat%ele_(0)%value(beam_energy$)
      ! other macro_init parameters will be as in init.tao
    elseif (lat%param%lattice_type .eq. linear_lattice$) then
      ! set beam energy to start of linac
      macro_init%E_0 = lat%ele_(0)%value(beam_energy$)
      ! set initial beam centroid
      macro_init%center = orb(0)%vec
    else
      call out_io (s_error$, r_name, &
                   "This lattice type not yet implemented for macroparticles!")
      call err_exit
    endif
    
    call init_macro_distribution (beam, u%beam%macro_init, .true.)
  endif


  ! Find if injecting into another lattice
  inject_loop: do jj = uni+1, size(s%u)
    if (s%u(jj)%coupling%coupled) then
      if (s%u(jj)%coupling%from_uni .eq. uni) then
	if (s%u(jj)%coupling%from_uni_ix_ele .ne. -1) then
	  extract_at_ix_ele = s%u(jj)%coupling%from_uni_ix_ele
	  exit inject_loop ! save jj for coupled universe
	else
	  call out_io (s_abort$, r_name, &
	       "Must specify an element when coupling lattices with macroparticles")
	  call err_exit
	endif
      endif
    endif
  enddo inject_loop
  
  ix1 = 0
  ! track through every element
  do j = 1, lat%n_ele_use
    ix2 = j
  
    ! track to the element
    call track_beam (lat, beam, ix1, ix2)

    ! Save beam at location if injecting into another lattice
    if (extract_at_ix_ele .eq. j) then
      s%u(jj)%coupling%injecting_beam = beam
    endif
    
    ! compute centroid orbit
    call tao_find_beam_centroid (beam, orb(ix2), ix2, u%beam)
    
    ix1 = ix2
  enddo
  ! now find the standard transfer matrices and twiss parameters
  call ring_make_mat6 (lat, -1, orb)
  call twiss_propagate_all (lat)

end subroutine tao_macro_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Find the centroid of all macroparticles in all slices in all bunches
! Also keep track of lost macroparticles if the optional arguments are passed

subroutine tao_find_beam_centroid (beam, orb, ix_ele, u_beam)

implicit none

type (beam_struct), target :: beam
type (coord_struct) :: orb, coord
integer, optional :: ix_ele
type (tao_beam_struct), optional :: u_beam
type (macro_struct), pointer :: macro

real charge, tot_charge

integer n_bunch, n_slice, n_macro

logical record_lost

character(20) :: r_name = "tao_find_beam_centroid"

  ! are we checking keep track of lost particles?
  if (present(ix_ele) .and. present(u_beam)) then
    record_lost = .true.
  else
    record_lost = .false.
  endif

  tot_charge = 0.0
  coord%vec = 0.0

  bunch_loop: do n_bunch = 1, size(beam%bunch)
    slice_loop: do n_slice = 1, size(beam%bunch(n_bunch)%slice)
      macro_loop: do n_macro = 1, size(beam%bunch(n_bunch)%slice(n_slice)%macro)
        ! only average over particles that haven't been lost
	if (record_lost) then
          if (u_beam%ix_lost(n_bunch*n_slice*n_macro) .ne. -1) cycle macro_loop
	else
	  if (macro%lost) cycle macro_loop
	endif
        macro => beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)
	! check for lost macro through this element
	if (macro%lost .and. record_lost) then
	  u_beam%ix_lost(n_bunch*n_slice*n_macro) = ix_ele
          call out_io (s_blank$, r_name, "Macroparticle lost at element \I\.", &
                                               ix_ele )
	  cycle macro_loop
	endif
        charge = macro%charge
	tot_charge = tot_charge + charge
	coord%vec = coord%vec + charge * macro%r%vec
      enddo macro_loop
    enddo slice_loop
  enddo bunch_loop
      
  ! average

  if (tot_charge .ne. 0) then
    orb%vec = coord%vec / tot_charge
  else 
    ! lost all macros
    orb%vec = 0.0
  endif
 
end subroutine tao_find_beam_centroid

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a particle from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated

subroutine tao_inject_particle (u, lat, orbit, from_where)

implicit none

type (tao_universe_struct) u
type (ring_struct) lat
type (coord_struct) orbit(0:)
integer from_where

type (ele_struct), save :: extract_ele
type (coord_struct) pos

character(20) :: r_name = "inject_particle"

  if (.not. u%coupling%coupled) return

  call init_ele (extract_ele)
  
  ! get particle perameters from previous universe at position s
  if (from_where .eq. from_model$) then
    call twiss_and_track_at_s (s%u(u%coupling%from_uni)%model, &
                   u%coupling%from_uni_s, extract_ele, &
		   s%u(u%coupling%from_uni)%model_orb, pos)
  elseif (from_where .eq. from_design$) then
    call twiss_and_track_at_s (s%u(u%coupling%from_uni)%design, &
                   u%coupling%from_uni_s, extract_ele, &
		   s%u(u%coupling%from_uni)%design_orb, pos)
  else
    call out_io (s_error$, r_name, &
                "Unknown source for particle injection!")
  endif

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

subroutine tao_inject_beam (u, lat, orbit, from_where)

implicit none
 
type (tao_universe_struct) u
type (ring_struct) lat
type (coord_struct) orbit(0:)
integer from_where
type (ele_struct), save :: extract_ele

character(20) :: r_name = "inject_beam"

  if (.not. u%coupling%coupled) return
  
  ! beam from previous universe at end of element should already be set
  ! but we still need the twiss parameters and everything else
  if (from_where .eq. from_model$) then
    extract_ele = s%u(u%coupling%from_uni)%model%ele_(u%coupling%from_uni_ix_ele)
  elseif (from_where .eq. from_design$) then
    extract_ele = s%u(u%coupling%from_uni)%design%ele_(u%coupling%from_uni_ix_ele)
  else
    call out_io (s_error$, r_name, &
                "Unknown source for beam injection!")
  endif
  
  !track through coupling element
  if (u%coupling%use_coupling_ele) then
    call make_mat6 (u%coupling%coupling_ele, s%u(u%coupling%from_uni)%design%param)
    call twiss_propagate1 (extract_ele, u%coupling%coupling_ele)
    call track1_beam (u%coupling%injecting_beam, u%coupling%coupling_ele, &
              s%u(u%coupling%from_uni)%design%param, u%coupling%injecting_beam)
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
  u%beam%beam = u%coupling%injecting_beam
  call tao_find_beam_centroid (u%coupling%injecting_beam, orbit(0))

end subroutine tao_inject_beam

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
type (beam_struct), pointer :: beam

character(20) :: r_name = "match_lats_init"

  if (.not. (u%coupling%coupled .and. u%coupling%use_coupling_ele)) return

  coupling_ele => u%coupling%coupling_ele
  call init_ele (extract_ele)
  call init_ele (inject_ele)
  
  ! set up coupling%injecting_beam 
  if (s%global%track_type .eq. 'macro') then
    beam => u%beam%beam
    call reallocate_beam (u%coupling%injecting_beam, size(beam%bunch), &
            size(beam%bunch(1)%slice), size(beam%bunch(1)%slice(1)%macro))
  endif

  ! right now, will only match design lattices
  if (u%coupling%match_to_design) then
    ! get twiss parameters from extracted lattice
    if (s%global%track_type .eq. 'single') then
      call twiss_and_track_at_s (s%u(u%coupling%from_uni)%design, &
                     u%coupling%from_uni_s, extract_ele, &
          	   s%u(u%coupling%from_uni)%design_orb, extract_pos)
    elseif (s%global%track_type .eq. 'macro') then
      extract_ele = s%u(u%coupling%from_uni)%design%ele_(u%coupling%from_uni_ix_ele)
    endif
    
    ! get twiss parameters for injected lattice
    ! This is performed before the standard lattice calculation so the design
    ! twiss parameters in ele_(0) will still be as set in the lattice file.
    if (s%global%track_type .eq. 'single') then
      inject_ele = u%design%ele_(0)
    elseif (s%global%track_type .eq. 'macro') then
      inject_ele = u%design%ele_(0)
      inject_ele%x%beta  = u%beam%macro_init%x%beta
      inject_ele%x%alpha = u%beam%macro_init%x%alpha
      inject_ele%y%beta  = u%beam%macro_init%y%beta
      inject_ele%y%alpha = u%beam%macro_init%y%alpha
    endif
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
