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
use random_mod
use rad_int_common

!

type (tao_lattice_struct), pointer :: this_bunch_track_lat ! for save_bunch_track routine

integer, parameter :: design$ = 1
integer, parameter :: model$ = 2

logical all_lost_already

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lattice_calc (calc_ok, init_design)
!
! Routine to calculate the lattice functions and TAO data. 
! Always tracks through the model lattice. 
! If initing the design lattices then the tracking will still be done
! through the model lattices. tao_init then transfers this into the
! design lattices.
! 
! Input:
!   init_design -- Logical, optional: To initialize design lattices
!
! Output:
!   calc_ok -- Logical: Set False if there was an error in the 
!                calculation like a particle was lost or a lat is unstable.
!-

subroutine tao_lattice_calc (calc_ok, init_design)

implicit none

logical, optional :: init_design

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_dat
type (tao_d1_data_struct), pointer :: d1_dat
type (coord_struct), allocatable, save :: orb(:)

integer i, j, ix, n_max, it, id
real(rp) :: delta_e = 0

character(20) :: r_name = "tao_lattice_calc"

logical initing_design
! hook_used refers to if a custom lattice calculation is performed
logical, automatic :: hook_used(size(s%u))
logical :: calc_ok
logical this_calc_ok

!

calc_ok = .true.

hook_used(:) = .false.
initing_design = logic_option (.false., init_design)

! do a custom lattice calculation if desired

if (tao_com%lattice_recalc) then
  do i = 1, size(s%u)
    u => s%u(i)
    if (.not. u%is_on) cycle
    call tao_lat_bookkeeper (u, u%model)
    call tao_hook_lattice_calc (u, u%model, hook_used(i), calc_ok)
  enddo
  if (all(hook_used)) tao_com%lattice_recalc = .false.
endif
    
! Closed orbit and Twiss calculation.
! This can be slow for large lattices so only do it if the lattice changed.

if (.not. tao_com%lattice_recalc) return

do i = 1, size(s%u)
  u => s%u(i)
  u%data(:)%good_model = .false. ! reset

  if (.not. u%is_on .or. hook_used(i)) cycle
  ! zero data array
  u%data%model_value = tiny(1.0_rp)
  do j = 1, 6
    u%model%orb%vec(j) = 0.0
  enddo
  u%model%orb(0) = u%model%lat%beam_start

  ! set up matching element
  if (initing_design) call tao_match_lats_init (u)

  select case (s%global%track_type)
  case ('single') 
    call tao_inject_particle (u, u%model%lat, u%model%orb)
    call tao_lat_bookkeeper (u, u%model)
    call tao_single_track (i, u%model, this_calc_ok)
  case ('beam') 
    call tao_inject_beam (u, u%model%lat, u%model%orb)
    call tao_lat_bookkeeper (u, u%model)
    call tao_beam_track (i, u%model, this_calc_ok)
  case ('macro')
    call tao_inject_macro_beam (u, u%model%lat, u%model%orb)
    call tao_lat_bookkeeper (u, u%model)
    call tao_macro_track (i, u%model, this_calc_ok)
  case default
    call out_io (s_error$, r_name, &
                   "This tracking type has yet to be implemented!")
    call out_io (s_blank$, r_name, &
                   "No tracking or twiss calculations will be perfomred.")
  end select

  if (this_calc_ok) then
    if (u%do_synch_rad_int_calc) then
      call radiation_integrals (u%model%lat, u%model%orb, u%model%modes, u%ix_rad_int_cache)
      call transfer_rad_int_struct (ric, u%model%rad_int)
    endif
    if (u%do_chrom_calc) call chrom_calc (u%model%lat, delta_e, &
                                    u%model%a%chrom, u%model%b%chrom, exit_on_error = .false.)
  else
    calc_ok = .false.
  endif

  call tao_load_data_array (u, -1)

  ! do multi-turn tracking if needed

  if (associated(u%ix_data(-2)%ix_datum)) then
    ix = u%ix_data(-2)%ix_datum(1)
    d2_dat => u%data(ix)%d1%d2
    n_max = 0
    do id = 1, size(d2_dat%d1)
      n_max = max(n_max, ubound(d2_dat%d1(id)%d, 1))
    enddo
    call reallocate_coord (orb, u%model%lat%n_ele_max)
    orb(0) = u%model%lat%beam_start
    do it = 0, n_max
      do id = 1, size(d2_dat%d1)
        d1_dat => d2_dat%d1(id)
        if (it >= lbound(d1_dat%d, 1) .and. it <= ubound(d1_dat%d, 1)) then
          select case (d1_dat%name)
          case ('x');   d1_dat%d(it)%model_value = orb(0)%vec(1)
          case ('p_x'); d1_dat%d(it)%model_value = orb(0)%vec(2)
          case ('y');   d1_dat%d(it)%model_value = orb(0)%vec(3)
          case ('p_y'); d1_dat%d(it)%model_value = orb(0)%vec(4)
          case ('z');   d1_dat%d(it)%model_value = orb(0)%vec(5)
          case ('p_z'); d1_dat%d(it)%model_value = orb(0)%vec(6)
          case default
            call out_io (s_fatal$, r_name, 'BAD MULTI_TURN_ORBIT D1_DATA%NAME: ' // d1_dat%name)
            call err_exit
          end select
        endif
      enddo
      call track_all (u%model%lat, orb)
      orb(0) = orb(u%model%lat%n_ele_track)
    enddo
  endif

enddo

! do any post-processing

call tao_hook_post_process_data ()
tao_com%lattice_recalc = .false.
tao_com%init_beam0 = .false.

if (.not. initing_design) then
  call tao_set_var_useit_opt
  call tao_set_data_useit_opt
endif

end subroutine tao_lattice_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_single_track (uni, tao_lat, calc_ok)

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat

integer uni, i, ii

character(20) :: r_name = "tao_single_track"

logical calc_ok

!

lat => tao_lat%lat
calc_ok = .true.

if (lat%param%lattice_type == circular_lattice$) then
  call closed_orbit_calc (lat, tao_lat%orb, 4, exit_on_error = .false.)
  if (.not. bmad_status%ok) then
    calc_ok = .false.
    do i = 0, ubound(tao_lat%orb, 1)
      tao_lat%orb(i)%vec = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    enddo
  endif
  if (s%global%matrix_recalc_on) call lat_make_mat6 (lat, -1, tao_lat%orb)
  call twiss_at_start (lat)
  if (.not. bmad_status%ok) then
    calc_ok = .false.
    return
  endif
endif

lat%param%ix_lost = not_lost$
lat%param%lost = .false.

call tao_load_data_array (s%u(uni), 0)

do i = 0, lat%n_ele_track

  if (i /= 0) then
    ! if doing linear tracking, first compute transfer matrix
    if (s%global%matrix_recalc_on .and. lat%ele(i)%tracking_method == linear$) &
                                             call make_mat6 (lat%ele(i), lat%param) 

    call track1 (tao_lat%orb(i-1), lat%ele(i), lat%param, tao_lat%orb(i))
  endif

  if (lat%param%lost) then
    lat%param%ix_lost = i
    calc_ok = .false.
    do ii = i+1, lat%n_ele_track
      tao_lat%orb(ii)%vec = 0
    enddo
    call out_io (s_blank$, r_name, "particle lost at element \I\.", lat%param%ix_lost)
    return
  endif

  call tao_calc_params (s%u(uni), i)
  call tao_load_data_array (s%u(uni), i)
    
enddo
  

end subroutine tao_single_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Right now, there is no beam tracking in lattices. If extracting from a
! lat then the beam is generated at the extraction point.

subroutine tao_beam_track (uni, tao_lat, calc_ok)

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct), pointer :: u
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (normal_modes_struct) :: modes
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve

integer uni, what_lat
integer j, i_uni, ip, ig, ic
integer n_bunch, n_part, i_uni_to
integer extract_at_ix_ele, n_lost
integer, allocatable, save :: ix_ele(:)

character(20) :: r_name = "tao_beam_track"

real(rp) :: value1, value2, f

logical post, calc_ok

! Initialize moment structure

this_bunch_track_lat => tao_lat
tao_lat%n_bunch_params2 = 0

!

calc_ok = .true.

call re_allocate (ix_ele,1)

u => s%u(uni)
beam => u%current_beam
beam = u%beam0
beam_init => u%beam_init
lat => tao_lat%lat

all_lost_already = .false.

! Find if injecting into another lattice

extract_at_ix_ele = -1
do i_uni_to = uni+1, size(s%u)
  if (.not. s%u(i_uni_to)%connect%connected) cycle
  if (s%u(i_uni_to)%connect%from_uni /= uni) cycle
  if (s%u(i_uni_to)%connect%from_uni_ix_ele .ne. -1) then
    extract_at_ix_ele = s%u(i_uni_to)%connect%from_uni_ix_ele
    exit ! save i_uni_to for connected universe
  else
    call out_io (s_abort$, r_name, &
         "Must specify an element when connecting lattices with a beam.")
    ! set initial beam centroid
    call err_exit
  endif
enddo

! If beam is injected into this lattice then no initialization wanted.

if (.not. u%connect%connected) then
  ! no beam tracking in lattices
  if (lat%param%lattice_type == circular_lattice$) then
    call tao_single_track (uni, tao_lat, calc_ok) 
    if (u%calc_beam_emittance) then
      call radiation_integrals (lat, tao_lat%orb, modes)
      if (extract_at_ix_ele .ne. -1) then
        f = lat%ele(extract_at_ix_ele)%value(E_TOT$) * &
                    (1+tao_lat%orb(extract_at_ix_ele)%vec(6)) / mass_of(lat%param%particle)
        beam_init%a_norm_emitt  = modes%a%emittance * f
        beam_init%b_norm_emitt  = modes%b%emittance * f
      endif
    endif
    ! transfer extracted particle info into macro_init
    if (extract_at_ix_ele .ne. -1) then
      beam_init%center  = tao_lat%lat%beam_start%vec
      ! other beam_init parameters will be as in tao.init, or as above
      call init_beam_distribution (lat%ele(extract_at_ix_ele), &
                               beam_init, s%u(i_uni_to)%connect%injecting_beam)
    endif
    return
  elseif (lat%param%lattice_type .ne. linear_lattice$) then
    call out_io (s_error$, r_name, &
                   "This lattice type not yet implemented for beam tracking!")
    call err_exit
  endif
endif

! track through every element

do j = 0, lat%n_ele_track

  ! track to the element and save for phase space plot

  if (tao_com%use_saved_beam_in_tracking) then
    beam = u%ele(j)%beam

  else
    if (j /= 0) then 
      ! if doing linear tracking, first compute transfer matrix
      if (s%global%matrix_recalc_on .and. lat%ele(j)%tracking_method == linear$) &
             call make_mat6 (lat%ele(j), lat%param) 
      call track_beam (lat, beam, j-1, j)
    endif

    if (u%ele(j)%save_beam) u%ele(j)%beam = beam
  endif
 
  ! Save beam at location if injecting into another lattice
  if (extract_at_ix_ele == j) then
    s%u(i_uni_to)%connect%injecting_beam = beam
  endif

  ! compute centroid orbit
  call tao_find_beam_centroid (beam, tao_lat%orb(j), uni, j, lat%ele(j))

  if (all_lost_already) exit

  ! Find lattice and beam parameters and load data

  call tao_calc_params (u, j)    
  call tao_load_data_array (u, j) 
enddo

! only post total lost if no extraction or extracting to a turned off lattice

post = .false.
if (extract_at_ix_ele == -1) post = .true.
if (uni < size(s%u)) then
  if (.not. s%u(uni+1)%is_on) post = .true.
endif

if (post) then
  n_lost = 0 
  do n_bunch = 1, size(beam%bunch)
    n_lost = n_lost + count(beam%bunch(n_bunch)%particle%ix_lost /= not_lost$)
  enddo
  if (n_lost .ne. 0) &
    call out_io (s_blank$, r_name, &
      "Total number of lost particles by the end of universe \I2\: \I5\.", &
                                  i_array = (/uni, n_lost /))
endif
  
end subroutine tao_beam_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Right now, there is no macroparticle tracking in lattices. If extracting from a
! lat then the macroparticle beam is generated at the extraction point.

subroutine tao_macro_track (uni, tao_lat, calc_ok)

use macro_utils_mod

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (tao_universe_struct), pointer :: u
type (macro_beam_struct), pointer :: beam
type (macro_init_struct), pointer :: macro_init
type (normal_modes_struct) :: modes

integer uni, what_lat
integer j, i_uni
integer n_bunch, n_slice, n_macro
integer extract_at_ix_ele, n_lost

character(20) :: r_name = "macro_track"

real(rp) :: value1, value2, f

logical post, calc_ok

!

calc_ok = .true.

u => s%u(uni)
beam => u%macro_beam%beam
macro_init => u%macro_beam%macro_init
lat => tao_lat%lat

all_lost_already = .false.

! Find if injecting into another lattice
extract_at_ix_ele = -1
do i_uni = uni+1, size(s%u)
  if (.not. s%u(i_uni)%connect%connected) cycle
  if (s%u(i_uni)%connect%from_uni /= uni) cycle
  if (s%u(i_uni)%connect%from_uni_ix_ele .ne. -1) then
    extract_at_ix_ele = s%u(i_uni)%connect%from_uni_ix_ele
    exit ! save i_uni for connected universe
  else
    call out_io (s_abort$, r_name, &
         "Must specify an element when connecting lattices with macroparticles")
    call err_exit
  endif
enddo 
  
! If beam is injected to here then no initialization wanted.
if (.not. u%connect%connected) then
  ! no macroparticle tracking in lattices
  if (lat%param%lattice_type == circular_lattice$ ) then
    call tao_single_track (uni, tao_lat, calc_ok) 
    if (u%macro_beam%calc_emittance) then
      call radiation_integrals (lat, tao_lat%orb, modes, u%ix_rad_int_cache)
      if (extract_at_ix_ele .ne. -1) then
        f = lat%ele(extract_at_ix_ele)%value(E_TOT$) * &
                    (1+tao_lat%orb(extract_at_ix_ele)%vec(6)) / mass_of(lat%param%particle)
        macro_init%a%norm_emit  = modes%a%emittance * f
        macro_init%b%norm_emit  = modes%b%emittance * f
      endif
    endif
    !transfer extracted particle info into macro_init
    if (extract_at_ix_ele .ne. -1) then
      macro_init%center  = tao_lat%lat%beam_start%vec
      ! other macro_init parameters will be as in init.tao, or as above
      call init_macro_distribution (s%u(i_uni)%connect%injecting_macro_beam, &
                                macro_init, lat%ele(extract_at_ix_ele), .true.)
    endif
    return
  elseif (lat%param%lattice_type .ne. linear_lattice$) then
    call out_io (s_error$, r_name, &
                   "This lattice type not yet implemented for macroparticles!")
    call err_exit
  endif
endif

! beginning element calculations
call tao_load_data_array (u, 0) 

! track through every element
do j = 1, lat%n_ele_track

  ! if doing linear tracking, first compute transfer matrix
  if (s%global%matrix_recalc_on .and. lat%ele(j)%tracking_method == linear$) &
           call make_mat6 (lat%ele(j), lat%param) 

  ! track to the element
  if (j /= 0) call track_macro_beam (lat, beam, j-1, j)
 
  ! Save beam at location if injecting into another lattice
  if (extract_at_ix_ele == j) then
    s%u(i_uni)%connect%injecting_macro_beam = beam
  endif
    
  ! compute centroid orbit
  call tao_find_macro_beam_centroid (beam, tao_lat%orb(j), uni, j, u%macro_beam)

  if (all_lost_already) exit

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
  if (n_lost .ne. 0) &
    call out_io (s_blank$, r_name, &
      "Total number of lost macroparticles by the end of universe \I2\: \I5\.", &
                                                      i_array = (/uni, n_lost /))
endif
  
end subroutine tao_macro_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Find the centroid of all particles in viewed bunches
! Also keep track of lost macroparticles if the optional arguments are passed

subroutine tao_find_beam_centroid (beam, orb, uni, ix_ele, ele)

implicit none

type (beam_struct), target :: beam
type (coord_struct) :: orb, coord
type (ele_struct), optional :: ele

integer, optional :: uni, ix_ele

integer n_bunch, n_part, n_lost, tot_part, i_ele

logical record_lost

character(100) line
character(20) :: r_name = "tao_find_beam_centroid"

!

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
if (record_lost .and. n_lost > 0) then
  line = "\I4\ particle(s) lost at element \I6\: " // ele%name
  if (size(s%u) > 1) line = trim(line) // " in universe \I3\ "
  call out_io (s_blank$, r_name, line, i_array = (/ n_lost, ix_ele, uni /) )
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

subroutine tao_inject_particle (u, lat, orb)

implicit none

type (tao_universe_struct) u
type (lat_struct) lat
type (coord_struct) orb(0:)

type (ele_struct), save :: extract_ele
type (coord_struct) pos
type (spin_polar_struct) :: polar

character(20) :: r_name = "inject_particle"

!

if (.not. u%connect%connected) then
  orb(0) = u%model%orb(0)
  polar%theta = u%beam_init%spin%theta
  polar%phi = u%beam_init%spin%phi

  call polar_to_spinor (polar, orb(0))
else
    
  if (.not. s%u(u%connect%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed")
    return
  endif
  
  call init_ele (extract_ele)
  
  ! get particle perameters from previous universe at position s
  call twiss_and_track_at_s (s%u(u%connect%from_uni)%model%lat, &
                             u%connect%from_uni_s, extract_ele, &
                             s%u(u%connect%from_uni)%model%orb, pos)

  ! track through connect element
  if (u%connect%use_connect_ele) then
    if (s%global%matrix_recalc_on) call make_mat6 (u%connect%connect_ele, &
                                             s%u(u%connect%from_uni)%model%lat%param)
    call twiss_propagate1 (extract_ele, u%connect%connect_ele)
    call track1 (pos, u%connect%connect_ele, &
                    s%u(u%connect%from_uni)%model%lat%param, pos)
    u%connect%connect_ele%value(E_TOT$) = extract_ele%value(E_TOT$)
    u%connect%connect_ele%floor = extract_ele%floor
    extract_ele = u%connect%connect_ele
  endif
  
  ! transfer to this lattice
  lat%ele(0)%a = extract_ele%a
  lat%ele(0)%b = extract_ele%b
  lat%ele(0)%x = extract_ele%x
  lat%ele(0)%y = extract_ele%y
  lat%ele(0)%z = extract_ele%z
  lat%ele(0)%value(E_TOT$) = extract_ele%value(E_TOT$)
  lat%ele(0)%c_mat   = extract_ele%c_mat
  lat%ele(0)%gamma_c = extract_ele%gamma_c
  lat%ele(0)%floor   = extract_ele%floor
  orb(0)      = pos
endif
        
end subroutine tao_inject_particle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a beam from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated.
!
! If there is no connection between lattice then this will initialize the beam

subroutine tao_inject_beam (u, lat, orb)

use tao_read_beam_mod

implicit none

type (lat_struct) lat
type (tao_universe_struct) u
type (coord_struct) orb(0:)
type (ele_struct), save :: extract_ele

type (lat_param_struct), pointer :: param

real(rp) v(6), bunch_charge
integer i, iu, ios, n_bunch, n_particle, n_in_file, n_in

character(20) :: r_name = "tao_inject_beam"
character(100) line

!

if (tao_com%use_saved_beam_in_tracking) return

! If there is an init file then read from the file

if (u%beam0_file /= "") then

  call open_beam_file (u%beam0_file)
  call read_beam_params (i, n_bunch, n_particle, bunch_charge)
  call set_beam_params (u%beam_init%n_bunch, u%beam_init%n_particle, u%beam_init%bunch_charge)
  call read_beam (u%beam0)
  call close_beam_file()
  call out_io (s_info$, r_name, 'Read initial beam distribution from: ' // u%beam0_file)

  call tao_find_beam_centroid (u%beam0, orb(0))  

elseif (.not. u%connect%connected) then
  if (tao_com%init_beam0 .or. .not. allocated(u%beam0%bunch)) then
    u%beam_init%center = u%model%lat%beam_start%vec
    if (u%beam_init%n_bunch < 1 .or. u%beam_init%n_particle < 1) then
      call out_io (s_fatal$, r_name, &
        'BEAM_INIT INITIAL BEAM PROPERTIES NOT SET FOR UNIVERSE: \i4\ ', u%ix_uni)
      call err_exit
    endif
    call init_beam_distribution (lat%ele(0), u%beam_init, u%beam0)
    call tao_find_beam_centroid (u%beam0, orb(0))
  endif

else
   
  if (.not. s%u(u%connect%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed.")
    return
  endif
  
  ! beam from previous universe at end of extracting element should already be set
  ! but we still need the twiss parameters and everything else

  extract_ele = s%u(u%connect%from_uni)%model%lat%ele(u%connect%from_uni_ix_ele)
  
  ! track through connection element

  if (u%connect%use_connect_ele) then
    param => s%u(u%connect%from_uni)%model%lat%param
    if (s%global%matrix_recalc_on) call make_mat6 (u%connect%connect_ele, param)
    call twiss_propagate1 (extract_ele, u%connect%connect_ele)
    call track1_beam_ele (u%connect%injecting_beam, u%connect%connect_ele, &
                        param, u%connect%injecting_beam)
    u%connect%connect_ele%value(E_TOT$) = extract_ele%value(E_TOT$)
    u%connect%connect_ele%floor = extract_ele%floor
    extract_ele = u%connect%connect_ele
  endif
    
  ! transfer to this lattice
  lat%ele(0)%a = extract_ele%a
  lat%ele(0)%b = extract_ele%b
  lat%ele(0)%z = extract_ele%z
  lat%ele(0)%value(E_TOT$) = extract_ele%value(E_TOT$)
  lat%ele(0)%c_mat   = extract_ele%c_mat
  lat%ele(0)%gamma_c = extract_ele%gamma_c
  lat%ele(0)%floor   = extract_ele%floor
  u%beam0 = u%connect%injecting_beam
  call tao_find_beam_centroid (u%connect%injecting_beam, orb(0))
endif

! record this beam
  
if (u%ele(0)%save_beam) u%ele(0)%beam = u%beam0

end subroutine tao_inject_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a beam from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated.
!
! If there is no connection between lattice then this will initialize the beam

subroutine tao_inject_macro_beam (u, lat, orb)

implicit none
 
type (tao_universe_struct) u
type (coord_struct) orb(0:)
type (lat_struct) lat
type (ele_struct), save :: extract_ele

type (lat_param_struct), pointer :: param

integer n_bunch, n_slice, n_macro

character(20) :: r_name = "tao_inject_macro_beam"

!

if (.not. u%connect%connected) then
  u%macro_beam%macro_init%center = u%model%lat%beam_start%vec
  call init_macro_distribution (u%macro_beam%beam, u%macro_beam%macro_init, &
                                lat%ele(0), .true.)
  u%macro_beam%ix_lost(:,:,:) = not_lost$
  call tao_find_macro_beam_centroid (u%macro_beam%beam, orb(0))
else
  if (.not. s%u(u%connect%from_uni)%is_on) then
    call out_io (s_error$, r_name, &
      "Injecting from a turned off universe! This will not do!")
    call out_io (s_blank$, r_name, &
      "No injection will be performed")
    return
  endif
  
  ! beam from previous universe at end of element should already be set
  ! but we still need the twiss parameters and everything else
  extract_ele = s%u(u%connect%from_uni)%model%lat%ele(u%connect%from_uni_ix_ele)
  
  !track through connection element
  if (u%connect%use_connect_ele) then
    param => s%u(u%connect%from_uni)%model%lat%param
    if (s%global%matrix_recalc_on) call make_mat6 (u%connect%connect_ele, param)
    call twiss_propagate1 (extract_ele, u%connect%connect_ele)
    call track1_macro_beam (u%connect%injecting_macro_beam, u%connect%connect_ele, &
                                              param, u%connect%injecting_macro_beam)
    u%connect%connect_ele%value(E_TOT$) = extract_ele%value(E_TOT$)
    u%connect%connect_ele%floor = extract_ele%floor
    extract_ele = u%connect%connect_ele
  endif
    
  ! transfer to this lattice
  lat%ele(0)%a = extract_ele%a
  lat%ele(0)%b = extract_ele%b
  lat%ele(0)%z = extract_ele%z
  lat%ele(0)%value(E_TOT$) = extract_ele%value(E_TOT$)
  lat%ele(0)%c_mat   = extract_ele%c_mat
  lat%ele(0)%gamma_c = extract_ele%gamma_c
  lat%ele(0)%floor   = extract_ele%floor
  u%macro_beam%beam = u%connect%injecting_macro_beam
  call tao_find_macro_beam_centroid (u%connect%injecting_macro_beam, orb(0))
  
  ! deterine if macroparticle already lost
  do n_bunch = 1, size(u%connect%injecting_macro_beam%bunch)
    do n_slice = 1, size(u%connect%injecting_macro_beam%bunch(n_bunch)%slice)
      do n_macro = 1, size(u%connect%injecting_macro_beam%bunch(n_bunch)%slice(n_slice)%macro)
  if (u%connect%injecting_macro_beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)%lost) &
    u%macro_beam%ix_lost(n_bunch, n_slice, n_macro) = 0
      enddo
    enddo
  enddo
endif
  
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
type (ele_struct), pointer :: connection_ele
type (beam_struct), pointer :: injecting_beam
type (macro_beam_struct), pointer :: injecting_macro_beam

character(20) :: r_name = "match_lats_init"

!

if (.not. (u%connect%connected .and. u%connect%use_connect_ele)) return

connection_ele => u%connect%connect_ele
call init_ele (extract_ele)
call init_ele (inject_ele)
  
! set up coupling%injecting_beam or macro_beam 
if (s%global%track_type == 'beam') then
  injecting_beam => s%u(u%connect%from_uni)%current_beam
  call reallocate_beam (u%connect%injecting_beam, size(injecting_beam%bunch), &
            size(injecting_beam%bunch(1)%particle))
elseif (s%global%track_type == 'macro') then
  injecting_macro_beam => s%u(u%connect%from_uni)%macro_beam%beam
  call reallocate_macro_beam (u%connect%injecting_macro_beam, &
                                size(injecting_macro_beam%bunch), &
                                size(injecting_macro_beam%bunch(1)%slice), &
                                size(injecting_macro_beam%bunch(1)%slice(1)%macro))
endif

! match design lattices
if (u%connect%match_to_design) then
  ! get twiss parameters from extracted lattice
  if (s%global%track_type == 'single') then
    call twiss_and_track_at_s (s%u(u%connect%from_uni)%design%lat, &
                     u%connect%from_uni_s, extract_ele, &
                     s%u(u%connect%from_uni)%design%orb, extract_pos)
  elseif (s%global%track_type == 'beam') then
    extract_ele = s%u(u%connect%from_uni)%design%lat%ele(u%connect%from_uni_ix_ele)
  elseif (s%global%track_type == 'macro') then
    extract_ele = s%u(u%connect%from_uni)%design%lat%ele(u%connect%from_uni_ix_ele)
  endif
    
  ! get twiss parameters for injected lattice
  ! This is performed before the standard lattice calculation so the design
  ! twiss parameters in ele(0) will still be as set in the lattice file.
  inject_ele = u%design%lat%ele(0)
else
  call out_io (s_warn$, r_name, &
                "The coupling element will only match the design lattices")
  return
endif

! set up matching element
connection_ele%value( beta_a0$) = extract_ele%a%beta
connection_ele%value(alpha_a0$) = extract_ele%a%alpha
connection_ele%value(  eta_a0$) = extract_ele%a%eta
connection_ele%value( etap_a0$) = extract_ele%a%etap
connection_ele%value( beta_b0$) = extract_ele%b%beta
connection_ele%value(alpha_b0$) = extract_ele%b%alpha
connection_ele%value(  eta_b0$) = extract_ele%b%eta
connection_ele%value( etap_b0$) = extract_ele%b%etap
  
connection_ele%value( beta_a1$) = inject_ele%a%beta
connection_ele%value(alpha_a1$) = inject_ele%a%alpha
connection_ele%value(  eta_a1$) = inject_ele%a%eta
connection_ele%value( etap_a1$) = inject_ele%a%etap
connection_ele%value( beta_b1$) = inject_ele%b%beta
connection_ele%value(alpha_b1$) = inject_ele%b%alpha
connection_ele%value(  eta_b1$) = inject_ele%b%eta
connection_ele%value( etap_b1$) = inject_ele%b%etap
  
connection_ele%value(dphi_a$)   = mod(inject_ele%a%phi - extract_ele%a%phi,twopi)
connection_ele%value(dphi_b$)   = mod(inject_ele%b%phi - extract_ele%b%phi,twopi)
  
! it's a linear element so no orbit need be passed
if (s%global%matrix_recalc_on) call make_mat6 (connection_ele, s%u(u%connect%from_uni)%design%lat%param)
  
end subroutine  tao_match_lats_init
 
end module tao_lattice_calc_mod

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine save_bunch_track (bunch, ele, s_travel)
!
! Custom routine to record bunch statistics.
! This routine is called by the bmad routine track1_bunch_csr.
!-

subroutine save_bunch_track (bunch, ele, s_travel)

use tao_lattice_calc_mod

implicit none

type (tao_lattice_struct), pointer :: t_lat
type (bunch_params_struct), allocatable :: tm(:)
type (bunch_struct) bunch
type (ele_struct) ele

real(rp) s_travel
integer n

! Make sure we have room in the array for the data

t_lat => this_bunch_track_lat
n = t_lat%n_bunch_params2

if (allocated(t_lat%bunch_params2)) then
  if (n+1 > size(t_lat%bunch_params2)) then
    allocate(tm(n))
    tm = t_lat%bunch_params2
    deallocate (t_lat%bunch_params2)
    allocate (t_lat%bunch_params2(max(10, 2*n)))
    t_lat%bunch_params2(1:n) = tm
    deallocate (tm)
  endif

else
  allocate(t_lat%bunch_params2(10))
endif

!

n = n + 1
t_lat%n_bunch_params2 = n
t_lat%bunch_params2(n)%s = ele%s - ele%value(l$) + s_travel
call calc_bunch_params (bunch, ele, t_lat%bunch_params2(n))

end subroutine
