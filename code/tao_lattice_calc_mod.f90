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

!

type (tao_lattice_struct), pointer :: this_bunch_track_lat ! for save_bunch_track routine

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_lattice_calc (calc_ok, ix_uni, init_design)
!
! Routine to calculate the lattice functions and TAO data. 
! Always tracks through the model lattice. 
! If initing the design lattices then the tracking will still be done
! through the model lattices. tao_init then transfers this into the
! design lattices.
! 
! Input:
!   ix_uni      -- Integer, optional: If present then calculation is restricted to
!                   the universe with this index.
!                   If present, tao_com%lattice_recalc will be ignored.
!   init_design -- Logical, optional: Set this True to initialize design lattices
!
! Output:
!   calc_ok -- Logical: Set False if there was an error in the 
!                calculation like a particle was lost or a lat is unstable.
!-

subroutine tao_lattice_calc (calc_ok, ix_uni, init_design)

implicit none

logical, optional :: init_design

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_dat
type (tao_d1_data_struct), pointer :: d1_dat
type (coord_struct), allocatable, save :: orb(:)
type (tao_lattice_struct), pointer :: tao_lat

integer, optional :: ix_uni
integer i, j, k, ix, n_max, it, id
real(rp) :: delta_e = 0

character(20) :: r_name = "tao_lattice_calc"

logical initing_design
logical :: calc_ok, recalc
logical this_calc_ok

!

calc_ok = .true.

initing_design = logic_option (.false., init_design)

! To save time, tao_com%lattice_recalc and s%u(:)%universe_recalc are used to
! determine what gets calculated. 
! If ix_uni is not present, the lattice functions of universe i are calculated
! if (tao_com%lattice_reclac = True and s%u(i)%universe_recalc = True. 

! do a custom lattice calculation if desired

if (.not. s%global%lattice_calc_on) return

recalc = tao_com%lattice_recalc .or. present(ix_uni)
if (.not. recalc) return

tao_com%ix0_taylor = -1   ! Reset taylor map
call tao_hook_lattice_calc (calc_ok)
    
! Closed orbit and Twiss calculation.
! This can be slow for large lattices so only do it if the lattice changed.

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (present(ix_uni)) then
    if (i /= ix_uni) cycle
  endif

  u => s%u(i)
  if (u%connect%connected) then
    if (s%u(u%connect%from_uni)%universe_recalc) u%universe_recalc = .true.
  endif

  if (.not. u%is_on .or. .not. u%universe_recalc) cycle

  ! Pointer to appropriate lattice and zero data array

  tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.
  u%data(:)%good_model = .false. ! reset
  u%data%model_value = tiny(1.0_rp)

  ! set up matching element

  if (initing_design) call tao_match_lats_init (u)

  call tao_lat_bookkeeper (u, tao_lat)

  ! 

  do k = 0, ubound(tao_lat%lat%branch, 1)
 
    do j = 1, 6
      tao_lat%lat_branch(k)%orbit%vec(j) = 0.0
    enddo

    select case (s%global%track_type)
    case ('single') 
      call tao_inject_particle (u, tao_lat, k)
      call tao_single_track (u, tao_lat, this_calc_ok, k)
    case ('beam')  ! Even when beam tracking we need to calculate the lattice parameters.
      call tao_inject_particle (u, tao_lat, k)
      call tao_single_track (u, tao_lat, this_calc_ok, k)
      call tao_inject_beam (u, tao_lat, k)
      call tao_beam_track (u, tao_lat, this_calc_ok, k)
    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING TYPE: ' // s%global%track_type)
      call err_exit
    end select

  enddo

  if (this_calc_ok) then
    if (u%do_synch_rad_int_calc) then
      call radiation_integrals (tao_lat%lat, tao_lat%lat_branch(0)%orbit, tao_lat%modes, u%ix_rad_int_cache)
      call transfer_rad_int_struct (ric, tao_lat%rad_int)
    endif
    if (u%do_chrom_calc) call chrom_calc (tao_lat%lat, delta_e, &
                             tao_lat%a%chrom, tao_lat%b%chrom, exit_on_error = .false.)
  else
    calc_ok = .false.
  endif

  call tao_load_data_array (u, -1, model$)

  ! do multi-turn tracking if needed. This is always the main lattice. 

  if (allocated(u%ix_data(-2)%ix_datum)) then
    ix = u%ix_data(-2)%ix_datum(1)
    d2_dat => u%data(ix)%d1%d2
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

  ! If calc is on common model then transfer data to base of all other universes

  if (tao_com%common_lattice .and. i == ix_common_uni$) then
    do j = lbound(s%u, 1), ubound(s%u, 1)
      if (j == ix_common_uni$) cycle
      s%u(j)%data(:)%base_value = u%data(:)%model_value
    enddo
  endif

enddo

! do any post-processing

call tao_hook_post_process_data ()
tao_com%lattice_recalc = .false.
s%u%universe_recalc = .true.

if (.not. initing_design) then
  call tao_set_var_useit_opt
  call tao_set_data_useit_opt
endif

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

integer i, ii, ix_lost, ix_branch

character(20) :: r_name = "tao_single_track"

logical calc_ok, err

!

orbit => tao_lat%lat_branch(ix_branch)%orbit

lat => tao_lat%lat
branch => tao_lat%lat%branch(ix_branch)
lat_branch => tao_lat%lat_branch(ix_branch)

calc_ok = .true.
lat%param%ix_lost = not_lost$
lat%param%lost = .false.

! Track

if (branch%param%lattice_type == circular_lattice$) then
  call closed_orbit_calc (lat, lat_branch%orbit, 4, exit_on_error = .false.)
  if (.not. bmad_status%ok) then
    calc_ok = .false.
    do i = 0, ubound(orbit, 1)
      orbit(i)%vec = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    enddo
  endif
  u%model_orb0 = orbit(0)   ! Save beginning orbit

else
  call track_all (lat, lat_branch%orbit, ix_branch)
  if (lat%param%lost) then
    calc_ok = .false.
    ix_lost = lat%param%ix_lost
    do ii = ix_lost+1, branch%n_ele_track
      orbit(ii)%vec = 0
    enddo
    call out_io (s_blank$, r_name, "particle lost at element \I0\: " // &
            trim(lat%ele(ix_lost)%name) // '  [s =\F9.2\]', &
            r_array = (/ lat%ele(ix_lost)%s /), i_array = (/ ix_lost /))
  endif
endif

! Twiss

if (u%mat6_recalc_on) then
  do i = 1, branch%n_ele_track
    if (lat%ele(i)%tracking_method == linear$) then
      call lat_make_mat6 (lat, i, ix_branch = ix_branch)
    else
      call lat_make_mat6 (lat, i, orbit, ix_branch)
    endif
  enddo
  if (branch%param%lattice_type == circular_lattice$) then
    call twiss_at_start (lat)
    if (.not. bmad_status%ok) then
      calc_ok = .false.
      return
    endif
  endif
  call twiss_propagate_all (lat, ix_branch)
endif

do i = 0, branch%n_ele_track
  call tao_load_data_array (u, i, model$)
enddo

end subroutine tao_single_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Right now, there is no beam tracking in circular lattices. 
! If extracting from a lat then the beam is generated at the extraction point.

subroutine tao_beam_track (u, tao_lat, calc_ok, ix_branch)

implicit none

type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
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

integer what_lat, n_alive, n_alive_old
integer i, j, i_uni, ip, ig, ic, ie1, ie2
integer n_bunch, n_part, i_uni_to, ix_lost
integer extract_at_ix_ele, n_lost, ix_branch
integer, allocatable, save :: ix_ele(:)

character(20) :: r_name = "tao_beam_track"

real(rp) :: value1, value2, f, time, old_time

logical post, calc_ok, too_many_lost, print_err, err, lost

! Initialize moment structure

this_bunch_track_lat => tao_lat
tao_lat%n_bunch_params2 = 0

!

calc_ok = .true.

call re_allocate (ix_ele,1)

beam => u%current_beam
branch => tao_lat%lat%branch(ix_branch)
lat_branch => tao_lat%lat_branch(ix_branch)
uni_branch => u%uni_branch(ix_branch)
uni_ele => uni_branch%ele
beam_init => uni_branch%beam_init
beam = uni_branch%beam0

lat => tao_lat%lat
lat%param%ix_lost = not_lost$  ! Needed by tao_evaluate_a_datum
ix_lost = not_lost$
lost = .false.
too_many_lost = .false.
orbit => tao_lat%lat_branch(ix_branch)%orbit

uni_ele(:)%n_particle_lost_here = 0

if (branch%param%lattice_type == circular_lattice$) then
  call out_io (s_fatal$, r_name, &
               'BEAM TRACKING WITH CIRCULAR LATTICE NOT YET IMPLEMENTED.', &
               'PLEASE SEE DCS IF YOU NEED THIS.')
  call err_exit
endif

! Transfer wakes from  design

do i = 1, branch%n_ele_max
  if (associated(lat%ele(i)%wake)) lat%ele(i)%wake%lr = u%design%lat%ele(i)%wake%lr
enddo

call zero_lr_wakes_in_lat (lat)

! Find if injecting into another lattice

extract_at_ix_ele = -1
do i_uni_to = u%ix_uni+1, ubound(s%u, 1)
  if (.not. s%u(i_uni_to)%connect%connected) cycle
  if (s%u(i_uni_to)%connect%from_uni /= u%ix_uni) cycle
  if (s%u(i_uni_to)%connect%from_uni_ix_ele == -1) then
    call out_io (s_abort$, r_name, &
         "Must specify an element when connecting lattices with a beam.")
    call err_exit
  endif
  extract_at_ix_ele = s%u(i_uni_to)%connect%from_uni_ix_ele
  exit ! save i_uni_to for connected universe
enddo

! If beam is injected into this lattice then no initialization wanted.

if (.not. u%connect%connected) then
  if (branch%param%lattice_type == circular_lattice$) then
    call tao_single_track (u, tao_lat, calc_ok, ix_branch) 
    if (extract_at_ix_ele /= -1) then
      if (u%calc_beam_emittance) then
        call radiation_integrals (lat, lat_branch%orbit, modes)
        f = lat%ele(extract_at_ix_ele)%value(E_tot$) * &
             (1+orbit(extract_at_ix_ele)%vec(6)) / mass_of(lat%param%particle)
        beam_init%a_norm_emitt  = modes%a%emittance * f
        beam_init%b_norm_emitt  = modes%b%emittance * f
      endif
      beam_init%center  = tao_lat%lat%beam_start%vec
      ! other beam_init parameters will be as in tao.init, or as above
      if (beam_init%a_norm_emitt == 0 .and. beam_init%b_norm_emitt == 0) then
        call out_io (s_abort$, r_name, &
                          'BOTH BEAM_INIT%A_NORM_EMITT AND %B_NORM_EMITT NOT SET!')
        call err_exit
      endif
      call init_beam_distribution (lat%ele(extract_at_ix_ele), &
                               beam_init, s%u(i_uni_to)%connect%injecting_beam)
    endif
    ! no beam tracking in circular lattices
    return
  endif
endif

! track through the lattice elements.
! The reference orbit for the transfer matrix calculation is taken to be the nominal 
! bunch center rather then the computed center to prevent jitter from changing things.

ie1 = uni_branch%ix_track_start
ie2 = branch%n_ele_track
if (uni_branch%ix_track_end > -1) ie2 = uni_branch%ix_track_end

print_err = .true.

if (s%global%beam_timer_on) then
  call run_timer ('START')
  old_time = 0
endif

n_alive_old = count(beam%bunch(s%global%bunch_to_plot)%particle(:)%ix_lost /= not_lost$)

do j = ie1, ie2

  ! track to the element and save for phase space plot

  if (tao_com%use_saved_beam_in_tracking) then
    beam = uni_ele(j)%beam

  else
    if (j /= ie1) then 
      call track_beam (lat, beam, j-1, j, too_many_lost)
    endif

    if (uni_ele(j)%save_beam) uni_ele(j)%beam = beam
  endif
 
  ! Save beam at location if injecting into another lattice
  if (extract_at_ix_ele == j) then
    s%u(i_uni_to)%connect%injecting_beam = beam
  endif

  ! compute centroid orbit

  call tao_find_beam_centroid (beam, orbit(j), too_many_lost, u, j, lat%ele(j))
  if (too_many_lost .and. .not. lost) then
    ix_lost = j
    lost = .true.
    call out_io (s_warn$, r_name, &
            "TOO MANY PARTICLES HAVE BEEN LOST AT ELEMENT #\i0\: " &
            // trim(lat%ele(j)%name), j)
  endif

  !

  n_alive = count(beam%bunch(s%global%bunch_to_plot)%particle(:)%ix_lost == not_lost$)
  uni_ele(j)%n_particle_lost_here = n_alive_old - n_alive
  n_alive_old = n_alive

  call calc_bunch_params (u%current_beam%bunch(s%global%bunch_to_plot), &
                             lat%ele(j), lat_branch%bunch_params(j), err, print_err)
  if (err) print_err = .false.  ! Only generate one message.
  call tao_load_data_array (u, j, model$) 

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

post = .false.
if (extract_at_ix_ele == -1) post = .true.
if (u%ix_uni < ubound(s%u, 1)) then
  if (.not. s%u(u%ix_uni+1)%is_on) post = .true.
endif

if (post) then
  n_lost = 0
  do n_bunch = 1, size(beam%bunch)
    n_lost = n_lost + count(beam%bunch(n_bunch)%particle%ix_lost /= not_lost$)
  enddo
  if (n_lost /= 0) &
    call out_io (s_blank$, r_name, &
      "Total number of lost particles by the end of universe \I2\: \I5\.", &
                                  i_array = (/u%ix_uni, n_lost /))
endif

lat%param%lost = lost
lat%param%ix_lost = ix_lost
  
end subroutine tao_beam_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Find the centroid of all particles in viewed bunches
! Also keep track of lost particles if the optional arguments are passed

subroutine tao_find_beam_centroid (beam, orb, too_many_lost, u, ix_ele, ele)

implicit none

type (beam_struct), target :: beam
type (coord_struct) :: orb, coord
type (ele_struct), optional :: ele
type (tao_universe_struct), optional :: u

integer, optional :: ix_ele

real(rp) charge_tot, charge

integer n_bunch, n_part, n_lost, i_ele

logical record_lost
logical too_many_lost

character(100) line
character(24) :: r_name = "tao_find_beam_centroid"

!

coord%vec = 0.0
n_lost = 0
charge_tot = 0
n_bunch = s%global%bunch_to_plot
  
! If optional arguments not present no verbose and
!  just check if particles lost

if (present(u) .or. present(ix_ele)) then
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
  elseif (beam%bunch(n_bunch)%particle(n_part)%ix_lost /= not_lost$) then
    cycle
  endif
  charge = beam%bunch(n_bunch)%particle(n_part)%charge
  coord%vec = coord%vec + beam%bunch(n_bunch)%particle(n_part)%r%vec * charge
  charge_tot = charge_tot + charge
enddo
      
! Post lost particles

if (record_lost .and. n_lost > 0) then
  line = "\I4\ particle(s) lost at element \I0\: " // ele%name
  if (size(s%u) > 1) line = trim(line) // " in universe \I3\ "
  call out_io (s_blank$, r_name, line, i_array = (/ n_lost, ix_ele, u%ix_uni /) )
endif
  
! average

if (charge_tot /= 0) then
  orb%vec = coord%vec / charge_tot
  too_many_lost = .false.
else 
  ! lost too many particles
  too_many_lost = .true.
  orb%vec = 0.0
endif
 
end subroutine tao_find_beam_centroid

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a particle from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated. If no injection then will set beginning orbit to
! whatever the user has specified. As always, tracking only occure in the model
! lattice.

subroutine tao_inject_particle (u, model, ix_branch)

implicit none

type (tao_universe_struct) u
type (tao_lattice_struct), target :: model
type (ele_struct), save :: extract_ele
type (ele_struct), pointer :: from_ele
type (coord_struct) pos
type (coord_struct), pointer :: orb0
type (spin_polar_struct) :: polar
type (branch_struct), pointer :: branch

integer ix_branch, i_ele_from, i_br_from

character(20) :: r_name = "inject_particle"

! Not main branch case

if (ix_branch /= 0) then
  branch => model%lat%branch(ix_branch)
  i_ele_from = branch%ix_from_ele
  i_br_from  = branch%ix_from_branch
  from_ele => model%lat%branch(i_br_from)%ele(i_ele_from)
  branch%ele(0)%x = from_ele%x
  branch%ele(0)%y = from_ele%y
  branch%ele(0)%a = from_ele%a
  branch%ele(0)%b = from_ele%b
  branch%ele(0)%z = from_ele%z
  branch%ele(0)%c_mat   = from_ele%c_mat
  branch%ele(0)%gamma_c = from_ele%gamma_c
  model%lat_branch(ix_branch)%orbit(0) = model%lat_branch(i_br_from)%orbit(i_ele_from)
  return
endif

! In u%model_orb0 is saved the last computed orbit. 
! This is important with common_lattice since tao_lat%lat_branch(0)%orbit(0) has been overwritten.

if (model%lat%branch(ix_branch)%param%lattice_type == linear_lattice$) then
  model%lat_branch(ix_branch)%orbit(0) = model%lat%beam_start
else
  model%lat_branch(ix_branch)%orbit(0) = u%model_orb0
endif

orb0 => model%lat_branch(ix_branch)%orbit(0)

! Not connected case is easy.

if (.not. u%connect%connected) then
  polar%theta = u%uni_branch(ix_branch)%beam_init%spin%theta
  polar%phi = u%uni_branch(ix_branch)%beam_init%spin%phi
  call polar_to_spinor (polar, orb0)
  return
endif

! Connected case

if (.not. s%u(u%connect%from_uni)%is_on) then
  call out_io (s_error$, r_name, &
                  "Injecting from a turned off universe! This will not do!", &
                  "No injection will be performed")
  return
endif
  
call init_ele (extract_ele)
  
! get particle perameters from previous universe at position s

call twiss_and_track_at_s (s%u(u%connect%from_uni)%model%lat, &
                             u%connect%from_uni_s, extract_ele, &
                             s%u(u%connect%from_uni)%model%lat_branch(0)%orbit, pos, ix_branch)

! track through connect element

if (u%connect%match_to_design) then
  if (u%mat6_recalc_on) call make_mat6 (u%connect%match_ele, &
                                             s%u(u%connect%from_uni)%model%lat%param)
  call twiss_propagate1 (extract_ele, u%connect%match_ele)
  call track1 (pos, u%connect%match_ele, &
                    s%u(u%connect%from_uni)%model%lat%param, pos)
  u%connect%match_ele%value(E_TOT$) = extract_ele%value(E_TOT$)
  u%connect%match_ele%floor = extract_ele%floor
  extract_ele = u%connect%match_ele
endif
  
! transfer to this lattice

model%lat%ele(0)%a = extract_ele%a
model%lat%ele(0)%b = extract_ele%b
model%lat%ele(0)%x = extract_ele%x
model%lat%ele(0)%y = extract_ele%y
model%lat%ele(0)%z = extract_ele%z
model%lat%ele(0)%value(E_TOT$) = extract_ele%value(E_TOT$)
model%lat%ele(0)%c_mat   = extract_ele%c_mat
model%lat%ele(0)%gamma_c = extract_ele%gamma_c
model%lat%ele(0)%floor   = extract_ele%floor
orb0      = pos
        
end subroutine tao_inject_particle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! This will inject a beam from a previous universe into this universe in
! preparation for tracking. The lattice where the extraction occurs will have
! already been calculated.
!
! If there is no connection between lattice then this will initialize the beam

subroutine tao_inject_beam (u, model, ix_branch)

use tao_read_beam_mod

implicit none

type (tao_universe_struct), target :: u
type (beam_init_struct) beam_init
type (tao_lattice_struct), target :: model
type (ele_struct), save :: extract_ele
type (coord_struct), pointer :: orb0
type (lat_param_struct), pointer :: param
type (tao_universe_branch_struct), pointer :: uni_branch

real(rp) v(6)
integer i, j, iu, ios, n_in_file, n_in, ix_branch

character(20) :: r_name = "tao_inject_beam"
character(100) line

logical too_many_lost, err

! Init

if (tao_com%use_saved_beam_in_tracking) return

uni_branch => u%uni_branch(ix_branch)
orb0 => model%lat_branch(0)%orbit(0)

model%lat%beam_start%vec = uni_branch%beam_init%center

! If there is an init file then read from the file

if (uni_branch%beam0_file /= "") then
  if (uni_branch%init_beam0 .or. .not. allocated(uni_branch%beam0%bunch)) then
    call tao_open_beam_file (uni_branch%beam0_file)
    call tao_set_beam_params (uni_branch%beam_init%n_bunch, uni_branch%beam_init%n_particle, &
                                                       uni_branch%beam_init%bunch_charge)
    call tao_read_beam (uni_branch%beam0, err)
    if (err) call err_exit
    call tao_close_beam_file()
    do i = 1, size(uni_branch%beam0%bunch)
      do j = 1, size(uni_branch%beam0%bunch(i)%particle)
        uni_branch%beam0%bunch(i)%particle(j)%r%vec = &
                  uni_branch%beam0%bunch(i)%particle(j)%r%vec + model%lat%beam_start%vec
      enddo
    enddo
    call out_io (s_info$, r_name, &
                  'Read initial beam distribution from: ' // uni_branch%beam0_file, &
                  'Centroid Offset: \6es12.3\ ', r_array = model%lat%beam_start%vec)
  endif
  call tao_find_beam_centroid (uni_branch%beam0, orb0, too_many_lost) 
  if (too_many_lost) then
    call out_io (s_warn$, r_name, "Not enough particles for beam init!")
    call err_exit
  endif
  if (uni_branch%ele(0)%save_beam) uni_branch%ele(0)%beam = uni_branch%beam0
  return
endif

! Not connected case

if (.not. u%connect%connected) then
  if (uni_branch%init_beam0 .or. .not. allocated(uni_branch%beam0%bunch)) then
    uni_branch%beam_init%center = model%lat%beam_start%vec
    if (uni_branch%beam_init%n_particle < 1) then
      call out_io (s_fatal$, r_name, &
        'BEAM_INIT INITIAL BEAM PROPERTIES NOT SET FOR UNIVERSE: \i4\ ', u%ix_uni)
      call err_exit
    endif
    beam_init = uni_branch%beam_init
    if (beam_init%n_bunch < 1) beam_init%n_bunch = 1   ! Default if not set.
    call init_beam_distribution (model%lat%ele(uni_branch%ix_track_start), beam_init, uni_branch%beam0)
    call tao_find_beam_centroid (uni_branch%beam0, orb0, too_many_lost)
    if (too_many_lost) then
      call out_io (s_warn$, r_name, "Not enough particles for beam init!")
      call err_exit
    endif
    uni_branch%init_beam0 = .false.
  endif
  if (uni_branch%ele(0)%save_beam) uni_branch%ele(0)%beam = uni_branch%beam0
  return
endif

! connected case...
   
if (.not. s%u(u%connect%from_uni)%is_on) then
  call out_io (s_error$, r_name, &
            "Injecting from a turned off universe! This will not do!", &
            "No injection will be performed.")
  return
endif
  
! beam from previous universe at end of extracting element should already be set
! but we still need the twiss parameters and everything else

extract_ele = s%u(u%connect%from_uni)%model%lat%ele(u%connect%from_uni_ix_ele)
  
! track through connection element

if (u%connect%match_to_design) then
  param => s%u(u%connect%from_uni)%model%lat%param
  if (u%mat6_recalc_on) call make_mat6 (u%connect%match_ele, param)
  call twiss_propagate1 (extract_ele, u%connect%match_ele)
  call track1_beam_simple (u%connect%injecting_beam, u%connect%match_ele, &
                      param, u%connect%injecting_beam)
  u%connect%match_ele%value(E_TOT$) = extract_ele%value(E_TOT$)
  u%connect%match_ele%floor = extract_ele%floor
  extract_ele = u%connect%match_ele
endif
    
! transfer to this lattice
model%lat%ele(0)%a = extract_ele%a
model%lat%ele(0)%b = extract_ele%b
model%lat%ele(0)%z = extract_ele%z
model%lat%ele(0)%value(E_TOT$) = extract_ele%value(E_TOT$)
model%lat%ele(0)%c_mat   = extract_ele%c_mat
model%lat%ele(0)%gamma_c = extract_ele%gamma_c
model%lat%ele(0)%floor   = extract_ele%floor
uni_branch%beam0 = u%connect%injecting_beam
call tao_find_beam_centroid (u%connect%injecting_beam, orb0, too_many_lost)
if (too_many_lost) then
  call out_io (s_warn$, r_name, "Not enough particles for beam init!")
  call err_exit
endif

if (uni_branch%ele(0)%save_beam) uni_branch%ele(0)%beam = uni_branch%beam0

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
type (ele_struct), pointer :: connection_ele
type (beam_struct), pointer :: injecting_beam

character(20) :: r_name = "match_lats_init"

!

if (.not. (u%connect%connected .and. u%connect%match_to_design)) return

connection_ele => u%connect%match_ele
call init_ele (extract_ele)
call init_ele (inject_ele)
  
! set up coupling%injecting_beam 

if (s%global%track_type == 'beam') then
  injecting_beam => s%u(u%connect%from_uni)%current_beam
  call reallocate_beam (u%connect%injecting_beam, size(injecting_beam%bunch), &
            size(injecting_beam%bunch(1)%particle))
endif

! match design lattices
if (u%connect%match_to_design) then
  ! get twiss parameters from extracted lattice
  if (s%global%track_type == 'single') then
    call twiss_and_track_at_s (s%u(u%connect%from_uni)%design%lat, &
                     u%connect%from_uni_s, extract_ele, &
                     s%u(u%connect%from_uni)%design%lat_branch(0)%orbit, extract_pos)
  elseif (s%global%track_type == 'beam') then
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
connection_ele%value(  eta_x0$) = extract_ele%x%eta
connection_ele%value( etap_x0$) = extract_ele%x%etap

connection_ele%value( beta_b0$) = extract_ele%b%beta
connection_ele%value(alpha_b0$) = extract_ele%b%alpha
connection_ele%value(  eta_y0$) = extract_ele%y%eta
connection_ele%value( etap_y0$) = extract_ele%y%etap
  
connection_ele%value( beta_a1$) = inject_ele%a%beta
connection_ele%value(alpha_a1$) = inject_ele%a%alpha
connection_ele%value(  eta_x1$) = inject_ele%x%eta
connection_ele%value( etap_x1$) = inject_ele%x%etap

connection_ele%value( beta_b1$) = inject_ele%b%beta
connection_ele%value(alpha_b1$) = inject_ele%b%alpha
connection_ele%value(  eta_y1$) = inject_ele%y%eta
connection_ele%value( etap_y1$) = inject_ele%y%etap
  
connection_ele%value(dphi_a$)   = mod(inject_ele%a%phi - extract_ele%a%phi,twopi)
connection_ele%value(dphi_b$)   = mod(inject_ele%b%phi - extract_ele%b%phi,twopi)
  
! it's a linear element so no orbit need be passed
if (u%mat6_recalc_on) call make_mat6 (connection_ele, &
                                             s%u(u%connect%from_uni)%design%lat%param)
  
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
type (bunch_params_struct), allocatable, save :: tm(:)
type (bunch_struct) bunch
type (ele_struct) ele

real(rp) s_travel
integer n

logical err

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
call calc_bunch_params (bunch, ele, t_lat%bunch_params2(n), err)

end subroutine
