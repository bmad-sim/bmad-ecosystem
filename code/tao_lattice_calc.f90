!+
! Subroutine tao_lattice_calc (universe, lattice, orbit)
!
! Routine to calculate the lattice functions.
! 
! Input:
!  universe  -- tao_universe_struct: OPTIONAL a speciffic universe to use
!                  if nothing specified then using s%u(:)
!  lattice   -- ring_struct: OPTIONAL, a speciffic lattice to use
!                  if nothing specified then using s%u(:)%model
!  orbit     -- coord_struct(0:): OPTIONAL, a speciffic orbit to place
!                  tracked orbit in
!                  if nothing specified then using s%u(:)%model_orb
!               NOTE: all three must be present or the default will be used
!
! Output:
!-

subroutine tao_lattice_calc (universe, lattice, orbit)

use tao_mod
use tao_data_mod
use macroparticle_mod

implicit none

type (tao_universe_struct), optional :: universe
type (ring_struct), optional :: lattice
type (coord_struct), optional :: orbit(0:)

type (coord_struct), allocatable :: orbit_temp(:)

integer i

character(20) :: r_name = "tao_lattice_calc"

! special means using speciffic lattice
logical special

! used refers to if a custom lattice calculation is performed
logical used_special
logical, automatic :: used(size(s%u))

special = .false.
used_special =.false.
used(:) = .false.

if (present(lattice) .and. present(orbit) .and. present(universe)) then
  special = .true.
elseif (present(universe)) then
  call out_io (s_warn$, r_name, "Must previde universe, lattice and orbit")
elseif (present(lattice)) then
  call out_io (s_warn$, r_name, "Must previde universe, lattice and orbit")
elseif (present(orbit)) then
  call out_io (s_warn$, r_name, "Must previde universe, lattice and orbit")
endif

! make sure useit is up-to-date
! only applicable for standard calculation
if (.not. special) then
  call tao_set_var_useit_opt
  call tao_set_data_useit_opt
endif

! do a custom lattice calculation if desired
if (s%global%lattice_recalc) then
  if (special) then
    call tao_hook_lattice_calc (universe, lattice, orbit, used_special)
    if (used_special) s%global%lattice_recalc = .false.
  else
    do i = 1, size(s%u)
      call tao_hook_lattice_calc (s%u(i), s%u(i)%model, s%u(i)%model_orb, &
                                                                   used(i))
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
    if (special) then
      call twiss_and_track (lattice, orbit_temp)
      call twiss_and_track (lattice, orbit_temp)
      orbit = orbit_temp
    else
      do i = 1, size(s%u)
        if (.not. used(i)) &
          call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
          call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
      enddo
    endif
    s%global%lattice_recalc = .false.
  elseif (s%global%track_type .eq. 'macro') then
    if (special) then
      call macro_track (universe, lattice, orbit)
      call macro_track (universe, lattice, orbit)
    else
      do i = 1, size(s%u)
        if (.not. used(i)) &
          call macro_track (s%u(i), s%u(i)%model, s%u(i)%model_orb)
          call macro_track (s%u(i), s%u(i)%model, s%u(i)%model_orb)
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
if (.not. special) &
  call tao_load_data_array ()


if (allocated(orbit_temp)) deallocate(orbit_temp)
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains

! right now this only does tracking through linacs

subroutine macro_track (u, lat, orb)

implicit none

type (tao_universe_struct), target :: u
type (ring_struct) :: lat
type (coord_struct) :: orb(0:)
type (beam_struct), pointer :: beam
integer, pointer :: to_ele(:)

integer j, jj, ix1, ix2

  beam => u%beam%beam
  to_ele => u%beam%to_ele

  ! set beam energy to start of linac
  u%beam%macro_init%E_0 = lat%ele_(0)%value(beam_energy$)

  call init_macro_distribution (beam, u%beam%macro_init, .true.)

  ! set initial beam centroid
  do j = 1, size(beam%bunch)
    do jj = 1, size(beam%bunch(j)%slice)
      beam%bunch(j)%slice(jj)%macro(:)%r = orb(0)
    enddo
  enddo
      
  if (u%beam%track_all) then
    ix1 = 0
    ! track through every element
    do j = 1, lat%n_ele_use
      ix2 = j
  
      ! track to the element
      call track_beam (lat, beam, ix1, ix2) 
      ! compute centroid orbit
      call find_beam_centroid (beam, orb(ix2))
      
      ix1 = ix2
    enddo
    ! now find the standard transfer matrices and twiss parameters
    call ring_make_mat6 (lat, -1, orb)
    call twiss_propagate_all (lat)
  endif

end subroutine macro_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!contains

subroutine find_beam_centroid (beam, orb)

implicit none

type (beam_struct), target :: beam
type (coord_struct) :: orb, coord
type (macro_struct), pointer :: macro

real charge, tot_charge

integer n_bunch, n_slice, n_macro, tot_macro

  tot_macro = 0
  tot_charge = 0.0
  coord%vec = 0.0

  do n_bunch = 1, size(beam%bunch)
    do n_slice = 1, size(beam%bunch(n_bunch)%slice)
      do n_macro = 1, size(beam%bunch(n_bunch)%slice(n_slice)%macro)
        macro => beam%bunch(n_bunch)%slice(n_slice)%macro(n_macro)
        charge = macro%charge
	coord%vec = coord%vec + charge * macro%r%vec
	tot_macro = tot_macro + 1
	tot_charge = tot_charge + charge
      enddo
    enddo
  enddo
      
  ! average thing

  orb%vec = coord%vec / (tot_macro * tot_charge)
 
end subroutine find_beam_centroid

end subroutine tao_lattice_calc
