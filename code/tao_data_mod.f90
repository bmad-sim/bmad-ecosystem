module tao_data_mod

use tao_mod

! These are data types specific to macroparticles

type macro_d2_data_names_struct
  character(16) :: name
  character(16), pointer :: d1_data(:)
end type

type this_coupling_struct
  real(rp) cbar(2,2)
  real(rp) coupling11, coupling12a, coupling12b, coupling22
  real(rp) f_11, f_12a, f_12b, f_22
  logical calc_done
end type

type (this_coupling_struct), save, allocatable, target :: cc(:)

type (macro_d2_data_names_struct), save, pointer :: macro_data_names(:)

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_load_data_array ()
!
! Routine to take data from the model lattice and model orbit
! and put that into the s%u(:)%data(:) arrays.
!
! Input:
!   s  -- Super_universe_struct:
!-

subroutine tao_load_data_array()

implicit none

type (tao_universe_struct), pointer :: u
integer i, j
character(20) :: r_name = 'tao_load_data_array'

! loop over all data in all universes 

do i = 1, size(s%u)
  u => s%u(i)
  call tao_data_coupling_init (u) 
  do j = 1, size(u%data)
    if (.not. u%data(j)%exists) cycle
    call tao_evaluate_a_datum (u%data(j), u%model, u%model_orb, &
                                                  u%data(j)%model_value)
    u%data(j)%s = u%model%ele_(u%data(j)%ix_ele)%s
  enddo
enddo

! do any post-processing
call tao_hook_post_process_data ()

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_data_coupling_init (u)
!
! Subroutine to initialize the coupling structure for a universe.
! This subroutine should be called before tao_evaluate_a_datum is
! called for a particular universe.
!
! Input:
!   u -- Tao_universe_struct: New universe.
!-

subroutine tao_data_coupling_init (u)

type (tao_universe_struct) u
integer m

! 

m = u%model%n_ele_max
if (.not. allocated(cc)) allocate (cc(0:m))
if (ubound(cc, 1) < m) then
  deallocate(cc)
  allocate(cc(0:m))
endif

cc%calc_done = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_evaluate_a_datum (datum, lattice, orb, datum_value)
!
! Subroutine to 
!-

subroutine tao_evaluate_a_datum (datum, lattice, orb, datum_value)

type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: data
type (tao_data_struct) datum
type (ring_struct) lattice
type (coord_struct) orb(0:)
type (modes_struct) mode

real(rp) datum_value
integer i, j, k, m, ix, ix1, ix2
character(20) :: r_name = 'tao_evaluate_a_datum'
logical found

!

call tao_hook_load_data_array (found, datum, lattice, orb, datum_value)
if (found) return

ix1 = datum%ix_ele
ix2 = datum%ix_ele2
ele => lattice%ele_(ix1)

select case (datum%data_type)

case ('orbit:x')
  call load_it (orb(:)%vec(1))
case ('orbit:y')
  call load_it (orb(:)%vec(3))
case ('orbit:z')
  call load_it (orb(:)%vec(5))

case ('orbit:x_p')
  call load_it (orb(:)%vec(2))
case ('orbit:y_p')
  call load_it (orb(:)%vec(4))
case ('orbit:z_p')
  call load_it (orb(:)%vec(6))

case ('phase:x')
  call relative_switch
  datum_value = lattice%ele_(ix2)%x%phi - lattice%ele_(ix1)%x%phi
case ('phase:y')
  call relative_switch
  datum_value = lattice%ele_(ix2)%y%phi - lattice%ele_(ix1)%y%phi

case ('beta:x')
  call load_it (lattice%ele_(:)%x%beta)
case ('beta:y')
  call load_it (lattice%ele_(:)%y%beta)

case ('alpha:x')
  call load_it (lattice%ele_(:)%x%alpha)
case ('alpha:y')
  call load_it (lattice%ele_(:)%y%alpha)

case ('eta:x')
  call load_it (lattice%ele_(:)%x%eta)
case ('eta:y')
  call load_it (lattice%ele_(:)%y%eta)

case ('etap:x')
  call load_it (lattice%ele_(:)%x%etap)
case ('etap:y')
  call load_it (lattice%ele_(:)%y%etap)

case ('beam_energy')
  call load_it (lattice%ele_(:)%value(beam_energy$))
  
case ('coupling:11b')
  call load_it (cc%coupling11, cc%f_11, coupling_here = .true.)
case ('coupling:12a')
  call load_it (cc%coupling12a, cc%f_12a, coupling_here = .true.)
case ('coupling:12b')
  call load_it (cc%coupling12b, cc%f_12b, coupling_here = .true.)
case ('coupling:22a')
  call load_it (cc%coupling22, cc%f_22, coupling_here = .true.)

case ('cbar:11')
  call load_it (cc%cbar(1,1), coupling_here = .true.)
case ('cbar:12')
  call load_it (cc%cbar(1,2), coupling_here = .true.)
case ('cbar:21')
  call load_it (cc%cbar(1,2), coupling_here = .true.)
case ('cbar:22')
  call load_it (cc%cbar(2,2), coupling_here = .true.)

case ('i5a_e6')
  call radiation_integrals (lattice, orb, mode)
  datum_value = mode%lin%i5a_e6
  datum%ix_ele_merit = lattice%n_ele_use

case ('i5b_e6')
  call radiation_integrals (lattice, orb, mode)
  datum_value = mode%lin%i5b_e6
  datum%ix_ele_merit = lattice%n_ele_use

case ('r56')
  call out_io (s_abort$, r_name, 'r56 constraint not yet implemented')
  call err_exit

case ('t566')
  call out_io (s_abort$, r_name, 't566 constraint not yet implemented')
  call err_exit

case ('floor:x')
  call relative_switch
  datum_value = lattice%ele_(ix2)%floor%x - lattice%ele_(ix1)%floor%x
case ('floor:y')
  call relative_switch
  datum_value = lattice%ele_(ix2)%floor%y - lattice%ele_(ix1)%floor%y
case ('floor:z')
  call relative_switch
  datum_value = lattice%ele_(ix2)%floor%z - lattice%ele_(ix1)%floor%z
case ('floor:theta')
  call relative_switch
  datum_value = lattice%ele_(ix2)%floor%theta - lattice%ele_(ix1)%floor%theta

case ('s_position')
  call relative_switch
  datum_value = lattice%ele_(ix2)%s - lattice%ele_(ix1)%s


case default
  call out_io (s_error$, r_name, 'UNKNOWN DATUM TYPE: ' // datum%data_type)
  return

end select

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine load_it (vec, f, coupling_here)

real(rp) vec(0:)
real(rp), optional :: f(0:)
integer ix_m, i
logical, optional :: coupling_here

!

if (ix2 < 0) then
  ix_m = ix1
  if (present(coupling_here)) call coupling_calc (ix_m)

!

else

  if (present(coupling_here)) then
    do i = ix1, ix2
      call coupling_calc (i)
    enddo
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec(ix1:ix2), 1) + ix1 - 1
  case ('max')
    ix_m = maxloc (vec(ix1:ix2), 1) + ix1 - 1
  case ('abs_min')
    ix_m = minloc (abs(vec(ix1:ix2)), 1) + ix1 - 1
  case ('abs_max')
    ix_m = maxloc (abs(vec(ix1:ix2)), 1) + ix1 - 1
  case default
    call out_io (s_abort$, r_name, 'BAD MERIT_TYPE: ' // datum%merit_type, &
                                   'FOR DATUM: ' // datum%data_type)
    call err_exit
  end select
endif

!

datum%ix_ele_merit = ix_m
datum_value = vec(ix_m)
if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ix_m))
if (present(f)) datum%conversion_factor = f(ix_m)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! contains

subroutine relative_switch

if (ix2 < 0) then
  ix2 = ix1
  ix1 = 0
endif  

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! contains

subroutine coupling_calc (ixd)

type (ele_struct), pointer :: ele
type (this_coupling_struct), pointer :: cc_p
integer ie, ixd
real(rp) f, f1, f2

!

if (cc(ixd)%calc_done) return

cc_p => cc(ixd)
ie = datum%ix_ele
ele => lattice%ele_(ie)

call c_to_cbar (ele, cc_p%cbar)
f = sqrt(ele%x%beta/ele%y%beta) 
f1 = f / ele%gamma_c
f2 = 1 / (f * ele%gamma_c)

cc_p%coupling11  = cc_p%cbar(1,1) * f1
cc_p%coupling12a = cc_p%cbar(1,2) * f2
cc_p%coupling12b = cc_p%cbar(1,2) * f1
cc_p%coupling22  = cc_p%cbar(2,2) * f2

cc_p%f_11  = f1
cc_p%f_12a = f2
cc_p%f_12b = f1
cc_p%f_22  = f2

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_macro_data (u, lat, orb, beam, i_ele)
!
! The macroparticle beam is not
! saved at every element due to memory constraints so any data must be
! calculated on the fly during the tracking.
!
! Obviously, data types that require a range doesn't work. That would require a
! bit more work to set up...
!
! Input:
!  u         -- tao_universe_struct
!  lat       -- ring_struct
!  orb       -- coord_struct(0:i_ele): beam centroid coords throguh to the last
!                tracked element
!  beam      -- beam_struct
!  e_ele     -- Integer: index of last tracked element
!
! Output:
!  u%data(:) -- tao_data_struct
!-

subroutine tao_macro_data (u, lat, orb, beam, i_ele)

use macroparticle_mod
use macro_utils_mod

implicit none

type (tao_universe_struct), intent(INOUT) :: u
type (ring_struct), intent(IN) :: lat
type (coord_struct), intent(IN) :: orb(0:)
type (beam_struct), intent(IN) :: beam
integer i_ele

type (tao_data_struct), pointer :: datum
real(rp), pointer :: datum_value
real(rp) dummy

integer i, ii, iii

! Find any datums associated with this element
! This is still not as optimized for speed as it could be
! Better optimized if I assume the element indices are in order
  do i = 1, size(u%beam%macro_data%d2)
    do ii = 1, size(u%beam%macro_data%d2(i)%d1)
      do iii = 1, size(u%beam%macro_data%d2(i)%d1(ii)%d)
												
	if (associated(u%beam%macro_data%d2(i)%d1(ii)%d(iii))) then
	  if (u%beam%macro_data%d2(i)%d1(ii)%d(iii)%ix_ele .eq. i_ele) then
														
	    datum => u%beam%macro_data%d2(i)%d1(ii)%d(iii)
	    datum_value => datum%model_value
	 
	    select case (datum%data_type)
           
	    case ('norm_emittance:x')
	      call calc_bunch_emittance (beam%bunch(1), lat%ele_(i_ele), &
	          	                 datum_value, dummy)
	      datum%s = lat%ele_(i_ele)%s
	 
	    case ('norm_emittance:y')
	      call calc_bunch_emittance (beam%bunch(1), lat%ele_(i_ele), &
	          	                 dummy, datum_value)
	      datum%s = lat%ele_(i_ele)%s
           
	    end select
	  endif
        endif

      enddo
    enddo
  enddo
		
end subroutine tao_macro_data

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_define_macro_data_types ()
!
! Defines which data types are specific to macroparticles
!-

subroutine tao_define_macro_data_types ()

implicit none

integer i

  if (associated(macro_data_names)) then
    do i = 1, size(macro_data_names)
      deallocate(macro_data_names(i)%d1_data)
    enddo
    deallocate(macro_data_names)
  endif
		
  ! total number of macro specifc data types
  allocate(macro_data_names(1))

  ! norm_emittance
  allocate(macro_data_names(1)%d1_data(2))
  macro_data_names(1)%name = 'norm_emittance'
  macro_data_names(1)%d1_data(1) = 'x'
  macro_data_names(1)%d1_data(2) = 'y'

  ! If Fortran was object oriented i could define the data functions here too!
		
end subroutine tao_define_macro_data_types

end module tao_data_mod
