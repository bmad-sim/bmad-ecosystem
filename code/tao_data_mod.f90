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

implicit none

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

implicit none

type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: data
type (tao_data_struct) datum
type (ring_struct) lattice
type (coord_struct) orb(0:)
type (modes_struct) mode
type (taylor_struct), save :: taylor(6)

real(rp) datum_value, mat6(6,6)

integer i, j, k, m, ix, ix1, ix2

character(20) :: r_name = 'tao_evaluate_a_datum'
character(16) data_type

logical found

!

call tao_hook_load_data_array (found, datum, lattice, orb, datum_value)
if (found) return

ix1 = datum%ix_ele
ix2 = datum%ix_ele2
ele => lattice%ele_(ix1)

data_type = datum%data_type
if (data_type(1:2) == 'r:') data_type = 'r:'
if (data_type(1:2) == 't:') data_type = 't:'

select case (data_type)

case ('orbit:x')
  call load_it (orb(:)%vec(1))
case ('orbit:y')
  call load_it (orb(:)%vec(3))
case ('orbit:z')
  call load_it (orb(:)%vec(5))

case ('orbit:p_x')
  call load_it (orb(:)%vec(2))
case ('orbit:p_y')
  call load_it (orb(:)%vec(4))
case ('orbit:p_z')
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

case ('r:')
  call relative_switch
  i = read_this_index (datum%data_type(3:3)); if (i == 0) return
  j = read_this_index (datum%data_type(4:4)); if (j == 0) return
  call transfer_matrix_calc (lattice, .true., mat6, ix1, ix2)
  datum_value = mat6(i, j)

case ('t:')
  call relative_switch
  i = read_this_index (datum%data_type(3:3)); if (i == 0) return
  j = read_this_index (datum%data_type(4:4)); if (j == 0) return
  k = read_this_index (datum%data_type(5:5)); if (k == 0) return
  call transfer_map_calc (lattice, taylor, ix1, ix2)
  datum_value = taylor_coef (taylor(i), j, k)

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

! dummy places for macroparticle data
case ('norm_emittance:x')
case ('norm_emittance:y')  
case ('dpz_dz')  
case ('bunch_sigma:x')  
case ('bunch_sigma:p_x')  
case ('bunch_sigma:y')  
case ('bunch_sigma:p_y')  
case ('bunch_sigma:z')  
case ('bunch_sigma:p_z')  
case ('bunch_alpha:x') 
case ('bunch_alpha:y') 
case ('bunch_beta:x') 
case ('bunch_beta:y') 
  
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

function read_this_index (char) result (ix)

  character(1) char
  integer ix

  ix = index('123456', char)
  if (ix == 0) then
    call out_io (s_abort$, r_name, 'BAD INDEX CONSTRAINT: ' // data%data_type)
    call err_exit
  endif

end function

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
type (bunch_params_struct) :: params

type (tao_data_struct), pointer :: datum
real(rp), pointer :: datum_value
real(rp) dummy, bunch(6)

integer i, ii, iii, i_ele
integer, save :: i_save

logical do_calc_bunch_params

character(20) :: r_name = "tao_macro_data"

  ! reset index counters
  if (i_ele .eq. 0) then
    do i = 1 , size(u%beam%macro_data%d2)
      u%beam%macro_data%d2(i)%d1%i_save = 1
    enddo
  endif
    

! Only calculate the bunch parameters once then fill in the datums 
! as they come up
  do_calc_bunch_params = .true.

! Find any datums associated with this element
! assume the element indices are in order and use the i_save index so that we
! only pass through the datums once per lattice calculation.
  d2_loop: do i = 1, size(u%beam%macro_data%d2)
    d1_loop: do ii = 1, size(u%beam%macro_data%d2(i)%d1)
      d_loop: do iii = u%beam%macro_data%d2(i)%d1(ii)%i_save, size(u%beam%macro_data%d2(i)%d1(ii)%d)
												
	if (associated(u%beam%macro_data%d2(i)%d1(ii)%d)) then
	  if (u%beam%macro_data%d2(i)%d1(ii)%d(iii)%ix_ele .eq. i_ele) then
            u%beam%macro_data%d2(i)%d1(ii)%i_save = iii
					    
	    if (do_calc_bunch_params) then
	      call calc_bunch_params (beam%bunch(s%global%bunch_to_plot), &
	                              lat%ele_(i_ele), params)
	      do_calc_bunch_params = .false.
	    endif

	    datum => u%beam%macro_data%d2(i)%d1(ii)%d(iii)
	    datum_value => datum%model_value
	 
	    select case (datum%data_type)
           
	      case ('norm_emittance:x')
	        datum_value = params%x%emitt
	     
	      case ('norm_emittance:y')
	        datum_value = params%y%emitt
             
	      case ('dpz_dz')
	        datum_value = params%dpz_dz
            
              case ('bunch_beta:x')
	        datum_value = params%x%beta
	        
              case ('bunch_beta:y')
	        datum_value = params%y%beta
	        
              case ('bunch_alpha:x')
	        datum_value = params%x%alpha
	        
	      case ('bunch_alpha:y')
	        datum_value = params%y%alpha
	        
	      case ('bunch_sigma:x')
	        datum_value = params%x%sigma
	        
	      case ('bunch_sigma:p_x')
	        datum_value = params%x%p_sigma
	        
	      case ('bunch_sigma:y')
	        datum_value = params%y%sigma
	        
	      case ('bunch_sigma:p_y')
	        datum_value = params%y%p_sigma
	        
	      case ('bunch_sigma:z')
	        datum_value = params%z_sigma
	        
	      case ('bunch_sigma:p_z')
	        datum_value = params%p_z_sigma
	      
	      case default
		call out_io (s_error$, r_name, &
		         "macro data type not found")
	    end select

	    datum%s = lat%ele_(i_ele)%s
	    cycle d1_loop

	  endif
        endif

      enddo d_loop
    enddo d1_loop
  enddo d2_loop
		
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

integer, parameter :: n_macro_data_types = 5

  if (associated(macro_data_names)) then
    do i = 1, size(macro_data_names)
      deallocate(macro_data_names(i)%d1_data)
    enddo
    deallocate(macro_data_names)
  endif
		
  ! total number of macro specifc data types
  allocate(macro_data_names(n_macro_data_types))

  ! norm_emittance
  allocate(macro_data_names(1)%d1_data(2))
  macro_data_names(1)%name = 'norm_emittance'
  macro_data_names(1)%d1_data(1) = 'x'
  macro_data_names(1)%d1_data(2) = 'y'

  ! dPz_dZ
  allocate(macro_data_names(2)%d1_data(1))
  macro_data_names(2)%name = 'dpz_dz'
  macro_data_names(2)%d1_data(1) = 'dpz_dz'

  ! Bunch Size
  allocate(macro_data_names(3)%d1_data(6))
  macro_data_names(3)%name = 'bunch_sigma'
  macro_data_names(3)%d1_data(1) = 'x'
  macro_data_names(3)%d1_data(2) = 'p_x'
  macro_data_names(3)%d1_data(3) = 'y'
  macro_data_names(3)%d1_data(4) = 'p_y'
  macro_data_names(3)%d1_data(5) = 'z'
  macro_data_names(3)%d1_data(6) = 'p_z'

  ! Beam beta
  allocate(macro_data_names(4)%d1_data(3))
  macro_data_names(4)%name = 'bunch_beta'
  macro_data_names(4)%d1_data(1) = 'x'
  macro_data_names(4)%d1_data(2) = 'y'
  macro_data_names(4)%d1_data(3) = 'z'

  ! Beam beta
  allocate(macro_data_names(5)%d1_data(3))
  macro_data_names(5)%name = 'bunch_alpha'
  macro_data_names(5)%d1_data(1) = 'x'
  macro_data_names(5)%d1_data(2) = 'y'
  macro_data_names(5)%d1_data(3) = 'z'

  ! If Fortran was object oriented i could define the data functions here too!
		
end subroutine tao_define_macro_data_types


end module tao_data_mod
