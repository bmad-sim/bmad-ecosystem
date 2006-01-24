module tao_data_mod

use tao_mod
use macroparticle_mod
use macro_utils_mod
use spin_mod

! These are data types specific to macroparticles

type this_coupling_struct
  real(rp) cbar(2,2)
  real(rp) coupling11, coupling12a, coupling12b, coupling22
  real(rp) f_11, f_12a, f_12b, f_22
  logical calc_done
end type

type (this_coupling_struct), save, allocatable, target :: cc(:)

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_data (data_value, data_weight, data_meas_value, dat_ix_dModel)
!
! Subroutine to get the values of the data used in optimization and put them
! in an array. The data is ordered starting with the first universe
!
! Input:
! 
! Output:
!   data_value(:)       -- Real, allocatable, optional: Data model values.
!   data_weight(:)      -- Real, allocatable, optional: Data  weights in the merit function.
!   data_meas_value(:)  -- Real, allocatable, optional: Data values when the data was taken.
!   data_ix_dModel(:)   -- Real, allocatable, optional: Data ix_dModel indices
!-

subroutine tao_get_data (data_value, data_weight, data_meas_value, data_ix_dModel)

implicit none

real(rp), allocatable, optional :: data_value(:), data_meas_value(:), data_weight(:)
real(rp), allocatable, optional :: data_ix_dModel(:)

integer i, j, iu
integer n_data

! 
  
  
n_data = 0
do iu = 1, size(s%u)
  if (.not. s%u(iu)%is_on) cycle
  n_data  = n_data + count(s%u(iu)%data(:)%useit_opt)
enddo
if (present(data_value))      call reallocate_real (data_value, n_data)
if (present(data_meas_value)) call reallocate_real (data_meas_value, n_data)
if (present(data_weight))     call reallocate_real (data_weight, n_data)
if (present(data_ix_dModel))  call reallocate_real (data_ix_dModel, n_data)

j = 0
do iu = 1, size(s%u)
  if (.not. s%u(iu)%is_on) cycle
  do i = 1, size(s%u(iu)%data)
    if (.not. s%u(iu)%data(i)%useit_opt) cycle
    j = j + 1
    if (present(data_value))        data_value(j)      = s%u(iu)%data(i)%model_value
    if (present(data_weight))       data_weight(j)     = s%u(iu)%data(i)%weight
    if (present(data_meas_value))   data_meas_value(j) = s%u(iu)%data(i)%meas_value
    if (present(data_ix_dModel))    data_ix_dModel(j)  = s%u(iu)%data(i)%ix_dModel
  enddo
enddo


end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_load_data_array (uni, ix_ele)
!
! Routine to take data from the model lattice and model orbit
! and put that into the s%u(:)%data(:) arrays.
!
! Input:
!   uni     -- Integer: universe where data resides
!   ix_ele  -- Integer: element to evaluate data at
!-

subroutine tao_load_data_array (u, ix_ele)

implicit none

type (tao_universe_struct) :: u
type (tao_ix_data_struct), pointer :: ix_data
integer uni, ix_ele
integer i
character(20) :: r_name = 'tao_load_data_array'

!

if (ix_ele .eq. 0) call tao_data_coupling_init (u) 
  
! find which datums to evaluate here
if (.not. associated(u%ix_data(ix_ele)%ix_datum)) return

ix_data => u%ix_data(ix_ele)
do i = 1, size(ix_data%ix_datum)
  call tao_evaluate_a_datum (u%data(ix_data%ix_datum(i)), u, u%model, &
              u%model_orb, u%data(ix_data%ix_datum(i))%model_value)
  u%data(ix_data%ix_datum(i))%s = u%model%ele_(ix_ele)%s
enddo
    

end subroutine tao_load_data_array

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
! Subroutine tao_evaluate_a_datum (datum, u, lattice, orb, datum_value, taylor_in)
!
! Subroutine to put the proper data in the specified datum
!
! Input:
!   datum        -- Tao_data_struct: What type of datum
!   u            -- Tao_universe_struct: Which universe to use.
!   lattice      -- Ring_struct: Lattice to use.
!   orb(0:)      -- Coord_struct: Orbit to use.
!   taylor_in(6) -- Taylor_struct, optional: Starting point for 
!                     tt: and t: constraints.
!     
! Output:
!   datum   -- Tao_data_struct: 
!     %ix_ele_merit -- For max/min type constraints: Place where value is max/min. 
!   datum_value -- Real(rp): Value of the datum.
!-

subroutine tao_evaluate_a_datum (datum, u, lattice, orb, datum_value, taylor_in)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (ring_struct) ::  lattice
type (coord_struct) :: orb(0:)
type (modes_struct) mode
type (taylor_struct), save :: taylor(6)
type (taylor_struct), optional :: taylor_in(6)
type (spin_polar_struct) polar

real(rp) datum_value, mat6(6,6), vec0(6)

integer, save :: ix_save = -1
integer i, j, k, m, ix, ix1, ix2, expnt(6)
!integer track_type

character(20) :: r_name = 'tao_evaluate_a_datum'
character(16) data_type

logical found

!

call tao_hook_evaluate_a_datum (found, datum, u, lattice, orb, datum_value)
if (found) return

ix1 = datum%ix_ele
ix2 = datum%ix_ele2
datum_value = 0  ! default

data_type = datum%data_type
if (data_type(1:2) == 'r:') data_type = 'r:'
if (data_type(1:2) == 't:') data_type = 't:'
if (data_type(1:3) == 'tt:') data_type = 'tt:'


select case (data_type)

case ('orbit:x')
  call load_it (orb(:)%vec(1))
case ('orbit:y')
  call load_it (orb(:)%vec(3))
case ('orbit:z')
  call load_it (orb(:)%vec(5))

case ('bpm:x')
  call tao_read_bpm (orb(ix1), lattice%ele_(ix1), x$, datum_value)
case ('bpm:y')
  call tao_read_bpm (orb(ix1), lattice%ele_(ix1), y$, datum_value)

case ('orbit:p_x')
  call load_it (orb(:)%vec(2))
case ('orbit:p_y')
  call load_it (orb(:)%vec(4))
case ('orbit:p_z')
  call load_it (orb(:)%vec(6))

case ('phase:x')
  datum_value = lattice%ele_(ix1)%x%phi - lattice%ele_(ix2)%x%phi
case ('phase:y')
  datum_value = lattice%ele_(ix1)%y%phi - lattice%ele_(ix2)%y%phi

case ('beta:x')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%beta
  endif
    
case ('beta:a')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%beta
  endif
    
case ('beta:y')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%beta
  endif

case ('beta:b')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%beta
  endif

case ('alpha:x')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%alpha
  endif
  
case ('alpha:a')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%alpha
  endif
  
case ('alpha:y')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%alpha
  endif

case ('alpha:b')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%alpha
  endif

case ('eta:x')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%eta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%eta
  else
    call load_it (lattice%ele_(:)%x%eta_lab)
  endif

case ('eta:y')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%eta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%eta
  else
    call load_it (lattice%ele_(:)%y%eta_lab)
  endif

case ('etap:x')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%etap
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%etap
  else
    call load_it (lattice%ele_(:)%x%etap_lab)
  endif

case ('etap:y')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%etap
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%etap
  else
    call load_it (lattice%ele_(:)%y%etap_lab)
  endif

case ('eta:a')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%eta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%eta
  else
    call load_it (lattice%ele_(:)%x%eta)
  endif

case ('eta:b')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%eta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%eta
  else
    call load_it (lattice%ele_(:)%y%eta)
  endif

case ('etap:a')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%etap
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%etap
  else
    call load_it (lattice%ele_(:)%x%etap)
  endif

case ('etap:b')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%etap
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%etap
  else
    call load_it (lattice%ele_(:)%y%etap)
  endif

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
  if (ix1 < ix2) return
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return
  call transfer_matrix_calc (lattice, .true., mat6, vec0, ix2, ix1)
  datum_value = mat6(i, j)

case ('t:')
  if (ix1 < ix2) return
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return
  k = tao_read_this_index (datum%data_type, 5); if (k == 0) return
  if (present(taylor_in)) then
    call tao_transfer_map_calc (lattice, taylor_in, ix2, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), j, k)
  else
    call tao_transfer_map_calc (lattice, taylor, ix2, ix1)
    datum_value = taylor_coef (taylor(i), j, k)
  endif

case ('tt:')
  if (ix1 < ix2) return
  expnt = 0
  i = tao_read_this_index (datum%data_type, 4); if (i == 0) return
  do j = 5, 15
    if (datum%data_type(j:j) == ' ') exit
    k = tao_read_this_index (datum%data_type, j); if (k == 0) return
    expnt(k) = expnt(k) + 1
  enddo
  if (present(taylor_in)) then
    call tao_transfer_map_calc (lattice, taylor_in, ix2, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), expnt)
  else
    call tao_transfer_map_calc (lattice, taylor, ix2, ix1)
    datum_value = taylor_coef (taylor(i), expnt)
  endif

case ('floor:x')
  if (ix2 >= 0) then
    datum_value = lattice%ele_(ix1)%floor%x - lattice%ele_(ix2)%floor%x
  else
    datum_value = lattice%ele_(ix1)%floor%x
  endif

case ('floor:y')
  if (ix2 >= 0) then
    datum_value = lattice%ele_(ix1)%floor%y - lattice%ele_(ix2)%floor%y
  else
    datum_value = lattice%ele_(ix1)%floor%y 
  endif

case ('floor:z')
  if (ix2 >= 0) then
    datum_value = lattice%ele_(ix1)%floor%z - lattice%ele_(ix2)%floor%z
  else
    datum_value = lattice%ele_(ix1)%floor%z 
  endif

case ('floor:theta')
  if (ix2 >= 0) then
    datum_value = lattice%ele_(ix1)%floor%theta - lattice%ele_(ix2)%floor%theta
  else
    datum_value = lattice%ele_(ix1)%floor%theta 
  endif

case ('s_position') 
  if (ix2 >= 0) then
    datum_value = lattice%ele_(ix1)%s - lattice%ele_(ix2)%s
  else
    datum_value = lattice%ele_(ix1)%s 
  endif

case ('norm_emittance:x')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%norm_emitt
  else
    datum_value = 0.0
  endif
  
case ('norm_emittance:y')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%norm_emitt
  else
    datum_value = 0.0
  endif
  
case ('norm_emittance:z')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%z%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif
  
case ('norm_emittance:a')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%norm_emitt
  else
    call orbit_amplitude_calc (lattice%ele_(ix1), orb(ix1), &
                               amp_na = datum_value, particle = electron$)
  endif
  
case ('norm_emittance:b')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%norm_emitt
  else
    call orbit_amplitude_calc (lattice%ele_(ix1), orb(ix1), &
                               amp_nb = datum_value, particle = electron$)
  endif
  
case ('dpx_dx') 
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%dpx_dx
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpy_dy') 
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%dpx_dx
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpz_dz') 
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%z%dpx_dx
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%z%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpa_da') 
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%dpx_dx
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpb_db') 
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%dpx_dx
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%dpx_dx
  else
    datum_value = 0.0
  endif

case ('sigma:x')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma:p_x')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%p_sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma:y')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma:p_y')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%p_sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma:z')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%z%sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%z%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma:p_z')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%z%p_sigma
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%z%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('spin:theta')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%spin%theta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = 0.0
  else
    call spinor_to_polar (orb(ix1), polar)
    datum_value = polar%theta
  endif
  
case ('spin:phi')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%spin%phi
  elseif (s%global%track_type .eq. "macro") then
    datum_value = 0.0
  else
    call spinor_to_polar (orb(ix1), polar)
    datum_value = polar%phi
  endif
  
case ('spin:polarity')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%spin%polarization
  elseif (s%global%track_type .eq. "macro") then
    datum_value = 0.0
  else
    datum_value = 1.0
  endif
  
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

if (datum%ele2_name == ' ') then
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

end subroutine coupling_calc

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_read_bpm (orb, ele, axis, reading)
!
! Find the reading on a bpm given the orbit at the BPM and the BPM offsets
!
! This routine will only give a nonzero reading for BMAD monitors and
! instruments.
!
! The BPM noise is in the ele%r(1,1) entry
!
! Input: 
!  orb      -- Coord_struct: Orbit position at BPM
!  ele      -- Ele_struct: the BPM
!  noise    -- real(rp): BPM gaussian noise (in meters)
!  axis     -- Integer: x$ or y$
!
! Output:
!  reading  -- Real(rp): BPM reading
!
!-

subroutine tao_read_bpm (orb, ele, axis, reading)

use random_mod

type (coord_struct) orb
type (ele_struct) ele

real(rp) reading
real(rp) ran_num

integer axis

character(20) :: r_name = "tao_read_bpm"

logical err

  if (.not. ele%is_on) then
    reading = 0.0
    return
  endif

  if (ele%key .ne. monitor$ .and. ele%key .ne. instrument$) then
    reading = 0.0
    return
  endif

  call ran_gauss (ran_num)
  
  if (axis .eq. x$) then
    reading = ele%r(1,1)*ran_num + &
               (orb%vec(1) - ele%value(x_offset_tot$) + ele%r(1,2)) * &
                                           cos(ele%value(tilt_tot$) + ele%r(1,4)) + &
               (orb%vec(3) - ele%value(y_offset_tot$) + ele%r(1,3)) * &
                                           sin(ele%value(tilt_tot$) + ele%r(1,4))              
  elseif (axis .eq. y$) then
    reading = ele%r(1,1)*ran_num &
              -(orb%vec(1) - ele%value(x_offset_tot$) + ele%r(1,2)) * &
                                           sin(ele%value(tilt_tot$) + ele%r(1,4)) + &
               (orb%vec(3) - ele%value(y_offset_tot$) + ele%r(1,3)) * &
                                           cos(ele%value(tilt_tot$) + ele%r(1,4))  
  else
    reading = 0.0
    call out_io (s_warn$, r_name, &
                 "This axis not supported for BPM reading!")
  endif
    
end subroutine tao_read_bpm

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_transfer_map_calc (lat, t_map, ix1, ix2, &
!                                         integrate, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between two elements.
!
! The transfer map is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn map is calculated.
! If ix2 < ix1 then the calculation will "wrap around" the lattice end.
! For example if ix1 = 900 and ix2 = 10 then the t_mat is the map from
! element 900 to the lattice end plus from 0 through 10. 
!
! If ix2 = ix1 then you get the unit map except if one_turn = True.
!
! Note: If integrate = False and if a taylor map does not exist for an 
! element this routine will make one and store it in the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Ring_struct: Lattice used in the calculation.
!   t_map(6)   -- Taylor_struct: Initial map (used when unit_start = False)
!   ix1        -- Integer, optional: Element start index for the calculation.
!                   Default is 0.
!   ix2        -- Integer, optional: Element end index for the calculation.
!                   Default is lat%n_ele_use.
!   integrate  -- Logical, optional: If present and True then do symplectic
!                   integration instead of concatenation. 
!                   Default = False.
!   one_turn   -- Logical, optional: If present and True then construct the
!                   one-turn map from ix1 back to ix1 (ignoring ix2).
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

subroutine tao_transfer_map_calc (lat, t_map, ix1, ix2, &
                                      integrate, one_turn, unit_start)

  use bmad_struct
  use bmad_interface
  use ptc_interface_mod, only: concat_taylor, ele_to_taylor, taylor_propagate1

  implicit none

  type (ring_struct) lat

  type (taylor_struct) :: t_map(:), taylor2(6)

  integer, intent(in), optional :: ix1, ix2
  integer i, i1, i2

  logical, optional :: integrate, one_turn, unit_start
  logical integrate_this, one_turn_this

!

  integrate_this = logic_option (.false., integrate)
  one_turn_this = logic_option (.false., one_turn)

  i1 = integer_option(0, ix1) 
  i2 = integer_option(lat%n_ele_use, ix2)
  if (one_turn_this) i2 = i1
 
  if (logic_option(.true., unit_start)) call taylor_make_unit (t_map)


  if (i2 < i1 .or. one_turn_this) then
    do i = i1+1, lat%n_ele_use
      call add_on_to_t_map
    enddo
    do i = 1, i2
      call add_on_to_t_map
    enddo

  else
    do i = i1+1, i2
      call add_on_to_t_map
    enddo
  endif

!--------------------------------------------------------
contains

subroutine add_on_to_t_map

  if (integrate_this) then
    call taylor_propagate1 (t_map, lat%ele_(i), lat%param)
  else
    if (.not. associated(lat%ele_(i)%taylor(1)%term)) then
      call ele_to_taylor (lat%ele_(i), lat%param)
    endif

    call concat_taylor (t_map, lat%ele_(i)%taylor, t_map)
  endif

end subroutine

end subroutine
                                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_transfer_map_calc_at_s (lat, t_map, s1, s2, &
!                                         integrate, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between two elements.
!
! The transfer map is from longitudinal position s1 to s2.
! If s1 and s2 are not present, the full 1-turn map is calculated.
! If s2 < s1 then the calculation will "wrap around" the lattice end.
! For example if s1 = 900 and s2 = 10 then the t_mat is the map from
! s = 900 to the lattice end plus from 0 through s = 10. 
!
! If s2 = s1 then you get the unit map except if one_turn = True.
!
! Note: If integrate = False and if a taylor map does not exist for an 
! element this routine will make one and store it in the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Ring_struct: Lattice used in the calculation.
!   t_map(6)   -- Taylor_struct: Initial map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start index for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end index for the calculation.
!                   Default is lat%param%total_length.
!   integrate  -- Logical, optional: If present and True then do symplectic
!                   integration instead of concatenation. 
!                   Default = False.
!   one_turn   -- Logical, optional: If present and True then construct the
!                   one-turn map from s1 back to s1 (ignoring s2).
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

subroutine tao_transfer_map_calc_at_s (lat, t_map, s1, s2, &
                                      integrate, one_turn, unit_start)

  use bmad_struct
  use bmad_interface
  use ptc_interface_mod, only: concat_taylor, ele_to_taylor, taylor_propagate1

  implicit none

  type (ring_struct) lat
  type (taylor_struct) :: t_map(:)

  real(rp), intent(in), optional :: s1, s2
  real(rp) ss1, ss2

  logical, optional :: integrate, one_turn, unit_start
  logical integrate_this, one_turn_this

!

  integrate_this = logic_option (.false., integrate)
  one_turn_this = logic_option (.false., one_turn)

  ss1 = 0;                       if (present(s1)) ss1 = s1
  ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
  if (one_turn_this) ss2 = ss1
 
  if (logic_option(.true., unit_start)) call taylor_make_unit (t_map)


  if (ss2 < ss1 .or. one_turn_this) then
    call transfer_this (ss1, lat%param%total_length)
    call transfer_this (0.0_rp, ss2)
  else
    call transfer_this (ss1, ss2)
  endif

!--------------------------------------------------------
! Known problems:
!   1) map type wigglers not treated properly.
!   2) need to reuse the taylor map? (is time really an issue?)

contains

subroutine transfer_this (s_1, s_2)

  type (ele_struct), save :: ele

  real(rp) s_1, s_2, s_now, s_end, ds
  real(rp), save :: ds_old = -1

  integer ix_ele
  integer, save :: ix_ele_old = -1

  logical kill_it

!

  call tao_ele_at_s (lat, s_1, ix_ele)

  if (ix_ele /= ix_ele_old) ele = lat%ele_(ix_ele)
  s_now = s_1
  kill_it = .false.

  do
    s_end = min(s_2, ele%s)
    if (ele%key == sbend$) then
      if (s_now /= lat%ele_(ix_ele-1)%s) ele%value(e1$) = 0
      if (s_end /= ele%s) ele%value(e2$) = 0
      if (s_now == lat%ele_(ix_ele-1)%s) kill_it = .true.
      if (s_end == ele%s) kill_it = .true.
    endif

    ds = s_end - s_now
    ele%value(l$) = ds
    if (lat%ele_(ix_ele)%value(l$) /= 0) ele%num_steps = &
                  1 + 0.999 * lat%ele_(ix_ele)%num_steps * ds / lat%ele_(ix_ele)%value(l$)

    if (ds /= ds_old .or. ix_ele /= ix_ele_old) kill_it = .true.

    if (kill_it) call kill_taylor (ele%taylor)

    if (integrate_this) then
      call taylor_propagate1 (t_map, ele, lat%param)
    else
      if (.not. associated(ele%taylor(1)%term)) then
        call ele_to_taylor (ele, lat%param)
      endif

      call concat_taylor (t_map, ele%taylor, t_map)
    endif

    if (s_end == s_2) then
      ix_ele_old = ix_ele
      ds_old = ds
      return
    endif

    s_now = s_end
    ix_ele = ix_ele + 1
    ele = lat%ele_(ix_ele)
  enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_mat6_calc_at_s (lat, mat6, s1, s2, one_turn, unit_start)
!
! Subroutine to calculate the transfer matrix between two elements.
!
! The transfer matrix is from longitudinal position s1 to s2.
! If s1 and s2 are not present, the full 1-turn map is calculated.
! If s2 < s1 then the calculation will "wrap around" the lattice end.
! For example if s1 = 900 and s2 = 10 then the t_mat is the map from
! s = 900 to the lattice end plus from 0 through s = 10. 
!
! If s2 = s1 then you get the unit matrix except if one_turn = True.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Ring_struct: Lattice used in the calculation.
!   mat6(6,6)  -- Taylor_struct: Initial map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start index for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end index for the calculation.
!                   Default is lat%param%total_length.
!   one_turn   -- Logical, optional: If present and True then construct the
!                   one-turn map from s1 back to s1 (ignoring s2).
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

subroutine tao_mat6_calc_at_s (lat, mat6, s1, s2, one_turn, unit_start)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) lat

  real(rp) mat6(:,:)
  real(rp), intent(in), optional :: s1, s2
  real(rp) ss1, ss2

  logical, optional :: one_turn, unit_start
  logical one_turn_this

!

  one_turn_this = logic_option (.false., one_turn)

  ss1 = 0;                       if (present(s1)) ss1 = s1
  ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
  if (one_turn_this) ss2 = ss1
 
  if (logic_option(.true., unit_start)) call mat_make_unit (mat6)


  if (ss2 < ss1 .or. one_turn_this) then
    call transfer_this (ss1, lat%param%total_length)
    call transfer_this (0.0_rp, ss2)
  else
    call transfer_this (ss1, ss2)
  endif

!--------------------------------------------------------
! Known problems:
!   1) map type wigglers not treated properly.
!   2) need to reuse mat6? (is time really an issue?)

contains

subroutine transfer_this (s_1, s_2)

  type (ele_struct), save :: ele
  real(rp) s_1, s_2, s_end, s_now, ds
  integer ix_ele

!

  call tao_ele_at_s (lat, s_1, ix_ele)
  ele = lat%ele_(ix_ele)
  s_now = s_1

  do
    s_end = min(s_2, ele%s)
    ds = s_end - s_now
    ele%value(l$) = ds
    if (ele%key == sbend$) then
      if (s_now /= lat%ele_(ix_ele-1)%s) ele%value(e1$) = 0
      if (s_end /= ele%s) ele%value(e2$) = 0
    endif

    call make_mat6 (ele, lat%param)
    mat6 = matmul (ele%mat6, mat6)

    if (s_end == s_2) return
    s_now = s_end
    ix_ele = ix_ele + 1
    ele = lat%ele_(ix_ele)
  enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_ele_at_s (lat, s, ix_ele)
!
! Subroutine to return the index of the element at position s.
! That is, ix_ele is choisen such that:
!     lat%ele_(ix_ele-1)%s < s <= lat%ele_(ix_ele)%s
!
! Note: s is evaluated modulo the lattice length:
!     s -> s - lat_length * floor(s/lat_length)
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- Ring_struct: Lattice of elements.
!   s   -- Real(rp): Longitudinal position.
!
! Output:
!   ix_ele -- Integer: Index of element at s.
!-

subroutine tao_ele_at_s (lat, s, ix_ele)

  use bmad

  implicit none

  type (ring_struct) lat
  real(rp) s, ss, ll
  integer ix_ele, n1, n2, n3

!

  ll = lat%param%total_length
  ss = s - ll * floor(s/ll)

  n1 = 0
  n3 = lat%n_ele_use

  do

    if (n3 == n1 + 1) then
      ix_ele = n3
      return
    endif

    n2 = (n1 + n3) / 2

    if (lat%ele_(n2)%s < ss) then
      n1 = n2
    else
      n3 = n2
    endif

  enddo

end subroutine


end module tao_data_mod
