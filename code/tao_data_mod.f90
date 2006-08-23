module tao_data_mod

use tao_mod
use macroparticle_mod
use macro_utils_mod
use spin_mod
use utilities_mod
use random_mod


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
!   data_weight(:)      -- Real, allocatable, optional: Data weights in the merit function.
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
! Subroutine tao_load_data_array (uni, ix_ele, track_type)
!
! Routine to take data from the model lattice and model orbit
! and put that into the s%u(:)%data(:) arrays.
!
! Input:
!   uni         -- Integer: universe where data resides
!   ix_ele      -- Integer: element to evaluate data at
!   track_type  -- Character(*): Type of particle tracking.
!-

subroutine tao_load_data_array (u, ix_ele, track_type)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum

integer, pointer :: ix_datum(:)
integer uni, ix_ele
integer i

character(*) track_type
character(20) :: r_name = 'tao_load_data_array'

!

if (ix_ele == 0) call tao_data_coupling_init (u) 
  
! find which datums to evaluate here
if (.not. associated(u%ix_data(ix_ele)%ix_datum)) return

ix_datum => u%ix_data(ix_ele)%ix_datum
do i = 1, size(ix_datum)
  datum => u%data(ix_datum(i))
  call tao_evaluate_a_datum (datum, u, u%model, track_type, datum%model_value)
  if (datum%ix_ele_merit > -1) datum%s = &
                                    u%model%lat%ele_(datum%ix_ele_merit)%s
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

m = u%model%lat%n_ele_max
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
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, track_type, &
!                                                 datum_value, taylor_in)
!
! Subroutine to put the proper data in the specified datum
!
! If datum results in NaN then datum_value = tiny(1.0_rp)
!
! Input:
!   datum        -- Tao_data_struct: What type of datum
!   u            -- Tao_universe_struct: Which universe to use.
!   tao_lat      -- Tao_lattice_struct: Lattice to use.
!   track_type   -- Character(*): Type of particle tracking.
!   taylor_in(6) -- Taylor_struct, optional: Starting point for 
!                     tt: and t: constraints.
!     
! Output:
!   datum   -- Tao_data_struct: 
!     %ix_ele_merit -- For max/min type constraints: Place where value is max/min. 
!   datum_value -- Real(rp): Value of the datum.
!-

subroutine tao_evaluate_a_datum (datum, u, tao_lat, track_type, &
                                                  datum_value, taylor_in)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat
type (ring_struct), pointer :: lat
type (modes_struct) mode
type (taylor_struct), save :: taylor(6)
type (taylor_struct), optional :: taylor_in(6)
type (spin_polar_struct) polar
type (ele_struct), pointer :: ele0
type (coord_struct), pointer :: orb0

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), one_pz
real(rp) gamma

integer, save :: ix_save = -1
integer i, j, k, m, n, ix, ix1, ix0, expnt(6), n_lat

character(*) track_type
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) data_type

logical found

!

call tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value)
if (found) return

ix0 = datum%ix_ele0
ix1 = datum%ix_ele
n_lat = tao_lat%lat%n_ele_use
datum_value = 0           ! default
datum%ix_ele_merit = ix1  ! default
lat => tao_lat%lat

data_type = datum%data_type
if (data_type(1:2) == 'r.') data_type = 'r.'
if (data_type(1:2) == 't.') data_type = 't.'
if (data_type(1:3) == 'tt.') data_type = 'tt.'
if (data_type(1:5) == 'wire.') data_type = 'wire.'

if (data_type == 'r.' .or. data_type == 't.' .or. data_type == 'tt.') then
  if (ix0 > ix1 .and. lat%param%lattice_type == linear_lattice$) then
    call out_io (s_error$, r_name, &
            'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
            'STARTING ELEMENT: ' // lat%ele_(ix0)%name, &
            'IS AFTER ENDING ELEMENT: ' // lat%ele_(ix1)%name)
    return
  endif
endif

!---------------------------------------------------

select case (data_type)

case ('orbit.x')
  call load_it (tao_lat%orb(:)%vec(1), ix0, ix1, datum_value, datum, lat)
case ('orbit.y')
  call load_it (tao_lat%orb(:)%vec(3), ix0, ix1, datum_value, datum, lat)
case ('orbit.z')
  call load_it (tao_lat%orb(:)%vec(5), ix0, ix1, datum_value, datum, lat)

case ('bpm.x')
  call tao_read_bpm (tao_lat%orb(ix1), lat%ele_(ix1), x_plane$, datum_value)
case ('bpm.y')
  call tao_read_bpm (tao_lat%orb(ix1), lat%ele_(ix1), y_plane$, datum_value)

case ('orbit.p_x')
  call load_it (tao_lat%orb(:)%vec(2), ix0, ix1, datum_value, datum, lat)
case ('orbit.p_y')
  call load_it (tao_lat%orb(:)%vec(4), ix0, ix1, datum_value, datum, lat)
case ('orbit.p_z')
  call load_it (tao_lat%orb(:)%vec(6), ix0, ix1, datum_value, datum, lat)

case ('phase.x')
  datum_value = lat%ele_(ix1)%x%phi - lat%ele_(ix0)%x%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele_(0)%x%phi + lat%ele_(n_lat)%x%phi 
case ('phase.y')
  datum_value = lat%ele_(ix1)%y%phi - lat%ele_(ix0)%y%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele_(0)%y%phi + lat%ele_(n_lat)%y%phi 

case ('phase_frac_diff')
  px = lat%ele_(ix1)%x%phi - lat%ele_(ix0)%x%phi
  if (ix0 > ix1) px = px - lat%ele_(0)%x%phi + lat%ele_(n_lat)%x%phi 
  py = lat%ele_(ix1)%y%phi - lat%ele_(ix0)%y%phi
  if (ix0 > ix1) py = py - lat%ele_(0)%y%phi + lat%ele_(n_lat)%y%phi 
  datum_value = modulo (px, twopi) - modulo (py, twopi)

case ('beta.x')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%x%beta, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%beta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%beta
  endif
    
case ('beta.y')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%y%beta, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%beta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%beta
  endif

case ('beta.z')
  if (track_type == "single") then
    datum_value = 0.0
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%beta
  elseif (track_type == "macro") then
    datum_value = 0.0
  endif

case ('beta.a')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%x%beta, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%beta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%beta
  endif
    
case ('beta.b')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%y%beta, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%beta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%beta
  endif

case ('beta.c')
  if (track_type == "single") then
    datum_value = 0.0
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%beta
  elseif (track_type == "macro") then
    datum_value = 0.0
  endif

case ('alpha.x')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%x%alpha, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%alpha
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%alpha
  endif
  
case ('alpha.y')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%y%alpha, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%alpha
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%alpha
  endif

case ('alpha.z')
  if (track_type == "single") then
    datum_value = 0.0
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%alpha
  elseif (track_type == "macro") then
    datum_value = 0.0
  endif

case ('alpha.a')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%x%alpha, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%alpha
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%alpha
  endif
  
case ('alpha.b')
  if (track_type == "single") then
    call load_it (lat%ele_(:)%y%alpha, ix0, ix1, datum_value, datum, lat)
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%alpha
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%alpha
  endif

case ('alpha.c')
  if (track_type == "single") then
    datum_value = 0.0
  elseif (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%alpha
  elseif (track_type == "macro") then
    datum_value = 0.0
  endif

case ('eta.x')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%eta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%eta
  else
    call load_it (lat%ele_(:)%x%eta_lab, ix0, ix1, datum_value, datum, lat)
  endif

case ('eta.y')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%eta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%eta
  else
    call load_it (lat%ele_(:)%y%eta_lab, ix0, ix1, datum_value, datum, lat)
  endif

case ('eta.z')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%eta
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif

case ('etap.x')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%etap
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%etap
  else
    call load_it (lat%ele_(:)%x%etap_lab, ix0, ix1, datum_value, datum, lat)
  endif

case ('etap.y')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%etap
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%etap
  else
    call load_it (lat%ele_(:)%y%etap_lab, ix0, ix1, datum_value, datum, lat)
  endif

case ('etap.z')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%etap
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif

case ('eta.a')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%eta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%eta
  else
    call load_it (lat%ele_(:)%x%eta, ix0, ix1, datum_value, datum, lat)
  endif

case ('eta.b')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%eta
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%eta
  else
    call load_it (lat%ele_(:)%y%eta, ix0, ix1, datum_value, datum, lat)
  endif

case ('eta.c')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%eta
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif

case ('etap.a')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%etap
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%etap
  else
    call load_it (lat%ele_(:)%x%etap, ix0, ix1, datum_value, datum, lat)
  endif

case ('etap.b')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%etap
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%etap
  else
    call load_it (lat%ele_(:)%y%etap, ix0, ix1, datum_value, datum, lat)
  endif

case ('etap.c')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%etap
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif

case ('beam_energy')
  call load_it (lat%ele_(0:n_lat)%value(beam_energy$) * &
                       (1+tao_lat%orb(0:n_lat)%vec(6)), ix0, ix1, datum_value, datum, lat)

case ('%beam_energy')
  call load_it (tao_lat%orb(0:n_lat)%vec(6), ix0, ix1, datum_value, datum, lat)
  
case ('coupling.11b')
  call load_it (cc%coupling11, ix0, ix1, datum_value, datum, lat, cc%f_11, coupling_here = .true.)
case ('coupling.12a')
  call load_it (cc%coupling12a, ix0, ix1, datum_value, datum, lat, cc%f_12a, coupling_here = .true.)
case ('coupling.12b')
  call load_it (cc%coupling12b, ix0, ix1, datum_value, datum, lat, cc%f_12b, coupling_here = .true.)
case ('coupling.22a')
  call load_it (cc%coupling22, ix0, ix1, datum_value, datum, lat, cc%f_22, coupling_here = .true.)

case ('cbar.11')
  call load_it (cc%cbar(1,1), ix0, ix1, datum_value, datum, lat, coupling_here = .true.)
case ('cbar.12')
  call load_it (cc%cbar(1,2), ix0, ix1, datum_value, datum, lat, coupling_here = .true.)
case ('cbar.21')
  call load_it (cc%cbar(2,1), ix0, ix1, datum_value, datum, lat, coupling_here = .true.)
case ('cbar.22')
  call load_it (cc%cbar(2,2), ix0, ix1, datum_value, datum, lat, coupling_here = .true.)

case ('i5a_e6')
  if (ix0 > 0 .or. ix1 > 0) then
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_use
    datum_value = sum(tao_lat%rad_int%lin_i5a_e6(ix0:ix1))
    datum%ix_ele_merit = ix1
  else
    datum_value = tao_lat%modes%lin%i5a_e6
    datum%ix_ele_merit = lat%n_ele_use
  endif

case ('i5b_e6')
  if (ix0 > 0 .or. ix1 > 0) then
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_use
    datum_value = sum(tao_lat%rad_int%lin_i5b_e6(ix0:ix1))
    datum%ix_ele_merit = ix1
  else
    datum_value = tao_lat%modes%lin%i5b_e6
    datum%ix_ele_merit = lat%n_ele_use
  endif

case ('r.')
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix0, ix1)
  datum_value = mat6(i, j)

case ('t.')
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return
  k = tao_read_this_index (datum%data_type, 5); if (k == 0) return
  if (present(taylor_in)) then
    call tao_transfer_map_calc (lat, taylor_in, ix0, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), j, k)
  else
    call tao_transfer_map_calc (lat, taylor, ix0, ix1)
    datum_value = taylor_coef (taylor(i), j, k)
  endif

case ('tt.')
  expnt = 0
  i = tao_read_this_index (datum%data_type, 4); if (i == 0) return
  do j = 5, 15
    if (datum%data_type(j:j) == ' ') exit
    k = tao_read_this_index (datum%data_type, j); if (k == 0) return
    expnt(k) = expnt(k) + 1
  enddo
  if (present(taylor_in)) then
    call tao_transfer_map_calc (lat, taylor_in, ix0, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), expnt)
  else
    call tao_transfer_map_calc (lat, taylor, ix0, ix1)
    datum_value = taylor_coef (taylor(i), expnt)
  endif

case ('floor.x')
  if (ix0 >= 0) then
    datum_value = lat%ele_(ix1)%floor%x - lat%ele_(ix0)%floor%x
  else
    datum_value = lat%ele_(ix1)%floor%x
  endif

case ('floor.y')
  if (ix0 >= 0) then
    datum_value = lat%ele_(ix1)%floor%y - lat%ele_(ix0)%floor%y
  else
    datum_value = lat%ele_(ix1)%floor%y 
  endif

case ('floor.z')
  if (ix0 >= 0) then
    datum_value = lat%ele_(ix1)%floor%z - lat%ele_(ix0)%floor%z
  else
    datum_value = lat%ele_(ix1)%floor%z 
  endif

case ('floor.theta')
  if (ix0 >= 0) then
    datum_value = lat%ele_(ix1)%floor%theta - lat%ele_(ix0)%floor%theta
  else
    datum_value = lat%ele_(ix1)%floor%theta 
  endif

case ('s_position') 
  if (ix0 >= 0) then
    datum_value = lat%ele_(ix1)%s - lat%ele_(ix0)%s
  else
    datum_value = lat%ele_(ix1)%s 
  endif

!---------------------------------------------------------
! Beam Emittance
  
case ('beam_emittance.x')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%norm_emitt
  else
    datum_value = 0.0
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma
  
case ('beam_emittance.y')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%norm_emitt
  else
    datum_value = 0.0
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma
  
case ('beam_emittance.z')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%norm_emitt
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma

case ('beam_emittance.a')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%norm_emitt
  else
    call orbit_amplitude_calc (lat%ele_(ix1), tao_lat%orb(ix1), &
                               amp_na = datum_value, particle = electron$)
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma
  
case ('beam_emittance.b')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%norm_emitt
  else
    call orbit_amplitude_calc (lat%ele_(ix1), tao_lat%orb(ix1), &
                               amp_nb = datum_value, particle = electron$)
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma
  
case ('beam_emittance.c')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%norm_emitt
  elseif (track_type == "macro") then
    datum_value = 0
  else
    datum_value = 0
  endif
  call convert_total_energy_to (lat%ele_(ix1)%value(beam_energy$), &
                                              lat%param%particle, gamma)
  datum_value = datum_value / gamma
  
case ('norm_beam_emittance.x')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%norm_emitt
  else
    datum_value = 0.0
  endif
  
case ('norm_beam_emittance.y')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%norm_emitt
  else
    datum_value = 0.0
  endif
  
case ('norm_beam_emittance.z')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%norm_emitt
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif

case ('norm_beam_emittance.a')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%norm_emitt
  else
    call orbit_amplitude_calc (lat%ele_(ix1), tao_lat%orb(ix1), &
                               amp_na = datum_value, particle = electron$)
  endif
  
case ('norm_beam_emittance.b')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%norm_emitt
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%norm_emitt
  else
    call orbit_amplitude_calc (lat%ele_(ix1), tao_lat%orb(ix1), &
                               amp_nb = datum_value, particle = electron$)
  endif
  
case ('norm_beam_emittance.c')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%c%norm_emitt
  elseif (track_type == "macro") then
    datum_value = 0
  else
    datum_value = 0
  endif
  
!---------------------------------------------------------
! Lattice Emittance
  
case ('lat_emittance.a')
  datum_value = tao_lat%modes%a%emittance  

case ('lat_emittance.b')
  datum_value = tao_lat%modes%b%emittance

case ('lat_emittance.c')
  datum_value = tao_lat%modes%z%emittance

!---------------------------------------------------------
case ('chrom.a')
  datum_value = tao_lat%a%chrom

case ('chrom.b')
  datum_value = tao_lat%b%chrom

case ('unstable_ring')
  datum_value = lat%param%growth_rate

case ('dpx_dx') 
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s12$) / tao_lat%bunch_params(ix1)%sigma(s11$)
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpy_dy') 
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s34$) / tao_lat%bunch_params(ix1)%sigma(s33$)
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpz_dz') 
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s56$) / tao_lat%bunch_params(ix1)%sigma(s55$)
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%z%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpa_da') 
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma_normal(s12$) / tao_lat%bunch_params(ix1)%sigma_normal(s11$)
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%a%dpx_dx
  else
    datum_value = 0.0
  endif

case ('dpb_db') 
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma_normal(s34$) / tao_lat%bunch_params(ix1)%sigma_normal(s33$)
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%b%dpx_dx
  else
    datum_value = 0.0
  endif

case ('sigma.x')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s11$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.p_x')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s22$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%x%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.y')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s33$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.p_y')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s44$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%y%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.z')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s55$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%z%sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.p_z')  
  if (track_type == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s66$))
  elseif (track_type == "macro") then
    datum_value = u%macro_beam%params%z%p_sigma
  else
    datum_value = 0.0
  endif
  
case ('sigma.xy')  
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s13$)
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif
  
case ('wire.')  
  if (track_type == "beam") then
    read (data_type(6:), '(a)') angle
    datum_value = tao_do_wire_scan (lat%ele_(ix1), angle, u%beam%beam)
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 0.0
  endif
  
case ('spin.theta')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%theta
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%theta
  endif
  
case ('spin.phi')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%phi
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%phi
  endif
  
case ('spin.polarity')
  if (track_type == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%polarization
  elseif (track_type == "macro") then
    datum_value = 0.0
  else
    datum_value = 1.0
  endif
  
case ('momentum_compaction')
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix0, ix1)
  ele0 => lat%ele_(ix0)
  orb0 => tao_lat%orb(ix0)
  call make_v_mats (ele0, v_mat, v_inv_mat)
  eta_vec = (/ ele0%x%eta, ele0%x%etap, ele0%y%eta, ele0%y%etap /)
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz
  datum_value = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)

case default
  call out_io (s_error$, r_name, 'UNKNOWN DATUM TYPE: ' // datum%data_type)
  call err_exit

end select

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine load_it (vec, ix0, ix1, datum_value, datum, lat, f, coupling_here)

implicit none

type (tao_data_struct) datum
type (ring_struct) lat

real(rp) vec(0:)
real(rp), optional :: f(0:)
real(rp) datum_value

character(20) :: r_name = 'tao_evaluate_a_datum'

integer ix_m, i, ix0, ix1, n_lat
logical, optional :: coupling_here

!

if (ix1 < ix0 .and. lat%param%lattice_type == linear_lattice$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // lat%ele_(ix0)%name, &
                'IS AFTER ENDING ELEMENT: ' // lat%ele_(ix1)%name)
  datum_value = 0
  return
endif
 
!

n_lat = lat%n_ele_use

if (datum%ele0_name == ' ') then
  ix_m = ix1
  if (present(coupling_here)) call coupling_calc (ix_m, datum, lat)

!

else

  if (ix1 < ix0) then   ! wrap around

    if (present(coupling_here)) then
      do i = ix0, n_lat
        call coupling_calc (i, datum, lat)
      enddo
      do i = 0, ix1
        call coupling_calc (i, datum, lat)
      enddo
    endif
  
    select case (datum%merit_type)
    case ('min')
      if (minval(vec(0:ix1)) < minval(vec(ix0:n_lat))) then
        ix_m = minloc (vec(0:ix1), 1) - 1
      else
        ix_m = minloc (vec(ix0:n_lat), 1) + ix0 - 1
      endif
    case ('max')
      if (maxval(vec(0:ix1)) > maxval(vec(ix0:n_lat))) then
        ix_m = maxloc (vec(0:ix1), 1) - 1
      else
        ix_m = maxloc (vec(ix0:n_lat), 1) + ix0 - 1
      endif
    case ('abs_min')
      if (minval(abs(vec(0:ix1))) < minval(abs(vec(ix0:n_lat)))) then
        ix_m = minloc (abs(vec(0:ix1)), 1) - 1
      else
        ix_m = minloc (abs(vec(ix0:n_lat)), 1) + ix0 - 1
      endif
    case ('abs_max')
      if (maxval(abs(vec(0:ix1))) > maxval(abs(vec(ix0:n_lat)))) then
        ix_m = maxloc (abs(vec(0:ix1)), 1) - 1
      else
        ix_m = maxloc (abs(vec(ix0:n_lat)), 1) + ix0 - 1
      endif
    case default
      call out_io (s_abort$, r_name, 'BAD MERIT_TYPE: ' // datum%merit_type, &
                                   'FOR DATUM: ' // tao_datum_name(datum))
      call err_exit
    end select

  else
    if (present(coupling_here)) then
      do i = ix0, ix1
        call coupling_calc (i, datum, lat)
      enddo
    endif

    select case (datum%merit_type)
    case ('min')
      ix_m = minloc (vec(ix0:ix1), 1) + ix0 - 1
    case ('max')
      ix_m = maxloc (vec(ix0:ix1), 1) + ix0 - 1
    case ('abs_min')
      ix_m = minloc (abs(vec(ix0:ix1)), 1) + ix0 - 1
    case ('abs_max')
      ix_m = maxloc (abs(vec(ix0:ix1)), 1) + ix0 - 1
    case default
      call out_io (s_abort$, r_name, &
                    'SINCE THIS DATUM: ' // datum%data_type, &
                    'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                    'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
      call err_exit
    end select

  endif

endif

!

datum%ix_ele_merit = ix_m
datum_value = vec(ix_m)
if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ix_m))
if (present(f)) datum%conversion_factor = f(ix_m)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine coupling_calc (ixd, datum, lat)

implicit none

type (ele_struct), pointer :: ele
type (this_coupling_struct), pointer :: cc_p
type (tao_data_struct) datum
type (ring_struct) lat

integer ie, ixd
real(rp) f, f1, f2

!

if (cc(ixd)%calc_done) return

cc_p => cc(ixd)
ie = datum%ix_ele
ele => lat%ele_(ie)

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
! The ele%r array is used for noise and steering calibration: 
!  r(1,1) = bpm noise (rms)
!  r(1,2) = horizontal alignment calibration
!  r(1,3) = vertical alignment calibration
!  r(1,4) = bpm tilt callibration
!  r(1,5) = BPM scal error (x_reading = a * x_real)
!
! Input: 
!  orb      -- Coord_struct: Orbit position at BPM
!  ele      -- Ele_struct: the BPM
!  noise    -- real(rp): BPM gaussian noise (in meters)
!  axis     -- Integer: x_plane$ or y_plane$
!
! Output:
!  reading  -- Real(rp): BPM reading
!
!-

subroutine tao_read_bpm (orb, ele, axis, reading)

type (coord_struct) orb
type (ele_struct) ele

real(rp) reading
real(rp) ran_num(2)

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
  
  if (axis == x_plane$) then
    reading = (1 + ele%r(1,5)*ran_num(1)) * (ele%r(1,1)*ran_num(2) + &
               (orb%vec(1) - ele%value(x_offset_tot$) + ele%r(1,2)) * &
                                           cos(ele%value(tilt_tot$) + ele%r(1,4)) + &
               (orb%vec(3) - ele%value(y_offset_tot$) + ele%r(1,3)) * &
                                           sin(ele%value(tilt_tot$) + ele%r(1,4)))
  elseif (axis == y_plane$) then
    reading = (1 + ele%r(1,5)*ran_num(1)) * (ele%r(1,1)*ran_num(2) &
              -(orb%vec(1) - ele%value(x_offset_tot$) + ele%r(1,2)) * &
                                           sin(ele%value(tilt_tot$) + ele%r(1,4)) + &
               (orb%vec(3) - ele%value(y_offset_tot$) + ele%r(1,3)) * &
                                           cos(ele%value(tilt_tot$) + ele%r(1,4))) 
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

  call ele_at_s (lat, s_1, ix_ele)

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

  call ele_at_s (lat, s_1, ix_ele)
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
! Subroutine tao_do_wire_scane (ele, theta, beam) result (moment)
!
! Returns the beam's second moment using the wire along the specified angle.
! Keep in mind that the actual correlation axis is 90 degrees off of the 
! wire angle
!
! This simulates a fast wire scanner that performs the scane over only one
! bunch. Obviously, this isn't realistic. Any dynamic effects will not be
! accounted for!
!
! The ele%r array is used for wire scanner resolution:
!  ele%r(1,1) = relative wire resolution RMS (i.e. 0.1 => 10%)
!  ele%r(1,4) = wire angle error in radians rms
!
! Input:
!  theta   -- Real(rp): wire angle wrt x axis (in degrees)
!  beam    -- Beam_struct: contains the beam distribution
!
! Output:
!  moment  -- Real(rp): second moment along axis specified by angle.
!-

function tao_do_wire_scan (ele, theta, beam) result (moment)

implicit none

type (ele_struct) ele
type (beam_struct) beam
real(rp), allocatable, save :: dist(:)

real(rp) theta, theta_rad, moment, ran_num(2)
real(rp) avg

  call ran_gauss (ran_num)
  
  ! angle in radians and correlation angle is 90 off from wire angle
  theta_rad = ele%r(1,4)*ran_num(1) + (theta - 90) * (2.0*pi / 360.0)

  if (.not. allocated (dist)) then
    allocate (dist(size(beam%bunch(1)%particle)))
  elseif (size(dist) .ne. size(beam%bunch(1)%particle)) then
    deallocate (dist)
    allocate (dist(size(beam%bunch(1)%particle)))
  endif

  ! Rotating the wire scanner is equivalent to rotating the beam by -theta

  dist =  beam%bunch(1)%particle%r%vec(1) * cos(-theta_rad ) &
        + beam%bunch(1)%particle%r%vec(3) * sin(-theta_rad)
  
  avg = sum (dist, mask = (beam%bunch(1)%particle%ix_lost == not_lost$)) &
          / count (beam%bunch(1)%particle%ix_lost == not_lost$)
        
  moment = (1 + ele%r(1,1)*ran_num(2)) * sum ((dist-avg)*(dist-avg), &
                 mask = (beam%bunch(1)%particle%ix_lost == not_lost$)) &
          / count (beam%bunch(1)%particle%ix_lost == not_lost$)

end function tao_do_wire_scan

end module tao_data_mod
