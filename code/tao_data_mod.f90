module tao_data_mod

use tao_mod
use tao_evaluate_mod
use spin_mod
use utilities_mod
use measurement_mod

type this_array_struct
  real(rp) cbar(2,2)
  real(rp) k_11a, k_12a, k_12b, k_22b
  real(rp) amp_a, amp_b, amp_na, amp_nb
  logical :: coupling_calc_done = .false.
  logical :: amp_calc_done = .false.
end type

type (this_array_struct), save, allocatable, target, private :: cc(:)

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
!
! Buffer routine for to_phase_and_coupling_reading.
!
! Input:
!   ele -- Ele_struct: The monitor.
!
! Output:
!   bpm_data     -- Bpm_phase_coupling_struct: Monitor values
!   valid_value  -- Logical: Valid data value?
!-

subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

implicit none

type (ele_struct) ele
type (bpm_phase_coupling_struct) bpm_data
type (bpm_phase_coupling_struct), save :: old_bpm_data

integer, save :: ix_ele_old = -1

logical valid_value
logical, save :: err

!

if (ix_ele_old /= ele%ix_ele) then
  call to_phase_and_coupling_reading (ele, old_bpm_data, err)
  ix_ele_old = ele%ix_ele
endif

bpm_data = old_bpm_data
valid_value = .not. err

end subroutine

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
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(iu)%is_on) cycle
  n_data  = n_data + count(s%u(iu)%data(:)%useit_opt)
enddo
if (present(data_value))      call re_allocate (data_value, n_data)
if (present(data_meas_value)) call re_allocate (data_meas_value, n_data)
if (present(data_weight))     call re_allocate (data_weight, n_data)
if (present(data_ix_dModel))  call re_allocate (data_ix_dModel, n_data)

j = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
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
! Subroutine tao_load_data_array (uni, ix_ele, who_dat)
!
! Routine to take data from the model lattice and model orbit
! and put that into the s%u(:)%data(:) arrays.
!
! Input:
!   uni         -- Integer: universe where data resides
!   ix_ele      -- Integer: element to evaluate data at
!   who_dat     -- Integer: Either: model$, base$, or design$
!-

subroutine tao_load_data_array (u, ix_ele, who_dat)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum

integer, pointer :: ix_datum(:)
integer uni, ix_ele, who_dat
integer i

character(20) data_source
character(20) :: r_name = 'tao_load_data_array'

logical good

!

if (ix_ele == 0) call tao_data_coupling_init (u) 
  
! find which datums to evaluate here
if (.not. allocated(u%ix_data(ix_ele)%ix_datum)) return

ix_datum => u%ix_data(ix_ele)%ix_datum
do i = 1, size(ix_datum)
  datum => u%data(ix_datum(i))

  select case (who_dat)
  case (model$) 
    call tao_evaluate_a_datum (datum, u, u%model, datum%model_value, datum%good_model)
    if (datum%ix_ele_merit > -1) datum%s = &
                                    u%model%lat%ele(datum%ix_ele_merit)%s
  case (design$)
    call tao_evaluate_a_datum (datum, u, u%design, datum%design_value, good)
  case (base$)
    call tao_evaluate_a_datum (datum, u, u%base, datum%base_value, good)
  end select
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

cc%coupling_calc_done = .false.
cc%amp_calc_done = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)
!
! Subroutine to put the proper data in the specified datum
!
! If datum results in NaN then datum_value = tiny(1.0_rp)
!
! Input:
!   datum         -- Tao_data_struct: What type of datum
!   u             -- Tao_universe_struct: Which universe to use.
!   tao_lat       -- Tao_lattice_struct: Lattice to use.
!     
! Output:
!   datum          -- Tao_data_struct: 
!     %ix_ele_merit   -- For max/min type constraints: Place where value is max/min. 
!   datum_value   -- Real(rp): Value of the datum.
!   valid_value   -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_data_struct), pointer :: dp
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (normal_modes_struct) mode
type (spin_polar_struct) polar
type (ele_struct), pointer :: ele, ele0
type (coord_struct), pointer :: orb0
type (bpm_phase_coupling_struct) bpm_data
type (taylor_struct), save :: taylor(6) ! Saved taylor map
type (floor_position_struct) floor

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), a_vec(4)
real(rp) gamma, one_pz, vec(2)
real(rp), allocatable, save ::value1(:)

integer, save :: ix_save = -1
integer i, j, k, m, n, ix, ix1, ix0, expnt(6), n_track, n_max

character(*), optional :: why_invalid
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) data_type, data_source, name

logical found, valid_value, err
logical, allocatable, save :: good1(:)

! See if there is a hook for this datum

call tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
if (found) return

! Check range

if (present(why_invalid)) why_invalid =  'ERROR IN TAO_EVALUATE_A_DATUM!' ! Generic

data_source = datum%data_source
data_type = datum%data_type
lat => tao_lat%lat

ix0 = datum%ix_ele0
ix1 = datum%ix_ele

if (datum%ele_name /= '') then
  ix1 = tao_valid_datum_index (lat, ix1, datum, valid_value)
  if (.not. valid_value) return
endif

if (datum%ele0_name /= '') then
  ix0 = tao_valid_datum_index (lat, ix0, datum, valid_value)
  if (.not. valid_value) return
endif

!

valid_value = .false.

n_track = tao_lat%lat%n_ele_track
n_max   = tao_lat%lat%n_ele_max

datum_value = 0           ! default
datum%ix_ele_merit = -1   ! default

if (data_type(1:11) == 'expression:')    data_type = 'expression:'
if (data_type(1:2)  == 'r.')             data_type = 'r.'
if (data_type(1:2)  == 't.')             data_type = 't.'
if (data_type(1:3)  == 'tt.')            data_type = 'tt.'
if (data_type(1:5)  == 'wire.')          data_type = 'wire.'
if (data_type(1:12) == 'periodic.tt.')   data_type = 'periodic.tt.'
if (data_type(1:14) == 'element_param.') data_type = 'element_param.'
if (data_type(1:4)  == 'emit') call convert_total_energy_to ( &
                    lat%ele(ix1)%value(E_tot$), lat%param%particle, gamma)

if (data_source /= "lattice" .and. data_source /= "beam") then
  call out_io (s_error$, r_name, &
          'UNKNOWN DATA_SOURCE: ' // data_source, &
          'FOR DATUM: ' // tao_datum_name(datum))
  call err_exit
endif


if (data_type(1:4) == 'bpm.') then
  if (ix0 /= ix1) then
    call out_io (s_error$, r_name, 'DATUM OVER A REGION NOT YET IMPLEMENTED FOR: ' // &
                                                                 tao_datum_name(datum))
    return
  endif
endif

if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) then
  if (data_source    == 'beam')  return
  if (data_type(1:3) == 'bpm')   return
  if (data_type(1:5) == 'orbit') return
endif

!---------------------------------------------------

select case (data_type)

case ('alpha.a')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%a%alpha, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%alpha
    valid_value = (tao_lat%bunch_params(ix1)%a%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
  
case ('alpha.b')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%b%alpha, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%alpha
    valid_value = (tao_lat%bunch_params(ix1)%b%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('alpha.z')
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%z%alpha
  valid_value = .true.

case ('beta.x')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%x%beta
    valid_value = (tao_lat%bunch_params(ix1)%x%norm_emitt /= 0)
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH data_source: ' // data_source)
    call err_exit
  endif
    
case ('beta.y')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%y%beta
    valid_value = (tao_lat%bunch_params(ix1)%y%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH data_source: ' // data_source)
    call err_exit
  endif

case ('beta.z')
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%z%beta
  valid_value = (tao_lat%bunch_params(ix1)%z%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'

case ('beta.a')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%a%beta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%beta
    valid_value = (tao_lat%bunch_params(ix1)%a%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
    
case ('beta.b')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%b%beta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%beta
    valid_value = (tao_lat%bunch_params(ix1)%b%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('bpm_orbit.x')
  if (data_source == 'beam') return ! bad
  call to_orbit_reading (tao_lat%orb(ix1), lat%ele(ix1), x_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_orbit.y')
  if (data_source == 'beam') return ! bad
  call to_orbit_reading (tao_lat%orb(ix1), lat%ele(ix1), y_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_eta.x')
  if (data_source == 'beam') return ! bad
  vec = (/ lat%ele(ix1)%x%eta, lat%ele(ix1)%y%eta /)
  call to_eta_reading (vec, lat%ele(ix1), x_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_eta.y')
  if (data_source == 'beam') return ! bad
  vec = (/ lat%ele(ix1)%x%eta, lat%ele(ix1)%y%eta /)
  call to_eta_reading (vec, lat%ele(ix1), y_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_phase.a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%phi_a

case ('bpm_phase.b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%phi_b

case ('bpm_k.22a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_22a

case ('bpm_k.12a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_12a

case ('bpm_k.11b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_11b

case ('bpm_k.12b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_12b

case ('bpm_cbar.22a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_22a

case ('bpm_cbar.12a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_12a

case ('bpm_cbar.11b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_11b

case ('bpm_cbar.12b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (lat%ele(ix1), bpm_data, valid_value)
  datum_value = bpm_data%k_12b

case ('cbar.11')
  if (data_source == 'beam') return
  call load_it (cc%cbar(1,1), ix0, ix1, datum_value, valid_value, &
                                     datum, tao_lat, calc_needed = .true.)
case ('cbar.12')
  if (data_source == 'beam') return
  call load_it (cc%cbar(1,2), ix0, ix1, datum_value, valid_value, &
                                     datum, tao_lat, calc_needed = .true.)
case ('cbar.21')
  if (data_source == 'beam') return
  call load_it (cc%cbar(2,1), ix0, ix1, datum_value, valid_value, &
                                     datum, tao_lat, calc_needed = .true.)
case ('cbar.22')
  if (data_source == 'beam') return
  call load_it (cc%cbar(2,2), ix0, ix1, datum_value, valid_value, &
                                     datum, tao_lat, calc_needed = .true.)

case ('chrom.a')
  if (data_source == 'beam') return
  datum_value = tao_lat%a%chrom
  valid_value = .true.

case ('chrom.b')
  if (data_source == 'beam') return
  datum_value = tao_lat%b%chrom
  valid_value = .true.

case ('dpx_dx') 
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%sigma(s12$) / tao_lat%bunch_params(ix1)%sigma(s11$)
  valid_value = .true.

case ('dpy_dy') 
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%sigma(s34$) / tao_lat%bunch_params(ix1)%sigma(s33$)
  valid_value = .true.

case ('dpz_dz') 
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%sigma(s56$) / tao_lat%bunch_params(ix1)%sigma(s55$)
  valid_value = .true.

case ('element_param.')
  if (data_source == 'beam') return ! bad
  call str_upcase (name, datum%data_type(15:))
  ix = attribute_index (lat%ele(ix1), name)
  if (ix < 1) then
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // datum%data_type)
    call err_exit
  endif
  call load_it (lat%ele(0:n_max)%value(ix), &
                                     ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('emit.x', 'norm_emit.x')
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%x%norm_emitt
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.

case ('emit.y', 'norm_emit.y')  
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%y%norm_emitt
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.

case ('emit.z', 'norm_emit.z')
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%z%norm_emitt
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.

case ('emit.a', 'norm_emit.a')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%norm_emitt
  elseif (data_source == 'lattice') then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_a)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_a, &
                              ix0, ix1, datum_value, valid_value, datum, tao_lat)
    else
      datum_value = gamma * tao_lat%modes%a%emittance  
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.
  
case ('emit.b', 'norm_emit.b')  
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%norm_emitt
  elseif (data_source == 'lattice') then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_b)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        valid_value = .false.
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_b, &
                              ix0, ix1, datum_value, valid_value, datum, tao_lat)
    else
      datum_value = gamma * tao_lat%modes%b%emittance
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.

case ('eta.x')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%x%eta
    valid_value = .true.
  else
    call load_it (lat%ele(:)%x%eta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('eta.y')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%y%eta
    valid_value = .true.
  else
    call load_it (lat%ele(:)%y%eta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('etap.x')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%x%etap
    valid_value = .true.
  else
    call load_it (lat%ele(:)%x%etap, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('etap.y')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%y%etap
    valid_value = .true.
  else
    call load_it (lat%ele(:)%y%etap, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('eta.a')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%eta
    valid_value = .true.
  else
    call load_it (lat%ele(:)%a%eta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('eta.b')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%eta
    valid_value = .true.
  else
    call load_it (lat%ele(:)%b%eta, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('etap.a')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%etap
    valid_value = .true.
  else
    call load_it (lat%ele(:)%a%etap, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('etap.b')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%etap
    valid_value = .true.
  else
    call load_it (lat%ele(:)%b%etap, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  endif

case ('e_tot')
  if (data_source == 'beam') return
  call load_it (lat%ele(0:n_track)%value(E_TOT$) * (1+tao_lat%orb(0:n_track)%vec(6)), &
                                     ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('%e_tot')
  if (data_source == 'beam') return
  call load_it (tao_lat%orb(0:n_track)%vec(6), ix0, ix1, &
                                            datum_value, valid_value, datum, tao_lat)
  
case ('expression:')
  if (data_source == 'beam') return ! bad

  ! The point here is that tao_evaluate_stack is much quicker than tao_to_real.
  ! So on the fist time through, construct datum%stack and for subsequent times, use
  ! datum%stack with tao_evaluate_stack.
  if (allocated (datum%stack)) then
    call tao_evaluate_stack (datum%stack, 1, .false., value1, good1, err)
    if (err) call err_exit
    datum_value = value1(1)
    valid_value = good1(1)

  else ! Only do this first time through...
    call tao_to_real (datum%data_type(12:), datum_value, err, .false., valid_value, datum%stack)
    ! Make sure that any datums used in the expression have already been evaluated.
    do i = 1, size(datum%stack)
      if (datum%stack(i)%name == '') cycle
      call tao_find_data (err, datum%stack(i)%name, d_array = d_array, print_err = .false.)
      if (err .or. size(d_array) == 0) cycle  ! Err -> This is not associated then not a datum.
      dp => d_array(1)%d
      if (dp%d1%d2%ix_uni < datum%d1%d2%ix_uni) cycle ! OK
      if (dp%d1%d2%ix_uni == datum%d1%d2%ix_uni .and. dp%ix_data < datum%d1%d2%ix_data) cycle
      call out_io (s_fatal$, r_name, 'DATUM: ' // tao_datum_name(datum), &
                                     'WHICH IS OF TYPE EXPRESSION:' // datum%data_type(12:), &
                                     'THE EXPRESSION HAS A COMPONENT: ' // datum%stack(i)%name, &
                                     'AND THIS COMPONENT IS EVALUATED AFTER THE EXPRESSION!')
      call err_exit
    enddo
  endif

case ('floor.x')
  datum_value = lat%ele(ix1)%floor%x - lat%ele(ix0)%floor%x
  valid_value = .true.

case ('floor.y')
  datum_value = lat%ele(ix1)%floor%y - lat%ele(ix0)%floor%y
  valid_value = .true.

case ('floor.z')
  datum_value = lat%ele(ix1)%floor%z - lat%ele(ix0)%floor%z
  valid_value = .true.

case ('floor.theta')
  datum_value = lat%ele(ix1)%floor%theta - lat%ele(ix0)%floor%theta
  valid_value = .true.

case ('floor.phi')
  datum_value = lat%ele(ix1)%floor%phi - lat%ele(ix0)%floor%phi
  valid_value = .true.

case ('gamma.a')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%a%gamma, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%a%gamma
    valid_value = (tao_lat%bunch_params(ix1)%a%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
  
case ('gamma.b')
  if (data_source == 'lattice') then
    call load_it (lat%ele(:)%b%gamma, ix0, ix1, datum_value, valid_value, datum, tao_lat)
  elseif (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%b%gamma
    valid_value = (tao_lat%bunch_params(ix1)%b%norm_emitt /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('gamma.z')
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%z%gamma
  valid_value = .true.

case ('i5a_e6')
  if (data_source == 'beam') return
  if (ix0 > 0 .or. ix1 > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5a_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      return
    endif
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5a_e6(ix0:ix1))
  else
    datum_value = tao_lat%modes%lin%i5a_e6
  endif
  valid_value = .true.

case ('i5b_e6')
  if (data_source == 'beam') return
  if (ix0 > 0 .or. ix1 > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5b_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      return
    endif
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5b_e6(ix0:ix1))
  else
    datum_value = tao_lat%modes%lin%i5b_e6
  endif
  valid_value = .true.

case ('k.11b')
  if (data_source == 'beam') return
  call load_it (cc%k_11a, ix0, ix1, datum_value, valid_value, &
                                  datum, tao_lat, calc_needed = .true.)
case ('k.12a')
  if (data_source == 'beam') return
  call load_it (cc%k_12a, ix0, ix1, datum_value, valid_value, &
                                  datum, tao_lat, calc_needed = .true.)
case ('k.12b')
  if (data_source == 'beam') return
  call load_it (cc%k_12b, ix0, ix1, datum_value, valid_value, &
                                 datum, tao_lat, calc_needed = .true.)
case ('k.22a')
  if (data_source == 'beam') return
  call load_it (cc%k_22b, ix0, ix1, datum_value, valid_value, &
                                 datum, tao_lat, calc_needed = .true.)

case ('momentum_compaction')
  if (data_source == 'beam') return
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix0, ix1)
  ele0 => lat%ele(ix0)
  orb0 => tao_lat%orb(ix0)
  call make_v_mats (ele0, v_mat, v_inv_mat)
  eta_vec = (/ ele0%a%eta, ele0%a%etap, ele0%b%eta, ele0%b%etap /)
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz
  datum_value = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)
  valid_value = .true.

case ('orbit.x')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(1), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.y')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(3), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.z')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(5), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.p_x')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(2), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.p_y')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(4), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.p_z')
  if (data_source == 'beam') return ! bad
  call load_it (tao_lat%orb(:)%vec(6), ix0, ix1, datum_value, valid_value, datum, tao_lat)

case ('orbit.amp_a')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_a, ix0, ix1, datum_value, valid_value, datum, tao_lat, calc_needed = .true.)

case ('orbit.amp_b')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_b, ix0, ix1, datum_value, valid_value, datum, tao_lat, calc_needed = .true.)

case ('orbit.norm_amp_a')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_na, ix0, ix1, datum_value, valid_value, datum, tao_lat, calc_needed = .true.)

case ('orbit.norm_amp_b')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_nb, ix0, ix1, datum_value, valid_value, datum, tao_lat, calc_needed = .true.)

case ('periodic.tt.')
  if (data_source == 'beam') return
  if (lat%param%lattice_type /= circular_lattice$) then
    call out_io (s_fatal$, r_name, 'LATTICE MUST BE CIRCULAR FOR A DATUM LIKE: ' // &
                                                                        datum%data_type)
    call err_exit
  endif
  
  call transfer_map_calc (lat, taylor, ix1, ix1, one_turn = .true.)
  do i = 1, 4
    call add_taylor_term (taylor(i), -1.0_rp, i)
  enddo
  call taylor_inverse (taylor, taylor)

  expnt = 0
  i = tao_read_this_index (datum%data_type, 13); if (i == 0) return
  do j = 14, 24
    if (datum%data_type(j:j) == ' ') exit
    k = tao_read_this_index (datum%data_type, j); if (k == 0) return
    expnt(k) = expnt(k) + 1
  enddo

  datum_value = taylor_coef (taylor(i), expnt)
  valid_value = .true.

case ('phase.a')
  if (data_source == 'beam') return ! bad
  datum_value = lat%ele(ix1)%a%phi - lat%ele(ix0)%a%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%a%phi + lat%ele(n_track)%a%phi 
  valid_value = .true.

case ('phase.b')
  if (data_source == 'beam') return ! bad
  datum_value = lat%ele(ix1)%b%phi - lat%ele(ix0)%b%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%b%phi + lat%ele(n_track)%b%phi 
  valid_value = .true.

case ('phase_frac.a')
  if (data_source == 'beam') return ! bad
  datum_value = lat%ele(ix1)%a%phi - lat%ele(ix0)%a%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%a%phi + lat%ele(n_track)%a%phi 
  datum_value = modulo2(datum_value, pi)
  valid_value = .true.

case ('phase_frac.b')
  if (data_source == 'beam') return ! bad
  datum_value = lat%ele(ix1)%b%phi - lat%ele(ix0)%b%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%b%phi + lat%ele(n_track)%b%phi 
  datum_value = modulo2(datum_value, pi)
  valid_value = .true.

case ('phase_frac_diff')
  if (data_source == 'beam') return ! bad
  px = lat%ele(ix1)%a%phi - lat%ele(ix0)%a%phi
  if (ix0 > ix1) px = px - lat%ele(0)%a%phi + lat%ele(n_track)%a%phi 
  py = lat%ele(ix1)%b%phi - lat%ele(ix0)%b%phi
  if (ix0 > ix1) py = py - lat%ele(0)%b%phi + lat%ele(n_track)%b%phi 
  datum_value = modulo2 (px - py, pi)
  valid_value = .true.

case ('r.')
  if (data_source == 'beam') return
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix0, ix1)
  datum_value = mat6(i, j)
  valid_value = .true.

case ('rel_floor.x', 'rel_floor.y', 'rel_floor.z', 'rel_floor.theta', 'rel_floor.phi')
  call init_floor (floor)

  if (ix1 > ix0) then
    do i = ix0+1, ix1
      call ele_geometry (floor, lat%ele(i), floor)
    enddo
  else
    do i = ix0+1, lat%n_ele_track
      call ele_geometry (floor, lat%ele(i), floor)
    enddo
    do i = 1, ix1
      call ele_geometry (floor, lat%ele(i), floor)
    enddo
  endif

  select case (data_type)
  case ('rel_floor.x')
    datum_value = lat%ele(ix1)%floor%x - lat%ele(ix0)%floor%x
  case ('rel_floor.y')
    datum_value = lat%ele(ix1)%floor%y - lat%ele(ix0)%floor%y
  case ('rel_floor.z')
    datum_value = lat%ele(ix1)%floor%z - lat%ele(ix0)%floor%z
  case ('rel_floor.theta')
    datum_value = lat%ele(ix1)%floor%theta - lat%ele(ix0)%floor%theta
  case ('rel_floor.phi')
    datum_value = lat%ele(ix1)%floor%phi - lat%ele(ix0)%floor%phi
  end select
  valid_value = .true.

case ('sigma.x')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s11$))
  valid_value = .true.
  
case ('sigma.p_x')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s22$))
  valid_value = .true.
  
case ('sigma.y')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s33$))
  valid_value = .true.
  
case ('sigma.p_y')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s44$))
  valid_value = .true.
  
case ('sigma.z')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s55$))
  valid_value = .true.
  
case ('sigma.p_z')  
  if (data_source == 'lattice') return
  datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s66$))
  valid_value = .true.
  
case ('sigma.xy')  
  if (data_source == 'lattice') return
  datum_value = tao_lat%bunch_params(ix1)%sigma(s13$)
  valid_value = .true.
  
case ('spin.theta')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%spin%theta
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%theta
  endif
  valid_value = .true.
  
case ('spin.phi')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%spin%phi
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%phi
  endif
  valid_value = .true.
  
case ('spin.polarity')
  if (data_source == 'beam') then
    datum_value = tao_lat%bunch_params(ix1)%spin%polarization
  else
    datum_value = 1.0
  endif
  valid_value = .true.
  
case ('s_position') 
  if (data_source == 'beam') return
  if (ix0 >= 0) then
    datum_value = lat%ele(ix1)%s - lat%ele(ix0)%s
  else
    datum_value = lat%ele(ix1)%s 
  endif
  valid_value = .true.

case ('tune.a')
  if (data_source == 'beam') return ! bad
  datum_value = lat%a%tune
  valid_value = .true.

case ('tune.b')
  if (data_source == 'beam') return ! bad
  datum_value = lat%b%tune
  valid_value = .true.

case ('t.', 'tt.')
  if (data_source == 'beam') return
  expnt = 0
  if (data_type == 't.') then
    i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
    do j = 4, 5
      k = tao_read_this_index (datum%data_type, j); if (k == 0) return
      expnt(k) = expnt(k) + 1
    enddo
  else
    i = tao_read_this_index (datum%data_type, 4); if (i == 0) return
    do j = 5, 15
      if (datum%data_type(j:j) == ' ') exit
      k = tao_read_this_index (datum%data_type, j); if (k == 0) return
      expnt(k) = expnt(k) + 1
    enddo
  endif

  if (tao_com%ix0_taylor /= ix0 .or. tao_com%ix1_taylor /= ix1) then
    if (tao_com%ix0_taylor == ix0 .and. ix1 > tao_com%ix1_taylor) then
      call transfer_map_calc (lat, taylor, tao_com%ix1_taylor, ix1, unit_start = .false.)
    else
      call transfer_map_calc (lat, taylor, ix0, ix1)
    endif
    tao_com%ix0_taylor = ix0
    tao_com%ix1_taylor = ix1
  endif
  datum_value = taylor_coef (taylor(i), expnt)
  valid_value = .true.

case ('unstable_orbit')
  valid_value = .true.
  if (lat%param%lattice_type /= linear_lattice$) return
  if (datum%ele_name == '') ix1 = lat%n_ele_track
  if (data_source == 'beam') then
    do i = 1, ix1
      datum_value = datum_value + (1 + ix1 - i) * u%ele(i)%n_lost_here
    enddo
    datum_value = datum_value / size(u%current_beam%bunch(s%global%bunch_to_plot)%particle)
  else
    if (lat%param%ix_lost == not_lost$) return
    datum_value = max(0, 1 + ix1 - lat%param%ix_lost)
  endif

case ('unstable_ring')
  if (data_source == 'beam') return
  datum_value = lat%param%growth_rate
  ! unstable_penalty is needed since at the meta stable borderline the growth rate is zero.
  if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty
  valid_value = .true.

case ('wall')
  if (data_source == 'beam') return
  print *, 'NOT YET IMPLEMENTED...'
  call err_exit

case ('wire.')  
  if (data_source == 'lattice') return
  read (data_type(6:), '(a)') angle
  datum_value = tao_do_wire_scan (lat%ele(ix1), angle, u%current_beam)
  valid_value = .true.
  
case default
  call out_io (s_error$, r_name, 'UNKNOWN DATUM TYPE: ' // datum%data_type)
  call err_exit

end select

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine load_it (vec, ix0, ix1, datum_value, valid_value, datum, tao_lat, calc_needed)

implicit none

type (tao_data_struct) datum
type (tao_lattice_struct) tao_lat

real(rp) vec(0:)
real(rp) datum_value

character(20) :: r_name = 'tao_evaluate_a_datum'

integer ix_m, i, ix0, ix1, n_track, ix_m2
logical valid_value
logical, optional :: calc_needed

!

if (ix1 < ix0 .and. tao_lat%lat%param%lattice_type == linear_lattice$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // tao_lat%lat%ele(ix0)%name, &
                'IS AFTER ENDING ELEMENT: ' // tao_lat%lat%ele(ix1)%name)
  valid_value = .false.
  return
endif
 
!

n_track = tao_lat%lat%n_ele_track

if (datum%ele0_name == ' ') then
  ix_m = ix1
  if (present(calc_needed)) call data_calc (ix_m, datum, tao_lat)
  datum_value = vec(ix_m)
  if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ix_m))

elseif (datum%merit_type == 'match') then

  ix_m = ix1
  if (present(calc_needed)) call data_calc (ix0, datum, tao_lat)
  if (present(calc_needed)) call data_calc (ix1, datum, tao_lat)
  datum_value = vec(ix1) - vec(ix0)

else

  if (ix1 < ix0) then   ! wrap around

    if (present(calc_needed)) then
      do i = ix0, n_track
        call data_calc (i, datum, tao_lat)
      enddo
      do i = 0, ix1
        call data_calc (i, datum, tao_lat)
      enddo
    endif
  
    select case (datum%merit_type)
    case ('min')
      ix_m = minloc (vec(0:ix1), 1) - 1
      ix_m2 = minloc (vec(ix0:n_track), 1) + ix0 - 1
      if (vec(ix_m2) < vec(ix_m2)) ix_m = ix_m2
      datum_value = vec(ix_m)

    case ('max')
      ix_m = maxloc (vec(0:ix1), 1) - 1
      ix_m2 = maxloc (vec(ix0:n_track), 1) + ix0 - 1
      if (vec(ix_m2) > vec(ix_m2)) ix_m = ix_m2
      datum_value = vec(ix_m)

    case ('abs_min')
      ix_m = minloc (abs(vec(0:ix1)), 1) - 1
      ix_m2 = minloc (abs(vec(ix0:n_track)), 1) + ix0 - 1
      if (abs(vec(ix_m2)) < abs(vec(ix_m2))) ix_m = ix_m2
      datum_value = abs(vec(ix_m))

    case ('abs_max')
      ix_m = maxloc (abs(vec(0:ix1)), 1) - 1
      ix_m2 = maxloc (abs(vec(ix0:n_track)), 1) + ix0 - 1
      if (abs(vec(ix_m2)) > abs(vec(ix_m2))) ix_m = ix_m2
      datum_value = abs(vec(ix_m))

    case ('int_min')
      datum_value = 0; ix_m = -1
      call integrate_min (ix1, tao_lat%lat%n_ele_track, datum_value, ix_m, tao_lat%lat, vec, datum)
      call integrate_min (0, ix0, datum_value, ix_m, tao_lat%lat, vec, datum)

    case ('int_max')
      datum_value = 0; ix_m = -1
      call integrate_max (ix1, tao_lat%lat%n_ele_track, datum_value, ix_m, tao_lat%lat, vec, datum)
      call integrate_max (0, ix0, datum_value, ix_m, tao_lat%lat, vec, datum)

    case default
      call out_io (s_abort$, r_name, 'BAD MERIT_TYPE: ' // datum%merit_type, &
                                   'FOR DATUM: ' // tao_datum_name(datum))
      call err_exit
    end select

  else
    if (present(calc_needed)) then
      do i = ix0, ix1
        call data_calc (i, datum, tao_lat)
      enddo
    endif

    select case (datum%merit_type)
    case ('min')
      ix_m = minloc (vec(ix0:ix1), 1) + ix0 - 1
      datum_value = vec(ix_m)

    case ('max')
      ix_m = maxloc (vec(ix0:ix1), 1) + ix0 - 1
      datum_value = vec(ix_m)

    case ('abs_min')
      ix_m = minloc (abs(vec(ix0:ix1)), 1) + ix0 - 1
      datum_value = abs(vec(ix_m))

    case ('abs_max')
      ix_m = maxloc (abs(vec(ix0:ix1)), 1) + ix0 - 1
      datum_value = abs(vec(ix_m))

    case ('int_min')
      datum_value = 0; ix_m = -1
      call integrate_min (ix0, ix1, datum_value, ix_m, tao_lat%lat, vec, datum)

    case ('int_max')
      datum_value = 0; ix_m = -1
      call integrate_max (ix0, ix1, datum_value, ix_m, tao_lat%lat, vec, datum)

    case default
      call out_io (s_abort$, r_name, &
                    'SINCE THIS DATUM: ' // tao_datum_name(datum), &
                    'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                    'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
      call err_exit
    end select

  endif

endif

!

datum%ix_ele_merit = ix_m
valid_value = .true.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_min (ix0, ix1, datum_value, ix_m, lat, vec, datum)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix0, ix1

!

val0 = datum%meas_value

do i = ix0, ix1
  if (ix_m < 0) ix_m = i
  if (vec(i) < vec(ix_m)) ix_m = i
  if (i == ix0) cycle
  dv0 = val0 - vec(i-1) 
  dv1 = val0 - vec(i)
  ds = lat%ele(i)%s - lat%ele(i-1)%s
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_max (ix0, ix1, datum_value, ix_m, lat, vec, datum)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix0, ix1

!

val0 = datum%meas_value

do i = ix0, ix1
  if (ix_m < 0) ix_m = i
  if (vec(i) > vec(ix_m)) ix_m = i
  if (i == ix0) cycle
  dv0 = vec(i-1) - val0
  dv1 = vec(i) - val0
  ds = lat%ele(i)%s - lat%ele(i-1)%s
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine data_calc (ixd, datum, tao_lat)

implicit none

type (ele_struct), pointer :: ele
type (this_array_struct), pointer :: cc_p
type (tao_data_struct) datum
type (tao_lattice_struct) tao_lat

integer ie, ixd
real(rp) f, f1, f2

!

cc_p => cc(ixd)
ie = datum%ix_ele
ele => tao_lat%lat%ele(ie)


select case (datum%data_type)
case ('k.11b', 'k.12a', 'k.12b', 'k.22a', 'cbar.11', 'cbar.12', 'cbar.21', 'cbar.22')
  if (cc_p%coupling_calc_done) return

  call c_to_cbar (ele, cc_p%cbar)
  f = sqrt(ele%a%beta/ele%b%beta) 
  f1 = f / ele%gamma_c
  f2 = 1 / (f * ele%gamma_c)

  cc_p%k_11a = cc_p%cbar(1,1) * f1
  cc_p%k_12a = cc_p%cbar(1,2) * f2
  cc_p%k_12b = cc_p%cbar(1,2) * f1
  cc_p%k_22b = cc_p%cbar(2,2) * f2
  cc_p%coupling_calc_done = .true.

! Amplitude calc

case ('orbit.amp_a', 'orbit.amp_b', 'orbit.norm_amp_a', 'orbit.norm_amp_b')
  if (cc_p%amp_calc_done) return
  call orbit_amplitude_calc (ele, tao_lat%orb(ie), cc_p%amp_a, cc_p%amp_b, &
                           cc_p%amp_na, cc_p%amp_nb, tao_lat%lat%param%particle)
  cc_p%amp_calc_done = .true.
end select

end subroutine data_calc

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_do_wire_scane (ele, wire_params, theta, beam) result (moment)
!
! Returns the beam's second moment using the wire along the specified angle.
! Keep in mind that the actual correlation axis is 90 degrees off of the 
! wire angle
!
! This simulates a fast wire scanner that performs the scane over only one
! bunch. Obviously, this isn't realistic. Any dynamic effects will not be
! accounted for!
!
! Input:
!  ele         -- Element_struct: 
!    %value(noise$) -- relative wire resolution RMS 
!    %value(tilt$)  -- wire angle error in radians rms.
!  theta       -- Real(rp): wire angle wrt x axis (in degrees)
!  beam        -- Beam_struct: contains the beam distribution
!
! Output:
!   moment  -- Real(rp): second moment along axis specified by angle.
!-

function tao_do_wire_scan (ele, theta, beam) result (moment)

use random_mod

implicit none

type (ele_struct) ele
type (beam_struct) beam

real(rp), allocatable, save :: dist(:)
real(rp) theta, theta_rad, moment, ran_num(2)
real(rp) avg

!

call ran_gauss (ran_num)

! angle in radians and correlation angle is 90 off from wire angle
theta_rad = ele%value(tilt_tot$)+ (theta - 90) * (2.0*pi / 360.0)

if (.not. allocated (dist)) then
  allocate (dist(size(beam%bunch(1)%particle)))
elseif (size(dist) /= size(beam%bunch(1)%particle)) then
  deallocate (dist)
  allocate (dist(size(beam%bunch(1)%particle)))
endif

! Rotating the wire scanner is equivalent to rotating the beam by -theta

dist =  beam%bunch(1)%particle%r%vec(1) * cos(-theta_rad ) &
        + beam%bunch(1)%particle%r%vec(3) * sin(-theta_rad)
  
avg = sum (dist, mask = (beam%bunch(1)%particle%ix_lost == not_lost$)) &
          / count (beam%bunch(1)%particle%ix_lost == not_lost$)
        
moment = (1 + ele%value(noise$)*ran_num(2)) * sum ((dist-avg)*(dist-avg), &
                 mask = (beam%bunch(1)%particle%ix_lost == not_lost$)) &
          / count (beam%bunch(1)%particle%ix_lost == not_lost$)

end function tao_do_wire_scan

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Function tao_valid_datum_index (lat, ix_ele, datum, valid) result (ix_loc)
! 
! Routine to see if an element index corresponds to an element with a definite 
! location such as an overlay or multipass element. 
!
! If the element is a super_lord then the super_slave element at the exit end
! will be returned. Otherwise ix_loc will be set to ix_ele.
!
! Input:
!   lat    -- Lat_struct: Lattice
!   ix_ele -- Integer: Index of element.
!   datum  -- Tao_data_struct: Used for error messages
!   valid  -- Logical: Set False if element does not have a definite location.
!               Set True otherwise
!   ix_loc -- Integer: Location of element in the tracking part of the lat%ele(:) array.
!-

function tao_valid_datum_index (lat, ix_ele, datum, valid) result (ix_loc)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum
integer ix_ele, ix_loc, ixc
logical valid
character(40) :: r_name = 'tao_valid_datum_index'

!

ix_loc = ix_ele
valid = .true.

if (ix_ele < 0 .or. ix_ele > lat%n_ele_max) then
  call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE! \i5\ ', ix_ele)
  valid = .false.
  return
endif

if (datum%data_type(1:14) == 'element_param.') return

if (ix_ele <= lat%n_ele_track) return

if (lat%ele(ix_ele)%control_type == super_lord$) then
  ixc = lat%ele(ix_ele)%ix2_slave
  ix_loc = lat%control(ixc)%ix_slave
  return
endif

valid = .false.
call out_io (s_error$, r_name, &
            'ELEMENT: ' // trim(lat%ele(ix_ele)%name) // &
            '    WHICH IS A: ' // control_name(lat%ele(ix_ele)%control_type), &
            'CANNOT BE USED IN DEFINING A DATUM SINCE IT DOES NOT HAVE ', &
            '   A DEFINITE LOCATION IN THE LATTICE.', &
            'FOR DATUM: ' // tao_datum_name(datum) )

end function

end module tao_data_mod
