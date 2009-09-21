module tao_data_and_eval_mod

use tao_mod
use spin_mod
use utilities_mod
use measurement_mod
use lat_geometry_mod

type this_array_struct
  real(rp) cbar(2,2)
  real(rp) k_11a, k_12a, k_12b, k_22b
  real(rp) amp_a, amp_b, amp_na, amp_nb
  real(rp) :: one = 1.0
  logical :: coupling_calc_done = .false.
  logical :: amp_calc_done = .false.
end type

type (this_array_struct), save, allocatable, target, private :: cc(:)

! used for parsing expressions

integer, parameter, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter, private :: numeric$ = 100, var$ = 101, var_on_the_fly$ = 102

integer, parameter, private :: eval_level(22) = (/ 1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 /)

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_data_on_the_fly (err, data_name, values, print_err)
!
! Routine to evaluate data "on-the-fly" of the form 
!     <universe>@dat:<source>:<data_type>[<ix_ele_start>&<ix_ele>]|<component>
!
! Input:
!   data_name -- Character(*): data name.
!   print_err -- Logical, optional :: Print error message? Default is True.
!
! Output:
!   err       -- Logical: True if there is an error. False otherwise
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_data_on_the_fly (err, data_name, values, print_err)

implicit none

type (tao_data_struct) datum
type (tao_universe_struct), pointer :: u

character(*) data_name
character(60) name, ele_name, component
character(40) :: r_name = 'tao_evaluate_data_on_the_fly'

real(rp), allocatable :: values(:)

integer i, j, num, ix, ix1, ios, n_tot
integer, allocatable, save :: ix_ele(:)

logical err, printit, valid
logical, optional :: print_err
logical, allocatable, save :: this_u(:)

!

call tao_pick_universe (data_name, name, this_u, err)
if (err) return

err = .true.
printit = logic_option (.true., print_err)
name = name(5:)  ! Strip off 'dat:'

! Get data source. Default is 'lattice'

ix1 = index(name, ':');  if (ix1 == 0) return
datum%data_source = name(1:ix1-1)
if (datum%data_source == '') datum%data_source = 'lattice'
name = name(ix1+1:)

! Get data type

ix1 = index(name, '[');  if (ix1 == 0) return
datum%data_type = name(1:ix1-1)
name = name(ix1+1:)

! Get ele_start

ix1 = index(name, ']');  if (ix1 == 0) return
ix = index(name, '&')
datum%ix_ele_start = -1
datum%ele_start_name = ''
if (ix /= 0 .and. ix < ix1) then
  if (is_integer(name(:ix-1))) then
    read (name(:ix-1), *, iostat = ios) datum%ix_ele_start
    if (ios /= 0) then
      if (printit) call out_io (s_error$, r_name, 'BAD ELE_START: ' // data_name)
      return
    endif
  else
    datum%ele_start_name = name(:ix-1)
  endif
  name = name(ix+1:)
endif

! Get ele

ix = index(name, ']')
ele_name = name(:ix-1)
name = name(ix+1:)

! Get component

if (name(1:1) == '|') then
  component = name(2:)
elseif (name(1:1) == ' ') then
  component = 'model'
else
  if (printit) call out_io (s_error$, r_name, 'COMPONENT: ' // data_name)
  return
endif

! Evaluate

n_tot = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(i)) cycle
  u => s%u(i)
  call elements_locator (ele_name, u%model%lat, ix_ele, err)
  if (err) return
  call re_allocate (values, n_tot + size(ix_ele))
  do j = 1, size(ix_ele)
    datum%ix_ele = ix_ele(j)
    select case (component)
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, values(n_tot+j), valid)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, values(n_tot+j), valid)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, values(n_tot+j), valid)
    case default
      call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // data_name)
      return
    end select
  enddo
  n_tot = n_tot + size(values)
enddo

err = .false.

end subroutine tao_evaluate_data_on_the_fly

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
!   data_value(:)       -- Real(Rp), allocatable, optional: Data model values.
!   data_weight(:)      -- Real(Rp), allocatable, optional: Data weights in the merit function.
!   data_meas_value(:)  -- Real(Rp), allocatable, optional: Data values when the data was taken.
!   data_ix_dModel(:)   -- Integer, allocatable, optional: Data ix_dModel indices
!-

subroutine tao_get_data (data_value, data_weight, data_meas_value, data_ix_dModel)

implicit none

real(rp), allocatable, optional :: data_value(:), data_meas_value(:), data_weight(:)
integer, allocatable, optional :: data_ix_dModel(:)

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
    if (datum%ix_ele_merit > -1) then
      datum%s = u%model%lat%ele(datum%ix_ele_merit)%s
    endif
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

recursive subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_data_struct), pointer :: dp
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (normal_modes_struct) mode
type (spin_polar_struct) polar
type (ele_struct), pointer :: ele, ele_start, ele_ref
type (coord_struct), pointer :: orb0
type (bpm_phase_coupling_struct) bpm_data
type (taylor_struct), save :: taylor(6) ! Saved taylor map
type (floor_position_struct) floor
type (coord_struct), pointer :: orbit(:)
type (branch_struct), pointer :: branch
type (lat_ele_loc_struct) loc_start, loc_ref, loc_ele
type (bunch_params_struct), pointer :: bunch_params(:)
type (tao_element_struct), pointer :: uni_ele(:)

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py, vec2(2)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), a_vec(4)
real(rp) gamma, one_pz, w0_mat(3,3), w_mat(3,3), vec3(3)
real(rp), allocatable, save ::value1(:), value_vec(:)
real(rp) theta, phi, psi

integer, save :: ix_save = -1
integer i, j, k, m, n, ix, ix_ele, ix_start, ix_ref, expnt(6), n_track, n_max, iz

character(*), optional :: why_invalid
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) data_type, data_source, name

logical found, valid_value, err
logical, allocatable, save :: good1(:)

! See if there is a hook for this datum

call tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
if (found) return

! Check range

data_source = datum%data_source
data_type = datum%data_type
lat => tao_lat%lat

loc_ele = tao_valid_datum_index (lat, datum%ele_name, datum%ix_ele, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ele   = loc_ele%ix_ele

loc_ref = tao_valid_datum_index (lat, datum%ele_ref_name, datum%ix_ele_ref, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ref   = loc_ref%ix_ele

loc_start = tao_valid_datum_index (lat, datum%ele_start_name, datum%ix_ele_start, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_start = loc_start%ix_ele

!

branch => lat%branch(datum%ix_branch)
orbit => tao_lat%lat_branch(datum%ix_branch)%orbit
bunch_params => tao_lat%lat_branch(datum%ix_branch)%bunch_params
uni_ele => u%uni_branch(datum%ix_branch)%ele

! ix_ref and ele_ref will point to the 0th element even if there is no refrence

nullify(ele_start, ele_ref, ele)
if (ix_ele >= 0)   ele       => branch%ele(ix_ele)
if (ix_start >= 0) ele_start => branch%ele(ix_start)
if (ix_ref >= 0)   ele_ref   => branch%ele(ix_ref)

valid_value = .false.

n_track = branch%n_ele_track
n_max   = branch%n_ele_max

datum_value = 0           ! default
datum%ix_ele_merit = -1   ! default

if (data_type(1:11) == 'expression:')    data_type = 'expression:'
if (data_type(1:2)  == 'r.')             data_type = 'r.'
if (data_type(1:2)  == 't.')             data_type = 't.'
if (data_type(1:3)  == 'tt.')            data_type = 'tt.'
if (data_type(1:5)  == 'wire.')          data_type = 'wire.'
if (data_type(1:12) == 'periodic.tt.')   data_type = 'periodic.tt.'
if (data_type(1:14) == 'element_param.') data_type = 'element_param.'

if (data_source /= "lattice" .and. data_source /= "beam") then
  call out_io (s_error$, r_name, &
          'UNKNOWN DATA_SOURCE: ' // data_source, &
          'FOR DATUM: ' // tao_datum_name(datum))
  why_invalid = 'UNKNOWN DATA_SOURCE: ' // data_source
  return
endif


if (data_type(1:4) == 'bpm.') then
  if (ix_start /= ix_ele) then
    call out_io (s_error$, r_name, 'DATUM OVER A REGION NOT YET IMPLEMENTED FOR: ' // &
                                                                 tao_datum_name(datum))
    if (present(why_invalid)) why_invalid = 'DATUM OVER A REGION NOT YET IMPLEMENTED FOR: ' // tao_datum_name(datum)
    return
  endif
endif

if (lat%param%ix_lost /= not_lost$ .and. ix_ele >= lat%param%ix_lost) then
  if (data_source == 'beam' .or. data_type(1:3) == 'bpm' .or. &
                                 data_type(1:5) == 'orbit') then
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE DUE TO PARTICLE LOSS.'
    return
  endif
endif

!---------------------------------------------------

if (present(why_invalid)) why_invalid = &
     'DATA_SOURCE = ' // trim (data_source) // ' NOT VALID FOR: ' // trim(tao_datum_name(datum)) ! Default

select case (data_type)

case ('alpha.a')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%a%alpha, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%a%alpha, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
  
case ('alpha.b')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%b%alpha, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%b%alpha, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('alpha.z')
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%z%alpha, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('beta.x')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%x%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%x%norm_emit /= 0)
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH DATA_SOURCE: ' // data_source)
    why_invalid = 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.x DATA_TYPE'
  endif
    
case ('beta.y')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%y%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%y%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH data_source: ' // data_source)
    why_invalid = 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.y DATA_TYPE'
  endif

case ('beta.z')
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%z%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  valid_value = valid_value .and. (bunch_params(ix_ele)%z%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'

case ('beta.a')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%a%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%a%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
    
case ('beta.b')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%b%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%b%beta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('bpm_orbit.x')
  if (data_source == 'beam') return ! bad
  call to_orbit_reading (orbit(ix_ele), ele, x_plane$, datum_value, err)
  valid_value = .not. (err .or. (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost))

case ('bpm_orbit.y')
  if (data_source == 'beam') return ! bad
  call to_orbit_reading (orbit(ix_ele), ele, y_plane$, datum_value, err)
  valid_value = .not. (err .or. (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost))

case ('bpm_eta.x')
  if (data_source == 'beam') return ! bad
  vec2 = (/ ele%x%eta, ele%y%eta /)
  call to_eta_reading (vec2, ele, x_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_eta.y')
  if (data_source == 'beam') return ! bad
  vec2 = (/ ele%x%eta, ele%y%eta /)
  call to_eta_reading (vec2, ele, y_plane$, datum_value, err)
  valid_value = .not. err

case ('bpm_phase.a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%phi_a

case ('bpm_phase.b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%phi_b

case ('bpm_k.22a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_22a

case ('bpm_k.12a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_12a

case ('bpm_k.11b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_11b

case ('bpm_k.12b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_12b

case ('bpm_cbar.22a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_22a

case ('bpm_cbar.12a')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_12a

case ('bpm_cbar.11b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_11b

case ('bpm_cbar.12b')
  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
  datum_value = bpm_data%k_12b

case ('cbar.11')
  if (data_source == 'beam') return
  call load_it (cc%cbar(1,1), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

case ('cbar.12')
  if (data_source == 'beam') return
  call load_it (cc%cbar(1,2), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

case ('cbar.21')
  if (data_source == 'beam') return
  call load_it (cc%cbar(2,1), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

case ('cbar.22')
  if (data_source == 'beam') return
  call load_it (cc%cbar(2,2), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

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
  datum_value = bunch_params(ix_ele)%sigma(s12$) / bunch_params(ix_ele)%sigma(s11$)
  valid_value = .true.

case ('dpy_dy') 
  if (data_source == 'lattice') return
  datum_value = bunch_params(ix_ele)%sigma(s34$) / bunch_params(ix_ele)%sigma(s33$)
  valid_value = .true.

case ('dpz_dz') 
  if (data_source == 'lattice') return
  datum_value = bunch_params(ix_ele)%sigma(s56$) / bunch_params(ix_ele)%sigma(s55$)
  valid_value = .true.

case ('e_tot')
  if (data_source == 'beam') return
  call load_it (branch%ele(0:n_track)%value(E_TOT$) * (1+orbit(0:n_track)%vec(6)), &
                            loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('%e_tot')
  if (data_source == 'beam') return
  call load_it (orbit(0:n_track)%vec(6), loc_ref, loc_start, loc_ele, &
                                     datum_value, valid_value, datum, lat, why_invalid)
  
case ('element_param.')
  if (data_source == 'beam') return ! bad
  call str_upcase (name, datum%data_type(15:))
  ix = attribute_index (ele, name)
  if (ix < 1) then
    call out_io (s_error$, r_name, 'BAD ELEMENT ATTRIBUTE : ' // datum%data_type)
    why_invalid = 'BAD ELEMENT ATTRIBUTE : ' // datum%data_type
  endif
  call load_it (branch%ele(0:n_max)%value(ix), &
                          loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('emit.x', 'norm_emit.x')
  call convert_total_energy_to (ele%value(E_tot$), lat%param%particle, gamma)
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%x%norm_emit, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.y', 'norm_emit.y')  
  call convert_total_energy_to (ele%value(E_tot$), lat%param%particle, gamma)
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%y%norm_emit, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.z', 'norm_emit.z')
  call convert_total_energy_to (ele%value(E_tot$), lat%param%particle, gamma)
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%z%norm_emit, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.a', 'norm_emit.a')
  call convert_total_energy_to (lat%ele(0)%value(E_tot$), lat%param%particle, gamma)
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%a%norm_emit, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'lattice') then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_a)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_a, &
                              loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    else
      datum_value = gamma * tao_lat%modes%a%emittance  
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.
  
case ('emit.b', 'norm_emit.b')  
  call convert_total_energy_to (lat%ele(0)%value(E_tot$), lat%param%particle, gamma)
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%b%norm_emit, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'lattice') then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_b)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        valid_value = .false.
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_b, &
                              loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    else
      datum_value = gamma * tao_lat%modes%b%emittance
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  valid_value = .true.

case ('eta.x')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%x%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%x%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('eta.y')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%y%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%y%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('etap.x')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%x%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%x%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('etap.y')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%y%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%y%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('eta.a')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%a%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%a%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('eta.b')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%b%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%b%eta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('etap.a')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%a%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%a%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('etap.b')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%b%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = .true.
  else
    call load_it (branch%ele(:)%b%etap, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  endif

case ('expression:')
  if (data_source == 'beam') return ! bad

  ! The point here is that tao_evaluate_stack is much quicker than tao_to_real.
  ! So on the fist time through, construct datum%stack and for subsequent times, use
  ! datum%stack with tao_evaluate_stack.
  if (allocated (datum%stack)) then
    call tao_evaluate_stack (datum%stack, 1, .false., value1, good1, err, .true.)
    if (err) return
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
      if (dp%d1%d2%ix_uni == datum%d1%d2%ix_uni .and. dp%ix_data < datum%ix_data) cycle
      call out_io (s_error$, r_name, 'DATUM: ' // tao_datum_name(datum), &
                                     'WHICH IS OF TYPE EXPRESSION:' // datum%data_type(12:), &
                                     'THE EXPRESSION HAS A COMPONENT: ' // datum%stack(i)%name, &
                                     'AND THIS COMPONENT IS EVALUATED AFTER THE EXPRESSION!')
      return
    enddo
    valid_value = .true.
  endif

case ('face1.offset')

case ('face2.offset')

case ('floor.x')
  call load_it (branch%ele(:)%floor%x, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('floor.y')
  call load_it (branch%ele(:)%floor%y, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('floor.z')
  call load_it (branch%ele(:)%floor%z, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('floor.theta')
  call load_it (branch%ele(:)%floor%theta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('floor.phi')
  call load_it (branch%ele(:)%floor%phi, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('gamma.a')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%a%gamma, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%a%gamma, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif
  
case ('gamma.b')
  if (data_source == 'lattice') then
    call load_it (branch%ele(:)%b%gamma, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  elseif (data_source == 'beam') then
    call load_it (bunch_params(:)%b%gamma, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
    valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!'
  endif

case ('gamma.z')
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%z%gamma, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('i5a_e6')
  if (data_source == 'beam') return
  if (ix_start > 0 .or. ix_ele > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5a_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      return
    endif
    ix_start = max(1, ix_start)
    if (ix_ele < 1) ix_ele = branch%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5a_e6(ix_start:ix_ele))
  else
    datum_value = tao_lat%modes%lin%i5a_e6
  endif
  valid_value = .true.

case ('i5b_e6')
  if (data_source == 'beam') return
  if (ix_start > 0 .or. ix_ele > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5b_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      return
    endif
    ix_start = max(1, ix_start)
    if (ix_ele < 1) ix_ele = branch%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5b_e6(ix_start:ix_ele))
  else
    datum_value = tao_lat%modes%lin%i5b_e6
  endif
  valid_value = .true.

case ('k.11b')
  if (data_source == 'beam') return
  call load_it (cc%k_11a, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
case ('k.12a')
  if (data_source == 'beam') return
  call load_it (cc%k_12a, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
case ('k.12b')
  if (data_source == 'beam') return
  call load_it (cc%k_12b, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
case ('k.22a')
  if (data_source == 'beam') return
  call load_it (cc%k_22b, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

case ('momentum_compaction')
  if (data_source == 'beam') return
  if (loc_ele%ix_branch /= 0 .or. loc_ele%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MATRIX CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(datum%ix_branch)%ele(ix_ref)
  endif

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = (/ ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap /)
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  if (ix_start < 0) ix_start = ix_ele
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix_ref, ix_start)

  call re_allocate2 (value_vec, 0, branch%n_ele_track, .false.)
  do i = ix_start, ix_ele
    value_vec(i) = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call load_it (value_vec, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('n_particle_loss')
  if (data_source /= 'beam') return
  if (ix_start == -1) then
    datum_value = uni_ele(ix_ele)%n_particle_lost_here
  else
    datum_value = sum(uni_ele(ix_start:ix_ele)%n_particle_lost_here)
  endif
  valid_value = .true.

case ('orbit.x')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(1), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.y')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(3), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.z')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(5), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.px')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(2), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.py')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(4), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.pz')
  if (data_source == 'beam') return ! bad
  call load_it (orbit(:)%vec(6), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.amp_a')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_a, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.amp_b')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_b, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.norm_amp_a')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_na, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('orbit.norm_amp_b')
  if (data_source == 'beam') return ! bad
  call load_it (cc%amp_nb, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)
  if (lat%param%ix_lost /= not_lost$ .and. ix_ele > lat%param%ix_lost) then
    valid_value = .false.
    why_invalid = 'Particle lost.'
  endif

case ('periodic.tt.')
  if (data_source == 'beam') return
  if (lat%param%lattice_type /= circular_lattice$) then
    call out_io (s_fatal$, r_name, 'LATTICE MUST BE CIRCULAR FOR A DATUM LIKE: ' // &
                                                                        datum%data_type)
    return
  endif
  
  call transfer_map_calc (lat, taylor, ix_ele, ix_ele, one_turn = .true.)
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

case ('phase.a', 'phase_frac.a')
  if (data_source == 'beam') return ! bad
  if (ix_ref < 0) then
    datum_value = ele%a%phi
  else
    datum_value = ele%a%phi - ele_ref%a%phi
    if (ix_ref > ix_ele) datum_value = datum_value - branch%ele(0)%a%phi + branch%ele(n_track)%a%phi 
  if (data_type == 'phase_frac.a') datum_value = modulo2(datum_value, pi)
  endif
  valid_value = .true.

case ('phase.b', 'phase_frac.b')
  if (data_source == 'beam') return ! bad
  if (ix_ref < 0) then
    datum_value = ele%b%phi
  else
    datum_value = ele%b%phi - ele_ref%b%phi
    if (ix_ref > ix_ele) datum_value = datum_value - branch%ele(0)%b%phi + branch%ele(n_track)%b%phi 
  endif
  if (data_type == 'phase_frac.b') datum_value = modulo2(datum_value, pi)
  valid_value = .true.

case ('phase_frac_diff')
  if (data_source == 'beam') return ! bad
  if (ix_ref < 0) then
    px = ele%a%phi 
    py = ele%b%phi 
  else
    px = ele%a%phi - ele_ref%a%phi
    py = ele%b%phi - ele_ref%b%phi
    if (ix_ref > ix_ele) px = px - branch%ele(0)%a%phi + branch%ele(n_track)%a%phi 
    if (ix_ref > ix_ele) py = py - branch%ele(0)%b%phi + branch%ele(n_track)%b%phi 
  endif
  datum_value = modulo2 (px - py, pi)
  valid_value = .true.

case ('r.')
  if (loc_ele%ix_branch /= 0 .or. loc_ele%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MATRIX CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif
  if (data_source == 'beam') return
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return

  if (ix_start < 0) ix_start = ix_ele
  if (ix_ref < 0) ix_ref = 0
  call transfer_matrix_calc (lat, .true., mat6, vec0, ix_ref, ix_start)

  call re_allocate2 (value_vec, 0, branch%n_ele_track, .false.)
  do k = ix_start, ix_ele
    value_vec(k) = mat6(i, j)
    if (k /= ix_ele) mat6 = matmul(branch%ele(k+1)%mat6, mat6)
  enddo
  call load_it (value_vec, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)

case ('rel_floor.x', 'rel_floor.y', 'rel_floor.z')

  if (ix_ref < 0) then
    ele_ref => lat%branch(datum%ix_branch)%ele(ix_ref)
  endif

  call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)
  vec3 = (/ ele%floor%x - ele_ref%floor%x, ele%floor%y - ele_ref%floor%y, ele%floor%z - ele_ref%floor%z /)
  vec3 = matmul (w0_mat, vec3)
  select case (data_type)
  case ('rel_floor.x')
    datum_value = vec3(1)
  case ('rel_floor.y')
    datum_value = vec3(2)
  case ('rel_floor.z')
    datum_value = vec3(3)
  end select
  valid_value = .true.

case ('rel_floor.theta', 'rel_floor.phi', 'rel_floor.psi')

  if (ix_ref < 0) then
    ele_ref => lat%branch(datum%ix_branch)%ele(ix_ref)
  endif

  call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)
  call floor_angles_to_w_mat (ele%floor%theta, ele%floor%phi, ele%floor%psi, w_mat)
  w_mat = matmul (w0_mat, w_mat)
  call floor_w_mat_to_angles (w_mat, 0.0_rp, theta, phi, psi)

  select case (data_type)
  case ('rel_floor.theta')
    datum_value = theta
  case ('rel_floor.phi')
    datum_value = phi
  case ('rel_floor.psi')
    datum_value = psi
  end select
  valid_value = .true.

case ('sigma.x')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s11$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)

case ('sigma.px')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s22$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)
  
case ('sigma.y')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s33$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)
  
case ('sigma.py')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s44$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)
  
case ('sigma.z')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s55$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)
  
case ('sigma.pz')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s66$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  datum_value = sqrt(datum_value)
  
case ('sigma.xy')  
  if (data_source == 'lattice') return
  call load_it (bunch_params(:)%sigma(s13$), loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  
case ('spin.theta')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%spin%theta, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  else
    call spinor_to_polar (orbit(ix_ele), polar)
    datum_value = polar%theta
  endif
  valid_value = .true.
  
case ('spin.phi')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%spin%phi, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  else
    call spinor_to_polar (orbit(ix_ele), polar)
    datum_value = polar%phi
  endif
  valid_value = .true.
  
case ('spin.polarity')
  if (data_source == 'beam') then
    call load_it (bunch_params(:)%spin%polarization, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid)
  else
    datum_value = 1.0
  endif
  valid_value = .true.
  
case ('s_position') 
  if (data_source == 'beam') return
  if (ix_ref >= 0) then
    datum_value = ele%s - ele_ref%s
  else
    datum_value = ele%s 
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
  if (loc_ele%ix_branch /= 0 .or. loc_ele%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MAP CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif
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

  if (ix_ref < 0) ix_ref = 0

  if (tao_com%ix_ref_taylor /= ix_ref .or. tao_com%ix_ele_taylor /= ix_ele) then
    if (tao_com%ix_ref_taylor == ix_ref .and. ix_ele > tao_com%ix_ele_taylor) then
      call transfer_map_calc (lat, taylor, tao_com%ix_ele_taylor, ix_ele, unit_start = .false.)
    else
      call transfer_map_calc (lat, taylor, ix_ref, ix_ele)
    endif
    tao_com%ix_ref_taylor = ix_ref
    tao_com%ix_ele_taylor = ix_ele
  endif
  datum_value = taylor_coef (taylor(i), expnt)
  valid_value = .true.

case ('unstable_orbit')
  if (lat%param%lattice_type /= linear_lattice$) return
  if (datum%ele_name == '') ix_ele = branch%n_ele_track
  if (data_source == 'beam') then
    do i = 1, ix_ele
      datum_value = datum_value + (1 + ix_ele - i) * uni_ele(i)%n_particle_lost_here
    enddo
    datum_value = datum_value / size(u%current_beam%bunch(s%global%bunch_to_plot)%particle)
  else
    iz = lat%param%ix_lost
    if (iz /= not_lost$) datum_value = max(0, 1 + ix_ele - iz)
  endif
  valid_value = .true.

case ('unstable_ring')
  if (data_source == 'beam') return
  datum_value = lat%param%growth_rate
  ! unstable_penalty is needed since at the meta stable borderline the growth rate is zero.
  if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty
  valid_value = .true.

case ('wall')
  if (data_source == 'beam') return
  print *, 'NOT YET IMPLEMENTED...'
  return

case ('wire.')  
  if (data_source == 'lattice') return
  read (data_type(6:), '(a)') angle
  datum_value = tao_do_wire_scan (ele, angle, u%current_beam)
  valid_value = .true.
  
case default
  call out_io (s_error$, r_name, 'UNKNOWN DATUM TYPE: ' // datum%data_type)
  return

end select

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine load_it (vec, loc_ref, loc_start, loc_ele, datum_value, valid_value, datum, lat, why_invalid, orbit)

implicit none

type (tao_data_struct) datum
type (lat_struct) lat
type (ele_struct), pointer :: ele, ele_start
type (lat_ele_loc_struct) loc_ref, loc_start, loc_ele
type (coord_struct), optional :: orbit(0:)

real(rp), target :: vec(0:)
real(rp) datum_value, ref_value
real(rp), pointer :: vec_ptr(:)

character(20) :: r_name = 'tao_evaluate_a_datum'
character(*), optional :: why_invalid

integer ix_m, i, n_track, ix_m2, ix_ref, ix_start, ix_ele

logical valid_value

!

ele_start => pointer_to_ele (lat, loc_start)
ele => pointer_to_ele (lat, loc_ele)
n_track = lat%branch(datum%ix_branch)%n_ele_track
ix_start = loc_start%ix_ele
ix_ref   = loc_ref%ix_ele
ix_ele   = loc_ele%ix_ele

if (ix_ele < ix_start .and. lat%param%lattice_type == linear_lattice$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // ele_start%name, &
                'IS AFTER ENDING ELEMENT: ' // ele%name)
  if (present(why_invalid)) why_invalid = 'ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum)
  valid_value = .false.
  return
endif

! Set up refrence value

if (ix_ref > -1) then
  if (present(orbit)) call data_calc (ix_ref, datum, lat, orbit)
  ref_value = vec(ix_ref)
else
  ref_value = 0
endif

 
! If ele_start does not exist

if (datum%ele_start_name == '' .or. ix_start == ix_ele) then
  if (present(orbit)) call data_calc (ix_ele, datum, lat, orbit)
  datum_value = vec(ix_ele) - ref_value
  if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(loc_ele%ix_ele))
  valid_value = .true.
  return
endif

! Set up the vector of values with the reference subtracted off

if (ref_value == 0) then
  vec_ptr => vec
else
  allocate(vec_ptr(0:ubound(vec,1)))
  vec_ptr = vec - ref_value
endif

! If there is a range

if (ix_ele < ix_start) then   ! wrap around

  if (present(orbit)) then
    do i = ix_start, n_track
      call data_calc (i, datum, lat, orbit)
    enddo
    do i = 0, ix_ele
      call data_calc (i, datum, lat, orbit)
    enddo
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = minloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) < vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = maxloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) > vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = minloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) < abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = maxloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) > abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))

  case ('int_min')
    datum_value = 0; ix_m = -1
    call integrate_min (ix_ele, n_track, datum_value, ix_m, lat, vec_ptr, datum)
    call integrate_min (0, ix_start, datum_value, ix_m, lat, vec_ptr, datum)

  case ('int_max')
    datum_value = 0; ix_m = -1
    call integrate_max (ix_ele, n_track, datum_value, ix_m, lat, vec_ptr, datum)
    call integrate_max (0, ix_start, datum_value, ix_m, lat, vec_ptr, datum)

  case default
    call out_io (s_abort$, r_name, 'BAD MERIT_TYPE: ' // datum%merit_type, &
                                 'FOR DATUM: ' // tao_datum_name(datum))
    return
  end select

else
  if (present(orbit)) then
    do i = ix_start, ix_ele
      call data_calc (i, datum, lat, orbit)
    enddo
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))

  case ('int_min')
    datum_value = 0; ix_m = -1
    call integrate_min (ix_start, ix_ele, datum_value, ix_m, lat, vec_ptr, datum)

  case ('int_max')
    datum_value = 0; ix_m = -1
    call integrate_max (ix_start, ix_ele, datum_value, ix_m, lat, vec_ptr, datum)

  case default
    call out_io (s_abort$, r_name, &
                  'SINCE THIS DATUM: ' // tao_datum_name(datum), &
                  'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                  'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
    return
  end select

endif

datum%ix_ele_merit = ix_m
valid_value = .true.
if (ref_value /= 0) deallocate (vec_ptr)

!

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_min (ix_start, ix_ele, datum_value, ix_m, lat, vec, datum)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) < vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
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

subroutine integrate_max (ix_start, ix_ele, datum_value, ix_m, lat, vec, datum)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) > vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
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

subroutine data_calc (ixd, datum, lat, orbit)

implicit none

type (ele_struct), pointer :: ele
type (this_array_struct), pointer :: cc_p
type (tao_data_struct) datum
type (lat_struct) lat
type (coord_struct) orbit(0:)

integer ie, ixd
real(rp) f, f1, f2

!

cc_p => cc(ixd)
ie = datum%ix_ele
ele => lat%ele(ie)

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
  call orbit_amplitude_calc (ele, orbit(ie), cc_p%amp_a, cc_p%amp_b, &
                           cc_p%amp_na, cc_p%amp_nb, lat%param%particle)
  cc_p%amp_calc_done = .true.
end select

end subroutine data_calc

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_do_wire_scan (ele, wire_params, theta, beam) result (moment)
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
!   datum  -- Tao_data_struct: Used for error messages and gives branch index.
!   valid  -- Logical: Set False if element does not have a definite location.
!               Set True otherwise
!   loc    -- Lat_ele_loc_struct: Location of element in the tracking part of the latticed.
!-

function tao_valid_datum_index (lat, ele_name, ix_ele, datum, valid, why_invalid) result (loc)

implicit none

type (lat_struct) lat
type (tao_data_struct) datum
type (ele_struct), pointer :: ele
type (lat_ele_loc_struct) :: loc

integer ix_ele, ixc, n_track, n_max

logical valid

character(*) ele_name
character(*), optional :: why_invalid
character(40) :: r_name = 'tao_valid_datum_index'

!

valid = .true.

loc%ix_ele = ix_ele
loc%ix_branch = datum%ix_branch

if (ele_name == '') return

n_track = lat%branch(datum%ix_branch)%n_ele_track
n_max   = lat%branch(datum%ix_branch)%n_ele_max

if (ix_ele < 0 .or. ix_ele > n_max) then
  call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE! \i5\ ', ix_ele)
  if (present(why_invalid)) why_invalid = 'ELEMENT INDEX OUT OF RANGE FOR' // tao_datum_name(datum)
  valid = .false.
  return
endif

if (datum%data_type(1:14) == 'element_param.') return
if (ix_ele <= n_track) return

ele => pointer_to_ele (lat, ix_ele, datum%ix_branch)

if (ele%lord_status == super_lord$) then
  ixc = ele%ix2_slave
  loc%ix_ele    = lat%control(ixc)%ix_slave
  loc%ix_branch = lat%control(ixc)%ix_branch
  return
endif

valid = .false.
call out_io (s_error$, r_name, &
            'ELEMENT: ' // trim(lat%ele(ix_ele)%name) // &
            '    WHICH IS A: ' // control_name(lat%ele(ix_ele)%lord_status), &
            'CANNOT BE USED IN DEFINING A DATUM SINCE IT DOES NOT HAVE ', &
            '   A DEFINITE LOCATION IN THE LATTICE.', &
            'FOR DATUM: ' // tao_datum_name(datum) )
if (present(why_invalid)) why_invalid = 'NO DEFINITE LOCATION IN LATTICE FOR: ' // tao_datum_name(datum)

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real (expression, value, err_flag, use_good_user, good, stack)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression    -- character(*): arithmetic expression
!   use_good_user -- Logical, optional: Use the good_user logical in evaluating good(:)
!                      Default is False.
!   print_err     -- Logical, optional: Print an error message? Default is True.
!
! Output:
!   value        -- real(rp): Value of arithmetic expression.
!   err_flag     -- Logical: TRUE on error.
!   good         -- Logical, optional: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   stack(:)  -- Tao_eval_stack1_struct, optional: Evaluation stack for the
!                  expression. This is useful to save if the same expression is
!                  to be evaluated repeatedly. 
!                  With this, tao_evaluate_stack can be called directly.
!-

subroutine tao_to_real (expression, value, err_flag, use_good_user, good, stack, print_err)

implicit none

type (tao_eval_stack1_struct), allocatable, optional :: stack(:)

character(*) :: expression

real(rp) value
real(rp), allocatable, save :: vec(:)

logical, allocatable, save :: ok(:)
logical err_flag, good_user
logical, optional :: good, use_good_user
logical, optional :: print_err

!

good_user = logic_option(.false., use_good_user)
call tao_evaluate_expression (expression, 1, good_user, &
                                            vec, ok, err_flag, print_err, stack)
if (err_flag) return
value = vec(1)
if (present(good)) good = ok(1)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real_vector (expression, n_size, use_good_user, value, good, 
!                                                                err_flag, print_err)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression    -- Character(*): Arithmetic expression.
!   n_size        -- Integer: Size of the value array. If the expression
!                               is a scaler then the value will be spread.
!                               If n_size = 0 then the natural size determined 
!                               by expression is used.
!   use_good_user -- Logical: Use the good_user logical in evaluating good(:)
!   print_err     -- Logical, optional: Print an error message? Default is True.
!  
! Output:
!   value(:)     -- Real(rp), allocatable: Value of arithmetic expression.
!   good(:)      -- Logical, allocatable: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   err_flag     -- Logical: True on error. False otherwise
!-

subroutine tao_to_real_vector (expression, n_size, use_good_user, value, good, &
                                                                      err_flag, print_err)

use random_mod

implicit none

real(rp), allocatable :: value(:)
logical, allocatable :: good(:)

integer n_size

character(*) :: expression
character(16) :: r_name = "tao_to_real_vector"

logical err_flag, err, use_good_user
logical, optional :: print_err

!

call tao_evaluate_expression (expression, n_size, use_good_user, &
                                                       value, good, err_flag, print_err)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_expression (expression, n_size, use_good_user, &
!                                   value, good, err_flag, print_err, stack)
!
! Mathematically evaluates a character expression.
!
! Input:
!   expression    -- Character(*): Arithmetic expression.
!   n_size        -- Integer: Size of the value array. If the expression
!                      is a scaler then the value will be spread.
!                      If n_size = 0, the natural size is determined by expression is used.
!   use_good_user -- Logical: Use the good_user logical in evaluating good(:)
!   print_err     -- Logical, optional: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!   
! Output:
!   value(:)  -- Real(rp), allocatable: Value of arithmetic expression.
!   good(:)   -- Logical, allocatable: Is the value valid? 
!                  Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                  orbit.x[23]|good_user is False.
!   err_flag  -- Logical: True on an error. EG: Invalid expression.
!                  A divide by zero is not an error but good(:) will be set to False.
!   stack(:)  -- Tao_eval_stack1_struct, optional: Evaluation stack for the
!                  expression. This is useful to save if the same expression is
!                  to be evaluated repeatedly. 
!                  With this, tao_evaluate_stack can be called directly.
!-

subroutine tao_evaluate_expression (expression, n_size, use_good_user, &
                                    value, good, err_flag, print_err, stack)

use random_mod

implicit none

type (tao_eval_stack1_struct), save :: stk(100)
type (tao_eval_stack1_struct), allocatable, optional :: stack(:)

integer i_lev, i_op, i, ios, n, n_size, n__size
integer op(200), ix_word, i_delim, i2, ix, ix_word2, ixb

real(rp), allocatable :: value(:)

character(*) :: expression
character(200) phrase
character(1) delim
character(40) word, word2
character(40) :: r_name = "tao_evaluate_expression"
character(40) saved_prefix

logical, allocatable :: good(:)
logical delim_found, split, ran_function_pending, use_good_user
logical err_flag, err, wild, printit
logical, optional :: print_err

! Don't destroy the input expression

err_flag = .true.
saved_prefix = ''
printit = logic_option(.true., print_err)

phrase = expression

! if phrase is blank then return 0.0

call string_trim (phrase, phrase, ios)
if (ios == 0) then
  call out_io (s_warn$, r_name, &
    "Expression is blank", len(phrase))
  value = 0.0
  return
endif
 
! General idea: Create a reverse polish stack that represents the expression.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

i_lev = 0
i_op = 0
ran_function_pending = .false.

do i = 1, size(stk)
  stk(i)%name = ''
  if (allocated(stk(i)%good)) deallocate (stk(i)%good)
  if (allocated(stk(i)%value_ptr)) deallocate (stk(i)%value_ptr)
enddo

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  call word_read (phrase, '+-*/()^,}[ ', word, ix_word, delim, &
                    delim_found, phrase)

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    call out_io (s_warn$, r_name, 'RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT')
    return
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
  ! get split at the "-" even though "-" is a delimiter

  split = .true.         ! assume initially that we have a split number
  if (ix_word == 0) then
    split = .false.
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'e' ) then
    split = .false.
  endif
  if (delim /= '-' .and. delim /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) == 0) split = .false.
  enddo

  ! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, &
                    delim_found, phrase)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

  ! Something like "lcav[lr(2).freq]" or "[2,4]@orbit.x[1,4] will get split on the "["

  do
    if (delim /= '[') exit

    call word_read (phrase, ']', word2, ix_word2, delim, delim_found, phrase)
    if (.not. delim_found) then
      call out_io (s_warn$, r_name, "NO MATCHING ']' FOR OPENING '[':" // expression)
      return
    endif
    word = word(:ix_word) // '[' // trim(word2) // ']'
    ix_word = ix_word + ix_word2 + 2
    if (phrase == ' ') then  
      delim_found = .false.
      delim = ' '
    elseif (phrase(1:1) == ' ') then  
      call string_trim (phrase, phrase, ix)
      if (index('+-*/()^,}[', phrase(1:1)) == 0) then
        delim = ' '
      else
        delim = phrase(1:1)
        phrase = phrase(2:)
      endif
    else          ! even more...
      call word_read (phrase, '[+-*/()^,}', word2, ix_word2, delim, &
                                                  delim_found, phrase)
      word = word(:ix_word) // trim(word2)       
      ix_word = ix_word + ix_word2 
    endif

  enddo

  ! If delim = "*" then see if this is being used as a wildcard
  ! Examples: "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|"
  ! If so, we have split in the wrong place and we need to correct this.

  if (delim == '*') then

    wild = .false.

    select case (phrase(1:1))
    case ( ']', '[', '|', '@')
      wild = .true.
    case ('.')
      if (index('0123456789', phrase(2:2)) /= 0) wild = .true.
    end select
    ixb = index(phrase, '|')
    if (ixb == 0) wild = .false.

    if (wild) then
      word = word(:ix_word) // '*' // phrase(1:ixb) 
      phrase = phrase(ixb+1:)
      call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, &
                                                delim_found, phrase)
      word = trim(word) // trim(word2)       
      ix_word = len_trim(word)
    endif

  endif

  !---------------------------
  ! Now see what we got...

  ! For a word ending in '|' then must be a construct like 'orbit.x|-model'.
  ! So store the 'orbit.x|' prefix

  if (ix_word /= 0) then
    if (word(ix_word:ix_word) == '|') then
      saved_prefix = word
      word = ''
      ix_word = 0
    endif
  endif

  ! For a "(" delim we must have a function

  if (delim == '(') then

    ran_function_pending = .false.
    if (ix_word /= 0) then
      word2 = word
      call downcase_string (word2)
      select case (word2)
      case ('sin')
        call pushit (op, i_op, sin$)
      case ('cos')
        call pushit (op, i_op, cos$)
      case ('tan') 
        call pushit (op, i_op, tan$)
      case ('asin') 
        call pushit (op, i_op, asin$)
      case ('acos') 
        call pushit (op, i_op, acos$)
      case ('atan') 
        call pushit (op, i_op, atan$)
      case ('abs') 
        call pushit (op, i_op, abs$)
      case ('sqrt') 
        call pushit (op, i_op, sqrt$)
      case ('log') 
        call pushit (op, i_op, log$)
      case ('exp') 
        call pushit (op, i_op, exp$)
      case ('ran') 
        call pushit (op, i_op, ran$)
        ran_function_pending = .true.
      case ('ran_gauss') 
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
      case default
        call out_io (s_warn$, r_name, 'UNEXPECTED CHARACTERS BEFORE "(": ', &
                                      'IN EXPRESSION: ' // expression)
        return
      end select
    endif

    call pushit (op, i_op, l_parens$)
    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call pushit (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (.not. ran_function_pending) then
        call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING BEFORE ")"', &
                                      'IN EXPRESSION: ' // expression)
        return
      endif
    else
      call pushit (stk%type, i_lev, numeric$)
      call tao_param_value_routine (word, saved_prefix, stk(i_lev), err, printit)
      if (err) then
        if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EVALUATING EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
        return
      endif
    endif

    ran_function_pending = .false.

    do
      do i = i_op, 1, -1       ! release pending ops
        if (op(i) == l_parens$) exit            ! break do loop
        call pushit (stk%type, i_lev, op(i))
      enddo

      if (i == 0) then
        call out_io (s_warn$, r_name, 'UNMATCHED ")" IN EXPRESSION: ' // expression)
        return
      endif

      i_op = i - 1

      call word_read (phrase, '+-*/()^,}', word, ix_word, delim, &
                    delim_found, phrase)
      if (ix_word /= 0) then
        call out_io (s_warn$, r_name, &
                'UNEXPECTED CHARACTERS AFTER ")" IN EXPRESSION: ' // expression)
        return
      endif

      if (delim /= ')') exit    ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call out_io (s_warn$, r_name, &
          '")(" CONSTRUCT DOES NOT MAKE SENSE IN EXPRESSION: ' // expression)
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call out_io (s_warn$, r_name, &
            'CONSTANT OR VARIABLE MISSING IN EXPRESSION: ' // expression)
      return
    endif
    call pushit (stk%type, i_lev, numeric$)
    call tao_param_value_routine (word, saved_prefix, stk(i_lev), err, printit)
    if (err) then
      if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
      return
    endif
  endif

  ! If we are here then we have an operation that is waiting to be identified


  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case (')')
    i_delim = r_parens$
  case ('^')
    i_delim = power$
  case (',', '}')
    i_delim = no_delim$
    call out_io (s_error$, r_name, &
                      'DELIMITOR FOUND OUT OF PLACE: ' // delim, &
                      'IN EXPRESSION: ' // expression)

      return
  case default
    if (delim_found) then
      call out_io (s_error$, r_name, 'INTERNAL ERROR')
      call err_exit
    endif
    i_delim = no_delim$
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (eval_level(op(i)) >= eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        call out_io (s_warn$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
        return
      endif
      call pushit (stk%type, i_lev, op(i))
    else
      exit
    endif
  enddo

  ! put the pending operation on the OP stack

  i_op = i
  if (i_delim == no_delim$) then
    exit parsing_loop
  else
    call pushit (op, i_op, i_delim)
  endif

enddo parsing_loop

!------------------------------------------------------------------
! Some error checks

if (i_op /= 0) then
  call out_io (s_warn$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
  return
endif

if (i_lev == 0) then
  call out_io (s_warn$, r_name, 'NO VALUE FOUND IN EXPRESSION: ' // expression)
  return
endif

n__size = 1
do i = 1, i_lev
  if (stk(i)%type /= numeric$ .and. stk(i)%type /= var$) cycle
  n = size(stk(i)%value)
  if (n == 1) cycle
  if (n__size == 1) n__size = n
  if (n /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
    return
  endif
enddo

if (n_size /= 0) then
  if (n__size /= 1 .and. n_size /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
    return
  endif
  n__size = n_size
endif

call tao_evaluate_stack (stk(1:i_lev), n__size, use_good_user, value, good, err_flag, printit)

! If the stack argument is present then copy stk to stack

if (present(stack)) then
  if (allocated(stack)) deallocate(stack)
  allocate (stack(i_lev))
  do i = 1, i_lev
    if (allocated (stk(i)%value)) then
      n = size(stk(i)%value)
      allocate (stack(i)%value(n), stack(i)%good(n))
      if (allocated (stack(i)%value_ptr)) allocate (stack(i)%value_ptr(n))
    endif
    stack(i) = stk(i)
  enddo
endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine pushit (stack, i_lev, value)

implicit none

integer stack(:), i_lev, value

character(6) :: r_name = "pushit"

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call out_io (s_warn$, r_name, 'STACK OVERFLOW.')
  call err_exit
endif

stack(i_lev) = value

end subroutine pushit
                       
end subroutine tao_evaluate_expression

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine tao_param_value_routine (str, saved_prefix, stack, err_flag, print_err)

implicit none

type (tao_eval_stack1_struct) stack
type (tao_real_pointer_struct), allocatable, save :: re_array(:)
type (tao_integer_array_struct), allocatable, save :: int_array(:)
type (tao_data_struct) datum

integer ios, i, n, ix, ix_ele, ix2

character(*) str, saved_prefix
character(60) name
character(40) :: r_name = 'tao_param_value_routine'

logical err_flag, print_err

! pi, etc

err_flag = .false.

select case (str)
case ('pi')
  call re_allocate(stack%value, 1)
  stack%value(1) = pi  
  return
case ('twopi')
  call re_allocate(stack%value, 1)
  stack%value(1) = twopi  
  return
case ('fourpi')
  call re_allocate(stack%value, 1)
  stack%value(1) = 4*pi  
  return
case ('sqrt_2')
  call re_allocate(stack%value, 1)
  stack%value(1) = sqrt(2.0_rp)  
  return
end select

! Case where str represents a number.

if (is_real(str)) then
  call re_allocate(stack%value, 1)
  read (str, *, iostat = ios) stack%value(1)
  if (ios /= 0) then
    if (print_err) call out_io (s_warn$, r_name, "This doesn't seem to be a number: " // str)
    err_flag = .true.
  endif
  return
endif

! Case where str is a variable name.
! Remember the last string so 'orbit.x|meas-ref' translates to 'orbit.x|meas - orbit.x|ref.'

ix = index(str, '|')
ix_ele = index(str, ':')
stack%name = str
if (ix == 0 .and. ix_ele == 0) stack%name = trim(saved_prefix) // str
if (ix /= 0) saved_prefix = str(1:ix)
name = stack%name

if (.not. allocated(re_array)) allocate (re_array(0))

! Look for a datum-on-the-fly.

if (name(1:4) == 'dat:' .or. index(name, '@dat:') /= 0) then
  call tao_evaluate_data_on_the_fly (err_flag, name, stack%value, print_err = .false.)
  stack%type = var_on_the_fly$
  return

! Look for a lattice element parameter 

elseif (name(1:4) == 'ele:' .or. index(name, '@ele:') /= 0) then
  call tao_evaluate_element_parameters (err_flag, name, stack%value, print_err = .false.)
  stack%type = var_on_the_fly$
  return

! Look for variable or data values

else
  call tao_find_var (err_flag, name, re_array = re_array, print_err = .false.)
  if (err_flag) call tao_find_data (err_flag, name, &
                     re_array = re_array, int_array = int_array, print_err = .false.)

  if (err_flag) return
endif

! Now transfer the information to the stack

if (size(re_array) /= 0) then
  stack%type = var$
  n = size(re_array)
  if (allocated(stack%value_ptr)) then
    if (size(stack%value_ptr) /= n) deallocate (stack%value_ptr)
  endif
  if (.not. allocated(stack%value_ptr)) allocate (stack%value_ptr(n))
  call re_allocate (stack%value, n)
  call re_allocate (stack%good, n)
  do i = 1, n
    stack%value(i) = re_array(i)%r
    stack%value_ptr(i)%r => re_array(i)%r
    ! good is only used with data and not variables
    if (associated(re_array(i)%good_value)) then
      stack%value_ptr(i)%good_user => re_array(i)%good_user
      stack%value_ptr(i)%good_value => re_array(i)%good_value
    else
      stack%good(i) = .true.
    endif
  enddo

elseif (size(int_array) /= 0) then
  n = size(int_array)
  call re_allocate (stack%value, n)
  call re_allocate (stack%good, n)
  do i = 1, n
    stack%value(i) = int_array(i)%i
    stack%good(i)  = .true.
  enddo

else
  if (print_err) call out_io (s_warn$, r_name, &
               "This doesn't seem to be datum value or variable value: " // name)
  err_flag = .true.
  return
endif

end subroutine tao_param_value_routine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_stack (stack, n_size, use_good_user, value, good, 
!                                                                err_flag, print_err)
!
! Routine to evaluate an expression stack.
!
! Input:
!   stack(:)      -- Tao_eval_stack1_struct: Expression stack
!   n_size        -- Integer: Result array size.
!   use_good_user -- Logical: Use the good_user logical in evaluating good(:)
!   print_err     -- Logical: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!
! Output:
!   value(:)     -- Real(rp), allocatable: Value of arithmetic expression.
!   good(:)      -- Logical, allocatable: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   err_flag     -- Logical: True on error. False otherwise
!-

subroutine tao_evaluate_stack (stack, n_size, use_good_user, value, good, err_flag, print_err)

implicit none

type (tao_eval_stack1_struct), target :: stack(:)
type (tao_eval_stack1_struct), pointer :: s(:)
type (tao_eval_stack1_struct) stk2(20)

real(rp), allocatable :: value(:)

integer i, i2, j, n
integer n_size

logical, allocatable :: good(:)
logical err_flag, use_good_user, print_err

character(20) :: r_name = 'tao_evaluate_stack'

! Calculate good

s => stack   ! For debugging purposes
err_flag = .true.

call re_allocate (good, n_size)
call re_allocate (value, n_size)

good = .true.
do i = 1, size(stack)
  if (allocated(stack(i)%value_ptr)) then
    if (associated(stack(i)%value_ptr(1)%good_value)) then    
      do j = 1, size(stack(i)%value_ptr)
        if (use_good_user) then
          stack(i)%good(j) = stack(i)%value_ptr(j)%good_value .and. stack(i)%value_ptr(j)%good_user
        else
          stack(i)%good(j) = stack(i)%value_ptr(j)%good_value
        endif
      enddo
    endif
  endif
  if (.not. allocated(stack(i)%good)) cycle
  if (size(stack(i)%good) == 1) then; good = good .and. stack(i)%good(1:1)
  else;                               good = good .and. stack(i)%good
  endif
enddo

! Go through the stack and perform the operations...

i2 = 0  ! stack pointer
do i = 1, size(stack)

  select case (stack(i)%type)
  case (numeric$) 
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (var_on_the_fly$)
    call tao_param_value_routine (stack(i)%name, '', stack(i), err_flag, print_err)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (var$)
    do j = 1, size(stack(i)%value)
      stack(i)%value(j) = stack(i)%value_ptr(j)%r
    enddo
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (unary_minus$) 
    stk2(i2)%value = -stk2(i2)%value

  case (unary_plus$) 
    ! Nothing to do

  case (plus$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) + stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value
    endif
    i2 = i2 - 1

  case (minus$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) - stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value
    endif
    i2 = i2 - 1

  case (times$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) * stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value
    endif
    i2 = i2 - 1

  case (divide$) 
    n = size(stk2(i2)%value)
    do j = 1, n
      if (stk2(i2)%value(j) == 0) then  ! Divide by 0 error!
        stk2(i2)%value(j) = 1
        if (n == 1) then
          good = .false.  ! All are false
        else
          good(j) = .false.
        endif
      endif
    enddo

    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) / stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value
    endif
    i2 = i2 - 1

  case (power$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) ** stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value
    endif
    i2 = i2 - 1

  case (sin$) 
    stk2(i2)%value = sin(stk2(i2)%value)

  case (cos$) 
    stk2(i2)%value = cos(stk2(i2)%value)

  case (tan$) 
    stk2(i2)%value = tan(stk2(i2)%value)

  case (asin$) 
    stk2(i2)%value = asin(stk2(i2)%value)

  case (acos$) 
    stk2(i2)%value = acos(stk2(i2)%value)

  case (atan$) 
    stk2(i2)%value = atan(stk2(i2)%value)

  case (abs$) 
    stk2(i2)%value = abs(stk2(i2)%value)

  case (sqrt$) 
    stk2(i2)%value = sqrt(stk2(i2)%value)

  case (log$) 
    stk2(i2)%value = log(stk2(i2)%value)

  case (exp$) 
    stk2(i2)%value = exp(stk2(i2)%value)

  case (ran$) 
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_uniform(stk2(i2)%value)

  case (ran_gauss$) 
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_gauss(stk2(i2)%value)

  case default
    call out_io (s_warn$, r_name, 'INTERNAL ERROR')
    call err_exit
  end select
enddo

if (i2 /= 1) then
  call out_io (s_warn$, r_name, 'INTERNAL ERROR')
  call err_exit
endif

if (size(stk2(1)%value) == 1 .and. n_size > 1) then
  value = stk2(1)%value(1)
else
  call value_transfer (value, stk2(1)%value)
endif

where (.not. good) value = 0

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine value_transfer (to_array, from_array)

real(rp), allocatable :: to_array(:)
real(rp) from_array(:)

!

call re_allocate (to_array, size(from_array))
to_array = from_array

end subroutine

end subroutine tao_evaluate_stack 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_int (str, i_int, err)
! 
! Converts a string to an integer
!
! If the string str is blank then i_int = 0
!-

subroutine tao_to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err
character(12) :: r_name = "tao_to_int"

!

call string_trim (str, str, ios)
if (ios .eq. 0) then
  i_int = 0
  return
endif
 
err = .false.
read (str, *, iostat = ios) i_int

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
  err = .true.
  return
endif

end subroutine

end module tao_data_and_eval_mod
