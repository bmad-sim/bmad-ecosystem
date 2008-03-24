module tao_data_mod

use tao_mod
use spin_mod
use utilities_mod
use random_mod

type this_coupling_struct
  real(rp) cbar(2,2)
  real(rp) coupling11, coupling12a, coupling12b, coupling22
  real(rp) f_11, f_12a, f_12b, f_22
  logical calc_done
end type

type (this_coupling_struct), save, allocatable, target, private :: cc(:)

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

cc%calc_done = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, taylor_in)
!
! Subroutine to put the proper data in the specified datum
!
! If datum results in NaN then datum_value = tiny(1.0_rp)
!
! Input:
!   datum        -- Tao_data_struct: What type of datum
!   u            -- Tao_universe_struct: Which universe to use.
!   tao_lat      -- Tao_lattice_struct: Lattice to use.
!   taylor_in(6) -- Taylor_struct, optional: Starting point for 
!                     tt: and t: constraints.
!     
! Output:
!   datum   -- Tao_data_struct: 
!     %ix_ele_merit -- For max/min type constraints: Place where value is max/min. 
!   datum_value -- Real(rp): Value of the datum.
!   valid_value -- Logical: Set false when there is a problem. Set true otherwise.
!-

subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, taylor_in)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat
type (lat_struct), pointer :: lat
type (normal_modes_struct) mode
type (taylor_struct), optional :: taylor_in(6)
type (spin_polar_struct) polar
type (ele_struct), pointer :: ele, ele0
type (coord_struct), pointer :: orb0

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), one_pz
real(rp) gamma

integer, save :: ix_save = -1
integer i, j, k, m, n, ix, ix1, ix0, expnt(6), n_lat

character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) data_type, data_source

logical found, valid_value

! See if there is a hook for this datum

call tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
if (found) return

! Check range

data_source = datum%data_source
data_type = datum%data_type
lat => tao_lat%lat
valid_value = .true.

ix0 = datum%ix_ele0
ix1 = datum%ix_ele

if (datum%ele_name /= '') then
  ix1 = tao_valid_datum_index (lat, ix1, valid_value)
  if (.not. valid_value) return
endif

if (datum%ele0_name /= '') then
  ix0 = tao_valid_datum_index (lat, ix0, valid_value)
  if (.not. valid_value) return
endif

!

n_lat = tao_lat%lat%n_ele_track
datum_value = 0           ! default
datum%ix_ele_merit = ix1  ! default

if (data_type(1:2) == 'r.') data_type = 'r.'
if (data_type(1:2) == 't.') data_type = 't.'
if (data_type(1:3) == 'tt.') data_type = 'tt.'
if (data_type(1:5) == 'wire.') data_type = 'wire.'

if (data_type(1:4) == 'emit') call convert_total_energy_to ( &
                    lat%ele(ix1)%value(E_tot$), lat%param%particle, gamma)

if (data_source /= "lattice" .and. data_source /= "beam") then
  call out_io (s_error$, r_name, &
          'UNKNOWN DATA_SOURCE: ' // data_source, &
          'FOR DATUM: ' // tao_datum_name(datum))
  call err_exit
endif

if (data_source == 'lattice') then
  select case (data_type)
  case ('beta.z', 'alpha.z', 'gamma.z', 'emit.x', 'norm_emit.x', &
        'emit.y', 'norm_emit.y', 'emit.z', 'norm_emit.z', &
        'dpx_dx', 'dpy_dy', 'dpz_dz', 'sigma.x', 'sigma.p_x', 'sigma.y', 'sigma.p_y', &
        'sigma.z', 'sigma.p_z', 'sigma.xy', 'wire.')
    call out_io (s_error$, r_name, 'DATA_SOURCE = "lattice" NOT VALID FOR DATUM: ' // &
                                   tao_datum_name(datum))
  end select
endif

if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost .and. &
                                      data_source == 'beam') valid_value = .false.

!---------------------------------------------------

select case (data_type)

case ('orbit.x')
  call load_it (tao_lat%orb(:)%vec(1), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.
case ('orbit.y')
  call load_it (tao_lat%orb(:)%vec(3), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.
case ('orbit.z')
  call load_it (tao_lat%orb(:)%vec(5), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.

case ('orbit.p_x')
  call load_it (tao_lat%orb(:)%vec(2), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.
case ('orbit.p_y')
  call load_it (tao_lat%orb(:)%vec(4), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.
case ('orbit.p_z')
  call load_it (tao_lat%orb(:)%vec(6), ix0, ix1, datum_value, valid_value, datum, lat)
  if (lat%param%ix_lost /= not_lost$ .and. ix1 >= lat%param%ix_lost) valid_value = .false.

case ('bpm.x')
  call orbit_to_bpm_reading (tao_lat%orb(ix1), lat%ele(ix1), x_plane$, datum_value)
case ('bpm.y')
  call orbit_to_bpm_reading (tao_lat%orb(ix1), lat%ele(ix1), y_plane$, datum_value)

case ('phase.a')
  datum_value = lat%ele(ix1)%a%phi - lat%ele(ix0)%a%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%a%phi + lat%ele(n_lat)%a%phi 
case ('phase.b')
  datum_value = lat%ele(ix1)%b%phi - lat%ele(ix0)%b%phi
  if (ix0 > ix1) datum_value = datum_value - lat%ele(0)%b%phi + lat%ele(n_lat)%b%phi 

case ('phase_frac_diff')
  px = lat%ele(ix1)%a%phi - lat%ele(ix0)%a%phi
  if (ix0 > ix1) px = px - lat%ele(0)%a%phi + lat%ele(n_lat)%a%phi 
  py = lat%ele(ix1)%b%phi - lat%ele(ix0)%b%phi
  if (ix0 > ix1) py = py - lat%ele(0)%b%phi + lat%ele(n_lat)%b%phi 
  datum_value = modulo (px, twopi) - modulo (py, twopi)

case ('beta.x')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%beta
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH data_source: ' // data_source)
    call err_exit
  endif
    
case ('beta.y')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%beta
  else
    call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // data_type, &
                                   'WITH data_source: ' // data_source)
    call err_exit
  endif

case ('beta.z')
  if (data_source == "lattice") then
    valid_value = .false.
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%beta
  endif

case ('beta.a')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%a%beta, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%beta
  endif
    
case ('beta.b')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%b%beta, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%beta
  endif

case ('alpha.a')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%a%alpha, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%alpha
  endif
  
case ('alpha.b')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%b%alpha, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%alpha
  endif

case ('alpha.z')
  if (data_source == "lattice") then
    valid_value = .false.
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%alpha
  endif

case ('gamma.a')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%a%gamma, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%gamma
  endif
  
case ('gamma.b')
  if (data_source == "lattice") then
    call load_it (lat%ele(:)%b%gamma, ix0, ix1, datum_value, valid_value, datum, lat)
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%gamma
  endif

case ('gamma.z')
  if (data_source == "lattice") then
    valid_value = .false.
  elseif (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%gamma
  endif

case ('eta.x')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%eta
  else
    call load_it (lat%ele(:)%x%eta, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('eta.y')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%eta
  else
    call load_it (lat%ele(:)%y%eta, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('etap.x')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%etap
  else
    call load_it (lat%ele(:)%x%etap, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('etap.y')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%etap
  else
    call load_it (lat%ele(:)%y%etap, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('eta.a')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%eta
  else
    call load_it (lat%ele(:)%a%eta, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('eta.b')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%eta
  else
    call load_it (lat%ele(:)%b%eta, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('etap.a')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%etap
  else
    call load_it (lat%ele(:)%a%etap, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('etap.b')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%etap
  else
    call load_it (lat%ele(:)%b%etap, ix0, ix1, datum_value, valid_value, datum, lat)
  endif

case ('e_tot')
  call load_it (lat%ele(0:n_lat)%value(E_TOT$) * (1+tao_lat%orb(0:n_lat)%vec(6)), &
                                     ix0, ix1, datum_value, valid_value, datum, lat)

case ('%e_tot')
  call load_it (tao_lat%orb(0:n_lat)%vec(6), ix0, ix1, &
                                            datum_value, valid_value, datum, lat)
  
case ('coupling.11b')
  call load_it (cc%coupling11, ix0, ix1, datum_value, valid_value, &
                                     datum, lat, cc%f_11, coupling_here = .true.)
case ('coupling.12a')
  call load_it (cc%coupling12a, ix0, ix1, datum_value, valid_value, &
                                     datum, lat, cc%f_12a, coupling_here = .true.)
case ('coupling.12b')
  call load_it (cc%coupling12b, ix0, ix1, datum_value, valid_value, &
                                     datum, lat, cc%f_12b, coupling_here = .true.)
case ('coupling.22a')
  call load_it (cc%coupling22, ix0, ix1, datum_value, valid_value, &
                                     datum, lat, cc%f_22, coupling_here = .true.)

case ('cbar.11')
  call load_it (cc%cbar(1,1), ix0, ix1, datum_value, valid_value, &
                                     datum, lat, coupling_here = .true.)
case ('cbar.12')
  call load_it (cc%cbar(1,2), ix0, ix1, datum_value, valid_value, &
                                     datum, lat, coupling_here = .true.)
case ('cbar.21')
  call load_it (cc%cbar(2,1), ix0, ix1, datum_value, valid_value, &
                                     datum, lat, coupling_here = .true.)
case ('cbar.22')
  call load_it (cc%cbar(2,2), ix0, ix1, datum_value, valid_value, &
                                     datum, lat, coupling_here = .true.)

case ('i5a_e6')
  if (ix0 > 0 .or. ix1 > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5a_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      valid_value = .false.
      return
    endif
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5a_e6(ix0:ix1))
    datum%ix_ele_merit = ix1
  else
    datum_value = tao_lat%modes%lin%i5a_e6
    datum%ix_ele_merit = lat%n_ele_track
  endif

case ('i5b_e6')
  if (ix0 > 0 .or. ix1 > 0) then
    if (.not. allocated(tao_lat%rad_int%lin_i5b_e6)) then
      call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
      valid_value = .false.
      return
    endif
    ix0 = max(1, ix0)
    if (ix1 < 1) ix1 = lat%n_ele_track
    datum_value = sum(tao_lat%rad_int%lin_i5b_e6(ix0:ix1))
    datum%ix_ele_merit = ix1
  else
    datum_value = tao_lat%modes%lin%i5b_e6
    datum%ix_ele_merit = lat%n_ele_track
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
    call transfer_map_calc (lat, taylor_in, ix0, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), j, k)
  else
    if (tao_com%ix0_taylor /= ix0 .or. tao_com%ix1_taylor /= ix1) then
      call transfer_map_calc (lat, tao_com%taylor, ix0, ix1)
      tao_com%ix0_taylor = ix0
      tao_com%ix1_taylor = ix1
    endif
    datum_value = taylor_coef (tao_com%taylor(i), j, k)
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
    call transfer_map_calc (lat, taylor_in, ix0, ix1, unit_start = .false.)
    datum_value = taylor_coef (taylor_in(i), expnt)
  else
    if (tao_com%ix0_taylor /= ix0 .or. tao_com%ix1_taylor /= ix1) then
      call transfer_map_calc (lat, tao_com%taylor, ix0, ix1)
      tao_com%ix0_taylor = ix0
      tao_com%ix1_taylor = ix1
    endif
    datum_value = taylor_coef (tao_com%taylor(i), expnt)
  endif

case ('floor.x')
  if (datum%ele0_name /= '') then
    datum_value = lat%ele(ix1)%floor%x - lat%ele(ix0)%floor%x
  else
    datum_value = lat%ele(ix1)%floor%x
  endif

case ('floor.y')
  if (datum%ele0_name /= '') then
    datum_value = lat%ele(ix1)%floor%y - lat%ele(ix0)%floor%y
  else
    datum_value = lat%ele(ix1)%floor%y 
  endif

case ('floor.z')
  if (datum%ele0_name /= '') then
    datum_value = lat%ele(ix1)%floor%z - lat%ele(ix0)%floor%z
  else
    datum_value = lat%ele(ix1)%floor%z 
  endif

case ('floor.theta')
  if (datum%ele0_name /= '') then
    datum_value = lat%ele(ix1)%floor%theta - lat%ele(ix0)%floor%theta
  else
    datum_value = lat%ele(ix1)%floor%theta 
  endif

case ('s_position') 
  if (ix0 >= 0) then
    datum_value = lat%ele(ix1)%s - lat%ele(ix0)%s
  else
    datum_value = lat%ele(ix1)%s 
  endif

case ('wall')
  print *, 'NOT YET IMPLEMENTED...'
  call err_exit

!---------------------------------------------------------
! Beam Emittance
  
case ('emit.x', 'norm_emit.x')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%x%norm_emitt
  else
    valid_value = .false.
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.y', 'norm_emit.y')  
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%y%norm_emitt
  else
    valid_value = .false.
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.z', 'norm_emit.z')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%z%norm_emitt
  else
    valid_value = .false.
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

case ('emit.a', 'norm_emit.a')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%a%norm_emitt
  elseif (data_source == "lattice") then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_a)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        valid_value = .false.
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_a, &
                              ix0, ix1, datum_value, valid_value, datum, lat)
    else
      datum_value = gamma * tao_lat%modes%a%emittance  
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma
  
case ('emit.b', 'norm_emit.b')  
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%b%norm_emitt
  elseif (data_source == "lattice") then
    if (lat%param%lattice_type == linear_lattice$) then
      if (.not. allocated(tao_lat%rad_int%lin_norm_emit_b)) then
        call out_io (s_error$, r_name, 'tao_lat%rad_int not allocated')
        valid_value = .false.
        return
      endif
      call load_it (tao_lat%rad_int%lin_norm_emit_b, &
                              ix0, ix1, datum_value, valid_value, datum, lat)
    else
      datum_value = gamma * tao_lat%modes%b%emittance
    endif
  endif
  if (data_type(1:4) == 'emit') datum_value = datum_value / gamma

!---------------------------------------------------------
case ('chrom.a')
  datum_value = tao_lat%a%chrom

case ('chrom.b')
  datum_value = tao_lat%b%chrom

case ('unstable_ring')
  datum_value = lat%param%growth_rate
  ! unstable_penalty is needed since at the meta stable borderline the growth rate is zero.
  if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty

case ('dpx_dx') 
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s12$) / tao_lat%bunch_params(ix1)%sigma(s11$)
  else
    valid_value = .false.
  endif

case ('dpy_dy') 
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s34$) / tao_lat%bunch_params(ix1)%sigma(s33$)
  else
    valid_value = .false.
  endif

case ('dpz_dz') 
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s56$) / tao_lat%bunch_params(ix1)%sigma(s55$)
  else
    valid_value = .false.
  endif

case ('sigma.x')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s11$))
  else
    valid_value = .false.
  endif
  
case ('sigma.p_x')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s22$))
  else
    valid_value = .false.
  endif
  
case ('sigma.y')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s33$))
  else
    valid_value = .false.
  endif
  
case ('sigma.p_y')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s44$))
  else
    valid_value = .false.
  endif
  
case ('sigma.z')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s55$))
  else
    valid_value = .false.
  endif
  
case ('sigma.p_z')  
  if (data_source == "beam") then
    datum_value = SQRT(tao_lat%bunch_params(ix1)%sigma(s66$))
  else
    valid_value = .false.
  endif
  
case ('sigma.xy')  
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%sigma(s13$)
  else
    valid_value = .false.
  endif
  
case ('wire.')  
  if (data_source == "beam") then
    read (data_type(6:), '(a)') angle
    datum_value = tao_do_wire_scan (lat%ele(ix1), angle, u%current_beam)
  else
    valid_value = .false.
  endif
  
case ('spin.theta')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%theta
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%theta
  endif
  
case ('spin.phi')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%phi
  else
    call spinor_to_polar (tao_lat%orb(ix1), polar)
    datum_value = polar%phi
  endif
  
case ('spin.polarity')
  if (data_source == "beam") then
    datum_value = tao_lat%bunch_params(ix1)%spin%polarization
  else
    datum_value = 1.0
  endif
  
case ('momentum_compaction')
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

case default
  call out_io (s_error$, r_name, 'UNKNOWN DATUM TYPE: ' // datum%data_type)
  call err_exit

end select

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine load_it (vec, ix0, ix1, datum_value, valid_value, datum, lat, f, coupling_here)

implicit none

type (tao_data_struct) datum
type (lat_struct) lat

real(rp) vec(0:)
real(rp), optional :: f(0:)
real(rp) datum_value

character(20) :: r_name = 'tao_evaluate_a_datum'

integer ix_m, i, ix0, ix1, n_lat, ix_m2
logical, optional :: coupling_here, valid_value

!

if (ix1 < ix0 .and. lat%param%lattice_type == linear_lattice$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // lat%ele(ix0)%name, &
                'IS AFTER ENDING ELEMENT: ' // lat%ele(ix1)%name)
  valid_value = .false.
  return
endif
 
!

n_lat = lat%n_ele_track

if (datum%ele0_name == ' ') then
  ix_m = ix1
  if (present(coupling_here)) call coupling_calc (ix_m, datum, lat)
  datum_value = vec(ix_m)
  if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ix_m))

elseif (datum%merit_type == 'match') then

  ix_m = ix1
  if (present(coupling_here)) call coupling_calc (ix0, datum, lat)
  if (present(coupling_here)) call coupling_calc (ix1, datum, lat)
  datum_value = vec(ix1) - vec(ix0)

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
      ix_m = minloc (vec(0:ix1), 1) - 1
      ix_m2 = minloc (vec(ix0:n_lat), 1) + ix0 - 1
      if (vec(ix_m2) < vec(ix_m2)) ix_m = ix_m2
      datum_value = vec(ix_m)

    case ('max')
      ix_m = maxloc (vec(0:ix1), 1) - 1
      ix_m2 = maxloc (vec(ix0:n_lat), 1) + ix0 - 1
      if (vec(ix_m2) > vec(ix_m2)) ix_m = ix_m2
      datum_value = vec(ix_m)

    case ('abs_min')
      ix_m = minloc (abs(vec(0:ix1)), 1) - 1
      ix_m2 = minloc (abs(vec(ix0:n_lat)), 1) + ix0 - 1
      if (abs(vec(ix_m2)) < abs(vec(ix_m2))) ix_m = ix_m2
      datum_value = abs(vec(ix_m))

    case ('abs_max')
      ix_m = maxloc (abs(vec(0:ix1)), 1) - 1
      ix_m2 = maxloc (abs(vec(ix0:n_lat)), 1) + ix0 - 1
      if (abs(vec(ix_m2)) > abs(vec(ix_m2))) ix_m = ix_m2
      datum_value = abs(vec(ix_m))

    case ('int_min')
      datum_value = 0; ix_m = -1
      call integrate_min (ix1, lat%n_ele_track, datum_value, ix_m, lat, vec, datum)
      call integrate_min (0, ix0, datum_value, ix_m, lat, vec, datum)

    case ('int_max')
      datum_value = 0; ix_m = -1
      call integrate_max (ix1, lat%n_ele_track, datum_value, ix_m, lat, vec, datum)
      call integrate_max (0, ix0, datum_value, ix_m, lat, vec, datum)

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
      call integrate_min (ix0, ix1, datum_value, ix_m, lat, vec, datum)

    case ('int_max')
      datum_value = 0; ix_m = -1
      call integrate_max (ix0, ix1, datum_value, ix_m, lat, vec, datum)

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
if (present(f)) datum%conversion_factor = f(ix_m)

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

subroutine coupling_calc (ixd, datum, lat)

implicit none

type (ele_struct), pointer :: ele
type (this_coupling_struct), pointer :: cc_p
type (tao_data_struct) datum
type (lat_struct) lat

integer ie, ixd
real(rp) f, f1, f2

!

if (cc(ixd)%calc_done) return

cc_p => cc(ixd)
ie = datum%ix_ele
ele => lat%ele(ie)

call c_to_cbar (ele, cc_p%cbar)
f = sqrt(ele%a%beta/ele%b%beta) 
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
! Subroutine tao_transfer_map_calc_at_s (lat, t_map, s1, s2, &
!                                         integrate, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between longitudinal positions
! s1 to s2.
!
! If s2 < s1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will "wrap around" the lattice end.
! For example, if s1 = 900 and s2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
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
!   lat        -- Lat_struct: Lattice used in the calculation.
!   t_map(6)   -- Taylor_struct: Initial map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start position for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end position for the calculation.
!                   Default is lat%param%total_length.
!   integrate  -- Logical, optional: If present and True then do symplectic
!                   integration instead of concatenation. 
!                   Default = False.
!   one_turn   -- Logical, optional: If present and True and s1 = s2 then 
!                   construct the one-turn map from s1 back to s1.
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

  type (lat_struct) lat
  type (taylor_struct) :: t_map(:)

  real(rp), intent(in), optional :: s1, s2
  real(rp) ss1, ss2

  logical, optional :: integrate, one_turn, unit_start
  logical integrate_this, one_turn_this, unit_start_this

  character(40) :: r_name = "tao_transfer_map_calc_at_s"

!

  integrate_this  = logic_option (.false., integrate)
  one_turn_this   = logic_option (.false., one_turn)
  unit_start_this = logic_option(.true., unit_start)

  ss1 = 0;                       if (present(s1)) ss1 = s1
  ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
 
  if (unit_start_this) call taylor_make_unit (t_map)

  if (ss1 == ss2 .and. .not. one_turn_this) return

! Normal case

if (ss1 < ss2) then 
  call transfer_this (ss1, ss2)

! For a circular lattice push through the origin.

elseif (lat%param%lattice_type == circular_lattice$) then
  call transfer_this (ss1, lat%param%total_length)
  call transfer_this (0.0_rp, ss2)

! For a linear lattice compute the backwards matrix

else
  if (.not. unit_start_this) then
    call out_io (s_fatal$, r_name, 'Backwards propagation with a non-unit starting map!')
    call err_exit
  endif

  call transfer_this (ss2, ss1)
  call taylor_inverse (t_map, t_map)
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

  if (ix_ele /= ix_ele_old) ele = lat%ele(ix_ele)
  s_now = s_1
  kill_it = .false.

  do
    s_end = min(s_2, ele%s)
    if (ele%key == sbend$) then
      if (s_now /= lat%ele(ix_ele-1)%s) ele%value(e1$) = 0
      if (s_end /= ele%s) ele%value(e2$) = 0
      if (s_now == lat%ele(ix_ele-1)%s) kill_it = .true.
      if (s_end == ele%s) kill_it = .true.
    endif

    ds = s_end - s_now
    ele%value(l$) = ds

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
    ele = lat%ele(ix_ele)
  enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_mat6_calc_at_s (lat, mat6, vec0, s1, s2, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between longitudinal positions
! s1 to s2.
!
! If s2 < s1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will "wrap around" the lattice end.
! For example, if s1 = 900 and s2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
!
! If s2 = s1 then you get the unit matrix except if one_turn = True.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Lat_struct: Lattice used in the calculation.
!   mat6(6,6)  -- Real(rp): Initial matrix (used when unit_start = False)
!   vec0(6)    -- Real(rp): Initial 0th order map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start index for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end index for the calculation.
!                   Default is lat%param%total_length.
!   one_turn   -- Logical, optional: If present and True then construct the
!                   one-turn map from s1 back to s1 (ignolat s2).
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    mat6(6,6) -- Real(rp): Transfer matrix.
!    vec0(6)   -- Real(rp): 0th order part of the map.
!-

subroutine tao_mat6_calc_at_s (lat, mat6, vec0, s1, s2, one_turn, unit_start)

use bmad_struct
use bmad_interface

implicit none

type (lat_struct) lat

real(rp) mat6(:,:), vec0(:)
real(rp), intent(in), optional :: s1, s2
real(rp) ss1, ss2

logical, optional :: one_turn, unit_start
logical one_turn_this

!

one_turn_this = logic_option (.false., one_turn)

ss1 = 0;                       if (present(s1)) ss1 = s1
ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
 
if (logic_option(.true., unit_start)) then
  call mat_make_unit (mat6)
  vec0 = 0
endif

! Normal case

if (ss1 < ss2 .or. (ss1 == ss2 .and. one_turn_this)) then
  call transfer_this (ss1, ss2)

! For a circular lattice push through the origin.

elseif (lat%param%lattice_type == circular_lattice$) then
  call transfer_this (ss1, lat%param%total_length)
  call transfer_this (0.0_rp, ss2)

! For a linear lattice compute the backwards matrix

else
  call transfer_this (ss2, ss1)
  call mat_inverse (mat6, mat6)
  vec0 = -matmul(mat6, vec0)

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
ele = lat%ele(ix_ele)
s_now = s_1

do
  s_end = min(s_2, ele%s)
  ds = s_end - s_now
  ele%value(l$) = ds
  if (ele%key == sbend$) then
    if (s_now /= lat%ele(ix_ele-1)%s) ele%value(e1$) = 0
    if (s_end /= ele%s) ele%value(e2$) = 0
  endif

  if (s%global%matrix_recalc_on) call make_mat6 (ele, lat%param)

  mat6 = matmul (ele%mat6, mat6)
  vec0 = matmul (ele%mat6, vec0) + ele%vec0

  if (s_end == s_2) return
  s_now = s_end
  ix_ele = ix_ele + 1
  ele = lat%ele(ix_ele)
enddo

end subroutine

end subroutine

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
! Function tao_valid_datum_index (lat, ix_ele, valid) result (ix_loc)
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
!
!   valid  -- Logical: Set false if element does not have a definite location.
!   ix_loc -- Integer: Location of element in the tracking part of the lat%ele(:) array.
!-

function tao_valid_datum_index (lat, ix_ele, valid) result (ix_loc)

implicit none

type (lat_struct) lat
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

if (ix_ele <= lat%n_ele_track) return

if (lat%ele(ix_ele)%control_type == super_lord$) then
  ixc = lat%ele(ix_ele)%ix2_slave
  ix_loc = lat%control(ixc)%ix_slave
  return
endif

valid = .false.
call out_io (s_error$, r_name, &
            'ELEMENT: ' // lat%ele(ix_ele)%name, &
            'WHICH IS A: ' // control_name(lat%ele(ix_ele)%control_type), &
            'CANNOT BE USED IN DEFINING A DATUM SINCE IT DOES NOT HAVE ', &
            '   A DEFINITE LOCATION IN THE LATTICE.')

end function

end module tao_data_mod
