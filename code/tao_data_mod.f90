module tao_data_mod

use tao_mod
use macroparticle_mod
use macro_utils_mod

! These are data types specific to macroparticles

type this_coupling_struct
  real(rp) cbar(2,2)
  real(rp) coupling11, coupling12a, coupling12b, coupling22
  real(rp) f_11, f_12a, f_12b, f_22
  logical calc_done
end type

type (this_coupling_struct), save, allocatable, target :: cc(:)

contains

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
!   ix_ele  -- Integer: element to valuate data at
!-

subroutine tao_load_data_array (u, ix_ele)

implicit none

type (tao_universe_struct) :: u
type (tao_ix_data_struct), pointer :: ix_data
integer uni, ix_ele
integer i
character(20) :: r_name = 'tao_load_data_array'

  if (ix_ele .eq. 0) call tao_data_coupling_init (u) 
  
  ! find which datums to evaluate here
  if (.not. associated(u%ix_data(ix_ele)%ix_datum)) return

  ix_data => u%ix_data(ix_ele)
  do i = 1, size(ix_data%ix_datum)
    call tao_evaluate_a_datum (u%data(ix_data%ix_datum(i)), u, u%model, &
              u%model_orb, ix_ele, u%data(ix_data%ix_datum(i))%model_value)
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
! Subroutine tao_evaluate_a_datum (datum, u, lattice, orb, ix_ele, datum_value)
!
! Subroutine to put the proper data in the specified datum
! Data types that require ranges should be evaulated at the ix2 element.
!-

subroutine tao_evaluate_a_datum (datum, u, lattice, orb, ix_ele, datum_value)

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele
type (tao_data_struct) datum
type (ring_struct) ::  lattice
type (coord_struct) :: orb(0:)
type (modes_struct) mode
type (taylor_struct), save :: taylor(6)
!type (macro_bunch_params_struct), save :: macro_params

real(rp) datum_value, mat6(6,6)

integer ix_ele
integer, save :: ix_save = -1
integer i, j, k, m, ix, ix1, ix2, expnt(6)
!integer track_type

character(20) :: r_name = 'tao_evaluate_a_datum'
character(16) data_type

logical found

!

call tao_hook_load_data_array (found, datum, u, lattice, orb, ix_ele, datum_value)
if (found) return

ix1 = datum%ix_ele
ix2 = datum%ix_ele2
ele => lattice%ele_(ix_ele)

data_type = datum%data_type
if (data_type(1:2) == 'r:') data_type = 'r:'
if (data_type(1:2) == 't:') data_type = 't:'
if (data_type(1:2) == 'tt:') data_type = 'tt:'


select case (data_type)

case ('orbit:x')
  call load_it (orb(:)%vec(1))
case ('orbit:y')
  call load_it (orb(:)%vec(3))
case ('orbit:z')
  call load_it (orb(:)%vec(5))

case ('bpm:x')
  call tao_read_bpm (orb(ix_ele), lattice%ele_(ix_ele), x$, datum_value)
case ('bpm:y')
  call tao_read_bpm (orb(ix_ele), lattice%ele_(ix_ele), y$, datum_value)

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
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%beta
  endif
    
case ('beta:y')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%beta)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%beta
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%beta
  endif

case ('alpha:x')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%x%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%x%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%x%alpha
  endif
  
case ('alpha:y')
  if (s%global%track_type .eq. "single") then
    call load_it (lattice%ele_(:)%y%alpha)
  elseif (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%y%alpha
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%y%alpha
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

case ('tt:')
  call relative_switch
  expnt = 0
  i = read_this_index (datum%data_type(4:4)); if (i == 0) return
  do j = 5, 15
    if (datum%data_type(j:j) == ' ') exit
    k = read_this_index (datum%data_type(j:j)); if (k == 0) return
    expnt(k) = expnt(k) + 1
  enddo
  call transfer_map_calc (lattice, taylor, ix1, ix2)
  datum_value = taylor_coef (taylor(i), expnt)

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
  
case ('norm_emittance:a')
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%a%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%a%norm_emitt
  else
    datum_value = 0.0
  endif
  
case ('norm_emittance:b')  
  if (s%global%track_type .eq. "beam") then
    datum_value = u%beam%params%b%norm_emitt
  elseif (s%global%track_type .eq. "macro") then
    datum_value = u%macro_beam%params%b%norm_emitt
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
    call out_io (s_abort$, r_name, 'BAD INDEX CONSTRAINT: ' // datum%data_type)
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

end subroutine coupling_calc

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_read_bpm
!
! Find the reading on a bpm given the orbit at the BPM and the BPM offsets
!
! Input: 
!  orb      -- Coord_struct: Orbit position at BPM
!  ele      -- Ele_struct: the BPM
!  axis     -- Integer: x$ or y$
!
! Output:
!  reading  -- Real(rp): BPM reading
!
!-

subroutine tao_read_bpm (orb, ele, axis, reading)

type (coord_struct) orb
type (ele_struct) ele
type (tao_d2_data_struct), pointer :: d2_ptr

real(rp) reading

integer axis

character(20) :: r_name = "tao_read_bpm"

logical err

  if (ele%key .ne. monitor$ .and. ele%key .ne. instrument$) then
    reading = 0.0
    return
  endif

  if (axis .eq. x$) then
    reading =  (orb%vec(1) - ele%value(x_offset_tot$)) * cos(ele%value(tilt_tot$)) + &
               (orb%vec(3) - ele%value(y_offset_tot$)) * sin(ele%value(tilt_tot$))              
  elseif (axis .eq. y$) then
    reading = -(orb%vec(1) - ele%value(x_offset_tot$)) * sin(ele%value(tilt_tot$)) + &
               (orb%vec(3) - ele%value(y_offset_tot$)) * cos(ele%value(tilt_tot$))  
  else
    call out_io (s_warn$, r_name, &
                 "This axis not supported for BPM reading!")
  endif
    
  call tao_find_data (err, s%u(3), 'bpm:x', d2_ptr)

end subroutine tao_read_bpm


end module tao_data_mod
