module tao_data_mod

use tao_mod

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

case ('s_length')
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

end module
