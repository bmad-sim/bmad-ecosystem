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

use tao_mod

implicit none

type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: data

integer i, j, k, m, ix, ix1, ix2
character(20) :: r_name = 'tao_load_data_array'
logical found

type this_coupling_struct
  real(rp) cbar(2,2)
  real(rp) coupling11, coupling12a, coupling12b, coupling22
  real(rp) f_11, f_12a, f_12b, f_22
  logical calc_done
end type
type (this_coupling_struct), save, allocatable :: cc(:)

! init setup for coupling calculation

  m = 0
  do i = 1, size(s%u)
    m = max(m, size(s%u(i)%model%ele_))
  enddo

  if (.not. allocated(cc)) allocate (cc(0:m))

  if (size(cc) < m+1) then
    deallocate(cc)
    allocate(cc(0:m))
  endif

! loop over all universes 

do i = 1, size(s%u)

  u => s%u(i)

! calc coupling if needed

  cc%calc_done = .false.
  do j = 1, size(u%data)
    if (u%data(j)%type(1:9) /= 'coupling:' .and. &
              u%data(j)%type(1:5) /= 'cbar:') cycle 
    ix1 = u%data(j)%ix_ele
    ix2 = u%data(j)%ix_ele2
    if (ix1 < 0) cycle
    if (ix2 < 0) then
      call coupling_calc (u%model%ele_(u%data(j)%ix_ele), cc(k))
    else
      do k = ix1, ix2
        call coupling_calc (u%model%ele_(k), cc(k))
      enddo
    endif
  enddo

! loop over all data structs

  do j = 1, size(u%data)
  
    data => u%data(j)
    if (.not. data%exists) cycle
    call tao_hook_load_data_array (data, found)
    if (found) cycle

    ix1 = data%ix_ele
    ix2 = data%ix_ele2
    ele => u%model%ele_(ix1)

! loop over all data 

    select case (data%type)

    case ('orbit:x')
      call load_it (u%model_orb(:)%vec(1))
    case ('orbit:y')
      call load_it (u%model_orb(:)%vec(3))
    case ('orbit:z')
      call load_it (u%model_orb(:)%vec(5))

    case ('phase:x')
      call load_it (u%model%ele_(:)%x%phi)
    case ('phase:y')
      call load_it (u%model%ele_(:)%y%phi)

    case ('beta:x')
      call load_it (u%model%ele_(:)%x%beta)
    case ('beta:y')
      call load_it (u%model%ele_(:)%y%beta)

    case ('eta:x')
      call load_it (u%model%ele_(:)%x%eta)
    case ('eta:y')
      call load_it (u%model%ele_(:)%y%eta)

    case ('coupling:11')
      call load_it (cc%coupling11, cc%f_11)
    case ('coupling:12a')
      call load_it (cc%coupling12a, cc%f_12a)
    case ('coupling:12b')
      call load_it (cc%coupling12b, cc%f_12b)
    case ('coupling:22')
      call load_it (cc%coupling22, cc%f_22)

    case ('cbar:11')
      call load_it (cc%cbar(1,1))
    case ('cbar:12')
      call load_it (cc%cbar(1,2))
    case ('cbar:21')
      call load_it (cc%cbar(1,2))
    case ('cbar:22')
      call load_it (cc%cbar(2,2))

    case default
      call out_io (s_error$, r_name, 'UNKNOWN DATA TYPE: ' // data%type)
      return

    end select

  enddo

enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine load_it (vec, f)

real(rp) vec(0:)
real(rp), optional :: f(0:)
integer ix_m

!

if (ix2 < 0) then
  ix_m = ix1
else
  select case (data%merit_type)
  case ('min')
    ix_m = minloc (vec(ix1:ix2), 1)
  case ('max')
    ix_m = maxloc (vec(ix1:ix2), 1)
  case ('abs_min')
    ix_m = minloc (abs(vec(ix1:ix2)), 1)
  case ('abs_max')
    ix_m = maxloc (abs(vec(ix1:ix2)), 1)
  case default
    call out_io (s_abort$, r_name, 'BAD MERIT_TYPE: ' // data%merit_type, &
                                   'FOR DATA: ' // data%type)
    call err_exit
  end select
endif

data%model_value = vec(ix_m)
if (data%merit_type(1:4) == 'abs_') data%model_value = abs(vec(ix_m))
if (present(f)) data%conversion_factor = f(ix_m)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! contains

subroutine coupling_calc (ele, cc)

type (ele_struct) ele
type (this_coupling_struct) cc

real(rp) f, f1, f2

!

if (cc%calc_done) return

call c_to_cbar (ele, cc%cbar)
f = sqrt(ele%x%beta/ele%y%beta) 
f1 = f / ele%gamma_c
f2 = 1 / (f * ele%gamma_c)

cc%coupling11  = cc%cbar(1,1) * f1
cc%coupling12a = cc%cbar(1,2) * f2
cc%coupling12b = cc%cbar(1,2) * f1
cc%coupling22  = cc%cbar(2,2) * f2

cc%f_11  = f1
cc%f_12a = f2
cc%f_12b = f1
cc%f_22  = f2

end subroutine

end subroutine
