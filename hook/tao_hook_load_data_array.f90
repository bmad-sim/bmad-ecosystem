!+
! Subroutine tao_hook_load_data_array (found, datum, lattice, orb, datum_value)
!
! Custom data types are defined here.
!-

subroutine tao_hook_load_data_array (found, datum, lattice, orb, datum_value)

use tao_mod

implicit none

type (tao_data_struct) datum
type (ring_struct) lattice
type (coord_struct) orb(0:)

real(rp) datum_value
logical found

integer ix1, ix2
type (ele_struct), pointer :: ele

character(20) :: r_name = 'tao_hook_load_data_array'
!



found = .true.

select case (datum%data_type)

  case default
    found = .false.

end select


contains
!-----------------------------------------------------------------------
!This internal subroutine is for use if there is a range of
!elemements associated with a datum. If there is only one then load_it
!just loads the data for element ix1, but if there are two lattice elements 
!associated with a datum load_it will find the proper model value
!for that datum based on the merit type. This subroutine doesn't have to be used
!and is included just out of convenience.

subroutine load_it (vec)

real(rp) vec(0:)
integer ix_m, i

!

if (ix2 < 0) then
  ix_m = ix1

!

else

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

end subroutine load_it

end subroutine
