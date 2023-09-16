program test

use bmad                 ! Define the structures we need to know about.
implicit none
type (lat_struct), target :: lat   ! This structure holds the lattice info
type (ele_struct), pointer :: ele, cleo
type (ele_pointer_struct), allocatable :: eles(:)
type (all_pointer_struct) a_ptr
integer i, ix, n_loc
logical err

! Programs should always implement "intelligent bookkeeping".
bmad_com%auto_bookkeeper = .false.

! Read in a lattice, and modify the ks solenoid strength of "cleo_sol".

call bmad_parser ("lat.bmad", lat)  ! Read in a lattice.

call lat_ele_locator ('CLEO_SOL', lat, eles, n_loc, err)  ! Find element
cleo => eles(1)%ele                        ! Point to cleo_sol element.
call pointer_to_attribute (cleo, 'KS', .true., a_ptr, err) ! Point to KS attribute.
a_ptr%r = a_ptr%r + 0.001         ! Modify KS component.
call set_flags_for_changed_attribute (cleo, cleo%value(ks$))
call lattice_bookkeeper (lat)
call lat_make_mat6 (lat, cleo%ix_ele)      ! Remake transfer matrix

! Calculate starting Twiss params if the lattice is closed, 
! and then propagate the Twiss parameters through the lattice.

if (lat%param%geometry == closed$) call twiss_at_start (lat)
call twiss_propagate_all (lat)      ! Propagate Twiss parameters

! Print info on the first 11 elements

print *, ' Ix  Name              Ele_type                   S      Beta_a'
do i = 0, 10
  ele => lat%ele(i)
  print '(i4,2x,a16,2x,a,2f12.4)', i, ele%name, key_name(ele%key), ele%s, ele%a%beta
enddo

! print information on the CLEO_SOL element.

print *
print *, '!---------------------------------------------------------'
print *, '! Information on element: CLEO_SOL'
print *
call type_ele (cleo, .false., 0, .false., 0, .true.)

deallocate (eles)

end program
