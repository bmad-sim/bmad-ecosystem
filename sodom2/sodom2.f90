program sodom2

use sodom2_mod
!use bmad
!!use sim_utils

implicit none

!type (lat_struct), target :: lat   ! This structure holds the lattice info
type (sodom2_params_struct) sodom
type (sodom2_com_struct), target :: sodom2_com

! Programs should always implement "intelligent bookkeeping".
bmad_com%auto_bookkeeper = .false.

!call bmad_parser ("lat.bmad", lat)  ! Read in a lattice.
call sodom2_read_params(sodom, sodom2_com)
call sodom2_init_params(sodom, sodom2_com)
print *, 'Initializing bunch...'
!call run_timer ('START')
call sodom2_init_bunch(sodom, sodom2_com)
print *, 'Tracking 1-turn...'
call sodom2_track_bunch(sodom, sodom2_com)
print *,  'Tracking complete. Calculating ADST and n-axis...'
call sodom2_construct_quaternions(sodom, sodom2_com)

call sodom2_construct_mat(sodom, sodom2_com)

call sodom2_eig(sodom2_com)

call sodom2_calc_ADST(sodom, sodom2_com)

call sodom2_calc_n(sodom, sodom2_com)

print *,  'Writing n-axis output to file...'
call sodom2_write(sodom, sodom2_com)

call sodom2_deallocate_memory(sodom2_com)

print *,  'Complete.'

!call lat_ele_locator ('CLEO_SOL', lat, eles, n_loc, err)  ! Find element
!cleo => eles(1)%ele                        ! Point to cleo_sol element.
!call pointer_to_attribute (cleo, 'KS', .true., a_ptr, err) ! Point to KS attribute.
!a_ptr%r = a_ptr%r + 0.001         ! Modify KS component.
!call set_flags_for_changed_attribute (cleo, cleo%value(ks$))
!call lattice_bookkeeper (lat)
!call lat_make_mat6 (lat, cleo%ix_ele)      ! Remake transfer matrix

! Calculate starting Twiss params if the lattice is closed, 
! and then propagate the Twiss parameters through the lattice.

!if (lat%param%geometry == closed$) call twiss_at_start (lat)
!call twiss_propagate_all (lat)      ! Propagate Twiss parameters

! Print info on the first 11 elements

!print *, ' Ix  Name              Ele_type                   S      Beta_a'
!do i = 0, 10
!  ele => lat%ele(i)
!  print '(i4,2x,a16,2x,a,2f12.4)', i, ele%name, key_name(ele%key), ele%s, ele%a%beta
!enddo

! print information on the CLEO_SOL element.

!print *
!print *, '!---------------------------------------------------------'
!print *, '! Information on element: CLEO_SOL'
!print *
!call type_ele (cleo, .false., 0, .false., 0, .true.)

!deallocate (eles)

end program
