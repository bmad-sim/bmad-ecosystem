program test

  use bmad                 ! Define the structures we need to know about.
  implicit none
  type (ring_struct) ring  ! This structure holds the lattice info
  integer i, ix

! Read in a lattice and calculate the twiss parameters.

  call bmad_parser ("bmad_test.lat", ring)    ! Read in a lattice.
  call twiss_at_start (ring)           ! Calculate starting Twiss params.
  call twiss_propagate_all (ring)      ! Propagate Twiss parameters

! Print info on the first 11 elements

  print *, ' Ix  Name              Ele_type                   S      Beta_x'
  do i = 0, 10
    print '(i4, 2x, a, 2x, a, 2f12.4)', i, ring%ele_(i)%name, &
                    key_name(ring%ele_(i)%key), &
                    ring%ele_(i)%s, ring%ele_(i)%x%beta
  enddo

! Find the CLEO_SOL element and print information on it.

  call element_locator ('CLEO_SOL', ring, ix)
  print *
  print *, '!---------------------------------------------------------'
  print *, '! Information on element: CLEO_SOL', ix
  print *
  call type_ele (ring%ele_(ix), .false., 0, .false., 0, .true., ring)

end program
