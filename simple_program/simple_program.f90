program test

  use bmad                 ! Define the structures we need to know about.
  implicit none
  type (lat_struct) ring  ! This structure holds the lattice info
  integer i, ix

! Read in a lattice and calculate the twiss parameters.

  call bmad_parser ("simple_program/lat.bmad", ring)    ! Read in a lattice.
  call twiss_at_start (ring)           ! Calculate starting Twiss params.
  call twiss_propagate_all (ring)      ! Propagate Twiss parameters

! Print info on the first 11 elements

  print *, ' Ix  Name              Ele_type                   S      Beta_a'
  do i = 0, 10
    print '(i4, 2x, a16, 2x, a, 2f12.4)', i, ring%ele(i)%name, &
                    key_name(ring%ele(i)%key), &
                    ring%ele(i)%s, ring%ele(i)%a%beta
  enddo

! Find the CLEO_SOL element and print information on it.

  call element_locator ('CLEO_SOL', ring, ix)
  print *
  print *, '!---------------------------------------------------------'
  print *, '! Information on element: CLEO_SOL', ix
  print *
  call type_ele (ring%ele(ix), .false., 0, .false., 0, .true., ring)

end program
