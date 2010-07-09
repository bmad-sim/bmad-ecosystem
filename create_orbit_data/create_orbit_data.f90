program create_orbit_data
  !+
  ! create_orbit_data
  !   create simulated orbit data by tracking through a lattice
  !
  ! mjf 2010.07.08
  !-

  use bmad
  use cesr_basic_mod
  use precision_def
  use nr

  implicit none

  include 'det_list.inc'

  integer i,i1,iturn,k,ibpm, af
  integer pntr
  integer bpm_num, bpm_ptr

  character*2 ew
  character*4 bpm_char
  character*12 label
  character(8) :: labels(140)
  character*3 new_label
  character*80 file
  real(8) gerr, noise
  real(8) :: xamp, yamp
  integer :: numturns
  character*10 :: formt
  character*20 turnschar
  character*80 :: prefix

  logical :: found, err

  integer :: bpm_to_lat(bpm_tot), n_loc
  character(14) :: outfilename
  character*40 :: lattice, cur_lattice
  character*80 :: lat_file
  type (lat_struct) :: lat
  type (coord_struct), allocatable :: orb(:), coord_by_turn(:,:)  ! Turn, ele_ix
  type (ele_pointer_struct), allocatable :: eles(:)

!  call mpm_init('BMAD')
!  call getlat( cur_lattice )

  call choose_cesr_lattice( lattice, lat_file, cur_lattice, lat )

!  call mpm_goodbye


  print *, "Number of turns (multiple of 2): "
  read *, numturns

  allocate(coord_by_turn(1:numturns,0:lat%n_ele_max))
  allocate(orb(0:lat%n_ele_max))

  write (turnschar,*) numturns
  turnschar = adjustl(turnschar)

  call init_random_seed()



  print*, "Enter level of noise to add to data in m"
  read *, noise

  print *, "X shaking amplitude in m: "
  read *, xamp
  print *, "Y shaking amplitude in m: "
  read *, yamp

  Print *, "Prefix for filename: "
  read *, prefix

  ! create a map of desired bpms' element index numbers
  call lat_ele_locator ("marker::det*", lat, eles, n_loc, err)

  do ibpm = 1,bpm_tot
    do i = 1, n_loc ! Loop over all elements found to match the search string.
      found = .false.
      if (eles(i)%ele%name == use_bpm(ibpm)) then
        bpm_to_lat(ibpm) = eles(i)%ele%ix_ele
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      print *, "ERROR, couldn't find ", use_bpm(ibpm)
    end if


  enddo
  deallocate(eles)


  do i1=1,2	! horz then vert
     if (i1 == 1) then
        write(file,'(a,a6)') trim(adjustl(prefix)), "-x.dat" 
     else
        write(file,'(a,a6)') trim(adjustl(prefix)), "-y.dat" 
     end if
     
     open(unit=21,name=file)
     write(21,8), bpm_tot
8    format(1x/28x,i3)
     
     ! Set starting amplitude
     if(i1==1) then
        coord_by_turn(1,0)%vec(1:6) = 0
        coord_by_turn(1,0)%vec(1) = xamp
     else
        coord_by_turn(1,0)%vec(:) = 0
        coord_by_turn(1,0)%vec(3) = yamp
     endif
     orb(0:) = coord_by_turn(1,0:)

     ! track through turns
     do iturn = 1, numturns
        call track_all(lat, orb)
        coord_by_turn(iturn,0:) = orb(0:)
        orb(0) = orb(lat%n_ele_track)
     enddo
!    call track(t,x, numturns)


    call init_random_seed()

    print *, "Number of turns: ", numturns
    do ibpm = 1, bpm_tot
      bpm_ptr = bpm_to_lat(ibpm)
      write (21, *) "x"

      write (21,12) "Location:", &
           use_bpm(ibpm)
      do af=1, 32
        write (21, *) "x"
      end do
      write (21, *) "# --- Bunch 1 ---"
12    format (2x,a9, 6x, a)
      do iturn = 1, numturns  ! turn number

        do k = 1, 2
          call gasdev_s(gerr)
          coord_by_turn(iturn,bpm_ptr)%vec(k*2-1) = &
               coord_by_turn(iturn,bpm_ptr)%vec(k*2-1) + &
               (gerr * noise).   ! m
        enddo

        write(21,7) (coord_by_turn(iturn,bpm_ptr)%vec(i),i=1,3,2) ! x, y position
7       format(55x,f9.6,2x,f9.6)

      enddo


    enddo

    close(unit=21) 

  enddo

  deallocate(coord_by_turn)
  deallocate(orb)
  stop
end program create_orbit_data



subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine init_random_seed

