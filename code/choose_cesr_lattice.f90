!+
! Subroutine choose_cesr_lattice (lattice, lat_file, current_lat, ring, choice)
!
! Subroutine to let the user choose a lattice. The subroutine will present a
! list to choose from.
!                                                               
! Modules Needed:
!   use bmad
!
! Input:
!   current_lat -- Character*40: Name of current lattice (will be stared in
!                       the list presented to the user).
!                       Use CALL GETLAT (CURRENT_LAT) to get the current name.
!                       Set CURRENT_LAT = ' ' if you do not want to use this
!                       feature.
!                       NOTE: You must be connected to the mpm to use GETLAT.
!   choice      -- Character*(*): [Optional] If present then this will be
!                       used as input instead of querying the user.
!
! Output:
!   lattice  -- Character*40: Lattice name choisen. If a file name is given
!                    and RING is not present then LATTICE = ""
!   lat_file -- Character*(*): Name of the lattice file. Typically:
!                    lat_file = 'U:[CESR.BMAD.LAT]BMAD_' // lattice // .LAT
!   ring     -- Ring_struct: OPTIONAL. If present then BMAD_PARSER is called
!               to load the RING structure.
!-

#include "CESR_platform.inc"
               
subroutine choose_cesr_lattice (lattice, lat_file, current_lat, ring, choice)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), optional :: ring

  character(len=*), optional :: choice
  character*(*) lat_file
  character*40 lattice, current_lat, lat_list(100)
  character*80 line
   
  integer i, num_lats, i_lat, ix, ios

  logical is_there, ask_for_lat, default_flag

!                   

  call get_lattice_list (lat_list, num_lats, 'BMAD_LAT:')

  ask_for_lat = .true.

  if (present(choice)) then
    line = choice
    call string_trim (line, line, ix)
    if (ix /= 0) ask_for_lat = .false.
  endif

! loop until we have a valid choice

  do

    if (ask_for_lat) then
      print *
      i_lat = 0
      do i = 1, num_lats
        if (lat_list(i) == current_lat) then
          print '(1x, a, i2, 3a)', '**', i, ') ', trim(lat_list(i)), &
                                           '   ! Current lattice in Data Base'
          i_lat = i
        else
          print '(i5, 2a)', i, ') ', lat_list(i)
        endif
      enddo
  
      print *, ' [Note: To be in this list a lattice file must have a name   ]'
      print *, ' [      of the form: U:[CESR.BMAD.LAT]bmad_<lattice_name>.lat]'

      print *
      print *, 'You can enter a Lattice number or a full file name.'
      if (i_lat == 0) then
        print '(a, $)', ' Choice: '
      else
        print '(a, i3, a, $)', ' Choice: <CR =', i_lat, '> '
      endif
      read (*, '(a)') line
    endif

    call string_trim (line, line, ix)
    line = line(:ix)

    if (ix == 0 .or. (ix == 1 .and. line == '*')) then
      default_flag = .true.
      do i_lat = 1, num_lats
        if (lat_list(i_lat) == current_lat) exit
      enddo
    else
      default_flag = .false.
      read (line, *, iostat = ios) i_lat
    endif

    if (default_flag .or. (ios == 0 .and. index('0123456789', line(1:1)) /= 0)) then
      if (i_lat < 1 .or. i_lat > num_lats) then
        print *, 'ERROR: WHICH LATTICE? TRY AGAIN...'
        ask_for_lat = .true.
        cycle  ! try again
      endif
      lattice = lat_list(i_lat)
      call lattice_to_bmad_file_name (lattice, lat_file)
    else
      lattice = ""
      lat_file = line
      inquire (file = lat_file, exist = is_there, name = lat_file)
      if (.not. is_there) then
        lattice = line
        lat_file = 'BMAD_LAT:bmad_' // lattice
        if (index(lattice, '.') == 0) lat_file = trim(lat_file) // '.lat' 
        ix = index(lattice, '.lat')
        if (ix /= 0) lattice = lattice(:ix-1)
        call FullFileName(lat_file, lat_file)
        inquire (file = lat_file, exist = is_there, name = lat_file)
        if (.not. is_there) then
          print *, 'READ ERROR OR FILE DOES NOT EXIST. TRY AGAIN...'
          ask_for_lat = .true.
          cycle
        endif
      endif
      ix = index(lat_file, ';')
      if (ix /= 0) lat_file = lat_file(:ix-1)
    endif
    exit

  enddo

! load ring if present

  if (present (ring)) then
    call bmad_parser (lat_file, ring)
    if (lattice /= "") then
      if (lattice /= ring%lattice) print *, &
           'WARNING FROM CHOOSE_CESR_LATTICE: LATTICE NAME IN RING DOES MATCH FILE NAME!'
    endif
    lattice = ring%lattice
  endif

end subroutine
