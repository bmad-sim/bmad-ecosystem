!+
! Subroutine CHOOSE_CESR_LATTICE (LATTICE, LAT_FILE, CURRENT_LAT, RING, CHOICE)
!
! Subroutine to let the user choose a lattice. The subroutine will present a
! list to choose from.
!                                                               
! Modules Needed:
!   use bmad
!
! Input:
!   CURRENT_LAT -- Character*40: Name of current lattice (will be stared in
!                       the list presented to the user).
!                       Use CALL GETLAT (CURRENT_LAT) to get the current name.
!                       Set CURRENT_LAT = ' ' if you do not want to use this
!                       feature.
!                       NOTE: You must be connected to the mpm to use GETLAT.
!   CHOICE      -- Character*(*): [Optional] If present then this will be
!                       used as input instead of querying the user.
!
! Output:
!   LATTICE  -- Character*40: Lattice name choisen. If a file name is given
!                    and RING is not present then LATTICE = ""
!   LAT_FILE -- Character*(*): Name of the lattice file. Typically:
!                    lat_file = 'U:[CESR.BMAD.LAT]BMAD_' // lattice // .LAT
!   RING     -- Ring_struct: OPTIONAL. If present then BMAD_PARSER is called
!               to load the RING structure.
!-

!$Id$
!$Log$
!Revision 1.9  2003/01/27 14:40:31  dcs
!bmad_version = 56
!
!Revision 1.8  2002/02/23 20:32:12  dcs
!Double/Single Real toggle added
!
!Revision 1.7  2002/01/11 17:06:59  dcs
!Fix Bug
!
!Revision 1.6  2002/01/08 21:44:38  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.5  2001/10/08 17:18:14  rwh24
!DCS changes to f90 files.
!Bug fixes to c file.
!
!Revision 1.4  2001/10/05 18:23:57  rwh24
!Bug Fixes
!
!Revision 1.3  2001/10/02 18:49:11  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

               
subroutine choose_cesr_lattice (lattice, lat_file, current_lat, ring, choice)

  use bmad_struct
  use bmad_interface
  use cesr_utils

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
      type *
      i_lat = 0
      do i = 1, num_lats
        if (lat_list(i) == current_lat) then
          type '(1x, a, i2, 3a)', '**', i, ') ', trim(lat_list(i)), &
                                           '   ! Current lattice in Data Base'
          i_lat = i
        else
          type '(i5, 2a)', i, ') ', lat_list(i)
        endif
      enddo
  
      type *, ' [Note: To be in this list a lattice file must have a name   ]'
      type *, ' [      of the form: U:[CESR.BMAD.LAT]bmad_<lattice_name>.lat]'

      type *
      type *, 'You can enter a Lattice number or a full file name.'
      if (i_lat == 0) then
        type '(a, $)', ' Choice: '
      else
        type '(a, i3, a, $)', ' Choice: <CR =', i_lat, '> '
      endif
      accept '(a)', line
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
        type *, 'ERROR: WHICH LATTICE? TRY AGAIN...'
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
        lat_file = 'BMAD_LAT:bmad_' // line 
        if (index(line, '.') == 0) lat_file = trim(lat_file) // '.lat' 
        call FullFileName(lat_file, lat_file)
        inquire (file = lat_file, exist = is_there, name = lat_file)
        if (.not. is_there) then
          type *, 'READ ERROR OR FILE DOES NOT EXIST. TRY AGAIN...'
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
    lattice = ring%lattice
  endif

end subroutine
