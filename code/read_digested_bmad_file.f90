!+
! Subroutine READ_DIGESTED_BMAD_FILE (IN_FILE_NAME, RING, VERSION)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! WRITE_DIGESTED_BMAD_FILE.F77]
!
! Note: This subroutine also reads in the common structures for BMAD_PARSER2
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     IN_FILE_NAME -- Character*(*): Name of the digested file
!
! Output:
!     RING      -- Ring_struct: Output structure
!     VERSION   -- Integer: Version number of RING.
!     STATUS    -- Common block status structure
!       .OK       -- Set .false. if read failure.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine read_digested_bmad_file (in_file_name, ring, version)

  use bmad_struct

  implicit none

  type (ring_struct)   ring

  integer d_unit, lunget, n_files, version, i, ix

  character*(*) in_file_name
  character*72 fname(3)

  logical found_it

!

  d_unit = lunget()
  bmad_status%ok = .true.

  open (unit = d_unit, file = in_file_name, status = 'old',  &
                     form = 'unformatted', readonly, shared, err = 9000)

  read (d_unit, err = 9100) n_files, version

  if (version < bmad_inc_version$) then
    if (bmad_status%type_out) type '(1x, a, i4, a, i4)',  &
            'READ_DIGESTED_BMAD_FILE: DIGESTED FILE VERSION OUT OF DATE',  &
            version, ' <', bmad_inc_version$
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

  if (version > bmad_inc_version$) then
    if (bmad_status%type_out) then
      type *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE HAS VERSION:',  &
                                                              version
      type *, '     GREATER THAN VERSION OF THIS PROGRAM:',  &
                                                  bmad_inc_version$
      type *, '     WILL NOT USE THE DIGESTED FILE.'
      type *, '     YOU SHOULD RECOMPILE THIS PROGRAM.'
    endif
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

  do i = 1, n_files
    read (d_unit, err = 9100) fname(1)
    ix = index(fname(1), ';')
    if (ix <= 1) then
      type *, 'READ_DIGESTED_BMAD_FILE: FILE LIST READ ERROR! ', fname(1)
      close (d_unit)
      bmad_status%ok = .false.
      return
    else
      fname(2) = fname(1)(:ix-1)
    endif
    inquire (file = fname(2), exist = found_it, name = fname(3))
    if (.not. found_It .or. fname(1) /= fname(3)) then
      if (bmad_status%type_out) then
        type *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE OUT OF DATE.'
      endif
      close (d_unit)
      bmad_status%ok = .false.
      return
    endif
  enddo

! we read (and write) the ring in pieces since it is
! too big to write in one piece

  read (d_unit, err = 9100) ring%parameters

  if (ring%n_ele_max > n_ele_maxx) then
    type *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                            ring%n_ele_max
    type *, '      IS GREATER THAN THE ELEMENT ARRAY SIZE:', n_ele_maxx
    type *, '      YOU NEED TO RECOMPILE'
    call err_exit
  endif

  do i = 0, ring%n_ele_max
    read (d_unit, err = 9100) ring%ele_(i)
  enddo

  if (ring%n_control_array > n_control_maxx) then
    type *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                      ring%n_control_array
    type *, '      IS GREATER THAN THE CONTROL ARRAY SIZE:', n_control_maxx
    type *, '      YOU NEED TO RECOMPILE'
    call err_exit
  endif

  do i = 1, ring%n_control_array
    read (d_unit, err = 9100) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_array
    read (d_unit, err = 9100) ring%ic_(i)
  enddo

  close (d_unit)
  return

!------------------

9000  continue
  if (bmad_status%type_out) then
    type *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE DOES NOT EXIST.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9100  continue
  if (bmad_status%type_out) then
    type *, 'READ_DIGESTED_BMAD_FILE: ERROR READING DIGESTED FILE.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

end subroutine
