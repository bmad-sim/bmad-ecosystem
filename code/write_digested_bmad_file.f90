!+
! Subroutine write_digested_bmad_file (digested_name, ring, n_files, file_names)
!
! Subroutine to write a digested file. The names of the original files used
! to create the RING structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     digested_name -- Character*72: Name for the digested file.
!     ring          -- Ring_struct: Input ring structure.
!     n_files       -- Number of original files
!     file_names(*) -- Character*72: Names of the original files used to create
!                       the ring structure.
!-

!$Id$
!$Log$
!Revision 1.3  2001/10/05 18:23:57  rwh24
!Bug Fixes
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine write_digested_bmad_file (digested_name, ring,  &
                                                  n_files, file_names)

  use bmad_struct
  use local_bmad_struct

  implicit none

  type (ring_struct)   ring

  integer d_unit, lunget, n_files, i, j

  character*(*) digested_name, file_names(*)
  character*200 fname

! write input file names to the digested file

  d_unit = lunget()
  open (unit = d_unit, file = digested_name, status = 'new',  &
                                          form = 'unformatted', err = 9000)
  write (d_unit) n_files, bmad_inc_version$
  do j = 1, n_files
    fname = file_names(j)
    write (d_unit) fname
  enddo

! write the ring structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.

  write (d_unit) ring%parameters
  do i = 0, ring%n_ele_max
    write (d_unit) ring%ele_(i)
  enddo

  do i = 1, ring%n_control_array
    write (d_unit) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_array
    write (d_unit) ring%ic_(i)
  enddo

! write the common strucures for bmad_parser and bmad_parser2
! ** this is obsolete ** 10/9/2000

  write (d_unit) pcom

  do i = 1, pcom%ivar_tot
    write (d_unit) var_(i)
  enddo

  close (d_unit)

  return

9000  type *
  type *, 'WRITE_DIGESTED_BMAD_FILE: ERROR OPENING FILE FOR OUTPUT: '
  type *, '    ', digested_name
  return

end subroutine
