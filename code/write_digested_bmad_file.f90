!+
! Subroutine write_digested_bmad_file (digested_name, ring, n_files, file_names)
!
! Subroutine to write a digested file. The names of the original files used
! to create the RING structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Modules Needed:
!   use bmad
!
! Input:
!     digested_name -- Character(*): Name for the digested file.
!     ring          -- Ring_struct: Input ring structure.
!     n_files       -- Number of original files
!     file_names(*) -- Character(*), optional: Names of the original 
!                       files used to create the ring structure.
!-

!$Id$
!$Log$
!Revision 1.8  2002/11/07 17:10:04  dcs
!Bug_fix
!
!Revision 1.7  2002/11/06 06:48:32  dcs
!Changed arg array
!
!Revision 1.6  2002/06/13 14:54:31  dcs
!Interfaced with FPP/PTC
!
!Revision 1.5  2002/02/23 20:32:30  dcs
!Double/Single Real toggle added
!
!Revision 1.4  2002/01/08 21:44:44  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.3  2001/10/05 18:23:57  rwh24
!Bug Fixes
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine write_digested_bmad_file (digested_name, ring,  &
                                                  n_files, file_names)

  use bmad

  implicit none

  type (ring_struct), target, intent(in) :: ring
  type (ele_digested_struct) :: u_ele
  type (ele_struct), pointer :: ele
  type (taylor_struct), pointer :: tt(:)
  
  integer d_unit, lunget, n_files, i, j, k, ix_w, ix_d, ix_m, ix_t(6)
  integer stat_b(12), stat, ierr

  character(*) digested_name
  character(*), optional :: file_names(:)
  character(200) fname

  external stat

! write input file names to the digested file

  d_unit = lunget()
  open (unit = d_unit, file = digested_name, form = 'unformatted', err = 9000)
  write (d_unit, err = 9010) n_files, bmad_inc_version$
  do j = 1, n_files
    fname = file_names(j)
    stat_b = 0
#ifdef CESR_UNIX
    ierr = stat(fname, stat_b)
#endif
    write (d_unit) fname, stat_b(10)  ! stat_b(10) = Modification date
  enddo

! write the ring structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.
! also: because of the pointers we need to disguise that we are writting ele's
!  by using the digested_ele_struct

  write (d_unit) ring%parameters
  
  do i = 0, ring%n_ele_max
  
    ele => ring%ele_(i)
    tt => ele%taylor
    u_ele%ele = ele
    
    ix_w = 0; ix_d = 0; ix_m = 0; ix_t = 0

    if (ele%pointer_init == has_been_inited$) then
      if (associated(ele%wig_term)) ix_w = size(ele%wig_term)
      if (associated(ele%descrip))  ix_d = 1
      if (associated(ele%a))        ix_m = 1
      if (associated(tt(1)%term))   ix_t = (/ (size(tt(j)%term), j = 1, 6) /)
    endif

    write (d_unit) u_ele%digested, ix_w, ix_d, ix_m, ix_t
    do j = 1, ix_w
      write (d_unit) ele%wig_term(j)
    enddo

    if (associated(ele%descrip))  write (d_unit) ele%descrip
    if (associated(ele%a))        write (d_unit) ele%a, ele%b
    
    do j = 1, 6
      do k = 1, ix_t(j)
        write (d_unit) tt(j)%term(k)
      enddo
    enddo

  enddo

! write the control info

  do i = 1, ring%n_control_array
    write (d_unit) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_array
    write (d_unit) ring%ic_(i)
  enddo

  close (d_unit)

  return

! Errors

9000  type *
  type *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT OPEN FILE FOR OUTPUT:'
  type *, '    ', trim(digested_name)
  type *, '     [This does not affect program operation]'
  return

9010  type *
  type *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT WRITE TO FILE FOR OUTPUT:'
  type *, '    ', trim(digested_name)
  type *, '     [This does not affect program operation]'
  close (d_unit)
  return

end subroutine
