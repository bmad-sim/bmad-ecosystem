!+
! Subroutine CESR_ELEMENTS_GET (NAME, N_FOUND, ELE)
!
! Subroutine to find the location of elements from [CESR.SURVEY]RING_MASTER.DAT
!
! Modules needed:
!     use bmad_struct
!     use bmad_interface
!
! Input:
!     NAME  -- character*(*). Name of element(s) to be found. NAME may contain
!              wild cards '*' and '%' with the normal wild card matching rules.
!              For example, NAME = 'Q*E' will find all the east quads.
!              The subroutine matches NAME with names found in RING_MASTER.DAT.
!
! Output:
!     N_FOUND -- Integer. Number of elements found.
!     ELE(:)  -- Ring_master_struct: Structure for the elements
!      %NAME     -- Character*6: Name of element. If NAME does not
!                    contain a wild card then, ELE(1).NAME = NAME.
!      %Z        -- Real*8: Distance from the IP to the center of the element.
!      %MAG_LEN  -- Real*8: Magnetic length of the element.
!      %PHYS_LEN -- Real*8: Physical length  of the element.
!      %COMMENT  -- Character*16. Comment for the element
!-

!$Id$
!$Log$
!Revision 1.1  2002/06/13 14:54:23  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:34:40  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.3  2001/10/02 18:50:51  rwh24
!More compatibility updates.
!
!Revision 1.2  2001/09/27 17:47:02  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine cesr_elements_get (name, n_found, ele)

  use bmad_struct
  use cesr_utils


  use precision_def

  implicit none

  type (ring_master_struct) :: ele(:)

  integer lun, n_found, ix, ixn, ix1, ix2, match_wild
  character name*(*), name_in*10, line*79
  logical finding

!

  n_found = 0
  call string_trim (name, name_in, ixn) ! strip off leading blanks
  call str_upcase (name_in, name_in)

  lun = lunget()
  open (unit = lun, file = '[cesr.survey]ring_master.dat', shared, &
                                             status = 'old', readonly)

  finding = .true.
  do while (finding)
    read (lun, '(a)')  line
    if (line(1:6) == '! Name') finding = .false.
  enddo

! Read a line from ring_master.dat and parse

  do while (.true.)

    read (lun, '(a)', end = 8000) line
    call string_trim (line, line, ix)
    if (ix /= 0) then
      

      if (match_wild (line(:ix), name_in(:ixn))) then  ! yes, a match
        n_found = n_found + 1
        ele(n_found)%name = line(:ix)
        read (line(ix+1:), *, err = 9000) ele(n_found)%z, &
                      ele(n_found)%mag_len, ele(n_found)%phys_len
        ix1 = index (line, '"')
        ix2 = index (line(ix1+1:), '"')
        if (ix1 == 0 .or. ix2 == 0) then
          type *, 'ERROR IN CESR_ELEMENTS_GET: NO COMMENT FOUND FOR', line
          call err_exit
        endif
        if (ix2 == 1) then
          ele(n_found)%comment = ' '
        else
          ele(n_found)%comment = line(ix1+1:ix1+ix2-1)
        endif
      endif
    endif

  enddo


8000  continue
  close (unit = lun)
  return

9000  continue
  close (unit = lun)
  type *
  type *, 'ERROR IN CESR_ELEMENTS_GET: CANNOT READ DATA IN RING_MASTER.DAT FOR ELEMENT:'
  type *, line
  call err_exit

  end
