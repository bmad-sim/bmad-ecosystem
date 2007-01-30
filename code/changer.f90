!........................................................................
!+
! Subroutine : subroutine changer (ring)
!
! Description:
!
! Arguments  : ring
!
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"
   subroutine changer (ring)

     use bmad_struct
     use bmad_interface

     implicit none

     type(lat_struct) ring

     character*80 line, last_line
     integer ix


     do while(.true.)
     print '(a, $)', ' BEAMBEAM: element change or GO> '
     read(5, '(a)') line
    
      ix = index(line, '!')
      if (ix /= 0) line = line(:ix-1)        ! strip off comments

      call str_upcase(line, line)
      call string_trim(line, line, ix)

      if (ix == 0) then       ! nothing typed. do the same thing
        line = last_line
      endif

      last_line = line

      if(line(1:1) .eq. 'G')exit

      call find_change( line, ring)
     end do
    return
 end

