!........................................................................
!+
!  module    : beambeam_interface 
!
! Description:
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
! Revision 1.1  2005/06/14 14:59:02  cesrulib
! Initial revision
!
!
!........................................................................
!
#include "CESR_platform.h"

module beambeam_interface

  interface
     subroutine size_beam(ring, end, scan_params, transmit, sib_j, n_typeout, orb_, phi_x, phi_y, past_params, past_lums, parameters)

       use bmad_struct
       use bmad_interface
       use bmadz_mod
       use bmadz_interface
       use scan_parameters
       use bsim_interface
       
       implicit none
       
       type(ring_struct) ring
       type(coord_struct) end(:)
       type(scan_params_struct) scan_params
       type(coord_struct) orb_(0:)
       
       logical, dimension(2) :: transmit
       
       integer sib_j, n_typeout
       
       real(rdef) phi_x, phi_y
       real(rdef) :: past_params(:,:), past_lums(:), parameters(:,:)
     end subroutine size_beam
  end interface

end module beambeam_interface
