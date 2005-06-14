!........................................................................
!+
! module scan_parameters
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

module scan_parameters

  use bmad

  type scan_params_struct
    character*80 lat_file
    character*80 file_name
    real(rdef) Q_z
    real(rdef) current  !mA/bunch
    real(rdef) sig_in(3)   !sigx, sigy, sigz, initial distribution
    real(rdef) sig_out(3)  !sigx, sigy, sigz, final distribution
    real(rdef) min_sig     !x/sigx + y/sigy + z/sigz > min_sig
    real(rdef) lum         !luminosity
    real(rdef) coupling_wb    !coupling of horizontal emittance into vertical weak beam
    real(rdef) coupling_sb    !coupling of horizontal emittance into vertical strong beam
    integer n_turn, particle, i_train, j_car, n_trains_tot, n_cars 
    integer n_part !initial distribution
    integer slices 
    integer n_part_out !final distribution
    logical lrbbi, beambeam_ip, close_pretz, close_vert, rec_taylor
    logical radiation
    character*5 fit
    type(coord_struct) final_pos_in ! gives final postion to which closed orbit will converge
    type(coord_struct) init(100)  !initial coordinate for single particle scan
    logical parallel !defines whether code is run in parallel
 end type scan_params_struct
end module
