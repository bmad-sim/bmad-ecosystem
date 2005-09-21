!........................................................................
!+
! module bsim_interface
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
! Revision 1.3  2005/09/21 20:19:06  dcs
! another try to get around compiler bug.
!
! Revision 1.2  2005/09/20 20:25:22  dcs
! cleaned up use statements.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

module bsim_interface
  
  interface
     subroutine close_pretzel(ring,i_dim,final_pos_in, final_pos_out)
       use bmad_struct, only: ring_struct, coord_struct
       implicit none
       type (ring_struct), intent(inout) :: ring
       type (coord_struct), optional, intent(in) :: final_pos_in
       type (coord_struct), optional, intent(out) :: final_pos_out
       integer, intent(in) :: i_dim
     end subroutine close_pretzel
  end interface
  
  interface
     subroutine close_vertical(ring,i_dim,final_pos_in, final_pos_out)
       use bmad_struct, only: ring_struct, coord_struct
       implicit none
       type(ring_struct), intent(inout) :: ring
       type(coord_struct), optional, intent(in) :: final_pos_in
       type(coord_struct), optional, intent(out) :: final_pos_out
       integer, intent(in) :: i_dim
     end subroutine close_vertical
  end interface
  
  
  interface
     subroutine read_turns()
       implicit none
     end subroutine read_turns
  end interface
  
  interface
     subroutine writefile(in_file, parameters)
       use precision_def
       implicit none
       real(rp), dimension(1:,1:), intent(in) ::  parameters
       character*60, intent(in) ::  in_file
     end subroutine writefile
  end interface
  
  interface
     subroutine gfit3d(phase_coord,parameters)
       use bmad_struct, only: coord_struct, rp
       implicit none
       type(coord_struct), dimension(1:), intent(in) :: phase_coord
       real(RP), dimension(1:,1:), intent(inout):: parameters
     end subroutine gfit3D
  end interface
  
  interface
     subroutine beambeam_setup (ring, particle,  current, scan_params, slices)
       use bmad_struct, only: ring_struct, rp
       use scan_parameters, only: scan_params_struct
       implicit none
       type ( ring_struct ) ring
       type (scan_params_struct) scan_params
       integer particle
       integer, optional, intent(in) :: slices
       real(rp) current
     end subroutine beambeam_setup
  end interface
  
  
  interface
     subroutine gaussian_dist (ele, mode, coupling, min_sig, coord_)
       use bmad_struct, only: modes_struct, coord_struct, ele_struct, rp
       implicit none
       type (modes_struct) mode
       type (coord_struct), allocatable :: coord_(:)
       type (ele_struct) ele
       real(rp) min_sig
       real(rp) coupling
     end subroutine gaussian_dist
  end interface
  
  interface
     subroutine histogram (ele, coord_, in_file, sig,a_out)
       use bmad_struct, only: coord_struct, ele_struct, rp
       implicit none
       type (coord_struct) coord_(:)
       type (ele_struct) ele
       real(rp) sig(3), a_out(3)
       character*60 in_file
     end subroutine histogram
  end interface
  
  interface
     subroutine histogram_new (ele, coord_, in_file, sig)
       use bmad_struct, only: coord_struct, ele_struct, rp
       implicit none
       type (coord_struct) coord_(:)
       type (ele_struct) ele
       real(rp) sig(3)
       character*60 in_file
     end subroutine histogram_new
  end interface
  
  interface
     subroutine luminosity_calc (ele, coord_, param, n_ok, lum)
       use bmad_struct, only: ele_struct, coord_struct, param_struct, rp
       implicit none
       type(ele_struct) ele
       type(coord_struct), allocatable :: coord_(:)
       type(param_struct) param
       real(rp) lum, f
       integer n_ok
     end subroutine luminosity_calc
  end interface
  
  interface
     subroutine lum_tracker(ring,n_part, start, end)
       use bmad_struct, only: ring_struct, coord_struct
       implicit none
       type(ring_struct) ring
       type(coord_struct) start(:),end(:)
       type(coord_struct), allocatable, save :: co_(:)
       integer n_part, i, j
     end subroutine lum_tracker
  end interface
  
  interface
     subroutine MARK_LRBBI_ONLY(master_ring, master_ring_oppos, ring, crossings)
       use bmad_struct, only: ring_struct, rp
       implicit none
       type (ring_struct), dimension(:) :: ring
       type (ring_struct) :: master_ring, master_ring_oppos
       real(rp), dimension(:,:) :: crossings
     end subroutine MARK_LRBBI_ONLY
  end interface
  
  interface
     subroutine save_last_pturns(pturns, end, n_part, turn)
       use bmad_struct, only: coord_struct
       type (coord_struct)  end(1:)
       integer pturns, n_part, turn
     end subroutine save_last_pturns
  end interface
  
end module bsim_interface
