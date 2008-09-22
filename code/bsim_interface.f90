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
! Revision 1.6  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.5.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.5  2005/09/21 20:59:07  dcs
! more changes to get around compiler bug.
!
! Revision 1.4  2005/09/21 20:35:58  dcs
! Put beambeam_setup in scan_parameters module.
!
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
#include "CESR_platform.inc"

module bsim_interface
  
  interface
     subroutine close_pretzel(ring,i_dim,final_pos_in, final_pos_out)
       use bmad_struct, only: lat_struct, coord_struct
       implicit none
       type (lat_struct), intent(inout) :: ring
       type (coord_struct), optional, intent(in) :: final_pos_in
       type (coord_struct), optional, intent(out) :: final_pos_out
       integer, intent(in) :: i_dim
     end subroutine close_pretzel
  end interface
  
  interface
     subroutine close_vertical(ring,i_dim,final_pos_in, final_pos_out)
       use bmad_struct, only: lat_struct, coord_struct
       implicit none
       type(lat_struct), intent(inout) :: ring
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
       character(*), intent(in) ::  in_file
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
     subroutine gaussian_dist (ele, mode, coupling, min_sig, coord)
       use bmad_struct, only: normal_modes_struct, coord_struct, ele_struct, rp
       implicit none
       type (normal_modes_struct) mode
       type (coord_struct), allocatable :: coord(:)
       type (ele_struct) ele
       real(rp) min_sig
       real(rp) coupling
     end subroutine gaussian_dist
  end interface
  
  interface
     subroutine histogram (ele, coord, in_file, sig,a_out)
       use bmad_struct, only: coord_struct, ele_struct, rp
       implicit none
       type (coord_struct) coord(:)
       type (ele_struct) ele
       real(rp) sig(3), a_out(3)
       character(*) in_file
     end subroutine histogram
  end interface
  
  interface
     subroutine histogram_new (ele, coord, in_file, sig)
       use bmad_struct, only: coord_struct, ele_struct, rp
       implicit none
       type (coord_struct) coord(:)
       type (ele_struct) ele
       real(rp) sig(3)
       character*60 in_file
     end subroutine histogram_new
  end interface
  
  interface
     subroutine luminosity_calc (ele, coord, param, n_ok, lum)
       use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, rp
       implicit none
       type(ele_struct) ele
       type(coord_struct), allocatable :: coord(:)
       type(lat_param_struct) param
       real(rp) lum, f
       integer n_ok
     end subroutine luminosity_calc
  end interface
  
  interface
     subroutine lum_tracker(ring,n_part, start, end)
       use bmad_struct, only: lat_struct, coord_struct
       implicit none
       type(lat_struct) ring
       type(coord_struct) start(:),end(:)
       type(coord_struct), allocatable, save :: co(:)
       integer n_part, i, j
     end subroutine lum_tracker
  end interface
  
  interface
     subroutine MARK_LRBBI_ONLY(master_ring, master_ring_oppos, ring, crossings)
       use bmad_struct, only: lat_struct, rp
       implicit none
       type (lat_struct), dimension(:) :: ring
       type (lat_struct) :: master_ring, master_ring_oppos
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
