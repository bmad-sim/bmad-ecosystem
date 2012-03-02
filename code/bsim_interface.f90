module bsim_interface
  
  interface
     subroutine writefile(in_file, parameters)
       use precision_def
       implicit none
       real(rp), dimension(1:,1:), intent(in) ::  parameters
       character(*), intent(in) ::  in_file
     end subroutine writefile
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
     subroutine save_last_pturns(pturns, end, n_part, turn)
       use bmad_struct, only: coord_struct
       type (coord_struct)  end(1:)
       integer pturns, n_part, turn
     end subroutine save_last_pturns
  end interface
  
 interface
  subroutine implement_pathlength_patch(path_length_patch,ring, delta_frf, frf) 
   use bmad_struct
   use bmad_interface
   implicit none
   type (lat_struct) ring
   real(rp), optional :: delta_frf, frf
   logical path_length_patch
  end subroutine implement_pathlength_patch
 end interface

integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module bsim_interface


