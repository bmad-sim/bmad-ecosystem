!+
! subroutine GAUSSIAN_DIST (ELE, MODE, coupling, MIN_SIG, COORD_)
!
! Subroutine to generate gaussian distribution in 6d phase space in coord_struct 
!  For each mode (x,y,z), get gaussian distribution of amplitudes and then
!  choose phase space angle with flat random distribution 0-2pi
! Input:
!   ELE  -- Element_struct
!   MODE -- Mode_struct : emittances
!   MIN_SIG -- Real: only keep phase space with (N_xsig + N_ysig + N_zsig > MIN_SIG)
!   COUPLING -- Real: fraction of horizontal emittance that mixes into vertical
! Output:
!    Coord_ -- Coord_struct: Array of phase space coordinates
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

 subroutine gaussian_dist (ele, mode, coupling, min_sig, coord_)

  use bmad_struct
  use bmad_interface
  use nr, ONLY : gasdev_s, ran

  implicit none

  type (modes_struct) mode
  type (coord_struct), allocatable :: coord_(:)
  type (ele_struct) ele

  real(rdef) a, b, z

  real(rdef) theta_a, theta_b, theta_z, dl, delta, a_amp, b_amp, z_amp
  real(rdef) z_emit, beta_z
!  real(rdef) ran
  real(rdef) v(4,4), temp_vec(4)
  real(rdef) mat(6,6)
  real(rdef) mu, cosmu
  real(rdef) min_sig
  real(rdef) coupling

  integer i, idum,p
  integer j

  p=size(coord_)
  call out_io(s_info$,"GAUSSIAN_DIST",' size=\i\ ',p)


  V(1:4,1:4)=0
  forall (j=1:4) V(j,j)=ele%gamma_c
  V(1:2,3:4) = ele%c_mat
  V(3,1) = -ele%c_mat(2,2)
  V(4,2) = -ele%c_mat(1,1)
  V(3,2) = ele%c_mat(1,2)
  V(4,1) = ele%c_mat(2,1)
  

  i = 0
  do while(i < size(coord_)-1)

  call gasdev_s(a)
  call gasdev_s(b)
  call gasdev_s(z)

  if(abs(a) + abs(b) + abs(z) < min_sig)cycle
     i = i+1


     theta_a = ran(idum)*twopi
     theta_b = ran(idum)*twopi
     theta_z = ran(idum)*twopi
     a_amp = abs(a) * mode%a%emittance
     b_amp = abs(b) * sqrt(mode%b%emittance **2 + (mode%a%emittance * coupling)**2)
     z_amp = abs(z) * mode%sig_z * mode%sige_e

     dl =  sqrt(z_amp * ele%z%beta)* cos(theta_z) 
     delta = sqrt(z_amp /ele%z%beta)* sin(theta_z)

     coord_(i)%vec(1) = sqrt(a_amp * ele%x%beta) * cos(theta_a) + &
                          delta * ele%x%eta  
     coord_(i)%vec(2) = -sqrt(a_amp/ele%x%beta)*(ele%x%alpha * cos(theta_a) +sin(theta_a)) + &
                          delta * sin(theta_z) * ele%x%etap

     coord_(i)%vec(3) = sqrt(b_amp * ele%y%beta) * cos(theta_b) + &
                          delta * ele%y%eta                          
     coord_(i)%vec(4) = -sqrt(b_amp/ele%y%beta)*(ele%y%alpha * cos(theta_b) +sin(theta_b)) + &
                          delta * ele%y%etap  

     coord_(i)%vec(5) = dl
     coord_(i)%vec(6) = delta

     temp_vec(1:4) = matmul(V,coord_(i)%vec(1:4))  !switch lab coordinates
     coord_(i)%vec(1:4) = temp_vec(1:4)

   end do

  return

end subroutine gaussian_dist
                                                            












