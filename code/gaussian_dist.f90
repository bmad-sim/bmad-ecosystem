!+
! subroutine GAUSSIAN_DIST (ELE, MODE, coupling, MIN_SIG, COORD)
!
! Subroutine to generate gaussian distribution in 6d phase space in coord_struct 
!  For each mode (x,y,z), get gaussian distribution of amplitudes and then
!  choose phase space angle with flat random distribution 0-2pi
! Input:
!   ELE  -- Element_struct
!   MODE -- normal_mode_struct : emittances
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
 subroutine gaussian_dist (ele, mode, coupling, min_sig, coord)

  use bmad_interface
  use nr, ONLY : gasdev_s, ran

  implicit none

  type (normal_modes_struct) mode
  type (coord_struct), allocatable :: coord(:)
  type (ele_struct) ele

  real(rp) a, b, z

  real(rp) theta_a, theta_b, theta_z, dl, delta, a_amp, b_amp, z_amp
  real(rp) z_emit, beta_z
!  real(rp) ran
  real(rp) v(4,4), temp_vec(4)
  real(rp) mat(6,6)
  real(rp) mu, cosmu
  real(rp) min_sig
  real(rp) coupling

  integer i, idum,p
  integer j

  p=size(coord)
  call out_io(s_info$,"GAUSSIAN_DIST",' size=\i\ ',p)


  V(1:4,1:4)=0
  forall (j=1:4) V(j,j)=ele%gamma_c
  V(1:2,3:4) = ele%c_mat
  V(3,1) = -ele%c_mat(2,2)
  V(4,2) = -ele%c_mat(1,1)
  V(3,2) = ele%c_mat(1,2)
  V(4,1) = ele%c_mat(2,1)
  

  i = 0
  do while(i < size(coord)-1)

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

     coord(i)%vec(1) = sqrt(a_amp * ele%a%beta) * cos(theta_a) + &
                          delta * ele%a%eta  
     coord(i)%vec(2) = -sqrt(a_amp/ele%a%beta)*(ele%a%alpha * cos(theta_a) +sin(theta_a)) + &
                          delta * sin(theta_z) * ele%a%etap

     coord(i)%vec(3) = sqrt(b_amp * ele%b%beta) * cos(theta_b) + &
                          delta * ele%b%eta                          
     coord(i)%vec(4) = -sqrt(b_amp/ele%b%beta)*(ele%b%alpha * cos(theta_b) +sin(theta_b)) + &
                          delta * ele%b%etap  

     coord(i)%vec(5) = dl
     coord(i)%vec(6) = delta

     temp_vec(1:4) = matmul(V,coord(i)%vec(1:4))  !switch lab coordinates
     coord(i)%vec(1:4) = temp_vec(1:4)

   end do

  return

end subroutine gaussian_dist
                                                            












