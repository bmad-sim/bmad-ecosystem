!+
! subroutine CLOSE_PRETZEL (ring, final_pos_in, final_pos_out, i_dim )
!
! Subroutine to eliminate differential displacement of strong and weak beam at IP
! by adjustment of PRETZING 13
!
! Input:
!   RING -- Ring_struct: Ring 
!   FINAL_POS_IN -- Coord_struct: optional, position to which closed orbit should converge
!   FINAL_POS_OUT -- Coord_struct: optional, actual position to which closed orbit has converged
! Output:
!    RING -- Ring_struct: Ring  with separators adjusted to bring beams into collision
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
! Revision 1.2  2006/04/18 22:16:15  ajl59
! added "use bmadz_interface" to stop fail out at beambeam_separation's new optional argument
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

subroutine close_pretzel (ring, i_dim, final_pos_in, final_pos_out)
  
  use bmad_interface
  use bmadz_mod
  use bmadz_interface

  implicit none
  
  type ( ring_struct ) ring
  type (ring_struct), save :: ring_oppos
  type (coord_struct), allocatable, save :: co_(:), co_oppos_(:)
  type (coord_struct) delta_ip, delta_ip_p
  type (coord_struct) delta_ip_0
  type (ele_struct) ele

  integer particle
  integer ierr
  integer i_dim
  integer i,ix_w/0/, ix_e/0/
  integer n
 

  real(rdef) dx, dv, dv_0
  real(rdef) kick_w, kick_e

  character(20) :: r_name='close_pretzel'
  character(200) write_line
 
!ANDREW
  type(coord_struct), optional:: final_pos_in, final_pos_out
  
  call reallocate_coord( co_, ring%n_ele_max )
  call reallocate_coord( co_oppos_, ring%n_ele_max )
  
  co_(0)%vec=0
  call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
  call track_all(ring, co_)

! find separators
  i=0
  do while ((ix_w == 0 .or. ix_e == 0) .and. i <= ring%n_ele_max)
     i=i+1
     if(ring%ele_(i)%name == 'H_SEP_08W') ix_w=i
     if(ring%ele_(i)%name == 'H_SEP_08E') ix_e=i
  enddo
  if(ix_w == 0 .or. ix_e == 0)then
     call out_io(s_abort$,r_name," Cannot find separators")
     stop
  endif
!    type *,ring%ele_(ix_w)%name,ring%ele_(ix_w)%value(hkick$) 
!    type *,ring%ele_(ix_e)%name,ring%ele_(ix_e)%value(hkick$) 

  

  
  call ring_make_mat6(ring, ix_w)
  call ring_make_mat6(ring, ix_e)
  call beambeam_separation(ring, delta_ip_0, i_dim)
  
  n=0
  
  write (write_line,'(1x,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)') n, &
       ' 8W(mr) = ',ring%ele_(ix_w)%value(hkick$)*1000., &
       ' 8E(mr) = ',ring%ele_(ix_e)%value(hkick$)*1000., &
       '   dx,dxp,dy,dyp (mm) = ',delta_ip_0%vec(1:4)*1000.  
  call out_io(s_info$,r_name, write_line)
  
  if(.not. present(final_pos_in)) then
     final_pos_in%vec(:) = 0.
     final_pos_out%vec(:) = 0
  end if
       
  do while (abs(delta_ip_0%vec(1)-final_pos_in%vec(1)) > 1.e-6 .and. n < 4)
     kick_w =  ring%ele_(ix_w)%value(hkick$)
     kick_e =  ring%ele_(ix_e)%value(hkick$)
     
     n = n+1
     
!  type *,' INITIAL SEPARATION :', delta_ip_0%vec(1)
     dv_0=-0.00001/n
     ring%ele_(ix_w)%value(hkick$) = kick_w + dv_0
     ring%ele_(ix_e)%value(hkick$) = kick_e + dv_0
     call ring_make_mat6(ring, -1)
     call beambeam_separation(ring, delta_ip_p, i_dim)
     
     dx = (delta_ip_p%vec(1) - (delta_ip_0%vec(1)))/dv_0
     dv = -(delta_ip_0%vec(1)-final_pos_in%vec(1))/dx
     
     ring%ele_(ix_w)%value(hkick$) = kick_w + dv
     ring%ele_(ix_e)%value(hkick$) = kick_e + dv
     
     call ring_make_mat6(ring, ix_w)
     call ring_make_mat6(ring, ix_e)
     call beambeam_separation(ring, delta_ip, i_dim)
     
     write (write_line,'(1x,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)') n, &
          ' 8W(mr) = ',ring%ele_(ix_w)%value(hkick$)*1000., &
          ' 8E(mr) = ',ring%ele_(ix_e)%value(hkick$)*1000., &
          '   dx,dxp,dy,dyp (mm) = ',delta_ip_0%vec(1:4)*1000.  
     call out_io(s_info$,r_name, write_line)
     
     
 

!  type *,' kick_w, kick_e ', kick_w , kick_e
!  type *,' INITIAL SEPARATION :', delta_ip_0%vec(1)
!  type *,' kick_w, kick_e ', kick_w + dv_0, kick_e +dv_0
!  type *,' derivative  separation: ', delta_ip_p%vec(1) - delta_ip_0%vec(1)
!  type *,' kick_w, kick_e ', kick_w + dv, kick_e +dv
!  type *,' Final separation: ', delta_ip%vec(1)

     delta_ip_0%vec = delta_ip%vec
  end do
  
  final_pos_out%vec(1:4)=delta_ip_0%vec(1:4)

  if(n>= 4)call out_io(s_warn$,r_name,' PRETZEL NOT CLOSED')
  
  
  
  return
  
end subroutine close_pretzel

