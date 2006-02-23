!+
! subroutine CLOSE_VERTICAL (ring, final_pos_in, final_pos_out, i_dim )
!
! Subroutine to eliminate differential displacement and angle
! of strong and weak beam at IP
! by adjustment of vertical separators
!
! Input:
!   RING -- Ring_struct: Ring 
!   FINAL_POS_IN -- Coord_struct: optional, desired position to which closed orbit should converge
!   FINAL_POS_OUT -- Coord_struct: optional, actual position to which closed orbit has converged
!   I_DIM -- 4 or 6 d tracking
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
! Revision 1.2  2006/02/23 16:22:49  ajl59
! added "use bmadz_interface" to stop fail out at beambeam_separation's new optional argument
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

subroutine close_vertical(ring, i_dim, final_pos_in, final_pos_out)
  
  use bmad
  use bmadz_mod
  use bmadz_interface
  
  implicit none
  
  type (ring_struct) ring
  type (ring_struct), save :: ring_oppos
  type (coord_struct), allocatable, save :: co_(:), co_oppos_(:)
  type (coord_struct) delta_ip_0, delta_ip, delta_ip_p
  type (ele_struct) ele

  integer particle
  integer ierr
  integer i,ix_w/0/, ix_e/0/, ix_7/0/
  integer n
  integer i_dim

  real(rdef) dy_5, dyp_5, dy_7, dyp_7, dv(2), dv_0, dvc7, vc7
  real(rdef) kick_w, kick_e
  real(rdef) mat_inv(2,2)
  real(rdef) det

  character(20) :: r_name='close_vertical'
  character(200) write_line

!ANDREW
!  integer, parameter :: direction=1   ! determines forward (1) or backward (-1) tracking in closed_orbit_calc, optional
  type (coord_struct), optional :: final_pos_in
  type (coord_struct), optional :: final_pos_out


  call reallocate_coord(co_,ring%n_ele_max)
  call reallocate_coord(co_oppos_,ring%n_ele_max)
  
  co_(0)%vec=0
  call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
! routine CLOSED_ORBIT_AT_START is obselete, replacing with:
!  call closed_orbit_calc (ring, co_(0), i_dim, direction)
  call track_all(ring, co_)
  
! find vertical separators
  i=0
  do while ((ix_w == 0 .or. ix_e == 0) .and. i <= ring%n_ele_max)
     i=i+1
     if(ring%ele_(i)%name == 'V_SEP_48W') ix_w=i
     if(ring%ele_(i)%name == 'V_SEP_48E') ix_e=i
  enddo
  if(ix_w == 0 .or. ix_e == 0)then
     call out_io(s_abort$,r_name," Cannot find separators")
     stop
  endif
  
! find vertical separators
  i = ring%n_ele_ring
  do while (ix_7 == 0 .and. i <= ring%n_ele_max)
     i=i+1
     if(ring%ele_(i)%name == 'RAW_VCROSING_7') ix_7=i
  enddo
  if(ix_7 == 0)then
     call out_io(s_abort$,r_name," Cannot find RAW_VCROSING7")      
     stop
  endif
!    type *,ring%ele_(ix_w)%name,ring%ele_(ix_w)%value(vkick$) 
!    type *,ring%ele_(ix_e)%name,ring%ele_(ix_e)%value(vkick$) 
  
  
  
  
  call ring_make_mat6(ring, ix_w)
  call ring_make_mat6(ring, ix_e)
  call beambeam_separation(ring, delta_ip_0, i_dim)
  
  n=0
  
  write (write_line,'(1x,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)') n, &
       ' 48W(mr) = ',ring%ele_(ix_w)%value(vkick$)*1000., &
       ' 48E(mr) = ',ring%ele_(ix_e)%value(vkick$)*1000., &
       '   dx,dxp,dy,dyp (mm) = ',delta_ip_0%vec(1:4)*1000.
  
  call out_io(s_info$,r_name, write_line)
  call out_io(s_info$,r_name,' RAW_VCROSING_7 = \f\ ' , ring%ele_(ix_7)%value(command$))
  
!  type '(a16,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)',' CLOSE_VERT: ',n, &
!         ' 48W(mr) = ',ring%ele_(ix_w)%value(vkick$)*1000., &
!         ' 48E(mr) = ',ring%ele_(ix_e)%value(vkick$)*1000., &
!          '   dx,dxp,dy,dyp (mm) = ',delta_ip_0%vec(1:4)*1000.
!  type '(a16,f)',' CLOSE_VERT: RAW_VCROSING_7 = ', ring%ele_(ix_7)%value(command$)  

  if(.not. present(final_pos_in)) then
     final_pos_in%vec(:) = 0.
     final_pos_out%vec(:) = 0
  end if

  do while ((abs(delta_ip_0%vec(3)-final_pos_in%vec(3)) > 1.e-6 .or. abs(delta_ip_0%vec(4)-final_pos_in%vec(4)) &
       > 1.e-6) .and. n < 10)
! Set kicks on separators
     kick_w =  ring%ele_(ix_w)%value(vkick$)
     kick_e =  ring%ele_(ix_e)%value(vkick$)
! Set knob
     vc7 = ring%ele_(ix_7)%value(command$)
     
     n = n+1
     
     !  type *,' INITIAL SEPARATION :', delta_ip_0%vec(1)
! Set increment for kick and knob
     dv_0=-0.00001
     dvc7=0.1
     
! Change kick
     ring%ele_(ix_w)%value(vkick$) = kick_w + dv_0
     ring%ele_(ix_e)%value(vkick$) = kick_e - dv_0
! Implement Changes
     call ring_make_mat6(ring, -1)
! Calculate final position of particle after orbit
     call beambeam_separation(ring, delta_ip_p, i_dim)
     
! Calculate derivatives
     dy_5 = (delta_ip_p%vec(3) - delta_ip_0%vec(3))/dv_0
     dyp_5 = (delta_ip_p%vec(4) - delta_ip_0%vec(4))/dv_0
     
! Reset kicks
     ring%ele_(ix_w)%value(vkick$) = kick_w
     ring%ele_(ix_e)%value(vkick$) = kick_e
! Adjust knob
     ring%ele_(ix_7)%value(command$) = vc7+dvc7
! Implement changes
     call ring_make_mat6(ring, -1)
! Calculate final position of particle after orbit
     call beambeam_separation(ring, delta_ip_p, i_dim)
     
     
! Calculate derivatives
     dy_7 = (delta_ip_p%vec(3) - delta_ip_0%vec(3))/dvc7
     dyp_7 = (delta_ip_p%vec(4) - delta_ip_0%vec(4))/dvc7
     
! Invert matrix of derivatives [dy/dv, dy'/dv; dy/dv7, dy'/dv7]
     det = dy_5*dyp_7 - dyp_5*dy_7
     mat_inv(1,1) = dyp_7/det
     mat_inv(1,2) = -dy_7/det
     mat_inv(2,1) = -dyp_5/det
     mat_inv(2,2) = dy_5/det
     
! Calculate changes that must be made to kick & knob
! (dv, dv7) = mat_inv*(y,y')  ; matrix inverse * initial position
     dv = -matmul(mat_inv, delta_ip_0%vec(3:4)-final_pos_in%vec(3:4) ) 
     
! Change kick and knob as necessary
     ring%ele_(ix_w)%value(vkick$) = kick_w + dv(1)
     ring%ele_(ix_e)%value(vkick$) = kick_e - dv(1)
     ring%ele_(ix_7)%value(command$) = vc7 + dv(2) 
     
! Implement changes
     call ring_make_mat6(ring, -1)
! Calculate final position of particle after orbit which should now be near final_pos
     call beambeam_separation(ring, delta_ip, i_dim)
     
     write (write_line,'(1x,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)') n, &
          ' 48W(mr) = ',ring%ele_(ix_w)%value(vkick$)*1000., &
          ' 48E(mr) = ',ring%ele_(ix_e)%value(vkick$)*1000., &
          '   dx,dxp,dy,dyp (mm) = ',delta_ip_0%vec(1:4)*1000.
     
     call out_io(s_info$,r_name, write_line)
     call out_io(s_info$,r_name," RAW_VCROSING_7 = \f\ ",ring%ele_(ix_7)%value(command$))
     
!  type '(a16,i2,1x,a10,f6.3,a10,f6.3,a24,4f8.4)',' CLOSE_VERT: ',n, &
!          ' 48W(mr) = ',ring%ele_(ix_w)%value(vkick$)*1000., &
!         ' 48E(mr) = ',ring%ele_(ix_e)%value(vkick$)*1000., &
!          '   dx,dxp,dy,dyp (mm) = ',delta_ip%vec(1:4)*1000.  
!  type '(a16,f)',' CLOSE_VERT: RAW_VCROSING_7 = ', ring%ele_(ix_7)%value(command$)  

!  type *,' kick_w, kick_e ', kick_w , kick_e
!  type *,' INITIAL SEPARATION :', delta_ip_0%vec(1)
!  type *,' kick_w, kick_e ', kick_w + dv_0, kick_e +dv_0
!  type *,' derivative  separation: ', delta_ip_p%vec(1) - delta_ip_0%vec(1)
!  type *,' kick_w, kick_e ', kick_w + dv, kick_e +dv
!  type *,' Final separation: ', delta_ip%vec(1)

! Reset initial position
     delta_ip_0%vec = delta_ip%vec
  end do

  final_pos_out%vec(1:4)=delta_ip_0%vec(1:4)
  
  if(n>= 10)call out_io(s_warn$,r_name,' VERT NOT CLOSED')
  
!  if(n>= 10)type *,' VERT NOT CLOSED'

  return
  
end subroutine close_vertical
                                                            










