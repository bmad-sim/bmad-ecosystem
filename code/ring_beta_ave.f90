!+
! Subroutine Ring_Beta_Ave(ring, cesr)
!
! Calculate the vertical and horizontal average betas for quads around the ring
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     ring -- record/ring_struct/
!
! Output:
!     cesr -- record/cesr_struct/
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine Ring_Beta_Ave(ring, cesr)
  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (cesr_struct)  cesr
  type (ele_struct)  elem
  integer rindex, qindex
  real betaxAve, betayAve

  do rindex = 1, ring%n_ele_symm
   elem = ring%ele_(rindex)
   if(elem%key==quadrupole$ .or. elem%key==sol_quad$) then
    if(elem%name(1:2)/='SK') then !since bmad doesn't know any better
      qindex = cesr%ix_cesr(rindex)
      call Beta_Ave(ring, rindex, betaxAve, betayAve)
      cesr%quad_(qindex)%x%beta_ave = betaxAve
      cesr%quad_(qindex)%y%beta_ave = betayAve
    endif
   endif
  enddo

  return
  end
