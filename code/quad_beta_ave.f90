!+
! Subroutine quad_beta_ave (ring, ix_ele, beta_x_ave, beta_y_ave)
!
! Subroutine to compute the average betas in a quad
!
! Modules Needed:
!   use bmad
!
! Input:
!     ring   -- Ring_struct: Ring structure
!     ix_ele -- Integer: Index of quadrupole
!
! Output:
!     beta_x_ave, beta_y_ave -- Real(rdef): Average betas in the quad.
!
! NOTE: This subroutine is only valid if there is no local coupling
!-

!$Id$
!$Log$
!Revision 1.6  2003/05/02 15:44:01  dcs
!F90 standard conforming changes.
!
!Revision 1.5  2003/01/27 14:40:41  dcs
!bmad_version = 56
!
!Revision 1.4  2002/06/13 14:54:28  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:22  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine quad_beta_ave (ring, ix_ele, beta_x_ave, beta_y_ave)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele

  integer ix_ele, ix

  real(rdef) beta_x_ave, beta_y_ave, k_quad

! Since the beta stored in the ELE is at the end we need
! to invert alpha for BEAVE routine

  ele = ring%ele_(ix_ele)

  if (ele%key /= quadrupole$ .and. ele%key /= sol_quad$ .and. &
        ele%key /= wiggler$) then
    print *, 'ERROR IN QUAD_BETA_AVE: ELEMENT NOT A QUAD, SOL_QUAD OR WIGGLER'
    print *, '      ', ele%name, '  ', key_name(ele%key)
    call err_exit
  endif

! if a controller element just use the value at the end of the first
! controlled element

  if (ix_ele > ring%n_ele_ring) then
    ix = ring%control_(ele%ix1_slave)%ix_slave
    beta_x_ave = ring%ele_(ix)%x%beta
    beta_y_ave = ring%ele_(ix)%y%beta
    return
  endif

! otherwise proceed as normal

  if (ele%x%beta == 0 .or. ele%y%beta == 0) then
    print *, 'ERROR IN QUAD_BETA_AVE: BETA IS ZERO AT END:',  &
                                             ele%x%beta, ele%y%beta
    print *, ele%name, '  ', key_name(ele%key)
    call err_exit
  endif

  k_quad = ele%value(k1$)
  
  if( ele%key == wiggler$)then
    beta_x_ave = b_ave(ele%x%beta, -ele%x%alpha,  0.0_rdef, ele%value(l$))
    beta_y_ave = b_ave(ele%y%beta, -ele%y%alpha, -k_quad, ele%value(l$))
  else
    beta_x_ave = b_ave(ele%x%beta, -ele%x%alpha,  k_quad, ele%value(l$))
    beta_y_ave = b_ave(ele%y%beta, -ele%y%alpha, -k_quad, ele%value(l$))
  endif

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

!+
!  real function b_ave(beta, alpha, kk, l)
!
!  returns average beta in  quad or drift given twiss parameters
!  beta and alpha at quad entrance, and quad strength and length k,l.
!-

function b_ave (beta, alpha, kk, l) result (this_ave)

  use precision_def

  implicit none

  real(rdef) beta,alpha,k,l,kk,a,x,g,gamma, this_ave
  
  k=sqrt(abs(kk))

  if (beta==0.0) then
    print *, 'ERROR IN B_AVE: INITIAL BETA IS ZERO'
    this_ave = -1
  endif
  
  gamma=(1.0+alpha**2)/beta

  if ((kk*l)==0.0) then
    this_ave=beta-alpha*l+gamma*l*l/3.  ! drift
    return
  endif

!

  a=alpha/k
  g=gamma/(k*k)
  x=k*l

  if(kk > 0.0) then ! focus case
    this_ave = 0.5*((beta+g)+(beta-g)*sin(2.*x)/(2.*x))
    this_ave = this_ave-a*(1.-cos(2.*x))/(2.*x)
  else
    this_ave = 0.5*((beta-g)+(beta+g)*sinh(2.*x)/(2.*x))
    this_ave = this_ave-a*(cosh(2.*x)-1.)/(2.*x)
  endif
  
end function

end subroutine
