!+
! Subroutine QUAD_BETA_AVE (RING, IX_ELE, BETA_X_AVE, BETA_Y_AVE)
!
! Subroutine to compute the average betas in a quad
!
! Modules Needed:
!   use bmad
!
! Input:
!     RING   -- Ring_struct: Ring structure
!     IX_ELE -- Integer: Index of quadrupole
!
! Output:
!     BETA_X_AVE, BETA_Y_AVE -- Real(rdef): Average betas in the quad.
!
! NOTE: This subroutine is only valid if there is no local coupling
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:22  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine quad_beta_ave (ring, ix_ele, beta_x_ave, beta_y_ave)

  use bmad
  implicit none
  type (ring_struct)  ring
  type (ele_struct)  ele

  integer ix_ele, ix

  real(rdef) beta_x_ave, beta_y_ave, beave, k_quad

! Since the beta stored in the ELE is at the end we need
! to invert alpha for BEAVE routine

  ele = ring%ele_(ix_ele)

  if (ele%key /= quadrupole$ .and. ele%key /= sol_quad$ .and. &
        ele%key /= wiggler$) then
    type *, 'ERROR IN QUAD_BETA_AVE: ELEMENT NOT A QUAD, SOL_QUAD OR WIGGLER'
    type *, '      ', ele%name, '  ', key_name(ele%key)
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
    type *, 'ERROR IN QUAD_BETA_AVE: BETA IS ZERO AT END:',  &
                                             ele%x%beta, ele%y%beta
    type *, ele%name, '  ', key_name(ele%key)
    call err_exit
  endif

  k_quad = ele%value(k1$)
  if( ele%key == wiggler$)then
    beta_x_ave = beave(ele%x%beta, -ele%x%alpha,  0., ele%value(l$))
    beta_y_ave = beave(ele%y%beta, -ele%y%alpha, -k_quad, ele%value(l$))
  else
    beta_x_ave = beave(ele%x%beta, -ele%x%alpha,  k_quad, ele%value(l$))
    beta_y_ave = beave(ele%y%beta, -ele%y%alpha, -k_quad, ele%value(l$))
  endif
  return
  end

