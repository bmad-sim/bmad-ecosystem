!+
! Subroutine quad_beta_ave (lat, ix_ele, beta_a_ave, beta_b_ave)
!
! Subroutine to compute the average betas in a quad
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lat structure
!   ix_ele -- Integer: Index of quadrupole
!
! Output:
!   beta_a_ave, beta_b_ave -- Real(rp): Average betas in the quad.
!
! NOTE: This subroutine is only valid if there is no local coupling
!-

#include "CESR_platform.inc"

subroutine quad_beta_ave (lat, ix_ele, beta_a_ave, beta_b_ave)

  use bmad_struct
  use bmad_interface, except_dummy => quad_beta_ave

  implicit none

  type (lat_struct)  lat
  type (ele_struct), save ::  ele

  integer ix_ele, ix

  real(rp) beta_a_ave, beta_b_ave, k_quad

! Since the beta stored in the ELE is at the end we need
! to invert alpha for BEAVE routine

  ele = lat%ele(ix_ele)

  if (ele%key /= quadrupole$ .and. ele%key /= sol_quad$ .and. &
        ele%key /= wiggler$) then
    print *, 'ERROR IN QUAD_BETA_AVE: ELEMENT NOT A QUAD, SOL_QUAD OR WIGGLER'
    print *, '      ', ele%name, '  ', key_name(ele%key)
    call err_exit
  endif

! if a controller element just use the value at the end of the first
! controlled element

  if (ix_ele > lat%n_ele_track) then
    ix = lat%control(ele%ix1_slave)%ix_slave
    beta_a_ave = lat%ele(ix)%a%beta
    beta_b_ave = lat%ele(ix)%b%beta
    return
  endif

! otherwise proceed as normal

  if (ele%a%beta == 0 .or. ele%b%beta == 0) then
    print *, 'ERROR IN QUAD_BETA_AVE: BETA IS ZERO AT END:',  &
                                             ele%a%beta, ele%b%beta
    print *, ele%name, '  ', key_name(ele%key)
    call err_exit
  endif

  k_quad = ele%value(k1$)
  
  if( ele%key == wiggler$)then
    beta_a_ave = b_ave(ele%a%beta, -ele%a%alpha,  0.0_rp, ele%value(l$))
    beta_b_ave = b_ave(ele%b%beta, -ele%b%alpha, -k_quad, ele%value(l$))
  else
    beta_a_ave = b_ave(ele%a%beta, -ele%a%alpha,  k_quad, ele%value(l$))
    beta_b_ave = b_ave(ele%b%beta, -ele%b%alpha, -k_quad, ele%value(l$))
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

  real(rp) beta,alpha,k,l,kk,a,x,g,gamma, this_ave
  
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
