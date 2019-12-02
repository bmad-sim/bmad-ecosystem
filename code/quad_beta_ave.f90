!+
! Subroutine quad_beta_ave (ele, beta_a_ave, beta_b_ave)
!
! Subroutine to compute the average betas in a quad
!
! Input:
!   ele   -- ele_struct: Element.
!
! Output:
!   beta_a_ave, beta_b_ave -- Real(rp): Average betas in the quad.
!
! NOTE: This subroutine is only valid if there is no local coupling
!-

subroutine quad_beta_ave (ele, beta_a_ave, beta_b_ave)

use bmad_interface, except_dummy => quad_beta_ave

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
type (branch_struct), pointer :: branch

integer ix

real(rp) beta_a_ave, beta_b_ave, k_quad

character(16), parameter :: r_name = 'quad_beta_ave'

! Since the beta stored in the ELE is at the end we need
! to invert alpha for BEAVE routine

if (ele%key /= quadrupole$ .and. ele%key /= sol_quad$ .and. &
      ele%key /= wiggler$ .and. ele%key /= undulator$) then
  call out_io (s_fatal$, r_name, 'ELEMENT NOT A QUAD, SOL_QUAD OR WIGGLER: [\i0\] ' // trim(ele%name) // '   ' // key_name(ele%key), ele%ix_ele)
  if (global_com%exit_on_error) call err_exit
endif

! if a controller element just use the value at the end of the first
! controlled element

branch => ele%branch
if (ele%ix_ele > branch%n_ele_track) then
  slave => pointer_to_slave(ele, 1)
  beta_a_ave = slave%a%beta
  beta_b_ave = slave%b%beta
  return
endif

! otherwise proceed as normal

if (ele%a%beta == 0 .or. ele%b%beta == 0) then
  call out_io (s_fatal$, r_name, 'BETA IS ZERO AT END:' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

k_quad = ele%value(k1$)

if (ele%key == wiggler$ .or. ele%key == undulator$)then
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
! function b_ave(beta, alpha, kk, l) result (this_ave)
!
! returns average beta in  quad or drift given twiss parameters
! beta and alpha at quad entrance, and quad strength and length k,l.
!-

function b_ave (beta, alpha, kk, l) result (this_ave)

use precision_def

implicit none

real(rp) beta,alpha,k,l,kk,a,x,g,gamma, this_ave

!

k=sqrt(abs(kk))

if (beta==0.0) then
  call out_io (s_error$, r_name, 'INITIAL BETA IS ZERO')
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
