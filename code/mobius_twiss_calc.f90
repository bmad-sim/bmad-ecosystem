!+
! Subroutine MOBIUS_TWISS_CALC (ELE, V_MAT)
!
! Subroutine calculate the mobius betas and etas which are effective
! projections of beta and eta in the X and Y planes. This is used for highly
! coupled lattices. For uncoupled lattices the mobius betas and etas
! are the same as the normal mode betas and etas.
!
! Modules Needed:
!   use bmad
!
! Input:
!     ELE        -- Ele_struct: Element
!     V_MAT(4,4) -- Real(rp): Normal mode to X-Y coords transformation
!
! Output:
!     ELE.X.MOBIUS_BETA -- Mobius betas
!     ELE.Y.MOBIUS_BETA
!     ELE.X.MOBIUS_ETA  -- Mobius etas
!     ELE.Y.MOBIUS_ETA
!
!-

!$Id$
!$Log$
!Revision 1.5  2003/07/09 01:38:16  dcs
!new bmad with allocatable ring%ele_(:)
!
!Revision 1.4  2003/01/27 14:40:39  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:19  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine mobius_twiss_calc (ele, v_mat)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct)  ele
  type (twiss_struct)  a, b

  real(rp) c11, c12, c21, c22, vec(4)
  real(rp) v_mat(4,4)

!

  a = ele%x
  b = ele%y

! beta calc

  c11 = ele%c_mat(1,1)
  c12 = ele%c_mat(1,2)
  c21 = ele%c_mat(2,1)
  c22 = ele%c_mat(2,2)

  ele%x%mobius_beta = ele%gamma_c**2 * a%beta + b%beta * c11**2 -  &
                              2 * b%alpha * c11 * c12  + b%gamma * c12**2

  ele%y%mobius_beta = ele%gamma_c**2 * b%beta + a%beta * c22**2 +  &
                              2 * a%alpha * c12 * c22  + a%gamma * c12**2

! eta calc

  vec(1) = a%eta
  vec(2) = a%etap
  vec(3) = b%eta
  vec(4) = b%etap

  vec = matmul (v_mat, vec)

  ele%x%mobius_eta = vec(1)
  ele%y%mobius_eta = vec(3)

  return
  end
