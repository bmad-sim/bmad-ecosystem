module four_mod

contains

subroutine xy_to_action(ring, ix, X, J, ok)

  use bmad
  use sim_utils_interface

  type(lat_struct) ring
  integer ix, i
  real(rp) X(6)
  real(rp) J(6)
  logical ok

  real(rp) Afour(6,6)
  real(rp) Afour_inv(6,6)
  type(ele_struct), pointer :: ele

  real(rp) sqrtabeta
  real(rp) sqrtbbeta

  ele => ring%ele(ix)

  sqrtabeta = sqrt(ele%a%beta)
  sqrtbbeta = sqrt(ele%b%beta)

  Afour = 0.0d0
  Afour(1,1) = sqrtabeta
  Afour(2,1) = -ele%a%alpha/sqrtabeta
  Afour(2,2) = 1.0d0/sqrtabeta
  Afour(1,6) = ele%x%eta
  Afour(2,6) = ele%x%etap
  Afour(3,3) = sqrtbbeta
  Afour(4,3) = -ele%b%alpha/sqrtbbeta
  Afour(4,4) = 1.0d0/sqrtbbeta
  Afour(3,6) = ele%y%eta
  Afour(4,6) = ele%y%etap
  Afour(5,5) = 1.0d0
  Afour(6,6) = 1.0d0

  ! do i=1,6
  !   write(*,'(6ES14.4)') Afour(i,:)
  ! enddo

  !call mat_inverse(Afour,Afour_inv, ok = ok)
  Afour_inv = mat_symp_conj(Afour)
  ok = .true.

  J = matmul(Afour_inv,X)

end subroutine

end module
