module mode3_mod

use bmad_struct
use bmad_interface
use eigen_mod

contains

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine normal_mode3_calc (mat, tune, G, V, synchrotron_motion)
!
! Decompose a 2n x 2n sympectic matrix into normal modes.
! For more details see:
!   ......
!
! Input: 
!   mat(2n:2n)  -- Real(rp): Symplectic matrix to decompose
!   syncrhotron_motion 
!               -- Logical, optional: if present and true then it is assumed that tune(n) 
!                    corresponds to the synchrotron motion and that tune(n) < 0
! Output: 
!   tune(n)  -- Real(rp): Eigenvalues are the tunes. 
!                 The ambiguity in tune is resolved so that det(G) = 1
!   G(2n:2n) -- Real(rp):  G is block diagonal. 
!                 Each 2X2 block is the matrix G for that mode 
!                 G(1,1) = 1/sqrt{beta}, G(2,2) = sqrt{beta}, etc.
!   V(2n:2n) -- Real(rp):  VUV^{-1} =  mat, where U is block diagonal and
!                 The on-diagonal V blocks are proportional to the unit matrix.
!-

subroutine normal_mode3_calc (mat, tune, G, V, synchrotron_motion)

implicit none

real(rp) mat(:,:), V(:,:), G(:,:), tune(:)

real(rp), allocatable :: eval_r(:), eval_i(:), evec_r(:,:), evec_i(:,:)
real(rp), allocatable ::  ReZ(:,:), ImZ(:,:)
real(rp), allocatable ::  ReZt(:,:), ImZt(:,:), Z(:,:)
real(rp), allocatable ::  ReKInv(:,:), ImKInv(:,:)
real(rp), allocatable :: Jtot(:,:)
real(rp), allocatable :: sum(:)
real(rp) Jzero(2,2), g2_mat(2,2), r_mat(2,2), h_mat(2,2)
real(rp) sdetZi, oroot2
real(rp) phi

integer n, i, j, l, p1, p2

logical error
logical, optional :: synchrotron_motion

! Error check

 n = size(mat, dim=1)
 p2 = size(mat, dim=2)
 if (p2 /= n .or. 2*(n/2) /= n) then
   print *,' The matrix must be square with even number of rows '
   stop
 endif

! Construct matrix for transforming to real basis

allocate(ReKInv(1:n,1:n), ImKInv(1:n,1:n))

oroot2 = 1/sqrt(2.)

ReKInv = 0
ImKInv = 0

do i = 2, n, 2
  ReKInv(i-1,i-1) = oroot2
  ReKInv(i,i-1) = oroot2
  ImKInv(i-1,i)= -oroot2
  ImKInv(i,i) = oroot2
end do

allocate(eval_r(n), eval_i(n), evec_r(n,n), evec_i(n,n))
call mat_eigen (mat,eval_r, eval_i, evec_r, evec_i, error)

do i = 2, n, 2
  tune(i/2) = atan2(eval_i(i-1),eval_r(i-1))
end do
  
! construct Z

allocate(ReZ(1:n,1:n),ImZ(1:n,1:n))
ReZ=0
ImZ=0

do i = 1,n
 do j=1,n
  ReZ(i,j) = evec_r(j,i)
  ImZ(i,j) = evec_i(j,i)
 end do
end do

! transform to real basis

allocate(ReZt(1:n,1:n), ImZt(1:n,1:n))
ReZt=0
ImZt=0

ReZt = matmul(ReZ,ReKInv) - matmul(ImZ,ImKInv)
ImZt = matmul(ImZ,ReKInv) + matmul(ReZ,ImKInv)

allocate(Z(1:n,1:n))

!normalize so that Z is symplectic and determinant of diagonal block is > 0

allocate(sum(n))
do i = 1, n/2

  sum(i) = 0.

  do j=1,n/2
    sum(i) = sum(i) + determinant(ReZt(2*j-1:2*j,2*i-1:2*i))
  end do

  if (determinant(ReZt(2*i-1:2*i,2*i-1:2*i)) > 0) then
    Z(1:n,2*i-1) = ReZt(1:n,2*i-1) /sqrt(abs(sum(i)))
    Z(1:n,2*i) = ReZt(1:n,2*i) /sqrt(abs(sum(i)))
   else !swap the columns
    Z(1:n,2*i) = ReZt(1:n,2*i-1) /sqrt(abs(sum(i)))
    Z(1:n,2*i-1) = ReZt(1:n,2*i) /sqrt(abs(sum(i)))
    tune(i) = twopi-tune(i)
   endif

end do

if (logic_option(.false., synchrotron_motion)) tune(n/2) = twopi-tune(n/2)

! Find G and combine into n X n matrix G
! create matrix J for transformation of H to G

allocate(Jtot(1:n,1:n))
Jtot = 0.
G = 0

do i = 2, n, 2
  l = i-1
  j = i
  sdetZi = abs(Z(l,l)*Z(j,j)-Z(l,j)*Z(j,l))
  H_mat = Z(l:j,l:j)/sqrt(sdetZi)

  phi = -atan2(h_mat(1,2),h_mat(1,1))
  if(h_mat(1,1)/cos(phi) < 0)phi = phi + twopi/2
  r_mat = cos(phi)
  r_mat(1,2)= sin(phi)
  r_mat(2,1) = -sin(phi)
  g2_mat = matmul(h_mat, r_mat)

  Jtot(i-1:i,i-1:i) = r_mat

  G(i-1,i-1) = g2_mat(2,2)
  G(i,i) = g2_mat(1,1)
  G(i,i-1) = -g2_mat(2,1)
end do

V = matmul(matmul(Z, Jtot), G)

deallocate(ReZ, ReZt, ImZt, Z)
deallocate(ReKInv, ImKInv)
deallocate(sum, Jtot)
deallocate(eval_r, eval_i, evec_r, evec_i)

end subroutine normal_mode3_calc

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine twiss3_propagate_all (lat)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
!-

subroutine twiss3_propagate_all (lat)

implicit none

type (lat_struct) lat

integer i

do i = 1, lat%n_ele_track
  call twiss3_propagate1 (lat%ele(i-1), lat%ele(i))
enddo

end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine twiss3_propagate1 (ele1, ele2)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
!-

subroutine twiss3_propagate1 (ele1, ele2)

implicit none

type mat2_struct
  real(rp) m(2,2)
end type

type (ele_struct) ele1, ele2
type (mat2_struct) w(3)

real(rp) gamma(3), tv(6,6), w_inv(2,2)

integer i, ik

!

if (.not. associated(ele2%mode3)) allocate(ele2%mode3)

tv = matmul (ele2%mat6, ele1%mode3%v)

do i = 1, 3
  ik = 2 * i - 1
  w(i)%m = tv(ik:ik+1,ik:ik+1)
  gamma(i) = sqrt(determinant (w(i)%m))
  w(i)%m = w(i)%m / gamma(i)
  call mat_symp_conj (w(i)%m, w_inv)
  ele2%mode3%v(1:6, ik:ik+1) = matmul(tv(1:6, ik:ik+1), w_inv)
enddo

call twiss1_propagate (ele1%mode3%a, w(1)%m,  ele2%value(l$), ele2%mode3%a)
call twiss1_propagate (ele1%mode3%b, w(2)%m,  ele2%value(l$), ele2%mode3%b)
call twiss1_propagate (ele1%mode3%c, w(3)%m,  0.0_rp,         ele2%mode3%c)

end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine twiss3_at_start (lat)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
!-

subroutine twiss3_at_start (lat)

implicit none

type (lat_struct) lat
real(rp) g(6,6), tune3(3)
integer n

!

if (.not. associated(lat%ele(0)%mode3)) allocate(lat%ele(0)%mode3)

call transfer_matrix_calc (lat, .true., lat%param%t1_with_RF)
call normal_mode3_calc (lat%param%t1_with_RF, tune3, g, lat%ele(0)%mode3%v)

call mode1_calc (g(1:2, 1:2), tune3(1), lat%ele(0)%mode3%a)
call mode1_calc (g(3:4, 3:4), tune3(2), lat%ele(0)%mode3%b)
call mode1_calc (g(5:6, 5:6), tune3(3), lat%ele(0)%mode3%c)

n = lat%n_ele_track
lat%ele(n)%mode3%a%phi = tune3(1)
lat%ele(n)%mode3%b%phi = tune3(2)
lat%ele(n)%mode3%c%phi = tune3(3)

!-------------------------------------------------------------------------------------
contains

subroutine mode1_calc (g, tune, twiss)

type (twiss_struct) twiss
real(rp) g(2,2), tune

!

twiss%beta = g(2,2)**2
twiss%alpha = g(2,1) * g(2,2)
twiss%gamma = (1 + twiss%alpha**2) / twiss%beta
twiss%phi = 0

end subroutine

end subroutine

end module
