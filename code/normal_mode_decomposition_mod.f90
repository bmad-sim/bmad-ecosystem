module normal_mode_decomposition_mod
     use bmad
     use eigen_mod

      type g_struct
        real(rp) mat(2,2)
      end type

      type decomposition_struct
        real(rp) mat_to_ele(6,6), mat(6,6), W(6,6), Y(6,6), tune(3)
      end type


!-----------------------
contains




   subroutine normal_mode_decomposition(mat,tune,Y, W, synchrotron_motion)
!
!  Decompose 2n X 2n sympectic matrix into normal modes.
!         First find eigenvalues and eigenvectors
!         Construct similarity transform V (ReV and ImV since V is complex) from diagonal basis to lab basis.
!         Transform ReV and ImV to V in real basis
!         Choose scale factor for each pair of columns of V so that V is symplectic 
!               (sum of determinant of 2X2 sub matrices is 1)
!         Choose relative phase of each pair of columns so that G(1,2)=0, and |G|=1
!
!     INPUT: 
!       Real -- mat(2n:2n): symplectic matrix to decompose
!       Logical (optional) -- syncrhotron_motion: if true then it is assumed that tune(n) corresponds to 
!                                 the synchrotron motion and that tune(n)<0
!     OUTPUT: 
!       Real -- tune(n):   eigenvalues are the tunes. The ambiguity in tune is resolved so that det(G)=1
!       Real -- Y(2n:2n):  Y is block diagonal. Each 2X2 block is the matrix G for that mode 
!                           G(1,1)=sqrt{beta}, G(1,2)=0, G(2,1)=alpha/sqrt{beta}, G(2,2)=1/sqrt{beta}
!       Real -- W(2n:2n):  WUW^{-1}= mat, where U is block diagonal with 2X2 blocks A_1,A_2,A_3,...
!                           A_1(1,1)=cos(tune(1))+alpha*sin(tune(1)), A_2(1,2)=beta*sin(tune(1)), etc. 

     implicit none

      type (g_struct) H(3),JJJ(3),G(3)

      real(rp) mat(:,:), W(:,:), Y(:,:), tune(:)
      real(rp), allocatable :: eval_r(:), eval_i(:), evec_r(:,:), evec_i(:,:)

      real(rp), allocatable ::  ReV(:,:), ImV(:,:), ReK(:,:), ImK(:,:)
      real(rp), allocatable ::  ReVt(:,:), ImVt(:,:), V(:,:)
      real(rp), allocatable ::  ReKInv(:,:), ImKInv(:,:)
      real(rp) sdetVi, oroot2
      real(rp), allocatable :: sum(:)
      real(rp) Jzero(2,2)
      real(rp) phi
      real(rp), allocatable :: YInv(:,:), Jtot(:,:)
      real(rp) deterrminant
      real(rp) M(2,2)

      integer n,i,j, l,p1, p2

      logical error
      logical, optional :: synchrotron_motion

       n = size(mat, dim=1)
       p2 = size(mat, dim=2)
       if(p2/= n .or. 2*(n/2)/= n)then
         print *,' The matrix must be square with even number of rows '
         stop
       endif


! construct matrix for transforming to real basis
  
       allocate(ReK(1:n,1:n), ImK(1:n,1:n))
       allocate(ReKInv(1:n,1:n), ImKInv(1:n,1:n))
       oroot2 = 1/sqrt(2.)
       ReK(1:n,1:n) = 0
       ImK(1:n,1:n) = 0
       ReKInv(1:n,1:n) = 0
       ImKInv(1:n,1:n) = 0
       do i =2,n,2
         ReK(i-1,i-1) = oroot2
         ReK(i-1,i) = oroot2
         ImK(i,i-1)= oroot2
         ImK(i,i) = -oroot2

         ReKInv(i-1,i-1) = oroot2
         ReKInv(i,i-1) = oroot2
         ImKInv(i-1,i)= -oroot2
         ImKInv(i,i) = oroot2
       end do

       allocate(eval_r(n), eval_i(n), evec_r(n,n), evec_i(n,n))
       call mat_eigen (mat,eval_r, eval_i, evec_r, evec_i, error)


!       do i=1,n
!        print '(3e12.4)', eval_r(i), eval_i(i), eval_r(i)**2+eval_i(i)**2
!       end do

       do i = 2,n,2
        tune(i/2) = atan2(eval_i(i-1),eval_r(i-1))
!        print '(e12.4)', tune(i/2)/twopi
       end do
       
  ! construct V

    allocate(ReV(1:n,1:n),ImV(1:n,1:n))
    ReV(1:n,1:n)=0
    ImV(1:n,1:n)=0

    do i = 1,n
      do j=1,n
       ReV(i,j) = evec_r(j,i)
       ImV(i,j) = evec_i(j,i)
      end do
    end do


! transform to real basis
    allocate(ReVt(1:n,1:n),ImVt(1:n,1:n))
    ReVt(1:n,1:n)=0
    ImVt(1:n,1:n)=0

    ReVt = matmul(ReV,ReKInv) - matmul(ImV,ImKInv)
    ImVt = matmul(ImV,ReKInv) + matmul(ReV,ImKInv)

  allocate(V(1:n,1:n))



  !normalize so that V is symplectic and determinant of diagonal block is > 0
 
     allocate(sum(n))
     do i=1,n/2
       sum(i) = 0.
        do j=1,n/2
          sum(i) = sum(i) + deterrminant(ReVt(2*j-1:2*j,2*i-1:2*i))
        end do

        if (deterrminant(ReVt(2*i-1:2*i,2*i-1:2*i)) > 0) then
          V(1:n,2*i-1) = ReVt(1:n,2*i-1) /sqrt(abs(sum(i)))
          V(1:n,2*i) = ReVt(1:n,2*i) /sqrt(abs(sum(i)))
         else !swap the columns
          V(1:n,2*i) = ReVt(1:n,2*i-1) /sqrt(abs(sum(i)))
          V(1:n,2*i-1) = ReVt(1:n,2*i) /sqrt(abs(sum(i)))
          tune(i) = twopi-tune(i)

        endif

     end do

     if(present(synchrotron_motion))then
       if(synchrotron_motion)tune(n/2)=twopi-tune(n/2)
     endif

  ! Find G and combine into n X n matrix Y

       Y(1:n, 1:n) = 0
       allocate(YInv(1:n,1:n))
       YInv(1:n,1:n) = 0

     do i=2,n,2
      l=i-1
      j=i
      sdetVi = abs(V(l,l)*V(j,j)-V(l,j)*V(j,l))
      H(i/2)%mat(1:2,1:2) = V(l:j,l:j)/sqrt(sdetVi)
     end do

  ! create matrix J for transformation of H to G
       allocate(Jtot(1:n,1:n))
       Jtot(1:n,1:n)=0.

    do i =2,n,2
     phi = -atan2(H(i/2)%mat(1,2),H(i/2)%mat(1,1))
     if(H(i/2)%mat(1,1)/cos(phi) < 0)phi = phi + twopi/2
     JJJ(i/2)%mat(1:2,1:2)=cos(phi)
     JJJ(i/2)%mat(1,2)= sin(phi)
     JJJ(i/2)%mat(2,1) = -sin(phi)
     G(i/2)%mat(1:2,1:2) = matmul(H(i/2)%mat(1:2,1:2),JJJ(i/2)%mat(1:2,1:2))

     Jtot(i-1:i,i-1:i)=JJJ(i/2)%mat(1:2,1:2)

     Y(i-1:i,i-1:i) = G(i/2)%mat(1:2,1:2)
     YInv(i-1,i-1) = G(i/2)%mat(2,2)
     YInv(i,i) = G(i/2)%mat(1,1)
     YInv(i,i-1) = -G(i/2)%mat(2,1)
    end do

     W=matmul(matmul(V(1:n,1:n),Jtot(1:n,1:n)),YInv(1:n,1:n))

      deallocate(ReV, ImV, ReK, ImK, ReVt, ImVt, V)
      deallocate( ReKInv, ImKInv)
      deallocate(sum)
      deallocate (YInv, Jtot)
      deallocate(eval_r, eval_i, evec_r, evec_i)

     return
     end subroutine normal_mode_decomposition

   function deterrminant(M)

    implicit none

     real(rp) M(2,2)
     real(rp) x, deterrminant

     x = M(1,1)*M(2,2)-M(1,2)*M(2,1)
     
     deterrminant = x

  return

  end function deterrminant

end module





