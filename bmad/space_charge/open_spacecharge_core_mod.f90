! OpenSC: a package for computing 3D space charge fields using convolution-based methods with integrated Green functions (IGFs).
! Boundary conditions in this version: (1) free-space, (2) rectangular pipe
! R.D. Ryne, February 2018

module open_spacecharge_core_mod

use, intrinsic :: iso_fortran_env
!uncomment the following for the version that uses the 1D FFT from Alan Miller's web page
!and compile with: gfortran test_mod.f90 fast_fourier_am.f90 open_spacecharge_mod.f90 test_opensc.f90
!use fast_fourier_am
!uncomment the following for the version that uses FFTW:
use fft_interface_mod
!$ use omp_lib

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128
!the following arrays are Green function arrays used by the rectangular pipe solver
complex(dp), allocatable, private,  dimension(:,:,:) :: cgrn1,cgrn2,cgrn3,cgrn4
integer, private :: g1ilo,g1jlo,g1klo,g2ilo,g2jlo,g2klo,g3ilo,g3jlo,g3klo,g4ilo,g4jlo,g4klo !lower bounds of grn1234 arrays

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine osc_freespace_solver(rho,gam0,delta,phi,efield,bfield,nlo,nhi,nlo_gbl,nhi_gbl,npad,idirectfieldcalc,igfflag)
! This will check the allocation of the global cgrn1. 
! OLD: note well: this assumes that osc_alloc_freespace_array has been called so that cgrn1 is already allocated!
! 
!use mpi
implicit none
integer, intent(in) :: igfflag,idirectfieldcalc
real(dp), intent(in), dimension(3) :: delta
integer, intent(in), dimension(3) :: nlo,nhi,nlo_gbl,nhi_gbl,npad
real(dp), intent(in) :: gam0
real(dp), intent(in), dimension(:,:,:) :: rho
real(dp), intent(out), dimension(:,:,:) :: phi
real(dp), intent(out), dimension(:,:,:,:) :: efield,bfield
real(dp) :: dx,dy,dz
integer :: ilo,ihi,jlo,jhi,klo,khi
integer :: ilo2,ihi2,jlo2,jhi2,klo2,khi2,iperiod,jperiod,kperiod
complex(dp), allocatable, dimension(:,:,:) :: crho

integer :: icomp,i,j,k,im1,ip1,jm1,jp1,km1,kp1
real(dp) :: beta0,xfac,yfac,zfac
integer :: mprocs,myrank,ierr
real(dp), parameter :: clight=299792458.d0
!
!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
myrank=0


if (igfflag == 0 .and. idirectfieldcalc.eq.1) then
  print *, 'Error: Direct field calc must use integrated Green function'
  print *, 'Aborting...'
  return
endif

dx=delta(1); dy=delta(2); dz=delta(3)
ilo=nlo(1);ihi=nhi(1);jlo=nlo(2);jhi=nhi(2);klo=nlo(3);khi=nhi(3)
ilo2=ilo; jlo2=jlo; klo2=klo

! Always call this. It will check the cgrn1 size, and re-allocate if things changed. 
call osc_alloc_freespace_array(nlo, nhi, npad)

iperiod=size(cgrn1,1); jperiod=size(cgrn1,2); kperiod=size(cgrn1,3)
ihi2=ilo2+iperiod-1; jhi2=jlo2+jperiod-1; khi2=klo2+kperiod-1
allocate(crho(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size complex array for the charge density
!! write(6,*)'iperiod,jperiod,kperiod=',iperiod,jperiod,kperiod

if(idirectfieldcalc.eq.0)then
  if(myrank.eq.0)write(6,*)'Solving for Phi'
  icomp=0
!compute/store the FFT of the charge density:
call getrhotilde(rho,crho,ilo,jlo,klo)
!compute/store the FFT of the green function:
call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
!compute the convolution:   (note that crho,cgrn1 are complex padded, phi is real not padded)
call conv3d(crho,cgrn1,phi,ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)

!original version did not handle the edge values correctly when differencing; here is a first order approx:
  do k=lbound(phi,3),ubound(phi,3)
    km1=k-1; if(km1.lt.lbound(phi,3))km1=k
    kp1=k+1; if(kp1.gt.ubound(phi,3))kp1=k
    zfac=merge(2.d0,1.d0,(km1.eq.k).or.(kp1.eq.k))
    do j=lbound(phi,2),ubound(phi,2)
      jm1=j-1; if(jm1.lt.lbound(phi,2))jm1=j
      jp1=j+1; if(jp1.gt.ubound(phi,2))jp1=j
      yfac=merge(2.d0,1.d0,(jm1.eq.j).or.(jp1.eq.j))
      do i=lbound(phi,1),ubound(phi,1)
        im1=i-1; if(im1.lt.lbound(phi,1))im1=i
        ip1=i+1; if(ip1.gt.ubound(phi,1))ip1=i
        xfac=merge(2.d0,1.d0,(im1.eq.i).or.(ip1.eq.i))
          efield(i,j,k,1)=-(phi(ip1,j,k)-phi(im1,j,k))/(2.d0*dx)*gam0*xfac
          efield(i,j,k,2)=-(phi(i,jp1,k)-phi(i,jm1,k))/(2.d0*dy)*gam0*yfac
          efield(i,j,k,3)=-(phi(i,j,kp1)-phi(i,j,km1))/(2.d0*dz)/gam0*zfac
      enddo
    enddo
  enddo
endif

if(idirectfieldcalc.eq.1)then
  if(myrank.eq.0)write(6,*)'Solving for Ex, Ey, Ez on mesh:', nhi
  if(.not.allocated(crho))allocate(crho(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
  call getrhotilde(rho,crho,ilo,jlo,klo)
  if(.not.allocated(cgrn1))then
    write(6,*)'something wrong. freespace green function should have been allocated during initialization'
    stop
  endif
  icomp=1
  call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
  call conv3d(crho,cgrn1,efield(:,:,:,1),ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  icomp=2
  call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
  call conv3d(crho,cgrn1,efield(:,:,:,2),ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  icomp=3
  call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
  call conv3d(crho,cgrn1,efield(:,:,:,3),ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
endif

! set the magnetic field:
beta0=sqrt(1-1/gam0**2)
do k=lbound(phi,3),ubound(phi,3)
  do j=lbound(phi,2),ubound(phi,2)
    do i=lbound(phi,1),ubound(phi,1)
      bfield(i,j,k,1)=-efield(i,j,k,2)*beta0/clight
      bfield(i,j,k,2)= efield(i,j,k,1)*beta0/clight
      bfield(i,j,k,3)=0.d0
    enddo
  enddo
enddo

end subroutine osc_freespace_solver

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!+
!
! Allocates the global cgrn1 double-sized mesh
!
! If it was already allocated, it will be resized
!
!
subroutine osc_alloc_freespace_array(nlo,nhi,npad)
implicit none
integer, intent(in), dimension(3) :: nlo, nhi, npad
integer :: i, glo(3), ghi(3)
logical :: good_sizes

! Set bounds
glo = nlo - nhi
ghi = nhi - nlo + npad

! NOTE: These are global, must be set
g1ilo = glo(1)
g1jlo = glo(2)
g1klo = glo(3)

! Check that mesh sizes haven't changed
if (allocated(cgrn1)) then
  good_sizes = .true. 
  do i=1, 3
    if (lbound(cgrn1, i) /= glo(i)) good_sizes = .false.
    if (ubound(cgrn1, i) /= ghi(i)) good_sizes = .false.
  enddo
  
  if (.not. good_sizes) then
    !! print *, 'green function array size changed, deallocating'
    deallocate(cgrn1)
  endif
endif

if(.not. allocated(cgrn1)) allocate(cgrn1( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)) )

end subroutine osc_alloc_freespace_array


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
 function coulombfun(u,v,w,gam) result(res)
 implicit none
 real(dp) :: res
 real(dp) :: u,v,w,gam
 if(u.eq.0.d0 .and. v.eq.0.d0 .and. w.eq.0.d0)then
   res=0.d0
   return
 endif
 res=1.d0/sqrt(u**2+v**2+(gam*w)**2)  !coulomb
!     res=u/(u**2+v**2+(gam*w)**2)**1.5d0  !x-electric field
 return
 end function coulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfcoulombfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
!-!         real(dp), external :: lafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
!     res=1.d0/sqrt(u**2+v**2+w**2)  !coulomb
!     res=u/(u**2+v**2+w**2)**1.5d0  !x-electric field
res=lafun(x2,y2,z2)-lafun(x1,y2,z2)-lafun(x2,y1,z2)-lafun(x2,y2,z1)-lafun(x1,y1,z1)+ &
    lafun(x1,y1,z2)+lafun(x1,y2,z1)+lafun(x2,y1,z1)
res=res/(dx*dy*dz*gam)

end function igfcoulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function lafun(x,y,z) result(res)
! lafun is the function involving log and atan in the PRSTAB paper (I should find a better name for this function)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
res=-0.5d0*z**2*atan(x*y/(z*r))-0.5d0*y**2*atan(x*z/(y*r))-0.5d0*x**2*atan(y*z/(x*r)) &
     +y*z*log(x+r)+x*z*log(y+r)+x*y*log(z+r)

end function lafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfexfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !-13 causes NaN's
real(dp), parameter :: em=-1.d-10
!-!         real(dp), external :: xlafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=0.d0
! res=xlafun(em,em,em) &
!    -xlafun(x1,em,em)-xlafun(em,y1,em)-xlafun(em,em,z1)-xlafun(x1,y1,z1) &
!    +xlafun(x1,y1,em)+xlafun(x1,em,z1)+xlafun(em,y1,z1)+xlafun(x2,em,em) &
!    -xlafun(ep,em,em)-xlafun(x2,y1,em)-xlafun(x2,em,z1)-xlafun(ep,y1,z1) &
!    +xlafun(ep,y1,em)+xlafun(ep,em,z1)+xlafun(x2,y1,z1)+xlafun(ep,y2,ep) &
!    -xlafun(x1,y2,ep)-xlafun(ep,ep,ep)-xlafun(ep,y2,z1)-xlafun(x1,ep,z1) &
!    +xlafun(x1,ep,ep)+xlafun(x1,y2,z1)+xlafun(ep,ep,z1)+xlafun(ep,ep,z2) &
!    -xlafun(x1,ep,z2)-xlafun(ep,y1,z2)-xlafun(ep,ep,ep)-xlafun(x1,y1,ep) &
!    +xlafun(x1,y1,z2)+xlafun(x1,ep,ep)+xlafun(ep,y1,ep)+xlafun(x2,y2,ep) &
!    -xlafun(ep,y2,ep)-xlafun(x2,ep,ep)-xlafun(x2,y2,z1)-xlafun(ep,ep,z1) &
!    +xlafun(ep,ep,ep)+xlafun(ep,y2,z1)+xlafun(x2,ep,z1)+xlafun(x2,ep,z2) &
!    -xlafun(ep,ep,z2)-xlafun(x2,y1,z2)-xlafun(x2,ep,ep)-xlafun(ep,y1,ep) &
!    +xlafun(ep,y1,z2)+xlafun(ep,ep,ep)+xlafun(x2,y1,ep)+xlafun(ep,y2,z2) &
!    -xlafun(x1,y2,z2)-xlafun(ep,ep,z2)-xlafun(ep,y2,ep)-xlafun(x1,ep,ep) &
!    +xlafun(x1,ep,z2)+xlafun(x1,y2,ep)+xlafun(ep,ep,ep)+xlafun(x2,y2,z2) &
!    -xlafun(ep,y2,z2)-xlafun(x2,ep,z2)-xlafun(x2,y2,ep)-xlafun(ep,ep,ep) &
!    +xlafun(ep,ep,z2)+xlafun(ep,y2,ep)+xlafun(x2,ep,ep)
else
  res=xlafun(x2,y2,z2)-xlafun(x1,y2,z2)-xlafun(x2,y1,z2)-xlafun(x2,y2,z1) &
      -xlafun(x1,y1,z1)+xlafun(x1,y1,z2)+xlafun(x1,y2,z1)+xlafun(x2,y1,z1)
endif
res=res/(dx*dy*dz)
return
end function igfexfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfeyfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !should include em also, but probably doesn't matter
!-!         real(dp), external :: ylafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=0.d0
! res=ylafun(ep,ep,ep) &
!    -ylafun(x1,ep,ep)-ylafun(ep,y1,ep)-ylafun(ep,ep,z1)-ylafun(x1,y1,z1) &
!    +ylafun(x1,y1,ep)+ylafun(x1,ep,z1)+ylafun(ep,y1,z1)+ylafun(x2,ep,ep) &
!    -ylafun(ep,ep,ep)-ylafun(x2,y1,ep)-ylafun(x2,ep,z1)-ylafun(ep,y1,z1) &
!    +ylafun(ep,y1,ep)+ylafun(ep,ep,z1)+ylafun(x2,y1,z1)+ylafun(ep,y2,ep) &
!    -ylafun(x1,y2,ep)-ylafun(ep,ep,ep)-ylafun(ep,y2,z1)-ylafun(x1,ep,z1) &
!    +ylafun(x1,ep,ep)+ylafun(x1,y2,z1)+ylafun(ep,ep,z1)+ylafun(ep,ep,z2) &
!    -ylafun(x1,ep,z2)-ylafun(ep,y1,z2)-ylafun(ep,ep,ep)-ylafun(x1,y1,ep) &
!    +ylafun(x1,y1,z2)+ylafun(x1,ep,ep)+ylafun(ep,y1,ep)+ylafun(x2,y2,ep) &
!    -ylafun(ep,y2,ep)-ylafun(x2,ep,ep)-ylafun(x2,y2,z1)-ylafun(ep,ep,z1) &
!    +ylafun(ep,ep,ep)+ylafun(ep,y2,z1)+ylafun(x2,ep,z1)+ylafun(x2,ep,z2) &
!    -ylafun(ep,ep,z2)-ylafun(x2,y1,z2)-ylafun(x2,ep,ep)-ylafun(ep,y1,ep) &
!    +ylafun(ep,y1,z2)+ylafun(ep,ep,ep)+ylafun(x2,y1,ep)+ylafun(ep,y2,z2) &
!    -ylafun(x1,y2,z2)-ylafun(ep,ep,z2)-ylafun(ep,y2,ep)-ylafun(x1,ep,ep) &
!    +ylafun(x1,ep,z2)+ylafun(x1,y2,ep)+ylafun(ep,ep,ep)+ylafun(x2,y2,z2) &
!    -ylafun(ep,y2,z2)-ylafun(x2,ep,z2)-ylafun(x2,y2,ep)-ylafun(ep,ep,ep) &
!    +ylafun(ep,ep,z2)+ylafun(ep,y2,ep)+ylafun(x2,ep,ep)
else
  res=ylafun(x2,y2,z2)-ylafun(x1,y2,z2)-ylafun(x2,y1,z2)-ylafun(x2,y2,z1) &
     -ylafun(x1,y1,z1)+ylafun(x1,y1,z2)+ylafun(x1,y2,z1)+ylafun(x2,y1,z1)
endif

res=res/(dx*dy*dz)

end function igfeyfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
function igfezfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !should include em also, but probably doesn't matter
!-!         real(dp), external :: zlafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=0.d0
! res=zlafun(ep,ep,ep) &
!    -zlafun(x1,ep,ep)-zlafun(ep,y1,ep)-zlafun(ep,ep,z1)-zlafun(x1,y1,z1) &
!    +zlafun(x1,y1,ep)+zlafun(x1,ep,z1)+zlafun(ep,y1,z1)+zlafun(x2,ep,ep) &
!    -zlafun(ep,ep,ep)-zlafun(x2,y1,ep)-zlafun(x2,ep,z1)-zlafun(ep,y1,z1) &
!    +zlafun(ep,y1,ep)+zlafun(ep,ep,z1)+zlafun(x2,y1,z1)+zlafun(ep,y2,ep) &
!    -zlafun(x1,y2,ep)-zlafun(ep,ep,ep)-zlafun(ep,y2,z1)-zlafun(x1,ep,z1) &
!    +zlafun(x1,ep,ep)+zlafun(x1,y2,z1)+zlafun(ep,ep,z1)+zlafun(ep,ep,z2) &
!    -zlafun(x1,ep,z2)-zlafun(ep,y1,z2)-zlafun(ep,ep,ep)-zlafun(x1,y1,ep) &
!    +zlafun(x1,y1,z2)+zlafun(x1,ep,ep)+zlafun(ep,y1,ep)+zlafun(x2,y2,ep) &
!    -zlafun(ep,y2,ep)-zlafun(x2,ep,ep)-zlafun(x2,y2,z1)-zlafun(ep,ep,z1) &
!    +zlafun(ep,ep,ep)+zlafun(ep,y2,z1)+zlafun(x2,ep,z1)+zlafun(x2,ep,z2) &
!    -zlafun(ep,ep,z2)-zlafun(x2,y1,z2)-zlafun(x2,ep,ep)-zlafun(ep,y1,ep) &
!    +zlafun(ep,y1,z2)+zlafun(ep,ep,ep)+zlafun(x2,y1,ep)+zlafun(ep,y2,z2) &
!    -zlafun(x1,y2,z2)-zlafun(ep,ep,z2)-zlafun(ep,y2,ep)-zlafun(x1,ep,ep) &
!    +zlafun(x1,ep,z2)+zlafun(x1,y2,ep)+zlafun(ep,ep,ep)+zlafun(x2,y2,z2) &
!    -zlafun(ep,y2,z2)-zlafun(x2,ep,z2)-zlafun(x2,y2,ep)-zlafun(ep,ep,ep) &
!    +zlafun(ep,ep,z2)+zlafun(ep,y2,ep)+zlafun(x2,ep,ep)
else
  res=zlafun(x2,y2,z2)-zlafun(x1,y2,z2)-zlafun(x2,y1,z2)-zlafun(x2,y2,z1) &
     -zlafun(x1,y1,z1)+zlafun(x1,y1,z2)+zlafun(x1,y2,z1)+zlafun(x2,y1,z1)
endif
res=res/(dx*dy*dz*gam) !note the factor of gam in the denominator here, as needed for Ez

end function igfezfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function xlafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
    !  if (z+r == 0 ) print *, 'ERROR:x, r, x+r', x, y, z,r, x+r
    !  if (y+r == 0 ) print *, 'ERROR:y, r, y+r', x, y, z,r, y+r
    !  if (z+r == 0 ) print *, 'ERROR:z, r, z+r', x, y, z,r, z+r
!     res=x*atan(y*z/(r*x))-z*atanh(r/y)-y*atanh(r/z)
!res=z-x*atan(z/x)+x*atan(y*z/(x*r))-z*log(y+r)-y*log(z+r)
res=z-x*atan(z/x)+x*atan(y*z/(x*r))
if (y+r /= 0) res = res -z*log(y+r)
if (z+r /= 0) res = res -y*log(z+r)
!     write(2,'(5(1pe12.5,1x))')x,y,z,r,res

end function xlafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function ylafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
!res=x-y*atan(x/y)+y*atan(z*x/(y*r))-x*log(z+r)-z*log(x+r)
res=x-y*atan(x/y)+y*atan(z*x/(y*r))
if (z+r /= 0) res = res -x*log(z+r)
if (x+r /= 0) res = res -z*log(x+r)   

end function ylafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function zlafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
!res=y-z*atan(y/z)+z*atan(x*y/(z*r))-y*log(x+r)-x*log(y+r)

res=y-z*atan(y/z)+z*atan(x*y/(z*r))
if (x+r /= 0) res = res -y*log(x+r)
if (y+r /= 0) res = res -x*log(y+r)

! TEST
!res = z*atan(x*y,(r*z)) - y*atanh(r/x) - x*atanh(r/y)      ! Bad
!res = z * atan(x*y/z*r) + y*log(1-x/r)/2 - y*log(1+x/r)/2  ! Works
!res = -z**2 * (x/(x**2 + z**2) - y/(y**2 + z**2)) + z*atan(x*y/(z*r)) -y*log(x+r) - x*log(y+r) ! Works

end function zlafun


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine getrhotilde(rho,crho,ilo,jlo,klo)
implicit none
integer :: ilo,jlo,klo
real(dp), dimension(ilo:,jlo:,klo:) :: rho
complex(dp), dimension(ilo:,jlo:,klo:) :: crho
integer :: ihi,jhi,khi,ihi2,jhi2,khi2,iperiod,jperiod,kperiod
ihi=ubound(rho,1); jhi=ubound(rho,2); khi=ubound(rho,3)
ihi2=ubound(crho,1); jhi2=ubound(crho,2); khi2=ubound(crho,3)
iperiod=size(crho,1); jperiod=size(crho,2); kperiod=size(crho,3)
crho(ilo:ihi2,jlo:jhi2,klo:khi2)=0.d0 !zero everything before storing rho in one octant
crho(ilo:ihi,jlo:jhi,klo:khi)=rho(ilo:ihi,jlo:jhi,klo:khi)
!fft the charge density:
call ccfft3d(crho,crho,(/1,1,1/),iperiod,jperiod,kperiod,0) !can I do the fft in place? problematic in parallel code
return
end subroutine
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_getgrnfree(cgrn,delta,gam,ilo_grn,jlo_grn,klo_grn,npad,icomp,igfflag)
implicit none
real(dp), dimension(3) :: delta
integer, dimension(3) :: npad
real(dp) :: dx,dy,dz,gam
integer :: ilo_grn,jlo_grn,klo_grn,igfflag,icomp
integer :: ilo_grn_gbl,jlo_grn_gbl,klo_grn_gbl !in a parallel code these would be in the arg list
integer :: ipad,jpad,kpad,ishift,jshift,kshift,iperiod,jperiod,kperiod
complex(dp), dimension(ilo_grn:,jlo_grn:,klo_grn:) :: cgrn
integer :: i,j,k,ip,jp,kp
real(dp) :: u,v,w
real(dp) :: gval
dx=delta(1); dy=delta(2); dz=delta(3)
ilo_grn_gbl=ilo_grn; jlo_grn_gbl=jlo_grn; klo_grn_gbl=klo_grn
iperiod=size(cgrn,1); jperiod=size(cgrn,2); kperiod=size(cgrn,3)
ipad=npad(1); jpad=npad(2); kpad=npad(3)
!this puts the Green function where it's needed so the convolution ends up in the correct location in the array
ishift=iperiod/2-(ipad+1)/2; jshift=jperiod/2-(jpad+1)/2; kshift=kperiod/2-(kpad+1)/2
!$ print *, 'OpenMP Green function calc'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k=klo_grn,ubound(cgrn,3)
  kp=klo_grn_gbl+mod(k-klo_grn_gbl+kshift ,kperiod) !in the serial code klo_grn_gbl = klo_grn, similarly for i,j
  w=kp*dz
  do j=jlo_grn,ubound(cgrn,2)
    jp=jlo_grn_gbl+mod(j-jlo_grn_gbl+jshift ,jperiod)
    v=jp*dy
   do i=ilo_grn,ubound(cgrn,1)
     ip=ilo_grn_gbl+mod(i-ilo_grn_gbl+ishift ,iperiod)
     u=ip*dx
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn(i,j,k)= cmplx(gval,0.d0, dp)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call ccfft3d(cgrn,cgrn,(/1,1,1/),iperiod,jperiod,kperiod,0)

end subroutine osc_getgrnfree

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine conv3d(crho,qgrn1,con,rilo,rjlo,rklo,g1ilo,g1jlo,g1klo,cilo,cjlo,cklo,iperiod,jperiod,kperiod)
implicit none
integer :: rilo,rjlo,rklo,g1ilo,g1jlo,g1klo,cilo,cjlo,cklo,iperiod,jperiod,kperiod
!input arrays:
complex(dp), dimension(rilo:,rjlo:,rklo:) :: crho
complex(dp), dimension(g1ilo:,g1jlo:,g1klo:) :: qgrn1
real(dp), dimension(cilo:,cjlo:,cklo:) :: con
! local:
complex(dp), allocatable, dimension(:,:,:) :: ccon
real(dp) :: fpei,qtot,factr
integer :: cihi,cjhi,ckhi

fpei=299792458.d0**2*1.d-7  ! this is 1/(4 pi eps0)
qtot=1.d0 !fix later: 1.d-9 ! 1 nC
allocate(ccon(cilo:cilo+iperiod-1,cjlo:cjlo+jperiod-1,cklo:cklo+kperiod-1))
ccon(:,:,:)=crho(:,:,:)*qgrn1(:,:,:)
call ccfft3d(ccon,ccon,(/-1,-1,-1/),iperiod,jperiod,kperiod,0)
!normalize:
!     factr=hx*hy*hz/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot
factr=    1.d0/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot
!store final result in original size (not double size) real array:
cihi=ubound(con,1)
cjhi=ubound(con,2)
ckhi=ubound(con,3)
con(cilo:cihi,cjlo:cjhi,cklo:ckhi)=factr*real(ccon(cilo:cihi,cjlo:cjhi,cklo:ckhi),dp)

end subroutine conv3d


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_rectpipe_solver(rho,a,b,gam0,delta,umin,phi,efield,bfield,nlo,nhi,nlo_gbl,nhi_gbl,idirectfieldcalc,igfflag)
!This assumes that cgrn1,...,cgrn2 have been allocated and calculated prior to entering this routine.
!This assumes that the rho and phi arrays have the same dimensions, called ilo,ihi,jlo,jhi,klo,khi below.
!The "gbl" variables are from the parallel code, so they are ignored in this serial code.
implicit none
integer, intent(in) :: igfflag,idirectfieldcalc
integer :: ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl
real(dp), intent(in) :: a,b,gam0
real(dp), intent(in), dimension(3) :: delta,umin
integer, intent(in), dimension(3) :: nlo,nhi,nlo_gbl,nhi_gbl
real(dp), intent(in), dimension(:,:,:) :: rho
real(dp), intent(out), dimension(:,:,:) :: phi
real(dp), intent(out), dimension(:,:,:,:) :: efield,bfield
complex(dp), allocatable, dimension(:,:,:) :: crho
integer :: iperiod,jperiod,kperiod
real(dp) :: dx,dy,dz,xmin,ymin,zmin
integer :: ilo2,ihi2,jlo2,jhi2,klo2,khi2

integer :: icomp
real(dp), parameter :: clight=299792458.d0
integer :: i,j,k,im1,ip1,jm1,jp1,km1,kp1
real(dp) :: beta0,xfac,yfac,zfac

dx=delta(1); dy=delta(2); dz=delta(3)
xmin=umin(1); ymin=umin(2); zmin=umin(3)
ilo=nlo(1);ihi=nhi(1);jlo=nlo(2);jhi=nhi(2);klo=nlo(3);khi=nhi(3)

ilo2=ilo; jlo2=jlo; klo2=klo
iperiod=size(cgrn1,1); jperiod=size(cgrn1,2); kperiod=size(cgrn1,3)
ihi2=ilo2+iperiod-1; jhi2=jlo2+jperiod-1; khi2=klo2+kperiod-1

write(6,*)'iperiod,jperiod,kperiod=',iperiod,jperiod,kperiod
 
if(idirectfieldcalc.eq.0)then
  write(6,*)'Solving for Phi'
  icomp=0

!compute/store the FFT of the charge density:
if(.not.allocated(crho))allocate(crho(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
call getrhotilde(rho,crho,ilo,jlo,klo)
write(6,*)'done computing FFT of the charge density'
!compute/store the FFT of the green function: this has been moved to the main program
if(.not.allocated(cgrn1))then
  write(6,*)'something wrong. rectpipe green functions should have been computed during initialization'
  stop
endif
call fftconvcorr3d(crho,cgrn1,cgrn2,cgrn3,cgrn4,phi,ilo,jlo,klo, &
     &g1ilo,g1jlo,g1klo,g2ilo,g2jlo,g2klo,g3ilo,g3jlo,g3klo,g4ilo,g4jlo,g4klo, &
     &ilo,jlo,klo,iperiod,jperiod,kperiod,dx,dy,dz,xmin,ymin,zmin)
write(6,*)'finished convolutions and correlations'

!original version did not handle the edge values correctly when differencing; here is a first order approx:
do k=lbound(phi,3),ubound(phi,3)
  km1=k-1; if(km1.lt.lbound(phi,3))km1=k
  kp1=k+1; if(kp1.gt.ubound(phi,3))kp1=k
  zfac=merge(2.d0,1.d0,(km1.eq.k).or.(kp1.eq.k))
  do j=lbound(phi,2),ubound(phi,2)
    jm1=j-1; if(jm1.lt.lbound(phi,2))jm1=j
    jp1=j+1; if(jp1.gt.ubound(phi,2))jp1=j
    yfac=merge(2.d0,1.d0,(jm1.eq.j).or.(jp1.eq.j))
      do i=lbound(phi,1),ubound(phi,1)
        im1=i-1; if(im1.lt.lbound(phi,1))im1=i
        ip1=i+1; if(ip1.gt.ubound(phi,1))ip1=i
        xfac=merge(2.d0,1.d0,(im1.eq.i).or.(ip1.eq.i))
        efield(i,j,k,1)=-(phi(ip1,j,k)-phi(im1,j,k))/(2.d0*dx)*gam0*xfac
        efield(i,j,k,2)=-(phi(i,jp1,k)-phi(i,jm1,k))/(2.d0*dy)*gam0*yfac
        efield(i,j,k,3)=-(phi(i,j,kp1)-phi(i,j,km1))/(2.d0*dz)/gam0*zfac
      enddo
    enddo
  enddo
endif

if(idirectfieldcalc.eq.1)then
  write(6,*)'**********ERROR********* rectpipe implementation presently requires directfieldcalc=0'
endif

! set the magnetic field:
beta0=sqrt(1-1/gam0**2)
do k=lbound(phi,3),ubound(phi,3)
  do j=lbound(phi,2),ubound(phi,2)
    do i=lbound(phi,1),ubound(phi,1)
      bfield(i,j,k,1)=-efield(i,j,k,2)*beta0/clight
      bfield(i,j,k,2)= efield(i,j,k,1)*beta0/clight
      bfield(i,j,k,3)=0.d0
    enddo
  enddo
enddo

end subroutine osc_rectpipe_solver



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_getgrnpipe(gam,a,b,delta,umin,npad)
!makes use of arrays cgrn1,cgrn2,cgrn3,cgrn4 declared at the beginning of this module
implicit none
real(dp), dimension(3) :: delta,umin
integer, dimension(3) :: npad
integer :: ipad,jpad,kpad,ishift,jshift,kshift,iperiod,jperiod,kperiod
real(dp) :: gam,a,b,hx,hy,hz,xmin,ymin,zmin
real(dp) :: u,v,w
integer :: i,j,k,ip,jp,kp
!     real(dp), external ::rfun

hx=delta(1); hy=delta(2); hz=delta(3)
xmin=umin(1); ymin=umin(2); zmin=umin(3)
!
iperiod=size(cgrn1,1); jperiod=size(cgrn1,2); kperiod=size(cgrn1,3)
ipad=npad(1); jpad=npad(2); kpad=npad(3)
!this puts the Green function where it's needed so the convolution ends up in the correct location in the array
      ishift=iperiod/2-(ipad+1)/2; jshift=jperiod/2-(jpad+1)/2; kshift=kperiod/2-(kpad+1)/2
!
!1:
do k=g1klo,ubound(cgrn1,3)
do j=g1jlo,ubound(cgrn1,2)
do i=g1ilo,ubound(cgrn1,1)
  ip=g1ilo+mod(i-g1ilo+ishift,iperiod)
  jp=g1jlo+mod(j-g1jlo+jshift,jperiod)
  kp=g1klo+mod(k-g1klo+kshift,kperiod)
  w=kp*hz
  v=jp*hy
  u=ip*hx
  cgrn1(i,j,k)=cmplx(0.d0,0.d0, dp)
  cgrn1(i,j,k)=rfun(u,v,w,gam,a,b,hz,1,1)
enddo
enddo
enddo
call ccfft3d(cgrn1,cgrn1,(/1,1,1/),iperiod,jperiod,kperiod,0)

!2:
do k=g2klo,ubound(cgrn2,3)
do j=g2jlo,ubound(cgrn2,2)
do i=g2ilo,ubound(cgrn2,1)
  ip=g2ilo+mod(i-g2ilo+ishift,iperiod)
  kp=g2klo+mod(k-g2klo+kshift,kperiod)
  w=kp*hz
  v=2.d0*ymin+(j-g2jlo)*hy ! v=2.d0*ymin+j*hy
  u=ip*hx
  cgrn2(i,j,k)=cmplx(0.d0,0.d0, dp)
  cgrn2(i,j,k)=rfun(u,v,w,gam,a,b,hz,1,-1)
enddo
enddo
enddo
call ccfft3d(cgrn2,cgrn2,(/1,-1,1/),iperiod,jperiod,kperiod,0)

!3:
do k=g3klo,ubound(cgrn3,3)
do j=g3jlo,ubound(cgrn3,2)
do i=g3ilo,ubound(cgrn3,1)
  jp=g3jlo+mod(j-g3jlo+jshift,jperiod)
  kp=g3klo+mod(k-g3klo+kshift,kperiod)
  w=kp*hz
  v=jp*hy
  u=2.d0*xmin+(i-g3ilo)*hx ! u=2.d0*xmin+i*hx
  cgrn3(i,j,k)=cmplx(0.d0,0.d0, dp)
  cgrn3(i,j,k)=rfun(u,v,w,gam,a,b,hz,-1,1)
enddo
enddo
enddo
call ccfft3d(cgrn3,cgrn3,(/-1,1,1/),iperiod,jperiod,kperiod,0)

!4:
do k=g4klo,ubound(cgrn4,3)
do j=g4jlo,ubound(cgrn4,2)
do i=g4ilo,ubound(cgrn4,1)
  kp=g4klo+mod(k-g4klo+kshift,kperiod)
  w=kp*hz
  v=2.d0*ymin+(j-g4jlo)*hy ! v=2.d0*ymin+j*hy
  u=2.d0*xmin+(i-g4ilo)*hx ! u=2.d0*xmin+i*hx
  cgrn4(i,j,k)=cmplx(0.d0,0.d0, dp)
  cgrn4(i,j,k)=rfun(u,v,w,gam,a,b,hz,-1,-1)
enddo
enddo
enddo
call ccfft3d(cgrn4,cgrn4,(/-1,-1,1/),iperiod,jperiod,kperiod,0)
return
end subroutine osc_getgrnpipe



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function rfun(u,v,w,gam,a,b,hz,i,j) result(res)
!this is the IGF version (non-IGF line is commented out)
!this version contains (i**m)*(j**n)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,a,b,hz
integer :: i,j
real(dp), parameter :: pi=acos(-1.0d0)
real(dp) :: ainv,binv,piainv,pibinv
real(dp) :: zfun,kapmn,term
integer :: m,n
ainv=1.0d0/a
binv=1.0d0/b
piainv=pi*ainv
pibinv=pi*binv
!     rfun=0.d0
res=0.d0
do m=1,5
do n=1,5
  kapmn=sqrt((m*pi*ainv)**2+(n*pi*binv)**2)
!!!!!   zfun=exp(-kapmn*abs(gam*w))
!!!!!   zfun=exp(-kapmn*abs(w))
zfun=(exp(-kapmn*abs(gam*w-hz))-2.d0*exp(-kapmn*abs(gam*w))+exp(-kapmn*abs(gam*w+hz)))/(hz**2*kapmn**2)
!       zfun=(exp(-kapmn*abs(w-hz))-2.d0*exp(-kapmn*abs(w))+exp(-kapmn*abs(w+hz)))/(hz**2*kapmn**2)
if(w.eq.0)zfun=zfun+2.d0/(hz*kapmn)
term=(i**m)*(j**n)*cos(m*u*piainv)*cos(n*v*pibinv)*zfun/kapmn
!       rfun=rfun+term
res=res+term
enddo
enddo
res=res*2.d0*pi*ainv*binv !here is the correct normalization

end function rfun



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine fftconvcorr3d(crho,qgrn1,qgrn2,qgrn3,qgrn4,con,rilo,rjlo,rklo, &
     g1ilo,g1jlo,g1klo,g2ilo,g2jlo,g2klo,g3ilo,g3jlo,g3klo,g4ilo,g4jlo,g4klo, &
     cilo,cjlo,cklo,iperiod,jperiod,kperiod,hx,hy,hz,xmin,ymin,zmin)
implicit none
integer :: rilo,rjlo,rklo,g1ilo,g1jlo,g1klo,g2ilo,g2jlo,g2klo,g3ilo,g3jlo,g3klo,g4ilo,g4jlo,g4klo
real(dp) :: hx,hy,hz,xmin,ymin,zmin
integer :: cilo,cjlo,cklo,iperiod,jperiod,kperiod
integer :: cihi,cjhi,ckhi
real(dp) :: fpei,qtot,factr
!input arrays:
complex(dp), dimension(rilo:,rjlo:,rklo:) :: crho
complex(dp), dimension(g1ilo:,g1jlo:,g1klo:) :: qgrn1
complex(dp), dimension(g2ilo:,g2jlo:,g2klo:) :: qgrn2
complex(dp), dimension(g3ilo:,g3jlo:,g3klo:) :: qgrn3
complex(dp), dimension(g4ilo:,g4jlo:,g4klo:) :: qgrn4
real(dp), dimension(cilo:,cjlo:,cklo:) :: con
!local arrays:
complex(dp), allocatable, dimension(:,:,:) :: ccon,ctmp


fpei=299792458.**2*1.e-7  ! this is 1/(4 pi eps0)
qtot=1.d0 !fix later: 1.d-9 ! 1 nC
!
allocate(ccon(cilo:cilo+iperiod-1,cjlo:cjlo+jperiod-1,cklo:cklo+kperiod-1))
allocate(ctmp(cilo:cilo+iperiod-1,cjlo:cjlo+jperiod-1,cklo:cklo+kperiod-1))
!
!1st term:
ccon(:,:,:)=crho(:,:,:)*qgrn1(:,:,:)
call ccfft3d(ccon,ccon,(/-1,-1,-1/),iperiod,jperiod,kperiod,0)

!2nd term:
ctmp(:,:,:)=crho(:,:,:)*qgrn2(:,:,:)
call ccfft3d(ctmp,ctmp,(/-1,1,-1/),iperiod,jperiod,kperiod,0)
ccon(:,:,:)=ccon(:,:,:)-ctmp(:,:,:)

!3rd term:
ctmp(:,:,:)=crho(:,:,:)*qgrn3(:,:,:)
call ccfft3d(ctmp,ctmp,(/1,-1,-1/),iperiod,jperiod,kperiod,0)
ccon(:,:,:)=ccon(:,:,:)-ctmp(:,:,:)

!4th term:
ctmp(:,:,:)=crho(:,:,:)*qgrn4(:,:,:)
call ccfft3d(ctmp,ctmp,(/1,1,-1/),iperiod,jperiod,kperiod,0)
ccon(:,:,:)=ccon(:,:,:)+ctmp(:,:,:)

!normalize:
!     factr=hx*hy*hz/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot
factr=    1.d0/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot

!store final result in original size (not double size) real array:
cihi=ubound(con,1)
cjhi=ubound(con,2)
ckhi=ubound(con,3)
con(cilo:cihi,cjlo:cjhi,cklo:ckhi)=factr*real(ccon(cilo:cihi,cjlo:cjhi,cklo:ckhi), dp)
return
end subroutine fftconvcorr3d


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_read_rectpipe_grn
implicit none
real(dp) :: apipeval,bpipeval,dxval,dyval,dzval,xminval,yminval,zminval,xmaxval,ymaxval,zmaxval,gamval
integer :: iloval,jloval,kloval,ihival,jhival,khival,iperval,jperval,kperval
open(unit=2,file='rectpipegrn.dat',form='unformatted',status='old')
read(2)apipeval,bpipeval,dxval,dyval,dzval,xminval,yminval,zminval,xmaxval,ymaxval,zmaxval, &
       iloval,jloval,kloval,ihival,jhival,khival,iperval,jperval,kperval,gamval
read(2)cgrn1,cgrn2,cgrn3,cgrn4
close(2)
write(6,*)'Done reading rectpipe transformed Green functions with these parameters:'
write(6,*)'apipe,bpipe=',apipeval,bpipeval
write(6,*)'dx,dy,dz=',dxval,dyval,dzval
write(6,*)'xmin,ymin,zmin=',xminval,yminval,zminval
write(6,*)'xmax,ymax,zmax=',xmaxval,ymaxval,zmaxval
write(6,*)'ilo,jlo,klo=',iloval,jloval,kloval
write(6,*)'ihi,jhi,khi=',ihival,jhival,khival
write(6,*)'convolution iperiod,jperiod,kperiod=',iperval,jperval,kperval
write(6,*)'gamma=',gamval

end subroutine osc_read_rectpipe_grn


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_write_rectpipe_grn(apipe,bpipe,delta,umin,umax,nlo,nhi,gamma)
implicit none
real(dp) :: apipe,bpipe,gamma
real(dp), dimension(3) :: delta,umin,umax
integer, dimension(3) :: nlo,nhi
real(dp) :: dx,dy,dz,xmin,ymin,zmin,xmax,ymax,zmax
integer :: nxlo,nylo,nzlo,nxhi,nyhi,nzhi
dx=delta(1); dy=delta(2); dz=delta(3)
xmin=umin(1); ymin=umin(2); zmin=umin(3)
xmax=umax(1); ymax=umax(2); zmax=umax(3)
nxlo=nlo(1); nylo=nlo(2); nzlo=nlo(3)
nxhi=nhi(1); nyhi=nhi(2); nzhi=nhi(3)
open(unit=2,file='rectpipegrn.dat',form='unformatted',status='unknown')
write(2)apipe,bpipe,dx,dy,dz,xmin,ymin,zmin,xmax,ymax,zmax,nxlo,nylo,nzlo,nxhi,nyhi,nzhi, &
        size(cgrn1,1),size(cgrn1,2),size(cgrn1,3),gamma
write(2)cgrn1,cgrn2,cgrn3,cgrn4
close(2)
write(6,*)'done writing rectpipe transformed Green functions'
return
end subroutine osc_write_rectpipe_grn

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_alloc_rectpipe_arrays(nlo,nhi,npad)
implicit none
integer, intent(in), dimension(3) :: nlo,nhi,npad
integer :: rilo,rihi,rjlo,rjhi,rklo,rkhi !rho dimensions
integer :: cilo,cihi,cjlo,cjhi,cklo,ckhi !phi dimensions
integer :: ipad,jpad,kpad
integer :: g1ihi,g1jhi,g1khi,g2ihi,g2jhi,g2khi,g3ihi,g3jhi,g3khi,g4ihi,g4jhi,g4khi

rilo=nlo(1); rihi=nhi(1); rjlo=nlo(2); rjhi=nhi(2); rklo=nlo(3); rkhi=nhi(3)
cilo=nlo(1); cihi=nhi(1); cjlo=nlo(2); cjhi=nhi(2); cklo=nlo(3); ckhi=nhi(3)
ipad=npad(1); jpad=npad(2); kpad=npad(3)

g1ilo=cilo-rihi; g1ihi=cihi-rilo; g1jlo=cjlo-rjhi; g1jhi=cjhi-rjlo; g1klo=cklo-rkhi; g1khi=ckhi-rklo
g2ilo=cilo-rihi; g2ihi=cihi-rilo; g2jlo=cjlo+rjlo; g2jhi=cjhi+rjhi; g2klo=cklo-rkhi; g2khi=ckhi-rklo
g3ilo=cilo+rilo; g3ihi=cihi+rihi; g3jlo=cjlo-rjhi; g3jhi=cjhi-rjlo; g3klo=cklo-rkhi; g3khi=ckhi-rklo
g4ilo=cilo+rilo; g4ihi=cihi+rihi; g4jlo=cjlo+rjlo; g4jhi=cjhi+rjhi; g4klo=cklo-rkhi; g4khi=ckhi-rklo

if(.not.allocated(cgrn1))allocate(cgrn1(g1ilo:g1ihi+ipad,g1jlo:g1jhi+jpad,g1klo:g1khi+kpad))
if(.not.allocated(cgrn2))allocate(cgrn2(g2ilo:g2ihi+ipad,g2jlo:g2jhi+jpad,g2klo:g2khi+kpad))
if(.not.allocated(cgrn3))allocate(cgrn3(g3ilo:g3ihi+ipad,g3jlo:g3jhi+jpad,g3klo:g3khi+kpad))
if(.not.allocated(cgrn4))allocate(cgrn4(g4ilo:g4ihi+ipad,g4jlo:g4jhi+jpad,g4klo:g4khi+kpad))

end subroutine osc_alloc_rectpipe_arrays

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine osc_cathodeimages_solver(rho,gam0,umin,delta,phi,efield,bfield,nlo,nhi,nlo_gbl,nhi_gbl,npad,idirectfieldcalc,igfflag,imethod)
! This will check the allocation of the global cgrn1. 
! OLD: note well: this assumes that osc_alloc_freespace_array has been called so that cgrn1 is already allocated!
! 
!use mpi
implicit none
integer, intent(in) :: igfflag,idirectfieldcalc,imethod
real(dp), intent(in), dimension(3) :: umin,delta
integer, intent(in), dimension(3) :: nlo,nhi,nlo_gbl,nhi_gbl,npad
real(dp), intent(in) :: gam0
!!!real(dp), intent(in), dimension(:,:,:) :: rho !routine changes rho when imethod=2; could be avoided w/ extra array if desired
real(dp), dimension(:,:,:) :: rho
real(dp), intent(out), dimension(:,:,:) :: phi
real(dp), intent(out), dimension(:,:,:,:) :: efield,bfield
real(dp) :: dx,dy,dz,r1,r2,zshft,umax3
integer :: ilo,ihi,jlo,jhi,klo,khi
integer :: ilo2,ihi2,jlo2,jhi2,klo2,khi2,iperiod,jperiod,kperiod,khalflen
real(dp), allocatable, dimension(:,:,:) :: phiimg
real(dp), allocatable, dimension(:,:,:,:) :: efieldimg
complex(dp), allocatable, dimension(:,:,:) :: crho

integer :: icomp,i,j,k,im1,ip1,jm1,jp1,km1,kp1
real(dp) :: beta0,xfac,yfac,zfac
integer :: mprocs,myrank,ierr
real(dp), parameter :: clight=299792458.d0
!
!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
myrank=0

if (igfflag == 0 .and. idirectfieldcalc.eq.1) then
  print *, 'Error: Direct field calc must use integrated Green function'
  print *, 'Aborting...'
  return
endif

dx=delta(1); dy=delta(2); dz=delta(3)
ilo=nlo(1);ihi=nhi(1);jlo=nlo(2);jhi=nhi(2);klo=nlo(3);khi=nhi(3)
ilo2=ilo; jlo2=jlo; klo2=klo

! Always call this. It will check the cgrn1 size, and re-allocate if things changed. 
call osc_alloc_freespace_array(nlo, nhi, npad)
call osc_alloc_image_array(nlo,nhi,npad)

iperiod=size(cgrn1,1); jperiod=size(cgrn1,2); kperiod=size(cgrn1,3)
ihi2=ilo2+iperiod-1; jhi2=jlo2+jperiod-1; khi2=klo2+kperiod-1
allocate(crho(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size complex array for the charge density
!! write(6,*)'iperiod,jperiod,kperiod=',iperiod,jperiod,kperiod

if(idirectfieldcalc.eq.0)then
  if(myrank.eq.0)write(6,*)'Solving for Phi including cathode images'
  icomp=0
!compute/store the FFT of the charge density:
call getrhotilde(rho,crho,ilo,jlo,klo)
!compute/store the FFT of the green function:
call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
!compute the convolution:   (note that crho,cgrn1 are complex padded, phi is real not padded)
call conv3d(crho,cgrn1,phi,ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)

!cathode image charges:
  allocate(phiimg(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3)))
if(imethod.eq.1)then ! convolution/correlation method
  call osc_getgrnimageconvcorr(cgrn2,umin,delta,gam0,g2ilo,g2jlo,g2klo,nlo,npad,icomp,igfflag) !cgrn2,cgrn1 are same except for domain
  call imageconvcorr3d(crho,cgrn2,phiimg,ilo,jlo,klo,g2ilo,g2jlo,g2klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  phi(:,:,:)=phi(:,:,:)-phiimg(:,:,:) !the image charges are negative
else
  write(6,*)'shift method...'
  khalflen=(nhi(3)-nlo(3)+1)/2
  do i=nlo(1),nhi(1)
  do j=nlo(2),nhi(2)
  do k=nlo(3),nlo(3)+khalflen-1
    r1=rho(i,j,k)
    r2=rho(i,j,nhi(3)+nlo(3)-k)
    rho(i,j,k)=r2                !used to say rhoflip
    rho(i,j,nhi(3)+nlo(3)-k)=r1  !used to say rhoflip
  enddo
  enddo
  enddo
  call getrhotilde(rho,crho,ilo,jlo,klo)
  umax3=umin(3)+(nhi(3)-nlo(3))*delta(3)
  zshft=umin(3)+umax3
  call osc_getgrnimageshift(cgrn1,delta,gam0,zshft,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
  call conv3d(crho,cgrn1,phiimg,ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  phi(:,:,:)=phi(:,:,:)-phiimg(:,:,:) !the image charges are negative

endif

!original version did not handle the edge values correctly when differencing; here is a first order approx:
  allocate(efieldimg(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),3))
  do k=lbound(phi,3),ubound(phi,3)
    km1=k-1; if(km1.lt.lbound(phi,3))km1=k
    kp1=k+1; if(kp1.gt.ubound(phi,3))kp1=k
    zfac=merge(2.d0,1.d0,(km1.eq.k).or.(kp1.eq.k))
    do j=lbound(phi,2),ubound(phi,2)
      jm1=j-1; if(jm1.lt.lbound(phi,2))jm1=j
      jp1=j+1; if(jp1.gt.ubound(phi,2))jp1=j
      yfac=merge(2.d0,1.d0,(jm1.eq.j).or.(jp1.eq.j))
      do i=lbound(phi,1),ubound(phi,1)
        im1=i-1; if(im1.lt.lbound(phi,1))im1=i
        ip1=i+1; if(ip1.gt.ubound(phi,1))ip1=i
        xfac=merge(2.d0,1.d0,(im1.eq.i).or.(ip1.eq.i))
          efield(i,j,k,1)=-(phi(ip1,j,k)-phi(im1,j,k))/(2.d0*dx)*gam0*xfac
          efield(i,j,k,2)=-(phi(i,jp1,k)-phi(i,jm1,k))/(2.d0*dy)*gam0*yfac
          efield(i,j,k,3)=-(phi(i,j,kp1)-phi(i,j,km1))/(2.d0*dz)/gam0*zfac
          efieldimg(i,j,k,1)=-(phiimg(ip1,j,k)-phiimg(im1,j,k))/(2.d0*dx)*gam0*xfac !needed below for B calc
          efieldimg(i,j,k,2)=-(phiimg(i,jp1,k)-phiimg(i,jm1,k))/(2.d0*dy)*gam0*yfac
          efieldimg(i,j,k,3)=-(phiimg(i,j,kp1)-phiimg(i,j,km1))/(2.d0*dz)/gam0*zfac
      enddo
    enddo
  enddo
endif

if(idirectfieldcalc.eq.1)then
  if(myrank.eq.0)write(6,*)'Solving for Ex, Ey, Ez including images on mesh:', nhi
  if(.not.allocated(crho))allocate(crho(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
  call getrhotilde(rho,crho,ilo,jlo,klo)
  if(.not.allocated(cgrn1))then
    write(6,*)'something wrong. freespace green function should have been allocated during initialization'
    stop
  endif
  do icomp=1,3
  call osc_getgrnfree(cgrn1,delta,gam0,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
  call conv3d(crho,cgrn1,efield(:,:,:,icomp),ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  enddo

!cathode image charges:
  allocate(efieldimg(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),3))
if(imethod.eq.1)then ! convolution/correlation method
  icomp=1
  call osc_getgrnimageconvcorr(cgrn2,umin,delta,gam0,g2ilo,g2jlo,g2klo,nlo,npad,icomp,igfflag) !cgrn2,cgrn1 are same except for domain
  call imageconvcorr3d(crho,cgrn2,efieldimg(:,:,:,1),ilo,jlo,klo,g2ilo,g2jlo,g2klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  efield(:,:,:,1)=efield(:,:,:,1)-efieldimg(:,:,:,1)
  icomp=2
  call osc_getgrnimageconvcorr(cgrn2,umin,delta,gam0,g2ilo,g2jlo,g2klo,nlo,npad,icomp,igfflag) !cgrn2,cgrn1 are same except for domain
  call imageconvcorr3d(crho,cgrn2,efieldimg(:,:,:,2),ilo,jlo,klo,g2ilo,g2jlo,g2klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  efield(:,:,:,2)=efield(:,:,:,2)-efieldimg(:,:,:,2)
  icomp=3
  call osc_getgrnimageconvcorr(cgrn2,umin,delta,gam0,g2ilo,g2jlo,g2klo,nlo,npad,icomp,igfflag) !cgrn2,cgrn1 are same except for domain
  call imageconvcorr3d(crho,cgrn2,efieldimg(:,:,:,3),ilo,jlo,klo,g2ilo,g2jlo,g2klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
  efield(:,:,:,3)=efield(:,:,:,3)-efieldimg(:,:,:,3)
else
  write(6,*)'direct field calc with shift method...'
  khalflen=(nhi(3)-nlo(3)+1)/2
  do i=nlo(1),nhi(1)
  do j=nlo(2),nhi(2)
  do k=nlo(3),nlo(3)+khalflen-1
    r1=rho(i,j,k)
    r2=rho(i,j,nhi(3)+nlo(3)-k)
    rho(i,j,k)=r2                !used to say rhoflip
    rho(i,j,nhi(3)+nlo(3)-k)=r1  !used to say rhoflip
  enddo
  enddo
  enddo
  call getrhotilde(rho,crho,ilo,jlo,klo)
  umax3=umin(3)+(nhi(3)-nlo(3))*delta(3)
  zshft=umin(3)+umax3
  do icomp=1,3
    call osc_getgrnimageshift(cgrn1,delta,gam0,zshft,g1ilo,g1jlo,g1klo,npad,icomp,igfflag)
    call conv3d(crho,cgrn1,efieldimg(:,:,:,icomp),ilo,jlo,klo,g1ilo,g1jlo,g1klo,ilo,jlo,klo,iperiod,jperiod,kperiod)
    efield(:,:,:,icomp)=efield(:,:,:,icomp)-efieldimg(:,:,:,icomp)
  enddo
endif ! imethod
endif ! idirectfieldcalc

! set the magnetic field:
beta0=sqrt(1-1/gam0**2)
do k=lbound(phi,3),ubound(phi,3)
  do j=lbound(phi,2),ubound(phi,2)
    do i=lbound(phi,1),ubound(phi,1)
!     bfield(i,j,k,1)=-efield(i,j,k,2)*beta0/clight
!     bfield(i,j,k,2)= efield(i,j,k,1)*beta0/clight
      bfield(i,j,k,1)=-(efield(i,j,k,2)+2.d0*efieldimg(i,j,k,2))*beta0/clight !I subtracted efieldimg above,
      bfield(i,j,k,2)= (efield(i,j,k,1)+2.d0*efieldimg(i,j,k,1))*beta0/clight !so add it back 2x to get sum
      bfield(i,j,k,3)=0.d0
    enddo
  enddo
enddo

end subroutine osc_cathodeimages_solver

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
! Allocates the global cgrn2 double-sized mesh for use by cathode image algorithm
! (cgrn2, originally envisioned for rectangular pipe, is being used here for images instead)
!
! If it was already allocated, it will be resized
!
!
subroutine osc_alloc_image_array(nlo,nhi,npad)
implicit none
integer, intent(in), dimension(3) :: nlo, nhi, npad
integer :: i, glo(3), ghi(3)
logical :: good_sizes

! Set bounds
glo(1:2) = nlo(1:2) - nhi(1:2)
ghi(1:2) = nhi(1:2) - nlo(1:2) + npad(1:2)
glo(3) =2*nlo(3)
ghi(3) =2*nhi(3) + npad(3)

! NOTE: These are global, must be set
g2ilo = glo(1)
g2jlo = glo(2)
g2klo = glo(3)

! Check that mesh sizes haven't changed
if (allocated(cgrn2)) then
  good_sizes = .true. 
  do i=1, 3
    if (lbound(cgrn2, i) /= glo(i)) good_sizes = .false.
    if (ubound(cgrn2, i) /= ghi(i)) good_sizes = .false.
  enddo
  
  if (.not. good_sizes) then
    !! print *, 'green function array size changed, deallocating'
    deallocate(cgrn2)
  endif
endif

if(.not. allocated(cgrn2)) allocate(cgrn2( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)) )

end subroutine osc_alloc_image_array

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine osc_getgrnimageshift(cgrn,delta,gam,zshft,ilo_grn,jlo_grn,klo_grn,npad,icomp,igfflag)
implicit none
real(dp), dimension(3) :: delta
integer, dimension(3) :: npad
real(dp) :: dx,dy,dz,gam,zshft
integer :: ilo_grn,jlo_grn,klo_grn,igfflag,icomp
integer :: ilo_grn_gbl,jlo_grn_gbl,klo_grn_gbl !in a parallel code these would be in the arg list
integer :: ipad,jpad,kpad,ishift,jshift,kshift,iperiod,jperiod,kperiod
complex(dp), dimension(ilo_grn:,jlo_grn:,klo_grn:) :: cgrn
integer :: i,j,k,ip,jp,kp
real(dp) :: u,v,w
real(dp) :: gval
dx=delta(1); dy=delta(2); dz=delta(3)
ilo_grn_gbl=ilo_grn; jlo_grn_gbl=jlo_grn; klo_grn_gbl=klo_grn
iperiod=size(cgrn,1); jperiod=size(cgrn,2); kperiod=size(cgrn,3)
ipad=npad(1); jpad=npad(2); kpad=npad(3)
!this puts the Green function where it's needed so the convolution ends up in the correct location in the array
ishift=iperiod/2-(ipad+1)/2; jshift=jperiod/2-(jpad+1)/2; kshift=kperiod/2-(kpad+1)/2
!$ print *, 'OpenMP Green function calc'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k=klo_grn,ubound(cgrn,3)
  kp=klo_grn_gbl+mod(k-klo_grn_gbl+kshift ,kperiod) !in the serial code klo_grn_gbl = klo_grn, similarly for i,j
  w=kp*dz+zshft
  do j=jlo_grn,ubound(cgrn,2)
    jp=jlo_grn_gbl+mod(j-jlo_grn_gbl+jshift ,jperiod)
    v=jp*dy
   do i=ilo_grn,ubound(cgrn,1)
     ip=ilo_grn_gbl+mod(i-ilo_grn_gbl+ishift ,iperiod)
     u=ip*dx
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn(i,j,k)= cmplx(gval,0.d0, dp)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call ccfft3d(cgrn,cgrn,(/1,1,1/),iperiod,jperiod,kperiod,0)

end subroutine osc_getgrnimageshift


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine osc_getgrnimageconvcorr(cgrn,umin,delta,gam,ilo_grn,jlo_grn,klo_grn,nlo,npad,icomp,igfflag)
implicit none
real(dp), dimension(3) :: umin,delta
integer, dimension(3) :: nlo,npad
real(dp) :: dx,dy,dz,gam
integer :: ilo_grn,jlo_grn,klo_grn,igfflag,icomp
integer :: ilo_grn_gbl,jlo_grn_gbl,klo_grn_gbl !in a parallel code these would be in the arg list
integer :: ipad,jpad,kpad,ishift,jshift,kshift,iperiod,jperiod,kperiod
complex(dp), dimension(ilo_grn:,jlo_grn:,klo_grn:) :: cgrn
integer :: i,j,k,ip,jp,kp
real(dp) :: u,v,w
real(dp) :: gval
dx=delta(1); dy=delta(2); dz=delta(3)
ilo_grn_gbl=ilo_grn; jlo_grn_gbl=jlo_grn; klo_grn_gbl=klo_grn
iperiod=size(cgrn,1); jperiod=size(cgrn,2); kperiod=size(cgrn,3)
ipad=npad(1); jpad=npad(2); kpad=npad(3)
!this puts the Green function where it's needed so the convolution ends up in the correct location in the array
ishift=iperiod/2-(ipad+1)/2; jshift=jperiod/2-(jpad+1)/2; kshift=kperiod/2-(kpad+1)/2
!$ print *, 'OpenMP Green function calc'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k=klo_grn,ubound(cgrn,3)
! kp=klo_grn_gbl+mod(k-klo_grn_gbl+kshift ,kperiod)
! w=kp*dz
  kp=k
  w=kp*dz+2.d0*umin(3)-2.d0*nlo(3)*delta(3)
  do j=jlo_grn,ubound(cgrn,2)
    jp=jlo_grn_gbl+mod(j-jlo_grn_gbl+jshift ,jperiod)
    v=jp*dy
   do i=ilo_grn,ubound(cgrn,1)
     ip=ilo_grn_gbl+mod(i-ilo_grn_gbl+ishift ,iperiod)
     u=ip*dx
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn(i,j,k)= cmplx(gval,0.d0, dp)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call ccfft3d(cgrn,cgrn,(/1,1,-1/),iperiod,jperiod,kperiod,0)

end subroutine osc_getgrnimageconvcorr

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine imageconvcorr3d(crho,qgrn1,con,rilo,rjlo,rklo,g1ilo,g1jlo,g1klo,cilo,cjlo,cklo,iperiod,jperiod,kperiod)
implicit none
integer :: rilo,rjlo,rklo,g1ilo,g1jlo,g1klo,cilo,cjlo,cklo,iperiod,jperiod,kperiod
!input arrays:
complex(dp), dimension(rilo:,rjlo:,rklo:) :: crho
complex(dp), dimension(g1ilo:,g1jlo:,g1klo:) :: qgrn1
real(dp), dimension(cilo:,cjlo:,cklo:) :: con
! local:
complex(dp), allocatable, dimension(:,:,:) :: ccon
real(dp) :: fpei,qtot,factr
integer :: cihi,cjhi,ckhi

print *, 'convolution/correlation'
fpei=299792458.d0**2*1.d-7  ! this is 1/(4 pi eps0)
qtot=1.d0 !fix later: 1.d-9 ! 1 nC
allocate(ccon(cilo:cilo+iperiod-1,cjlo:cjlo+jperiod-1,cklo:cklo+kperiod-1))
ccon(:,:,:)=crho(:,:,:)*qgrn1(:,:,:)
call ccfft3d(ccon,ccon,(/-1,-1,+1/),iperiod,jperiod,kperiod,0)
!normalize:
!     factr=hx*hy*hz/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot
factr=    1.d0/( (1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod) )*fpei*qtot
!store final result in original size (not double size) real array:
cihi=ubound(con,1)
cjhi=ubound(con,2)
ckhi=ubound(con,3)
con(cilo:cihi,cjlo:cjhi,cklo:ckhi)=factr*real(ccon(cilo:cihi,cjlo:cjhi,cklo:ckhi),dp)

end subroutine imageconvcorr3d

end module open_spacecharge_core_mod
