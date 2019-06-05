module open_spacecharge_mod

use, intrinsic :: iso_fortran_env
use open_spacecharge_core_mod

implicit none


! Fortran 2008
!integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
!integer, parameter, private :: qp = REAL128

type mesh3d_struct
  integer :: nlo(3) = [ 1,  1,  1]       ! Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: nhi(3) = [64, 64, 64]       ! Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: npad(3) = [ 1,  1,  1]      ! Array padding for cyclic convolution
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma
  real(dp) :: charge                    ! Total charge on mesh
  real(dp), allocatable, dimension(:,:,:) :: rho        ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: phi        ! electric potential grid
  real(dp), allocatable, dimension(:,:,:,:) :: efield   ! electric field grid
  real(dp), allocatable, dimension(:,:,:,:) :: bfield   ! magnetic field grid
end type

contains


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine print_mesh3d(mesh3d)
type(mesh3d_struct) :: mesh3d
print *, '------------------------'
print *, 'Mesh: '
print *, 'nlo: ', mesh3d%nlo
print *, 'nhi: ', mesh3d%nhi
print *, 'min: ', mesh3d%min
print *, 'max: ', mesh3d%max
print *, 'delta: ', mesh3d%delta
print *, 'gamma: ', mesh3d%gamma
print *, 'charge: ', mesh3d%charge
if (allocated(mesh3d%rho)) print *, 'rho allocated'
if (allocated(mesh3d%phi)) print *, 'phi allocated'
if (allocated(mesh3d%efield)) print *, 'efield allocated'
if (allocated(mesh3d%bfield)) print *, 'bfield allocated'
print *, '------------------------'
end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine space_charge_freespace(mesh3d, direct_field_calc, integrated_green_function)
type(mesh3d_struct) :: mesh3d
integer :: idirectfieldcalc=1
integer :: igfflag=1
logical, optional :: direct_field_calc, integrated_green_function

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0

call osc_freespace_solver(mesh3d%rho, mesh3d%gamma, &
  mesh3d%delta(1), mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, mesh3d%npad, idirectfieldcalc,igfflag)
end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine space_charge_rectpipe(mesh3d, apipe, bpipe, direct_field_calc, integrated_green_function)
type(mesh3d_struct) :: mesh3d
real(dp) :: apipe, bpipe
integer :: idirectfieldcalc=1
integer :: igfflag=1
logical, optional :: direct_field_calc, integrated_green_function

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0

call osc_rectpipe_solver(mesh3d%rho, apipe, bpipe, mesh3d%gamma, &
  mesh3d%delta(1), mesh3d%min, mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, idirectfieldcalc,igfflag)

end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine deposit_particles(xa, ya, za, mesh3d, total_charge, qa, resize_mesh)
!
! Deposits particle arrays onto mesh
!
! Input:
!   xa           -- REAL64: x coordinate array
!   ya           -- REAL64: y coordinate array
!   za           -- REAL64: z coordinate array
!   qa           -- REAL64, optional: charge coordinate array
!   total_charge -- REAL64, optional: total charge of particles, used only if qa is not present
!   resize_mesh  -- logical, optional : Set mesh bounds to fit bunch. 
!                            default  : .true.
!
! Output:
!   mesh3d      -- mesh3d_struct:   
!                     %rho(:,:,:) 
!                     %charge
!
!-

!routine for charge deposition
subroutine deposit_particles(xa, ya, za, mesh3d, qa, total_charge, resize_mesh)

real(dp) :: xa(:),ya(:),za(:)
type(mesh3d_struct) :: mesh3d
real(dp), optional :: qa(:), total_charge
logical, optional :: resize_mesh
logical :: resize

real(dp) :: dx,dy,dz,xmin,ymin,zmin, xmax,ymax,zmax, charge1, min(3), max(3), delta(3), pad(3)
real(dp) :: dxi,dyi,dzi,ab,de,gh
integer :: ilo, jlo, klo, ihi, jhi, khi, nlo(3), nhi(3)
integer :: n,ip,jp,kp, n_particles

if (present(resize_mesh)) then
    resize = resize_mesh
else
    ! Default
    resize = .true.
endif

if (resize) then
  min = [minval(xa), minval(ya), minval(za)]
  max = [maxval(xa), maxval(ya), maxval(za)] 
  delta =(max(:) - min(:) ) / (mesh3d%nhi(:) - mesh3d%nlo(:) + 1)
  ! Pad by by 1.1 bins
  pad = delta*1.1
  min = min !- pad
  max = max + pad
  delta =(max(:) - min(:) ) / (mesh3d%nhi(:) - mesh3d%nlo(:) + 1)

  mesh3d%min = min
  mesh3d%max = max
  mesh3d%delta = delta
endif

if (.not. allocated(mesh3d%rho)) then
  nlo = mesh3d%nlo
  nhi = mesh3d%nhi  
  allocate(mesh3d%rho(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3)))
  allocate(mesh3d%phi(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3)))
  allocate(mesh3d%efield(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3), 3))
  allocate(mesh3d%bfield(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3), 3))
endif

dx = mesh3d%delta(1)
dy = mesh3d%delta(2)
dz = mesh3d%delta(3)

dxi=1.d0/dx
dyi=1.d0/dy
dzi=1.d0/dz

ilo = mesh3d%nlo(1)
jlo = mesh3d%nlo(2)
klo = mesh3d%nlo(3)
ihi = mesh3d%nhi(1)
jhi = mesh3d%nhi(2)
khi = mesh3d%nhi(3)

xmin = mesh3d%min(1)
ymin = mesh3d%min(2)
zmin = mesh3d%min(3)

n_particles = size(xa)



if (present(qa)) then
    mesh3d%charge = sum(qa)
elseif (present(total_charge) .and. .not. present(qa)) then
    mesh3d%charge = total_charge
    charge1 = total_charge/n_particles
else
    write(*, *) 'Error, must specify qa or total_charge'
    return
endif

! Clear 
mesh3d%rho(ilo:ihi,jlo:jhi,klo:khi) = 0

do n=1, n_particles
  !if(lostflag(n).ne.0.d0)cycle
  ip=floor((xa(n)-xmin)*dxi+1) !this is 1-based; use ilogbl for the general case
  jp=floor((ya(n)-ymin)*dyi+1)
  kp=floor((za(n)-zmin)*dzi+1)
  ab=((xmin-xa(n))+ip*dx)*dxi
  de=((ymin-ya(n))+jp*dy)*dyi
  gh=((zmin-za(n))+kp*dz)*dzi
! this "if" statement slows things down, but I'm using it for debugging purposes

!  if(ip<ilo .or. jp<jlo .or. kp<klo .or. ip>ihi .or. jp>jhi .or. kp>khi)then
!    ifail=ifail+1
!    cycle
!  else

  if (present(qa)) charge1 = qa(n)

  mesh3d%rho(ip,jp,kp)=mesh3d%rho(ip,jp,kp)+ab*de*gh*charge1
  mesh3d%rho(ip,jp+1,kp)=mesh3d%rho(ip,jp+1,kp)+ab*(1.-de)*gh*charge1
  mesh3d%rho(ip,jp+1,kp+1)=mesh3d%rho(ip,jp+1,kp+1)+ab*(1.-de)*(1.-gh)*charge1
  mesh3d%rho(ip,jp,kp+1)=mesh3d%rho(ip,jp,kp+1)+ab*de*(1.-gh)*charge1
  mesh3d%rho(ip+1,jp,kp+1)=mesh3d%rho(ip+1,jp,kp+1)+(1.-ab)*de*(1.-gh)*charge1
  mesh3d%rho(ip+1,jp+1,kp+1)=mesh3d%rho(ip+1,jp+1,kp+1)+(1.-ab)*(1.-de)*(1.-gh)*charge1
  mesh3d%rho(ip+1,jp+1,kp)=mesh3d%rho(ip+1,jp+1,kp)+(1.-ab)*(1.-de)*gh*charge1
  mesh3d%rho(ip+1,jp,kp)=mesh3d%rho(ip+1,jp,kp)+(1.-ab)*de*gh*charge1
!  endif
enddo

!if(ifail.ne.0)write(6,*)'(depose_rho_scalar) ifail=',ifail


end subroutine deposit_particles




!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine interpolate_field(x, y, z, mesh3d, E, B)
!
! Interpolate field on mesh
!
! Input:
!   x, y, z   -- REAL64: coordinates to interpolate
!   mesh3d    -- mesh3d_struct:  contains efield, bfield
!
! Output:
!   E(3)      -- REAL64, optional: interpolated electric field at x, y, z   
!   B(3)      -- REAL64, optional: interpolated magnetic field at x, y, z   
!
!-

subroutine interpolate_field(x, y, z, mesh3d, E, B)
type(mesh3d_struct) ::  mesh3d
real(dp) :: x, y, z
real(dp), optional :: E(3), B(3)
real(dp) :: hxi,hyi,hzi,ab,de,gh
integer :: ip,jp,kp,ip1,jp1,kp1
integer :: nflag
nflag=0
hxi=1.d0/mesh3d%delta(1); hyi=1.d0/mesh3d%delta(2); hzi=1.d0/mesh3d%delta(3)      
ip=floor((x-mesh3d%min(1))*hxi+1)
jp=floor((y-mesh3d%min(2))*hyi+1)
kp=floor((z-mesh3d%min(3))*hzi+1)

if(ip<1 .or. ip>mesh3d%nhi(1)-1)then
  nflag=1
  write(6,*)'ierror: ip=', ip, (x-mesh3d%min(1)), (x-mesh3d%min(1))/mesh3d%delta(1)
  if(ip<1)then
    ip=1
  else
    ip=mesh3d%nhi(1)-1
  endif
endif
if(jp<1 .or. jp>mesh3d%nhi(2)-1)then
  nflag=1
  write(6,*)'jerror: jp=', jp, (y-mesh3d%min(2)), (y-mesh3d%min(2))/mesh3d%delta(2)
  write(6,*)ab,de,gh
  if(jp<1)then
    jp=1
  else
    jp=mesh3d%nhi(2)-1
  endif
endif
  
if(kp<1 .or. kp>mesh3d%nhi(3)-1)then
    nflag=1
    write(6,*)'kerror:  kp=',kp
!!!!!!!!!!write(6,*)ab,de,gh
  if(kp<1)then
    kp=1
  else
    kp=mesh3d%nhi(3)-1
  endif
endif
ab=((mesh3d%min(1)-x)+ip*mesh3d%delta(1))*hxi
de=((mesh3d%min(2)-y)+jp*mesh3d%delta(2))*hyi
gh=((mesh3d%min(3)-z)+kp*mesh3d%delta(3))*hzi
if(nflag.eq.1)then
  write(6,*)ab,de,gh
  nflag=0
endif

ip1=ip+1
jp1=jp+1
kp1=kp+1

if (present(E)) then 
  E=mesh3d%efield(ip, jp,  kp,  :)*ab*de*gh                 &
   +mesh3d%efield(ip, jp1, kp,  :)*ab*(1.-de)*gh            &
   +mesh3d%efield(ip, jp1, kp1, :)*ab*(1.-de)*(1.-gh)       &
   +mesh3d%efield(ip, jp,  kp1, :)*ab*de*(1.-gh)            &
   +mesh3d%efield(ip1,jp,  kp1, :)*(1.-ab)*de*(1.-gh)       &
   +mesh3d%efield(ip1,jp1, kp1, :)*(1.-ab)*(1.-de)*(1.-gh)  &
   +mesh3d%efield(ip1,jp1, kp,  :)*(1.-ab)*(1.-de)*gh       &
   +mesh3d%efield(ip1,jp,  kp,  :)*(1.-ab)*de*gh
endif

if (present(B)) then 
  B=mesh3d%bfield(ip, jp,  kp,  :)*ab*de*gh                 &
   +mesh3d%bfield(ip, jp1, kp,  :)*ab*(1.-de)*gh            &
   +mesh3d%bfield(ip, jp1, kp1, :)*ab*(1.-de)*(1.-gh)       &
   +mesh3d%bfield(ip, jp,  kp1, :)*ab*de*(1.-gh)            &
   +mesh3d%bfield(ip1,jp,  kp1, :)*(1.-ab)*de*(1.-gh)       &
   +mesh3d%bfield(ip1,jp1, kp1, :)*(1.-ab)*(1.-de)*(1.-gh)  &
   +mesh3d%bfield(ip1,jp1, kp,  :)*(1.-ab)*(1.-de)*gh       &
   +mesh3d%bfield(ip1,jp,  kp,  :)*(1.-ab)*de*gh
endif

end subroutine


end module
