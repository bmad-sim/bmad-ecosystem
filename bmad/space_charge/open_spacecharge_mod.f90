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


!+
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

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


!+
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

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
  mesh3d%delta, mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, mesh3d%npad, idirectfieldcalc,igfflag)
end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine space_charge_cathodeimages(mesh3d, direct_field_calc, integrated_green_function, image_method)
type(mesh3d_struct) :: mesh3d
integer :: idirectfieldcalc=1
integer :: igfflag=1
logical, optional :: direct_field_calc, integrated_green_function
integer, optional :: image_method
integer :: imethod

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0

imethod=1 ! =1 for convolution/correlation, =2 for shift method
if(present(image_method) .and. image_method.eq.2)imethod=2

call osc_cathodeimages_solver(mesh3d%rho, mesh3d%gamma, &
  mesh3d%min,mesh3d%delta, mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, mesh3d%npad, idirectfieldcalc,igfflag,imethod)
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
logical :: resize, good_mesh_sizes

real(dp) :: dx,dy,dz,xmin,ymin,zmin, xmax,ymax,zmax, charge1, min(3), max(3), delta(3), pad(3)
real(dp) :: dxi,dyi,dzi,ab,de,gh
integer :: i, ilo, jlo, klo, ihi, jhi, khi, nlo(3), nhi(3)
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
  delta =(max(:) - min(:) ) / (mesh3d%nhi(:) - mesh3d%nlo(:) )
  
  ! Small padding to protect against indexing errors
  min = min - 1.0e-6_dp*delta
  max = max + 1.0e-6_dp*delta
  delta =(max(:) - min(:) ) / (mesh3d%nhi(:) - mesh3d%nlo(:) )
  
  mesh3d%min = min
  mesh3d%max = max
  mesh3d%delta = delta
endif


! If allocated, check that the mesh sizes agree. Otherwise deallocate and below will allocate. 
if (allocated(mesh3d%rho)) then
  nlo = mesh3d%nlo
  nhi = mesh3d%nhi  
  good_mesh_sizes= .true. 
  do i=1, 3
    if (lbound(mesh3d%rho, i) /= nlo(i) ) good_mesh_sizes = .false.
    if (ubound(mesh3d%rho, i) /= nhi(i) ) good_mesh_sizes = .false.
  enddo
  
  if (.not. good_mesh_sizes) then
   ! print *, 'mesh size changed, deallocating'
    deallocate(mesh3d%rho)
    deallocate(mesh3d%phi)
    deallocate(mesh3d%efield)
    deallocate(mesh3d%bfield)
  endif
endif

if (.not. allocated(mesh3d%rho)) then
  nlo = mesh3d%nlo
  nhi = mesh3d%nhi 
  !print *, 'Allocating new meshes, nlo, nhi: ', nlo, nhi
  allocate(mesh3d%rho(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3)))
  allocate(mesh3d%phi(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3)))
  allocate(mesh3d%efield(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3), 3))
  allocate(mesh3d%bfield(nlo(1):nhi(1),nlo(2):nhi(2), nlo(3):nhi(3), 3))
endif


if (mesh3d%delta(1) == 0) mesh3d%delta(1) = 1d-10
if (mesh3d%delta(2) == 0) mesh3d%delta(2) = 1d-10
if (mesh3d%delta(3) == 0) mesh3d%delta(3) = 1d-10

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
    write(6,*)'z=',z
    write(6,*)'mesh3d%min(3)=',mesh3d%min(3)
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


!------------------------------------------------------------------------
!+
! Subroutine space_charge_3d(mesh3d, offset, at_cathode, calc_bfield, image_efield)
!
! Performs the space charge calculation using the integrated Green function method
! and FFT-based convolutions. 
!
! Input:
!   mesh3d        -- mesh3d_struct: populated with %rho
!
!   offset        -- real(3), optional: Offset coordinates x0, y0, z0 to evaluate the field,
!                    relative to rho. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to rho. 
!
!   at_cathode    -- logical, optional: Maintain constant voltage at the cathode 
!                       using image charges. Default is False. 
!
!   calc_bfield   -- logical, optional: Calculate the magnetic field mesh3d%bfield
!
!                     Default: False
!
!   image_efield(:,:,:,:) -- real(dp), allocatable, optional: Must be present when at_cathode if true.
!                       
!
! Output:
!   mesh3d        -- mesh3d_struct: populated with %efield, and optionally %bfield                             
!
!
!
! Notes: 
!   The magnetic field components can be calculated by:
!     Bx = -(beta/c) * Ey
!     By =  (beta/c) * Ex
!     Bz = 0
!   The image charges move in the opposite direction, so the signs are flipped. 
!     
!
!-
subroutine space_charge_3d(mesh3d, offset, at_cathode, calc_bfield, image_efield)
type(mesh3d_struct) :: mesh3d
real(dp), allocatable, dimension(:,:,:,:), optional :: image_efield   ! electric field grid
real(dp), optional :: offset(3)
real(dp) :: offset0(3), beta
real(dp), parameter :: c_light = 299792458.0
logical, optional :: at_cathode, calc_bfield
logical :: bcalc

if (present(calc_bfield)) then
  bcalc = calc_bfield
else
  bcalc = .false.
endif

if (.not. present(offset)) then
  offset0 = 0
else
  offset0 = offset
endif

! Free space field
call osc_freespace_solver2(mesh3d%rho, mesh3d%gamma, mesh3d%delta, efield=mesh3d%efield, offset=offset0)

! Optional B field
if (bcalc) then
  beta = sqrt(1-1/mesh3d%gamma**2)
  mesh3d%bfield = 0
  mesh3d%bfield(:,:,:,1) = -(beta/c_light)*mesh3d%efield(:,:,:,2)
  mesh3d%bfield(:,:,:,2) =  (beta/c_light)*mesh3d%efield(:,:,:,1)
endif

! Cathode calc
if (.not. present(at_cathode)) return
if (.not. at_cathode) return

! Allocate scratch array for the image field
if (present(image_efield)) then
  if (.not. allocated(image_efield)) allocate(image_efield, mold=mesh3d%efield)
endif

! Image field, with an offset assuming the cathode is at z=0.
! The offset is the z width of the mesh, plus 2 times the distance of the mesh from the cathode.
offset0(3) = offset0(3) + 2*mesh3d%min(3) + (mesh3d%max(3)-mesh3d%min(3))

! Flip the charge mesh in z, with opposite charge sign
call osc_freespace_solver2(-mesh3d%rho(:,:,size(mesh3d%rho,3):1:-1), &
  mesh3d%gamma, mesh3d%delta, efield=image_efield, offset=offset0)  
  
  
! Finally add fields
if (bcalc) then
  ! Opposite sign for beta, because image charges are moving in the negative z direction
  mesh3d%bfield(:,:,:,1) = mesh3d%bfield(:,:,:,1) + (beta/c_light)*image_efield(:,:,:,2)
  mesh3d%bfield(:,:,:,2) = mesh3d%bfield(:,:,:,2) - (beta/c_light)*image_efield(:,:,:,1)
endif  
    
mesh3d%efield = mesh3d%efield + image_efield

end subroutine



!------------------------------------------------------------------------
!+
! Subroutine osc_freespace_solver2(rho, gamma, delta, efield, phi, offset)
!
! Deposits particle arrays onto mesh
!
! Input:
!   rho          -- REAL64(:,:,:): charge density array in x, y, z
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!   gamma        -- REAL64: relativistic gamma
!   icomp        -- integer: Field component requested:
!                        0: phi (scalar potential)
!                       
!
!   efield        -- REAL64(:,:,:,3), optional: allocated electric field array to populate.
!                      
!                                     The final index corresponds to components
!                                     1: Ex
!                                     2: Ey
!                                     3: Ez                                   
!                                     If present, all components will be computed.    
!
!   phi           -- REAL64(:,:,:), optional: allocated potential array to populate
!
!   offset        -- real(3), optional: Offset coordinates x0, y0, z0 to evaluate the field,
!                    relative to rho. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to rho. 
!
! Output:
!   efield        -- REAL64(:,:,:,:) : electric field                                 
!   phi           -- REAL64(:,:,:)   : potential
!
!
! Notes: 
!   The magnetic field components can be calculated by:
!     Bx = -(beta/c) * Ey
!     By =  (beta/c) * Ex
!     Bz = 0
!
!-
subroutine osc_freespace_solver2(rho, gamma, delta, efield, phi, offset)


real(dp), intent(in), dimension(:,:,:) :: rho
real(dp), intent(in) :: gamma, delta(3)
real(dp), optional, intent(out), dimension(:,:,:,:) :: efield
real(dp), optional, intent(out), dimension(:,:,:) :: phi
real(dp), intent(in), optional :: offset(3)
! internal arrays
complex(dp), allocatable, dimension(:,:,:) :: crho, cgrn
real(dp) :: factr, offset0=0
real(dp), parameter :: clight=299792458.0
real(dp), parameter :: fpei=299792458.0**2*1.00000000055d-7  ! this is 1/(4 pi eps0) after the 2019 SI changes

integer :: nx, ny, nz, nx2, ny2, nz2
integer :: icomp, ishift, jshift, kshift

! Sizes
nx = size(rho, 1); ny = size(rho, 2); nz = size(rho, 3)
nx2 = 2*nx; ny2 = 2*ny; nz2 = 2*nz; 

! Allocate complex scratch arrays
allocate(crho(nx2, ny2, nz2))
allocate(cgrn(nx2, ny2, nz2))

! rho -> crho -> FFT(crho)
crho = 0
crho(1:nx, 1:ny, 1:nz) = rho ! Place in one octant
call ccfft3d(crho, crho, [1,1,1], nx2, ny2, nz2, 0) 

! Loop over phi, Ex, Ey, Ez
do icomp=0, 3
  if ((icomp == 0) .and. (.not. present(phi))) cycle
  if ((icomp == 1) .and. (.not. present(efield))) exit

  call osc_get_cgrn_freespace(cgrn, delta, gamma, icomp, offset=offset)
  
  !  cgrn -> FFT(cgrn)
  call ccfft3d(cgrn, cgrn, [1,1,1], nx2, ny2, nz2, 0)  
  
  ! Multiply FFT'd arrays, re-use cgrn
  cgrn=crho*cgrn

  ! Inverse FFT
  call ccfft3d(cgrn, cgrn, [-1,-1,-1], nx2, ny2, nz2, 0)  
  
  ! This is where the output is shifted to
  ishift = nx-1
  jshift = ny-1
  kshift = nz-1
  
  ! Extract field
  factr = fpei/(nx2*ny2*nz2)
  
  if (icomp == 0) then
    phi(:,:,:) = factr * real(cgrn(1+ishift:nx+ishift, 1+jshift:ny+jshift, 1+kshift:nz+kshift), dp)  
  else
    efield(:,:,:,icomp) = factr * real(cgrn(1+ishift:nx+ishift, 1+jshift:ny+jshift, 1+kshift:nz+kshift), dp)
  endif
    
enddo

end subroutine osc_freespace_solver2


!------------------------------------------------------------------------
!+
! Subroutine osc_get_cgrn_freespace(cgrn, delta, gamma, icomp, offset)
!
! Computes the free space Green function on a mesh with given spacings in the lab frame.
! The computation is performed in the rest fram by boosting the coordinates by gamma.
!
!
! Input:
!   cgrn         -- COMPLEX128(:,:,:): pre-allocated array 
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!   gamma        -- REAL64: relativistic gamma
!   icomp        -- integer: Field component requested:
!                        0: phi (scalar potential)
!                        1: Ex
!                        2: Ey
!                        3: Ez
!   offset        -- real(3), optional: Offset coordinates for the center of the grid in [m]. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to the rho_mesh center. 
!                              
! Output:
!   cgrn         -- COMPLEX128(:,:,:): Green function array
!   
!                
! Notes:
!   Internally, dz -> dz*gamma.             
!   For efficients, the indefinite functions lafun2, xlafun2 are actually evaluated 
!   on a grid slightly offset by -dx/2, -dy/2, -dz/2,
!   and these points are used to evaluate the integral with array addition and subtraction. 
!
!-
subroutine osc_get_cgrn_freespace(cgrn, delta, gamma, icomp, offset)

complex(dp), intent(out), dimension(:,:,:) :: cgrn 
real(dp), intent(in), dimension(3) :: delta
integer, intent(in) :: icomp
real(dp), intent(in) :: gamma
real(dp), intent(in), optional :: offset(3)
! Local
real(dp) :: dx,dy,dz
real(dp) :: u,v,w, umin, vmin, wmin
real(dp) :: gval, factor
integer :: imin, imax, jmin, jmax, kmin, kmax
integer :: i,j,k, isize, jsize, ksize

! Mesh spacings. dz is stretched in the rest frame.
dx=delta(1); dy=delta(2); dz=delta(3)*gamma

! Special cases: Ex and Ey are enhanced by gamma
if ((icomp==1) .or. (icomp==2)) then
  factor = gamma /(dx*dy*dz)
else
  factor = 1.0 /(dx*dy*dz)
endif

! Evaluate on an offset grid, for use in the indefinite integral evaluation below. 
isize = size(cgrn,1); jsize=size(cgrn,2); ksize=size(cgrn,3)
umin = (0.5-isize/2) *dx
vmin = (0.5-jsize/2) *dy
wmin = (0.5-ksize/2) *dz

! Add optional offset
if (present(offset)) then
  umin = umin + offset(1)
  vmin = vmin + offset(2) 
  wmin = wmin + offset(3)*gamma ! Don't forget this!
endif

! !$ print *, 'OpenMP Green function calc osc_get_cgrn_freespace'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k = 1, ksize
  w = (k-1)*dz + wmin

  do j=1, jsize
    v=(j-1)*dy + vmin

   do i=1, isize
     u = (i-1)*dx + umin
     
     if(icomp == 0) gval=lafun2(u,v,w)*factor
     if(icomp == 1) gval=xlafun2(u,v,w)*factor
     if(icomp == 2) gval=xlafun2(v,w,u)*factor
     if(icomp == 3) gval=xlafun2(w,u,v)*factor
     cgrn(i,j,k)= cmplx(gval, 0, dp)
     
    enddo
  enddo
enddo
!$OMP END PARALLEL DO

! Evaluate the indefinite integral over the cube     
!cgrn = cgrn(1:,1:,1:) - cgrn(:-1,1:,1:) - cgrn(1:,:-1,1:) - cgrn(1:,1:,:-1) - cgrn(:-1,:-1,:-1) + cgrn(:-1,:-1,1:) + cgrn(:-1,1:,:-1)  + cgrn(1:,:-1,:-1)
!       (x2,y2,z2) -      (x1,y2,z2) -     (x2,y1,z2) -    (x2,y2,z1) -    (x1,y1,z1)   +    (x1,y1,z2)   +    (x1,y2,z1)  +    (x2,y1,z1)
cgrn(1:isize-1, 1:jsize-1, 1:ksize-1) = &
    + cgrn(2:isize,   2:jsize,   2:ksize) &
    - cgrn(1:isize-1, 2:jsize,   2:ksize) &
    - cgrn(2:isize,   1:jsize-1, 2:ksize) &
    - cgrn(2:isize,   2:jsize,   1:ksize-1) &
    - cgrn(1:isize-1, 1:jsize-1, 1:ksize-1) &
    + cgrn(1:isize-1, 1:jsize-1, 2:ksize) &
    + cgrn(1:isize-1, 2:jsize,   1:ksize-1)  &
    + cgrn(2:isize,   1:jsize-1, 1:ksize-1)

end subroutine osc_get_cgrn_freespace


!------------------------------------------------------------------------
!+
! elemental real(dp) function xlafun2(x, y, z)
!
! The indefinite integral:
! \int x/r^3 dx dy dz = x*atan((y*z)/(r*x)) -z*log(r+y) + y*log((r-z)/(r+z))/2
!
! This corresponds to the electric field component Ex.
! Other components can be computed by permuting the arguments
!
!-
elemental real(dp) function xlafun2(x, y, z)
  real(dp), intent(in) :: x, y, z
  real(dp) :: r
  r=sqrt(x**2+y**2+z**2)
  xlafun2 = x*atan((y*z)/(r*x)) -z*log(r+y) + y*log((r-z)/(r+z))/2
end function

!------------------------------------------------------------------------
!+
! elemental real(dp) function lafun2(x,y,z)
!
! The indefinite integral:
! \int 1/r dx dy dz = 
!          -z**2*atan(x*y/(z*r))/2 - y**2*atan(x*z/(y*r))/2 -x**2*atan(y*z/(x*r))/2 
!           +y*z*log(x+r) + x*z*log(y+r) + x*y*log(z+r)
!
! This corresponds to the scalar potential.
! Other components can be computed by permuting the arguments
!
!-
elemental real(dp) function lafun2(x,y,z)
  real(dp), intent(in) :: x, y, z
  real(dp) :: r
  r=sqrt(x**2+y**2+z**2)
  lafun2 = -z**2*atan(x*y/(z*r))/2 - y**2*atan(x*z/(y*r))/2 -x**2*atan(y*z/(x*r))/2 &
           +y*z*log(x+r) + x*z*log(y+r) + x*y*log(z+r)
end function 


end module
