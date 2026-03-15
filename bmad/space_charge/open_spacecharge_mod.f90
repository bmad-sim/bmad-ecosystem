module open_spacecharge_mod

use, intrinsic :: iso_fortran_env
use open_spacecharge_core_mod

implicit none


! Fortran 2008
!integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
!integer, parameter, private :: qp = REAL128

! Saved mesh geometry for lazy resizing across calls to deposit_particles.
! Persists between calls so delta stays constant when the bunch barely changes.
real(dp), save, private :: saved_mesh_min(3) = 0
real(dp), save, private :: saved_mesh_max(3) = 0
logical, save, private :: saved_mesh_valid = .false.

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
! Subroutine deposit_particles(xa, ya, za, mesh3d, total_charge, qa, resize_mesh, mesh_growth_factor, mesh_shrink_factor)
!
! Deposits particle arrays onto mesh
!
! Input:
!   xa                  -- REAL64: x coordinate array
!   ya                  -- REAL64: y coordinate array
!   za                  -- REAL64: z coordinate array
!   qa                  -- REAL64, optional: charge coordinate array
!   total_charge        -- REAL64, optional: total charge of particles, used only if qa is not present
!   resize_mesh         -- logical, optional : Set mesh bounds to fit bunch. 
!                                   default  : .true.
!   mesh_growth_factor  -- REAL64, optional: Fractional padding when growing mesh (default: 0 = tight fit).
!   mesh_shrink_factor  -- REAL64, optional: Fractional threshold for shrinking mesh (default: 0 = tight fit).
!
! Output:
!   mesh3d      -- mesh3d_struct:   
!                     %rho(:,:,:) 
!                     %charge
!
!-

!routine for charge deposition
subroutine deposit_particles(xa, ya, za, mesh3d, qa, total_charge, resize_mesh, mesh_growth_factor, mesh_shrink_factor)

real(dp) :: xa(:),ya(:),za(:)
type(mesh3d_struct) :: mesh3d
real(dp), optional :: qa(:), total_charge
logical, optional :: resize_mesh
real(dp), intent(in), optional :: mesh_growth_factor, mesh_shrink_factor
logical :: resize, good_mesh_sizes, needs_resize
real(dp) :: growth_factor, shrink_factor

real(dp) :: dx,dy,dz,xmin,ymin,zmin, xmax,ymax,zmax, charge1, min(3), max(3), delta(3), pad(3)
real(dp) :: dxi,dyi,dzi,ab,de,gh
real(dp) :: oab  ! 1-ab
integer :: i, ilo, jlo, klo, ihi, jhi, khi, nlo(3), nhi(3)
integer :: n,ip,jp,kp, n_particles
real(dp), allocatable :: rho_priv(:,:,:)

if (present(resize_mesh)) then
    resize = resize_mesh
else
    ! Default
    resize = .true.
endif

growth_factor = 0
if (present(mesh_growth_factor)) growth_factor = mesh_growth_factor
shrink_factor = 0
if (present(mesh_shrink_factor)) shrink_factor = mesh_shrink_factor

if (resize) then
  ! Fused min/max: single pass over all coordinate arrays instead of 6 separate minval/maxval calls.
  min(1) = xa(1); max(1) = xa(1)
  min(2) = ya(1); max(2) = ya(1)
  min(3) = za(1); max(3) = za(1)
  do n = 2, size(xa)
    if (xa(n) < min(1)) then; min(1) = xa(n); else if (xa(n) > max(1)) then; max(1) = xa(n); endif
    if (ya(n) < min(2)) then; min(2) = ya(n); else if (ya(n) > max(2)) then; max(2) = ya(n); endif
    if (za(n) < min(3)) then; min(3) = za(n); else if (za(n) > max(3)) then; max(3) = za(n); endif
  enddo

  if (growth_factor == 0 .and. shrink_factor == 0) then
    ! Original tight-fit behavior: compute delta, add tiny relative padding.
    delta = (max - min) / (mesh3d%nhi - mesh3d%nlo)
    min = min - 1.0e-6_dp * delta
    max = max + 1.0e-6_dp * delta
  else
    ! Lazy resize: only change mesh bounds when the bunch no longer fits or is much smaller.
    ! This keeps delta stable across calls, enabling Green function FFT caching.
    needs_resize = .not. saved_mesh_valid
    if (.not. needs_resize) then
      ! Grow: bunch exceeds saved mesh bounds
      if (any(min < saved_mesh_min) .or. any(max > saved_mesh_max)) then
        needs_resize = .true.
      else
        ! Shrink: bunch range significantly smaller than mesh range in any dimension
        delta = saved_mesh_max - saved_mesh_min  ! reuse delta temporarily for mesh_range
        if (any(max - min < delta * (1.0_dp - shrink_factor))) needs_resize = .true.
      endif
    endif

    if (needs_resize) then
      pad = (max - min) * (growth_factor * 0.5_dp)
      ! Ensure minimum padding for numerical safety (handles zero-range dimensions)
      where (pad < 1.0e-10_dp) pad = 1.0e-6_dp
      min = min - pad
      max = max + pad
      saved_mesh_min = min
      saved_mesh_max = max
      saved_mesh_valid = .true.
    else
      min = saved_mesh_min
      max = saved_mesh_max
    endif
  endif

  delta = (max - min) / (mesh3d%nhi - mesh3d%nlo)
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

! Deposit particles onto mesh using trilinear weighting.
if (n_particles >= 50000) then
  ! OpenMP path: each thread accumulates into a private rho copy, then sum into shared mesh.
  !$OMP PARALLEL PRIVATE(rho_priv, n, ip, jp, kp, ab, de, gh, oab, charge1) &
  !$OMP& SHARED(mesh3d, xa, ya, za, qa, n_particles, dxi, dyi, dzi, dx, dy, dz, xmin, ymin, zmin, ilo, ihi, jlo, jhi, klo, khi)
  allocate(rho_priv(ilo:ihi+1, jlo:jhi+1, klo:khi+1))
  rho_priv = 0

  if (present(qa)) then
    !$OMP DO SCHEDULE(STATIC)
    do n = 1, n_particles
      ip = floor((xa(n)-xmin)*dxi+1)
      jp = floor((ya(n)-ymin)*dyi+1)
      kp = floor((za(n)-zmin)*dzi+1)
      ab = ((xmin-xa(n))+ip*dx)*dxi
      de = ((ymin-ya(n))+jp*dy)*dyi
      gh = ((zmin-za(n))+kp*dz)*dzi
      oab = 1.0_dp - ab
      charge1 = qa(n)
      rho_priv(ip,  jp,  kp)   = rho_priv(ip,  jp,  kp)   + ab *de     *gh     *charge1
      rho_priv(ip,  jp+1,kp)   = rho_priv(ip,  jp+1,kp)   + ab *(1-de) *gh     *charge1
      rho_priv(ip,  jp+1,kp+1) = rho_priv(ip,  jp+1,kp+1) + ab *(1-de) *(1-gh) *charge1
      rho_priv(ip,  jp,  kp+1) = rho_priv(ip,  jp,  kp+1) + ab *de     *(1-gh) *charge1
      rho_priv(ip+1,jp,  kp+1) = rho_priv(ip+1,jp,  kp+1) + oab*de     *(1-gh) *charge1
      rho_priv(ip+1,jp+1,kp+1) = rho_priv(ip+1,jp+1,kp+1) + oab*(1-de) *(1-gh) *charge1
      rho_priv(ip+1,jp+1,kp)   = rho_priv(ip+1,jp+1,kp)   + oab*(1-de) *gh     *charge1
      rho_priv(ip+1,jp,  kp)   = rho_priv(ip+1,jp,  kp)   + oab*de     *gh     *charge1
    enddo
    !$OMP END DO
  else
    !$OMP DO SCHEDULE(STATIC)
    do n = 1, n_particles
      ip = floor((xa(n)-xmin)*dxi+1)
      jp = floor((ya(n)-ymin)*dyi+1)
      kp = floor((za(n)-zmin)*dzi+1)
      ab = ((xmin-xa(n))+ip*dx)*dxi
      de = ((ymin-ya(n))+jp*dy)*dyi
      gh = ((zmin-za(n))+kp*dz)*dzi
      oab = 1.0_dp - ab
      rho_priv(ip,  jp,  kp)   = rho_priv(ip,  jp,  kp)   + ab *de     *gh     *charge1
      rho_priv(ip,  jp+1,kp)   = rho_priv(ip,  jp+1,kp)   + ab *(1-de) *gh     *charge1
      rho_priv(ip,  jp+1,kp+1) = rho_priv(ip,  jp+1,kp+1) + ab *(1-de) *(1-gh) *charge1
      rho_priv(ip,  jp,  kp+1) = rho_priv(ip,  jp,  kp+1) + ab *de     *(1-gh) *charge1
      rho_priv(ip+1,jp,  kp+1) = rho_priv(ip+1,jp,  kp+1) + oab*de     *(1-gh) *charge1
      rho_priv(ip+1,jp+1,kp+1) = rho_priv(ip+1,jp+1,kp+1) + oab*(1-de) *(1-gh) *charge1
      rho_priv(ip+1,jp+1,kp)   = rho_priv(ip+1,jp+1,kp)   + oab*(1-de) *gh     *charge1
      rho_priv(ip+1,jp,  kp)   = rho_priv(ip+1,jp,  kp)   + oab*de     *gh     *charge1
    enddo
    !$OMP END DO
  endif

  !$OMP CRITICAL (deposit_rho_reduce)
  mesh3d%rho(ilo:ihi,jlo:jhi,klo:khi) = mesh3d%rho(ilo:ihi,jlo:jhi,klo:khi) + rho_priv(ilo:ihi,jlo:jhi,klo:khi)
  !$OMP END CRITICAL (deposit_rho_reduce)
  deallocate(rho_priv)
  !$OMP END PARALLEL

else
  ! Serial path for small particle counts (avoids OpenMP overhead)
  do n = 1, n_particles
    ip = floor((xa(n)-xmin)*dxi+1)
    jp = floor((ya(n)-ymin)*dyi+1)
    kp = floor((za(n)-zmin)*dzi+1)
    ab = ((xmin-xa(n))+ip*dx)*dxi
    de = ((ymin-ya(n))+jp*dy)*dyi
    gh = ((zmin-za(n))+kp*dz)*dzi
    oab = 1.0_dp - ab
    if (present(qa)) charge1 = qa(n)
    mesh3d%rho(ip,  jp,  kp)   = mesh3d%rho(ip,  jp,  kp)   + ab *de     *gh     *charge1
    mesh3d%rho(ip,  jp+1,kp)   = mesh3d%rho(ip,  jp+1,kp)   + ab *(1-de) *gh     *charge1
    mesh3d%rho(ip,  jp+1,kp+1) = mesh3d%rho(ip,  jp+1,kp+1) + ab *(1-de) *(1-gh) *charge1
    mesh3d%rho(ip,  jp,  kp+1) = mesh3d%rho(ip,  jp,  kp+1) + ab *de     *(1-gh) *charge1
    mesh3d%rho(ip+1,jp,  kp+1) = mesh3d%rho(ip+1,jp,  kp+1) + oab*de     *(1-gh) *charge1
    mesh3d%rho(ip+1,jp+1,kp+1) = mesh3d%rho(ip+1,jp+1,kp+1) + oab*(1-de) *(1-gh) *charge1
    mesh3d%rho(ip+1,jp+1,kp)   = mesh3d%rho(ip+1,jp+1,kp)   + oab*(1-de) *gh     *charge1
    mesh3d%rho(ip+1,jp,  kp)   = mesh3d%rho(ip+1,jp,  kp)   + oab*de     *gh     *charge1
  enddo
endif


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
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine interpolate_field_batch(x, y, z, mesh3d, n_pts, E, B)
!
! Vectorized batch interpolation of fields at multiple points.
! OpenMP parallelized over the particle dimension.
!
! Input:
!   x(n_pts), y(n_pts), z(n_pts) -- REAL64: coordinate arrays
!   mesh3d                        -- mesh3d_struct: contains efield, bfield
!   n_pts                         -- integer: number of points to interpolate
!
! Output:
!   E(3,n_pts) -- REAL64, optional: interpolated electric field
!   B(3,n_pts) -- REAL64, optional: interpolated magnetic field
!-
subroutine interpolate_field_batch(x, y, z, mesh3d, n_pts, E, B)
type(mesh3d_struct), intent(in) :: mesh3d
integer, intent(in) :: n_pts
real(dp), intent(in) :: x(n_pts), y(n_pts), z(n_pts)
real(dp), optional, intent(out) :: E(3, n_pts), B(3, n_pts)

real(dp) :: hxi, hyi, hzi, ab, de, gh, oab, ode, ogh
real(dp) :: w000, w010, w011, w001, w101, w111, w110, w100
integer :: n, ip, jp, kp, ip1, jp1, kp1

hxi = 1.0_dp / mesh3d%delta(1)
hyi = 1.0_dp / mesh3d%delta(2)
hzi = 1.0_dp / mesh3d%delta(3)

!$OMP PARALLEL DO PRIVATE(n, ip, jp, kp, ip1, jp1, kp1, ab, de, gh, oab, ode, ogh, &
!$OMP& w000, w010, w011, w001, w101, w111, w110, w100) &
!$OMP& SHARED(mesh3d, x, y, z, E, B, n_pts, hxi, hyi, hzi) SCHEDULE(STATIC)
do n = 1, n_pts
  ip = floor((x(n) - mesh3d%min(1)) * hxi + 1)
  jp = floor((y(n) - mesh3d%min(2)) * hyi + 1)
  kp = floor((z(n) - mesh3d%min(3)) * hzi + 1)

  ! Clamp to valid range
  ip = max(1, min(ip, mesh3d%nhi(1) - 1))
  jp = max(1, min(jp, mesh3d%nhi(2) - 1))
  kp = max(1, min(kp, mesh3d%nhi(3) - 1))

  ab  = ((mesh3d%min(1) - x(n)) + ip * mesh3d%delta(1)) * hxi
  de  = ((mesh3d%min(2) - y(n)) + jp * mesh3d%delta(2)) * hyi
  gh  = ((mesh3d%min(3) - z(n)) + kp * mesh3d%delta(3)) * hzi
  oab = 1.0_dp - ab
  ode = 1.0_dp - de
  ogh = 1.0_dp - gh

  ip1 = ip + 1; jp1 = jp + 1; kp1 = kp + 1

  ! Precompute trilinear weights
  w000 = ab  * de  * gh;   w010 = ab  * ode * gh
  w011 = ab  * ode * ogh;  w001 = ab  * de  * ogh
  w101 = oab * de  * ogh;  w111 = oab * ode * ogh
  w110 = oab * ode * gh;   w100 = oab * de  * gh

  if (present(E)) then
    E(1,n) = mesh3d%efield(ip,jp,kp,1)*w000 + mesh3d%efield(ip,jp1,kp,1)*w010 &
           + mesh3d%efield(ip,jp1,kp1,1)*w011 + mesh3d%efield(ip,jp,kp1,1)*w001 &
           + mesh3d%efield(ip1,jp,kp1,1)*w101 + mesh3d%efield(ip1,jp1,kp1,1)*w111 &
           + mesh3d%efield(ip1,jp1,kp,1)*w110 + mesh3d%efield(ip1,jp,kp,1)*w100
    E(2,n) = mesh3d%efield(ip,jp,kp,2)*w000 + mesh3d%efield(ip,jp1,kp,2)*w010 &
           + mesh3d%efield(ip,jp1,kp1,2)*w011 + mesh3d%efield(ip,jp,kp1,2)*w001 &
           + mesh3d%efield(ip1,jp,kp1,2)*w101 + mesh3d%efield(ip1,jp1,kp1,2)*w111 &
           + mesh3d%efield(ip1,jp1,kp,2)*w110 + mesh3d%efield(ip1,jp,kp,2)*w100
    E(3,n) = mesh3d%efield(ip,jp,kp,3)*w000 + mesh3d%efield(ip,jp1,kp,3)*w010 &
           + mesh3d%efield(ip,jp1,kp1,3)*w011 + mesh3d%efield(ip,jp,kp1,3)*w001 &
           + mesh3d%efield(ip1,jp,kp1,3)*w101 + mesh3d%efield(ip1,jp1,kp1,3)*w111 &
           + mesh3d%efield(ip1,jp1,kp,3)*w110 + mesh3d%efield(ip1,jp,kp,3)*w100
  endif

  if (present(B)) then
    B(1,n) = mesh3d%bfield(ip,jp,kp,1)*w000 + mesh3d%bfield(ip,jp1,kp,1)*w010 &
           + mesh3d%bfield(ip,jp1,kp1,1)*w011 + mesh3d%bfield(ip,jp,kp1,1)*w001 &
           + mesh3d%bfield(ip1,jp,kp1,1)*w101 + mesh3d%bfield(ip1,jp1,kp1,1)*w111 &
           + mesh3d%bfield(ip1,jp1,kp,1)*w110 + mesh3d%bfield(ip1,jp,kp,1)*w100
    B(2,n) = mesh3d%bfield(ip,jp,kp,2)*w000 + mesh3d%bfield(ip,jp1,kp,2)*w010 &
           + mesh3d%bfield(ip,jp1,kp1,2)*w011 + mesh3d%bfield(ip,jp,kp1,2)*w001 &
           + mesh3d%bfield(ip1,jp,kp1,2)*w101 + mesh3d%bfield(ip1,jp1,kp1,2)*w111 &
           + mesh3d%bfield(ip1,jp1,kp,2)*w110 + mesh3d%bfield(ip1,jp,kp,2)*w100
    B(3,n) = mesh3d%bfield(ip,jp,kp,3)*w000 + mesh3d%bfield(ip,jp1,kp,3)*w010 &
           + mesh3d%bfield(ip,jp1,kp1,3)*w011 + mesh3d%bfield(ip,jp,kp1,3)*w001 &
           + mesh3d%bfield(ip1,jp,kp1,3)*w101 + mesh3d%bfield(ip1,jp1,kp1,3)*w111 &
           + mesh3d%bfield(ip1,jp1,kp,3)*w110 + mesh3d%bfield(ip1,jp,kp,3)*w100
  endif
enddo
!$OMP END PARALLEL DO

end subroutine interpolate_field_batch


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
! Persistent scratch arrays — only reallocated when doubled mesh size changes
complex(dp), allocatable, save, dimension(:,:,:) :: crho, cgrn
integer, save :: alloc_nx2 = 0, alloc_ny2 = 0, alloc_nz2 = 0
! Green function FFT cache — avoids recomputing when (delta, gamma, offset, mesh size) unchanged
complex(dp), allocatable, save, dimension(:,:,:,:) :: grn_fft_cache  ! (nx2, ny2, nz2, 0:3)
real(dp), save :: grn_cache_delta(3) = 0
real(dp), save :: grn_cache_gamma = -1
real(dp), save :: grn_cache_offset(3) = 0
integer, save :: grn_cache_nx2 = 0, grn_cache_ny2 = 0, grn_cache_nz2 = 0
logical, save :: grn_cache_valid(0:3) = .false.
real(dp) :: factr, offset_eff(3)
real(dp), parameter :: clight=299792458.0
real(dp), parameter :: fpei=299792458.0**2*1.00000000055d-7  ! this is 1/(4 pi eps0) after the 2019 SI changes

integer :: nx, ny, nz, nx2, ny2, nz2
integer :: icomp, ishift, jshift, kshift
logical :: cache_hit

! Sizes
nx = size(rho, 1); ny = size(rho, 2); nz = size(rho, 3)
nx2 = 2*nx; ny2 = 2*ny; nz2 = 2*nz; 

! Allocate complex scratch arrays only when size changes.
! Thread safety: guard with OMP CRITICAL in case this is ever called from a parallel region.
!$OMP CRITICAL (solver2_scratch_lock)
if (nx2 /= alloc_nx2 .or. ny2 /= alloc_ny2 .or. nz2 /= alloc_nz2) then
  if (allocated(crho)) deallocate(crho)
  if (allocated(cgrn)) deallocate(cgrn)
  allocate(crho(nx2, ny2, nz2))
  allocate(cgrn(nx2, ny2, nz2))
  alloc_nx2 = nx2; alloc_ny2 = ny2; alloc_nz2 = nz2
endif
!$OMP END CRITICAL (solver2_scratch_lock)

! Effective offset (treat absent as zero for cache key)
if (present(offset)) then
  offset_eff = offset
else
  offset_eff = 0
endif

! Check if Green function FFT cache is valid for current parameters
cache_hit = (nx2 == grn_cache_nx2 .and. ny2 == grn_cache_ny2 .and. nz2 == grn_cache_nz2 .and. &
             all(delta == grn_cache_delta) .and. gamma == grn_cache_gamma .and. &
             all(offset_eff == grn_cache_offset))
if (.not. cache_hit) then
  grn_cache_valid = .false.
  grn_cache_delta = delta
  grn_cache_gamma = gamma
  grn_cache_offset = offset_eff
  grn_cache_nx2 = nx2; grn_cache_ny2 = ny2; grn_cache_nz2 = nz2
  if (allocated(grn_fft_cache)) then
    if (size(grn_fft_cache, 1) /= nx2 .or. size(grn_fft_cache, 2) /= ny2 .or. &
        size(grn_fft_cache, 3) /= nz2) then
      deallocate(grn_fft_cache)
      allocate(grn_fft_cache(nx2, ny2, nz2, 0:3))
    endif
  else
    allocate(grn_fft_cache(nx2, ny2, nz2, 0:3))
  endif
endif

! rho -> crho -> FFT(crho)
! Zero only the padding region (octants beyond 1:nx, 1:ny, 1:nz), then fill the data octant.
crho(1:nx, 1:ny, 1:nz) = rho
crho(nx+1:nx2, :, :) = 0
crho(1:nx, ny+1:ny2, :) = 0
crho(1:nx, 1:ny, nz+1:nz2) = 0
call ccfft3d(crho, crho, [1,1,1], nx2, ny2, nz2, 0) 

! Loop over phi, Ex, Ey, Ez
do icomp=0, 3
  if ((icomp == 0) .and. (.not. present(phi))) cycle
  if ((icomp == 1) .and. (.not. present(efield))) exit

  if (grn_cache_valid(icomp)) then
    ! Cache hit: multiply directly from cached FFT'd Green function
    cgrn = crho * grn_fft_cache(:,:,:,icomp)
  else
    ! Cache miss: compute Green function, forward FFT, cache, multiply
    call osc_get_cgrn_freespace(cgrn, delta, gamma, icomp, offset=offset_eff)
    call ccfft3d(cgrn, cgrn, [1,1,1], nx2, ny2, nz2, 0)
    grn_fft_cache(:,:,:,icomp) = cgrn
    grn_cache_valid(icomp) = .true.
    cgrn = crho * cgrn
  endif

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
