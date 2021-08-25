! This code is part of the OpenCSR package developed at:
! https://github.com/ChristopherMayes/OpenCSR


module csr3d_mod

use, intrinsic :: iso_fortran_env
use elliptic_integral_mod, only : ellipinc 
use fft_interface_mod, only : ccfft3d
!$ use omp_lib, only: omp_get_max_threads

implicit none

! Fortran 2008
!integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
!integer, parameter, private :: qp = REAL128


!------------------------------------------------------------------------
!------------------------------------------------------------------------
contains


!------------------------------------------------------------------------
!+
! Subroutine csr3d_steady_state_solver(density, gamma, rho, delta, wake, offset, normalize)
!
! Deposits particle arrays onto mesh
!
! Input:
!   density      -- REAL64(:,:,:): charge density array in x, y, z
!                                  This should be normalized so that:
!                                  sum(density)*product(delta) = 1
!
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!
!   gamma        -- REAL64: relativistic gamma
!
!   rho          -- REAL64: bending radius. Can be positive or negative.
!
!   icomp        -- integer: Field component requested:
!                        0: phi (scalar potential)
!                       
!   wake        -- REAL64(:,:,:,3): allocated wakefield array to populate.
!                      
!                                     The final index corresponds to components
!                                     1: Wx
!                                     2: Wy
!                                     3: Ws                                   
!                                     If present, all components will be computed.    
!
!   phi           -- REAL64(:,:,:), optional: allocated potential array to populate
!
!   normalize     -- REAL64, optional: Return normalized units (1/m^2)
!                                      Otherwise further multiply by 1/(4*pi*epsilon_0 * rho)
!                                      to return V/m
!                    Default: True
!
!   offset        -- real(3), optional: Offset coordinates x0, y0, z0 to evaluate the field,
!                    relative to rho. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to rho. 
!
! Output:
!   wake        -- REAL64(:,:,:,3) :  3D CSR wake                         
!
!
! Notes: 
!
!   This algorithm is ultimately memory limited. For each component, two complex,
!   double-sized arrays are needed, requiring the following:
!   Bytes per component calc = nx*ny*nz * 2 (arrays of complex numbers) *
!            2^3 (double-sized 3D arrays) * 128 (bits/complex number) / 8 (bits/Bytes)
!         = 256 Bytes * nx*ny*nz
!   Some examples:
!      32x32x32    =>   8 MB
!      64x64x64    =>  67 MB
!      128x128x128 => 538 MB
!      256x256x256 => 4.3 GB
!      512x512x512 =>  34 GB
!   For smaller grid sizes, the algorithm could be improved by calculating all three
!   Green function arrays at the same time. 
!
!-
subroutine csr3d_steady_state_solver(density, gamma, rho, delta, wake, offset, normalize)


real(dp), intent(in), dimension(:,:,:) :: density
real(dp), intent(in) :: gamma, rho, delta(3)
real(dp), optional, intent(out), dimension(:,:,:,:) :: wake
real(dp), intent(in), optional :: offset(3)
logical, optional :: normalize
logical :: normalize0
! internal arrays
real(dp), allocatable, dimension(:,:,:) :: density_prime
complex(dp), allocatable, dimension(:,:,:) :: complex_density_prime, cgrn
real(dp) :: dx, dy, dz
real(dp) :: factor, offset0=0
real(dp), parameter :: clight=299792458.0
real(dp), parameter :: fpei=299792458.0**2*1.00000000055d-7  ! this is 1/(4 pi eps0) after the 2019 SI changes

integer :: nx, ny, nz, nx2, ny2, nz2
integer :: i, icomp, ishift, jshift, kshift

! optional normalization
normalize0 = .true.
if (present(normalize)) normalize0 = normalize

! Sizes
nx = size(density, 1); ny = size(density, 2); nz = size(density, 3)
nx2 = 2*nx; ny2 = 2*ny; nz2 = 2*nz; 

! Grid spacings (conveniences)
dx = delta(1); dy = delta(2); dz = delta(3)




! Allocate complex scratch arrays
!allocate(density_prime(nx, ny, nz))
allocate(cgrn(nx2, ny2, nz2))
allocate(complex_density_prime(nx2, ny2, nz2))


! density -> d/dz density -> FFT 
! This will calculate d/dz density and place in the complex array
call calc_density_derivative_complex(density, complex_density_prime, dz)
!call calc_density_derivative(density, density_prime, delta(3))
!complex_density_prime = 0
!complex_density_prime(1:nx, 1:ny, 1:nz) = density_prime

! Debug:
!call write_2d('density', density(:, ny/2, :))
!call write_2d('density_prime', density_prime(:, ny/2, :))
!call write_2d('complex_density_prime', real(complex_density_prime(1:nx, ny/2, 1:nz)))

! FFT complex density 
call ccfft3d(complex_density_prime, complex_density_prime, [1,1,1], nx2, ny2, nz2, 0) 

! Loop over Gx, Gy, Gs
do icomp=1, 3
  
  ! Get the green array
  call get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset=offset)

  ! This step fails because of NaN
  !  cgrn -> FFT(cgrn)
  call ccfft3d(cgrn, cgrn, [1,1,1], nx2, ny2, nz2, 0)  
   
  print *,  "Multiply FFT'd arrays, re-use cgrn"
  ! Multiply FFT'd arrays, re-use cgrn
  cgrn=complex_density_prime*cgrn

  ! Inverse FFT
  call ccfft3d(cgrn, cgrn, [-1,-1,-1], nx2, ny2, nz2, 0)  


  ! Factor for normalized 1/m^2 units, or V/m
  if (normalize0) then
    factor = dx*dy*dz/(nx2*ny2*nz2)
  else
    ! Green function handles negative rho
    factor = fpei/abs(rho) * dx*dy*dz/(nx2*ny2*nz2)
  endif
  
  ! This is where the output is shifted to
  ishift = nx-1
  jshift = ny-1
  kshift = nz-1

 ! Extract field
  wake(:,:,:,icomp) = factor*real(cgrn(1+ishift:nx+ishift, 1+jshift:ny+jshift, 1+kshift:nz+kshift), dp)

enddo

end subroutine csr3d_steady_state_solver



subroutine calc_density_derivative_complex(density, density_prime, dz)
real(dp), intent(in), dimension(:,:,:) :: density
real(dp), intent(in) :: dz
complex(dp), dimension(:,:,:), intent(inout) :: density_prime
integer :: nx, ny, nz

print *, 'calc_density_derivative...'

! Sizes
nx = size(density, 1)
ny = size(density, 2)
nz = size(density, 3)

! Allocate double-sized arrays. This array will be placed in the first octant
!allocate( density_prime( 2*size(density, 1), 2*size(density, 2), 2*size(density, 3) ))
density_prime = 0 ! Zero out everything

! Endpoints: second order Forward and backwards stencils
density_prime(1:nx,1:ny,1)  = (-(3/2)*density(:,:,1)    + 2*density(:,:,2)    -(1/2)*density(:,:,3) )/dz
density_prime(1:nx,1:ny,nz) = ( (1/2)*density(:,:,nz-2) - 2*density(:,:,nz-1) +(3/2)*density(:,:,nz) )/dz

! Second order central differences
density_prime(1:nx,1:ny,2:nz-1) =  (density(:,:,3:nz) - density(:,:,1:nz-2))/(2*dz)

print *, 'calc_density_derivative...Done'

end subroutine



!------------------------------------------------------------------------
!+
! Subroutine get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset)
!
! Computes the free space Green function on a mesh with given spacings in the lab frame.
! The computation is performed in the rest fram by boosting the coordinates by gamma.
!
!
! Input:
!   cgrn         -- COMPLEX128(:,:,:): pre-allocated array 
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!   gamma        -- REAL64: relativistic gamma
!   rho          -- REAL64: bending radius. Can be positive or negative.
!   icomp        -- integer: Wake component requested:
!
!                        1: Wx
!                        2: Wy
!                        3: Ws
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
!   Internal formulas are for rho > 0
!   
!   Internally, mesh spacings in scaled units
!   x - > x/rho
!   y ->  y/rho
!   z ->  z/(2*rho)   
!
!
!-
subroutine get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset)

complex(dp), intent(out), dimension(0:,0:,0:) :: cgrn ! Convenient indexing 
real(dp), intent(in), dimension(3) :: delta
integer, intent(in) :: icomp
real(dp), intent(in) :: gamma, rho
real(dp), intent(in), optional :: offset(3)
! Local
real(dp) :: dx,dy,dz
real(dp) :: x, y, z, xmin, ymin, zmin
real(dp) :: gval
integer :: imin, imax, jmin, jmax, kmin, kmax
integer :: i,j,k, isize, jsize, ksize
integer :: n_threads

! Mesh spacings in scaled units
! x - > x/rho
! y ->  y/rho
! z ->  z/(2*rho)


dx=delta(1)/rho; dy=delta(2)/abs(rho); dz=delta(3)/abs(2*rho) ! This will handle negative rho

! Grid min
isize = size(cgrn,1); jsize=size(cgrn,2); ksize=size(cgrn,3)

xmin = ( -isize/2 +1 ) *dx
ymin = ( -jsize/2 +1 ) *dy
zmin = ( -ksize/2 +1 ) *dz

! Add optional offset
if (present(offset)) then
  xmin = xmin + offset(1)/abs(rho)  ! Do not flip the sign of the offset. 
  ymin = ymin + offset(2)/abs(rho)
  zmin = zmin + offset(3)/abs(2*rho)
endif



!$ n_threads = omp_get_max_threads()
!$ if (n_threads > 1) write(*,'(a, i0)') 'OpenMP Green function calc get_cgrn_csr3d, n_threads = ', n_threads
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k = 0, ksize-1
  z = k*dz + zmin

  do j=0, jsize-1
    y=j*dy + ymin

    do i=0, isize-1
      x = i*dx + xmin
      cgrn(i,j,k) =  psi0(x, y, z, gamma, icomp, abs(dx), dy, dz)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO

! Flip psi_x sign for negative rho
if (icomp == 1 .and. rho < 0) cgrn = -cgrn

!call write_2d('psi', real(cgrn(:, jsize/2-1 +1, :)))

end subroutine get_cgrn_csr3d



!------------------------------------------------------------------------
!+
! elemental real(dp) function psi_s(x, y, z, gamma)
!
!
! Eq. 24 from Ref[X] without the prefactor e beta^2 / (2 rho^2)
!-
elemental real(dp) function psi_s(x, y, z, gamma)
real(dp), intent(in) :: x, y, z, gamma
real(dp) :: beta, beta2, kap, alp
    
  
    
    
! Check for the origin
if ((x == 0) .and. (y == 0) .and. (z == 0)) then
    psi_s = 0
else
    beta2 = 1-1/gamma**2
    beta = sqrt(beta2)  
    
    alp = alpha(x, y, z, gamma)
    
    kap = 2*(alp - z)/beta ! Simpler form of kappa
    !kap = sqrt(x**2 + y**2 + 4*(1+x) * sin(alp)**2) 
    
    psi_s =  (cos(2*alp) - 1/(1+x)) / (kap - beta*(1+x)*sin(2*alp) )
    
endif
end function psi_s









!------------------------------------------------------------------------
!+
! function psi0(x, y, z, gamma, icomp, dx, dy, dz)
!
! Same as psi, but handles singularities along the coordinate axes.
! This should be used when creating a grid of points. 
!

elemental real(dp) function psi0(x, y, z, gamma, icomp, dx, dy, dz)
real(dp), intent(in) :: x, y, z, gamma, dx, dy, dz
integer, intent(in) :: icomp
select case (icomp)

    case(1, 2)
    ! x, y
    ! There are singularities along these axes. Do simple average
    ! TODO: tanh-sinh quadrature
    if ((x == 0) .and. (y == 0)) then
        ! Along z axis
        psi0 = (psi(-dx/2, y, z, gamma, icomp) + psi( dx/2, y, z, gamma, icomp))/2
 
    elseif ((y == 0) .and. (z == 0)) then
        ! Along x axis
        psi0 = (psi(x, -dy/2, z, gamma, icomp) + psi(x, dy/2, z, gamma, icomp))/2
    else    
        psi0 = psi(x, y, z, gamma, icomp)
    endif
    
    case(3)
    ! s
    if ((x == 0) .and. (y == 0) .and. (z == 0)) then
        ! Origin = 0
        psi0 = 0
    else
        psi0 = psi(x, y, z, gamma, icomp)
    endif
    
    
end select    
end function



! Conveniences

real(dp) function psi_x0(x, y, z, gamma, dx, dy, dz)
real(dp), intent(in) :: x, y, z, gamma, dx, dy, dz
psi_x0 = psi0(x, y, z, gamma, 1, dx, dy, dz)
end function

real(dp) function psi_y0(x, y, z, gamma, dx, dy, dz)
real(dp), intent(in) :: x, y, z, gamma, dx, dy, dz
psi_y0 = psi0(x, y, z, gamma, 2, dx, dy, dz)
end function


real(dp) function psi_s0(x, y, z, gamma)
real(dp), intent(in) :: x, y, z, gamma
real(dp) :: dx=0, dy=0, dz=0
    psi_s0 = psi0(x, y, z, gamma, 3, dx, dy, dz)
end function





!------------------------------------------------------------------------
!+
! elemental real(dp) function psi(x, y, z, gamma, icomp)
!
! Returns a component of the 3D CSR potential psi from Ref. [X]
!
! Eq. 24 from Ref[X] without the prefactor e beta^2 / (2 rho^2)
! This is actually psi_x_hat
!
!
! Inputs
! ------
! Same as alpha plus:
!
! icomp : integer
!          1: Wx
!          2: Wy
!          3: Ws
!
! Returns:
! psi : real
!      psi_x if icomp == 1
!      psi_y if icomp == 2
!      psi_s if icomp == 3
!
!-
elemental real(dp) function psi(x, y, z, gamma, icomp)
real(dp), intent(in) :: x, y, z, gamma
integer, intent(in) :: icomp
real(dp) :: beta, beta2, kap, alp
real(dp) :: sin2a, cos2a, kap2, sin2a2, x2, y2, y4, xp, xp2, xy, xy2, f1, f2, arg2, f, e

beta2 = 1-1/gamma**2
beta = sqrt(beta2)

alp = alpha(x, y, z, gamma)
kap = 2*(alp - z)/beta ! Simpler form of kappa
!kap = sqrt(x**2 + y**2 + 4*(1+x) * sin(alp)**2) 

if (icomp == 3) then
    if ((x == 0) .and. (y == 0) .and. (z == 0)) then
      psi = 0
    else
      psi = (cos(2*alp) - 1/(1+x)) / (kap - beta*(1+x)*sin(2*alp) )
    endif
    return
endif

! Common patterns
sin2a = sin(2*alp)
cos2a = cos(2*alp)

kap2 = kap**2
sin2a2 = sin2a**2

x2 = x**2 
y2 = y**2
y4 = y2**2
xp = x + 1
xp2 = xp**2
xy2 = x2 + y2
xy = sqrt(xy2)

! More complicated pattens
f1 = 2 + 2*x +x2
f2 = (2+x)**2
arg2 = -4 * xp / xy2 

! Elliptic integrals of the first and second kind F(phi|m), E(phi|m)
call ellipinc(alp, arg2, F, E)

if (icomp == 1) then
! psi_x    
! psi_x is actually psi_x_hat that includes the psi_phi term
! There is an extra ] in the numerator of the second term. All terms should multiply E. 

    psi = f1*F / (xp*xy) - (x2*f2 + y2*f1)*E / (xp*(y2+f2)*xy)  &
        + ( kap2 - 2*beta2*xp2 + beta2*xp*f1*cos2a  ) / (beta *xp*(kap2 - beta2*xp2*sin2a2)) &
        + kap*( y4 - x2*f2 - 2*beta2*y2*xp2 )*sin2a / ( xy2*(y2 + f2)*(kap2-beta2*xp2*sin2a2)  ) &
        + kap*beta2*xp*( x2*f2 + y2*f1 )*sin2a*cos2a / ( xy2*(y2+f2)*(kap2-beta2*xp2*sin2a2)  ) &
        - (2/beta2)* F/xy ! Include the phi term      

elseif (icomp == 2) then
! psi_y
    psi = y * ( &
            F/xy - (x*(2+x)+y2)*E / ((y2+f2)*xy) &
            - beta*(1-xp*cos2a) / (kap2-beta2*xp2*sin2a2) &
            + kap*xp*( -(2+beta2)*y2 + (-2+beta2)*x*(2+x) ) * sin2a / ( (y4 + x2*f2 + 2*y2*f1)*( kap2-beta2*xp2*sin2a2 ) ) &
            + kap*beta2*xp2*(y2 + x*(2+x))*sin2a*cos2a / ( ( y4 + x2*f2 + 2*y2*f1)*(kap2 -beta2*xp2*sin2a2)  ) &
            )      
endif


end function psi





!------------------------------------------------------------------------
!+
! elemental real(dp) function alpha(x, y, z, gamma)
! 
! alpha angle from Eq. 6 in Ref. [X]
!
! When z = 0, this uses the the quadratic solution. 
! Otherwise, it uses the quartic solution in Appendix A. 
!
!
  
! Inputs
! ------
!   x  : real
!        horizontal position divided by rho
!        This is chi in Ref. [X]
!   y  : real
!        vertical position divided by rho 
!        This is zeta in Ref. [X]
!   z  : real
!        longitudinal position divided by 2*rho
!        This is xi in Ref. [X]
!
!   gamma : real
!           realativistic gamma
!
! Returns
! -------
!  alpha : real
!
!-
elemental real(dp) function alpha(x, y, z, gamma)
real(dp), intent(in) :: x, y, z, gamma
real(dp) :: beta2, b, c, eta, nu, zeta,  omega, omega3, m
real(dp) :: arg1, arg2, arg3, temp

beta2 = 1-1/gamma**2

if (z == 0.0) then
    ! Quadratic solution 
    
    b = 3 * (1 - beta2 - beta2*x) / beta2 / (1+x)    
    c = -3*(x**2 + y**2)/(4*(1+x))

    alpha = sqrt( (-b + sqrt(b**2 - 4*c))/2 )
    
else    
    !Quartic solution 
                          
    ! Terms of the depressed quartic equation
    eta = -6 * z / (beta2 * (1+x))
    nu = 3 * (1/beta2 - 1 - x) / (1+x)
    zeta = (3.0_dp/4.0_dp) * (4* z**2 /beta2 - x**2 - y**2) / (1+x)
     
    ! Omega calc and cube root
    temp = (eta**2/16 - zeta * nu/6 + nu**3/216)  
    omega = temp + sqrt(temp**2 - (zeta/3 + nu**2/36)**3)
    omega3 = omega**(1.0_dp/3.0_dp)   ! Do not forget the 1, 3 literals must be defined to double precision!!  

    ! Eq. (A2) from Ref[1]
    m = -nu/3 + (zeta/3 + nu**2/36) /omega3 + omega3
     
    arg1 = sqrt(2 * abs(m))
    arg2 = -2 * (m + nu)
    arg3 = 2 * eta / arg1
        
    if (z >= 0) then
        alpha =  (arg1 + sqrt(abs(arg2 - arg3)))/2
    else
        alpha = (-arg1 + sqrt(abs(arg2 + arg3)))/2
    endif

endif
end function alpha




!------------------------------------------------------------------------
! Debugging utils


! -----------------
! -----------------
! Write grid utility for debugging
!
subroutine write_2d(fname, grid)
real(dp) :: grid(:,:)
integer :: i, j, outfile
character(*) :: fname

open(newunit=outfile, file = trim(fname))

write(outfile, *)  size(grid, 1)
write(outfile, *)  size(grid, 2)
do i = 1, size(grid, 1)
  do j = 1, size(grid, 2)
    write(outfile, *) grid(i, j)
  enddo
enddo  
close(outfile)
write(*,*) 'Grid written to ', trim(fname)  
end subroutine
! -----------------
! -----------------



subroutine ellipinc_test()
real(dp) :: phi, m, mc, elb, eld
real(dp) :: ellipkinc, ellipeinc
phi = 0.1
m = 0.5

call ellipinc(phi, m, ellipkinc, ellipeinc)
print *, 'ellipinc(phi, m, ellipkinc, ellipeinc)', phi, m, ellipkinc, ellipeinc

print *, ' ------ negative m test ---------' 
phi = 0.1761732772710074
m = -2001999.9999999998

call ellipinc(phi, m, ellipkinc, ellipeinc)
print *, 'ellipinc(phi, m, ellipkinc, ellipeinc)', phi, m, ellipkinc, ellipeinc


!     Inputs: phi = amplitude, mc = complementary parameter, 0 <= m < 1
!
!     Outputs: elb = B(phi|m), eld = D(phi|m)

end subroutine



end module


