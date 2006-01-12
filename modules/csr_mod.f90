#include "CESR_platform.inc"

module csr_mod

use make_mat6_mod
use beam_def_struct
use bookkeeper_mod, only: attribute_bookkeeper

!

type csr_values_struct
  real(rp) g, v1, v3, w2, theta
end type

! Binning structures for the csr calculation.

type csr_particle_bin_struct   ! Sub-Structure for binning particle averages
  real(rp) x0, y0      ! Transverse center of the particle distrubution
  real(rp) z0_edge     ! Left (min z) edge of bin
  real(rp) z1_edge     ! Right (max z) edge of bin
  real(rp) z_center    ! z at center of bin.
  real(rp) sig_x       ! particle's RMS width
  real(rp) sig_y       ! particle's RMS width
  real(rp) lsc_d0
  real(rp) lsc_d1
  real(rp) charge      ! charge of the particles
  real(rp) dcharge_ds  ! charge derivative
end type

type csr_kick_bin_struct ! Sub-structure for csr calculation cache
  real(rp) I_csr
  real(rp) phi
  real(rp) d             ! Distance between source point and end of element.
  real(rp) dz_particles  ! Distance between source and kicked particles.
  real(rp) s_prime       ! Source point location.
end type

type csr_bin_struct             ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) :: dz_bin = 0        ! Bin width
  real(rp) ds_track_step        ! True step size
  type (csr_particle_bin_struct), allocatable :: bunch1(:)  
  type (csr_kick_bin_struct), allocatable :: kick1(:) ! Array of caches
end type

!

type tsc_struct   ! transverse space charge structure
  type (coord_struct) closed_orb
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
  real(rp) sig_z
endtype    

!

type csr_common_struct                  ! Common block for csr calc
  real(rp) :: ds_track_step = 0         ! Tracking step size
  integer :: n_bin = 0                  ! Number of bins used
  integer :: particle_bin_span = 2      ! Longitudinal particle length / dz_bin
  logical :: lcsr_component_on = .true. ! Longitudinal csr component
  logical :: lsc_component_on = .true.  ! Longitudinal space charge component
  logical :: tsc_component_on = .true.  ! Transverse space charge component
end type

type (csr_common_struct), save, target :: csr_com

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)
!
! Routine to track a bunch of particles through the element lat%ele_(ix_ele)
! with csr radiation effects.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bunch_start -- Bunch_struct: Starting bunch position.
!   lat         -- Ring_struct: Lattice.
!   ix_ele      -- Integer: lat%ele_(ix_ele) is the element to track through.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

subroutine track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)

implicit none

type (ring_struct) lat
type (bunch_struct), target :: bunch_start, bunch_end
type (particle_struct), pointer :: pt
type (ele_struct), pointer :: ele
type (ele_struct), save :: runt
type (csr_bin_struct) bin

real(rp) s_travel
integer i, j, nb, ix_ele, n_step

character(20) :: r_name = 'track1_bunch_sc'

! Init

ele => lat%ele_(ix_ele)
bunch_end = bunch_start

! n_step is the number of steps to take when tracking through the element.
! bin%ds_step is the true step length.

if (csr_com%n_bin == 0) then
  call out_io (s_fatal$, r_name, 'CSR_COM%N_BIN NOT SET!')
  call err_exit
endif

if (csr_com%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'CSR_COM%DS_TRACK_STEP NOT SET!')
  call err_exit
endif

n_step = ele%value(l$) / csr_com%ds_track_step
bin%ds_track_step = ele%value(l$) / n_step  

! runt is the element that is tracked through at each step.

runt = ele
runt%value(l$) = bin%ds_track_step
call attribute_bookkeeper (runt, lat%param)
bmad_com%auto_bookkeeper = .false.   ! make things go faster

! Loop over the tracking steps

do i = 1, n_step

  if (ele%key == sbend$) then
    if (i == 1) then
      runt%value(e2$) = 0
    elseif (i == 2) then
      runt%value(e1$) = 0
    elseif (i == n_step) then
      runt%value(e2$) = ele%value(e2$)
    endif
  endif

  s_travel = bin%ds_track_step * (i - 1)

  if (csr_com%lcsr_component_on .or. csr_com%lsc_component_on) then
    call csr_bin_particles (bunch_end%particle, bin)    
    call csr_bin_kicks (lat, ix_ele, s_travel, bin)
  endif

  ! loop over all particles

  do j = 1, size(bunch_end%particle)
    pt => bunch_end%particle(j)

    ! csr kick

    if (csr_com%lsc_component_on) call lsc_kick (bin, pt)
    if (csr_com%lcsr_component_on) call lcsr_kick (bin, pt)

    ! track through the runt

    call track1_particle (pt, runt, lat%param, pt)

  enddo

enddo

bmad_com%auto_bookkeeper = .false.   

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_particles (particle, bin)
!
! Routine to bin the particles longitudinally in s. 
!
! To avoid noise in the cacluation, every particle is considered to have a 
! triangular distribution with a base length  given by 
!   csr_com%particle_bin_span * bin%dz_bin. 
! That is, particles will, in general, overlap multiple bins. 
!
! Modules needed:
!   use csr_mod
!
! Input:
!   particle(:) -- Particle_struct: Array of particles
!   csr_com     -- Csr_common_struct: CSR common block (not an argument).
!     %n_bin               -- Number of bins.
!     %particle_bin_span   -- Particle length / dz_bin. 
!                               Default is particle_bin_span = 2.
!
! Output:
!   bin     -- Csr_bin_struct: The bin structure.
!     %dz_bin    -- Bin longitudinal length
!     %bunch1(1:) -- Array of bins.
!-

subroutine csr_bin_particles (particle, bin)

implicit none

type this_local_struct   ! Temporary structure 
  real(rp) :: charge     ! how much charge of particle in bin
  real(rp) :: x0, y0     ! particle center
  integer ib             ! bin index
end type

type (particle_struct), target :: particle(:)
type (particle_struct), pointer :: p
type (csr_bin_struct), target :: bin
type (this_local_struct), save, allocatable :: tloc(:)
type (csr_particle_bin_struct), pointer :: bunch1

real(rp) z_min, z_max, f, dz_particle, dz
real(rp) zp_center, zp0, zp1, dz_bin2, zb0, zb1, charge, dcharge_ds

integer i, j, n, ix0, ib, ic

! Init bins...
! The left edge of bin%bunch1(1) is at z_min
! The right edge of bin%bunch1(n_bin) is at z_max

dz = maxval(particle(:)%r%vec(5)) - minval(particle(:)%r%vec(5)) 
bin%dz_bin = dz / (csr_com%n_bin - (csr_com%particle_bin_span + 1))
bin%dz_bin = 1.0000001 * bin%dz_bin     ! to prevent round off problems
dz_bin2 = bin%dz_bin / 2
z_min = (maxval(particle(:)%r%vec(5)) + minval(particle(:)%r%vec(5))) / 2
z_min = z_min - csr_com%n_bin * bin%dz_bin / 2
z_max = z_min + csr_com%n_bin * bin%dz_bin / 2
dz_particle = csr_com%particle_bin_span * bin%dz_bin

! allocate memeory for the bins

if (allocated(bin%bunch1)) then
  if (size(bin%bunch1, 1) < csr_com%n_bin) deallocate (bin%bunch1)
endif

if (.not. allocated(bin%bunch1)) &
    allocate (bin%bunch1(csr_com%n_bin), bin%kick1(-1:csr_com%n_bin))

! Fill in some z information

do i = 1, csr_com%n_bin
  bin%bunch1(i)%z0_edge  = z_min + (i - 1) * bin%dz_bin
  bin%bunch1(i)%z_center = bin%bunch1(i)%z0_edge + bin%dz_bin / 2
  bin%bunch1(i)%z1_edge  = bin%bunch1(i)%z0_edge + bin%dz_bin
enddo

! Init the tloc structure...
! Each particle is distributed longitudinally in a triangular fashion.
! The tloc records how much charge for a given particle is in a bin.

n = size(particle) * (csr_com%particle_bin_span + 2)
if (allocated(tloc)) then
  if (size(tloc) < n) deallocate (tloc)
endif

if (.not. allocated(tloc)) allocate (tloc(n))
tloc%ib = -1

! Compute the particle distribution center in each bin

bin%bunch1(:)%charge = 0
bin%bunch1(:)%dcharge_ds = 0
bin%bunch1(:)%x0 = 0
bin%bunch1(:)%y0 = 0
bin%bunch1(:)%sig_x = 0
bin%bunch1(:)%sig_y = 0


f = 2.0 / dz_particle**2

ic = 0
do i = 1, size(particle)
  p => particle(i)
  if (.not. (p%ix_lost == not_lost$)) cycle
  zp_center = p%r%vec(5) ! center of particle
  zp0 = zp_center - dz_particle / 2       ! particle left edge 
  zp1 = zp_center + dz_particle / 2       ! particle right edge 
  ix0 = nint((zp0 - z_min) / bin%dz_bin)  ! left most bin index
  do j = 0, csr_com%particle_bin_span+1
    ib = j + ix0
    bunch1 => bin%bunch1(ib)
    zb0 = bin%bunch1(ib)%z0_edge
    zb1 = bin%bunch1(ib)%z1_edge   ! edges of the bin
    charge = charge_in_bin (zb0, zb1)
    dcharge_ds = (charge_in_bin (zb0+dz_bin2, zb1+dz_bin2) - &
                    charge_in_bin (zb0-dz_bin2, zb1-dz_bin2)) / bin%dz_bin
    bunch1%charge = bunch1%charge + charge
    bunch1%dcharge_ds = bunch1%dcharge_ds + dcharge_ds
    bunch1%x0 = bunch1%x0 + p%r%vec(1) * charge
    bunch1%y0 = bunch1%y0 + p%r%vec(3) * charge
    ic = ic + 1
    tloc(ic)%charge = charge
    tloc(ic)%x0 = p%r%vec(1)
    tloc(ic)%y0 = p%r%vec(3)
    tloc(ic)%ib = ib
  enddo
enddo

do ib = 1, csr_com%n_bin
  if (bin%bunch1(ib)%charge == 0) cycle
  bin%bunch1(ib)%x0 = bin%bunch1(ib)%x0 / bin%bunch1(ib)%charge
  bin%bunch1(ib)%y0 = bin%bunch1(ib)%y0 / bin%bunch1(ib)%charge
enddo

! Compute the particle distribution sigmas in each bin
! Abs is used instead of the usual formula to lessen the effect
! of non-Gaussian tails

do ic = 1, size(tloc)
  if (tloc(ic)%ib < 0) cycle
  bunch1 => bin%bunch1(tloc(ic)%ib)
  bunch1%sig_x = bunch1%sig_x + abs(tloc(ic)%x0 - bunch1%x0) * tloc(ic)%charge
  bunch1%sig_y = bunch1%sig_y + abs(tloc(ic)%y0 - bunch1%y0) * tloc(ic)%charge
enddo

f = sqrt(pi/2)
do ib = 1, csr_com%n_bin
  bunch1 => bin%bunch1(ib)
  if (bunch1%charge == 0) cycle
  bunch1%sig_x = f * bunch1%sig_x / bunch1%charge
  bunch1%sig_y = f * bunch1%sig_y / bunch1%charge
  bunch1%lsc_d0 = bunch1%sig_x * bunch1%sig_y
  if (bunch1%sig_x == 0 .and. bunch1%sig_y == 0) then
    bunch1%lsc_d1 = 0
  else
    bunch1%lsc_d1 = (bunch1%sig_x**2 + bunch1%sig_y**2) / &
                                        (bunch1%sig_x + bunch1%sig_y)
  endif
enddo

!-------------------------------------------
contains

function charge_in_bin (z0_bin, z1_bin) result (charge)

real(rp) z0_bin, z1_bin, charge, z1, z2

! Integrate over left triangular half of particle distribution

z1 = max(zp0, z0_bin)        ! left integration edge
z2 = min(zp_center, z1_bin)  ! right integration edge
if (z2 > z1) then            ! If left particle half is in bin ...
  charge = f * p%charge * ((z2 - zp0)**2 - (z1 - zp0)**2)
else
  charge = 0
endif

! Integrate over right triangular half of particle distribution

z1 = max(zp_center, z0_bin)  ! left integration edge
z2 = min(zp1, z1_bin)        ! right integration edge
if (z2 > z1) then         ! If right particle half is in bin ...
  charge = charge + f * p%charge * ((z1 - zp1)**2 - (z2 - zp1)**2)
endif

end function

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine subroutine csr_bin_kicks (lat, ix_ele, s_travel, bin)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   lat       -- Ring_struct: Lattice.
!   ix_ele    -- Integer: lat%ele_(ix_ele) is the element to set up cache for.
!   s_travel  -- Real(rp): Distance between the beginning of the element we are
!          tracking through and the kick point (which is within this element). 
!   bin       -- Csr_bin_struct: Binned particle averages.
!     %bunch1(:)   -- bin array of particle averages.
!
! Output:
!   bin     -- Csr_bin_struct: Binned particle averages.
!     %kick1(:) -- CSR kick calculation bin array. 
!-

subroutine csr_bin_kicks (lat, ix_ele, s_travel, bin)

implicit none

type (csr_bin_struct), target :: bin
type (ring_struct), target :: lat
type (csr_kick_bin_struct), pointer :: kick1
type (csr_values_struct) val

real(rp) s_travel, s_kick, s0_kick_ele

integer i, ix_ele, n_ele_pp

character(16) :: r_name = 'csr_cache_setup'

! The kick point P is fixed.
! The source point P' varies from bin to bin.
! n_ele_pp is the number of elements between P' and P excluding the 
! elements containing P and P'. n_ele_pp will be incremented by
! next_element_params_calc as we go from one bin to the next.

n_ele_pp = -1
call next_element_params_calc (n_ele_pp, s0_kick_ele, val)
s_kick = s0_kick_ele + s_travel  ! absolute s value at point P.

call convert_total_energy_to (lat%ele_(ix_ele)%value(beam_energy$), &
                                               lat%param%particle, bin%gamma)
bin%gamma2 = bin%gamma**2

! special case

bin%kick1(-1)%dz_particles = -1 * bin%dz_bin
bin%kick1(-1)%I_csr = 0

! Loop over all kick1 bins

do i = 0, csr_com%n_bin

  kick1 => bin%kick1(i)
  kick1%dz_particles = i * bin%dz_bin

  ! Calculate what element the kick point is in.

  do
    kick1%d = d_calc_csr(kick1%dz_particles, val, bin)
    kick1%s_prime = s_kick - (kick1%d + val%v1)     ! s value at P'
    if (kick1%s_prime > s0_kick_ele) exit       ! If in element exit loop
    call next_element_params_calc (n_ele_pp, s0_kick_ele, val)
  enddo

  ! calculate csr
  
  kick1%I_csr = I_csr (kick1%dz_particles, kick1%d, val, bin)

enddo

!-----------------------------------------------
contains

subroutine next_element_params_calc (n_ele_pp, s0_kick_ele, val)

type (csr_values_struct) val
type (ele_struct), pointer :: ele

integer i, n_ele_pp, ix_source

real(rp) dphi, s0_kick_ele
real(rp) g_i(100), d_i(100)

! n_ele_pp is the number of elements between P' and P excluding the 
! P and P' elements
! Assume a drift before the first element if needed.

n_ele_pp = n_ele_pp + 1
ix_source = ix_ele - n_ele_pp  ! Index of current P' element

ele => lat%ele_(ix_source)  ! Pointer to the P' element

! Assume a drift before the first element if needed.

if (ix_source == 0) then
  s0_kick_ele = -1e20  ! something large and negative
else
  s0_kick_ele = lat%ele_(ix_source-1)%s  ! s value at beginning edge of the P' element
endif

! calculate new values for d and g for this element and store in arrays
! d_i(i) is the length of the ith element
! g_i(i) is the bending radius of the ith element

val%g = 0
if (ele%key == sbend$) val%g = ele%value(g$)
g_i(n_ele_pp+1) = val%g

d_i(n_ele_pp+1) = ele%value(l$)

if (n_ele_pp == 0) d_i(1) = s_travel

! calculate new v1, v3, etc.

val%v1 = 0; val%v3 = 0
val%w2 = 0
val%theta = 0

do i = n_ele_pp, 1, -1
  dphi = d_i(i) * g_i(i)
  val%v1 = val%v1 + d_i(i)
  val%v3 = val%v3 + d_i(i) * (val%theta**2 + val%theta*dphi + dphi**2 / 3) / 2
  val%w2 = val%w2 + d_i(i) * (val%theta + dphi/2)
  val%theta = val%theta + dphi
enddo

end subroutine

end subroutine
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function I_csr (z, d, val, bin) result (I_this)
!
! Routine to calculate the CSR kick integral.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   z   -- Real(rp): Distance between source and kicked particles.
!   d   -- Real(rp): Distance between source and kick points.
!   val -- Csr_values_struct: Other parameters needed in the calculation.
!   bin -- Csr_bin_struct:
!
! Output:
!   I_this -- Real(rp): CSR kick integral.
!-

function I_csr (z, d, val, bin) result (i_this)

implicit none

type (csr_values_struct) val
type (csr_bin_struct) bin

real(rp) z, d
real(rp) phi, t, a, k, i_this

!

  if (z <= 0) then
    I_this = 0
    return
  endif

  phi = val%g * d
  t = bin%gamma * (d + val%v1)
  a = bin%gamma2 * (val%w2 + phi * val%v1 + phi * d / 2)
  k = bin%gamma * (val%theta + phi)
  I_this = -2 * bin%gamma * (t + a * k) / (t**2 + a**2) + 1 / (bin%gamma2 * z) 

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function z_calc_csr (d, val, bin, dz_dd) result (z_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   d   -- Real(rp): Distance between the source point and the kick point.
!   val -- Csr_values_struct: Other parameters needed in the calculation.
!   bin -- Csr_bin_struct:
!
! Output:
!   z_this -- Real(rp), Distance between source and kick particles.
!   dz_dd  -- Real(rp), optional: Derivative: dz/dp.
!-

function z_calc_csr (d, val, bin, dz_dd) result (z_this)

implicit none

type (csr_values_struct) val
type (csr_bin_struct) bin

real(rp) d, phi, z_this, v1d, w2d
real(rp), optional :: dz_dd

! Special case

if (val%v1 == 0 .and. d == 0) then
  z_this = 0
  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2)
  return
endif

! General case

phi = val%g * d
v1d = val%v1 + d
w2d = 2*val%w2 - phi*d
z_this = v1d / (2 * bin%gamma2) + (val%v3 + phi**2 * d / 6 - w2d**2/(8*v1d))

if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2) + &
                      (phi**2/2 + phi*w2d/(2*v1d) + w2d**2/(8*v1d**2))

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function d_calc_csr (dz_particles, val, bin) result (d_this)
!
! Routine to calculate the distance between source and kick points.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   dz_particles -- Real(rp): Distance between source and kicked particles
!   val          -- Csr_valuse_struct: 
!   bin          -- Csr_bin_struct:
!
! Output:
!   d_this -- Real(rp): Distance between source and kick points.
!-

function d_calc_csr (dz_particles, val, bin) result (d_this)

implicit none

type (csr_values_struct) val
type (csr_bin_struct) bin

real(rp) dz_particles, z1, z2, dz_dd, dz1_dd, dz2_dd
real(rp) eps, d_this, delta_d_old, delta_d, z_this, d_old
real(rp) d1, d2, a, b, c

integer j

character(8) :: r_name = 'd_calc_csr'

! If not in a bend then can do the calculation exactly

if (val%g == 0) then
  a = 1 / bin%gamma2
  b = 2 * (val%v3 - dz_particles)
  c = -val%w2**2
  d_this = (-b + sqrt(b**2 - 4 * a * c)) / (2 * a) - val%v1
  return
endif

! Use Newton's method until root is bracketed.
! Initial d1 is just an approximate guess to get in the ball park.

eps = 1e-10 * abs(dz_particles) + 1e-14
z2 = 0
d1 = min(2 * bin%gamma2 * dz_particles, (6 * dz_particles / val%g**2) ** (0.3333))

do 

  z1 = z_calc_csr (d1, val, bin, dz1_dd) - dz_particles
  if (abs(z1) < eps) then
    d_this = d1
    return
  endif
  if (z1 * z2 < 0) exit  ! have bracket

  d2 = d1 - z1 / dz1_dd 
  z2 = z_calc_csr (d2, val, bin, dz2_dd) - dz_particles
  if (abs(z2) < eps) then
    d_this = d2
    return
  endif
  if (z1 * z2 < 0) exit  ! have bracket

  d1 = d2 - z2 / dz2_dd 

enddo

! Use rtsafe from Numerical Recipes

if (z1 > 0) then
  d_this = d1
  d1 = d2
  d2 = d_this
endif

d_this = (d1 + d2) / 2
delta_d_old = abs (d2 - d1)
delta_d     = delta_d_old
z_this = z_calc_csr(d_this, val, bin, dz_dd) - dz_particles

do j = 1, 100

  ! Bisect if Newton out of range or not decreasing fast enough

  if (((d_this-d2)*dz_dd-z_this) * ((d_this-d1)*dz_dd-z_this) > 0 .or. &
                              abs(2*z_this) > abs(delta_d_old * dz_dd)) then
    delta_d_old = delta_d
    delta_d = (d2 - d1) / 2
    d_this = d1 + delta_d
    if (d1 == d_this) return

  else
    delta_d_old = delta_d
    delta_d = z_this / dz_dd
    d_old = d_this
    d_this = d_this - delta_d
    if (d_this == d_old) return
  endif

  if (abs(z_this) < eps) return
  z_this = z_calc_csr(d_this, val, bin, dz_dd) - dz_particles

  if (z_this < 0) then
    d1 = d_this
  else
    d2 = d_this
  endif

enddo

call out_io (s_abort$, r_name, 'ALGORITHM NOT CONVERGING.')
call err_exit

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_kick (bin, particle)
!
! Routine to calculate the longitudinal space charge kick.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bin -- Csr_particle_bin_struct: Binned beam
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Final particle energy.
!-

subroutine lsc_kick (bin, particle)

implicit none

type (csr_bin_struct), target :: bin
type (csr_particle_bin_struct), pointer :: bunch1
type (particle_struct) particle

real(rp) x, y, sx, sy, a, b, dz, sg, f

integer i

!

do i = 1, csr_com%n_bin

  bunch1 => bin%bunch1(i)
  if (bunch1%charge == 0) cycle

  x = particle%r%vec(1) - bunch1%x0
  y = particle%r%vec(3) - bunch1%y0
  dz = abs(particle%r%vec(5) - bunch1%z_center) 
  sg = bin%gamma * abs(dz)
  sx = bunch1%sig_x
  sy = bunch1%sig_y
  if (sx == 0 .or. sy == 0) then
    a = 0
    b = 0
  else
    a = sx * sy * exp(((x/sx)**2 + (y/sy)**2) / 2)
    b = bin%gamma * (sx**2 + sy**2) / (sx + sy)
  endif

  f = bunch1%charge * bin%ds_track_step / (a + b * sg + bin%gamma2 * sg**2)   
  f = f * r_e / (e_charge * bin%gamma)

! If the particle position is within the bin then use a linear interpolation.

  if (particle%r%vec(5) < bunch1%z0_edge) then
    particle%r%vec(6) = particle%r%vec(6) - f
   
  elseif (particle%r%vec(5) > bunch1%z1_edge) then
    particle%r%vec(6) = particle%r%vec(6) + f

  else
    particle%r%vec(6) = particle%r%vec(6) + f * (2 * dz / bin%dz_bin)

  endif

end do

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lcsr_kick (bin, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
! Modules needed:
!   use csr_mod
!
!   bin -- Csr_particle_bin_struct: Binned beam
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Final particle energy.
!-

subroutine lcsr_kick (bin, particle)

implicit none

type (csr_bin_struct) bin
type (particle_struct) particle

real(rp) dz_bin, zp, r1, r2, dz, I_csr, dpz
integer i, i0, i_del

! the longitudinal csr kick is computed via integrating:
!    dcharge_ds * I_csr 
! The trick here is to compute for a given bin%bunch1(i)%dchearge_ds the
! corresponding j for bin%kick1(j)%I_csr.
! We actually use a weighted average between %kick1(j)%I_csr and %kick1(j+1)%I_csr
! so that the integral varies smoothly as a function of particle%r%vec(5).

zp = particle%r%vec(5)
i0 = int((zp - bin%bunch1(1)%z0_edge) / bin%dz_bin) + 1
r2 = modulo2 ((zp - bin%bunch1(1)%z_center) / bin%dz_bin, 0.5_rp)

if (r2 > 0) then
  i_del = i0
  r1 = 1 - r2
else
  i_del = i0 - 1
  r1 = -r2
  r2 = 1 - r1
endif

! Loop over all bins and integrate. We take advantage of the fact that 
! I_csr = 0 for source particles ahead of the kicked particle.

dpz = 0
do i = 1, i0
  I_csr = bin%kick1(i_del-i)%I_csr * r1 + bin%kick1(i_del-i+1)%I_csr * r2
  dpz = dpz + I_csr * bin%bunch1(i)%dcharge_ds
enddo

particle%r%vec(6) = particle%r%vec(6) + &
      dpz * bin%ds_track_step * r_e / (e_charge * bin%gamma)

end subroutine

end module
