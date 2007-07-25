#include "CESR_platform.inc"

module csr_mod

use make_mat6_mod
use beam_utils
use bookkeeper_mod, only: attribute_bookkeeper

!+
! See the paper:
!   "An Efficient Formalism for Simulating Coherent Synchrotorn Radiaion"
!   D. Sagan
!-

type csr_kick_factor_struct
  real(rp) g            ! Bending strength = 1 / R at source.
  real(rp) v, v1, v3
  real(rp) w2
  real(rp) theta        ! Kicked particle angle.
  real(rp) L, L_vec(3)  ! Vector between source and kick locations.
end type

! Structure for a single bin.

type csr_bin1_struct   ! Structure for a single particle bin.
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
  real(rp) kick_csr    ! CSR kick
  real(rp) kick_lsc    ! LSC Kick.
end type

! Kicks in this structure are a function of the particle separation

type csr_kick_bin_struct ! Sub-structure for csr calculation cache
  real(rp) I_csr         ! Kick integral.
  real(rp) kick_csr      ! Kick.
  real(rp) kick_lsc      ! Kick.
  real(rp) phi           ! Source point angle.
  real(rp) d             ! Distance between source point and end of element.
  real(rp) dz_particles  ! Distance between source and kicked particles.
  real(rp) s_prime       ! Source point location.
end type

type csr_bin_struct             ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) beta                 ! Relativistic beta factor.
  real(rp) :: dz_bin = 0        ! Bin width
  real(rp) ds_track_step        ! True step size
  real(rp) y2                   ! Height of source particle.
  real(rp) kick_factor           ! Coefficient to scale the kick
  type (csr_bin1_struct), allocatable :: bin1(:)  
  type (csr_kick_bin_struct), allocatable :: kick1(:) ! Array of caches
end type

!

type tsc_struct   ! transverse space charge structure
  type (coord_struct) closed_orb
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) sig_z
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
end type    

!+
! Note: Shielding is simulated via the image current due to the 
! top and bottom walls. The side walls are neglected.
!-

type csr_common_struct                   ! Common block for csr calc
  real(rp) :: ds_track_step = 0          ! Tracking step size
  real(rp) :: beam_chamber_height = 0    ! Used in shielding calculation.
  real(rp) :: sigma_cutoff = 0.1         ! Cutoff for the lsc calc. If a bin sigma
                                         !  is < cutoff * sigma_ave then ignore.
  integer :: n_bin = 0                   ! Number of bins used
  integer :: particle_bin_span = 2       ! Longitudinal particle length / dz_bin
  integer :: n_shield_images = 0         ! Chamber wall shielding. 0 = no shielding.
  logical :: lcsr_component_on = .true.  ! Longitudinal csr component
  logical :: lsc_component_on = .true.   ! Longitudinal space charge component
  logical :: tsc_component_on = .true.   ! Transverse space charge component
  logical :: small_angle_approx = .true. ! Use lcsr small angle approximation?
end type

type (csr_common_struct), save, target :: csr_com

interface 
  subroutine save_bunch_track (bunch, ele, s_travel)
    use beam_def_struct, only: bunch_struct, ele_struct, rp
    implicit none
    type (bunch_struct) bunch
    type (ele_struct) ele
    real(rp) s_travel
  end subroutine
end interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)
!
! Routine to track a bunch of particles through the element lat%ele(ix_ele)
! with csr radiation effects.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bunch_start -- Bunch_struct: Starting bunch position.
!   lat         -- lat_struct: Lattice.
!   ix_ele      -- Integer: lat%ele(ix_ele) is the element to track through.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

subroutine track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)

implicit none

type (lat_struct) lat
type (bunch_struct), target :: bunch_start, bunch_end
type (particle_struct), pointer :: pt
type (ele_struct), pointer :: ele
type (ele_struct), save :: runt
type (csr_bin_struct), save :: bin

real(rp) s_start
integer i, j, ns, nb, ix_ele, n_step

character(20) :: r_name = 'track1_bunch_sc'
logical auto_bookkeeper

! Init

ele => lat%ele(ix_ele)
bunch_end = bunch_start

! n_step is the number of steps to take when tracking through the element.
! bin%ds_step is the true step length.

if (csr_com%n_bin <= csr_com%particle_bin_span + 1) then
  call out_io (s_fatal$, r_name, &
            'CSR_COM%N_BIN MUST BE GREATER THAN CSR_COM%PARTICLE_BIN_SPAN+1!')
  call err_exit
endif

if (csr_com%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'CSR_COM%DS_TRACK_STEP NOT SET!')
  call err_exit
endif

! make sure that ele_len / track_step is an integer.

n_step = max (1, nint(ele%value(l$) / csr_com%ds_track_step))
bin%ds_track_step = ele%value(l$) / n_step

auto_bookkeeper = bmad_com%auto_bookkeeper ! save state
bmad_com%auto_bookkeeper = .false.   ! make things go faster

! Loop over the tracking steps
! runt is the element that is tracked through at each step.

do i = 0, n_step

  call slice_ele_calc (ele, lat%param, i, n_step, runt)

  s_start = i * bin%ds_track_step

  call csr_bin_particles (bunch_end%particle, bin)    

  ! ns = 0 is the unshielded kick.
  ! For the shielding image currents never use the small angle approximation

  do ns = 0, csr_com%n_shield_images

    ! %kick_factor takes into account that at the endpoints we are only putting in a half kick.
    bin%kick_factor = 1
    if (i == 0 .or. i == n_step) bin%kick_factor = 0.5

    bin%y2 = ns * csr_com%beam_chamber_height

    if (ns == 0) then
      call csr_bin_kicks (lat, ix_ele, s_start, bin, csr_com%small_angle_approx)

    else
      ! The factor of two is due to there being image currents both above and below.
      ! The factor of -1^ns accounts for the sign of the image currents
      bin%kick_factor = bin%kick_factor * 2 * (-1)**ns
      call csr_bin_kicks (lat, ix_ele, s_start, bin, .false.)
    endif

    ! loop over all particles and give them a kick

    do j = 1, size(bunch_end%particle)
      pt => bunch_end%particle(j)
      call csr_kick_calc (bin, pt)
    enddo

  enddo

  ! track through the runt

  if (i /= n_step) then
    do j = 1, size(bunch_end%particle)
      pt => bunch_end%particle(j)
      call track1_particle (pt, runt, lat%param, pt)
    enddo
  endif

  call save_bunch_track (bunch_end, ele, s_start)

enddo

bmad_com%auto_bookkeeper = auto_bookkeeper  ! restore state

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
!   particle(:)          -- Particle_struct: Array of particles
!   csr_com              -- Csr_common_struct: CSR common block (not an argument).
!     %n_bin             -- Number of bins.
!     %particle_bin_span -- Particle length / dz_bin. 
!
! Output:
!   bin           -- Csr_bin_struct: The bin structure.
!     %dz_bin     -- Bin longitudinal length
!     %bin1(1:) -- Array of bins.
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
type (csr_bin1_struct), pointer :: bin1

real(rp) z_min, z_max, f, dz_particle, dz
real(rp) zp_center, zp0, zp1, dz_bin2, zb0, zb1, charge, dcharge_ds

integer i, j, n, ix0, ib, ic

character(20) :: r_name = 'csr_bin_particles'

! Init bins...
! The left edge of bin%bin1(1) is at z_min
! The right edge of bin%bin1(n_bin) is at z_max

if (.not. csr_com%lcsr_component_on .and. .not. csr_com%lsc_component_on) return

dz = maxval(particle(:)%r%vec(5)) - minval(particle(:)%r%vec(5)) 
bin%dz_bin = dz / (csr_com%n_bin - (csr_com%particle_bin_span + 1))
bin%dz_bin = 1.0000001 * bin%dz_bin     ! to prevent round off problems
dz_bin2 = bin%dz_bin / 2
z_min = (maxval(particle(:)%r%vec(5)) + minval(particle(:)%r%vec(5))) / 2
z_min = z_min - csr_com%n_bin * bin%dz_bin / 2
z_max = z_min + csr_com%n_bin * bin%dz_bin / 2
dz_particle = csr_com%particle_bin_span * bin%dz_bin

! allocate memeory for the bins

if (allocated(bin%bin1)) then
  if (size(bin%bin1, 1) < csr_com%n_bin) deallocate (bin%bin1)
endif

if (.not. allocated(bin%bin1)) &
    allocate (bin%bin1(csr_com%n_bin), bin%kick1(-csr_com%n_bin:csr_com%n_bin))

! Fill in some z information

do i = 1, csr_com%n_bin
  bin%bin1(i)%z0_edge  = z_min + (i - 1) * bin%dz_bin
  bin%bin1(i)%z_center = bin%bin1(i)%z0_edge + bin%dz_bin / 2
  bin%bin1(i)%z1_edge  = bin%bin1(i)%z0_edge + bin%dz_bin
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

bin%bin1(:)%charge = 0
bin%bin1(:)%dcharge_ds = 0
bin%bin1(:)%x0 = 0
bin%bin1(:)%y0 = 0
bin%bin1(:)%sig_x = 0
bin%bin1(:)%sig_y = 0


f = 2.0 / dz_particle**2

! The contribution to the charge in a bin from a particle is computed from the overlap
! between the particle and the bin.
 
! The gradient of the charge, %dcharge_ds, is computed by using two bins shifted 
! by one half the bin width and assuming a linear gradiant inbetween.

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
    bin1 => bin%bin1(ib)
    zb0 = bin%bin1(ib)%z0_edge
    zb1 = bin%bin1(ib)%z1_edge   ! edges of the bin
    charge = charge_in_bin (zb0, zb1)
    dcharge_ds = (charge_in_bin (zb0+dz_bin2, zb1+dz_bin2) - &
                    charge_in_bin (zb0-dz_bin2, zb1-dz_bin2)) / bin%dz_bin
    bin1%charge = bin1%charge + charge
    bin1%dcharge_ds = bin1%dcharge_ds + dcharge_ds
    bin1%x0 = bin1%x0 + p%r%vec(1) * charge
    bin1%y0 = bin1%y0 + p%r%vec(3) * charge
    ic = ic + 1
    tloc(ic)%charge = charge
    tloc(ic)%x0 = p%r%vec(1)
    tloc(ic)%y0 = p%r%vec(3)
    tloc(ic)%ib = ib
  enddo
enddo

do ib = 1, csr_com%n_bin
  if (bin%bin1(ib)%charge == 0) cycle
  bin%bin1(ib)%x0 = bin%bin1(ib)%x0 / bin%bin1(ib)%charge
  bin%bin1(ib)%y0 = bin%bin1(ib)%y0 / bin%bin1(ib)%charge
enddo

! Compute the particle distribution sigmas in each bin
! Abs is used instead of the usual formula to lessen the effect
! of non-Gaussian tails

do ic = 1, size(tloc)
  if (tloc(ic)%ib < 0) cycle
  bin1 => bin%bin1(tloc(ic)%ib)
  bin1%sig_x = bin1%sig_x + abs(tloc(ic)%x0 - bin1%x0) * tloc(ic)%charge
  bin1%sig_y = bin1%sig_y + abs(tloc(ic)%y0 - bin1%y0) * tloc(ic)%charge
enddo

f = sqrt(pi/2)
do ib = 1, csr_com%n_bin
  bin1 => bin%bin1(ib)
  if (bin1%charge == 0) cycle
  bin1%sig_x = f * bin1%sig_x / bin1%charge
  bin1%sig_y = f * bin1%sig_y / bin1%charge
  bin1%lsc_d0 = bin1%sig_x * bin1%sig_y
  if (bin1%sig_x == 0 .and. bin1%sig_y == 0) then
    bin1%lsc_d1 = 0
  else
    bin1%lsc_d1 = (bin1%sig_x**2 + bin1%sig_y**2) / &
                                        (bin1%sig_x + bin1%sig_y)
  endif
enddo

!---------------------------------------------------------------------------
contains

! computes the contribution to the charge in a bin from
! a given particle.
! z0_bin, z1_bin are the edge positions of the bin

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
! Subroutine csr_bin_kicks (lat, ix_ele, s_travel, bin, small_anlge_approx)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   lat       -- lat_struct: Lattice.
!   ix_ele    -- Integer: lat%ele(ix_ele) is the element to set up cache for.
!   s_travel  -- Real(rp): Distance between the beginning of the element we are
!                  tracking through and the kick point (which is within this element).
!   bin       -- Csr_bin_struct: Binned particle averages.
!     %bin1(:)  -- bin array of particle averages.
!   small_angle_approx -- Logical: If True then use a small angle approximation.
!
! Output:
!   bin     -- Csr_bin_struct: Binned particle averages.
!     %kick1(:) -- CSR kick calculation bin array. 
!     %bin1(:)%kick_csr -- Integrated kick
!-

subroutine csr_bin_kicks (lat, ix_ele, s_travel, bin, small_angle_approx)

implicit none

type (csr_bin_struct), target :: bin
type (lat_struct), target :: lat
type (csr_kick_bin_struct), pointer :: kick1
type (csr_kick_factor_struct) k_factor

real(rp) s_travel, s_kick, s0_kick_ele, coef

integer i, ix_ele, n_ele_pp, n_bin

logical small_angle_approx

character(16) :: r_name = 'csr_bin_kicks'

! The kick point P is fixed.
! The source point P' varies from bin to bin.
! n_ele_pp is the number of elements between P' and P excluding the 
! elements containing P and P'. n_ele_pp will be incremented by
! next_element_params_calc as we go from one bin to the next.

if (.not. csr_com%lcsr_component_on) return

n_ele_pp = -1
call next_element_params_calc (n_ele_pp, s0_kick_ele, k_factor)
s_kick = s0_kick_ele + s_travel  ! absolute s value at point P.

call convert_total_energy_to (lat%ele(ix_ele)%value(E_TOT$), &
                                   lat%param%particle, bin%gamma, beta = bin%beta)
bin%gamma2 = bin%gamma**2

! Loop over all kick1 bins and compute the kick or kick integral.
! The loop steps in increasing dz since that is what next_element_params_calc expects.

do i = lbound(bin%kick1, 1), ubound(bin%kick1, 1)

  kick1 => bin%kick1(i)
  kick1%dz_particles = i * bin%dz_bin   ! Notice negative sign.

  ! Calculate what element the kick point is in.

  do
    kick1%d = d_calc_csr(kick1%dz_particles, k_factor, bin, small_angle_approx)
    kick1%s_prime = s_kick - (kick1%d + k_factor%v1)     ! s value at P'
    if (kick1%s_prime > s0_kick_ele) exit       ! If in element exit loop
    call next_element_params_calc (n_ele_pp, s0_kick_ele, k_factor)
  enddo

  ! calculate csr.
  ! I_csr is only calculated for particles with y = 0 and not for image currents.

  if (bin%y2 == 0) then
    call I_csr (kick1, k_factor, bin)
  else
    call kick_csr (kick1, k_factor, bin)
  endif

enddo

! Now calculate the kick for a particle at the center of a bin.

coef = bin%ds_track_step * r_e / (e_charge * bin%gamma)
n_bin = csr_com%n_bin

if (bin%y2 == 0) then
  call lsc_bin_kicks_y0 (bin)
  do i = 1, n_bin
    bin%bin1(i)%kick_csr = coef * &
                  dot_product(bin%kick1(i-1:0:-1)%I_csr, bin%bin1(1:i)%dcharge_ds)
  enddo

else
  do i = 1, n_bin
    bin%bin1(i)%kick_csr = coef * &
                  dot_product(bin%kick1(i-1:i-n_bin:-1)%kick_csr, bin%bin1(1:n_bin)%charge)
    bin%bin1(i)%kick_lsc = coef * &
                  dot_product(bin%kick1(i-1:i-n_bin:-1)%kick_lsc, bin%bin1(1:n_bin)%charge)
  enddo
endif

!----------------------------------------------------------------------------
contains

subroutine next_element_params_calc (n_ele_pp, s0_kick_ele, k_factor)

type (csr_kick_factor_struct) k_factor
type (ele_struct), pointer :: source_ele

integer i, n_ele_pp, ix_source

real(rp) phi, dphi, s0_kick_ele
real(rp), allocatable, save :: g_i(:), d_i(:)

! n_ele_pp is the number of elements between P' and P excluding the 
! P and P' elements
! Assume a drift before the first element if needed.

n_ele_pp = n_ele_pp + 1
ix_source = ix_ele - n_ele_pp  ! Index of current P' element

source_ele => lat%ele(ix_source)  ! Pointer to the P' element

! Assume a drift before the first element if needed.

if (ix_source == 0) then
  s0_kick_ele = -1e20  ! something large and negative
else
  s0_kick_ele = lat%ele(ix_source-1)%s  ! s value at beginning edge of the P' element
endif

! calculate new values for d and g for this element and store in arrays
! d_i(i) is the length of the ith element
! g_i(i) is the bending radius of the ith element

k_factor%g = 0
if (source_ele%key == sbend$) k_factor%g = source_ele%value(g$)

if (.not. allocated(g_i)) allocate(g_i(100), d_i(100))
if (size(g_i) <= n_ele_pp) then
  call re_allocate (g_i, 2*size(g_i))
  call re_allocate (d_i, 2*size(d_i))
endif

g_i(n_ele_pp+1) = k_factor%g
d_i(n_ele_pp+1) = source_ele%value(l$)

if (n_ele_pp == 0) d_i(1) = s_travel

! calculate new v1, v3, etc.

k_factor%v = 0
k_factor%v1 = 0; k_factor%v3 = 0
k_factor%w2 = 0
k_factor%theta = 0

do i = n_ele_pp, 1, -1

  dphi = d_i(i) * g_i(i)
  k_factor%v1 = k_factor%v1 + d_i(i)

  if (small_angle_approx) then
    k_factor%v3 = k_factor%v3 + d_i(i) * &
                    (k_factor%theta**2 + k_factor%theta*dphi + dphi**2 / 3) / 2
    k_factor%w2 = k_factor%w2 + d_i(i) * (k_factor%theta + dphi/2)
  else
    phi = k_factor%theta 
    if (g_i(i) == 0) then
      k_factor%v = k_factor%v + d_i(i) * cos(phi)
      k_factor%w2 = k_factor%w2 + d_i(i) * sin(phi)
    else
      k_factor%v = k_factor%v + (cos(phi + dphi) - cos(phi)) / g_i(i)
      k_factor%w2 = k_factor%w2 + (sin(phi + dphi) - sin(phi)) / g_i(i)
    endif
    k_factor%v3 = k_factor%v1 - k_factor%v
  endif

  k_factor%theta = k_factor%theta + dphi

enddo

end subroutine

end subroutine
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_bin_kicks_y0 (bin)
!
! Routine to cache intermediate values needed for the lsc calculation.
! This routine is not for image currents.
!
! Modules needed:
!   use lsc_mod
!
! Input:
!   bin       -- lsc_bin_struct: Binned particle averages.
!     %bin1(:)  -- bin array of particle averages.
!
! Output:
!   bin     -- lsc_bin_struct: Binned particle averages.
!     %bin1(:)%kick_lsc -- Integrated kick.
!-

subroutine lsc_bin_kicks_y0 (bin)

implicit none

type (csr_bin_struct), target :: bin
type (csr_bin1_struct), pointer :: bin1

real(rp) sx, sy, a, b, c, dz, factor, sig_x_ave, sig_y_ave

integer i, j

character(16) :: r_name = 'lsc_bin_kicks'

! If there are too few particles in a bin the sigma calc may give a bad value.
! This can be a problem if the computed sigma is small to ignore any bins
! with a very small sigma. 
! To know what is "small" compute the average sigma

sig_x_ave = dot_product(bin%bin1(:)%sig_x, bin%bin1(:)%charge) / sum(bin%bin1(:)%charge)
sig_y_ave = dot_product(bin%bin1(:)%sig_y, bin%bin1(:)%charge) / sum(bin%bin1(:)%charge)

! Compute the kick at the center of each bin

bin%bin1(:)%kick_lsc = 0
if (.not. csr_com%lcsr_component_on) return

do i = 1, csr_com%n_bin
  bin1 => bin%bin1(i)
  sx = bin1%sig_x
  sy = bin1%sig_y
  if (sx < sig_x_ave * csr_com%sigma_cutoff .or. &
                          sy < sig_y_ave * csr_com%sigma_cutoff) cycle
  a = sx * sy
  b = bin%gamma * (sx**2 + sy**2) / (sx + sy)
  c = bin%gamma**2

  do j = 1, csr_com%n_bin
    if (i == j) cycle
    dz = bin%bin1(j)%z_center - bin%bin1(i)%z_center
    bin%bin1(j)%kick_lsc = bin%bin1(j)%kick_lsc + &
                   bin1%charge * sign(1.0_rp, dz) / (a + b * abs(dz) + c * dz**2)
  enddo

enddo

factor = bin%kick_factor * bin%ds_track_step * r_e / (e_charge * bin%gamma)
bin%bin1(:)%kick_lsc = factor * bin%bin1(:)%kick_lsc

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine I_csr (kick1, k_factor, bin) 
!
! Routine to calculate the CSR kick integral.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   kick1 -- Csr_kick_bin_struct: 
!     %dz_particles -- Real(rp): Distance between source and kicked particles.
!     %d            -- Real(rp): Distance between source and kick points.
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!
! Output:
!   kick1 -- Csr_kick_bin_struct: 
!     %I_csr -- Real(rp): CSR kick integral.
!-

subroutine I_csr (kick1, k_factor, bin)

implicit none

type (csr_kick_bin_struct) kick1
type (csr_kick_factor_struct) k_factor
type (csr_bin_struct) bin

real(rp) z, d
real(rp) phi, t, a, k

! 

z = kick1%dz_particles
d = kick1%d

if (z <= 0) then
  kick1%I_csr = 0
  return
endif

phi = k_factor%g * d
t = bin%gamma * (d + k_factor%v1)
a = bin%gamma2 * (k_factor%w2 + phi * k_factor%v1 + phi * d / 2)
k = bin%gamma * (k_factor%theta + phi)
kick1%I_csr = -2 * bin%kick_factor * &
                  bin%gamma * (t + a * k) / (t**2 + a**2) + 1 / (bin%gamma2 * z) 

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kick_csr (kick1, k_factor, bin) result (kick_this)
!
! Routine to calculate the CSR kick integral.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   kick1 -- Csr_kick_bin_struct: 
!     %dz_particles -- Real(rp): Distance between source and kicked particles.
!     %d            -- Real(rp): Distance between source and kick points.
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!
! Output:
!   kick1 -- Csr_kick_bin_struct: 
!     %kick_csr -- Real(rp): CSR kick.
!-

subroutine kick_csr (kick1, k_factor, bin)

implicit none

type (csr_kick_bin_struct) kick1
type (csr_kick_factor_struct), target :: k_factor
type (csr_kick_factor_struct), pointer :: kf
type (csr_bin_struct) bin

real(rp) z, d
real(rp) N_vec(3), G_vec(3), B_vec(3), Bp_vec(3), NBp_vec(3), NBpG_vec(3), rad_cross_vec(3)
real(rp) phi, sin_phi, cos_phi, OneNBp, OneNBp3, radiate, coulomb1, coulomb2

!

kick1%kick_csr = 0
kick1%kick_lsc = 0

z = kick1%dz_particles
d = kick1%d

kf => k_factor
phi = kf%g * d
sin_phi = sin(phi)
cos_phi = cos(phi)

N_vec = k_factor%L_vec / kf%L
B_vec = (/ cos(kf%theta), sin(kf%theta), 0.0_rp /)
Bp_vec = bin%beta * (/ cos_phi, -sin_phi, 0.0_rp /)
G_vec = bin%beta**2 * kf%g * (/ sin_phi, cos_phi, 0.0_rp /)

OneNBp = 1 - sum(N_vec * Bp_vec)
OneNBp3 = OneNBp**3

NBp_vec = N_vec - Bp_vec
NBpG_vec = cross(NBp_vec, G_vec)
rad_cross_vec = cross(N_vec, NBpG_vec)

coulomb2 = bin%gamma * z / ((bin%gamma * z)**2 + bin%y2**2)**(1.5)

if (csr_com%lcsr_component_on) then
  radiate  = dot_product (B_vec, rad_cross_vec) / (kf%L * OneNBp3)
  coulomb1 = dot_product (B_vec, NBp_vec) / (bin%gamma2 * kf%L**2 * OneNBp3)
  kick1%kick_csr = bin%kick_factor * (radiate + coulomb1 - coulomb2)
endif

if (csr_com%lsc_component_on) then
  kick1%kick_lsc = bin%kick_factor * coulomb2
endif

!-----------------------------------------------------------------------------
contains
function cross (a_vec, b_vec) result (c_vec)

real(rp), intent(in) :: a_vec(3), b_vec(3)
real(rp) :: c_vec(3)

c_vec(1) = a_vec(2) * b_vec(3) - a_vec(3) * b_vec(2)
c_vec(2) = a_vec(3) * b_vec(1) - a_vec(1) * b_vec(3)
c_vec(3) = a_vec(1) * b_vec(2) - a_vec(2) * b_vec(1)

end function

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function d_calc_csr (dz_particles, k_factor, bin, small_angle_approx) result (d_this)
!
! Routine to calculate the distance between source and kick points.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   dz_particles -- Real(rp): Distance between source and kicked particles
!   k_factor     -- csr_kick_factor_struct: Value
!   bin          -- Csr_bin_struct:
!   small_angle_approx -- Logical: If True then use a small angle approximation.
!
! Output:
!   d_this -- Real(rp): Distance between source and kick points.
!-

function d_calc_csr (dz_particles, k_factor, bin, small_angle_approx) result (d_this)

implicit none

type (csr_kick_factor_struct) k_factor
type (csr_bin_struct) bin

real(rp) dz_particles, z1, z2, dz_dd, dz1_dd, dz2_dd
real(rp) eps, d_this, delta_d_old, delta_d, dz_calc, d_old
real(rp) d1, d2, a, b, c

integer j

logical small_angle_approx

character(8) :: r_name = 'd_calc_csr'

! If not in a bend then can do the calculation exactly

if (k_factor%g == 0) then
  a = 1 / bin%gamma2
  b = 2 * (k_factor%v3 - dz_particles)
  c = dz_particles**2 - (k_factor%w2**2 + bin%y2**2)
  d_this = (-b + sqrt(b**2 - 4 * a * c)) / (2 * a) - k_factor%v1
  return
endif

! Use Newton's method until root is bracketed.
! Initial d1 is just an approximate guess to get in the ball park.

eps = 1e-10 * abs(dz_particles) + 1e-14
z2 = 0
if (bin%y2 /= 0) then
  d1 = sqrt(3 * bin%y2 / k_factor%g)
elseif (dz_particles >= 0) then
  d1 = min(2 * bin%gamma2 * dz_particles, (6 * dz_particles / k_factor%g**2) ** (0.3333))
else
  d1 = (bin%y2**2 + dz_particles**2) / (2 * dz_particles)
endif

! Now bracket the root

do 

  z1 = z_calc_csr (d1, k_factor, bin, small_angle_approx, dz1_dd) - dz_particles
  if (abs(z1) < eps) then
    d_this = d1
    return
  endif
  if (z1 * z2 < 0) exit  ! have bracket

  d2 = d1 - z1 / dz1_dd 
  z2 = z_calc_csr (d2, k_factor, bin, small_angle_approx, dz2_dd) - dz_particles
  if (abs(z2) < eps) then
    d_this = d2
    return
  endif
  if (z1 * z2 < 0) exit  ! have bracket

  d1 = d2 - z2 / dz2_dd 

enddo

! Find the bracketed root.
! This follows rtsafe from Numerical Recipes.

if (z1 > 0) then
  d_this = d1
  d1 = d2
  d2 = d_this
endif

d_this = (d1 + d2) / 2
delta_d_old = abs (d2 - d1)
delta_d     = delta_d_old
dz_calc = z_calc_csr(d_this, k_factor, bin, small_angle_approx, dz_dd) - dz_particles

do j = 1, 100

  ! Bisect if Newton out of range or not decreasing fast enough

  if (((d_this-d2)*dz_dd-dz_calc) * ((d_this-d1)*dz_dd-dz_calc) > 0 .or. &
                              abs(2*dz_calc) > abs(delta_d_old * dz_dd)) then
    delta_d_old = delta_d
    delta_d = (d2 - d1) / 2
    d_this = d1 + delta_d
    if (d1 == d_this) return

  else
    delta_d_old = delta_d
    delta_d = dz_calc / dz_dd
    d_old = d_this
    d_this = d_this - delta_d
    if (d_this == d_old) return
  endif

  if (abs(dz_calc) < eps) return
  dz_calc = z_calc_csr(d_this, k_factor, bin, small_angle_approx, dz_dd) - dz_particles

  if (dz_calc < 0) then
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
! Function z_calc_csr (d, k_factor, bin, small_angle_approx, dz_dd) result (z_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   d        -- Real(rp): Distance between the source point and the kick point.
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!   small_angle_approx -- Logical: If True then use a small angle approximation.
!
! Output:
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!     %L        -- Only calculated in the non small angle approx case.
!     %L_vec    -- Only calculated in the non small angle approx case.
!   z_this -- Real(rp), Distance between source and kick particles.
!   dz_dd  -- Real(rp), optional: Derivative: dz/dp.
!-

function z_calc_csr (d, k_factor, bin, small_angle_approx, dz_dd) result (z_this)

implicit none

type (csr_kick_factor_struct), target :: k_factor
type (csr_kick_factor_struct), pointer :: kf
type (csr_bin_struct) bin

real(rp) d, phi, phi2, onec, z_this, v1d, w2d, y22
real(rp) RoneMCos, Rsin
real(rp), optional :: dz_dd

logical small_angle_approx

! Special case

if (k_factor%v1 == 0 .and. d == 0 .and. bin%y2 == 0) then
  z_this = 0
  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2)
  return
endif

!

phi = k_factor%g * d
v1d = k_factor%v1 + d

! General case with small angle approx

if (small_angle_approx) then
  w2d = 2*k_factor%w2 - phi*d
  y22 = 4 * bin%y2**2
  z_this = v1d / (2 * bin%gamma2) + &
                      (k_factor%v3 + phi**2 * d / 6 - (w2d**2 + y22)/(8*v1d))

  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2) + &
                      (phi**2/2 + phi*w2d/(2*v1d) + (w2d**2 + y22)/(8*v1d**2))

! General case without small angle approx

else
  kf => k_factor
  if (abs(phi) < 1e-2) then
    phi2 = phi**2
    RoneMCos = phi * d * (1.0/2 - phi2/24 + phi2**2/720)
    Rsin = d * (1 - phi2 / 6 + phi2**2 / 120)
  else
    RoneMCos = (1 - cos(phi)) / kf%g
    Rsin = sin(phi) / kf%g
  endif

  kf%L_vec(1) = Rsin + kf%v
  kf%L_vec(2) = kf%w2 - RoneMCos 
  kf%L_vec(3) = bin%y2
  kf%L = sqrt(sum(kf%L_vec * kf%L_vec))

  z_this = v1d / (2 * bin%gamma2) + (v1d - kf%L)

  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2) + &
               1 - (kf%L_vec(1) * cos(phi) + kf%L_vec(2) * sin(phi)) / kf%L
                

endif

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_kick_calc (bin, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
! Modules needed:
!   use csr_mod
!
!   bin -- Csr_bin1_struct: Binned beam
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Particle_struct: Particle to kick.
!     %r%vec(6) -- Final particle energy.
!-

subroutine csr_kick_calc (bin, particle)

implicit none

type (csr_bin_struct) bin
type (particle_struct) particle

real(rp) zp, r1, r2, dz, dpz
integer i, i0, i_del

! We use a weighted average between %kick1(j)%I_csr and %kick1(j+1)%I_csr
! so that the integral varies smoothly as a function of particle%r%vec(5).

zp = particle%r%vec(5)
i0 = int((zp - bin%bin1(1)%z_center) / bin%dz_bin) + 1
r1 = (zp - bin%bin1(i0)%z_center) / bin%dz_bin
r2 = 1 - r1

if (r1 < 0 .or. r1 > 1 .or. i0 < 1 .or. i0 >= csr_com%n_bin) then
  print *, 'CSR INTERNAL ERROR!'
  call err_exit
endif

if (csr_com%lcsr_component_on) then
  particle%r%vec(6) = particle%r%vec(6) + &
              r2 * bin%bin1(i0)%kick_csr + r1 * bin%bin1(i0+1)%kick_csr
endif

if (csr_com%lsc_component_on) then
  particle%r%vec(6) = particle%r%vec(6) + &
              r2 * bin%bin1(i0)%kick_lsc + r1 * bin%bin1(i0+1)%kick_lsc
endif

end subroutine

end module
