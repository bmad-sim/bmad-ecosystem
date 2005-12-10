#include "CESR_platform.inc"

module csr_mod

use make_mat6_mod
use beam_def_struct
use bookkeeper_mod, only: attribute_bookkeeper

type csr_local_cache_struct
  real(rp) gamma, gamma2
  real(rp) g, v1, v3, w2, theta, c_const
end type

type csr_bin1_struct   ! structure for a single bin
  real(rp) x0, y0      ! Center of the particle distrubution
  real(rp) z1_edge     ! Left (min z) edge of bin
  real(rp) z2_edge     ! Right (max z) edge of bin
  real(rp) z           ! z at center of bin.
  real(rp) sig_x       ! particle's RMS width
  real(rp) sig_y       ! particle's RMS width
  real(rp) lsc_d0
  real(rp) lsc_d1
  real(rp) charge      ! charge of the particles
  real(rp) dcharge_ds  ! charge derivative
end type

type csr_bin_struct             ! Structure for the bins
  real(rp) dz_bin                 ! Bin width
  integer :: n_bin_particle = 2   ! Longitudinal particle length / dz_bin
  integer n_bin                   ! Number of bins used
  real(rp) gamma, gamma2
  type (csr_bin1_struct), allocatable :: bin1(:)  ! Array of bins
end type

type csr_cache1_struct
  real(rp) I_csr
  real(rp) phi
  real(rp) z        ! Distance between source and kicked particles.
  real(rp) s_prime  ! Source point
end type

type csr_cache_struct
  type (csr_cache1_struct), allocatable :: cache1(:) ! Array of caches
end type

type csr_particle_bin_struct ! Structure for dividing a particle into bins
  integer ix0_bin                     ! index of first bin
  real(rp), allocatable :: charge(:)  ! how much charge in n^th bin
  real(rp), allocatable :: dcharge_ds(:)   
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

! sc_com%v(i) holds the parameters for the transverse space charge calculation
! at lattice element i.

type csr_common_struct
  real(rp) ds_track_step                   ! Tracking step size
  logical :: lcsr_component_on = .true.    ! Longitudinal csr component
  logical :: lsc_component_on = .true.     ! Longitudinal space charge component
  logical :: tsc_component_on = .true.     ! Transverse space charge component
end type

type (csr_common_struct), save, target :: csr_com

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)
!
! Routine to track a bunch of particles through and element with
! csr radiation effects.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   ele         -- Ele_struct: The element to track through.
!   param       -- Param_struct: General parameters.
!   bin         -- Csr_bin_struct: Binning information.
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
type (csr_cache_struct) cache
type (csr_bin_struct) bin

real(rp) ds_step, gamma, s_travel

integer i, j, nb, ix_ele, n_step

character(20) :: r_name = 'track1_bunch_sc'

! Init

ele => lat%ele_(ix_ele)
call convert_total_energy_to (ele%value(beam_energy$), lat%param%particle, gamma)

bunch_end = bunch_start

if (csr_com%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'DS_TRACK_STEP NOT SET!')
  call err_exit
else
  n_step = ele%value(l$) / csr_com%ds_track_step
endif

ds_step = ele%value(l$) / n_step

runt = ele
runt%value(l$) = ds_step
call attribute_bookkeeper (runt, lat%param)
bmad_com%auto_bookkeeper = .false.   

! Loop over tracking steps

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

  s_travel = ds_step * (i - 1)

  if (csr_com%lcsr_component_on .or. csr_com%lcsr_component_on) then
    call csr_bin_particles (bunch_end%particle, bin)    
    call csr_cache_setup (lat, ix_ele, s_travel, bin, cache)
  endif

  ! kick the particles

  do j = 1, size(bunch_end%particle)

    pt => bunch_end%particle(j)
    ! csr kick

    if (csr_com%lsc_component_on) call lsc_kick (bin, pt)
    if (csr_com%lcsr_component_on) call lcsr_kick (bin, pt, cache)

    ! track through element

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
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!   particle(:) -- Particle_struct: Array of particles
!
! Output:
!   bin(:)  -- Csr_bin_struct: The bins.
!   dbin(:) -- Csr_bin_struct: The bin differntial.
!-

subroutine csr_bin_particles (particle, bin)

implicit none

type (particle_struct), target :: particle(:)
type (particle_struct), pointer :: p
type (csr_bin_struct), target :: bin
type (csr_particle_bin_struct), save, allocatable :: pbin(:)
type (csr_bin1_struct), pointer :: b1

real(rp) z_min, z_max, f, dz_particle, dz, z_center
real(rp) zp_center, zp0, zp1, dz_bin2, zb0, zb1

integer i, j, np, np_bin, ix0, ib

! Init bins...
! z_min is the minimum z value at which there is charge.

np = size(particle)
dz_particle = bin%n_bin_particle * bin%dz_bin
z_min = minval(particle(:)%r%vec(5)) - dz_particle / 2.0
z_max = maxval(particle(:)%r%vec(5)) + dz_particle / 2.0

! renormalize so z_max - z_min = integer
! The left edge of bin%bin1(1) is at z_min
! The right edge of bin%bin1(n_bin) is at z_max

dz = int((z_max - z_min + 1) / bin%dz_bin) - (z_max - z_min) / bin%dz_bin
z_min = z_min - dz / 2
z_max = z_max + dz / 2

bin%n_bin = nint((z_max - z_min) / bin%dz_bin)

if (allocated(bin%bin1)) then
  if (size(bin%bin1, 1) < bin%n_bin) deallocate (bin%bin1)
endif

if (.not. allocated(bin%bin1)) allocate (bin%bin1(bin%n_bin))

do i = 1, bin%n_bin
  bin%bin1(i)%z1_edge = z_min + (i - 1) * bin%dz_bin
  bin%bin1(i)%z = bin%bin1(i)%z1_edge + bin%dz_bin / 2
  bin%bin1(i)%z2_edge = bin%bin1(i)%z1_edge + bin%dz_bin
enddo

! Init pbin structure...
! For each particle(i) there is a pbin(i) holding the information of
! how the particle is distributed in the bins.
! Each particle is distributed longitudinally in a triangular fashion

if (allocated(pbin)) then
  if (size(pbin) < size(particle) .or. &
                    size(pbin(1)%charge) < bin%n_bin_particle + 1) then
    do i = 1, size(pbin)
      deallocate (pbin(i)%charge, pbin(i)%dcharge_ds)
    enddo
    deallocate (pbin)
  endif
endif

if (.not. allocated(pbin)) then
  allocate (pbin(size(particle)))
  do i = 1, size(pbin)
    allocate (pbin(i)%charge(bin%n_bin_particle+1))
  enddo  
endif

do i = 1, np
  p => particle(i)
  if (.not. (p%ix_lost == not_lost$)) cycle
  zp_center = p%r%vec(5) ! center of particle
  zp0 = z_center - dz_particle / 2.0     ! particle left edge 
  zp1 = z_center + dz_particle / 2.0     ! particle right edge 
  ix0 = int(zp0 - z_min)                 ! left most bin index
  pbin(i)%ix0_bin = ix0
  f = 2.0 / dz_particle**2
  dz_bin2 = bin%dz_bin / 2
  do j = 0, bin%n_bin_particle
    zb0 = bin%bin1(j+ix0)%z1_edge
    zb1 = bin%bin1(j+ix0)%z2_edge   ! edges of the bin
    pbin(i)%charge(j+1) = charge_in_bin (zb0, zb1)
    pbin(i)%dcharge_ds = (charge_in_bin (zb0+dz_bin2, zb1+dz_bin2) - &
                    charge_in_bin (zb0-dz_bin2, zb1-dz_bin2)) / bin%dz_bin
  enddo
enddo

! Compute the particle distribution center in each bin

bin%bin1(:)%charge = 0
bin%bin1(:)%dcharge_ds = 0
bin%bin1(:)%x0 = 0
bin%bin1(:)%y0 = 0
bin%bin1(:)%sig_x = 0
bin%bin1(:)%sig_y = 0

do i = 1, np
  p => particle(i)
  if (.not. (p%ix_lost == not_lost$)) cycle
  do j = 0, bin%n_bin_particle
    ib = j + pbin(i)%ix0_bin
    b1 => bin%bin1(ib)
    b1%charge = b1%charge + pbin(i)%charge(j+1)
    b1%dcharge_ds = b1%dcharge_ds + pbin(i)%dcharge_ds(j+1)
    b1%x0 = b1%x0 + p%r%vec(1) * pbin(i)%charge(j+1)
    b1%y0 = b1%y0 + p%r%vec(3) * pbin(i)%charge(j+1)
  enddo
enddo

do ib = 1, bin%n_bin
  if (bin%bin1(ib)%charge == 0) cycle
  bin%bin1(ib)%x0 = bin%bin1(ib)%x0 / bin%bin1(ib)%charge
  bin%bin1(ib)%y0 = bin%bin1(ib)%y0 / bin%bin1(ib)%charge
enddo

! Compute the particle distribution sigmas in each bin
! Abs is used instead of the usual formula to lessen the effect
! of non-Gaussian tails

do i = 1, np
  p => particle(i)
  if (.not. (p%ix_lost == not_lost$)) cycle
  do j = 0, bin%n_bin_particle
    ib = j + pbin(i)%ix0_bin
    b1 => bin%bin1(ib)
    b1%sig_x = b1%sig_x + abs(p%r%vec(1) - b1%x0) * pbin(i)%charge(j+1)
    b1%sig_y = b1%sig_y + abs(p%r%vec(3) - b1%y0) * pbin(i)%charge(j+1)
  enddo
enddo

f = sqrt(pi/2)
do ib = 1, bin%n_bin
  b1 => bin%bin1(ib)
  if (b1%charge == 0) cycle
  b1%sig_x = f * b1%sig_x / b1%charge
  b1%sig_y = f * b1%sig_y / b1%charge
  b1%lsc_d0 = b1%sig_x * b1%sig_y
  b1%lsc_d1 = (b1%sig_x**2 + b1%sig_y**2) / (b1%sig_x + b1%sig_y)
enddo

!-------------------------------------------
contains

function charge_in_bin (z1_bin, z2_bin) result (charge)

real(rp) z1_bin, z2_bin, charge, z0, z1

! Integrate over left triangular half of particle distribution

z0 = max(zp0, zb0)        ! left integration edge
z1 = min(zp_center, zb1)  ! right integration edge
if (z1 > z0) then         ! If left particle half is in bin ...
  charge = f * p%charge * ((z1 - zb0)**2 - (z0 - zb0)**2)
else
  charge = 0
endif

! Integrate over right triangular half of particle distribution
z0 = max(zp_center, zb0)  ! left integration edge
z1 = min(zp1, zb1)        ! right integration edge
if (z1 > z0) then         ! If right particle half is in bin ...
  charge = charge + f * p%charge * ((zb1 - z0)**2 - (zb1 - z1)**2)
endif

end function

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine subroutine csr_cache_setup (lat, ix_ele, s_travel, bin, cache)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

subroutine csr_cache_setup (lat, ix_ele, s_travel, bin, cache)

implicit none

type (csr_bin_struct) bin
type (csr_cache_struct), target :: cache
type (ring_struct), target :: lat
type (ele_struct), pointer :: ele
type (csr_cache1_struct), pointer :: c1
type (csr_local_cache_struct) p

real(rp) c_const, s0, s_travel, s, s_boundary, d

integer i, ix_ele, n_ele

character(16) :: r_name = 'csr_cache_setup'

! The "P' element" is the element containing the source point P'.
! The "P element" is the element containing the kick point P.

n_ele = 0   
call next_element_calc (n_ele, p)
s = s_boundary + s_travel  ! s value at point P.

do i = 0, bin%n_bin

  c1 => cache%cache1(i)
  c1%z = i * bin%dz_bin

  ! first point is easy

  if (i == 0) then
    c1%phi = 0 
    c1%I_csr = 0
    c1%s_prime = s0 
    cycle
  endif

  ! calculate where we are

  do
    d = d_calc_csr(c1%z, p)
    c1%s_prime = s0 - d
    if (c1%s_prime > s_boundary) exit  ! If in element exit loop
    call next_element_calc (n_ele, p)
  enddo

  ! calculate csr
  
  c1%I_csr = I_csr (c1%z, d, p)

enddo

!-----------------------------------------------
contains

subroutine next_element_calc(n_ele, p)

type (csr_local_cache_struct) p

integer i, n_ele, ix_e

real(rp), save :: I0_csr
real(rp) dphi, zz, dd
real(rp) g_i(100), d_i(100)

!

n_ele = n_ele + 1          ! Number of elements including the P and P' elements
ix_e = ix_ele + 1 - n_ele  ! Index of current P' element
if (ix_e == 0) then
  call out_io (s_fatal$, r_name, 'CSR CALC EXTENDS BACK BEFORE BEGINNING OF LATTICE!')
  call err_exit
endif

s_boundary = lat%ele_(ix_e-1)%s  ! s value at beginning edge of the P' element

ele => lat%ele_(ix_e)  ! Pointer to the P' element

call convert_total_energy_to (ele%value(beam_energy$), lat%param%particle, p%gamma)
p%gamma2 = p%gamma**2

! calculate new values for d and g for this element and store in arrays

p%g = 0
if (ele%key == sbend$) p%g = ele%value(g$)
g_i(n_ele) = p%g

d_i(n_ele) = ele%value(l$)
if (n_ele == 1) d_i(n_ele) = s_travel

! calculate new v1, v3, etc.

p%v1 = 0; p%v3 = 0
p%w2 = 0
p%theta = 0

do i = 1, n_ele-1
  dphi = d_i(i) * g_i(i)
  p%v1 = p%v1 + d_i(i)
  p%v3 = p%v3 - d_i(i) * (3*p%theta**2 + 3*p%theta*dphi + dphi**2) / 6
  p%w2 = p%w2 + d_i(i) * (p%theta + dphi/2)
  p%theta = p%theta + dphi
enddo

! Calculate change in p%c_const at boundary

if (n_ele == 1) then
  p%c_const = 0
else
  zz = z_calc_csr (0.0_rp, p)
  p%c_const = p%c_const + (I0_csr - I_csr(zz, 0.0_rp, p))
endif

dd = d_i(n_ele)
zz = z_calc_csr (dd, p)
I0_csr = I_csr(zz, dd, p)   ! calculate for next time.

end subroutine

end subroutine
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function I_csr
!
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

function I_csr (z, d, p) result (i_this)

implicit none

type (csr_local_cache_struct) p

real(rp) z, d
real(rp) phi, t, a, k, i_this

!

  phi = p%g * d
  t = p%gamma * (d + p%v1)
  a = p%gamma2 * (p%w2 + phi * p%v1 + phi * d / 2)
  k = p%gamma * (p%theta + phi)
  I_this = -2 * p%gamma * (t + a * k) / (t**2 + a**2) + &
                                          1 / (p%gamma2 * z) + p%c_const

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function z_calc_csr
!
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

function z_calc_csr (d, p, dz_dd) result (z_this)

implicit none

type (csr_local_cache_struct) p
real(rp) d, phi, z_this, v1d, w2d
real(rp), optional :: dz_dd

!

phi = p%g * d
v1d = p%v1 + d
w2d = 2*p%w2 - phi*d
z_this = v1d / (2 * p%gamma2) - (p%v3 - phi**2 / 6 + w2d**2/(8*v1d))

if (present(dz_dd)) then
  dz_dd = 1 / (2 * p%gamma2) - (3*phi*p%g + phi*w2d/(4*v1d) - phi*w2d**2/(8*v1d**2))
endif

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function d_calc_csr
!
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

function d_calc_csr (z_target, p) result (d_this)

implicit none

type (csr_local_cache_struct) p
real(rp) z_target, z1, z2, dz_dd, dz1_dd, dz2_dd
real(rp) eps, d_this, delta_d_old, delta_d, z_this, d
real(rp), save :: d1 = 0, d2

integer j

character(8) :: r_name = 'd_calc_csr'

! Use Newton's method until root is bracketed.

eps = 1e-10 * abs(z_target) + 1e-14
z2 = 0

do 

  z1 = z_calc_csr (d1, p, dz1_dd) - z_target
  if (z1 * z2 < 0) exit  ! have bracket
  if (abs(z1) < eps) then
    d_this = d1
    return
  endif

  d2 = d1 - z1 / dz1_dd 
  z2 = z_calc_csr (d, p, dz2_dd) - z_target
  if (z1 * z2 < 0) exit  ! have bracket
  if (abs(z2) < eps) then
    d_this = d2
    return
  endif

  d1 = d2 - z2 / dz2_dd 

enddo

! Use rtsafe from Numerical Recipes

if (z1 > 0) then
  d_this = d1
  d1 = d2
  d2 = d_this
endif

d_this = (d1 + d2) / 2
delta_d_old = abs (d1 - d2)
delta_d     = delta_d_old
z_this = z_calc_csr(d_this, p, dz_dd) - z_target

do j = 1, 100

  ! Bisect if Newton out of range or not decre4asing fast enough

  if (((d_this-d2)*dz_dd-z_this) * ((d_this-d1)*dz_dd-z_this) > 0 .or. &
                              abs(2*z_this) > abs(delta_d_old * dz_dd)) then
    delta_d_old = delta_d
    delta_d = (d2 - d1) / 2
    d_this = d1 + delta_d
    if (d1 == d_this) return

  else
    delta_d_old = delta_d
    delta_d = d_this / dz_dd
    d = d_this
    d_this = d_this - delta_d
    if (d_this == d) return
  endif

  if (abs(z_this) < eps) return
  z_this = z_calc_csr(d_this, p, dz_dd) - z_target

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
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

subroutine lsc_kick (bin, particle)

implicit none

type (csr_bin_struct), target :: bin
type (csr_bin1_struct), pointer :: b1
type (particle_struct) particle

real(rp) x, y, sx, sy, a, b, c, zb_center, abs_z, sg, zz

integer i

!

do i = 1, bin%n_bin

  b1 => bin%bin1(i)

  x = particle%r%vec(1) - b1%x0
  y = particle%r%vec(3) - b1%y0
  abs_z = abs(particle%r%vec(5) - zb_center) 
  sg = bin%gamma * abs_z
  sx = b1%sig_x
  sy = b1%sig_y
  a = sx * sy * exp(((x/sx)**2 + (y/sy)**2) / 2)
  b = bin%gamma * (sx**2 + sy**2) / (sx + sy)
  c = bin%gamma2

!

  if (particle%r%vec(5) < b1%z1_edge) then
    particle%r%vec(6) = particle%r%vec(6) - b1%charge / (a + b * sg + c * sg**2)   
   
  elseif (particle%r%vec(5) > b1%z2_edge) then
    particle%r%vec(6) = particle%r%vec(6) + b1%charge / (a + b * sg + c * sg**2) 

  elseif (particle%r%vec(5) > zb_center) then
    zz = b1%z2_edge - particle%r%vec(5)

  else  ! particle%r%vec(5) < zb_center

  endif

end do

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lcsr_kick (bin, particle, cache)
!
! Routine to
!
! Modules needed:
!   use csr_mod
!
! Input:
!
! Output:
!
!-

subroutine lcsr_kick (bin, particle, cache)

implicit none

type (csr_bin_struct) bin
type (particle_struct) particle
type (csr_cache_struct) cache

real(rp) dz_bin, zp, r1, r2, dz, I_csr

integer i, ix_bin, i0, i_del, id

!

zp = particle%r%vec(5)
i0 = int((zp - bin%bin1(1)%z1_edge) / bin%dz_bin) + 1
i_del = nint((zp - bin%bin1(1)%z) / bin%dz_bin)
r2 = (zp - bin%bin1(1)%z) / bin%dz_bin - i_del

if (r2 > 0) then
  id = 1
else
  id = -1
  r2 = -r2
endif

r1 = 1 - r2

! loop over all bins

do i = i0, bin%n_bin
  dz = zp - bin%bin1(ix_bin)%z 
  I_csr = cache%cache1(i+i_del)%I_csr * r1 + cache%cache1(i+i_del+id)%I_csr * r2
  particle%r%vec(6) = particle%r%vec(6) + I_csr * bin%bin1(i)%dcharge_ds
enddo

end subroutine

end module
