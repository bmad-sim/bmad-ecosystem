!+
! Old version of CSR tracking.
!
! See the paper:
!   "An Efficient Formalism for Simulating Coherent Synchrotron Radiation"
!   D. Sagan
!-

module csr_old_mod

use beam_utils

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
  real(rp) dcharge_density_dz ! gradiant between this and preceeding bin
  real(rp) kick_csr    ! CSR kick
  real(rp) kick_lsc    ! LSC Kick.
end type

! Kicks in this structure are a function of the particle separation

type csr_kick1_struct ! Sub-structure for csr calculation cache
  real(rp) I_csr         ! Kick integral.
  real(rp) I_int_csr     ! Integrated Kick integral.
  real(rp) k_csr         ! Kick.
  real(rp) phi           ! Source point angle.
  real(rp) d             ! Distance between source point and end of element.
  real(rp) dz_particles  ! Distance between source and kicked particles.
  real(rp) s_prime       ! Source point location.
end type

type csr_bin_struct             ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) rel_mass             ! m_particle / m_electron
  real(rp) beta                 ! Relativistic beta factor.
  real(rp) :: dz_bin = 0        ! Bin width
  real(rp) ds_track_step        ! True step size
  real(rp) y2                   ! Height of source particle.
  real(rp) kick_factor          ! Coefficient to scale the kick
  integer particle              ! Particle type
  type (csr_bin1_struct), allocatable :: bin1(:)  
  type (csr_kick1_struct), allocatable :: kick1(:) ! Array of caches
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr_old (bunch_start, lat, ele, bunch_end, err, s_start, s_end)
!
! Routine to track a bunch of particles through an element with csr radiation effects.
!
! Input:
!   bunch_start -- Bunch_struct: Starting bunch position.
!   lat         -- lat_struct: Lattice.
!   ele         -- Ele_struct: The element to track through.
!   s_start     -- real(rp), optional: Starting position relative to ele. Default = 0
!   s_end       -- real(rp), optional: Ending position. Default is ele length.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. EG: Too many particles lost.
!-

subroutine track1_bunch_csr_old (bunch_start, lat, ele, bunch_end, err, s_start, s_end)

implicit none

type (lat_struct) lat
type (bunch_struct), target :: bunch_start, bunch_end
type (coord_struct), pointer :: pt
type (ele_struct) :: ele
type (ele_struct), save :: runt
type (ele_struct), pointer :: ele0
type (csr_bin_struct), save :: bin

real(rp), optional :: s_start, s_end
real(rp) s0_step
integer i, j, ns, nb, n_step, n_live, iu_wake

character(*), parameter :: r_name = 'track1_bunch_csr'
logical err, auto_bookkeeper

! Init

err = .true.

! No CSR for a zero length element.
! And taylor elements get ignored.

if (ele%value(l$) == 0 .or. ele%key == taylor$) then
  ele%csr_method = off$
  ele%space_charge_method = off$
  call track1_bunch_hom (bunch_end, ele, lat%param, bunch_end)
  err = .false.
  if (ele%key == taylor$ .and. ele%value(l$) == 0 .and. csr_param%print_taylor_warning) then
    call out_io (s_warn$, r_name, 'CSR calc for taylor element not done: ' // ele%name)
  endif
  return
endif

! n_step is the number of steps to take when tracking through the element.
! bin%ds_step is the true step length.

bunch_end = bunch_start

if (csr_param%n_bin <= csr_param%particle_bin_span + 1) then
  call out_io (s_fatal$, r_name, &
            'CSR_PARAM%N_BIN MUST BE GREATER THAN CSR_PARAM%PARTICLE_BIN_SPAN+1!')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (csr_param%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'CSR_PARAM%DS_TRACK_STEP NOT SET!')
  if (global_com%exit_on_error) call err_exit
  return
endif

! make sure that ele_len / track_step is an integer.

n_step = max (1, nint(ele%value(l$) / csr_param%ds_track_step))
bin%ds_track_step = ele%value(l$) / n_step

auto_bookkeeper = bmad_com%auto_bookkeeper ! save state
bmad_com%auto_bookkeeper = .false.   ! make things go faster

! Loop over the tracking steps
! runt is the element that is tracked through at each step.

do i = 0, n_step

  ! track through the runt

  if (i /= 0) then
    call create_uniform_element_slice (ele, lat%param, i, n_step, runt, s_start, s_end)
    runt%csr_method = off$
    runt%space_charge_method = off$
    call track1_bunch_hom (bunch_end, runt, lat%param, bunch_end)
  endif

  s0_step = i * bin%ds_track_step
  if (present(s_start)) s0_step = s0_step + s_start

  ! Cannot do a realistic calculation if there are less particles than bins

  n_live = count(bunch_end%particle%state == alive$)
  if (n_live < csr_param%n_bin) then
    call out_io (s_error$, r_name, 'NUMBER OF LIVE PARTICLES: \i0\ ', &
                          'LESS THAN NUMBER OF BINS FOR CSR CALC.', &
                          'AT ELEMENT: ' // trim(ele%name) // '  [# \i0\] ', &
                          i_array = [n_live, ele%ix_ele ])
    return
  endif

  call csr_bin_particles (ele, bunch_end%particle, bin)    

  ! ns = 0 is the unshielded kick.
  ! For the shielding image currents never use the small angle approximation

  bin%bin1(:)%kick_csr = 0
  bin%bin1(:)%kick_lsc = 0

  do ns = 0, csr_param%n_shield_images

    ! %kick_factor takes into account that at the endpoints we are only putting in a half kick.

    bin%kick_factor = 1
    if (i == 0 .or. i == n_step) bin%kick_factor = 0.5

    bin%y2 = ns * csr_param%beam_chamber_height

    if (ns == 0) then
      call csr_bin_kicks (lat, ele, s0_step, bin, csr_param%small_angle_approx)

    else
      ! The factor of two is due to there being image currents both above and below.
      ! The factor of -1^ns accounts for the sign of the image currents
      bin%kick_factor = bin%kick_factor * 2 * (-1)**ns
      call csr_bin_kicks (lat, ele, s0_step, bin, .false.)
    endif

  enddo

  ! loop over all particles and give them a kick

  do j = 1, size(bunch_end%particle)
    if (bunch_end%particle(j)%state /= alive$) cycle
    call csr_kick_calc (ele, bin, bunch_end%particle(j))
  enddo

  call save_bunch_track (bunch_end, ele, s0_step)

  ! Record to file?

  if (csr_param%write_csr_wake) then
    iu_wake = lunget()
    open (iu_wake, file = 'old_csr_wake.dat', access = 'append')
    if (i == 0) then
      write (iu_wake, '(a)') '!------------------------------------------------------------'
      write (iu_wake, '(a, i6, 2x, a)') '! ', ele%ix_ele, trim(ele%name)
    endif
    write (iu_wake, '(a)') '!#-----------------------------'
    write (iu_wake, '(a, i4, f12.6)') '! ', i, ele%s_start + s0_step
    write (iu_wake, '(a)') '!         Z   Charge/Meter     CSR_Wake' 
    do j = 1, csr_param%n_bin
      write (iu_wake, '(f12.6, 2es14.6)') bin%bin1(j)%z_center, &
                                  bin%bin1(j)%charge/bin%dz_bin, bin%bin1(j)%kick_csr
    enddo
    close (iu_wake)
  endif

enddo

bmad_com%auto_bookkeeper = auto_bookkeeper  ! restore state
err = .false.

end subroutine track1_bunch_csr_old

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_particles (ele, particle, bin)
!
! Routine to bin the particles longitudinally in s. 
!
! To avoid noise in the cacluation, every particle is considered to have a 
! triangular distribution with a base length  given by 
!   csr_param%particle_bin_span * bin%dz_bin. 
! That is, particles will, in general, overlap multiple bins. 
!
! Input:
!   ele                  -- ele_struct: Element being tracked through.
!   particle(:)          -- Coord_struct: Array of particles
!   csr_param            -- Csr_parameter_struct: CSR common block (not an argument).
!     %n_bin             -- Number of bins.
!     %particle_bin_span -- Particle length / dz_bin. 
!
! Output:
!   bin           -- Csr_bin_struct: The bin structure.
!     %dz_bin     -- Bin longitudinal length
!     %bin1(1:) -- Array of bins.
!-

subroutine csr_bin_particles (ele, particle, bin)

implicit none

type this_local_struct   ! Temporary structure 
  real(rp) :: charge     ! how much charge of particle in bin
  real(rp) :: x0, y0     ! particle center
  integer ib             ! bin index
end type

type (ele_struct) ele
type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p
type (csr_bin_struct), target :: bin
type (this_local_struct), allocatable :: tloc(:)
type (csr_bin1_struct), pointer :: bin1

real(rp) z_center, z_min, z_max, f, dz_particle, dz, z_maxval, z_minval
real(rp) zp_center, zp0, zp1, zb0, zb1, charge

integer i, j, n, ix0, ib, ic

character(20) :: r_name = 'csr_bin_particles'

! Init bins...
! The left edge of bin%bin1(1) is at z_min
! The right edge of bin%bin1(n_bin) is at z_max
! The first and last bins are empty.

if (ele%space_charge_method /= slice$ .and. ele%csr_method /= one_dim$) return

z_maxval = maxval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
z_minval = minval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
dz = z_maxval - z_minval
bin%dz_bin = dz / (csr_param%n_bin - 2 - (csr_param%particle_bin_span + 1))
bin%dz_bin = 1.0000001 * bin%dz_bin     ! to prevent round off problems
z_center = (z_maxval + z_minval) / 2
z_min = z_center - csr_param%n_bin * bin%dz_bin / 2
z_max = z_center + csr_param%n_bin * bin%dz_bin / 2
dz_particle = csr_param%particle_bin_span * bin%dz_bin

! allocate memeory for the bins

if (allocated(bin%bin1)) then
  if (size(bin%bin1, 1) < csr_param%n_bin) deallocate (bin%bin1, bin%kick1)
endif

if (.not. allocated(bin%bin1)) allocate (bin%bin1(csr_param%n_bin), bin%kick1(-csr_param%n_bin:csr_param%n_bin))

! Fill in some z information

do i = 1, csr_param%n_bin
  bin%bin1(i)%z0_edge  = z_min + (i - 1) * bin%dz_bin
  bin%bin1(i)%z_center = bin%bin1(i)%z0_edge + bin%dz_bin / 2
  bin%bin1(i)%z1_edge  = bin%bin1(i)%z0_edge + bin%dz_bin
enddo

! Init the tloc structure...
! Each particle is distributed longitudinally in a triangular fashion.
! The tloc records how much charge for a given particle is in a bin.

n = size(particle) * (csr_param%particle_bin_span + 2)
allocate (tloc(n))
tloc%ib = -1

! Compute the particle distribution center in each bin

bin%bin1(:)%charge = 0
bin%bin1(:)%x0 = 0
bin%bin1(:)%y0 = 0
bin%bin1(:)%sig_x = 0
bin%bin1(:)%sig_y = 0
bin%bin1(:)%dcharge_density_dz = 0

f = 2.0 / dz_particle**2

! The contribution to the charge in a bin from a particle is computed from the overlap
! between the particle and the bin.
 
ic = 0
do i = 1, size(particle)
  p => particle(i)
  if (p%state /= alive$) cycle
  zp_center = p%vec(5) ! center of particle
  zp0 = zp_center - dz_particle / 2       ! particle left edge 
  zp1 = zp_center + dz_particle / 2       ! particle right edge 
  ix0 = nint((zp0 - z_min) / bin%dz_bin)  ! left most bin index
  do j = 0, csr_param%particle_bin_span+1
    ib = j + ix0
    bin1 => bin%bin1(ib)
    zb0 = bin%bin1(ib)%z0_edge
    zb1 = bin%bin1(ib)%z1_edge   ! edges of the bin
    charge = charge_in_bin (zb0, zb1)
    bin1%charge = bin1%charge + charge
    bin1%x0 = bin1%x0 + p%vec(1) * charge
    bin1%y0 = bin1%y0 + p%vec(3) * charge
    ic = ic + 1
    tloc(ic)%charge = charge
    tloc(ic)%x0 = p%vec(1)
    tloc(ic)%y0 = p%vec(3)
    tloc(ic)%ib = ib
  enddo
enddo

do ib = 1, csr_param%n_bin
  if (ib /= 1) bin%bin1(ib)%dcharge_density_dz = (bin%bin1(ib)%charge - bin%bin1(ib-1)%charge) / bin%dz_bin**2
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
do ib = 1, csr_param%n_bin
  bin1 => bin%bin1(ib)
  if (bin1%charge == 0) cycle
  bin1%sig_x = f * bin1%sig_x / bin1%charge
  bin1%sig_y = f * bin1%sig_y / bin1%charge
  bin1%lsc_d0 = bin1%sig_x * bin1%sig_y
  if (bin1%sig_x == 0 .and. bin1%sig_y == 0) then
    bin1%lsc_d1 = 0
  else
    bin1%lsc_d1 = (bin1%sig_x**2 + bin1%sig_y**2) / (bin1%sig_x + bin1%sig_y)
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

end subroutine csr_bin_particles

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_kicks (lat, ele, s_travel, bin, small_anlge_approx)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Input:
!   lat       -- lat_struct: Lattice.
!   ele       -- Element_struct: Element to set up cache for.
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

subroutine csr_bin_kicks (lat, ele, s_travel, bin, small_angle_approx)

implicit none

type (csr_bin_struct), target :: bin
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct) ele
type (csr_kick1_struct), pointer :: kick1
type (csr_kick_factor_struct) k_factor

real(rp) s_travel, s_kick, s0_kick_ele, coef, e_tot, f1
real(rp), allocatable :: g_i(:), d_i(:)

integer i, n_ele_pp, n_bin

logical small_angle_approx

character(16) :: r_name = 'csr_bin_kicks'

! Assume a linear energy gain

branch => lat%branch(ele%ix_branch)
f1 = s_travel / ele%value(l$)
e_tot = f1 * branch%ele(ele%ix_ele-1)%value(e_tot$) + (1 - f1) * ele%value(e_tot$)
call convert_total_energy_to (e_tot, branch%param%particle, bin%gamma, beta = bin%beta)
bin%gamma2 = bin%gamma**2
bin%rel_mass = mass_of(lat%param%particle) / m_electron 
bin%particle = lat%param%particle

! The kick point P is fixed.
! The source point P' varies from bin to bin.
! n_ele_pp is the number of elements between P' and P excluding the 
! elements containing P and P'. n_ele_pp will be incremented by
! next_element_params_calc as we go from one bin to the next.

n_ele_pp = -1
call next_element_params_calc (n_ele_pp, s0_kick_ele, k_factor)
s_kick = s0_kick_ele + s_travel  ! absolute s value at point P.

! Loop over all kick1 bins and compute the kick or kick integral.
! The loop steps in increasing dz since that is what next_element_params_calc expects.

do i = lbound(bin%kick1, 1), ubound(bin%kick1, 1)

  kick1 => bin%kick1(i)
  kick1%dz_particles = i * bin%dz_bin

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
    call I_csr (kick1, i, k_factor, bin)
    if (bin%kick1(i)%I_int_csr == 0 .and. i /= lbound(bin%kick1, 1)) then
      bin%kick1(i)%I_int_csr = (bin%kick1(i)%I_csr + bin%kick1(i-1)%I_csr) * bin%dz_bin / 2
    endif
  else
    call kick_image_charge (kick1, k_factor, bin)
  endif

enddo

! 

coef = bin%ds_track_step * r_e / (bin%rel_mass * e_charge * abs(charge_of(lat%param%particle)) * bin%gamma)
n_bin = csr_param%n_bin

! CSR & Image charge kick

if (bin%y2 == 0) then
  if (ele%space_charge_method == slice$) then  ! Longitudinal component
    do i = 1, n_bin
      bin%bin1(i)%kick_csr = coef * dot_product(bin%kick1(i:1:-1)%I_int_csr, bin%bin1(1:i)%dcharge_density_dz)
    enddo
  endif

else  ! Image charge
  do i = 1, n_bin
    bin%bin1(i)%kick_csr = bin%bin1(i)%kick_csr + coef * &
                  dot_product(bin%kick1(i-1:i-n_bin:-1)%k_csr, bin%bin1(1:n_bin)%charge)
  enddo
endif

! Longitudinal space charge kick

if (ele%space_charge_method == slice$) then  
  if (bin%y2 == 0) then
    call lsc_y0_kick_calc (ele, bin)
  endif
endif

!----------------------------------------------------------------------------
contains

subroutine next_element_params_calc (n_ele_pp, s0_kick_ele, k_factor)

type (csr_kick_factor_struct) k_factor
type (ele_struct), pointer :: source_ele

integer i, n_ele_pp, ix_source

real(rp) phi, dphi, s0_kick_ele

! n_ele_pp is the number of elements between P' and P excluding the 
! P and P' elements
! Assume a drift before the first element if needed.

n_ele_pp = n_ele_pp + 1
ix_source = ele%ix_ele - n_ele_pp  ! Index of current P' element

source_ele => branch%ele(ix_source)  ! Pointer to the P' element

! Assume a drift before the first element if needed.

if (ix_source == 0) then
  s0_kick_ele = -1d20  ! something large and negative
else
  s0_kick_ele = branch%ele(ix_source-1)%s  ! s value at beginning edge of the P' element
endif

! calculate new values for d and g for this element and store in arrays
! d_i(i) is the length of the ith element
! g_i(i) is the bending radius of the ith element

k_factor%g = 0
if (source_ele%key == sbend$) k_factor%g = source_ele%value(g$)

if (.not. allocated(g_i)) allocate(g_i(n_ele_pp+100), d_i(n_ele_pp+100))
if (size(g_i) <= n_ele_pp) then
  call re_allocate (g_i, n_ele_pp + 100)
  call re_allocate (d_i, n_ele_pp + 100)
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
    k_factor%v3 = k_factor%v3 + d_i(i) * (k_factor%theta**2 + k_factor%theta*dphi + dphi**2 / 3) / 2
    k_factor%w2 = k_factor%w2 + d_i(i) * (k_factor%theta + dphi/2)
  else
    phi = k_factor%theta 
    if (g_i(i) == 0) then
      k_factor%v = k_factor%v + d_i(i) * cos(phi)
      k_factor%w2 = k_factor%w2 + d_i(i) * sin(phi)
    else
      k_factor%v = k_factor%v + (sin(phi + dphi) - sin(phi)) / g_i(i)
      k_factor%w2 = k_factor%w2 + (cos(phi) - cos(phi + dphi)) / g_i(i)
    endif
    k_factor%v3 = k_factor%v1 - k_factor%v
  endif

  k_factor%theta = k_factor%theta + dphi

enddo

end subroutine next_element_params_calc

end subroutine csr_bin_kicks
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_y0_kick_calc (ele, bin)
!
! Routine to cache intermediate values needed for the lsc calculation.
! This routine is not for image currents.
!
! Input:
!   ele       -- ele_struct: Element being tracked through.
!   bin       -- csr_bin_struct: Binned particle averages.
!     %bin1(:)  -- bin array of particle averages.
!
! Output:
!   bin     -- csr_bin_struct: Binned particle averages.
!     %bin1(:)%kick_lsc -- Integrated kick.
!-

subroutine lsc_y0_kick_calc (ele, bin)

implicit none

type (ele_struct) ele
type (csr_bin_struct), target :: bin
type (csr_bin1_struct), pointer :: bin1

real(rp) sx, sy, a, b, c, dz, factor, sig_x_ave, sig_y_ave, charge_tot

integer i, j

character(16) :: r_name = 'lsc_bin_kicks'

! If there are too few particles in a bin the sigma calc may give a bad value.
! This can be a problem if the computed sigma is small.
! Therefore  ignore any bins with a very small sigma. 
! To know what is "small" is, compute the average sigma

charge_tot = sum(bin%bin1(:)%charge)
if (charge_tot == 0) return

sig_x_ave = dot_product(bin%bin1(:)%sig_x, bin%bin1(:)%charge) / charge_tot
sig_y_ave = dot_product(bin%bin1(:)%sig_y, bin%bin1(:)%charge) / charge_tot

if (sig_y_ave == 0 .or. sig_x_ave == 0) return  ! Symptom of not enough particles.

! Compute the kick at the center of each bin

bin%bin1(:)%kick_lsc = 0
if (ele%space_charge_method /= slice$) return

do i = 1, csr_param%n_bin
  bin1 => bin%bin1(i)
  sx = bin1%sig_x
  sy = bin1%sig_y
  if (sx < sig_x_ave * csr_param%sigma_cutoff .or. sy < sig_y_ave * csr_param%sigma_cutoff) then
    bin1%sig_x = 0  ! Mark for tsc calc.
    bin1%sig_y = 0
    cycle
  endif
  a = sx * sy
  b = bin%gamma * (sx**2 + sy**2) / (sx + sy)
  c = bin%gamma**2

  do j = 1, csr_param%n_bin
    if (i == j) cycle
    dz = bin%bin1(j)%z_center - bin%bin1(i)%z_center
    bin%bin1(j)%kick_lsc = bin%bin1(j)%kick_lsc + bin1%charge * sign(1.0_rp, dz) / (a + b * abs(dz) + c * dz**2)
  enddo

enddo

factor = bin%kick_factor * bin%ds_track_step * r_e / &
          (bin%rel_mass * e_charge * abs(charge_of(bin%particle)) * bin%gamma)
bin%bin1(:)%kick_lsc = factor * bin%bin1(:)%kick_lsc

end subroutine lsc_y0_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine I_csr (kick1, i_bin, k_factor, bin) 
!
! Routine to calculate the CSR kick integral.
!
! Input:
!   kick1    -- Csr_kick1_struct: 
!     %dz_particles -- Real(rp): Distance between source and kicked particles.
!     %d            -- Real(rp): Distance between source and kick points.
!   i_bin    -- Integer: Bin index.
!   k_factor -- Csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!
! Output:
!   kick1 -- Csr_kick1_struct: 
!     %I_csr     -- Real(rp): CSR kick integral.
!     %I_int_csr -- Real(rp): Integral of I_csr. Only calculated for i_bin = 1 since
!                     it is not needed otherwise.
!-

subroutine I_csr (kick1, i_bin, k_factor, bin)

implicit none

type (csr_kick1_struct) kick1
type (csr_kick_factor_struct) k_factor
type (csr_bin_struct) bin

real(rp) z, d
real(rp) phi, t, a, k, gam, gam2, g

integer i_bin

! 

z = kick1%dz_particles
d = kick1%d

kick1%I_int_csr = 0

if (z <= 0) then
  kick1%I_csr = 0
  return
endif

gam = bin%gamma
gam2 = bin%gamma2
g = k_factor%g 
phi = g * d
t = d + k_factor%v1
a = gam * (k_factor%w2 + phi * k_factor%v1 + phi * d / 2)
k = gam * (k_factor%theta + phi)

kick1%I_csr = bin%kick_factor * (-2 * (t + a * k) / (t**2 + a**2) + 1 / (gam2 * z))
if (i_bin == 1 .and.k_factor%v1 == 0) then
  kick1%I_int_csr = bin%kick_factor * (-(d*g)**2 / 4 + log(2*gam2 * z / d) / gam2)
endif

end subroutine I_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kick_image_charge (kick1, k_factor, bin) 
!
! Routine to calculate the image charge kick.
!
! Input:
!   kick1    -- Csr_kick1_struct: 
!     %dz_particles -- Real(rp): Distance between source and kicked particles.
!     %d            -- Real(rp): Distance between source and kick points.
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!
! Output:
!   kick1    -- Csr_kick1_struct: 
!     %k_csr -- Real(rp): Image charge kick.
!-

subroutine kick_image_charge (kick1, k_factor, bin)

implicit none

type (csr_kick1_struct) kick1
type (csr_kick_factor_struct), target :: k_factor
type (csr_kick_factor_struct), pointer :: kf
type (csr_bin_struct) bin

real(rp) z, d
real(rp) N_vec(3), G_vec(3), B_vec(3), Bp_vec(3), NBp_vec(3), NBpG_vec(3), rad_cross_vec(3)
real(rp) phi, sin_phi, cos_phi, OneNBp, OneNBp3, radiate, coulomb1

!

kick1%k_csr = 0

z = kick1%dz_particles
d = kick1%d

kf => k_factor
phi = kf%g * d
sin_phi = sin(phi)
cos_phi = cos(phi)

N_vec = k_factor%L_vec / kf%L
B_vec = [cos(kf%theta), sin(kf%theta), 0.0_rp ]
Bp_vec = bin%beta * [cos_phi, -sin_phi, 0.0_rp ]
G_vec = bin%beta**2 * kf%g * [sin_phi, cos_phi, 0.0_rp ]

OneNBp = 1 - sum(N_vec * Bp_vec)
OneNBp3 = OneNBp**3

NBp_vec = N_vec - Bp_vec
NBpG_vec = cross(NBp_vec, G_vec)
rad_cross_vec = cross(N_vec, NBpG_vec)

radiate  = dot_product (B_vec, rad_cross_vec) / (kf%L * OneNBp3)
coulomb1 = dot_product (B_vec, NBp_vec) / (bin%gamma2 * kf%L**2 * OneNBp3)
kick1%k_csr = bin%kick_factor * (radiate + coulomb1)

!-----------------------------------------------------------------------------
contains
function cross (a_vec, b_vec) result (c_vec)

real(rp), intent(in) :: a_vec(3), b_vec(3)
real(rp) :: c_vec(3)

c_vec(1) = a_vec(2) * b_vec(3) - a_vec(3) * b_vec(2)
c_vec(2) = a_vec(3) * b_vec(1) - a_vec(1) * b_vec(3)
c_vec(3) = a_vec(1) * b_vec(2) - a_vec(2) * b_vec(1)

end function

end subroutine kick_image_charge

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function d_calc_csr (dz_particles, k_factor, bin, small_angle_approx) result (d_this)
!
! Routine to calculate the distance between source and kick points.
!
! Input:
!   dz_particles -- Real(rp): Distance between source and kicked particles
!   k_factor     -- csr_kick_factor_struct: Geometric values
!   bin          -- Csr_bin_struct:
!   small_angle_approx -- Logical: If True then use a small angle approximation.
!
! Output:
!   d_this   -- Real(rp): Distance between source and kick points.
!   k_factor -- csr_kick_factor_struct: Geometric values
!     %l        -- L vector length.
!     %l_vec(3) -- L vector components.
!-

function d_calc_csr (dz_particles, k_factor, bin, small_angle_approx) result (d_this)

implicit none

type (csr_kick_factor_struct) k_factor
type (csr_bin_struct) bin

real(rp) dz_particles, z1, z2, dz_dd, dz1_dd, dz2_dd
real(rp) eps_z, eps_d, d_this, delta_d_old, delta_d, dz_calc, d_old
real(rp) d1, d2, a, b, c, d_plus_v1

integer j

logical small_angle_approx

character(8) :: r_name = 'd_calc_csr'

! If not in a bend then can do the calculation exactly

if (k_factor%g == 0) then
  a = 1 / bin%gamma2
  b = 2 * (k_factor%v3 - dz_particles)
  c = - (k_factor%w2**2 + bin%y2**2)
  d_this = (-b + sqrt(b**2 - 4 * a * c)) / (2 * a) - k_factor%v1
  k_factor%l_vec(1) = d_this + k_factor%v
  k_factor%l_vec(2) = k_factor%w2
  k_factor%l_vec(3) = bin%y2
  k_factor%L = sqrt(dot_product(k_factor%L_vec, k_factor%L_vec))
  return
endif

! Use Newton's method until root is bracketed.
! Initial d1 is just an approximate guess to get in the ball park.

eps_z = 1d-10 * abs(dz_particles) + 1d-14
eps_d = 1d-10 * (1 + abs(k_factor%v))

z2 = 0
if (bin%y2 /= 0) then
  d1 = sqrt(3 * bin%y2 / abs(k_factor%g))
elseif (dz_particles >= 0) then
  d1 = min(2 * bin%gamma2 * dz_particles, (6 * dz_particles / k_factor%g**2) ** (0.3333))
else
  d1 = (bin%y2**2 + dz_particles**2) / (2 * dz_particles)
endif

! Now bracket the root

do 

  z1 = z_calc_csr (d1, k_factor, bin, small_angle_approx, dz1_dd) - dz_particles
  if (abs(z1) < eps_z) then
    d_this = d1
    return
  endif
  if (z1 * z2 < 0) exit  ! have bracket

  d2 = d1 - z1 / dz1_dd 
  z2 = z_calc_csr (d2, k_factor, bin, small_angle_approx, dz2_dd) - dz_particles
  if (abs(z2) < eps_z) then
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

  !! if (abs(dz_calc) < eps_z) return
  if (abs(d2-d1) < eps_d) return
  dz_calc = z_calc_csr(d_this, k_factor, bin, small_angle_approx, dz_dd) - dz_particles

  if (dz_calc < 0) then
    d1 = d_this
  else
    d2 = d_this
  endif

enddo

call out_io (s_abort$, r_name, 'ALGORITHM NOT CONVERGING.')
if (global_com%exit_on_error) call err_exit

end function d_calc_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function z_calc_csr (d, k_factor, bin, small_angle_approx, dz_dd) result (z_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle.
!
! Input:
!   d        -- Real(rp): Distance between the source point and the kick point.
!   k_factor -- csr_kick_factor_struct: Other parameters needed in the calculation.
!   bin      -- Csr_bin_struct:
!   small_angle_approx -- Logical: If True then use a small angle approximation.
!
! Output:
!   k_factor -- csr_kick_factor_struct: Geometric values
!     %l        -- L vector length.
!     %l_vec(3) -- L vector components.
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

character(*), parameter :: r_name = 'z_calc_csr'

! Special cases

if (k_factor%v1 == 0 .and. d == 0 .and. bin%y2 == 0) then
  z_this = 0
  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2)
  return
endif

if (bin%y2 == 0 .and. d < 0) then
  z_this = 2 * d
  if (present(dz_dd)) dz_dd = 2
  return
endif

!

phi = k_factor%g * d
v1d = k_factor%v1 + d
kf => k_factor

! General case with small angle approx

if (small_angle_approx) then
  if (v1d == 0 .or. bin%gamma2 == 0) then
    call out_io (s_fatal$, r_name, 'ERROR IN CSR CALC.') 
    if (global_com%exit_on_error) call err_exit
  endif

  w2d = 2*k_factor%w2 - phi*d
  y22 = 4 * bin%y2**2
  z_this = v1d / (2 * bin%gamma2) + (k_factor%v3 + phi**2 * d / 6 - (w2d**2 + y22)/(8*v1d))

  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2) + (phi**2/2 + phi*w2d/(2*v1d) + (w2d**2 + y22)/(8*v1d**2))

! General case without small angle approx

else
  if (abs(phi) < 1d-2) then
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
  kf%L = sqrt(dot_product(kf%L_vec, kf%L_vec))

  z_this = v1d / (2 * bin%gamma2) + (v1d - kf%L)

  if (present(dz_dd)) dz_dd = 1 / (2 * bin%gamma2) + 1 - (kf%L_vec(1) * cos(phi) - kf%L_vec(2) * sin(phi)) / kf%L
                

endif

end function z_calc_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_kick_calc (ele, bin, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
!   ele      -- ele_struct: Element being tracked through.
!   bin      -- Csr_bin1_struct: Binned beam
!   particle -- Coord_struct: Particle to kick.
!     %vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Coord_struct: Particle to kick.
!     %vec(6) -- Final particle energy.
!-

subroutine csr_kick_calc (ele, bin, particle)

implicit none

type (ele_struct) ele
type (csr_bin_struct), target :: bin
type (coord_struct), target :: particle
type (csr_bin1_struct), pointer :: bin1

real(rp) zp, r1, r0, dz, dpz, kx, ky, f0, f, beta0
real(rp), pointer :: vec(:)
integer i, i0, i_del

! We use a weighted average between %kick1(j)%I_csr and %kick1(j+1)%I_csr
! so that the integral varies smoothly as a function of particle%vec(5).

if (.not. allocated(bin%bin1)) return  ! True if kicks are turned off.

zp = particle%vec(5)
i0 = int((zp - bin%bin1(1)%z_center) / bin%dz_bin) + 1
r1 = (zp - bin%bin1(i0)%z_center) / bin%dz_bin
r0 = 1 - r1
vec => particle%vec

if (r1 < 0 .or. r1 > 1 .or. i0 < 1 .or. i0 >= csr_param%n_bin) then
  print *, 'CSR INTERNAL ERROR!'
  if (global_com%exit_on_error) call err_exit
endif

vec(6) = vec(6) + r0 * bin%bin1(i0)%kick_csr + r1 * bin%bin1(i0+1)%kick_csr

! Longitudinal space charge

if (ele%space_charge_method == slice$) then
  vec(6) = vec(6) + r0 * bin%bin1(i0)%kick_lsc + r1 * bin%bin1(i0+1)%kick_lsc
endif

! Must update beta and z due to the energy change

beta0 = particle%beta
call convert_pc_to ((1+vec(6))* particle%p0c, particle%species, beta = particle%beta)
vec(5) = vec(5) * particle%beta / beta0

! Transverse space charge.

if (ele%space_charge_method == slice$) then
  f0 = bin%kick_factor * bin%ds_track_step * r_e / (twopi * &
           bin%dz_bin * bin%rel_mass * e_charge * abs(charge_of(bin%particle)) * bin%gamma**3)

  bin1 => bin%bin1(i0)
  if (bin1%sig_x /= 0) then
    call bbi_kick ((vec(1)-bin1%x0)/bin1%sig_x, (vec(3)-bin1%y0)/bin1%sig_y, bin1%sig_y/bin1%sig_x, kx, ky)
    f = f0 * r0 * bin1%charge / (bin1%sig_x + bin1%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    vec(2) = vec(2) - kx * f
    vec(4) = vec(4) - ky * f
  endif

  bin1 => bin%bin1(i0+1)
  if (bin1%sig_x /= 0) then
    call bbi_kick ((vec(1)-bin1%x0)/bin1%sig_x, (vec(3)-bin1%y0)/bin1%sig_y, bin1%sig_y/bin1%sig_x, kx, ky)
    f = f0 * r1 * bin1%charge / (bin1%sig_x + bin1%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    vec(2) = vec(2) - kx * f   
    vec(4) = vec(4) - ky * f
  endif

endif

end subroutine csr_kick_calc

end module
