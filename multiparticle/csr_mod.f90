!+
! See the paper:
!   "An Efficient Formalism for Simulating Coherent Synchrotron Radiation"
!   D. Sagan
!-

module csr_mod

use make_mat6_mod
use beam_utils
use spline_mod
use nr, only: zbrent

type csr_source_ele_struct
  type (ele_struct), pointer :: ele     ! Source ele
  type (coord_struct) orbit_c           ! centroid Orbit at exit end
  type (floor_position_struct) floor_c  ! Floor position of centroid at exit end
  type (spline_struct) spline           ! spline for centroid orbit. spline%x = s-position with x = 0 the beginning of the element.
end type

type csr_bunch_slice_struct   ! Structure for a single particle bin.
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
  real(rp) L             ! Distance between source and kick points.
  real(rp) dz_particles  ! Distance between source and kicked particles at constant time.
  real(rp) s_source      ! source point location.
  real(rp) z_source      ! distance between source and referece particle.
  real(rp) theta_L       ! Angle of L vector
  real(rp) theta_sl      ! Angle between velocity of particle at source pt and L
  real(rp) theta_tot     ! Angle between velocity of kicked particle and velocoty of source particle
  real(rp) g_bend        ! Source point 1/bending_radius
  integer ix_ele_source  ! source element index.
  type (floor_position_struct) floor_s  ! Floor position of source pt
end type

type csr_top_level_struct             ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) rel_mass             ! m_particle / m_electron
  real(rp) beta                 ! Relativistic beta factor.
  real(rp) :: dz_slice = 0      ! Bin width
  real(rp) ds_track_step        ! True step size
  real(rp) s_kick               ! Kick point longitudinal location
  real(rp) z_kick               ! distance between kick and referece particle.
  real(rp) y2                   ! Height of source particle.
  real(rp) kick_factor          ! Coefficient to scale the kick
  logical small_angle_approx
  type(floor_position_struct) floor_k   ! Floor coords at kick point
  integer ix_ele_kick           ! kicked element index.
  integer particle              ! Particle type
  type (csr_bunch_slice_struct), allocatable :: slice(:)    ! slice(i) refers to the i^th bunch slice.
  type (csr_kick1_struct), allocatable :: kick1(:)          ! kick1(i) referes to the kick between two slices i bins apart.
  type (csr_source_ele_struct), allocatable :: source_ele(:)
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch_start, ele, bunch_end, err, s_start, s_end, centroid)
!
! Routine to track a bunch of particles through an element with csr radiation effects.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bunch_start -- Bunch_struct: Starting bunch position.
!   ele         -- Ele_struct: The element to track through. Must be part of a lattice.
!   s_start     -- real(rp), optional: Starting position relative to ele. Default = 0
!   s_end       -- real(rp), optional: Ending position. Default is ele length.
!   centroid(0:) -- coord_struct, optional: Centroid orbit. Only needed if the
!                     central orbit is far from the zero orbit.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. EG: Too many particles lost.
!-

subroutine track1_bunch_csr (bunch_start, ele, bunch_end, err, s_start, s_end, centroid)

implicit none

type (bunch_struct), target :: bunch_start, bunch_end
type (coord_struct), pointer :: pt
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (ele_struct), save :: runt
type (ele_struct), pointer :: ele0, s_ele
type (csr_top_level_struct), target :: csr_top
type (csr_source_ele_struct), pointer :: source_ele
type (coord_struct), optional :: centroid(0:)
type (floor_position_struct) floor

real(rp), optional :: s_start, s_end
real(rp) s0_step, vec0(6), vec(6)

integer i, j, ie, ns, nb, n_step, n_live

character(*), parameter :: r_name = 'track1_bunch_csr'
logical err, auto_bookkeeper

! Init

err = .true.
branch => ele%branch
csr_top%small_angle_approx = csr_param%small_angle_approx

! No CSR for a zero length element.
! And taylor elements get ignored.

if (ele%value(l$) == 0 .or. ele%key == taylor$) then
  ele%csr_calc_on = .false.
  call track1_bunch_hom (bunch_end, ele, branch%param, bunch_end)
  err = .false.
  ! Only do warning if previous element needed
  if (ele%key == taylor$ .and. csr_param%print_taylor_warning) then
    ele0 => pointer_to_next_ele (ele, -1)
    if (ele0%csr_calc_on) call out_io (s_warn$, r_name, &
                        'CSR calc for taylor element not done: ' // ele%name)
  endif
  return
endif

! n_step is the number of steps to take when tracking through the element.
! csr_top%ds_step is the true step length.

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

! Calculate beam centroid info at element edges

allocate (csr_top%source_ele(0:ele%ix_ele))

do i = 0, ele%ix_ele
  source_ele => csr_top%source_ele(i)
  source_ele%ele => branch%ele(i)  ! Pointer to the P' element
  s_ele => source_ele%ele

  if (present(centroid)) then
    source_ele%orbit_c = centroid(i)
    vec = source_ele%orbit_c%vec
    floor%r = [vec(1), vec(3), s_ele%value(l$)]
    source_ele%floor_c = coords_local_curvilinear_to_floor (floor, s_ele)
    source_ele%floor_c%theta = s_ele%floor%theta + asin(vec(2) / sqrt((1+vec(6)**2 - vec(2)**2)))
  else
    call init_coord (source_ele%orbit_c, ele = s_ele, element_end = downstream_end$)
    source_ele%floor_c = s_ele%floor
  endif

  vec = source_ele%orbit_c%vec
  source_ele%floor_c%theta = s_ele%floor%theta - asin(vec(2) / sqrt((1 + vec(6))**2 - vec(2)**2 - vec(4)**2))

  if (s_ele%value(l$) /= 0) then
    vec0 = csr_top%source_ele(i-1)%orbit_c%vec
    vec = source_ele%orbit_c%vec
    call create_a_spline (source_ele%spline, [0.0_rp, vec0(1)], [s_ele%value(l$), vec(1)], &
                          vec0(2) / sqrt((1+vec0(6))**2 - vec(2)**2), vec(2) / sqrt((1+vec(6))**2 - vec(2)**2))
  endif
enddo

! make sure that ele_len / track_step is an integer.

n_step = max (1, nint(ele%value(l$) / csr_param%ds_track_step))
csr_top%ds_track_step = ele%value(l$) / n_step

auto_bookkeeper = bmad_com%auto_bookkeeper ! save state
bmad_com%auto_bookkeeper = .false.   ! make things go faster

! Loop over the tracking steps
! runt is the element that is tracked through at each step.

do i = 0, n_step

  ! track through the runt

  if (i /= 0) then
    call create_uniform_element_slice (ele, branch%param, i, n_step, runt, s_start, s_end)
    runt%csr_calc_on = .false.
    call track1_bunch_hom (bunch_end, runt, branch%param, bunch_end)
  endif

  s0_step = i * csr_top%ds_track_step
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

  ie = ele%ix_ele
  csr_top%s_kick = s0_step + ele%s - ele%value(l$)
  csr_top%z_kick = (csr_top%source_ele(ie-1)%orbit_c%vec(5) + &
                      dz_from_ele_beginning(s0_step, ele, csr_top%source_ele(ie)%spline)) 

  call csr_bin_particles (bunch_end%particle, csr_top)

  ! ns = 0 is the unshielded kick.
  ! For the shielding image currents never use the small angle approximation

  csr_top%slice(:)%kick_csr = 0
  csr_top%slice(:)%kick_lsc = 0

  do ns = 0, csr_param%n_shield_images
    ! The factor of -1^ns accounts for the sign of the image currents
    ! Take into account that at the endpoints we are only putting in a half kick.
    ! The factor of two is due to there being image currents both above and below.

    csr_top%kick_factor = (-1)**ns
    if (i == 0 .or. i == n_step) csr_top%kick_factor = csr_top%kick_factor / 2
    if (ns /= 0) csr_top%kick_factor = 2* csr_top%kick_factor

    csr_top%y2 = ns * csr_param%beam_chamber_height

    call csr_bin_kicks (ele, s0_step, csr_top)
  enddo

  ! loop over all particles and give them a kick

  do j = 1, size(bunch_end%particle)
    if (bunch_end%particle(j)%state /= alive$) cycle
    call csr_kick_calc (csr_top, bunch_end%particle(j))
  enddo

  call save_bunch_track (bunch_end, ele, s0_step)
enddo

bmad_com%auto_bookkeeper = auto_bookkeeper  ! restore state
err = .false.

end subroutine track1_bunch_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_parcticles (particle, csr_top)
!
! Routine to bin the particles longitudinally in s. 
!
! To avoid noise in the cacluation, every particle is considered to have a 
! triangular distribution with a base length  given by 
!   csr_param%particle_bin_span * csr_top%dz_slice. 
! That is, particles will, in general, overlap multiple bins. 
!
! Input:
!   particle(:)          -- Coord_struct: Array of particles
!   csr_param            -- Csr_parameter_struct: CSR common block (not an argument).
!     %n_bin             -- Number of bins.
!     %particle_bin_span -- Particle length / dz_slice. 
!
! Output:
!   csr_top         -- Csr_top_level_struct: The bin structure.
!     %dz_slice     -- Bin longitudinal length
!     %slice(1:) -- Array of bins.
!-

subroutine csr_bin_particles (particle, csr_top)

implicit none

type this_local_struct   ! Temporary structure 
  real(rp) :: charge     ! how much charge of particle in bin
  real(rp) :: x0, y0     ! particle center
  integer ib             ! bin index
end type

type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p
type (csr_top_level_struct), target :: csr_top
type (this_local_struct), allocatable :: tloc(:)
type (csr_bunch_slice_struct), pointer :: slice

real(rp) z_center, z_min, z_max, f, dz_particle, dz, z_maxval, z_minval
real(rp) zp_center, zp0, zp1, zb0, zb1, charge

integer i, j, n, ix0, ib, ic

character(20) :: r_name = 'csr_bin_particles'

! Init bins...
! The left edge of csr_top%slice(1) is at z_min
! The right edge of csr_top%slice(n_bin) is at z_max
! The first and last bins are empty.

if (.not. csr_param%lcsr_component_on .and. .not. csr_param%lsc_component_on .and. &
    .not. csr_param%tsc_component_on .and. csr_param%n_shield_images == 0) return

z_maxval = maxval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
z_minval = minval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
dz = z_maxval - z_minval
csr_top%dz_slice = dz / (csr_param%n_bin - 2 - (csr_param%particle_bin_span + 1))
csr_top%dz_slice = 1.0000001 * csr_top%dz_slice     ! to prevent round off problems
z_center = (z_maxval + z_minval) / 2
z_min = z_center - csr_param%n_bin * csr_top%dz_slice / 2
z_max = z_center + csr_param%n_bin * csr_top%dz_slice / 2
dz_particle = csr_param%particle_bin_span * csr_top%dz_slice

! allocate memeory for the bins

if (allocated(csr_top%slice)) then
  if (size(csr_top%slice, 1) < csr_param%n_bin) deallocate (csr_top%slice)
endif

if (.not. allocated(csr_top%slice)) &
    allocate (csr_top%slice(csr_param%n_bin), csr_top%kick1(-csr_param%n_bin:csr_param%n_bin))

! Fill in some z information

do i = 1, csr_param%n_bin
  csr_top%slice(i)%z0_edge  = z_min + (i - 1) * csr_top%dz_slice
  csr_top%slice(i)%z_center = csr_top%slice(i)%z0_edge + csr_top%dz_slice / 2
  csr_top%slice(i)%z1_edge  = csr_top%slice(i)%z0_edge + csr_top%dz_slice
enddo

! Init the tloc structure...
! Each particle is distributed longitudinally in a triangular fashion.
! The tloc records how much charge for a given particle is in a bin.

n = size(particle) * (csr_param%particle_bin_span + 2)
allocate (tloc(n))
tloc%ib = -1

! Compute the particle distribution center in each bin

csr_top%slice(:)%charge = 0
csr_top%slice(:)%x0 = 0
csr_top%slice(:)%y0 = 0
csr_top%slice(:)%sig_x = 0
csr_top%slice(:)%sig_y = 0
csr_top%slice(:)%dcharge_density_dz = 0

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
  ix0 = nint((zp0 - z_min) / csr_top%dz_slice)  ! left most bin index
  do j = 0, csr_param%particle_bin_span+1
    ib = j + ix0
    slice => csr_top%slice(ib)
    zb0 = csr_top%slice(ib)%z0_edge
    zb1 = csr_top%slice(ib)%z1_edge   ! edges of the bin
    charge = charge_in_bin (zb0, zb1)
    slice%charge = slice%charge + charge
    slice%x0 = slice%x0 + p%vec(1) * charge
    slice%y0 = slice%y0 + p%vec(3) * charge
    ic = ic + 1
    tloc(ic)%charge = charge
    tloc(ic)%x0 = p%vec(1)
    tloc(ic)%y0 = p%vec(3)
    tloc(ic)%ib = ib
  enddo
enddo

do ib = 1, csr_param%n_bin
  if (ib /= 1) csr_top%slice(ib)%dcharge_density_dz = &
                  (csr_top%slice(ib)%charge - csr_top%slice(ib-1)%charge) / csr_top%dz_slice**2
  if (csr_top%slice(ib)%charge == 0) cycle
  csr_top%slice(ib)%x0 = csr_top%slice(ib)%x0 / csr_top%slice(ib)%charge
  csr_top%slice(ib)%y0 = csr_top%slice(ib)%y0 / csr_top%slice(ib)%charge
enddo

! Compute the particle distribution sigmas in each bin
! Abs is used instead of the usual formula to lessen the effect
! of non-Gaussian tails

do ic = 1, size(tloc)
  if (tloc(ic)%ib < 0) cycle
  slice => csr_top%slice(tloc(ic)%ib)
  slice%sig_x = slice%sig_x + abs(tloc(ic)%x0 - slice%x0) * tloc(ic)%charge
  slice%sig_y = slice%sig_y + abs(tloc(ic)%y0 - slice%y0) * tloc(ic)%charge
enddo

f = sqrt(pi/2)
do ib = 1, csr_param%n_bin
  slice => csr_top%slice(ib)
  if (slice%charge == 0) cycle
  slice%sig_x = f * slice%sig_x / slice%charge
  slice%sig_y = f * slice%sig_y / slice%charge
  slice%lsc_d0 = slice%sig_x * slice%sig_y
  if (slice%sig_x == 0 .and. slice%sig_y == 0) then
    slice%lsc_d1 = 0
  else
    slice%lsc_d1 = (slice%sig_x**2 + slice%sig_y**2) / (slice%sig_x + slice%sig_y)
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
! Subroutine csr_bin_kicks (ele, ds_kick_pt, csr_top)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Input:
!   ele          -- element_struct: Element being tracked through.
!   ds_kick_pt   -- real(rp): Distance between the beginning of the element we are
!                    tracking through and the kick point (which is within this element).
!   csr_top      -- csr_top_level_struct: 
!
! Output:
!   csr_top         -- csr_top_level_struct: 
!     %kick1(:)          -- CSR kick calculation bin array. 
!     %slice(:)%kick_csr -- Integrated kick
!-

subroutine csr_bin_kicks (ele, ds_kick_pt, csr_top)

implicit none

type (csr_top_level_struct), target :: csr_top
type (branch_struct), pointer :: branch
type (ele_struct) ele
type (csr_kick1_struct), pointer :: kick1

real(rp) ds_kick_pt, s_kick, coef, e_tot, f1

integer i, n_bin, ix_source

character(16) :: r_name = 'csr_bin_kicks'

! Assume a linear energy gain in a cavity

branch => ele%branch
f1 = ds_kick_pt / ele%value(l$)
e_tot = f1 * branch%ele(ele%ix_ele-1)%value(e_tot$) + (1 - f1) * ele%value(e_tot$)
call convert_total_energy_to (e_tot, branch%param%particle, csr_top%gamma, beta = csr_top%beta)
csr_top%gamma2 = csr_top%gamma**2
csr_top%rel_mass = mass_of(branch%param%particle) / m_electron 
csr_top%particle = branch%param%particle
s_kick = ele%s + ds_kick_pt - ele%value(l$) ! absolute s value at point P.
csr_top%s_kick = s_kick

! The kick point P is fixed.
! The source point P' varies from bin to bin.
! As we go from bin to bin, P' will move backward
! Initially source point assumed in same element that contains the kick point.

ix_source = ele%ix_ele

! Loop over all kick1 bins and compute the kick or kick integral.
! The loop steps in increasing dz since that is what next_element_params_calc expects.

do i = lbound(csr_top%kick1, 1), ubound(csr_top%kick1, 1)

  kick1 => csr_top%kick1(i)
  kick1%dz_particles = i * csr_top%dz_slice
  kick1%ix_ele_source = ix_source

  ! Calculate what element the kick point is in.

  do
    kick1%s_source = s_source_calc(kick1, csr_top)
    if (kick1%s_source /= real_garbage$) exit       ! If in source element exit loop
    kick1%ix_ele_source = kick1%ix_ele_source - 1
  enddo

  ! calculate csr.
  ! I_csr is only calculated for particles with y = 0 and not for image currents.

  if (csr_top%y2 == 0) then
    call I_csr (kick1, i, csr_top)
    if (csr_top%kick1(i)%I_int_csr == 0 .and. i /= lbound(csr_top%kick1, 1)) then
      csr_top%kick1(i)%I_int_csr = (csr_top%kick1(i)%I_csr + csr_top%kick1(i-1)%I_csr) * csr_top%dz_slice / 2
    endif
  else
    call kick_image_charge (kick1, csr_top)
  endif

enddo

! 

coef = csr_top%ds_track_step * r_e / (csr_top%rel_mass * e_charge * abs(charge_of(branch%param%particle)) * csr_top%gamma)
n_bin = csr_param%n_bin

! CSR & Image charge kick

if (csr_top%y2 == 0) then
  if (csr_param%lcsr_component_on) then
    do i = 1, n_bin
      csr_top%slice(i)%kick_csr = coef * dot_product(csr_top%kick1(i:1:-1)%I_int_csr, csr_top%slice(1:i)%dcharge_density_dz)
    enddo
  endif

else  ! Image charge
  do i = 1, n_bin
    csr_top%slice(i)%kick_csr = csr_top%slice(i)%kick_csr + coef * &
                  dot_product(csr_top%kick1(i-1:i-n_bin:-1)%k_csr, csr_top%slice(1:n_bin)%charge)
  enddo
endif

! Space charge kick

if (csr_param%lsc_component_on) then
  if (csr_top%y2 == 0) then
    call lsc_y0_kick_calc (csr_top)
  endif
endif

end subroutine csr_bin_kicks
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_y0_kick_calc (csr_top)
!
! Routine to cache intermediate values needed for the lsc calculation.
! This routine is not for image currents.
!
! Modules needed:
!   use lsc_mod
!
! Input:
!   csr_top       -- csr_top_level_struct: 
!     %slice(:)   -- bin array of particle averages.
!
! Output:
!   bin     -- csr_top_level_struct: Binned particle averages.
!     %slice(:)%kick_lsc -- Integrated kick.
!-

subroutine lsc_y0_kick_calc (csr_top)

implicit none

type (csr_top_level_struct), target :: csr_top
type (csr_bunch_slice_struct), pointer :: slice

real(rp) sx, sy, a, b, c, dz, factor, sig_x_ave, sig_y_ave, charge_tot

integer i, j

character(*), parameter :: r_name = 'lsc_y0_kick_calc'

! If there are too few particles in a bin the sigma calc may give a bad value.
! This can be a problem if the computed sigma is small.
! Therefore  ignore any bins with a very small sigma. 
! To know what is "small" is, compute the average sigma

charge_tot = sum(csr_top%slice(:)%charge)
if (charge_tot == 0) return

sig_x_ave = dot_product(csr_top%slice(:)%sig_x, csr_top%slice(:)%charge) / charge_tot
sig_y_ave = dot_product(csr_top%slice(:)%sig_y, csr_top%slice(:)%charge) / charge_tot

if (sig_y_ave == 0 .or. sig_x_ave == 0) return  ! Symptom of not enough particles.

! Compute the kick at the center of each bin

csr_top%slice(:)%kick_lsc = 0
if (.not. csr_param%lsc_component_on) return

do i = 1, csr_param%n_bin
  slice => csr_top%slice(i)
  sx = slice%sig_x
  sy = slice%sig_y
  if (sx < sig_x_ave * csr_param%sigma_cutoff .or. sy < sig_y_ave * csr_param%sigma_cutoff) then
    slice%sig_x = 0  ! Mark for tsc calc.
    slice%sig_y = 0
    cycle
  endif
  a = sx * sy
  b = csr_top%gamma * (sx**2 + sy**2) / (sx + sy)
  c = csr_top%gamma**2

  do j = 1, csr_param%n_bin
    if (i == j) cycle
    dz = csr_top%slice(j)%z_center - csr_top%slice(i)%z_center
    csr_top%slice(j)%kick_lsc = csr_top%slice(j)%kick_lsc + &
                   slice%charge * sign(1.0_rp, dz) / (a + b * abs(dz) + c * dz**2)
  enddo

enddo

factor = csr_top%kick_factor * csr_top%ds_track_step * r_e / &
          (csr_top%rel_mass * e_charge * abs(charge_of(csr_top%particle)) * csr_top%gamma)
csr_top%slice(:)%kick_lsc = factor * csr_top%slice(:)%kick_lsc

end subroutine lsc_y0_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine I_csr (kick1, i_bin, csr_top) 
!
! Routine to calculate the CSR kick integral.
!
! Input:
!   kick1      -- csr_kick1_struct: 
!   i_bin      -- integer: Bin index.
!   csr_top    -- csr_top_level_struct:
!
! Output:
!   kick1     -- csr_kick1_struct: 
!     %I_csr     -- real(rp): CSR kick integral.
!     %I_int_csr -- real(rp): Integral of I_csr. Only calculated for i_bin = 1 since it is not needed otherwise.
!-

subroutine I_csr (kick1, i_bin, csr_top)

implicit none

type (csr_kick1_struct) kick1
type (csr_top_level_struct) csr_top

real(rp) z, gam2, g, ds
integer i_bin

! 

kick1%I_int_csr = 0

z = kick1%dz_particles
if (z <= 0) then
  kick1%I_csr = 0
  return
endif

gam2 = csr_top%gamma2

kick1%I_csr = csr_top%kick_factor * (1 / (gam2 * z) - &
    2 * (1 + gam2 * kick1%theta_tot * kick1%theta_sl) /(kick1%L * (1 + gam2 * kick1%theta_sl**2)))
if (i_bin == 1) then
  ds = (csr_top%s_kick + csr_top%z_kick) - (kick1%s_source + kick1%z_source)
  kick1%I_int_csr = csr_top%kick_factor * (-kick1%theta_tot**2 / 4 + log(2*gam2 * z / ds) / gam2)
endif

end subroutine I_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kick_image_charge (kick1, csr_top) 
!
! Routine to calculate the image charge kick.
!
! Input:
!   kick1    -- Csr_kick1_struct: 
!   csr_top      -- Csr_top_level_struct:
!
! Output:
!   kick1    -- Csr_kick1_struct: 
!     %k_csr -- Real(rp): Image charge kick.
!-

subroutine kick_image_charge (kick1, csr_top)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_kick1_struct), pointer :: k
type (csr_top_level_struct) csr_top

real(rp) N_vec(3), G_vec(3), B_vec(3), Bp_vec(3), NBp_vec(3), NBpG_vec(3), rad_cross_vec(3)
real(rp) z, sin_phi, cos_phi, OneNBp, OneNBp3, radiate, coulomb1, theta

!

k => kick1

k%k_csr = 0
z = k%dz_particles

N_vec = (csr_top%floor_k%r - k%floor_s%r) / k%L

theta = csr_top%floor_k%theta
B_vec = [sin(theta), 0.0_rp, cos(theta)]        ! Beta vector at kicked point

theta = k%floor_s%theta
Bp_vec = csr_top%beta * [sin(theta), 0.0_rp, cos(theta) ]  ! beta vector at source point
G_vec = csr_top%beta**2 * k%g_bend * [-cos(theta), 0.0_rp, sin(theta) ] ! Acceleration vector

OneNBp = 1 - sum(N_vec * Bp_vec)
OneNBp3 = OneNBp**3

NBp_vec = N_vec - Bp_vec
NBpG_vec = cross_product(NBp_vec, G_vec)
rad_cross_vec = cross_product(N_vec, NBpG_vec)

radiate  = dot_product (B_vec, rad_cross_vec) / (k%L * OneNBp3)
coulomb1 = dot_product (B_vec, NBp_vec) / (csr_top%gamma2 * k%L**2 * OneNBp3)
kick1%k_csr = csr_top%kick_factor * (radiate + coulomb1)

end subroutine kick_image_charge

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function s_source_calc (kick1, csr_top) result (s_source)
!
! Routine to calculate the distance between source and kick points.
!
! Input:
!   kick1              -- csr_kick1_struct:
!   csr_top            -- csr_top_level_struct:
!
! Output:
!   s_source    -- real(rp): source s-position
!   source_ele -- source_ele_struct: Geometric values
!-

function s_source_calc (kick1, csr_top) result (s_source)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_top_level_struct), target :: csr_top
type (csr_source_ele_struct), pointer :: source_ele
type (ele_struct), pointer :: s_ele
type (floor_position_struct), pointer :: fk, f0, fs

real(rp) a, b, c, dz, s_source, beta2, L0, Lz, s0
real(rp) z0, z1

character(*), parameter :: r_name = 's_sourcecsr'

! If at beginning of lattice assume an infinite drift.
! s_source will be negative

dz = kick1%dz_particles   ! Target distance.
source_ele => csr_top%source_ele(kick1%ix_ele_source)
s_ele => source_ele%ele
beta2 = csr_top%beta**2
f0 => source_ele%floor_c
fk => csr_top%floor_k
fs => kick1%floor_s

if (s_ele%ix_ele == 0) then
  L0 = sqrt((fk%r(1) - f0%r(1))**2 + (fk%r(3) - f0%r(3))**2 + csr_top%y2**2)
  Lz = L0 * cos(f0%theta)
  a = 1 - beta2 * Lz**2
  b = 2 * (csr_top%s_kick + csr_top%z_kick - dz - beta2 * Lz)
  c = ((csr_top%s_kick + csr_top%z_kick) - (source_ele%orbit_c%s + source_ele%orbit_c%vec(5)) - dz)**2 - beta2 * L0
  s_source = (-b + sqrt(b**2 - 4 * a * c)) / (2 * a)

  fs%r = [f0%r(1) + s_source * sin(f0%theta), csr_top%y2, f0%r(3) + s_source * cos(f0%theta)]
  kick1%L = sqrt(dot_product(fk%r-fs%r, fk%r-fs%r))
  return
endif

! Look at ends of the source element to make sure that we are within the element

s_source = real_garbage$

if (s_ele%value(l$) == 0) return  ! Source point cannot be here.
if (0 < ddz_calc_csr(s_ele%value(l$)) .or. ddz_calc_csr(0.0_rp) < 0) return

s_source = zbrent (ddz_calc_csr, 0.0_rp, s_ele%value(l$), 1d-8)

!----------------------------------------------------------------------------
contains

!+
! Function ddz_calc_csr (s) result (ddz_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle at constant time minus the target distance.
!
! Input:
!   s        -- Real(rp): Distance from start of element.
!
! Output:
!   ddz_this -- Real(rp): Distance between source and kick particles: Calculated - Wanted.
!-

function ddz_calc_csr (s) result (ddz_this)

implicit none

type (floor_position_struct) fs
real(rp), intent(in) :: s
real(rp) :: ddz_this, x, z, l_vec(3)

character(*), parameter :: r_name = 'ddz_calc_csr'

! 

call spline1_evaluate(source_ele%spline, s, x)
fs%r = [x, csr_top%y2, s]
fs = coords_local_curvilinear_to_floor (fs, s_ele)



l_vec(1) = fk%r(1) - fs%r(1)
l_vec(2) = csr_top%y2
l_vec(3) = fk%r(3) - fs%r(3)
kick1%L = sqrt(dot_product(L_vec, L_vec))

z = (csr_top%s_kick + csr_top%z_kick) - (s + s_ele%s - s_ele%value(l$) + &
                    dz_from_ele_beginning(s, s_ele, source_ele%spline)) - csr_top%beta * kick1%L
ddz_this = z - kick1%dz_particles

end function ddz_calc_csr

end function s_source_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_kick_calc (csr_top, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
!   csr_top  -- csr_top_level_struct: 
!   particle -- coord_struct: Particle to kick.
!     %vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Coord_struct: Particle to kick.
!     %vec(6) -- Final particle energy.
!-

subroutine csr_kick_calc (csr_top, particle)

implicit none

type (csr_top_level_struct), target :: csr_top
type (coord_struct), target :: particle
type (csr_bunch_slice_struct), pointer :: slice

real(rp) zp, r1, r0, dz, dpz, kx, ky, f0, f, beta0
real(rp), pointer :: vec(:)
integer i, i0, i_del

! We use a weighted average between %kick1(j)%I_csr and %kick1(j+1)%I_csr
! so that the integral varies smoothly as a function of particle%vec(5).

if (.not. allocated(csr_top%slice)) return  ! True if kicks are turned off.

zp = particle%vec(5)
i0 = int((zp - csr_top%slice(1)%z_center) / csr_top%dz_slice) + 1
r1 = (zp - csr_top%slice(i0)%z_center) / csr_top%dz_slice
r0 = 1 - r1
vec => particle%vec

if (r1 < 0 .or. r1 > 1 .or. i0 < 1 .or. i0 >= csr_param%n_bin) then
  print *, 'CSR INTERNAL ERROR!'
  if (global_com%exit_on_error) call err_exit
endif

vec(6) = vec(6) + r0 * csr_top%slice(i0)%kick_csr + r1 * csr_top%slice(i0+1)%kick_csr

! Longitudinal space charge

if (csr_param%lsc_component_on) then
  vec(6) = vec(6) + r0 * csr_top%slice(i0)%kick_lsc + r1 * csr_top%slice(i0+1)%kick_lsc
endif

! Must update beta and z due to the energy change

beta0 = particle%beta
call convert_pc_to ((1+vec(6))* particle%p0c, particle%species, beta = particle%beta)
vec(5) = vec(5) * particle%beta / beta0

! Transverse space charge.

if (csr_param%tsc_component_on) then
  f0 = csr_top%kick_factor * csr_top%ds_track_step * r_e / (twopi * &
           csr_top%dz_slice * csr_top%rel_mass * e_charge * abs(charge_of(csr_top%particle)) * csr_top%gamma**3)

  slice => csr_top%slice(i0)
  if (slice%sig_x /= 0) then
    call bbi_kick ((vec(1)-slice%x0)/slice%sig_x, (vec(3)-slice%y0)/slice%sig_y, &
                                                       slice%sig_y/slice%sig_x, kx, ky)
    f = f0 * r0 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    vec(2) = vec(2) - kx * f
    vec(4) = vec(4) - ky * f
  endif

  slice => csr_top%slice(i0+1)
  if (slice%sig_x /= 0) then
    call bbi_kick ((vec(1)-slice%x0)/slice%sig_x, (vec(3)-slice%y0)/slice%sig_y, &
                                                       slice%sig_y/slice%sig_x, kx, ky)
    f = f0 * r1 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    vec(2) = vec(2) - kx * f   
    vec(4) = vec(4) - ky * f
  endif

endif

end subroutine csr_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function dz_from_ele_beginning(s, ele, spline) result (this_dz)
!
! Routine to calculate the change in particle z from the beginning of an element to a particular point.
!
! Input:
!   s       -- real(rp): s-position from beginning of element.
!   ele     -- ele_struct: Element tracked through
!   spline  -- spline_struct: Spline of x-position as a function of s.
!
! Output:
!   this_dz -- real(rp): change in particle z.
!-

function dz_from_ele_beginning(s, ele, spline) result (this_dz)

implicit none

type (ele_struct) ele
type (spline_struct) spline

real(rp) s, this_dz, c(0:3)

! x' = c(1) + 2*c2*s + 3*c3*s^2
! dz = -Integral: x'^2/2 ds

c = spline%coef
this_dz = -(c(1)**2 * s + 2*c(1)*c(2) * s**2 + (6*c(1)*c(3) + 4*c(2)**2) * s**3/3 + 3*c(2)*c(3) * s**4 + 9*c(3)**2 * s**5/5) / 2

! For bends add: dz = -Integral: x/rho ds

if (ele%key /= sbend$) return
this_dz = this_dz - ele%value(g$) * (c(0)*s + c(1)*s**2/2 + c(2)*s**3/3 + c(3)* s**4/4)

end function dz_from_ele_beginning

end module
