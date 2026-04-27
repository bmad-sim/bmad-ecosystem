!+
! See the paper:
!   "Coherent Synchrotron Radiation Simulations for Off-Axis Beams Using the Bmad Toolkit"
!   D. Sagan & C. Mayes
!   Proceedings of IPAC2017, Copenhagen, Denmark THPAB076
!   https://accelconf.web.cern.ch/ipac2017/papers/thpab076.pdf
!-

module csr_and_space_charge_mod

use beam_utils
use bmad_interface
use spline_mod
use super_recipes_mod, only: super_zbrent
use open_spacecharge_mod
use csr3d_mod, only: csr3d_steady_state_solver

! csr_ele_info_struct holds info for a particular lattice element
! The centroid "chord" is the line from the centroid position at the element entrance to
! the centroid position at the element exit end.

type csr_ele_info_struct
  type (ele_struct), pointer :: ele            ! lattice element
  type (coord_struct) orbit0, orbit1           ! centroid orbit at entrance/exit ends
  type (floor_position_struct) floor0, floor1  ! Floor position of centroid at entrance/exit ends
  type (floor_position_struct) ref_floor0, ref_floor1  ! Floor position of element ref coords at entrance/exit ends
  type (spline_struct) spline                  ! Spline for centroid orbit. spline%x = distance along chord.
                                               !   The spline is zero at the ends by construction.
  real(rp) theta_chord                         ! Reference angle of chord in z-x plane
  real(rp) L_chord                             ! Chord Length. Negative if bunch moves backwards in element.
  real(rp) dL_s                                ! L_s(of element) - L_chord
end type

type csr_bunch_slice_struct  ! Structure for a single particle bin.
  real(rp) :: x0 = 0, y0 = 0 ! Transverse center of the particle distrubution
  real(rp) :: z0_edge = 0    ! Left (min z) edge of bin
  real(rp) :: z1_edge = 0    ! Right (max z) edge of bin
  real(rp) :: z_center = 0   ! z at center of bin.
  real(rp) :: sig_x = 0      ! particle's RMS width
  real(rp) :: sig_y = 0      ! particle's RMS width
  real(rp) :: charge = 0     ! charge of the particles
  real(rp) :: dcharge_density_dz = 0      ! Charge density gradient 
  real(rp) :: edge_dcharge_density_dz = 0 ! gradient between this and preceeding bin. [Evaluated at bin edge.]
  real(rp) :: kick_csr = 0                ! CSR kick
  real(rp) :: coef_lsc_plus(0:2,0:2) = 0  ! LSC Kick coefs.
  real(rp) :: coef_lsc_minus(0:2,0:2) = 0 ! LSC Kick coefs.
  real(rp) :: kick_lsc = 0
  real(rp) :: n_particle = 0 ! Number of particles in slice can be a fraction since particles span multiple bins.
end type

! csr_kick1_struct stores the CSR kick, kick integral etc. for a given source and kick positions.
! This structure also holds info on parameters that go into the kick calculation.
! Since an integration step involves one kick position and many source positions,
! the information that is only kick position dependent is held in the csr_struct and
! the csr_struct holds an array of csr_kick1_structs, one for each dz.

type csr_kick1_struct ! Sub-structure for csr calculation cache
  real(rp) I_csr            ! Kick integral.
  real(rp) I_int_csr        ! Integrated Kick integral.
  real(rp) image_kick_csr   ! kick.
  real(rp) L_vec(3)         ! L vector in global coordinates.
  real(rp) L                ! Distance between source and kick points.
  real(rp) dL               ! = epsilon_L = Ls - L
  real(rp) dz_particles     ! Kicked particle - source particle position at constant time.
  real(rp) s_chord_source   ! Source point coordinate along chord.
  real(rp) theta_L          ! Angle of L vector
  real(rp) theta_sl         ! Angle between velocity of particle at source pt and L
  real(rp) theta_lk         ! Angle between L and velocity of kicked particle
  integer ix_ele_source     ! Source element index.
  type (floor_position_struct) floor_s  ! Floor position of source pt
end type

type csr_particle_position_struct
  real(rp) :: r(3)       ! particle position
  real(rp) :: charge     ! particle charge
end type

type csr_struct           ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) rel_mass             ! m_particle / m_electron
  real(rp) beta                 ! Relativistic beta factor.
  real(rp) :: dz_slice = 0      ! Bin width
  real(rp) ds_track_step        ! True step size
  real(rp) s_kick               ! Kick point longitudinal location (element ref coords) from entrance end
  real(rp) s_chord_kick         ! Kick point along beam centroid line
  real(rp) y_source             ! Height of source particle.
  real(rp) kick_factor          ! Coefficient to scale the kick
  real(rp) actual_track_step    ! ds_track_step scalled by Length_centroid_chord / Length_element ratio
  real(rp) x0_bunch, y0_bunch   ! Bunch centroid
  type(floor_position_struct) floor_k   ! Floor coords at kick point
  integer species                       ! Particle type
  integer ix_ele_kick                   ! Same as element being tracked through.
  type (csr_bunch_slice_struct), allocatable :: slice(:)    ! slice(i) refers to the i^th bunch slice.
  type (csr_kick1_struct), allocatable :: kick1(:)          ! kick1(i) referes to the kick between two slices i bins apart.
  type (csr_ele_info_struct), allocatable :: eleinfo(:)     ! Element-by-element information.
  type (ele_struct), pointer :: kick_ele                    ! Element where the kick pt is == ele tracked through.
  type (mesh3d_struct) :: mesh3d
  type (csr_particle_position_struct), allocatable :: position(:)
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch, ele, centroid, err, s_start, s_end, bunch_track)
!
! Routine to track a bunch of particles through an element with csr radiation effects.
!
! Input:
!   bunch         -- Bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: The element to track through. Must be part of a lattice.
!   centroid(0:)  -- coord_struct, Approximate beam centroid orbit for the lattice branch.
!                      Calculate this before beam tracking by tracking a single particle.
!   s_start       -- real(rp), optional: Starting position relative to ele. Default = 0
!   s_end         -- real(rp), optional: Ending position. Default is ele length.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- Bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. EG: Too many particles lost.
!   bunch_track -- bunch_track_struct, optional: track information if the tracking method does
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!-

subroutine track1_bunch_csr (bunch, ele, centroid, err, s_start, s_end, bunch_track)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: pt, c0
type (ele_struct), target :: ele
type (bunch_track_struct), optional :: bunch_track
type (branch_struct), pointer :: branch
type (ele_struct) :: runt
type (ele_struct), pointer :: ele0, s_ele
type (csr_struct), target :: csr
type (csr_ele_info_struct), pointer :: eleinfo
type (coord_struct), target :: centroid(0:)
type (floor_position_struct) floor

real(rp), optional :: s_start, s_end
real(rp) s0_step, vec0(6), vec(6), theta_chord, theta0, theta1, L
real(rp) e_tot, f1, x, z, ds_step, s_save_last

integer i, j, n, ie, ns, nb, n_step, n_live, i_step
integer :: iu_wake

character(*), parameter :: r_name = 'track1_bunch_csr'
logical err, err_flag, parallel0, parallel1

! Init

err = .true.
csr%kick_ele => ele    ! Element where the kick pt is == ele tracked through.
branch => pointer_to_branch(ele)

if (ele%space_charge_method == fft_3d$) then
  c0 => centroid(ele%ix_ele)
  call convert_pc_to((1+c0%vec(6)) * c0%p0c, c0%species, gamma = csr%mesh3d%gamma)
  csr%mesh3d%nhi = space_charge_com%space_charge_mesh_size
endif

! No CSR for a zero length element.
! And taylor elements get ignored.

if (ele%value(l$) == 0 .or. ele%key == taylor$) then
  call track1_bunch_hom (bunch, ele)
  err = .false.
  return
endif

! n_step is the number of steps to take when tracking through the element.
! csr%ds_step is the true step length.

if (ele%csr_method == one_dim$ .and. (ele%key == wiggler$ .or. ele%key == undulator$) .and. &
                                (ele%field_calc == planar_model$ .or. ele%field_calc == helical_model$)) then
  call out_io (s_warn$, r_name, 'CALCULATION OF CSR EFFECTS IN PLANAR OR HELICAL MODEL WIGGLERS MAY BE INVALID SINCE', &
                                'CSR TRACKING USES A SPLINE FIT FOR THE CENTROID ORBIT WHICH IS INACCURATE FOR A WIGGLING ORBIT.', &
                                'FOR: ' // ele%name)
endif

if (space_charge_com%n_bin <= space_charge_com%particle_bin_span + 1) then
  if (space_charge_com%n_bin == 0) then
    call out_io (s_error$, r_name, &
            'SPACE_CHARGE_COM%N_BIN IS ZERO WHICH INDICATES THAT THE SPACE_CHARGE_COM STRUCTURE HAS NOT BEEN SET BY THE USER.', &
            'ALL PARTICLES IN THE BUNCH WILL BE MARKED AS LOST.')
  else
    call out_io (s_error$, r_name, &
            'SPACE_CHARGE_COM%N_BIN (= \i0\) MUST BE GREATER THAN SPACE_CHARGE_COM%PARTICLE_BIN_SPAN+1 (= \i0\+1)!', &
            'ALL PARTICLES IN THE BUNCH WILL BE MARKED AS LOST.', &
             i_array = [space_charge_com%n_bin, space_charge_com%particle_bin_span])
  endif
  bunch%particle%state = lost$
  return
endif

! make sure that ele_len / track_step is an integer.

csr%ds_track_step = ele%value(csr_ds_step$)
if (csr%ds_track_step == 0) csr%ds_track_step = space_charge_com%ds_track_step
if (csr%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'NEITHER SPACE_CHARGE_COM%DS_TRACK_STEP NOR CSR_TRACK_STEP FOR THIS ELEMENT ARE SET! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

n_step = max (1, nint(ele%value(l$) / csr%ds_track_step))
csr%ds_track_step = ele%value(l$) / n_step

! Calculate beam centroid info at element edges, etc.

n = min(ele%ix_ele+10, branch%n_ele_track)
allocate (csr%eleinfo(0:n))

do i = 0, n
  eleinfo => csr%eleinfo(i)
  eleinfo%ele => branch%ele(i)  ! Pointer to the P' element
  s_ele => eleinfo%ele
  eleinfo%ref_floor1 = branch%ele(i)%floor
  eleinfo%ref_floor1%r(2) = 0  ! Make sure in horizontal plane

  eleinfo%orbit1 = centroid(i)
  vec = eleinfo%orbit1%vec
  floor%r = [vec(1), vec(3), s_ele%value(l$)]
  eleinfo%floor1 = coords_local_curvilinear_to_floor (floor, s_ele)
  eleinfo%floor1%r(2) = 0  ! Make sure in horizontal plane
  ! The 1d-14 is inserted to avoid 0/0 error at the beginning of an e_gun.
  eleinfo%floor1%theta = s_ele%floor%theta + asin(vec(2) / (1.0_rp + 1d-14 + vec(6)))

  if (i == 0) then
    eleinfo%ref_floor0 = csr%eleinfo(i)%ref_floor1
    eleinfo%floor0   = csr%eleinfo(i)%floor1
    eleinfo%orbit0   = csr%eleinfo(i)%orbit1
  else
    eleinfo%ref_floor0 = csr%eleinfo(i-1)%ref_floor1
    eleinfo%floor0   = csr%eleinfo(i-1)%floor1
    eleinfo%orbit0   = csr%eleinfo(i-1)%orbit1
  endif

  vec0 = eleinfo%orbit0%vec
  vec = eleinfo%orbit1%vec
  theta_chord = atan2(eleinfo%floor1%r(1)-eleinfo%floor0%r(1), eleinfo%floor1%r(3)-eleinfo%floor0%r(3))
  eleinfo%theta_chord = theta_chord

  ! An element with a negative step length (EG off-center beam in a patch which just has an x-pitch), can
  ! have floor%theta and theta_chord 180 degress off. Hence, restrict theta0 & theta1 to be within [-pi/2,pi/2]
  theta0 = modulo2(eleinfo%floor0%theta - theta_chord, pi/2)
  theta1 = modulo2(eleinfo%floor1%theta - theta_chord, pi/2)
  eleinfo%L_chord = sqrt((eleinfo%floor1%r(1)-eleinfo%floor0%r(1))**2 + (eleinfo%floor1%r(3)-eleinfo%floor0%r(3))**2)
  if (abs(eleinfo%L_chord) < 1d-8) then   ! 1d-8 is rather arbitrary.
    eleinfo%spline = spline_struct()
    eleinfo%dL_s = 0
    cycle
  endif

  if (s_ele%key == match$) cycle  ! Match elements offsets are ignored so essentially they are like markers.

  ! With a negative step length, %L_chord is negative
  parallel0 = (abs(modulo2(eleinfo%floor0%theta - theta_chord, pi)) < pi/2)
  parallel1 = (abs(modulo2(eleinfo%floor1%theta - theta_chord, pi)) < pi/2)
  if (parallel0 .neqv. parallel1) then
    call out_io (s_fatal$, r_name, 'VERY CONFUSED CSR CALCULATION! PLEASE SEEK HELP ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
  if (.not. parallel0) eleinfo%L_chord = -eleinfo%L_chord

  eleinfo%spline = create_a_spline ([0.0_rp, 0.0_rp], [eleinfo%L_chord, 0.0_rp], theta0, theta1)
  eleinfo%dL_s = dspline_len(0.0_rp, eleinfo%L_chord, eleinfo%spline)
enddo

!

csr%species = bunch%particle(1)%species
csr%ix_ele_kick = ele%ix_ele
csr%actual_track_step = csr%ds_track_step * (csr%eleinfo(ele%ix_ele)%L_chord / ele%value(l$))

call save_a_bunch_step (ele, bunch, bunch_track, s_start)

!----------------------------------------------------------------------------------------
! Loop over the tracking steps
! runt is the element that is tracked through at each step.

do i_step = 0, n_step

  ! track through the runt

  if (i_step /= 0) then
    call element_slice_iterator (ele, branch%param, i_step, n_step, runt, s_start, s_end)
    call track1_bunch_hom (bunch, runt)
  endif

  s0_step = i_step * csr%ds_track_step
  if (present(s_start)) s0_step = s0_step + s_start

  ! Cannot do a realistic calculation if there are less particles than bins

  n_live = count(bunch%particle%state == alive$)
  if (n_live < space_charge_com%n_bin) then
    call out_io (s_error$, r_name, 'NUMBER OF LIVE PARTICLES: \i0\ ', &
                          'LESS THAN NUMBER OF BINS FOR CSR CALC.', &
                          'AT ELEMENT: ' // trim(ele%name) // '  [# \i0\] ', &
                          i_array = [n_live, ele%ix_ele ])
    return
  endif

  ! Assume a linear energy gain in a cavity

  f1 = s0_step / ele%value(l$)
  e_tot = f1 * branch%ele(ele%ix_ele-1)%value(e_tot$) + (1 - f1) * ele%value(e_tot$)
  call convert_total_energy_to (e_tot, branch%param%particle, csr%gamma, beta = csr%beta)
  csr%gamma2 = csr%gamma**2
  csr%rel_mass = mass_of(branch%param%particle) / m_electron 

  call csr_bin_particles (ele, bunch%particle, csr, err_flag); if (err_flag) return

  csr%s_kick = s0_step
  csr%s_chord_kick = s_ref_to_s_chord (s0_step, csr%eleinfo(ele%ix_ele))
  z = csr%s_chord_kick
  x = spline1(csr%eleinfo(ele%ix_ele)%spline, z)
  theta_chord = csr%eleinfo(ele%ix_ele)%theta_chord
  csr%floor_k%r = [x*cos(theta_chord)+z*sin(theta_chord), 0.0_rp, -x*sin(theta_chord)+z*cos(theta_chord)] + &
                      csr%eleinfo(ele%ix_ele)%floor0%r
  csr%floor_k%theta = theta_chord + spline1(csr%eleinfo(ele%ix_ele)%spline, z, 1)

  ! ns = 0 is the unshielded kick.

  if (ele%space_charge_method == slice$ .or. ele%csr_method == one_dim$) then
    do ns = 0, space_charge_com%n_shield_images
      ! The factor of -1^ns accounts for the sign of the image currents
      ! Take into account that at the endpoints we are only putting in a half kick.
      ! The factor of two is due to there being image currents both above and below.

      csr%kick_factor = (-1)**ns
      if (i_step == 0 .or. i_step == n_step) csr%kick_factor = csr%kick_factor / 2
      if (ns /= 0) csr%kick_factor = 2 * csr%kick_factor

      csr%y_source = ns * space_charge_com%beam_chamber_height

      call csr_bin_kicks (ele, s0_step, csr, err_flag)
      if (err_flag) return
    enddo
  endif

  ! Give particles a kick

  call csr_and_sc_apply_kicks (ele, csr, bunch%particle)

  call save_a_bunch_step (ele, bunch, bunch_track, s0_step)

  ! Record wake to file?

  if (space_charge_com%diagnostic_output_file /= '') then
    iu_wake = lunget()
    open (iu_wake, file = space_charge_com%diagnostic_output_file, access = 'append')
    if (i_step == 0) then
      write (iu_wake, '(a)') '!------------------------------------------------------------'
      write (iu_wake, '(a, i6, 2x, a)') '! ', ele%ix_ele, trim(ele%name)
    endif
    write (iu_wake, '(a)') '!#-----------------------------'
    write (iu_wake, '(a, i4, f12.6)') '! Step index:', i_step
    write (iu_wake, '(a, f12.6)') '! S-position:', s0_step
    write (iu_wake, '(a)') '!         Z   Charge/Meter    CSR_Kick/m       I_CSR/m      S_Source' 
    if (allocated(csr%kick1)) then
      ds_step = csr%kick_factor * csr%actual_track_step
      do j = 1, space_charge_com%n_bin
        ele0 => branch%ele(csr%kick1(j)%ix_ele_source)
        write (iu_wake, '(f14.10, 4es14.6)') csr%slice(j)%z_center, &
                    csr%slice(j)%charge/csr%dz_slice, csr%slice(j)%kick_csr/ds_step, &
                    csr%kick1(j)%I_csr/ds_step, ele0%s_start + csr%kick1(j)%s_chord_source
      enddo
    elseif (i_step == 0) then
      write (iu_wake, '(a)') 'Note: CSR wake not calculated for element: ' // ele%name
      write (iu_wake, '(a)') '      Check settings of ele%space_charge_method and ele%csr_method.'
    endif
    close (iu_wake)
  endif

enddo

err = .false.

end subroutine track1_bunch_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_particles (ele, particle, csr)
!
! Routine to bin the particles longitudinally in s. 
!
! To avoid noise in the calculation, every particle is considered to have a 
! triangular distribution with a base length  given by
!   space_charge_com%particle_bin_span * csr%dz_slice.
! That is, particles will, in general, overlap multiple bins. 
!
! Input:
!   ele                  -- ele_struct: Element being tracked through.
!   particle(:)          -- coord_struct: Array of particles
!
! Other:
!   space_charge_com     -- space_charge_common_struct: Space charge/CSR common block
!     %n_bin                 -- Number of bins.
!     %particle_bin_span     -- Particle length / dz_slice. 
!
! Output:
!   csr                  -- csr_struct: The bin structure.
!     %dz_slice             -- Bin longitudinal length
!     %slice(1:)            -- Array of bins.
!-

subroutine csr_bin_particles (ele, particle, csr, err_flag)

implicit none

type (ele_struct) ele
type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p
type (csr_struct), target :: csr
type (csr_bunch_slice_struct), pointer :: slice

real(rp) z_center, z_min, z_max, dz_particle, dz, z_maxval, z_minval, c_tot, n_tot
real(rp) zp_center, zp0, zp1, zb0, zb1, charge, overlap_fraction, charge_tot, sig_x_ave, sig_y_ave, f
integer i, j, n, ix0, ib, ib2, ib_center, n_bin_eff, n_bin, pbs

! Thread-private accumulation arrays for OpenMP binning
real(rp), allocatable :: bin_n_particle(:), bin_charge(:), bin_x0(:), bin_y0(:)
real(rp), allocatable :: bin_sig_x(:), bin_sig_y(:)
real(rp) :: z1_over, z2_over, overlap, p_charge, p_x, p_y, inv_dz2

logical err_flag

character(*), parameter :: r_name = 'csr_bin_particles'

! Init bins...
! The left edge of csr%slice(1) is at z_min
! The right edge of csr%slice(n_bin) is at z_max
! The first and last bins are empty.

err_flag = .false.

if (ele%space_charge_method /= slice$ .and. ele%csr_method /= one_dim$) return

n_bin_eff = space_charge_com%n_bin - 2 - (space_charge_com%particle_bin_span + 1)
if (n_bin_eff < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF CSR BINS TOO SMALL: \i0\ ', &
              'MUST BE GREATER THAN 3 + PARTICLE_BIN_SPAN.', i_array = [space_charge_com%n_bin])
  if (global_com%exit_on_error) call err_exit
  particle%state = lost$
  err_flag = .true.
  return
endif

z_maxval = maxval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
z_minval = minval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
dz = z_maxval - z_minval
csr%dz_slice = dz / n_bin_eff
csr%dz_slice = 1.0000001 * csr%dz_slice     ! to prevent round off problems
z_center = (z_maxval + z_minval) / 2
z_min = z_center - space_charge_com%n_bin * csr%dz_slice / 2
z_max = z_center + space_charge_com%n_bin * csr%dz_slice / 2
dz_particle = space_charge_com%particle_bin_span * csr%dz_slice

if (dz == 0) then
  call out_io (s_fatal$, r_name, 'LONGITUDINAL WIDTH OF BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
endif

! allocate memeory for the bins

if (allocated(csr%slice)) then
  if (size(csr%slice, 1) < space_charge_com%n_bin) deallocate (csr%slice)
endif

if (.not. allocated(csr%slice)) &
    allocate (csr%slice(space_charge_com%n_bin), csr%kick1(-space_charge_com%n_bin:space_charge_com%n_bin))

csr%slice(:) = csr_bunch_slice_struct()  ! Zero everything

! Fill in some z information

do i = 1, space_charge_com%n_bin
  csr%slice(i)%z0_edge  = z_min + (i - 1) * csr%dz_slice
  csr%slice(i)%z_center = csr%slice(i)%z0_edge + csr%dz_slice / 2
  csr%slice(i)%z1_edge  = csr%slice(i)%z0_edge + csr%dz_slice
enddo

! Compute the particle distribution center in each bin

! The contribution to the charge in a bin from a particle is computed from the overlap
! between the particle and the bin.

n_bin = space_charge_com%n_bin
pbs = space_charge_com%particle_bin_span
inv_dz2 = 1.0_rp / dz_particle**2

!$OMP parallel private(bin_n_particle, bin_charge, bin_x0, bin_y0, zp_center, zp0, zp1, ix0, j, ib, zb0, zb1, z1_over, z2_over, overlap, charge, p_charge, p_x, p_y)
allocate(bin_n_particle(n_bin), bin_charge(n_bin), bin_x0(n_bin), bin_y0(n_bin))
bin_n_particle = 0; bin_charge = 0; bin_x0 = 0; bin_y0 = 0

!$OMP do
do i = 1, size(particle)
  if (particle(i)%state /= alive$) cycle
  zp_center = particle(i)%vec(5)
  zp0 = zp_center - dz_particle / 2
  zp1 = zp_center + dz_particle / 2
  p_charge = particle(i)%charge
  p_x = particle(i)%vec(1)
  p_y = particle(i)%vec(3)
  ix0 = nint((zp0 - z_min) / csr%dz_slice)
  do j = 0, pbs+1
    ib = j + ix0
    zb0 = csr%slice(ib)%z0_edge
    zb1 = csr%slice(ib)%z1_edge
    ! Inline particle_overlap_in_bin for thread safety (avoids host-scope variable sharing)
    overlap = 0
    z1_over = max(zp0, zb0)
    z2_over = min(zp_center, zb1)
    if (z2_over > z1_over) overlap = 2 * real(((z2_over - zp0)**2 - (z1_over - zp0)**2), rp) * inv_dz2
    z1_over = max(zp_center, zb0)
    z2_over = min(zp1, zb1)
    if (z2_over > z1_over) overlap = overlap + 2 * real(((z1_over - zp1)**2 - (z2_over - zp1)**2), rp) * inv_dz2
    bin_n_particle(ib) = bin_n_particle(ib) + overlap
    charge = overlap * p_charge
    bin_charge(ib) = bin_charge(ib) + charge
    bin_x0(ib) = bin_x0(ib) + p_x * charge
    bin_y0(ib) = bin_y0(ib) + p_y * charge
  enddo
enddo
!$OMP end do

!$OMP critical
do ib = 1, n_bin
  csr%slice(ib)%n_particle = csr%slice(ib)%n_particle + bin_n_particle(ib)
  csr%slice(ib)%charge = csr%slice(ib)%charge + bin_charge(ib)
  csr%slice(ib)%x0 = csr%slice(ib)%x0 + bin_x0(ib)
  csr%slice(ib)%y0 = csr%slice(ib)%y0 + bin_y0(ib)
enddo
!$OMP end critical
deallocate(bin_n_particle, bin_charge, bin_x0, bin_y0)
!$OMP end parallel

do ib = 1, space_charge_com%n_bin
  if (ib /= 1) csr%slice(ib)%edge_dcharge_density_dz = (csr%slice(ib)%charge - csr%slice(ib-1)%charge) / csr%dz_slice**2
  if (ib == 1) then
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib+1)%charge - csr%slice(ib)%charge) / csr%dz_slice**2
  elseif (ib == space_charge_com%n_bin) then
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib)%charge - csr%slice(ib-1)%charge) / csr%dz_slice**2
  else
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib+1)%charge - csr%slice(ib-1)%charge) / (2 * csr%dz_slice**2)
  endif

  if (csr%slice(ib)%charge == 0) cycle
  csr%slice(ib)%x0 = csr%slice(ib)%x0 / csr%slice(ib)%charge
  csr%slice(ib)%y0 = csr%slice(ib)%y0 / csr%slice(ib)%charge
enddo

csr%x0_bunch = sum(csr%slice%x0 * csr%slice%charge) / sum(csr%slice%charge)
csr%y0_bunch = sum(csr%slice%y0 * csr%slice%charge) / sum(csr%slice%charge)

! Compute the particle distribution sigmas in each bin.
! Abs(x-x0) is used instead of the usual formula involving (x-x0)^2 to lessen the effect
! of non-Gaussian tails.

!$OMP parallel private(bin_sig_x, bin_sig_y, zp_center, zp0, zp1, ix0, j, ib, zb0, zb1, z1_over, z2_over, overlap, charge, p_charge, p_x, p_y)
allocate(bin_sig_x(n_bin), bin_sig_y(n_bin))
bin_sig_x = 0; bin_sig_y = 0

!$OMP do
do i = 1, size(particle)
  if (particle(i)%state /= alive$) cycle
  zp_center = particle(i)%vec(5)
  zp0 = zp_center - dz_particle / 2
  zp1 = zp_center + dz_particle / 2
  p_charge = particle(i)%charge
  p_x = particle(i)%vec(1)
  p_y = particle(i)%vec(3)
  ix0 = nint((zp0 - z_min) / csr%dz_slice)
  do j = 0, pbs+1
    ib = j + ix0
    zb0 = csr%slice(ib)%z0_edge
    zb1 = csr%slice(ib)%z1_edge
    ! Inline particle_overlap_in_bin
    overlap = 0
    z1_over = max(zp0, zb0)
    z2_over = min(zp_center, zb1)
    if (z2_over > z1_over) overlap = 2 * real(((z2_over - zp0)**2 - (z1_over - zp0)**2), rp) * inv_dz2
    z1_over = max(zp_center, zb0)
    z2_over = min(zp1, zb1)
    if (z2_over > z1_over) overlap = overlap + 2 * real(((z1_over - zp1)**2 - (z2_over - zp1)**2), rp) * inv_dz2
    charge = overlap * p_charge
    bin_sig_x(ib) = bin_sig_x(ib) + abs(p_x - csr%slice(ib)%x0) * charge
    bin_sig_y(ib) = bin_sig_y(ib) + abs(p_y - csr%slice(ib)%y0) * charge
  enddo
enddo
!$OMP end do

!$OMP critical
do ib = 1, n_bin
  csr%slice(ib)%sig_x = csr%slice(ib)%sig_x + bin_sig_x(ib)
  csr%slice(ib)%sig_y = csr%slice(ib)%sig_y + bin_sig_y(ib)
enddo
!$OMP end critical
deallocate(bin_sig_x, bin_sig_y)
!$OMP end parallel

charge_tot = 0;  sig_x_ave = 0;  sig_y_ave = 0
f = sqrt(pi/2)  ! This corrects for the fact that |x - x0| is used instead of (x - x0)^2 to compute the sigma.
do ib = 1, space_charge_com%n_bin
  slice => csr%slice(ib)
  if (slice%n_particle < space_charge_com%sc_min_in_bin) cycle
  slice%sig_x = f * slice%sig_x / slice%charge
  slice%sig_y = f * slice%sig_y / slice%charge
  charge_tot = charge_tot + slice%charge
  sig_x_ave = slice%sig_x * slice%charge
  sig_y_ave = slice%sig_y * slice%charge
enddo

if (charge_tot == 0) then   ! Not enough particles for calc
  csr%slice%sig_x = 1  ! Something large
  csr%slice%sig_y = 1  ! Something large
  return
endif

sig_x_ave = sig_x_ave / charge_tot
sig_y_ave = sig_y_ave / charge_tot

! At the ends, for bins where there are not enough particles to calculate sigma, use
! the sigmas of nearest bin that has a valid sigmas

do ib = space_charge_com%n_bin/2, space_charge_com%n_bin
  slice => csr%slice(ib)
  if (slice%n_particle < space_charge_com%sc_min_in_bin) then
    slice%sig_x = csr%slice(ib-1)%sig_x
    slice%sig_y = csr%slice(ib-1)%sig_y
  else
    if (slice%sig_x < sig_x_ave * space_charge_com%lsc_sigma_cutoff) slice%sig_x = sig_x_ave * space_charge_com%lsc_sigma_cutoff
    if (slice%sig_y < sig_y_ave * space_charge_com%lsc_sigma_cutoff) slice%sig_y = sig_y_ave * space_charge_com%lsc_sigma_cutoff
  endif
enddo

do ib = space_charge_com%n_bin/2, 1, -1
  slice => csr%slice(ib)
  if (slice%n_particle < space_charge_com%sc_min_in_bin) then
    slice%sig_x = csr%slice(ib+1)%sig_x
    slice%sig_y = csr%slice(ib+1)%sig_y
  else
    if (slice%sig_x < sig_x_ave * space_charge_com%lsc_sigma_cutoff) slice%sig_x = sig_x_ave * space_charge_com%lsc_sigma_cutoff
    if (slice%sig_y < sig_y_ave * space_charge_com%lsc_sigma_cutoff) slice%sig_y = sig_y_ave * space_charge_com%lsc_sigma_cutoff
  endif
enddo


end subroutine csr_bin_particles

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_kicks (ele, ds_kick_pt, csr, err_flag)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Input:
!   ele          -- ele_struct: Element being tracked through.
!   ds_kick_pt   -- real(rp): Distance between the beginning of the element we are
!                    tracking through and the kick point (which is within this element).
!   csr          -- csr_struct: 
!
! Output:
!   csr          -- csr_struct: 
!     %kick1(:)          -- CSR kick calculation bin array. 
!     %slice(:)%kick_csr -- Integrated kick
!   err_flag     -- logical: Set True if there is an error. False otherwise
!-

subroutine csr_bin_kicks (ele, ds_kick_pt, csr, err_flag)

implicit none

type (ele_struct) ele
type (csr_struct), target :: csr
type (branch_struct), pointer :: branch
type (csr_kick1_struct), pointer :: kick1

real(rp) ds_kick_pt, coef, dr_match(3)

integer i, n_bin, m_fft
logical err_flag

! Contiguous arrays for vectorized convolution
real(rp), allocatable :: I_int_arr(:), edge_dcharge_arr(:), image_kick_arr(:), charge_arr(:)

! FFT workspace for O(n_bin * log(n_bin)) convolution
complex(rp), allocatable :: fft_a(:), fft_b(:)

character(16) :: r_name = 'csr_bin_kicks'

! The kick point P is fixed.
! Loop over all kick1 bins and compute the kick.
! When y_source == 0 (no image charges), only positive bins contribute to CSR kick
! (source behind kicked particle). Negative bins have dz <= 0 and I_csr returns 0.

err_flag = .false.

if (csr%y_source == 0) then
  ! CSR only: skip negative bins since I_csr = 0 for dz_particles <= 0
  do i = 1, ubound(csr%kick1, 1)
    kick1 => csr%kick1(i)
    kick1%dz_particles = i * csr%dz_slice

    if (i == 1) then
      kick1%ix_ele_source = csr%ix_ele_kick
      dr_match = 0
    else
      kick1%ix_ele_source = csr%kick1(i-1)%ix_ele_source
    endif

    kick1%s_chord_source = s_source_calc(kick1, csr, err_flag, dr_match)
    if (err_flag) return
    call I_csr (kick1, i, csr)
  enddo

else
  ! Image charges: need full range for both positive and negative separations
  do i = lbound(csr%kick1, 1), ubound(csr%kick1, 1)
    kick1 => csr%kick1(i)
    kick1%dz_particles = i * csr%dz_slice

    if (i == lbound(csr%kick1, 1)) then
      kick1%ix_ele_source = csr%ix_ele_kick
      dr_match = 0
    else
      kick1%ix_ele_source = csr%kick1(i-1)%ix_ele_source
    endif

    kick1%s_chord_source = s_source_calc(kick1, csr, err_flag, dr_match)
    if (err_flag) return
    call image_charge_kick_calc (kick1, csr)
  enddo
endif

!

coef = csr%actual_track_step * classical_radius(csr%species) / &
                              (csr%rel_mass * e_charge * abs(charge_of(csr%species)) * csr%gamma)
n_bin = space_charge_com%n_bin

! CSR & Image charge kick
! Use contiguous arrays for vectorized convolution instead of strided struct access.

if (csr%y_source == 0) then
  if (ele%csr_method == one_dim$) then
    allocate(I_int_arr(n_bin), edge_dcharge_arr(n_bin))
    do i = 1, n_bin
      I_int_arr(i) = csr%kick1(i)%I_int_csr
      edge_dcharge_arr(i) = csr%slice(i)%edge_dcharge_density_dz
    enddo

    ! FFT-based linear convolution: O(n_bin * log(n_bin)) instead of O(n_bin^2).
    ! kick_csr(i) = coef * sum_{k=1}^{i} I_int_arr(k) * edge_dcharge_arr(i+1-k)
    ! This equals the first n_bin elements of the linear convolution of I_int_arr with edge_dcharge_arr.
    ! Zero-pad to avoid circular convolution artifacts.
    m_fft = 1
    do while (m_fft < 2 * n_bin)
      m_fft = m_fft * 2
    enddo

    allocate(fft_a(m_fft), fft_b(m_fft))
    fft_a = (0.0_rp, 0.0_rp)
    fft_b = (0.0_rp, 0.0_rp)
    do i = 1, n_bin
      fft_a(i) = cmplx(I_int_arr(i), 0.0_rp, rp)
      fft_b(i) = cmplx(edge_dcharge_arr(i), 0.0_rp, rp)
    enddo

    call fft_1d(fft_a, -1)
    call fft_1d(fft_b, -1)
    fft_a = fft_a * fft_b
    call fft_1d(fft_a, 1)

    do i = 1, n_bin
      csr%slice(i)%kick_csr = coef * real(fft_a(i), rp) / m_fft
    enddo

    deallocate(fft_a, fft_b, I_int_arr, edge_dcharge_arr)
  endif

else  ! Image charge
  allocate(image_kick_arr(-n_bin:n_bin), charge_arr(n_bin))
  do i = -n_bin, n_bin
    image_kick_arr(i) = csr%kick1(i)%image_kick_csr
  enddo
  do i = 1, n_bin
    charge_arr(i) = csr%slice(i)%charge
  enddo
  do i = 1, n_bin
    csr%slice(i)%kick_csr = csr%slice(i)%kick_csr + coef * &
                  dot_product(image_kick_arr(i-1:i-n_bin:-1), charge_arr(1:n_bin))
  enddo
  deallocate(image_kick_arr, charge_arr)
endif

! Longitudinal space charge kick
! Note: image charges do not contribute to LSC.

if (ele%space_charge_method == slice$ .and. csr%y_source == 0) then
  call lsc_kick_params_calc (ele, csr)
endif

end subroutine csr_bin_kicks
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function s_source_calc (kick1, csr, err_flag, dr_match) result (s_source)
!
! Routine to calculate the distance between source and kick points.
!
! Input:
!   kick1         -- csr_kick1_struct:
!   csr           -- csr_struct:
!   dr_match(3)   -- real(rp): Discontinuity factor if there is a match element between source and kick elements.
!
! Output:
!   s_source      -- real(rp): source s-position.
!   kick1         -- csr_kick1_struct:
!   err_flag      -- logical: Set True if there is an error. Untouched otherwise.
!   dr_match(3)   -- real(rp): Discontinuity factor if there is a match element between source and kick elements.
!-

function s_source_calc (kick1, csr, err_flag, dr_match) result (s_source)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_ele_info_struct), pointer :: einfo_s, einfo_k
type (floor_position_struct), pointer :: fk, f0, fs
type (ele_struct), pointer :: ele

real(rp) a, b, c, dz, s_source, beta2, L0, Lz, ds_source
real(rp) z0, z1, sz_kick, sz0, Lsz0, ddz0, ddz1, dr_match(3)

integer i, last_step, status
logical err_flag

character(*), parameter :: r_name = 's_source_calc'

! Each interation of the loop looks for a possible source point in the lattice element with 
! index kick1%ix_ele_source. If found, return. If not, move on to another element

dz = kick1%dz_particles   ! Target distance.
beta2 = csr%beta**2
last_step = 0
s_source = 0        ! To prevent uninitalized complaints from the compiler.

do

  einfo_s => csr%eleinfo(kick1%ix_ele_source)

  ! If at beginning of lattice assume an infinite drift.
  ! s_source will be negative

  if (kick1%ix_ele_source == 0) then
    fk => csr%floor_k
    fs => kick1%floor_s
    f0 => einfo_s%floor1

    L0 = sqrt((fk%r(1) - f0%r(1))**2 + (fk%r(3) - f0%r(3))**2 + csr%y_source**2)
    ! L_z is the z-component from lat start to the kick point.
    Lz = (fk%r(1) - f0%r(1)) * sin(f0%theta) + (fk%r(3) - f0%r(3)) * cos(f0%theta) 
    ! Lsz0 is Ls from the lat start to the kick point
    Lsz0 = dspline_len(0.0_rp, csr%s_chord_kick, csr%eleinfo(csr%ix_ele_kick)%spline) + csr%s_chord_kick
    do i = 1, csr%ix_ele_kick - 1
      Lsz0 = Lsz0 + csr%eleinfo(i)%dL_s + csr%eleinfo(i)%L_chord
    enddo

    a = 1/csr%gamma2
    b = 2 * (Lsz0 - dz - beta2 * Lz)
    c = (Lsz0 - dz)**2 - beta2 * L0**2
    ds_source = -(-b + sqrt(b**2 - 4 * a * c)) / (2 * a)
    s_source = einfo_s%ele%s + ds_source

    fs%r = [f0%r(1) + ds_source * sin(f0%theta), csr%y_source, f0%r(3) + ds_source * cos(f0%theta)]
    fs%theta = f0%theta
    kick1%L_vec = fk%r - fs%r
    kick1%L = sqrt(dot_product(kick1%L_vec, kick1%L_vec))
    kick1%dL = lsz0 - ds_source - kick1%L  ! Remember s_source is negative
    kick1%theta_L = atan2(kick1%L_vec(1), kick1%L_vec(3))
    kick1%theta_sl = f0%theta - kick1%theta_L
    einfo_k => csr%eleinfo(csr%ix_ele_kick)
    kick1%theta_lk = kick1%theta_L - (spline1(einfo_k%spline, csr%s_chord_kick, 1) + einfo_k%theta_chord)
    return
  endif

  ! Match elements can have an orbit discontinuity which is non-physical.
  ! So if there is a match element then shift the orbit using dr_match to remove the discontinuity.
  ! Any angle discontinuity is ignored.

  ele => einfo_s%ele

  if (ele%key == floor_shift$) then
    call out_io (s_fatal$, r_name, &
        'CSR CALC IS NOT ABLE TO HANDLE A FLOOR_SHIFT ELEMENT: ' // ele%name, &
        'WHILE TRACKING THROUGH ELEMENT: ' // csr%kick_ele%name)
    if (global_com%exit_on_error) call err_exit
    err_flag = .true.
    return
  endif

  ! Look at ends of the source element and check if the source point is within the element or not.
  ! Generally dz decreases with increasing s but this may not be true for patch elements.

  ddz0 = ddz_calc_csr(0.0_rp, status)
  ddz1 = ddz_calc_csr(einfo_s%L_chord, status)

  if (last_step == -1 .and. ddz1 > 0) then  ! Roundoff error is causing ddz1 to be positive.
    s_source = ddz_calc_csr(einfo_s%L_chord, status)
    return
  endif

  if (last_step == 1 .and. ddz0 < 0) then  ! Roundoff error is causing ddz0 to be negative.
    s_source = ddz_calc_csr(0.0_rp, status)
    return
  endif

  if (ddz0 < 0 .and. ddz1 < 0) then
    if (last_step == 1) exit       ! Round off error can cause problems
    last_step = -1
    kick1%ix_ele_source = kick1%ix_ele_source - 1

    einfo_s => csr%eleinfo(kick1%ix_ele_source)
    if (einfo_s%ele%key == match$) then
      dr_match = einfo_s%floor1%r - einfo_s%floor0%r ! discontinuity in x.
      kick1%ix_ele_source = kick1%ix_ele_source - 1
    endif

    cycle
  endif

  if (ddz0 > 0 .and. ddz1 > 0) then
    if (kick1%ix_ele_source == csr%ix_ele_kick) return  ! Source ahead of kick pt => ignore.
    if (last_step == -1) exit      ! Round off error can cause problems
    last_step = 1
    kick1%ix_ele_source = kick1%ix_ele_source + 1

    if (csr%eleinfo(kick1%ix_ele_source)%ele%key == match$) then
      dr_match = 0
      kick1%ix_ele_source = kick1%ix_ele_source + 1
    endif

    cycle
  endif

  ! Only possibility left is that root is bracketed.

  s_source = super_zbrent (ddz_calc_csr, 0.0_rp, einfo_s%L_chord, 1e-12_rp, 1e-8_rp, status)
  return
    
enddo

call out_io (s_fatal$, r_name, &
    'CSR CALCULATION ERROR. PLEASE REPORT THIS...', &
    'WHILE TRACKING THROUGH ELEMENT: ', csr%kick_ele%name)
if (global_com%exit_on_error) call err_exit
err_flag = .true.

!----------------------------------------------------------------------------
contains

!+
! Function ddz_calc_csr (s_chord_source, status) result (ddz_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle at constant time minus the target distance.
!
! Input:
!   s_chord_source  -- real(rp): Chord distance from start of element.
!
! Output:
!   ddz_this        -- real(rp): Distance between source and kick particles: Calculated - Wanted.
!-

function ddz_calc_csr (s_chord_source, status) result (ddz_this)

implicit none

type (csr_ele_info_struct), pointer :: ce
real(rp), intent(in) :: s_chord_source
real(rp) ddz_this, x, z, c, s, dtheta_L
real(rp) s0, s1, ds, theta_L, dL

integer status
integer i

character(*), parameter :: r_name = 'ddz_calc_csr'

! 

x = spline1(einfo_s%spline, s_chord_source)
c = cos(einfo_s%theta_chord)
s = sin(einfo_s%theta_chord)
kick1%floor_s%r = [x*c + s_chord_source*s, csr%y_source, -x*s + s_chord_source*c] + einfo_s%floor0%r + dr_match   ! Floor coordinates

kick1%L_vec = csr%floor_k%r - kick1%floor_s%r
kick1%L = sqrt(dot_product(kick1%L_vec, kick1%L_vec))
kick1%theta_L = atan2(kick1%L_vec(1), kick1%L_vec(3))

s0 = s_chord_source
s1 = csr%s_chord_kick

if (kick1%ix_ele_source == csr%ix_ele_kick) then
  ds = s1 - s0
  ! dtheta_L = angle of L line in centroid chord ref frame
  dtheta_L = einfo_s%spline%coef(1) + einfo_s%spline%coef(2) * (2*s0 + ds) + einfo_s%spline%coef(3) * (3*s0**2 + 3*s0*ds + ds**2)
  dL = dspline_len(s0, s1, einfo_s%spline, dtheta_L) ! = Ls - L
  ! Ls is negative if the source pt is ahead of the kick pt (ds < 0). But L is always positive. 
  if (ds < 0) dL = dL + 2 * ds    ! Correct for L always being positive.
  kick1%theta_sl = spline1(einfo_s%spline, s0, 1) - dtheta_L
  kick1%theta_lk = dtheta_L - spline1(einfo_s%spline, s1, 1)

else
  ! In an element where the beam centroid takes a backstep, %theta_chord will be anti-parallel to
  ! other angles. This is why pi/2 is used with modulo2 to make sure angles are in the range [-pi/2, pi/2]
  theta_L = kick1%theta_L
  dL = dspline_len(s0, einfo_s%L_chord, einfo_s%spline, modulo2(theta_L-einfo_s%theta_chord, pi/2))
  do i = kick1%ix_ele_source+1, csr%ix_ele_kick-1
    ce => csr%eleinfo(i)
    if (ce%ele%key == match$) cycle  ! Match elements are adjusted to give zero displacement.
    dL = dL + dspline_len(0.0_rp, ce%L_chord, ce%spline, modulo2(theta_L-ce%theta_chord, pi/2))
  enddo
  ce => csr%eleinfo(csr%ix_ele_kick)
  dL = dL + dspline_len(0.0_rp, s1, ce%spline, modulo2(theta_L-ce%theta_chord, pi/2))
  kick1%theta_sl = modulo2((spline1(einfo_s%spline, s0, 1) + einfo_s%theta_chord) - theta_L, pi/2)
  kick1%theta_lk = modulo2(theta_L - (spline1(ce%spline, s1, 1) + ce%theta_chord), pi/2)
endif

kick1%floor_s%theta = kick1%theta_sl + kick1%theta_L

! The above calc for dL neglected csr%y_source. So must correct for this.

if (csr%y_source /= 0) dL = dL - (kick1%L - sqrt(kick1%L_vec(1)**2 + kick1%L_vec(3)**2))
kick1%dL = dL
ddz_this = kick1%L / (2 * csr%gamma2) + dL
ddz_this = ddz_this - kick1%dz_particles

end function ddz_calc_csr

end function s_source_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_kick_params_calc (ele, csr)
!
! Routine to cache intermediate values needed for the lsc calculation.
! This routine is not for image currents.
!
! Input:
!   ele       -- Element_struct: Element to set up cache for.
!   csr       -- csr_struct: 
!     %slice(:)   -- bin array of particle averages.
!
! Output:
!   csr     -- csr_struct: Binned particle averages.
!     %slice(:)   -- bin array of particle averages.
!-

subroutine lsc_kick_params_calc (ele, csr)

use da2_mod

implicit none

type (ele_struct) ele
type (csr_struct), target :: csr

real(rp), pointer :: f(:,:)
real(rp) factor, f00, c_val, dz_half, z_center_i
real(rp) sx2, sy2, g_z1_s, g_z2_s, h_z1_s, h_z2_s, alph, bet, ss, f1(0:2,0:2)

integer i, j, n_bin

! Shared read-only arrays (source slice properties, independent of kick index i)
real(rp), allocatable :: sx_v(:), sy_v(:), a_v(:), b_v(:)
real(rp), allocatable :: charge_v(:), dcdz_v(:), z_center_v(:)
real(rp), allocatable :: radix_v(:), sr_v(:)

! Thread-private work arrays (recomputed each iteration of the i loop)
real(rp), allocatable :: abs_z_v(:), z1_v(:), z2_v(:), rho0_v(:), drho_v(:)
real(rp), allocatable :: b2cz1_v(:), b2cz2_v(:), abcz1_v(:), abcz2_v(:)
real(rp), allocatable :: atz1_v(:), atz2_v(:), bcd_v(:), dk0_v(:)
integer, allocatable :: sign_v(:)

character(*), parameter :: r_name = 'lsc_kick_params_calc'

!

if (ele%space_charge_method /= slice$) return

n_bin = space_charge_com%n_bin

factor = csr%kick_factor * csr%actual_track_step * classical_radius(csr%species) / &
                              (csr%rel_mass * e_charge * abs(charge_of(csr%species)) * csr%gamma)

c_val = csr%gamma**2
dz_half = csr%dz_slice / 2

! Precompute j-dependent arrays (source slice properties, independent of kick index i)
! These are shared read-only across threads.

allocate(sx_v(n_bin), sy_v(n_bin), a_v(n_bin), b_v(n_bin))
allocate(charge_v(n_bin), dcdz_v(n_bin), z_center_v(n_bin))
allocate(radix_v(n_bin), sr_v(n_bin))

do j = 1, n_bin
  sx_v(j) = csr%slice(j)%sig_x
  sy_v(j) = csr%slice(j)%sig_y
  charge_v(j) = csr%slice(j)%charge
  dcdz_v(j) = csr%slice(j)%dcharge_density_dz
  z_center_v(j) = csr%slice(j)%z_center
enddo

! Quantities that depend only on the source slice (j), not on kick slice (i)
a_v = sx_v * sy_v
b_v = csr%gamma * (sx_v**2 + sy_v**2) / (sx_v + sy_v)
radix_v = -b_v**2 + 4 * a_v * c_val
sr_v = sqrt(abs(radix_v))

! Compute the kick at the center of each bin (parallelized).
! i = index of slice where kick is computed.
! Each iteration writes only to csr%slice(i)%kick_lsc, so the outer loop is parallelizable.
! Thread-private work arrays are allocated once per thread (inside the parallel region, outside the do loop).

!$OMP parallel default(none) &
!$OMP   shared(n_bin, csr, z_center_v, charge_v, dcdz_v, a_v, b_v, radix_v, sr_v, c_val, dz_half, factor) &
!$OMP   private(i, j, z_center_i, abs_z_v, z1_v, z2_v, rho0_v, drho_v, bcd_v, abcz1_v, abcz2_v, &
!$OMP           b2cz1_v, b2cz2_v, atz1_v, atz2_v, dk0_v, sign_v)
allocate(abs_z_v(n_bin), z1_v(n_bin), z2_v(n_bin), rho0_v(n_bin), drho_v(n_bin))
allocate(b2cz1_v(n_bin), b2cz2_v(n_bin), abcz1_v(n_bin), abcz2_v(n_bin))
allocate(atz1_v(n_bin), atz2_v(n_bin), bcd_v(n_bin), dk0_v(n_bin))
allocate(sign_v(n_bin))

!$OMP do
do i = 1, n_bin
  z_center_i = csr%slice(i)%z_center
  abs_z_v = abs(z_center_i - z_center_v)
  do j = 1, n_bin
    if (z_center_i > z_center_v(j)) then
      sign_v(j) = 1
    elseif (z_center_i < z_center_v(j)) then
      sign_v(j) = -1
    else
      sign_v(j) = 0
    endif
  enddo
  z1_v = abs_z_v - dz_half
  z2_v = abs_z_v + dz_half
  drho_v = dcdz_v * sign_v
  rho0_v = charge_v / csr%dz_slice - drho_v * abs_z_v

  ! Self-slice override (diagonal: i == j)
  drho_v(i) = dcdz_v(i)
  rho0_v(i) = 0
  z1_v(i) = 0
  z2_v(i) = dz_half
  sign_v(i) = -2        ! Factor of 2 accounts for 1/2 we did not integrate over.

  bcd_v = 2 * c_val * rho0_v - b_v * drho_v
  abcz1_v = a_v + b_v * z1_v + c_val * z1_v**2
  abcz2_v = a_v + b_v * z2_v + c_val * z2_v**2
  b2cz1_v = b_v + 2 * c_val * z1_v
  b2cz2_v = b_v + 2 * c_val * z2_v
  atz1_v = atan(b2cz1_v / sr_v) / sr_v
  atz2_v = atan(b2cz2_v / sr_v) / sr_v
  do j = 1, n_bin
    if (radix_v(j) <= 0) then
      atz1_v(j) = log((b2cz1_v(j) - sr_v(j)) / (b2cz1_v(j) + sr_v(j))) / (2 * sr_v(j))
      atz2_v(j) = log((b2cz2_v(j) - sr_v(j)) / (b2cz2_v(j) + sr_v(j))) / (2 * sr_v(j))
    endif
  enddo
  dk0_v = factor * ((2 * atz2_v * bcd_v + drho_v * log(abcz2_v)) - &
                     (2 * atz1_v * bcd_v + drho_v * log(abcz1_v))) / (2 * c_val)
  csr%slice(i)%kick_lsc = csr%slice(i)%kick_lsc + sum(sign_v * dk0_v)
enddo
!$OMP end do
deallocate(abs_z_v, z1_v, z2_v, rho0_v, drho_v)
deallocate(b2cz1_v, b2cz2_v, abcz1_v, abcz2_v, atz1_v, atz2_v, bcd_v, dk0_v)
deallocate(sign_v)
!$OMP end parallel

! Transverse dependence DA2 coefficients.
! Each outer iteration i writes only to csr%slice(i)%coef_lsc_plus/minus, so i is parallelizable.
! The inner j loop has a serial dependency (harmonic-mean accumulation), but different i are independent.

if (space_charge_com%lsc_kick_transverse_dependence) then
  !$OMP parallel default(none) &
  !$OMP   shared(n_bin, csr, z_center_v, charge_v, dcdz_v, a_v, b_v, radix_v, sr_v, sx_v, sy_v, c_val, dz_half, factor) &
  !$OMP   private(i, j, z_center_i, abs_z_v, z1_v, z2_v, rho0_v, drho_v, bcd_v, abcz1_v, abcz2_v, &
  !$OMP           b2cz1_v, b2cz2_v, atz1_v, atz2_v, dk0_v, sign_v, &
  !$OMP           f, f1, f00, sx2, sy2, g_z1_s, g_z2_s, h_z1_s, h_z2_s, alph, bet, ss)
  allocate(abs_z_v(n_bin), z1_v(n_bin), z2_v(n_bin), rho0_v(n_bin), drho_v(n_bin))
  allocate(b2cz1_v(n_bin), b2cz2_v(n_bin), abcz1_v(n_bin), abcz2_v(n_bin))
  allocate(atz1_v(n_bin), atz2_v(n_bin), bcd_v(n_bin), dk0_v(n_bin))
  allocate(sign_v(n_bin))

  !$OMP do
  do i = 1, n_bin
    z_center_i = csr%slice(i)%z_center

    abs_z_v = abs(z_center_i - z_center_v)

    do j = 1, n_bin
      if (z_center_i > z_center_v(j)) then
        sign_v(j) = 1
      elseif (z_center_i < z_center_v(j)) then
        sign_v(j) = -1
      else
        sign_v(j) = 0
      endif
    enddo

    z1_v = abs_z_v - dz_half
    z2_v = abs_z_v + dz_half
    drho_v = dcdz_v * sign_v
    rho0_v = charge_v / csr%dz_slice - drho_v * abs_z_v

    drho_v(i) = dcdz_v(i)
    rho0_v(i) = 0
    z1_v(i) = 0
    z2_v(i) = dz_half
    sign_v(i) = -2

    bcd_v = 2 * c_val * rho0_v - b_v * drho_v
    abcz1_v = a_v + b_v * z1_v + c_val * z1_v**2
    abcz2_v = a_v + b_v * z2_v + c_val * z2_v**2
    b2cz1_v = b_v + 2 * c_val * z1_v
    b2cz2_v = b_v + 2 * c_val * z2_v

    atz1_v = atan(b2cz1_v / sr_v) / sr_v
    atz2_v = atan(b2cz2_v / sr_v) / sr_v

    do j = 1, n_bin
      if (radix_v(j) <= 0) then
        atz1_v(j) = log((b2cz1_v(j) - sr_v(j)) / (b2cz1_v(j) + sr_v(j))) / (2 * sr_v(j))
        atz2_v(j) = log((b2cz2_v(j) - sr_v(j)) / (b2cz2_v(j) + sr_v(j))) / (2 * sr_v(j))
      endif
    enddo

    dk0_v = factor * ((2 * atz2_v * bcd_v + drho_v * log(abcz2_v)) - &
                       (2 * atz1_v * bcd_v + drho_v * log(abcz1_v))) / (2 * c_val)

    do j = 1, n_bin
      if (dk0_v(j) == 0) cycle

      f00 = 1 / dk0_v(j)
      sx2 = sx_v(j)**2
      sy2 = sy_v(j)**2

      g_z1_s = a_v(j) * (b2cz1_v(j)*bcd_v(j) + 4*abcz1_v(j)*atz1_v(j)*bcd_v(j)*c_val - drho_v(j)*radix_v(j)) / &
               (2*abcz1_v(j)*c_val*radix_v(j))
      g_z2_s = a_v(j) * (b2cz2_v(j)*bcd_v(j) + 4*abcz2_v(j)*atz2_v(j)*bcd_v(j)*c_val - drho_v(j)*radix_v(j)) / &
               (2*abcz2_v(j)*c_val*radix_v(j))
      h_z1_s = (b_v(j)*b_v(j)*b2cz1_v(j)*bcd_v(j) - 6*a_v(j)*b2cz1_v(j)*bcd_v(j)*c_val - &
                4*abcz1_v(j)*b2cz1_v(j)*bcd_v(j)*c_val - 24*abcz1_v(j)**2*atz1_v(j)*bcd_v(j)*c_val**2 + &
                drho_v(j)*radix_v(j)**2 - 2*b_v(j)*b2cz1_v(j)*bcd_v(j)*c_val*z1_v(j) - &
                2*b2cz1_v(j)*bcd_v(j)*(c_val*z1_v(j))**2)
      h_z2_s = (b_v(j)*b_v(j)*b2cz2_v(j)*bcd_v(j) - 6*a_v(j)*b2cz2_v(j)*bcd_v(j)*c_val - &
                4*abcz2_v(j)*b2cz2_v(j)*bcd_v(j)*c_val - 24*abcz2_v(j)**2*atz2_v(j)*bcd_v(j)*c_val**2 + &
                drho_v(j)*radix_v(j)**2 - 2*b_v(j)*b2cz2_v(j)*bcd_v(j)*c_val*z2_v(j) - &
                2*b2cz2_v(j)*bcd_v(j)*(c_val*z2_v(j))**2)

      alph = f00**2 * (g_z2_s - g_z1_s)
      bet = f00**3 * (g_z2_s - g_z1_s)**2 + (f00 * a_v(j))**2 * &
            (h_z2_s/abcz2_v(j)**2 - h_z1_s/abcz1_v(j)**2) / (4 * c_val * radix_v(j)**2)

      ss = 1.0_rp / sign_v(j)
      f1(0,0) = ss * f00
      f1(0,1) = ss * alph / (2 * sy2)
      f1(0,2) = ss * (alph + 2*bet) / (8 * sy2**2)
      f1(1,0) = ss * alph / (2 * sx2)
      f1(1,1) = ss * (alph + 2*bet) / (4 * sx2 * sy2)
      f1(2,0) = ss * (alph + 2*bet) / (8 * sx2**2)

      if (f1(0,0) > 0) then
        f => csr%slice(i)%coef_lsc_plus
      else
        f => csr%slice(i)%coef_lsc_minus
      endif

      if (f(0,0) == 0) then  ! If first time
        f = f1
      else
        f = da2_div(da2_mult(f1, f), f + f1)
      endif
    enddo

  enddo
  !$OMP end do
  deallocate(abs_z_v, z1_v, z2_v, rho0_v, drho_v)
  deallocate(b2cz1_v, b2cz2_v, abcz1_v, abcz2_v, atz1_v, atz2_v, bcd_v, dk0_v)
  deallocate(sign_v)
  !$OMP end parallel
endif

deallocate(sx_v, sy_v, a_v, b_v, charge_v, dcdz_v, z_center_v)
deallocate(radix_v, sr_v)

end subroutine lsc_kick_params_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine I_csr (kick1, i_bin, csr) 
!
! Routine to calculate the CSR kick integral (at y = 0)
!
! Input:
!   kick1      -- csr_kick1_struct: 
!   i_bin      -- integer: Bin index.
!   csr    -- csr_struct:
!
! Output:
!   kick1     -- csr_kick1_struct: 
!     %I_csr     -- real(rp): CSR kick integral.
!     %I_int_csr -- real(rp): Integral of I_csr. Only calculated for i_bin = 1 since it is not needed otherwise.
!-

subroutine I_csr (kick1, i_bin, csr)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_kick1_struct), pointer :: k
type (csr_ele_info_struct), pointer :: eleinfo
type (spline_struct), pointer :: spl

real(rp) g_bend, z, zz, Ls, L, dtheta_L, dL, s_chord_kick
integer i_bin, ix_ele_kick

!

kick1%I_int_csr = 0
kick1%I_csr = 0

! No kick when source particle is ahead of the kicked particle

z = kick1%dz_particles
if (z <= 0) return

!

k => kick1
k%I_csr = -csr%kick_factor * 2 * (k%dL / z + csr%gamma2 * k%theta_sl * k%theta_lk / (1 + csr%gamma2 * k%theta_sl**2)) / k%L

! I_csr Integral 

if (i_bin == 1) then
  ix_ele_kick = csr%ix_ele_kick
  s_chord_kick = csr%s_chord_kick
  do
    if (s_chord_kick /= 0) exit  ! Not at edge of element and element has finite length
    ix_ele_kick = ix_ele_kick - 1
    if (ix_ele_kick == 0) return  ! No kick from before beginning of lattice
    s_chord_kick = csr%eleinfo(ix_ele_kick)%L_chord
  enddo

  spl => csr%eleinfo(ix_ele_kick)%spline
  g_bend = -spline1(spl, s_chord_kick, 2) / sqrt(1 + spline1(spl, s_chord_kick, 1)**2)**3

  if (k%ix_ele_source == ix_ele_kick) then
    Ls = k%L + k%dL
    k%I_int_csr = -csr%kick_factor * ((g_bend * Ls/2)**2 - log(2 * csr%gamma2 * z / Ls) / csr%gamma2)  
  else
    ! Since source pt is in another element, split integral into pieces.
    ! First integrate over element containing the kick point.
    L = s_chord_kick
    eleinfo => csr%eleinfo(ix_ele_kick)
    dtheta_L = eleinfo%spline%coef(1) + eleinfo%spline%coef(2) * L + eleinfo%spline%coef(3) * L**2
    dL = dspline_len(0.0_rp, L, eleinfo%spline, dtheta_L) ! = Ls - L
    Ls = L + dL
    zz = L / (2 * csr%gamma2) + dL
    k%I_int_csr = -csr%kick_factor * ((g_bend * Ls/2)**2 - log(2 * csr%gamma2 * zz / Ls) / csr%gamma2)
    ! Now add on rest of integral
    k%I_int_csr = k%I_int_csr + k%I_csr * (z - zz)
  endif

else
  kick1%I_int_csr = (kick1%I_csr + csr%kick1(i_bin-1)%I_csr) * csr%dz_slice / 2
endif

end subroutine I_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine image_charge_kick_calc (kick1, csr) 
!
! Routine to calculate the image charge kick.
!
! Input:
!   kick1    -- csr_kick1_struct: 
!   csr      -- csr_struct:
!
! Output:
!   kick1             -- csr_kick1_struct: 
!     %image_kick_csr -- real(rp): Image charge kick.
!-

subroutine image_charge_kick_calc (kick1, csr)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_kick1_struct), pointer :: k
type (spline_struct), pointer :: sp

real(rp) N_vec(3), G_vec(3), B_vec(3), Bp_vec(3), NBp_vec(3), NBpG_vec(3), rad_cross_vec(3)
real(rp) z, sin_phi, cos_phi, OneNBp, OneNBp3, radiate, coulomb1, theta, g_bend

!

k => kick1
sp => csr%eleinfo(k%ix_ele_source)%spline

g_bend = -spline1(sp, k%s_chord_source, 2) / sqrt(1 + spline1(sp, k%s_chord_source, 1)**2)**3
theta = k%floor_s%theta

Bp_vec = csr%beta * [sin(theta), 0.0_rp, cos(theta) ]             ! beta vector at source point
G_vec = csr%beta**2 * g_bend * [-cos(theta), 0.0_rp, sin(theta) ] ! Acceleration vector

k%image_kick_csr = 0
z = k%dz_particles

N_vec = (csr%floor_k%r - k%floor_s%r) / k%L

theta = csr%floor_k%theta
B_vec = [sin(theta), 0.0_rp, cos(theta)]        ! Beta vector at kicked point


OneNBp = 1 - sum(N_vec * Bp_vec)
OneNBp3 = OneNBp**3

NBp_vec = N_vec - Bp_vec
NBpG_vec = cross_product(NBp_vec, G_vec)
rad_cross_vec = cross_product(N_vec, NBpG_vec)

radiate  = dot_product (B_vec, rad_cross_vec) / (k%L * OneNBp3)
coulomb1 = dot_product (B_vec, NBp_vec) / (csr%gamma2 * k%L**2 * OneNBp3)
kick1%image_kick_csr = csr%kick_factor * (radiate + coulomb1)

end subroutine image_charge_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_and_sc_apply_kicks (ele, csr, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
! Input:
!   ele         -- ele_struct: Element being tracked through.
!   csr         -- csr_struct: 
!   particle(:) -- coord_struct: Particles to kick.
!
! Output:
!   particle(:) -- Coord_struct: Particles with kick applied.
!-

subroutine csr_and_sc_apply_kicks (ele, csr, particle)

use da2_mod

implicit none

type (ele_struct) ele
type (csr_struct), target :: csr
type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p
type (csr_bunch_slice_struct), pointer :: slice

real(rp) zp, r1, r0, dz, dpz, nk(2), dnk(2,2), f0, f, beta0, dpx, dpy, x, y, f_tot
real(rp) Evec(3), factor, pz0, dct_ave, new_beta, ef
integer i, n, i0, i_del, ip

! CSR kick and Slice space charge kick
! We use a weighted average so that the integral varies smoothly as a function of particle%vec(5).

if (ele%csr_method == one_dim$ .or. ele%space_charge_method == slice$) then

  ! CSR-only kick (no space charge): each particle reads shared slice data and updates only its own vec(6).
  if (ele%csr_method == one_dim$ .and. ele%space_charge_method /= slice$) then
    !$OMP parallel do default(none) shared(particle, csr, space_charge_com, global_com) private(ip, zp, i0, r1, r0)
    do ip = 1, size(particle)
      if (particle(ip)%state /= alive$) cycle
      zp = particle(ip)%vec(5)
      i0 = int((zp - csr%slice(1)%z_center) / csr%dz_slice) + 1
      r1 = (zp - csr%slice(i0)%z_center) / csr%dz_slice
      r0 = 1 - r1
      if (r1 < -0.01_rp .or. r1 > 1.01_rp .or. i0 < 1 .or. i0 >= space_charge_com%n_bin) then
        !$OMP critical
        call out_io (s_error$, 'csr_and_sc_apply_kicks', 'CSR INTERNAL ERROR!')
        !$OMP end critical
        if (global_com%exit_on_error) call err_exit
      endif
      particle(ip)%vec(6) = particle(ip)%vec(6) + r0 * csr%slice(i0)%kick_csr + r1 * csr%slice(i0+1)%kick_csr
    enddo
    !$OMP end parallel do

  else
  ! General case with possible space charge
  !$OMP parallel do default(none) &
  !$OMP   shared(particle, csr, ele, space_charge_com, global_com) &
  !$OMP   private(ip, p, slice, zp, i0, r1, r0, dpz, x, y, beta0, f0, f, dpx, dpy, nk, dnk)
  do ip = 1, size(particle)
    p => particle(ip)
    if (p%state /= alive$) cycle
    zp = p%vec(5)
    i0 = int((zp - csr%slice(1)%z_center) / csr%dz_slice) + 1
    r1 = (zp - csr%slice(i0)%z_center) / csr%dz_slice
    r0 = 1 - r1

    ! r1 should be in [0,1] but allow for some round-off error
    if (r1 < -0.01_rp .or. r1 > 1.01_rp .or. i0 < 1 .or. i0 >= space_charge_com%n_bin) then
      !$OMP critical
      call out_io (s_error$, 'csr_and_sc_apply_kicks', 'CSR INTERNAL ERROR!')
      !$OMP end critical
      if (global_com%exit_on_error) call err_exit
    endif

    ! CSR kick
    ! We use a weighted average so that the integral varies smoothly as a function of particle%vec(5).

    if (ele%csr_method == one_dim$) then
      p%vec(6) = p%vec(6) + r0 * csr%slice(i0)%kick_csr + r1 * csr%slice(i0+1)%kick_csr
    endif

  ! Slice space charge kick

  if (ele%space_charge_method == slice$) then
    if (space_charge_com%lsc_kick_transverse_dependence) then
      x = p%vec(1)
      y = p%vec(3)
      dpz = 0

      slice => csr%slice(i0)
      if (slice%coef_lsc_plus(0,0) /= 0)  dpz = dpz + r0 / da2_evaluate(slice%coef_lsc_plus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)
      if (slice%coef_lsc_minus(0,0) /= 0) dpz = dpz + r0 / da2_evaluate(slice%coef_lsc_minus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)

      slice => csr%slice(i0+1)
      if (slice%coef_lsc_plus(0,0) /= 0)  dpz = dpz + r1 / da2_evaluate(slice%coef_lsc_plus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)
      if (slice%coef_lsc_minus(0,0) /= 0) dpz = dpz + r1 / da2_evaluate(slice%coef_lsc_minus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)

      p%vec(6) = p%vec(6) + dpz

    else
      p%vec(6) = p%vec(6) + r0 * csr%slice(i0)%kick_lsc + r1 * csr%slice(i0+1)%kick_lsc
    endif

    ! Must update beta and z due to the energy change

    beta0 = p%beta
    call convert_pc_to ((1+p%vec(6))* p%p0c, p%species, beta = p%beta)
    p%vec(5) = p%vec(5) * p%beta / beta0

    ! Transverse space charge. Like the beam-beam interaction but in this case the particles are going
    ! in the same direction longitudinally. In this case, the kicks due to the electric and magnetic fields
    ! tend to cancel instead of add as in the bbi case.

    f0 = csr%kick_factor * csr%actual_track_step * classical_radius(csr%species) / (twopi * &
             csr%dz_slice * csr%rel_mass * e_charge * abs(charge_of(p%species)) * (1 + p%vec(6)) * csr%gamma**3)

    dpx = 0;  dpy = 0

    ! Interpolate between slices. If one of the slices does not have 

    slice => csr%slice(i0)
    call bbi_kick (p%vec(1)-slice%x0, p%vec(3)-slice%y0, [slice%sig_x, slice%sig_y], nk, dnk)
    f = f0 * r0 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    dpx = -nk(1) * f
    dpy = -nk(2) * f

    slice => csr%slice(i0+1)
    call bbi_kick (p%vec(1)-slice%x0, p%vec(3)-slice%y0, [slice%sig_x, slice%sig_y], nk, dnk)
    f = f0 * r1 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    dpx = dpx - nk(1) * f   
    dpy = dpy - nk(2) * f

    if (slice%sig_x == 0 .or. slice%sig_y == 0) call err_exit

    ! 1/2 of the kick is electric and this leads to an energy change.
    ! The formula for p%vec(6) assumes that dpx and dpy are small so only the linear terms are kept.

    p%vec(2) = p%vec(2) + dpx
    p%vec(4) = p%vec(4) + dpy
    p%vec(6) = p%vec(6) + 0.5_rp * (dpx*p%vec(2) + dpy*p%vec(4))
  endif


  enddo
  !$OMP end parallel do
  endif  ! end of general case (else branch)
endif

! Mesh Space charge kick

if (ele%space_charge_method == fft_3d$) then
  call apply_fft_3d_kicks(csr, particle)
endif

end subroutine csr_and_sc_apply_kicks

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine apply_fft_3d_kicks (csr, particle)
!
! Routine to apply FFT-based 3D space charge kicks to particles.
! Deposits charge on a mesh, solves for the field, interpolates back, and applies kicks.
!
! Input:
!   csr         -- csr_struct: Contains mesh, position arrays, and tracking parameters.
!   particle(:) -- coord_struct: Particles to kick.
!
! Output:
!   particle(:) -- coord_struct: Particles with kicks applied.
!-

subroutine apply_fft_3d_kicks (csr, particle)

implicit none

type (csr_struct), target :: csr
type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p

real(rp) :: dct_ave, factor, pz0, ef, dpz, new_beta
integer :: i, n, n_interp
real(rp), allocatable :: x_interp(:), y_interp(:), z_interp(:), E_batch(:,:)
integer, allocatable :: interp_idx(:)

!

if (.not. allocated(csr%position)) allocate(csr%position(size(particle)))
if (size(csr%position) < size(particle)) then
  deallocate(csr%position)
  allocate(csr%position(size(particle)))
endif

! Do the calculation with respect to the average of (time - time_ref) so that adding a constant time offset
! will not affect the calculation.

dct_ave = sum(particle%vec(5)/particle%beta, particle%state == alive$) / count(particle%state == alive$)

n = 0
do i = 1, size(particle)
  p => particle(i)
  if (p%state /= alive$) cycle
  n = n + 1
  csr%position(n)%r = [p%vec(1), p%vec(3), p%vec(5) - dct_ave * p%beta]
  csr%position(n)%charge = p%charge
enddo

call deposit_particles (csr%position(1:n)%r(1), csr%position(1:n)%r(2), csr%position(1:n)%r(3), csr%mesh3d, qa=csr%position(1:n)%charge, &
  mesh_growth_factor=space_charge_com%mesh_growth_factor, mesh_shrink_factor=space_charge_com%mesh_shrink_factor)
call space_charge_3d(csr%mesh3d)

! Gather coordinates for alive particles
allocate(x_interp(size(particle)), y_interp(size(particle)), z_interp(size(particle)))
allocate(interp_idx(size(particle)), E_batch(3, size(particle)))
n_interp = 0
do i = 1, size(particle)
  if (particle(i)%state /= alive$) cycle
  n_interp = n_interp + 1
  x_interp(n_interp) = particle(i)%vec(1)
  y_interp(n_interp) = particle(i)%vec(3)
  z_interp(n_interp) = particle(i)%vec(5) - dct_ave * particle(i)%beta
  interp_idx(n_interp) = i
enddo

! Batch interpolation and kick application
if (n_interp > 0) then
  call interpolate_field_batch(x_interp, y_interp, z_interp, csr%mesh3d, n_interp, E=E_batch)

  !$OMP PARALLEL DO PRIVATE(n, i, p, factor, pz0, ef, dpz, new_beta)
  do n = 1, n_interp
    i = interp_idx(n)
    p => particle(i)
    factor = csr%actual_track_step / (p%p0c * p%beta)
    pz0 = sqrt((1.0_rp + p%vec(6))**2 - p%vec(2)**2 - p%vec(4)**2)
    ! Considering magnetic field also, effectively reduces this force by 1/gamma^2
    p%vec(2) = p%vec(2) + E_batch(1,n) * factor / csr%mesh3d%gamma**2
    p%vec(4) = p%vec(4) + E_batch(2,n) * factor / csr%mesh3d%gamma**2
    ef = E_batch(3,n) * factor
    dpz = sqrt_alpha(1 + p%vec(6), ef*ef + 2 * ef * pz0)
    p%vec(6) = p%vec(6) + dpz
    call convert_pc_to (p%p0c * (1 + p%vec(6)), p%species, beta = new_beta)
    p%vec(5) = p%vec(5) * new_beta / p%beta
    p%beta = new_beta
  enddo
  !$OMP END PARALLEL DO
endif

deallocate(x_interp, y_interp, z_interp, interp_idx, E_batch)

end subroutine apply_fft_3d_kicks

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function dspline_len (s_chord0, s_chord1, spline, dtheta_ref) result (dlen)
!
! Routine to calculate the difference in length between the spline curve length and a referece line.
! Referece line is centroid chord (referece system of the spline) rotated by dtheta_ref.
!
! Input:
!   s_chord0    -- real(rp): Start position along centroid chord.
!   s_chord1    -- real(rp): Stop position along central_chord.
!   spline      -- spline_struct: Spline of x-position as a function of s.
!   dtheta_ref  -- real(rp), optional: angle to rotate the reference line from the centroid chord.
!                    Default is 0.
!
! Output:
!   dlen -- real(rp): L_spline - L_chord
!-

function dspline_len (s_chord0, s_chord1, spline, dtheta_ref) result (dlen)

implicit none

type (spline_struct) spline

real(rp) s_chord0, s_chord1, dlen
real(rp), optional :: dtheta_ref
real(rp) c(0:3), s0, ds

! x' = c(1) + 2*c2*s + 3*c3*s^2
! dlen = Integral: x'^2/2 ds

c = spline%coef
s0 = s_chord0
ds = s_chord1 - s_chord0

if (present(dtheta_ref)) then
  c(1) = c(1) - dtheta_ref
endif

dlen = (ds / 2) * ( &
        c(1)**2 + &
        (2*s0 + ds) * 2*c(1)*c(2) + &
        (3*s0*s0 + 3*s0*ds + ds*ds) * (6*c(1)*c(3) + 4*c(2)**2) / 3 + &
        (4*s0**3 + 6*s0*s0*ds + 4*s0*ds*ds + ds**3) * 3*c(2)*c(3) + &
        (5*s0**4 + 10*s0**3*ds + 10*s0*s0*ds*ds + 5*s0*ds**3 + ds**4) * 9*c(3)**2 &
       )

end function dspline_len

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function s_ref_to_s_chord (s_ref, eleinfo) result (s_chord)
!
! Routine to calculate s_chord given s_ref.
!
! Input:
!   s_ref     -- real(rp): s-position along element ref coords.
!   eleinfo     -- csr_ele_info_struct: Element info
!
! Output:
!   s_chord   -- real(rp): s-position along centroid chord.
!-

function s_ref_to_s_chord (s_ref, eleinfo) result (s_chord)

implicit none

type (csr_ele_info_struct), target :: eleinfo
type (ele_struct), pointer :: ele

real(rp) s_ref, s_chord, dtheta, dr(3), x, g, t

!

ele => eleinfo%ele
g = ele%value(g$)

if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. abs(g) > 1d-5) then
  dtheta = eleinfo%floor0%theta - eleinfo%ref_floor0%theta
  dr = eleinfo%ref_floor0%r - eleinfo%floor0%r
  t = eleinfo%ref_floor0%theta + pi/2
  x = dr(1) * cos(t) + dr(3) * sin(t)
  s_chord = abs(ele%value(rho$)) * atan2(s_ref * cos(dtheta), abs(ele%value(rho$) + x + s_ref * sin(dtheta)))

else
  s_chord = s_ref * eleinfo%L_chord /ele%value(l$)
endif

end function s_ref_to_s_chord







!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr3d (bunch, ele, centroid, err, bunch_track)
!
! EXPERIMENTAL. NOT CURRENTLY OPERATIONAL!
!
! Routine to track a bunch of particles through an element using
! steady-state 3D CSR.
!
!
! Input:
!   bunch         -- Bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: The element to track through. Must be part of a lattice.
!   centroid(0:)  -- coord_struct, Approximate beam centroid orbit for the lattice branch.
!                      Calculate this before beam tracking by tracking a single particle.
!   s_start       -- real(rp), optional: Starting position relative to ele. Default = 0
!   s_end         -- real(rp), optional: Ending position. Default is ele length.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- Bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. EG: Too many particles lost.
!   bunch_track -- bunch_track_struct, optional: track information if the tracking method does
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!
!
! Notes:
!   The core routines are from the OpenCSR package developed at:
!   https://github.com/ChristopherMayes/OpenCSR
!
!-

subroutine track1_bunch_csr3d (bunch, ele, centroid, err, s_start, s_end, bunch_track)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: c0, p
type (ele_struct), target :: ele
type (bunch_track_struct), optional :: bunch_track
type (branch_struct), pointer :: branch
type (ele_struct) :: runt
type (csr_struct), target :: csr
type (coord_struct), target :: centroid(0:)
type (coord_struct), pointer :: particle(:)

real(rp), optional :: s_start, s_end
real(rp) s0_step
real(rp) e_tot, ds_step
real(rp) Evec(3), factor, pz0
integer i, j, n, ie, ns, nb, n_step, n_live, i_step

character(*), parameter :: r_name = 'track1_bunch_csr3d'
logical err, err_flag, parallel0, parallel1

! EXPERIMENTAL. NOT CURRENTLY OPERATIONAL!

! Init

err = .true.
branch => pointer_to_branch(ele)

! This calc is only valid in SBEND elements, nonzero length
! No CSR for a zero length element.
! And taylor elements get ignored.
if (ele%value(l$) == 0 .or. ele%key /= sbend$) then
  call track1_bunch_hom (bunch, ele)
  err = .false.
  return
endif

! Set gamma, mesh size
c0 => centroid(ele%ix_ele)
call convert_pc_to((1+c0%vec(6)) * c0%p0c, c0%species, gamma = csr%mesh3d%gamma)
csr%mesh3d%nhi = space_charge_com%csr3d_mesh_size

! n_step is the number of steps to take when tracking through the element.
! csr%ds_step is the true step length.

particle => bunch%particle
! make sure that ele_len / track_step is an integer.

csr%ds_track_step = ele%value(csr_ds_step$)
if (csr%ds_track_step == 0) csr%ds_track_step = space_charge_com%ds_track_step
if (csr%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'NEITHER SPACE_CHARGE_COM%DS_TRACK_STEP NOR CSR_TRACK_STEP FOR THIS ELEMENT ARE SET! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

n_step = max (1, nint(ele%value(l$) / csr%ds_track_step))
csr%ds_track_step = ele%value(l$) / n_step
csr%species = bunch%particle(1)%species

if (.not. allocated(csr%position)) allocate(csr%position(size(particle)))
if (size(csr%position) < size(particle)) then
  deallocate(csr%position)
  allocate(csr%position(size(particle)))
endif

!----------------------------------------------------------------------------------------
! Loop over the tracking steps
! runt is the element that is tracked through at each step.

call save_a_bunch_step (ele, bunch, bunch_track, s_start)

do i_step = 0, n_step

  ! track through the runt

  if (i_step /= 0) then
    call element_slice_iterator (ele, branch%param, i_step, n_step, runt, s_start, s_end)
    call track1_bunch_hom (bunch, runt)
  endif
  particle => bunch%particle
  
  s0_step = i_step * csr%ds_track_step
  if (present(s_start)) s0_step = s0_step + s_start

  e_tot = ele%value(e_tot$)
  call convert_total_energy_to (e_tot, branch%param%particle, csr%gamma, beta = csr%beta)
  csr%gamma2 = csr%gamma**2

  ! Bin particles
  ! TODO: simplify this into a function, and use the function above for fft_3d
  n = 0
  do i = 1, size(particle)
    p => particle(i)
    if (p%state /= alive$) cycle
    n = n + 1
    csr%position(n)%r = p%vec(1:5:2)
    csr%position(n)%charge = p%charge
  enddo
  call deposit_particles (csr%position(1:n)%r(1), csr%position(1:n)%r(2), csr%position(1:n)%r(3), csr%mesh3d, qa=csr%position(1:n)%charge, &
    mesh_growth_factor=space_charge_com%mesh_growth_factor, mesh_shrink_factor=space_charge_com%mesh_shrink_factor)  

  ! Give particles a kick
  ! TODO: simplify with fft_3d
  print *, '---------- CSR Steady_State_3D ----------'
  
  call csr3d_steady_state_solver(csr%mesh3d%rho/product(csr%mesh3d%delta), &
              csr%mesh3d%gamma, ele%value(rho$), csr%mesh3d%delta, csr%mesh3d%efield, normalize=.false.)
  
  ! Handle first and last steps for half-kicks
  csr%kick_factor = csr%ds_track_step
  if (i_step == 0 .or. i_step == n_step) csr%kick_factor = csr%kick_factor / 2
  print *, 'kick by ds', csr%kick_factor, i_step
  
  
  print *, 'Interpolating field and kicking particles'
  do i = 1, size(particle)
    p => particle(i)
    if (p%state /= alive$) cycle
    call interpolate_field(p%vec(1), p%vec(3), p%vec(5), csr%mesh3d, E=Evec)

    factor = csr%kick_factor / (p%p0c  * p%beta) 

    pz0 = sqrt( (1.0_rp + p%vec(6))**2 - p%vec(2)**2 - p%vec(4)**2 ) ! * p0 
    p%vec(2) = p%vec(2) + Evec(1)*factor
    p%vec(4) = p%vec(4) + Evec(2)*factor
    p%vec(6) = sqrt(p%vec(2)**2 + p%vec(4)**2 + (Evec(3)*factor + pz0)**2) -1.0_rp
    ! Set beta
    call convert_pc_to (p%p0c * (1 + p%vec(6)), p%species, beta = p%beta)
  enddo

  call save_a_bunch_step (ele, bunch, bunch_track, s0_step+s_start)

enddo

err = .false.

end subroutine track1_bunch_csr3d

end module
