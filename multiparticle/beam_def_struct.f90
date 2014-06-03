module beam_def_struct

use bmad_struct

type beam_spin_struct
  real(rp) :: polarization = 1.0 ! 100% polarized is default.
  real(rp) :: theta = 0.0        ! polarization direction in polar coordinates
  real(rp) :: phi = 0.0          ! polarization direction in polar coordinates
end type

type bunch_struct
  type (coord_struct), allocatable :: particle(:)
  integer, allocatable :: ix_z(:)  ! bunch%ix_z(1) is index of head particle, etc.
  real(rp) charge_tot   ! Total charge in a bunch (Coul).
  real(rp) charge_live  ! Charge of live particles (Coul).
  real(rp) z_center     ! Longitudinal center of bunch (m). Note: Generally, z_center of 
                        !   bunch #1 is 0 and z_center of the other bunches is negative.
  real(rp) t_center     ! Center of bunch creation time relative to head bunch.
  integer ix_ele        ! Element this bunch is at.
  integer ix_bunch      ! Bunch index. Head bunch = 1, etc.
end type

type beam_struct
  type (bunch_struct), allocatable :: bunch(:)
end type

type ellipse_beam_init_struct
  integer :: part_per_ellipse = 0  ! number of particles per ellipse
  integer :: n_ellipse = 1         ! number of ellipses (>= 1)
  real(rp) :: sigma_cutoff = 0     ! sigma cutoff of the representation
end type

type kv_beam_init_struct
  integer :: part_per_phi(2) = 0    ! number of particles per angle variable.
  integer :: n_I2 = 0               ! number of I2
  real(rp) :: A = 0                 ! A = I1/e
end type

type grid_beam_init_struct
  integer :: n_x = 0         ! Number of columns.
  integer :: n_px = 0        ! Number of rows.
  real(rp) :: x_min = 0      ! Lower x limit.
  real(rp) :: x_max = 0      ! Upper x limit.
  real(rp) :: px_min = 0     ! Lower px limit.
  real(rp) :: px_max = 0     ! Upper px limit,
end type

type beam_init_struct
  character(200) :: file_name = ''           ! For distribution_type(1) = 'FILE'
  character(16) :: distribution_type(3) = '' ! distribution type (in x-px, y-py, and z-pz planes)
                                             ! "ELLIPSE", "KV", "GRID", "FILE", "", or "RAN_GAUSS" 
  type (ellipse_beam_init_struct) ellipse(3) ! Parameters for ellipse beam distribution
  type (kv_beam_init_struct) KV              ! Parameters for KV beam distribution
  type (grid_beam_init_struct) grid(3)       ! Parameters for grid beam distribution
  !!! The following are for Random distributions
  real(rp) :: center_jitter(6) = 0.0  ! Bunch center rms jitter
  real(rp) :: emit_jitter(2)   = 0.0  ! a and b bunch emittance rms jitter normalized to emittance
  real(rp) :: sig_z_jitter     = 0.0  ! bunch length RMS jitter 
  real(rp) :: sig_e_jitter     = 0.0  ! energy spread RMS jitter 
  integer :: n_particle = 0           ! Number of random particles per bunch.
  logical :: renorm_center = .true.   ! Renormalize centroid?
  logical :: renorm_sigma = .true.    ! Renormalize sigma?
  character(16) :: random_engine = 'pseudo' ! Or 'quasi'. Random number engine to use. 
  character(16) :: random_gauss_converter = 'exact'  
                                            ! Or 'quick'. Uniform to gauss conversion method.
  real(rp) :: random_sigma_cutoff = -1      ! Cut-off in sigmas.
  !!! The following are used  by all distribution types
  type(beam_spin_struct)  spin        ! Initialize the spin
  real(rp) :: a_norm_emit = 0         ! a-mode normalized emittance (emit * gamma)
  real(rp) :: b_norm_emit = 0         ! b-mode normalized emittance (emit * gamma)
  real(rp) :: a_emit = 0              ! a-mode emittance
  real(rp) :: b_emit = 0              ! b-mode emittance
  real(rp) :: dPz_dz = 0              ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0           ! Bench center offset relative to reference.
  real(rp) :: dt_bunch = 0            ! Time between bunches.
  real(rp) :: sig_z = 0               ! Z sigma in m.
  real(rp) :: sig_e = 0               ! e_sigma in dE/E.
  real(rp) :: bunch_charge = 1        ! charge in a bunch. Default is something non-zero.
  integer :: n_bunch = 1              ! Number of bunches.
  integer :: species = not_set$       ! positron$, etc. not_set$ => use referece particle.
  logical :: init_spin     = .false.  ! initialize beam spinors
  logical :: full_6D_coupling_calc = .false. ! Use V from 6x6 1-turn mat to match distribution?  
                                             !   Else use 4x4 1-turn mat used.
  logical :: use_t_coords = .false.   ! If true, the distributions will be taken as in t-coordinates  
  logical :: use_z_as_t   = .false.   ! Only used if  use_t_coords = .true.
                                      !   If true,  z describes the t distribution 
                                      !   If false, z describes the s distribution
end type

type bunch_params_struct
  type (twiss_struct) :: x, y, z, a, b, c
  type (coord_struct) :: centroid  ! Lab frame
  type (beam_spin_struct) :: spin  ! polarization
  real(rp) sigma(6,6)              ! beam size matrix
  real(rp) s                       ! Longitudinal position.
  real(rp) charge_live             ! Charge of all non-lost particle
  integer n_particle_tot
  integer n_particle_live          ! all non-lost particles
  integer n_particle_lost_in_ele
end type

!-------------------------------------------------------------------------
! CSR parameters

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

!+
! Note: Shielding is simulated via the image current due to the 
! top and bottom walls. The side walls are neglected.
!-

type csr_parameter_struct                  ! Common block for csr calc
  real(rp) :: ds_track_step = 0            ! Tracking step size
  real(rp) :: beam_chamber_height = 0      ! Used in shielding calculation.
  real(rp) :: sigma_cutoff = 0.1           ! Cutoff for the lsc calc. If a bin sigma
                                           !  is < cutoff * sigma_ave then ignore.
  integer :: n_bin = 0                     ! Number of bins used
  integer :: particle_bin_span = 2         ! Longitudinal particle length / dz_bin
  integer :: n_shield_images = 0           ! Chamber wall shielding. 0 = no shielding.
  integer :: ix1_ele_csr = -1              ! Start index for csr tracking
  integer :: ix2_ele_csr = -1              ! Stop index for csr tracking
  logical :: lcsr_component_on = .true.    ! Longitudinal csr component
  logical :: lsc_component_on = .true.     ! Longitudinal space charge component
  logical :: tsc_component_on = .false.    ! Transverse space charge component
  logical :: small_angle_approx = .true.   ! Use lcsr small angle approximation?
  logical :: print_taylor_warning = .true. ! Print warning if Taylor element is present?
end type

type (csr_parameter_struct), save, target :: csr_param

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If n_bunch = 0 then all macro beam pointers will be deallocated.
! Rule: If beam%bunch(:) is allocated, beam%bunch(i)%particle(:) will be allocated.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   beam -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle)

implicit none

type (beam_struct) beam

integer i, n_bunch, n_particle

! Deallocate if needed

if (allocated(beam%bunch)) then
  if (n_bunch == 0 .or. size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

do i = 1, n_bunch
  call reallocate_bunch (beam%bunch(i), n_particle)
enddo

end subroutine reallocate_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reallocate_bunch (bunch, n_particle)
! 
! Subroutine to reallocate particles within a bunch_struct.
!
! Modules needed:
!   use bunch_mod
!
! Input:
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   bunch -- bunch_struct: Allocated bunch_struct structure.
!-

subroutine reallocate_bunch (bunch, n_particle)

implicit none

type (bunch_struct) bunch

integer i, n_particle

! Deallocate if needed

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= n_particle) deallocate (bunch%particle, bunch%ix_z)
endif

if (.not. allocated(bunch%particle)) then
  allocate (bunch%particle(n_particle), bunch%ix_z(n_particle))
  bunch%ix_z = 0
endif

end subroutine reallocate_bunch

end module
