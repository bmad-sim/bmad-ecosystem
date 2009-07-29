module beam_def_struct

use bmad_struct

! Sigma matrix elements (21 total)
integer, parameter :: s11$ = 1, s12$ = 2, s13$ = 3, s14$ =  4, s15$ =  5
integer, parameter :: s16$ = 6, s22$ = 7, s23$ = 8, s24$ = 9
integer, parameter :: s25$ = 10, s26$ = 11, s33$ = 12, s34$ = 13, s35$ = 14
integer, parameter :: s36$ = 15, s44$ = 16, s45$ = 17, s46$ = 18
integer, parameter :: s55$ = 19, s56$ = 20, s66$ = 21

! Phase space distribution types
integer, parameter :: grid$ = 1       ! uniform rectangular grid
integer, parameter :: ellipse$ = 2    ! ellipses representing a Gaussian dist.
integer, parameter :: KV$ = 3         ! Kapchinsky-Vladimirsky distribution

type beam_spin_struct
  real(rp) :: polarization = 1.0 ! i.e. 80% polarized
  real(rp) :: theta = 0.0  ! polar coordinates
  real(rp) :: phi = 0.0    ! polar coordinates
end type

type particle_struct
  type (coord_struct) r   ! Center of the particle
  real(rp) charge         ! charge in a particle (Coul).
  integer :: ix_z = 0     ! Index for ordering the particles longitudinally.
                          !   particle(1)%ix_z is index of head particle.
  integer :: ix_lost = not_lost$  ! When the particle been lost in tracking
                                  !   ix_lost set to index of element where lost.
end type

type bunch_struct
  type (particle_struct), allocatable :: particle(:)
  real(rp) charge   ! Total charge in a bunch (Coul).
  real(rp) z_center ! Longitudinal center of bunch (m). Note: Generally, z_center of 
                    !   bunch #1 is 0 and z_center of the other bunches is negative.
  real(rp) t_center ! Longitudinal center of bunch (sec).
  integer ix_ele    ! Element this bunch is at.
  integer ix_bunch  ! Bunch index. Head bunch = 1, etc.
end type

type beam_struct
  type (bunch_struct), allocatable :: bunch(:)
end type

type tail_weighted_beam_init_struct
   integer :: type(3) = 0              ! distribution type (in x-px, y-py, and z-pz planes)
   
                                       ! Params for ellipse$ and KV$
   integer :: part_per_ellipse(3) = 0    ! number of particles per ellipse
   
                                       ! Params for ellipse$
   integer :: n_ellipse(3) = 1           ! number of ellipses (>= 1)
   real(rp) :: sigma_cutoff(3) = 0       ! sigma cutoff of the representation

                                       ! Params for KV$
   integer :: n_I2 = 0                   ! number of I2
   real(rp) :: A = 0                     ! A = I1/e

                                       ! Params for grid$
   integer :: n_x(3) = 0                 ! number of columns
   integer :: n_px(3) = 0                ! number of rows
   real(rp) :: minima(6) = 0             ! upper and lower limits in (x,px,y,py,z,pz)
   real(rp) :: maxima(6) = 0
end type

type beam_init_struct
  real(rp) a_norm_emitt     ! a-mode emittance
  real(rp) b_norm_emitt     ! b-mode emittance
  real(rp) :: dPz_dz = 0    ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0 ! Bench center offset relative to reference.
  real(rp) dt_bunch         ! Time between bunches.
  real(rp) sig_z            ! Z sigma in m.
  real(rp) sig_e            ! e_sigma in dE/E.
  real(rp) bunch_charge     ! charge in a bunch.
  real(rp) :: center_jitter(6) = 0.0 ! Bunch center rms jitter
  real(rp) :: emitt_jitter(2)  = 0.0 ! a and b bunch emittance rms jitter normalized to emittance
  real(rp) :: sig_z_jitter     = 0.0 ! bunch length RMS jitter 
  real(rp) :: sig_e_jitter     = 0.0 ! energy spread RMS jitter 
  type(beam_spin_struct)  spin       ! Initialize the spin
  integer :: n_particle = 0          ! Number of simulated particles per bunch.
  integer :: n_bunch = 1             ! Number of bunches.
  logical :: renorm_center = .true.  ! Renormalize centroid?
  logical :: renorm_sigma = .true.   ! Renormalize sigma?
  logical :: init_spin     = .false. ! initialize beam spinors
  logical :: is_random = .true.      ! Random distribution?  Or a tail-weighted distribution?
  character(16) :: random_engine = 'pseudo' ! Or 'quasi'. Random number engine to use. 
  character(16) :: random_gauss_converter = 'exact'  
                                            ! Or 'limited'. Uniform to gauss conversion method.
  real(rp) :: random_sigma_cutoff = 4.0     ! Used with 'limited' converter. Cut-off in sigmas.

  type (tail_weighted_beam_init_struct) tw_beam_init ! Parameters for a tail-weighted beam distribution
end type

type bunch_lat_param_struct
   real(rp) beta, alpha, gamma
   real(rp) eta, etap
   real(rp) norm_emitt ! normalized emittance
end type

type bunch_params_struct
  type (bunch_lat_param_struct) :: x, y, z, a, b, c
  type (coord_struct) :: centroid  ! Lab frame
  type (beam_spin_struct) :: spin  ! polarization
  real(rp) sigma(21)               ! projected sigma matrix
  real(rp) s                       ! Longitudinal position.
  real(rp) charge_live             ! Charge of all non-lost particle
  integer n_live_particle          ! all non-lost particles
end type

! How close to polarization vector for particle to be polarized?	

 real(rp), parameter, private ::  sigma_theta = 1e-3 ! 1 milliradian
 real(rp), parameter, private ::  sigma_phi = 1e-3

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
