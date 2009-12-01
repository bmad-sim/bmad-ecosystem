module ping_anal_mod

use bmad

! Structure for a given bpm

type bpm_ping_struct
  character(8) name       ! bpm name
  integer index           ! bpm index
  real(rp), allocatable :: x(:), y(:)  ! orbit data
  ! Fit parameters...
  real(rp) x0, y0                   ! Baseline
  real(rp) amp_a, amp_b, amp_z      ! Osc Amp
  real(rp) beta_a, beta_b           ! Beta functions
  real(rp) phi_a, phi_b             ! Phase advance
  real(rp) Cbar11, Cbar12, Cbar22   ! Coupling
  real(rp) eta_x, eta_y             ! Dispersions
end type

! Input parameters

type ping_param_struct
  character(100) lattice_file
  character(100) data_file
  integer ix_phase_meas       ! Phase measurement index.
  integer turn_start          ! Start turn for analysis window.
  integer turn_stop           ! Stop turn for analysis window.
end type

! Structure for everything

type ping_universe_struct
  type (ping_param_struct) param
  type (bpm_ping_struct), allocatable :: bpm(:)
  ! Fit parameters...
  real(rp) n_damp     ! Damping number of turns.
  real(rp) n_damp_z   ! Longitudinal damping.
  real(rp) mu_oct     ! Octupole decoherence.
  real(rp) sig_x      ! Horizontal sigma (needed for octupole decoherence)
  real(rp) sig_y      ! Vertical sigma (needed for octupole decoherence)
  real(rp) Q_a, Q_b   ! Tunes.
end type


end module
