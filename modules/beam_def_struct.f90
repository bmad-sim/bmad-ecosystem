module beam_def_struct

use bmad_struct
use bmad_interface

! Sigma matrix elements (21 total)
integer, parameter :: s11$ = 1, s12$ = 2, s13$ = 3, s14$ =  4, s15$ =  5
integer, parameter :: s16$ = 6, s22$ = 7, s23$ = 8, s24$ = 9
integer, parameter :: s25$ = 10, s26$ = 11, s33$ = 12, s34$ = 13, s35$ = 14
integer, parameter :: s36$ = 15, s44$ = 16, s45$ = 17, s46$ = 18
integer, parameter :: s55$ = 19, s56$ = 20, s66$ = 21

integer, parameter :: not_lost$ = -1

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
  real(rp) charge   ! total charge in a bunch (Coul).
  real(rp) s_center ! longitudinal center of bunch (m).
end type

type beam_struct
  type (bunch_struct), allocatable :: bunch(:)
end type

type beam_init_struct
  real(rp) a_norm_emitt     ! a-mode emittance
  real(rp) b_norm_emitt     ! b-mode emittance
  real(rp) :: dPz_dz = 0    ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0 ! Bench center offset relative to reference.
  real(rp) ds_bunch         ! Distance between bunches.
  real(rp) sig_z            ! Z sigma in m.
  real(rp) sig_e            ! e_sigma in dE/E.
  real(rp) bunch_charge     ! charge in a bunch.
  real(rp) :: center_jitter(6) = 0.0 ! Bunch center rms jitter
  real(rp) :: emitt_jitter(2)  = 0.0 ! a and b bunch emittance rms jitter normalized to emittance
  real(rp) :: sig_z_jitter     = 0.0 ! bunch length RMS jitter 
  real(rp) :: sig_e_jitter     = 0.0 ! energy spread RMS jitter 
  type(beam_spin_struct)  spin       ! Initialize the spin
  integer :: n_particle = 0          ! Number of simulated particles per bunch.
  integer :: n_bunch = 0             ! Number of bunches.
  logical :: renorm_center = .true.  ! Renormalize centroid?
  logical :: renorm_sigma = .true.   ! Renormalize sigma?
  logical :: preserve_dist = .false. ! use the same grid distributon each time
  logical :: init_spin     = .false. ! initialize beam spinors
end type

type bunch_param_struct
  real(rp) beta, alpha, gamma
  real(rp) eta, etap
  real(rp) norm_emitt ! normalized emittance
end type

type bunch_params_struct
  type (bunch_param_struct) :: x, y, z, a, b, c
  type (coord_struct) :: centroid  ! Lab frame
  type (beam_spin_struct) :: spin  ! polarization
  real(rp) sigma(21)               ! projected sigma matrix
  real(rp) sigma_normal(21)        ! normal mode sigma matrix
  real(rp) s                       ! Longitudinal position.
  integer n_particle               ! all non-lost particles
end type

! How close to polarization vector for particle to be polarized?	

 real(rp), parameter, private ::  sigma_theta = 1e-3 ! 1 milliradian
 real(rp), parameter, private ::  sigma_phi = 1e-3

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_particle (start, ele, param, end)
!
! Subroutine to track a particle through an element.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   start  -- struct: Starting coords.
!   ele    -- Ele_struct: Element to track through.
!   param  -- Param_struct: Global parameters.
!
! Output:
!   end    -- struct: Ending coords.
!-

subroutine track1_particle (start, ele, param, end)

  implicit none

  type (particle_struct) :: start
  type (particle_struct) :: end
  type (ele_struct) :: ele
  type (param_struct), intent(inout) :: param

! transfer z-order index, charge, etc

  end = start
  if (start%ix_lost /= not_lost$) return
  if (ele%key == marker$) return

  call track1 (start%r, ele, param, end%r)
  if (param%lost) end%ix_lost = ele%ix_ele

  if (end%ix_lost /= not_lost$) then
    end%r%vec = 0
    end%charge = 0
    return
  endif

end subroutine

end module
