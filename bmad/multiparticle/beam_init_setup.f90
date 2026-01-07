!+
! Function beam_init_setup (beam_init_in, ele, species, modes, err_flag) result (beam_init_set)
!
! Routine to setup a beam_init_struct instance:
!   * Error checks
!   * Comput norm_emit from emit or vice versa depending upon what is set.
!   * If modes is present: Merge modes info into beam_init_in if parameters in beam_init_in are set.
!     For example, if beam_init_in%sig_z is set negative, beam_init_set%sig_z will be set to the value in modes.
!
! Input:
!   beam_init_in  -- beam_init_struct: Input parameters
!   ele           -- ele_struct:
!   species       -- integer: Beam particle species.
!   modes         -- normal_modes_struct, optional: Normal mode parameters.
!
! Ouput:
!   err_flag      -- logical, optional: Set true if there is an error. False otherwise.
!   beam_init_set -- beam_init_struct: See above.
!-

function beam_init_setup (beam_init_in, ele, species, modes, err_flag) result (beam_init_set)

use bmad_routine_interface, dummy => beam_init_setup

implicit none

type (beam_init_struct), target :: beam_init_in, beam_init_set
type (beam_init_struct), pointer :: bis
type (ele_struct) ele
type (normal_modes_struct), optional :: modes

real(rp) ran_g(2), beta_gamma
integer species
logical, optional :: err_flag

character(*), parameter :: r_name = 'beam_init_setup'

! Check

if (present(err_flag)) err_flag = .true.
bis => beam_init_set
bis = beam_init_in
if (ele%value(p0c$) == 0) return

! Sanity checks.

if (bis%a_emit /= 0 .and. bis%a_norm_emit /= 0 .and. .not. (bis%a_emit < 0 .and. bis%a_norm_emit < 0)) then
  call out_io (s_error$, r_name, 'I am confused! Both a_emit and a_norm_emit are set non-zero in the beam_init_struct structure.', &
                                   'Please set one or the other to zero.')
  return
endif

if (bis%b_emit /= 0 .and. bis%b_norm_emit /= 0 .and. .not. (bis%b_emit < 0 .and. bis%b_norm_emit < 0)) then
  call out_io (s_error$, r_name, 'I am confused! Both b_emit and b_norm_emit are set non-zero in the beam_init_struct structure.', &
                                   'Please set one or the other to zero.')
  return
endif

! a-emit

beta_gamma = ele%value(p0c$) / mass_of(species)

if (bis%a_emit < 0 .or. bis%a_norm_emit < 0) then
  if (present(modes)) then
    bis%a_emit = modes%a%emittance
    bis%a_norm_emit = modes%a%emittance * beta_gamma
  else
    call out_io (s_warn$, r_name, &
                    'a_emit and/or a_norm_emit is set negative but the calling program has not provided lattice emittance numbers!', &
                    'Will use a value of zero...')
    bis%a_emit = 0
    bis%a_norm_emit = 0
  endif
elseif (bis%a_norm_emit /= 0) then
  bis%a_emit = bis%a_norm_emit / beta_gamma
else
  bis%a_norm_emit = bis%a_emit * beta_gamma
endif

! b-emit

if (bis%b_emit < 0 .or. bis%b_norm_emit < 0) then
  if (present(modes)) then
    bis%b_emit = modes%b%emittance
    bis%b_norm_emit = modes%b%emittance * beta_gamma
  else
    call out_io (s_warn$, r_name, &
                    'b_emit and/or b_norm_emit is set negative but the calling program has not provided lattice emittance numbers!', &
                    'Will use a value of zero...')
    bis%b_emit = 0
    bis%b_norm_emit = 0
  endif
elseif (bis%b_norm_emit /= 0) then
  bis%b_emit = bis%b_norm_emit / beta_gamma
else
  bis%b_norm_emit = bis%b_emit * beta_gamma
endif

! sig_z

if (bis%sig_z < 0) then
  if (present(modes)) then
    bis%sig_z = modes%sig_z
  else
    call out_io (s_warn$, r_name, &
                    'sig_z is set negative but the calling program has not provided lattice emittance numbers!', &
                    'Will use a value of zero...')
    bis%sig_z = 0
  endif
endif

! sig_pz

if (bis%sig_pz < 0) then
  if (present(modes)) then
    bis%sig_pz = modes%sigE_E
  else
    call out_io (s_warn$, r_name, &
                    'sig_pz is set negative but the calling program has not provided lattice emittance numbers!', &
                    'Will use a value of zero...')
    bis%sig_pz = 0
  endif
endif

! Checking that |bis%dpz_dz| < mode%sigE_E / mode%sig_z

if (abs(bis%dPz_dz * bis%sig_z) > bis%sig_pz) then
  call out_io (s_error$, r_name, "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  return
endif

! Add jitter if needed

if (any(bis%emit_jitter /= 0)) then
  call ran_gauss(ran_g) ! ran(3:4) for z and e jitter used below
  bis%a_emit      = bis%a_emit      * (1 + bis%emit_jitter(1) * ran_g(1))
  bis%b_emit      = bis%b_emit      * (1 + bis%emit_jitter(2) * ran_g(2))
  bis%a_norm_emit = bis%a_norm_emit * (1 + bis%emit_jitter(1) * ran_g(1))
  bis%b_norm_emit = bis%b_norm_emit * (1 + bis%emit_jitter(2) * ran_g(2))
endif

if (present(err_flag)) err_flag = .false.

end function beam_init_setup 

