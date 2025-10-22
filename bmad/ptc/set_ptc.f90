!+
! Subroutine set_ptc (e_tot, particle, taylor_order, integ_order, n_step, &
!                                             no_cavity, force_init)
!
! Subroutine to initialize PTC.
!
! Note: At some point before you use PTC to compute Taylor maps etc.
!   you have to call set_ptc with both e_tot and particle args present. 
!   Always supply both of these args together or not at all. 
!
! Note: If you just want to use FPP without PTC then call the FPP routine init directly.
!
! Note: This subroutine cannot be used if you want to have "knobs" (in the PTC sense).
!
! Note: Call this routine to transfer the value of the electric dipole moment from 
!   bmad_com%electric_dipole_moment to PTC.
!
! Input:
!   e_tot          -- Real(rp), optional: Energy in eV.
!   particle       -- Integer, optional: Type of particle:
!                       electron$, proton$, etc.
!   taylor_order   -- Integer, optional: Maximum order of the taylor polynomials.
!                       0 => Use default.
!   integ_order    -- Integer, optional: Default Order for the drift-kick-drift 
!                       sympletic integrator. Possibilities are: 2, 4, or 6
!                       Default = 2
!   n_step         -- Integer, optional: Default Number of integration steps.
!                       Default = 1
!   no_cavity      -- Logical, optional: No RF Cavity exists? 
!                       Default = False.
!                       Corresponds to the nocavity option of the PTC init routine.
!                       no_cavity = .true. will turn any cavity into a drift.
!   force_init     -- logical, optional: If present and True then force a PTC init.
!-

subroutine set_ptc (e_tot, particle, taylor_order, integ_order, n_step, no_cavity, force_init) 

use ptc_interface_mod, dummy => set_ptc
use mad_like, only: make_states, pmaMUON, pmaE, PHASE0, &
              assignment(=), nocavity0, operator(+), in_bmad_units, &
              berz, init, set_madx, lp, superkill, TIME0, init_all, SPIN0, C_WATCH_USER
use madx_ptc_module, only: ptc_ini_no_append, append_empty_layout, m_u, bmadl, use_info, &
              use_info_m, check_longitudinal, bmad_automatic, OLD_SURVEY
use c_tpsa, only: c_verbose, E_MUON, USE_QUATERNION

implicit none

integer, optional :: integ_order, particle, n_step, taylor_order
integer this_method, this_steps, t_order

real(rp), optional :: e_tot
real(dp) this_energy

logical, optional :: no_cavity, force_init
logical params_present, c_verbose_save

character(16) :: r_name = 'set_ptc'

! ptc cannot be used with photons

if (logic_option(.false., force_init)) ptc_private%init_ptc_needed = .true.

if (present(particle)) then
  if (particle == photon$) return
endif

! Some init

OLD_SURVEY = .false.
!!BMAD_AUTOMATIC = .true. ! For c_normal calc. This enables automatic testing of whether RF is on or off.
USE_QUATERNION = .true.
E_MUON = bmad_com%electric_dipole_moment
CHECK_LONGITUDINAL = .false. ! MAD-X uses the True setting.
call in_bmad_units

! do not call set_mad

params_present = present(e_tot) .and. present(particle)

if (ptc_private%init_ptc_needed .and. params_present) then
  if (particle == muon$ .or. particle == antimuon$) then
    call make_states (pmaMUON/pmaE)
  elseif (particle == positron$ .or. particle == electron$) then
    call make_states(.true._lp)
  elseif (particle == proton$ .or. particle == antiproton$) then
    if (particle /= proton$ .and. particle /= antiproton$) then
      call out_io (s_error$, r_name, 'PTC IS NOT ABLE TO HANDLE PARTICLES OF TYPE: ' // species_name(particle), 'USING PROTON/ANTIPROTON')
    endif
    call make_states(.false._lp)
  else
    call make_states (mass_of(particle)/mass_of(electron$), anomalous_moment_of(particle), real(charge_of(particle), dp))
    call out_io (s_warn$, r_name, 'Note: Radiation calculation in PTC not correct for particles of type: ' // species_name(particle))
  endif

  ! Use PTC time tracking
  ptc_private%base_state = ptc_private%base_state + TIME0
  PHASE0 = 0
endif

if (present(no_cavity))           ptc_private%base_state = ptc_private%base_state + NOCAVITY0
if (bmad_com%spin_tracking_on)    ptc_private%base_state = ptc_private%base_state + SPIN0

if (present (integ_order)) then
  this_method = integ_order
  bmad_com%default_integ_order = integ_order
else
  this_method = bmad_com%default_integ_order
endif

if (present (n_step)) then
  this_steps = n_step
else
  this_steps = 10
endif

if (present(taylor_order)) then
  t_order = taylor_order
  if (t_order == 0) t_order = ptc_private%taylor_order_saved
  ptc_private%taylor_order_saved = t_order
endif

if (params_present) then
  if (ptc_private%init_ptc_needed .or. ptc_private%e_tot_set /= e_tot .or. present(integ_order) .or. present(n_step)) then
    this_energy = 1d-9 * e_tot
    if (this_energy == 0) then
      call out_io (s_fatal$, r_name, 'E_TOT IS 0.')
      if (global_com%exit_on_error) call err_exit
    endif
    call set_madx (energy = this_energy, method = this_method, step = this_steps)
    ptc_private%e_tot_set  = e_tot
    ! Only do this once
    if (ptc_private%init_ptc_needed) call ptc_ini_no_append 
    ptc_private%init_ptc_needed = .false.
  endif
endif

! Do not call init before the call to make_states.

if (.not. ptc_private%init_ptc_needed .or. logic_option(.false., force_init)) then  ! If make_states has been called
  t_order = 0
  if (present(taylor_order)) t_order = taylor_order
  if (t_order == 0) t_order = bmad_com%taylor_order
  if (t_order == 0) t_order = ptc_private%taylor_order_saved
  if (ptc_private%taylor_order_ptc /= t_order) then
    ! Due to Bmad vs PTC units bug, call init with nocavity
    call init (ptc_private%base_state+NOCAVITY0, t_order, 0)
    ptc_private%init_spin_needed = .true.
    c_verbose_save = c_verbose
    c_verbose = .false.
    c_verbose = c_verbose_save
    ptc_private%taylor_order_ptc = t_order
  endif

  if (ptc_private%init_spin_needed) then
    call init_all (ptc_private%base_state, t_order, 0)
    ptc_private%init_spin_needed = .false.
  endif
endif

! Superkill tells PTC to do a through cleanup when killing a fibre.

SUPERKILL = .false.
C_WATCH_USER = .false.       ! Suppress some print statements 

!

use_info   = .true.    ! Enable storage space in fibre%i
use_info_m = .true.    ! Enable matrix storage in fibre%i%m

end subroutine set_ptc
