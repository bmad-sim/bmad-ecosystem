!+
! Subroutine set_ptc (param, taylor_order, integ_order, &
!                               num_steps, no_cavity, exact_calc)
!
! Subroutine to initialize ptc.
! Note: This subroutine cannot be used if there are "knobs".
! This subroutine replaces:
!     make_states
!     set_mad
!     init
!
! Modules needed:
!   use accelerator
!
! Input:
!   param        -- Param_struct, optional: BMAD parameters:
!     %energy      -- Energy.
!     %particle    -- Type of particle.
!   taylor_order -- Integer, optional: Maximum order of the taylor polynomials.
!   integ_order  -- Integer, optional: Default Order for the drift-kick-drift 
!                     sympletic integrator. Possibilities are: 2, 4, or 6
!                     Default = 2
!   num_steps    -- Integer, optional: Default Number of integration steps.
!                     Default = 1
!   no_cavity    -- Logical, optional: RF Cavity exists? 
!                     Default = False.
!   exact_calc   -- logical, optional: Exact Model? 
!                     Default = False.
!-

subroutine set_ptc (param, taylor_order, integ_&
                        order, num_steps, no_cavity, exact_calc) 

  use accelerator

  implicit none

  type (param_struct), optional :: param

  integer, optional :: integ_order, num_steps, taylor_order
  integer this_method, this_steps
  integer nd2, npara

  real(dp) this_energy

  logical, optional :: no_cavity, exact_calc
  logical, save :: init_needed = .true.
              
! do not call set_mad

  if (init_needed .and. present(param)) then
    if (param%particle == positron$ .or. param%particle == electron$) then
      call make_states(.true.)
    else
      call make_states(.false.)
    endif
    EXACT_MODEL = .false.
    ALWAYS_EXACTMIS = .false.
    init_needed = .false.
  endif

  if (present(exact_calc)) then
    EXACT_MODEL = exact_calc
    ALWAYS_EXACTMIS = exact_calc
  endif
    
  if (present(no_cavity)) default = default+nocavity
  
  if (present (integ_order)) then
    this_method = integ_order
    bmad_com%default_integ_order = integ_order
  else
    this_method = bmad_com%default_integ_order
  endif

  if (present (num_steps)) then
    this_steps = num_steps
    bmad_com%default_num_steps = num_steps
  else
    this_steps = bmad_com%default_num_steps
  endif

  if (present(param)) then
    if (bmad_com%energy /= param%energy .or. &
                        present(integ_order) .or. present(num_steps)) then
      this_energy = param%energy
      call set_mad (energy = this_energy, method = this_method, &
                                                       step = this_steps)
      bmad_com%energy  = param%energy
    endif
  endif

  if (present(taylor_order)) then  
    if (bmad_com%taylor_order_ptc /= taylor_order) then
      call init (default, taylor_order, 0, berz, nd2, &
                                               bmad_com%real_8_map_init)
      bmad_com%taylor_order_ptc = taylor_order
    endif
  endif
  
end subroutine  
