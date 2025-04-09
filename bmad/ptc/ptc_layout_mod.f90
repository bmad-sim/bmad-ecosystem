!+
! Module ptc_layout_mod
!
! Module of PTC layout interface routines.
! Also see: ptc_interface_mod
!-

module ptc_layout_mod

use ptc_interface_mod
use changed_attribute_bookkeeper

implicit none

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine type_ptc_layout (lay)
!
! Subroutine to print the global information in a layout
!
! Input:
!   lay - layout: layout to use.
!+

subroutine type_ptc_layout (lay)

use s_def_all_kinds, only: layout

type (layout) lay

!

if (.not. associated(lay%start)) then
  print *, 'Warning from TYPE_LAYOUT: Layout NOT Associated'
  return
endif

print *, 'Name:         ', lay%name
print *, 'N:            ', lay%N,        '  ! Number of Elements'
print *, 'LatPos:       ', lay%lastpos,  '  ! Last position'

end subroutine type_ptc_layout

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine branch_to_ptc_m_u (branch)
!
! Subroutine to create a PTC layout from a Bmad lattice branch.
! Note: If lat_to_ptc_layout has already been setup, you should first do a 
!           call kill_ptc_layouts(lat)
! This deallocates the pointers in PTC
!
! Note: If not already done, before you call this routine you need to first call:
!    call set_ptc (...)
! [This is normally done in bmad_parser.]
!
! Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there 
! will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
! In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the 
! exit end of ele will correspond to the exit end of ele%ptc_fibre.
!
! Input:
!   branch                -- branch_struct: Input branch.
!
! Output:
!   branch(:)%ptc              -- Pointers to generated layouts.
!   branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine branch_to_ptc_m_u (branch)

use s_fibre_bundle, only: ring_l, append, lp, layout, fibre
use mad_like, only: set_up, kill, lielib_print
use madx_ptc_module, only: m_u, m_t, append_empty_layout, survey, make_node_layout

type (branch_struct) :: branch
type (ele_struct) drift_ele
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: chain_ele(:)

integer n, ib, ie, ix_pass, ix_pass0, n_links
logical err_flag, doneit, ele_inserted_in_layout

! Transfer elements.

ix_pass0 = 0
ele_inserted_in_layout = .false.

call layout_init_stuff (branch)
call add_ptc_layout_to_list (branch%ptc,  m_u%end)   ! Save layout

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)

  call multipass_chain(ele, ix_pass, n_links, chain_ele)
  if (ix_pass > 1) then
    ele%ptc_fibre => chain_ele(1)%ele%ptc_fibre
    ix_pass0 = ix_pass
    cycle
  endif

  ! Use new layout for multipass regions.

  if (ix_pass0 /= ix_pass .and. ele_inserted_in_layout) then
    call layout_end_stuff (branch, ie, .false.)
    call layout_init_stuff (branch)
    call add_ptc_layout_to_list (branch%ptc, m_u%end)   ! Save layout
    ele_inserted_in_layout = .false.
  endif

  call ele_to_fibre (ele, ele%ptc_fibre, .true., err_flag, for_layout = .true.)

  ele_inserted_in_layout = .true.
  ix_pass0 = ix_pass

enddo

call layout_end_stuff (branch, ie, .true.)

! Set bookkeeping state

do ie = 0, branch%n_ele_max
  branch%ele(ie)%bookkeeping_state%ptc = ok$
enddo
branch%param%bookkeeping_state%ptc = ok$

!-----------------------------------------------------------------------------
contains

subroutine layout_init_stuff (branch)

type (branch_struct) branch

!

call append_empty_layout(m_u)
call set_up(m_u%end)

write (m_u%end%name, '(a, i4)') 'm_u Bmad branch:', branch%ix_branch ! Diagnostic

end subroutine layout_init_stuff

!-----------------------------------------------------------------------------
! contains

subroutine layout_end_stuff (branch, ie, at_end)

type (branch_struct) branch
logical at_end
integer ie

!

if (branch%param%geometry == closed$ .and. ie == branch%n_ele_track) then
  m_u%end%closed = .true.
else
  m_u%end%closed = .false.
endif

call set_ptc_quiet(12, set$, n)

m_u%end%closed = .true.
doneit = .true.
call ring_l (m_u%end, doneit)
call make_node_layout (m_u%end)
call survey (m_u%end)

call set_ptc_quiet(12, unset$, n)

end subroutine layout_end_stuff

end subroutine branch_to_ptc_m_u

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   
! 
! Routine to add a layout the a list of layouts.
!
! Input:
!   branch_ptc_info  -- ptc_branch1_struct: List of layouts
!   layout_end       -- layout: ptc layout
!
! Output:
!   branch_ptc_info  -- ptc_branch1_struct: Updated list.
!-

subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   

type (ptc_branch1_struct) branch_ptc_info, temp_info
type (layout), target :: layout_end
integer n

!

if (allocated(branch_ptc_info%m_u_layout)) then
  n = size(branch_ptc_info%m_u_layout)
  call move_alloc(branch_ptc_info%m_u_layout, temp_info%m_u_layout)
  allocate(branch_ptc_info%m_u_layout(n+1))
  branch_ptc_info%m_u_layout(1:n) = temp_info%m_u_layout
  branch_ptc_info%m_u_layout(n+1)%ptr => layout_end
  deallocate(temp_info%m_u_layout)

else
  allocate(branch_ptc_info%m_u_layout(1))
  branch_ptc_info%m_u_layout(1)%ptr => layout_end  
endif 

end subroutine add_ptc_layout_to_list

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_setup_tracking_with_damping_and_excitation (branch, include_damping, include_excitation, ptc_state, closed_orb)
!
! Routine to setup PTC tracking that includes radiation damping and stochastic excitation.
! To track use the routine track_probe with ptc_state as the state argument.
!
! Input:
!   branch              -- branch_struct: Lattice branch to setup.
!   include_damping     -- logical: Include radiation damping?
!   include_excitation  -- logical: Include radiation excitation?
!
! Output:
!   ptc_state           -- internal_state: Use with track_probe when tracking.
!   closed_orb(6)       -- real(dp): Closed orbit at start of branch. Only set if excitation is included.
!-

subroutine ptc_setup_tracking_with_damping_and_excitation (branch, include_damping, include_excitation, ptc_state, closed_orb)

use s_fitting_new

type (branch_struct) branch
type (internal_state) ptc_state
type (probe) prb
type (probe_8) prb8
type (c_damap) cda

real(dp) closed_orb(6)
logical include_damping, include_excitation

! If including excitation then need to do a setup calculation with envelope on.

if (include_excitation) then
  call alloc(prb8)
  call alloc(cda) 

  ptc_state = ptc_private%base_state + envelope0 + radiation0 
  call find_orbit_x (branch%ptc%m_t_layout, closed_orb, ptc_state, 1d-8, fibre1 = 1)

  cda = 1
  prb = closed_orb
  prb8 = prb + cda
  compute_stoch_kick = .true.
  call track_probe (prb8, ptc_state, fibre1 = branch%ele(1)%ptc_fibre)
  compute_stoch_kick = .false.

  call kill(prb8)
  call kill(cda)

else
  ptc_state = ptc_private%base_state 
  if (include_damping) ptc_state = ptc_state + radiation0
  call find_orbit_x (branch%ptc%m_t_layout, closed_orb, ptc_state, 1d-8, fibre1 = 1)
endif

! Now set ptc_state without envelope on.

ptc_state = ptc_private%base_state
if (include_damping)    ptc_state = ptc_state + radiation0
if (include_excitation) ptc_state = ptc_state + stochastic0

end subroutine ptc_setup_tracking_with_damping_and_excitation

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, pz)
!
! Routine to compute the transfer matrices for the individual elements and closed orbit 
! for a lattice branch with closed geometry.
!
! Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.
!
! Input:
!   branch            -- branch_struct: Lattice branch.
!   pz                -- real(rp), optional: energy offset around which to 
!                          calculate the matrices if there is no RF.
! Output:
!   branch            -- branch_struct: Lattice branch containing the matrices.
!     %ele(i)%fibre%i%m(6,6)    -- matrices.
!     %ele(i)%fibre%i%fix0(6)   -- Closed orbit at entrance.
!     %ele(i)%fibre%i%fix(6)    -- Closed orbit at exit.
!-

subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, pz)

use s_fitting_new, only: compute_linear_one_magnet_maps, fibre, layout

type (branch_struct) branch
type (fibre), pointer :: ptc_fibre
type (layout), pointer :: ptc_layout

real(rp), optional :: pz
real(rp) c_mat(6,6), c_mat_inv(6,6)
integer i

!

ptc_fibre => branch%ele(1)%ptc_fibre
call compute_linear_one_magnet_maps (ptc_fibre, ptc_private%base_state, pz)

end subroutine ptc_one_turn_mat_and_closed_orbit_calc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)
!
! Routine to calculate emittances, etc.
!
! Note: This routine calls the PTC init_all routine.
!
! Input: 
!   ele -- ele_struct: Element at which to evaluate the parameters.
!
! Output:
!   norm_mode       -- Normal_modes_struct
!     %a%tune, %b%tune, %z%tune
!     %a%alpha_damp, etc.
!     %a%emittance, etc.
!   sigma_map(6,6)  -- real(rp): Sigma matrix (Bmad coordinates).
!   closed_orb      -- coord_struct: Closed orbit at ele (Bmad coordinates).
!                        Notice: This closed orbit is the closed orbit with radiation on.
!-

subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)

use pointer_lattice

type (ele_struct) ele
type (normal_modes_struct) norm_mode
type (coord_struct) closed_orb
type (fibre), pointer :: ptc_fibre
type (branch_struct), pointer :: branch
type (internal_state) ptc_state
type (layout), pointer :: ptc_layout
type (c_damap) id
type (probe_8) xs
type (probe) xs0
type (c_normal_form) cc_norm

real(rp) sigma_mat(6,6), emit(3), ptc_sigma_mat(6,6), tune(3), damp(3), energy_loss, dp_loss, epsc
complex(rp) cmplx_sigma_mat(6,6)

character(*), parameter :: r_name = 'ptc_emit_calc'

!

branch => pointer_to_branch(ele)
if (.not. rf_is_on(branch)) then
  call out_io (s_error$, r_name, 'RF is not on! Cannot compute emittances in PTC.')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (.not. bmad_com%radiation_damping_on) then
  call out_io (s_info$, r_name, 'Note!! Radiation damping is always on for PTC emittance calc.')
endif

!

ptc_state = ptc_private%base_state - NOCAVITY0 + RADIATION0

ptc_fibre => pointer_to_fibre(ele)
ptc_layout => ptc_fibre%parent_layout
call find_orbit_x (closed_orb%vec, ptc_state, 1e-8_rp, fibre1 = ptc_fibre) 

ptc_state = ptc_state + ENVELOPE0

call alloc (id)
call alloc (xs)
call alloc (cc_norm)

id=1
xs0=closed_orb%vec
xs=xs0+id
call track_probe(xs, ptc_state, fibre1 = ptc_fibre)
id=xs
call GET_loss(ptc_layout, energy_loss, dp_loss)
norm_mode%e_loss = energy_loss * 1e9_rp

lielib_print(16) = 0    ! Do not print eigenvalue info.
epsc = EPS_EIGENVALUES_OFF_UNIT_CIRCLE
EPS_EIGENVALUES_OFF_UNIT_CIRCLE = max(1d-1, epsc)  ! Set larger since rad damping is on.

call ptc_set_rf_state_for_c_normal(ptc_state%nocavity)
call c_normal(id, cc_norm)
lielib_print(16) = 1
EPS_EIGENVALUES_OFF_UNIT_CIRCLE = epsc

call init_coord(closed_orb, closed_orb, ele, downstream_end$)

norm_mode%a%tune = cc_norm%tune(1)   ! Fractional tune with damping
norm_mode%b%tune = cc_norm%tune(2)
norm_mode%z%tune = cc_norm%tune(3)

norm_mode%a%alpha_damp = cc_norm%damping(1)
norm_mode%b%alpha_damp = cc_norm%damping(2)
norm_mode%z%alpha_damp = cc_norm%damping(3)

norm_mode%a%emittance = cc_norm%emittance(1)
norm_mode%b%emittance = cc_norm%emittance(2)
norm_mode%z%emittance = cc_norm%emittance(3)

sigma_mat = cc_norm%s_ij0

!

call kill (id)
call kill (xs)
call kill (cc_norm)

end subroutine ptc_emit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_spin_calc (ele, norm_mode, sigma_mat, closed_orb)
!
! Routine to equilibrium polarizations, etc.
!
! Input: 
!   ele -- ele_struct: Element at which to evaluate the parameters.
!
! Output:
!   norm_mode       -- Normal_modes_struct
!     %a%tune, %b%tune, %z%tune
!     %a%alpha_damp, etc.
!     %a%emittance, etc.
!   sigma_map(6,6)  -- real(rp): Sigma matrix (Bmad coordinates).
!   closed_orb      -- coord_struct: Closed orbit at ele (Bmad coordinates).
!                        Notice: This closed orbit is the closed orbit with radiation on.
!-

subroutine ptc_spin_calc (ele, norm_mode, sigma_mat, closed_orb)

use pointer_lattice

type (ele_struct) ele
type (normal_modes_struct) norm_mode
type (coord_struct) closed_orb
type (fibre), pointer :: ptc_fibre
type (branch_struct), pointer :: branch
type (internal_state) ptc_state
type (layout), pointer :: ptc_layout
type (c_damap) id, U_1, as, a0, a1, a2
type (probe_8) xs
type (probe) xs0
type (c_normal_form) cc_norm
type (c_spinor) isf

real(rp) sigma_mat(6,6), emit(3), ptc_sigma_mat(6,6), tune(3), damp(3), energy_loss, dp_loss
real(rp) depol, n0(3), dn_dpz(3)
complex(rp) cmplx_sigma_mat(6,6)

integer i1, i2, i3

character(*), parameter :: r_name = 'ptc_spin_calc'

!

branch => pointer_to_branch(ele)
if (.not. rf_is_on(branch)) then
  call out_io (s_error$, r_name, 'RF is not on! Cannot compute emittances in PTC.')
  if (global_com%exit_on_error) call err_exit
  return
endif

!

ptc_state = ptc_private%base_state - NOCAVITY0 + RADIATION0 + SPIN0

ptc_fibre => pointer_to_fibre(ele)
ptc_layout => ptc_fibre%parent_layout
call find_orbit_x (closed_orb%vec, ptc_state, 1e-8_rp, fibre1 = ptc_fibre) 

ptc_state = ptc_state + ENVELOPE0

call alloc (id, U_1, as, a0, a1, a2)
call alloc (xs)
call alloc (cc_norm)
call alloc (isf)


id=1
xs0=closed_orb%vec
xs=xs0+id
call track_probe(xs, ptc_state, fibre1 = ptc_fibre)
id=xs
call GET_loss(ptc_layout, energy_loss, dp_loss)
norm_mode%e_loss = energy_loss * 1e9_rp

call ptc_set_rf_state_for_c_normal(ptc_state%nocavity)
call c_normal(id, cc_norm, dospin = ptc_state%spin)
call c_full_canonise(cc_norm%atot, U_1, as, a0, a1, a2)

isf = 2   ! isf = (0,1,0)
call makeSO3(as)
isf = as%s * isf 
!m_spinor=1
!l_spinor=3
!l_spinor=as%s*l_spinor
!m_spinor=as%s*m_spinor 

depol=0    ! tao_dep^-1. See Bmad manual Eq. 20.21
do i1=1,3
do i2=1,6
do i3=1,6
  depol =  depol + (ISF%v(i1).d.i2) *  id%e_ij(i2,i3) * (ISF%v(i1).d.i3)
enddo
enddo
enddo
 
depol = -depol/2

do i1 = 1, 3
  dn_dpz(i1) = isf%v(i1).d.6
  n0(i1) = isf%v(i1)
enddo

call init_coord(closed_orb, closed_orb, ele, downstream_end$)

norm_mode%a%tune = cc_norm%tune(1)   ! Fractional tune with damping
norm_mode%b%tune = cc_norm%tune(2)
norm_mode%z%tune = cc_norm%tune(3)

norm_mode%a%alpha_damp = cc_norm%damping(1)
norm_mode%b%alpha_damp = cc_norm%damping(2)
norm_mode%z%alpha_damp = cc_norm%damping(3)

norm_mode%a%emittance = cc_norm%emittance(1)
norm_mode%b%emittance = cc_norm%emittance(2)
norm_mode%z%emittance = cc_norm%emittance(3)

sigma_mat = cc_norm%s_ij0

!

call kill (id, U_1, as, a0, a1, a2)
call kill (xs)
call kill (cc_norm)
call kill (isf)

end subroutine ptc_spin_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_track_all (branch, orbit, track_state, err_flag)
!
! Routine to track from the start to the end of a lattice branch. 
! 
! Input:
!   branch      -- lat_struct: Lat to track through.
!   orbit(0)    -- Coord_struct: Coordinates at beginning of branch.
!
! Output:
!   orbit(0:*)   -- Coord_struct: Orbit array.
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!   err_flag     -- Logical, optional: Set true if particle lost or error. False otherwise
!-

subroutine ptc_track_all (branch, orbit, track_state, err_flag)

use madx_ptc_module

type (branch_struct), target :: branch
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: ele
type (fibre), pointer :: fib

real(rp) x(6)

integer i
integer, optional :: track_state

logical, optional :: err_flag

! Init

if (present(err_flag)) err_flag = .true.
if (present(track_state)) track_state = moving_forward$

!

if (orbit(0)%state == not_set$) call init_coord(orbit(0), orbit(0)%vec, branch%ele(0), downstream_end$) 
x = orbit(0)%vec

do i = 1, branch%n_ele_track
  ele => branch%ele(i)

  orbit(i) = orbit(i-1)
  call check_aperture_limit (orbit(i), ele, first_track_edge$, branch%param)
  if (orbit(i)%state /= alive$) then
    if (present(err_flag)) err_flag = .false.
    return
  endif

  fib => branch%ele(i)%ptc_fibre%next
  call track_probe_x (x, ptc_private%base_state, branch%ele(i-1)%ptc_fibre%next, fib)
  call init_coord (orbit(i), x, ele, downstream_end$)

  call check_aperture_limit (orbit(i), ele, second_track_edge$, branch%param)
  if (orbit(i)%state /= alive$) then
    if (present(err_flag)) err_flag = .false.
    return
  endif

enddo

end subroutine ptc_track_all

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)
!
! Routine to calculate the closed orbit of a lattice branch using PTC.
! This routine assumes the associated PTC layout has been crated 
! with lat_to_ptc_layout.
!
! Input:
!   branch          -- branch_struct: Branch of a lattice.
!   radiation_damping_on
!                   -- logical, optional: If True, radiation dampling is included in the calculation. 
!                        Default is the setting of bmad_com%%radiation_damping_on.
!
! Output:
!   closed_orbit(:) -- coord_struct, allocatable: closed_orbit
!-

subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)

use madx_ptc_module

type (branch_struct), target :: branch
type (coord_struct), allocatable :: closed_orbit(:)
type (fibre), pointer :: fib
type (internal_state) ptc_state

real(rp) x(6)
integer i

logical, optional :: radiation_damping_on

!

if (logic_option(bmad_com%radiation_damping_on, radiation_damping_on)) then
  ptc_state = ptc_private%base_state + radiation0
else
  ptc_state = ptc_private%base_state - radiation0
endif

call reallocate_coord(closed_orbit, branch%n_ele_max)

x = 0
fib => branch%ele(0)%ptc_fibre%next
call find_orbit_x (x, ptc_state, 1.0d-7, fibre1 = fib)  ! find closed orbit
call init_coord (closed_orbit(0), x, branch%ele(0), downstream_end$)

do i = 1, branch%n_ele_track
  fib => branch%ele(i)%ptc_fibre%next
  call track_probe_x (x, ptc_state, branch%ele(i-1)%ptc_fibre%next, fib)
  call init_coord (closed_orbit(i), x, branch%ele(i), downstream_end$)
enddo

end subroutine ptc_closed_orbit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_one_turn_map_at_ele (ele, orb0, map, ptc_state, pz, rf_on)
!
! Routine to calculate the PTC one turn map for a ring.
!
! Input:
!   ele     -- ele_struct: Element determining start/end position for one turn map.
!   map     -- probe_8: Will be allocated.
!   pz      -- real(rp), optional: momentum deviation of closed orbit. 
!                                  Default = 0
!   rf_on   -- integer, optional: RF state for calculation. yes$, no$, or maybe$ (default)
!                   maybe$ means that RF state in branch is used.
!
! Output:
!   orb0(6)     -- real(rp): Closed orbit around which map is made.
!   map         -- probe_8: Map.
!     %x           -- Orbital part.
!     %q%x         -- Spin part.
!   ptc_state   -- internal_state: PTC state used for tracking.
!-

subroutine ptc_one_turn_map_at_ele (ele, orb0, map, ptc_state, pz, rf_on)

use madx_ptc_module

type (ele_struct), target :: ele
type (fibre), pointer :: fib
type (c_damap) da_map
type (probe_8) map
type (internal_state) ptc_state
type (probe) p0

real(rp), optional :: pz
real(dp) orb0(6)

integer, optional :: rf_on

logical rf_on_state, spin_on

!

spin_on = bmad_com%spin_tracking_on

select case (integer_option(maybe$, rf_on))
case (yes$)
  ptc_state = ptc_private%base_state - nocavity0
case (no$)
  ptc_state = ptc_private%base_state + nocavity0
case (maybe$)
  rf_on_state = rf_is_on(ele%branch)
  if (rf_on_state) then
    ptc_state = ptc_private%base_state - nocavity0
  else
    ptc_state = ptc_private%base_state + nocavity0
  endif
case default
  stop
end select

if (spin_on) ptc_state = ptc_state + spin0

! Find closed orbit

orb0 = 0
if (present(pz)) orb0(6) = pz
fib => ele%ptc_fibre%next
call find_orbit_x (orb0, ptc_state, 1.0d-5, fibre1 = fib)  ! find closed orbit

! Construct map.

call alloc(da_map)
call alloc(map)

p0 = orb0
da_map = 1
map = da_map + p0
call track_probe (map, ptc_state, fibre1 = fib)

call kill(da_map)

end subroutine ptc_one_turn_map_at_ele

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_map_to_normal_form (one_turn_map, ptc_nocavity, normal_form, phase_map, spin_tune)
!
! Routine to do normal form analysis on a map.
! Note: All output quantities must be allocated prior to calling this routine.
!
! Input:
!   one_turn_map    -- probe_8: One turn map.
!   ptc_nocavity    -- logical: PTC nocavity parameter setting when map was made.
!
! Output:
!   normal_form     -- c_normal_form: Normal form decomposition.
!   phase_map(3)    -- c_taylor: Phase Taylor maps.
!   spin_tune       -- c_taylor, optional: Spin tune Taylor map.
!-

subroutine ptc_map_to_normal_form (one_turn_map, ptc_nocavity, normal_form, phase_map, spin_tune)

use pointer_lattice

implicit none

type (probe_8) one_turn_map
type (c_normal_form) normal_form
type (c_taylor) phase_map(3)
type (c_taylor), optional :: spin_tune
type (c_damap) c_map

logical ptc_nocavity

!! type (c_damap) as, a2, a1, a0
!! integer n(6)

!

call alloc(c_map)
c_map = one_turn_map

call ptc_set_rf_state_for_c_normal(ptc_nocavity)
call c_normal (c_map, normal_form, dospin = bmad_com%spin_tracking_on, phase = phase_map, nu_spin = spin_tune)

call kill(c_map)


!! call c_full_factor_map (normal_form%atot, as, a0, a1, a2)  ! as = spin, a0 = closed orbit, a1 = linear, a2 = non-linear
!! call clean(a0, a0, prec = 1d-10)
!! call clean(a1, a1, prec = 1d-10)
!! call clean(a2, a2, prec = 1d-10)
!! call clean(as, as, prec = 1d-10)
!! 
!! call print (a0)   ! a_tot = a_0 o a_1 o a_2 o a_s
!! print *, '===================================='
!! call print (a1)
!! n = 0
!! n(6) = 3
!! print *, (a0%v(1).sub.n)
!! print *, (normal_form%atot%v(1).sub.n)
!! call kill(c_map, as, a2, a1, a0)

!! call print (normal_form%atot)

end subroutine ptc_map_to_normal_form

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_form_taylors(one_turn_taylor, rf_on, dhdj, A, A_inverse)
!
! Do a normal form decomposition on a one-turn taylor map M:
!   M = A o R o A_inverse
! where A maps Floquet (fully normalized) coordinates to lab coordinates. 
! In Floquet coordinates, the amplitudes are defined as J_i = (1/2) (x_i^2 + p_i^2).
! The map R = exp(:h:) is a pure rotation with h = h(J) is a function of the amplitudes only.
! The angles (phase advances) are given by phi_i = 2pi*dh/dJ_i.
! The taylor terms of dhdj are therefore the tunes, chromaticities, amplitude dependent tune shifts, etc.
!
! The mapping procedure for one turn is:
!  z_Floquet_in = A_inverse o z_Lab_in
!  [phi_a, phi_b, phi_c] = 2 pi * dhdj o z_Floquet_in
!  z_Floquet_out = RotationMatrix(phi_a, phi_b, phi_c) . z_Floquet_in
!  z_Lab_out = A o z_Floquet_out
!
!
! Input: 
!   one_turn_taylor -- taylor_struct      : one turn taylor map
!   rf_on           -- logical            : Was the map calculated with RF on?
!
! Output:
!   A             -- taylor_struct, optional: Map from Floquet coordinates to Lab coordinates
!   A_inverse     -- taylor_struct, optional: Map from Lab coordinates to Floquet coordinates
!   dhdj            -- taylor_struct, optional: Map from Floquet coordinates to phase advances
!-
subroutine normal_form_taylors (one_turn_taylor, rf_on, dhdj, A, A_inverse)

use madx_ptc_module

type (taylor_struct) :: one_turn_taylor(6)
type (taylor_struct), optional :: A_inverse(6), dhdj(6), A(6)
type (damap) :: da_map
type (real_8) :: map8(6)
type (normalform) :: normal
type (internal_state) :: state
logical :: rf_on
!

if (rf_on) then
  state = ptc_private%base_state - nocavity0
else
  state = ptc_private%base_state + nocavity0
endif

call alloc(map8)
call alloc(da_map)
call alloc(normal)

! Convert to real_8, then a damap
map8 = one_turn_taylor

! Get the normal form
da_map = map8
normal = da_map

! Convert to taylor_structs

! A
if (present(A)) then
  A(1:6) = normal%A_t%v(1:6)
endif

! A_inverse
if (present(A_inverse)) then
  map8 = normal%A_t**(-1)
  A_inverse = map8
endif

! dhdj 
if (present(dhdj)) then
  dhdj(1:6) = normal%dhdj%v(1:6)
endif

! Cleanup

call kill(map8)
call kill(da_map)
call kill(normal)

end subroutine normal_form_taylors

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_form_complex_taylors
!
! UNDER DEVELOPMENT
!-

subroutine normal_form_complex_taylors (one_turn_taylor, rf_on, F, L, A, A_inverse, order)

use madx_ptc_module

type (taylor_struct) :: one_turn_taylor(6)
type (complex_taylor_struct), optional :: F(6), L(6)
type (taylor_struct), optional ::  A(6), A_inverse(6)
type (c_damap) :: cda, cdaLinear
type (damap) :: da
type (real_8) :: map8(6)
type (c_normal_form) :: complex_normal_form
type(c_vector_field) :: fvecfield
type (internal_state) :: state
integer :: i
integer, optional :: order
logical :: rf_on, c_verbose_save

!

if (rf_on) then
  state = ptc_private%base_state - nocavity0
else
  state = ptc_private%base_state + nocavity0
endif

!

c_verbose_save = c_verbose
c_verbose = .false.

call alloc(map8)
call alloc(da)
call alloc(cda)
call alloc(cdaLinear)
call alloc(fvecfield)
call alloc(complex_normal_form)

! Convert to real_8, then a da map, then complex da map
map8 = one_turn_taylor
da = map8
cda = da

! Complex normal form in phasor basis
! See: fpp-ptc-read-only/build_book_example_g95/the_fpp_on_line_glossary/complex_normal.htm
! M = A o N o A_inverse.

call ptc_set_rf_state_for_c_normal(state%nocavity)
call c_normal(cda, complex_normal_form, dospin=my_false, no_used=integer_option(1, order))

if (present(F) .or. present(L)) then
  cda = complex_normal_form%N

  ! Move to the phasor basis
  cda=from_phasor(-1)*cda*from_phasor(1)

  !   N = L exp(F.grad), where L is the linear (rotation) map.
  call c_factor_map(cda, cdaLinear, fvecfield, 1) ! 1 => Dragt-Finn direction

  ! Zero out small coefficients
  call c_clean_vector_field(fvecfield, fvecfield, 1.d-8 )
  ! Output
endif

if (present(L)) then
  L(1:6) = cdaLinear%v(1:6)
endif

if (present(F)) then
  F(1:6) = fvecfield%v(1:6)
endif

if(present(A)) then
  da = complex_normal_form%a_t
  A(1:6) = da%v(1:6)
endif

if(present(A_inverse)) then
  da = complex_normal_form%a_t**(-1)
  A_inverse(1:6) = da%v(1:6)
endif

! Cleanup
call kill(map8)
call kill(da)
call kill(cda)
call kill(cdaLinear)
call kill(fvecfield)
call kill(complex_normal_form)

! Reset PTC state
use_complex_in_ptc=my_false
c_verbose = c_verbose_save

end subroutine normal_form_complex_taylors

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_form_rd_terms(one_turn_map, normal_form, rf_on, order)
!
! Calculate resonance driving terms using PTC.  The contents of
! normal_form%rd_term(:)%c_val are directly comparable to equations
! 97 of "The Sextupole Scheme for the Swiss Light Source (SLS): An
! An Analytic Approach" SLS Note 9/97, Johan Bengtsson.
!
! Input:
!   one_turn_map   -- probe_8: One-turn map.
!   rf_on          -- logical: perform calculation with RF on?
!   order          -- integer, optional: order for normal_form_calculation.
!
! Output:
!   normal_form%h(:)%c_val -- bmad_normal_form_struct: complex values for one-turn driving terms.
!-
subroutine normal_form_rd_terms(one_turn_map, normal_form, rf_on, order)

use madx_ptc_module

type (probe_8) one_turn_map
type (bmad_normal_form_struct) normal_form
integer i, order_for_normal_form
integer, optional :: order
logical :: rf_on, c_verbose_save

type (c_vector_field) F
type (c_taylor) vb
type (c_damap) cda
type (c_damap) a2, a_step1
type (c_damap) acs1_a
type (c_normal_form) complex_normal_form
type (internal_state) :: state

!

if (.not. allocated(normal_form%h)) return

order_for_normal_form = integer_option(ptc_private%taylor_order_ptc, order)

if (rf_on) then
  state = ptc_private%base_state - nocavity0
else
  state = ptc_private%base_state + nocavity0
endif

! Set PTC state

c_verbose_save = c_verbose
c_verbose = .false.

call alloc(vb)
call alloc(F)
call alloc(cda)
call alloc(complex_normal_form)
call alloc(a2)
call alloc(a_step1)
call alloc(acs1_a)

! Complex normal form in phasor basis
! See: fpp-ptc-read-only/build_book_example_g95/the_fpp_on_line_glossary/complex_normal.htm
! M = A o N o A_inverse.

cda = one_turn_map

call ptc_set_rf_state_for_c_normal(state%nocavity)
call c_normal(cda, complex_normal_form, dospin=my_false)  !, no_used=1) 
call c_fast_canonise(complex_normal_form%a_t,acs1_a)

a_step1 = acs1_a**(-1)*cda*acs1_a
a2 = (a_step1.sub.1)**(-1)*a_step1

F=log(a2) ! c_vector_field <- c_damap
vb = getpb(F)*c_phasor()

do i = 1, size(normal_form%h(:))
  if (normal_form%h(i)%id /= '') then
    normal_form%h(i)%c_val = vb .par. normal_form%h(i)%id
  endif
enddo

! Cleanup
call kill(vb)
call kill(F)
call kill(cda)
call kill(complex_normal_form)
call kill(a2)
call kill(a_step1)
call kill(acs1_a)

! Reset PTC state

use_complex_in_ptc=my_false
c_verbose = c_verbose_save

end subroutine normal_form_rd_terms

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+

subroutine set_ptc_verbose(on)

use madx_ptc_module

logical :: on

!

c_verbose = on

end subroutine set_ptc_verbose

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine update_ele_from_fibre (ele)
!
! Routine to update a bmad lattice element when the associated PTC fibre has been modified.
! Remember to call lattice_bookkeeper after calling this routine.
!
! Input:
!   ele           -- ele_struct: Element with corresponding ele%ptc_fibre fibre.
!
! Output:
!   ele       -- ele_struct: Modified element. 
!-

subroutine update_ele_from_fibre (ele)

use precision_constants, only: volt_c

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (fibre), pointer :: fib

real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx), tilt, kick
integer nmul, i, ix
character(40) name

!

fib => ele%ptc_fibre
branch => pointer_to_branch(ele)

call update_this_real (ele%value(l$), fib%mag%p%ld)
call update_this_real (ele%value(p0c$), 1e9_rp * fib%mag%p%p0c)

if (attribute_name(ele, num_steps$) == 'NUM_STEPS') then
  call update_this_real (ele%value(num_steps$), real(fib%mag%p%nst, rp))
  if (ele%value(num_steps$) == 0) then
    ele%value(ds_step$) = 0
  else
    ele%value(ds_step$) = ele%value(l$) / ele%value(num_steps$)
  endif
endif

! If integrator_order is defined for this element then update

name = attribute_name(ele, integrator_order$)
if (name(1:1) /= '!') call update_this_real (ele%value(integrator_order$), real(fib%mag%p%method, rp))

! Multipole

a_pole = 0
b_pole = 0
a_pole_elec = 0
b_pole_elec = 0

if (associated(fib%mag%an)) then
  nmul = size(fib%mag%an)
  a_pole(0:nmul-1) = fib%mag%an
  b_pole(0:nmul-1) = fib%mag%bn
endif

call multipole_ab_to_kt (a_pole, b_pole, knl, tn)

! Electric Multipole

if (associated(fib%mag%tp10)) then
  if (associated(fib%mag%tp10%ae)) then
    if (.not. associated(ele%a_pole_elec)) call multipole_init(ele, electric$)
    nmul = size(fib%mag%tp10%ae)
    a_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%ae(1:n_pole_maxx+1)
    b_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%be(1:n_pole_maxx+1)
  endif
endif

!

if (ele%key == sbend$) then
  call update_this_real (ele%value(ref_tilt_tot$), fib%mag%p%tiltd)
else
  call update_this_real (ele%value(tilt_tot$), fib%mag%p%tiltd)
endif

!

select case (ele%key)

case (octupole$)
  call update_this_real (ele%value(k3$), knl(3))
  call update_this_real (ele%value(tilt$), tn(3))
  a_pole(3) = 0
  b_pole(3) = 0

case (quadrupole$)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  a_pole(1) = 0
  b_pole(1) = 0

case (hkicker$)
  call update_this_real (ele%value(kick$), knl(1))
  a_pole(1) = 0
  b_pole(1) = 0

case (lcavity$, rfcavity$)
  call update_this_real (ele%value(rf_frequency$), fib%mag%freq)

  select case (nint(ele%value(cavity_type$)))
  case (traveling_wave$);     call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  case (standing_wave$);      call update_this_real (ele%value(voltage$), fib%mag%volt*0.5d6)
  case (ptc_standard$);       call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  end select

  if (ele%key == lcavity$) then
    call update_this_real (ele%value(phi0$), pi/2 - fib%mag%lag/twopi)
  else
    call update_this_real (ele%value(phi0$), fib%mag%lag/twopi)
  endif

case (multipole$)
  ele%a_pole = knl
  ele%b_pole = tn
  a_pole = 0
  b_pole = 0

case (sbend$)
  call update_this_real (ele%value(g$), fib%mag%p%b0)
  call update_this_real (ele%value(angle$), ele%value(g$) * ele%value(l$))
  call update_this_real (ele%value(hgap$), fib%mag%hgap(1))
  call update_this_real (ele%value(fint$), fib%mag%fint(1))
  call update_this_real (ele%value(hgapx$), fib%mag%hgap(2))
  call update_this_real (ele%value(fintx$), fib%mag%fint(2))

  if (nint(ele%value(ptc_field_geometry$)) == straight$) then
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1) + ele%value(angle$)/2)
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2) + ele%value(angle$)/2)
  else
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1))
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2))
  endif

case (sextupole$)
  call update_this_real (ele%value(k2$), knl(2))
  call update_this_real (ele%value(tilt$), tn(2))
  a_pole(2) = 0
  b_pole(2) = 0

case (solenoid$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)

case (sol_quad$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  a_pole(1) = 0
  b_pole(1) = 0

case (wiggler$, undulator$)

case default
end select

! multipoles

if (any(a_pole /= 0) .or. any(b_pole /= 0)) then
  if (.not. associated(ele%a_pole)) call multipole_init(ele, magnetic$)
  ele%a_pole = a_pole
  ele%b_pole = b_pole
endif

if (any(a_pole_elec /= 0) .or. any(b_pole_elec /= 0)) then
  if (.not. associated(ele%a_pole_elec)) call multipole_init(ele, electric$)
  ele%a_pole_elec = a_pole_elec
  ele%b_pole_elec = b_pole_elec
endif

! Fringes

!------------------------------------------------------------------------
contains

subroutine update_this_real (var, value)

real(rp) var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_real

!------------------------------------------------------------------------
! contains

subroutine update_this_integer (var, value)

integer var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_integer

end subroutine update_ele_from_fibre

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max, 
!                                 even_steps, r_typical, dx_tol_bend, use_2nd_order)
!
! Routine to calculate the optimum number of tracking steps and order
! of the integrator (2, 4, or 6) for each fibre in a layout.
!
! See the Bmad manual chapter on PTC for more details.
!
! Input:
!   ptc_layout    -- layout:
!   kl_max        -- real(rp): Maximum K1*L per tracking step.
!   ds_max        -- real(rp): Maximum ds for any step. 
!                      Useful when including other physicas like space charge.
!   even_steps    -- logical, optional: Always use an even number of steps for a fibre?
!                      Useful if need to evaluate at the center of fibres.
!   r_typical     -- real(rp), optional: Typical transverse offset. Used for computing the
!                      effective contribution of K1*L due to sextupoles.
!   dx_tol_bend   -- real(rp): Tolerable residual orbit in a bend.
!   use_2nd_order -- logical, optional: If present and True then force the use of 2nd order
!                       integrator.
!   crossover(2)  -- integer, optional: crossover points between orders for all
!                       elements except wigglers. Default is [4, 18].
!   crossover_wiggler(2)
!                 -- integer, optional: crossover points for wigglers. Default is [30, 60].
!
! Output:
!   ptc_layout -- layout: Lattice with the optimum number of tracking steps and integrator order.
!-

subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max, &
                    even_steps, r_typical, dx_tol_bend, use_2nd_order, crossover, crossover_wiggler)

use madx_ptc_module

type (layout) ptc_layout

real(rp) kl_max
real(rp), optional :: ds_max, dx_tol_bend, r_typical

integer, optional :: crossover(2), crossover_wiggler(2)
integer :: limit_int(2), lp14

logical, optional :: even_steps(2), use_2nd_order

!

resplit_cutting = 0   ! Ignore ds_max
if (present(ds_max)) then
  if (ds_max > 0) resplit_cutting = 2
endif

limit_int = [4, 18]

if (logic_option(.false., use_2nd_order)) limit_int = [10000, 10001] ! Something big.

call set_ptc_quiet(14, set$, lp14)

call thin_lens_resplit (ptc_layout, kl_max, lim = crossover, &
            limit_wiggler = crossover_wiggler, lmax0 = ds_max, sexr = r_typical, xbend = dx_tol_bend)

call set_ptc_quiet(14, unset$, lp14)

end subroutine ptc_calculate_tracking_step_size

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx, 
!                                                           even, crossover, crossover_wiggler)
!
! Routine to resplit (that is, recalculate the number of integration steps for an element)
! For the fibres in all layouts. After doing a resplit, the tune (and any other relavent
! "adjustable" parameters) should be adjusted to the correct values.
!
! Input:
!   dKL_max      -- real(rp): Maximum K1 * L quadrupole strength allowed for an integration step.
!                     Reasonable value would be something like 0.04.
!   l_max        -- real(rp): Maximum step length. Ignored if set to 0.
!   l_max_drift_only
!                -- logical: If True then l_max is only used for splitting drifts.
!   bend_dorb    -- real(rp): Residual bend orbit error. With some integration methods a zero
!                     orbit at the start of the bend will not be zero at the end. In this case,
!                     bend_dorb sets a maximum allowable orbit deviation. If set to zero, 
!                     this argument will be ignored. A resonable value is 10d-7. Note that the 
!                     actual orbit deviation is not simply related to bend_dorb and can be larger.
!                     In any case, lowering bend_dorb (without making it zero) will lower the 
!                     orbit deviation.
!   sex_dx       -- real(rp): To split sextupoles, sex_dx is used as the reference position 
!                     about which the quadrupole strength is calculated. This quadrupole strength
!                     is then used with dKL_max to calculate the number of integration steps.
!                     Set to zero to ignore.
!   even         -- logical, optional: If True then each fibre  will have an even number of steps.
!                     If False then the number of steps will be odd. If not present then number
!                     of steps is not constrained to be even or odd.
!   crossover(2) -- integer, optional: crossover(1) sets the maximum number of 2nd order integration
!                     steps to use. If the number of steps would exceed crossover(1) then integration
!                     is switched to 4th order. crossover(2) sets the maximum number of 4th order
!                     integration steps. If this number is exceeded, 6th order integration is used.
!                     Currently the default in PTC is [4, 18].
!   crossover_wiggler(2)
!                 -- integer, optional: crossover for wiggler elements.
!-


subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx, &
                                                            even, crossover, crossover_wiggler)

use s_fitting, only: thin_lens_restart, thin_lens_resplit
use madx_ptc_module, only: m_u

type(layout), pointer :: r

real(rp) dKL_max, bend_dorb, sex_dx, l_max
integer, optional :: crossover(2), crossover_wiggler(2)
integer lp14
logical l_max_drift_only
logical, optional :: even

!

r => m_u%start

call set_ptc_quiet(14, set$, lp14)

call thin_lens_restart(r, universe=.true.)
call thin_lens_resplit(r, dKL_max, even, crossover, crossover_wiggler, l_max, bend_dorb, sex_dx, universe=.true.)

call set_ptc_quiet(14, unset$, lp14)

end subroutine ptc_layouts_resplit

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_check_for_lost_particle (state, ptc_fibre, do_reset)
!
! Routine to check if a particle has been lost when tracking with PTC.
!
! Input:
!   do_reset -- logical: If True then reset ptc flags.
!
! Output:
!   state     -- integer: Same as coord_struct%state. alive$, lost$, lost_neg_x$, etc.
!   ptc_fibre -- fibre, pointer: Pointer to fibre where particle lost. Nullified if particle alive.
!-

subroutine ptc_check_for_lost_particle (state, ptc_fibre, do_reset)

use definition, only: check_stable, lost_node, lost_fibre, xlost, reset_aperture_flag

type (fibre), pointer :: ptc_fibre
integer state
logical do_reset

!

if (check_stable) then
  state = alive$
  return
endif

!

ptc_fibre => lost_fibre
state = lost$

!! ptc_fibre%cylindrical_map%p%aperture%...
!! xlost(1:6) is lost position

!

if (do_reset) call reset_aperture_flag()

end subroutine ptc_check_for_lost_particle

end module
