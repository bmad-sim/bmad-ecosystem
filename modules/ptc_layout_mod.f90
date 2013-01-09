!+
! Module ptc_layout_mod
!
! Module of PTC layout interface routines.
! Also see: ptc_interface_mod
!-

module ptc_layout_mod

use ptc_interface_mod
use multipass_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine type_ptc_layout (lay)
!
! Subroutine to print the global information in a layout
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   lay - layout: layout to use.
!+

subroutine type_ptc_layout (lay)

use s_def_all_kinds, only: layout

implicit none

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
! Subroutine lat_to_ptc_layout (lat)
!
! Subroutine to create a PTC layout from a Bmad lattice branch.
! Note: If ptc_layout has been already used then you should first do a 
!           call kill(ptc_layout)
! This deallocates the pointers in the layout
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
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   lat -- lat_struct: Input lattice
!
! Output:
!   lat%branch(:)%ptc              -- Pointers to generated layouts.
!   lat%branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine lat_to_ptc_layout (lat)

use madx_ptc_module, only: m_u, m_t, fibre, append_empty_layout, survey, make_node_layout, append_point, &
                           set_up, ring_l

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: chain_ele(:)
type (fibre), pointer :: fib
type (layout), pointer :: lay

real(rp) ptc_orientation(3,3), ang(3)

integer i, j, ix_pass, n_links
logical logic

character(20), parameter :: r_name = 'lat_to_ptc_layout'

! Setup m_u

do i = 0, ubound(lat%branch, 1)
  call branch_to_ptc_layout (lat%branch(i))
enddo

! setup m_t

do i = 0, ubound(lat%branch, 1)

  call append_empty_layout(m_t)
  call set_up (m_t%end)
  write (m_t%end%name, '(a, i4)') 'm_t Bmad branch:', i  ! For diagnostic purposes

  branch => lat%branch(i)
  call add_ptc_layout_to_list(branch%ptc,  m_t%end)   ! Save layout

  lay => m_t%end

  do j = 0, branch%n_ele_track
    ele => branch%ele(j)

    if (.not. associated(ele%ptc_fibre)) then
      call out_io (s_fatal$, r_name, 'NO FIBRE ASSOCIATED WITH ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    endif

    if (tracking_uses_end_drifts(ele)) then
      call append_this_fibre(ele%ptc_fibre%previous%previous)
      call append_this_fibre(ele%ptc_fibre%previous)
    endif

    call append_this_fibre(ele%ptc_fibre, .true.)

  enddo

  lay%closed = .true.

  logic = .true.
  call ring_l(lay, logic)

  ele => branch%ele(0)
  ang = [ele%floor%phi, ele%floor%theta, ele%floor%psi]
  call bmad_patch_parameters_to_ptc (ang, ptc_orientation)

  call survey (m_t%end, ptc_orientation, ele%floor%r)

enddo

!-----------------------------------------------------------------------------
contains

subroutine append_this_fibre(ele_fib, do_point)

type (fibre), pointer :: ele_fib
logical, optional :: do_point

!

call append_point(lay, ele_fib)
fib => lay%end

if (ele%key == patch$ .or. ele%key == floor_shift$) then
  fib%dir = ele%value(ptc_dir$)
else
  fib%dir = ele%orientation
endif

fib%charge = ele%branch%param%rel_tracking_charge

if (logic_option(.false., do_point)) ele%ptc_fibre => fib

end subroutine append_this_fibre

end subroutine lat_to_ptc_layout 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine branch_to_ptc_layout (lat)
!
! Subroutine to create a PTC layout from a Bmad lattice branch..
! Note: If ptc_layout has been already used then you should first do a 
!           call kill(ptc_layout)
! This deallocates the pointers in the layout
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
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   branch -- branch_struct: Input branch.
!
! Output:
!   branch(:)%ptc              -- Pointers to generated layouts.
!   branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine branch_to_ptc_layout (branch)

use s_fibre_bundle, only: ring_l, append, lp, layout, fibre
use mad_like, only: set_up, kill, lielib_print
use madx_ptc_module, only: m_u, m_t, append_empty_layout, survey, make_node_layout

implicit none

type (branch_struct) :: branch
type (ele_struct) drift_ele
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: chain_ele(:)

integer n, ib, ie, ix_pass, ix_pass0
logical doneit, ele_inserted_in_layout

! Transfer elements.

ix_pass0 = 0
ele_inserted_in_layout = .false.

call layout_init_stuff
call add_ptc_layout_to_list (branch%ptc,  m_u%end)   ! Save layout

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)

  call multipass_chain(ele, ix_pass, chain_ele = chain_ele)
  if (ix_pass > 1) then
    ele%ptc_fibre => chain_ele(1)%ele%ptc_fibre
    ix_pass0 = ix_pass
    cycle
  endif

  ! Use new layout for multipass regions.

  if (ix_pass0 /= ix_pass .and. ele_inserted_in_layout) then
    call layout_end_stuff
    call layout_init_stuff
    call add_ptc_layout_to_list (branch%ptc, m_u%end)   ! Save layout
    ele_inserted_in_layout = .false.
  endif

  !

  if (tracking_uses_end_drifts(ele)) then
    call create_hard_edge_drift (ele, upstream_end$, drift_ele)
    call ele_to_fibre (drift_ele, drift_ele%ptc_fibre, branch%param, .true., for_layout = .true.)
  endif

  call ele_to_fibre (ele, ele%ptc_fibre, branch%param, .true., for_layout = .true.)

  if (tracking_uses_end_drifts(ele)) then
    call create_hard_edge_drift (ele, downstream_end$, drift_ele)
    ! ele%ptc_fibre points to last PTC fibre.
    call ele_to_fibre (drift_ele, ele%ptc_fibre, branch%param, .true., for_layout = .true.)
  endif

  ele_inserted_in_layout = .true.
  ix_pass0 = ix_pass

enddo

call layout_end_stuff

!-----------------------------------------------------------------------------
contains

subroutine layout_init_stuff ()

call append_empty_layout(m_u)
call set_up(m_u%end)

write (m_u%end%name, '(a, i4)') 'm_u Bmad branch:', branch%ix_branch ! Diagnostic

end subroutine layout_init_stuff

!-----------------------------------------------------------------------------
! contains

subroutine layout_end_stuff ()

if (branch%param%geometry == closed$ .and. ie == branch%n_ele_track) then
  m_u%end%closed = .true.
else
  m_u%end%closed = .false.
endif

n = lielib_print(12)
lielib_print(12) = 0  ! No printing info messages

m_u%end%closed = .true.
doneit = .true.
call ring_l (m_u%end, doneit)
call survey (m_u%end)
call make_node_layout (m_u%end)

lielib_print(12) = n

end subroutine layout_end_stuff

end subroutine branch_to_ptc_layout

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   
! 
! Routine to add a layout the a list of layouts.
!
! Module needed:
!   use ptc_layout_mod
!
! Input:
!   branch_ptc_info  -- ptc_branch1_info_struct: List of layouts
!   layout_end       -- layout: ptc layout
!
! Output:
!   branch_ptc_info  -- ptc_branch1_info_struct: Updated list.
!-

subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   

implicit none

type (ptc_branch1_info_struct) branch_ptc_info, temp_info
type (layout), target :: layout_end
integer n

!

if (allocated(branch_ptc_info%layout)) then
  n = size(branch_ptc_info%layout)
  call move_alloc(branch_ptc_info%layout, temp_info%layout)
  allocate(branch_ptc_info%layout(n+1))
  branch_ptc_info%layout(1:n) = temp_info%layout
  branch_ptc_info%layout(n)%ptr => layout_end
  deallocate(temp_info%layout)

else
  allocate(branch_ptc_info%layout(1))
  branch_ptc_info%layout(1)%ptr => layout_end  
endif 

end subroutine add_ptc_layout_to_list

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine kill_ptc_layouts (lat)
!
! Routine to kill the layouts associated with a Bmad lattice.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input: 
!   lat  -- lat_struct: Bmad lattice with associated layouts.
!-

subroutine kill_ptc_layouts (lat)

use madx_ptc_module

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer ib, il

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do il = 1, size(branch%ptc%layout)
    call kill_layout_in_universe(branch%ptc%layout(il)%ptr)
  enddo
  deallocate(branch%ptc%layout)
enddo

end subroutine kill_ptc_layouts

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)
!
! Routine to calculate emittances, etc.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input: 
!   ele -- ele_struct: Element at which to evaluate the parameters.
!
! Output:
!   norm_mode       -- Normal_modes_struct
!     %e_loss
!     %a%tune, %b%tune, %z%tune
!     %a%alpha_damp, etc.
!     %a%emittance, etc.
!   sigma_map(6,6)  -- real(rp): Sigma matrix (Bmad coordinates).
!   closed_orb      -- coord_struct: Closed orbit at ele (Bmad coordinates).
!-

subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)

use madx_ptc_module

implicit none

type (ele_struct) ele
type (internal_state) state
type (normal_modes_struct) norm_mode
type (normal_spin) normal
type (damapspin) da_map
type (probe) x_probe
type (probe_8) x_probe8  
type (coord_struct) closed_orb

real(rp) sigma_mat(6,6)
real(dp) x(6), energy, deltap

!

check_krein = .false.

state = (default - nocavity0) + radiation0  ! Make sure have RF + radiation on.

x = 0
call find_orbit_x (x, state, 1.0d-5, fibre1 = ele%ptc_fibre%next)  ! find closed orbit
call vec_ptc_to_bmad (x, ele%ptc_fibre%next%beta0, closed_orb%vec)

call get_loss (ele%ptc_fibre%parent_layout, energy, deltap)
norm_mode%e_loss = 1d9 * energy
norm_mode%z%alpha_damp = deltap

call init (state, 1, 0)  ! First order DA
call alloc(normal)
call alloc(da_map)
call alloc(x_probe8)

normal%stochastic = .false. ! Normalization of the stochastic kick not needed.

x_probe = x
da_map = 1
x_probe8 = x_probe + da_map

! Remember: ptc calculates things referenced to the beginning of a fibre while
! Bmad references things at the exit end.

state = state+envelope0
call track_probe (x_probe8, state, fibre1 = ele%ptc_fibre%next)
da_map = x_probe8
normal = da_map

norm_mode%a%tune = normal%tune(1)   ! Fractional tune with damping
norm_mode%b%tune = normal%tune(2)
norm_mode%z%tune = normal%tune(3)

norm_mode%a%alpha_damp = normal%damping(1)
norm_mode%b%alpha_damp = normal%damping(2)
norm_mode%z%alpha_damp = normal%damping(3)

norm_mode%a%emittance = normal%emittance(1)
norm_mode%b%emittance = normal%emittance(2)
norm_mode%z%emittance = normal%emittance(3)

call sigma_mat_ptc_to_bmad (normal%s_ij0, ele%ptc_fibre%next%beta0, sigma_mat)

call kill(normal)
call kill(da_map)
call kill(x_probe8)

end subroutine ptc_emit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_closed_orbit_calc (branch, closed_orbit)
!
! Routine to calculate the closed orbit of a lattice branch using PTC.
! This routine will use the associate PTC layout if it exists and will
! create one if it does not.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   branch          -- branch_struct: Branch of a lattice.
!
! Output:
!   closed_orbit(:) -- coord_struct, allocatable: closed_orbit
!-

subroutine ptc_closed_orbit_calc (branch, closed_orbit)

use madx_ptc_module

implicit none

type (branch_struct), target :: branch
type (coord_struct), allocatable :: closed_orbit(:)
type (fibre), pointer :: fib

real(dp) x(6)
real(rp) vec(6)

integer i

!

x = 0
fib => branch%ele(0)%ptc_fibre%next
call find_orbit_x (x, default, 1.0d-5, fibre1 = fib)  ! find closed orbit
call vec_ptc_to_bmad (x, fib%beta0, vec)
call init_coord (closed_orbit(0), vec, branch%ele(0), .true., branch%param%particle)

do i = 1, branch%n_ele_track
  fib => branch%ele(i)%ptc_fibre%next
  call track_probe_x (x, default, branch%ele(i-1)%ptc_fibre%next, fib)
  call vec_ptc_to_bmad (x, fib%beta0, vec)
  call init_coord (closed_orbit(i), vec, branch%ele(i), .true., branch%param%particle)
enddo

end subroutine ptc_closed_orbit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_one_turn_map_at_ele (ele, rf_on, pz, map)
!
! Routine to calculate the one turn map for a ring.
! Note: Use set_ptc(no_cavity = True/False) set turn on/off the RF cavities.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele     -- ele_struct: Element determining start/end position for one turn map.
!   order   -- integer, optional: Order of the map. If not given then default order is used.
!
! Output:
!   map(6)  -- taylor_struct: Bmad taylor map
!-

subroutine ptc_one_turn_map_at_ele (ele, rf_on, pz, map)

use madx_ptc_module

implicit none

type (ele_struct), target :: ele
type (taylor_struct) map(6)
type (internal_state) state
type (fibre), pointer :: fib
type (damap) da_map
type (real_8) ray(6)

real(rp) pz
real(dp) x(6)

logical rf_on

!

if (rf_on) then
  state = default - nocavity0 
else
  state = default + nocavity0 
endif

! Find closed orbit

x = 0
x(5) = pz
fib => ele%ptc_fibre%next
call find_orbit_x (x, state, 1.0d-5, fibre1 = fib)  ! find closed orbit

! Construct map.

call alloc(da_map)
call alloc(ray)
da_map = 1   ! Identity
ray = da_map + x
call track_probe_x (ray, state, fibre1 = fib)

call real_8_to_taylor(ray, fib%beta0, map)

call kill(ray)
call kill(da_map)

end Subroutine ptc_one_turn_map_at_ele

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine write_ptc_flat_file_lattice (file_name, branch)
!
! Routine to create a PTC flat file lattice from a Bmad branch.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   file_name     -- character(*): Flat file name.
!   branch        -- branch_struct: Branch containing a layout.
!-

subroutine write_ptc_flat_file_lattice (file_name, branch)

use pointer_lattice

implicit none

type (branch_struct) branch
character(*) file_name

!

call print_complex_single_structure (branch%ptc%layout(1)%ptr, file_name)

end subroutine write_ptc_flat_file_lattice

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine modify_ptc_fibre (ele, param)
!
! Routine to modify an existing PTC fibre. 
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele           -- ele_struct: Element with corresponding PTC fibre.
!   param         -- lat_param_struct:
!
! Output:
!   ele%ptc_fibre 
!-

subroutine modify_ptc_fibre_attribute (ele, param)

use madx_ptc_module

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (keywords) ptc_key

real(rp) value

integer i

character(32), parameter :: r_name = 'modify_ptc_fibre_attribute'

!

call ele_to_an_bn (ele, param, ptc_key%list%k, ptc_key%list%ks, ptc_key%list%nmul)

do i = ptc_key%list%nmul, 1, -1
  call add (ele%ptc_fibre,  i, 0, ptc_key%list%k(i))
  call add (ele%ptc_fibre, -i, 0, ptc_key%list%ks(i))
enddo

end subroutine modify_ptc_fibre_attribute 

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
! Module needed:
!   use ptc_interface_mod
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
!   dx_tol_bend   -- real(rp): Tollerable residual orbit in a bend.
!   use_2nd_order -- logical, optional: If present and True then force the use of 2nd order
!                       integrator.
!
! Output:
!   ptc_layout -- layout: Lattice with the optimum number of tracking steps.
!-

subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max, &
                                    even_steps, r_typical, dx_tol_bend, use_2nd_order)

use madx_ptc_module

implicit none

type (layout) ptc_layout

real(rp) kl_max
real(rp), optional :: ds_max, dx_tol_bend, r_typical

integer :: limit_int(2)

logical, optional :: even_steps, use_2nd_order

!

resplit_cutting = 0   ! Ignore ds_max
if (present(ds_max)) then
  if (ds_max > 0) resplit_cutting = 2
endif

limit_int = [4, 18]
if (logic_option(.false., use_2nd_order)) limit_int = [10000, 10001] ! Something big.

call thin_lens_resplit (ptc_layout, kl_max, lim = limit_int, &
                        lmax0 = ds_max, sexr = r_typical, xbend = dx_tol_bend)

end subroutine ptc_calculate_tracking_step_size

end module
