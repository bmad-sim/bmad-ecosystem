!+
! Module ptc_layout_mod
!
! Module of PTC layout interface routines.
! Also see: ptc_interface_mod
!-

module ptc_layout_mod

use ptc_interface_mod

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
! Subroutine to create a PTC layout from a Bmad lat.
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

use s_fibre_bundle, only: ring_l, append, lp, layout, fibre
use mad_like, only: set_up, kill, lielib_print
use madx_ptc_module, only: m_u, append_empty_layout, survey, make_node_layout

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct) drift_ele
type (ele_struct), pointer :: ele

integer n, ib, ie
logical doneit

! transfer elements.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  call append_empty_layout(m_u)
  call set_up(m_u%end)

  allocate(branch%ptc%layout(1))
  branch%ptc%layout(1)%ptr => m_u%end   ! Save layout

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (tracking_uses_hard_edge_model(ele)) then
      call create_hard_edge_drift (ele, entrance_end$, drift_ele)
      call ele_to_fibre (drift_ele, drift_ele%ptc_fibre, branch%param%particle, .true., for_layout = .true.)
    endif

    call ele_to_fibre (ele, ele%ptc_fibre, branch%param%particle, .true., for_layout = .true.)

    if (tracking_uses_hard_edge_model(ele)) then
      call create_hard_edge_drift (ele, exit_end$, drift_ele)
      ! ele%ptc_fibre points to last PTC fibre.
      call ele_to_fibre (drift_ele, ele%ptc_fibre, branch%param%particle, .true., for_layout = .true.)
    endif
  enddo

  ! End stuff

  if (branch%param%lattice_type == circular_lattice$) then
    m_u%end%closed = .true.
  else
    m_u%end%closed = .false.
  endif


  n = lielib_print(12)
  lielib_print(12) = 0  ! No printing info messages

  doneit = .true.
  call ring_l (m_u%end, doneit)
  call survey (m_u%end)
  call make_node_layout (m_u%end)

  lielib_print(12) = n

  branch%ele(0)%ptc_fibre => branch%ele(1)%ptc_fibre%previous

enddo

end subroutine lat_to_ptc_layout

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
type (layout), pointer :: ptc_layout
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

ptc_layout => ele%branch%ptc%layout(1)%ptr

state = (default - nocavity0) + radiation0  ! Set state flags

x = 0
call find_orbit_x (x, state, 1.0d-5, fibre1 = ele%ptc_fibre%next)  ! find_orbit == find closed orbit
call vec_ptc_to_bmad (x, ele%ptc_fibre%next%beta0, closed_orb%vec)

call get_loss (ptc_layout, energy, deltap)
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
! Subroutine one_turn_map_at_ele (ele, map, order)
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

subroutine one_turn_map_at_ele (ele, map, order)

implicit none

type (ele_struct), target :: ele
type (taylor_struct) map(6)


integer, optional :: order

!

!ptc_layout => ele%branch%ptc%layout(1)%ptr

!x = 0
!call find_orbit_x (x, default, 1.0d-5, fibre1 = ele%ptc_fibre%next)  ! find_orbit == find closed orbit
!closed_orb%vec = x

end Subroutine one_turn_map_at_ele

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
! Subroutine modify_ptc_fibre (ele)
!
! Routine to modify an existing PTC fibre. 
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele           -- ele_struct: Element with corresponding PTC fibre.
!
! Output:
!   ele%ptc_fibre 
!-

subroutine modify_ptc_fibre_attribute (ele, attribute, value)

implicit none

type (ele_struct), target :: ele

real(rp) value

character(*) attribute
character(32), parameter :: r_name = 'modify_ptc_fibre_attribute'

!

select case (ele%key)

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN ELEMENT TYPE: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit
end select

!

end subroutine modify_ptc_fibre_attribute 

end module
