module ptc_track_with_radiation_mod

! Etienne wanted the "zhe" stuff to be standalone and so duplicated structures in 

use ptc_layout_mod
use duan_zhe_map, only: tree_element_zhe => tree_element, probe_zhe => probe, track_tree_probe_complex_zhe, zhe_ini

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_setup_map_with_radiation (map_with_rad, ele1, ele2, orbit1, map_order, err_flag)
!
! Routine to construct a map including radiation damping and excitation.
!
! ele1/ele2 must have an associated PTC layout (which can be constructed by calling lat_to_ptc_layout).
!
! To track after calling this routine track by calling ptc_track_with_radiation.
!
! Input:
!   ele1            -- ele_struct: The map starts at the exit end of ele1.
!   ele2            -- ele_struct, optional: The map ends at the exit end of ele2. If not present, the 
!                       1-turn map will be constructed.
!   orbit1          -- coord_struct, optional: Orbit at ele1 about which the map is constructed.
!   map_order       -- integer, optional: Order of the map. If not set the currently set order is used.
!                         
!
! Output:
!   map_with_rad(3) -- tree_element_zhe: Transport map.
!   err_flag        -- logical, optional: Set True if there is an error such as not associated PTC layout.
!-

subroutine ptc_setup_map_with_radiation (map_with_rad, ele1, ele2, orbit1, map_order, err_flag)

use pointer_lattice

implicit none

type (tree_element_zhe) map_with_rad(3)
type (ele_struct) ele1
type (ele_struct), optional :: ele2
type (coord_struct), optional :: orbit1
type (layout), pointer :: ptc_layout
type (internal_state) state
type (branch_struct), pointer :: branch
type (fibre), pointer :: f1, f2
type (tree_element) tree_map(3)

real(rp) orb(6)

integer, optional :: map_order
integer order, val_save

logical, optional :: err_flag

character(*), parameter :: r_name = 'ptc_setup_map_with_radiation'

!

if (present(err_flag)) err_flag = .true.

call zhe_ini(bmad_com%spin_tracking_on)
use_bmad_units = .true.

state = default0 + radiation0 + envelope
if (.not. rf_is_on(ele1%branch)) state = state + nocavity0

order = integer_option(ptc_com%taylor_order_ptc, map_order)
call init_all(state, order, 0)

branch => pointer_to_branch(ele1)
ptc_layout => branch%ptc%m_t_layout

if (.not. associated(ptc_layout)) then
  call out_io (s_fatal$, r_name, 'NO ASSOCIATED PTC LAYOUT PRESENT!')
  if (global_com%exit_on_error) call err_exit
  return
endif

f1 => pointer_to_ptc_ref_fibre(ele1)
if (present(ele2)) then
  f2 => pointer_to_ptc_ref_fibre(ele2)
else
  f2 => f1
endif

if (present(orbit1)) then
  orb = orbit1%vec
else
  orb = 0
  call find_orbit_x(orb, STATE, 1.0d-8, fibre1 = f1)
endif

call set_ptc_quiet(0, set$, val_save)
call fill_tree_element_line_zhe(state, f1, f2, order, orb, stochprec = 1d-10, sagan_tree = tree_map)
call set_ptc_quiet(0, unset$, val_save)

call copy_this_tree (tree_map, map_with_rad)

use_bmad_units = .false. ! Since Zhe stuff is standalone this will not affect the use of map_with_rad.
if (present(err_flag)) err_flag = .false.

!----------------------------------------------------------------------------------
contains

subroutine copy_this_tree(t,u)

implicit none
type(tree_element) :: t(3)
type(tree_element_zhe) :: u(3)
integer i

do i = 1, 3
  u(i)%cc         =>t(i)%cc
  u(i)%jl         =>t(i)%jl
  u(i)%jv         =>t(i)%jv
  u(i)%n          =>t(i)%n
  u(i)%np         =>t(i)%np
  u(i)%no         =>t(i)%no
  u(i)%fixr       =>t(i)%fixr
  u(i)%ds         =>t(i)%ds
  u(i)%beta0      =>t(i)%beta0
  u(i)%fix        =>t(i)%fix
  u(i)%fix0       =>t(i)%fix0
  u(i)%e_ij       =>t(i)%e_ij
  u(i)%rad        =>t(i)%rad
  u(i)%eps        =>t(i)%eps
  u(i)%symptrack  =>t(i)%symptrack
  u(i)%usenonsymp =>t(i)%usenonsymp
  u(i)%factored   =>t(i)%factored
enddo

end subroutine copy_this_tree

end subroutine ptc_setup_map_with_radiation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_track_with_radiation (orbit, map_with_rad, rad_damp, rad_fluct)
!
! Routine to track through a map that includes radiation.
! Use the routine ptc_setup_map_with_radiation to construct the map.
!
! Input:
!   orbit         -- coord_struct: Starting orbit.
!   map_with_rad  -- tree_element_zhe: Map with radiation included.
!   rad_damp      -- logical, optional: Override the setting of bmad_com%radiation_damping_on
!   rad_fluct     -- logical, optional: Override the setting of bmad_com%radiation_fluctuations_on
!   
! Output:
!   orbit         -- coord_struct: Ending orbit after tracking through the map.
!     %state        -- Set to lost$ if there is a problem.
!-

subroutine ptc_track_with_radiation (orbit, map_with_rad, rad_damp, rad_fluct)

use rotation_3d_mod
use duan_zhe_map, only: assignment(=)

implicit none

type (coord_struct) orbit
type (tree_element_zhe) map_with_rad(3)
type (probe_zhe) z_probe

logical, optional :: rad_damp, rad_fluct
logical damp, fluct

!

damp   = logic_option(bmad_com%radiation_damping_on, rad_damp)
fluct = logic_option(bmad_com%radiation_fluctuations_on, rad_fluct)

!

z_probe = orbit%vec
z_probe%q%x = [1, 0, 0, 0]

call track_tree_probe_complex_zhe (map_with_rad, z_probe, bmad_com%spin_tracking_on, rad_damp, rad_fluct)
if (z_probe%u) orbit%state = lost$   ! %u = T => "unstable".

orbit%vec = z_probe%x
if (bmad_com%spin_tracking_on) then
  orbit%spin = rotate_vec_given_quat(z_probe%q%x, orbit%spin)
endif

end subroutine ptc_track_with_radiation

end module
