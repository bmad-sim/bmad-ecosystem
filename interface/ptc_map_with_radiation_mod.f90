module ptc_map_with_radiation_mod

! Etienne wanted the "zhe" stuff to be standalone and so duplicated structures in 

use ptc_layout_mod
use duan_zhe_map, only: tree_element_zhe => tree_element, probe_zhe => probe, track_tree_probe_complex_zhe, &
                        zhe_ini

type ptc_map_with_rad_struct
  type (tree_element_zhe) sub_map(3)    ! Type tree_element in PTC
  character(200) lattice_file     ! Name of the lattice file
  integer map_order
  logical radiation_damping_on
  integer ix_branch
  integer ix_ele_start            ! Start point for making the map
  integer ix_ele_end              ! End point for making the map
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_setup_map_with_radiation (map_with_rad, ele1, ele2, map_order, include_damping, orbit1, err_flag)
!
! Routine to construct a map including radiation damping and excitation.
! Note: The setting of bmad_com%radiation_damping_on will determine if damping is included in the map.
!
! ele1/ele2 must have an associated PTC layout (which can be constructed by calling lat_to_ptc_layout).
!
! To track after calling this routine track by calling ptc_track_with_radiation.
! To cleanup memory after using, call ptc_kill_map_with_radiation.
! To save a map call ptc_write_map_with_radiation.
! To read a saved map call ptc_read_map_with_radiation.
! To set the random number seed call: ptc_ran_seed_put.
!
! Input:
!   ele1            -- ele_struct: The map starts at the exit end of ele1.
!   ele2            -- ele_struct, optional: The map ends at the exit end of ele2. If not present, the 
!                        1-turn map will be constructed.
!   map_order       -- integer, optional: Order of the map. 
!                        If not present or less than 1, the currently set order is used.
!   include_damping -- logical, optional: If True (default), the map will be constructed with radiation damping included.
!                        If False, the map will not be constructed with radiation dampling included. 
!                        Since radiation damping can always be turned off when tracking, if you are only concerned about
!                        the orbital motion, there is no reason to create a map without damping. 
!                        However, the spin map is constructed about the closed orbit so the spin map will be affected
!                        by whether damping is on or not. 
!                        To the extent that the damping is small the shift in the spin map will be small.
!   orbit1          -- coord_struct, optional: Orbit at ele1 about which the map is constructed.
!                        If not present then the orbit will be computed using PTC tracking.
!
! Output:
!   map_with_rad    -- ptc_map_with_rad_struct: Transport map.
!   err_flag        -- logical, optional: Set True if there is an error such as not associated PTC layout.
!-

subroutine ptc_setup_map_with_radiation (map_with_rad, ele1, ele2, map_order, include_damping, orbit1, err_flag)

use pointer_lattice

implicit none

type (ptc_map_with_rad_struct) map_with_rad
type (ele_struct) ele1
type (ele_struct), optional :: ele2
type (coord_struct), optional :: orbit1
type (layout), pointer :: ptc_layout
type (internal_state) state
type (branch_struct), pointer :: branch
type (fibre), pointer :: f1, f2
type (tree_element) tree_map(3)
type (c_damap) c_map1, c_ident
type (probe_8) pb8
type (probe) pb

real(rp) orb(6), orb0(6), e_ij(6,6)

integer, optional :: map_order
integer val_save

logical, optional :: err_flag, include_damping

character(*), parameter :: r_name = 'ptc_setup_map_with_radiation'

!

if (present(err_flag)) err_flag = .true.

call zhe_ini(bmad_com%spin_tracking_on)

if (logic_option(.true., include_damping)) then
  state = ptc_com%base_state + radiation0 + envelope0
  map_with_rad%radiation_damping_on = .true.
else
  state = ptc_com%base_state + envelope0
  map_with_rad%radiation_damping_on = .false.
endif

if (bmad_com%spin_tracking_on) state = state + spin0
if (.not. rf_is_on(ele1%branch)) state = state + nocavity0

map_with_rad%map_order = integer_option(ptc_com%taylor_order_ptc, map_order)
if (map_with_rad%map_order < 1) map_with_rad%map_order = ptc_com%taylor_order_ptc

call init_all(state, map_with_rad%map_order, 0)

branch => pointer_to_branch(ele1)
ptc_layout => branch%ptc%m_t_layout
map_with_rad%ix_branch = branch%ix_branch
map_with_rad%lattice_file = branch%lat%input_file_name

if (.not. associated(ptc_layout)) then
  call out_io (s_fatal$, r_name, 'NO ASSOCIATED PTC LAYOUT PRESENT!')
  if (global_com%exit_on_error) call err_exit
  return
endif

f1 => pointer_to_fibre(ele1)
map_with_rad%ix_ele_start = ele1%ix_ele

if (present(ele2)) then
  f2 => pointer_to_fibre(ele2)
  map_with_rad%ix_ele_end = ele2%ix_ele
else
  f2 => f1
  map_with_rad%ix_ele_end = ele1%ix_ele
endif

if (present(orbit1)) then
  orb = orbit1%vec
else
  orb = 0
  call find_orbit_x(orb, STATE, 1.0d-8, fibre1 = f1)
  if (.not. check_stable) then
    call out_io (s_error$, r_name, 'CANNOT FIND CLOSED ORBIT WHEN TRCKING WITH RADIATION IN PTC!')
    return
  endif
endif

call set_ptc_quiet(0, set$, val_save)

!!call radia_full(ptc_layout, f1 = f1, estate = state0, e_ij=e_ij, ngen=10, bunch_zhe=bunch_zhe)

call alloc(pb8)
call alloc(c_map1, c_ident)
pb = orb
c_ident = 1
pb8 = pb + c_ident

! If f1 == f2 then want one turn map instead of unit map.

if (associated(f1, f2)) then
  call propagate(pb8, state, fibre1 = f1)
else
  call propagate(pb8, state, fibre1 = f1, fibre2 = f2)
endif

c_map1 = pb8

call fill_tree_element_line_zhe_outside_map(c_map1, as_is=.false., stochprec=1.d-10, tree_zhe=map_with_rad%sub_map) 

call set_ptc_quiet(0, unset$, val_save)
call kill (pb8)
call kill (c_map1, c_ident)

if (present(err_flag)) err_flag = .false.

end subroutine ptc_setup_map_with_radiation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_track_map_with_radiation (orbit, map_with_rad, rad_damp, rad_fluct)
!
! Routine to track through a map that includes radiation.
! To construct the map, use the routine ptc_setup_map_with_radiation.
! To cleanup memory after using, call ptc_kill_map_with_radiation.
! To save a map call ptc_write_map_with_radiation.
! To read a saved map call ptc_read_map_with_radiation.
! To set the random number seed call: ptc_ran_seed_put.
!
! Input:
!   orbit            -- coord_struct: Starting orbit.
!   map_with_rad     -- ptc_map_with_rad_struct: Map with radiation included.
!   rad_damp         -- logical, optional: Override the setting of bmad_com%radiation_damping_on
!   rad_fluct        -- logical, optional: Override the setting of bmad_com%radiation_fluctuations_on
!   
! Output:
!   orbit            -- coord_struct: Ending orbit after tracking through the map.
!     %state            -- Set to lost$ if there is a problem.
!-

subroutine ptc_track_map_with_radiation (orbit, map_with_rad, rad_damp, rad_fluct)

use rotation_3d_mod
use duan_zhe_map, only: assignment(=), C_VERBOSE_ZHE

implicit none

type (coord_struct) orbit
type (ptc_map_with_rad_struct) map_with_rad
type (probe_zhe) z_probe

logical, optional :: rad_damp, rad_fluct
logical damp, fluct

!

damp   = logic_option(bmad_com%radiation_damping_on, rad_damp)
fluct = logic_option(bmad_com%radiation_fluctuations_on, rad_fluct)
C_VERBOSE_ZHE = .false.

!

z_probe = orbit%vec
z_probe%q%x = [1, 0, 0, 0]

call track_tree_probe_complex_zhe (map_with_rad%sub_map, z_probe, bmad_com%spin_tracking_on, damp, fluct)
if (z_probe%u) orbit%state = lost$   ! %u = T => "unstable".

orbit%vec = z_probe%x
if (bmad_com%spin_tracking_on) then
  orbit%spin = rotate_vec_given_quat(z_probe%q%x, orbit%spin)
endif

end subroutine ptc_track_map_with_radiation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_write_map_with_radiation(map_with_rad, file_name, file_unit)
!
! Routine to create or append to a binary file containing a ptc_map_with_rad_struct map.
!
! Either file_name or file_unit must be present but not both.
! If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
! and to close the file afterwards.
!
! Input:
!   map_with_rad    -- ptc_map_with_rad_struct: Map with radiation included.
!   file_name       -- character(*), optional: Name of binary file to create.
!   file_unit       -- integer, optional: File unit number to append to.
!-

subroutine ptc_write_map_with_radiation(map_with_rad, file_name, file_unit)

type (ptc_map_with_rad_struct), target :: map_with_rad
type (tree_element_zhe), pointer :: t

integer i, j, k, iu
integer, optional :: file_unit
character(*), optional :: file_name

!

if (present(file_name)) then
  iu = lunget()
  open (iu, file = file_name, form = 'unformatted')
else
  iu = file_unit
endif

write (iu) map_with_rad%lattice_file
write (iu) map_with_rad%map_order, map_with_rad%radiation_damping_on, &
            map_with_rad%ix_branch, map_with_rad%ix_ele_start, map_with_rad%ix_ele_end

do k = 1, 3
  t => map_with_rad%sub_map(k)
  call write_real1(t%cc)
  call write_real1(t%fixr)
  call write_real1(t%fix)
  call write_real1(t%fix0)
  call write_int1(t%jl)
  call write_int1(t%jv)
  call write_int0(t%n)
  call write_int0(t%np)
  call write_int0(t%no)
  call write_real2(t%e_ij)
  call write_real2(t%rad)
  call write_real0(t%ds)
  call write_real0(t%beta0)
  call write_real0(t%eps)
  call write_logic0(t%symptrack)
  call write_logic0(t%usenonsymp)
  call write_logic0(t%factored)
enddo

if (present(file_name)) close(iu)

!-------------------------------------------------------------------
contains

subroutine write_real0(rr)
real(rp), pointer :: rr
write (iu) rr
end subroutine write_real0

subroutine write_real1(rr)
real(rp), pointer :: rr(:)
write (iu) size(rr)
write (iu) rr
end subroutine write_real1

subroutine write_real2(rr)
real(rp), pointer :: rr(:,:)
write (iu) size(rr, 1), size(rr, 2)
write (iu) rr
end subroutine write_real2

subroutine write_int0(rr)
integer, pointer :: rr
write (iu) rr
end subroutine write_int0

subroutine write_int1(rr)
integer, pointer :: rr(:)
write (iu) size(rr)
write (iu) rr
end subroutine write_int1

subroutine write_logic0(rr)
logical, pointer :: rr
write (iu) rr
end subroutine write_logic0

end subroutine ptc_write_map_with_radiation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_read_map_with_radiation(map_with_rad, err_flag, file_name, file_unit)
!
! Routine to read a binary file containing a ptc_map_with_rad_struct map
!
! Either file_name or file_unit must be present but not both.
! File_unit is used when there are multiple maps in a file.
! If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
! and to close the file afterwards.
!
! Input:
!   file_name       -- character(*), optional: Name of binary file.
!   file_unit       -- integer, optional: File unit number read from.
!
! Output:
!   map_with_rad    -- ptc_map_with_rad_struct: Map with radiation included.
!   err_flag        -- Logical: Set True if there is a read error.
!-

subroutine ptc_read_map_with_radiation(map_with_rad, err_flag, file_name, file_unit)

type (ptc_map_with_rad_struct), target :: map_with_rad
type (tree_element_zhe), pointer :: t

integer i, j, k, iu
integer, optional :: file_unit
logical err_flag
character(*), optional :: file_name

!

err_flag = .true.

if (present(file_name)) then
  iu = lunget()
  open (iu, file = file_name, form = 'unformatted', status = 'old', err = 9100)
else
  iu = file_unit
endif

read (iu, err = 9000, end = 9000) map_with_rad%lattice_file
read (iu, err = 9000, end = 9000) map_with_rad%map_order, map_with_rad%radiation_damping_on, &
            map_with_rad%ix_branch, map_with_rad%ix_ele_start, map_with_rad%ix_ele_end

do k = 1, 3
  t => map_with_rad%sub_map(k)
  if (.not. read_real1(t%cc)) goto 9000
  if (.not. read_real1(t%fixr)) goto 9000
  if (.not. read_real1(t%fix)) goto 9000
  if (.not. read_real1(t%fix0)) goto 9000
  if (.not. read_int1(t%jl)) goto 9000
  if (.not. read_int1(t%jv)) goto 9000
  if (.not. read_int0(t%n)) goto 9000
  if (.not. read_int0(t%np)) goto 9000
  if (.not. read_int0(t%no)) goto 9000
  if (.not. read_real2(t%e_ij)) goto 9000
  if (.not. read_real2(t%rad)) goto 9000
  if (.not. read_real0(t%ds)) goto 9000
  if (.not. read_real0(t%beta0)) goto 9000
  if (.not. read_real0(t%eps)) goto 9000
  if (.not. read_logic0(t%symptrack)) goto 9000
  if (.not. read_logic0(t%usenonsymp)) goto 9000
  if (.not. read_logic0(t%factored)) goto 9000
enddo

if (present(file_name)) close(iu)

call zhe_ini(bmad_com%spin_tracking_on)

err_flag = .false.

!

9000 continue
if (present(file_name)) close(iu)
9100 continue

!-------------------------------------------------------------------
contains

function read_real0(rr) result (ok)
real(rp), pointer :: rr
integer ios
logical ok
allocate(rr)
read (iu, iostat = ios) rr
ok = (ios == 0)
end function read_real0

function read_real1(rr) result (ok)
real(rp), pointer :: rr(:)
logical ok
integer n, ios1, ios2
read (iu, iostat = ios1) n
allocate(rr(n))
read (iu, iostat = ios2) rr
ok = (ios1 == 0 .and. ios2 == 0)
end function read_real1

function read_real2(rr) result (ok)
real(rp), pointer :: rr(:,:)
logical ok
integer n1, n2, ios1, ios2
read (iu, iostat = ios1) n1, n2
allocate(rr(n1,n2))
read (iu, iostat = ios2) rr
ok = (ios1 == 0 .and. ios2 == 0)
end function read_real2

function read_int0(rr) result (ok)
integer, pointer :: rr
integer ios
logical ok
allocate(rr)
read (iu, iostat = ios) rr
ok = (ios == 0)
end function read_int0

function read_int1(rr) result (ok)
integer, pointer :: rr(:)
logical ok
integer n, ios1, ios2
read (iu, iostat = ios1) n
allocate(rr(n))
read (iu, iostat = ios2) rr
ok = (ios1 == 0 .and. ios2 == 0)
end function read_int1

function read_logic0(rr) result (ok)
integer ios
logical ok
logical, pointer :: rr
allocate(rr)
read (iu, iostat = ios) rr
ok = (ios == 0)
end function read_logic0

end subroutine ptc_read_map_with_radiation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ptc_kill_map_with_radiation(map_with_rad)
!
! Routine to kill a binary file containing a ptc_map_with_rad_struct map
!
! Input:
!   map_with_rad     -- ptc_map_with_rad_struct: Map with radiation included.
!
! Output:
!   map_with_rad     -- ptc_map_with_rad_struct: Deallocated map.
!-

subroutine ptc_kill_map_with_radiation(map_with_rad)

type (ptc_map_with_rad_struct), target :: map_with_rad
type (tree_element_zhe), pointer :: t
integer k

!

do k = 1, 3
  t => map_with_rad%sub_map(k)
  deallocate(t%cc)
  deallocate(t%fixr)
  deallocate(t%fix)
  deallocate(t%fix0)
  deallocate(t%jl)
  deallocate(t%jv)
  deallocate(t%n)
  deallocate(t%np)
  deallocate(t%no)
  deallocate(t%e_ij)
  deallocate(t%rad)
  deallocate(t%ds)
  deallocate(t%beta0)
  deallocate(t%eps)
  deallocate(t%symptrack)
  deallocate(t%usenonsymp)
  deallocate(t%factored)
enddo

end subroutine ptc_kill_map_with_radiation

end module
