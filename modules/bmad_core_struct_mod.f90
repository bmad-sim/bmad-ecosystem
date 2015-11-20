!+
! Module bmad_core_struct_mod
!
! Collection of routines for initializing, allocating, and deallocating bmad structures.
! Also included are routines that set struct1 = struct12 without using 
! the overloaded equal sign.
!
! NOTE: NO ROUTINES IN THIS MODULE HAVE ACCESS TO THE OVERLOADED
! EQUAL SIGN USED TO SET ELE1 = ELE2, LAT1 = LAT2 ETC.
!-

module bmad_core_struct_mod

use basic_bmad_mod
use basic_bmad_interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord (...)
!
! Routine to initialize a coord_struct. 
!
! This routine is an overloaded name for:
!   Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)
!   Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)
!
! Note: Unless shift_vec6 is set to False, if ele is an init_ele or e_gun, orb%vec(6) is shifted
! so that the particle's energy is maintained at ele%value(p0c_start$).
!
! Note: For a photon, orb%vec(5) is set depending upon where the photon is relative to the element.
! Note: If the particle is initialized with element_end = inside$, orb%s will not be set.
! Note: If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle 
!       except if direction is set.
!
! Modules needed:
!   use bmad
!
! Input:
!   orb_in       -- Coord_struct: Input orbit.
!   vec(6)       -- real(rp), optional: Coordinate vector. If not present then taken to be zero.
!   ele          -- ele_struct, optional: Particle is initialized to start from the entrance end of ele
!   element_end  -- integer, optional: upstream_end$, downstream_end$, inside$, or start_end$.
!                     Must be present if ele argument is present.
!                     start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1.
!                     Default is upstream_end$.
!   particle     -- Integer, optional: Particle type (electron$, etc.). 
!                     If particle = not_set$ and orb_in is present, use orb_in%species instead.
!   dirction     -- Integer, optional: +1 -> moving downstream +s direciton, -1 -> moving upstream.
!                     0 -> Ignore. Default is to not change orb%direction except for photons which get set
!                     according to orb%vec(6).
!   E_photon     -- real(rp), optional: Photon energy if particle is a photon. Ignored otherwise.
!   t_ref_offset -- real(rp), optional: Offset of the reference time. This is non-zero when
!                     there are multiple bunches and the reference time for a particular particle
!                     is pegged to the time of the center of the bunch.
!   shift_vec6   -- Logical, optional: If present and False, prevent the shift of orb%vec(6).
!
! Output:
!   orb -- Coord_struct: Initialized coordinate.
!                 Note: For photons, orb%vec(6) is computed as sqrt(1 - vec(2)^2 - vec(4)^2) if needed.
!-

interface init_coord
  module procedure init_coord1
  module procedure init_coord2
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine reallocate_coord (...)
!
! Routine to allocate or reallocate at allocatable coord_struct array.
! reallocate_coord is an overloaded name for:
!   reallocate_coord_n (coord, n_coord)
!   reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable coord_struct array to at least:
!     coord(0:n_coord)                            if n_coord arg is used.
!     coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.
!
! The old coordinates are saved
! If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
! In any case, coord(n)%vec for n > 0 is set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:)  -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!   lat       -- lat_struct: Lattice 
!   ix_branch -- Integer, optional: Branch to use. Default is 0 (main branch).
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

interface reallocate_coord
  module procedure reallocate_coord_n
  module procedure reallocate_coord_lat
end interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine transfer_twiss (ele_in, ele_out)
!
! Routine to transfer the twiss parameters from one element to another.
!
! Moduels needed:
!   use bmad
!
! Input:
!   ele_in   -- Ele_struct: Element with existing Twiss parameters.
!
! Output:
!   ele_out  -- Ele_struct: Element receiving the Twiss parameters.
!-

subroutine transfer_twiss (ele_in, ele_out)

implicit none

type (ele_struct) ele_in, ele_out

!

ele_out%x       = ele_in%x
ele_out%y       = ele_in%y
ele_out%a       = ele_in%a
ele_out%b       = ele_in%b
ele_out%z       = ele_in%z
ele_out%c_mat   = ele_in%c_mat
ele_out%gamma_c = ele_in%gamma_c

end subroutine transfer_twiss

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_lat_parameters (lat_in, lat_out)
!
! Subroutine to transfer the lat parameters (such as lat%name, lat%param, etc.)
! from one lat to another. The only stuff that is not transfered are things
! that are (or have) pointers or arrays
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat with parameters set.
!-

subroutine transfer_lat_parameters (lat_in, lat_out)

implicit none

type (lat_struct), intent(in) :: lat_in
type (lat_struct) :: lat_out

!

lat_out%use_name                  = lat_in%use_name
lat_out%lattice                   = lat_in%lattice
lat_out%input_file_name           = lat_in%input_file_name
lat_out%title                     = lat_in%title
lat_out%a                         = lat_in%a
lat_out%b                         = lat_in%b
lat_out%z                         = lat_in%z
lat_out%param                     = lat_in%param
lat_out%lord_state                = lat_in%lord_state
lat_out%beam_start                = lat_in%beam_start
lat_out%pre_tracker               = lat_in%pre_tracker
lat_out%version                   = lat_in%version
lat_out%n_ele_track               = lat_in%n_ele_track
lat_out%n_ele_max                 = lat_in%n_ele_max
lat_out%n_control_max             = lat_in%n_control_max
lat_out%n_ic_max                  = lat_in%n_ic_max
lat_out%input_taylor_order        = lat_in%input_taylor_order
lat_out%absolute_time_tracking    = lat_in%absolute_time_tracking
lat_out%ptc_uses_hard_edge_drifts = lat_in%ptc_uses_hard_edge_drifts
lat_out%surface                   => lat_in%surface  

end subroutine transfer_lat_parameters

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
!
! Subroutine to transfer a Taylor map from one element to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in       -- Ele_struct: Element with the Taylor map.
!   taylor_order -- Integer, optional: Order to truncate the Taylor map at.
!
! Output:
!   ele_out      -- Ele_struct: Element receiving the Taylor map truncated to
!                     order taylor_order.
!-

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)

implicit none

type (ele_struct) ele_in, ele_out
integer, optional :: taylor_order
integer it, ix, k 

!

do it = 1, 6

  if (present(taylor_order)) then
    ix = 0
    do k = 1, size(ele_in%taylor(it)%term)
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
      ix = ix + 1
    enddo
  else
    ix = size(ele_in%taylor(it)%term)
  endif

  if (.not. associated(ele_out%taylor(it)%term)) allocate (ele_out%taylor(it)%term(ix))
  if (size(ele_out%taylor(it)%term) /= ix) allocate (ele_out%taylor(it)%term(ix))

  ix = 0
  do k = 1, size(ele_in%taylor(it)%term)
    if (present(taylor_order)) then
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
    endif
    ix = ix + 1
    ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
  enddo

enddo

ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

end subroutine transfer_ele_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ele (ele1, ele2, nullify_pointers)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal operator.
! The result is that ele2's pointers will point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1             -- Ele_struct:
!   nullify_pointers -- Logical, optional: If present and True then nullify the 
!                         pointers in ele2 except for the ele2%lat and ele2%lord pointers. 
!                         This gives a "bare bones" copy where one does not have to 
!                         worry about deallocating allocated structure components later.
!
! Output:
!   ele2 -- Ele_struct:
!-

subroutine transfer_ele (ele1, ele2, nullify_pointers)

type (ele_struct), target :: ele1
type (ele_struct) :: ele2
logical, optional :: nullify_pointers

!

ele2 = ele1

if (logic_option (.false., nullify_pointers)) then
  call deallocate_ele_pointers (ele2, .true.)
  ele2%branch => ele1%branch  ! Reinstate
  ele2%lord   => ele1%lord    ! Reinstate
endif

end subroutine transfer_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_eles (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1(:) -- Ele_struct:
!
! Output:
!   ele2(:) -- Ele_struct:
!-

subroutine transfer_eles (ele1, ele2)

type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(inout) :: ele2(:)

ele2 = ele1

end subroutine transfer_eles

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branch (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1 -- Branch_struct:
!
! Output:
!   branch2 -- Branch_struct:
!-

subroutine transfer_branch (branch1, branch2)

type (branch_struct) :: branch1
type (branch_struct) :: branch2

!

branch2 = branch1

end subroutine transfer_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branches (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1(:) -- Branch_struct:
!
! Output:
!   branch2(:) -- Branch_struct:
!-

subroutine transfer_branches (branch1, branch2)

type (branch_struct) :: branch1(:)
type (branch_struct) :: branch2(:)

branch2 = branch1

end subroutine transfer_branches

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_lat (lat1, lat2)
!
! Subroutine to set lat2 = lat1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end lat2's pointers point to the same memory as lat1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   lat1 -- lat_struct:
!
! Output:
!   lat2 -- lat_struct:
!-

subroutine transfer_lat (lat1, lat2)

type (lat_struct), intent(in) :: lat1
type (lat_struct), intent(out) :: lat2

lat2 = lat1

end subroutine transfer_lat

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wig (wig_in, wig_out)
!
! Subroutine to point wig_out => wig_in
!
! Modules needed:
!   use bmad
!
! Input:
!   wig_in  -- Wig_struct, pointer: Input wiggler field.
!
! Output:
!   wig_out -- Wig_struct, pointer: Output wiggler field.
!-

subroutine transfer_wig (wig_in, wig_out)

implicit none

type (wig_struct), pointer :: wig_in, wig_out

!

if (.not. associated(wig_in) .and. .not. associated(wig_out)) return
if (associated(wig_in, wig_out)) return

! If both associated must be pointing to different memory locations

if (associated(wig_in) .and. associated(wig_out)) then
  wig_out%n_link = wig_out%n_link - 1
  if (wig_out%n_link == 0) then
    deallocate (wig_out%term)
    deallocate (wig_out)
  endif
  wig_out => wig_in
  wig_out%n_link = wig_out%n_link + 1

elseif (associated(wig_out)) then 
  wig_out%n_link = wig_out%n_link - 1
  if (wig_out%n_link == 0) then
    deallocate (wig_out%term)
    deallocate (wig_out)
  else
    nullify (wig_out)
  endif

elseif (associated(wig_in)) then 
  wig_out => wig_in
  wig_out%n_link = wig_out%n_link + 1
endif

end subroutine transfer_wig

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wall3d (wall3d_in, wall3d_out)
!
! Subroutine to point wall3d_out => wall3d_in
!
! Modules needed:
!   use bmad
!
! Input:
!   wall3d_in(:)  -- Wall3d_struct, pointer: Input wall3dgler field.
!
! Output:
!   wall3d_out(:) -- Wall3d_struct, pointer: Output wall3dgler field.
!-

subroutine transfer_wall3d (wall3d_in, wall3d_out)

implicit none

type (wall3d_struct), pointer :: wall3d_in(:), wall3d_out(:)

!

if (.not. associated(wall3d_in) .and. .not. associated(wall3d_out)) return
if (associated(wall3d_in, wall3d_out)) return

! If both associated must be pointing to different memory locations

if (associated(wall3d_out)) call deallocate_wall3d_pointer(wall3d_out)

if (associated(wall3d_in)) then 
  wall3d_out => wall3d_in
  wall3d_out%n_link = wall3d_out%n_link + 1
endif

end subroutine transfer_wall3d

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_em_field (field_in, field_out)
!
! Subroutine to transfer the field info from one struct to another.
! In the end will have:
!     field_out%map  => field_in%map
!     field_out%grid => field_in%grid
!
! Modules needed:
!   use bmad
!
! Input:
!   field_in -- Field_struct, pointer: Input RF field.
!
! Output:
!   field_out -- Field_struct, pointer: Output RF field.
!-

subroutine transfer_em_field (field_in, field_out)

implicit none

type (em_fields_struct), pointer :: field_in, field_out
type (em_field_mode_struct), pointer :: mode, mode_in, mode_out

integer i

! Rule: If field_in or field_out is associated then %mode must be allocated

if (.not. associated(field_in) .and. .not. associated(field_out)) return

! field_in exists and field_out does not exist: Create field_out.

if (.not. associated(field_out)) then
  call init_em_field (field_out, size(field_in%mode))
  field_out%mode = field_in%mode
  do i = 1, size(field_out%mode)
    mode => field_out%mode(i)
    if (associated(mode%map)) mode%map%n_link = mode%map%n_link + 1
    if (associated(mode%grid)) mode%grid%n_link = mode%grid%n_link + 1
  enddo
  return
endif

! field_in does not exist and field_out exists: Deallocate field_out.

if (.not. associated(field_in)) then
  call init_em_field (field_out, 0)
  return
endif

! Both field_in and field_out exist: If both point to the same memory then need
! to do nothing. Otherwise need to transfer the data.

call init_em_field (field_out, size(field_in%mode))

do i = 1, size(field_out%mode)

  mode_in => field_in%mode(i)
  mode_out => field_out%mode(i)

  if (associated(mode_in%map) .and. associated(mode_out%map)) then
    if (.not. associated(mode_in%map, mode_out%map)) then
      mode_out%map%n_link = mode_out%map%n_link - 1
      if (mode_out%map%n_link == 0) deallocate (mode_out%map)
      mode_out%map => mode_in%map
      mode_out%map%n_link = mode_out%map%n_link + 1
    endif
  elseif (associated(mode_out%map) .and. .not. associated(mode_in%map)) then 
    mode_out%map%n_link = mode_out%map%n_link - 1
    if (mode_out%map%n_link == 0) deallocate (mode_out%map)
  elseif (associated(mode_in%map) .and. .not. associated(mode_out%map)) then 
    mode_out%map => mode_in%map
    mode_out%map%n_link = mode_out%map%n_link + 1
  endif

  if (associated(mode_in%grid) .and. associated(mode_out%grid)) then
    if (.not. associated(mode_in%grid, mode_out%grid)) then
      mode_out%grid%n_link = mode_out%grid%n_link - 1
      if (mode_out%grid%n_link == 0) deallocate (mode_out%grid)
      mode_out%grid => mode_in%grid
      mode_out%grid%n_link = mode_out%grid%n_link + 1
    endif
  elseif (associated(mode_out%grid) .and. .not. associated(mode_in%grid)) then 
    mode_out%grid%n_link = mode_out%grid%n_link - 1
    if (mode_out%grid%n_link == 0) deallocate (mode_out%grid)
  elseif (associated(mode_in%grid) .and. .not. associated(mode_out%grid)) then 
    mode_out%grid => mode_in%grid
    mode_out%grid%n_link = mode_out%grid%n_link + 1
  endif

  mode_out = mode_in

enddo

end subroutine transfer_em_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_wake (wake_in, wake_out)

implicit none

type (wake_struct), pointer :: wake_in, wake_out
integer n_sr_long, n_sr_trans, n_lr

!

if (associated (wake_in)) then
  n_sr_long   = size(wake_in%sr_long%mode)
  n_sr_trans  = size(wake_in%sr_trans%mode)
  n_lr        = size(wake_in%lr)
  call init_wake (wake_out, n_sr_long, n_sr_trans, n_lr)
  wake_out    = wake_in
else
  if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0)
endif

end subroutine transfer_wake

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine deallocate_ele_pointers (ele, nullify_only, nullify_branch, dealloc_poles)
!
! Subroutine to deallocate the pointers in an element.
!
! Note: ele%branch is always nullified. 
!
! Note: For Taylor elements that are slice_slaves: The ele%taylor(:)%term pointers are nullified and 
! not deallocated since these pointers always just point to the lord's corresponding components.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele            -- ele_struct: Element with pointers.
!   nullify_only   -- Logical, optional: If present and True: Nullify & do not deallocate.
!   nullify_branch -- Logical, optional: Nullify ele%branch? Default is True.
!   dealloc_poles  -- Logical, optional: Dealloc ele%a/b_pole, ele%a/b_pole_elec? Default is True.
!
! Output:
!   ele -- Ele_struct: Element with deallocated pointers.
!-

subroutine deallocate_ele_pointers (ele, nullify_only, nullify_branch, dealloc_poles)

implicit none

type (ele_struct), target :: ele
type (em_field_mode_struct), pointer :: mode
logical, optional, intent(in) :: nullify_only, nullify_branch, dealloc_poles
integer i

! %lord and %lat never point to something that has been allocated for the element
! so just nullify these pointers.

if (logic_option(.true., nullify_branch)) nullify (ele%branch)
nullify (ele%lord)

! nullify

if (logic_option (.false., nullify_only)) then
  nullify (ele%wig)
  nullify (ele%rad_int_cache)
  nullify (ele%r)
  nullify (ele%descrip)
  nullify (ele%a_pole, ele%b_pole)
  nullify (ele%a_pole_elec, ele%b_pole_elec)
  nullify (ele%wake)
  forall (i = 1:size(ele%taylor)) ele%taylor(i)%term => null()
  nullify (ele%ptc_genfield%field)
  nullify (ele%ptc_fibre)
  nullify (ele%mode3)
  nullify (ele%wall3d)
  nullify (ele%em_field)
  return
endif

! Normal deallocate.

if (associated (ele%a_pole) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole, ele%b_pole)
endif

if (associated (ele%a_pole_elec) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole_elec, ele%b_pole_elec)
endif

if (associated (ele%rad_int_cache))  deallocate (ele%rad_int_cache)
if (associated (ele%r))              deallocate (ele%r)
if (associated (ele%descrip))        deallocate (ele%descrip)
if (associated (ele%mode3))          deallocate (ele%mode3)
if (associated (ele%wake))        deallocate (ele%wake)

call deallocate_wall3d_pointer (ele%wall3d)

if (associated (ele%em_field)) then
  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    if (associated (mode%map)) then
      mode%map%n_link = mode%map%n_link - 1
      if (mode%map%n_link == 0) deallocate (ele%em_field%mode(i)%map)
    endif
    if (associated (mode%grid)) then
      mode%grid%n_link = mode%grid%n_link - 1
      if (mode%grid%n_link == 0) deallocate (ele%em_field%mode(i)%grid)
    endif
  enddo
  deallocate (ele%em_field)
endif

if (associated(ele%wig)) then
  ele%wig%n_link = ele%wig%n_link - 1
  if (ele%wig%n_link == 0) then
    deallocate (ele%wig)
  else
    nullify (ele%wig)
  endif
endif

if (associated (ele%taylor(1)%term)) then
  if (ele%slave_status == slice_slave$ .and. ele%key == taylor$) then
    forall (i = 1:size(ele%taylor)) ele%taylor(i)%term => null()
  else
    do i = 1, size(ele%taylor)
      deallocate (ele%taylor(i)%term)
    enddo
  endif
endif

call kill_ptc_genfield (ele%ptc_genfield%field)

end subroutine deallocate_ele_pointers

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine kill_ptc_genfield (ptc_genfield)
!
! Subroutine to kill a ptc_genfield.
!
! Modules needed:
!   use bmad
!
! Input:
!   ptc_genfield -- Genfield, pointer: ptc_genfield to kill.
!
! Output:
!   ptc_genfield -- Genfield, pointer: Killed ptc_genfield.
!-

subroutine kill_ptc_genfield (ptc_genfield)

use tpsalie_analysis, only: kill 

implicit none

type (genfield), pointer :: ptc_genfield

!

if (associated(ptc_genfield)) then
  call kill (ptc_genfield)
  deallocate (ptc_genfield)
endif

end subroutine kill_ptc_genfield

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine deallocate_lat_pointers (lat)
!
! Subroutine to deallocate the pointers in a lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lat with pointers.
!
! Output:
!   lat -- lat_struct: Lat with deallocated pointers.
!-

subroutine deallocate_lat_pointers (lat)

implicit none

type (lat_struct) lat
integer i

!

if (associated (lat%ele)) then
  call deallocate_ele_array_pointers (lat%ele)
  call deallocate_ele_pointers (lat%ele_init)
endif

if (allocated(lat%control))  deallocate (lat%control)
if (allocated(lat%ic))       deallocate (lat%ic)

if (allocated(lat%attribute_alias)) deallocate(lat%attribute_alias)

! Do not need to deallocate stuff in lat%branch(0) since
! these pointers have been deallocated above.

if (allocated (lat%branch)) then
  call deallocate_wall3d_pointer (lat%branch(0)%wall3d)

  do i = 1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers (lat%branch(i)%ele)
    deallocate (lat%branch(i)%param, lat%branch(i)%a, lat%branch(i)%b, lat%branch(i)%z)
    call deallocate_wall3d_pointer (lat%branch(i)%wall3d)
  enddo
  deallocate (lat%branch)
endif

!

lat%n_ele_track  = -1
lat%n_ele_max  = -1

end subroutine deallocate_lat_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_ele_array_pointers (eles)
!
! Routine to deallocate the pointers of all the elements in an 
! element array and the array itself.
!
! Modules needed:
!   use bmad
!
! Input:
!   eles(:) -- Ele_struct, pointer: Array of elements.
!
! Output:
!   eles(:) -- Ele_struct, pointer: Deallocated array.
!-

subroutine deallocate_ele_array_pointers (eles)

implicit none

type (ele_struct), pointer :: eles(:)
integer i

!

do i = lbound(eles, 1), ubound(eles, 1)
  call deallocate_ele_pointers (eles(i))
enddo

deallocate (eles)

end subroutine deallocate_ele_array_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_wall3d_pointer (wall3d)
!
! Routine to deallocate a wall3d pointer.
!
! Input:
!   wall3d(:) -- wall3d_struct, pointer: Pointer to wall3d structure.
!
! Output:
!   wall3d(:) -- wall3d_struct, pointer: deallocated
!-

subroutine deallocate_wall3d_pointer (wall3d)

implicit none

type (wall3d_struct), pointer :: wall3d(:)
integer i

!

if (associated (wall3d)) then
  wall3d%n_link = wall3d%n_link - 1
  if (wall3d(1)%n_link == 0) then
    do i = 1, size(wall3d)
      deallocate (wall3d(i)%section)
    enddo
    deallocate (wall3d)
  else
    nullify(wall3d)
  endif
endif

end subroutine deallocate_wall3d_pointer

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)

implicit none

type (coord_struct) orb, orb_temp
type (ele_struct), optional :: ele
real(rp), optional :: vec(:), t_ref_offset, E_photon
integer, optional :: element_end, particle, direction
logical, optional :: shift_vec6

!

orb_temp = coord_struct()
if (present(vec)) orb_temp%vec = vec

call init_coord2 (orb, orb_temp, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)


end subroutine init_coord1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_ref_offset, shift_vec6)

implicit none

type (coord_struct) orb, orb_in, orb_in_save
type (ele_struct), optional, target :: ele

real(rp), optional :: E_photon, t_ref_offset
real(rp) p0c, e_tot, ref_time

integer, optional :: element_end, particle, direction
integer species, dir

logical, optional :: shift_vec6

character(16), parameter :: r_name = 'init_coord1'

! Use temporary orb_in_save so if actual arg for vec, particle, or E_photon
! is part of the orb actual arg things do not get overwriten.

orb_in_save = orb_in  ! Needed if actual args orb and orb_in are the same.
orb = orb_in

species = orb_in_save%species
if (present(particle)) then
  if (particle /= not_set$) species = particle
endif

if (present(particle)) orb%species = particle

if (orb%species == not_set$ .and. present(ele) .and. associated (ele%branch)) then
  orb%species = default_tracking_species(ele%branch%param)
endif

if (orb%species == not_set$) then
  call out_io (s_warn$, r_name, 'NO PARTICLE SPECIES GIVEN. USING POSITRONS AS DEFAULT!')
  orb%species = positron$
endif

orb%state = alive$

dir = integer_option(0, direction)
if (dir == 0) then
  if (orb%species == photon$ .and. orb%vec(6) > 0) orb%direction = 1
  if (orb%species == photon$ .and. orb%vec(6) < 0) orb%direction = -1
else
  orb%direction = dir
  if (orb%species == photon$) orb%vec(6) = orb%direction * abs(orb%vec(6))
endif

! Set location and species

if (present(element_end)) orb%location = element_end

if (orb%location == start_end$) then
  if (orb%direction == 1) then
    orb%location = upstream_end$
  else
    orb%location = downstream_end$
  endif
endif

! Energy values

if (present(ele)) then
  if (.not. present(element_end)) then
    call out_io (s_fatal$, r_name, 'Rule: "element_end" argument must be present if "ele" argument is.')
    call err_exit
  endif
  if (orb%location /= upstream_end$ .or. ele%key == beginning_ele$) then
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    ref_time = ele%ref_time
    if (orb%location == downstream_end$) orb%s = ele%s
  else
    p0c = ele%value(p0c_start$)
    e_tot = ele%value(e_tot_start$)
    ref_time = ele%value(ref_time_start$)
    orb%s = ele%s - ele%value(l$)
  endif
endif

! Photon

if (orb%species == photon$) then

  orb%phase = orb%phase
  orb%field = orb%field
  orb%path_len = 0
  orb%beta = 1

  if (present(ele)) orb%p0c = p0c

  if (present(E_photon)) then
    if (E_photon /= 0) orb%p0c = E_photon
  endif

  ! If original particle is not a photon, photon is launched in same direciton as particle 
  if (orb_in_save%species /= photon$ .and. orb_in_save%species /= not_set$) then
    orb%vec(2:4:2) = orb%vec(2:4:2) / (1 + orb%vec(6))
    if (dir == 0) orb%direction = orb_in_save%direction
  endif
  orb%vec(6) = orb%direction * sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  if (orb%location == downstream_end$) then
    orb%vec(5) = ele%value(l$)
  else
    orb%vec(5) = 0
  endif

else
  orb%spin = orb%spin

endif

! If ele is present...

orb%ix_ele = -1

if (present(ele)) then

  orb%ix_ele = ele%ix_ele
  if (ele%slave_status == slice_slave$ .and. associated(ele%lord)) orb%ix_ele = ele%lord%ix_ele

  if (ele%key == beginning_ele$) orb%location = downstream_end$

  if (orb%species /= photon$) then

    orb%p0c = p0c

    ! E_gun shift. Only time p0c_start /= p0c for an init_ele is when there is an e_gun present in the branch.
    if (logic_option(.true., shift_vec6)) then
      if (ele%key == beginning_ele$) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_start$) - ele%value(p0c$)) / ele%value(p0c$)
      elseif (ele%key == e_gun$ .or. (ele%key == marker$ .and. ele%value(e_tot_ref_init$) /= 0)) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_ref_init$) - ele%value(p0c$)) / ele%value(p0c$)
      endif
    endif

    if (orb%vec(6) == 0) then
      orb%beta = p0c / e_tot
    else
      call convert_pc_to (p0c * (1 + orb%vec(6)), orb%species, beta = orb%beta)
    endif

    ! Do not set %t if %beta = 0 since %t may be a good value.

    if (orb%beta == 0) then
      if (orb%vec(5) /= 0) then
        call out_io (s_error$, r_name, 'Z-POSITION IS NONZERO WITH BETA = 0.', &
                                       'THIS IS NONSENSE SO SETTING Z TO ZERO.')
        orb%vec(5) = 0
      endif

    else
      orb%t = ref_time - orb%vec(5) / (orb%beta * c_light)
      if (present(t_ref_offset)) orb%t = orb%t + t_ref_offset
    endif
  endif

endif

end subroutine init_coord2

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function default_tracking_species (param) result (species)
!
! Routine to return the default species for tracking through a lattice branch.
!
! Input:
!   param  -- lat_param_struct: Parameters for a lattice branch.
!
! Output:
!   species -- Integer: Default species to be used for tracking.
!-

function default_tracking_species (param) result (species)

implicit none

type (lat_param_struct) param
integer species

!

if (param%default_tracking_species == ref_particle$) then
  species = param%particle
elseif (param%default_tracking_species == anti_ref_particle$) then
  species = antiparticle(param%particle)
else
  species = param%default_tracking_species
endif

end function default_tracking_species

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_lat (lat, n)
!
! Subroutine to initialize a Bmad lat.
! 
! Modules needed:
!   use bmad
!
! Input:
!   n    -- Integer, optional: Upper bound lat%ele(0:) array is initialized to.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n)

implicit none

type (lat_struct)  lat

integer, optional :: n

!

call deallocate_lat_pointers (lat)
if (present(n)) call allocate_lat_ele_array(lat, n)
call init_ele (lat%ele_init)
call init_ele (lat%beam_start_ele, def_beam_start$)

call reallocate_control (lat, 100)

lat%title = ' '
lat%use_name = ' '
lat%lattice = ' '
lat%input_file_name = ' '

lat%param = lat_param_struct()
call set_status_flags (lat%param%bookkeeping_state, stale$)

call init_mode_info (lat%a)
call init_mode_info (lat%b)
call init_mode_info (lat%z)

lat%n_ele_track = 0
lat%n_ele_max = 0
lat%n_control_max = 0
lat%n_ic_max = 0
lat%input_taylor_order = 0
lat%version = -1
lat%absolute_time_tracking   = bmad_com%absolute_time_tracking_default

call allocate_branch_array (lat, 0)

!----------------------------------------
contains

subroutine init_mode_info (t)
type (mode_info_struct) t
t%tune = 0
t%emit = 0
t%chrom = 0
end subroutine init_mode_info

end subroutine init_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine init_ele (ele, key, sub_key, ix_ele, branch)
!
! Subroutine to initialize a Bmad element.
!
! Modules needed:
!   use bmad
!
! Input:
!   key       -- Integer, optional: Key to initialize to. EG: quadrupole$, etc.
!   sub_key   -- Integer, optional: Sub-key to initialize to.
!   ix_ele    -- Integer, optional: ix_ele index to initalize to. Default = -1.
!   branch    -- branch_struct: Branch to point ele%branch and ele%ix_branch to.
!                  Otherwise ele%branch is nullified and ele%ix_branch = 0.
!
! Output:
!   ele -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele, key, sub_key, ix_ele, branch)

implicit none

type (ele_struct)  ele
type (branch_struct), optional, target :: branch
integer, optional :: key, sub_key
integer, optional :: ix_ele

!

call deallocate_ele_pointers (ele)

ele = ele_struct()

if (present(branch)) then
  ele%branch => branch
  ele%ix_branch = branch%ix_branch
endif

if (present(ix_ele)) ele%ix_ele = ix_ele

if (present(key)) ele%key = key
if (present(key)) call set_ele_defaults(ele)

if (present(sub_key)) ele%sub_key = sub_key

end subroutine init_ele

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine elec_multipole_init (ele, zero)
!
! Subroutine to allocate memory for the ele%a_pole_elec and ele%b_pole_elec 
! electric multipole vectors.
!
! Modules needed:
!   use bmad
!
! Input:
!   zero -- Logical, optional: If present and True then zero the arrays
!             even if they already exist when this routine is called. 
!             Default is False which means that if the arrays already 
!             exist then this routine will do nothing.
!
! Output:
!   ele -- Ele_struct: Element holding the elec_multipoles.
!     %a_pole_elec(0:n_pole_maxx) -- Elec_multipole An array 
!     %b_pole_elec(0:n_pole_maxx) -- Multipole Bn array
!-

subroutine elec_multipole_init (ele, zero)

implicit none

type (ele_struct) ele
logical, optional :: zero

! If %a_pole_elec and %b_pole_elec already exist then zero them if zero argument present 
! and True.

if (associated (ele%a_pole_elec)) then
  if (logic_option(.false., zero)) then
    ele%a_pole_elec = 0
    ele%b_pole_elec = 0
  endif

! If memory not allocated then allocate and zero.

else
  allocate (ele%a_pole_elec(0:n_pole_maxx), ele%b_pole_elec(0:n_pole_maxx))
  ele%a_pole_elec = 0
  ele%b_pole_elec = 0
endif

end subroutine elec_multipole_init

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_init (ele, zero)
!
! Subroutine to allocate memory for the ele%a_pole and ele%b_pole multipole 
! vectors.
!
! Modules needed:
!   use bmad
!
! Input:
!   zero -- Logical, optional: If present and True then zero the arrays
!             even if they already exist when this routine is called. 
!             Default is False which means that if the arrays already 
!             exist then this routine will do nothing.
!
! Output:
!   ele -- Ele_struct: Element holding the multipoles.
!     %a_pole(0:n_pole_maxx) -- Multipole An array 
!     %b_pole(0:n_pole_maxx) -- Multipole Bn array
!-

subroutine multipole_init (ele, zero)

implicit none

type (ele_struct) ele
logical, optional :: zero

! If %a_pole and %b_pole already exist then zero them if zero argument present 
! and True.

if (associated (ele%a_pole)) then
  if (logic_option(.false., zero)) then
    ele%a_pole = 0
    ele%b_pole = 0
  endif

! If memory not allocated then allocate and zero.

else
  allocate (ele%a_pole(0:n_pole_maxx), ele%b_pole(0:n_pole_maxx))
  ele%a_pole = 0
  ele%b_pole = 0
endif

end subroutine multipole_init

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_em_field (em_field, n_mode)
!
! Subroutine to initialize a em_field_struct pointer.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_mode   -- Integer: Size of %modes(:) to create. If 0, deallocate em_field
!
! Output:
!   em_field -- em_field_struct, pointer: Initialized structure.
!-

subroutine init_em_field (em_field, n_mode)

type (em_fields_struct), pointer :: em_field
type (em_field_mode_struct), pointer :: mode

integer n_mode

integer i

! Cases where nothing is to be done

if (n_mode < 1 .and. .not. associated(em_field)) return

if (n_mode > 0 .and. associated(em_field)) then
  if (size(em_field%mode) == n_mode) return
endif

! Must deallocate existing.

if (associated(em_field)) then
  do i = 1, size(em_field%mode)
    mode => em_field%mode(i)
    if (associated(mode%map)) then
      mode%map%n_link = mode%map%n_link - 1
      if (mode%map%n_link == 0) deallocate (mode%map)
    endif
    if (associated(mode%grid)) then
      mode%grid%n_link = mode%grid%n_link - 1
      if (mode%grid%n_link == 0) deallocate (mode%grid)
    endif
  enddo
  deallocate(em_field)
endif
  
if (n_mode < 1) return

! n_mode > 0 case.

allocate(em_field)
allocate(em_field%mode(n_mode))

end subroutine init_em_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr)
!
! Subroutine to initialize a wake struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_sr_long  -- Integer: Number of terms: wake%nr(n_sr_long).
!   n_sr_trans -- Integer: Number of terms: wake%nr(n_sr_trans).
!   n_lr            -- Integer: Number of terms: wake%nr(n_lr)
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!               If all inputs are 0 then wake is deallocated.
!-

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr)

implicit none

type (wake_struct), pointer :: wake
integer n_sr_long, n_sr_trans, n_lr

! Deallocate wake if all inputs are zero.

if (n_sr_long == 0 .and. n_sr_trans == 0 .and. n_lr == 0) then
  if (associated(wake)) deallocate (wake)
  return
endif

!

if (associated (wake)) then
  if (size(wake%sr_long%mode) /= n_sr_long) then
    deallocate (wake%sr_long%mode)
    allocate (wake%sr_long%mode(n_sr_long))
  endif
  if (size(wake%sr_trans%mode) /= n_sr_trans) then
    deallocate (wake%sr_trans%mode)
    allocate (wake%sr_trans%mode(n_sr_trans))
  endif
  if (size(wake%lr) /= n_lr) then
    deallocate (wake%lr)
    allocate (wake%lr(n_lr))
  endif

else
  allocate (wake)
  allocate (wake%sr_long%mode(n_sr_long))
  allocate (wake%sr_trans%mode(n_sr_trans))
  allocate (wake%lr(n_lr))
endif

end subroutine init_wake

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:)  -- Element array to reallocate.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!   ix_branch   -- Integer, optional: Branch index. Default is 0.
!
! Output:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:) -- Ele_struct, pointer: Resized element array.
!-

subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)

implicit none

type (lat_struct), target :: lat
integer, optional :: upper_bound
integer, optional :: ix_branch
integer ix_br, i

!

ix_br = integer_option (0, ix_branch)

if (ix_br == 0) then
  call allocate_element_array (lat%ele, upper_bound, .true.)
  if (allocated(lat%branch)) then
    do i = 0, ubound(lat%ele, 1)
      lat%ele(i)%branch => lat%branch(0)
    enddo
    lat%branch(0)%ele => lat%ele
  endif

else
  call allocate_element_array (lat%branch(ix_br)%ele, upper_bound, .true.)
  do i = 0, ubound(lat%branch(ix_br)%ele, 1)
    lat%branch(ix_br)%ele(i)%branch => lat%branch(ix_br)
  enddo
  lat%branch(ix_br)%ele%ix_branch = ix_br
endif


end subroutine allocate_lat_ele_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_element_array (ele, upper_bound, init_ele0)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Note: Use allocate_lat_ele_array instead for all ele(:) arrays that
!       are part of a lattice.
!   
!
! Modules needed:
!   use bmad
!
! Input:
!   ele(:)      -- Ele_struct, pointer: Element array.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!   init_ele0   -- Logical, optional: If present and True and ele(:) array has not been allocated then set:
!                     ele(0)%name = 'BEGINNING'
!                     ele(0)%key = beginning_ele$
!                     ele(0)%mat6 = unit matrix
!
! Output:
!   ele(:)      -- Ele_struct, pointer: Allocated element array.
!-

subroutine allocate_element_array (ele, upper_bound, init_ele0)

implicit none

type (ele_struct), pointer :: ele(:)
type (ele_struct), pointer :: temp_ele(:)

integer, optional :: upper_bound
integer curr_ub, ub, i

logical, optional :: init_ele0

! get new size

ub = 10
if (associated (ele)) ub = max (int(1.3*size(ele)), ub)
if (present(upper_bound))  ub = upper_bound

!  save ele if present

if (associated (ele)) then
  if (ub == ubound(ele, 1)) return
  curr_ub = min(ub, ubound(ele, 1))
  do i = curr_ub+1, ubound(ele, 1)
    call deallocate_ele_pointers(ele(i))
  enddo
  temp_ele => ele
  allocate(ele(0:ub))
  call transfer_eles (temp_ele(0:curr_ub), ele(0:curr_ub))
  deallocate (temp_ele)
else
  curr_ub = -1
  allocate(ele(0:ub))
endif

! 

do i = curr_ub+1, ub
  call init_ele (ele(i))
  ele(i)%ix_ele = i
end do

if (logic_option(.false., init_ele0) .and. curr_ub == -1) then
  ele(0)%name = 'BEGINNING'
  ele(0)%key = beginning_ele$
  call mat_make_unit (ele(0)%mat6)
  call set_ele_defaults(ele(0))
endif

end subroutine allocate_element_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_branch_array (lat, upper_bound, lat)
!
! Subroutine to allocate or re-allocate an branch array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- Lat_struct: 
!     %branch(:)  -- Branch array to be allocated.
!   upper_bound -- Integer: Desired upper bound.
! 
! Output:
!   lat         -- Lat_struct: 
!     %branch(:)  -- Allocated branch array.
!-

subroutine allocate_branch_array (lat, upper_bound)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (branch_struct), pointer :: temp_branch(:)

integer :: upper_bound
integer curr_ub, ub, i, j

character(20) :: r_name = 'allocate_branch_array'

!  save branch if present

ub = upper_bound
if (allocated (lat%branch)) then
  if (ub == ubound(lat%branch, 1)) return
  curr_ub = min(ub, ubound(lat%branch, 1))
  allocate (temp_branch(0:curr_ub))
  call transfer_branches (lat%branch(0:curr_ub), temp_branch)
  do i = curr_ub+1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers(lat%branch(i)%ele)
    deallocate(lat%branch(i)%n_ele_track)
    deallocate(lat%branch(i)%n_ele_max)
  enddo
  deallocate (lat%branch)
  allocate(lat%branch(0:ub))
  call transfer_branches (temp_branch(0:curr_ub), lat%branch(0:curr_ub))
  deallocate (temp_branch)
else
  curr_ub = -1
  allocate(lat%branch(0:ub))
  lat%branch(0)%ele            => lat%ele
  lat%branch(0)%param          => lat%param
  lat%branch(0)%a              => lat%a
  lat%branch(0)%b              => lat%b
  lat%branch(0)%z              => lat%z
  lat%branch(0)%n_ele_track    => lat%n_ele_track
  lat%branch(0)%n_ele_max      => lat%n_ele_max
  if (associated(lat%ele)) then
    do i = 0, ubound(lat%ele, 1)
      lat%ele(i)%branch => lat%branch(0)
    enddo
  endif
endif

! 

do i = curr_ub+1, ub
  branch => lat%branch(i)
  branch%lat => lat
  branch%name = ''
  branch%ix_branch = i
  branch%ix_from_branch = -1
  branch%ix_from_ele = -1
  if (i == 0) cycle
  allocate(branch%n_ele_track)
  allocate(branch%n_ele_max)
  allocate(branch%param)
  allocate(branch%a, branch%b, branch%z)
  !!!! branch%param = lat%param
  call set_status_flags (branch%param%bookkeeping_state, stale$)
end do

do i = 0, ub
  branch => lat%branch(i)
  if (.not. associated (branch%ele)) cycle
  do j = 0, ubound(branch%ele, 1)
    branch%ele(j)%branch => branch
  enddo
enddo

end subroutine allocate_branch_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_n (coord, n_coord)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_n (coord, n_coord)

type (coord_struct), allocatable :: coord(:)
type (coord_struct), allocatable :: old(:)

integer, intent(in) :: n_coord
integer i, n_old

!

if (allocated (coord)) then

  n_old = ubound(coord, 1)
  if (n_old >= n_coord) return
  allocate(old(0:n_old))

  do i = 0, n_old
    old(i) = coord(i)
  enddo

  deallocate (coord)
  allocate (coord(0:n_coord))

  do i = 0, n_old
    coord(i) = old(i)
  enddo

  deallocate(old)

else
  allocate (coord(0:n_coord))
endif

end subroutine reallocate_coord_n

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_lat (coord, lat, ix_branch)

type (coord_struct), allocatable :: coord(:)
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer, optional :: ix_branch

!

branch => lat%branch(integer_option(0, ix_branch))

if (allocated(coord)) then
  call reallocate_coord_n (coord, branch%n_ele_max)
else
  allocate (coord(0:branch%n_ele_max))
endif

end subroutine reallocate_coord_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_array (coord_array, lat)
!
! Subroutine to allocate an allocatable coord_array_struct array to
! the proper size for a lattice.
!
! Note: Any old coordinates are not saved except for coord_array(:)%orbit(0).
! If, at input, coord_array is not allocated, coord_array(:)%orbit(0)%vec is set to zero.
! In any case, all other %vec components are set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:) -- Coord_struct, allocatable: Allocatable array.
!   lat      -- lat_struct: 
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord_array (coord_array, lat)

implicit none

type (coord_array_struct), allocatable :: coord_array(:)
type (lat_struct) lat
type (coord_struct), allocatable :: start(:)

integer i, j, nb

!

if (.not. allocated(lat%branch)) return
nb = ubound(lat%branch, 1)

if (allocated (coord_array)) then
  if (size(coord_array) /= nb + 1) then
    call reallocate_coord(start, nb)
    do i = 0, nb
      start(i) = coord_array(i)%orbit(0)
    enddo
    deallocate (coord_array)
    allocate (coord_array(0:nb))
    do i = 0, nb
      call reallocate_coord (coord_array(i)%orbit, lat%branch(i)%n_ele_max)
      coord_array(i)%orbit(0) = start(i)
    enddo
  endif
else
  allocate (coord_array(0:nb))
  do i = 0, nb
    call reallocate_coord (coord_array(i)%orbit, lat%branch(i)%n_ele_max)
  enddo
endif

end subroutine reallocate_coord_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_control(lat, n) 
!
! Function to reallocate the lat%control(:) and lat%ic(:) arrays.
! The old data in the arrays will be saved.
! 
! Modules needed:
!   use bmad
!
! Input:
!   lat  -- Lat_struct: Lattice.
!   n    -- Integer: Array size for lat%control(:) and lat%ic(:).
!
! Output:
!   lat  -- Lat_struct: Lattice.
!     %control(:) -- Control Array with size at least n.
!     %ic(:)      -- Control Array.
!-

subroutine reallocate_control (lat, n)

implicit none

type (lat_struct) lat
type (control_struct), allocatable :: control(:)
integer, intent(in) :: n
integer i, n_old

!

if (.not. allocated(lat%control)) then
  allocate (lat%control(n), lat%ic(n))
  return
endif

n_old = size(lat%control)
if (n_old >= n) return

call move_alloc (lat%control, control)

allocate (lat%control(n))
do i = 1, n_old
  call move_alloc(control(i)%stack, lat%control(i)%stack)
  lat%control(i)%ix_lord   = control(i)%ix_lord
  lat%control(i)%slave     = control(i)%slave
  lat%control(i)%ix_attrib = control(i)%ix_attrib
enddo

call re_allocate(lat%ic, max(n, size(lat%ic) + n - n_old))
lat%ic(n_old+1:) = 0

end subroutine reallocate_control

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_status_flags (bookkeeping_state, stat)
!
! Routine to set the bookkeeping status block.
!
! Input:
!   stat          -- Integer: bookkeeping status. ok$, stale$, etc.
!
! Output:
!   bookkeeping_state -- bookkeeping_state_struct: 
!-

subroutine set_status_flags (bookkeeping_state, stat)

implicit none

type (bookkeeping_state_struct) bookkeeping_state
integer stat

!

bookkeeping_state%control        = stat
bookkeeping_state%s_position     = stat
bookkeeping_state%floor_position = stat
bookkeeping_state%ref_energy     = stat
bookkeeping_state%attributes     = stat
bookkeeping_state%mat6           = stat
bookkeeping_state%rad_int        = stat
bookkeeping_state%ptc            = stat

end subroutine set_status_flags

end module bmad_core_struct_mod
