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

use basic_bmad_interface

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

private reallocate_coord_n, reallocate_coord_lat

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

ele_out%x         = ele_in%x
ele_out%y         = ele_in%y
ele_out%a         = ele_in%a
ele_out%b         = ele_in%b
ele_out%z         = ele_in%z
ele_out%c_mat     = ele_in%c_mat
ele_out%gamma_c   = ele_in%gamma_c
ele_out%mode_flip = ele_in%mode_flip

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
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, ele2's pointers point to the same memory as ele1's.
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
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, branch2's pointers point to the same memory as branch1's.
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
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, branch2's pointers point to the same memory as branch1's.
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
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, lat2's pointers point to the same memory as lat1's.
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

if (associated(wall3d_out)) call unlink_wall3d(wall3d_out)

if (associated(wall3d_in)) then 
  wall3d_out => wall3d_in
  wall3d_out%n_link = wall3d_out%n_link + 1
endif

end subroutine transfer_wall3d

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_exact_bend (exact_bend_in, exact_bend_out)
!
! Subroutine to point exact_bend_out => exact_bend_in
!
! Modules needed:
!   use bmad
!
! Input:
!   exact_bend_in  -- Exact_bend_struct, pointer: Input exact_bendgler field.
!
! Output:
!   exact_bend_out -- Exact_bend_struct, pointer: Output exact_bendgler field.
!-

subroutine transfer_exact_bend (exact_bend_in, exact_bend_out)

implicit none

type (exact_bend_struct), pointer :: exact_bend_in, exact_bend_out

!

if (.not. associated(exact_bend_in) .and. .not. associated(exact_bend_out)) return
if (associated(exact_bend_in, exact_bend_out)) return

! If both associated must be pointing to different memory locations

if (associated(exact_bend_out)) call unlink_exact_bend(exact_bend_out)

if (associated(exact_bend_in)) then 
  exact_bend_out => exact_bend_in
  exact_bend_out%n_link = exact_bend_out%n_link + 1
endif

end subroutine transfer_exact_bend

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_fieldmap (ele_in, ele_out, who)
!
! Subroutine to transfer the field info from one element to another.
! In the end will have:
!     ele_out%cartesian_map    => ele_in%cartesian_map
!     ele_out%cylindrical_map  => ele_in%cylindrical_map
!     ele_out%grid_field       => ele_in%grid_field
!     ele_out%taylor_field     => ele_in%taylor_field
!
! Input:
!   ele_in -- Ele_struct, pointer: Input element.
!   who    -- integer: Possibilities are: all$, cartesian_map$, cylindrical_map$, grid_field$, or taylor_field$
!
! Output:
!   ele_out -- Ele_struct, pointer: Output element.
!-

subroutine transfer_fieldmap (ele_in, ele_out, who)

implicit none

type (ele_struct) :: ele_in, ele_out

integer who
integer i, nm

! Cartesian_map

if (who == all$ .or. who == cartesian_map$) then
  if (associated(ele_in%cartesian_map) .and. associated(ele_out%cartesian_map)) then
    if (size(ele_in%cartesian_map) /= size(ele_out%cartesian_map)) then
      call unlink_fieldmap (cartesian_map = ele_out%cartesian_map)
      nm = size(ele_in%cartesian_map)
      allocate (ele_out%cartesian_map(nm))
      do i = 1, nm
        ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
        ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%cartesian_map)
        if (associated(ele_out%cartesian_map(i)%ptr, ele_in%cartesian_map(i)%ptr)) then
          ele_out%cartesian_map(i) = ele_in%cartesian_map(i) ! Make sure same info
        else
          ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link - 1
          if (ele_out%cartesian_map(i)%ptr%n_link == 0) deallocate (ele_out%cartesian_map(i)%ptr)
          ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
          ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%cartesian_map) .and. .not. associated(ele_out%cartesian_map)) then
    nm = size(ele_in%cartesian_map)
    allocate (ele_out%cartesian_map(nm))
    do i = 1, nm
      ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
      ele_out%cartesian_map(i)%ptr%n_link = ele_in%cartesian_map(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%cartesian_map) .and. associated(ele_out%cartesian_map)) then
    call unlink_fieldmap (cartesian_map = ele_out%cartesian_map)
  endif
endif

! Cylindrical_map

if (who == all$ .or. who == cylindrical_map$) then
  if (associated(ele_in%cylindrical_map) .and. associated(ele_out%cylindrical_map)) then
    if (size(ele_in%cylindrical_map) /= size(ele_out%cylindrical_map)) then
      call unlink_fieldmap (cylindrical_map = ele_out%cylindrical_map)
      nm = size(ele_in%cylindrical_map)
      allocate (ele_out%cylindrical_map(nm))
      do i = 1, nm
        ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
        ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%cylindrical_map)
        if (associated(ele_out%cylindrical_map(i)%ptr, ele_in%cylindrical_map(i)%ptr)) then
          ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
        else
          ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link - 1
          if (ele_out%cylindrical_map(i)%ptr%n_link == 0) deallocate (ele_out%cylindrical_map(i)%ptr)
          ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
          ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%cylindrical_map) .and. .not. associated(ele_out%cylindrical_map)) then
    nm = size(ele_in%cylindrical_map)
    allocate (ele_out%cylindrical_map(nm))
    do i = 1, nm
      ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
      ele_out%cylindrical_map(i)%ptr%n_link = ele_in%cylindrical_map(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%cylindrical_map) .and. associated(ele_out%cylindrical_map)) then
    call unlink_fieldmap (cylindrical_map = ele_out%cylindrical_map)
  endif
endif

! Grid_field

if (who == all$ .or. who == grid_field$) then
  if (associated(ele_in%grid_field) .and. associated(ele_out%grid_field)) then
    if (size(ele_in%grid_field) /= size(ele_out%grid_field)) then
      call unlink_fieldmap (grid_field = ele_out%grid_field)
      nm = size(ele_in%grid_field)
      allocate (ele_out%grid_field(nm))
      do i = 1, nm
        ele_out%grid_field(i) = ele_in%grid_field(i)
        ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%grid_field)
        if (associated(ele_out%grid_field(i)%ptr, ele_in%grid_field(i)%ptr)) then
          ele_out%grid_field(i) = ele_in%grid_field(i)
        else
          ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link - 1
          if (ele_out%grid_field(i)%ptr%n_link == 0) deallocate (ele_out%grid_field(i)%ptr)
          ele_out%grid_field(i) = ele_in%grid_field(i)
          ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%grid_field) .and. .not. associated(ele_out%grid_field)) then
    nm = size(ele_in%grid_field)
    allocate (ele_out%grid_field(nm))
    do i = 1, nm
      ele_out%grid_field(i) = ele_in%grid_field(i)
      ele_out%grid_field(i)%ptr%n_link = ele_in%grid_field(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%grid_field) .and. associated(ele_out%grid_field)) then
    call unlink_fieldmap (grid_field = ele_out%grid_field)
  endif
endif

! Taylor_field

if (who == all$ .or. who == taylor_field$) then
  if (associated(ele_in%taylor_field) .and. associated(ele_out%taylor_field)) then
    if (size(ele_in%taylor_field) /= size(ele_out%taylor_field)) then
      call unlink_fieldmap (taylor_field = ele_out%taylor_field)
      nm = size(ele_in%taylor_field)
      allocate (ele_out%taylor_field(nm))
      do i = 1, nm
        ele_out%taylor_field(i) = ele_in%taylor_field(i)
        ele_out%taylor_field(i)%ptr%n_link = ele_out%taylor_field(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%taylor_field)
        if (associated(ele_out%taylor_field(i)%ptr, ele_in%taylor_field(i)%ptr)) then
          ele_out%taylor_field(i) = ele_in%taylor_field(i)
        else
          ele_out%taylor_field(i)%ptr%n_link = ele_out%taylor_field(i)%ptr%n_link - 1
          if (ele_out%taylor_field(i)%ptr%n_link == 0) deallocate (ele_out%taylor_field(i)%ptr)
          ele_out%taylor_field(i) = ele_in%taylor_field(i)
          ele_out%taylor_field(i)%ptr%n_link = ele_out%taylor_field(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%taylor_field) .and. .not. associated(ele_out%taylor_field)) then
    nm = size(ele_in%taylor_field)
    allocate (ele_out%taylor_field(nm))
    do i = 1, nm
      ele_out%taylor_field(i) = ele_in%taylor_field(i)
      ele_out%taylor_field(i)%ptr%n_link = ele_in%taylor_field(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%taylor_field) .and. associated(ele_out%taylor_field)) then
    call unlink_fieldmap (taylor_field = ele_out%taylor_field)
  endif
endif

end subroutine transfer_fieldmap

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
logical, optional, intent(in) :: nullify_only, nullify_branch, dealloc_poles
integer i

! %lord and %lat never point to something that has been allocated for the element
! so just nullify these pointers.

if (logic_option(.true., nullify_branch)) nullify (ele%branch)
nullify (ele%lord)

! nullify

if (logic_option (.false., nullify_only)) then
  nullify (ele%descrip)
  nullify (ele%control_var)
  nullify (ele%cartesian_map)
  nullify (ele%cylindrical_map)
  nullify (ele%exact_bend)
  nullify (ele%taylor_field)
  nullify (ele%grid_field)
  nullify (ele%ptc_fibre)
  nullify (ele%mode3)
  nullify (ele%photon)
  nullify (ele%rad_int_cache)
  nullify (ele%space_charge)
  nullify (ele%wake)
  nullify (ele%wall3d)
  nullify (ele%r)
  nullify (ele%a_pole, ele%b_pole)
  nullify (ele%a_pole_elec, ele%b_pole_elec)
  forall (i = 1:size(ele%taylor)) ele%taylor(i)%term => null()
  nullify (ele%ptc_genfield%field)
  return
endif

! Normal deallocate.

if (associated (ele%a_pole) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole, ele%b_pole)
endif

if (associated (ele%a_pole_elec) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole_elec, ele%b_pole_elec)
endif

if (associated (ele%descrip))        deallocate (ele%descrip)
if (associated (ele%control_var))    deallocate (ele%control_var)
if (associated (ele%rad_int_cache))  deallocate (ele%rad_int_cache)
if (associated (ele%r))              deallocate (ele%r)
if (associated (ele%photon))         deallocate (ele%photon)
if (associated (ele%mode3))          deallocate (ele%mode3)
if (associated (ele%wake))           deallocate (ele%wake)
if (associated (ele%space_charge))   deallocate (ele%space_charge)

call unlink_wall3d (ele%wall3d)

if (associated (ele%cartesian_map)) then
  call unlink_fieldmap (cartesian_map = ele%cartesian_map)
endif

if (associated (ele%cylindrical_map)) then
  call unlink_fieldmap (cylindrical_map = ele%cylindrical_map)
endif

call unlink_exact_bend (ele%exact_bend)

if (associated (ele%taylor_field)) then
  call unlink_fieldmap (taylor_field = ele%taylor_field)
endif

if (associated (ele%grid_field)) then
  call unlink_fieldmap (grid_field = ele%grid_field)
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
  call unlink_wall3d (lat%branch(0)%wall3d)

  do i = 1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers (lat%branch(i)%ele)
    deallocate (lat%branch(i)%param, lat%branch(i)%a, lat%branch(i)%b, lat%branch(i)%z)
    call unlink_wall3d (lat%branch(i)%wall3d)
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

if (.not. associated(eles)) return

do i = lbound(eles, 1), ubound(eles, 1)
  call deallocate_ele_pointers (eles(i))
enddo

deallocate (eles)

end subroutine deallocate_ele_array_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine unlink_wall3d (wall3d)
!
! Routine to deallocate a wall3d pointer.
!
! Input:
!   wall3d(:) -- wall3d_struct, pointer: Pointer to wall3d structure.
!
! Output:
!   wall3d(:) -- wall3d_struct, pointer: deallocated
!-

subroutine unlink_wall3d (wall3d)

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

end subroutine unlink_wall3d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine unlink_exact_bend (exact_bend)
!
! Routine to deallocate a exact_bend pointer.
!
! Input:
!   exact_bend -- exact_bend_struct, pointer: Pointer to exact_bend structure.
!
! Output:
!   exact_bend -- exact_bend_struct, pointer: deallocated
!-

subroutine unlink_exact_bend (exact_bend)

implicit none

type (exact_bend_struct), pointer :: exact_bend
integer i

!

if (associated (exact_bend)) then
  exact_bend%n_link = exact_bend%n_link - 1
  if (exact_bend%n_link == 0) then
    deallocate (exact_bend)
  else
    nullify(exact_bend)
  endif
endif

end subroutine unlink_exact_bend

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
! Subroutine unlink_fieldmap (cartesian_map, cylindrical_map, taylor_field, grid_field)
!
! Subroutine to unlink the field components of an element.
!
! Input:
!   cartesian_map(:)   -- cartesian_map_struct, pointer, optional: cartesian_map component.
!   cylindrical_map(:) -- cylindrical_map_struct, pointer, optional: cylindrical_map component.
!   taylor_field(:)    -- taylor_field_struct, pointer, optional: taylor_field component.
!   grid_field(:)      -- grid_field_struct, pointer, optional: grid_field component.
!
! Output:
!   cartesian_map(:)   -- cartesian_map_struct, pointer, optional: cartesian_map component.
!   cylindrical_map(:) -- cylindrical_map_struct, pointer, optional: cylindrical_map component.
!   taylor_field(:)    -- taylor_field_struct, pointer, optional: taylor_field component.
!   grid_field(:)      -- grid_field_struct, pointer, optional: grid_field component.
!-

subroutine unlink_fieldmap (cartesian_map, cylindrical_map, taylor_field, grid_field)

type (cartesian_map_struct), pointer, optional :: cartesian_map(:)
type (cylindrical_map_struct), pointer, optional :: cylindrical_map(:)
type (taylor_field_struct), pointer, optional :: taylor_field(:)
type (grid_field_struct), pointer, optional :: grid_field(:)

integer i

! In theory, %prt should always be associated but a bad digested file can leave %ptr unassociated.

if (present(cartesian_map)) then
  do i = 1, size(cartesian_map)
    if (.not. associated(cartesian_map(i)%ptr)) cycle
    cartesian_map(i)%ptr%n_link = cartesian_map(i)%ptr%n_link - 1
    if (cartesian_map(i)%ptr%n_link == 0) deallocate (cartesian_map(i)%ptr)
  enddo
  deallocate (cartesian_map)
endif

if (present(cylindrical_map)) then
  do i = 1, size(cylindrical_map)
    if (.not. associated(cylindrical_map(i)%ptr)) cycle
    cylindrical_map(i)%ptr%n_link = cylindrical_map(i)%ptr%n_link - 1
    if (cylindrical_map(i)%ptr%n_link == 0) deallocate (cylindrical_map(i)%ptr)
  enddo
  deallocate (cylindrical_map)
endif

if (present(grid_field)) then
  do i = 1, size(grid_field)
    if (.not. associated(grid_field(i)%ptr)) cycle
    grid_field(i)%ptr%n_link = grid_field(i)%ptr%n_link - 1
    if (grid_field(i)%ptr%n_link == 0) deallocate (grid_field(i)%ptr)
  enddo
  deallocate (grid_field)
endif

if (present(taylor_field)) then
  do i = 1, size(taylor_field)
    if (.not. associated(taylor_field(i)%ptr)) cycle
    taylor_field(i)%ptr%n_link = taylor_field(i)%ptr%n_link - 1
    if (taylor_field(i)%ptr%n_link == 0) deallocate (taylor_field(i)%ptr)
  enddo
  deallocate (taylor_field)
endif

end subroutine unlink_fieldmap

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

character(*), parameter :: r_name = 'reallocate_coord_n'

!

if (allocated (coord)) then

  if (lbound(coord, 1) /= 0) then
    call out_io (s_fatal$, r_name, 'ORBIT ARRAY LOWER BOUND NOT EQUAL TO ZERO!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

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
! Routine to reallocate the lat%control(:) and lat%ic(:) arrays.
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
  lat%control(i)%lord      = control(i)%lord
  lat%control(i)%slave     = control(i)%slave
  lat%control(i)%ix_attrib = control(i)%ix_attrib
enddo

call re_allocate(lat%ic, max(n, size(lat%ic) + n - n_old))
lat%ic(n_old+1:) = 0

end subroutine reallocate_control

end module bmad_core_struct_mod
