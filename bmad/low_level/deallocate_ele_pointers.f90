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

use equal_mod, dummy => deallocate_ele_pointers

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
  nullify (ele%ac_kick)
  nullify (ele%control)
  nullify (ele%converter)
  nullify (ele%cartesian_map)
  nullify (ele%cylindrical_map)
  nullify (ele%gen_grad_map)
  nullify (ele%grid_field)
  nullify (ele%ptc_fibre)
  nullify (ele%mode3)
  nullify (ele%photon)
  nullify (ele%rad_map)
  nullify (ele%high_energy_space_charge)
  nullify (ele%wake)
  nullify (ele%wall3d)
  nullify (ele%r)
  nullify (ele%custom)
  nullify (ele%a_pole, ele%b_pole)
  nullify (ele%a_pole_elec, ele%b_pole_elec)
  forall (i = 1:size(ele%taylor)) ele%taylor(i)%term => null()
  forall (i = 0:3) ele%spin_taylor(i)%term => null()
  return
endif

! Normal deallocate.

if (associated (ele%a_pole) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole, ele%b_pole)
endif

if (associated (ele%a_pole_elec) .and. logic_option(.true., dealloc_poles)) then
  deallocate (ele%a_pole_elec, ele%b_pole_elec)
endif

if (associated (ele%descrip))                   deallocate (ele%descrip)
if (associated (ele%control))                   deallocate (ele%control)
if (associated (ele%converter))                 deallocate (ele%converter)
if (associated (ele%rad_map))                   deallocate (ele%rad_map)
if (associated (ele%r))                         deallocate (ele%r)
if (associated (ele%custom))                    deallocate (ele%custom)
if (associated (ele%photon))                    deallocate (ele%photon)
if (associated (ele%mode3))                     deallocate (ele%mode3)
if (associated (ele%wake))                      deallocate (ele%wake)
if (associated (ele%high_energy_space_charge))  deallocate (ele%high_energy_space_charge)
if (associated (ele%gen_grad_map))              deallocate (ele%gen_grad_map)

if (allocated (ele%multipole_cache))            deallocate (ele%multipole_cache)

call unlink_wall3d (ele%wall3d)

if (associated (ele%cartesian_map)) then
  call unlink_fieldmap (cartesian_map = ele%cartesian_map)
endif

if (associated (ele%cylindrical_map)) then
  call unlink_fieldmap (cylindrical_map = ele%cylindrical_map)
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

if (associated (ele%spin_taylor(0)%term)) then
  if (ele%slave_status == slice_slave$ .and. ele%key == taylor$) then
    forall (i = 0:3) ele%spin_taylor(i)%term => null()
  else
    do i = 0, 3
      deallocate (ele%spin_taylor(i)%term)
    enddo
  endif
endif

end subroutine deallocate_ele_pointers

