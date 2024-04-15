module equal_mod

use bmad_routine_interface

interface assignment (=)
  module procedure ele_equal_ele
  module procedure ele_vec_equal_ele_vec
  module procedure lat_equal_lat 
  module procedure lat_vec_equal_lat_vec 
  module procedure branch_equal_branch
  module procedure taylor_equal_taylor
  module procedure taylors_equal_taylors
  module procedure em_taylor_equal_em_taylor
  module procedure em_taylors_equal_em_taylors
  module procedure complex_taylor_equal_complex_taylor
  module procedure complex_taylors_equal_complex_taylors
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

interface operator (+)
  module procedure em_field_plus_em_field
end interface

interface operator (*)
  module procedure map1_times_map1
end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function map1_times_map1 (map2, map1) result (map_out)
!
! Routine to concatenate two spin orbital linear maps.
!     map_out = map2(map1)
! Order is like applying matrices. map1 is before map2.
!
! Input:
!   map2      -- spin_orbit_map1_struct: Second map.
!   map1      -- spin_orbit_map1_struct: First map.
!
! Output:
!   map_out   -- spin_orbit_map1_struct: Concatenated map.
!-

function map1_times_map1 (map2, map1) result (map_out)

type (spin_orbit_map1_struct), intent(in) :: map1, map2
type (spin_orbit_map1_struct) map_out
real(rp) m2q1(0:3,6)
integer i

!

map_out%orb_mat = matmul(map2%orb_mat, map1%orb_mat)
map_out%vec0 = matmul(map2%orb_mat, map1%vec0) + map2%vec0

map_out%spin_q(:,0) = quat_mul(map2%spin_q(:,0), map1%spin_q(:,0))

m2q1 = matmul(map2%spin_q(:,1:6), map1%orb_mat)
do i = 1, 6
  map_out%spin_q(:,i) = quat_mul(map2%spin_q(:,0), map1%spin_q(:,i)) + quat_mul(m2q1(:,i), map1%spin_q(:,0))
enddo

end function map1_times_map1

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function em_field_plus_em_field (field1, field2) result (field_tot)
!
! Routine to add fields.
!
! Note: This subroutine is called by the overloaded plus sign:
!		field_tot = field1 + field2 
!
! Input:
!   field1 -- em_field_struct: Input field
!   field2 -- em_field_struct: Input field
!
! Output:
!   field_tot -- em_field_struct: Combined field.
!-

function em_field_plus_em_field (field1, field2) result (field_tot)

type (em_field_struct), intent(in) :: field1, field2
type (em_field_struct) field_tot

!

field_tot%e = field1%e + field2%e
field_tot%b = field1%b + field2%b

field_tot%de = field1%de + field2%de
field_tot%db = field1%db + field2%db

end function em_field_plus_em_field 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equal_ele (ele_out, ele_in)
!
! Subroutine that is used to set one element equal to another. 
! This routine takes care of the pointers in ele_out. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele_out = ele_in
!
! Input:
!   ele_in  -- Ele_struct: Input element.
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine ele_equal_ele (ele_out, ele_in)

implicit none
	
type (ele_struct), intent(inout), target :: ele_out
type (ele_struct), intent(in), target :: ele_in

call ele_equals_ele(ele_out, ele_in, .true.)

end subroutine ele_equal_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equals_ele (ele_out, ele_in, update_nametable)
!
! Subroutine that is used to set an element equal to another.
! Note: Use ele_equal_ele instead unless you know what you are doing.
! 
!
! Input:
!   ele_in            -- Ele_struct: Input element.
!   update_nametable  -- logical: If true, update the nametable. If false, do not.
!                         Note: nametable updates can take time if this routine
!                         is called a many times. See remove_eles_from_lat as an example.
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine ele_equals_ele (ele_out, ele_in, update_nametable)

implicit none
	
type (ele_struct), intent(inout), target :: ele_out
type (ele_struct), intent(in), target :: ele_in
type (ele_struct) ele_save
type (nametable_struct), pointer :: nt
type (converter_sub_distribution_struct), pointer :: sd_in, sd_out

integer i, j, n, n1, n2, ub(2), ub1
logical update_nametable, comensurate

! 1) Save ele_out pointers in ele_save
! 2) Set ele_out = ele_in.

call transfer_ele (ele_out, ele_save)
call transfer_ele (ele_in, ele_out)

! ele_out%ix_ele and ele_out%ix_branch should not change.
! ele_out%branch should not change if ele_out is a component of a lat_struct.
!   Otherwise ele_out%lat should point to ele_in%lat (For cases where ele_out 
!   represents a sliced piece of ele_in)

ele_out%ix_ele    = ele_save%ix_ele    ! This should not change.
ele_out%ix_branch = ele_save%ix_branch ! This should not change.

if (update_nametable .and. ele_out%ix_ele > -1) then          ! If part of a lattice...
  ele_out%branch => ele_save%branch    !   then ele_out%branch should not change.
  if (associated(ele_out%branch)) then
    n = ele_nametable_index(ele_out)
    nt => ele_out%branch%lat%nametable
    ! During parsing the nametable may not have yet been updated so do not modify the
    ! nametable if this is the case.
    if (n <= nt%n_max .and. allocated(nt%name)) then
      if (nt%name(n) /= ele_out%name) call nametable_change1(ele_out%branch%lat%nametable, ele_out%name, n)
    endif
  endif
endif

! Transfer pointer info.
! When finished ele_out's pointers will be pointing to a different memory
! location from ele_in's so that the elements are separate.
! Exceptions: %em_field%mode%cylindrical_map, %em_field%mode%grid.

! %cartesian_map

ele_out%cartesian_map => ele_save%cartesian_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, cartesian_map$)

! %ac_kick

ele_out%ac_kick => ele_save%ac_kick  ! reinstate
call transfer_ac_kick (ele_in%ac_kick, ele_out%ac_kick)

! %cylindrical_map, etc.

ele_out%cylindrical_map => ele_save%cylindrical_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, cylindrical_map$)

ele_out%gen_grad_map => ele_save%gen_grad_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, gen_grad_map$)

ele_out%grid_field => ele_save%grid_field ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, grid_field$)

! %rad_map

if (associated(ele_in%rad_map)) then
  if (associated (ele_save%rad_map)) then
      ele_out%rad_map => ele_save%rad_map
  else
    allocate (ele_out%rad_map)
  endif
  ele_out%rad_map = ele_in%rad_map
else
  if (associated (ele_save%rad_map)) deallocate (ele_save%rad_map)
endif

! %r

if (associated(ele_in%r)) then
  if (associated (ele_save%r)) then
    if (all(lbound(ele_save%r) == lbound(ele_in%r)) .and. &
        all(ubound(ele_save%r) == ubound(ele_in%r)) ) then
      ele_out%r => ele_save%r
    else
      deallocate (ele_save%r)
      allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                       lbound(ele_in%r,2):ubound(ele_in%r,2), &
                       lbound(ele_in%r,3):ubound(ele_in%r,3)))
    endif
  else
    allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                     lbound(ele_in%r,2):ubound(ele_in%r,2), &
                     lbound(ele_in%r,3):ubound(ele_in%r,3)))
  endif
  ele_out%r = ele_in%r
else
  if (associated (ele_save%r)) deallocate (ele_save%r)
endif

! %photon

if (associated(ele_in%photon)) then
  ele_out%photon => ele_save%photon  ! reinstate
  if (.not. associated(ele_out%photon)) allocate(ele_out%photon)

  if (allocated (ele_in%photon%grid%pt)) then
    ub = ubound(ele_in%photon%grid%pt)
    if (allocated (ele_out%photon%grid%pt)) then
      if (any(ub /= ubound(ele_out%photon%grid%pt))) deallocate (ele_out%photon%grid%pt)
    endif
    if (.not. allocated (ele_out%photon%grid%pt)) allocate (ele_out%photon%grid%pt(0:ub(1), 0:ub(2)))
  else
    if (allocated(ele_out%photon%grid%pt)) deallocate (ele_out%photon%grid%pt)
  endif

  if (allocated (ele_in%photon%pixel%pt)) then
    ub = ubound(ele_in%photon%pixel%pt)
    if (allocated (ele_out%photon%pixel%pt)) then
      if (any(ub /= ubound(ele_out%photon%pixel%pt))) deallocate (ele_out%photon%pixel%pt)
    endif
    if (.not. allocated (ele_out%photon%pixel%pt)) allocate (ele_out%photon%pixel%pt(0:ub(1), 0:ub(2)))
  else
    if (allocated(ele_out%photon%pixel%pt)) deallocate (ele_out%photon%pixel%pt)
  endif

  ele_out%photon = ele_in%photon
else
  if (associated (ele_save%photon)) deallocate (ele_save%photon)
endif

! %control

if (associated(ele_in%control)) then
  n1 = size(ele_in%control%var)
  n2 = -1
  if (allocated(ele_in%control%ramp)) n2 = size(ele_in%control%ramp)

  ele_out%control => ele_save%control   ! reinstate

  if (associated(ele_out%control)) then
    if (size(ele_out%control%var) /= n1)  deallocate(ele_out%control%var)
    if (allocated(ele_out%control%ramp)) then
      if (size(ele_out%control%ramp) /= n2) deallocate(ele_out%control%ramp)
    endif
  endif

  if (.not. associated(ele_out%control)) allocate(ele_out%control)
  if (.not. allocated(ele_out%control%var)) allocate(ele_out%control%var(n1))
  if (.not. allocated(ele_out%control%ramp) .and. n2 > -1) allocate(ele_out%control%ramp(n2))

  ele_out%control = ele_in%control

else
  if (associated (ele_save%control)) deallocate(ele_save%control)
endif

! %converter

if (associated(ele_in%converter)) then
  n = size(ele_in%converter%dist)
  ele_out%converter => ele_save%converter   ! reinstate
  if (associated(ele_out%converter)) then
    comensurate = .false.
    if (size(ele_out%converter%dist) /= n) goto 100
    do i = 1, n
      if (size(ele_in%converter%dist(i)%sub_dist) /= size(ele_out%converter%dist(i)%sub_dist)) goto 100
      do j = 1, size(ele_in%converter%dist(i)%sub_dist)
        sd_in => ele_in%converter%dist(i)%sub_dist(j)
        sd_out => ele_out%converter%dist(i)%sub_dist(j)
        if (.not. all(sd_in%prob_pc_r%prob == sd_out%prob_pc_r%prob)) goto 100
        if (size(sd_in%dir_out%beta%fit_1d_r) /= size(sd_out%dir_out%beta%fit_1d_r)) goto 100
        if (size(sd_in%dir_out%alpha_x%fit_1d_r) /= size(sd_out%dir_out%alpha_x%fit_1d_r)) goto 100
        if (size(sd_in%dir_out%alpha_y%fit_1d_r) /= size(sd_out%dir_out%alpha_y%fit_1d_r)) goto 100
      enddo
    enddo
    comensurate = .true.
    100 continue
    if (.not. comensurate) deallocate(ele_out%converter)
  endif

  if (.not. associated(ele_out%converter)) then
    allocate(ele_out%converter)
    allocate(ele_out%converter%dist(n))
  endif
  ele_out%converter = ele_in%converter

else
  if (associated (ele_save%converter)) deallocate(ele_save%converter)
endif

! %foil

if (associated(ele_in%foil)) then
  n = size(ele_in%foil%material)
  ele_out%foil => ele_save%foil   ! reinstate
  if (associated(ele_out%foil)) then
    comensurate = .false.
    if (size(ele_out%foil%material) /= n) deallocate(ele_out%foil)
  endif

  if (.not. associated(ele_out%foil)) then
    allocate(ele_out%foil)
    allocate(ele_out%foil%material(n))
  endif
  ele_out%foil = ele_in%foil

else
  if (associated (ele_save%foil)) deallocate(ele_save%foil)
endif

! %taylor

do i = 1, 6
  ele_out%taylor(i)%term => ele_save%taylor(i)%term ! reinstate
  ele_out%taylor(i) = ele_in%taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %spin_taylor

do i = 0, 3
  ele_out%spin_taylor(i)%term => ele_save%spin_taylor(i)%term ! reinstate
  ele_out%spin_taylor(i) = ele_in%spin_taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %wall3d

ele_out%wall3d => ele_save%wall3d        ! reinstate
call transfer_wall3d (ele_in%wall3d, ele_out%wall3d)

! %a_pole, and %b_pole

if (associated(ele_in%a_pole)) then
  ele_out%a_pole => ele_save%a_pole   ! reinstate
  ele_out%b_pole => ele_save%b_pole   ! reinstate
  call multipole_init (ele_out, magnetic$)
  ele_out%a_pole = ele_in%a_pole
  ele_out%b_pole = ele_in%b_pole
else
  if (associated (ele_save%a_pole)) deallocate (ele_save%a_pole, ele_save%b_pole)
endif

! %a_pole_elec, and %b_pole_elec

if (associated(ele_in%a_pole_elec)) then
  ele_out%a_pole_elec => ele_save%a_pole_elec   ! reinstate
  ele_out%b_pole_elec => ele_save%b_pole_elec   ! reinstate
  call multipole_init (ele_out, electric$)
  ele_out%a_pole_elec = ele_in%a_pole_elec
  ele_out%b_pole_elec = ele_in%b_pole_elec
else
  if (associated (ele_save%a_pole_elec)) deallocate (ele_save%a_pole_elec, ele_save%b_pole_elec)
endif

! %custom

if (associated(ele_in%custom)) then
  if (associated (ele_save%custom)) then
    if (size(ele_save%custom) == size(ele_in%custom)) then
      ele_out%custom => ele_save%custom
    else
      deallocate (ele_save%custom)
      allocate (ele_out%custom(size(ele_in%custom)))
    endif
  else
    allocate (ele_out%custom(size(ele_in%custom)))
  endif
  ele_out%custom = ele_in%custom
else
  if (associated (ele_save%custom)) deallocate (ele_save%custom)
endif

! %descrip

if (associated(ele_in%descrip)) then
  if (associated (ele_save%descrip)) then
    ele_out%descrip => ele_save%descrip
  else
    allocate (ele_out%descrip)
  endif
  ele_out%descrip = ele_in%descrip
else
  if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
endif

! %mode3

if (associated(ele_in%mode3)) then
  if (associated (ele_save%mode3)) then
    ele_out%mode3 => ele_save%mode3
  else
    allocate (ele_out%mode3)
  endif
  ele_out%mode3 = ele_in%mode3
else
  if (associated (ele_save%mode3)) deallocate (ele_save%mode3)
endif

! %high_energy_space_charge

if (associated(ele_in%high_energy_space_charge)) then
  if (associated (ele_save%high_energy_space_charge)) then
    ele_out%high_energy_space_charge => ele_save%high_energy_space_charge
  else
    allocate (ele_out%high_energy_space_charge)
  endif
  ele_out%high_energy_space_charge = ele_in%high_energy_space_charge
else
  if (associated (ele_save%high_energy_space_charge)) deallocate (ele_save%high_energy_space_charge)
endif

! %wake

ele_out%wake => ele_save%wake  ! reinstate
call transfer_wake (ele_in%wake, ele_out%wake)

end subroutine ele_equals_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_vec_equal_ele_vec (ele1, ele2)
!
! Subroutine that is used to set one ele vector equal to another. 
! This routine takes care of the pointers in ele1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele1(:) = ele2(:)
!
! Input:
!   ele2(:) -- ele_struct: Input ele vector.
!
! Output:
!   ele1(:) -- ele_struct: Output ele vector.
!-

subroutine ele_vec_equal_ele_vec (ele1, ele2)

implicit none
	
type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(in) :: ele2(:)

integer i

! error check

if (size(ele1) /= size(ele2)) then
  print *, 'ERROR IN ele_vec_equal_ele_vec: ARRAY SIZES ARE NOT THE SAME!'
  if (global_com%exit_on_error) call err_exit
endif

! transfer

do i = 1, size(ele1)
  call ele_equal_ele (ele1(i), ele2(i))
enddo

end subroutine ele_vec_equal_ele_vec 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_equal_lat (lat_out, lat_in)
!
! Subroutine that is used to set one lat equal to another. 
! This routine takes care of the pointers in lat_in. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat_out = lat_in
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat.
!-

subroutine lat_equal_lat (lat_out, lat_in)

implicit none

type (lat_struct), intent(inout), target :: lat_out
type (lat_struct), intent(in), target :: lat_in
type (branch_struct), pointer :: branch_out
type (control_struct), pointer :: c_in, c_out
integer i, n, nb, ne, ie, n_out, n_in
logical do_alloc

! Kill allociated PTC layouts if they exist

call kill_ptc_layouts(lat_out)

! If the element arrays have not been initialized in lat_in then deallocate lat_out.

if (.not. allocated (lat_in%branch) .or. .not. associated(lat_in%ele)) then
  call deallocate_lat_pointers (lat_out)
  return
endif

! Care must be taken here since lat%ele points to the same memory as lat%branch(0).
! First take care of the branch lines.

nb = ubound(lat_in%branch, 1)
call allocate_branch_array (lat_out, nb)

do i = 0, nb
  do_alloc = .true.
  if (associated(lat_out%branch(i)%ele)) then
    if (ubound(lat_out%branch(i)%ele, 1) >= lat_in%branch(i)%n_ele_max) do_alloc = .false.
  endif
  if (do_alloc) then
    ne = min(lat_in%branch(i)%n_ele_max+10, ubound(lat_in%branch(i)%ele, 1))
    call allocate_lat_ele_array (lat_out, ne, i)
  endif
  branch_out => lat_out%branch(i)
  branch_out = lat_in%branch(i)
  branch_out%lat => lat_out
  do ie = 0, ubound(branch_out%ele, 1)
    branch_out%ele(ie)%ix_ele = ie
    branch_out%ele(ie)%ix_branch = i
    branch_out%ele(ie)%branch => branch_out
  enddo
enddo

lat_out%ele_init = lat_in%ele_init
nullify(lat_out%ele_init%branch)

! non-pointer transfer

call transfer_lat_parameters (lat_in, lat_out)
call ramper_slave_setup(lat_out, .true.)


end subroutine lat_equal_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_vec_equal_lat_vec (lat1, lat2)
!
! Subroutine that is used to set one lat vector equal to another. 
! This routine takes care of the pointers in lat1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat1(:) = lat2(:)
!
! Input:
!   lat2(:) -- lat_struct: Input lat vector.
!
! Output:
!   lat1(:) -- lat_struct: Output lat vector.
!-

subroutine lat_vec_equal_lat_vec (lat1, lat2)

implicit none
	
type (lat_struct), intent(inout) :: lat1(:)
type (lat_struct), intent(in) :: lat2(:)

integer i

! error check

if (size(lat1) /= size(lat2)) then
  print *, 'ERROR IN lat_vec_equal_lat_vec: ARRAY SIZES ARE NOT THE SAME!'
  if (global_com%exit_on_error) call err_exit
endif

! transfer

do i = 1, size(lat1)
  call lat_equal_lat (lat1(i), lat2(i))
enddo

end subroutine lat_vec_equal_lat_vec 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine branch_equal_branch (branch1, branch2)
!
! Subroutine that is used to set one branch equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		branch1 = branch2
!
! Input:
!   branch2 -- branch_struct: Input branch.
!
! Output:
!   branch1 -- branch_struct: Output branch.
!-

subroutine branch_equal_branch (branch1, branch2)

implicit none
	
type (branch_struct), intent(inout) :: branch1
type (branch_struct), intent(in) :: branch2
integer i, n
logical do_alloc

!

do_alloc = .true.
if (associated(branch1%ele)) then
  if (ubound(branch1%ele, 1) >= branch2%n_ele_max) do_alloc = .false.
endif

if (do_alloc) then
  n = min(branch2%n_ele_max+10, ubound(branch2%ele, 1))
  call allocate_element_array (branch1%ele, n)
endif

do i = 0, branch2%n_ele_max
  branch1%ele(i)  = branch2%ele(i)
enddo

call transfer_wall3d (branch2%wall3d, branch1%wall3d)

call transfer_branch_parameters(branch2, branch1)

end subroutine branch_equal_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine coord_equal_coord (coord1, coord2)
!
! Subroutine that is used to set one coord equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		coord1 = coord2
!
! Input:
!   coord2 -- coord_struct: Input coord.
!
! Output:
!   coord1 -- coord_struct: Output coord.
!-

elemental subroutine coord_equal_coord (coord1, coord2)

implicit none
	
type (coord_struct), intent(inout) :: coord1
type (coord_struct), intent(in) :: coord2

!

coord1%vec = coord2%vec
coord1%spin = coord2%spin
 
end subroutine coord_equal_coord

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine taylor_equal_taylor (taylor1, taylor2)
!
! Subroutine that is used to set one taylor equal to another. 
! This routine takes care of the pointers in taylor1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		taylor1 = taylor2
!
! Input:
!   taylor2 -- Taylor_struct: Input taylor.
!
! Output:
!   taylor1 -- Taylor_struct: Output taylor.
!-

subroutine taylor_equal_taylor (taylor1, taylor2)

implicit none
	
type (taylor_struct), intent(inout) :: taylor1
type (taylor_struct), intent(in) :: taylor2

!

taylor1%ref = taylor2%ref

if (associated(taylor2%term)) then
  call init_taylor_series (taylor1, size(taylor2%term))
  taylor1%term = taylor2%term
  taylor1%ref = taylor2%ref
else
  if (associated (taylor1%term)) deallocate (taylor1%term)
endif

end subroutine taylor_equal_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine taylors_equal_taylors (taylor1, taylor2)
!
! Subroutine to transfer the values from one taylor map to another:
!     Taylor1 <= Taylor2
!
! Input:
!   taylor2(:) -- Taylor_struct: Taylor map.
!
! Output:
!   taylor1(:) -- Taylor_struct: Taylor map. 
!-

subroutine taylors_equal_taylors (taylor1, taylor2)

implicit none

type (taylor_struct), intent(inout) :: taylor1(:)
type (taylor_struct), intent(in)    :: taylor2(:)

integer i

!

do i = 1, size(taylor1)
  taylor1(i) = taylor2(i)
enddo

end subroutine taylors_equal_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_taylor_series (bmad_taylor, n_term, save_old)
!
! Subroutine to initialize or extend a Bmad Taylor series (6 of these series make
! a Taylor map). Note: This routine does not zero the terms.
!
! Input:
!   bmad_taylor -- taylor_struct: Old structure.
!   n_term      -- integer: Number of terms to allocate. 
!                    n_term < 0 => bmad_taylor%term pointer will be disassociated.
!   save_old    -- logical, optional: If True then save any old terms and ref orbit when
!                    bmad_taylor is resized. If False zero the ref orbit. Default is False.
!
! Output:
!   bmad_taylor -- Taylor_struct: Initalized structure.
!-

subroutine init_taylor_series (bmad_taylor, n_term, save_old)

implicit none

type (taylor_struct) bmad_taylor
type (taylor_term_struct), pointer :: term(:)
integer n_term
integer n
logical, optional :: save_old

!

if (.not. logic_option (.false., save_old)) bmad_taylor%ref = 0

if (n_term < 0) then
  if (associated(bmad_taylor%term)) deallocate(bmad_taylor%term)
  return
endif

if (.not. associated (bmad_taylor%term)) then
  allocate (bmad_taylor%term(n_term))
  return
endif

if (size(bmad_taylor%term) == n_term) return

!

if (logic_option (.false., save_old) .and. n_term > 0 .and. size(bmad_taylor%term) > 0) then
  n = min (n_term, size(bmad_taylor%term))
  term => bmad_taylor%term
  allocate (bmad_taylor%term(n_term))
  bmad_taylor%term(1:n) = term(1:n)
  deallocate (term)

else
  deallocate (bmad_taylor%term)
  allocate (bmad_taylor%term(n_term))
endif

end subroutine init_taylor_series

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine em_taylor_equal_em_taylor (em_taylor1, em_taylor2)
!
! Subroutine that is used to set one em_taylor equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		em_taylor1 = em_taylor2
!
! Input:
!   em_taylor2 -- Em_taylor_struct: Input em_taylor.
!
! Output:
!   em_taylor1 -- Em_taylor_struct: Output em_taylor.
!-

subroutine em_taylor_equal_em_taylor (em_taylor1, em_taylor2)

implicit none
	
type (em_taylor_struct), intent(inout) :: em_taylor1
type (em_taylor_struct), intent(in) :: em_taylor2
integer n2

!

em_taylor1%ref = em_taylor2%ref

if (allocated(em_taylor2%term)) then
  n2 = size(em_taylor2%term)
  if (allocated(em_taylor1%term)) then
    if (size(em_taylor1%term) /= n2) then
      deallocate(em_taylor1%term)
      allocate (em_taylor1%term(n2))
    endif
  else
    allocate (em_taylor1%term(n2))
  endif
  em_taylor1%term = em_taylor2%term

else
  if (allocated(em_taylor1%term)) deallocate (em_taylor1%term)
endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine em_taylors_equal_em_taylors (em_taylor1, em_taylor2)
!
! Subroutine to transfer the values from one em_taylor map to another:
!     Em_taylor1 <= Em_taylor2
!
! Input:
!   em_taylor2(:) -- Em_taylor_struct: Em_taylor map.
!
! Output:
!   em_taylor1(:) -- Em_taylor_struct: Em_taylor map. 
!-

subroutine em_taylors_equal_em_taylors (em_taylor1, em_taylor2)

implicit none

type (em_taylor_struct), intent(inout) :: em_taylor1(:)
type (em_taylor_struct), intent(in)    :: em_taylor2(:)

integer i

!

do i = 1, size(em_taylor1)
  em_taylor1(i) = em_taylor2(i)
enddo

end subroutine em_taylors_equal_em_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_em_taylor_series (em_taylor, n_term, save_old)
!
! Subroutine to initialize a Bmad Em_taylor series (6 of these series make
! a Em_taylor map). Note: This routine does not zero the structure. The calling
! routine is responsible for setting all values.
!
! Input:
!   em_taylor   -- Em_taylor_struct: Old structure.
!   n_term      -- Integer: Number of terms to allocate. 
!                   n_term < 0 => em_taylor%term pointer will be disassociated.
!   save_old    -- Logical, optional: If True then save any old terms when
!                   em_taylor is resized. Default is False.
!
! Output:
!   em_taylor -- Em_taylor_struct: Initalized structure.
!-

subroutine init_em_taylor_series (em_taylor, n_term, save_old)

implicit none

type (em_taylor_struct) em_taylor
type (em_taylor_term_struct), allocatable :: term(:)
integer n_term
integer n
logical, optional :: save_old

!

if (n_term < 0) then
  if (allocated(em_taylor%term)) deallocate(em_taylor%term)
  return
endif

if (.not. allocated (em_taylor%term)) then
  allocate (em_taylor%term(n_term))
  return
endif

if (size(em_taylor%term) == n_term) return

!

if (logic_option (.false., save_old) .and. n_term > 0 .and. size(em_taylor%term) > 0) then
  n = min (n_term, size(em_taylor%term))
  call move_alloc(em_taylor%term, term)
  allocate (em_taylor%term(n_term))
  em_taylor%term(1:n) = term(1:n)
  deallocate (term)

else
  deallocate (em_taylor%term)
  allocate (em_taylor%term(n_term))
endif

end subroutine init_em_taylor_series

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)
!
! Subroutine that is used to set one complex_taylor equal to another. 
! This routine takes care of the pointers in complex_taylor1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		complex_taylor1 = complex_taylor2
!
! Input:
!   complex_taylor2 -- complex_taylor_struct: Input complex_taylor.
!
! Output:
!   complex_taylor1 -- complex_taylor_struct: Output complex_taylor.
!-

subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)

implicit none
	
type (complex_taylor_struct), intent(inout) :: complex_taylor1
type (complex_taylor_struct), intent(in) :: complex_taylor2

!

complex_taylor1%ref = complex_taylor2%ref

if (associated(complex_taylor2%term)) then
  call init_complex_taylor_series (complex_taylor1, size(complex_taylor2%term))
  complex_taylor1%term = complex_taylor2%term
else
  if (associated (complex_taylor1%term)) deallocate (complex_taylor1%term)
endif

end subroutine complex_taylor_equal_complex_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)
!
! Subroutine to transfer the values from one complex_taylor map to another:
!     complex_taylor1 <= complex_taylor2
!
! Input:
!   complex_taylor2(:) -- complex_taylor_struct: complex_taylor map.
!
! Output:
!   complex_taylor1(:) -- complex_taylor_struct: complex_taylor map. 
!-

subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)

implicit none

type (complex_taylor_struct), intent(inout) :: complex_taylor1(:)
type (complex_taylor_struct), intent(in)    :: complex_taylor2(:)

integer i

!

do i = 1, size(complex_taylor1)
  complex_taylor1(i) = complex_taylor2(i)
enddo

end subroutine complex_taylors_equal_complex_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_complex_taylor_series (complex_taylor, n_term, save)
!
! Subroutine to initialize a Bmad complex_taylor series (6 of these series make
! a complex_taylor map). Note: This routine does not zero the structure. The calling
! routine is responsible for setting all values.
!
! Input:
!   complex_taylor -- complex_taylor_struct: Old structure.
!   n_term      -- Integer: Number of terms to allocate. 
!                   n_term < 1 => complex_taylor%term pointer will be disassociated.
!   save        -- Logical, optional: If True then save any old terms when
!                   complex_taylor is resized. Default is False.
!
! Output:
!   complex_taylor -- complex_taylor_struct: Initalized structure.
!-

subroutine init_complex_taylor_series (complex_taylor, n_term, save)

implicit none

type (complex_taylor_struct) complex_taylor
type (complex_taylor_term_struct), allocatable :: term(:)
integer n_term
integer n
logical, optional :: save

!

if (n_term < 1) then
  if (associated(complex_taylor%term)) deallocate(complex_taylor%term)
  return
endif

if (.not. associated (complex_taylor%term)) then
  allocate (complex_taylor%term(n_term))
  return
endif

if (size(complex_taylor%term) == n_term) return

!

if (logic_option (.false., save) .and. n_term > 0 .and. size(complex_taylor%term) > 0) then
  n = min (n_term, size(complex_taylor%term))
  allocate (term(n))
  term = complex_taylor%term(1:n)
  deallocate (complex_taylor%term)
  allocate (complex_taylor%term(n_term))
  complex_taylor%term(1:n) = term
  deallocate (term)

else
  deallocate (complex_taylor%term)
  allocate (complex_taylor%term(n_term))
endif

end subroutine init_complex_taylor_series

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine bunch_equal_bunch (bunch1, bunch2)
!
! Subroutine to set one particle bunch equal to another.
!
! Note: This subroutine is called by the overloaded equal sign:
!    bunch1 = bunch2
!
! Input: 
!   bunch2 -- bunch_struct: Input bunch
!
! Output
!   bunch1 -- bunch_struct: Output bunch
!-

subroutine bunch_equal_bunch (bunch1, bunch2)

implicit none

type (bunch_struct), intent(inout) :: bunch1
type (bunch_struct), intent(in)    :: bunch2

integer i, np

!

if (.not. allocated(bunch2%particle)) then
  if (allocated(bunch1%particle)) deallocate (bunch1%particle, bunch1%ix_z)
else
  np = size(bunch2%particle)
  if (.not. allocated(bunch1%particle)) allocate(bunch1%particle(np), bunch1%ix_z(np))
  if (size(bunch1%particle) /= size(bunch2%particle)) then
    deallocate (bunch1%particle, bunch1%ix_z)
    allocate (bunch1%particle(np), bunch1%ix_z(np))
  endif
  bunch1%particle = bunch2%particle
  bunch1%ix_z     = bunch2%ix_z
endif

bunch1%charge_tot  = bunch2%charge_tot
bunch1%charge_live = bunch2%charge_live
bunch1%t_center    = bunch2%t_center
bunch1%ix_ele      = bunch2%ix_ele
bunch1%ix_bunch    = bunch2%ix_bunch
bunch1%n_live      = bunch2%n_live

end subroutine bunch_equal_bunch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine beam_equal_beam (beam1, beam2)
!
! Subroutine to set one particle beam equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!    beam1 = beam2
!
! Input: 
!  beam2 -- beam_struct: Input beam
!
! Output
!  beam1 -- beam_struct: Output beam
!
!-

subroutine beam_equal_beam (beam1, beam2)

implicit none

type (beam_struct), intent(inout) :: beam1
type (beam_struct), intent(in)    :: beam2

integer i, j, n_bun, n_particle
logical allocate_this

! The following rule must be observed: If beam%bunch is allocated then
! beam%bunch%particle must be also.

n_bun = size(beam2%bunch)

allocate_this = .true.
if (allocated(beam1%bunch)) then
  if (size(beam1%bunch) /= size(beam2%bunch)) then
    do i = 1, size(beam1%bunch)
      deallocate (beam1%bunch(i)%particle)
    enddo
    deallocate (beam1%bunch)
  else
    allocate_this = .false.
  endif
endif

if (allocate_this) allocate (beam1%bunch(n_bun))

do i = 1, n_bun
  beam1%bunch(i) = beam2%bunch(i)
enddo

end subroutine beam_equal_beam

end module

