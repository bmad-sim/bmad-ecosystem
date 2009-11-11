module transfer_map_mod

use bmad_struct
use bmad_interface

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine transfer_map_calc_at_s (lat, t_map, s1, s2, &
!                                         integrate, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between longitudinal positions
! s1 to s2.
!
! If s2 < s1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will 'wrap around' the lattice end.
! For example, if s1 = 900 and s2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
!
! If s2 = s1 then you get the unit map except if one_turn = True.
!
! Note: If integrate = False and if a taylor map does not exist for an 
! element this routine will make one and store it in the element.
!
! Modules Needed:
!   use transfer_map_mod
!
! Input:
!   lat        -- Lat_struct: Lattice used in the calculation.
!   t_map(6)   -- Taylor_struct: Initial map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start position for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end position for the calculation.
!                   Default is lat%param%total_length.
!   integrate  -- Logical, optional: If present and True then do symplectic
!                   integration instead of concatenation. 
!                   Default = False.
!   one_turn   -- Logical, optional: If present and True, and s1 = s2, and the lattice
!                   is circular: Construct the one-turn map from s1 back to s1.
!                   Otherwise t_map is unchanged or the unit map if unit_start = T.
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then t_map will be
!                   used as the starting map instead of the unit map.
!                   Default = True
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

subroutine transfer_map_calc_at_s (lat, t_map, s1, s2, &
                                      integrate, one_turn, unit_start)

use ptc_interface_mod, only: concat_taylor, ele_to_taylor, taylor_propagate1, taylor_inverse
use bookkeeper_mod, only: makeup_super_slave1, attribute_bookkeeper

implicit none

type (lat_struct) lat
type (taylor_struct) :: t_map(:)
type (taylor_struct) a_map(6)

real(rp), intent(in), optional :: s1, s2
real(rp) ss1, ss2

logical, optional :: integrate, one_turn, unit_start
logical integrate_this, one_turn_this, unit_start_this

character(40) :: r_name = 'transfer_map_calc_at_s'

!

integrate_this  = logic_option (.false., integrate)
one_turn_this   = logic_option (.false., one_turn)
unit_start_this = logic_option(.true., unit_start)

ss1 = 0;                       if (present(s1)) ss1 = s1
ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
 
if (unit_start_this) call taylor_make_unit (t_map)

! One turn or not calc?

if (ss1 == ss2 .and. (.not. one_turn_this .or. &
                lat%param%lattice_type == linear_lattice$)) return

! Normal case

if (ss1 < ss2) then 
  call transfer_this (t_map, ss1, ss2)

! For a circular lattice push through the origin.

elseif (lat%param%lattice_type == circular_lattice$) then
  call transfer_this (t_map, ss1, lat%param%total_length)
  call transfer_this (t_map, 0.0_rp, ss2)

! For a linear lattice compute the backwards matrix

else
  if (unit_start_this) then
    call transfer_this (t_map, ss2, ss1)
    call taylor_inverse (t_map, t_map)
  else
    call taylor_make_unit (a_map)
    call transfer_this (a_map, ss2, ss1)
    call taylor_inverse (a_map, a_map)
    call concat_taylor (t_map, a_map, t_map)
    call kill_taylor (a_map)
  endif

endif

!------------------------------------------------------------------------
contains

subroutine transfer_this (map, s_1, s_2)

type (taylor_struct) :: map(:)
type (ele_struct), pointer, save :: ele
type (ele_struct), pointer, save :: runt
type (ele_struct), target, save :: runt_save

real(rp) s_1, s_2, s_now, s_end, ds
real(rp), save :: ds_old = -1

integer i, ix_ele
integer, save :: ix_ele_old = -1

logical kill_it, track_entrance, track_exit, track_entire_ele
logical, save :: old_track_end = .false.

! Init

call ele_at_s (lat, s_1, ix_ele)
s_now = s_1

! Loop over all the element to track through.

do

  ele => lat%ele(ix_ele)
  s_end = min(s_2, ele%s)

  track_entrance   = (s_now == lat%ele(ix_ele-1)%s) 
  track_exit       = (s_end == ele%s)
  track_entire_ele = (track_entrance .and. track_exit)

  ds = s_end - s_now

  ! runt points at ele if we are tracking through the entire element and
  ! at runt_save if only a partial track. We do this since we will mangle
  ! the element with a partial track.

  if (ix_ele /= ix_ele_old) then
    if (track_entire_ele) then
      runt => ele
    else  ! partial track
      runt_save = ele
      runt => runt_save
      old_track_end = .false.
    endif
  endif

  ! We only need to do the "split" bookkeeping if we are only partially tracking
  ! through the element and we have not done the bookkeeping before.

  if (.not. track_entire_ele) then

    ! Kill the saved taylor map if it does not apply to the present integration step.

    kill_it = .false.
    if (ds /= ds_old .or. ix_ele /= ix_ele_old) then
      kill_it = .true.
    elseif (ele%key == sbend$) then
      if (track_entrance .or. track_exit .or. old_track_end) kill_it = .true.
    elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then
      kill_it = .true.
    elseif (ele%key == lcavity$) then
      kill_it = .true.
    endif

    if (kill_it) then
      call kill_taylor (runt%taylor)
      runt%value(l$) = ds
      call makeup_super_slave1 (runt, ele, s_now-lat%ele(ix_ele-1)%s, &
                                              lat%param, track_entrance, track_exit)
      call attribute_bookkeeper (runt, lat%param)
    endif

  endif

  ! Now for the integration step

  if (integrate_this) then
    call taylor_propagate1 (map, runt, lat%param)
  else
    if (.not. associated(runt%taylor(1)%term)) then
      call ele_to_taylor (runt, lat%param)
    endif
    call concat_taylor (map, runt%taylor, map)
  endif

  ! Save the present integration step parameters so that if this routine
  ! is called in the future we can tell if the saved taylor map is still valid.

  ds_old = ds
  ix_ele_old = ix_ele
  old_track_end = track_entrance .or. track_exit

  ! Are we done?

  if (s_end == s_2) return

  ! We are not done so move to the next element.

  s_now = s_end
  ix_ele = ix_ele + 1

enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine mat6_calc_at_s (lat, mat6, vec0, s1, s2, one_turn, unit_start)
!
! Subroutine to calculate the transfer map between longitudinal positions
! s1 to s2.
!
! If s2 < s1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will "wrap around" the lattice end.
! For example, if s1 = 900 and s2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
!
! If s2 = s1 then you get the unit matrix except if one_turn = True.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Lat_struct: Lattice used in the calculation.
!   mat6(6,6)  -- Real(rp): Initial matrix (used when unit_start = False)
!   vec0(6)    -- Real(rp): Initial 0th order map (used when unit_start = False)
!   s1         -- Real(rp), optional: Element start index for the calculation.
!                   Default is 0.
!   s2         -- Real(rp), optional: Element end index for the calculation.
!                   Default is lat%param%total_length.
!   one_turn   -- Logical, optional: If present and True, and s1 = s2, and the lattice
!                   is circular: Construct the one-turn matrix from s1 back to s1.
!                   Otherwise mat6 is unchanged or the unit map if unit_start = T.
!                   Default = False.
!   unit_start -- Logical, optional: If present and False then mat6 will be
!                   used as the starting matrix instead of the unit matrix.
!                   Default = True
!
! Output:
!    mat6(6,6) -- Real(rp): Transfer matrix.
!    vec0(6)   -- Real(rp): 0th order part of the map.
!-

subroutine mat6_calc_at_s (lat, mat6, vec0, s1, s2, one_turn, unit_start)

use bookkeeper_mod, only: makeup_super_slave1, attribute_bookkeeper

implicit none

type (lat_struct) lat

real(rp) mat6(:,:), vec0(:)
real(rp), intent(in), optional :: s1, s2
real(rp) ss1, ss2

logical, optional :: one_turn, unit_start
logical one_turn_this

!

one_turn_this = logic_option (.false., one_turn)

ss1 = 0;                       if (present(s1)) ss1 = s1
ss2 = lat%param%total_length;  if (present(s2)) ss2 = s2
 
if (logic_option(.true., unit_start)) then
  call mat_make_unit (mat6)
  vec0 = 0
endif

! One turn or not calc?

if (ss1 == ss2 .and. (.not. one_turn_this .or. &
                lat%param%lattice_type == linear_lattice$)) return

! Normal case

if (ss1 < ss2) then
  call transfer_this (ss1, ss2)

! For a circular lattice push through the origin.

elseif (lat%param%lattice_type == circular_lattice$) then
  call transfer_this (ss1, lat%param%total_length)
  call transfer_this (0.0_rp, ss2)

! For a linear lattice compute the backwards matrix

else
  call transfer_this (ss2, ss1)
  call mat_inverse (mat6, mat6)
  vec0 = -matmul(mat6, vec0)

endif


!--------------------------------------------------------

contains

subroutine transfer_this (s_1, s_2)

type (ele_struct), pointer, save :: ele
type (ele_struct), pointer, save :: runt
type (ele_struct), target, save :: runt_save

real(rp) s_1, s_2, s_end, s_now, ds
real(rp), save :: ds_old = -1

integer ix_ele
integer, save :: ix_ele_old = -1

logical track_entrance, track_exit, track_entire_ele, kill_it
logical, save :: old_track_end = .false.

! Init

call ele_at_s (lat, s_1, ix_ele)
s_now = s_1

! Loop over all the element to track through.

do

  ele => lat%ele(ix_ele)
  s_end = min(s_2, ele%s)

  track_entrance   = (s_now == lat%ele(ix_ele-1)%s) 
  track_exit       = (s_end == ele%s)
  track_entire_ele = (track_entrance .and. track_exit)

  ds = s_end - s_now

  ! runt points at ele if we are tracking through the entire element and
  ! at runt_save if only a partial track. We do this since we will mangle
  ! the element with a partial track.

  if (ix_ele /= ix_ele_old) then
    if (track_entire_ele) then
      runt => ele
    else  ! partial track
      runt_save = ele
      runt => runt_save
      old_track_end = .false.
    endif
  endif

  ! We only need to do the "split" bookkeeping if we are only partially tracking
  ! through the element and we have not done the bookkeeping before.

  if (.not. track_entire_ele) then

    ! Kill the saved matrix if it does not apply to the present integration step.

    kill_it = .false.
    if (ds /= ds_old .or. ix_ele /= ix_ele_old) then
      kill_it = .true.
    elseif (ele%key == sbend$) then
      if (track_entrance .or. track_exit .or. old_track_end) kill_it = .true.
    elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then
      kill_it = .true.
    elseif (ele%key == lcavity$) then
      kill_it = .true.
    endif

    if (kill_it) then
      runt%value(l$) = ds
      call makeup_super_slave1 (runt, ele, s_now-lat%ele(ix_ele-1)%s, &
                                            lat%param, track_entrance, track_exit)
      call attribute_bookkeeper (runt, lat%param)

      call make_mat6 (runt, lat%param)
    endif

  endif

  ! Now for the integration step

  mat6 = matmul (runt%mat6, mat6)
  vec0 = matmul (runt%mat6, vec0) + runt%vec0

  ! Save the present integration step parameters so that if this routine
  ! is called in the future we can tell if the saved taylor map is still valid.

  ds_old = ds
  ix_ele_old = ix_ele
  old_track_end = track_entrance .or. track_exit

  ! Are we done?

  if (s_end == s_2) return

  ! We are not done so move to the next element.

  s_now = s_end
  ix_ele = ix_ele + 1

enddo

end subroutine

end subroutine

end module
