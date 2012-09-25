module bookkeeper_mod

use wall3d_mod
use bmad_utils_mod
use multipole_mod
use lat_geometry_mod
use equality_mod
use em_field_mod
use crystal_param_mod

integer, parameter :: off$ = 1, on$ = 2
integer, parameter :: save_state$ = 3, restore_state$ = 4

private control_bookkeeper1, makeup_overlay_and_girder_slave, super_lord_length_bookkeeper 
private makeup_group_lord, makeup_super_slave1, makeup_super_slave
private compute_slave_aperture 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_attribute (lat, ele, attrib)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
! Also will do some dependent variable bookkeeping when a particular attribute has 
! been altered. Look at this routine's code for more details.
!
! set_flages_for_changed_attribute is an overloaded name for:
!   set_flages_for_changed_real_attribute 
!   set_flages_for_changed_inteter_attribute 
!   set_flages_for_changed_logical_attribute 
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: Lattice with the changed attribute.
!   ele    -- ele_struct, optional: Element being modified.
!               If not present, mark the entire lattice as being modified.
!   attrib -- Real(rp), integer, or logical; optional: Attribute that has been changed.
!
! Output:
!   lat  -- lat_struct: Lattice with appropriate changes.
!-

interface set_flags_for_changed_attribute
  module procedure set_flags_for_changed_real_attribute 
  module procedure set_flags_for_changed_integer_attribute 
  module procedure set_flags_for_changed_logical_attribute 
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lattice_bookkeeper (lat, err_flag)
!
! Subroutine to do a "complete" bookkeeping job on a lattice:
!   lord/slave control
!   reference energy calc
!   s-position calc
!   geometry (floor position) calc
!
! Not done are Twiss, transfer matrices, and orbit calculations.  
!
! Note: This this routine does a complete job of bookking
! and could be unacceptably slow if lat%auto_bookkeeper = True.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice needing bookkeeping.
!
! Output:
!   lat      -- lat_struct: Lattice with bookkeeping done.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine lattice_bookkeeper (lat, err_flag)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (bookkeeping_state_struct), pointer :: stat

integer i, j

logical, optional :: err_flag
logical found, err

character(20), parameter :: r_name = 'lattice_bookkeeper'

! Control bookkeeper is called twice to make sure, for example, that multipass bend
! are correctly computed.

if (present(err_flag)) err_flag = .true.

call control_bookkeeper (lat, do_free_eles = .true.)
call lat_compute_ref_energy_and_time (lat, err)
if (err) return
call control_bookkeeper (lat, super_and_multipass_only = .true., do_free_eles = .true.)

! Global geometry

call s_calc (lat)
call lat_geometry (lat)

! multipass slaves with ref_orbit set may depend upon the geometry so recalc.
! Also free elements may have had their control status set but this is bogus so just reset to ok$

found = .false.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  do j = 1, branch%n_ele_track
    ele => branch%ele(j)
    if (ele%slave_status == free$ .and. ele%bookkeeping_state%control == stale$) ele%bookkeeping_state%control = ok$
    if (.not. bmad_com%auto_bookkeeper .and. ele%bookkeeping_state%control /= stale$) cycle
    if (ele%slave_status == multipass_slave$ .and. ele%ref_orbit /= 0) then
      call makeup_multipass_slave (lat, ele)
      call attribute_bookkeeper (ele, branch%param)
      found = .true.
    endif
  enddo
  branch%param%bookkeeping_state%attributes = ok$
enddo

if (found) then
  call s_calc (lat)
  call lat_geometry (lat)
endif

! See if all status flags have been properly reset.
! Exception is mat6 flag since the bookkeeping routines do not touch this.

if (.not. bmad_com%auto_bookkeeper) then

  stat => lat%lord_state
  if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
      stat%s_position == stale$ .or. stat%ref_energy == stale$) then
    call out_io (s_info$, r_name, 'Stale bookkeeping lord_status flags detected.', &
                                  'Please contact DCS!', 'Status: \5i6\ ', &
            i_array = [stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
  endif
  call reset_status_flags_to_ok(stat)

  do i = 0, ubound(lat%branch, 1)

    branch => lat%branch(0)
    stat => branch%param%bookkeeping_state
    if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
        stat%s_position == stale$ .or. stat%ref_energy == stale$) then
      call out_io (s_info$, r_name, 'Stale bookkeeping status flags detected at: \i0\.', &
                                    'Please contact DCS!', 'Status: \5i6\ ', &
              i_array = [i, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
    endif
    call reset_status_flags_to_ok(stat)

    do j = 0, branch%n_ele_max
      if (branch%ele(j)%key == null_ele$) cycle 
      stat => branch%ele(j)%bookkeeping_state
      if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
          stat%s_position == stale$ .or. stat%ref_energy == stale$) then
        call out_io (s_info$, r_name, 'Stale bookkeeping status flags detected at: \i0\, \i0\.', &
                                      'Please contact DCS!', 'Status: \5i6\ ', &
              i_array = [i, j, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
      endif
      call reset_status_flags_to_ok(stat)
    enddo

  enddo
endif

if (present(err_flag)) err_flag = .false.

!----------------------------------------------------------
contains

subroutine reset_status_flags_to_ok (stat)
  type (bookkeeping_state_struct) stat

  if (stat%control /= ok$        .and. stat%control /= super_ok$)        stat%control = ok$
  if (stat%attributes /= ok$     .and. stat%attributes /= super_ok$)     stat%attributes = ok$
  if (stat%floor_position /= ok$ .and. stat%floor_position /= super_ok$) stat%floor_position = ok$
  if (stat%s_position /= ok$         .and. stat%s_position /= super_ok$)         stat%s_position = ok$
  if (stat%ref_energy /= ok$     .and. stat%ref_energy /= super_ok$)     stat%ref_energy = ok$

end subroutine reset_status_flags_to_ok

end subroutine lattice_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (lat, ele, super_and_multipass_only, do_free_eles)
!
! Subroutine to transfer attibute information from lord to slave elements.
! This subroutine will call attribute_bookkeeper.
! Note: To do a complete bookkeeping job on a lattice use:
!   lattice_bookkeeper
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: lattice to be used
!   ele    -- Ele_struct, optional: Element whose attribute values 
!               have been changed. If not present bookkeeping will be done 
!               for all elements.
!   super_and_multipass_only 
!          -- Logical, optional: Default is False. If True then only do bookkeeping 
!               for superposition and multipass elements only. Default is False. 
!               This argument is used by lattice_bookkeeper and should not
!               to set unless you know what you are doing.
!   do_free_eles
!         -- Logical, optional: Default is False. If True then call 
!               attribute_bookkeeper for elements in the tracking part of the lattice.
!-

recursive subroutine control_bookkeeper (lat, ele, super_and_multipass_only, do_free_eles)

implicit none

type (lat_struct), target :: lat
type (ele_struct), optional :: ele
type (ele_struct), pointer :: slave, lord, branch_ele, ele2
type (branch_struct), pointer :: branch

integer i, j, ie, ib, n1, n2

logical, optional :: super_and_multipass_only, do_free_eles
logical sm_only, all_bookkeeping_done

character(20), parameter :: r_name = 'control_bookkeeper'

! Check that super_slave lengths add up to the super_lord_length.
! When super_and_multipass_only is present, this check has already been
! done so we don't need to do it again.

sm_only = logic_option (.false., super_and_multipass_only)

if (.not. sm_only) then
  if (present(ele)) then
    if (ele%lord_status == super_lord$) call super_lord_length_bookkeeper (lat, ele) 
  else
    call super_lord_length_bookkeeper (lat) 
  endif
endif

! If ele is present we only do bookkeeping for this one element and its slaves

if (present(ele)) then
  call control_bookkeeper1 (lat, ele, sm_only)
  do ie = 1, ele%n_slave
    slave => pointer_to_slave(ele, ie)
    call control_bookkeeper (lat, slave, sm_only)
  enddo
  return
endif

! Else we need to make up all the lords...
! First mark all the elements needing bookkeeping

n1 = lat%n_ele_track+1
n2 = lat%n_ele_max

if (bmad_com%auto_bookkeeper) then
  lat%ele(n1:n2)%bookkeeping_state%control = stale$  ! Bookkeeping done on this element yet?
endif

! Now do the control bookkeeping.
! Need to do this from the top level down.
! The top level are those lord elements that have no lords.

do
  all_bookkeeping_done = .true.
  ie_loop: do ie = n1, n2
    ele2 => lat%ele(ie)
    if (ele2%key == null_ele$) cycle
    if (ele2%bookkeeping_state%control /= stale$ .and. ele2%bookkeeping_state%attributes /= stale$) cycle
    do j = 1, ele2%n_lord
      lord => pointer_to_lord(ele2, j)
      if (lord%bookkeeping_state%control /= stale$ .and. lord%bookkeeping_state%attributes /= stale$) cycle
      all_bookkeeping_done = .false.  ! This element remains to be done.
      cycle ie_loop ! Do not do bookkeeping yet if lord not done yet.
    enddo
    call control_bookkeeper1 (lat, ele2, sm_only)
    ele2%bookkeeping_state%control = ok$  ! Done with this element
  enddo ie_loop
  if (all_bookkeeping_done) exit  ! And we are done
enddo

! and now the slaves in the tracking lattice

do ib = 0, ubound(lat%branch, 1)
  if (.not. bmad_com%auto_bookkeeper .and. lat%branch(ib)%param%bookkeeping_state%control /= stale$) cycle
  do ie = 0, lat%branch(ib)%n_ele_track
    ele2 => lat%branch(ib)%ele(ie)
    if (ele2%key == null_ele$) cycle
    if (ele2%slave_status == free$) cycle
    if (.not. bmad_com%auto_bookkeeper .and. ele2%bookkeeping_state%control /= stale$) cycle
    call control_bookkeeper1 (lat, ele2, sm_only)
    ele2%bookkeeping_state%control = ok$
  enddo
  lat%branch(ib)%param%bookkeeping_state%control = ok$
enddo

! Update attributes for elements in the tracking part of the lattice.

if (logic_option(.false., do_free_eles)) then
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    if (.not. bmad_com%auto_bookkeeper .and. branch%param%bookkeeping_state%attributes /= stale$) cycle
    do j = 0, branch%n_ele_track
      call attribute_bookkeeper (branch%ele(j), branch%param)
    enddo
    branch%param%bookkeeping_state%attributes = ok$
  enddo
endif

lat%lord_state%control = ok$
lat%lord_state%attributes = ok$

end subroutine control_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper1 (lat, ele, sm_only)
!
! This routine is for control bookkeeping for a single element.
! This subroutine is only to be called from control_bookkeeper and is
! not meant for general use.
!-

subroutine control_bookkeeper1 (lat, ele, sm_only)

type (lat_struct), target :: lat
type (ele_struct) ele

logical sm_only, call_a_bookkeeper

! Init

if (sm_only) then
  if (ele%lord_status /= super_lord$ .and. ele%slave_status /= super_slave$ .and. &
      ele%lord_status /= multipass_lord$ .and. ele%slave_status /= multipass_slave$) return
endif

! First make sure the attribute bookkeeping for this element is correct since
! the makeup_*_slave routines may need it.

if (ele%key /= overlay$ .and. ele%key /= group$) call attribute_bookkeeper (ele, lat%branch(ele%ix_branch)%param)

! Slave bookkeeping

call_a_bookkeeper = .false.

if (ele%slave_status == super_slave$) then
  ! Attrubute bookkeeping is done in the makeup_super_slave
  call makeup_super_slave (lat, ele)

elseif (ele%slave_status == overlay_slave$) then
  call makeup_overlay_and_girder_slave (lat, ele)
  call_a_bookkeeper = .true.

elseif (ele%slave_status == multipass_slave$) then
  call makeup_multipass_slave (lat, ele)
  if (ele%n_lord > 1) call makeup_overlay_and_girder_slave (lat, ele)
  call_a_bookkeeper = .true.
endif

! Lord bookkeeping

if (ele%lord_status == group_lord$) then
  call makeup_group_lord (lat, ele)
  call_a_bookkeeper = .true.

elseif (ele%lord_status == super_lord$) then
  call adjust_super_lord_s_position (lat, ele)
  call_a_bookkeeper = .true.
endif

! If bookkeeping has been done by a makeup_*_salve routine then
! attribute_bookkeeper must be called again.
! This is true even if the lattice is static since a slave element
! can have its lord's dependent attribute values.
! Example: super_slave will, at this point, have its lord's num_steps value but 
! num_steps in the slave is different from the lord due to differences in length.

if (call_a_bookkeeper) call attribute_bookkeeper (ele, lat%branch(ele%ix_branch)%param)

end subroutine control_bookkeeper1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine super_lord_length_bookkeeper (lat, ele)
!
! Subroutine to make sure the length of the slaves of a super_lord add up to the
! length of the lord. If not, make an adjustment to the slave length.
!
! Note: This routine is private and is called by control_bookkeeper. 
!
! Modules needed:
!   use bookkeeper_mod
!
! Input:
!   lat   -- Lat_struct: Lattice.
!   ele   -- Ele_struct, optional: Index of super_lord element to check.
!                  If not present, bookkeeping will be done for all super_lords.
!
! Output:
!   lat  -- Lat_struct: Lattice with adjustments made.
!-

subroutine super_lord_length_bookkeeper (lat, ele)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: lord0, lord2, slave, slave2
type (branch_struct), pointer :: branch
type (ele_struct), optional :: ele

real(rp) sum_len_slaves, coef, vary_sublength, length_pos, length_neg
real(rp) d_length, d_length_pos, d_length_neg
real(rp) dl_tol

integer j, k, ie, ixa_lord0, ixb_lord0, ixa_lord2, ixb_lord2
integer ix_pos_edge, ix_neg_edge, ixa, ixb

logical pos_extension_lord_exists, neg_extension_lord_exists, all_extension_lord_exists
logical length_adjustment_made, overlap_a, overlap_b, overlap_all
logical, allocatable :: can_vary_length(:)

character(40) :: r_name = 'super_lord_length_bookkeeper'

!

if (.not. bmad_com%auto_bookkeeper) then
  if (lat%branch(0)%param%bookkeeping_state%s_position /= stale$) return
endif

dl_tol = 10 * bmad_com%significant_length

length_adjustment_made = .false.

do ie = lat%n_ele_track+1, lat%n_ele_max

  lord0 => lat%ele(ie)

  if (lord0%key == null_ele$) cycle
  if (lord0%lord_status /= super_lord$) cycle

  if (present(ele)) then
    if (ele%ix_ele /= ie) cycle
  endif

  if (.not. bmad_com%auto_bookkeeper .and. lord0%bookkeeping_state%s_position /= stale$) cycle

  sum_len_slaves = 0
  do j = 1, lord0%n_slave
    slave => pointer_to_slave(lord0, j)
    sum_len_slaves = sum_len_slaves + slave%value(l$)
  enddo

  ! Nothing to be done if the lengths add up.

  if (abs(sum_len_slaves - lord0%value(l$)) < dl_tol * (1 + lord0%value(l$))) cycle

  ! Now we need to adjust some super_slave lengths.
  ! We try varying the length of all the slaves except
  ! those that are part of a "contained" super_lord. A "contained" super_lord
  ! is a super_lord that the present lord (lord0) completely overlaps.
  ! This is necessary since otherwise the length of the contained super_lord
  ! would not be consistant with the lengths of its slaves.

  ! The complication here is that we will need to adjust the lengths
  ! of the elements to either side of the lord0 to keep the lengths of other
  ! super_lords consistant with their slaves. We need to know if these other
  ! super_lords extend in the positive, negative or both directions past lord0.

  length_adjustment_made = .true.

  slave => pointer_to_slave(lord0, 1)
  ixa_lord0 = slave%ix_ele  ! Index at entrance end of lord0

  slave => pointer_to_slave(lord0, lord0%n_slave)
  ixb_lord0 = slave%ix_ele  ! Index at exit end of lord0

  pos_extension_lord_exists = .false.
  neg_extension_lord_exists = .false.
  all_extension_lord_exists = .false.

  ix_pos_edge = lat%n_ele_max
  ix_neg_edge = 0

  vary_sublength = 0
  call re_allocate(can_vary_length, lord0%n_slave, .false.)

  slave_loop: do j = 1, lord0%n_slave

    slave => pointer_to_slave(lord0, j)
    can_vary_length(j) = .true. ! Can be varied.

    do k = 1, slave%n_lord
      lord2 => pointer_to_lord(slave, k)
      if (lord0%ix_ele == lord2%ix_ele) cycle  ! Ignore self
      slave2 => pointer_to_slave(lord2, 1)            ! Slave at entrance end
      ixa_lord2 = slave2%ix_ele

      slave2 => pointer_to_slave(lord2, lord2%n_slave) ! Slave at exit end
      ixb_lord2 = slave2%ix_ele

      overlap_a = .false. ! entrance end of lord2 overlaps lord0?
      overlap_b = .false. ! exit end of lord2 overlaps lord0?
      overlap_all = .false.

      ! Case where lord0 does not wrap around the IP

      if (ixa_lord0 <= ixb_lord0) then
        if (ixa_lord2 >= ixa_lord0 .and. ixa_lord2 <= ixb_lord0) overlap_a = .true.
        if (ixb_lord2 >= ixa_lord0 .and. ixb_lord2 <= ixb_lord0) overlap_b = .true.
        if (.not. overlap_a .and. .not. overlap_b) then
          if (ixa_lord2 <= ixb_lord2) then           ! If lord2 does not wrap
            if (ixa_lord2 < ixa_lord0 .and. ixb_lord2 > ixb_lord0) overlap_all = .true.
          else                                       ! If lord2 does wrap
            if (ixa_lord2 < ixa_lord0 .or. ixb_lord2 > ixb_lord0) overlap_all = .true.
          endif
        endif

      ! Case where lord0 does wrap around the IP

      else
        if (ixa_lord2 >= ixa_lord0 .or. ixa_lord2 <= ixb_lord0) overlap_a = .true.
        if (ixb_lord2 >= ixa_lord0 .or. ixb_lord2 <= ixb_lord0) overlap_b = .true.
        if (.not. overlap_a .and. .not. overlap_b) then
          if (ixa_lord2 > ixb_lord2) then ! and lord2 wraps also
            if (ixa_lord2 < ixa_lord0 .and. ixb_lord2 > ixb_lord0) overlap_all = .true.
          endif
        endif
      endif

      ! Contained?

      if (overlap_a .and. overlap_b) then  ! If contained
        can_vary_length(j) = .false.
      elseif (overlap_a) then
        pos_extension_lord_exists = .true.
        ix_pos_edge = min (ix_pos_edge, ixa_lord2)
      elseif (overlap_b) then
        neg_extension_lord_exists = .true.
        ix_neg_edge = max (ix_neg_edge, ixb_lord2)
      elseif (overlap_all) then
        all_extension_lord_exists = .true.
      endif

    enddo

    if (can_vary_length(j)) vary_sublength = vary_sublength + slave%value(l$)
  enddo slave_loop

  ! If we have not found any slaves to vary we are in trouble

  if (vary_sublength == 0) then
    call out_io (s_fatal$, r_name, 'CANNOT VARY LENGTH OF SUPER_LORD: ' // lord0%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ! Calculate positive and negative extension length changes

  coef = (lord0%value(l$) - sum_len_slaves) / vary_sublength
  d_length = lord0%value(l$) - sum_len_slaves

  if (pos_extension_lord_exists) then
    length_pos = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lord0, j)
      if (slave%ix_ele < ix_pos_edge) cycle
      if (.not. can_vary_length(j)) cycle
      length_pos = length_pos + slave%value(l$)
    enddo
    d_length_pos = length_pos * coef
  endif

  if (neg_extension_lord_exists) then
    length_neg = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lord0, j)
      if (slave%ix_ele > ix_neg_edge) cycle
      if (.not. can_vary_length(j)) cycle
      length_neg = length_neg + slave%value(l$)
    enddo    
    d_length_neg = length_neg * coef
  endif

  ! Vary the slave lengths

  do j = 1, lord0%n_slave
    slave => pointer_to_slave(lord0, j)
    if (.not. can_vary_length(j)) cycle
    slave%value(l$) = slave%value(l$) * (1 + coef)
    call set_ele_status_stale (slave, attribute_group$)
  enddo

  ! Now to make the adjustments to either side of lord0.

  branch => lat%branch(slave%ix_branch)

  ixa = ixa_lord0 - 1
  if (ixa == 0) ixa = branch%n_ele_track

  ixb = ixb_lord0 + 1
  if (ixb == branch%n_ele_track + 1) ixb = 1 

  if (all_extension_lord_exists) then
    if (pos_extension_lord_exists .and. neg_extension_lord_exists) then
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length_neg
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length_pos
    elseif (pos_extension_lord_exists) then
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - (d_length - d_length_pos)
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length_pos
    elseif (neg_extension_lord_exists) then
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length_neg
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - (d_length - d_length_neg)
    else
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length / 2
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length / 2
    endif

  else ! An all_extension_lord does not exist
    if (pos_extension_lord_exists) &
            branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length_pos
    if (neg_extension_lord_exists) &
            branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length_neg
  endif

  call set_ele_status_stale (branch%ele(ixa), attribute_group$)
  call set_ele_status_stale (branch%ele(ixb), attribute_group$)

enddo

! If there has been a length adjustment then we need to make sure everything is ok.

if (length_adjustment_made) then
  do ie = lat%n_ele_track+1, lat%n_ele_max
    lord0 => lat%ele(ie)
    if (lord0%lord_status /= super_lord$) cycle
    sum_len_slaves = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lord0, j)
      sum_len_slaves = sum_len_slaves + slave%value(l$)
    enddo
    if (abs(sum_len_slaves - lord0%value(l$)) > dl_tol * abs(lord0%value(l$))) then
      call out_io (s_fatal$, r_name, &
              'INCONSISTANT SUPER_LORD/SUPER_SLAVE LENGTHS!', &
              'LORD: ' // lord0%name, &
              'LENGTH: \es16.9\ ', &
              'SUM OF SLAVE LENGTHS: \es16.9\ ', r_array = [lord0%value(l$), sum_len_slaves] )
      if (global_com%exit_on_error) call err_exit
      return
    endif
  enddo
endif

end subroutine super_lord_length_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine adjust_super_lord_s_position (lat, lord)
!
! Subroutine to adjust the positions of the slaves of a super_lord due
! to changes in the lord's s_offset.
!-

subroutine adjust_super_lord_s_position (lat, lord)

implicit none

type (lat_struct), target :: lat
type (ele_struct) lord
type (ele_struct), pointer :: slave 

real(rp) s_start_lord, s_start_slave, s_end_lord, tot_len, sig_l

character(40) :: r_name = 'adjust_super_lord_s_position'

!

if (lord%lord_status /= super_lord$) then
  call out_io (s_abort$, r_name, 'ELEMENT IS NOT A LORD! ' // lord%name)
  if (global_com%exit_on_error) call err_exit
  return 
endif

! If a super lord is moved then we just need to adjust the start and end edges.
! Since we don't want to kill taylor maps due to round-off errors we only change
! the length if the percentage or absolute change is more than 10^-10

s_end_lord = lord%s + lord%value(s_offset_tot$)
slave => pointer_to_slave(lord, lord%n_slave)
sig_l = bmad_com%significant_length

if (abs(s_end_lord - slave%s) < sig_l * (1 + abs(slave%value(l$)))) return

slave%value(l$) = slave%value(l$) + (s_end_lord - slave%s)
slave%s = s_end_lord

! Adjust start position of the first slave

slave => pointer_to_slave(lord, 1)
tot_len = lat%branch(slave%ix_branch)%param%total_length
s_start_lord = s_end_lord - lord%value(l$)

s_start_slave = slave%s - slave%value(l$)
if (s_start_slave > s_start_lord + tot_len / 2) s_start_slave = s_start_slave - tot_len
slave%value(l$) = slave%value(l$) + (s_start_slave - s_start_lord)

end subroutine adjust_super_lord_s_position

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_lord (lat, lord)
!
! Subroutine to calculate the attributes of group slave elements.
! This routine is private to bookkeeper_mod.
!-

Subroutine makeup_group_lord (lat, lord)   

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: lord
type (ele_struct), pointer :: slave

real(rp) delta, coef
real(rp), pointer :: r_ptr

integer ix, iv, i

logical moved, err_flag

character(20) :: r_name = 'makeup_group_lord'

!

delta = lord%value(command$) - lord%value(old_command$)    ! change
lord%value(old_command$) = lord%value(command$) ! save old

moved = .false.   ! have we longitudinally moved an element?

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i, ix)
  iv = lat%control(ix)%ix_attrib
  if (iv == l$) moved = .true.
  coef = lat%control(ix)%coef
  call pointer_to_indexed_attribute (slave, iv, .false., r_ptr, err_flag)
  if (err_flag) call err_exit
  r_ptr = r_ptr + delta * coef
  call set_flags_for_changed_attribute (lat, slave, r_ptr)
enddo

if (moved) then
  call s_calc (lat)       ! recompute s distances
  call lat_geometry (lat)
endif

lord%bookkeeping_state%control = ok$

end subroutine makeup_group_lord

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_multipass_slave (lat, slave)
!
! Subroutine to calcualte the attributes of multipass slave elements.
! This routine is not meant for general use.
!-

subroutine makeup_multipass_slave (lat, slave)

implicit none

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord, patch_in_lord, patch_in_slave
type (branch_struct), pointer :: branch
type (floor_position_struct), pointer :: f0, f1
type (coord_struct) start, end

real(rp) s, slave_val(num_ele_attrib$), arg
real(rp) d, e, r_lord, r_slave, cos_lord, cos_e, sin_lord, sin_lorde
real(rp) ang_slave, ang_lord, ang_slave_old, d1, d2
real(rp) cos_e2, d_theta, ang_dlord, cos_lorde1, cos_dlord
real(rp) w0_mat(3,3), w1_mat(3,3), w1_inv_mat(3,3), offset(3), dw_mat(3,3)
real(rp) theta, phi, psi, w0_inv_mat(3,3)

integer i, j, ix_slave, ic, ix_s0, ix_patch_in_slave, n_pass
character(40) :: r_name = 'makeup_multipass_slave'

!

branch => lat%branch(slave%ix_branch)
call set_ele_status_stale (slave, attribute_group$)
slave%bookkeeping_state%control = ok$

ix_slave = slave%ix_ele
j =  lat%ic(slave%ic1_lord)
lord => lat%ele(lat%control(j)%ix_lord)
n_pass = j - lord%ix1_slave + 1  ! pass number for slave

slave_val = slave%value  ! save

slave%value = lord%value
if (lord%key == lcavity$ .or. lord%key == rfcavity$) then
  slave%value(dphi0$)        = slave_val(dphi0$)
endif

! A slave's field_master = T irregardless of the lord's setting.
! This is to make attribute_bookkeeper compute the correct normalized field strength.

slave%value(E_tot_start$)    = slave_val(E_tot_start$)
slave%value(p0c_start$)      = slave_val(p0c_start$)
slave%value(e_tot$)          = slave_val(e_tot$)
slave%value(p0c$)            = slave_val(p0c$)
slave%value(delta_ref_time$) = slave_val(delta_ref_time$)
slave%value(ref_time_start$) = slave_val(ref_time_start$)
slave%value(n_ref_pass$)     = 0
if (attribute_index(slave, 'FIELD_MASTER') /= 0) slave%field_master = .true.

! A match element with match_end$: Restore initial Twiss parameters (which
! are calculated in twiss_propagate1).

if (lord%key == match$) then
  if (lord%value(match_end$) /= 0) then
    slave%value(beta_a0$)    = slave_val(beta_a0$)
    slave%value(beta_b0$)    = slave_val(beta_b0$)
    slave%value(alpha_a0$)   = slave_val(alpha_a0$)
    slave%value(alpha_b0$)   = slave_val(alpha_b0$)
    slave%value(eta_x0$)     = slave_val(eta_x0$)
    slave%value(eta_y0$)     = slave_val(eta_y0$)
    slave%value(etap_x0$)    = slave_val(etap_x0$)
    slave%value(etap_y0$)    = slave_val(etap_y0$)
    slave%value(c_11$:c_22$) = slave_val(c_11$:c_22$)
    slave%value(gamma_c$)    = slave_val(gamma_c$)
  endif

  if (lord%value(match_end_orbit$) /= 0) then
    slave%value(x0$)  = slave_val(x0$)
    slave%value(px0$) = slave_val(px0$)
    slave%value(y0$)  = slave_val(y0$)
    slave%value(py0$) = slave_val(py0$)
    slave%value(z0$)  = slave_val(z0$)
    slave%value(pz0$) = slave_val(pz0$)
  endif
endif

! multipoles

if (associated (slave%a_pole)) then
  slave%a_pole           = lord%a_pole
  slave%b_pole           = lord%b_pole
  slave%multipoles_on    = lord%multipoles_on
  slave%scale_multipoles = lord%scale_multipoles
endif

! RF wakes

call transfer_rf_wake (lord%rf_wake, slave%rf_wake)

if (associated (slave%rf_wake)) then
  do i = 1, size(lord%rf_wake%lr)
    slave%rf_wake%lr(i)%t_ref = lord%rf_wake%lr(i)%t_ref - slave%ref_time
  enddo
endif

! EM fields
! Note: %em_field%mode%map, %em_field_mode%grid are in common between lord and slave.

if (associated(slave%em_field)) then
  do i = 1, size(lord%em_field%mode)
    slave%em_field%mode(i) = lord%em_field%mode(i)
  enddo
endif

! methods

slave%mat6_calc_method = lord%mat6_calc_method
slave%tracking_method  = lord%tracking_method
slave%map_with_offsets = lord%map_with_offsets
slave%is_on            = lord%is_on
slave%aperture_at      = lord%aperture_at
slave%aperture_type    = lord%aperture_type

! patch element.
! The reference energy may be zero while parsing in a lattice file so only do
! the computation if we have a non-zero energy.

if (lord%key == patch$ .and. slave%value(p0c$) /= 0) then

  select case (lord%ref_orbit)

  case (patch_in$)
    ic = lord%ix1_slave + nint(lord%value(n_ref_pass$)) - 1
    ix_s0 = lat%control(ic)%ix_slave  ! Index of slave element on the reference pass

    if (.not. has_orientation_attributes(slave)) call err_exit ! This should not be

    ! ref pass slave parameters are fixed.
    if (ix_s0 == ix_slave) then
      slave%value(x_offset$) = slave_val(x_offset$)
      slave%value(y_offset$) = slave_val(y_offset$)
      slave%value(z_offset$) = slave_val(z_offset$)
      slave%value(x_pitch$)  = slave_val(x_pitch$)
      slave%value(y_pitch$)  = slave_val(y_pitch$)
      slave%value(tilt$)     = slave_val(tilt$)
      return     
    endif

    f0 => branch%ele(ix_s0)%floor      ! Coords at ref slave exit end.
    f1 => branch%ele(ix_slave-1)%floor ! Coords at this slave entrance end.
    call floor_angles_to_w_mat (f0%theta, f0%phi, f0%psi, w0_mat)
    call floor_angles_to_w_mat (f1%theta, f1%phi, f1%psi, w1_mat)

    call mat_inverse (w1_mat, w1_inv_mat)
    dw_mat = matmul (w1_inv_mat, w0_mat) 
    call floor_w_mat_to_angles (dw_mat, 0.0_rp, theta, phi, psi)
    slave%value(x_pitch$) = theta
    slave%value(y_pitch$) = phi
    slave%value(tilt$) = psi

    offset = [f0%x-f1%x, f0%y-f1%y, f0%z-f1%z]
    if (slave%value(translate_after$) == 0) then
      offset = matmul(w1_inv_mat, offset)
    else
      call mat_inverse (w0_mat, w0_inv_mat)
      offset = matmul(w0_inv_mat, offset)
    endif
    slave%value(x_offset$) = offset(1)
    slave%value(y_offset$) = offset(2)
    slave%value(z_offset$) = offset(3)

  case (patch_out$)
    ic = lord%ix1_slave + nint(lord%value(n_ref_pass$)) - 1
    ix_s0 = lat%control(ic)%ix_slave  ! Index of slave element on the reference pass

    patch_in_lord => pointer_to_lord(lord, 1)
    patch_in_slave => pointer_to_slave(patch_in_lord, n_pass)
    start%vec = 0
    do i = patch_in_slave%ix_ele, ix_slave - 1
      call track1 (start, branch%ele(i), branch%param, end)
      start = end
    enddo
    slave%value(x_offset$) = end%vec(1) 
    slave%value(y_offset$) = end%vec(3)
    slave%value(z_offset$) = end%vec(5)
    slave%value(x_pitch$) = end%vec(2)
    slave%value(y_pitch$) = end%vec(4)

  end select
endif

! An sbend is tricky since the reference orbit changes with energy.

if (lord%key == sbend$ .and. slave%value(p0c$) /= 0 .and. lord%value(g$) /= 0) then

  if (lord%ref_orbit /= single_ref$) then
    if (lord%value(k1$) /= 0 .or. lord%value(k2$) /= 0 .or. associated(lord%a_pole)) then
      call out_io (s_fatal$, r_name, &
            'MULTIPASS BEND ELEMENT: ' // lord%name, &
            'WITH THE REF_ORBIT ATTRIBUTE SET TO: ' // ref_orbit_name(lord%ref_orbit), &
            'HAS A NONZERO HIGHER ORDER MULTIPOLE!', &
            'THIS IS NOT ALLOWED. SEE THE BMAD MANUAL FOR MORE DETAILS.')
      if (global_com%exit_on_error) call err_exit
    endif
  endif

  ! 

  select case (lord%ref_orbit)

  ! %ref_orbit = single_ref$

  case (single_ref$)
    slave%value(b_field$)     = lord%value(b_field$) * slave%value(p0c$) / lord%value(p0c$) 
    slave%value(b_field_err$) = lord%value(b_field$) + lord%value(b_field_err$) - &
                                                                   slave%value(b_field$)
    slave%value(g_err$) = (lord%value(g$) + lord%value(g_err$)) * &
                                    lord%value(p0c$) / slave%value(p0c$) - lord%value(g$)

  ! %ref_orbit = match_global_coords$

  case (match_global_coords$)

    slave%value(g$) = lord%value(g$) * lord%value(p0c$) / slave%value(p0c$)

    ! e1 and e2 and l are determined by the reference orbit of this pass with respect
    ! to the reference orbit of the reference pass. 
    ! Assumption: the slave element lies in the (x, z) plane.

    if (slave%floor%phi /= 0 .or. slave%floor%psi /= 0) then
       call out_io (s_fatal$, r_name, 'MULTIPASS ELEMENT: ' // lord%name, &
                     'WHICH HAS REF_ORBIT = MATCH_GLOBAL_COORDS DOES NOT LIE IN THE (X, Z) PLANE!')
      if (global_com%exit_on_error) call err_exit
    endif

    ic = lord%ix1_slave + nint(lord%value(n_ref_pass$)) - 1
    ix_s0 = lat%control(ic)%ix_slave  ! Index of slave element on the reference pass
    if (ix_s0 == ix_slave) return     ! Do not need calculation for ref slave.

    f0 => branch%ele(ix_s0-1)%floor    ! Coords at ref slave entrance end.
    f1 => branch%ele(ix_slave-1)%floor ! Coords at this slave entrance end.

    d_theta = modulo2(f1%theta - f0%theta, pi)
    !! if (abs(d_theta) > pi/4) return  ! Stop calc if too unphysical.

    ! d1 is the distance between the reference trajectory entrance points between 
    ! the slave and the lord.
    d1 = ((f1%x - f0%x) * cos(f1%theta) - (f1%y - f0%y) * cos(f1%theta)) / &
                                        cos(d_theta + lord%value(e1$))
    !! if (abs(d1 * slave%value(g$)) > 0.1) return  ! Stop calc if too unphysical.

    ! Iterate to converge to a solution.

    r_lord  = 1 / lord%value(g$)
    r_slave = 1 / slave%value(g$)
    ang_lord = lord%value(angle$)
    ang_dlord = d_theta + ang_lord
    cos_lord = cos(ang_lord);   cos_lorde1 = cos(ang_lord - lord%value(e1$))
    sin_lord = sin(ang_lord);   cos_dlord = cos(ang_dlord)
    cos_e2 = cos(lord%value(e2$))
    ang_slave     = ang_lord   ! Init guess
    ang_slave_old = ang_slave  
    do i = 1, 10  ! limit interations in case of nonconvergance
      d2 = (r_lord * (cos_lord - 1) + d1 * cos_lorde1 + r_slave * &
              (cos(ang_dlord - ang_slave) - cos_dlord)) / cos_e2
      arg = (r_lord * (sin(ang_dlord) - sin(d_theta)) - &
                        d1 * sin(d_theta + lord%value(e1$)) + &
                        d2 * sin(ang_dlord - lord%value(e2$))) / r_slave
      if (abs(arg) > 1) then
        call out_io (s_error$, r_name, &
              'MULTIPASS MATCH_GLOBAL_COORDS CALC ERROR FOR: ' // lord%name, &
              'MATCHING ABORTED. CURRENT PARAMETERS ARE NOT CORRECT.')
        if (global_com%exit_on_error) call err_exit
      endif
      ang_slave = asin(arg)
      if (abs(ang_slave - ang_slave_old) < 1e-6 * abs(ang_slave)) exit
      ang_slave_old = ang_slave
    enddo

    slave%value(angle$) = ang_slave
    slave%value(l$) = ang_slave * r_slave
    slave%value(rho$) = r_slave
    slave%value(e1$) = lord%value(e1$) + d_theta
    slave%value(e2$) = lord%value(e2$) + ang_slave - ang_lord - d_theta

  ! %ref_orbit = match_at_entrance$ or match_at_exit$

  case (match_at_entrance$, match_at_exit$)
    slave%value(g$) = lord%value(g$) * lord%value(p0c$) / slave%value(p0c$)
    ! Iterate to converge to a solution
    if (slave%value(g$) /= lord%value(g$)) then
      if (lord%ref_orbit == match_at_entrance$) then
        e = lord%value(e2$)
      else
        e = lord%value(e1$)
      endif
      r_lord  = 1 / lord%value(g$)
      r_slave = 1 / slave%value(g$)
      ang_lord = lord%value(angle$)
      cos_lord = cos(ang_lord); cos_e = cos(e)
      sin_lord = sin(ang_lord); sin_lorde = sin(ang_lord - e)
      ang_slave     = ang_lord
      ang_slave_old = ang_slave
      ! d is the distance between the reference trajectory end points between the slave
      ! and the lord at the opposite end of the match end.
      do
        d = (r_lord * (cos_lord - 1) + r_slave * (cos(ang_lord - ang_slave) - cos_lord) ) / cos_e
        arg = (r_lord * sin_lord + d * sin_lorde) / r_slave
        if (abs(arg) > 1) then
          call out_io (s_error$, r_name, &
                'MULTIPASS MATCH_AT_ENTRANCE/EXIT CALC ERROR FOR: ' // lord%name, &
                'MATCHING ABORTED. CURRENT PARAMETERS ARE NOT CORRECT.')
          if (global_com%exit_on_error) call err_exit
        endif
        ang_slave = asin(arg)
        if (abs(ang_slave - ang_slave_old) < 1e-6 * abs(ang_slave)) exit
        ang_slave_old = ang_slave
      enddo
      slave%value(angle$) = ang_slave
      slave%value(l$) = ang_slave * r_slave
      slave%value(rho$) = r_slave
      if (lord%ref_orbit == match_at_entrance$) then
        slave%value(e2$) = e + ang_slave - ang_lord 
      else
        slave%value(e1$) = e + ang_slave - ang_lord 
      endif
    endif

  case default
    call out_io (s_fatal$, r_name, 'BAD REF_ORBIT VALUE: \i0\ ', &
                           'FOR: ' // lord%name, i_array = [lord%ref_orbit] )
    if (global_com%exit_on_error) call err_exit

  end select
endif

end subroutine makeup_multipass_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (lat, slave)
!
! Subroutine to calcualte the attributes of superposition slave elements.
! This routine is not meant for general use.
!-
       
subroutine makeup_super_slave (lat, slave)

implicit none

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord, slave0, lord1
type (ele_struct) :: sol_quad
type (branch_struct), pointer :: branch

integer i, j, ix_con, ix, ix_slave, ix_lord, ix_order

real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3, c, s
real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(num_ele_attrib$)
real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4), leng, offset
real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

logical is_first, is_last, err_flag, major_method_set_done 

character(20) :: r_name = 'makeup_super_slave'

! Super_slave:

branch => lat%branch(slave%ix_branch)
ix_slave = slave%ix_ele

slave%bookkeeping_state%control = ok$
call set_ele_status_stale (slave, attribute_group$)

if (slave%slave_status /= super_slave$) then
  call out_io(s_abort$, r_name, "ELEMENT IS NOT A SUPER SLAVE: " // slave%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

!

slave%field_calc = refer_to_lords$

!-----------------------------------------------------------------------
! 1 super_lord for this super_slave: just transfer attributes except length

if (slave%n_lord == 1) then

  lord => pointer_to_lord(slave, 1, ix_con, ix_order)
  is_first = (ix_order == 1)
  is_last  = (ix_order == lord%n_slave)

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (.not. is_first) then
    if (.not. all(slave%map_ref_orb_in == branch%ele(ix_slave-1)%map_ref_orb_out)) then
      slave0 => pointer_to_slave(lord, ix_order-1)
      slave%map_ref_orb_in = slave0%map_ref_orb_out
      if (associated(slave%rad_int_cache)) slave%rad_int_cache%stale = .true. ! Forces recalc
    endif
  endif

  ! Find the offset from the longitudinal start of the lord to the start of the slave

  offset = 0 ! length of all slaves before this one
  do i = 1, ix_con - lord%ix1_slave
    slave0 => pointer_to_slave(lord, i)
    offset = offset + slave0%value(l$)
  enddo

  ! If this is the last slave, adjust it's length to be consistant with
  ! The lord length. Then do the rest of the bookkeeping

  if (is_last) slave%value(l$) = lord%value(l$) - offset

  call makeup_super_slave1 (slave, lord, offset, lat%branch(slave%ix_branch)%param, &
                                  is_first, is_last, err_flag)
  if (err_flag) return

  return

endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.
                                         
k_x = 0
k_y = 0
x_kick = 0
y_kick = 0
sum_1 = 0
sum_2 = 0
sum_3 = 0
sum_4 = 0
ks_sum = 0
ks_xp_sum = 0
ks_xo_sum = 0
ks_yp_sum = 0
ks_yo_sum = 0

!

value = 0
value(l$) = slave%value(l$)
value(E_tot_start$)    = slave%value(E_tot_start$)
value(p0c_start$)      = slave%value(p0c_start$)
value(E_tot$)          = slave%value(E_tot$)
value(p0c$)            = slave%value(p0c$)
value(delta_ref_time$) = slave%value(delta_ref_time$)
value(ref_time_start$) = slave%value(ref_time_start$)

s_slave = slave%s - value(l$)/2  ! center of slave
slave%is_on = .false.

major_method_set_done = .false.

! sum over all lords...

do j = 1, slave%n_lord

  lord => pointer_to_lord(slave, j, ix_con)
  is_first = (ix_con == lord%ix1_slave)
  is_last  = (ix_con == lord%ix2_slave)

  ! Do some error checking.

  if (lord%lord_status /= super_lord$) then
    call out_io (s_abort$, r_name, &
          "SUPER_SLAVE HAS A CONTROL ELEMENT THAT IS NOT A SUPER_LORD", &
          'SLAVE: ' //  slave%name // '  \i\ ', &
          'LORD:  ' //  lord%name  // '  \i\ ', i_array = [ix_slave, lord%ix_ele] )
    if (global_com%exit_on_error) call err_exit
  endif

  if (associated(lord%rf_wake)) then
    call out_io (s_abort$, r_name, &
            'SUPERPOSITION OF ELEMENTS WITH WAKES NOT YET IMPLEMENTED!', &
            'SUPER_LORD: ' // lord%name)
    if (global_com%exit_on_error) call err_exit
  endif

  ! Physically, the lord length cannot be less than the slave length.
  ! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

  if (abs(slave%value(l$)) > abs(lord%value(l$))) then
    coef = 1
  else
    coef = slave%value(l$) / lord%value(l$) 
  endif

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (.not. is_first) then
    if (.not. all(slave%map_ref_orb_in == branch%ele(ix_slave-1)%map_ref_orb_out)) then
      slave%map_ref_orb_in = branch%ele(ix_slave-1)%map_ref_orb_out
      if (associated(slave%rad_int_cache)) slave%rad_int_cache%stale = .true. ! Forces recalc
    endif
  endif

  ! Choose the smallest ds_step of all the lords.

  if (value(ds_step$) == 0 .or. lord%value(ds_step$) < value(ds_step$)) &
                                        value(ds_step$) = lord%value(ds_step$)

  ! Coupler and aperture calc.

  call compute_slave_aperture (value, slave, lord, is_first, is_last)

  if (slave%key == lcavity$) call compute_slave_coupler (value, slave, lord, is_first, is_last)

  ! Methods.
  ! Major_method_set_done = T means the slave methods have been set using a "major" element.
  ! A "major" element is something other than a pipe, monitor, etc.

  if (j == 1 .or. .not. major_method_set_done) then
    slave%mat6_calc_method = lord%mat6_calc_method
    slave%tracking_method  = lord%tracking_method
    slave%map_with_offsets = lord%map_with_offsets
  endif

  select case (lord%key)
  case (hkicker$, vkicker$, kicker$, instrument$, monitor$, pipe$, rcollimator$, ecollimator$)
  case default
    if (.not. major_method_set_done) then
      if (slave%mat6_calc_method /= lord%mat6_calc_method) then
        lord1 => pointer_to_lord(slave, 1)
        call out_io(s_abort$, r_name, 'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(lord1%name))
        if (global_com%exit_on_error) call err_exit
      endif
      if (slave%tracking_method /= lord%tracking_method) then
        lord1 => pointer_to_lord(slave, 1)
        call out_io(s_abort$, r_name, ' TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(lord1%name))
        if (global_com%exit_on_error) call err_exit
      endif
      if (slave%map_with_offsets .neqv. lord%map_with_offsets) then
        lord1 => pointer_to_lord(slave, 1)
        call out_io(s_abort$, r_name, 'MAP_WITH_OFFSETS DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(lord1%name))
        if (global_com%exit_on_error) call err_exit
      endif
    endif

    major_method_set_done = .true.
  end select

  ! descriptive strings.

  if (associated(lord%descrip)) then
    if (.not. associated(slave%descrip)) allocate (slave%descrip)
    slave%descrip = lord%descrip
  endif

  if (lord%type /= '') slave%type = lord%type
  if (lord%alias /= '') slave%alias = lord%alias

  !----------------------------------------------------
  ! kicks, etc.

  if (.not. lord%is_on) cycle
  slave%is_on = .true.  ! Slave is on if at least one lord is on

  if (slave%key == em_field$) cycle  ! Field info is stored in the lord elements.

  tilt = lord%value(tilt_tot$)

  if (lord%key == hkicker$) then
    x_kick = x_kick + lord%value(kick$) * cos(tilt) * coef
    y_kick = y_kick + lord%value(kick$) * sin(tilt) * coef
  elseif (lord%key == vkicker$) then
    x_kick = x_kick - lord%value(kick$) * sin(tilt) * coef
    y_kick = y_kick + lord%value(kick$) * cos(tilt) * coef
  elseif (lord%key == kicker$) then
    c = cos(tilt) * coef
    s = sin(tilt) * coef
    x_kick = x_kick + c * lord%value(hkick$) - s * lord%value(vkick$)
    y_kick = y_kick + s * lord%value(hkick$) + c * lord%value(vkick$)
  else
    x_kick = x_kick + lord%value(hkick$) * coef
    y_kick = y_kick + lord%value(vkick$) * coef
  endif

  !------

  select case (slave%key)

  ! sextupole

  case (sextupole$) 

    cos_n = lord%value(k2$) * cos(3 * tilt)
    sin_n = lord%value(k2$) * sin(3 * tilt)            
    
    k_x = k_x + cos_n
    k_y = k_y + sin_n

  ! octupole

  case (octupole$)

    cos_n = lord%value(k3$) * cos(4 * tilt)
    sin_n = lord%value(k3$) * sin(4 * tilt)        
    
    k_x = k_x + cos_n
    k_y = k_y + sin_n

  ! solenoid/quadrupole combo.

  case (solenoid$, sol_quad$, quadrupole$)

    x_p = lord%value(x_pitch_tot$);  x_o = lord%value(x_offset_tot$)
    y_p = lord%value(y_pitch_tot$);  y_o = lord%value(y_offset_tot$)

    s_del = s_slave - (lord%s + lord%value(s_offset_tot$) - lord%value(l$)/2)
    s_del = modulo2 (s_del, branch%param%total_length/2)

    ks = lord%value(ks$)

    ks_sum = ks_sum + ks

    ks_xp_sum = ks_xp_sum + ks * x_p
    ks_yp_sum = ks_yp_sum + ks * y_p

    ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
    ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

    cos_n = lord%value(k1$) * cos(2 * tilt)
    sin_n = lord%value(k1$) * sin(2 * tilt)

    k_x = k_x + cos_n
    k_y = k_y + sin_n

    sum_1 = sum_1 + cos_n * x_p + sin_n * y_p
    sum_2 = sum_2 + sin_n * x_p - cos_n * y_p

    sum_3 = sum_3 + cos_n * (x_o + x_p * s_del) + sin_n * (y_o + y_p * s_del)
    sum_4 = sum_4 + sin_n * (x_o + x_p * s_del) - cos_n * (y_o + y_p * s_del)

  ! bend_sol_quad

  case (bend_sol_quad$)
    call out_io (s_abort$, r_name, 'CODING NOT YET IMPLEMENTED FOR: ' // key_name(slave%key))
    if (global_com%exit_on_error) call err_exit

  ! Everything else

  case default
    ! Everything else has no special needs. 

  end select

enddo

if (slave%tracking_method == bmad_standard$ .and. slave%key == em_field$) slave%tracking_method = runge_kutta$
if (slave%mat6_calc_method == bmad_standard$ .and. slave%key == em_field$) slave%mat6_calc_method = tracking$

!-------------------------------------------------------------------------------
! stuff sums into slave element

if (slave%key == em_field$) then
  slave%value = value
  return  ! Field info is stored in the lord elements.
endif

if (x_kick == 0 .and. y_kick == 0) then
  if (slave%key == hkicker$ .or. slave%key == vkicker$) then
    value(kick$) = 0
  else
    value(hkick$) = 0
    value(vkick$) = 0
  endif
elseif (slave%key == hkicker$) then
  value(kick$) = sqrt(x_kick**2 + y_kick**2)
  value(tilt$) = atan2(y_kick, x_kick)
elseif (slave%key == vkicker$) then
  value(kick$) = sqrt(x_kick**2 + y_kick**2)
  value(tilt$) = atan2(-x_kick, y_kick)
elseif (slave%key == kicker$) then
  value(tilt$) = 0
  value(hkick$) = x_kick
  value(vkick$) = y_kick
else
  value(hkick$) = x_kick
  value(vkick$) = y_kick
endif

slave%value = value

!-----------------------------

select case (slave%key)

case (sextupole$) 

  if (k_x == 0 .and. k_y == 0) goto 8000

  k2 = sqrt(k_x**2 + k_y**2)
  tilt = atan2(k_y, k_x) / 3

  if (tilt > pi/6) then
    k2 = -k2
    tilt = tilt - pi/3
  elseif (tilt < -pi/6) then
    k2 = -k2
    tilt = tilt + pi/3
  endif

  slave%value(k2$) = k2
  slave%value(tilt$) = tilt

! octupole

case (octupole$)

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) goto 8000

  k3 = sqrt(k_x**2 + k_y**2)
  tilt = atan2(k_y, k_x) / 4

  if (tilt > pi/8) then
    k3 = -k3
    tilt = tilt - pi/4
  elseif (tilt < -pi/8) then
    k3 = -k3
    tilt = tilt + pi/4
  endif

  slave%value(k3$) = k3
  slave%value(tilt$) = tilt

! sol_quad, etc.

case (solenoid$, sol_quad$, quadrupole$)

  ks = ks_sum
  slave%value(ks$) = ks

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) goto 8000

  if (ks /= 0) then
    x_o_sol = ks_xo_sum / ks
    x_p_sol = ks_xp_sum / ks
    y_o_sol = ks_yo_sum / ks
    y_p_sol = ks_yp_sum / ks
  endif

  if (k_x == 0 .and. k_y == 0) then  ! pure solenoid
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    slave%value(x_offset$) = x_o_sol
    slave%value(y_offset$) = y_o_sol
    slave%value(x_pitch$)  = x_p_sol
    slave%value(y_pitch$)  = y_p_sol
  endif   

  ! here if have quadrupole component

  if (k_x /= 0 .or. k_y /= 0) then
    k1 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 2

    if (tilt > pi/4) then
      k1 = -k1
      tilt = tilt - pi/2
    elseif (tilt < -pi/4) then
      k1 = -k1
      tilt = tilt + pi/2
    endif

    slave%value(k1$) = k1
    slave%value(tilt$) = tilt

    cos_n = k_x / (k_x**2 + k_y**2)
    sin_n = k_y / (k_x**2 + k_y**2)

    slave%value(x_pitch$)  = cos_n * sum_1 + sin_n * sum_2
    slave%value(y_pitch$)  = sin_n * sum_1 - cos_n * sum_2
    slave%value(x_offset$) = cos_n * sum_3 + sin_n * sum_4
    slave%value(y_offset$) = sin_n * sum_3 - cos_n * sum_4
  endif

  ! if ks /= 0 then we have to recalculate the offsets and pitches.

  if (ks /= 0 .and. (k_x /= 0 .or. k_y /= 0)) then

    x_p = slave%value(x_pitch$) - x_p_sol; x_o = slave%value(x_offset$) - x_o_sol
    y_p = slave%value(y_pitch$) - y_p_sol; y_o = slave%value(y_offset$) - y_o_sol

    if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) goto 8000

    t_2 = [x_o, x_p, y_o, y_p]
    call tilt_coords (tilt, t_2)

    l_slave = slave%value(l$)

    t_1 = [t_2(2), 0.0_rp, t_2(4), 0.0_rp]
    t_2(1) = t_2(1) + ks * t_2(4) / k1 
    t_2(3) = t_2(3) + ks * t_2(2) / k1
             
    call mat_make_unit (T_end)
    T_end(4,1) =  ks / 2
    T_end(2,3) = -ks / 2

    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    sol_quad%value(ks$) = ks
    sol_quad%value(k1$) = k1
    sol_quad%value(l$)  = l_slave
    call make_mat6 (sol_quad, branch%param)
    T_tot = sol_quad%mat6(1:4,1:4)

    r_off = matmul (T_end, l_slave * t_1 / 2 - t_2) 
    r_off = matmul (T_tot, r_off) + matmul (T_end, l_slave * t_1 / 2 + t_2)

    call mat_make_unit (mat4)
    mat4(:,2) = mat4(:,2) + l_slave * T_tot(:,1) / 2
    mat4(:,4) = mat4(:,4) + l_slave * T_tot(:,3) / 2
    mat4(1,2) = mat4(1,2) + l_slave / 2
    mat4(3,4) = mat4(3,4) + l_slave / 2
    mat4 = mat4 - T_tot

    call mat_inverse (mat4, mat4_inv)
    beta = matmul (mat4_inv, r_off)

    call tilt_coords (-tilt, beta)

    slave%value(x_offset$) = beta(1) + x_o_sol
    slave%value(x_pitch$)  = beta(2) + x_p_sol
    slave%value(y_offset$) = beta(3) + y_o_sol
    slave%value(y_pitch$)  = beta(4) + y_p_sol
  endif

! bend_sol_quad

case (bend_sol_quad$)
  call out_io (s_abort$, r_name, &
               'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
  if (global_com%exit_on_error) call err_exit

end select

! If the slave has %field_master = T then we need to convert k1, etc values to field quantities.

8000 continue

if (slave%field_master) then
  slave%field_master = .false.   ! So attribute_bookkeeper will do the right thing.
  call attribute_bookkeeper (slave, branch%param)
  slave%field_master = .true.
else
  call attribute_bookkeeper (slave, branch%param)
endif

end subroutine makeup_super_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset,
!                                    param, at_entrance_end, at_exit_end, err_flag, old_slice)
!
! Routine to create an element that represents a longitudinal slice of the original element.
! Note: This routine assumes that the following call has been made before hand:
!    call transfer_ele (ele_in, sliced_ele, .true.)
! This routine only has to be done once for the life of the sliced_ele variable.
!
! Note: To save tracking computation time, if ele_in has taylor, symp_lie_ptc, or symp_map 
! for tracking_method or mat6_calc_method, then this will be changed to symp_lie_bmad 
! for wigglers and bmad_standard for everything else.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in          -- Ele_struct: Original element to slice
!   l_slice         -- Real(rp): Length of the slice
!   offset          -- Real(rp): Offset of entrance end of sliced_ele from entrance end of ele_in.
!   param           -- Lat_param_struct: lattice paramters.
!   at_entrance_end -- Logical: Sliced_ele contains the ele's entrance end?
!   at_exit_end     -- Logical: Sliced_ele contains the ele's exit end?
!   old_slice       -- Logical, optional: Previous slice. If present this saves computation
!                        time of the refernece energy and time at the start of the present slice.
!
! Output:
!   sliced_ele -- Ele_struct: Sliced_ele element with appropriate values set.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!-

recursive subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset, &
                                       param, at_entrance_end, at_exit_end, err_flag, old_slice)

implicit none

type (ele_struct), target :: sliced_ele, ele_in
type (ele_struct), optional :: old_slice
type (ele_struct) :: ele2
type (lat_param_struct) param

real(rp) l_slice, offset, e_len, ref_time_start, p0c_start, e_tot_start
real(rp) time_ref_orb_out(6)

logical at_entrance_end, at_exit_end, err_flag, err2_flag

character(24) :: r_name = 'create_element_slice'

! Err check. Remember: the element length may be negative

err_flag = .true.
sliced_ele%ix_ele = -2  ! Indicate sliced ele is not an element in the lattice.

if (ele_in%key == taylor$ .or. ele_in%key == hybrid$) then
  call out_io (s_fatal$, r_name, &
        'CANNOT SLICE ELEMENT OF TYPE: ' // key_name(ele_in%key), &
        'CANNOT SLICE: ' // ele_in%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

e_len = ele_in%value(l$)
if (l_slice*e_len < 0 .or. abs(l_slice) > abs(e_len) + bmad_com%significant_length) then
  call out_io (s_fatal$, r_name, &
        'SLICE LENGTH IS OUT OF RANGE FOR ELEMENT: ' // ele_in%name, &
        'LENGTH: \2es12.3\ ', r_array = [l_slice, e_len])
  if (global_com%exit_on_error) call err_exit
  return
endif

! Simple case where ele length is zero

if (e_len == 0) then
  sliced_ele = ele_in
  err_flag = .false.
  return
endif

! Save values from old_slice if present

if (present(old_slice)) then
  p0c_start       = old_slice%value(p0c$)
  e_tot_start     = old_slice%value(e_tot$)
  ref_time_start  = old_slice%ref_time
  time_ref_orb_out = old_slice%time_ref_orb_out
endif

! The sliced element is treated as a super_slave to the original element except
! if that element is a super_slave in which case the sliced element has the same lords
! as the original element.

sliced_ele%value(l$) = l_slice

if (ele_in%slave_status /= super_slave$) then
  sliced_ele%slave_status = slice_slave$
  sliced_ele%n_lord = 1
  sliced_ele%lord   => ele_in
endif

call makeup_super_slave1 (sliced_ele, ele_in, offset, param, at_entrance_end, at_exit_end, err2_flag)
if (err2_flag) return

sliced_ele%s = ele_in%s - e_len + offset + sliced_ele%value(l$)

! Use a speedier tracking method.

select case (sliced_ele%tracking_method)
case (taylor$, symp_map$, symp_lie_ptc$)
  select case (sliced_ele%key)
  case (wiggler$); sliced_ele%tracking_method = symp_lie_bmad$
  case default;    sliced_ele%tracking_method = bmad_standard$
  end select
end select

select case (sliced_ele%mat6_calc_method)
case (taylor$, symp_map$, symp_lie_ptc$)
  select case (sliced_ele%key)
  case (wiggler$); sliced_ele%mat6_calc_method = symp_lie_bmad$
  case default;    sliced_ele%mat6_calc_method = bmad_standard$
  end select
end select

sliced_ele%field_calc = refer_to_lords$

! Makeup_super_slave1 does not compute reference energy or time so need to do it here.

if (offset == 0) then
  p0c_start      = ele_in%value(p0c_start$)
  e_tot_start    = ele_in%value(e_tot_start$)
  ref_time_start = ele_in%value(ref_time_start$)
  sliced_ele%time_ref_orb_in = ele_in%time_ref_orb_in
elseif (present(old_slice)) then
  sliced_ele%time_ref_orb_in = time_ref_orb_out
elseif (ele_has_constant_ds_dt_ref(ele_in)) then
  p0c_start      = ele_in%value(p0c$)
  e_tot_start    = ele_in%value(e_tot$)
  ref_time_start = ele_in%ref_time - ele_in%value(delta_ref_time$) * (ele_in%value(l$) - offset) / ele_in%value(l$)
  sliced_ele%time_ref_orb_in = 0
else
  call transfer_ele (sliced_ele, ele2)
  call create_element_slice (ele2, ele_in, offset, 0.0_rp, param, .true., .false., err2_flag)
  if (err2_flag) return
  p0c_start      = ele2%value(p0c$)
  e_tot_start    = ele2%value(e_tot$)
  ref_time_start = ele2%ref_time
  sliced_ele%time_ref_orb_in = ele2%time_ref_orb_out
endif

call ele_compute_ref_energy_and_time (sliced_ele, param, e_tot_start, p0c_start, ref_time_start, err2_flag)
if (err2_flag) return

err_flag = .false.

end subroutine create_element_slice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave1 (slave, lord, offset, param, at_entrance_end, at_exit_end)
!
! Routine to construct a super_slave from a super_lord when the slave has only one lord.
! Note: Reference energy and times are not computed in this routine.
!
! Modules needed:
!   use bmad
!
! Input:
!   slave  -- Ele_struct: Slave element.
!   lord   -- Ele_struct: Lord element.
!   offset -- Real(rp): offset of entrance end of slave from entrance end of the lord.
!   param  -- Lat_param_struct: lattice paramters.
!   at_entrance_end -- Logical: Slave contains the lord's entrance end?
!   at_exit_end     -- Logical: Slave contains the lord's exit end?
!
! Output:
!   slave    -- Ele_struct: Slave element with appropriate values set.
!   err_flag -- Logical: Set true if there is an error. False otherwise.
!-

subroutine makeup_super_slave1 (slave, lord, offset, param, at_entrance_end, at_exit_end, err_flag)

implicit none

type (ele_struct), target :: slave, lord
type (lat_param_struct) param

real(rp) offset, s_del, coef
real(rp) value(num_ele_attrib$)
integer i
logical at_entrance_end, at_exit_end, err_flag
character(24) :: r_name = 'makeup_super_slave1'

! Physically, the lord length cannot be less than the slave length.
! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

err_flag = .true.

if (lord%value(l$) == 0) then
  call out_io (s_fatal$, r_name, 'LORD HAS ZERO LENGTH!')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (abs(slave%value(l$)) > abs(lord%value(l$))) then
  coef = 1
else
  coef = slave%value(l$) / lord%value(l$) 
endif

! Reference energy and time computed in ele_compute_ref_energy_and_time.

value = lord%value
value(l$)              = slave%value(l$)                ! do not change slave length, etc.
value(delta_ref_time$) = slave%value(delta_ref_time$)
value(ref_time_start$) = slave%value(ref_time_start$)
value(E_tot_start$)    = slave%value(E_tot_start$)
value(p0c_start$)      = slave%value(p0c_start$)
value(E_tot$)          = slave%value(E_tot$)
value(p0c$)            = slave%value(p0c$)
value(num_steps$)      = slave%value(num_steps$)

!

if (lord%key == hkicker$ .or. lord%key == vkicker$) then
  value(kick$)    = lord%value(kick$) * coef
  value(bl_kick$) = lord%value(bl_kick$) * coef
else
  value(hkick$)    = lord%value(hkick$) * coef
  value(vkick$)    = lord%value(vkick$) * coef
  value(bl_hkick$) = lord%value(bl_hkick$) * coef
  value(bl_vkick$) = lord%value(bl_vkick$) * coef
endif

if (slave%key == rfcavity$) value(voltage$) = lord%value(voltage$) * coef

call compute_slave_aperture (value, slave, lord, at_entrance_end, at_exit_end)

if (slave%key == lcavity$) then
  slave%value(coupler_at$) = no_end$
  call compute_slave_coupler (value, slave, lord, at_entrance_end, at_exit_end)
endif

! s_del is the distance between lord and slave centers

if (has_orientation_attributes(slave)) then
  s_del = offset + slave%value(l$)/2 - lord%value(l$)/2
  value(x_pitch$)  = value(x_pitch_tot$)
  value(y_pitch$)  = value(y_pitch_tot$)
  value(x_offset$) = value(x_offset_tot$) + s_del * value(x_pitch_tot$)
  value(y_offset$) = value(y_offset_tot$) + s_del * value(y_pitch_tot$)
  value(tilt$)     = value(tilt_tot$)
endif

slave%value = value
slave%is_on = lord%is_on
slave%mat6_calc_method = lord%mat6_calc_method
slave%tracking_method  = lord%tracking_method
slave%map_with_offsets = lord%map_with_offsets

if (slave%tracking_method == bmad_standard$ .and. slave%key == em_field$) slave%tracking_method = runge_kutta$
if (slave%mat6_calc_method == bmad_standard$ .and. slave%key == em_field$) slave%mat6_calc_method = tracking$

! wiggler fields and electro-magnetic fields

if (slave%key == wiggler$) slave%value(n_pole$) = lord%value(n_pole$) * coef

! If an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

if (slave%key == sbend$) then
  if (.not. at_entrance_end) then 
    slave%value(e1$)    = 0
    slave%value(h1$)    = 0
    slave%value(fint$)  = 0
    slave%value(hgap$)  = 0
  endif
  if (.not. at_exit_end) then   ! first slave bend
    slave%value(e2$)    = 0
    slave%value(h2$)    = 0
    slave%value(fintx$) = 0
    slave%value(hgapx$) = 0
  endif
endif                       

! If there are long range wakes they must be scaled.

if (associated (slave%rf_wake)) then
  slave%rf_wake%lr%freq_in   = lord%rf_wake%lr%freq_in
  slave%rf_wake%lr%freq      = lord%rf_wake%lr%freq
  slave%rf_wake%lr%Q         = lord%rf_wake%lr%Q
  slave%rf_wake%lr%angle     = lord%rf_wake%lr%angle
  slave%rf_wake%lr%m         = lord%rf_wake%lr%m
  slave%rf_wake%lr%polarized = lord%rf_wake%lr%polarized
  slave%rf_wake%lr%r_over_q  = lord%rf_wake%lr%r_over_q * coef
endif

!

if (slave%key == lcavity$) then
  slave%value(e_loss$) = lord%value(e_loss$) * coef
endif

slave%bookkeeping_state%attributes = stale$
call attribute_bookkeeper (slave, param)

err_flag = .false.

end subroutine makeup_super_slave1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_aperture (value, slave, lord, at_entrance_end, at_exit_end)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_aperture (value, slave, lord, at_entrance_end, at_exit_end)

implicit none

type (ele_struct) slave, lord
real(rp) value(num_ele_attrib$)
logical at_entrance_end, at_exit_end

! 

slave%aperture_at = no_end$

select case (lord%aperture_at)
case (exit_end$) 
  if (at_exit_end) slave%aperture_at = exit_end$
case (entrance_end$)
  if (at_entrance_end) slave%aperture_at = entrance_end$
case (both_ends$)
  if (at_entrance_end .and. at_exit_end) then
    slave%aperture_at = both_ends$
  elseif (at_entrance_end) then
    slave%aperture_at = entrance_end$
  elseif (at_exit_end) then 
    slave%aperture_at = exit_end$
  endif
case (continuous$)
  slave%aperture_at = continuous$
end select

if (slave%aperture_at == no_end$) then
  value(x1_limit$) = 0
  value(x2_limit$) = 0
  value(y1_limit$) = 0
  value(y2_limit$) = 0
endif

end subroutine compute_slave_aperture

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_coupler (value, slave, lord, at_entrance_end, at_exit_end)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_coupler (value, slave, lord, at_entrance_end, at_exit_end)

implicit none

type (ele_struct) slave, lord
real(rp) value(num_ele_attrib$)
logical at_entrance_end, at_exit_end

!

select case (nint(lord%value(coupler_at$)))
case (exit_end$) 
  if (at_exit_end) slave%value(coupler_at$) = exit_end$
case (entrance_end$)
  if (at_entrance_end) slave%value(coupler_at$) = entrance_end$
case (both_ends$)
  if (at_entrance_end .and. at_exit_end) then
    slave%value(coupler_at$) = both_ends$
  elseif (at_entrance_end) then
    slave%value(coupler_at$) = entrance_end$
  elseif (at_exit_end) then 
    slave%value(coupler_at$) = exit_end$
  endif
end select

if (nint(slave%value(coupler_at$)) == no_end$) then
  value(coupler_strength$) = 0
endif

end subroutine compute_slave_coupler

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_overlay_and_girder_slave (lat, slave)
!
! This routine is not meant for general use.
!-

subroutine makeup_overlay_and_girder_slave (lat, slave)

implicit none

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (branch_struct), pointer :: branch

real(rp) value(num_ele_attrib_extended$), coef, ds, s_slave
real(rp) t, x_off, y_off, x_pitch, y_pitch
real(rp), pointer :: r_ptr
integer i, ix_con, ix, iv, ix_slave, icom, l_stat
logical used(num_ele_attrib_extended$), multipole_set, err_flag

character(40) :: r_name = 'makeup_overlay_and_girder_slave'

!
                             
branch => lat%branch(slave%ix_branch)

slave%bookkeeping_state%control = ok$
call set_ele_status_stale (slave, attribute_group$)

l_stat = slave%lord_status
ix_slave = slave%ix_ele

value = 0
used = .false.
multipole_set = .false.
slave%on_a_girder = .false.

do i = 1, slave%n_lord
  lord => pointer_to_lord(slave, i, ix_con)

  if (lord%lord_status == multipass_lord$) cycle
  if (lord%lord_status == group_lord$) cycle

  if (lord%lord_status == girder_lord$ .and. has_orientation_attributes(slave)) then
    s_slave = slave%s - slave%value(l$)/2
    if (s_slave > lord%value(s_max$)) s_slave = s_slave - lat%branch(slave%ix_branch)%param%total_length
    ds = s_slave - lord%value(s_center$) 
    if (lord%value(tilt$) == 0) then
      x_off = slave%value(x_offset$)
      y_off = slave%value(y_offset$)
      x_pitch = slave%value(x_pitch$)
      y_pitch = slave%value(y_pitch$)
    else
      t = lord%value(tilt$)
      x_off = slave%value(x_offset$) * cos(t) - slave%value(y_offset$) * sin(t)
      y_off = slave%value(x_offset$) * sin(t) + slave%value(y_offset$) * cos(t)
      x_pitch = slave%value(x_pitch$) * cos(t) - slave%value(y_pitch$) * sin(t)
      y_pitch = slave%value(x_pitch$) * sin(t) + slave%value(y_pitch$) * cos(t)
    endif
    slave%value(x_offset_tot$) = x_off + ds * lord%value(x_pitch$) + lord%value(x_offset$)
    slave%value(y_offset_tot$) = y_off + ds * lord%value(y_pitch$) + lord%value(y_offset$)
    slave%value(s_offset_tot$) = slave%value(s_offset$) + lord%value(s_offset$)
    slave%value(x_pitch_tot$)  = x_pitch  + lord%value(x_pitch$)
    slave%value(y_pitch_tot$)  = y_pitch  + lord%value(y_pitch$)
    slave%value(tilt_tot$)     = slave%value(tilt$)     + lord%value(tilt$)
    slave%on_a_girder = .true.
    cycle
  endif

  if (lord%lord_status /= overlay_lord$) then
    call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY_LORD \i\ ', ix_slave)
    call type_ele (slave, .true., 0, .false., 0, .true., lat)
    if (global_com%exit_on_error) call err_exit
  endif     

  coef = lat%control(ix_con)%coef
  iv = lat%control(ix_con)%ix_attrib
  call pointer_to_indexed_attribute (lord, lord%ix_value, .false., r_ptr, err_flag)
  if (err_flag) call err_exit
  value(iv) = value(iv) + r_ptr * coef
  used(iv) = .true.
  if (iv > num_ele_attrib$) multipole_set = .true.
enddo

where (used(1:num_ele_attrib$)) slave%value = value(1:num_ele_attrib$)
if (multipole_set) then
  where (used(a0$:a20$)) slave%a_pole = value(a0$:a20$)
  where (used(b0$:b20$)) slave%b_pole = value(b0$:b20$)
endif

! If no girder then simply transfer tilt to tilt_tot, etc.

if (.not. slave%on_a_girder .and. has_orientation_attributes(slave)) then
  slave%value(tilt_tot$)     = slave%value(tilt$)
  slave%value(x_offset_tot$) = slave%value(x_offset$)
  slave%value(y_offset_tot$) = slave%value(y_offset$)
  slave%value(s_offset_tot$) = slave%value(s_offset$)
  slave%value(x_pitch_tot$)  = slave%value(x_pitch$)
  slave%value(y_pitch_tot$)  = slave%value(y_pitch$)
endif

end subroutine makeup_overlay_and_girder_slave 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine attribute_bookkeeper (ele, param, force_bookkeeping)
!
! Routine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
!
! Note: This routine does not do any other bookkeeping. Consider using
! control_bookkeeper or lattice_bookkeeper instead.
! 
! BEAMBEAM:   
!     bbi_const$ = param%n_part * charge$ * classical_radius_factor /
!                           (2 * pi * p0c$ * (sig_x$ + sig_y$)
!
! CRYSTAL:
!     graze_angle_in$
!     graze_angle_out$
!     tilt_corr$
!
! ELSEPARATOR:
!     e_field$ = sqrt(hkick$**2 + vkick$**2) * p0c$ / l$
!     voltage$ = e_field$ * gap$ 
!
! LCAVITY:    
!     voltage$ = gradient$ * L$ 
! 
! RFCAVITY:   
!     rf_frequency$ = harmon$ * c_light / param%total_length (only if harmon$ /= 0)
!
! SBEND:      
!     angle$   = L$ * G$
!     l_chord$ = 2 * sin(Angle$/2) / G$
!     rho$     = 1 / G$
!
! WIGGLER (map_type):
!     B_MAX$    
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!
! WIGGLER (periodic_type):
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!     n_pole$ = L$ / l_pole$
!
! Modules needed:
!   use bmad
!
! Input:
!   ele            -- Ele_struct: Element with attributes 
!   param          -- lat_param_struct: 
!   force_bookkeeping 
!                  -- Logical, optional: If present and True then force
!                       attribute bookkeeping to be done independent of
!                       the state of ele%bookkeeping_stat%attributes.
! Output:
!   ele            -- Ele_struct: Element with self-consistant attributes.
!
! Programming Note: If the dependent attributes are changed then 
!       the attribute_free routine must be modified.
!-

subroutine attribute_bookkeeper (ele, param, force_bookkeeping)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) start, end
type (em_field_struct) field

real(rp) factor, gc, f2, phase, E_tot, polarity, dval(num_ele_attrib$)
real(rp), pointer :: val(:)

integer i, n

character(20) ::  r_name = 'attribute_bookkeeper'

logical, optional :: force_bookkeeping
logical err_flag
logical non_offset_changed, offset_changed, offset_nonzero, is_on
logical, save :: v_mask(num_ele_attrib$), offset_mask(num_ele_attrib$)
logical :: init_needed = .true.

! Some init

val => ele%value

! Intelligent bookkeeping

if (.not. bmad_com%auto_bookkeeper) then
  if (ele%bookkeeping_state%attributes /= stale$ .and. &
      .not. logic_option(.false., force_bookkeeping)) return

  if (ele%lord_status /= not_a_lord$) then
    call set_ele_status_stale (ele, control_group$)
  endif

  if (ele%old_value(l$) /= val(l$)) then
    call set_ele_status_stale (ele, s_position_group$)
  endif

  if (ele%lord_status /= overlay_lord$ .and. ele%lord_status /= group_lord$ .and. &
      ele%lord_status /= multipass_lord$) then
    call set_ele_status_stale (ele, mat6_group$)
  endif

endif

ele%bookkeeping_state%attributes = ok$
ele%bookkeeping_state%rad_int = stale$

! For auto bookkeeping if no change then we don't need to do anything

if (bmad_com%auto_bookkeeper) then

  val(check_sum$) = 0
  if (associated(ele%a_pole)) val(check_sum$) = sum(ele%a_pole) + sum(ele%b_pole)

  dval = val - ele%old_value
  dval(x1_limit$:y2_limit$) = 0  ! Limit changes do not need bookkeeping
  dval(scratch$) = 0
  if (all(dval == 0) .and. ele%key /= capillary$) return
endif

! Transfer tilt to tilt_tot, etc.

if (.not. ele%on_a_girder .and. has_orientation_attributes(ele)) then
  val(tilt_tot$)     = val(tilt$)
  val(x_offset_tot$) = val(x_offset$)
  val(y_offset_tot$) = val(y_offset$)
  val(s_offset_tot$) = val(s_offset$)
  val(x_pitch_tot$)  = val(x_pitch$)
  val(y_pitch_tot$)  = val(y_pitch$)
endif

! Taylor elements need no more bookkeeping

if (ele%key == taylor$) then
  ele%old_value = val
  return
endif

! Field_master...

if (ele%field_master) then

  if (val(p0c$) == 0) then
    factor = 0
  else
    factor = c_light / val(p0c$)
  endif

  select case (ele%key)
  case (quadrupole$)
    val(k1$) = factor * val(B1_gradient$)
  case (sextupole$)
    val(k2$) = factor * val(B2_gradient$)
  case (octupole$)
    val(k3$) = factor * val(B3_gradient$)
  case (solenoid$)
    val(ks$) = factor * val(Bs_field$)
  case (sol_quad$)
    val(ks$) = factor * val(Bs_field$)
    val(k1$) = factor * val(B1_gradient$)
  case (sbend$)
    val(g$)     = factor * val(B_field$)
    val(g_err$) = factor * val(B_field_err$)
    val(k1$)    = factor * val(B1_gradient$)
    val(k2$)    = factor * val(B2_gradient$)
  case (hkicker$)
    val(kick$) = factor * val(BL_kick$)
  case (vkicker$)
    val(kick$) = factor * val(BL_kick$)
  case (bend_sol_quad$)
    val(g$)  = factor * val(B_field$)
    val(k1$) = factor * val(B1_gradient$)
    val(ks$) = factor * val(Bs_field$)
  end select

  if (has_hkick_attributes(ele%key)) then
    val(hkick$) = factor * val(BL_hkick$)
    val(vkick$) = factor * val(BL_vkick$)
  endif

else

  factor = val(p0c$) / c_light

  select case (ele%key)
  case (quadrupole$)
    val(B1_gradient$) = factor * val(k1$)
  case (sextupole$)
    val(B2_gradient$) = factor * val(k2$)
  case (octupole$)
    val(B3_gradient$) = factor * val(k3$)
  case (solenoid$)
    val(Bs_field$)    = factor * val(ks$)
  case (sol_quad$)
    val(Bs_field$)    = factor * val(ks$)
    val(B1_gradient$) = factor * val(k1$)
  case (sbend$)
    val(B_field$)     = factor * val(g$)
    val(B_field_err$) = factor * val(g_err$)
    val(B1_gradient$) = factor * val(k1$)
    val(B2_gradient$) = factor * val(k2$)
  case (hkicker$)
    val(BL_kick$) = factor * val(kick$)
  case (vkicker$) 
    val(BL_kick$) = factor * val(kick$)
  case (bend_sol_quad$)
    val(B_field$)     = factor * val(g$)
    val(B1_gradient$) = factor * val(k1$)
    val(Bs_field$)    = factor * val(ks$)
  end select

  if (has_hkick_attributes(ele%key)) then
    val(BL_hkick$) = factor * val(hkick$)
    val(BL_vkick$) = factor * val(vkick$)
  endif

endif

! Dependent attribute bookkeeping.
! Note: If the dependent attributes are changed then attribute_free 
!       must be modified.

! num_steps

if (val(ds_step$) <= 0) then
  if (val(num_steps$) <= 0 .or. abs(val(l$)) == 0) then
    val(ds_step$) = bmad_com%default_ds_step
  else
    val(ds_step$) = val(num_steps$) * abs(val(l$))
  endif
endif
 
val(num_steps$) = max(1, nint(abs(val(l$) / val(ds_step$))))

!----------------------------------
! General bookkeeping...

select case (ele%key)

! Bends

case (sbend$)

  val(angle$) = val(l$) * val(g$)

  if (val(l$) == 0 .or. val(g$) == 0) then
    val(l_chord$) = 0
  else
    val(l_chord$) = 2 * sin(val(angle$)/2) / val(g$)
  endif

  if (val(g$) == 0) then
    val(rho$) = 0
  else
    val(rho$) = 1 / val(g$)
  endif

  if (ele%old_value(g$) /= val(g$)) then
    call set_ele_status_stale (ele, floor_position_group$)
  endif

! Lcavity

case (lcavity$)
  if (ele%lord_status /= multipass_lord$) then
    if (val(phi0$) /= ele%old_value(phi0$) .or. val(dphi0$) /= ele%old_value(dphi0$) .or. &
        val(gradient$) /= ele%old_value(gradient$) .or. val(e_loss$) /= ele%old_value(e_loss$) .or. &
        val(l$) /= ele%old_value(l$)) then
      call set_ele_status_stale (ele, ref_energy_group$)
    endif
  endif

  val(voltage$) = val(gradient$) * val(l$)
  if (val(rf_frequency$) == 0) then
    val(l_hard_edge$) = 0
  else
    val(l_hard_edge$) = c_light * nint(val(n_cell$)) / (2 * val(rf_frequency$))
  endif

! E_Gun

case (e_gun$)
  if (ele%lord_status /= multipass_lord$) then
    if (val(gradient$) /= ele%old_value(gradient$) .or. val(l$) /= ele%old_value(l$)) then
      call set_ele_status_stale (ele, ref_energy_group$)
      val(voltage$) = val(gradient$) * val(l$)
      val(voltage_err$) = val(gradient_err$) * val(l$)
    endif
    val(l_hard_edge$) = val(l$)
  endif


! RFcavity

case (rfcavity$)
  if (val(harmon$) /= 0) val(rf_frequency$) =  val(harmon$) * c_light / param%total_length 

  if (val(rf_frequency$) == 0) then
    val(l_hard_edge$) = 0
  else
    val(l_hard_edge$) = c_light * nint(val(n_cell$)) / (2 * val(rf_frequency$))
  endif

! Solenoid

case (solenoid$)
  val(l_hard_edge$) = val(l$)

! BeamBeam

case (beambeam$)

  if (val(n_slice$) == 0) val(n_slice$) = 1.0 ! revert to default

  if (val(charge$) == 0 .or. param%n_part == 0) then
    val(bbi_const$) = 0

  else

    if (val(sig_x$) == 0 .or. val(sig_y$) == 0) then
      call out_io(s_abort$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!')
      call type_ele(ele, .true., 0, .false., 0, .false.)
      if (global_com%exit_on_error) call err_exit
    endif

    val(bbi_const$) = -param%n_part * val(charge$) * classical_radius_factor /  &
                             (2 * pi * val(p0c$) * (val(sig_x$) + val(sig_y$)))

  endif

! Crystal

case (crystal$, multilayer_mirror$)

  if (ele%key == crystal$) then
    call crystal_type_to_crystal_params (ele, err_flag)
    call crystal_attribute_bookkeeper (ele)
  else
    call multilayer_type_to_multilayer_params (ele, err_flag)
  endif

  gc = 0
  if (val(d_source$) /= 0) gc = val(graze_angle_in$) / (2 * val(d_source$))
  if (val(d_detec$) /= 0)  gc = gc + val(graze_angle_out$) / (2 * val(d_detec$))

  val(c2_curve_tot$) = val(c2_curve$) + gc / 2
  val(c3_curve_tot$) = val(c3_curve$)
  val(c4_curve_tot$) = val(c4_curve$) + gc**3 / 4

! Elseparator

case (elseparator$)

  if (val(l$) == 0 .or. val(gap$) == 0) then
    val(e_field$) = 0
    val(voltage$) = 0
  else
    val(e_field$) = sqrt(val(hkick$)**2 + val(vkick$)**2) * val(p0c$) / val(l$)
    val(voltage$) = val(e_field$) * val(gap$) 
  endif


! Wiggler

case (wiggler$) 

  ! Calculate b_max for map_type wigglers. 

  if (ele%sub_key == map_type$ .and. val(b_max$) == 0) then
    is_on = ele%is_on  ! Save
    polarity = val(polarity$)
    ele%is_on = .true.
    val(polarity$) = 1
    start%vec = 0
    val(b_max$) = 0
    n = nint(val(num_steps$))
    do i = 0, n
      call em_field_calc (ele, param, i * val(l$) / n, 0.0_rp, start, .true., field)
      val(b_max$) = max(val(b_max$), sqrt(sum(field%b**2)))
    enddo
    ele%is_on = is_on
    val(polarity$) = polarity
  endif

  if (val(p0c$) == 0) then
    val(k1$) = 0
  else
    val(k1$) = -0.5 * (c_light * val(b_max$) / val(p0c$))**2
  endif

  if (val(b_max$) == 0) then
    val(rho$) = 0
  else
    val(rho$) = val(p0c$) / (c_light * val(b_max$))
  endif

  if (val(l_pole$) == 0) then
    val(n_pole$) = 0
  else
    val(n_pole$) = val(l$) / val(l_pole$)
  endif

  ! Periodic_type wigglers have a single %wig%term for use with tracking, etc.
  ! The phase of this term is set so that tracking with a particle starting
  ! on-axis ends on-axis. For this to be true, there must be an integer number
  ! of poles.

  ! For super_slave and sliced elements, the phi_z is set by the position with respect to the lord in
  ! the routine makeup_super_slave1 and so should not be touched here.

  if (ele%sub_key == periodic_type$ .and. ele%slave_status /= super_slave$ .and. &
      ele%slave_status /= multipass_slave$ .and. ele%slave_status /= slice_slave$) then
    if (.not. associated(ele%wig)) then
      allocate (ele%wig)
      allocate (ele%wig%term(1))
    endif

    if (val(l$) == 0) then
      ele%wig%term(1)%ky = 0
    else
      ele%wig%term(1)%ky = pi * val(n_pole$) / val(l$)
    endif
    ele%wig%term(1)%coef   = val(b_max$)
    ele%wig%term(1)%kx     = 0
    ele%wig%term(1)%kz     = ele%wig%term(1)%ky
    ele%wig%term(1)%phi_z  = -ele%wig%term(1)%kz * val(l$) / 2 
    ele%wig%term(1)%type   = hyper_y$
  endif

end select

! It might be true that value == old_value here even though this was not true at the start
! of this routine. For example: A super_slave, via makeup_super_slave (which will be
! called by control_bookkeeper1 before calling this routine), will inherit its lords 
! num_steps value but then this routine will reset num_steps to be the correct value.

! So stop here if nothing has truely changed.

if (bmad_com%auto_bookkeeper) then
  if (all(val == ele%old_value)) return
endif

! Since things have changed we need to kill the Taylor Map and ptc_genfield.

if (init_needed) then
  v_mask = .true.
  v_mask([x_offset$, y_offset$, s_offset$, &
            tilt$, x_pitch$, y_pitch$, x_offset_tot$, y_offset_tot$, s_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$]) = .false.
  offset_mask = .not. v_mask
  v_mask( [x1_limit$, x2_limit$, y1_limit$, y2_limit$] ) = .false.
  init_needed = .false.
endif

dval = val - ele%old_value
dval(scratch$) = 0

if (has_orientation_attributes(ele)) then
  non_offset_changed = (any(dval /= 0 .and. v_mask))
  offset_changed =  (any(dval /= 0 .and. offset_mask))
  offset_nonzero = (any(val /= 0 .and. offset_mask))
else
  non_offset_changed = (any(dval /= 0))
  offset_changed = .false.
  offset_nonzero = .false.
endif

! If an element has just been offset and bmad_com%conserve_taylor_map = T then 
! conserve the taylor map.

if (associated(ele%taylor(1)%term) .and. ele%map_with_offsets .and. &
        offset_nonzero .and. offset_changed .and. .not. non_offset_changed .and. &
        bmad_com%conserve_taylor_maps .and. ele%key /= patch$) then
  ele%map_with_offsets = .false.
  if (associated(ele%branch) .and. ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
      lord%map_with_offsets = .false.
    enddo
  endif
  if (any(ele%old_value /= 0 .and. offset_mask)) non_offset_changed = .true.  ! To trigger kill_taylor below
  call out_io (s_info$, r_name, &
      'Note: bmad_com%conserve_taylor_maps = True (this is the default)', &
      'But: Element has just been offset: ' // ele%name, &
      "To conserve the element's Taylor map, I will set ele%map_with_offsets = False.")
endif

! Kill the taylor map and ptc_genfield if necessary.

if (non_offset_changed .or. (offset_changed .and. ele%map_with_offsets)) then
  if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)
  if (associated(ele%ptc_genfield)) call kill_ptc_genfield(ele%ptc_genfield)
endif

! Make stale ele%rad_int_cache if allocated

if (associated(ele%rad_int_cache)) ele%rad_int_cache%stale = .true.  ! Forces recalc

! Set old_value = value

ele%old_value = val

end subroutine attribute_bookkeeper

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_integer_attribute (lat, ele, attrib)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_integer_attribute (lat, ele, attrib)

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch

integer, target :: attrib
integer, pointer :: a_ptr
real(rp) dummy

!

call set_flags_for_changed_real_attribute (lat, ele, dummy)

!

branch => lat%branch(ele%ix_branch)
a_ptr => attrib

if (ele%value(p0c$) /= ele%value(p0c_start$)) then
  if (associated(a_ptr, ele%tracking_method) .or. associated(a_ptr, ele%field_calc)) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif
endif

end subroutine set_flags_for_changed_integer_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_logical_attribute (lat, ele, attrib)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_logical_attribute (lat, ele, attrib)

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch

logical, target :: attrib
logical, pointer :: a_ptr
real(rp) dummy

!

call set_flags_for_changed_real_attribute (lat, ele, dummy)

end subroutine set_flags_for_changed_logical_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_real_attribute (lat, ele, attrib)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_real_attribute (lat, ele, attrib)

implicit none

type (lat_struct), target :: lat
type (ele_struct), optional, target :: ele
type (branch_struct), pointer :: branch
type (em_field_mode_struct), pointer :: mode

real(rp), optional, target :: attrib
real(rp), pointer :: a_ptr
real(rp) v_mat(4,4), v_inv_mat(4,4), eta_vec(4), eta_xy_vec(4)
real(rp), target :: unknown_attrib

integer i, j, ib

logical coupling_change, found

! If ele is not present then must reinit everything.

if (.not. present(ele)) then
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    call set_status_flags (branch%param%bookkeeping_state, stale$)
    do j = 0, ubound(branch%ele, 1)
      call set_status_flags (branch%ele(j)%bookkeeping_state, stale$)
    enddo
  enddo
  return
endif

!-------------------
! For a particular elemement...

branch => lat%branch(ele%ix_branch)

! If a lord then set the control flag stale

if (ele%lord_status /= not_a_lord$) call set_ele_status_stale (ele, control_group$)

! Groups and overlays do not have any dependent attributes. 
! For all others set the attributes flag stale.

if (ele%key /= group$ .and. ele%key /= overlay$) then
  call set_ele_status_stale (ele, attribute_group$)
endif

! Transfer matrix calc needs to be flagged

if (ele%lord_status /= overlay_lord$ .and. ele%lord_status /= group_lord$ .and. &
    ele%lord_status /= multipass_lord$) then
  call set_ele_status_stale (ele, mat6_group$)
endif

! If attrib is not present then point to a dummy location which will not match when 
! the associated() function is used below.

if (.not. present(attrib)) then
  call set_ele_status_stale (ele, all_groups$)
endif

! Use a_ptr with the associated function to see which attribute has been changed.

a_ptr => attrib

! A limit change does not need any bookkeeping

if (associated(a_ptr, ele%value(x1_limit$)) .or. associated(a_ptr, ele%value(x2_limit$)) .or. &
    associated(a_ptr, ele%value(y1_limit$)) .or. associated(a_ptr, ele%value(y2_limit$))) return

! A length change involves changes in the floor position.

if (associated(a_ptr, ele%value(l$))) then
  if (ele%lord_status /= overlay_lord$ .and. ele%lord_status /= group_lord$) then
    call set_ele_status_stale (ele, s_position_group$)
    call set_ele_status_stale (ele, floor_position_group$)
  endif
  if (ele%value(p0c$) /= ele%value(p0c_start$)) call set_ele_status_stale (ele, ref_energy_group$)
endif

! E_tot and p0c can be varied in an init_ele or a multipass lord with n_ref_pass = 0.
! In addition, for an init_ele, must also set e_tot_start and p0c_start. This is important
! for lattices with an e_gun element

if (associated(a_ptr, ele%value(e_tot$))) then
  call convert_total_energy_to (ele%value(e_tot$), branch%param%particle, pc = ele%value(p0c$))
  call set_ele_status_stale (ele, ref_energy_group$)
  if (ele%key == init_ele$) then
    ele%value(e_tot_start$) = ele%value(e_tot$)
    ele%value(p0c_start$) = ele%value(p0c$)
  endif
  return
endif

if (associated(a_ptr, ele%value(p0c$))) then
  call convert_pc_to (ele%value(p0c$), branch%param%particle, e_tot = ele%value(e_tot$))
  call set_ele_status_stale (ele, ref_energy_group$)
  if (ele%key == init_ele$) then
    ele%value(e_tot_start$) = ele%value(e_tot$)
    ele%value(p0c_start$) = ele%value(p0c$)
  endif
  return
endif

if (associated(a_ptr, ele%value(e_tot_start$))) then
  call convert_total_energy_to (ele%value(e_tot_start$), branch%param%particle, pc = ele%value(p0c_start$))
  call set_ele_status_stale (ele, ref_energy_group$)
  return
endif

if (associated(a_ptr, ele%value(p0c_start$))) then
  call convert_pc_to (ele%value(p0c_start$), branch%param%particle, e_tot = ele%value(e_tot_start$))
  call set_ele_status_stale (ele, ref_energy_group$)
  return
endif

!

select case (ele%key)
case (init_ele$) 
  coupling_change = .false.

  if (associated(a_ptr, ele%a%beta) .or. associated(a_ptr, ele%a%alpha)) then
    if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
    return
  endif

  if (associated(a_ptr, ele%b%beta) .or. associated(a_ptr, ele%b%alpha)) then
    if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta
    return
  endif

  if (associated(a_ptr, ele%c_mat(1,1)) .or. associated(a_ptr, ele%c_mat(1,2)) .or. & 
          associated(a_ptr, ele%c_mat(2,1)) .or. associated(a_ptr, ele%c_mat(2,2))) then
    ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + ele%c_mat(1,2)*ele%c_mat(2,1))
    coupling_change = .true.
  endif

  if (associated(a_ptr, ele%x%eta) .or. associated(a_ptr, ele%x%etap) .or. &
      associated(a_ptr, ele%y%eta) .or. associated(a_ptr, ele%y%etap) .or. &
      coupling_change) then 
    call make_v_mats (ele, v_mat, v_inv_mat)
    eta_xy_vec = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap]
    eta_vec = matmul (v_inv_mat, eta_xy_vec)
    ele%a%eta  = eta_vec(1)
    ele%a%etap = eta_vec(2)
    ele%b%eta  = eta_vec(3)
    ele%b%etap = eta_vec(4)
    return
  endif

  if (associated(a_ptr, ele%a%eta) .or. associated(a_ptr, ele%a%etap) .or. &
      associated(a_ptr, ele%b%eta) .or. associated(a_ptr, ele%b%etap)) then 
    call make_v_mats (ele, v_mat, v_inv_mat)
    eta_vec = [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap]
    eta_xy_vec = matmul (v_mat, eta_vec)
    ele%x%eta  = eta_xy_vec(1)
    ele%x%etap = eta_xy_vec(2)
    ele%y%eta  = eta_xy_vec(3)
    ele%y%etap = eta_xy_vec(4)
    return
  endif

  if (associated(a_ptr, ele%floor%x) .or. associated(a_ptr, ele%floor%y) .or. &
      associated(a_ptr, ele%floor%z) .or. associated(a_ptr, ele%floor%theta) .or. &
      associated(a_ptr, ele%floor%phi) .or. associated(a_ptr, ele%floor%psi)) then
    call set_ele_status_stale (ele, floor_position_group$)
    return
  endif

case (sbend$)
  if (associated(a_ptr, ele%value(angle$)) .or. associated(a_ptr, ele%value(g$)) .or. &
      associated(a_ptr, ele%value(rho$)) .or. associated(a_ptr, ele%value(b_field$))) then
    call set_ele_status_stale (ele, floor_position_group$)
  endif

case (branch$, photon_branch$)
  if (associated(a_ptr, ele%value(direction$))) then
    branch%param%bookkeeping_state%floor_position = stale$
    branch%ele(0)%bookkeeping_state%floor_position = stale$
  endif

case (lcavity$, e_gun$)
  if (associated(a_ptr, ele%value(gradient$)) .or. associated(a_ptr, ele%value(phi0$)) .or. &
      associated(a_ptr, ele%value(dphi0$)) .or. associated(a_ptr, ele%value(e_loss$))) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

  if (associated(ele%em_field)) then
    found = .false.
    do i = 1, size(ele%em_field%mode)
      mode => ele%em_field%mode(i)
      if (associated(a_ptr, mode%dphi0_ref)) found = .true.
      if (associated(a_ptr, mode%field_scale)) found = .true.
      if (mode%master_scale > 0) found = found .or. (associated(a_ptr, ele%value(mode%master_scale)))
      if (associated(a_ptr, mode%dphi0_ref)) found = .true.
    enddo
    if (found) call set_ele_status_stale (ele, ref_energy_group$)
  endif

case (patch$)
  if (associated(a_ptr, ele%value(e_tot_offset$))) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

end select

!

end subroutine set_flags_for_changed_real_attribute

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a lattice. An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6.
!
! NOTE: When set_on_off is used with switch = save_state$ the state
! of an element, ele%is_on, is saved in ele%old_is_on. That is, the saved
! state "history" is only one level deep. A problem arises when a routine 
! uses set_on_off with switch = save_state$ and then calls a routine which 
! also calls set_on_off switch = save_state$. For example:
!     call set_on_off (rfcavity$, lat, save_state$)
!     call set_on_off (rfcavity$, lat, off$)
!     call some_subroutine (lat)
!     call set_on_off (rfcavity$, lat, restore_state$)
! If some_subroutine also calls set_on_off with switch = save_state$ then
! ele%old_is_on will be set to False erasing the saved state. The subsequent
! call to set_on_off with switch = restore_state$ will switch all cavities off.
! It is thus important to avoid this situation. To help avoid this situation,
! All Bmad routines in the Bmad library will make sure to preserve ele%old_is_on.
!
! Modules needed:
!   use bmad
!
! Input:
!   key          -- Integer: Key name of elements to be turned on or off.
!                      [Key = quadrupole$, etc.]
!   lat          -- lat_struct: lattice structure holding the elements
!   switch       -- Integer: 
!                     on$            => Turn elements on.  
!                     off$           => Turn elements off. 
!                     save_state$    => Save present on/off state. 
!                                         No turning on or off is done.
!                     restore_state$ => Restore saved on/off state.
!   orb(0:)     -- Coord_struct, optional: Needed for lat_make_mat6
!   use_ref_orb -- Logical, optional: If present and true then use the
!                    present ele%map_ref_orb for the reference orbit for
!                    constructing %mat6. Default is false.
!   ix_branch   -- integer, optional: If present then only set for 
!                    this lattice branch.
!
! Output:
!   lat -- lat_struct: Modified lattice.
!-

subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), optional :: orb(0:)
type (coord_struct) ref_orb

integer i, ib, key
integer, intent(in) :: switch
integer, optional :: ix_branch

logical, optional :: use_ref_orb
logical old_state

character(20) :: r_name = 'set_on_off'

!

do ib = 0, ubound(lat%branch, 1)

  if (present(ix_branch)) then
    if (ix_branch /= ib) cycle
  endif

  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max

    if (branch%ele(i)%key /= key) cycle

    old_state = branch%ele(i)%is_on

    select case (switch)
    case (on$)
      branch%ele(i)%is_on = .true.
    case (off$)
      branch%ele(i)%is_on = .false.
    case (save_state$)
      branch%ele(i)%old_is_on = branch%ele(i)%is_on
      cycle
    case (restore_state$)
      branch%ele(i)%is_on = branch%ele(i)%old_is_on
    case default
      call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
      if (global_com%exit_on_error) call err_exit
    end select

    if (old_state .neqv. branch%ele(i)%is_on) then
      if (logic_option (.false., use_ref_orb)) then
        ref_orb%vec = branch%ele(i)%map_ref_orb_in
        call make_mat6(branch%ele(i), branch%param, ref_orb)
      else
        call set_ele_status_stale (branch%ele(i), mat6_group$)
        call lat_make_mat6(lat, i, orb, ib)
      endif
    endif

  enddo
enddo

end subroutine set_on_off

end module
