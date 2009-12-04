module bookkeeper_mod

use bmad_interface
use bmad_utils_mod
use multipole_mod
use lat_geometry_mod

integer, parameter :: off$ = 1, on$ = 2
integer, parameter :: save_state$ = 3, restore_state$ = 4

private control_bookkeeper1, makeup_overlay_and_girder_slave 
        
contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lattice_bookkeeper (lat)
!
! Subroutine to do a complete bookkeeping job on a lattice.
!
! This this routine does a complete job of bookking and could be unacceptably
! slow if used, for example, in the inner loop of an optimizer. In this case
! consider using only control_bookkeeper instead.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice needing bookkeeping.
!
! Output:
!   lat   -- lat_struct: Lattice with bookkeeping done.
!-

subroutine lattice_bookkeeper (lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
integer i, j
logical found

! Control bookkeeper is called twice to make sure, for example, that the z_patch for a 
! wiggler super_lord is computed. Other reasons include multipass bends.

call control_bookkeeper (lat)
call compute_reference_energy (lat)
call control_bookkeeper (lat, super_and_multipass_only = .true.)

! Make sure attributes are updated for all elements in the tracking part of the lattice.
! We do this since control_bookkeeper does not take care of free elements.
! Also overlay slaves with field_master = T need doing...

do i = 0, ubound(lat%branch, 1)
  do j = 1, lat%branch(i)%n_ele_track
    call attribute_bookkeeper (lat%branch(i)%ele(j), lat%param)
  enddo
enddo

! Global geometry

call s_calc (lat)
call lat_geometry (lat)

! multipass slaves with ref_orbit set may may depend upon the geometry so recalc.

found = .false.

do i = 0, ubound(lat%branch, 1)
  do j = 1, lat%branch(i)%n_ele_track
    ele => lat%branch(i)%ele(j)
    if (ele%slave_status == multipass_slave$ .and. ele%ref_orbit /= 0) then
      call makeup_multipass_slave (lat, ele)
      call attribute_bookkeeper (ele, lat%param)
      found = .true.
    endif
  enddo
enddo

if (found) then
  call s_calc (lat)
  call lat_geometry (lat)
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (lat, ix_ele, ix_branch, super_and_multipass_only)
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
!   lat       -- lat_struct: lattice to be used
!   ix_ele    -- Integer, optional: Index of element whose attribute values 
!                  have been changed. If not present bookkeeping will be done 
!                  for all elements.
!   ix_branch -- Integer, optional: Branch index. Default is 0.
!   super_and_multipass_only 
!             -- Logical, optional: If True then only do bookkeeping for 
!                  superposition and multipass elements only. Default is False. 
!                  This argument is used by lattice_bookkeeper and should not
!                  to set unless you know what you are doing.
!-

recursive subroutine control_bookkeeper (lat, ix_ele, ix_branch, super_and_multipass_only)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, lord

integer, optional :: ix_ele, ix_branch
integer ie, ib, ix0, j, n1, n2

logical, optional :: super_and_multipass_only
logical sm_only, did_bookkeeping

! Check that super_slave lengths add up to the super_lord_length.
! With super_only (used by lattice_bookkeeper) this check has already been
! done so we don't need to do it again.

ix0 = integer_option(0, ix_ele)
ib = integer_option (0, ix_branch)
sm_only = logic_option (.false., super_and_multipass_only)
ele => lat%branch(ib)%ele(ix0)

if (.not. sm_only) then
  if (.not. present(ix_ele) .or. ele%lord_status == super_lord$) &
                                    call super_lord_length_bookkeeper (lat, ix_ele)
endif

! If ix_ele is present we only do bookkeeping for this one element and its slaves

if (present(ix_ele)) then
  call control_bookkeeper1 (lat, ele, sm_only)
  do ie = 1, ele%n_slave
    slave => pointer_to_slave (lat, ele, ie)
    call control_bookkeeper (lat, slave%ix_ele, slave%ix_branch, sm_only)
  enddo
  return
endif

! Else we need to make up all the lords. 
! Need to do this from the top level down.
! The top level are those lord elements that have no lords.

n1 = lat%n_ele_track+1
n2 = lat%n_ele_max
lat%ele(n1:n2)%bmad_logic = .false.  ! Bookkeeping done on this element yet?

do
  did_bookkeeping = .true.
  ie_loop: do ie = n1, n2
    if (lat%ele(ie)%bmad_logic) cycle
    do j = 1, lat%ele(ie)%n_lord
      lord => pointer_to_lord (lat, lat%ele(ie), j)
      if (.not. lord%bmad_logic) then
        did_bookkeeping = .false.  ! This element remains to be done.
        cycle ie_loop ! Do not do bookkeeping yet if lord not done yet.
      endif
    enddo
    call control_bookkeeper1 (lat, lat%ele(ie), sm_only)
    lat%ele(ie)%bmad_logic = .true.  ! Done this element
  enddo ie_loop
  if (did_bookkeeping) exit  ! And we are done
enddo

! and now the slaves in the tracking lattice

do ib = 0, ubound(lat%branch, 1)
  do ie = 1, lat%branch(ib)%n_ele_track
    if (lat%branch(ib)%ele(ie)%slave_status /= free$) then
      call control_bookkeeper1 (lat, lat%ele(ie), sm_only)
    endif
  enddo
enddo

end subroutine

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

logical sm_only, called_a_bookkeeper

! Init

if (sm_only) then
  if (ele%lord_status /= super_lord$ .and. ele%slave_status /= super_slave$ .and. &
      ele%lord_status /= multipass_lord$ .and. ele%slave_status /= multipass_slave$) return
endif

! Make sure the bookkeeping for this element is correct.

call attribute_bookkeeper (ele, lat%param)

! Slave bookkeeping

called_a_bookkeeper = .false.

if (ele%slave_status == super_slave$) then
  call makeup_super_slave (lat, ele)
  called_a_bookkeeper = .true.

elseif (ele%slave_status == overlay_slave$) then
  call makeup_overlay_and_girder_slave (lat, ele)
  called_a_bookkeeper = .true.

elseif (ele%slave_status == multipass_slave$) then
  call makeup_multipass_slave (lat, ele)
  if (ele%n_lord > 1) call makeup_overlay_and_girder_slave (lat, ele)
  called_a_bookkeeper = .true.
endif

! Lord bookkeeping

if (ele%lord_status == group_lord$) then
  call makeup_group_lord (lat, ele)
  called_a_bookkeeper = .true.

elseif (ele%lord_status == super_lord$) then
  call adjust_super_lord_s_position (lat, ele)
  called_a_bookkeeper = .true.

endif

! 

if (called_a_bookkeeper) call attribute_bookkeeper (ele, lat%param)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine super_lord_length_bookkeeper (lat, ix_ele)
!
! Subroutine to make sure the length of the slaves of a super_lord add up to the
! length of the lord. If not, make an adjustment to the slave length.
!
! Note: This routine is called by control_bookkeeper. Generally this
! routine does not have to be called directly.
!
! Modules needed:
!   use bookkeeper_mod
!
! Input:
!   lat    -- Lat_struct: Lattice.
!   ix_ele -- Integer, optional: Index of super_lord element to check.
!                  If not present bookkeeping will be done for all super_lords.
!
! Output:
!   lat  -- Lat_struct: Lattice with adjustments made.
!-

subroutine super_lord_length_bookkeeper (lat, ix_ele)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: lord0, lord2, slave, slave2
type (branch_struct), pointer :: branch

real(rp) length, coef, length_start, length_pos, length_neg
real(rp) d_length, d_length_pos, d_length_neg
real(rp) dl_tol

integer, optional :: ix_ele
integer j, k, ie, ix0, ixa_lord0, ixb_lord0, ixa_lord2, ixb_lord2
integer ix_pos_edge, ix_neg_edge, ixa, ixb

logical pos_extension_lord_exists, neg_extension_lord_exists, all_extension_lord_exists
logical length_adjustment_made, overlap_a, overlap_b, overlap_all

character(40) :: r_name = 'super_lord_length_bookkeeper'

!

dl_tol = 10 * bmad_com%significant_longitudinal_length

length_adjustment_made = .false.
ix0 = integer_option(0, ix_ele)

do ie = lat%n_ele_track+1, lat%n_ele_max

  lord0 => lat%ele(ie)

  if (lord0%lord_status /= super_lord$) cycle
  if (ix0 /= 0 .and. ix0 /= ie) cycle

  length = 0
  do j = 1, lord0%n_slave
    slave => pointer_to_slave(lat, lord0, j)
    length = length + slave%value(l$)
  enddo

  ! Nothing to be done if the lengths add up.

  if (abs(length - lord0%value(l$)) < dl_tol * (1 + lord0%value(l$))) cycle

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

  slave = pointer_to_slave(lat, lord0, 1)
  ixa_lord0 = slave%ix_ele  ! Index at entrance end of lord0

  slave = pointer_to_slave(lat, lord0, lord0%n_slave)
  ixb_lord0 = slave%ix_ele  ! Index at exit end of lord0

  pos_extension_lord_exists = .false.
  neg_extension_lord_exists = .false.
  all_extension_lord_exists = .false.

  ix_pos_edge = lat%n_ele_max
  ix_neg_edge = 0

  length_start = 0
  slave_loop: do j = 1, lord0%n_slave

    slave => pointer_to_slave(lat, lord0, j)
    slave%bmad_logic = .true. ! Can be varied.

    do k = 1, slave%n_lord
      lord2 => pointer_to_lord(lat, slave, k)
      if (lord0%ix_ele == lord2%ix_ele) cycle  ! Ignore self
      slave2 => pointer_to_slave(lat, lord2, 1)            ! Slave at entrance end
      ixa_lord2 = slave2%ix_ele

      slave2 => pointer_to_slave(lat, lord2, lord2%n_slave) ! Slave at exit end
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
        slave%bmad_logic = .false.
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

    if (slave%bmad_logic) length_start = length_start + slave%value(l$)
  enddo slave_loop

  ! If we have not found any slaves to vary we are in trouble

  if (length_start == 0) then
    call out_io (s_fatal$, r_name, 'CANNOT VARY LENGTH OF SUPER_LORD: ' // lord0%name)
    call err_exit
  endif

  ! Calculate positive and negative extension length changes

  coef = lord0%value(l$) / length_start
  d_length = lord0%value(l$) - length_start

  if (pos_extension_lord_exists) then
    length_pos = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lat, lord0, j)
      if (slave%ix_ele < ix_pos_edge) cycle
      if (.not. slave%bmad_logic) cycle
      length_pos = length_pos + slave%value(l$)
    enddo
    d_length_pos = length_pos * (coef - 1)
  endif

  if (neg_extension_lord_exists) then
    length_neg = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lat, lord0, j)
      if (slave%ix_ele > ix_neg_edge) cycle
      if (.not. slave%bmad_logic) cycle
      length_neg = length_neg + slave%value(l$)
    enddo    
    d_length_neg = length_neg * (coef - 1)
  endif

  ! Vary the slave lengths

  do j = 1, lord0%n_slave
    slave => pointer_to_slave(lat, lord0, j)
    if (.not. slave%bmad_logic) cycle
    slave%value(l$) = slave%value(l$) * coef
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
    if (pos_extension_lord_exists .and. neg_extension_lord_exists) then
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length_pos
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length_neg
    elseif (pos_extension_lord_exists) then
      branch%ele(ixb)%value(l$) = branch%ele(ixb)%value(l$) - d_length_pos
    elseif (neg_extension_lord_exists) then
      branch%ele(ixa)%value(l$) = branch%ele(ixa)%value(l$) - d_length_neg
    else  
      ! No extensions -> Nothing to be done.
    endif
  endif

enddo

! If there has been a length adjustment then we need to make sure everything is ok.

if (length_adjustment_made) then
  do ie = lat%n_ele_track+1, lat%n_ele_max
    lord0 => lat%ele(ie)
    if (lord0%lord_status /= super_lord$) cycle
    length = 0
    do j = 1, lord0%n_slave
      slave => pointer_to_slave(lat, lord0, j)
      length = length + slave%value(l$)
    enddo
    if (abs(length - lord0%value(l$)) > dl_tol * abs(lord0%value(l$))) then
      call out_io (s_fatal$, r_name, &
              'INCONSISTANT SUPER_LORD/SUPER_SLAVE LENGTHS!', &
              'LORD: ' // lord0%name, &
              'LENGTH: \es16.9\ ', &
              'SUM OF SLAVE LENGTHS: \es16.9\ ', r_array = (/ lord0%value(l$), length /) )
      call err_exit
    endif
  enddo
endif

end subroutine

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

real(rp) s_start, s_start2, s_end, tot_len

character(40) :: r_name = 'adjust_super_lord_s_position'

!

if (lord%lord_status /= super_lord$) then
   call out_io (s_abort$, r_name, 'ELEMENT IS NOT A LORD! ' // lord%name)
   call err_exit
endif

! If a super lord is moved then we just need to adjust the start and end edges.
! Since we don't want to kill taylor maps due to round-off errors we only change
! the length if the percentage or absolute change is more than 10^-10

s_end = lord%s + lord%value(s_offset_tot$)
slave => pointer_to_slave (lat, lord, lord%n_slave)
s_start = slave%s - slave%value(l$)
if (abs(s_end - s_start - slave%value(l$)) > bmad_com%significant_longitudinal_length * &
                              (1 + abs(slave%value(l$)))) slave%value(l$) = s_end - s_start
slave%s = s_end

! Adjust start position

slave => pointer_to_slave (lat, lord, 1)
tot_len = lat%branch(slave%ix_branch)%param%total_length
s_start = s_end - lord%value(l$)
if (s_start < 0) s_start = s_start + tot_len
s_start2 = slave%s - slave%value(l$)
if (s_start < 0) s_start = s_start + tot_len
slave%value(l$) = slave%value(l$) + s_start2 - s_start

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_lord (lat, lord)
!
! Subroutine to calculate the attributes of group slave elements
!-

Subroutine makeup_group_lord (lat, lord)   

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: lord
type (ele_struct), pointer :: slave

real(rp) delta, coef
real(rp), pointer :: r_ptr

integer ix, iv, slave_stat, i

logical moved, err_flag

character(20) :: r_name = 'makeup_group_lord'

!

delta = lord%value(command$) - lord%value(old_command$)    ! change
lord%value(old_command$) = lord%value(command$) ! save old

moved = .false.   ! have we longitudinally moved an element?

do i = lord%ix1_slave, lord%ix2_slave

  ix = lat%control(i)%ix_slave
  iv = lat%control(i)%ix_attrib
  slave_stat = lat%ele(ix)%slave_status
  slave => lat%ele(ix)

  if (iv == l$) then
    moved = .true.
    if (slave_stat /= free$ .and. slave_stat /= super_slave$) then
      call out_io (s_abort$, r_name, "A GROUP: " // lord%name, &
                    "CONTROLS THE LENGTH OF A LORD ELEMENT: " // slave%name)
      call err_exit
    endif
  endif
  coef = lat%control(i)%coef
  call pointer_to_indexed_attribute (slave, iv, .false., r_ptr, err_flag)
  if (err_flag) call err_exit
  r_ptr = r_ptr + delta * coef
enddo

if (moved) call s_calc (lat)       ! recompute s distances

end subroutine

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

real(rp) s, slave_val(n_attrib_maxx)
real(rp) d, e, r_lord, r_slave, cos_lord, cos_e, sin_lord, sin_lorde
real(rp) ang_slave, ang_lord, ang_slave_old, d1, d2
real(rp) cos_e2, d_theta, ang_dlord, cos_lorde1, cos_dlord
real(rp) w0_mat(3,3), w1_mat(3,3), w1_inv_mat(3,3), offset(3), dw_mat(3,3)
real(rp) theta, phi, psi, w0_inv_mat(3,3)

integer i, j, ix_slave, ic, ix_s0, ix_patch_in_slave, n_pass
character(40) :: r_name = 'makeup_multipass_slave'

!

ix_slave = slave%ix_ele
branch => lat%branch(slave%ix_branch)
j =  lat%ic(slave%ic1_lord)
lord => lat%ele(lat%control(j)%ix_lord)
n_pass = j - lord%ix1_slave + 1  ! pass number for slave

slave_val = slave%value  ! save

slave%value = lord%value
if (lord%key == lcavity$ .or. lord%key == rfcavity$) then
  slave%value(dphi0$)        = slave_val(dphi0$)
  slave%value(E_tot_start$)  = slave_val(E_tot_start$)
  slave%value(p0c_start$)    = slave_val(p0c_start$)
endif

! A slave's field_master = T irregardless of the lord's setting.
! This is to make attribute_bookkeeper compute the correct normalized field strength.

slave%value(e_tot$) = slave_val(e_tot$)
slave%value(p0c$)   = slave_val(p0c$)
slave%value(n_ref_pass$)    = 0
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
  slave%a_pole = lord%a_pole
  slave%b_pole = lord%b_pole
endif

! wakes

if (associated (slave%r)) slave%r = lord%r
if (associated (slave%const)) slave%const = lord%const
if (associated (slave%wake)) then
  slave%wake%sr_table      = lord%wake%sr_table
  slave%wake%sr_mode_long  = lord%wake%sr_mode_long
  slave%wake%sr_mode_trans = lord%wake%sr_mode_trans
  slave%wake%lr            = lord%wake%lr
  do i = 1, size(lord%wake%lr)
    slave%wake%lr(i)%t_ref = lord%wake%lr(i)%t_ref - slave%ref_time
  enddo
endif

! methods

slave%mat6_calc_method = lord%mat6_calc_method
slave%tracking_method  = lord%tracking_method
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

    offset = (/ f0%x-f1%x, f0%y-f1%y, f0%z-f1%z /)
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

    patch_in_lord => pointer_to_lord(lat, lord, 1)
    patch_in_slave => pointer_to_slave (lat, patch_in_lord, n_pass)
    start%vec = 0
    do i = patch_in_slave%ix_ele, ix_slave - 1
      call track1 (start, branch%ele(i), lat%param, end)
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

  select case (lord%ref_orbit)

  case (single_ref$)
    slave%value(b_field$)     = lord%value(b_field$) * slave%value(p0c$) / lord%value(p0c$) 
    slave%value(b_field_err$) = lord%value(b_field$) + lord%value(b_field_err$) - &
                                                                   slave%value(b_field$)
    slave%value(g_err$) = (lord%value(g$) + lord%value(g_err$)) * &
                                    lord%value(p0c$) / slave%value(p0c$) - lord%value(g$)


  case (match_global_coords$)

    slave%value(g$) = lord%value(g$) * lord%value(p0c$) / slave%value(p0c$)

    ! e1 and e2 and l are determined by the reference orbit of this pass with respect
    ! to the reference orbit of the reference pass. 
    ! Assumption: the slave element lies in the (x, z) plane.

    if (slave%floor%phi /= 0 .or. slave%floor%psi /= 0) then
       call out_io (s_fatal$, r_name, 'MULTIPASS ELEMENT: ' // lord%name, &
                     'WHICH HAS REF_ORBIT = MATCH_GLOBAL_COORDS DOES NOT LIE IN THE (X, Z) PLANE!')
      call err_exit
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
      ang_slave = asin((r_lord * (sin(ang_dlord) - sin(d_theta)) - &
                        d1 * sin(d_theta + lord%value(e1$)) + &
                        d2 * sin(ang_dlord - lord%value(e2$))) / r_slave)
      if (abs(ang_slave - ang_slave_old) < 1e-6 * abs(ang_slave)) exit
      ang_slave_old = ang_slave
    enddo

    slave%value(angle$) = ang_slave
    slave%value(l$) = ang_slave * r_slave
    slave%value(rho$) = r_slave
    slave%value(e1$) = lord%value(e1$) + d_theta
    slave%value(e2$) = lord%value(e2$) + ang_slave - ang_lord - d_theta


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
        ang_slave = asin((r_lord * sin_lord + d * sin_lorde) / r_slave)
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

    if (lord%value(k1$) /= 0 .or. lord%value(k2$) /= 0 .or. associated(lord%a_pole)) then
      call out_io (s_fatal$, r_name, &
            'MULTIPASS BEND ELEMENT: ' // lord%name, &
            'WITH THE REF_ORBIT ATTRIBUTE SET TO: ' // ref_orbit_name(lord%ref_orbit), &
            'HAS A NONZERO HIGHER ORDER MULTIPOLE!', &
            'THIS IS NOT ALLOWED. SEE THE BMAD MANUAL FOR MORE DETAILS.')
      call err_exit
    endif

  case default
    call out_io (s_fatal$, r_name, 'BAD REF_ORBIT VALUE: \i0\ ', &
                           'FOR: ' // lord%name, i_array = (/ lord%ref_orbit /) )
    call err_exit

  end select
endif

end subroutine

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
type (ele_struct), pointer :: lord, slave0
type (ele_struct), save :: sol_quad
type (branch_struct), pointer :: branch

integer i, j, ix_con, ix, ix_slave, ix_lord

real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3, c, s
real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
real(rp) a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4), leng, offset
real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

logical, save :: init_needed = .true.

character(20) :: r_name = 'makeup_super_slave'

! init

if (init_needed) then
  call init_ele (sol_quad)
  sol_quad%key = sol_quad$
  init_needed = .false.
endif

! Super_slave:

if (slave%slave_status /= super_slave$) then
   call out_io(s_abort$, r_name, "ELEMENT IS NOT AN SUPER SLAVE: " // slave%name)
  call err_exit
endif

branch => lat%branch(slave%ix_branch)
ix_slave = slave%ix_ele

! If this slave is the last slave for some lord (so that then longitudinal 
! end of the slave matches the end of the lord) then the limits of the lord
! are transfered to the slave.

!-----------------------------------------------------------------------
! 1 super_lord for this super_slave: just transfer attributes except length

if (slave%n_lord == 1) then

  ix_con = lat%ic(slave%ic1_lord)  
  ix = lat%control(ix_con)%ix_lord
  lord => lat%ele(ix)

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (ix_con /= lord%ix1_slave) then
    slave%map_ref_orb_in = branch%ele(ix_slave-1)%map_ref_orb_out
  endif

  ! Find the offset from the longitudinal start of the lord to the start of the slave

  offset = 0 ! length of all slaves before this one
  do i = lord%ix1_slave, ix_con-1
    j = lat%control(i)%ix_slave
    offset = offset + branch%ele(j)%value(l$)
  enddo

  ! If this is the last slave, adjust it's length to be consistant with
  ! The lord length. Then do the rest of the bookkeeping

  if (ix_con == lord%ix2_slave) slave%value(l$) = lord%value(l$) - offset

  call makeup_super_slave1 (slave, lord, offset, lat%param, &
                             (ix_con == lord%ix1_slave), (ix_con == lord%ix2_slave))

  return

endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.
                                         
k_x = 0
k_y = 0
x_kick = 0
y_kick = 0
a_tot = 0
b_tot = 0
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
value(E_tot$) = slave%value(E_tot$)
value(p0c$) = slave%value(p0c$)

s_slave = slave%s - value(l$)/2  ! center of slave
slave%is_on = .false.

! sum over all lords...

do j = 1, slave%n_lord

  lord => pointer_to_lord (lat, slave, j, ix_con)

  ! Physically, the lord length cannot be less than the slave length.
  ! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

  if (abs(slave%value(l$)) > abs(lord%value(l$))) then
    coef = 1
  else
    coef = slave%value(l$) / lord%value(l$) 
  endif

  if (lord%lord_status /= super_lord$) then
    call out_io (s_abort$, r_name, &
          "SUPER_SLAVE HAS A CONTROL ELEMENT THAT IS NOT A SUPER_LORD", &
          'SLAVE: ' //  slave%name // '  \i\ ', &
          'LORD:  ' //  lord%name  // '  \i\ ', i_array = (/ ix_slave, ix /) )
    call err_exit
  endif

  if (associated(lord%wake)) then
    call out_io (s_abort$, r_name, &
            'SUPERPOSITION OF MULTIPLE ELEMENTS WITH WAKES NOT YET IMPLEMENTED!', &
            'SUPER_LORD: ' // lord%name)
    call err_exit
  endif

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (ix_con /= lord%ix1_slave) then
    slave%map_ref_orb_in = branch%ele(ix_slave-1)%map_ref_orb_out
  endif

  ! Choose the smallest ds_step of all the lords.

  if (value(ds_step$) == 0 .or. lord%value(ds_step$) < value(ds_step$)) &
                                        value(ds_step$) = lord%value(ds_step$)

  ! Coupler and aperture calc.

  call compute_slave_aperture (value, slave, lord, &
                               ix_con == lord%ix1_slave, ix_con == lord%ix2_slave)

  if (slave%key == lcavity$) call compute_slave_coupler (value, slave, lord, &
                               ix_con == lord%ix1_slave, ix_con == lord%ix2_slave)

  ! Methods

  if (j == slave%ic1_lord) then
    slave%mat6_calc_method = lord%mat6_calc_method
    slave%tracking_method  = lord%tracking_method
  else
    if (slave%mat6_calc_method /= lord%mat6_calc_method) then
      ix = lat%control(lat%ic(slave%ic1_lord))%ix_lord
      call out_io(s_abort$, r_name, 'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT', &
           'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(branch%ele(ix)%name))
      call err_exit
    endif
    if (slave%tracking_method /= lord%tracking_method) then
      ix = lat%control(lat%ic(slave%ic1_lord))%ix_lord
      call out_io(s_abort$, r_name, ' TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT', &
           'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(branch%ele(ix)%name))
      call err_exit
    endif
  endif

  ! kicks, etc.

  if (.not. lord%is_on) cycle
  slave%is_on = .true.  ! Slave is on if at least one lord is on

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

  if (associated(lord%a_pole)) then
    call multipole_ele_to_kt (lord, +1, knl, t, .true.)
    call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
    a_tot = a_tot + a
    b_tot = b_tot + b
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
    s_del = modulo2 (s_del, lat%param%total_length/2)

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
    call out_io (s_abort$, r_name, &
               'CODING NOT YET IMPLEMENTED FOR: ' // key_name(slave%key))
    call err_exit

  ! hkicker, vkicker, kicker, etc. looks like drifts

  case (hkicker$, vkicker$, kicker$, instrument$, monitor$, pipe$, rcollimator$)

  ! default

  case default
    call out_io (s_abort$, r_name, &
               'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
    call err_exit

  end select

enddo

!------------------------------
! stuff sums into slave element

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

if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
  call multipole_init(slave)
  call multipole_ab_to_kt(a_tot, b_tot, knl, t)
  call multipole_kt_to_ab(knl*slave%value(l$), t-tilt, a, b)
  slave%a_pole = a
  slave%b_pole = b
  slave%value(radius$) = 1
elseif (associated(slave%a_pole)) then
  deallocate (slave%a_pole, slave%b_pole)
endif

!-----------------------------

select case (slave%key)

case (sextupole$) 

  if (k_x == 0 .and. k_y == 0) return

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

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

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

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

  if (ks /= 0) then
    x_o_sol = ks_xo_sum / ks
    x_p_sol = ks_xp_sum / ks
    y_o_sol = ks_yo_sum / ks
    y_p_sol = ks_yp_sum / ks
  endif

  if (k_x == 0 .and. k_y == 0) then  ! pure solenoid
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    deallocate (slave%a_pole, slave%b_pole, stat = ix)
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

    if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) return

    t_2 = (/ x_o, x_p, y_o, y_p /)
    call tilt_coords (tilt, t_2, .true.)  ! set

    l_slave = slave%value(l$)

    t_1 = (/ t_2(2), 0.0_rp, t_2(4), 0.0_rp /)
    t_2(1) = t_2(1) + ks * t_2(4) / k1 
    t_2(3) = t_2(3) + ks * t_2(2) / k1
             
    call mat_make_unit (T_end)
    T_end(4,1) =  ks / 2
    T_end(2,3) = -ks / 2

    sol_quad%value(ks$) = ks
    sol_quad%value(k1$) = k1
    sol_quad%value(l$)  = l_slave
    call make_mat6 (sol_quad, lat%param)
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

    call tilt_coords (tilt, beta, .false.)  ! unset

    slave%value(x_offset$) = beta(1) + x_o_sol
    slave%value(x_pitch$)  = beta(2) + x_p_sol
    slave%value(y_offset$) = beta(3) + y_o_sol
    slave%value(y_pitch$)  = beta(4) + y_p_sol
  endif

! bend_sol_quad

case (bend_sol_quad$)
  call out_io (s_abort$, r_name, &
               'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
  call err_exit

end select

! If the slave has %field_master = T then we need to convert k1, etc values to field quantities.

if (slave%field_master) then
  slave%field_master = .false.   ! So attribute_bookkeeper will do the right thing.
  call attribute_bookkeeper (slave, lat%param)
  slave%field_master = .true.
endif


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave1 (slave, lord, offset, param, at_entrance_end, at_exit_end)
!
! Routine to transfer the %value, %wig_term, and %wake%lr information from a 
! superposition lord to a slave when the slave has only one lord.
!
! Note: attribute_bookkeeper needs to be called after calling this routine.
!
! Modules needed:
!   use bmad
!
! Input:
!   slave  -- Ele_struct: Slave element.
!     %value(l$) -- Length of slave.
!   lord   -- Ele_struct: Lord element.
!   offset -- Real(rp): offset of entrance end of slave from entrance end of the lord.
!   param  -- Lat_param_struct: lattice paramters.
!   at_entrance_end -- Logical: Slave contains the lord's entrance end?
!   at_exit_end     -- Logical: Slave contains the lord's exit end?
!
! Output:
!   slave -- Ele_struct: Slave element with appropriate values set.
!-

subroutine makeup_super_slave1 (slave, lord, offset, param, at_entrance_end, at_exit_end)

implicit none

type (ele_struct), target :: slave, lord
type (lat_param_struct) param

real(rp) offset, s_del, coef, r, e_tot_start
real(rp) value(n_attrib_maxx)
integer i
logical at_entrance_end, at_exit_end
character(24) :: r_name = 'makeup_super_slave1'

! Physically, the lord length cannot be less than the slave length.
! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

if (lord%value(l$) == 0) then
  call out_io (s_fatal$, r_name, 'LORD HAS ZERO LENGTH!')
  call err_exit
endif

if (abs(slave%value(l$)) > abs(lord%value(l$))) then
  coef = 1
else
  coef = slave%value(l$) / lord%value(l$) 
endif

value = lord%value
value(l$) = slave%value(l$)                 ! do not change slave length

if (lord%key == wiggler$) then
  value(z_patch$) = slave%value(z_patch$)
endif

if (lord%key == hkicker$ .or. lord%key == vkicker$) then
  value(kick$) = lord%value(kick$) * coef
else
  value(hkick$) = lord%value(hkick$) * coef
  value(vkick$) = lord%value(vkick$) * coef
endif

if (slave%key == rfcavity$) value(voltage$) = lord%value(voltage$) * coef

slave%aperture_at = no_end$
call compute_slave_aperture (value, slave, lord, at_entrance_end, at_exit_end)

if (slave%key == lcavity$) then
  slave%value(coupler_at$) = no_end$
  call compute_slave_coupler (value, slave, lord, at_entrance_end, at_exit_end)
endif

! s_del is the distance between lord and slave centers

s_del = offset + slave%value(l$)/2 - lord%value(l$)/2
value(x_pitch$)  = value(x_pitch_tot$)
value(y_pitch$)  = value(y_pitch_tot$)
value(x_offset$) = value(x_offset_tot$) + s_del * value(x_pitch_tot$)
value(y_offset$) = value(y_offset_tot$) + s_del * value(y_pitch_tot$)
value(tilt$)     = value(tilt_tot$)

slave%value = value
slave%is_on = lord%is_on
slave%mat6_calc_method = lord%mat6_calc_method
slave%tracking_method  = lord%tracking_method

! If a wiggler: 
! must keep track of where we are in terms of the unsplit wiggler.
! This is for anything which does not try to make a homogeneous approximation.
! l_original is the length of the unsplit original wiggler.
! l_start is the starting point with respect to the original wiggler.
! l_end is the ending point with respect to the original wiggler.

if (slave%key == wiggler$) then
  slave%value(n_pole$) = lord%value(n_pole$) * coef
  slave%value(l_original$) = lord%value(l$)

  slave%value(l_start$) = offset
  slave%value(l_end$)   = slave%value(l_start$) + slave%value(l$)

  if (associated(lord%wig_term)) then
    if (.not. associated (slave%wig_term) .or. &
            size(slave%wig_term) /= size(lord%wig_term)) then
      if (associated (slave%wig_term)) deallocate (slave%wig_term)
      allocate (slave%wig_term(size(lord%wig_term)))
    endif
    do i = 1, size(lord%wig_term)
      slave%wig_term(i) = lord%wig_term(i)
      slave%wig_term(i)%phi_z = lord%wig_term(i)%phi_z + &
                             lord%wig_term(i)%kz * slave%value(l_start$)
    enddo
  else
    if (associated (slave%wig_term)) deallocate (slave%wig_term)
  endif

endif

! If a custom element: 
! Must keep track of where we are in terms of the unsplit element.
! See wiggler above for more details.

if (slave%key == custom$) then
  slave%value(l_original$) = lord%value(l$)
  slave%value(l_start$)    = offset
  slave%value(l_end$)      = slave%value(l_start$) + slave%value(l$)
endif

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

if (associated (slave%wake)) then
  slave%wake%lr%freq_in   = lord%wake%lr%freq_in
  slave%wake%lr%freq      = lord%wake%lr%freq
  slave%wake%lr%Q         = lord%wake%lr%Q
  slave%wake%lr%angle     = lord%wake%lr%angle
  slave%wake%lr%m         = lord%wake%lr%m
  slave%wake%lr%polarized = lord%wake%lr%polarized
  slave%wake%lr%r_over_q  = lord%wake%lr%r_over_q * coef
endif

! lcavity energy bookkeeping

if (slave%key == lcavity$) then
  slave%value(e_loss$) = lord%value(e_loss$) * coef
  if (any(slave%value /= slave%old_value)) then ! Only do this if necessary
    r = offset / lord%value(l$)
    e_tot_start = (1-r) * lord%value(e_tot_start$) + r * lord%value(e_tot$) 
    call attribute_bookkeeper (slave, param)
  endif
endif

end subroutine

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
real(rp) value(n_attrib_maxx)
logical at_entrance_end, at_exit_end

!

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
end select

if (slave%aperture_at == no_end$) then
  value(x1_limit$) = 0
  value(x2_limit$) = 0
  value(y1_limit$) = 0
  value(y2_limit$) = 0
endif

end subroutine

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
real(rp) value(n_attrib_maxx)
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

end subroutine

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

real(rp) value(n_attrib_special_maxx), coef, ds
real(rp), pointer :: r_ptr
integer i, ix_con, ix, iv, ix_slave, icom, l_stat
logical used(n_attrib_special_maxx), multipole_set, err_flag

character(40) :: r_name = 'makeup_overlay_and_girder_slave'

!
                             
branch => lat%branch(slave%ix_branch)
l_stat = slave%lord_status
ix_slave = slave%ix_ele

value = 0
used = .false.
multipole_set = .false.
slave%on_a_girder = .false.

do i = 1, slave%n_lord
  lord => pointer_to_lord (lat, slave, i, ix_con)

  if (lord%lord_status == multipass_lord$) cycle
  if (lord%lord_status == group_lord$) cycle

  if (lord%lord_status == girder_lord$) then
    ds = (slave%s - slave%value(l$)/2) - lord%value(s_center$) 
    slave%value(x_offset_tot$) = slave%value(x_offset$) + &
                   ds * lord%value(x_pitch$) + lord%value(x_offset$)
    slave%value(y_offset_tot$) = slave%value(y_offset$) + &
                   ds * lord%value(y_pitch$) + lord%value(y_offset$)
    slave%value(s_offset_tot$) = slave%value(s_offset$) + lord%value(s_offset$)
    slave%value(x_pitch_tot$)  = slave%value(x_pitch$)  + lord%value(x_pitch$)
    slave%value(y_pitch_tot$)  = slave%value(y_pitch$)  + lord%value(y_pitch$)
    slave%value(tilt_tot$)     = slave%value(tilt$)     + lord%value(tilt$)
    slave%on_a_girder = .true.
    cycle
  endif

  if (lord%lord_status /= overlay_lord$) then
    call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY_LORD \i\ ', ix_slave)
    call type_ele (slave, .true., 0, .false., 0, .true., lat)
    call err_exit
  endif     

  coef = lat%control(ix_con)%coef
  iv = lat%control(ix_con)%ix_attrib
  call pointer_to_indexed_attribute (lord, lord%ix_value, .false., r_ptr, err_flag)
  if (err_flag) call err_exit
  value(iv) = value(iv) + r_ptr * coef
  used(iv) = .true.
  if (iv > n_attrib_maxx) multipole_set = .true.
enddo

where (used(1:n_attrib_maxx)) slave%value = value(1:n_attrib_maxx)
if (multipole_set) then
  where (used(a0$:a20$)) slave%a_pole = value(a0$:a20$)
  where (used(b0$:b20$)) slave%b_pole = value(b0$:b20$)
endif

! If no girder then simply transfer tilt to tilt_tot, etc.

if (.not. slave%on_a_girder) then
  slave%value(tilt_tot$)     = slave%value(tilt$)
  slave%value(x_offset_tot$) = slave%value(x_offset$)
  slave%value(y_offset_tot$) = slave%value(y_offset$)
  slave%value(s_offset_tot$) = slave%value(s_offset$)
  slave%value(x_pitch_tot$)  = slave%value(x_pitch$)
  slave%value(y_pitch_tot$)  = slave%value(y_pitch$)
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
! Note: This routine does not do any other bookkeeping. Consider using
! control_bookkeeper or lattice_bookkeeper instead.
! 
! BEAMBEAM:   
!     bbi_const$ = param%n_part * mass_of(param%particle) * charge$ * r_e /
!                           (2 * pi * p0c$ * (sig_x$ + sig_y$)
!
! ELSEPARATOR:
!     e_field$ = sqrt(hkick$**2 + vkick$**2) * p0c$ / l$
!     voltage$ = e_field$ * gap$ 
!
! LCAVITY:    
!     delta_e$ = gradient$ * L$ 
!     E_tot$   = E_tot_start$ + gradient$ * l$ * cos(phase)
!     p0c$     = sqrt(E_tot$**2 - mc2^2)
! 
! RFCAVITY:   
!     rf_frequency$ = harmon$ * c_light / param%total_length (only if harmon$ /= 0)
!
! SBEND:      
!     angle$   = L$ * G$
!     l_chord$ = 2 * sin(Angle$/2) / G$
!     rho$     = 1 / G$
!
! WIGGLER:    
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!     n_pole$ = L$ / l_pole$
!     z_patch$
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- Ele_struct: Element with attributes 
!   param      -- lat_param_struct: 
!
! Output:
!   ele            -- Ele_struct: Element with self-consistant attributes.
!     %map_ref_orb_out -- Reference orbit to be used for the next
!                         super_slave wiggler z_patch calculation. 
!                         This is to be only used by the control_bookkeeper routine.
!
! Programming Note: If the dependent attributes are changed then 
!       the attribute_free routine must be modified.
!-

subroutine attribute_bookkeeper (ele, param)

use symp_lie_mod, only: symp_lie_bmad

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) start, end

real(rp) factor, f2, phase, E_tot, dval(n_attrib_maxx)
real(rp), pointer :: val(:)

character(20) ::  r_name = 'attribute_bookkeeper'

logical non_offset_changed, offset_changed, offset_nonzero, z_patch_calc_needed
logical v_mask(n_attrib_maxx), offset_mask(n_attrib_maxx)
logical :: init_needed = .true.
logical :: debug = .false.  ! For debugging purposes

! If no change then we don't need to do anything

if (ele%key == taylor$) return

val => ele%value
val(check_sum$) = 0
if (associated(ele%a_pole)) val(check_sum$) = sum(ele%a_pole) + sum(ele%b_pole)
z_patch_calc_needed = (ele%key == wiggler$ .and. val(z_patch$) == 0 .and. val(p0c$) /= 0)

if (all(val == ele%old_value) .and. .not. z_patch_calc_needed) return
if (debug) dval = val - ele%old_value

! Transfer tilt to tilt_tot, etc.

if (.not. ele%on_a_girder .and. ele%key /= match$) then
  val(tilt_tot$)     = val(tilt$)
  val(x_offset_tot$) = val(x_offset$)
  val(y_offset_tot$) = val(y_offset$)
  val(s_offset_tot$) = val(s_offset$)
  val(x_pitch_tot$)  = val(x_pitch$)
  val(y_pitch_tot$)  = val(y_pitch$)
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

  val(hkick$) = factor * val(BL_hkick$)
  val(vkick$) = factor * val(BL_vkick$)

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

  val(BL_hkick$) = factor * val(hkick$)
  val(BL_vkick$) = factor * val(vkick$)

endif

! Dependent attribute bookkeeping.
! Note: If the dependent attributes are changed then attribute_free 
!       must be modified.

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

! Lcavity
! Only do the calculation if the starting energy is not zero since 
! attribute_bookkeeper can be called before the attributes are set.

case (lcavity$)
  if (val(E_tot_start$) /= 0) then
    val(delta_e$) = val(gradient$) * val(L$) 
    phase = twopi * (val(phi0$) + val(dphi0$)) 
    E_tot = val(E_tot_start$) + val(gradient$) * val(l$) * cos(phase)
    E_tot = E_tot - val(e_loss$) * param%n_part * e_charge
    if (e_tot /= val(e_tot$)) then ! Only do this if necessary
      val(E_tot$) = E_tot
      call convert_total_energy_to (E_tot, param%particle, pc = val(p0c$))
      call convert_total_energy_to (val(e_tot_start$), param%particle, pc = val(p0c_start$))
    endif
  endif

! RFcavity

case (rfcavity$)
  if (val(harmon$) /= 0) val(rf_frequency$) =  val(harmon$) * c_light / param%total_length 

! BeamBeam

case (beambeam$)

  if (val(n_slice$) == 0) val(n_slice$) = 1.0 ! revert to default

  if (val(charge$) == 0 .or. param%n_part == 0) then
    val(bbi_const$) = 0

  else

    if (val(sig_x$) == 0 .or. val(sig_y$) == 0) then
      call out_io(s_abort$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!')
      call type_ele(ele, .true., 0, .false., 0, .false.)
      call err_exit
    endif

    val(bbi_const$) = -param%n_part * mass_of(param%particle) * val(charge$) * r_e /  &
                             (2 * pi * val(p0c$) * (val(sig_x$) + val(sig_y$)))

  endif

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

  ! Periodic_type wigglers have a single %wig_term for use with tracking, etc.
  ! The phase of this term is set so that tracking with a particle starting
  ! on-axis ends on-axis. For this to be true, there must be an integer number
  ! of poles.

  if (ele%sub_key == periodic_type$) then
    if (.not. associated(ele%wig_term)) allocate (ele%wig_term(1))

    if (val(l_pole$) == 0) then
      ele%wig_term(1)%ky = 0
    else
      ! Use an integer number of poles in calculating ky.
      ele%wig_term(1)%ky = pi * val(n_pole$) / val(l$)
    endif
    ele%wig_term(1)%coef   = val(b_max$)
    ele%wig_term(1)%kx     = 0
    ele%wig_term(1)%kz     = ele%wig_term(1)%ky
    ele%wig_term(1)%phi_z  = -ele%wig_term(1)%ky * val(l$) / 2 
    ele%wig_term(1)%type   = hyper_y$
  endif

end select

! num_steps

if (val(ds_step$) /= 0) ele%num_steps = abs(nint(val(l$) / val(ds_step$)))
if (val(ds_step$) == 0 .or. ele%num_steps == 0) ele%num_steps = 1

! If things have changed we need to kill the Taylor Map and gen_field.
! The old_value array tells us this.

if (all(val == ele%old_value) .and. .not. z_patch_calc_needed) return

if (init_needed) then
  v_mask = .true.
  v_mask( (/ x_offset$, y_offset$, s_offset$, tilt$, x_pitch$, &
            y_pitch$, x_offset_tot$, y_offset_tot$, s_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$, z_patch$ /) ) = .false.
  offset_mask = .not. v_mask
  offset_mask(z_patch$) = .false.
  v_mask( (/ x1_limit$, x2_limit$, y1_limit$, y2_limit$ /) ) = .false.
  init_needed = .false.
endif

non_offset_changed = (any(val /= ele%old_value .and. v_mask))
offset_changed =  (any(val /= ele%old_value .and. offset_mask))
offset_nonzero = (any(val /= 0 .and. offset_mask))

! If an element has just been offset and bmad_com%conserve_taylor_map = T then 
! conserve the taylor map.

if (associated(ele%taylor(1)%term) .and. ele%map_with_offsets .and. &
        offset_nonzero .and. bmad_com%conserve_taylor_maps .and. &
        .not. non_offset_changed .and. ele%key /= patch$) then
  ele%map_with_offsets = .false.
  call out_io (s_info$, r_name, &
      'Note: bmad_com%conserve_taylor_maps = True (this is the default)', &
      'But: Element has just been offset: ' // ele%name, &
      "To conserve the element's Taylor map, I will set ele%map_with_offsets = False.")
endif

! Kill the taylor map and gen_field if necessary.

if (non_offset_changed .or. (offset_changed .and. ele%map_with_offsets)) then
  if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)
  if (associated(ele%gen_field)) call kill_gen_field(ele%gen_field)
  if (ele%key == wiggler$) then
    val(z_patch$) = 0
    z_patch_calc_needed = (ele%key == wiggler$ .and. val(p0c$) /= 0)
  endif
endif

! compute the z_patch for a wiggler if needed.
! This is normally zero except for split wiggler sections.

if (z_patch_calc_needed) then
  start = ele%map_ref_orb_in
  call symp_lie_bmad (ele, param, start, end, .false., offset_ele = .false.)
  val(z_patch$) = end%vec(5)
  end%vec(5) = 0
  if (val(z_patch$) == 0) val(z_patch$) = 1e-30 ! something non-zero.
  ele%map_ref_orb_out = end             ! save for next super_slave
endif

ele%old_value = val

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine changed_attribute_bookkeeper (lat, ele, a_ptr)
!
! Subroutine to do bookkeeping when a particular attribute has been altered.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: Lattice with the changed attribute.
!   ele    -- ele_struct: Element being modified.
!   a_ptr  -- Real(rp), pointer: Pointer to the changed attribute.
!
! Output:
!   lat  -- lat_struct: Lattice with appropriate changes.
!-

subroutine changed_attribute_bookkeeper (lat, ele, a_ptr)

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: ele

real(rp), pointer :: a_ptr
real(rp) v_mat(4,4), v_inv_mat(4,4), eta_vec(4), eta_xy_vec(4)

logical coupling_change

!

if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)

if (ele%key == init_ele$) then
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
    ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + &
                                                ele%c_mat(1,2)*ele%c_mat(2,1))
    coupling_change = .true.
  endif

  if (associated(a_ptr, ele%x%eta) .or. associated(a_ptr, ele%x%etap) .or. &
      associated(a_ptr, ele%y%eta) .or. associated(a_ptr, ele%y%etap) .or. &
      coupling_change) then 
    call make_v_mats (ele, v_mat, v_inv_mat)
    eta_xy_vec = (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)
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
    eta_vec = (/ ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap /)
    eta_xy_vec = matmul (v_mat, eta_vec)
    ele%x%eta  = eta_xy_vec(1)
    ele%x%etap = eta_xy_vec(2)
    ele%y%eta  = eta_xy_vec(3)
    ele%y%etap = eta_xy_vec(4)
    return
  endif

endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_lat_taylors (lat_in, lat_out, type_out, transfered_all)
!
! Subroutine to transfer the taylor maps from the elements of one lattice to
! the elements of another. The elements are matched between the lattices so 
! that the appropriate element in lattice_out will get the correct Taylor map
! even if the order of the elements is different in the 2 lattices.
!
! Note: The transfered Taylor map will be truncated to bmad_com%taylor_order.
! Note: If the taylor_order of an element in lattice_in is less than 
!   bmad_com%taylor_order then it will not be used.  
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in    -- lat_struct: Input lattice with Taylor maps.
!   type_out  -- Logical: If True then print a message for each Taylor map
!                 transfered.
!
! Output:
!   lat_out   -- lat_struct: lattice to receive the Taylor maps.
!   transfered_all 
!             -- Logical, optional: Set True if a Taylor map is found
!                 for all elements in lattice_out that need one. False otherwise.
!-

subroutine transfer_lat_taylors (lat_in, lat_out, type_out, transfered_all)

implicit none

type (lat_struct), target, intent(in) :: lat_in
type (lat_struct), target, intent(inout) :: lat_out
type (ele_struct), pointer :: ele_in, ele_out

integer i, j
integer n_in, ix_in(ubound(lat_in%ele, 1))
 
logical, intent(in)  :: type_out
logical, optional :: transfered_all

character(25) :: r_name = 'transfer_lat_taylors'

! check global parameters

if (present(transfered_all)) transfered_all = .true.

if (lat_in%ele(0)%value(E_tot$) /= lat_out%ele(0)%value(E_tot$)) then
  if (type_out) then
     call out_io (s_warn$, r_name, &
            'THE LATTICE ENERGIES ARE DIFFERENT. TAYLOR MAPS NOT TRANSFERED.')
  endif
  if (present(transfered_all)) transfered_all = .false.
  return
endif

! Find the taylor series in the first lattice.

n_in = 0
do i = 1, lat_in%n_ele_max
  if (associated(lat_in%ele(i)%taylor(1)%term)) then
    if (bmad_com%taylor_order > lat_in%ele(i)%taylor_order) cycle
    n_in = n_in + 1
    ix_in(n_in) = i
  endif
enddo

! Go through lattice_out and match elements.
! If we have a match transfer the Taylor map.
! Call attribute_bookkeeper before transfering the taylor map to make sure
! the ele_out%value(:) arry is correct. 

out_loop: do i = 1, lat_out%n_ele_max

  ele_out => lat_out%ele(i)

  do j = 1, n_in

    ele_in => lat_in%ele(ix_in(j))

    if (equivalent_taylor_attributes (ele_in, ele_out)) then
      if (type_out) call out_io (s_info$, r_name, &
          ' Reusing Taylor from: ' // trim(ele_in%name) // '  to: ' //  ele_out%name)
      call attribute_bookkeeper (ele_out, lat_out%param)
      call transfer_ele_taylor (ele_in, ele_out, bmad_com%taylor_order)
      cycle out_loop
    endif

  enddo

  if (ele_out%tracking_method == taylor$ .or. &
                  ele_out%mat6_calc_method == taylor$ .and. type_out) then
    call out_io (s_warn$, r_name, ' NO TAYLOR FOR: ' // ele_out%name)
    if (present(transfered_all)) transfered_all = .false.
  endif

enddo out_loop

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_on_off (key, lat, switch, orb, use_ref_orb)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a lattice. An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6.
!
! Modules needed:
!   use bmad
!
! Input:
!   key          -- Integer: Key name of elements to be turned on or off.
!                      [Key = quadrupole$, etc.]
!   lat             -- lat_struct: lattice structure holding the elements
!   switch       -- Integer: 
!                     on$            => Turn elements on.  
!                     off$           => Turn elements off. 
!                     save_state$    => Save present on/off state. 
!                                         No turning on or off is done.
!                     restore_state$ => Restore saved on/off state.
!   orb(0:)     -- Coord_struct, optional: Needed for lat_make_mat6
!   use_ref_orb -- Logical, optional: If present and true then use the
!                    present ele%ref_orb. Default is false.
!
! Output:
!   lat -- lat_struct: Modified lattice.
!-

subroutine set_on_off (key, lat, switch, orb, use_ref_orb)

implicit none

type (lat_struct) lat
type (coord_struct), optional :: orb(0:)
type (coord_struct) ref_orb

integer i, key               
integer, intent(in) :: switch

logical, optional :: use_ref_orb
logical old_state

character(20) :: r_name = 'set_on_off'

!

do i = 1, lat%n_ele_max

  if (lat%ele(i)%key /= key) cycle

  old_state = lat%ele(i)%is_on

  select case (switch)
  case (on$) 
    lat%ele(i)%is_on = .true.
  case (off$)
    lat%ele(i)%is_on = .false.
  case (save_state$)
    lat%ele(i)%old_is_on = lat%ele(i)%is_on
    cycle
  case (restore_state$)
    lat%ele(i)%is_on = lat%ele(i)%old_is_on
  case default
    call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
    call err_exit
  end select

  if (old_state .neqv. lat%ele(i)%is_on) then
    if (logic_option (.false., use_ref_orb)) then
      ref_orb = lat%ele(i)%map_ref_orb_in
      call make_mat6(lat%ele(i), lat%param, ref_orb)
    else
      call lat_make_mat6(lat, i, orb)
    endif
  endif

enddo

end subroutine

end module
