!+
! Subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset,
!                          param, include_upstream_end, include_downstream_end, err_flag, old_slice, orb_in)
!
! Routine to create an element that represents a longitudinal slice of the original element.
!
! Note: To save tracking computation time, if ele_in has taylor, or symp_lie_ptc
! for tracking_method or mat6_calc_method, then this will be changed to symp_lie_bmad 
! for wigglers and bmad_standard for everything else.
!
! Input:
!   ele_in            -- Ele_struct: Original element to slice
!   l_slice           -- Real(rp): Length of the slice
!   offset            -- Real(rp): Offset of entrance end of sliced_ele from entrance end of ele_in.
!   param             -- Lat_param_struct: lattice paramters.
!   include_upstream_end   -- Logical: Sliced_ele contains the ele's entrance end?
!   include_downstream_end -- Logical: Sliced_ele contains the ele's exit end?
!   old_slice         -- Logical, optional: Previous slice. If present this saves computation
!                          time of the reference energy and time at the start of the present slice.
!                          Also makes the ref energy continuous (there can be some small differences when
!                          using, say, runge_kutta tracking due to tracking tolerances).
!   orb_in            -- coord_struct, optional: Incoming orbit if calling routine is doing tracking through the slice.
!                         This is used when old_slice is not present and there may be an adjustment needed to
!                         the orbit ref energy (EG space charge tracking does not keep track of ref energy through an lcavity).
!
! Output:
!   sliced_ele -- Ele_struct: Sliced_ele element with appropriate values set.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!-

recursive subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset, &
                         param, include_upstream_end, include_downstream_end, err_flag, old_slice, orb_in)

use bookkeeper_mod, dummy_except => create_element_slice

implicit none

type (ele_struct), target :: sliced_ele, ele_in, ele0
type (ele_struct), optional :: old_slice
type (coord_struct), optional :: orb_in
type (ele_struct) :: ele2
type (lat_param_struct) param
type (coord_struct) time_ref_orb_out

real(rp) l_slice, offset, in_len, r, p0c_set
integer i
logical include_upstream_end, include_downstream_end, err_flag, err2_flag, rad_map_stale

character(*), parameter :: r_name = 'create_element_slice'

! Init

err_flag = .true.
in_len = ele_in%value(l$)

! Save values from old_slice if present in case the old_slice actual arg is same as sliced_ele.

rad_map_stale = .true.
ele0%value(l$) = 0
if (present(old_slice)) then
  ele0%value(p0c$)   = old_slice%value(p0c$)
  ele0%value(e_tot$) = old_slice%value(e_tot$)
  ele0%ref_time      = old_slice%ref_time
  ele0%value(l$)     = old_slice%value(l$)
  time_ref_orb_out   = old_slice%time_ref_orb_out
  if (associated(old_slice%rad_map)) rad_map_stale = old_slice%rad_map%stale
endif

!

if (.not. associated(sliced_ele%lord, ele_in) .or. sliced_ele%ix_ele /= ix_slice_slave$) then
  call transfer_ele(ele_in, sliced_ele, .true.)
endif

! A sad_mult with zero length is treated differently since edge fields are ignored
! Therefore, make sure a sad_mult has a finite length.

if (ele_in%key == sad_mult$ .and. l_slice == 0) then
  sliced_ele%value(l$) = sign(bmad_com%significant_length/100, in_len)
else
  sliced_ele%value(l$) = l_slice
endif
sliced_ele%lord_status = not_a_lord$
sliced_ele%slave_status = slice_slave$
sliced_ele%n_slave = 0
sliced_ele%n_slave_field = 0
sliced_ele%ix_ele = ix_slice_slave$  ! Indicate sliced ele is not an element in the lattice.
sliced_ele%s_start = ele_in%s - in_len + offset
sliced_ele%s = sliced_ele%s_start + sliced_ele%value(l$)
sliced_ele%map_ref_orb_in  = coord_struct()
sliced_ele%map_ref_orb_out = coord_struct()

! The sliced element is treated as a super_slave to the original element except if
! the original element is itself a super_slave. In this case the sliced element is a super_slave
! of the original elements lords.

if (ele_in%slave_status == slice_slave$) then
  sliced_ele%lord => ele_in%lord
else
  sliced_ele%lord => ele_in
endif

sliced_ele%n_lord = 1
if (ele_in%slave_status == super_slave$) sliced_ele%n_lord = ele_in%n_lord

! Err check. Remember: the element length may be negative

if (l_slice*in_len < 0 .or. abs(l_slice) > abs(in_len) + bmad_com%significant_length) then
  call out_io (s_fatal$, r_name, &
        'SLICE LENGTH IS OUT OF RANGE FOR ELEMENT: ' // ele_in%name, &
        'LENGTH: \2es12.3\ ', r_array = [l_slice, in_len])
  if (global_com%exit_on_error) call err_exit
  return
endif

if (ele_in%key == hybrid$) then
  call out_io (s_fatal$, r_name, &
        'CANNOT SLICE ELEMENT OF TYPE: ' // key_name(ele_in%key), &
        'CANNOT SLICE: ' // ele_in%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

! Simple case where ele length is zero

if (in_len == 0) then
  err_flag = .false.
  return
endif

!

call makeup_super_slave1 (sliced_ele, ele_in, offset, param, include_upstream_end, include_downstream_end, err2_flag)
if (err2_flag) return

! See if %rad_map can be saved. Can only do this if the old_slice is the same as the new one.
! Therefore cannot save with wigglers and other elements which are not uniform longitudinally

if (present(old_slice) .and. .not. rad_map_stale .and. .not. include_upstream_end .and. &
    .not. include_downstream_end .and. ele0%value(l$) == sliced_ele%value(l$) .and. &
    (ele_in%tracking_method == bmad_standard$ .or. ele_in%field_calc == bmad_standard$)) then
  select case (ele_in%key)
  case (sbend$, quadrupole$, sextupole$, octupole$, thick_multipole$)
    if (associated(sliced_ele%rad_map)) sliced_ele%rad_map%stale = .false.
  end select
endif

! For a sliced taylor element the %taylor%term components point to the lord components. 
! The routine deallocate_ele_pointers will only nullify and never deallocate these components of a slice_slave.

if (ele_in%key == taylor$) then
  do i = 1, 6
    sliced_ele%taylor(i)%term => ele_in%taylor(i)%term
  enddo
endif

! Use a speedier tracking method.

select case (sliced_ele%tracking_method)
case (taylor$, symp_lie_ptc$)
  if (sliced_ele%field_calc == fieldmap$) then
    select case (sliced_ele%key)
    case (wiggler$, undulator$); sliced_ele%tracking_method = symp_lie_bmad$
    case default;                sliced_ele%tracking_method = symp_lie_ptc$
    end select
  else
    select case (sliced_ele%key)
    case (wiggler$, undulator$); sliced_ele%tracking_method = symp_lie_bmad$
    case (em_field$);            sliced_ele%tracking_method = symp_lie_ptc$
    case default;                sliced_ele%tracking_method = bmad_standard$
    end select
  endif
end select

select case (sliced_ele%mat6_calc_method)
case (taylor$, symp_lie_ptc$)
  if (sliced_ele%field_calc == fieldmap$) then
    select case (sliced_ele%key)
    case (wiggler$, undulator$); sliced_ele%mat6_calc_method = symp_lie_bmad$
    case default;                sliced_ele%mat6_calc_method = symp_lie_ptc$
    end select
  else
    select case (sliced_ele%key)
    case (wiggler$, undulator$); sliced_ele%mat6_calc_method = symp_lie_bmad$
    case (em_field$);            sliced_ele%mat6_calc_method = symp_lie_ptc$
    case default;                sliced_ele%mat6_calc_method = bmad_standard$
    end select
  endif
end select

sliced_ele%field_calc = refer_to_lords$

! Makeup_super_slave1 does not compute reference energy or time so need to do it here.

if (offset == 0) then
  ele0%value(p0c$)      = ele_in%value(p0c_start$)
  ele0%value(e_tot$)    = ele_in%value(e_tot_start$)
  ele0%ref_time         = ele_in%value(ref_time_start$)
  sliced_ele%time_ref_orb_in = ele_in%time_ref_orb_in

elseif (present(old_slice)) then
  sliced_ele%time_ref_orb_in = time_ref_orb_out

elseif (ele_has_constant_ds_dt_ref(ele_in)) then
  ele0%value(p0c$)      = ele_in%value(p0c$)
  ele0%value(e_tot$)    = ele_in%value(e_tot$)
  ele0%ref_time = ele_in%ref_time - ele_in%value(delta_ref_time$) * (ele_in%value(l$) - offset) / ele_in%value(l$)
  sliced_ele%time_ref_orb_in%vec = 0

else
  call transfer_ele (sliced_ele, ele2)
  call create_element_slice (ele2, ele_in, offset, 0.0_rp, param, .true., .false., err2_flag)
  if (err2_flag) return
  ele0%value(p0c$)      = ele2%value(p0c$)
  ele0%value(e_tot$)    = ele2%value(e_tot$)
  ele0%ref_time         = ele2%ref_time
  sliced_ele%time_ref_orb_in = ele2%time_ref_orb_out
endif

ele0%ref_species = ele_in%ref_species
call ele_compute_ref_energy_and_time (ele0, sliced_ele, param, err2_flag)
if (err2_flag) return

if (.not. include_upstream_end) sliced_ele%time_ref_orb_in%location = inside$
if (.not. include_downstream_end) sliced_ele%time_ref_orb_out%location = inside$

! Round off can throw off the ending ref energy.
! This can be problematic when include_downstream_end = T since a particle traveling through and
! into the next element will have an energy mismatch and track1 flag this as an error.
! Solution is to just set the end ref energy to what it should be.

if (include_downstream_end) then
  sliced_ele%value(p0c$)      = ele_in%value(p0c$)
  sliced_ele%value(e_tot$)    = ele_in%value(e_tot$)
endif

! Round off can throw off the beginning ref energy. 
! Adjust if necessary.

if (present(orb_in)) then
  if (orb_in%time_dir * orb_in%direction == 1) then
    p0c_set = sliced_ele%value(p0c_start$)
  else
    p0c_set = sliced_ele%value(p0c$)
  endif

  r = orb_in%p0c / p0c_set
  orb_in%vec(2) = orb_in%vec(2) * r
  orb_in%vec(4) = orb_in%vec(4) * r
  orb_in%vec(6) = orb_in%vec(6) * r + (r - 1.0_rp)
  orb_in%p0c = p0c_set
endif

err_flag = .false.

end subroutine create_element_slice


