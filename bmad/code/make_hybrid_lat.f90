!+
! Subroutine make_hybrid_lat (lat_in, lat_out, use_taylor, orb0_arr)
!
! Routine to concatinate together groups of elements in the tracking part of a lattice to 
! make a lattice with fewer elements to track through.
! This is used to speed up computation times.
! The concatinated elements in the new lattice are known as hybrid elements.
!
! The %select component of an element in lat_in is used to decide how an element is treated.
! [That is, the calling routine needs to set the %select component before calling this routine.]
! In general, %select = True means keep the element (do not hybridize) and False means hybridize.
! Additional rules:
!   * If %select = True for a lord, all the slave elements will be set with %select = True.
!   * If a super_slave has %select = True, all super_lords of this slave will be set with %select = True. 
!      And all super_slaves of these super_lords will also be so marked. This is done to prevent the 
!      hybridized lattice (lat_out) from being in an ill defined state.
!   * If two overlay elements control the same attribute of a given slave, and if one of the overlay elements
!      has %select = True, %select of the other overlay will be set to True.
!   * If the elements to either side of an element are not being hybridized then the element itself will
!      be retained (there is no advantage to hybridization in this case). Exception: This rule is not
!      applied to a super_slave elements.
!   * Any lord element where one or more slaves is hybridized will be removed from the lattice.
!      Exception: multipass_lords will be retained if one or more slaves are retained.
!
! With use_taylor = False, before this routine is called, the lat_in%ele(:)%mat6 matrices must be calculated.
!
! The hybrid elements will have the following parameters defined:
!   L                 Length of the combinded elements.
!   E_TOT_START       Starting energy.
!   DELTA_E_REF       Change in reference energy through the hybrid.
!   DELTA_REF_TIME    Time needed for the reference particle to transverse the hybrid.
!
! Note: Lat_out must not be the same actual argument as lat_in.
!
! Input:
!   lat_in         -- lat_struct: Input lattice.
!     %branch(:)%ele(:)%select  -- Roughly: Set True to keep and False to hybridize. See above.
!   use_taylor     -- Logical, optional: If present and True then the
!                        hybrid elements will have a taylor series 
!                        instead of a simple linear matrix. If an element to
!                        be concatenated has a taylor series then this taylor
!                        series will be concatenated with the other elements
!                        in the hybrid element. 
!   orb0_arr(0:)   -- Coord_array_struct, optional: Central orbit for taylor stuff.
!                       Each orb0_arr(i)%orbit(:) holds the orbit for the i^th lattice branch
!
! Output:
!   lat_out   -- lat_struct: Lattice with hybrid elements. 
!                  Note: Lat_out must not be the same actual argument as lat_in.
!-

subroutine make_hybrid_lat (lat_in, lat_out, use_taylor, orb0_arr)

use ptc_interface_mod, except_dummy => make_hybrid_lat

implicit none

type (lat_struct), target :: lat_in, lat_out
type (coord_array_struct), optional :: orb0_arr(0:)
type (coord_struct) c0, c2
type (ele_struct), pointer :: ele, ele_in, ele_out, lord, lord2, slave
type (branch_struct), pointer :: branch, b_in, b_out
type (control_struct), pointer :: ctl1, ctl2

real(rp) ref_time0

integer j_in, i_out, k, n, ix_hyb
integer j, ix, ic, o_key, n_con, n_ic, n_lord, ib, ie, il, il2, is, i
integer, allocatable :: ica(:)

logical init_hybrid_needed
logical do_taylor, err_flag, something_has_changed
logical, optional :: use_taylor

character(*), parameter :: r_name = 'make_hybrid_lat'

! Init

lat_out = lat_in
do_taylor = logic_option(.false., use_taylor)

! If a lord has %select = True then each of its slaves must have %select = True.

do ie = lat_out%n_ele_track+1, lat_out%n_ele_max
  lord => lat_out%ele(ie)
  if (.not. lord%select) cycle
  call select_this_and_slaves (lord)
enddo

! If a super_slave has %select = T then %select of all of its super_lords must be set to T.

do
  something_has_changed = .false.

  do ib = 0, ubound(lat_out%branch, 1)
    branch => lat_out%branch(ib)
    branch%ele(0)%select = .true.  ! The beginning marker ele always gets preserved
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (.not. ele%select) cycle
      if (ele%slave_status /= super_slave$) cycle
      do il = 1, ele%n_lord
        lord => pointer_to_lord(ele, il, ctl1)
        if (lord%select) cycle
        if (lord%lord_status /= super_lord$) cycle
        call select_this_and_slaves(lord)
        something_has_changed = .true.
      enddo
    enddo
  enddo

  if (.not. something_has_changed) exit
enddo  

! If a overlay with %select = T controls the same attribute of some slave as some other overlay,
! %select of this overlay must be set to True.

do
  something_has_changed = .false.

  do ib = 0, ubound(lat_out%branch, 1)
    branch => lat_out%branch(ib)
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (.not. ele%select) cycle

      do il = 1, ele%n_lord
        lord => pointer_to_lord(ele, il, ctl1)
        if (.not. lord%select) cycle
        if (lord%key /= overlay$) cycle
        do il2 = 1, ele%n_lord
          if (il == il2) cycle
          lord2 => pointer_to_lord(ele, il, ctl2)
          if (lord2%select) cycle
          if (ctl1%ix_attrib /= ctl2%ix_attrib) cycle
          lord2%select = .true.
          call select_this_and_slaves(lord2)
          something_has_changed = .true.
        enddo
      enddo
    enddo
  enddo

  if (.not. something_has_changed) exit
enddo  

! Keep all tracking elements that do not have elements around them to be hybridized.

do ib = 0, ubound(lat_out%branch, 1)
  branch => lat_out%branch(ib)
  do ie = 1, branch%n_ele_track-1
    ele => branch%ele(ie)
    if (ele%slave_status == super_slave$) cycle
    if (branch%ele(ie-1)%select .and. branch%ele(ie+1)%select) ele%select = .true.
  enddo
  lat_in%branch(ib)%ele(0:branch%n_ele_max)%select = branch%ele(0:branch%n_ele_max)%select
enddo

! Keep a lord only if all slaves have %select = T.
! Exception: multipass_lord will be kept if it has one slave that is kept

do ie = lat_out%n_ele_track+1, lat_out%n_ele_max
  lord => lat_out%ele(ie)
  if (lord%select) cycle
  if (lord%lord_status == multipass_lord$) then
    do i = 1, lord%n_slave
      slave => pointer_to_slave(lord, i, ctl1)
      if (slave%select) then
        lord%select = .true.
      else
        ctl1%attribute = 'REMOVE'
      endif
    enddo
  else
    if (lord_has_slave_to_be_hybridized(lord)) cycle
    lord%select = .true.
  endif
enddo

! Mark elements for deletion and then delete them. 

do ib = 0, ubound(lat_out%branch, 1)
  branch => lat_out%branch(ib)
  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%select) cycle
    ele%ix_ele = -1   ! mark for delection
  enddo
enddo

call remove_eles_from_lat (lat_out)

! Now hybridize elements as needed

ix_hyb = 0

do ib = 0, ubound(lat_out%branch, 1)
  b_out => lat_out%branch(ib)
  b_in  => lat_in%branch(ib)

  i_out = 0                        ! index for current out lat
  init_hybrid_needed = .true.      ! we need to init out lat element

  call init_coord (c0, particle = b_in%param%particle)
  call init_coord (c2, particle = b_in%param%particle)

  ! loop over all lat_in elements

  do j_in = 1, b_in%n_ele_track

    ele_in => b_in%ele(j_in)

    ! if not hybridizing...

    if (ele_in%select) then

      ! On to the next out-element which is a simple element

      i_out = i_out + 1                     ! starting next element
      ele_out => b_out%ele(i_out)
      ele_in%ixx = i_out

      if (ele_out%name /= ele_in%name) then
        call out_io (s_fatal$, r_name, 'HYBRIDIZATION INTERNAL ERROR. PLASE CONTACT A BMAD MAINTAINER.')
        if (global_com%exit_on_error) call err_exit
      endif

      init_hybrid_needed = .true.                ! need to init next ele

    ! here if no match found...
    ! If this is the first element after a matched element then just transfer in
    ! to out.

    else

      ele_in%ixx = 0          ! point to nothing

      if (init_hybrid_needed) then
        i_out = i_out + 1                       ! starting next element
        call insert_element(lat_out, ele_in, i_out, ib)
        ele_out => b_out%ele(i_out)
        ele_out%lord_status = not_a_lord$
        ele_out%slave_status = free$
        ele_out%n_slave = 0
        ele_out%n_slave_field = 0
        ele_out%ix1_slave = 0
        ele_out%n_lord = 0
        ele_out%n_lord_field = 0
        ele_out%ic1_lord = 0
        ele_out%value(e_tot_start$) = ele_in%value(e_tot_start$)
        ele_out%value(p0c_start$)   = ele_in%value(p0c_start$)
        ref_time0                   = b_in%ele(j_in-1)%ref_time
        ele_out%value(ref_time_start$) = ref_time0
        ele_out%value(delta_e_ref$)    = ele_in%value(e_tot$) - ele_out%value(e_tot_start$)
        ele_out%value(delta_ref_time$) = ele_in%ref_time - ref_time0 
        ele_out%value(e_tot$)    = ele_in%value(e_tot$)
        ele_out%value(p0c$)      = ele_in%value(p0c$)
        ele_out%tracking_method  = taylor$
        ele_out%mat6_calc_method = taylor$
        ele_out%field_calc       = bmad_standard$

        ix_hyb = ix_hyb + 1

        if (do_taylor) then
          if (present (orb0_arr)) then
            c0 = orb0_arr(ib)%orbit(j_in-1)
          else
            c0%vec = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
          endif
          if (.not. associated(ele_out%taylor(1)%term)) then ! construct taylor
            call ele_to_taylor (ele_out, c0)
          endif
        endif

        init_hybrid_needed = .false.

      else
        if (ele_in%key == marker$ .or. ele_in%key == photon_fork$ .or. ele_in%key == fork$) cycle

        if (do_taylor) then
          if (associated(ele_in%taylor(1)%term)) then
            call concat_ele_taylor (ele_out%taylor, ele_in, err_flag)
          else
            call taylor_propagate1 (ele_out%taylor, ele_in, b_in%param, err_flag)
          endif
        else
          ele_out%mat6 = matmul(ele_in%mat6, ele_out%mat6)
          ele_out%vec0 = ele_in%vec0 + matmul(ele_in%mat6, ele_out%vec0)
        endif

        ele_out%s = ele_in%s
        ele_out%ref_time = ele_in%ref_time
        ele_out%value(l$) = ele_out%value(l$) + ele_in%value(l$)
        ele_out%s_start = ele_out%s - ele_out%value(l$)
        if (ele_in%value(hkick$) /= 0 .or. ele_in%value(vkick$) /= 0) then
          c2%vec = 0
          call offset_particle (ele_in, set$, c2)
          call offset_particle (ele_in, unset$, c2)
          ele_out%value(hkick$) = ele_out%value(hkick$) + c2%vec(2)
          ele_out%value(vkick$) = ele_out%value(vkick$) + c2%vec(4)
        endif

        ele_out%value(x1_limit$)       = ele_in%value(x1_limit$)
        ele_out%value(x2_limit$)       = ele_in%value(x2_limit$)
        ele_out%value(y1_limit$)       = ele_in%value(y1_limit$)
        ele_out%value(y2_limit$)       = ele_in%value(y2_limit$)

        ele_out%value(delta_e_ref$)    = ele_in%value(E_tot$) - ele_out%value(e_tot_start$)
        ele_out%value(delta_ref_time$) = ele_in%ref_time - ref_time0
        ele_out%value(e_tot$)          = ele_in%value(e_tot$)
        ele_out%value(p0c$)            = ele_in%value(p0c$)

        o_key = ele_out%key 

        if (ele_in%key == drift$ .and. (o_key == drift$ .or. o_key == marker$ .or. o_key == photon_fork$ .or. o_key == fork$)) then
          ele_out%key = drift$
          write (ele_out%name, '(a, i0)') 'DRIFT_HYBRID', ix_hyb
        else
          if (do_taylor) then
            ele_out%key = taylor$
          else
            ele_out%key = hybrid$
          endif
          write (ele_out%name, '(a, i0)') 'HYBRID', ix_hyb
        endif

        ele_out%a       = ele_in%a
        ele_out%b       = ele_in%b
        ele_out%c_mat   = ele_in%c_mat
        ele_out%gamma_c = ele_in%gamma_c
      endif

      if (.not. do_taylor) then
        call mat6_to_taylor(ele_out%vec0, ele_out%mat6, ele_out%taylor)
      endif

    endif ! ele%select

  enddo ! ele loop

  if (b_out%n_ele_track == 0) then
    call out_io (s_fatal$, r_name, 'OUTPUT LATTICE HAS ZERO TRACKING ELEMENTS!')
    if (global_com%exit_on_error) call err_exit
  endif

enddo ! branch loop

! End

call create_lat_ele_nametable (lat_out, lat_out%nametable)
call lat_sanity_check (lat_out, err_flag)

!------------------------------------------------------------------------------------
contains

recursive subroutine select_this_and_slaves (lord)

type (ele_struct) lord
type (ele_struct), pointer :: slave

integer i

!

lord%select = .true.

do i = 1, lord%n_slave
  slave => pointer_to_slave (lord, i)
  slave%select = .true.
  call select_this_and_slaves(slave)
enddo

end subroutine select_this_and_slaves

!------------------------------------------------------------------------------------
! contains

recursive function lord_has_slave_to_be_hybridized (lord) result (has_hy_slave)

type (ele_struct) lord
type (ele_struct), pointer :: slave

logical has_hy_slave
integer i

!

has_hy_slave = .false.

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  if (slave%n_slave == 0) then
    has_hy_slave = (.not. slave%select)
  else
    has_hy_slave = lord_has_slave_to_be_hybridized(slave)
  endif
  if (has_hy_slave) return
enddo

end function lord_has_slave_to_be_hybridized

end subroutine make_hybrid_lat
