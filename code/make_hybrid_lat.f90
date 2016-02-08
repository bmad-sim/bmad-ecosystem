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
!
! For elements in the tracking part of the lat_in lattice (that is, for non-lord elements):
!    If %select is True, the element is retained in lat_out and not hybridized.
!    If %select is False, the element will be hybridized except if the element is needed to preserve
!      lattice bookkeeping (discussed below) or if the elements to either side of it are not being 
!      hybridized so that there are no neighboring elements to hybridize with.
!
! For lord elements in the lat_in lattice:
!    If %select is True, the element will be retained in lat_out. In this case, all slave elements of this
!      element (and slaves of slaves, etc) will also be retained.
!    If %select is False, the element may be retained or removed in lat_out using the folowing rules:
!      If the element can be retained without changing what elements will be retained or removed/hybridized then
!        the element will be retained.
!      If removal of the element will cause lattice bookkeeping problemsexcept if a sub-slave is being hybridized.
!      If a sub-slave is being hybridized,  in which case the element will be removed to preserve lattice bookkeeping.
!
! For hybrid elements lat_out%ele(i)%tracking_method and 
! lat_out%ele(i)%mat6_calc_method are set as follows:
!
!   use_taylor    tracking_method     mat6_calc_method
!   ----------    ---------------     ----------------  
!   False         linear$             static$
!   True          taylor$             taylor$
!
! With use_taylor = False, before this routine is called, the lat_in%ele(:)%mat6 matrices must be calculated.
!
! The hybrid elements will have the following parameters defined:
!   L                 Length of the combinded elements.
!   E_TOT_START       Starting energy.
!   DELTA_E           Change in energy through the hybrid.
!   DELTA_REF_TIME    Time needed for the reference particle to transverse the hybrid.
!
! Note: Lat_out must not be the same actual argument as lat_in.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in         -- lat_struct: Input lattice.
!     %branch(:)%ele(:)%select  -- Set this element component to enable hybridization. See above
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
type (ele_struct), pointer :: ele_in, ele_out, lord, lord2, slave
type (branch_struct), pointer :: b_in, b_out
type (control_struct), pointer :: ctl1, ctl2

real(rp) ref_time0

integer j_in, i_out, k, n, ix_hyb
integer j, ix, ic, o_key, n_con, n_ic, n_lord, ib, ie, il, il2, is
integer, allocatable :: ica(:)

logical init_hybrid_needed
logical do_taylor, err_flag, something_has_changed
logical, optional :: use_taylor

character(*), parameter :: r_name = 'make_hybrid_lat'

! Init

lat_out = lat_in
do_taylor = logic_option(.false., use_taylor)

! Must keep all subslaves of a kept lord

do ie = lat_out%n_ele_track+1, lat_out%n_ele_max
  call select_slaves (lat_out%ele(ie))
enddo

! If a super_lord has %select = True then each of its slaves must have all of their lords havd %select = True.

do ie = lat_out%n_ele_track+1, lat_out%n_ele_max
  lord => lat_out%ele(ie)
  if (.not. lord%select) cycle
  if (lord%lord_status /= super_lord$) cycle

  do is = 1, lord%n_slave
    slave => pointer_to_slave(lord, is)
    do il = 1, slave%n_lord
      lord2 => pointer_to_lord(slave, il)
      if (lord2%lord_status /= super_lord$) cycle
      if (lord2%select) cycle
      lord2%select = .true.
      call select_slaves (lord2)
    enddo
  enddo
enddo

! If a controller (group or overlay) controls the same attribute of some slave as some other controller with 
! %select = True, then %select of this controller must be set to True.

do
  something_has_changed = .false.

  do ib = 0, ubound(lat_out%branch, 1)
    b_out => lat_out%branch(ib)
    b_out%ele(0)%select = .true.  ! The beginning marker ele always gets preserved
    do ie = 1, b_out%n_ele_track
      ele_out => b_out%ele(ie)
      do il = 1, ele_out%n_lord
        lord => pointer_to_lord(ele_out, il, ctl1)
        if (.not. lord%select) cycle
        if (lord%key /= overlay$ .and. lord%key /= group$) cycle
        do il2 = 1, ele_out%n_lord
          lord2 => pointer_to_lord(ele_out, il, ctl2)
          if (lord2%select) cycle
          if (ctl1%ix_attrib /= ctl2%ix_attrib) cycle
          lord2%select = .true.
          call select_slaves(lord2)
          something_has_changed = .true.
        enddo
      enddo
    enddo
  enddo

  if (.not. something_has_changed) exit
enddo  

! Keep all lords that do not have any slaves that will be hybridized.

do ie = lat_out%n_ele_track+1, lat_out%n_ele_max
  lord => lat_out%ele(ie)
  if (lord%select) cycle
  if (lord_has_slave_to_be_hybridized(lord)) cycle
  lord%select = .true.
enddo

! Keep all tracking elements that do not have elements around them to be hybridized.

do ib = 0, ubound(lat_out%branch, 1)
  b_out => lat_out%branch(ib)
  do ie = 2, b_out%n_ele_track-1
    ele_out => b_out%ele(ie)
    if (b_out%ele(ie-1)%select .and. b_out%ele(ie+1)%select) ele_out%select = .true.
  enddo
  lat_in%branch(ib)%ele(:)%select = b_out%ele(:)%select
enddo

! Mark elements for deletion and then delete them. 
! The first element in a group of elements to be hybridized does not get deleted and
! this element will then be transformed into a hybrid element.

do ib = 0, ubound(lat_out%branch, 1)
  b_out => lat_out%branch(ib)
  do ie = 1, b_out%n_ele_max
    ele_out => b_out%ele(ie)
    if (ele_out%select) cycle
    if (ie <= b_out%n_ele_track .and. b_out%ele(ie-1)%select) cycle
    ele_out%key = -1   ! mark for delection
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

  ! loop over all in lat elements

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
    ! to out. Else modify the ele_out%mat6 matrix using the ele_in%mat6 matrix.

    else

      ele_in%ixx = 0          ! point to nothing

      if (init_hybrid_needed) then
        i_out = i_out + 1                       ! starting next element
        ele_out => b_out%ele(i_out)
        ele_out = ele_in
        ele_out%lord_status = not_a_lord$
        ele_out%slave_status = free$
        ele_out%n_slave = 0
        ele_out%ix1_slave = 0
        ele_out%ix2_slave = -1
        ele_out%n_lord = 0
        ele_out%ic1_lord = 0
        ele_out%ic2_lord = -1
        ele_out%tracking_method = linear$
        ele_out%mat6_calc_method = static$
        ele_out%value(e_tot_start$) = b_in%ele(j_in-1)%value(e_tot$)
        ref_time0                   = b_in%ele(j_in-1)%ref_time
        ele_out%value(ref_time_start$) = ref_time0
        ele_out%value(delta_e$)        = ele_in%value(e_tot$) - ele_out%value(e_tot_start$)
        ele_out%value(delta_ref_time$) = ele_in%ref_time - ref_time0 

        ix_hyb = ix_hyb + 1

        if (do_taylor) then
          if (present (orb0_arr)) then
            c0 = orb0_arr(ib)%orbit(j_in)
          else
            c0%vec = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
          endif
          ele_out%tracking_method = taylor$
          ele_out%mat6_calc_method = taylor$
          if (.not. associated(ele_out%taylor(1)%term)) then ! construct taylor
            call ele_to_taylor (ele_out, b_in%param, ele_out%taylor, c0)
          endif
        endif

        init_hybrid_needed = .false.

      else
        if (ele_in%key == marker$ .or. ele_in%key == photon_fork$ .or. ele_in%key == fork$) cycle

        if (do_taylor) then
          if (associated(ele_in%taylor(1)%term)) then
            call concat_ele_taylor (ele_out%taylor, ele_in, ele_out%taylor)
          else
            call taylor_propagate1 (ele_out%taylor, ele_in, b_in%param)
          endif
        else
          ele_out%mat6 = matmul(ele_in%mat6, ele_out%mat6)
          ele_out%vec0 = ele_in%vec0 + matmul(ele_in%mat6, ele_out%vec0)
        endif

        ele_out%s = ele_in%s
        ele_out%ref_time = ele_in%ref_time
        ele_out%value(l$) = ele_out%value(l$) + ele_in%value(l$)
        if (ele_in%value(hkick$) /= 0 .or. ele_in%value(vkick$) /= 0) then
          c2%vec = 0
          call offset_particle (ele_in, b_in%param, set$, c2, set_multipoles = .false.)
          call offset_particle (ele_in, b_in%param, unset$, c2, set_multipoles = .false.)
          ele_out%value(hkick$) = ele_out%value(hkick$) + c2%vec(2)
          ele_out%value(vkick$) = ele_out%value(vkick$) + c2%vec(4)
        endif

        ele_out%value(x1_limit$)       = ele_in%value(x1_limit$)
        ele_out%value(x2_limit$)       = ele_in%value(x2_limit$)
        ele_out%value(y1_limit$)       = ele_in%value(y1_limit$)
        ele_out%value(y2_limit$)       = ele_in%value(y2_limit$)
        ele_out%value(delta_e$)        = ele_in%value(e_tot$) - ele_out%value(e_tot_start$)
        ele_out%value(delta_ref_time$) = ele_in%ref_time - ref_time0

        o_key = ele_out%key 

        if (ele_in%key == drift$ .and. (o_key == drift$ .or. o_key == marker$ .or. o_key == photon_fork$ .or. o_key == fork$)) then
          ele_out%key = drift$
          write (ele_out%name, '(a, i0)') 'DRIFT_HYBRID', ix_hyb
        else
          ele_out%key = hybrid$
          write (ele_out%name, '(a, i0)') 'HYBRID', ix_hyb
        endif

        ele_out%a       = ele_in%a
        ele_out%b       = ele_in%b
        ele_out%c_mat   = ele_in%c_mat
        ele_out%gamma_c = ele_in%gamma_c

      endif

    endif ! ele%select

  enddo ! ele loop

  if (b_out%n_ele_track == 0) then
    call out_io (s_fatal$, r_name, 'OUTPUT LATTICE HAS ZERO TRACKING ELEMENTS!')
    if (global_com%exit_on_error) call err_exit
  endif

enddo ! branch loop

! End

call lat_sanity_check (lat_out, err_flag)

!------------------------------------------------------------------------------------
contains

recursive subroutine select_slaves (lord)

type (ele_struct) lord
type (ele_struct), pointer :: slave

integer i

!

if (.not. lord%select) return

do i = 1, lord%n_slave
  slave => pointer_to_slave (lord, i)
  slave%select = .true.
  call select_slaves(slave)
enddo

end subroutine select_slaves

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
