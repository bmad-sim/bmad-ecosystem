!+
! Subroutine make_hybrid_lat (lat_in, keep_ele,
!                    remove_markers, lat_out, ix_out, use_taylor, orb0)
!
! Subroutine to concatinate together the elements in a lat to make
! a lat with fewer elements. This is used to speed up computation times.
! The concatinated elements in the new lat are known as hybrid elements.
!
! Note: Lat_out must not be the same actual argument as lat_in.
!
! Note: The routine update_hybrid_list can be useful in constructing the keep_ele array.
!
! Note: For hybrid elements lat_out%ele(i)%tracking_method and 
! lat_out%ele(i)%mat6_calc_method are set as follows:
!
!   use_taylor    tracking_method     mat6_calc_method
!   ----------    ---------------     ----------------  
!   False         linear$             no_method$
!   True          taylor$             taylor$
!
! The hybrid elements will have the following parameters defined:
!   L                 Length of the combinded elements.
!   E_TOT_START       Starting energy.
!   DELTA_E           Change in energy through the hybrid.
!   DELTA_REF_TIME    Time needed for the reference particle to transverse the hybrid.
!
! Note: For use_taylor = .false. You need to have made the 
! lat_in%ele(:)%mat6 matrices before you call this routine.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in         -- lat_struct: Input lattice.
!   keep_ele(lat_in%n_ele_max) 
!                  -- Logical array: keep_ele(I) = True indicates an element
!                        that is NOT to be concatenated. That is, there will be
!                        a corresponding element in lat_out.
!   remove_markers -- Logical: If .true. then marker elements in keep_ele
!                        are removed from lat_out. in this case ix_out
!                        points to the element before the marker
!   use_taylor     -- Logical, optional: If present and True then the
!                        hybrid elements will have a taylor series 
!                        instead of a simple linear matrix. If an element to
!                        be concatenated has a taylor series then this taylor
!                        series will be concatenated with the other elements
!                        in the hybrid element. 
!   orb0(0:)       -- Coord_struct, optional: Central orbit for taylor stuff.
!
! Output:
!   lat_out   -- lat_struct: Lattice with hybrid elements. 
!                  Note: Lat_out must not be the same actual argument as lat_in.
!   ix_out(lat_in%n_ele_max) 
!             -- Integer, optional: Ix_out(i) is the index for lat_in%ele(i) 
!                  of the corresponding element in lat_out%ele(). 
!                  ix_out(i) set to 0 if lat_in%ele(i) is concatenated.
!-

#include "CESR_platform.inc"

subroutine make_hybrid_lat (r_in, keep_ele, remove_markers, &
                                       r_out, ix_out, use_taylor, orb0)

use ptc_interface_mod, except_dummy => make_hybrid_lat
use bmad_utils_mod

implicit none

type (lat_struct), target :: r_in, r_out
type (coord_struct), optional, volatile :: orb0(0:)
type (coord_struct) c0, c2
type (ele_struct), pointer :: ele_in, ele_out

real(rp) e_vec(4), ref_time0

integer j_in, i_out, k, n
integer j, ix, ic, o_key, n_con, n_ic, n_lord
integer, allocatable, save :: ica(:)
integer, optional :: ix_out(:)

logical init_hybrid_needed, remove_markers, keep_ele(:)
logical z_decoupled, do_taylor
logical, optional :: use_taylor

! Init

n = count(keep_ele(1:r_in%n_ele_max))
call init_lat (r_out, n+100)

if (present(use_taylor)) then
  do_taylor = use_taylor
else
  do_taylor = .false.
endif

if (all(r_in%ele(1:r_in%n_ele_track)%mat6(6,5) == 0) .and. .not. do_taylor) then
  z_decoupled = .true.
else
  z_decoupled = .false.
endif

i_out = 0                        ! index for current out lat
r_out%ele(0) = r_in%ele(0)     !
init_hybrid_needed = .true.         ! we need to init out lat element

if (r_in%n_ele_track == 0) then
  print *, 'ERROR IN make_hybrid_lat: LAT_IN%n_ele_track = 0!'
  call err_exit
endif

! loop over all in lat elements

do j_in = 1, r_in%n_ele_track

  ele_in => r_in%ele(j_in)

  ! if a match...

  if (keep_ele(j_in)) then

    ! if current out-element is a hybrid then calculate dispersion part of mat6

    if (i_out /= 0) then
      if (ele_out%key == hybrid$ .and. z_decoupled) &
                      call mat6_dispersion (e_vec, ele_out%mat6)
    endif

    ! on to the next out-element which is a simple element

    if (remove_markers .and. (ele_in%key == marker$ .or. &
              ele_in%key == photon_branch$ .or. ele_in%key == branch$)) then
      ele_in%ixx = i_out
    else
      i_out = i_out + 1                     ! starting next element
      if (i_out > ubound(r_out%ele, 1)) call allocate_lat_ele_array(r_out)
      ele_out => r_out%ele(i_out)
      ele_out = ele_in   ! single element
      ele_in%ixx = i_out
    endif

    init_hybrid_needed = .true.                ! need to init next ele

  ! here if no match found...
  ! If this is the first element after a matched element then just transfer in
  ! to out. Else modify the out mat6 transfer matrix using the in MAT6 matrix.

  else

    ele_in%ixx = 0          ! point to nothing

    if (init_hybrid_needed) then
      i_out = i_out + 1                       ! starting next element
      if (i_out > ubound(r_out%ele, 1)) call allocate_lat_ele_array(r_out)
      ele_out => r_out%ele(i_out)
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
      ele_out%mat6_calc_method = no_method$
      ele_out%value(e_tot_start$) = r_in%ele(j_in-1)%value(e_tot$)
      ref_time0                   = r_in%ele(j_in-1)%ref_time
      ele_out%value(delta_e$)        = ele_in%value(e_tot$) - ele_out%value(e_tot_start$)
      ele_out%value(delta_ref_time$) = ele_in%ref_time - ref_time0 

      if (present (orb0)) then
        c0 = orb0(j_in)
      else
        c0%vec = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      endif

      if (z_decoupled) then
        e_vec(1:4) = ele_out%mat6(1:4, 6)
      elseif (do_taylor) then
        ele_out%tracking_method = taylor$
        ele_out%mat6_calc_method = taylor$
        if (.not. associated(ele_out%taylor(1)%term)) then ! construct taylor
          call ele_to_taylor (ele_out, r_in%param, c0)
        endif
      endif

      init_hybrid_needed = .false.

    else
      if (ele_in%key == marker$ .or. ele_in%key == photon_branch$ .or. &
                                                     ele_in%key == branch$) cycle

      if (z_decoupled) then
        e_vec = matmul(ele_in%mat6(1:4,1:4), e_vec)
        e_vec(1:4) = e_vec(1:4) + ele_in%mat6(1:4, 6)
        ele_out%mat6(1:4,1:4) = matmul(ele_in%mat6(1:4,1:4), &
                                                    ele_out%mat6(1:4,1:4))
      elseif (do_taylor) then
        if (associated(ele_in%taylor(1)%term)) then
          call concat_ele_taylor (ele_out%taylor, ele_in, ele_out%taylor)
        else
          call taylor_propagate1 (ele_out%taylor, ele_in, r_in%param)
        endif
      else
        ele_out%mat6 = matmul(ele_in%mat6, ele_out%mat6)
      endif

      ele_out%s = ele_in%s
      ele_out%ref_time = ele_in%ref_time
      ele_out%value(l$) = ele_out%value(l$) + ele_in%value(l$)
      if (ele_in%value(hkick$) /= 0 .or. ele_in%value(vkick$) /= 0) then
        c2%vec = 0
        call offset_particle (ele_in, r_in%param, c2, set$, &
                        set_canonical = .false., set_multipoles = .false.)
        call offset_particle (ele_in, r_in%param, c2, unset$, &
                        set_canonical = .false., set_multipoles = .false.)
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
      if (ele_in%key == drift$ .and. (o_key == drift$ .or. &
           o_key == marker$ .or. o_key == photon_branch$ .or. o_key == branch$)) then
        ele_out%name = 'DRIFT_HYBRID' 
        ele_out%key = drift$
      else
        ele_out%name = 'HYBRID'
        ele_out%key = hybrid$
      endif

      ele_out%a       = ele_in%a
      ele_out%b       = ele_in%b
      ele_out%c_mat   = ele_in%c_mat
      ele_out%gamma_c = ele_in%gamma_c

    endif

    if (ele_out%key == hybrid$ .and. .not. do_taylor) then
      if (present(orb0)) then
        ele_out%vec0 = orb0(j_in)%vec - matmul(ele_out%mat6, c0%vec)
      else
        ele_out%vec0 = 0
      endif
    endif

  endif ! keep_ele

enddo

! end cleanup

if (ele_out%key == hybrid$ .and. z_decoupled)  &
                        call mat6_dispersion (e_vec, ele_out%mat6)

call transfer_lat_parameters (r_in, r_out)
r_out%n_ele_track  = i_out
r_out%n_ele_track = i_out

! put control elements in

do j_in = r_in%n_ele_track+1, r_in%n_ele_max
  ele_in => r_in%ele(j_in)    
  if (keep_ele(j_in)) then
    i_out = i_out + 1
    if (i_out > ubound(r_out%ele, 1)) call allocate_lat_ele_array(r_out)
    ele_out => r_out%ele(i_out)
    ele_in%ixx = i_out
    ele_out = ele_in
  endif
enddo
r_out%n_ele_max = i_out

!--------------------------------
! Construct control info...
! First construct r_out%control(:) array

allocate (ica(size(r_in%control)))

n_con = 0
do i_out = 1, r_out%n_ele_max

  ele_out => r_out%ele(i_out)
  if (ele_out%n_slave == 0) cycle

  n_con = n_con + ele_out%n_slave
  if (n_con > size(r_out%control)) call reallocate_control(r_out, n_con + 500)

  do j = ele_out%ix1_slave, ele_out%ix2_slave
    k = n_con + j - ele_out%ix2_slave
    ica(j) = k
    ix = r_in%control(j)%ix_slave
    if (r_in%ele(ix)%ixx == 0) then
      r_out%control(k)%ix_lord = i_out
      r_out%control(k)%ix_slave = ubound(r_out%ele, 1) ! point to dummy ele
      r_out%control(k)%ix_attrib = -1
      r_out%control(k)%coef = 0
    else
      r_out%control(k)%ix_lord = i_out
      r_out%control(k)%ix_slave = r_in%ele(ix)%ixx
      r_out%control(k)%ix_attrib = r_in%control(j)%ix_attrib
      r_out%control(k)%coef = r_in%control(j)%coef
    endif
  enddo

  ele_out%ix1_slave = n_con - ele_out%n_slave + 1
  ele_out%ix2_slave = n_con

enddo

r_out%n_control_max = n_con

! Second construct r_out%ic array
! Some lords might not exist so we need to be carful.

n_ic = 0

do i_out = 1, r_out%n_ele_max
  
  ele_out => r_out%ele(i_out)
  if (ele_out%n_lord == 0) cycle

  n_lord = 0  ! Actual number of lords
  do j = ele_out%ic1_lord, ele_out%ic2_lord
    ix = r_in%control(r_in%ic(j))%ix_lord
    if (.not. keep_ele(ix)) cycle
    n_ic = n_ic + 1
    n_lord = n_lord + 1
    r_out%ic(n_ic) = ica(r_in%ic(j))
    ic = r_out%ic(n_ic)
    ix = r_in%control(ic)%ix_lord
    if (ix == 0) then
      print *, 'WARNING IN make_hybrid_lat: LORD ELEMENT ',  &
                  'LIST NOT COMPLETE FOR: ', ele_out%name
    endif
  enddo

  ele_out%n_lord = n_lord
  ele_out%ic1_lord = n_ic - n_lord + 1
  ele_out%ic2_lord = n_ic
  if (ele_out%n_lord == 0) ele_out%slave_status = free$

enddo

r_out%n_ic_max = n_ic
deallocate (ica)

! end

if (r_out%n_ele_track == 0) then
  print *, 'ERROR IN make_hybrid_lat: OUTPUT LAT HAS 0 ELEMENTS!'
  call err_exit
endif

r_out%ele_init = r_in%ele_init

call check_lat_controls (r_out, .false.)
if (present (ix_out)) ix_out(1:r_in%n_ele_max) = r_in%ele(1:r_in%n_ele_max)%ixx

end subroutine
