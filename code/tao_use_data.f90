!+
! subroutine tao_use_data (action, data_type, locations)
!
! Veto, restore or use specified datums. range syntax is just like
!    indexing in fortran: 1:34, 46, 58:78
!
! Input:
!   action	     -- character(*): veto, use or restore
!   data_type        -- character(*): the selected data name or all
!   locations        -- character(*): the index location expression
!
! Output:
!-

subroutine tao_use_data (action, data_type, locations)

use tao_mod

implicit none

type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr

character(*) :: action
character(*) :: data_type
character(*) :: locations

logical, allocatable :: action_logic(:) !which elements do we take action on?
logical, automatic :: picked(size(s%u))
logical err, all_selected

integer which, i, iu, n1, n2, ix1, ix2, j, ix
integer err_num

character(3) select_all
character(12) :: r_name = "tao_use_data"
character(40) d_name, d1_name
character(200) line

! decipher action

call match_word (action, name$%use_veto_restore, which)

! Are we selecting all data?
call str_upcase(select_all, data_type(1:3))
all_selected = .false.
if (index(select_all, 'ALL') .ne. 0) all_selected = .true.

! loop over the universes to do.

call tao_pick_universe (data_type, d_name, picked, err)
if (err) return

do iu = 1, size(s%u)

  ! find data name and name

  if (.not. picked(iu)) cycle
  if (all_selected) then
    if (locations .eq. ' ') locations = "all"
    do j = 1, s%u(iu)%n_d2_data_used
      d2_ptr => s%u(iu)%d2_data(j)
      ix = index(d_name, ':')
      if (ix .ne. 0) then !do only specified dimension
        d1_name = data_type(ix+1:)
        d1_loop: do i = 1, size(d2_ptr%d1)
          if (d1_name == d2_ptr%d1(i)%name) then
            d1_ptr => d2_ptr%d1(i)
            exit d1_loop
          endif
          if (i .eq. size(d2_ptr%d1)) then
            call out_io (s_error$, r_name, "Couldn't find d1_data name: " // d1_name)
            cycle 
          endif
        enddo d1_loop
      else ! we want all dimensions
        nullify(d1_ptr)
      endif
      call use_d2_data () ! with d1_ptr and d2_ptr set
      if (err_num == -1) return
    enddo
  else
    call tao_find_data (err, s%u(iu), d_name, d2_ptr, d1_ptr)
    if (err) return
    call use_d2_data () ! with d1_ptr and d2_ptr set
    if (err_num == -1) return
  endif
enddo

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

!make sure d2_ptr and d1_ptr is set properly before calling this!
subroutine use_d2_data ()


! find locations
 
! set d%good_user based on action and action_logic

if (associated(d1_ptr)) then
  n1 = lbound(d1_ptr%d, 1)
  n2 = ubound(d1_ptr%d, 1)
  allocate(action_logic(n1:n2))
  line = locations
  call location_decode (line, action_logic, n1, err_num) 
  if (err_num == -1) then
    deallocate(action_logic)
    return
  endif
  call use (d1_ptr)
  deallocate(action_logic)

else
  do i = 1, size(d2_ptr%d1)
    n1 = lbound(d2_ptr%d1(i)%d, 1)
    n2 = ubound(d2_ptr%d1(i)%d, 1)
    allocate(action_logic(n1:n2))
    line = locations
    call location_decode (line, action_logic, n1, err_num) 
    if (err_num == -1) then
      deallocate(action_logic)
      return
    endif
    call use (d2_ptr%d1(i))
    deallocate(action_logic)
    if (err) return
  enddo

endif


! Optimizer bookkeeping and Print out changes.
call tao_set_data_useit_opt()
call tao_data_show_use (d2_ptr)

end subroutine use_d2_data

!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine use (d1)

type (tao_d1_data_struct) d1
integer i

!

select case (which)
case (use$)
  d1%d(:)%good_user = action_logic
case (veto$)
  where (action_logic) d1%d(:)%good_user = .not. action_logic
case (restore$)
  where (action_logic) d1%d(:)%good_user = action_logic
case default
  call out_io (s_error$, r_name, &
                    "Internal error picking name$%use_veto_restore")
  err = .true.
end select

end subroutine

end subroutine tao_use_data
