!+
! subroutine tao_use_data (action, data_name)
!
! Veto, restore or use specified datums. range syntax is just like
!    indexing in fortran: 1:34,46,58:78
!
! Input:
!   action	     -- character(*): veto, use or restore
!   data_name    -- character(*): the selected data
!
! Output:
!-

subroutine tao_use_data (action, data_name)

use tao_mod

implicit none

type (tao_d1_data_array_struct), allocatable, save :: d1_dat(:)
type (tao_data_array_struct), allocatable, save :: d_dat(:)

character(*) :: action
character(*) :: data_name

integer i, which

logical err

character(12) :: r_name = "tao_use_data"

! decipher action

call match_word (action, name$%use_veto_restore, which)

! If "use" is choisen then must veto everything first.

call tao_find_data (err, data_name, d1_array = d1_dat, d_array = d_dat)

if (which == use$) then
  if (err) return
  do i = 1, size(d1_dat)
    d1_dat(i)%d1%d%good_user = .false.
  enddo
endif

! now do the set

do i = 1, size(d_dat)
  select case (which)
  case (use$, restore$)
    d_dat(i)%d%good_user = .true. 
  case (veto$)
    d_dat(i)%d%good_user = .false.
  case default
    call out_io (s_error$, r_name, "Internal error picking name$%use_veto_restore")
    err = .true.
  end select
enddo

! Optimizer bookkeeping and Print out changes.

call tao_set_data_useit_opt()

do i = 1, size(d1_dat)
  if (i > 1) then
    if (associated(d1_dat(i)%d1%d2, d1_dat(i-1)%d1%d2)) cycle
  endif
  call tao_data_show_use (d1_dat(i)%d1%d2)
enddo

end subroutine tao_use_data
