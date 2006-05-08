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

type (tao_data_array_struct), allocatable, save :: d_dat(:)

character(*) :: action
character(*) :: data_name

integer i, which

logical err

character(12) :: r_name = "tao_use_data"

! decipher action

call match_word (action, name$%use_veto_restore, which)

! If "use" is choisen then must veto everything first.

if (which == use$) then
  call tao_find_data (err, data_name, d_array = d_dat, all_elements = .true.)
  if (err) return
  do i = 1, size(d_dat)
    d_dat(i)%d%good_user = .false.
  enddo
endif

! now do the set

call tao_find_data (err, data_name, d_array = d_dat)
if (err) return

do i = 1, size(d_dat)
  select case (which)
  case (use$, restore$)
    d_dat(i)%d%good_user = .true. 
  case (veto$)
    d_dat(i)%d%good_user = .false.
  case default
    call out_io (s_error$, r_name, &
                    "Internal error picking name$%use_veto_restore")
    err = .true.
  end select
enddo

! Optimizer bookkeeping and Print out changes.

call tao_set_data_useit_opt()

call tao_data_show_use (d_dat(1)%d%d1%d2)
do i = 2, size(d_dat)
  if (.not. associated(d_dat(i)%d%d1%d2, d_dat(i-1)%d%d1%d2)) &
        call tao_data_show_use (d_dat(i)%d%d1%d2)
enddo

end subroutine tao_use_data
