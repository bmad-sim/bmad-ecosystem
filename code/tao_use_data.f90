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

use tao_interface, dummy => tao_use_data

implicit none

type (tao_d1_data_array_struct), allocatable :: d1_dat(:)
type (tao_data_array_struct), allocatable :: d_dat(:)

character(*) :: action
character(*) :: data_name
character(16) match

integer i, which

logical err

character(12) :: r_name = "tao_use_data"

! decipher action

call match_word (action, ["use    ", "veto   ", "restore"], which, matched_name = match)

! If "use" is choisen then must veto everything first.

call tao_find_data (err, data_name, d1_array = d1_dat, d_array = d_dat)
if (err) return

if (match == 'use') then
  do i = 1, size(d1_dat)
    d1_dat(i)%d1%d%good_user = .false.
  enddo
endif

! Now do the set

do i = 1, size(d_dat)
  select case (match)
  case ('use', 'restore')
    d_dat(i)%d%good_user = .true. 
  case ('veto')
    d_dat(i)%d%good_user = .false.
  case default
    call out_io (s_error$, r_name, "Internal error picking use/veto/restore")
    return
  end select
enddo

! Optimizer bookkeeping 

do i = 1, size(d1_dat)
  call tao_set_data_useit_opt(d1_dat(i)%d1%d)
enddo

! And print changes.

do i = 1, size(d1_dat)
  if (i > 1) then
    if (associated(d1_dat(i)%d1%d2, d1_dat(i-1)%d1%d2)) cycle
  endif
  call tao_data_show_use (d1_dat(i)%d1%d2)
enddo

end subroutine tao_use_data
