!+
! subroutine tao_data_show_use(d2_data)
!
! Displays what data is used by the optimizer for the specied d2_data
!
! Input:
!   d2_data	-- tao_d2_data_struct: 
!-

subroutine tao_data_show_use (d2_data)

use tao_mod

implicit none

type (tao_d2_data_struct), intent(in) :: d2_data

character(17) :: r_name = "tao_data_show_use"
character(200) line

integer i

! print out used data for each name

do i = 1, size(d2_data%d1)
  call location_encode (line, d2_data%d1(i)%d%useit_opt, &
                      d2_data%d1(i)%d%exists, lbound(d2_data%d1(i)%d, 1))
  write (line, '(2x, 3a, t20, 2a)') trim(d2_data%name), &
                  ':',  trim(d2_data%d1(i)%name), "Using: ", line(1:170)
  call out_io (s_blank$, r_name, line)
enddo

end subroutine tao_data_show_use
