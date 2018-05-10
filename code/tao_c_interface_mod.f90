module tao_c_interface_mod
use iso_c_binding
use tao_interface
use fortran_cpp_utils

implicit none

contains

!------------------------------------------------------------------------
!+ 
! function tao_c_out_io_buffer_num_lines() bind(c) result(n_lines)
!
! Function to access out_io_buffer%n_lines from C
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   n_lines -- integer(c_int): C access to out_io_buffer%n_lines
!-

function tao_c_out_io_buffer_num_lines() bind(c) result(n_lines)
integer(c_int) ::  n_lines
n_lines = out_io_buffer_num_lines()
end function

!------------------------------------------------------------------------
!+ 
! function tao_c_out_io_buffer_get_line(n) bind(c) result(c_string_ptr)
!
! Function to access out_io_buffer%line(n) from C
!
! Modules needed:
!   use tao_c_interface_mod
!
! Input:
!   n            -- integer(c_int): line number of out_io_buffer%lines      
!
! Output:
!   c_string_ptr -- type(c_ptr): C pointer to access out_io_buffer%lines(n)
!-

function tao_c_out_io_buffer_get_line(n) bind(c) result(c_string_ptr)
type(c_ptr) :: c_string_ptr
 character(c_char) :: c_line(n_char_show+1) 
integer(c_int), value :: n
integer :: i
i = n
c_line = c_string(out_io_buffer_get_line(i))
c_string_ptr = c_loc(c_line(1)) ! must point to (1) to avoid gfortran compiler error
end function

!------------------------------------------------------------------------
!+ 
! function tao_c_command (c_str) bind(c) result(err)
!
! Function to send a command to Tao from C. 
!
! Modules needed:
!   use tao_c_interface_mod
!
! Input:
!   c_str(*) -- character(c_char): Tao command string   
!
! Output:
!   err      -- integer(c_int): Error code from tao_top_level (0 => no errors)
!-

function tao_c_command (c_str) bind(c) result(err)
character(len=1,kind=c_char), dimension(*), intent(in) ::  c_str
character(200) :: f_str
integer :: errcode
integer(c_int) :: err
call to_f_str (c_str, f_str)
call tao_top_level(command = trim(f_str), errcode = errcode)
err = errcode
end function

end module tao_c_interface_mod
