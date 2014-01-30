module tao_c_interface_mod
use iso_c_binding
use tao_mod
use fortran_cpp_utils

implicit none

contains

function tao_c_scratch_n_lines() bind(c) result(n)
integer(c_int) ::  n
n = scratch%n_lines
end function

function tao_c_scratch_line(n) bind(c) result(c_string_ptr)
type(c_ptr) :: c_string_ptr
integer(c_int), value :: n
integer :: i
i = n
scratch%c_line = c_string(scratch%lines(i)) 
c_string_ptr = c_loc(scratch%c_line)

end function

function tao_c_set_init_args(c_str) bind(c) result(error)
use tao_command_mod
character(len=1,kind=c_char), dimension(*), intent(in) ::  c_str
character(2000) :: f_str
character(200) :: cmd_word(10)
integer(c_int) :: error
integer :: n
logical:: err
call to_f_str(c_str, f_str)
call tao_cmd_split (f_str, 10, cmd_word, .false., err)
!count words actually found
!TODO: generalize 10
do n=1, 10
  if (cmd_word(n) == '') exit
enddo
call tao_parse_command_args (err, cmd_word(1:n-1) )
error = merge(1, 0, err)
end function


function tao_c_command(c_str) bind(c) result(err)
character(len=1,kind=c_char), dimension(*), intent(in) ::  c_str
character(200) :: f_str
integer :: errcode
integer(c_int) :: err
call to_f_str (c_str, f_str)
call tao_top_level(command = trim(f_str), errcode = errcode)
err = errcode
end function

end module tao_c_interface_mod