!+
! Module tao_c_interface_mod
!
! Module containing a set of routines primarily used to communicate with Python via ctypes.
! The python side code is in:
!   tao/python/pytao/tao_ctypes/core.py
!-


module tao_c_interface_mod

use tao_interface
use fortran_cpp_utils

implicit none

type tao_c_interface_common_struct
  real(c_double), allocatable :: c_real(:)
  integer(c_int), allocatable :: c_integer(:)
  character(c_char) :: c_line(n_char_show+1) 
  integer(c_int) :: n_real = 0, n_int = 0
end type

type (tao_c_interface_common_struct), target, save :: tao_c_interface_com

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_init_tao (c_str) bind(c) result (err)
!
! Function to initialize Tao.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Input:
!   c_str(*) -- character(c_char): Tao startup string.
!
! Output:
!   err      -- integer(c_int): Error code from tao_top_level (0 => no errors)
!-

function tao_c_init_tao (c_str) bind(c) result (err)
character(c_char), intent (in) ::  c_str(*)
character(200) :: f_str
integer :: errcode
integer(c_int) :: err

! For keep terminal printing on for debugging.

call out_io_print_and_capture_setup (print_on = .false., capture_state = 'BUFFERED', capture_add_null = .true.)
call to_f_str (c_str, f_str)
call tao_top_level(command = trim(f_str), errcode = errcode)
err = errcode

end function tao_c_init_tao

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_command (c_str) bind(c) result (err)
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

function tao_c_command (c_str) bind(c) result (err)
character(c_char), intent (in) ::  c_str(*)
character(200) :: f_str
integer :: errcode
integer(c_int) :: err

!

call to_f_str (c_str, f_str)
call tao_top_level(command = trim(f_str), errcode = errcode)
err = errcode

end function tao_c_command

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_out_io_buffer_num_lines() bind(c) result (n_lines)
!
! Function to access string buffer out_io_buffer%n_lines number of lines from C.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   n_lines -- integer(c_int): C access to out_io_buffer%n_lines
!-

function tao_c_out_io_buffer_num_lines() bind(c) result (n_lines)
integer(c_int) ::  n_lines

!

n_lines = out_io_buffer_num_lines()

end function tao_c_out_io_buffer_num_lines

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_real_array_size() bind(c) result (n_size)
!
! Function to access tao_c_interface_com%n_real number of items in the real array from C.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   n_size -- integer(c_int): C access to tao_c_interface_com%n_real.
!-

function tao_c_real_array_size() bind(c) result (n_size)
integer(c_int) ::  n_size

!

n_size = tao_c_interface_com%n_real

end function tao_c_real_array_size

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_integer_array_size() bind(c) result (n_size)
!
! Function to access tao_c_interface_com%n_int number of items in the integer array from C.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   n_size -- integer(c_int): C access to out_int_array%n_size
!-

function tao_c_integer_array_size() bind(c) result (n_size)
integer(c_int) ::  n_size

!

n_size = tao_c_interface_com%n_int

end function tao_c_integer_array_size

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_out_io_buffer_get_line(n) bind(c) result (c_string_ptr)
!
! Function to access string buffer out_io_buffer%line(n) from C.
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

function tao_c_out_io_buffer_get_line(n) bind(c) result (c_string_ptr)
type(c_ptr) :: c_string_ptr
integer(c_int), value :: n
integer :: i

!

i = n
call to_c_str (out_io_buffer_get_line(i), tao_c_interface_com%c_line)
c_string_ptr = c_loc(tao_c_interface_com%c_line(1)) ! must point to (1) to avoid gfortran compiler error

end function tao_c_out_io_buffer_get_line

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_get_real_array() bind(c) result (array_ptr)
!
! Function to get the buffered real array from C.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   array_ptr -- type(c_ptr): C pointer to the real array.
!-

function tao_c_get_real_array() bind(c) result (array_ptr)
type(c_ptr) :: array_ptr
integer :: i

!

array_ptr = c_loc(tao_c_interface_com%c_real(1))

end function tao_c_get_real_array

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_get_integer_array() bind(c) result (array_ptr)
!
! Function to get the buffered integer array from C.
!
! Modules needed:
!   use tao_c_interface_mod
!
! Output:
!   array_ptr -- type(c_ptr): C pointer to the real array.
!-

function tao_c_get_integer_array() bind(c) result (array_ptr)
type(c_ptr) :: array_ptr
integer :: i

!

array_ptr = c_loc(tao_c_interface_com%c_integer(1))

end function tao_c_get_integer_array

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine tao_c_out_io_buffer_reset() bind(c)
!
! Routine to reset the buffer.
!-

subroutine tao_c_out_io_buffer_reset() bind(c)

call out_io_buffer_reset()
if (allocated(tao_c_interface_com%c_real))  deallocate (tao_c_interface_com%c_real)
if (allocated(tao_c_interface_com%c_integer)) deallocate (tao_c_interface_com%c_integer)
tao_c_interface_com%n_real = 0
tao_c_interface_com%n_int = 0

end subroutine tao_c_out_io_buffer_reset

end module tao_c_interface_mod
