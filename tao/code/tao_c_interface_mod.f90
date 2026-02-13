!+
! Module tao_c_interface_mod
!
! Module containing a set of routines to interface with C.
! This is primarily used to communicate with Python via ctypes.
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

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_c_double (re, n, exact, init_val)
!
! Routine to reallocate an array of c_double reals.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   re(:)  -- Real(c_double), Allocatable: Real array.
!   n      -- Integer: Size wanted.
!   exact  -- Logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   re(:)  -- Real(c_double), Allocatable: Allocated array with size(re) >= n.
!-

subroutine re_allocate_c_double (re, n, exact, init_val)

implicit none

real(c_double), allocatable :: re(:), temp_re(:)
real(c_double), optional :: init_val

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (allocated(re)) then
  n_old = size(re)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (re, temp_re)
  allocate (re(n))
  if (present(init_val)) re = init_val
  n_save = min(n, n_old)
  re(1:n_save) = temp_re(1:n_save)
  deallocate (temp_re)  
else
  allocate (re(n))
  if (present(init_val)) re = init_val
endif

end subroutine re_allocate_c_double

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function tao_c_init_tao (c_str) bind(c) result (err)
!
! Function, which can be called from C, to initialize Tao.
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

s%init   = tao_init_struct()
s%initialized = .false.

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

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine tao_c_get_beam_track_active_element() bind(c) result (ix_ele)
!
! Get the active element index being tracked through `tao_beam_track`.
!-

function tao_c_get_beam_track_element() bind(c) result (ix_ele)

integer(c_int) :: ix_ele
ix_ele = s%com%ix_beam_track_active_element

end function tao_c_get_beam_track_element

end module tao_c_interface_mod
