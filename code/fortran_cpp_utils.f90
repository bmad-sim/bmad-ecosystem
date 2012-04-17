module fortran_and_cpp

use precision_def

type c_dummy_struct
  real(rp) dummy
end type

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine remove_null_in_string 
! 
! This is an overloaded routine for:
!   remove_null_in_string_char (str_char, str_out)
!   remove_null_in_string_arr (str_arr, str_out)
!
! Routine to convert a null character in a string to a blank.
! All characters thereafter are similarly converted.
! This is useful for converting a C style string to Fortran.
! If there is no null character then str_out = str_in.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   str_char   -- Character(*): Input string with null character.
!   str_arr(*) -- Character(1): Input array of null terminated character(1) characters.
!
! Output:
!   str_out -- Character(*): String with null character converted.
!-

interface remove_null_in_string
  module procedure remove_null_in_string_char
  module procedure remove_null_in_string_arr
end interface

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic (logic) result (c_log)
!
! Function to convert from a fortran logical to a C logical.
!
! Modules needed:
!   use fortran_and_cpp
!
! Input:
!   logic -- Logical: Fortran logical.
!
! Output:
!   c_log -- Integer: C logical.
!-

function c_logic (logic) result (c_log)

implicit none

logical logic
integer c_log

!

if (logic) then
  c_log = 1
else
  c_log = 0
endif

end function c_logic

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
!
! Modules needed:
!   use fortran_and_cpp
!
! Input:
!   logic -- Integer: C logical.
!
! Output:
!   f_log -- Logical: Fortran logical.
!-

function f_logic (logic) result (f_log)

implicit none

integer logic
logical f_log

!

if (logic == 0) then
  f_log = .false.
else
  f_log = .true.
endif

end function f_logic

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function r_size (ptr) result (this_size)
!
! Function to return the size of a real pointer.
! If the pointer is not associated then 0 is returned.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   ptr(:) -- Real(rp), pointer: Pointer to an array.
!
! Output:
!   this_size -- Integer: Size of array. 0 if not associated.
!-

function r_size (ptr) result (this_size)

implicit none

real(rp), pointer :: ptr(:)
integer this_size

this_size = 0
if (associated(ptr)) this_size = size(ptr)

end function r_size

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function i_size (ptr) result (this_size)
!
! Function to return the size of an integer pointer.
! If the pointer is not associated then 0 is returned.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   ptr(:) -- Integer, pointer: Pointer to an array.
!
! Output:
!   this_size -- Integer: Size of array. 0 if not associated.
!-


function i_size (ptr) result (this_size)

implicit none

integer, pointer :: ptr(:)
integer this_size

this_size = 0
if (associated(ptr)) this_size = size(ptr)

end function i_size

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine remove_null_in_string_arr (str_in, str_out)
! 
! This routine overloaded by:
!        remove_null_in_string
! See remove_null_in_string for more details.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   str_in(*) -- Character(1): Input character array. Null terminated.
!
! Output:
!   str_out -- Character(*): String with null character converted.
!-

subroutine remove_null_in_string_arr (str_in, str_out)

implicit none

character(1) str_in(*)
character(*) str_out
integer ix

!

str_out = ''
do ix = 1, 32000
  if (str_in(ix) == char(0)) return
  str_out(ix:ix) = str_in(ix)
enddo

end subroutine remove_null_in_string_arr

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine remove_null_in_string_char (str_in, str_out)
! 
! This routine overloaded by:
!        remove_null_in_string
! See remove_null_in_string for more details.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   str_in -- Character(*): Input string with null character.
!
! Output:
!   str_out -- Character(*): String with null character converted.
!-

subroutine remove_null_in_string_char (str_in, str_out)

implicit none

character(*) str_in, str_out
integer ix

!

ix = index(str_in, char(0))
if (ix == 0) then
  str_out = str_in
else
  str_out = str_in(1:ix-1)
endif

end subroutine remove_null_in_string_char

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_str (str) result (c_string)
!
! Function to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   str   -- Character(*): Input character string
!
! Output:
!   c_str -- Character(*): String with a null put just after the last
!             non-blank character.
!-

function c_str (str) result (c_string)

implicit none

character(*) str
character(len_trim(str)+1) c_string

c_string = trim(str) // char(0)

end function c_str

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function mat2arr (mat) result (arr)
!
! Function to take a matrix and turn it into an array:
!   arr(n2*(i-1) + j) = mat(i,j)
! where n2 = size(mat,2).
! This is used for passing matrices to C++ routines.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   mat(:,:)  -- Real(rp): Input matrix
!
! Output:
!   arr(:)   -- Real(rp): Output array 
!-

function mat2arr (mat) result (arr)

implicit none

real(rp) mat(:,:)
real(rp) arr(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) arr(n2*(i-1) + j) = mat(i,j)
 
end function mat2arr

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tensor2arr (tensor) result (arr)
!
! Function to take a tensorrix and turn it into an array:
!   arr(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! where n2 = size(tensor,2).
! This is used for passing tensorrices to C++ routines.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   tensor(:,:)  -- Real(rp): Input tensorrix
!
! Output:
!   arr(:)   -- Real(rp): Output array 
!-

function tensor2arr (tensor) result (arr)

implicit none

real(rp) tensor(:,:,:)
real(rp) arr(size(tensor))
integer i, j, k, n1, n2, n3

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) arr(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function tensor2arr

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function arr2mat (arr, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = arr(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   arr(:)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- Real(rp): Output matrix
!-

function arr2mat (arr, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
real(rp) arr(:)
real(rp) mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = arr(n2*(i-1) + j) 
 
end function arr2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function arr2tensor (arr, n1, n2, n3) result (tensor)
!
! Function to take a an array and turn it into a tensor:
!   tensor(i,j) = arr(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Modules needed:
!  use fortran_and_cpp
!
! Input:
!   arr(:)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- Real(rp): Output tensor.
!-

function arr2tensor (arr, n1, n2, n3) result (tensor)

implicit none

integer i, j, k, n1, n2, n3
real(rp) arr(:)
real(rp) tensor(n1,n2,n3)

forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = arr(n3*n2*(i-1) + n3*(j-1) + k) 
 
end function arr2tensor

end module

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine set_string (string_in, n_in, string_out)
!
! This routine is called by a C/C++ routine to set:
!   string_out = string_in
! Where:
!   string_in is a string generated on the C/C++ side.
!   string_out is a string originally defined in a fortran routine.
!
! Input:
!   string_in -- Character(n_in): Input string.
!   n_in      -- Integer: Length of string_in.
!   n_out     -- Integer: Length of string_out.
!
! Output:
!   string_out -- Character(n_out): Set equal to string_in.
!-

subroutine set_string (string_in, n_in, string_out, n_out)

implicit none

integer n_in, n_out
character(n_in) string_in
character(n_out) string_out

!

string_out = string_in

end subroutine set_string
