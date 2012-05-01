module fortran_cpp_utils

use precision_def
use, intrinsic :: iso_c_binding

type c_dummy_struct
  real(rp) dummy
end type

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function mat2vec (mat) result (vec)
!
! Overloaded function to take a matrix and turn it into an array in 
! C standard row-major order:
!   vec(n2*(i-1) + j) = mat(i,j)
! where n2 = size(mat,2).
! This is used for passing matrices to C++ routines.
!
! Overloaded functions:
!   real_mat2vec   (real_mat)   result (real_vec)
!   int_mat2vec    (int_mat)    result (int_vec)
!   complx_mat2vec (complx_mat) result (complx_vec)
!   bool_mat2vec   (bool_mat)   result (bool_vec)
!   ptr_mat2vec (ptr_mat) result (ptr_vec)
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   real_mat(:,:)   -- Real(rp): Input matrix
!   int_mat(:,:)    -- Integer: Input matrix
!   complx_mat(:,:) -- Çomplex(rp): Input matrix
!   bool_mat(:,:)   -- Logical: Input matrix
!   ptr_mat(:,:)    -- type(c_ptr): Input matrix
!
! Output:
!   real_vec(:)   -- Real(rp): Output array 
!   int_vec(:)    -- Integer: Output array 
!   complx_vec(:) -- Complex(rp): Output array 
!   bool_vec(:)   -- Logical(c_bool): Output array 
!   ptr_vec(:)    -- type(c_ptr): Output array 
!-

interface mat2vec
  module procedure real_mat2vec
  module procedure int_mat2vec
  module procedure cmplx_mat2vec
  module procedure bool_mat2vec
  module procedure ptr_mat2vec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tensor2arr (tensor) result (arr)
!
! Function to take a tensor and turn it into an array in 
! C standard row-major order::
!   arr(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j,k)
! where n2 = size(tensor,2).
! This is used for passing tensorrices to C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   tensor(:,:)  -- Real(rp): Input tensorrix
!
! Output:
!   arr(:)   -- Real(rp): Output array 
!-

interface tensor2vec
  module procedure real_tensor2vec
  module procedure cmplx_tensor2vec
  module procedure ptr_tensor2vec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array in C standard row-major 
! order and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- Real(rp): Output matrix
!-

interface vec2mat
  module procedure real_vec2mat
  module procedure int_vec2mat
  module procedure cmplx_vec2mat
  module procedure bool_vec2mat
  module procedure ptr_vec2mat
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function vec2tensor (vec, n1, n2, n3) result (tensor)
!
! Function to take a an array in C standard row-major 
! order and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- Real(rp): Output tensor.
!-

interface vec2tensor
  module procedure real_vec2tensor
  module procedure cmplx_vec2tensor
  module procedure ptr_vec2tensor
end interface

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
!  use fortran_cpp_utils
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

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
! This function overloads:
!   f_logic_int  (int_logic) result (f_log)
!   f_logic_bool (bool_logic) result (f_log)
!
! Modules needed:
!   use fortran_cpp_utils
!
! Input:
!   int_logic  -- Integer: C logical.
!   bool_logic -- Logical(c_bool): C logical.
!
! Output:
!   f_log -- Logical: Fortran logical.
!-

interface f_logic
  module procedure f_logic_bool1
  module procedure f_logic_bool_vec
  module procedure f_logic_int
end interface

interface c_logic
  module procedure c_logic1
  module procedure c_logic_vec
end interface

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic1 (logic) result (c_log)
!
! Function to convert from a fortran logical to a C logical.
! See c_logic for more details.
!
! Modules needed:
!   use fortran_cpp_utils
!
! Input:
!   logic -- Logical: Fortran logical.
!
! Output:
!   c_log -- Integer: C logical.
!-

function c_logic1 (logic) result (c_log)

implicit none

logical logic
logical(c_bool) c_log

!

if (logic) then
  c_log = 1
else
  c_log = 0
endif

end function c_logic1

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic_vec (logic) result (c_log)
!
! Function to convert from a fortran logical to a C logical.
! See c_logic for more details.
!
! Modules needed:
!   use fortran_cpp_utils
!
! Input:
!   logic -- Logical: Fortran logical.
!
! Output:
!   c_log -- Integer: C logical.
!-

function c_logic_vec (logic) result (c_log)

implicit none

logical logic(:)
logical(c_bool) c_log(size(logic))
integer i

!

do i = 1, size(logic)
  c_log(i) = c_logic1(logic(i))
enddo

end function c_logic_vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic_int (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
! This function is overloaded by f_logic.
! See f_logic for more details.
!-

function f_logic_int (logic) result (f_log)

implicit none

integer logic
logical f_log

!

if (logic == 0) then
  f_log = .false.
else
  f_log = .true.
endif

end function f_logic_int

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic_bool1 (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
! This function is overloaded by f_logic.
! See f_logic for more details.
!-

function f_logic_bool1 (logic) result (f_log)

implicit none

logical(c_bool) logic
logical f_log
integer int_logic

interface
  subroutine bool_to_int (logic, int_logic) bind(c)
    import c_bool, c_int
    logical(c_bool) logic
    integer(c_int) int_logic
  end subroutine
end interface

!

call bool_to_int (logic, int_logic)

if (int_logic == 0) then
  f_log = .false.
else
  f_log = .true.
endif

end function f_logic_bool1

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic_bool_vec (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
! This function is overloaded by f_logic.
! See f_logic for more details.
!-

function f_logic_bool_vec (logic) result (f_log)

implicit none

logical(c_bool) logic(:)
logical f_log(size(logic))
integer i

!

do i = 1, size(logic)
  f_log(i) = f_logic_bool1(logic(i))
enddo

end function f_logic_bool_vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function r_size (ptr) result (this_size)
!
! Function to return the size of a real pointer.
! If the pointer is not associated then 0 is returned.
!
! Modules needed:
!  use fortran_cpp_utils
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
!  use fortran_cpp_utils
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
!  use fortran_cpp_utils
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
!  use fortran_cpp_utils
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
! Subroutine to_c_str (f_string, c_string)
!
! Subroutine to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   f_string   -- Character(*): Input character string
!
! Output:
!   c_string(*) -- Character(kind=c_char): String with a null put just after the last
!                    non-blank character.
!-

subroutine to_c_str (f_string, c_string)

implicit none

character(*) f_string
character(kind=c_char) c_string(*)
integer i

!

do i = 1, len_trim(f_string)
  c_string(i) = f_string(i:i)
enddo

c_string(i) = char(0)

end subroutine to_c_str

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine to_f_str (c_string, f_string)
!
! Subroutine to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   c_string(*) -- Character(kind=c_char): C-style string.
!
! Output:
!   f_string -- Character(*): Output character string.
!-

subroutine to_f_str (c_string, f_string)

implicit none

character(*) f_string
character(kind=c_char) c_string(*)
integer i

!

do i = 1, len(f_string)
  if (c_string(i) == char(0)) then
    f_string(i:) = ''
    return
  endif
  f_string(i:i) = c_string(i)
enddo

end subroutine to_f_str

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_mat2vec (mat) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- Real(rp): Input matrix
!
! Output:
!   vec(:)   -- Real(rp): Output array 
!-

function real_mat2vec (mat) result (vec)

implicit none

real(rp) mat(:,:)
real(rp) vec(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function real_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function int_mat2vec (mat) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- integer: Input matrix
!
! Output:
!   vec(:)   -- integer: Output array 
!-

function int_mat2vec (mat) result (vec)

implicit none

integer mat(:,:)
integer vec(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function int_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_mat2vec (mat) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- logical: Input matrix
!
! Output:
!   vec(:)   -- logical: Output array 
!-

function bool_mat2vec (mat) result (vec)

implicit none

logical mat(:,:)
logical(c_bool) vec(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
do i = 1, n1
do j = 1, n2 
  vec(n2*(i-1) + j) = c_logic(mat(i,j))
enddo
enddo

end function bool_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function ptr_mat2vec (mat) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- type(c_ptr): Input matrix
!
! Output:
!   vec(:)   -- type(c_ptr): Output array 
!-

function ptr_mat2vec (mat) result (vec)

implicit none

type(c_ptr) mat(:,:)
type(c_ptr) vec(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
do i = 1, n1
do j = 1, n2 
  vec(n2*(i-1) + j) = mat(i,j)
enddo
enddo

end function ptr_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_mat2vec (mat) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- complex(rp): Input matrix
!
! Output:
!   vec(:)   -- complex(rp): Output array 
!-

function cmplx_mat2vec (mat) result (vec)

implicit none

complex(rp) mat(:,:)
complex(rp) vec(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function cmplx_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_tensor2vec (tensor) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:)  -- Real(rp): Input tensorrix
!
! Output:
!   vec(:)   -- Real(rp): Output array 
!-

function real_tensor2vec (tensor) result (vec)

implicit none

real(rp) tensor(:,:,:)
real(rp) vec(size(tensor))
integer i, j, k, n1, n2, n3

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function real_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_tensor2vec (tensor) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:)  -- complex(rp): Input tensorrix
!
! Output:
!   vec(:)   -- complex(rp): Output array 
!-

function cmplx_tensor2vec (tensor) result (vec)

implicit none

complex(rp) tensor(:,:,:)
complex(rp) vec(size(tensor))
integer i, j, k, n1, n2, n3

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function cmplx_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function ptr_tensor2vec (tensor) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:)  -- type(c_ptr): Input tensorrix
!
! Output:
!   vec(:)   -- type(c_ptr): Output array 
!-

function ptr_tensor2vec (tensor) result (vec)

implicit none

type(c_ptr) tensor(:,:,:)
type(c_ptr) vec(size(tensor))
integer i, j, k, n1, n2, n3

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function ptr_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- Real(rp): Output matrix
!-

function real_vec2mat (vec, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
real(rp) vec(*)
real(rp) mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end function real_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function int_vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- integer: Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- integer: Output matrix
!-

function int_vec2mat (vec, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
integer vec(*)
integer mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end function int_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- logical: Input array.
!   n1       -- integer: Size of first mat index.
!   n2       -- integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- logical: Output matrix
!-

function bool_vec2mat (vec, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
logical(c_bool) vec(*)
logical mat(n1,n2)

do i = 1,n1
do j = 1,n2
  mat(i,j) = f_logic(vec(n2*(i-1) + j))
enddo
enddo

end function bool_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- complex(rp): Input array.
!   n1       -- integer: Size of first mat index.
!   n2       -- integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- complex(rp): Output matrix
!-

function cmplx_vec2mat (vec, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
complex(rp) vec(*)
complex(rp) mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end function cmplx_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function ptr_vec2mat (vec, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- type(c_ptr): Input array.
!   n1       -- integer: Size of first mat index.
!   n2       -- integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- type(c_ptr): Output matrix
!-

function ptr_vec2mat (vec, n1, n2) result (mat)

implicit none

integer i, j, n1, n2
type(c_ptr) vec(*)
type(c_ptr) mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end function ptr_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_vec2tensor (vec, n1, n2, n3) result (tensor)
!
! Function to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- Real(rp): Output tensor.
!-

function real_vec2tensor (vec, n1, n2, n3) result (tensor)

implicit none

integer i, j, k, n1, n2, n3
real(rp) vec(*)
real(rp) tensor(n1,n2,n3)

forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end function real_vec2tensor

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_vec2tensor (vec, n1, n2, n3) result (tensor)
!
! Function to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- complex(rp): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- complex(rp): Output tensor.
!-

function cmplx_vec2tensor (vec, n1, n2, n3) result (tensor)

implicit none

integer i, j, k, n1, n2, n3
complex(rp) vec(*)
complex(rp) tensor(n1,n2,n3)

forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end function cmplx_vec2tensor

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function ptr_vec2tensor (vec, n1, n2, n3) result (tensor)
!
! Function to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Modules needed:
!  use fortran_cpp_utils
!
! Input:
!   vec(*)   -- type(c_ptr): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- type(c_ptr): Output tensor.
!-

function ptr_vec2tensor (vec, n1, n2, n3) result (tensor)

implicit none

integer i, j, k, n1, n2, n3
type(c_ptr) vec(*)
type(c_ptr) tensor(n1,n2,n3)

forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end function ptr_vec2tensor

!-----------------------------------------------------------------------------

end module

!-----------------------------------------------------------------------------


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
