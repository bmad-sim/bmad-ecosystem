module fortran_cpp_utils

use precision_def
use, intrinsic :: iso_c_binding

type c_dummy_struct
  real(rp) dummy
end type

interface
  pure subroutine bool_to_int (logic, int_logic) bind(c)
    import c_bool, c_int
    logical(c_bool), intent(in) :: logic
    integer(c_int), intent(inout) :: int_logic
  end subroutine
end interface

interface
  pure subroutine int_to_bool (int_logic, logic) bind(c)
    import c_bool, c_int
    integer(c_int), intent(in) :: int_logic
    logical(c_bool), intent(inout) :: logic
  end subroutine
end interface

interface vec2fvec
  module procedure bool_vec2fvec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function fscalar2scalar (f_scalar, n) result (c_scalar)
!
! Overloaded function to translate a scalar from Fortran form to C form.
! Overloads:
!   bool_fscalar2scalar (bool_f_scalar, n) result (bool_c_scalar)
!
! Input:
!   bool_f_scalar  -- Logical: Input scalar.
!   n              -- Integer: 0 if actual scalar is not allocated, 1 otherwise.
!
! Output:
!   bool_c_scalar  -- Logical(c_bool): Output scalar
!-

interface fscalar2scalar
  module procedure bool_fscalar2scalar
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function fvec2vec (f_vec, n) result (c_vec)
!
! Overloaded function to translate a vector from Fortran form to C form.
! Overloads:
!   bool_fvec2vec  (bool_f_vec, n)  result (bool_c_vec)
!   int_fvec2vec   (int_f_vec, n)   result (int_c_vec)
!   cmplx_fvec2vec (cmplx_f_vec, n) result (cmplx_c_vec)
!   real_fvec2vec  (real_f_vec, n)  result (real_c_vec)
!
! Input:
!   bool_f_vec(:)  -- logical: Input vector
!   int_f_vec(:)   -- integer: Input vector
!   cmplx_f_vec(:) -- complex: Input vector
!   real_f_vec(:)  -- real(rp): Input vector
!
! Output:
!   bool_c_vec(:)  -- logical(c_bool): Output array 
!   int_c_vec(:)   -- integer(c_int): Output array 
!   cmplx_c_vec(:) -- complex(c_double_complex): Output array 
!   real_c_vec(:)  -- real(c_double): Output array 
!-

interface fvec2vec
  module procedure real_fvec2vec
  module procedure int_fvec2vec
  module procedure cmplx_fvec2vec
  module procedure bool_fvec2vec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function mat2vec (mat, n) result (vec)
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
!   cmplx_mat2vec  (cmplx_mat)  result (cmplx_vec)
!   bool_mat2vec   (bool_mat)   result (bool_vec)
!
! Input:
!   real_mat(:,:)   -- Real(rp): Input matrix
!   int_mat(:,:)    -- Integer: Input matrix
!   cmplx_mat(:,:)  -- Complex(rp): Input matrix
!   bool_mat(:,:)   -- Logical: Input matrix
!   n               -- Integer: Number of elements. Normally this is size(mat). 
!                        Set to 0 if actual mat arg is not allocated.
!
! Output:
!   real_vec(*)   -- Real(c_double): Output array 
!   int_vec(*)    -- Integer(c_int): Output array 
!   cmplx_vec(*) -- complex(c_double_complex): Output array 
!   bool_vec(*)   -- Logical(c_bool): Output array 
!-

interface mat2vec
  module procedure real_mat2vec
  module procedure int_mat2vec
  module procedure cmplx_mat2vec
  module procedure bool_mat2vec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tensor2vec (tensor, n) result (vec)
!
! Function to take a tensor and turn it into an array in 
! C standard row-major order::
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j,k)
! where n2 = size(tensor,2).
! This is used for passing tensorrices to C++ routines.
!
! Input:
!   tensor(:,:)     -- Real(rp): Input tensorrix
!   n               -- Integer: Number of elements. Normally this 
!                        is size(tensor), 0 if actual tensor arg is not allocated.
!
! Output:
!   vec(*)   -- Real(c_double): Output array 
!-

interface tensor2vec
  module procedure real_tensor2vec
  module procedure int_tensor2vec
  module procedure cmplx_tensor2vec
  module procedure bool_tensor2vec
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine vec2mat (vec, mat)
!
! Overloaded routine to take a an array in C standard row-major 
! order and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Overloaded functions:
!   real_vec2mat
!   int_vec2mat
!   cmplx_vec2mat
!   bool_vec2mat
!
! Input:
!   vec(*)   -- Real(c_double): Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(:,:)  -- Real(rp): Output matrix
!-

interface vec2mat
  module procedure real_vec2mat
  module procedure int_vec2mat
  module procedure cmplx_vec2mat
  module procedure bool_vec2mat
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine vec2tensor (vec, tensor)
!
! Routine to take a an array in C standard row-major 
! order and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Input:
!   vec(*)   -- Real(c_double): Input array.
!   n1       -- Integer: Size of first tensor index.
!   n2       -- Integer: Size of second tensor index.
!   n3       -- Integer: Size of third tensor index.
!
! Output:
!   tensor(n1,n2,n3)  -- Real(rp): Output tensor.
!-

interface vec2tensor
  module procedure real_vec2tensor
  module procedure int_vec2tensor
  module procedure cmplx_vec2tensor
  module procedure bool_vec2tensor
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
! Input:
!   int_logic  -- Integer: C logical.
!   bool_logic -- Logical(c_bool): C logical.
!
! Output:
!   f_log -- Logical: Fortran logical.
!-

interface f_logic
  module procedure f_logic_bool
  module procedure f_logic_int
end interface

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic
!
! Function to convert from a C logical to a Fortran logical.
! This function overloads:
!   c_logic1  (logic) result (c_logic)
!   c_logic_vec (logic_vec) result (c_logic_vec)
!
! Input:
!   logic        -- logical: Fortran logical.
!   logic_vec(:) -- logical: Fortran logical vector.
!
! Output:
!   c_log          -- logical(c_bool): C boolien. 0 if False, 1 if True.
!   c_logic_vec(:) -- logical(c_bool): C boolien vector
!-

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
! Input:
!   logic -- Logical: Fortran logical.
!
! Output:
!   c_log -- Integer: C logical.
!-

pure function c_logic1 (logic) result (c_log)

implicit none

logical, intent(in) :: logic
logical(c_bool) c_log

!

if (logic) then
  call int_to_bool (1, c_log)
else
  call int_to_bool (0, c_log)
endif

end function c_logic1

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic_vec (logic_vec) result (c_bool_vec)
!
! Function to convert from a fortran logical to a C logical.
! See c_logic for more details.
!
! Input:
!   logic_vec(:) -- logical: Fortran logical vector.
!
! Output:
!   c_bool_vec(:) -- logical(c_bool): C logical vector.
!-

pure function c_logic_vec (logic_vec) result (c_bool_vec)

implicit none

logical, intent(in) :: logic_vec(:)
logical(c_bool) c_bool_vec(size(logic_vec))
integer i

!

do i = 1, size(logic_vec)
  c_bool_vec(i) = c_logic1(logic_vec(i))
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

integer, intent(in) :: logic
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
! Function f_logic_bool (logic) result (f_log)
!
! Function to convert from a C logical to a Fortran logical.
! This function is overloaded by f_logic.
! See f_logic for more details.
!-

function f_logic_bool (logic) result (f_log)

implicit none

logical(c_bool), intent(in) :: logic
logical f_log
integer int_logic

!

call bool_to_int (logic, int_logic)

if (int_logic == 0) then
  f_log = .false.
else
  f_log = .true.
endif

end function f_logic_bool

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function r_size (ptr) result (this_size)
!
! Function to return the size of a real pointer.
! If the pointer is not associated then 0 is returned.
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
  if (str_in(ix) == c_null_char) return
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

ix = index(str_in, c_null_char)
if (ix == 0) then
  str_out = str_in
else
  str_out = str_in(1:ix-1)
endif

end subroutine remove_null_in_string_char

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine to_c_str (f_str, c_str)
!
! Subroutine to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Input:
!   f_str   -- Character(*): Input character string
!
! Output:
!   c_str(*) -- Character(kind=c_char): String with a null put just after the last
!                    non-blank character.
!-

pure subroutine to_c_str (f_str, c_str)

implicit none

character(*), intent(in) :: f_str
character(kind=c_char), intent(out) ::  c_str(*)
integer i

!

do i = 1, len_trim(f_str)
  c_str(i) = f_str(i:i)
enddo

c_str(i) = c_null_char

end subroutine to_c_str

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! function c_string(f_str) result (c_str)
!
! Functional form of subroutine to_c_str (f_str, c_str)
!
! IMPRTANT NOTE: This routine is meant to be use "in place" like:
!   call open_dir (c_string(my_dir), ...)
! Do not use with the equal operator:
!   my_c_str = c_string(my_f_str)   ! DO NOT USE THIS SYNTAX !!!
! The problem is that the output of c_string is an array and my_c_str is an array so the above set will
! fail if the array size of my_c_str and the array size of the output of c_string are not the same.
!
! Input:
!   f_str   -- Character(*): Input character string
!
! Output:
!   c_str(len_trim(f_str)+1) 
!           -- Character(kind=c_char): String with a null at the end.
!-

pure function c_string(f_str) result (c_str)

implicit none
character(*), intent(in) :: f_str
character(kind=c_char), dimension(len_trim(f_str)+1) :: c_str

call to_c_str (f_str, c_str)

end function c_string

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine to_f_str (c_str, f_str)
!
! Subroutine to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Input:
!   c_str(*) -- Character(kind=c_char): C-style string.
!
! Output:
!   f_str -- Character(*): Output character string.
!-

subroutine to_f_str (c_str, f_str)

implicit none

character(*) f_str
character(kind=c_char) c_str(*)
integer i

!

do i = 1, len(f_str)
  if (c_str(i) == c_null_char) then
    f_str(i:) = ''
    return
  endif
  f_str(i:i) = c_str(i)
enddo

end subroutine to_f_str


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_fscalar2scalar (f_scalar, n) result (c_scalar)
!
! Function transform from Fortran to C.
! See fscalar2scalar for more details
!-

function bool_fscalar2scalar (f_scalar, n) result (c_scalar)

implicit none

integer n
logical f_scalar
logical(c_bool) c_scalar


!

if (n == 0) then
  call int_to_bool (0, c_scalar)
else
  c_scalar = c_logic(f_scalar)
endif

end function bool_fscalar2scalar

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_fvec2vec (f_vec, n) result (c_vec)
!
! Function transform from Fortran to C.
! See fvec2vec for more details
!-

function real_fvec2vec (f_vec, n) result (c_vec)

implicit none

integer n, i
real(rp) f_vec(:)
real(c_double), target :: c_vec(n)

!

forall (i = 1:n) c_vec(i) = f_vec(i)

end function real_fvec2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function int_fvec2vec (f_vec, n) result (c_vec)
!
! Function transform from Fortran to C.
! See fvec2vec for more details
!-

function int_fvec2vec (f_vec, n) result (c_vec)

implicit none

integer n, i
integer f_vec(:)
integer(c_int), target :: c_vec(n)

forall (i = 1:n) c_vec(i) = f_vec(i)
 
end function int_fvec2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_fvec2vec (f_vec, n) result (c_vec)
!
! Function transform from Fortran to C.
! See fvec2vec for more details
!-

function cmplx_fvec2vec (f_vec, n) result (c_vec)

implicit none

integer n, i
complex(rp) f_vec(:)
complex(c_double_complex), target :: c_vec(n)

forall (i = 1:n) c_vec(i) = f_vec(i)
 
end function cmplx_fvec2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_fvec2vec (f_vec, n) result (c_vec)
!
! Function transform from Fortran to C.
! See fvec2vec for more details
!-

function bool_fvec2vec (f_vec, n) result (c_vec)

implicit none

integer n, i
logical f_vec(:)
logical(c_bool), target :: c_vec(n)

forall (i = 1:n) c_vec(i) = c_logic(f_vec(i))
 
end function bool_fvec2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_mat2vec (mat, n) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- Real(rp): Input matrix
!
! Output:
!   vec(:)   -- Real(c_double): Output array 
!-

function real_mat2vec (mat, n) result (vec)

implicit none

integer n
real(rp) mat(:,:)
real(c_double) vec(n)
integer i, j, n1, n2

if (n == 0) return ! Real arg not allocated
n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function real_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function int_mat2vec (mat, n) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- integer: Input matrix
!
! Output:
!   vec(:)   -- integer(c_int): Output array 
!-

function int_mat2vec (mat, n) result (vec)

implicit none

integer n
integer mat(:,:)
integer(c_int) vec(n)
integer i, j, n1, n2

if (n == 0) return ! Real arg not allocated
n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function int_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_mat2vec (mat, n) result (vec)
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

function bool_mat2vec (mat, n) result (vec)

implicit none

integer n
logical mat(:,:)
logical(c_bool) vec(n)
integer i, j, n1, n2

if (n == 0) return ! Real arg not allocated

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
! Function cmplx_mat2vec (mat, n) result (vec)
!
! Function to take a matrix and turn it into an array:
!   vec(n2*(i-1) + j) = mat(i,j)
! See mat2vec for more details
!
! Input:
!   mat(:,:)  -- complex(rp): Input matrix
!
! Output:
!   vec(:)   -- complex(c_double_complex): Output array 
!-

function cmplx_mat2vec (mat, n) result (vec)

implicit none

integer n
complex(rp) mat(:,:)
complex(c_double_complex) vec(n)
integer i, j, n1, n2

if (n == 0) return ! Real arg not allocated
n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) vec(n2*(i-1) + j) = mat(i,j)
 
end function cmplx_mat2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function real_tensor2vec (tensor, n) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:,:)  -- Real(rp): Input tensorrix
!
! Output:
!   vec(:)   -- Real(c_double): Output array 
!-

function real_tensor2vec (tensor, n) result (vec)

implicit none

integer n
real(rp) tensor(:,:,:)
real(c_double), target :: vec(n)
integer i, j, k, n1, n2, n3

if (n == 0) return ! Real arg not allocated
n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function real_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function int_tensor2vec (tensor, n) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:,:)  -- Integer: Input tensorrix
!
! Output:
!   vec(:)   -- Integer(c_int): Output array 
!-

function int_tensor2vec (tensor, n) result (vec)

implicit none

integer n
integer tensor(:,:,:)
integer(c_int), target :: vec(n)
integer i, j, k, n1, n2, n3

if (n == 0) return ! Real arg not allocated
n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function int_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function cmplx_tensor2vec (tensor, n) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:,:)  -- complex(rp): Input tensorrix
!
! Output:
!   vec(:)   -- complex(c_double_complex): Output array 
!-

function cmplx_tensor2vec (tensor, n) result (vec)

implicit none

integer n
complex(rp) tensor(:,:,:)
complex(c_double_complex), target :: vec(n)
integer i, j, k, n1, n2, n3

if (n == 0) return ! Real arg not allocated
n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = tensor(i,j,k)
 
end function cmplx_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function bool_tensor2vec (tensor, n) result (vec)
!
! Function to take a tensor and turn it into an array:
!   vec(n3*n2*(i-1) + n3*(j - 1) + k) = tensor(i,j, k)
! See tensor2vec for more details
!
! Input:
!   tensor(:,:,:)  -- logical: Input tensorrix
!
! Output:
!   vec(:)   -- logical(c_bool): Output array 
!-

function bool_tensor2vec (tensor, n) result (vec)

implicit none

integer n
logical tensor(:,:,:)
logical(c_bool), target :: vec(n)
integer i, j, k, n1, n2, n3

if (n == 0) return ! Real arg not allocated
n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor, 3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) vec(n3*n2*(i-1) + n3*(j-1) + k) = c_logic(tensor(i,j,k))
 
end function bool_tensor2vec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bool_vec2fvec (c_vec, f_vec)
!
!
! Input:
!   c_vec(*)   -- Logical(c_bool): Input array.
!
! Output:
!   f_vec(n1,n2)  -- Logical: Output f_vec
!-

subroutine bool_vec2fvec (c_vec, f_vec)

implicit none

integer i
logical(c_bool) c_vec(*)
logical f_vec(:)


do i = 1, size(f_vec)
  f_vec(i) = f_logic(c_vec(i))
enddo
 
end subroutine bool_vec2fvec

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine real_vec2mat (vec, mat)
!
! Subroutine to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Input:
!   vec(*)   -- Real(c_double): Input array.
!
! Output:
!   mat(n1,n2)  -- Real(rp): Output matrix
!-

subroutine real_vec2mat (vec, mat)

implicit none

integer i, j, n1, n2
real(c_double) vec(*)
real(rp) mat(:,:)

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end subroutine real_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine int_vec2mat (vec, mat)
!
! Subroutine to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Input:
!   vec(*)   -- integer: Input array.
!
! Output:
!   mat(:,:)  -- integer: Output matrix
!-

subroutine int_vec2mat (vec, mat)

implicit none

integer i, j, n1, n2
integer(c_int) vec(*)
integer mat(:,:)

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end subroutine int_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bool_vec2mat (vec, mat)
!
! Subroutine to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Input:
!   vec(*)   -- logical: Input array.
!
! Output:
!   mat(:,:)  -- logical: Output matrix
!-

subroutine bool_vec2mat (vec, mat)

implicit none

integer i, j, n1, n2
logical(c_bool) vec(*)
logical mat(:,:)

n1 = size(mat, 1); n2 = size(mat, 2)
do i = 1,n1
do j = 1,n2
  mat(i,j) = f_logic(vec(n2*(i-1) + j))
enddo
enddo

end subroutine bool_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine cmplx_vec2mat (vec, mat)
!
! Subroutine to take a an array and turn it into a matrix:
!   mat(i,j) = vec(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Input:
!   vec(*)   -- complex(c_double_complex): Input array.
!
! Output:
!   mat(:,:)  -- complex(rp): Output matrix
!-

subroutine cmplx_vec2mat (vec, mat)

implicit none

integer i, j, n1, n2
complex(c_double_complex) vec(*)
complex(rp) mat(:,:)

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) mat(i,j) = vec(n2*(i-1) + j) 
 
end subroutine cmplx_vec2mat

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine real_vec2tensor (vec, tensor)
!
! Subroutine to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Input:
!   vec(*)   -- Real(rp): Input array.
!
! Output:
!   tensor(:,:,:)  -- Real(rp): Output tensor.
!-

subroutine real_vec2tensor (vec, tensor)

implicit none

integer i, j, k, n1, n2, n3
real(c_double) vec(*)
real(rp) tensor(:,:,:)

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor,3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end subroutine real_vec2tensor

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine int_vec2tensor (vec, tensor)
!
! Subroutine to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Input:
!   vec(*)   -- integer: Input array.
!
! Output:
!   tensor(:,:,:)  -- integer(c_int): Output tensor.
!-

subroutine int_vec2tensor (vec, tensor)

implicit none

integer i, j, k, n1, n2, n3
integer(c_int) vec(*)
integer tensor(:,:,:)

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor,3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end subroutine int_vec2tensor

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine cmplx_vec2tensor (vec, tensor)
!
! Subroutine to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Input:
!   vec(*)   -- complex(c_double_complex): Input array.
!
! Output:
!   tensor(:,:,:)  -- complex(rp): Output tensor.
!-

subroutine cmplx_vec2tensor (vec, tensor)

implicit none

integer i, j, k, n1, n2, n3
complex(c_double_complex) vec(*)
complex(rp) tensor(:,:,:)

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor,3)
forall (i = 1:n1, j = 1:n2, k = 1:n3) tensor(i,j,k) = vec(n3*n2*(i-1) + n3*(j-1) + k) 
 
end subroutine cmplx_vec2tensor

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bool_vec2tensor (vec, tensor)
!
! Subroutine to take a an array and turn it into a tensor:
!   tensor(i,j) = vec(n3*n2*(i-1) + n3*j + k) 
! This is used for getting tensorrices from C++ routines.
!
! Input:
!   vec(*)   -- logical(c_bool): Input array.
!
! Output:
!   tensor(:,:,:)  -- logical: Output tensor.
!-

subroutine bool_vec2tensor (vec, tensor)

implicit none

integer i, j, k, n1, n2, n3
logical(c_bool) vec(*)
logical tensor(:,:,:)

n1 = size(tensor, 1); n2 = size(tensor, 2); n3 = size(tensor,3)
do i = 1, n1
do j = 1, n2
do k = 1, n3
  tensor(i,j,k) = f_logic(vec(n3*n2*(i-1) + n3*(j-1) + k))
enddo
enddo
enddo

end subroutine bool_vec2tensor

!-----------------------------------------------------------------------------

end module


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine set_fortran_string (string_in, n_in, string_out)
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

subroutine set_fortran_string (string_in, n_in, string_out, n_out)

implicit none

integer n_in, n_out
character(n_in) string_in
character(n_out) string_out

!

string_out = string_in

end subroutine set_fortran_string
