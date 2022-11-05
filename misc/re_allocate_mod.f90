module re_allocate_mod

use utilities_mod
use sim_utils_struct

!+
! Subroutine re_allocate (...)
!
! Subroutine to reallocate a 1-dim allocatable array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
! 
! Note: For pointers to an array use the re_associate routine.
! Note: Also see the re_allocate2 and re_allocate2d routines.
! 
! Note: using exact = False can increase computation speed by
! preventing unneccessary deallocations/reallocations.
!
! Note: re_allocate is an overloaded name for: 
!   Subroutine re_allocate_string (str, n, exact, init_val)
!   Subroutine re_allocate_var_string1 (vstr, n, exact, init_val)
!   Subroutine re_allocate_var_string (var_str, n, exact, init_val)
!   Subroutine re_allocate_integer (inte, n, exact, init_val)
!   Subroutine re_allocate_real (re, n, exact, init_val)
!   Subroutine re_allocate_all_pointer (a_ptr, n, exact)
!   Subroutine re_allocate_complex (cmpl, n, exact, init_val)
!   Subroutine re_allocate_logical (logic, n, exact, init_val)
!
! Input:
!   str(:)      -- character(*), allocatable: String array.
!   vstr        -- character(:), allocatable: Variable length string
!   var_str(:)  -- var_length_string_struct, allocatable: Variable length string structure.
!   inte(:)     -- integer, allocatable: Integer array.
!   re(:)       -- real(rp), Allocatable: Real array.
!   a_ptr(:)    -- all_pointer_struct: array of all_pointer_structs.
!   cmpl(:)     -- complex(rp), Allocatable: Complex array.
!   logic(:)    -- logical, allocatable: Logical array.
!   n           -- integer: Minimum size needed for 1-dimensional arrays.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than n. 
!                    Default is True.
!   init_val    -- optional: If present, init created array space to this value.
!                    If init_val is not present, no init will be done.
!
! Output:
!   str(:)      -- character(*), allocatable: Allocated array. 
!   var_str(:)  -- var_length_string_struct, allocatable: Allocated array.
!   inte(:)     -- integer, allocatable: Allocated array. 
!   re(:)       -- real(rp), Allocatable: Allocated array. 
!   a_ptr(:)    -- all_pointer_struct: Real pointer array.
!   cmpl(:)     -- complex(rp), Allocatable: Allocated Array.
!   logic(:)    -- logical, allocatable: Allocated array.
!-

interface re_allocate
  module procedure re_allocate_string
  module procedure re_allocate_var_string1
  module procedure re_allocate_var_string
  module procedure re_allocate_integer
  module procedure re_allocate_logical
  module procedure re_allocate_real
  module procedure re_allocate_all_pointer
  module procedure re_allocate_complex
end interface

!+
! Subroutine re_allocate2 (...)
!
! Subroutine to reallocate a 1-dim allocatable array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Note: For pointers to an array use the re_associate routine.
! Also see the re_allocate and re_allocate2d routines.
!
! Note: using exact = False can increase computation speed by
! preventing unneccessary deallocations/reallocations.
!
! Note: re_allocate2 is an overloaded name for: 
!   Subroutine re_allocate2_string (str, n_min, n_max, exact, init_val)
!   Subroutine re_allocate2_var_string (var_str, n_min, n_max, exact, init_val)
!   Subroutine re_allocate2_integer (inte, n_min, n_max, exact, init_val)
!   Subroutine re_allocate2_real (re, n_min, n_max, exact, init_val)
!   Subroutine re_allocate2_all_pointer (a_ptr, n_min, n_max, exact)
!   Subroutine re_allocate2_complex (cmpl, n_min, n_max, exact, init_val)
!   Subroutine re_allocate2_logical (logic, n_min, n_max, exact, init_val)
!
! Input:
!   str(:)        -- character(*), allocatable: String array.
!   var_str(:)    -- var_length_string_struct, allocatable: Variable length string structure.
!   inte(:)       -- integer, allocatable: Integer array.
!   re(:)         -- real(rp), Allocatable: Real array.
!   a_ptr(:)      -- all_pointer_struct: array of all_pointer_structs.
!   cmpl(:)       -- complex(rp), Allocatable: Complex array.
!   logic(:)      -- logical, allocatable: Logical array.
!   n_min         -- integer: Desired lower bound.
!   n_max         -- integer: Desired upper bound.
!   exact         -- logical, optional: If present and False then the size of 
!                      the output array is permitted to be larger than n. 
!                      Default is True.
!   init_val      -- optional: If present, init created array space to this value.
!                      If init_val is not present, no init will be done.
!
! Output:
!   str(:)        -- character(*), allocatable: Allocated array. 
!   var_str(:)    -- var_length_string_struct, allocatable: Allocated array.
!   inte(:)       -- integer, allocatable: Allocated array. 
!   re(:)         -- real(rp), Allocatable: Allocated array. 
!   a_ptr(:)      -- all_pointer_struct: Real pointer array.
!   cmpl(:)       -- complex(rp), Allocatable: Allocated Array.
!   logic(:)      -- logical, allocatable: Allocated array.
!-

interface re_allocate2
  module procedure re_allocate2_string
  module procedure re_allocate2_var_string
  module procedure re_allocate2_integer
  module procedure re_allocate2_logical
  module procedure re_allocate2_real
  module procedure re_allocate2_all_pointer
  module procedure re_allocate2_complex
end interface

!+
! Subroutine re_allocate2d (...)
!
! Subroutine to reallocate a 2-dim allocatable array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
! 
! Note: For pointers to an array use the re_associate routine.
! 
! Note: using exact = False can increase computation speed by
! preventing unneccessary deallocations/reallocations.
!
! Note: re_allocate2d is an overloaded name for:
!   Subroutine re_allocate_string2d (str2, n1, n2, exact, init_val, lb1, lb2)
!   Subroutine re_allocate_var_string2d (var_str2, n1, n2, exact, init_val, lb1, lb2)
!   Subroutine re_allocate_integer2d (inte2, n1, n2, exact, init_val, lb1, lb2)
!   Subroutine re_allocate_real2d (re2, n1, n2, exact, init_val, lb1, lb2)
!   Subroutine re_allocate_logical2d (logic2, n1, n2, exact, init_val, lb1, lb2)
!
! Input:
!   str2(:,:)       -- character(*), allocatable: String array.
!   var_str2(:,:)   -- var_length_string_struct, allocatable: Variable length string structure.
!   inte2(:,:)      -- integer, allocatable: Integer array.
!   re2(:,:)        -- real(rp), allocatable: Real array.
!   logic2(:,:)     -- logical, allocatable: Logical array.
!   n1, n2          -- integer: Minimum size needed for 2-dimensional arrays.
!   exact           -- logical, optional: If present and False then the size of 
!                        the output array is permitted to be larger than n. 
!                        Default is True.
!   init_val        -- optional: If present, init created array space to this value.
!                        If init_val is not present, no init will be done.
!   lb1, lb2        -- integer, optional: Array lower bounds. If not present, lower bounds are 1.
!
! Output:
!   str2(:,:)       -- character(*), allocatable: Allocated array.
!   var_str2(:,:)   -- var_length_string_struct, allocatable: Allocated array.
!   inte2(:,:)      -- integer, allocatable: Allocated array.
!   re2(:,:)        -- real(rp), allocatable: Allocated array.
!   logic2(:,:)     -- logical, allocatable: Allocated array. 
!-

interface re_allocate2d
  module procedure re_allocate_string2d
  module procedure re_allocate_var_string2d
  module procedure re_allocate_integer2d
  module procedure re_allocate_logical2d
  module procedure re_allocate_real2d
end interface

!-----------------------------------------
!+
! Subroutine re_associate (...)
!
! Subroutine to reassociate a 1-dim allocatable array.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Note: For allocatable arrays use the re_allocate routine.
! 
! Note: using exact = False can increase computation speed by
! preventing unneccessary deallocations/reallocations.
!
! Note: re_associate is an overloaded name for: 
!   Subroutine re_associate_string (str, n, exact, init_val)
!   Subroutine re_associate_var_string (var_str, n, exact, init_val)
!   Subroutine re_associate_integer (inte, n, exact, init_val)
!   Subroutine re_associate_real (re, n, exact, init_val)
!   Subroutine re_associate_logical (logic, n, exact, init_val)
!
! Input:
!   str(:)      -- Character(*), pointer: String array.
!   var_str(:)  -- var_length_string_struct, pointer: Variable length string structure.
!   inte(:)     -- integer, pointer: Integer array.
!   re(:)       -- real(rp), Pointer: Real array.
!   logic(:)    -- logical, pointer: Logical array.
!   n           -- integer: Minimum size needed.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than n. 
!                    Default is True.
!   init_val    -- optional: If present, init created array space to this value.
!                    If init_val is not present, no init will be done.
!
! Output:
!   str(:)      -- Character(*), pointer: Associated array with size(str) >= n.
!   var_str(:)  -- var_length_string_struct, pointer: Associated array with size(inte) >= n.
!   inte(:)     -- integer, pointer: Associated array with size(inte) >= n.
!   re(:)       -- real(rp), Pointer: Associated array with size(re) >= n.
!   logic(:)    -- logical, pointer: Associated array with size(logic) >= n.
!-

interface re_associate
  module procedure re_associate_string
  module procedure re_associate_var_string
  module procedure re_associate_integer
  module procedure re_associate_logical
  module procedure re_associate_real
end interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_string (str, n, exact, init_val)
!
! Routine to reallocate an array of strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   str(:) -- Character(*), allocatable: String array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   str(:) -- Character(*), allocatable: Allocated array with size(str) >= n.
!-

subroutine re_allocate_string (str, n, exact, init_val)

implicit none

integer n, n_old, n_save
character(*), allocatable :: str(:)
character(len(str)), allocatable :: temp_str(:)
character(*), optional :: init_val

logical, optional :: exact

!

if (allocated(str)) then
  n_old = size(str)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (str, temp_str)
  allocate (str(n))
  if (present(init_val)) str = init_val
  n_save = min(n, n_old)
  str(1:n_save) = temp_str(1:n_save)
  deallocate (temp_str)  
else
  allocate (str(n))
  if (present(init_val)) str = init_val
endif

end subroutine re_allocate_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_var_string (var_str, n, exact, init_val)
!
! Routine to reallocate an array of variable length strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   var_str(:)  -- var_length_string_struct, allocatable: String array.
!   n           -- integer: Size wanted.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than n. 
!                    Default is True.
!
! Output:
!   var_str(:) -- var_length_string_struct, allocatable: Allocated array with size(var_str) >= n.
!-

subroutine re_allocate_var_string (var_str, n, exact, init_val)

implicit none

integer n, n_old, n_save
type(var_length_string_struct), allocatable :: var_str(:)
type(var_length_string_struct), allocatable :: temp_var_str(:)
type(var_length_string_struct), optional :: init_val

logical, optional :: exact

!

if (allocated(var_str)) then
  n_old = size(var_str)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (var_str, temp_var_str)
  allocate (var_str(n))
  if (present(init_val)) var_str = init_val
  n_save = min(n, n_old)
  var_str(1:n_save) = temp_var_str(1:n_save)
  deallocate (temp_var_str)  
else
  allocate (var_str(n))
  if (present(init_val)) var_str = init_val
endif

end subroutine re_allocate_var_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_var_string1 (vstr, n, exact, init_val)
!
! Routine to reallocate a variable length string
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   vstr        -- character(:), allocatable: Variable length string.
!   n           -- integer: Size wanted.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than n. 
!                    Default is True.
!
! Output:
!   vstr(:)     -- character(:), allocatable: String with size(vstr) >= n.
!-

subroutine re_allocate_var_string1 (vstr, n, exact, init_val)

implicit none

integer n, n_old, n_save
character(:), allocatable :: vstr
character(:), allocatable :: temp_vstr
character, optional :: init_val

logical, optional :: exact

!

if (allocated(vstr)) then
  n_old = len(vstr)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (vstr, temp_vstr)
  allocate (character(n):: vstr)
  if (present(init_val)) vstr = init_val
  n_save = min(n, n_old)
  vstr(1:n_save) = temp_vstr(1:n_save)
  deallocate (temp_vstr)  
else
  allocate (character(n):: vstr)
  if (present(init_val)) vstr = init_val
endif

end subroutine re_allocate_var_string1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_integer (inte, n, exact, init_val)
!
! Routine to reallocate an array of integers.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   inte(:) -- integer, allocatable: Integer array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   inte(:) -- integer, allocatable: Allocated array with size(inte) >= n.
!-

subroutine re_allocate_integer (inte, n, exact, init_val)

implicit none

integer, allocatable :: inte(:), temp_inte(:)

integer, intent(in) :: n
integer n_save, n_old
integer, optional :: init_val

logical, optional :: exact

!

if (allocated(inte)) then
  n_old = size(inte)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (inte, temp_inte)
  allocate (inte(n))
  if (present(init_val)) inte = init_val
  n_save = min(n, n_old)
  inte(1:n_save) = temp_inte(1:n_save)
  deallocate (temp_inte)  
else
  allocate (inte(n))
  if (present(init_val)) inte = init_val
endif

end subroutine re_allocate_integer

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_complex (cmpl, n, exact, init_val)
!
! Routine to reallocate an array of complex numbers.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   cmpl(:)  -- Complex(rp), Allocatable: Complex array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   cmpl(:)  -- Complex(rp), Allocatable: Allocated array with size(cmpl) >= n.
!-

subroutine re_allocate_complex (cmpl, n, exact, init_val)

implicit none

complex(rp), allocatable :: cmpl(:), temp_cmpl(:)
complex(rp), optional :: init_val

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (allocated(cmpl)) then
  n_old = size(cmpl)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (cmpl, temp_cmpl)
  allocate (cmpl(n))
  if (present(init_val)) cmpl = init_val
  n_save = min(n, n_old)
  cmpl(1:n_save) = temp_cmpl(1:n_save)
  deallocate (temp_cmpl)  
else
  allocate (cmpl(n))
  if (present(init_val)) cmpl = init_val
endif

end subroutine re_allocate_complex

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_real (re, n, exact, init_val)
!
! Routine to reallocate an array of reals.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   re(:)  -- real(rp), Allocatable: Real array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   re(:)  -- real(rp), Allocatable: Allocated array with size(re) >= n.
!-

subroutine re_allocate_real (re, n, exact, init_val)

implicit none

real(rp), allocatable :: re(:), temp_re(:)
real(rp), optional :: init_val

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

end subroutine re_allocate_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_all_pointer (a_ptr, n, exact)
!
! Routine to reallocate an array of all_pointer_structs.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The pointers of the array are preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   a_ptr(:)  -- All_pointer_struct, Allocatable: All_pointer array.
!   n         -- integer: Size wanted.
!   exact     -- logical, optional: If present and False then the size of 
!                  the output array is permitted to be larger than n. 
!                  Default is True.
!
! Output:
!   a_ptr(:)  -- All_pointer_struct, Allocatable: Allocated array with size(re) >= n.
!-

subroutine re_allocate_all_pointer (a_ptr, n, exact)

implicit none

type(all_pointer_struct), allocatable :: a_ptr(:), temp_a(:)

integer, intent(in) :: n
integer n_save, n_old, i

logical, optional :: exact

!

if (allocated(a_ptr)) then
  n_old = size(a_ptr)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc(a_ptr, temp_a)
  allocate (a_ptr(n))
  n_save = min(n, n_old)
  a_ptr(1:n_save) = temp_a(1:n_save) 
  deallocate (temp_a)
  do i = n_save+1, n
    a_ptr(i)%r => null()
  enddo
else
  allocate (a_ptr(n))
endif

end subroutine re_allocate_all_pointer

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_logical (logic, n, exact, init_val)
!
! Routine to reallocate a string array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the size of the array
!
! Input:
!   logic(:) -- logical, allocatable: Logical array.
!   n      -- logical: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   logic(:) -- logical, allocatable: Allocated array with size(logic) >= n.
!-

subroutine re_allocate_logical (logic, n, exact, init_val)

implicit none

logical, allocatable :: logic(:), temp_logic(:)

integer, intent(in) :: n
integer i, n_save, n_old

logical, optional :: exact
logical, optional :: init_val

!

if (allocated(logic)) then
  n_old = size(logic)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc(logic, temp_logic)
  allocate (logic(n))
  if (present(init_val)) logic = init_val
  n_save = min(n, n_old)
  do i = 1, n_save
    call transfer_logical (temp_logic(i), logic(i)) 
  enddo
  deallocate (temp_logic)  
else
  allocate (logic(n))
  if (present(init_val)) logic = init_val
endif

end subroutine re_allocate_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_string (str, n1, n2, exact, init_val)
!
! Routine to reallocate an array of strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   str(:) -- Character(*), allocatable: String array.
!   n1     -- integer: Desired lower bound.
!   n2     -- integer: Desired upper bound.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   str(:) -- Character(*), allocatable: Allocated array with 
!               bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_string (str, n1, n2, exact, init_val)

implicit none

integer n1, n2, n1_old, n2_old, n1_save, n2_save
character(*), allocatable :: str(:)
character(len(str)), allocatable :: temp_str(:)
character(*), optional :: init_val

logical, optional :: exact

!

if (allocated(str)) then
  n1_old = lbound(str, 1); n2_old = ubound(str, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(str, temp_str)
  allocate (str(n1:n2))
  if (present(init_val)) str = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  str(n1_save:n2_save) = temp_str(n1_save:n2_save)
  deallocate (temp_str)  
else
  allocate (str(n1:n2))
  if (present(init_val)) str = init_val
endif

end subroutine re_allocate2_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_var_string (var_str, n1, n2, exact, init_val)
!
! Routine to reallocate an array of var_strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   var_str(:) -- var_length_string_struct, allocatable: Var_string array.
!   n1         -- integer: Desired lower bound.
!   n2         -- integer: Desired upper bound.
!   exact      -- logical, optional: If present and False then the size of 
!                   the output array is permitted to be larger than [n1, n2]. 
!                   Default is True.
!
! Output:
!   var_str(:) -- var_length_string_struct, allocatable: Allocated array with 
!                   bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_var_string (var_str, n1, n2, exact, init_val)

implicit none

integer n1, n2, n1_old, n2_old, n1_save, n2_save
type(var_length_string_struct), allocatable :: var_str(:)
type(var_length_string_struct), allocatable :: temp_var_str(:)
type(var_length_string_struct), optional :: init_val

logical, optional :: exact

!

if (allocated(var_str)) then
  n1_old = lbound(var_str, 1); n2_old = ubound(var_str, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(var_str, temp_var_str)
  allocate (var_str(n1:n2))
  if (present(init_val)) var_str = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  var_str(n1_save:n2_save) = temp_var_str(n1_save:n2_save)
  deallocate (temp_var_str)  
else
  allocate (var_str(n1:n2))
  if (present(init_val)) var_str = init_val
endif

end subroutine re_allocate2_var_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_integer (inte, n1, n2, exact, init_val)
!
! Routine to reallocate an array of integers.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   inte(:) -- integer, allocatable: Integer array.
!   n1     -- integer: Desired lower bound.
!   n2     -- integer: Desired upper bound.
!   exact   -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   inte(:) -- integer, allocatable: Allocated array with 
!                bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_integer (inte, n1, n2, exact, init_val)

implicit none

integer, allocatable :: inte(:), temp_inte(:)

integer, intent(in) :: n1, n2
integer n1_save, n2_save, n1_old, n2_old
integer, optional :: init_val

logical, optional :: exact

!

if (allocated(inte)) then
  n1_old = lbound(inte, 1); n2_old = ubound(inte, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(inte, temp_inte)
  allocate (inte(n1:n2))
  if (present(init_val)) inte = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  inte(n1_save:n2_save) = temp_inte(n1_save:n2_save)
  deallocate (temp_inte)  
else
  allocate (inte(n1:n2))
  if (present(init_val)) inte = init_val
endif

end subroutine re_allocate2_integer

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_complex (cmpl, n1, n2, exact, init_val)
!
! Routine to reallocate an array of complex numbers.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   cmpl(:)  -- Complex(rp), Allocatable: Complex array.
!   n1     -- integer: Desired lower bound.
!   n2     -- integer: Desired upper bound.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   cmpl(:)  -- Complex(rp), Allocatable: Allocated array with 
!                bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_complex (cmpl, n1, n2, exact, init_val)

implicit none

complex(rp), allocatable :: cmpl(:), temp_cmpl(:)
complex(rp), optional :: init_val

integer, intent(in) :: n1, n2
integer n1_save, n2_save, n1_old, n2_old

logical, optional :: exact

!

if (allocated(cmpl)) then
  n1_old = lbound(cmpl, 1); n2_old = ubound(cmpl, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(cmpl, temp_cmpl)
  allocate (cmpl(n1:n2))
  if (present(init_val)) cmpl = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  cmpl(n1_save:n2_save) = temp_cmpl(n1_save:n2_save)
  deallocate (temp_cmpl)  
else
  allocate (cmpl(n1:n2))
  if (present(init_val)) cmpl = init_val
endif

end subroutine re_allocate2_complex

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_real (re, n1, n2, exact, init_val)
!
! Routine to reallocate an array of reals.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   re(:)  -- real(rp), Allocatable: Real array.
!   n1     -- integer: Desired lower bound.
!   n2     -- integer: Desired upper bound.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   re(:)  -- real(rp), Allocatable: Allocated array with 
!                bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_real (re, n1, n2, exact, init_val)

implicit none

real(rp), allocatable :: re(:), temp_re(:)
real(rp), optional :: init_val

integer, intent(in) :: n1, n2
integer n1_save, n2_save, n1_old, n2_old

logical, optional :: exact

!

if (allocated(re)) then
  n1_old = lbound(re, 1); n2_old = ubound(re, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(re, temp_re)
  allocate (re(n1:n2))
  if (present(init_val)) re = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  re(n1_save:n2_save) = temp_re(n1_save:n2_save)
  deallocate (temp_re)  
else
  allocate (re(n1:n2))
  if (present(init_val)) re = init_val
endif

end subroutine re_allocate2_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_all_pointer (a_ptr, n1, n2, exact)
!
! Routine to reallocate an array of all_pointer_structs
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   a_ptr(:) -- All_pointer_struct, allocatable: Real array.
!   n1       -- integer: Desired lower bound.
!   n2       -- integer: Desired upper bound.
!   exact    -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   a_ptr(:)  -- All_pointer_struct, allocatable: Allocated array with 
!                  bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_all_pointer (a_ptr, n1, n2, exact)

implicit none

type(all_pointer_struct), allocatable :: a_ptr(:), temp_a_ptr(:)

integer, intent(in) :: n1, n2
integer n1_save, n2_save, n1_old, n2_old, i

logical, optional :: exact

!

if (allocated(a_ptr)) then
  n1_old = lbound(a_ptr, 1); n2_old = ubound(a_ptr, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(a_ptr, temp_a_ptr)
  allocate (a_ptr(n1:n2))
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  a_ptr(n1_save:n2_save) = temp_a_ptr(n1_save:n2_save)
  deallocate (temp_a_ptr)  

else
  allocate (a_ptr(n1:n2))
endif

end subroutine re_allocate2_all_pointer

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate2_logical (logic, n1, n2, exact, init_val)
!
! Routine to reallocate a logical array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the size of the array
!
! Input:
!   logic(:) -- logical, allocatable: Logical array.
!   n1     -- integer: Desired lower bound.
!   n2     -- integer: Desired upper bound.
!   exact    -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   logic(:) -- logical, allocatable: Allocated array with 
!                 bounds spanning at least [n1, n2]
!-

subroutine re_allocate2_logical (logic, n1, n2, exact, init_val)

implicit none

logical, allocatable :: logic(:), temp_logic(:)
logical, optional :: init_val

integer, intent(in) :: n1, n2
integer n1_save, n2_save, n1_old, n2_old

logical, optional :: exact

!

if (allocated(logic)) then
  n1_old = lbound(logic, 1); n2_old = ubound(logic, 1)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old) return
  call move_alloc(logic, temp_logic)
  allocate (logic(n1:n2))
  if (present(init_val)) logic = init_val
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  logic(n1_save:n2_save) = temp_logic(n1_save:n2_save)
  deallocate (temp_logic)  
else
  allocate (logic(n1:n2))
  if (present(init_val)) logic = init_val
endif

end subroutine  re_allocate2_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_string2d (str2, n1, n2, exact, init_val, lb1, lb2)
!
! Routine to reallocate an array of strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   str2(:,:) -- Character(*), allocatable: String array.
!   n1, n2    -- integer: Size wanted.
!   exact     -- logical, optional: If present and False then the size of 
!                  the output array is permitted to be larger than [n1, n2]. 
!                  Default is True.
!
! Output:
!   str2(:,:) -- Character(*), allocatable: Allocated array with 
!                  bounds spanning at least [n1, n2]
!-

subroutine re_allocate_string2d (str2, n1, n2, exact, init_val, lb1, lb2)

implicit none

integer n1, n2, n1_old, n2_old, n1_save, n2_save, l1, l2
integer, optional :: lb1, lb2
character(*), allocatable :: str2(:,:)
character(len(str2)), allocatable :: temp_str2(:,:)
character(*), optional :: init_val

logical, optional :: exact

!

l1 = integer_option(1, lb1)
l2 = integer_option(1, lb2)

if (allocated(str2)) then
  n1_old = ubound(str2, 1); n2_old = ubound(str2, 2)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2_old <= n2) return
  call move_alloc(str2, temp_str2)
  allocate (str2(l1:n1, l2:n2))
  if (present(init_val)) str2 = init_val
  n1_save = min(n1, n1_old); n2_save = min(n2, n2_old)
  str2(l1:n1_save,l2:n2_save) = temp_str2(l1:n1_save,l2:n2_save)
  deallocate (temp_str2)  
else
  allocate (str2(l1:n1, l2:n2))
  if (present(init_val)) str2 = init_val
endif

end subroutine re_allocate_string2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_var_string2d (var_str2, n1, n2, exact, init_val, lb1, lb2)
!
! Routine to reallocate an array of var_strings.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   var_str2(:,:) -- var_length_string_struct, allocatable: Var_string array.
!   n1, n2        -- integer: Size wanted.
!   exact         -- logical, optional: If present and False then the size of 
!                      the output array is permitted to be larger than [n1, n2]. 
!                      Default is True.
!
! Output:
!   var_str2(:,:) -- var_length_string_struct, allocatable: Allocated array with 
!                      bounds spanning at least [n1, n2]
!-

subroutine re_allocate_var_string2d (var_str2, n1, n2, exact, init_val, lb1, lb2)

implicit none

integer n1, n2, n1_old, n2_old, n1_save, n2_save, l1, l2
integer, optional :: lb1, lb2
type(var_length_string_struct), allocatable :: var_str2(:,:)
type(var_length_string_struct), allocatable :: temp_var_str2(:,:)
type(var_length_string_struct), optional :: init_val

logical, optional :: exact

!

l1 = integer_option(1, lb1)
l2 = integer_option(1, lb2)

if (allocated(var_str2)) then
  n1_old = ubound(var_str2, 1); n2_old = ubound(var_str2, 2)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2_old <= n2) return
  call move_alloc(var_str2, temp_var_str2)
  allocate (var_str2(l1:n1, l2:n2))
  if (present(init_val)) var_str2 = init_val
  n1_save = min(n1, n1_old); n2_save = min(n2, n2_old)
  var_str2(l1:n1_save,l2:n2_save) = temp_var_str2(l1:n1_save,l2:n2_save)
  deallocate (temp_var_str2)  
else
  allocate (var_str2(l1:n1, l2:n2))
  if (present(init_val)) var_str2 = init_val
endif

end subroutine re_allocate_var_string2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_integer2d (inte2, n1, n2, exact, init_val, lb1, lb2)
!
! Routine to reallocate an array of integers.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   inte2(:,:) -- integer, allocatable: Integer array.
!   n1, n2     -- integer: Size wanted.
!   exact      -- logical, optional: If present and False then the size of 
!                   the output array is permitted to be larger than [n1, n2]. 
!                   Default is True.
!
! Output:
!   inte2(:,:) -- integer, allocatable: Allocated array with 
!                   bounds spanning at least [n1, n2]
!-

subroutine re_allocate_integer2d (inte2, n1, n2, exact, init_val, lb1, lb2)

implicit none

integer, allocatable :: inte2(:,:), temp_inte2(:,:)
integer, optional :: init_val

integer, intent(in) :: n1, n2
integer, optional :: lb1, lb2
integer n1_save, n2_save, l1, l2, n1_old, n2_old

logical, optional :: exact

!

l1 = integer_option(1, lb1)
l2 = integer_option(1, lb2)

if (allocated(inte2)) then
  n1_old = ubound(inte2, 1); n2_old = ubound(inte2, 2)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2_old <= n2) return
  call move_alloc(inte2, temp_inte2)
  allocate (inte2(l1:n1, l2:n2))
  if (present(init_val)) inte2 = init_val
  n1_save = min(n1, n1_old); n2_save = min(n2, n2_old)
  inte2(l1:n1_save,l2:n2_save) = temp_inte2(l1:n1_save,l2:n2_save)
  deallocate (temp_inte2)  
else
  allocate (inte2(l1:n1, l2:n2))
  if (present(init_val)) inte2 = init_val
endif

end subroutine re_allocate_integer2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_real2d (re2, n1, n2, exact, init_val, lb1, lb2)
!
! Routine to reallocate an array of reals.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the original size of the array
!
! Input:
!   re2(:,:) -- real(rp), Allocatable: Real array.
!   n1, n2   -- integer: Size wanted.
!   exact    -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than [n1, n2]. 
!                 Default is True.
!
! Output:
!   re2(:,:) -- real(rp), Allocatable: Allocated array with 
!                 bounds spanning at least [n1, n2]
!-

subroutine re_allocate_real2d (re2, n1, n2, exact, init_val, lb1, lb2)

implicit none

real(rp), allocatable :: re2(:,:), temp_re2(:,:)
real(rp), optional :: init_val

integer, intent(in) :: n1, n2
integer, optional :: lb1, lb2
integer n1_save, n2_save, l1, l2, n1_old, n2_old

logical, optional :: exact

!

l1 = integer_option(1, lb1)
l2 = integer_option(1, lb2)

if (allocated(re2)) then
  n1_old = ubound(re2, 1); n2_old = ubound(re2, 2)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2_old <= n2) return
  n1_save = min(n1, n1_old); n2_save = min(n2, n2_old)
  allocate (temp_re2(n1_save,n2_save))
  temp_re2 = re2(l1:n1_save,l2:n2_save)
  deallocate (re2)
  allocate (re2(l1:n1, l2:n2))
  if (present(init_val)) re2 = init_val
  re2(l1:n1_save,l2:n2_save) = temp_re2
  deallocate (temp_re2)  
else
  allocate (re2(l1:n1, l2:n2))
  if (present(init_val)) re2 = init_val
endif

end subroutine re_allocate_real2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_allocate_logical2d (logic2, n1, n2, exact, init_val, lb1, lb2)
!
! Routine to reallocate a logical array.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if [n1, n2] is less than the size of the array
!
! Input:
!   logic2(:,:) -- logical, allocatable: Logical array.
!   n1, n2      -- integer: Size wanted.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than [n1, n2]. 
!                    Default is True.
!
! Output:
!   logic2(:,:) -- logical, allocatable: Allocated array with 
!                    bounds spanning at least [n1, n2]
!-

subroutine re_allocate_logical2d (logic2, n1, n2, exact, init_val, lb1, lb2)

implicit none

logical, allocatable :: logic2(:,:), temp_logic2(:,:)
logical, optional :: init_val

integer, intent(in) :: n1, n2
integer, optional :: lb1, lb2
integer n1_save, n2_save, l1, l2, n1_old, n2_old

logical, optional :: exact

!

l1 = integer_option(1, lb1)
l2 = integer_option(1, lb2)

if (allocated(logic2)) then
  n1_old = ubound(logic2, 1); n2_old = ubound(logic2, 2)
  if (n1 == n1_old .and. n2 == n2_old) return
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2_old <= n2) return
  call move_alloc(logic2, temp_logic2)
  allocate (logic2(l1:n1, l2:n2))
  if (present(init_val)) logic2 = init_val
  n1_save = min(n1, n1_old); n2_save = min(n2, n2_old)
  logic2(l1:n1_save,l2:n2_save) = temp_logic2(l1:n1_save,l2:n2_save)
  deallocate (temp_logic2)  
else
  allocate (logic2(l1:n1, l2:n2))
  if (present(init_val)) logic2 = init_val
endif

end subroutine re_allocate_logical2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_string (str, n, exact, init_val)
!
! Routine to reassociate an array of strings.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   str(:) -- Character(*), pointer: String array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   str(:) -- Character(*), pointer: Allocated array with size(str) >= n.
!-

subroutine re_associate_string (str, n, exact, init_val)

implicit none

character(*), pointer :: str(:)
character(len(str)), pointer :: temp_str(:)
character(*), optional :: init_val

integer n, n_old, n_save

logical, optional :: exact

!

if (associated(str)) then
  n_old = size(str)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_str => str
  allocate (str(n))
  if (present(init_val)) str = init_val
  str(1:n_save) = temp_str
  deallocate (temp_str)  
else
  allocate (str(n))
  if (present(init_val)) str = init_val
endif

end subroutine re_associate_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_var_string (var_str, n, exact, init_val)
!
! Routine to reassociate an array of var_strings.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   var_str(:) -- var_length_string_struct, pointer: Var_string array.
!   n          -- integer: Size wanted.
!   exact      -- logical, optional: If present and False then the size of 
!                   the output array is permitted to be larger than n. 
!                   Default is True.
!
! Output:
!   var_str(:) -- var_length_string_struct, pointer: Allocated array with size(var_str) >= n.
!-

subroutine re_associate_var_string (var_str, n, exact, init_val)

implicit none

type(var_length_string_struct), pointer :: var_str(:)
type(var_length_string_struct), pointer :: temp_var_str(:)
type(var_length_string_struct), optional :: init_val

integer n, n_old, n_save

logical, optional :: exact

!

if (associated(var_str)) then
  n_old = size(var_str)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_var_str => var_str
  allocate (var_str(n))
  if (present(init_val)) var_str = init_val
  var_str(1:n_save) = temp_var_str
  deallocate (temp_var_str)  
else
  allocate (var_str(n))
  if (present(init_val)) var_str = init_val
endif

end subroutine re_associate_var_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_integer (inte, n, exact, init_val)
!
! Routine to reassociate an array of integers.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   inte(:) -- integer, pointer: Integer array.
!   n       -- integer: Size wanted.
!   exact   -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   inte(:) -- integer, pointer: Allocated array with size(inte) >= n.
!-

subroutine re_associate_integer (inte, n, exact, init_val)

implicit none

integer, pointer :: inte(:), temp_inte(:)

integer, intent(in) :: n
integer n_save, n_old
integer, optional :: init_val

logical, optional :: exact

!

if (associated(inte)) then
  n_old = size(inte)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_inte => inte
  allocate (inte(n))
  if (present(init_val)) inte = init_val
  inte(1:n_save) = temp_inte
  deallocate (temp_inte)  
else
  allocate (inte(n))
  if (present(init_val)) inte = init_val
endif

end subroutine re_associate_integer

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_real (re, n, exact, init_val)
!
! Routine to reassociate an array of reals.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   re(:)  -- real(rp), Pointer: Real array.
!   n      -- integer: Size wanted.
!   exact  -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   re(:)  -- real(rp), Pointer: Allocated array with size(re) >= n.
!-

subroutine re_associate_real (re, n, exact, init_val)

implicit none

real(rp), pointer :: re(:), temp_re(:)
real(rp), optional :: init_val

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (associated(re)) then
  n_old = size(re)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_re => re
  allocate (re(n))
  if (present(init_val)) re = init_val
  re(1:n_save) = temp_re
  deallocate (temp_re)  
else
  allocate (re(n))
  if (present(init_val)) re = init_val
endif

end subroutine re_associate_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_logical (logic, n, exact, init_val)
!
! Routine to reassociate a string array.
! This is modeled after the reassociate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the size of the array
!
! Input:
!   logic(:) -- logical, pointer: Logical array.
!   n        -- logical: Size wanted.
!   exact    -- logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   logic(:) -- logical, pointer: Allocated array with size(logic) >= n.
!-

subroutine re_associate_logical (logic, n, exact, init_val)

implicit none

logical, pointer :: logic(:), temp_logic(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact
logical, optional :: init_val

!

if (associated(logic)) then
  n_old = size(logic)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_logic => logic
  allocate (logic(n))
  if (present(init_val)) logic = init_val
  logic(1:n_save) = temp_logic
  deallocate (temp_logic)  
else
  allocate (logic(n))
  if (present(init_val)) logic = init_val
endif

end subroutine re_associate_logical

end module
