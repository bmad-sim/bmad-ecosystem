!+
! Module complex_taylor_mod
!
! Note: This module is an exact copy of taylor_mod, with:
!   taylor -> complex_taylor
!   real -> complex
!   format modifications to: type_complex_taylors 
!-

module complex_taylor_mod

use sim_utils

! Note: the complex_taylor_struct uses the Bmad standard (x, p_x, y, p_y, z, p_z) 
! the universal_complex_taylor in Etienne's PTC uses (x, p_x, y, p_y, p_z, -c*t)
! %ref is the reference point about which the complex_taylor expansion was made

type complex_taylor_term_struct
  complex(rp) :: coef
  integer :: expn(6)  
end type

type complex_taylor_struct
  complex (rp) :: ref = 0
  type (complex_taylor_term_struct), pointer :: term(:) => null()
end type

!

interface assignment (=)
  module procedure complex_taylor_equal_complex_taylor
  module procedure complex_taylors_equal_complex_taylors
end interface

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function complex_taylor_coef (bmad_complex_taylor, exp)
! Function complex_taylor_coef (bmad_complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series.
!
! Note: complex_taylor_coef is overloaded by:
!   complex_taylor_coef1 (bmad_complex_taylor, exp)
!   complex_taylor_coef2 (bmad_complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
! Using the complex_taylor_coef2 form limits obtaining coefficients to 9th order
! or less. Also: complex_taylor_coef2 does not check that all i1, ..., i9 are between
! 1 and 6.
!
! For example: To get the 2nd order term corresponding to 
!   y(out) = Coef * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (complex_taylor_struct) bmad_complex_taylor(6)      ! complex_taylor Map
!   ...
!   coef = complex_taylor_coef (bmad_complex_taylor(3), 6, 6)  ! 1st possibility or ...
!   coef = complex_taylor_coef (bmad_complex_taylor(3), [0, 0, 0, 0, 0, 2 ])  
!
! Modules needed:
!   use bmad
!
! Input (complex_taylor_coef1):
!   bmad_complex_taylor -- complex_taylor_struct: complex_taylor series.
!   exp(6)      -- Integer: Array of exponent indices.
!
! Input (complex_taylor_coef2):
!   bmad_complex_taylor -- complex_taylor_struct: complex_taylor series.
!   i1, ..., i9 -- Integer, optional: indexes (each between 1 and 6).
!
! Output:
!   complex_taylor_coef -- complex(rp): Coefficient.
!-

interface complex_taylor_coef
  module procedure complex_taylor_coef1
  module procedure complex_taylor_coef2
end interface

private complex_taylor_coef1, complex_taylor_coef2

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine add_complex_taylor_term (bmad_complex_taylor, coef, exp)
! Subroutine add_complex_taylor_term (bmad_complex_taylor, coef, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! If bmad_complex_taylor does not have a term with the same exponents as expn then a new
! term is added to bmad_complex_taylor so the total number of terms is increased by one.
!
! If bmad_complex_taylor already has a term with the same exponents then
! the "replace" argument determines what happens:
!   If replace = False (default) then
!      coef is added to the coefficient of the old term.
!   If replace = True then
!      coef replaces the coefficient of the old term.
! In both these cases, the number of terms in bmad_complex_taylor remains the same.
!
! Note: add_complex_taylor_term is overloaded by:
!   add_complex_taylor_term1 (bmad_complex_taylor, exp, replace)
!   add_complex_taylor_term2 (bmad_complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)
! Using the add_complex_taylor_term2 form limits obtaining coefficients to 9th order
! or less. Also: add_complex_taylor_term2 does not check that all i1, ..., i9 are between
! 1 and 6.
!
! For example: To add the 2nd order term corresponding to:
!   y(out) = 1.34 * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (complex_taylor_struct) bmad_complex_taylor(6)      ! complex_taylor Map
!   ...
!   coef = add_complex_taylor_term (bmad_complex_taylor(3), 1.34_rp, 6, 6)  ! 1st possibility or ...
!   coef = add_complex_taylor_term (bmad_complex_taylor(3), 1.34_rp, [0, 0, 0, 0, 0, 2 ])  
!
! Modules needed:
!   use bmad
!
! Input (add_complex_taylor_term1):
!   bmad_complex_taylor -- complex_taylor_struct: complex_taylor series.
!   coef        -- complex(rp): Coefficient.
!   exp(6)      -- Integer: Array of exponent indices.
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Input (add_complex_taylor_term2):
!   bmad_complex_taylor -- complex_taylor_struct: complex_taylor series.
!   coef        -- complex(rp): Coefficient.
!   i1, ..., i9 -- Integer, optional: Exponent indexes (each between 1 and 6).
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Output:
!   bmad_complex_taylor -- complex_taylor_struct: New series with term added
!-

interface add_complex_taylor_term
  module procedure add_complex_taylor_term1
  module procedure add_complex_taylor_term2
end interface

private add_complex_taylor_term1, add_complex_taylor_term2

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)
!
! Subroutine that is used to set one complex_taylor equal to another. 
! This routine takes care of the pointers in complex_taylor1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		complex_taylor1 = complex_taylor2
!
! Input:
!   complex_taylor2 -- complex_taylor_struct: Input complex_taylor.
!
! Output:
!   complex_taylor1 -- complex_taylor_struct: Output complex_taylor.
!-

subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)

implicit none
	
type (complex_taylor_struct), intent(inout) :: complex_taylor1
type (complex_taylor_struct), intent(in) :: complex_taylor2

!

complex_taylor1%ref = complex_taylor2%ref

if (associated(complex_taylor2%term)) then
  call init_complex_taylor_series (complex_taylor1, size(complex_taylor2%term))
  complex_taylor1%term = complex_taylor2%term
else
  if (associated (complex_taylor1%term)) deallocate (complex_taylor1%term)
endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)
!
! Subroutine to transfer the values from one complex_taylor map to another:
!     complex_taylor1 <= complex_taylor2
!
! Modules needed:
!   use bmad
!
! Input:
!   complex_taylor2(:) -- complex_taylor_struct: complex_taylor map.
!
! Output:
!   complex_taylor1(:) -- complex_taylor_struct: complex_taylor map. 
!-

subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)

implicit none

type (complex_taylor_struct), intent(inout) :: complex_taylor1(:)
type (complex_taylor_struct), intent(in)    :: complex_taylor2(:)

integer i

!

do i = 1, size(complex_taylor1)
  complex_taylor1(i) = complex_taylor2(i)
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function complex_taylor_coef1 (bmad_complex_taylor, exp)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series. This routine is used by the overloaded function
! complex_taylor_coef. See complex_taylor_coef for more details.
!-

function complex_taylor_coef1 (bmad_complex_taylor, exp) result (coef)

implicit none

type (complex_taylor_struct), intent(in) :: bmad_complex_taylor

complex(rp) coef

integer, intent(in) :: exp(:)
integer i

!

coef = 0

do i = 1, size(bmad_complex_taylor%term)
  if (all(bmad_complex_taylor%term(i)%expn == exp)) then
    coef = bmad_complex_taylor%term(i)%coef
    return
  endif
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function complex_taylor_coef2 (bmad_complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series. This routine is used by the overloaded function
! complex_taylor_coef. See complex_taylor_coef for more details.
!-

function complex_taylor_coef2 (bmad_complex_taylor, i1, i2, i3, &
                          i4, i5, i6, i7, i8, i9) result (coef)

implicit none

type (complex_taylor_struct), intent(in) :: bmad_complex_taylor

complex(rp) coef

integer, intent(in), optional :: i1, i2, i3, i4, i5, i6, i7, i8, i9
integer i, exp(6)

!

exp = 0
if (present (i1)) exp(i1) = exp(i1) + 1
if (present (i2)) exp(i2) = exp(i2) + 1
if (present (i3)) exp(i3) = exp(i3) + 1
if (present (i4)) exp(i4) = exp(i4) + 1
if (present (i5)) exp(i5) = exp(i5) + 1
if (present (i6)) exp(i6) = exp(i6) + 1
if (present (i7)) exp(i7) = exp(i7) + 1
if (present (i8)) exp(i8) = exp(i8) + 1
if (present (i9)) exp(i9) = exp(i9) + 1

coef = 0

do i = 1, size(bmad_complex_taylor%term)
  if (all(bmad_complex_taylor%term(i)%expn == exp)) then
    coef = bmad_complex_taylor%term(i)%coef
    return
  endif
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine type_complex_taylors (bmad_complex_taylor, max_order, lines, n_lines)
!
! Subroutine to print or put in a string array a Bmad complex taylor map.
! If the lines(:) argument is not present, the element information is printed to the terminal.
!
! Moudles needed:
!   use bmad
!
! Input:
!   bmad_complex_taylor(6) -- complex_taylor_struct: Array of complex_taylors.
!   max_order      -- Integer, optional: Maximum order to print.
!
! Output:
!   lines(:)     -- Character(100), allocatable, optional :: Character array to hold the output. 
!                     If not present, the information is printed to the terminal.
!   n_lines      -- Integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_complex_taylors (bmad_complex_taylor, max_order, lines, n_lines)

use re_allocate_mod

implicit none

type (complex_taylor_struct), intent(in), target :: bmad_complex_taylor(6)
type (complex_taylor_term_struct), pointer :: tt
type (complex_taylor_struct) tlr

integer, optional, intent(out) :: n_lines
integer, optional :: max_order
integer i, j, k, nl, ix

character(*), optional, allocatable :: lines(:)
character(100), allocatable :: li(:)
character(100) fmt1, fmt2, fmt

! If not allocated then not much to do

if (.not. associated(bmad_complex_taylor(1)%term)) then
  nl = 2
  allocate (li(nl))
  li(1) = '---------------------------------------------------'
  li(2) = 'A complex_taylor Map Does Not Exist.' 

! Normal case

else
  nl = 14 + sum( [(size(bmad_complex_taylor(i)%term), i = 1, 6) ])
  allocate(li(nl))

  write (li(1), *) 'complex_taylor Terms:'
  write (li(2), *) &
         'Out       Re(coef)            Im(coef)        Exponents           Order        Reference'
  nl = 2


  fmt1 = '(i4, a, 2f20.12, 6i3, i9)'
  fmt2 = '(i4, a, 1p, 2d20.11, 0p, 6i3, i9)'

  do i = 1, 6
    nl=nl+1; li(nl) = ' -----------------------------------------------------------------------'

    nullify (tlr%term)
    
    if (.not. associated(bmad_complex_taylor(i)%term) ) cycle
    
    call sort_complex_taylor_terms (bmad_complex_taylor(i), tlr)

    do j = 1, size(bmad_complex_taylor(i)%term)

      tt => tlr%term(j)

      if (present(max_order)) then
        if (sum(tt%expn) > max_order) cycle
      endif

      if (abs(tt%coef) < 1d5) then
        fmt = fmt1
      else
        fmt = fmt2
      endif

      !if (j == 1) then
      !  nl=nl+1; write (li(nl), fmt) i, ':', tt%coef, &
      !              (tt%expn(k), k = 1, 6), sum(tt%expn), bmad_complex_taylor(i)%ref
      !else
      if (j==1) then
        nl=nl+1
        write(li(nl), '(a, 2f18.9)') '   Reference: ', bmad_complex_taylor(i)%ref
      endif
        nl=nl+1; write (li(nl), fmt) i, ':', tt%coef, &
                    (tt%expn(k), k = 1, 6), sum(tt%expn)
      !endif
    enddo

    deallocate (tlr%term)

  enddo
endif

! Finish

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print '(1x, a)', trim(li(i))
  enddo
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine complex_taylor_make_unit (bmad_complex_taylor)
!
! Subroutine to make the unit complex_taylor map:
!       r(out) = Map * r(in) = r(in)
!
! Modules needed:
!   use bmad
!
! Output:
!   bmad_complex_taylor(6) -- complex_taylor_struct: Unit complex_taylor map .
!-

subroutine complex_taylor_make_unit (bmad_complex_taylor)

implicit none

type (complex_taylor_struct) bmad_complex_taylor(:)
integer i

do i = 1, size(bmad_complex_taylor)
  call init_complex_taylor_series (bmad_complex_taylor(i), 1)
  bmad_complex_taylor(i)%term(1)%coef = 1.0
  bmad_complex_taylor(i)%term(1)%expn = 0
  bmad_complex_taylor(i)%term(1)%expn(i) = 1
  bmad_complex_taylor(i)%ref = 0
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_complex_taylor_term1 (bmad_complex_taylor, coef, expn, replace)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! This routine is used by the overloaded function add_complex_taylor_term. 
! See add_complex_taylor_term for more details.
!-

subroutine add_complex_taylor_term1 (bmad_complex_taylor, coef, expn, replace)

implicit none

type (complex_taylor_struct) bmad_complex_taylor

complex(rp) coef
integer expn(:), i, n

logical, optional :: replace

! Search for an existing term of the same type

n = size(bmad_complex_taylor%term)

do i = 1, n
  if (all(bmad_complex_taylor%term(i)%expn == expn)) then
    if (logic_option(.false., replace)) then
      bmad_complex_taylor%term(i)%coef = coef
    else
      bmad_complex_taylor%term(i)%coef = coef + bmad_complex_taylor%term(i)%coef
    endif
    if (bmad_complex_taylor%term(i)%coef == 0) then  ! Kill this term
      bmad_complex_taylor%term(i:n-1) = bmad_complex_taylor%term(i+1:n)
      call init_complex_taylor_series (bmad_complex_taylor, n-1, .true.)
    endif
    return
  endif
enddo

! new term

call init_complex_taylor_series (bmad_complex_taylor, n+1, .true.)
bmad_complex_taylor%term(n+1)%coef = coef
bmad_complex_taylor%term(n+1)%expn = expn

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_complex_taylor_term2 (bmad_complex_taylor, coef, 
!                        i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! This routine is used by the overloaded function add_complex_taylor_term. 
! See add_complex_taylor_term for more details.
!-

subroutine add_complex_taylor_term2 (bmad_complex_taylor, coef, &
                        i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)

implicit none

type (complex_taylor_struct) bmad_complex_taylor

complex(rp) coef

integer, intent(in), optional :: i1, i2, i3, i4, i5, i6, i7, i8, i9
integer i, n, expn(6)

logical, optional :: replace

! 

expn = 0
if (present (i1)) expn(i1) = expn(i1) + 1
if (present (i2)) expn(i2) = expn(i2) + 1
if (present (i3)) expn(i3) = expn(i3) + 1
if (present (i4)) expn(i4) = expn(i4) + 1
if (present (i5)) expn(i5) = expn(i5) + 1
if (present (i6)) expn(i6) = expn(i6) + 1
if (present (i7)) expn(i7) = expn(i7) + 1
if (present (i8)) expn(i8) = expn(i8) + 1
if (present (i9)) expn(i9) = expn(i9) + 1

call add_complex_taylor_term1 (bmad_complex_taylor, coef, expn, replace)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_complex_taylor_series (bmad_complex_taylor, n_term, save)
!
! Subroutine to initialize a Bmad complex_taylor series (6 of these series make
! a complex_taylor map). Note: This routine does not zero the structure. The calling
! routine is responsible for setting all values.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_complex_taylor -- complex_taylor_struct: Old structure.
!   n_term      -- Integer: Number of terms to allocate. 
!                   n_term < 1 => bmad_complex_taylor%term pointer will be disassociated.
!   save        -- Logical, optional: If True then save any old terms when
!                   bmad_complex_taylor is resized. Default is False.
!
! Output:
!   bmad_complex_taylor -- complex_taylor_struct: Initalized structure.
!-

subroutine init_complex_taylor_series (bmad_complex_taylor, n_term, save)

implicit none

type (complex_taylor_struct) bmad_complex_taylor
type (complex_taylor_term_struct), allocatable :: term(:)
integer n_term
integer n
logical, optional :: save

!

if (n_term < 1) then
  if (associated(bmad_complex_taylor%term)) deallocate(bmad_complex_taylor%term)
  return
endif

if (.not. associated (bmad_complex_taylor%term)) then
  allocate (bmad_complex_taylor%term(n_term))
  return
endif

if (size(bmad_complex_taylor%term) == n_term) return

!

if (logic_option (.false., save) .and. n_term > 0 .and. size(bmad_complex_taylor%term) > 0) then
  n = min (n_term, size(bmad_complex_taylor%term))
  allocate (term(n))
  term = bmad_complex_taylor%term(1:n)
  deallocate (bmad_complex_taylor%term)
  allocate (bmad_complex_taylor%term(n_term))
  bmad_complex_taylor%term(1:n) = term
  deallocate (term)

else
  deallocate (bmad_complex_taylor%term)
  allocate (bmad_complex_taylor%term(n_term))
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kill_complex_taylor (bmad_complex_taylor)
!
! Subroutine to deallocate a Bmad complex_taylor map.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_complex_taylor(:) -- complex_taylor_struct: complex_taylor to be deallocated. It is OK
!                       if bmad_complex_taylor has already been deallocated.
!
! Output:
!   bmad_complex_taylor(:) -- complex_taylor_struct: deallocated complex_taylor structure.
!-

subroutine kill_complex_taylor (bmad_complex_taylor)

implicit none

type (complex_taylor_struct) bmad_complex_taylor(:)

integer i

!

do i = 1, size(bmad_complex_taylor)
  if (associated(bmad_complex_taylor(i)%term)) deallocate (bmad_complex_taylor(i)%term)
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! subroutine sort_complex_taylor_terms (complex_taylor_in, complex_taylor_sorted)
!
! Subroutine to sort the complex_taylor terms from "lowest" to "highest" of
! a complex_taylor series.
! This subroutine is needed because what comes out of PTC is not sorted.
!
! Uses function complex_taylor_exponent_index to sort.
!
! Note: complex_taylor_sorted needs to have been initialized.
! Note: complex_taylor_sorted cannot be complex_taylor_in. That is it is not legal to write:
!           call sort_complex_taylor_terms (this_complex_taylor, this_complex_taylor)
!
! Modules needed:
!   use bmad
!
! Input:
!   complex_taylor_in     -- complex_taylor_struct: Unsorted complex_taylor series.
!
! Output:
!   complex_taylor_sorted -- complex_taylor_struct: Sorted complex_taylor series.
!-

subroutine sort_complex_taylor_terms (complex_taylor_in, complex_taylor_sorted)

use nr

implicit none

type (complex_taylor_struct), intent(in)  :: complex_taylor_in
type (complex_taylor_struct) :: complex_taylor_sorted
type (complex_taylor_term_struct), allocatable :: tt(:)

integer, allocatable :: ord(:), ix(:)

integer i, n, expn(6)

! init

n = size(complex_taylor_in%term)
if (associated(complex_taylor_sorted%term)) deallocate(complex_taylor_sorted%term)
allocate(complex_taylor_sorted%term(n), ix(n), ord(n), tt(n))

!

tt = complex_taylor_in%term

do i = 1, n
  expn = tt(i)%expn
  ord(i) = complex_taylor_exponent_index(expn)
enddo

call indexx (ord, ix)

do i = 1, n
  complex_taylor_sorted%term(i)= tt(ix(i))
enddo

deallocate(ord, ix, tt)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function complex_taylor_exponent_index(expn) result(index)
!
! Function to associate a unique number with a complex_taylor exponent.
!
! The number associated with a complex_taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
! 
! Input:
!   expn(6)  -- integer: complex_taylor exponent
!
! Output:
!   index    -- integer: Sorted complex_taylor series.
!-
function complex_taylor_exponent_index(expn) result(index)
implicit none
integer :: expn(6), index
!
index = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
              expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
end function


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine complex_taylor_to_mat6 (a_complex_taylor, r_in, vec0, mat6, r_out)
!
! Subroutine to calculate, from a complex_taylor map and about some trajectory:
!   The 1st order (Jacobian) transfer matrix.
!
! Modules needed:
!   use bmad
!
! Input:
!   a_complex_taylor(6) -- complex_taylor_struct: complex_taylor map.
!   r_in(6)     -- complex(rp): Coordinates at the input. 
!
! Output:
!   vec0(6)   -- complex(rp): 0th order tranfsfer map
!   mat6(6,6) -- complex(rp): 1st order transfer map (6x6 matrix).
!   r_out(6)  -- complex(rp), optional: Coordinates at output.
!-

subroutine complex_taylor_to_mat6 (a_complex_taylor, r_in, vec0, mat6, r_out)

implicit none

type (complex_taylor_struct), target, intent(in) :: a_complex_taylor(6)
complex(rp), intent(in) :: r_in(:)
complex(rp), optional :: r_out(:)
type (complex_taylor_term_struct), pointer :: term

complex(rp), intent(out) :: mat6(6,6), vec0(6)
complex(rp) prod, t(6), t_out, out(6), r_ref

integer i, j, k, l

! mat6 calc

mat6 = 0
vec0 = 0
out = 0

do i = 1, 6
  r_ref = a_complex_taylor(i)%ref

  terms: do k = 1, size(a_complex_taylor(i)%term)
    term => a_complex_taylor(i)%term(k)
 
    t_out = term%coef
    do l = 1, 6
      if (term%expn(l) == 0) cycle
      t(l) = (r_in(l) - r_ref)**term%expn(l)
      t_out = t_out * t(l)
    enddo

    out(i) = out(i) + t_out
 
    do j = 1, 6
 
      if (term%expn(j) == 0) cycle
      if (term%expn(j) > 1 .and. r_in(j) == r_ref) cycle

      if (term%expn(j) > 1)then
        prod = term%coef * term%expn(j) * (r_in(j) - r_ref)**(term%expn(j)-1)
      else  ! term%expn(j) == 1
        prod = term%coef
      endif

      do l = 1, 6
        if (term%expn(l) == 0) cycle
        if (l == j) cycle
        prod = prod * t(l)
      enddo

      mat6(i,j) = mat6(i,j) + prod

    enddo

  enddo terms

  vec0(i) = out(i) - sum(mat6(i,:) * r_in)

enddo

if (present(r_out)) r_out = out

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine mat6_to_complex_taylor (vec0, mat6, bmad_complex_taylor)
!
! Subroutine to form a first order complex_taylor map from the 6x6 transfer
! matrix and the 0th order transfer vector.
!
! Modules needed:
!   use bmad
!
! Input:
!   vec0(6)   -- 0th order transfer vector.
!   mat6(6,6) -- 6x6 transfer matrix.
!
! Output:
!   bmad_complex_taylor(6) -- complex_taylor_struct: first order complex_taylor map.
!-

subroutine mat6_to_complex_taylor (vec0, mat6, bmad_complex_taylor)

implicit none

type (complex_taylor_struct) bmad_complex_taylor(6)

complex(rp), intent(in) :: mat6(6,6)
complex(rp), intent(in) :: vec0(6)

integer i, j, n
!

call kill_complex_taylor(bmad_complex_taylor)

do i = 1, 6
  n = count(mat6(i,1:6) /= 0)
  if (vec0(i) /= 0) n = n + 1
  allocate (bmad_complex_taylor(i)%term(n))

  n = 0

  if (vec0(i) /= 0) then
    n = n + 1
    bmad_complex_taylor(i)%term(1)%coef = vec0(i)
    bmad_complex_taylor(i)%term(1)%expn = 0
  endif

  do j = 1, 6
    if (mat6(i,j) /= 0) then
      n = n + 1
      bmad_complex_taylor(i)%term(n)%coef = mat6(i,j)
      bmad_complex_taylor(i)%term(n)%expn = 0
      bmad_complex_taylor(i)%term(n)%expn(j) = 1
    endif
  enddo

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track_complex_taylor (start_orb, bmad_complex_taylor, end_orb)
!
! Subroutine to track using a complex_taylor map.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_complex_taylor(6) -- complex_taylor_struct: complex_taylor map.
!   start_orb      -- complex(rp): Starting coords.
!
! Output:
!   end_orb        -- complex(rp): Ending coords.
!-

subroutine track_complex_taylor (start_orb, bmad_complex_taylor, end_orb)

implicit none

type (complex_taylor_struct), intent(in) :: bmad_complex_taylor(:)

complex(rp), intent(in) :: start_orb(:)
complex(rp), intent(out) :: end_orb(:)
complex(rp) diff_orb(6)
complex(rp), allocatable :: expn(:,:)

integer i, j, k, ie, e_max, i_max

!

diff_orb = start_orb - end_orb

! size cache matrix

e_max = 0
i_max = size(bmad_complex_taylor)

do i = 1, i_max
  do j = 1, size(bmad_complex_taylor(i)%term)
    e_max = max (e_max, maxval(bmad_complex_taylor(i)%term(j)%expn)) 
  enddo
enddo

allocate (expn(0:e_max, i_max))

! Fill in cache matrix

expn(0,:) = 1.0d0  !for when ie=0
expn(1,:) = start_orb(:)
do j = 2, e_max
  expn(j,:) = expn(j-1,:) * diff_orb(:)
enddo

! compute complex_taylor map

end_orb = 0

do i = 1, i_max
  do j = 1, size(bmad_complex_taylor(i)%term)
    end_orb(i) = end_orb(i) + bmad_complex_taylor(i)%term(j)%coef * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(1), 1) * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(2), 2) * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(3), 3) * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(4), 4) * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(5), 5) * &
                       expn(bmad_complex_taylor(i)%term(j)%expn(6), 6)
  enddo
enddo

end subroutine


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine truncate_complex_taylor_to_order (complex_taylor_in, order, complex_taylor_out)
!
! Subroutine to throw out all terms in a complex_taylor map that are above a certain order.
!
! Modules needed:
!   use bmad
!
! Input:
!   complex_taylor_in(:)   -- complex_taylor_struct: Input complex_taylor map.
!   order          -- Integer: Order above which terms are dropped.
!
! Output:
!   complex_taylor_out(:)  -- complex_taylor_struct: Truncated complex_taylor map.
!-

subroutine truncate_complex_taylor_to_order (complex_taylor_in, order, complex_taylor_out)

implicit none

type (complex_taylor_struct) :: complex_taylor_in(:), complex_taylor_out(:)
type (complex_taylor_struct) :: complex_taylor

integer order
integer i, j, n

!

do i = 1, size(complex_taylor_in)

  complex_taylor = complex_taylor_in(i) ! in case actual args complex_taylor_out and complex_taylor_in are the same

  ! Count the number of terms

  n = 0
  do j = 1, size(complex_taylor%term)
    if (sum(complex_taylor%term(j)%expn) <= order) n = n + 1
  enddo

  call init_complex_taylor_series (complex_taylor_out(i), n)

  n = 0
  do j = 1, size(complex_taylor%term)
    if (sum(complex_taylor%term(j)%expn) > order) cycle
    n = n + 1
    complex_taylor_out(i)%term(n) = complex_taylor%term(j)
  enddo

  complex_taylor_out(i)%ref = complex_taylor_in(i)%ref

enddo

end subroutine truncate_complex_taylor_to_order

end module
