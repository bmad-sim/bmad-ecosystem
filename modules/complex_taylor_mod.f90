!+
! Module complex_taylor_mod
!
! Note: This module is essentually a copy of taylor_mod, with:
!   taylor -> complex_taylor
!   real -> complex
!   format modifications to: type_complex_taylors 
!-

module complex_taylor_mod

use equal_mod

!

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function complex_taylor_coef (complex_taylor, exp)
! Function complex_taylor_coef (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series.
!
! Note: complex_taylor_coef is overloaded by:
!   complex_taylor_coef1 (complex_taylor, exp)
!   complex_taylor_coef2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
! Using the complex_taylor_coef2 form limits obtaining coefficients to 9th order
! or less. Also: complex_taylor_coef2 does not check that all i1, ..., i9 are between
! 1 and 6.
!
! For example: To get the 2nd order term corresponding to 
!   y(out) = Coef * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (complex_taylor_struct) complex_taylor(6)      ! complex_taylor Map
!   ...
!   coef = complex_taylor_coef (complex_taylor(3), 6, 6)  ! 1st possibility or ...
!   coef = complex_taylor_coef (complex_taylor(3), [0, 0, 0, 0, 0, 2 ])  
!
! Input (complex_taylor_coef1):
!   complex_taylor -- complex_taylor_struct: complex_taylor series.
!   exp(6)      -- Integer: Array of exponent indices.
!
! Input (complex_taylor_coef2):
!   complex_taylor -- complex_taylor_struct: complex_taylor series.
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
! Subroutine add_complex_taylor_term (complex_taylor, coef, exp)
! Subroutine add_complex_taylor_term (complex_taylor, coef, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! If complex_taylor does not have a term with the same exponents as expn then a new
! term is added to complex_taylor so the total number of terms is increased by one.
!
! If complex_taylor already has a term with the same exponents then
! the "replace" argument determines what happens:
!   If replace = False (default) then
!      coef is added to the coefficient of the old term.
!   If replace = True then
!      coef replaces the coefficient of the old term.
! In both these cases, the number of terms in complex_taylor remains the same.
!
! Note: add_complex_taylor_term is overloaded by:
!   add_complex_taylor_term1 (complex_taylor, exp, replace)
!   add_complex_taylor_term2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)
! Using the add_complex_taylor_term2 form limits obtaining coefficients to 9th order
! or less. Also: add_complex_taylor_term2 does not check that all i1, ..., i9 are between
! 1 and 6.
!
! For example: To add the 2nd order term corresponding to:
!   y(out) = 1.34 * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (complex_taylor_struct) complex_taylor(6)      ! complex_taylor Map
!   ...
!   coef = add_complex_taylor_term (complex_taylor(3), 1.34_rp, 6, 6)  ! 1st possibility or ...
!   coef = add_complex_taylor_term (complex_taylor(3), 1.34_rp, [0, 0, 0, 0, 0, 2 ])  
!
! Input (add_complex_taylor_term1):
!   complex_taylor -- complex_taylor_struct: complex_taylor series.
!   coef        -- complex(rp): Coefficient.
!   exp(6)      -- Integer: Array of exponent indices.
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Input (add_complex_taylor_term2):
!   complex_taylor -- complex_taylor_struct: complex_taylor series.
!   coef        -- complex(rp): Coefficient.
!   i1, ..., i9 -- Integer, optional: Exponent indexes (each between 1 and 6).
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Output:
!   complex_taylor -- complex_taylor_struct: New series with term added
!-

private add_complex_taylor_term1, add_complex_taylor_term2

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Elemental subroutine complex_taylor_clean (complex_taylor)
!
! Routine to "cleanup" a complex_taylor map by removing terms whose coefficients are small.
!
! Input:
!   complex_taylor -- Taylor_struct: Taylor series or array (since routine is elemental).
!
! Output:
!  complex_taylor -- Taylor_struct: Cleaned Taylor series/map
!-

elemental subroutine complex_taylor_clean (complex_taylor)

implicit none

type (complex_taylor_struct), intent(inout) :: complex_taylor
type (complex_taylor_struct) t2
real(rp), parameter :: eps = 1d-10, ps_bound = 0.1
integer i, n


! ps_bound is an approximate upper bound to the phase space coordinates

if (.not. associated(complex_taylor%term)) return

n = 0 
do i = 1, size(complex_taylor%term)
  if (abs(complex_taylor%term(i)%coef) * ps_bound**sum(complex_taylor%term(i)%expn) > eps) n = n + 1
enddo

if (n == size(complex_taylor%term)) return

t2%term => complex_taylor%term
allocate(complex_taylor%term(n))

n = 0 
do i = 1, size(t2%term)
  if (abs(t2%term(i)%coef) * ps_bound**sum(t2%term(i)%expn) > eps) then
    n = n + 1
    complex_taylor%term(n) = t2%term(i)
  endif
enddo

deallocate(t2%term)

end subroutine complex_taylor_clean

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function complex_taylor_coef1 (complex_taylor, exp)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series. This routine is used by the overloaded function
! complex_taylor_coef. See complex_taylor_coef for more details.
!-

function complex_taylor_coef1 (complex_taylor, exp) result (coef)

implicit none

type (complex_taylor_struct), intent(in) :: complex_taylor

complex(rp) coef

integer, intent(in) :: exp(:)
integer i

!

coef = 0

do i = 1, size(complex_taylor%term)
  if (all(complex_taylor%term(i)%expn == exp)) then
    coef = complex_taylor%term(i)%coef
    return
  endif
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function complex_taylor_coef2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
!
! Function to return the coefficient for a particular complex_taylor term
! from a complex_taylor Series. This routine is used by the overloaded function
! complex_taylor_coef. See complex_taylor_coef for more details.
!-

function complex_taylor_coef2 (complex_taylor, i1, i2, i3, &
                          i4, i5, i6, i7, i8, i9) result (coef)

implicit none

type (complex_taylor_struct), intent(in) :: complex_taylor

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

do i = 1, size(complex_taylor%term)
  if (all(complex_taylor%term(i)%expn == exp)) then
    coef = complex_taylor%term(i)%coef
    return
  endif
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine type_complex_taylors (complex_taylor, max_order, lines, n_lines, file_id, out_type, clean)
!
! Subroutine to print or put in a string array a Bmad taylor map.
! If the lines(:) argument is not present, the element information is printed to the terminal.
!
! Input:
!   complex_taylor(:) -- complex_taylor_struct: Array of taylors.
!   max_order      -- integer, optional: Maximum order to print.
!   file_id        -- integer, optional: If present, write output to a file with handle file_id.
!   out_type       -- character(*), optional: Determins the string to be used for the output type column.
!                       '' (default)  -> '1', '2', '3', etc.
!                       'PHASE'       -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz'
!                       'SPIN'        -> 'S1', 'Sx', 'Sy', 'Sz' (quaternion representation)
!                       'NONE'        -> No out column
!                       Anything else -> Use this for the output column.
!   clean          -- logical, optional: If True then do not include terms whose coefficients
!                       are negligible. Default is false
!
! Output:
!   lines(:)     -- character(*), allocatable, optional :: Character array to hold the output. 
!                     If not present, the information is printed to the terminal. 
!                     Char width should be 120 or above for out_type = 'PHASE' but can be less for other out_types.
!   n_lines      -- integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_complex_taylors (complex_taylor, max_order, lines, n_lines, file_id, out_type, clean)

use re_allocate_mod

implicit none

type (complex_taylor_struct), intent(in), target :: complex_taylor(:)
type (complex_taylor_term_struct), pointer :: tt
type (complex_taylor_struct) tlr, t2

integer, optional, intent(out) :: n_lines
integer, optional :: max_order, file_id
integer i, j, k, nl, ix, nt

logical, optional :: clean

character(*), optional :: out_type
character(*), optional, allocatable :: lines(:)
character(120), allocatable :: li(:)
character(40) fmt1, fmt2, fmt, o_type
character(20) out_str
character(2), parameter :: spin_out(4) = ['S1', 'Sx', 'Sy', 'Sz']
character(2), parameter :: phase_out(6) = ['X ', 'Px', 'Y ', 'Py', 'Z ', 'Pz']

! If not allocated then not much to do

nt = size(complex_taylor)

if (.not. associated(complex_taylor(1)%term)) then
  nl = 2
  allocate (li(nl))
  li(1) = '---------------------------------------------------'
  li(2) = 'A Taylor Map Does Not Exist.' 

! Normal case

else

  fmt1 = '(a, 2f20.12,  1x, 6i3, i8, 2f17.9)'
  fmt2 = '(a, 2es20.11, 1x, 6i3, i8, 2f17.9)'

  o_type = ''
  if (present(out_type)) o_type = out_type

  nl = 20 + nt + sum( [(size(complex_taylor(i)%term), i = 1, nt) ])
  allocate(li(nl))

  nl = 0

  select case (o_type)
  case ('PHASE')
    nl=nl+1; write (li(nl), '(a)') ' Out      Re Coef             Im Coef         Exponents          Order      Reference'
    nl=nl+1; li(nl) =              ' --------------------------------------------------------------------------------------------------------'
  case ('NONE')
    nl=nl+1; write (li(nl), '(a)') '       Re Coef             Im Coef         Exponents          Order'
    nl=nl+1; li(nl) =              '   ----------------------------------------------------------------'
    fmt1 = '(' // fmt1(5:);  fmt2 = '(' // fmt2(5:)
  case default
    nl=nl+1; write (li(nl), '(a)') ' Out      Re Coef             Im Coef         Exponents          Order'
    nl=nl+1; li(nl) =              ' ---------------------------------------------------------------------'
  end select

  do i = 1, nt

    out_str = ' ??:'
    select case(o_type)
    case ('')
      write (out_str, '(i3, a)') i, ':'
    case ('PHASE')
      if (i <=6) out_str = ' ' // phase_out(i) // ':'
    case ('SPIN')
      if (i <=4) out_str = ' ' // spin_out(i) // ':'
    case default
      out_str = o_type // ':'
    end select

    if (logic_option(.false., clean)) then
      t2 = complex_taylor(i)
      call complex_taylor_clean(t2)
    else
      t2%term => complex_taylor(i)%term
    endif

    if (size(t2%term) == 0) then
      nl=nl+1; write (li(nl), '(a, 6x, a)') trim(out_str), 'No Terms. Always evaluates to zero.'

    else
      nullify (tlr%term)
      call sort_complex_taylor_terms (t2, tlr)

      do j = 1, size(t2%term)
        tt => tlr%term(j)

        if (present(max_order)) then
          if (sum(tt%expn) > max_order) cycle
        endif

        if (abs(tt%coef) < 1d5) then
          fmt = fmt1
        else
          fmt = fmt2
        endif

        if (j == 1 .and. o_type == 'PHASE') then
          nl=nl+1; write (li(nl), fmt) trim(out_str), tt%coef, (tt%expn(k), k = 1, 6), sum(tt%expn), complex_taylor(i)%ref
        elseif (o_type == 'NONE') then
          nl=nl+1; write (li(nl), fmt) tt%coef, (tt%expn(k), k = 1, 6), sum(tt%expn)
        else
          nl=nl+1; write (li(nl), fmt) trim(out_str), tt%coef, (tt%expn(k), k = 1, 6), sum(tt%expn)
        endif
      enddo

      deallocate (tlr%term)
    endif

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

if (present(file_id)) then
  do i = 1, nl
    write (file_id, '(a)') trim(li(i))
  enddo
endif

end subroutine type_complex_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine complex_taylor_make_unit (complex_taylor)
!
! Subroutine to make the unit complex_taylor map:
!       r(out) = Map * r(in) = r(in)
!
! Output:
!   complex_taylor(6) -- complex_taylor_struct: Unit complex_taylor map .
!-

subroutine complex_taylor_make_unit (complex_taylor)

implicit none

type (complex_taylor_struct) complex_taylor(:)
integer i

do i = 1, size(complex_taylor)
  call init_complex_taylor_series (complex_taylor(i), 1)
  complex_taylor(i)%term(1)%coef = 1.0
  complex_taylor(i)%term(1)%expn = 0
  complex_taylor(i)%term(1)%expn(i) = 1
  complex_taylor(i)%ref = 0
enddo

end subroutine complex_taylor_make_unit

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_complex_taylor_term1 (complex_taylor, coef, expn, replace)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! This routine is used by the overloaded function add_complex_taylor_term. 
! See add_complex_taylor_term for more details.
!-

subroutine add_complex_taylor_term1 (complex_taylor, coef, expn, replace)

implicit none

type (complex_taylor_struct) complex_taylor

complex(rp) coef
integer expn(:), i, n

logical, optional :: replace

! Search for an existing term of the same type

n = size(complex_taylor%term)

do i = 1, n
  if (all(complex_taylor%term(i)%expn == expn)) then
    if (logic_option(.false., replace)) then
      complex_taylor%term(i)%coef = coef
    else
      complex_taylor%term(i)%coef = coef + complex_taylor%term(i)%coef
    endif
    if (complex_taylor%term(i)%coef == 0) then  ! Kill this term
      complex_taylor%term(i:n-1) = complex_taylor%term(i+1:n)
      call init_complex_taylor_series (complex_taylor, n-1, .true.)
    endif
    return
  endif
enddo

! new term

call init_complex_taylor_series (complex_taylor, n+1, .true.)
complex_taylor%term(n+1)%coef = coef
complex_taylor%term(n+1)%expn = expn

end subroutine add_complex_taylor_term1

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_complex_taylor_term2 (complex_taylor, coef, i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)
!
! Routine to add a complex_taylor term to a complex_taylor series.
!
! This routine is used by the overloaded function add_complex_taylor_term. 
! See add_complex_taylor_term for more details.
!-

subroutine add_complex_taylor_term2 (complex_taylor, coef, i1, i2, i3, i4, i5, i6, i7, i8, i9, replace)

implicit none

type (complex_taylor_struct) complex_taylor

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

call add_complex_taylor_term1 (complex_taylor, coef, expn, replace)

end subroutine add_complex_taylor_term2

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kill_complex_taylor (complex_taylor)
!
! Subroutine to deallocate a Bmad complex_taylor map.
!
! Input:
!   complex_taylor(:) -- complex_taylor_struct: complex_taylor to be deallocated. It is OK
!                       if complex_taylor has already been deallocated.
!
! Output:
!   complex_taylor(:) -- complex_taylor_struct: deallocated complex_taylor structure.
!-

subroutine kill_complex_taylor (complex_taylor)

implicit none

type (complex_taylor_struct) complex_taylor(:)

integer i

!

do i = 1, size(complex_taylor)
  if (associated(complex_taylor(i)%term)) deallocate (complex_taylor(i)%term)
enddo

end subroutine kill_complex_taylor

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
! Input:
!   complex_taylor_in     -- complex_taylor_struct: Unsorted complex_taylor series.
!
! Output:
!   complex_taylor_sorted -- complex_taylor_struct: Sorted complex_taylor series.
!-

subroutine sort_complex_taylor_terms (complex_taylor_in, complex_taylor_sorted)


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

call indexer (ord, ix)

do i = 1, n
  complex_taylor_sorted%term(i)= tt(ix(i))
enddo

deallocate(ord, ix, tt)

end subroutine sort_complex_taylor_terms

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

end function complex_taylor_exponent_index

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine complex_taylor_to_mat6 (a_complex_taylor, r_in, vec0, mat6, r_out)
!
! Subroutine to calculate, from a complex_taylor map and about some trajectory:
!   The 1st order (Jacobian) transfer matrix.
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

end subroutine complex_taylor_to_mat6

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine mat6_to_complex_taylor (vec0, mat6, complex_taylor)
!
! Subroutine to form a first order complex_taylor map from the 6x6 transfer
! matrix and the 0th order transfer vector.
!
! Input:
!   vec0(6)   -- 0th order transfer vector.
!   mat6(6,6) -- 6x6 transfer matrix.
!
! Output:
!   complex_taylor(6) -- complex_taylor_struct: first order complex_taylor map.
!-

subroutine mat6_to_complex_taylor (vec0, mat6, complex_taylor)

implicit none

type (complex_taylor_struct) complex_taylor(6)

complex(rp), intent(in) :: mat6(6,6)
complex(rp), intent(in) :: vec0(6)

integer i, j, n
!

call kill_complex_taylor(complex_taylor)

do i = 1, 6
  n = count(mat6(i,1:6) /= 0)
  if (vec0(i) /= 0) n = n + 1
  allocate (complex_taylor(i)%term(n))

  n = 0

  if (vec0(i) /= 0) then
    n = n + 1
    complex_taylor(i)%term(1)%coef = vec0(i)
    complex_taylor(i)%term(1)%expn = 0
  endif

  do j = 1, 6
    if (mat6(i,j) /= 0) then
      n = n + 1
      complex_taylor(i)%term(n)%coef = mat6(i,j)
      complex_taylor(i)%term(n)%expn = 0
      complex_taylor(i)%term(n)%expn(j) = 1
    endif
  enddo

enddo

end subroutine mat6_to_complex_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track_complex_taylor (start_orb, complex_taylor, end_orb)
!
! Subroutine to track using a complex_taylor map.
!
! Input:
!   complex_taylor(6) -- complex_taylor_struct: complex_taylor map.
!   start_orb      -- complex(rp): Starting coords.
!
! Output:
!   end_orb        -- complex(rp): Ending coords.
!-

subroutine track_complex_taylor (start_orb, complex_taylor, end_orb)

implicit none

type (complex_taylor_struct), intent(in) :: complex_taylor(:)

complex(rp), intent(in) :: start_orb(:)
complex(rp), intent(out) :: end_orb(:)
complex(rp) diff_orb(6)
complex(rp), allocatable :: expn(:,:)

integer i, j, k, ie, e_max, i_max

!

diff_orb = start_orb - end_orb

! size cache matrix

e_max = 0
i_max = size(complex_taylor)

do i = 1, i_max
  do j = 1, size(complex_taylor(i)%term)
    e_max = max (e_max, maxval(complex_taylor(i)%term(j)%expn)) 
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
  do j = 1, size(complex_taylor(i)%term)
    end_orb(i) = end_orb(i) + complex_taylor(i)%term(j)%coef * &
                       expn(complex_taylor(i)%term(j)%expn(1), 1) * &
                       expn(complex_taylor(i)%term(j)%expn(2), 2) * &
                       expn(complex_taylor(i)%term(j)%expn(3), 3) * &
                       expn(complex_taylor(i)%term(j)%expn(4), 4) * &
                       expn(complex_taylor(i)%term(j)%expn(5), 5) * &
                       expn(complex_taylor(i)%term(j)%expn(6), 6)
  enddo
enddo

end subroutine track_complex_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine truncate_complex_taylor_to_order (complex_taylor_in, order, complex_taylor_out)
!
! Subroutine to throw out all terms in a complex_taylor map that are above a certain order.
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
