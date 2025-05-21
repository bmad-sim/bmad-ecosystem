!+
! Module taylor_mod
!
! Note: The companion complex_taylor_mod is the same as this one, with 
!   taylor -> complex_taylor
!   real -> complex
! When editing this file, please update bmad_complex_taylor_mod
!-

module taylor_mod

use bmad_routine_interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Elemental subroutine taylor_clean (bmad_taylor)
!
! Routine to "cleanup" a taylor map by removing terms whose coefficients are small.
!
! Input:
!   bmad_taylor -- Taylor_struct: Taylor series or array (since routine is elemental).
!
! Output:
!    bmad_taylor -- Taylor_struct: Cleaned Taylor series/map
!-

elemental subroutine taylor_clean (bmad_taylor)

implicit none

type (taylor_struct), intent(inout) :: bmad_taylor
type (taylor_struct) t2
real(rp), parameter :: eps = 1d-12, ps_bound = 0.2
integer i, n


! ps_bound is an approximate upper bound to the phase space coordinates

if (.not. associated(bmad_taylor%term)) return

n = 0 
do i = 1, size(bmad_taylor%term)
  if (abs(bmad_taylor%term(i)%coef) * ps_bound**sum(bmad_taylor%term(i)%expn) > eps) n = n + 1
enddo

if (n == size(bmad_taylor%term)) return

t2%term => bmad_taylor%term
allocate(bmad_taylor%term(n))

n = 0 
do i = 1, size(t2%term)
  if (abs(t2%term(i)%coef) * ps_bound**sum(t2%term(i)%expn) > eps) then
    n = n + 1
    bmad_taylor%term(n) = t2%term(i)
  endif
enddo

deallocate(t2%term)

end subroutine taylor_clean

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function taylor_coef (bmad_taylor, expn) result (taylor_coef)
!
! Function to return the coefficient for a particular taylor term from a Taylor Series.
!
! For example: To get the 2nd order term corresponding to 
!   y(out) = Coef * p_z^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (taylor_struct) bmad_taylor(6)      ! Taylor Map
!   ...
!   coef = taylor_coef (bmad_taylor(3), [0, 0, 0, 0, 0, 2 ]) or  
!   coef = taylor_coef (bmad_taylor(3), taylor_expn([6,6]))
!
! Input (taylor_coef1):
!   bmad_taylor -- Taylor_struct: Taylor series.
!   exp(6)      -- Integer: Array of exponent indices.
!
! Output:
!   taylor_coef -- Real(rp): Coefficient. Set to zero if not found.
!-

function taylor_coef (bmad_taylor, expn) result (coef)

implicit none

type (taylor_struct) :: bmad_taylor

real(rp) coef

integer, intent(in) :: expn(:)
integer i

!

coef = 0

do i = 1, size(bmad_taylor%term)
  if (all(bmad_taylor%term(i)%expn == expn)) then
    coef = bmad_taylor%term(i)%coef
    return
  endif
enddo

end function taylor_coef

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function taylor_expn (coord_array) result (expn_array)
!
! Function to return the exponantial 6-vector for a taylor term given the "coordinate array".
! For example:
!   taylor_expn([3,4]) = [0, 0, 1, 1, 0, 0]
!   taylor_expn([2,5,5]) = [0, 1, 0, 0, 2, 0]
!
! Input:
!   coord_array(:)    -- integer, series of coordinates.
!
! Output:
!   expn_array(6)     -- integer, array of exponents.
!-

function taylor_expn (coord_array) result (expn_array)

implicit none

integer coord_array(:), expn_array(6)
integer i, j, k

!

expn_array = 0

do i = 1, size(coord_array)
  j = coord_array(i)
  expn_array(j) = expn_array(j) + 1
enddo

end function taylor_expn

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine taylor_make_quaternion_unit (bmad_taylor)
!
! Subroutine to make the "unit" quaternion Taylor map. Used for spin maps.
!
! Output:
!   bmad_taylor(0:3)  -- taylor_struct: Unit quaternion Taylor map.
!-

subroutine taylor_make_quaternion_unit (bmad_taylor)

implicit none

type (taylor_struct) bmad_taylor(0:3)
integer i

!

do i = 0, 3
  call init_taylor_series (bmad_taylor(i), 1)
  if (i == 0) bmad_taylor(i)%term(1)%coef = 1.0
  bmad_taylor(i)%term(1)%expn = 0
enddo

end subroutine taylor_make_quaternion_unit

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine taylor_make_unit (bmad_taylor, ref_orbit)
!
! Subroutine to make the "unit" Taylor map around the reference orbit:
!   r(out) = Map * dr(in) = dr(in) + ref_orbit = r(in) 
! where
!   dr(in) = r(in) - ref_orbit
!
! Input:
!   ref_orbit(6)    -- real(rp), optional: Reference orbit. Taken to be zero if not present.
!
! Output:
!   bmad_taylor(6)  -- taylor_struct: Unit Taylor map.
!     %ref            -- Set to the reference orbit.
!-

subroutine taylor_make_unit (bmad_taylor, ref_orbit)

implicit none

type (taylor_struct) bmad_taylor(:)
real(rp), optional :: ref_orbit(:)
integer i

!

if (present(ref_orbit)) then
  bmad_taylor%ref = ref_orbit
else
  bmad_taylor%ref = 0
endif

do i = 1, size(bmad_taylor)
  if (bmad_taylor(i)%ref == 0) then
    call init_taylor_series (bmad_taylor(i), 1)
    bmad_taylor(i)%term(1)%coef = 1.0
    bmad_taylor(i)%term(1)%expn = 0
    bmad_taylor(i)%term(1)%expn(i) = 1
  else
    call init_taylor_series (bmad_taylor(i), 2, .true.)
    bmad_taylor(i)%term(1)%coef = bmad_taylor(i)%ref
    bmad_taylor(i)%term(1)%expn = 0
    bmad_taylor(i)%term(2)%coef = 1.0
    bmad_taylor(i)%term(2)%expn = 0
    bmad_taylor(i)%term(2)%expn(i) = 1
  endif
enddo

end subroutine taylor_make_unit

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function taylor_term_index (bmad_taylor, expn, create) result (ix_term)
!
! Routine to return the index of a particular taylor term.
!
! Input:
!   bmad_taylor   -- taylor_struct: Taylor series to search.
!   expn(6)       -- integer: Exponents to match to.
!   create        -- logical: If True, crate taylor term if it does not exist.
!
! Output:
!   ix_term       -- integer: Index of taylor term in bmad_taylor%term(:). Will be 0 if 
!                      create = False and the term does not exist.
!-

function taylor_term_index (bmad_taylor, expn, create) result (ix_term)

implicit none

type (taylor_struct) bmad_taylor

integer expn(6), ix_term
integer i, n

logical create 

! First search for existing term.

if (associated(bmad_taylor%term)) then
  n = size(bmad_taylor%term)
else
  n = 0
endif

do i = 1, n
  if (.not. all(bmad_taylor%term(i)%expn == expn)) cycle
  ix_term = i
  return
enddo

! Not found.

if (create) then
  call init_taylor_series (bmad_taylor, n+1, .true.)
  bmad_taylor%term(n+1)%expn = expn
  ix_term = n+1
else
  ix_term = 0
endif

end function taylor_term_index

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_taylor_term (bmad_taylor, coef, expn, replace)
!
! Routine to add a Taylor term to a Taylor series.
!
! If bmad_taylor does not have a term with the same exponents then a new
! term is added to bmad_taylor so the total number of terms is increased by one.
!
! If bmad_taylor already has a term with the same exponents then
! the "replace" argument determines what happens:
!   If replace = False (default) then
!      coef is added to the coefficient of the old term.
!   If replace = True then
!      coef replaces the coefficient of the old term.
! In both these cases, the number of terms in bmad_taylor remains the same.
!
! For example: To add the 2nd order term corresponding to:
!   y(out) = 1.34 * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (taylor_struct) bmad_taylor(6)      ! Taylor Map
!   ...
!   coef = add_taylor_term (bmad_taylor(3), 1.34_rp, [0, 0, 0, 0, 0, 2 ])    ! or
!   coef = add_taylor_term (bmad_taylor(3), 1.34_rp, taylor_expn([6,6]))  
!
! Input:
!   bmad_taylor -- Taylor_struct: Taylor series.
!   coef        -- Real(rp): Coefficient.
!   expn(6)     -- Integer: Array of exponent indices.
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Output:
!   bmad_taylor -- Taylor_struct: New series with term added.
!-

subroutine add_taylor_term (bmad_taylor, coef, expn, replace)

implicit none

type (taylor_struct) bmad_taylor

real(rp) coef
integer expn(:), n, nn

logical, optional :: replace

! 

n = taylor_term_index(bmad_taylor, expn, .true.)

if (logic_option(.false., replace)) then
  bmad_taylor%term(n)%coef = coef
else
  bmad_taylor%term(n)%coef = coef + bmad_taylor%term(n)%coef
endif

nn = size(bmad_taylor%term)
if (bmad_taylor%term(n)%coef == 0) then  ! Kill this term
  bmad_taylor%term(n:nn-1) = bmad_taylor%term(n+1:nn)
  call init_taylor_series (bmad_taylor, nn-1, .true.)
endif

end subroutine add_taylor_term

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine remove_taylor_term (bmad_taylor, expn)
!
! Routine to remove a Taylor term to a Taylor series.
!
! If bmad_taylor does not have a term with the same exponents then nothing is done.
!
! For example: To remove the 2nd order term corresponding to:
!   y(out) = xxx * p_z(in)^2 
! [This is somtimes refered to as the T_366 term]
! The call would be:
!   type (taylor_struct) bmad_taylor(6)      ! Taylor Map
!   ...
!   coef = remove_taylor_term (bmad_taylor(3), [0, 0, 0, 0, 0, 2 ])  ! or
!   coef = remove_taylor_term (bmad_taylor(3), taylor_expn([6, 6]))  
!
! Input
!   bmad_taylor -- Taylor_struct: Taylor series.
!   expn(6)     -- Integer: Array of exponent indices.
!
! Output:
!   bmad_taylor -- Taylor_struct: New series with term removed
!-
!-

subroutine remove_taylor_term (bmad_taylor, expn)

implicit none

type (taylor_struct) bmad_taylor
type (taylor_term_struct), pointer :: term(:)

integer expn(:), i, j, k, n

! Search for an existing term of the same type

n = size(bmad_taylor%term)

do i = 1, n
  if (all(bmad_taylor%term(i)%expn /= expn)) cycle
  term => bmad_taylor%term
  allocate (bmad_taylor%term(n-1))
  k = 0
  do j = 1, n
    if (j == i) cycle
    k = k + 1
    bmad_taylor%term(k) = term(j)
  enddo
  deallocate(term)
  return
enddo

end subroutine remove_taylor_term

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine taylor_extract_zeroth_order_part (bmad_taylor, zeroth_part)
!
! Routine to return the zeroth order ("constant") part of a taylor map.
!
! Input:
!   bmad_taylor(:)    -- taylor_struct: Taylor map.
!
! Output:
!   zeroth_part(:)    -- real(rp): vector of constant terms in the map.
!-

subroutine taylor_extract_zeroth_order_part (bmad_taylor, zeroth_part)

implicit none

type (taylor_struct) bmad_taylor(:)
real(rp) zeroth_part(:)
integer i, j

!

zeroth_part = 0

do i = 1, size(bmad_taylor)
  if (.not. associated(bmad_taylor(i)%term)) cycle

  do j = 1, size(bmad_taylor(i)%term)
    if (any(bmad_taylor(i)%term(j)%expn /= 0)) cycle
    zeroth_part(i) = bmad_taylor(i)%term(j)%coef
    exit
  enddo
enddo

end subroutine taylor_extract_zeroth_order_part

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine sort_taylor_terms (taylor_in, taylor_sorted, min_val)
!
! Subroutine to sort the taylor terms from "lowest" to "highest" of
! a taylor series.
! This subroutine is needed because what comes out of PTC is not sorted.
!
! Uses function taylor_exponent_index to sort.
!
! Note: It is OK for taylor_sorted to be the same actual argument as taylor_in. That is, this is OK:
!           call sort_taylor_terms (this_taylor, this_taylor)
!
! Input:
!   taylor_in     -- Taylor_struct: Unsorted taylor series.
!   min_val       -- real(rp), optional: If present then any terms lower than this value will be dropped.
!
! Output:
!   taylor_sorted -- Taylor_struct: Sorted taylor series.
!-

subroutine sort_taylor_terms (taylor_in, taylor_sorted, min_val)


implicit none

type (taylor_struct), intent(in)  :: taylor_in
type (taylor_struct) :: taylor_sorted
type (taylor_term_struct), allocatable :: tt(:)

real(rp), optional :: min_val

integer, allocatable :: ord(:), ix(:)

integer i, n, j

! init

n = size(taylor_in%term)
if (associated(taylor_sorted%term)) deallocate(taylor_sorted%term)
allocate(taylor_sorted%term(n), ix(n), ord(n), tt(n))

!

tt = taylor_in%term

do i = 1, n
  ord(i) = taylor_exponent_index(tt(i)%expn)
enddo

call indexer (ord, ix)

j = 0
do i = 1, n
  if (present(min_val)) then
    if (abs(tt(ix(i))%coef) <= min_val) cycle
  endif
  j = j + 1
  taylor_sorted%term(j)= tt(ix(i))
enddo

if (j < n) call init_taylor_series (taylor_sorted, j, .true.)
taylor_sorted%ref = taylor_in%ref

end subroutine sort_taylor_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function taylor_exponent_index (expn) result(index)
!
! Function to associate a unique number with a taylor exponent.
!
! The number associated with a taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
! 
! Input:
!   expn(6)  -- integer: taylor exponent
!
! Output:
!   index    -- integer: Sorted taylor series.
!-

function taylor_exponent_index (expn) result(index)

implicit none

integer :: expn(6), index

!
index = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
              expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)

end function taylor_exponent_index

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine taylor_to_mat6 (a_taylor, r_in, vec0, mat6, r_out)
!
! Routine to calculate from a Taylor map the 1st order (Jacobian) transfer matrix.
! Output transfer map is:
!   x_out = mat6 * x_in + vec0
! Where x_out and x_in are phase space coordinates. Notice that vec0 is *not* the
! zeroth order part of the map around r_out.
!
! Input:
!   a_taylor(6) -- Taylor_struct: Taylor map.
!   r_in(6)     -- Real(rp): Coordinates at the input about which the Jacobian is evaluated. 
!
! Output:
!   vec0(6)   -- Real(rp): See above
!   mat6(6,6) -- Real(rp): 1st order transfer map (6x6 matrix).
!   r_out(6)  -- Real(rp), optional: Coordinates at output with r_in as the starting coords.
!-

subroutine taylor_to_mat6 (a_taylor, r_in, vec0, mat6, r_out)

implicit none

type (taylor_struct), target, intent(in) :: a_taylor(6)
real(rp), intent(in) :: r_in(:)
real(rp), optional :: r_out(:)
type (taylor_term_struct), pointer :: term

real(rp), intent(out) :: mat6(6,6), vec0(6)
real(rp) prod, t(6), t_out, out(6), r_diff(6)

integer i, j, k, l

! mat6 calc

mat6 = 0
vec0 = 0
out = 0

r_diff = r_in - a_taylor(:)%ref

!

do i = 1, 6
  terms: do k = 1, size(a_taylor(i)%term)
    term => a_taylor(i)%term(k)
 
    t_out = term%coef
    do l = 1, 6
      if (term%expn(l) == 0) cycle
      t(l) = r_diff(l)**term%expn(l)
      t_out = t_out * t(l)
    enddo
    out(i) = out(i) + t_out
 
    do j = 1, 6
      if (term%expn(j) == 0) cycle
      if (term%expn(j) > 1 .and. r_diff(j) == 0) cycle

      if (term%expn(j) > 1)then
        prod = term%coef * term%expn(j) * r_diff(j)**(term%expn(j)-1)
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

end subroutine taylor_to_mat6

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine mat6_to_taylor (vec0, mat6, bmad_taylor, ref_orb)
!
! Subroutine to form a first order Taylor map from the 6x6 transfer
! matrix and the 0th order transfer vector.
!
! Input:
!   vec0(6)     -- real(rp): 0th order transfer vector.
!   mat6(6,6)   -- real(rp): 6x6 transfer matrix.
!   ref_orb(6)  -- real(rp), optional: Reference orbit at entrance to map. 
!                   Default is zero orbit.
!
! Output:
!   bmad_taylor(6) -- taylor_struct: first order taylor map.
!-

subroutine mat6_to_taylor (vec0, mat6, bmad_taylor, ref_orb)

implicit none

type (taylor_struct) bmad_taylor(6)

real(rp), intent(in) :: mat6(6,6)
real(rp), intent(in) :: vec0(6)
real(rp), optional :: ref_orb(6)

integer i, j, n

!

call kill_taylor(bmad_taylor)
if (present(ref_orb)) bmad_taylor%ref = ref_orb

do i = 1, 6
  n = count(mat6(i,1:6) /= 0)
  if (vec0(i) /= 0) n = n + 1
  allocate (bmad_taylor(i)%term(n))

  n = 0

  if (vec0(i) /= 0) then
    n = n + 1
    bmad_taylor(i)%term(1)%coef = vec0(i)
    bmad_taylor(i)%term(1)%expn = 0
  endif

  do j = 1, 6
    if (mat6(i,j) /= 0) then
      n = n + 1
      bmad_taylor(i)%term(n)%coef = mat6(i,j)
      bmad_taylor(i)%term(n)%expn = 0
      bmad_taylor(i)%term(n)%expn(j) = 1
    endif
  enddo
enddo

end subroutine mat6_to_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function track_taylor (start_orb, bmad_taylor, ref_orb) result (end_track)
!
! Routine to track using a Taylor map.
!
! This routine can be used for both orbital phase space and (spin) quaternion tracking.
! In the case of quaternion tracking, the result is the rotation rotation map which
! then must be applied to the starting spin coords. Note: In general, the quaternion
! returned by this routine will not have norm = 1 so it must be properly normalized.
!
! Input:
!   start_orb(6)   -- real(rp): Starting **phase space** coords.
!   bmad_taylor(:) -- taylor_struct: Taylor map. Array size is 6 for phase space 
!                       tracking and 4 for quaternion spin tracking.
!   ref_orb(6)     -- real(rp), optional: Phase space reference orbit. 
!                       Needed if doing spin tracking since in this case the bmad_taylor 
!                       argument will not contain the reference orbit.
!
! Output:
!   end_track(:)   -- real(rp): Ending coords for phase space tracking and a quaternion 
!                       for quaternion tracking. Array size will be 6 for phase space 
!                       tracking and 4 for quaternion tracking (same as bmad_taylor(:)).
!-

function track_taylor (start_orb, bmad_taylor, ref_orb) result (end_track)

implicit none

type (taylor_struct) :: bmad_taylor(:)

real(rp) :: start_orb(:)
real(rp) :: end_track(size(bmad_taylor))
real(rp), allocatable :: expn(:, :)
real(rp), optional :: ref_orb(:)
real(rp) diff_orb(size(start_orb))

integer i, j, k, ie, e_max, n_size

! 

if (present(ref_orb)) then
  diff_orb = start_orb - ref_orb
else
  diff_orb = start_orb - bmad_taylor%ref
endif

! size cache matrix

e_max = 0
n_size = size(end_track)

do i = 1, n_size
  do j = 1, size(bmad_taylor(i)%term)
    e_max = max (e_max, maxval(bmad_taylor(i)%term(j)%expn)) 
  enddo
enddo

allocate (expn(0:e_max, 6))

! Fill in cache matrix

expn(0,:) = 1.0d0
do j = 1, e_max
  expn(j,:) = expn(j-1,:) * diff_orb(:)
enddo

! Compute taylor map

end_track = 0

do i = 1, n_size
  do j = 1, size(bmad_taylor(i)%term)
    end_track(i) = end_track(i) + bmad_taylor(i)%term(j)%coef * &
                       expn(bmad_taylor(i)%term(j)%expn(1), 1) * &
                       expn(bmad_taylor(i)%term(j)%expn(2), 2) * &
                       expn(bmad_taylor(i)%term(j)%expn(3), 3) * &
                       expn(bmad_taylor(i)%term(j)%expn(4), 4) * &
                       expn(bmad_taylor(i)%term(j)%expn(5), 5) * &
                       expn(bmad_taylor(i)%term(j)%expn(6), 6)
  enddo
enddo

end function track_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine truncate_taylor_to_order (taylor_in, order, taylor_out)
!
! Subroutine to throw out all terms in a taylor map that are above a certain order.
!
! Input:
!   taylor_in(:)   -- Taylor_struct: Input Taylor map.
!   order          -- Integer: Order above which terms are dropped.
!
! Output:
!   taylor_out(:)  -- Taylor_struct: Truncated Taylor map.
!-

subroutine truncate_taylor_to_order (taylor_in, order, taylor_out)

implicit none

type (taylor_struct) :: taylor_in(:), taylor_out(:)
type (taylor_struct) :: taylor

integer order
integer i, j, n

!

do i = 1, size(taylor_in)

  taylor = taylor_in(i) ! in case actual args taylor_out and taylor_in are the same

  ! Count the number of terms

  n = 0
  do j = 1, size(taylor%term)
    if (sum(taylor%term(j)%expn) <= order) n = n + 1
  enddo

  call init_taylor_series (taylor_out(i), n, .true.)

  n = 0
  do j = 1, size(taylor%term)
    if (sum(taylor%term(j)%expn) > order) cycle
    n = n + 1
    taylor_out(i)%term(n) = taylor%term(j)
  enddo

  taylor_out(i)%ref = taylor_in(i)%ref

enddo

end subroutine truncate_taylor_to_order

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine evaluate_em_taylor (r_pos, em_taylor, field, dfield)
!
! Routine to evaluate the field using an em_taylor map.
! Note: dField/dz values will be calculated using Maxwell's Equations for a static field.
!
! Input:
!   r_pos(2)     -- real(rp): (x, y) position
!   em_taylor(3) -- em_taylor_struct: em_taylor map.
!
! Output:
!   field(3)     -- real(rp): Field.
!   dfield(3,3)  -- real(rp), optional: Field derivatives.
!-

subroutine evaluate_em_taylor (r_pos, em_taylor, field, dfield)

implicit none

type (em_taylor_struct) :: em_taylor(3)

real(rp) r_pos(2), field(3)
real(rp), optional :: dfield(3,3)

real(rp), allocatable :: expn(:, :)

integer i, j, ie_max, iex, iey

! size cache matrix

ie_max = 0

do i = 1, 3
  do j = 1, size(em_taylor(i)%term)
    ie_max = max (ie_max, maxval(em_taylor(i)%term(j)%expn)) 
  enddo
enddo

allocate (expn(0:ie_max, 2))

! Fill in cache matrix

expn(0,:) = 1

do j = 1, ie_max
  expn(j,:) = expn(j-1,:) * r_pos(:)
enddo

! Compute taylor map

field = 0

do i = 1, 3
  do j = 1, size(em_taylor(i)%term)
    field(i) = field(i) + em_taylor(i)%term(j)%coef * &
                            expn(em_taylor(i)%term(j)%expn(1), 1) * &
                            expn(em_taylor(i)%term(j)%expn(2), 2)
  enddo
enddo

if (present(dfield)) then
  dfield = 0
  do i = 1, 3
    do j = 1, size(em_taylor(i)%term)
      iex = em_taylor(i)%term(j)%expn(1)
      iey = em_taylor(i)%term(j)%expn(2)
      if (iex > 0) dfield(i,1) = dfield(i,1) + iex * em_taylor(i)%term(j)%coef * expn(iex-1, 1) * expn(iey, 2)
      if (iey > 0) dfield(i,2) = dfield(i,2) + iey * em_taylor(i)%term(j)%coef * expn(iex, 1) * expn(iey-1, 2)
    enddo
  enddo

  dfield(1,3) = dfield(3,1)   ! Curl is zero
  dfield(2,3) = dfield(3,2)   ! Curl is zero
  dfield(3,3) = -(dfield(1,1) + dfield(2,2))  ! Divergence is zero
endif

end subroutine evaluate_em_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine add_em_taylor_term (em_taylor, coef, expn, replace)
!
! Routine to add a Em_taylor term to a Em_taylor series.
!
! If em_taylor does not have a term with the same exponents then a new
! term is added to em_taylor so the total number of terms is increased by one.
!
! If em_taylor already has a term with the same exponents then
! the "replace" argument determines what happens:
!   If replace = False (default) then
!      coef is added to the coefficient of the old term.
!   If replace = True then
!      coef replaces the coefficient of the old term.
! In both these cases, the number of terms in em_taylor remains the same.
!
! Input:
!   em_taylor   -- Em_taylor_struct: Em_taylor series.
!   coef        -- Real(rp): Coefficient.
!   expn(2)     -- Integer: Array of exponent indices.
!   replace     -- Logical, optional: Replace existing term? Default is False.
!
! Output:
!   em_taylor -- Em_taylor_struct: New series with term added
!-
!-

subroutine add_em_taylor_term (em_taylor, coef, expn, replace)

implicit none

type (em_taylor_struct) em_taylor

real(rp) coef
integer expn(:), i, n

logical, optional :: replace

! Search for an existing term of the same type

n = size(em_taylor%term)

do i = 1, n
  if (all(em_taylor%term(i)%expn == expn)) then
    if (logic_option(.false., replace)) then
      em_taylor%term(i)%coef = coef
    else
      em_taylor%term(i)%coef = coef + em_taylor%term(i)%coef
    endif
    if (em_taylor%term(i)%coef == 0) then  ! Kill this term
      em_taylor%term(i:n-1) = em_taylor%term(i+1:n)
      call init_em_taylor_series (em_taylor, n-1, .true.)
    endif
    return
  endif
enddo

! new term

call init_em_taylor_series (em_taylor, n+1, .true.)
em_taylor%term(n+1)%coef = coef
em_taylor%term(n+1)%expn = expn

end subroutine add_em_taylor_term

end module
