#include "CESR_platform.inc"

module bmad_taylor_mod

  use precision_def

! Note: the taylor_struct uses the bmad standard (x, p_x, y, p_y, z, p_z) 
! the universal_taylor in Etienne's PTC uses (x, p_x, y, p_y, p_z, -c*t)
! %ref is the reference point about which the taylor expansion was made

  type taylor_term_struct
    real(rp) :: coef
    integer :: exp(6)  
  end type

  type taylor_struct
    real (rp) ref             
    type (taylor_term_struct), pointer :: term(:) => null()
  end type

!

  interface assignment (=)
    module procedure taylor_equal_taylor
  end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine taylor_equal_taylor (taylor1, taylor2)
!
! Subroutine to transfer the values from one taylor map to another:
!     Taylor1 <= Taylor2
!
! Modules needed:
!   use bmad
!
! Input:
!   taylor2(:) -- Taylor_struct: Taylor map.
!
! Output:
!   taylor1(:) -- Taylor_struct: Taylor map. 
!-

subroutine taylor_equal_taylor (taylor1, taylor2)

  implicit none

  type (taylor_struct), intent(inout) :: taylor1(:)
  type (taylor_struct), intent(in)    :: taylor2(:)

  integer i

!

  if (associated (taylor1(1)%term)) then
    do i = 1, size(taylor1)
      deallocate (taylor1(i)%term)
    enddo
  endif

  do i = 1, size(taylor1)
    allocate (taylor1(i)%term(size(taylor2(i)%term)))
    taylor1(i)%term = taylor2(i)%term
  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function taylor_coef (bmad_taylor, exp)
!
! Function to return the coefficient for a particular taylor term
! from a Taylor Series.
! For example: To get the term corresponding to 
!   y_out = Coef * p_z^2 
! [This is somtimes refered to as the T_366 term]
! The call would be
!   type (taylor_struct) bmad_taylor(6)      ! Taylor Map
!   ...
!   coef = taylor_coef (bmad_taylor(3), (/ 0, 0, 0, 0, 0, 2 /) )
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor -- Taylor_struct: Taylor series.
!   exp(6)      -- Integer: Array of exponent indices.
!
! Output:
!   taylor_coef -- Real(rp): Coefficient.
!-

function taylor_coef (bmad_taylor, exp) result (coef)

  implicit none

  type (taylor_struct), intent(in) :: bmad_taylor

  real(rp) coef

  integer, intent(in) :: exp(:)
  integer i

!

  coef = 0

  do i = 1, size(bmad_taylor%term)
    if (all(bmad_taylor%term(i)%exp == exp)) then
      coef = bmad_taylor%term(i)%coef
      return
    endif
  enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine type_taylors (bmad_taylor)
!
! Subroutine to print in a nice format a BMAD Taylor Map at the terminal.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: 6 taylor series: (x, P_x, y, P_y, z, P_z) 
!-

subroutine type_taylors (bmad_taylor)

  implicit none

  type (taylor_struct), target :: bmad_taylor(:)
  integer i, n_lines
  character*80, pointer :: lines(:)

!

  nullify (lines)

  call type2_taylors (bmad_taylor, lines, n_lines)

  do i = 1, n_lines
    print *, trim(lines(i))
  enddo

  deallocate(lines)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine type2_taylors (bmad_taylor, lines, n_lines)
!
! Subroutine to write a BMAD taylor map in a nice format to a character array.
!
! Moudles needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Array of taylors.
!
! Output:
!   lines(:)     -- Character*80, allocatable: Character array to hold the 
!                     output. The array size of lines(:) will be set by
!                     this subroutine.
!   n_lines      -- Number of lines in lines(:).
!-

subroutine type2_taylors (bmad_taylor, lines, n_lines)

  implicit none

  type (taylor_struct), intent(in), target :: bmad_taylor(6)
  type (taylor_term_struct), pointer :: tt
  type (taylor_struct) tlr

  integer, intent(out) :: n_lines
  integer i, j, k, nl, ix

  character*80, pointer :: lines(:)
  character*40, fmt1, fmt2, fmt

! If not allocated then not much to do

  if (.not. associated(bmad_taylor(1)%term)) then
    n_lines = 2
    allocate (lines(n_lines))
    lines(1) = '---------------------------------------------------'
    lines(2) = 'A Taylor Map Does Not Exist.' 
    return
  endif

! Normal case

  deallocate (lines, stat = ix)
  n_lines = 8 + sum( (/ (size(bmad_taylor(i)%term), i = 1, 6) /) )
  allocate(lines(n_lines))
  

  write (lines(1), *) 'Taylor Terms:'
  write (lines(2), *) &
        'Out     Coef              Exponents           Order        Reference'
  nl = 2


  fmt1 = '(i4, a, f20.12, 6i3, i9, f18.9)'
  fmt2 = '(i4, a, 1p, e20.11, 0p, 6i3, i9, f18.9)'

  do i = 1, 6
    lines(nl+1) = '---------------------------------------------------'
    nl = nl + 1

    nullify (tlr%term)
    call sort_taylor_terms (bmad_taylor(i), tlr)

    do j = 1, size(bmad_taylor(i)%term)

      tt => tlr%term(j)

      if (abs(tt%coef) < 1e5) then
        fmt = fmt1
      else
        fmt = fmt2
      endif

      if (j == 1) then
        write (lines(nl+j), fmt) i, ':', tt%coef, &
                    (tt%exp(k), k = 1, 6), sum(tt%exp), bmad_taylor(i)%ref
      else
        write (lines(nl+j), fmt) i, ':', tt%coef, &
                    (tt%exp(k), k = 1, 6), sum(tt%exp)
      endif
    enddo

    nl = nl + size(bmad_taylor(i)%term)
    deallocate (tlr%term)

  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_taylor (bmad_taylor)
!
! Subroutine to initialize (nullify) a BMAD Taylor map.
! Note: You must be sure that bmad_taylor is not allocated before using this
! routine otherwise you will have a memory leak!
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: New structure.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Initalized structure.
!-

subroutine init_taylor (bmad_taylor)

  implicit none

  type (taylor_struct) bmad_taylor(:)

  integer i

!

  do i = 1, size(bmad_taylor)
    nullify (bmad_taylor(i)%term)
  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine kill_taylor (bmad_taylor)
!
! Subroutine to deallocate a BMAD taylor map.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Taylor to be deallocated. It is OK
!                       if bmad_taylor has already been deallocated.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: deallocated Taylor structure.
!-

subroutine kill_taylor (bmad_taylor)

  implicit none

  type (taylor_struct) bmad_taylor(:)

  integer i, ix

!

  do i = 1, size(bmad_taylor)
    deallocate (bmad_taylor(i)%term, stat = ix)
  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! subroutine sort_taylor_terms (taylor_in, taylor_sorted)
!
! Subroutine to sort the taylor terms from "lowest" to "highest" of
! a taylor series.
! This subroutine is needed since what comes out of PTC is not sorted.
!
! The number associated with a taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! Where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
!
! Note: taylor_sorted needs to have been initialized.
! Note: taylor_sorted cannot be taylor_in. That is it is not legal to write:
!           call sort_taylor_terms (this_taylor, this_taylor)
!
! Modules needed:
!   use bmad
!
! Input:
!   taylor_in     -- Taylor_struct: Unsorted taylor series.
!
! Output:
!   taylor_sorted -- Taylor_struct: Sorted taylor series.
!-

subroutine sort_taylor_terms (taylor_in, taylor_sorted)

  use nr

  implicit none

  type (taylor_struct), intent(in)  :: taylor_in
  type (taylor_struct) :: taylor_sorted
  type (taylor_term_struct), allocatable :: tt(:)
  
  integer, allocatable :: ord_(:), ix_(:)

  integer i, j, n, expn(6)

! init

  n = size(taylor_in%term)
  if (associated(taylor_sorted%term)) deallocate(taylor_sorted%term)
  allocate(taylor_sorted%term(n), ix_(n), ord_(n), tt(n))

!

  tt = taylor_in%term

  do i = 1, n
    expn = tt(i)%exp
    ord_(i) = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
                expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
  enddo

  call indexx (ord_, ix_)

  do i = 1, n
    taylor_sorted%term(i)= tt(ix_(i))
  enddo

  deallocate(ord_, ix_, tt)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine taylor_to_mat6 (a_taylor, c0, mat6, c1)
!
! Subroutine to calculate the linear (Jacobian) matrix about 
! some trajectory from a Taylor map.
!
! Modules needed:
!   use bmad
!
! Input:
!   a_taylor(6) -- Taylor_struct: Taylor map.
!   c0          -- Coord_struct: Coordinates at the input. 
!
! Output:
!   mat6(6,6) -- Real(rp): Jacobian.
!   c1(6)     -- Real(rp): Oth order transport vector.
!-

subroutine taylor_to_mat6 (a_taylor, c0, mat6, c1)

  implicit none

  type (taylor_struct), target, intent(in) :: a_taylor(6)
  real(rp), intent(in) :: c0(:)
  real(rp), intent(out) :: c1(:)
  type (taylor_term_struct), pointer :: term

  real(rp), intent(out) :: mat6(6,6)
  real(rp) prod

  integer i, j, k, l

! mat6 calc

  mat6 = 0
  
  do i = 1, 6
    do j = 1, 6
      do k = 1, size(a_taylor(i)%term)

        term => a_taylor(i)%term(k)

        if (term%exp(j) == 0) cycle
        if (term%exp(j) > 1 .and. c0(j) == 0) cycle

        if (term%exp(j) > 1)then
          prod = term%coef * term%exp(j) * c0(j) ** (term%exp(j)-1)
        else
          prod = term%coef
        endif

        do l = 1, 6
          if (term%exp(l) == 0) cycle
          if (l == j) cycle
          prod = prod * c0(l) ** term%exp(l)
        enddo

        mat6(i,j) = mat6(i,j) + prod

      enddo
    enddo
  enddo

! c1 calc

  call track_taylor (c0, a_taylor, c1)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine mat6_to_taylor (mat6, vec0, bmad_taylor)
!
! Subroutine to form a first order Taylor map from the 6x6 transfer
! matrix and the 0th order transfer vector.
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6) -- 6x6 transfer matrix.
!   vec0(6)   -- 0th order transfer vector.
!
! Output:
!   bmad_taylor(6) -- Taylor_struct: first order taylor map.
!-

subroutine mat6_to_taylor (mat6, vec0, bmad_taylor)

  implicit none

  type (taylor_struct) bmad_taylor(6)

  real(rp), intent(in) :: mat6(6,6)
  real(rp), intent(in) :: vec0(6)

  integer i, j, n
!

  call kill_taylor(bmad_taylor)

  do i = 1, 6
    n = count(mat6(i,1:6) /= 0)
    if (vec0(i) /= 0) n = n + 1
    allocate (bmad_taylor(i)%term(n))

    n = 0

    if (vec0(i) /= 0) then
      n = n + 1
      bmad_taylor(i)%term(1)%coef = vec0(i)
      bmad_taylor(i)%term(1)%exp = 0
    endif

    do j = 1, 6
      if (mat6(i,j) /= 0) then
        n = n + 1
        bmad_taylor(i)%term(n)%coef = mat6(i,j)
        bmad_taylor(i)%term(n)%exp = 0
        bmad_taylor(i)%term(n)%exp(j) = 1
      endif
    enddo

  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track_taylor (start, bmad_taylor, end)
!
! Subroutine to track using a Taylor map.
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Taylor map.
!   start          -- Real(rp): Starting coords.
!
! Output:
!   end            -- Real(rp): Ending coords.
!-

subroutine track_taylor (start, bmad_taylor, end)

  implicit none
  
  type (taylor_struct), intent(in) :: bmad_taylor(6)
  real(rp), intent(in) :: start(6)
  real(rp), intent(out) :: end(6)
  
  real(rp) delta
  
  integer i, j, k, ie
  
!

  end = 0
  do i = 1, 6
    j_loop: do j = 1, size(bmad_taylor(i)%term)
      delta =  bmad_taylor(i)%term(j)%coef
      do k = 1, 6
        ie = bmad_taylor(i)%term(j)%exp(k) 
        if (ie == 0) cycle
        if (start(k) == 0) cycle j_loop  ! delta = 0 in this case 
        delta = delta * start(k) ** ie
      enddo
      end(i) = end(i) + delta
    enddo j_loop
  enddo

end subroutine


end module
