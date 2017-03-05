!+
! Function c_multi (n, m, no_n_fact, c_full) result (c_out)
!
! Subroutine to compute multipole factors:
!          c_multi(n, m) =  +/- ("n choose m")/n!
! This is used in calculating multipoles.
!
! Input:
!   n,m       -- Integer: For n choose m
!   no_n_fact -- Logical, optional: If present and true then
!                 c_out = +/- "n choose m".
!   c_full(:,:) --    real(rp), optional:  If present, will be populated with
!                     return entire c(n_pole_maxx,n_pole_maxx) matrix
!
! Output:
!   c_out  -- Real(rp): Multipole factor.
!-

function c_multi (n, m, no_n_fact, c_full) result (c_out)

use bmad_struct

implicit none

integer, intent(in) :: n, m
integer in, im

real(rp) c_out
real(rp), save :: n_factorial(0:n_pole_maxx)
real(rp), save :: c(0:n_pole_maxx, 0:n_pole_maxx)
real(rp), optional :: c_full(0:n_pole_maxx, 0:n_pole_maxx)

logical, save :: init_needed = .true.
logical, optional :: no_n_fact

! The magnitude of c(n, m) is number of combinations normalized by n!

if (init_needed) then

  c(0, 0) = 1

  do in = 1, n_pole_maxx
    c(in, 0) = 1
    c(in, in) = 1
    do im = 1, in-1
      c(in, im) = c(in-1, im-1) + c(in-1, im)
    enddo
  enddo

  n_factorial(0) = 1

  do in = 0, n_pole_maxx
    if (in > 0) n_factorial(in) = in * n_factorial(in-1)
    do im = 0, in
      c(in, im) = c(in, im) / n_factorial(in)
      if (mod(im, 4) == 0) c(in, im) = -c(in, im)
      if (mod(im, 4) == 3) c(in, im) = -c(in, im)
    enddo
  enddo

  init_needed = .false.

endif

!

if (logic_option (.false., no_n_fact)) then
  c_out = c(n, m) * n_factorial(n)
else
  c_out = c(n, m)
endif

if (present(c_full)) c_full = c

end function c_multi
