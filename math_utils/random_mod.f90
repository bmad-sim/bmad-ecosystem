!+
! Module random_mod
!
! Module for random number generation.
!-

#include "CESR_platform.inc"

module random_mod

use precision_def
use physical_constants
use output_mod
use dcslib_interface

!+
! Subroutine ran_gauss (harvest)
!
! Subroutine to return a gaussian distributed random number with unit sigma.
! This routine uses the same algorithm as gasdev from Numerical Recipes.
!
! Note: ran_gauss is an overloaded name for:
!     ran_gauss_scaler   ! harvest is a scaler
!     ran_gauss_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
! Note: Use ran_gauss_converter to set which conversion routine to use.
!
! Modules needed:
!   use random_mod
!
! Output:
!   harvest    -- Real(rp): Random number. 
! Or
!   harvest(:) -- Real(rp): Random number array. 
!                  For quasi_random$ numbers, the array size must be less than 6.
!-

interface ran_gauss
  module procedure ran_gauss_scaler
  module procedure ran_gauss_vector
end interface

!+
! Subroutine ran_uniform (harvest)
!
! Subroutine to return a random number uniformly distributed in the 
! interval [0, 1]. This routine uses the same algorithm as ran or sobseq
! from Numberical Recipes in Fortran90.
! See ran_engine.
!
! Note: ran_uniform is an overloaded name for:
!     ran_uniform_scaler   ! harvest is a scaler
!     ran_uniform_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
!
! Modules needed:
!   use random_mod
!
! Output:
!   harvest    -- Real(rp): Random number. 
! Or
!   harvest(:) -- Real(rp): Random number array. 
!                  For quasi_random$ numbers the array size must be less than 6.
!-

interface ran_uniform
  module procedure ran_uniform_scaler
  module procedure ran_uniform_vector
end interface

! common variables for random number generator.

integer, private, save :: my_seed = 0
integer, private, parameter :: k4b=selected_int_kind(9)
integer(k4b), private, save :: ix = -1, iy = -1, k
real(sp), private, save :: am
integer(k4b), private, parameter :: ia = 16807, im = 2147483647
integer(k4b), private, parameter :: iq = 127773, ir = 2836
logical, private, save :: number_stored = .false.

! index_quasi is used so that the scaler routines know which component of the
! quasi-random vector to return.

integer, private, save :: index_quasi = 1

integer, parameter :: pseudo_random$ = 1, quasi_random$ = 2
integer, private, save :: engine = pseudo_random$

integer, parameter :: limited_gaussian$ = 3, exact_gaussian$ = 4
integer, private, save :: gauss_converter = exact_gaussian$
real(rp), save :: gauss_sigma_cut = 4.0
integer :: n_pts_per_sigma = 20

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_scaler (harvest)
!
! Subroutine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!-

subroutine ran_gauss_scaler (harvest)

use nr, only: erf_s

implicit none

real(rp), intent(out) :: harvest
real(rp) a(2), v1, v2, r
real(rp), save :: h_saved
real(rp), allocatable, save :: g(:)

integer i, ss, ix
integer, save :: n_g

logical, save :: init_needed = .true.

! quasi-random must use the limited_gaussian.
! g is the normalized error function and maps from the 
! interval [0, 0.5] to [0, infinity].

if (engine == quasi_random$ .or. gauss_converter == limited_gaussian$) then
  if (init_needed) then
    n_g = gauss_sigma_cut*n_pts_per_sigma
    allocate (g(0:n_g))
    do i = 0, n_g-1
      g(i) = 0.5 * erf_s (gauss_sigma_cut * i / (sqrt_2 * n_g))
    enddo
    g(n_g) = 0.50000001
    init_needed = .false.
  endif
  call ran_uniform(r)
  if (r > 0.5) then
    r = r - 0.5
    ss = 1
  else
    r = 0.5 - r
    ss = -1
  endif
  call bracket_index(g, 0, n_g, r, ix)    
  harvest = ss * (ix + (g(ix+1) - r) / (g(ix+1) - g(ix))) / n_pts_per_sigma
  return
endif

! If we have a stored value then just use it

if (number_stored) then
  harvest = h_saved
  number_stored = .false.
  return
endif

! else we generate a number

do
  call ran_uniform(a)
  v1 = 2*a(1) - 1
  v2 = 2*a(2) - 1
  r = v1**2 + v2**2
  if (r > 0 .and. r < 1) exit   ! In unit circle
enddo

r = sqrt(-2*log(r)/r)
harvest = v1 * r
h_saved = v2 * r
number_stored = .true.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_vector (harvest)
!
! Subroutine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!-

subroutine ran_gauss_vector (harvest)

implicit none

real(rp), intent(out) :: harvest(:)
integer i

!

do i = 1, size(harvest)
  index_quasi = i
  call ran_gauss_scaler (harvest(i))
enddo

index_quasi = 1

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_engine (set, get)
!
! Subroutine to set what random number generator algorithm is used.
! If this routine is never called then pseudo_random$ is used.
! With sobseq quasi-random numbers the maximum dimension is 6.
!
! Modules needed:
!   use random_mod
! 
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!                'pseudo' -> Uses ran from Numerical Recipies (F90).
!                'quasi'  -> Uses sobseq from Numerical Recipes.
!   get -- Character, optional: Get the current (before any set) random number engine. 
!-

subroutine ran_engine (set, get)

implicit none

character(*), optional :: set, get
character(16) :: r_name = 'ran_engine'

! get

if (present (get)) then
  select case (engine)
  case (pseudo_random$)
    get = 'pseudo'
  case (quasi_random$)
    get = 'quasi'
  end select
endif

! set

if (present(set)) then
  select case (set)
  case ('pseudo')
    engine = pseudo_random$
  case ('quasi')
    engine = quasi_random$
    number_stored = .false.
    index_quasi = 1
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER ENGINE NAME: ' // set)
  end select
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut)
!
! Subroutine to set what conversion routine is used for converting
! uniformly distributed random numbers to Gaussian distributed random numbers.
!
! If this routine is not called then exact_gaussian$ is used.
!
! exact_gaussian$ is a straight forward converter as explained in Numerical recipes.
!
! limited_gaussian$ is a quick a dirty approximation with a cutoff so that no 
! numbers will be generated beyound what is set for sigma_cut. 
!
! Note: Because of technical issues, when using the quasi_random$ number generator
! (see the ran_engine routine), the limited_gaussian$ method will automatically be 
! used independent of what was set with this routine.
!
! Modules needed:
!   use random_mod
! 
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!             'exact'
!             'limited'
!   set_sigma_cut -- Real(rp), optional: Sigma cutoff used with the 
!                 limited_gaussian$ converter. Initially: sigma_cut = 4.0.
!
! Output:
!   get -- Character(*), optional: Get the current (before any set) gaussian converter.
!   get_sigma_cut -- Real(rp), optional: Get the current (before andy set) sigma cutoff.
!-

subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut)

implicit none

real(rp), optional :: set_sigma_cut, get_sigma_cut
character(*), optional :: set, get
character(16) :: r_name = 'ran_gauss_converter'

! get converter

if (present (get)) then
  select case (gauss_converter)
  case (limited_gaussian$)
    get = 'limited'
  case (exact_gaussian$)
    get = 'exact'
  end select
endif

! get sigma_cut

if (present(get_sigma_cut)) get_sigma_cut = gauss_sigma_cut

! set converter

if (present(set)) then
  select case (set)
  case ('limited')
    gauss_converter = limited_gaussian$
  case ('exact')
    gauss_converter = exact_gaussian$
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER GAUSS_CONVERTER NAME: ' // set)
  end select
endif

! set sigma_cut

if (present(set_sigma_cut)) gauss_sigma_cut = set_sigma_cut

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_put (seed)
!
! Subroutine to seed the random number generator. 
!
! If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
! the system clock will be used to generate the seed.
!
! The seed is only used with the pseudo_random$ engine.
!
! Use the subroutine ran_seed_get(seed) to get the seed used.
!
! Modules needed:
!   use random_mod
!
! Intput:
!   seed -- Integer: Seed number. If seed = 0 then a 
!             seed will be choosen based upon the system clock.
!-

subroutine ran_seed_put (seed)

use nr, only: sobseq

implicit none

integer seed
integer v(10)
real(rp) dum(2)

! Random number generator init
! Only init if seed >= 0

if (seed > -1) then

  if (seed == 0) then
    call date_and_time (values = v)
    my_seed = v(2) + 11*v(3) + 111*v(5) + 1111*v(6) + 11111*v(7) + 111111*v(8)
  else
    my_seed = seed
  endif

  am = nearest(1.0,-1.0)/IM
  iy = ior(ieor(888889999, abs(my_seed)), 1)
  ix = ieor(777755555, abs(my_seed))
  number_stored = .false.

endif

! Quasi-random number generator init

call sobseq (dum, 0)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_get (seed)
! 
! Subroutine to return the seed used for the random number generator.
!
! Modules needed:
!   use random_mod
!
! Output:
!   seed -- Integer: Random number seed used.
!-

subroutine ran_seed_get (seed)

implicit none
integer seed

seed = my_seed

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_scaler (harvest)
!
! Subroutine to return a random number uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!-

subroutine ran_uniform_scaler (harvest)

use nr, only: sobseq

implicit none

real(rp), intent(out) :: harvest
real(rp), save :: r(6)
character :: r_name = 'ran_uniform_scaler'

!

if (iy < 0) call ran_seed_put(my_seed)

! quasi-random

if (engine == quasi_random$) then
  if (index_quasi == 1) call sobseq (r)
  if (index_quasi > 6) then
    call out_io (s_error$, r_name, 'NUMBER OF DIMENSIONS WANTED IS TOO LARGE!')
    call err_exit
  endif
  harvest = r(index_quasi)
  return
endif

! Pseudo-random

ix = ieor(ix, ishft(ix, 13)) !Marsaglia shift sequence with period 2^32 - 1.
ix = ieor(ix, ishft(ix, -17))
ix = ieor(ix, ishft(ix, 5))
k = iy/iq                  ! Park-Miller sequence by Schrage's method,
iy = ia*(iy-k*iq)-ir*k     ! period 2^31-2.

if (iy < 0) iy = iy+im

! Combine the two generators with masking to ensure nonzero value.

harvest = am * ior(iand(im, ieor(ix,iy)), 1) 

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_vector (harvest)
!
! Subroutine to return a vector of random numbers uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!-

subroutine ran_uniform_vector (harvest)

implicit none

real(rp), intent(out) :: harvest(:)
integer i

!

do i = 1, size(harvest)
  index_quasi = i
  call ran_uniform_scaler (harvest(i))
enddo

index_quasi = 1  ! reset

end subroutine

end module

