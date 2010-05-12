!+
! Module random_mod
!
! Module for random number generation.
!-

module random_mod

use precision_def
use physical_constants
use output_mod
use sim_utils_interface

integer, parameter :: kr4b = selected_int_kind(9)

type random_state_struct
  integer(kr4b) :: ix = -1, iy = -1
  logical :: number_stored = .false.
  real(rp) :: h_saved = 0
end type

type (random_state_struct), private, save :: r_state

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
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

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
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
real(sp), private, save :: am
integer(kr4b), parameter :: im = 2147483647

integer, parameter :: pseudo_random$ = 1, quasi_random$ = 2
integer, private, save :: engine = pseudo_random$

integer, parameter :: quick_gaussian$ = 3, exact_gaussian$ = 4
integer, private, save :: gauss_converter = exact_gaussian$
real(rp), save :: gauss_sigma_cut = -1
integer, private :: n_pts_per_sigma = 20

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_scaler (harvest, index_quasi)
!
! Subroutine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!-

subroutine ran_gauss_scaler (harvest, index_quasi)

use nr, only: erf_s, erf

implicit none

real(rp), intent(out) :: harvest
real(rp) a(2), v1, v2, r, sigma_cut, fac
real(rp), allocatable, save :: g(:)

integer, optional :: index_quasi
integer i, ss, ix, n_g_new
integer, save :: n_g = -1

! quasi-random must use the quick_gaussian since the exact_gaussian can
! use several uniform random numbers to generate a single Gaussian random number.
! This invalidates the algorithm used to generate a quasi-random Gaussian vector.

! g is the normalized error function and maps from the 
! interval [0, 0.5] to [0, infinity].

if (engine == quasi_random$ .or. gauss_converter == quick_gaussian$) then

  if (gauss_sigma_cut > 0) then
    sigma_cut = gauss_sigma_cut
  else 
    sigma_cut = 10
  endif
  n_g_new = max(nint(sigma_cut * n_pts_per_sigma), 10)

  if (n_g_new /= n_g) then
    n_g = n_g_new
    if (allocated(g)) deallocate(g)
    allocate (g(0:n_g))
    fac = 2 * erf_s (sigma_cut/sqrt_2)
    do i = 0, n_g-1
      g(i) = erf_s (sigma_cut * i / (sqrt_2 * n_g)) / fac
    enddo
    g(n_g) = 0.50000000001_rp
  endif

  call ran_uniform_scaler (r, index_quasi)
  if (r > 0.5) then
    r = r - 0.5
    ss = 1
  else
    r = 0.5 - r
    ss = -1
  endif
  call bracket_index(g, 0, n_g, r, ix)    
  harvest = ss * (ix + (g(ix+1) - r) / (g(ix+1) - g(ix))) * sigma_cut / n_g
  return

endif

! Loop until we get an acceptable number

do 

  ! If we have a stored value then just use it

  if (r_state%number_stored) then
    r_state%number_stored = .false.
    harvest = r_state%h_saved
    if (gauss_sigma_cut < 0 .or. abs(harvest) < gauss_sigma_cut) return
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
  r_state%h_saved = v2 * r
  r_state%number_stored = .true.

  harvest = v1 * r
  if (gauss_sigma_cut < 0 .or. abs(harvest) < gauss_sigma_cut) return

enddo

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
  call ran_gauss_scaler (harvest(i), i)
enddo

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
    r_state%number_stored = .false.
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
! quick_gaussian$ is a quick a dirty approximation with a cutoff so that no 
! numbers will be generated beyound what is set for sigma_cut. 
!
! A negative sigma_cut means that the exact_gaussian$ will not be limited
! and the quick_gaussian$ will use a default of 10.0
!
! Note: Because of technical issues, when using the quasi_random$ number generator
! (see the ran_engine routine), the quick_gaussian$ method will automatically be 
! used independent of what was set with this routine.
!
! Modules needed:
!   use random_mod
! 
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!             'exact'
!             'quick'  (Old deprecated: 'limited')
!   set_sigma_cut -- Real(rp), optional: Sigma cutoff. Initially: sigma_cut = -1.
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
  case (quick_gaussian$)
    get = 'quick'
  case (exact_gaussian$)
    get = 'exact'
  end select
endif

! get sigma_cut

if (present(get_sigma_cut)) get_sigma_cut = gauss_sigma_cut

! set converter

if (present(set)) then
  select case (set)
  case ('quick', 'limited')
    gauss_converter = quick_gaussian$
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
! Subroutine ran_seed_put (seed, state)
!
! Subroutine to seed the random number generator. 
!
! If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
! the system clock will be used to generate the seed.
!
! The seed is only used with the pseudo_random$ engine.
!
! If state is given, that will be used instead of seed
! 
! Use the subroutine ran_seed_get(seed) to get the seed used.
!
! Modules needed:
!   use random_mod
!
! Intput:
!   seed  -- Integer, optional: Seed number. If seed = 0 then a 
!              seed will be choosen based upon the system clock.
!   state -- random_state_struct, optional: Internal state of the random engine.
!-

subroutine ran_seed_put (seed, state)

use nr, only: sobseq

implicit none

type (random_state_struct), optional :: state
integer, optional :: seed
integer v(10)
real(rp) dum(2)

! Quasi-random number generator init

call sobseq (dum, 0)

am = nearest(1.0,-1.0) / im

! Init Random number generator with state

if (present(state)) then
  r_state = state
  return
endif

! init

if (seed == 0) then
  call date_and_time (values = v)
  my_seed = v(2) + 11*v(3) + 111*v(5) + 1111*v(6) + 11111*v(7) + 111111*v(8)
else
  my_seed = seed
endif

r_state%iy = ior(ieor(888889999, abs(my_seed)), 1)
r_state%ix = ieor(777755555, abs(my_seed))

r_state%number_stored = .false.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_get (seed, state)
! 
! Subroutine to return the seed used for the random number generator.
!
! Note: The internal state can be used to put the pseudo-random
! number generator into a known state. See ran_seed_put
!
! Modules needed:
!   use random_mod
!
! Output:
!   seed  -- Integer, optional: Random number seed used.
!   state -- random_state_struct, optional: Internal state of the random engine.
!-

subroutine ran_seed_get (seed, state)

implicit none

type (random_state_struct), optional :: state
integer, optional :: seed

!

if (present(seed)) seed = my_seed
if (present(state)) state = r_state

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_scaler (harvest, index_quasi)
!
! Subroutine to return a random number uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!-

subroutine ran_uniform_scaler (harvest, index_quasi)

use nr, only: sobseq

implicit none

real(rp), intent(out) :: harvest
real(rp), save :: r(6)

integer(kr4b) k, ix_q
integer, optional :: index_quasi

integer(kr4b), parameter :: ia = 16807
integer(kr4b), parameter :: iq = 127773, ir = 2836

character :: r_name = 'ran_uniform_scaler'

! If r_state%iy < 0 then the random number generator has never bee initialized.

if (r_state%iy < 0) call ran_seed_put(my_seed)

! quasi-random

if (engine == quasi_random$) then
  ix_q = integer_option(1, index_quasi)
  if (ix_q == 1) call sobseq (r)
  if (ix_q > 6) then
    call out_io (s_error$, r_name, 'NUMBER OF DIMENSIONS WANTED IS TOO LARGE!')
    call err_exit
  endif
  harvest = r(ix_q)
  return
endif

! Pseudo-random
! Marsaglia shift sequence with period 2^32 - 1.

r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 13)) 
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, -17))
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 5))
k = r_state%iy/iq         ! Park-Miller sequence by Schrage's method,
r_state%iy = ia*(r_state%iy - k*iq) - ir * k       ! period 2^31-2.

if (r_state%iy < 0) r_state%iy = r_state%iy + im

! Combine the two generators with masking to ensure nonzero value.

harvest = am * ior(iand(im, ieor(r_state%ix, r_state%iy)), 1) 

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
  call ran_uniform_scaler (harvest(i), i)
enddo

end subroutine

end module

