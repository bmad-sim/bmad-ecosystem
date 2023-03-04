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

!

integer, private, parameter :: kr4b = selected_int_kind(9)
integer(kr4b), private, parameter :: im_nr_ran = 2147483647
integer(i4_b), private, parameter :: sobseq_maxbit = 30, sobseq_maxdim = 6

! common variables for random number generator.

character(8), parameter :: ran_engine_name(2) = [character(8):: 'pseudo', 'quasi']
character(8), parameter :: ran_gauss_converter_name(4) = [character(8):: '', '', 'quick', 'exact']

integer, parameter :: pseudo_random$ = 1, quasi_random$ = 2
integer, parameter :: quick_gaussian$ = 3, exact_gaussian$ = 4

type random_state_struct
  integer(kr4b) :: ix = -1, iy = -1
  logical :: number_stored = .false.
  real(rp) :: h_saved = 0
  integer :: engine = pseudo_random$
  ! Params
  integer :: seed = 0
  real(sp) :: am = 0
  integer :: gauss_converter = exact_gaussian$
  real(rp) :: gauss_sigma_cut = -1  ! Only used if positive.
  integer(i4_b) :: in_sobseq = 0
  integer(i4_b) :: ix_sobseq(sobseq_maxdim) = 0
  real(rp) :: x_sobseq(sobseq_maxdim) = 0
end type

type (random_state_struct), private, target, save :: ran_state_save
type (random_state_struct), private, target, save, allocatable :: thread_ran_state(:)
logical, private, save :: thread_state_allocated = .false.

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine ran_gauss (harvest, ran_state, sigma_cut)
!
! Routine to return a gaussian distributed random number with unit sigma.
! This routine uses the same algorithm as gasdev from Numerical Recipes.
!
! Note: ran_gauss is an overloaded name for:
!     ran_gauss_scalar   ! harvest is a scalar
!     ran_gauss_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
! Note: Use ran_gauss_converter to set which conversion routine to use.
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                   See the ran_seed_put documentation for more details.
!   sigma_cut -- real(rp), optional: If present and positive will override setting 
!                   of ran_state%gauss_sigma_cut.
!
! Output:
!   harvest    -- real(rp): Random number. 
! Or
!   harvest(:) -- real(rp): Random number array. 
!                  For quasi_random$ numbers, the array size must be less than 6.
!-

interface ran_gauss
  module procedure ran_gauss_scalar
  module procedure ran_gauss_vector
end interface

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine ran_uniform (harvest, ran_state)
!
! Routine to return a random number uniformly distributed in the 
! interval [0, 1]. This routine uses the same algorithm as ran or sobseq
! from Numberical Recipes in Fortran90.
! See ran_engine.
!
! Note: ran_uniform is an overloaded name for:
!     ran_uniform_scalar   ! harvest is a scalar
!     ran_uniform_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   harvest    -- Real(rp): Random number. 
! Or
!   harvest(:) -- Real(rp): Random number array. 
!                  For quasi_random$ numbers the array size must be less than 6.
!-

interface ran_uniform
  module procedure ran_uniform_scalar
  module procedure ran_uniform_vector
end interface

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_scalar (harvest, ran_state, sigma_cut, index_quasi, sigma_cut_ued)
!
! Routine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!
! Note: The index_quasi argument is used internally for the quasi-random number generator.
!-

subroutine ran_gauss_scalar (harvest, ran_state, sigma_cut, index_quasi, sigma_cut_used)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

real(rp), intent(out) :: harvest
real(rp), optional :: sigma_cut
real(rp), optional :: sigma_cut_used
real(rp) a(2), v1, v2, r, sig_cut, fac
real(rp), parameter :: sigma_max = 8

integer, parameter :: n_pts_per_sigma = 25
integer, parameter :: max_g = sigma_max * n_pts_per_sigma
integer, optional :: index_quasi
integer i, ss, ix

real(rp), save :: erf_array(0:max_g) = 0

! quasi-random must use the quick_gaussian since the exact_gaussian can
! use several uniform random numbers to generate a single Gaussian random number.
! This invalidates the algorithm used to generate a quasi-random Gaussian vector.

! ran_state%g is the normalized error function and maps from the 
! interval [0, 0.5] to [0, infinity].

r_state => pointer_to_ran_state(ran_state)

sig_cut = 1000
if (r_state%gauss_sigma_cut > 0) sig_cut = r_state%gauss_sigma_cut
if (present(sigma_cut)) then
  if (sigma_cut > 0) sig_cut = sigma_cut
endif

!

if (r_state%engine == quasi_random$ .or. r_state%gauss_converter == quick_gaussian$) then
  ! Init g

  sig_cut = min(sigma_max, sig_cut)

  if (erf_array(1) == 0) then
    fac = 2 * erf (sigma_max/sqrt_2)
    do i = 0, max_g-1
      erf_array(i) = erf (i / (n_pts_per_sigma * sqrt_2)) / fac
    enddo
    erf_array(max_g) = 0.50000000001_rp
  endif

  !

  call ran_uniform_scalar (r, ran_state, index_quasi)
  if (r > 0.5) then
    r = r - 0.5
    ss = 1
  else
    r = 0.5 - r
    ss = -1
  endif

  if (sig_cut < sigma_max) then
    fac = n_pts_per_sigma * sig_cut
    ix = int(fac)
    fac = fac - ix
    fac = 2.0_rp * (erf_array(ix) * (1 - fac) + erf_array(ix+1) * fac)
    r = fac * r  
  endif

  ix = bracket_index(r, erf_array, 0)
  harvest = (ix + (r - erf_array(ix)) / (erf_array(ix+1) - erf_array(ix))) * ss / n_pts_per_sigma
  if (present(sigma_cut_used)) sigma_cut_used = sig_cut
  return
endif

! Loop until we get an acceptable number

do 
  ! If we have a stored value then just use it

  if (r_state%number_stored) then
    r_state%number_stored = .false.
    harvest = r_state%h_saved
    if (sig_cut < 0 .or. abs(harvest) < sig_cut) exit
  endif

  ! else we generate a number

  do
    call ran_uniform(a, ran_state)
    v1 = 2*a(1) - 1
    v2 = 2*a(2) - 1
    r = v1**2 + v2**2
    if (r > 0 .and. r < 1) exit   ! In unit circle
  enddo

  r = sqrt(-2*log(r)/r)
  r_state%h_saved = v2 * r
  r_state%number_stored = .true.

  harvest = v1 * r
  if (abs(harvest) < sig_cut) exit
enddo

if (present(sigma_cut_used)) sigma_cut_used = sig_cut

end subroutine ran_gauss_scalar

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_vector (harvest, ran_state, sigma_cut)
!
! Routine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!-

subroutine ran_gauss_vector (harvest, ran_state, sigma_cut)

implicit none

type (random_state_struct), optional, target :: ran_state

real(rp), optional :: sigma_cut
real(rp), intent(out) :: harvest(:)
real(rp) sigma_cut_actual
integer i

!

do
  do i = 1, size(harvest)
    call ran_gauss_scalar (harvest(i), ran_state, sigma_cut, i, sigma_cut_actual)
  enddo

  if (sum(harvest*harvest) < sigma_cut_actual**2) return
enddo

end subroutine ran_gauss_vector

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_engine (set, get, ran_state)
!
! Routine to set what random number generator algorithm is used.
! If this routine is never called then pseudo_random$ is used.
! With sobseq quasi-random numbers the maximum dimension is 6.
!
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!                'pseudo' -> Uses ran from Numerical Recipies (F90).
!                'quasi'  -> Uses sobseq from Numerical Recipes.
!                ''       -> Do nothing.
!   get -- Character, optional: Get the current (before any set) random number engine. 
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!-

subroutine ran_engine (set, get, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

character(*), optional :: set, get
character(16) :: r_name = 'ran_engine'

! Get

r_state => pointer_to_ran_state(ran_state)

if (present (get)) then
  select case (r_state%engine)
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
    r_state%engine = pseudo_random$
  case ('quasi')
    r_state%engine = quasi_random$
    r_state%number_stored = .false.
  case ('')
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER ENGINE NAME: ' // set)
  end select
endif

end subroutine ran_engine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state)
!
! Routine to set what conversion routine is used for converting
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
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!             'exact'
!             'quick'  ! Old deprecated: 'limited'
!             ''       ! Do nothing
!   set_sigma_cut -- Real(rp), optional: Sigma cutoff. Initially: sigma_cut = -1.
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   get -- Character(*), optional: Get the current (before any set) gaussian converter.
!   get_sigma_cut -- Real(rp), optional: Get the current (before any set) sigma cutoff.
!-

subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

real(rp), optional :: set_sigma_cut, get_sigma_cut

character(*), optional :: set, get
character(16) :: r_name = 'ran_gauss_converter'

! Get converter

r_state => pointer_to_ran_state(ran_state)

if (present (get)) then
  select case (r_state%gauss_converter)
  case (quick_gaussian$)
    get = 'quick'
  case (exact_gaussian$)
    get = 'exact'
  end select
endif

! get sigma_cut

if (present(get_sigma_cut)) then
  get_sigma_cut = r_state%gauss_sigma_cut
endif

! set converter

if (present(set)) then
  select case (set)
  case ('quick', 'limited')
    r_state%gauss_converter = quick_gaussian$
  case ('exact')
    r_state%gauss_converter = exact_gaussian$
  case ('')
    ! Do nothing
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER GAUSS_CONVERTER NAME: ' // set)
  end select
endif

! set sigma_cut

if (present(set_sigma_cut)) then
  r_state%gauss_sigma_cut = set_sigma_cut
endif

end subroutine ran_gauss_converter

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_put (seed, mpi_offset)
!
! Routine to seed a random number generator. 
!
! If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
! the system clock will be used to generate the seed.
!
! Note: The seed is only used with the pseudo_random$ engine.
! Note: Use the subroutine ran_seed_get(seed) to get the seed used.
! Note: Use pointer_to_ran_state() to access the ran state directly.
!
! Input:
!   seed        -- integer, optional: Seed number. If seed = 0 then a 
!                   seed will be choosen based upon the system clock.
!   mpi_offset  -- integer, optional: Offset added to seed. Default is zero.
!                   Used with MPI processes ensure different threads use different random numbers.
!-

subroutine ran_seed_put (seed, mpi_offset)

!$ use omp_lib, only: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

implicit none

integer :: seed
integer nt, max_t, n
integer, optional :: mpi_offset

real(rp) dum(2)

! OpenMP put

!$  if (.not. thread_state_allocated) call allocate_thread_states()
!$  nt = OMP_GET_THREAD_NUM()
!$  max_t = OMP_GET_MAX_THREADS()
!$  if (nt == 0) then
!$    do n = 0, max_t-1
!$      call this_seed_put(seed, n)
!$    enddo
!$  endif
!$OMP BARRIER
!$  return

! Non-OpenMP put

call this_seed_put (seed, integer_option(0, mpi_offset))

!---------------------------------------------
contains

subroutine this_seed_put (seed, nt)

type (random_state_struct), pointer :: r_state

integer :: seed
integer v(10), nt

!

r_state => pointer_to_ran_state(ix_thread = nt)

r_state%in_sobseq = 0
r_state%ix_sobseq = 0

r_state%am = nearest(1.0,-1.0) / im_nr_ran

if (seed == 0) then
  call date_and_time (values = v)
  r_state%seed = v(1) + v(2) + 11*v(3) + 111*v(5) + 1111*v(6) + 11111*v(7) + 111111*v(8) + nt
else
  r_state%seed = seed + nt
endif

r_state%iy = ior(ieor(888889999, abs(r_state%seed)), 1)
r_state%ix = ieor(777755555, abs(r_state%seed))

r_state%number_stored = .false.

end subroutine this_seed_put

end subroutine ran_seed_put

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_get (seed)
! 
! Routine to return the seed used for the random number generator.
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   seed      -- Integer, optional: Random number seed used.
!-

subroutine ran_seed_get (seed)

implicit none

type (random_state_struct), pointer :: r_state

integer :: seed

!

r_state => pointer_to_ran_state()
seed = r_state%seed

end subroutine ran_seed_get 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_default_state (set_state, get_state)
!
! Routine to set or get the state of the default random number generator.
! See the ran_seed_put documentation for more details
!
! Input:
!   set_state -- random_state_struct, optional: State to set the default generator to.
!
! Output:
!   get_state -- random_state_struct, optional: Returns the state of the default generator.
!-

subroutine ran_default_state (set_state, get_state)

implicit none

type (random_state_struct), optional :: set_state, get_state
type (random_state_struct), pointer :: state_ptr

!

state_ptr => pointer_to_ran_state()

if (present(set_state)) state_ptr = set_state
if (present(get_state)) get_state = state_ptr

end subroutine ran_default_state

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_scalar (harvest, ran_state, index_quasi)
!
! Routine to return a random number uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!
! Note: The index_quasi argument is used internally for the quasi-random number generator.
!-

subroutine ran_uniform_scalar (harvest, ran_state, index_quasi)

implicit none

type (random_state_struct), pointer :: r_state
type (random_state_struct), optional, target :: ran_state

real(rp), intent(out) :: harvest

integer(kr4b) k, ix_q
integer, optional :: index_quasi


integer(kr4b), parameter :: ia = 16807
integer(kr4b), parameter :: iq = 127773, ir = 2836

character :: r_name = 'ran_uniform_scalar'

! If r_state%iy < 0 then the random number generator has never been initialized.

r_state => pointer_to_ran_state(ran_state)
if (r_state%iy < 0) call ran_seed_put(r_state%seed)

! quasi-random

if (r_state%engine == quasi_random$) then
  ix_q = integer_option(1, index_quasi)
  if (ix_q == 1) call super_sobseq (r_state%x_sobseq, ran_state)
  if (ix_q > sobseq_maxdim) then
    call out_io (s_error$, r_name, 'NUMBER OF DIMENSIONS WANTED IS TOO LARGE!')
    if (global_com%exit_on_error) call err_exit
  endif
  harvest = r_state%x_sobseq(ix_q)
  return
endif

! Pseudo-random
! Marsaglia shift sequence with period 2^32 - 1.

r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 13)) 
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, -17))
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 5))
k = r_state%iy/iq         ! Park-Miller sequence by Schrage's method,
r_state%iy = ia*(r_state%iy - k*iq) - ir * k       ! period 2^31-2.

if (r_state%iy < 0) r_state%iy = r_state%iy + im_nr_ran

! Combine the two generators with masking to ensure nonzero value.

harvest = r_state%am * ior(iand(im_nr_ran, ieor(r_state%ix, r_state%iy)), 1) 

end subroutine ran_uniform_scalar

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_vector (harvest, ran_state)
!
! Routine to return a vector of random numbers uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!-

subroutine ran_uniform_vector (harvest, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state

real(rp), intent(out) :: harvest(:)
integer i

!

do i = 1, size(harvest)
  call ran_uniform_scalar (harvest(i), ran_state, i)
enddo

end subroutine ran_uniform_vector

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine super_sobseq (x, ran_state)
!
! Routine patterened after sobseq in Numerical Recipes.
! Difference is that this version has an argument for the internal state.
!
! Input:
!   ran_state -- random_state_struct, optional: Generator state.
!                     See the ran_seed_put documentation for more details.
!
! Output: 
!   x(:)      -- real(dp): Random vector.
!   ran_state -- random_state_struct, optional: Generator state.
!-

subroutine super_sobseq (x, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

real(rp), dimension(:), intent(out) :: x
real(dp), parameter :: fac=1.0_dp/2.0_dp**sobseq_maxbit

integer(i4_b) :: im, j
integer(i4_b), parameter :: ip(sobseq_maxdim) = [0,1,1,2,1,4], mdeg(sobseq_maxdim) = [1,2,3,3,4,4] 
! Note: If sobseq_maxbit or sobseq_maxdim are changed, iv needs to be recomputed by
! running the original sobseq routine.
integer(i4_b), parameter :: iv(sobseq_maxdim*sobseq_maxbit) = [ &
                      536870912, 536870912, 536870912, 536870912, 536870912, &
                      536870912, 805306368, 268435456, 805306368, 805306368, &
                      268435456, 268435456, 671088640, 939524096, 939524096, &
                      402653184, 402653184, 671088640, 1006632960, 738197504, &
                      335544320, 1006632960, 872415232, 603979776, 570425344, &
                      436207616, 234881024, 167772160, 838860800, 100663296, &
                      855638016, 1023410176, 721420288, 285212672, 150994944, &
                      385875968, 713031680, 562036736, 411041792, 713031680, &
                      763363328, 1031798784, 1069547520, 331350016, 616562688, &
                      566231040, 88080384, 465567744, 538968064, 975175680, &
                      920649728, 853540864, 941621248, 497025024, 808452096, &
                      756023296, 1062207488, 489684992, 605028352, 198180864, &
                      673710080, 431489024, 381157376, 952631296, 706215936, &
                      898105344, 1010565120, 1072431104, 258736128, 208928768, &
                      1026818048, 804519936, 572653568, 540672000, 771883008, &
                      316801024, 531759104, 864944128, 858980352, 271384576, &
                      453181440, 758317056, 206110720, 954400768, 715816960, &
                      941195264, 545488896, 550076416, 361594880, 238256128, &
                      1073725440, 742375424, 817971200, 813154304, 559235072, &
                      590921728, 536879104, 438312960, 954261504, 417505280, &
                      301998080, 328081408, 805318656, 1024462848, 340963328, &
                      1009913856, 419434496, 685969408, 671098880, 565721088, &
                      238651392, 172697600, 897587200, 640919552, 1006648320, &
                      334244864, 732843008, 297131008, 826291200, 121168896, &
                      570434048, 976886272, 417426944, 704744960, 169882112, &
                      361637376, 855651072, 757939456, 609285376, 553894656, &
                      756025600, 1071843072, 713042560, 429498752, 909831040, &
                      847291520, 127413632, 464754048, 1069563840, 1073730496, &
                      1068349120, 499194688, 947127616, 486080448, 538976288, &
                      536877600, 383778848, 954376224, 663894048, 137247648, &
                      808464432, 268451088, 256901168, 204607536, 676930576, &
                      875760848, 673720360, 939532728, 783810616, 306915352, &
                      1066773016, 775659528, 1010580540, 738202604, 460062740, &
                      766893116, 476151092, 856477748, 572662306, 436222522, &
                      537001998, 536972810, 229785522, 1000343598, 858993459, &
                      1023421741, 805503019, 805552913, 357112905, 214971443]

character(16), parameter :: r_name = 'super_sobseq'

! Calc

r_state => pointer_to_ran_state(ran_state)

im = r_state%in_sobseq
do j = 1, sobseq_maxbit
  if (.not. btest(im,0)) exit
  im = im/2
end do

if (j > sobseq_maxbit) then
  call out_io (s_fatal$, r_name, 'SOBSEQ_MAXBIT TOO SMALL')
  if (global_com%exit_on_error) call err_exit
  return
endif

im = (j-1) * sobseq_maxdim
j = min(size(x), sobseq_maxdim)
r_state%ix_sobseq(1:j) = ieor(r_state%ix_sobseq(1:j),iv(1+im:j+im))
x(1:j) = r_state%ix_sobseq(1:j) * fac
r_state%in_sobseq = r_state%in_sobseq + 1

end subroutine super_sobseq

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function pointer_to_ran_state(ran_state, ix_thread) result (ran_state_ptr)
!
! Routine to point to the appropriate state structure for generating random numbers
!
! Input:
!   ran_state     -- random_state_struct, optional: Point to this if present.
!                      Otherwise point to the global saved state.
!   ix_thread     -- integer, optional: Thread index.
!
! Output:
!   ran_state_ptr -- random_state_struct, pointer: Pointer to the appropriate state.
!-

function pointer_to_ran_state(ran_state, ix_thread) result (ran_state_ptr)

!$ use omp_lib, only: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: ran_state_ptr

integer, optional :: ix_thread
integer nt

!

if (present(ran_state)) then
  ran_state_ptr => ran_state
  return
endif

!$ if (.not. thread_state_allocated) call allocate_thread_states()
!$ nt = integer_option(OMP_GET_THREAD_NUM(), ix_thread)
!$ ran_state_ptr => thread_ran_state(nt)
!$ return

ran_state_ptr => ran_state_save

end function pointer_to_ran_state

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine allocate_thread_states()
!
! Routine to allocate random number state structures when openMP is used.
!-

subroutine allocate_thread_states()

!$ use omp_lib, only: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

implicit none

integer nt, max_t

!

!$ if (thread_state_allocated) return
!$ nt = OMP_GET_THREAD_NUM()
!$ max_t = OMP_GET_MAX_THREADS()
!$ if (nt == 0) allocate (thread_ran_state(0:max_t-1))
!$OMP BARRIER
!$ if (nt == 0) thread_state_allocated = .true.
!$OMP BARRIER

end subroutine allocate_thread_states

end module
