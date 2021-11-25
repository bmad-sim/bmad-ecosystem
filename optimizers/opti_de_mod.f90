module opti_de_mod

use sim_utils

type opti_de_param_struct
  real(rp) :: CR = 0.8                    ! Crossover probability
  real(rp) :: F  = 0.8                    ! Mixing number
  real(rp) :: l_best = 0.0                ! Percentage of best vector.
  logical :: use_2nd_diff   = .false.     ! use F * (x_4 - x_5) term
  logical :: binomial_cross = .false.
  logical :: randomize_F    = .false.
  logical :: minimize_merit = .true.      ! Alternative is to maximize.
end type

type (opti_de_param_struct), save :: opti_de_param

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function opti_de (v_best, generations, population, merit_func, 
!                                           v_del, status) result (best_merit)
!    
! Differential Evolution for Optimal Control Problems.
! This optimizer is based upon the work of Storn and Price: 
!   R. Storn, and K. V. Price, 
!   "Minimizing the real function of the ICEC'96 contest by differential evolution"
!   IEEE conf. on Evolutionary Computation, 842-844 (1996).
!
! Both arrays: v_best(:) and v_del(:) need to have the same size.
!
! The "perturbed vector" is
!   v = x_1 + l_best * (x_best - x_1) + F * (x_2 - x_3) + F * (x_4 - x_5)
! The last term F * (x_4 - x_5) is only used if opti_de_param%use_2nd_diff = T.
!
! The crossover can be either "Exponential" or "Binary". 
! Exponential crossover is what is described in the paper.
! With Exponential crossover the crossover parameters form a contiguous block
! and the average number of crossover parameters is approximately
!     average crossovers ~ min(D, CR / (1 - CR))
! where D is the total number of parameters.
!
! With Binary crossover the probability of crossover of a parameter is 
! uncorrelated with the probability of crossover of any other parameter and
! the average number of crossovers is
!     average crossovers = D * CR
!
! The parameters used for the DE are in opti_de_param:
!                                 Default
!   opti_de_param%CR               0.8    ! Crossover Probability.
!   opti_de_param%F                0.8    !
!   opti_de_param%l_best           0.0    ! Percentage of best solution used.
!   opti_de_param%binomial_cross   False  ! IE: Default = Exponential.
!   opti_de_param%use_2nd_diff     False  ! use F * (x_4 - x_5) term
!   opti_de_param%randomize_F      False  ! 
!   opti_de_param%minimize_merit   True   ! F => maximize the Merit func.
!
! opti_de_param%randomize_F = True means that the F that is used for a given 
! generation  is randomly choisen to be within the range 
! [0, 2*F] with average F.
!
! The Merit function prototype is:
!   function merit_func (vec, status, iter_count) result (merit)
!     use precision_def
!     real(rp) vec(:)     ! Input: trial solution.
!     integer status      ! Output: Set non-zero to terminate opti_de.
!     integer iter_count  ! Input/Output: Initially set to zero by opti_de.
!                         !   To be used by merit_func for display purposes, etc.
!     real(rp) merit      ! Output: Merit value corresponting to vec.
!   end function
!
! Input:
!   generations -- Integer: Number of generations to evolve.
!   population  -- Integer: Number of trial solutions in any generation.
!   merit_func  -- Function: See the prototype above.
!   v_best(:)   -- Real(rp): Starting solution.
!   v_del(:)    -- Real(rp): Deltas used for starting population.
!
! Output:
!   v_best(:)  -- Real(rp): Final solution with best merit value.
!   opti_de    -- Real(rp): Best merit value with final solution. 
!   status     -- Integer:  Setting of the status argument from merit_func.
!   best_merit -- Real(rp): Merit at minimum.
!-

function opti_de (v_best, generations, population, merit_func, v_del, status) result (best_merit)

use indexer_mod

implicit none

integer, intent(in) :: generations, population
real(rp), intent(in) :: v_del(:)
real(rp), intent(out) :: v_best(:)
real(rp) :: best_merit

interface
  function merit_func (vec, status, iter_count) result (merit)
    use precision_def
    real(rp) vec(:)               
    integer status 
    integer iter_count
    real(rp) merit
  end function
end interface

type solution_struct
  real(rp), pointer :: vec(:)
  real(rp) merit
end type

type (solution_struct) :: trial(population) 

real(rp) :: v1(size(v_best)), v2(size(v_best))
real(rp) :: F, merit, r_ran, r_pop(population), r_vec(size(v_best))

integer status, iter_count
integer :: i, ii, nd, ng, i_best, k, kk, n_vec
integer :: i1(population), i2(population), i3(population)
integer :: i4(population), i5(population)
integer p1, p2, p3, p4, p5

logical this_better_merit, this_best_merit

character(*), parameter :: r_name = 'opti_de'

! Error check

nd = size(v_best)

if (size(v_del) /= nd) then
  call out_io (s_error$, r_name, 'ARRAY SIZES NOT THE SAME!', &
                               'FOR V_DEL, AND V_BEST: \2i6\ ', i_array=[size(v_del), nd])
  if (global_com%exit_on_error) call err_exit
endif

if (population < 4) then
  call out_io (s_error$, r_name, 'POPULATION MUST BE AT LEAST 4! IT IS: \i0\ ', population)
  if (global_com%exit_on_error) call err_exit
endif

if (population < 6 .and. opti_de_param%use_2nd_diff) then
  call out_io (s_error$, r_name, 'POPULATION MUST BE AT LEAST 6 WITH USE_2ND_DIFF! IT IS: \i0\ ', population)
  if (global_com%exit_on_error) call err_exit
endif

! Initialize

status = 0
iter_count = 0

n_vec = 3
if (opti_de_param%use_2nd_diff) n_vec = 5

! Find the best of the initial population

do i = 1, population

  allocate (trial(i)%vec(nd))

  if (i == 1) then
    trial(1)%vec = v_best
  else
    call random_number (v1)
    trial(i)%vec = v_best + v_del * (2*v1 - 1)
  endif

  trial(i)%merit = merit_func(trial(i)%vec, status, iter_count)

  if (i == 1) then
    i_best = 1
  elseif (opti_de_param%minimize_merit) then
    if (trial(i)%merit < best_merit) i_best = i
  else
    if (trial(i)%merit > best_merit) i_best = i
  endif

  best_merit = trial(i_best)%merit
  v_best = trial(i_best)%vec

  if (status /= 0) return

enddo

! Loop over all generations

do ng = 1, generations

  if (opti_de_param%randomize_F) then
    call random_number (r_ran)
    F = 2 * opti_de_param%f * r_ran
  else
    F = opti_de_param%f
  endif

  ! i1, ... i5 gives the indexes for constructing the perturbed vector.

  call random_number (r_pop)
  call indexer (r_pop, i1)   ! i1 is a random permutation

  call random_number (r_pop)
  call indexer (r_pop, i2)   ! i2 is a random permutation

  call random_number (r_pop)
  call indexer (r_pop, i3)   ! i3 is a random permutation

  if (opti_de_param%use_2nd_diff) then
    call random_number (r_pop)
    call indexer (r_pop, i4)   ! i4 is a random permutation

    call random_number (r_pop)
    call indexer (r_pop, i5)   ! i5 is a random permutation
  endif

  ! Loop over the entire population

  do i = 1, population
    do k = 1, population
      kk = mod (k + i, population) + 1
      p1 = i1(kk) 
      if (p1 /= i) exit
    enddo

    do k = 1, population
      kk = mod (k + i, population) + 1
      p2 = i2(kk) 
      if (p2 /= i .and. p2 /= p1) exit
    enddo

    do k = 1, population
      kk = mod (k + i, population) + 1
      p3 = i3(kk) 
      if (p3 /= i .and. p3 /= p1 .and. p3 /= p2) exit
    enddo

    if (opti_de_param%use_2nd_diff) then
      do k = 1, population
        kk = mod (k + i, population) + 1
        p4 = i4(kk) 
        if (p4 /= i .and. p4 /= p1 .and. p4 /= p2 .and. p4 /= p3) exit
      enddo

      do k = 1, population
        kk = mod (k + i, population) + 1
        p5 = i5(kk) 
        if (p5 /= i .and. p5 /= p1 .and. p5 /= p2 .and. p5 /= p3 .and. p5 /= p4) exit
      enddo
    endif

    ! Find random vectors

    ! Construct the perturbed vector

    v1 = trial(p1)%vec + opti_de_param%l_best * (v_best - trial(p1)%vec) + F * (trial(p2)%vec - trial(p3)%vec)
    if (opti_de_param%use_2nd_diff) v1 = v1 + F * (trial(p4)%vec - trial(p5)%vec)

    ! Perform crossover

    v2 = trial(i)%vec

    if (opti_de_param%binomial_cross) then
      call random_number (r_vec)
      where (r_vec < opti_de_param%CR) v2 = v1

    else  ! exponential crossover
      call random_number(r_ran) 
      ii = r_ran * nd + 1
      do k = 1, nd
        call random_number(r_ran) 
        if (r_ran > opti_de_param%CR) exit
        ii = mod(ii, nd) + 1
        v2(ii) = v1(ii)
        ii = ii + 1
      enddo
    endif

    ! Compare

    merit = merit_func (v2, status, iter_count)

    this_better_merit = .false.
    this_best_merit = .false.

    if (opti_de_param%minimize_merit) then
      if (merit < trial(i)%merit) this_better_merit = .true.
      if (merit < best_merit) this_best_merit = .true.
    else
      if (merit > trial(i)%merit) this_better_merit = .true.
      if (merit > best_merit) this_best_merit = .true.
    endif

    if (this_better_merit) then
      trial(i)%vec = v2
      trial(i)%merit = merit
    endif

    if (this_best_merit) then
      v_best = v2
      best_merit = merit
    endif

    if (status /= 0) return
  enddo

enddo

! Cleanup

do i = 1, population
  deallocate (trial(i)%vec)
enddo

end function opti_de

end module
