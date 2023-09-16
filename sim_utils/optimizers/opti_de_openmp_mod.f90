module opti_de_openmp_mod

use opti_de_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function opti_de_openmp (v_best, generations, population, merit_func, 
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

function opti_de_openmp (v_best, generations, population, merit_func, v_del, status) result (best_merit)

!$ use omp_lib

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
  real(rp), allocatable :: vec(:)
  real(rp) merit
  integer status
end type

type (solution_struct) :: solution(population) 
type (solution_struct) :: try(population)

real(rp) :: v1(size(v_best)), r_vec(size(v_best))
real(rp) :: F, r_ran, r_pop(population)

integer status, iter_count
integer :: i, ii, ng, nd, i_best, k, kk, n_vec
integer :: i1(population), i2(population), i3(population)
integer :: i4(population), i5(population)
integer p1, p2, p3, p4, p5

logical this_better_merit, this_best_merit

character(*), parameter :: r_name = 'opti_de_openmp'

! Initialize

nd = size(v_best)
status = 0
iter_count = 0

n_vec = 3
if (opti_de_param%use_2nd_diff) n_vec = 5

! Error check

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

! Find the best of the initial population

do i = 1, population
  allocate (solution(i)%vec(nd), try(i)%vec(nd))

  if (i == 1) then
    solution(1)%vec = v_best
  else
    call random_number (v1)
    solution(i)%vec = v_best + v_del * (2*v1 - 1)
  endif
enddo

!$OMP parallel do
do i = 1, population
  solution(i)%merit = merit_func(solution(i)%vec, solution(i)%status, iter_count)
enddo
!$OMP end parallel do

do i = 1, population
  if (i == 1) then
    i_best = 1
  elseif (opti_de_param%minimize_merit) then
    if (solution(i)%merit < best_merit) i_best = i
  else
    if (solution(i)%merit > best_merit) i_best = i
  endif

  best_merit = solution(i_best)%merit
  v_best = solution(i_best)%vec

  if (solution(i)%status /= 0) return
enddo

!---------------------------------
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

    ! Construct the perturbed vector

    v1 = solution(p1)%vec + opti_de_param%l_best * (v_best - solution(p1)%vec) + F * (solution(p2)%vec - solution(p3)%vec)
    if (opti_de_param%use_2nd_diff) v1 = v1 + F * (solution(p4)%vec - solution(p5)%vec)

    ! Perform crossover

    try(i)%vec = solution(i)%vec

    if (opti_de_param%binomial_cross) then
      call random_number (r_vec)
      where (r_vec < opti_de_param%CR) try(i)%vec = v1

    else  ! exponential crossover
      call random_number(r_ran) 
      ii = r_ran * nd + 1
      do k = 1, nd
        call random_number(r_ran) 
        if (r_ran > opti_de_param%CR) exit
        ii = mod(ii, nd) + 1
        try(i)%vec(ii) = v1(ii)
        ii = ii + 1
      enddo
    endif
  enddo

  ! Compute merit

  !$OMP parallel do
  do i = 1, population
    try(i)%merit = merit_func (try(i)%vec, try(i)%status, iter_count)
  enddo
  !$OMP end parallel do

  ! Compare solutions

  do i = 1, population
    if (opti_de_param%minimize_merit) then
      this_better_merit = (try(i)%merit < solution(i)%merit)
      this_best_merit   = (try(i)%merit < best_merit)
    else
      this_better_merit = (try(i)%merit > solution(i)%merit)
      this_best_merit   = (try(i)%merit > best_merit)
    endif

    if (this_better_merit) solution(i) = try(i)

    if (this_best_merit) then
      v_best = try(i)%vec
      best_merit = try(i)%merit
    endif

    if (try(i)%status /= 0) return
  enddo
enddo

end function opti_de_openmp

end module
