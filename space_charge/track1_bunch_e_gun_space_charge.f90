!+
! Subroutine track1_bunch_e_gun_space_charge (bunch, ele, err)
!
! Subroutine to track a bunch of particles through an e_gun.
!
! Input:
!   bunch     -- bunch_struct: Starting bunch position.
!   ele       -- Ele_struct: E_gun element to track through. Must be part of a lattice.
!
! Output:
!   bunch     -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch_e_gun_space_charge (bunch, ele, err)

use bmad, dummy => track1_bunch_e_gun_space_charge

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (coord_struct), pointer :: p

real(rp) s0, dt_step, t_now, t_end, dt_max
real(rp) :: t_emit(size(bunch%particle))
integer i, j, n, n_pre_born, n_emit_max
integer :: ix_t_emit(size(bunch%particle))

logical err, finished, radiation_included
logical :: is_tracked(size(bunch%particle))

character(*), parameter :: r_name = 'track1_bunch_e_gun_space_charge'

!

branch => pointer_to_branch(ele)
t_emit = 1e30_rp  ! Something large
n_emit_max = max(1, nint(ele%value(emit_fraction$) * size(bunch%particle)))

! If a particle is marked as still being in a preceding element then update. 
! Also if starting from the cathode, init all elements to be pre_born. The element to track through 
! could be a slice or super_slave so must check if The beginning of the element is at the cathode.

s0 = branch%ele(0)%s
t_now = 1e30_rp  ! Something large
n_pre_born = 0

do i = 1, size(bunch%particle)
  p => bunch%particle(i)

  if (p%ix_ele < ele%ix_ele) then
    p%ix_ele = ele%ix_ele
    p%location = upstream_end$
  endif

  if (ele%s_start == s0 .and. p%location == upstream_end$) p%state = pre_born$

  if (p%state == pre_born$) then
    t_emit(i) = p%t
    n_pre_born = n_pre_born + 1
  endif

  call track1_preprocess (p, ele, branch%param, err, finished, radiation_included)

  t_now = min(p%t, t_now)
enddo

call indexer(t_emit, ix_t_emit)

dt_max = ele%value(dt_max$)
if (dt_max == 0) then
  dt_max = 1e-10   !!! Arbitrary for testing!!
  call out_io (s_warn$, r_name, 'Element: ' // ele%name, 'Does not have dt_max set!')
endif

! Track

if (.false.) then
do
  t_end = t_now + dt_max

  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    if (p%state == pre_born$ .and. p%t <= t_end) p%state = alive$
  enddo 

  call track_bunch_time(branch%lat, bunch, t_end, ele%s)

  ! Apply SC kick
  ! Need to apply SC kick to newly born particles proportional to the time from birth to the end of the time step.

  t_now = t_end

  finished = .true.
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    if (p%state == pre_born$ .or. (p%s < ele%s .and. p%state == alive$)) then
      finished = .false.
      exit
    endif
  enddo 
  if (finished) exit
enddo
endif

!

do i = 1, size(bunch%particle)
  call track1_postprocess (bunch%particle(i), ele, branch%param, bunch%particle(i))
enddo

call out_io (s_error$, r_name, 'E_GUN TRACKING WITH CATHODE NOT YET IMPLEMENTED!')

end subroutine track1_bunch_e_gun_space_charge


! Idea: Can apply SC kick to 
! Question: How to determine dt_step?
! Question: When to end tracking?

