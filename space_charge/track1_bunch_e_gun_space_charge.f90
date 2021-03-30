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
use nr, only: indexx

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (coord_struct), pointer :: p

real(rp) s0, dt_step, t_now, t_end
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

! If starting from the cathode, init all elements to be pre_born. The element to track through 
! could be a slice or super_slave so must check if The beginning of the element is at the cathode.

s0 = branch%ele(0)%s
if (ele%s_start == s0) then
  where (bunch%particle%location == upstream_end$) bunch%particle%state = pre_born$ 
else
  ! Track particles from current position to equal time.
endif

n_pre_born = 0

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  call track1_preprocess (p, ele, branch%param, err, finished, radiation_included)
  if (finished) call err_exit   ! I don't know what to do with this!
  if (p%state == pre_born$) then
    t_emit(i) = p%t
    n_pre_born = n_pre_born + 1
  endif
enddo

call indexx(t_emit, ix_t_emit)

! Track

do
  t_end = t_end + ele%value(dt_max$)

  if (n_pre_born == 0) then
    call track_bunch_time(branch%lat, bunch, t_end)
    cycle
  endif

  do n = 1, min(n_pre_born, n_emit_max)
    if (t_emit(ix_t_emit(n)) > t_end) exit
  enddo

  n = n - 1

  if (n < 1) then ! No particle emitted here
    call track_bunch_time(branch%lat, bunch, t_end)
    cycle
  endif

  is_tracked = .false.

  do i = 1, n
    j = ix_t_emit(i)
    p => bunch%particle(j)
    p%state = alive$
    call track1_time_runge_kutta (p, ele, branch%param, p, err, t_end = t_end)
    is_tracked(j) = .true.
  enddo

  ix_t_emit(1:n_pre_born-n) = ix_t_emit(n+1:n_pre_born)  
  n_pre_born = n_pre_born - n

  do i = 1, size(bunch%particle)
    if (is_tracked(i)) cycle
    p => bunch%particle(i)
    call track1_time_runge_kutta (p, ele, branch%param, p, err, t_end = t_end)
  enddo    

  ! Apply SC kick
  ! Need to apply SC kick to newly born particles proportional to the time from birth to the end of the time step.
  call err_exit  ! Place holder

enddo

!

do i = 1, size(bunch%particle)
  call track1_postprocess (bunch%particle(i), ele, branch%param, bunch%particle(i))
  if (finished) call err_exit   ! I don't know what to do with this.
enddo

end subroutine track1_bunch_e_gun_space_charge


! Idea: Can apply SC kick to 
! Question: How to determine dt_step?
! Question: When to end tracking?

