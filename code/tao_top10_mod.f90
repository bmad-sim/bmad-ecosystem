module tao_top10_mod


use tao_struct
use tao_interface
use cesr_utils
use tao_dmerit_mod

! structure for making lists of the biggest contributors to the merit function.

type tao_top10_struct
  character(16) name   ! name of contributor
  real(rp) value       ! contribution to the merit function
  integer index        ! index of contributor.
  logical valid        ! valid entry?
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_top10_print ()
!
! Routine to print out the top10 contributors to the merit function.
!
! Input:
!-

subroutine tao_top10_print ()

implicit none

type (tao_top10_struct) top_merit(10)
type (tao_top10_struct) top_dmerit(10)
type (tao_top10_struct) top_delta(10)

real(rp) delta, a_max, merit
integer i, j, n, nl, nu

character(18) name
character(100) fmt, lines(20)
character(20) :: r_name = 'tao_top10_print'

! tao_merit also calculates the contrribution of the individual
! variables and data to the merit function.

merit = tao_merit()
call tao_dmerit_calc ()

! top_merit stores the top contributors to the merit function.
! top_dmerit stores the top dmerit/dvar values
! top_delta stores the top |var_model - var_design| 

top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '
top_merit(:)%value = 0; top_merit(:)%index = 0
top_dmerit(:)%valid = .false.; top_dmerit(:)%name = ' '
top_dmerit(:)%value = 0; top_dmerit(:)%index = 0
top_delta(:)%valid  = .false.; top_delta(:)%name  = ' '
top_delta(:)%value = 0; top_delta(:)%index = 0

nu = size(s%u)
do i = 1, nu
  do j = 1, size(s%u(i)%data)
    if (.not. s%u(i)%data(j)%useit_opt) cycle
    name = s%u(i)%data(j)%data_type
    if (nu > 1) write (name, '(2a, i0)') trim(name), ';', i
    call tao_to_top10 (top_merit, s%u(i)%data(j)%merit, name, &
                                                s%u(i)%data(j)%ix_d1, 'max')
  enddo
enddo


do j = 1, size(s%var)
  if (.not. s%var(j)%useit_opt) cycle
  name = s%var(j)%v1%name
  call tao_to_top10 (top_merit, s%var(j)%merit, name, s%var(j)%ix_v1, 'max')
  call tao_to_top10 (top_dmerit, s%var(j)%dmerit_dvar, name, &
                                                      s%var(j)%ix_v1, 'max')
  delta = s%var(j)%model_value - s%var(j)%design_value
  call tao_to_top10 (top_delta, delta, name, s%var(j)%ix_v1, 'max')
enddo

! write results


a_max = max(1.1, maxval(abs(top_delta(:)%value)))
n = max(0, 6 - int(log10(a_max)))

write (fmt, '(a, i1, a)') &
   '((1x, a10, i5, f11.1, 3x), (a8, i5, 1pe12.3, 3x), (a8, i5, 0pf11.', n, '))'


nl = 0
lines(nl+1) = ' '
lines(nl+2) = &
  '       Top10 merit          |     Top10 derivative      |      Top10 delta'
lines(nl+3) = &
  ' Name         ix      Value | Name       ix  Derivative | Name       ix     delta'
nl = nl + 3

do i = 1, 10
  nl = nl + 1
  write (lines(nl), fmt) &
      top_merit(i)%name,  top_merit(i)%index,  top_merit(i)%value, &
      top_dmerit(i)%name, top_dmerit(i)%index, top_dmerit(i)%value,  &
      top_delta(i)%name,  top_delta(i)%index,  top_delta(i)%value
enddo

nl = nl + 1
write (lines(nl), *) 'Merit:  ', merit

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_top10 (top10, value, name, c_index, order)
!
! Routine to order the largest contributors to the merit function in
! a list. Call this routine for each contributor.
!
! Note: Before first calling this routine set:
!   top10(:)%valid = .false.
!
! Input:
!   value   -- Real(rp): value of the contributor.
!   name    -- Character(16): Name of the contributor..
!   c_index -- Integer: Index of the contributor.
!   order   -- Character(16): Ordering of the list. Possibilities are:
!                 'max'     -- #1 has the maximum value.
!                 'min'     -- #1 has the minimum value.
!                 'abs_max' -- #1 has the maximum aplitude.
!                 'abs_min' -- #1 has the maximum aplitude.
!
! Output:
!   top10(:) -- Tao_top10_struct: List of top contributors.
!                 Note that the list is not limited to 10 entries.
!-

subroutine tao_to_top10 (top10, value, name, c_index, order)

implicit none

type (tao_top10_struct) top10(:)

integer c_index, ix, n
real(rp) value

character(*) name, order
character(20) :: r_name = 'tao_to_top10'

! Find where in list the current contributor is.

n = size(top10)
do ix = n, 1, -1
  if (.not. top10(ix)%valid) cycle
  select case (order)
  case ('max')
    if (value < top10(ix)%value) exit
  case ('min')
    if (value > top10(ix)%value) exit
  case ('abs_max')  
    if (abs(value) < abs(top10(ix)%value)) exit
  case ('abs_min')  
    if (abs(value) > abs(top10(ix)%value)) exit
  case default
    call out_io (s_abort$, r_name, 'BAD "ORDER" ARGUMENT: ' // order)
  end select
enddo

ix = ix + 1          ! place to put current contributor.
if (ix > n) return   ! not big enough to be in list.

! Move the people below the current contributor down to make room and
! then put the contributor in.

top10(ix+1:n) = top10(ix:n-1) 
top10(ix) = tao_top10_struct(name, value, c_index, .true.)

end subroutine


end module
