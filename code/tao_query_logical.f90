!+
! Subroutine tao_query_logical (true_sym, false_sym, prompt, logic)
!
! Subroutine to prompt for a logical answer.
! Default is always original state of logical.
!
! Input:
!   true_sym  -- Character(*): 
!   false_sym -- Character(*):
!   prompt    -- Character(*):
!   logic     -- Logical: Default.
!
! Output
!   logic     -- Logical: 
!-

subroutine tao_query_logical (true_sym, false_sym, prompt, logic)

implicit none

logical logic
character(*) true_sym, false_sym, prompt
character(40) ans, str

!

do

  if (logic) then
    str = true_sym
  else
    str = false_sym
  endif

  write (*, '(10a)', advance = 'NO'), trim(prompt), ' (', trim(true_sym), &
      '/', trim(false_sym), ') <', trim(str), '> ' 
  read  (*, '(a)') ans

  if (len_trim(ans) == 0) return
  if (ans == true_sym) then
    logic = .true.
    return
  elseif (ans == false_sym) then
    logic = .false.
    return
  endif

enddo

end subroutine