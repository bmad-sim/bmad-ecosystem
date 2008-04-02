subroutine logic_get( true_ans, false_ans, question, ans )

  implicit none
  character(*) :: true_ans, false_ans
  character(1) :: true_ans_a, false_ans_a
  character(*) :: question
  logical :: ans
  character(80) :: temp
  integer :: len

  write (*, "(a)") question
  accept "(a)", temp

  true_ans_a = true_ans(1:1)
  false_ans_a = false_ans(1:1)
  call string_trim(temp, temp, len)
  call str_upcase(true_ans_a, true_ans_a)
  call str_upcase(false_ans_a, false_ans_a)
  call str_upcase(temp, temp)

  if (temp(1:1) == true_ans_a) then
    ans = .true.
!  else if (temp(1:1) == false_ans_a) then
!    ans = .false.
  else
! Accept anything else as false
!    write (*, "(a)") "Unknown answer!!!  Assuming false for now!"
    ans = .false.
  end if

end subroutine logic_get
