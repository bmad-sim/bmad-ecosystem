subroutine logic_get( true_ans, false_ans, question, ans )

  implicit none
  character(1) :: true_ans, false_ans
  character(*) :: question
  logical :: ans
  character(80) :: temp
  integer :: len

  write (*, "(a)") question
  accept "(a)", temp

  call string_trim(temp, temp, len)
  call str_upcase(true_ans, true_ans)
  call str_upcase(false_ans, false_ans)
  call str_upcase(temp, temp)

  if (temp(1:1) == true_ans(1)) then
    ans = .true.
  else if (temp(1:1) == false_ans(1)) then
    ans = .false.
  else
    write (*, "(a)") "Unknown answer!!!  Assuming false for now!"
    ans = .false.
  end if

end subroutine logic_get
