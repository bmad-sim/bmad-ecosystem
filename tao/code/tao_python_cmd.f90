!+
! Subroutine tao_python_cmd (input_str)
!
! This routine has been replaced by tao_pipe_cmd.
! This routine is keep around for backwards compatibility.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_python_cmd (input_str)

implicit none
character(*) input_str

!

call tao_pipe_cmd(input_str)

end subroutine
