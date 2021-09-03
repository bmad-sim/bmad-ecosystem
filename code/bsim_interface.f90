module bsim_interface
  
  interface
     subroutine set_tune3(lat, target_tunes, everything_ok, verbose_in, write_out_in, regex_in)
       use bmad_struct, only: lat_struct, rp
       implicit none
       type(lat_struct), target :: lat
       real(rp) target_tunes(3)
       logical everything_ok
       logical, optional :: verbose_in
       integer, optional :: write_out_in
       character(40), optional :: regex_in
       
     end subroutine set_tune3
  end interface
  
integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module bsim_interface


