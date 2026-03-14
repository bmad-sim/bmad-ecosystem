subroutine set_env(env_name, env_value, err_flag)

use iso_c_binding
implicit none

interface
  function setenv(name, value, overwrite) bind(C, name="setenv")
    use iso_c_binding
    integer(c_int) :: setenv
    character(kind=c_char), intent(in) :: name(*), value(*)
    integer(c_int), value :: overwrite
  end function
end interface

character(*) env_name, env_value
logical err_flag
integer(c_int) :: ret

!

ret = setenv(trim(env_name) // c_null_char, trim(env_value) // c_null_char, 1_c_int)
err_flag = (ret /= 0)

end subroutine
