!+
! Module tao_common
!
! module for the common block for the optimizers.
!-

module tao_common

use tao_mod

type tao_alias_struct
  character(24) :: name
  character(100) :: string
end type

type tao_common_struct
  type (tao_alias_struct) alias(100)
  logical opti_init        ! init needed?
  logical opti_at_limit    ! Variable at limit?
  character(40) cmd_arg(9) ! Command file arguments.
  integer :: n_alias = 0
end type

type (tao_super_universe_struct), pointer :: s_com
type (tao_common_struct) tao_com

end module
