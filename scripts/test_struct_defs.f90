type wake_lr_struct
  integer(8) :: aaa
  integer :: bbb
  character(200) :: file = ''
  real(rp) :: t_ref = 0             ! time reference value for computing the wake amplitude.
                                    !  This is used to prevent value overflow with long trains.
  real(rp) :: freq_spread = 0       ! Random frequency spread of long range modes.
end type

type wake_struct
  type (wake_lr_struct) :: sr = wake_sr_struct('', null(), null(), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, .true.) ! Short-range wake
end type

