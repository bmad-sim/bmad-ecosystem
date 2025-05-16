module modulo2_mod

use precision_def

!+
! Function modulo2 (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Input:
!   x    -- Real(sp), Real(rp), or Integer
!   amp  -- Real(sp), Real(rp), or Integer: Must be positive.
!
! Output:
!   mod2 -- Real(sp), Real(rp), or Integer: Result
!-

interface modulo2
  module procedure modulo2_sp 
  module procedure modulo2_dp
  module procedure modulo2_qp
  module procedure modulo2_int
end interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_sp (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Input:
!   x    -- Real(sp): 
!   amp  -- Real(sp): Must be positive.
!
! Output:
!   mod2 -- Real(sp): Result
!-

elemental function modulo2_sp (x, amp) result (mod2)


  use precision_def

  implicit none

  real(sp), intent(in) :: x, amp
  real(sp) mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_dp (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Input:
!   x    -- Real(rp): 
!   amp  -- Real(rp): Must be positive.
!
! Output:
!   mod2 -- Real(rp): Result
!-

elemental function modulo2_dp (x, amp) result (mod2)


  use precision_def

  implicit none

  real(rp), intent(in) :: x, amp
  real(rp) mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_qp (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Input:
!   x    -- Real(rp): 
!   amp  -- Real(rp): Must be positive.
!
! Output:
!   mod2 -- Real(rp): Result
!-

elemental function modulo2_qp (x, amp) result (mod2)


  use precision_def

  implicit none

  real(quadp), intent(in) :: x, amp
  real(quadp) mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_int (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Input:
!   x    -- Integer: 
!   amp  -- Integer: Must be positive.
!
! Output:
!   mod2 -- Integer: Result
!-

elemental function modulo2_int (x, amp) result (mod2)


  use precision_def

  implicit none

  integer, intent(in) :: x, amp
  integer mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function


end module
