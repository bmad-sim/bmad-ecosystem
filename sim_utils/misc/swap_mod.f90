module swap_mod

use precision_def

implicit none

interface swap
  module procedure swap_i
  module procedure swap_r
  module procedure swap_rv
  module procedure swap_c
  module procedure swap_cv
  module procedure swap_cm
end interface

contains

subroutine swap_i(a,b)
integer, intent(inout) :: a,b
integer :: dum
dum=a
a=b
b=dum
end subroutine swap_i


subroutine swap_r(a,b)
real(rp), intent(inout) :: a,b
real(rp) :: dum
dum=a
a=b
b=dum
end subroutine swap_r


subroutine swap_rv(a,b)
real(rp), intent(inout) :: a(:),b(:)
real(rp), dimension(size(a)) :: dum
dum=a
a=b
b=dum
end subroutine swap_rv


subroutine swap_c(a,b)
complex(rp), intent(inout) :: a,b
complex(rp) :: dum
dum=a
a=b
b=dum
end subroutine swap_c


subroutine swap_cv(a,b)
complex(rp), intent(inout) :: a(:),b(:)
complex(rp), dimension(size(a)) :: dum
dum=a
a=b
b=dum
end subroutine swap_cv


subroutine swap_cm(a,b)
complex(rp), intent(inout) :: a(:,:),b(:,:)
complex(rp), dimension(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
end subroutine swap_cm

end module
