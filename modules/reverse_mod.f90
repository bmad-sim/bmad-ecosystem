#include "CESR_platform.inc"

module reverse_mod

  use bmad_struct
  use bmad_interface
  use ptc_interface_mod, only: taylor_inverse

contains

!+
! Subroutine ring_reverse (ring_in, ring_rev)
!
! Subroutine to construct a ring structure with the elements in reversed order.
! This may be used for backward tracking through the ring. 
!
! The correspondence between elements in the two rings is as follows:
!     ring_rev%ele_(ring%n_ele_ring+1-i) = ring_in%ele_(i)  
!                                                for 0 < i <= ring%n_ele_ring
!     ring_rev%ele_(i)                   = ring_in%ele_(i)   
!                                                for ring%n_ele_ring < i 
!
! All longitudial quantities (for example, the value of ks for a solenoid) 
! are flipped in sign for the reversed ring. 
! This means that the appropriate transformation from particle coordinates
! in one ring the corrisponding coordinates in the other is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Note: The Twiss parameters will not be correct for the reversed ring.
! You will need to compute them.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring_in -- Ring_struct: Input ring.
!
! Output:
!   ring_rev -- Ring_struct: Ring with the elements in reversed order.
!-

subroutine ring_reverse (ring_in, ring_rev)

  implicit none

  type (ring_struct), intent(in) :: ring_in
  type (ring_struct), intent(out), target :: ring_rev
  type (ele_struct), pointer :: ele
  type (control_struct), pointer :: con

  integer i, n, i1, i2, nr, n_con
  integer :: ix_con(size(ring_in%control_))

! transfer

  n_con = size(ring_in%control_)

  ring_rev = ring_in

  nr = ring_rev%n_ele_ring
  ring_rev%ele_(1:nr) = ring_rev%ele_(nr:1:-1)

! flip longitudinal stuff, maps

  do i = 1, ring_rev%n_ele_max
    call reverse_ele (ring_rev%ele_(i))
  enddo

! correct control information

  do i = 1, ring_rev%n_control_max
    con => ring_rev%control_(i)
    if (con%ix_slave <= nr) con%ix_slave = nr+1-con%ix_slave
    if (con%ix_lord <= nr)  con%ix_lord  = nr+1-con%ix_lord
  enddo

! slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.


  forall (i = 1:n_con) ix_con(i) = i 

  do i = nr+1, ring_rev%n_ele_max
    ele => ring_rev%ele_(i)
    if (ele%control_type /= super_lord$) cycle
    i1 = ele%ix1_slave
    i2 = ele%ix2_slave
    ring_rev%control_(i1:i2) = ring_rev%control_(i2:i1:-1)
    ix_con(i1:i2) = ix_con(i2:i1:-1)
    ele%s = ring_rev%param%total_length - ele%s + ele%value(l$)
  enddo

  n = ring_rev%n_ic_max
  ring_rev%ic_(1:n) = ix_con(ring_rev%ic_(1:n))

! Cleanup

  call s_calc (ring_rev) 
  call check_ring_controls (ring_rev, .true.)

end subroutine

!--------------------------------------------------------------------------
!+
! Subroutine reverse_ele (ele)
!
! Subroutine to "reverse" an element for backward tracking.
!
! All longitudial quantities (for example, the value of ks for a solenoid) 
! are flipped in sign for the reversed ring. 
! This means that the appropriate transformation that corresponds to the 
! reverse transformation is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Input element.
!
! Output:
!   ele -- Ele_struct: Reversed element.
!-

subroutine reverse_ele (ele)

  use ptc_interface_mod, only: kill

  implicit none

  type (ele_struct) ele

  integer i, j, sum245

  real(rp) tempp

! Flip longitudinal attributes

  ele%value(x_pitch$) = -ele%value(x_pitch$)
  ele%value(y_pitch$) = -ele%value(y_pitch$)

  select case (ele%key)

  case (solenoid$, sol_quad$)
    ele%value(ks$) = -ele%value(ks$)

  case (rfcavity$)
    ele%value(phi0$) = -ele%value(phi0$)

  case (sbend$)
    tempp = ele%value(e1$)
    ele%value(e1$) = ele%value(e2$)
    ele%value(e2$) = tempp

! For wigglers:
!       phi -> -phi -  k_z * Length
! This transforms:
!       (B_x, B_y, B_z) -> (B_x, B_y, -B_z)

  case (wiggler$)
    if (associated(ele%wig_term)) then
      do i = 1, size(ele%wig_term)
        ele%wig_term(i)%phi_z = -ele%wig_term(i)%phi_z - &
                                        ele%wig_term(i)%kz * ele%value(l$)
      enddo
    endif

  end select

! The inverse of the Taylor map is: 
!   M * T^(-1) * M
! where: 
!   T = the input taylor series
!   M = The map (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)

  if (associated(ele%taylor(1)%term)) then

    call taylor_inverse (ele%taylor, ele%taylor)

! Apply M to the right

    do i = 1, 6
      do j = 1, size(ele%taylor(i)%term)
        sum245 = ele%taylor(i)%term(j)%exp(2) + ele%taylor(i)%term(j)%exp(4) + &
                                                  ele%taylor(i)%term(j)%exp(5)
        if (mod(sum245, 2) == 1) ele%taylor(i)%term(j)%coef = &
                                                    -ele%taylor(i)%term(j)%coef
      end do
    end do

! Apply M to the left

    ele%taylor(2)%term(:)%coef = -ele%taylor(2)%term(:)%coef
    ele%taylor(4)%term(:)%coef = -ele%taylor(4)%term(:)%coef
    ele%taylor(5)%term(:)%coef = -ele%taylor(5)%term(:)%coef

  endif

! kill any gen_field

  if (associated(ele%gen_field)) call kill (ele%gen_field)

! reverse mat6

  call mat_symp_conj (ele%mat6, ele%mat6, 6, 6)
  ele%vec0 = -matmul(ele%mat6, ele%vec0)

  ele%mat6(2,:) = -ele%mat6(2,:)
  ele%mat6(4,:) = -ele%mat6(4,:)
  ele%mat6(5,:) = -ele%mat6(5,:)

  ele%vec0(2) = -ele%vec0(2)
  ele%vec0(4) = -ele%vec0(4)
  ele%vec0(5) = -ele%vec0(5)

  ele%mat6(:,2) = -ele%mat6(:,2)
  ele%mat6(:,4) = -ele%mat6(:,4)
  ele%mat6(:,5) = -ele%mat6(:,5)

end subroutine

end module


