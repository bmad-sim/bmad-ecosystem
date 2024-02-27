!+
! Function spin_taylor_to_linear (spin_taylor, normalize, dref_orb, is_on) result (spin_map1)
!
! Routine to truncate a Taylor spin map to order 1.
!
! Linear map normalization involves adjusting the 0th order quaternion to have unit magnitude
! and the first order quaternions to be perpendicular to the 0th order quaternion.
!
! Input:
!   spin_taylor(0:3)    -- taylor_struct: Taylor spin map.
!   normalize           -- logical: If True, normalize the linear map.
!   dref_orb(6)         -- real(rp): Change in Reference orbit: output_map1_ref - input_taylor_ref.
!   is_on               -- logical: Is map turned on? If not spin_map1 will be the unit map.
!
! Output:
!   spin_map1(0:3,0:6)  -- real(rp): First order spin map.
!-

function spin_taylor_to_linear (spin_taylor, normalize, dref_orb, is_on) result (spin_map1)

use bmad_routine_interface, dummy => spin_taylor_to_linear

implicit none

type (taylor_struct), target :: spin_taylor(0:3)
type (taylor_struct), pointer :: st
type (taylor_term_struct), pointer :: term

real(rp) dref_orb(6), spin_map1(0:3,0:6)
real(rp) t(6), t_out, prod, f_norm

integer i, j, k, l, n, p
logical normalize, is_on

!

spin_map1 = 0

if (.not. is_on) then
  spin_map1(0,0) = 1.0_rp
  return
endif

!

if (all(dref_orb == 0)) then
  do i = 0, 3
    st => spin_taylor(i)

    do k = 1, size(st%term)
      term => st%term(k)
      n = sum(st%term(k)%expn)

      select case (n)
      case (0)
        spin_map1(i,0) = st%term(k)%coef
      case (1)
        do p = 1, 6
          if (st%term(k)%expn(p) == 0) cycle
          spin_map1(i,p) = st%term(k)%coef
          exit
        enddo
      end select
    enddo
  enddo

!

else
  do i = 0, 3
    st => spin_taylor(i)

    do k = 1, size(st%term)
      term => st%term(k)

      t_out = term%coef
      do l = 1, 6
        if (term%expn(l) == 0) cycle
        t(l) = dref_orb(l)**term%expn(l)
        t_out = t_out * t(l)
      enddo
      spin_map1(i,0) = spin_map1(i,0) + t_out

      do j = 1, 6
        if (term%expn(j) == 0) cycle
        if (term%expn(j) > 1 .and. dref_orb(j) == 0) cycle

        if (term%expn(j) > 1)then
          prod = term%coef * term%expn(j) * dref_orb(j)**(term%expn(j)-1)
        else  ! term%expn(j) == 1
          prod = term%coef
        endif

        do l = 1, 6
          if (term%expn(l) == 0) cycle
          if (l == j) cycle
          prod = prod * t(l)
        enddo

        spin_map1(i,j) = spin_map1(i,j) + prod
      enddo
    enddo
  enddo
endif

! Normalize

if (.not. normalize) return

! Renormalize to make 0th order qaternion have unit length
f_norm = 1.0_rp / norm2(spin_map1(:,0))
spin_map1(:,:) = spin_map1(:,:)  * f_norm

! Now make first order quaternions perpendicular to the 0th order quaternions.
do j = 1, 6
  f_norm = dot_product(spin_map1(:,0), spin_map1(:,j))
  spin_map1(:,j) = spin_map1(:,j) - f_norm * spin_map1(:,0)
enddo

end function
