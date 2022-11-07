!+
! Subroutine gen_grad_to_em_taylor (gen_grad, iz, em_taylor)
!
! Routine to take a gen_grad coefs at a given z-position and return the equivalent Taylo field map.
!
! Input:
!   gen_grad      -- gen_grad_map_struct: Gen_grad map.
!   iz            -- integer: z-plane index to evaluate.
!   em_taylor(3)  -- em_taylor_struct: Map for (Bx, By, Bz) or (Ex, Ey, Ez) field.
!-

subroutine gen_grad_to_em_taylor (gen_grad, iz, em_taylor)

use bmad_interface, dummy => gen_grad_to_em_taylor

implicit none

type (gen_grad_map_struct), target :: gen_grad
type (em_taylor_struct), target :: em_taylor(3)
type (gen_grad1_struct), pointer :: gg

real(rp) coef
real(rp), allocatable :: coef_x(:,:), coef_y(:,:), coef_z(:,:)
real(rp), allocatable ::xy_plus(:), xy_zero(:), xy_minus(:)

integer iz
integer i, j, k, d, n, m, io, ix, m_max, iord, it

logical is_even

! Find largest order

io = 0
m_max = 0
do i = 1, size(gen_grad%gg)
  gg => gen_grad%gg(i)
  io = max(io, gg%m + ubound(gg%deriv,2) - 1)
  m_max = max(m_max, gg%m)
enddo

allocate (coef_x(0:io,0:io), coef_y(0:io,0:io), coef_z(0:io,0:io))  ! (poly order, x^n power)
allocate (xy_plus(m_max+1), xy_zero(m_max), xy_minus(m_max-1))

coef_x = 0; coef_y = 0; coef_z = 0

!

do i = 1, size(gen_grad%gg)
  gg => gen_grad%gg(i)
  m = gg%m

  do j = 0, m/2
    xy_plus(m-2*j) = (-1)**j * n_choose_k(m+1, 2*j+1)
  enddo

  do j = 0, (m+1)/2
    xy_plus(m-2*j+1) = (-1)**j * n_choose_k(m+1, 2*j)
  enddo

  do j = 0, (m-1)/2
    xy_zero(m-2*j-1) = (-1)**j * n_choose_k(m, 2*j+1)
  enddo

  do j = 0, m/2
    xy_zero(m-2*j) = (-1)**j * n_choose_k(m, 2*j)
  enddo

  do j = 0, (m-2)/2
    xy_minus(m-2*j-2) = (-1)**j * n_choose_k(m-1, 2*j+1)
  enddo

  do j = 0, (m-1)/2
    xy_minus(m-2*j-1) = (-1)**j * n_choose_k(m-1, 2*j)
  enddo

  is_even = .false.
  do d = 0, ubound(gg%deriv,2)
    is_even = (.not. is_even)
    coef = gg%deriv(iz,n) * gen_grad%field_scale
    if (coef == 0) cycle

    if (is_even) then
      n = d / 2
      coef = coef * (-1)**n * factorial(m) / (4**n * factorial(n) * factorial(n+m))
      iord = m + d - 1
      if (gg%sincos == sin$) then
        do k = 0, n
          do j = 0, (m-2)/2
            it = 2*k+m-2*j-2
            coef_x(iord,it) = coef_x(iord,it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(m-2*j-2)  ! Sin
          enddo
          do j = 0, (m-1)/2
            it = 2*k+m-2*j-1
            coef_y(iord,it) = coef_y(iord,it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(m-2*j-1)  ! Cos
          enddo
        enddo

        do k = 0, n-1
          do j = 0, m/2
            it = 2*k+m-2*j
            coef_x(iord,it) = coef_x(iord,it) + coef * n * n_choose_k(n-1,k) * xy_plus(m-2*j)  ! Sin
          enddo
          do j = 0, (m-1)/2
            it = 2*k+m-2*j-1
            coef_y(iord,it) = coef_y(iord,it) - coef * n * n_choose_k(n-1,k) * xy_plus(m-2*j+1)  ! Cos
          enddo
        enddo

      else  ! sincos = cos$
        do k = 0, n
          do j = 0, (m-1)/2
            it = 2*k+m-2*j-1
            coef_x(iord,it) = coef_x(iord,it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(m-2*j-1)  ! Cos
          enddo
          do j = 0, (m-2)/2
            it = 2*k+m-2*j-2
            coef_y(iord,it) = coef_y(iord,it) - coef * (n+m) * n_choose_k(n,k) * xy_minus(m-2*j-2)  ! Sin
          enddo
        enddo

        do k = 0, n-1
          do j = 0, (m-1)/2
            it = 2*k+m-2*j-1
            coef_x(iord,it) = coef_x(iord,it) - coef * n * n_choose_k(n-1,k) * xy_plus(m-2*j+1)  ! Cos
          enddo
          do j = 0, m/2
            it = 2*k+m-2*j
            coef_y(iord,it) = coef_y(iord,it) + coef * n * n_choose_k(n-1,k) * xy_plus(m-2*j)  ! Sin
          enddo
        enddo
      endif

    else  ! is odd
      n = (d - 1) / 2
      coef = coef * (-1)**n * factorial(m) / (4**n * factorial(n) * factorial(n+m))
      iord = m + d - 1
      if (gg%sincos == sin$) then
        do k = 0, n
          do j = 0, (m-1)/2
            it = 2*k+m-2*j-1
            coef_y(iord,it) = coef_y(iord,it) + coef * n_choose_k(n,k) * xy_zero(m-2*j-1)  ! Sin
          enddo
        enddo

      else
        do k = 0, n
          do j = 0, (m-1)/2
            it = 2*k+m-2*j
            coef_y(iord,it) = coef_y(iord,it) + coef * n_choose_k(n,k) * xy_zero(m-2*j)  ! Cos
          enddo
        enddo
      endif

    endif  ! is_even
  enddo  ! deriv index
enddo  ! gg

!

n = count(coef_x /= 0)
call init_em_taylor_series(em_taylor(1), n)

n = 0
do io = 0, ubound(coef_x,1)
  do ix = 0, io
    if (coef_x(io,ix) == 0) cycle
    n = n + 1
    em_taylor(1)%term(n)%coef = coef_x(io,ix)
    !em_taylor(1)%expn = 
  enddo
enddo

end subroutine
