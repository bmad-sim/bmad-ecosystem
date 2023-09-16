!+
! Subroutine gen_grad_at_s_to_em_taylor (ele, gen_grad, s_pos, em_taylor)
!
! Routine to return the equivalent Taylor field map at a point s_pos.
!
! Input:
!   ele           -- ele: Element containing the map.
!   gen_grad      -- gen_grad_map_struct: Gen_grad map.
!   s_pos         -- real(rp): Position to evaluate em_taylor at.
!
! Output:
!   em_taylor(3)  -- em_taylor_struct: Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
!-

subroutine gen_grad_at_s_to_em_taylor (ele, gen_grad, s_pos, em_taylor)

use bmad_interface, dummy => gen_grad_at_s_to_em_taylor

implicit none

type em_taylor_coef_struct
  real(rp), allocatable :: c(:,:) ! (deriv-order, x_power)  Note: x_power + y_power = deriv_order
end type

type (ele_struct) ele
type (gen_grad_map_struct), target :: gen_grad
type (em_taylor_struct), target :: em_taylor(3)
type (gen_grad1_struct), pointer :: gg
type (em_taylor_coef_struct) em_coef(3)

real(rp) s_pos, s0, coef, scale, z_rel, s_here
real(rp), allocatable ::xy_plus(:), xy_zero(:), xy_minus(:)
real(rp), allocatable :: der(:)

integer iz0, nd
integer i, j, k, d, n, m, io, ix, m_max, iord, it, ig

logical is_even
character(*), parameter :: r_name = 'gen_grad_at_s_to_em_taylor'

! Find where to interpolate

select case (gen_grad%ele_anchor_pt)
case (anchor_beginning$); s0 = 0
case (anchor_center$);    s0 = ele%value(l$) / 2
case (anchor_end$);       s0 = ele%value(l$)
end select

s_here = s_pos - s0 - gen_grad%r0(3)

iz0 = floor(s_here / gen_grad%dz)
if (iz0 < gen_grad%iz0) iz0 = iz0 + 1 ! Allow one dz width out-of-bounds.
if (iz0 >= gen_grad%iz1) iz0 = iz0 - 1 ! Allow one dz width out-of-bounds.

if (iz0 < gen_grad%iz0 .or. iz0 >= gen_grad%iz1) then
  !if (.not. logic_option(.false., grid_allow_s_out_of_bounds)) then
    call out_io (s_error$, r_name, 'PARTICLE Z  \F10.3\ POSITION OUT OF BOUNDS.', &
                                   'FOR GEN_GRAD_MAP IN ELEMENT: ' // ele%name, r_array = [s_pos])
  !endif
  return
endif

!

z_rel = s_here - iz0 * gen_grad%dz

! Find largest order

io = 0
m_max = 0
do i = 1, size(gen_grad%gg)
  gg => gen_grad%gg(i)
  io = max(io, max(1, gg%m-1) + gg%n_deriv_max)
  m_max = max(m_max, gg%m)
enddo

allocate (em_coef(1)%c(0:io,0:io), em_coef(2)%c(0:io,0:io), em_coef(3)%c(0:io,0:io))  ! (poly order, x^n power)
allocate (xy_plus(0:m_max+1), xy_zero(0:m_max), xy_minus(0:m_max-1))

em_coef(1)%c = 0; em_coef(2)%c = 0; em_coef(3)%c = 0

!

do ig = 1, size(gen_grad%gg)
  gg => gen_grad%gg(ig)
  nd = gg%n_deriv_max
  m = gg%m
  is_even = .false.

  do j = 0, m/2
    xy_plus(m-2*j) = (-1)**j * n_choose_k(m+1, 2*j+1)  ! S_xy(m+1) coefs
  enddo

  do j = 0, (m+1)/2
    xy_plus(m-2*j+1) = (-1)**j * n_choose_k(m+1, 2*j)  ! C_xy(m+1) coefs
  enddo

  do j = 0, divide_by_two(m-1)
    xy_zero(m-2*j-1) = (-1)**j * n_choose_k(m, 2*j+1)  ! S_xy(m) coefs
  enddo

  do j = 0, m/2
    xy_zero(m-2*j) = (-1)**j * n_choose_k(m, 2*j)  ! C_xy(m) coefs
  enddo

  do j = 0, divide_by_two(m-2)
    xy_minus(m-2*j-2) = (-1)**j * n_choose_k(m-1, 2*j+1)  ! S_xy(m-1) coefs
  enddo

  do j = 0, divide_by_two(m-1)
    xy_minus(m-2*j-1) = (-1)**j * n_choose_k(m-1, 2*j)  ! C_xy(m-1) coefs
  enddo

  call re_allocate2(der, 0, nd, .false.)

  do d = 0, nd
    is_even = (.not. is_even)
    coef = poly_eval(gg%deriv(iz0, d:), z_rel, diff_coef=.true.)
    if (coef == 0) cycle

    if (is_even) then
      n = d / 2
      coef = coef * (-1)**n * factorial(m) / (4**n * factorial(n) * factorial(n+m))
      iord = m + d - 1

      if (gg%sincos == sin$) then
        if (m == 0) cycle

        do k = 0, n
          do j = 0, divide_by_two(m-2)
            it = m-2*j-2
            em_coef(1)%c(iord,2*k+it) = em_coef(1)%c(iord,2*k+it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(it)  ! Sin
          enddo
          do j = 0, divide_by_two(m-1)
            it = m-2*j-1
            em_coef(2)%c(iord,2*k+it) = em_coef(2)%c(iord,2*k+it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(it)  ! Cos
          enddo
        enddo

        do k = 0, n-1
          do j = 0, m/2
            it = m-2*j
            em_coef(1)%c(iord,2*k+it) = em_coef(1)%c(iord,2*k+it) + coef * n * n_choose_k(n-1,k) * xy_plus(it)  ! Sin
          enddo
          do j = 0, (m+1)/2
            it = m-2*j+1
            em_coef(2)%c(iord,2*k+it) = em_coef(2)%c(iord,2*k+it) - coef * n * n_choose_k(n-1,k) * xy_plus(it)  ! Cos
          enddo
        enddo

      else  ! sincos = cos$
        if (m == 0) then
          if (d == 0) cycle

          ! For m = 0 use the mapping: (sin(m-1), cos(m-1)) -> (-sin(1-m), cos(1-m))
          do k = 0, n
            it = 1
            em_coef(1)%c(iord,2*k+it) = em_coef(1)%c(iord,2*k+it) + coef * (n) * n_choose_k(n,k)  ! Cos
            it = 0
            em_coef(2)%c(iord,2*k+it) = em_coef(2)%c(iord,2*k+it) + coef * (n) * n_choose_k(n,k)  ! Sin
          enddo

        else  ! m /= 0
          do k = 0, n
            do j = 0, divide_by_two(m-1)
              it = m-2*j-1
              em_coef(1)%c(iord,2*k+it) = em_coef(1)%c(iord,2*k+it) + coef * (n+m) * n_choose_k(n,k) * xy_minus(it)  ! Cos
            enddo
            do j = 0, divide_by_two(m-2)
              it = m-2*j-2
              em_coef(2)%c(iord,2*k+it) = em_coef(2)%c(iord,2*k+it) - coef * (n+m) * n_choose_k(n,k) * xy_minus(it)  ! Sin
            enddo
          enddo
        endif

        do k = 0, n-1
          do j = 0, (m+1)/2
            it = m-2*j+1
            em_coef(1)%c(iord,2*k+it) = em_coef(1)%c(iord,2*k+it) + coef * n * n_choose_k(n-1,k) * xy_plus(it)  ! Cos
          enddo
          do j = 0, m/2
            it = m-2*j
            em_coef(2)%c(iord,2*k+it) = em_coef(2)%c(iord,2*k+it) + coef * n * n_choose_k(n-1,k) * xy_plus(it)  ! Sin
          enddo
        enddo
      endif

    else  ! is odd
      n = (d - 1) / 2
      coef = coef * (-1)**n * factorial(m) / (4**n * factorial(n) * factorial(n+m))
      iord = m + d - 1
      if (gg%sincos == sin$) then
        do k = 0, n
          do j = 0, divide_by_two(m-1)
            it = m-2*j-1
            em_coef(3)%c(iord,2*k+it) = em_coef(3)%c(iord,2*k+it) + coef * n_choose_k(n,k) * xy_zero(it)  ! Sin
          enddo
        enddo

      else
        do k = 0, n
          do j = 0, m/2
            it = m-2*j
            em_coef(3)%c(iord,2*k+it) = em_coef(3)%c(iord,2*k+it) + coef * n_choose_k(n,k) * xy_zero(it)  ! Cos
          enddo
        enddo
      endif

    endif  ! is_even
  enddo  ! deriv index
enddo  ! gg

!

scale = gen_grad%field_scale * master_parameter_value(gen_grad%master_parameter, ele)

do i = 1, 3
  n = count(em_coef(i)%c /= 0)
  call init_em_taylor_series(em_taylor(i), n)

  n = 0
  do io = 0, ubound(em_coef(i)%c,1)
    do ix = 0, io
      if (em_coef(i)%c(io,ix) == 0) cycle
      n = n + 1
      em_taylor(i)%term(n)%coef = em_coef(i)%c(io,ix) * scale
      em_taylor(i)%term(n)%expn = [ix, io-ix]
    enddo
  enddo
enddo

!-------------------------------------------------------------
contains

! This routine is needed since "(-1)/2" evaluate to 0 but what is wanted is -1.

function divide_by_two(mm) result (m2)
integer mm, m2
if (mm < 0) then
  m2 = -1
else
  m2 = mm/2
endif
end function divide_by_two

end subroutine
