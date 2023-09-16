!+
! Subroutine complex_error_function (wr, wi, zr, zi)
!
! This routine evaluates the function w(z) in the first quadrant of
! the complex plane. That is:
!    zr > 0 and zi > 0
! Three different expressions, pade1, pade2, and
! asymp, are used, depending on where z lies in the quadrant.
!
! Input:
!   zr -- Real(rp): Real part of z.
!   zi -- Real(rp): Imaginary part of z.
!
! Output:
!   wr -- Real(rp): Real part of w.
!   wi -- Real(rp): Imaginary part of w.
!-

subroutine complex_error_function (wr, wi, zr, zi)

  use precision_def

  implicit none

  real(rp) x1, x2, x3, x4, x5, y1, y2, r2, wr, wi, zr, zi, eps1, &
       yc, yc1, yc2

  parameter (x1 = 4.1,   x2 = 3.6, x3 = 3.5, x4 = 2.7, x5 = 2.2)
  parameter (y1 = 1.275, y2 = 1.095)
  parameter (r2 = 8.7025)

! data for pade1

  real(rp) c1, c2, c3, c4, c5, c6, &
       d1, d2, d3, d4, d5, d6, d7, &
       u2r, u2i, u3r, u3i, u4r, u4i, u5r, u5i, u6r, u6i, u7r, u7i, &
       fr, fi, dr, di, de

  parameter (c1 = -1.25647718d0,  c2 = 8.25059158d-1, &
             c3 = -3.19300157d-1, c4 = 7.63191605d-2, &
             c5 = -1.04697938d-2, c6 = 6.44878652d-4)

  parameter (d1 = -2.38485635d0,  d2 = 2.51608137d0, &
             d3 = -1.52579040d0,  d4 = 5.75922693d-1, &
             d5 = -1.35740709d-1, d6 = 1.85678083d-2, &
             d7 = -1.14243694d-3)


! data for pade2

  real(rp) c0r, c0i, c1r, c1i, c2r, c2i, c3r, c3i, &
       d1r, d1i, d2r, d2i, d3r, d3i, d4r, d4i, &
       z2r, z2i, z3r, z3i, z4r, z4i, &
       zr3

  parameter (c0r = 1.23409804d-4, c0i = 2.01157318d-1, &
             c1r = 2.33746715d-1, c1i = 1.61133338d-1, &
             c2r = 1.25689814d-1, c2i = -4.0422725d-2, &
             c3r = 8.92089179d-3, c3i = -1.81293213d-2)

  parameter (d1r = 1.19230984d0,   d1i = -1.16495901d0, &
             d2r = 8.9401545d-2,   d2i = -1.07372867d0, &
             d3r = -1.68547429d-1, d3i = -2.70096451d-1, &
             d4r = -3.20997564d-2, d4i = -1.58578639d-2)

! data for asymp

  real(rp) a1p, a2p, a3p, a4p, a5p, &
       b1, b2, b3, b4, b5, &
       dd1e, dde1, dd2e, dde2, dd3e, dde3, dd4e, dde4, dd5e, dde5, &
       dd1r, ddr1, dd2r, ddr2, dd3r, ddr3, dd4r, ddr4, dd5r, ddr5, &
       wi0, eps, pi2, &
       xx1, xx2

  parameter (a1p = 1.94443615d-1, a2p = 7.64384940d-2, &
             a3p = 1.07825546d-2, a4p = 4.27695730d-4, a5p = 2.43202531d-6)

  parameter (b1 = 3.42901327d-1, b2 = 1.036610830d0, b3 = 1.756683649d0, &
             b4 = 2.532731674d0, b5 = 3.436159119d0)

  parameter (pi2 = 1.12837917d0)

  parameter (xx1 = 3.5d0, xx2 = 4.2d0)

  parameter (eps = 0.01)

!

  eps1 = 0.0625d0 * (zr-x3)
!
  if (zr < x5) then
    if (zr*zr+zi*zi >= r2) then
      go to 30
    else
      go to 10
    endif
  endif
!
  if (zr <= x4) then
    yc1 = -1.40 * (zr-x4) + y1
    yc2 =  1.75 * (zr-x4) + y1
    if (zi >= yc1) go to 30
    if (zi >= yc2) then
      go to 10
    else
       go to 20
    endif
  endif
!
  if (zr < x2) then
    yc = -0.2 * (zr-x4) + y1
    if (zi >= yc) goto 30
    if (zr >= x3 .and. zr*zi < eps1) then
      go to 30
    else
      go to 20
    endif
  endif
!
  if (zr <= x1) then
    yc = -1.4 * (zr-x2) + y2
    if (zi >= yc) go to 30
    if (zr*zi < eps1) then
      go to 30
    else
      go to 20
    endif
  endif
!
  goto 30     ! since zr > x1

!+
! This Section of code formally called pade1(wr, wi, zr, zi)
!
! This program calculates a pade approximation of w(z)
! around the origin.
!-

10    continue

  u2r = zi * zi - zr * zr
  u2i =  - 2.d0 * zr * zi
  u3r =  - u2r * zi - u2i * zr
  u3i = u2r * zr - u2i * zi
  u4r =  - u3r * zi - u3i * zr
  u4i = u3r * zr - u3i * zi
  u5r =  - u4r * zi - u4i * zr
  u5i = u4r * zr - u4i * zi
  u6r =  - u5r * zi - u5i * zr
  u6i = u5r * zr - u5i * zi
  u7r =  - u6r * zi - u6i * zr
  u7i = u6r * zr - u6i * zi

  fr = 1.d0 - c1 * zi + c2 * u2r + c3 * u3r + c4 * u4r + c5 * u5r + c6 * u6r
  fi = c1 * zr + c2 * u2i + c3 * u3i + c4 * u4i + c5 * u5i + c6 * u6i

  dr = 1.d0 - d1 * zi + d2 * u2r + d3 * u3r + d4 * u4r + d5 * u5r + d6 * u6r + d7 * u7r
  di = d1 * zr + d2 * u2i + d3 * u3i + d4 * u4i + d5 * u5i + d6 * u6i + d7 * u7i
  de = dr * dr + di * di

  wr = (fr * dr + fi * di) / de
  wi = (fi * dr - fr * di) / de

  return


!+
! This section of code formally called pade2(wr, wi, zr, zi)
!
! This program calculates a pade approximation of w(z)
! around the point z = 3.
!-

20    continue

  zr3 = zr - 3.d0

  z2r = zr3 * zr3 - zi * zi
  z2i = 2.d0 * zr3 * zi
  z3r = z2r * zr3 - z2i * zi
  z3i = z2r * zi + z2i * zr3
  z4r = z3r * zr3 - z3i * zi
  z4i = z3r * zi + z3i * zr3

  fr = c0r + c1r * zr3 - c1i * zi + c2r * z2r - c2i * z2i + c3r * z3r - c3i * z3i
  fi = c0i + c1r * zi + c1i * zr3 + c2r * z2i + c2i * z2r + c3r * z3i + c3i * z3r

  dr = 1.d0 + d1r * zr3 - d1i * zi + d2r * z2r - d2i * z2i + d3r * z3r - d3i * z3i + &
        d4r * z4r - d4i * z4i
  di = d1r * zi + d1i * zr3 + d2r * z2i + d2i * z2r + d3r * z3i + d3i * z3r + d4r * z4i + &
        d4i * z4r
  de = dr * dr + di * di

  wr = (fr * dr + fi * di) / de
  wi = (fi * dr - fr * di) / de

  return


!+
! This section of code formally called asymp(wr, wi, zr, zi)
!
! This program calculates an asymptotic expression of
! w(z) valid away from the origin.
!-

30    continue

  ddr1 = zr + b1
  dd1r = zr - b1
  ddr2 = zr + b2
  dd2r = zr - b2
  ddr3 = zr + b3
  dd3r = zr - b3
  ddr4 = zr + b4
  dd4r = zr - b4
  ddr5 = zr + b5
  dd5r = zr - b5
  dde1 = ddr1**2 + zi**2
  dd1e = dd1r**2 + zi**2
  dde2 = ddr2**2 + zi**2
  dd2e = dd2r**2 + zi**2
  dde3 = ddr3**2 + zi**2
  dd3e = dd3r**2 + zi**2
  dde4 = ddr4**2 + zi**2
  dd4e = dd4r**2 + zi**2
  dde5 = ddr5**2 + zi**2
  dd5e = dd5r**2 + zi**2

  wi = a1p * (ddr1/dde1 + dd1r/dd1e) + a2p * (ddr2/dde2 + dd2r/dd2e) + &
        a3p * (ddr3/dde3 + dd3r/dd3e) + a4p * (ddr4/dde4 + dd4r/dd4e) + &
        a5p * (ddr5/dde5 + dd5r/dd5e)

  if (zr >= xx1) then
    eps1 = .04d0 / (zr - 3.29d0) - .034d0
    if ((zr*zi < eps1) .or. ((zr >= xx2) .and. (zr * zi < eps))) then
      dde1 = ddr1 * ddr1
      dd1e = dd1r * dd1r
      dde2 = ddr2 * ddr2
      dd2e = dd2r * dd2r
      dde3 = ddr3 * ddr3
      dd3e = dd3r * dd3r
      dde4 = ddr4 * ddr4
      dd4e = dd4r * dd4r
      dde5 = ddr5 * ddr5
      dd5e = dd5r * dd5r

      wi0 = a1p * (ddr1/dde1 + dd1r/dd1e) + a2p * (ddr2/dde2 + dd2r/dd2e) + &
           a3p * (ddr3/dde3 + dd3r/dd3e) + a4p * (ddr4/dde4 + dd4r/dd4e) + &
           a5p * (ddr5/dde5 + dd5r/dd5e)

      wr = exp(- zr * zr) + 2.d0 * wi0 * zr * zi - pi2 * zi

      return

    endif
  endif

  wr = (a1p * (1.d0/dde1 + 1.d0/dd1e) + a2p * (1.d0/dde2 + 1.d0/dd2e) + &
        a3p * (1.d0/dde3 + 1.d0/dd3e) + a4p * (1.d0/dde4 + 1.d0/dd4e) + &
        a5p * (1.d0/dde5 + 1.d0/dd5e)) * zi

end subroutine
