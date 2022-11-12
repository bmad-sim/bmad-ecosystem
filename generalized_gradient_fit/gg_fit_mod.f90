module gg_fit_mod

use lmdif_mod
use bmad_struct
use expression_mod, only: sin$, cos$
use super_recipes_mod

implicit none

!

type gg1_struct
  integer :: sincos = -1    ! sin$ or cos$
  integer :: m = -1         ! Azimuthal index
  integer :: sym_score = 1  ! Symmetry score. 0 => do not use in fit.
  real(rp), allocatable :: deriv(:,:)    ! (iz, deriv)
  integer, allocatable :: ix_var(:)      ! Index to var_vec(:) array.
end type

type (gg1_struct), allocatable, target :: gg1(:)

type fit_info_struct
  real(rp) :: rms0 = 0      ! Initial field rms
  real(rp) :: rms_fit = 0   ! Field_fit - Field_data rms
end type

type (fit_info_struct), allocatable, target :: fit(:)

integer, parameter :: m_max = 10
integer :: Nx_min, Nx_max, Ny_min, Ny_max, Nz_min, Nz_max
integer :: n_grid_pts
integer :: n_cycles = 100000
integer :: every_n_th_plane = 1, n_planes_add = 0, n_deriv_keep = -1
integer :: n_deriv_max
integer :: iz_min = int_garbage$, iz_max = int_garbage$    ! Compute range
integer :: m_cos(m_max) = -1
integer :: m_sin(m_max) = -1
integer :: sym_x = 0, sym_y = 0
integer n_var0, n_var, n_merit

real(rp), allocatable :: Bx_dat(:,:,:), By_dat(:,:,:), Bz_dat(:,:,:)
real(rp), allocatable :: Bx_fit(:,:), By_fit(:,:), Bz_fit(:,:)
real(rp), allocatable :: Bx_diff(:,:), By_diff(:,:), Bz_diff(:,:)
real(rp), allocatable :: dBx_dvar(:,:), dBy_dvar(:,:), dBz_dvar(:,:)

! del_grid is the difference between (x,y,z) points in the units of the grid file.
! length_scale is table file (x,y,z) units in meters. length_scale = table_units/meters
! field_scale is the field table field in Tesla or V/m.

real(rp) :: del_grid(3), r0_grid(3), field_scale, length_scale, r_max
real(rp) :: z_min = real_garbage$, z_max = real_garbage$
real(rp) :: lmdif_eps = 1d-12

real(rp), allocatable :: var_vec(:), merit_vec(:), dB_dvar_vec(:,:)
real(rp) merit

logical printit

character(10) :: optimizer = 'lmdif'
character(100) :: field_file, out_file = ''

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine read_field_table (field_file)
! 
! Routine to read in a field table. Units are: cm, Gauss.
! Look at the code for the format.
!
! Input:
!   field_file    -- character(*): Name of file.
!-

subroutine read_field_table (field_file)

implicit none

real(rp) xx, yy, zz, Bx, By, Bz

integer ios, ios1, ios2
integer :: i, j, k, ix, iy, iz

logical, allocatable :: valid_field(:,:,:)

character(*) :: field_file
character(200) line

! Binary
! Notice that the binary table stores lengths in meters and fields in Tesla independent of length_scale and field_scale.

if (index(field_file, '.binary') /= 0) then

  open (1, file = field_file, form = 'unformatted', status = 'OLD', action = 'READ')
  read (1) length_scale, field_scale
  read (1) Nx_min, Nx_max
  read (1) Ny_min, Ny_max
  read (1) Nz_min, Nz_max
  read (1) del_grid
  read (1) r0_grid

  allocate (Bx_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), By_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))
  allocate (Bz_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), valid_field(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))

  valid_field = .false.
  do
    read (1, iostat = ios) i, j, k, Bx_dat(i,j,k), By_dat(i,j,k), Bz_dat(i,j,k)
    if (ios /= 0) exit
    valid_field(i,j,k) = .true.
  enddo

! ASCII

else
  open (1, file = field_file, status = 'OLD', action = 'READ')

  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) length_scale
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('length_scale', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) field_scale
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('field_scale', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) del_grid
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('del_grid', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) r0_grid
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('r0_grid', line)

  Nx_min = 100000; Nx_max = -100000
  Ny_min = 100000; Ny_max = -100000
  Nz_min = 100000; Nz_max = -100000

  do
    read (1, '(a)', iostat = ios1) line
    if (line(1:1) == '!' .or. line == '') cycle
    read (line, *, iostat = ios2) xx, yy, zz, Bx, By, Bz
    if (ios1 < 0) exit
    if (ios1 /= 0 .or. ios2 /= 0) call read_err('xx, yy, zz, Bx, By, Bz', line)
    Nx_min = min(Nx_min, nint (xx/del_grid(1))); Nx_max = max(Nx_max, nint (xx/del_grid(1)))
    Ny_min = min(Ny_min, nint (yy/del_grid(2))); Ny_max = max(Ny_max, nint (yy/del_grid(2)))
    Nz_min = min(Nz_min, nint (zz/del_grid(3))); Nz_max = max(Nz_max, nint (zz/del_grid(3)))
  enddo

  allocate (Bx_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), By_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))
  allocate (Bz_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), valid_field(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))

  rewind (1)
  read (1, '(a)', iostat = ios1) line
  read (1, '(a)', iostat = ios1) line
  read (1, '(a)', iostat = ios1) line
  read (1, '(a)', iostat = ios1) line

  do
    read (1, '(a)', iostat = ios1) line
    if (line(1:1) == '!' .or. line == '') cycle
    read (line, *, iostat = ios2) xx, yy, zz, Bx, By, Bz
    if (ios1 < 0) exit
    if (ios1 /= 0 .or. ios2 /= 0) call read_err('xx, yy, zz, Bx, By, Bz', line)
    i = nint (xx/del_grid(1))
    j = nint (yy/del_grid(2))
    k = nint (zz/del_grid(3))

    Bx_dat(i,j,k) = Bx * field_scale
    By_dat(i,j,k) = By * field_scale
    Bz_dat(i,j,k) = Bz * field_scale
    valid_field(i,j,k) = .true.
  end do

  del_grid = del_grid * length_scale
  r0_grid = r0_grid * length_scale
endif

!

close (1)

print '(a, f12.4)', 'Length scale:', length_scale
print '(a, es11.3)', 'Field scale:', field_scale
print '(a, 2(f10.5, a, i5, a, 4x))', 'x table range (meters, index):', Nx_min*del_grid(1), ' (', Nx_min, ')', Nx_max*del_grid(1), ' (', Nx_max, ')'
print '(a, 2(f10.5, a, i5, a, 4x))', 'y table range (meters, index)):', Ny_min*del_grid(2), ' (', Ny_min, ')', Ny_max*del_grid(2), ' (', Ny_max, ')'
print '(a, 2(f10.5, a, i5, a, 4x))', 'z table range (meters, index)):', Nz_min*del_grid(3), ' (', Nz_min, ')', Nz_max*del_grid(3), ' (', Nz_max, ')'

n_grid_pts = (Nx_max-Nx_min+1) * (Ny_max-Ny_min+1) * (Nz_max-Nz_min+1)

print *, 'Field table read: ', trim(field_file)
print *, 'Number of grid field points:       ', n_grid_pts

if (n_grid_pts /= count(valid_field)) then
  print *, 'Number points missing in table: ', n_grid_pts - count(valid_field)
  print *, 'Will stop here!'
  do i = Nx_min, Nx_max
  do j = Ny_min, Ny_max
  do k = Nz_min, Nz_max
    if (valid_field(i,j,k)) cycle
    print '(a, 3i5, 3es12.4)', 'No value at:', i, j, k, [i,j,k]*del_grid/length_scale
  enddo
  enddo
  enddo
  stop
endif

!

if (z_min == real_garbage$) then
  iz_min = Nz_min
else
  iz_min = max(Nz_min, nint(z_min/del_grid(3)))
endif

if (z_max == real_garbage$) then
  iz_max = Nz_max
else
  iz_max = min(Nz_max, nint(z_max/del_grid(3)))
endif

r_max = norm2([max(abs(Nx_min), abs(Nx_max)) * del_grid(1), max(abs(Ny_min), abs(Ny_max)) * del_grid(2)])
print *, 'r_max radius: ', r_max

!--------------------------------------
contains

subroutine read_err(str, line)

character(*) str, line

print *, 'ERROR READING FIELD TABLE FILE: ', trim(field_file)
print *, 'ERROR OBTAINING VALUE FOR: ', trim(str)
print *, 'ERROR PARSING: ', trim(line)
stop

end subroutine read_err

end subroutine read_field_table

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine write_binary_field_table (field_file)
! 
! Routine to write in a field table. Units are: Tesla (or V/m) independent of field_scale value.
!
! Input:
!   field_file    -- character(*): Name of file
!-

subroutine write_binary_field_table (field_file)

integer ix, iy, iz
character(*) field_file

! Write binary table

open (1, file = trim(field_file) // '.binary', form = 'unformatted')

write (1) length_scale, field_scale
write (1) Nx_min, Nx_max
write (1) Ny_min, Ny_max
write (1) Nz_min, Nz_max
write (1) del_grid
write (1) r0_grid

do ix = Nx_min, Nx_max
do iy = Ny_min, Ny_max
do iz = Nz_min, Nz_max
  write (1) ix, iy, iz, Bx_dat(ix,iy,iz), By_dat(ix,iy,iz), Bz_dat(ix,iy,iz)
enddo
enddo
enddo

close (1)

end subroutine write_binary_field_table

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine write_ascii_field_table (field_file)
! 
! Routine to write a subset in z of the field table. 
!
! Input:
!   field_file    -- character(*): Name of file.
!-

subroutine write_ascii_field_table (field_file)

integer ix, iy, iz, n
character(*) field_file
character(40) fmt

! Write ascii table

ix = index(field_file, '.binary')
if (ix == 0) then
  open (1, file = trim(field_file) // '.partial')
else
  open (1, file = field_file(:ix-1) // '.partial')
endif

n = nint(log10(length_scale))
write (fmt, '(a,i0,a)') '(3f12.', 6+n, ', 3es20.11)'

write (1, '(es14.6, a)') length_scale, '  ! length_scale in meters'
write (1, '(es14.6, a)') field_scale,  '  ! field_scale in Tesla or V/m'
write (1, '(2i8, a)') Nx_min, Nx_max,  '  ! Nx_min, Nx_max'
write (1, '(2i8, a)') Ny_min, Ny_max,  '  ! Ny_min, Ny_max'
write (1, '(2i8, a)') iz_min, iz_max,  '  ! Nz_min, Nz_max'
write (1, '(3f10.6, a)') del_grid/length_scale, '  ! del_grid'
write (1, '(3f10.6, a)') r0_grid/length_scale,  '  ! r0_grid'


do ix = Nx_min, Nx_max
do iy = Ny_min, Ny_max
do iz = iz_min, iz_max
  write (1, fmt) ix*del_grid(1)/length_scale, iy*del_grid(2)/length_scale, iz*del_grid(3)/length_scale, &
              Bx_dat(ix,iy,iz)/field_scale, By_dat(ix,iy,iz)/field_scale, Bz_dat(ix,iy,iz)/field_scale
enddo
enddo
enddo

close (1)

end subroutine write_ascii_field_table

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

subroutine fit_field()

type (super_mrqmin_storage_struct) storage
type (gg1_struct), pointer :: gg

real(rp), allocatable :: vec0(:), weight(:), merit0_vec(:)
real(rp) merit0, x, y, v, chisq, a_lambda

integer iloop, im, id, n_gg, n, ig
integer ix, iy, iz, status, iz0, iz1

logical at_end

!

allocate (Bx_diff(Nx_min:Nx_max, Ny_min:Ny_max), By_diff(Nx_min:Nx_max, Ny_min:Ny_max), &
                                                            Bz_diff(Nx_min:Nx_max, Ny_min:Ny_max))
allocate (Bx_fit(Nx_min:Nx_max, Ny_min:Ny_max), By_fit(Nx_min:Nx_max, Ny_min:Ny_max), &
                                                            Bz_fit(Nx_min:Nx_max, Ny_min:Ny_max))

!

n_gg = (count(m_cos /= -1) + count(m_sin /= -1))
n_var0 = n_gg * (n_deriv_max + 1)
do im = 1, m_max
  if (m_cos(im) == 0) n_var0 = n_var0 - 1
  if (m_sin(im) == 0) then
    print *, 'M = 0 sin coefficients do not contribute to the field. Please remove m_sin set.'
    stop
  endif
enddo

n_merit = 3 * (Nx_max-Nx_min+1) * (Ny_max-Ny_min+1) * (2*n_planes_add+1)
allocate (merit_vec(n_merit), merit0_vec(n_merit), dB_dvar_vec(n_merit, n_var0))
allocate (dBx_dvar(Nx_min:Nx_max, Ny_min:Ny_max), dBy_dvar(Nx_min:Nx_max, Ny_min:Ny_max))
allocate (dBz_dvar(Nx_min:Nx_max, Ny_min:Ny_max))
allocate (gg1(n_gg))

n_var = 0
n_gg = 0

do im = 1, m_max
  if (m_cos(im) /= -1) then
    n_gg = n_gg + 1
    gg => gg1(n_gg)
    gg%sincos = cos$
    gg%m = m_cos(im)
    allocate (gg%deriv(iz_min:iz_max, 0:n_deriv_max), gg%ix_var(0:n_deriv_max))
    gg%deriv = 0; gg%ix_var = -1
    gg%sym_score = sym_score_calc (gg%sincos, gg%m)
    do id = 0, n_deriv_max
      if (gg%m == 0 .and. id == 0) cycle
      n_var = n_var + 1
      gg%ix_var(id) = n_var
    enddo
  endif

  if (m_sin(im) /= -1) then
    n_gg = n_gg + 1
    gg => gg1(n_gg)
    gg%sincos = sin$
    gg%m = m_sin(im)
    allocate (gg%deriv(iz_min:iz_max, 0:n_deriv_max), gg%ix_var(0:n_deriv_max))
    gg%deriv = 0; gg%ix_var = -1
    gg%sym_score = sym_score_calc (gg%sincos, gg%m)
    do id = 0, n_deriv_max
      n_var = n_var + 1
      gg%ix_var(id) = n_var
    enddo
  endif
enddo

!

allocate (fit(iz_min:iz_max))
allocate (var_vec(n_var), vec0(n_var), weight(n_merit))
vec0 = 0
var_vec = 0

do iz = iz_min, iz_max, every_n_th_plane
  print *, '===================================='
  print '(a, f10.4, a, i0, a)', ' Plane:', iz*del_grid(3), ' (', iz, ')'
  print *, 'Cycle   Merit (rms)'

  iz0 = max(Nz_min, iz-n_planes_add)
  iz1 = min(Nz_max, iz+n_planes_add)

  call initial_lmdif
  var_vec = vec0
  call merit_calc(vec0, iz, .true.)
  merit0_vec = merit_vec
  print '(i5, es14.6)', 0, merit
  merit0 = merit
  fit(iz)%rms0 = merit0 

  ! 

  select case (optimizer)
  case ('lmdif')
    do iloop = 1, n_cycles
      call suggest_lmdif (var_vec, merit_vec, lmdif_eps, n_cycles, at_end)
      call merit_calc(var_vec, iz, .false.)
      if (merit < 0.99*merit0) then
        print '(i5, es14.6)', iloop, merit
        merit0 = merit
      endif
      if (at_end) exit
    enddo

  case ('lm')
    a_lambda = -1
    weight = 1
    do iloop = 1, n_cycles
      call super_mrqmin (merit0_vec, weight, var_vec, chisq, mrqmin_func, storage, a_lambda, status)
      if (merit < 0.99*merit0) then
        print '(i5, es14.6, es12.3)', iloop, merit, a_lambda
        merit0 = merit
      endif
      if (a_lambda > 1e20_rp) exit
    enddo

  case default
    print *, 'Unknown optimizer: ' // optimizer
    print *, 'Possibilities are "lmdif", "lm"'
    stop
  end select

  !

  fit(iz)%rms_fit = merit

  do ig = 1, size(gg1)
    gg => gg1(ig)
    do id = 0, n_deriv_max
      if (gg%ix_var(id) < 1) cycle
      gg%deriv(iz,id) = var_vec(gg%ix_var(id))
    enddo
  enddo

  print '(a)', '  Ix    m   nd  sym             Deriv       Deriv*r_max^(m+d-1)*(d+m)*m!/((d/2)!*(d/2+m)!)'

  do ig = 1, size(gg1)
    gg => gg1(ig)
    do id = 0, n_deriv_max
      if (gg%ix_var(id) > 0) then
        v = gg%deriv(iz,id)
        print '(i4, 3i5, 2x, a, es16.6, es12.2)', gg%ix_var(id), gg%m, id, &
                      gg%sym_score, sincos_name(gg%sincos), v, &
                      v * r_max**(gg%m+id-1) * (id+gg%m) * factorial(gg%m) / (factorial(id/2) * factorial(id/2+gg%m))
      else
        print '(i4, 3i5, 2x, a, es16.6)', gg%ix_var(id), gg%m, id, gg%sym_score, sincos_name(gg%sincos)
      endif
    enddo
  enddo

  if (printit) then
    do ix = Nx_min, Nx_max
      x = ix * del_grid(1)
      do iy = Ny_min, Ny_max
        y = iy * del_grid(2)
        print '(2i4, 2f8.4, 3(4x, 2f8.4))', ix, iy, x, y, &
                      Bx_dat(ix,iy,iz), Bx_fit(ix,iy), &
                      By_dat(ix,iy,iz), By_fit(ix,iy), &
                      Bz_dat(ix,iy,iz), Bz_fit(ix,iy)
      enddo
    enddo
  endif
enddo

!-------------------------------------------------------------------
contains

subroutine merit_calc(var_vec, iz, dB_calc, b_fit_vec)

type (gg1_struct), pointer :: gg

real(rp) var_vec(:)
real(rp), optional :: b_fit_vec(:)
real(rp) x, y, rho, theta, Bx, By, Bz, B_rho, B_theta, f, ff, p_rho, coef

integer i, iz, ix, iy, iv, m, id, nn, izz, ig, n0, n3, ix_merit
logical :: dB_calc, is_even

!

ix_merit = 0

do izz = iz0, iz1

  n3 = size(Bx_fit)
  n0 = ix_merit * n3

  Bx_fit = 0
  By_fit = 0
  Bz_fit = 0

  do ig = 1, size(gg1)
    gg => gg1(ig)
    do id = 0, n_deriv_max
      iv = gg%ix_var(id)
      if (iv < 1) cycle

      dBx_dvar = 0
      dBy_dvar = 0
      dBz_dvar = 0

      m = gg%m
      is_even = (modulo(id,2) == 0)

      if (is_even) then
        nn = id / 2
      else
        nn = (id - 1) / 2
      endif

      coef = poly_eval(var_vec(iv:iv+n_deriv_max-id), (izz-iz)*del_grid(3), .true.)

      f = (-0.25_rp)**nn * factorial(m) / (factorial(nn) * factorial(nn+m))

      do ix = Nx_min, Nx_max
        x = ix * del_grid(1)
        do iy = Ny_min, Ny_max
          y = iy * del_grid(2)
          rho = sqrt(x*x + y*y)
          theta = atan2(y,x)

          if (id+m-1 == 0) then  ! Covers case where rho = 0
            p_rho = 1
          else
            p_rho = rho**(id+m-1)
          endif

          ff = f * p_rho 

          if (is_even) then
            if (gg%sincos == sin$) then
              B_rho   = ff * (2*nn+m) * sin(m*theta)
              B_theta = ff * m * cos(m*theta)
            else
              B_rho   = ff * (2*nn+m) * cos(m*theta)
              B_theta = -ff * m * sin(m*theta)
            endif

            Bx = B_rho * cos(theta) - B_theta * sin(theta)
            By = B_rho * sin(theta) + B_theta * cos(theta)

            Bx_fit(ix,iy) = Bx_fit(ix,iy) + Bx * coef
            By_fit(ix,iy) = By_fit(ix,iy) + By * coef

            if (dB_calc) then
              dBx_dvar(ix,iy) = Bx
              dBy_dvar(ix,iy) = By
            endif

          else
            if (gg%sincos == sin$) then
              Bz = ff * sin(m*theta)
            else
              Bz = ff * cos(m*theta)
            endif

            Bz_fit(ix,iy) = Bz_fit(ix,iy) + Bz * coef
            if (dB_calc) then
              dBz_dvar(ix,iy) = Bz
            endif
          endif
        enddo   ! ix = Nx_min, Nx_max
      enddo   ! iy = Ny_min, Ny_max

      if (dB_calc) then
        dB_dvar_vec(n0+1:n0+n3, iv)        = reshape(dBx_dvar(:,:), [n3])
        dB_dvar_vec(n0+n3+1:n0+2*n3, iv)   = reshape(dBy_dvar(:,:), [n3])
        dB_dvar_vec(n0+2*n3+1:n0+3*n3, iv) = reshape(dBz_dvar(:,:), [n3])
      endif

    enddo   ! id = 0, n_deriv_max
  enddo   ! ig = 1, size(gg1)

  !

  merit_vec(n0+1:n0+n3)        = reshape(Bx_dat(:,:,izz)-Bx_fit, [n3])
  merit_vec(n0+n3+1:n0+2*n3)   = reshape(By_dat(:,:,izz)-By_fit, [n3])
  merit_vec(n0+2*n3+1:n0+3*n3) = reshape(Bz_dat(:,:,izz)-Bz_fit, [n3])

  if (present(b_fit_vec)) then
    b_fit_vec(n0+1:n0+n3)        = reshape(Bx_fit, [n3])
    b_fit_vec(n0+n3+1:n0+2*n3)   = reshape(By_fit, [n3])
    b_fit_vec(n0+2*n3+1:n0+3*n3) = reshape(Bz_fit, [n3])
  endif

  ix_merit = ix_merit + 3

enddo   ! izz = iz0, iz1

merit = 3 * norm2(merit_vec) / (n3 * ix_merit)

end subroutine merit_calc

!---------------------------------------------------------------------------
! contains

subroutine mrqmin_func (var, yfit, dy_dvar, status)

real(rp), intent(in) :: var(:)
real(rp), intent(out) :: yfit(:)
real(rp), intent(out) :: dy_dvar(:,:)
integer status

!

call merit_calc(var, iz, .false., yfit)
dy_dvar = dB_dvar_vec

end subroutine mrqmin_func

end subroutine fit_field

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

function sincos_name(sincos) result (name)

integer sincos
character(3) name

!

select case (sincos)
case (sin$);  name = 'sin'
case (cos$);  name = 'cos'
case default; name = '???'
end select

end function sincos_name

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

function sym_score_calc (sincos, m) result (sym_score)

integer sincos, m, sym_score

!

if (sym_x == 0 .and. sym_y == 0) then
  sym_score = 1

elseif (sym_x == 0) then
  if (sincos == sin$) then
    sym_score = 1 - sym_y
  else
    sym_score = 1 + sym_y
  endif

elseif (sym_y == 0) then
  if (sincos == sin$) then
    sym_score = 1 - (-1)**m * sym_x
  else
    sym_score = 1 + (-1)**m * sym_x
  endif

else  ! sym in both planes
  if (sincos == sin$) then
    sym_score = 1 - sym_y - (-1)**m * sym_x + (-1)**m * sym_x * sym_y
  else
    sym_score = 1 + sym_y + (-1)**m * sym_x + (-1)**m * sym_x * sym_y
  endif
endif

end function sym_score_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

subroutine write_gg()

type (gg1_struct), pointer :: gg
integer ig, iz, m, n, nmax
character(40) fmt

!

if (out_file == '') out_file = 'gg'

!--------------------------

open (1, file = trim(out_file) // '.plot')

write (1, '(4(a, f8.4))')  '# del_grid = [', del_grid(1), ',', del_grid(2), ',', del_grid(3), ']'
write (1, '(4(a, f8.4))')  '# r0_grid  = [', r0_grid(1), ',', r0_grid(2), ',', r0_grid(3), ']'
write (1, '(a, f12.6, a)') '# r_max = ', r_max, '  ! Max transverse radius'

do ig = 1, size(gg1)
  gg => gg1(ig)
  if (ig > 1) then
    write (1, *)    ! Two blank lines is to separate the data sets for gnuplot plotting
    write (1, *)
  endif
  write (1, '(a, i2)') '# m    =', gg%m
  write (1, '(2a)')    '# type = ', sincos_name(gg%sincos)
  write (1, '(a, i2)') '# Iz     z_pos    Derivs'

  do iz = iz_min, iz_max
    write (1, '(i4, f10.4, 99es16.8)') iz, iz*del_grid(3)+r0_grid(3), fit(iz)%rms0, gg%deriv(iz,:)
  enddo
enddo

close (1)

!--------------------------

open (1, file = trim(out_file) // '.info')

write (1, '(4(a, f8.4))')  '# del_grid = [', del_grid(1), ',', del_grid(2), ',', del_grid(3), ']'
write (1, '(4(a, f8.4))')  '# r0_grid  = [', r0_grid(1), ',', r0_grid(2), ',', r0_grid(3), ']'
write (1, '(a, f12.6, a)') '# r_max = ', r_max, '  ! Max transverse radius'

do ig = 1, size(gg1)
  gg => gg1(ig)
  m = gg%m
  if (ig > 1) then
    write (1, *)    ! Two blank lines is to separate the data sets for gnuplot plotting
    write (1, *)
  endif
  write (1, '(a, i2)') '# m    =', gg%m
  write (1, '(2a)')    '# type = ', sincos_name(gg%sincos)
  write (1, '(a, i2)') '# Iz     z_pos    Init_RMS     RMS/RMS0     Deriv * r_max^(m+d-1) * (d+m) * m! / ((d/2)! * (d/2+m)!)'

  nmax = ubound(gg%deriv,2)
  do iz = iz_min, iz_max
    write (1, '(i4, f10.4, 2es13.4, 99es16.8)') iz, iz*del_grid(3)+r0_grid(3), &
                  fit(iz)%rms0, fit(iz)%rms_fit/fit(iz)%rms0, &
                  (gg%deriv(iz,n) * r_max**(max(0,m+n-1)) * (n+m) * factorial(m) / (factorial(n/2) * factorial(n/2+m)), n = 0, nmax)
  enddo

enddo

close (1)

!--------------------------

open (1, file = trim(out_file) // '.bmad', recl = 500)

write (1, '(a)')           '  gen_grad_map = {'
write (1, '(a, f10.6, a)') '    dz =', del_grid(3), ','

do ig = 1, size(gg1)
  gg => gg1(ig)
  write (1, '(a)')        '    curve = {'
  write (1, '(a, i2, a)') '      m    =', gg%m, ','
  write (1, '(3a)')       '      kind = ', sincos_name(gg%sincos), ','
  write (1, '(a, i2)')    '      derivs = {'

  write (fmt, '(a, i0, a)') '(f13.4, a, ', size(gg%deriv,2), 'es20.12, a)' 
  do iz = iz_min, iz_max-1
    write (1, fmt) iz*del_grid(3), ':', gg%deriv(iz,:), ','
  enddo
  write (1, fmt) iz_max*del_grid(3), ':', gg%deriv(iz_max,:)
  write (1, '(a)') '      }'
  if (ig == size(gg1)) then
    write (1, '(a)') '    }'
  else
    write (1, '(a)') '    },'
  endif
enddo

write (1, '(a)') '  }'

close (1)

end subroutine write_gg

end module
