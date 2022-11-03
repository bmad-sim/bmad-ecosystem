module gg_fit_mod

use lmdif_mod
use bmad_struct
use expression_mod, only: sin$, cos$

implicit none

!

type gg1_struct
  integer :: sincos = -1    ! sin$ or cos$
  integer :: m = -1         ! Azimuthal index
  real(rp), allocatable :: deriv(:,:)
end type

type (gg1_struct), allocatable, target :: gg1(:)

type var_info_struct
  integer :: sincos = -1    ! sin$ or cos$
  integer :: m = -1         ! Azimuthal index
  integer :: id = -1        ! derivative index
  integer :: sym_score = 1  ! Symmetry score. 0 => do not use in fit.
  integer :: ix_var = -1    ! Index to var_vec(:) array.
  integer :: ix_gg = -1     ! Index to gg1 array
end type

type (var_info_struct), allocatable, target :: var_info(:)

integer, parameter :: m_max = 10
integer :: Nx_min, Nx_max, Ny_min, Ny_max, Nz_min, Nz_max
integer :: n_grid_pts
integer :: n_cycles
integer :: every_n_th_plane
integer :: n_deriv_max
integer :: ix_z_min = -1, ix_z_max = -1    ! Compute range
integer :: m_cos(m_max) = -1
integer :: m_sin(m_max) = -1
integer :: sym_x = 0, sym_y = 0
integer n_var0, n_var, n_merit

real(rp), allocatable :: Bx_dat(:,:,:), By_dat(:,:,:), Bz_dat(:,:,:)
real(rp), allocatable :: Bx_fit(:,:), By_fit(:,:), Bz_fit(:,:)
real(rp), allocatable :: Bx_diff(:,:), By_diff(:,:), Bz_diff(:,:)

! del_grid is the difference between (x,y,z) points in the units of the grid file.
! length_scale is table file (x,y,z) units in meters. length_scale = table_units/meters
! field_scale is the field table field in Tesla or V/m.

real(rp) :: del_grid(3), r0_grid(3), field_scale, length_scale
real(rp) :: lmdif_eps = 1d-12

real(rp), allocatable :: var_vec(:), merit_vec(:)
real(rp) merit

logical printit

character(100) :: field_file , out_file = ''

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
!   field_file    -- character(*): Name of file.Y
!                      If file starts with the string "binary_' this is a binary file.
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

  close (1)

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
  read (line, *, iostat = ios2) Nx_min, Nx_max
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('Nx_min, Nx_max', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) Ny_min, Ny_max
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('Ny_min, Ny_max', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) Nz_min, Nz_max
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('Nz_min, Nz_max', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) del_grid
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('del_grid', line)
  read (1, '(a)', iostat = ios1) line
  read (line, *, iostat = ios2) r0_grid
  if (ios1 /= 0 .or. ios2 /= 0) call read_err('r0_grid', line)

  allocate (Bx_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), By_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))
  allocate (Bz_dat(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), valid_field(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))

  do
    read (1, '(a)', iostat = ios1) line
    if (line(1:1) == '!' .or. line == '') cycle
    read (line, *, iostat = ios2) xx, yy, zz, Bx, By, Bz
    if (ios1 < 0) exit
    if (ios1 /= 0 .or. ios2 /= 0) call read_err('xx, yy, zz, Bx, By, Bz', line)
    i = nint (xx/del_grid(1))
    j = nint (yy/del_grid(2))
    k = nint (zz/del_grid(3))

    if (i < Nx_min .or. i > Nx_max .or. j < Ny_min .or. j > Ny_max .or. k < Nz_min .or. k > Nz_max) then
      print *, 'COMPUTED INDEX OUT OF RANGE FOR LINE: ', trim(line)
      stop
    endif
    Bx_dat(i,j,k) = Bx * field_scale
    By_dat(i,j,k) = By * field_scale
    Bz_dat(i,j,k) = Bz * field_scale
    valid_field(i,j,k) = .true.
  end do

  close (1)

  del_grid = del_grid * length_scale
  r0_grid = r0_grid * length_scale
endif

!

n_grid_pts = (Nx_max-Nx_min+1) * (Ny_max-Ny_min+1) * (Nz_max-Nz_min+1)

print *, 'Field table read: ', trim(field_file)
print *, '  Number of grid field points:       ', n_grid_pts

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
! Routine to write in a field table. Units are: cm, Gauss.
! Look at the code for the format.
!
! Input:
!   field_file    -- character(*): Name of file.Y
!                      If file starts with the string "binary_' this is a binary file.
!-

subroutine write_binary_field_table (field_file)

integer ix, iy, iz
character(*) field_file

! Write binary table
! Notice that the binary table stores lengths in meters and fields in Tesla independent of length_scale and field_scale.

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

end subroutine write_binary_field_table

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

subroutine fit_field()

type (var_info_struct), pointer :: vinfo
real(rp), allocatable :: vec0(:)
real(rp) merit, merit0, x, y
integer ixz, ixz0, ixz1, iloop, iv0, im, id, n_gg
integer ix, iy, iz

logical at_end

!

n_merit = 3 * (Nx_max-Nx_min+1) * (Ny_max-Ny_min+1)
allocate (merit_vec(n_merit))

allocate (Bx_diff(Nx_min:Nx_max, Ny_min:Ny_max), By_diff(Nx_min:Nx_max, Ny_min:Ny_max), &
                                                            Bz_diff(Nx_min:Nx_max, Ny_min:Ny_max))
allocate (Bx_fit(Nx_min:Nx_max, Ny_min:Ny_max), By_fit(Nx_min:Nx_max, Ny_min:Ny_max), &
                                                            Bz_fit(Nx_min:Nx_max, Ny_min:Ny_max))

!

ixz0 = ix_z_min
if (ixz0 == -1) ixz0 = Nz_min

ixz1 = ix_z_max
if (ixz1 == -1) ixz1 = Nz_max

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

allocate (var_info(n_var0))
allocate (gg1(n_gg))

iv0 = 0
n_var = 0
n_gg = 0

do im = 1, m_max
  if (m_cos(im) /= -1) then
    n_gg = n_gg + 1
    gg1(n_gg)%sincos = cos$
    gg1(n_gg)%m = m_cos(im)
    allocate (gg1(n_gg)%deriv(ixz0:ixz1, 0:n_deriv_max))
    gg1(n_gg)%deriv = 0

    do id = 0, n_deriv_max
      if (m_cos(im) == 0 .and. id == 0) cycle
      var_info(iv0) = var_info_setup(cos$, m_cos(im), id, iv0, n_var, n_gg)
    enddo
  endif

  if (m_sin(im) /= -1) then
    n_gg = n_gg + 1
    gg1(n_gg)%sincos = sin$
    gg1(n_gg)%m = m_sin(im)
    allocate (gg1(n_gg)%deriv(ixz0:ixz1, 0:n_deriv_max))
    gg1(n_gg)%deriv = 0

    do id = 0, n_deriv_max
      var_info(iv0) = var_info_setup(sin$, m_sin(im), id, iv0, n_var, n_gg)
    enddo
  endif
enddo

!

allocate (var_vec(n_var), vec0(n_var))
vec0 = 0
var_vec = 0

do ixz = ixz0, ixz1
  call initial_lmdif
  call merit_calc(vec0, ixz)
  print '(i5, es14.6)', 0, merit
  merit0 = 2*merit

  do iloop = 1, n_cycles
    call suggest_lmdif (var_vec, merit_vec, lmdif_eps, n_cycles, at_end)
    call merit_calc(var_vec, ixz)
    if (merit < 0.99*merit0) then
      print '(i5, es14.6)', iloop, merit
      merit0 = merit
    endif
    if (at_end) exit
  enddo

  do iv0 = 1, size(var_info)
    vinfo => var_info(iv0)
    if (vinfo%ix_var < 1) cycle
    gg1(vinfo%ix_gg)%deriv(ixz,vinfo%id) = var_vec(vinfo%ix_var)
  enddo

  print *
  print *, 'Plane:', ixz
  print *, ' Ix    m  der  sym            Coef'
  do iv0 = 1, size(var_info)
    vinfo => var_info(iv0)
    if (vinfo%ix_var > 0) then
      print '(i4, 3i5, 2x, a, es16.6)', iv0, vinfo%m, vinfo%id, &
                              vinfo%sym_score, sincos_name(vinfo%sincos), var_vec(vinfo%ix_var)
    else
      print '(i4, 3i5, 2x, a, es16.6)', iv0, vinfo%m, vinfo%id, &
                              vinfo%sym_score, sincos_name(vinfo%sincos)
    endif
  enddo

  if (printit) then
    do ix = Nx_min, Nx_max
      x = ix * del_grid(1)
      do iy = Ny_min, Ny_max
        y = iy * del_grid(2)
        print '(2i4, 2f8.4, 3(4x, 2f8.4))', ix, iy, x, y, &
                      Bx_dat(ix,iy,ixz), Bx_fit(ix,iy), &
                      By_dat(ix,iy,ixz), By_fit(ix,iy), &
                      Bz_dat(ix,iy,ixz), Bz_fit(ix,iy)
      enddo
    enddo
  endif
enddo

!----------------------------------
contains

subroutine merit_calc(var_vec, ixz)

type (var_info_struct) info
real(rp) var_vec(:)
real(rp) x, y, rho, theta, Bx, By, Bz, B_rho, B_theta, f, ff, p_rho
integer ixz, ix, iy, iv0, iv, m, id, nn, n_merit, n3
logical is_even

!

Bx_fit = 0
By_fit = 0
Bz_fit = 0

do iv0 = 1, size(var_info)
  info = var_info(iv0)
  iv = info%ix_var
  if (iv < 1) cycle
  if (var_vec(iv) == 0) cycle

  id = info%id
  m = info%m
  is_even = (modulo(id,2) == 0)

  if (is_even) then
    nn = id / 2
  else
    nn = (id - 1) / 2
  endif

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

      ff = f * p_rho * var_vec(iv)

      if (is_even) then
        if (info%sincos == sin$) then
          B_rho   = ff * (2*nn+m) * sin(m*theta)
          B_theta = ff * m * cos(m*theta)
        else
          B_rho   = ff * (2*nn+m) * cos(m*theta)
          B_theta = -ff * m * sin(m*theta)
        endif
        Bx = B_rho * cos(theta) - B_theta * sin(theta)
        By = B_rho * sin(theta) + B_theta * cos(theta)
        Bx_fit(ix,iy) = Bx_fit(ix,iy) + Bx
        By_fit(ix,iy) = By_fit(ix,iy) + By

      else
        if (info%sincos == sin$) then
          Bz = ff * sin(m*theta)
        else
          Bz = ff * cos(m*theta)
        endif
        Bz_fit(ix,iy) = Bz_fit(ix,iy) + Bz
      endif

    enddo
  enddo
enddo

n_merit = size(merit_vec)
n3 = n_merit/3
merit_vec(1:n3)        = reshape(Bx_dat(:,:,ixz)-Bx_fit, [n3])
merit_vec(n3+1:2*n3)   = reshape(By_dat(:,:,ixz)-By_fit, [n3])
merit_vec(2*n3+1:)     = reshape(Bz_dat(:,:,ixz)-Bz_fit, [n3])

merit = sum(sqrt(merit_vec*merit_vec)) / n3

end subroutine merit_calc

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

function var_info_setup (sincos, m, id, iv0, n_var, n_gg) result (var_info)

type (var_info_struct) var_info
integer sincos, m, id, iv0, n_var, n_gg, sym_score
logical is_even

!

is_even = (modulo(id,2) == 0)
iv0 = iv0 + 1

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

if (sym_score == 0) then 
  var_info = var_info_struct(sincos, m, id, sym_score, -1, -1)
else
  n_var = n_var + 1
  var_info = var_info_struct(sincos, m, id, sym_score, n_var, n_gg)
endif

end function var_info_setup

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

subroutine write_gg()

type (gg1_struct), pointer :: gg
integer ig, iz

if (out_file == '') out_file = 'gg.dat'
open (1, file = out_file, status = 'new')

write (1, '(4(a, f8.4))') '# del_grid = [', del_grid(1), ',', del_grid(2), ',', del_grid(3), ']'
write (1, '(4(a, f8.4))') '# r0_grid  = [', r0_grid(1), ',', r0_grid(2), ',', r0_grid(3), ']'

do ig = 1, size(gg1)
  gg => gg1(ig)
  if (ig > 1) then
    write (1, *)    ! Two blank lines is to separate the data sets for gnuplot plotting
    write (1, *)
  endif
  write (1, '(a, i2)') '#  m    =', gg%m
  write (1, '(a, i2)') '#  type =', sincos_name(gg%sincos)
  write (1, '(a, i2)') '# Iz     z_pos   Derivs'

  do iz = Nz_min, Nz_max
    write (1, '(i4, f10.4, 99es20.12)') iz, iz*del_grid(3)+r0_grid(3), gg%deriv(iz,:)
  enddo
enddo

close (1)

end subroutine write_gg

end module
