module gg_fit_mod

use lmdif_mod
implicit none

!

type var_info_struct
  integer cossin
  integer m
  integer nd    ! deriv
end type

type (var_info_struct), allocatable :: var_info(:)
type (lat_struct) lat

integer, parameter :: m_max = 10
integer :: Nx_min, Nx_max, Ny_min, Ny_max, Nz_min, Nz_max
integer :: n_grid_pts
integer :: n_cycles
integer :: every_n_th_plane
integer :: n_deriv_max
integer :: ix_z_min = -1, ix_z_max = -1    ! Compute range
integer :: m_cos(m_max) = -1
integer :: m_sin(m_max) = -1
integer :: Bx_sym_x = 0, Bx_sym_y = 0, By_sym_x = 0, By_sym_y = 0, Bz_sym_x = 0, Bz_sym_y = 0
integer n_var, n_merit

real(rp), allocatable :: Bx_fit(:,:,:), By_fit(:,:,:), Bz_fit(:,:,:)
real(rp), allocatable :: Bx_dat(:,:,:), By_dat(:,:,:), Bz_dat(:,:,:)
real(rp), allocatable :: Bx_diff(:,:), By_diff(:,:), Bz_diff(:,:)

real(rp) :: del_grid(3), field_scale, length_scale
real(rp) :: lmdif_eps = 1d-12

real(rp), allocatable :: var_vec(:), merit_vec(:)

real(rp), allocatable :: div(:,:,:), curl_x(:,:,:), curl_y(:,:,:), curl_z(:,:,:)
real(rp) div_max, div_scale_max, div_scale, r0_grid(3)
real(rp) curl_x_max, curl_x_scale_max, curl_x_scale
real(rp) curl_y_max, curl_y_scale_max, curl_y_scale
real(rp) curl_z_max, curl_z_scale_max, curl_z_scale
real(rp) B_diff, merit_coef, merit_k, merit_data, dB_rms
real(rp) merit_tot, coef_weight, k_weight
real(rp) merit_x, merit_y, merit_z, x_offset_map, y_offset_map, z_offset_map

character(100) :: field_file

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

real(rp), allocatable :: valid_field(:,:,:)
real(rp) xx, yy, zz, Bx, By, Bz

integer ios, ios1, ios2
integer :: i, j, k, ix, iy, iz

character(*) :: field_file
character(200) line


! Binary
! Notice that the binary table stores lengths in meters and fields in Tesla independent of length_scale and field_scale.

if (field_file(1:8) == 'binary::') then

  open (1, file = field_file(9:), form = 'unformatted', status = 'OLD', action = 'READ')
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
    read (1, iostat = ios) i, j, k, Bx_dat(i,j,k), By_dat(i,j,k), Bz_dat(i,j,k), valid_field(i,j,k)
    if (ios /= 0) exit
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
n_opti_grid_pts = count(opti_field)

print *, 'Field table read: ', trim(field_file)
print *, '  Number of grid field points:       ', n_grid_pts
print *, '  Number points used in optimization:', n_opti_grid_pts

if (n_grid_pts /= count(valid_field)) then
  print *, 'Number points missing in table: ', n_grid_pts - count(valid_field)
  print *, 'Will stop here!'
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

character(*) field_file

! Write binary table
! Notice that the binary table stores lengths in meters and fields in Tesla independent of length_scale and field_scale.

open (1, file = trim(field_file) // '.binary')

write (1) length_scale, field_scale
write (1) Nx_min, Nx_max
write (1) Ny_min, Ny_max
write (1) Nz_min, Nz_max
write (1) del_grid
write (1) r0_grid

do ix = Nx_min, Nx_max, Nx_max-Nx_min
do iy = Ny_min, Ny_max
do iz = Nz_min, Nz_max
  write (1) i, j, k, Bx_dat(i,j,k), By_dat(i,j,k), Bz_dat(i,j,k), valid_field(i,j,k)
enddo
enddo
enddo

end subroutine write_binary_field_table

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
!-

subroutine fit_field_table()

integer ixz, ixz0, ixz1, iloop, iv, im, id
logical at_end
!

n_var = (count(m_cos /= -1) + count(m_sin /= -1)) * (n_deriv_max + 1)
do im = 1, m_max
  if (m_cos(im) == 0) n_var = n_var - 1
  if (m_sin(im) == 0) then
    print *, 'M = 0 sin coefficients do not contribute to the field. Please remove m_sin set.'
    stop
  endif
enddo



n_merit = 3*(Nx_max-Nx_min+1) * (Ny_max-Ny_min+1)

allocate (var_vec(n_var), var_info(n_var)
allocate (merit_vec(n_merit), diff(Nx_min:Nx_max, Ny_min:Ny_max))

allocate (Bx_diff(Nx_max:Nx_min, Ny_max:Ny_min), By_diff(Nx_max:Nx_min, Ny_max:Ny_min), Bz_diff(Nx_max:Nx_min, Ny_max:Ny_min))

!

ixz0 = ix_z_min
if (ixz0 == -1) ixz0 = lbound(Bx_fit,3)

ixz1 = ix_z_max
if (ixz1 == -1) ixz1 = ubound(Bx_fit,3)

var_vec = 0
iv = 0
do im = 1, m_max
  if (m_cos(im) /= -1) then
    do id = 0, n_deriv_max
      if (m_cos(im) == 0 .and. id == 0) cycle
      iv = iv + 1
      var_info(iv) = var_info_struct(cos$, m_cos(im), id)
    enddo
  endif

  if (m_sin(im) /= -1) then
    do id = 0, n_deriv_max
      iv = iv + 1
      var_info(iv) = var_info_struct(sin$, m_sin(im), id)
    enddo
  endif
enddo

do ixz = ixz0, ixz1
  call initial_lmdif
  call merit_calc()

  do iloop = 1, 10000
    call suggest_lmdif (var_vec, merit_vec, lmdif_eps, 10000, at_end
    call merit_calc()
    print '(i5, es14.4)', iloop, merit
    if (at_end) exit
  enddo
enddo


!----------------------------------
contains

subroutine merit_calc()

type (var_info_struct) info
real(rp) x, y, rho, theta, Bxy(2), Bz
integer ix, iy, iv, m, nd, d2
logical is_even

!

Bx_diff = Bx_data(:,:,ixz)
By_diff = By_data(:,:,ixz)
Bz_diff = Bz_data(:,:,ixz)

do iv = 1, n_var
  info = var_info(iv)
  nd = info%nd
  m = info%m
  is_even = (modulo(nd,2) == 0)

  if (is_even) then
    d2 = nd2 / 2
    f = (-0.25_rp)**d2 * factorial(m) * (2*d2+m) / (factorial(d2) * factorial(d2+m))
  else
    d2 = (nd2 - 1) / 2
    f = (-0.25_rp)**d2 * factorial(m) / (factorial(d2) * factorial(d2+m))
  endif

  do ix = Nx_min, Nx_max
    x = ix * del_grid(1)
    do iy = Ny_min, Ny_max
      y = iy * del_grid(2)
      rho = sqrt(x*x + y*y)
      theta = atan2(y,x)

      if (is_even) then
        if (m == 0) then  ! Must be cossin = cos$ and nd > 0
          Bxy = (f * rho**(nd+m-1) * var_vec(iv)) * [cos(theta), sin(theta)]
        elseif (info%cossin == sin$) then
          Bxy = (f * rho**(nd+m-1) * var_vec(iv)) * [sin((m-1)*theta), cos((m-1)*theta)]
        else
          Bxy = (f * rho**(nd+m-1) * var_vec(iv)) * [cos((m-1)*theta), sin((m-1)*theta)]
        endif
        Bx_diff(ix,iy) = Bx_diff(ix,iy) - Bxy(1)
        By_diff(ix,iy) = By_diff(ix,iy) - Bxy(2)
      else
        if (info%cossin == sin$) then
          Bz = f * rho**(nd+m-1) * var_vec(iv) * sin(m*theta)
        else
          Bz = f * rho**(nd+m-1) * var_vec(iv) * cos(m*theta)
        endif
        Bz_diff(ix,iy) = Bz_diff(ix,iy) - Bz
      endif

    enddo
  enddo
enddo

nn = n_merit/3
merit_vec(1:nn)        = reshape(Bx_diff, [nn])
merit_vec(nn+1:2*nn)   = reshape(By_diff, [nn])
merit_vec(2*nn+1:)     = reshape(Bz_diff, [nn])

merit = norm2(merit_vec) / n_merit

end subroutine merit_calc

end subroutine fit_field_table()

end module
