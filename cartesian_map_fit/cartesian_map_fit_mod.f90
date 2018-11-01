module wig_mod

use bmad

integer :: Nx, Ny, Nz, Nz_ref
integer :: n_pts, n_term, n_data, n_var, n_var_term
integer :: n_loops, n_cycles
integer :: iz1, iz2  ! z-edges of the fit region


type bf_struct
  type (cartesian_map_term1_struct), tm
  real(rp) dky_dkx, dky_dkz
  real(rp) coef_x, coef_y, coef_z
  real(rp), pointer :: s_x(:), c_x(:), s_y(:), c_y(:), s_z(:), c_z(:)
  real(rp), pointer :: ds_y(:), dc_y(:)
  integer sign_x, sign_dc_x
end type     


type (bf_struct) bf(500)
type (term_struct) term(500), old_term(500)

real(rp), allocatable :: Bx_fit(:,:,:), By_fit(:,:,:), Bz_fit(:,:,:)
real(rp), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)

real(rp) :: del(3), sumB

real(rp), allocatable :: div(:,:,:), curl_x(:,:,:), curl_y(:,:,:), curl_z(:,:,:)
real(rp) div_max, div_scale_max, div_scale
real(rp) curl_x_max, curl_x_scale_max, curl_x_scale
real(rp) curl_y_max, curl_y_scale_max, curl_y_scale
real(rp) curl_z_max, curl_z_scale_max, curl_z_scale
real(rp) B_int, B_diff, merit_coef, merit_data, dB_rms
real(rp) :: merit_tot, weight_coef_all = 1.0, weight_coef(400)
real(rp) merit_x, merit_y, merit_z

real(rp), allocatable :: alpha_wf(:,:), beta_wf(:)

logical :: dyda_calc = .true.
logical :: fft_write = .false.
logical :: mask_x0 = .false., mask_y0 = .false., mask_phi_z = .false.
logical :: mask_kx = .false., mask_ky = .false., mask_kz = .false.

character(80) :: field_file
character(16) :: symmetry

real(rp), allocatable :: Bx_in(:,:,:), By_in(:,:,:), Bz_in(:,:,:)
real(rp), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:) 
integer Nx_in, Ny_in, Nz_in

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine read_field_table_file (field_file)
! 
! Routine to read in a field table. Units are: cm, Gauss.
! Look at the code for the format.
!
! Input:
!   field_file    -- character(*): Name of file.Y
!                      If file starts with the string "binary_' this is a binary file.
!-

subroutine read_field_table_file (field_file)

implicit none

integer :: i, j, k, ix, iy, iz

character(*) :: field_file
character(80) line

character(8) :: columns(7) = (/ '1 X   ', '2 Y   ', '3 Z   ', &
                                '4 BX  ', '5 BY  ', '6 BZ  ', '0 [CM]' /)


if (index(field_file, 'binary_') /= 0) then

  open (1, file = field_file, status = 'old', form = 'unformatted', readonly)
  read (1) Nx_in, Ny_in, Nz_in, del
  allocate (Bx_in(Nx_in,Ny_in,Nz_in), By_in(Nx_in,Ny_in,Nz_in), Bz_in(Nx_in,Ny_in,Nz_in))
  allocate (X(Nx_in,Ny_in,Nz_in), Y(Nx_in,Ny_in,Nz_in), Z(Nx_in,Ny_in,Nz_in))

  do k = 1, Nz_in
    read (1) Bx_in(:,:,k)
    read (1) By_in(:,:,k)
    read (1) Bz_in(:,:,k)
  end do

  close (1)

else

  open (1, file = field_file, STATUS='OLD', readonly, shared)
  read (1, *) Nx_in, Ny_in, Nz_in
  allocate (Bx_in(Nx_in,Ny_in,Nz_in), By_in(Nx_in,Ny_in,Nz_in), Bz_in(Nx_in,Ny_in,Nz_in))
  allocate (X(Nx_in,Ny_in,Nz_in), Y(Nx_in,Ny_in,Nz_in), Z(Nx_in,Ny_in,Nz_in))
             
  do i = 1, size(columns)
    read (1, '(a)') line
    if (index(line, trim(columns(i))) == 0) then
      print *, 'ERROR: COLUMNS IN FIELD TABLE NOT IN CORRECT ORDER. EXPECTED: ', columns(i)
      print *, '       GOT: ', trim(line)
      call err_exit
    endif
  enddo

  do i = 1, Nx_in
  do j = 1, Ny_in
  do k = 1, Nz_in
    read (1, *) X(i,j,k), Y(i,j,k), Z(i,j,k), Bx_in(i,j,k), By_in(i,j,k), Bz_in(i,j,k)
  end do
  end do
  end do

  close (1)

  del(1) = (X(Nx_in,1,1) - X(1,1,1)) / (Nx_in - 1)
  del(2) = (Y(1,Ny_in,1) - Y(1,1,1)) / (Ny_in - 1)
  del(3) = (Z(1,1,Nz_in) - Z(1,1,1)) / (Nz_in - 1)

  do i = 1, Nx_in
  do j = 1, Ny_in
  do k = 1, Nz_in
    if (abs(X(i,j,k) - (i-1)*del(1)) > 1e-5 * del(1)) then
      print *, 'ERROR IN ORDERING OF FIELD TABLE'
      print *, 'ERROR: X WRONG FOR:', i, j, k
      print *, '       AT:      ', X(i,j,k)-X(1,1,1), X(i,j,k), X(1,1,1)
      print *, '       COMPUTED:', del(1)*(i-1)
      call err_exit
    endif
    if (abs(Y(i,j,k) - (j-1)*del(2)) > 1e-5 * del(2)) then
      print *, 'ERROR IN ORDERING OF FIELD TABLE'
      print *, 'ERROR: Y WRONG FOR:', i, j, k
      print *, '       AT:      ', Y(i,j,k)-Y(1,1,1), Y(i,j,k), Y(1,1,1)
      print *, '       COMPUTED:', del(2)*(j-1)
      call err_exit
    endif
    if (abs(Z(i,j,k) - (k-1)*del(3)) > 1e-5 * del(3)) then
      print *, 'ERROR IN ORDERING OF FIELD TABLE'
      print *, 'ERROR: Z WRONG FOR:', i, j, k
      print *, '       AT:      ', Z(i,j,k)-Z(1,1,1), Z(i,j,k), Z(1,1,1)
      print *, '       COMPUTED:', del(3)*(k-1)
      call err_exit
    endif
  end do
  end do
  end do

  del = 0.01 * del  ! convert to meters

  Bx_in = Bx_in / 1e4  ! convert Gauss to Tesla
  By_in = By_in / 1e4
  Bz_in = Bz_in / 1e4

endif

end subroutine read_field_file

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine read_wiggler_fit_param_file (param_file, write_curl)
!
! Routine to read in the parameters for the program
!
! Input:
!   param_file    -- character(*): Name of the parameter file.
!   write_curl    -- logical: Write computed curl of field? Should be zero.
!-

subroutine read_wiggler_fit_param_file (param_file, write_curl)

implicit none

real(rp) ave_x, ave_z

integer :: i, j, k, ix, iy, iz

logical write_curl

! sorted_field_trimmed.txt'

character(*) :: param_file

namelist / wght / weight_coef_all, n_loops, n_cycles, &
    weight_coef, vary_phi_z, fft_write, Nx, Ny, iz1, iz2, symmetry

namelist / parameters / field_file, Nz_ref, term

! Read in parameters and starting fit

weight_coef = -1
symmetry = 'NONE'
open (1, file = 'fit.weight', status = 'old', readonly, shared)
read (1, nml = wght)
close (1)
where (weight_coef < 0) weight_coef = weight_coef_all

term%coef = 0
open (1, file = param_file, status = 'old', readonly, shared)
read (1, nml = parameters)
close (1)

! Read in field data file

call read_field_table_file (field_file)

if (Nx_in < Nx .or. Ny_in < Ny) then
  print *, 'ERROR: DATA FILE SIZE MISMATCH'
  print *, '       FROM PROGRAM Nx, Ny:  ', Nx, Ny
  print *, '       FROM DATA FILE:', Nx_in, Ny_in
  call err_exit
endif

!-------------------------------------
! read in starting point
! first throw out any term with a zero coef.

n_term = 0
do i = 1, size(term)
  if (term(i)%coef /= 0) then
    n_term = n_term + 1
    term(n_term) = term(i) 
  endif
enddo

! allocate

if (Nx == 0) Nx = Nx_in
if (Ny == 0) Ny = Ny_in

if (iz1 == 0) iz1 = 1
if (iz2 == 0) iz2 = Nz_in
Nz = iz2 - iz1 + 1

if (vary_phi_z) then
  n_var_term = 4  ! number of vars per term
else
  n_var_term = 3
endif

n_var = n_var_term * n_term
n_pts = Nx * Ny * Nz
n_data = 3*n_pts + n_term

allocate (alpha_wf(n_var, n_var), beta_wf(n_var))

do i = 1, size(bf)
  allocate(bf(i)%s_x(Nx), bf(i)%c_x(Nx))
  allocate(bf(i)%s_y(Ny), bf(i)%c_y(Ny))
  allocate(bf(i)%ds_y(Ny), bf(i)%dc_y(Ny))
  allocate(bf(i)%s_z(Nz), bf(i)%c_z(Nz))
enddo

allocate (curl_x(Nx, Ny, Nz), curl_y(Nx, Ny, Nz), curl_z(Nx, Ny, Nz))
allocate (div(Nx, Ny, Nz))

allocate (Bx_fit(Nx, Ny, Nz), By_fit(Nx, Ny, Nz), Bz_fit(Nx, Ny, Nz))
allocate (Bx(Nx, Ny, Nz), By(Nx, Ny, Nz), Bz(Nx, Ny, Nz))

Bx = Bx_in(1:Nx, 1:Ny, iz1:iz2) 
By = By_in(1:Nx, 1:Ny, iz1:iz2)
Bz = Bz_in(1:Nx, 1:Ny, iz1:iz2)

!

! This is for seeing how well the curl of the field is zero.

If (write_curl) then
  iy = 1
  do ix = 1, Nx
    ave_x = sum(abs(Bx(ix, iy, 1:Nz)))
    ave_z = sum(abs(Bz(ix, iy, 1:Nz))) 
    print '(a, i4, 2f10.2)', ' Bx (y = 0):', ix, &
                              100*(ix-1)*del(1), 1e4*ave_x/Nz
    print '(a, i4, 2f10.2)', ' Bz (y = 0):', ix, &
                              100*(ix-1)*del(1), 1e4*ave_z/Nz
  enddo

  ix = 1
  do iy = 1, Ny
    ave_x = sum(abs(Bx(ix, iy, 1:Nz)))
    print '(a, i4, 2f10.2)', ' Bx (x = 0):', iy, &
                                100*(iy-1)*del(2), 1e4*ave_x/Nz
  enddo

  call curl_calc (Bx, By, Bz)

endif

!----------------------------------------------------------------------

sumB = sum(abs(Bx) + abs(By) + abs(Bz))

end subroutine read_wiggler_fit_param_file
               
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine curl_calc (B_x, B_y, B_z)

implicit none

integer ix0, iy0, iz0
integer ix, iy, iz

real(rp) B_x(:,:,:), B_y(:,:,:), B_z(:,:,:)
real(rp) dx, dy, dz

div_max = 0; div_scale_max = 0
curl_x_max = 0; curl_x_scale_max = 0
curl_y_max = 0; curl_y_scale_max = 0
curl_z_max = 0; curl_z_scale_max = 0

do ix = 1, Nx-1
do iy = 1, Ny-1
do iz = 1, Nz-1 ! , Nz-1

  dx = (sum(B_x(ix+1,iy:iy+1,iz:iz+1)) - sum(B_x(ix,iy:iy+1,iz:iz+1)))/del(1)
  dy = (sum(B_y(ix:ix+1,iy+1,iz:iz+1)) - sum(B_y(ix:ix+1,iy,iz:iz+1)))/del(2)
  dz = (sum(B_z(ix:ix+1,iy:iy+1,iz+1)) - sum(B_z(ix:ix+1,iy:iy+1,iz)))/del(3)

  div_scale = (abs(dx) + abs(dy) + abs(dz)) / 4
  div(ix,iy,iz) = (dx + dy + dz) / 4

!

  dy = (sum(B_z(ix,iy+1,iz:iz+1)) - sum(B_z(ix,iy,iz:iz+1)))/del(2)
  dz = (sum(B_y(ix,iy:iy+1,iz+1)) - sum(B_y(ix,iy:iy+1,iz)))/del(3)

  curl_x_scale = (abs(dy) + abs(dz)) / 2
  curl_x(ix,iy,iz) = (dy - dz) / 2

!

  dx = (sum(B_z(ix+1,iy,iz:iz+1)) - sum(B_z(ix,iy,iz:iz+1)))/del(1)
  dz = (sum(B_x(ix:ix+1,iy,iz+1)) - sum(B_x(ix:ix+1,iy,iz)))/del(3)

  curl_y_scale = (abs(dx) + abs(dz)) / 2
  curl_y(ix,iy,iz) = (dx - dz) / 2

!

  dx = (sum(B_y(ix+1,iy:iy+1,iz)) - sum(B_y(ix,iy:iy+1,iz)))/del(1)
  dy = (sum(B_x(ix:ix+1,iy+1,iz)) - sum(B_x(ix:ix+1,iy,iz)))/del(2)

  curl_z_scale = (abs(dx) + abs(dy)) / 2
  curl_z(ix,iy,iz) = (dx - dy) / 2

!

  div_max = max(div_max, abs(div(ix,iy,iz)))
  div_scale_max = max(div_scale_max, abs(div_scale))

  curl_x_max = max(curl_x_max, abs(curl_x(ix,iy,iz)))
  curl_x_scale_max = max(curl_x_scale_max, abs(curl_x_scale))

  curl_y_max = max(curl_y_max, abs(curl_y(ix,iy,iz)))
  curl_y_scale_max = max(curl_y_scale_max, abs(curl_y_scale))

  curl_z_max = max(curl_z_max, abs(curl_z(ix,iy,iz)))
  curl_z_scale_max = max(curl_z_scale_max, abs(curl_z_scale))

enddo
enddo
enddo

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine print_stuff

implicit none

call db_rms_calc

print *, 'Nx, Ny, Nz:', Nx, Ny, Nz
print *, 'iz1, iz2:', iz1, iz2
print *, 'N_term:', n_term
print *, 'Chi2:     ', merit_tot
print *, 'Chi2_coef:', merit_coef
print *, 'B_diff (G):', 1e4*B_diff / (3 * n_pts)
print *, 'dB_rms (G):', dB_rms
print *, 'B_Merit   :', B_diff / sumB
print *, 'B_int:  ', B_int
print '(a, 4f10.2)', ' B_diff [x, y, z, tot]:', &
                  1e4*sum(abs(Bx_fit(:,:,:)-Bx(:,:,:)))/n_pts, &
                  1e4*sum(abs(By_fit(:,:,:)-By(:,:,:)))/n_pts, &
                  1e4*sum(abs(Bz_fit(:,:,:)-Bz(:,:,:)))/n_pts, &
                  1e4*B_diff/n_pts
print '(a, 4f10.2)', ' B_rms [x, y, z, tot]: ', &
                  1e4*sqrt(merit_x/n_pts), 1e4*sqrt(merit_y/n_pts), &
                  1e4*sqrt(merit_z/n_pts), dB_rms
print '(a, 4f10.2)', ' B_dat [x, y, z, tot]: ', &
                  1e4*sum(abs(Bx(:,:,:)))/n_pts, &
                  1e4*sum(abs(By(:,:,:)))/n_pts, &
                  1e4*sum(abs(Bz(:,:,:)))/n_pts, &
                  1e4*sumB/n_pts

end subroutine  

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine db_rms_calc

db_rms = 1e4 * sqrt(merit_data/n_pts) 

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine funcs_wf (a, deriv_calc)


implicit none

type (cartesian_map_term1_struct), pointer :: tm

real(rp), intent(in) :: a(:)

real(rp), dimension(3,n_var_term) :: dyda_wf
real(rp) kx_i, ky_i, kz_i
real(rp) dBx_i, dBy_i, dBz_i, f_c_i(3), f_kx_i(3), f_kz_i(3), f_phi_i(3)
real(rp) dBx_j, dBy_j, dBz_j, f_c_j(3), f_kx_j(3), f_kz_j(3), f_phi_j(3)
real(rp) df_kz_i(3), df_kz_j(3), dB(3)
real(rp) c_s_c, s_s_s, s_c_c, c_c_s, s_ds_c, c_dc_c, c_ds_s
real(rp) s_s_c, c_c_c, c_s_s

integer ix, iy, iz, i, j, jj, icoef, ikx, ikz, iphi, jcoef, jkx, jkz, jphi
integer ik, jk

logical, optional :: deriv_calc

! Compute the fit

Bx_fit = 0
By_fit = 0
Bz_fit = 0

do i = 1, n_term

  tm => bf(i)%tm

  jj = 1
  icoef = n_var_term*(i-1) + 1

  if (.not. mask_kx) then
    jj = jj + 1
    ikx   = n_var_term*(i-1) + jj
    tm%kx = a(ikx)
  endif

  if (.not. mask_ky) then
    jj = jj + 1
    iky   = n_var_term*(i-1) + jj
    tm%ky = a(iky)
  endif

  if (.not. mask_kz) then
    jj = jj + 1
    ikz   = n_var_term*(i-1) + jj
    tm%kz = a(ikz)
  endif

  if (.not. mask_x0) then
    jj = jj + 1
    ix0   = n_var_term*(i-1) + jj
    tm%x0 = a(ix0)
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    iy0   = n_var_term*(i-1) + jj
    tm%y0 = a(iy0)
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    iphi_z   = n_var_term*(i-1) + jj
    tm%phi_z = a(iphi_z)
  endif

  if (abs(tm%kx) > 1000 .or. abs(tm%ky) > 1000 .or. abs(tm%kz) > 1000) then
    print *, 'FUNCS_WF: |K| is > 1000', i, tm%kx, tm%ky, tm%kz
    merit_tot = 1e100_rp
    return
  endif


  if (bf(i)%kx > 0) then
    bf(i)%switch = hhyper_y$
    bf(i)%ky = sqrt(bf(i)%kz**2 + bf(i)%kx**2)
    bf(i)%dky_dkx =  bf(i)%kx
    bf(i)%dky_dkz =  bf(i)%kz
  elseif (bf(i)%kx > -abs(bf(i)%kz)) then
    bf(i)%switch = hhyper_xy$
    bf(i)%ky = sqrt(bf(i)%kz**2 - bf(i)%kx**2)
    bf(i)%dky_dkx = -bf(i)%kx
    bf(i)%dky_dkz =  bf(i)%kz
  else
    bf(i)%switch = hhyper_x$
    bf(i)%ky = sqrt(bf(i)%kx**2 - bf(i)%kz**2)
    bf(i)%dky_dkx =  bf(i)%kx
    bf(i)%dky_dkz = -bf(i)%kz
  endif

  kx_i = bf(i)%kx * del(1)
  ky_i = bf(i)%ky * del(2)
  kz_i = bf(i)%kz * del(3)

  if (bf(i)%switch == hhyper_x$ .or. bf(i)%switch == hhyper_xy$) then
    bf(i)%sign_x = +1
    bf(i)%sign_dc_x = +1
    bf(i)%coef_x = bf(i)%coef * bf(i)%kx 
    do ix = 1, Nx
      bf(i)%s_x(ix) =  sinh(kx_i * (ix-1)) 
      bf(i)%c_x(ix) =  cosh(kx_i * (ix-1)) 
    enddo
  else
    bf(i)%sign_x = -1
    bf(i)%sign_dc_x = -1
    bf(i)%coef_x = -bf(i)%coef * bf(i)%kx 
    do ix = 1, Nx
      bf(i)%s_x(ix) = sin(kx_i * (ix-1)) 
      bf(i)%c_x(ix) = cos(kx_i * (ix-1)) 
    enddo
  endif

  bf(i)%coef_y = bf(i)%coef
  if (bf(i)%switch == hhyper_y$ .or. bf(i)%switch == hhyper_xy$) then
    do iy = 1, Ny
      bf(i)%s_y(iy) = sinh_k (bf(i)%ky, del(2) * (iy-1)) 
      bf(i)%c_y(iy) = cosh (ky_i * (iy-1))
      bf(i)%ds_y(iy) = dsinh_k (bf(i)%ky, del(2) * (iy-1))
      bf(i)%dc_y(iy) = dcosh_k (bf(i)%ky, del(2) * (iy-1))
    enddo
  else
    do iy = 1, Ny
      bf(i)%s_y(iy) = sin_k (bf(i)%ky, del(2) * (iy-1)) 
      bf(i)%c_y(iy) = cos (ky_i * (iy-1)) 
      bf(i)%ds_y(iy) = dsin_k (bf(i)%ky, del(2) * (iy-1))
      bf(i)%dc_y(iy) = dcos_k (bf(i)%ky, del(2) * (iy-1))
    enddo
  endif

  bf(i)%coef_z = -bf(i)%coef * bf(i)%kz 
  do iz = iz1, iz2 ! 1, Nz
    bf(i)%s_z(iz) = sin(kz_i * (iz-Nz_ref) + bf(i)%phi_z) 
    bf(i)%c_z(iz) = cos(kz_i * (iz-Nz_ref) + bf(i)%phi_z) 
  enddo

  do ix = 1, Nx
  do iy = 1, Ny
  do iz = iz1, iz2 ! 1, Nz

    Bx_fit(ix,iy,iz) = Bx_fit(ix,iy,iz) + bf(i)%coef_x * &
                                bf(i)%s_x(ix) * bf(i)%s_y(iy) * bf(i)%c_z(iz)
    By_fit(ix,iy,iz) = By_fit(ix,iy,iz) + bf(i)%coef_y * &
                                bf(i)%c_x(ix) * bf(i)%c_y(iy) * bf(i)%c_z(iz)
    Bz_fit(ix,iy,iz) = Bz_fit(ix,iy,iz) + bf(i)%coef_z * &
                                bf(i)%c_x(ix) * bf(i)%s_y(iy) * bf(i)%s_z(iz)

  enddo
  enddo
  enddo

enddo

!

B_int = sum (By_fit(1, 1, :)) / Nz

B_diff = sum(abs(Bx_fit(:,:,:)-Bx(:,:,:))) + &
         sum(abs(By_fit(:,:,:)-By(:,:,:))) + &
         sum(abs(Bz_fit(:,:,:)-Bz(:,:,:))) 

merit_coef = sum(abs((term(1:n_term)%coef) * weight_coef(1:n_term)))

merit_x = sum((Bx_fit(:,:,:)-Bx(:,:,:))**2) 
merit_y = sum((By_fit(:,:,:)-By(:,:,:))**2) 
merit_z = sum((Bz_fit(:,:,:)-Bz(:,:,:))**2) 

merit_data = merit_x + merit_y + merit_z

merit_tot = merit_data + merit_coef

! Compute derivitives

if (present(deriv_calc)) then
  if (.not. deriv_calc) return
endif

if (.not. dyda_calc) return

beta_wf = 0

do i = 1, n_term

  icoef = n_var_term*(i-1) + 1
  ikx   = n_var_term*(i-1) + 2
  ikz   = n_var_term*(i-1) + 3
  if (vary_phi_z) iphi  = n_var_term*(i-1) + 4

  do ix = 1, Nx
  do iy = 1, Ny
  do iz = iz1, iz2 ! 1, Nz

    s_s_c = bf(i)%s_x(ix) * bf(i)%s_y(iy) * bf(i)%c_z(iz)
    c_c_c = bf(i)%c_x(ix) * bf(i)%c_y(iy) * bf(i)%c_z(iz)
    c_s_s = bf(i)%c_x(ix) * bf(i)%s_y(iy) * bf(i)%s_z(iz) 

    c_s_c = bf(i)%c_x(ix) * bf(i)%s_y(iy) * bf(i)%c_z(iz) 
    s_s_s = bf(i)%s_x(ix) * bf(i)%s_y(iy) * bf(i)%s_z(iz)
    s_c_c = bf(i)%s_x(ix) * bf(i)%c_y(iy) * bf(i)%c_z(iz) 
    c_c_s = bf(i)%c_x(ix) * bf(i)%c_y(iy) * bf(i)%s_z(iz)
    s_ds_c = bf(i)%s_x(ix) * bf(i)%ds_y(iy) * bf(i)%c_z(iz) 
    c_dc_c = bf(i)%c_x(ix) * bf(i)%dc_y(iy) * bf(i)%c_z(iz)
    c_ds_s = bf(i)%c_x(ix) * bf(i)%ds_y(iy) * bf(i)%s_z(iz)

    dyda_wf(1, 1) =  bf(i)%sign_x * bf(i)%kx * s_s_c
    dyda_wf(2, 1) =  c_c_c
    dyda_wf(3, 1) = -bf(i)%kz * c_s_s

    dyda_wf(1, 2) =  bf(i)%sign_x * bf(i)%coef * s_s_c + &
                       bf(i)%coef_x * (ix-1)*del(1) * c_s_c + &
                       bf(i)%coef_x *bf(i)%dky_dkx * s_ds_c
    dyda_wf(2, 2) =  bf(i)%coef_y * bf(i)%sign_dc_x * (ix-1)*del(1) * s_c_c + &
                          bf(i)%coef_y *bf(i)%dky_dkx * c_dc_c
    dyda_wf(3, 2) =  bf(i)%coef_z * bf(i)%sign_dc_x * (ix-1)*del(1) * s_s_s + &
                          bf(i)%coef_z *bf(i)%dky_dkx * c_ds_s

    dyda_wf(1, 3) = -bf(i)%coef_x * (iz-Nz_ref)*del(3) * s_s_s + &
                          bf(i)%coef_x *bf(i)%dky_dkz * s_ds_c
    dyda_wf(2, 3) = -bf(i)%coef_y * (iz-Nz_ref)*del(3) * c_c_s + &
                          bf(i)%coef_y *bf(i)%dky_dkz * c_dc_c
    dyda_wf(3, 3) = -bf(i)%coef * c_s_s + &
                          bf(i)%coef_z * (iz-Nz_ref)*del(3) * c_s_c + &
                          bf(i)%coef_z *bf(i)%dky_dkz * c_ds_s

    if (vary_phi_z) then
      dyda_wf(1, 4) = -bf(i)%coef_x * s_s_s
      dyda_wf(2, 4) = -bf(i)%coef_y * c_c_s
      dyda_wf(3, 4) =  bf(i)%coef_z * c_s_c
    endif

    dB(1) = Bx_fit(ix,iy,iz) - Bx(ix,iy,iz)
    dB(2) = By_fit(ix,iy,iz) - By(ix,iy,iz)
    dB(3) = Bz_fit(ix,iy,iz) - Bz(ix,iy,iz)

    do j = 1, n_var_term
      jj = n_var_term*(i-1) + j
      beta_wf(jj) = beta_wf(jj) - dot_product(dyda_wf(:,j), dB)
    enddo

  enddo
  enddo
  enddo

  beta_wf(icoef) = beta_wf(icoef) - sign(weight_coef(i), term(i)%coef) / 2

enddo

! alpha_wf calc

alpha_wf = 0

do i = 1, n_term

  icoef = n_var_term*(i-1) + 1
  ikx   = n_var_term*(i-1) + 2
  ikz   = n_var_term*(i-1) + 3
  if (vary_phi_z) iphi  = n_var_term*(i-1) + 4

  alpha_wf(icoef, icoef) = 0           ! was: weight_coef(i)

  do j = 1, n_term
    jcoef = n_var_term*(j-1) + 1
    jkx   = n_var_term*(j-1) + 2
    jkz   = n_var_term*(j-1) + 3
    if (vary_phi_z) jphi  = n_var_term*(j-1) + 4

    if (j < i) then
      ik = n_var_term*(i-1) + n_var_term
      jk = n_var_term*(j-1) + n_var_term
      alpha_wf(icoef:ik,jcoef:jk) = transpose(alpha_wf(jcoef:jk, icoef:ik))
      cycle
    endif

    do ix = 1, Nx
      do iy = 1, Ny

        dBx_i = bf(i)%sign_x * bf(i)%s_x(ix) * bf(i)%s_y(iy)
        dBy_i =                bf(i)%c_x(ix) * bf(i)%c_y(iy)
        dBz_i =                bf(i)%c_x(ix) * bf(i)%s_y(iy)

        f_c_i(1)  =  dBx_i * bf(i)%kx
        f_c_i(2)  =  dBy_i
        f_c_i(3)  = -dBz_i * bf(i)%kz

        f_kx_i(1) =  dBx_i * bf(i)%coef + &
               bf(i)%coef_x * (ix-1)*del(1) * bf(i)%c_x(ix) * bf(i)%s_y(iy) + &
               bf(i)%coef_x * bf(i)%dky_dkx * bf(i)%s_x(ix) * bf(i)%ds_y(iy) 
        f_kx_i(2) =  bf(i)%coef_y * bf(i)%sign_dc_x * (ix-1)*del(1) * bf(i)%s_x(ix) * bf(i)%c_y(iy) + &
               bf(i)%coef_y * bf(i)%dky_dkx * bf(i)%c_x(ix) * bf(i)%dc_y(iy) 
        f_kx_i(3) =  bf(i)%coef_z * bf(i)%sign_dc_x * (ix-1)*del(1) * bf(i)%s_x(ix) * bf(i)%s_y(iy) + &
               bf(i)%coef_z * bf(i)%dky_dkx * bf(i)%c_x(ix) * bf(i)%ds_y(iy) 

        f_kz_i(1) = bf(i)%coef_x * bf(i)%dky_dkz * bf(i)%s_x(ix) * bf(i)%ds_y(iy) 
        f_kz_i(2) = bf(i)%coef_y * bf(i)%dky_dkz * bf(i)%c_x(ix) * bf(i)%dc_y(iy) 
        f_kz_i(3) = bf(i)%coef_z * bf(i)%dky_dkz * bf(i)%c_x(ix) * bf(i)%ds_y(iy) - &
               bf(i)%coef * dBz_i 

        df_kz_i(1) =  bf(i)%coef_x * bf(i)%s_x(ix) * bf(i)%s_y(iy)
        df_kz_i(2) =  bf(i)%coef_y * bf(i)%c_x(ix) * bf(i)%c_y(iy) 
        df_kz_i(3) =  bf(i)%coef_z * bf(i)%c_x(ix) * bf(i)%s_y(iy) 

        f_phi_i(1) = -bf(i)%coef_x * bf(i)%s_x(ix) * bf(i)%s_y(iy)
        f_phi_i(2) = -bf(i)%coef_y * bf(i)%c_x(ix) * bf(i)%c_y(iy)
        f_phi_i(3) =  bf(i)%coef_z * bf(i)%c_x(ix) * bf(i)%s_y(iy)


        dBx_j = bf(j)%sign_x * bf(j)%s_x(ix) * bf(j)%s_y(iy)
        dBy_j =                bf(j)%c_x(ix) * bf(j)%c_y(iy)
        dBz_j =                bf(j)%c_x(ix) * bf(j)%s_y(iy)

        f_c_j(1)  =  dBx_j * bf(j)%kx
        f_c_j(2)  =  dBy_j
        f_c_j(3)  = -dBz_j * bf(j)%kz

        f_kx_j(1) =  dBx_j * bf(j)%coef + &
               bf(j)%coef_x * (ix-1)*del(1) * bf(j)%c_x(ix) * bf(j)%s_y(iy) + &
               bf(j)%coef_x * bf(j)%dky_dkx * bf(j)%s_x(ix) * bf(j)%ds_y(iy) 
        f_kx_j(2) =  bf(j)%coef_y * bf(j)%sign_dc_x * (ix-1)*del(1) * bf(j)%s_x(ix) * bf(j)%c_y(iy) + &
               bf(j)%coef_y * bf(j)%dky_dkx * bf(j)%c_x(ix) * bf(j)%dc_y(iy) 
        f_kx_j(3) =  bf(j)%coef_z * bf(j)%sign_dc_x * (ix-1)*del(1) * bf(j)%s_x(ix) * bf(j)%s_y(iy) + &
               bf(j)%coef_z * bf(j)%dky_dkx * bf(j)%c_x(ix) * bf(j)%ds_y(iy) 

        f_kz_j(1) = bf(j)%coef_x * bf(j)%dky_dkz * bf(j)%s_x(ix) * bf(j)%ds_y(iy) 
        f_kz_j(2) = bf(j)%coef_y * bf(j)%dky_dkz * bf(j)%c_x(ix) * bf(j)%dc_y(iy) 
        f_kz_j(3) = bf(j)%coef_z * bf(j)%dky_dkz * bf(j)%c_x(ix) * bf(j)%ds_y(iy) - bf(j)%coef * dBz_j 

        df_kz_j(1) =  bf(j)%coef_x * bf(j)%s_x(ix) * bf(j)%s_y(iy)
        df_kz_j(2) =  bf(j)%coef_y * bf(j)%c_x(ix) * bf(j)%c_y(iy) 
        df_kz_j(3) =  bf(j)%coef_z * bf(j)%c_x(ix) * bf(j)%s_y(iy) 
        

        alpha_wf(icoef, jcoef) = alpha_wf(icoef, jcoef) +  sum_z(i, j, f_c_i, f_c_j,   0, 0)
        alpha_wf(icoef, jkx)   = alpha_wf(icoef, jkx)   +  sum_z(i, j, f_c_i, f_kx_j,  0, 0)
        alpha_wf(icoef, jkz)   = alpha_wf(icoef, jkz)   +  sum_z(i, j, f_c_i, f_kz_j,  0, 0)
        alpha_wf(icoef, jkz)   = alpha_wf(icoef, jkz)   +  sum_z(i, j, f_c_i, df_kz_j, 2, 0)

        alpha_wf(ikx, jcoef) = alpha_wf(ikx, jcoef) +  sum_z(i, j, f_kx_i, f_c_j,   0, 0)
        alpha_wf(ikx, jkx)   = alpha_wf(ikx, jkx)   +  sum_z(i, j, f_kx_i, f_kx_j,  0, 0)
        alpha_wf(ikx, jkz)   = alpha_wf(ikx, jkz)   +  sum_z(i, j, f_kx_i, f_kz_j,  0, 0)
        alpha_wf(ikx, jkz)   = alpha_wf(ikx, jkz)   +  sum_z(i, j, f_kx_i, df_kz_j, 2, 0)

        alpha_wf(ikz, jcoef) = alpha_wf(ikz, jcoef) +  sum_z(i, j,  f_kz_i,  f_c_j,  0, 0)
        alpha_wf(ikz, jcoef) = alpha_wf(ikz, jcoef) +  sum_z(i, j, df_kz_i,  f_c_j,  1, 0)
        alpha_wf(ikz, jkx)   = alpha_wf(ikz, jkx)   +  sum_z(i, j,  f_kz_i,  f_kx_j, 0, 0)
        alpha_wf(ikz, jkx)   = alpha_wf(ikz, jkx)   +  sum_z(i, j, df_kz_i,  f_kx_j, 1, 0)
        alpha_wf(ikz, jkz)   = alpha_wf(ikz, jkz)   +  sum_z(i, j,  f_kz_i,  f_kz_j, 0, 0)
        alpha_wf(ikz, jkz)   = alpha_wf(ikz, jkz)   +  sum_z(i, j, df_kz_i,  f_kz_j, 1, 0)
        alpha_wf(ikz, jkz)   = alpha_wf(ikz, jkz)   +  sum_z(i, j,  f_kz_i, df_kz_j, 2, 0)
        alpha_wf(ikz, jkz)   = alpha_wf(ikz, jkz)   +  sum_z(i, j, df_kz_i, df_kz_j, 3, 0)

        if (vary_phi_z) then
          f_phi_j(1) = -bf(j)%coef_x * bf(j)%s_x(ix) * bf(j)%s_y(iy)
          f_phi_j(2) = -bf(j)%coef_y * bf(j)%c_x(ix) * bf(j)%c_y(iy)
          f_phi_j(3) =  bf(j)%coef_z * bf(j)%c_x(ix) * bf(j)%s_y(iy)

          alpha_wf(iphi, jcoef) = alpha_wf(iphi, jcoef) +  sum_z(i, j, f_phi_i, f_c_j,   0, 1)
          alpha_wf(iphi, jkx)   = alpha_wf(iphi, jkx)   +  sum_z(i, j, f_phi_i, f_kx_j,  0, 1)
          alpha_wf(iphi, jkz)   = alpha_wf(iphi, jkz)   +  sum_z(i, j, f_phi_i, f_kz_j,  0, 1)
          alpha_wf(iphi, jkz)   = alpha_wf(iphi, jkz)   +  sum_z(i, j, f_phi_i, df_kz_j, 2, 1)

          alpha_wf(icoef, jphi) = alpha_wf(icoef, jphi) +  sum_z(i, j, f_c_i,   f_phi_j, 0, 2)
          alpha_wf(ikx, jphi)   = alpha_wf(ikx,   jphi) +  sum_z(i, j, f_kx_i,  f_phi_j, 0, 2)
          alpha_wf(ikz, jphi)   = alpha_wf(ikz,   jphi) +  sum_z(i, j,  f_kz_i, f_phi_j, 0, 2)
          alpha_wf(ikz, jphi)   = alpha_wf(ikz,   jphi) +  sum_z(i, j, df_kz_i, f_phi_j, 1, 2)
          alpha_wf(iphi, jphi)  = alpha_wf(iphi,  jphi) +  sum_z(i, j, f_phi_i, f_phi_j, 0, 3)
        endif

      enddo
    enddo

  enddo
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function sum_z (ii, jj, f_i, f_j, n_deriv, n_phi) result (sumz)

real(rp) f_i(3), f_j(3), sumz
real(rp) a_cc(3), a_ss(3), a_sum, a_dif
real(rp) k_ave, k_dif, d

integer ii, jj, n_deriv, n_phi


! 

select case (n_phi)

case(0)
  a_cc(1:2) = cos(bf(ii)%phi_z) * cos(bf(jj)%phi_z)
  a_ss(1:2) = sin(bf(ii)%phi_z) * sin(bf(jj)%phi_z)
  a_cc(3) = a_ss(1)
  a_ss(3) = a_cc(1)

case(1)
  a_cc(1:2) =  sin(bf(ii)%phi_z) * cos(bf(jj)%phi_z)
  a_ss(1:2) = -cos(bf(ii)%phi_z) * sin(bf(jj)%phi_z)
  a_cc(3) = -a_ss(1)
  a_ss(3) = -a_cc(1)
  
case(2)
  a_cc(1:2) =  cos(bf(ii)%phi_z) * sin(bf(jj)%phi_z)
  a_ss(1:2) = -sin(bf(ii)%phi_z) * cos(bf(jj)%phi_z)
  a_cc(3) = -a_ss(1)
  a_ss(3) = -a_cc(1)

case(3)
  a_cc(1:2) =  sin(bf(ii)%phi_z) * sin(bf(jj)%phi_z)
  a_ss(1:2) =  cos(bf(ii)%phi_z) * cos(bf(jj)%phi_z)
  a_cc(3) =  a_ss(1)
  a_ss(3) =  a_cc(1)

case default
  print *, 'SUM_Z: PHI ERROR!', n_phi
  call err_exit

end select

a_sum = sum(f_i*f_j*a_cc) + sum(f_i*f_j*a_ss)
a_dif = sum(f_i*f_j*a_cc) - sum(f_i*f_j*a_ss)

!

d = del(3)
k_ave = d * (bf(jj)%kz + bf(ii)%kz) / 2
k_dif = d * (bf(jj)%kz - bf(ii)%kz) / 2

! no derivatives

select case (n_deriv)

case(0)

  if (abs(k_dif) < 1e-6) then   ! equal case
    if (abs(k_ave) < 1e-10) then
      sumz = a_sum * Nz/2 + a_dif * Nz/2
    else
      sumz = a_sum * Nz/2 + a_dif * Csc(k_ave)*Sin(k_ave*Nz)/2 
    endif
  else
    sumz = a_sum * Csc(k_dif)*Sin(k_dif*Nz)/2 + a_dif * Csc(k_ave)*Sin(k_ave*Nz)/2
  endif

! single derivative.

case(1, 2)

  if (n_deriv == 2) k_dif = -k_dif

  if (abs(k_dif) < 1e-6) then   ! equal case
    if (abs(k_ave) < 1e-10) then
      sumz = a_dif * d*k_ave*(1 - Nz**2) / 12
    else
      sumz = a_dif * d*(Nz*Cos(k_ave*Nz)*Csc(k_ave) - &
                       Cot(k_ave)*Csc(k_ave)*Sin(k_ave*Nz))/4 
    endif
  else
    sumz = a_sum * d*(-Nz*Cos(k_dif*Nz)*Csc(k_dif) + Cot(k_dif)*Csc(k_dif)*Sin(k_dif*Nz))/ 4 + &
           a_dif * d*(Nz*Cos(k_ave*Nz)*Csc(k_ave) - Cot(k_ave)*Csc(k_ave)*Sin(k_ave*Nz)) / 4 
  endif

!

case(3)

  if (abs(k_dif) < 1e-6) then   ! equal case
    if (abs(k_ave) < 1e-10) then
      sumz = 0
    else
      sumz = a_sum * d**2 * (Nz**3 - Nz) / 24 + &
           a_dif * d**2*(-2*Nz*Cos(k_ave*Nz)*Cot(k_ave)*Csc(k_ave) - &
              Nz**2*Csc(k_ave)*Sin(k_ave*Nz) + Cot(k_ave)**2*Csc(k_ave)*Sin(k_ave*Nz) + &
              Csc(k_ave)**3*Sin(k_ave*Nz)) / 8
    endif
  else
    sumz = a_sum * d**2*(2*Nz*Cos(k_dif*Nz)*Cot(k_dif)*Csc(k_dif) + &
             Nz**2*Csc(k_dif)*Sin(k_dif*Nz) - Cot(k_dif)**2*Csc(k_dif)*Sin(k_dif*Nz) - &
             Csc(k_dif)**3*Sin(k_dif*Nz)) / 8 + &
           a_dif * d**2*(-2*Nz*Cos(k_ave*Nz)*Cot(k_ave)*Csc(k_ave) - &
             Nz**2*Csc(k_ave)*Sin(k_ave*Nz) + Cot(k_ave)**2*Csc(k_ave)*Sin(k_ave*Nz) + &
             Csc(k_ave)**3*Sin(k_ave*Nz)) / 8
  endif

!

case default

  print *, 'Error in sum_z: bad n_deriv', n_deriv
  call err_exit

end select

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function csc(x) result (csc_trig)
implicit none
real(rp) x, csc_trig
csc_trig = 1 / sin(x)
end function

function cot(x) result (cot_trig)
implicit none
real(rp) x, cot_trig
cot_trig = cos(x) / sin(x)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function sinh_k (k, y)

implicit none

real (rp) sinh_k, k, y, arg

! sinh(k*y) / k

arg = k * y

if (abs(arg) < 1e-8) then
  sinh_k = y + arg**2 * y / 6
else
  sinh_k = sinh(arg) / k 
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function sin_k (k, y)

implicit none

real (rp) sin_k, k, y, arg

! sin(k*y) / k

arg = k * y

if (abs(arg) < 1e-8) then
  sin_k = y - arg**2 * y / 6
else
  sin_k = sin(arg) / k 
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function dsinh_k (k, y)

implicit none

real (rp) dsinh_k, k, y, arg

! d/dk[sinh(k*y) / k] / k

arg = k * y

if (abs(arg) < 1e-8) then
  dsinh_k = (1 + arg**2/10) * y**3 / 3 
else
  dsinh_k = (arg*cosh(arg) - sinh(arg)) / k**3 
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function dsin_k (k, y)

implicit none

real (rp) dsin_k, k, y, arg

! d/dk[sin(k*y) / k] / k

arg = k * y

if (abs(arg) < 1e-8) then
  dsin_k = (-1 + arg**2/10) * y**3 / 3 
else
  dsin_k = (arg*cos(arg) - sin(arg)) / k**3 
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function dcosh_k (k, y)

implicit none

real (rp) dcosh_k, k, y, arg

! d/dk[cosh(k*y)] / k

arg = k * y

if (abs(arg) < 1e-8) then
  dcosh_k = y**2 * (1 + arg**2/6)
else
  dcosh_k = y * sinh(arg) / k
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function dcos_k (k, y)

implicit none

real (rp) dcos_k, k, y, arg

! d/dk[cos(k*y)] / k

arg = k * y

if (abs(arg) < 1e-8) then
  dcos_k = y**2 * (1 - arg**2/6)
else
  dcos_k = y * sin(arg) / k
endif

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

	SUBROUTINE mrqmin_wf(a,maska,covar,alpha,chisq,alamda)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagmult
	USE nr, ONLY : covsrt,gaussj
	IMPLICIT NONE
	REAL(RP), DIMENSION(:), INTENT(INOUT) :: a
	REAL(RP), DIMENSION(:,:), INTENT(OUT) :: covar, alpha
	REAL(RP), INTENT(OUT) :: chisq
	REAL(RP), INTENT(INOUT) :: alamda
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTEGER(I4B) :: ma,ndata
	INTEGER(I4B), SAVE :: mfit
  integer, save :: a_count = 0
	call mrqmin_private_wf

!-------------------------------------------
	CONTAINS

	SUBROUTINE mrqmin_private_wf
	REAL(RP), SAVE :: ochisq
	REAL(RP), DIMENSION(:), ALLOCATABLE, SAVE :: atry, beta
	REAL(RP), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
	ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
		size(alpha,1),size(alpha,2)/),'mrqmin: ma')

	mfit=count(maska)

	if (alamda < 0.0) then
		allocate(atry(ma),beta(ma),da(ma,1))
		alamda=1.0_rp
		call mrqcof_wf(a, alpha, beta)
		ochisq=chisq
		atry=a
	end if

	covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
	call diagmult(covar(1:mfit,1:mfit),1.0_rp+alamda)
	da(1:mfit,1)=beta(1:mfit)
	call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))

	if (alamda == 0.0) then
		call covsrt(covar,maska)
		call covsrt(alpha_wf,maska)
		deallocate(atry,beta,da)
		RETURN
	end if

	atry=a+unpack(da(1:mfit,1),maska,0.0_rp)
	call mrqcof_wf(atry,covar,da(1:mfit,1))

	if (chisq < ochisq) then
    if (a_count >= 3) then
  		alamda=0.1_rp*alamda
      a_count = 0
    else
      a_count = a_count + 1
    endif
		ochisq=chisq
		alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
		beta(1:mfit)=da(1:mfit,1)
		a=atry
	else
		alamda=10.0_rp*alamda
    a_count = 0
		chisq=ochisq
	end if
	END SUBROUTINE mrqmin_private_wf

!----------------------------------------

	SUBROUTINE mrqcof_wf(a, alpha, beta)
	REAL(RP), DIMENSION(:), INTENT(IN) :: a
  REAL(RP), DIMENSION(:), INTENT(OUT) :: beta
  REAL(RP), DIMENSION(:,:), INTENT(OUT) :: alpha
	INTEGER(I4B) :: j,k,l,m

	call funcs_wf(a)
  alpha = alpha_wf
  beta = beta_wf
  chisq = merit_tot

	END SUBROUTINE mrqcof_wf
	END SUBROUTINE mrqmin_wf

end module
