module wig_mod

use bmad

integer :: Nx_min, Nx_max, Ny_min, Ny_max, Nz_min, Nz_max, Nz_ref
integer :: n_pts, n_term, n_data, n_var, n_var_per_term
integer :: n_loops, n_cycles

type term_struct
  type (cartesian_map_term1_struct), cmt
  real(rp), pointer :: s_x(:), c_x(:), s_y(:), c_y(:), s_z(:), c_z(:)
end type     

type (term_struct) term(500)

real(rp), allocatable :: Bx_fit(:,:,:), By_fit(:,:,:), Bz_fit(:,:,:)
real(rp), allocatable :: Bx_in(:,:,:), By_in(:,:,:), Bz_in(:,:,:)
real(rp), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:) 

real(rp) :: del(3), sumB

real(rp), allocatable :: div(:,:,:), curl_x(:,:,:), curl_y(:,:,:), curl_z(:,:,:)
real(rp) div_max, div_scale_max, div_scale
real(rp) curl_x_max, curl_x_scale_max, curl_x_scale
real(rp) curl_y_max, curl_y_scale_max, curl_y_scale
real(rp) curl_z_max, curl_z_scale_max, curl_z_scale
real(rp) B_int, B_diff, merit_coef, merit_data, dB_rms
real(rp) :: merit_tot, weight_coef_all = 1.0, weight_coef(400)
real(rp) merit_x, merit_y, merit_z

logical :: dyda_calc = .true.
logical :: fft_write = .false.
logical :: mask_x0 = .true., mask_y0 = .true., mask_phi_z = .false.

character(80) :: field_file

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
  read (1) Nx_min, Nx_max, del(1)
  read (1) Ny_min, Ny_max, del(2)
  read (1) Nz_min, Nz_max, del(3)
  allocate (Bx_in(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), By_in(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), Bz_in(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))
  allocate (X(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), Y(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max), Z(Nx_min:Nx_max,Ny_min:Ny_max,Nz_min:Nz_max))

  do j = Ny_min, Ny_max
  do k = Nz_min, Nz_max
    read (1) Bx_in(:,j,k)
    read (1) By_in(:,j,k)
    read (1) Bz_in(:,j,k)
  end do
  enddo

  close (1)

else

  open (1, file = field_file, STATUS='OLD', readonly, shared)
  read (1, *) Nx_min, Nx_max, del(1)
  read (1, *) Ny_min, Ny_max, del(2)
  read (1, *) Nz_min, Nz_max, del(3)
  allocate (Bx_in(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max), By_in(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max), Bz_in(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max))
  allocate (X(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max), Y(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max), Z(Nx_min:Nx_max, Ny_min:Ny_max, Nz_min:Nz_max))
             
  do i = 1, size(columns)
    read (1, '(a)') line
    if (index(line, trim(columns(i))) == 0) then
      print *, 'ERROR: COLUMNS IN FIELD TABLE NOT IN CORRECT ORDER. EXPECTED: ', columns(i)
      print *, '       GOT: ', trim(line)
      call err_exit
    endif
  enddo

  do i = Nx_min, Nx_max
  do j = Ny_min, Ny_max
  do k = Nz_min, Nz_max
    read (1, *) X(i,j,k), Y(i,j,k), Z(i,j,k), Bx_in(i,j,k), By_in(i,j,k), Bz_in(i,j,k)
  end do
  end do
  end do

  close (1)

  do i = Nx_min, Nx_max
  do j = Ny_min, Ny_max
  do k = Nz_min, Nz_max
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

type term_input_struct
  real(rp) :: coef = 0
  real(rp) :: kx = 0, ky = 0, kz = 0
  real(rp) :: phi_z = 0
  character(2) :: family = 0        ! family_x$, etc.
  real(rp) :: x0 = 0, y0 = 0
end type

type (term_input_struct), cmt(500)

real(rp) ave_x, ave_z

integer :: i, j, k, ix, iy, iz

logical write_curl

! sorted_field_trimmed.txt

character(*) :: param_file

namelist / wght / weight_coef_all, n_loops, n_cycles, weight_coef, fft_write
namelist / parameters / field_file, Nz_ref, cmt

! Read in parameters and starting fit

weight_coef = -1
open (1, file = 'fit.weight', status = 'old', readonly, shared)
read (1, nml = wght)
close (1)
where (weight_coef < 0) weight_coef = weight_coef_all

coef = 0
open (1, file = param_file, status = 'old', readonly, shared)
read (1, nml = parameters)
close (1)

! Read in field data file

call read_field_table_file (field_file)

!-------------------------------------
! read in starting point
! first throw out any term with a zero coef.

n_term = 0
do i = 1, size(term)
  if (term(i)%coef == 0) cycle
  n_term = n_term + 1
  term(n_term)%cmt%coef  = term(i)%coef
  term(n_term)%cmt%kx    = term(i)%kx
  term(n_term)%cmt%ky    = term(i)%ky
  term(n_term)%cmt%kz    = term(i)%kz
  term(n_term)%cmt%phi_z = term(i)%phi_z
  term(n_term)%cmt%x0    = term(i)%x0
  term(n_term)%cmt%y0    = term(i)%y0
  do k = 1, ubound(cartesian_map_family_name, 1)
    if (upcase(term(i)%family) == cartesian_map_family_name(k)) then
      term(n_term)%cmt%family = k
      exit
    endif
  enddo
  if (k == ubound(cartesian_map_family_name, 1)+1) then
    print *, 'Bad family name for term: ', term(i)%family
    stop
  endif
enddo

n_var_per_term = 6  ! number of vars per term
if (mask_x0) n_var_per_term = n_var_per_term - 1
if (mask_y0) n_var_per_term = n_var_per_term - 1
if (mask_phi_z) n_var_per_term = n_var_per_term - 1

n_var = n_var_per_term * n_term
n_pts = (Nx_max-Nx_min+1) * (Ny_max-Ny_min+1) * (Nz_max-Nz_min+1)
n_data = 3*n_pts + n_term

do i = 1, size(term)
  allocate(term(i)%s_x(Nx_min:Nx_max), term(i)%c_x(Nx_min:Nx_max))
  allocate(term(i)%s_y(Ny_min:Ny_max), term(i)%c_y(Ny_min:Ny_max))
  allocate(term(i)%s_z(Nz_min:Nz_max), term(i)%c_z(Nz_min:Nz_max))
  term(i)%coef = coef
enddo

allocate (curl_x(Nx, Ny, Nz), curl_y(Nx, Ny, Nz), curl_z(Nx, Ny, Nz))
allocate (div(Nx, Ny, Nz))

allocate (Bx_fit(Nx, Ny, Nz), By_fit(Nx, Ny, Nz), Bz_fit(Nx, Ny, Nz))

! This is for seeing how well the curl of the field is zero.

If (write_curl) then
  call curl_calc (Bx_in, By_in, Bz_in)
endif

sumB = sum(abs(Bx_in) + abs(By_in) + abs(Bz_in))

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

do ix = Nx_min, Nx_max-1
do iy = Ny_min, Ny_max-1
do iz = Nz_min, Nz_max-1

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

end subroutine curl_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine print_stuff

implicit none

call db_rms_calc

print *, 'Nx_min, Nx_max:', Nx_min, Nx_max
print *, 'Ny_min, Ny_max:', Ny_min, Ny_max
print *, 'Nz_min, Nz_max:', Nz_min, Nz_max
print *, 'N_term:', n_term
print *, 'Chi2:     ', merit_tot
print *, 'Chi2_coef:', merit_coef
print *, 'B_diff (G):', B_diff / (3 * n_pts)
print *, 'dB_rms (G):', dB_rms
print *, 'B_Merit   :', B_diff / sumB
print *, 'B_int:  ', B_int
print '(a, 4f10.2)', ' B_diff [x, y, z, tot]:', &
                  sum(abs(Bx_fit(:,:,:)-Bx(:,:,:)))/n_pts, &
                  sum(abs(By_fit(:,:,:)-By(:,:,:)))/n_pts, &
                  sum(abs(Bz_fit(:,:,:)-Bz(:,:,:)))/n_pts, &
                  B_diff/n_pts
print '(a, 4f10.2)', ' B_rms [x, y, z, tot]: ', &
                  sqrt(merit_x/n_pts), sqrt(merit_y/n_pts), &
                  sqrt(merit_z/n_pts), dB_rms
print '(a, 4f10.2)', ' B_dat [x, y, z, tot]: ', &
                  sum(abs(Bx(:,:,:)))/n_pts, &
                  sum(abs(By(:,:,:)))/n_pts, &
                  sum(abs(Bz(:,:,:)))/n_pts, &
                  sumB/n_pts

end subroutine  

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine db_rms_calc

db_rms = sqrt(merit_data/n_pts) 

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine funcs_cf (a, yfit, dyda, status)

implicit none

type (cartesian_map_term1_struct), pointer :: cmt
type (term_struct), pointer :: tt

real(rp), intent(in) :: a(:)
real(rp), intent(out) :: yfit(:)
real(rp), intent(out) :: dyda(:, :)
real(rp) norm_x, norm_y, norm_z
real(rp) dnormx_dkx, dnormx_dky, dnormx_dkz, dnormy_dkx, dnormy_dky, dnormy_dkz, dnormz_dkx, dnormz_dky, dnormz_dkz

integer sgn_x, sgn_y, sgn_z, dsgn_x, dsgn_y, dsgn_z
integer status
integer ix, iy, iz, i, j, j_term, jj, jb

! Compute the fit

Bx_fit = 0
By_fit = 0
Bz_fit = 0

jj = 0
j_term = 0

do i = 1, n_term
  tt => term(i)
  cmt => tt%cmt
  jj = j_term

  jj = jj + 1
  cmt%coef = a(jj)

  jj = jj + 1
  cmt%kx = a(jj)

  jj = jj + 1
  cmt%kz = a(jj)

  if (.not. mask_x0) then
    jj = jj + 1
    cmt%x0 = a(jj)
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    cmt%y0 = a(jj)
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    cmt%phi_z = a(jj)
  endif

  if (abs(cmt%kx) > 1000 .or. abs(cmt%ky) > 1000 .or. abs(cmt%kz) > 1000) then
    print *, 'FUNCS_WF: |K| is > 1000', i, cmt%kx, cmt%ky, cmt%kz
    merit_tot = 1e100_rp
    return
  endif

  ! The forms are constructed to make the field a continuous function of (kx, ky, kz) when going
  ! between hyper_y and hyper_xy and between hyper_xy and hyper_x. Note that the field is *not*
  ! continuous between hyper_y and hyper_x.

  if (cmt%kx > 0) then
    cmt%form = hyper_y$
    cmt%ky = sqrt(cmt%kz**2 + cmt%kx**2)
  elseif (cmt%kx > -abs(cmt%kz)) then
    cmt%form = hyper_xy$
    cmt%ky = sqrt(cmt%kz**2 - cmt%kx**2)
  else
    cmt%form = hyper_x$
    cmt%ky = sqrt(cmt%kx**2 - cmt%kz**2)
  endif

  !

  do iz = Nz_min, Nz_max
    tt%s_z(iz) = sin(cmt%kz * del(3) * (iz-Nz_ref) + tt%phi_z)
    tt%c_z(iz) = cos(cmt%kz * del(3) * (iz-Nz_ref) + tt%phi_z)
  enddo

  sgn_x = 1; sgn_y = 1; sgn_z = 1

  select case (cmt%form)
  case (hyper_y$)
    do ix = Nx_min, Nx_max
      tt%s_x(ix) = sin(cmt%kx * (del(1) * ix + cmt%x0))
      tt%c_x(ix) = cos(cmt%kx * (del(1) * ix + cmt%x0))
    enddo
    do iy = Ny_min, Ny_max
      tt%s_y(iy) = sinh (cmt%ky * (del(2) * iy + cmt%y0))
      tt%c_y(iy) = cosh (cmt%ky * (del(2) * iy + cmt%y0))
    enddo

    select case (cmt%family)
    case (family_x$);  sgn_z = -1
    case (family_y$);  sgn_x = -1; sgn_z = -1
    case (family_qu$); sgn_z = -1
    case (family_sq$); sgn_x = -1; sgn_z = -1
    end select

    sgn_dky_dkx = 1; sgn_dky_dkz = 1
    norm_x = cmt%kx / cmt%ky
    norm_y = 1
    norm_z = cmt%kz / cmt%ky
    dnormx_dkx = 1 / cmt%ky; dnromx_dky = -norm_x / cmt%ky; dnormx_dkz = 0
    dnormy_dkx = 0; dnromy_dky = 0; dnormy_dkz = 0
    dnormz_dkx = 0; dnromz_dky = -norm_z / cmt%ky; dnormz_dkz = 1 / cmt%ky

    dsgn_x = -1; dsgn_y = 1; dsgn_z = -1

  case (hyper_xy$)
    do ix = Nx_min, Nx_max
      tt%s_x(ix) =  sinh(cmt%kx * (del(1) * ix + cmt%x0))
      tt%c_x(ix) =  cosh(cmt%kx * (del(1) * ix + cmt%x0))
    enddo
    do iy = Ny_min, Ny_max
      tt%s_y(iy) = sinh (cmt%ky * (del(2) * iy + cmt%y0))
      tt%c_y(iy) = cosh (cmt%ky * (del(2) * iy + cmt%y0))
    enddo

    norm = 1 / cmt%kz
    select case (cmt%family)
    case (family_x$);  sgn_z = -1
    case (family_y$);  sgn_z = -1
    case (family_qu$); sgn_z = -1
    case (family_sq$); sgn_z = -1
    end select

    sgn_dky_dkx = -1; sgn_dky_dkz = 1
    norm_x = cmt%kx / cmt%kz
    norm_y = cmt%ky / cmt%kz
    norm_z = 1
    dnormx_dkx = 1 / cmt%kz; dnromx_dky = 0; dnormx_dkz = -norm_x / cmt%kz
    dnormy_dkx = 0; dnromy_dky = 1 / cmt%kz; dnormy_dkz = -norm / cmt%kz
    dnormz_dkx = 0; dnromz_dky = 0; dnormz_dkz = 0

    dsgn_x = 1; dsgn_y = 1; dsgn_z = -1

  case (hyper_x$)
    do ix = Nx_min, Nx_max
      tt%s_x(ix) =  sinh(cmt%kx * (del(1) * ix + cmt%x0))
      tt%c_x(ix) =  cosh(cmt%kx * (del(1) * ix + cmt%x0))
    enddo
    do iy = Ny_min, Ny_max
      tt%s_y(iy) = sin (cmt%ky * (del(2) * iy + cmt%y0))
      tt%c_y(iy) = cos (cmt%ky * (del(2) * iy + cmt%y0))
    enddo

    norm = 1 / cmt%kx
    select case (cmt%family)
    case (family_x$);  sgn_y = -1; sgn_z = -1
    case (family_y$);  sgn_z = -1
    case (family_qu$); sgn_z = -1
    case (family_sq$); sgn_x = -1
    end select

    sgn_dky_dkx = 1; sgn_dky_dkz = -1
    norm_x = 1
    norm_y = cmt%ky / cmt%kx
    norm_z = cmt%kz / cmt%kx
    dnormx_dkx = 0; dnomx_dky = 0; dnormx_dkz = 0
    dnormy_dkx = -norm_y / cmt%kx; dnromy_dky = 1 / cmt%kx; dnormy_dkz = 0
    dnormz_dkx = -norm_z / cmt%kx; dnromz_dky = 0; dnormz_dkz = 1 / cmt%kx

    dsgn_x = 1; dsgn_y = -1; dsgn_z = -1
  end select

  !

  Bx_fit = 0
  By_fit = 0
  Bz_fit = 0
  jk = j_term
  jb = 0

  if (optimizer_running) then
    do ix = Nx_min, Nx_max, Nx_max-Nx_min
    do iy = Ny_min, Ny_max
    do iz = Nz_min, Nz_max
      call add_to_field()
    enddo
    enddo
    enddo

    do ix = Nx_min, Nx_max
    do iy = Ny_min, Ny_max, Ny_max-Ny_min
    do iz = Nz_min, Nz_max
      call add_to_field()
    enddo
    enddo
    enddo

    do ix = Nx_min, Nx_max
    do iy = Ny_min, Ny_max
    do iz = Nz_min, Nz_max, Nz_max-Nz_min
      call add_to_field()
    enddo
    enddo
    enddo

  else
    do ix = Nx_min, Nx_max
    do iy = Ny_min, Ny_max
    do iz = Nz_min, Nz_max
      call add_to_field()
    enddo
    enddo
    enddo
  endif

  j_term = jj

enddo

!

merit_x = 0
merit_y = 0
merit_z = 0

if (optimizer_running) then
  do ix = Nx_min, Nx_max, Nx_max-Nx_min
  do iy = Ny_min, Ny_max
  do iz = Nz_min, Nz_max
    jj = jj + 1
    yfit(jj) = Bx_fit(ix,iy,iz) - Bx_in(ix,iy,iz)
    merit_x = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = By_fit(ix,iy,iz) - By_in(ix,iy,iz)
    merit_y = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = Bz_fit(ix,iy,iz) - Bz_in(ix,iy,iz)
    merit_z = merit_x + yfit(jj)**2
  enddo
  enddo
  enddo

  do ix = Nx_min, Nx_max
  do iy = Ny_min, Ny_max, Ny_max-Ny_min
  do iz = Nz_min, Nz_max
    jj = jj + 1
    yfit(jj) = Bx_fit(ix,iy,iz) - Bx_in(ix,iy,iz)
    merit_x = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = By_fit(ix,iy,iz) - By_in(ix,iy,iz)
    merit_y = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = Bz_fit(ix,iy,iz) - Bz_in(ix,iy,iz)
    merit_z = merit_x + yfit(jj)**2
  enddo
  enddo
  enddo

  do ix = Nx_min, Nx_max
  do iy = Ny_min, Ny_max
  do iz = Nz_min, Nz_max, Nz_max-Nz_min
    jj = jj + 1
    yfit(jj) = Bx_fit(ix,iy,iz) - Bx_in(ix,iy,iz)
    merit_x = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = By_fit(ix,iy,iz) - By_in(ix,iy,iz)
    merit_y = merit_x + yfit(jj)**2
    jj = jj + 1
    yfit(jj) = Bz_fit(ix,iy,iz) - Bz_in(ix,iy,iz)
    merit_z = merit_x + yfit(jj)**2
  enddo
  enddo
  enddo

  return
endif

!

merit_coef = sum(abs((term(1:n_term)%coef) * weight_coef(1:n_term)))
merit_data = merit_x + merit_y + merit_z
merit_tot = merit_data + merit_coef

!---------------------------------------------
contains

subroutine add_to_field()

real(rp) 
integer jj

!

x_pos = del(1) * ix + cmt%x0
y_pos = del(2) * iy + cmt%y0
z_pos = del(3) * iz

!---

select case (cmt%family)
case (family_x$)
  Bx_fit(ix,iy,iz) = Bx_fit(ix,iy,iz) + cmt%coef * norm_x * c_x * c_y * c_z * sgn_x
  By_fit(ix,iy,iz) = By_fit(ix,iy,iz) + cmt%coef * norm_y * s_x * s_y * c_z * sgn_y
  Bz_fit(ix,iy,iz) = Bz_fit(ix,iy,iz) + cmt%coef * norm_z * s_x * c_y * s_z * sgn_z

  dyda(jb+1,j_term+1) = norm_x * c_x * c_y * c_z * sgn_x   ! dBx_fit/dCoef
  dyda(jb+2,j_term+1) = norm_y * s_x * s_y * c_z * sgn_y
  dyda(jb+3,j_term+1) = norm_z * s_x * c_y * s_z * sgn_z

  dyda(jb+1,j_term+2) = cmt%coef * c_y * c_z * sgn_x * (dnormx_dkx * c_x + norm_x * s_x * x_pos * dsgn_x)
  dyda(jb+2,j_term+2) = cmt%coef * s_y * c_z * sgn_y * (dnormy_dkx * s_x + norm_y * c_x * x_pos)
  dyda(jb+3,j_term+2) = cmt%coef * c_y * s_z * sgn_z * (dnormz_dkx * s_x + norm_z * c_x * x_pos)

  dyda(jb+1,j_term+3) = cmt%coef * c_x * c_y * sgn_x * (dnormx_dkz * c_z + norm_x * s_z * z_pos * dsgn_z)
  dyda(jb+2,j_term+3) = cmt%coef * s_x * s_y * sgn_y * (dnormy_dkz * c_z + norm_y * s_z * z_pos * dsgn_z)
  dyda(jb+3,j_term+3) = cmt%coef * s_x * c_y * sgn_z * (dnormz_dkz * s_z + norm_z * c_z * z_pos)

  jj = j_term+3
  if (.not. mask_x0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * s_x * c_y * c_z * sgn_x * cmt%kx * dsgn_x
    dyda(jb+2,jj) = cmt%coef * norm_y * c_x * s_y * c_z * sgn_y * cmt%kx 
    dyda(jb+3,jj) = cmt%coef * norm_z * c_x * c_y * s_z * sgn_z * cmt%kx 
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * s_y * c_z * sgn_x * cmt%ky * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * c_y * c_z * sgn_y * cmt%ky 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * s_y * s_z * sgn_z * cmt%ky * dsgn_y
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * c_y * s_z * sgn_x * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * s_y * s_z * sgn_y * dsgn_y 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * c_y * c_z * sgn_z
  endif

  ! Contribution from fact that ky is dependent upon kx and kx

  select case (cmt%form)
  case (hyper_y$)
    fx = cmt%coef * c_x * c_z * sgn_x * (dnormx_dky * c_y + norm_x * s_y * y_pos * dsgn_y) / cmt%ky
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y * y_pos) / cmt%ky
    fz = cmt%coef * s_x * s_z * sgn_z * (dnormz_dky * c_y + norm_z * s_y * y_pos * dsgn_y) / cmt%ky
  case (hyper_xy$)
    fx = cmt%coef * c_x * c_z * sgn_x * (norm_x * sinhx_x(cmt%ky*y_pos) * y_pos**2)
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y + sinhx_x(cmt%ky*y_pos)) y_pos / cmt%kz
    fz = cmt%coef * s_x * s_z * sgn_z * (norm_z * sinhx_x(cmt%ky*y_pos) * y_pos**2)
  case (hyper_x$)
    fx = -cmt%coef * c_x * c_z * sgn_x * (norm_x * sinx_x(cmt%ky*y_pos) * y_pos**2)
    fy =  cmt%coef * s_x * c_z * sgn_y * (c_y - sinx_x(cmt%ky*y_pos)) * y_pos / cmt%kx
    fz = -cmt%coef * s_x * s_z * sgn_z * (norm_z * sinx_x(cmt%ky*y_pos) * y_pos**2)
  end select

  dyda(jb+1,j_term+2) = dyda(jb+1,j_term+2) + fx * cmt%kx * sgn_dky_dkx
  dyda(jb+2,j_term+2) = dyda(jb+2,j_term+2) + fy * cmt%kx * sgn_dky_dkx
  dyda(jb+3,j_term+2) = dyda(jb+3,j_term+2) + fz * cmt%kx * sgn_dky_dkx
  dyda(jb+1,j_term+3) = dyda(jb+1,j_term+3) + fx * cmt%kz * sgn_dky_dkz
  dyda(jb+2,j_term+3) = dyda(jb+2,j_term+3) + fy * cmt%kz * sgn_dky_dkz
  dyda(jb+3,j_term+3) = dyda(jb+3,j_term+3) + fz * cmt%kz * sgn_dky_dkz


case (family_y$)
  Bx_fit(ix,iy,iz) = Bx_fit(ix,iy,iz) + cmt%coef * norm_x * s_x * s_y * c_z * sgn_x
  By_fit(ix,iy,iz) = By_fit(ix,iy,iz) + cmt%coef * norm_y * c_x * c_y * c_z * sgn_y
  Bz_fit(ix,iy,iz) = Bz_fit(ix,iy,iz) + cmt%coef * norm_z * c_x * s_y * s_z * sgn_z

  dyda(jb+1,j_term+1) = norm_x * s_x * s_y * c_z * sgn_x
  dyda(jb+2,j_term+1) = norm_y * c_x * c_y * c_z * sgn_y
  dyda(jb+3,j_term+1) = norm_z * c_x * s_y * s_z * sgn_z

  dyda(jb+1,j_term+2) = cmt%coef * s_y * c_z * sgn_x * (dnormx_dkx * s_x + norm_x * c_x * x_pos)
  dyda(jb+2,j_term+2) = cmt%coef * c_y * c_z * sgn_y * (dnormy_dkx * c_x + norm_y * s_x * x_pos * dsgn_x)
  dyda(jb+3,j_term+2) = cmt%coef * s_y * s_z * sgn_z * (dnormz_dkx * c_x + norm_z * s_x * x_pos * dsgn_x)

  dyda(jb+1,j_term+3) = cmt%coef * s_x * s_y * sgn_x * (dnormx_dkz * c_z + norm_x * s_z * z_pos * dsgn_z)
  dyda(jb+2,j_term+3) = cmt%coef * c_x * c_y * sgn_y * (dnormy_dkz * c_z + norm_y * s_z * z_pos * dsgn_z)
  dyda(jb+3,j_term+3) = cmt%coef * c_x * s_y * sgn_z * (dnormz_dkz * s_z + norm_z * c_z * z_pos)

  jj = j_term+3
  if (.not. mask_x0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * c_y * c_z * sgn_x * cmt%kx
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * c_y * c_z * sgn_y * cmt%kx * dsgn_x 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * s_y * s_z * sgn_z * cmt%kx * dsgn_x
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * s_x * c_y * c_z * sgn_x * cmt%ky * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * c_x * s_y * c_z * sgn_y * cmt%ky 
    dyda(jb+3,jj) = cmt%coef * norm_z * c_x * c_y * s_z * sgn_z * cmt%ky
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * s_x * s_y * s_z * sgn_x * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * c_x * c_y * s_z * sgn_y * dsgn_y 
    dyda(jb+3,jj) = cmt%coef * norm_z * c_x * s_y * c_z * sgn_z
  endif

  ! Contribution from fact that ky is dependent upon kx and kx

  select case (cmt%form)
  case (hyper_y$)
    fx = cmt%coef * c_x * c_z * sgn_x * (dnormx_dky * s_y + norm_x * c_y * y_pos) / cmt%ky
    fy = cmt%coef * s_x * c_z * sgn_y * (s_y * y_pos * dsgn_y) / cmt%ky
    fz = cmt%coef * s_x * s_z * sgn_z * (dnormz_dky * s_y + norm_z * c_y * y_pos) / cmt%ky
  case (hyper_xy$)
    fx = cmt%coef * c_x * c_z * sgn_x * (norm_x * sinhx_x(cmt%ky*y_pos) * y_pos**2)
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y + sinhx_x(cmt%ky*y_pos)) y_pos / cmt%kz
    fz = cmt%coef * s_x * s_z * sgn_z * (norm_z * sinhx_x(cmt%ky*y_pos) * y_pos**2)
  case (hyper_x$)
    fx = -cmt%coef * c_x * c_z * sgn_x * (norm_x * sinx_x(cmt%ky*y_pos) * y_pos**2)
    fy =  cmt%coef * s_x * c_z * sgn_y * (c_y - sinx_x(cmt%ky*y_pos)) * y_pos / cmt%kx
    fz = -cmt%coef * s_x * s_z * sgn_z * (norm_z * sinx_x(cmt%ky*y_pos) * y_pos**2)
  end select

  dyda(jb+1,j_term+2) = dyda(jb+1,j_term+2) + fx * cmt%kx * sgn_dky_dkx
  dyda(jb+2,j_term+2) = dyda(jb+2,j_term+2) + fy * cmt%kx * sgn_dky_dkx
  dyda(jb+3,j_term+2) = dyda(jb+3,j_term+2) + fz * cmt%kx * sgn_dky_dkx
  dyda(jb+1,j_term+3) = dyda(jb+1,j_term+3) + fx * cmt%kz * sgn_dky_dkz
  dyda(jb+2,j_term+3) = dyda(jb+2,j_term+3) + fy * cmt%kz * sgn_dky_dkz
  dyda(jb+3,j_term+3) = dyda(jb+3,j_term+3) + fz * cmt%kz * sgn_dky_dkz

!---

case (family_qu$)
  Bx_fit(ix,iy,iz) = Bx_fit(ix,iy,iz) + cmt%coef * norm * c_x * s_y * c_z * sgn_x
  By_fit(ix,iy,iz) = By_fit(ix,iy,iz) + cmt%coef * norm * s_x * c_y * c_z * sgn_y
  Bz_fit(ix,iy,iz) = Bz_fit(ix,iy,iz) + cmt%coef * norm * s_x * s_y * s_z * sgn_z

  dyda(jb+1,j_term+1) = norm_x * c_x * s_y * c_z * sgn_x
  dyda(jb+2,j_term+1) = norm_y * s_x * c_y * c_z * sgn_y
  dyda(jb+3,j_term+1) = norm_z * s_x * s_y * s_z * sgn_z

  dyda(jb+1,j_term+2) = cmt%coef * c_y * c_z * sgn_x * (dnormx_dkx * c_x + norm_x * s_x * x_pos * dsgn_x)
  dyda(jb+2,j_term+2) = cmt%coef * s_y * c_z * sgn_y * (dnormy_dkx * s_x + norm_y * c_x * x_pos)
  dyda(jb+3,j_term+2) = cmt%coef * c_y * s_z * sgn_z * (dnormz_dkx * s_x + norm_z * c_x * x_pos)

  dyda(jb+1,j_term+3) = cmt%coef * c_x * c_y * sgn_x * (dnormx_dkz * c_z + norm_x * s_z * z_pos * dsgn_z)
  dyda(jb+2,j_term+3) = cmt%coef * s_x * s_y * sgn_y * (dnormy_dkz * c_z + norm_y * s_z * z_pos * dsgn_z)
  dyda(jb+3,j_term+3) = cmt%coef * s_x * c_y * sgn_z * (dnormz_dkz * s_z + norm_z * c_z * z_pos)

  jj = j_term+3
  if (.not. mask_x0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * s_x * c_y * c_z * sgn_x * cmt%kx * dsgn_x
    dyda(jb+2,jj) = cmt%coef * norm_y * c_x * s_y * c_z * sgn_y * cmt%kx 
    dyda(jb+3,jj) = cmt%coef * norm_z * c_x * c_y * s_z * sgn_z * cmt%kx 
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * s_y * c_z * sgn_x * cmt%ky * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * c_y * c_z * sgn_y * cmt%ky 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * s_y * s_z * sgn_z * cmt%ky * dsgn_y
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * c_y * s_z * sgn_x * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * s_y * s_z * sgn_y * dsgn_y 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * c_y * c_z * sgn_z
  endif

  ! Contribution from fact that ky is dependent upon kx and kx

  select case (cmt%form)
  case (hyper_y$)
    fx = cmt%coef * c_x * c_z * sgn_x * (dnormx_dky * c_y + norm_x * s_y * y_pos * dsgn_y) / cmt%ky
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y * y_pos) / cmt%ky
    fz = cmt%coef * s_x * s_z * sgn_z * (dnormz_dky * c_y + norm_z * s_y * y_pos * dsgn_y) / cmt%ky
  case (hyper_xy$)
    fx = cmt%coef * c_x * c_z * sgn_x * (norm_x * sinhx_x(cmt%ky*y_pos) * y_pos**2)
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y + sinhx_x(cmt%ky*y_pos)) y_pos / cmt%kz
    fz = cmt%coef * s_x * s_z * sgn_z * (norm_z * sinhx_x(cmt%ky*y_pos) * y_pos**2)
  case (hyper_x$)
    fx = -cmt%coef * c_x * c_z * sgn_x * (norm_x * sinx_x(cmt%ky*y_pos) * y_pos**2)
    fy =  cmt%coef * s_x * c_z * sgn_y * (c_y - sinx_x(cmt%ky*y_pos)) * y_pos / cmt%kx
    fz = -cmt%coef * s_x * s_z * sgn_z * (norm_z * sinx_x(cmt%ky*y_pos) * y_pos**2)
  end select

  dyda(jb+1,j_term+2) = dyda(jb+1,j_term+2) + fx * cmt%kx * sgn_dky_dkx
  dyda(jb+2,j_term+2) = dyda(jb+2,j_term+2) + fy * cmt%kx * sgn_dky_dkx
  dyda(jb+3,j_term+2) = dyda(jb+3,j_term+2) + fz * cmt%kx * sgn_dky_dkx
  dyda(jb+1,j_term+3) = dyda(jb+1,j_term+3) + fx * cmt%kz * sgn_dky_dkz
  dyda(jb+2,j_term+3) = dyda(jb+2,j_term+3) + fy * cmt%kz * sgn_dky_dkz
  dyda(jb+3,j_term+3) = dyda(jb+3,j_term+3) + fz * cmt%kz * sgn_dky_dkz

!---

case (family_sq$)
  Bx_fit(ix,iy,iz) = Bx_fit(ix,iy,iz) + cmt%coef * norm_x * s_x * c_y * c_z * sgn_x
  By_fit(ix,iy,iz) = By_fit(ix,iy,iz) + cmt%coef * norm_y * c_x * s_y * c_z * sgn_y
  Bz_fit(ix,iy,iz) = Bz_fit(ix,iy,iz) + cmt%coef * norm_z * c_x * c_y * s_z * sgn_z

  dyda(jb+1,j_term+1) = norm_x * s_x * c_y * c_z * sgn_x
  dyda(jb+2,j_term+1) = norm_y * c_x * s_y * c_z * sgn_y
  dyda(jb+3,j_term+1) = norm_z * c_x * c_y * s_z * sgn_z

  dyda(jb+1,j_term+2) = cmt%coef * c_y * c_z * sgn_x * (dnormx_dkx * c_x + norm_x * s_x * x_pos * dsgn_x)
  dyda(jb+2,j_term+2) = cmt%coef * s_y * c_z * sgn_y * (dnormy_dkx * s_x + norm_y * c_x * x_pos)
  dyda(jb+3,j_term+2) = cmt%coef * c_y * s_z * sgn_z * (dnormz_dkx * s_x + norm_z * c_x * x_pos)

  dyda(jb+1,j_term+3) = cmt%coef * c_x * c_y * sgn_x * (dnormx_dkz * c_z + norm_x * s_z * z_pos * dsgn_z)
  dyda(jb+2,j_term+3) = cmt%coef * s_x * s_y * sgn_y * (dnormy_dkz * c_z + norm_y * s_z * z_pos * dsgn_z)
  dyda(jb+3,j_term+3) = cmt%coef * s_x * c_y * sgn_z * (dnormz_dkz * s_z + norm_z * c_z * z_pos)

  jj = j_term+3
  if (.not. mask_x0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * s_x * c_y * c_z * sgn_x * cmt%kx * dsgn_x
    dyda(jb+2,jj) = cmt%coef * norm_y * c_x * s_y * c_z * sgn_y * cmt%kx 
    dyda(jb+3,jj) = cmt%coef * norm_z * c_x * c_y * s_z * sgn_z * cmt%kx 
  endif

  if (.not. mask_y0) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * s_y * c_z * sgn_x * cmt%ky * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * c_y * c_z * sgn_y * cmt%ky 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * s_y * s_z * sgn_z * cmt%ky * dsgn_y
  endif

  if (.not. mask_phi_z) then
    jj = jj + 1
    dyda(jb+1,jj) = cmt%coef * norm_x * c_x * c_y * s_z * sgn_x * dsgn_y
    dyda(jb+2,jj) = cmt%coef * norm_y * s_x * s_y * s_z * sgn_y * dsgn_y 
    dyda(jb+3,jj) = cmt%coef * norm_z * s_x * c_y * c_z * sgn_z
  endif

  ! Contribution from fact that ky is dependent upon kx and kx

  select case (cmt%form)
  case (hyper_y$)
    fx = cmt%coef * c_x * c_z * sgn_x * (dnormx_dky * c_y + norm_x * s_y * y_pos * dsgn_y) / cmt%ky
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y * y_pos) / cmt%ky
    fz = cmt%coef * s_x * s_z * sgn_z * (dnormz_dky * c_y + norm_z * s_y * y_pos * dsgn_y) / cmt%ky
  case (hyper_xy$)
    fx = cmt%coef * c_x * c_z * sgn_x * (norm_x * sinhx_x(cmt%ky*y_pos) * y_pos**2)
    fy = cmt%coef * s_x * c_z * sgn_y * (c_y + sinhx_x(cmt%ky*y_pos)) y_pos / cmt%kz
    fz = cmt%coef * s_x * s_z * sgn_z * (norm_z * sinhx_x(cmt%ky*y_pos) * y_pos**2)
  case (hyper_x$)
    fx = -cmt%coef * c_x * c_z * sgn_x * (norm_x * sinx_x(cmt%ky*y_pos) * y_pos**2)
    fy =  cmt%coef * s_x * c_z * sgn_y * (c_y - sinx_x(cmt%ky*y_pos)) * y_pos / cmt%kx
    fz = -cmt%coef * s_x * s_z * sgn_z * (norm_z * sinx_x(cmt%ky*y_pos) * y_pos**2)
  end select

  dyda(jb+1,j_term+2) = dyda(jb+1,j_term+2) + fx * cmt%kx * sgn_dky_dkx
  dyda(jb+2,j_term+2) = dyda(jb+2,j_term+2) + fy * cmt%kx * sgn_dky_dkx
  dyda(jb+3,j_term+2) = dyda(jb+3,j_term+2) + fz * cmt%kx * sgn_dky_dkx
  dyda(jb+1,j_term+3) = dyda(jb+1,j_term+3) + fx * cmt%kz * sgn_dky_dkz
  dyda(jb+2,j_term+3) = dyda(jb+2,j_term+3) + fy * cmt%kz * sgn_dky_dkz
  dyda(jb+3,j_term+3) = dyda(jb+3,j_term+3) + fz * cmt%kz * sgn_dky_dkz
end select

! coef derivative

jb = jb + 3

end subroutine add_to_field

end subroutine funcs_wf

end module
