!+
! Program cartesian_map_fit
!
! Program to fit a cartesian map to a field table.
!-

program cartesian_map_fit

use cartesian_map_fit_mod
use super_recipes_mod
use opti_de_mod, only: opti_de !MPI , opti_de_mpi

implicit none

type (super_mrqmin_storage_struct) storage
type (coord_struct) orbit
type (em_field_struct) field
type (cartesian_map_term1_struct), pointer :: c1, cmt
type (term_struct), pointer :: tt

real(rp), allocatable :: a(:), a0(:), y_target(:), y0(:), y1(:)
real(rp), allocatable :: covar(:,:), alpha(:,:), weight(:), dyda(:,:), var_step(:), dyda2(:,:)
real(rp) :: chisq, chisq1, chisq2, alamda, da, dy_max, dy, fret, tol, z, delta
real(rp) dfield(3), merit, dyda_err, dyda_this, rel_err

logical, allocatable :: maska(:)

integer i, j, jj, nn, k, ijk, status, population
integer :: iter, nret, lib$get_foreign, istat
integer ix, iy, iz, cur_num

character(16) :: mode, mode_in, num_str, arg_in
character(100) file_name
character(100) fname

!

mode_in = 'STD'
num_str = ''

if (command_argument_count() > 0) then
  call get_command_argument(1, arg_in)
  if (is_integer(arg_in)) then
    num_str = arg_in
    if (command_argument_count() > 1) call get_command_argument(2, mode_in)
  else
    mode_in = arg_in
  endif
endif

!

call match_word (mode_in, [character(16):: 'STD', 'fft', 'debug', 'binary', 'fit_table'], ix, .false., .true., mode)

if (ix == 0) then
  print *, 'I DO NOT UNDERSTAND MODE: ', mode_in
  stop
endif

if (num_str == '') then
  open (1, file = 'number.in', status = 'OLD', action = 'READ')
  read (1, *) cur_num
  close(1)
else
  if (.not. is_integer(num_str)) then
    print *, 'I DO NOT UNDERSTAND NUMBER: ', num_str
    stop
  endif
  read (num_str, *) cur_num
endif

file_name = 'fit'
write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'

call read_cartesian_map_fit_param_file (fname, .false.)

!

allocate (a(n_var), maska(n_var), y_target(n_data_tot), weight(n_data_tot))
allocate (covar(n_var,n_var), alpha(n_var,n_var))
allocate (y_fit(n_data_tot), dyda(n_data_tot,n_var), var_step(n_var))

y_target = 0

do i = 1, n_term
  jj = (i - 1) * n_var_per_term 
  jj=jj+1; a(jj) = term(i)%cmt%coef; var_step(jj) = de_coef_step
  jj=jj+1; a(jj) = term(i)%kxy;      var_step(jj) = de_k_step
  jj=jj+1; a(jj) = term(i)%cmt%kz;   var_step(jj) = de_k_step

  if (.not. mask_x0) then
    jj=jj+1; a(jj) = term(i)%cmt%x0; var_step(jj) = de_x0_y0_step
  endif

  if (.not. mask_y0) then
    jj=jj+1; a(jj) = term(i)%cmt%y0; var_step(jj) = de_x0_y0_step
  endif

  if (.not. mask_phi_z) then
    jj=jj+1; a(jj) = term(i)%cmt%phi_z; var_step(jj) = de_phi_z_step
  endif
enddo

print *, '------------------------------------------------------'
print *, 'File: ', trim(fname)

call print_stuff(0, a)

! Write binary table
! Notice that the binary table stores lengths in meters and fields in Tesla independent of length_scale and field_scale.

if (mode == 'binary') then
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
  close (1)
  stop
endif

! Write a fit table

if (mode == 'fit_table') then
  open (1, file = 'fit.table')
  write (1, '(2a)') real_to_string(length_scale, 16, 4),   '! length_scale'
  write (1, '(2a)') real_to_string(field_scale, 16, 4),    '! field_scale'
  write (1, '(2i6, 4x, a)') Nx_min, Nx_max, '! Nx_min, Nx_max'
  write (1, '(2i6, 4x, a)') Ny_min, Ny_max, '! Ny_min, Ny_max'
  write (1, '(2i6, 4x, a)') Nz_min, Nz_max, '! Nz_min, Nz_max'
  write (1, '(a, 4x, a)') reals_to_table_row(del_grid/length_scale, 16, 6, 2),       '! del_grid'
  write (1, '(a, 4x, a)') reals_to_table_row(r0_grid/length_scale, 16, 6, 2),        '! r0_grid'

  do ix = Nx_min, Nx_max
  do iy = Ny_min, Ny_max
  do iz = Nz_min, Nz_max
    write (1, '(3f13.6, 3es16.7)') [ix*del_grid(1), iy*del_grid(2), iz*del_grid(3)] / length_scale, &
                     [Bx_fit(ix, iy, iz), By_fit(ix, iy, iz), Bz_fit(ix, iy, iz)] / field_scale
  enddo
  enddo
  enddo

  print *, 'Wrote fit table file: fit.table'
  stop
endif

! Debug. 
! Note: When debugging use a map with only one term.

if (mode == 'debug') then
  calc_at_surface_grid_points_only = .true.
  call funcs_lm(a, y_fit, dyda, status) 

  allocate (a0(n_var), y0(n_data_tot), y1(n_data_tot), dyda2(n_data_tot, n_var))
  dyda_err = 0
  a0 = a
  delta = 1e-5
  do i = 1, n_var
    a(i) = a0(i) - delta
    call funcs_lm(a, y0, dyda2, status) 
    a(i) = a0(i) + delta
    call funcs_lm(a, y1, dyda2, status) 
    a(i) = a0(i)

    do j = 1, size(y_fit)
      dyda_this = (y1(j) - y0(j)) / (2 * delta)
      if (dyda_this == 0 .and. dyda(j,i) == 0) cycle
      rel_err = abs(dyda_this - dyda(j,i)) / (abs(dyda_this) + abs(dyda(j,i)))
      if (rel_err > dyda_err) then
        dyda_err = rel_err
        print *, 'Max err:', i, j, dyda_err
      endif
    enddo
  enddo
  stop

  call init_coord (orbit, orbit, lat%ele(0), downstream_end$)

  c1 => c_map%ptr%term(1)
  tt => term(1)
  cmt => tt%cmt

  c1%coef = cmt%coef   ! There may be a sign flip
  c1%kx   = cmt%kx
  c1%ky   = cmt%ky

  jj = 0
  do i = Nx_min, Nx_max
  do j = Ny_min, Ny_max
  do k = Nz_min, Nz_max
    orbit%vec(1) = i * del_grid(1) - r0_grid(1)
    orbit%vec(3) = j * del_grid(2) - r0_grid(2)
    z = k * del_grid(3) - r0_grid(3)
    call em_field_calc (lat%ele(1), lat%param, z, orbit, .true., field)
    print *
    print *, '!--------------------------------------------------'
    print '(3i4, 3es14.6)', i, j, k, field%B
    print '(3i4, 3es14.6)', i, j, k, Bx_fit(i,j,k), By_fit(i,j,k), Bz_fit(i,j,k)

    call get_this_derivative ('dCoef: ', c1%coef,  1d-5, dyda(jj+1:jj+3, 1), dfield, 0)
    call get_this_derivative ('dkxy:  ', tt%kxy,   1d-5, dyda(jj+1:jj+3, 2), dfield, 1)
    call get_this_derivative ('dkz:   ', c1%kz,    1d-5, dyda(jj+1:jj+3, 3), dfield, 1)
    call get_this_derivative ('dx0:   ', c1%x0,    1d-5, dyda(jj+1:jj+3, 4), dfield, 0)
    call get_this_derivative ('dy0:   ', c1%y0,    1d-5, dyda(jj+1:jj+3, 5), dfield, 0)
    call get_this_derivative ('dphi_z:', c1%phi_z, 1d-5, dyda(jj+1:jj+3, 6), dfield, 0)

    jj = jj + 3
  enddo
  enddo
  enddo

  stop
endif

! FFT

if (mode == 'fft') then
  calc_at_surface_grid_points_only = .false.
  merit = funcs_de(a, status, jj)
  call fft_field()
  stop
endif

!-----------------------------------------------------------------------------
! Optimize

alamda = -1
maska = .true.
covar = 0
weight = 1
weight(n_data_grid+1:n_data_tot:3) = coef_weight
weight(n_data_grid+2:n_data_tot:3) = k_weight
weight(n_data_grid+3:n_data_tot:3) = k_weight
population = de_var_to_population_factor * n_var

calc_at_surface_grid_points_only = .true.
call funcs_lm(a, y_fit, dyda, status) 

do ijk = 1, n_loops

  select case (optimizer)
  case ('de')
    dyda_calc = .false.
    merit_tot = opti_de (a, n_cycles, population, funcs_de, var_step, status)


  case ('lm')
    dyda_calc = .true.
    do i = 1, n_cycles
      call super_mrqmin (y_target, weight, a, chisq, funcs_lm, storage, alamda, status, maska)
      print *, i, chisq, alamda
      if (alamda > 1e20 .or. alamda < 1e-20) exit
    enddo

    call funcs_lm(a, y_fit, dyda, status) 

  case default
    print *, 'Unknown optimizer: ', optimizer
    stop
  end select

  !
 
  cur_num = cur_num + 1
  open (1, file = 'number.in')
  write (1, *) cur_num
  close(1)
  write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'
 
  print *, 'File: ', fname
 
  call write_bmad_lattice_file (fname, lat, output_form = one_file$)

  open (1, file = fname, status = 'old', position = 'append')
  write (1, *)
  write (1, *) 'end_file'
  write (1, *) 
  write (1, *) '&Parameters'
  write (1, *) '  field_file   = "', trim(field_file), '"'
  write (1, *) '  optimizer    = "', trim(optimizer), '"'
  write (1, *) '  coef_weight  = ', coef_weight
  write (1, *) '  k_weight     = ', k_weight
  write (1, *) '  n_loops      = ', n_loops
  write (1, *) '  n_cycles     = ', n_cycles
  write (1, *) '  mask_x0      = ', mask_x0
  write (1, *) '  mask_y0      = ', mask_y0
  write (1, *) '  mask_phi_z   = ', mask_phi_z
  write (1, *) '  de_var_to_population_factor =', de_var_to_population_factor
  write (1, *) '  de_coef_step       =', de_coef_step
  write (1, *) '  de_k_step          =', de_k_step
  write (1, *) '  de_x0_y0_step      =', de_x0_y0_step
  write (1, *) '  de_phi_z_step      =', de_phi_z_step
  write (1, *) '/end'
  write (1, *)
  call print_stuff(1, a)
  close (1)

  if (alamda > 1e20 .or. alamda < 1e-20) exit

enddo

!----------------------------------------------------------
contains

subroutine get_this_derivative (name, var, delta, dBda, dfield, vtype)

type (em_field_struct) f0, f1
real(rp) var, delta, dfield(3), dBda(3)
integer vtype
character(*) name

call vary_var(var, -delta, vtype) 
call em_field_calc (lat%ele(1), lat%param, z, orbit, .true., f0)

call vary_var(var, 2*delta, vtype) 
call em_field_calc (lat%ele(1), lat%param, z, orbit, .true., f1)

call vary_var(var, -delta, vtype) 
dfield = (f1%B - f0%B) / (2 * delta)

print *
print '(a, 3es14.6)', name, dfield
print '(a, 3es14.6)', name, dBda
print '(a, 3es14.6)', name, dfield - dBda

end subroutine

!---------------------------------------------------------
! contains

subroutine vary_var (var, delta, vtype)

real(rp) var, delta, kxy
integer vtype

!

var = var + delta
if (vtype == 0) return

select case (tt%cmt%form)
case (hyper_y$)
  cmt%ky = (tt%kxy**2 + cmt%kz**2) / (2 * tt%kxy)
  cmt%kx = sqrt(cmt%ky**2 - cmt%kz**2)

case (hyper_x$)
  cmt%kx = (tt%kxy**2 + cmt%kz**2) / (2 * tt%kxy)
  cmt%ky = -sqrt(cmt%kx**2 - cmt%kz**2)

case (hyper_xy$)
  cmt%kx = (tt%kxy - sqrt(2 * cmt%kz**2 - tt%kxy**2)) / 2
  cmt%ky = sqrt(cmt%kz**2 - cmt%kx**2)
end select

c1%kx = cmt%kx
c1%ky = cmt%ky

end subroutine vary_var

!---------------------------------------------------------
! contains

subroutine fft_field()

complex(rp), allocatable :: fx(:), fy(:), fz(:)
complex(rp) fc
real(rp), allocatable :: amp(:)
real(rp) ff, norm

integer ix, iy, N, i, j, if
integer, allocatable :: indx(:)

!

N = Nz_max - Nz_min + 1
allocate (fx(0:N/2), fy(0:N/2), fz(0:N/2), amp(0:N/2), indx(0:N/2))

ix = -nint((x_offset_map - r0_grid(1)) / del_grid(1))
iy = -nint((y_offset_map - r0_grid(2)) / del_grid(2))

if (ix < Nx_min .or. ix > Nx_max) then
  print *, 'X = 0 LIES OUTSIDE OF THE GRID! CANNOT DO FFT ANALYSIS!'
  stop
endif

if (iy < Ny_min .or. iy > Ny_max) then
  print *, 'Y = 0 LIES OUTSIDE OF THE GRID! CANNOT DO FFT ANALYSIS!'
  stop
endif

fx = 0; fy = 0; fz = 0

do iz = Nz_min, Nz_max
  do if = 0, N/2
    fc = exp(-twopi * I_imag * if * iz / N)
    fx(if) = fx(if) + (Bx_fit(ix,iy,iz) - Bx_dat(ix,iy,iz)) * fc
    fy(if) = fy(if) + (By_fit(ix,iy,iz) - By_dat(ix,iy,iz)) * fc
    fz(if) = fz(if) + (Bz_fit(ix,iy,iz) - Bz_dat(ix,iy,iz)) * fc
  enddo
enddo

norm = 2.0 / (Nz_max - Nz_min + 1)

amp = abs(fx)
call indexer(-amp, indx)
print *, 'B_x FFT (X-family):'
print *,   '    kz            Amp           Phase'
do i = 1, min(5, size(indx))
  j = indx(i-1) - 1
  ff = twopi * j / (N * del_grid(3))
  print '(3es14.4)', ff, norm*amp(j), atan2(aimag(fx(j)), real(fx(j)))
enddo

amp = abs(fy)
call indexer(-amp, indx)
print *
print *, 'B_y FFT (Y-family):'
print *,   '    kz            Amp           Phase'
do i = 1, min(5, size(indx))
  j = indx(i-1) - 1
  ff = twopi * j / (N * del_grid(3))
  print '(3es14.4)', ff, norm*amp(j), atan2(aimag(fy(j)), real(fy(j)))
enddo

amp = abs(fz)
call indexer(-amp, indx)
print *
print *, 'B_z FFT (SQ-family):'
print *,   '    kz            Amp           Phase'
do i = 1, min(5, size(indx))
  j = indx(i-1) - 1
  ff = twopi * j / (N * del_grid(3))
  print '(3es14.4)', ff, norm*amp(j), atan2(aimag(fz(j)), real(fz(j)))
enddo

end subroutine fft_field

end program
