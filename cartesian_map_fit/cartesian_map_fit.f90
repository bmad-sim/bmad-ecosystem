!+
! Program cartesian_map_fit
!
! Program to fit a cartesian map to a field table.
!-

program cartesian_map_fit

use cartesian_map_fit_mod
use super_recipes_mod
use write_lat_file_mod

implicit none

type (super_mrqmin_storage_struct) storage
type (coord_struct) orbit
type (em_field_struct) field
type (cartesian_map_term1_struct), pointer :: c1, cmt
type (term_struct), pointer :: tt

real(rp), allocatable :: a(:), y0(:)
real(rp), allocatable :: covar(:,:), alpha(:,:), weight(:), yfit(:), dyda(:,:)
real(rp) :: chisq, chisq1, chisq2, alamda, da, dy_max, dy, fret, tol, z, delta
real(rp) dfield(3)

logical, allocatable :: maska(:)

integer i, j, jj, nn, k, ijk, status
integer :: iter, nret, lib$get_foreign, istat
integer ix, iy, iz, cur_num

character(32) file_name
character(40) fname

!

open (1, file = 'number.in', readonly, shared)
read (1, *) cur_num
close(1)

file_name = 'fit'
write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'

call read_cartesian_map_fit_param_file (fname, .false.)

!

allocate (a(n_var), maska(n_var), y0(n_data), weight(n_var))
allocate (covar(n_var,n_var), alpha(n_var,n_var))
allocate (yfit(n_data), dyda(n_data,n_var))

y0 = 0

do i = 1, n_term
  jj = (i - 1) * n_var_per_term 
  jj=jj+1; a(jj) = term(i)%cmt%coef
  jj=jj+1; a(jj) = term(i)%kxy
  jj=jj+1; a(jj) = term(i)%cmt%kz

  if (.not. mask_x0) then
    jj=jj+1; a(jj) = term(i)%cmt%x0
  endif

  if (.not. mask_y0) then
    jj=jj+1; a(jj) = term(i)%cmt%y0
  endif

  if (.not. mask_phi_z) then
    jj=jj+1; a(jj) = term(i)%cmt%phi_z
  endif
enddo

call funcs_cf(a, yfit, dyda, status) 

print *, '------------------------------------------------------'
print *, 'File: ', fname
call print_stuff

! Debug

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
  orbit%vec(1) = i * del(1)
  orbit%vec(3) = j * del(2)
  z = k * del(3)
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

!

alamda = -1
maska = .true.
covar = 0

do ijk = 1, n_loops

  do i = 1, n_cycles
    call super_mrqmin (y0, weight, a, chisq, funcs_cf, storage, alamda, status, maska)
    type *, i, chisq, alamda
    if (alamda > 1e20 .or. alamda < 1e-20) exit
  enddo
 
  call funcs_cf(a, yfit, dyda, status) 
 
  !
 
  cur_num = cur_num + 1
  open (1, file = 'number.in')
  write (1, *) cur_num
  close(1)
  write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'
 
  print *, 'File: ', fname
  call print_stuff
  call db_rms_calc
 
  call write_bmad_lattice_file (fname, lat)
  open (1, file = fname, status = 'old')

  write (1, *)
  write (1, *) 'end_file'
  write (1, *)
  write (1, *) 'Chi2:     ', merit_tot
  write (1, *) 'Chi2_coef:', merit_coef
  write (1, *) 'B_diff (G):', 1e4*B_diff / (3*n_grid_pts)
  write (1, *) 'dB_rms (G):', dB_rms
  write (1, *) 'B_Merit   :', B_diff / sumB_in
  write (1, *) 'B_int:', B_int
  write (1, *) 
  write (1, *) '&Parameters'
  write (1, *) '  field_file   = "', trim(field_file), '"'
  write (1, *) '  coef_weight  = ', coef_weight
  write (1, *) '  n_loops      = ', n_loops
  write (1, *) '  n_cycles     = ', n_cycles
  write (1, *) '  fft_write    = ', fft_write
  write (1, *) '  z0_phase_ref = ', z0_phase_ref
  write (1, *) '  mask_x0      = ', mask_x0
  write (1, *) '  mask_y0      = ', mask_y0
  write (1, *) '  mask_phi_z   = ', mask_phi_z
  write (1, *) '/end'
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

!--------------------------------------------
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

end program
