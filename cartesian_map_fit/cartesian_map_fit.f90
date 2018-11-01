!+
! Program wiggler_fit
!
! Program to fit a field table for a wiggler field to the Bmad standard
! wiggler fit function.
!-

program wiggler_fit

use wig_mod
use nr

implicit none

real(rp), allocatable :: a(:), a1(:)
real(rp), allocatable :: covar(:,:), alpha(:,:)
real(rp) :: chisq, chisq1, chisq2, alamda, da, dy_max, dy, fret, tol

logical, allocatable :: maska(:)

integer i, j, nn, k, ijk
integer :: iter, nret, lib$get_foreign, istat
integer ix, iy, iz, cur_num
integer, allocatable :: ix_(:)

character(32) file_name
character(40) fname

!

open (1, file = 'number.in', readonly, shared)
read (1, *) cur_num
close(1)

file_name = 'fit'
write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'

call read_wiggler_fit_param_file (fname, .false.)

!

allocate (a(n_var), maska(n_var), a1(n_var))
allocate (covar(n_var,n_var), alpha(n_var,n_var))
allocate (ix_(n_term))

do i = 1, n_term
  if (vary_phi_z) then
    a(4*i-3:4*i) = (/ term(i)%coef, term(i)%kx, term(i)%kz, term(i)%phi_z /)
  else
    a(3*i-2:3*i) = (/ term(i)%coef, term(i)%kx, term(i)%kz /)
  endif
enddo

call funcs_wf(a, .false.) 

print *, '------------------------------------------------------'
print *, 'File: ', fname
call print_stuff

!

alamda = -1
maska = .true.
covar = 0

do ijk = 1, n_loops

  do i = 1, n_cycles
    call mrqmin_wf (a, maska, covar, alpha, chisq, alamda)
    type *, i, chisq, alamda
    if (alamda > 1e20 .or. alamda < 1e-20) exit
  enddo
 
  call funcs_wf(a, .false.) 
 
  !
 
  cur_num = cur_num + 1
  open (1, file = 'number.in')
  write (1, *) cur_num
  close(1)
  write (fname, '(a, i4.4, a)') trim(file_name), cur_num, '.in'
 
  print *, 'File: ', fname
  call print_stuff
  call db_rms_calc
 
  open (1, file = fname, carriagecontrol = 'list')
  write (1, *) 'Chi2:     ', merit_tot
  write (1, *) 'Chi2_coef:', merit_coef
  write (1, *) 'B_diff (G):', 1e4*B_diff / (3*n_pts)
  write (1, *) 'dB_rms (G):', dB_rms
  write (1, *) 'B_Merit   :', B_diff / sumB
  write (1, *) 'B_int:', B_int
  write (1, *) 'Weight_coef: ', weight_coef_all
  write (1, *) 
  write (1, *) '&Parameters'
  write (1, *) '  field_file = "', trim(field_file), '"'
  write (1, *) '  Nz_ref     =', Nz_ref
  call indexx(term(1:n_term)%kz, ix_)
  do i = 1, n_term
    ix = ix_(i)
    write (1, '(1x, a, i3, a, 3f13.7, f13.8)') '  term(', i, ') =', &
                term(ix)%coef, term(ix)%kx, term(ix)%kz, term(ix)%phi_z
  enddo
  write (1, *) '/end'
  close (1)

  if (alamda > 1e20 .or. alamda < 1e-20) exit

enddo


end program
