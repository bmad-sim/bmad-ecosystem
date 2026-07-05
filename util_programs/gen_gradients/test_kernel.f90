program test_kernel
use gen_gradients_mod
implicit none

real(rp) :: g_ref, x, y
real(rp) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
real(rp) :: B(3), A(3), dA(3,3)
real(rp) :: Bref(3), Aref(3), dAref(3,3)
real(rp) :: v, maxerr, e
integer :: nent, i, j, knd, n, m, u

open(newunit=u, file='gg_ref_case.txt', status='old')
read(u,*) g_ref, x, y
read(u,*) nent
gval = 0
do i = 1, nent
  read(u,*) knd, n, m, v
  gval(knd, n, m) = v
enddo
read(u,*) Bref(1), Bref(2), Bref(3)
read(u,*) Aref(1), Aref(2), Aref(3)
do i = 1, 3
  read(u,*) dAref(i,1), dAref(i,2), dAref(i,3)
enddo
close(u)

call gg_field_potential_calc(g_ref, gval, x, y, .true., .true., B, A, dA)

maxerr = 0
call chk('B ', B, Bref, 3, maxerr)
call chk('A ', A, Aref, 3, maxerr)
do i = 1, 3
  call chk('dA', dA(i,:), dAref(i,:), 3, maxerr)
enddo

print '(a, es12.4)', 'MAX RELATIVE ERROR = ', maxerr
if (maxerr < 1e-10_rp) then
  print '(a)', 'KERNEL VALIDATION PASSED'
else
  print '(a)', 'KERNEL VALIDATION FAILED'
  call exit(1)
endif

contains
subroutine chk(lbl, got, ref, nn, mx)
character(*) lbl
integer nn, k
real(rp) got(nn), ref(nn), mx, r
do k = 1, nn
  r = abs(got(k)-ref(k)) / max(abs(ref(k)), 1e-30_rp)
  if (r > mx) mx = r
  print '(2a, i2, 3(a, es22.14))', lbl, ' idx=', k, '  got=', got(k), ' ref=', ref(k), ' relerr=', r
enddo
end subroutine
end program
