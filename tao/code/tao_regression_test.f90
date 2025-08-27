!+
! Subroutine tao_regression_test()
!
! This routine generates output for the tao_test regression test in the regression_tests/tao_test directory.
!-

subroutine tao_regression_test()

use tao_set_mod, dummy => tao_regression_test

implicit none

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch

real(rp) r
real(rp), allocatable :: val(:)
integer iu, ii
logical err

character(200) excite_zero(3), veto
character(60) :: expr(10) = [character(60):: &
                  '[anomalous_moment_of(proton), mass_of(electron)]', &
                  '(46.5/anomalous_moment_of(proton))^2-pi', &
                  '[3,4] * [1]@ele::q1[k1]', &
                  'pi + 10**@data::*|model*7', &
                  '100*var::*[*]|model*10', &
                  'mass_of(#3He+2) / charge_of(Si++)', &
                  'bbb*charge_of(aaa)', &
                  'max([1,2,3]) + min(4,2,3)', &
                  'atan2(1,2)', &
                  '1' &
                ]

!

iu = lunget()
open (iu, file = 'output.now')

!

s%global%expression_tree_on = .true.
call tao_set_symbolic_number_cmd ('aaa', 'species(Li+5)')
call tao_set_symbolic_number_cmd ('bbb', '34*2')

do ii = 1, size(expr)
  call tao_evaluate_expression(expr(ii), 0, .false., val, err)
  write (iu, '(2a, 9es18.10)') quote(expr(ii)), ' REL 1E-9', val
enddo


!

u => s%u(1)
branch => u%model%lat%branch(0)
tao_branch => u%model%tao_branch(0)
excite_zero = ''
veto = ''

call tao_spin_polarization_calc (branch, tao_branch, excite_zero, veto)

write(iu, '(a, es18.10)')  '"spin-tune" ABS 4E-9', tao_branch%spin%tune / twopi
write(iu, '(a, es18.10)')  '"Polarization Limit ST" REL 4E-9                  ', tao_branch%spin%pol_limit_st
write(iu, '(a, es18.10)')  '"Polarization Limit DK" REL 4E-9                   ', tao_branch%spin%pol_limit_dk
write(iu, '(a, 3es18.10)') '"Polarization Limits DK (a,b,c-modes)" REL 4E-9    ', tao_branch%spin%pol_limit_dk_partial
write(iu, '(a, 3es18.10)') '"Polarization Limits DK (bc,ac,ab-modes)" REL 4E-9 ', tao_branch%spin%pol_limit_dk_partial2
write(iu, '(a, es18.10)')  '"Polarization Rate BKS" REL 4E-9', tao_branch%spin%pol_rate_bks
write(iu, '(a, es18.10)')  '"Depolarization Rate" REL 4E-9',   tao_branch%spin%depol_rate
write(iu, '(a, 3es18.10)') '"Depolarization Rate" REL 4E-9',   tao_branch%spin%depol_rate_partial
write(iu, '(a, 3es18.10)') '"Depolarization Rate" REL 4E-9',   tao_branch%spin%depol_rate_partial2
write(iu, '(a, es18.10)')  '"Integral g^3 * b_hat * n_0" REL 4E-9         ', tao_branch%spin%integral_bn
write(iu, '(a, es18.10)')  '"Integral g^3 * b_hat * dn/ddelta" REL 4E-9   ', tao_branch%spin%integral_bdn
write(iu, '(a, es18.10)')  '"Integral g^3 (1 - 2(n * s_hat)/9)" REL 4E-9  ', tao_branch%spin%integral_1ns
write(iu, '(a, es18.10)')  '"Integral g^3 * 11 (dn/ddelta)^2 / 9" REL 4E-9', tao_branch%spin%integral_dn2

!

close(iu)
end subroutine
