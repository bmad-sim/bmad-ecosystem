!+
! Subroutine tao_regression_test(cmd_str)
!
! This routine generates output for the tao_test regression test in the regression_tests/tao_test directory.
!-

subroutine tao_regression_test(cmd_str)

use tao_set_mod, dummy => tao_regression_test
use pointer_lattice, only: operator(.sub.), operator(**), operator(*), alloc, kill, print, ci_phasor, assignment(=)

implicit none

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ptc_normal_form_struct), pointer :: ptc_nf
type (tao_lattice_struct), pointer :: tao_lat

real(rp) r, z1, z2
real(rp), allocatable :: val(:)
integer i, iu, ii, expo(6)

logical err, dflt_source(16)

character(*) cmd_str
character(200) excite_zero(3), veto
character(60) :: expr(16) = [character(60):: &
                  'q*d[2:4]|design', &
                  '[anomalous_moment_of(proton), mass_of(electron)]', &
                  '(46.5/anomalous_moment_of(proton))^2-pi', &
                  '[3,4] * [1]@ele::q1[k1]', &
                  'pi + 10**@data::*|model*7', &
                  '100*var::*[*]|model*10', &
                  'mass_of(#3He+2) / charge_of(Si++)', &
                  'bbb*charge_of(aaa)', &
                  'max([-1,-2,-3]) + min(4,2,3) + sum(-3,-1.0e-1)', &
                  'atan2(1,2)', &
                  'data::twiss.end[1]|model-design+1e-10', &
                  '[1,2] - [3, lat::r.11[beginning&end->-0.5*l]]', &
                  'sum(abs(ele::sbend::b*[angle]))', &
                  'ran_gauss() + ran_gauss(0.1) + ran()', &
                  'ele::q1-1[l]', &
                  'ele::q1+1[l]' &
                        ]

!

iu = lunget()
select case (cmd_str)
case ('append')
  open (iu, file = 'output.now', access = 'append')
case default
  open (iu, file = 'output.now')
end select

!

u => s%u(1)
branch => u%model%lat%branch(0)
tao_branch => u%model%tao_branch(0)

! Closed lattice regressions
! For some reason the regression tests on Mac shows large value deviations so the tollerance is set to 1e-2.

if (branch%param%geometry == closed$) then
  tao_lat => tao_pointer_to_tao_lat (u, model$)
  if (.not. u%calc%one_turn_map) call tao_ptc_normal_form (.true., tao_lat, branch%ix_branch)
  ptc_nf  => tao_branch%ptc_normal_form

  do i = 0, ptc_private%taylor_order_ptc-1
    expo = [0, 0, 0, 0, 0, i]
    if (ptc_nf%state%nocavity) then
      z1 =  real(ptc_nf%phase(1) .sub. expo)
      z2 =  real(ptc_nf%phase(2) .sub. expo)
      write (iu, '(a, i0, a, 2es18.7)') '"Chrom-noRF-', i, '" REL 1e-2', z1, z2
    else
      z1 =  real(ptc_nf%u_phase(1) .sub. expo)
      z2 =  real(ptc_nf%u_phase(2) .sub. expo)
      write (iu, '(a, i0, a, 2es18.7)') '"Chrom-RF-', i, '" REL 1e-2', z1, z2
    endif
  enddo

  do i = 1, ptc_private%taylor_order_ptc-1
    expo = [0, 0, 0, 0, 0, i]
    if (ptc_nf%state%nocavity) then
      z1 = -real(ptc_nf%phase(3) .sub. expo) / branch%param%total_length
      z2 =  real(ptc_nf%path_length .sub. expo) / branch%param%total_length
      write (iu, '(a, i0, a, 2es18.7)') '"Slip-noRF-', i, '" REL 1e-2', z1, z2
    else
      z1 = -real(ptc_nf%u_phase(3) .sub. expo) / branch%param%total_length
      z2 =  real(ptc_nf%u_path_length .sub. expo) / branch%param%total_length
      write (iu, '(a, i0, a, 2es18.7)') '"Slip-RF-', i, '" REL 1e-2', z1, z2
    endif
  enddo

  return
endif

! Open lattice regressions

s%global%expression_tree_on = .true.
call tao_set_symbolic_number_cmd ('aaa', 'species(Li+5)')
call tao_set_symbolic_number_cmd ('bbb', '34*2')
call ran_seed_put(1234)
dflt_source = .false.
dflt_source(1) = .true.

do ii = 1, size(expr)
  if (dflt_source(ii)) then
    call tao_evaluate_expression(expr(ii), 0, .false., val, err, dflt_source = 'var')
  else
    call tao_evaluate_expression(expr(ii), 0, .false., val, err)
  endif
  write (iu, '(2a, 9es18.10)') quote(expr(ii)), ' REL 1E-9', val
enddo


!

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
