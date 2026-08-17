!+
! Program tao_nml_test
!
! Equivalence test for the Tao init file reader in tao_nml_mod.
!
! Back compatibility for Tao init files means the new reader must accept exactly what a Fortran
! namelist read accepted and produce the same values. Rather than assert that from reasoning,
! this test reads the *same file* twice: once with a real namelist read and once with the new
! reader, and compares the results component by component. Any divergence shows up as a
! mismatch in output.now.
!
! The test file exercises the syntax Tao init files actually use: comments, quoted and
! unquoted strings, embedded quote marks, multiple values per item, values continued across
! lines, repeat counts, short value lists for an array, subscripted items, structure
! components, several groups in one file, and both the "/" and "&end" terminators.
!-

program tao_nml_test

use tao_nml_mod
use tao_attrib_resolve_mod

implicit none

type (tao_nml_group_struct) group
type (tao_ptr_struct) ptr
type (tao_datum_input), target :: dat
type (tao_nml_ref_struct) item_ref

real(rp) :: rv, arr(6), rng(8), rep(4), tabbed(3), null_mid(3)
integer :: iv, n_mismatch, iu, i, isub
logical :: lv, err, eof, has_sub
character(40) :: sv
character(200) :: why
character(:), allocatable :: head, rest

! Saved results of the namelist pass.

real(rp) :: rv_nml, arr_nml(6), rng_nml(8), rep_nml(4), tabbed_nml(3), null_mid_nml(3)
integer :: iv_nml
logical :: lv_nml
character(40) :: sv_nml
type (tao_datum_input) :: dat_nml

namelist / testgrp / rv, iv, lv, sv, arr, dat, rng, rep, tabbed, null_mid

!

call write_test_file()

n_mismatch = 0
open (newunit = iu, file = 'output.now', recl = 300)

! Pass 1: the reference. A real Fortran namelist read.

call reset_vars()
block
  integer iu2
  open (newunit = iu2, file = 'nml_test_input.nml', status = 'old')
  read (iu2, nml = testgrp)
  close (iu2)
end block

rv_nml = rv;  iv_nml = iv;  lv_nml = lv;  sv_nml = sv;  arr_nml = arr;  dat_nml = dat
rng_nml = rng;  rep_nml = rep;  tabbed_nml = tabbed;  null_mid_nml = null_mid

! Pass 2: the new reader.

call reset_vars()
block
  integer iu2
  open (newunit = iu2, file = 'nml_test_input.nml', status = 'old')
  call tao_nml_group_read (iu2, 'nml_test_input.nml', 'testgrp', group, eof, err, why)
  close (iu2)

  if (err .or. eof) then
    write (iu, '(3a)') '"reader-status"  STR  ', quote('FAILED: ' // trim(why)), ' '
    close (iu)
    stop
  endif

  do i = 1, group%n_item
    call tao_nml_item_ref (group%item(i), item_ref, err)
    head = item_ref%head
    rest = item_ref%rest
    if (err) cycle

    select case (head)
    case ('rv');   call tao_nml_value_set (group%item(i), rv, err)
    case ('iv');   call tao_nml_value_set (group%item(i), iv, err)
    case ('lv');   call tao_nml_value_set (group%item(i), lv, err)
    case ('sv');   call tao_nml_value_set (group%item(i), sv, err)

    case ('arr');       call set_arr (group%item(i), arr)
    case ('rng');       call set_arr (group%item(i), rng)
    case ('rep');       call set_arr (group%item(i), rep)
    case ('tabbed');    call set_arr (group%item(i), tabbed)
    case ('null_mid');  call set_arr (group%item(i), null_mid)

    case ('dat')
      call tao_res_tao_datum_input (dat, rest, ptr, err, why)
      if (err) then
        call tao_nml_err (group%item(i), why)
      else
        call tao_nml_value_set (group%item(i), ptr, err)
      endif

    case default
      call tao_nml_unknown (group%item(i), 'testgrp', err)
    end select
  enddo
end block

! Compare.

call cmp_real ('rv', rv_nml, rv)
call cmp_int ('iv', iv_nml, iv)
call cmp_logic ('lv', lv_nml, lv)
call cmp_str ('sv', sv_nml, sv)
do i = 1, 6
  call cmp_real ('arr' // int_str(i), arr_nml(i), arr(i))
enddo
do i = 1, 8
  call cmp_real ('rng' // int_str(i), rng_nml(i), rng(i))
enddo
do i = 1, 4
  call cmp_real ('rep' // int_str(i), rep_nml(i), rep(i))
enddo
do i = 1, 3
  call cmp_real ('tab' // int_str(i), tabbed_nml(i), tabbed(i))
enddo
do i = 1, 3
  call cmp_real ('nul' // int_str(i), null_mid_nml(i), null_mid(i))
enddo
call cmp_str ('dat%ele_name', dat_nml%ele_name, dat%ele_name)
call cmp_str ('dat%data_type', dat_nml%data_type, dat%data_type)
call cmp_str ('dat%merit_type', dat_nml%merit_type, dat%merit_type)
call cmp_real ('dat%weight', dat_nml%weight, dat%weight)
call cmp_logic ('dat%good_user', dat_nml%good_user, dat%good_user)

write (iu, '(a, i0)') '"n-mismatch"        ABS 0     ', n_mismatch
write (iu, '(a, i0)') '"n-items"           ABS 0     ', group%n_item

! Multiple groups in one file, and skipping to a wanted group.

block
  integer iu2
  type (tao_nml_group_struct) g2
  open (newunit = iu2, file = 'nml_test_input.nml', status = 'old')
  call tao_nml_group_read (iu2, 'nml_test_input.nml', 'second_group', g2, eof, err, why)
  close (iu2)
  write (iu, '(3a)') '"skip-to-group"     STR       ', quote(trim(g2%name)), ''
  write (iu, '(a, l1, a)') '"skip-err"          STR       "', err, '"'
  if (g2%n_item >= 1) then
    write (iu, '(3a)') '"second-item-1"     STR       ', quote(trim(g2%item(1)%name)), ''
    write (iu, '(3a)') '"second-value-1"    STR       ', quote(trim(g2%item(1)%value)), ''
  endif
end block

! A group with no terminating "/" is ended by end of file rather than being an error.
!
! This is a deliberate difference from a namelist read, not a match: a namelist read of an
! unterminated group aborts, or with iostat= reports end of file and discards whatever it had
! read. Neither is useful, and files like this do occur (Tao's own "write namelist -variable"
! can crash partway and leave one behind). No well formed file is affected, since a well formed
! file terminates its groups, so nothing that used to load loads differently.
!
! There is no namelist reference to compare against here for that reason, so the expected
! values are asserted directly.

block
  integer iu2
  type (tao_nml_group_struct) g3

  open (newunit = iu2, file = 'nml_test_bad.nml', status = 'replace')
  write (iu2, '(a)') '&testgrp'
  write (iu2, '(a)') '  rv = 8.25'
  write (iu2, '(a)') '  iv = 77'
  close (iu2)

  call reset_vars()
  open (newunit = iu2, file = 'nml_test_bad.nml', status = 'old')
  call tao_nml_group_read (iu2, 'nml_test_bad.nml', 'testgrp', g3, eof, err, why)
  close (iu2)

  do i = 1, g3%n_item
    call tao_nml_item_ref (g3%item(i), item_ref, err)
    if (err) cycle
    select case (item_ref%head)
    case ('rv');  call tao_nml_value_set (g3%item(i), rv, err)
    case ('iv');  call tao_nml_value_set (g3%item(i), iv, err)
    end select
  enddo

  write (iu, '(a, l1, a)') '"unterminated-err"  STR       "', err, '"'
  write (iu, '(a, l1, a)') '"unterminated-eof"  STR       "', eof, '"'
  write (iu, '(a, es14.6)') '"unterminated-rv"   ABS 1e-14 ', rv
  write (iu, '(a, i0)')     '"unterminated-iv"   ABS 0     ', iv
end block

! Expressions. A numeric scalar may be given an expression rather than a plain number, so the
! named constants that tao_evaluate_expression knows about, such as pi and c_light, can be used.
!
! This is the second deliberate difference from a namelist read, and again there is nothing to
! compare against since a namelist read rejects all of these. The expected values are asserted
! directly.
!
! Note the absence of a division: "/" ends the group, so "c_light / 1e9" would not reach the
! evaluator. That is not new behavior, it is the namelist file syntax being kept as is.

block
  integer iu2
  type (tao_nml_group_struct) g4

  open (newunit = iu2, file = 'nml_test_expr.nml', status = 'replace')
  write (iu2, '(a)') '&exprgrp'
  write (iu2, '(a)') '  rv = 2*pi              ! named constant, and "*" is a multiply here'
  write (iu2, '(a)') '  arr(1) = c_light * 1e-9'
  write (iu2, '(a)') '  arr(2) = sqrt(m_electron) * 1e-3'
  write (iu2, '(a)') '  iv = 40 + 60           ! would read as just 40 without this'
  write (iu2, '(a)') '  rep = 3*7.5            ! an array, so still a repeat count'
  write (iu2, '(a)') '/'
  close (iu2)

  call reset_vars()
  open (newunit = iu2, file = 'nml_test_expr.nml', status = 'old')
  call tao_nml_group_read (iu2, 'nml_test_expr.nml', 'exprgrp', g4, eof, err, why)
  close (iu2)

  do i = 1, g4%n_item
    call tao_nml_item_ref (g4%item(i), item_ref, err)
    if (err) cycle
    select case (item_ref%head)
    case ('rv');   call tao_nml_value_set (g4%item(i), rv, err)
    case ('iv');   call tao_nml_value_set (g4%item(i), iv, err)
    case ('rep');  call tao_nml_value_set (g4%item(i), rep, err)
    case ('arr');  call set_arr (g4%item(i), arr)
    end select
    if (err) exit
  enddo

  write (iu, '(a, l1, a)')  '"expr-err"          STR       "', err, '"'
  write (iu, '(a, es14.6)') '"expr-twopi"        ABS 1e-14 ', rv
  write (iu, '(a, es14.6)') '"expr-c-light"      ABS 1e-14 ', arr(1)
  write (iu, '(a, es14.6)') '"expr-sqrt"         ABS 1e-14 ', arr(2)
  write (iu, '(a, i0)')     '"expr-integer"      ABS 0     ', iv
  write (iu, '(a, 3es14.6)') '"expr-repeat-kept"  ABS 1e-14 ', rep(1:3)

  ! A value that is neither a number nor an expression is still an error.

  block
    type (tao_nml_item_struct) bad_item
    bad_item%name = 'rv'
    bad_item%value = 'not_a_constant'
    bad_item%file = 'nml_test_expr.nml'
    bad_item%i_line = 1
    call tao_nml_value_set (bad_item, rv, err)
    write (iu, '(a, l1, a)') '"expr-bad-err"      STR       "', err, '"'
  end block
end block

write (iu, '(a, i0)') '"n-mismatch-final"  ABS 0     ', n_mismatch

close (iu)

! Remove the scratch input files so the test leaves nothing behind.

block
  integer iu2
  open (newunit = iu2, file = 'nml_test_input.nml', status = 'old')
  close (iu2, status = 'delete')
  open (newunit = iu2, file = 'nml_test_bad.nml', status = 'old')
  close (iu2, status = 'delete')
  open (newunit = iu2, file = 'nml_test_expr.nml', status = 'old')
  close (iu2, status = 'delete')
end block

!--------------------------------------------------------------------------
contains

subroutine set_arr (item, a)
type (tao_nml_item_struct) item
real(rp), target :: a(:)
type (tao_nml_ref_struct) r
type (tao_ptr_struct) p
character(:), allocatable :: v(:)
character(200) w
integer i1, i2, k, nv, iv
logical e

call tao_nml_item_ref (item, r, e)
if (e) return

if (.not. r%has_sub) then
  call tao_nml_value_set (item, a, e)
  return
endif

call tao_nml_ref_bounds (item, r, 1, size(a), i1, i2, e)
if (e) return

if (i1 == i2) then
  call tao_nml_value_set (item, a(i1), e)
  return
endif

call tao_nml_split_values (item%value, v, nv, e, w)
if (e) return
iv = 0
do k = i1, i2
  iv = iv + 1
  if (iv > nv) exit
  if (v(iv) == '') cycle
  p = tao_ptr_struct()
  p%r => a(k)
  call tao_set_ptr_value (p, v(iv), e, w)
  if (e) return
enddo
end subroutine set_arr

!--------------------------------------------------------------------------

subroutine reset_vars ()
rv = 0
iv = 0
lv = .false.
sv = ''
arr = -1
rng = -1
rep = -1
tabbed = -1
null_mid = -1
dat = tao_datum_input()
end subroutine reset_vars

!--------------------------------------------------------------------------

subroutine cmp_real (name, a, b)
character(*) name
real(rp) a, b
if (a /= b) n_mismatch = n_mismatch + 1
write (iu, '(3a, es14.6, es14.6)') '"', trim(name), '" ABS 1e-14 ', a, b
end subroutine cmp_real

subroutine cmp_int (name, a, b)
character(*) name
integer a, b
if (a /= b) n_mismatch = n_mismatch + 1
write (iu, '(3a, i0, 1x, i0)') '"', trim(name), '" ABS 0 ', a, b
end subroutine cmp_int

subroutine cmp_logic (name, a, b)
character(*) name
logical a, b
if (a .neqv. b) n_mismatch = n_mismatch + 1
write (iu, '(6a)') '"', trim(name), '" STR ', quote(logic_str(a)), ' ', quote(logic_str(b))
end subroutine cmp_logic

subroutine cmp_str (name, a, b)
character(*) name
character(*) a, b
if (a /= b) n_mismatch = n_mismatch + 1
write (iu, '(6a)') '"', trim(name), '" STR ', quote(trim(a)), ' ', quote(trim(b))
end subroutine cmp_str

!--------------------------------------------------------------------------

subroutine write_test_file ()
integer iu2
character(1), parameter :: tb = achar(9)
open (newunit = iu2, file = 'nml_test_input.nml', status = 'replace')
write (iu2, '(a)') '! A leading comment outside any group.'
write (iu2, '(a)') '&testgrp'
write (iu2, '(a)') '  rv = 1.25d0        ! a real with a comment after it'
write (iu2, '(a)') '  iv = 42'
write (iu2, '(a)') '  lv = T'
write (iu2, '(a)') "  sv = 'a string, with a comma'"
write (iu2, '(a)') '  arr = 1.0, 2.0,'
write (iu2, '(a)') '        3.0, 4.0     ! value list continued across lines'
write (iu2, '(a)') '  arr(6) = 9.5       ! subscripted item'
write (iu2, '(a)') '  rng(1:3) = 11.0, 12.0, 13.0   ! subscript range'
write (iu2, '(a)') '  rng(5:) = 25.0                ! open ended range'
write (iu2, '(a)') '  rep = 3*7.5                   ! repeat count'
write (iu2, '(a)') '  tabbed = 1.0'//tb//'2.0'//tb//tb//'3.0   ! tab separated, with a comment'
write (iu2, '(a)') '  null_mid = 1.0,,3.0           ! null value in the middle'
write (iu2, '(a)') "  dat%ele_name = 'Q01W'"
write (iu2, '(a)') "  dat%data_type = 'orbit.x'"
write (iu2, '(a)') "  dat%merit_type = 'target'"
write (iu2, '(a)') '  dat%weight = 3.5'
write (iu2, '(a)') '  dat%good_user = F'
write (iu2, '(a)') '/'
write (iu2, '(a)') ''
write (iu2, '(a)') '&second_group'
write (iu2, '(a)') '  alpha = 7'
write (iu2, '(a)') '  beta = 8'
write (iu2, '(a)') '&end'
close (iu2)
end subroutine write_test_file

end program tao_nml_test
