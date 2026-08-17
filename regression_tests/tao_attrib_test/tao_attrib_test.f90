!+
! Program tao_attrib_test
!
! Test of the generated structure component resolvers in tao_attrib_resolve_mod.
!
! The resolvers map a component name string onto a pointer to that component. They replace the
! use of a Fortran namelist read on a scratch file for setting structure components from a
! string. This test checks that:
!   * Scalars, array elements, whole arrays, and nested structures all resolve.
!   * Values written through the returned pointer land in the structure.
!   * Names are case insensitive and tolerate surrounding whitespace.
!   * Malformed and unknown references produce an error with a useful explanation
!     rather than being silently ignored (which is what a namelist read does).
!
! Note: For the positive tests, a resolution failure would leave the component at its default
! value and so would be caught by the value check. The "no-errors" line is a further check that
! none of the positive cases set the error flag.
!-

program tao_attrib_test

use tao_attrib_resolve_mod

implicit none

type (bmad_common_struct) bc
type (beam_init_struct) bi
type (tao_datum_input) datum
type (tao_ptr_struct) ptr

integer iu
logical err, any_err
character(200) why

!

open (newunit = iu, file = 'output.now', recl = 300)
any_err = .false.

! Scalar components of each intrinsic type.

call bc_test ('abs_tol_adaptive_tracking');  if (.not. err) ptr%r = 1.5_rp
write (iu, '(a, es14.6)')  '"real-scalar"       ABS 1e-14 ', bc%abs_tol_adaptive_tracking

call bc_test ('taylor_order');               if (.not. err) ptr%i = 7
write (iu, '(a, i0)')      '"integer-scalar"    ABS 0     ', bc%taylor_order

call bc_test ('sr_wakes_on');                if (.not. err) ptr%l = .false.
write (iu, '(3a)')         '"logical-scalar"    STR       ', quote(logic_str(bc%sr_wakes_on))

call bi_test ('species');                    if (.not. err) ptr%str = 'proton'
write (iu, '(2a)')         '"char-scalar"       STR       ', quote(trim(bi%species))

! Array element and whole array.

call bc_test ('d_orb(3)');                   if (.not. err) ptr%r = 2.5_rp
write (iu, '(a, es14.6)')  '"array-element"     ABS 1e-14 ', bc%d_orb(3)

call bc_test ('d_orb')
write (iu, '(a, i0)')      '"whole-array-size"  ABS 0     ', size_or_zero(ptr)

call bi_test ('distribution_type(2)');       if (.not. err) ptr%str = 'GRID'
write (iu, '(2a)')         '"char-array-elem"   STR       ', quote(trim(bi%distribution_type(2)))

! Nested structures: scalar parent and array parent.

call bi_test ('kv%n_i2');                    if (.not. err) ptr%i = 11
write (iu, '(a, i0)')      '"nested-scalar"     ABS 0     ', bi%kv%n_I2

call bi_test ('kv%part_per_phi(2)');         if (.not. err) ptr%i = 13
write (iu, '(a, i0)')      '"nested-array"      ABS 0     ', bi%kv%part_per_phi(2)

call bi_test ('ellipse(2)%sigma_cutoff');    if (.not. err) ptr%r = 3.5_rp
write (iu, '(a, es14.6)')  '"struct-array"      ABS 1e-14 ', bi%ellipse(2)%sigma_cutoff

! Case insensitivity and surrounding whitespace.

call bc_test ('  ABS_TOL_TRACKING  ');       if (.not. err) ptr%r = 4.5_rp
write (iu, '(a, es14.6)')  '"case-and-space"    ABS 1e-14 ', bc%abs_tol_tracking

! A Tao init file input structure.

ptr = tao_ptr_struct()
call tao_res_tao_datum_input (datum, 'ele_name', ptr, err, why)
any_err = (any_err .or. err)
if (.not. err) ptr%str = 'Q01W'
write (iu, '(2a)')         '"input-struct"      STR       ', quote(trim(datum%ele_name))

write (iu, '(3a)')         '"no-errors"         STR       ', quote(logic_str(any_err))

! Error cases. Each must set err and give an explanation.
! From here on, err being True is the expected result so any_err is not accumulated.

call bc_test ('no_such_component')
write (iu, '(4a)') '"err-unknown"       STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bc_test ('taylor_order(2)')
write (iu, '(4a)') '"err-sub-on-scalar" STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bc_test ('d_orb(99)')
write (iu, '(4a)') '"err-sub-range"     STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bi_test ('kv')
write (iu, '(4a)') '"err-struct-bare"   STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bc_test ('taylor_order%foo')
write (iu, '(4a)') '"err-sub-of-scalar" STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bi_test ('ellipse%sigma_cutoff')
write (iu, '(4a)') '"err-no-subscript"  STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bc_test ('')
write (iu, '(4a)') '"err-blank"         STR ', quote(logic_str(err)), ' ', quote(trim(why))

call bc_test ('d_orb(x)')
write (iu, '(4a)') '"err-bad-subscript" STR ', quote(logic_str(err)), ' ', quote(trim(why))

close (iu)

!--------------------------------------------------------------------------
contains

! err is accumulated into any_err. The "no-errors" line is written before the error cases
! begin, so only the positive cases contribute to it.

subroutine bc_test (name)
character(*) name
ptr = tao_ptr_struct()
call tao_res_bmad_common_struct (bc, name, ptr, err, why)
any_err = (any_err .or. err)
end subroutine bc_test

subroutine bi_test (name)
character(*) name
ptr = tao_ptr_struct()
call tao_res_beam_init_struct (bi, name, ptr, err, why)
any_err = (any_err .or. err)
end subroutine bi_test

function size_or_zero (ptr) result (n)
type (tao_ptr_struct) ptr
integer n
n = 0
if (associated(ptr%r1)) n = size(ptr%r1)
end function size_or_zero

end program tao_attrib_test
