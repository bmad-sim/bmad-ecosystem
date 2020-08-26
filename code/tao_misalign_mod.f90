module tao_misalign_mod

use tao_data_and_eval_mod

type (ele_struct), target, save :: zero_ele

! with respect to what?
  integer, parameter :: wrt_beam$ = 1, wrt_model$ = 2, wrt_design$ = 3
  integer, parameter ::  wrt_survey$ = 4, wrt_i_beam$ = 5

contains

!------------------------------------------------------------------------------
!+
! Subroutine tao_misalign (ele_type, where, ele_attribute, misalign_value, wrt_what)
!
! Routine to introduce gaussian distributed random misalignments to 
! specified elements of a given type in the model lattice in the specified universe.
!
! If ele_type begins with "*@" choose all universes. If ele_type
! begin with "n@" then choose universe n. Otherwise the viewed universe is used.
!
! If misalign_value is prepended by 'x' then the misalignment value will be a
! relative misalignment. Otherwise, it's an absolute rms value.
!
! To deal with Woodleyfied lattices this routine will check if the element to
! misalign is butted up against another of the same, if so, both are treated as
! a single element.
!
! In the special case where you are changing a sbend strength then use
! dg. However, if a relative error is specified it will be relative to 'g'.
!
!  Input:
!      ele_type        -- character(*) : Element type to be misaligned
!      where           -- character(*) : Ele index to search
!                                        '*' means search all elements
!                                        This isn't used with wrt_beam$
!      ele_attrib      -- character(*) : Attribute to be misaligned
!      misalign_value  -- character(*) : Spread of gaussian errors. 
!      wrt             -- character(*): misalign wrt what?
!                         If wrt_i_beam then only misalign if on an I-Beam
!                         If wrt_survey then only misalign if NOT on an I-Beam
!
!
!  Output:
!      s%u(i)%model    -- lat_struct  : Model lattice with gaussian 
!                                        errors introduced on current values.
!                                        i being specified universe
!-
!------------------------------------------------------------------------------

subroutine tao_misalign (wrt, ele_type, where, ele_attrib, misalign_value)

use random_mod
use attribute_mod, only: attribute_free
use bookkeeper_mod, only: set_flags_for_changed_attribute

implicit none

type (ele_struct), pointer :: ele => null()
type (ele_struct), pointer :: ele2 => null()
type (ele_struct), pointer :: ref_ele => null()
type (coord_struct), allocatable, save :: orb(:)
type (lat_struct), pointer :: lat
type (all_pointer_struct) a_ptr

character(*) :: ele_type, ele_attrib, misalign_value, where, wrt
character(20) :: r_name = 'tao_misalign'

real(rp) misalign_value_num, gauss_err, pitch_offset, ref_value
real(rp) theta, length, Ecav, E0
real(rp), pointer :: value

integer i, j, k, ix_attrib, universe, ix_current, wrt_what
integer time_values(8), seed(1)
integer :: ix_key(2), key_next, num_loc

logical err, found, rel_error, found_double, rel_sbend_error_flag
logical, allocatable, save :: which_univ(:)
logical, allocatable, save :: action_logic(:)

! random number seeded dulat tao initialization

select case (wrt)
case ('wrt_model')
  wrt_what = wrt_model$
case ('wrt_design')
  wrt_what = wrt_design$
case ('wrt_survey')
  wrt_what = wrt_survey$
  call init_ele (zero_ele)
end select

call upcase_string(ele_type)
ele_type = ele_type(1:len(trim(ele_type)))
call tao_pick_universe(ele_type, ele_type, which_univ, err)
if (err) return

call match_word (ele_type, key_name, ix_key(1))
if (ix_key(1) .le. 0) then
  call out_io (s_error$, r_name, "Error matching key name")
  return
endif

!Find an element in the lattice of this key type then find the attribute in
!this element
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. which_univ(i)) cycle
  do j = 1, s%u(i)%model%lat%n_ele_max
    if (s%u(i)%model%lat%ele(j)%key .eq. ix_key(1)) then
      ele => s%u(i)%model%lat%ele(j)
      exit
    endif
  enddo
enddo
if (.not. associated(ele)) then
  call out_io (s_error$, r_name, "Could not find such an element in the specified lattices")
  return
endif

call upcase_string(ele_attrib)
call pointer_to_attribute (ele, ele_attrib, .false., a_ptr, err, .true., ix_attrib)
if (err .or. .not. associated(a_ptr%r)) then
  call out_io (s_error$, r_name, "Error matching attribute name")
  return
endif

rel_error = .false.
rel_sbend_error_flag = .false.
if (index(misalign_value, 'x') .eq. 1) then
  ! realitive error
  rel_error = .true.
  call tao_to_real (misalign_value(2:), misalign_value_num, err)
  if (err) return
  if (ix_key(1) .eq. sbend$ .and. ix_attrib .eq. dg$) rel_sbend_error_flag = .true.
else
  call tao_to_real (misalign_value, misalign_value_num, err)
  if (err) return
endif

num_loc = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (num_loc .lt. s%u(i)%model%lat%n_ele_max) &
       num_loc = s%u(i)%model%lat%n_ele_max
enddo
allocate (action_logic(0:num_loc))
call location_decode (where, action_logic, 0, num_loc)
if (num_loc .eq. -1) goto 999

found_double = .false.

! set misalignment
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. which_univ(i)) cycle
  lat => s%u(i)%model%lat
  do j = 1, lat%n_ele_max
    if (found_double) then
      found_double = .false.
      cycle
    endif
    if (lat%ele(j)%key .eq. ix_key(1) .and. action_logic(j)) then
      call find_ele_and_ref_ele (wrt_what, i, j, ref_ele, ele)
      if (.not. attribute_free (j, 0, ele_attrib, lat)) cycle
      if (rel_error) then
        if (rel_sbend_error_flag) then
          ele%value(ix_attrib) = ref_ele%value(g$) * gauss_err() * misalign_value_num
        else
          ele%value(ix_attrib) = ref_ele%value(ix_attrib) * (1 + gauss_err() * misalign_value_num)
        endif
      else
        ele%value(ix_attrib) = ref_ele%value(ix_attrib) + gauss_err() * misalign_value_num
      endif
      call set_flags_for_changed_attribute (ele, ele%value(ix_attrib))

      if (j .eq. lat%n_ele_max) exit

      if (lat%ele(j+1)%key .eq. ix_key(1) .and. action_logic(j)) then
        call find_ele_and_ref_ele (wrt_what, i, j+1, ref_ele, ele2)
        ! Found Woodley double
        ele2%value(ix_attrib) = ele%value(ix_attrib)
        found_double = .true.
        call set_flags_for_changed_attribute (ele2, ele2%value(ix_attrib))
      endif
      
    end if
  end do

  s%u(i)%calc%lattice = .true.
end do


999 continue

deallocate (action_logic)

end subroutine tao_misalign

!------------------------------------------------------------------------------
subroutine find_ele_and_ref_ele (wrt_what, ix_uni, ix_ele, ref_ele, ele)

implicit none

type (ele_struct), pointer :: ref_ele
type (ele_struct), pointer :: ele

integer wrt_what, ix_ele, ix_uni

select case (wrt_what)
case (wrt_model$)
  ref_ele => s%u(ix_uni)%model%lat%ele(ix_ele)
  ele => s%u(ix_uni)%model%lat%ele(ix_ele)
case (wrt_design$)
  ref_ele => s%u(ix_uni)%design%lat%ele(ix_ele)
  ele => s%u(ix_uni)%model%lat%ele(ix_ele)
case (wrt_survey$)
  ref_ele => zero_ele
  ele => s%u(ix_uni)%model%lat%ele(ix_ele)
end select

end subroutine find_ele_and_ref_ele

end module tao_misalign_mod

! This must be outide the module for the Lahey compiler not to complain
!------------------------------------------------------------------------------
! Returns random numbers of a Gaussian distribution with RMS=1.

real(rp) function gauss_err () result (error)

use precision_def
use random_mod

implicit none
 
call ran_gauss(error)

end function gauss_err

