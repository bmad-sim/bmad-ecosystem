#include "CESR_platform.inc"

module taylor_mod

  use accelerator_struct

  interface assignment (=)
    module procedure equal_real_8_taylor
    module procedure equal_taylor_real_8
    module procedure equal_universal_universal
  end interface

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine equal_real_8_taylor (y8, bmad_taylor)
!
! Subroutine to overload "=" in expressions
!       y8 = bmad_taylor
!
! Modules needed:
!   use accelerator
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Input taylor series array.
!
! Output:
!   y8(:) -- real_8: PTC Taylor series array.
!-

subroutine equal_real_8_taylor (y8, bmad_taylor)

!  use s_tracking

  implicit none

  type (real_8), intent(out) :: y8(:)
  type (taylor_struct), intent(in) :: bmad_taylor(:)

  call taylor_to_real_8 (bmad_taylor, y8, .true.)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine equal_taylor_real_8 (bmad_taylor, y8)
!
! Subroutine to overload "=" in expressions
!       bmad_taylor = y8
!
! Modules needed:
!   use accelerator
!
! Input:
!   y8(:) -- real_8: PTC Taylor series array.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Input taylor series array.
!-

subroutine equal_taylor_real_8 (bmad_taylor, y8)

!  use s_tracking

  implicit none

  type (real_8), intent(in) :: y8(:)
  type (taylor_struct), intent(out) :: bmad_taylor(:)

  call real_8_to_taylor (y8, bmad_taylor, .true.)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_to_taylor (y8, bmad_taylor, switch_z)
!
! Subroutine to convert from a real_8 taylor map in Etienne's PTC 
! to a taylor map in BMAD.
! The conversion can also convert from the PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! to the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Modules needed:
!   use accelerator
!
! Input:
!   y8(6)       -- Real_8: Taylor map.
!   switch_z    -- Logical, optional: If True then switch coordinate 
!                    conventions. Default is True.
!
! Output:
!   bmad_taylor(6) -- Taylor_struct:
!-

subroutine real_8_to_taylor (y8, bmad_taylor, switch_z)

!  use s_tracking

  implicit none

  type (real_8), intent(in) :: y8(:)
  type (taylor_struct), intent(out) :: bmad_taylor(:)
  type (universal_taylor) :: u_t(6)

  integer i

  logical, optional :: switch_z

!

  do i = 1, 6
    u_t(i) = 0  ! nullify
    u_t(i) = y8(i)%t
  enddo

  call universal_to_bmad_taylor (u_t, bmad_taylor, switch_z)

  do i = 1, 6
    u_t(i) = -1  ! deallocate
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_real_8 (bmad_taylor, y8, switch_z)
!
! Subroutine to convert from a taylor map in BMAD to a
! real_8 taylor map in Etienne's PTC.
! The conversion can also convert from the the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
! to PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
!
! Modules needed:
!   use accelerator
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Taylor map.
!   switch_z       -- Logical, optional: If True then switch coordinate 
!                       conventions. Default is True.
!
! Output:
!   y8(6)       -- Real_8: Taylor map.
!-

subroutine taylor_to_real_8 (bmad_taylor, y8, switch_z)

!  use s_tracking

  implicit none

  type (real_8), intent(out) :: y8(:)
  type (taylor_struct), intent(in) :: bmad_taylor(:)
  type (universal_taylor) :: u_t

  integer i, j, ii, n

  logical, optional :: switch_z
  logical switch

! init

  call kill (y8)
  call real_8_init (y8, .true.)

!

  do i = 1, 6

    switch = .true.
    if (present(switch_z)) switch = switch_z

    ii = i
    if (switch) then
      if (i == 5) ii = 6
      if (i == 6) ii = 5
    endif

    n = size(bmad_taylor(i)%term)
    allocate (u_t%n, u_t%nv, u_t%c(n), u_t%j(n,6))
    u_t%n = n
    u_t%nv = 6

    do j = 1, n
      if (switch) then
        u_t%j(j,:) = bmad_taylor(i)%term(j)%exp((/1,2,3,4,6,5/))
        u_t%c(j) = bmad_taylor(i)%term(j)%coef * &
                                      (-1)**bmad_taylor(i)%term(j)%exp(5)
        if (i == 5) u_t%c(j) = -u_t%c(j)
      else
        u_t%j(j,:) = bmad_taylor(i)%term(j)%exp(:)
        u_t%c(j) = bmad_taylor(i)%term(j)%coef
      endif
    enddo

    y8(ii) = u_t
    u_t = -1   ! deallocate
        
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_bmad_to_ptc (vec_bmad, vec_ptc)
!
! Subroutine to convert between BMAD and PTC coordinates.
! PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Input:
!   vec_bmad(6) -- Real(rdef): Input BMAD vector.
!
! Output:
!   vec_ptc(6)  -- Real(dp): Output PTC vector.
!-

subroutine vec_bmad_to_ptc (vec_bmad, vec_ptc)

!  use s_tracking

  implicit none
  
  real(rdef), intent(in)  :: vec_bmad(:)
  real(dp), intent(out)   :: vec_ptc(:)
  real(dp) temp_vec(6)

  temp_vec = vec_bmad((/1,2,3,4,6,5/))
  vec_ptc = temp_vec
  vec_ptc(6) = -vec_ptc(6)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_ptc_to_bmad (vec_ptc, vec_bmad)
!
! Subroutine to convert between BMAD and PTC coordinates.
! PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Input:
!   vec_ptc(6)  -- Real(rdef): Input PTC vector.
!
! Output:
!   vec_bmad(6) -- Real(rdef): Output BMAD vector.
!-

subroutine vec_ptc_to_bmad (vec_ptc, vec_bmad)

  implicit none
  
  real(dp), intent(in)    :: vec_ptc(:)
  real(rdef), intent(out) :: vec_bmad(:)
  real(rdef) temp(6)

  temp = vec_ptc((/1,2,3,4,6,5/))
  vec_bmad = temp
  vec_bmad(5) = -vec_bmad(5)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Elemental subroutine equal_universal_universal (ut1, ut2)
!
! Subroutine to transfer the values in one universal taylor variable to
! another. Note: ut1 needs to have been initialized.
!
! Modules needed:
!   use accelerator
!
! Input:
!   ut2 -- Universal_taylor:
!
! Output:
!   ut1 -- Universal_taylor:
!-

elemental subroutine equal_universal_universal (ut1, ut2)

!  use definition

  implicit none

  type (universal_taylor), intent(inout) :: ut1
  type (universal_taylor), intent(in)    :: ut2

!

  if (associated (ut1%n)) deallocate (ut1%n, ut1%nv, ut1%c, ut1%j)
  allocate (ut1%n, ut1%nv, ut1%c(ut2%n), ut1%j(ut2%n, ut2%nv))

  ut1%n  = ut2%n
  ut1%nv = ut2%nv
  ut1%c  = ut2%c
  ut1%j  = ut2%j

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_init (y, set_taylor)
!
! Subroutine to allocate a PTC real_8 variable.
! The internal kind parameter will be set to 0.
!
! Note: If this variable has been used before make sure you have 
! deallocated using:
!   call kill(y)
!
! Modules needed:
!   use accelerator
!
! Input:
!   y(:)       -- Real_8: 
!   set_taylor -- Logical, optional :: If present and True then make
!                   y the identity taylor series (kind = 2).
!
! Output:
!   y(:) -- Real_8:
!-

subroutine real_8_init (y, set_taylor)

!  use s_tracking

  implicit none
  
  type (real_8) :: y(:)
  real(dp) :: x(6) = (/ 0, 0, 0, 0, 0, 0 /)

  logical, optional :: set_taylor

!

  call alloc(y)
  y = bmad_com%real_8_map_init

  if (present(set_taylor)) then
    if (set_taylor) y = x   ! converts y to taylor (kind = 2)
  endif

end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor, switch_z)
!
! Subroutine to convert from a universal_taylor map in Etienne's PTC 
! to a taylor map in BMAD.
! The conversion can also convert from the PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! to the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Modules needed:
!   use accelerator
!
! Input:
!   u_taylor(6) -- Universal_taylor: Universal_taylor map.
!   switch_z    -- Logical, optional: If True then switch coordinate 
!                    conventions. Default is True.
!
! Output:
!   bmad_taylor(6)   -- Taylor_struct:
!-

Subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor, switch_z)

!  use definition

  implicit none

  type (universal_taylor), intent(in) :: u_taylor(:)
  type (taylor_struct) :: bmad_taylor(:)

  integer i, j, k, ii, n

  logical, optional :: switch_z
  logical switch

! Remember to suppress any terms that have a zero coef.  

  do i = 1, 6

    switch = .true.
    if (present(switch_z)) switch = switch_z

    ii = i
    if (switch) then
      if (i == 5) ii = 6
      if (i == 6) ii = 5
    endif

    if (associated(bmad_taylor(i)%term)) deallocate(bmad_taylor(i)%term)

    n = count(u_taylor(ii)%c(:) /= 0)
    allocate(bmad_taylor(i)%term(n))

    k = 0
    do j = 1, u_taylor(ii)%n
      if (u_taylor(ii)%c(j) == 0) cycle
      k = k + 1
      if (switch) then
        bmad_taylor(i)%term(k)%exp  = u_taylor(ii)%j(j, (/1,2,3,4,6,5/))
        bmad_taylor(i)%term(k)%coef = u_taylor(ii)%c(j)
        bmad_taylor(i)%term(k)%coef = bmad_taylor(i)%term(k)%coef * &
                                      (-1)**bmad_taylor(i)%term(k)%exp(5)
        if (i == 5) bmad_taylor(i)%term(k)%coef = -bmad_taylor(i)%term(k)%coef
      else
        bmad_taylor(i)%term(k)%exp  = u_taylor(ii)%j(j,:)
        bmad_taylor(i)%term(k)%coef = u_taylor(ii)%c(j)
      endif
    enddo

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_real_8 (y1, y2, y3)
!
! Subroutine to concatinate two real_8 taylor series.
! This subroutine assumes that y1, y2, and y3 have been allocated.
!
! Modules needed:
!   use accelerator
!
! Input:
!   y1(6) -- real_8: Input.
!   y2(6) -- real_8: Input.
!
! Output
!   y3(6) -- real_8: Concatinated output.
!-

subroutine concat_real_8 (y1, y2, y3)

  implicit none

  type (real_8), intent(in) :: y1(:), y2(:)
  type (real_8), intent(out) :: y3(:)
  type (damap) da1, da2, da3

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

  call alloc(da1)
  call alloc(da2)
  call alloc(da3)

! concat

  da1 = y1
  da2 = y2

  da3 = da1 .o. da2  ! concat with constant terms
  
  y3 = da3

! kill temp vars

  call kill (da1)
  call kill (da2)
  call kill (da3)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_genfield (bmad_taylor, gen_field, r0)
!
! Subroutine to construct a genfield (partially inverted map) from a taylor
! map.
! Note: The constant terms of the taylor map are removed in the process.
! Note: The genfield uses PTC coordinates.
!
! Moudules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Input taylor map.
!
! Output:
!   gen_field      -- Genfield: Output partially inverted map.
!   r0(6)          -- Real(rdef): The constant part of the bmad_taylor map
!-

subroutine taylor_to_genfield (bmad_taylor, gen_field, r0)

  implicit none

  type (taylor_struct), intent(in) :: bmad_taylor(6)
  type (genfield), intent(inout) :: gen_field
  type (taylor_struct) taylor_(6)
  type (damap) da_map
  type (real_8) y(6)

  real(rdef), intent(out) :: r0(6)

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! Remove constant terms from the taylor map first. This is probably
! not needed but we do it to make sure everything is alright.
! Also remove terms that have higher order then bmad_com%taylor_order

  call remove_constant_taylor (bmad_taylor, taylor_, r0, .true.)

! allocate pointers

  call alloc (gen_field)
  call alloc (da_map)
  call alloc (y)

! calculate the gen_field

  y = taylor_
  da_map = y
  gen_field = da_map

! cleanup

  call kill (da_map)
  call kill (y)
  call kill_taylor (taylor_)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine remove_constant_taylor (taylor_in, taylor_out, r0, &
!                                                 remove_higher_order_terms)
!
! Subroutine to remove the constant part of a taylor series.
! Optionally terms that are higher order than bmad_com%taylor_order can
! be removed.
!
! Note: It is assumed that taylor_out has been deallocated before the call to
! this routine. Calling this routine with the first two actual arguments the
! same is prohibited.
!
! Moudules needed:
!   use bmad
!
! Input:
!  taylor_in(6) -- Taylor_struct: Input taylor map.
!  remove_higher_order_terms -- Logical: If True then terms that are higher
!                               order than bmad_com%taylor_order are removed.
!
! Output:
!   taylor_out(6)  -- Taylor_struct: Taylor with constant terms removed.
!   r0(6)          -- Real(rdef): The constant part of the taylor map
!-

subroutine remove_constant_taylor (taylor_in, taylor_out, r0, &
                                                 remove_higher_order_terms)

  implicit none

  type (taylor_struct), intent(in) :: taylor_in(:)
  type (taylor_struct) taylor_out(:)

  real(rdef), intent(out) :: r0(:)

  integer i, j, n, nn, ss

  logical, intent(in) :: remove_higher_order_terms

!

  r0 = 0

  do i = 1, 6

    n = size(taylor_in(i)%term)

    do j = 1, size(taylor_in(i)%term)
      if (all(taylor_in(i)%term(j)%exp == 0)) then
        n = n - 1
        r0(i) = taylor_in(i)%term(j)%coef
      endif
      if (remove_higher_order_terms) then
        if (sum(taylor_in(i)%term(j)%exp) > bmad_com%taylor_order) n = n - 1
      endif
    enddo

    allocate (taylor_out(i)%term(n))

    nn = 0
    do j = 1, size(taylor_in(i)%term)
      ss = sum(taylor_in(i)%term(j)%exp)
      if (ss == 0 .or. (remove_higher_order_terms .and. &
                                            ss > bmad_com%taylor_order)) cycle
      nn = nn + 1
      taylor_out(i)%term(nn) = taylor_in(i)%term(j)
    enddo

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_inverse (taylor, taylor_inv)
!
! Subroutine to invert a taylor map.
!
! Moudules needed:
!   use bmad
!
! Input:
!   taylor(6)     -- Taylor_struct: Input taylor map.
!
! Output:
!   taylor_inv(6) -- Taylor_struct: Inverted taylor map.
!-

subroutine taylor_inverse (taylor, taylor_inv)

  implicit none

  type (taylor_struct), intent(in) :: taylor(:)
  type (taylor_struct), intent(out) :: taylor_inv(:)
  type (taylor_struct) tlr(6)
  type (real_8) y(6), yc(6)
  type (damap) da

  real(rdef) r0(6)
  real(8) r8(6)

! The inverse operation of PTC ignores constant terms so we have to take
! them out and then put them back in.

  call remove_constant_taylor (taylor, tlr, r0, .true.)

  call alloc(da)
  call alloc(y)

! compute inverse

  y = tlr
  da = y
  da = da**(-1)
  y = da

! put constant terms back in

  if (any(r0 /= 0)) then
    call real_8_init(yc)
    yc = real(r0, 8)
    call concat_real_8 (y, yc, y)
    call kill (yc)
  endif

! transfer inverse to taylor_inv

  taylor_inv = y

! clean up

  call kill (da)
  call kill (y)
  call kill_taylor (tlr)

end subroutine
end module
