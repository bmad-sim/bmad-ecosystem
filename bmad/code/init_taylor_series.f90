!+
! Subroutine init_taylor_series (bmad_taylor, n_term, save_old)
!
! Subroutine to initialize or extend a Bmad Taylor series (6 of these series make
! a Taylor map). Note: This routine does not zero the terms.
!
! Input:
!   bmad_taylor -- taylor_struct: Old structure.
!   n_term      -- integer: Number of terms to allocate. 
!                    n_term < 0 => bmad_taylor%term pointer will be disassociated.
!   save_old    -- logical, optional: If True then save any old terms and ref orbit when
!                    bmad_taylor is resized. If False zero the ref orbit. Default is False.
!
! Output:
!   bmad_taylor -- Taylor_struct: Initalized structure.
!-

subroutine init_taylor_series (bmad_taylor, n_term, save_old)

use bmad_routine_interface, dummy => init_taylor_series

implicit none

type (taylor_struct) bmad_taylor
type (taylor_term_struct), pointer :: term(:)
integer n_term
integer n
logical, optional :: save_old

!

if (.not. logic_option (.false., save_old)) bmad_taylor%ref = 0

if (n_term < 0) then
  if (associated(bmad_taylor%term)) deallocate(bmad_taylor%term)
  return
endif

if (.not. associated (bmad_taylor%term)) then
  allocate (bmad_taylor%term(n_term))
  return
endif

if (size(bmad_taylor%term) == n_term) return

!

if (logic_option (.false., save_old) .and. n_term > 0 .and. size(bmad_taylor%term) > 0) then
  n = min (n_term, size(bmad_taylor%term))
  term => bmad_taylor%term
  allocate (bmad_taylor%term(n_term))
  bmad_taylor%term(1:n) = term(1:n)
  deallocate (term)

else
  deallocate (bmad_taylor%term)
  allocate (bmad_taylor%term(n_term))
endif

end subroutine init_taylor_series
