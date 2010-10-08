!+
! subroutine mat_type (mat, nunit, header, num_form)
!
! subroutine to output matraces to the terminal or to a file
!
! Modules needed:
!   use sim_utils
!
! input: 
!   mat(:,:)  -- Real(rp): Matrix to print out.
!   nunit     -- Integer, optional: Unit for writing:
!                    > 0 output to file only with unit = nunit
!                    = 0 output to terminal only (default)
!                    < 0 output to terminal and file with unit = abs(nunit).
!   header   -- Character(*), optional: Title to print above matrix.
!   num_form -- Character(*), optional: Format for the numbers. Default
!                    is "es13.5" if any |term| > 100 and "f11.6" otherwise.
!-

subroutine mat_type (mat, nunit, header, num_form)

use precision_def

implicit none

integer, optional, intent(in) :: nunit
integer size1, size2, munit, i, j, iu

real(rp) mat(:,:)

character(*), optional :: header, num_form
character(24) format1

!

iu = 0
if (present(nunit)) iu = nunit

size1 = size(mat, 1)
size2 = size(mat, 2)

if (present(num_form)) then
  format1 = num_form
elseif (any(abs(mat) > 100)) then
  write (format1, '(a, i2.2, a)') '(3x, 1p, ', size2, 'es13.5)'
else
  write (format1, '(a, i2.2, a)') '(3x, ', size2, 'f11.6)'
endif

if (iu <= 0) then
  if (present(header)) print *, header
  do i = 1, size1
    print format1, (mat(i, j), j = 1, size2)
  enddo
endif

if (iu /= 0) then
  munit = abs(iu)
  if (present(header)) write (munit, *) header
  do i = 1, size1
    write (munit, format1) (mat(i, j), j = 1, size2)
  enddo
endif

end subroutine
