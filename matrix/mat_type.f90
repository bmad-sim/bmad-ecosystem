!+
! subroutine mat_type (mat, nunit, header, num_form, lines, n_lines)
!
! subroutine to output matraces to the terminal or to a file
!
! Modules needed:
!   use sim_utils
!
! input: 
!   mat(:,:)  -- real(rp): Matrix to print out.
!   nunit     -- integer, optional: Unit for writing:
!                    > 0 output to file only with unit = nunit
!                    = 0 output to terminal only (default)
!                    < 0 output to terminal and file with unit = abs(nunit).
!   header    -- character(*), optional: Title to print above matrix.
!   num_form  -- character(*), optional: Format for the numbers. 
!                    Default is "(3x, NNNf11.6)" where NNN is the matrix row size. 
!                    if any |term| > 100 then "(3x, NNNes14.5)" will be the default.
!   lines(:)  -- character(*), optional: Character array to hold the output.
!                   If present, output to the terminal will not be done.
!   n_lines   -- integer, optional: Number of lines in lines(:) that hold valid output.
!                   n_lines bust be present if lines(:) is.
!-

subroutine mat_type (mat, nunit, header, num_form, lines, n_lines)

use precision_def

implicit none

integer, optional :: nunit, n_lines
integer size1, size2, munit, i, j, iu

real(rp) mat(:,:)

character(*), optional :: header, num_form, lines(:)
character(24) format1

!

iu = 0
if (present(nunit)) iu = nunit
if (present(n_lines)) n_lines = 0

size1 = size(mat, 1)
size2 = size(mat, 2)

if (present(num_form)) then
  format1 = num_form
elseif (any(abs(mat) > 100)) then
  write (format1, '(a, i2.2, a)') '(3x, 1p, ', size2, 'es15.5)'
else
  write (format1, '(a, i2.2, a)') '(3x, ', size2, 'f11.6)'
endif

if (iu <= 0) then
  if (present(lines)) then
    if (present(header)) then
      n_lines = n_lines + 1; write (lines(n_lines), '(a)') header
    do i = 1, size1
      n_lines = n_lines + 1; write (lines(n_lines), format1) (mat(i, j), j = 1, size2)
    enddo
  endif

  else
    if (present(header)) print '(a)', header
    do i = 1, size1
      print format1, (mat(i, j), j = 1, size2)
    enddo
  endif
endif

if (iu /= 0) then
  munit = abs(iu)
  if (present(header)) write (munit, *) header
  do i = 1, size1
    write (munit, format1) (mat(i, j), j = 1, size2)
  enddo
endif

end subroutine
