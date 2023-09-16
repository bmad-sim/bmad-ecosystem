!+
! subroutine mat_type (mat, nunit, header, num_form, lines, n_lines)
!
! subroutine to output matraces to the terminal or to a file
! When writing to a file, it is assumed that the file has been opened.
!
! input: 
!   mat(:,:)  -- real(rp): Matrix to print out.
!   nunit     -- integer, optional: Unit for writing:
!                    > 0 output to file only with unit = nunit.
!                    = 0 output to terminal only (default)
!                    < 0 output to terminal and file with unit = abs(nunit).
!   header    -- character(*), optional: Title to print above matrix.
!   num_form  -- character(*), optional: Format for the numbers. 
!                    Default is "3x, NNNf11.6" where NNN is the matrix row size. 
!                    if any |term| > 100 then "3x, NNNes14.5" will be the default.
!   lines(:)  -- character(*), optional: Character array to hold the output.
!                   If present, output to the terminal will not be done.
!   n_lines   -- integer, optional: Number of lines in lines(:) that hold valid output.
!                   n_lines bust be present if lines(:) is.
!-

subroutine mat_type (mat, nunit, header, num_form, lines, n_lines)

use sim_utils, dummy => mat_type

implicit none

integer, optional :: nunit, n_lines
integer size1, size2, munit, i, j, iu, ios

real(rp) mat(:,:)

character(*), optional :: header, num_form, lines(:)
character(400) line
character(24) format1, form

!

iu = 0
if (present(nunit)) iu = nunit
if (present(n_lines)) n_lines = 0

size1 = size(mat, 1)
size2 = size(mat, 2)

if (present(num_form)) then
  form = '(' // trim(num_form) // ')'
  do i = 1, size1
    write (line, form, iostat = ios) mat(i,:)
    if (ios /= 0) exit
  enddo

  if (ios == 0) then
    format1 = form
  else
    format1 = '(3x, 100es15.5)'
  endif

elseif (any(abs(mat) > 100)) then
  format1 = '(3x, 100es15.5)'

else
  format1 = '(3x, 100f11.6)'
endif

if (iu <= 0) then
  if (present(lines)) then
    if (present(header)) then
      n_lines=n_lines+1; write (lines(n_lines), '(a)') header
    endif

    do i = 1, size1
      n_lines=n_lines+1; call mat_line(lines(n_lines), format1, mat(i,:))
    enddo

  else
    if (present(header)) print '(a)', header
    do i = 1, size1
      call mat_line(line, format1, mat(i,:))
      print '(a)', trim(line)
    enddo
  endif
endif

if (iu /= 0) then
  munit = abs(iu)
  if (present(header)) write (munit, *) header
  do i = 1, size1
    call mat_line(line, format1, mat(i,:))
    write (munit, '(a)') trim(line)
  enddo
endif

!------------------------------------------
contains

subroutine mat_line(line, format1, row)

character(*) line, format1
real(rp) row(:)
integer ix, ix2, nw, p, ios

!

write (line, format1, iostat = ios) row

end subroutine mat_line

end subroutine
