module plot_frequency_mod

use bmad

type column_struct 
  real row(15)
end type

contains

!---------------------------------------------------

subroutine read_scan (file_name, column, size) 

implicit none

type(column_struct), allocatable, intent(out) :: column(:)
character(200) file_name
character(200) line
character(20) word
integer i, n, ix, j, size
integer k
integer ios

open(unit=2, file=file_name, status='OLD')
ios=0
i=0
do while(ios >= 0)
  read(2,'(a200)', iostat=ios)line
  if(ios /= 0)cycle
  if(index(line,'!') /= 0) cycle
  i = i+1
end do

if(.not. allocated(column))allocate(column(1:i))
size = i

rewind(2)
i=0
ios=0
do while(ios >= 0)
  read(2,'(a200)', iostat=ios)line
  if(ios /= 0)cycle
  if(index(line,'!') /= 0) cycle
  i = i+1

  ix=0
  do j=1,12
    call string_trim(line(ix+1:), line, ix)
    word=line(1:ix)
    if(index(word,'NaN') /= 0)then
      column(i)%row(j) = 0.
     else
      read(word,*) column(i)%row(j)
    endif
  end do
end do

close(unit=2)

end subroutine
   
end module
