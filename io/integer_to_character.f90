subroutine integer_to_character(i,c)

  implicit none
  integer i
  integer n
  character*(*) c

  n = len(c)
  c(1:n)=' '
  if(abs(i)<10)write(c(n:n),'(i1)')i
  if(abs(i)>=10 .and.abs(i)<100)write(c(n-1:n),'(i2)')i
  if(abs(i)>=100 .and.abs(i)<1000)write(c(n-2:n),'(i3)')i
  if(abs(i) >= 1000)print '(a,i,a)',' SUBROUTINE INTEGER_TO_CHARACTER: ',i,' >= 1000'
  return
 end
