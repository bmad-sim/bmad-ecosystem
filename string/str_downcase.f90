!+
! Subroutine str_downcase (destination, source)
!
! Subroutine to convert a string to down case.
!
! Input:
!   source -- Character(*): Source string.
!
! Output:
!   destination -- Character(*): Down cased string.
!-

subroutine str_downcase(dst, src)

implicit none

character*(*) dst,src
integer i,s,dlen,slen

dlen=len(dst)
slen=len(src)
i=1
do i=1,dlen
   if (i.le.slen) then
      s=ichar(src(i:i))
      if (s .ge. 65 .and. s .le. 90) then
         s=s+32
      endif
      dst(i:i)=char(s)
   else
      dst(i:i)=' '
   endif
enddo

end subroutine str_downcase


