!   match=str_find_first_substring(str,where,which,'str1','str2'..)
!   find where one of str1,str2... found LEFTMOST in 'str'. 
!   'which' tells which substr matched.
!   match=true. if some substring matches.
!   this implementation limited to 13 substrings
!    Need module for optional arguments to be tested present
module str_find_first_substring_module
 contains

   function str_find_first_substring &
(line,where,which,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13) result (is_match)
  implicit none
  character(*), intent(in):: line
  character(*), intent(in), optional:: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
  integer ws(13)
  integer, intent(out):: which,where
  logical is_match,pres(13)
  integer i
  is_match=.false. ; pres(13)=present(s13)
  pres(1)=present(s1) ; pres(2)=present(s2) ; pres(3)=present(s3)
  pres(4)=present(s4) ; pres(5)=present(s5) ; pres(6)=present(s6)
  pres(7)=present(s7) ; pres(8)=present(s8) ; pres(9)=present(s9)
  pres(10)=present(s10) ; pres(11)=present(s11) ; pres(12)=present(s12)
  where=9998 ; which=0
  ws=9999    !clear locs
  if(pres(1)) ws(1)=index(line,trim(s1))
  if(pres(2)) ws(2)=index(line,trim(s2))
  if(pres(3)) ws(3)=index(line,trim(s3))
  if(pres(4)) ws(4)=index(line,trim(s4))
  if(pres(5)) ws(5)=index(line,trim(s5))
  if(pres(6)) ws(6)=index(line,trim(s6))
  if(pres(7)) ws(7)=index(line,trim(s7))
  if(pres(8)) ws(8)=index(line,trim(s8))
  if(pres(9)) ws(9)=index(line,trim(s9))
  if(pres(10)) ws(10)=index(line,trim(s10))
  if(pres(11)) ws(11)=index(line,trim(s11))
  if(pres(12)) ws(12)=index(line,trim(s12))
  if(pres(13)) ws(13)=index(line,trim(s13))
  do i=1,13
   if((ws(i).gt.0).and.(ws(i).lt.where)) then
    is_match=.true. ; where=ws(i) ; which=i
   endif
  enddo 
  if (.not.is_match) where=0
  return
end function str_find_first_substring
end module
