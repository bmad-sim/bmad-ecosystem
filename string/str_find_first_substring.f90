!   match=str_find_first_substring(str,where,which,'str1','str2'..)
!   find where one of str1,str2... found in 'str'. 
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
  character(*), intent(in):: s1
  character(*), intent(in), optional:: s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
  integer, intent(out):: which,where
  logical is_match,pres
  integer i
  is_match=.false.
  do i=1,13
   pres=.false.
   select case(i)
    case(1)
      pres=.true. ;where = index(line,trim( s1))
    case(2) 
     if(present(s2)) pres=.true.; if(pres) where=index(line,trim(s2))
    case(3) 
     if(present(s3)) pres=.true.; if(pres) where=index(line,trim(s3))
    case(4) 
     if(present(s4)) pres=.true.; if(pres) where=index(line,trim(s4))
    case(5) 
     if(present(s5)) pres=.true.; if(pres) where=index(line,trim(s5))
    case(6) 
     if(present(s6)) pres=.true.; if(pres) where=index(line,trim(s6))
    case(7) 
     if(present(s7)) pres=.true.; if(pres) where=index(line,trim(s7))
    case(8) 
     if(present(s8)) pres=.true.; if(pres) where=index(line,trim(s8))
    case(9) 
     if(present(s9)) pres=.true.; if(pres) where=index(line,trim(s9))
    case(10) 
     if(present(s10)) pres=.true.; if(pres) where=index(line,trim(s10))
    case(11) 
     if(present(s11)) pres=.true.; if(pres) where=index(line,trim(s11))
    case(12) 
     if(present(s12)) pres=.true.; if(pres) where=index(line,trim(s12))
    case(13) 
     if(present(s13)) pres=.true.; if(pres) where=index(line,trim(s13))
   end select
   if(.not.pres)then
    where=0 ; return
   elseif(where.gt.0) then
    which=i ; is_match=.true.;  return
   endif 
  enddo 
end function str_find_first_substring
end module
