integer function numtest(string)
implicit none
integer str_find_first_in_set,str_find_first_not_in_set
character octset*8,intset*12,hexset*22,fltset*13,expset*15
parameter (octset='01234567') !use these to test legal
parameter (intset='-+0123456789')
parameter (hexset='0123456789ABCDEFabcdef')
parameter (fltset='-+0123456789.')
parameter (expset='-+0123456789.Ee')
character*(*) string
integer length,bad,first
length=len(string)  !length of user (sub)string
numtest=-1    !init to bad num
if(string(1:1) == '"') then !might be octal
 bad=str_find_first_not_in_set(string(2:length),octset)
 if(bad == 0) numtest=2 !flag as good attempt at oct
elseif((string(1:1) == 'H').or.(string(1:1) == 'h')) then
 bad=str_find_first_not_in_set(string(2:length),hexset)
 numtest=2    !assume hex (integer) unless
 if(bad /= 0) numtest=1 !flag as plain string (non-hex chars)
else                    !decimal or ascii
 bad=str_find_first_not_in_set(string(1:length),intset)
 numtest=2    !integer unless
 if(bad > 0) then  !last chance is as real (float)
  numtest=3             !presume float
  bad=str_find_first_not_in_set(string(1:length),fltset)
  if(bad > 0)  then  !some non-numeric, check first char
   numtest=4
   bad=str_find_first_not_in_set(string(1:length),expset)
    if(bad > 0) then
     numtest=1    !assume ascii with some nums
     first=str_find_first_in_set(string(1:length),intset)
     if(first == 1) numtest=-1  !started with num, cannot change
    endif
  endif
 endif
endif
return
end
