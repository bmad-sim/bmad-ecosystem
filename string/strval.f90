real function strval(string,type,ng)
!  format convert ascii string to real num using supplied type
implicit none
character*(*) string
integer length,type,ival
logical ng
ng=.true.      !if conversion fails
length=len(string)  !length of user (sub)string
strval=0.0
if(type < 0)  return  !not a number
if(type == 2) then  !integer
 if((string(1:1) == 'H').or.(string(1:1) == 'h')) then
  read(string(2:length),16,err=99) ival !hex
16    format(z)
 elseif(string(1:1) == '"') then
  read(string(2:length),8,err=99) ival !octal
8    format(o)
 else
  read(string(1:length),10,err=99) ival !base ten integer
10    format(i)
 endif
 ng=.false.
 strval=ival            !return uniformly as float (=>use nint)
elseif(type == 3) then  !real
 read(string(1:length),1,err=99) strval
1  format(f)
elseif(type == 4) then  !real
 read(string(1:length),11,err=99) strval
11  format(e<length>.7)
endif
ng=.false.
99  return
end
