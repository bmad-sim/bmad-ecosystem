real function strval(string,type,ng)
!  format convert ascii string to real num using supplied type
implicit none
character*(*):: string
character (len=20) fmt   !change to gnu.org suggested form
integer:: length,type,ival
logical:: ng
ng=.true.      !if conversion fails
length=len_trim(string)  !length of user (sub)string
strval=0.0
if(type < 0)  return  !not a number
if(type == 2) then  !integer
  write(fmt,*) length-1   !for Hex, Octal 
 if((string(1:1) == 'H').or.(string(1:1) == 'h')) then
  read(string(2:length),"(z"//adjustl(fmt)//")",err=99) ival !hex
 elseif(string(1:1) == '"') then
  read(string(2:length),"(o"//adjustl(fmt)//")",err=99) ival !octal
 else
  write(fmt,*) length
  read(string(1:length),"(i"//adjustl(fmt)//")",err=99) ival !base ten integer
 endif
 ng=.false.
 strval=ival            !return uniformly as float (=>use nint)
elseif(type>2) then  !real
 read(string(1:length),*,err=99) strval
endif
ng=.false.
99  return
end
