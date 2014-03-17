integer function csr_sparse(string,begs,ends,typ,val)
!  CSR_SPARSE  FUNCTION    SOFT        C.SOFT      SBP         95.00.00
!tpl  Ntokens=csr_sparse(string,begs,ends,typs,vals) !breakdown string
!  !!! to avoid name conflict with numerical recipies routine,
!  !! old "sparse" changed to csr_sparse !! sbp 1 sept 1998
!  Parse a string, returning start and stop of tokens,
!  as well as their type and value (floating).
!  Default delimiters are spaces tabs, and / .
!  Using subr sparse_delim, may define alternative delimiter set.
!  begs,ends = integer vectors of token start,stop char #.
!  typs = integer vector with 'type' of each substring: where
!  1=notanumber -1=illegal attemped number  2=integer  3=float 4=exp
!   5 = hyperlink
!  vals = floating point value of each token, if number. else 0.0
use csr_sparse_mod
implicit none
integer iadd,i,numtest,begs(80),ends(80),typ(80),init
real val(80),strval
integer length,iorg,ibeg,iend !temp char pointers
character string*(*),ic*1,hypstar*1,hypend*1
data hypstar,hypend /'<','>'/
logical ng,hyper
integer str_find_first_in_set,str_find_first_not_in_set
csr_sparse=0      !init # of tokens found
length=len(string)              !# chars in user string
ibeg=1        !allow first pass
iend=0
iadd=0

do while ((ibeg /= 0).and.(iend < length)) !just go till non left
 iorg=iend+1+iadd        !start after end of prev
 if(iorg > length) goto 666  !ended on "
 ibeg=str_find_first_not_in_set(string(iorg:length),whitesp)
 if(ibeg /= 0) then     !found a start
  csr_sparse=csr_sparse+1    !incr token count
  ibeg=ibeg+iorg-1              !correct to loc in full string
  ic=string(ibeg:ibeg)     !1st char of new token
  typ(csr_sparse)=0    !clear old type
  if((strd /= char(0)).and.(ic == strd)) then  !test for explicit charstring
   ibeg=ibeg+1    !skip charstring delim
   iend=str_find_first_in_set(string(ibeg:length),strd)-1
   typ(csr_sparse)=1
   val(csr_sparse)=0.0
   iadd=1    !1 extra char to skip
  elseif(ic == hypstar) then
   typ(csr_sparse)=5    !hypertext
   iend=str_find_first_in_set(string(ibeg:length),hypend)
   iadd=0
  else
   iend=str_find_first_in_set(string(ibeg:length),whitesp)-1
   typ(csr_sparse)=0
   iadd=0
  endif
  if(iend <= 0) iend=length     !no whts => all is token
  if(iend /= length) iend=iend+ibeg-1 !correct loc as before
  begs(csr_sparse)=ibeg             !returned locators
  ends(csr_sparse)=iend
!    if(ic == strd) print *,ibeg,iend
  if(typ(csr_sparse) == 0) typ(csr_sparse)=numtest(string(ibeg:iend))
  val(csr_sparse)=strval(string(ibeg:iend),typ(csr_sparse),ng)
  if(typ(csr_sparse) /= 5)   then
   if(ng) typ(csr_sparse)=1      !not a num
  endif
 endif
enddo
666  return
end
