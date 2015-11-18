!      call index_vec(string,sub,locs,nfound)
!      Get locs in 'string'  where 'sub' occurs. If sub=' ', counts blanks.
!      if sub='xxx' will find 4 (not 2) occurances in string 'hello xxxxxx'
       subroutine index_vec(str,sub,locs,nf)
       implicit none
       character*(*) str,sub
       character*200 temp,clr*1
       integer:: locs(*),nf,sl,i,j
       logical:: got
       nf=0 ; sl=max(1,len_trim(sub))
       temp=str ; got=.true. ;clr=' ' ;if(sub(1:1).eq.' ') clr='*'
       do while(got)
        i=index(temp,sub(1:sl)) ; got=i.gt.0
        if(got) then
         nf=nf+1 ; locs(nf)=i ; temp(j:j)=clr
        endif
       enddo
       return
       end
        
