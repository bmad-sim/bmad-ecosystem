         integer function mode_find(vec,n)
!        find most frequent value
         implicit none
         integer:: i,j,n,ival,used,maxb,maxin,vec(n),vstack(2,1500)   !allow days worth of distict vals:: 1440

         maxb=0 ; maxin=0 ; used=0
         do i=1,n
          ival=vec(i)
          do j=1,used
           if(vstack(1,j).eq.ival) then
            vstack(2,j)=vstack(2,j)+1
            if(vstack(2,j).gt.maxin) then
             maxb=j ; maxin=vstack(2,j)
            endif
            goto 1
           endif
          enddo
          !If arrive here, make new entry
          used=used+1 ; if(used.eq.1) maxb=1 ; vstack(1,used)=ival ; vstack(2,used)=1
1         continue
         enddo
         mode_find=vstack(1,maxb)
         end function mode_find
