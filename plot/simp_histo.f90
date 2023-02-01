         subroutine simp_histo(n,m,vals,minv,maxv,nb,xw,yw,xb,yb,ixo,iyo,ind)
         ! create 2 plot vectors showing distrib of data in vals
         ! n= data pts, m=max counts/bin, minv,maxv = range of valuesto use for histo
         ! nb=bins, xw,xy= plot region pixels ; ixo,iyo= result vecs
         ! ind=total points in output line  
         implicit none
         integer:: i,nn,ind,ic=6  
         integer:: borg,ix0,ix1,iy,n,m,nb,xw,yw,xb,yb,ixo(*),iyo(*),b,bins(nb)
         real:: v,vals(*),minv,perb,maxv

         perb=(maxv-minv)/float(nb)
         bins(1:nb)=0  !clear
         do i=1,n
          b=1+(vals(i)-minv)/perb  ; borg=b ;  b=max(1,min(nb,b))
          bins(b)=bins(b)+1
         enddo
         
         ind=0 ; ixo(1:nb*2+2)=xb ; iyo(1:nb*2+2)=yb
         do nn=1,nb       !generate plot lines
          ix0=((nn-1)*xw)/nb ; ix1=(nn*xw)/nb  ; iy=(bins(nn)*(yw-yb))/(m+1)
          ixo(ind+1)=ix0+xb ; ixo(ind+2)=ix1+xb ; iyo(ind+1)=iy+yb ; iyo(ind+2)=iy+yb ; ind=ind+2
          v=minv+perb*nn
          if(mod(nn,2).eq.1) call gi_putreals(ic,ix0+xb,yb-16,v,1,5,2)
         enddo 

         do i=5,m-5,5
          iy=yb+(i*yw)/m
          call  gg_linedots(1,xb,iy,xb+xw,iy,5,1)  !int*4 icol,ixy12(nl)
          call gi_integer(1,xb+xw+5,iy,i,1,2)
         enddo
         call gi_putstr(1,xb+xw+5,yb+yw-10,' #')
         
         ind=ind+1 ; ixo(ind)=ix1 ; iyo(ind)=yb
         return
          
         end subroutine 
