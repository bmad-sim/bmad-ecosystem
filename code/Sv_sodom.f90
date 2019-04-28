module sodom_like
use pointer_lattice
implicit none
 private equal_sodom_fourier


type c_sodom_fourier
 type(complex_quaternion), pointer :: qd(:,:)=>null()   ! 0:nphix,0:nphiy Location of quaternions
 type(complex_quaternion), pointer ::  qout(:,:)=>null() !  -mx:my, -my:my  Fourier terms
 type(complex_quaternion), pointer ::  atot(:,:)=>null() !  -mx:my, -my:my  Fourier terms
 integer nphix,nphiy 
  real(dp)  rx,ry,mux,muy,phix,phiy   !  rx=sqrt(Ix), ry = sqrt(Iy), mux, muy invariant phase advance
  real(dp) dphix,dphiy ! phase spacing  twopi/nphix
 integer  mx,my ! number of Fourier modes  -mx:mx , -my:my
 real(dp) closed_orbit(6),isf(3)
 real(dp) a(4,4),ai(4,4),delta,nu_res 
 complex(dp) epsilon
 integer no,nrx,nry,nrs
 logical nr
 type(complex_quaternion) logres
 type(internal_state) state
 type(integration_node), pointer:: t
end type c_sodom_fourier
 

   INTERFACE assignment (=)
     MODULE PROCEDURE equal_sodom_fourier
  end  INTERFACE


contains



 subroutine alloc_sodom_fourier(f,g)
 implicit none
 type(c_sodom_fourier), intent(inout):: f
 type(c_sodom_fourier), optional, intent(in):: g
 integer i,j

 if(present(g)) then
   f%nphix=g%nphix
   f%nphiy=g%nphiy
   f%closed_orbit= g%closed_orbit
  f%phix=g%phix
  f%phiy=g%phiy
  f%t=>g%t
  f%no=g%no
  f%a= g%a
  f%ai= g%ai
  f%mux=g%mux
  f%muy=g%muy
  f%rx=g%rx
  f%ry=g%ry
  f%state=g%state

 f%nr= g%nr
 f%nrx= g%nrx
 f%nry= g%nry
 f%nrs= g%nrs
  f%delta=g%delta
  f%epsilon=g%epsilon
  f%nu_res=g%nu_res
  f%logres=g%logres
  f%isf=g%isf
 allocate(f%qd(0:f%nphix,0:f%nphiy)) 
 allocate(f%qout(-f%mx:f%mx,-f%my:f%my))
 allocate(f%atot(-f%mx:f%mx,-f%my:f%my))
 f%qd=g%qd
 f%qout=g%qout
 f%atot=g%atot

 f%nrx=g%nrx
 f%nry=g%nry
 f%nrs=g%nrs
 f%nr=g%nr
else
 f%closed_orbit=0.0_dp
 f%phix=0.0_dp
 f%phiy=0.0_dp
  f%t=>null()
  f%no=1
  f%a=0.0_dp
  f%ai=0.0_dp
  f%mux=0.0_dp
  f%muy=0.0_dp
  f%rx=0
  f%ry=0
  f%state=only_4d0
 f%nr=.false.
 f%nrx=0
 f%nry=0
 f%nrs=0
  f%delta=0
  f%epsilon=0
  f%nu_res=0
  f%logres=0.0_dp
  f%isf=0
  f%isf(2)=1.0_dp
 allocate(f%qd(0:f%nphix,0:f%nphiy)) 
  allocate(f%atot(-f%mx:f%mx,-f%my:f%my))
 allocate(f%qout(-f%mx:f%mx,-f%my:f%my))
do i=0,3
 f%qd%x(i)=0.0_dp
 f%qout%x(i)=0.0_dp
 f%atot%x(i)=0.0_dp
enddo

 f%qout(0,0)%x(0)=1.0_dp
 f%atot(0,0)%x(0)=1.0_dp
 f%qd%x(0)=1.0_dp
f%dphix=twopi/f%nphix
f%dphiy=twopi/f%nphiy
endif

 



 end subroutine alloc_sodom_fourier

 subroutine equal_sodom_fourier(f,g)
 implicit none
 type(c_sodom_fourier), intent(inout) :: f
 type(c_sodom_fourier), intent(in) :: g
 integer i,j
 type(c_sodom_fourier) t

 
 f%closed_orbit= g%closed_orbit
 f%phix= g%phix
 f%phiy= g%phiy
 f%t=> g%t

  f%no=g%no
  f%a=g%a
  f%ai=g%ai

  f%ai=g%ai
  f%ai=g%ai
  f%ai=g%ai
  f%ai=g%ai

 f%nr= g%nr
 f%nrx= g%nrx
 f%nry= g%nry
 f%nrs= g%nrs
  f%delta=g%delta
  f%epsilon=g%epsilon
  f%nu_res=g%nu_res
  f%logres=g%logres
  f%isf=g%isf

  f%qd=g%qd
  f%qout=g%qout
  f%atot=g%atot
 
 end subroutine equal_sodom_fourier

 subroutine data_c_sodom_fourier(fq,xi)
 implicit none
 type(layout), target :: r
 type(c_sodom_fourier) fq
 real(dp)  xi(6)
 type(c_damap) c_map,id_s
 type(c_normal_form) c_n
 type(probe_8) xs
 type(probe) xs0
 type(internal_state) state
 
 integer i,j,n(2),count
 real(dp) x(6),phi(2) ,xr(6)
type(c_ray) cray

 
CALL FIND_ORBIT_x(fq%closed_orbit,fq%STATE,c_1d_5,node1=fq%t)  ! (3)
 
XS0=fq%closed_orbit
 
write(6,'(6(1x,g20.13))')fq%closed_orbit

state=fq%state+spin0

fq%phix=0
fq%phiy=0
call init_all(STATE,fq%no,0)

call alloc(c_map,id_s)
call alloc(c_n)
call alloc(xs)

id_s=1

xs=xs0 + id_s

 CALL propagate(XS,STATE,node1=fq%t)  ! (4)

c_map=xs
call  c_normal(c_map,c_n,dospin=my_true) 
write(6,*) "c_n%tune(1:2),c_n%spin_tune"
write(6,*) c_n%tune(1:2),c_n%spin_tune
fq%mux=c_n%tune(1)*twopi 
fq%muy=c_n%tune(2)*twopi 


fq%a=c_n%atot
fq%ai=c_n%atot**(-1)
 


!!! etienne cornell !!!!

xr=xi-fq%closed_orbit



if(fq%no==1) then
!!! etienne cornell !!!!
xr(1:4)=matmul(fq%ai,xr(1:4))
 fq%rx=sqrt(xr(1)**2+xr(2)**2)
 fq%ry=sqrt(xr(3)**2+xr(4)**2)
 
 fq%phix=atan2(-xr(2),xr(1))
 fq%phiy=atan2(-xr(4),xr(3))
 
write(6,*) "initial phase ",fq%phix,fq%phiy
 
!!! etienne cornell !!!!
! x(1:4)=matmul(fq%a,x(1:4)) 
else
 cray%x=0
 cray%x(1:6)=x(1:6)

 cray=c_n%a_t.o.cray

  x = cray%x(1:6)

endif

write(6,*) " closed "
write(6,'(6(1x,g20.13))')fq%closed_orbit

!write(6,*) " start 0"



do i=0,fq%nphix
do j=0,fq%nphiy
x=0

x(1)= fq%rx*cos(i*fq%dphix+fq%phix)
x(2)=-fq%rx*sin(i*fq%dphix+fq%phix)
x(3)= fq%ry*cos(j*fq%dphiy+fq%phiy)
x(4)=-fq%ry*sin(j*fq%dphiy+fq%phiy)

if(fq%no==1) then
 x(1:4)=matmul(fq%a,x(1:4)) 
else
 cray%x=0
 cray%x(1:6)=x(1:6)

 cray=c_n%a_t.o.cray

  x = cray%x(1:6)

endif




x=x+fq%closed_orbit

XS0%x=x
XS0%q=1.0_dp

 CALL propagate(XS0,STATE,node1=fq%t) 
 
fq%qd(i,j)=XS0%q
 
!write(6,'(6(1x,g20.13))') xs0%x(1:6)
enddo
enddo
!write(6,*) " end 0"

call kill(c_map,id_s)
call kill(c_n)
call kill(xs)

 end subroutine data_c_sodom_fourier

 subroutine fourier_c_sodom(fq)
 implicit none
  type(c_sodom_fourier) fq
 real(dp) norm,phi1,phi2
 integer i,j,k1,k2,i1,j1,k11,k22


norm=fq%nphix*fq%nphiy

do i=-fq%mx,fq%mx
do j=-fq%my,fq%my

i1=lbound(fq%qout,1)+i+fq%mx
j1=lbound(fq%qout,2)+j+fq%my
 fq%qout(i1,j1)=0.0_dp

do k1=0,fq%nphix-1
do k2=0,fq%nphiy-1

phi1=k1*fq%dphix+fq%phix
phi2=k2*fq%dphiy+fq%phiy


k11=lbound(fq%qd,1)+k1
k22=lbound(fq%qd,2)+k2


 fq%qout(i1,j1) = (exp(-i_*(i*phi1+j*phi2))/norm)  * fq%qd(k11,k22)  + fq%qout(i1,j1)

enddo
enddo

enddo
enddo

end subroutine fourier_c_sodom

 subroutine log_one_res_c_sodom_qout(fq,phi)
 implicit none
  type(c_sodom_fourier) fq
  real(dp),optional :: phi(2) 
   type(complex_quaternion) q1,dq1,m1,m2,q,dq
   integer i
   complex(dp) ex
   real(dp) cpsi,spsi,epsilon,co

  ex=1.0_dp
   if(present(phi)) then
     ex=exp(i_*(fq%nrx*phi(1)+fq%nry*phi(2)))
    endif
  q1=fq%qout(0,0)
  m1=fq%qout(fq%nrx,fq%nry)*ex
  m2=fq%qout(-fq%nrx,-fq%nry)*conjg(ex)

  dq1=q1+m1+m2
  dq1=(1.0_dp/abs(dq1))*dq1
  q1=dq1
  dq1%x(0)=dq1%x(0)-1.0_dp
  q=dq1
  dq=dq1
  do i=2,1000
   co=i
   dq=-dq1*dq
   q=q+(1.0_dp/co)*dq
  enddo

   fq%logres=q
q=(1.0_dp/twopi)*q

   fq%delta=q%x(2)
   fq%epsilon=(q%x(1)-i_*q%x(3)) !/2.0_dp    ! book definition
epsilon=abs(fq%epsilon)
   fq%nu_res= sqrt(fq%delta**2+epsilon**2)
cpsi=(ex+conjg(ex))/2.0_dp
spsi=(ex-conjg(ex))/2.0_dp/i_
   fq%isf(1)=epsilon*cpsi
   fq%isf(2)=fq%delta
   fq%isf(3)=-epsilon*spsi
   fq%isf=   fq%isf/fq%nu_res
write(6,*) fq%isf
pause 555
 end subroutine log_one_res_c_sodom_qout


 subroutine c_linear_quaternion_fourier(m_in,m_out,as) 
!#restricted: normal
!# This routine normalises the constant part of the spin matrix. 
!# m_out=as**(-1)*m_in*as
  implicit none
  type(complex_quaternion), intent(inout) :: m_in,m_out,as
  type(complex_quaternion) q0,q1,e_y,q3,qs
  real(dp) alpha,cosalpha,sinalpha
  integer i

q0=m_in 
alpha=1.0d0/sqrt(q0%x(0)**2+q0%x(1)**2+q0%x(2)**2+q0%x(3)**2)

q0=alpha*q0   ! normalised


         as=1

q1=q0
q1%x(0)=0.0_dp
qs=0.0_dp
qs%x(0)=1.0_dp/sqrt(q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
q1=q1*qs   ! q1=n

e_y=0.0_dp
e_y%x(2)=1.0_dp
 

q3=q1*e_y

 ! q3 =-n.j + n x j . l

cosalpha=-q3%x(0)

sinalpha=sqrt(q3%x(1)**2+q3%x(2)**2+q3%x(3)**2)



alpha= atan2(sinalpha,cosalpha)



if(alpha==0.and.cosalpha/=-1.0_dp) then
! write(6,*)sinalpha,cosalpha
! write(6,*) "weird in c_normal_spin_linear_quaternion "
! pause 123 
 q3=1.0_dp


else

if(cosalpha==-1.0_dp)  then
 q3=-1.0_dp
else 
 q3%x(0)=cos(alpha/2)
 q3%x(1:3)=-sin(alpha/2)*q3%x(1:3)/sinalpha 

endif


endif

  

as=q3   

     m_out=as**(-1)*m_in*as
!        m_out=AS**(-1)*m_in*as
q0=m_out

alpha=2*atan2(real(q0%x(2)),real(q0%x(0)))
 
 end  subroutine c_linear_quaternion_fourier


 subroutine normal_fourier_c_sodom(fq)
 implicit none
  type(c_sodom_fourier) fq
type(complex_quaternion), allocatable :: at(:,:), ai(:,:)
 integer k1,k2,j,i1,j1,nr(2)
  complex(dp) q,al,qt,alp,qp,qtp
 type(complex_quaternion)m_out,a0

call c_linear_quaternion_fourier(fq%qout(0,0),m_out,a0)  

call constant_simil_fourier_c_sodom(fq,a0)



allocate(at(-fq%mx:fq%mx,-fq%my:fq%my),ai(-fq%mx:fq%mx,-fq%my:fq%my))

if(.not.fq%nr) then  ! no resonance
  do k1=-fq%mx,fq%mx
  do k2=-fq%my,fq%my
  
  if(k1==0.and.k2==0) then 
   ai(0,0)=1.0_dp
  else
   q=(fq%qout(k1,k2)%x(1)-i_*fq%qout(k1,k2)%x(3))/2.0_dp
   qt=fq%qout(0,0)%x(0)+i_*fq%qout(0,0)%x(2)-(fq%qout(0,0)%x(0)-i_*fq%qout(0,0)%x(2))*exp(i_*(k1*fq%mux+k2*fq%muy))
   al=q/qt
   qp=(fq%qout(k1,k2)%x(1)+i_*fq%qout(k1,k2)%x(3))/2.0_dp
   qtp=fq%qout(0,0)%x(0)-i_*fq%qout(0,0)%x(2)-(fq%qout(0,0)%x(0)+i_*fq%qout(0,0)%x(2))*exp(i_*(k1*fq%mux+k2*fq%muy))
   alp=qp/qtp
  
  ai(k1,k2)%x(0) = 0.0_dp
  ai(k1,k2)%x(1) = al+alp
  ai(k1,k2)%x(3) = i_*(al-alp)
  ai(k1,k2)%x(2) = fq%qout(k1,k2)%x(2)/(1.0_dp-exp(i_*(k1*fq%mux+k2*fq%muy)))/fq%qout(0,0)%x(0) 
  endif
  enddo
  enddo
else
  do k1=-fq%mx,fq%mx
  do k2=-fq%my,fq%my
  
  if(k1==0.and.k2==0) then 
   ai(0,0)=1.0_dp
  else
   nr(1)=abs(fq%nrx-k1)+abs(fq%nry-k2)+abs(fq%nrs+1)
   nr(2)=abs(fq%nrx+k1)+abs(fq%nry+k2)+abs(fq%nrs-1)
   if(nr(1)/=0.and.nr(2)/=0) then
   q=(fq%qout(k1,k2)%x(1)-i_*fq%qout(k1,k2)%x(3))/2.0_dp
   qt=fq%qout(0,0)%x(0)+i_*fq%qout(0,0)%x(2)-(fq%qout(0,0)%x(0)-i_*fq%qout(0,0)%x(2))*exp(i_*(k1*fq%mux+k2*fq%muy))
   al=q/qt
   else
    al=0
   endif
   nr(1)=abs(fq%nrx-k1)+abs(fq%nry-k2)+abs(fq%nrs-1)
   nr(2)=abs(fq%nrx+k1)+abs(fq%nry+k2)+abs(fq%nrs+1)
   if(nr(1)/=0.and.nr(2)/=0) then
   qp=(fq%qout(k1,k2)%x(1)+i_*fq%qout(k1,k2)%x(3))/2.0_dp
   qtp=fq%qout(0,0)%x(0)-i_*fq%qout(0,0)%x(2)-(fq%qout(0,0)%x(0)+i_*fq%qout(0,0)%x(2))*exp(i_*(k1*fq%mux+k2*fq%muy))
   alp=qp/qtp
   else
   alp=0
   endif  
  ai(k1,k2)%x(0) = 0.0_dp
  ai(k1,k2)%x(1) = al+alp
  ai(k1,k2)%x(3) = i_*(al-alp)
  ai(k1,k2)%x(2) = fq%qout(k1,k2)%x(2)/(1.0_dp-exp(i_*(k1*fq%mux+k2*fq%muy)))/fq%qout(0,0)%x(0) 
  endif
  enddo
  enddo
endif
 
 call unit_fourier_c_sodom(fq,ai,ai)
 call phase_advance_fourier_c_sodom(fq,ai,at)

 call inv_fourier_c_sodom(fq,ai,ai)


call mul_fourier_c_sodom(fq,at,fq%qout,at)
call mul_fourier_c_sodom(fq,at,ai,fq%qout)



!!!! sum total

call constant_mul_fourier_c_sodom_r_atot(fq,a0)
call mul_fourier_c_sodom(fq,fq%atot(:,:),ai,fq%atot(:,:))

call print(ai(0,0))
call print(ai(1,1))
call print(ai(fq%mx,fq%my))

deallocate(at,ai)


end subroutine normal_fourier_c_sodom



 subroutine constant_mul_fourier_c_sodom_r_atot(fq,a0)
 implicit none
  type(c_sodom_fourier) fq
  type(complex_quaternion) a0
 integer k1,k2,i,j
 

do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my

 fq%atot(k1,k2) =  fq%atot(k1,k2)*a0


enddo
enddo
 
end subroutine constant_mul_fourier_c_sodom_r_atot


 subroutine constant_mul_fourier_c_sodom_r_qout(fq,a0)
 implicit none
  type(c_sodom_fourier) fq
  type(complex_quaternion) a0
 integer k1,k2,i,j
 

do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my

 fq%qout(k1,k2) =  fq%qout(k1,k2)*a0


enddo
enddo
 
end subroutine constant_mul_fourier_c_sodom_r_qout

 subroutine constant_simil_fourier_c_sodom(fq,a0)
 implicit none
  type(c_sodom_fourier) fq
  type(complex_quaternion) a0
 integer k1,k2,i,j
 

do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my

 fq%qout(k1,k2) =  a0**(-1)*fq%qout(k1,k2)*a0 

enddo
enddo
 
end subroutine constant_simil_fourier_c_sodom

 subroutine unit_fourier_c_sodom(fq,f1,ft)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) f1(:,:),ft(:,:)
type(complex_quaternion) q1
type(complex_quaternion), allocatable :: qd(:,:)
 integer  ij(2)
 integer i,j,k
real(dp) norm

  allocate(qd(0:fq%nphix,0:fq%nphiy))

do i=0,fq%nphix
do j=0,fq%nphiy
ij(1)=i
ij(2)=j
 call evaluate_fourier_c_sodom(fq,f1,q1,ij=ij)
 norm=1.d0/sqrt(q1%x(0)**2+q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
  q1=norm*q1
 qd(i,j)=q1
enddo
enddo
call fourier_c_quaternion_sodom(fq,qd,ft)
 
deallocate(qd)


end subroutine unit_fourier_c_sodom


 subroutine evaluate_fourier_c_sodom(fq,f,q,phi,ij)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) q,f(:,:)
 real(dp), optional :: phi(2)
 real(dp) ph(2)
 integer, optional :: ij(2)
 integer k1,k2,k11,k22
 

 q=0.0_dp
 if(present(phi)) then
  ph=phi
 else
  ph(1)=ij(1)*fq%dphix+fq%phix
  ph(2)=ij(2)*fq%dphiy+fq%phiy
 endif
do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my
k11=k1+fq%mx+lbound(f,1)
k22=k2+fq%my+lbound(f,2)
 q  =  q  + exp(i_*(k1*ph(1)+k2*ph(2)))*f(k11,k22)
enddo
enddo
 
 end subroutine evaluate_fourier_c_sodom

 subroutine fourier_c_quaternion_sodom(fq,qd,q)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) q(:,:),qd(:,:)
 real(dp) norm,phi1,phi2
 integer i,j,k1,k2,i1,j1,k11,k22


 
norm=fq%nphix*fq%nphiy

do i=-fq%mx,fq%mx
do j=-fq%my,fq%my

i1=lbound(q,1)+i+fq%mx
j1=lbound(q,2)+j+fq%my
 q(i1,j1)=0.0_dp

do k1=0,fq%nphix-1
do k2=0,fq%nphiy-1

phi1=k1*fq%dphix+fq%phix
phi2=k2*fq%dphiy+fq%phiy


k11=lbound(qd,1)+k1
k22=lbound(qd,2)+k2


 q(i1,j1) = (exp(-i_*(i*phi1+j*phi2))/norm)  * qd(k11,k22)  + q(i1,j1)

enddo
enddo

enddo
enddo

 
end subroutine fourier_c_quaternion_sodom

 subroutine phase_advance_fourier_c_sodom(fq,f,q,phi)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) q(:,:),f(:,:)
 real(dp), optional :: phi(2)
 real(dp) ph(2)
 integer k1,k2,k11,k22
 
 if(present(phi)) then
  ph=phi
 else
  ph(1)=fq%mux
  ph(2)=fq%muy
 endif
do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my
k11=k1+fq%mx+lbound(f,1)
k22=k2+fq%my+lbound(f,2)
 q(k11,k22)  =   exp(i_*(k1*ph(1)+k2*ph(2)))*f(k11,k22)
enddo
enddo
 
end subroutine phase_advance_fourier_c_sodom

 subroutine inv_fourier_c_sodom(fq,f1,ft)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) f1(:,:),ft(:,:)
type(complex_quaternion) q1
type(complex_quaternion), allocatable :: qd(:,:)
 integer  ij(2)
 integer i,j
real(dp) norm
  allocate(qd(0:fq%nphix,0:fq%nphiy))
 do i=0,fq%nphix
do j=0,fq%nphiy
ij(1)=i
ij(2)=j
 call evaluate_fourier_c_sodom(fq,f1,q1,ij=ij)
 q1=q1**(-1)
 norm=1.d0/sqrt(q1%x(0)**2+q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
  q1=norm*q1
 qd(i,j)=q1
enddo
enddo
 call fourier_c_quaternion_sodom(fq,qd,ft)

deallocate(qd)

end subroutine inv_fourier_c_sodom

 subroutine mul_fourier_c_sodom(fq,f1,f2,ft)
 implicit none
  type(c_sodom_fourier) fq
 type(complex_quaternion) f1(:,:),f2(:,:),ft(:,:)
type(complex_quaternion) q1,q2
type(complex_quaternion), allocatable :: qd(:,:)
 integer  ij(2)
 integer i,j
 real(dp) norm

  allocate(qd(0:fq%nphix,0:fq%nphiy))
 do i=0,fq%nphix
do j=0,fq%nphiy
ij(1)=i
ij(2)=j
 call evaluate_fourier_c_sodom(fq,f1,q1,ij=ij)
 call evaluate_fourier_c_sodom(fq,f2,q2,ij=ij)
 q1=q1*q2
 norm=1.d0/sqrt(q1%x(0)**2+q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
  q1=norm*q1
 qd(i,j)=q1
enddo
enddo
 call fourier_c_quaternion_sodom(fq,qd,ft)

deallocate(qd)

end subroutine mul_fourier_c_sodom


 subroutine print_fourier_sodom(fq,q,mfile,eps)
 implicit none
  type(c_sodom_fourier) fq
    integer, optional :: mfile
 type(complex_quaternion) q(:,:) 
 real(dp), optional :: eps
logical pr
    integer k1,k2,i1,j1,mf

 

      mf=6
     if(present(mfile)) mf=mfile


do k1=-fq%mx,fq%mx
do k2=-fq%my,fq%my

i1=lbound(q,1)+k1+fq%mx
j1=lbound(q,2)+k2+fq%my
 call print(q(i1,j1),0,prec=eps,pr=pr)
if(pr) write(mf,*) k1,k2
 call print(q(i1,j1),mf,prec=eps,pr=pr)
enddo
enddo



 end subroutine print_fourier_sodom

end module sodom_like