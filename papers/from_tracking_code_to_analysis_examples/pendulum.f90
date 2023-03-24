program pendulum
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none
type(internal_state),target :: state
type(c_damap)  m,Id,r,a1
type(c_vector_field) vf
type(c_taylor) K,theta,p,phase,I,a_op,k_io
type(c_normal_form) N
real(dp) prec,F,dt,c4,c6,tpf
complex(dp) v
integer map_order,io,km,ji(2),ns

prec=1.d-13 ! for printing
longprint=.false. 
 
state=only_2d0

map_order=6
call init_all(state,map_order,0)
call alloc(m,id,r,a1)
call alloc(vf);call alloc(K,theta,p,phase,I,a_op,k_io);call alloc(N);

f=1/2.d0 !  This pendulum has a period of tau=2
dt=5.d-2
theta=1.d0.cmono.1                     !(1a)
p=1.d0.cmono.2                         !(1b)

K=p**2/2+2*(2*pi*f)**2*sin(theta/2.d0)**2  !(1c)

a1%v(1)=theta/sqrt(2*pi*f)             !(2a)
a1%v(2)=p*sqrt(2*pi*f)                 !(2b)

write(6,*);write(6,*) " This is Hamiltonian in linear phasors "; write(6,*);  

K=K*a1*from_phasor()                   !(2c) 
 
       call print(K,6,prec)       
       
      tpf=K.sub.'11'
       
       do io=3,map_order
      a_op=0.0_dp
       k_io=K.sub.io
       call c_taylor_cycle(K_io,size=ns)  !(3a)

!!!   Cycling over all monomials 
    do km=1,ns
       call c_taylor_cycle(K_io,ii=km,value=v,j=ji) !(3b)
       if(ji(1)/=ji(2)) then                       !(3c)
         v=v/((ji(1)-ji(2))*(2*i_)*tpf)    !(3d)
         a_op=a_op+(v.cmono.ji)                !(3e)
       endif
    enddo
    if(io==4) then
      write(6,*); write(6,*) " Eq 3.73 of the book is printed "
      call print(a_op,6)
    endif   
       vf=cgetvectorfield( a_op )                  !(4a)
       K=exp(vf,K)                                 !(4b)

       enddo
 
write(6,*);write(6,*) " This is K_new directly normalised "; write(6,*);  

       call print(K,6,prec) 
       
       write(6,*);write(6,*) " From analytical calculations "; 
write(6,*);       
write(6,99) pi*f, "(phi+ phi-) + ",-1/16.0_dp/4," (phi+ phi-)^2 + ", &
 -1/512.0_dp/pi/f/8," (phi+ phi-)^3  "
write(6,*); 
99   format(' ',(g23.16,a15,g23.16,a17,g23.16,a16))

!!!!!! Now Map Methods !!!!!!

K=p**2/2+2*(2*pi*f)**2*sin(theta/2.d0)**2  !(5)

!!! vf is the force field of the pendulum
vf=getvectorfield( -dt*K )                  !(6)
!vf%v(1)=dt*(theta.pb.K)                    !(6a)
!vf%v(2)=dt*(p.pb.K)                        !(6b)
!vf%v(1)=-dt*(K.d.2)                        !(6c)
!vf%v(2)= dt*(K.d.1)                        !(6d)

id=1                      ! (7)
m=exp(vf,id)              ! (8)

write(6,*);  
 write(6,*);write(6,*) " Normalising the map ";
write(6,*);  
call  c_normal(m,n)       ! (9)

r=n%a_t**(-1)*m*n%a_t     ! (10)

r=from_phasor(-1)*r*from_phasor()          ! (11)
call print(r,6,prec)

       phase=-i_*log(r%v(2).k.(2)).cut.map_order  ! (12a)
       phase=phase/dt                             ! (12b)
       id%v(1)=id%v(1)*2.d0                       ! (12c)
       id%v(2)=1.d0                               ! (12d)
       phase=phase.o.id                           ! (12e)
       
write(6,*);write(6,*) " Analytical results from dK_infinity/dJ Eq.(3.22) ";
write(6,*);  
       
       write(6,100) 2*pi*f,"+ ",-1.d0/8," J + ",-3.d0/512/pi/f," J^2   " ! (13a)
       
       call print(phase,6,prec)                   ! (13b)
 

       write(6,*); write(6,102) pi*f," theta^2 + ", 1/(4*pi*f)," p^2"
       c4=-1536*pi**3*f**3
       write(6,104) 80*pi**4*f**4/c4," theta^4 + ",-24*pi**2*f**2/c4, &
       " theta^2 p^2 + ",-3/c4 ," p^4"       
       
       c6=2949120*pi**5*f**5
       write(6,106) 1472*pi**6*f**6/c6," theta^6 + ",2640*pi**4*f**4/c6, &
       " theta^4 p^2 + ",1620*pi**2*f**2/c6, & 
       " theta^2 p^4 + ",135/c6 ," p^6"; write(6,*);    ! (14a)   
    
 I=((((1.d0.cmono.1)**2+(1.d0.cmono.2)**2)/2.d0)*n%a_t**(-1)).cut.map_order  ! (14b)   
  
       call print(I,6,prec)                             ! (14c)

       
102   format(' ',(g23.16,a11,g23.16,a4))
104   format(' ',(g23.16,a11,g23.16,a15,g23.16,a4))
106   format(' ',(g23.16,a11,g23.16,a15,g23.16,a15,g23.16,a4))
100    format(' ',(g23.16,a2,g23.16,a5,g23.16,a7))

    end program pendulum
    
