program standard_map
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none
type(internal_state),target :: state
type(c_damap)  m,Id,r,a1,m_n,c,rho1,n_np1
type(c_vector_field) vf,g_io
type(c_factored_lie) f_op,t_op
type(c_taylor) K,theta,p,phase,I,f_pb,k_io,t_pb,a_op
type(c_normal_form) N
real(dp) prec,F,dt,c4,c6,tpf,beta
complex(dp) v,expmu,om
integer map_order,io,km,ji(2),ns,mf,a
logical t_o

prec=1.d-13 ! for printing
longprint=.false. 
c_verbose=.false.
call kanalnummer(mf,'result_of_standard_map.txt')
 
state=only_2d0

map_order=8
call init_all(state,map_order,0)
call alloc(m,id,r,a1,m_n,c,rho1,n_np1)
call alloc(vf);call alloc(K,theta,p,phase,I,f_pb,k_io,t_pb,a_op);call alloc(N);
call alloc(f_op);call alloc(t_op);call alloc(g_io);



f=1/2.d0 !  This pendulum has a period of tau=2
dt=5.d-2
theta=1.d0.cmono.1                     !(1a)
p=1.d0.cmono.2                         !(1b)
 
! One quadratic step of integration
!  K=p**2/2+2*(2*pi*f)**2*sin(theta/2.d0)**2  

 theta = theta + (dt/2) * p                                   ! (2a)
 p = p - dt * 2*(2*pi*f)**2* sin(theta/2.d0)*cos(theta/2.d0)  ! (2b)
 theta = theta + (dt/2) * p                                   ! (2c)
 
M%v(1)=theta    ! (3a)
M%v(2)=p        ! (3b)

 write(mf,*);
 write(mf,*) " Normalising the standard map using FPP software "; 
  
call  c_normal(m,n)       ! (4)
write(mf,*) " The tune of the map is ", n%tune(1)
write(mf,*); write(mf,*) " Normalised map  ";write(mf,*); 

        id=exp(n%ker)         ! (5)
        call print(id,mf,prec)

!!!!!!!!!!!!!! Normalised Pseudo-Hamiltonian !!!!!!!!!!!!!!!        
        vf=0
        
        do io=1,size(n%ker%f)
         vf%v(1)=vf%v(1)+n%ker%f(io)%v(1)  ! (6a)
         vf%v(2)=vf%v(2)+n%ker%f(io)%v(2)  ! (6b)
        enddo    
        
        k_io=-cgetpb(vf)/dt                ! (6c)

 write(mf,*); write(mf,*) " Normalised Pseudo-Hamiltonian  ";write(mf,*); 
   
        call print(k_io,mf,prec)
        
!!!!!!!!!!!!!!!!! Map and Vector Brackets !!!!!!!!!!!!!!!!!!
write(mf,*);
write(mf,*) " Normalising the map using our software using Vector Brackets ";
write(mf,*);
  
        t_op=0 ;f_op=0;

        beta=(1-(pi*f*dt)**2)**(1.0_dp/2)/(2*pi*f)      ! (ia)
        
        a1%v(1)=sqrt(beta).cmono.1                      ! (ib)
        a1%v(2)=(1.0_dp/sqrt(beta)).cmono.2             ! (ic)
        
        c=from_phasor()                                 ! (iia)
        m_n=c**(-1)*a1**(-1)*m*a1*c                     ! (iib)
       
        rho1=m_n.sub.1                   ! (iiia)
        expmu=m_n%v(1).sub.'1'           ! (iiib)

        t_op%dir=1   ! (iva)
        f_op%dir=-1  ! (ivb)
        
 do io=2,map_order
          
         n_np1= m_n*rho1**(-1); n_np1= exp_inv(T_op,n_np1)           ! (v)
         g_io%v(1)=n_np1%v(1).sub.io      ! (via)
         g_io%v(2)=n_np1%v(2).sub.io      ! (vib)
do a=1,2
          call c_taylor_cycle(g_io%v(a),size=ns)  
  
!!!   Cycling over all monomials 
    do km=1,ns
       call c_taylor_cycle(g_io%v(a),ii=km,value=v,j=ji)   ! (vii)
       if(ji(1)-ji(2)+(-1)**a/=0) then                   ! (viiia)
         v=v/(1-expmu**(ji(2)-ji(1)-(-1)**a))            ! (viiib)
         f_op%f(io)%v(a)=f_op%f(io)%v(a)+(v.cmono.ji)   ! (viiic)
       else
         t_op%f(io)%v(a)=t_op%f(io)%v(a)+(v.cmono.ji)   ! (viiid)
       endif
    enddo
 enddo   
        m_n=exp(-f_op%f(io))*m_n ;m_n=exp(f_op%f(io),m_n) ! (ix)
 enddo
 
        t_op%f(1)%v(1)= log(expmu).cmono.1       ! (xa)
        t_op%f(1)%v(2)=-log(expmu).cmono.2       ! (xb)       
        
write(mf,*);write(mf,*) " Normalised map  "; write(mf,*);
   
        id=exp(t_op)  ! (xi)
        call print(id,mf,prec)

!!!!!!!!!!! Computing the normalised Pseudo-Hamiltonian !!!!!!!!!! 

        vf=0
        do io=1,size(t_op%f)
         vf%v(1)=vf%v(1)+t_op%f(io)%v(1)     ! (xiia)
         vf%v(2)=vf%v(2)+t_op%f(io)%v(2)     ! (xiib)
        enddo    
     k_io=-cgetpb(vf)/dt                     ! (xiii)

write(mf,*);  
write(mf,*) " Normalised Pseudo-Hamiltonian  ";
write(mf,*);  
        call print(k_io,mf,prec)
 
!!!!!!!!!!!!!! Map and Poisson Brackets !!!!!!!!!!!!!!  
write(mf,*);
write(mf,*) " Normalising the map using our software using Poisson Brackets ";
write(mf,*);
        t_op=0
        f_op=0
        beta=(1-(pi*f*dt)**2)**(1.0_dp/2)/(2*pi*f)      ! (Ia)
        
        a1%v(1)=sqrt(beta).cmono.1                      ! (Ib)
        a1%v(2)=(1.0_dp/sqrt(beta)).cmono.2             ! (Ic)
        
        c=from_phasor()                                 ! (IIa)
        m_n=c**(-1)*a1**(-1)*m*a1*c                     ! (IIb)
        
        rho1=m_n.sub.1                   ! (IIIa)
        expmu=m_n%v(1).sub.'1'           ! (IIIb)
        
        t_op%dir=1  ! same as above
        f_op%dir=-1   
        
 do io=2,map_order
         t_op%f(io)=0                    ! (IVa)
         f_op%f(io)=0                    ! (IVb)
         n_np1=m_n*rho1**(-1); n_np1=exp_inv(t_op,n_np1)           ! (V)
         
         g_io%v(1)=n_np1%v(1).sub.io      ! (VIa)
         g_io%v(2)=n_np1%v(2).sub.io      ! (VIb)
         
         k_io=cgetpb( g_io )              ! (VIc)
         
         f_pb=0.0_dp                    ! (VIIa)
         t_pb=0.0_dp                    ! (VIIb)

          call c_taylor_cycle(K_io,size=ns)   ! same as above
  
!!!   Cycling over all monomials 
    do km=1,ns
       call c_taylor_cycle(K_io,ii=km,value=v,j=ji)  
       if(ji(1)/=ji(2)) then                     !(VIIIa)
         v=v/(1-expmu**(ji(2)-ji(1)))            !(VIIIb)
         f_pb=f_pb+(v.cmono.ji)                  !(VIIIc)
       else
         t_pb=t_pb+(v.cmono.ji)                  !(VIIId)
       endif
       
    enddo
        t_op%f(io)=cgetvectorfield(t_pb)         ! (IXa)
        f_op%f(io)=cgetvectorfield(f_pb)!        ! (IXb)
        m_n=exp(-f_op%f(io))*m_n; m_n=exp(f_op%f(io),m_n) ! (X)

 enddo  
        t_op%f(1)%v(1)= log(expmu).cmono.1       ! (XIa)
        t_op%f(1)%v(2)=-log(expmu).cmono.2       ! (XIb)
       
 write(mf,*);write(mf,*) " Normalised map  "; write(mf,*);

        id=exp(t_op)                           ! (XII)
        call print(id,mf,prec)

!!!!!!!!!!!! Computing the pseudo-Hamiltonian !!!!!!!!!!!!  
write(mf,*);  
  write(mf,*) " Normalising the standard map using our own little software ";
  write(mf,*) " by first computing the 'logarithm' of the map and then  ";
  write(mf,*) " normalising it using Poisson Brackets ";
write(mf,*);  

        vf=log(M)                                   ! (1a)
        k=-getpb(vf)/dt                             ! (1b)

        beta=(1-(pi*f*dt)**2)**(1.0_dp/2)/(2*pi*f)  ! (2a)
        
        a1%v(1)=sqrt(beta).cmono.1                  ! (2b)
        a1%v(2)=(1.0_dp/sqrt(beta)).cmono.2         ! (2c)

!! This is Hamiltonian in linear phasors 

K=K*a1*from_phasor()                   ! (3a) 
 
      tpf=K.sub.'11'                   ! (3b)
       
       do io=3,map_order
       a_op=0.0_dp
       k_io=K.sub.io                      ! (4a)
       call c_taylor_cycle(K_io,size=ns)  ! (4b)

!!!   Cycling over all monomials 
    do km=1,ns
       call c_taylor_cycle(K_io,ii=km,value=v,j=ji) ! (5a)
       if(ji(1)/=ji(2)) then                        ! (5b)
         v=v/((ji(1)-ji(2))*(2*i_)*tpf)             ! (5c)
         a_op=a_op+(v.cmono.ji)                     ! (5d)
       endif
    enddo
    
       vf=cgetvectorfield( a_op )                  ! (6a)
       K=exp(vf,K)                                 ! (6b)

       enddo
write(mf,*);write(mf,*) " This is K_new directly normalised "; write(mf,*);  

       call print(K,mf,prec)  
       vf=cgetvectorfield( K )      ! (7a)
       vf=-dt*vf                    ! (7b)
       id=exp(vf)
       

  
 write(mf,*);write(mf,*) " Normalised map  "; write(mf,*);
  
        call print(id,mf,prec)       

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  write(mf,*) 
  write(mf,*) " Normalising the standard map using our own little software ";
  write(mf,*) " by first computing the 'logarithm' of the map and then  ";
  write(mf,*) " normalising it using Vector Brackets ";
    

Write(6,*)"Transform vector field with formula (F.dot.a_k^-1) o a ---> type t"
Write(6,*) " Or with exp(Lie Bracket A) F ---> type f "
read(5,*) t_o
  if(t_o) then
    Write(mf,*) " Transforming vector field with formula (F.dot.a_k^-1) o a "
  else
    Write(mf,*) " Transforming vector field with exp(Lie Bracket A) F "            
  endif
        f_op=0
        t_op=0
           
        beta=(1-(pi*f*dt)**2)**(1.0_dp/2)/(2*pi*f)      ! (1A)
        
        a1%v(1)=sqrt(beta).cmono.1                      ! (1B)
        a1%v(2)=(1.0_dp/sqrt(beta)).cmono.2             ! (1C)
        
        M=from_phasor(-1)*a1**(-1)*m*a1*from_phasor()   ! (2)

        g_io=log(M)  !                                   ! (3A)
        g_io=-(1.0_dp/dt)*g_io                           ! (3B)
        om= g_io%v(1).sub.'1'                            ! (3C)
        
       do io=2,map_order
!!!   Cycling over all monomials 
do a=1,2
          m_n%v(a)=g_io%v(a).sub.io   ! (4)
          call c_taylor_cycle(m_n%v(a),size=ns)  
 
!!!   Cycling over all monomials 
    do km=1,ns
       call c_taylor_cycle(m_n%v(a),ii=km,value=v,j=ji)  
       if(ji(1)-ji(2)+(-1)**a/=0) then                     ! (5A)
         v=v/(om*(ji(1)-ji(2)+(-1)**a))                    ! (5B)
         f_op%f(io)%v(a)=f_op%f(io)%v(a)+(v.cmono.ji)      ! (5C)
       endif 
    enddo
enddo   

    if(t_o) then    ! if true, use (F.dot.a_k^-1) o a
          id=exp(f_op%f(io))    ! (6A)
          a1=id**(-1)           ! (6B)
          r=0
         
          do ns=1,2
          do a=1,2
           r%v(ns)=g_io%v(a)*(a1%v(ns).d.a)+r%v(ns)  ! (7)
          enddo
          r%v(ns)=r%v(ns)*id                         ! (8)
          enddo
         
          do ns=1,2
            g_io%v(ns)=r%v(ns)  
          enddo    
    else       !  use a Lie Bracket operator
          g_io=exp_ad(f_op%f(io),g_io)      ! (6')
    endif
 
 enddo
   
 write(mf,*);write(mf,*) " Normalised map  "; write(mf,*);
   
       id=exp(g_io)                      ! (9)
       call print(id,mf,prec)

 write(mf,*);write(mf,*) " This is K_new directly normalised "; write(mf,*);  
 
       k=cgetpb(g_io)   ! (10)
       call print(K,mf,prec)  
       
       close(mf)

!!!!!!!!!!!!!!!!!!!   Illustrating the Logarithm    !!!!!!!!!!!!!!!!!!
!
                               
          lielib_print(3)=1   ! printing iterates
write(6,*);write(6,*) "Testing the logarithm of a map ";write(6,*);
write(6,*);write(6,*) "First Case ";write(6,*);

          vf=log(M) 


!!!!  forcing linear behaviour until 10^-15 is reached
write(6,*); write(6,*) "Second Case ";write(6,*);   
       
          vf=log(M,epso=1.d-14)

write(6,*); write(6,*) "Third Case ";write(6,*);  
extra_terms_log=.true.
          vf=log(M) 

end program standard_map