program very_damped_map
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none
type(internal_state),target :: state
type(c_damap)  m,a1,m_n,R_c
type(c_vector_field) T_prime,vf
type(c_taylor) K,x,p 
type(c_normal_form) N,N_N
real(dp) prec,alpha,nu,radius,z(2),rad0,r1,r1p
complex(dp) zz(2)
integer map_order,io,km,ji(2),ns,mf,a,i,mf1
real(dp), allocatable :: co(:)
logical :: normalise=.false.

prec=1.d-13 ! for printing
longprint=.false. 
c_verbose=.false.
 
state=only_2d0

write(6,*) " map order and alpha "
read(5,*) map_order, alpha
allocate(co(0:map_order/2))
call init_all(state,map_order,0)
call alloc(m,a1,m_n,R_c)
call alloc(T_prime);call alloc(vf);call alloc(K,x,p);
call alloc(N);call alloc(N_N)
  
nu=0.4433d0
x=1.d0.cmono.1                     !(1a)
p=1.d0.cmono.2                     !(1b)
 
 k = (1.0_dp+alpha-x**2)*(cos(twopi*nu)*x + sin(twopi*nu)*p)  ! (2a)
 p = cos(twopi*nu)*p - sin(twopi*nu)*x                         ! (2b)  
 x = k                                             ! (2c)                      
 
M%v(1)=x                   ! (3a)
M%v(2)=p                   ! (3b)

call  c_normal(m,n)        ! (4)

if(normalise) then         ! (5a)
  remove_tune_shift=.true.
  m_n=n%a_t**(-1)*m*n%a_t
  call  c_normal(m_n,n_n)
  m_n=to_phasor()*m_n*from_phasor()
  a1= to_phasor()*n_n%a_t**(-1)*from_phasor()
  call  flatten_c_factored_lie(n_n%ker,vf)
  T_prime=a1*vf
else                       ! (5b)
  m_n=to_phasor()*n%a_t**(-1)*m*n%a_t*from_phasor()
  vf%v(1)=+i_*(1.d0.cmono.1)*(twopi*n%tune(1))
  vf%v(2)=-i_*(1.d0.cmono.2)*(twopi*n%tune(1))
  R_c=exp(vf,m_n)
  T_prime=log(R_c)
endif 

radius=sqrt(4*alpha/(2*cos(twopi*nu)**2+1))  ! (6)
rad0=radius; write(6,*) " naive average radius = ",radius;

a=int(map_order/2)
do i=0,a; ji(1)=i+1;ji(2)=i;co(i)=T_prime%v(1).sub.ji;enddo; ! (7)
  
 radius=rad0**2    ! Newton search for equilibrium radius
do i=1,10
  r1=0;r1p=0;
 do io=0,a; r1 =co(io)*radius**io+r1; enddo;
 do io=1,a; r1p =io*co(io)*radius**(io-1)+r1p; enddo;              ! (8)
 radius=radius-r1/r1p; write(6,*) sqrt(radius); 
enddo
radius=sqrt(radius)
write(6,*) map_order,"th order radius for alpha = ",alpha,"tune =",nu 
write(6,*) "radius =" ,radius

call kanalnummer(mf,"plot.dat")
call kanalnummer(mf1,"naive.dat")
ns=1000
do i=1,ns
m%v(1)=radius*cos(twopi*i/ns)   ! results of perturbation theory            ! (9a)
m%v(2)=radius*sin(twopi*i/ns)
m=n%a_t.o.m !A_t(z) where z=(radius*cos(twopi*i/ns),radius*sin(twopi*i/ns)) ! (9b)
z(1)=m%v(1).sub.'0'
z(2)=m%v(2).sub.'0'
write(mf,*) z
 z(1)=rad0*cos(twopi*i/ns)   ! naive result 
 z(2)=rad0*sin(twopi*i/ns)
write(mf1,*) z
enddo     
close(mf)
close(mf1)
 
end program very_damped_map