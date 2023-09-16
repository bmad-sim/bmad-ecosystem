program Guignard_Hamiltonian_cs
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none

    interface
       subroutine build_lattice_als(ALS,MIS,error,exact,sl,thin,onecell)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: ALS
         logical(lp) mis
         real(dp),optional :: error(6)
         logical, optional :: exact,sl,thin,onecell
       end subroutine build_lattice_als
    end interface

type(layout), pointer:: ALS
real(dp) prec,closed_orbit(6) 
real(dp) circ,ds,s,ds_ave,intp(2)
 
type(internal_state),target :: state 
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map,id_s,U_c,D,f,A,b,R,U
type(c_normal_form) normal_form
type(c_vector_field) Gh
type(c_taylor) h
type(c_taylor), allocatable :: hn(:)
integer :: pos =1
integer i,map_order,mf,k,ns,n_mode
type(probe) ray_closed
type(probe_8) ray
type(fibre), pointer :: p
type(integration_node), pointer :: t
character*48 :: command_gino
logical int_step,used_ds_ave,asym
integer icase
real(dp), allocatable:: theta(:),bet(:),ht(:)
real(dp) hv
!!!!!!!!!!!!!!!!!!!!!


c_verbose=.false.
prec=1.d-6 ! for printing 
use_info = .true.
longprint=.false.

Write(6,*) " Random errors -> t, no errors -> f"
read(5,*) mis 

1 write(6,*) " Choose the state "
if(.not.mis) write(6,*) "only_2d0  -> 1-d-f map if no errors -> type 1"
write(6,*) "only_4d0  -> 2-d-f map              -> type 2"
write(6,*) "delta0    -> 2-d-f map + delta      -> type 3"
write(6,*) "nocavity0 -> 3-d-f map if no cavity -> type 4"
write(6,*) "default0  -> 3-d-f map if    cavity -> type 5   <- SLOW "
read(5,*) i

    select case(i)
    case(1)
      state=only_2d0
    case(2)
      state=only_4d0  
    case(3)
      state=delta0  
    case(4)
      state=nocavity0  
    case(5)
      state=default0 
    case default 
    write(6,*) "Choose between 1 to 5"
    goto 1
    end select



call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start
call build_lattice_als(ALS,mis,exact=.false.) 

!!!!Fitting the tune and controlling the step size  !!!! 
 call kanalnummer(mf,file="fit_tune.txt")
  write(mf,*) "select layout"                  
  write(mf,*) "  1  "
  write(mf,*) "L MAX  "
  write(mf,*) "0.3 "     ! (1) all ds of the order  0.3 metre
  write(mf,*) "CUTTING ALGORITHM "   
  write(mf,*) "2 "                ! Drifts are cut as well
  write(mf,*) "LIMIT FOR CUTTING "
  write(mf,*) "10000 10000 "
  write(mf,*) "THIN LENS "
  write(mf,*) "1000.d0 "
  write(mf,*) "set families"
  write(mf,*)  "2 "
  write(mf,*) "1 NO "
  write(mf,*) "QF"
  write(mf,*) " 2, 1 " 
  write(mf,*) "1 NO "
  write(mf,*) "QD"
  write(mf,*) " 2, 2 "
  write(mf,*) "FITTUNE"
  write(mf,*) " 0.0000000001 "
  write(mf,*) " 0.3678 , 0.2712345 "
  write(mf,*) "deallocate families"
  write(mf,*) " return" 
 close(mf)

write(6,*) " Do you want to break the symmetry of the lattice "
write(6,*) " by mispowering a single sextupole ?  "
write(6,*) " Yes -> t    No -> f  "
read(5,*) asym

if(asym) then
 p=>als%start; call move_to(ALS,p,"SF",pos);  ! (2a) 
 call add(p,3,1,10d0); !call add(p,2,1,.01d0);    ! (2b)
endif 

 call read_ptc_command77("fit_tune.txt")  

 call MAKE_NODE_LAYOUT(als)
 
call kanalnummer(mf,"guignard_hamiltonian.txt")

write(6,*) " int_step t or f "
read(5,*) int_step      ! (3)
write(6,*) " Uniform ds -> t "
read(5,*) used_ds_ave   ! (4)
 
write(6,*) " case 1,2  "
write(6,*) " case = 1 -> linear transformation "
write(6,*) " case = 2 -> full nonlinear transformation: not interesting "
read(5,*) icase         ! (5)

write(6,*) " Number of Fourier modes: make it 12 or more "
read(5,*) n_mode
allocate(hn(0:n_mode))     ! (6)

map_order=4                ! (7)    

call init_all(state,map_order,0)

call alloc(one_turn_map, id_s,U_c,D,f,A,b,R,U)
call alloc(normal_form); call alloc(ray); 
call alloc(Gh);
call alloc(hn); call alloc(h);


closed_orbit=0.d0;                                            
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=1)   
 
ray_closed=closed_orbit     
id_s=1;   
ray=id_s+ray_closed;          
                          
call propagate(als,RAY,state,fibre1=1)   ! (8a)

one_turn_map=ray                         ! (8b)      

call c_normal(one_turn_map,normal_form)  ! (8c)

call  c_canonise(normal_form%A_t ,U_c,f,A,b) ! (9a)

id_s=U_c                                     ! (9b)


if(used_ds_ave) then 
 circ=twopi
else
 call GET_LENGTH(als,circ) 
endif

if(int_step) then           ! (10)
 ns=als%t%n                           
else
 ns=als%n
endif

allocate(bet(0:ns),theta(0:ns),ht(0:ns))
bet=0.d0
theta=0.d0 
ds_ave=twopi/ns             ! (11)
ds=0.d0

if(icase==1) U_c=U_c.cut.2

s=0.d0 
ray=U_c+ray_closed;        ! (12)

p=>als%start               ! (13a)
t=>p%t1                    ! (13b)

do i=1,ns

if(mod(i,100)==0) then
  write(6,*) ns-i, " steps remaining "
endif

if(used_ds_ave) then
   ds=ds_ave                ! (14a)
else
  if(int_step) then
   ds=twopi*t%s(5)/circ     ! (14b)
  else
   ds=twopi*p%mag%p%ld/circ ! (14c)
  endif 
endif


if(int_step) then
 call propagate(als,RAY,state,node1=i,node2=i+1)   ! (15a)
else
 call propagate(als,RAY,state,fibre1=i,fibre2=i+1) ! (15b)  
endif

  ray_closed=ray ! Saving orbit   
  U=ray ! copying in map  

  U_c=U

bet(i)=1.d0/((u_c%v(1).sub.'10')**2+(u_c%v(1).sub.'01')**2) ! (15c)


if(icase==1) U_c=U_c.cut.2      ! (16a)
call  c_canonise(U_c,U_c,f,A,b) ! (16b)
U=U_c**(-1)*U  ! (16c)

  Gh=log(U)      ! (17a)
  h=getpb(Gh)    ! (17b)


! Checking convergence of the logarithm
a=exp(-(Gh.cut.2),(U.sub.1))  ! (18)
do k=1,c_%nd2
 if(abs(full_abs(a%v(k))-1)>1.d-5) then
    call print(a,6)
    write(6,*);write(6,*) "Log failed at element ",i, p%mag%name
    stop
 endif
enddo

s=s+ds;  theta(i)=s
do k=0,n_mode
 hn(k)=hn(k)-exp(-i_*k*theta(i))*h/circ ! (19)
enddo 

ray=ray_closed+U_c  

if(int_step) then   ! (20)
 t=>t%next
 p=>t%parent_fibre
else
 p=>p%next
endif

enddo

if(asym) then
 write(mf,*) " Results for an asymmetric ring "
else
 write(mf,*) " Results for 12-fold symmetric ring "
endif
if(.not.used_ds_ave) then
write(mf,*) " (Circumference/pi) x Hamitonian of the ring Guignard style" 
else
 write(mf,*) " 1/pi x Hamitonian of the ring Guignard style"
endif
write(mf,*) "  "

do k=0,n_mode
write(mf,*);write(mf,*) k;write(mf,*);
hn(k)=hn(k)*from_phasor()               ! (21a)
call print((circ/pi)*hn(k),mf,prec)     ! (21b)

enddo

do i=0,ns

hv=(hn(0).sub.'11')
do k=1,n_mode
 hv=hv + 2*(hn(k).sub.'11')*exp(i_*k*theta(i))  ! (22)
enddo
 ht(i) = hv
enddo




write(mf,*) " Tracked Canonical Transformation"
write(mf,*) "  "
call print(U_c,mf,prec)
write(mf,*) " Original Canonical Transformation "
write(mf,*) "  "
call print(id_s,mf,prec)

write(mf,*) " Products with inverse   "
write(mf,*) "  "

u_c=u_c*id_s**(-1)           ! (23)
 
call print(u_c,mf,prec)

close(mf)

if(.not.used_ds_ave) then
 call kanalnummer(mf,"plot_linear_H.dat")
 
 do i=0,ns
  write(mf,*) theta(i),bet(i),ht(i)*2  ! 2J= (x + i p)*(x-i p) (24)
 enddo

close(mf)
endif

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program Guignard_Hamiltonian_cs

 