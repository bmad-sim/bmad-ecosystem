program Guignard_normal_form
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
real(dp) circ,ds,s,ds_ave,intp(3),ts
 
type(internal_state),target :: state 
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map,id_s,U_c,U,A,fi,b  
type(c_normal_form) normal_form
type(c_vector_field) logN,f_lin,f_non,h_left 
integer :: pos =1
integer i,map_order,mf,ns,n_mode,km 
type(probe) ray_closed
type(probe_8) ray
type(fibre), pointer :: p
type(integration_node), pointer :: t
character*48 :: command_gino
logical int_step,used_ds_ave,asym
integer icase
type(c_vector_field_fourier) G,F,K
!!!!!!!!!!!!!!!!!!!!!


c_verbose=.false.
prec=1.d-6 ! for printing 
use_info = .true.
longprint=.false.

mis=.false.
state=only_2d0


call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

call build_lattice_als(ALS,mis,exact=.false.,onecell=.true.) ! (0)

!!!!Fitting the tune and controlling the step size  !!!! 
 call kanalnummer(mf,file="fit_tune.txt")
  write(mf,*) "select layout"                  
  write(mf,*) "  1  "
  write(mf,*) "L MAX  "
  write(mf,*) "0.3 "     ! (1)    ! all ds < 0.3
  write(mf,*) "CUTTING ALGORITHM "   
  write(mf,*) "2 "                ! Drifts are cut as well
  write(mf,*) "LIMIT FOR CUTTING "
  write(mf,*) "10000 10000 "
  write(mf,*) "THIN LENS "
  write(mf,*) "1000.d0 "
  write(mf,*) " return" 
 close(mf)

 p=>als%start; call move_to(ALS,p,"SF",pos);  ! (2) 
 call add(p,3,1,10d0); 
 

 call read_ptc_command77("fit_tune.txt")  

 call MAKE_NODE_LAYOUT(als)
 
call kanalnummer(mf,"guignard_hamiltonian.txt")

int_step=.true.
used_ds_ave=.true.
icase=1

write(6,*) " Enter number of Fourier modes"
read(5,*) n_mode
n_fourier=n_mode

map_order=4                ! (3)    
call init_all(state,map_order,0)

call alloc(one_turn_map, id_s,U_c,A,U,fi,b) 
call alloc(normal_form); call alloc(ray); 
call alloc(f_non);call alloc(f_lin)
call alloc(logN); call alloc(h_left)

call alloc(G);call alloc(F);call alloc(K);


closed_orbit=0.d0;                                            
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=1)   
 
ray_closed=closed_orbit     
id_s=1;   
ray=id_s+ray_closed;          
                          
call propagate(als,RAY,state,fibre1=1)   ! (4)

one_turn_map=ray                         ! (5)      

call c_normal(one_turn_map,normal_form)  ! (6)
U=normal_form%A_t    ! att=c_n%A_t*c_n%As
 
intp=0 
intp(1)=1  

! id_s is a rotation
id_s=U**(-1)*one_turn_map*U ! (7)

call extract_linear_from_normalised(id_s,b,a,f_lin,f_non,intp) ! (8)

call  c_canonise(normal_form%A_t ,U_c,fi,A,b) ! (9)
 

h_left=f_lin
h_left=(1.d0/twopi)*h_left
f_lin= (1.d0/twopi)*(from_phasor()*f_lin)
f_non= (1.d0/twopi)*(from_phasor()*f_non)   ! (10)
prec=1.d-10
write(mf,*); write(mf,*) " Results of one-turn map Normalisation" 
write(mf,*) " ____________________________lin__________________________________"
call print(f_lin,mf,prec)
write(mf,*) " ____________________________non__________________________________"
call print(f_non,mf,prec)
write(mf,*) " ______________________________________________________________";
write(mf,*);

circ=twopi
ns=als%t%n                           
ds_ave=twopi/ns          ! (11)
ds=ds_ave                ! 

U_c=U_c.cut.2

s=0.d0 
ray=U_c+ray_closed;        ! (12)

p=>als%start               ! (13a)
t=>p%t1                    ! (13b)

do i=1,ns

 if(mod(i,1000)==0) then
   write(6,*) ns-i, " steps remaining "
 endif
 
 call propagate(als,RAY,state,node1=i,node2=i+1)   ! (14)


  ray_closed=ray ! Saving orbit  ! (15a)
  U=ray ! copying in map  ! (15b)
 
  U_c=exp(-ds*h_left,U)        !  (15c)
  U_c=U_c.cut.2                !  (15d)
  U=U_c**(-1)*U                !  (15e)

  logN=log(U)      ! (16)

! Checking convergence of the logarithm
a=exp(-(logN.cut.2),(U.sub.1))  ! (17)
do km=1,c_%nd2
 if(abs(full_abs(a%v(km))-1)>1.d-5) then
    call print(a,6)
    write(6,*);write(6,*) "Log failed at element ",i, p%mag%name
    stop
 endif
enddo

s=s+ds; 
do km=-n_mode,n_mode
 G%f(km)=G%f(km)+(exp(-i_*km*s)/circ)*logN ! (18)
enddo

 ray=ray_closed+U_c                     ! (19)

 t=>t%next
 p=>t%parent_fibre
enddo

U=from_phasor()
call transform_vector_field_fourier_by_map(G,G,u) ! (20)

prec=1.d-5
write(mf,*); write(mf,*) " Results of Guignard Normalisation" 
write(mf,*) " One exponent k=0 with",n_fourier, "modes"
call normalise_vector_field_fourier(G,F,K)     ! (21)
call c_clean_vector_field_fourier(K,K,prec)
call print(K%f(0),mf)

write(mf,*) " Factored k=0 with",n_fourier, "modes"
K=G
call normalise_vector_field_fourier_factored(k)    ! (22)
call c_clean_vector_field_fourier(K,K,prec)
call print(K%f(0),mf)

close(mf)

! call kanalnummer(mf,file="plot.dat")
! do i=0,n_mode
!  write(6,*) " doing mode",i
! k=g
! n_fourier=i
! call normalise_vector_field_fourier_factored(K) 
! ts=imag(K%f(0)%v(1).sub.'21')
! write(mf,*) i, ts
! enddo
!close(mf)

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program Guignard_normal_form

 