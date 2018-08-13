program Guignard_normal_form_average_x
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
type(c_taylor) fonction,fonction_FLOQUET,phase(3)
type(c_damap)  one_turn_map,id_s,U_c,U,A,fi,b,a_cs  
type(c_normal_form) normal_form
type(c_vector_field) logN,f_lin,f_non,h_left 
integer :: pos =1
integer i,map_order,mf,mf1,ns,n_mode,km 
type(probe) ray_closed
type(probe_8) ray,ray_cs_twiss
type(fibre), pointer :: p
type(integration_node), pointer :: t
character*120 :: dsc
logical :: doit=.true., used_ds_ave 
type(c_vector_field_fourier) G,F,K,f1
real(dp) DX_AVERAGE_DCS,betax_1,theta
complex(dp) coe
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

call build_lattice_als(ALS,mis,exact=.false.,thin=.true.,onecell=.true.) 

write(6,*) " Give integration step ds "
write(6,*) " > 3 and nothing is cut; each step is a full magnet"
write(6,*) "  real fun starts around ds=0.5 "
read(5,*) ts
write(dsc,*) ts

!!!!Fitting the tune and controlling the step size  !!!! 
 call kanalnummer(mf,file="fit_tune.txt")
  write(mf,*) "select layout"                  
  write(mf,*) "  1  "
  write(mf,*) "L MAX  "
  write(mf,'(a120)')dsc  !  all ds < 0.9
  write(mf,*) "CUTTING ALGORITHM "   
  write(mf,*) "2 "                ! Drifts are cut as well
  write(mf,*) "LIMIT FOR CUTTING "
  write(mf,*) "10000 10000 "
  write(mf,*) "THIN LENS "
  write(mf,*) "1000.d0 "
!  write(mf,*) "MISALIGN EVERYTHING"
!  write(mf,*) "0 0 0 0 0.1   0  5 "
  write(mf,*) "return"
 close(mf)


 call read_ptc_command77("fit_tune.txt")  


courant_snyder_teng_edwards=.true.
time_lie_choice=.true.

p=>als%start

 call MAKE_NODE_LAYOUT(als)
 
call kanalnummer(mf,"guignard_hamiltonian.txt")
call kanalnummer(mf1,"guignard_canonical_transformation.txt")

Write(6,*) " Constant phase advance per step        ---> t "
write(6,*) " Courant Snyder phase advance per step ----> f "
read(5,*) used_ds_ave


write(6,*) " Enter number of Fourier modes"
read(5,*) n_mode

n_fourier=n_mode

map_order=4                ! (3)    
call init_all(state,map_order,0)

call alloc(one_turn_map, id_s,U_c,A,U,fi,b,a_cs) 
call alloc(normal_form); call alloc(ray); 
call alloc(f_non);call alloc(f_lin)
call alloc(logN); call alloc(h_left)
call alloc(fonction,fonction_FLOQUET)
call alloc(G);call alloc(F);call alloc(K);
call alloc(phase);call alloc(ray_cs_twiss)

DX_AVERAGE_DCS=0.d0
closed_orbit=0.d0;                                            
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=1)   
 
ray_closed=closed_orbit     
id_s=1;   
ray=id_s+ray_closed;          
                          
call propagate(als,RAY,state,fibre1=1)   ! (4)

one_turn_map=ray                         ! (5)      

call c_normal(one_turn_map,normal_form)  ! (6)
U=normal_form%A_t    
 

fonction=1.0_dp.cmono.1
  call C_AVERAGE(fonction,U,fonction_FLOQUET) ! (7)
a_cs=U.sub.1
call c_canonise(a_cs,a_cs)  

 
intp=0 
intp(1)=1  

! id_s is a rotation
id_s=U**(-1)*one_turn_map*U 
call extract_linear_from_normalised(id_s,b,a,f_lin,f_non,intp) ! (8)

call  c_canonise(normal_form%A_t ,U_c)       ! (9a)
id_s=U_c.sub.1                                     ! (9b)
h_left=(1.d0/twopi)*f_lin


if(used_ds_ave) then
  circ=twopi
else
  call GET_LENGTH(als,circ) 
endif

ns=als%t%n                           
ds_ave=twopi/ns          ! (11)
ds=ds_ave                ! 
if(used_ds_ave) then
 Write(mf1,*) " Constant phase advance per step "
 write(mf1,*) "dtheta =",ds_ave, ' radians '
 write(mf1,*) "ds <= ",ts, " metres "
else
 write(mf1,*) " Approximate Courant Snyder phase advance per step "
 write(mf1,*) "ds <= ",ts
endif

U_c=U_c.cut.2

s=0.d0 
ray=U_c+ray_closed;        ! (12)
ray_cs_twiss=ray_closed+a_cs

p=>als%start               ! (13a)
t=>p%t1                    ! (13b)

do i=1,ns

 if(mod(i,ns/10)==0) then
   write(6,*) ns-i, " steps remaining "
 endif
 
 call propagate(als,ray,state,node1=i,node2=i+1)   ! (14a)
 
 call propagate(als,ray_cs_twiss,state,node1=i,node2=i+1)   ! (14b)


  ray_closed=ray ! Saving orbit  ! (15a)
  a_cs=ray_cs_twiss
  a_cs=a_cs.sub.1
  call c_full_canonise(a_cs,a_cs,phase=phase) 

  U=ray ! copying in map  ! (15b)

if(used_ds_ave) then
 U_c=exp(-ds*h_left,U)        !  (16a)
 U_c=U_c.cut.2                !  (16b)
 U=U_c**(-1)*U                !  (16c)
else 
 ds=twopi*t%s(5)/circ 
 U_c=U.cut.2                  !  (16d)
 call  c_canonise(U_c,U_c)    !  (16e)
 U=U_c**(-1)*U                !  (16f)
endif
  logN=log(U)                 !  (17)


! Checking convergence of the logarithm
a=exp(-(logN.cut.2),(U.sub.1))  ! (18)
do km=1,c_%nd2
 if(abs(full_abs(a%v(km))-1)>1.d-5) then
    call print(a,6)
    write(6,*);write(6,*) "Log failed at element ",i, p%mag%name
    stop
 endif
enddo

s=s+ds; 
do km=-n_mode,n_mode
 G%f(km)=G%f(km)+(exp(-i_*km*s)/twopi)*logN ! (19)
enddo

  betax_1=(a_cs%v(1).sub.'10')**2+(a_cs%v(1).sub.'01')**2  ! (20a)
 if((p%mag%name(1:2)=="SF".or.p%mag%name(1:2)=="SD").and.t%cas==case0) then
  DX_AVERAGE_DCS=(betax_1)**1.5_DP*p%mag%BN(3)/4.0_DP &    ! (20b) 
    *(-SIN(PHASE(1)*TWOPI)+SIN((PHASE(1)-normal_form%TUNE(1))*TWOPI)) &
    /(1.0_DP-COS(normal_form%TUNE(1)*TWOPI)) + DX_AVERAGE_DCS
 endif

 ray=ray_closed+U_c                     ! (21)
 ray_cs_twiss=ray_closed+a_cs
 t=>t%next
 p=>t%parent_fibre
enddo
 
U=from_phasor()
call transform_vector_field_fourier_by_map(G,G,u) ! (22)

 
prec=1.d-5



call alloc(f1)
if(used_ds_ave.and.doit) then
 call normalise_vector_field_fourier(G,F,K)     ! (23a)
else
 n_extra=50
 call normalise_vector_field_fourier(G,F,K,F1)  ! (23b)
endif
 
prec=1.d-11

call c_clean_vector_field_fourier(K,K,prec)
call c_clean_vector_field_fourier(g,g,prec)

write(mf1,*) " Original force "
 call print_vector_field_fourier(g,mf1)
write(mf1,*) " Normalised force "
 call print_vector_field_fourier(k,mf1)

  write(mf,*);write(mf,*)"Results for <x> ";write(mf,*)n_fourier,"modes"; write(mf,*);

  betax_1=(a_cs%v(1).sub.'10')**2+(a_cs%v(1).sub.'01')**2
 DX_AVERAGE_DCS=DX_AVERAGE_DCS*SQRT(betax_1)            ! (26) or (20c)
 write(mf,'(a11,F20.13,a20)')'  d<x>/dCS ', dx_average_dCS, " < ---- analytical"
  write(mf,*); write(mf,*) "Result of map normal form "
  call print(fonction_FLOQUET,mf)          ! (27)

theta=0.d0
call c_evaluate_vector_field_fourier(f1,theta,f_non)

f_non=to_phasor()*f_non   !  Turns f_non in Cartesian basis
id_s=exp(f_non,a_cs)
call c_evaluate_vector_field_fourier(f,theta,f_non)

f_non=to_phasor()*f_non !  Turns f_non in Cartesian basis
A=exp(f_non,id_s)


  if(used_ds_ave) then
   write(mf,*);write(mf,*) "Constant linear phase advance"
   write(mf,*) "dtheta =",ds_ave, ' radians '
   write(mf,*) "ds <= ",ts, " metres "
  else
   write(mf,*);write(mf,*) "Courant-Snyder linear phase advance "
   write(mf,*) "ds <= ",ts, " metres "
  endif
 write(mf,*);write(mf,*) "Result of Guignard normal form "
 fonction=1.0_dp.cmono.1
  call C_AVERAGE(fonction,A,fonction_FLOQUET)   ! (28)
  call print(fonction_FLOQUET,mf)                  ! (29)
close(mf)

 write(mf1,*) " Linear Courant-Snyder transformation at point of calculation for <x> "
call print(a_cs,mf1)
if(.not.used_ds_ave.or.(.not.doit)) then
 write(mf1,*) " F1 : removing the s-dependence of the phase advance"
 call print_vector_field_fourier(f1,mf1)
endif
 write(mf1,*) " F : nonlinear transformation"
 call print_vector_field_fourier(f,mf1)
close(mf1)
call ptc_end(graphics_maybe=1,flat_file=.false.)

end program Guignard_normal_form_average_x

 