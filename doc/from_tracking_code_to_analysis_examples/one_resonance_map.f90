program one_resonance_map
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
real(dp) dmu,mdotmu,adotmu
 
type(internal_state),target :: state 
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, three_turn_map,id,co_moving_map
type(c_normal_form) normal_form
type(c_vector_field) Fh,F2
type(c_taylor) h3,h3t,h2,h3c,jm,ja,jx,jy
integer :: pos =1, nind(11)
integer i,map_order,mf
type(probe) ray_closed
type(probe_8) ray
type(fibre), pointer :: p
character*48 :: command_gino
!!!!!!!!!!!!!!!!!!!!!


c_verbose=.false.
prec=1.d-8 ! for printing 
use_info = .true.
longprint=.false.
 
state=only_2d0

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

call build_lattice_als(ALS,mis,exact=.false.) 
 
!!!!Fitting the tune to nu_x=0.334 !!!! 
 call kanalnummer(mf,file="fit_tune.txt")
  write(mf,*) "select layout"                  
  write(mf,*) "  1  "
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
  write(mf,*) " 0.334 , 0.2712345 "
  write(mf,*) "deallocate families"
  write(mf,*) " return "
 close(mf)
 p=>als%start; call move_to(ALS,p,"SF",pos);write(6,*) pos
 call add(p,3,0,1.d0);
 call read_ptc_command77("fit_tune.txt")
 
call kanalnummer(mf,"result_of_3nu_x.txt")

map_order=6
call init_all(state,map_order,0)

call alloc(one_turn_map, three_turn_map,id,co_moving_map)
call alloc(normal_form); call alloc(ray)
call alloc(Fh);call alloc(F2);call alloc(h3,h3t,h3c,h2);

closed_orbit=0.d0;                                            
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   
 
ray_closed=closed_orbit     
id=1;   
! ray= closed orbit + identity map  
ray=id+ray_closed;          
                          
call propagate(als,RAY,state,fibre1=pos)  

one_turn_map=ray                                      ! (1)
write(mf,*);write(mf,*)" Map produced by code ";write(mf,*);   
call print(one_turn_map,mf,prec)

three_turn_map= one_turn_map**3        !3-turn map    ! (2)   

call  c_normal(three_turn_map,normal_form)            ! (3)
write(mf,*);write(mf,*)"tune ",normal_form%tune(1);write(mf,*); 

Fh=0
do i=1,normal_form%ker%n     ! Rotation so all exponents commute
 Fh=normal_form%ker%f(i)+Fh                           ! (4)
enddo

h3=(cgetpb(Fh)*to_phasor())*normal_form%a_t**(-1)     ! (5)

write(mf,*);write(mf,*)" Invariant of the 3-turn map";write(mf,*);   
call print(h3,mf,prec)

Write(mf,*); Write(mf,*) " Checking that it is indeed invariant"
Write(mf,*) " h3*one_turn_map - h3  ";write(mf,*);
h3t=h3*one_turn_map
h3t=h3t-h3

call print(h3t,mf,prec)

write(mf,*);write(mf,*)" Invariant using the Logarithm ";write(mf,*);   

Fh=log(three_turn_map)                                ! (6)
h3=getpb(Fh)

call print(h3,mf,prec)

h3t=h3*one_turn_map
h3t=h3t-h3

Write(mf,*)
Write(mf,*) " Checking that it is indeed invariant"
Write(mf,*) " h3*one_turn_map - h3  ";write(mf,*);
call print(h3t,mf,prec)

normal_form%nres=0;normal_form%m=0; 
do i=1,map_order+1                                    ! (7a)
  if(mod(i,3)==0) then
   normal_form%nres=normal_form%nres+1               
   normal_form%m(1,normal_form%nres)=i
  endif
enddo 
 
call  c_normal(one_turn_map,normal_form)              ! (7b)

h2=(pi/3.d0)*((1.d0.cmono.1)**2+(1.d0.cmono.2)**2)    ! (8a)
F2=getvectorfield(h2)                                 ! (8b)                       

id=normal_form%a_t**(-1)*one_turn_map*normal_form%a_t ! (8c)
co_moving_map=exp(F2)*id                              ! (8d)
Fh=log(co_moving_map)                                 ! (8e)
h3c=getpb(Fh)*normal_form%a_t**(-1)                   ! (8f)

write(mf,*); write(mf,*) " Invariant of the co-moving map "; write(mf,*);   

call print(h3c,mf,prec)

h3t=h3c*one_turn_map
h3t=h3t-h3c

Write(mf,*)
Write(mf,*) " Checking that it is indeed invariant   "
Write(mf,*) " h3*one_turn_map - h3  ";write(mf,*);
call print(h3t,mf,prec)
call kill(one_turn_map, three_turn_map,id,co_moving_map)
call kill(normal_form); 
call kill(ray)
call kill(Fh);
call kill(F2);
call kill(h3,h3t,h3c,h2);

!!!!Fitting the tune to nu_x=0.37123 nu_y=0.3135  !!!!
!!!!   nu_x +2 nu_y = 0.99823
 call kanalnummer(mf,file="fit_tune.txt")
  write(mf,*) "select layout"                  
  write(mf,*) "  1  "
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
  write(mf,*) " 0.37123 , 0.3135 "
  write(mf,*) "deallocate families"
  write(mf,*) " return "
 close(mf)
 call read_ptc_command77("fit_tune.txt") 
close(mf)

call kanalnummer(mf,"result_of_nu_x+2nu_y.txt")
write(mf,*)" Results of the nu_x +2 nu_y = 1 resonance";write(mf,*);

map_order=4
state=only_4d0

call init_all(state,map_order,0)


call alloc(one_turn_map, three_turn_map,id,co_moving_map)
call alloc(normal_form); call alloc(ray)
call alloc(Fh);call alloc(F2);
call alloc(h3,h3t,h3c,h2,jm,ja,jx,jy);

closed_orbit=0.d0;                                            
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   
 
ray_closed=closed_orbit     
id=1;   
! ray= closed orbit + identity map  
ray=id+ray_closed;          
                          
call propagate(als,RAY,state,fibre1=pos)  

one_turn_map=ray                                      ! (A)

 normal_form%m=0; 
 normal_form%nres=2
 normal_form%m(1,1)=1 ; normal_form%m(2,1)=2;         ! (B1)
 normal_form%m(1,2)=2 ; normal_form%m(2,2)=4;         ! (B2)
 normal_form%positive=.false.
 
call  c_normal(one_turn_map,normal_form)              ! (C)

mdotmu=2*pi/5.d0                                             ! (D1)
adotmu=2*pi*(2*normal_form%tune(1)-normal_form%tune(2))/5.d0 ! (D2)

jx=(0.5d0.cmono.'2')+(0.5d0.cmono.'02')             ! (E1)
jy=(0.5d0.cmono.'002')+(0.5d0.cmono.'0002')         ! (E2)
jm=(jx+2*jy)                                        ! (E3)
ja=(2*jx-jy)                                        ! (E4)
 
h2=mdotmu*jm+adotmu*ja                                ! (F1)
F2=getvectorfield(h2)                                 ! (F2)                           

id=normal_form%a_t**(-1)*one_turn_map*normal_form%a_t ! (F3)
co_moving_map=exp(F2,id)                              ! (F4)

Fh=log(co_moving_map)                                 ! (F5)
h3c=getpb(Fh)*normal_form%a_t**(-1)                   ! (F6)

h3t=h3c*one_turn_map
h3t=h3t-h3c                                           ! (G)

call print(h3c,mf,prec)

Write(mf,*)
Write(mf,*) " Checking that it is indeed invariant   "
Write(mf,*) " h3*one_turn_map - h3  ";write(mf,*);
call print(h3t,mf,prec)

close(mf)


call ptc_end(graphics_maybe=1,flat_file=.false.)

end program one_resonance_map

