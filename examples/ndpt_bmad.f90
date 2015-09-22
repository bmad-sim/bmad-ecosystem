program one_turn_orbital_map_phase_ad
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none

    interface
       subroutine build_lattice_als0(ALS,MIS,beta)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: ALS
         real(dp) beta
         logical(lp) mis
       end subroutine build_lattice_als0
    end interface

type(layout), pointer:: ALS
real(dp) prec,closed_orbit(6),mat(6,6),a(6,6),L,closed_bmad(6)
complex(dp) ac(6,6),w(6)
real(dp) beta,gamma,alpha
type(internal_state),target :: state
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, Id,a_1,a_2,b_2
type(damap) ido,fl
type(normalform) no
type(real_8) y(6),t
type(c_normal_form) normal_form
type(c_taylor) e2,r2,z1,z2,z1_new,z2_new,e2t,e1,phase(3)
integer i,map_order,pos,mf
c_mess_up_vector=.true.; b_mess=-1.0_dp;
c_verbose=.false.
prec=1.d-6 ! for printing
longprint=.false. 
call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start
mis=.false.

call build_lattice_als0(ALS,mis, 0.5d0) 

state=nocavity0 +time0
!state=delta0

!write(6,*) "Write 't' for Courant-Snyder "
!!write(6,*) "Write 'f' for Anti-Courant-Snyder "
!read(5,*) courant_snyder_teng_edwards
ndpt_bmad=0
use_bmad_units=.false.

 call kanalnummer(mf,"junk_ptc_units.txt")
 



 
map_order=2
call init_all(state,map_order,0)


call alloc(one_turn_map,id,a_1,a_2,b_2)
call alloc(phase)
call alloc(y) 
call alloc(t)
call alloc(normal_form)
call alloc(e2,r2,z1,z2,z1_new,z2_new,e2t,e1)

                                      ! (1c)
closed_orbit=0.d0
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (2)
write(mf,*) closed_orbit

 
 
id=1   ! map is set to identity                                      ! (3)

! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (4)
 
 
 
call propagate(als,y(1:6),state,fibre1=1)                            ! (5)
 

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (6)

 
                                                 ! (7)

call  c_normal(one_turn_map,normal_form)                             ! (8a)

call print(normal_form%KER,mf)

call kill(one_turn_map,id,a_1,a_2,b_2)
call kill(phase)
call kill(y) 
call kill(t)
call kill(normal_form)
call kill(e2,r2,z1,z2,z1_new,z2_new,e2t,e1)

call init(state,map_order,0)

closed_orbit=0.d0
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (2)
write(mf,*) closed_orbit

 

call alloc(ido,fl)
call alloc(no)
call alloc(y)

ido=1
y(1:6)=closed_orbit(1:6)+ido                                          ! (4)
 
 
call propagate(als,y(1:6),state,fibre1=1)      
 




ido=y

!ido=fl*ido*fl

no=ido

do i=1,2
 no%dhdj%v(i) = no%dhdj%v(i)*twopi
enddo
do i=1,c_%nd
 call print(no%dhdj%v(i),mf)
enddo

close(mf)

call kill(ido,fl)
call kill(no)
call kill(y)


!!!!!!!!!!!!!!!!!!!!!!!!!   Using BMAD UNITS in PTC   !!!!!!!!!!!!!!!!!!!!!!!!!!!
ndpt_bmad=1
use_bmad_units=.true.
 state=nocavity0 +time0
call init_all(state,map_order,0)

call kanalnummer(mf,"junk_bmad_units.txt")
 
call alloc(one_turn_map,id,a_1,a_2,b_2)
call alloc(phase)
call alloc(y) 
call alloc(t)
call alloc(normal_form)
call alloc(e2,r2,z1,z2,z1_new,z2_new,e2t,e1)

                                      ! (1c)
closed_orbit=0.d0
!call FIND_ORBIT_LAYOUT_noda_bmad(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1) 
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (2)
write(mf,*) closed_orbit

 

id=1   ! map is set to identity                                      ! (3)

! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (4)
 
 
 
call propagate(als,y(1:6),state,fibre1=1)                            ! (5)
 

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (6)

 
                                                 ! (7)

call  c_normal(one_turn_map,normal_form)                             ! (8a)

call print(normal_form%KER,mf)

call kill(one_turn_map,id,a_1,a_2,b_2)
call kill(phase)
call kill(y) 
call kill(t)
call kill(normal_form)
call kill(e2,r2,z1,z2,z1_new,z2_new,e2t,e1)

call init(state,map_order,0)

closed_orbit=0.d0
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (2)
write(mf,*) closed_orbit

 
 

call alloc(ido,fl)
call alloc(no)
call alloc(y)

ido=1
y(1:6)=closed_orbit(1:6)+ido                                          ! (4)
 
 
 
 
call propagate(als,y(1:6),state,fibre1=1)      
 




ido=y

 

no=ido

do i=1,2
 no%dhdj%v(i) = no%dhdj%v(i)*twopi
enddo
do i=1,c_%nd
 call print(no%dhdj%v(i),mf)
enddo

close(mf)
call kill(ido,fl)
call kill(no)
call kill(y)
call ptc_end

end program one_turn_orbital_map_phase_ad

subroutine  build_lattice_als0(ALS,mis,beta)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: ALS

real(dp) :: alpha,lbend, cut, ksd, ksf ,beta
type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10
type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20
type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1
type(layout) :: sfline,sdline,sup1,supb
logical(lp) mis
type(work) w
!-----------------------------------


call make_states(.true.)
exact_model = .false.
call update_states
madlength = .false.


call set_mad(beta = beta, method = 6, step = 3)

madkind2 = drift_kick_drift


L1  = drift("L1 ",  2.832695d0);L2  = drift("L2 ",  0.45698d0);
L3  = drift("L3 ",  0.08902d0);L4  = drift("L4 ",  0.2155d0);
L5  = drift("L5 ",  0.219d0);L6  = drift("L6 ",  0.107078d0);
L7  = drift("L7 ",  0.105716d0);L8  = drift("L8 ",  0.135904d0);
L9  = drift("L9 ",  0.2156993d0);L10 = drift("L10",  0.089084d0);
L11= drift("L11",  0.235416d0);L12= drift("L12",  0.1245d0);
L13= drift("L13",  0.511844d0);L14= drift("L14",  0.1788541d0);
L15= drift("L15",  0.1788483d0);L16= drift("L16",  0.511849d0);
L17= drift("L17",  0.1245d0);L18= drift("L18",  0.235405d0);
L19= drift("L19",  0.089095d0);L20= drift("L20",  0.2157007d0);
L21= drift("L21",  0.177716d0);L22= drift("L22",  0.170981d0);
L23= drift("L23",  0.218997d0);L24 = drift ("L24",  0.215503d0);
L25 = drift ("L25",  0.0890187d0);L26 = drift ("L26",  0.45698d0);
L27 = drift ("L27",  2.832696d0);L27a  = drift (" L27a",  0.8596d0);
L27b  = drift (" L27b",  0.1524d0);L27c  = drift (" L27c",  0.04445d0);
L27d  = drift (" L27d",  1.776246d0);ds  = drift (" DS  ", 0.1015d0);

QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+6.447435260914397D-03)
QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.3368D0-2.593018157427161D-02); 
QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  

!!! 1/2 mad-x value
ksf=-41.3355516397069748d0;
ksd=56.2564709584745489d0;

sf=sextupole ("sf",2.d0*0.1015d0, K2= ksf);
sd= sextupole("sd", 2.d0*0.1015d0, K2= ksd);

 VC5=marker("vc5");
ALPHA=0.17453292519943295769236907684886d0;
 
LBEND=0.86621d0;
 
BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
 
CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1.0d0,REV_FREQ=500.0d6)

sfline=1*sf;
sdline=1*sd;

SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27;

SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND1+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27;

ALS = 11*sup1+supb+cav;
 
ALS = .ring.ALS

call survey(ALS)

w=als%start

write(6,*) " beta0 = ", w%beta0
if(mis) then
 sig=1.d-5; cut=4.d0; call MESS_UP_ALIGNMENT(ALS,SIG,cut);
endif
end subroutine build_lattice_als0
