program program_ALS
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none

    interface
       subroutine build_ALS(ALS,MIS)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: ALS
         logical(lp) mis
       end subroutine build_ALS
    end interface

type(layout), pointer:: ALS
type(damapspin) id,m
real(dp) closed_orbit(6), circ,prec,beta
type(probe) ray_closed
type(probe_8) ray
type(internal_state) state
type(normal_spin) normal
integer :: no1=0,np1=0,number_of_ac_plane=0
integer  sta,i,mf,j,mfp
logical(lp) mis,tot
type(vecresonance) diagonal_field
!!!!!!  Complex FPP Objects !!!!!!
type(c_damap) cmap
type(c_normal_form) cn,cnn
type(c_taylor) x2
type(taylor) x2r
symp=.false.
prec=1.d-6 ! for printing
call ptc_ini

ALS=>m_u%start

Write(6,*) " small misalignments and rotations in lattice ? input t or f "
read(5,*) mis

call build_ALS(ALS,mis) 

1 continue
!Write(6,*) " Map with  6d without cavity -> 1 "
!Write(6,*) " Map with only 4d + delta  -> 2 "
!Write(6,*) " Map with only 6d + cavity  -> 3 "
!Write(6,*) " Map with only 6d + radiation + cavity  -> 4 "
Write(6,*) " Map with only 6d + radiation + cavity + stochasticity  -> 5 "
!read(5,*) sta
sta=5
select case(sta)
case(1)
 state=nocavity0
case(2)
 state=delta0
case(3)
 state=default0
case(4)
 state=default0+radiation0
 Write(6,*) " Do you want a total nonsymplectic normalisation?"
 Write(6,*) " Write 't' otherwise 'f' "
 read(5,*) tot
case(5)
 state=default0+radiation0+envelope0
case default
Write(6,*) "try again, wrong input "
goto 1
end select
if(sta/=5) then
 Write(6,*) " Modulate the magnet Bend1 ? yes -> t, no -> f "
 read(5,*) mis
else
 mis=my_false
endif
!!!! modulate.txt sets the magnet BEND1 as a modulated magnet !!!! 
if(mis) then
  number_of_ac_plane=1
 call kanalnummer(mf,file="modulate.txt")
  write(mf,*) "select layout"
  write(mf,*) 1
  write(mf,*) " MODULATE"
  write(mf,*) " BEND1" 
  write(mf,*) "1.d0 1.d-3 0.00d0   !   DC_ac,A_ac,theta_ac"
  write(mf,*) "1.d0   1      ! 0  !D_ac,n_ac, n_coeff "
  write(mf,*) "1 0 0.000003d0 "
  write(mf,*) "0  0 0 "
  write(mf,*) " return "
 close(mf)
 call read_ptc_command77("modulate.txt")
 state=state+modulation0
!!!! circ is the circumference of the ring !!!! 
 call get_length(als,circ)
!!!! set a modulation clock !!!!!!
 RAY_CLOSED%ac%om=0.1234d0*twopi/circ
 RAY_CLOSED%ac%x(1)=0.d0
 RAY_CLOSED%ac%x(2)=0.d0
 write(6,*) " Modulation tune =",circ*RAY_CLOSED%ac%om/twopi
endif
write(6,*) " bmad units ? "
read(5,*) i

if(i==1) then
 use_bmad_units=.true.
 ndpt_bmad=1
endif
 if(use_bmad_units) then
   call kanalnummer(mf,file="results_bmad.txt")
 else
    call kanalnummer(mf,file="results.txt")
endif
if(sta==5) then
 no1=1
  write(6,*) " Order set to 1 for stochastic radiation calculation "
else
 do while(no1<1) 
  write(6,*) " input an order > 0 (not too big first)! "
  read(5,*) no1
 enddo
endif

np1=0    !!!  No system parameters in the map
call init_all(state,no1,np1) 
!call init(state,no1,np1) ! PTC and FPP are properly initialized
!!! Complex TPSA/FPP initialization
write(6,*)  " c_%nd,np1=c_%np,ndpt1=c_%ndpt,ac_rf=number_of_ac_plane,ptc=my_true"
write(6,*)  c_%nd,c_%np,c_%ndpt,number_of_ac_plane
!call c_init(c_%NO,c_%nd,np1=c_%np,ndpt1=c_%ndpt,ac_rf=number_of_ac_plane,ptc=my_true)  


 
call alloc(cmap); call alloc(cn);

call alloc(id,m)  
call alloc(ray)
call alloc(normal)
call alloc(diagonal_field)
if(tot) normal%n%jtune=1

closed_orbit=0.d0 ! initial guess for closed orbit

call FIND_ORBIT_probe_x(ALS,closed_orbit,state,1.d-5,fibre1=1)

write(6,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)
write(mf,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)
 
RAY_CLOSED=closed_orbit

id=1;   ray=RAY_CLOSED+id;   ! ray= closed orbit + identity map

CALL TRACK_PROBE(ALS,RAY,STATE,FIBRE1=1) ! One turn map is computed via the ray
!if(use_bmad_units) call convert_ptc_to_bmad_stochastic(ray,als%start%beta0,my_false,display=my_true)
call print(ray,mf)


!!! The ray contains a truncated power series algebra
!!! The map we are interested to compute is around the closed orbit: 
!!!  it is really part of a diffential algebra
m=ray   !   The ray is "officially" turned into a "damapspin"  (DA in Berz's honour)

normal = m  ! The map is normalised Orbital+Modulation

 ! The normalised vector field is put in phasors' basis
diagonal_field=normal%n%normal%nonlinear


write(mf,*) " "
write(mf,*) " Real Tunes (damping or momentum compaction effects )"
write(mf,*) " "
select case(sta)
case(1)
write(mf,*) " 6th plane should have path length information  "
case(2)
write(mf,*) " Should be zero: no damping and no momentum compaction effects"
case(3)
write(mf,*) " Should be zero: no damping and no momentum compaction effects"
case(4)
write(mf,*) " damping decrements "
end select
write(mf,*) " "

call print(diagonal_field%cos,mf,prec)

write(mf,*) " "
write(mf,*) "Ordinary Tunes in Radians (or imaginary part of vector field)"
write(mf,*) " "
call print(diagonal_field%sin,mf,prec)


!!!!  Here is the complex FPP calculation !!!!
call kanalnummer(mfp,"phasors.txt")
 cmap=to_phasor()
 write(mfp,*) " to_phasor "
 call print(cmap,mfp)
 cmap=from_phasor()   
 write(mfp,*) "  "
 write(mfp,*) " from_phasor "
call print(cmap,mfp)
close(mfp)

cmap=m


call c_normal(cmap,cn)

write(mf,*) " "
write(mf,*) " complex map tunes "
call print(cn%ker,mf,prec)

if(sta/=4.and.c_%no>1) then
call alloc(x2)

call alloc(x2r)

x2r=1.0_dp.mono.'2'  ! Creates x^2

call AVERAGE(x2r,normal%n%a_t,x2r)

write(mf,*) "  "
write(mf,*) " Average of x^2 using standard FPP "
write(mf,*) "  "
call print(x2r,mf)

x2=1.0_dp.cmono.'2' ! Creates x^2

write(mf,*) "  "
write(mf,*) " Average of x^2 using Complex FPP "
write(mf,*) "  "
call  C_AVERAGE(x2,cn%a_t,x2) 

call print(x2,mf)

write(mf,*) "  "
write(mf,*) " Half the horizontal beta functions using the matrix A "
write(mf,*) "  "

beta=(cn%a_t%v(1).sub.'1')**2+(cn%a_t%v(1).sub.'01')**2

write(mf,*) " Half the beta = (A_11)^2+(A_12)^2 = ", beta/2

elseif(c_%no<2.and.(sta/=4.and.sta/=5)) then
 write(mf,*) "  "
 write(mf,*) " Average of x^2 cannot be computed with no<2 "
 write(mf,*) "  "


endif




if(sta==5) then 

 write(mf,*) "  " 
 write(mf,*) "Equilibrium Beam Sizes " 
 write(mf,*) "  " 
 do i=1,6
  write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11)') "    Real Package -> ", "<x_",i," x_",i,"> = ",normal%s_ij0(i,i)
  write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11)') " Complex Package -> ", "<x_",i," x_",i,"> = ",real(cn%s_ij0(i,i))
 enddo 
 do i=1,6
 do j=i,6 
    if(i/=j) write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11)') "    Real Package -> ","<x_",i," x_",j,"> = ",normal%s_ij0(i,j)
    if(i/=j) write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11)') " Complex Package -> ","<x_",i," x_",j,"> = ",real(cn%s_ij0(i,j))
 enddo
 enddo 

 write(mf,*) "  " 
 write(mf,*) "Equilibrium Beam Sizes in Phasors Basis " 
 write(mf,*) "  " 
 do i=1,6
 do j=i,6 
  write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11,a5,D18.11)') "    Real Package -> ", "<h_",i," h_",j,"> = ",real(normal%s_ijr(i,j)), &
  ' + i ',aimag(normal%s_ijr(i,j))
  write(mf,'(a20,1X,a3,I1,a3,i1,a4,D18.11,a5,D18.11)') " Complex Package -> ", "<h_",i," h_",j,"> = ",real(cn%s_ijr(i,j)),' + i ' &
  ,aimag(cn%s_ijr(i,j))
 enddo  
 enddo
 


 write(mf,*) "  " 
 write(mf,*) "Equilibrium Chao emittances " 
 write(mf,*) "  " 
 do i=1,3
 write(mf,'(a20,1x,a2,i1,a5,I1,a3,i1,a6,D18.11)')"    Real Package -> ", "e_",i,"= <h_",2*i-1," h_",2*i,">/2 = " & 
  ,normal%emittance(i)  
 write(mf,'(a20,1x,a2,i1,a5,I1,a3,i1,a6,D18.11)')" Complex Package -> ", "e_",i,"= <h_",2*i-1," h_",2*i,">/2 = ",cn%emittance(i)  
 enddo

 
endif
 
close(mf)

write(6,*) "   "
write(6,*) " hit return to terminate program "
write(6,*) "   "
pause 

end program program_ALS


!=================================================================


subroutine  build_ALS(ALS,mis)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: ALS

real(dp) :: alpha,lbend, cut, ksd, ksf 
type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10
type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20
type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
 type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1
type(layout) :: sfline,sdline,sup1,supb
logical(lp) mis
!-----------------------------------

call make_states(.true.)
exact_model = .false.
!default = default + nocavity  
call update_states
madlength = .false.


call set_mad(energy = 1.5d0, method = 6, step = 3)

madkind2 = drift_kick_drift


  L1  = drift("L1 ",  2.832695d0)
  L2  = drift("L2 ",  0.45698d0)
  L3  = drift("L3 ",  0.08902d0)
  L4  = drift("L4 ",  0.2155d0)
  L5  = drift("L5 ",  0.219d0)
  L6  = drift("L6 ",  0.107078d0)
  L7  = drift("L7 ",  0.105716d0)
  L8  = drift("L8 ",  0.135904d0)
  L9  = drift("L9 ",  0.2156993d0)
  L10 = drift("L10",  0.089084d0)
   L11= drift("L11",  0.235416d0)
   L12= drift("L12",  0.1245d0)
   L13= drift("L13",  0.511844d0)
   L14= drift("L14",  0.1788541d0)
   L15= drift("L15",  0.1788483d0)
   L16= drift("L16",  0.511849d0)
   L17= drift("L17",  0.1245d0)
   L18= drift("L18",  0.235405d0)
   L19= drift("L19",  0.089095d0)
   L20= drift("L20",  0.2157007d0)
   L21= drift("L21",  0.177716d0)
   L22= drift("L22",  0.170981d0)
   L23= drift("L23",  0.218997d0)
 L24 = drift ("L24",  0.215503d0)
 L25 = drift ("L25",  0.0890187d0)
 L26 = drift ("L26",  0.45698d0)
 L27 = drift ("L27",  2.832696d0)
 L27a  = drift (" L27a",  0.8596d0)
 L27b  = drift (" L27b",  0.1524d0)
 L27c  = drift (" L27c",  0.04445d0)
 L27d  = drift (" L27d",  1.776246d0)
 ds  = drift (" DS  ", 0.1015d0)

  QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+6.447435260914397D-03)
  QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
  QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.3368D0-2.593018157427161D-02); 
  QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
  QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
  QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  

!!! 1/2 mad-x value
ksf=(-41.67478927130080d0+0.3392376315938252d0);ksd= (56.36083889436033d0-0.1043679358857811d0);
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


if(mis) then
 sig=1.d-5
 cut=4.d0
 call MESS_UP_ALIGNMENT(ALS,SIG,cut)
endif
end subroutine build_ALS
