program one_turn_orbital_map_phase_ad
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
real(dp) prec,closed_orbit(6),mat(6,6),a(6,6),L
complex(dp) ac(6,6),w(6)
real(dp) beta,gamma,alpha
type(internal_state),target :: state
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, Id,a_1,a_2,b_2
type(real_8) y(6)
type(c_normal_form) normal_form
type(c_taylor) e2,r2,z1,z2,z1_new,z2_new,e2t,e1,phase(3)
integer i,map_order,pos
c_mess_up_vector=.true.; b_mess=-1.0_dp;
c_verbose=.false.
prec=1.d-6 ! for printing
longprint=.false. 
call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

call build_lattice_als(ALS,mis) 

state=only_2d0

write(6,*) "Write 't' for Courant-Snyder "
write(6,*) "Write 'f' for Anti-Courant-Snyder "
read(5,*) courant_snyder_teng_edwards

pos=2
map_order=2
call init_all(state,map_order,0)

call alloc(one_turn_map,id,a_1,a_2,b_2)
call alloc(phase)
call alloc(y) 
call alloc(normal_form)
call alloc(e2,r2,z1,z2,z1_new,z2_new,e2t,e1)
z1=1.e0_dp.cmono.'10'                                                ! (1a)
z2=1.e0_dp.cmono.'01'                                                ! (1b)
r2=z1**2+z2**2                                                       ! (1c)
closed_orbit=0.0d0
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (2)

id=1   ! map is set to identity                                      ! (3)

! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (4)

call propagate(als,y(1:6),state,fibre1=1)                            ! (5)

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (6)
closed_orbit=y                                                       ! (7)

call  c_normal(one_turn_map,normal_form)                             ! (8a)

write(6,'(1/,a50,1/)') " Canonical Transformation coming from Normal Form "
call print(normal_form%a_t,6)                                        ! (8b)

call c_canonise(normal_form%a_t,a_1,phase=phase);phase(1)=0.0_dp;    ! (9a) 

if(courant_snyder_teng_edwards) then
 write(6,'(1/,a50,1/)') " Courant-Snyder Canonical Transformation  "
  else
 write(6,'(1/,a50,1/)') " Anti-Courant-Snyder Canonical Transformation  "
endif
call print(a_1,6,prec)                                               ! (9b) 
                                                           
! map is added to closed orbit and put into the 6 polymorphs

y(1:6)=closed_orbit(1:6)+a_1                                         ! (10)

call propagate(als,y(1:6),state,fibre1=1,fibre2=pos)                 ! (11)

b_2=y(1:6)                                                           ! (12a)
call c_canonise(b_2,a_2,phase=phase)                                 ! (12b)
!!!!!!!!!!!!!!!!!!!!!!!   One turn map at pos = 2 !!!!!!!!!!!!!!!!!!!!
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=pos)   ! (13)

y(1:6)=closed_orbit(1:6)+id                                          ! (14)

call propagate(als,y(1:6),state,fibre1=pos)                          ! (15)
  
one_turn_map=y(1:6)                                                  ! (16)

call  c_normal(one_turn_map,normal_form)                             ! (17)
call  c_canonise(normal_form%a_t,a_2);                               ! (18) 

!!!!!!!!!!!!!!!!!!!   Computes invariants at 1 and 2 !!!!!!!!!!!!!!!!!

e1 = r2*a_1**(-1)                                                    ! (19a)
e2t = r2*b_2**(-1)                                                   ! (19b)
e2 = r2*a_2**(-1)                                                    ! (19c)

write(6,'(1/,a54,1/)')" Invariant at position = 1 computed from one turn map "
call print(e1,6)                                                     ! (20a)
write(6,'(1/,a54,1/)')" Invariant at position = 2  tracked from position = 1 "
call print(e2t,6)                                                    ! (20b)
write(6,'(1/,a54,1/)')" Invariant at position = 2  computed from one turn map"
call print(e2,6)                                                     ! (20c)



write(6,'(1/,a54,1/)') " Phase advance from position = 1 to position = 2      "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Analytic results !!!!!!!!!!!!!!!!!!!!!!
if(pos==2.and.courant_snyder_teng_edwards) then
 L=2.832695d0 ! Length of first element L1
 beta=e1.sub.'02'
 alpha=(e1.sub.'11')/2

 write(6,"(a48,G20.13)") " Based on theory, the phase advance should be = ", &
                                                 atan(L/(beta-L*alpha))/twopi
else
 write(6,*) " Based on theory, the phase advance should be zero "
endif
!!!!!!!!!!!!!!!!!!!!!!!!! The code's results !!!!!!!!!!!!!!!!!!!!!!!!!

call print(phase(1),6,prec)                                          ! (21)

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program one_turn_orbital_map_phase_ad
