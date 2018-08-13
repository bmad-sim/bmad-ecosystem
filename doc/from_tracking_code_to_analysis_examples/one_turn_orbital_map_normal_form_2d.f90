program one_turn_orbital_map_normal_form_2d
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
real(dp) prec,closed_orbit(6),mat(6,6),a(6,6)
complex(dp) ac(6,6),w(6)
real(dp) beta,gamma,alpha
type(internal_state),target :: state
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, Id
type(real_8) y(6)
type(c_normal_form) normal_form
type(c_taylor) e1,r2,z1,z2,z1_new,z2_new,e2
integer i,map_order
c_mess_up_vector=.true.; b_mess=-1.0_dp;
c_verbose=.false.
prec=1.d-6 ! for printing
longprint=.false. 
 
 closed_orbit=0.d0
call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

call build_lattice_als(ALS,mis) 

state=only_2d0

map_order=2
call init_all(state,map_order,0)

call alloc(one_turn_map,id)
call alloc(y) 
call alloc(normal_form)
call alloc(e1,r2,z1,z2,z1_new,z2_new,e2)

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (1)

id=1   ! map is set to identity                                      ! (2)

! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (3)

call propagate(als,y(1:6),state,fibre1=1)                            ! (4)

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (5)
closed_orbit=y                                                       ! (6)

one_turn_map=one_turn_map.sub.1                                      ! (7a)
call print(one_turn_map,6,prec)                                      ! (7b)


call  c_normal(one_turn_map,normal_form)                             ! (8)

write(6,*) " "
write(6,*) " tune = ",normal_form%tune(1)                            ! (8a)

write(6,*) " ";write(6,*) "  Constant part of the map";
 write(6,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)   ! (8b)

mat=one_turn_map
write(6,*) " ";write(6,*) " One-turn Map ";
do i=1,c_%nd2
 write(6,'(a5,i1,a5,6(1x,g12.5))') " row ",i," --> ",mat(i,1:c_%nd2) ! (8c)
enddo

a=normal_form%a_t
write(6,*) " ";write(6,*) "  Canonical transformation A";
do i=1,c_%nd2
 write(6,'(a5,i1,a5,6(1x,g12.5))') " row ",i," --> ",a(i,1:c_%nd2)  ! (8d)
enddo

mat=normal_form%a_t**(-1)*one_turn_map*normal_form%a_t
write(6,*) " ";write(6,*) "  Normal Form ";
do i=1,c_%nd2
 write(6,'(a5,i1,a5,6(1x,g12.5))') " row ",i," --> ",mat(i,1:c_%nd2) ! (8e)
enddo

z1=1.e0_dp.cmono.'10'                                                ! (9a)
z2=1.e0_dp.cmono.'01'                                                ! (9b)
z1_new=z1*normal_form%a_t**(-1)                                      ! (9c)
z2_new=z2*normal_form%a_t**(-1)                                      ! (9d)

r2=z1**2+z2**2                                                       ! (10a)
e1 = r2*normal_form%a_t**(-1)                                        ! (10b)
e2 = z1_new**2+z2_new**2                                             ! (10c)

beta = a(1,1)**2+a(1,2)**2                                           ! (11a)
gamma = a(2,1)**2+a(2,2)**2                                          ! (11b)
alpha = -a(1,1)*a(2,1)-a(1,2)*a(2,2)                                 ! (11c)
write(6,*)  " "
write(6,*) "beta = ", beta
write(6,*) "gamma = ", gamma
write(6,*) " 2 x alpha = ", 2*alpha
write(6,*)  " "
write(6,*)  " Courant-Snyder Invariant: r^2 o a^(-1) "
call print(e1,6)
write(6,*)  " Courant-Snyder Invariant : z1_new**2 + z2_new**2 "
call print(e2,6)

write(6,*) " ";write(6,*) " w_1 =",(i_-alpha)/sqrt(beta)      ! (12a)
               write(6,*) " w_2 =",-sqrt(beta)                ! (12b)

ac=from_phasor(-1) * normal_form%a_t**(-1)                    ! (13a)
ac=transpose(ac)                                              ! (13b)
 
write(6,*)" ";write(6,*)"Complex Canonical transformation A";write(6,*)" ";
write(6,'(a8,17x,a1,18x,a1,1/)') " column ","1","2"
do i=1,c_%nd2
 write(6,'(a5,i1,a5,6(1x,g12.5,1x,g12.5))') " row ",i," --> ",ac(i,1:c_%nd2)
enddo
 
call ptc_end(graphics_maybe=1,flat_file=.false.)

end program one_turn_orbital_map_normal_form_2d
