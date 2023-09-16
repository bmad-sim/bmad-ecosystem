program program_one_turn_map
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
real(dp) prec,closed_orbit(6),mat(6,6)
type(internal_state),target :: state
logical(lp) :: mis=.true. 
type(c_damap)  one_turn_map, Id
type(real_8) y(6)
integer i,map_order

prec=1.d-6 ! for printing
longprint=.false. 
 
call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

call build_lattice_als(ALS,mis) 

state=nocavity0

map_order=1
call init_all(state,map_order,0)
call alloc(one_turn_map,id)
call alloc(y) 
 
call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-5_dp,fibre1=1)     ! (1)

id=1   ! map is set to identity                                      ! (2)
write(6,*) " id%v(1) "
call print(id%v(1),6)
! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (3)
write(6,*) " Y(1) = closed_orbit(1)+id%v(1) "
call print(y(1),6)
call propagate(als,y(1:6),state,fibre1=1)                            ! (4)

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (5)

call print(one_turn_map,6,prec)

mat=one_turn_map                                                     ! (6)
closed_orbit=y                                                       ! (7)

write(6,*) " ";
write(6,*) "  Constant part of the map and linear part (matrix) ";
 write(6,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)
do i=1,6
 write(6,'(a5,i1,a5,6(1x,g12.5))') " row ",i," --> ",mat(i,1:6)
enddo

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program program_one_turn_map
