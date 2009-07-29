!+
! Module bmad_utils_mod
!
! Module for subroutines that use bmad_struct structures but do not
! call other routines in bmad_interface.
!
! ALSO: THESE ROUTINES DO NOT HAVE ACCESS TO THE OVERLOADED
! ELE1 = ELE2 AND LAT1 = LAT2.
!-

#include "CESR_platform.inc"

module bmad_utils_mod

use bmad_struct
use make_mat6_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_controller_controls (contrl, name, err)
!
! Routine to check for problems when setting up group or overlay controllers.
!
! Modules needed:
!   use bmad
!
! Input:
!   contrl(:)   -- Control_struct: control info. 1 element for each slave.
!   name        -- Character(*): Lord name. Used for error reporting.
!
! Output:
!   err         -- Logical: Set true if there is a problem. False otherwise.
!-

subroutine check_controller_controls (contrl, name, err)

implicit none

type (control_struct) contrl(:)
integer i, j
logical err
character(*) name
character(40) :: r_name = 'check_controller_controls'

!

err = .true.

do i = 1, size(contrl)
  do j = i+1, size(contrl)
    if (contrl(i)%ix_slave == contrl(j)%ix_slave .and. &
              contrl(i)%ix_attrib == contrl(j)%ix_attrib) then
      call out_io (s_error$, r_name, 'DUPLICATE SLAVE CONTROL FOR LORD: ' // name)
      return
    endif
  enddo
enddo

err = .false.

end subroutine check_controller_controls

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord (orb, vec)
! 
! Subroutine to initialize a coord_struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   vec(6) -- real(rp): Coordinate vector. If not present then taken to be zero.
!
! Output:
!   orb -- Coord_struct: Initialized coordinate.
!-

subroutine init_coord (orb, vec)

implicit none

type (coord_struct) orb
real(rp), optional :: vec(:)

!

orb%vec = 0
if (present(vec)) orb%vec = vec

orb%spin = 0

end subroutine init_coord

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)
!
! Function to convert a character string  (eg: "drift") to an index (eg: drift$).
!
! Modules needed:
!   use bmad
!
! Input:
!   key_str        -- Character(*): Name of the key. Result is case insensitive.
!   abbrev_allowed -- Logical, optional: Abbreviations (eg: "quad") allowed?
!                       Default is False. At least 3 characters are needed 
!                       (except for rfcavity elements) if True.
!
! Output:
!   key_index -- Integer: Index of the key. Set to -1 if key_name not recognized.
!-

function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)

implicit none

character(*) key_str
character(16) name

logical, optional :: abbrev_allowed
logical abbrev

integer key_index
integer i, n_name, n_match

!

n_match = 0
key_index = -1
if (key_str == '') return

call str_upcase (name, key_str)
call string_trim (name, name, n_name)

abbrev = logic_option(.false., abbrev_allowed)

do i = 1, n_key
  if (abbrev .and. (n_name > 2 .or. name(1:2) == "RF")) then
    if (name(:n_name) == key_name(i)(1:n_name)) then
      key_index = i
      n_match = n_match + 1
    endif
  else
    if (name == key_name(i)) then
      key_index = i
      return
    endif
  endif
enddo

if (abbrev .and. n_match > 1) key_index = -1  ! Multiple matches are not valid

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine zero_ele_offsets (ele)
!
! Subroutine to zero the offsets, pitches and tilt of an element.
!
! Modules needed:
!   use bmad
!
! Output:
!   ele -- Ele_struct: Element with no (mis)orientation.
!-

subroutine zero_ele_offsets (ele)

implicit none

type (ele_struct) ele

ele%value(tilt$) = 0
ele%value(x_pitch$) = 0
ele%value(y_pitch$) = 0
ele%value(x_offset$) = 0
ele%value(y_offset$) = 0
ele%value(s_offset$) = 0

ele%value(tilt_tot$) = 0
ele%value(x_pitch_tot$) = 0
ele%value(y_pitch_tot$) = 0
ele%value(x_offset_tot$) = 0
ele%value(y_offset_tot$) = 0
ele%value(s_offset_tot$) = 0

end subroutine zero_ele_offsets

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_add_pitch (ele, mat6)
!
! Subroutine to modify a first order transfer matrix to include the affect
! of an element pitch. Note that this routine does not correct the 0th order
! part of the map. It is assumed that on input the transfer map
! does not include the affect of any pitches.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element with pitches
!     %value(x_pitch_tot$) -- Horizontal pitch
!     %value(y_pitch_tot$) -- Vertical pitch
!   mat6(6,6) -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   mat6(6,6) -- Real(rp): 1st order xfer map with pitches.
!-

subroutine mat6_add_pitch (ele, mat6)

implicit none

type (ele_struct) ele
real(rp) mat6(:,:), x_pitch, y_pitch

!

if (ele%value(x_pitch_tot$) == 0 .and. ele%value(y_pitch_tot$) == 0) return

x_pitch = ele%value(x_pitch_tot$)
y_pitch = ele%value(y_pitch_tot$)

mat6(5,6) = mat6(5,6) - mat6(5,2) * x_pitch - mat6(5,4) * y_pitch

mat6(5,1) = mat6(5,1) - x_pitch * (mat6(1,1) - 1) 
mat6(5,2) = mat6(5,2) - x_pitch *  mat6(1,2)
mat6(5,3) = mat6(5,3) - x_pitch *  mat6(1,3)
mat6(5,4) = mat6(5,4) - x_pitch *  mat6(1,4)

mat6(5,1) = mat6(5,1) - y_pitch *  mat6(3,1)
mat6(5,2) = mat6(5,2) - y_pitch *  mat6(3,2)
mat6(5,3) = mat6(5,3) - y_pitch * (mat6(3,3) - 1)
mat6(5,4) = mat6(5,4) - y_pitch *  mat6(3,4)

mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2) + &
                    mat6(5,4) * mat6(1,3) - mat6(5,3) * mat6(1,4)
mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2) + &
                    mat6(5,4) * mat6(2,3) - mat6(5,3) * mat6(2,4)
mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4) + &
                    mat6(5,2) * mat6(3,1) - mat6(5,1) * mat6(3,2)
mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4) + &
                    mat6(5,2) * mat6(4,1) - mat6(5,1) * mat6(4,2)

end subroutine mat6_add_pitch

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_total_energy_to (E_tot, particle, 
!                                         gamma, kinetic, beta, pc, brho)
!
! Routine to calculate the momentum, etc. from a particle's total energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   E_tot    -- Real(rp): Total energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   pc      -- Real(rp), optional: Particle momentum
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine convert_total_energy_to (E_tot, particle, &
                                         gamma, kinetic, beta, pc, brho)

implicit none

real(rp), intent(in) :: E_tot
real(rp), intent(out), optional :: kinetic, beta, pc, brho, gamma
real(rp) pc_new, mc2

integer, intent(in) :: particle
character(20) :: r_name = 'convert_total_energy_to'

!

mc2 = mass_of(particle)
if (E_tot < mc2) then
  call out_io (s_abort$, r_name, &
            'ERROR: TOTAL ENERGY IS LESS THAN REST MASS:\f10.0\ ', E_tot)
  call err_exit
endif

pc_new = E_tot * sqrt(1.0 - (mc2/E_tot)**2)
if (present(pc))     pc     = pc_new
if (present(beta))    beta    = pc_new / E_tot  
if (present(kinetic)) kinetic = E_tot - mc2
if (present(brho))    brho    = pc_new / c_light
if (present(gamma))   gamma   = E_tot / mc2

end subroutine convert_total_energy_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho)
!
! Routine to calculate the energy, etc. from a particle's momentum.
!
! Modules needed:
!   use bmad
!
! Input:
!   pc       -- Real(rp): Particle momentum
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   E_tot   -- Real(rp), optional: Total energy of the particle.
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho)

implicit none

real(rp), intent(in) :: pc
real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma
real(rp) E_tot_new, mc2

integer, intent(in) :: particle
character(20) :: r_name = 'convert_pc_to'

!

mc2 = mass_of(particle)
E_tot_new = sqrt(pc**2 + mc2**2)

if (present(E_tot))   E_tot   = E_tot_new
if (present(beta))    beta    = pc / E_tot_new
if (present(kinetic)) kinetic = E_tot_new - mc2
if (present(brho))    brho    = pc / c_light
if (present(gamma))   gamma   = E_tot_new / mc2

end subroutine convert_pc_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine wiggler_vec_potential (ele, here, vec_pot)
!
! Subroutine to calculate the normalized vector potential at 
! a point for a wiggler. The normalized potental a_norm is defined by:
!      p_cononical = p_mv - a_norm
! The Gauge used here is the same one as used in PTC and has A_x = 0.
! 
! Modules needed:
!   use bmad
!
! Input:
!   ele     -- Ele_struct: wiggler element.
!   here    -- Coord_struct: Coordinates for calculating the vector pot.
!
! Output:
!   vec_pot(3) -- Real(rp): Normalized vector potential
!-

subroutine wiggler_vec_potential (ele, here, vec_pot)

implicit none

type (ele_struct), target, intent(in) :: ele
type (coord_struct), intent(in) :: here
real(rp), intent(out) :: vec_pot(3)

type (wig_term_struct), pointer :: t

real(rp) c_x, s_x, c_y, s_y, c_z, s_z
real(rp) x, y, s, coef

integer i

!

if (ele%key /= wiggler$) then
  print *, 'ERROR IN WIGGLER_VEC_POTENTIAL. ELEMENT NOT A WIGGLER: ', ele%name
  call err_exit
endif

!

x = here%vec(1)
y = here%vec(3)
s = here%vec(5)

vec_pot = 0

do i = 1, size(ele%wig_term)
  t => ele%wig_term(i)

    if (t%type == hyper_y$) then
      c_x = cos(t%kx * x)
      s_x = sin(t%kx * x)
    elseif (t%type == hyper_x$ .or. t%type == hyper_xy$) then
      c_x = cosh(t%kx * x)
      s_x = sinh(t%kx * x)
    else
      print *, 'ERROR IN WIGGLER_VEC_POTENTIAL: UNKNOWN TERM TYPE!'
      call err_exit
    endif

    if (t%type == hyper_y$ .or. t%type == hyper_xy$) then
      c_y = cosh (t%ky * y)
      s_y = sinh (t%ky * y)
    else
      c_y = cos (t%ky * y)
      s_y = sin (t%ky * y)
    endif

    c_z = cos (t%kz * s + t%phi_z)
    s_z = sin (t%kz * s + t%phi_z)

    coef = ele%value(polarity$) * t%coef

    vec_pot(2) = vec_pot(2) - coef  * (t%kz / (t%kx * t%ky)) * s_x * s_y * s_z
    vec_pot(3) = vec_pot(3) - coef  * (1 / t%kx)             * s_x * c_y * c_z
  enddo


end subroutine wiggler_vec_potential


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_lat_parameters (lat_in, lat_out)
!
! Subroutine to transfer the lat parameters (such as lat%name, lat%param, etc.)
! from one lat to another. The only stuff that is not transfered are the
! arrays and lat%ele_init.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat with parameters set.
!-

subroutine transfer_lat_parameters (lat_in, lat_out)

implicit none

type (lat_struct), intent(in) :: lat_in
type (lat_struct) :: lat_out

!

lat_out%name =                 lat_in%name
lat_out%lattice =              lat_in%lattice
lat_out%input_file_name =      lat_in%input_file_name
lat_out%title =                lat_in%title
lat_out%a =                    lat_in%a
lat_out%b =                    lat_in%b
lat_out%z =                    lat_in%z
lat_out%param =                lat_in%param
lat_out%version =              lat_in%version
lat_out%n_ele_track =          lat_in%n_ele_track
lat_out%n_ele_max =            lat_in%n_ele_max
lat_out%n_control_max =        lat_in%n_control_max
lat_out%n_ic_max =             lat_in%n_ic_max
lat_out%input_taylor_order =   lat_in%input_taylor_order
lat_out%beam_start =           lat_in%beam_start

end subroutine transfer_lat_parameters

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
!
! Subroutine to transfer a Taylor map from one element to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in       -- Ele_struct: Element with the Taylor map.
!   taylor_order -- Integer: Order to truncate the Taylor map at.
!
! Output:
!   ele_out      -- Ele_struct: Element receiving the Taylor map truncated to
!                     order taylor_order.
!-

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)

implicit none

type (ele_struct) ele_in, ele_out
integer, intent(in) :: taylor_order
integer it, ix, k 

!

do it = 1, 6
  ix = 0
  do k = 1, size(ele_in%taylor(it)%term)
   if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= taylor_order) ix = ix + 1
  enddo
  if (.not. associated(ele_out%taylor(it)%term)) &
                          allocate (ele_out%taylor(it)%term(ix))
  if (size(ele_out%taylor(it)%term) /= ix) &
                          allocate (ele_out%taylor(it)%term(ix))
  ix = 0
  do k = 1, size(ele_in%taylor(it)%term)
   if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= taylor_order) then
      ix = ix + 1
      ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
    endif      
  enddo
enddo

ele_out%taylor_order = taylor_order
ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

if (ele_in%key == wiggler$) ele_out%value(z_patch$) = ele_in%value(z_patch$)

end subroutine transfer_ele_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_lat (lat, n)
!
! Subroutine to initialize a BMAD lat.
! 
! Modules needed:
!   use bmad
!
! Input:
!   n    -- Integer: Upper bound lat%ele(0:) array is initialized to.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n)

implicit none

type (lat_struct)  lat
integer n

!

call deallocate_lat_pointers (lat)
call allocate_lat_ele_array(lat, n)
call init_ele (lat%ele_init)

allocate (lat%control(100))
allocate (lat%ic(100))

lat%title = ' '
lat%name = ' '
lat%lattice = ' '
lat%input_file_name = ' '
lat%param%stable = .true.
lat%param%particle = positron$

call init_coord(lat%beam_start)

call init_mode_info (lat%a)
call init_mode_info (lat%b)
call init_mode_info (lat%z)

lat%n_ele_track = 0
lat%n_ele_max = 0
lat%n_control_max = 0
lat%n_ic_max = 0
lat%input_taylor_order = 0
lat%version = -1

call allocate_branch_array (lat%branch, 0, lat)  
lat%branch(0)%name = 'MAIN'

!----------------------------------------
contains

subroutine init_mode_info (t)
type (mode_info_struct) t
t%tune = 0
t%emit = 0
t%chrom = 0
end subroutine 

end subroutine init_lat

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+ 
! Function equivalent_taylor_attributes (ele1, ele2) result (equiv)
!
! Subroutine to see if two elements are equivalent in terms of attributes so
! that their Taylor Maps would be the same. 
! If the reference orbit about which the Taylor map is made is zero then
! two elements can be equivalent even if the names are different.
!
! This routine is used to see if a taylor map from one element may be 
! used for another and thus save some computation time. Taylor map elements
! Are considered *never* to be equivalent since their maps are never computed.
!
! Modules needed:
!   use bmad
!
! Input: 
!   ele1 -- Ele_struct: Element with a Taylor map
!   ele2 -- Ele_struct: Element that might receive the Taylor map from ele1.
!
! Output:
!   equiv -- logical: True if elements are equivalent.
!-

function equivalent_taylor_attributes (ele1, ele2) result (equiv)

implicit none

type (ele_struct) :: ele1, ele2

integer it

logical equiv
logical vmask(n_attrib_maxx)

!

equiv = .false.

if (ele1%key /= ele2%key) return
if (ele1%sub_key /= ele2%sub_key) return
if (ele1%map_with_offsets .neqv. ele2%map_with_offsets) return
if (ele1%integrator_order /= ele2%integrator_order) return
if (ele1%name /= ele2%name .and. any(ele1%taylor%ref /= 0)) return

vmask = .true.
if (ele1%key == wiggler$ .and. ele1%sub_key == map_type$) then
  vmask( (/ k1$, rho$, b_max$, z_patch$, p0c$ /) ) = .false.
endif
if (.not. ele1%map_with_offsets) then
  vmask( (/ x_offset$, y_offset$, s_offset$, tilt$, x_pitch$, &
            y_pitch$, x_offset_tot$, y_offset_tot$, s_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$/) ) = .false.
endif

if (any(ele1%value /= ele2%value .and. vmask)) return

if (associated(ele1%wig_term) .neqv. associated(ele2%wig_term)) return
if (associated(ele1%wig_term)) then
  if (size(ele1%wig_term) /= size(ele2%wig_term)) return
  do it = 1, size(ele1%wig_term)
    if (ele1%wig_term(it)%coef  /= ele2%wig_term(it)%coef)  cycle
    if (ele1%wig_term(it)%kx    /= ele2%wig_term(it)%kx)    cycle
    if (ele1%wig_term(it)%ky    /= ele2%wig_term(it)%ky)    cycle
    if (ele1%wig_term(it)%kz    /= ele2%wig_term(it)%kz)    cycle
    if (ele1%wig_term(it)%phi_z /= ele2%wig_term(it)%phi_z) cycle
  enddo
endif

if (ele1%key == taylor$) return  

equiv = .true.


end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine clear_lat_1turn_mats (lat)
!
! Subroutine to clear the 1-turn matrices in the lat structure:
!   lat%param%t1_no_RF
!   lat%param%t1_with_RF
! This will force any routine dependent upon these to do a remake.
!
! Modules needed:
!   use bmad
!
! Output:
!   lat -- lat_struct: Lat with 1-turn matrices cleared.
!-

subroutine clear_lat_1turn_mats (lat)

implicit none

type (lat_struct) lat

lat%param%t1_no_RF = 0
lat%param%t1_with_RF = 0

end subroutine clear_lat_1turn_mats


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ele (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct:
!
! Output:
!   ele2 -- Ele_struct:
!-

subroutine transfer_ele (ele1, ele2)

type (ele_struct) :: ele1
type (ele_struct) :: ele2

!

ele2 = ele1

end subroutine transfer_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_eles (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1(:) -- Ele_struct:
!
! Output:
!   ele2(:) -- Ele_struct:
!-

subroutine transfer_eles (ele1, ele2)

type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(inout) :: ele2(:)

ele2 = ele1

end subroutine transfer_eles

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branch (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1 -- Branch_struct:
!
! Output:
!   branch2 -- Branch_struct:
!-

subroutine transfer_branch (branch1, branch2)

type (branch_struct) :: branch1
type (branch_struct) :: branch2

!

branch2 = branch1

end subroutine transfer_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branches (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1(:) -- Branch_struct:
!
! Output:
!   branch2(:) -- Branch_struct:
!-

subroutine transfer_branches (branch1, branch2)

type (branch_struct) :: branch1(:)
type (branch_struct) :: branch2(:)

branch2 = branch1

end subroutine transfer_branches

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_lat (lat1, lat2)
!
! Subroutine to set lat2 = lat1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end lat2's pointers point to the same memory as lat1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   lat1 -- lat_struct:
!
! Output:
!   lat2 -- lat_struct:
!-

subroutine transfer_lat (lat1, lat2)

type (lat_struct), intent(in) :: lat1
type (lat_struct), intent(out) :: lat2

lat2 = lat1

end subroutine transfer_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord (coord, n_coord)
!
! Subroutine to reallocate an allocatable  coord_struct array to at least:
!     coord(0:n_coord)
! Note: The old coordinates are not saved except for coord(0).
! If, at input, coord(:) is not allocated then coord(0)%vec is set to zero.
! In any case, coord(n)%vec for n > 0 is set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:) -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord (coord, n_coord)

type (coord_struct), allocatable :: coord(:)
type (coord_struct) start

integer, intent(in) :: n_coord
integer i

!

if (allocated (coord)) then
  if (size(coord) < n_coord + 1) then
    start = coord(0)
    deallocate (coord)
    allocate (coord(0:n_coord))
    coord(0) = start
    do i = 1, n_coord
      call init_coord (coord(i))
    enddo
  endif
else
  allocate (coord(0:n_coord))
  do i = 0, n_coord
    call init_coord (coord(i))
  enddo
endif

end subroutine reallocate_coord

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_control(lat, n) 
!
! Function to reallocate the lat%control(:) and lat%ic(:) arrays.
! The old data in the arrays will be saved.
! 
! Modules needed:
!   use bmad
!
! Input:
!   lat  -- Lat_struct: Lattice.
!   n    -- Integer: Array size for lat%control(:) and lat%ic(:).
!
! Output:
!   lat  -- Lat_struct: Lattice.
!     %control(:) -- Control Array with size at least n.
!     %ic(:)      -- Control Array with size at least n.
!-

subroutine reallocate_control (lat, n)

implicit none

type (lat_struct) lat
type (control_struct), allocatable :: control(:)
real(rp), allocatable :: ic(:)
integer, intent(in) :: n
integer n_old

!

if (.not. allocated(lat%control)) then
  allocate (lat%control(n), lat%ic(n))
  return
endif

n_old = size(lat%control)
if (n_old >= n) return

allocate (control(n_old), ic(n_old))
control = lat%control
ic = lat%ic
deallocate (lat%control, lat%ic)
allocate (lat%control(n), lat%ic(n))
lat%control(1:n_old) = control
lat%ic(1:n_old) = ic

deallocate (control, ic)

end subroutine reallocate_control

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine deallocate_ele_pointers (ele, nullify_only)
!
! Subroutine to deallocate the pointers in an element.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- ele_struct: Element with pointers.
!   nullify_only -- Logical, optional: If present and True then
!               Just nullify. Do not deallocate.
!
! Output:
!   ele -- Ele_struct: Element with deallocated pointers.
!-

subroutine deallocate_ele_pointers (ele, nullify_only)

implicit none

type (ele_struct) ele
logical, optional, intent(in) :: nullify_only

! nullify only

if (present (nullify_only)) then
  if (nullify_only) then
    nullify (ele%wig_term)
    nullify (ele%const)
    nullify (ele%r)
    nullify (ele%descrip)
    nullify (ele%a_pole, ele%b_pole)
    nullify (ele%wake)
    nullify (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
              ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)
    nullify (ele%gen_field)
    nullify (ele%mode3)
    return
  endif
endif

! Normal deallocate

if (associated (ele%wig_term)) deallocate (ele%wig_term)
if (associated (ele%const))    deallocate (ele%const)
if (associated (ele%r))        deallocate (ele%r)
if (associated (ele%descrip))  deallocate (ele%descrip)
if (associated (ele%a_pole))   deallocate (ele%a_pole, ele%b_pole)
if (associated (ele%mode3))    deallocate (ele%mode3)

if (associated (ele%wake)) then
  if (associated (ele%wake%sr_table))       deallocate (ele%wake%sr_table)
  if (associated (ele%wake%sr_mode_long))  deallocate (ele%wake%sr_mode_long)
  if (associated (ele%wake%sr_mode_trans)) deallocate (ele%wake%sr_mode_trans)
  if (associated (ele%wake%lr))        deallocate (ele%wake%lr)
  deallocate (ele%wake)
endif

if (associated (ele%taylor(1)%term)) deallocate &
         (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
         ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)

call kill_gen_field (ele%gen_field)

end subroutine deallocate_ele_pointers

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine kill_gen_field (gen_field)
!
! Subroutine to kill a gen_field.
!
! Modules needed:
!   use bmad
!
! Input:
!   gen_field -- Genfield, pointer: gen_field to kill.
!
! Output:
!   gen_field -- Genfield, pointer: Killed gen_field.
!-

subroutine kill_gen_field (gen_field)

use tpsalie_analysis, only: kill 

implicit none

type (genfield), pointer :: gen_field

!

if (associated(gen_field)) then
  call kill (gen_field)
  deallocate (gen_field)
endif

end subroutine kill_gen_field

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine init_ele (ele, key, sub_key, ix_ele, ix_branch)
!
! Subroutine to initialize a Bmad element. Element is initialized to be free
! (not a lord or slave) and all %values set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   key     -- Integer, optional: Key to initialize to. EG: quadrupole$, etc.
!   sub_key -- Integer, optional: Sub-key to initialize to.
!   ix_ele     -- Integer, optional: ix_ele index to initalize to. Default = -1.
!   ix_branch  -- Integer, optional: Branch index to initalize to. Default = 0.
!
! Output:
!   ele -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele, key, sub_key, ix_ele, ix_branch)

implicit none

type (ele_struct)  ele
integer, optional :: key, sub_key
integer, optional :: ix_branch, ix_ele

!

ele%type = ' '
ele%alias = ' '
ele%name = '<Initialized>'
ele%attribute_name = ' '

ele%key = integer_option (0, key)
ele%sub_key = integer_option (0, sub_key)

ele%value(:) = 0
ele%old_value(:) = 0
call init_coord (ele%map_ref_orb_in)
call init_coord (ele%map_ref_orb_out)

ele%lord_status = free$
ele%slave_status = free$
ele%ix_value = 0
ele%ic1_lord = 0
ele%ic2_lord = -1
ele%n_lord = 0
ele%ix1_slave = 0
ele%ix2_slave = -1
ele%n_slave = 0
ele%ix_pointer = 0
ele%s = 0
ele%ref_time = 0
ele%ix_branch = 0
ele%ix_ele = -1

if (present(ix_branch)) ele%ix_branch = ix_branch
if (present(ix_ele)) ele%ix_ele = ix_ele

call init_floor (ele%floor)

ele%mat6_calc_method = bmad_standard$
ele%tracking_method  = bmad_standard$
ele%field_calc       = bmad_standard$
ele%integrator_order = bmad_com%default_integ_order
ele%ref_orbit  = 0
ele%num_steps = 0

ele%is_on             = .true.
ele%multipoles_on     = .true.
ele%symplectify       = .false.
ele%map_with_offsets  = .true.
ele%on_a_girder       = .false.
ele%csr_calc_on       = .true.

ele%field_master = .false.
ele%coupler_at  = exit_end$
ele%aperture_at = exit_end$
ele%offset_moves_aperture = .false.

call deallocate_ele_pointers (ele)

! init Twiss

ele%c_mat = 0
ele%gamma_c = 1.0

ele%x%eta  = 0
ele%x%etap = 0

ele%y%eta  = 0
ele%y%etap = 0

ele%a%beta     = 0
ele%a%alpha    = 0
ele%a%gamma    = 0
ele%a%eta      = 0
ele%a%etap     = 0
ele%a%phi      = 0
ele%a%sigma    = 0
ele%a%emit     = 0

ele%b%beta     = 0
ele%b%alpha    = 0
ele%b%gamma    = 0
ele%b%eta      = 0
ele%b%etap     = 0
ele%b%phi      = 0
ele%b%sigma    = 0
ele%b%emit     = 0

ele%z%beta     = 0
ele%z%alpha    = 0
ele%z%gamma    = 0
ele%z%eta      = 0
ele%z%etap     = 0
ele%z%phi      = 0
ele%z%sigma    = 0
ele%z%emit     = 0

! This is needed because of a compiler and/or totalview bug

allocate (ele%r(1,1))
ele%r = 0.0

end subroutine init_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+ 
! Subroutine init_floor (floor)
!
! Routine to initialize a floor_position_struct to zero.
!
! Output:
!   floor -- Floor_position_struct: Floor coordinates to init.
!-

subroutine init_floor (floor)

implicit none

type (floor_position_struct) floor

!

floor%x = 0
floor%y = 0
floor%z = 0
floor%theta = 0
floor%phi   = 0
floor%psi   = 0

end subroutine init_floor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- Lat_struct: Lattice with element array.
!     %ele(:)      -- Ele_struct, pointer: Element array to reallocate.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!   ix_branch   -- Integer, optional: Branch index. Default is 0.
!
! Output:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:) -- Ele_struct, pointer: Resized element array.
!-

subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)

implicit none

type (lat_struct), target :: lat
integer, optional :: upper_bound
integer, optional :: ix_branch
integer i_branch

!

i_branch = integer_option (0, ix_branch)

if (i_branch == 0) then
  call allocate_ele_array (lat%ele, upper_bound)
  if (allocated(lat%branch)) lat%branch(0)%ele => lat%ele
else
  call allocate_ele_array (lat%branch(i_branch)%ele, upper_bound)
endif

end subroutine allocate_lat_ele_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_ele_array (ele, upper_bound)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
! Note: Do not use this routine for the lat%ele(:) array!
!   Use allocate_lat_ele_array instead.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele(:)      -- Ele_struct, pointer: Element array.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!
! Output:
!   ele(:)      -- Ele_struct, pointer: Element array.
!-

subroutine allocate_ele_array (ele, upper_bound)

implicit none

type (ele_struct), pointer :: ele(:)
type (ele_struct), pointer :: temp_ele(:)

integer, optional :: upper_bound
integer curr_ub, ub, i

! get new size

ub = 10
if (associated (ele)) ub = max (int(1.3*size(ele)), ub)
if (present(upper_bound))  ub = upper_bound

!  save ele if present

if (associated (ele)) then
  if (ub == ubound(ele, 1)) return
  curr_ub = min(ub, ubound(ele, 1))
  allocate (temp_ele(0:curr_ub))
  call transfer_eles (ele(0:curr_ub), temp_ele)
  do i = curr_ub+1, ubound(ele, 1)
    call deallocate_ele_pointers(ele(i))
  enddo
  deallocate (ele)
  allocate(ele(0:ub))
  call transfer_eles (temp_ele(0:curr_ub), ele(0:curr_ub))
  deallocate (temp_ele)
else
  curr_ub = -1
  allocate(ele(0:ub))
endif

! 

do i = curr_ub+1, ub
  call init_ele (ele(i))
  ele(i)%ix_ele = i
end do

end subroutine allocate_ele_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_branch_array (branch, upper_bound, lat)
!
! Subroutine to allocate or re-allocate an branch array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   branch(:)   -- Branch_struct, pointer: Branch array.
!   upper_bound -- Integer: Desired upper bound.
!   lat         -- Lat_struct, optional: Needed if branch is not allocated
! 
! Output:
!   branch(:)   -- Branch_struct, pointer: Branch array.
!-

subroutine allocate_branch_array (branch, upper_bound, lat)

implicit none

type (branch_struct), allocatable :: branch(:)
type (branch_struct), pointer :: temp_branch(:)
type (lat_struct), optional, target :: lat

integer :: upper_bound
integer curr_ub, ub, i

character(20) :: r_name = 'allocate_branch_array'

!  save branch if present

ub = upper_bound
if (allocated (branch)) then
  if (ub == ubound(branch, 1)) return
  curr_ub = min(ub, ubound(branch, 1))
  allocate (temp_branch(0:curr_ub))
  call transfer_branches (branch(0:curr_ub), temp_branch)
  do i = curr_ub+1, ubound(branch, 1)
    call deallocate_ele_array_pointers(branch(i)%ele)
    deallocate(branch(i)%n_ele_track)
    deallocate(branch(i)%n_ele_max)
  enddo
  deallocate (branch)
  allocate(branch(0:ub))
  call transfer_branches (temp_branch(0:curr_ub), branch(0:curr_ub))
  deallocate (temp_branch)
else
  curr_ub = -1
  allocate(branch(0:ub))
  if (.not. present(lat)) then
    call out_io (s_fatal$, r_name, 'LAT ARGUMENT MISSING.')
    call err_exit
  endif
  lat%branch(0)%ele => lat%ele
  lat%branch(0)%param => lat%param
  lat%branch(0)%n_ele_track => lat%n_ele_track
  lat%branch(0)%n_ele_max => lat%n_ele_max
endif

! 

do i = curr_ub+1, ub
  branch(i)%ix_branch = i
  if (i == 0) cycle
  allocate(branch(i)%n_ele_track)
  allocate(branch(i)%n_ele_max)
  allocate(branch(i)%param)
end do

end subroutine allocate_branch_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine deallocate_lat_pointers (lat)
!
! Subroutine to deallocate the pointers in a lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lat with pointers.
!
! Output:
!   lat -- lat_struct: Lat with deallocated pointers.
!-

subroutine deallocate_lat_pointers (lat)

implicit none

type (lat_struct) lat
integer i

!

if (associated (lat%ele)) then
  call deallocate_ele_array_pointers (lat%ele)
  call deallocate_ele_pointers (lat%ele_init)
  deallocate (lat%control)
  deallocate (lat%ic)
endif

call deallocate_branch (lat%branch)

lat%n_ele_track  = -1
lat%n_ele_max  = -1

end subroutine deallocate_lat_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_branch (branch)
!
! Subroutine to deallocate a branch array and everything in it.
!
! Modules needed:
!   use bmad
!
! Input:
!   branch(:) -- Branch_struct, allocatable: Array of lines.
!
! Output:
!   branch(:) -- Branch_struct, allocatable: Deallocated array.
!-

subroutine deallocate_branch (branch)

implicit none

type (branch_struct), allocatable :: branch(:)
integer i

!

if (allocated (branch)) then
  do i = 1, ubound(branch, 1)
    call deallocate_ele_array_pointers (branch(i)%ele)
    deallocate (branch(i)%param)
  enddo
  deallocate (branch)
endif

end subroutine deallocate_branch

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_ele_array_pointers (eles)
!
! Routine to deallocate the pointers of all the elements in an 
! element array and the array itself.
!
! Modules needed:
!   use bmad
!
! Input:
!   eles(:) -- Ele_struct, pointer: Array of elements.
!
! Output:
!   eles(:) -- Ele_struct, pointer: Deallocated array.
!-

subroutine deallocate_ele_array_pointers (eles)

implicit none

type (ele_struct), pointer :: eles(:)
integer i

!

do i = lbound(eles, 1), ubound(eles, 1)
  call deallocate_ele_pointers (eles(i))
enddo

deallocate (eles)

end subroutine deallocate_ele_array_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine transfer_mat_from_twiss (ele1, ele2, m)
!
! Subroutine to make a 6 x 6 transfer matrix from the twiss parameters
! at two points.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct: Element with twiss parameters for the starting point.
!     %a, %b -- a-mode and b-mode Twiss paramters
!       %beta   -- Beta parameter.
!       %alpha  -- Alpha parameter.
!       %phi    -- Phase at initial point.
!     %x  %y -- dispersion values
!       %eta    -- Dispersion at initial point.
!       %etap   -- Dispersion derivative at initial point.
!     %c_mat(2,2) -- Coupling matrix
!   ele2 -- Ele_struct: Element with twiss parameters for the ending point.
!
! Output:
!   m(6,6) -- Real(rp): Transfer matrix between the two points.
!-

subroutine transfer_mat_from_twiss (ele1, ele2, m)

implicit none

type (ele_struct) ele1, ele2

real(rp) m(6,6), v_mat(4,4), v_inv_mat(4,4), det
character(20) :: r_name = 'transfer_mat_from_twiss'

! Error check

if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele1%name)
  call err_exit
endif

if (ele2%a%beta == 0 .or. ele2%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele2%name)
  call err_exit
endif

! Transfer matrices without coupling or dispersion

call mat_make_unit (m)
call transfer_mat2_from_twiss (ele1%a, ele2%a, m(1:2,1:2))
call transfer_mat2_from_twiss (ele1%b, ele2%b, m(3:4,3:4))

! Add in coupling

if (any(ele1%c_mat /= 0)) then
  call mat_det (ele1%c_mat, det)
  ele1%gamma_c = sqrt(1-det)
  call make_v_mats (ele1, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (m(1:4,1:4), v_inv_mat)
endif

if (any(ele2%c_mat /= 0)) then
  call mat_det (ele2%c_mat, det)
  ele2%gamma_c = sqrt(1-det)
  call make_v_mats (ele2, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (v_mat, m(1:4,1:4))
endif

! Add in dispersion.

m(1:4,6) = (/ ele2%x%eta, ele2%x%etap, ele2%y%eta, ele2%y%etap /) - &
        matmul (m(1:4,1:4), (/ ele1%x%eta, ele1%x%etap, ele1%y%eta, ele1%y%etap /)) 

! The m(5,x) terms follow from the symplectic condition.

m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)


end subroutine transfer_mat_from_twiss

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine match_ele_to_mat6 (ele, vec0, mat6, err_flag)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Match element.
!     %value(beta_a0$) -- Beta_a at the start
!
! Output:
!   vec0(6)   -- Real(rp): Currently just set to zero.
!   mat6(6,6) -- Real(rp): Transfer matrix.
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

subroutine match_ele_to_mat6 (ele, vec0, mat6, err_flag)

implicit none

type (ele_struct), target :: ele, ele0, ele1

real(rp) mat6(6,6), vec0(6)
real(rp), pointer :: v(:)

logical err_flag

! Special case where match_end is set but there is no beginning beta value yet.
! In this case, just return the unit matrix and set the err_flag.

if (ele%value(match_end$) /= 0 .and. ele%value(beta_a0$) == 0) then
  call mat_make_unit (mat6)
  vec0 = 0
  err_flag = .true.
  return
endif

!

err_flag = .false.

v => ele%value

ele0%a%beta   = v(beta_a0$)
ele0%a%alpha  = v(alpha_a0$)
ele0%a%phi    = 0
ele0%x%eta    = v(eta_x0$)
ele0%x%etap   = v(etap_x0$)

ele0%b%beta   = v(beta_b0$)
ele0%b%alpha  = v(alpha_b0$)
ele0%b%phi    = 0
ele0%y%eta    = v(eta_y0$)
ele0%y%etap   = v(etap_y0$)

ele1%a%beta   = v(beta_a1$)
ele1%a%alpha  = v(alpha_a1$)
ele1%a%phi    = v(dphi_a$)
ele1%x%eta    = v(eta_x1$)
ele1%x%etap   = v(etap_x1$)

ele1%b%beta   = v(beta_b1$)
ele1%b%alpha  = v(alpha_b1$)
ele1%b%phi    = v(dphi_b$)
ele1%y%eta    = v(eta_y1$)
ele1%y%etap   = v(etap_y1$)

ele0%c_mat(1,:) = (/ v(c_11$), v(c_12$) /)
ele0%c_mat(2,:) = (/ v(c_21$), v(c_22$) /)
ele0%gamma_c    = v(gamma_c$)

ele1%c_mat = 0 
ele1%gamma_c = 1

ele0%name = ele%name
ele1%name = ele%name

call transfer_mat_from_twiss (ele0, ele1, mat6)

! Kick part

vec0 = (/ v(x1$), v(px1$), v(y1$), v(py1$), v(z1$), v(pz1$) /) - &
       matmul (mat6, (/ v(x0$), v(px0$), v(y0$), v(py0$), v(z0$), v(pz0$) /))

end subroutine match_ele_to_mat6

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_wake (wake_in, wake_out)

implicit none

type (wake_struct), pointer :: wake_in, wake_out
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

!

if (associated (wake_in)) then
  n_sr_table       = size(wake_in%sr_table)
  n_sr_mode_long  = size(wake_in%sr_mode_long)
  n_sr_mode_trans = size(wake_in%sr_mode_trans)
  n_lr        = size(wake_in%lr)
  call init_wake (wake_out, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
  wake_out%sr_file   = wake_in%sr_file
  wake_out%lr_file   = wake_in%lr_file
  wake_out%z_sr_mode_max  = wake_in%z_sr_mode_max
  wake_out%sr_table       = wake_in%sr_table
  wake_out%sr_mode_long  = wake_in%sr_mode_long
  wake_out%sr_mode_trans = wake_in%sr_mode_trans
  wake_out%lr        = wake_in%lr
else
  if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0, 0)
endif

end subroutine transfer_wake

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
!
! Subroutine to initialize a wake struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_sr_table       -- Integer: Number of terms: wake%sr_table(0:n_sr-1).
!   n_sr_mode_long  -- Integer: Number of terms: wake%nr(n_sr_mode_long).
!   n_sr_mode_trans -- Integer: Number of terms: wake%nr(n_sr_mode_trans).
!   n_lr        -- Integer: Number of terms: wake%nr(n_lr)
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!               If all inputs are 0 then wake is deallocated.
!-

subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

implicit none

type (wake_struct), pointer :: wake
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

! Deallocate wake if all inputs are zero.

if (n_sr_table == 0 .and. n_sr_mode_long == 0 .and. n_sr_mode_trans == 0 .and. n_lr == 0) then
  if (associated(wake)) then
    deallocate (wake%sr_table)
    deallocate (wake%sr_mode_long)
    deallocate (wake%sr_mode_trans)
    deallocate (wake%lr)
    deallocate (wake)
  endif
  return
endif

!

if (associated (wake)) then
  if (size(wake%sr_table) /= n_sr_table) then
    deallocate (wake%sr_table)
    allocate (wake%sr_table(0:n_sr_table-1))
  endif
  if (size(wake%sr_mode_long) /= n_sr_mode_long) then
    deallocate (wake%sr_mode_long)
    allocate (wake%sr_mode_long(n_sr_mode_long))
  endif
  if (size(wake%sr_mode_trans) /= n_sr_mode_trans) then
    deallocate (wake%sr_mode_trans)
    allocate (wake%sr_mode_trans(n_sr_mode_trans))
  endif
  if (size(wake%lr) /= n_lr) then
    deallocate (wake%lr)
    allocate (wake%lr(n_lr))
  endif

else
  allocate (wake)
  allocate (wake%sr_table(0:n_sr_table-1))
  allocate (wake%sr_mode_long(n_sr_mode_long))
  allocate (wake%sr_mode_trans(n_sr_mode_trans))
  allocate (wake%lr(n_lr))
endif

end subroutine init_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine calc_superimpose_key (ele1, ele2) result (ele3)
!
! Function to decide what ele3%key and ele3%sub_key should be
! when two elements, ele1, and ele2, are superimposed.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct:
!     %key
!     %sub_key
!   ele2 -- Ele_struct:
!     %key
!     %sub_key
!
! Output:
!   ele3 -- Ele_struct:
!     %key
!     %sub_key
!-

subroutine calc_superimpose_key (ele1, ele2, ele3)

implicit none

type (ele_struct), target :: ele1, ele2, ele3
integer key1, key2
integer, pointer :: key3

!

key1 = ele1%key
key2 = ele2%key
key3 => ele3%key

key3 = -1  ! Default if no superimpse possible
ele3%sub_key = 0

if (key1 == key2) then
  if (key1 == wiggler$ .and. ele1%sub_key /= ele2%sub_key) return  ! Bad combo
  key3 = key1
  if (key1 == wiggler$) ele3%sub_key = ele1%sub_key
  return
endif

if (key1 == drift$) then
  key3 = key2
  ele3%sub_key = ele2%sub_key
  return
endif

if (key2 == drift$) then
  key3 = key1
  ele3%sub_key = ele1%sub_key
  return
endif

if (any(key1 == (/ rcollimator$, monitor$, instrument$ /))) then
  key3 = key2
  return
endif

if (any(key2 == (/ rcollimator$, monitor$, instrument$ /))) then
  key3 = key1
  return
endif

if (any(key1 == (/ kicker$, hkicker$, vkicker$ /))) then
  if (any(key2 == (/ kicker$, hkicker$, vkicker$ /))) then
    key3 = kicker$
  else
    key3 = key2
  endif
  return
endif

if (any(key2 == (/ kicker$, hkicker$, vkicker$ /))) then
  key3 = key1
  return
endif

!

select case (key1)

case (quadrupole$,  solenoid$, sol_quad$) 
  select case (key2)
  case (quadrupole$);    key3 = sol_quad$
  case (solenoid$);      key3 = sol_quad$
  case (sol_quad$);      key3 = sol_quad$
  case (bend_sol_quad$); key3 = bend_sol_quad$
  case (sbend$);         key3 = bend_sol_quad$
  end select

case (bend_sol_quad$)
  select case (key2)
  case (quadrupole$);    key3 = bend_sol_quad$
  case (solenoid$);      key3 = bend_sol_quad$
  case (sol_quad$);      key3 = bend_sol_quad$
  case (sbend$);         key3 = bend_sol_quad$
  end select
end select

end subroutine calc_superimpose_key

end module
