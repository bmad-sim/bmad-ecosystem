!+
! Subroutine write_lattice_scibmad_format(scibmad_file, lat, err_flag)
!
! Routine to create a SciBmad lattice file.
!
! Input:
!   scibmad_file  -- character(*): SciBmad lattice file name.
!   lat           -- lat_struct: Lattice
!
! Output:
!   err_flag      -- logical, optional: Set True if there is a problem. That is, if the file could
!                     not be opened or if the lattice contains something that could not be
!                     translated. Set False otherwise.
!-

subroutine write_lattice_scibmad_format(scibmad_file, lat, err_flag)

use write_lattice_file_mod, dummy => write_lattice_scibmad_format
use bmad_routine_interface, dummy2 => write_lattice_scibmad_format
use expression_mod
use taylor_mod, only: mat6_to_taylor
use super_recipes_mod, only: super_sort

implicit none

! There is one this_expr_struct for each controlled (SciBmad element name, SciBmad attribute) pair.
! Expressions are collected by SciBmad attribute since a single Bmad attribute (EG: HKICK of a tilted
! element) may map to multiple SciBmad attributes and multiple Bmad attributes may map to a single
! SciBmad attribute.
!
! Control_type component of this_expr_struct:
!   control_lord$   - Element is a group or overlay
!   group$          - Non-controller element with group control of attribute
!   overlay$        - Non-controller element with overlay control of attribute

type this_expr_struct
  integer :: control_type = not_set$
  type (ele_struct), pointer :: bmad_ele => null()
  real(rp) :: base_value = 0          ! Present value of the SciBmad attribute. Used with groups.
  real(rp) :: sum_ctl = 0             ! Present value of the overlay controlled part of the SciBmad attribute.
  integer :: ix_attrib_counted(20) = 0 ! Bmad attributes already counted in sum_ctl.
  integer :: n_attrib_counted = 0
  character(100) :: sort_name = ''
  character(100) :: scibmad_ele = ''
  character(40) :: scibmad_attrib = ''
  character(5000) :: expr = ''
  character(100), allocatable :: group_var_names(:)
  real(rp), allocatable :: group_var_values(:)
end type

type (lat_struct), target :: lat, lat2
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele2, lord, slave, slave2, multi_lord
type (coord_struct), pointer :: orb
type (multipass_region_lat_struct), target :: mult_lat
type (multipass_all_info_struct), target :: m_info
type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele
type (multipass_ele_info_struct), pointer :: e_info
type (ele_pointer_struct), allocatable :: named_eles_ptr(:)  ! List of unique element names 
type (lat_ele_order_struct) order
type (ele_attribute_struct) info
type (taylor_struct) taylor(6), spin_taylor(0:3)
type (nametable_struct) var_nametab
type (control_struct), pointer :: ctl
type (this_expr_struct), allocatable, target :: expr(:)
type (this_expr_struct), pointer :: e_ptr

real(rp) factor(2), tot, resid

integer n, i, j, k, ix, ib, ie, iu, is, it, iv, n_names, ix_match, ix_pass, ix_r, ios, n_expr, ix_expr
integer ix_lord, ix_super, ie1, ib1, eles_not_translated(20)
integer n_done
integer, allocatable :: an_indexx(:), expr_index(:), done_index(:)

logical has_been_added, in_multi_region, is_group, is_mult
logical has_planar_wiggler, is_added
logical xlate_err    ! Set True if something in the lattice cannot be translated.
logical, optional :: err_flag

character(*) scibmad_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(200) name2
character(100) name, look_for, ele_name, sort_name
character(40) sci_attrib(2)
character(40), allocatable :: scibmad_names(:)
character(200), allocatable :: done_names(:)
character(240) fname
character(4000) line
character(*), parameter :: r_name = 'write_lattice_scibmad_format'
character(20) :: scibmad_ele_type(n_key$)

! err_flag is set False only after the file has been successfully written.

if (present(err_flag)) err_flag = .true.
xlate_err = .false.

scibmad_ele_type(drift$)                = 'Drift'
scibmad_ele_type(sbend$)                = 'SBend'
scibmad_ele_type(quadrupole$)           = 'Quadrupole'
scibmad_ele_type(group$)                = '??Group'
scibmad_ele_type(sextupole$)            = 'Sextupole'
scibmad_ele_type(overlay$)              = '??Overlay'
scibmad_ele_type(custom$)               = 'Marker'
scibmad_ele_type(taylor$)               = 'LineElement'
scibmad_ele_type(rfcavity$)             = 'RFCavity'
scibmad_ele_type(elseparator$)          = 'Drift'           !!! Not translated
scibmad_ele_type(beambeam$)             = 'BeamBeam'
scibmad_ele_type(wiggler$)              = 'LineElement'   ! SciBmad has no Wiggler constructor. Uses kind = "Wiggler".
scibmad_ele_type(sol_quad$)             = 'Solenoid'
scibmad_ele_type(marker$)               = 'Marker'
scibmad_ele_type(kicker$)               = 'Kicker'
scibmad_ele_type(hybrid$)               = 'Marker'
scibmad_ele_type(octupole$)             = 'Octupole'
scibmad_ele_type(rbend$)                = 'SBend'
scibmad_ele_type(multipole$)            = 'Multipole'
scibmad_ele_type(ab_multipole$)         = 'Multipole'
scibmad_ele_type(solenoid$)             = 'Solenoid'
scibmad_ele_type(patch$)                = 'Patch'
scibmad_ele_type(lcavity$)              = 'RFCavity'
scibmad_ele_type(null_ele$)             = 'NullEle'
scibmad_ele_type(beginning_ele$)        = 'Marker'
scibmad_ele_type(match$)                = 'LineElement'
scibmad_ele_type(monitor$)              = 'Drift'
scibmad_ele_type(instrument$)           = 'Drift'
scibmad_ele_type(hkicker$)              = 'Kicker'
scibmad_ele_type(vkicker$)              = 'Kicker'
scibmad_ele_type(rcollimator$)          = 'Drift'
scibmad_ele_type(ecollimator$)          = 'Drift'
scibmad_ele_type(girder$)               = 'Nothing!'
scibmad_ele_type(converter$)            = 'Converter'
scibmad_ele_type(photon_fork$)          = 'Marker'      !!! Not translated
scibmad_ele_type(fork$)                 = 'Marker'      !!! Not translated
scibmad_ele_type(mirror$)               = 'Marker'      !!! Not translated
scibmad_ele_type(crystal$)              = 'Marker'      !!! Not translated
scibmad_ele_type(pipe$)                 = 'Drift'
scibmad_ele_type(capillary$)            = 'Drift'       !!! Not translated
scibmad_ele_type(multilayer_mirror$)    = 'Drift'       !!! Not translated
scibmad_ele_type(e_gun$)                = 'EGun'
scibmad_ele_type(em_field$)             = 'EMField'
scibmad_ele_type(floor_shift$)          = 'FloorShift'
scibmad_ele_type(fiducial$)             = 'Fiducial'
scibmad_ele_type(undulator$)            = 'LineElement'   ! SciBmad has no Undulator constructor. Uses kind = "Wiggler".
scibmad_ele_type(diffraction_plate$)    = 'Marker'
scibmad_ele_type(photon_init$)          = 'Marker'
scibmad_ele_type(sample$)               = 'Marker'
scibmad_ele_type(detector$)             = 'Marker'
scibmad_ele_type(sad_mult$)             = 'Marker'
scibmad_ele_type(mask$)                 = 'Marker'
scibmad_ele_type(ac_kicker$)            = 'Marker'
scibmad_ele_type(lens$)                 = 'Marker'
scibmad_ele_type(crab_cavity$)          = 'CrabCavity'
scibmad_ele_type(ramper$)               = 'Ramper'
scibmad_ele_type(rf_bend$)              = 'RFBend'
scibmad_ele_type(gkicker$)              = 'Kicker'
scibmad_ele_type(foil$)                 = 'Marker'
scibmad_ele_type(thick_multipole$)      = 'ThickMultipole'
scibmad_ele_type(pickup$)               = 'Drift'
scibmad_ele_type(feedback$)             = 'Drift'
scibmad_ele_type(fixer$)                = 'Fixer'

eles_not_translated = -1
eles_not_translated(1:18) = [elseparator$, photon_fork$, fork$, mirror$, crystal$, diffraction_plate$, photon_init$, &
                           sample$, detector$, sad_mult$, mask$, ac_kicker$, lens$, foil$, pickup$, feedback$, hybrid$, custom$]

! Elements with the same name but different SciBmad definitions are given unique names.

lat2 = lat
call this_create_unique_ele_names(lat2, '_n?')

! Open file

call fullfilename(scibmad_file, fname)
iu = lunget()
open (iu, file = fname, status = 'unknown', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE FOR WRITING: ' // trim(fname))
  return
endif

write (iu, '(4a)') '# Translated using Bmad based Bmad-to-SciBmad translation code.'
write (iu, '(4a)') '# Translated from Bmad lattice file: ', trim(lat%input_file_name)
write (iu, '(a)')
write (iu, '(a)')  'using Beamlines'
write (iu, '(a)')  'using BeamTracking'

! Write the four-potential function used by wiggler and undulator elements.

has_planar_wiggler = .false.

do ib = 0, ubound(lat2%branch, 1)
  branch => lat2%branch(ib)
  do ie = 1, branch%n_ele_track
    if (is_planar_wiggler(branch%ele(ie))) has_planar_wiggler = .true.
  enddo
enddo

if (has_planar_wiggler) call write_planar_wiggler_four_potential(iu)

! Write functions for Taylor elements

do ib = 0, ubound(lat2%branch, 1)
  branch => lat2%branch(ib)
  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)

    if (any(ele%key == eles_not_translated)) then
      call out_io(s_warn$, r_name, 'Element translation problem for ' // ele_full_name(ele) // ' of type: ' // key_name(ele%key), &
                                  '   Will translate to a: ' // scibmad_ele_type(ele%key))
    endif

    select case (ele%key)
    case (match$)
      call mat6_to_taylor(ele%vec0, ele%mat6, taylor)
      call write_this_taylor(iu, ele, taylor)
      cycle

    case (taylor$)
      call write_this_taylor(iu, ele, ele%taylor)
      cycle
    end select

  enddo
enddo

! Write element defs

! Note: Beamlines cannot currently handle multipass nor superimpose so ignore.
! Stuff that is commented out due to this is marked by "!!!"

n_names = 0
n = lat2%n_ele_max
allocate (scibmad_names(n), an_indexx(n), named_eles_ptr(n))

write (iu, '(a)')
write (iu, '(a)') '@elements begin'

do ib = 0, ubound(lat2%branch, 1)
  branch => lat2%branch(ib)
  ele_loop: do ie = 0, branch%n_ele_track   !!! Note: Not n_ele_max since superimpose/multipass not handled
    ele => branch%ele(ie)
    ele_name = scibmad_ele_name(ele%name, ib)

    if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$ .or. ele%key == girder$) cycle   ! Not currently handled
    if (ele%key == null_ele$) cycle

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, scibmad_names, an_indexx, n_names, ix_match, has_been_added, &
                                                                  named_eles_ptr, scibmad_ele_name(ele%name, ib))
    if (.not. has_been_added) cycle

    ! Write element def

    call ele_def_line(ele, ele_name, line, .true.)
    call write_lat_line(line, iu, .true., ampersand_at_ends = .false.)

  enddo ele_loop
enddo

write (iu, '(a)') 'end    # @elements'

!------------------------------------------------------------------------------------------------------
! Overlay and group elements....

write (iu, '(a)') '#---------------------------------------------------------------------------------------'
write (iu, '(a)') '# Overlay and Group elements'
write (iu, '(a)')

! Make a list of controlled attributes.
! Note: A single Bmad attribute (EG: HKICK of a tilted element) may map to multiple SciBmad
! attributes and multiple Bmad attributes may map to a single SciBmad attribute. So expressions
! are collected by (SciBmad element, SciBmad attribute).

allocate (expr(2*lat2%n_control_max), expr_index(2*lat2%n_control_max))
allocate (done_names(lat2%n_control_max), done_index(lat2%n_control_max))
n_expr = 0; n_done = 0

do ie = lat2%n_ele_track+1, lat2%n_ele_max
  lord => lat2%ele(ie)

  if (lord%key == girder$) then
    xlate_err = .true.
    cycle
  endif

  if (lord%key /= overlay$ .and. lord%key /= group$) cycle

  do is = 1, lord%n_slave
    slave => pointer_to_slave(lord, is, ctl)
    if (.not. allocated(ctl%stack)) cycle   ! Knot point control. Message given below.

    ! If the lord controls multiple elements with the same name (EG: a family of quadrupoles that all
    ! have the same definition), there is only one SciBmad element so only count the control once.

    name2 = int_str(ie) // ':' // trim(slave%name) // ':' // ctl%attribute
    call find_index(name2, done_names, done_index, n_done, ix, add_to_list = .true., has_been_added = is_added)
    if (.not. is_added) cycle

    call scibmad_attrib_name(ctl%attribute, slave, n, sci_attrib, factor)
    is_group = (lord%key == group$)

    do j = 1, n
      sort_name = trim(slave%name) // ':' // sci_attrib(j)
      call find_index(sort_name, expr%sort_name, expr_index, n_expr, ix_expr, add_to_list = .true., has_been_added = is_added)
      e_ptr => expr(ix_expr)
      if (is_added) then
        e_ptr%bmad_ele => slave
        e_ptr%scibmad_ele = scibmad_ele_name(slave%name)
        e_ptr%scibmad_attrib = sci_attrib(j)
        if (slave%key == group$ .or. slave%key == overlay$) then
          e_ptr%control_type = control_lord$
        elseif (is_group) then
          e_ptr%control_type = group$
        else
          e_ptr%control_type = overlay$
        endif
        e_ptr%base_value = scibmad_multipole_value(slave, sci_attrib(j), is_mult)
        if (.not. is_mult) e_ptr%base_value = factor(j) * value_of_attribute(slave, ctl%attribute)
      endif

      ! Present value of the overlay controlled part. If multiple lords control a given Bmad attribute,
      ! the attribute value is the sum of the contributions of all the lords so only count it once.

      if (.not. is_group .and. ctl%ix_attrib > 0 .and. ctl%ix_attrib <= num_ele_attrib$) then
        if (all(e_ptr%ix_attrib_counted(1:e_ptr%n_attrib_counted) /= ctl%ix_attrib) .and. &
                                  e_ptr%n_attrib_counted < size(e_ptr%ix_attrib_counted)) then
          e_ptr%n_attrib_counted = e_ptr%n_attrib_counted + 1
          e_ptr%ix_attrib_counted(e_ptr%n_attrib_counted) = ctl%ix_attrib
          e_ptr%sum_ctl = e_ptr%sum_ctl + factor(j) * slave%value(ctl%ix_attrib)
        endif
      endif

      if (is_group) then
        line = '((' // trim(expression_kernel(ctl%stack, lord, .false.)) // ') - (' // &
                       trim(expression_kernel(ctl%stack, lord, .true., e_ptr)) // '))'
      else
        line = '(' // trim(expression_kernel(ctl%stack, lord, .false.)) // ')'
      endif

      if (factor(j) /= 1.0_rp) line = re_str(factor(j)) // ' * ' // trim(line)
      e_ptr%expr = trim(e_ptr%expr) // ' + ' // trim(line)
    enddo
  enddo
enddo

! First print constants used in expressions.

call nametable_init(var_nametab)
write (iu, '(a)') 'c1 = Context('

do ie = lat2%n_ele_track+1, lat2%n_ele_max
  lord => lat2%ele(ie)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle

  do is = 1, lord%n_slave
    slave => pointer_to_slave(lord, is, ctl)
    if (.not. allocated(ctl%stack)) then
      call out_io(s_warn$, r_name, ele_full_name(lord) // ' Uses knot points for the control curve. This cannot yet be translated!')
      xlate_err = .true.
      exit
    endif

    do k = 1, size(ctl%stack)
      if (ctl%stack(k)%type /= variable$) cycle
      call find_index(ctl%stack(k)%name, var_nametab, ix_match, add_to_list = .true., has_been_added = is_added)
      if (is_added) write (iu, '(6x, 2a, es24.17, a)') trim(scibmad_ele_name(ctl%stack(k)%name)), ' = ', ctl%stack(k)%value, ','
    enddo
  enddo
enddo

! Output controller vars 

do ie = lat2%n_ele_track+1, lat2%n_ele_max
  lord => lat2%ele(ie)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle
  do iv = 1, size(lord%control%var)
    sort_name = trim(lord%name) // ':' // lord%control%var(iv)%name
    call find_index(sort_name, expr%sort_name, expr_index, n_expr, ix_expr)
    name = trim(scibmad_ele_name(lord%name)) // '_' // downcase(lord%control%var(iv)%name)
    if (ix_expr == 0) then
      write (iu, '(6x, 2a, es24.16, a)') trim(name), ' = ', lord%control%var(iv)%value, ','
    else
      write (iu, '(6x, 6a)') trim(name), ' = ', trim(def_expr(expr(ix_expr), .true.)), ','
    endif
  enddo
enddo

! Output associated (old) group controller parameters

do iv = 1, n_expr
  e_ptr => expr(iv)
  if (e_ptr%control_type /= group$) cycle

  write (iu, '(6x, a)') trim(e_ptr%scibmad_ele) // '_' // trim(e_ptr%scibmad_attrib) // ' = ' // re_str(e_ptr%base_value) // ','
  if (.not. allocated(e_ptr%group_var_names)) cycle

  do i = 1, size(e_ptr%group_var_names)
    write (iu, '(6x, 7a)') 'old_', trim(e_ptr%group_var_names(i)), '__', trim(e_ptr%scibmad_ele), '_', &
                              trim(e_ptr%scibmad_attrib), ' = ' // re_str(e_ptr%group_var_values(i)) // ','
  enddo
enddo

write (iu, '(6x, a)') ')'    ! End context construct
write (iu, '(a)')


! Now output deferred expressions.
! First: do overlay controlled parameters.
! A SciBmad multipole component may have contributions from Bmad attributes that are not controlled
! (EG: The bend angle contribution to Kn0). Such contributions are constant so just add them in.

do iv = 1, n_expr
  e_ptr => expr(iv)
  if (e_ptr%control_type /= overlay$) cycle

  line = def_expr(e_ptr, .false.)
  tot = scibmad_multipole_value(e_ptr%bmad_ele, e_ptr%scibmad_attrib, is_mult)
  if (is_mult) then
    resid = tot - e_ptr%sum_ctl
    if (abs(resid) > 1e-14_rp * max(abs(tot), abs(e_ptr%sum_ctl))) line = re_str(resid) // ' + ' // trim(line)
  endif

  write (iu, '(5a)') trim(e_ptr%scibmad_ele), '.', trim(e_ptr%scibmad_attrib), ' = DefExpr(c -> ', trim(line) // ')'
enddo

! Second: do group controlled parameters.
! A group varies an attribute incrementally: When a group variable is changed, the change in the
! control expression is added to the present attribute value. To do this, the present attribute value
! and the group variable values at the time of the last evaluation ("old" values) are stored in the context.

do iv = 1, n_expr
  e_ptr => expr(iv)
  if (e_ptr%control_type /= group$) cycle

  write (iu, '(5a)') trim(e_ptr%scibmad_ele), '.', trim(e_ptr%scibmad_attrib), ' = DefExpr(c ->'
  write (iu, '(12x, 1a)') 'begin'
  write (iu, '(14x, 7a)') 'result = c.', trim(e_ptr%scibmad_ele), '_', trim(e_ptr%scibmad_attrib), ' + ', trim(def_expr(e_ptr, .false.))
  if (allocated(e_ptr%group_var_names)) then
    do k = 1, size(e_ptr%group_var_names)
      name = 'c.old_' // trim(e_ptr%group_var_names(k)) // '__' // trim(e_ptr%scibmad_ele) // '_' // trim(e_ptr%scibmad_attrib)
      write (iu, '(14x, 3a)') trim(name), ' = c.', trim(e_ptr%group_var_names(k))
    enddo
  endif
  write (iu, '(14x, 7a)') 'c.', trim(e_ptr%scibmad_ele), '_', trim(e_ptr%scibmad_attrib), ' = result'
  write (iu, '(14x, 1a)') 'return result'
  write (iu, '(12x, 1a)') 'end)'
enddo

!------------------------------
! Define Branches.

do ib = 0, ubound(lat2%branch, 1)
  branch => lat2%branch(ib)

  write (iu, '(a)')
  branch%name = downcase(branch%name)
  if (branch%name == '') branch%name = 'branch' // int_str(ib+1)
  line = trim(branch%name) // ' = Branch(['

  do ie = 0, branch%n_ele_track
    ele => branch%ele(ie)
    call write_scibmad_element (line, iu, ele)
  enddo
 
  line = line(:len_trim(line)-1) // ']; name = ' // quote(branch%name) // ')'
  call write_lat_line (line, iu, .true., ampersand_at_ends = .false.)
enddo

! Define Lattice

line = 'lat = Lattice(['
do ib = 0, ubound(lat2%branch, 1)
  branch => lat2%branch(ib)
  if (branch%ix_from_branch > -1) cycle
  line = trim(line) // ', ' // branch%name
enddo

ix = index(line, '[, ')
line = line(:ix) // trim(line(ix+3:)) // '], context = c1)'
write (iu, '(a)')
write (iu, '(a)') trim(line)

! cleanup

close(iu)
deallocate (scibmad_names, an_indexx)
!!! deallocate (mult_lat%branch)

if (present(err_flag)) err_flag = xlate_err

!----------------------------------------------------------------------------------------------
contains

function jbool(logic) result (bool_str)
logical logic
character(5) bool_str

if (logic) then
  bool_str = 'true'
else
  bool_str = 'false'
endif

end function jbool

!----------------------------------------------------------------------------------------------
! contains
!
! Returns True if ele is a wiggler or undulator that uses Bmad's periodic planar model with kx = 0.
! This is the only wiggler model that can currently be translated.

function is_planar_wiggler(ele) result (is_planar)

type (ele_struct) ele
type (ele_struct), pointer :: ele2
logical is_planar

!

is_planar = .false.
if (ele%key /= wiggler$ .and. ele%key /= undulator$) return
if (ele%value(kx$) /= 0) return
if (ele%value(l_period$) == 0 .or. ele%value(l$) == 0) return
ele2 => pointer_to_field_ele(ele, 1)
if (ele2%field_calc /= planar_model$) return
is_planar = .true.

end function is_planar_wiggler

!----------------------------------------------------------------------------------------------
! contains
!
! Write the Julia function that gives the four-potential of Bmad's periodic planar wiggler model.
! This function must be present in the translated lattice file since it is not part of SciBmad.

subroutine write_planar_wiggler_four_potential (iu)

integer iu

!

write (iu, '(a)')
write (iu, '(a)') '# Four-potential, and its derivatives, of the Bmad periodic planar wiggler/undulator model with kx = 0.'
write (iu, '(a)') '# params = (B_max [Tesla], k_w [1/m], phase [rad]). The gauge used is phi = Ay = As = 0 with'
write (iu, '(a)') '#   Ax = (B_max / k_w) * cosh(k_w*y) * sin(k_w*s + phase)'
write (iu, '(a)') '# which gives the Bmad field'
write (iu, '(a)') '#   Bx = 0,  By = B_max*cosh(k_w*y)*cos(k_w*s + phase),  Bs = -B_max*sinh(k_w*y)*sin(k_w*s + phase).'
write (iu, '(a)') '# The potential is in physical units so four_potential_normalized must be false.'
write (iu, '(a)') '# The derivative tuple order is:'
write (iu, '(a)') '#   (dphi/dx, dphi/dy, dphi/ds, dphi/dt, dAx/dx, dAx/dy, dAx/ds, dAx/dt,'
write (iu, '(a)') '#    dAy/dx,  dAy/dy,  dAy/ds,  dAy/dt,  dAs/dx, dAs/dy, dAs/ds, dAs/dt)'
write (iu, '(a)') '# Note: Bmad measures z with respect to a reference particle that follows the wiggling on-axis'
write (iu, '(a)') '# trajectory while SciBmad uses a straight line reference. The z of the two codes will therefore'
write (iu, '(a)') '# differ by the on-axis path lengthening of the wiggler.'
write (iu, '(a)')
write (iu, '(a)') '@inline function planar_wiggler_four_potential(x, y, s, t, params)'
write (iu, '(a)') '  B_max, k_w, phase = params'
write (iu, '(a)') '  theta = k_w * s + phase'
write (iu, '(a)') '  A_x = (B_max / k_w) * cosh(k_w * y) * sin(theta)'
write (iu, '(a)') '  dA_x_dy = B_max * sinh(k_w * y) * sin(theta)'
write (iu, '(a)') '  dA_x_ds = B_max * cosh(k_w * y) * cos(theta)'
write (iu, '(a)') '  z = zero(A_x)'
write (iu, '(a)') '  potential = (z, A_x, z, z)'
write (iu, '(a)') '  derivatives = (z, z, z, z,'
write (iu, '(a)') '                 z, dA_x_dy, dA_x_ds, z,'
write (iu, '(a)') '                 z, z, z, z,'
write (iu, '(a)') '                 z, z, z, z)'
write (iu, '(a)') '  return potential, derivatives'
write (iu, '(a)') 'end'

end subroutine write_planar_wiggler_four_potential

!----------------------------------------------------------------------------------------------
! contains
!
! Construct the SciBmad element definition line for ele.
! If warn = False, no warnings are issued and xlate_err is not touched. This is used when checking
! if elements with the same name have the same definition.

subroutine ele_def_line(ele, ele_name, line, warn)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele2

real(rp) f, length, k_wig, n_per, phase
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer i, j, n, ix, n_step, i_order, n_wig

character(*) ele_name, line
character(1) prefix
logical warn

!

length = ele%value(l$)

line = '  ' // trim(ele_name) // ' = ' // trim(scibmad_ele_type(ele%key)) // '('

if (ele%ix_ele == 0) then
  line = trim(line) // ', pc_ref = ' // re_str(ele%value(p0c$))
  line = trim(line) // ', species_ref = Species(' // quote(openpmd_species_name(ele%ref_species)) // ')'
  !! if (ele%a%beta /= 0) line = trim(line) // ', beta_a = ' // re_str(ele%a%beta)
  !! if (ele%b%beta /= 0) line = trim(line) // ', beta_b = ' // re_str(ele%b%beta)
  !! if (ele%a%alpha /= 0) line = trim(line) // ', alpha_a = ' // re_str(ele%a%alpha)
  !! if (ele%b%alpha /= 0) line = trim(line) // ', alpha_b = ' // re_str(ele%b%alpha)
  !! if (ele%x%eta /= 0) line = trim(line) // ', eta_x = ' // re_str(ele%x%eta)
  !! if (ele%y%eta /= 0) line = trim(line) // ', eta_y = ' // re_str(ele%y%eta)
  !! if (ele%x%etap /= 0) line = trim(line) // ', etap_x = ' // re_str(ele%x%etap)
  !! if (ele%y%etap /= 0) line = trim(line) // ', etap_y = ' // re_str(ele%y%etap)
  !! if (any(ele%c_mat /= 0)) line = trim(line) // ', c_mat = [' // re_str(ele%c_mat(1,1)) // ', ' // re_str(ele%c_mat(1,2)) // &
  !!                                                             '; ' // re_str(ele%c_mat(2,1)) // ', ' // re_str(ele%c_mat(2,2)) // ']'
  !! orb => lat2%particle_start
  !! if (any(orb%vec /= 0)) line = trim(line) // ', particle.orbit = [' // re_str(orb%vec(1)) // ', ' // re_str(orb%vec(2)) // ', ' // &
  !!                   re_str(orb%vec(3)) // ', ' // re_str(orb%vec(4)) // ', ' // re_str(orb%vec(5)) // ', ' // re_str(orb%vec(6)) // ']'
  !! if (any(orb%spin /= 0)) line = trim(line) // ', particle.spin = [' // &
  !!                                        re_str(orb%spin(1)) // ', ' // re_str(orb%spin(2)) // ', ' //re_str(orb%spin(3)) // ']'

endif

if (.not. ele%is_on) write (line, '(3a)') trim(line), ', is_on = ', jbool(ele%is_on)

!

if (ele%key == sbend$) then
  line = trim(line) // ', L = ' // re_str(length)
  if (ele%value(e1$) /= 0) line = trim(line) // ', e1 = ' // re_str(ele%value(e1$))
  if (ele%value(e2$) /= 0) line = trim(line) // ', e2 = ' // re_str(ele%value(e2$))

  if (ele%value(g$) /= 0)  line = trim(line) // ', g_ref = ' // re_str(ele%value(g$))
  if (ele%value(ref_tilt$) /= 0)  line = trim(line) // ', tilt_ref = ' // re_str(ele%value(ref_tilt$))
  if (ele%value(roll$) /= 0)  line = trim(line) // ', roll = ' // re_str(ele%value(roll$))
  !!! if (ele%value(fint$)*ele%value(hgap$) /= 0)    line = trim(line) // ', edge_int1 = ' // re_str(ele%value(fint$)*ele%value(hgap$))
  !!! if (ele%value(fintx$)*ele%value(hgapx$) /= 0)  line = trim(line) // ', edge_int2 = ' // re_str(ele%value(fintx$)*ele%value(hgapx$))
  if (ele%value(fint$)*ele%value(hgap$) /= 0 .or. ele%value(fintx$)*ele%value(hgapx$) /= 0) then
    call out_io(s_warn$, r_name, 'BEND EDGE_INT PARAMETER CANNOT YET BE TRANSLATED!')
    xlate_err = .true.
  endif

elseif (has_attribute(ele, 'L')) then
  if (length /= 0) line = trim(line) // ', L = ' // re_str(length)
endif

! Magnetic multipoles

call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, magnetic$, include_kicks$)
if (ele%key == sbend$) then 
  b_pole(0) = b_pole(0) + ele%value(angle$)
  ix = max(0, ix)
endif

if (ele%field_master) then
  f = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  prefix = 'B'
else
  f = 1
  prefix = 'K'
endif

if (length /= 0) f = f / length

do j = 0, ix
  if (length == 0) then
    if (a_pole(j) /= 0) line = trim(line) // ', ' // prefix // 's' // int_str(j) // 'L = ' // re_str(f * factorial(j) * a_pole(j))
    if (b_pole(j) /= 0) line = trim(line) // ', ' // prefix // 'n' // int_str(j) // 'L = ' // re_str(f * factorial(j) * b_pole(j))
  else
    if (a_pole(j) /= 0) line = trim(line) // ', ' // prefix // 's' // int_str(j) // ' = ' // re_str(f * factorial(j) * a_pole(j))
    if (b_pole(j) /= 0) line = trim(line) // ', ' // prefix // 'n' // int_str(j) // ' = ' // re_str(f * factorial(j) * b_pole(j))
  endif
enddo

! Electric multipoles

call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, electric$, include_kicks$)

do j = 0, ix
  if (a_pole(j) /= 0) line = trim(line) // ', Es' // int_str(j) // ' = ' // re_str(factorial(j) * a_pole(j))
  if (b_pole(j) /= 0) line = trim(line) // ', En' // int_str(j) // ' = ' // re_str(factorial(j) * b_pole(j))
enddo

!

if (has_attribute(ele, 'X1_LIMIT')) then
  if (ele%value(x1_limit$) /= 0) line = trim(line) // ', x1_limit = ' // trim(aper_str(-ele%value(x1_limit$)))
  if (ele%value(x2_limit$) /= 0) line = trim(line) // ', x2_limit = ' // trim(aper_str(ele%value(x2_limit$)))
  if (ele%value(y1_limit$) /= 0) line = trim(line) // ', y1_limit = ' // trim(aper_str(-ele%value(y1_limit$)))
  if (ele%value(y2_limit$) /= 0) line = trim(line) // ', y2_limit = ' // trim(aper_str(ele%value(y2_limit$)))

  if (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0 .or. &
      ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) then
    if (ele%aperture_type == elliptical$) then
      line = trim(line) // ', aperture_shape = ApertureShape.Elliptical'
    else
      line = trim(line) // ', aperture_shape = ApertureShape.Rectangular'
    endif
  endif
endif

!


select case (ele%key)
case (match$, taylor$)
  line = trim(line) // ', transport_map = map_' // trim(ele_name)

! Only the periodic planar model (with kx = 0) can be translated. The field is defined by the
! four-potential written by write_planar_wiggler_four_potential and is integrated by the
! Yoshida integrator using the Bmad step size.

case (wiggler$, undulator$)
  if (is_planar_wiggler(ele)) then

    ! With an integer number of periods the vector potential vanishes at both ends of the element
    ! so the canonical momenta used by SciBmad are the same as the momenta used by Bmad there.
    ! If the number of periods is not an integer, adjust the period so that it is. This is an
    ! approximation and not an error.

    ele2 => pointer_to_field_ele(ele, 1)   ! In case ele is a super_slave. If not, ele2 == ele.
    n_per = ele2%value(l$) / ele2%value(l_period$)
    n_wig = max(1, nint(n_per))
    if (abs(n_per - n_wig) > 1e-8_rp * n_wig .and. warn) then
      call out_io(s_warn$, r_name, ele_full_name(ele2) // ' does not have an integer number of periods.', &
                            '     L_PERIOD will be adjusted to make an integer number of periods.')
    endif

    k_wig = twopi * n_wig / ele2%value(l$)
    n_step = max(1, nint(ele2%value(num_steps$)))
    i_order = nint(ele2%value(integrator_order$))
    phase = k_wig * ((ele%s_start - ele2%s_start) - 0.5_rp * ele2%value(l$))
    if (all(i_order /= [2, 4, 6, 8])) i_order = 4   ! Yoshida only accepts these orders.

    line = trim(line) // ', kind = ' // quote('Wiggler')
    line = trim(line) // ', four_potential = planar_wiggler_four_potential'
    line = trim(line) // ', four_potential_params = (' // re_str(ele%value(b_max$)) // ', ' // &
                                      re_str(k_wig) // ', ' // re_str(phase) // ')'
    line = trim(line) // ', four_potential_normalized = false'
    line = trim(line) // ', tracking_method = Yoshida(order = ' // int_str(i_order) // &
                                                   ', n_steps = ' // int_str(n_step) // ')'
  else
    if (warn) call out_io(s_warn$, r_name, ele_full_name(ele) // ' does not use the periodic planar model. This cannot yet be translated!')
    if (warn) xlate_err = .true.
  endif
end select

!

if (ele%key == patch$) then
  if (ele%value(t_offset$) /= 0)      line = trim(line) // ', dt = ' // re_str(ele%value(t_offset$))
  if (ele%value(x_offset$) /= 0)      line = trim(line) // ', dx = ' // re_str(ele%value(x_offset$))
  if (ele%value(y_offset$) /= 0)      line = trim(line) // ', dy = ' // re_str(ele%value(y_offset$))
  if (ele%value(z_offset$) /= 0)      line = trim(line) // ', dz = ' // re_str(ele%value(z_offset$))
  if (ele%value(y_pitch$) /= 0)       line = trim(line) // ', dx_rot = ' // re_str(-ele%value(y_pitch$))
  if (ele%value(x_pitch$) /= 0)       line = trim(line) // ', dy_rot = ' // re_str(ele%value(x_pitch$))
  if (ele%value(tilt$) /= 0)          line = trim(line) // ', dz_rot = ' // re_str(ele%value(tilt$))
  if (ele%value(E_tot_offset$) /= 0)  line = trim(line) // ', dE_ref = ' // re_str(ele%value(E_tot_offset$))
  if (ele%value(E_tot_set$) /= 0)     line = trim(line) // ', E_ref = ' // re_str(ele%value(E_tot_set$))

else
  if (has_attribute(ele, 'X_PITCH')) then
    if (ele%value(x_offset$) /= 0)  line = trim(line) // ', x_offset = ' // re_str(ele%value(x_offset$))
    if (ele%value(y_offset$) /= 0)  line = trim(line) // ', y_offset = ' // re_str(ele%value(y_offset$))
    if (ele%value(z_offset$) /= 0)  line = trim(line) // ', z_offset = ' // re_str(ele%value(z_offset$))
    if (ele%value(y_pitch$) /= 0)  line = trim(line) // ', x_rot = ' // re_str(-ele%value(y_pitch$))
    if (ele%value(x_pitch$) /= 0)  line = trim(line) // ', y_rot = ' // re_str(ele%value(x_pitch$))
  endif

  if (has_attribute(ele, 'TILT')) then
    if (ele%value(tilt$) /= 0)  line = trim(line) // ', tilt = ' // re_str(ele%value(tilt$))
  endif
endif

!

if (has_attribute(ele, 'KS')) then
  if (ele%field_master) then
    if (ele%value(bs_field$) /= 0)  line = trim(line) // ', bsol_field = ' // re_str(ele%value(bs_field$))
  else
    if (ele%value(ks$) /= 0)  line = trim(line) // ', Ksol = ' // re_str(ele%value(ks$))
  endif
endif

!

if (has_attribute(ele, 'RF_FREQUENCY')) then
  if (is_true(ele%value(harmon_master$))) then
    if (ele%value(harmon$) /= 0)  line = trim(line) // ', harmon = ' // re_str(ele%value(harmon$))
  else
    if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
  endif
endif

if (ele%key == lcavity$) then
  if (ele%value(voltage$)+ele%value(voltage_err$) /= 0)  line = trim(line) // ', voltage = ' // re_str((ele%value(voltage$) + ele%value(voltage_err$)))
  if (ele%value(phi0$) /= 0)  line = trim(line) // ', phi0 = ' // re_str(ele%value(phi0$) + ele%value(phi0_err$))
  ! Note: SaganCavity wants n_cells to be an integer.
  line = trim(line) // ', tracking_method = SaganCavity(n_cells = ' // int_str(nint(ele%value(n_rf_steps$))) // &
                                                    ', L_active = ' // re_str(ele%value(L_active$)) // ')'

elseif (has_attribute(ele, 'RF_FREQUENCY')) then
  if (ele%key == rfcavity$) line = trim(line) // ', zero_phase = PhaseRef.AboveTransition'
  if (ele%value(voltage$) /= 0)  line = trim(line) // ', voltage = ' // re_str(ele%value(voltage$)/abs(charge_of(lat2%branch(ele%ix_branch)%param%particle)))
  if (ele%value(phi0$) /= 0)  line = trim(line) // ', phi0 = ' // re_str(ele%value(phi0$))
endif

if (has_attribute(ele, 'CAVITY_TYPE')) then
  if (nint(ele%value(cavity_type$)) == standing_wave$) then
    line = trim(line) // ', traveling_wave = false'
  else
    line = trim(line) // ', traveling_wave = true'
  endif
endif

!

if (ele%type /= ' ') line = trim(line) // ', label = ' // quote(ele%type)
if (ele%alias /= ' ') line = trim(line) // ', alias = ' // quote(ele%alias)
if (associated(ele%descrip)) line = trim(line) // ', description = ' // quote(ele%descrip)

!

if (ele%key == fork$ .or. ele%key == photon_fork$) then
  n = nint(ele%value(ix_to_branch$))
!!!      line = trim(line) // ', to_line = ' // quote(downcase(lat2%branch(n)%name))
  if (ele%value(ix_to_element$) > 0) then
    i = nint(ele%value(ix_to_element$))
!!!        line = trim(line) // ', to_element = ' // quote(scibmad_ele_name(lat2%branch(n)%ele(i)))
  endif
endif

!

ix = index(line, '(, ')
if (ix == 0) then
  line = trim(line) // ')'
else
  line = line(1:ix) // trim(line(ix+3:)) // ')'
endif

end subroutine ele_def_line

!----------------------------------------------------------------------------------------------
! contains

function aper_str (limit) result (ap_str)

real(rp) limit
character(24) ap_str

!

if (limit == 0) then
  ap_str = 'NaN'
else
  ap_str = re_str(limit)
endif

end function aper_str

!----------------------------------------------------------------------------------------------
! contains

subroutine write_scibmad_element (line, iu, ele)

type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu, ix

!

if (ele%orientation == 1) then
  write (line, '(4a)') trim(line), ' ', trim(scibmad_ele_name(ele%name, ele%ix_branch)), ','
else
  write (line, '(4a)') trim(line), ' reverse(', trim(scibmad_ele_name(ele%name, ele%ix_branch)), '),'
endif

if (len_trim(line) > 100) call write_lat_line(line, iu, .false., ampersand_at_ends = .false.)

end subroutine write_scibmad_element

!--------------------------------------------------------------------------------
! contains

function scibmad_ele_name(name, ix_branch) result (name_out)

character(40) name, name_out
integer, optional :: ix_branch
integer ix, ib

!

ib = integer_option(0, ix_branch) + 1
name_out = downcase(name)
if (name_out == 'end') name_out = 'end_b' // int_str(ib)
if (name_out == 'beginning') name_out = 'begin_b' // int_str(ib)


ix = index(name_out, '#')
if (ix /= 0) name_out = name_out(1:ix-1) // '_s' // name_out(ix+1:)

ix = index(name_out, '\')     !'
if (ix /= 0) name_out = name_out(1:ix-1) // '_m' // name_out(ix+1:)

call str_substitute(name_out, '.', '_')

end function scibmad_ele_name

!--------------------------------------------------------------------------------
! contains

subroutine write_this_taylor(iu, ele, taylor)

type (ele_struct) ele
type (taylor_struct), target :: taylor(6)
type (taylor_term_struct) term

integer iu, i, j, k
integer e_max(6)
character(200) line

!

write (iu, '(a)')
write (iu, '(9a)') 'function map_', trim(scibmad_ele_name(ele%name)), '(v, q)'

e_max = 0
do i = 1, 6
  do j = 1, 6
    e_max(j) = max(e_max(j), maxval(taylor(i)%term(:)%expn(j)))
  enddo
enddo

do i = 0, 3
  if (.not. associated(ele%spin_taylor(i)%term)) cycle
  if (size(ele%spin_taylor(i)%term) == 0) cycle
  do j = 1, 6
    e_max(j) = max(e_max(j), maxval(ele%spin_taylor(i)%term(:)%expn(j)))
  enddo
enddo

!

do i = 1, 6
  write (iu, '(2(a, i0))') '  v_out', i, '= '
  do j = 1, size(taylor(i)%term)
    term = taylor(i)%term(j)
    if (write_lat_debug_flag) then  ! Used for regression tests
      write (line, '(4x, es12.4)') term%coef
    else
      write (line, '(4x, es24.16)') term%coef
    endif

    do k = 1, 6
      if (term%expn(k) == 0) cycle
      if (term%expn(k) == 1) then
        write (line, '(a, 3(a, i0))') trim(line), '*v[', k, ']'
      else
        write (line, '(a, 3(a, i0))') trim(line), '*v[', k, ']^', term%expn(k)
      endif
    enddo
    if (j  == size(taylor(i)%term)) then
      write (iu, '(a)') line 
    else
      write (iu, '(a)') trim(line) // ' +' 
    endif
  enddo
enddo

!

write (iu, '(a)') 

do i = 0, 3
  if (.not. associated(ele%spin_taylor(i)%term)) then
    write (iu, '(a, i0, 2a)') ' q_out', i, ' = ', unit_spin_map(i)
    cycle
  elseif (size(ele%spin_taylor(i)%term) == 0) then
    write (iu, '(a, i0, 2a)') '  q_out', i, ' = ', unit_spin_map(i)
    cycle
  endif

  write (iu, '(2(a, i0))') '  q_out', i, ' = '
  do j = 1, size(ele%spin_taylor(i)%term)
    term = ele%spin_taylor(i)%term(j)
    if (write_lat_debug_flag) then  ! Used for regression tests
      write (line, '(4x, es13.5)') term%coef
    else
      write (line, '(4x, es24.16)') term%coef
    endif

    do k = 1, 6
      if (term%expn(k) == 0) cycle
      if (term%expn(k) == 1) then 
        write (line, '(a, 3(a, i0))') trim(line), '*q[', k, ']'
      else
        write (line, '(a, 3(a, i0))') trim(line), '*q[', k, ']^', term%expn(k)
      endif
    enddo
    if (j  == size(ele%spin_taylor(i)%term)) then
      write (iu, '(a)') line 
    else
      write (iu, '(a)') trim(line) // ' +' 
    endif
  enddo
enddo

write (iu, '(a)') 
write (iu, '(a)') '  return (v_out1, v_out2, v_out3, v_out4, v_out5, v_out6), (q_out1, q_out2, q_out3, q_out4)'
write (iu, '(a)') 'end'

end subroutine write_this_taylor

!------------------------------------------------------
! contains

recursive function expression_kernel(stack, lord, use_old_names, e_ptr) result(expr_str)

type (expression_atom_struct), target :: stack(:)
type (expression_atom_struct) :: stack2(size(stack))
type (ele_struct) lord
type (this_expr_struct), optional :: e_ptr

integer ix_match, i, n
character(1000) expr_str
character(100) var_name
logical use_old_names, is_new

!

stack2 = stack

do i = 1, size(stack2)
  select case (downcase(stack2(i)%name))
  case ('c_light', 'm_electron', 'm_proton', 'm_neutron', 'm_muon', 'm_pion_0', 'm_pion_charged', &
        'm_deuteron', 'm_helion', 'h_planck')
    stack2(i)%name = upcase(stack2(i)%name)
  case ('pi', 'sqrt', 'log', 'exp', 'sin', 'cos', 'tan', 'cot', 'asin', 'acos', 'atan', 'sinh', 'cosh', &
        'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'acoth', 'abs', 'factorial', 'sign')
    stack2(i)%name = downcase(stack2(i)%name)
  case ('twopi')
    stack2(i)%name = '2*pi'
  case ('fourpi')
    stack2(i)%name = '4*pi'
  case ('e', 'e_log')
    stack2(i)%name = 'exp(1.0)'
  case ('sqrt_2')
    stack2(i)%name = 'sqrt(2.0)'
  case ('degrad')
    stack2(i)%name = '(180 / pi)'
  case ('degrees', 'raddeg')
    stack2(i)%name = '(pi / 180)'
  case ('r_e')
    stack2(i)%name = 'R_ELECTRON'
  case ('r_p')
    stack2(i)%name = 'R_PROTON'
  case ('h_bar_planck')
    stack2(i)%name = 'H_BAR'
  case ('e_charge')
    stack2(i)%name = 'E_CHARGE'
  case ('fine_struct_const')
    stack2(i)%name = 'FINE_STRUCTURE'
  case ('emass')
    stack2(i)%name = '(1e-9 * M_ELECTRON)'
  case ('pmass')
    stack2(i)%name = '(1e-9 * M_PROTON)'
  case ('anom_moment_electron')
    stack2(i)%name = 'ANOMALY_ELECTRON'
  case ('anom_moment_muon')
    stack2(i)%name = 'ANOMALY_MUON'
  case ('anom_moment_proton')
    stack2(i)%name = 'gyromagnetic_anomaly(Species("proton"))'
  case ('anom_moment_deuteron')
    stack2(i)%name = 'gyromagnetic_anomaly(Species("deuteron"))'

  case ('atan2')
    stack2(i)%name = 'atan'
  case ('modulo')
    stack2(i)%name = 'mod'
  case ('sinc')
    stack2(i)%name = 'sincu'
  case ('ran')
    stack2(i)%name = 'rand'
  case ('ran_gauss')
    stack2(i)%name = 'randn'
  case ('int')
    stack2(i)%name = 'trunc'
  case ('nint')
    stack2(i)%name = 'round'
  case ('floor')
    stack2(i)%name = 'floor'
  case ('ceiling')
    stack2(i)%name = 'ceil'
  case ('mass_of')
    stack2(i)%name = 'massof'
  case ('charge_of')
    stack2(i)%name = 'chargeof'
  case ('anomalous_modment_of')
    stack2(i)%name = ''
  case ('species')
    stack2(i)%name = 'Species'

  case default
    select case (stack2(i)%type)
    case (constant$)      ! Something like "c_light"
      stack2(i)%name = upcase(stack2(i)%name)

    case (variable$)
      stack2(i)%name = 'c.' // downcase(stack2(i)%name)

    case default
      if (stack2(i)%type > var_offset$ .and. stack2(i)%type < var_offset$ + n_var_max$) then
        if (use_old_names) then
          ! The "old" value of a group variable is its present value.
          var_name = trim(scibmad_ele_name(lord%name)) // '_' // trim(downcase(stack2(i)%name))
          is_new = .true.
          n = 0
          if (allocated(e_ptr%group_var_names)) then
            n = size(e_ptr%group_var_names)
            is_new = (.not. any(e_ptr%group_var_names == var_name))
          endif
          if (is_new) then
            call re_allocate(e_ptr%group_var_names, n+1)
            call re_allocate(e_ptr%group_var_values, n+1)
            e_ptr%group_var_names(n+1) = var_name
            e_ptr%group_var_values(n+1) = lord%control%var(stack2(i)%type - var_offset$)%value
          endif
          stack2(i)%name = 'c.old_' // trim(var_name) // '__' // trim(e_ptr%scibmad_ele) // '_' // trim(e_ptr%scibmad_attrib)

        else
          stack2(i)%name = 'c.' // trim(scibmad_ele_name(lord%name)) // '_' // downcase(stack2(i)%name)
        endif
      endif
    end select
  end select
enddo

expr_str = expression_stack_to_string(stack2)

end function expression_kernel

!------------------------------------------------------
! contains

! Return SciBmad attribute name(s) given Bmad attribute name.
! Since a Bmad attribute may map to more than one SciBmad attribute (EG: HKICK of a tilted element
! maps to both the normal and skew n = 0 multipole components), up to two names are returned.
! n_sci = 0 => Attribute cannot be translated.
! The value of the SciBmad attribute sci_attrib_name(i) is factor(i) times the value of the Bmad attribute.

subroutine scibmad_attrib_name(bmad_attrib_name, ele, n_sci, sci_attrib_name, factor)

type (ele_struct) ele

integer n_sci
character(*) bmad_attrib_name
character(*) sci_attrib_name(2)
real(rp) factor(2)

!

n_sci = 1
sci_attrib_name = ''
factor = 1.0_rp

if ((bmad_attrib_name(1:1) == 'A' .or. bmad_attrib_name(1:1) == 'B') .and. is_integer(bmad_attrib_name(2:), ix)) then
  if (ele%field_master) then
    sci_attrib_name(1) = 'B'
    factor(1) = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    sci_attrib_name(1) = 'K'
    factor(1) = 1
  endif

  if (bmad_attrib_name(1:1) == 'A') then
    sci_attrib_name(1) = sci_attrib_name(1)(1:1) // 's'
  else
    sci_attrib_name(1) = sci_attrib_name(1)(1:1) // 'n'
  endif

  sci_attrib_name(1) = trim(sci_attrib_name(1)) // bmad_attrib_name(2:)
  factor(1) = factor(1) * factorial(ix)

  if (ele%value(l$) == 0) then
    sci_attrib_name(1) = trim(sci_attrib_name(1)) // 'L'
  else
    factor(1) = factor(1) / ele%value(l$)
  endif

  return
endif

! Kick attributes translate to n = 0 multipole components.

select case (bmad_attrib_name)
case ('KICK', 'HKICK', 'VKICK', 'BL_KICK', 'BL_HKICK', 'BL_VKICK')
  call scibmad_kick_attrib_name(bmad_attrib_name, ele, n_sci, sci_attrib_name, factor)
  return
end select

select case (bmad_attrib_name)
case ('B1_GRADIENT');   sci_attrib_name(1) = 'Bn1'
case ('B2_GRADIENT');   sci_attrib_name(1) = 'Bn2'
case ('B3_GRADIENT');   sci_attrib_name(1) = 'Bn3'
case ('K1');            sci_attrib_name(1) = 'Kn1'
case ('K2');            sci_attrib_name(1) = 'Kn2'
case ('K3');            sci_attrib_name(1) = 'Kn3'
case ('E1');            sci_attrib_name(1) = 'e1'
case ('E2');            sci_attrib_name(1) = 'e2'
case ('G');             sci_attrib_name(1) = 'g_ref'
case ('ANGLE');         sci_attrib_name(1) = 'g_ref'
case ('L');             sci_attrib_name(1) = 'L'
case ('X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH', 'TILT')
  if (ele%key == patch$) then
    select case (bmad_attrib_name)
    case ('X_OFFSET');      sci_attrib_name(1) = 'dx'
    case ('Y_OFFSET');      sci_attrib_name(1) = 'dy'
    case ('Z_OFFSET');      sci_attrib_name(1) = 'dz'
    case ('X_PITCH');       sci_attrib_name(1) = 'dy_rot'
    case ('Y_PITCH');       sci_attrib_name(1) = 'dx_rot'; factor(1) = -1
    case ('TILT');          sci_attrib_name(1) = 'dz_rot'
    end select
  else
    select case (bmad_attrib_name)
    case ('X_OFFSET');      sci_attrib_name(1) = 'x_offset'
    case ('Y_OFFSET');      sci_attrib_name(1) = 'y_offset'
    case ('Z_OFFSET');      sci_attrib_name(1) = 'z_offset'
    case ('X_PITCH');       sci_attrib_name(1) = 'y_rot'
    case ('Y_PITCH');       sci_attrib_name(1) = 'x_rot'; factor(1) = -1
    case ('TILT');          sci_attrib_name(1) = 'z_rot'
    end select
  endif

case ('T_OFFSET');      sci_attrib_name(1) = 't_offset'
case ('KS');            sci_attrib_name(1) = 'Ksol'
case ('BS_FIELD');      sci_attrib_name(1) = 'Bsol'

! These group specific attributes vary the lengths of neighboring elements. There is no
! SciBmad equivalent.

case ('START_EDGE', 'END_EDGE', 'ACCORDION_EDGE', 'S_POSITION', 'LORD_PAD1', 'LORD_PAD2')
  n_sci = 0
  xlate_err = .true.
  call out_io(s_warn$, r_name, 'Group control of the ' // trim(bmad_attrib_name) // ' attribute of element ' // &
                                                trim(ele%name) // ' not yet coded for translation.')

case default
  if (ele%key == group$ .or. ele%key == overlay$) then
    sci_attrib_name(1) = bmad_attrib_name
    return   ! No problem translating controller vars.
  endif
  n_sci = 0
  xlate_err = .true.
  call out_io(s_warn$, r_name, 'Attribute not yet coded for translation: ' // trim(bmad_attrib_name), &
                               'Please report this.')
end select

end subroutine scibmad_attrib_name

!------------------------------------------------------
! contains

! Return SciBmad attribute name(s) for the Bmad kick attributes KICK, HKICK, VKICK and the
! corresponding integrated field attributes BL_KICK, BL_HKICK, BL_VKICK.
! In SciBmad a kick is represented by the n = 0 multipole components so a kick attribute of a
! tilted element must be distributed between the normal and skew components.
! The conversion mirrors what multipole_ele_to_ab does with the kick attributes (which is what is
! used when writing the element definitions) so that overlay/group controlled values are consistent
! with the element definition values.

subroutine scibmad_kick_attrib_name(bmad_attrib_name, ele, n_sci, sci_attrib_name, factor)

type (ele_struct) ele

integer n_sci, key, i
real(rp) factor(2), f0, tilt, coef(2)
character(*) bmad_attrib_name
character(*) sci_attrib_name(2)
character(1) prefix
logical is_hkick

! is_hkick = True if the attribute gives a kick in the horizontal plane (in the element body frame).

key = ele%key
is_hkick = (index(bmad_attrib_name, 'VKICK') == 0)
if (key == vkicker$) is_hkick = .false.
if (key == hkicker$) is_hkick = .true.

! coef(1) is the normal (Kn0/Bn0) coefficient and coef(2) is the skew (Ks0/Bs0) coefficient.
! Note: For kicker type elements the kick is defined in the element body frame so there is no
! rotation by the element tilt.

select case (key)
case (hkicker$, vkicker$, kicker$, ac_kicker$)
  if (is_hkick) then
    coef = [-1.0_rp, 0.0_rp]
  else
    coef = [0.0_rp, 1.0_rp]
  endif

case (elseparator$)   ! Kick is electric
  if (ele%value(l$) == 0) then
    n_sci = 0
    return
  endif

  if (is_hkick) then
    sci_attrib_name(1) = 'En0'
    factor(1) = -ele%value(p0c$) / ele%value(l$)
  else
    sci_attrib_name(1) = 'Es0'
    factor(1) = ele%value(p0c$) / ele%value(l$)
  endif
  n_sci = 1
  return

case default
  if (key == sbend$ .or. key == rf_bend$) then
    tilt = ele%value(ref_tilt_tot$)
  else
    tilt = ele%value(tilt_tot$)
  endif

  if (is_hkick) then
    coef = [-cos(tilt), -sin(tilt)]
  else
    coef = [-sin(tilt), cos(tilt)]
  endif
end select

! BL_KICK, BL_HKICK and BL_VKICK are integrated field values so no scaling by the reference momentum.

if (bmad_attrib_name(1:3) == 'BL_') then
  prefix = 'B'
  f0 = 1
elseif (ele%field_master) then
  prefix = 'B'
  f0 = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
else
  prefix = 'K'
  f0 = 1
endif

if (ele%value(l$) /= 0) f0 = f0 / ele%value(l$)

!

n_sci = 0

if (coef(1) /= 0) then
  n_sci = n_sci + 1
  sci_attrib_name(n_sci) = prefix // 'n0'
  factor(n_sci) = f0 * coef(1)
endif

if (coef(2) /= 0) then
  n_sci = n_sci + 1
  sci_attrib_name(n_sci) = prefix // 's0'
  factor(n_sci) = f0 * coef(2)
endif

if (ele%value(l$) == 0) then
  do i = 1, n_sci
    sci_attrib_name(i) = trim(sci_attrib_name(i)) // 'L'
  enddo
endif

end subroutine scibmad_kick_attrib_name

!------------------------------------------------------
! contains

! Return the value of the SciBmad multipole attribute sci_attrib_name as computed when writing the
! element definition. is_multipole is set False if sci_attrib_name is not a multipole attribute.
! This is needed since a given SciBmad multipole component may get contributions from several
! Bmad attributes (EG: Kn0 of a bend gets contributions from HKICK, VKICK, DG and the bend angle)
! and the part not controlled by an overlay must be added in when writing a controlled value.

function scibmad_multipole_value(ele, sci_attrib_name, is_multipole) result (value)

type (ele_struct) ele

real(rp) value, ff, a_p(0:n_pole_maxx), b_p(0:n_pole_maxx)
integer nlen, nord, ixp
character(*) sci_attrib_name
character(40) nam
logical is_multipole

!

value = 0
is_multipole = .false.

nam = sci_attrib_name
if (nam(2:2) /= 'n' .and. nam(2:2) /= 's') return

! Electric multipoles are not scaled by the element length.

if (nam(1:1) == 'E') then
  if (.not. is_integer(nam(3:), nord)) return
  if (nord > n_pole_maxx) return
  call multipole_ele_to_ab(ele, .false., ixp, a_p, b_p, electric$, include_kicks$)
  ff = 1

else
  if ((nam(1:1) == 'B') .neqv. ele%field_master) return   ! Prefix is 'B' if and only if field_master.

  nlen = len_trim(nam)
  if (nam(nlen:nlen) == 'L') then
    if (ele%value(l$) /= 0) return
    nam = nam(1:nlen-1)
  else
    if (ele%value(l$) == 0) return
  endif

  if (.not. is_integer(nam(3:), nord)) return
  if (nord > n_pole_maxx) return

  call multipole_ele_to_ab(ele, .false., ixp, a_p, b_p, magnetic$, include_kicks$)
  if (ele%key == sbend$) b_p(0) = b_p(0) + ele%value(angle$)

  if (ele%field_master) then
    ff = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    ff = 1
  endif

  if (ele%value(l$) /= 0) ff = ff / ele%value(l$)
endif

!

if (nam(2:2) == 's') then
  value = ff * factorial(nord) * a_p(nord)
else
  value = ff * factorial(nord) * b_p(nord)
endif

is_multipole = .true.

end function scibmad_multipole_value

!------------------------------------------------------
! contains

function def_expr(e_ptr, add_def_prefix) result (expr_str)

type (this_expr_struct) e_ptr
character(5000) expr_str
logical add_def_prefix

! First three characters of %expr are " + " which can be dropped

expr_str = trim(e_ptr%expr(4:))
if (add_def_prefix) expr_str = 'DefExpr(c -> ' // trim(expr_str) // ')'

end function def_expr

!------------------------------------------------------
! contains
!
! For each set of tracking elements that share a name, if the SciBmad definitions of the elements
! are not all the same, append the suffix to the names of all the elements in the set.
! The "?" in the suffix is replaced by an index which numbers the elements in lattice order.
! Since the comparison uses the SciBmad definition line, any difference that would show up in the
! translation (parameter values, multipoles, is_on, aperture type, etc.) is detected.
! Note: The beginning and end elements are not renamed since scibmad_ele_name makes their names unique.

subroutine this_create_unique_ele_names (lat, suffix)

type (lat_struct), target :: lat
type (nametable_struct), pointer :: ntab
type (ele_struct), pointer :: ele0, ele

integer i_nt, i_end, i2, ix_p, nn
integer, allocatable :: indx(:)

logical all_same

character(*) suffix
character(40) suff, name0
character(4000) line0, line2

! Find '?' character

ix_p = index(suffix, '?')
if (ix_p == 0) then
  call out_io (s_error$, r_name, 'SUFFIX DOES NOT HAVE A "?" CHARACTER: ' // suffix)
  return
endif

suff = suffix
call str_upcase (suff, suff)

!

ntab => lat%nametable
allocate(indx(ntab%n_max - ntab%n_min + 1))
i_nt = ntab%n_min

do
  if (i_nt > ntab%n_max) exit

  ! Find the range [i_nt, i_end] of nametable entries that share the same name.

  name0 = ntab%name(ntab%index(i_nt))
  i_end = i_nt
  do
    if (i_end == ntab%n_max) exit
    if (ntab%name(ntab%index(i_end+1)) /= name0) exit
    i_end = i_end + 1
  enddo

  ! Collect the tracking elements with this name. Lord elements are not written as elements
  ! so are ignored.

  nn = 0
  if (i_end > i_nt .and. name0 /= 'BEGINNING' .and. name0 /= 'END') then
    do i2 = i_nt, i_end
      ele => pointer_to_ele(lat, ntab%index(i2))
      if (ele%ix_ele > lat%branch(ele%ix_branch)%n_ele_track) cycle
      nn = nn + 1
      indx(nn) = ntab%index(i2)
    enddo
  endif

  i_nt = i_end + 1
  if (nn < 2) cycle

  ! Compare definitions. The name is left out of the definition line since that is what is being decided.

  ele0 => pointer_to_ele(lat, indx(1))
  call ele_def_line(ele0, '', line0, .false.)
  all_same = .true.
  do i2 = 2, nn
    ele => pointer_to_ele(lat, indx(i2))
    call ele_def_line(ele, '', line2, .false.)
    if (line2 /= line0 .or. .not. same_transport_map(ele0, ele)) then
      all_same = .false.
      exit
    endif
  enddo

  if (all_same) cycle

  ! The nametable index array does not have any sort order with respect to the order in the lattice
  ! So do a sort

  call super_sort(indx(1:nn))

  do i2 = 1, nn
    ele => pointer_to_ele(lat, indx(i2))
    ele%name = trim(ele%name) // suff(1:ix_p-1) // int_str(i2) // suff(ix_p+1:)
  enddo
enddo

end subroutine this_create_unique_ele_names

!------------------------------------------------------
! contains
!
! Match and Taylor elements reference a transport map function which is not part of the definition line.
! Return True if the maps of ele1 and ele2 are the same (or if neither element has a map).

function same_transport_map(ele1, ele2) result (is_same)

type (ele_struct) ele1, ele2
logical is_same
integer i

!

is_same = .true.

select case (ele1%key)
case (match$)
  is_same = (all(ele1%vec0 == ele2%vec0) .and. all(ele1%mat6 == ele2%mat6))

case (taylor$)
  do i = 1, 6
    if (.not. same_taylor(ele1%taylor(i)%term, ele2%taylor(i)%term)) is_same = .false.
  enddo
  do i = 0, 3
    if (.not. same_taylor(ele1%spin_taylor(i)%term, ele2%spin_taylor(i)%term)) is_same = .false.
  enddo
end select

end function same_transport_map

!------------------------------------------------------
! contains

function same_taylor(term1, term2) result (is_same)

type (taylor_term_struct), pointer :: term1(:), term2(:)
logical is_same
integer j

!

is_same = .false.
if (associated(term1) .neqv. associated(term2)) return

if (associated(term1)) then
  if (size(term1) /= size(term2)) return
  do j = 1, size(term1)
    if (term1(j)%coef /= term2(j)%coef) return
    if (any(term1(j)%expn /= term2(j)%expn)) return
  enddo
endif

is_same = .true.

end function same_taylor

end subroutine write_lattice_scibmad_format
