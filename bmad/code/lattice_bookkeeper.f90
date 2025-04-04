!+
! Subroutine lattice_bookkeeper (lat, err_flag)
!
! Subroutine to do basic bookkeeping for a lattice:
!   lord/slave control
!   reference energy calc
!   s-position calc
!   geometry (floor position) calc
!
! Not done are any higher level calculations. Twiss, transfer matrices, 
! orbit calculations, anything involving tracking, etc.
!
! Note: This this routine does a complete job of bookking
! and could be unacceptably slow if lat%auto_bookkeeper = True.
!
! Note: The documentation for a routine should say if a call to lattice_bookkeeper is needed
! afterwards. If it is not mentioned in the documentation for a routine, a call to 
! lattice_bookkeeper is not needed afterwards.
!
! Input:
!   lat   -- lat_struct: Lattice needing bookkeeping.
!
! Output:
!   lat      -- lat_struct: Lattice with bookkeeping done.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine lattice_bookkeeper (lat, err_flag)

use bookkeeper_mod, dummy => lattice_bookkeeper
use precision_constants, only: e_muon  ! PTC
use bmad_parser_struct, only: bp_com

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (bookkeeping_state_struct), pointer :: stat

real(rp) dval(num_ele_attrib$), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
integer i, j, ix

logical, optional :: err_flag
logical found, err, auto_saved

character(20), parameter :: r_name = 'lattice_bookkeeper'

! Set PTC E_MUON just to make sure it has the same value as bmad_com%electric_dipole_moment

E_MUON = bmad_com%electric_dipole_moment

! Turn on intelligent bookkeeping while this routine is running
! If bookkeeping has not been intelligent then mark everything as stale.

auto_saved = bmad_com%auto_bookkeeper
if (auto_saved) then
  bmad_com%auto_bookkeeper = .false.
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 0, branch%n_ele_max
      call set_ele_status_stale (branch%ele(j), all_groups$, .false.) 
      call attributes_need_bookkeeping(branch%ele(j), dval)

      if (any(dval /= 0) .and. bp_com%parser_name == '') then  ! If not parsing then error
        call out_io (s_warn$, r_name, &
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
          '!!!!! Using intelligent bookkeeping is now mandated for all Bmad based programs.               !!!!!', &
          '!!!!! See the "Intelligent Bookkeeping" section in the Bmad manual.                            !!!!!', &
          '!!!!! This program will run now but if this program modifies any lattice parameters, the       !!!!!', &
          '!!!!! correctness of the results is questionable.                                              !!!!!', &
          '!!!!! Contact the maintainer of this program with this information.                            !!!!!', &
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      endif
    enddo
  enddo
endif

if (present(err_flag)) err_flag = .true.

! Cases where bookkeeping routines have to be called multiple times:
!   * Control_bookkeeper to make sure that girders are correctly adjusted.
!       Especially when lattice has bends with field_master = True.
!   * lat_geometry in case the energy changes and there is a bend with field_master = T.
!   * Recompute ref energy for cases where a flexible patch has changed its geometry and this
!       affects the reference energy due to the presence of lcavity elements.

call ramper_slave_setup(lat)

do i = 1, 3
  call control_bookkeeper (lat, err_flag = err);     if (err) return
  call lat_geometry (lat)
  call s_calc (lat)
  call lat_compute_ref_energy_and_time (lat, err);   if (err) return
enddo

call ptc_bookkeeper (lat)

! Make sure the multipole cache is initialized

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  do j = 0, branch%n_ele_max
    ele => branch%ele(j)
    if (.not. allocated(ele%multipole_cache)) cycle
    if (.not. ele%multipole_cache%mag_valid) call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole)
    if (.not. ele%multipole_cache%elec_valid) call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, electric$)
  enddo
enddo

! See if all status flags have been properly reset.
! Exception is mat6 flag since the bookkeeping routines do not touch this.

stat => lat%lord_state
if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
    stat%s_position == stale$ .or. stat%ref_energy == stale$) then
  call out_io (s_info$, r_name, 'Stale bookkeeping lord_status flags detected.', &
                                'Please contact DCS!', 'Status: \5i6\ ', &
          i_array = [stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
endif
call reset_status_flags_to_ok(stat)

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  stat => branch%param%bookkeeping_state
  if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
      stat%s_position == stale$ .or. stat%ref_energy == stale$) then
    call out_io (s_info$, r_name, 'Stale bookkeeping status flags detected at branch: \i0\.', &
                                  'Please contact DCS!', 'Status: \5i6\ ', &
            i_array = [i, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
  endif
  call reset_status_flags_to_ok(stat)

  do j = 0, branch%n_ele_max
    ele => branch%ele(j)
    if (ele%key == null_ele$) cycle 
    stat => ele%bookkeeping_state
    if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
        stat%s_position == stale$ .or. stat%ref_energy == stale$) then
      call out_io (s_info$, r_name, &
            'Stale bookkeeping status flags detected at element: ' // trim(ele%name) // ' (\i0\>>\i0\).', &
            'Please contact DCS!', 'Status: \5i6\ ', &
            i_array = [i, j, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
    endif
    call reset_status_flags_to_ok(stat)
  enddo

enddo

!

if (present(err_flag)) err_flag = .false.
bmad_com%auto_bookkeeper = auto_saved

!----------------------------------------------------------
contains

subroutine reset_status_flags_to_ok (stat)
  type (bookkeeping_state_struct) stat

  if (stat%control /= ok$        .and. stat%control /= super_ok$)        stat%control = ok$
  if (stat%attributes /= ok$     .and. stat%attributes /= super_ok$)     stat%attributes = ok$
  if (stat%floor_position /= ok$ .and. stat%floor_position /= super_ok$) stat%floor_position = ok$
  if (stat%s_position /= ok$     .and. stat%s_position /= super_ok$)     stat%s_position = ok$
  if (stat%ref_energy /= ok$     .and. stat%ref_energy /= super_ok$)     stat%ref_energy = ok$

end subroutine reset_status_flags_to_ok

end subroutine lattice_bookkeeper

