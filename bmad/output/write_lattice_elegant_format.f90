!+ 
! Subroutine write_lattice_elegant_format (out_file_name, lat, ref_orbit, &
!        use_matrix_model, include_apertures, dr12_drift_max, ix_branch, err)
!
! Subroutine to write an Elegant lattice file using the 
! information in a lat_struct. Optionally, only part of the lattice can be generated.
!
! To write a Bmad lattice file, use: write_bmad_lattice_file
!
! Note: sol_quad elements are replaced by a drift-matrix-drift model.
!
! Note: A wiggler or undulator is translated to an Elegant CWIGGLER element along with an SDDS
! file, or a pair of files, holding the harmonic expansion of the field. The Cartesian harmonic
! expansion used by CWIGGLER is the same one used by a Bmad cartesian_map so the translation is
! exact. See the cwiggler_translate routine below for the details and for the list of Bmad field
! maps that cannot be represented this way. A wiggler that cannot be translated to a CWIGGLER is
! replaced by a drift-matrix-drift model.
!
! Input:
!   out_file_name     -- character(*): Name of the mad output lattice file.
!   lat               -- lat_struct: Holds the lattice information.
!   ref_orbit(0:)     -- coord_struct, allocatable, optional: Referece orbit for sad_mult and patch elements.
!                          This argument must be present if the lattice has sad_mult or patch elements and is
!                          being translated to MAD-8 or SAD.
!   use_matrix_model  -- logical, optional: Use a drift-matrix_drift model for wigglers/undulators?
!                           [A MAD "matrix" is a 2nd order Taylor map.] This switch is ignored for SAD conversion.
!                           Default is False -> Use a CWIGGLER element.
!                           Note: sol_quad elements always use a drift-matrix-drift model.
!   include_apertures -- logical, optional: If True (the default), add to the output lattice a zero length
!                           collimator element next to any non-collimator element that has an aperture.
!                           Note: MADX translations for non-drift elements can handle non-collimator elements 
!                           with an aperture so in this case this argument is ignored.
!   dr12_drift_max    -- real(rp), optional: Max deviation for drifts allowed before a correction matrix element
!                           is added. Default value is 1d-5. A negative number means use default.
!   ix_branch         -- Integer, optional: Index of lattice branch to use. Default = 0.
!
! Output:
!   err               -- logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_lattice_elegant_format (out_file_name, lat, ref_orbit, &
                    use_matrix_model, include_apertures, dr12_drift_max, ix_branch, err)

use taylor_mod, dummy => write_lattice_elegant_format
use write_lattice_file_mod, dummy3 => write_lattice_elegant_format
use ptc_interface_mod, only: taylor_inverse, concat_taylor

implicit none

type (lat_struct), target :: lat, lat_out
type (ele_struct), pointer :: ele, ele1, ele2, lord, sol_ele, first_sol_edge
type (ele_struct) :: drift_ele, ab_ele, taylor_ele, col_ele, kicker_ele, null_ele, bend_ele, quad_ele
type (ele_struct) :: wig_ele
type (cartesian_map_struct), target :: wig_cmap    ! Scratch map for periodic type wigglers.
type (coord_struct) orb_start, orb_end, orb_center
type (coord_struct), allocatable, optional :: ref_orbit(:)
type (coord_struct), allocatable :: orbit_out(:)
type (taylor_term_struct) :: term
type (branch_struct), pointer :: branch, branch_out 
type (taylor_struct) taylor_a(6), taylor_b(6)
type (taylor_struct), pointer :: taylor_ptr(:)
type (all_pointer_struct) a_ptr

real(rp), optional :: dr12_drift_max
real(rp) field, hk, vk, limit(2), length, a, b, f, e2, beta, r_max, r0, dr12_max
real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) tilt, x_pitch, y_pitch, etilt, epitch, eyaw, offset(3), w_mat(3,3)

integer, optional :: ix_branch
integer, allocatable :: n_repeat(:), an_indexx(:)
integer i, j, ib, j2, k, n, ix, i_unique, i_line, iout, iu, n_names, j_count, f_count, ix_ele
integer ie, ios, a_count, ix_lord, ix_match, iv, ifa, ix_pole_max
integer ix1, ix2, n_lord, aperture_at, n_name_change_warn, n_elsep_warn, n_taylor_order_saved
integer :: ix_line_min, ix_line_max, n_warn_max, n_wig_model_err, print_wig_model_err_max

character(*), parameter :: r_name = "write_lattice_elegant_format"
character(*) out_file_name
character(300) line, knl_str, ksl_str
character(200) sdds_base
character(40) orig_name, str, bmad_params(20), elegant_params(20)
character(40), allocatable :: names(:)
character(4000) line_out   ! Can be this large for taylor maps.
character(2) continue_char, eol_char, comment_char, separator_char
character(1), parameter :: num(9) = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

logical, optional :: use_matrix_model, include_apertures, err
logical init_needed, err_flag, monopole
logical parsing, warn_printed, converted, ptc_exact_model
logical print_err, wig_note_printed, wig_ok

! Use ptc exact_model = True since this is needed to get the drift nonlinear terms

ptc_exact_model = ptc_com%exact_model
ptc_com%exact_model = .true.
dr12_max = real_option(1d-5, dr12_drift_max)
if (dr12_max < 0) dr12_max = 1d-5

! Init

n_warn_max = 10
n_wig_model_err = 0
print_wig_model_err_max = 5
wig_note_printed = .false.

ix = integer_option(0, ix_branch)
if (ix < 0 .or. ix > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: /i0/ ', i_array = [ix])
  return
endif

branch => lat%branch(ix)

comment_char = '!'
continue_char = ' &'
eol_char = ''
separator_char = ','
ix_line_max = 80

call out_io (s_info$, r_name, '! NOTE: ELEGANT TRANSLATION IN DEVELOPMENT. PLEASE USE WITH CAUTION!.')

ix_line_min = ix_line_max - 20

call init_ele (col_ele)
call init_ele (drift_ele, drift$)
call init_ele (taylor_ele, taylor$)
call init_ele (ab_ele, ab_multipole$)
call init_ele (kicker_ele, kicker$) 
call init_ele (quad_ele, quadrupole$)
call init_ele (bend_ele, sbend$)
call multipole_init (ab_ele, magnetic$)
null_ele%key = null_ele$

allocate (names(branch%n_ele_max+10), an_indexx(branch%n_ele_max+10)) ! list of element names

call out_io (s_info$, r_name, &
      'Note: In general, Bmad lattice elements can have attributes that cannot be translated. ', &
      '      For example, higher order terms in a Taylor element.', &
      '      Please use caution when using a translated lattice.')


! open file

if (present(err)) err = .true.
n_taylor_order_saved = ptc_private%taylor_order_ptc

iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

! Any CWIGGLER harmonic files are named after the lattice file so that they land in the same
! place and so that a relative file name in the lattice resolves the same way the lattice does.

sdds_base = out_file_name
n = len_trim(sdds_base)
if (n > 4) then
  if (sdds_base(n-3:n) == '.lte') sdds_base = sdds_base(1:n-4)
endif

!-----------------------------------------------------------------------------
! Translation is a two step process:
!   1) Create a new lattice called lat_out making substitutions for sol_quad and wiggler elements, etc..
!   2) Use lat_out to create the lattice file.

lat_out = lat
call allocate_lat_ele_array(lat_out, 2*branch%n_ele_max, branch%ix_branch)
branch_out => lat_out%branch(branch%ix_branch)

if (present(ref_orbit)) then
  call reallocate_coord(orbit_out, size(ref_orbit))
  orbit_out = ref_orbit
else
  call reallocate_coord(orbit_out, branch%n_ele_max)
endif

f_count = 0    ! fringe around bends and quads. Also drift nonlinearities.
j_count = 0    ! drift around solenoid or sol_quad index. Also z shift count.
a_count = 0    ! Aperture count
i_unique = 1000

! Loop over all input elements

nullify(first_sol_edge)
n_name_change_warn = 0
n_elsep_warn = 0
ix_ele = 0

do
  ix_ele = ix_ele + 1
  if (ix_ele > branch_out%n_ele_track) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle  ! Has been marked for delection

  val => ele%value

  ! If there is an aperture with an element that is not an ecoll or rcoll then need to make a separate

  if ((val(x1_limit$) /= 0 .or. val(x2_limit$) /= 0 .or. val(y1_limit$) /= 0 .or. val(y2_limit$) /= 0) .and. &
      ele%key /= ecollimator$ .and. ele%key /= rcollimator$ .and. logic_option(.true., include_apertures) .and. &
      (ele%key == drift$)) then

    if (val(x1_limit$) /= val(x2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric x_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(x1_limit$) = max(val(x1_limit$), val(x2_limit$))
    endif

    if (val(y1_limit$) /= val(y2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric y_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(y1_limit$) = max(val(y1_limit$), val(y2_limit$))
    endif

    ! create ecoll and rcoll elements.

    if (ele%aperture_type == rectangular$) then
      col_ele%key = rcollimator$
    else
      col_ele%key = ecollimator$
    endif
    a_count = a_count + 1
    write (col_ele%name, '(a, i0)')  'COLLIMATOR_N', a_count
    col_ele%value = val
    col_ele%value(l$) = 0
    val(x1_limit$) = 0; val(x2_limit$) = 0; val(y1_limit$) = 0; val(y2_limit$) = 0; 
    aperture_at = ele%aperture_at  ! Save since ele pointer will be invalid after the insert
    if (aperture_at == both_ends$ .or. aperture_at == downstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
    endif
    if (aperture_at == both_ends$ .or. aperture_at == upstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele, branch_out%ix_branch, orbit_out)
    endif
    ix_ele = ix_ele - 1 ! Want to process the element again on the next loop.

    cycle ! cycle since ele pointer is invalid

  endif

  ! If the bend has a roll then put kicker elements just before and just after

  if (ele%key == sbend$ .and. val(roll$) /= 0) then
    j_count = j_count + 1
    write (kicker_ele%name,   '(a, i0)') 'ROLL_Z', j_count
    kicker_ele%value(hkick$) =  val(angle$) * (1 - cos(val(roll$))) / 2
    kicker_ele%value(vkick$) = -val(angle$) * sin(val(roll$)) / 2
    val(roll$) = 0   ! So on next iteration will not create extra kickers.
    call insert_element (lat_out, kicker_ele, ix_ele, branch_out%ix_branch, orbit_out)
    call insert_element (lat_out, kicker_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
    cycle
  endif

  ! If there is a multipole component then put multipole elements at half strength 
  ! just before and just after the element. Exception: With Elegant if there is only one multipole for a quad, sextupole, or octupole

  monopole = .false.
  select case (ele%key)
  case (quadrupole$, sextupole$, octupole$, thick_multipole$)    ! Elegant
    call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, magnetic$, include_kicks$)
    if (count(knl /= 0) == 1 .and. all(knl(0:3) == 0)) monopole = .true.
  end select

  if (.not. monopole .and. ele%key /= multipole$ .and. ele%key /= ab_multipole$ .and. ele%key /= null_ele$ .and. ele%key /= sad_mult$) then
    call multipole_ele_to_ab (ele, .true., ix_pole_max, ab_ele%a_pole, ab_ele%b_pole)
    if (ix_pole_max > -1) then
      ab_ele%a_pole = ab_ele%a_pole / 2
      ab_ele%b_pole = ab_ele%b_pole / 2
      if (associated(ele%a_pole)) then
        deallocate (ele%a_pole, ele%b_pole)
        call attribute_bookkeeper(ele, .true.)
      endif
      j_count = j_count + 1
      write (ab_ele%name, '(a1, a, i0)') key_name(ele%key), 'MULTIPOLE_', j_count
      call insert_element (lat_out, ab_ele, ix_ele, branch_out%ix_branch, orbit_out)
      call insert_element (lat_out, ab_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
      cycle
    endif
  endif

  ! If there are nonzero kick values and this is not a kick type element then put
  ! kicker elements at half strength just before and just after the element.
  ! Also add a matrix element to get the change in z correct.
  ! A sad_mult gets translated to a matrix element which has kick components so no extra kickers needed here.

  if (has_hkick_attributes(ele%key)) then
    if (ele%key /= kicker$ .and. ele%key /= hkicker$ .and. ele%key /= vkicker$ .and. ele%key /= sad_mult$) then
      if (val(hkick$) /= 0 .or. val(vkick$) /= 0) then
        j_count = j_count + 1
        write (kicker_ele%name, '(a1, a, i0)') key_name(ele%key), '_KICKER_', j_count
        kicker_ele%value(hkick$) = val(hkick$) / 2
        kicker_ele%value(vkick$) = val(vkick$) / 2
        val(hkick$) = 0; val(vkick$) = 0
        if (ele%key == sbend$) then
          f = val(dg$) * val(l$) / 2
          kicker_ele%value(hkick$) = kicker_ele%value(hkick$) - cos(ele%value(ref_tilt_tot$)) * f
          kicker_ele%value(vkick$) = kicker_ele%value(vkick$) - sin(ele%value(ref_tilt_tot$)) * f
          val(dg$) = 0
        endif
        !!! write (taylor_ele%name, '(a, i0)') 'Z_SHIFTER', j_count 
        taylor_ele%name = ele%name
        call taylor_make_unit(taylor_ele%taylor)
        orb_start = orbit_out(ix_ele-1)
        orb_start%vec(2) = orb_start%vec(2) - kicker_ele%value(hkick$)
        orb_start%vec(4) = orb_start%vec(4) - kicker_ele%value(vkick$)
        call track1 (orb_start, ele, branch_out%param, orb_end) 
        ele%key = -1  ! Mark to ignore
        f = (ele%map_ref_orb_out%vec(5) - ele%map_ref_orb_in%vec(5)) - (orb_end%vec(5) - orb_start%vec(5))
        call add_taylor_term (taylor_ele%taylor(5), f, [0, 0, 0, 0, 0, 0])
        call insert_element (lat_out, kicker_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        call insert_element (lat_out, taylor_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
        call insert_element (lat_out, kicker_ele, ix_ele+3, branch_out%ix_branch, orbit_out)
        cycle
      endif
    endif
  endif

  ! A wiggler/undulator is translated to an Elegant CWIGGLER element if possible. Otherwise, and
  ! always for a sol_quad, a drift-matrix-drift model is used.

  if (ele%key == wiggler$ .or. ele%key == undulator$ .or. ele%key == sol_quad$) then

    ! Try for a CWIGGLER. If the wiggler has been sliced due to superposition, the CWIGGLER is
    ! constructed from the super_lord since a slice is not an integer number of wiggler periods.
    ! In this case the slices, and any markers in between, are replaced by the super_lord. This
    ! is only done if nothing else has been superimposed on the wiggler since such an element
    ! would be lost.

    if (ele%key /= sol_quad$ .and. .not. logic_option(.false., use_matrix_model)) then
      ix1 = -1; ix2 = -1
      wig_ok = .true.

      if (ele%slave_status == super_slave$) then
        lord => pointer_to_lord(ele, 1)
        call find_element_ends (lord, ele1, ele2)
        ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
        ! If the wiggler wraps around the origin we are in trouble.
        if (ix2 < ix1) then
          call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
          if (global_com%exit_on_error) call err_exit
        endif

        do i = ix1+1, ix2
          ele1 => branch_out%ele(i)
          if (ele1%key == -1 .or. ele1%key == marker$) cycle
          if (ele1%slave_status == super_slave$ .and. ele1%n_lord == 1) then
            ele2 => pointer_to_lord(ele1, 1)
            if (ele2%ix_ele == lord%ix_ele) cycle
          endif
          wig_ok = .false.
          exit
        enddo

        if (.not. wig_ok) call out_io (s_warn$, r_name, &
            'Wiggler has other elements superimposed on it so it cannot be translated to a CWIGGLER: ' // lord%name, &
            'Will use a drift-matrix-drift model instead.')
      else
        lord => ele
      endif

      if (wig_ok .and. (lord%key == wiggler$ .or. lord%key == undulator$)) then
        if (cwiggler_translate(lord, line_out, .false.)) then
          if (ele%slave_status == super_slave$) then
            call out_io (s_warn$, r_name, &
                'Note: Not translating to Elegant the markers within wiggler: ' // lord%name)
            wig_ele = lord
            wig_ele%slave_status = free$
            wig_ele%lord_status  = not_a_lord$
            wig_ele%n_lord = 0; wig_ele%n_slave = 0; wig_ele%n_lord_field = 0; wig_ele%n_slave_field = 0
            lord%key = -1  ! Mark for deletion
            do i = ix1+1, ix2
              branch_out%ele(i)%key = -1  ! Mark for deletion
            enddo
            call insert_element (lat_out, wig_ele, ix2+1, branch_out%ix_branch, orbit_out)
            ix_ele = ix2 + 1
          endif
          cycle
        endif
      endif
    endif

    ! Drift-matrix-drift model.

    call out_io (s_warn$, r_name, 'Converting element to drift-matrix-drift model: ' // ele%name)
    drift_ele%value = ele%value
    drift_ele%value(l$) = -val(l$) / 2
    call make_mat6 (drift_ele, branch_out%param)
    taylor_ele%mat6 = matmul(matmul(drift_ele%mat6, ele%mat6), drift_ele%mat6)
    call mat6_to_taylor (taylor_ele%vec0, taylor_ele%mat6, taylor_ele%taylor)

    ! Add drifts before and after wigglers and sol_quads so total length is invariant
    j_count = j_count + 1
    write (drift_ele%name, '(a, i0)') 'DRIFT_Z', j_count
    taylor_ele%name = ele%name
    drift_ele%value(l$) = val(l$) / 2
    ele%key = -1 ! Mark to ignore
    call insert_element (lat_out, drift_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
    call insert_element (lat_out, taylor_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
    call insert_element (lat_out, drift_ele, ix_ele+3, branch_out%ix_branch, orbit_out)
    ix_ele = ix_ele + 2
    cycle
  endif

enddo

! For a patch that is *not* associated with the edge of a solenoid: A z_offset must be split into a drift + patch

ix_ele = 0

do
  ix_ele = ix_ele + 1
  if (ix_ele > branch_out%n_ele_track) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle

  ! If the name has more than 16 characters then replace the name by something shorter and unique.

  orig_name = ele%name

  if (len_trim(ele%name) > 16) then
    i_unique = i_unique + 1
    write (ele%name(12:), '(i0)') i_unique
  endif

  ! Replace element name containing "/" or "#" with "_"

  do
    j = max(index(ele%name, '\'), index(ele%name, '#'))        ! '
    if (j == 0) exit
    ele%name(j:j) = '_'
  enddo

  if (ele%name /= orig_name .and. n_name_change_warn <= n_warn_max) then
    call out_io (s_info$, r_name, 'Element name changed from: ' // trim(orig_name) // ' to: ' // ele%name)
    if (n_name_change_warn == n_warn_max) call out_io (s_info$, r_name, &
                           'Enough name change warnings. Will stop issuing them now.')
    n_name_change_warn = n_name_change_warn + 1
  endif

  !

  val => ele%value

  if (ele%key == patch$ .and. ele%value(z_offset$) /= 0) then
    drift_ele%name = 'DRIFT_' // ele%name
    drift_ele%value(l$) = val(z_offset$)
    call insert_element (lat_out, drift_ele, ix_ele, branch_out%ix_branch, orbit_out)
    ix_ele = ix_ele + 1
    ele => branch_out%ele(ix_ele)
    val => ele%value
    val(z_offset$) = 0
  endif
enddo

!-------------------------------------------------------------------------------------------------
! Now write info to the output file...
! lat lattice name

write (iu, '(3a)') comment_char, ' File generated by: write_lattice_foreign_format', trim(eol_char)
write (iu, '(4a)') comment_char, ' Bmad lattice file: ', trim(lat%input_file_name), trim(eol_char)
if (lat%lattice /= '') write (iu, '(4a)') comment_char, ' Bmad lattice name: ', trim(lat%lattice), trim(eol_char)
write (iu, '(a)')

! write element parameters

n_names = 0                          ! number of names stored in the list
ix_ele = 0

do
  ix_ele = ix_ele + 1
  if (ix_ele > branch_out%n_ele_track) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle

  val => ele%value

  if (ele%key == elseparator$) then 
    n_elsep_warn = n_elsep_warn + 1
    ele%key = drift$
    call out_io (s_info$, r_name, 'Elseparator being converted into a drift for ELEGANT conversion: ' // ele%name)  
  endif

  ! Do not make duplicate specs

  call find_index (ele%name, names, an_indexx, n_names, ix_match)
  if (ix_match > 0) cycle

  ! Add to the list of elements

  if (size(names) < n_names + 10) then
    call re_allocate(names, 2*size(names))
    call re_allocate(an_indexx, 2*size(names))
  endif
 
  call find_index (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true.)

  !-------------------------------------------------------------------------------------------

  bmad_params = ''
  elegant_params = ''

  ! Special case where element is a pure thick N-pole of order greater than octupole.
  select case (ele%key)
  case (quadrupole$, sextupole$, octupole$, thick_multipole$)    ! Elegant
    call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, magnetic$, include_kicks$)
    if (count(knl /= 0) == 1 .and. all(knl(0:3) == 0)) then
      n = find_location(knl /= 0, .true.) - 1
      write (line_out, '(2a)') trim(ele%name) // ': mult'
      call value_to_line(line_out, knl(n), 'knl', 'R')
      write (line_out, '(2a, i0)') trim(line_out), ', order = ', n
      bmad_params(:5) = [character(40):: 'l', 'tilt', 'x_offset', 'y_offset', 'z_offset']
      elegant_params(:5) = [character(40):: 'l', 'tilt', 'dx', 'dy', 'dz']

      do i = 1, size(bmad_params)
        if (bmad_params(i) == '') exit
        call pointer_to_attribute (ele, upcase(bmad_params(i)), .true., a_ptr, err_flag)
        call value_to_line (line_out, a_ptr%r, elegant_params(i), 'R')
      enddo
      call write_line(line_out)
      cycle
    endif
  end select

  !

  select case (ele%key)

  case (instrument$, detector$, monitor$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': moni'
    bmad_params(:4) = [character(40):: 'l', 'tilt', 'x_offset', 'y_offset']
    elegant_params(:4) = [character(40):: 'l', 'tilt', 'dx', 'dy']


  case (drift$, pipe$)   ! Elegant
    if (ele%csr_method == off$) then
      write (line_out, '(2a)') trim(ele%name) // ': edrift'
    else
      write (line_out, '(2a)') trim(ele%name) // ': csredrift'
    endif

    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']

  case (gkicker$)
    write (line_out, '(2a)') trim(ele%name) // ': malign'
    bmad_params(:6) = [character(40):: 'x_kick', 'y_kick', 'z_kick', 'px_kick', 'py_kick', 'pz_kick']
    elegant_params(:6) = [character(40):: 'dx', 'dy', 'dz', 'dxp', 'dyp', 'dp']

  case (hkicker$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': ehkick'
    bmad_params(:6) = [character(40):: 'l', 'kick', 'tilt', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:6) = [character(40):: 'l', 'kick', 'tilt', 'dx', 'dy', 'dz']

  case (vkicker$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': evkick'
    bmad_params(:6) = [character(40):: 'l', 'kick', 'tilt', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:6) = [character(40):: 'l', 'kick', 'tilt', 'dx', 'dy', 'dz']

  case (kicker$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': ekicker'
    bmad_params(:7) = [character(40):: 'l', 'hkick', 'vkick', 'tilt', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:7) = [character(40):: 'l', 'hkick', 'vkick', 'tilt', 'dx', 'dy', 'dz']

  case (sbend$)   ! Elegant
    if (ele%csr_method == off$) then
      write (line_out, '(2a)') trim(ele%name) // ': csbend'
    else
      write (line_out, '(2a)') trim(ele%name) // ': csrcsbend'
    endif

    if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) line_out = trim(line_out) // ', malign_method = 2'

    select case (nint(ele%value(fringe_at$)))
    case (entrance_end$); line_out = trim(line_out) // ', edge2_effects = 0'
    case (exit_end$);     line_out = trim(line_out) // ', edge1_effects = 0'
    case (no_end$);       line_out = trim(line_out) // ', edge1_effects = 0, edge2_effects = 0'
    end select

    call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, magnetic$, include_kicks$)
    call value_to_line (line_out, ele%value(dg$)*ele%value(rho$), 'fse_dipole', 'R')
    call value_to_line (line_out, b_pole(0) - ele%value(dg$)*ele%value(l$), 'xkick', 'R')
    do n = 1, 8
      call value_to_line (line_out, b_pole(n)*factorial(n)/ele%value(l$), 'k' // int_str(n), 'R')
    enddo


    if (ele%value(fint$) == ele%value(fintx$)) then
      if (ele%value(fint$) /= 0.5_rp) call value_to_line (line_out, ele%value(fint$), 'fint', 'R', .false.)
    else
      if (ele%value(fint$) /= 0.5_rp) call value_to_line (line_out, ele%value(fint$), 'fint1', 'R', .false.)
      if (ele%value(fintx$) /= 0.5_rp) call value_to_line (line_out, ele%value(fintx$), 'fint2', 'R', .false.)
    endif

    bmad_params(:12) = [character(40):: 'l', 'angle', 'e1', 'e2', 'ref_tilt', 'roll', 'h1', 'h2', &
                                                             'vkick', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:12) = [character(40):: 'l', 'angle', 'e1', 'e2', 'tilt', 'etilt', 'h1', 'h2', 'ykick', 'dx', 'dy', 'dz']

  case (quadrupole$)   ! Elegant
    knl = knl / ele%value(l$)
    if (knl(2) == 0) then
      write (line_out, '(2a)') trim(ele%name) // ': kquad'
      if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) line_out = trim(line_out) // ', malign_method = 2'
    else
      write (line_out, '(2a)') trim(ele%name) // ': kquse'
      call value_to_line (line_out, 0.5_rp*knl(2)*cos(3*(tilts(2)-tilts(1))), 'k2', 'R')
    endif

    tilt = tilts(1)
    call value_to_line (line_out, knl(1), 'k1', 'R')

    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']

  case (sextupole$)   ! Elegant
    knl = knl / ele%value(l$)
    write (line_out, '(2a)') trim(ele%name) // ': ksext'
    if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) line_out = trim(line_out) // ', malign_method = 2'
    call value_to_line (line_out, knl(2), 'k2', 'R')
    call value_to_line (line_out, knl(1)*cos(0.5_rp*(tilts(1)-tilts(2))), 'k1', 'R')
    call value_to_line (line_out, knl(1)*sin(0.5_rp*(tilts(1)-tilts(2))), 'j1', 'R')
    call value_to_line (line_out, knl(0)*cos(tilts(0)), 'hkick', 'R')
    call value_to_line (line_out, knl(0)*sin(tilts(0)), 'vkick', 'R')

    tilt = tilts(2)
    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']

  case (octupole$)   ! Elegant
    knl = knl / ele%value(l$)
    write (line_out, '(2a)') trim(ele%name) // ': koct'
    call value_to_line (line_out, knl(3), 'k3', 'R')
    call value_to_line (line_out, knl(0)*cos(tilts(0)), 'hkick', 'R')
    call value_to_line (line_out, knl(0)*sin(tilts(0)), 'vkick', 'R')

    tilt = tilts(3)
    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']

  case (solenoid$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': sole'
    bmad_params(:5) = [character(40):: 'l', 'ks', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:5) = [character(40):: 'l', 'ks', 'dx', 'dy', 'dz']

  case (taylor$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': ematrix'
    do i = 1, 6
      f = taylor_coef(ele%taylor(i), [0,0,0,0,0,0])
      call value_to_line (line_out, f, 'c' // int_str(i), 'R')

      do j = 1, 6
        f = taylor_coef(ele%taylor(i), taylor_expn([j]))
        call value_to_line (line_out, f, 'r' // int_str(i) // int_str(j), 'R')

        do k = 1, j
          f = taylor_coef(ele%taylor(i), taylor_expn([j,k]))
          call value_to_line (line_out, f, 'r' // int_str(i) // int_str(j), 'R')
        enddo
      enddo
    enddo

    tilt = ele%value(tilt$)
    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']

  case (beambeam$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': beambeam'
    call value_to_line (line_out, strong_beam_strength(ele)*e_charge, 'charge', 'R')
    bmad_params(:4) = [character(40):: 'x_offset', 'y_offset', 'sig_x', 'sig_y']
    elegant_params(:4) = [character(40):: 'xcenter', 'ycenter', 'xsize', 'ysize']

  case (marker$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': mark'
    bmad_params(:2) = [character(40):: 'x_offset', 'y_offset']
    elegant_params(:2) = [character(40):: 'dx', 'dy']

  case (ab_multipole$, multipole$)   ! Elegant
    call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, include_kicks$)
    orig_name = ele%name
    ab_ele = ele
    ele%key = -1

    do i = 1, ix_pole_max
      if (knl(i) == 0) cycle
      if (count(knl /= 0) > 1) ab_ele%name = trim(orig_name) // '__' // int_str(i)
      write (line_out, '(2a)') trim(ab_ele%name) // ': mult'
      call insert_element(lat_out, ab_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
      ix_ele = ix_ele + 1
      call value_to_line (line_out, knl(i), 'knl', 'R')
      call value_to_line (line_out, tilts(i), 'tilt', 'R')
      line_out = trim(line_out) // ', order = ' // int_str(i)
      call value_to_line (line_out, ab_ele%value(x_offset$), 'dx', 'R')
      call value_to_line (line_out, ab_ele%value(y_offset$), 'dy', 'R')
      call value_to_line (line_out, ab_ele%value(z_offset$), 'dz', 'R')
      call write_line (line_out)
    enddo
    cycle  

  case (ecollimator$, rcollimator$)   ! Elegant
    if (ele%key == ecollimator$) then
      write (line_out, '(2a)') trim(ele%name) // ': ecol'
    else
      write (line_out, '(2a)') trim(ele%name) // ': rcol'
    endif
    call value_to_line (line_out, ab_ele%value(l$), 'l', 'R')

    r_max = (ele%value(x2_limit$) + ele%value(x1_limit$)) / 2
    r0 = (ele%value(x2_limit$) - ele%value(x1_limit$)) / 2
    if (ele%offset_moves_aperture) r0 = r0 + ele%value(x_offset$)
    call value_to_line (line_out, r_max, 'x_max', 'R')
    call value_to_line (line_out, r0, 'dx', 'R')

    r_max = (ele%value(y2_limit$) + ele%value(y1_limit$)) / 2
    r0 = (ele%value(y2_limit$) - ele%value(y1_limit$)) / 2
    if (ele%offset_moves_aperture) r0 = r0 + ele%value(y_offset$)
    call value_to_line (line_out, r_max, 'y_max', 'R')
    call value_to_line (line_out, r0, 'dy', 'R')

  case (wiggler$, undulator$)   ! Elegant
    ! The substitution loop above has already verified that the translation is possible.
    if (.not. cwiggler_translate(ele, line_out, .true.)) then
      call out_io (s_error$, r_name, 'INTERNAL ERROR: CANNOT TRANSLATE WIGGLER TO A CWIGGLER: ' // trim(ele%name), &
                                     'PLEASE REPORT THIS. CONVERTING TO A DRIFT.')
      write (line_out, '(2a)') trim(ele%name) // ': edrift'
      bmad_params(:1) = [character(40):: 'l']
      elegant_params(:1) = [character(40):: 'l']
    endif

  case (rfcavity$, lcavity$)   ! Elegant
    if (ele%key == rfcavity$) then
      write (line_out, '(2a)') trim(ele%name) // ': rfca'
      call value_to_line (line_out, 360.0_rp*(ele%value(phi0$)+ele%value(phi0_multipass$)), 'phase', 'R')
    else
      write (line_out, '(2a)') trim(ele%name) // ': rfca, change_p0 = 1'
      call value_to_line (line_out, 360.0_rp*(ele%value(phi0$)+ele%value(phi0_multipass$))+90.0_rp, 'phase', 'R')
    endif

    if (nint(ele%value(cavity_type$)) == standing_wave$) then
      line_out = trim(line_out) // ', body_focus_model="SRS", standing_wave = 1, end1_focus=1, end2_focus=1'
    endif

    bmad_params(:3) = [character(40):: 'l', 'voltage', 'rf_frequency']
    elegant_params(:3) = [character(40):: 'l', 'volt', 'freq']

  case (crab_cavity$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': rfdf'
    call value_to_line (line_out, 360.0_rp*(ele%value(phi0$)+ele%value(phi0_multipass$)), 'phase', 'R')
    bmad_params(:7) = [character(40):: 'l', 'voltage', 'rf_frequency', 'tilt', 'x_offset', 'y_offset', 'z_offset']
    elegant_params(:7) = [character(40):: 'l', 'voltage', 'frequency', 'tilt', 'dx', 'dy', 'dz']

  case (patch$)   ! Elegant
    if (ele%value(tilt$) /= 0) then
      write (line_out, '(2a)') trim(ele%name) // '_rot: rotate'
      call value_to_line (line_out, ele%value(tilt$),  'tilt', 'R')
    endif

    write (line_out, '(2a)') trim(ele%name) // ': malign'
    bmad_params(:7) = [character(40):: 'x_offset', 'y_offset', 'z_offset', 't_offset', 'e_tot_offset', 'x_pitch', 'y_pitch']
    elegant_params(:7) = [character(40):: 'dx', 'dy', 'dz', 'dt', 'de', 'dxp', 'dyp']

  case (floor_shift$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': floor'
    call value_to_line (line_out, ele%floor%r(1),  'x', 'R')
    call value_to_line (line_out, ele%floor%r(2),  'y', 'R')
    call value_to_line (line_out, ele%floor%r(3),  'z', 'R')
    call value_to_line (line_out, ele%floor%theta, 'theta', 'R')
    call value_to_line (line_out, ele%floor%phi,   'phi', 'R')
    call value_to_line (line_out, ele%floor%psi,   'psi', 'R')

  case (match$)   ! Elegant
    write (line_out, '(2a)') trim(ele%name) // ': ematrix'
    call value_to_line (line_out, ele%vec0(1),  'C1', 'R')
    call value_to_line (line_out, ele%vec0(2),  'C2', 'R')
    call value_to_line (line_out, ele%vec0(3),  'C3', 'R')
    call value_to_line (line_out, ele%vec0(4),  'C4', 'R')
    call value_to_line (line_out, ele%vec0(5),  'C5', 'R')
    call value_to_line (line_out, ele%vec0(6),  'C6', 'R')
    do i = 1, 6; do j = 1, 6
      call value_to_line (line_out, ele%mat6(i,j),  'R' // num(i) // num(j), 'R')
    enddo; enddo

  case default
    call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO TRANSLATE ELEMENT: ' // ele%name, &
                                   'WHICH IS OF TYPE: ' // key_name(ele%key), &
                                   'CONVERTING TO DRIFT')
    write (line_out, '(2a)') trim(ele%name) // ': drift'
    bmad_params(:1) = [character(40):: 'l']
    elegant_params(:1) = [character(40):: 'l']
  end select

  !--------------------------------------------------------

  select case (ele%key)
  case (sbend$, patch$, drift$)   ! Elegant
    ! Pass

  case (quadrupole$, sextupole$, octupole$, taylor$)    ! Elegant
    x_pitch = ele%value(x_pitch$)
    y_pitch = ele%value(y_pitch$)
    call floor_angles_to_w_mat(x_pitch, y_pitch, tilt, w_mat)

    if (x_pitch == 0 .or. y_pitch == 0) then
      epitch = -y_pitch  ! alpha_x
      eyaw = x_pitch     ! alpha_y
      etilt = tilt       ! alpha_z
    else
      epitch = -atan2(w_mat(2,3), w_mat(3,3))
      etilt = -atan2(w_mat(1,2), w_mat(1,1))
      eyaw = -atan2(w_mat(1,3), w_mat(2,3)/sin(epitch))
    endif

    offset = matmul(w_mat, [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)])
    call value_to_line (line_out, etilt, 'tilt', 'R')
    call value_to_line (line_out, epitch, 'pitch', 'R')
    call value_to_line (line_out, eyaw, 'yaw', 'R')
    call value_to_line (line_out, offset(1), 'dx', 'R')
    call value_to_line (line_out, offset(2), 'dy', 'R')
    call value_to_line (line_out, offset(3), 'dz', 'R')

  ! Elegant
  
  case (instrument$, detector$, monitor$, hkicker$, vkicker$, kicker$, wiggler$, undulator$)  ! Has tilt but not pitches.
    if (has_orientation_attributes(ele) .and. (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0)) then
      call out_io (s_warn$, r_name, 'X_PITCH OR Y_PITCH PARAMETERS OF A ' // trim(key_name(ele%key)) // ' CANNOT BE TRANSLATED TO ELEGANT: ' // ele%name)
    endif
  case default
    if (has_orientation_attributes(ele) .and. (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0)) then
      call out_io (s_warn$, r_name, 'TILT, X_PITCH OR Y_PITCH PARAMETERS OF A ' // trim(key_name(ele%key)) // ' CANNOT BE TRANSLATED TO ELEGANT: ' // ele%name)
    endif
  end select

  !--------------------------------------------------------

  do i = 1, size(bmad_params)
    if (bmad_params(i) == '') exit
    call pointer_to_attribute (ele, upcase(bmad_params(i)), .true., a_ptr, err_flag)
    call value_to_line (line_out, a_ptr%r, elegant_params(i), 'R')
  enddo

  call write_line(line_out)
  cycle
enddo

!---------------------------------------------------------------------------------------
! Write the lattice line

i_unique = 1000
i_line = 0
init_needed = .true.
line = ' '

do n = 1, branch_out%n_ele_track
  ele => branch_out%ele(n)
  if (ele%key == -1) cycle

  if (init_needed) then
    write (iu, '(a)')
    write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
    write (iu, '(a)')
    i_line = i_line + 1
    write (line_out, '(a, i0, 2a)') 'line_', i_line, ': line = (', ele%name
    iout = 0
    init_needed = .false.

  else
    ix = len_trim(line_out) + len_trim(ele%name)

    if (ix > 75) then
      write (iu, '(3a)') trim(line_out), trim(separator_char), trim(continue_char)
      iout = iout + 1
      line_out = '   ' // ele%name
    else
      line_out = trim(line_out) // trim(separator_char) // ' ' // ele%name
    endif
  endif

  if (ele%key == patch$ .and. ele%value(tilt$) /= 0) then
    line_out = trim(line_out) // ', ' //  trim(ele%name) // '_rot'
  endif

  ! Output line if long enough or at end

  if (n == branch_out%n_ele_track .or. iout > 48) then
    line_out = trim(line_out) // ')'
    write (iu, '(2a)') trim(line_out), trim(eol_char)
    line_out = ' '
    init_needed = .true.
  endif

enddo

!------------------------------------------
! Use statement

write (iu, '(a)')
write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
write (iu, '(a)')

line_out = 'lat: line = (line_1'

do i = 2, i_line
  write (line_out, '(3a, i0)') trim(line_out), trim(separator_char), ' line_', i
enddo

line_out = trim(line_out) // ')'
call write_line (line_out)

! Write twiss parameters for a non-closed lattice.

! End stuff

call out_io (s_info$, r_name, 'Written ELEGANT lattice file: ' // trim(out_file_name))

deallocate (names)
if (present(err)) err = .false.

call deallocate_lat_pointers (lat_out)
if (associated(wig_cmap%ptr)) deallocate (wig_cmap%ptr)

! Restore ptc settings

if (n_taylor_order_saved /= ptc_private%taylor_order_ptc) call set_ptc (taylor_order = n_taylor_order_saved) 
ptc_com%exact_model = ptc_exact_model

close(iu)

!------------------------------------------------------------------------
contains

subroutine write_line (line_out)

implicit none

character(*) line_out
integer ix, ix1, ix2, ix3

! Prefer to breakup a line after a comma

do
  if (len_trim(line_out) < ix_line_max) exit
  ix1 = index(line_out(ix_line_min+1:), ',')
  ix2 = index(line_out(ix_line_min+1:), '=')
  ix3 = index(line_out(ix_line_min+1:), ' ')

  if (ix1 /= 0 .and. ix1+ix_line_min < ix_line_max) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0 .and. ix2+ix_line_min < ix_line_max) then
    ix = ix2 + ix_line_min
  elseif (ix3 /= 0 .and. ix3+ix_line_min < ix_line_max) then
    ix = ix3 + ix_line_min
  elseif (ix1 /= 0) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0) then
    ix = ix2 + ix_line_min
  else
    ix = ix3 + ix_line_min
  endif

  write (iu, '(2a)') line_out(:ix), trim(continue_char)
  line_out = '    ' // line_out(ix+1:)
enddo

write (iu, '(2a)') trim(line_out), trim(eol_char)

end subroutine write_line

!------------------------------------------------------------------------
! contains
!+
! Function cwiggler_translate (ele, line_out, do_write) result (is_ok)
!
! Translate a wiggler or undulator to an Elegant CWIGGLER element.
!
! Elegant's CWIGGLER uses the same Cartesian harmonic field expansion that a Bmad cartesian_map
! uses so the translation is exact. With (x, y, z) measured from the entrance end of the element,
! a Bmad term of the family_y/hyper_y variety gives (see em_field_calc):
!   Bx = -coef * (kx/ky) * sin(kx*(x+x0)) * sinh(ky*(y+y0)) * cos(kz*z + phi_z)
!   By =  coef *          cos(kx*(x+x0)) * cosh(ky*(y+y0)) * cos(kz*z + phi_z)
! while Elegant's horizontal wiggler with normal poles gives (see routine GWigB in the Elegant
! file gwig.c):
!   Bx =  B_MAX * Cmn * (kx/ky) * sin(kx*x) * sinh(ky*y) * cos(kz*z + Phase)
!   By = -B_MAX * Cmn *          cos(kx*x) * cosh(ky*y) * cos(kz*z + Phase)
! So the two are the same expansion with B_MAX * Cmn = -coef. Doing this for all four of the
! Elegant expansions gives the table:
!
!   Bmad family  Bmad form  Elegant expansion                   B_MAX * Cmn   k constraint
!   -----------  ---------  ---------------------------------   -----------   ---------------
!   family_y     hyper_y    BY_FILE with BY_SPLIT_POLE = 0       -coef         ky^2 = kx^2+kz^2
!   family_y     hyper_x    BY_FILE with BY_SPLIT_POLE = 1       -coef*ky/kx   kx^2 = ky^2+kz^2
!   family_x     hyper_x    BX_FILE with BX_SPLIT_POLE = 0        coef         kx^2 = ky^2+kz^2
!   family_x     hyper_y    BX_FILE with BX_SPLIT_POLE = 1        coef*kx/ky   ky^2 = kx^2+kz^2
!
! B_MAX is set to 1 Tesla so that the Cmn column of a harmonic file is simply the field in Tesla.
! This avoids any normalization ambiguity and makes Elegant's radiation integral estimate, which
! uses B_MAX * sum(Cmn), come out right.
!
! Elegant has no analog of a Bmad hyper_xy form term nor of a family_qu or family_sq term so a
! map using any of these cannot be translated. Nor can a cylindrical_map, grid_field, or
! gen_gradients wiggler. In these cases is_ok is set False and the calling code falls back to a
! drift-matrix-drift model.
!
! Periodic type (planar_model and helical_model) wigglers do not store a cartesian_map but they
! do have an exactly equivalent one which create_wiggler_cartesian_map constructs. This is the
! same map that the PTC and symp_lie_bmad tracking use.
!
! Input:
!   ele       -- ele_struct: Wiggler or undulator element.
!   do_write  -- logical: If True, write the harmonic file(s) and construct the element line.
!                  If False, only test whether the translation is possible.
!
! Output:
!   line_out  -- character(*): Elegant element definition. Only set if do_write is True.
!   is_ok     -- logical: True if the element can be translated to a CWIGGLER.
!-

function cwiggler_translate (ele, line_out, do_write) result (is_ok)

type (ele_struct), target :: ele
type (ele_struct), pointer :: f_ele
type (cartesian_map_struct), pointer :: cm
type (cartesian_map_term1_struct), pointer :: ct

real(rp) scale, s_anchor, z0, kw, kz_min, kz_max, len_ele, ph, dn, ff, x_off, y_off
real(rp), allocatable :: c_amp(:), phase(:)

integer it, jt, iw, ix_swap, n_term, n_row, n_period, n_step, n_by, n_bx, form_by, form_bx, iu2, ios2
integer, allocatable :: ib(:), indx(:)

character(*) line_out
character(200) by_file, bx_file, fname, full_name

logical do_write, is_ok

! Get the element that holds the field. For a super_slave this is the lord and so the CWIGGLER
! would not correspond to the slave. The calling code resolves this before calling here.

is_ok = .false.
len_ele = ele%value(l$)
if (len_ele == 0) return

f_ele => pointer_to_field_ele(ele, 1)
if (f_ele%value(l$) /= len_ele) then
  call cwig_err (ele, 'GETS ITS FIELD FROM AN ELEMENT OF A DIFFERENT LENGTH.', do_write)
  return
endif

! Get the Cartesian map.

select case (f_ele%field_calc)
case (planar_model$, helical_model$)
  if (f_ele%value(l_period$) == 0) then
    call cwig_err (ele, 'IS A PERIODIC TYPE WIGGLER WITH L_PERIOD = 0.', do_write)
    return
  endif
  call create_wiggler_cartesian_map (f_ele, wig_cmap)
  cm => wig_cmap

case (fieldmap$)
  if (associated(f_ele%cylindrical_map) .or. associated(f_ele%grid_field) .or. associated(f_ele%gen_gradients)) then
    call cwig_err (ele, 'HAS A FIELD MAP THAT IS NOT A CARTESIAN_MAP.', do_write)
    return
  endif
  if (.not. associated(f_ele%cartesian_map)) then
    call cwig_err (ele, 'HAS FIELD_CALC = FIELDMAP BUT NO FIELD MAP.', do_write)
    return
  endif
  if (size(f_ele%cartesian_map) /= 1) then
    call cwig_err (ele, 'HAS MORE THAN ONE CARTESIAN_MAP.', do_write)
    return
  endif
  cm => f_ele%cartesian_map(1)

case default
  call cwig_err (ele, 'HAS A FIELD_CALC SETTING THAT CANNOT BE TRANSLATED: ' // int_str(f_ele%field_calc), do_write)
  return
end select

if (cm%field_type /= magnetic$) then
  call cwig_err (ele, 'HAS AN ELECTRIC CARTESIAN_MAP.', do_write)
  return
endif

if (.not. associated(cm%ptr)) return
if (.not. allocated(cm%ptr%term)) return
if (size(cm%ptr%term) == 0) return

! Sort the terms into the By (family_y) and Bx (family_x) sets and compute the Cmn and Phase
! columns. Elegant has one split-pole switch per set so all terms of a set must have the same form.

n_term = size(cm%ptr%term)
allocate (c_amp(n_term), phase(n_term), ib(n_term), indx(n_term))

scale = cm%field_scale * master_parameter_value(cm%master_parameter, ele)

select case (cm%ele_anchor_pt)
case (anchor_center$); s_anchor = len_ele / 2
case (anchor_end$);    s_anchor = len_ele
case default;          s_anchor = 0
end select
z0 = s_anchor + cm%r0(3)   ! Bmad map z = (distance from entrance) - z0.

n_by = 0; n_bx = 0
form_by = 0; form_bx = 0
kz_min = 1e30_rp; kz_max = 0

do it = 1, n_term
  ct => cm%ptr%term(it)

  if (ct%kz <= 0) then
    call cwig_err (ele, 'HAS A CARTESIAN_MAP TERM WITH KZ <= 0 WHICH ELEGANT DOES NOT ALLOW.', do_write)
    return
  endif

  if (ct%x0 /= cm%ptr%term(1)%x0 .or. ct%y0 /= cm%ptr%term(1)%y0) then
    call cwig_err (ele, 'HAS CARTESIAN_MAP TERMS WITH DIFFERING X0/Y0 OFFSETS.', do_write)
    return
  endif

  ! Elegant checks the Maxwell constraint on the wave numbers to a relative accuracy of 1d-6.

  select case (ct%form)
  case (hyper_y$); ff = sqrt(ct%kx**2 + ct%kz**2) / ct%ky
  case (hyper_x$); ff = sqrt(ct%ky**2 + ct%kz**2) / ct%kx
  case default
    call cwig_err (ele, 'HAS A HYPER_XY FORM CARTESIAN_MAP TERM WHICH ELEGANT CANNOT MODEL.', do_write)
    return
  end select

  if (abs(ff - 1) > 1d-8) then
    call cwig_err (ele, 'HAS A CARTESIAN_MAP TERM WHOSE WAVE NUMBERS DO NOT SATISFY THE MAXWELL CONSTRAINT.', do_write)
    return
  endif

  select case (ct%family)
  case (family_y$)
    n_by = n_by + 1
    ib(it) = 1
    if (form_by == 0) form_by = ct%form
    if (ct%form /= form_by) then
      call cwig_err (ele, 'HAS BOTH HYPER_Y AND HYPER_X FORM FAMILY_Y TERMS. ELEGANT HAS ONE BY_SPLIT_POLE SWITCH.', do_write)
      return
    endif
    if (ct%form == hyper_y$) then
      c_amp(it) = -scale * ct%coef
    else
      if (ct%ky == 0) then
        call cwig_err (ele, 'HAS A SPLIT POLE HORIZONTAL WIGGLER TERM WITH KY = 0 WHICH ELEGANT DOES NOT ALLOW.', do_write)
        return
      endif
      c_amp(it) = -scale * ct%coef * ct%ky / ct%kx
    endif

  case (family_x$)
    n_bx = n_bx + 1
    ib(it) = 2
    if (form_bx == 0) form_bx = ct%form
    if (ct%form /= form_bx) then
      call cwig_err (ele, 'HAS BOTH HYPER_X AND HYPER_Y FORM FAMILY_X TERMS. ELEGANT HAS ONE BX_SPLIT_POLE SWITCH.', do_write)
      return
    endif
    if (ct%form == hyper_x$) then
      c_amp(it) = scale * ct%coef
    else
      if (ct%kx == 0) then
        call cwig_err (ele, 'HAS A SPLIT POLE VERTICAL WIGGLER TERM WITH KX = 0 WHICH ELEGANT DOES NOT ALLOW.', do_write)
        return
      endif
      c_amp(it) = scale * ct%coef * ct%kx / ct%ky
    endif

  case default
    call cwig_err (ele, 'HAS A QU OR SQ FAMILY CARTESIAN_MAP TERM WHICH ELEGANT CANNOT MODEL.', do_write)
    return
  end select

  ph = modulo(ct%phi_z - ct%kz * z0, twopi)
  if (ph >= twopi) ph = 0   ! Guard against roundoff. Elegant demands a non-negative first phase.
  phase(it) = ph

  kz_min = min(kz_min, ct%kz)
  kz_max = max(kz_max, ct%kz)
enddo

! Elegant normalizes the wave numbers to kw = twopi/(L/PERIODS) with PERIODS an integer. Also the
! integration step size must resolve the highest harmonic: Elegant demands
! pi * kz_max / kw <= STEPS_PER_PERIOD and that STEPS_PER_PERIOD be a multiple of 4.

dn = len_ele * kz_min / twopi
n_period = max(1, nint(dn))
kw = twopi * n_period / len_ele

if (do_write .and. abs(dn - n_period) > 1d-6 * dn) then
  call out_io (s_warn$, r_name, 'Wiggler length is not an integer number of field periods: ' // trim(ele%name), &
                                'The Elegant CWIGGLER PERIODS parameter will be set to ' // int_str(n_period) // '.')
endif

n_step = max(12, 4 * ceiling(pi * kz_max / (4 * kw)))

if (.not. do_write) then
  is_ok = .true.
  return
endif

! Write the harmonic file(s). The rows are ordered by increasing kz since Elegant treats the first
! row as the fundamental.

indx = [(it, it = 1, n_term)]
do it = 2, n_term
  do jt = it, 2, -1
    if (cm%ptr%term(indx(jt))%kz >= cm%ptr%term(indx(jt-1))%kz) exit
    ix_swap = indx(jt); indx(jt) = indx(jt-1); indx(jt-1) = ix_swap
  enddo
enddo

by_file = ''; bx_file = ''

do iw = 1, 2
  if (iw == 1) then
    if (n_by == 0) cycle
    n_row = n_by
    fname = trim(sdds_base) // '.' // trim(ele%name) // '.by.sdds'
  else
    if (n_bx == 0) cycle
    n_row = n_bx
    fname = trim(sdds_base) // '.' // trim(ele%name) // '.bx.sdds'
  endif

  call fullfilename (fname, full_name)
  iu2 = lunget()
  open (iu2, file = full_name, iostat = ios2)
  if (ios2 /= 0) then
    call out_io (s_error$, r_name, 'CANNOT OPEN CWIGGLER HARMONIC FILE: ' // trim(fname))
    return
  endif

  write (iu2, '(a)') 'SDDS1'
  write (iu2, '(4a)') '&description text="CWIGGLER harmonics from the Bmad cartesian_map of element ', &
                                                            trim(ele%name), '", contents="cwiggler harmonics" &end'
  write (iu2, '(a)') '&column name=Cmn, units=T, type=double &end'
  write (iu2, '(a)') '&column name=KxOverKw, type=double &end'
  write (iu2, '(a)') '&column name=KyOverKw, type=double &end'
  write (iu2, '(a)') '&column name=KzOverKw, type=double &end'
  write (iu2, '(a)') '&column name=Phase, units=rad, type=double &end'
  write (iu2, '(a)') '&data mode=ascii, &end'
  write (iu2, '(i0)') n_row

  do jt = 1, n_term
    it = indx(jt)
    if (ib(it) /= iw) cycle
    ct => cm%ptr%term(it)
    write (iu2, '(5(a, 1x))') re_str(c_amp(it)), re_str(ct%kx/kw), re_str(ct%ky/kw), re_str(ct%kz/kw), re_str(phase(it))
  enddo

  close (iu2)

  if (iw == 1) then
    by_file = fname
  else
    bx_file = fname
  endif
enddo

! Construct the element definition. FORCE_MATCHED is turned off since it would inset the field
! from the ends of the element and so give a field different from the Bmad one.

if (.not. wig_note_printed) then
  write (iu, '(3a)') comment_char, &
      ' CWIGGLER B_MAX is set to 1 so that the Cmn column of a harmonic file is the field in Tesla.', trim(eol_char)
  wig_note_printed = .true.
endif

write (line_out, '(2a)') trim(ele%name), ': cwiggler, b_max = 1, integration_order = 4, force_matched = 0'
call value_to_line (line_out, len_ele, 'l', 'R')
line_out = trim(line_out) // ', periods = ' // int_str(n_period) // ', steps_per_period = ' // int_str(n_step)

if (by_file /= '') then
  line_out = trim(line_out) // ', by_file = "' // trim(by_file) // '"'
  if (form_by == hyper_x$) line_out = trim(line_out) // ', by_split_pole = 1'
endif

if (bx_file /= '') then
  line_out = trim(line_out) // ', bx_file = "' // trim(bx_file) // '"'
  if (form_bx == hyper_y$) line_out = trim(line_out) // ', bx_split_pole = 1'
endif

! The map r0 and term x0/y0 shift the field within the element so they add to the element offset.

x_off = ele%value(x_offset$) + cm%r0(1) - cm%ptr%term(1)%x0
y_off = ele%value(y_offset$) + cm%r0(2) - cm%ptr%term(1)%y0

call value_to_line (line_out, x_off, 'dx', 'R')
call value_to_line (line_out, y_off, 'dy', 'R')
call value_to_line (line_out, ele%value(z_offset$), 'dz', 'R')
call value_to_line (line_out, ele%value(tilt$), 'tilt', 'R')

is_ok = .true.

end function cwiggler_translate

!------------------------------------------------------------------------
! contains
!
! Print a rate limited explanation of why a wiggler cannot be translated to a CWIGGLER.
! Only printed on the testing pass so that the message appears once per element instance.

subroutine cwig_err (ele, why, do_write)

type (ele_struct) ele
character(*) why
logical do_write

!

if (do_write) return

n_wig_model_err = n_wig_model_err + 1
if (n_wig_model_err > print_wig_model_err_max) return

call out_io (s_warn$, r_name, 'Cannot translate to an Elegant CWIGGLER the element: ' // trim(ele%name), &
                              'Since this element ' // trim(why), &
                              'Will use a drift-matrix-drift model instead.')

if (n_wig_model_err == print_wig_model_err_max) call out_io (s_warn$, r_name, &
                              'Enough wiggler translation warnings. Will stop issuing them now.')

end subroutine cwig_err

end subroutine write_lattice_elegant_format
