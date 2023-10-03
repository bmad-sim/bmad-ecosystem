!+ 
! Subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
!        use_matrix_model, include_apertures, dr12_drift_max, ix_start, ix_end, ix_branch, converted_lat, err)
!
! Subroutine to write a Elegant, MAD-8, MAD-X, OPAL, SAD, or XSIF lattice file using the 
! information in a lat_struct. Optionally, only part of the lattice can be generated.
! [XSIF is a variant of MAD8 used by SLAC.]
!
!=========================================================================================
! NOTE: ELEGANT TRANSLATION IN DEVELOPMENT. PLEASE CONTACT DAVID SAGAN IF YOU WANT TO USE!
!=========================================================================================
!
! To write a Bmad lattice file, use: write_bmad_lattice_file
!
! Note: When translating to XSIF or MAD8: sad_mult and patch element are translated
!  to a XSIF/MAD8 matrix element (which is a 2nd order map). In this case, the ref_orbit orbit is
!  used as the reference orbit for construction of the 2nd order map.
!
! If a sad_mult or patch element is translated to a matrix element, and the referece orbit
! is non-zero, the calculation must use 2nd order maps thourghout in order to avoid "feed down".
! If the PTC map order is different from 2, PTC will be temperarily switched to 2. 
!
! The MAD drift model is approximate and this can be a problem if the reference orbit is large.
! For a drift, the value of transfer matrix element R12 is equal to L/(1+pz) for small
! deviations of the ref_orbit from zero. dr12_drift_max sets the maximum deviation of R12 beyound 
! which an extra matrix element is inserted to make the MAD model better agree with Bmad.
!
! Note: sol_quad elements are replaced by a drift-matrix-drift or solenoid-quad model.
! Note: wiggler elements are replaced by a drift-matrix-drift or drift-bend model.
!
! Input:
!   out_type          -- character(*): Either 'ELEGANT', 'XSIF', 'MAD-8', 'MAD-X', 'SAD', or 'OPAL-T'.
!   out_file_name     -- character(*): Name of the mad output lattice file.
!   lat               -- lat_struct: Holds the lattice information.
!   ref_orbit(0:)     -- coord_struct, allocatable, optional: Referece orbit for sad_mult and patch elements.
!                          This argument must be present if the lattice has sad_mult or patch elements and is
!                          being translated to MAD-8 or SAD.
!   use_matrix_model  -- logical, optional: Use a drift-matrix_drift model for wigglers/undulators?
!                           [A MAD "matrix" is a 2nd order Taylor map.] This switch is ignored for SAD conversion.
!                           Default is False -> Use a bend-drift-bend model. 
!                           Note: sol_quad elements always use a drift-matrix-drift model.
!   include_apertures -- logical, optional: If True (the default), add to the output lattice a zero length
!                           collimator element next to any non-collimator element that has an aperture.
!                           Note: MADX translations for non-drift elements can handle non-collimator elements 
!                           with an aperture so in this case this argument is ignored.
!   dr12_drift_max    -- real(rp), optional: Max deviation for drifts allowed before a correction matrix element
!                           is added. Default value is 1d-5. A negative number means use default.
!   ix_start          -- integer, optional: Starting index of lat%ele(i)
!                           used for output.
!   ix_end            -- integer, optional: Ending index of lat%ele(i)
!                           used for output.
!   ix_branch         -- Integer, optional: Index of lattice branch to use. Default = 0.
!
! Output:
!   converted_lat     -- lat_struct, optional: Equivalent Bmad lattice with wiggler and 
!                           sol_quad elements replaced by their respective models.
!   err               -- logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
      use_matrix_model, include_apertures, dr12_drift_max, ix_start, ix_end, ix_branch, converted_lat, err)

use mad_mod, dummy2 => write_lattice_in_foreign_format
use bmad, dummy => write_lattice_in_foreign_format
use write_lat_file_mod, dummy3 => write_lattice_in_foreign_format
use ptc_interface_mod, only: taylor_inverse, concat_taylor

implicit none

type (lat_struct), target :: lat, lat_model, lat_out
type (lat_struct), optional, target :: converted_lat
type (ele_struct), pointer :: ele, ele1, ele2, lord, sol_ele, first_sol_edge
type (ele_struct) :: drift_ele, ab_ele, taylor_ele, col_ele, kicker_ele, null_ele, bend_ele, quad_ele
type (coord_struct) orb_start, orb_end, orb_center
type (coord_struct), allocatable, optional :: ref_orbit(:)
type (coord_struct), allocatable :: orbit_out(:)
type (taylor_term_struct) :: term
type (branch_struct), pointer :: branch, branch_out
type (mad_energy_struct) energy
type (mad_map_struct) mad_map
type (taylor_struct) taylor_a(6), taylor_b(6)
type (taylor_struct), pointer :: taylor_ptr(:)
type (all_pointer_struct) a_ptr

real(rp), optional :: dr12_drift_max
real(rp) field, hk, vk, limit(2), length, a, b, f, e2, beta, r_max, r0, dr12_max
real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) tilt, x_pitch, y_pitch, etilt, epitch, eyaw, offset(3), w_mat(3,3)

integer, optional :: ix_start, ix_end, ix_branch
integer, allocatable :: n_repeat(:), an_indexx(:)
integer i, j, ib, j2, k, n, ix, i_unique, i_line, iout, iu, n_names, j_count, f_count, ix_ele
integer ie, ie1, ie2, ios, a_count, ix_lord, ix_match, iv, ifa, ix_pole_max
integer ix1, ix2, n_lord, aperture_at, n_name_change_warn, n_elsep_warn, n_taylor_order_saved
integer :: ix_line_min, ix_line_max, n_warn_max = 10

character(*), parameter :: r_name = "write_lattice_in_foreign_format"
character(*) out_type, out_file_name
character(300) line, knl_str, ksl_str
character(40) orig_name, str, bmad_params(20), elegant_params(20)
character(40), allocatable :: names(:)
character(4000) line_out   ! Can be this large for taylor maps.
character(2) continue_char, eol_char, comment_char, separator_char

logical, optional :: use_matrix_model, include_apertures, err
logical init_needed, mad_out, err_flag
logical parsing, warn_printed, converted, ptc_exact_model

! SAD translation

if (out_type == 'SAD') then
  call write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)
  return
endif

! Use ptc exact_model = True since this is needed to get the drift nonlinear terms

ptc_exact_model = ptc_com%exact_model
ptc_com%exact_model = .true.
dr12_max = real_option(1d-5, dr12_drift_max)
if (dr12_max < 0) dr12_max = 1d-5

! Init

ix = integer_option(0, ix_branch)
if (ix < 0 .or. ix > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: /i0/ ', i_array = [ix])
  return
endif

branch => lat%branch(ix)

if (out_type == 'MAD-X' .or. out_type == 'OPAL-T') then
  comment_char = '//'
  continue_char = ''
  eol_char = ';'
  separator_char = ','
  ix_line_max = 100

elseif (out_type == 'MAD-8' .or. out_type == 'XSIF' .or. out_type == 'ELEGANT') then
  comment_char = '!'
  continue_char = ' &'
  eol_char = ''
  separator_char = ','
  ix_line_max = 80

else
  call out_io (s_error$, r_name, 'BAD OUT_TYPE: ' // out_type)
  return
endif

if (out_type == 'ELEGANT') call out_io (s_warn$, r_name, '! NOTE: ELEGANT TRANSLATION IN DEVELOPMENT. PLEASE CONTACT DAVID SAGAN IF YOU WANT TO USE!')

mad_out = .false.
if (out_type == 'MAD-X' .or. out_type == 'MAD-8') mad_out = .true.

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

ie1 = integer_option(1, ix_start)
ie2 = integer_option(branch%n_ele_track, ix_end)

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
ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle  ! Has been marked for delection

  val => ele%value

  ! If there is an aperture with an element that is not an ecoll or rcoll then need to make a separate
  ! element with the aperture info. Exception: MAD-X can handle apertures on non-collimator elements.

  if ((val(x1_limit$) /= 0 .or. val(x2_limit$) /= 0 .or. val(y1_limit$) /= 0 .or. val(y2_limit$) /= 0) .and. &
      ele%key /= ecollimator$ .and. ele%key /= rcollimator$ .and. logic_option(.true., include_apertures) .and. &
      (ele%key == drift$ .or. out_type /= 'MAD-X')) then

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
      ie2 = ie2 + 1
    endif
    if (aperture_at == both_ends$ .or. aperture_at == upstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 1
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
    ie2 = ie2 + 2
    cycle
  endif

  ! If there is a multipole component then put multipole elements at half strength 
  ! just before and just after the element.

  if (ele%key /= multipole$ .and. ele%key /= ab_multipole$ .and. ele%key /= null_ele$ .and. ele%key /= sad_mult$) then
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
      ie2 = ie2 + 2
      cycle
    endif
  endif

  ! If there are nonzero kick values and this is not a kick type element then put
  ! kicker elements at half strength just before and just after the element.
  ! Also add a matrix element to get the change in z correct.
  ! A sad_mult gets translated to a matrix element which has kick components so no extra kickers needed here.
  ! Exception: MAD-X sbend has K0 and K0S attributes.

  if (has_hkick_attributes(ele%key) .and. .not. (ele%key == sbend$  .and. out_type == 'MAD-X')) then
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
        ie2 = ie2 + 3
        cycle
      endif
    endif
  endif

  ! A quadrupole with fringe = full or soft_edge_only has its fringe kicks modeled as a 2nd order map.

  iv = nint(ele%value(fringe_type$))
  if (mad_out .and. ele%key == quadrupole$ .and. (iv == full$ .or. iv == soft_edge_only$)) then
    quad_ele = ele
    ele%value(fringe_type$) = none$

    if (ptc_private%taylor_order_ptc /= 2) call set_ptc (taylor_order = 2) 

    f_count = f_count + 1
    ie = ix_ele

    ifa = nint(ele%value(fringe_at$))
    if (ifa == entrance_end$ .or. ifa == both_ends$) then
      quad_ele%value(fringe_at$) = entrance_end$
      quad_ele%value(l$) = 1d-30
      call ele_to_taylor (quad_ele, branch_out%param, orbit_out(ie-1), orbital_taylor = taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'Q_FRINGE_IN', f_count
      call insert_element (lat_out, taylor_ele, ie, branch_out%ix_branch, orbit_out)
      ie = ie + 1
      ie2 = ie2 + 1
    endif

    if (ifa == exit_end$ .or. ifa == both_ends$) then
      quad_ele%value(fringe_at$) = exit_end$
      quad_ele%value(l$) = 1d-30
      call ele_to_taylor (quad_ele, branch_out%param, orbit_out(ie), orbital_taylor = taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'Q_FRINGE_OUT', f_count
      call insert_element (lat_out, taylor_ele, ie+1, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 1
    endif

    cycle
  endif

  ! A bend with fringe = sad_full or has non-zero dg has its fringe kicks modeled as a 1st order map.

  iv = nint(ele%value(fringe_type$))
  if (ele%key == sbend$ .and. ((mad_out .and. iv == sad_full$) .or. (out_type == 'MAD-8' .and. ele%value(dg$) /= 0))) then

    if (ptc_private%taylor_order_ptc /= 1) call set_ptc (taylor_order = 1)

    f_count = f_count + 1
    ie = ix_ele

    bend_ele = ele
    bend_ele%value(l$) = ele%value(l$)/2
    bend_ele%value(angle$) = ele%value(angle$)/2
    bend_ele%value(e2$) = 0
    call set_fringe_on_off (bend_ele%value(fringe_at$), exit_end$, off$)
    call track1 (orbit_out(ie-1), bend_ele, branch_out%param, orb_center)

    if (at_this_ele_end(entrance_end$, nint(ele%value(fringe_at$))) .or. ele%value(dg$) /= 0) then
      call ele_to_taylor (bend_ele, branch_out%param, orbit_out(ie-1), orbital_taylor = taylor_a)

      bend_ele%value(fringe_type$) = basic_bend$
      bend_ele%value(dg$) = 0
      orb_start = orb_center
      orb_start%direction = -1
      orb_start%species = antiparticle(orb_center%species)
      call track1 (orb_start, bend_ele, branch_out%param, orb_start)  ! bactrack to entrance end
      call ele_to_taylor (bend_ele, branch_out%param, orb_start, orbital_taylor = taylor_b)

      call taylor_inverse (taylor_b, taylor_b)
      call concat_taylor (taylor_a, taylor_b, taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'B_FRINGE_IN', f_count
      call insert_element (lat_out, taylor_ele, ie, branch_out%ix_branch, orbit_out)
      ele => branch_out%ele(ix_ele+1)
      call kill_taylor (taylor_a)
      call kill_taylor (taylor_b)
      ie = ie + 1
      ie2 = ie2 + 1
    endif

    if (at_this_ele_end(exit_end$, nint(ele%value(fringe_at$))) .or. ele%value(dg$) /= 0) then
      bend_ele = ele
      bend_ele%value(l$) = ele%value(l$)/2
      bend_ele%value(angle$) = ele%value(angle$)/2
      bend_ele%value(e1$) = 0
      call set_fringe_on_off (bend_ele%value(fringe_at$), entrance_end$, off$)

      call ele_to_taylor (bend_ele, branch_out%param, orb_center, orbital_taylor = taylor_a)

      bend_ele%value(fringe_type$) = basic_bend$
      bend_ele%value(dg$) = 0
      call ele_to_taylor (bend_ele, branch_out%param, orb_center, orbital_taylor = taylor_b)
      call taylor_inverse (taylor_b, taylor_b)

      call concat_taylor (taylor_b, taylor_a, taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'B_FRINGE_OUT', f_count
      call insert_element (lat_out, taylor_ele, ie+1, branch_out%ix_branch, orbit_out)
      call kill_taylor (taylor_a)
      call kill_taylor (taylor_b)
      ie2 = ie2 + 1
    endif

    ele%value(fringe_type$) = basic_bend$
    ele%value(dg$) = 0
    cycle
  endif

  ! A drift where the ref orbit is too large needs an added 1st order matrix element 

  f = ele%value(l$) / (1 + orbit_out(ele%ix_ele)%vec(6))
  if (mad_out .and. ele%key == drift$ .and. ele%name(1:7) /= 'DRIFT_Z' .and. abs(ele%mat6(1,2) - f) > dr12_max) then
    if (ptc_private%taylor_order_ptc /= 1) call set_ptc (taylor_order = 1) 

    drift_ele = ele
    drift_ele%value(l$) = -ele%value(l$)
    call make_mat6_mad (drift_ele, branch_out%param, orbit_out(ix_ele), orb_end)
    call mat6_to_taylor (drift_ele%vec0, drift_ele%mat6, taylor_a)

    drift_ele%value(l$) = ele%value(l$)
    call ele_to_taylor (drift_ele, branch_out%param, orbit_out(ix_ele-1), orbital_taylor = taylor_b)
    call concat_taylor (taylor_a, taylor_b, taylor_ele%taylor)
    call kill_taylor (taylor_a)
    call kill_taylor (taylor_b)

    taylor_ele%name = 'TAYLOR_' // ele%name
    call insert_element (lat_out, taylor_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
    ie2 = ie2 + 1
    ix_ele = ix_ele + 1
    cycle
  endif

  ! Convert sol_quad_and wiggler elements to an "equivalent" set of elements.
  ! NOTE: FOR NOW, SOL_QUAD USES DRIFT-MATRIX-DRIFT MODEL!

  if (ele%key == wiggler$ .or. ele%key == undulator$ .or. ele%key == sol_quad$) then
    if (logic_option(.false., use_matrix_model) .or. ele%key == sol_quad$) then
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
      ie2 = ie2 + 2
      cycle

    ! Non matrix model...
    ! If the wiggler has been sliced due to superposition, throw 
    ! out the markers that caused the slicing.

    else
      if (ele%key == wiggler$ .or. ele%key == undulator$) then  ! Not a sol_quad
        if (ele%slave_status == super_slave$) then
          ! Create the wiggler model using the super_lord
          lord => pointer_to_lord(ele, 1)
          call out_io (s_warn$, r_name, 'Converting element to drift-bend-drift model: ' // lord%name)
          call create_planar_wiggler_model (lord, lat_model)
          ! Remove all the slave elements and markers in between.
          call out_io (s_warn$, r_name, &
              'Note: Not translating to MAD/XSIF the markers within wiggler: ' // lord%name)
          call find_element_ends (lord, ele1, ele2)
          ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
          lord%key = -1 ! mark for deletion
          ! If the wiggler wraps around the origin we are in trouble.
          if (ix2 < ix1) then 
            call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
            if (global_com%exit_on_error) call err_exit
          endif
          do i = ix1+1, ix2
            branch_out%ele(i)%key = -1  ! mark for deletion
          enddo
          ix_ele = ix_ele + (ix2 - ix1 - 1)
        else
          call out_io (s_warn$, r_name, 'Converting element to drift-bend-drift model: ' // ele%name)
          call create_planar_wiggler_model (ele, lat_model)
          ele%key = -1 ! Mark to ignore
        endif

      else   ! sol_quad
        call create_sol_quad_model (ele, lat_model)  ! NOT YET IMPLEMENTED!
        ele%key = -1 ! Mark to ignore
      endif

      do j = 1, lat_model%n_ele_track
        call insert_element (lat_out, lat_model%ele(j), ix_ele+j, branch_out%ix_branch, orbit_out)
      enddo

      ie2 = ie2 + lat_model%n_ele_track - 1
      cycle
    endif
  endif

enddo

! For a patch that is *not* associated with the edge of a solenoid: A z_offset must be split into a drift + patch

ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle

  ! If the name has more than 16 characters then replace the name by something shorter and unique.

  orig_name = ele%name

  if (len_trim(ele%name) > 16) then
    i_unique = i_unique + 1
    write (ele%name, '(a, i0)') ele%name(1:11), i_unique
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

write (iu, '(3a)') comment_char, ' File generated by: write_lattice_in_foreign_format', trim(eol_char)
write (iu, '(4a)') comment_char, ' Bmad Lattice File: ', trim(lat%input_file_name), trim(eol_char)
if (lat%lattice /= '') write (iu, '(4a)') comment_char, ' Bmad Lattice: ', trim(lat%lattice), trim(eol_char)
write (iu, '(a)')

! beam definition

select case (out_type)
case ('MAD-8', 'MAD-X', 'XSIF')
  ele => branch_out%ele(ie1-1)

  write (line_out, '(7a)') 'beam_def: Beam, Particle = ', trim(species_name(branch_out%param%particle)),  &
        ', Energy = ', re_str(1d-9*ele%value(E_TOT$)), ', Npart = ', re_str(branch_out%param%n_part), trim(eol_char)
  call write_line (line_out)
  write (iu, '(a)')
end select

! write element parameters

n_names = 0                          ! number of names stored in the list
ix_ele = ie1 - 1

do   ! ix_ele = 1e1, ie2
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  if (ele%key == -1) cycle

  val => ele%value

  if (out_type == 'XSIF' .or. out_type == 'ELEGANT') then
    if (ele%key == elseparator$) then 
      n_elsep_warn = n_elsep_warn + 1
      ele%key = drift$  ! XSIF does not have elsep elements.
      call out_io (s_info$, r_name, 'Elseparator being converted into a drift for ' //out_type // ' conversion: ' // ele%name)  
    endif
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

  !------------------------------------
  ! ELEGANT conversion

  if (out_type == 'ELEGANT') then

    bmad_params = ''
    elegant_params = ''

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
      call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, magnetic$, include_kicks$)
      knl = knl / ele%value(l$)

      if (knl(2) == 0) then
        write (line_out, '(2a)') trim(ele%name) // ': kquad'
        if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) line_out = trim(line_out) // ', malign_method = 2'
      else
        write (line_out, '(2a)') trim(ele%name) // ': kquse'
        call value_to_line (line_out, knl(2)*cos(3*(tilts(2)-tilts(1)))/2, 'k2', 'R')
      endif

      tilt = tilts(1)
      call value_to_line (line_out, knl(1), 'k1', 'R')

      bmad_params(:1) = [character(40):: 'l']
      elegant_params(:1) = [character(40):: 'l']

    case (sextupole$)   ! Elegant
      call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, magnetic$, include_kicks$)
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
      call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilts, magnetic$, include_kicks$)
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
      do i = 1, ix_pole_max
        if (knl(i) == 0) cycle
        ab_ele%name = trim(orig_name) // '__' // int_str(i)
        write (line_out, '(2a)') trim(ab_ele%name) // ': mult'
        call insert_element(lat_out, ab_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 1;  ix_ele = ix_ele + 1
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
      write (line_out, '(2a)') trim(ele%name) // ': wiggler'
      bmad_params(:7) = [character(40):: 'l', 'b_max', 'x_offset', 'y_offset', 'z_offset', 'tilt', 'n_pole']
      elegant_params(:7) = [character(40):: 'l', 'b', 'dx', 'dy', 'dz', 'tilt', 'poles']

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
      if (all([ele%value(x_pitch$), ele%value(y_pitch$), ele%value(x_offset$), ele%value(x_offset$), ele%value(x_offset$)] == 0)) then
        write (line_out, '(2a)') trim(ele%name) // ': rotate'
        call value_to_line (line_out, ele%value(tilt$),  'tilt', 'R')
      else
        write (line_out, '(2a)') trim(ele%name) // ': malign'
        bmad_params(:5) = [character(40):: 'x_offset', 'y_offset', 'z_offset', 't_offset', 'e_tot_offset']
        elegant_params(:5) = [character(40):: 'dx', 'dy', 'dz', 'dt', 'de']
        if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
          call out_io (s_warn$, r_name, 'PITCH OR TILT PARAMETERS OF A PATCH CANNOT BE TRANSLATED TO ELEGANT: ' // ele%name)
        endif
      endif

    case (floor_shift$)   ! Elegant
      write (line_out, '(2a)') trim(ele%name) // ': floor'
      call value_to_line (line_out, ele%floor%r(1),  'x', 'R')
      call value_to_line (line_out, ele%floor%r(2),  'y', 'R')
      call value_to_line (line_out, ele%floor%r(3),  'z', 'R')
      call value_to_line (line_out, ele%floor%theta, 'theta', 'R')
      call value_to_line (line_out, ele%floor%phi,   'phi', 'R')
      call value_to_line (line_out, ele%floor%psi,   'psi', 'R')

!    case (match$)   ! Elegant
!      write (line_out, '(2a)') trim(ele%name) // ': ilmatrix'

    case default
      call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO TRANSLATE ELEMENT: ' // ele%name, &
                                     'WHICH IS OF TYPE: ' // key_name(ele%key), &
                                     'CONVERTING TO DRIFT')
      write (line_out, '(2a)') trim(ele%name) // ': drift'
      bmad_params(:1) = [character(40):: 'l']
      elegant_params(:1) = [character(40):: 'l']
    end select

    !------

    select case (ele%key)
    case (sbend$, patch$, drift$)
      ! Pass

    case (quadrupole$, sextupole$, octupole$, taylor$)
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

    case (instrument$, detector$, monitor$, hkicker$, vkicker$, kicker$)  ! Has tilt but not pitches.
      if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) then
        call out_io (s_warn$, r_name, 'X_PITCH OR Y_PITCH PARAMETERS OF A ' // trim(key_name(ele%key)) // ' CANNOT BE TRANSLATED TO ELEGANT: ' // ele%name)
      endif

    case default
      if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
        call out_io (s_warn$, r_name, 'TILT, X_PITCH OR Y_PITCH PARAMETERS OF A ' // trim(key_name(ele%key)) // ' CANNOT BE TRANSLATED TO ELEGANT: ' // ele%name)
      endif
    end select

    do i = 1, size(bmad_params)
      if (bmad_params(i) == '') exit
      call pointer_to_attribute (ele, upcase(bmad_params(i)), .true., a_ptr, err_flag)
      call value_to_line (line_out, a_ptr%r, elegant_params(i), 'R')
    enddo

    call write_line(line_out)
    cycle
  endif

  !------------------------------------
  ! OPAL conversion
  
  if (out_type == 'OPAL-T') then

    select case (ele%key)

    ! OPAL-T
    case (marker$)
      write (line_out, '(a)') trim(ele%name) // ': marker'
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)


    ! OPAL-T
    case (drift$, instrument$, pipe$, detector$, monitor$)
      write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case (sbend$)
      write (line_out, '(2a)') trim(ele%name) // ': sbend, l = ', re_str(val(l$))
      call value_to_line (line_out, val(b_field$), 'k0', 'R')
      call value_to_line (line_out, val(e_tot$), 'designenergy', 'R')
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case (quadrupole$)
      write (line_out, '(2a)') trim(ele%name) // ': quadrupole, l = ', re_str(val(l$))
      !Note that OPAL-T has k1 = dBy/dx, and that bmad needs a -1 sign for electrons
      call value_to_line (line_out, -1*val(b1_gradient$), 'k1', 'R')
      !elemedge The edge of the field is specifieda bsolute (floor space co-ordinates) in m.
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case default
      call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO TRANSLATE ELEMENT: ' // ele%name, &
                                     'WHICH IS OF TYPE: ' // key_name(ele%key), &
                                     'CONVERTING TO DRIFT')
      write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    end select

    call write_line(line_out)
    cycle
  endif

  !-----------------------------------
  ! For anything else but OPAL and ELEGANT

  select case (ele%key)

  ! drift MAD

  case (drift$, instrument$, pipe$, detector$, monitor$)

    write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
  
  ! beambeam MAD

  case (beambeam$)

    line_out = trim(ele%name) // ': beambeam'
    call value_to_line (line_out, val(sig_x$), 'sigx', 'R')
    call value_to_line (line_out, val(sig_y$), 'sigy', 'R')
    call value_to_line (line_out, val(x_offset$), 'xma', 'R')
    call value_to_line (line_out, val(y_offset$), 'yma', 'R')
    call value_to_line (line_out, val(charge$), 'charge', 'R')


  ! r/ecollimator MAD

  case (ecollimator$, rcollimator$)

    if (out_type == 'MAD-X') then
      write (line_out, '(2a)') trim(ele%name) // ': collimator, l = ', re_str(val(l$))
    else
      write (line_out, '(2a)') trim(ele%name) // ': ' // trim(key_name(ele%key)) // ', l = ', re_str(val(l$))
      call value_to_line (line_out, val(x1_limit$), 'xsize', 'R')
      call value_to_line (line_out, val(y1_limit$), 'ysize', 'R')
    endif

  ! elseparator MAD

  case (elseparator$)

    write (line_out, '(2a)') trim(ele%name) // ': elseparator, l = ', re_str(val(l$))
    hk = val(hkick$)
    vk = val(vkick$)

    if (hk /= 0 .or. vk /= 0) then

      ix = len_trim(line_out) + 1
      field = 1.0d3 * sqrt(hk**2 + vk**2) * val(E_TOT$) / val(l$)
      if (out_type == 'MAD-X') then
        write (line_out(ix:), '(2a)') ', ey = ', re_str(field)
      else
        write (line_out(ix:), '(2a)') ', e = ',re_str(field)
      endif

      if (branch_out%param%particle == positron$) then
        tilt = -atan2(hk, vk) + val(tilt$)
      else
        tilt = -atan2(hk, vk) + val(tilt$) + pi
      endif
      ix = len_trim(line_out) + 1
      write (line_out(ix:), '(2a)') ', tilt = ', re_str(tilt)

    endif

  ! hkicker MAD

  case (hkicker$)

    write (line_out, '(2a)') trim(ele%name) // ': hkicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(kick$), 'kick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! kicker MAD

  case (kicker$)

    write (line_out, '(2a)') trim(ele%name) // ': kicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(hkick$), 'hkick', 'R')
    call value_to_line (line_out, val(vkick$), 'vkick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! vkicker MAD

  case (vkicker$)

    write (line_out, '(2a)') trim(ele%name) // ': vkicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(kick$), 'kick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! marker MAD

  case (marker$, fork$, photon_fork$)

    line_out = trim(ele%name) // ': marker'

  ! octupole MAD

  case (octupole$)

    write (line_out, '(2a)') trim(ele%name) // ': octupole, l = ', re_str(val(l$))

    call value_to_line (line_out, val(k3$), 'k3', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! quadrupole MAD

  case (quadrupole$)

    write (line_out, '(2a)') trim(ele%name) // ': quadrupole, l = ', re_str(val(l$))
    call value_to_line (line_out, val(k1$), 'k1', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! sbend MAD

  case (sbend$)

    write (line_out, '(2a)') trim(ele%name) // ': sbend, l = ', re_str(val(l$))

    call value_to_line (line_out, val(angle$), 'angle', 'R')
    call value_to_line (line_out, val(e1$), 'e1', 'R')
    call value_to_line (line_out, val(e2$), 'e2', 'R')
    call value_to_line (line_out, val(k1$), 'k1', 'R')
    call value_to_line (line_out, val(ref_tilt$), 'tilt', 'R')
    if (out_type == 'MAD-X') then
      call value_to_line (line_out, val(fint$), 'fint', 'R')
      call value_to_line (line_out, val(fintx$), 'fintx', 'R')
      call value_to_line (line_out, val(hgap$), 'hgap', 'R')
    else
      if (val(fintx$) /= val(fint$)) then
        call out_io (s_info$, r_name, 'FINTX != FINT FOR BEND' // ele%name, 'CANNOT TRANSLATE FINTX')
      endif
      call value_to_line (line_out, val(fint$), 'fint', 'R')
      call value_to_line (line_out, val(hgap$), 'hgap', 'R')
    endif

    ! MAD-X sbend kick fields. MAD-8 conversion uses matrix elements to either side (see above).

    if (out_type == 'MAD-X' .and. ele%value(l$) /= 0) then
      call multipole_ele_to_ab (ele, .false., ix, a_pole, b_pole, magnetic$, include_kicks$)
      call value_to_line (line_out, val(dg$) + b_pole(0)/val(l$), 'k0', 'R')
      call value_to_line (line_out, a_pole(0)/val(l$), 'k0s', 'R')
    endif

  ! sextupole MAD

  case (sextupole$)

    write (line_out, '(2a)') trim(ele%name) // ': sextupole, l = ', re_str(val(l$))
    call value_to_line (line_out, val(k2$), 'k2', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! taylor MAD

  case (taylor$, sad_mult$, patch$, match$)

    if (ele%key == patch$ .and. out_type == 'MAD-X') then
      ele%key = null_ele$
      orig_name = ele%name
      if (val(x_offset$) /= 0 .or. val(y_offset$) /= 0 .or. val(z_offset$) /= 0) then
        drift_ele%name = trim(orig_name) // '__t'
        call insert_element(lat_out, drift_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 1;  ix_ele = ix_ele + 1
        line_out = trim(drift_ele%name) // ': translation'
        call value_to_line (line_out, val(x_offset$), 'dx', 'R')
        call value_to_line (line_out, val(y_offset$), 'dy', 'R')
        call value_to_line (line_out, val(z_offset$), 'ds', 'R')
        call write_line(line_out)
      endif

      if (val(x_pitch$) /= 0) then
        drift_ele%name = trim(orig_name) // '__y'
        call insert_element(lat_out, drift_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 1;  ix_ele = ix_ele + 1
        call write_line(trim(drift_ele%name) // ': yrotation, angle = ' // re_str(-val(x_pitch$)))
      endif

      if (val(y_pitch$) /= 0) then
        drift_ele%name = trim(orig_name) // '__x'
        call insert_element(lat_out, drift_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 1;  ix_ele = ix_ele + 1
        call write_line(trim(drift_ele%name) // ': xrotation, angle = ' // re_str(-val(y_pitch$)))
      endif

      if (val(tilt$) /= 0) then
        drift_ele%name = trim(orig_name) // '__s'
        call insert_element(lat_out, drift_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 1;  ix_ele = ix_ele + 1
        call write_line(trim(drift_ele%name) // ': srotation, angle = ' // re_str(val(tilt$)))
      endif

      cycle
    endif

    if (associated (ele%taylor(1)%term)) then
      taylor_ptr => ele%taylor
    elseif (ele%key == match$) then
      allocate(taylor_ptr(6))
      call ele_to_taylor (ele, branch%param, orbital_taylor = taylor_ptr)
    else
      allocate(taylor_ptr(6))
      if (.not. present(ref_orbit)) then
        call out_io (s_error$, r_name, &
                      'ORBIT ARGUMENT NEEDS TO BE PRESENT WHEN TRANSLATING', &
                      'A LATTICE WITH A SAD_MULT OR PATCH ELEMENT')           
        cycle
      endif
      if (ptc_private%taylor_order_ptc /= 2) call set_ptc (taylor_order = 2) 
      call ele_to_taylor (ele, branch%param, orbit_out(ix_ele-1), .true., orbital_taylor = taylor_ptr)
    endif

    line_out = trim(ele%name) // ': matrix'
    warn_printed = .false.
    call value_to_line (line_out, val(l$), 'l', 'R')

    do i = 1, 6
      do k = 1, size(taylor_ptr(i)%term)
        term = taylor_ptr(i)%term(k)

        select case (sum(term%expn))
        case (0)
          select case (out_type)
          case ('MAD-8') 
            write (str, '(a, i0, a)') 'kick(', i, ')'
          case ('MAD-X') 
            write (str, '(a, i0)') 'kick', i
          case ('XSIF') 
            call out_io (s_error$, r_name, 'XSIF DOES NOT HAVE A CONSTRUCT FOR ZEROTH ORDER TAYLOR TERMS NEEDED FOR: ' // ele%name)
            cycle
          end select
          call value_to_line (line_out, term%coef, str, 'R')

        case (1)
          j = maxloc(term%expn, 1)
          select case (out_type)
          case ('MAD-8')
            write (str, '(a, i0, a, i0, a)') 'rm(', i, ',', j, ')'
          case ('MAD-X')
            write (str, '(a, 2i0)') 'rm', i, j
          case ('XSIF')
            write (str, '(a, 2i0)') 'r', i, j
          end select

          if (j == i) then
            call value_to_line (line_out, term%coef, str, 'R', .false.)
          else
            call value_to_line (line_out, term%coef, str, 'R')
          endif

        case (2)
          j = maxloc(term%expn, 1)
          term%expn(j) = term%expn(j) - 1
          j2 = maxloc(term%expn, 1)
          select case (out_type)
          case ('MAD-8')
            write (str, '(a, 3(i0, a))') 'tm(', i, ',', j, ',', j2, ')'
          case ('MAD-X')
            write (str, '(a, 3i0)') 'tm', i, j, j2
          case ('XSIF')
            write (str, '(a, 3i0)') 't', i, j, j2
          end select
          call value_to_line (line_out, term%coef, str, 'R')

        case default
          if (.not. warn_printed .and. ele%key == taylor$) then
            call out_io (s_warn$, r_name, &
                  'Higher order taylor term(s) in: ' // trim(ele%name) // &
                  ' cannot be converted to mad matrix term')
            warn_printed = .true.
          endif  
        end select
      enddo

    enddo

    if (.not. associated(ele%taylor(1)%term)) deallocate(taylor_ptr)

  ! rfcavity MAD

  case (rfcavity$)

    write (line_out, '(2a)') trim(ele%name) // ': rfcavity, l = ', re_str(val(l$))
    call value_to_line (line_out, val(voltage$)/1E6, 'volt', 'R')
    call value_to_line (line_out, val(phi0$)+val(phi0_multipass$)+0.5, 'lag', 'R')
    call value_to_line (line_out, val(harmon$), 'harmon', 'I')

  ! lcavity MAD

  case (lcavity$)

    write (line_out, '(2a)') trim(ele%name) // ': lcavity, l = ', re_str(val(l$))
    call value_to_line (line_out, val(gradient$)*val(l$)/1d6, 'deltae', 'R')
    call value_to_line (line_out, val(rf_frequency$)/1d6, 'freq', 'R')
    call value_to_line (line_out, val(phi0$)+val(phi0_multipass$), 'phi0', 'R')

  ! solenoid MAD

  case (solenoid$)

    write (line_out, '(2a)') trim(ele%name) // ': solenoid, l = ', re_str(val(l$))
    call value_to_line (line_out, val(ks$), 'ks', 'R')

  ! multipole MAD

  case (multipole$, ab_multipole$)

    knl = 0; tilts = 0
    call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilts)
    write (line_out, '(2a)') trim(ele%name) // ': multipole'  

    if (out_type == 'MAD-X') then
      knl_str = ''; ksl_str = ''
      call multipole_ele_to_ab (ele, .true., ix_pole_max, a_pole, b_pole)
      do i = 0, 9
        if (all(knl(i:) == 0)) exit
        if (abs(a_pole(i)) < 1d-12 * abs(b_pole(i))) a_pole(i) = 0  ! Round to zero insignificant value
        if (abs(b_pole(i)) < 1d-12 * abs(a_pole(i))) b_pole(i) = 0  ! Round to zero insignificant value
        call value_to_line (knl_str,  b_pole(i) * factorial(i), '', 'R', .false.)
        call value_to_line (ksl_str, -a_pole(i) * factorial(i), '', 'R', .false.)
      enddo
      if (any(b_pole /= 0)) line_out = trim(line_out) // ', knl = {' // trim(knl_str(3:)) // '}'
      if (any(a_pole /= 0)) line_out = trim(line_out) // ', ksl = {' // trim(ksl_str(3:)) // '}'

    else
      do i = 0, 9
        write (str, '(a, i0, a)') 'K', i, 'L'
        call value_to_line (line_out, knl(i), str, 'R')
        write (str, '(a, i0)') 'T', i
        call value_to_line (line_out, tilts(i), str, 'R')
      enddo
    endif

  ! unknown MAD

  case default

    call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO TRANSLATE ELEMENT: ' // ele%name, &
                                   'WHICH IS OF TYPE: ' // key_name(ele%key), &
                                   'CONVERTING TO DRIFT')
    line_out = trim(ele%name) // ': drift, l = ' // re_str(val(l$))

  end select

  ! Add apertures for mad-x. Use 1 meter for unset apertures

  if (out_type == 'MAD-X' .and. logic_option(.true., include_apertures)) then
    if (val(x1_limit$) /= 0 .or. val(y1_limit$) /= 0) then
      limit = [val(x1_limit$), val(y1_limit$)]
      where (limit == 0) limit = 1
      if (ele%aperture_type == rectangular$) then
        line_out = trim(line_out) // ', apertype = rectangle'
      else
        line_out = trim(line_out) // ', apertype = ellipse'
      endif
      write (line_out, '(6a)') trim(line_out), ', aperture = {', re_str(limit(1)), ', ', re_str(limit(2)), '}'
    endif
  endif

  ! write element spec to file

  call write_line(line_out)

enddo

!---------------------------------------------------------------------------------------
! Write the lattice line
! MAD has a limit of 4000 characters so we may need to break the lat into pieces.

i_unique = 1000
i_line = 0
init_needed = .true.
line = ' '

do n = ie1, ie2
  ele => branch_out%ele(n)
  if (ele%key == null_ele$) cycle  ! Will happen with patch elements translated to MAD-X
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

  ! Output line if long enough or at end

  if (n == ie2 .or. iout > 48) then
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

if (out_type == 'MAD-X') then
  write (iu, '(a)') 'use, period = lat;'
elseif (out_type /= 'OPAL-T') then
  write (iu, '(a)') 'use, lat'
endif

!---------------------------------------------------
! Element offsets for MAD.
! This must come after use statement.

if (out_type(1:3) == 'MAD') then

  write (iu, '(a)')
  write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
  write (iu, '(a)')

  allocate (n_repeat(n_names))
  n_repeat = 0

  do ix_ele = ie1, ie2

    ele => branch_out%ele(ix_ele)
    val => ele%value

    ! sad_mult and patch elements are translated to a matrix which does not have offsets.
    ! And marker like elements also do not have offsets

    if (ele%key == sad_mult$ .or. ele%key == patch$) cycle
    if (ele%key == marker$ .or. ele%key == fork$ .or. ele%key == photon_fork$) cycle

    !

    call find_index (ele%name, names, an_indexx, n_names, ix_match)
    if (ix_match == 0) cycle ! Happens for translated to MADX patch elements.
    n_repeat(ix_match) = n_repeat(ix_match) + 1
    
    if (val(x_pitch$) == 0 .and. val(y_pitch$) == 0 .and. &
        val(x_offset_tot$) == 0 .and. val(y_offset_tot$) == 0 .and. val(z_offset_tot$) == 0) cycle

    write (iu, '(3a)') 'select, flag = error, clear', trim(eol_char)
    write (iu, '(3a, i0, 2a)') 'select, flag = error, range = ', trim(ele%name), &
                                    '[', n_repeat(ix_match), ']', trim(eol_char)

    line_out = 'ealign'
    call value_to_line (line_out,  val(x_pitch$), 'dtheta', 'R')
    call value_to_line (line_out, -val(y_pitch$), 'dphi', 'R')
    call value_to_line (line_out, val(x_offset$) - val(x_pitch$) * val(l$) / 2, 'dx', 'R')
    call value_to_line (line_out, val(y_offset$) - val(y_pitch$) * val(l$) / 2, 'dy', 'R')
    call value_to_line (line_out, val(z_offset$), 'ds', 'R')
    call write_line (line_out)

  enddo

  deallocate (n_repeat)

endif

! Write twiss parameters for a non-closed lattice.

if (branch_out%param%geometry == open$ .and. (out_type == 'MAD-8' .or. out_type == 'MAD-X' .or. out_type == 'XSIF')) then
  ele => branch_out%ele(ie1-1)
  orb_start = lat%particle_start
  beta = ele%value(p0c$) / ele%value(E_tot$)
  write (iu, '(a)')
  write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
  write (iu, '(a)')
  write (iu, '(12a)') 'initial: beta0, betx = ', re_str(ele%a%beta), ', bety = ', re_str(ele%b%beta), &
                      ', alfx = ', re_str(ele%a%alpha), ', alfy = ', re_str(ele%b%alpha), ', ', trim(continue_char)
  write (iu, '(5x, 12a)') 'dx = ', re_str(ele%a%eta), ', dpx = ', re_str(ele%a%etap), & 
                        ', dy = ', re_str(ele%b%eta), ', dpy = ', re_str(ele%b%etap), ', ', trim(continue_char)
  write (iu, '(5x, 12a)') 'x = ', re_str(orb_start%vec(1)), ', px = ', re_str(orb_start%vec(2)), &
                        ', y = ', re_str(orb_start%vec(3)), ', py = ', re_str(orb_start%vec(4)), &
                        ', t = ', re_str(orb_start%vec(5)*beta), ', pt = ', re_str(orb_start%vec(6)/beta), trim(eol_char)



  if (ele%a%beta /= 0 .and. ele%b%beta /= 0) then
    write (iu, '(a)') 'twiss, beta0 = initial;'
  endif
endif

! End stuff

call out_io (s_info$, r_name, 'Written ' // trim(out_type) // ' lattice file: ' // trim(out_file_name))

deallocate (names)
if (present(err)) err = .false.

if (present(converted_lat)) then
  converted_lat = lat
  converted_lat%branch(branch%ix_branch) = branch_out
  converted_lat%n_ele_max = converted_lat%n_ele_track
  do ib = 0, ubound(converted_lat%branch, 1)
    branch => converted_lat%branch(ib)
    do i = 1, branch%n_ele_track
      branch%ele(i)%slave_status = free$
      branch%ele(i)%n_lord = 0
    enddo
  enddo
  converted_lat%n_control_max = 0
  converted_lat%n_ic_max = 0
endif

call deallocate_lat_pointers (lat_out)
call deallocate_lat_pointers (lat_model)

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

end subroutine write_lattice_in_foreign_format
