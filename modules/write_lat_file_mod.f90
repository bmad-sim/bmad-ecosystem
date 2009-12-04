module write_lat_file_mod

use bmad_struct
use bmad_interface
use multipole_mod
use multipass_mod
use element_modeling_mod
use lat_ele_loc_mod

private str, rchomp, write_out

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_bmad_lattice_file (bmad_file, lat, err)
!
! Subroutine to write a Bmad lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
!
! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   bmad_file     -- Character(*): Name of the output lattice file.
!   lat           -- lat_struct: Holds the lattice information.
!   ix_start      -- Integer, optional: Starting index of lat%ele(i)
!                       used for output.
!   ix_end        -- Integer, optional: Ending index of lat%ele(i)
!                       used for output.
!
! Output:
!   err    -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_bmad_lattice_file (bmad_file, lat, err)

implicit none

type multipass_info_struct
  integer ix_region
  logical region_start_pt
  logical region_stop_pt
end type

type (multipass_info_struct), allocatable :: multipass(:)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, super, slave, lord, s1, s2, multi_lord
type (ele_struct), save :: ele_init
type (wig_term_struct) wt
type (control_struct) ctl
type (taylor_term_struct) tm
type (multipass_all_info_struct), target :: m_info
type (lr_wake_struct), pointer :: lr
type (ele_pointer_struct), pointer :: ss1(:), ss2(:)

real(rp) s0, x_lim, y_lim, val

character(*) bmad_file
character(4000) line
character(4) last
character(40) name, look_for, attrib_name
character(200) wake_name, file_name
character(40), allocatable :: names(:)
character(200), allocatable, save :: sr_wake_name(:), lr_wake_name(:)
character(40) :: r_name = 'write_bmad_lattice_file'
character(10) angle

integer j, k, n, ix, iu, iuw, ios, ixs, n_sr, n_lr, ix1, ie, ib
integer unit(6), ix_names, ix_match
integer ix_slave, ix_ss, ix_l, ix_r, ix_pass
integer ix_top, ix_super, default_val

logical, optional :: err
logical unit_found, write_term, match_found, found, in_multi_region, expand_lat_out
logical is_multi_sup, x_lim_good, y_lim_good, is_default

! Init...
! Count the number of foreign wake files

if (present(err)) err = .true.
call init_ele (ele_init)

n_sr = 0
n_lr = 0
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (.not. associated(ele%wake)) cycle
  if (ele%wake%sr_file(1:6) == 'xsif::') n_sr = n_sr + 1 
  if (ele%wake%lr_file(1:6) == 'xsif::') n_lr = n_lr + 1  
enddo
call re_allocate(sr_wake_name, n_sr)
call re_allocate(lr_wake_name, n_lr)

n_sr = 0
n_lr = 0

! Open the file

iu = lunget()
call fullfilename (bmad_file, file_name)
open (iu, file = file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(bmad_file))
  return
endif

! Non-elemental stuff

if (lat%title /= ' ') &
          write (iu, *) 'title, "', trim(lat%title), '"'
if (lat%lattice /= ' ') &
          write (iu, *) 'parameter[lattice] = "', trim(lat%lattice), '"'
write (iu, *) 'parameter[lattice_type] = ', lattice_type(lat%param%lattice_type)
write (iu, *) 'parameter[taylor_order] =', lat%input_taylor_order

write (iu, *)
write (iu, *) 'parameter[e_tot] =', &
                    trim(str(lat%ele(0)%value(e_tot$)))
write (iu, *) 'parameter[particle] = ', particle_name(lat%param%particle)
if (lat%param%n_part /= 0) &
        write (iu, *) 'parameter[n_part] = ', lat%param%n_part

ele => lat%ele(0) 

if (ele%floor%x /= 0) &
      write (iu, *) 'beginning[x_position] = ', trim(str(ele%floor%x))
if (ele%floor%y /= 0) &
      write (iu, *) 'beginning[y_position] = ', trim(str(ele%floor%y))
if (ele%floor%z /= 0) &
      write (iu, *) 'beginning[z_position] = ', trim(str(ele%floor%z))
if (ele%floor%theta /= 0) &
  write (iu, *) 'beginning[theta_position] = ', trim(str(ele%floor%theta))
if (ele%floor%phi /= 0) &
      write (iu, *) 'beginning[phi_position] = ', trim(str(ele%floor%phi))
if (ele%floor%psi /= 0) &
      write (iu, *) 'beginning[psi_position] = ', trim(str(ele%floor%psi))

if (lat%param%lattice_type /= circular_lattice$) then
  write (iu, *)
  if (ele%a%beta /= 0) &
      write (iu, *) 'beginning[beta_a] = ', trim(str(ele%a%beta))
  if (ele%a%alpha /= 0) &
      write (iu, *) 'beginning[alpha_a] = ', trim(str(ele%a%alpha))
  if (ele%a%phi /= 0) &
      write (iu, *) 'beginning[phi_a] = ', trim(str(ele%a%phi))
  if (ele%a%eta /= 0) &
      write (iu, *) 'beginning[eta_a] = ', trim(str(ele%a%eta))
  if (ele%a%etap /= 0) &
      write (iu, *) 'beginning[etap_a] = ', trim(str(ele%a%etap))
  if (ele%b%beta /= 0) &
      write (iu, *) 'beginning[beta_b] = ', trim(str(ele%b%beta))
  if (ele%b%alpha /= 0) &
      write (iu, *) 'beginning[alpha_b] = ', trim(str(ele%b%alpha))
  if (ele%b%phi /= 0) &
      write (iu, *) 'beginning[phi_b] = ', trim(str(ele%b%phi))
  if (ele%b%eta /= 0) &
      write (iu, *) 'beginning[eta_b] = ', trim(str(ele%b%eta))
  if (ele%b%etap /= 0) &
      write (iu, *) 'beginning[etap_b] = ', trim(str(ele%b%etap))
  if (ele%c_mat(1,1) /= 0) &
      write (iu, *) 'beginning[c11] = ', trim(str(ele%c_mat(1,1)))
  if (ele%c_mat(1,2) /= 0) &
      write (iu, *) 'beginning[c12] = ', trim(str(ele%c_mat(1,2)))
  if (ele%c_mat(2,1) /= 0) &
      write (iu, *) 'beginning[c21] = ', trim(str(ele%c_mat(2,1)))
  if (ele%c_mat(2,2) /= 0) &
      write (iu, *) 'beginning[c22] = ', trim(str(ele%c_mat(2,2)))
endif

! Element stuff

write (iu, *)
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, *)

ixs = 0
ix_names = 0
allocate (names(lat%n_ele_max))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 1, branch%n_ele_max

    ele => branch%ele(ie)

    multi_lord = pointer_to_multipass_lord (ele, lat, ix_pass) 

    if (ele%key == null_ele$) cycle
    if (ele%slave_status == multipass_slave$) cycle ! Ignore for now
    if (ele%lord_status == super_lord$ .and. multi_lord%ix_ele > 0) cycle
    if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle

    if (ie == lat%n_ele_track+1) then
      write (iu, *)
      write (iu, '(a)') '!-------------------------------------------------------'
      write (iu, '(a)') '! Overlays, groups, etc.'
      write (iu, *)
    endif

    ! For a super_slave just create a dummy drift. 

    if (ele%slave_status == super_slave$) then
      ixs = ixs + 1
      ele%ixx = ixs
      write (iu, '(a, i3.3, 2a)') 'slave_drift_', ixs, &
                                        ': drift, l = ', trim(str(ele%value(l$)))
      cycle
    endif

    ! Do not write anything for elements that have a duplicate name.

    call find1_indexx (ele%name, names, ix_names, ix_match, match_found)
    if (match_found) cycle

    names(ix_match+1:ix_names+1) = names(ix_match:ix_names)
    names(ix_match) = ele%name
    ix_names = ix_names + 1

    ! Overlays and groups

    if (ele%lord_status == overlay_lord$ .or. ele%lord_status == group_lord$) then
      if (ele%lord_status == overlay_lord$) then
        write (line, '(2a)') trim(ele%name), ': overlay = {'
      else
        write (line, '(2a)') trim(ele%name), ': group = {'
      endif
      j_loop: do j = ele%ix1_slave, ele%ix2_slave
        ctl = lat%control(j)
        slave => pointer_to_ele (lat, ctl%ix_slave, ctl%ix_branch)
        do k = ele%ix1_slave, j-1 ! do not use elements w/ duplicate names
          if (lat%ele(lat%control(k)%ix_slave)%name == slave%name) cycle j_loop
        enddo
        if (j == ele%ix1_slave) then
          write (line, '(3a)') trim(line), trim(slave%name)
        else
          write (line, '(3a)') trim(line), ', ', trim(slave%name)
        endif
        name = attribute_name(slave, ctl%ix_attrib)  
        if (name /= ele%component_name) &
                line = trim(line) // '[' // trim(name) // ']'
        if (ctl%coef /= 1) write (line, '(3a)') trim(line), '/', trim(str(ctl%coef))
      enddo j_loop
      line = trim(line) // '}'
      if (ele%component_name == ' ') then
        line = trim(line) // ', command'
      else
        line = trim(line) // ', ' // ele%component_name
      endif
      if (ele%lord_status == overlay_lord$) then
        ix = ele%ix_value
        if (ele%value(ix) /= 0) write (line, '(3a)') &
                            trim(line), ' = ', str(ele%value(ix))
      endif
      if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
      if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
      if (associated(ele%descrip)) line = trim(line) // &
                              ', descrip = "' // trim(ele%descrip) // '"'
      call write_out (line, iu, .true.)
      cycle
    endif

    ! Girder

    if (ele%lord_status == girder_lord$) then
      write (line, '(2a)') trim(ele%name), ': girder = {'
      do j = ele%ix1_slave, ele%ix2_slave
        ix1 = lat%control(j)%ix_slave
        if (j == ele%ix2_slave) then
          write (line, '(3a)') trim(line), trim(lat%ele(ix1)%name), '}'
        else
          write (line, '(3a)') trim(line), trim(lat%ele(ix1)%name), ', '
        endif
      enddo
    else
      line = trim(ele%name) // ': ' // key_name(ele%key)
    endif

    ! Branch

    if (ele%key == branch$ .or. ele%key == photon_branch$) then
      n = nint(ele%value(ix_branch_to$))
      line = trim(line) // ', to = ' // trim(lat%branch(n)%name)
    endif

    ! Other elements

    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
    if (associated(ele%descrip)) line = trim(line) // &
                              ', descrip = "' // trim(ele%descrip) // '"'

    ! patch_in_slave

    if (ele%slave_status == patch_in_slave$) then
      lord => pointer_to_lord (lat, ele, 1)
      line = trim(line) // ', ref_patch = ' // trim(lord%name)
    endif

    ! Create a null_ele element for a superposition and fill in the superposition
    ! information.

    is_multi_sup = .false.
    if (ele%lord_status == multipass_lord$) then
      ix1 = lat%control(ele%ix1_slave)%ix_slave
      if (lat%ele(ix1)%lord_status == super_lord$) is_multi_sup = .true.
    endif

    if (ele%lord_status == super_lord$ .or. is_multi_sup) then
      write (iu, '(a)') "x__" // trim(ele%name) // ": null_ele"
      line = trim(line) // ', superimpose, ele_beginning, ref = x__' // trim(ele%name)
    endif

    ! If the wake file is not BMAD Format (Eg: XSIF format) then create a new wake file.
    ! If first three characters of the file name are '...' then it is a foreign file.

    if (associated(ele%wake)) then

      ! Short-range

      if (ele%wake%sr_file /= ' ') then

        wake_name = ele%wake%sr_file

        if (wake_name(1:3) == '...') then
          found = .false.
          do n = 1, n_sr
            if (wake_name == sr_wake_name(n)) then
              found = .true.
              exit
            endif
          enddo
          if (.not. found) then
            n = n_sr + 1
            n_sr = n
            sr_wake_name(n_sr) = wake_name
          endif
          write (wake_name, '(a, i0, a)') 'sr_wake_file_', n, '.bmad'
          if (.not. found) then
            call out_io (s_info$, r_name, 'Creating SR Wake file: ' // trim(wake_name))
            iuw = lunget()
            open (iuw, file = wake_name)
            write (iuw, *) '!      z           Wz             Wt'
            write (iuw, *) '!     [m]       [V/C/m]       [V/C/m^2]'
            do n = lbound(ele%wake%sr_table, 1), ubound(ele%wake%sr_table, 1)
              write (iuw, '(3es14.5)') ele%wake%sr_table(n)%z, &
                                    ele%wake%sr_table(n)%long, ele%wake%sr_table(n)%trans
            enddo
            close(iuw)
          endif
        endif

        line = trim(line) // ',  sr_wake_file = "' // trim(wake_name) // '"'

      endif

      ! Long-range

      if (ele%wake%lr_file /= ' ') then

        wake_name = ele%wake%lr_file

        if (wake_name(1:3) == '...') then
          found = .false.
          do n = 1, n_lr
            if (wake_name == lr_wake_name(n)) then
              found = .true.
              exit
            endif
          enddo
          if (.not. found) then
            n = n_lr + 1
            n_lr = n
            lr_wake_name(n_lr) = wake_name
          endif
          write (wake_name, '(a, i0, a)') 'lr_wake_file_', n, '.bmad'
          if (.not. found) then
            call out_io (s_info$, r_name, 'Creating LR Wake file: ' // trim(wake_name))
            iuw = lunget()
            open (iuw, file = wake_name)
            write (iuw, '(14x, a)') &
   'Freq       R/Q      Q       m  Polarization   b_sin       b_cos       a_sin       a_cos       t_ref'
            write (iuw, '(14x, a)') &
              '[Hz]  [Ohm/m^(2m)]             [Rad/2pi]'
            do n = lbound(ele%wake%lr, 1), ubound(ele%wake%lr, 1)
              lr => ele%wake%lr(n)
              if (lr%polarized) then
                write (angle, '(f10.6)') lr%angle
              else
                angle = '     unpol'
              endif
              if (any ( (/ lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref /) /= 0)) then
                write (iuw, '(a, i0, a, 3es14.5, i6, a, 5es12.2)') 'lr(', n, ') =', lr%freq_in, &
                      lr%R_over_Q, lr%Q, lr%m, angle, lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref
              else
                write (iuw, '(a, i0, a, 3es14.5, i6, a)') 'lr(', n, ') =', &
                      lr%freq_in, lr%R_over_Q, lr%Q, lr%m, angle
              endif
            enddo
            close(iuw)
          endif
        endif

        line = trim(line) // ',  lr_wake_file = "' // trim(wake_name) // '"'

      endif

    endif

    ! Decide if x1_limit, etc. are to be output directly or combined. 

    x_lim = ele%value(x1_limit$) 
    x_lim_good = .false.
    if (x_lim /=0 .and. ele%value(x2_limit$) == x_lim) x_lim_good = .true.

    y_lim = ele%value(y1_limit$) 
    y_lim_good = .false.
    if (y_lim /=0 .and. ele%value(y2_limit$) == y_lim) y_lim_good = .true.

    ! Print the element attributes.

    do j = 1, n_attrib_maxx
      attrib_name = attribute_name(ele, j)
      val = ele%value(j)
      if (val == 0) cycle
      if (j == check_sum$) cycle
      if (x_lim_good .and. (j == x1_limit$ .or. j == x2_limit$)) cycle
      if (y_lim_good .and. (j == y1_limit$ .or. j == y2_limit$)) cycle
      if (.not. attribute_free (ele, attrib_name, lat, .false., .true.)) cycle
      if (attrib_name == 'DS_STEP' .and. val == bmad_com%default_ds_step) cycle
      if (attrib_name == null_name$) then
        print *, 'ERROR IN WRITE_BMAD_LATTICE_FILE:'
        print *, '      ELEMENT: ', ele%name
        print *, '      HAS AN UNKNOWN ATTRIBUTE INDEX:', j
        stop
      endif

      if (attrib_name == 'COUPLER_AT') then
        if (nint(val) /= exit_end$) then
          line = trim(line) // ', coupler_at = ' // element_end_name(nint(val))
        endif
        cycle
      endif

      select case (attribute_type(attrib_name))
      case (is_logical$)
        write (line, '(4a, l1)') trim(line), ', ', trim(attrib_name), ' = ', (val /= 0)
      case (is_integer$)
        write (line, '(4a, i0)') trim(line), ', ', trim(attrib_name), ' = ', int(val)
      case (is_real$)
        line = trim(line) // ', ' // trim(attrib_name) // ' = ' // str(val)
      case (is_name$)
        name = attribute_value_name (attrib_name, val, ele, is_default)
          if (.not. is_default) then
            line = trim(line) // ', ' // trim(attrib_name) // ' = ' // name
          endif
      end select

    enddo ! attribute loop

    ! Print the combined limits if needed.

    if (x_lim_good .and. y_lim_good .and. x_lim == y_lim) then
      line = trim(line) // ', aperture = ' // str(x_lim)
    else
      if (x_lim_good) line = trim(line) // ', x_limit = ' // str(x_lim)
      if (y_lim_good) line = trim(line) // ', y_limit = ' // str(y_lim)
    endif

    ! Encode methods, etc.

    if (ele%mat6_calc_method /= bmad_standard$) line = trim(line) // &
                ', mat6_calc_method = ' // calc_method_name(ele%mat6_calc_method)
    if (ele%tracking_method /= bmad_standard$) line = trim(line) // &
                ', tracking_method = ' // calc_method_name(ele%tracking_method)
    if (ele%field_calc /= bmad_standard$) line = trim(line) // &
                ', field_calc = ' // calc_method_name(ele%field_calc)
    if (ele%symplectify) line = trim(line) // ', symplectify'
    if (attribute_index(ele, 'FIELD_MASTER') /= 0 .and. ele%field_master) &
                line = trim(line) // ', field_master = True'
    if (.not. ele%is_on) line = trim(line) // ', is_on = False'
    if (.not. ele%map_with_offsets) line = trim(line) // ', map_with_offsets = False'
    if (.not. ele%csr_calc_on) line = trim(line) // ', csr_calc_on = False'
    if (ele%offset_moves_aperture) line = trim(line) // ', offset_moves_aperture = True'
    if (ele%aperture_at /= exit_end$) line = trim(line) // ', aperture_at = ' // & 
                                                         element_end_name(ele%aperture_at)

    default_val = rectangular$
    if (ele%key == ecollimator$) default_val = elliptical$
    if (ele%aperture_type /= default_val) line = trim(line) // &
                                ', aperture_type = ' // shape_name(ele%aperture_type)

    if (ele%integrator_order /= bmad_com%default_integ_order) write (line, '(a, i0)') &
                   trim(line) // ', integrator_order = ', ele%integrator_order
    if (ele%ref_orbit /= 0) line = trim(line) // ', ref_orbit = ' // ref_orbit_name(ele%ref_orbit)

    ! Multipass lord 

    if (ele%lord_status == multipass_lord$ .and. .not. ele%field_master .and. ele%value(n_ref_pass$) == 0) then
      select case (ele%key)
        case (quadrupole$, sextupole$, octupole$, solenoid$, sol_quad$, sbend$, &
              hkicker$, vkicker$, kicker$, elseparator$, bend_sol_quad$)
        line = trim(line) // ', e_tot = ' // str(ele%value(e_tot$))
      end select
    endif


    call write_out (line, iu, .false.)  

    ! Encode taylor

    if (ele%key == taylor$) then
      do j = 1, 6
        unit_found = .false.
        unit = 0
        unit(j:j) = 1
        do k = 1, size(ele%taylor(j)%term)
          tm = ele%taylor(j)%term(k)
          write_term = .false.
          if (all(tm%exp == unit)) then
            unit_found = .true.
            if (tm%coef /= 1) write_term = .true.
          else
            write_term = .true.
          endif
          if (write_term) write (line, '(2a, i1, 3a, 6i2, a)') &
                trim(line), ', {', j, ': ', trim(str(tm%coef)), ',', tm%exp, '}'
        enddo
        if (.not. unit_found) write (line, '(2a, i1, a, 6i2, a)') &
                trim(line), ', {', j, ': 0,', tm%exp, '}'
      enddo
    endif

    if (associated(ele%a_pole)) then
      do j = 0, ubound(ele%a_pole, 1)
        if (ele%a_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0$)) // ' = ' // str(ele%a_pole(j))
        if (ele%b_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0$)) // ' = ' // str(ele%b_pole(j))
      enddo
    endif
    
    if (ele%key == wiggler$ .and. ele%sub_key == map_type$) then
      line = trim(line) // ', &'
      call write_out (line, iu, .true.)  
      do j = 1, size(ele%wig_term)
        wt = ele%wig_term(j)
        last = '}, &'
        if (j == size(ele%wig_term)) last = '}'
        write (iu, '(a, i3, 11a)') ' term(', j, ')={', trim(str(wt%coef)), ', ', &
          trim(str(wt%kx)), ', ', trim(str(wt%ky)), ', ', trim(str(wt%kz)), &
          ', ', trim(str(wt%phi_z)), trim(last)  
      enddo
    else
      call write_out (line, iu, .true.)  
    endif

  enddo ele_loop
enddo  ! branch loop

!----------------------------------------------------------
! Lattice Layout...

! Multipass stuff...

allocate (multipass(lat%n_ele_max))
multipass(:)%ix_region = 0
multipass(:)%region_start_pt = .false.
multipass(:)%region_stop_pt   = .false.

call multipass_all_info (lat, m_info)

if (size(m_info%top) /= 0) then

  ! Go through and mark all 1st pass regions
  ! In theory the original lattice file could have something like:
  !   lat: line = (..., m1, m2, ..., m1, -m2, ...)
  ! where m1 and m2 are multipass lines. The first pass region (m1, m2) looks 
  ! like this is one big region but the later (m1, -m2) signals that this 
  ! is not so.
  ! We thus go through all the first pass regions and compare them to the
  ! corresponding higher pass regions. If we find two elements that are contiguous
  ! in the first pass region but not contiguous in some higher pass region, 
  ! we need to break the first pass region into two.

  ix_r = 0
  in_multi_region = .false.

  do ie = 1, lat%n_ele_track+1
    ele => lat%ele(ie)
    ix_pass = m_info%bottom(ie)%ix_pass
    if (ix_pass /= 1 .or. ie == lat%n_ele_track+1) then  ! Not a first pass region
      if (in_multi_region) multipass(ie-1)%region_stop_pt = .true.
      in_multi_region = .false.
      cycle
    endif
    ! If start of a new region...
    if (.not. in_multi_region) then  
      ix_r = ix_r + 1
      multipass(ie)%ix_region = ix_r
      multipass(ie)%region_start_pt = .true.
      in_multi_region = .true.
      ix_top = m_info%bottom(ie)%ix_top(1)
      ix_super = m_info%bottom(ie)%ix_super(1)
      ss1 => m_info%top(ix_top)%slave(:,ix_super)
      cycle
    endif
    ix_top = m_info%bottom(ie)%ix_top(1)
    ix_super = m_info%bottom(ie)%ix_super(1)
    ss2 => m_info%top(ix_top)%slave(:, ix_super)
    do ix_pass = 2, size(ss1)
      s1 => ss1(ix_pass)%ele
      s2 => ss2(ix_pass)%ele
      if (abs(s1%ix_ele - s2%ix_ele) /= 1) then  ! If not contiguous then need a new region
        ix_r = ix_r + 1
        multipass(ie-1)%region_stop_pt = .true.
        multipass(ie)%region_start_pt = .true.
        exit
      endif
    enddo
    ss1 => ss2
    multipass(ie)%ix_region = ix_r
  enddo

  ! Each 1st pass region is now a valid multipass line.
  ! Write out this info.

  write (iu, *)
  write (iu, '(a)') '!-------------------------------------------------------'

  ix_r = 0
  in_multi_region = .false.

  do ie = 1, lat%n_ele_track

    ix_pass = m_info%bottom(ie)%ix_pass
    if (ix_pass /= 1) cycle 

    if (multipass(ie)%region_start_pt) then
      if (ix_r > 0) then
        line = line(:len_trim(line)-1) // ')'
        call write_out (line, iu, .true.)
      endif
      ix_r = ix_r + 1
      write (iu, *)
      write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ': line[multipass] = ('
    endif

    call write_line_element (line, iu, lat%ele(ie), lat)

  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_out (line, iu, .true.)

end if

! Main line.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

write (iu, *)
line = 'main_line: line = ('

in_multi_region = .false.
do ie = 1, lat%n_ele_track

  if (.not. m_info%bottom(ie)%multipass) then
    call write_line_element (line, iu, lat%ele(ie), lat)
    cycle
  endif

  ix_top = m_info%bottom(ie)%ix_top(1)
  ix_super = m_info%bottom(ie)%ix_super(1)
  ix1 = m_info%top(ix_top)%slave(1,ix_super)%ele%ix_ele
  ix_r = multipass(ix1)%ix_region

  ! If entering new multipass region
  if (.not. in_multi_region) then
    in_multi_region = .true.
    if (multipass(ix1)%region_start_pt) then
      write (line, '(2a, i2.2, a)') trim(line), ' multi_line_', ix_r, ','
      look_for = 'stop'
    else
      write (line, '(2a, i2.2, a)') trim(line), ' -multi_line_', ix_r, ','
      look_for = 'start'
    endif
  endif

  if (look_for == 'start' .and. multipass(ix1)%region_start_pt .or. &
      look_for == 'stop' .and. multipass(ix1)%region_stop_pt) then 
    in_multi_region = .false.
  endif

enddo

line = line(:len_trim(line)-1) // ')'
call write_out (line, iu, .true.)

write (iu, *)
write (iu, *) 'use, main_line'

! Branch lines

do ib = 1, ubound(lat%branch, 1)

  branch => lat%branch(ib)

  write (iu, *)
  line = trim(branch%name) // ': line = ('

  do ie = 1, branch%n_ele_track
    call write_line_element (line, iu, branch%ele(ie), lat) 
  enddo

enddo

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

expand_lat_out = .false.
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (ele%slave_status == super_slave$) cycle

  if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
    if (ele%value(dphi0$) == 0) cycle
    if (.not. expand_lat_out) call write_expand_lat_header
    write (iu, *)
    write (iu, '(3a)') trim(ele%name), '[dphi0] = ', trim(str(ele%value(dphi0$)))
  endif

  if (ele%key == patch$ .and. ele%ref_orbit /= 0) then
    if (.not. expand_lat_out) call write_expand_lat_header
    write (iu, *)
    if (ele%value(x_offset$) /= 0) write (iu, '(3a)') trim(ele%name), '[x_offset] = ', trim(str(ele%value(x_offset$)))
    if (ele%value(y_offset$) /= 0) write (iu, '(3a)') trim(ele%name), '[y_offset] = ', trim(str(ele%value(y_offset$)))
    if (ele%value(z_offset$) /= 0) write (iu, '(3a)') trim(ele%name), '[z_offset] = ', trim(str(ele%value(z_offset$)))
    if (ele%value(x_pitch$) /= 0)  write (iu, '(3a)') trim(ele%name), '[x_pitch] = ', trim(str(ele%value(x_pitch$)))
    if (ele%value(y_pitch$) /= 0)  write (iu, '(3a)') trim(ele%name), '[y_pitch] = ', trim(str(ele%value(y_pitch$)))
    if (ele%value(tilt$) /= 0)     write (iu, '(3a)') trim(ele%name), '[tilt]    = ', trim(str(ele%value(tilt$)))
  endif

enddo

! cleanup

close(iu)
deallocate (names)
deallocate (multipass)
deallocate (m_info%top, m_info%bottom)
if (present(err)) err = .false.

!--------------------------------------------------------------------------------
contains

subroutine write_expand_lat_header ()

write (iu, *)
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, *)
write (iu, '(a)') 'expand_lattice'
expand_lat_out = .true.

end subroutine

end subroutine

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine write_line_element (line, iu, ele, lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu
integer j, ix

!

if (ele%slave_status == super_slave$) then
  ! If a super_lord element starts at the beginning of this slave element,
  !  put in the null_ele marker 'x__' + lord_name for the superposition.
  do j = 1, ele%n_lord
    lord => pointer_to_lord (lat, ele, j)
    lord_name = lord%name
    m_lord => pointer_to_multipass_lord (lord, lat)
    if (associated(m_lord)) lord_name = lord%name
    slave => pointer_to_slave (lat, lord, 1) 
    if (slave%ix_ele == ele%ix_ele) then
      write (line, '(4a)') trim(line), ' x__', trim(lord_name), ',' 
    endif
  enddo
  write (line, '(2a, i3.3, a)') trim(line), ' slave_drift_', ele%ixx, ','

elseif (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord (lat, ele, 1)
  write (line, '(4a)') trim(line), ' ', trim(lord%name), ','

else
  write (line, '(4a)') trim(line), ' ', trim(ele%name), ','
endif

if (len_trim(line) > 80) call write_out(line, iu, .false.)

end subroutine

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function str(rel) result (str_out)

implicit none

real(rp) rel
integer pl
character(20) str_out
character(16) fmt

!

if (rel == 0) then
  str_out = '0'
  return
endif

pl = floor(log10(abs(rel)))

if (pl > 5) then
  fmt = '(2a, i1)'
  if (pl > 9) fmt = '(2a, i2)'
  write (str_out, fmt) trim(rchomp(rel/10.0**pl, 0)), 'E', pl

elseif (pl > -3) then
  str_out = rchomp(rel, pl)

else
  fmt = '(2a, i2)'
  if (pl < -9)  fmt = '(2a, i3)'
  write (str_out, fmt) trim(rchomp(rel*10.0**(-pl), 0)), 'E', pl

endif

end function

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function rchomp (rel, plc) result (out)

implicit none

real(rp) rel
character(16) out
character(8) :: fmt = '(f16.xx)'
integer it, plc, ix

!

write (fmt(6:7), '(i2.2)') 8-plc
write (out, fmt) rel
do it = 16, 1, -1
  if (out(it:it) == ' ') cycle
  if (out(it:it) == '0') then
    out(it:it) = ' '
    cycle
  endif
  if (out(it:it) == '.') out(it:it) = ' '
  call string_trim(out, out, ix)
  return
enddo

end function

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
! Input:
!   end_is_neigh -- Logical: If true then write out everything.
!                     Otherwise wait for a full line of 76 characters or so.

subroutine write_out (line, iu, end_is_neigh)

implicit none

character(*) line
integer i, iu
logical end_is_neigh
logical, save :: init = .true.

!

outer_loop: do 

  if (len_trim(line) < 76) then
    if (end_is_neigh) then
      call write_this (line)
      init = .true.
    endif
    return
  endif
      
  do i = 74, 1, -1
    if (line(i:i) == ',') then
      call write_this (line(:i) // ' &')
      line = line(i+1:)
      cycle outer_loop
    endif
  enddo

  do i = 75, len_trim(line)
    if (line(i:i) == ',') then
      call write_this (line(:i) // ' &')
      line = line(i+1:)
      cycle outer_loop
    endif
  enddo

  if (end_is_neigh) then
    call write_this (line)
    init = .true.
    return
  endif

enddo outer_loop


contains

subroutine write_this (line2)

character(*) line2

!

if (init) then
  init = .false.
  write (iu, '(a)') trim(line2)
else
  write (iu, '(2x, a)') trim(line2)
endif

end subroutine

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find1_indexx (name, names, n_max, ix_match, match_found)
!
! Subroutine to find a matching name in a list of names sorted in increasing
! alphabetical order.
!
! Input:
!   name     -- Character(40): Name to match to.
!   names(:) -- Character(40): Array of sorted names.
!   n_max    -- Integer Only names(1:n_max) are used.   
!
! Output:
!   ix_match  -- Integer: 
!                  If a match is found then:
!                      names(ix_match) = name
!                      names(ix_match-1) /= name
!                  If no match is found then:
!                      names(ix_match) > name  ! ix_match may be > size(names)
!                      names(ix_match-1) < name
!   match_found -- Logical: Set True if a match is found. False otherwise
!-

subroutine find1_indexx (name, names, n_max, ix_match, match_found)

implicit none

integer ix1, ix2, ix3, n_max, ix_match

character(40) name, names(:)
character(40) this_name

logical match_found

! simple case

match_found = .false.

if (n_max == 0) then
  ix_match = 1
  return
endif

!

ix1 = 1
ix3 = n_max

do

  ix2 = (ix1 + ix3) / 2 
  this_name = names(ix2)

  if (this_name == name) then
    do ! if there are duplicate names in the list choose the first one
      if (ix2 == 1) exit
      if (names(ix2-1) /= this_name) exit
      ix2 = ix2 - 1
    enddo
    ix_match = ix2
    match_found = .true.
    return
  elseif (this_name < name) then
    ix1 = ix2 + 1
  else
    ix3 = ix2 - 1
  endif
                     
  if (ix1 > ix3) then
    if (this_name < name) then
      ix_match = ix2 + 1
    else
      ix_match = ix2
    endif
    return
  endif

enddo

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+ 
! Subroutine bmad_to_mad_or_xsif (out_type, out_file_name, lat, &
!                           use_matrix_model, ix_start, ix_end, converted_lat, err)
!
! Subroutine to write a MAD-8, MAD-X, or XSIF lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
!
! Note: sol_quad elements are replaced by a drift-matrix-drift or solenoid-quad model.
! Note: wiggler elements are replaced by a drift-matrix-drift or drift-bend model.

! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   out_type      -- Character(*): Either 'XSIF', 'MAD-8' or 'MAD-X'.
!   out_file_name -- Character(*): Name of the mad output lattice file.
!   lat           -- lat_struct: Holds the lattice information.
!   use_matrix_model
!                 -- Logical, optional: Use a drift-matrix_drift model for wigglers
!                       and sol_quad elements? Default is False.
!   ix_start      -- Integer, optional: Starting index of lat%ele(i)
!                       used for output.
!   ix_end        -- Integer, optional: Ending index of lat%ele(i)
!                       used for output.
!
! Output:
!   converted_lat -- Lat_struct, optional: Equivalent Bmad lattice with wiggler and 
!                       sol_quad elements replaced by their respective models.
!   err           -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine bmad_to_mad_or_xsif (out_type, out_file_name, lat, &
                          use_matrix_model, ix_start, ix_end, converted_lat, err)

type (lat_struct), target :: lat, lat_out, lat_model
type (lat_struct), optional :: converted_lat
type (ele_struct), pointer :: ele, ele1, ele2, lord
type (ele_struct), save :: drift_ele, ab_ele, taylor_ele
type (taylor_term_struct) :: term

real(rp) field, hk, vk, tilt
real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx)

integer, optional :: ix_start, ix_end
integer i, j, j2, k, n, ix, i_unique, i_line, iout, iu, n_list, j_count, ix_ele
integer ie1, ie2, ios, t_count, ix_lord, ix1, ix2, n_lord

character(*) out_type, out_file_name
character(300) line
character(40), allocatable :: name_list(:)
character(1000) line_out
character(8) str
character(20) :: r_name = "bmad_to_mad_or_xsif"
character(2) continue_char, eol_char

logical init_needed
logical parsing
logical, optional :: use_matrix_model, err

! open file

if (present(err)) err = .true.
iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

! Init

if (out_type == 'MAD-X') then
  continue_char = ''
  eol_char = ';'
elseif (out_type == 'MAD-8' .or. out_type == 'XSIF') then
  continue_char = ' &'
  eol_char = ''
else
  call out_io (s_error$, r_name, 'BAD OUT_TYPE: ' // out_type)
  return
endif

call init_ele (drift_ele, drift$)
call init_ele (taylor_ele, taylor$)
call init_ele (ab_ele, ab_multipole$)
call multipole_init (ab_ele)

ie1 = 1    ! element index
if (present(ix_start)) ie1 = ix_start

ie2 = lat%n_ele_track
if (present(ix_end)) ie2 = ix_end

allocate (name_list(3*(ie2-ie1+1))) ! list of element names

! Transfer info to lat_out and make substitutions for sol_quad and wiggler elements.

lat_out = lat
j_count = 0    ! drift around solenoid or sol_quad index
t_count = 0    ! taylor element count.
i_unique = 1000

ix_ele = ie1 
do 

  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => lat_out%ele(ix_ele)

  ! If the name has more than 16 characters then replace the name by something shorter and unique.

  if (len_trim(ele%name) > 16) then
    call out_io (s_warn$, r_name, 'Shortening element name: ' // ele%name)
    i_unique = i_unique + 1
    write (ele%name, '(a, i0)') ele%name(1:11), i_unique
  endif

  ! Replace element name containing "/" or "#" with "_"

  do
    j = index (ele%name, '\')         ! '
    j = index (ele%name, '#')   
    if (j == 0) exit
    ele%name(j:j) = '_'
  enddo

  ! If there is a multipole component then put half at the beginning and half at the end

  if (associated(ele%a_pole) .and. ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
    call multipole_ele_to_ab (ele, lat%param%particle, ab_ele%a_pole, ab_ele%b_pole, .true.)
    ab_ele%a_pole = ab_ele%a_pole / 2
    ab_ele%b_pole = ab_ele%b_pole / 2
    deallocate (ele%a_pole, ele%b_pole)
    j_count = j_count + 1
    write (ab_ele%name,   '(a, i3.3)') 'MULTIPOLE_Z', j_count
    call insert_element (lat_out, ab_ele, ix_ele)
    call insert_element (lat_out, ab_ele, ix_ele+2)
    ie2 = ie2 + 2
    cycle
  endif

  ! Convert sol_quad_and wiggler elements.
  ! NOTE: FOR NOW SOL_QUAD  USES DRIFT-MATRIX-DRIFT MODEL!

  if (ele%key == wiggler$ .or. ele%key == sol_quad$) then
    if (logic_option(.false., use_matrix_model) .or. ele%key == sol_quad$) then

      drift_ele%value(l$) = -ele%value(l$) / 2
      call make_mat6 (drift_ele, lat_out%param)
      taylor_ele%mat6 = matmul(matmul(drift_ele%mat6, ele%mat6), drift_ele%mat6)
      call mat6_to_taylor (taylor_ele%vec0, taylor_ele%mat6, taylor_ele%taylor)

      ! Add drifts before and after wigglers and sol_quads so total length is invariant
      j_count = j_count + 1
      t_count = t_count + 1
      write (drift_ele%name, '(a, i3.3)') 'DRIFT_Z', j_count
      write (taylor_ele%name, '(a, i3.3)') 'SOL_QUAD', j_count
      drift_ele%value(l$) = ele%value(l$) / 2
      ele%key = -1 ! Mark for deletion
      call remove_eles_from_lat (lat_out)
      call insert_element (lat_out, drift_ele, ix_ele)
      call insert_element (lat_out, taylor_ele, ix_ele+1)
      call insert_element (lat_out, drift_ele, ix_ele+2)
      ie2 = ie2 + 2
      cycle

    ! Non matrix model...
    ! If the wiggler has been sliced due to superposition, throw 
    ! out the markers that caused the slicing.

    else
      if (ele%key == wiggler$) then
        if (ele%slave_status == super_slave$) then
          ! Create the wiggler model using the super_lord
          lord => pointer_to_lord(lat, ele, 1)
          call create_wiggler_model (lord, lat_model)
          ! Remove all the slave elements and markers in between.
          call out_io (s_warn$, r_name, &
              'Note: Not translating to MAD/XSIF the markers within wiggler: ' // lord%name)
          lord%key = -1 ! mark for deletion
          call find_element_ends (lat_out, lord, ele1, ele2)
          ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
          ! If the wiggler wraps around the origin we are in trouble.
          if (ix2 < ix1) then 
            call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
            call err_exit
          endif
          do i = ix1+1, ix2
            lat_out%ele(i)%key = -1  ! mark for deletion
          enddo
          ie2 = ie2 - (ix2 - ix1 - 1)
        else
          call create_wiggler_model (ele, lat_model)
        endif
      else
        call create_sol_quad_model (ele, lat_model)  ! NOT YET IMPLEMENTED!
      endif
      ele%key = -1 ! Mark for deletion
      call remove_eles_from_lat (lat_out)
      do j = 1, lat_model%n_ele_track
        call insert_element (lat_out, lat_model%ele(j), ix_ele+j-1)
      enddo
      ie2 = ie2 + lat_model%n_ele_track - 1
      cycle
    endif
  endif

enddo

! lat lattice name

write (iu, '(2a)') '! File generated by BMAD_TO_MAD/BMAD_TO_XSIF', trim(eol_char)
write (iu, '(3a)') '! Bmad Lattice File: ', trim(lat_out%input_file_name), trim(eol_char)
write (iu, '(3a)') '! Bmad Lattice: ', trim(lat_out%lattice), trim(eol_char)
write (iu, *)

! beam definition

ele => lat_out%ele(ie1-1)

write (iu, '(2a, 2(a, es12.5), a)')  &
      'beam_def: Beam, Particle = ', trim(particle_name(lat_out%param%particle)),  &
      ', Energy =', 1e-9*ele%value(E_TOT$), ', Npart =', lat_out%param%n_part, trim(eol_char)

write (iu, *)

! write element parameters

n_list = 0                          ! number of names stored in the list

do ix_ele = ie1, ie2

  ele => lat_out%ele(ix_ele)

  if (out_type == 'XSIF') then
    if (ele%key == elseparator$) ele%key = drift$  ! XSIF does not have elsep elements.
  endif

  ! do not make duplicate specs

  if (any (ele%name == name_list(1:n_list))) cycle

  ! Add to the list of elements

  n_list = n_list + 1
  name_list(n_list) = ele%name

  ! select key

  val => ele%value

  select case (ele%key)

  ! drift

  case (drift$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': drift, l =', val(l$)
  ! beambeam

  case (beambeam$)

    line_out = trim(ele%name) // ': beambeam'
    call value_to_line (line_out, val(sig_x$), 'sigx', 'es13.5', 'R')
    call value_to_line (line_out, val(sig_y$), 'sigy', 'es13.5', 'R')
    call value_to_line (line_out, val(x_offset$), 'xma', 'es13.5', 'R')
    call value_to_line (line_out, val(y_offset$), 'yma', 'es13.5', 'R')
    call value_to_line (line_out, val(charge$), 'charge', 'es13.5', 'R')

  ! elseparator

  case (elseparator$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': elseparator, l =', val(l$)
    hk = val(hkick$)
    vk = val(vkick$)

    if (hk /= 0 .or. vk /= 0) then

      ix = len_trim(line_out) + 1
      field = 1.0e3 * sqrt(hk**2 + vk**2) * val(E_TOT$) / val(l$)
      if (out_type == 'MAD-X') then
        write (line_out(ix:), '(a, es13.5)') ', ey =', field
      else
        write (line_out(ix:), '(a, es13.5)') ', e =', field
      endif

      if (lat_out%param%particle == positron$) then
        tilt = -atan2(hk, vk) + val(tilt$)
      else
        tilt = -atan2(hk, vk) + val(tilt$) + pi
      endif
      ix = len_trim(line_out) + 1
      write (line_out(ix:), '(a, es13.5)') ', tilt =', tilt

    endif

  ! kicker

  case (kicker$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': kicker, l =', val(l$)

    call value_to_line (line_out, val(hkick$), 'hkick', 'es13.5', 'R')
    call value_to_line (line_out, val(vkick$), 'vkick', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! hkicker

  case (hkicker$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': hkicker, l =', val(l$)

    call value_to_line (line_out, val(hkick$), 'kick', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! vkicker

  case (vkicker$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': vkicker, l =', val(l$)

    call value_to_line (line_out, val(vkick$), 'kick', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! marker

  case (marker$, branch$, photon_branch$)

    line_out = trim(ele%name) // ': marker'

  ! octupole

  case (octupole$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': octupole, l =', val(l$)

    call value_to_line (line_out, val(k3$), 'k3', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! quadrupole

  case (quadrupole$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': quadrupole, l =', val(l$)
    call value_to_line (line_out, val(k1$), 'k1', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! sbend

  case (sbend$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': sbend, l =', val(l$)

    call value_to_line (line_out, val(angle$), 'angle', 'es13.5', 'R')
    call value_to_line (line_out, val(e1$), 'e1', 'es13.5', 'R')
    call value_to_line (line_out, val(e2$), 'e2', 'es13.5', 'R')
    call value_to_line (line_out, val(k1$), 'k1', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! sextupole

  case (sextupole$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': sextupole, l =', val(l$)
    call value_to_line (line_out, val(k2$), 'k2', 'es13.5', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'es13.5', 'R')

  ! taylor

  case (taylor$)

    line_out = trim(ele%name) // ': matrix'

    do i = 1, 6
      do k = 1, size(ele%taylor(i)%term)
        term = ele%taylor(i)%term(k)

        select case (sum(term%exp))
        case (1)
          j = maxloc(term%exp, 1)
          if (out_type == 'MAD-8') then
            write (str, '(a, i0, a, i0, a)') 'rm(', i, ',', j, ')'
          elseif (out_type == 'MAD-X') then
            write (str, '(a, 2i0)') 'rm', i, j
          elseif (out_type == 'XSIF') then
            write (str, '(a, 2i0)') 'r', i, j
          endif
          call value_to_line (line_out, term%coef, str, 'es13.5', 'R')
          
        case (2)
          j = maxloc(term%exp, 1)
          term%exp(j) = term%exp(j) - 1
          j2 = maxloc(term%exp, 1)
          if (out_type == 'MAD-8') then
            write (str, '(a, 3(i0, a))') 'tm(', i, ',', j, ',', j2, ')'
          elseif (out_type == 'MAD-X') then
            write (str, '(a, 3i0)') 'tm', i, j, j2
          elseif (out_type == 'XSIF') then
            write (str, '(a, 3i0)') 't', i, j, j2
          endif
          call value_to_line (line_out, term%coef, str, 'es13.5', 'R')


        case default
          call out_io (s_error$, r_name, &
                  'TAYLOR TERM: \es12.2\ : \6i3\ ', &
                  'IN ELEMENT: ' // ele%name, &
                  'CANNOT BE CONVERTED TO MAD MATRIX TERM', &
                  r_array = (/ term%coef /), i_array = term%exp)
        end select
      enddo

      if (ele%mat6(i,i) == 0) then
        if (out_type == 'MAD-8') then
          write (str, '(a, i0, a, i0, a)') 'rm(', i, ',', i, ')'
        elseif (out_type == 'MAD-X') then
          write (str, '(a, 2i0)') 'rm', i, i
        elseif (out_type == 'XSIF') then
          write (str, '(a, 2i0)') 'r', i, i
        endif
        call value_to_line (line_out, 0.0_rp, str, 'es13.5', 'R')
      endif

    enddo

  ! rfcavity

  case (rfcavity$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': rfcavity, l =', val(l$)
    call value_to_line (line_out, val(voltage$)/1E6, 'volt', 'es13.5', 'R')
    call value_to_line (line_out, val(phi0$)+val(dphi0$), 'lag', 'es13.5', 'R')
    call value_to_line (line_out, val(harmon$), 'harmon', 'i8', 'I')

  ! lcavity

  case (lcavity$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': lcavity, l =', val(l$)
    call value_to_line (line_out, val(gradient$)*val(l$)/1e6, 'deltae', 'f11.4', 'R')
    call value_to_line (line_out, val(rf_frequency$)/1e6, 'freq', 'es13.5', 'R')
    call value_to_line (line_out, val(phi0$)+val(dphi0$), 'phi0', 'es13.5', 'R')

  ! solenoid

  case (solenoid$)

    write (line_out, '(a, es13.5)') trim(ele%name) // ': solenoid, l =', val(l$)
    call value_to_line (line_out, val(ks$), 'ks', 'es13.5', 'R')

  ! multipole

  case (multipole$, ab_multipole$)

    call multipole_ele_to_kt (ele, lat_out%param%particle, knl, tilts, .true.)
    write (line_out, '(a, es13.5)') trim(ele%name) // ': multipole'  
    do i = 0, 9
      write (str, '(a, i1, a)') 'K', i, 'L'
      call value_to_line (line_out, knl(i), str, 'es13.5', 'R')
      write (str, '(a, i1)') 'T', i
      call value_to_line (line_out, tilts(i), str, 'es13.5', 'R')
    enddo

  ! unknown

  case default

    call out_io (s_error$, r_name, 'UNKNOWN ELEMENT TYPE: ' // key_name(ele%key), &
                                  'CONVERTING TO MARKER')
    line_out = trim(ele%name) // ': marker'

  end select

  ! write element spec to file

  do
    if (len_trim(line_out) < 76) exit
    ix = index(line_out(61:), ' ')
    write (iu, '(2a)') line_out(:60+ix), trim(continue_char)
    line_out = '    ' // line_out(60+ix:)
  enddo

  write (iu, '(2a)') trim(line_out(:75)), trim(eol_char)

enddo

! Write the lattice line
! bmad has a limit of 4000 characters so we may need to break the lat into pieces.

i_unique = 1000
i_line = 0
init_needed = .true.
line = ' '

do n = ie1, ie2

  ele => lat_out%ele(n)

  if (init_needed) then
    write (iu, *)
    write (iu, *) '!---------------------------------', trim(eol_char)
    write (iu, *)
    i_line = i_line + 1
    write (line, '(a, i0, 2a)') 'line_', i_line, ': line = (', ele%name
    iout = 0
    init_needed = .false.
    cycle
  endif

  ix = len_trim(line) + len_trim(ele%name)

  if (ix > 75) then
    write (iu, '(3a)') trim(line), ',', trim(continue_char)
    iout = iout + 1
    line = '   ' // ele%name
  else
    line = trim(line) // ', ' // ele%name
  endif

  if (n == ie2 .or. iout > 48) then
    line = trim(line) // ')'
    write (iu, '(2a)') trim(line), trim(eol_char)
    line = ' '
    init_needed = .true.
  endif

enddo

write (iu, *)
write (iu, *) '!---------------------------------', trim(eol_char)
write (iu, *)
line = 'lat: line = (line_1'
do i = 2, i_line
  write (line, '(2a, i0)') trim(line), ', line_', i
enddo
line = trim(line) // ')'
write (iu, *) trim(line), trim(eol_char)
if (out_type == 'MAD-X') then
  write (iu, *) 'use, period=lat;'
else
  write (iu, *) 'use, lat'
endif

! Write twiss parameters for a linear lattice.

ele => lat_out%ele(ie1-1)
if (lat_out%param%lattice_type /= circular_lattice$) then
  write (iu, *)
  write (iu, *) '!---------------------------------', trim(eol_char)
  write (iu, *)
  write (iu, '(2(a, es12.5), 2a)') 'TWISS, betx =', ele%a%beta, ', bety =', ele%b%beta, ',', trim(continue_char)
  write (iu, '(5x, 2(a, es12.5), 2a)') 'alfx =', ele%a%alpha, ', alfy =', ele%b%alpha, ',', trim(continue_char)
  write (iu, '(5x, 2(a, es12.5), 2a)') 'dx =', ele%a%eta, ', dpx = ', ele%a%etap, ',', trim(continue_char)
  write (iu, '(5x, 2(a, es12.5), a)') 'dy =', ele%b%eta, ', dpy = ', ele%b%etap, trim(eol_char)
endif

!

call out_io (s_info$, r_name, 'Written ' // trim(out_type) // &
                                ' lattice file: ' // trim(out_file_name))

deallocate (name_list)
if (present(err)) err = .false.

if (present(converted_lat)) then
  converted_lat = lat_out
  converted_lat%n_ele_max = converted_lat%n_ele_track
  do i = 1, converted_lat%n_ele_track
    converted_lat%ele(i)%slave_status = free$
    converted_lat%ele(i)%n_lord = 0
    converted_lat%ele(i)%ic2_lord = converted_lat%ele(i)%ic1_lord - 1
  enddo
  converted_lat%n_control_max = 0
  converted_lat%n_ic_max = 0
endif

call deallocate_lat_pointers (lat_out)
call deallocate_lat_pointers (lat_model)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine value_to_line (line, value, str, fmt, typ)

use precision_def

implicit none

character(*) line, str, fmt
character fmt2*40, typ*1

real(rp) value

!

if (value == 0) return
fmt2 = '(4a,' // trim(fmt) // ')'
if (typ == 'R') then
  write (line, fmt2) trim(line), ', ', trim(str), ' =', value
elseif (typ == 'I') then
  write (line, fmt2) trim(line), ', ', trim(str), ' =', nint(value)
else
  print *, 'ERROR IN VALUE_TO_LINE. BAD "TYP": ', typ 
  call err_exit
endif
end subroutine

end module
