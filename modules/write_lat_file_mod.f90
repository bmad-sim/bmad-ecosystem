module write_lat_file_mod

use bmad_struct
use bmad_interface
use multipole_mod

private str, rchomp, write_out, element_out, bmad_to_mad_or_xsif

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
  integer ix_slave_series
  integer ix_pass
  integer ix_region
  logical region_start_pt
  logical region_stop_pt
end type

type (multipass_info_struct), allocatable :: multipass(:)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, super, slave, lord
type (ele_struct), save :: ele_init
type (wig_term_struct) wt
type (control_struct) ctl
type (taylor_term_struct) tm

real(rp) s0

character(*) bmad_file
character(4000) line
character(4) last
character(40) name, look_for
character(200) wake_name, file_name
character(40), allocatable :: names(:)
character(200), allocatable, save :: sr_wake_name(:), lr_wake_name(:)
character(40) :: r_name = 'write_bmad_lattice_file'

integer i, j, k, n, ix, iu, iuw, ios, ixs, n_sr, n_lr, ix1
integer unit(6), ix_names, ix_match, n_pass_max, n_1st_pass
integer ix_slave, ix_ss, ix_l, ixs1, ixs2, ix_r, ix_pass
integer ix_ss1, ix_ss2, ix_multi_lord
integer, allocatable :: ix_slave_series(:,:)

logical, optional :: err
logical unit_found, write_term, match_found, found, in_multi_region, expand_lat_out
logical is_multi_sup

! Init...
! Count the number of foreign wake files

if (present(err)) err = .true.
call init_ele (ele_init)

n_sr = 0
n_lr = 0
do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (.not. associated(ele%wake)) cycle
  if (ele%wake%sr_file(1:6) == 'xsif::') n_sr = n_sr + 1 
  if (ele%wake%lr_file(1:6) == 'xsif::') n_lr = n_lr + 1  
enddo
call re_allocate(sr_wake_name, n_sr, 200)
call re_allocate(lr_wake_name, n_lr, 200)

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

ele_loop: do i = 1, lat%n_ele_max

  ele => lat%ele(i)

  ix_multi_lord = multipass_lord_index(i, lat, ix_pass) 

  if (ele%key == null_ele$) cycle
  if (ele%control_type == multipass_slave$) cycle ! Ignore for now
  if (ele%control_type == super_lord$ .and. ix_multi_lord > 0) cycle
  if (ele%control_type == super_slave$ .and. ix_pass > 1) cycle

  if (i == lat%n_ele_track+1) then
    write (iu, *)
    write (iu, '(a)') '!-------------------------------------------------------'
    write (iu, '(a)') '! Overlays, groups, etc.'
    write (iu, *)
  endif

  ! For a super_slave just create a dummy drift. 

  if (ele%control_type == super_slave$) then
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

  if (ele%control_type == overlay_lord$ .or. ele%control_type == group_lord$) then
    if (ele%control_type == overlay_lord$) then
      write (line, '(2a)') trim(ele%name), ': overlay = {'
    else
      write (line, '(2a)') trim(ele%name), ': group = {'
    endif
    j_loop: do j = ele%ix1_slave, ele%ix2_slave
      ctl = lat%control(j)
      ix = ctl%ix_slave
      slave => lat%ele(ix)
      do k = ele%ix1_slave, j-1 ! do not use elements w/ duplicate names
        if (lat%ele(lat%control(k)%ix_slave)%name == slave%name) cycle j_loop
      enddo
      if (j == ele%ix1_slave) then
        write (line, '(3a)') trim(line), trim(slave%name)
      else
        write (line, '(3a)') trim(line), ', ', trim(slave%name)
      endif
      name = attribute_name(slave, ctl%ix_attrib)  
      if (name /= ele%attribute_name) &
              line = trim(line) // '[' // trim(name) // ']'
      if (ctl%coef /= 1) write (line, '(3a)') trim(line), '/', trim(str(ctl%coef))
    enddo j_loop
    line = trim(line) // '}'
    if (ele%attribute_name == ' ') then
      line = trim(line) // ', command'
    else
      line = trim(line) // ', ' // ele%attribute_name
    endif
    if (ele%control_type == overlay_lord$) then
      ix = ele%ix_value
      if (ele%value(ix) /= 0) write (line, '(3a)') &
                          trim(line), ' = ', str(ele%value(ix))
    endif
    call write_out (line, iu, .true.)
    cycle
  endif

  ! Girder

  if (ele%control_type == girder$) then
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

  ! other elements

  if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
  if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
  if (associated(ele%descrip)) line = trim(line) // &
                            ', descrip = "' // trim(ele%descrip) // '"'

  ! Create a null_ele element for a superposition and fill in the superposition
  ! information.

  is_multi_sup = .false.
  if (ele%control_type == multipass_lord$) then
    ix1 = lat%control(ele%ix1_slave)%ix_slave
    if (lat%ele(ix1)%control_type == super_lord$) is_multi_sup = .true.
  endif

  if (ele%control_type == super_lord$ .or. is_multi_sup) then
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
          write (iuw, *) '              Freq       R/Q      Q       m  Polarization_Angle'
          write (iuw, *) '              [Hz]  [Ohm/m^(2m)]             [Radians/2pi]'
          do n = lbound(ele%wake%lr, 1), ubound(ele%wake%lr, 1)
            if (ele%wake%lr(n)%polarized) then
              write (iuw, '(a, i0, a, 3es14.5, i6, f10.6)') 'lr(', n, ') =', &
                    ele%wake%lr(n)%freq_in, ele%wake%lr(n)%R_over_Q, &
                     ele%wake%lr(n)%Q, ele%wake%lr(n)%m, ele%wake%lr(n)%angle
            else
              write (iuw, '(a, i0, a, 3es14.5, i6, f10.6)') 'lr(', n, ') =', &
                    ele%wake%lr(n)%freq_in, ele%wake%lr(n)%R_over_Q, &
                    ele%wake%lr(n)%Q, ele%wake%lr(n)%m
            endif
          enddo
          close(iuw)
        endif
      endif

      line = trim(line) // ',  lr_wake_file = "' // trim(wake_name) // '"'

    endif

  endif

  ! Now for the rest of the element attributes.

  do j = 1, n_attrib_maxx

    ! Exclude dependent variables

    if (j == check_sum$ .and. ele%key /= patch$) cycle
    if (j == E_TOT$) cycle
    if (j == p0c$) cycle
    if (j == tilt_tot$) cycle
    if (j == x_pitch_tot$) cycle
    if (j == y_pitch_tot$) cycle
    if (j == x_offset_tot$) cycle
    if (j == y_offset_tot$) cycle
    if (j == s_offset_tot$) cycle


    select case (ele%key)
    case (beambeam$)
      if (j == bbi_const$) cycle
    case (elseparator$)
      if (j == e_field$) cycle
      if (j == voltage$) cycle
    case (lcavity$)
      if (j == e_loss$) cycle
      if (j == delta_e$) cycle
      if (j == p0c_start$) cycle
      if (j == E_TOT_START$) cycle
    case (wiggler$)
      if (j == k1$) cycle
      if (j == rho$) cycle
    case (sbend$)
      if (j == l_chord$) cycle
      if (j == angle$) cycle
      if (j == rho$) cycle
    end select

    if (ele%field_master) then
      select case (ele%key)
      case (quadrupole$)
        if (j == k1$) cycle
      case (sextupole$)
        if (j == k2$) cycle
      case (octupole$)
        if (j == k3$) cycle
      case (solenoid$)
        if (j == ks$) cycle
      case (sol_quad$) 
        if (j == ks$) cycle
        if (j == k1$) cycle
      case (sbend$)
        if (j == g$) cycle
        if (j == g_err$) cycle
      case (hkicker$)
        if (j == kick$) cycle
      case (vkicker$)
        if (j == kick$) cycle
      end select

      if (j == hkick$) cycle
      if (j == vkick$) cycle

    else
      select case (ele%key)
      case (quadrupole$)
        if (j == b1_gradient$) cycle
      case (sextupole$)
        if (j == b2_gradient$) cycle
      case (octupole$)
        if (j == b3_gradient$) cycle
      case (solenoid$)
        if (j == bs_field$) cycle
      case (sol_quad$) 
        if (j == bs_field$) cycle
        if (j == b1_gradient$) cycle
      case (sbend$)
        if (j == b_field$) cycle
        if (j == b_field_err$) cycle
      case (hkicker$)
        if (j == bl_kick$) cycle
      case (vkicker$)
        if (j == bl_kick$) cycle
      end select

      if (j == bl_hkick$) cycle
      if (j == bl_vkick$) cycle

    endif
      
    !

    if (ele%value(j) == 0) cycle
    line = trim(line) // ', ' // trim(attribute_name(ele, j)) // &
                                                  ' = ' // str(ele%value(j))

    if (attribute_name(ele, j) == null_name) then
      print *, 'ERROR IN WRITE_BMAD_LATTICE_FILE:'
      print *, '      ELEMENT: ', ele%name
      print *, '      HAS AN UNKNOWN ATTRIBUTE INDEX:', j
      stop
    endif

  enddo ! attribute loop

  if (ele%mat6_calc_method /= bmad_standard$) line = trim(line) // &
          ', mat6_calc_method = ' // calc_method_name(ele%mat6_calc_method)
  if (ele%tracking_method /= bmad_standard$) line = trim(line) // &
          ', tracking_method = ' // calc_method_name(ele%tracking_method)
  if (ele%symplectify) line = trim(line) // ', symplectify'
  if (.not. ele%is_on) line = trim(line) // ', is_on = False'
  call write_out (line, iu, .false.)  

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
              trim(line), ', {', j, ': ', trim(str(tm%coef)), ',', unit, '}'
        if (.not. unit_found) write (line, '(2a, i1, a, 6i2, a)') &
              trim(line), ', {', j, ': 0,', unit, '}'
      enddo
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

!----------------------------------------------------------
! Lattice Layout...

! Multipass stuff...
! First get an upper bound on the number of 1st pass slaves in the tracking lattice 
! and the maximum number of passes. 

n_pass_max = 0
n_1st_pass = 0

do i = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(i)
  if (lord%control_type /= multipass_lord$) cycle
  n_pass_max = max(n_pass_max, lord%n_slave) 
  ix_slave = lat%control(lord%ix1_slave)%ix_slave
  if (lat%ele(ix_slave)%control_type == super_lord$) then
    n_1st_pass = n_1st_pass + lat%ele(ix_slave)%n_slave
  else
    n_1st_pass = n_1st_pass + 1
  endif
enddo

allocate (multipass(lat%n_ele_max))
multipass(:)%ix_pass = 0
multipass(:)%ix_region = 0
multipass(:)%region_start_pt = .false.
multipass(:)%region_stop_pt   = .false.

if (n_1st_pass > 0) then

  ! Mark all the corresponding slaves

  ix_ss = 0
  allocate (ix_slave_series(n_1st_pass, n_pass_max))

  do i = lat%n_ele_track+1, lat%n_ele_max
    lord => lat%ele(i)
    if (lord%control_type /= multipass_lord$) cycle 
    ixs = lat%control(lord%ix1_slave)%ix_slave
    if (lat%ele(ixs)%control_type == super_lord$) then
      do j = lord%ix1_slave, lord%ix2_slave
        ix_pass = j + 1 - lord%ix1_slave
        ixs = lat%control(j)%ix_slave
        super => lat%ele(ixs)
        do k = super%ix1_slave, super%ix2_slave
          ix_ss2 = ix_ss + k + 1 - super%ix1_slave
          ixs = lat%control(k)%ix_slave
          ix_slave_series(ix_ss2, ix_pass) = ixs
          multipass(ixs)%ix_slave_series = ix_ss2
          multipass(ixs)%ix_pass = ix_pass
        enddo
      enddo
      ix_ss = ix_ss + super%n_slave
    else
      ix_ss = ix_ss + 1
      do j = lord%ix1_slave, lord%ix2_slave
        ix_pass = j + 1 - lord%ix1_slave
        ixs = lat%control(j)%ix_slave
        ix_slave_series(ix_ss, ix_pass) = ixs
        multipass(ixs)%ix_slave_series = ix_ss
        multipass(ixs)%ix_pass = ix_pass
      enddo
    endif
  enddo

  ! Now go through and mark all 1st pass regions
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

  do i = 1, lat%n_ele_track+1
    ele => lat%ele(i)
    ix_pass = multipass(i)%ix_pass
    if (ix_pass /= 1) then  ! Not a first pass region
      if (in_multi_region) multipass(i-1)%region_stop_pt = .true.
      in_multi_region = .false.
      cycle
    endif
    ! If start of a new region...
    if (.not. in_multi_region) then  
      ix_r = ix_r + 1
      multipass(i)%ix_region = ix_r
      multipass(i)%region_start_pt = .true.
      in_multi_region = .true.
      ix_ss1 = multipass(i)%ix_slave_series
      cycle
    endif
    ix_ss2 = multipass(i)%ix_slave_series
    do ix_pass = 2, size(ix_slave_series, 2)
      ixs1 = ix_slave_series(ix_ss1, ix_pass)
      ixs2 = ix_slave_series(ix_ss2, ix_pass)
      if (abs(ixs1 - ixs2) /= 1) then  ! If not contiguous then need a new region
        ix_r = ix_r + 1
        multipass(i-1)%region_stop_pt = .true.
        multipass(i)%region_start_pt = .true.
        exit
      endif
    enddo
    ix_ss1 = ix_ss2
    multipass(i)%ix_region = ix_r
  enddo

  ! Each 1st pass region is now a valid multipass line.
  ! Write out this info.

  write (iu, *)
  write (iu, '(a)') '!-------------------------------------------------------'

  ix_r = 0
  in_multi_region = .false.

  do i = 1, lat%n_ele_track

    ix_pass = multipass(i)%ix_pass
    if (ix_pass /= 1) cycle 

    if (multipass(i)%region_start_pt) then
      if (ix_r > 0) then
        line = line(:len_trim(line)-1) // ')'
        call write_out (line, iu, .true.)
      endif
      ix_r = ix_r + 1
      write (iu, *)
      write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ': line[multipass] = ('
    endif

    call write_line_element (line, iu, i, lat)

  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_out (line, iu, .true.)

  deallocate (ix_slave_series)
end if

! Main line.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

write (iu, *)
line = 'main_line: line = ('

in_multi_region = .false.
do i = 1, lat%n_ele_track

  if (multipass(i)%ix_pass == 0) then
    call write_line_element (line, iu, i, lat)
    cycle
  endif

  ix_ss = multipass(i)%ix_slave_series
  ix1 = ix_slave_series(ix_ss, 1)
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

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

if (n_pass_max > 0) then
  expand_lat_out = .false.
  do i = 1, lat%n_ele_max
    ele => lat%ele(i)
    if (ele%control_type == super_slave$) cycle
    if (ele%key /= lcavity$ .and. ele%key /= rfcavity$) cycle
    if (ele%value(dphi0$) == 0) cycle
    if (.not. expand_lat_out) then
      write (iu, *)
      write (iu, '(a)') '!-------------------------------------------------------'
      write (iu, *)
      write (iu, '(a)') 'expand_lattice'
      write (iu, *)
      expand_lat_out = .true.
    endif
    write (iu, '(3a)') trim(ele%name), '[dphi0] = ', trim(str(ele%value(dphi0$)))
  enddo
endif

! cleanup

close(iu)
deallocate (names)
deallocate (multipass)
if (present(err)) err = .false.

end subroutine

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine write_line_element (line, iu, ix_ele, lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord

character(*) line
character(40) lord_name

integer iu, ix_ele, ixm
integer j, ix, ic

!

ele => lat%ele(ix_ele)

if (ele%control_type == super_slave$) then
  ! If a super_lord element starts at the beginning of this slave element,
  !  put in the null_ele marker 'x__' + lord_name for the superposition.
  do j = ele%ic1_lord, ele%ic2_lord
    ix = lat%control(lat%ic(j))%ix_lord
    lord => lat%ele(ix)
    lord_name = lord%name
    ixm = multipass_lord_index(ix, lat)
    if (ixm > 0) lord_name = lat%ele(ixm)%name
    if (lat%control(lord%ix1_slave)%ix_slave == ix_ele) then
      write (line, '(4a)') trim(line), ' x__', trim(lord_name), ',' 
    endif
  enddo
  write (line, '(2a, i3.3, a)') trim(line), ' slave_drift_', ele%ixx, ','

elseif (ele%control_type == multipass_slave$) then
  ic = lat%ic(ele%ic1_lord)
  ix = lat%control(ic)%ix_lord
  lord => lat%ele(ix)
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

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine bmad_to_mad (mad_file, lat, ix_start, ix_end, err)
!
! Subroutine to write a mad lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
!
! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   mad_file    -- Character(*): Name of the mad output lattice file.
!   lat         -- lat_struct: Holds the lattice information.
!   ix_start    -- Integer, optional: Starting index of lat%ele(i)
!                       used for output.
!   ix_end      -- Integer, optional: Ending index of lat%ele(i)
!                       used for output.
!
! Output:
!   err    -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine bmad_to_mad (mad_file, lat, ix_start, ix_end, err)

implicit none

type (lat_struct)  lat
integer, optional :: ix_start, ix_end
character(*) mad_file 
logical, optional :: err

!

call bmad_to_mad_or_xsif ('MAD', mad_file, lat, ix_start, ix_end, err)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine bmad_to_xsif (xsif_file, lat, ix_start, ix_end, err)
!
! Subroutine to write a xsif lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
!
! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   xsif_file   -- Character(*): Name of the xsif output lattice file.
!   lat         -- lat_struct: Holds the lattice information.
!   ix_start    -- Integer, optional: Starting index of lat%ele(i)
!                       used for output.
!   ix_end      -- Integer, optional: Ending index of lat%ele(i)
!                       used for output.
!
! Output:
!   err    -- Logical: Set True if, say a file could not be opened.
!-

subroutine bmad_to_xsif (xsif_file, lat, ix_start, ix_end, err)

implicit none

type (lat_struct)  lat
integer, optional :: ix_start, ix_end
character(*) xsif_file 
logical, optional :: err

! 

call bmad_to_mad_or_xsif ('XSIF', xsif_file, lat, ix_start, ix_end, err)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine bmad_to_mad_or_xsif (out_type, out_file_name, lat, ix_start, ix_end, err)

type (lat_struct), target :: lat
type (ele_struct), save :: ele
type (ele_struct), save :: drift_ele, ab_ele

integer, optional :: ix_start, ix_end
integer i, j, n, ix, i_unique, i_line, iout, iu, n_list, j_count
integer ie1, ie2, ios

character(*) out_type, out_file_name
character(40) ele_name
character(300) line
character(20) :: r_name = "write_lattice_line"
character(40), allocatable :: name_list(:)

logical init_needed
logical, optional :: err

! open file

if (present(err)) err = .true.
iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

! lat lattice name

write (iu, '(2a)') '! File generated by BMAD_TO_', out_type
write (iu, '(2a)') '! Bmad Lattice File: ', trim(lat%input_file_name)
write (iu, '(2a)') '! Bmad Lattice: ', trim(lat%lattice)
write (iu, *)

! Init

call init_ele (drift_ele)
call init_ele (ab_ele)
call init_ele (ele)

drift_ele%key = drift$
ab_ele%key = ab_multipole$

j_count = 0    ! drift around solenoid or sol_quad index

ie1 = 1    ! element index
if (present(ix_start)) ie1 = ix_start

ie2 = lat%n_ele_track
if (present(ix_end)) ie2 = ix_end

allocate (name_list(3*(ie2-ie1+1))) ! list of element names
n_list = 0                          ! number of names stored in the list

! beam definition

ele = lat%ele(ie1-1)

write (iu, '(2a, 2(a, es12.5))')  &
      'beam_def: Beam, Particle = ', trim(particle_name(lat%param%particle)),  &
      ', Energy =', 1e-9*ele%value(E_TOT$), ', Npart =', lat%param%n_part

write (iu, *)

! Write element definitions...

i_unique = 1000

do i = ie1, ie2

  ele = lat%ele(i)

  ! If the name has more than 16 characters then replace the name by something shorter and unique.

  ele_name = ele%name
  if (len_trim(ele_name) > 16) then
    call out_io (s_warn$, r_name, 'Shortening element name: ' // ele%name)
    do 
      i_unique = i_unique + 1
      write (ele_name, '(a4, i0)') key_name(ele%key), i_unique
      if (all(ele_name /= name_list(1:n_list))) exit
    enddo      
  endif

  ! Replace element name containing "/" with "_"

  do
    j = index (ele_name, '\')         ! '
    if (j == 0) exit
    ele_name(j:j) = '_'
  enddo

  ! If there is a multipole component then put half at the beginning and half at the end

  if (associated(ele%a_pole) .and. ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
    ab_ele%a_pole = ele%a_pole / 2
    ab_ele%b_pole = ele%b_pole / 2
    write (ab_ele%name, '(a, i3.3)') 'MULTIPOLE_Z', j_count
    call element_out (out_type, ab_ele, ab_ele%name, lat, name_list, n_list, iu)
  endif

  ! Add drifts before and after wigglers and sol_quads so total length is invariant

  if (ele%key == wiggler$ .or. ele%key == sol_quad$) then
    j_count = j_count + 1
    write (drift_ele%name, '(a, i3.3)') 'DRIFT_Z', j_count
    drift_ele%value(l$) = ele%value(l$) / 2
    call element_out (out_type, drift_ele, drift_ele%name, lat, name_list, n_list, iu)

    drift_ele%value(l$) = -drift_ele%value(l$)
    call make_mat6 (drift_ele, lat%param)
    ele%mat6 = matmul(matmul(drift_ele%mat6, ele%mat6), drift_ele%mat6)
    call element_out (out_type, ele, ele_name, lat, name_list, n_list, iu)

    drift_ele%value(l$) = ele%value(l$) / 2
    call element_out (out_type, drift_ele, drift_ele%name, lat, name_list, n_list, iu)

  else
    call element_out (out_type, ele, ele_name, lat, name_list, n_list, iu)
  endif

  if (associated(ele%a_pole) .and. ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
    call element_out(out_type, ab_ele, ab_ele%name, lat, name_list, n_list, iu)
  endif

enddo

! Write the lattice line
! bmad has a limit of 4000 characters so we may need to break the lat into pieces.

i_unique = 1000
i_line = 0
init_needed = .true.
line = ' '

do n = 1, n_list

  ele_name = name_list(n)

  ix = len_trim(line)

  if (ix > 60 .or. n == n_list) then
    if (iout >= 50 .or. n == n_list) then
      i = len_trim(ele_name)
      line(ix+1:) = ', ' // ele_name(:i) // ')'
      write (iu, '(2a)') trim(line)
      line = ' '
      init_needed = .true.
    else
      write (iu, '(2a)') trim(line(:ix)), ', &'
      iout = iout + 1
      line = '   ' // ele_name
    endif

  elseif (init_needed) then
    write (iu, *)
    write (iu, *) '!---------------------------------'
    write (iu, *)
    i_line = i_line + 1
    write (line, '(a, i0, 2a)') 'line_', i_line, ': line = (', ele_name
    iout = 0
    init_needed = .false.

  else
    line(ix+1:) = ', ' // ele_name
  endif

enddo

write (iu, *)
write (iu, *) '!---------------------------------'
write (iu, *)
line = 'lat: line = (line_1'
do i = 2, i_line
  write (line, '(2a, i0)') trim(line), ', line_', i
enddo
line = trim(line) // ')'
write (iu, *) trim(line)
write (iu, *) 'use, lat'

! Write twiss parameters for a linear lattice.

ele = lat%ele(ie1-1)
if (lat%param%lattice_type /= circular_lattice$) then
  write (iu, *)
  write (iu, *) '!---------------------------------'
  write (iu, *)
  write (iu, '(3(a, es12.5))') &
      'TWISS, betx =', ele%a%beta, ', bety =', ele%b%beta, ', &'
  write (iu, '(5x, 3(a, es12.5))') &
      'alfx =', ele%a%alpha, ', alfy =', ele%b%alpha, ', &'
  write (iu, '(5x, 3(a, es12.5))') &
      'dx =', ele%a%eta, ', dpx = ', ele%a%etap, ', &'
  write (iu, '(5x, 3(a, es12.5))') &
      'dy =', ele%b%eta, ', dpy = ', ele%b%etap
endif

!

call out_io (s_info$, r_name, 'Written ' // trim(out_type) // &
                                ' lattice file: ' // trim(out_file_name))

deallocate (name_list)
if (present(err)) err = .false.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine element_out (out_type, ele, ele_name, lat, name_list, n_list, iu)

type (ele_struct), target :: ele
type (lat_struct) lat

integer i, j, n_list, ix, iu

character(*) out_type
character(1000) line
character(40) ele_name
character(40), allocatable :: name_list(:)
character(8) str

real(rp) field, hk, vk, tilt
real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx)

logical parsing

! Add to the list of elements

n_list = n_list + 1
name_list(n_list) = ele_name

! do not make duplicate specs

do i = 1, n_list - 1
  if (ele_name == name_list(i)) return
enddo

! select key

val => ele%value

select case (ele%key)

! drift

case (drift$)

  write (line, '(a, es13.5)') trim(ele_name) // ': drift, l =', val(l$)
! beambeam

case (beambeam$)

  line = trim(ele_name) // ': beambeam'
  call value_to_line (line, val(sig_x$), 'sigx', 'es13.5', 'R')
  call value_to_line (line, val(sig_y$), 'sigy', 'es13.5', 'R')
  call value_to_line (line, val(x_offset$), 'xma', 'es13.5', 'R')
  call value_to_line (line, val(y_offset$), 'yma', 'es13.5', 'R')
  call value_to_line (line, val(charge$), 'charge', 'es13.5', 'R')

! elseparator

case (elseparator$)

  write (line, '(a, es13.5)') trim(ele_name) // ': elseparator, l =', val(l$)
  hk = val(hkick$)
  vk = val(vkick$)

  if (hk /= 0 .or. vk /= 0) then

    ix = len_trim(line) + 1
    field = 1.0e3 * sqrt(hk**2 + vk**2) * val(E_TOT$) / val(l$)
    write (line(ix:), '(a, es13.5)') ', e =', field

    if (lat%param%particle == positron$) then
      tilt = -atan2(hk, vk) + val(tilt$)
    else
      tilt = -atan2(hk, vk) + val(tilt$) + pi
    endif
    ix = len_trim(line) + 1
    write (line(ix:), '(a, es13.5)') ', tilt =', tilt

  endif

! kicker

case (kicker$)

  write (line, '(a, es13.5)') trim(ele_name) // ': kicker, l =', val(l$)

  call value_to_line (line, val(hkick$), 'hkick', 'es13.5', 'R')
  call value_to_line (line, val(vkick$), 'vkick', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! hkicker

case (hkicker$)

  write (line, '(a, es13.5)') trim(ele_name) // ': hkicker, l =', val(l$)

  call value_to_line (line, val(hkick$), 'kick', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! vkicker

case (vkicker$)

  write (line, '(a, es13.5)') trim(ele_name) // ': vkicker, l =', val(l$)

  call value_to_line (line, val(vkick$), 'kick', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! marker

case (marker$, branch$, photon_branch$)

  line = trim(ele_name) // ': marker'

! octupole

case (octupole$)

  write (line, '(a, es13.5)') trim(ele_name) // ': octupole, l =', val(l$)

  call value_to_line (line, val(k3$), 'k3', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! quadrupole

case (quadrupole$)

  write (line, '(a, es13.5)') trim(ele_name) // ': quad, l =', val(l$)
  call value_to_line (line, val(k1$), 'k1', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! sbend

case (sbend$)

  write (line, '(a, es13.5)') trim(ele_name) // ': sbend, l =', val(l$)

  call value_to_line (line, val(angle$), 'angle', 'es13.5', 'R')
  call value_to_line (line, val(e1$), 'e1', 'es13.5', 'R')
  call value_to_line (line, val(e2$), 'e2', 'es13.5', 'R')
  call value_to_line (line, val(k1$), 'k1', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! sextupole

case (sextupole$)

  write (line, '(a, es13.5)') trim(ele_name) // ': sextupole, l =', val(l$)
  call value_to_line (line, val(k2$), 'k2', 'es13.5', 'R')
  call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! rfcavity

case (rfcavity$)

  write (line, '(a, es13.5)') trim(ele_name) // ': rfcavity, l =', val(l$)
  call value_to_line (line, val(voltage$)/1E6, 'volt', 'es13.5', 'R')
  call value_to_line (line, val(phi0$)+val(dphi0$), 'lag', 'es13.5', 'R')
  call value_to_line (line, val(harmon$), 'harmon', 'i8', 'I')

! lcavity

case (lcavity$)

  write (line, '(a, es13.5)') trim(ele_name) // ': lcavity, l =', val(l$)
  call value_to_line (line, val(gradient$)*val(l$)/1e6, 'deltae', 'f11.4', 'R')
  call value_to_line (line, val(rf_frequency$)/1e6, 'freq', 'es13.5', 'R')
  call value_to_line (line, val(phi0$)+val(dphi0$), 'phi0', 'es13.5', 'R')

! solenoid

case (solenoid$)

  write (line, '(a, es13.5)') trim(ele_name) // ': solenoid, l =', val(l$)
  call value_to_line (line, val(ks$), 'ks', 'es13.5', 'R')

! wiggler or solquad must be a matrix

case (wiggler$, sol_quad$, match$)

  line = trim(ele_name) // ': matrix'
  do i = 1, 6
    do j = 1, 6
      if (out_type == 'MAD') then
        write (str, '(a, i0, a, i0, a)') 'rm(', i, ',', j, ')'
      elseif (out_type == 'XSIF') then
        write (str, '(a, i0, i0)') 'r', i, j
      else
        print *, 'ERROR: BMAD OUT_TYPE: ', out_type
        call err_exit
      endif
      call value_to_line (line, ele%mat6(i,j), str, 'es13.5', 'R')
    enddo
  enddo

! multipole

case (multipole$, ab_multipole$)

  call multipole_ele_to_kt (ele, lat%param%particle, knl, tilts, .true.)
  write (line, '(a, es13.5)') trim(ele_name) // ': multipole'  
  do i = 0, 9
    write (str, '(a, i1, a)') 'K', i, 'L'
    call value_to_line (line, knl(i), str, 'es13.5', 'R')
    write (str, '(a, i1)') 'T', i
    call value_to_line (line, tilts(i), str, 'es13.5', 'R')
  enddo

! unknown

case default

  print *, 'ERROR: UNKNOWN ELEMENT: ', key_name(ele%key), ele%key
  print *, '       CONVERTING TO MARKER'

  line = trim(ele_name) // ': marker'

end select

! write element spec to file

do
  if (len_trim(line) < 76) exit
  ix = index(line(61:), ' ')
  write (iu, '(2a)') line(:60+ix), '&'
  line = '    ' // line(60+ix:)
enddo

write (iu, '(a)') trim(line(:75))

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
