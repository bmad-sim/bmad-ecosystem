!+
! Subroutine tao_show_cmd (what, stuff)
!
! Show information on variable, parameters, elements, etc.
!
! Input:
!   what  -- Character(*): What to show.
!   stuff -- Character(*): Particular stuff to show.
!-

module tao_show_mod

contains

!--------------------------------------------------------------------

recursive subroutine tao_show_cmd (what, stuff)

use tao_mod
use tao_top10_mod
use tao_command_mod, only: tao_cmd_split
use random_mod
use csr_mod, only: csr_com

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, save, target :: v_array(:)
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: c
type (tao_plot_region_struct), pointer :: region
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_ele_shape_struct), pointer :: shape

type (lr_wake_struct), pointer :: lr
type (ele_struct), pointer :: ele
type (coord_struct) orb
type (ele_struct) ele3

type show_lat_column_struct
  character(16) name
  character(16) format
  integer field_width
end type

type (show_lat_column_struct) column(40)

real(rp) f_phi, s_pos, l_lat
real(rp) :: delta_e = 0

character(*) :: what, stuff
character(24) :: var_name
character(24)  :: plane, imt, lmt, amt, rmt, irmt, iimt
character(80) :: word(2), fmt, fmt2, fmt3
character(8) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name
character(100), pointer :: ptr_lines(:)
character(100) file_name
character(200) stuff2
character(40) ele_name, name, sub_name
character(60) nam

character(16) :: show_names(16) = (/ &
   'data        ', 'var         ', 'global      ', 'alias       ', 'top10       ', &
   'optimizer   ', 'ele         ', 'lattice     ', 'constraints ', 'plot        ', &
   'write       ', 'hom         ', 'opt_vars    ', 'universe    ', 'taylor      ', &
   'beam        ' /)

character(200), allocatable, save :: lines(:)
character(200) line, line1, line2, line3
character(9) angle

integer :: data_number, ix_plane
integer nl, loc, ixl, iu, nc, n_size, ix_u, ios, ie
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju
integer num_locations
integer, allocatable, save :: ix_ele(:)

logical err, found, at_ends, first_time
logical show_all, name_found
logical, automatic :: picked(size(s%u))
logical, allocatable :: show_here(:)

namelist / custom_show_list / column

!

call re_allocate (ix_ele,1)
call re_allocate (lines, 200, 500)

err = .false.

lines = " "
nl = 0

rmt  = '(a, 9es16.8)'
irmt = '(a, i0, a, es16.8)'
imt  = '(a, 9i8)'
iimt = '(a, i0, a, i8)'
lmt  = '(a, 9l)'
amt  = '(9a)'

u => s%u(s%global%u_view)

if (s%global%phase_units == radians$) f_phi = 1
if (s%global%phase_units == degrees$) f_phi = 180 / pi
if (s%global%phase_units == cycles$)  f_phi = 1 / twopi

! find what to show

if (what == ' ') then
  call out_io (s_error$, r_name, 'SHOW WHAT?')
  return
endif

call match_word (what, show_names, ix)
if (ix == 0) then
  call out_io (s_error$, r_name, 'SHOW WHAT? WORD NOT RECOGNIZED: ' // what)
  return
endif

if (ix < 0) then
  call out_io (s_error$, r_name, 'SHOW WHAT? AMBIGUOUS: ' // what)
  return
endif

call tao_cmd_split (stuff, 2, word, .false., err)


select case (show_names(ix))

!----------------------------------------------------------------------
! beam

case ('beam')

  u => s%u(s%global%u_view)

  nl=nl+1; write(lines(nl), '(a, i3)') 'Universe: ', s%global%u_view
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), rmt) 'beam_init%a_norm_emitt      = ', u%beam_init%a_norm_emitt
  nl=nl+1; write(lines(nl), rmt) 'beam_init%b_norm_emitt      = ', u%beam_init%b_norm_emitt
  nl=nl+1; write(lines(nl), rmt) 'beam_init%dPz_dz            = ', u%beam_init%dPz_dz
  nl=nl+1; write(lines(nl), rmt) 'beam_init%center            = ', u%beam_init%center
  nl=nl+1; write(lines(nl), rmt) 'beam_init%ds_bunch          = ', u%beam_init%ds_bunch
  nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z             = ', u%beam_init%sig_z
  nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e             = ', u%beam_init%sig_e
  nl=nl+1; write(lines(nl), rmt) 'beam_init%bunch_charge      = ', u%beam_init%bunch_charge
  nl=nl+1; write(lines(nl), rmt) 'beam_init%center_jitter     = ', u%beam_init%center_jitter
  nl=nl+1; write(lines(nl), rmt) 'beam_init%emitt_jitter      = ', u%beam_init%emitt_jitter
  nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z_jitter      = ', u%beam_init%sig_z_jitter
  nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e_jitter      = ', u%beam_init%sig_e_jitter
  nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%polarization = ', u%beam_init%spin%polarization
  nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%theta        = ', u%beam_init%spin%theta
  nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%phi          = ', u%beam_init%spin%phi
  nl=nl+1; write(lines(nl), imt) 'beam_init%n_particle        = ', u%beam_init%n_particle
  nl=nl+1; write(lines(nl), imt) 'beam_init%n_bunch           = ', u%beam_init%n_bunch
  nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_center     = ', u%beam_init%renorm_center
  nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_sigma      = ', u%beam_init%renorm_sigma
  nl=nl+1; write(lines(nl), lmt) 'beam_init%preserve_dist     = ', u%beam_init%preserve_dist
  nl=nl+1; write(lines(nl), lmt) 'beam_init%init_spin         = ', u%beam_init%init_spin
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%sr_wakes_on               = ', bmad_com%sr_wakes_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%lr_wakes_on               = ', bmad_com%lr_wakes_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%trans_space_charge_on     = ', bmad_com%trans_space_charge_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%coherent_synch_rad_on     = ', bmad_com%coherent_synch_rad_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%spin_tracking_on          = ', bmad_com%spin_tracking_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%radiation_damping_on      = ', bmad_com%radiation_damping_on
  nl=nl+1; write(lines(nl), lmt) 'bmad_com%radiation_fluctuations_on = ', bmad_com%radiation_fluctuations_on
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), rmt) 'csr_com%ds_track_step      = ', csr_com%ds_track_step
  nl=nl+1; write(lines(nl), imt) 'csr_com%n_bin              = ', csr_com%n_bin
  nl=nl+1; write(lines(nl), imt) 'csr_com%particle_bin_span  = ', csr_com%particle_bin_span
  nl=nl+1; write(lines(nl), lmt) 'csr_com%lcsr_component_on  = ', csr_com%lcsr_component_on
  nl=nl+1; write(lines(nl), lmt) 'csr_com%lsc_component_on   = ', csr_com%lsc_component_on
  nl=nl+1; write(lines(nl), lmt) 'csr_com%tsc_component_on   = ', csr_com%tsc_component_on
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), amt) 'global%track_type           = ', s%global%track_type
  nl=nl+1; write(lines(nl), lmt) 'global%save_beam_everywhere = ', s%global%save_beam_everywhere

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! optimized_vars

case ('opt_vars')

  call tao_var_write (' ')

!----------------------------------------------------------------------
! hom

case ('hom')

  nl=nl+1; lines(nl) = &
        '       #        Freq         R/Q           Q   m  Polarization_Angle'
  do i = 1, size(u%model%lat%ele)
    ele => u%model%lat%ele(i)
    if (ele%key /= lcavity$) cycle
    if (ele%control_type == multipass_slave$) cycle
    nl=nl+1; write (lines(nl), '(a, i6)') ele%name, i
    do j = 1, size(ele%wake%lr)
      lr => ele%wake%lr(j)
      angle = '-'
      if (lr%polarized) write (angle, '(f9.4)') lr%angle
      nl=nl+1; write (lines(nl), '(i8, 3es12.4, i4, a)') j, &
                  lr%freq, lr%R_over_Q, lr%Q, lr%m, angle
    enddo
    nl=nl+1; lines(nl) = ' '
  enddo
  nl=nl+1; lines(nl) = '       #        Freq         R/Q           Q   m  Polarization_Angle'

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! write

case ('write')

  iu = lunget()
  file_name = s%global%write_file
  ix = index(file_name, '*')
  if (ix /= 0) then
    s%global%n_write_file = s%global%n_write_file + 1
    write (file_name, '(a, i3.3, a)') file_name(1:ix-1), &
                      s%global%n_write_file, trim(file_name(ix+1:))
  endif

  open (iu, file = file_name, position = 'APPEND', status = 'UNKNOWN')
  call output_direct (iu)  ! tell out_io to write to a file

  call out_io (s_blank$, r_name, ' ', 'Tao> show ' // stuff, ' ')
  call tao_show_cmd (word(1), word(2))  ! recursive

  call output_direct (0)  ! reset to not write to a file
  close (iu)
  call out_io (s_blank$, r_name, 'Written to file: ' // file_name)

  return

!----------------------------------------------------------------------
! alias

case ('alias')

  call re_allocate (lines, len(lines(1)), tao_com%n_alias+10)
  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, tao_com%n_alias
    nl=nl+1; lines(nl) = trim(tao_com%alias(i)%name) // ' = "' // &
                                    trim(tao_com%alias(i)%string) // '"'
  enddo
  
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, 'ALL')
  call tao_show_constraints (0, 'TOP10')

!----------------------------------------------------------------------
! data

case ('data')


! If just "show data" then show all names

  call tao_pick_universe (word(1), line1, picked, err)
  if (err) return

  if (line1 == ' ') then  ! just specified a universe

    do iu = 1, size(s%u)

      if (.not. picked(iu)) cycle

      u => s%u(s%global%u_view)

      nl=nl+1; write(lines(nl), *) ' '
      if (size(s%u) > 1) then
        nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', iu
      endif

      do i = 1, size(u%d2_data)
        d2_ptr => u%d2_data(i)
        if (d2_ptr%name == ' ') cycle
        nl=nl+1; lines(nl) = ' '
        do j = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
          d1_ptr => d2_ptr%d1(j)
          nl=nl+1; write (lines(nl), '(i5, 2x, 4a, i0, a, i0, a)') j, &
                      trim(d2_ptr%name), '.', trim(d1_ptr%name), &
                      '[', lbound(d1_ptr%d, 1), ':', ubound(d1_ptr%d, 1), ']'
        enddo
      enddo
    enddo

    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the data

  call tao_find_data (err, word(1), d2_ptr, d1_ptr, d_array)
  if (err) return

  n_size = 0
  if (allocated(d_array)) n_size = size(d_array)

! If d_ptr points to something then show the datum info.

  if (n_size == 1) then
    d_ptr => d_array(1)%d
    nl=nl+1; write(lines(nl), *) ' '
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d_ptr%d1%d2%ix_uni
    endif
    nl=nl+1; write(lines(nl), amt)  '%name              = ', d_ptr%name
    nl=nl+1; write(lines(nl), amt)  '%ele0_name         = ', d_ptr%ele0_name
    nl=nl+1; write(lines(nl), amt)  '%ele_name          = ', d_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%data_type         = ', d_ptr%data_type
    nl=nl+1; write(lines(nl), amt)  '%data_source       = ', d_ptr%data_source
    nl=nl+1; write(lines(nl), imt)  '%ix_ele0           = ', d_ptr%ix_ele0
    nl=nl+1; write(lines(nl), imt)  '%ix_ele            = ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_merit      = ', d_ptr%ix_ele_merit
    nl=nl+1; write(lines(nl), imt)  '%ix_dmodel         = ', d_ptr%ix_dModel
    nl=nl+1; write(lines(nl), imt)  '%ix_d1             = ', d_ptr%ix_d1
    nl=nl+1; write(lines(nl), imt)  '%ix_data           = ', d_ptr%ix_data
    nl=nl+1; write(lines(nl), imt)  '%ix_bunch          = ', d_ptr%ix_bunch
    nl=nl+1; write(lines(nl), rmt)  '%model             = ', d_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%design            = ', d_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%meas              = ', d_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%ref               = ', d_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%base              = ', d_ptr%base_value
    nl=nl+1; write(lines(nl), rmt)  '%old               = ', d_ptr%old_value   
    nl=nl+1; write(lines(nl), rmt)  '%fit               = ', d_ptr%fit_value
    nl=nl+1; write(lines(nl), rmt)  '%conversion_factor = ', d_ptr%conversion_factor
    nl=nl+1; write(lines(nl), rmt)  '%s                 = ', d_ptr%s
    nl=nl+1; write(lines(nl), amt)  '%merit_type        = ', d_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit             = ', d_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit       = ', d_ptr%delta_merit
    nl=nl+1; write(lines(nl), rmt)  '%weight            = ', d_ptr%weight
    nl=nl+1; write(lines(nl), lmt)  '%exists            = ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_model        = ', d_ptr%good_model
    nl=nl+1; write(lines(nl), lmt)  '%good_meas         = ', d_ptr%good_meas
    nl=nl+1; write(lines(nl), lmt)  '%good_ref          = ', d_ptr%good_ref
    nl=nl+1; write(lines(nl), lmt)  '%good_user         = ', d_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%good_opt          = ', d_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%good_plot         = ', d_ptr%good_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_plot        = ', d_ptr%useit_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_opt         = ', d_ptr%useit_opt

! Else show the d1_data info.

  elseif (associated(d1_ptr)) then

    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d1_ptr%d2%ix_uni
    endif
    
    nl=nl+1; write(lines(nl), '(2a)') 'Data name: ', trim(d2_ptr%name) // '.' // d1_ptr%name

    line1 = '                                                                      |   Useit'
    line2 = '     Name                         Meas         Model        Design    | Opt  Plot'
    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2

! if a range is specified, show the data range   

    call re_allocate (lines, len(lines(1)), nl+100+size(d1_ptr%d))

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      if (size(lines) > nl + 50) call re_allocate (lines, len(lines(1)), nl+100)
      nl=nl+1; write(lines(nl), '(i5, 2x, a20, 3es14.4, 2l6)') d_ptr%ix_d1, &
                     d_ptr%name, d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
    enddo

    nl=nl+1; lines(nl) = line2
    nl=nl+1; lines(nl) = line1

! else we must have a valid d2_ptr.

  elseif (associated(d2_ptr)) then

    call re_allocate (lines, len(lines(1)), nl+100+size(d2_ptr%d1))

    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d2_ptr%ix_uni
    endif
    nl=nl+1; write(lines(nl), '(2a)') 'D2_Data type:    ', d2_ptr%name
    nl=nl+1; write(lines(nl), '(5x, a)') '                   Bounds'
    nl=nl+1; write(lines(nl), '(5x, a)') 'D1_Data name    lower: Upper' 

    do i = 1, size(d2_ptr%d1)
      if (size(lines) > nl + 50) call re_allocate (lines, len(lines(1)), nl+100)
      nl=nl+1; write(lines(nl), '(5x, a, i5, a, i5)') d2_ptr%d1(i)%name, &
                  lbound(d2_ptr%d1(i)%d, 1), '.', ubound(d2_ptr%d1(i)%d, 1)
    enddo

    if (any(d2_ptr%descrip /= ' ')) then
      call re_allocate (lines, len(lines(1)), nl+100+size(d2_ptr%descrip))
      nl=nl+1; write (lines(nl), *)
      nl=nl+1; write (lines(nl), '(a)') 'Descrip:'
      do i = 1, size(d2_ptr%descrip)
        if (d2_ptr%descrip(i) /= ' ') then
          nl=nl+1; write (lines(nl), '(i4, 2a)') i, ': ', d2_ptr%descrip(i)
        endif
      enddo
    endif

! error

  else
    lines(1) = 'TRY BEING MORE SPECIFIC.'
    nl = 1
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! ele

case ('ele', 'taylor')

  call str_upcase (ele_name, word(1))

  if (index(ele_name, '*') /= 0 .or. index(ele_name, '%') /= 0) then
    write (lines(1), *) 'Matches to name:'
    nl = 1
    do loc = 1, u%model%lat%n_ele_max
      if (.not. match_wild(u%model%lat%ele(loc)%name, ele_name)) cycle
      if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
      nl=nl+1; write (lines(nl), '(i8, 2x, a)') loc, u%model%lat%ele(loc)%name
      name_found = .true.
    enddo
    if (.not. name_found) then
      nl=nl+1; write (lines(nl), *) '   *** No Matches to Name Found ***'
    endif

! else no wild cards

  else  

    call tao_locate_element (ele_name, s%global%u_view, ix_ele)
    loc = ix_ele(1)
    if (loc < 0) return

    write (lines(nl+1), *) 'Element #', loc
    nl = nl + 1

    ! Show the element info
    if (show_names(ix) == 'ele') then
      call type2_ele (u%model%lat%ele(loc), ptr_lines, n, .true., 6, .false., &
              s%global%phase_units, .true., u%model%lat, .true., .true., &
              s%global%show_ele_wig_terms)
    else
      call type2_ele (u%model%lat%ele(loc), ptr_lines, n, .true., 6, .true., &
              s%global%phase_units, .true., u%model%lat, .false., .false., &
              s%global%show_ele_wig_terms)
    endif
    if (size(lines) < nl+n+100) call re_allocate (lines, len(lines(1)), nl+n+100)
    lines(nl+1:nl+n) = ptr_lines(1:n)
    nl = nl + n
    deallocate (ptr_lines)
    if (show_names(ix) == 'ele') then
      nl=nl+1; lines(nl) = '[Conversion from Global to Screen: (Z, X) -> (-X, -Y)]'
    endif

    orb = u%model%orb(loc)
    fmt = '(2x, a, 3p2f11.4)'
    write (lines(nl+1), *) ' '
    write (lines(nl+2), *)   'Orbit: [mm, mrad]'
    write (lines(nl+3), fmt) "X  X':", orb%vec(1:2)
    write (lines(nl+4), fmt) "Y  Y':", orb%vec(3:4)
    write (lines(nl+5), fmt) "Z  Z':", orb%vec(5:6)
    nl = nl + 5

    ! Show data associated with this element
    call show_ele_data (u, loc, lines, nl)

    found = .false.
    do i = loc + 1, u%model%lat%n_ele_max
      if (u%model%lat%ele(i)%name /= ele_name) cycle
      if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
      if (found) then
        nl=nl+1; write (lines(nl), *)
        found = .true.
      endif 
      nl=nl+1;  write (lines(nl), *) &
                'Note: Found another element with same name at:', i
    enddo

  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! global

case ('global')

  nl=nl+1; write (lines(nl), imt) 'n_universes          = ', size(s%u)
  nl=nl+1; write (lines(nl), rmt) 'de_lm_step_ratio     = ', s%global%de_lm_step_ratio
  nl=nl+1; write (lines(nl), rmt) 'lm_opt_deriv_reinit  = ', s%global%lm_opt_deriv_reinit
  nl=nl+1; write (lines(nl), rmt) 'lmdif_eps            = ', s%global%lmdif_eps
  nl=nl+1; write (lines(nl), rmt) 'y_axis_plot_dmin     = ', s%global%y_axis_plot_dmin
  nl=nl+1; write (lines(nl), imt) 'u_view               = ', s%global%u_view
  nl=nl+1; write (lines(nl), imt) 'n_opti_loops         = ', s%global%n_opti_loops
  nl=nl+1; write (lines(nl), imt) 'n_opti_cycles        = ', s%global%n_opti_cycles
  nl=nl+1; write (lines(nl), imt) 'bunch_to_plot        = ', s%global%bunch_to_plot
  nl=nl+1; write (lines(nl), imt) 'random_seed          = ', s%global%random_seed
  if (s%global%random_seed == 0) then
    call ran_seed_get(ix)
    nl=nl+1; write (lines(nl), imt) 'random_seed (generated) = ', ix
  endif
  nl=nl+1; write (lines(nl), imt) 'n_curve_pts          = ', s%global%n_curve_pts
  nl=nl+1; write (lines(nl), amt) 'track_type           = ', s%global%track_type
  nl=nl+1; write (lines(nl), amt) 'phase_units          = ', &
                              frequency_units_name(s%global%phase_units)
  nl=nl+1; write (lines(nl), amt) 'optimizer            = ', s%global%optimizer
  nl=nl+1; write (lines(nl), amt) 'prompt_string        = ', s%global%prompt_string
  nl=nl+1; write (lines(nl), amt) 'var_out_file         = ', s%global%var_out_file
  nl=nl+1; write (lines(nl), amt) 'print_command        = ', s%global%print_command
  nl=nl+1; write (lines(nl), amt) 'init_file            = ', s%global%init_file
  nl=nl+1; write (lines(nl), amt) 'beam_file            = ', s%global%beam_file
  nl=nl+1; write (lines(nl), lmt) 'auto_scale           = ', s%global%auto_scale
  nl=nl+1; write (lines(nl), lmt) 'label_lattice_elements = ', s%global%label_lattice_elements
  nl=nl+1; write (lines(nl), lmt) 'label_keys           = ', s%global%label_keys
  nl=nl+1; write (lines(nl), lmt) 'derivative_recalc    = ', s%global%derivative_recalc
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on        = ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on        = ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), lmt) 'opt_with_ref         = ', s%global%opt_with_ref 
  nl=nl+1; write (lines(nl), lmt) 'opt_with_base        = ', s%global%opt_with_base
  nl=nl+1; write (lines(nl), lmt) 'plot_on              = ', s%global%plot_on
  nl=nl+1; write (lines(nl), lmt) 'matrix_recalc_on     = ', s%global%matrix_recalc_on
  nl=nl+1; write (lines(nl), lmt) 'var_limits_on        = ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), lmt) 'save_beam_everywhere = ', s%global%save_beam_everywhere
  nl=nl+1; write (lines(nl), lmt) 'use_saved_beam_in_tracking = ', &
                                              s%global%use_saved_beam_in_tracking
  nl=nl+1; write (lines(nl), lmt) 'show_ele_wig_terms   = ', s%global%show_ele_wig_terms

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! lattice

case ('lattice')
  
  if (word(1) .eq. ' ') then
    nl=nl+1; write (lines(nl), '(a, i3)') 'Universe: ', s%global%u_view
    nl=nl+1; write (lines(nl), '(a, i5, a, i5)') 'Number of elements to track through:', &
                                          1, '  through', u%model%lat%n_ele_track
    if (u%model%lat%n_ele_max .gt. u%model%lat%n_ele_track) then
      nl=nl+1; write (lines(nl), '(a, i5, a, i5)') 'Lord elements:   ', &
                        u%model%lat%n_ele_track+1, '  through', u%model%lat%n_ele_max
    else
      nl=nl+1; write (lines(nl), '(a)') "there are NO Lord elements"
    endif

    nl=nl+1; write (lines(nl), '(a, f0.3)') "Lattice length: ", u%model%lat%param%total_length

    if (u%is_on) then
      nl=nl+1; write (lines(nl), '(a)') 'This universe is turned ON'
    else
      nl=nl+1; write (lines(nl), '(a)') 'This universe is turned OFF'
    endif

    if (.not. u%model%lat%param%stable .or. .not. u%model%lat%param%stable) then
      nl=nl+1; write (lines(nl), '(a, l)') 'Model lattice stability: ', &
                                                            u%model%lat%param%stable
      nl=nl+1; write (lines(nl), '(a, l)') 'Design lattice stability:', &
                                                            u%design%lat%param%stable
      call out_io (s_blank$, r_name, lines(1:nl))
      return
    endif
 
    call radiation_integrals (u%model%lat, &
                                  u%model%orb, u%model%modes, u%ix_rad_int_cache)
    call radiation_integrals (u%design%lat, &
                                  u%design%orb, u%design%modes, u%ix_rad_int_cache)
    if (u%model%lat%param%lattice_type == circular_lattice$) then
      call chrom_calc (u%model%lat, delta_e, &
                          u%model%a%chrom, u%model%b%chrom, exit_on_error = .false.)
      call chrom_calc (u%design%lat, delta_e, &
                          u%design%a%chrom, u%design%b%chrom, exit_on_error = .false.)
    endif

    nl=nl+1; write (lines(nl), *)
    nl=nl+1; write (lines(nl), '(17x, a)') '       X          |            Y'
    nl=nl+1; write (lines(nl), '(17x, a)') 'Model     Design  |     Model     Design'
    fmt = '(1x, a10, 1p 2e11.3, 2x, 2e11.3, 2x, a)'
    fmt2 = '(1x, a10, 2f11.3, 2x, 2f11.3, 2x, a)'
    fmt3 = '(1x, a10, 2f11.4, 2x, 2f11.4, 2x, a)'
    f_phi = 1 / twopi
    l_lat = u%model%lat%param%total_length
    n = u%model%lat%n_ele_track
    if (u%model%lat%param%lattice_type == circular_lattice$) then
      nl=nl+1; write (lines(nl), fmt2) 'Q', f_phi*u%model%lat%ele(n)%a%phi, &
            f_phi*u%design%lat%ele(n)%a%phi, f_phi*u%model%lat%ele(n)%b%phi, &
            f_phi*u%design%lat%ele(n)%b%phi,  '! Tune'
      nl=nl+1; write (lines(nl), fmt2) 'Chrom', u%model%a%chrom, & 
            u%design%a%chrom, u%model%b%chrom, u%design%b%chrom, '! dQ/(dE/E)'
      nl=nl+1; write (lines(nl), fmt2) 'J_damp', u%model%modes%a%j_damp, &
          u%design%modes%a%j_damp, u%model%modes%b%j_damp, &
          u%design%modes%b%j_damp, '! Damping Partition #'
      nl=nl+1; write (lines(nl), fmt) 'Emittance', u%model%modes%a%emittance, &
          u%design%modes%a%emittance, u%model%modes%b%emittance, &
          u%design%modes%b%emittance, '! Meters'
    endif
    nl=nl+1; write (lines(nl), fmt) 'Alpha_damp', u%model%modes%a%alpha_damp, &
          u%design%modes%a%alpha_damp, u%model%modes%b%alpha_damp, &
          u%design%modes%b%alpha_damp, '! Damping per turn'
    nl=nl+1; write (lines(nl), fmt) 'I4', u%model%modes%a%synch_int(4), &
          u%design%modes%a%synch_int(4), u%model%modes%b%synch_int(4), &
          u%design%modes%b%synch_int(4), '! Radiation Integral'
    nl=nl+1; write (lines(nl), fmt) 'I5', u%model%modes%a%synch_int(5), &
          u%design%modes%a%synch_int(5), u%model%modes%b%synch_int(5), &
          u%design%modes%b%synch_int(5), '! Radiation Integral'

    nl=nl+1; write (lines(nl), *)
    nl=nl+1; write (lines(nl), '(19x, a)') 'Model     Design'
    fmt = '(1x, a12, 1p2e11.3, 3x, a)'
    nl=nl+1; write (lines(nl), fmt) 'Sig_E/E:', u%model%modes%sigE_E, &
              u%design%modes%sigE_E
    nl=nl+1; write (lines(nl), fmt) 'Energy Loss:', u%model%modes%e_loss, &
              u%design%modes%e_loss, '! Energy_Loss (eV / Turn)'
    nl=nl+1; write (lines(nl), fmt) 'J_damp:', u%model%modes%z%j_damp, &
          u%design%modes%z%j_damp, '! Longitudinal Damping Partition #'
    nl=nl+1; write (lines(nl), fmt) 'Alpha_damp:', u%model%modes%z%alpha_damp, &
          u%design%modes%z%alpha_damp, '! Longitudinal Damping per turn'
    nl=nl+1; write (lines(nl), fmt) 'Alpha_p:', u%model%modes%synch_int(1)/l_lat, &
                 u%design%modes%synch_int(1)/l_lat, '! Momentum Compaction'
    nl=nl+1; write (lines(nl), fmt) 'I1:', u%model%modes%synch_int(1), &
                 u%design%modes%synch_int(1), '! Radiation Integral'
    nl=nl+1; write (lines(nl), fmt) 'I2:', u%model%modes%synch_int(2), &
                 u%design%modes%synch_int(2), '! Radiation Integral'
    nl=nl+1; write (lines(nl), fmt) 'I3:', u%model%modes%synch_int(3), &
                 u%design%modes%synch_int(3), '! Radiation Integral'

    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! Here to show info on particular elements

  column(:)%name = ""
  column(1)  = show_lat_column_struct("index",   "i6",     7)
  column(2)  = show_lat_column_struct("name",    "a24",   25)
  column(3)  = show_lat_column_struct("key",     "a16",   16)
  column(4)  = show_lat_column_struct("s",       "f10.3", 10)
  column(5)  = show_lat_column_struct("beta_a",  "f7.2",   7)
  column(6)  = show_lat_column_struct("phi_a",   "f8.3",   8)
  column(7)  = show_lat_column_struct("eta_a",   "f5.1",   5)
  column(8)  = show_lat_column_struct("orbit_x", "3p, f8.3",   8)
  column(9)  = show_lat_column_struct("beta_b",  "f7.2",   7)
  column(10) = show_lat_column_struct("phi_b",   "f8.3",   8)
  column(11) = show_lat_column_struct("eta_b",   "f5.1",   5)
  column(12) = show_lat_column_struct("orbit_y", "3p, f8.3",   8)

  call string_trim(stuff, stuff2, ix)
  at_ends = .true.
  allocate (show_here(0:u%model%lat%n_ele_max))

  do
    if (stuff2(1:ix) == 'middle') then
      call string_trim(stuff2(ix+1:), stuff2, ix)
      at_ends = .false.
    elseif (stuff2(1:ix) == 'custom') then
      call string_trim(stuff2(ix+1:), stuff2, ix)
      file_name = stuff2(1:ix)
      call string_trim(stuff2(ix+1:), stuff2, ix)
      iu = lunget()
      open (iu, file = file_name, status = 'old', iostat = ios)
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name)
        return
      endif
      column(:)%name = ""
      read (iu, nml = custom_show_list, iostat = ios)
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT READ "CUSTOM_SHOW_LIST" NAMELIST IN FILE: ' // file_name)
        return
      endif
    else
      exit
    endif
  enddo
  
  if (ix == 0 .or. stuff2(1:ix) == 'all') then
    show_here = .true.
  else
    call location_decode (stuff2, show_here, 0, num_locations)
    if (num_locations .eq. -1) then
      call out_io (s_error$, r_name, "Syntax error in range list!")
      deallocate(show_here)
      return
    endif
  endif

  if (at_ends) then
    at_ends = .true.
    write (line1, '(6x, a)') 'Model values at End of Element:'
  else
    at_ends = .false.
    write (line1, '(6x, a)') 'Model values at Center of Element:'
  endif

  ix = 1
  line2 = ""
  line3 = ""
  do i = 1, size(column)
    if (column(i)%name == "") cycle
    ix2 = ix + column(i)%field_width
    select case (column(i)%name)
    case ("index")
      line2(ix:) = "Ix" 
    case ("name")
      line2(ix:) = "Name" 
    case ("key")
      line2(ix:) = "key" 
    case ("s")
      line2(ix2-3:) = "S" 
    case ("beta_a")
      line2(ix2-4:) = "Beta" 
      line3(ix2-4:) = "  A "
    case ("beta_b")
      line2(ix2-4:) = "Beta" 
      line3(ix2-4:) = "  B "
    case ("alpha_a")
      line2(ix2-5:) = "Alpha" 
      line3(ix2-5:) = "   A "
    case ("alpha_b")
      line2(ix2-5:) = "Alpha" 
      line3(ix2-5:) = "   B "
    case ("phi_a")
      line2(ix2-3:) = "Phi" 
      line3(ix2-3:) = " A "
    case ("phi_b")
      line2(ix2-3:) = "Phi" 
      line3(ix2-3:) = " B "
    case ("eta_a")
      line2(ix2-3:) = "Eta" 
      line3(ix2-3:) = " A "
    case ("eta_b")
      line2(ix2-3:) = "Eta" 
      line3(ix2-3:) = " B "
    case ("etap_a")
      line2(ix2-4:) = "Etap" 
      line3(ix2-4:) = "  A "
    case ("etap_b")
      line2(ix2-4:) = "Etap" 
      line3(ix2-4:) = "  B "
    case ("orbit_x")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "   X "
    case ("orbit_px")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "  Px "
    case ("orbit_y")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "   Y "
    case ("orbit_py")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "  Py "
    case ("orbit_z")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "   Z "
    case ("orbit_pz")
      line2(ix2-5:) = "Orbit" 
      line3(ix2-5:) = "  Pz "
    case ("x")
      
    case default
      call out_io (s_error$, r_name, 'BAD NAME FOUND IN COLUMN SPEC: ' // column(i)%name)
      return
    end select
    column(i)%format = "(" // trim(column(i)%format) // ")"
    ix = ix2
  enddo

  lines(nl+1) = line1
  lines(nl+2) = line2
  lines(nl+3) = line3
  nl=nl+3

  do ie = 0, u%model%lat%n_ele_track
    if (.not. show_here(ie)) cycle
    if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
    ele => u%model%lat%ele(ie)
    if (ie == 0 .or. at_ends) then
      ele3 = ele
      orb = u%model%orb(ie)
      s_pos = ele3%s
    else
      call twiss_and_track_partial (u%model%lat%ele(ie-1), ele, &
                u%model%lat%param, ele%value(l$)/2, ele3, u%model%orb(ie-1), orb)
      s_pos = ele%s-ele%value(l$)/2
    endif

    line = ""
    ix = 1
    do i = 1, size(column)
      if (column(i)%name == "") cycle
      select case (column(i)%name)
      case ("index")
        write (line(ix:), column(i)%format, iostat = ios) ie
      case ("name")
        write (line(ix:), column(i)%format, iostat = ios) ele%name
      case ("key")
        write (line(ix:), column(i)%format, iostat = ios) key_name(ele%key)
      case ("s")
        write (line(ix:), column(i)%format, iostat = ios) s_pos
      case ("beta_a")
        write (line(ix:), column(i)%format, iostat = ios) ele3%a%beta
      case ("beta_b")
        write (line(ix:), column(i)%format, iostat = ios) ele3%b%beta
      case ("alpha_a")
        write (line(ix:), column(i)%format, iostat = ios) ele3%a%alpha
      case ("alpha_b")
        write (line(ix:), column(i)%format, iostat = ios) ele3%b%alpha
      case ("phi_a")
        write (line(ix:), column(i)%format, iostat = ios) ele3%a%phi
      case ("phi_b") 
        write (line(ix:), column(i)%format, iostat = ios) ele3%b%phi
      case ("eta_a")
        write (line(ix:), column(i)%format, iostat = ios) ele3%a%eta
      case ("eta_b")
        write (line(ix:), column(i)%format, iostat = ios) ele3%b%eta
      case ("etap_a")
        write (line(ix:), column(i)%format, iostat = ios) ele3%a%etap
      case ("etap_b")
        write (line(ix:), column(i)%format, iostat = ios) ele3%b%etap
      case ("orbit_x")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(1)
      case ("orbit_px")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(2)
      case ("orbit_y")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(3)
      case ("orbit_py")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(4)
      case ("orbit_z")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(5)
      case ("orbit_pz")
        write (line(ix:), column(i)%format, iostat = ios) orb%vec(6)
      case ("x")

      end select
      ix  = ix + column(i)%field_width
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'BAD FORMAT: ' // column(i)%format, &
                                       'FOR DISPLAYING: ' // column(i)%name)
        return
      endif
    enddo

    nl=nl+1; lines(nl) = line

  enddo

  lines(nl+1) = line2
  lines(nl+2) = line3
  lines(nl+3) = line1
  nl=nl+3

  first_time = .true.  
  do ie = u%model%lat%n_ele_track+1, u%model%lat%n_ele_max
    if (.not. show_here(ie)) cycle
    if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
    ele => u%model%lat%ele(ie)
    if (first_time) then
      nl=nl+1; lines(nl) = ' '
      nl=nl+1; lines(nl) = 'Lord Elements:'
      first_time = .false.
    endif
    nl=nl+1
    write (lines(nl), '(i6, 1x, a24, 1x, a16, f10.3, 2(f7.2, f8.3, f5.1, f8.3))') &
          ie, ele%name, key_name(ele%key)
  enddo



  call out_io (s_blank$, r_name, lines(1:nl))

  deallocate(show_here)

!----------------------------------------------------------------------
! optimizer

case ('optimizer')

  do i = 1, size(s%u)
    u => s%u(i)
    call out_io (s_blank$, r_name, ' ', 'Data Used:')
    write (lines(1), '(a, i4)') 'Universe: ', i
    if (size(s%u) > 1) call out_io (s_blank$, r_name, lines(1))
    do j = 1, size(u%d2_data)
      if (u%d2_data(j)%name == ' ') cycle
      call tao_data_show_use (u%d2_data(j))
    enddo
  enddo

  call out_io (s_blank$, r_name, ' ', 'Variables Used:')
  do j = 1, size(s%v1_var)
    if (s%v1_var(j)%name == ' ') cycle
    call tao_var_show_use (s%v1_var(j))
  enddo

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; write (lines(nl), amt) 'optimizer:        ', s%global%optimizer
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! plots

case ('plot')

! word(1) is blank => print overall info

  if (word(1) == ' ') then

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Templates:        Plot.Graph'
    nl=nl+1; lines(nl) = '             --------- ----------'
    do i = 1, size(s%template_plot)
      p => s%template_plot(i)
      if (p%name == ' ') cycle
      ix = 21 - len_trim(p%name)
      name = ' '
      name(ix:) = trim(p%name)
      if (allocated(p%graph)) then
        do j = 1, size(p%graph)
          nl=nl+1; write (lines(nl), '(2x, 3a)') name(1:20), '.', p%graph(j)%name
          name = ' '
        enddo
      else
        nl=nl+1; write (lines(nl), '(2x, 2a)') name(1:20), '.'
      endif
      nl=nl+1; lines(nl) = ' '
    enddo

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = '[Visible]     Plot Region         <-->  Template' 
    nl=nl+1; lines(nl) = '---------     -----------               ------------'
    do i = 1, size(s%plot_page%region)
      region => s%plot_page%region(i)
      nl=nl+1; write (lines(nl), '(3x l1, 10x, a20, 2a)') region%visible, &
                                    region%name, '<-->  ', region%plot%name
    enddo

    ! shapes

    if (s%plot_page%ele_shape(1)%key > 0) then
      nl=nl+1; lines(nl) = ' '
      nl=nl+1; lines(nl) = 'Element Shapes:'
      nl=nl+1; lines(nl) = &
            'Key             Ele_Name        Shape         Color        dy_pix   Draw_Name?'
      nl=nl+1; lines(nl) = &
            '-----------     ------------    --------      -----        -------  ---------'

      do i = 1, size(s%plot_page%ele_shape)
        shape => s%plot_page%ele_shape(i)
        if (shape%key < 1) cycle
        nl=nl+1; write (lines(nl), '(4a, f10.4, l3)') shape%key_name(1:16), &
                  shape%ele_name(1:16), shape%shape(1:14), shape%color(1:10), &
                  shape%dy_pix, shape%draw_name
      enddo
    endif

    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! Find particular plot

  call tao_find_plots (err, word(1), 'BOTH', plot, graph, curve, print_flag = .false.)
  if (err) return

! print info on particular plot, graph, or curve

  if (allocated(curve)) then
    c => curve(1)%c
    g => c%g
    p => g%p
    if (associated(p%r)) then
      nl=nl+1; lines(nl) = 'Region.Graph.Curve: ' // trim(p%r%name) // '.' // &
                                                  trim(g%name) // '.' // c%name
    endif
    nl=nl+1; lines(nl) = 'Plot.Graph.Curve:   ' // trim(p%name) // '.' // &
                                                  trim(g%name) // '.' // c%name
    nl=nl+1; write (lines(nl), amt) 'name                    = ', c%name
    nl=nl+1; write (lines(nl), amt) 'data_source             = ', c%data_source
    nl=nl+1; write (lines(nl), amt) 'data_type               = ', c%data_type
    nl=nl+1; write (lines(nl), amt) 'ele_ref_name            = ', c%ele_ref_name
    nl=nl+1; write (lines(nl), imt) 'ix_ele_ref              = ', c%ix_ele_ref
    nl=nl+1; write (lines(nl), imt) 'ix_ele_ref_track        = ', c%ix_ele_ref_track
    nl=nl+1; write (lines(nl), imt) 'ix_bunch                = ', c%ix_bunch
    nl=nl+1; write (lines(nl), imt) 'ix_universe             = ', c%ix_universe
    nl=nl+1; write (lines(nl), imt) 'symbol_every            = ', c%symbol_every
    nl=nl+1; write (lines(nl), rmt) 'x_axis_scale_factor     = ', c%x_axis_scale_factor
    nl=nl+1; write (lines(nl), rmt) 'y_axis_scale_factor     = ', c%y_axis_scale_factor
    nl=nl+1; write (lines(nl), lmt) 'use_y2                  = ', c%use_y2
    nl=nl+1; write (lines(nl), lmt) 'draw_line               = ', c%draw_line
    nl=nl+1; write (lines(nl), lmt) 'draw_symbols            = ', c%draw_symbols
    nl=nl+1; write (lines(nl), lmt) 'convert                 = ', c%convert
    nl=nl+1; write (lines(nl), lmt) 'draw_interpolated_curve = ', c%draw_interpolated_curve
    

  elseif (allocated(graph)) then
    g => graph(1)%g
    p => g%p
    if (associated(p%r)) then
      nl=nl+1; lines(nl) = 'Region.Graph: ' // trim(p%r%name) // '.' // trim(g%name)
    endif
    nl=nl+1; lines(nl) = 'Plot.Graph:   ' // trim(p%name) // '.' // trim(g%name)
    nl=nl+1; write (lines(nl), amt) 'name                 = ', g%name
    nl=nl+1; write (lines(nl), amt) 'type                 = ', g%type
    nl=nl+1; write (lines(nl), amt) 'title                = ', g%title
    nl=nl+1; write (lines(nl), amt) 'title_suffix         = ', g%title_suffix
    nl=nl+1; write (lines(nl), '(a, 4f10.2, 2x, a)') &
                                    'margin               = ', g%margin
    nl=nl+1; write (lines(nl), imt) 'box                  = ', g%box
    nl=nl+1; write (lines(nl), imt) 'ix_universe          = ', g%ix_universe
    nl=nl+1; write (lines(nl), imt) 'box                  = ', g%box
    nl=nl+1; write (lines(nl), lmt) 'valid                = ', g%valid
    nl=nl+1; write (lines(nl), lmt) 'y2_mirrors_y         = ', g%y2_mirrors_y
    nl=nl+1; write (lines(nl), rmt) 'y%max                = ', g%y%max
    nl=nl+1; write (lines(nl), rmt) 'y%min                = ', g%y%min
    nl=nl+1; write (lines(nl), imt) 'y%major_div          = ', g%y%major_div
    nl=nl+1; write (lines(nl), imt) 'y%places             = ', g%y%places
    nl=nl+1; write (lines(nl), lmt) 'y%draw_label         = ', g%y%draw_label
    nl=nl+1; write (lines(nl), lmt) 'y%draw_numbers       = ', g%y%draw_numbers
    nl=nl+1; write (lines(nl), rmt) 'y2%max               = ', g%y2%max
    nl=nl+1; write (lines(nl), rmt) 'y2%min               = ', g%y2%min
    nl=nl+1; write (lines(nl), imt) 'y2%major_div         = ', g%y2%major_div
    nl=nl+1; write (lines(nl), imt) 'y2%places            = ', g%y2%places
    nl=nl+1; write (lines(nl), lmt) 'y2%draw_label        = ', g%y2%draw_label
    nl=nl+1; write (lines(nl), lmt) 'y2%draw_numbers      = ', g%y2%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'limited              = ', g%limited
    nl=nl+1; write (lines(nl), lmt) 'clip                 = ', g%clip
    nl=nl+1; write (lines(nl), lmt) 'draw_axes            = ', g%draw_axes
    nl=nl+1; lines(nl) = 'Curves:'
    do i = 1, size(g%curve)
      nl=nl+1; write (lines(nl), amt) '   ', g%curve(i)%name
    enddo

  elseif (allocated(plot)) then
    p => plot(1)%p
    if (associated(p%r)) then
      nl=nl+1; lines(nl) = 'Region:  ' // trim(p%r%name)
    endif
    nl=nl+1; lines(nl) = 'Plot:  ' // p%name
    nl=nl+1; write (lines(nl), amt) 'x_axis_type          = ', p%x_axis_type
    nl=nl+1; write (lines(nl), rmt) 'x_divisions          = ', p%x_divisions
    nl=nl+1; write (lines(nl), rmt) 'x%max                = ', p%x%max
    nl=nl+1; write (lines(nl), rmt) 'x%min                = ', p%x%min
    nl=nl+1; write (lines(nl), imt) 'x%major_div          = ', p%x%major_div
    nl=nl+1; write (lines(nl), imt) 'x%places             = ', p%x%places
    nl=nl+1; write (lines(nl), lmt) 'x%draw_label         = ', p%x%draw_label
    nl=nl+1; write (lines(nl), lmt) 'x%draw_numbers       = ', p%x%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'independent_graphs   = ', p%independent_graphs
    
    nl=nl+1; write (lines(nl), *) 'Graphs:'
    do i = 1, size(p%graph)
      nl=nl+1; write (lines(nl), amt) '   ', p%graph(i)%name
    enddo

  else
    call out_io (s_error$, r_name, 'This is not a graph')
    return
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! top10

case ('top10')

  call tao_top10_print ()

!----------------------------------------------------------------------
! universe
    
case ('universe')

  if (word(1) == ' ') then
    ix_u = s%global%u_view
  else
    read (word(1), *, iostat = ios) ix_u
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER')
      return
    endif
    if (ix_u < 1 .or. ix_u > size(s%u)) then
      call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE')
      return
    endif
  endif

  u => s%u(ix_u)

  nl = 0
  nl=nl+1; write(lines(nl), imt) '%ix_uni                = ', u%ix_uni
  nl=nl+1; write(lines(nl), imt) '%n_d2_data_used        = ', u%n_d2_data_used
  nl=nl+1; write(lines(nl), imt) '%n_data_used           = ', u%n_data_used
  nl=nl+1; write(lines(nl), lmt) '%do_synch_rad_int_calc = ', u%do_synch_rad_int_calc
  nl=nl+1; write(lines(nl), lmt) '%do_chrom_calc         = ', u%do_chrom_calc
  nl=nl+1; write(lines(nl), lmt) '%calc_beam_emittance   = ', u%calc_beam_emittance
  nl=nl+1; write(lines(nl), lmt) '%is_on                 = ', u%is_on
  nl=nl+1; write(lines(nl), amt) '%beam_init_file        = ', trim(u%beam_init_file)

  call out_io (s_blank$, r_name, lines(1:nl)) 

!----------------------------------------------------------------------
! variable
    
case ('var')

  if (.not. associated (s%v1_var)) then
    call out_io (s_error$, r_name, 'NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!')
    return 
  endif

! If 'n@' is present then write out stuff for universe n

  ix = index(word(1), '@')
  if (ix /= 0) then
    if (ix == 1) then
      ix_u = s%global%u_view
    else
      read (word(1)(:ix-1), *, iostat = ios) ix_u
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER')
        return
      endif
      if (ix_u == 0) ix_u = s%global%u_view
      if (ix_u < 1 .or. ix_u > size(s%u)) then
        call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE')
        return
      endif
    endif
    write (lines(1), '(a, i4)') 'Variables controlling universe:', ix_u
    write (lines(2), '(5x, a)') '                    '
    write (lines(3), '(5x, a)') 'Name                '
    nl = 3
    do i = 1, size(s%var)
      if (.not. s%var(i)%exists) cycle
      found = .false.
      do j = 1, size(s%var(i)%this)
        if (s%var(i)%this(j)%ix_uni == ix_u) found = .true.
      enddo
      if (.not. found) cycle
      nam = tao_var1_name(s%var(i))
      nl=nl+1; write(lines(nl), '(5x, a, a40)') nam(1:25), s%var(i)%name
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! If just "show var" then show all namees

  if (word(1) == '*') then
    call tao_var_write (' ')
    return
  endif

  if (word(1) == ' ') then
    write (lines(1), '(5x, a)') '                      Bounds'
    write (lines(2), '(5x, a)') 'Name                Lower  Upper'
    nl = 2
    do i = 1, size(s%v1_var)
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
      nl=nl+1; write(lines(nl), '(5x, a20, i5, i7)') v1_ptr%name, &
                                       lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the variables

  call string_trim (word(2), word(2), ix)
! are we looking at a range of locations?

  call tao_find_var(err, word(1), v1_ptr, v_array) 
  if (err) return
  n_size = 0
  if (allocated(v_array)) n_size = size(v_array)

! v_ptr is valid then show the variable info.

  if (n_size == 1) then

    v_ptr => v_array(1)%v

    nl=nl+1; write(lines(nl), amt)  'Name          = ', v_ptr%name        
    nl=nl+1; write(lines(nl), amt)  'Alias         = ', v_ptr%alias       
    nl=nl+1; write(lines(nl), amt)  'Ele_name      = ', v_ptr%ele_name    
    nl=nl+1; write(lines(nl), amt)  'Attrib_name   = ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  'Ix_var        = ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  'Ix_dvar       = ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  'Ix_v1         = ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), rmt)  'Model_value   = ', v_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  'Base_value    = ', v_ptr%base_value

    if (.not. allocated (v_ptr%this)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      do i = 1, size(v_ptr%this)
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%Ix_uni:        ', &
                                                            v_ptr%this(i)%ix_uni
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%Ix_ele:        ', v_ptr%this(i)%ix_ele
        if (associated (v_ptr%this(i)%model_ptr)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_ptr:   ', &
                                                            v_ptr%this(i)%model_ptr
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_ptr:   <not associated>'
        endif
        if (associated (v_ptr%this(i)%base_ptr)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_ptr:    ', &
                                                            v_ptr%this(i)%base_ptr
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_ptr:    <not associated>'
        endif
      enddo
    endif

    nl=nl+1; write(lines(nl), rmt)  '%Design_value     = ', v_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%Old_value        = ', v_ptr%old_value
    nl=nl+1; write(lines(nl), rmt)  '%Meas_value       = ', v_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%Ref_value        = ', v_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%Correction_value = ', v_ptr%correction_value
    nl=nl+1; write(lines(nl), rmt)  '%High_lim         = ', v_ptr%high_lim
    nl=nl+1; write(lines(nl), rmt)  '%Low_lim          = ', v_ptr%low_lim
    nl=nl+1; write(lines(nl), rmt)  '%Step             = ', v_ptr%step
    nl=nl+1; write(lines(nl), rmt)  '%Weight           = ', v_ptr%weight
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit      = ', v_ptr%delta_merit
    nl=nl+1; write(lines(nl), amt)  '%Merit_type       = ', v_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%Merit            = ', v_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%dMerit_dVar      = ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), lmt)  '%Exists           = ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%Good_var         = ', v_ptr%good_var
    nl=nl+1; write(lines(nl), lmt)  '%Good_user        = ', v_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%Good_opt         = ', v_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_opt        = ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_plot       = ', v_ptr%useit_plot

! check if there is a variable number
! if no variable number requested, show a range

  elseif (associated(v1_ptr)) then

    nc = 0
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      nc = max(nc, len_trim(v_ptr%name))
    enddo

    write(lines(1), '(2a)') 'Variable name:  ', v1_ptr%name
    lines(2) = ' '
    line1 = '       Name'
    line1(nc+17:) = 'Meas         Model        Design  Useit_opt'
    write (lines(3), *) line1
    nl = 3
    ! if a range is specified, show the variable range   
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200)
      nl=nl+1
      write(lines(nl), '(i6, 2x, a)') v_ptr%ix_v1, v_ptr%name
      write(lines(nl)(nc+9:), '(3es14.4, 7x, l)') v_ptr%meas_value, &
                 v_ptr%model_value, v_ptr%design_value, v_ptr%useit_opt
    enddo
    nl=nl+1
    write (lines(nl), *) line1

  else
    lines(1) = '???'
    nl = 1
  endif

! print out results

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "INTERNAL ERROR, SHOULDN'T BE HERE!")
  return

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine show_ele_data (u, i_ele, lines, nl)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum
character(*) :: lines(:)
character(100) l1
integer i_ele, nl, i

logical found_one

!

nl=nl+1; write (lines(nl), '(a)') "  "
write (l1, '(a, 20x, a)') "Data Name", &
          "Data Type             |  Model Value  |  Design Value |  Base Value"
nl=nl+1; lines(nl) = l1

found_one = .false.
do i = 1, size(u%data)
  if (u%data(i)%ix_ele .eq. i_ele) then
    found_one = .true.
    datum => u%data(i)
    nl=nl+1; write (lines(nl), "(a, t30, a20, 3(1x, es15.5))") &
                trim(tao_datum_name(datum)),  datum%data_type, datum%model_value, &
                datum%design_value, datum%base_value 
  endif
enddo

if (found_one) then
  nl=nl+1; lines(nl) = l1
else
  nl=nl+1; write (lines(nl), '(a)') "No data associated with this element."
endif

end subroutine show_ele_data

end subroutine tao_show_cmd

end module
