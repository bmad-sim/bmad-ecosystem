module synrad_write_power_mod

use synrad_mod

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine write_power_header (iu, file, gen_params, synrad_mode)

  implicit none

  type (synrad_param_struct) gen_params
  type (synrad_mode_struct) synrad_mode
  character(*) file
  integer iu

!

  open (unit = iu, file = file)
  write (iu, *) 'Lattice: ', gen_params%lat_file
  write (iu, *) 'I_beam    =', gen_params%i_beam,    ' ! Amps/beam'
  write (iu, *) 'Epsilon_y =', gen_params%epsilon_y, ' ! Vertical emittance'

  write (iu, '(3(/,2x,a))') &
'          Segment                                  ', &
'  Ix  Name          S_seg      X_seg     P/len      P/Area     P_tot     Phot/sec      A Beta    B Beta    A Eta     Ele Type       Relevant               Ele',&
'                     (m)        (m)      (W/m)     (W/mm^2)      (W)      (1/s)         (m)        (m)      (m)      at s_mid       Attribute              Name'
!' S_seg      P/len      P/Area     Phot/sec',&
!'  (m)       (W/m)     (W/mm^2)     (1/s)'

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_results (wall, branch, gen_params, ix_ele1, ix_ele2, synrad_mode)

  implicit none

  type (wall_struct), target :: wall
  type (wall_seg_struct), pointer :: seg
  type (seg_power_struct), pointer :: ep
  type (seg_power_struct), target :: zero_power
  type (synrad_param_struct) gen_params
  type (synrad_mode_struct) :: synrad_mode
  type (branch_struct) branch

  real(rp) value

  integer ix_ele1, ix_ele2, i
  integer key

  character*16 seg_name, ep_source_name, e_source_name, p_source_name
  character fmt*100, ep_name*2
  character file1*50, ele_num*6
  character ele_at_seg*16, attrib*10

!

  call set_wall_eles (wall, branch)

  file1 = 'synch_power_' // trim(wall_name(wall%side))
  if (ix_ele1 == ix_ele2) then
    write (ele_num, '(i0)') ix_ele1
    file1 = trim(file1) // '_' // trim(ele_num)
  endif
  file1 = trim(file1) // '.dat'
  call downcase_string (file1)
  call write_power_header (1, file1, gen_params, synrad_mode)


  do i = 1, wall%n_seg_max
    seg => wall%seg(i)
    key = branch%ele(seg%ix_ele)%key
    attrib = ' '
    value = 0

    if (key == quadrupole$) then
      attrib = 'k1 = '
      value = branch%ele(seg%ix_ele)%value(k1$)
    elseif (key == sol_quad$) then
      attrib = 'k1 = '
      value = branch%ele(seg%ix_ele)%value(k1$)
    elseif (key == solenoid$) then
      attrib = 'ks = '
      value = branch%ele(seg%ix_ele)%value(ks$)
    elseif (key == sbend$) then
      attrib = 'G = '
      value = branch%ele(seg%ix_ele)%value(g$)
    elseif (key == rbend$) then
      attrib = 'G = '
      value = branch%ele(seg%ix_ele)%value(g$)
    elseif (key == sextupole$) then
      attrib = 'k2 = '
      value = branch%ele(seg%ix_ele)%value(k2$)
    elseif (key == wiggler$) then
      attrib = 'B_max = '
      value = branch%ele(seg%ix_ele)%value(b_max$)
!      call type_ele(branch%ele(seg%ix_ele))
    end if

    if (gen_params%filter_phantom_photons .and. wall%pt(seg%ix_pt)%phantom) then
      ep => zero_power
    else
      ep => seg%power
    endif

    if (ep%main_source%ix_ele == 0) then
      ep_source_name = '--------'
      ep_name = '--'
    else
      ep_source_name = branch%ele(ep%main_source%ix_ele)%name
    endif

    seg_name = wall%pt(seg%ix_pt)%name
    call str_substitute (seg_name, ' ', '_', .true.)

    fmt = '(i6, 1x, a10, f10.4, 2es11.3, 3es12.4, 3f10.3, 2x, a16, 1x, a10, es12.4, 1x, a)'
    write (1, fmt) &
              i, seg_name, seg%s, seg%x, &
              ep%power_per_len, &
              1.e-6 * (ep%power_per_area), &
              ep%power_tot, ep%photons_per_sec, &
              seg%a%beta, seg%b%beta,seg%a%eta, key_name(key), &
              attrib, value, trim(branch%ele(seg%ix_ele)%name)

!    fmt = '(f10.4, 3es11.3)'
!    write (1, fmt) &
!              seg%s, &
!              ep%power_per_len, &
!              1.e-6 * (ep%power_per_area), &
!              ep%photons_per_sec

  enddo
  close (unit = 1)
  print *, 'Written: ', file1

end subroutine write_power_results

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine write_header (iu, file, gen_params, synrad_mode)

  implicit none

  type (synrad_param_struct) gen_params
  type (synrad_mode_struct) synrad_mode
  character(*) file
  integer iu

!

  open (unit = iu, file = file)
  write (iu, *) 'Lattice: ', gen_params%lat_file
  write (iu, *) 'I_beam    =', gen_params%i_beam,    ' ! Amps/beam'
  write (iu, *) 'Input eps_y =', 1.e6*gen_params%epsilon_y, ' ! mm-mrad'
  write (iu, *) 'filter_phantom_photons =', gen_params%filter_phantom_photons
  write (iu, *) 'Positron eps_x =', 1.e6*synrad_mode%pos_mode%a%emittance, ' ! mm-mrad'
  write (iu, *) 'Positron eps_y =', 1.e6*synrad_mode%pos_mode%b%emittance, ' ! mm-mrad'
  write (iu, *) 'Positron sig_z =', synrad_mode%pos_mode%sig_z, ' ! m'
  write (iu, *) 'Positron sigE/E =',1.e2*synrad_mode%pos_mode%sigE_E, ' ! %'

  write (iu, *) 'Electron eps_x =', 1.e6*synrad_mode%ele_mode%a%emittance, ' ! mm-mrad'
  write (iu, *) 'Electron eps_y =', 1.e6*synrad_mode%ele_mode%b%emittance, ' ! mm-mrad'
  write (iu, *) 'Electron sig_z =', synrad_mode%ele_mode%sig_z, ' ! m'
  write (iu, *) 'Electron sigE/E =',1.e2*synrad_mode%ele_mode%sigE_E, ' ! %'

  write (iu, *)
  write (iu, '(a)') 'Attribute values are G (bends), k1 (quads), k2 (sextupoles), ks (solenoids), B_max (wigglers)'

  write (iu, '(3(/,2x,a))') &
'Segment                                  ', &
'  Nr      S_seg     L_seg      X_seg     P/len      P/Area     P_tot     Phot/sec      A Beta    B Beta    A Eta     Ele Type     Attribute       Element       Source          Primary        S_source    Nr ',&
'           (m)       (m)        (m)      (W/m)     (W/mm^2)      (W)      (1/s)         (m)        (m)      (m)      at s_mid       Value          Name        Ele Type         Source           (m)     Sources'

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_results (wall, branch, gen_params, ix_ele1, ix_ele2, synrad_mode)

  implicit none

  type (wall_struct), target :: wall
  type (wall_seg_struct), pointer :: seg
  type (seg_power_struct), pointer :: ep
  type (synrad_param_struct) gen_params
  type (synrad_mode_struct) :: synrad_mode
  type (branch_struct) branch
  type (seg_power_struct), target :: zero_power

  integer ix_ele1, ix_ele2, i
  character*16 seg_name, ep_source_name, ep_source_key_name
  character fmt*130, ep_name*2
  character file1*50, ele_num*6
  character ele_at_seg*16, attrib*10
  integer key, ep_source_key
  real(rp) value

!

  call set_wall_eles (wall, branch)

  file1 = 'synrad_' // trim(wall_name(wall%side))
  if (ix_ele1 == ix_ele2) then
    write (ele_num, '(i0)') ix_ele1
    file1 = trim(file1) // '_' // trim(ele_num)
  endif
  file1 = trim(file1) // '.txt'
  call downcase_string (file1)
  call write_header (1, file1, gen_params, synrad_mode)


  do i = 1, wall%n_seg_max
    seg => wall%seg(i)
    key = branch%ele(seg%ix_ele)%key
    attrib = ' '
    value = 0

    if (key == quadrupole$) then
      attrib = 'k1 = '
      value = branch%ele(seg%ix_ele)%value(k1$)
    elseif (key == sol_quad$) then
      attrib = 'k1 = '
      value = branch%ele(seg%ix_ele)%value(k1$)
    elseif (key == solenoid$) then
      attrib = 'ks = '
      value = branch%ele(seg%ix_ele)%value(ks$)
    elseif (key == sbend$) then
      attrib = 'G = '
      value = branch%ele(seg%ix_ele)%value(g$)
    elseif (key == rbend$) then
      attrib = 'G = '
      value = branch%ele(seg%ix_ele)%value(g$)
    elseif (key == sextupole$) then
      attrib = 'k2 = '
      value = branch%ele(seg%ix_ele)%value(k2$)
    elseif (key == wiggler$) then
      attrib = 'B_max = '
      value = branch%ele(seg%ix_ele)%value(b_max$)
!      call type_ele(branch%ele(seg%ix_ele))
    end if

    if (gen_params%filter_phantom_photons .and. wall%pt(seg%ix_pt)%phantom) then
      ep => zero_power
    else
      ep => seg%power
    endif

    if (ep%main_source%ix_ele == 0) then
      ep_source_name = '------------------'
      ep_name = '--'
      ep_source_key = 0
      ep_source_key_name = '----------'
    else
      ep_source_name = branch%ele(ep%main_source%ix_ele)%name
      ep_source_key = branch%ele(ep%main_source%ix_ele)%key
      ep_source_key_name = key_name(ep_source_key)
    endif

    seg_name = wall%pt(seg%ix_pt)%name
    call str_substitute (seg_name, ' ', '_', .true.)

    fmt = '(i6, 1x, f10.4, f10.4, 2es11.3, 3es12.4, 3f10.3, 2x, a10, 1x, es12.4, 1x, a18, 1x, a10, 1x, a18, f11.4, i6)'
    write (1, fmt) &
              i, seg%s, seg%len, seg%x, &
              ep%power_per_len, &
              1.e-6 * (ep%power_per_area), &
              ep%power_tot, ep%photons_per_sec, &
              seg%a%beta, seg%b%beta,seg%a%eta, key_name(key), &
              value, trim(branch%ele(seg%ix_ele)%name), &
              ep_source_key_name, trim(ep_source_name), ep%main_source%s, ep%n_source

  enddo
  close (unit = 1)
  print *, 'Written: ', file1

end subroutine write_results

end module
