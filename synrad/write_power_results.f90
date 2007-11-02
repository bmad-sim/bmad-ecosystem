!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_results (wall, lat, gen_params, use_ele_ix)

  use sr_struct
  use sr_interface

  implicit none

  type (wall_struct), target :: wall
  type (wall_seg_struct), pointer :: seg
  type (sr_power_struct), pointer :: ep
  type (synrad_param_struct) gen_params
  type (lat_struct) lat

  integer use_ele_ix, i
  character*16 seg_name, ep_source_name, e_source_name, p_source_name
  character fmt*80, ep_name*2
  character file1*50, ele_num*6

!

  file1 = 'synch_power_' // trim(wall_name(wall%side))
  if (use_ele_ix .ne. 0) then
    write (ele_num, '(i6.6)') use_ele_ix
    file1 = trim(file1) // '_' // trim(ele_num)
  endif
  file1 = trim(file1) // '.dat'
  call downcase_string (file1)
  call write_power_header (1, file1, gen_params)


  do i = 1, wall%n_seg_tot
    seg => wall%seg(i)

    ep => seg%sr_power

    if (ep%ix_ele_source == 0) then
      ep_source_name = '--------'
      ep_name = '--'
    else
      ep_source_name = lat%ele(ep%ix_ele_source)%name
    endif

    call convert_blanks_to_underscore (wall%pt(seg%ix_pt)%name, seg_name)

    fmt = '(i6,1x,a10,f10.3,e11.3,e11.3,e12.4,e12.4,e12.4)'
    write (1, fmt) &
              i, seg_name, seg%s, seg%x, &
              ep%power_per_len, &
              1.e-6 * (ep%power_per_area), &
              ep%power, ep%photons_per_sec

  enddo
  close (unit = 1)
  type *, 'Written: ', file1

end subroutine write_power_results
