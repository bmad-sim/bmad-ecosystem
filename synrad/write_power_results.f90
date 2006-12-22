!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_results (wall, ring, gen_params)

  use sr_struct
  use sr_interface

  implicit none

  type (wall_struct), target :: wall
  type (wall_seg_struct), pointer :: seg
  type (sr_power_struct), pointer :: ep
  type (general_param_struct) gen_params
  type (ring_struct) ring

  integer i
  character*16 seg_name, ep_source_name, e_source_name, p_source_name
  character fmt*80, ep_name*2
  character file1*40, file2*40, file3*40

!

  file1 = 'synch_power_' // trim(wall_name(wall%side)) // '.dat'
  call write_power_header (1, file1, gen_params)


  do i = 1, wall%n_seg_tot
    seg => wall%seg(i)

    ep => seg%sr_power

    if (ep%ix_ele_source == 0) then
      ep_source_name = '--------'
      ep_name = '--'
    else
      ep_source_name = ring%ele_(ep%ix_ele_source)%name
    endif

    call convert_blanks_to_underscore (wall%pt(seg%ix_pt)%name, seg_name)

    fmt = '(i6,1x,a10,f8.3,f8.3,e12.2,e12.4,e12.4,e12.4)'
    write (1, fmt) &
              i, seg_name, seg%s, seg%x, &
              ep%power_per_len, &
              1.e-6 * (ep%power_per_area), &
              ep%power, ep%photons_per_sec

  enddo
  close (unit = 1)
  type *, 'Written: ', file1

end subroutine write_power_results
