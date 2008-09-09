module synrad_output_mod

use sr_mod
use cesrv_struct
use cesrv_interface
use synrad_window_mod

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Types out the crotch window hit information

subroutine ray_output ( window, ring )

  implicit none

  type (crotch_window_struct), target :: window(:)
  type (lat_struct) ring

  integer i, j, lun, lunget, iw(n_windows$)
  real(rp), pointer :: sigma
  type (coord_struct), pointer :: coord
  character*16 line
  logical target

  call get_window_numbers( window, iw )

  print *, 'Do you want ray information at the projected target '
  call get_input_string ('or at the window (default) ?  (Enter t or w)', line)
  call str_upcase(line, line)
  line = adjustl(line)
  target = .false.
  if (line(1:1) == 'T') target = .true.

  lun = lunget()
  open (lun, file = 'windows.out')

  if (target) then
    print *, 'Crotch Window Data'
    write (lun, '(1x, a)') 'Crotch Window Data'
  else
    print *, 'Crotch Window Data'
    write (lun, '(1x, a)') 'Crotch Window Data'
  endif

  write (lun, '(1x, 6a13)') 'x pos   ', 'y pos   ', 'sigma_y   ', &
           'dx/ds   ', 'dy/ds   ', 'Source Ele '
  write (*, '(1x, 6a13)') 'x pos   ', 'y pos   ', 'sigma_y   ', &
           'dx/ds   ', 'dy/ds   ', 'Source Ele '

  do i=1,n_windows$

    if (iw(i) == 0) exit
    write (lun, '(1x a)') window(iw(i))%name
    print *, window(iw(i))%name," has ",window(iw(i))%n_ray_hit,"ray hits."
    do j=1,window(iw(i))%n_ray_hit
      
      if (target) then
        sigma => window(iw(i))%ray_hit_(j)%sig_y_eff 
        coord => window(iw(i))%ray_hit_(j)%target_coord
      else
        sigma => window(iw(i))%ray_hit_(j)%window_sig_y
        coord => window(iw(i))%ray_hit_(j)%hit_coord
      endif

      write (*, '(1x, 5f13.9, 1x, a16)') coord%vec(1), coord%vec(3), &
            sigma, coord%vec(2), coord%vec(4), &
            ring%ele(window(iw(i))%ray_hit_(j)%ray%ix_source)%name
      write (lun, '(1x, 5f13.9, 1x, a16)') coord%vec(1), coord%vec(3), &
            sigma, coord%vec(2), coord%vec(4), &
            ring%ele(window(iw(i))%ray_hit_(j)%ray%ix_source)%name
    
    enddo

    write (lun, '(1x)')

  enddo
             
  close (lun)
  print *, "Windows.out written"

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine sextupole_output ( u )

  implicit none

  type (universe_struct), target :: u
  type (var_struct), pointer :: avar

  integer i, lun, lunget, ix
  real(rp) kl, bxtotal, bytotal, bxlocal, bylocal

  !

  bxtotal = 0
  bytotal = 0

  lun = lunget()
  open (lun, file = 'sexts.out')
  write (lun, '(1x, 7a12)') 'ix_ele  ', 'name   ', 'kl      ', &
           'bxlocal   ', 'bxtotal  ', 'bylocal  ', 'bytotal  '
  write (*, '(1x, 7a12)') 'ix_ele  ', 'name   ', 'kl      ', &
           'bxlocal   ', 'bxtotal  ', 'bylocal  ', 'bytotal  '

  do i = lbound(u%sex_k2%v, 1), ubound(u%sex_k2%v,1)
    avar => u%sex_k2%v(i)
    if (.not. avar%exists) cycle
    ix = avar%ix_ele
    kl = avar%model * u%ring%ele(ix)%value(l$)
    bylocal = kl / 2 * ( u%orb(ix)%vec(1)**2 + u%orb(ix)%vec(3)**2)
    bxlocal = kl * ( u%orb(ix)%vec(1) * u%orb(ix)%vec(3))
    bxtotal = bxtotal + abs(bxlocal)
    bytotal = bytotal + abs(bylocal)
    write (lun, '(1x,i6,a16,5ES11.3)') avar%ix_ele, trim(avar%name), kl, bxlocal,&
             bxtotal, bylocal, bytotal
    write (*, '(1x,i6,a16,5ES11.3)') avar%ix_ele, trim(avar%name), kl, bxlocal,&
             bxtotal, bylocal, bytotal
  enddo
  close (lun)
  print *, "sexts.out written"

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine write_power_header (iu, file, gen_params)

  implicit none

  type (synrad_param_struct) gen_params
  character(*) file
  integer iu

!

  open (unit = iu, file = file)
  write (iu, *) 'Lattice: ', gen_params%lat_file
  write (iu, *) 'I_beam    =', gen_params%i_beam,    ' ! Amps/beam'
  write (iu, *) 'Epsilon_y =', gen_params%epsilon_y, ' ! Vertical emittance'

  write (iu, '(3(/,2x,a))') &
'          Segment                                  ', &
'  Ix  Name          S_seg      X_seg     P/len      P/Area     P_tot     Phot/sec      A Beta    B Beta    Ele Type       Relevant               Ele',&
'                     (m)        (m)      (W/m)     (W/mm^2)      (W)      (1/s)         (m)        (m)     at s_mid       Attribute              Name'

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_results (wall, lat, gen_params, use_ele_ix)

  implicit none

  type (wall_struct), target :: wall
  type (wall_seg_struct), pointer :: seg
  type (sr_power_struct), pointer :: ep
  type (synrad_param_struct) gen_params
  type (lat_struct) lat

  integer use_ele_ix, i
  character*16 seg_name, ep_source_name, e_source_name, p_source_name
  character fmt*100, ep_name*2
  character file1*50, ele_num*6
  character ele_at_seg*16, attrib*10
  integer key
  real(rp) value

!

  call set_wall_eles (wall, lat)

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
    key = lat%ele(seg%ix_ele)%key
    attrib = ' '
    value = 0

    if (key == quadrupole$) then
      attrib = 'k1 = '
      value = lat%ele(seg%ix_ele)%value(k1$)
    elseif (key == sol_quad$) then
      attrib = 'k1 = '
      value = lat%ele(seg%ix_ele)%value(k1$)
    elseif (key == solenoid$) then
      attrib = 'ks = '
      value = lat%ele(seg%ix_ele)%value(ks$)
    elseif (key == sbend$) then
      attrib = 'G = '
      value = lat%ele(seg%ix_ele)%value(g$)
    elseif (key == rbend$) then
      attrib = 'G = '
      value = lat%ele(seg%ix_ele)%value(g$)
    elseif (key == sextupole$) then
      attrib = 'k2 = '
      value = lat%ele(seg%ix_ele)%value(k2$)
    elseif (key == wiggler$) then
      attrib = 'B_max = '
      value = lat%ele(seg%ix_ele)%value(b_max$)
!      call type_ele(lat%ele(seg%ix_ele))
    end if

    ep => seg%sr_power

    if (ep%ix_ele_source == 0) then
      ep_source_name = '--------'
      ep_name = '--'
    else
      ep_source_name = lat%ele(ep%ix_ele_source)%name
    endif

    call convert_blanks_to_underscore (wall%pt(seg%ix_pt)%name, seg_name)

    fmt = '(i6,1x,a10,f10.3,e11.3,e11.3,e12.4,e12.4,e12.4,f10.3,f10.3,2x,a16,1x,a10,e12.4,1x,a40)'
    write (1, fmt) &
              i, seg_name, seg%s, seg%x, &
              ep%power_per_len, &
              1.e-6 * (ep%power_per_area), &
              ep%power, ep%photons_per_sec, &
              seg%a%beta, seg%b%beta, key_name(key), &
              attrib, value, lat%ele(seg%ix_ele)%name

  enddo
  close (unit = 1)
  type *, 'Written: ', file1

end subroutine write_power_results

end module
