module io_mod

  use bmad_struct
  use bmad_interface

  private str, rchomp, write_out

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_bmad_lattice_file (lattice_name, ring)
!
! Subroutine to write a Bmad lattice file using the information in
! a ring_struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice_name -- Character(*): Name of the lattice file
!   ring         -- Ring_struct: Holds the lattice information.
!-

subroutine write_bmad_lattice_file (lattice_name, ring)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (wig_term_struct) wt
  type (control_struct) ctl
  type (taylor_term_struct) tm

  real(rp) s0

  character(*) lattice_name
  character(4000) line
  character(4) last
  character(16) name

  integer i, j, k, ix, iu, ios, ixs
  integer unit(6)

  logical slave_here, unit_found, write_term

! Open the file

  iu = lunget()
  open (iu, file = lattice_name, iostat = ios)
  if (ios /= 0) then
    print *, 'ERROR IN WRITE_BMAD_LATTICE_FILE: CANNOT OPEN FILE: ', &
                                                          trim(lattice_name)
    return
  endif

! Non-elemental stuff

  if (ring%title /= ' ') &
            write (iu, *) 'title, "', trim(ring%title), '"'
  if (ring%lattice /= ' ') &
            write (iu, *) 'parameter[lattice] = "', trim(ring%lattice), '"'
  write (iu, *) 'parameter[lattice_type] = ', lattice_type(ring%param%lattice_type)
  write (iu, *) 'parameter[taylor_order] =', ring%input_taylor_order

  write (iu, *)
  write (iu, *) 'parameter[beam_energy] =', &
                      trim(str(ring%ele_(0)%value(beam_energy$)))
  write (iu, *) 'beam, particle = ', particle_name(ring%param%particle)
  if (ring%param%n_part /= 0) &
          write (iu, *) 'beam, n_part = ', particle_name(ring%param%n_part)

  ele => ring%ele_(0) 

  if (ele%position%x /= 0) &
        write (iu, *) 'beginning[x_position] = ', trim(str(ele%position%x))
  if (ele%position%y /= 0) &
        write (iu, *) 'beginning[y_position] = ', trim(str(ele%position%y))
  if (ele%position%z /= 0) &
        write (iu, *) 'beginning[z_position] = ', trim(str(ele%position%z))
  if (ele%position%theta /= 0) &
    write (iu, *) 'beginning[theta_position] = ', trim(str(ele%position%theta))
  if (ele%position%phi /= 0) &
        write (iu, *) 'beginning[phi_position] = ', trim(str(ele%position%phi))
  if (ele%position%psi /= 0) &
        write (iu, *) 'beginning[psi_position] = ', trim(str(ele%position%psi))

  if (ring%param%lattice_type /= circular_lattice$) then
    write (iu, *)
    if (ele%x%beta /= 0) &
        write (iu, *) 'beginning[beta_x] = ', trim(str(ele%x%beta))
    if (ele%x%alpha /= 0) &
        write (iu, *) 'beginning[alpha_x] = ', trim(str(ele%x%alpha))
    if (ele%x%phi /= 0) &
        write (iu, *) 'beginning[phi_x] = ', trim(str(ele%x%phi))
    if (ele%x%eta /= 0) &
        write (iu, *) 'beginning[eta_x] = ', trim(str(ele%x%eta))
    if (ele%x%etap /= 0) &
        write (iu, *) 'beginning[etap_x] = ', trim(str(ele%x%etap))
    if (ele%y%beta /= 0) &
        write (iu, *) 'beginning[beta_y] = ', trim(str(ele%y%beta))
    if (ele%y%alpha /= 0) &
        write (iu, *) 'beginning[alpha_y] = ', trim(str(ele%y%alpha))
    if (ele%y%phi /= 0) &
        write (iu, *) 'beginning[phi_y] = ', trim(str(ele%y%phi))
    if (ele%y%eta /= 0) &
        write (iu, *) 'beginning[eta_y] = ', trim(str(ele%y%eta))
    if (ele%y%etap /= 0) &
        write (iu, *) 'beginning[etap_y] = ', trim(str(ele%y%etap))
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
  write (iu, *) '!-------------------------------------------------------'
  write (iu, *)

  slave_here = .false.
  ixs = 0

  do i = 1, ring%n_ele_max

    ele => ring%ele_(i)

    if (i == ring%n_ele_ring+1) then
      write (iu, *)
      write (iu, *) '!-------------------------------------------------------'
      write (iu, *)
    endif

    if (ele%control_type == super_slave$) then
      ele%ixx = 0
      if (slave_here) cycle
      s0 = ring%ele_(i-1)%s
      ixs = ixs + 1
      ele%ixx = ixs
      slave_here = .true.
      cycle
    endif

    if (slave_here) then
      s0 = ring%ele_(i-1)%s - s0
      write (iu, '(a, i3.3, 2a)') ' slave_drift_', ixs, ': drift, l = ', trim(str(s0))
      slave_here = .false.
    endif

! Overlays and groups

    if (ele%control_type == overlay_lord$ .or. ele%control_type == group_lord$) then
      if (ele%control_type == overlay_lord$) then
        write (line, '(2a)') trim(ele%name), ': overlay = {'
      else
        write (line, '(2a)') trim(ele%name), ': group = {'
      endif
      do j = ele%ix1_slave, ele%ix2_slave
        ctl = ring%control_(j)
        ix = ctl%ix_slave
        if (j == ele%ix1_slave) then
          write (line, '(3a)') trim(line), trim(ring%ele_(ix)%name)
        else
          write (line, '(3a)') trim(line), ', ', trim(ring%ele_(ix)%name)
        endif
        name = attribute_name(ring%ele_(ix), ctl%ix_attrib)  
        if (name /= ele%attribute_name) &
                line = trim(line) // '[' // trim(name) // ']'
        if (ctl%coef /= 1) write (line, '(3a)') trim(line), '/', trim(str(ctl%coef))
      enddo
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

!

    line = trim(ele%name) // ': ' // key_name(ele%key)
    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'

    if (ele%control_type == super_lord$) then
      line = trim(line) // ', superimpose, offset = ' // &
                                             trim(str(ele%s-ele%value(l$)/2))
    endif

    if (associated(ele%descrip)) line = trim(line) // &
                              ', descrip = "' // trim(ele%descrip) // '"'
    if (associated(ele%wake%sr_file)) line = trim(line) // &
                              ',  sr_file = "' // trim(ele%wake%sr_file) // '"'
    if (associated(ele%wake%lr_file)) line = trim(line) // &
                              ',  lr_file = "' // trim(ele%wake%lr_file) // '"'

    do j = 1, n_attrib_maxx

      if (j == check_sum$ .and. ele%key /= patch$) cycle
      if (j == beam_energy$) cycle

      select case (ele%key)
      case (beambeam$)
        if (j == bbi_const$) cycle
      case (elseparator$)
        if (j == e_field$) cycle
        if (j == voltage$) cycle
      case (lcavity$)
        if (j == e_loss$) cycle
        if (j == delta_e$) cycle
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
        case (sbend$)
          if (j == g$) cycle
        end select
      else
        if (j == b_field$) cycle
        if (j == e_field$) cycle
      endif
      
      if (ele%value(j) == 0) cycle
      line = trim(line) // ', ' // trim(attribute_name(ele, j)) // &
                                                  ' = ' // str(ele%value(j))

    enddo

    if (ele%mat6_calc_method /= bmad_standard$) line = trim(line) // &
            ', mat6_calc_method = ' // calc_method_name(ele%mat6_calc_method)
    if (ele%tracking_method /= bmad_standard$) line = trim(line) // &
            ', tracking_method = ' // calc_method_name(ele%tracking_method)
    if (ele%num_steps > 1) write (line, '(2a, i3)') trim(line), &
            ', num_steps =', ele%num_steps
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

    if (associated(ele%a)) then
      do j = 0, ubound(ele%a, 1)
        if (ele%a(j) /= 0) write (line, '(2a, i2.2, 2a)') &
                    trim(line), ', a', j, ' = ', trim(str(ele%a(j)))
        if (ele%b(j) /= 0) write (line, '(2a, i2.2, 2a)') &
                    trim(line), ', b', j, ' = ', trim(str(ele%b(j)))
      enddo
    endif
    
    if (associated(ele%wig_term)) then
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

  enddo


! Lattice Layout

  write (iu, *)
  write (iu, *) '!-------------------------------------------------------'
  write (iu, *)

  line = 'main_line: line = ('

  do i = 1, ring%n_ele_ring
    ele => ring%ele_(i)
    if (ele%control_type == super_slave$) then
      if (ele%ixx == 0) cycle
      write (line, '(2a, i3.3, a)') trim(line), ' slave_drift_', ele%ixx, ','
    else
      write (line, '(4a)') trim(line), ' ', trim(ele%name), ','
    endif
    if (mod(i, 20) == 0) call write_out(line, iu, .false.)
  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_out (line, iu, .true.)

  write (iu, *)
  write (iu, *) 'use, main_line'

  close(iu)

end subroutine

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
    write (str_out, fmt) trim(rchomp(rel/10**pl, 0)), 'E', pl

  elseif (pl > -3) then
    str_out = rchomp(rel, pl)

  else
    fmt = '(2a, i2)'
    if (pl < -9)  fmt = '(2a, i3)'
    write (str_out, fmt) trim(rchomp(rel*10**(-pl), 0)), 'E', pl

  endif

end function

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

end module
