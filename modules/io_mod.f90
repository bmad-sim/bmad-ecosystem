module io_mod

  use bmad_struct
  use bmad_interface
  use multipole_mod

  private str, rchomp, write_out

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_bmad_lattice_file (bmad_file, ring)
!
! Subroutine to write a Bmad lattice file using the information in
! a ring_struct. Optionally only part of the lattice can be generated.
!
! Modules needed:
!   use io_mod
!
! Input:
!   bmad_file     -- Character(*): Name of the output lattice file.
!   ring          -- Ring_struct: Holds the lattice information.
!   ix_start      -- Integer, optional: Starting index of ring%ele_(i)
!                       used for output.
!   ix_end        -- Integer, optional: Ending index of ring%ele_(i)
!                       used for output.
!-

subroutine write_bmad_lattice_file (bmad_file, ring)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele, slave
  type (ele_struct) ele_init
  type (wig_term_struct) wt
  type (control_struct) ctl
  type (taylor_term_struct) tm

  real(rp) s0

  character(*) bmad_file
  character(4000) line
  character(4) last
  character(16) name
  character(16), allocatable :: names(:)
 
  integer i, j, k, ix, iu, ios, ixs, ix1
  integer unit(6), ix_names, ix_match

  logical slave_here, unit_found, write_term, match_found

! Init

  call init_ele (ele_init)

! Open the file

  iu = lunget()
  open (iu, file = bmad_file, iostat = ios)
  if (ios /= 0) then
    print *, 'ERROR IN WRITE_BMAD_LATTICE_FILE: CANNOT OPEN FILE: ', &
                                                          trim(bmad_file)
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
          write (iu, *) 'beam, n_part = ', ring%param%n_part

  ele => ring%ele_(0) 

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
  ix_names = 0
  allocate (names(ring%n_ele_max))

  ele_loop: do i = 1, ring%n_ele_max

    ele => ring%ele_(i)

    if (ele%key == null_ele$) cycle

    if (i == ring%n_ele_use+1) then
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
        ctl = ring%control_(j)
        ix = ctl%ix_slave
        slave => ring%ele_(ix)
        do k = ele%ix1_slave, j-1 ! do not use elements w/ duplicate names
          if (ring%ele_(ring%control_(k)%ix_slave)%name == slave%name) exit j_loop
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

! I_beam

    if (ele%control_type == i_beam$) then
      write (line, '(2a)') trim(ele%name), ': i_beam = {'
      do j = ele%ix1_slave, ele%ix2_slave
        ix1 = ring%control_(j)%ix_slave
        if (j == ele%ix2_slave) then
          write (line, '(3a)') trim(line), trim(ring%ele_(ix1)%name), '}'
        else
          write (line, '(3a)') trim(line), trim(ring%ele_(ix1)%name), ', '
        endif
      enddo
    else
      line = trim(ele%name) // ': ' // key_name(ele%key)
    endif

! other elements

    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'

    if (ele%control_type == super_lord$) then
      line = trim(line) // ', superimpose, offset = ' // &
                                             trim(str(ele%s-ele%value(l$)/2))
    endif

    if (associated(ele%descrip)) line = trim(line) // &
                              ', descrip = "' // trim(ele%descrip) // '"'

    if (associated(ele%wake)) then
      if (ele%wake%sr_file /= ' ') line = &
            trim(line) // ',  sr_file = "' // trim(ele%wake%sr_file) // '"'
      if (ele%wake%lr_file /= ' ') line = &
            trim(line) // ',  lr_file = "' // trim(ele%wake%lr_file) // '"'
    endif

    do j = 1, n_attrib_maxx

! Exclude dependent variables

      if (j == check_sum$ .and. ele%key /= patch$) cycle
      if (j == beam_energy$) cycle
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
        if (j == energy_start$) cycle
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
          if (j == b_gradient$) cycle
        case (sextupole$)
          if (j == b_gradient$) cycle
        case (octupole$)
          if (j == b_gradient$) cycle
        case (solenoid$)
          if (j == b_field$) cycle
        case (sol_quad$) 
          if (j == b_field$) cycle
          if (j == b_gradient$) cycle
        case (sbend$)
          if (j == b_field$) cycle
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

      if (attribute_name(ele, j) == null_name) print *, 'Null: ', ele%name, j

    enddo

    if (ele%mat6_calc_method /= bmad_standard$) line = trim(line) // &
            ', mat6_calc_method = ' // calc_method_name(ele%mat6_calc_method)
    if (ele%tracking_method /= bmad_standard$) line = trim(line) // &
            ', tracking_method = ' // calc_method_name(ele%tracking_method)
    if (ele%num_steps /=ele_init%num_steps) write (line, '(2a, i3)') trim(line), &
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
        if (ele%a(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0$)) // ' = ' // str(ele%a(j))
        if (ele%b(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0$)) // ' = ' // str(ele%b(j))
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

  enddo ele_loop


! Lattice Layout

  write (iu, *)
  write (iu, *) '!-------------------------------------------------------'
  write (iu, *)

  line = 'main_line: line = ('

  do i = 1, ring%n_ele_use
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

! cleanup

  close(iu)
  deallocate (names)

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
!   name     -- Character(16): Name to match to.
!   names(:) -- Character(16): Array of sorted names.
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

  character(16) name, names(:)
  character(16) this_name

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
! Subroutine bmad_to_mad (mad_file, ring, ix_start, ix_end)
!
! Subroutine to write a mad lattice file using the information in
! a ring_struct. Optionally only part of the lattice can be generated.
!
! Modules needed:
!   use io_mod
!
! Input:
!   mad_file      -- Character(*): Name of the mad output lattice file.
!   ring          -- Ring_struct: Holds the lattice information.
!   ix_start      -- Integer, optional: Starting index of ring%ele_(i)
!                       used for output.
!   ix_end        -- Integer, optional: Ending index of ring%ele_(i)
!                       used for output.
!-

subroutine bmad_to_mad (mad_file, ring, ix_start, ix_end)

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele, drift_ele

  integer, optional :: ix_start, ix_end
  integer i, j, il, ix, iout, n, ix1, ix2, iu, ie1, ie2, n_list

  real(rp) field, hk, vk, tilt
  real(rp) knl(0:n_pole_maxx), tilt_(0:n_pole_maxx)

  character(*) mad_file 
  character(300) line
  character(6) line_name(3)
  character(16) name
  character(16), allocatable :: name_list(:)
  character(8) str

  logical init_needed, parsing

! open file

  iu = lunget()
  open (iu, file = mad_file)

! ring lattice name

  write (iu, '(a)')  '! File generated by BMAD_TO_MAD'
  write (iu, '(2a)') '! BMAD file: ', trim(mad_file)
  write (iu, '(2a)') '! Lattice: ', ring%lattice
  write (iu, *)

! beam definition

  ie1 = 1    ! element index
  if (present(ix_start)) ie1 = ix_start

  ie2 = ring%n_ele_use
  if (present(ix_end)) ie2 = ix_end

  allocate (name_list(ie2-ie1+1)) ! list of unique names
  n_list = 0                     ! number of names stored in the list


  ele = ring%ele_(ie1-1)

  write (iu, '(2a, 2(a, es12.5))')  &
        'BEAM, Particle = ', trim(particle_name(ring%param%particle)),  &
        ', Energy =', 1e-9*ele%value(beam_energy$), ', Npart =', ring%param%n_part

  write (iu, *)

! add drifts before and after wigglers and sol_quads so total length is invariant

  call init_ele (drift_ele)
  drift_ele%key = drift$

  j = 0    ! drift around solenoid or sol_quad index

  do i = ie1, ie2
    ele = ring%ele_(i)
    if (ele%key == wiggler$ .or. ele%key == sol_quad$) then
      j = j + 1
      write (drift_ele%name, '(a, i3.3)') 'DZ', j
      drift_ele%value(l$) = ele%value(l$) / 2
      call element_out (drift_ele)

      drift_ele%value(l$) = -drift_ele%value(l$)
      call make_mat6 (drift_ele, ring%param)
      ele%mat6 = matmul(matmul(drift_ele%mat6, ele%mat6), drift_ele%mat6)
      call element_out (ele)

      drift_ele%value(l$) = ele%value(l$) / 2
      call element_out (drift_ele)

    else
      call element_out (ele)
    endif
  enddo

  deallocate (name_list)

! Write lattice line
! mad has a limit of 5000 characters so we may need to break the ring into
! pieces

  line_name(1) = 'line_1'
  line_name(2) = 'line_2'
  line_name(3) = 'line_3'

  il = 0
  init_needed = .true.
  line = ' '

  do n = ie1, ie2

    ele = ring%ele_(n)
    ix = len_trim(line)

! replace element name containing "/" with "_"

    parsing = .true.
    name = ele%name
    do while (parsing)
      i = index (name, '\')         ! '
      if (i /= 0) then
        name(i:i) = '_'
      else
        parsing = .false.
      endif
    enddo

!

    if (ix > 60 .or. n == ie2) then
      if (iout >= 75 .or. n == ie2) then
        i = len_trim(ele%name)
        line(ix+1:) = ', ' // name(:i) // ')'
        write (iu, '(2a)') trim(line)
        line = ' '
        init_needed = .true.
      else
        write (iu, '(2a)') trim(line(:ix)), ', &'
        iout = iout + 1
        line = '   ' // name
      endif

    elseif (init_needed) then
      write (iu, *)
      write (iu, *) '!---------------------------------'
      write (iu, *)
      il = il + 1
      line = line_name(il) // ': line = (' // name
      iout = 0
      init_needed = .false.

    else
      line(ix+1:) = ', ' // name
    endif

  enddo

  write (iu, *)
  write (iu, *) '!---------------------------------'
  write (iu, *)
  write (iu, *) 'ring: line = (line_1', (', ', line_name(i), i = 2, il), ')'

  if (ring%param%lattice_type /= circular_lattice$) then
    write (iu, *)
    write (iu, *) '!---------------------------------'
    write (iu, *)
    write (iu, '(3(a, es12.5))') &
        'TWISS, betx =', ele%x%beta, ', bety =', ele%y%beta, ', &'
    write (iu, '(5x, 3(a, es12.5))') &
        'alfx =', ele%x%alpha, ', alfy =', ele%y%alpha, ', &'
    write (iu, '(5x, 3(a, es12.5))') &
        'dx =', ele%x%eta, ', dpx = ', ele%x%etap, ', &'
    write (iu, '(5x, 3(a, es12.5))') &
        'dy =', ele%y%eta, ', dpy = ', ele%y%etap
  endif

  print *, 'Written: ', trim(mad_file)

!-----------------------------------------------------------------------------
contains

subroutine element_out (ele)

  type (ele_struct), target :: ele
  integer i
  real(rp), pointer :: val(:)

! replace element name containing "/" with "_"

  name = ele%name
  do while (parsing)
    ix = index (name, '\')  ! '
    if (ix /= 0) then
      name(ix:ix) = '_'
    else
      parsing = .false.
    endif
  enddo
                                                        
! do not make duplicate specs

  do i = 1, n_list
    if (ele%name == name_list(i)) return
  enddo

  n_list = n_list + 1
  name_list(n_list) = ele%name
  val => ele%value

! select key

  select case (ele%key)

! drift

  case (drift$)

    write (line, '(a, es13.5)') trim(name) // ': drift, l =', val(l$)
! beambeam

  case (beambeam$)

    line = trim(name) // ': beambeam'
    call value_to_line (line, val(sig_x$), 'sigx', 'es13.5', 'R')
    call value_to_line (line, val(sig_y$), 'sigy', 'es13.5', 'R')
    call value_to_line (line, val(x_offset$), 'xma', 'es13.5', 'R')
    call value_to_line (line, val(y_offset$), 'yma', 'es13.5', 'R')
    call value_to_line (line, val(charge$), 'charge', 'es13.5', 'R')

! elseparator

  case (elseparator$)

    write (line, '(a, es13.5)') trim(name) // ': elseparator, l =', val(l$)
    hk = val(hkick$)
    vk = val(vkick$)

    if (hk /= 0 .or. vk /= 0) then

      ix = len_trim(line) + 1
      field = 1.0e3 * sqrt(hk**2 + vk**2) * val(beam_energy$) / val(l$)
      write (line(ix:), '(a, es13.5)') ', e =', field

      if (ring%param%particle == positron$) then
        tilt = -atan2(hk, vk) + val(tilt$)
      else
        tilt = -atan2(hk, vk) + val(tilt$) + pi
      endif
      ix = len_trim(line) + 1
      write (line(ix:), '(a, es13.5)') ', tilt =', tilt

    endif

! kicker

  case (kicker$)

    write (line, '(a, es13.5)') trim(name) // ': kicker, l =', val(l$)

    call value_to_line (line, val(hkick$), 'hkick', 'es13.5', 'R')
    call value_to_line (line, val(vkick$), 'vkick', 'es13.5', 'R')
    call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! marker

  case (marker$)

    line = trim(name) // ': marker'

! octupole

  case (octupole$)

    write (line, '(a, es13.5)') trim(name) // ': octupole, l =', val(l$)

    call value_to_line (line, val(k3$), 'k3', 'es13.5', 'R')
    call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! quadrupole

  case (quadrupole$)

    write (line, '(a, es13.5)') trim(name) // ': quad, l =', val(l$)
    call value_to_line (line, val(k1$), 'k1', 'es13.5', 'R')
    call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! sbend

  case (sbend$)

    write (line, '(a, es13.5)') trim(name) // ': sbend, l =', val(l$)

    call value_to_line (line, val(angle$), 'angle', 'es13.5', 'R')
    call value_to_line (line, val(e1$), 'e1', 'es13.5', 'R')
    call value_to_line (line, val(e2$), 'e2', 'es13.5', 'R')
    call value_to_line (line, val(k1$), 'k1', 'es13.5', 'R')
    call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! sextupole

  case (sextupole$)

    write (line, '(a, es13.5)') trim(name) // ': sextupole, l =', val(l$)
    call value_to_line (line, val(k2$), 'k2', 'es13.5', 'R')
    call value_to_line (line, val(tilt$), 'tilt', 'es13.5', 'R')

! rfcavity

  case (rfcavity$)

    write (line, '(a, es13.5)') trim(name) // ': rfcavity, l =', val(l$)
    call value_to_line (line, val(voltage$)/1E6, 'volt', 'es13.5', 'R')
    call value_to_line (line, val(phi0$)+val(dphi0$), 'lag', 'es13.5', 'R')
    call value_to_line (line, val(harmon$), 'harmon', 'i8', 'I')

! lcavity

  case (lcavity$)

    write (line, '(a, es13.5)') trim(name) // ': lcavity, l =', val(l$)
    call value_to_line (line, val(gradient$)*val(l$)/1e6, 'deltae', 'f11.4', 'R')
    call value_to_line (line, val(rf_frequency$)/1e6, 'freq', 'es13.5', 'R')
    call value_to_line (line, val(phi0$)+val(dphi0$), 'phi0', 'es13.5', 'R')

! solenoid

  case (solenoid$)

    write (line, '(a, es13.5)') trim(name) // ': solenoid, l =', val(l$)
    call value_to_line (line, val(ks$), 'ks', 'es13.5', 'R')

! wiggler or solquad must be a matrix

  case (wiggler$, sol_quad$)

    write (iu, '(2a)') trim(name), ': matrix, &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(1,1) =', ele%mat6(1,1),  &
          ', rm(1,2) =', ele%mat6(1,2), ', rm(1,3) =', ele%mat6(1,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(1,4) =', ele%mat6(1,4),  &
          ', rm(1,5) =', ele%mat6(1,5), ', rm(1,6) =', ele%mat6(1,6), ', &'

    write (iu, '(5x, 6(a, es13.5))') 'rm(2,1) =', ele%mat6(2,1),  &
          ', rm(2,2) =', ele%mat6(2,2), ', rm(2,3) =', ele%mat6(2,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(2,4) =', ele%mat6(2,4),  &
          ', rm(2,5) =', ele%mat6(2,5), ', rm(2,6) =', ele%mat6(2,6), ', &'

    write (iu, '(5x, 6(a, es13.5))') 'rm(3,1) =', ele%mat6(3,1),  &
          ', rm(3,2) =', ele%mat6(3,2), ', rm(3,3) =', ele%mat6(3,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(3,4) =', ele%mat6(3,4),  &
          ', rm(3,5) =', ele%mat6(3,5), ', rm(3,6) =', ele%mat6(3,6), ', &'

    write (iu, '(5x, 6(a, es13.5))') 'rm(4,1) =', ele%mat6(4,1),  &
          ', rm(4,2) =', ele%mat6(4,2), ', rm(4,3) =', ele%mat6(4,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(4,4) =', ele%mat6(4,4),  &
          ', rm(4,5) =', ele%mat6(4,5), ', rm(4,6) =', ele%mat6(4,6), ', &'

    write (iu, '(5x, 6(a, es13.5))') 'rm(5,1) =', ele%mat6(5,1),  &
          ', rm(5,2) =', ele%mat6(5,2), ', rm(5,3) =', ele%mat6(5,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(5,4) =', ele%mat6(5,4),  &
          ', rm(5,5) =', ele%mat6(5,5), ', rm(5,6) =', ele%mat6(5,6), ', &'

    write (iu, '(5x, 6(a, es13.5))') 'rm(6,1) =', ele%mat6(6,1),  &
          ', rm(6,2) =', ele%mat6(6,2), ', rm(6,3) =', ele%mat6(6,3), ', &'
    write (iu, '(5x, 6(a, es13.5))') 'rm(6,4) =', ele%mat6(6,4),  &
          ', rm(6,5) =', ele%mat6(6,5), ', rm(6,6) =', ele%mat6(6,6)

! multipole

  case (multipole$, ab_multipole$)

    call multipole_ele_to_kt (ele, ring%param%particle, knl, tilt_, .true.)
    write (line, '(a, es13.5)') trim(name) // ': multipole'  
    do i = 0, 9
      write (str, '(a, i1, a)') 'K', i, 'L'
      call value_to_line (line, knl(i), str, 'es13.5', 'R')
      write (str, '(a, i1)') 'T', i
      call value_to_line (line, tilt_(i), str, 'es13.5', 'R')
    enddo

! unknown

  case default

    print *, 'ERROR: UNKNOWN ELEMENT: ', key_name(ele%key), ele%key
    print *, '       CONVERTING TO MARKER'

    line = trim(name) // ': marker'

  end select

! write element spec to mad file

  do
    if (len_trim(line) < 76) exit
    ix = index(line(61:), ' ')
    write (iu, '(2a)') line(:60+ix), '&'
    line = '    ' // line(60+ix:)
  enddo

  write (iu, '(a)') trim(line(:75))

end subroutine

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
