!+
! Subroutine READ_BUTNS_FILE (butns_num, butns, db, read_ok, type_err)
!
! Subroutine to read in the information in a BUTNS.nnnnn file.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   butns_num -- Integer: Number in BUTNS.nnnnn file name.
!   type_err  -- Logical: If True then open error message will be printed.
!                         If False then open error message will be supressed.
!
! Output:
!   butns -- Butns_struct: Orbit information.
!       %lattice    -- Character*40: Lattice name.
!       %save_set   -- Integer: Save set number.
!       %date       -- Character*20: Date orbit was taken
!       %file_num   -- Integer: Equal to butns_num.
!       %turn       -- Integer: Turn number for injection orbits. 0 otherwise.
!       %comment(5) -- Character*72: Comment.
!       %det(0:99)%amp(4) -- Integer: raw button numbers.
!       %det(0:99)%x_orb  -- Real: Horizontal orbit in meters.
!       %det(0:99)%y_orb  -- Real: Horizontal orbit in meters.
!   db    -- Db_struct: Structure holding the steering settings.
!     %csr_horz_cur(i)%cu_now -- CU settings for CSR HORZ CUR
!     %csr_hbnd_cur(i)%cu_now -- CU settings for CSR HBND CUR
!     %csr_vert_cur(i)%cu_now -- CU settings for CSR VERT CUR
!     %csr_hsp_volt(i)%cu_now -- CU settings for CSR HSP VVAL
!     %csr_vsp_volt(i)%cu_now -- CU settings for CSR VSP VOLT
!     %scir_vertcur(i)%cu_now -- CU settings for SCIR VERTCUR
!     %scir_pos_stp(i)%cu_now -- CU settings for SCIR POS STP
!     %scir_enc_cnt(i)%cu_now -- CU settings for SCIR ENC CNT
!     %scir_pos_rd(i)%cu_now  -- CU settings for SCIR POS RD
!   read_ok -- Logical: Set True if butns file was successfuly parsed.
!
! Note: orbit numbers are from 0 to 99
!
! Note: db%csr_hsp_volt is actually obtained from the node CSR HSP VVAL which
! records the actual voltage as opposed to the command. Since CSR HSP VVAL
! is a readback node the values put in db%csr_hsp_volt will not be exactly the
! actual commands (and if the separators have been turned off they will not
! even be approximately the same).
!-

subroutine read_butns_file (butns_num, butns, db, read_ok, type_err)

  use bmad_struct
  use cesr_mod

  implicit none

  type (db_struct) db
  type (butns_struct) butns

  integer vec(120), det_type(120)
  integer i, ix, j, butns_num, iu, lunget, ios, raw(4, 120)
  integer n_node, n_ele

  real(rp) x_orbit(120), y_orbit(120), rdummy

  character line_in*130, file_in*40

  logical read_ok, type_err, err_flag, is_ok(120)

! init comments

  butns%comment = ' '

! compute filename

  read_ok = .false.

  call form_file_name_with_number ('ORBIT', butns_num, file_in, err_flag)
  if (err_flag) return

  butns%file_num = butns_num

! read header line in the data file (to get lattice name)
! and sterring strengths

  iu = lunget()
  open(unit = iu, file = file_in, status = 'old', action = 'READ', iostat = ios)
  if (ios /= 0) then
    if (type_err) print *, &
          'ERROR IN READ_BUTNS_FILE: ERROR OPENING: ', trim(file_in)
    return
  endif

  read (iu, '(a)') line_in
  butns%lattice = line_in(61:)
  butns%date = line_in(30:)                     ! get date
  read (line_in(54:), *) butns%save_set

  do
    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERINGS IN ORBIT FILE'
      goto 1000
    endif
    ix = index(line_in, 'END BUTNS')
    if (ix /= 0) exit
    ix = index(line_in, 'TURN=')
    if (ix /= 0) then
      read (line_in(ix+5:), '(i)', iostat = ios) butns%turn
      if (ios /= 0) then
        print *, 'ERROR IN READ_BUTNS_FILE: ERROR READING TURN NUMBER.'
        print *, '      IN FILE: ', trim(file_in)
      endif
    endif
  enddo

  read (line_in(ix+10:), *, iostat = ios) n_node
  if (ios /= 0) then
    print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE NUMBER IN ORBIT FILE'
    goto 1000
  endif

! read data base valuse stored in the orbit file

  do i = 1, n_node

    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE IN ORBIT FILE'
      exit
    endif

    read (line_in(14:), *, iostat = ios) n_ele

    if (line_in(2:13) == 'CSR COMMENTS') then
      do j = 1, n_ele
        read (iu, '(1x, a)', iostat = ios) butns%comment(j)
        if (ios /= 0) then
          print *, 'ERROR IN ORBIT_READ: ERROR READING COMMENT #', j
          exit
        endif
      enddo
      cycle
    endif

    read (iu, '(12x, 10i6)') vec(1:n_ele)

    if (line_in(2:13) == 'CSR HORZ CUR') then
      db%csr_horz_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR VERT CUR') then
      db%csr_vert_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR HBND CUR') then
      db%csr_hbnd_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR HSP VVAL') then
      call hsp_vval_to_volt (vec, vec)
      n_ele = size(db%csr_hsp_volt)
      db%csr_hsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR VSP VOLT') then
      db%csr_vsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR VERTCUR') then
      db%scir_vertcur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR POS STP') then
      db%scir_pos_stp(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR ENC CNT') then
      db%scir_enc_cnt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR POS RD ') then
      db%scir_pos_rd(1:n_ele)%cu_now = vec(1:n_ele)
    else
      print *, 'ERROR IN ORBIT_READ: UNKNOWN NODE IN ORBIT FILE: ', line_in(2:13)
      goto 1000
    endif

  enddo

! close

  1000 continue
  close (unit = iu)

! read in the raw data

  call butfilget (raw, file_in, rdummy, det_type)      ! read in raw data
  call butcon (raw, 1, 100, y_orbit, x_orbit)
  is_ok = .false.
  call det_ok (raw, 1, 100, det_type, is_ok)

  do i = 1, 100
    j = i
    if (i == 100) j = 0
    butns%det(j)%amp = raw(1:4, i)
    butns%det(j)%x_orb = x_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%y_orb = y_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%ok = is_ok(i)
  end do

  read_ok = .true.

end subroutine

