!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine read_outline (outside, inside, ring, type_warning)

  use sr_struct
  use sr_interface

  implicit none

  type (wall_struct) inside, outside
  type (ring_struct) ring

  type (outline_struct) outline_(100), outline, outline_in, z
  type (wall_list_struct) vac_ele(2000)
  type (concat_struct) concat

  integer n, i, j, n1, n2, ixx
  integer n_vac_parts, n_outline, ix_in, ix_out, n_concat_parts
  integer lun, ix, ios

  real(rp) s_ave, f, del_s, factor, s_overlay, s_fudge

  character string*80, tab, name*16, pt_name*16, units*16
  parameter (tab = char(9))

  real(rp) fake$ / -9.999e-31 /

  logical found, in_equals_out, in_overlay, type_warning
  logical was_joint_or_flange, is_joint_or_flange

  type input_struct
    real(rp) s, x
    logical phantom
  end type

  type (input_struct) in(200), out(200)

  namelist / element / z, in, out, units, in_equals_out
  namelist / concat_element / concat

! init

  outside%side = outside$
  inside%side = inside$

! read in list of elements

  lun = lunget()
  open (lun, file = 'vac_db_sr.txt', status = 'old', readonly, iostat = ios)
  if (ios == 0) then
    type *, 'Note: Local Vacuum Element File used: VAC_DB_SR.TXT'
  else
    open (lun, file = '$BMAD_LAYOUT:vac_db_sr.txt',  &
                                                     status = 'old', readonly)
  endif

  i = 0
  do
    read(lun, '(a)',end=100) string
    i = i + 1
    ix = index(string, tab)
    vac_ele(i)%name = string(:ix-1)
    call str_upcase (vac_ele(i)%name, vac_ele(i)%name)
    read (string(ix+1:), *) vac_ele(i)%s, vac_ele(i)%len
  end do

  100 continue
  n_vac_parts = i
  close (unit = lun)

  do i = 2, n_vac_parts
    if (vac_ele(i)%s < vac_ele(i-1)%s) then
      type *, 'ERROR: VACUUM ELEMENT LIST FROM VAC_DB_SR.TXT'
      type *, '       NOT IN CORRECT ORDER:'
      type '(7x, i4, f10.2, 2x, a)', i-1, vac_ele(i-1)%s, vac_ele(i-1)%name
      type '(7x, i4, f10.2, 2x, a)', i, vac_ele(i)%s, vac_ele(i)%name
    endif
  end do

! read in outlines of different parts
! use a local file if there is one

  lun = lunget()
  open (lun, file = 'beam_pipe.outline',  status = 'old', &
                                                     readonly, iostat = ios)
  if (ios == 0) then
    type *, 'Note: Local Outline File used: BEAM_PIPE.OUTLINE'
  else
    open (unit = lun, file = '$BMAD_LAYOUT:beam_pipe.outline',  &
                                                     status = 'old', readonly)
  endif

  ixx = 0
  do

    z%name = ' '
    z%blueprint = ' '
    z%n_out = -1
    z%n_in = -1
    z%ix_in_slide = 0
    z%ix_out_slide = 0
    z%has_alley = .false.
    z%zero_is_center = .false.
    z%s_center = 0
    z%overlay = .false.

    in%x = fake$
    out%x = fake$
    in%phantom = .false.
    out%phantom = .false.
    z%in%name = ' '
    z%out%name = ' '
    units = 'METRIC'
    in_equals_out = .false.

    read (lun, nml = element, end = 200)

    z%out%s = out%s
    z%out%x = out%x
    z%out%phantom = out%phantom

    z%in%s = in%s
    z%in%x = in%x
    z%in%phantom = in%phantom

    ixx = ixx + 1

    name = z%name
    do i = 1, size(z%out)
      if (z%out(i)%x == fake$) then
        z%n_out = i - 1
        exit
      endif
      if (z%out(i)%name == ' ') then
        z%out(i)%name = name
      else
        name = z%out(i)%name
      endif
      z%out(i)%blueprint = z%blueprint
    enddo

    if (z%n_out == -1) then
      type *, 'ERROR IN READ_OUTLINE: OUT ARRAY OVERFLOW FOR: ', z%name
      call err_exit
    endif


    if (in_equals_out) then
      z%in = z%out
      z%n_in = z%n_out
    else
      name = z%name
      do i = 1, size(z%in)
        if (z%in(i)%x == fake$) then
          z%n_in = i - 1
          exit
        endif
        if (z%in(i)%name == ' ') then
          z%in(i)%name = name
        else
          name = z%in(i)%name
        endif
        z%in(i)%blueprint = z%blueprint
      enddo

      if (z%n_in == -1) then
        type *, 'ERROR IN READ_OUTLINE: OUT ARRAY OVERFLOW FOR: ', z%name
        call err_exit
      endif
    endif

    do j = 2, z%n_in
      if (z%in(j)%s < z%in(j-1)%s .and. .not. z%has_alley) then
        type *, 'ERROR INSIDE OUTLINE HAS BACKSTEP IN S: ', z%name
        type '(10x, i10, f10.4)', j-1, z%in(j-1)%s
        type '(10x, i10, f10.4)', j, z%in(j)%s
      endif
    enddo

    do j = 2, z%n_out
      if (z%out(j)%s < z%out(j-1)%s .and. .not. z%has_alley) then
        type *, 'ERROR OUTSIDE OUTLINE HAS BACKSTEP IN S: ', z%name
        type '(10x, i10, f10.4)', j-1, z%out(j-1)%s
        type '(10x, i10, f10.4)', j, z%out(j)%s
      endif
    enddo

    if (z%n_out == 0 .xor. z%n_in == 0) then
      type *, 'ERROR: OUTLINE HAS ONE SIDE BUT NOT THE OTHER:', z%name
      type *, '       N_OUT, N_IN:', z%n_out, z%n_in
      z%n_out = 0; z%n_in = 0
    endif

    call str_upcase (units, units)
    if (units == 'ENGLISH') then
      factor = 0.0254
      z%in%s  = z%in%s * factor
      z%in%x  = z%in%x * factor
      z%out%s = z%out%s * factor
      z%out%x = z%out%x * factor
      z%s_center = z%s_center * factor
    elseif (units /= 'METRIC') then
      type *, 'ERROR IN READ_OUTLINE: UNKNOWN UNITS: ', units
      type *, '      FOR ELEMENT: ', name
      call err_exit
    endif

    outline_(ixx) = z

    call str_upcase (name, z%name)

    if (name == 'S FLANGE' .and. z%n_out == 0) then
      if (type_warning) type *, 'WARNING: STANDARD FLANGE ("S FLANGE") IGNORED'
      cycle
    endif

!    if (name == 'D SLD JNT') then
!      if (type_warning) type *, 'WARNING: FUDGING "? SLD JNT" TO BE "D SLD JNT"'
!      outline_(ixx+1) = outline_(ixx)
!      ixx = ixx + 1
!      outline_(ixx)%name = '? SLD JNT'
!      outline_(ixx)%in%name = '? SLD JNT'
!      outline_(ixx)%out%name = '? SLD JNT'
!    endif

  end do

  200 continue
  n_outline = ixx

! form concatenated elements

  rewind (unit = lun)

  do

    concat%part(:)%direction = 1
    concat%overlay = .false.
    concat%part%name = 'End-of-List'

    read (lun, nml = concat_element, end = 300)

    do j = 1, size(concat%part)

      if (concat%part(j)%name == 'End-of-List') then
        n_concat_parts = j - 1
        exit
      endif

      found = .false.
      do i = 1, n_outline
        if (concat%part(j)%name == outline_(i)%name) then
          found = .true.
          outline_in = outline_(i)
          exit
        endif
      enddo

      if (.not. found) then
        type *, 'ERROR: CANNOT FIND OUTLINE FOR CONCATINATION'
        type *, '       CONCATINATION:  ', concat%name
        type *, '       OUTLINE WANTED: ', concat%part(j)%name
        cycle
      endif

      if (concat%name == outline_in%name) then
        type *, 'ERROR: CONCATINATION NAME CANNOT BE THE SAME AS AN EXISTING'
        type *, '       OUTLINE NAME: ', concat%name
      endif

      if (concat%part(j)%direction == -1) then
        call outline_reverse (outline_in, outline_in)
      elseif (concat%part(j)%direction /= 1) then
        type *, 'ERROR: REVERSE/FORWARD DIRECTION NOT SPECIFIED FOR CONCATINATION'
        type *, '       CONCATINATION: ', concat%name
        type *, '       OUTLINE: ', concat%part(j)%name
        cycle
      endif

      if (j == 1) then
        outline = outline_in
      else
        call outline_concat (outline, outline_in, outline)
      endif

    enddo

    if (j == 1) then
      where (outline%in%name == concat%part(1)%name) &
                                       outline%in%name = concat%name
      where (outline%out%name == concat%part(1)%name) &
                                       outline%out%name = concat%name
    endif


    outline%name = concat%name
    n_outline = n_outline + 1
    outline_(n_outline) = outline
    outline_(n_outline)%overlay = concat%overlay

  enddo

  300 continue
  close(unit = lun)

! match outlines with element names

  outline_loop: do j = 1, n_vac_parts

    do i = 1, n_outline
      call str_upcase(name, outline_(i)%name)
      if (name == vac_ele(j)%name) then
        vac_ele(j)%ix_outline = i
        if (outline_(i)%n_out == 0) cycle outline_loop
        if (j .ne. 1) then
          call check_end (outline_(i)%out, 1, 'OUTSIDE BEGINNING')
          call check_end (outline_(i)%in, 1, 'INSIDE BEGINNING')
        endif
        if (j .ne. n_vac_parts) then
          call check_end (outline_(i)%out, outline_(i)%n_out, 'OUTSIDE ENDING')
          call check_end (outline_(i)%in, outline_(i)%n_in, 'INSIDE ENDING')
        endif
        cycle outline_loop
      endif
    enddo
    if (type_warning) type *, 'WARNING: NO OUTLINE FOR: ', vac_ele(j)%name
  enddo outline_loop

! init wall

  outside%pt(0)%s = 0.0
  outside%pt%name = '????????'
  ix_out = -1

  inside%pt(0)%s = 0.0
  inside%pt%name = '????????'
  ix_in = -1

!--------------------------------
! go through all the vacuum elements and calculate what the wall should be

  in_overlay = .false.
  was_joint_or_flange = .false.

  do j = 1, n_vac_parts

    i = vac_ele(j)%ix_outline
    if (i == 0) cycle        ! cannot do anything if we do not have an outline
    if (outline_(i)%n_out == 0 .or. outline_(i)%n_in == 0) cycle
    outline = outline_(i)

! if in an overlay then skip this part if it is within the overlayed region

    if (in_overlay .and. &
              vac_ele(j)%s + outline%out(1)%s < outside%pt(ix_out)%s) cycle

! if current part is an overlay then delete points in the overlayed region

    if (outline%overlay) then

      s_overlay = vac_ele(j)%s + outline%out(1)%s  ! beginning edge of region

      if (ix_out /= -1) then
        do n = ix_out, 0, -1
          if (outside%pt(n)%s .le. s_overlay) exit
        enddo
        ix_out = n
      endif

      if (ix_in /= -1) then
        do n = ix_in, 0, -1
          if (inside%pt(n)%s .le. s_overlay) exit
        enddo
        ix_in = n
      endif

      s_overlay = vac_ele(j)%s + outline%out(outline%n_out)%s ! end edge
      in_overlay = .true.

    else
      in_overlay = .false.
    endif

! stretch if a sliding joint

    if (vac_ele(j)%name(3:) == 'SLD JNT') then
      n1 = outline%ix_out_slide
      n2 = outline%n_out
      if (vac_ele(j)%len == 0) then
        if (type_warning) &
              type *, 'WARNING: NO LENGTH FOR SLIDING JOINT AT:', vac_ele(j)%s
        del_s = 6.00*0.0254 - (outline%out(n2)%s - outline%out(1)%s)
      else
        del_s = (4.37*0.0254 + vac_ele(j)%len) - &
                                  (outline%out(n2)%s - outline%out(1)%s)
      endif
      outline%out(n1:n2)%s = outline%out(n1:n2)%s + del_s
      n1 = outline%ix_in_slide
      n2 = outline%n_in
      outline%in(n1:n2)%s = outline%in(n1:n2)%s + del_s
    endif

    if (outline%s_center /= 0) then
      outline%out%s = outline%out%s - outline%s_center
      outline%in%s  = outline%in%s  - outline%s_center
    elseif (.not. outline%zero_is_center) then
      s_ave = (outline%out(outline%n_out)%s + outline%out(1)%s) / 2
      outline%out%s = outline%out%s - s_ave
      outline%in%s  = outline%in%s  - s_ave
    endif

! fudge if overlap between sliding joint and Standard flange

    is_joint_or_flange = .false.
    if (vac_ele(j)%name(3:) == 'SLD JNT' .or. &
                    vac_ele(j)%name == 'S FLANGE') is_joint_or_flange = .true.

    if (was_joint_or_flange .and. is_joint_or_flange) then
      s_fudge = outside%pt(ix_out)%s - (vac_ele(j)%s + outline%out(1)%s)
      if (s_fudge < 0) s_fudge = 0
      if (type_warning .and. s_fudge .gt. 0.002) then
        type *, 'WARNING: S_FUDGE BETWEEN JOINT AND FLANGE > 2 mm: ', s_fudge
        type *, '         NEAR S =',  outside%pt(ix_out)%s
      endif
    else
      s_fudge = 0
    endif

    was_joint_or_flange = is_joint_or_flange

! calculate outside wall points

    if (outline%n_out > 0) then
      do n = 1, outline%n_out
        ix_out = ix_out + 1
        outside%pt(ix_out)%s = vac_ele(j)%s + outline%out(n)%s + s_fudge
        outside%pt(ix_out)%x = outline%out(n)%x
        outside%pt(ix_out)%name    = outline%out(n)%name
        outside%pt(ix_out)%phantom = outline%out(n)%phantom
        if (n == 1) outside%pt(ix_out)%name = 'ARC'
        if (outline%has_alley) then
          outside%pt(ix_out)%type = possible_alley$
        else
          outside%pt(ix_out)%type = no_alley$
        endif
      end do
    endif

    in_overlay = outline%overlay

! calculate inside wall points

    if (outline%n_in > 0) then
      do n = 1, outline%n_in
        ix_in = ix_in + 1
        inside%pt(ix_in)%s = vac_ele(j)%s + outline%in(n)%s + s_fudge
        inside%pt(ix_in)%x = -outline%in(n)%x
        inside%pt(ix_in)%name    = outline%in(n)%name
        inside%pt(ix_in)%phantom = outline%in(n)%phantom
        if (n == 1) inside%pt(ix_in)%name = 'ARC'
        if (outline%has_alley) then
          inside%pt(ix_in)%type = possible_alley$
        else
          inside%pt(ix_in)%type = no_alley$
        endif
      end do
    endif

  end do

! cleanup

  outside%n_pt_tot = ix_out
  inside%n_pt_tot = ix_in

  call delete_overlaping_wall_points (outside)
  call delete_overlaping_wall_points (inside)

  forall (i = 1:ix_out) outside%pt(i)%ix_pt = i
  forall (i = 1:ix_in)  inside%pt(i)%ix_pt = i

  call create_alley (inside)
  call create_alley (outside)

! check that endpoints are correct

  if (abs(outside%pt(outside%n_pt_tot)%s - ring%ele_(ring%n_ele_ring)%s) > &
                                                                    0.01) then
    type *, 'ERROR: OUTSIDE WALL ENDS AT:', outside%pt(outside%n_pt_tot)%s
    type *, '       AND NOT AT RING END OF:', ring%ele_(ring%n_ele_ring)%s
    call err_exit
  else
    outside%pt(outside%n_pt_tot)%s = ring%ele_(ring%n_ele_ring)%s
  endif

  if (abs(inside%pt(inside%n_pt_tot)%s - ring%ele_(ring%n_ele_ring)%s) > &
                                                                    0.01) then
    type *, 'ERROR: INSIDE WALL ENDS AT:', inside%pt(inside%n_pt_tot)%s
    type *, '       AND NOT AT RING END OF:', ring%ele_(ring%n_ele_ring)%s
    call err_exit
  else
    inside%pt(inside%n_pt_tot)%s = ring%ele_(ring%n_ele_ring)%s
  endif

! write to file

  f = 1/0.0254
  open (unit = 11, file = 'outside_wall.dat', carriagecontrol = 'list')
  write (11, *) &
   '   Ix        S(m)       X(m)    Name                   S(in)       X(in) '
  do n=0,outside%n_pt_tot
    call convert_blanks_to_underscore (outside%pt(n)%name, pt_name)
    write(11, '(1x, i4, f12.3, f12.5, 3x, a16, 2f12.3)') n, &
            outside%pt(n)%s, outside%pt(n)%x, pt_name, &
            f*outside%pt(n)%s, f*outside%pt(n)%x
  end do
  close (unit = 11)

!

  open (unit = 11, file = 'inside_wall.dat', carriagecontrol = 'list')
  write (11, *) &
   '   Ix        S(m)       X(m)    Name                   S(in)       X(in) '
  do n=0,inside%n_pt_tot
    call convert_blanks_to_underscore (inside%pt(n)%name, pt_name)
    write(11, '(1x, i4, f12.3, f12.5, 3x, a16, 2f12.3)') n, &
            inside%pt(n)%s, inside%pt(n)%x, pt_name, &
            f*inside%pt(n)%s, f*inside%pt(n)%x
  end do
  close (unit = 11)

! do some checking

  call check_wall (inside, ring)
  call check_wall (outside, ring)

! segment wall

  call break_wall_into_segments (inside, .1_rp)
  call break_wall_into_segments (outside, .1_rp)

end subroutine
