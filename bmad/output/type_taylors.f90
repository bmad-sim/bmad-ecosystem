!+
! Subroutine type_taylors (bmad_taylor, max_order, lines, n_lines, file_id, out_style, clean, out_var_suffix)
!
! Subroutine to print or put in a string array a Bmad taylor map.
! If the lines(:) argument is not present, the element information is printed to the terminal.
!
! Input:
!   bmad_taylor(:)  -- taylor_struct: Array of taylors.
!   max_order       -- integer, optional: Maximum order to print.
!   file_id         -- integer, optional: If present, write output to a file with handle file_id.
!   out_style       -- character(*), optional: Determins the string to be used for the output type column.
!                        '' (default) -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz' If size(bmad_taylor) = 6
!                                     -> 'S1', 'Sx', 'Sy', 'Sz' (Spin quaternions) If size(bmad_taylor) = 4
!                                     -> '1', '2', '3', etc. Otherwiase
!                        'NUMBER'     -> '1', '2', '3', etc. 
!                        'BMAD'       -> Bmad lattice file format. 
!                                          If size(bmad_taylor) = 6 --> orbital map assumed.
!                                          If size(bmad_taylor) = 4 --> spin map assumed.
!                        'SCIBMAD'    -> SciBmad format which is Julia/GTPSA compatible.
!   clean           -- logical, optional: If True then do not include terms whose coefficients
!                       are negligible. Default is false.
!   out_var_suffix  -- character(*), optional: If out_style = 'SCIBMAD', out_var_suffix is used as the
!                       suffix of the variable holding the taylor map. Default is "z". 
!                       For example, if "z" is the suffix then:
!                         Descriptor = "d_z", orbital map name = "v_z", ref orbit name = v0_z, and spin map name = "q_z".
!
! Output:
!   lines(:)     -- character(*), allocatable, optional :: Character array to hold the output. 
!                     If not present, the information is printed to the terminal.
!                     For out_style = 'BMAD', Suggested length of lines characters.
!   n_lines      -- integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_taylors (bmad_taylor, max_order, lines, n_lines, file_id, out_style, clean, out_var_suffix)

use taylor_mod, dummy => type_taylors

implicit none

type (taylor_struct), intent(in), target :: bmad_taylor(:)
type (taylor_term_struct), pointer :: tt
type (taylor_struct) tlr, t2

integer, optional, intent(out) :: n_lines
integer, optional :: max_order, file_id
integer i, j, k, n, ie, nl, ix, nt, max_ord

logical, optional :: clean

character(*), optional :: out_style, out_var_suffix
character(*), optional, allocatable :: lines(:)
character(200), allocatable :: li(:)
character(80) suffix
character(40) fmt1, fmt2, fmt, style
character(4) out_str
character(2), parameter :: spin_out(4) = ['S1', 'Sx', 'Sy', 'Sz']
character(2), parameter :: phase_out(6) = ['X ', 'Px', 'Y ', 'Py', 'Z ', 'Pz']
character(2), parameter :: q_out(4) = ['q0', 'q1', 'q2', 'q3']

! If not allocated then not much to do

nt = size(bmad_taylor)

style = ''
if (present(out_style)) style = out_style

if (.not. associated(bmad_taylor(1)%term)) then
  nl = 2
  allocate (li(nl))
  li(1) = '---------------------------------------------------'
  li(2) = 'A Taylor Map Does Not Exist.' 

! SciBmad format

elseif (style == 'SCIBMAD') then
  suffix = ''
  if (present(out_var_suffix)) suffix = out_var_suffix
  if (suffix == '') suffix = 'z'

  allocate(li(20))
  li = ''
  nl = 0

  max_ord = 0
  do i = 1, nt
    t2 = bmad_taylor(i)
    if (logic_option(.false., clean)) call taylor_clean(t2)
    do j = 1, size(t2%term)
      max_ord = max(max_ord, sum(t2%term(j)%expn))
    enddo
  enddo
  if (present(max_order)) max_ord = min(max_ord, max_order)

  if (nt == 6) then
    nl=nl+1; li(nl) = 'using GTPSA'
    nl=nl+1; write (li(nl), '(3a, i0, a, i0, a)') 'd_', trim(suffix), ' = Descriptor(', nt, ', ', max_ord, ')'
    nl=nl+1; write (li(nl), '(3a, 6(es24.16, a))') 'v0_', trim(suffix), ' = [', bmad_taylor(1)%ref, ','
    do j = 2, 6
      nl=nl+1; write (li(nl), '(8x, es24.16, a)') bmad_taylor(j)%ref, ','
    enddo
    n = len_trim(li(nl))
    li(nl) = li(nl)(1:n-1) // ']'
    nl=nl+1; write (li(nl), '(5a)') 'v_', trim(suffix), ' = zeros(TPS64{d_', trim(suffix), '}, 6)'
  else
    nl=nl+1; li(nl) = 'using ReferenceFrameRotations'
    nl=nl+1; write(li(nl), '(3a)') 'q_', trim(suffix), ' = Quaternion{TPS64{d_', trim(suffix), '}}(0,0,0,0)'
  endif

  do i = 1, nt
    t2 = bmad_taylor(i)
    if (logic_option(.false., clean)) call taylor_clean(t2)
    nullify (tlr%term)
    call sort_taylor_terms (t2, tlr)

    do j = 1, size(tlr%term)
      if (nl >= size(li)) call re_allocate(li, 2*nl, .true., '')

      tt => tlr%term(j)

      if (present(max_order)) then
        if (sum(tt%expn) > max_order) cycle
      endif

      if (nt == 4) then
        nl=nl+1; write (li(nl), '(5a, 7(i0, a), es24.16)') 'q_', trim(suffix), '.', q_out(i), '[', i, '][[', (tt%expn(k), ',', k = 1, 5), tt%expn(6), ']] =', tt%coef  
      else
        nl=nl+1; write (li(nl), '(3a, 7(i0, a), es24.16)') 'v_', trim(suffix), '[', i, '][[', (tt%expn(k), ',', k = 1, 5), tt%expn(6), ']] =', tt%coef  
      endif
    enddo
  enddo

  call re_allocate(li, nl)

! Bmad lattice format

elseif (style == 'BMAD') then
  allocate(li(10))
  li = ''
  nl = 1

  do i = 1, nt
    t2 = bmad_taylor(i)

    if (logic_option(.false., clean)) call taylor_clean(t2)
    nullify (tlr%term)
    call sort_taylor_terms (t2, tlr)


    if (nt == 6 .and. taylor_coef(t2, taylor_expn([i])) == 0) then
      li(nl) = trim(li(nl)) // '{' // int_str(i) // ': 1.0 | ' // int_str(i) // '},'
    endif

    do j = 1, size(tlr%term)
      if (len_trim(li(nl)) > min(100, len(li(1))-30)) then
        nl = nl + 1
        if (nl > size(li)) call re_allocate(li, nl+10, .true., '')
      endif

      tt => tlr%term(j)

      if (present(max_order)) then
        if (sum(tt%expn) > max_order) cycle
      endif

      li(nl) = trim(li(nl)) // ' {'

      if (nt == 4) then
        li(nl) = trim(li(nl)) // spin_out(i) // ': ' // real_str(tt%coef, 16)
      else
        li(nl) = trim(li(nl)) // int_str(i) // ': ' // real_str(tt%coef, 16)
      endif

      if (any(tt%expn > 6) .or. sum(tt%expn) > 10) then
        write (li(nl), '(a, 6(a, i0), a)') trim(li(nl)) // ',', (' ', tt%expn(k), k = 1, 6), '},'
      else
        li(nl) = trim(li(nl)) // ' |'
        ix = len_trim(li(nl)) + 1
        do ie = 1, 6
          do k = 1, tt%expn(ie)
            ix = ix + 1
            li(nl)(ix:ix) = int_str(ie)
          enddo
        enddo
        li(nl) = trim(li(nl)) // '},'
      endif

    enddo
  enddo

  call re_allocate(li, nl)
  do i = 1, nl
    li(i) = '    ' // li(i)
  enddo
  n = len_trim(li(nl))
  li(nl) = li(nl)(:n-1)   ! Strip trailing comma

! Normal case

else
  nl = 8 + nt + sum( [(size(bmad_taylor(i)%term), i = 1, nt) ])
  allocate(li(nl))

  write (li(1), '(a)') ' Taylor Terms:'
  if (style == 'SPIN') then
    write (li(2), '(a)') ' Out      Coef               Exponents           Order'
  else
    write (li(2), '(a)') ' Out      Coef               Exponents           Order       Reference'
  endif
  nl = 2

  fmt1 = '(a, f22.14,  1x, 6i3, i9, f18.9)'
  fmt2 = '(a, es22.13, 1x, 6i3, i9, f18.9)'

  do i = 1, nt
    nl=nl+1; li(nl) = ' ---------------------------------------------------'

    if (style == '') then
      select case (nt)
      case (4);     style = 'SPIN'
      case (6);     style = 'PHASE'
      case default; style = 'NUMBER'
      end select
    endif

    select case(style)
    case ('NUMBER')
      write (out_str, '(i3, a)') i, ':'
    case ('PHASE')
      if (i <=6) out_str = ' ' // phase_out(i) // ':'
    case ('SPIN')
      if (i <=4) out_str = ' ' // spin_out(i) // ':'
    end select

    t2 = bmad_taylor(i)
    if (logic_option(.false., clean)) call taylor_clean(t2)

    if (size(t2%term) == 0) then
      nl=nl+1; write (li(nl), '(a, 6x, a)') out_str, 'No Terms. Always evaluates to zero.'

    else
      nullify (tlr%term)
      call sort_taylor_terms (t2, tlr)

      do j = 1, size(t2%term)
        tt => tlr%term(j)

        if (present(max_order)) then
          if (sum(tt%expn) > max_order) cycle
        endif

        if (1d-6 < abs(tt%coef) .and. abs(tt%coef) < 1d5) then
          fmt = fmt1
        else
          fmt = fmt2
        endif

        if (j == 1 .and. style /= 'SPIN') then
          nl=nl+1; write (li(nl), fmt) out_str, tt%coef, (tt%expn(k), k = 1, 6), sum(tt%expn), bmad_taylor(i)%ref
        else
          nl=nl+1; write (li(nl), fmt) out_str, tt%coef, (tt%expn(k), k = 1, 6), sum(tt%expn)
        endif
      enddo

      deallocate (tlr%term)
    endif

  enddo
endif

! Finish

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print '(1x, a)', trim(li(i))
  enddo
endif

if (present(file_id)) then
  do i = 1, nl
    write (file_id, '(a)') trim(li(i))
  enddo
endif

end subroutine type_taylors

