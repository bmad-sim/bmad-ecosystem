module tune_plane_res_mod

  use quick_plot

! Structure for holding the input parameters

  type res_params_struct
    real(rp) x_min, x_max, y_min, y_max
    real(rp) q_s
    real(rp) length
    real(rp) :: x_win_len = 776, y_win_len = 600
    real(rp) scale
    integer pqr_max
    integer r_max
    integer p_max, q_max, pq_max
    integer p_restrict, q_restrict
    integer sum_diff
    integer x_div, y_div
    integer places
    character(16) units
    logical plot_all_1st_and_2nd
    logical show_labels
  end type

! Structure for a single resonance line

  type res_line_struct
    integer p, q, r, n
    integer tag
    logical overlap
    real(rp) x(2)           ! x-coords of end points of resonance line on graph
    real(rp) y(2)           ! y-coords of end points of resonance line on graph
    real(rp) x_lab, y_lab   ! label position
    real(rp) angle_lab      ! label angle in degrees
  end type

! Combined structure

  type res_struct
    type (res_line_struct) line(500)
    integer num_line
  end type

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine res_line_calc (param, res)
!
! Subroutine to calculate the resonance lines for a given region in the tune plane.
!
! Modules needed:
!   use tune_plane_res_mod
!
! Input:
!   param -- Res_params_struct: Input parameters.
!
! Output:
!   res -- Res_struct: List of resonance lines.
!-

subroutine res_line_calc (param, res)

  type (res_params_struct) param
  type (res_struct), target :: res
  type (res_line_struct), pointer :: line

  integer p, q, r, n, i, j, k, ix, kmax, nmin, nmax
  integer nl, ix_pix, iy_pix
  integer plim, rlim, res_tag
  integer pj, qj, rj, nj, ix_x(2), ix_y(2)
  integer fac1

  real(rp) e, res_offset, min_space, strength
  real(rp) n12, n14, n23, n34
  real(rp) x1, y1, x2, y2, xi, yi, xj, yj
  real(rp) fac2, lab_height, x_labj, y_labj
  real(rp) distance2, xi2, yi2, xj2, yj2, fac, orderj, nj_equiv
  real(rp) spacing, old_spacing

! 

  if (param%y_min == 0 .and. param%y_max == 0) then
    param%y_min = param%x_min
    param%y_max = param%x_max
  endif

  if (param%x_max <= param%x_min) then
    print *, 'ERROR param%X_MAX LESS THAN OR EQUAL TO param%X_MIN'
    print *, 'param%X_MAX =', param%x_max
    print *, 'param%X_MIN =', param%x_min
    stop
  endif

  if (param%x_max <= param%x_min) then
    print *, 'ERROR param%X_MAX LESS THAN OR EQUAL TO param%X_MIN'
    print *, 'param%X_MAX =', param%x_max
    print *, 'param%X_MIN =', param%x_min
    stop
  endif

  res%num_line = 0

! write the parameters to the data file

  if (param%units == 'kHz') then
    param%scale = 390.32
    param%places = 0
  elseif (param%units == '1') then
    param%scale = 1.0
    param%places = 2
  else
    print *, 'ERROR: BAD "param%UNITS" PARAMETER'
    stop
  endif

  param%x_min = param%x_min / param%scale
  param%x_max = param%x_max / param%scale
  param%y_min = param%y_min / param%scale
  param%y_max = param%y_max / param%scale
  param%q_s  = param%q_s  / param%scale

  if (max(param%x_min-param%x_max, param%y_min-param%y_max) > 2 .or. param%q_s > 2) then
    print *,  &
      'ERROR: param%X_MIN, param%X_MAX, param%Y_MIN, param%Y_MAX, OR param%Q_S NOT MATCHED TO AXIS_TOGGLE'
    stop
  endif

  res_offset = .002 * max((param%x_max-param%x_min), (param%y_max-param%y_min))
  min_space  = .06 * max((param%x_max-param%x_min), (param%y_max-param%y_min))
  lab_height = .03 * max((param%x_max-param%x_min), (param%y_max-param%y_min))

! search the tune plane for the resonances and
! store the results in array RES%LINE(:)
!
! from the symmetry (p, q, r, n) -> (-p, -q, -r, -n) we can restrict
! ourselves to:
!            0 >= q  and  (p >=0 when q = 0)
!
! resonances are taged by a label depending on which of the sides of
! the plot boundry they cross:
!                 1  =>  cross x = param%x_min
!                 2  =>  cross y = param%y_min
!                 3  =>  cross x = param%x_max
!                 4  =>  cross y = param%y_max
! since every resonance line crosses two sides the two labels are combined
! to form a 2 digit number with the lowest number first.  e.g.
!          23  =>  cross y = param%y_min and x = param%x_max, etc.

  q_loop: do q = 0, param%pqr_max

    plim = param%pqr_max - q

    p_loop: do p = -plim, plim

      rlim = min(param%r_max, param%pqr_max - (abs(p) + abs(q)))

! make sure that one of the restrictions is not violated. if a violation occurs
! go to the end of the do loop

      if (p < 0 .and. q == 0) cycle p_loop  ! do not duplicate resonances

      if (param%p_max >= 0 .and. abs(p) > param%p_max) cycle p_loop
      if (param%q_max >= 0 .and. abs(q) > param%q_max) cycle p_loop
      if (param%pq_max >= 0 .and. abs(p)+abs(q) > param%pq_max) cycle p_loop

      if (.not. param%plot_all_1st_and_2nd .or. abs(p) + abs(q) > 2) then
        if (param%p_restrict == 0 .and. mod(abs(p), 2) == 1) cycle p_loop
        if (param%p_restrict == 1 .and. mod(abs(p), 2) == 0) cycle p_loop
        if (param%q_restrict == 0 .and. mod(abs(q), 2) == 1) cycle p_loop
        if (param%q_restrict == 1 .and. mod(abs(q), 2) == 0) cycle p_loop
        if (param%sum_diff == 0 .and. p < 0 .and. q > 0) cycle p_loop
        if (param%sum_diff == 1 .and. p > 0 .and. q > 0) cycle p_loop
      endif

! now loop over r

      r_loop: do r = -rlim, rlim

! N12 is the (non-integer) value of N for the line going through
! the corner formed by the intersection of the #1 and #2 sides, etc.
! E is added in order to suppress resonance lines that just nick the
! edge of the boundry

        e = .0001

        n12 = p*(param%x_min+e) + q*(param%y_min+e) + r*param%q_s
        n23 = p*(param%x_max-e) + q*(param%y_min+e) + r*param%q_s
        n34 = p*(param%x_max-e) + q*(param%y_max-e) + r*param%q_s
        n14 = p*(param%x_min+e) + q*(param%y_max-e) + r*param%q_s

        nmin = floor(min(n12, n23, n34, n14)) + 1
        nmax = floor(max(n12, n23, n34, n14))

        n_loop: do n = nmin, nmax

          res%num_line = res%num_line + 1
          line => res%line(res%num_line)

          line%angle_lab = atan2(float(-p), float(q)) * 180 / pi
          if (q == 0) line%angle_lab = 90.0   ! instead of -90

          line%p = p
          line%q = q
          line%r = r
          line%n = n

! which sides does the line go through?
! (x1, y1) and (x2, y2) are the points at which the res line goes through the
! boundry

          if (n <= n12 .and. n < n34) then
            line%tag = 23
            line%x(1) = (n-r*param%q_s-q*param%y_min)/p
            line%y(1) = param%y_min
            line%x(2) = param%x_max
            line%y(2) = (n-r*param%q_s-p*param%x_max)/q

          else if (n <= n14 .and. n < n23) then
            line%tag = 12
            line%x(1) = param%x_min
            line%y(1) = (n-r*param%q_s-p*param%x_min)/q
            line%x(2) = (n-r*param%q_s-q*param%y_min)/p
            line%y(2) = param%y_min

          else if (n >= n12 .and. n > n34) then
            line%tag = 14
            line%x(1) = param%x_min
            line%y(1) = (n-r*param%q_s-p*param%x_min)/q
            line%x(2) = (n-r*param%q_s-q*param%y_max)/p
            line%y(2) = param%y_max

          else if (n >= n14 .and. n > n23) then
            line%tag = 34
            line%x(1) = (n-r*param%q_s-q*param%y_max)/p
            line%y(1) = param%y_max
            line%x(2) = param%x_max
            line%y(2) = (n-r*param%q_s-p*param%x_max)/q

          else if (n >= n12 .and. n <= n14) then
            line%tag = 13
            line%x(1) = param%x_min
            line%y(1) = (n-r*param%q_s-p*param%x_min)/q
            line%x(2) = param%x_max
            line%y(2) = (n-r*param%q_s-p*param%x_max)/q

          else
            line%tag = 24
            line%x(1) = (n-r*param%q_s-q*param%y_min)/p
            line%y(1) = param%y_min
            line%x(2) = (n-r*param%q_s-q*param%y_max)/p
            line%y(2) = param%y_max

          endif

! mark all higher order resonances that overlaps a lower order one
! e.g. (2, 2, 0, 4) overlaps (1, 1, 0, 2).

          line%overlap = .false.      ! assume no overlap

          do j = 2, (max(abs(p),abs(q)))
            if (mod(p, j) == 0  .and.  mod(q, j) == 0 .and.  &
                     mod(r, j) == 0  .and.  mod(n, j) == 0) then
              line%overlap = .true.      ! resonance overlap
              exit
            endif
          enddo    ! j

        enddo n_loop

      enddo r_loop

    enddo p_loop

  enddo q_loop  

! the border lines are treated as resonance lines for the label placement 
! calculation

  res%line(res%num_line+1)%p = 1
  res%line(res%num_line+1)%q = 0
  res%line(res%num_line+1)%x(1) = param%x_min

  res%line(res%num_line+2)%p = 0
  res%line(res%num_line+2)%q = 1
  res%line(res%num_line+2)%x(1) = param%y_min

  res%line(res%num_line+3)%p = 1
  res%line(res%num_line+3)%q = 0
  res%line(res%num_line+3)%x(1) = param%x_max

  res%line(res%num_line+4)%p = 0
  res%line(res%num_line+4)%q = 1
  res%line(res%num_line+4)%x(1) = param%y_max

!--------------------------------------------------------------------------
! find label positions
                                               
! we look for a clear place on the tune plane to put the label. The goodness
! of a possible label point is measured by SPACING. 0 <= SPACING <= 2 with
! 0. being the worst and 2 being the best case. SPACING is essentually a measure
! of the distance from the possible label point to the nearest res line. 

  i_loop: do i = 1, res%num_line

    if (res%line(i)%overlap) cycle

    p = res%line(i)%p
    q = res%line(i)%q
    r = res%line(i)%r
    n = res%line(i)%n
    res_tag = res%line(i)%tag

    x1 = res%line(i)%x(1)
    y1 = res%line(i)%y(1)
    x2 = res%line(i)%x(2)
    y2 = res%line(i)%y(2)

! KMAX is the number of possible label points along the res line

    kmax = max(3*sqrt((x2 - x1)**2 + (y2 - y1)**2)/min_space, 1.0)
    old_spacing = -1   ! init

    k_loop: do k = 1, kmax

      xi = x1 + (x2-x1)*mod((2.0*k+1+2*(kmax/4))/(2.0*kmax), 1.0)
      yi = y1 + (y2-y1)*mod((2.0*k+1+2*(kmax/4))/(2.0*kmax), 1.0)

      fac = sqrt(1.0*(p**2+q**2))
      xi2 = xi + (lab_height*p)/fac
      yi2 = yi + (lab_height*q)/fac

      spacing = 2.0

! check the other res lines to see if they are too close to the possible
! label point

      j_loop: do j = 1, res%num_line+4

        pj = res%line(j)%p
        qj = res%line(j)%q
        rj = res%line(j)%r
        nj = res%line(j)%n

        x_labj = res%line(j)%x_lab
        y_labj = res%line(j)%y_lab
        orderj = abs(pj) + abs(qj) + abs(rj)

        if (res%line(j)%overlap) cycle      ! skip overlapping resonances
        if ((p*qj - pj*q) .eq. 0) cycle     ! skip parallel lines

        if (j .gt. res%num_line) then   ! boundry line
          nj_equiv = res%line(j)%x(1)
          orderj = 0.
        else
          nj_equiv = nj - rj*param%q_s
        endif

        xj = 1.0*((n-r*param%q_s)*qj - (nj_equiv)*q)/(p*qj - pj*q)   ! intersection point
        yj = 1.0*((nj_equiv)*p - (n-r*param%q_s)*pj)/(p*qj - pj*q)   ! intersection point

        xj2 = xj + qj*fac*lab_height/(p*qj - pj*q)
        yj2 = yj - pj*fac*lab_height/(p*qj - pj*q)

        distance2 = min((xi-xj)**2+(yi-yj)**2, &
                              (xi2-xj2)**2+(yi2-yj2)**2) / min_space**2

        if (sqrt((x_labj-xj)**2+(y_labj-yj)**2) .lt. min_space .or. &
              sqrt((x_labj-xj2)**2+(y_labj-yj2)**2) .lt. min_space) then 
          distance2 = distance2 / 4.
          fac1 = param%pqr_max
          fac2 = 0.
        else 
          fac1 = (param%pqr_max + 1. - orderj)
          fac2 = (orderj/param%pqr_max) ** 6
        endif

        spacing = min(spacing, distance2**fac1 + fac2)

! quit this point if worse than before

        if (spacing < old_spacing) cycle k_loop

      enddo j_loop

      res%line(i)%x_lab = xi
      res%line(i)%y_lab = yi

      old_spacing = spacing

      if (spacing .gt. 1.9) exit ! can't do much better so next label

    enddo k_loop

  enddo i_loop

! Rescale


  do i = 1, res%num_line
    line => res%line(i)
    line%x = line%x * param%scale
    line%y = line%y * param%scale
    line%x_lab = line%x_lab * param%scale
    line%y_lab = line%y_lab * param%scale
  enddo

  param%x_min = param%x_min * param%scale
  param%x_max = param%x_max * param%scale
  param%y_min = param%y_min * param%scale
  param%y_max = param%y_max * param%scale

end subroutine


!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine res_line_write_list (i_unit, param, res)

  implicit none

  type (res_params_struct) param
  type (res_struct), target :: res
  type (res_line_struct), pointer :: line

  integer i_unit, i

!

  write (i_unit, '(a)')        '&res_params'
  write (i_unit, '(a, f12.6)') '  param%x_min = ', param%x_min
  write (i_unit, '(a, f12.6)') '  param%x_max = ', param%x_max
  write (i_unit, '(a, f12.6)') '  param%y_min = ', param%y_min
  write (i_unit, '(a, f12.6)') '  param%y_max = ', param%y_max
  write (i_unit, '(a, i0)')    '  param%x_div = ', param%x_div
  write (i_unit, '(a, i0)')    '  param%y_div = ', param%y_div
  write (i_unit, '(a, f12.6)') '  param%Q_s   = ', param%Q_s
  write (i_unit, '(a, a)')     '  param%units = ', param%units
  write (i_unit, '(a, i0)')    '  param%pqr_max = ', param%pqr_max
  write (i_unit, '(a, i0)')    '  param%pq_max  = ', param%pq_max
  write (i_unit, '(a, i0)')    '  param%p_max   = ', param%p_max  
  write (i_unit, '(a, i0)')    '  param%q_max   = ', param%q_max  
  write (i_unit, '(a, i0)')    '  param%r_max   = ', param%r_max  
  write (i_unit, '(a, i0)')    '  param%p_restrict = ', param%p_restrict
  write (i_unit, '(a, i0)')    '  param%q_restrict = ', param%q_restrict
  write (i_unit, '(a, i0)')    '  param%sum_diff   = ', param%sum_diff
  write (i_unit, '(a, l1)')    '  param%plot_all_1st_and_2nd = ', param%plot_all_1st_and_2nd
  write (i_unit, '(a, l1)')    '  param%show_labels          = ', param%show_labels
  write (i_unit, '(a)')        '/'
  write (i_unit, *)

  write (i_unit, *) '    P     Q     R     N'
  do i = 1, res%num_line
    line => res%line(i)
    write (i_unit, '(4i6)') line%p, line%q, line%r, line%n
  enddo

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine res_line_plot (plot_type, param, res)

  implicit none

  type (res_params_struct) param
  type (res_struct), target :: res
  type (res_line_struct), pointer :: line
  
  real(rp) x_off, y_off, ang, dQx, dQy
  real(rp) x1, x_page_len, y1, y_page_len, x_graph_len, y_graph_len

  integer i, nl, wid, pqr_max
  integer id, places1

  character(*) plot_type
  character(32) str, fmt
  character(64) lines(50)


! Open the page

  call qp_open_page (plot_type, id, param%x_win_len, param%y_win_len, "POINTS")
  call qp_get_layout_attrib ("PAGE", x1, x_page_len, y1, y_page_len, "POINTS")

! Place the graph. Make sure that the aspect ratio is correct.

  call qp_set_page_border (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, "POINTS")
  call qp_set_box (1, 1, 1, 1)

  x_graph_len = x_page_len - 260
  y_graph_len = y_page_len - 130
  dQx = param%x_max - param%x_min
  dQy = param%y_max - param%y_min
  if (x_graph_len / y_graph_len > dQx / dQy) then
    x_graph_len = y_graph_len * dQx / dQy
  else
    y_graph_len = x_graph_len * dQy / dQx
  endif

  call qp_set_graph_placement (60.0_rp, x_graph_len, 50.0_rp, y_graph_len, "POINTS")

! plot 

  call qp_set_axis ('X', param%x_min, param%x_max, param%x_div, param%places, 'Fractional Horizontal Tune, Q\dx\u')
  call qp_set_axis ('Y', param%y_min, param%y_max, param%y_div, param%places, 'Fractional Vertical Tune, Q\dy\u')
  call qp_set_graph_attrib (draw_grid = .false.)
  call qp_draw_axes 
  call qp_draw_graph_title('Resonance Lines')

  pqr_max = 0
  do i = 1, res%num_line
    line => res%line(i)
    pqr_max = max(pqr_max, abs(line%p) + abs(line%q) + abs(line%r)) 
  enddo

!

  do i = 1, res%num_line

    line => res%line(i)

    if (line%overlap) cycle

    wid = 2 * (pqr_max - abs(line%p) - abs(line%q) - abs(line%r)) + 1

    if (line%r == 0) then
      call qp_draw_line (line%x(1), line%x(2), line%y(1), line%y(2), width = wid, line_pattern = 'solid')
    else
      call qp_draw_line (line%x(1), line%x(2), line%y(1), line%y(2), width = wid, line_pattern = 'dashed')
    endif

    call qp_draw_symbol (line%x_lab, line%y_lab, type = 'circle')

    write (str, '(5(a, i0))') '(', line%p, ',', line%q, ',', line%r, ',', line%n, ')'

    ! x_off, y_off are to keep the text away from the resonance line.

    ang = line%angle_lab
    call qp_from_inch_rel (-0.05 * sin(ang*pi/180), 0.05 * cos(ang*pi/180), x_off, y_off, "DATA")
    call qp_draw_text (str, line%x_lab+x_off, line%y_lab+y_off, &
                                      justify = "CB", height = 10.0_rp, angle = ang)

  enddo

!

  lines(1) = 'Resonance line restrictions:'
  lines(2) = '  [pQ\dx\u + qQ\dy\u + rQ\ds\u = n]'
  write (lines(3), '(a, i0)') '  |p| + |q| + |r| <= ', param%pqr_max
  write (lines(4), '(a, i0)') '  |r| <= ', param%r_max
  nl = 4

  if (param%p_max >= 0) then
    nl = nl + 1
    write (lines(nl), '(a, i0)') '  |p| <= ', param%p_max
  endif

  if (param%q_max >= 0) then
    nl = nl + 1
    write (lines(nl), '(a, i0)') '  |q| <= ', param%q_max
  endif

  if (param%p_restrict == 0) then
    nl = nl + 1
    lines(nl) = '  p even'
  elseif (param%p_restrict == 1) then
    nl = nl + 1
    lines(nl) = '  p odd'
  endif

  if (param%q_restrict == 0) then
    nl = nl + 1
    lines(nl) = '  q even'
  elseif (param%q_restrict == 1) then
    nl = nl + 1
    lines(nl) = '  q odd'
  endif

  if (param%sum_diff == 0) then
    nl = nl + 1
    lines(nl) = '  Sum Resonances Only'
  elseif (param%sum_diff == 1) then
    nl = nl + 1
    lines(nl) = '  Difference Resonances Only'
  endif

  nl = nl + 1
  lines(nl) = ' '

  if (param%r_max /= 0) then
    nl = nl + 1
    places1 = param%places + 1
    write (fmt, '(a, i0, a)') '(a, f0.', places1, ')'
    write (lines(nl), fmt) '  Q\ds\u = ', param%q_s * param%scale
  endif

  call qp_set_text_attrib ("LEGEND", height = 15.0_rp)
  call qp_draw_text_legend (lines(1:nl))

end subroutine

end module
