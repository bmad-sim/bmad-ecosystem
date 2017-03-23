!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, always_allocate)
!
! Routine to find a plot using the region or plot name.
! A region or plot name is something like: name = "top"
! A graph name is something like: name  = "top.x"
! A curve name is something like: name  = "top.x.c1"
! The wild card "*" can be used so name = "top.*.c1" could 
!   return "top.x.c1" and "top.y.c1".
!
! Unless always_allocate = T, The graph(:) array will only be allocated if 
! the graph portion of name is not blank.
!   For example: name = "top" will leave graph(:) unallocated.
! Unless always_allocate = T, The curve(:) array will only be allocated if 
! the curve portion of name is not blank.
!   For example: name = "top.x" will leave curve(:) unallocated.
!
! Input:
!   name       -- Character(*): Name of plot or region.
!   where      -- Character(*): Where to look: 'TEMPLATE', 'REGION', 'BOTH', 'COMPLETE'
!                   For where = 'BOTH', if something is found in a plot region,
!                   then the templates will not be searched
!                   where = 'COMPLETE' is used by the python command and does not allow abbreviations.
!                   'COMPLETE' should not otherwise be used.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!   always_allocate 
!              -- Logical, optional: If present and True then always allocate 
!                   graph(:) and curve(:) arrays except if there is an error. 
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   plot(:)  -- Tao_plot_array_struct, allocatable, optional: Array of plots.
!   graph(:) -- Tao_graph_array_struct, allocatable, optional: Array of graphs.
!   curve(:) -- Tao_curve_array_struct, allocatable, optional: Array of curves.
!-

subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, always_allocate)

use tao_set_mod, dummy => tao_find_plots

implicit none

type (tao_plot_array_struct), allocatable, optional :: plot(:)
type (tao_graph_array_struct), allocatable, optional :: graph(:)
type (tao_curve_array_struct), allocatable, optional :: curve(:)
type (tao_plot_array_struct), allocatable :: p(:)
type (tao_graph_array_struct), allocatable :: g(:)
type (tao_curve_array_struct), allocatable :: c(:)
type (tao_plot_region_struct), pointer :: region
type (tao_plot_array_struct), allocatable :: p_temp(:)

integer i, j, k, ix, np, ng, nc, n_exact

character(*) name, where
character(40) plot_name, graph_name, curve_name
character(28) :: r_name = 'tao_find_plots'

logical, optional :: print_flag, always_allocate
logical err, have_exact_match

! Init

if (present(plot)) then
  if (allocated(plot)) deallocate(plot)
endif

if (present(graph)) then
  if (allocated(graph)) deallocate(graph)
endif

if (present(curve)) then
  if (allocated(curve)) deallocate(curve)
endif

if (allocated(p)) deallocate(p)
if (allocated(g)) deallocate(g)
if (allocated(c)) deallocate(c)

! Parse name argument

err = .false.

if (name == "") then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BLANK "WHERE" LOCATION')
  err = .true.
  return
endif

ix = index(name, '.')
if (ix == 0) then
  plot_name = name
  graph_name = ' '
else
  plot_name = name(1:ix-1)
  graph_name = name(ix+1:)
endif

np = 0
allocate (p(10))

!

if (plot_name(1:1) == '@') then
  if ((plot_name(1:2) == '@R' .and. where == 'TEMPLATE') .or. &
      (plot_name(1:2) == '@T' .and. where == 'REGION')) then
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'CONFUSED PLOT INDEX: ' // plot_name)
    err = .true.
    return
  endif

  select case (plot_name(1:2))
  case ('@R')
    call tao_find_plot_region (err, plot_name, region, print_flag)
    if (err) return
    call point_to_plot (region%plot)

  case ('@T')
    call tao_set_integer_value (ix, '', plot_name(3:), err, 1, size(s%plot_page%template))
    if (err) then
      if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BAD PLOT TEMPLATE INDEX: ' // plot_name)
      return
    endif
    call point_to_plot (s%plot_page%template(ix))

  case default
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'CONFUSED PLOT INDEX: ' // plot_name)
    err = .true.
    return
  end select

!

else
  select case (where)
  case ('REGION')
    n_exact = count(s%plot_page%region%name == plot_name) + &
              count(s%plot_page%region%plot%name == plot_name)
  case ('BOTH')
    n_exact = count(s%plot_page%region%name == plot_name) + &
              count(s%plot_page%region%plot%name == plot_name) + &
              count(s%plot_page%template%name == plot_name)
  case ('TEMPLATE')
    n_exact = count(s%plot_page%template%name == plot_name)
  case ('COMPLETE')
    n_exact = count(s%plot_page%region%name == plot_name) + &
              count(s%plot_page%template%name == plot_name)
  case default
    if (logic_option(.true., print_flag)) call out_io (s_fatal$, r_name, 'BAD "WHERE" LOCATION: ' // where)
    call err_exit
  end select

  have_exact_match = (n_exact /= 0)

  ! Match name to region or plot 

  if (where == 'REGION' .or. where == 'BOTH') then
    do i = 1, size(s%plot_page%region)
      if (plot_name /= '*') then 
        if (have_exact_match) then
          if (s%plot_page%region(i)%name /= plot_name .and. s%plot_page%region(i)%plot%name /= plot_name) cycle
        else
          if (index(s%plot_page%region(i)%name, trim(plot_name)) /= 1 .and. &
              index(s%plot_page%region(i)%plot%name, trim(plot_name)) /= 1) cycle
        endif
      endif

      call point_to_plot(s%plot_page%region(i)%plot)
    enddo
  endif

  if (where == 'TEMPLATE' .or. (where == 'BOTH' .and. np == 0)) then
    do i = 1, size(s%plot_page%template)
      if (s%plot_page%template(i)%phantom) cycle
      if (plot_name /= '*') then 
        if (have_exact_match) then
          if (s%plot_page%template(i)%name /= plot_name) cycle
        else
          if (index(s%plot_page%template(i)%name, trim(plot_name)) /= 1) cycle
        endif
      endif

      call point_to_plot(s%plot_page%template(i))
    enddo
  endif

  ! Complete

  if (where == 'COMPLETE') then
    do i = 1, size(s%plot_page%region)
      if (s%plot_page%region(i)%name /= plot_name) cycle
      call point_to_plot(s%plot_page%region(i)%plot)
    enddo

    do i = 1, size(s%plot_page%template)
      if (s%plot_page%template(i)%name /= plot_name) cycle
      call point_to_plot(s%plot_page%template(i))
    enddo
  endif

  ! Allocate space

  if (np == 0) then
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'PLOT NOT FOUND: ' // plot_name)
    err = .true.
    return
  endif

endif

! Find the number of graphs and allocate.

if (present(plot)) allocate(plot(np))
if (present(plot)) plot = p(1:np)

if (.not. present(graph) .and. .not. present(curve)) return

ix = index(graph_name, '.')
if (ix == 0) then
  curve_name = ' '
else
  curve_name = graph_name(ix+1:)
  graph_name = graph_name(1:ix-1)
endif

if (logic_option(.false., always_allocate) .and. graph_name == ' ') graph_name = '*'
if (graph_name == ' ') return

ng = 0
do i = 1, np
  do j = 1, size(p(i)%p%graph)
    if (p(i)%p%graph(j)%name == graph_name .or. graph_name == '*') ng = ng + 1
  enddo
enddo

if (ng == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                  'GRAPH NOT FOUND: ' // trim(plot_name) // '.' // graph_name)
  err = .true.
  return
endif

allocate (g(ng))
if (present(graph)) allocate (graph(ng))

! Now that we have counted and allocated the matches set the graph pointers

ng = 0
do i = 1, np
  do j = 1, size(p(i)%p%graph)
    if (p(i)%p%graph(j)%name == graph_name .or. graph_name == '*') then
      ng = ng + 1
      g(ng)%g => p(i)%p%graph(j)
    endif
  enddo
enddo

if (present(graph)) graph = g

! Find number of curves that match

if (.not. present(curve)) return

if (logic_option(.false., always_allocate) .and. curve_name == ' ') curve_name = '*'
if (curve_name == ' ') return

nc = 0
do j = 1, ng
  do k = 1, size(g(j)%g%curve)
    if (g(j)%g%curve(k)%name == curve_name .or. curve_name == '*') nc = nc + 1
  enddo
enddo

if (nc == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                  'CURVE NOT FOUND: ' // name)
  err = .true.
  return
endif

allocate (c(nc))
if (present(curve)) allocate (curve(nc))

! Now that we have counted and allocated the matches set the curve pointers

nc = 0
do j = 1, ng
  do k = 1, size(g(j)%g%curve)
    if (g(j)%g%curve(k)%name == curve_name .or. curve_name == '*') then
      nc = nc + 1
      c(nc)%c => g(j)%g%curve(k)
    endif
  enddo
enddo

if (present(curve)) curve = c

!------------------------------------------------------------
contains

subroutine point_to_plot (this_plot)

type (tao_plot_struct), target :: this_plot
type (tao_plot_array_struct), allocatable :: p_temp(:)

!

if (np+1 > size(p)) then
  call move_alloc (p, p_temp)
  allocate (p(2*np))
  p(1:np) = p_temp(1:np)
  deallocate (p_temp)
endif

np = np + 1
p(np)%p => this_plot

end subroutine point_to_plot

end subroutine tao_find_plots
