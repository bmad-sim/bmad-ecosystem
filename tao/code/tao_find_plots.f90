!+
! Subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, blank_means_all, only_visible)
!
! Routine to find a plot using the region or plot name.
! A region or plot name is something like: name = "top"
! A graph name is something like: name  = "top.x"
! A curve name is something like: name  = "top.x.c1"
! The wild card "*" can be used so name = "top.*.c1" could 
!   return "top.x.c1" and "top.y.c1".
!
! If blank_means_all = F (the default), something name = "orbit" would not return any graphs or curves.
! If blank_means_all = T, blank graph or curve fields get interpreted as "*". So something line name = "orbit" 
! is interpreted as "orbit.*.*" and all graphs and all curves  of orbit will be returned. 
! If name = "orbit.g.x" then the setting of blank_means_all will be irrelavent.
!
! Use "@Rnnn" to specify a region by index where "nnn" is the index.
! Use "@Tnnn" to specify a region by index where "nnn" is the index.
! Use "T::<name>" to specify a template plot that would be masked if where = "BOTH".
!
! Input:
!   name            -- Character(*): Name of plot or region.
!   where           -- Character(*): Where to look: 'TEMPLATE', 'REGION', 'BOTH' 
!                        For where = 'BOTH', if something is found in a plot region,
!                        then the templates will not be searched
!   print_flag      -- Logical, optional: If present and False then surpress error
!                        messages. Default is True.
!   blank_means_all -- Logical, optional: If present and True then blank graph or curve fields get  interpreted as "*".
!   only_visible    -- logical, optional: Default is True. If True and s%global%plot_on = True then only include
!                        visible plots. If False then plot visible setting is ignored.
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   plot(:)  -- Tao_plot_array_struct, allocatable, optional: Array of plots. If error => size set to 0.
!   graph(:) -- Tao_graph_array_struct, allocatable, optional: Array of graphs. If error => size set to 0.
!   curve(:) -- Tao_curve_array_struct, allocatable, optional: Array of curves. If error => size set to 0.
!-

subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, blank_means_all, only_visible)

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
character(20) where_str
character(*), parameter :: r_name = 'tao_find_plots'

logical, optional :: print_flag, blank_means_all, only_visible
logical err, have_exact_match, only_vis

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

! Parse name argument

where_str = where
err = .false.
only_vis = logic_option(.true., only_visible) .and. s%global%plot_on

if (name == "") then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BLANK "WHERE" LOCATION')
  call end_stuff(4)
  return
endif

ix = index(name, '.')
if (ix == 0) then
  plot_name = name
  graph_name = ''
else
  plot_name = name(1:ix-1)
  graph_name = name(ix+1:)
endif

np = 0
allocate (p(10))

if (plot_name(1:3) == 'T::') then
  where_str = 'TEMPLATE'
  plot_name = plot_name(4:)
endif

! Old name translation

if (len_trim(plot_name) > 1 .and. index('ddispersion', trim(plot_name)) == 1) plot_name = 'etap_dispersion'

! And find plots

if (plot_name(1:1) == '@') then
  select case (plot_name(1:2))
  case ('@R')
    call tao_find_plot_region (err, plot_name, region, print_flag)
    if (err) then
      if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BAD PLOT REGION INDEX: ' // plot_name)
      call end_stuff(4)
      return
    endif
    call point_to_plot (region%plot, p, np)

  case ('@T')
    call tao_set_integer_value (ix, '', plot_name(3:), err, 1, size(s%plot_page%template))
    if (err) then
      if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BAD PLOT TEMPLATE INDEX: ' // plot_name)
      call end_stuff(4)
      return
    endif
    call point_to_plot (s%plot_page%template(ix), p, np)

  case default
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'CONFUSED PLOT INDEX: ' // plot_name)
    call end_stuff(4)
    return
  end select

!

else
  select case (where_str)
  case ('REGION')
    if (only_vis) then
      n_exact = count(s%plot_page%region%name == plot_name .and. s%plot_page%region%visible) + &
                count(s%plot_page%region%plot%name == plot_name .and. s%plot_page%region%visible)
    else
      n_exact = count(s%plot_page%region%name == plot_name) + count(s%plot_page%region%plot%name == plot_name)
    endif
  case ('BOTH')
    if (only_vis) then
      n_exact = count(s%plot_page%region%name == plot_name .and. s%plot_page%region%visible) + &
                count(s%plot_page%region%plot%name == plot_name .and. s%plot_page%region%visible) + &
                count(s%plot_page%template%name == plot_name)
    else
      n_exact = count(s%plot_page%region%name == plot_name) + &
                count(s%plot_page%region%plot%name == plot_name) + count(s%plot_page%template%name == plot_name)
    endif
  case ('TEMPLATE')
    n_exact = count(s%plot_page%template%name == plot_name)
  case default
    if (logic_option(.true., print_flag)) call out_io (s_fatal$, r_name, 'BAD "WHERE" LOCATION: ' // where)
    call end_stuff(4)
    return
  end select

  have_exact_match = (n_exact /= 0)

  ! Match name to region or plot 

  if (where_str == 'REGION' .or. where_str == 'BOTH') then
    do i = 1, size(s%plot_page%region)
      region => s%plot_page%region(i)
      if (only_vis .and. .not. region%visible) cycle
      if (region%plot%name == '') cycle
      if (plot_name /= '*') then 
        if (have_exact_match) then
          if (region%name /= plot_name .and. region%plot%name /= plot_name) cycle
        else
          if (index(region%name, trim(plot_name)) /= 1 .and. &
              index(region%plot%name, trim(plot_name)) /= 1) cycle
        endif
      endif

      call point_to_plot(region%plot, p, np)
    enddo
  endif

  if (where_str == 'TEMPLATE' .or. (where_str == 'BOTH' .and. np == 0)) then
    do i = 1, size(s%plot_page%template)
      if (s%plot_page%template(i)%phantom) cycle
      if (s%plot_page%template(i)%name == '') cycle
      if (plot_name /= '*') then 
        if (have_exact_match) then
          if (s%plot_page%template(i)%name /= plot_name) cycle
        else
          if (index(s%plot_page%template(i)%name, trim(plot_name)) /= 1) cycle
        endif
      endif

      call point_to_plot(s%plot_page%template(i), p, np)
    enddo
  endif

  ! Allocate space

  if (np == 0) then
    select case (where_str)
    case ('REGION')
      if (n_exact == 0) then
        if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
            'PLOT REGION NOT FOUND: ' // plot_name, 'USE THE COMMAND "show plot" TO SEE A LIST OF REGIONS.')
      else
        if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
            'PLOT REGION DOES NOT HAVE AN ASSOCIATED PLOT: ' // plot_name, 'USE THE COMMAND "show plot" TO SEE A LIST OF REGIONS AND ASSOCIATED PLOTS.')
      endif
    case ('TEMPLATE')
      if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
          'PLOT TEMPLATE NOT FOUND: ' // plot_name, 'USE THE COMMAND "show plot -templates" TO SEE A LIST OF PLOTS.')
    case ('BOTH')
      if (n_exact == 0) then
        if (logic_option(.true., print_flag)) then
          if (only_vis) then
            call out_io (s_error$, r_name, &
                '"' // trim(plot_name) // '" IS NOT THE NAME OF A PLOT REGION THAT CONTAINS A PLOT WITH NOR A PLOT TEMPLATE', &
                'USE THE COMMAND "show plot" AND "show plot -templates" TO SEE A LIST OF REGIONS AND TEMPLATES.')
          else
            if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                '"' // trim(plot_name) // '" IS NOT THE NAME OF A PLOT REGION WITH NOR A PLOT TEMPLATE', &
                'USE THE COMMAND "show plot" AND "show plot -templates" TO SEE A LIST OF REGIONS AND TEMPLATES.')
          endif
        endif
      else
        if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
            'PLOT REGION DOES NOT HAVE AN ASSOCIATED PLOT: ' // plot_name, 'USE THE COMMAND "show plot" TO SEE A LIST OF REGIONS AND ASSOCIATED PLOTS.')
      endif
    end select
    call end_stuff(4)
    return
  endif

endif

! Find the number of graphs and allocate.

if (present(plot)) allocate(plot(np))
if (present(plot)) plot = p(1:np)

if (.not. present(graph) .and. .not. present(curve)) then
  return
endif

ix = index(graph_name, '.')
if (ix == 0) then
  curve_name = ''
else
  curve_name = graph_name(ix+1:)
  graph_name = graph_name(1:ix-1)
endif

if (logic_option(.false., blank_means_all) .and. graph_name == '') graph_name = '*'
if (graph_name == '') then
  call end_stuff(2)
  return
endif

ng = 0
do i = 1, np
  if (p(i)%p%name == '') cycle  ! Blank name means plot is invalid.
  do j = 1, size(p(i)%p%graph)
    if (p(i)%p%graph(j)%name == graph_name .or. graph_name == '*') ng = ng + 1
  enddo
enddo

if (ng == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                  'GRAPH NOT FOUND: ' // trim(plot_name) // '.' // graph_name)
  call end_stuff(4)
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

if (logic_option(.false., blank_means_all) .and. curve_name == '') curve_name = '*'
if (curve_name == '') then
  call end_stuff(1)
  return
endif

nc = 0
do j = 1, ng
  if (.not. allocated(g(j)%g%curve)) cycle
  do k = 1, size(g(j)%g%curve)
    if (g(j)%g%curve(k)%name == curve_name .or. curve_name == '*') nc = nc + 1
  enddo
enddo

if (nc == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'CURVE NOT FOUND: ' // name)
  call end_stuff(4)
  return
endif

allocate (c(nc))
if (present(curve)) allocate (curve(nc))

! Now that we have counted and allocated the matches set the curve pointers

nc = 0
do j = 1, ng
  if (.not. allocated(g(j)%g%curve)) cycle
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

subroutine point_to_plot (this_plot, p, np)

type (tao_plot_struct), target :: this_plot
type (tao_plot_array_struct), allocatable :: p(:), p_temp(:)
integer np

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

!------------------------------------------------------------
! contains

subroutine end_stuff(status)
integer status
!
if (status >= 1 .and. present(curve)) then
  if (allocated(curve)) deallocate(curve)
  allocate(curve(0))
endif

if (status >= 2 .and. present(graph)) then
  if (allocated(graph)) deallocate(graph)
  allocate(graph(0))
endif

if (status >= 3 .and. present(plot)) then
  if (allocated(plot)) deallocate(plot)
  allocate(plot(0))
endif

if (status == 4) err = .true.

end subroutine end_stuff

end subroutine tao_find_plots
