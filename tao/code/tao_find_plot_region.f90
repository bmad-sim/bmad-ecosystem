!+
! Subroutine tao_find_plot_region (err, where, region, print_flag)
!
! Routine to find a region using the region name.
!
! Input:
!   where      -- Character(*): Region name.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   region   -- Tao_plot_region_struct, pointer: Region found.
!-

subroutine tao_find_plot_region (err, where, region, print_flag)

use tao_set_mod, dummy => tao_find_plot_region

implicit none

type (tao_plot_region_struct), pointer :: region

integer i, ix

character(*) where
character(40) this_region, graph_name
character(28) :: r_name = 'tao_find_plot_region'

logical, optional :: print_flag
logical err

! Parse where argument. where could be something like "r12.g.c" in which case we 
! need to strip off the graph.curve part.

ix = index(where, '.')
if (ix == 0) then
  this_region = where
else
  this_region = where(1:ix-1)
endif

! Region index?

if (this_region(1:2) == '@R') then
  call tao_set_integer_value (ix, '', this_region(3:), err, 1, size(s%plot_page%region))
  if (err) then
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, 'BAD PLOT REGION INDEX: ' // where)
    return
  endif
  region => s%plot_page%region(ix)
  return
endif

! Match plot name to region

err = .false.

call match_word (this_region, s%plot_page%region%name, ix, .true.)

if (ix < 1) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                    'PLOT REGION NOT FOUND: ' // this_region)
  err = .true.
else
  region => s%plot_page%region(ix)  
endif

end subroutine tao_find_plot_region

