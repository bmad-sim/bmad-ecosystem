!+
! Subroutine tao_plot_struct_transfer (plot_in, plot_out)
!
! Subroutine to transfer the information from one plot_struct to another.
! This routine handles keeping the pointers separate.
!
! Input:
!   plot_in  -- Tao_plot_struct: Input structure.
!
! Output:
!   plot_out -- Tao_plot_struct: Output struture.
!-

subroutine tao_plot_struct_transfer (plot_in, plot_out)

use tao_interface, dummy => tao_plot_struct_transfer

implicit none

type (tao_plot_struct), target :: plot_in, plot_out
type (tao_curve_struct), pointer :: c_in, c_out
type (tao_plot_region_struct), pointer :: region
integer i, j, n, ix_plot

! Deallocate plot_out allocatable arrays.

if (allocated(plot_out%graph)) deallocate (plot_out%graph)

! Set plot_out = plot_in. 

region => plot_out%r
ix_plot = plot_out%ix_plot
plot_out = plot_in
plot_out%r => region
plot_out%ix_plot = ix_plot

if (allocated(plot_out%graph)) then
  do i = 1, size(plot_out%graph)
    plot_out%graph(i)%p => plot_out
    if (allocated (plot_out%graph(i)%curve)) then
      do j = 1, size(plot_out%graph(i)%curve)
        plot_out%graph(i)%curve(j)%g => plot_out%graph(i)
      enddo
    endif
  enddo
endif

end subroutine
