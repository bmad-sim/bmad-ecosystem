!+
! Subroutine tao_plot_struct_transfer (plot_in, plot_out, preserve_region)
!
! Subroutine to transfer the information from one plot_struct to another.
! This routine handles keeping the pointers separate.
!
! Input:
!   plot_in  -- Tao_plot_struct: Input structure.
!   preserve_region -- Logical: If True then plot_out%region is not overwritten.
!
! Output:
!   plot_out -- Tao_plot_struct: Output struture.
!-

subroutine tao_plot_struct_transfer (plot_in, plot_out, preserve_region)

use tao_mod

implicit none

type (tao_plot_struct) plot_in, plot_out
type (tao_plot_region_struct) region

integer i, j, n
logical preserve_region

! Deallocate plot_out pointers

if (associated(plot_out%graph)) then

  do i = 1, size(plot_out%graph)
    if (.not. associated(plot_out%graph(i)%curve)) cycle
    do j = 1, size(plot_out%graph(i)%curve)
      if (.not. associated(plot_out%graph(i)%curve(j)%x_dat)) cycle
      deallocate (plot_out%graph(i)%curve(j)%x_dat)
      deallocate (plot_out%graph(i)%curve(j)%y_dat)
    enddo
    deallocate (plot_out%graph(i)%curve)
  enddo

  deallocate (plot_out%graph)

endif

! Set plot_out = plot_in. 
! This makes the pointers of plot_out point to the same memory as plot_in

region = plot_out%region
plot_out = plot_in
if (preserve_region) plot_out%region = region

! Now allocate storage for the plot_out pointers and 
! copy the data from plot_in to plot_out

if (associated(plot_in%graph)) then

  allocate (plot_out%graph(size(plot_in%graph)))

  do i = 1, size(plot_out%graph)

    plot_out%graph(i) = plot_in%graph(i)

    if (.not. associated (plot_in%graph(i)%curve)) cycle
    allocate (plot_out%graph(i)%curve(size(plot_in%graph(i)%curve)))

    do j = 1, size(plot_out%graph(i)%curve)
      plot_out%graph(i)%curve(j) = plot_in%graph(i)%curve(j)
      if (.not. associated (plot_in%graph(i)%curve(j)%x_dat)) cycle
      n = size(plot_in%graph(i)%curve(j)%x_dat)
      allocate (plot_out%graph(i)%curve(j)%x_dat(n))
      allocate (plot_out%graph(i)%curve(j)%y_dat(n))
      plot_out%graph(i)%curve(j)%x_dat = plot_in%graph(i)%curve(j)%x_dat 
      plot_out%graph(i)%curve(j)%y_dat = plot_in%graph(i)%curve(j)%y_dat 
    enddo

  enddo

endif

end subroutine
