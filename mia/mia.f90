program mia

  use orbit_mia
  use mia_plot
  use mia_input
  use mia_matrixops

  implicit none

  character* 1 char
  integer :: nset, istat, i_set
  logical :: good_input

  type(data_set), allocatable :: data(:)
  type(cbpm_analysis)data_struct

  !Check for bad input.
!  good_input = .false.

  !MIA gives errors if you want to use a number other than 2,
  !so there's no reason to ask for a number.
!  do while (good_input == .false.)
!     istat=in4get1('Number of files to read?',nset) 
!     write (*, "(a)") " Number of files to read? "        
!     accept "(i)", nset
!    if (nset < 1 .or. nset > 2) then
!        print *, "Input must be 1 or 2."
!     else
!        good_input = .true.
!     endif
!  enddo

  nset = 2

  allocate (data(nset))
  allocate (file(nset))

  do i_set = 1, nset             !Do for first input data 
     call read_bpm_data(data(i_set),i_set,nset)
     if (i_set==1) then
        call initialize_structures(data(1), nset)
     endif

     call svd_fft(data(i_set))
!     call logic_get( 'Y', 'N', 'Plot data? (Y/N) ', plot_more)
!     do while(plot_more)
        call plots(data, 1, i_set)      !Use plot_it
!   call logic_get( 'P','C',' (P)lot more to screen or (C)ontinue?', plot_more)
!     end do
  end do

  do i_set =1, nset
     call bpm_ops(data,i_set, nset)
  end do
  call match_tau_column(nset,data)
  call convert_data_from_pi_matrix(data)
  call calculate_with_known_spacing(data)
!  call logic_get( 'Y', 'N', 'Plot data? (Y/N) ', plot_more)
!  do while(plot_more)
    call plots(data, 2, 2)         !Use plot_it2
!    call logic_get( 'P','C',' (P)lot more to screen or (C)ontinue?', plot_more)
!  end do
  call output (data(1))

end program mia





