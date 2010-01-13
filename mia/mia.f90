program mia

  use orbit_mia
  use mia_plot
  use mia_input
  use mia_matrixops
  use mia_parse
  use mia_veto

  implicit none

  character* 1 char
  integer :: istat, i_set, narg, b
  logical :: good_input, more_files, first_run

  type(data_set), allocatable :: data(:)
  integer :: time(8)
  real :: startTime, endTime

  call date_and_time(values=time)
  startTime = time(5) * 3600 + time(6) * 60 &
       + time(7) + 0.001 * time(8)

  nset = 2
  more_files = .false.
  first_run = .true.
  do while (more_files .or. first_run)
     allocate (data(nset))

     do i_set = 1, nset             !Do for first input data 
        !See if there are any command line arguments:
        if (i_set==1) then
           narg = IARGC()
           if (narg>0) then
              call get_args(data)
           endif
        endif
        
        call read_bpm_data(data(i_set),i_set,nset)
        
        if (i_set==1) then
           call initialize_structures(data(1), nset)
        endif
       
        call svd_fft(data(i_set))
        !Plot after analysis; for debugging, single file analysis, or
        !SVD/FFT data only uncomment this.
!        call plots(data, 1, i_set)      !Use plot_it

     end do


871  do i_set =1, nset
        call bpm_ops(data,i_set, nset, first_run)
     end do
     first_run = .false.
!     call deall_file
     call match_tau_column(nset,data)
     call convert_data_from_pi_matrix(data)
     call calculate_with_known_spacing(data)

     if (.not. silentMode) then
        call plots(data, 2)         !Use plot_it2
     end if
     if (vetoBPM) then
        !***Call subroutine to remove a bpm,
        !***redo analysis with new poshis

        call veto(data, inputPasser)
        call clean(.true.)
        call initialize_structures(data(1), nset)

        call svd_fft(data(1))
        call svd_fft(data(2))
        goto 871

     end if

     call output (data)

!Disabled to make batch analysis easier; can uncomment this
!to give the user the option of running MIA with a new set of files
!without exiting
!     if (.not. silentMode) then
!        call logic_get('Y', 'N', 'Repeat calculations with different files?',&
!             more_files)
!     else
!        more_files = .false.
!     endif
     if (.not. more_files) then
        deallocate (data)
     end if
     call clean (more_files)
  enddo
  call date_and_time(values=time)
  endTime = time(5) * 3600 + time(6) * 60 &
       + time(7) + 0.001 * time(8) - startTime
  Print *, "MIA ran in ", endTime, "seconds."

end program mia





