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

  !nset is number of files to analyze--this code wants one file with
  !horizontal shaking and one with vertical shaking.
  nset = 2
  more_files = .false.
  first_run = .true.




  !Whole program can be looped asking for more files.
  !It's disabled for testing but code using variable more_files 
  !can be uncommented.
  do while (more_files .or. first_run)
     allocate (data(nset))

     !Check for and parse command line arguments:
     narg = IARGC()
     if (narg>0) then
        call get_args(data)
     endif


     !Read data and perform SVD and FFT:
     do i_set = 1, nset            
        
        call read_data(data(i_set),i_set,nset)



        !Initialize data structures after reading number of
        !turns, BPMs, etc from first file (second file is expected
        !to match)
        if (i_set == 1) then
           call initialize_struct(data(i_set), .false.)
       end if
       !Parse detector names, get S positions, etc:
       call bpm_ops(data,i_set, nset, first_run)

!        Print *, "Num turns: ", data(i_set)%numturns
!        Print *, "Num dets: ", data(i_set)%bpmproc
       
        call svd_fft(data(i_set))

!**For single-file SVD analysis, uncomment the following line:
!        call plots(data, 1)

     end do


     first_run = .false.
     call match_tau_column(nset,data)


871  continue
     call autoCut(data)
!     call autoVeto(data)
567  continue
     if (vetoBPM) then
        
        call clean(.false.)
        call initialize_struct(data(1), .true.)
        call findPair(nset,data)

        call svd_fft(data(1))
        call svd_fft(data(2))
        call match_tau_column(nset,data)
        vetoBPM = .false.

        goto 871

     end if


     call tune(data)
     call convert_data_from_pi_matrix(data)
     call calculate_with_known_spacing(data)

     if (.not. silentMode) then
        call plots(data, 2)         !Use plot_it2
     end if
     !Check if user has vetoed a detector:
     if (vetoBPM) then
        !Parse user input and remove vetoed detector(s)
        call veto(data, inputPasser)
        goto 567

     end if

!Write the output files
     call output (data)

     call fftOut(data)

     if (.not. noSPos) then
        call cesrv_out()
     end if

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





