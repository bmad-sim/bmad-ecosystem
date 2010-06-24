module mia_plot

  !This module contains the plotting subroutines for MIA.
  use quick_plot
  use mia_types
  use mia_matrixops

  logical :: windowOpen = .false.     !Whether a plotting window is open
contains
  subroutine plots(data, plot)     
    !Routine can be used to plot data from different  
    !matrices. Can use digital display but also writes 
    !a postscript that can be printed for a hard copy

    implicit none
    type(data_set) :: data(*)     !Data set
    character(1) ans              !User input (any character to continue)
    INTEGER :: id = 0             !ID for page--arbitrary number
    Integer :: plot               !To use plot_it or plot_it2
    logical :: printit = .false.  !To print or not to print
    logical :: plot_more          !To plot or not to plot
    character(50) :: input        !User input for plotting options
    integer :: col                !Keeps track of columns to plot
    integer :: graph_set, &       !Which set of graphs to use for debug_plots
         old_set                  !Previous value of graph_set
    logical :: gif, gifOpen       !gif is true if a GIF file is to be plotted
    character(15) :: gifName      !Filename for a gif file
    integer ::numGif              !Number of GIF files (used for filenames)
    if (postScript) then
       gif = .true.
    else
       gif = .false.
    end if
    col = 1
    plot_more = .true.
    windowOpen = .false.    
    graph_set = 5
    old_set = 3
    gifOpen = .false.
    numGif=1
54  format (i)

    do while(plot_more)
       if (.not. windowOpen) then
          !Following calls open the digital displays
          if (gif) then
             write (gifName,54), numGif
             if (numGif < 10) then
                gifName = "0" // adjustl(gifName)
             end if
             gifName = "mia" // adjustl(gifName) // ".gif"
             numGif = numGif+1
             call qp_open_page ("GIF-L", id, 600.0_rp, 550.0_rp, "POINTS",&
                  plot_file=gifName)
             gifOpen = .true.
          else
             call qp_open_page ("X", id, 600.0_rp, 550.0_rp, "POINTS")
          end if
          call qp_set_page_border (0.02_rp,0.02_rp,0.035_rp,0.035_rp,"%PAGE")
          call qp_set_margin (0.07_rp, 0.01_rp, 0.01_rp, 0.05_rp, "%PAGE")
          call qp_draw_text (" Orbit MIA ", 0.5_rp, 1.0_rp, "%PAGE", "CT") 
          windowOpen = .true.
       endif


       if (debugMode) then
          !Plot first, then ask questions.
          call qp_draw_text (" Orbit MIA ", 0.5_rp, 1.0_rp, "%PAGE", "CT") 
          if (graph_set<7) then
             call debug_plots(data,graph_set, col)
          else
             call plot_it2(data)
          end if
176       write (*,*), "Plotting options:"
          write (*,'(4x,a)')  "1 - Beta Ratios (A Mode)", &
               "2 - Beta Ratios (B Mode)",&
               "3 - Pi eigenvectors and SVD eigenvalues (A Mode file)",&
               "4 - Pi eigenvectors and SVD eigenvalues (B Mode file)",&
               "5 - Tau eigenvectors and their fft (A Mode file)",&
               "6 - Tau eigenvectors and their fft (B Mode file)",&
               "7 - Beta and coupling",&
               "8 - Phase advance and coupling",&
               "col n - Plot column n of the current set",&
               " help - Print this message", &
             " print - Toggle generating gif file of the current plot on/off",&
           " veto [#] - Remove detector # (ex 12W) from data and redo analysis"
          write (*,*), "Type anything else to continue"
          accept "(a)", input

          old_set = graph_set
          call lowerCase(input)
          select case(input)
          case('help')
             goto 176
          case('q')
             plot_more = .false.
          case('quit')
             plot_more = .false.
          case('print')
             !Generate a GIF file
             if (gif) then
                gif=.false.
             else
                gif = .true.
             end if
          case ('1')
             graph_set = 1
          case ('2')
             graph_set = 2
          case ('3')
             graph_set = 3
          case ('4')
             graph_set = 4
          case ('5')
             graph_set = 5
          case ('6')
             graph_set = 6
          case ('7')
             phase = .false.
             graph_set = 7
          case ('8')
             phase = .true.
             graph_set = 8             
          case default
             select case(input(1:4))
             case('col ')
12              format (i)
                !Read column number
                read (input(5:len_trim(input)),*), col
             case('veto')
                vetoBPM = .true.
                inputPasser = input
                plot_more = .false.
             case default
                plot_more = .false.
             end select
          end select
          !Reset column to 1 when plotting something different
          !Also reset it once the last columns are reached
          if (old_set .ne. graph_set) then
             col = 1
          else if (col+1 > 2*NUM_BPMS .and. graph_set>4 .or. &
               col+3 > 2*NUM_BPMS .and. graph_set <5) then
             col = 1
          end if
          !*End debug mode case
       else
          !Default plots
          call plot_it2 (data)
177       write (*,*), "Plotting options:"
          write (*,'(2x,a)')  " phase - Plot phase advance", &
               " beta - Plot betas"," help - Print this message", &
               " print - Toggle making gif file of the current plot on/off", &
               " veto [#] - Remove detector # (ex 12W) from data and redo analysis"
          write (*,*), "Type anything else to continue"
          accept "(a)", input
         
          call lowerCase(input)
          
          Select case(input(1:5))
          case ('phase')
             phase = .true.
          case('print')
             if (gif) then
                gif = .false.
             else
                gif = .true.
             end if
          case('beta')
             phase = .false.
          case('help')
             goto 177     !Go reprint the options
          case default
             select case(input(1:4))
             case('veto')
                !Just pass input along;
                !subroutines to remove dets are called in the main program.
                vetoBPM = .true.
                inputPasser = input
             end select
             plot_more = .false.
          end select

       end if

       !GIF files are written when the page is closed
       !If an X window is being used instead, then the page is
       !cleared for the next set of plots.
       if (gifOpen) then
          call qp_close_page
          gifOpen = .false.
          windowOpen = .false.
       else
          if (gif) then
             call qp_close_page
             windowOpen = .false.
             gifOpen = .false.
          else
             call qp_clear_page
          end if
       end if
    end do

    !Close plot window
    if (windowOpen) then
       call qp_close_page
       windowOpen = .false.
    endif

  end subroutine plots

  subroutine debug_plots(data,graph_set,col)
    !
    !This subroutine provides plots of SVD results and other intermediates
    !for use in diagnosing problems with MIA or input files. For example,
    !looking at the graphs of the eigenvalues and tau matrix will tell you
    !if noise has drowned out your signal.
    !
    !*Programmers beware: there is danger ahead. Please wear eye protection.
    !

    type(data_set) :: data(*)          !Data from file
    integer :: i, &                    !Counter
         xdiv, &                       !Divisions of x axis
         graph, &                      !Choice to graph
         arr_length, icolumn, count, nset, n_graphs,&
         graph_set, &                  !User choice for sets of data to plot
         format_num, &
         col, &                        !Column of pi or tau matrix to plot
         A_file, B_file                !Which file number is A or B mode
    real(rp) xlength, &                !Length of x axis (# values)
         miny, maxy                    !Min and max y values
    real(rp), allocatable :: xcoord(:), & !X coordinates to plot
         ycoord(:,:), &                  !Y coordinates to plot
         sPos(:)                       !Array of S positions for plotting
    character(40) :: title
    character(40),allocatable :: titl(:)    !Titles
    logical :: graph_more              !To graph more
    integer :: ix, iy, ix_tot, iy_tot  !Number of plots in X and y directions
    real(rp), allocatable :: b(:), &   !Number of BPMs (for plotting)
         nt(:), &                      !Number of turns
         lam_log(:)                    !Log of lambda (for plotting)
    real(rp) :: xmin, xmax             !X min and max for qplot
    character(10) :: tempChar(4)       !Holds column # to use in plot titles
    logical :: sPoss
    allocate (b(2*NUM_BPMS))
    allocate (nt(NUM_TURNS))
    allocate (sPos(NUM_BPMS))
    nset = 2                 !Make nset a global variable
!    A_file = data_struc%set_num_a
!    B_file = data_struc%set_num_b

A_file = 1
B_file = 2

    do i=1,2*NUM_BPMS
       b(i) = i
    end do
    do i=1,NUM_TURNS
       nt(i) = i
    end do

    ix_tot = 1
    ix = 1                         
    graph_more = .true.

    graph = 1
    DO WHILE (graph_more)
       IF (ALLOCATED(xcoord)) DEALLOCATE(xcoord)  
       IF (ALLOCATED(ycoord)) DEALLOCATE(ycoord)
       if (ALLOCATED(titl)) DEALLOCATE(titl)
       sPoss=.false.
       if (graph_set > 4) then
          n_graphs = 4
          arr_length = NUM_TURNS
          allocate(xcoord(arr_length))
          xcoord = nt(1:NUM_TURNS)
          allocate (titl(4))
       else if (graph_set > 2) then
          n_graphs  = 5 
          arr_length = 2*NUM_BPMS
          allocate(xcoord(arr_length))
          xcoord = b
          allocate (titl(5))
       else
          n_graphs  = 5
          arr_length = NUM_BPMS
          allocate(xcoord(arr_length))
          xcoord = b(1:NUM_BPMS)
          allocate (titl(5))
!          sPoss = .true.
       end if
       allocate(ycoord(n_graphs,arr_length))
       xlength = arr_length
       iy_tot = n_graphs
       iy = iy_tot                           !Graph position 1 is the bottom 

       call qp_set_box(ix, iy, ix_tot, iy_tot)

                      !!!!ZOMBIES INVADE!!!!!
             !!!They bring with them six plotting options!!!
       select case (graph_set)
       case(1)
          titl=(/"\gg\u2\d \gb Ratio - A Mode","\gg\u2\d\gb Ratio - B Mode",&
               "<J\dt\u>\gD\u2\d - A Mode","<J\dt\u>\gD\u2\d - B Mode",&
               "Angle - A Mode"/)
          do i=1, NUM_BPMS
             ycoord(1,i) = data_struc%loc(i)%a%gam2_beta_ratio
             ycoord(2,i) = data_struc%loc(i)%b%gam2_beta_ratio
             ycoord(3,i) = data_struc%loc(i)%a%magnitude2(1)
             ycoord(4,i) = data_struc%loc(i)%b%magnitude2(1)
             ycoord(5,i) = atan(data_struc%loc(i)%a%ratio(1))
          end do
       case (2)
          titl=(/'\gg\u2\d\gb Ratio - A Mode','\gg\u2\d\gb Ratio - B Mode',&
               '<J\dt\u>\gD\u2\d - A Mode','<J\dt\u>\gD\u2\d - B Mode',&
               'Angle - B Mode'/)
          do i=1, NUM_BPMS
             ycoord(1,i) = data_struc%loc(i)%a%gam2_beta_ratio
             ycoord(2,i) = data_struc%loc(i)%b%gam2_beta_ratio
             ycoord(3,i) = data_struc%loc(i)%a%magnitude2(2)
             ycoord(4,i) = data_struc%loc(i)%b%magnitude2(2)
             ycoord(5,i) = atan(data_struc%loc(i)%b%ratio(2))
          end do
       case (3)
!          Print *, 'Welcome to case 3. Watch out for frogs.
          do i=1,4
             call intToChar((col+i-1),tempChar(i))
          end do
          titl = (/('A Mode \gP Column       '//trim(tempChar(1))), &
               ('A Mode \gP Column       '//trim(tempChar(2))),&
               ('A Mode \gP Column       '//trim(tempChar(3))),&
               ('A Mode \gP Column       '//trim(tempChar(4))),&
               'Log of Eignenvalues, A Mode'/)
          do i=1, 2*NUM_BPMS
             ycoord(1,i) = data(A_file)%pi_mat(i,col)
             ycoord(2,i) = data(A_file)%pi_mat(i,col+1)
             ycoord(3,i) = data(A_file)%pi_mat(i,col+2)
             ycoord(4,i) = data(A_file)%pi_mat(i,col+3)
             ycoord(5,i) = log10(data(A_file)%lambda(i))
!             Print *, "Lambda: ", data(A_file)%lambda(i)
          end do
          col = col+4
       case(4)
          do i=1,4
             call intToChar((col+i-1),tempChar(i))
          end do
          titl = (/('B Mode \gP Column      '//trim(tempChar(1))), &
               ('B Mode \gP Column      '//trim(tempChar(2))),&
               ('B Mode \gP Column      '//trim(tempChar(3))),&
               ('B Mode \gP Column      '//trim(tempChar(4))),&
               'Log of Eignenvalues, B Mode'/)
          do i=1, 2*NUM_BPMS
             ycoord(1,i) = data(B_file)%pi_mat(i,col)
             ycoord(2,i) = data(B_file)%pi_mat(i,col+1)
             ycoord(3,i) = data(B_file)%pi_mat(i,col+2)
             ycoord(4,i) = data(B_file)%pi_mat(i,col+3)
             ycoord(5,i) = log10(data(B_file)%lambda(i))
          end do
          col = col+4
       case(5)
          do i=1,2
             call intToChar((col+i-1),tempChar(i))
          end do
          titl = (/('\gt Column         '//trim(tempChar(1))), &
               ('FFT Spectrum Column'//trim(tempChar(1))),&
               ('\gt Column         '//trim(tempChar(2))),&
               ('FFT Spectrum Column'//trim(tempChar(2))),&
               ('he   l0000lo'//trim(tempChar(2)))/)
          do i=1, NUM_TURNS
             ycoord(1,i) = data(A_file)%tau_mat(i,col)
             ycoord(2,i) = data(A_file)%spectrum(i,col)
             ycoord(3,i) = data(A_file)%tau_mat(i,col+1)
             ycoord(4,i) = data(A_file)%spectrum(i,col+1)
          end do
          col = col+2
       case(6)
          do i=1,2
             call intToChar((col+i-1),tempChar(i))
          end do
          titl = (/('\gt Column         '//trim(tempChar(1))), &
               ('FFT Spectrum Column'//trim(tempChar(1))),&
               ('\gt Column         '//trim(tempChar(2))),&
               ('FFT Spectrum Column'//trim(tempChar(2))),&
               ('y0000000 00o'//trim(tempChar(2)))/)
          do i=1, NUM_TURNS
             ycoord(1,i) = data(B_file)%tau_mat(i,col)
             ycoord(2,i) = data(B_file)%spectrum(i,col)
             ycoord(3,i) = data(B_file)%tau_mat(i,col+1)
             ycoord(4,i) = data(B_file)%spectrum(i,col+1)
          end do
          col = col+2
       end select


       do graph = 1,n_graphs
!          tempTitle = trim(titl(graph))
          call min_max_y(ycoord(graph,:), miny, maxy, ycoord(graph,:),&
               arr_length)
          if (sPoss) then
             call sortsPos(xcoord, ycoord(graph,:))
          end if
          xmin = minval(xcoord)
          xmax = maxval(xcoord)
          call calcAxes(xmin,xmax,xdiv)
          !Draws a single graph
          call qp_save_state(.true.)
          call qp_set_box (ix, iy, ix_tot, iy_tot)
          call qp_set_symbol_attrib (star5_filled$, height = 5.0_rp)
          call qp_set_axis ("X", xmin, xmax, div=xdiv)
          call qp_calc_and_set_axis("Y", miny, maxy, 4, 6, "GENERAL")
          call qp_draw_graph (xcoord,ycoord(graph,:),title=titl(graph))
          if (sPoss) then
             !Line at IP (s=0):
             call qp_draw_line(0.0_rp,0.0_rp, miny,maxy,width=5, color = 2,&
                  style=4)
             !Line at L3 (if plot goes that far):
             if (maxval(xcoord) > ip_L3) then
                call qp_draw_line(ip_L3,ip_L3,miny,maxy,width=5,&
                     color=3,style=4)
             end if
          end if
          call qp_restore_state
          iy = iy-1
          ix = 1
       end do
       graph_more = .false.
    end do

    deallocate (b)
    deallocate (nt)
    if (allocated(xcoord)) deallocate (xcoord)
    if (allocated(ycoord)) deallocate (ycoord)
  end subroutine debug_plots

  subroutine min_max_y(vector, miny, maxy, ycoord, arr_length)

    !
    !Finds min and max y values from a given vector to be plotted.
    !Also assigns the vector to ycoord.
    !Used in plotting graphs.
    !

    real(rp) :: vector(:), &       !Some vector passed into the subroutine
         miny, maxy, &             !Find the min and max values of the vector
         max_factor, min_factor, & !Used to find order of magnitude to
                                   !round min, max properly
         ycoord(:)                 !Takes the values of vector
    integer :: arr_length

    !Change factors to integers?
    ycoord(:) = vector(1:arr_length)
    miny = minval(ycoord(1:arr_length))
    maxy = maxval(ycoord(1:arr_length))
    !Rounding
    !Finds the power to multiply maxy by to get 4 figures
    !Inintalize factors to 0 in case if condition is not met. 
    max_factor = 0
    min_factor = 0
    !Use absolute value to avoid neg. numbers
    if (maxy /= 0)    max_factor = log10(abs(maxy)) 
    if (miny /= 0)    min_factor = log10(abs(miny))

    min_factor = floor(min_factor)
    call rounding(max_factor)
!Adjust this to apply to a wider range (ex 10 vs. 1000)
    miny = floor(miny*10/(10**min_factor))*(10**min_factor)/10
    maxy = ceiling(maxy*10/(10**max_factor))*(10**max_factor)/10

  end subroutine min_max_y

  subroutine rounding(factor)

    !Rounds some number.
    !Used in min_max_y

    real(rp) :: factor   !Some number, called factor because it is only used
    !for min_max_y's factors
    if (factor < 0) then
       factor = ceiling(factor)
    else
       factor = floor(factor)
    endif

  end subroutine rounding

  subroutine plot_it2 (data)
    !
    !Depending upon settings in the plots() subroutine, plots will either be
    !generated on the display or in a postscript file. By default, 
    !beta and cbar are plotted. There is an option to replace beta
    !with phase advance and replot. 
    !Other plotting options are available through the debug option in 
    !subroutine debug_plots().
    !

    type(data_set) :: data(*)          !Data from file
    integer :: i, &                    !Counter
         xdiv, &                       !Divisions of x axis
         graph, &                      !Choice to graph
         arr_length, icolumn, count, nset
    real(rp) xlength, &                !Length of x axis (# values)
         miny, maxy                    !Min and max y values
    real(rp), allocatable :: xcoord(:), &  !X coordinates to plot
         ycoord(:), &                  !Y coordinates to plot
         sPos(:)                       !Array of S positions for plotting
    real(rp) :: sf2                    !Factor
    character(80) title, titl          !Titles
    logical :: graph_more              !To graph more
    integer :: ix, iy, ix_tot, iy_tot  !Number of plots in X and y directions
    real(rp), allocatable :: b(:), &   !Number of BPMs (for plotting)
         phi(:), &                     !Phi (for x-value plotting)
         nt(:), &                      !Number of turns
         lam_log(:)                    !Log of lambda (for plotting)
    real(rp) :: xmin, xmax             !X min and max for qplot
    logical :: badFile, badColumn      !Errors in user input
    real (rp):: sMax   !Should be global?
    logical :: sPoss
    sPoss = .true.
    sMax = endLoc/2

    allocate (b(NUM_BPMS))
    allocate (phi(NUM_BPMS))
    allocate (nt(NUM_TURNS))
    allocate (sPos(NUM_BPMS))

    do i=1,NUM_BPMS
       b(i) = i
    end do
    do i=1,NUM_TURNS
       nt(i) = i
    end do

    ix_tot = 1
    ix = 1                         

!    graph_more = .true.

    iy_tot = 5        !Graph 5 plots by default
    iy = iy_tot                             !Graph position 1 is the bottom 
    DO graph = 1, 5
       IF (ALLOCATED(xcoord)) DEALLOCATE(xcoord)  
       IF (ALLOCATED(ycoord)) DEALLOCATE(ycoord)

       arr_length = NUM_BPMS
       allocate(xcoord(NUM_BPMS))
       allocate(ycoord(NUM_BPMS))
       xlength = NUM_BPMS
       xdiv = 16
       call qp_set_box(ix, iy, ix_tot, iy_tot)

       !
       !Formats
       !
36     format('Phase Advance - A Mode')
37     format('Phase Advance - B Mode')
38     format('Gamma**2 Beta - A Mode')
39     format('Gamma**2 Beta - B Mode')
43     format('Inv Gamma Cbar (1,1)')  
44     format('Inv Gamma Cbar (1,2)')
45     format('Inv Gamma Cbar (2,2)')
  
       if (graph == 1 .and. phase) then
          ycoord(:) = data_struc%loc(:)%a%phi
          call phase_bound(ycoord,miny,maxy)
          write(titl,36)
       else if (graph == 2 .and. phase) then
          ycoord(:) = data_struc%loc(:)%b%phi
          call phase_bound(ycoord,miny,maxy)
          write(titl,37)
       else if (graph == 1 .and. .not. phase) then
          call min_max_y(data_struc%loc(:)%a%gam2_beta, miny, maxy, &
               ycoord(:), arr_length)
          write(titl,38)
       else if (graph == 2 .and. .not. phase) then
          call min_max_y(data_struc%loc(:)%b%gam2_beta, miny, maxy, &
               ycoord(:), arr_length)
          write(titl,39)
       else if (graph == 3) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(1,1), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,43)
       else if (graph == 4) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(1,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,44)
       else if (graph == 5) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(2,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,45)
       end if
       if (sPoss) then
          call sortsPos(xcoord, ycoord)
       else
          xcoord = b(1:NUM_BPMS)
       end if
       xmax = maxval(xcoord)
       xmin = minval(xcoord)
       call calcAxes(xmin,xmax,xdiv)
       !Draws graph
       call qp_save_state(.true.)
       call qp_set_box (ix, iy, ix_tot, iy_tot)
       call qp_set_symbol_attrib (star5_filled$, height = 5.0_rp)
       call qp_set_axis ("X", xmin, xmax, div = xdiv)
       call qp_calc_and_set_axis("Y", miny, maxy, &
            4, 6, "GENERAL")
       call qp_draw_graph (xcoord, ycoord(:),title = titl)
       !Removed line at IP because sPos no longer goes negative.
       !Line at IP (s=0):
!       call qp_draw_line(0.0_rp,0.0_rp, miny,maxy,width=5, color = 2,&
!            style=4)
       !Line at L3 (if plot goes that far):
       call qp_draw_line(ip_L3,ip_L3, miny,maxy,width=5, color = 3,&
            style=4)

       call qp_restore_state

       IF (iy <= 1) THEN
          graph_more = .false.
       ENDIF

       iy = iy-1
       ix = 1

    enddo

    deallocate (b)
    deallocate (phi)
    deallocate (nt)
    if (allocated(xcoord)) deallocate (xcoord)
    if (allocated(ycoord)) deallocate (ycoord)
    if (allocated(lam_log)) deallocate (lam_log)
  end subroutine plot_it2

  subroutine arrange_phi(phi, xcoord)

    integer :: i, numEast
    real(rp) :: phi(:), xcoord(:)
    real(rp), allocatable :: phi_old(:)
    allocate(phi_old(NUM_BPMS))
    phi_old(:) = phi(:)

    do i=1, NUM_BPMS
       if (.not. data_struc%proc((NUM_BPMS-i)+1)%is_west) then
          !Place values in phi in reverse order to put
          !east BPMs first
          numEast = i
          phi(i) = phi_old(NUM_BPMS-i+1)
  !        xcoord(i) = data_struc%proc(NUM_BPMS-i+1)%number
       else
          phi(i) = phi_old(i-numEast)
 !         xcoord(i) = data_struc%proc(i-numEast)%number
       endif
    enddo

    !Reverse order of east BPMS (they are in ascending order)
    phi_old(:) = phi(:)
    do i=1, numEast
       phi(i) = phi_old(numEast-i+1)
!       xcoord(i) = data_struc%proc(NUM_BPMS-numEast+i)%number
    enddo

    deallocate(phi_old)
  end subroutine arrange_phi

  subroutine sortsPos(xcoord,ycoord)
    !
    !Set S position for the x coordinates when plotting.
    !Also sets a division point for a large gap between east and west
    !bpms. Otherwise, plots by S position can be difficult to read.
    !Coordinates are reorganized so east BPMs come before west.
    !

    integer :: i,j, break, numEast
    real(rp) :: xcoord(:),ycoord(:), div
    real(rp), allocatable :: tempy(:),sPos(:)
    logical :: east !Becomes true when east BPMs are found
    logical :: twice !True if xcoord contains x*NUM_BPMS elements
    allocate (tempy(num_bpms))
    allocate(sPos(num_bpms))

    numEast=0
    do i=1, num_bpms
       xcoord(i) = data_struc%proc(i)%sPos
!       sPos(i) = data_struc%proc(i)%sPos
!       if ( .not. data_struc%proc(i)%is_west) then
!          if (.not.east) then
!             break = i-1 !Index of last west detector
!             east = .true.
!          end if
!          numEast = numEast + 1
!          sPos(i) = sPos(i) - endLoc
!          xcoord(i-break) = sPos(i)
!          tempy(i-break) = ycoord(i)
!       end if
    end do

!    do i=1,break
!       xcoord(numEast+i) = sPos(i)
!       tempy(numEast+i) = ycoord(i)
!    end do

!    ycoord = tempy
!    call quickSortpos(xcoord, ycoord)
    call slowSort(xcoord, ycoord)
  end subroutine sortsPos

  recursive subroutine quickSortpos(xcoord, ycoord)
    !Uses the recursive quicksort algorithm. Doens't work
    !quite as intended... Not in use until the bugs have been
    !corrected.

    real(rp) :: xcoord(:),ycoord(:) !The x and y coordinate vectors to sort
    integer :: i, j                 !Indices for quicksort algorithm
    real(rp) :: pivotVal            !Value of the pivot (xcoord)
    logical :: lowerFound, higherFound ! Has it found a higher and lower value
                                      ! to flip across the pivot?
    real(rp) :: tempx, tempy
    logical :: keepGoing

    pivotVal = xcoord(size(xcoord)/2)
    i=1
    j=size(xcoord)
    keepGoing = .false.

    Print *, "Pivot: ", pivotVal, "i", i,"j", j

    do while (i<j)
       lowerFound = .false.
       higherFound = .false.

       do while (i<j .and. .not. lowerFound)
          if (xcoord(i) > pivotVal) then
             lowerFound = .true.
          else
             i=i+1
          end if
       end do

       do while (j> i .and. .not. higherFound)
          if (xcoord(j) < pivotVal) then
             higherFound = .true.
          else
             j = j-1
          end if
       end do
       if (lowerFound .and. higherFound) then
!          Print *, "Switching ", ycoord(i), " with ", ycoord(j)
          call swap(xcoord(i), xcoord(j))
          call swap(ycoord(i), ycoord(j))
          keepGoing = .true.
       end if
    end do

1211 continue
    if (keepGoing) then
       Print *, "Left side"
       call quickSortpos(xcoord(1:i), ycoord(1:i))
       Print *, "Right side"
       call quickSortpos(xcoord(i+1:size(xcoord)), ycoord(i+1:size(ycoord)))
    end if

  end subroutine quickSortpos

  subroutine slowSort(xArr, yArr)
    real(rp) :: xArr(:), yArr(:)
    integer :: i,j

    do i=1, size(xArr)
       do j=1, size(xArr)
          if (xArr(j) > xArr(i) ) then
             call swap(xArr(i),xArr(j))
             call swap(yArr(i),yArr(j))
          end if


       end do
    end do

  end subroutine slowSort


  subroutine swap (a,b)
    real(rp) :: a, b, temp

    temp = a
    a = b
    b = temp
  end subroutine swap

  subroutine rearrange(arr,break)

    real(rp) :: arr(:)
    real(rp),allocatable :: temp(:)
    integer :: break, numEast, i

    allocate(temp(size(arr)))
    numEast = num_bpms-break

    do i=break+1,num_bpms
       temp(i-break+1) = arr(i)
    end do


    do i=1,break
       temp(i+numEast) = arr(i)
    end do

    arr = temp
  end subroutine rearrange

  subroutine superSize(half, double)

    real(rp) :: half(:),double(:)
    integer :: i
    do i=1,num_bpms
       Print *, "I am:", i
       double(2*i-1) = half(i)
       double(2*i) = half(i)
    end do

  end subroutine superSize

  subroutine lowerCase(string)
    !Modified from Fortran 90/95 for Scientists and Engineers
    !Written by S.J. Chapman
    
    !Shifts a string to all lowercase letters
    
    character(5) :: string
    integer :: i
    
    do i=1, 5
       if ( LGE(string(i:i), 'A') .and. LLE(string(i:i),'Z')) then
          string(i:i) = ACHAR( IACHAR ( string(i:i)) - 32)
       end if
    end do
  end subroutine lowerCase

  subroutine intToChar(int,char)
    !
    !Converts an integer into a character string.
    !
    integer :: int
    character(3) :: char

    write (char,'(i3)'),int
    char = trim(char)

  end subroutine intToChar


  subroutine phase_bound(array, miny, maxy)
    !
    !Adjusts phase range from -pi to +pi for plotting
    !
    real(rp) :: array(:), miny, maxy
    integer :: i

    do i=1, NUM_BPMS
       array(i) = modulo(array(i)+pi,2*pi)-pi
    end do
    miny = minval(array)
    maxy = maxval(array)
  end subroutine phase_bound

  subroutine calcAxes(xmin,xmax,xdiv)
    !
    !Calculates axes labels so that they are multiples of
    !5,25,50,100, etc.
    !
    real (rp) :: xmin,xmax
    integer :: xdiv, range,modterm,modx, interv
    real :: interval

    interval = 6 !Try six divisions first
    if (xmax > 300) then
       xmin = floor(xmin/100)*100
       xmax = ceiling(xmax/100)*100
    else  
       xmin = floor(xmin/10)*10
       xmax = ceiling(xmax/10)*10
    end if

   !Move xmin from 1 to 0 for markers to be at 0,5,etc instead of 1,6,...
    if (mod(xmin,10.0_rp) == 1) then
       xmin = xmin-1
    end if
    range = ceiling(xmax-xmin)
    xdiv=ceiling(range/interval)
!    if (mod(xdiv,5) == 0) then
!       goto 511
!    end if

    !For xmax < 5, this algorithm must be modified.
    if (range .ge. 500) then
       interv = 100
    else if (range .ge. 150) then
       interv = 50
    else if(range .ge. 100) then
       interv = 25
    else if (range .ge. 50) then
       interv = 10
    else
       interv = 5
    end if
    xdiv = range / interv

!    modx = mod(xdiv,modterm)
!    if (modx>0) then
!       if (modx .ge. (modterm/2)) then
!          xdiv = xdiv + (modterm-modx)
!       else
!          xdiv = xdiv - modx
!          interval = interval + 1
!       end if
!       xmax = xdiv * interval+xmin
!       xdiv=interval
!    end if

511 continue
!    xdiv = interval
  end subroutine calcAxes

end module mia_plot
