module mia_plot

  !This module contains the plotting subroutines for MIA.
  use quick_plot
  use mia_types

  logical :: windowOpen = .false.     !Whether a plotting window is open
contains
  subroutine plots(data, plot, iset)     
    !Routine can be used to plot data from different  
    !matrices. Can use digital display but also writes 
    !a postscript that can be printed for a hard copy

    implicit none
    type(data_set) :: data(*)     !Data set
    character(1) ans              !User input (any character to continue)
    INTEGER :: id = 0             !ID for page--arbitrary number
    Integer :: plot               !To use plot_it or plot_it2
    Integer :: iset               !File number
    logical :: printit = .false.  !To print or not to print
    logical :: plot_more          !To plot or not to plot

    call logic_get( 'Y', 'N', 'Plot data? (Y/N) ', plot_more)

    !Only opens a new page on the first call of the subroutine.
    if (plot_more .and. windowOpen == .false.) then
       !Following calls open the digital displays
       call qp_open_page ("X", id, 600.0_rp, 470.0_rp, "POINTS")
       call qp_set_page_border (0.02_rp, 0.02_rp, 0.035_rp, 0.035_rp, "%PAGE")
       call qp_set_margin (0.07_rp, 0.01_rp, 0.01_rp, 0.05_rp, "%PAGE")
       windowOpen = .true.
    endif

    do while(plot_more)
       call qp_draw_text (" Orbit MIA ", 0.5_rp, 1.0_rp, "%PAGE", "CT") 
       if (plot == 2) then
          call plot_it2 (data)
       else
          !Calls plot_it by default
          call plot_it(data(iset), iset)
       endif             

       !    write (*, "(a)") " Hit any key to continue: "
       !    accept "(a)", ans
       !    call qp_close_page

       call logic_get('Y', 'N', 'Plot more data? (Y/N) ', plot_more) 

       if (plot_more) then
          call logic_get('P','C',' (P)rint a plot or (C)ontinue',printit)
          if (printit) then
             !       call qp_close_page 
             call qp_open_page ("PS-L")    !Generates a PS file
             if (plot == 2) then
                call plot_it2(data)
             else
                call plot_it(data(iset), iset)
                !Generates the plot in ps file
             endif
             call qp_close_page    !Closes and names PS file as quick_plot.ps
          end if
       endif
       call qp_clear_page
    end do

    if (windowOpen .and. plot == 2) then
       call qp_close_page
       windowOpen = .false.   
    endif

  end subroutine plots

  subroutine plot_it (data, iset)

    IMPLICIT NONE

    type(data_set)data              !Data from file
    integer ::  xdiv, &             !# Divisions in x
         graph,&                    !What to plot
         istat,in4get1, &           !Status and get input
         ix, iy, &                  !Place in x and y in the window
         ix_tot, iy_tot,&           !Total number of divisions in x and y
         count, i, &                !Counters
         divisions, &
         icolumn, arr_length, iset
    real(rp) xlength, &             !Length of X axis (# data points)
         miny, maxy                 !Min and max y values
    real(rp), allocatable :: xcoord(:), &  !Xcoordinates to be plotted
         lam_log(:), &              !Log of lambda (for plotting)
         ycoord(:)                  !Y coordinates to be plotted
    character(80) title, titl       !Titles of graph
    logical graph_more              !To graph more
    real(rp), allocatable ::  b(:), & !X values for plotting (num. BPMs)
          nt(:)                     !X values for plotting (num. turns)

    allocate (b(2*NUM_BPMS))
    allocate (nt(NUM_TURNS))


    do i=1,NUM_TURNS
       nt(i) = i
    end do

    do i=1,2*NUM_BPMS
       b(i) = i
    end do

    graph_more = .true.
    ix_tot = 1                           !Graphs are stacked vertically,
    ix = ix_tot                          !Change ix_tot to plot side-by-side

    Print *, "How many graphs do you want?"
    PRINT *, "Enter a number between 0 and 5"
    Read *, iy_tot

    If (iy_tot>5) then
       print *, "Changing number of graphs to 5 (maximum)."
       iy_tot = 5
    endif
    if (iy_tot<=0) then
       print *, "Number of graphs 0 or less--nothing is being plotted."
       graph_more = .false.
    endif

    iy = iy_tot                             !Graph position 1 is the bottom 

    DO WHILE (graph_more)        !Plot graphs as long as the user wants more.

       IF (ALLOCATED(xcoord)) DEALLOCATE(xcoord)  
       IF (ALLOCATED(ycoord)) DEALLOCATE(ycoord)

98     write (*,'(1x,a)') "For Pi, enter 1","For Frequency Peak, &
            enter 2",  "For Lambda, enter 3", "For Spectrum, enter 4", &
            "For Tau, enter 5"
       write (*, "(a)") " Enter choice "        
       accept "(i)", graph

       if (graph<=0 .or. graph>5) then
          PRINT *, "Enter a number between 1 and 5."
          goto 98
       else if (graph ==1 .or. graph > 3) then
          !Column is needed for choices 1, 4, and 5
99        write (*, "(a)") " Which column should be plotted? "        
          accept "(i)", icolumn
          if (icolumn < 1 .or. icolumn > 2*NUM_BPMS)then
             Print *, "Column should be between 1 and ", 2*NUM_BPMS
             goto 99
          endif
       endif

       call qp_set_box (ix, iy, ix_tot, iy_tot)

       !Allocate or assign variables common to all options if input is valid:
       if (graph<=3) then
          arr_length = 2*NUM_BPMS
          allocate(xcoord(arr_length))
          allocate(ycoord(arr_length))
          xlength = 2.0*NUM_BPMS
          xdiv = 2*NUM_BPMS
          xcoord = b
       else
          arr_length = NUM_TURNS
          allocate(xcoord(arr_length))
          allocate(ycoord(arr_length))
          xdiv = 16
          xcoord = nt
          xlength = NUM_TURNS
       endif

       !Graph stuff
       if (graph == 1) then
          call min_max_y(data%pi_mat(:,icolumn), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,22) iset, icolumn, data%shortName
22        format('PI - Iset =',i2,'  Column = ',i3, '  File = ', a)
       else if (graph == 2) then
          !nmax is an integer; all others are real(rp)
          !Can't use a generic statement because of the type of nmax??
          !call min_max_y(data%fr_peak(:))
          miny = minval(data%fr_peak(:))
          maxy = maxval(data%fr_peak(:))
          ycoord(:) = data%fr_peak(:)
          write(titl,23) iset, data%shortName     
23        format('FR_PEAK - Iset =',i2,'  File = ', a)
       else if (graph == 3) then
          !Logarithmic scale option on qplot is not available yet
          !(See documentation)
          !Plots lambda on a logarithmic scale
          IF (ALLOCATED(lam_log)) DEALLOCATE(lam_log)    
          allocate (lam_log(arr_length))
          count = 1
          do while (count<=arr_length)
             lam_log(count) = log10(data%lambda(count))
             count = count+1
          enddo
          !LAPACK sorts lambdas
          !          call sort_l(lam_log, 2*NUM_BPMS)
          call min_max_y(lam_log(:), miny, maxy, ycoord(:), arr_length)   
          write(titl,24) iset,  data%shortName
24        format('Lambda (log10) - Iset =',i2,'  File = ', a)
       else if (graph == 4) then
          call min_max_y(data%spectrum(:,icolumn), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,25) iset,icolumn, data%shortName
25        format('Spectrum - Iset =',i2,'  Column = ',i2, '  File = ', a)
       else if (graph == 5) then
          call min_max_y(data%tau_mat(:,icolumn), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,21) iset, icolumn, data%shortName
21        format('TAU - Iset =',i2,'  Column = ',i3, '  File = ', a)
       end if


       divisions = ceiling(10.0 - (iy_tot/10.0))
       !Draws graph
       call qp_save_state(.true.)
       call qp_set_box (ix, iy, ix_tot, iy_tot)
       !  call qp_set_graph_attrib (draw_grid = .true.)
       call qp_set_symbol_attrib (star5_filled$, height = 5.0_rp)
       call qp_set_axis ("X", 0.0_rp, xlength, div = xdiv)
       !call qp_set_axis ("Y", miny, maxy, places = 4, draw_numbers = .true., &
       !div = divisions)
       call qp_calc_and_set_axis("Y", miny, maxy, 4, 6, "GENERAL")

       !**********************************************************
       !There are issues with overlapping text. This may be an
       !issue in quick_plot, not this program.
       !**********************************************************
       !place = 1-(1.0/iy_tot*iy)
       !  call qp_draw_text (titl, 1.0_rp, place, "%PAGE", "LT", 12.0_rp)

       call qp_draw_graph (xcoord, ycoord(:),title = titl)
       call qp_restore_state

       IF (iy == 1) THEN
          graph_more = .false.
       ENDIF

       iy = iy-1

    ENDDO
    !Formats

    deallocate (b)
    deallocate (nt)
    if (allocated(ycoord)) deallocate (ycoord)
    if (allocated(xcoord)) deallocate (xcoord)
    if (allocated(lam_log)) deallocate (lam_log)
  end subroutine plot_it


  subroutine min_max_y(vector, miny, maxy, ycoord, arr_length)

    !
    !Finds min and max y values from a given vector to be plotted.
    !Also assigns the vector to ycoord.
    !Used in plotting graphs.
    !

    real(rp) :: vector(:), &        !Some vector passed into the subroutine
         miny, maxy, &              !Find the min and max values of the vector
         max_factor, min_factor, &  !Used to find order of magnitude to
         !round min, max properly
         ycoord(:)                  !Takes the values of vector
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

    miny = floor(miny*100/(10**min_factor))*(10**min_factor)/100
    maxy = ceiling(maxy*100/(10**max_factor))*(10**max_factor)/100

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
    !routine can be used to plot data from different matrices
    !can use digital display but also writes a postscript
    !that can be printed for a hard copy

    type(data_set) :: data(*)          !Data from file
    integer :: i, &                    !Counter
         xdiv, &                       !Divisions of x axis
         graph, &                      !Choice to graph
         istat,&                       !Status
         in4get1,reaget1 , &           !Get input
         arr_length, iset, icolumn, count, nset
    real(rp) xlength, &                !Length of x axis (# values)
         miny, maxy                    !Min and max y values
    real(rp), allocatable :: xcoord(:), &  !X coordinates to plot
         ycoord(:)                     !Y coordinates to plot
    real(rp) :: sf2                    !Factor
    character(80) title, titl          !Titles
    logical :: graph_more              !To graph more
    integer :: ix, iy, ix_tot, iy_tot  !Number of plots in X and y directions
    real(rp), allocatable :: b(:), &   !Number of BPMs (for plotting)
         phi(:), &                     !Phi (for x-value plotting)
         nt(:), &                      !Number of turns
         lam_log(:)                    !Log of lambda (for plotting)
    real(rp) :: xmin, xmax
    logical :: badFile, badColumn

    allocate (b(2*NUM_BPMS))
    allocate (phi(2*NUM_BPMS))
    allocate (nt(NUM_TURNS))
    nset = 2                 !Make nset a global variable

    do i=1,2*NUM_BPMS
       b(i) = i
    end do
    do i=1,NUM_TURNS
       nt(i) = i
    end do

    ix_tot = 1                           !Graphs are stacked vertically,
    ix = ix_tot                          !Change ix_tot to plot side-by-side
    graph_more = .true.

    Print *, "How many graphs do you want?"
    PRINT *, "Enter a number between 0 and 5"
    Read *, iy_tot

    If (iy_tot>5) then
       print *, "Changing number of graphs to 5 (maximum)."
       iy_tot = 5
    endif
    if (iy_tot<=0) then
       print *, "Number of graphs 0 or less--nothing is being plotted."
       graph_more = .false.
    endif

    iy = iy_tot                             !Graph position 1 is the bottom 

    DO WHILE (graph_more)
       IF (ALLOCATED(xcoord)) DEALLOCATE(xcoord)  
       IF (ALLOCATED(ycoord)) DEALLOCATE(ycoord)
       !Moved inside loop because of changed x coordinates for phi
       write (*,'(1x,a)') "For A mode Beta Ratio (X), enter 1", &
            "For A mode Beta Ratio (Y), enter 2", &
            "For B mode Beta Ratio (X), enter 3", &
            "For B mode Beta Ratio (Y), enter 4", &
            "For A mode Magnitude**2 (X), enter 5", &
            "For A mode Magnitude**2 (Y), enter 6", &
            "For B mode Magnitude**2 (X), enter 7", &
            "For B mode Magnitude**2 (Y), enter 8", &
            "For A mode Gam**2 Beta Ratio, enter 9", &
            "For B mode Gam**2 Beta Ratio, enter 10", & 
            "For A mode Phase Advance, enter 11", &
            "For B mode Phase Advance, enter 12", &
            "For A mode Gamma**2 Beta, enter 13", &
            "For B mode Gamma**2 Beta, enter 14", &
            "For Sqrt(Beta) Cbar (1,1), enter 15", &
            "For Sqrt(Beta) Cbar (1,2), enter 16", &
            "For Sqrt(Beta) Cbar (2,2), enter 17", &
            "For Inv Gamma Cbar (1,1), enter 18", &
            "For Inv Gamma Cbar (1,2), enter 19", &
            "For Inv Gamma Cbar (2,2), enter 20", & 
            "For Pi, enter 21", &
            "For Frequency Peak,enter 22", &
            "For Lambda, enter 23", &
            "For Spectrum, enter 24", &
            "For Tau, enter 25"
       write (*, "(a)") " Enter choice "        
       accept "(i)", graph

       if (graph <= 8 .and. graph >= 5) then
          Write (*, '(1x,a)') "Type Magnitude**2 Scale Factor"
          write (*, "(a)") " Factor = "        
          accept "(f)", sf2
       end if

       if (graph >=21 .and. graph<=23) then
          arr_length = 2*NUM_BPMS
          allocate(xcoord(arr_length))
          allocate(ycoord(arr_length))
          xlength = 2.0*NUM_BPMS
          xdiv = 2*NUM_BPMS
          xcoord = b
       else if (graph>23) then
          arr_length = NUM_TURNS
          allocate(xcoord(arr_length))
          allocate(ycoord(arr_length))
          xdiv = 16
          xcoord = nt
          xlength = NUM_TURNS
       else
          arr_length = NUM_BPMS
          allocate (xcoord(arr_length))
          allocate (ycoord(arr_length))
          xlength = NUM_BPMS
          xdiv = NUM_BPMS
          xcoord = b
       endif
       xmin = minval(xcoord)
       xmax = maxval(xcoord)

       badFile = .true.
       if (graph > 20) then
          do while (badFile)
             write (*, "(a)") "Which file should be used?"
             accept "(i)", iset
             
             if (graph ==21 .or. graph > 23) then
                badColumn = .true.
                do while (badColumn)
                   !Column is needed for choices 21, 24, and 25
                   write (*, "(a)") " Which column should be plotted? "        
                   accept "(i)", icolumn
                   if (icolumn > 0 .and. icolumn < 2*NUM_BPMS) then
                      badColumn = .false.
                   else
                      Print *, "Bad column--choose a number between 1 and", &
                           2*NUM_BPMS
                   endif
                enddo
             endif
             if (iset >0 .and. iset <= nset) then
                badFile = .false.
             else
                print *, "Bad input--choose 1 or 2."
             endif
          enddo
       endif


       call qp_set_box(ix, iy, ix_tot, iy_tot)

       !
       !Formats
       !
26     format('(X) Ratio - A Mode')
27     format('(Y) Ratio - A Mode')
28     format('(X) Ratio - B Mode')
29     format('(Y) Ratio - B Mode')
30     format('(X) Magnitude**2 - A Mode')
31     format('(Y) Magnitude**2 - A Mode')
32     format('(X) Magnitude**2 - B Mode')
33     format('(Y) Magnitude**2 - B Mode')
34     format('Gam**2 Beta Ratio - A Mode')
35     format('Gam**2 Beta Ratio - B Mode')
36     format('Phase Advance - A Mode')
37     format('Phase Advance - B Mode')
38     format('Gamma**2 Beta - A Mode')
39     format('Gamma**2 Beta - B Mode')
40     format('Sqrt(Beta-A) Cbar (1,1)')
41     format('Sqrt(Aver-Beta) Cbar (1,2)')
42     format('Sqrt(Beta-B) Cbar(2,2)')
43     format('Inv Gamma Cbar (1,1)')  
44     format('Inv Gamma Cbar (1,2)')
45     format('Inv Gamma Cbar (2,2)')

       if (graph == 1) then
          call min_max_y(data_struc%loc(:)%a%ratio(1), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,26)
       else if (graph == 2) then
          call min_max_y(data_struc%loc(:)%a%ratio(2), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,27)
       else if (graph == 3) then
          call min_max_y(data_struc%loc(:)%b%ratio(1), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,28)
       else if (graph == 4) then
          call min_max_y(data_struc%loc(:)%b%ratio(2), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,29)
       else if (graph == 5) then
          call min_max_y(data_struc%loc(:)%a%magnitude2(1), &
               miny, maxy, ycoord(:), arr_length)
          ycoord(:) = data_struc%loc(:)%a%magnitude2(1) / maxy *sf2
          maxy = sf2
          write(titl,30)
       else if (graph == 6) then
          call min_max_y(data_struc%loc(:)%a%magnitude2(2), &
               miny, maxy, ycoord(:), arr_length)
          ycoord(:) = data_struc%loc(:)%a%magnitude2(2) / maxy * sf2
          maxy = sf2
          write(titl,31)
       else if (graph == 7) then
          call min_max_y(data_struc%loc(:)%b%magnitude2(1), &
               miny, maxy, ycoord(:), arr_length)
          ycoord(:) = data_struc%loc(:)%b%magnitude2(1) / maxy * sf2
          maxy = sf2
          write(titl,32)
       else if (graph == 8) then
          call min_max_y(data_struc%loc(:)%b%magnitude2(2), &
               miny, maxy, ycoord(:), arr_length)
          ycoord(:) = data_struc%loc(:)%b%magnitude2(2) / maxy * sf2
          maxy = sf2
          write(titl,33)
       else if (graph == 9) then
          call min_max_y(data_struc%loc(:)%a%gam2_beta_ratio, &
               miny, maxy, ycoord(:), arr_length)
          write(titl,34)
       else if (graph == 10) then
          call min_max_y(data_struc%loc(:)%b%gam2_beta_ratio, &
               miny, maxy, ycoord(:), arr_length)
          write(titl,35)
       else if (graph == 11) then
          phi(:) = data_struc%loc(:)%a%phi
          call arrange_phi(phi(:), xcoord(:))
          call min_max_y(phi(:), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,36)
          xmax = maxval(xcoord)
          xmin = minval(xcoord)
       else if (graph == 12) then
          phi(:) = data_struc%loc(:)%b%phi
          call arrange_phi(phi(:), xcoord(:))
          call min_max_y(phi(:), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,37)
          xmax = maxval(xcoord)
          xmin = minval(xcoord)
       else if (graph == 13) then
          call min_max_y(data_struc%loc(:)%a%gam2_beta, miny, maxy, &
               ycoord(:), arr_length)
          write(titl,38)
       else if (graph == 14) then
          call min_max_y(data_struc%loc(:)%b%gam2_beta, miny, maxy, &
               ycoord(:), arr_length)
          write(titl,39)
       else if (graph == 15) then
          call min_max_y(data_struc%loc(:)%sqrt_beta_cbar(1,1),&
               miny, maxy, ycoord(:), arr_length)
          write(titl,40)
       else if (graph == 16) then
          call min_max_y(data_struc%loc(:)%sqrt_beta_cbar(1,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,41)
       else if (graph == 17) then
          call min_max_y(data_struc%loc(:)%sqrt_beta_cbar(2,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,42)
       else if (graph == 18) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(1,1), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,43)
       else if (graph == 19) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(1,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,44)
       else if (graph == 20) then
          call min_max_y(data_struc%loc(:)%inv_gamma_cbar(2,2), &
               miny, maxy, ycoord(:), arr_length)
          write(titl,45)
       else if (graph == 21) then
          call min_max_y(data(iset)%pi_mat(:,icolumn), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,22) iset, icolumn, data(iset)%shortName
22        format('PI - Iset =',i2,'  Column = ',i3, '  File = ', a)
       else if (graph == 22) then
          !nmax is an integer; all others are real(rp)
          !Can't use a generic statement because of the type of nmax??
          !call min_max_y(data%fr_peak(:))
          miny = minval(data(iset)%fr_peak(:))
          maxy = maxval(data(iset)%fr_peak(:))
          ycoord(:) = data(iset)%fr_peak(:)
          write(titl,23) iset, data(iset)%shortName     
23        format('FR_PEAK - Iset =',i2,'  File = ', a)
       else if (graph == 23) then
          !Logarithmic scale option on qplot is not available yet
          !(See documentation)
          !Plots lambda on a logarithmic scale
          IF (ALLOCATED(lam_log)) DEALLOCATE(lam_log)    
          allocate (lam_log(arr_length))
          do count=1, arr_length
             lam_log(count) = log10(data(iset)%lambda(count))
          enddo
!          call sort_l(lam_log, 2*NUM_BPMS)
          call min_max_y(lam_log(:), miny, maxy, ycoord(:), arr_length)   
          write(titl,24) iset,  data(iset)%shortName
24        format('Lambda (log10) - Iset =',i2,'  File = ', a)
       else if (graph == 24) then
          call min_max_y(data(iset)%spectrum(:,icolumn), miny, maxy,&
               ycoord(:), arr_length)
          write(titl,25) iset,icolumn, data(iset)%shortName
25        format('Spectrum - Iset =',i2,'  Column = ',i2, '  File = ', a)
       else if (graph == 25) then
          call min_max_y(data(iset)%tau_mat(:,icolumn), miny, maxy, ycoord(:),&
               arr_length)
          write(titl,21) iset, icolumn, data(iset)%shortName
21        format('TAU - Iset =',i2,'  Column = ',i3, '  File = ', a)
       end if

       !Draws graph
       call qp_save_state(.true.)
       call qp_set_box (ix, iy, ix_tot, iy_tot)
       !call qp_set_graph_attrib (draw_grid = .true.)
       call qp_set_symbol_attrib (star5_filled$, height = 5.0_rp)
!       call qp_set_axis ("X", 0.0_rp, xlength, div = xdiv)
       !call qp_set_axis ("Y", miny, maxy, places = 4, draw_numbers = .true., &
       !div = 10)
       call qp_set_axis ("X", 0.0_rp, xlength, div = xdiv)
!       call qp_calc_and_set_axis("X", xmin, xmax, 4, 6, "GENERAL")
       call qp_calc_and_set_axis("Y", miny, maxy, 4, 6, "GENERAL")
       call qp_draw_graph (xcoord, ycoord(:),title = titl)
       call qp_restore_state

       IF (iy <= 1) THEN
          graph_more = .false.
       ENDIF

       iy = iy-1

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

end module mia_plot
