module quick_plot_struct

use precision_def

integer, parameter :: White$ = 0, Black$ = 1, Red$ = 2, Green$ = 3
integer, parameter :: Blue$ = 4, Cyan$ = 5, Magenta$ = 6, Yellow$ = 7
integer, parameter :: Orange$ = 8, Yellow_Green$ = 9, Light_Green$ = 10
integer, parameter :: Navy_Blue$ = 11, Purple$ = 12, Redish_Purple$ = 13
integer, parameter :: Dark_Grey$ = 14, Light_Grey$ = 15

character(16), parameter :: qp_color_name(0:15) = (/ 'White        ', &
    'Black        ', 'Red          ', 'Green        ', 'Blue         ', &
    'Cyan         ', 'Magenta      ', 'Yellow       ', 'Orange       ', &
    'Yellow_Green ', 'Light_Green  ', 'Navy_Blue    ', 'Purple       ', &
    'Redish_Purple', 'Dark_Grey    ', 'Light_Grey   ' /)

integer, parameter :: solid$ = 1, dashed$ = 2, dash_dot$ = 3
integer, parameter :: dotted$ = 4, dash_dot3$ = 5

character(16) :: qp_line_style_name(5) = (/ 'solid    ', &
    'dashed   ', 'dash_dot ', 'dotted   ', 'dash_dot3' /)

integer, parameter :: solid_fill$ = 1, no_fill$ = 2
integer, parameter :: hatched$ = 3, cross_hatched$ = 4

character(16) :: qp_fill_name(4) = (/ 'solid_fill   ', 'no_fill      ', &
                                      'hatched      ', 'cross_hatched' /)

integer, parameter :: square$ = 0, dot$ = 1, plus$ = 2, times$ = 3
integer, parameter :: circle$ = 4, x_symbol$ = 5, triangle$ = 7
integer, parameter :: circle_plus$ = 8, circle_dot$ = 9
integer, parameter :: square_concave$ = 10, diamond$ = 11
integer, parameter :: star5$ = 12, triangle_filled$ = 13, red_cross$ = 14
integer, parameter :: star_of_david$ = 15, square_filled$ = 16
integer, parameter :: circle_filled$ = 17, star5_filled$ = 18

character(16) :: qp_symbol_type_name(0:18) = (/ 'square         ', &
    'dot            ', 'plus           ', 'times          ', 'circle         ', &
    'x_symbol       ', '---------------', 'triangle       ', 'circle_plus    ', &
    'circle_dot     ', 'square_concave ', 'diamond        ', 'star5          ', &
    'triangle_filled', 'red_cross      ', 'star_of_david  ', 'square_filled  ', &
    'circle_filled  ', 'star5_filled   ' /)

integer, parameter :: dflt_draw$ = 1, dflt_set$ = 2

end module
