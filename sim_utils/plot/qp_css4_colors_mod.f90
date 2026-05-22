module qp_css4_colors_mod

implicit none

! CSS4 (CSS Color Level 4) named colors for quick_plot.
!
! Index scheme:
!   0-15:    Original quick_plot named colors (RGB values match CSS4 spec)
!   16:      Transparent (special)
!   17-154:  CSS4 named colors (138 colors, excluding 10 that duplicate existing names)
!   155+:    Continuous color mapping
!
! Colors are ordered by visual distinctiveness (farthest-point-first algorithm)
! so that devices with limited color indices get the most useful colors first.
!
! Excluded CSS4 names (same name and RGB as original quick_plot colors):
!   black, blue, cyan, green, magenta, orange, purple, red, white, yellow
!
! Historical note: green, orange, and purple were updated to match CSS4 values.
!   Old green  (0, 255, 0)   -> now available as "lime"  (index 17)
!   Old orange (255, 127, 0) -> closest: "darkorange" (255, 140, 0)
!   Old purple (127, 0, 255) -> closest: "blueviolet" (138, 43, 226)

integer, parameter :: qp_n_css4_colors = 138
integer, parameter :: qp_max_color_index$ = 154
integer, parameter :: qp_continuous_color_start$ = 155
integer, parameter :: qp_custom_color$ = -2   ! Sentinel: use qp_custom_rgb instead of palette lookup

! Storage for arbitrary RGB color specified via hex (#RRGGBB) or RGB(r,g,b) string.
! Set by qp_string_to_enum when it returns qp_custom_color$.
integer, save :: qp_custom_rgb(3) = [0, 0, 0]  ! R, G, B values (0-255)

! RGB values for all discrete colors (indices 0:154)
! Ordered: original 16 + transparent, then CSS4 by distinctiveness.
integer, parameter :: qp_color_red(0:154) = [ &
  255,   0, 255,   0,   0,   0, 255, 255, 255, 127, &  ! 0-9
    0,   0, 128, 255,  85, 170,   0,   0, 128,   0, &  ! 10-19
    0, 205, 106, 128, 238, 127, 240,  50, 148, 244, &  ! 20-29
   95, 154, 186, 255, 220, 144,  72, 255, 100,   0, &  ! 30-39
   46, 160, 255, 173, 102, 199, 189,  32, 184, 128, &  ! 40-49
  255, 173,  75, 178, 255, 205,   0,  60, 219,  34, &  ! 50-59
  153,   0, 147,  65, 220, 218, 102, 135,  85, 210, &  ! 60-69
  222, 107,  70, 143, 188, 250, 255, 245,  47,  25, &  ! 70-79
  221, 192, 218, 139, 105,  72,  30, 230, 123, 255, &  ! 80-89
  216, 255, 233, 255,   0, 119, 138,   0, 224, 255, &  ! 90-99
  255, 175, 176, 240, 255, 255, 255,  64, 165, 152, &  ! 100-109
  240, 240, 211, 139,   0, 238, 255, 135, 245, 255, &  ! 110-119
  112, 255, 210, 255, 250, 255, 245,   0, 139, 255, &  ! 120-129
  245, 255, 255, 176, 255, 248, 250, 255, 240, 250, &  ! 130-139
  253, 255, 124, 169,   0, 127, 169,  47, 105, 255, &  ! 140-149
  128, 211, 119, 112,   0 ]

integer, parameter :: qp_color_green(0:154) = [ &
  255,   0,   0, 128,   0, 255,   0, 255, 165, 255, &  ! 0-9
  255, 127,   0,   0,  85, 170,   0, 255,   0,   0, &  ! 10-19
  139,  92,  90, 128, 130, 255, 230, 205,   0, 164, &  ! 20-29
  158, 205,  85, 192,  20, 238, 209,  69, 149, 206, &  ! 30-39
  139,  82, 105, 216,  51,  21, 183, 178, 134, 128, &  ! 40-49
   99, 255,   0,  34, 250, 133,   0, 179, 112, 139, &  ! 50-59
   50, 191, 112, 105, 220, 165, 205, 206, 107, 105, &  ! 60-69
  184, 142, 130, 188, 143, 128, 215, 222,  79,  25, &  ! 70-79
  160, 192, 112,  69, 105,  61, 144, 230, 104, 228, &  ! 80-89
  191, 127, 150,  20, 100, 136,  43, 250, 255, 140, &  ! 90-99
  160, 238, 196, 255, 228, 255, 240, 224,  42, 251, &  ! 100-109
  128, 248, 211,   0, 128, 232, 255, 206, 245, 182, &  ! 110-119
  128, 239, 180, 218, 240, 222, 245,   0,   0, 228, &  ! 120-129
  255, 235, 245, 224, 248, 248, 250, 250, 255, 235, &  ! 130-139
  245, 250, 252, 169, 255, 255, 169,  79, 105,   0, &  ! 140-149
  128, 211, 136, 128, 255 ]

integer, parameter :: qp_color_blue(0:154) = [ &
  255,   0,   0,   0, 255, 255, 255,   0,   0,   0, &  ! 0-9
  127, 255, 128, 127,  85, 170,   0,   0,   0, 128, &  ! 10-19
  139,  92, 205,   0, 238, 212, 140,  50, 211,  96, &  ! 20-29
  160,  50, 211, 203,  60, 144, 204,   0, 237, 209, &  ! 30-39
   87,  45, 180, 230, 153, 133, 107, 170,  11, 128, &  ! 40-49
   71,  47, 130,  34, 205,  63, 205, 113, 147,  34, &  ! 50-59
  204, 255, 219, 225, 220,  32, 170, 250,  47,  30, &  ! 60-69
  135,  35, 180, 143, 143, 114,   0, 179,  79, 112, &  ! 70-79
  221, 192, 214,  19, 105, 139, 255, 250, 238, 225, &  ! 80-89
  216,  80, 122, 147,   0, 153, 226, 154, 255,   0, &  ! 90-99
  122, 238, 222, 240, 196, 224, 245, 208,  42, 152, &  ! 100-109
  128, 255, 211, 139, 128, 170, 240, 235, 220, 193, &  ! 110-119
  144, 213, 140, 185, 230, 173, 245, 139,   0, 181, &  ! 120-129
  250, 205, 238, 230, 220, 255, 210, 250, 255, 215, &  ! 130-139
  230, 240,   0, 169, 255,   0, 169,  79, 105, 255, &  ! 140-149
  128, 211, 153, 144, 127 ]

character(24), parameter :: qp_css4_color_name(17:154) = [ character(24) :: &
  'lime', &
  'maroon', &
  'navy', &
  'darkcyan', &
  'indianred', &  ! index 21
  'slateblue', &
  'olive', &
  'violet', &
  'aquamarine', &
  'khaki', &  ! index 26
  'limegreen', &
  'darkviolet', &
  'sandybrown', &
  'cadetblue', &
  'yellowgreen', &  ! index 31
  'mediumorchid', &
  'pink', &
  'crimson', &
  'lightgreen', &
  'mediumturquoise', &  ! index 36
  'orangered', &
  'cornflowerblue', &
  'darkturquoise', &
  'seagreen', &
  'sienna', &  ! index 41
  'hotpink', &
  'lightblue', &
  'rebeccapurple', &
  'mediumvioletred', &
  'darkkhaki', &  ! index 46
  'lightseagreen', &
  'darkgoldenrod', &
  'gray', &
  'tomato', &
  'greenyellow', &  ! index 51
  'indigo', &
  'firebrick', &
  'lemonchiffon', &
  'peru', &
  'mediumblue', &  ! index 56
  'mediumseagreen', &
  'palevioletred', &
  'forestgreen', &
  'darkorchid', &
  'deepskyblue', &  ! index 61
  'mediumpurple', &
  'royalblue', &
  'gainsboro', &
  'goldenrod', &
  'mediumaquamarine', &  ! index 66
  'lightskyblue', &
  'darkolivegreen', &
  'chocolate', &
  'burlywood', &
  'olivedrab', &  ! index 71
  'steelblue', &
  'darkseagreen', &
  'rosybrown', &
  'salmon', &
  'gold', &  ! index 76
  'wheat', &
  'darkslategray', &
  'midnightblue', &
  'plum', &
  'silver', &  ! index 81
  'orchid', &
  'saddlebrown', &
  'dimgray', &
  'darkslateblue', &
  'dodgerblue', &  ! index 86
  'lavender', &
  'mediumslateblue', &
  'mistyrose', &
  'thistle', &
  'coral', &  ! index 91
  'darksalmon', &
  'deeppink', &
  'darkgreen', &
  'lightslategray', &
  'blueviolet', &  ! index 96
  'mediumspringgreen', &
  'lightcyan', &
  'darkorange', &
  'lightsalmon', &
  'paleturquoise', &  ! index 101
  'lightsteelblue', &
  'honeydew', &
  'bisque', &
  'lightyellow', &
  'lavenderblush', &  ! index 106
  'turquoise', &
  'brown', &
  'palegreen', &
  'lightcoral', &
  'aliceblue', &  ! index 111
  'lightgray', &
  'darkmagenta', &
  'teal', &
  'palegoldenrod', &
  'ivory', &  ! index 116
  'skyblue', &
  'beige', &
  'lightpink', &
  'slategray', &
  'papayawhip', &  ! index 121
  'tan', &
  'peachpuff', &
  'linen', &
  'navajowhite', &
  'whitesmoke', &  ! index 126
  'darkblue', &
  'darkred', &
  'moccasin', &
  'mintcream', &
  'blanchedalmond', &  ! index 131
  'seashell', &
  'powderblue', &
  'cornsilk', &
  'ghostwhite', &
  'lightgoldenrodyellow', &  ! index 136
  'snow', &
  'azure', &
  'antiquewhite', &
  'oldlace', &
  'floralwhite', &  ! index 141
  'lawngreen', &
  'darkgray', &
  'aqua', &
  'chartreuse', &
  'darkgrey', &  ! index 146
  'darkslategrey', &
  'dimgrey', &
  'fuchsia', &
  'grey', &
  'lightgrey', &  ! index 151
  'lightslategrey', &
  'slategrey', &
  'springgreen' ]

! Integer constants for CSS4 color indices.
! Note: "tan$" is omitted because it conflicts with the tangent operator (tan$)
! defined in bmad_struct.f90. Use the string "tan" with qp_string_to_enum instead.

integer, parameter :: lime$ = 17, maroon$ = 18, navy$ = 19, darkcyan$ = 20
integer, parameter :: indianred$ = 21, slateblue$ = 22, olive$ = 23, violet$ = 24
integer, parameter :: aquamarine$ = 25, khaki$ = 26, limegreen$ = 27, darkviolet$ = 28
integer, parameter :: sandybrown$ = 29, cadetblue$ = 30, yellowgreen$ = 31, mediumorchid$ = 32
integer, parameter :: pink$ = 33, crimson$ = 34, lightgreen$ = 35, mediumturquoise$ = 36
integer, parameter :: orangered$ = 37, cornflowerblue$ = 38, darkturquoise$ = 39, seagreen$ = 40
integer, parameter :: sienna$ = 41, hotpink$ = 42, lightblue$ = 43, rebeccapurple$ = 44
integer, parameter :: mediumvioletred$ = 45, darkkhaki$ = 46, lightseagreen$ = 47
integer, parameter :: darkgoldenrod$ = 48, gray$ = 49, tomato$ = 50, greenyellow$ = 51
integer, parameter :: indigo$ = 52, firebrick$ = 53, lemonchiffon$ = 54, peru$ = 55
integer, parameter :: mediumblue$ = 56, mediumseagreen$ = 57, palevioletred$ = 58
integer, parameter :: forestgreen$ = 59, darkorchid$ = 60, deepskyblue$ = 61
integer, parameter :: mediumpurple$ = 62, royalblue$ = 63, gainsboro$ = 64, goldenrod$ = 65
integer, parameter :: mediumaquamarine$ = 66, lightskyblue$ = 67, darkolivegreen$ = 68
integer, parameter :: chocolate$ = 69, burlywood$ = 70, olivedrab$ = 71, steelblue$ = 72
integer, parameter :: darkseagreen$ = 73, rosybrown$ = 74, salmon$ = 75, gold$ = 76
integer, parameter :: wheat$ = 77, darkslategray$ = 78, midnightblue$ = 79, plum$ = 80
integer, parameter :: silver$ = 81, orchid$ = 82, saddlebrown$ = 83, dimgray$ = 84
integer, parameter :: darkslateblue$ = 85, dodgerblue$ = 86, lavender$ = 87
integer, parameter :: mediumslateblue$ = 88, mistyrose$ = 89, thistle$ = 90, coral$ = 91
integer, parameter :: darksalmon$ = 92, deeppink$ = 93, darkgreen$ = 94, lightslategray$ = 95
integer, parameter :: blueviolet$ = 96, mediumspringgreen$ = 97, lightcyan$ = 98
integer, parameter :: darkorange$ = 99, lightsalmon$ = 100, paleturquoise$ = 101
integer, parameter :: lightsteelblue$ = 102, honeydew$ = 103, bisque$ = 104
integer, parameter :: lightyellow$ = 105, lavenderblush$ = 106, turquoise$ = 107
integer, parameter :: brown$ = 108, palegreen$ = 109, lightcoral$ = 110, aliceblue$ = 111
integer, parameter :: lightgray$ = 112, darkmagenta$ = 113, teal$ = 114
integer, parameter :: palegoldenrod$ = 115, ivory$ = 116, skyblue$ = 117, beige$ = 118
integer, parameter :: lightpink$ = 119, slategray$ = 120, papayawhip$ = 121
! tan$ omitted (conflict with tangent operator)
integer, parameter :: peachpuff$ = 123, linen$ = 124, navajowhite$ = 125
integer, parameter :: whitesmoke$ = 126, darkblue$ = 127, darkred$ = 128, moccasin$ = 129
integer, parameter :: mintcream$ = 130, blanchedalmond$ = 131, seashell$ = 132
integer, parameter :: powderblue$ = 133, cornsilk$ = 134, ghostwhite$ = 135
integer, parameter :: lightgoldenrodyellow$ = 136, snow$ = 137, azure$ = 138
integer, parameter :: antiquewhite$ = 139, oldlace$ = 140, floralwhite$ = 141
integer, parameter :: lawngreen$ = 142, darkgray$ = 143, aqua$ = 144, chartreuse$ = 145
integer, parameter :: darkgrey$ = 146, darkslategrey$ = 147, dimgrey$ = 148, fuchsia$ = 149
integer, parameter :: grey$ = 150, lightgrey$ = 151, lightslategrey$ = 152
integer, parameter :: slategrey$ = 153, springgreen$ = 154

end module qp_css4_colors_mod
