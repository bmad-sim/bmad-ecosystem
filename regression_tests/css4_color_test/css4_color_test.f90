program css4_color_test
!
! Test program to verify CSS4 color name lookup and RGB values.
! Outputs each color name, its resolved index, and the RGB values that
! would be sent to the plotting backend.
!
! Usage: Run this program and compare output against css4_reference_swatch_expected.txt
!

use quick_plot_struct

implicit none

integer :: ix_enum, i
logical :: error
character(24) :: test_names(138)
character(24) :: round_trip_name

! CSS4 color names to test (case variations to exercise case-insensitivity)
test_names = [character(24) :: &
  'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', &
  'blanchedalmond', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', &
  'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 'crimson', 'darkblue', 'darkcyan', &
  'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', &
  'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', &
  'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', &
  'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', &
  'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', &
  'goldenrod', 'gray', 'greenyellow', 'grey', 'honeydew', 'hotpink', 'indianred', &
  'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon', &
  'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgray', &
  'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', &
  'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', &
  'limegreen', 'linen', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', &
  'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', &
  'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', &
  'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orangered', &
  'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', &
  'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'rebeccapurple', 'rosybrown', &
  'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', &
  'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', &
  'steelblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', &
  'whitesmoke', 'yellowgreen' ]

! Header
write(*, '(a)') '# CSS4 Color Test Results'
write(*, '(a)') '# index  name                     R    G    B  round_trip_name'

! Test original colors first
write(*, '(a)') '# --- Original quick_plot colors ---'
do i = 0, 15
  round_trip_name = qp_enum_to_string(i, 'color')
  write(*, '(i5, 2x, a24, 3i5, 2x, a24)') i, trim(qp_color_name(i)), &
    qp_color_red(i), qp_color_green(i), qp_color_blue(i), trim(round_trip_name)
enddo

! Test CSS4 colors
write(*, '(a)') '# --- CSS4 colors ---'
do i = 1, 138
  ix_enum = qp_string_to_enum(test_names(i), 'color', error=error)
  round_trip_name = qp_enum_to_string(ix_enum, 'color')

  if (error) then
    write(*, '(a, 2x, a24, a)') 'ERROR', trim(test_names(i)), ' - name not found!'
  else if (ix_enum < 0 .or. ix_enum > qp_max_color_index$) then
    write(*, '(a, 2x, a24, a, i8)') 'ERROR', trim(test_names(i)), ' - bad index:', ix_enum
  else
    write(*, '(i5, 2x, a24, 3i5, 2x, a24)') ix_enum, trim(test_names(i)), &
      qp_color_red(ix_enum), qp_color_green(ix_enum), qp_color_blue(ix_enum), &
      trim(round_trip_name)
  endif
enddo

! Test that original color names still resolve correctly
write(*, '(a)') '# --- Verify original names still work ---'
ix_enum = qp_string_to_enum('red', 'color', error=error)
write(*, '(a, i5, l3)') 'red -> ', ix_enum, .not. error
ix_enum = qp_string_to_enum('Salmon', 'color', error=error)
write(*, '(a, i5, l3)') 'Salmon -> ', ix_enum, .not. error
ix_enum = qp_string_to_enum('DodgerBlue', 'color', error=error)
write(*, '(a, i5, l3)') 'DodgerBlue -> ', ix_enum, .not. error

! Test continuous color still works
ix_enum = qp_string_to_enum('Z0.5000000000000', 'color', error=error)
write(*, '(a, i12, l3)') 'Z0.5 -> ', ix_enum, .not. error
write(*, '(a, l3)') 'Z0.5 >= continuous_start: ', ix_enum >= qp_continuous_color_start$

end program css4_color_test
