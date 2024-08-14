&tao_plot_page
  !plot_page%size = 500, 600
  !plot_page%size = 1000, 500
  plot_page%size = 1000, 650
  plot_page%draw_graph_title_suffix = F
  plot_page%text_height = 18
  place(1) = 'r11', 'floor'
  !place(2) = 'r22', 'multiorbit'
/

!!!------------  Multi floor plan ------------

&floor_plan_drawing
 include_default_shapes = T
 !ele_shape(N) = "<ele_id>" "<shape>" "<color>" "<size>" "<label>" <draw> <multi> <line_width>
 ele_shape(1) = "building_wall::overhead_view_floor_left" "solid_line" "orange" 
 ele_shape(2) = "building_wall::overhead_view_floor_right" "solid_line" "orange" 
 !ele_shape(3) = "building_wall::overhead_view_floor_right_wall" "solid_line" "black"
 ele_shape(3) = "sbend::*B*" "box" "blue" 0.25 "none"
 ele_shape(4) = "quadrupole::*" "xbox" "red" 0.1552 "none"
 ele_shape(5) = "sbend::*E" "xbox" "orange" 0.15 "none"
 ele_shape(6) = "marker::fq*" "box" "red" 0.15 "none" 
 ele_shape(7) = "marker::fe*" "box" "orange" 0.15 "none"
/

&tao_template_plot
 plot%name = "floor"
 plot%n_graph = 1
/

&tao_template_graph
 graph_index = 1
 graph%name = "overhead_view_floor"
 graph%type = "floor_plan"
 graph%box = 1, 1, 1, 1
 graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
 graph%ix_universe = -2 ! Draw all universes.
 !graph%y%major_div_nominal = 4
 !graph%x%major_div_nominal = 5
 graph%y%max = 2 !104
 graph%y%min = -2 !99
 !graph%y%max = 25
 graph%x%label_offset=-0.05
 graph%floor_plan%correct_distortion = F
 graph%floor_plan%size_is_absolute = T
 graph%floor_plan%view='zx'
 graph%floor_plan%orbit_color = "green"
 graph%floor_plan%orbit_width = 5
 graph%floor_plan%orbit_scale = 1
 graph%x%label = "Z (m)"
 graph%y%label = "X (m)"
 graph%title = "Floor Plan"
/


!!!------------  Multi beta a ------------

&tao_template_plot
  plot%name = 'multiba'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'ba'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gb\dA\u (m)'
  graph%y%label_offset=.2
  graph%x%label_offset=-0.05
  graph%n_curve = 6
  graph%title = "Horizontal Beta"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'beta.a'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'beta.a'
  curve(1)%y_axis_scale_factor = 1

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'beta.a'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'beta.a'
  curve(2)%y_axis_scale_factor = 1

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'beta.a'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'beta.a'
  curve(3)%y_axis_scale_factor = 1

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'beta.a'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'beta.a'
  curve(4)%y_axis_scale_factor = 1

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'beta.a'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'beta.a'
  curve(5)%y_axis_scale_factor = 1

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'beta.a'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'beta.a'
  curve(6)%y_axis_scale_factor = 1
/

!!!------------  Multi beta b ------------

&tao_template_plot
  plot%name = 'multibb'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'bb'

 !graph%title = ''
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gb\dB\u (m)'
  graph%y%label_offset=.2
  graph%x%label_offset=-0.05
  graph%n_curve = 6
  graph%title = "Vertical Beta"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'beta.b'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'beta.b'
  curve(1)%y_axis_scale_factor = 1

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'beta.b'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'beta.b'
  curve(2)%y_axis_scale_factor = 1

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'beta.b'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'beta.b'
  curve(3)%y_axis_scale_factor = 1

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'beta.b'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'beta.b'
  curve(4)%y_axis_scale_factor = 1

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'beta.b'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'beta.b'
  curve(5)%y_axis_scale_factor = 1

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'beta.b'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'beta.b'
  curve(6)%y_axis_scale_factor = 1
/

!!!------------  Multi Alpha a ------------

&tao_template_plot
  plot%name = 'multiaa'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'aa'

  graph_index = 1

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\ga\dA\u (m)'
  graph%y%label_offset=.2
  graph%x%label_offset=-0.05
  graph%n_curve = 6
  graph%title = "Horizontal Alpha"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'alpha.a'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'alpha.a'
  curve(1)%y_axis_scale_factor = 1

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'alpha.a'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'alpha.a'
  curve(2)%y_axis_scale_factor = 1

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'alpha.a'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'alpha.a'
  curve(3)%y_axis_scale_factor = 1

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'alpha.a'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'alpha.a'
  curve(4)%y_axis_scale_factor = 1

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'alpha.a'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'alpha.a'
  curve(5)%y_axis_scale_factor = 1

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'alpha.a'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'alpha.a'
  curve(6)%y_axis_scale_factor = 1
/

!!!------------  Multi alpha b ------------

&tao_template_plot
  plot%name = 'multiab'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'ab'

 !graph%title = ''
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\ga\dB\u (m)'
  graph%y%label_offset=.2
  graph%x%label_offset=-0.05
  graph%n_curve = 6
  graph%title = "Vertical Alpha"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'alpha.b'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'alpha.b'
  curve(1)%y_axis_scale_factor = 1

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'alpha.b'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'alpha.b'
  curve(2)%y_axis_scale_factor = 1

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'alpha.b'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'alpha.b'
  curve(3)%y_axis_scale_factor = 1

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'alpha.b'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'alpha.b'
  curve(4)%y_axis_scale_factor = 1

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'alpha.b'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'alpha.b'
  curve(5)%y_axis_scale_factor = 1

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'alpha.b'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'alpha.b'
  curve(6)%y_axis_scale_factor = 1
/

!!!------------  Multi orbit ------------

&tao_template_plot
  plot%name = 'multiorbitx'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'orbitx'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'x (mm)'
  graph%y%label_offset=.2
  graph%n_curve = 6

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'orbit.x'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'orbit.x'
  curve(1)%y_axis_scale_factor = 1000

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'orbit.x'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'orbit.x'
  curve(2)%y_axis_scale_factor = 1000

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'orbit.x'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'orbit.x'
  curve(3)%y_axis_scale_factor = 1000

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'orbit.x'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'orbit.x'
  curve(4)%y_axis_scale_factor = 1000

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'orbit.x'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'orbit.x'
  curve(5)%y_axis_scale_factor = 1000

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'orbit.x'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'orbit.x'
  curve(6)%y_axis_scale_factor = 1000
/

&tao_template_plot
  plot%name = 'multiorbity'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'orbity'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'y (mm)'
  graph%y%label_offset=.2
  graph%n_curve = 6

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'orbit.y'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'orbit.y'
  curve(1)%y_axis_scale_factor = 1000

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'orbit.y'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'orbit.y'
  curve(2)%y_axis_scale_factor = 1000

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'orbit.y'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'orbit.y'
  curve(3)%y_axis_scale_factor = 1000

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'orbit.y'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'orbit.y'
  curve(4)%y_axis_scale_factor = 1000

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'orbit.y'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'orbit.y'
  curve(5)%y_axis_scale_factor = 1000

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'orbit.y'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'orbit.y'
  curve(6)%y_axis_scale_factor = 1000
/


!!!------------  Multi eta ------------

&tao_template_plot
  plot%name = 'multietax'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etax'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gy\dx\u (cm)'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Horizontal Dispersion"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'eta.x'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'eta.x'
  curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'eta.x'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'eta.x'
  curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'eta.x'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'eta.x'
  curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'eta.x'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'eta.x'
  curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'eta.x'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'eta.x'
  curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'eta.x'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'eta.x'
  curve(6)%y_axis_scale_factor = 100
/

&tao_template_plot
  plot%name = 'multietay'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etay'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gy\dx\u (cm)'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Vertical Dispersion"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'eta.y'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'eta.y'
  curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'eta.y'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'eta.y'
  curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'eta.y'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'eta.y'
  curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'eta.y'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'eta.y'
  curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'eta.y'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'eta.y'
  curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'eta.y'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'eta.y'
  curve(6)%y_axis_scale_factor = 100
/

!!!------------  Multi etap ------------

&tao_template_plot
  plot%name = 'multietapx'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etax'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'd\gy\dx\u'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Horizontal Dispersion Prime"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'etap.x'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'etap.x'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'etap.x'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'etap.x'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'etap.x'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'etap.x'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'etap.x'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'etap.x'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'etap.x'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'etap.x'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'etap.x'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'etap.x'
  !curve(6)%y_axis_scale_factor = 100
/

&tao_template_plot
  plot%name = 'multietay'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etay'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'd\gy\dx\u'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Vertical Dispersion Prime"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'etap.y'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'etap.y'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'etap.y'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'etap.y'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'etap.y'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'etap.y'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'etap.y'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'etap.y'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'etap.y'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'etap.y'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'etap.y'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'etap.y'
  !curve(6)%y_axis_scale_factor = 100
/

!!!------------  Multi time ------------

&tao_template_plot
  plot%name = 'multitime'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'time'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'Time (s)'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Time"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'time'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'time'
  curve(1)%y_axis_scale_factor = 0.00000001

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'time'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'time'
  curve(2)%y_axis_scale_factor = 0.00000001

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'time'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'time'
  curve(3)%y_axis_scale_factor = 0.00000001

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'time'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'time'
  curve(4)%y_axis_scale_factor = 0.00000001

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'time'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'time'
  curve(5)%y_axis_scale_factor = 0.00000001

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'time'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'time'
  curve(6)%y_axis_scale_factor = 0.00000001
/

!!!------------  Multi phase ------------

&tao_template_plot
  plot%name = 'multiphasex'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etax'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gf\dx\u (rad)'
  graph%y%label_offset=.2
  graph%n_curve = 6

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'phase.a'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'phase.a'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'phase.a'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'phase.a'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'phase.a'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'phase.a'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'phase.a'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'phase.a'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'phase.a'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'phase.a'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'phase.a'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'phase.a'
  !curve(6)%y_axis_scale_factor = 100
/

&tao_template_plot
  plot%name = 'multiphasey'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'etay'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = '\gf\dy\u (rad)'
  graph%y%label_offset=.2
  graph%n_curve = 6

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'phase.b'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'phase.b'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'phase.b'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'phase.b'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'phase.b'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'phase.b'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'phase.b'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'phase.b'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'phase.b'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'phase.b'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'phase.b'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'phase.b'
  !curve(6)%y_axis_scale_factor = 100
/

!!!------------  Multi R56 ------------

&tao_template_plot
  plot%name = 'multir56'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'r56'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'R\d56\u (m)'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "R56"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'r.56'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'r.56'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'r.56'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'r.56'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'r.56'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'r.56'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'r.56'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'r.56'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'r.56'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'r.56'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'r.56'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'r.56'
  !curve(6)%y_axis_scale_factor = 100
/

!!!------------  Multi Momentum Compaction ------------

&tao_template_plot
  plot%name = 'multiMomentumCompaction'
  plot%x_axis_type = 's'
  plot%x%label = 's (m)'
  !graph%ix_universe = -2
  plot%n_graph = 1
/

&tao_template_graph
  graph%name = 'MomComp'

 !graph%title = 'Lattice orbit'
 !graph%x%draw_numbers = .false.
 !graph%x%draw_label = .false.
  graph_index = 1
 
 !graph%y%min = 0
 !graph%y%max =  15

  graph%margin = 0.06, 0.06, 0.05, 0.1, "%BOX"
  graph%y%label = 'Momentum Compaction (m)'
  graph%y%label_offset=.2
  graph%n_curve = 6
  graph%title = "Momentum Compaction"

  curve(1)%ix_universe = 1
  curve(1)%data_source = 'lattice'
  curve(1)%data_type   = 'momentum_compaction'
  curve(1)%y_axis_scale_factor = 1
  curve(1)%line%color = 1
  curve(1)%line%width=2
  curve(1)%draw_symbols=.false.
  curve(1)%legend_text = 'momentum_compaction'
  !curve(1)%y_axis_scale_factor = 100

  curve(2)%ix_universe = 2
  curve(2)%data_source = 'lattice'
  curve(2)%data_type   = 'momentum_compaction'
  curve(2)%y_axis_scale_factor = 1
  curve(2)%line%color = 2
  curve(2)%line%width=2
  curve(2)%draw_symbols=.false.
  curve(2)%legend_text = 'momentum_compaction'
  !curve(2)%y_axis_scale_factor = 100

  curve(3)%ix_universe = 3
  curve(3)%data_source = 'lattice'
  curve(3)%data_type   = 'momentum_compaction'
  curve(3)%y_axis_scale_factor = 1
  curve(3)%line%color = 3
  curve(3)%line%width=2
  curve(3)%draw_symbols=.false.
  curve(3)%legend_text = 'momentum_compaction'
  !curve(3)%y_axis_scale_factor = 100

  curve(4)%ix_universe = 4
  curve(4)%data_source = 'lattice'
  curve(4)%data_type   = 'momentum_compaction'
  curve(4)%y_axis_scale_factor = 1
  curve(4)%line%color = 4
  curve(4)%line%width=2
  curve(4)%draw_symbols=.false.
  curve(4)%legend_text = 'momentum_compaction'
  !curve(4)%y_axis_scale_factor = 100

  curve(5)%ix_universe = 5
  curve(5)%data_source = 'lattice'
  curve(5)%data_type   = 'momentum_compaction'
  curve(5)%y_axis_scale_factor = 1
  curve(5)%line%color = 5
  curve(5)%line%width=2
  curve(5)%draw_symbols=.false.
  curve(5)%legend_text = 'momentum_compaction'
  !curve(5)%y_axis_scale_factor = 100

  curve(6)%ix_universe = 6
  curve(6)%data_source = 'lattice'
  curve(6)%data_type   = 'momentum_compaction'
  curve(6)%y_axis_scale_factor = 1
  curve(6)%line%color = 6
  curve(6)%line%width=2
  curve(6)%draw_symbols=.false.
  curve(6)%legend_text = 'momentum_compaction'
  !curve(6)%y_axis_scale_factor = 100
/
