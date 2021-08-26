wall_generator
 
                          __   ___       ___  __       ___  __   __  
 |  |  /\  |    |        / _` |__  |\ | |__  |__)  /\   |  /  \ |__) 
 |/\| /~~\ |___ |___ ___ \__> |___ | \| |___ |  \ /~~\  |  \__/ |  \ 
 
 Please provide a lattice file. Format:
 wall_generator <lattice>
   Example: wall_generator lat.bmad
 wall_generator <lattice> <n_angles> <ds>
   Example: wall_generator lat.bmad 8


This will produce a file like:

Wall for:lat.bmad                                                                                                                                                                                                
                 x          normal_x                 y          normal_y                 z          normal_z            ix_ele       angle_index                 s
                 m                 1                 m                 1                 m                 1                 1                 1                 m
  2.0000000000E-01  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00                 1                 0  0.0000000000E+00
 -2.0000000000E-01  0.0000000000E+00  2.4492935983E-17  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00                 1                 1  0.0000000000E+00
  2.0000000000E-01  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00  1.0000000000E+00  0.0000000000E+00                 1                 0  1.0000000000E+00
 -2.0000000000E-01  0.0000000000E+00  2.4492935983E-17  0.0000000000E+00  1.0000000000E+00  0.0000000000E+00                 1                 1  1.0000000000E+00


which can be read by another program to plot. See parse_wall_generator.ipynb for an example of how to do this. 

