# Use this file with the command:
#       load 'surface_plot.pg'
# Note: If your data file is not named "surface.dat" please edit the last line below.

set print "-"
print "Note: Plotted is -Z (not Z) so the view is from above."

set dgrid3d splines
set grid
set view 50
set style data lines
set contour base
set hidden3d
set xlabel "X"
set ylabel "Y"
set zlabel "-Z" offset -3
splot 'surface.dat' u 3:4:(-$5)
