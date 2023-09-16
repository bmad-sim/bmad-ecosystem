

#!/usr/bin/env gnuplot

emins = "30 600 1400 1600"

set xlabel 'Energy [eV]'
set ylabel 'Reflectivity' offset 1,0
set key right out reverse Left title 'Angle' samplen 0.5 box
unset colorbox
set xtics rotate by -45

outname = 'Al_0nmAl.pdf'
set out outname
set term pdf size 20,5 fontscale 0.7
#set term png size 1500,400 fontscale 0.7
set multiplot layout 1,words(emins)

do for [emin in emins] {

fs = system('ls -1 Al_0nmAl/'.emin.'*raw.dat')

# pull out the angle from the filename
ang(f) = sprintf("%g",1*f[strstrt(f,"_raw")-4:strstrt(f,"_raw")-1])

plot for [i=1:words(fs)] word(fs,i) u 1:2:(log(1*ang(word(fs,i)))) w lp palette pt 7 ps 0.5 t ang(word(fs,i))

}
unset multi
set out
print 'Wrote to '.outname
