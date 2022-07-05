f(x) = exp(-x**2/2) / sqrt(2*pi)

a1 = 0.1290322580645162
s1 = 1.7
b1 = 0.3629032258064516
c1 = 1.2
cur1(x) = a1 * exp(-x**2/(2*s1**2)) / (sqrt(2*pi) * s1) + 3*b1 / (pi * (1+(x/c1)**6))

a2 = 0.3149433686518084
s2 = 0.5
b2 = 0.2088587290695706
c2 = 1.64
cur2(x) = a2 * exp(-x**2/(2*s2**2)) / (sqrt(2*pi) * s2) + 3*b2 / (pi * (1+(x/c2)**6))

set xrange [-4:4]
set yrange [0:0.5]

set border linewidth 2
set tics font ", 12"
set xlabel font ", 14"
set ylabel font ", 14"
set size ratio 0.5
set xlabel "X"
set ylabel 'Probability Density' 

plot f(x) lw 3 title "Gaussian", cur1(x) dashtype 2 lw 3 title "", cur2(x) dashtype 2 lw 3 title ""
