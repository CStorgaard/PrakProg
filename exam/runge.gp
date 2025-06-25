# runge.gp
set terminal svg size 600,400 enhanced
set output 'runge.svg'

set title "Runge test: Berrut r₁ vs Berrut r₂ vs truth"
set xlabel "x"
set ylabel "y"
set key outside right

# styles
set style line 1 lc rgb 'red'    lw 2       # r1
set style line 2 lc rgb 'blue'   lw 2 dt 2  # r2
set style line 3 lc rgb 'black'  lw 1       # true f

plot \
  'Out_runge.txt' index 0 using 1:2 with points pt 7 ps 1.2 title 'nodes', \
  'Out_runge.txt' index 1 using 1:2 with lines ls 1 title 'r1 (const-exact)', \
  'Out_runge.txt' index 1 using 1:3 with lines ls 2 title 'r2 (linear-exact)', \
  'Out_runge.txt' index 1 using 1:4 with lines ls 3 title 'f(x)'
