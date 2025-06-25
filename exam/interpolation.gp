# interpolation.gp

set terminal svg size 600,400 dynamic enhanced
set output 'interpolation.svg'

set title "Berrut Interpolation: first (r1) vs second (r2)"
set xlabel "x"
set ylabel "y"
set key outside right

# semi-transparent red for r1, semi-transparent blue for r2
set style line 1 lc rgb '#FF0000' lw 2        # r1: red
set style line 2 lc rgb '#0000FF' lw 2 dt 2   # r2: blue, dashed

plot \
  'Out.txt' index 0 using 1:2 with points pt 7 ps 1.2 title 'data', \
  'Out.txt' index 1 using 1:2 with lines ls 1 title 'r1 (constant-exact)', \
  'Out.txt' index 1 using 1:3 with lines ls 2 title 'r2 (linear-exact)'
