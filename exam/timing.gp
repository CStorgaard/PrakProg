# timings.gp
set terminal svg size 600,400 enhanced
set output 'timing.svg'

set title "Berrut Evaluation Time vs. n"
set xlabel "n (number of nodes)"
set ylabel "avg eval time (ns)"
set key left top

# log–log to check O(n) linear slope
set logscale x 2
set logscale y 10

plot \
  'Out_times.txt' using 1:2 with linespoints lw 2 pt 4 title 'avg time'
