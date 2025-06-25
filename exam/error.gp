# error.gp
set terminal svg size 600,400 enhanced
set output 'error.svg'

set title "Interpolation Absolute Error of f(x) = x + 0.2sin(5x)"
set xlabel "x"
set ylabel "error=|r(x) - f(x)|"
set key outside right

plot \
  'Out_error.txt' using 1:5 with lines lw 2 lc rgb 'red'   title 'err1 (const-exact)', \
  'Out_error.txt' using 1:6 with lines lw 2 lc rgb 'blue'  title 'err2 (linear-exact)'
