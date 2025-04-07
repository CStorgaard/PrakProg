set terminal svg enhanced size 800,600
set output "Convergence_rmin.svg"

set title "Convergence study: (E₀ + 0.5) vs r_{min}"
set xlabel "r_{min}"
set ylabel "E₀ + 0.5"
set logscale x
# Reverse the x-axis; assuming rmin spans from 1e-5 to 0.01:
set xrange [8e-6:0.03] reverse
set grid

plot "< grep '^rmin' ConvergenceStudy.txt" using \
     (real(trim(substr(strcol(3), 1, strlen(strcol(3))-1)))):(real(strcol(6)) + 0.5) \
     with linespoints title "E₀ + 0.5 vs r_{min}"
