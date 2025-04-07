set terminal svg enhanced size 800,600
set output "Convergence_max.svg"

set title "Convergence of energy with varying r_{max} (r_{min} fixed) "
set xlabel "r_{max}"
set ylabel "E0"
set key top right
set grid

plot "< sed -n '2,8p' ConvergenceStudy.txt" using \
  (real(substr(strcol(3), 1, strlen(strcol(3))-1))):(real(strcol(6))) with linespoints title "E0 vs r_{max}"
