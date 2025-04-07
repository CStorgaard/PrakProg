set terminal svg enhanced size 800,600
set output "Convergence_acc.svg"

set title "Convergence study: E₀ vs acc"
set xlabel "acc"
set ylabel "E₀"
set logscale x
set grid
set xrange [8e-8:0.03] reverse

plot "< grep '^acc' ConvergenceStudy.txt" using \
     (real(substr(strcol(3), 1, strlen(strcol(3))-1))):(real(strcol(6))) \
     with linespoints title "E₀ vs acc"
