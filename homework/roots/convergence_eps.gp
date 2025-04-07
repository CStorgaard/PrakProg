set terminal svg enhanced size 800,600
set output "Convergence_eps.svg"

set title "Convergence study: E₀ vs eps"
set xlabel "eps"
set ylabel "E₀"
set logscale x
set grid
set xrange [8e-8:0.03] reverse

# The eps lines in ConvergenceStudy.txt look like:
# eps = 0.0100000, E0 = -0.499978722027009
# Field 3 is "0.0100000," (with a trailing comma) and field 6 is the energy.
plot "< grep '^eps' ConvergenceStudy.txt" using \
     (real(substr(strcol(3), 1, strlen(strcol(3))-1))):(real(strcol(6))) \
     with linespoints title "E₀ vs eps"
