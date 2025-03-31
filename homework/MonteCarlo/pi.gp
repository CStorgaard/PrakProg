# pi.gp
set terminal svg size 800,600
set output "pi.svg"

set title "Pi estimated and actual error"
set xlabel "N"
set ylabel "Error"
set key top right
set logscale x
set logscale y
set grid

plot "pi_results.txt" using 1:3 with lines title "Estimated error", \
     "pi_results.txt" using 1:4 with lines title "Actual error"