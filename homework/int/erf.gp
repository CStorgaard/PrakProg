# plot.gp
set terminal svg size 800,600
set output "erf.svg"

set title "Error function"
set xlabel "x"
set ylabel "erf(x)"
set key top right

# Add labels: adjust the text and position as needed.
plot "erf.txt" using 1:2:3 with yerrorbars title "Error function", \