# plot.gp
set terminal svg size 800,600
set output "decayFit.svg"

set title "Error function"
set xlabel "x"
set ylabel "y"
set key top right

# Add labels: adjust the text and position as needed.
plot "erf.txt" using 1:2:3 with yerrorbars title "Error function", \