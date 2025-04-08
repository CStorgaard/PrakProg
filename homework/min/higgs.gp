set terminal svg size 800,600
set output "Higgs.svg"

set title "Higgs data and fit"
set xlabel "Energy [GeV]"
set ylabel "Signal"
set key top right
set datafile commentschars "#"

# Plot the experimental data (block index 0) as dots with error bars.
plot "plot.dat" index 0 using 1:2:3 with yerrorbars pointtype 7 pointsize 1 title "Experimental Data", \
     "< sed -n '/^# Fitted/,/^$/p' plot.dat | sed '1d'" using 1:2 with lines lw 2 title "Fitted Curve"
