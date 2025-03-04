# plot.gp
set terminal svg size 800,600
set output "decayFit.svg"

set title "Radioactive Decay Fit"
set xlabel "Time (days)"
set ylabel "Activity"
set key top right

# Add labels: adjust the text and position as needed.
set label 1 "Calculated half-life: 4.123 days" at graph 0.05, 0.90
set label 2 "True half-life: 3.6319 days" at graph 0.05, 0.85

set logscale y

plot "dataPoints.txt" using 1:2:3 with yerrorbars title "Data", \
     "fitPoints.txt" using 1:2 with lines title "Best Fit", \
     "fitEnvelope.txt" using 1:3 with lines title "1-sigma upper", \
     "fitEnvelope.txt" using 1:4 with lines title "1-sigma lower"