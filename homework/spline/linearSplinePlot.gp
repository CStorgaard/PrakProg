set terminal svg size 800,600
set output "linearSplinePlot.svg"

set title "Linear Spline, Integral, and Data Points"
set xlabel "x"
set ylabel "y"

plot "dataPoints.txt" using 1:2 with points pt 7 ps 1.5 title "Data Points", \
     "linearSplineResults.txt" using 1:2 with lines title "Linear Spline S(x)", \
     "linearSplineResults.txt" using 1:4 with lines title "Linear Spline Integral"
