set terminal svg size 800,600
set output "quadraticSplinePlot.svg"

set title "Quadratic Spline, Integral, and Data Points"
set xlabel "x"
set ylabel "y"

plot "dataPoints.txt" using 1:2 with points pt 7 ps 1.5 title "Data Points", \
     "quadraticSplineResults.txt" using 1:2 with lines title "Quadratic Spline S(x)", \
     "quadraticSplineResults.txt" using 1:4 with lines title "Quadratic Spline Integral"
