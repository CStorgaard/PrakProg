set terminal svg size 800,600
set output "cubicSplinePlot.svg"

set title "Cubic Spline, Integral, and Data Points"
set xlabel "x"
set ylabel "y"

plot "dataPoints.txt" using 1:2 with points pt 7 ps 1.5 title "Data Points", \
     "cubicSplineResults.txt" using 1:2 with lines title "Cubic Spline S(x)", \
     "cubicSplineResults.txt" using 1:4 with lines title "Cubic Spline Integral"
