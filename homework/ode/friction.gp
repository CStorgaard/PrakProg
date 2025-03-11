set terminal svg size 800,600 fname 'Helvetica'  # Use SVG output
set output "friction.svg"                             # Output file

set multiplot layout 1,2 title "Harmonic Oscillator"  # 1 row, 2 columns

# First subplot: y(x)
set xlabel "x"
set ylabel "y"
set title "y(x)"
set xrange [0:10]
set yrange [-1.5:1.5]
set size square
plot "friction.txt" using 1:2 with linespoints title "y(x)"


# Second subplot: Phase portrait: y' vs. y
set xlabel "y"
set ylabel "y'"
set title "Phase Portrait"
set xrange [-1.5:1.5]
set yrange [-1.8:1.5]
set size square
plot "friction.txt" using 2:3 with linespoints title "Phase Portrait"

unset multiplot
set output
