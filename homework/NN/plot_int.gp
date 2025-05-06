set terminal svg size 800,600
set output 'int.svg'
set title "Net Antiderivative A(x)"
set xlabel 'x'
set ylabel 'A(x)'
plot \
  'Out.txt' using 1:6 title "A(x)" with lines
