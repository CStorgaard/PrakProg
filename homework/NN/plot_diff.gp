set terminal svg size 800,600
set output 'diff.svg'
set title "Net' and Net''"
set xlabel 'x'
set ylabel 'derivatives'
plot \
  'Out.txt' using 1:4 title "net'(x)" with lines, \
  'Out.txt' using 1:5 title "net''(x)" with lines
