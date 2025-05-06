# plot_func.gp
set terminal svg size 800,600
set output 'func.svg'
set title 'g(x) vs. Neural Net Response'
set xlabel 'x'
set ylabel 'y'
plot \
  'Out.txt' using 1:2 with lines title 'g(x) + 0.1', \
  'Out.txt' using 1:3 with lines title 'net(x)'