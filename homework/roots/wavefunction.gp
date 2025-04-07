# Wavefunction.gp
set terminal svg size 800,600
set output "Wavefunction.svg"

set title "Calculated and exact wavefunction"
set xlabel "r"
set ylabel "\psi (r)"
set key top right
set grid

plot "WavefunctionOutput.txt" using 1:2 with lines title "Numeric wavefunction", \
     "WavefunctionOutput.txt" using 1:3 with lines title "Actual Wavefunction"