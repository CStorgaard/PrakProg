# planet_orbit_gr_compare.gp

set terminal svg size 1000,400 fname 'Helvetica'
set output "planet_orbit_gr_compare.svg"

set multiplot layout 1,3 title "Planet Orbit GR: Three Scenarios"

###############################################################################
# Left Subplot: u(phi) for scenarios 1, 2, 3
###############################################################################
set xlabel "phi"
set ylabel "u"
set title "u(phi)"
plot "planet_orbit_gr.txt" using 1:2 with lines title "Scenario 1", \
     "planet_orbit_gr.txt" using 1:4 with lines title "Scenario 2", \
     "planet_orbit_gr.txt" using 1:6 with lines title "Scenario 3"

###############################################################################
# Middle Subplot: u'(phi) for scenarios 1, 2, 3
###############################################################################
set xlabel "phi"
set ylabel "u'"
set title "u'(phi)"
plot "planet_orbit_gr.txt" using 1:3 with lines title "Scenario 1", \
     "planet_orbit_gr.txt" using 1:5 with lines title "Scenario 2", \
     "planet_orbit_gr.txt" using 1:7 with lines title "Scenario 3"

###############################################################################
# Right Subplot: Phase portrait for each scenario, plotting u' vs. u
###############################################################################
set xlabel "u"
set ylabel "u'"
set title "Phase Portrait"
plot "planet_orbit_gr.txt" using 2:3 with dot title "Scenario 1", \
     "planet_orbit_gr.txt" using 4:5 with lines title "Scenario 2", \
     "planet_orbit_gr.txt" using 6:7 with lines title "Scenario 3"

unset multiplot
set output
