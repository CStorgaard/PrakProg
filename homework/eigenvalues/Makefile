##############################################################################
# Paths
##############################################################################
LIB_DIR = /home/fyn/repos/PrakProg/homework/lib

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.dll
	cp $(LIB_DIR)/vec.dll .
	cp $(LIB_DIR)/matrix.dll .
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/vec.dll -reference:$(LIB_DIR)/matrix.dll main.cs

$(LIB_DIR)/vec.dll: $(LIB_DIR)/vec.cs
	mcs -target:library -out:$(LIB_DIR)/vec.dll $(LIB_DIR)/vec.cs

$(LIB_DIR)/matrix.dll: $(LIB_DIR)/matrix.cs $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.cs

##############################################################################
# 1) Generate data for Δr convergence study (fix rmax=10)
##############################################################################
dr_data.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono main.exe -study:dr -rmax:10 > dr_data.txt

dr_plot.svg: dr_data.txt
	gnuplot -e "\
	set terminal svg; \
	set output 'dr_plot.svg'; \
	set title 'Convergence with Δr (rmax fixed)'; \
	set xlabel 'Δr'; \
	set ylabel 'ε0'; \
	plot 'dr_data.txt' using 1:2 with linespoints title 'ε0 vs Δr'"

##############################################################################
# 2) Generate data for rmax convergence study (fix dr=0.3)
##############################################################################
rmax_data.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono main.exe -study:rmax -dr:0.3 > rmax_data.txt

rmax_plot.svg: rmax_data.txt
	gnuplot -e "\
	set terminal svg; \
	set output 'rmax_plot.svg'; \
	set title 'Convergence with rmax (Δr fixed)'; \
	set xlabel 'rmax'; \
	set ylabel 'ε0'; \
	plot 'rmax_data.txt' using 1:2 with linespoints title 'ε0 vs rmax'"

##############################################################################
# 3) Generate wavefunction data (assuming main.exe -study:wavefunction prints
#    (r[i], f0, f1, f2, ...) in columns)
##############################################################################
wavefunction_data.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono main.exe -study:wavefunction > wavefunction_data.txt

wavefunction_plot.svg: wavefunction_data.txt
	gnuplot -e "\
	set terminal svg; \
	set output 'wavefunction_plot.svg'; \
	set title 'Hydrogen s-wave radial wavefunctions'; \
	set xlabel 'r'; \
	set ylabel 'f(r)'; \
	plot 'wavefunction_data.txt' using 2:3 with lines title 'f0 numeric', \
	     'wavefunction_data.txt' using 2:4 with lines title 'f0 analytic', \
	     'wavefunction_data.txt' using 2:5 with lines title 'f1 numeric', \
	     'wavefunction_data.txt' using 2:6 with lines title 'f1 analytic', \
	     'wavefunction_data.txt' using 2:7 with lines title 'f2 numeric', \
	     'wavefunction_data.txt' using 2:8 with lines title 'f2 analytic'"


##############################################################################
# Combined targets
##############################################################################
# Make both Δr and rmax plots:
plot: dr_plot.svg rmax_plot.svg
	@echo "Plots generated: dr_plot.svg and rmax_plot.svg"

# Make wavefunction plot
wavefunction: wavefunction_plot.svg
	@echo "Wavefunction plot generated: wavefunction_plot.svg"

# Make everything (studies + wavefunction)
all: main.exe dr_plot.svg rmax_plot.svg wavefunction_plot.svg

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe dr_data.txt rmax_data.txt wavefunction_data.txt dr_plot.svg rmax_plot.svg wavefunction_plot.svg
