##############################################################################
# Paths
##############################################################################
LIB_DIR = ../lib

##############################################################################
# Output
##############################################################################
Out.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/MonteCarlo.dll
	# Copy the integrator library to the current directory
	cp $(LIB_DIR)/MonteCarlo.dll .
	# Compile main.cs referencing only integrator.dll
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/MonteCarlo.dll main.cs

##############################################################################
# Build integrator.dll
##############################################################################
$(LIB_DIR)/MonteCarlo.dll: $(LIB_DIR)/montecarlo.cs
	mcs -target:library -out:$(LIB_DIR)/MonteCarlo.dll $(LIB_DIR)/montecarlo.cs

##############################################################################
# Plotting
##############################################################################
pi.svg: pi_results.txt pi.gp
	gnuplot pi.gp

quasi_pi.svg: pi_results_quasi.txt quasi_pi.gp
	gnuplot quasi_pi.gp

all: pi.svg quasi_pi.svg
	@echo "Plots generated: pi.svg and quasi_pi.svg"
##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe integrator.dll
