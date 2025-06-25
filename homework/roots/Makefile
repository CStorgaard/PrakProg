##############################################################################
# Default target
##############################################################################
.PHONY: all plot-all clean test
all: Out.txt WavefunctionOutput.txt ConvergenceStudy.txt wavefunction.svg plot-all

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
# Data generation
##############################################################################
WavefunctionOutput.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< wavefunction > $@

ConvergenceStudy.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< convergence > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/roots.dll
	cp $(LIB_DIR)/roots.dll .
	mcs -target:exe \
	    -out:main.exe \
	    -reference:$(LIB_DIR)/roots.dll \
	    main.cs

##############################################################################
# Build roots.dll
##############################################################################
$(LIB_DIR)/roots.dll: $(LIB_DIR)/roots.cs \
                     $(LIB_DIR)/vec.cs \
                     $(LIB_DIR)/matrix.cs \
                     $(LIB_DIR)/ode.cs
	mcs -target:library \
	    -out:$(LIB_DIR)/roots.dll \
	    $(LIB_DIR)/roots.cs \
	    $(LIB_DIR)/vec.cs \
	    $(LIB_DIR)/matrix.cs \
	    $(LIB_DIR)/ode.cs

##############################################################################
# Plotting
##############################################################################
wavefunction.svg: WavefunctionOutput.txt wavefunction.gp
	gnuplot wavefunction.gp

convergence_min.svg: ConvergenceStudy.txt convergence_min.gp
	gnuplot convergence_min.gp

convergence_max.svg: ConvergenceStudy.txt convergence_max.gp
	gnuplot convergence_max.gp

convergence_eps.svg: ConvergenceStudy.txt convergence_eps.gp
	gnuplot convergence_eps.gp

convergence_acc.svg: ConvergenceStudy.txt convergence_acc.gp
	gnuplot convergence_acc.gp

plot-all: wavefunction.svg \
          convergence_min.svg \
          convergence_max.svg \
          convergence_eps.svg \
          convergence_acc.svg
	@echo "Plots generated: wavefunction.svg, convergence_*.svg"

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe roots.dll