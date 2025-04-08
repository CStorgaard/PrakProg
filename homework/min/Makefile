##############################################################################
# Default target: build data, main.exe, generate plot, and output file.
##############################################################################
all: data main.exe higgs.svg Out.txt

##############################################################################
# Paths
##############################################################################
LIB_DIR = ../lib

##############################################################################
# Output: Run main.exe and save its output.
##############################################################################
Out.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/min.dll
	# Copy the integrator library to the current directory
	cp $(LIB_DIR)/min.dll .
	# Compile main.cs referencing only $(LIB_DIR)/min.dll
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/min.dll main.cs

##############################################################################
# Data file: generate the experimental data file
##############################################################################
data: higgs.data.txt

higgs.data.txt:
	@echo "# energy E[GeV], signal σ(E) [certain units], experimental uncertainty Δσ [same units]" > higgs.data.txt
	@echo "101 -0.25 2.0" >> higgs.data.txt
	@echo "103 -0.30 2.0" >> higgs.data.txt
	@echo "105 -0.15 1.9" >> higgs.data.txt
	@echo "107 -1.71 1.9" >> higgs.data.txt
	@echo "109 0.81 1.9" >> higgs.data.txt
	@echo "111 0.65 1.9" >> higgs.data.txt
	@echo "113 -0.91 1.9" >> higgs.data.txt
	@echo "115 0.91 1.9" >> higgs.data.txt
	@echo "117 0.96 1.6" >> higgs.data.txt
	@echo "119 -2.52 1.6" >> higgs.data.txt
	@echo "121 -1.01 1.6" >> higgs.data.txt
	@echo "123 2.01 1.6" >> higgs.data.txt
	@echo "125 4.83 1.6" >> higgs.data.txt
	@echo "127 4.58 1.6" >> higgs.data.txt
	@echo "129 1.26 1.3" >> higgs.data.txt
	@echo "131 1.01 1.3" >> higgs.data.txt
	@echo "133 -1.26 1.3" >> higgs.data.txt
	@echo "135 0.45 1.3" >> higgs.data.txt
	@echo "137 0.15 1.3" >> higgs.data.txt
	@echo "139 -0.91 1.3" >> higgs.data.txt
	@echo "141 -0.81 1.1" >> higgs.data.txt
	@echo "143 -1.41 1.1" >> higgs.data.txt
	@echo "145 1.36 1.1" >> higgs.data.txt
	@echo "147 0.50 1.1" >> higgs.data.txt
	@echo "149 -0.45 1.1" >> higgs.data.txt
	@echo "151 1.61 1.1" >> higgs.data.txt
	@echo "153 -2.21 1.1" >> higgs.data.txt
	@echo "155 -1.86 0.9" >> higgs.data.txt
	@echo "157 1.76 0.9" >> higgs.data.txt
	@echo "159 -0.50 0.9" >> higgs.data.txt

##############################################################################
# Build min.dll
##############################################################################
$(LIB_DIR)/min.dll: $(LIB_DIR)/min.cs $(LIB_DIR)/vec.cs $(LIB_DIR)/matrix.cs $(LIB_DIR)/ode.cs
	mcs -target:library -out:$(LIB_DIR)/min.dll $(LIB_DIR)/min.cs $(LIB_DIR)/vec.cs $(LIB_DIR)/matrix.cs $(LIB_DIR)/ode.cs

##############################################################################
# Plotting: generate the Higgs.svg plot using gnuplot and plot.dat
##############################################################################
higgs.svg: plot.dat higgs.gp
	gnuplot higgs.gp

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe $(LIB_DIR)/min.dll
