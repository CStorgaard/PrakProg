##############################################################################
# Paths
##############################################################################
LIB_DIR = ../lib

.PHONY: all clean test
all: main.exe Out.txt erf.svg

##############################################################################
# Run main.exe to capture all console output (integrals, log messages, etc.)
##############################################################################
Out.txt: main.exe
	@echo "Running main.exe → Out.txt"
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/integrator.dll
	@echo "Compiling main.exe"
	cp $(LIB_DIR)/integrator.dll .
	mcs -target:exe \
	    -out:main.exe \
	    -reference:$(LIB_DIR)/integrator.dll \
	    main.cs

##############################################################################
# Build integrator.dll
##############################################################################
$(LIB_DIR)/integrator.dll: $(LIB_DIR)/integrator.cs
	@echo "Compiling integrator.dll"
	mcs -target:library \
	    -out:$(LIB_DIR)/integrator.dll \
	    $(LIB_DIR)/integrator.cs

##############################################################################
# Plotting the error function
##############################################################################
# Your program will have already created erf.txt
erf.svg: erf.txt erf.gp
	@echo "Plotting erf.svg"
	gnuplot erf.gp

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe integrator.dll
