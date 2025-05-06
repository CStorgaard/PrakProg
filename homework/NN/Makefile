##############################################################################
# Default target: build data, main.exe, generate plot, and output file.
##############################################################################

.PHONY: all clean test

all: main.exe Out.txt func.svg diff.svg int.svg

##############################################################################
# Paths
##############################################################################
LIB_DIR = ../lib

##############################################################################
# Library sources
##############################################################################
LIB_SRCS :=	$(LIB_DIR)/GradientDescent.cs \
            $(LIB_DIR)/vec.cs \
            $(LIB_DIR)/matrix.cs \
            $(LIB_DIR)/ode.cs \
            $(LIB_DIR)/ANN.cs

##############################################################################
# Build ANN.dll
##############################################################################
$(LIB_DIR)/ANN.dll: $(LIB_SRCS)
	mcs -target:library -out:$@ $^

##############################################################################
# Build main.exe 
##############################################################################
main.exe: main.cs $(LIB_DIR)/ANN.dll
	cp $(LIB_DIR)/ANN.dll .
	mcs -target:exe -out:$@ -reference:ANN.dll $<

##############################################################################
# Output: Run main.exe and save its output.
##############################################################################
Out.txt: main.exe
	# run the program under Mono, capture output
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Plot function vs. network
##############################################################################
func.svg: Out.txt plot_func.gp
	@echo "-> Generating func.svg via Gnuplot…"
	gnuplot plot_func.gp

##############################################################################
# Plot first derivatives
##############################################################################
diff.svg: Out.txt plot_diff.gp
	@echo "-> Generating diff.svg via Gnuplot…"
	gnuplot plot_diff.gp

##############################################################################
# Plot antiderivatives
##############################################################################
int.svg: Out.txt plot_int.gp
	@echo "-> Generating int.svg via Gnuplot…"
	gnuplot plot_int.gp


##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	rm -f main.exe $(LIB_DIR)/ANN.dll
