##############################################################################
# Default target
##############################################################################
.PHONY: all plot-all clean test
all: plot-all

##############################################################################
# Paths
##############################################################################
LIB_DIR = ../lib

##############################################################################
# Output (unused by plot targets but available if needed)
##############################################################################
Out.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.dll $(LIB_DIR)/spline.dll
	cp $(LIB_DIR)/vec.dll .
	cp $(LIB_DIR)/matrix.dll .
	cp $(LIB_DIR)/spline.dll .
	mcs -target:exe \
	    -out:main.exe \
	    -reference:$(LIB_DIR)/vec.dll \
	    -reference:$(LIB_DIR)/matrix.dll \
	    -reference:$(LIB_DIR)/spline.dll \
	    main.cs

$(LIB_DIR)/vec.dll: $(LIB_DIR)/vec.cs
	mcs -target:library \
	    -out:$(LIB_DIR)/vec.dll \
	    $(LIB_DIR)/vec.cs

$(LIB_DIR)/matrix.dll: $(LIB_DIR)/matrix.cs $(LIB_DIR)/vec.dll
	mcs -target:library \
	    -out:$(LIB_DIR)/matrix.dll \
	    -reference:$(LIB_DIR)/vec.dll \
	    $(LIB_DIR)/matrix.cs

$(LIB_DIR)/spline.dll: $(LIB_DIR)/spline.cs $(LIB_DIR)/matrix.dll $(LIB_DIR)/vec.dll
	mcs -target:library \
	    -out:$(LIB_DIR)/spline.dll \
	    -reference:$(LIB_DIR)/matrix.dll \
	    -reference:$(LIB_DIR)/vec.dll \
	    $(LIB_DIR)/spline.cs

##############################################################################
# Plot targets
##############################################################################
plot-linear: main.exe
	mono main.exe
	gnuplot linearSplinePlot.gp
	@echo "Linear spline plot generated as linearSplinePlot.svg"

plot-quadratic: main.exe
	mono main.exe
	gnuplot quadraticSplinePlot.gp
	@echo "Quadratic spline plot generated as quadraticSplinePlot.svg"

plot-cubic: main.exe
	mono main.exe
	gnuplot cubicSplinePlot.gp
	@echo "Cubic spline plot generated as cubicSplinePlot.svg"

plot-gnuplot: main.exe
	mono main.exe
	gnuplot gnuplotSplinePlot.gp
	@echo "Gnuplot spline plot generated as gnuplotSplinePlot.svg"

plot-all: plot-linear plot-quadratic plot-cubic plot-gnuplot

# Alias for backward compatibility
plot: plot-all

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe vec.dll matrix.dll spline.dll
