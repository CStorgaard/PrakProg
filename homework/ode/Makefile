##############################################################################
# Paths
##############################################################################
LIB_DIR = /home/fyn/repos/PrakProg/homework/lib

##############################################################################
# Output
##############################################################################
Out.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< > $@

##############################################################################
# Build main.exe
##############################################################################
main.exe: main.cs $(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.dll $(LIB_DIR)/spline.dll $(LIB_DIR)/ode.dll
	cp $(LIB_DIR)/vec.dll .
	cp $(LIB_DIR)/matrix.dll .
	cp $(LIB_DIR)/spline.dll .
	cp $(LIB_DIR)/ode.dll .
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/vec.dll -reference:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/spline.dll -reference:$(LIB_DIR)/ode.dll main.cs

$(LIB_DIR)/vec.dll: $(LIB_DIR)/vec.cs
	mcs -target:library -out:$(LIB_DIR)/vec.dll $(LIB_DIR)/vec.cs

$(LIB_DIR)/matrix.dll: $(LIB_DIR)/matrix.cs $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.cs

$(LIB_DIR)/spline.dll: $(LIB_DIR)/spline.cs $(LIB_DIR)/matrix.dll $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/spline.dll -reference:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/spline.cs

$(LIB_DIR)/ode.dll: $(LIB_DIR)/ode.cs $(LIB_DIR)/matrix.dll $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/ode.dll -reference:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/ode.cs

##############################################################################
# Data Generation Targets
##############################################################################
harmonic.txt: main.exe
	mono main.exe harmonic

friction.txt: main.exe
	mono main.exe friction

planet_orbit_gr.txt: main.exe
	mono main.exe planet_orbit_gr

data: harmonic.txt friction.txt planet_orbit_gr.txt

##############################################################################
# Plot target
##############################################################################
plot-harmonic: harmonic.txt
	gnuplot harmonic.gp

plot-friction: friction.txt
	gnuplot friction.gp

plot-planet_orbit_gr: planet_orbit_gr.txt
	gnuplot planet_orbit_gr.gp

plot-all: plot-harmonic plot-friction plot-planet_orbit_gr


##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe vec.dll matrix.dll
