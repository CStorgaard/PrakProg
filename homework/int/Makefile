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
main.exe: main.cs $(LIB_DIR)/integrator.dll
	# Copy the integrator library to the current directory
	cp $(LIB_DIR)/integrator.dll .
	# Compile main.cs referencing only integrator.dll
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/integrator.dll main.cs

##############################################################################
# Build integrator.dll
##############################################################################
$(LIB_DIR)/integrator.dll: $(LIB_DIR)/integrator.cs
	mcs -target:library -out:$(LIB_DIR)/integrator.dll $(LIB_DIR)/integrator.cs

##############################################################################
# Plotting
##############################################################################
erf.svg: erf.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set title "Error function";\
	set xlabel "x";\
	set ylabel "erf(x)";\
	set key top left;\
	set grid xtics ytics;\
	plot "erf.txt" using 1:2:3 with yerrorlines title "erf(x)"\
	' | gnuplot > $@

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe integrator.dll
