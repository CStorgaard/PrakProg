##############################################################################
# Default target: build library, all executables, run them, and generate plots.
##############################################################################

.PHONY: all clean

all: \
    Berrut.dll \
    main.exe Out.txt interpolation.svg \
    runge.exe Out_runge.txt runge.svg \
    time.exe Out_times.txt timing.svg \
    error.exe Out_error.txt error.svg

##############################################################################
# Library source
##############################################################################
LIB_SRCS := Berrut.cs
BERRUT_DLL := Berrut.dll

##############################################################################
# Build Berrut.dll
##############################################################################
$(BERRUT_DLL): $(LIB_SRCS)
	@mcs -target:library -out:$@ $^

##############################################################################
# Build all executables
##############################################################################
main.exe:      main.cs       $(BERRUT_DLL)
	@mcs -target:exe -out:$@ -reference:$(BERRUT_DLL) $<

runge.exe:     main_runge.cs $(BERRUT_DLL)
	@mcs -target:exe -out:$@ -reference:$(BERRUT_DLL) $<

time.exe:      main_time.cs  $(BERRUT_DLL)
	@mcs -target:exe -out:$@ -reference:$(BERRUT_DLL) $<

error.exe:     main_error.cs $(BERRUT_DLL)
	@mcs -target:exe -out:$@ -reference:$(BERRUT_DLL) $<

##############################################################################
# Run the programs under Mono, capture outputs
##############################################################################
Out.txt:         main.exe
	@MONO_PATH=. mono $< > $@

Out_runge.txt:   runge.exe
	@MONO_PATH=. mono $< > $@

Out_times.txt:   time.exe
	@MONO_PATH=. mono $< > $@

Out_error.txt:   error.exe
	@MONO_PATH=. mono $< > $@

##############################################################################
# Generate plots via Gnuplot
##############################################################################
interpolation.svg: Out.txt        interpolation.gp
	@echo "-> Generating $@ …"
	@gnuplot interpolation.gp

runge.svg:         Out_runge.txt  runge.gp
	@echo "-> Generating $@ …"
	@gnuplot runge.gp

timing.svg:       Out_times.txt  timing.gp
	@echo "-> Generating $@ …"
	@gnuplot timing.gp

error.svg:         Out_error.txt   error.gp
	@echo "-> Generating $@ …"
	@gnuplot error.gp

##############################################################################
# Clean up all build artifacts and outputs
##############################################################################
clean:
	@rm -f *.dll *.exe *.txt *.svg
