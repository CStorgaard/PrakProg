CXXFLAGS += -std=c++23 -Wall
LDLIBS += -lstdc++
CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: plot.test.svg

plot.erf.svg: data.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set key out;\
	set tics out;\
	set title "Error function";\
	set xlabel "x";\
	set ylabel "y";\
	plot "$<" using 1:2 index 0 with lines title "Error function"\
	,"erf.data" using 1:2 index 0 with points title "Error function data"\
	' | tee log.erf.gpi | gnuplot

plot.gamma.svg: data.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set key out;\
	set tics out;\
	set title "Γ function";\
	set xlabel "x";\
	set ylabel "";\
	plot "$<" using 1:2 index 1 with lines title "Γ function",\
	     "$<" using 1:2 index 3 with points title "factorials"\
	' | tee log.gamma.gpi | gnuplot

plot.lngamma.svg: data.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set key out;\
	set tics out;\
	set title "LogΓ function";\
	set xlabel "x";\
	set ylabel "";\
	plot "$<" using 1:2 index 2 with lines title "LogΓ function",\
	     "$<" using 1:2 index 4 with points title " log factorials"\
	' | tee log.lngamma.gpi | gnuplot

data.txt: main.exe
	mono $< > $@

plot.test.svg: plot.gpi
	gnuplot $< > $@  

sfuns.dll: sfuns.cs 
	$(MKLIB)  

main.exe: main.cs sfuns.dll
	$(MKEXE)  

clean: 
	$(RM) *.dll *.exe Out.txt [lL]og* *.svg data.txt