.PHONY: all clean

all: Out.txt

Out.txt: main.exe
	mono $< > $@

main.exe: main.cs 
	mcs -reference:System.Numerics.dll -out:$@ $<

clean:
	rm -f *.dll *.exe Out.txt