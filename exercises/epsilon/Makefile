Out.txt: main.exe
	mono $< > $@
	# mono main.exe > Out.txt

main.exe: main.cs doublecompare.dll
	mcs -target:exe -out:main.exe -reference:doublecompare.dll main.cs


doublecompare.dll : doublecompare.cs
	mcs -target:library -out:doublecompare.dll doublecompare.cs

clean:
	rm -f main.exe Out.txt

test:
	echo rm -f main.exe Out.txt
