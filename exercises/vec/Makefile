Out.txt: main.exe
	mono $< > $@
	# mono main.exe > Out.txt

main.exe: main.cs vec.dll
	mcs -target:exe -out:main.exe -reference:vec.dll main.cs

vec.dll: vec.cs
	mcs -target:library -out:vec.dll vec.cs

clean:
	rm -f main.exe Out.txt

test:
	echo rm -f main.exe Out.txt
