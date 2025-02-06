Out.txt : main.exe
	mono $< >$@
	# mono main.exe > Out.txt

main.exe : main.cs sfuns.dll
	mcs -target:exe -reference:sfuns.dll -out:main.exe main.cs

sfuns.dll : sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs

clean:
	$(RM) sfuns.dll main.exe Out.txt

test:
	echo $(RM)
