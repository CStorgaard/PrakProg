LIB_DIR = ../lib

Out.txt: main.exe
	MONO_PATH=$(LIB_DIR) mono $< > $@

main.exe: main.cs $(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.dll
	cp $(LIB_DIR)/vec.dll .
	cp $(LIB_DIR)/matrix.dll .
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/vec.dll -reference:$(LIB_DIR)/matrix.dll main.cs

$(LIB_DIR)/vec.dll: $(LIB_DIR)/vec.cs
	mcs -target:library -out:$(LIB_DIR)/vec.dll $(LIB_DIR)/vec.cs

$(LIB_DIR)/matrix.dll: $(LIB_DIR)/matrix.cs $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.cs

clean:
	rm -f *.txt *.dll *.exe

test:
	echo rm -f main.exe Out.txt
