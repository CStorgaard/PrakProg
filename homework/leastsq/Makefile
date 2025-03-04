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
main.exe: main.cs $(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.dll
	cp $(LIB_DIR)/vec.dll .
	cp $(LIB_DIR)/matrix.dll .
	mcs -target:exe -out:main.exe -reference:$(LIB_DIR)/vec.dll -reference:$(LIB_DIR)/matrix.dll main.cs

$(LIB_DIR)/vec.dll: $(LIB_DIR)/vec.cs
	mcs -target:library -out:$(LIB_DIR)/vec.dll $(LIB_DIR)/vec.cs

$(LIB_DIR)/matrix.dll: $(LIB_DIR)/matrix.cs $(LIB_DIR)/vec.dll
	mcs -target:library -out:$(LIB_DIR)/matrix.dll -reference:$(LIB_DIR)/vec.dll $(LIB_DIR)/matrix.cs

##############################################################################
# Plot target
##############################################################################
plot: main.exe
    #1) Run the program to produce dataPoints.txt and fitPoints.txt
	mono main.exe
    # 2) Call Gnuplot with the script
	gnuplot plot.gp
	@echo "Plot generated as decayFit.svg"

##############################################################################
# Housekeeping
##############################################################################
clean:
	rm -f *.txt *.dll *.exe *.svg

test:
	@echo rm -f main.exe vec.dll matrix.dll
