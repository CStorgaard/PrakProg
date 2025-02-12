.PHONY: all clean

# The default target builds main.exe and then runs the program to create Out.txt
all: Out.txt

# Out.txt depends on main.exe and input.txt
Out.txt: main.exe input.txt
	@echo "Running main.exe with arguments from input.txt..."
	@mono main.exe $(shell cat input.txt) > Out.txt

# Compile main.exe from main.cs
main.exe: main.cs
	mcs -target:exe -out:main.exe main.cs

# Create an input file with the command-line argument
input.txt:
	@echo "-numbers:1,2,3,4,5" > input.txt

# Clean up generated files
clean:
	rm -f main.exe Out.txt input.txt